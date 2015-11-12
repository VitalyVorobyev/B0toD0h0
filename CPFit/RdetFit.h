#ifndef RDETFIT_H
#define RDETFIT_H

#include "RkRdetRnpPdf.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/FunctionMinimum.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>

#include "TLine.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

#include "RkRdetRnpPdf.h"

RkRdetRnpPdf* m_pdf;
TChain* m_tree;
TTree* m_good_tree;
Int_t m_exp,m_bin,m_flv;
Int_t m_b0f;
Double_t m_dt;
Double_t m_z_sig;
Double_t m_z_sig_mc;
Double_t m_A,m_B;
Double_t m_costhBcms;
Double_t m_chi2_vtx_d0;
int m_svd = 2;

// * Event-by-event parameters
int ntrk_rec;
int ntrk_asc;
int ndf_rec;
int ndf_asc;
int good_icpv;
double sz_rec;
double sz_asc;
double chisq_rec;
double chisq_asc;
// *

using namespace std;
using namespace ROOT;
using namespace Minuit2;

const double sol = 2.99792458;
bool draw_plots = true;
const double cm2ps = 78.48566945838871754705;
int xi;
const int NDots = 250;
const int NBins = NDots;
double dtmax = 70;
double dtmin = -dtmax;
//const double dzmax =  dtmax/78.4857;
//const double dzmin = -dzmax;

int m_mode = 4;
bool make_Rdet_fit = false;
bool make_Olr_fit  = false;
bool sgl_trk       = false;
bool mlt_trk       = false;
bool include_s0    = false;
bool include_dt0   = false;
bool pipi_fit      = false;
bool d0_fit        = false;
bool ddalitz       = false;

int init(const int _mode){
  cout << "init: mode " << _mode << endl;
  m_pdf->SetRange(dtmax);
  return 0;
}

void GetEvent(const int i){
  m_good_tree->GetEvent(i);
  m_dt = (m_z_sig-m_z_sig_mc)*0.1*cm2ps;
  sz_rec /= 10.;
  return;
}

void SetPDFParams(const vector<double>& par){
  int i=0;
  // * Outlayer * //
  m_pdf->Set_f_ol_sgl(par[i++]);
  m_pdf->Set_f_ol_mlt(par[i++]);
  m_pdf->Set_sigma_ol(par[i++]);

  // * multiple Rdet * //
  m_pdf->Set_Srec(par[i],par[i+1]); i+=2;

  // * single Rdet * //
  m_pdf->Set_Smn_rec(par[i++]);
  m_pdf->Set_Stl_rec(par[i++]);
  m_pdf->Set_ftl_rec(par[i++]);
  return;
}

TTree* GetGoodTTree(TTree* tree, const int bin = 0, const int flv = 0){
  stringstream out;
  out.str("");
  out << "z_sig_mc>-1";
  if(m_svd == 2) out << " && exp>30";
  else           out << " && exp<30";
  if(bin)        out << " && bin_mc == " << bin;
  if(flv)        out << " && flv_mc == " << flv;
  cout << out.str() << endl;
  return tree->CopyTree(out.str().c_str());
}

class pdfFcn : public FCNBase{
public:
  pdfFcn(void){
    theErrorDef = 1;
    NTot = m_good_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcn() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double loglh = 0;
    double sigpdf;
    SetPDFParams(par);
    for(int i=0; i<NTot; i++){
      GetEvent(i);
      if(include_s0){ sz_rec += par[8];}
      if(include_dt0){ m_dt  += par[9]*chisq_rec/ndf_rec;}
      sigpdf = m_pdf->PdfRrec(m_dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec);
      if(!std::isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf);
      else{ cout << "pdfFcn: pdf = " << sigpdf << endl;}
    }
    cout << "loglh: " << loglh;
    if(mlt_trk) cout << ", Srec0: " << par[0+3] << ", Srec1: " << par[1+3];
    else        cout << ", Smn: " << par[2+3] << ", Stl: " << par[3+3] << ", ftl: " << par[4+3];
    if(include_s0) cout << ", s0 = " << par[8];
    if(include_dt0) cout << ", dt0 = " << par[9];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

int Mode(const int mode){
  switch (mode) {
  case 1:  return 1; // pi0
  case 10: return 10;// D*0 pi0
  case 20: return 20;// D*0 eta
  case 4:  return 3; // omega
  case 5:  return 5; // eta'
  default: return 2; // eta
  }
}

int h0Mode(const int mode){
  switch (mode) {
  case 1:  return 10; // pi0
  case 10: return 10; // D*0 pi0
  case 2:  return 10; // eta->gg
  case 20: return 10; // D*0 eta
  case 3:  return 20; // eta->ppp
  case 4:  return 20; // omega
  case 5:  return 10; // eta'
  }
}

bool is_sgl_vertex(const int mode){
  switch (mode) {
  case 1:  return true;  // pi0
  case 10: return true;  // D*0 pi0
  case 2:  return true;  // eta->gg
  case 20: return true;  // D*0 eta
  case 3:  return false; // eta->ppp
  case 4:  return false; // omega
  case 5:  return false; // eta'
  }
}

string GenFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/mixtree_m" << Mode(mode) << "_mh0" << h0Mode(mode) << ".root";
  return out.str();
}

#endif
