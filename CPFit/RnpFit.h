#ifndef RNPFIT_H
#define RNPFIT_H

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
Double_t m_z_asc;
Double_t m_z_asc_mc;
Double_t m_A,m_B;
Double_t m_costhBcms;
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
bool make_Rnp_fit  = false;
bool make_Rdet_fit = false;
bool make_Olr_fit  = false;
bool make_sgl_fit  = false;
bool make_mlt_fit  = true;
bool include_s0    = false;

int init(const int _mode){
  cout << "init: mode " << _mode << endl;
  m_pdf->SetRange(dtmax);
  return 0;
}

void GetEvent(const int i){
  m_good_tree->GetEvent(i);
  m_dt = (m_z_asc-m_z_asc_mc)*0.1*cm2ps;
//  m_dt *= 0.1*cm2ps;
  sz_asc /= 10.;
  return;
}

void SetPDFParams(const vector<double>& par){
  int i=0;
  // * Outlayer * //
  m_pdf->Set_f_ol_sgl(par[i++]);
  m_pdf->Set_f_ol_mlt(par[i++]);
  m_pdf->Set_sigma_ol(par[i++]);

  // * multiple Rdet * //
  m_pdf->Set_Sasc(par[i],par[i+1]); i+=2;

  // * single Rdet * //
  m_pdf->Set_Smn_asc(par[i++]);
  m_pdf->Set_Stl_asc(par[i++]);
  m_pdf->Set_ftl_asc(par[i++]);

  // * single Rnp * //
  m_pdf->Set_fd_np_sgl(par[i],par[i+1]); i+=2;
  m_pdf->Set_fp_np_sgl(par[i++]);
  m_pdf->Set_tau_np_p_sgl(par[i],par[i+1]); i+=2;
  m_pdf->Set_tau_np_n_sgl(par[i],par[i+1]); i+=2;

  // * multiple Rnp */
  m_pdf->Set_fd_np_mlt(par[i],par[i+1]); i+=2;
  m_pdf->Set_fd_np_st_mlt(par[i++]);
  m_pdf->Set_fd_np_xi_mlt(par[i++]);
  m_pdf->Set_fd_np_stxi_mlt(par[i++]);
  m_pdf->Set_fp_np_mlt(par[i++]);
  m_pdf->Set_fn_np_mlt(par[i++]);
  m_pdf->Set_tau_np_p_mlt(par[i],par[i+1]); i+=2;
  m_pdf->Set_tau_np_p_xi_mlt(par[i++]);
  m_pdf->Set_tau_np_p_stxi_mlt(par[i++]);
  m_pdf->Set_tau_np_n_mlt(par[i],par[i+1]); i+=2;
  m_pdf->Set_tau_np_n_xi_mlt(par[i++]);
  m_pdf->Set_tau_np_n_stxi_mlt(par[i++]);

  return;
}

TTree* GetGoodTTree(TTree* tree, const int bin = 0, const int flv = 0){
  stringstream out;
  out.str("");
  out << "z_asc_mc>-1";
  if(m_svd == 2)   out << " && exp>30";
  else             out << " && exp<30";
  if(bin)          out << " && bin_mc == " << bin;
  if(flv)          out << " && flv_mc == " << flv;
  if(make_sgl_fit) out << " && ndf_asc==0";
  if(make_mlt_fit) out << " && ndf_asc>0";
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
//    double norm = 0;
//    double Norm = 0;
//    double pdfval;
//    const double ddt = (dtmax - dtmin)/100.;
//    double dt;
    SetPDFParams(par);
    for(int i=0; i<NTot; i++){
      GetEvent(i);
      if(include_s0) { sz_asc += par[par.size()-1];}
//      for(int j=0; j<100; j++){
//        dt = dtmin + (j+0.5)*ddt;
//        pdfval = m_pdf->RascRnp(dt,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
//        if(!std::isnan(pdfval)) norm += pdfval;
//      }
//      norm *= ddt;
//      Norm += norm;
      sigpdf = m_pdf->PdfRascRnp(m_dt,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
//      if(!std::isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf/norm);
      if(!std::isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf);
      else{ cout << "pdfFcn: pdf = " << sigpdf << endl;}
    }
    cout << "loglh: " << loglh;// << ", Norm: " << Norm;// << ", tau: " << par[0] << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3];
    if(make_mlt_fit) cout << ", Srec0: " << par[0+3] << ", Srec1: " << par[1+3];
    else             cout << ", Smn: " << par[2+3] << ", Stl: " << par[3+3] << ", ftl: " << par[4+3];
    if(include_s0) cout << ", s0: " << par[par.size()-1];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

double calc_norm_single(const int ntrk_rec,const double& sz_rec,const double& chisq_z_rec,const int ndf_z_rec,const int ntrk_asc,const double& sz_asc,const double& chisq_z_asc,const int ndf_z_asc, const bool otlr = true){
  double norm = 0;
  const double ddt = (dtmax - dtmin)/100.;
  double dt;
  for(int i=0; i<100; i++){
    dt = dtmin + (i+0.5)*ddt;
    double pdf = m_pdf->RascRnp(m_dt,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
    if(pdf>0) norm += pdf;
  }
  return norm*ddt;
}

class pdfFcnSingle : public FCNBase{
public:
  pdfFcnSingle(void){
    theErrorDef = 1;
    NTot = m_good_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnSingle() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    SetPDFParams(par);
    double loglh = 0;
    double sigpdf;
    double sum_norm = 0;
    for(int i=0; i<NTot; i++){
//      m_good_tree->GetEvent(i);
//      m_dt *= 10./sol*cm2ps;
      GetEvent(i);
      const double norm = 1;//calc_norm_single(ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,false);
      sigpdf = m_pdf->RascRnp(m_dt,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
      sum_norm += norm/NTot;
      if(!std::isnan(sigpdf) && sigpdf>0){ loglh += -2*TMath::Log(sigpdf/norm);}
      else{ cout << "pdfFcnSingle: pdf = " << sigpdf << endl;}
    }
    cout << "loglh: " << loglh << ", norm = " << sum_norm << ", tau: " << par[0] << ", dm = " << par[1] << ", A: " << par[2] << ", B: " << par[3] << endl;
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
