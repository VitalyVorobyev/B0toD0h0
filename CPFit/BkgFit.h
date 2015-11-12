#ifndef BKGFIT_H
#define BKGFIT_H

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

#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"

#include "RbkgPdf.h"

const double separate_plots = false;
const double fix_all   = false;
const double only_mlt  = false;
const double only_sgl  = false;

const double fix_delta = false;
const double fix_lt    = false;
const double fix_main  = false;
const double fix_tail  = false;

const double fix_mlt   = false || only_sgl || fix_all;
const double fix_sgl   = false || only_mlt || fix_all;

RbkgPdf* m_pdf;
TTree* m_tree;
TTree* m_good_tree;
Int_t m_exp,m_bin,m_flv;
Double_t m_tag_LH;
Double_t m_dt;
Double_t m_z_sig, m_z_asc;
const int m_svd = 2;

// * Event-by-event parameters
int ntrk_rec;
int ntrk_asc;
int ndf_rec;
int ndf_asc;
double sz_rec;
double sz_asc;
// *

void SetPdfParams(const vector<double>& par){
  m_pdf->SetTau(par[0]);
  m_pdf->Set_mu(par[1]);
  m_pdf->Set_mu_delta(par[2]);

  m_pdf->Set_f_delta_mlt(par[3]);
  m_pdf->Set_f_tail_mlt(par[4]);
  m_pdf->Set_S_main_mlt(par[5]);
  m_pdf->Set_S_tail_mlt(par[6]);

  m_pdf->Set_f_delta_sgl(par[7]);
  m_pdf->Set_f_tail_sgl(par[8]);
  m_pdf->Set_S_main_sgl(par[9]);
  m_pdf->Set_S_tail_sgl(par[10]);
  return;
}

using namespace std;
using namespace ROOT;
using namespace Minuit2;

bool draw_plots = true;
const double sol = 2.99792458;
const double cm2ps = 78.48566945838871754705;

const int NDots = 250;
const int NBins = 250;
const double dtmax = 70;
const double dtmin = -dtmax;

int bin(const int j){
  if(j<8) return j-8;
  else    return j-7;
}
int flv_ind(const int flv){
  if(flv == 1) return 0;
  else         return 1;
}
int bin_ind(const int bin){
  if(bin<0) return bin+8;
  else      return bin+7;
}
int flv(const int k){
  if(k) return -1;
  else  return  1;
}

double sum_sigma(const double& s1, const double& s2){
  return sqrt(s1*s1+s2*s2);
}

int init(const int _mode){
  cout << "init: mode " << _mode << endl;
  m_pdf->SetRange(dtmax);
  return 0;
}

void GetEvent(const int i){
  m_good_tree->GetEvent(i);
  m_tag_LH>0 ? m_flv = 1 : m_flv = -1;
  m_dt = (m_z_sig - m_z_asc)*0.1*cm2ps;
  sz_rec /= 10.;
  sz_asc /= 10.;
  return;
}

TTree* GetGoodTTree(TTree* tree, const int bin = 0, const int flv = 0){
  stringstream out;
  out.str("");
  out << "z_sig_mc<-98";
  if(m_svd == 2) out << " && exp>30";
  else           out << " && exp<30";
  if(bin)      out << " && bin == " << bin;
  if(flv)      out << " && flv_mc == " << flv;
  if(only_mlt) out << " && ndf_asc!=0";
  if(only_sgl) out << " && ndf_asc==0";
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
    SetPdfParams(par);

    double loglh = 0;
    double pdf;
    double norm = 0;

    for(int i=0; i<NTot; i++){
      GetEvent(i);
      const double s = sum_sigma(sz_asc,sz_rec);
      pdf = m_pdf->Pdf(m_dt,s,ndf_asc);
      loglh += -2*TMath::Log(pdf);
    }
    cout << "loglh: " << loglh;
    for(int i=0; i<10; i++) cout << " " << par[i];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

int draw_fit_results(const vector< vector<double> >& vals,const vector< vector<double> >& errs, const vector<double> Aref, const vector<double> Bref){
  const int NFits = vals.size();
  double Fits[NFits],FitsErr[NFits];
  double TAU[NFits],TAUErr[NFits];
  double DM[NFits],DMErr[NFits];
  double A[NFits],AErr[NFits];
  double B[NFits],BErr[NFits];
  double ARef[NFits],BRef[NFits];

  for(int i=0; i<NFits; i++){
    if(NFits == 16) { Fits[i] = bin(i); FitsErr[i] = 0;}
    else            { Fits[i] = i+1; FitsErr[i] = 0;}
    TAU[i] = vals[i][0]; TAUErr[i] = errs[i][0];
    DM[i]  = vals[i][1]; DMErr[i]  = errs[i][1];
    A[i]   = vals[i][2]; AErr[i]   = errs[i][2];
    B[i]   = vals[i][3]; BErr[i]   = errs[i][3];
    ARef[i]= Aref[i];    BRef[i]   = Bref[i];
  }

  TMultiGraph* mgA = new TMultiGraph("mgA","mgA");
  TMultiGraph* mgB = new TMultiGraph("mgB","mgB");
  TGraphErrors* grTAU = new TGraphErrors(NFits,Fits,TAU,FitsErr,TAUErr);
  grTAU->SetMarkerStyle(20);
  grTAU->SetMarkerSize(1);
  grTAU->SetLineStyle(0);
  grTAU->SetMarkerColor(kBlue);
  TGraphErrors* grDM = new TGraphErrors(NFits,Fits,DM,FitsErr,DMErr);
  grDM->SetMarkerStyle(20);
  grDM->SetMarkerSize(1);
  grDM->SetMarkerColor(kBlue);
  TGraphErrors* grA = new TGraphErrors(NFits,Fits,A,FitsErr,AErr);
  grA->SetMarkerStyle(20);
  grA->SetMarkerSize(1.);
  grA->SetMarkerColor(kBlue);
  TGraph* grAref = new TGraphErrors(NFits,Fits,ARef);
  grAref->SetMarkerStyle(20);
  grAref->SetMarkerSize(1.);
  grAref->SetMarkerColor(kRed);
  TGraphErrors* grB = new TGraphErrors(NFits,Fits,B,FitsErr,BErr);
  grB->SetMarkerStyle(20);
  grB->SetMarkerSize(1.);
  grB->SetMarkerColor(kBlue);
  mgA->Add(grA);
  mgA->Add(grAref);
  TGraph* grBref = new TGraphErrors(NFits,Fits,BRef);
  grBref->SetMarkerStyle(20);
  grBref->SetMarkerSize(1.);
  grBref->SetMarkerColor(kRed);
  mgB->Add(grB);
  mgB->Add(grBref);

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  c2->Draw();
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.49);
  pad1->Draw();
  TPad* pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.49);
  pad2->Draw();
  TPad* pad3 = new TPad("pad3","pad3",0.01,0.51,0.49,0.99);
  pad3->Draw();
  TPad* pad4 = new TPad("pad4","pad4",0.51,0.51,0.99,0.99);
  pad4->Draw();

  pad1->cd();
  grTAU->Draw("ap");
  pad2->cd();
  grDM->Draw("ap");
  pad3->cd();
//  pad3->GetYaxis()->SetRangeUser(-1.,1.);
//  mgA->GetYaxis()->SetRangeUser(-1.,1.);
  mgA->Draw("ap");
  pad4->cd();
//  mgB->GetYaxis()->SetRangeUser(-1.,1.);
  mgB->Draw("ap");

  c2->Update();
  c2->Print("fit1.png");
  c2->Print("fit1.root");
  system("display fit1.png &");

  return 0;
}

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
