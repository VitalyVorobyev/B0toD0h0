#ifndef PSEUDOTOYFIT_H
#define PSEUDOTOYFIT_H

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

#include "RkRdetRnpPdf.h"

RkRdetRnpPdf* m_pdf;
TTree* m_tree;
TTree* m_good_tree;
Int_t m_exp,m_bin,m_flv;
Double_t m_dt;
Double_t m_A,m_B;
Double_t m_costhBcms;
const int m_svd = 2;

// * Event-by-event parameters
int ntrk_rec;
int ntrk_asc;
int ndf_rec;
int ndf_asc;
double sz_rec;
double sz_asc;
double chisq_rec;
double chisq_asc;
// *

const int EventsInASample = 1000;
int SampleNumber = 0;
int NSamples = 0;

using namespace std;
using namespace ROOT;
using namespace Minuit2;

const double sol = 2.99792458;
bool draw_plots = true;
const double cm2ps = 78.48566945838871754705;
const double beta = 23.;
const double m_btau     = 1.534;//1.520;//0.460/sol*cm2ps;//0.507;// ps
const double m_dm       = 0.459;//0.510/sol;//0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
const double m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
const double m_cos2beta = TMath::Cos(2.*beta/180.*TMath::Pi());
int xi;
const int NDots = 50;
const int NBins = 50;
const double dtmax = 10;
const double dtmin = -dtmax;
const double dzmax =  dtmax/78.4857;
const double dzmin = -dzmax;

const bool no_interf     = true;
const bool no_np         = false;
const bool fix_btau      = false;
const bool fix_dm        = false || no_interf;
const bool fix_sin2beta  = false || no_interf;
const bool fix_cos2beta  = false || no_interf;

double Carr[8], Sarr[8], Karr[8], Kbarr[8];
// CLEO measurements
const double Carr_CLEO[8] = { 0.365, 0.710, 0.481,0.008,-0.757,-0.884,-0.462, 0.106};
const double Sarr_CLEO[8] = {-0.179,-0.013,-0.147,0.938, 0.386,-0.162,-0.616,-1.063};
const double Karr_CLEO[8] = { 0.067, 0.093, 0.030,0.089, 0.079, 0.102, 0.123, 0.159};
const double Kbarr_CLEO[8]= { 0.016, 0.027, 0.013,0.040, 0.015, 0.012, 0.026, 0.109};

//Model integrals
//const double Carr_model[8] = {0.594511, -0.329868,-0.61714, -0.758994,-0.368921,  0.133975, 0.484073, 0.742303};
//const double Sarr_model[8] = {-0.357055,-0.680259,-0.277483,-0.259506, 0.440863,  0.623397, 0.358263,-0.0558073};
//const double Karr_model[8] = {0.0665352,0.0835056, 0.0360653,0.0961057,0.0722859, 0.0980402,0.123951, 0.198887};
//const double Kbarr_model[8]= {0.0304473,0.00571654,0.0022683,0.0264409,0.00969127,0.0105694,0.0311109,0.10838};

const double Carr_model[8] = {0.554166,-0.00454239,-0.608498,-0.934161,-0.565509,0.0800234,0.480612,0.68653};
const double Sarr_model[8] = {0.446559,0.825851,0.678567,-0.000976537,-0.572031,-0.734359,-0.372541,0.0279944};
const double Karr_model[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
const double Kbarr_model[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};

// rho Ks
const double Carr_rhoKs[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
const double Sarr_rhoKs[8] = {0,0,0,0,0,0,0,0};
const double Karr_rhoKs[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_rhoKs[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

// no interfer
const double Carr_noint[8] = {0,0,0,0,0,0,0,0};
const double Sarr_noint[8] = {0,0,0,0,0,0,0,0};
const double Karr_noint[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_noint[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

inline double C(const int bin){ return Carr[abs(bin)-1];}
inline double S(const int bin){ return bin>0 ? Sarr[abs(bin)-1] : -Sarr[abs(bin)-1];}
inline double K(const int bin){ return bin>0 ? Karr[abs(bin)-1] : Kbarr[abs(bin)-1];}
inline double A(const int flv,const int bin){ return flv*(K(bin)-K(-bin))/(K(bin)+K(-bin));}
inline double B(const int flv,const int bin){ return 2.*flv*xi*sqrt(K(bin)*K(-bin))/(K(bin)+K(-bin))*(C(bin)*m_sin2beta+S(bin)*m_cos2beta);}
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

int init(const int _mode){
  cout << "init: mode " << _mode << endl;
  m_pdf->SetRange(dtmax);
  if(no_interf){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_noint[i];
      Sarr[i] = Sarr_noint[i];
      Karr[i] = Karr_noint[i];
      Kbarr[i] = Kbarr_noint[i];
    }
  } else if(_mode == 5){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_rhoKs[i];
      Sarr[i] = Sarr_rhoKs[i];
      Karr[i] = Karr_rhoKs[i];
      Kbarr[i] = Kbarr_rhoKs[i];
    }
  } else if(_mode == 4){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_model[i];
      Sarr[i] = Sarr_model[i];
      Karr[i] = Karr_model[i];
      Kbarr[i] = Kbarr_model[i];
    }
  } else{
    for(int i=0; i<8; i++){
      Carr[i] = Carr_CLEO[i];
      Sarr[i] = Sarr_CLEO[i];
      Karr[i] = Karr_CLEO[i];
      Kbarr[i] = Kbarr_CLEO[i];
    }
  }
  cout << "Initialized arrays:" << endl;
  for(int i=0;i<8;i++){
    cout << i+1 << " " << Karr[i] << " " << Kbarr[i] << " " << Carr[i] << " " << Sarr[i] << endl;
  }
  return 0;
}

void GetEvent(const int i){
  m_good_tree->GetEvent(i);
  m_dt *= 0.1*cm2ps;
  sz_rec /= 10.;
  sz_asc /= 10.;
  return;
}

void SetPDFParams(const vector<double>& par){
  m_pdf->SetTauDm(par[0],par[1]);
  m_pdf->SetSinCos(par[2],par[3]);
  return;
}

TTree* GetGoodTTree(TTree* tree, const int _mode, const int bin = 0, const int flv = 0){
  double demin, demax, mbcmin, mbcmax, chisq_sig_max, chisq_asc_max, sz_sig_max, sz_asc_max;
  double bdtg_cut;
  int mode, h0mode;
  mbcmin = 5.272;
  mbcmax = 5.286;
  chisq_sig_max = 10;
  chisq_asc_max = 10;
  sz_sig_max = 0.1;
  sz_asc_max = 0.2;
  switch(_mode){
  case 1:// D0 pi0
    demin = -0.08;
    demax =  0.07;
    mode = 1;
    h0mode = 10;
    bdtg_cut = 0.910;
    break;
  case 2:// D0 eta->gg
    demin = -0.08;
    demax =  0.07;
    mode = 2;
    h0mode = 10;
    bdtg_cut = 0.690;
    break;
  case 3:// D0 eta->pi+pi-pi0
    demin = -0.05;
    demax =  0.04;
    mode = 2;
    h0mode = 20;
    bdtg_cut = 0.500;
    break;
  case 4:// D0 omega->pi+pi-pi0
    demin = -0.06;
    demax =  0.05;
    mode = 3;
    h0mode = 20;
    bdtg_cut = 0.125;
  case 5:// D0 omega->pi+pi-pi0 (Ks rho model)
    demin = -0.06;
    demax =  0.05;
    mode = 3;
    h0mode = 20;
    bdtg_cut = 0.125;
    break;
  }
  stringstream out;
  out.str("");
  if(m_svd == 2) out << "exp>30 && ";
  else         out << "exp<30 && ";
  if(bin)      out << "bin_mc == " << bin << " && ";
  if(flv)      out << "flv_mc == " << flv << " && ";
  out << "de<" << demax << " && " << "de>" << demin << " && ";
  out << "mbc<" << mbcmax << " && " << "mbc>" << mbcmin << " && ";
  out << "chisq_z_sig>0 && chisq_z_sig<" << chisq_sig_max << " && ";
  out << "chisq_z_asc>0 && chisq_z_asc<" << chisq_asc_max << " && ";
  out << "sz_sig>0 && sz_sig<" << sz_sig_max << " && ";
  out << "sz_asc>0 && sz_asc<" << sz_asc_max << " && ";
  out << "bdtg>" << bdtg_cut << " && ";
//  out << "ndf_z_asc>0 && ";
  if(no_np) out << "!nptag && ";
  out << "dz*0.1*78.48>" << dtmin << " && dz*0.1*78.48<" << dtmax;
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
    double norm = 0;

    SetPDFParams(par);
    for(int i=SampleNumber*EventsInASample; i<(SampleNumber+1)*EventsInASample && i<NTot; i++){
      GetEvent(i);
      m_pdf->SetAkCk(m_costhBcms,0.5*10.58);
      m_pdf->SetKKCS(K(m_bin),K(-m_bin),C(m_bin),S(m_bin));
      m_pdf->SetFlvXi(m_flv,xi);
      if(!no_np){
        sigpdf = m_pdf->Pdf(m_dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,no_interf);
      } else{
        sigpdf = m_pdf->NoNPPdf(m_dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,no_interf);
      }
      if(!std::isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf);
      else{ cout << "pdfFcn: pdf = " << sigpdf << endl;}
    }
    cout << "loglh: " << loglh << ", norm: " << norm/32. << ", tau: " << par[0] << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3];
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
    double pdf = m_pdf->PdfAB(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
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
    const double& btau = par[0];
    const double& dm   = par[1];
    const double& A    = par[2];
    const double& B    = par[3];
    double loglh = 0;
    double sigpdf;
    m_pdf->SetTauDm(btau,dm);
    m_pdf->SetAB(A,B);
    double sum_norm = 0;
    for(int i=SampleNumber*EventsInASample; i<(SampleNumber+1)*EventsInASample && i<NTot; i++){
      GetEvent(i);
      m_pdf->SetAkCk(m_costhBcms,0.5*10.580);
      const double norm = 1;//calc_norm_single(ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,false);
      sigpdf = m_pdf->PdfAB(m_dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
      sum_norm += norm/EventsInASample;
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

int draw_fit_results2(const vector< vector<double> >& vals,const vector< vector<double> >& errs){
  const int NFits = vals.size();

  TH1I* tau_hist = new TH1I("tau_hist","tau_hist",50,1.4,1.6);
  TH1I* dm_hist = new TH1I("dm_hist","dm_hist",50,0.2,0.7);
  TH1I* cos_hist = new TH1I("cos_hist","cos_hist",50,0.2,1.);
  TH1I* sin_hist = new TH1I("sin_hist","sin_hist",50,0.2,1.);

  TH1I* pull_tau_hist = new TH1I("pull_tau_hist","pull_tau_hist",50,-5.,5.);
  TH1I* pull_dm_hist = new TH1I("pull_dm_hist","pull_dm_hist",50,-5.,5.);
  TH1I* pull_cos_hist = new TH1I("pull_cos_hist","pull_cos_hist",50,-5.,5.);
  TH1I* pull_sin_hist = new TH1I("pull_sin_hist","pull_sin_hist",50,-5.,5.);

  for(int i=0; i<NFits; i++){
    tau_hist->Fill(vals[i][0]); pull_tau_hist->Fill(vals[i][0]-m_btau/errs[i][0]);
    dm_hist->Fill(vals[i][1]);  pull_dm_hist->Fill(vals[i][1]-m_dm/errs[i][1]);
    cos_hist->Fill(vals[i][3]); pull_cos_hist->Fill(vals[i][3]-m_cos2beta/errs[i][3]);
    sin_hist->Fill(vals[i][2]); pull_sin_hist->Fill(vals[i][2]-m_sin2beta/errs[i][2]);
  }

  TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
  c1->cd();
  c1->Draw();
  TPad* pad1a = new TPad("pad1a","pad1a",0.01,0.01,0.49,0.49);
  pad1a->Draw();
  TPad* pad2a = new TPad("pad2a","pad2a",0.51,0.01,0.99,0.49);
  pad2a->Draw();
  TPad* pad3a = new TPad("pad3a","pad3a",0.01,0.51,0.49,0.99);
  pad3a->Draw();
  TPad* pad4a = new TPad("pad4a","pad4a",0.51,0.51,0.99,0.99);
  pad4a->Draw();

  pad1a->cd();
  pull_tau_hist->Draw("e");
  pad2a->cd();
  pull_dm_hist->Draw("e");
  pad3a->cd();
  pull_cos_hist->Draw("e");
  pad4a->cd();
  pull_sin_hist->Draw("e");

  c1->Update();
  c1->Print("ptoy_pulls.png");
  c1->Print("ptoy_pulls.root");
  system("display ptoy_pulls.png &");


  TCanvas* c2 = new TCanvas("c2","c2",1200,1200);
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
  tau_hist->Draw("e");
  pad2->cd();
  dm_hist->Draw("e");
  pad3->cd();
  cos_hist->Draw("e");
  pad4->cd();
  sin_hist->Draw("e");

  c2->Update();
  c2->Print("ptoy.png");
  c2->Print("ptoy.root");
  system("display ptoy.png &");

  return 0;
}

#endif // PSEUDOTOYFIT_H
