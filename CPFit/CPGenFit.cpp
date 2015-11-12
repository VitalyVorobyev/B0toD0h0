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

#include "Rk.h"

TFile *ifile;
TTree* m_tree;
TTree* m_raw_tree;
Int_t m_exp,m_bin,m_flv,m_good_icpv;
Double_t m_dz;
Double_t m_A,m_B;
Double_t m_costhBcms;

using namespace std;
using namespace ROOT;
using namespace Minuit2;

const bool RkFit = true;

const double sol = 2.99792458;
bool draw_plots = true;
const double cm2ps = 78.4857;
const double beta = 67.;//67.;
const double m_btau     = 1.534;//(RkFit ? 0.0198462 : 0.0459);//0.0454193;
const double m_dm       = 0.510;//0.452;//(RkFit ? 40 : 17.);//0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
const double m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
const double m_cos2beta = TMath::Cos(2.*beta/180.*TMath::Pi());
int xi;
const int NDots = 100;
const int NBins = 50;
//const double dtmax = (RkFit ? 10 : 30);
//const double dtmin = -dtmax;
const double dzmax = 10.;// dtmax/78.4857;
const double dzmin = -10.;//-dzmax;

double btau = m_btau;
double dm = m_dm;
double sin2beta = m_sin2beta;
double cos2beta = m_cos2beta;

const bool make_k_fit   = false;
const bool no_interf    = false;
const bool draw_b_bbar  = false;

const bool fix_btau     = true;
const bool fix_dm       = true  || no_interf;
const bool fix_sin2beta = false || no_interf;
const bool fix_cos2beta = false || no_interf;

bool is_good_event(){
//  if(flv != -1) return false;
//  if(abs(bin) == 1 || abs(bin) == 8) return false;
  if(abs(m_dz)>dzmax)        return false;
  if(!m_bin || abs(m_bin)>8) return false;
  return true;
}

double Carr[8], Sarr[8], Karr[8], Kbarr[8];
// CLEO measurements
const double Carr_CLEO[8] = { 0.365, 0.710, 0.481,0.008,-0.757,-0.884,-0.462, 0.106};
const double Sarr_CLEO[8] = {-0.179,-0.013,-0.147,0.938, 0.386,-0.162,-0.616,-1.063};
const double Karr_CLEO[8] = {0.067,0.093,0.030,0.089,0.079,0.102,0.123,0.159};
const double Kbarr_CLEO[8]= {0.016,0.027,0.013,0.040,0.015,0.012,0.026,0.109};

//Model integrals
//const double Carr_model[8] = {0.594511, -0.329868,-0.61714, -0.758994,-0.368921,  0.133975, 0.484073, 0.742303};
//const double Sarr_model[8] = {-0.357055,-0.680259,-0.277483,-0.259506, 0.440863,  0.623397, 0.358263,-0.0558073};
//const double Karr_model[8] = {0.0665352,0.0835056, 0.0360653,0.0961057,0.0722859, 0.0980402,0.123951, 0.198887};
//const double Kbarr_model[8]= {0.0304473,0.00571654,0.0022683,0.0264409,0.00969127,0.0105694,0.0311109,0.10838};

const double Carr_model[8] = {0.554166,-0.00454239,-0.608498,-0.934161,-0.565509,0.0800234,0.480612,0.68653};
const double Sarr_model[8] = {0.446559,0.825851,0.678567,-0.000976537,-0.572031,-0.734359,-0.372541,0.0279944};
const double Karr_model[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
const double Kbarr_model[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};

//const double Karr_fit[8] = {0.0618,};
//const double Kbarr_fit[8]= {};

// rho Ks
const double Carr_rhoKs[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
const double Sarr_rhoKs[8] = {0,0,0,0,0,0,0,0};
const double Karr_rhoKs[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_rhoKs[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

inline double C(const int bin){ return Carr[abs(bin)-1];}
inline double S(const int bin){ return bin>0 ? Sarr[abs(bin)-1] : -Sarr[abs(bin)-1];}
inline double K(const int bin){ return bin>0 ? Karr[abs(bin)-1] : Kbarr[abs(bin)-1];}
inline double A(const int flv,const int bin){ return flv*(K(bin)-K(-bin))/(K(bin)+K(-bin));}
inline double B(const int flv,const int bin){ return 2.*flv*xi*sqrt(K(bin)*K(-bin))/(K(bin)+K(-bin))*(C(bin)*sin2beta+S(bin)*cos2beta);}
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

void GetEvent(const int i){
  m_tree->GetEvent(i);
  RkFit ? m_dz *= 0.1*cm2ps : m_dz *= 100/sol;
  return;
}

int SetGoodTTree(const int _mode, const int bin = 0, const int flv = 0){
  double demin, demax, mbcmin, mbcmax;
  delete m_tree;
//          chisq_sig_max, chisq_asc_max, sz_sig_max, sz_asc_max;
  double bdtg_cut;
  int mode, h0mode;
  mbcmin = 5.272;
  mbcmax = 5.286;
//  chisq_sig_max = 10;
//  chisq_asc_max = 10;
//  sz_sig_max = 0.1;
//  sz_asc_max = 0.2;
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
//  if(m_svd == 2) out << "exp>30 && ";
//  else         out << "exp<30 && ";
  if(bin)      out << "bin_mc == " << bin << " && ";
  if(flv)      out << "flv_mc == " << flv << " && ";
  out << "de<" << demax << " && " << "de>" << demin << " && ";
  out << "mbc<" << mbcmax << " && " << "mbc>" << mbcmin << " && ";
//  out << "chisq_z_sig>0 && chisq_z_sig<" << chisq_sig_max << " && ";
//  out << "chisq_z_asc>0 && chisq_z_asc<" << chisq_asc_max << " && ";
//  out << "sz_sig>0 && sz_sig<" << sz_sig_max << " && ";
//  out << "sz_asc>0 && sz_asc<" << sz_asc_max << " && ";
  out << "good_icpv == 1 && ";
  out << "bdtg>" << bdtg_cut << " && ";
//  out << "ndf_z_asc>0 && ";
//  if(no_np) out << "!nptag && ";
  out << "dz_mc*0.1*78.48>" << dzmin << " && dz_mc*0.1*78.48<" << dzmax;
  m_tree = m_raw_tree->CopyTree(out.str().c_str());
  return m_tree->GetEntries();
}

void calc_K(const vector<double>& A,const vector<double>& B){
  if(A.size() != 16 || B.size() != 16){
    cout << "calc_K: wrong vectors size: " << A.size() << ", " << B.size() << endl;
    return;
  }
  cout << "sin2beta = " << sin2beta << ", cos2beta = " << cos2beta << endl;
  double Norm1 = 0, Norm2 = 0;
  double K[16], Kb[16];
  cout << "calc K:" << endl;
  for(int i=0; i<16; i++){
    const double a = A[i]*(C(bin(i))*sin2beta+S(bin(i))*cos2beta);// K - Kb
    const double b = B[i];
    const double c = sqrt(a*a + b*b);// K + Kb
    K[i]  = 0.5*(c+a);
    Kb[i] = 0.5*(c-a);
    if(i<8) Norm1 += K[i]+Kb[i];
    else    Norm2 += K[i]+Kb[i];
//    cout << "a = " << a << ", b = " << b << ", c = " << c << ", C = " << C(bin(i)) << ", S = " << S(bin(i)) << endl;
  }
  cout << " Norm1 = " << Norm1 << ", Norm2 = " << Norm2 << endl;
  for(int i=0; i<8; i++){
    K[i]   /= Norm1; Kb[i]   /= Norm1;
    K[i+8] /= Norm2; Kb[i+8] /= Norm2;
  }
  for(int i=0; i<16; i++){
    cout << "  K[" << bin(i) << "] = " << K[i] << ", Kb[" << bin(i) << "] = " << Kb[i] << endl;
  }
  return;
}


int init_arrs(const int mode){
  if(mode == 5){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_rhoKs[i];
      Sarr[i] = Sarr_rhoKs[i];
      Karr[i] = Karr_rhoKs[i];
      Kbarr[i] = Kbarr_rhoKs[i];
    }
  } else if(mode == 4 || mode == 3 || mode == 2){
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

double pdfSig(void){
  const double dt = m_dz;
  if(no_interf){
    return TMath::Exp(-TMath::Abs(dt)/btau)/(2.*btau);
  }
  const int flv = m_flv;
  const int bin = m_bin;
  double pdf = K(bin)+K(-bin) - flv*(K(bin)-K(-bin))*cos(dm*dt) + 2.*flv*xi*sqrt(K(bin)*K(-bin))*sin(dm*dt)*(C(bin)*sin2beta+S(bin)*cos2beta);
  double pdf_norm = K(bin)+K(-bin) - flv*(K(bin)-K(-bin))/(dm*dm*btau*btau+1);
  if(pdf<0 || pdf_norm<=0) return 0;
  return TMath::Exp(-TMath::Abs(dt)/btau)*pdf/(2.*btau*pdf_norm);
}

double pdfSigSingle(void){
  const double dt = m_dz;
  const double A  = m_A;
  const double B  = m_B;
  const double pdf = 1. - A*cos(dm*dt) + B*sin(dm*dt);
  const double pdf_norm = 1. - A/(dm*dm*btau*btau+1);
  if(pdf<0 || pdf_norm<=0) return 0;
  return TMath::Exp(-TMath::Abs(dt)/btau)*pdf/(2.*btau*pdf_norm);
}

double calc_norm_single(void){
  double norm = 0;
  const double ddz = (dzmax - dzmin)/100.;
  for(int i=0; i<100; i++){
    m_dz = dzmin + (i+0.5)*ddz;
    double pdf = pdfSigSingle();
    if(pdf>0) norm += pdf;
  }
  return norm*ddz;
}

double calc_norms(double norms[2][16]){
  double norm = 0;
  const double ddz = (dzmax - dzmin)/100.;
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      norms[k][j] = 0;
    }
  }
  for(int k=0; k<2; k++){
    m_flv = flv(k);
    for(int j=0; j<16; j++){
      m_bin = bin(j);
      for(int i=0; i<100; i++){
        m_dz = dzmin + (i+0.5)*ddz;
        double pdf = pdfSig();
        if(pdf>0){
          norms[k][j] += pdf;
          norm        += pdf;
        }
      }
    }
  }
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      norms[k][j] *= ddz;
    }
  }
  return norm*ddz;
}

class pdfFcn : public FCNBase{
public:
  pdfFcn(void){
    theErrorDef = 1;
    NTot = m_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcn() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    btau     = par[0];
    dm       = par[1];
    sin2beta = par[2];
    cos2beta = par[3];
    double loglh = 0;
    double sigpdf;
    double norm = 0;

    if(make_k_fit){
      double Ksum = 0;
      for(int i=0; i<8; i++){
        Karr[i] = par[4+2*i]; Ksum += Karr[i];
        if(i != 7){ Kbarr[i] = par[5+2*i]; Ksum += Kbarr[i];}
      }
      Kbarr[7] = 1. - Ksum;
      cout << "Kbarr[7] = " << Kbarr[7] << endl;
    }

    if(RkFit){
      RkPdf rkpdf;
      rkpdf.SetRange(dzmax);
      rkpdf.SetTauDm(btau,dm);
      rkpdf.SetSinCos(sin2beta,cos2beta);
      for(int i=0; i<NTot; i++){
        GetEvent(i);
        if(!is_good_event()) continue;
        rkpdf.SetAkCk(m_costhBcms,0.5*10.58);
        rkpdf.SetKKCS(K(m_bin),K(-m_bin),C(m_bin),S(m_bin));
        rkpdf.SetFlvXi(m_flv,xi);
        sigpdf = rkpdf.Pdf(m_dz,no_interf);
        if(!isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf);
      }
    } else{
      double norms[2][16];
      norm = calc_norms(norms);
      for(int i=0; i<NTot; i++){
        GetEvent(i);
        if(!is_good_event()) continue;
        sigpdf = pdfSig();
        if(!isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf/norms[flv_ind(m_flv)][bin_ind(m_bin)]);
      }
    }
    cout << "loglh: " << loglh << ", norm: " << norm/32. << ", tau: " << par[0] << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3] << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

class pdfFcnSingle : public FCNBase{
public:
  pdfFcnSingle(vector<double>* _vec, vector<double>* _cosBvec){
    theErrorDef = 1;
    vec = _vec;
    cosBvec = _cosBvec;
    NTot = vec->size();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnSingle() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    btau = par[0];
    dm   = par[1];
    m_A  = par[2];
    m_B  = par[3];
    double loglh = 0;
    double sigpdf;
    const double norm = calc_norm_single();
    if(RkFit){
      RkPdf rkpdf;
      rkpdf.SetTauDm(btau,dm);
      rkpdf.SetAB(m_A,m_B);
      rkpdf.SetRange(dzmax);
      for(int i=0; i<NTot; i++){
        m_dz = vec->at(i);
        if(abs(m_dz)>dzmax) continue;
        rkpdf.SetAkCk(cosBvec->at(i),0.5*10.580);
        sigpdf = rkpdf.PdfAB(m_dz,no_interf);
        if(!isnan(sigpdf) && sigpdf>0){
            loglh += -2*TMath::Log(sigpdf);
//            cout << sigpdf << endl;
        }
        else{ cout << "pdfFcnSingle: pdf = " << sigpdf << endl;}
      }
    } else{
      for(int i=0; i<NTot; i++){
        m_dz = vec->at(i);
        if(abs(m_dz)>dzmax) continue;
        sigpdf = pdfSigSingle();
        if(!isnan(sigpdf) && sigpdf>0) loglh += -2*TMath::Log(sigpdf/norm);
      }
    }
    cout << "loglh: " << loglh << ", norm = " << norm << ", tau: " << par[0] << ", dm = " << par[1] << ", A: " << par[2] << ", B: " << par[3] << endl;
    return loglh;
  }
private:
  double theErrorDef;
//  double A,B;
  vector<double>* vec;
  vector<double>* cosBvec;
  int NTot;
};

double calc_chisq(const int n,TH1I* hist,const TGraph* gr){
  double chisq = 0;
  double* xarr = gr->GetX();
  double* yarr = gr->GetY();
  for(int i=0; i<n; i++){
    double val = hist->GetBinContent(hist->FindBin(xarr[i]));
    if(val>0)chisq += (val-yarr[i])*(val-yarr[i])/TMath::Sqrt(val);
  }
  return chisq/(n-1);
}

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

  TCanvas* c2 = new TCanvas("c2","c2",800,400);
  c2->cd();
  c2->Draw();
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.99);
  pad1->SetGrid();
  pad1->Draw();
  TPad* pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.99);
  pad2->SetGrid();
  pad2->Draw();

  pad1->cd();
  grTAU->Draw("ap");
  pad2->cd();
  grDM->Draw("ap");

  c2->Update();
  c2->Print("fit1.png");
  c2->Print("fit1.root");

  system("display fit1.png &");

  TCanvas* c1 = new TCanvas("c1","c1",800,400);
  c1->cd();
  c1->Draw();
  TPad* pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.99);
  pad3->Draw();
  TPad* pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.99);
  pad4->Draw();

  pad3->cd();
  pad3->SetGrid();
//  pad3->GetYaxis()->SetRangeUser(-1.,1.);
//  mgA->GetYaxis()->SetRangeUser(-1.,1.);
  mgA->Draw("ap");
  pad4->cd();
  pad4->SetGrid();
//  mgB->GetYaxis()->SetRangeUser(-1.,1.);
  mgB->Draw("ap");

  c1->Update();
  c1->Print("fitAB.png");
  c1->Print("fitAB.root");

  system("display fitAB.png &");
  return 0;
}

MnUserParameterState make_single_fit(vector<double>& vec,vector<double>& cosBvec){
  btau     = m_btau;
  dm       = m_dm;
  int NPar = 0;
  MnUserParameters upar;
  upar.Add(string("btau"),btau,0.2*btau,0,10*btau); NPar++;
  upar.Add(string("dm"),dm,0.1*dm,0.,5.*dm);        NPar++;
  upar.Add(string("A"),m_A,0.1,-3.,3.);             NPar++;
  upar.Add(string("B"),m_B,0.1,-3.,3.);             NPar++;

  cout << "Single fit: " << vec.size() << "=" << cosBvec.size() << " events, A0 = " << m_A << ", B0 = " << m_B << endl;

  pdfFcnSingle* theFCN = new pdfFcnSingle(&vec,&cosBvec);
  MnMigrad migrad(*theFCN,upar);
  migrad.Fix("btau");
  migrad.Fix("dm");
//  migrad.Fix("A");
//  migrad.Fix("B");
  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    cout << pstate.Value(i) << " +- " << pstate.Error(i);
    cout << " (" << (pstate.Value(i)-upar.Value(i))/pstate.Error(i) << ")" << endl;
    upar.SetValue(i,pstate.Value(i));
  }

  if(draw_plots){
    btau = upar.Value(0);
    dm   = upar.Value(1);
    m_A  = upar.Value(2);
    m_B  = upar.Value(3);

    TH1I dh("dh","dh",NBins,dzmin,dzmax);
    dh.SetMarkerStyle(20);
    dh.SetMarkerColor(kBlue);
    dh.SetMarkerSize(1.2);
    const int NTot = vec.size();
    for(int i=0; i<NTot; i++){
      if(abs(vec[i])<dzmax) dh.Fill(vec[i]);
    }
    double pdf_arr[NDots],dz_arr[NDots];
    double norm = 0;

    const double ddz = (dzmax-dzmin)/(double)NDots;
    const double Nddz = ddz*NTot*NDots/NBins;

    if(RkFit){
      double pdf = 0;
      RkPdf rkpdf;
      rkpdf.SetRange(dzmax);
      rkpdf.SetTauDm(btau,dm);
      rkpdf.SetAB(m_A,m_B);
      for(int i=0; i<NDots; i++){
        pdf_arr[i] = 0;
        dz_arr[i] = dzmin+(i+0.5)*ddz;
        for(int j=0; j<vec.size(); j++){
          rkpdf.SetAkCk(cosBvec[j],0.5*10.58);
          pdf = rkpdf.PdfAB(dz_arr[i],no_interf);
          pdf_arr[i] += pdf;
          norm       += pdf;
        }
        pdf_arr[i] *= Nddz/vec.size();
        }
        norm *= Nddz/vec.size();
      } else{
        for(int i=0; i<NDots; i++){
        dz_arr[i]  = dzmin+(i+0.5)*ddz; m_dz = dz_arr[i];
        pdf_arr[i] = pdfSigSingle()*Nddz;
        norm += pdf_arr[i];
      }
    }

    cout << "norm = " << norm << endl;

    TGraph* gr = new TGraph(NDots,dz_arr,pdf_arr);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.);
    gr->SetLineWidth(2.);
    gr->SetMarkerColor(kRed);

    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    c1->SetLogy();
    dh.Draw("e");
    gr->Draw("same");

    c1->Print("fit.png");
    system("display fit.png &");
  }

  return pstate;
}

int make_full_fit(void){
  stringstream out;
  int NPar = 0;
  MnUserParameters upar;
  cout << " Init parameters:" << endl;
  upar.Add(string("btau"),btau,0.2*btau,0,10*btau); NPar++;
  upar.Add(string("dm"),dm,0.1*dm,0.7*dm,1.3*dm);   NPar++;
  upar.Add(string("sin2beta"),sin2beta,0.1,-1.,1.); NPar++;
  upar.Add(string("cos2beta"),cos2beta,0.1,-1.,1.); NPar++;
  if(make_k_fit){
    for(int i=0; i<8; i++){
      out.str("");
      out << "K" << i+1;
      upar.Add(out.str(),K(i+1),0.1,0.,1.); NPar++;
      if(i != 7){
        out.str("");
        out << "Kb" << i+1;
        upar.Add(out.str(),K(-(i+1)),0.1,0.,1.); NPar++;
      }
    }
  }

  cout << " Start MIGRAD" << endl;
  pdfFcn* theFCN = new pdfFcn();
  MnMigrad migrad(*theFCN,upar);

  if(fix_btau                 ) migrad.Fix("btau");
  if(fix_dm       || no_interf) migrad.Fix("dm");
  if(fix_sin2beta || no_interf) migrad.Fix("sin2beta");
  if(fix_cos2beta || no_interf) migrad.Fix("cos2beta");

  cout << "Starting minimization" << endl;
  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  cout << "after migrad" << endl;

  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << pstate.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  const int NTot = m_tree->GetEntries();
  btau     = upar.Value(0);
  dm       = upar.Value(1);
  sin2beta = upar.Value(2);
  cos2beta = upar.Value(3);

  double dz_arr[NDots];
  const double ddz = (dzmax-dzmin)/(double)NDots;
  cout << "Init hists and fill graphs" << endl;
  for(int j=0; j<NDots; j++){
    dz_arr[j] = dzmin+(j+0.5)*ddz;
    cout << dz_arr[j] << " ";
  }
  cout << endl;

//  if(make_k_fit){
//    double Ksum = 0;
//    for(int i=0; i<8; i++){
//      Karr[i]  = upar.Value(4+2*i); Ksum += Karr[i];
//      if(i != 7){ Kbarr[i] = upar.Value(5+2*i); Ksum += Kbarr[i];}
//    }
//    Kbarr[7] = 1. - Ksum;
//    cout << "Kb[8] fit value = " << Kbarr[7] << endl;
//  }

//  if(RkFit){


    RkPdf rkpdf;
    rkpdf.SetRange(dzmax);
    rkpdf.SetTauDm(btau,dm);
    rkpdf.SetSinCos(sin2beta,cos2beta);

    if(no_interf){
      rkpdf.SetKKCS(K(1),K(1),C(1),S(1));
      int Nev = 0;
      TH1I* dh = new TH1I("dh","dh",NBins,dzmin,dzmax);
      dh->SetMarkerStyle(20);
      dh->SetMarkerColor(kBlue);
      dh->SetMarkerSize(1.2);
      double pdf_arr[NDots];
      double norm = 0, norm1 = 0;
      for(int j=0; j<NDots; j++){
        pdf_arr[j] = 0;
        for(int l=0; l<NTot; l++){
          GetEvent(l);
          rkpdf.SetFlvXi(m_flv,xi);
          rkpdf.SetAkCk(m_costhBcms,0.5*10.58);
          pdf_arr[j] += rkpdf.Pdf(dz_arr[j],no_interf);
        }
        cout << pdf_arr[j] << " ";
//        pdf_arr[j] *= ddz/NTot;
        norm += pdf_arr[j];
      }
      cout << endl;
      cout << "Norm = " << norm << endl;
      for(int i=0; i<NDots; i++) pdf_arr[i] /= norm;

      for(int i=0; i<NTot; i++){
        GetEvent(i);
        if(!is_good_event()) continue;
        dh->Fill(m_dz);
        Nev++;
      }
      cout << "Nev = " << Nev << endl;

      cout << "Set norm" << endl;
      for(int j=0; j<NDots; j++){
        pdf_arr[j]  *= Nev*NDots/(double)NBins;
        norm1       += pdf_arr[j];
      }
      cout << "Norm: " << norm1 << endl;

      TGraph* gr = new TGraph(NDots,dz_arr,pdf_arr);
      gr->SetMarkerStyle(kDot);
      gr->SetMarkerSize(1.5);
      gr->SetLineWidth(2);

      TCanvas* c1 = new TCanvas("c1","c1",800,600);
      c1->cd();
      dh->Draw("e");
      gr->Draw("same");

      c1->Print("genfit_no_int.png");
      c1->Print("genfit_no_int.root");

      system("display genfit_no_int.png &");
      return 0;
    }

    if(draw_b_bbar){
        rkpdf.SetKKCS(K(1),K(1),C(1),S(1));
        int Nev[2];
        TH1I* dh[2];
        for(int i=0; i<2; i++){
          dh[i] = !i ? new TH1I("dh","dh",NBins,dzmin,dzmax) : new TH1I("dhb","dhb",NBins,dzmin,dzmax);
          Nev[i] = 0;
          dh[i]->SetMarkerStyle(21);
          !i ? dh[i]->SetMarkerColor(kBlue) : dh[i]->SetMarkerColor(kRed);
          dh[i]->SetMarkerSize(1.1);
        }
        double pdf_arr[2][NDots];
        double norm[2], norm1[2];
        for(int k=0; k<2; k++){
          norm[k] = 0; norm1[k] = 0;
          rkpdf.SetFlvXi(flv(k),xi);
          for(int j=0; j<NDots; j++){
            pdf_arr[k][j] = 0;
            for(int l=0; l<NTot; l++){
              GetEvent(l);
              rkpdf.SetAkCk(m_costhBcms,0.5*10.58);
              pdf_arr[k][j] += rkpdf.Pdf(dz_arr[j],no_interf);
            }
            pdf_arr[k][j] *= ddz/NTot;
            norm[k] += pdf_arr[k][j];
          }
        }
        for(int k=0; k<2; k++){
          for(int i=0; i<NDots; i++){
              pdf_arr[k][i] /= norm[k];
          }
        }

        for(int i=0; i<NTot; i++){
          GetEvent(i);
          if(!is_good_event()) continue;
          dh[flv_ind(m_flv)]->Fill(m_dz);
          Nev[flv_ind(m_flv)]++;
        }
        cout << "Nev = " << Nev << endl;

        cout << "Set norm" << endl;
        for(int k=0; k<2; k++){
          for(int j=0; j<NDots; j++){
            pdf_arr[k][j]  *= Nev[k]*NDots/(double)NBins;
            norm1[k]       += pdf_arr[k][j];
          }
        }
//        cout << "Norm: " << norm1 << endl;

        TGraph* gr = new TGraph(NDots,dz_arr,pdf_arr[0]);
        gr->SetMarkerStyle(kDot);
        gr->SetMarkerSize(1.5);
        gr->SetLineWidth(2);

        TGraph* grb = new TGraph(NDots,dz_arr,pdf_arr[1]);
        grb->SetMarkerStyle(kDot);
        grb->SetMarkerSize(1.5);
        grb->SetLineWidth(2);

        TCanvas* c1 = new TCanvas("c1","c1",800,600);
        c1->cd();
        dh[0]->Draw("e");
        dh[1]->Draw("e,same");
        gr->Draw("same");
        grb->Draw("same");

        c1->Print("genfit_no_int.png");
        c1->Print("genfit_no_int.root");

        system("display genfit_no_int.png &");
        return 0;
    }

    int Nev[2][16];
    TH1I* dh[2][16];
    double pdf_arr[2][16][NDots];
    double norm[2][16];
    double norm1[2][16];
    for(int k=0; k<2; k++){
      m_flv = flv(k);
      cout << "flv = " << m_flv << endl;
      rkpdf.SetFlvXi(flv(k),xi);
      for(int i=0; i<16; i++){
        Nev[k][i]    = 0;
        norm[k][i]   = 0;
        norm1[k][i]  = 0;
        out.str("");
        out << "dh" << bin(i) << "_" << k;
        dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dzmin,dzmax);
        dh[k][i]->SetMarkerStyle(20);
        dh[k][i]->SetMarkerColor(kBlue);
        dh[k][i]->SetMarkerSize(1.1);
        m_bin = bin(i);
        cout << "  bin = " << m_bin << endl;
        rkpdf.SetKKCS(K(m_bin),K(-m_bin),C(m_bin),S(m_bin));
        for(int j=0; j<NDots; j++){
//          m_dz = dz_arr[j];
//          cout << "    dz = " << m_dz << endl;
          pdf_arr[k][i][j] = 0;
          for(int l=0; l<NTot; l++){
            GetEvent(l);
            rkpdf.SetAkCk(m_costhBcms,0.5*10.58);
            pdf_arr[k][i][j] += rkpdf.Pdf(dz_arr[j],no_interf);
//            cout << m_dz << " ";
          }
//          cout << endl;
          cout << pdf_arr[k][i][j] << " ";
          pdf_arr[k][i][j] *= ddz/NTot;
          norm[k][i] += pdf_arr[k][i][j];
        }
        cout << endl;
      }
    }

//    RkPdf rkpdf;
//    rkpdf.SetRange(dzmax);
//    rkpdf.SetTauDm(btau,dm);
//    rkpdf.SetSinCos(sin2beta,cos2beta);

//  } else{
//    for(int k=0; k<2; k++){
//      m_flv = flv(k);
//      for(int i=0; i<16; i++){
//        Nev[k][i]    = 0;
//        norm[k][i]   = 0;
//        norm1[k][i]  = 0;
//        out.str("");
//        out << "dh" << bin(i) << "_" << k;
//        dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dzmin,dzmax);
//        for(int j=0; j<NDots; j++){
//          m_dz = dz_arr[j];
//          m_bin = bin(i);
//          pdf_arr[k][i][j] = pdfSig()*ddz;
//          norm[k][i] += pdf_arr[k][i][j];
//        }
//      }
//    }
//  }

  int NEveCounter = 0;
  for(int i=0; i<NTot; i++){
//    m_tree->GetEvent(i);
    GetEvent(i);
    if(!is_good_event()) continue;
    const int k = flv_ind(m_flv);
    const int j = bin_ind(m_bin);
    dh[k][j]->Fill(m_dz);
    Nev[k][j]++;
    NEveCounter++;
  }
  cout << "NEveCounter = " << NEveCounter << endl;

  cout << "Set norm" << endl;
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      for(int j=0; j<NDots; j++){
        pdf_arr[k][i][j]  *= Nev[k][i]*NDots/(double)NBins;
        norm1[k][i]       += pdf_arr[k][i][j];
      }
    }
  }

  TGraph* gr[2][16];
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
//        cout << k << " " << i+1 << " " << norm[k][i] << " " << normb[k][i] << " ";
//        cout << k << " " << i+1 << " " << norm[k][i] << " " << norm[k][i+8] << " ";
//        cout << Nev[k][i] << " " << Nevb[k][i] << " ";
//        cout << Nev[k][i] << " " << Nev[k][i+8] << " ";
//        cout << dh[k][i]->GetEntries() << " " << dhb[k][i]->GetEntries() << " ";
//        cout << norm1[k][i] << " " << normb1[k][i] << endl;
      gr[k][i]  = new TGraph(NDots,dz_arr,pdf_arr[k][i]);
      gr[k][i]->SetMarkerStyle(kDot);
      gr[k][i]->SetMarkerSize(1.5);
      gr[k][i]->SetLineWidth(2);
      k == 0 ? gr[k][i]->SetMarkerColor(kRed) : gr[k][i]->SetMarkerColor(kRed);
    }
  }

//  double chisq[2][16];//,chisqb[2][8];
//  for(int k=0; k<2; k++){
//    for(int i=0; i<8; i++){
//      chisq[k][i] = calc_chisq(NBins,dh[k][i],gr[k][i]);
//    }
//  }

  cout << "Draw" << endl;
  TCanvas* c1 = new TCanvas("c1","c1",2600,800);
  TPad *pad[8];
  TPad *padb[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad" << i;
    pad[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.99);
    pad[i]->Draw();
    out.str("");
    out << "padb" << i;
    padb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.49);
    padb[i]->Draw();
  }
  int k = 0;
  for(int i=0; i<8; i++){
    pad[i]->cd();
    dh[k][i]->Draw("e");
    gr[k][i]->Draw("same");

    padb[i]->cd();
    dh[k][i+8]->Draw("e");
    gr[k][i+8]->Draw("same");
  }
  c1->Print("gen_fit.root");
  c1->Print("gen_fit.png");

  system("display gen_fit.png &");

  TCanvas* c2 = new TCanvas("c2","c2",2600,800);
  c2->cd();
  TPad *pad1[8];
  TPad *padb1[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad1" << i;
    pad1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.99);
    pad1[i]->Draw();
    out.str("");
    out << "padb1" << i;
    padb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.49);
    padb1[i]->Draw();
  }
  k = 1;
  for(int i=0; i<8; i++){
    pad1[i]->cd();
    dh[k][i]->Draw("e");
    gr[k][i]->Draw("same");

    padb1[i]->cd();
    dh[k][i+8]->Draw("e");
    gr[k][i+8]->Draw("same");
  }
  c2->Print("genb_fit.root");
  c2->Print("genb_fit.png");

  system("display genb_fit.png &");

//  for(int k=0; k<2; k++){
//    for(int i=0; i<8; i++){
//      cout << "chisq: " << k << " " << i+1 << " " << chisq[k][i] << " " << chisq[k][i+8] << endl;
//    }
//  }
  cout << "NEveCounter = " << NEveCounter << endl;
  return 0;
}

int main(const int argn, const char** argv){
  int _mode = 4;
  int _binfit = -1;
  int _flv = 1;
  if(argn == 2) sscanf(argv[1],"%d",&_mode);
  if(argn == 3){
    sscanf(argv[1],"%d",&_mode);
    sscanf(argv[2],"%d",&_binfit);
    if(_binfit) return -2;
  }
  if(argn == 4){
    sscanf(argv[1],"%d",&_mode);
    sscanf(argv[2],"%d",&_binfit);
    sscanf(argv[3],"%d",&_flv);
    if(abs(_flv) != 1) return -3;
  }
  switch(_mode){
  case 1:
    cout << "Mode: pi0" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s5_full.root");
    xi = 1;
    break;
  case 2:
    cout << "Mode: eta -> gg" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    xi = -1;
    break;
  case 3:
    cout << "Mode: eta -> pi+pi-pi0" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    xi = -1;
    break;
  case 4:
    cout << "Mode: omega" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s5_full.root");
    xi = -1;
    break;
  case 5:
    cout << "Mode: omega (rho Ks0)" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s3_full.root");
    xi = -1;
    break;
  default:
    cout << "Wrong mode" << endl;
    return -1;
  }
  init_arrs(_mode);
  m_raw_tree = (TTree*)ifile->Get("TEvent");
  m_raw_tree->SetBranchAddress("exp",   &m_exp);
  m_raw_tree->SetBranchAddress("bin_mc",&m_bin);
  m_raw_tree->SetBranchAddress("flv_mc",&m_flv);
  m_raw_tree->SetBranchAddress("good_icpv",&m_good_icpv);
  if(!RkFit) m_raw_tree->SetBranchAddress("dt_mc", &m_dz);
  else{
    m_raw_tree->SetBranchAddress("costhBcms", &m_costhBcms);
    m_raw_tree->SetBranchAddress("dz_mc", &m_dz);
  }

  SetGoodTTree(_mode,0,0);

  const int NTot = m_tree->GetEntries();
  int j,k;

  vector< vector<double> > vals_vec, errs_vec;
  vector<double> Aref,Bref;

  if(argn == 4){
    if(abs(_binfit)>8) return -1;
    if(_binfit){
      vector<double> dzvec,cosBvec;
      for(int i=0; i<NTot; i++){
        GetEvent(i);
        if(m_bin == _binfit && m_flv == _flv){
          dzvec.push_back(m_dz);
          cosBvec.push_back(m_costhBcms);
        }
      }
      m_flv = _flv; m_bin = _binfit;
      Aref.push_back(A(m_flv,m_bin));
      Bref.push_back(B(m_flv,m_bin));
      m_A = Aref[0]; m_B = Bref[0];
      MnUserParameterState upar = make_single_fit(dzvec,cosBvec);
      vals_vec.push_back(upar.Params());
      errs_vec.push_back(upar.Errors());
    } else{
      vector<double> Afit_vec, Bfit_vec;
      draw_plots = false;
      vector<double> dzvec[16],cosBvec[16];
      for(int i=0; i<NTot; i++){
        GetEvent(i);
        if(m_flv == _flv){
          if(!m_bin) continue;
          j = bin_ind(m_bin);
          dzvec[j].push_back(m_dz);
          cosBvec[j].push_back(m_costhBcms);
        }
      }
      m_flv = _flv;
      for(int i=0; i<16; i++){
        m_bin = bin(i);
        Aref.push_back(A(m_flv,m_bin));
        Bref.push_back(B(m_flv,m_bin));
        m_A = Aref[i]; m_B = Bref[i];
        MnUserParameterState upar = make_single_fit(dzvec[i],cosBvec[i]);
        Afit_vec.push_back(upar.Value(2));
        Bfit_vec.push_back(upar.Value(3));
        vals_vec.push_back(upar.Params());
        errs_vec.push_back(upar.Errors());
      }
      calc_K(Afit_vec,Bfit_vec);
    }
  } else if(argn == 3){
    draw_plots = false;
    vector<double> dzvec[2][16],cosBvec[2][16];
    for(int i=0; i<NTot; i++){
//      m_tree->GetEvent(i);
      GetEvent(i);
      j = bin_ind(m_bin);
      k = flv_ind(m_flv);
      dzvec[k][j].push_back(m_dz);
      cosBvec[k][j].push_back(m_costhBcms);
//      cosBvec[k][j].push_back(0);
     }
     for(k=0; k<2; k++){
       for(j=0; j<16; j++){
         m_flv = flv(k);
         m_bin = bin(j);
         Aref.push_back(A(m_flv,m_bin));
         Bref.push_back(B(m_flv,m_bin));
         m_A = Aref[16*k+j]; m_B = Bref[16*k+j];
         MnUserParameterState upar = make_single_fit(dzvec[k][j],cosBvec[k][j]);
         vals_vec.push_back(upar.Params());
         errs_vec.push_back(upar.Errors());
       }
     }
  } else {
    cout << "Make full fit" << endl;
    make_full_fit();
    return 0;
  }

  ofstream ofile;
  ofile.open("fit.txt",std::ofstream::out);
  ofile << "argn = " << argn << ", bin = " << _binfit << ", flv = " << _flv << endl;
  for(int i=0; i<vals_vec.size(); i++){
    cout << "Fit procedure " << i+1 << ":" << endl;
    cout << "  Aref = " << Aref[i] << ", Bref = " << Bref[i] << endl;
    for(int j=0; j<vals_vec[i].size(); j++){
      cout << "  " << vals_vec[i][j] << " +- " << errs_vec[i][j] << endl;
      ofile << vals_vec[i][j] << " +- " << errs_vec[i][j] << " ";
    }
    ofile << endl;
  }

  ofile.close();
  draw_fit_results(vals_vec,errs_vec,Aref,Bref);
  return 0;
}
