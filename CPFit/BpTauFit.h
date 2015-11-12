#ifndef BPTAUFIT_H
#define BPTAUFIT_H

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
#include <string>

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
#include "TLine.h"
#include "TPaveText.h"

#include "RbkgPdf.h"
#include "RkRdetRnpPdf.h"

int m_svd = 0;
int good_icpv_mlt;
int good_icpv_sgl;

bool generic = false;
bool signal  = false;
bool mlt_asc = false;
bool sgl_asc = false;

bool single_d0 = false;
const double m_btau = 1.638;
double m_fbkg_svd1 = 0;
double m_fbkg_svd2 = 0;

RbkgPdf* m_pdf_back_svd1;
RbkgPdf* m_pdf_back_svd2;
RkRdetRnpPdf* m_pdf_sig_svd1;
RkRdetRnpPdf* m_pdf_sig_svd2;
TChain* m_tree;
TTree* m_good_tree;
Int_t m_exp;
Double_t m_dt;

// * Event-by-event parameters
int ntrk_asc;
int ndf_rec;
int ndf_asc;
double sz_rec;
double sz_asc;
double chisq_rec;
double chisq_asc;
double m_costhBcms;
// *

void SetBkgPdfParams(const vector<double>& par){
  m_pdf_back_svd1->SetTau(par[0]);
  m_pdf_back_svd1->Set_mu(par[1]);
  m_pdf_back_svd1->Set_mu_delta(par[2]);

  m_pdf_back_svd1->Set_f_delta_mlt(par[3]);
  m_pdf_back_svd1->Set_f_tail_mlt(par[4]);
  m_pdf_back_svd1->Set_S_main_mlt(par[5]);
  m_pdf_back_svd1->Set_S_tail_mlt(par[6]);

  m_pdf_back_svd1->Set_f_delta_sgl(par[7]);
  m_pdf_back_svd1->Set_f_tail_sgl(par[8]);
  m_pdf_back_svd1->Set_S_main_sgl(par[9]);
  m_pdf_back_svd1->Set_S_tail_sgl(par[10]);

  m_pdf_back_svd2->SetTau(par[11]);
  m_pdf_back_svd2->Set_mu(par[12]);
  m_pdf_back_svd2->Set_mu_delta(par[13]);

  m_pdf_back_svd2->Set_f_delta_mlt(par[14]);
  m_pdf_back_svd2->Set_f_tail_mlt(par[15]);
  m_pdf_back_svd2->Set_S_main_mlt(par[16]);
  m_pdf_back_svd2->Set_S_tail_mlt(par[17]);

  m_pdf_back_svd2->Set_f_delta_sgl(par[18]);
  m_pdf_back_svd2->Set_f_tail_sgl(par[19]);
  m_pdf_back_svd2->Set_S_main_sgl(par[20]);
  m_pdf_back_svd2->Set_S_tail_sgl(par[21]);
  return;
}

using namespace std;
using namespace ROOT;
using namespace Minuit2;

bool draw_plots = true;
const double sol = 2.99792458;
const double cm2ps = 78.48566945838871754705;

const int NDots = 250;
const int NBins = NDots;
const double dtmax = 70;
const double dtmin = -dtmax;

double sum_sigma(const double& s1, const double& s2){
  return sqrt(s1*s1+s2*s2);
}

int init(){
  m_pdf_back_svd1->SetRange(dtmax);
  m_pdf_back_svd2->SetRange(dtmax);
  m_pdf_sig_svd1->SetRange(dtmax);
  m_pdf_sig_svd2->SetRange(dtmax);
  if(signal){
    m_fbkg_svd1 = 0;
    m_fbkg_svd2 = 0;
  } else if(generic){
    m_fbkg_svd1 = 0.0596;
    m_fbkg_svd2 = 0.0596;
  } else{// data
    if(single_d0){
      m_fbkg_svd1 = 0.0826;// +- 0.098
      m_fbkg_svd2 = 0.0911;// +- 0.098
    } else{
      m_fbkg_svd1 = 0.0820;// +- 0.093
      m_fbkg_svd2 = 0.0912;// +- 0.093
    }
  }
  return 0;
}

void GetEvent(const int i){
  m_good_tree->GetEvent(i);
  m_dt *= cm2ps;
  if(signal) m_dt *= 0.1;
  sz_rec /= 10.;// mm2cm
  sz_asc /= 10.;// mm2cm
  return;
}

double Pdf(const double dt){
  const double fbkg = m_exp > 30 ? m_fbkg_svd2 : m_fbkg_svd1;
  m_exp > 30 ? m_pdf_sig_svd2->SetAkCk(m_costhBcms,0.5*10.58) : m_pdf_sig_svd1->SetAkCk(m_costhBcms,0.5*10.58);
  double pdf = 0;
  if(m_exp > 30){// svd 2
    if(single_d0) pdf = (1-fbkg)*m_pdf_sig_svd2->Pdf(dt,1,sz_rec,chisq_rec,0,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,true);
    else          pdf = (1-fbkg)*m_pdf_sig_svd2->Pdf(dt,2,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,true);
    if(fbkg>0.001)pdf += fbkg*m_pdf_back_svd2->Pdf(dt,sum_sigma(sz_asc,sz_rec),ndf_asc);
  } else{// svd 1
    if(single_d0) pdf = (1-fbkg)*m_pdf_sig_svd1->Pdf(dt,1,sz_rec,chisq_rec,0,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,true);
    else          pdf = (1-fbkg)*m_pdf_sig_svd1->Pdf(dt,2,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,true);
    if(fbkg>0.001)pdf += fbkg*m_pdf_back_svd1->Pdf(dt,sum_sigma(sz_asc,sz_rec),ndf_asc);
  }
  return pdf;
}

TTree* GetGoodTTree(TTree* tree, const bool sideband){
  const double demin = -0.03;
  const double demax = 0.04;
  const double bdtg_cut = -0.44;
  const double mbcmin = sideband ? 5.200 : 5.272;
  const double mbcmax = sideband ? 5.260 : 5.289;
//  const double chisq_sig_max = single_d0 ? 1000 : 10;
//  const double chisq_asc_max = 10;
//  const double sz_sig_max = 0.2;// cm
//  const double sz_asc_max = 0.2;// cm

  stringstream out;
  out.str("");
  if(signal) out << "(bpf == 1 || bpf == 5 || bpf == 10) && ";
  if(m_svd == 2) out << "exp>30 && ";
  if(m_svd == 1) out << "exp<30 && ";
  out << "de<" << demax << " && " << "de>" << demin << " && ";
  out << "mbc<" << mbcmax << " && " << "mbc>" << mbcmin << " && ";
//  out << "chisq_z_sig>0 && chisq_z_sig/ndf_z_sig<" << chisq_sig_max << " && ";
//  out << "chisq_z_asc>0 && chisq_z_asc/ndf_z_asc<" << chisq_asc_max << " && ";
  if(single_d0) out << "good_icpv_sgl == 1 && ";
  else          out << "good_icpv_mlt == 1 && ";
  out << "bdtg>" << bdtg_cut << " && ";
  if(mlt_asc) out << "ndf_z_asc!=0 && ";
  if(sgl_asc) out << "ndf_z_asc==0 && ";
  if(single_d0){
    if(!signal){
      out << "abs(dz_d0)<" << dtmax/cm2ps;
    } else{
      out << "abs(dz_d0)*0.1<" << dtmax/cm2ps;
    }
  } else{
    if(!signal){
      out << "abs(dz)<" << dtmax/cm2ps;
    } else{
      out << "abs(dz)*0.1<" << dtmax/cm2ps;
    }
  }
  return tree->CopyTree(out.str().c_str());
}

void SetPDFParams(const vector<double>& par){
  m_pdf_sig_svd1->SetTauDm(par[0],0);
  m_pdf_sig_svd2->SetTauDm(par[0],0);

  m_pdf_sig_svd1->Set_f_ol_sgl(par[1]);
  m_pdf_sig_svd1->Set_f_ol_mlt(par[2]);

  m_pdf_sig_svd2->Set_f_ol_sgl(par[3]);
  m_pdf_sig_svd2->Set_f_ol_mlt(par[4]);
  return;
}

class pdfFcnBkg : public FCNBase{
public:
  pdfFcnBkg(void){
    theErrorDef = 1;
    NTot = m_good_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnBkg() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    SetBkgPdfParams(par);
    double loglh = 0;
    double pdf;
    for(int i=0; i<NTot; i++){
      GetEvent(i);
      const double s = sum_sigma(sz_asc,sz_rec);
      pdf = m_exp>30 ? m_pdf_back_svd2->Pdf(m_dt,s,ndf_asc) : m_pdf_back_svd1->Pdf(m_dt,s,ndf_asc);
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
    double pdf;
    SetPDFParams(par);
    for(int i=0; i<NTot; i++){
      GetEvent(i);
      pdf = Pdf(m_dt);
      if(!std::isnan(pdf) && pdf>0) loglh += -2*TMath::Log(pdf);
      else{ cout << "pdfFcn: pdf = " << pdf << endl;}
    }
    cout << "loglh: " << loglh << ", tau: " << par[0];// << ", f_olr_sgl = " << par[1] << ", f_olr_mlt = " << par[2];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

void Draw_NoTag(const int NFreePar){
  const int NTot = m_good_tree->GetEntries();
  const double ddt = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];

  for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }

  double Norm = 0;
  double pdf_array[NDots];

  stringstream out;
  out.str("");
  out << "Lifetime fit D^{0}#pi^{+}";
  if(single_d0)   out << " (only D^{0})";
  if(m_svd == 2)  out << ", SVD2";
  if(m_svd == 1)  out << ", SVD1";
  if(generic)     out << ", GenMC";
  else if(signal) out << ", SigMC";
  else            out << ", Data";
  if(sgl_asc) out << " (sgl)";
  if(mlt_asc) out << " (mlt)";

  TH1I* DH = new TH1I("DH",out.str().c_str(),NBins,dtmin,dtmax);
  DH->SetMarkerStyle(21);
  DH->SetMarkerColor(kBlue);
  DH->SetMarkerSize(1.1);

  double pdfval;
  for(int i=0; i<NDots; i++) pdf_array[i] = 0;
  for(int j=0; j<NTot; j++){
    GetEvent(j);
    DH->Fill(m_dt);
    for(int i=0; i<NDots; i++){
      pdfval = Pdf(dt_arr[i]);
      if(!std::isnan(pdfval)) pdf_array[i] += pdfval;
      else cout << "pdf is nan: " << dt_arr[i] << " " << sz_rec << " " << chisq_rec << " " << ndf_rec << " " << ntrk_asc << " " << sz_asc << " " << chisq_asc << " " << ndf_asc << endl;
    }
  }

  double pull_array[NBins];
  const int nBD = NDots/NBins;
  double chisq = 0;
  for(int i=0; i<NDots; i++){
    pdf_array[i] *= ddtb;
    Norm += pdf_array[i];
    if(!(i%nBD)){
      const int bin = i/nBD;
      const int bin_content = DH->GetBinContent(bin+1);
      if(bin_content){
        pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
        chisq += pull_array[bin]*pull_array[bin];
      }
    }
  }
  cout << "Norm = " << Norm/nBD << ", Nevents = " << NTot << endl;
  chisq /= NBins;
  cout << "chi2/n.d.f. = " << chisq << endl;

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetLineWidth(2);
  GR->SetMarkerColor(kBlue);

  TGraphErrors* GRP = new TGraphErrors(NBins,dt_barr,pull_array,dt_barr_err,dt_pull_err);
  GRP->SetMarkerStyle(21);
  GRP->SetMarkerSize(1.2);
  GRP->SetMarkerColor(kBlue);
  GRP->SetLineColor(0);

  TCanvas* c1 = new TCanvas("c1","c1",600,800);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.19);
  pad1->Draw(); pad2->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  DH->GetXaxis()->SetTitle("#Deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.06);
  DH->GetXaxis()->SetTitleOffset(0.75);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->Draw("e");

  TPaveText *pt = new TPaveText(0.75,0.65,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt->AddText(out.str().c_str());
  pt->Draw();

  GR->Draw("same");

  pad2->cd();
  GRP->GetYaxis()->SetRangeUser(-5,5);
  GRP->GetXaxis()->SetRangeUser(dtmin,dtmax);
  GRP->Draw();
  TLine* zeroline = new TLine(dtmin,0,dtmax,0);
  zeroline->SetLineWidth(2);
  zeroline->Draw("AP");

  TLine* plusline = new TLine(dtmin,3,dtmax,3);
  plusline->SetLineWidth(1);
  plusline->SetLineStyle(kDashed);
  plusline->Draw();

  TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
  minusline->SetLineWidth(1);
  minusline->SetLineStyle(kDashed);
  minusline->Draw();

  out.str("");
  out << "lifetime_d0pip";
  if(single_d0)    out << "_d0";
  if(mlt_asc)      out << "_mlt_";
  if(sgl_asc)      out << "_sgl_";
  if(m_svd == 2)   out << "svd2";
  if(m_svd == 1)   out << "svd1";
  if(signal)       out << "_sig";
  else if(generic) out << "_gen";
  else             out << "data";
  if(!NFreePar)    out << "_def";
  string rootname = out.str() + string(".root");

  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());
  return;
}

#endif // BPTAUFIT_H
