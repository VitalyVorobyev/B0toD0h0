#include "dalitzmodel.h"

#include "TH2F.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>

using namespace std;

PhaseSpace phsp;
const int GridSize = 2500;
//const double M_D = 1.86483;
//const double M_D_sq = M_D*M_D;
//const double M_Ks = 0.497614;
//const double M_Ks_sq = M_Ks*M_Ks;
//const double M_pi = 0.13957018;
//const double M_pi_sq = M_pi*M_pi;
//const double M_min = (M_pi+M_Ks)*(M_pi+M_Ks);
//const double M_max = (M_D-M_pi)*(M_D-M_pi);

const double Carr_CLEO[8] = { 0.710, 0.365, 0.106,-0.462,-0.884,-0.757,0.008, 0.481};
const double Sarr_CLEO[8] = {-0.013,-0.179,-1.063,-0.616,-0.162, 0.386,0.938,-0.147};
double C(const int bin){return Carr_CLEO[8-bin];}
double S(const int bin){return Sarr_CLEO[8-bin];}
//const double Karr_Belle[8] = { 0.067, 0.093, 0.030,0.089, 0.079, 0.102, 0.123, 0.159};
//const double Kbarr_Belle[8]= { 0.016, 0.027, 0.013,0.040, 0.015, 0.012, 0.026, 0.109};
const double Karr_CLEO[8]     = { 0.077, 0.098, 0.030, 0.080, 0.071, 0.099, 0.124, 0.165};
const double Karr_CLEO_err[8] = { 0.004, 0.004, 0.002, 0.004, 0.003, 0.004, 0.004, 0.004};
const double Kbarr_CLEO[8]    = { 0.020, 0.032, 0.013, 0.040, 0.018, 0.016, 0.029, 0.088};
const double Kbarr_CLEO_err[8]= { 0.002, 0.002, 0.001, 0.003, 0.002, 0.002, 0.002, 0.004};
double K(const int bin){return Karr_CLEO[8-bin];}
double Kb(const int bin){return Kbarr_CLEO[8-bin];}

const double Carr_CLEO_stat_err[8] = {0.034,0.071,0.105,0.100,0.056,0.099,0.080,0.080};
const double Carr_CLEO_sist_err[8] = {0.038,0.078,0.100,0.082,0.054,0.065,0.087,0.070};
const double Sarr_CLEO_stat_err[8] = {0.097,0.166,0.174,0.188,0.130,0.208,0.120,0.177};
const double Sarr_CLEO_sist_err[8] = {0.031,0.048,0.066,0.052,0.041,0.067,0.047,0.107};

double intC[8], intS[8], intK[8], intKb[8];

//EvtComplex amp_Belle2010(EvtVector4R p4_p,EvtVector4R moms1,EvtVector4R moms2, EvtVector4R moms3){
//  ** A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010 **
//  EvtResonance2 Res01(p4_p,moms1,moms2, 1.638, 133.2, 0.0508, .89166, 1);//K*(892)
//  EvtResonance2 Res02(p4_p,moms1,moms2, 2.210, 358.9, 0.294 , 1.412 , 0);//K0*(1430)
//  EvtResonance2 Res03(p4_p,moms1,moms2, 0.890, 314.8, 0.0985, 1.4256, 2);//K2*(1430)
//  EvtResonance2 Res04(p4_p,moms1,moms2, 0.880,  82.0, 0.322 , 1.717 , 1);//K*(1680)
//  EvtResonance2 Res05(p4_p,moms1,moms2, 0.650, 120.0, 0.232 , 1.414 , 1);//K*(1410)
//  EvtResonance2 Res06(p4_p,moms1,moms3, 0.149, 325.4, 0.0508, .89166, 1);//DCS K*(892)
//  EvtResonance2 Res07(p4_p,moms1,moms3, 0.360,  87.0, 0.294 , 1.412 , 0);//DCS K0*(1430)
//  EvtResonance2 Res08(p4_p,moms1,moms3, 0.230, 275.0, 0.0985, 1.4256, 2);//DCS K2*(1430)
//  EvtResonance2 Res09(p4_p,moms1,moms3, 2.100, 130.0, 0.322 , 1.717 , 1);//DCS K*(1680)
//  EvtResonance2 Res10(p4_p,moms1,moms3, 0.420, 253.0, 0.232 , 1.414 , 1);//DCS K*(1410)
//  EvtResonance2 Res11(p4_p,moms3,moms2, 1.000,   0.0, 0.1490, 0.7717, 1);//Rho
//  EvtResonance2 Res12(p4_p,moms3,moms2, .0343, 112.0, .00849, .78265, 1);//Omega
//  EvtResonance2 Res13(p4_p,moms3,moms2, 0.385, 207.3, 0.05  , 0.977,  0);//f0(980)
//  EvtResonance2 Res14(p4_p,moms3,moms2, 1.250,  69.0, 0.272 , 1.31 ,  0);//f0(1370)
//  EvtResonance2 Res15(p4_p,moms3,moms2, 1.440, 342.9, 0.1851, 1.2754, 2);//f2(1270)
//  EvtResonance2 Res16(p4_p,moms3,moms2, 0.490,  64.0, 0.400 , 1.465,  1);//Rho(1450)
//  EvtResonance2 Res17(p4_p,moms3,moms2, 1.560, 214.0, 0.453 , 0.522,  0);//sigma1
//  EvtResonance2 Res18(p4_p,moms3,moms2, 0.200, 212.0, 0.088 , 1.033,  0);//sigma2
//  return EvtComplex(-2.537, 0.923) + Res01.resAmpl() + Res02.resAmpl() +
//    Res03.resAmpl() + Res04.resAmpl() + Res05.resAmpl() + Res06.resAmpl() +
//    Res07.resAmpl() + Res08.resAmpl() + Res09.resAmpl() + Res10.resAmpl() +
//    Res11.resAmpl() + Res12.resAmpl() + Res13.resAmpl() + Res14.resAmpl() +
//    Res15.resAmpl() + Res16.resAmpl() + Res17.resAmpl() + Res18.resAmpl();
//}

//double E_Ks_star(const double& mp){return (mp-M_pi_sq+M_Ks_sq)*0.5/sqrt(mp);}
//double E_pi_star(const double& mp){return (M_D_sq-mp-M_pi_sq)*0.5/sqrt(mp);}

//double mm_max(const double& mp){
//  if(mp>M_max || mp<M_min) return -1;
//  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
//  const double sum = EKs+Epi;
//  const double PKs2 = EKs*EKs-M_Ks_sq;
//  const double Ppi2 = Epi*Epi-M_pi_sq;
//  const double dif = sqrt(PKs2)-sqrt(Ppi2);
//  return (sum-dif)*(sum+dif);
//}
//double mm_min(const double& mp){
//  if(mp>M_max || mp<M_min) return -1;
//  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
//  const double sum = EKs+Epi;
//  const double PKs2 = EKs*EKs-M_Ks_sq;
//  const double Ppi2 = Epi*Epi-M_pi_sq;
//  const double sumP = sqrt(PKs2)+sqrt(Ppi2);
//  return (sum-sumP)*(sum+sumP);
//}

//bool IsInPlot(const double& mp,const double& mm){
//  if(mm >= mm_min(mp) && mm <= mm_max(mp)) return true;
//  else return false;
//}

void changeAB(double& A, double& B){
    A = A + B;
    B = A - B;
    A = A - B;
    return;
}

//int GetModelBin(const double& mp,const double& mm){
//  EvtVector4R pd,pks,ppip,ppim;
//  Get4Vs(mp,mm,pd,pks,ppip,ppim);
//  EvtComplex A,Abar;
//  A     = amp_Belle2010(pd,pks,ppip,ppim);
//  Abar  = amp_Belle2010(pd,pks,ppim,ppip);
//  delta = arg(A) - arg(Abar);
//  if(delta<   -PI/8) delta += 2.*PI;
//  if(delta>15.*PI/8) delta -= 2.*PI;
//  for(int i=1; i<=8; i++){
//    if(PI*(i-1.5)/4 < delta && delta < PI*(i-0.5)/4){
//      if(mp>mm) return  i;
//      else      return -i;
//    }
//  }
//  return 0;
//}

//TH2F* hist;

//int GetBin(const double& mp,const double& mm){
//  int gbin = hist->FindBin(mp,mm);
//  return (int)floor(hist->GetBinContent(gbin)+0.01);
//}

void calc_integrals(const bool model_binning){
//  const int BINS[9] = {0,8,1,2,3,4,5,6,7};
//  const int BINS[9] = {0,7,6,5,4,3,2,1,8};
//  const int BINS[9] = {0,1,2,3,4,5,6,7,8};
  double mmMax;
  double mmMin;
  const double M_min = phsp.m_min();
  const double M_max = phsp.m_max();
  double mp,mm;
  const double dm = (M_max - M_min)/GridSize;
  double delta, P, Pbar;
  PhaseSpace phsp;
  DalitzModel dmodel;
//  Int_t colors[9] = {0,1,2,3,4,5,6,7,8}; // #colors >= #levels - 1
//  Int_t colors[9] = {0,1,9,2,3,6,4,5,7}; // #colors >= #levels - 1
  Int_t colors[9]   = {0,7,1,0,2,3,46,4,5};// // #colors >= #levels - 1
  gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  TH2I* bin_hist = new TH2I("bin_hist","bin_hist",500,M_min,M_max,500,M_min,M_max);
  bin_hist->SetStats(0);
  bin_hist->SetContour(9);
  bin_hist->SetMarkerSize(1.7);
  bin_hist->GetXaxis()->SetTitle("m_{+}^{2}, GeV^{2}/c^{4}");
  bin_hist->GetXaxis()->SetLabelSize(0.05);
  bin_hist->GetXaxis()->SetTitleSize(0.05);
  bin_hist->GetXaxis()->SetTitleOffset(0.8);
  bin_hist->GetYaxis()->SetTitle("m_{-}^{2}, GeV^{2}/c^{4}");
  bin_hist->GetYaxis()->SetLabelSize(0.05);
  bin_hist->GetYaxis()->SetTitleSize(0.05);
  bin_hist->GetYaxis()->SetTitleOffset(0.8);
  double Majorant = 0;
  int bin;//,binhist;
  for(int i=0; i<8; i++){ intC[i] = 0; intS[i] = 0; intK[i] = 0; intKb[i] = 0;}
  for(mp = M_min; mp<=M_max; mp += dm){
    mmMax = phsp.mm_max(mp);
    mmMin = phsp.mm_min(mp);
    mmMin = mmMin > mp ? mmMin : mp;
    for(mm = mmMin; mm<=mmMax; mm += dm){
      if(!phsp.IsInPlot(mp,mm)) continue;
      bin = model_binning ? abs(dmodel.GetBin(mp,mm)) : phsp.GetBin(mp,mm);
      if(!bin || abs(bin)>8) continue;
      bin_hist->SetBinContent(bin_hist->FindBin(mp,mm),abs(bin));
      bin_hist->SetBinContent(bin_hist->FindBin(mm,mp),abs(bin));
      dmodel.PPbarDelta(mp,mm,P,Pbar,delta);
      if(std::isnan(delta)) continue;
      if(P+Pbar > Majorant) Majorant = P+Pbar;
      if(bin < 0){
        changeAB(P,Pbar);
        delta = -delta;
      }
      bin = abs(bin) - 1;
      intC[bin]  += sqrt(P*Pbar)*TMath::Cos(delta);
      intS[bin]  += sqrt(P*Pbar)*TMath::Sin(delta);
      intK[bin]  += P;
      intKb[bin] += Pbar;
    }
  }
  double norm = 0;
  for(int i=0; i<8; i++){
    intC[i] /= sqrt(intK[i]*intKb[i]);
    intS[i] /= sqrt(intK[i]*intKb[i]);
    norm += intK[i] + intKb[i];
  }

  cout << "Norm     = " << norm << endl;
  cout << "Majorant = " << Majorant << endl;
  for(int i=0; i<8; i++){
    intK[i]  /= norm;
    intKb[i] /= norm;
    cout << i+1 << ": C = " << intC[i] << ", S = " << intS[i] << ", K = " << intK[i] << ", Kbar = " << intKb[i] << ", Q = " << intC[i]*intC[i] +intS[i]*intS[i] << endl;
  }

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  bin_hist->SetMinimum(-0.5);
  bin_hist->SetMaximum(8.5);
  bin_hist->Draw("ColZ");
  if(model_binning){
    c1->Print("model_binning.png");
    c1->Print("model_binning.eps");
    c1->Print("model_binning.root");
    system("display model_binning.png &");
  } else{
    c1->Print("cleo_binning.png");
    c1->Print("cleo_binning.eps");
    c1->Print("cleo_binning.root");
    system("display cleo_binning.png &");
  }

  return;
}

void draw_CS(const bool model_binning){
  double Cerr[8],Serr[8];
  double Carr[8],Sarr[8];
  for(int i=0; i<8; i++){
    Cerr[i] = sqrt(Carr_CLEO_sist_err[i]*Carr_CLEO_sist_err[i]+Carr_CLEO_stat_err[i]*Carr_CLEO_stat_err[i]);
    Serr[i] = sqrt(Sarr_CLEO_sist_err[i]*Sarr_CLEO_sist_err[i]+Sarr_CLEO_stat_err[i]*Sarr_CLEO_stat_err[i]);
    Carr[i] = C(i+1); Sarr[i] = S(i+1);
  }
  TGraph* grCS = new TGraph(8,intC,intS);
  grCS->SetMarkerStyle(20);
  grCS->SetMarkerSize(1);
  grCS->SetLineStyle(0);
  grCS->SetMarkerColor(kBlue);
  TGraphErrors* grCS_CLEO = new TGraphErrors(8,Carr,Sarr,Cerr,Serr);
  grCS_CLEO->SetMarkerStyle(20);
  grCS_CLEO->SetMarkerSize(1);
  grCS_CLEO->SetLineStyle(0);
  grCS_CLEO->SetMarkerColor(kRed);

  double circ_x[100],circ_y[100];
  for(int i=0; i<100; i++){
    circ_x[i] = cos(0.01*i*EvtConst::twoPi);
    circ_y[i] = sin(0.01*i*EvtConst::twoPi);
  }
  TGraph* gr_circ = new TGraph(100,circ_x,circ_y);
  gr_circ->SetMarkerStyle(20);
  gr_circ->SetMarkerSize(0.3);
  gr_circ->SetLineStyle(1);
  gr_circ->SetMarkerColor(kBlack);
  TMultiGraph* mgCS = new TMultiGraph("mgCS","mgCS");
  mgCS->Add(grCS);
  mgCS->Add(grCS_CLEO);
  mgCS->Add(gr_circ);

  TCanvas* c1 = new TCanvas("c1","c1",600,700);
  c1->cd();
  mgCS->Draw("ap");
  c1->Update();
  if(model_binning){
    c1->Print("CS_model.eps");
    c1->Print("CS_model.root");
    system("evince CS_model.eps &");
  } else{
    c1->Print("CS.eps");
    c1->Print("CS.root");
    system("evince CS.eps &");
  }
  return;
}

void draw_K(const bool model_binning){
  double bin[8]     = {1,2,3,4,5,6,7,8};
  double bin_err[8] = {0,0,0,0,0,0,0,0};
  double Karr[8],Kbarr[8];
  for(int i=0; i<8; i++){ Karr[i] = K(i+1); Kbarr[i] = Kb(i+1);}
  TGraph* grK  = new TGraph(8,bin,intK);
  grK->SetMarkerStyle(20);
  grK->SetMarkerSize(1);
  grK->SetLineStyle(0);
  grK->SetMarkerColor(kBlue);
  TGraph* grKb = new TGraph(8,bin,intKb);
  grKb->SetMarkerStyle(21);
  grKb->SetMarkerSize(1);
  grKb->SetLineStyle(0);
  grKb->SetMarkerColor(kBlue);
  TGraphErrors* grK_CLEO  = new TGraphErrors(8,bin,Karr,bin_err,Karr_CLEO_err);
  grK_CLEO->SetMarkerStyle(20);
  grK_CLEO->SetMarkerSize(1);
  grK_CLEO->SetLineStyle(0);
  grK_CLEO->SetMarkerColor(kRed);
  TGraphErrors* grKb_CLEO = new TGraphErrors(8,bin,Kbarr,bin_err,Kbarr_CLEO_err);
  grKb_CLEO->SetMarkerStyle(21);
  grKb_CLEO->SetMarkerSize(1);
  grKb_CLEO->SetLineStyle(0);
  grKb_CLEO->SetMarkerColor(kRed);

  TMultiGraph* mgK = new TMultiGraph("mgK","mgK");
  mgK->Add(grK);
  mgK->Add(grKb);
  mgK->Add(grK_CLEO);
  mgK->Add(grKb_CLEO);

  TCanvas* c2 = new TCanvas("c2","c2",600,700);
  c2->cd();
  c2->SetGrid();
  mgK->Draw("ap");
  c2->Update();
  if(model_binning){
    c2->Print("K_model.eps");
    c2->Print("K_model.root");
    system("evince K_model.eps &");
  } else{
    c2->Print("K.eps");
    c2->Print("K.root");
    system("evince K.eps &");
  }
  return;
}

int main(int argc, char** argv){
//  TFile* binning = TFile::Open("../Tuples/dkpp_belle_ddd.root");
//  hist = (TH2F*)binning->Get("dkpp_bin_h");
  bool model_binning = false;
  if(argc == 2){
    if(string(argv[1]) == string("model")) model_binning = true;
  }

  calc_integrals(model_binning);
  draw_CS(model_binning);
  draw_K(model_binning);
  return 0;
}
