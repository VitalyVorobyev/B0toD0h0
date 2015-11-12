//#include "cuts.h"
using namespace RooFit;

void draw_dz1(){
  TFile *pi0file = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigPi0_s3new_full.root");
  TTree *pi0tree = (TTree*)pi0file->Get("TEvent");

  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigOmega_s1new_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const double dzmin = -0.5;
  const double dzmax =  0.5;
  RooArgSet argset;
  RooRealVar dz_mc_sig1("dz_mc_sig1","#deltaz_{sig}",dzmin,dzmax,"cm");
  RooRealVar dz_mc_asc1("dz_mc_asc1","#deltaz_{asc}",dzmin,dzmax,"cm");
  argset.add(dz_mc_sig1); argset.add(dz_mc_asc1);

  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  argset.add(b0f);

  RooDataSet ds_pi0("ds_pi0","ds_pi0",pi0tree,argset);
  ds_pi0.Print();
  RooDataSet* dsr_pi0_sig = ds_pi0.reduce(RooArgSet(dz_mc_sig1));
  RooDataHist* dh_pi0_sig = dsr_pi0_sig->binnedClone();
  RooHistPdf* hpdf_pi0_sig = new RooHistPdf("hpdf_pi0_sig","hpdf_pi0_sig",RooArgSet(dz_mc_sig1),*dh_pi0_sig);

  RooDataSet* dsr_pi0_asc = ds_pi0.reduce(RooArgSet(dz_mc_asc1));
  RooDataHist* dh_pi0_asc = dsr_pi0_asc->binnedClone();
  RooHistPdf* hpdf_pi0_asc = new RooHistPdf("hpdf_pi0_asc","hpdf_pi0_asc",RooArgSet(dz_mc_asc1),*dh_pi0_asc);

  RooDataSet ds_omega("ds_omega","ds_omega",omegatree,argset);
  ds_omega.Print();
  RooDataSet* dsr_omega_sig = ds_omega.reduce(RooArgSet(dz_mc_sig1));
  RooDataHist* dh_omega_sig = dsr_omega_sig->binnedClone();
  RooHistPdf* hpdf_omega_sig = new RooHistPdf("hpdf_omega_sig","hpdf_omega_sig",RooArgSet(dz_mc_sig1),*dh_omega_sig);

  RooDataSet* dsr_omega_asc = ds_omega.reduce(RooArgSet(dz_mc_asc1));
  RooDataHist* dh_omega_asc = dsr_omega_asc->binnedClone();
  RooHistPdf* hpdf_omega_asc = new RooHistPdf("hpdf_omega_asc","hpdf_omega_asc",RooArgSet(dz_mc_asc1),*dh_omega_asc);

///////////
// Plots //
///////////
  RooPlot* dzFrameSig = dz_mc_sig1.frame();
  hpdf_pi0_sig->plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
  hpdf_omega_sig->plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kRed));
  TCanvas* dz_cm_sig = new TCanvas("#deltaz_{sig}, cm","#deltaz_{sig}, cm",700,500);
  dz_cm_sig->cd();
  dsr_pi0_sig->statOn(dzFrameSig,Label("#pi^{0}"),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  dsr_omega_sig->statOn(dzFrameSig,Layout(0.6,0.98,0.7),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  dzFrameSig->GetXaxis()->SetTitleOffset(0.75);
  dzFrameSig->GetXaxis()->SetLabelSize(0.05);
  dzFrameSig->GetXaxis()->SetTitleSize(0.06);
  dzFrameSig->GetYaxis()->SetTitleOffset(1.6);
  dzFrameSig->Draw();
  dz_cm_sig->Update();
  dz_cm_sig->Print("dz_sig.png");
  dz_cm_sig->Print("dz_sig.root");


  RooPlot* dzFrameAsc = dz_mc_asc1.frame();
  hpdf_pi0_asc->plotOn(dzFrameAsc,DataError(RooAbsData::SumW2),MarkerSize(1));
  hpdf_omega_asc->plotOn(dzFrameAsc,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kRed));
  TCanvas* dz_cm_asc = new TCanvas("#deltaz_{asc}, mm","#deltaz_{asc}, mm",700,500);
  dz_cm_asc->cd();
  dsr_pi0_asc->statOn(dzFrameAsc,Label("#pi^{0}"),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  dsr_omega_asc->statOn(dzFrameAsc,Layout(0.6,0.98,0.7),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  TPad *pad2 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad2->Draw();
  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  dzFrameAsc->GetXaxis()->SetTitleOffset(0.75);
  dzFrameAsc->GetXaxis()->SetLabelSize(0.05);
  dzFrameAsc->GetXaxis()->SetTitleSize(0.06);
  dzFrameAsc->GetYaxis()->SetTitleOffset(1.6);
  dzFrameAsc->Draw();
  dz_cm_asc->Update();
  dz_cm_asc->Print("dz_asc.png");
  dz_cm_asc->Print("dz_asc.root");

}
