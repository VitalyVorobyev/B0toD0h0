#include "cuts.h"
using namespace RooFit;

void draw_dz(){
  TFile *pi0file = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigPi0_s3new_full.root");
  TTree *pi0tree = (TTree*)pi0file->Get("TEvent");

//  TFile *etafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigEta_full.root");
//  TTree *etatree = (TTree*)etafile->Get("TEvent");

  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigOmega_s1new_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  TFile *genfile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_gen_0_2new.root");
  TTree *gentree = (TTree*)genfile->Get("TEvent");

  RooArgSet argset;
//  RooRealVar mbc("mbc","M_{bc}",5.2,5.23,"GeV");
//  argsetpi0.add(mbc); argsetpeta.add(mbc); argsetomega.add(mbc);
//  RooRealVar de("de","#DeltaE",-0.15,0.3,"GeV");
//  argsetpi0.add(de); argsetpeta.add(de); argsetomega.add(de);
  RooRealVar dz("dz","#Deltaz",-1,1,"mm");
  argset.add(dz);

  RooDataSet ds_gen("ds_gen","ds_gen",gentree,argset);
  ds_gen.Print();
  RooDataHist* dh_gen = ds_gen.binnedClone();
  RooHistPdf* hpdf_gen = new RooHistPdf("hpdf_gen","hpdf_gen",RooArgSet(dz),*dh_gen);

  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  h0mode.defineType("ppp",20);
  argset.add(h0mode);

  RooDataSet ds_pi0("ds_pi0","ds_pi0",pi0tree,argset);
  ds_pi0.Print();
  RooDataSet* dsr_pi0 = ds_pi0.reduce(RooArgSet(dz));
  RooDataHist* dh_pi0 = dsr_pi0->binnedClone();
  RooHistPdf* hpdf_pi0 = new RooHistPdf("hpdf_pi0","hpdf_pi0",RooArgSet(dz),*dh_pi0);

//  RooDataSet ds_etagg("ds_etagg","ds_etagg",etatree,argset,"h0mode == 10");
//  ds_etagg.Print();
//  RooDataSet* dsr_etagg = ds_etagg.reduce(RooArgSet(dz));
//  RooDataHist* dh_etagg = dsr_etagg->binnedClone();
//  RooHistPdf* hpdf_etagg = new RooHistPdf("hpdf_etagg","hpdf_etagg",RooArgSet(dz),*dh_etagg);

//  RooDataSet ds_etappp("ds_etappp","ds_etappp",etatree,argset,"h0mode == 20");
//  ds_etappp.Print();
//  RooDataSet* dsr_etappp = ds_etappp.reduce(RooArgSet(dz));
//  RooDataHist* dh_etappp = dsr_etappp->binnedClone();
//  RooHistPdf* hpdf_etappp = new RooHistPdf("hpdf_etappp","hpdf_etappp",RooArgSet(dz),*dh_etappp);

 RooDataSet ds_omega("ds_omega","ds_omega",omegatree,argset);
  ds_omega.Print();
  RooDataSet* dsr_omega = ds_omega.reduce(RooArgSet(dz));
  RooDataHist* dh_omega = dsr_omega->binnedClone();
  RooHistPdf* hpdf_omega = new RooHistPdf("hpdf_omega","hpdf_omega",RooArgSet(dz),*dh_omega);

///////////
// Plots //
///////////
  RooPlot* dzFrame = dz.frame();
  hpdf_pi0.plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
//  hpdf_etagg.plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kRed));
//  hpdf_etappp.plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kBlack));
  hpdf_omega.plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kRed));
  hpdf_gen.plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1),LineColor(kBlack));
  TCanvas* dz_cm = new TCanvas("#Deltaz, mm","#Deltaz, mm",700,500);
  dz_cm->cd();
  ds_pi0.statOn(dzFrame,Label("#pi^{0}"),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  ds_omega.statOn(dzFrame,Layout(0.6,0.98,0.7),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  ds_gen.statOn(dzFrame,Layout(0.15,0.49,0.9),What("MR"),Label("background"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  dzFrame->GetXaxis()->SetTitleOffset(0.75);
  dzFrame->GetXaxis()->SetLabelSize(0.05);
  dzFrame->GetXaxis()->SetTitleSize(0.06);
  dzFrame->GetYaxis()->SetTitleOffset(1.6);
  dzFrame->Draw();
  dz_cm->Update();
  dz_cm->Print("dz.png");
  dz_cm->Print("dz.root");

}
