#include "cuts.h"
using namespace RooFit;

void bdt_in_modes(void){
  TFile *pifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0.root");
  TTree *pitree = (TTree*)pifile->Get("TEvent");
  TFile *etafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigEta_full.root");
  TTree *etatree = (TTree*)etafile->Get("TEvent");
  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigOmega_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const bool save_flag = true;
  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("badpi",5);
  b0f.defineType("fsr",10);
  argset.add(b0f);

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin = -0.3;
  const double deMax = 0.3;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
//   mbc.setRange("Signal",mbc_min,mbc_max);
//   mbc.setRange("mbcSignal",mbc_min,mbc_max);
//   mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
//   de.setRange("Signal",de_min,de_max);
//   de.setRange("mbcSignal",deMin,deMax);
//   de.setRange("deSignal",de_min,de_max);
  RooRealVar bdtgsl("bdtgsl","bdtgsl",-1.,1.); argset.add(bdtgsl);
  bdtgsl.setRange("Cut",bdtgsl_cut,1.);
  RooDataSet pids("pids","pids",pitree,argset,"mbc>0||mbc<=0");
  pids.Print();
  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  h0mode.defineType("ppp",20);
  argset.add(h0mode);

  RooDataSet etads("etads","etads",etatree,argset,"mbc>0||mbc<=0");
  etads.Print();
  RooDataSet omegads("omegads","omegads",omegatree,argset,"mbc>0||mbc<=0");
  omegads.Print();
  stringstream out;

  /////////////////
  // BDTG Signal //
  /////////////////
  // pi //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min << " && de<" << de_max;
  RooDataSet* ds_pi = pids.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dh_pi = ds_pi->binnedClone();
  RooHistPdf* hpdf_pi = new RooHistPdf("hpdf_pi","hpdf_pi",RooArgSet(bdtgsl),*dh_pi);
  const double eff_pi = 100*hpdf_pi->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();

  // eta -> gg //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min << " && de<" << de_max << " && h0mode == 10";
  RooDataSet* ds_etagg = etads.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dh_etagg = ds_etagg->binnedClone();
  RooHistPdf* hpdf_etagg = new RooHistPdf("hpdf_etagg","hpdf_etagg",RooArgSet(bdtgsl),*dh_etagg);
  const double eff_etagg = 100*hpdf_etagg->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();

  // eta -> ppp //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min_etappp << " && de<" << de_max_etappp << " && h0mode == 20";
  RooDataSet* ds_etappp = etads.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dh_etappp = ds_etappp->binnedClone();
  RooHistPdf* hpdf_etappp = new RooHistPdf("hpdf_etappp","hpdf_etappp",RooArgSet(bdtgsl),*dh_etappp);
  const double eff_etappp = 100*hpdf_etappp->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();

  // eta -> omega //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min_omega << " && de<" << de_max_omega;
  RooDataSet* ds_omega = omegads.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dh_omega = ds_omega->binnedClone();
  RooHistPdf* hpdf_omega = new RooHistPdf("hpdf_omega","hpdf_omega",RooArgSet(bdtgsl),*dh_omega);
  const double eff_omega = 100*hpdf_omega->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();

  ///////////
  // Plots //
  ///////////
  // signal sl //
  RooPlot* bdtgSigslFrame = bdtgsl.frame();
  hpdf_pi->plotOn(bdtgSigslFrame,LineWidth(2),LineColor(kRed));
  hpdf_etagg->plotOn(bdtgSigslFrame,LineWidth(2));
  hpdf_etappp->plotOn(bdtgSigslFrame,LineWidth(2),LineColor(kRed),LineStyle(kDashed));
  hpdf_omega->plotOn(bdtgSigslFrame,LineWidth(2),LineColor(kBlack));

  TCanvas* cbdtgs_sl = new TCanvas("bdtgs_sl","BDTg",700,500);
  cbdtgs_sl->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetLogy();

  bdtgSigslFrame->GetXaxis()->SetTitleSize(0.05);
  bdtgSigslFrame->GetXaxis()->SetTitleOffset(0.85);
  bdtgSigslFrame->GetXaxis()->SetLabelSize(0.04);
  bdtgSigslFrame->GetYaxis()->SetTitleOffset(1.6);
  bdtgSigslFrame->Draw();

  out << setprecision(4);
  TPaveText *pt = new TPaveText(0.4,0.84,0.6,0.91,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "pi0->gg: " << eff_pi;
  pt->AddText(out.str().c_str());
  pt->SetTextColor(kRed);
  pt->Draw();
  TPaveText *pt2 = new TPaveText(0.4,0.77,0.6,0.84,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "eta->gg: " << eff_etagg;
  pt2->AddText(out.str().c_str());
  pt2->SetTextColor(kBlue);
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(0.4,0.70,0.6,0.77,"brNDC");
  pt3->SetFillColor(0);
  pt3->SetTextAlign(12);
  out.str("");
  out << "eta->ppp: " << eff_etappp;
  pt3->AddText(out.str().c_str());
  pt3->SetTextColor(kBlack);
  pt3->Draw();
  TPaveText *pt4 = new TPaveText(0.4,0.63,0.6,0.70,"brNDC");
  pt4->SetFillColor(0);
  pt4->SetTextAlign(12);
  out.str("");
  out << "omega->ppp: " << eff_omega;
  pt4->AddText(out.str().c_str());
  pt4->SetTextColor(kBlack);
  pt4->Draw();

  TLine *line_s = new TLine(bdtgsl_cut,0,bdtgsl_cut,0.1);
  line_s->SetLineColor(kBlue);
  line_s->SetLineStyle(1);
  line_s->SetLineWidth((Width_t)2.);
  line_s->Draw();
  if(save_flag) cbdtgs_sl->Print("../Reports/pics/bdtg_modes.png");

  return;
}
