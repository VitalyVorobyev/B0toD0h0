#include "cuts.h"
using namespace RooFit;

void bdt_in_regions(void){
  TFile *contfile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_cont.root");
  TTree *conttree = (TTree*)contfile->Get("TEvent");
  TFile *sigfile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0.root");
  TTree *sigtree = (TTree*)sigfile->Get("TEvent");

  const bool save_flag = false;
  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
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
  RooRealVar bdtgfr("bdtgfr","bdtgfr",-1.,1.); argset.add(bdtgfr);
  bdtgfr.setRange("Cut",bdtgfr_cut,1.);
  RooDataSet sigds("sigds","sigds",sigtree,argset,"mbc>0||mbc<=0 && b0f != -1");
  sigds.Print();
  RooDataSet contds("contds","contds",conttree,argset,"mbc>0||mbc<=0 && b0f == -1");
  contds.Print();
  stringstream out;

  /////////////////
  // BDTG Signal //
  /////////////////
  // soft cut //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min << " && de<" << de_max;
  RooDataSet* dssl_s_sig = sigds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_s_sig = dssl_s_sig->binnedClone();
  RooHistPdf* hslpdf_s_sig = new RooHistPdf("hslpdf_s_sig","hslpdf_s_sig",RooArgSet(bdtgsl),*dhsl_s_sig);
  const double sl_eff_s = 100*hslpdf_s_sig->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();
//   delete dhsl_s_sig;
//   delete dssl_s_sig;
  RooDataSet* dsfr_s_sig = sigds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_s_sig = dsfr_s_sig->binnedClone();
  RooHistPdf* hfrpdf_s_sig = new RooHistPdf("hfrpdf_s_sig","hfrpdf_s_sig",RooArgSet(bdtgfr),*dhfr_s_sig);
  const double fr_eff_s = 100*hfrpdf_s_sig->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal();
//   delete dhfr_s_sig;
//   delete dsfr_s_sig;

  // soft cut //
  out.str("");
  out << "mbc<5.289 && mbc>5.250 && de<0.1 && de>-0.1";
  RooDataSet* dssl_side_sig = sigds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_side_sig = dssl_side_sig->binnedClone();
  RooHistPdf* hslpdf_side_sig = new RooHistPdf("hslpdf_side_sig","hslpdf_side_sig",RooArgSet(bdtgsl),*dhsl_side_sig);
  const double sl_eff_sl = 100*hslpdf_side_sig->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();
//   delete dhsl_side_sig;
//   delete dssl_side_sig;
  RooDataSet* dsfr_side_sig = sigds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_side_sig = dsfr_side_sig->binnedClone();
  RooHistPdf* hfrpdf_side_sig = new RooHistPdf("hfrpdf_side_sig","hfrpdf_side_sig",RooArgSet(bdtgfr),*dhfr_side_sig);
  const double fr_eff_sl = 100*hfrpdf_side_sig->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal();
//   delete dhfr_side_sig;
//   delete dsfr_side_sig;

  // full range //
  out.str("");
  RooDataSet* dssl_full_sig = sigds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_full_sig = dssl_full_sig->binnedClone();
  RooHistPdf* hslpdf_full_sig = new RooHistPdf("hslpdf_full_sig","hslpdf_full_sig",RooArgSet(bdtgsl),*dhsl_full_sig);
  const double sl_eff_fr = 100*hslpdf_full_sig->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal();
//   delete dhsl_full_sig;
//   delete dssl_full_sig;
  RooDataSet* dsfr_full_sig = sigds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_full_sig = dsfr_full_sig->binnedClone();
  RooHistPdf* hfrpdf_full_sig = new RooHistPdf("hfrpdf_full_sig","hfrpdf_full_sig",RooArgSet(bdtgfr),*dhfr_full_sig);
  const double fr_eff_fr = 100*hfrpdf_full_sig->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal();
//   delete dhfr_full_sig;
//   delete dsfr_full_sig;

  ////////////////////
  // BDTG Continuum //
  ////////////////////
  // soft cut //
  out.str("");
  out << "mbc>" << mbc_min << " && mbc<" << mbc_max << " && de>" << de_min << " && de<" << de_max;
  RooDataSet* dssl_s_cont = contds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_s_cont = dssl_s_cont->binnedClone();
  RooHistPdf* hslpdf_s_cont = new RooHistPdf("hslpdf_s_cont","hslpdf_s_cont",RooArgSet(bdtgsl),*dhsl_s_cont);
  const double sl_reg_s = 100*(1-hslpdf_s_cont->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal());
//   delete dhsl_s_cont;
//   delete dssl_s_cont;
  RooDataSet* dsfr_s_cont = contds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_s_cont = dsfr_s_cont->binnedClone();
  RooHistPdf* hfrpdf_s_cont = new RooHistPdf("hfrpdf_s_cont","hfrpdf_s_cont",RooArgSet(bdtgfr),*dhfr_s_cont);
  const double fr_reg_s = 100*(1-hfrpdf_s_cont->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal());
//   delete dhfr_s_cont;
//   delete dsfr_s_cont;

  // soft cut //
  out.str("");
  out << "mbc<5.289 && mbc>5.250 && de<0.1 && de>-0.1";
  RooDataSet* dssl_side_cont = contds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_side_cont = dssl_side_cont->binnedClone();
  RooHistPdf* hslpdf_side_cont = new RooHistPdf("hslpdf_side_cont","hslpdf_side_cont",RooArgSet(bdtgsl),*dhsl_side_cont);
  const double sl_reg_sl = 100*(1-hslpdf_side_cont->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal());
//   delete dhsl_side_cont;
//   delete dssl_side_cont;
  RooDataSet* dsfr_side_cont = contds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_side_cont = dsfr_side_cont->binnedClone();
  RooHistPdf* hfrpdf_side_cont = new RooHistPdf("hfrpdf_side_cont","hfrpdf_side_cont",RooArgSet(bdtgfr),*dhfr_side_cont);
  const double fr_reg_sl = 100*(1-hfrpdf_side_cont->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal());
//   delete dhfr_side_cont;
//   delete dsfr_side_cont;

  // full range //
  out.str("");
  RooDataSet* dssl_full_cont = contds.reduce(RooArgSet(bdtgsl),out.str().c_str());
  RooDataHist* dhsl_full_cont = dssl_full_cont->binnedClone();
  RooHistPdf* hslpdf_full_cont = new RooHistPdf("hslpdf_full_cont","hslpdf_full_cont",RooArgSet(bdtgsl),*dhsl_full_cont);
  const double sl_reg_fr = 100*(1-hslpdf_full_cont->createIntegral(RooArgSet(bdtgsl),NormSet(RooArgSet(bdtgsl)),Range("Cut"))->getVal());
//   delete dhsl_full_cont;
//   delete dssl_full_cont;
  RooDataSet* dsfr_full_cont = contds.reduce(RooArgSet(bdtgfr),out.str().c_str());
  RooDataHist* dhfr_full_cont = dsfr_full_cont->binnedClone();
  RooHistPdf* hfrpdf_full_cont = new RooHistPdf("hfrpdf_full_cont","hfrpdf_full_cont",RooArgSet(bdtgfr),*dhfr_full_cont);
  const double fr_reg_fr = 100*(1-hfrpdf_full_cont->createIntegral(RooArgSet(bdtgfr),NormSet(RooArgSet(bdtgfr)),Range("Cut"))->getVal());
//   delete dhfr_full_cont;
//   delete dsfr_full_cont;

  ///////////
  // Plots //
  ///////////
  // signal sl //
  RooPlot* bdtgSigslFrame = bdtgsl.frame();
  hslpdf_s_sig->plotOn(bdtgSigslFrame,LineWidth(2),LineColor(kRed));
  hslpdf_side_sig->plotOn(bdtgSigslFrame,LineWidth(2));
  hslpdf_full_sig->plotOn(bdtgSigslFrame,LineWidth(2),LineColor(kBlack));

  TCanvas* cbdtgs_sl = new TCanvas("bdtgs_sl","BDTg signal soft cut",700,500);
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
  out << "SR: " << sl_eff_s;
  pt->AddText(out.str().c_str());
  pt->SetTextColor(kRed);
  pt->Draw();
  TPaveText *pt2 = new TPaveText(0.4,0.77,0.6,0.84,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "SC: " << sl_eff_sl;
  pt2->AddText(out.str().c_str());
  pt2->SetTextColor(kBlue);
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(0.4,0.70,0.6,0.77,"brNDC");
  pt3->SetFillColor(0);
  pt3->SetTextAlign(12);
  out.str("");
  out << "FR: " << sl_eff_fr;
  pt3->AddText(out.str().c_str());
  pt3->SetTextColor(kBlack);
  pt3->Draw();

  TLine *line_s = new TLine(bdtgsl_cut,0,bdtgsl_cut,0.1);
  line_s->SetLineColor(kBlue);
  line_s->SetLineStyle(1);
  line_s->SetLineWidth((Width_t)2.);
  line_s->Draw();
  if(save_flag) cbdtgs_sl->Print("../Reports/pics/bdtg_diff_sl_sig.png");

  // signal fr //
  RooPlot* bdtgSigfrFrame = bdtgfr.frame();
  hfrpdf_s_sig->plotOn(bdtgSigfrFrame,LineWidth(2),LineColor(kRed));
  hfrpdf_side_sig->plotOn(bdtgSigfrFrame,LineWidth(2));
  hfrpdf_full_sig->plotOn(bdtgSigfrFrame,LineWidth(2),LineColor(kBlack));

  TCanvas* cbdtgs_fr = new TCanvas("bdtgs_fr","BDTg signal full range",700,500);
  cbdtgs_fr->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetLogy();

  bdtgSigfrFrame->GetXaxis()->SetTitleSize(0.05);
  bdtgSigfrFrame->GetXaxis()->SetTitleOffset(0.85);
  bdtgSigfrFrame->GetXaxis()->SetLabelSize(0.04);
  bdtgSigfrFrame->GetYaxis()->SetTitleOffset(1.6);
  bdtgSigfrFrame->Draw();
  
  TPaveText *pt = new TPaveText(0.4,0.84,0.6,0.91,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "SR: " << fr_eff_s;
  pt->AddText(out.str().c_str());
  pt->SetTextColor(kRed);
  pt->Draw();
  TPaveText *pt2 = new TPaveText(0.4,0.77,0.6,0.84,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "SC: " << fr_eff_sl;
  pt2->AddText(out.str().c_str());
  pt2->SetTextColor(kBlue);
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(0.4,0.70,0.6,0.77,"brNDC");
  pt3->SetFillColor(0);
  pt3->SetTextAlign(12);
  out.str("");
  out << "FR: " << fr_eff_fr;
  pt3->AddText(out.str().c_str());
  pt3->SetTextColor(kBlack);
  pt3->Draw();

  TLine *line_s = new TLine(bdtgfr_cut,0,bdtgfr_cut,0.04);
  line_s->SetLineColor(kBlue);
  line_s->SetLineStyle(1);
  line_s->SetLineWidth((Width_t)2.);
  line_s->Draw();
  if(save_flag) cbdtgs_fr->Print("../Reports/pics/bdtg_diff_fr_sig.png");

  // continuum sl //
  RooPlot* bdtgContslFrame = bdtgsl.frame();
  hslpdf_s_cont->plotOn(bdtgContslFrame,LineWidth(2),LineColor(kRed));
  hslpdf_side_cont->plotOn(bdtgContslFrame,LineWidth(2));
  hslpdf_full_cont->plotOn(bdtgContslFrame,LineWidth(2),LineColor(kBlack));

  TCanvas* cbdtgc_sl = new TCanvas("cbdtgc_sl","BDTg continuum soft cut",700,500);
  cbdtgc_sl->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetLogy();

  bdtgContslFrame->GetXaxis()->SetTitleSize(0.05);
  bdtgContslFrame->GetXaxis()->SetTitleOffset(0.85);
  bdtgContslFrame->GetXaxis()->SetLabelSize(0.04);
  bdtgContslFrame->GetYaxis()->SetTitleOffset(1.6);
  bdtgContslFrame->Draw();
  
  TPaveText *pt = new TPaveText(0.4,0.84,0.6,0.91,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "SR: " << sl_reg_s;
  pt->AddText(out.str().c_str());
  pt->SetTextColor(kRed);
  pt->Draw();
  TPaveText *pt2 = new TPaveText(0.4,0.77,0.6,0.84,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "SC: " << sl_reg_sl;
  pt2->AddText(out.str().c_str());
  pt2->SetTextColor(kBlue);
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(0.4,0.70,0.6,0.77,"brNDC");
  pt3->SetFillColor(0);
  pt3->SetTextAlign(12);
  out.str("");
  out << "FR: " << sl_reg_fr;
  pt3->AddText(out.str().c_str());
  pt3->SetTextColor(kBlack);
  pt3->Draw();

  TLine *line_s = new TLine(bdtgsl_cut,0,bdtgsl_cut,0.04);
  line_s->SetLineColor(kBlue);
  line_s->SetLineStyle(1);
  line_s->SetLineWidth((Width_t)2.);
  line_s->Draw();
  if(save_flag) cbdtgc_sl->Print("../Reports/pics/bdtg_diff_sl_cont.png");

  // continuum fr //
  RooPlot* bdtgContfrFrame = bdtgfr.frame();
  hfrpdf_s_cont->plotOn(bdtgContfrFrame,LineWidth(2),LineColor(kRed));
  hfrpdf_side_cont->plotOn(bdtgContfrFrame,LineWidth(2));
  hfrpdf_full_cont->plotOn(bdtgContfrFrame,LineWidth(2),LineColor(kBlack));

  TCanvas* cbdtgc_fr = new TCanvas("cbdtgc_fr","BDTg continuum full region",700,500);
  cbdtgc_fr->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetLogy();

  bdtgContfrFrame->GetXaxis()->SetTitleSize(0.05);
  bdtgContfrFrame->GetXaxis()->SetTitleOffset(0.85);
  bdtgContfrFrame->GetXaxis()->SetLabelSize(0.04);
  bdtgContfrFrame->GetYaxis()->SetTitleOffset(1.6);
  bdtgContfrFrame->Draw();

  TPaveText *pt = new TPaveText(0.4,0.84,0.6,0.91,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "SR: " << fr_reg_s;
  pt->AddText(out.str().c_str());
  pt->SetTextColor(kRed);
  pt->Draw();
  TPaveText *pt2 = new TPaveText(0.4,0.77,0.6,0.84,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "SC: " << fr_reg_sl;
  pt2->AddText(out.str().c_str());
  pt2->SetTextColor(kBlue);
  pt2->Draw();
  TPaveText *pt3 = new TPaveText(0.4,0.70,0.6,0.77,"brNDC");
  pt3->SetFillColor(0);
  pt3->SetTextAlign(12);
  out.str("");
  out << "FR: " << fr_reg_fr;
  pt3->AddText(out.str().c_str());
  pt3->SetTextColor(kBlack);
  pt3->Draw();

  TLine *line_s = new TLine(bdtgfr_cut,0,bdtgfr_cut,0.1);
  line_s->SetLineColor(kBlue);
  line_s->SetLineStyle(1);
  line_s->SetLineWidth((Width_t)2.);
  line_s->Draw();
  if(save_flag) cbdtgc_fr->Print("../Reports/pics/bdtg_diff_fr_cont.png");

  return;
}
