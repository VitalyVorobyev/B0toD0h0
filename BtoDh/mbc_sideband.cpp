#include "cuts.h"
using namespace RooFit;

void mbc_sideband(){
  TChain* tree_gen = new TChain("TEvent");
  tree_gen->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-2.root");
  tree_gen->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_3-5.root");

  TChain* tree_data = new TChain("TEvent");
  tree_data->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_data.root");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
//  b0f.defineType("rho",3);

  stringstream out;
  RooArgSet argset;

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin = -0.3;
  const double deMax = 0.3;

  RooRealVar mbc("mbc","M_{bc}",0.5*(mbc_min+mbc_max),mbcMin,mbcMax,"GeV");
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("Sideband",mbcMin,5.25);
  argset.add(mbc);

  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV");
  de.setRange("Signal",deMin,deMax);
  de.setRange("Sideband",deMin,deMax);
  argset.add(de);

  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtgs("bdtgs","bdtgs",bdtgs_cut,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  RooDataSet dsdat("dsdat","dsdat",tree_data,argset,"mbc>0||mbc<=0");
  
  RooDataSet* dsdatsb = dsdat.reduce(CutRange("Sideband"));
  dsdatsb->Print();
  
  argset.add(b0f);
  RooDataSet dsgen("dsgen","dsgen",tree_gen,argset,"mbc>0||mbc<=0");
  
  RooDataSet* dsgensb = dsgen.reduce(CutRange("Sideband"));
  dsgensb->Print();
  RooDataSet* dsgensr = dsgen.reduce(CutRange("Signal"));
  dsgensr->Print();

  /////////
  // pdf //
  /////////
  RooRealVar c1("c1","c1",mc_c1,-10.,10.);
  RooRealVar c2("c2","c2",mc_c2,-10.,10.);
  RooChebychev pdf("pdf","pdf",de,RooArgSet(c1,c2));
  
  pdf.fitTo(*dsdatsb);
  pdf.fitTo(*dsgensr);
  pdf.fitTo(*dsgensb,Timer(true));
   
  /////////////
  //  Plots  //
  /////////////
  // Signal //
  RooPlot* deFrameSig = de.frame();
  dsgensr->plotOn(deFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrameSig,LineWidth(2));

  RooHist* hdepull = deFrameSig->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E Generic SR","#Delta E Generic SR",600,700);
  cm->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrameSig->GetXaxis()->SetTitleSize(0.05);
  deFrameSig->GetXaxis()->SetTitleOffset(0.85);
  deFrameSig->GetXaxis()->SetLabelSize(0.04);
  deFrameSig->GetYaxis()->SetTitleOffset(1.6);
  deFrameSig->Draw();
  
  TPaveText *ptmbc = new TPaveText(0.6,0.83,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrameSig->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  TLine *de_line_RIGHT;
  de_line_RIGHT = new TLine(de_max,0,de_max,30);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  de_line_LEFT = new TLine(de_min,0,de_min,30);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(deMin,3,deMax,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(deMin,0,deMax,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(deMin,-3,deMax,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
  
  
  // Generic sideband //
  RooPlot* deFrameGenSB = de.frame();
  dsgensb->plotOn(deFrameGenSB,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrameGenSB,LineWidth(2));

  RooHist* hdepull = deFrameGenSB->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmGCB = new TCanvas("#Delta E Generic SB","#Delta E Generic SB",600,700);
  cmGCB->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrameGenSB->GetXaxis()->SetTitleSize(0.05);
  deFrameGenSB->GetXaxis()->SetTitleOffset(0.85);
  deFrameGenSB->GetXaxis()->SetLabelSize(0.04);
  deFrameGenSB->GetYaxis()->SetTitleOffset(1.6);
  deFrameGenSB->Draw();

  TLine *de_line_RIGHT;
  de_line_RIGHT = new TLine(de_max,0,de_max,30);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  de_line_LEFT = new TLine(de_min,0,de_min,30);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  TPaveText *ptmbc = new TPaveText(0.6,0.83,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrameGenSB->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(deMin,3,deMax,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(deMin,0,deMax,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(deMin,-3,deMax,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cmGCB->Update();
  
  // Data sidedbnd //
  RooPlot* deFrameData = de.frame();
  dsdatsb->plotOn(deFrameData,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrameData,LineWidth(2));

  RooHist* hdepull = deFrameData->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmD = new TCanvas("#Delta E Data SB","#Delta E Data SB",600,700);
  cmD->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrameData->GetXaxis()->SetTitleSize(0.05);
  deFrameData->GetXaxis()->SetTitleOffset(0.85);
  deFrameData->GetXaxis()->SetLabelSize(0.04);
  deFrameData->GetYaxis()->SetTitleOffset(1.6);
  deFrameData->Draw();

  TLine *de_line_RIGHT;
  de_line_RIGHT = new TLine(de_max,0,de_max,30);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  de_line_LEFT = new TLine(de_min,0,de_min,30);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();
  
  TPaveText *ptmbc = new TPaveText(0.6,0.83,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrameData->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(deMin,3,deMax,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(deMin,0,deMax,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(deMin,-3,deMax,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cmD->Update();
}
