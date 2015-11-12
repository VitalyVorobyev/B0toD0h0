#include "cuts.h"
using namespace RooFit;

void Signal_de_fit(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  b0f.defineType("fsr",10);
  b0f.defineType("bad_pi0",5);

  RooArgSet argset;

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;

  RooRealVar mbc("mbc","mbc",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.3,0.3,"GeV"); argset.add(de);
  de.setRange("signal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtgs("bdtgs","bdtgs",0.98,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  argset.add(b0f);

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  RooRealVar de0("de0","de0",0.018,-0.1,0.1);
  RooRealVar s1("s1","s1",0.021,0.,0.5);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar deCBl("deCBl","deCBl",-0.014,-0.1,0.1);
  RooRealVar sCBl("sCBl","sCBl",0.032,0.,0.5);
  RooRealVar nl("nl","nl",12.1,0.,100.);
  RooRealVar alphal("alphal","alphal",0.72,-10.,10.);

  RooRealVar deCBr("deCBr","deCBr",-0.015,-0.1,0.1);
  RooRealVar sCBr("sCBr","sCBr",0.066,0.,0.5);
  RooRealVar nr("nr","nr",19.7,0.,100.);
  RooRealVar alphar("alphar","alphar",-0.75,-10.,10.);
  
  RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);
  
  RooRealVar fCBl("fCBl","fCBl",0.20,0.,1.);
  RooRealVar fCBr("fCBr","fCBr",0.20,0.,1.);
  
  RooAddPdf pdf("pdf","pdf",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
//  RooAddPdf pdf("pdf","pdf",RooArgList(CBr,g1),RooArgSet(fCBr));
  
  pdf.fitTo(ds,Verbose(),Timer(true));
  
  /////////////
  //  Plots  //
  /////////////
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
//  pdf.plotOn(deFrame,Components(CBl),LineStyle(kDashed));
//  pdf.plotOn(deFrame,Components(CBr),LineStyle(kDashed));
//  pdf.plotOn(deFrame,Components(g1),LineStyle(kDashed));
  pdf.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E, Signal","#Delta E, Signal",600,700);
  cm->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.05);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  TPaveText *pt = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(out.str().c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,1500);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,1500);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();
  
  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(-0.3,3,0.3,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(-0.3,0,0.3,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(-0.3,-3,0.3,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
}
