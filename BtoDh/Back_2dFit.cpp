#include "cuts.h"
using namespace RooFit;

void Back_2dFit(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_gen.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);

  RooArgSet argset;
  
  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin = -0.15;
  const double deMax = 0.3;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtgs("bdtgs","bdtgs",bdtgs_cut,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  argset.add(b0f);

  RooDataSet ds("ds","ds",tree,argset);
  RooDataSet* ds0 = ds.reduce(RooArgSet(de,mbc));
  RooDataHist* dh = ds0->binnedClone();
  ds.Print();
  dh->Print();

  ////////////
  // de pdf //
  ////////////
//  RooRealVar c1("c1","c1",mc_c1,-1.,0.); if(cComb) c1.setConstant(kTRUE);
//  RooRealVar c2("c2","c2",mc_c2,0.,0.1); if(cComb) c2.setConstant(kTRUE);
  RooRealVar c3("c3","c3",mc_c1,-1.,0.); if(cComb) c3.setConstant(kTRUE);
  RooRealVar c4("c4","c4",mc_c2,0.,0.1); if(cComb) c4.setConstant(kTRUE);
  RooChebychev pdf_de("pdf_de","pdf_de",de,RooArgSet(c3,c4));

  /////////////
  // mbc pdf //
  /////////////
  RooRealVar argpar("argpar","argus shape parameter",mc_argpar,-100.,-1.); if(cComb)
  argpar.setConstant(kTRUE);
  RooRealVar argedge("argedge","argedge",mc_argedge,5.285,5.292); if(cComb) argedge.setConstant(kTRUE);
  RooArgusBG pdf_mbc("pdf_mbc","Argus PDF",mbc,argedge,argpar);

  /////////
  // pdf //
  /////////
  RooProdPdf pdf("pdf","pdf",RooArgList(pdf_de,pdf_mbc));

  pdf.fitTo(*dh,Verbose(),Timer(true));

  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf.plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  pdf.plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

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
  deFrame->GetXaxis()->SetLabelSize(0.04);
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

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,200);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,200);
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
  
  // mbc //
  RooPlot* mbcFrame = mbc.frame();
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf.plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc.frame(Title("#Delta E pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("M_{bc}, Signal","M_{bc}, Signal",600,700);
  cmmbc->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);

  mbcFrame->GetXaxis()->SetTitleSize(0.05);
  mbcFrame->GetXaxis()->SetTitleOffset(0.85);
  mbcFrame->GetXaxis()->SetLabelSize(0.04);
  mbcFrame->GetYaxis()->SetTitleOffset(1.6);
  mbcFrame->Draw();

  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,150);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,150);
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();
  
  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  mbcPull->SetMarkerSize(0.05);
  mbcPull->Draw();
  TLine *mbc_lineUP = new TLine(mbcMin,3,mbcMax,3);
  mbc_lineUP->SetLineColor(kBlue);
  mbc_lineUP->SetLineStyle(2);
  mbc_lineUP->Draw();
  TLine *mbc_line = new TLine(mbcMin,0,mbcMax,0);
  mbc_line->SetLineColor(kBlue);
  mbc_line->SetLineStyle(1);
  mbc_line->SetLineWidth((Width_t)2.);
  mbc_line->Draw();
  TLine *mbc_lineDOWN = new TLine(mbcMin,-3,mbcMax,-3);
  mbc_lineDOWN->SetLineColor(kBlue);
  mbc_lineDOWN->SetLineStyle(2);
  mbc_lineDOWN->Draw();
  
  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.15,0.1),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");

  cmmbc->Update();
}
