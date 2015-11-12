#include "cuts.h"
using namespace RooFit;

void BB_comb_2d_fit(const int _mode = 4){
  TChain* tree = new TChain("TEvent");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_bb_0-1_full.root");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
  b0f.defineType("peak1",3);
  b0f.defineType("peak2",4);
  b0f.defineType("peak3",11);
  b0f.defineType("peak4",20);
  RooCategory d0f("d0f","d0f");
  d0f.defineType("true",1);

  RooCategory mode("mode","mode");
  RooCategory h0mode("h0mode","h0mode");

  string MODE,Mode;
  double mbcMin = 5.2;
//  if(_mode == 3 || _mode == 4) mbcMin = 5.215;
  const double mbcMax = 5.29;
  const double deMin = -0.15;
  const double deMax = 0.3;
  double BDTG_MIN = 0;
  double BDTG_MAX = 1;
  switch(_mode){
  case 1: BDTG_MIN = bdtg_cut_pi0;
          mode.defineType("pi0",1);
          h0mode.defineType("gg",10);
          MODE = string("pi0");
          Mode = string("#pi0");
          break;
  case 2: BDTG_MIN = bdtg_cut_etagg;
          mode.defineType("eta",2);
          h0mode.defineType("gg",10);
          MODE = string("etagg");
          Mode = string("#eta#to#gamma#gamma");
          break;
  case 3: BDTG_MIN = bdtg_cut_etappp;
          mode.defineType("eta",2);
          h0mode.defineType("ppp",20);
          MODE = string("etappp");
          Mode = string("#eta#to#pi^{+}#pi^{-}#pi^{0}");
          break;
  case 4: BDTG_MAX = bdtg_cut_omega;
          mode.defineType("omega",3);
          h0mode.defineType("ppp",20);
          MODE = string("omega");
          Mode = string("#omega");
          break;
  default:
          return;
  }

  RooArgSet argset;
  argset.add(mode);
  argset.add(h0mode);

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
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut);// argset.add(atckpi_max);

  argset.add(b0f); argset.add(d0f);

  RooDataSet ds("ds","ds",tree,argset,"(mbc>0||mbc<=0) && (de>0||de<=0)");// to eliminate NANs

  ds.Print();

  /////////////////
  // BB comb PDF //
  /////////////////
  /////////////
  // mbc pdf //
  /////////////
  RooRealVar edge("edge","edge",5.29,5.28,5.30,"GeV");// edge.setConstant(kTRUE);
  RooRealVar mbctau("mbctau","mbctau",-4.,-100,0.,"GeV");
  RooRealVar pow1("pow1","pow1",1.,0.1,10);
  RooRealVar pow2("pow2","pow2",1.,0.1,10);

  RooGenericPdf pdf_mbc_comb("pdf_mbc_comb","pow(@1-@0,@2)*exp(@3*pow(@1-@0,@4))",RooArgList(mbc,edge,pow1,mbctau,pow2));

  ////////////
  // de pdf //
  ////////////
  RooRealVar c3("c3","c3",-2.221e+04,-2.4e4,-2.0e4);// c3.setConstant(kTRUE);
  RooRealVar c31("c31","c31",8496,8000.,10000);     // c31.setConstant(kTRUE);
  RooRealVar c32("c32","c32",-812.5,-900.,-700.);   // c32.setConstant(kTRUE);
  RooFormulaVar _c3("_c3","@0+@1*@3+@2*@3*@3",RooArgSet(c3,c31,c32,mbc));
  RooExponential pdf_de_comb("pdf_de_comb","pdf_de_comb",de,_c3);

  RooProdPdf pdf("pdf","pdf",pdf_mbc_comb,Conditional(pdf_de_comb,de));

  RooRealVar de0("de0","de0",0.,-0.1,0.1,"GeV");
  RooRealVar deWidth("deWidth","deWidth",0.01,0.,0.02,"GeV");
  RooGaussian deGaus("deGaus","deGaus",de,de0,deWidth);

  RooRealVar mbc0("mbc0","mbc0",5.275,5.27,5.285,"GeV");
  RooRealVar mbcWidth("mbcWidth","mbcWidth",0.004,0.,0.1,"GeV");
  RooGaussian mbcGaus("mbcGaus","mbcGaus",mbc,mbc0,mbcWidth);

  RooProdPdf pdfPeak("pdfPeak","pdfPeak",RooArgList(deGaus,mbcGaus));

  RooRealVar NPeak("NPeak","NPeak",0.,0.,10000.);
  RooRealVar NBack("NBack","NBack",20000,0.,100000.);

//  RooAddPdf pdf("pdf","pdf",RooArgSet(pdfBack,pdfPeak),RooArgList(NBack,NPeak));

//  pdf_mbc_comb.fitTo(ds,Verbose(),Timer(true));
  pdf.fitTo(ds,Verbose(),Timer(true));

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
  TPaveText *pt = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "mode: " << Mode;
  pt->AddText(out.str().c_str());
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,120);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,120);
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
  out.str("");
  out << "../Reports/pics/de_2dfit_comb_" << MODE << ".png";
  cm->Print(out.str().c_str());
  out.str("");
  out << "../Reports/pics/de_2dfit_comb_" << MODE << ".root";
  cm->Print(out.str().c_str());
  
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

  TPaveText *ptmbc = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "mode: " << Mode;
  pt->AddText(out.str().c_str());
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,70);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,70);
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

  cmmbc->Update();
  out.str("");
  out << "../Reports/pics/mbc_bb_2dfit_comb_" << MODE << ".png";
  cmmbc->Print(out.str().c_str());
  out.str("");
  out << "../Reports/pics/mbc_bb_2dfit_comb_" << MODE << ".root";
  cmmbc->Print(out.str().c_str());
 
  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.15,0.1),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");
//  Chi2->Print();
  return;
}

