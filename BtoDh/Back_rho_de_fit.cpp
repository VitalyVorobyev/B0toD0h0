#include "cuts.h"
using namespace RooFit;

void Back_rho_de_fit(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_gen.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  RooCategory b0f("b0f","b0f");
//  b0f.defineType("comb",-1);
  b0f.defineType("rho",3);
//  b0f.defineType("unknown",0);
  
  RooCategory d0f("d0f","d0f");
  d0f.defineType("true",1);
//  d0f.defineType("fake",-1);

  RooArgSet argset;
  RooRealVar mbc("mbc","mbc",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.15,0.3,"GeV"); argset.add(de);
  de.setRange("signal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
//  RooRealVar bdt("bdt","bdt",bdt_cut,1.); argset.add(bdt);
  RooRealVar bdtgs("bdtgs","bdtgs",0.98,1.); argset.add(bdtgs);
//  RooRealVar bdtg("bdtg","bdtg",0.95,1.); argset.add(bdtg);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  argset.add(b0f);
//  argset.add(d0f);
  
  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();
  
//   TFile kestfile("RhoHist.root","RECREATE");
//   RooKeysPdf kest("kest","kest",de,ds,RooKeysPdf::MirrorBoth,1);
// //  RooKeysPdf kest1("kest1","kest1",de,ds,RooKeysPdf::NoMirror);
//   de.setBins(1000);
//   RooDataHist* hist = new RooDataHist("hist","hist",RooArgSet(de));
//   hist->Print();
//   RooArgSet* deset = new RooArgSet(de);
//   deset->Print();
//   kest.fillDataHist(hist,deset,1);
//   hist->Print();
//   RooHistPdf histpdf("histpdf","histpdf",RooArgSet(de),*hist);//????
//   histpdf.Write();
//   kestfile.Close();
//   de.setBins(60);
//   ifile->cd();
  
  stringstream out;
  out.str("");
  out << "de<" << de_max << "&&de>" << de_min;
  const int gen_back = ds.sumEntries(out.str().c_str());
  cout << "Gen background: " << gen_back << endl;

  RooRealVar atanpar1("atanpar1","atanpar1",-200,-1000.,-80.); if(cRho) atanpar1.setConstant(kTRUE);
  RooRealVar atanpar2("atanpar2","atanpar2",-23.,-35.,-10.); if(cRho) atanpar2.setConstant(kTRUE);
  RooGenericPdf pdf("pdf","TMath::ATan(@0*@1+@2)+0.5*TMath::Pi()",RooArgSet(de,atanpar1,atanpar2));

//   RooRealVar atanpar1("atanpar1","atanpar1",85,80.,100.); if(cRho) atanpar1.setConstant(kTRUE);
//   RooRealVar atanpar2("atanpar2","atanpar2",-11.,-12.,-10.); if(cRho) atanpar2.setConstant(kTRUE);
//   RooGenericPdf pdf("pdf","TMath::Erfc(@0*@1-@2)",RooArgSet(de,atanpar1,atanpar2));
  
  pdf.fitTo(ds,Verbose(),Timer(true));
  
  /////////////
  //  Plots  //
  /////////////
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,LineWidth(2));
//  kest1.plotOn(deFrame,LineWidth(2),LineColor(kMagenta));

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

  TPaveText *pt = new TPaveText(0.5,0.7,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
  out.str("");
  out << "Signal region: (" << de_min << "," << de_max << ")";
  pt->AddText(out.str().c_str());
  out.str("");
  out << "Events in the SR: " << gen_back;
  pt->AddText(out.str().c_str());
  pt->Draw();

  TLine *de_lineLEFT = new TLine(de_min,0,de_min,1.5*gen_back);
  de_lineLEFT->SetLineColor(kRed);
  de_lineLEFT->SetLineStyle(1);
  de_lineLEFT->Draw();
  
  TLine *de_lineRIGHT = new TLine(de_max,0,de_max,1.5*gen_back);
  de_lineRIGHT->SetLineColor(kRed);
  de_lineRIGHT->SetLineStyle(1);
  de_lineRIGHT->Draw();

  TLine *de_lineFITEDGE = new TLine(-0.15,0,-0.15,400);
  de_lineFITEDGE->SetLineColor(kGreen);
  de_lineFITEDGE->SetLineStyle(1);
  de_lineFITEDGE->Draw();
  
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
