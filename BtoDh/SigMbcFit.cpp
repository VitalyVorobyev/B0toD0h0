#include "cuts.h"
using namespace RooFit;

void SigMbcFit(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");
  
  const double mbc_min = 5.26;
  const double mbc_max = 5.29;
  
  RooCategory b0f("b0f","b0f");
  b0f.defineType("good",1);
  b0f.defineType("bad",5);
  b0f.defineType("fsr",10);

  RooArgSet argset;
  RooRealVar mbc("mbc","M_{bc}",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.3,0.3,"GeV");
  argset.add(de);
//   de.setRange("signal",de_min,de_max);
//   de.setRange("fit",-0.15,0.3);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtgs("bdtgs","bdtgs",0.98,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  argset.add(b0f);

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  RooRealVar mbc0("mbc0","mbc0",5.2846,5.26,5.30);
  RooRealVar sl("sl","sl",0.011,0.,0.5);
  RooRealVar sr("sr","sr",0.0016,0.,0.5);
  RooBifurGauss bg("bg","bg",mbc,mbc0,sl,sr);

  RooRealVar mbc00("mbc00","mbc00",5.280,5.26,5.30);
  RooRealVar sll("sll","sll",0.0032,0.,0.5);
  RooRealVar srr("srr","srr",0.0023,0.,0.5);
  RooBifurGauss bgg("bgg","bgg",mbc,mbc00,sll,srr);
 
  RooRealVar f("f","f",0.2,0.,1.);   
  RooAddPdf pdf("pdf","pdf",RooArgList(bg,bgg),RooArgSet(f));

  stringstream out;  
//   out.str("");
//   out << "de<" << de_max << "&&de>" << de_min;
//   const int gen_back = ds.sumEntries(out.str().c_str());
//   cout << "Gen background: " << gen_back << endl;
  
  pdf.fitTo(ds,Verbose(),Timer(true));
  
  /////////////
  //  Plots  //
  /////////////
  RooPlot* deFrame = mbc.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(bgg),LineWidth(1),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(bg),LineWidth(1),LineStyle(kDashed));
  pdf.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = mbc.frame(Title("#Delta E pull distribution"));
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

  TPaveText *pt = new TPaveText(0.55,0.87,0.95,0.93,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
//   out.str("");
//   out << "Signal region: (" << de_min << "," << de_max << ")";
//   pt->AddText(out.str().c_str());
//   out.str("");
//   out << "Events in the SR: " << gen_back;
//   pt->AddText(out.str().c_str());
  pt->Draw();

//   TLine *de_lineLEFT = new TLine(de_min,0,de_min,0.3*gen_back);
// //  TLine *de_lineLEFT = new TLine(de_min,0,de_min,18);
//   de_lineLEFT->SetLineColor(kRed);
//   de_lineLEFT->SetLineStyle(1);
//   de_lineLEFT->Draw();
//   
//   TLine *de_lineRIGHT = new TLine(de_max,0,de_max,0.3*gen_back);
// //  TLine *de_lineRIGHT = new TLine(de_max,0,de_max,18);
//   de_lineRIGHT->SetLineColor(kRed);
//   de_lineRIGHT->SetLineStyle(1);
//   de_lineRIGHT->Draw();
  
  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(mbc_min,3,mbc_max,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(mbc_min,0,mbc_max,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(mbc_min,-3,mbc_max,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
}
