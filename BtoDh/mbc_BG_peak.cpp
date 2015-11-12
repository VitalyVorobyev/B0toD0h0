#include "cuts.h"
using namespace RooFit;

void mbc_BG_peak(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig.root");
//  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/Tuples/b2dh_charged_3.root");
//  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");
  
  const int N = tree->GetEntries();
  double _mbc,_de,_md,_mk,_mpi0,_bdtgs;
  int _b0f;
  tree->SetBranchAddress("mbc",&_mbc);
  tree->SetBranchAddress("de",&_de);
  tree->SetBranchAddress("md",&_md);
  tree->SetBranchAddress("mk",&_mk);
  tree->SetBranchAddress("mpi0",&_mpi0);
  tree->SetBranchAddress("bdtgs",&_bdtgs);
  tree->SetBranchAddress("b0f",&_b0f);
  
  TTree* newtree = new TTree("newtree","newtree");
  newtree->Branch("mbc",&_mbc,"mbc/D");
//  newtree->Branch("de",&_de,"de/D");
//  newtree->Branch("md",&_md,"md/D");
//  newtree->Branch("mk",&_mk,"mk/D");
//  newtree->Branch("mpi0",&_mpi0,"mpi0/D");
//  newtree->Branch("bdtgs",&_bdtgs,"bdtgs/D");
//  newtree->Branch("b0f",&_b0f,"b0f/I");

  for(int i=0; i<N; i++){
    tree->GetEvent(i);
    //if(_b0f != 3) continue;
//     if(_bdtgs<0.98) continue;
//     if(_de<-0.15 || _de>0.3) continue;
//     if(_md<(DMass-md_cut) || _md>(DMass+md_cut)) continue;
//     if(_mk<(KMass-mk_cut) || _mk>(KMass+mk_cut)) continue;
//     if(_mpi0<(Pi0Mass-mpi0_cut) || _mpi0>(Pi0Mass+mpi0_cut)) continue;
    newtree->Fill();
  }
  newtree->Print();
  
   RooArgSet argset;
//   RooCategory b0f("b0f","b0f");
//   b0f.defineType("rho",3);
//   argset.add(b0f);
// 
//   RooRealVar de("de","#DeltaE",-0.15,0.3,"GeV");
//   argset.add(de);
//   RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
//   RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//   RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
//  RooRealVar bdtgs("bdtgs","bdtgs",0.98,1.); argset.add(bdtgs);

  RooRealVar mbc("mbc","mbc",5.26,5.29,"GeV"); argset.add(mbc);
  RooDataSet ds("ds","ds",newtree,argset);
  ds.Print();

  RooRealVar mean("mean","mean",5.278,5.27,5.29);
  RooRealVar sigma("sigma","sigma",0.005,0.001,0.5);
  RooGaussian gaus("gaus","gaus",mbc,mean,sigma);
  
  gaus.fitTo(ds,Verbose());

  /////////////
  //  Plots  //
  /////////////
  RooPlot* deFrame = mbc.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  gaus.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = mbc.frame(Title("M_{bc} pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("M_{bc} peaking BG","M_{bc} peaking BG",600,700);
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
  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
//   out.str("");
//   out << "Signal region: (" << de_min << "," << de_max << ")";
//   pt->AddText(out.str().c_str());
//   out.str("");
//   out << "Events in the SR: " << gen_back;
//   pt->AddText(out.str().c_str());
//   pt->Draw();

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
  TLine *de_lineUP = new TLine(5.20,3,5.29,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(5.20,0,5.29,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(5.20,-3,5.29,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
}
