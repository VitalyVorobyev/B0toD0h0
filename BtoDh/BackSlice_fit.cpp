#include "cuts.h"
using namespace RooFit;

void BackSlice_fit(const int m_mode = 1, const bool bbflag = false){
  // mode 1 -> pi0
  // mode 2 -> eta -> gg
  // mode 3 -> eta -> pi+pi-pi0
  // mode 4 -> omega
//  const bool projection_flag = true;
//  const bool save_flag       = true;

  double BDTG_MIN = 0;
  double BDTG_MAX = 1;

  RooCategory mode("mode","mode");
  RooCategory h0mode("h0mode","h0mode");
  RooArgSet argset;

  switch(m_mode){
  case 1:
    mode.defineType("pi0",1);
    h0mode.defineType("gg",10);
    BDTG_MIN = bdtg_cut_pi0;
    break;
  case 2:
    mode.defineType("eta",2);
    h0mode.defineType("gg",10);
    BDTG_MIN = bdtg_cut_etagg;
    break;
  case 3:
    mode.defineType("eta",2);
    h0mode.defineType("ppp",20);
    BDTG_MIN = bdtg_cut_etappp;
    break;
  case 4:
    mode.defineType("omega",3);
    h0mode.defineType("ppp",20);
    BDTG_MIN = bdtg_cut_omega;
    break;
  }

  argset.add(mode); argset.add(h0mode);

  TFile *ifile;
  if(!bbflag){
    //ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_cont_0-1_full.root");
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1_full.root");
  } else{
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_bb_0-1_full.root");
  }
  TTree *tree = (TTree*)ifile->Get("TEvent");

//  string cuts;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
  argset.add(b0f);

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin  = -0.3;
  const double deMax  =  0.3;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);

//  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
//  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);
//  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  const int NSlices = 10;
  double c1_arr[NSlices], c1_arr_err[NSlices];
  double c2_arr[NSlices], c2_arr_err[NSlices];
  double mbc_arr[NSlices], mbc_arr_err[NSlices];
  double chisq_arr[NSlices];
  stringstream out;
  double dmbc = (mbcMax-mbcMin)/NSlices;

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar c1("c1","c1",mc_c1,-10.,10.);// if(cComb && false) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",0,-10.,10.);// c2.setConstant(kTRUE);
  RooChebychev pdf("pdf","pdf",de,RooArgSet(c1,c2));

  RooDataSet* ds0;
  for(int i=0; i<NSlices; i++){
    out.str("");
    out << "mbc>" << mbcMin+i*dmbc << " && mbc<" << mbcMin+(i+1)*dmbc;
    ds0 = (RooDataSet*)ds.reduce(RooArgSet(de),out.str().c_str());

    pdf.fitTo(*ds0,Verbose(),Timer(true));

    c1_arr[i] = c1.getVal(); c1_arr_err[i] = c1.getError();
    c2_arr[i] = c2.getVal(); c2_arr_err[i] = c2.getError();
    mbc_arr[i] = mbcMin+(0.5+i)*dmbc; mbc_arr_err[i] = 0;

    RooPlot* deFrame = de.frame();
    ds0->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(deFrame,LineWidth(2));

//    out.str("");
//    out << "dE, Signal " << i;
//    TCanvas* cmmbc = new TCanvas(out.str().c_str(),out.str().c_str(),600,600);
//    cmmbc->cd();

//    TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
//    pad1->Draw();

//    pad1->cd();
//    pad1->SetLeftMargin(0.15);
//    pad1->SetFillColor(0);

//    deFrame->GetXaxis()->SetTitleSize(0.05);
//    deFrame->GetXaxis()->SetTitleOffset(0.85);
//    deFrame->GetXaxis()->SetLabelSize(0.04);
//    deFrame->GetYaxis()->SetTitleOffset(1.6);
//    deFrame->Draw();

    chisq_arr[i] = deFrame->chiSquare();
  }

  for(int i=0; i<NSlices; i++){
    cout << c1_arr[i] << " " << c2_arr[i] << " " << chisq_arr[i] << endl;
  }

  TGraphErrors* gr_c1 = new TGraphErrors(NSlices,mbc_arr,c1_arr,mbc_arr_err,c1_arr_err);
  gr_c1->SetTitle("c_{1}");
  gr_c1->SetMarkerSize(1);
  gr_c1->SetMarkerColor(kBlue);
  gr_c1->SetMarkerStyle(21);
  TGraphErrors* gr_c2 = new TGraphErrors(NSlices,mbc_arr,c2_arr,mbc_arr_err,c2_arr_err);
  gr_c2->SetTitle("c_{2}");
  gr_c2->SetMarkerSize(1);
  gr_c2->SetMarkerColor(kBlue);
  gr_c2->SetMarkerStyle(21);
  TGraphErrors* gr_chisq = new TGraphErrors(NSlices,mbc_arr,chisq_arr,mbc_arr_err,mbc_arr_err);
  gr_chisq->SetTitle("Chi2");
  gr_chisq->SetMarkerSize(1);
  gr_chisq->SetMarkerColor(kBlue);
  gr_chisq->SetMarkerStyle(21);

  TCanvas* cm = new TCanvas("de","#Delta E, Cmb background",800,400);
  cm->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.99);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  gr_c1->Draw("AP");
  gr_c1->Fit("pol2");

  pad2->cd();
  pad2->SetGrid();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  gr_c2->Draw("AP");
  gr_c2->Fit("pol2");

  cm->Update();

  TCanvas* cm_chi2 = new TCanvas("chi2","#chi^{2}",400,400);
  cm_chi2->cd();
  cm_chi2->SetLeftMargin(0.15);
  cm_chi2->SetFillColor(0);
  gr_chisq->Draw("AP");

  cm_chi2->Update();

  return;
}
