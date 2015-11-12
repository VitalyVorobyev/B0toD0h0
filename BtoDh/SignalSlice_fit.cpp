#include "cuts.h"
using namespace RooFit;

void SignalSlice_fit(const int m_mode = 1, const int _b0f=-1, const bool SF = false){
  // mode 1 -> pi0
  // mode 2 -> eta
  // mode 3 -> omega
//  const bool projection_flag = true;
//  const bool save_flag       = true;
  const bool bg_fit          = false;

  double de_sig_min,de_sig_max;
  double BDTG_MIN = 0;
  double BDTG_MAX = 1;

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  double deMin  = SF ? -0.1 : -0.25;
  double deMax  = -deMin;
  if(_b0f == 3){
    deMin = -0.15;
    deMax = -0.10;
  }

  int _mode,_h0mode;
  bool gg_flag = true;
  TFile *ifile;
  switch(m_mode){
    case 1:
      _mode = 1;
      _h0mode = 10;
      de_sig_min = de_min;
      de_sig_max = de_max;
      BDTG_MIN = bdtg_cut_pi0;
//      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_full.root");
      if(_b0f != 3) ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
      else          ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-2.root");
      break;
    case 2:
      _mode = 2;
      _h0mode = 10;
      de_sig_min = de_min;      de_sig_max = de_max;
//      BDTG_MIN = bdtg_cut_etagg;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
      break;
    case 3:
      _mode = 2;
      _h0mode = 20;
      gg_flag = false;
      de_sig_min = de_min_etappp;
      de_sig_max = de_max_etappp;
//      BDTG_MIN = bdtg_cut_etappp;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
      break;
    case 4:
      _mode = 3;
      _h0mode = 20;
      gg_flag = false;
      de_sig_min = de_min_omega;
      de_sig_max = de_max_omega;
//      BDTG_MIN = bdtg_cut_omega;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s5_full.root");
      break;
    default:
      return;
  }

  bool remove_right_CB_flag  = false;
  if(_b0f == 5) remove_right_CB_flag = true;

  TTree *tree = (TTree*)ifile->Get("TEvent");
  RooArgSet argset;

//  string cuts;
  RooCategory b0f("b0f","b0f");
  if(_b0f == 1 || _b0f == -1) b0f.defineType("signal",1);
  if(_b0f == 1 || _b0f == -1) b0f.defineType("fsr",10);
  if(_b0f == 5 || _b0f == -1){
    b0f.defineType("bad_pi0",5);
//    b0f.defineType("bad_pi0_m11",-11);
    b0f.defineType("bad_pi0_4",4);
  }
  if(_b0f == 3) b0f.defineType("rho",3);

  argset.add(b0f);

  RooCategory mode("mode","mode");
  mode.defineType("sig",_mode);
  if(_b0f != 3) argset.add(mode);

  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("h0sig",_h0mode);
  if(_b0f != 3) argset.add(h0mode);

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good",1);
  if(_b0f != 3) argset.add(good_icpv);

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_sig_min,de_sig_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_sig_min,de_sig_max);

  if(_b0f != 3){ RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);}
  if(_b0f != 3){ RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);}
  if(_b0f != 3){ RooRealVar chi2_vtx_d0("chi2_vtx_d0","chi2_vtx_d0",0,50); argset.add(chi2_vtx_d0);}
  if(_b0f != 3){ RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);}
//  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  const int NSlices = _b0f == 3 ? 10 : 30;
  double de_arr[NSlices], de_arr_err[NSlices];
  double chisq_arr[NSlices];
  double mbc0_arr[NSlices], mbc0_arr_err[NSlices];

  double sl_arr[NSlices], sl_arr_err[NSlices];
  double sr_arr[NSlices], sr_arr_err[NSlices];
  double s_arr[NSlices],  s_arr_err[NSlices];
//  double sll_arr[NSlices], sll_arr_err[NSlices];
//  double srr_arr[NSlices], srr_arr_err[NSlices];
//  double mbc00_arr[NSlices], mbc00_arr_err[NSlices];
//  double alpha_arr[NSlices], alpha_arr_err[NSlices];
//  int N_arr[NSlices];
  stringstream out;
  double dde = (deMax-deMin)/NSlices;

  RooDataSet ds("ds","ds",tree,argset);
//  RooDataHist* dh = ds0->binnedClone();  
  ds.Print();
//  ds0->Print();

  /////////////
  // mbc pdf //
  /////////////
  ///  mbc0
  if(bg_fit){
    RooRealVar mbc0("mbc0","mbc0",5.27940,5.27,5.29);
    RooRealVar sl("sl","sl",0.00255,0.,0.05);
    RooRealVar sr("sr","sr",0.002,0.,0.5);
    RooBifurGauss pdf("pdf","pdf",mbc,mbc0,sl,sr);
  } else{
    RooRealVar mbc0("mbc0","mbc0",5.27940,5.27,5.29);
    RooRealVar s("s","s",0.00255,0.,0.05);
    RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
    RooNovosibirsk pdf("pdf","pdf",mbc,mbc0,s,alpha);
  }
//  RooGaussian g("g","g",mbc,mbc0,sr);
//  RooRealVar mbc00("mbc00","mbc00",5.27940,5.27,5.29);
//  RooRealVar sll("sll","sll",0.006,0.001,0.5);
//  RooRealVar srr("srr","srr",0.004,0.,0.5);
//  RooBifurGauss pdf("pdf","pdf",mbc,mbc0,sl,sr);
//  RooRealVar deCBl("deCBl","deCBl",get_deCBl(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG) deCBl.setConstant(kTRUE);
//  RooRealVar sCBl("sCBl","sCBl",,0.,0.5); if(cSIG) sCBl.setConstant(kTRUE);
//  RooRealVar n("n","n",1.,0.,100.);
//  RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
//  RooRealVar alpha("alpha","alpha",0.2,0.01,2.); alpha.setConstant(kTRUE);
//  RooCBShape CB("CB","CB",mbc,mbc00,sll,alpha,n);
//  RooRealVar f("f","f",0.1,0.,1.);// f.setConstant(kTRUE);
//  RooRealVar fCB("fCB","fCB",0.1,0.,1.);// f.setConstant(kTRUE);
//  RooNovosibirsk pdf("pdf","pdf",mbc,mbc0,sl,alpha);
//  RooAddPdf pdf("pdf","pdf",RooArgSet(Sib,bgg),RooArgList(f));

  RooDataSet* ds0;
  for(int i=0; i<NSlices; i++){
    out.str("");
    out << "de>" << deMin+i*dde << " && de<" << deMin+(i+1)*dde;
    de_arr[i] = deMin+(0.5+i)*dde; de_arr_err[i] = 0;

    ds0 = (RooDataSet*)ds.reduce(RooArgSet(mbc),out.str().c_str());
    pdf.fitTo(*ds0,Verbose(),Timer(true));
    if(bg_fit){
      sl_arr[i] = sl.getVal(); sl_arr_err[i] = sl.getError();
      sr_arr[i] = sr.getVal(); sr_arr_err[i] = sr.getError();
      mbc0_arr[i] = mbc0.getVal(); mbc0_arr_err[i] = mbc0.getError();
    } else{
      mbc0_arr[i] = mbc0.getVal(); mbc0_arr_err[i] = mbc0.getError();
      s_arr[i]    = s.getVal();    s_arr_err[i]    = s.getError();
    }

//    sll_arr[i] = sll.getVal(); sll_arr_err[i] = sll.getError();
//    srr_arr[i] = srr.getVal(); srr_arr_err[i] = srr.getError();
//    mbc00_arr[i] = mbc00.getVal(); mbc00_arr_err[i] = mbc00.getError();
//    alpha_arr[i] = alpha.getVal(); alpha_arr_err[i] = alpha.getError();

    RooPlot* mbcFrame = mbc.frame();
    ds0->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(mbcFrame,LineWidth(2));

    if(_b0f == 3){
    out.str("");
    out << "M_{bc}, Signal " << i;
    TCanvas* cmmbc = new TCanvas(out.str().c_str(),out.str().c_str(),600,600);
    cmmbc->cd();

    TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
    pad1->Draw();

    pad1->cd();
    pad1->SetLeftMargin(0.15);
    pad1->SetFillColor(0);

    mbcFrame->GetXaxis()->SetTitleSize(0.05);
    mbcFrame->GetXaxis()->SetTitleOffset(0.85);
    mbcFrame->GetXaxis()->SetLabelSize(0.04);
    mbcFrame->GetYaxis()->SetTitleOffset(1.6);
    mbcFrame->Draw();
    }

    chisq_arr[i] = mbcFrame->chiSquare();
  }

  for(int i=0; i<NSlices; i++){
    if(bg_fit) cout << mbc0_arr[i] << " " << sl_arr[i] << " " << sr_arr[i] << " " << chisq_arr[i] << endl;
    else       cout << mbc0_arr[i] << " " << s_arr[i] << " " << " " << chisq_arr[i] << endl;
  }

  TGraphErrors* gr_mbc0 = new TGraphErrors(NSlices,de_arr,mbc0_arr,de_arr_err,mbc0_arr_err);
  gr_mbc0->SetTitle("M^{0}_{bc}");
  gr_mbc0->SetMarkerSize(1);
  gr_mbc0->SetMarkerColor(kBlue);
  gr_mbc0->SetMarkerStyle(21);
//  TGraphErrors* gr_mbc00 = new TGraphErrors(NSlices,de_arr,mbc00_arr,de_arr_err,mbc00_arr_err);
//  gr_mbc00->SetMarkerSize(1);
//  gr_mbc00->SetMarkerColor(kRed);
//  gr_mbc00->SetMarkerStyle(21);
//  TMultiGraph* mg_mbc0 = new TMultiGraph();
//  mg_mbc0->SetTitle("Mbc0");
//  mg_mbc0->Add(gr_mbc0);
//  mg_mbc0->Add(gr_mbc00);

//  TGraphErrors* gr_alpha = new TGraphErrors(NSlices,de_arr,alpha_arr,de_arr_err,alpha_arr_err);
//  gr_alpha->SetTitle("#alpha");
//  gr_alpha->SetMarkerSize(1);
//  gr_alpha->SetMarkerColor(kBlue);
//  gr_alpha->SetMarkerStyle(21);

  if(bg_fit){
    TGraphErrors* gr_sl = new TGraphErrors(NSlices,de_arr,sl_arr,de_arr_err,sl_arr_err);
    gr_sl->SetTitle("#sigma_{l}");
//    gr_sl->SetTitleSize(40);
    gr_sl->SetMarkerSize(1);
    gr_sl->SetMarkerColor(kBlue);
    gr_sl->SetMarkerStyle(21);
    TGraphErrors* gr_sr = new TGraphErrors(NSlices,de_arr,sr_arr,de_arr_err,sr_arr_err);
    gr_sr->SetTitle("#sigma_{r}");
//    gr_sr->SetTitleSize(40);
    gr_sr->SetMarkerSize(1);
    gr_sr->SetMarkerColor(kRed);
    gr_sr->SetMarkerStyle(21);
    TMultiGraph* mg_s = new TMultiGraph();
    mg_s->SetTitle("Width");
//    mg_s->SetTitleSize(40);
    mg_s->Add(gr_sl);
    mg_s->Add(gr_sr);
  } else{
    TGraphErrors* gr_s = new TGraphErrors(NSlices,de_arr,s_arr,de_arr_err,s_arr_err);
    gr_s->SetTitle("#sigma_{Mbc}");
//    gr_s->SetTitleSize(40);
    gr_s->SetMarkerSize(1);
    gr_s->SetMarkerColor(kBlue);
    gr_s->SetMarkerStyle(21);
  }

//  TGraphErrors* gr_sll = new TGraphErrors(NSlices,de_arr,sll_arr,de_arr_err,sll_arr_err);
//  gr_sll->SetMarkerSize(1);
//  gr_sll->SetMarkerColor(kBlue);
//  gr_sll->SetMarkerStyle(21);
//  TGraphErrors* gr_srr = new TGraphErrors(NSlices,de_arr,srr_arr,de_arr_err,srr_arr_err);
//  gr_srr->SetMarkerSize(1);
//  gr_srr->SetMarkerColor(kRed);
//  gr_srr->SetMarkerStyle(21);
//  TMultiGraph* mg_ss = new TMultiGraph();
//  mg_ss->SetTitle("Width 2");
//  mg_ss->Add(gr_sll);
//  mg_ss->Add(gr_srr);

  TGraphErrors* gr_chisq = new TGraphErrors(NSlices,de_arr,chisq_arr,de_arr_err,de_arr_err);
  gr_chisq->SetTitle("Chi2");
  gr_chisq->SetMarkerSize(1);
  gr_chisq->SetMarkerColor(kBlue);
  gr_chisq->SetMarkerStyle(21);

  TCanvas* cm = new TCanvas("de_slice","#Delta E, Signal",800,400);
  cm->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.99);
//  TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.49);
//  TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.49);
  pad1->Draw();
  pad2->Draw();
//  pad3->Draw();
//  pad4->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetGrid();
  gr_mbc0->GetXaxis()->SetTitle("#DeltaE (GeV)");
  gr_mbc0->GetXaxis()->SetTitleSize(0.05);
  gr_mbc0->GetXaxis()->SetLabelSize(0.05);
  gr_mbc0->GetYaxis()->SetLabelSize(0.05);
  gr_mbc0->GetXaxis()->SetRangeUser(-0.3,0.3);
  gr_mbc0->Draw("AP");
  if(!(m_mode>2 && _b0f == 1)){ gr_mbc0->Fit("pol2");}
  else{
    TF1 *vsfit = new TF1("vsfit","[0]*TMath::Erf(([1]-x)/[2])+[3]",deMin,deMax);
    vsfit->SetParameters(-1,0,1,5.28);
    gr_mbc0->Fit(vsfit,"qr");
    double vspar[4];
    vsfit->GetParameters(vspar);
    cout << "Par: " << vspar[0] << " " << vspar[1] << " " << vspar[2] << " " << vspar[3] << endl;
  }
  pad1->Update();

  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  pad2->SetGrid();
  if(bg_fit) mg_s->Draw("AP");
  else{
    gr_s->GetXaxis()->SetTitle("#DeltaE (GeV)");
    gr_s->GetXaxis()->SetTitleSize(0.05);
    gr_s->GetXaxis()->SetLabelSize(0.05);
    gr_s->GetYaxis()->SetLabelSize(0.05);
    gr_s->GetXaxis()->SetRangeUser(-0.3,0.3);
    gr_s->Draw("AP");
    gr_s->Fit("pol2");
  }

  cm->Update();

//  pad3->cd();
//  pad3->SetLeftMargin(0.15);
//  pad3->SetFillColor(0);
//  pad3->SetGrid();
//  gr_alpha->Draw("AP");

  TCanvas* c2 = new TCanvas("de_chi2","#Delta E, Signal",400,400);
  c2->cd();
  c2->SetGrid();
  gr_chisq->Draw("AP");
  c2->Update();

  return;
}
