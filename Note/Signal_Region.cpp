#include "cuts.h"
using namespace RooFit;

void Signal_Region(void){
  const bool fix_m_pi0      = false;
  const bool fix_m_d0       = false;
  const bool fix_m_etagg    = false;
  const bool fix_m_etappp   = false;
  const bool fix_m_omega    = false;
  const bool fix_mbc_pi0    = false;
  const bool fix_de_pi0     = false;
  const bool fix_mbc_etagg  = false;
  const bool fix_de_etagg   = false;
  const bool fix_mbc_etappp = false;
  const bool fix_de_etappp  = false;
  const bool fix_mbc_omega  = false;
  const bool fix_de_omega   = false;
  const bool fix_mbc_etapgg = false;
  const bool fix_de_etapgg  = false;
  const bool fix_mbc_etapppp= false;
  const bool fix_de_etapppp = false;
  const bool fix_dm_etapgg  = false;
  const bool fix_dm_etapppp = false;

  const bool skip_m_pi0      = false;
  const bool skip_m_d0       = false;
  const bool skip_m_etagg    = false;
  const bool skip_m_etappp   = false;
  const bool skip_m_omega    = false;
  const bool skip_mbc_pi0    = false;
  const bool skip_de_pi0     = false;
  const bool skip_mbc_etagg  = false;
  const bool skip_de_etagg   = false;
  const bool skip_mbc_etappp = false;
  const bool skip_de_etappp  = false;
  const bool skip_mbc_omega  = false;
  const bool skip_de_omega   = false;
  const bool skip_mbc_etapgg = false;
  const bool skip_de_etapgg  = false;
  const bool skip_mbc_etapppp= false;
  const bool skip_de_etapppp = false;
  const bool skip_dm_etapgg  = false;
  const bool skip_dm_etapppp = false;

  double int_dm_etapppp_sig, int_dm_etapppp_plot;
  double int_dm_etapgg_sig, int_dm_etapgg_plot;
//  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/Tuples/b2dh_data.root");
//  TTree *tree = (TTree*)ifile->Get("TEvent");

  TFile *file_pi0 = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigPi0_s7.root");
  TTree *tree_pi0 = (TTree*)file_pi0->Get("TEvent");

  TFile *file_eta = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigEta_s2.root");
  TTree *tree_eta = (TTree*)file_eta->Get("TEvent");

  TFile *file_etap = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigETAP_s1.root");
  TTree *tree_etap = (TTree*)file_etap->Get("TEvent");

  TFile *file_omega = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigOmega_s5.root");
  TTree *tree_omega = (TTree*)file_omega->Get("TEvent");

  RooArgSet argset;

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good_icpv",1);
  argset.add(good_icpv);
  RooCategory mode("mode","mode");
  mode.defineType("pi0",1);
  mode.defineType("eta",2);
  mode.defineType("omega",3);
  mode.defineType("etap",5);
  argset.add(mode);
  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  h0mode.defineType("ppp",20);
  argset.add(h0mode);
  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  b0f.defineType("wrph",5);
  b0f.defineType("fsr",10);
  argset.add(b0f);

  RooRealVar chi2_mass_d0("chi2_mass_d0","chi2_mass_d0",0.,50.); argset.add(chi2_mass_d0);
  RooRealVar mh0("mh0","m(h^{0})",0.,1.,"GeV");
  mh0.setRange("pi0_fit",Pi0Mass-0.01,Pi0Mass+0.01);
  mh0.setRange("pi0_plot",0.11,0.16);
  mh0.setRange("eta_fit",EtaMass-0.01,EtaMass+0.01);
  mh0.setRange("eta_plot",EtaMass-0.05,EtaMass+0.048);
  mh0.setRange("omega_fit",OmegaMass-0.01,OmegaMass+0.012);
  mh0.setRange("omega_plot",OmegaMass-0.05,OmegaMass+0.05);
  argset.add(mh0);
  RooRealVar mbc("mbc","M_{bc}",5.24,5.291,"GeV");
  mbc.setRange("mbc_fit",5.28-0.01,5.28+0.01);
  mbc.setRange("mbc_plot",5.24,5.291);
  argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.3,0.3,"GeV");
  de.setRange("de_fit",-0.03,0.03);
  de.setRange("de_plot",-0.3,0.3);
  argset.add(de);
  RooRealVar md("md_raw","m(D^{0})",1.8,1.92,"GeV");
  md.setRange("md_fit",DMass-0.01,DMass+0.01);
  md.setRange("md_plot",1.8,1.92);
  argset.add(md);
  RooRealVar mk("mk","m(K_{S}^{0})",0.48,0.52,"GeV");   argset.add(mk);
  RooRealVar mpi0("mpi0","m(#pi^{0})",0.08,0.18,"GeV");// argset.add(mpi0);
  RooRealVar bdt("bdt","BDT",0.,1.); argset.add(bdt);

  RooRealVar dmetap("dmetap","#Deltam(#eta`)",0.395,0.425);
  dmetap.setRange("dmetap_fit",0.410-0.003,0.410+0.003);
  dmetap.setRange("dmetap_plot",0.395,0.425);


  stringstream out;
  // ** Definition of the cuts ** //
  const string pi0modecut    = string("mode == 1 && h0mode == 10");
  const string etaggmodecut  = string("mode == 2 && h0mode == 10");
  const string etapppmodecut = string("mode == 2 && h0mode == 20");
  const string omegamodecut  = string("mode == 3 && h0mode == 20");
  const string etapggmodecut  = string("mode == 5 && h0mode == 10");
  const string etappppmodecut  = string("mode == 5 && h0mode == 20");
  out.str("");
//  out << "abs(md_raw-" << DMass << ")<" << md_cut;
  out << "md_raw>" << md_min << " && md_raw<" << md_max;
  const string mdcut = out.str();
  out.str("");
  out << "abs(mk-" <<  KMass << ")<" << mk_cut;
  const string mkcut = out.str();
  out.str("");
//  out << "de>" << de_min << " && de<" << de_max;
  out << "de>" << de_min_pi0 << " && de<" << de_max_pi0;
  const string depi0cut = out.str();
  out.str("");
//  out << "de>" << de_min << " && de<" << de_max;
  out << "de>" << de_min_etagg << " && de<" << de_max_etagg;
  const string deetaggcut = out.str();
  out.str("");
  out << "de>" << de_min_etappp << " && de<" << de_max_etappp;
  const string deetapppcut = out.str();
  out.str("");
  out << "de>" << de_min_etapgg << " && de<" << de_max_etapgg;
  const string deetapggcut = out.str();
  out.str("");
  out << "de>" << de_min_etapppp << " && de<" << de_max_etapppp;
  const string deetappppcut = out.str();
  out.str("");
  out << "de>" << de_min_omega << " && de<" << de_max_omega;
  const string deomegacut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_pi0 << " && mbc<" << mbc_max_pi0;
  const string mbcpi0cut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_etagg << " && mbc<" << mbc_max_etagg;
  const string mbcetaggcut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_etappp << " && mbc<" << mbc_max_etappp;
  const string mbcetapppcut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_etapgg << " && mbc<" << mbc_max_etapgg;
  const string mbcetapggcut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_etapppp << " && mbc<" << mbc_max_etapppp;
  const string mbcetappppcut = out.str();
  out.str("");
  out << "mbc>" << mbc_min_omega << " && mbc<" << mbc_max_omega;
  const string mbcomegacut = out.str();
  out.str("");
  out << "mpi0>" << mpi0_min << " && mpi0<" << mpi0_max;
  const string mpi0cut = out.str();
  out.str("");
  out << "mh0>" << mpi0_min << " && mh0<" << mpi0_max;
  const string mh0pi0cut = out.str();
  out.str("");
  out << "mh0>" << metagg_min << " && mh0<" << metagg_max;
  const string mh0etaggcut = out.str();
  out.str("");
  out << "mh0>" << metappp_min << " && mh0<" << metappp_max;
  const string mh0etapppcut = out.str();
  out.str("");
  out << "dmetap>" << dmetapgg_min << " && dmetap<" << dmetapgg_max;
  const string dmetapggcut = out.str();
  out.str("");
  out << "dmetap>" << dmetapppp_min << " && dmetap<" << dmetapppp_max;
  const string dmetappppcut = out.str();
  out.str("");
//  out << "abs(mh0-" << OmegaMass << ")<" << momega_cut;
  out << "mh0>" << momega_min << " && mh0<" << momega_max;
  const string mh0omegacut = out.str();
  out.str("");
  out << "bdt>" << bdt_cut_pi0;
  const string bdtpi0cut = out.str();
  out.str("");
  out << "bdt>" << bdt_cut_etagg;
  const string bdtetaggcut = out.str();
  out.str("");
  out << "bdt>" << bdt_cut_etappp;
  const string bdtetapppcut = out.str();
  out.str("");
  out << "bdt>" << bdt_cut_omega;
  const string bdtomegacut = out.str();

  const string ps(" && ");
  // ** End ** //

  RooRealVar mean("mean","mean",1.,-0.1,10.);
  RooRealVar sl("sl","sl",0.1,0.,10.);
  RooRealVar sr("sr","sr",0.1,0.,10.);
  RooBifurGauss pdf_mh0("pdf_mh0","pdf_mh0",mh0,mean,sl,sr);
  RooBifurGauss pdf_md0("pdf_md0","pdf_md0",md,mean,sl,sr);
  RooBifurGauss pdf_mbc("pdf_mbc","pdf_mbc",mbc,mean,sl,sr);
  RooBifurGauss pdf_de("pdf_de","pdf_de",de,mean,sl,sr);
  RooBifurGauss pdf_dmetap("pdf_dmetap","pdf_dmetap",dmetap,mean,sl,sr);

  RooRealVar  s1("s1","s1",0.005,0.,0.1);
  RooGaussian g1("g1","g1",mh0,mean,s1);
  RooRealVar fg1("fg1","fg1",0.9,0.,1.);
  RooGaussian g1d("g1d","g1d",md,mean,s1);
  RooGaussian g1de("g1de","g1de",de,mean,s1);
  RooGaussian g1dm("g1dm","g1dm",dmetap,mean,s1);

  RooRealVar alphal("alphal","alphal",1., 0.,10.); alphal.setConstant(kTRUE);
  RooRealVar nl("nl","nl",5.8,0.,100.);// nl.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",-1.,-10.,0.); alphar.setConstant(kTRUE);
  RooRealVar nr("nr","nr",2.2,0.,100.);// nr.setConstant(kTRUE);
  RooRealVar alpha_mbc("alpha_mbc","alpha_mbc",0.139,0.01,2.);// alpha.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",mh0,mean,sl,alphal,nl);
  RooCBShape CBr("CBr","CBr",mh0,mean,sr,alphar,nr);
  RooCBShape CBld("CBld","CBld",md,mean,sl,alphal,nl);
  RooCBShape CBrd("CBrd","CBrd",md,mean,sr,alphar,nr);

  RooRealVar deCBl("deCBl","deCBl",0.,-0.2,0.1);
  RooRealVar deCBr("deCBr","deCBr",0.,-0.2,0.1);
  RooCBShape CBlde("CBlde","CBlde",de,deCBl,sl,alphal,nl);
  RooCBShape CBrde("CBrde","CBrde",de,deCBr,sr,alphar,nr);
  RooCBShape CBlde0("CBlde0","CBlde0",de,deCBl,s1,alphal,nl);
//  RooCBShape CBrde0("CBrde0","CBrde0",de,mean,sr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",0.46,0.,1.);
  RooRealVar fCBr("fCBr","fCBr",0.31,0.,1.);

  RooAddPdf pdf_m_h0("pdf_m_h0","pdf_m_h0",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  RooAddPdf pdf_m_d0("pdf_m_d0","pdf_m_d0",RooArgList(CBld,CBrd,g1d),RooArgSet(fCBl,fCBr));
  RooNovosibirsk pdf_m_mbc("pdf_m_mbc","pdf_m_mbc",mbc,mean,s1,alpha_mbc);
  RooAddPdf pdf_m_h0_dg("pdf_m_h0_dg","pdf_m_h0_dg",RooArgList(g1,pdf_mh0),RooArgSet(fg1));
  RooAddPdf pdf_dm_etap_dg("pdf_dm_etap_dg","pdf_dm_etap_dg",RooArgList(g1dm,pdf_dmetap),RooArgSet(fg1));
  RooAddPdf pdf_m_de("pdf_m_de","pdf_m_de",RooArgList(CBlde,CBrde,g1de),RooArgSet(fCBl,fCBr));
//  RooAddPdf pdf_m_de1("pdf_m_de1","pdf_m_de1",RooArgList(g1de,pdf_de),RooArgSet(fg1));
  RooAddPdf pdf_m_de0("pdf_m_de0","pdf_m_de0",RooArgList(CBlde0,pdf_de),RooArgSet(fCBl));
  // * * //

  if(!skip_m_omega){
  // m_omega //
  const string m_omega_cuts = omegamodecut + ps + mdcut + ps + mkcut + ps + mbcomegacut + ps + deomegacut + ps + bdtomegacut + ps + mpi0cut;
  RooDataSet ds_m_omega("ds_m_omega","ds_m_omega",tree_omega,RooArgSet(argset,mpi0),m_omega_cuts.c_str());

  mean.setVal(7.81904e-01);
  sl.setVal(7.15783e-03);
  sr.setVal(7.32162e-03);
  if(fix_m_omega){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mh0.fitTo(ds_m_omega,Timer(true),Range("omega_fit"));
  const double m_omega_min = mean.getVal() - 3.*sl.getVal();
  const double m_omega_max = mean.getVal() + 3.*sr.getVal();

  fCBl.setVal(4.68894e-01);
  fCBr.setVal(3.23399e-01);
  mean.setVal(7.81940e-01);
  nl.setVal(5.27939e+00);
  nr.setVal(2.48841e+00);
  s1.setVal(4.56829e-03);
  sl.setVal(9.31832e-03);
  sr.setVal(9.77746e-03);
  if(fix_m_omega){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_h0.fitTo(ds_m_omega,Timer(true),Range("omega_plot"),Save());
  mh0.setRange("omega_sig",m_omega_min,m_omega_max);
  const double int_m_omega_sig = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("omega_sig"))->getVal();
  const double int_m_omega_plot = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("omega_plot"))->getVal();

  RooPlot* m_omegaFrame = mh0.frame(Range("omega_plot"),Title("Signal range of #omega mass"));
  ds_m_omega.plotOn(m_omegaFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_h0.plotOn(m_omegaFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* m_omega_cm = new TCanvas("m_omega_cm","m_omega_cm",600,400);
  m_omega_cm->cd();
  ds_m_omega.statOn(m_omegaFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_omegaFrame->GetXaxis()->SetTitleSize(0.05);
  m_omegaFrame->GetXaxis()->SetTitleOffset(0.85);
  m_omegaFrame->GetXaxis()->SetLabelSize(0.05);
  m_omegaFrame->GetYaxis()->SetTitleOffset(1.6);
  m_omegaFrame->Draw();
  TLine *momegalineLEFT = new TLine(m_omega_min,0,m_omega_min,3200);
  momegalineLEFT->SetLineColor(kRed);
  momegalineLEFT->SetLineWidth(2);
  momegalineLEFT->SetLineStyle(1);
  momegalineLEFT->Draw();
  TLine *momegalineRIGHT = new TLine(m_omega_max,0,m_omega_max,3200);
  momegalineRIGHT->SetLineColor(kRed);
  momegalineRIGHT->SetLineWidth(2);
  momegalineRIGHT->SetLineStyle(1);
  momegalineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_m_omega_sig/int_m_omega_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_omega_cm->Update();
  m_omega_cm->Print("pics/m_omega_cut.eps");
  m_omega_cm->Print("pics/m_omega_cut.root");
  // ** end of m(omega) ** //
  }// skip

  if(!skip_m_etagg){
  // m(eta->gamma gamma) //
  const string m_etagg_cuts = etaggmodecut + ps + mdcut + ps + mkcut + ps + mbcetaggcut + ps + deetaggcut + ps + bdtetaggcut;
  RooDataSet ds_m_etagg("ds_m_etagg","ds_m_etagg",tree_eta,argset,m_etagg_cuts.c_str());
  ds_m_etagg.Print();

  mean.setVal(5.46546e-01);
  sl.setVal(9.64747e-03);
  sr.setVal(9.03598e-03);
  if(fix_m_etagg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mh0.fitTo(ds_m_etagg,Timer(true),Range("eta_fit"));
  const double m_etagg_min = mean.getVal() - 3.*sl.getVal();
  const double m_etagg_max = mean.getVal() + 3.*sr.getVal();

  fCBl.setVal(5.63992e-01);
  fCBr.setVal(2.03950e-02);
  mean.setVal(5.45632e-01);
  nl.setVal(9.99990e+01);
  nr.setVal(8.96532e-01);
  s1.setVal(1.19950e-02);
  sl.setVal(9.30707e-03);
  sr.setVal(1.55787e-03);
  if(fix_m_etagg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_h0.fitTo(ds_m_etagg,Timer(true),Range("eta_plot"));
  mh0.setRange("etagg_sig",m_etagg_min,m_etagg_max);
  const double int_m_etagg_sig = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("etagg_sig"))->getVal();
  const double int_m_etagg_plot = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("eta_plot"))->getVal();

  RooPlot* m_etaggFrame = mh0.frame(Range("eta_plot"),Title("Signal range of #eta#rightarrow#gamma#gamma mass"));
  ds_m_etagg.plotOn(m_etaggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_h0.plotOn(m_etaggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* m_etagg_cm = new TCanvas("m_etagg_cm","m_etagg_cm",600,400);
  m_etagg_cm->cd();
  ds_m_etagg.statOn(m_etaggFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_etaggFrame->GetXaxis()->SetTitleSize(0.05);
  m_etaggFrame->GetXaxis()->SetTitleOffset(0.85);
  m_etaggFrame->GetXaxis()->SetLabelSize(0.05);
  m_etaggFrame->GetYaxis()->SetTitleOffset(1.6);
  m_etaggFrame->Draw();
  TLine *metagglineLEFT = new TLine(m_etagg_min,0,m_etagg_min,3200);
  metagglineLEFT->SetLineColor(kRed);
  metagglineLEFT->SetLineWidth(2);
  metagglineLEFT->SetLineStyle(1);
  metagglineLEFT->Draw();
  TLine *metagglineRIGHT = new TLine(m_etagg_max,0,m_etagg_max,3200);
  metagglineRIGHT->SetLineColor(kRed);
  metagglineRIGHT->SetLineWidth(2);
  metagglineRIGHT->SetLineStyle(1);
  metagglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_m_etagg_sig/int_m_etagg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_etagg_cm->Update();
  m_etagg_cm->Print("pics/m_etagg_cut.eps");
  m_etagg_cm->Print("pics/m_etagg_cut.root");
  // ** end of m(eta->gamma gamma) ** //
  }//skip

  if(!skip_dm_etapgg){
  // dm(eta') (eta->gamma gamma) //
  const string dm_etapgg_cuts = etapggmodecut + ps + mdcut + ps + mkcut + ps + mbcetapggcut + ps + deetapggcut + ps + bdtetaggcut + ps + mh0etaggcut;
  RooDataSet ds_dm_etapgg("ds_dm_etapgg","ds_dm_etapgg",tree_etap,RooArgSet(argset,dmetap),dm_etapgg_cuts.c_str());
  ds_dm_etapgg.Print();

  mean.setVal(4.10197e-01);
  sl.setVal(2.83874e-03);
  sr.setVal(2.50381e-03);
  if(fix_dm_etapgg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_dmetap.fitTo(ds_dm_etapgg,Timer(true),Range("dmetap_fit"));
  const double dm_etapgg_min = mean.getVal() - 3.*sl.getVal();
  const double dm_etapgg_max = mean.getVal() + 3.*sr.getVal();

  fg1.setVal(3.55901e-01);
  mean.setVal(4.10007e-01);
  s1.setVal(4.66009e-03);
  sl.setVal(2.24283e-03);
  sr.setVal(2.60270e-03);
  if(fix_dm_etapgg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fg1.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fg1.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_dm_etap_dg.fitTo(ds_dm_etapgg,Timer(true),Range("dmetap_plot"));
  dmetap.setRange("dmetapgg_sig",dm_etapgg_min,dm_etapgg_max);
  int_dm_etapgg_sig = pdf_dm_etap_dg.createIntegral(RooArgSet(dmetap),NormSet(RooArgSet(dmetap)),Range("dmetapgg_sig"))->getVal();
  int_dm_etapgg_plot = pdf_dm_etap_dg.createIntegral(RooArgSet(dmetap),NormSet(RooArgSet(dmetap)),Range("dmetap_plot"))->getVal();

  RooPlot* dm_etapggFrame = dmetap.frame(Range("detap_plot"),Title("Signal range of #Deltam(#eta`) for #eta#rightarrow#gamma#gamma"));
  ds_dm_etapgg.plotOn(dm_etapggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));

  pdf_dm_etap_dg.plotOn(dm_etapggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* dm_etapgg_cm = new TCanvas("dm_etapgg_cm","dm_etapgg_cm",600,400);
  dm_etapgg_cm->cd();
  ds_dm_etapgg.statOn(dm_etapggFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  dm_etapggFrame->GetXaxis()->SetTitleSize(0.05);
  dm_etapggFrame->GetXaxis()->SetTitleOffset(0.85);
  dm_etapggFrame->GetXaxis()->SetLabelSize(0.05);
  dm_etapggFrame->GetYaxis()->SetTitleOffset(1.6);
  dm_etapggFrame->Draw();
  TLine *dmetapgglineLEFT = new TLine(dm_etapgg_min,0,dm_etapgg_min,150);
  dmetapgglineLEFT->SetLineColor(kRed);
  dmetapgglineLEFT->SetLineWidth(2);
  dmetapgglineLEFT->SetLineStyle(1);
  dmetapgglineLEFT->Draw();
  TLine *dmetapgglineRIGHT = new TLine(dm_etapgg_max,0,dm_etapgg_max,150);
  dmetapgglineRIGHT->SetLineColor(kRed);
  dmetapgglineRIGHT->SetLineWidth(2);
  dmetapgglineRIGHT->SetLineStyle(1);
  dmetapgglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_dm_etapgg_sig/int_dm_etapgg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  dm_etapgg_cm->Update();
  dm_etapgg_cm->Print("pics/dm_etapgg_cut.eps");
  dm_etapgg_cm->Print("pics/dm_etapgg_cut.root");
  // ** end of dm(eta') (eta->gamma gamma) ** //
  }//skip

  if(!skip_m_etappp){
  // m(eta->p+p-p0) //
  const string m_etappp_cuts = etapppmodecut + ps + mdcut + ps + mkcut + ps + mbcetapppcut + ps + deetapppcut + ps + bdtetapppcut + ps + mpi0cut;
  RooDataSet ds_m_etappp("ds_m_etappp","ds_m_etappp",tree_eta,RooArgSet(argset,mpi0),m_etappp_cuts.c_str());

  mean.setVal(5.47517e-01);
  sl.setVal(3.29716e-03);
  sr.setVal(3.30856e-03);
  if(fix_m_etappp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mh0.fitTo(ds_m_etappp,Timer(true),Range("eta_fit"));
  const double m_etappp_min = mean.getVal() - 3.*sl.getVal();
  const double m_etappp_max = mean.getVal() + 3.*sr.getVal();

  fCBl.setVal(4.01138e-01);
  fCBr.setVal(2.95008e-01);
  mean.setVal(5.47548e-01);
  nl.setVal(4.10093e+00);
  nr.setVal(2.13769e+00);
  s1.setVal(3.61594e-03);
  sl.setVal(2.37987e-03);
  sr.setVal(2.31083e-03);
  if(fix_m_etappp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_h0.fitTo(ds_m_etappp,Timer(true),Range("eta_plot"));
  mh0.setRange("etappp_sig",m_etappp_min,m_etappp_max);
  const double int_m_etappp_sig = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("etappp_sig"))->getVal();
  const double int_m_etappp_plot = pdf_m_h0.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("eta_plot"))->getVal();

  RooPlot* m_etapppFrame = mh0.frame(Range("eta_plot"),Title("Signal range of #eta#rightarrow#pi^{+}#pi^{-}#pi^{0} mass"));
  ds_m_etappp.plotOn(m_etapppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_h0.plotOn(m_etapppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* m_etappp_cm = new TCanvas("m_etappp_cm","m_etappp_cm",600,400);
  m_etappp_cm->cd();
  ds_m_etappp.statOn(m_etapppFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_etapppFrame->GetXaxis()->SetTitleSize(0.05);
  m_etapppFrame->GetXaxis()->SetTitleOffset(0.85);
  m_etapppFrame->GetXaxis()->SetLabelSize(0.05);
  m_etapppFrame->GetYaxis()->SetTitleOffset(1.6);
  m_etapppFrame->Draw();
  TLine *metappplineLEFT = new TLine(m_etappp_min,0,m_etappp_min,3200);
  metappplineLEFT->SetLineColor(kRed);
  metappplineLEFT->SetLineWidth(2);
  metappplineLEFT->SetLineStyle(1);
  metappplineLEFT->Draw();
  TLine *metappplineRIGHT = new TLine(m_etappp_max,0,m_etappp_max,3200);
  metappplineRIGHT->SetLineColor(kRed);
  metappplineRIGHT->SetLineWidth(2);
  metappplineRIGHT->SetLineStyle(1);
  metappplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_m_etappp_sig/int_m_etappp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_etappp_cm->Update();
  m_etappp_cm->Print("pics/m_etappp_cut.eps");
  m_etappp_cm->Print("pics/m_etappp_cut.root");
  // ** end of m(eta->p+p-p0) ** //
  }// skip

  if(!skip_dm_etapppp){
  // dm(eta') (eta->pi+ pi- pi0) //
  RooArgSet newargset(argset);
  newargset.add(dmetap);
  newargset.add(mpi0);
  const string dm_etapppp_cuts = etappppmodecut + ps + mdcut + ps + mkcut + ps + mbcetappppcut + ps + deetappppcut + ps + bdtetapppcut + ps + mh0etapppcut + ps + mpi0cut;
  RooDataSet ds_dm_etapppp("ds_dm_etapppp","ds_dm_etapppp",tree_etap,newargset,dm_etapppp_cuts.c_str());
  ds_dm_etapppp.Print();

  mean.setVal(4.10385e-01);
  sl.setVal(2.66981e-03);
  sr.setVal(2.22528e-03);
  if(fix_dm_etapppp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_dmetap.fitTo(ds_dm_etapppp,Timer(true),Range("dmetap_fit"));
  const double dm_etapppp_min = mean.getVal() - 3.*sl.getVal();
  const double dm_etapppp_max = mean.getVal() + 3.*sr.getVal();

  fg1.setVal(2.85758e-01);
  mean.setVal(4.10251e-01);
  s1.setVal(7.99697e-03);
  sl.setVal(2.17622e-03);
  sr.setVal(2.02066e-03);
  if(fix_dm_etapppp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fg1.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fg1.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_dm_etap_dg.fitTo(ds_dm_etapppp,Timer(true),Range("dmetap_plot"));
  dmetap.setRange("dmetapppp_sig",dm_etapppp_min,dm_etapppp_max);
  int_dm_etapppp_sig = pdf_dm_etap_dg.createIntegral(RooArgSet(dmetap),NormSet(RooArgSet(dmetap)),Range("dmetapppp_sig"))->getVal();
  int_dm_etapppp_plot = pdf_dm_etap_dg.createIntegral(RooArgSet(dmetap),NormSet(RooArgSet(dmetap)),Range("dmetap_plot"))->getVal();

  RooPlot* dm_etappppFrame = dmetap.frame(Range("detap_plot"),Title("Signal range of #Deltam(#eta`) for #eta#rightarrow#pi^{+}#pi^{-}#pi^{0}"));
  ds_dm_etapppp.plotOn(dm_etappppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));

  pdf_dm_etap_dg.plotOn(dm_etappppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* dm_etapppp_cm = new TCanvas("dm_etapppp_cm","dm_etapppp_cm",600,400);
  dm_etapppp_cm->cd();
  ds_dm_etapppp.statOn(dm_etappppFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  dm_etappppFrame->GetXaxis()->SetTitleSize(0.05);
  dm_etappppFrame->GetXaxis()->SetTitleOffset(0.85);
  dm_etappppFrame->GetXaxis()->SetLabelSize(0.05);
  dm_etappppFrame->GetYaxis()->SetTitleOffset(1.6);
  dm_etappppFrame->Draw();
  TLine *dmetapppplineLEFT = new TLine(dm_etapppp_min,0,dm_etapppp_min,60);
  dmetapppplineLEFT->SetLineColor(kRed);
  dmetapppplineLEFT->SetLineWidth(2);
  dmetapppplineLEFT->SetLineStyle(1);
  dmetapppplineLEFT->Draw();
  TLine *dmetapppplineRIGHT = new TLine(dm_etapppp_max,0,dm_etapppp_max,60);
  dmetapppplineRIGHT->SetLineColor(kRed);
  dmetapppplineRIGHT->SetLineWidth(2);
  dmetapppplineRIGHT->SetLineStyle(1);
  dmetapppplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_dm_etapppp_sig/int_dm_etapppp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  dm_etapppp_cm->Update();
  dm_etapppp_cm->Print("pics/dm_etapppp_cut.eps");
  dm_etapppp_cm->Print("pics/dm_etapppp_cut.root");
  // ** end of dm(eta') (eta->pi+ pi- pi0) ** //
  }//skip

  if(!skip_m_pi0){
  // m(pi0) //
  const string m_pi0_cuts = pi0modecut + ps + mdcut + ps + mkcut + ps + mbcpi0cut + ps + depi0cut + ps + bdtpi0cut;
  RooDataSet ds_m_pi0("ds_m_pi0","ds_m_pi0",tree_pi0,RooArgSet(argset,mpi0),m_pi0_cuts.c_str());

  mean.setVal(1.34330e-01);
  sl.setVal(6.19746e-03);
  sr.setVal(6.46114e-03);
  if(fix_m_pi0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mh0.fitTo(ds_m_pi0,Timer(true),Range("pi0_fit"));
  const double m_pi0_min = mean.getVal() - 3.*sl.getVal();
  const double m_pi0_max = mean.getVal() + 3.*sr.getVal();

  fg1.setVal(5.84785e-01);
  mean.setVal(1.34535e-01);
  s1.setVal(5.45215e-03);
  sl.setVal(9.95957e-03);
  sr.setVal(9.12687e-03);
  if(fix_m_pi0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fg1.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fg1.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_h0_dg.fitTo(ds_m_pi0,Timer(true),Range("pi0_plot"));
  mh0.setRange("pi0_sig",m_pi0_min,m_pi0_max);
  const double int_m_pi0_sig  = pdf_m_h0_dg.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("pi0_sig"))->getVal();
  const double int_m_pi0_plot = pdf_m_h0_dg.createIntegral(RooArgSet(mh0),NormSet(RooArgSet(mh0)),Range("pi0_plot"))->getVal();

  RooPlot* m_pi0Frame = mh0.frame(Range("pi0_plot"),Title("Signal range of #pi^{0} mass"));
  ds_m_pi0.plotOn(m_pi0Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_h0_dg.plotOn(m_pi0Frame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* m_pi0_cm = new TCanvas("m_pi0_cm","m_pi0_cm",600,400);
  m_pi0_cm->cd();
  ds_m_pi0.statOn(m_pi0Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_pi0Frame->GetXaxis()->SetTitleSize(0.05);
  m_pi0Frame->GetXaxis()->SetTitleOffset(0.85);
  m_pi0Frame->GetXaxis()->SetLabelSize(0.05);
  m_pi0Frame->GetYaxis()->SetTitleOffset(1.6);
  m_pi0Frame->Draw();
  TLine *mpi0lineLEFT = new TLine(m_pi0_min,0,m_pi0_min,3200);
  mpi0lineLEFT->SetLineColor(kRed);
  mpi0lineLEFT->SetLineWidth(2);
  mpi0lineLEFT->SetLineStyle(1);
  mpi0lineLEFT->Draw();
  TLine *mpi0lineRIGHT = new TLine(m_pi0_max,0,m_pi0_max,3200);
  mpi0lineRIGHT->SetLineColor(kRed);
  mpi0lineRIGHT->SetLineWidth(2);
  mpi0lineRIGHT->SetLineStyle(1);
  mpi0lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_m_pi0_sig/int_m_pi0_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_pi0_cm->Update();
  m_pi0_cm->Print("pics/m_pi0_cut.eps");
  m_pi0_cm->Print("pics/m_pi0_cut.root");
  // ** end of m(pi^0) ** //
  }//skip

  if(!skip_m_d0){
  // m(D0) //
  const string m_d0_cuts = pi0modecut + ps + mkcut + ps + mbcpi0cut + ps + depi0cut + ps + bdtpi0cut + ps + mh0pi0cut;
  RooDataSet ds_m_d0("ds_m_d0","ds_m_d0",tree_pi0,RooArgSet(argset,mpi0),m_d0_cuts.c_str());

  mean.setVal(1.86452e+00);
  sl.setVal(4.32366e-03);
  sr.setVal(4.60337e-03);
  if(fix_m_d0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_md0.fitTo(ds_m_d0,Timer(true),Range("md_fit"));
  const double m_d0_min = mean.getVal() - 3.*sl.getVal();
  const double m_d0_max = mean.getVal() + 3.*sr.getVal();

  alphal.setVal(1.81613e+00);
  alphar.setVal(-2.00119e+00);
  fCBl.setVal(6.21124e-02);
  fCBr.setVal(5.79150e-01);
  mean.setVal(1.86466e+00);
  nl.setVal(6.78465e+00);
  nr.setVal(2.15390e+00);
  s1.setVal(6.07718e-03);
  sl.setVal(1.35035e-02);
  sr.setVal(3.59869e-03);
  if(fix_m_d0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    alphal.setConstant(kTRUE);
    alphar.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    alphal.setConstant(kFALSE);
    alphar.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_d0.fitTo(ds_m_d0,Timer(true),Range("md_plot"));
  md.setRange("md_sig",m_d0_min,m_d0_max);
  const double int_m_d0_sig = pdf_m_d0.createIntegral(RooArgSet(md),NormSet(RooArgSet(md)),Range("md_sig"))->getVal();
  const double int_m_d0_plot = pdf_m_d0.createIntegral(RooArgSet(md),NormSet(RooArgSet(md)),Range("md_plot"))->getVal();

  RooPlot* m_d0Frame = md.frame(Range("md_plot"),Title("Signal range of D^{0} mass"));
  ds_m_d0.plotOn(m_d0Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_d0.plotOn(m_d0Frame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* m_d0_cm = new TCanvas("m_d0_cm","m_d0_cm",600,400);
  m_d0_cm->cd();
  ds_m_d0.statOn(m_d0Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_d0Frame->GetXaxis()->SetTitleSize(0.05);
  m_d0Frame->GetXaxis()->SetTitleOffset(0.85);
  m_d0Frame->GetXaxis()->SetLabelSize(0.05);
  m_d0Frame->GetYaxis()->SetTitleOffset(1.6);
  m_d0Frame->Draw();
  TLine *md0lineLEFT = new TLine(m_d0_min,0,m_d0_min,3200);
  md0lineLEFT->SetLineColor(kRed);
  md0lineLEFT->SetLineWidth(2);
  md0lineLEFT->SetLineStyle(1);
  md0lineLEFT->Draw();
  TLine *md0lineRIGHT = new TLine(m_d0_max,0,m_d0_max,3200);
  md0lineRIGHT->SetLineColor(kRed);
  md0lineRIGHT->SetLineWidth(2);
  md0lineRIGHT->SetLineStyle(1);
  md0lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_m_d0_sig/int_m_d0_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_d0_cm->Update();
  m_d0_cm->Print("pics/m_d0_cut.eps");
  m_d0_cm->Print("pics/m_d0_cut.root");

  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of m(D0) ** //
  }//skip

  if(!skip_mbc_pi0){
  // Mbc for D0 pi0 //
  const string mbc_pi0_cuts = pi0modecut + ps + mkcut + ps + depi0cut + ps + bdtpi0cut + ps + mh0pi0cut + ps + mdcut;
  RooDataSet ds_mbc_pi0("ds_mbc_pi0","ds_mbc_pi0",tree_pi0,RooArgSet(argset,mpi0),mbc_pi0_cuts.c_str());

  mean.setVal(5.28014e+00);
  sl.setVal(3.68288e-03);
  sr.setVal(2.58289e-03);
  if(fix_mbc_pi0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_pi0,Timer(true),Range("mbc_fit"));
  const double mbc_pi0_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_pi0_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(1.40466e-01);
  mean.setVal(5.27987e+00);
  s1.setVal(3.15806e-03);
  if(fix_mbc_pi0){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_pi0,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_pi0",mbc_pi0_min,mbc_pi0_max);
  const double int_mbc_pi0_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_pi0"))->getVal();
  const double int_mbc_pi0_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_pi0Frame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#pi^{0}"));
  ds_mbc_pi0.plotOn(mbc_pi0Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_pi0Frame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_pi0_cm = new TCanvas("mbc_pi0_cm","mbc_pi0_cm",600,400);
  mbc_pi0_cm->cd();
  ds_mbc_pi0.statOn(mbc_pi0Frame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_pi0Frame->GetXaxis()->SetTitleSize(0.05);
  mbc_pi0Frame->GetXaxis()->SetTitleOffset(0.85);
  mbc_pi0Frame->GetXaxis()->SetLabelSize(0.05);
  mbc_pi0Frame->GetYaxis()->SetTitleOffset(1.6);
  mbc_pi0Frame->Draw();
  TLine *mbc_pi0lineLEFT = new TLine(mbc_pi0_min,0,mbc_pi0_min,3200);
  mbc_pi0lineLEFT->SetLineColor(kRed);
  mbc_pi0lineLEFT->SetLineWidth(2);
  mbc_pi0lineLEFT->SetLineStyle(1);
  mbc_pi0lineLEFT->Draw();
  TLine *mbc_pi0lineRIGHT = new TLine(mbc_pi0_max,0,mbc_pi0_max,3200);
  mbc_pi0lineRIGHT->SetLineColor(kRed);
  mbc_pi0lineRIGHT->SetLineWidth(2);
  mbc_pi0lineRIGHT->SetLineStyle(1);
  mbc_pi0lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_pi0_sig/int_mbc_pi0_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_pi0_cm->Update();
  mbc_pi0_cm->Print("pics/mbc_pi0_cut.eps");
  mbc_pi0_cm->Print("pics/mbc_pi0_cut.root");
  // ** end of mbc for D pi0 ** //
  }//skip

  if(!skip_de_pi0){
  // dE for D0 pi0 //
  const string de_pi0_cuts = pi0modecut + ps + mdcut + ps + mkcut + ps + mbcpi0cut + ps + bdtpi0cut + ps + mpi0cut;
  RooDataSet ds_de_pi0("ds_de_pi0","ds_de_pi0",tree_pi0,RooArgSet(argset,mpi0),de_pi0_cuts.c_str());

  mean.setVal(2.71374e-03);
  sl.setVal(4.36853e-02);
  sr.setVal(2.98068e-02);
  if(fix_de_pi0){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_pi0,Timer(true),Range("de_fit"));
  const double de_pi0_min = mean.getVal() - 3.*sl.getVal();
  const double de_pi0_max = mean.getVal() + 3.*sr.getVal();

  alphal.setVal(5.67008e-01);
  alphar.setVal(-1.01331e+00);
  deCBl.setVal(-1.47561e-02);
  deCBr.setVal(-2.37767e-02);
  fCBl.setVal(7.07311e-01);
  fCBr.setVal(1.08329e-01);
  mean.setVal(1.85549e-02);
  nl.setVal(7.57483e+00);
  nr.setVal(8.57788e+00);
  s1.setVal(2.03879e-02);
  sl.setVal(3.20359e-02);
  sr.setVal(3.79709e-02);

  if(fix_de_pi0){
    deCBl.setConstant(kTRUE);
    deCBr.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    alphal.setConstant(kTRUE);
    alphar.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    deCBr.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    alphal.setConstant(kFALSE);
    alphar.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de.fitTo(ds_de_pi0,Timer(true),Range("de_plot"));
  de.setRange("de_pi0_sig",de_pi0_min,de_pi0_max);
  const double int_de_pi0_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_pi0_sig"))->getVal();
  const double int_de_pi0_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_pi0Frame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#pi^{0}"));
  ds_de_pi0.plotOn(de_pi0Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(de_pi0Frame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_pi0_cm = new TCanvas("de_pi0_cm","de_pi0_cm",600,400);
  de_pi0_cm->cd();
  ds_de_pi0.statOn(de_pi0Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_pi0Frame->GetXaxis()->SetTitleSize(0.05);
  de_pi0Frame->GetXaxis()->SetTitleOffset(0.85);
  de_pi0Frame->GetXaxis()->SetLabelSize(0.05);
  de_pi0Frame->GetYaxis()->SetTitleOffset(1.6);
  de_pi0Frame->Draw();
  TLine *de_pi0lineLEFT = new TLine(de_pi0_min,0,de_pi0_min,3200);
  de_pi0lineLEFT->SetLineColor(kRed);
  de_pi0lineLEFT->SetLineWidth(2);
  de_pi0lineLEFT->SetLineStyle(1);
  de_pi0lineLEFT->Draw();
  TLine *de_pi0lineRIGHT = new TLine(de_pi0_max,0,de_pi0_max,3200);
  de_pi0lineRIGHT->SetLineColor(kRed);
  de_pi0lineRIGHT->SetLineWidth(2);
  de_pi0lineRIGHT->SetLineStyle(1);
  de_pi0lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_pi0_sig/int_de_pi0_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_pi0_cm->Update();
  de_pi0_cm->Print("pics/de_pi0_cut.eps");
  de_pi0_cm->Print("pics/de_pi0_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 pi0 ** //
  }//skip

  if(!skip_mbc_etagg){
  // Mbc for D0 eta->gg //
  const string mbc_etagg_cuts = etaggmodecut + ps + mkcut + ps + deetaggcut + ps + bdtetaggcut + ps + mh0etaggcut + ps + mdcut;
  RooDataSet ds_mbc_etagg("ds_mbc_etagg","ds_mbc_etagg",tree_eta,argset,mbc_etagg_cuts.c_str());

  mean.setVal(5.28005e+00);
  sl.setVal(3.43095e-03);
  sr.setVal(2.52900e-03);
  if(fix_mbc_etagg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_etagg,Timer(true),Range("mbc_fit"));
  const double mbc_etagg_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_etagg_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(1.08683e-01);
  mean.setVal(5.27979e+00);
  s1.setVal(2.98633e-03);
  if(fix_mbc_etagg){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_etagg,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_etagg",mbc_etagg_min,mbc_etagg_max);
  const double int_mbc_etagg_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_etagg"))->getVal();
  const double int_mbc_etagg_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_etaggFrame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#eta(#rightarrow#gamma#gamma)"));
  ds_mbc_etagg.plotOn(mbc_etaggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_etaggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_etagg_cm = new TCanvas("mbc_etagg_cm","mbc_etagg_cm",600,400);
  mbc_etagg_cm->cd();
  ds_mbc_etagg.statOn(mbc_etaggFrame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_etaggFrame->GetXaxis()->SetTitleSize(0.05);
  mbc_etaggFrame->GetXaxis()->SetTitleOffset(0.85);
  mbc_etaggFrame->GetXaxis()->SetLabelSize(0.05);
  mbc_etaggFrame->GetYaxis()->SetTitleOffset(1.6);
  mbc_etaggFrame->Draw();
  TLine *mbc_etagglineLEFT = new TLine(mbc_etagg_min,0,mbc_etagg_min,3200);
  mbc_etagglineLEFT->SetLineColor(kRed);
  mbc_etagglineLEFT->SetLineWidth(2);
  mbc_etagglineLEFT->SetLineStyle(1);
  mbc_etagglineLEFT->Draw();
  TLine *mbc_etagglineRIGHT = new TLine(mbc_etagg_max,0,mbc_etagg_max,3200);
  mbc_etagglineRIGHT->SetLineColor(kRed);
  mbc_etagglineRIGHT->SetLineWidth(2);
  mbc_etagglineRIGHT->SetLineStyle(1);
  mbc_etagglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_etagg_sig/int_mbc_etagg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_etagg_cm->Update();
  mbc_etagg_cm->Print("pics/mbc_etagg_cut.eps");
  mbc_etagg_cm->Print("pics/mbc_etagg_cut.root");
  // ** end of mbc for D eta->gg ** //
  }//skip

  if(!skip_de_etagg){
  // dE for D0 eta->gg //
  const string de_etagg_cuts = etaggmodecut + ps + mdcut + ps + mkcut + ps + mbcetaggcut + ps + bdtetaggcut + ps + mh0etaggcut;
  RooDataSet ds_de_etagg("ds_de_etagg","ds_de_etagg",tree_eta,argset,de_etagg_cuts.c_str());

  mean.setVal(1.30334e-03);
  sl.setVal(3.35854e-02);
  sr.setVal(2.48597e-02);
  if(fix_de_etagg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_etagg,Timer(true),Range("de_fit"));
  const double de_etagg_min = mean.getVal() - 3.*sl.getVal();
  const double de_etagg_max = mean.getVal() + 3.*sr.getVal();

  alphal.setVal(9.02616e-01);
  alphar.setVal(-1.61665e+00);
  deCBl.setVal(-2.04410e-02);
  deCBr.setVal(-2.05836e-02);
  fCBl.setVal(4.92314e-01);
  fCBr.setVal(1.34440e-01);
  mean.setVal(7.62352e-03);
  nl.setVal(1.00000e+02);
  nr.setVal(1.21212e+01);
  s1.setVal(2.05598e-02);
  sl.setVal(3.08165e-02);
  sr.setVal(4.10312e-02);
  if(fix_de_etagg){
    deCBl.setConstant(kTRUE);
    deCBr.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    alphal.setConstant(kTRUE);
    alphar.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    deCBr.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    alphal.setConstant(kFALSE);
    alphar.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de.fitTo(ds_de_etagg,Timer(true),Range("de_plot"));
  de.setRange("de_etagg_sig",de_etagg_min,de_etagg_max);
  const double int_de_etagg_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_etagg_sig"))->getVal();
  const double int_de_etagg_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_etaggFrame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#eta(#rightarrow#gamma#gamma)"));
  ds_de_etagg.plotOn(de_etaggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(de_etaggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_etagg_cm = new TCanvas("de_etagg_cm","de_etagg_cm",600,400);
  de_etagg_cm->cd();
  ds_de_etagg.statOn(de_etaggFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_etaggFrame->GetXaxis()->SetTitleSize(0.05);
  de_etaggFrame->GetXaxis()->SetTitleOffset(0.85);
  de_etaggFrame->GetXaxis()->SetLabelSize(0.05);
  de_etaggFrame->GetYaxis()->SetTitleOffset(1.6);
  de_etaggFrame->Draw();
  TLine *de_etagglineLEFT = new TLine(de_etagg_min,0,de_etagg_min,3200);
  de_etagglineLEFT->SetLineColor(kRed);
  de_etagglineLEFT->SetLineWidth(2);
  de_etagglineLEFT->SetLineStyle(1);
  de_etagglineLEFT->Draw();
  TLine *de_etagglineRIGHT = new TLine(de_etagg_max,0,de_etagg_max,3200);
  de_etagglineRIGHT->SetLineColor(kRed);
  de_etagglineRIGHT->SetLineWidth(2);
  de_etagglineRIGHT->SetLineStyle(1);
  de_etagglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_etagg_sig/int_de_etagg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_etagg_cm->Update();
  de_etagg_cm->Print("pics/de_etagg_cut.eps");
  de_etagg_cm->Print("pics/de_etagg_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 eta->gg ** //
  }//skip

  if(!skip_mbc_etappp){
  // Mbc for D0 eta->ppp //
  const string mbc_etappp_cuts = etapppmodecut + ps + mkcut + ps + deetapppcut + ps + bdtetapppcut + ps + mh0etapppcut + ps + mdcut + ps + mpi0cut;
  RooDataSet ds_mbc_etappp("ds_mbc_etappp","ds_mbc_etappp",tree_eta,RooArgSet(argset,mpi0),mbc_etappp_cuts.c_str());

  mean.setVal(5.27985e+00);
  sl.setVal(3.12139e-03);
  sr.setVal(2.49806e-03);
  if(fix_mbc_etappp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_etappp,Timer(true),Range("mbc_fit"));
  const double mbc_etappp_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_etappp_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(9.11142e-02);
  mean.setVal(5.27972e+00);
  s1.setVal(2.82721e-03);
  if(fix_mbc_etappp){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_etappp,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_etappp",mbc_etappp_min,mbc_etappp_max);
  const double int_mbc_etappp_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_etappp"))->getVal();
  const double int_mbc_etappp_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_etapppFrame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  ds_mbc_etappp.plotOn(mbc_etapppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_etapppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_etappp_cm = new TCanvas("mbc_etappp_cm","mbc_etappp_cm",600,400);
  mbc_etappp_cm->cd();
  ds_mbc_etappp.statOn(mbc_etapppFrame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_etapppFrame->GetXaxis()->SetTitleSize(0.05);
  mbc_etapppFrame->GetXaxis()->SetTitleOffset(0.85);
  mbc_etapppFrame->GetXaxis()->SetLabelSize(0.05);
  mbc_etapppFrame->GetYaxis()->SetTitleOffset(1.6);
  mbc_etapppFrame->Draw();
  TLine *mbc_etappplineLEFT = new TLine(mbc_etappp_min,0,mbc_etappp_min,3200);
  mbc_etappplineLEFT->SetLineColor(kRed);
  mbc_etappplineLEFT->SetLineWidth(2);
  mbc_etappplineLEFT->SetLineStyle(1);
  mbc_etappplineLEFT->Draw();
  TLine *mbc_etappplineRIGHT = new TLine(mbc_etappp_max,0,mbc_etappp_max,3200);
  mbc_etappplineRIGHT->SetLineColor(kRed);
  mbc_etappplineRIGHT->SetLineWidth(2);
  mbc_etappplineRIGHT->SetLineStyle(1);
  mbc_etappplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_etappp_sig/int_mbc_etappp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_etappp_cm->Update();
  mbc_etappp_cm->Print("pics/mbc_etappp_cut.eps");
  mbc_etappp_cm->Print("pics/mbc_etappp_cut.root");
  // ** end of mbc for D eta->ppp ** //
  }//skip

  if(!skip_de_etappp){
  // dE for D0 eta->ppp //
  const string de_etappp_cuts = etapppmodecut + ps + mdcut + ps + mkcut + ps + mbcetapppcut + ps + bdtetapppcut + ps + mh0etapppcut + ps + mpi0cut;
  RooDataSet ds_de_etappp("ds_de_etappp","ds_de_etappp",tree_eta,RooArgSet(argset,mpi0),de_etappp_cuts.c_str());

  mean.setVal(1.08281e-03);
  sl.setVal(1.81824e-02);
  sr.setVal(1.47628e-02);
  if(fix_de_etappp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_etappp,Timer(true),Range("de_fit"));
  const double de_etappp_min = mean.getVal() - 3.*sl.getVal();
  const double de_etappp_max = mean.getVal() + 3.*sr.getVal();

  alphal.setVal(9.57376e-01);
  alphar.setVal(-1.29682e+00);
  deCBl.setVal(-1.59062e-02);
  deCBr.setVal(8.58390e-04);
  fCBl.setVal(4.22100e-01);
  fCBr.setVal(5.33223e-01);
  mean.setVal(1.57449e-02);
  nl.setVal(1.85938e+00);
  nr.setVal(3.72925e+00);
  s1.setVal(9.62635e-03);
  sl.setVal(2.25108e-02);
  sr.setVal(1.25048e-02);
  if(fix_de_etappp){
    deCBl.setConstant(kTRUE);
    deCBr.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    alphal.setConstant(kTRUE);
    alphar.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    deCBr.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    alphal.setConstant(kFALSE);
    alphar.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de.fitTo(ds_de_etappp,Timer(true),Range("de_plot"));
  de.setRange("de_etappp_sig",de_etappp_min,de_etappp_max);
  const double int_de_etappp_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_etappp_sig"))->getVal();
  const double int_de_etappp_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_etapppFrame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  ds_de_etappp.plotOn(de_etapppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(de_etapppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_etappp_cm = new TCanvas("de_etappp_cm","de_etappp_cm",600,400);
  de_etappp_cm->cd();
  ds_de_etappp.statOn(de_etapppFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_etapppFrame->GetXaxis()->SetTitleSize(0.05);
  de_etapppFrame->GetXaxis()->SetTitleOffset(0.85);
  de_etapppFrame->GetXaxis()->SetLabelSize(0.05);
  de_etapppFrame->GetYaxis()->SetTitleOffset(1.6);
  de_etapppFrame->Draw();
  TLine *de_etappplineLEFT = new TLine(de_etappp_min,0,de_etappp_min,3200);
  de_etappplineLEFT->SetLineColor(kRed);
  de_etappplineLEFT->SetLineWidth(2);
  de_etappplineLEFT->SetLineStyle(1);
  de_etappplineLEFT->Draw();
  TLine *de_etappplineRIGHT = new TLine(de_etappp_max,0,de_etappp_max,3200);
  de_etappplineRIGHT->SetLineColor(kRed);
  de_etappplineRIGHT->SetLineWidth(2);
  de_etappplineRIGHT->SetLineStyle(1);
  de_etappplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_etappp_sig/int_de_etappp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_etappp_cm->Update();
  de_etappp_cm->Print("pics/de_etappp_cut.eps");
  de_etappp_cm->Print("pics/de_etappp_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 eta->ppp ** //
  }//skip

  if(!skip_mbc_omega){
  // Mbc for D0 omega //
  const string mbc_omega_cuts = omegamodecut + ps + mkcut + ps + deomegacut + ps + bdtomegacut + ps + mh0omegacut + ps + mdcut + ps + mpi0cut;
  RooDataSet ds_mbc_omega("ds_mbc_omega","ds_mbc_omega",tree_omega,RooArgSet(argset,mpi0),mbc_omega_cuts.c_str());

  mean.setVal(5.27987e+00);
  sl.setVal(3.15573e-03);
  sr.setVal(2.51960e-03);
  if(fix_mbc_omega){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_omega,Timer(true),Range("mbc_fit"));
  const double mbc_omega_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_omega_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(9.74025e-02);
  mean.setVal(5.27975e+00);
  s1.setVal(2.85890e-03);
  if(fix_mbc_omega){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_omega,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_omega",mbc_omega_min,mbc_omega_max);
  const double int_mbc_omega_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_omega"))->getVal();
  const double int_mbc_omega_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_omegaFrame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#omega"));
  ds_mbc_omega.plotOn(mbc_omegaFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_omegaFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_omega_cm = new TCanvas("mbc_omega_cm","mbc_omega_cm",600,400);
  mbc_omega_cm->cd();
  ds_mbc_omega.statOn(mbc_omegaFrame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_omegaFrame->GetXaxis()->SetTitleSize(0.05);
  mbc_omegaFrame->GetXaxis()->SetTitleOffset(0.85);
  mbc_omegaFrame->GetXaxis()->SetLabelSize(0.05);
  mbc_omegaFrame->GetYaxis()->SetTitleOffset(1.6);
  mbc_omegaFrame->Draw();
  TLine *mbc_omegalineLEFT = new TLine(mbc_omega_min,0,mbc_omega_min,3200);
  mbc_omegalineLEFT->SetLineColor(kRed);
  mbc_omegalineLEFT->SetLineWidth(2);
  mbc_omegalineLEFT->SetLineStyle(1);
  mbc_omegalineLEFT->Draw();
  TLine *mbc_omegalineRIGHT = new TLine(mbc_omega_max,0,mbc_omega_max,3200);
  mbc_omegalineRIGHT->SetLineColor(kRed);
  mbc_omegalineRIGHT->SetLineWidth(2);
  mbc_omegalineRIGHT->SetLineStyle(1);
  mbc_omegalineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_omega_sig/int_mbc_omega_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_omega_cm->Update();
  mbc_omega_cm->Print("pics/mbc_omega_cut.eps");
  mbc_omega_cm->Print("pics/mbc_omega_cut.root");
  // ** end of mbc for D omega ** //
  }//skip

  if(!skip_de_omega){
  // dE for D0 omega //
  const string de_omega_cuts = omegamodecut + ps + mdcut + ps + mkcut + ps + mbcomegacut + ps + bdtomegacut + ps + mh0omegacut;
  RooDataSet ds_de_omega("ds_de_omega","ds_de_omega",tree_omega,RooArgSet(argset,mpi0),de_omega_cuts.c_str());

  mean.setVal(1.25827e-03);
  sl.setVal(1.92832e-02);
  sr.setVal(1.51722e-02);
  if(fix_de_omega){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_omega,Timer(true),Range("de_fit"));
  const double de_omega_min = mean.getVal() - 3.*sl.getVal();
  const double de_omega_max = mean.getVal() + 3.*sr.getVal();

  alphal.setVal(4.25000e-01);
  alphar.setVal(-1.07006e+00);
  deCBl.setVal(-3.56010e-02);
  deCBr.setVal(2.94853e-03);
  fCBl.setVal(2.35851e-01);
  fCBr.setVal(3.36197e-01);
  mean.setVal(-2.65321e-03);
  nl.setVal(4.42277e+00);
  nr.setVal(4.21416e+00);
  s1.setVal(1.79495e-02);
  sl.setVal(1.96292e-02);
  sr.setVal(1.21546e-02);
  if(fix_de_omega){
    deCBl.setConstant(kTRUE);
    deCBr.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    fCBr.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    nr.setConstant(kTRUE);
    alphal.setConstant(kTRUE);
    alphar.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    deCBr.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    fCBr.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    nr.setConstant(kFALSE);
    alphal.setConstant(kFALSE);
    alphar.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de.fitTo(ds_de_omega,Timer(true),Range("de_plot"));
  de.setRange("de_omega_sig",de_omega_min,de_omega_max);
  const double int_de_omega_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_omega_sig"))->getVal();
  const double int_de_omega_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_omegaFrame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#omega"));
  ds_de_omega.plotOn(de_omegaFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(de_omegaFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_omega_cm = new TCanvas("de_omega_cm","de_omega_cm",600,400);
  de_omega_cm->cd();
  ds_de_omega.statOn(de_omegaFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_omegaFrame->GetXaxis()->SetTitleSize(0.05);
  de_omegaFrame->GetXaxis()->SetTitleOffset(0.85);
  de_omegaFrame->GetXaxis()->SetLabelSize(0.05);
  de_omegaFrame->GetYaxis()->SetTitleOffset(1.6);
  de_omegaFrame->Draw();
  TLine *de_omegalineLEFT = new TLine(de_omega_min,0,de_omega_min,3200);
  de_omegalineLEFT->SetLineColor(kRed);
  de_omegalineLEFT->SetLineWidth(2);
  de_omegalineLEFT->SetLineStyle(1);
  de_omegalineLEFT->Draw();
  TLine *de_omegalineRIGHT = new TLine(de_omega_max,0,de_omega_max,3200);
  de_omegalineRIGHT->SetLineColor(kRed);
  de_omegalineRIGHT->SetLineWidth(2);
  de_omegalineRIGHT->SetLineStyle(1);
  de_omegalineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_omega_sig/int_de_omega_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_omega_cm->Update();
  de_omega_cm->Print("pics/de_omega_cut.eps");
  de_omega_cm->Print("pics/de_omega_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 omega ** //
  }//skip


  if(!skip_mbc_etapgg){
  // Mbc for D0 eta' (eta->gg) //
  RooArgSet argset3(argset);
  argset3.add(dmetap);
  const string mbc_etapgg_cuts = etapggmodecut + ps + mkcut + ps + deetapggcut + ps + bdtetaggcut + ps + mh0etaggcut + ps + mdcut + ps + dmetapggcut;
  RooDataSet ds_mbc_etapgg("ds_mbc_etapgg","ds_mbc_etapgg",tree_etap,argset3,mbc_etapgg_cuts.c_str());

  mean.setVal(5.27985e+00);
  sl.setVal(3.12139e-03);
  sr.setVal(2.49806e-03);
  if(fix_mbc_etapgg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_etapgg,Timer(true),Range("mbc_fit"));
  const double mbc_etapgg_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_etapgg_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(9.11142e-02);
  mean.setVal(5.27972e+00);
  s1.setVal(2.82721e-03);
  if(fix_mbc_etapgg){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_etapgg,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_etapgg",mbc_etapgg_min,mbc_etapgg_max);
  const double int_mbc_etapgg_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_etapgg"))->getVal();
  const double int_mbc_etapgg_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_etapggFrame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#eta` (#eta#rightarrow#gamma#gamma)"));
  ds_mbc_etapgg.plotOn(mbc_etapggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_etapggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_etapgg_cm = new TCanvas("mbc_etapgg_cm","mbc_etapgg_cm",600,400);
  mbc_etapgg_cm->cd();
  ds_mbc_etapgg.statOn(mbc_etapggFrame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_etapggFrame->GetXaxis()->SetTitleSize(0.05);
  mbc_etapggFrame->GetXaxis()->SetTitleOffset(0.85);
  mbc_etapggFrame->GetXaxis()->SetLabelSize(0.05);
  mbc_etapggFrame->GetYaxis()->SetTitleOffset(1.6);
  mbc_etapggFrame->Draw();
  TLine *mbc_etapgglineLEFT = new TLine(mbc_etapgg_min,0,mbc_etapgg_min,100);
  mbc_etapgglineLEFT->SetLineColor(kRed);
  mbc_etapgglineLEFT->SetLineWidth(2);
  mbc_etapgglineLEFT->SetLineStyle(1);
  mbc_etapgglineLEFT->Draw();
  TLine *mbc_etapgglineRIGHT = new TLine(mbc_etapgg_max,0,mbc_etapgg_max,100);
  mbc_etapgglineRIGHT->SetLineColor(kRed);
  mbc_etapgglineRIGHT->SetLineWidth(2);
  mbc_etapgglineRIGHT->SetLineStyle(1);
  mbc_etapgglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_etapgg_sig/int_mbc_etapgg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_etapgg_cm->Update();
  mbc_etapgg_cm->Print("pics/mbc_etapgg_cut.eps");
  mbc_etapgg_cm->Print("pics/mbc_etapgg_cut.root");
  // ** end of mbc for D eta' (eta->gg) ** //
  }//skip

  if(!skip_de_etapgg){
  // dE for D0 eta'(->gg) //
  RooArgSet argset2(argset);
  argset2.add(dmetap);
  const string de_etapgg_cuts = etapggmodecut + ps + mdcut + ps + mkcut + ps + mbcetapggcut + ps + bdtetaggcut + ps + mh0etaggcut + ps + dmetapggcut;
  RooDataSet ds_de_etapgg("ds_de_etapgg","ds_de_etapgg",tree_etap,argset2,de_etapgg_cuts.c_str());

  mean.setVal(6.58973e-03);
  sl.setVal(3.13538e-02);
  sr.setVal(1.91869e-02);
  if(fix_de_etapgg){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_etapgg,Timer(true),Range("de_fit"));
  const double de_etapgg_min = mean.getVal() - 3.*sl.getVal();
  const double de_etapgg_max = mean.getVal() + 3.*sr.getVal();

  deCBl.setVal(-2.47968e-05);
  fCBl.setVal(4.17267e-01);
  mean.setVal(1.13805e-02);
  nl.setVal(9.99987e+01);
  s1.setVal(2.83275e-02);
  sl.setVal(3.49433e-02);
  sr.setVal(1.45921e-02);
  if(fix_de_etapgg){
    deCBl.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de0.fitTo(ds_de_etapgg,Timer(true),Range("de_plot"));
  de.setRange("de_etapgg_sig",de_etapgg_min,de_etapgg_max);
  const double int_de_etapgg_sig  = pdf_m_de0.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_etapgg_sig"))->getVal();
  const double int_de_etapgg_plot = pdf_m_de0.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_etapggFrame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#eta'(#eta#rightarrow#gamma#gamma)"));
  ds_de_etapgg.plotOn(de_etapggFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de0.plotOn(de_etapggFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_etapgg_cm = new TCanvas("de_etapgg_cm","de_etapgg_cm",600,400);
  de_etapgg_cm->cd();
  ds_de_etapgg.statOn(de_etapggFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_etapggFrame->GetXaxis()->SetTitleSize(0.05);
  de_etapggFrame->GetXaxis()->SetTitleOffset(0.85);
  de_etapggFrame->GetXaxis()->SetLabelSize(0.05);
  de_etapggFrame->GetYaxis()->SetTitleOffset(1.6);
  de_etapggFrame->Draw();
  TLine *de_etapgglineLEFT = new TLine(de_etapgg_min,0,de_etapgg_min,200);
  de_etapgglineLEFT->SetLineColor(kRed);
  de_etapgglineLEFT->SetLineWidth(2);
  de_etapgglineLEFT->SetLineStyle(1);
  de_etapgglineLEFT->Draw();
  TLine *de_etapgglineRIGHT = new TLine(de_etapgg_max,0,de_etapgg_max,200);
  de_etapgglineRIGHT->SetLineColor(kRed);
  de_etapgglineRIGHT->SetLineWidth(2);
  de_etapgglineRIGHT->SetLineStyle(1);
  de_etapgglineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_etapgg_sig/int_de_etapgg_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_etapgg_cm->Update();
  de_etapgg_cm->Print("pics/de_etapgg_cut.eps");
  de_etapgg_cm->Print("pics/de_etapgg_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 eta' (eta->gg) ** //
  }//skip

  if(!skip_mbc_etapppp){
  // Mbc for D0 eta' (eta->ppp) //
  RooArgSet argset2(argset);
  argset2.add(mpi0);
  argset2.add(dmetap);
  const string mbc_etapppp_cuts = etappppmodecut + ps + mkcut + ps + deetappppcut + ps + bdtetapppcut + ps + mh0etapppcut + ps + mdcut + ps + mpi0cut + ps + dmetappppcut;
  RooDataSet ds_mbc_etapppp("ds_mbc_etapppp","ds_mbc_etapppp",tree_etap,argset2,mbc_etapppp_cuts.c_str());

  mean.setVal(5.27985e+00);
  sl.setVal(3.12139e-03);
  sr.setVal(2.49806e-03);
  if(fix_mbc_etapppp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_mbc.fitTo(ds_mbc_etapppp,Timer(true),Range("mbc_fit"));
  const double mbc_etapppp_min = mean.getVal() - 2.5*sl.getVal();
  const double mbc_etapppp_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(3.99910e-02);
  mean.setVal(5.27957e+00);
  s1.setVal(2.58833e-03);
  if(fix_mbc_etapppp){
    mean.setConstant(kTRUE);
    s1.setConstant(kTRUE);
    alpha_mbc.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    s1.setConstant(kFALSE);
    alpha_mbc.setConstant(kFALSE);
  }

  RooFitResult* r = pdf_m_mbc.fitTo(ds_mbc_etapppp,Timer(true),Range("mbc_plot"));
  mbc.setRange("mbc_sig_etapppp",mbc_etapppp_min,mbc_etapppp_max);
  const double int_mbc_etapppp_sig = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_etapppp"))->getVal();
  const double int_mbc_etapppp_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* mbc_etappppFrame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D^{0}#eta` (#eta#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  ds_mbc_etapppp.plotOn(mbc_etappppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(mbc_etappppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* mbc_etapppp_cm = new TCanvas("mbc_etapppp_cm","mbc_etapppp_cm",600,400);
  mbc_etapppp_cm->cd();
  ds_mbc_etapppp.statOn(mbc_etappppFrame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  mbc_etappppFrame->GetXaxis()->SetTitleSize(0.05);
  mbc_etappppFrame->GetXaxis()->SetTitleOffset(0.85);
  mbc_etappppFrame->GetXaxis()->SetLabelSize(0.05);
  mbc_etappppFrame->GetYaxis()->SetTitleOffset(1.6);
  mbc_etappppFrame->Draw();
  TLine *mbc_etapppplineLEFT = new TLine(mbc_etapppp_min,0,mbc_etapppp_min,100);
  mbc_etapppplineLEFT->SetLineColor(kRed);
  mbc_etapppplineLEFT->SetLineWidth(2);
  mbc_etapppplineLEFT->SetLineStyle(1);
  mbc_etapppplineLEFT->Draw();
  TLine *mbc_etapppplineRIGHT = new TLine(mbc_etapppp_max,0,mbc_etapppp_max,100);
  mbc_etapppplineRIGHT->SetLineColor(kRed);
  mbc_etapppplineRIGHT->SetLineWidth(2);
  mbc_etapppplineRIGHT->SetLineStyle(1);
  mbc_etapppplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_etapppp_sig/int_mbc_etapppp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  mbc_etapppp_cm->Update();
  mbc_etapppp_cm->Print("pics/mbc_etapppp_cut.eps");
  mbc_etapppp_cm->Print("pics/mbc_etapppp_cut.root");
  // ** end of mbc for D eta' (eta->ppp) ** //
  }//skip

  if(!skip_de_etapppp){
  // dE for D0 eta'(->ppp) //
  RooArgSet argset1(argset);
  argset1.add(mpi0);
  argset1.add(dmetap);
  const string de_etapppp_cuts = etappppmodecut + ps + mdcut + ps + mkcut + ps + mbcetappppcut + ps + bdtetapppcut + ps + mh0etapppcut + ps + mpi0cut + ps + dmetappppcut;
  RooDataSet ds_de_etapppp("ds_de_etapppp","ds_de_etapppp",tree_etap,argset1,de_etapppp_cuts.c_str());

  mean.setVal(-1.71340e-03);
  sl.setVal(1.25908e-02);
  sr.setVal(1.24989e-02);
  if(fix_de_etapppp){
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
  } else{
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
  }
  pdf_de.fitTo(ds_de_etapppp,Timer(true),Range("de_fit"));
  const double de_etapppp_min = mean.getVal() - 3.*sl.getVal();
  const double de_etapppp_max = mean.getVal() + 3.*sr.getVal();

  deCBl.setVal(-1.18815e-02);
  fCBl.setVal(2.78244e-01);
  mean.setVal(-1.40767e-03);
  nl.setVal(3.80990e+00);
  s1.setVal(3.52310e-02);
  sl.setVal(1.11162e-02);
  sr.setVal(1.09681e-02);
  if(fix_de_etapppp){
    deCBl.setConstant(kTRUE);
    mean.setConstant(kTRUE);
    sl.setConstant(kTRUE);
    sr.setConstant(kTRUE);
    fCBl.setConstant(kTRUE);
    nl.setConstant(kTRUE);
    s1.setConstant(kTRUE);
  } else{
    deCBl.setConstant(kFALSE);
    mean.setConstant(kFALSE);
    sl.setConstant(kFALSE);
    sr.setConstant(kFALSE);
    fCBl.setConstant(kFALSE);
    nl.setConstant(kFALSE);
    s1.setConstant(kFALSE);
  }
  RooFitResult* r = pdf_m_de0.fitTo(ds_de_etapppp,Timer(true),Range("de_plot"));
  de.setRange("de_etapppp_sig",de_etapppp_min,de_etapppp_max);
  const double int_de_etapppp_sig  = pdf_m_de0.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_etapppp_sig"))->getVal();
  const double int_de_etapppp_plot = pdf_m_de0.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* de_etappppFrame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D^{0}#eta'(#eta#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  ds_de_etapppp.plotOn(de_etappppFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de0.plotOn(de_etappppFrame,VisualizeError(*r,1),LineWidth(2));

  TCanvas* de_etapppp_cm = new TCanvas("de_etapppp_cm","de_etapppp_cm",600,400);
  de_etapppp_cm->cd();
  ds_de_etapppp.statOn(de_etappppFrame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  de_etappppFrame->GetXaxis()->SetTitleSize(0.05);
  de_etappppFrame->GetXaxis()->SetTitleOffset(0.85);
  de_etappppFrame->GetXaxis()->SetLabelSize(0.05);
  de_etappppFrame->GetYaxis()->SetTitleOffset(1.6);
  de_etappppFrame->Draw();
  TLine *de_etapppplineLEFT = new TLine(de_etapppp_min,0,de_etapppp_min,60);
  de_etapppplineLEFT->SetLineColor(kRed);
  de_etapppplineLEFT->SetLineWidth(2);
  de_etapppplineLEFT->SetLineStyle(1);
  de_etapppplineLEFT->Draw();
  TLine *de_etapppplineRIGHT = new TLine(de_etapppp_max,0,de_etapppp_max,60);
  de_etapppplineRIGHT->SetLineColor(kRed);
  de_etapppplineRIGHT->SetLineWidth(2);
  de_etapppplineRIGHT->SetLineStyle(1);
  de_etapppplineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_etapppp_sig/int_de_etapppp_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  de_etapppp_cm->Update();
  de_etapppp_cm->Print("pics/de_etapppp_cut.eps");
  de_etapppp_cm->Print("pics/de_etapppp_cut.root");
  alphal.setVal( 1.); alphal.setConstant(kTRUE);
  alphar.setVal(-1.); alphar.setConstant(kTRUE);
  // ** end of dE for D0 eta' (eta->ppp) ** //
  }//skip

  if(!skip_m_omega){
  cout << "m(omega)      cuts: " << m_omega_min << " " << m_omega_max << endl;
  cout << " Integral: " << int_m_omega_sig << "/" << int_m_omega_plot << " = " << int_m_omega_sig/int_m_omega_plot << endl;
  }

  if(!skip_m_etagg){
  cout << "m(etagg)     cuts: " << m_etagg_min << " " << m_etagg_max << endl;
  cout << " Integral: " << int_m_etagg_sig << "/" << int_m_etagg_plot << " = " << int_m_etagg_sig/int_m_etagg_plot << endl;
  }

  if(!skip_dm_etapgg){
  cout << "dm(etapgg)   cuts: " << dm_etapgg_min << " " << dm_etapgg_max << endl;
  cout << " Integral: " << int_dm_etapgg_sig << "/" << int_dm_etapgg_plot << " = " << int_dm_etapgg_sig/int_dm_etapgg_plot << endl;
  }

  if(!skip_m_etappp){
  cout << "m(etappp)      cuts: " << m_etappp_min << " " << m_etappp_max << endl;
  cout << " Integral: " << int_m_etappp_sig << "/" << int_m_etappp_plot << " = " << int_m_etappp_sig/int_m_etappp_plot << endl;
  }

  if(!skip_dm_etapppp){
  cout << "dm(etapppp)    cuts: " << dm_etapppp_min << " " << dm_etapppp_max << endl;
  cout << " Integral: " << int_dm_etapppp_sig << "/" << int_dm_etapppp_plot << " = " << int_dm_etapppp_sig/int_dm_etapppp_plot << endl;
  }

  if(!skip_m_pi0){
  cout << "m(pi0)        cuts: " << m_pi0_min << " " << m_pi0_max << endl;
  cout << " Integral: " << int_m_pi0_sig << "/" << int_m_pi0_plot << " = " << int_m_pi0_sig/int_m_pi0_plot << endl;
  }

  if(!skip_m_d0){
  cout << "m(D0)         cuts: " << m_d0_min << " " << m_d0_max << endl;
  cout << " Integral: " << int_m_d0_sig << "/" << int_m_d0_plot << " = " << int_m_d0_sig/int_m_d0_plot << endl;
  }

  if(!skip_mbc_pi0){
  cout << "Mbc(pi0)      cuts: " << mbc_pi0_min << " " << mbc_pi0_max << endl;
  cout << " Integral: " << int_mbc_pi0_sig << "/" << int_mbc_pi0_plot << " = " << int_mbc_pi0_sig/int_mbc_pi0_plot << endl;
  }

  if(!skip_de_pi0){
  cout << "dE(pi0)       cuts: " << de_pi0_min << " " << de_pi0_max << endl;
  cout << " Integral: " << int_de_pi0_sig << "/" << int_de_pi0_plot << " = " << int_de_pi0_sig/int_de_pi0_plot << endl;
  }

  if(!skip_mbc_etagg){
  cout << "Mbc(eta->gg)  cuts: " << mbc_etagg_min << " " << mbc_etagg_max << endl;
  cout << " Integral: " << int_mbc_etagg_sig << "/" << int_mbc_etagg_plot << " = " << int_mbc_etagg_sig/int_mbc_etagg_plot << endl;
  }

  if(!skip_de_etagg){
  cout << "dE(eta->gg)   cuts: " << de_etagg_min << " " << de_etagg_max << endl;
  cout << " Integral: " << int_de_etagg_sig << "/" << int_de_etagg_plot << " = " << int_de_etagg_sig/int_de_etagg_plot << endl;
  }

  if(!skip_mbc_etappp){
  cout << "Mbc(eta->ppp) cuts: " << mbc_etappp_min << " " << mbc_etappp_max << endl;
  cout << " Integral: " << int_mbc_etappp_sig << "/" << int_mbc_etappp_plot << " = " << int_mbc_etappp_sig/int_mbc_etappp_plot << endl;
  }

  if(!skip_de_etappp){
  cout << "dE(eta->ppp)  cuts: " << de_etappp_min << " " << de_etappp_max << endl;
  cout << " Integral: " << int_de_etappp_sig << "/" << int_de_etappp_plot << " = " << int_de_etappp_sig/int_de_etappp_plot << endl;
  }

  if(!skip_mbc_omega){
  cout << "Mbc(omega)    cuts: " << mbc_omega_min << " " << mbc_omega_max << endl;
  cout << " Integral: " << int_mbc_omega_sig << "/" << int_mbc_omega_plot << " = " << int_mbc_omega_sig/int_mbc_omega_plot << endl;
  }

  if(!skip_de_omega){
  cout << "dE(omega)     cuts: " << de_omega_min << " " << de_omega_max << endl;
  cout << " Integral: " << int_de_omega_sig << "/" << int_de_omega_plot << " = " << int_de_omega_sig/int_de_omega_plot << endl;
  }

  if(!skip_mbc_etapgg){
    cout << "mbc(eta')(->gg)  cuts: " << mbc_etapgg_min << " " << mbc_etapgg_max << endl;
    cout << " Integral: " << int_mbc_etapgg_sig << "/" << int_mbc_etapgg_plot << " = " << int_mbc_etapgg_sig/int_mbc_etapgg_plot << endl;
  }

  if(!skip_de_etapgg){
    cout << "dE(eta')(->gg)  cuts: " << de_etapgg_min << " " << de_etapgg_max << endl;
    cout << " Integral: " << int_de_etapgg_sig << "/" << int_de_etapgg_plot << " = " << int_de_etapgg_sig/int_de_etapgg_plot << endl;
  }

  if(!skip_mbc_etapppp){
    cout << "mbc(eta')(->ppp)  cuts: " << mbc_etapppp_min << " " << mbc_etapppp_max << endl;
    cout << " Integral: " << int_mbc_etapppp_sig << "/" << int_mbc_etapppp_plot << " = " << int_mbc_etapppp_sig/int_mbc_etapppp_plot << endl;
  }

  if(!skip_de_etapppp){
    cout << "dE(eta')(->ppp)   cuts: " << de_etapppp_min << " " << de_etapppp_max << endl;
    cout << " Integral: " << int_de_etapppp_sig << "/" << int_de_etapppp_plot << " = " << int_de_etapppp_sig/int_de_etapppp_plot << endl;
  }
  return;
}
