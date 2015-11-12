#include "cuts.h"
using namespace RooFit;

void PurityFit(const int _mode){
  TChain* tree = new TChain("TEvent");
//  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1.root");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_uds_2_12.root");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_charm_2_12.root");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_charged_2_12.root");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_mixed_2_12.root");

  gROOT->ProcessLine(".L pdfs/RooRhoDeltaEPdf.cxx+");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  b0f.defineType("fsr",10);
  b0f.defineType("bad_pi0",5);
  b0f.defineType("rho2",2);
  b0f.defineType("rho3",3);
  b0f.defineType("rho4",4);
  b0f.defineType("rho11",11);
  b0f.defineType("comb",-1);

  RooCategory mode("mode","mode");
  RooCategory h0mode("h0mode","h0mode");

  double BDTG_MIN = 0;
  double BDTG_MAX = 1;
  bool gg_flag = true;
  double Mbc_min;
  double Mbc_max;
  double dE_min;
  double dE_max;
  int m_mode,m_h0mode;
  double mh0_min, mh0_max;
  string label;
  switch(_mode){
  case 1:
    label = string("#pi^{0}");
    BDTG_MIN = bdt_cut_pi0;
    mode.defineType("pi0",1);
    h0mode.defineType("gg",10);
    Mbc_min = mbc_min_pi0;
    Mbc_max = mbc_max_pi0;
    dE_min  = de_min_pi0;
    dE_max  = de_max_pi0;
    m_mode = 1;
    m_h0mode = 10;
    mh0_min = mpi0_min;
    mh0_max = mpi0_max;
    break;
  case 2:
    label = string("#eta#rightarrow#gamma#gamma");
    BDTG_MIN = bdtg_cut_etagg;
    mode.defineType("eta",2);
    h0mode.defineType("gg",10);
    Mbc_min = mbc_min;
    Mbc_max = mbc_max;
    dE_min  = de_min;
    dE_max  = de_max;
    m_mode = 2;
    m_h0mode = 10;
    mh0_min = EtaGGMass-metagg_cut;
    mh0_max = EtaGGMass+metagg_cut;
    break;
  case 3:
    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    BDTG_MIN = bdtg_cut_etappp;
    gg_flag = false;
    mode.defineType("eta",2);
    h0mode.defineType("ppp",20);
    Mbc_min = mbc_min;
    Mbc_max = mbc_max;
    dE_min  = de_min_etappp;
    dE_max  = de_max_etappp;
    m_mode = 2;
    m_h0mode = 20;
    mh0_min = EtaMass-metappp_cut;
    mh0_max = EtaMass+metappp_cut;
    break;
  case 4:
    label = string("#omega");
    BDTG_MIN = bdtg_cut_omega;
    gg_flag = false;
    mode.defineType("omega",3);
    h0mode.defineType("ppp",20);
    Mbc_min = mbc_min_omega;
    Mbc_max = mbc_max_omega;
    dE_min  = de_min_omega;
    dE_max  = de_max_omega;
    m_mode = 3;
    m_h0mode = 20;
    mh0_min = OmegaMass-momega_cut;
    mh0_max = OmegaMass+momega_cut;
    break;
  default:
    return;
  }

  RooArgSet argset;
  argset.add(mode);
  argset.add(h0mode);
  argset.add(b0f);

  RooCategory flv("flv_mc","flv_mc");
  flv.defineType("B0",1);
  flv.defineType("anti-B0",-1);
  argset.add(flv);

  RooCategory bin("bin","bin");
  bin.defineType("1",1); bin.defineType("-1",-1);
  bin.defineType("2",2); bin.defineType("-2",-2);
  bin.defineType("3",3); bin.defineType("-3",-3);
  bin.defineType("4",4); bin.defineType("-4",-4);
  bin.defineType("5",5); bin.defineType("-5",-5);
  bin.defineType("6",6); bin.defineType("-6",-6);
  bin.defineType("7",7); bin.defineType("-7",-7);
  bin.defineType("8",8); bin.defineType("-8",-8);
  argset.add(bin);

  RooSuperCategory binflv("binflv","binflv",RooArgSet(bin,flv));

  const double mbcMin = 5.20;
  const double mbcMax = 5.2885;
  const double deMin = -0.15;
  const double deMax = 0.3;
  const double elliscaleDe  = TMath::Sqrt(4./TMath::Pi());
  const double elliscaleMbc = TMath::Sqrt(4./TMath::Pi());

  RooRealVar mbc_center("mbc_center","mbc_center",0.5*(Mbc_min+Mbc_max),Mbc_min,Mbc_max); mbc_center.setConstant(kTRUE);
  RooRealVar mbc_center_eq("mbc_center_eq","mbc_center_eq",mr_argedge_3-0.5*(Mbc_max-Mbc_min)*elliscaleMbc,Mbc_min,Mbc_max); mbc_center_eq.setConstant(kTRUE);
  RooRealVar de_center("de_center","de_center",0.5*(dE_min+dE_max),dE_min,dE_max); de_center.setConstant(kTRUE);
  RooRealVar mbc_radius("mbc_radius","mbc_radius",0.5*(Mbc_max-Mbc_min)*elliscaleMbc,0,0.5*(mbcMax-mbcMin)); mbc_radius.setConstant(kTRUE);
  RooRealVar de_radius("de_radius","de_radius",0.5*(dE_max-dE_min)*elliscaleDe,0.,0.5*(deMax-deMin)); de_radius.setConstant(kTRUE);
  RooRealVar mbc_radius1("mbc_radius1","mbc_radius1",0.5*(Mbc_max-Mbc_min),0,0.5*(mbcMax-mbcMin)); mbc_radius1.setConstant(kTRUE);
  RooRealVar de_radius1("de_radius1","de_radius1",0.5*(dE_max-dE_min),0.,0.5*(deMax-deMin)); de_radius1.setConstant(kTRUE);

  cout << 0.5*(Mbc_min+Mbc_max) << " " << 0.5*(Mbc_max-Mbc_min) << endl;
  cout << 0.5*(dE_min+dE_max) << " " << 0.5*(dE_max-dE_min) << endl;

  mbc_center.Print();
  mbc_center_eq.Print();

  RooRealVar mbc("mbc","M_{bc}",0.5*(Mbc_min+Mbc_max),mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",Mbc_min,Mbc_max);
  mbc.setRange("mbcSignal",Mbc_min,Mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);

  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",dE_min,dE_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",dE_min,dE_max);

//   de.setRange("Ellips",dE_min,dE_max);
//   RooFormulaVar mbclo("mbclo","@1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.00001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
//   RooFormulaVar mbchi("mbchi","@1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.00001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
//   mbc.setRange("Ellips",mbclo,mbchi);

  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mh0("mh0","mh0",mh0_min,mh0_max,"GeV"); argset.add(mh0);
  RooRealVar mpi0("mpi0","mpi0",mpi0_min,mpi0_max,"GeV"); if(_mode!=2) argset.add(mpi0);
  RooRealVar bdt("bdt","bdt",BDTG_MIN,BDTG_MAX); argset.add(bdt);

  argset.add(b0f);
  RooDataSet ds_sig("ds_sig","ds_sig",tree,argset,"mbc>0||mbc<=0 && (b0f == 1 || b0f == 5 || b0f == 10)");
  RooDataSet ds_bkg("ds_bkg","ds_bkg",tree,argset,"mbc>0||mbc<=0 && !(b0f == 1 || b0f == 5 || b0f == 10)");
  RooDataHist dh("dh","dh");
  dh.add(ds_sig,"",1./0.563);
  dh.add(ds_bkg,"",1./0.949);

  stringstream out;
  out.str("");
  out << "de<" << dE_max << " && de>" << dE_min;
  out << " && mbc>" << Mbc_min << " && mbc<" << Mbc_max;
  Roo1DTable* sigtable = ds.table(b0f,out.str().c_str());
  sigtable->Print();
  sigtable->Print("v");

  Roo1DTable* fulltable = ds.table(b0f);
  fulltable->Print();
  fulltable->Print("v");

//  RooDataHist* dh = ds0->binnedClone();

  ds.Print();
  int _b0f = -1;
  ////////////////
  // Signal PDF //
  ////////////////
  ////////////
  // de pdf //
  ////////////
  if(gg_flag){
    RooRealVar  de0("de0","de0",get_de0(m_mode,m_h0mode,_b0f),-0.2,0.1); if(cSig) de0.setConstant(kTRUE);
    RooRealVar  s1("s1","s1",get_s1(m_mode,m_h0mode,_b0f),0.,0.5);       if(cSig) s1.setConstant(kTRUE);
    RooGaussian g1("g1","g1",de,de0,s1);

    RooRealVar deCBl("deCBl","deCBl",get_deCBl(m_mode,m_h0mode,_b0f),-0.2,0.1); if(cSig) deCBl.setConstant(kTRUE);
    RooRealVar sCBl("sCBl","sCBl",get_sCBl(m_mode,m_h0mode,_b0f),0.,0.5);       if(cSig) sCBl.setConstant(kTRUE);
    RooRealVar alphal("alphal","alphal", get_alphal(m_mode,m_h0mode,_b0f), 0.,10.);                          if(cSIG) alphal.setConstant(kTRUE);
    RooRealVar nl("nl","nl",2.,0.,100.); nl.setConstant(kTRUE);

    RooRealVar deCBr("deCBr","deCBr",get_deCBr(m_mode,m_h0mode,_b0f),-0.2,0.1); if(cSig) deCBr.setConstant(kTRUE);
    RooRealVar sCBr("sCBr","sCBr",get_sCBr(m_mode,m_h0mode,_b0f),0.,0.5);       if(cSig) sCBr.setConstant(kTRUE);
    RooRealVar alphar("alphar","alphar",get_alphar(m_mode,m_h0mode,_b0f),-10.,0.);                          if(cSig) alphar.setConstant(kTRUE);
    RooRealVar nr("nr","nr",2,0.,100.); nr.setConstant(kTRUE);

    RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
    RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

    RooRealVar fCBl("fCBl","fCBl",get_fCBl(m_mode,m_h0mode,_b0f),0.,1.); if(cSig) fCBl.setConstant(kTRUE);
    RooRealVar fCBr("fCBr","fCBr",get_fCBr(m_mode,m_h0mode,_b0f),0.,1.); if(cSig) fCBr.setConstant(kTRUE);

    RooAddPdf pdf_de_sig("pdf_de_sig","pdf_de_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  } else{
    RooRealVar de0_201("de0_201","de0_201",get_de0(m_mode,m_h0mode,1),-0.1,0.1); if(cSig) de0_201.setConstant(kTRUE);
    RooRealVar  s1_201("s1_201","s1_201",get_s1(m_mode,m_h0mode,1),0.,0.5); if(cSig) s1_201.setConstant(kTRUE);
    RooGaussian g1_201("g1_201","g1_201",de,de0_201,s1_201);

    RooRealVar deCBl_201("deCBl_201","deCBl_201",get_deCBl(m_mode,m_h0mode,1),-0.1,0.1);     if(cSig) deCBl_201.setConstant(kTRUE);
    RooRealVar sCBl_201("sCBl_201","sCBl_201",get_sCBl(m_mode,m_h0mode,1),0.,0.5);           if(cSig) sCBl_201.setConstant(kTRUE);
    RooRealVar nl_201("nl_201","nl_201",2.,0.,100.); nl_201.setConstant(kTRUE);
    RooRealVar alphal_201("alphal_201","alphal_201",get_alphal(m_mode,m_h0mode,1),-10.,10.); if(cSig) alphal_201.setConstant(kTRUE);
    RooRealVar deCBr_201("deCBr_201","deCBr_201",get_deCBr(m_mode,m_h0mode,1),-0.1,0.1);     if(cSig) deCBr_201.setConstant(kTRUE);
    RooRealVar sCBr_201("sCBr_201","sCBr_201",get_sCBr(m_mode,m_h0mode,1),0.,0.5);           if(cSig) sCBr_201.setConstant(kTRUE);
    RooRealVar nr_201("nr_201","nr_201",2.,0.,100.); nr_201.setConstant(kTRUE);
    RooRealVar alphar_201("alphar_201","alphar_201",get_alphar(m_mode,m_h0mode,1),-10.,10.); if(cSig) alphar_201.setConstant(kTRUE);

    RooCBShape CBl_201("CBl_201","CBl_201",de,deCBl_201,sCBl_201,alphal_201,nl_201);
    RooCBShape CBr_201("CBr_201","CBr_201",de,deCBr_201,sCBr_201,alphar_201,nr_201);

    RooRealVar fCBl_201("fCBl_201","fCBl_201",get_fCBl(m_mode,m_h0mode,1),0.,1.); if(cSig) fCBl_201.setConstant(kTRUE);
    if(_mode == 3){
      fCBl_201.setVal(0.);
      fCBl_201.setConstant(kTRUE);
      alphal_201.setConstant(kTRUE);
    }
    RooRealVar fCBr_201("fCBr_201","fCBr_201",get_fCBr(m_mode,m_h0mode,1),0.,1.); if(cSig) fCBr_201.setConstant(kTRUE);

    RooAddPdf pdf_de1("pdf_de1","pdf_de1",RooArgList(CBl_201,CBr_201,g1_201),RooArgSet(fCBl_201,fCBr_201));

    RooRealVar  de0_205("de0_205","de0_205",get_de0(m_mode,m_h0mode,5),-0.2,0.1); if(cSig) de0_205.setConstant(kTRUE);
    RooRealVar  s1_205("s1_205","s1_205",get_s1(m_mode,m_h0mode,5),0.,0.5);       if(cSig) s1_205.setConstant(kTRUE);
    RooGaussian g1_205("g1_205","g1_205",de,de0_205,s1_205);

    RooRealVar deCBl_205("deCBl_205","deCBl_205",get_deCBl(m_mode,m_h0mode,5),-0.1,0.1); if(cSig) deCBl_205.setConstant(kTRUE);
    RooRealVar sCBl_205("sCBl_205","sCBl_205",get_sCBl(m_mode,m_h0mode,5),0.,0.5);       if(cSig) sCBl_205.setConstant(kTRUE);
    RooRealVar nl_205("nl_205","nl_205",2,0.,100.); nl_205.setConstant(kTRUE);
    RooRealVar alphal_205("alphal_205","alphal_205",get_alphal(m_mode,m_h0mode,5),-10.,10.);  if(cSig) alphal_205.setConstant(kTRUE);
    RooCBShape CBl_205("CBl_205","CBl_205",de,deCBl_205,sCBl_205,alphal_205,nl_205);

    RooRealVar fCBl_205("fCBl_205","fCBl_205",get_fCBl(m_mode,m_h0mode,5),0.,1.); if(cSig) fCBl_205.setConstant(kTRUE);

    RooAddPdf pdf_de5("pdf_de5","pdf_de5",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));
  }

  /////////////
  // mbc pdf //
  /////////////
  if(gg_flag){
    RooRealVar a_s("a_s","a_s",get_a_s(_mode)); if(cSig) a_s.setConstant(kTRUE);
    RooRealVar b_s("b_s","b_s",get_b_s(_mode)); if(cSig) b_s.setConstant(kTRUE);
    RooRealVar c_s("c_s","c_s",get_c_s(_mode),0.0015,0.0035);// if(cSig) c_s.setConstant(kTRUE);
    RooFormulaVar S("S","S","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s,b_s,a_s));

    RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);

    RooRealVar a_mbc0("a_mbc0","a_mbc0",get_a_mbc0(_mode)); if(cSig) a_mbc0.setConstant(kTRUE);
    RooRealVar b_mbc0("b_mbc0","b_mbc0",get_b_mbc0(_mode)); if(cSig) b_mbc0.setConstant(kTRUE);
    RooRealVar c_mbc0("c_mbc0","c_mbc0",get_c_mbc0(_mode),5.277,5.285);// if(cSig) c_mbc0.setConstant(kTRUE);
    RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));
    RooNovosibirsk pdf_mbc_sig("pdf_mbc_sig","pdf_mbc_sig",mbc,MBC0,S,alpha);
  } else{
    RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
    RooRealVar c0("c0","c0",get_c0(_mode)); if(cSig) c0.setConstant(kTRUE);
    RooRealVar c1("c1","c1",get_c1(_mode)); if(cSig) c1.setConstant(kTRUE);
    RooRealVar c2("c2","c2",get_c2(_mode)); if(cSig) c2.setConstant(kTRUE);
    RooRealVar mbc0("mbc0","mbc0",5.284,5.277,5.29);// if(cSig) mbc0.setConstant(kTRUE);
    RooFormulaVar MBC("MBC","MBC","@0+@1*TMath::Erf((@2-@3))/@4",RooArgList(mbc0,c0,c1,de,c2));

    RooRealVar a_s1("a_s1","a_s1",get_a_s(_mode),0.15,0.45); if(cSig) a_s1.setConstant(kTRUE);
    RooRealVar b_s1("b_s1","b_s1",get_b_s(_mode),-0.05,0.05); if(cSig) b_s1.setConstant(kTRUE);
    RooRealVar c_s1("c_s1","c_s1",get_c_s(_mode),0.0015,0.0035);// if(cSig) c_s1.setConstant(kTRUE);
    RooFormulaVar S1("S1","S1","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s1,b_s1,a_s1));
    RooNovosibirsk pdf_mbc1("pdf_mbc1","pdf_mbc1",mbc,MBC,S1,alpha);

    RooRealVar a_s5("a_s5","a_s5",get_a5_s(_mode)); if(cSig) a_s5.setConstant(kTRUE);
    RooRealVar b_s5("b_s5","b_s5",get_b5_s(_mode)); if(cSig) b_s5.setConstant(kTRUE);
    RooRealVar c_s5("c_s5","c_s5",get_c5_s(_mode),0.0015,0.0055); if(cSig) c_s5.setConstant(kTRUE);
    RooFormulaVar S5("S5","S5","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s5,b_s5,a_s5));

    RooRealVar a_mbc0("a_mbc0","a_mbc0",get_a5_mbc0(_mode)); if(cSig) a_mbc0.setConstant(kTRUE);
    RooRealVar b_mbc0("b_mbc0","b_mbc0",get_b5_mbc0(_mode)); if(cSig) b_mbc0.setConstant(kTRUE);
    RooRealVar c_mbc0("c_mbc0","c_mbc0",get_c5_mbc0(_mode),5.27,5.29); if(cSig) c_mbc0.setConstant(kTRUE);
    RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));
    RooNovosibirsk pdf_mbc5("pdf_mbc5","pdf_mbc5",mbc,MBC0,S5,alpha);
  }

  /////////
  // pdf //
  /////////
  if(gg_flag){
    RooProdPdf pdf_sig("pdf_sig","pdf_sig",pdf_de_sig,Conditional(pdf_mbc_sig,mbc));
  } else{
    RooRealVar f_201("f_201","f_201",get_f201(m_mode,m_h0mode),0.,1.); if(cSig) f_201.setConstant(kTRUE);
    RooProdPdf pdf1_sig("pdf1_sig","pdf1_sig",pdf_de1,Conditional(pdf_mbc1,mbc));
    RooProdPdf pdf5_sig("pdf5_sig","pdf5_sig",pdf_de5,Conditional(pdf_mbc5,mbc));
    RooAddPdf  pdf_sig("pdf_sig","pdf_sig",RooArgList(pdf1_sig,pdf5_sig),RooArgSet(f_201));
  }

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar c10("c10","c10",get_cmb_c10(_mode),-10,50.); if(cComb) c10.setConstant(kTRUE);
  RooRealVar c11("c11","c11",get_cmb_c11(_mode),-50,0.);  if(cComb) c11.setConstant(kTRUE);
  RooFormulaVar c1_cmb("c1_cmb","@0+@1*@2",RooArgSet(c10,c11,mbc));
  RooRealVar c2_cmb("c2_cmb","c2_cmb",get_cmb_c20(_mode),-0.1,1);     if(cComb) c2_cmb.setConstant(kTRUE);
  RooChebychev pdf_de_comb_bb("pdf_de_comb_bb","pdf_de_comb_bb",de,RooArgSet(c1_cmb,c2_cmb));

  RooRealVar C1("C1","C1",get_cmb_c1(_mode),-10,50.); if(cComb) C1.setConstant(kTRUE);
  RooRealVar C2("C2","C2",get_cmb_c2(_mode),-0.1,1);  if(cComb) C2.setConstant(kTRUE);
  RooChebychev pdf_de_comb_qq("pdf_de_comb_qq","pdf_de_comb_qq",de,RooArgSet(C1,C2));
  /////////////
  // mbc pdf //
  /////////////
  RooRealVar argedge("argedge","argedge",5.288,5.285,5.29); //argedge.setConstant(kTRUE);
  RooRealVar argpar_cmb_bb("argpar_cmb_bb","argpar_cmb_bb",get_argpar_bb(_mode),-300,-10.); if(cComb) argpar_cmb_bb.setConstant(kTRUE);
  RooArgusBG pdf_mbc_comb_ar("pdf_mbc_comb_ar","Argus PDF",mbc,argedge,argpar_cmb_bb);

  RooRealVar mbc0_cmb_bb("mbc0_cmb_bb","mbc0_cmb_bb",get_mbc0_cmb_bb(_mode),5.25,5.29,"GeV");// if(cComb) mbc0_cmb_bb.setConstant(kTRUE);
  RooRealVar mbcWidth_cmb_bb("mbcWidth","mbcWidth",get_mbcw_cmb_bb(_mode),0.,0.1,"GeV"); if(cComb) mbcWidth_cmb_bb.setConstant(kTRUE);
  RooGaussian mbcGaus_cmb_bb("mbcGaus","mbcGaus",mbc,mbc0_cmb_bb,mbcWidth_cmb_bb);

  RooRealVar f_g("f_g","f_g",get_f_g_cmb_bb(_mode),0.4,0.7);if(_mode == 2 || !gg_flag){ f_g.setConstant(kTRUE);}
  RooAddPdf pdf_mbc_cmb_bb("pdf_mbc_cmb_bb","pdf_mbc_cmb_bb",RooArgList(mbcGaus_cmb_bb,pdf_mbc_comb_ar),RooArgSet(f_g));

  RooRealVar argpar_cmb_qq("argpar_cmb_qq","argpar_cmb_qq",get_argpar_qq(_mode),-300,-10.); if(cComb) argpar_cmb_qq.setConstant(kTRUE);
  RooArgusBG pdf_mbc_cmb_qq("pdf_mbc_cmb_qq","pdf_mbc_cmb_qq",mbc,argedge,argpar_cmb_qq);
  
  /////////
  // pdf //
  /////////
  RooRealVar f_bb("f_bb","f_bb",0.3,0.,1.);
  RooProdPdf pdf_cmb_bb("pdf_cmb_bb","pdf_cmb_bb",pdf_mbc_cmb_bb,Conditional(pdf_de_comb_bb,de));
  RooProdPdf pdf_cmb_qq("pdf_cmb_qq","pdf_cmb_qq",pdf_mbc_cmb_qq,Conditional(pdf_de_comb_qq,de));
  RooAddPdf pdf_comb("pdf_comb","pdf_comb",RooArgSet(pdf_cmb_bb,pdf_cmb_qq),RooArgList(f_bb));

  /////////////////////
  // Peaking bkg PDF //
  /////////////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar de0r("de0r","de0r",get_de0r(_mode),-0.2,0.12);         if(cPeak) de0r.setConstant(kTRUE);
  RooRealVar slopel("slopel","slopel",get_slopel(_mode),-1.e5,0.);  if(cPeak) slopel.setConstant(kTRUE);
  RooRealVar sloper("sloper","sloper",get_sloper(_mode),-10000,0.); if(cPeak) sloper.setConstant(kTRUE);
  RooRealVar steep("steep","steep",get_steep(_mode),0.,1000.);      if(cPeak) steep.setConstant(kTRUE);
  RooRealVar p5("p5","p5",get_p5(_mode),0.01,1000.);                if(cPeak) p5.setConstant(kTRUE);
  RooRhoDeltaEPdf pdf_de_peak("pdf_de_peak","pdf_de_peak",de,de0r,slopel,sloper,steep,p5);
//  RooGenericPdf pdf_de_peak("pdf_de_peak","1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) > 0 ? 1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) : 0.001",RooArgSet(de,de0r,slopel,sloper,steep,p5));
  /////////////
  // mbc pdf //
  /////////////
  if(gg_flag){
    RooRealVar b_peak_s("b_peak_s","b_peak_s",get_peak_b_s(_mode),-0.1,0.1); if(cPeak) b_peak_s.setConstant(kTRUE);
    RooRealVar k_peak_s("k_peak_s","k_peak_s",get_peak_k_s(_mode),-0.1,0.1); if(cPeak) k_peak_s.setConstant(kTRUE);
    RooFormulaVar S_peak("S_peak","S_peak","@0+@1*@2",RooArgList(b_peak_s,de,k_peak_s));
    RooRealVar alpha_peak("alpha_peak","alpha_peak",0.139,0.01,2.); alpha_peak.setConstant(kTRUE);
    RooRealVar b_peak_mbc0("b_peak_mbc0","b_peak_mbc0",get_peak_b_mbc0(_mode),5.25,5.29); if(cPeak) b_peak_mbc0.setConstant(kTRUE);
    RooRealVar k_peak_mbc0("k_peak_mbc0","k_peak_mbc0",get_peak_k_mbc0(_mode),-0.1,0.1);  if(cPeak) k_peak_mbc0.setConstant(kTRUE);
    RooFormulaVar MBC0_peak("MBC0_peak","MBC0_peak","@0+@1*@2",RooArgList(b_peak_mbc0,de,k_peak_mbc0));
    RooNovosibirsk pdf_mbc_peak("pdf_mbc_peak","pdf_mbc_peak",mbc,MBC0_peak,S_peak,alpha_peak);
  } else{
//    RooRealVar argedge("argedge","argedge",5.288,5.285,5.29); //argedge.setConstant(kTRUE);
    RooRealVar argpar_peak_bb("argpar_peak_bb","argpar_peak_bb",get_argpar_bb(_mode),-300,-10.); if(cPeak) argpar_peak_bb.setConstant(kTRUE);
    RooArgusBG pdf_mbc_peak_ar("pdf_mbc_peak_ar","Argus PDF",mbc,argedge,argpar_peak_bb);

    RooRealVar mbc0_peak("mbc0_peak","mbc0_peak",get_peak_b_mbc0(_mode),5.25,5.291,"GeV"); if(cPeak) mbc0_peak.setConstant(kTRUE);
    RooRealVar mbcWidth_peak("mbcWidth_peak","mbcWidth_peak",get_peak_b_s(_mode),0.,0.1,"GeV"); if(cPeak) mbcWidth_peak.setConstant(kTRUE);
    RooGaussian mbcGaus_peak("mbcGaus_peak","mbcGaus_peak",mbc,mbc0_peak,mbcWidth_peak);
    RooRealVar f_g_peak("f_g_peak","f_g_peak",get_f_g_cmb_bb(_mode),0.,1.); if(cPeak) f_g_peak.setConstant(kTRUE);

    RooAddPdf pdf_mbc_peak("pdf_mbc_peak","pdf_mbc_peak",RooArgList(mbcGaus_peak,pdf_mbc_peak_ar),RooArgSet(f_g_peak));
  }
  /////////
  // pdf //
  /////////
  RooProdPdf pdf_peak("pdf_peak","pdf_peak",pdf_de_peak,Conditional(pdf_mbc_peak,mbc));

  //////////////////
  // Complete PDF //
  //////////////////
  RooRealVar Nsig("Nsig","Nsig",1150,0.,10000.);
//  RooRealVar Npbg("Npbg","Npbg",100,0,100000.);
  RooRealVar Ncmb("Ncmb","Ncmb",2288,0,100000);
  switch(_mode){
  case 1:
    RooRealVar Npbg("Npbg","Npbg",100,0,100000.);
    break;
  case 2:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0051);
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  case 3:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0081);
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  case 4:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0031);
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  default:
    return -1;
  }

  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_peak,pdf_comb),RooArgList(Nsig,Npbg,Ncmb));

  RooArgSet* params = pdf.getParameters(RooArgSet(de,mbc));
//  RooArgset* initParams = (RooArgSet*) params->snapshot();

  pdf.fitTo(ds,Verbose(),Timer(true));

  params->printLatex(OutputFile("PurityFit.tex"));

   RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intRho  = pdf_peak->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   const double nsig = intSig->getVal()*Nsig.getVal();
   const double nsig_err = intSig->getVal()*Nsig.getError();
   const double nsig_err_npq = TMath::Sqrt(nsig*(Nsig.getVal()-nsig)/Nsig.getVal());
   const double nsig_err_total = TMath::Sqrt(nsig_err*nsig_err+nsig_err_npq*nsig_err_npq);
   const double nrho = intRho->getVal()*Npbg.getVal();
   const double nrho_err = _mode == 1 ? intRho->getVal()*Npbg.getError() : intRho->getVal()*f_bb.getError()*Ncmb.getVal()*f_p_f_bbc.getVal();
   const double nrho_err_npq = TMath::Sqrt(nrho*(Npbg.getVal()-nrho)/Npbg.getVal());
   const double nrho_err_total = TMath::Sqrt(nrho_err*nrho_err+nrho_err_npq*nrho_err_npq);
   const double ncmb = intCmb->getVal()*Ncmb.getVal();
   const double ncmb_err = intCmb->getVal()*Ncmb.getError();
   const double ncmb_err_npq = TMath::Sqrt(ncmb*(Ncmb.getVal()-ncmb)/Ncmb.getVal());
   const double ncmb_err_total = TMath::Sqrt(ncmb_err*ncmb_err+ncmb_err_npq*ncmb_err_npq);
   const double purity = nsig/(nsig+nrho+ncmb);
   const double purity_err = nsig_err_total/(nsig+nrho+ncmb);

   de.setRange("Ellips",dE_min,dE_max);
   RooFormulaVar mbclo("mbclo","@1-@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   RooFormulaVar mbchi("mbchi","@1+@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   mbc.setRange("Ellips",mbclo,mbchi);

   de.setRange("Elli",dE_min,dE_max);
   RooFormulaVar mbclo1("mbclo1","@1-@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   RooFormulaVar mbchi1("mbchi1","@1+@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   mbc.setRange("Elli",mbclo1,mbchi1);

   RooAbsReal* intSigEl = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intRhoEl = pdf_peak->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intCmbEl = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   const double nsigEl = intSigEl->getVal()*Nsig.getVal();
   const double nsig_errEl = intSigEl->getVal()*Nsig.getError();
   const double nsig_errEl_npq = TMath::Sqrt(fabs(nsigEl*(Nsig.getVal()-nsigEl)/Nsig.getVal()));
   const double nsig_errEl_total = TMath::Sqrt(fabs(nsig_errEl*nsig_errEl+nsig_errEl_npq*nsig_errEl_npq));
   const double nrhoEl = intRhoEl->getVal()*Npbg.getVal();
   const double nrho_errEl = _mode == 1 ? intRhoEl->getVal()*Npbg.getError() : intRhoEl->getVal()*f_bb.getError()*Ncmb.getVal()*f_p_f_bbc.getVal();
   const double nrho_errEl_npq = TMath::Sqrt(fabs(nrhoEl*(Npbg.getVal()-nrhoEl)/Npbg.getVal()));
   const double nrho_errEl_total = TMath::Sqrt(fabs(nrho_errEl*nrho_errEl+nrho_errEl_npq*nrho_errEl_npq));
   const double ncmbEl = intCmbEl->getVal()*Ncmb.getVal();
   const double ncmb_errEl = intCmbEl->getVal()*Ncmb.getError();
   const double ncmb_errEl_npq = TMath::Sqrt(fabs(ncmbEl*(Ncmb.getVal()-ncmbEl)/Ncmb.getVal()));
   const double ncmb_errEl_total = TMath::Sqrt(ncmb_errEl*ncmb_errEl+ncmb_errEl_npq*ncmb_errEl_npq);
   const double purityEl = nsigEl/(nsigEl+nrhoEl+ncmbEl);
   const double purity_errEl = nsig_errEl_total/(nsigEl+nrhoEl+ncmbEl);


   RooAbsReal* intSigEl1 = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   const double intElli = intSigEl1->getVal();
   RooAbsReal* intRhoEl1 = pdf_peak->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   RooAbsReal* intCmbEl1 = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   const double nsigEl1 = intSigEl1->getVal()*Nsig.getVal();
   const double nsig_errEl1 = intSigEl1->getVal()*Nsig.getError();
   const double nsig_errEl1_npq = TMath::Sqrt(fabs(nsigEl1*(Nsig.getVal()-nsigEl1)/Nsig.getVal()));
   const double nsig_errEl1_total = TMath::Sqrt(fabs(nsig_errEl1*nsig_errEl1+nsig_errEl1_npq*nsig_errEl1_npq));
   const double nrhoEl1 = intRhoEl1->getVal()*Npbg.getVal();
   const double nrho_errEl1 = _mode == 1 ? intRhoEl1->getVal()*Npbg.getError() : intRhoEl1->getVal()*f_bb.getError()*Ncmb.getVal()*f_p_f_bbc.getVal();
   const double nrho_errEl1_npq = TMath::Sqrt(fabs(nrhoEl1*(Npbg.getVal()-nrhoEl1)/Npbg.getVal()));
   const double nrho_errEl1_total = TMath::Sqrt(fabs(nrho_errEl1*nrho_errEl1+nrho_errEl1_npq*nrho_errEl1_npq));
   const double ncmbEl1 = intCmbEl1->getVal()*Ncmb.getVal();
   const double ncmb_errEl1 = intCmbEl1->getVal()*Ncmb.getError();
   const double ncmb_errEl1_npq = TMath::Sqrt(ncmbEl1*(Ncmb.getVal()-ncmbEl1)/Ncmb.getVal());
   const double ncmb_errEl1_total = TMath::Sqrt(ncmb_errEl1*ncmb_errEl1+ncmb_errEl1_npq*ncmb_errEl1_npq);
   const double purityEl1 = nsigEl1/(nsigEl1+nrhoEl1+ncmbEl1);
   const double purity_errEl1 = nsig_errEl1_total/(nsigEl1+nrhoEl1+ncmbEl1);

  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(pdf_peak),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("mbcSignal"));
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

  out.str("");
  out << "(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius1.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius1.getVal() << "<1";
  Roo1DTable* ellitable1 = ds.table(b0f,out.str().c_str());

  out.str("");
  out << "(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius.getVal() << "<1";
  Roo1DTable* ellitable = ds.table(b0f,out.str().c_str());

  const int NSIGNAL_ELLI   = ellitable1->get("signal") + ellitable1->get("fsr") + ellitable1->get("bad_pi0");
//  const int NSIGNAL_ELLIPS = ellitable->get("signal")  + ellitable->get("fsr")  + ellitable->get("bad_pi0");

  stringstream out1;
  TPaveText *pt = new TPaveText(0.5,0.6,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out1.str("");
  out1 << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out1.str().c_str());
  out1.str("");
  out1 << "S: " << (int)(nsigEl1+0.5) << " #pm " << (int)(nsig_errEl1_total+0.5) << " (" << NSIGNAL_ELLI << ")";
//  out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  pt->AddText(out1.str().c_str());
//  out1.str("");
//  out1 << "S_{2}: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5) << " (" << NSIGNAL_ELLIPS << ")";
//  pt->AddText(out1.str().c_str());
  out1.str("");
//  out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  out1 << "P: " << std::fixed << std::setprecision(2) << purityEl1*100. << " #pm " << purity_errEl1*100;
  pt->AddText(out1.str().c_str());
  pt->AddText(label.c_str());
  pt->Draw();

  TLine *de_line_RIGHT;
  de_line_RIGHT = new TLine(dE_max,0,dE_max,120);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  de_line_LEFT = new TLine(dE_min,0,dE_min,120);
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
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_sig),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_peak),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull  = mbc.frame(Title("#Delta E pull distribution"));
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

  TPaveText *ptmbc = new TPaveText(0.2,0.6,0.7,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out1.str("");
  out1 << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out1.str().c_str());
  out1.str("");
  out1 << "S: " << (int)(nsigEl1+0.5) << " #pm " << (int)(nsig_errEl1_total+0.5) << " (" << NSIGNAL_ELLI << ")";
  ptmbc->AddText(out1.str().c_str());
//  out1.str("");
//  out1 << "S_{2}: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5) << " (" << NSIGNAL_ELLIPS << ")";
//  ptmbc->AddText(out1.str().c_str());
  out1.str("");
  out1 << "P: " << std::fixed << std::setprecision(2) << purityEl1*100. << " #pm " << purity_errEl1*100;
  ptmbc->AddText(out1.str().c_str());
  ptmbc->AddText(label.c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT;
  mbc_line_RIGHT = new TLine(Mbc_max,0,Mbc_max,40);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT;
  mbc_line_LEFT = new TLine(Mbc_min,0,Mbc_min,40);
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

  double DEMIN = -0.15;
  if(keysflag) DEMIN = -0.3;
  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,DEMIN,0.1),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");

  // Show signal ranges
  TEllipse* elli = new TEllipse(de_center.getVal(),mbc_center.getVal(),de_radius.getVal(),mbc_radius.getVal());
  elli->SetFillColor(0);
  elli->SetFillStyle(0);
  elli->SetLineColor(kBlue);
  elli->SetLineWidth(2);
  TEllipse* elli1 = new TEllipse(de_center.getVal(),mbc_center.getVal(),de_radius1.getVal(),mbc_radius1.getVal());
  elli1->SetFillColor(0);
  elli1->SetFillStyle(0);
//  elli1->SetLineColor(kBlue);
  elli1->SetLineColor(kRed);
  elli1->SetLineWidth(2);
  TLine* l1 = new TLine(dE_min,Mbc_min,dE_max,Mbc_min);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->SetLineWidth(2);
  TLine* l2 = new TLine(dE_min,Mbc_max,dE_max,Mbc_max);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(1);
  l2->SetLineWidth(2);
  TLine* l3 = new TLine(dE_min,Mbc_min,dE_min,Mbc_max);
  l3->SetLineColor(kRed);
  l3->SetLineStyle(1);
  l3->SetLineWidth(2);
  TLine* l4 = new TLine(dE_max,Mbc_min,dE_max,Mbc_max);
  l4->SetLineColor(kRed);
  l4->SetLineStyle(1);
  l4->SetLineWidth(2);

  TCanvas* ellican = new TCanvas("ellican","ellican",400,400);
  ellican->cd();
  out.str("");
  out << "bdtg>" << BDTG_MIN << " && de>-0.15 && de<0.20 && mbc>5.265 && b0f != 1 && b0f != 5 && b0f != 10 && b0f != 0";
  tree->Draw("mbc:de",out.str().c_str());
  tree->SetMarkerStyle(6);
  tree->SetMarkerColor(kBlue);
  out.str("");
  out << "bdtg>" << BDTG_MIN << " && de>-0.15 && de<0.20 && mbc>5.265 && (b0f == 1 || b0f == 5 || b0f == 10)";
  tree->Draw("mbc:de",out.str().c_str(),"same");
//  elli->Draw();
  elli1->Draw();
//  ellican->Pad().GetXaxis()->SetTitle("#DeltaE (GeV)");
//  l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();

//  TCanvas* sigcan = new TCanvas("sigcan","sigcan",400,400);
//  sigcan->cd();
//  out <<
//  tree->Draw("mbc:de","bdtg>0.98 && de>-0.15 && de<0.20 && mbc>5.265 && (b0f == 1 || b0f == 5 || b0f == 10)");
//  elli->Draw(); elli1->Draw(); l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();

//  TCanvas* backcan = new TCanvas("backcan","backcan",400,400);
//  backcan->cd();
//  tree->Draw("mbc:de","bdtg>0.98 && de>-0.15 && de<0.20 && mbc>5.265 && !(b0f == 1 || b0f == 5 || b0f == 10)");
//  elli->Draw(); elli1->Draw(); l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();

  cout << "Rectangle:" << endl;
  out.str("");
  out << "de<" << dE_max << " && de>" << dE_min;
  out << " && mbc>" << Mbc_min << " && mbc<" << Mbc_max;
  Roo1DTable* recttable = ds.table(b0f,out.str().c_str());
  recttable->Print();
  recttable->Print("v");

  cout << "Ellips:" << endl;
//  out.str("");
//  out << "(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius.getVal() << "<1";
//  cout << out.str() << endl;
//  Roo1DTable* ellitable = ds.table(b0f,out.str().c_str());
  ellitable->Print();
  ellitable->Print("v");

  cout << "Elli:" << endl;
//  out.str("");
//  out << "(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius1.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius1.getVal() << "<1";
//  cout << out.str() << endl;
//  Roo1DTable* ellitable1 = ds.table(b0f,out.str().c_str());
  ellitable1->Print();
  ellitable1->Print("v");

  Roo1DTable* fulltable = ds.table(b0f);
  fulltable->Print();
  fulltable->Print("v");
  const int NSigTotal = fulltable->get("signal") + fulltable->get("fsr") + fulltable->get("bad_pi0");
  const double TruePur = ((double)NSIGNAL_ELLI)/(NSIGNAL_ELLI+ellitable1->get("comb")+ellitable1->get("rho2")+ellitable1->get("rho3")+ellitable1->get("rho4")+ellitable1->get("rho11"));

  cout << "Rectangle:" << endl;
  cout << "Nsig = " << nsig <<" +- " << nsig_err << " +- " << nsig_err_npq << " (" << nsig_err_total << ")" << endl;
  cout << "Npbg = " << nrho <<" +- " << nrho_err << " +- " << nrho_err_npq << " (" << nrho_err_total << ")" << endl;
  cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << " +- " << ncmb_err_npq << " (" << ncmb_err_total << ")" << endl;
  cout << "Pury = " << purity << " +- " << purity_err << endl;

  cout << "Ellips:" << endl;
  cout << "Nsig = " << nsigEl <<" +- " << nsig_errEl << " +- " << nsig_errEl_npq << " (" << nsig_errEl_total << ")" << endl;
  cout << "Npbg = " << nrhoEl <<" +- " << nrho_errEl << " +- " << nrho_errEl_npq << " (" << nrho_errEl_total << ")" << endl;
  cout << "Ncmb = " << ncmbEl <<" +- " << ncmb_errEl << " +- " << ncmb_errEl_npq << " (" << ncmb_errEl_total << ")" << endl;
  cout << "Pury = " << purityEl << " +- " << purity_errEl << endl;

  cout << "Elli:" << endl;
  cout << "Nsig = " << nsigEl1 <<" +- " << nsig_errEl1 << " +- " << nsig_errEl1_npq << " (" << nsig_errEl1_total << ")" << endl;
  cout << "Npbg = " << nrhoEl1 <<" +- " << nrho_errEl1 << " +- " << nrho_errEl1_npq << " (" << nrho_errEl1_total << ")" << endl;
  cout << "Ncmb = " << ncmbEl1 <<" +- " << ncmb_errEl1 << " +- " << ncmb_errEl1_npq << " (" << ncmb_errEl1_total << ")" << endl;
  cout << "Pury = " << purityEl1 << " +- " << purity_errEl1 << endl;

  cout << "Elli signal integral: " << intElli << endl;
  cout << "Nsig (full range): " << Nsig.getVal() << " +- " << Nsig.getError() << " (" << NSigTotal << ")" << endl;
  cout << "True purity: " << TruePur << endl;
}

