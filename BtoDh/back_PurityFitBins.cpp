#include "cuts.h"
#include "MyTools.h"

void PurityFitBins(const int _mode,const bool gen_pseudotoy = false){
  TChain* tree = new TChain("TEvent");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1.root");

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
  string label;
  switch(_mode){
  case 1:
    label = string("#pi^{0}");
    BDTG_MIN = bdtg_cut_pi0;
    mode.defineType("pi0",1);
    h0mode.defineType("gg",10);
    Mbc_min = mbc_min;
    Mbc_max = mbc_max;
    dE_min  = de_min;
    dE_max  = de_max;
    m_mode = 1;
    m_h0mode = 10;
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
    break;
  case 4:
    label = string("#omega");
    BDTG_MIN = bdtg_cut_omega;
    gg_flag = false;
    mode.defineType("omega",3);
    h0mode.defineType("ppp",20);
    Mbc_min = mbc_min;
    Mbc_max = mbc_max;
    dE_min  = de_min_omega;
    dE_max  = de_max_omega;
    m_mode = 3;
    m_h0mode = 20;
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
  const double mbcMax = 5.287;
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
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);
  RooRealVar tag_LH("tag_LH","tag_LH",-1.,1.); argset.add(tag_LH);
//  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  argset.add(b0f);
  RooDataSet* ds = new RooDataSet("ds","ds",tree,argset,"mbc>0||mbc<=0");

  stringstream out;
  out.str("");
  out << "de<" << dE_max << " && de>" << dE_min;
  out << " && mbc>" << Mbc_min << " && mbc<" << Mbc_max;
  Roo1DTable* sigtable = ds->table(b0f,out.str().c_str());
  sigtable->Print();
  sigtable->Print("v");

  Roo1DTable* fulltable = ds->table(b0f);
  fulltable->Print();
  fulltable->Print("v");

  ds->Print();
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

    RooRealVar deCBl("deCBl","deCBl",get_deCBl(m_mode,m_h0mode,_b0f),-0.2,0.1);      if(cSig) deCBl.setConstant(kTRUE);
    RooRealVar sCBl("sCBl","sCBl",get_sCBl(m_mode,m_h0mode,_b0f),0.,0.5);            if(cSig) sCBl.setConstant(kTRUE);
    RooRealVar alphal("alphal","alphal", get_alphal(m_mode,m_h0mode,_b0f), 0.,10.);  if(cSIG) alphal.setConstant(kTRUE);
    RooRealVar nl("nl","nl",2.,0.,100.); nl.setConstant(kTRUE);

    RooRealVar deCBr("deCBr","deCBr",get_deCBr(m_mode,m_h0mode,_b0f),-0.2,0.1);    if(cSig) deCBr.setConstant(kTRUE);
    RooRealVar sCBr("sCBr","sCBr",get_sCBr(m_mode,m_h0mode,_b0f),0.,0.5);          if(cSig) sCBr.setConstant(kTRUE);
    RooRealVar alphar("alphar","alphar",get_alphar(m_mode,m_h0mode,_b0f),-10.,0.); if(cSig) alphar.setConstant(kTRUE);
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
    RooRealVar c_s("c_s","c_s",get_c_s(_mode),0.0015,0.0035); if(cSig) c_s.setConstant(kTRUE);
    RooFormulaVar S("S","S","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s,b_s,a_s));

    RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);

    RooRealVar a_mbc0("a_mbc0","a_mbc0",get_a_mbc0(_mode)); if(cSig) a_mbc0.setConstant(kTRUE);
    RooRealVar b_mbc0("b_mbc0","b_mbc0",get_b_mbc0(_mode)); if(cSig) b_mbc0.setConstant(kTRUE);
    RooRealVar c_mbc0("c_mbc0","c_mbc0",get_c_mbc0(_mode),5.27,5.285); if(cSig) c_mbc0.setConstant(kTRUE);
    RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));
    RooNovosibirsk pdf_mbc_sig("pdf_mbc_sig","pdf_mbc_sig",mbc,MBC0,S,alpha);
  } else{
    RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
    RooRealVar c0("c0","c0",get_c0(_mode)); if(cSig) c0.setConstant(kTRUE);
    RooRealVar c1("c1","c1",get_c1(_mode)); if(cSig) c1.setConstant(kTRUE);
    RooRealVar c2("c2","c2",get_c2(_mode)); if(cSig) c2.setConstant(kTRUE);
    RooRealVar mbc0("mbc0","mbc0",5.284,5.28,5.29); if(cSig) mbc0.setConstant(kTRUE);
    RooFormulaVar MBC("MBC","MBC","@0+@1*TMath::Erf((@2-@3))/@4",RooArgList(mbc0,c0,c1,de,c2));

    RooRealVar a_s1("a_s1","a_s1",get_a_s(_mode),0.15,0.45); if(cSig) a_s1.setConstant(kTRUE);
    RooRealVar b_s1("b_s1","b_s1",get_b_s(_mode),-0.05,0.05); if(cSig) b_s1.setConstant(kTRUE);
    RooRealVar c_s1("c_s1","c_s1",get_c_s(_mode),0.0015,0.0035); if(cSig) c_s1.setConstant(kTRUE);
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

  RooRealVar mbc0_cmb_bb("mbc0_cmb_bb","mbc0_cmb_bb",get_mbc0_cmb_bb(_mode),5.25,5.29,"GeV"); if(cComb) mbc0_cmb_bb.setConstant(kTRUE);
  RooRealVar mbcWidth_cmb_bb("mbcWidth","mbcWidth",get_mbcw_cmb_bb(_mode),0.,0.1,"GeV"); if(cComb) mbcWidth_cmb_bb.setConstant(kTRUE);
  RooGaussian mbcGaus_cmb_bb("mbcGaus","mbcGaus",mbc,mbc0_cmb_bb,mbcWidth_cmb_bb);

  RooRealVar f_g("f_g","f_g",get_f_g_cmb_bb(_mode),0.,1.); if(cComb) f_g.setConstant(kTRUE);
  RooAddPdf pdf_mbc_cmb_bb("pdf_mbc_cmb_bb","pdf_mbc_cmb_bb",RooArgList(mbcGaus_cmb_bb,pdf_mbc_comb_ar),RooArgSet(f_g));

  RooRealVar argpar_cmb_qq("argpar_cmb_qq","argpar_cmb_qq",get_argpar_qq(_mode),-300,-10.); if(cComb) argpar_cmb_qq.setConstant(kTRUE);
  RooArgusBG pdf_mbc_cmb_qq("pdf_mbc_cmb_qq","pdf_mbc_cmb_qq",mbc,argedge,argpar_cmb_qq);
  
  /////////
  // pdf //
  /////////
  RooRealVar f_bb("f_bb","f_bb",0.2,0.,1.);
  RooProdPdf pdf_cmb_bb("pdf_cmb_bb","pdf_cmb_bb",pdf_mbc_cmb_bb,Conditional(pdf_de_comb_bb,de));
  RooProdPdf pdf_cmb_qq("pdf_cmb_qq","pdf_cmb_qq",pdf_mbc_cmb_qq,Conditional(pdf_de_comb_qq,de));
  RooAddPdf pdf_comb("pdf_comb","pdf_comb",RooArgSet(pdf_cmb_bb,pdf_cmb_qq),RooArgList(f_bb));

  /////////////////////
  // Peaking bkg PDF //
  /////////////////////
  /////////////////
  // peak de pdf //
  /////////////////
  RooRealVar de0r("de0r","de0r",get_de0r(_mode),-0.2,0.12);         if(cPeak) de0r.setConstant(kTRUE);
  RooRealVar slopel("slopel","slopel",get_slopel(_mode),-1.e5,0.);  if(cPeak) slopel.setConstant(kTRUE);
  RooRealVar sloper("sloper","sloper",get_sloper(_mode),-10000,0.); if(cPeak) sloper.setConstant(kTRUE);
  RooRealVar steep("steep","steep",get_steep(_mode),0.,1000.);      if(cPeak) steep.setConstant(kTRUE);
  RooRealVar p5("p5","p5",get_p5(_mode),0.01,1000.);                if(cPeak) p5.setConstant(kTRUE);
  RooRhoDeltaEPdf pdf_de_peak("pdf_de_peak","pdf_de_peak",de,de0r,slopel,sloper,steep,p5);
//  RooGenericPdf pdf_de_peak("pdf_de_peak","1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) > 0 ? 1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) : 0.001",RooArgSet(de,de0r,slopel,sloper,steep,p5));
  //////////////////
  // peak mbc pdf //
  //////////////////
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
    RooRealVar argpar_peak_bb("argpar_peak_bb","argpar_peak_bb",get_argpar_bb(_mode),-300,-10.); if(cPeak) argpar_peak_bb.setConstant(kTRUE);
    RooArgusBG pdf_mbc_peak_ar("pdf_mbc_peak_ar","Argus PDF",mbc,argedge,argpar_peak_bb);

    RooRealVar mbc0_peak("mbc0_peak","mbc0_peak",get_peak_b_mbc0(_mode),5.25,5.291,"GeV"); if(cPeak) mbc0_peak.setConstant(kTRUE);
    RooRealVar mbcWidth_peak("mbcWidth_peak","mbcWidth_peak",get_peak_b_s(_mode),0.,0.1,"GeV"); if(cPeak) mbcWidth_peak.setConstant(kTRUE);
    RooGaussian mbcGaus_peak("mbcGaus_peak","mbcGaus_peak",mbc,mbc0_peak,mbcWidth_peak);
    RooRealVar f_g_peak("f_g_peak","f_g_peak",get_f_g_cmb_bb(_mode),0.,1.); if(cPeak) f_g_peak.setConstant(kTRUE);

    RooAddPdf pdf_mbc_peak("pdf_mbc_peak","pdf_mbc_peak",RooArgList(mbcGaus_peak,pdf_mbc_peak_ar),RooArgSet(f_g_peak));
  }
  //////////////
  // peak pdf //
  //////////////
  RooProdPdf pdf_peak("pdf_peak","pdf_peak",pdf_de_peak,Conditional(pdf_mbc_peak,mbc));

  //////////////////
  // Complete PDF //
  //////////////////
  RooRealVar Nsig("Nsig","Nsig",1150,0.,10000.);
  RooRealVar Ncmb("Ncmb","Ncmb",2288,0,100000);
  double F_p_f_bbc = 0;
  switch(_mode){
  case 1:
    RooRealVar Npbg("Npbg","Npbg",100,0,100000.);
    break;
  case 2:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0051);
    F_p_f_bbc = f_p_f_bbc.getVal();
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  case 3:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0081);
    F_p_f_bbc = f_p_f_bbc.getVal();
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  case 4:
    RooConstVar f_p_f_bbc("f_p_f_bbc","f_p_f_bbc",0.0031);
    F_p_f_bbc = f_p_f_bbc.getVal();
    RooFormulaVar Npbg("Npbg","Npbg","@0*@1*@2",RooArgList(Ncmb,f_bb,f_p_f_bbc));
    break;
  default:
    return -1;
  }

  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_peak,pdf_comb),RooArgList(Nsig,Npbg,Ncmb));
  pdf.fitTo(*ds,Verbose(),Timer(true));

  if(_mode == 1) F_p_f_bbc = Npbg.getVal()/(Ncmb.getVal()*f_bb.getVal());

  de.setRange("Elli",dE_min,dE_max);
  RooFormulaVar mbclo1("mbclo1","@1-@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
  RooFormulaVar mbchi1("mbchi1","@1+@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
  mbc.setRange("Elli",mbclo1,mbchi1);

  int Flv,Bin;
  RooAbsReal* m_intSigEl1 = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
  const double m_sigint = m_intSigEl1->getVal();
  RooAbsReal* m_intRhoEl1 = pdf_peak.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
  const double m_peakint = m_intRhoEl1->getVal();
  RooAbsReal* m_intCmbEl1 = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
  const double m_combint = m_intCmbEl1->getVal();

  cout << "Signal region signal  fraction: " << m_sigint << endl;
  cout << "Signal region combin  fraction: " << m_combint << endl;
  cout << "Signal region peaking fraction: " << m_peakint << endl;

  double NSPredicted[2][16];
  double NSPredictedTot[2][16];
  double NEvTot[2][16];
  double NSPredictedErr[2][16];

  double NSPredictedTagV[2][16];
  double NSPredictedTotTagV[2][16];
  double NEvTotTagV[2][16];
  double NSPredictedErrTagV[2][16];

  double NS[2][16];
  double NC[2][16];
  double NP[2][16];
  double Pur[2][16];

  double NSTagV[2][16];
  double NCTagV[2][16];
  double NPTagV[2][16];
  double PurTagV[2][16];

  double NSPrTot = 0;
  double NSPrTotTagV = 0;
  double NSTot = 0;
  double NCTot = 0;
  double NPTot = 0;
  double NGoodEventsTot = 0;
  double NGoodEvents[2][16];
//  double NSFit[2][16];
//  double NSFitErr[2][16];
//  double NCFit[2][16];
//  double NCFitErr[2][16];
//  double NPFit[2][16];
//  double NPFitErr[2][16];

//  TTree* new_tree = ds->tree();
//  new_tree->Print();

  for(int i=0; i<2; i++){
    for(int j=0; j<16; j++){
      NGoodEvents[i][j] = 0;
      NSPredicted[i][j] = 0;
      NSPredictedTot[i][j] = 0;
      NEvTot[i][j] = 0;
      NSPredictedErr[i][j] = 0;
      NSPredictedTagV[i][j] = 0;
      NSPredictedTotTagV[i][j] = 0;
      NEvTotTagV[i][j] = 0;
      NSPredictedErrTagV[i][j] = 0;
      NS[i][j] = 0;
      NC[i][j] = 0;
      NP[i][j] = 0;
      Pur[i][j] = 0;
      NSTagV[i][j] = 0;
      NCTagV[i][j] = 0;
      NPTagV[i][j] = 0;
      PurTagV[i][j] = 0;
    }
  }

  const int NTot = tree->GetEntries();
  Int_t l_exp, l_flv, l_bin;
  Double_t l_tag;

  Int_t l_b0f, l_mode, l_h0mode, l_good_icpv;
  Double_t l_bdtg, l_mbc, l_de, l_mk, l_md;

  tree->SetBranchAddress("exp",&l_exp);
  tree->SetBranchAddress("flv_mc",&l_flv);
  tree->SetBranchAddress("bin_mc",&l_bin);

  tree->SetBranchAddress("b0f",&l_b0f);
  tree->SetBranchAddress("mode",&l_mode);
  tree->SetBranchAddress("h0mode",&l_h0mode);
//  tree->SetBranchAddress("good_icpv",&l_good_icpv);

  tree->SetBranchAddress("tag_LH",&l_tag);

  tree->SetBranchAddress("bdtg",&l_bdtg);
  tree->SetBranchAddress("mbc",&l_mbc);
  tree->SetBranchAddress("de",&l_de);
  tree->SetBranchAddress("mk",&l_mk);
  tree->SetBranchAddress("md",&l_md);

  cout << "Filling the arrays..." << endl;
  double average_wrong_tag[16] = 0;
  for(int i=0; i<16; i++){average_wrong_tag[i] = 0;}
  for(int i=0; i<NTot; i++){
    if(!(i%100000)) cout << i << endl;
    tree->GetEvent(i);
//    l_tag *= -1;
//    l_flv *= -1;
    if(l_mode != m_mode || l_h0mode != m_h0mode) continue;
    if(l_b0f<-1 || l_b0f == 0 || l_b0f>11) continue;
    if(l_bin == 0) continue;
    if(l_bdtg<BDTG_MIN) continue;
    if(abs(l_md-DMass-md_cut)>md_cut) continue;
    if(abs(l_mk-KMass-mk_cut)>md_cut) continue;
    const double ElR2 = EllipsR2(l_de,l_mbc,mbc_center.getVal(),de_center.getVal(),mbc_radius1.getVal(),de_radius1.getVal());
    if(ElR2>1) continue;

    NGoodEventsTot++;
    const int bin_bin = bin_ind(l_bin);
    const int flv_bin = flv_ind(l_flv);
    const int flv_tagv = l_tag>0 ? 1 : -1;
    const int flv_bin_tagv = flv_ind(flv_tagv);
    NGoodEvents[flv_bin][bin_bin]++;
//    cout << "flv_tagv: " << flv_tagv << ", flv: " << l_flv << endl;
    average_wrong_tag[bin_bin] += get_wtag_prob(l_tag,l_exp);
//    cout << l_tag << " " << l_exp << " " << get_wbin(l_tag) << " " << get_wtag_prob(l_tag,l_exp) << endl;
//    cout << get_wtag_prob(l_tag,l_exp) << " ";
//    cout << bin_bin << " " << flv_bin << " " << flv_tagv << " " << flv_bin_tagv << endl;

    NEvTot[flv_bin][bin_bin]++;
    NEvTotTagV[flv_bin_tagv][bin_bin]++;

    if(l_b0f == 1 || l_b0f == 5 || l_b0f == 10){
      NS[flv_bin][bin_bin]++;
      NSTagV[flv_bin_tagv][bin_bin]++;
      NSTot++;
    } else if(l_b0f == -1){
      NC[flv_bin][bin_bin]++;
      NCTagV[flv_bin_tagv][bin_bin]++;
      NCTot++;
    } else{
      NP[flv_bin][bin_bin]++;
      NPTagV[flv_bin_tagv][bin_bin]++;
      NPTot++;
    }
  }
  for(int i=0; i<16; i++){ average_wrong_tag[i] /= (NGoodEvents[0][i]+NGoodEvents[1][i]);}
  cout << "Average wrong tag probability: " << endl;
  for(int i=0; i<8; i++){
    cout << (8-i) << " bin :" << average_wrong_tag[i] << ", " << -(8-i) << " bin :" << average_wrong_tag[15-i] << endl;
  }

  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      const int l_bin = bin(j);
      const int l_flv = flv(k);
      double frac      = 0.5*N(l_bin,l_flv)*Nsig.getVal()*m_sigint;
      double frac_tagv = 0.5*N(l_bin,l_flv,average_wrong_tag[j])*Nsig.getVal()*m_sigint;
      NSPredicted[k][j] += frac;
      NSPredictedTot[k][j] += frac;
      NSPrTot += frac;
      NSPredictedTagV[k][j] += frac_tagv;
      NSPredictedTotTagV[k][j] += frac_tagv;
      NSPrTotTagV += frac_tagv;

      Pur[k][j] = NS[k][j] / NEvTot[k][j];
      PurTagV[k][j] = NSTagV[k][j] / NEvTotTagV[k][j];

      frac /= Nsig.getVal();
      const double sigma1 = frac*Nsig.getError();
      const double snpq2  = Nsig.getVal() * frac * (1 - frac);
      NSPredictedErr[k][j] = sqrt((sigma1*sigma1 + snpq2));

      frac_tagv /= Nsig.getVal();
      const double sigma1_tagv = frac_tagv * Nsig.getError();
      const double snpq2_tagv  = Nsig.getVal() * frac_tagv * (1 - frac_tagv);
      NSPredictedErrTagV[k][j] = sqrt((sigma1_tagv*sigma1_tagv + snpq2_tagv));
    }
  }

//  RooRealVar* nsig_vec[2][16];
//  RooRealVar* ncomb_vec[2][16];
//  RooRealVar* ncomb_bb_vec[2][16];
//  RooRealVar* npeak_vec[2][16];
//  RooAddPdf*  pdf_vec[2][16];
//  RooSimultaneous PDF("PDF","PDF",binflv);

//  f_bb.setConstant(kTRUE);
//  argedge.setConstant(kTRUE);
//  for(int i=0; i<2; i++){
//    Flv = flv(i);
//    for(int j=0; j<16; j++){
//      Bin = bin(j); Flv = flv(i);
//      out.str("");
//      out << "bin == " << Bin << " && flv_mc == " << Flv;
//      const double m_nsig = 0.5*Nsig.getVal()*N(Bin,Flv);

//      out.str("");
//      out << "(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius1.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius1.getVal() << "<1";
//      Roo1DTable* m_ellitable1 = ds->table(RooArgSet(b0f,flv,bin),out.str().c_str());
//      Roo1DTable* m_totaltable = ds->table(RooArgSet(b0f,flv,bin));

//      string bf_str; out.str("");
//      if(Flv == 1) out << string(";B0;");
//      else         out << string(";anti-B0;");
//      out << Bin << "}";
//      bf_str = string("{signal")  + out.str();
//      NS[i][j]  = m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] = m_totaltable->get(bf_str.c_str());
//      bf_str = string("{fsr")     + out.str();
//      NS[i][j] += m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      bf_str = string("{bad_pi0") + out.str();
//      NS[i][j] += m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      NSTot    += NS[i][j];
//      bf_str = string("{comb")     + out.str();
//      NC[i][j]  = m_ellitable1->get(bf_str.c_str());//*m_combint;
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      bf_str = string("{rho2")     + out.str();
//      NP[i][j]  = m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      bf_str = string("{rho3")     + out.str();
//      NP[i][j] += m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      bf_str = string("{rho4")     + out.str();
//      NP[i][j] += m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      bf_str = string("{rho11")    + out.str();
//      NP[i][j] += m_ellitable1->get(bf_str.c_str());
//      NEvTot[i][j] += m_totaltable->get(bf_str.c_str());
//      Pur[i][j] = NS[i][j]/(NS[i][j] + NC[i][j] + NP[i][j])*100;

//      const double frac = 0.5*N(Bin,Flv)*m_sigint;
//      const double sigma1   = frac * Nsig.getError();
//      const double snpq2 = Nsig.getVal() * frac * (1 - frac);
//      NSPredicted[i][j]    = frac * Nsig.getVal(); NSPrTot += NSPredicted[i][j];
//      NSPredictedTot[i][j] = Nsig.getVal()*0.5*N(Bin,Flv);
//      NSPredictedErr[i][j] = sqrt((sigma1*sigma1 + snpq2));

//      cout << "Define fractions";
//      out.str("");
//      out << "Nsig" << i << "_" << j;
//      nsig_vec[i][j] = new RooRealVar(out.str().c_str(),out.str().c_str(),m_nsig,0.,Nsig.getVal());
//      out.str("");
//      out << "Ncmb" << i << "_" << j;
//      ncomb_vec[i][j] = new RooRealVar(out.str().c_str(),out.str().c_str(),Ncmb.getVal()/32*(1.-f_bb.getVal()),0.,Ncmb.getVal());
//      out.str("");
//      out << "Ncmb_bb" << i << "_" << j;
//      ncomb_bb_vec[i][j] = new RooRealVar(out.str().c_str(),out.str().c_str(),Ncmb.getVal()/32*f_bb.getVal(),0.,Ncmb.getVal());
//      out.str("");
//      out << "Npbg" << i << "_" << j;
//      npeak_vec[i][j] = new RooRealVar(out.str().c_str(),out.str().c_str(),Npbg.getVal()/32,0.,Npbg.getVal());
//      if(_mode != 1) npeak_vec[i][j]->setConstant(kTRUE);
//      cout << ", Define pdf";
//      out.str("");
//      out << "pdf" << i << "_" << j;
//      pdf_vec[i][j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(pdf_sig,pdf_peak,pdf_de_comb_bb,pdf_de_comb_qq),RooArgList(*nsig_vec[i][j],*npeak_vec[i][j],*ncomb_bb_vec[i][j],*ncomb_vec[i][j]));

//      cout << ", Add pdf" << endl;
//      out.str("");
//      if(i == 0) out << "{" << Bin << ";anti-B0}";
//      else       out << "{" << Bin << ";B0}";
//      PDF.addPdf(*pdf_vec[i][j],out.str().c_str());
//    }
//  }

//  PDF.fitTo(*ds,Verbose(),Timer(true));

  double bins[32];
  double bins_err[32];
  double sig_true[32];
  double sig_true_tagv[32];
//  double sig_fit[32], sig_fit_err[32], sig_fit_diff[32];
  double sig_pre[32], sig_pre_err[32], sig_pre_diff[32];
  double sig_pre_tagv[32], sig_pre_err_tagv[32], sig_pre_diff_tagv[32];
//  double chi2_fit = 0;
  double chi2_pre = 0;
  double chi2_pre_tagv = 0;

  double NPPredicted[2][16];
  double NCPredicted[2][16];
  double NBPredicted[2][16];

  double NPPredictedTagV[2][16];
  double NCPredictedTagV[2][16];
  double NBPredictedTagV[2][16];

//  double frac_pre[2][16];
//  double frac_pre_tagv[2][16];

  for(int i=0; i<2; i++){
    for(int j=0; j<16; j++){
      NBPredicted[i][j] = NEvTot[i][j] - NSPredicted[i][j];
      if(NBPredicted[i][j]<0){
        NSPredicted[i][j] = NEvTot[i][j];
        NBPredicted[i][j] = 0;
      }

      NBPredictedTagV[i][j] = NEvTotTagV[i][j] - NSPredictedTagV[i][j];
      if(NBPredictedTagV[i][j]<0){
        NSPredictedTagV[i][j] = NEvTotTagV[i][j];
        NBPredictedTagV[i][j] = 0;
      }
      const double denominator = m_combint+F_p_f_bbc*f_bb.getVal()*m_peakint;
      const double comb_coeff = (F_p_f_bbc*f_bb.getVal()*m_peakint)/denominator;
      const double peak_coeff = m_combint/denominator;
      NPPredicted[i][j] = NBPredicted[i][j] * comb_coeff;
      NCPredicted[i][j] = NBPredicted[i][j] * peak_coeff;

      NPPredictedTagV[i][j] = NBPredictedTagV[i][j] * comb_coeff;
      NCPredictedTagV[i][j] = NBPredictedTagV[i][j] * peak_coeff;

//      NSFit[i][j] = nsig_vec[i][j]->getVal()*m_sigint;   NSFitErr[i][j] = nsig_vec[i][j]->getError()*m_sigint;
//      NCFit[i][j] = ncomb_vec[i][j]->getVal()*m_combint; NCFitErr[i][j] = ncomb_vec[i][j]->getError()*m_combint;
//      NPFit[i][j] = npeak_vec[i][j]->getVal()*m_peakint; NPFitErr[i][j] = npeak_vec[i][j]->getError()*m_peakint;
//      const double snpqsq_sig = nsig_vec[i][j]->getVal()*m_sigint*(1-m_sigint);
//      NSFitErr[i][j] = sqrt(NSFitErr[i][j]*NSFitErr[i][j]+snpqsq_sig);

      cout << "Flv: " << flv(i) << ", Bin: " << bin(j) << endl;
//      cout << " Perfect tagging:" << endl;
//      cout << ", NS(model): " << NSPredicted[i][j] << " +- " << NSPredictedErr[i][j] << ", NC: " << NCPredicted[i][j] << ", NP: " << NPPredicted[i][j] << endl;
//      cout << "  NS       : " << NS[i][j] << ", NC: " << NC[i][j] << ", NP: " << NP[i][j] << ", Pur: " << Pur[i][j] << endl;
//      cout << " TagV:" << endl;
      cout << ", NS(model): " << NSPredictedTagV[i][j] << " +- " << NSPredictedErrTagV[i][j] << ", NC: " << NCPredictedTagV[i][j] << ", NP: " << NPPredictedTagV[i][j] << endl;
      cout << "  NS       : " << NSTagV[i][j] << ", NC: " << NCTagV[i][j] << ", NP: " << NPTagV[i][j] << ", Pur: " << PurTagV[i][j] << endl;
//      cout << "  NS(fit)  : " << NSFit[i][j] << " +- " << NSFitErr[i][j] << ", NCFit = " << NCFit[i][j] << " +- " << NCFitErr[i][j];
//      cout << ", NPFit = " << NPFit[i][j] << " +- " << NPFitErr[i][j] << endl << endl;

      int k = i*16 + j;
      bins[k] = k;
      bins_err[k] = 0;
      sig_true[k] = NS[i][j];
//      sig_fit[k] = NSFit[i][j]; sig_fit_err[k] = NSFitErr[i][j];
      sig_pre[k] = NSPredicted[i][j]; sig_pre_err[k] = NSPredictedErr[i][j];
//      sig_fit_diff[k] = sig_fit[k] - sig_true[k];
      sig_pre_diff[k] = sig_pre[k] - sig_true[k];
//      chi2_fit += sig_fit_diff[k]*sig_fit_diff[k]/(sig_fit_err[i]*sig_fit_err[i]);
      chi2_pre += sig_pre_diff[k]*sig_pre_diff[k]/(sig_pre_err[i]*sig_pre_err[i]);

      sig_true_tagv[k] = NSTagV[i][j];
//      sig_fit_tagv[k] = NSFit[i][j]; sig_fit_err[k] = NSFitErr[i][j];
      sig_pre_tagv[k] = NSPredictedTagV[i][j]; sig_pre_err_tagv[k] = NSPredictedErrTagV[i][j];
//      sig_fit_diff_tagv[k] = sig_fit[k] - sig_true[k];
      sig_pre_diff_tagv[k] = sig_pre_tagv[k] - sig_true_tagv[k];
//      chi2_fit += sig_fit_diff[k]*sig_fit_diff[k]/(sig_fit_err[i]*sig_fit_err[i]);
      chi2_pre_tagv += sig_pre_diff_tagv[k]*sig_pre_diff_tagv[k]/(sig_pre_err_tagv[i]*sig_pre_err_tagv[i]);
    }
  }
//  chi2_fit /= 32.;
  chi2_pre /= 32.;
  chi2_pre_tagv /= 32.;

//  cout << "NSPrTot     = " << NSPrTot << " (" << Nsig.getVal()*m_sigint << ")" << endl;
  cout << "NSPrTotTagV = " << NSPrTotTagV << " (" << Nsig.getVal()*m_sigint << ")" << endl;
  cout << "NSTot: " << NSTot << ", NCTot: " << NCTot << ", NPTot: " << NPTot << endl;

  TMultiGraph* mg_sig = new TMultiGraph("mg_sig","N_{sig} vs. Dalitz bin & Flavor");
  TMultiGraph* mg_sig_diff = new TMultiGraph("mg_sig_diff","N_{sig}-N_{sig}^{true} vs. Dalitz bin & Flavor");
  TGraph* gr_sig_true = new TGraph(32,bins,sig_true);
  gr_sig_true->SetMarkerStyle(20);
  gr_sig_true->SetMarkerColor(kBlue);
  gr_sig_true->SetMarkerSize(1.);

  TGraph* gr_sig_true_tagv = new TGraph(32,bins,sig_true_tagv);
  gr_sig_true_tagv->SetMarkerStyle(20);
  gr_sig_true_tagv->SetMarkerColor(kBlue);
  gr_sig_true_tagv->SetMarkerSize(1.);
//  TGraphErrors* gr_sig_fit = new TGraphErrors(32,bins,sig_fit,bins_err,sig_fit_err);
//  gr_sig_fit->SetMarkerStyle(20);
//  gr_sig_fit->SetMarkerColor(kRed);
//  gr_sig_fit->SetMarkerSize(1.);
  TGraphErrors* gr_sig_pre = new TGraphErrors(32,bins,sig_pre,bins_err,sig_pre_err);
  gr_sig_pre->SetMarkerStyle(21);
  gr_sig_pre->SetMarkerColor(kBlue);
  gr_sig_pre->SetMarkerSize(1.);

  TGraphErrors* gr_sig_pre_tagv = new TGraphErrors(32,bins,sig_pre_tagv,bins_err,sig_pre_err_tagv);
  gr_sig_pre_tagv->SetMarkerStyle(20);
  gr_sig_pre_tagv->SetMarkerColor(kBlack);
  gr_sig_pre_tagv->SetMarkerSize(1.);
//  TGraphErrors* gr_sig_fit_diff = new TGraphErrors(32,bins,sig_fit_diff,bins_err,sig_fit_err);
//  gr_sig_fit_diff->SetMarkerStyle(20);
//  gr_sig_fit_diff->SetMarkerColor(kRed);
//  gr_sig_fit_diff->SetMarkerSize(1.);
  TGraphErrors* gr_sig_pre_diff = new TGraphErrors(32,bins,sig_pre_diff,bins_err,sig_pre_err);
  gr_sig_pre_diff->SetMarkerStyle(20);
  gr_sig_pre_diff->SetMarkerColor(kBlue);
  gr_sig_pre_diff->SetMarkerSize(1.);

  TGraphErrors* gr_sig_pre_diff_tagv = new TGraphErrors(32,bins,sig_pre_diff_tagv,bins_err,sig_pre_err_tagv);
  gr_sig_pre_diff_tagv->SetMarkerStyle(20);
  gr_sig_pre_diff_tagv->SetMarkerColor(kBlack);
  gr_sig_pre_diff_tagv->SetMarkerSize(1.);

  TCanvas* c_sig = new TCanvas("c_sig","c_sig",800,800);
  c_sig->Divide(1,2);
  c_sig->Draw();
  c_sig->cd(1);
  c_sig->Pad()->SetGrid();
//  mg_sig->Add(gr_sig_true);
  mg_sig->Add(gr_sig_true_tagv);
//  mg_sig->Add(gr_sig_fit);
//  mg_sig->Add(gr_sig_pre);
  mg_sig->Add(gr_sig_pre_tagv);
  mg_sig->Draw("ap");

  c_sig->cd(2);
  c_sig->Pad()->SetGrid();
//  mg_sig_diff->Add(gr_sig_fit_diff);
//  mg_sig_diff->Add(gr_sig_pre_diff);
  mg_sig_diff->Add(gr_sig_pre_diff_tagv);
  mg_sig_diff->Draw("ap");
  pt = new TPaveText(0.65,0.80,0.98,0.98,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/nbins = " << chi2_pre_tagv;
  pt->AddText(out.str().c_str());
  pt->Draw();
  c_sig->Update();

//  cout << "Chi2/nbins(fit) = " << chi2_fit << endl;
  cout << "Chi2/nbins(pre)      = " << chi2_pre << endl;
  cout << "Chi2/nbins(pre) TagV = " << chi2_pre_tagv << endl;

//  return;
  if(!gen_pseudotoy) return;

  delete tree;
  delete ds;

  // * Background dt pdf * //
  RooRealVar dt("dt","dt",-70,70);
  RooRealVar f_bkg_delta("f_delta","f_delta",get_f_bkg_delta(_mode,55,1));
  RooRealVar mu_bkg_delta("mu_bkg_delta","mu_bkg_delta",get_mu_bkg_delta(_mode,55));
  RooRealVar mu_bkg("mu_bkg","mu_bkg",get_mu_bkg(_mode,55));
  RooRealVar sigma_mn_bkg("sigma_mn_bkg","sigma_mn_bkg",get_sigma_mn_bkg(_mode,55,1));
  RooRealVar sigma_tl_bkg("sigma_tl_bkg","sigma_tl_bkg",get_sigma_tl_bkg(_mode,55,1));
  RooRealVar tau_bkg("tau_bkg","tau_bkg",get_tau_bkg(_mode,55,1));

  RooFormulaVar dt_delta("dt_delta","dt_delta","@0-@1",RooArgSet(dt,mu_bkg_delta));
  RooGaussModel gauss_mn_model("gauss_mn_model","gauss_mn_model",dt,mu_bkg,sigma_mn_bkg);
  RooGaussModel gauss_tl_model("gauss_tl_model","gauss_tl_model",dt,mu_bkg,sigma_tl_bkg);
  RooRealVar f_bkg_tail("f_bkg_tail","f_bkg_tail",get_f_tl_bkg(_mode,55,1));
  RooAddModel gauss_model("gauss_model","gauss_model",RooArgSet(gauss_tl_model,gauss_mn_model),RooArgList(f_bkg_tail));

  RooTruthModel delta_model("delta_model","delta_model",dt);
  RooDecay dec_delta_pdf("dec_delta_pdf","dec_delta_pdf",dt,tau_bkg,delta_model,RooDecay::DoubleSided);
  RooDecay dec_gauss_pdf("dec_gauss_pdf","dec_gauss_pdf",dt,tau_bkg,gauss_model,RooDecay::DoubleSided);
//  RooDecay bkg_dt_pdf("bkg_dt_pdf","bkg_dt_pdf",dt,tau_bkg,delta_model,RooDecay::DoubleSided);
  RooAddPdf bkg_dt_pdf("bkg_dt_pdf","bkg_dt_pdf",RooArgSet(dec_delta_pdf,dec_gauss_pdf),RooArgList(f_bkg_delta));

  TChain* tree_sig = new TChain("TEvent");
  switch (_mode) {
  case 1:
    tree_sig->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
    break;
  case 2:
    tree_sig->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
    break;
  case 3:
    tree_sig->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
    break;
  case 4:
    tree_sig->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s5_full.root");
    break;
  default:
      break;
  }

  const int NSigTot = tree_sig->GetEntries();
  out.str("");
  out << "toy_tree_m" << _mode << "_v2" << ".root";
  TFile* toy_file = new TFile(out.str().c_str(),"RECREATE");
  TTree* toy_tree = new TTree("TEvent","TEvent");

  Int_t m_exp, m_bin, m_flv, m_good_icpv, m_ntrk_rec, m_ntrk_asc, m_ndf_rec, m_ndf_asc, m_b0f;
  Double_t m_chi2_vtx_d0, m_costhBcms, m_tag, m_z_sig, m_z_asc, m_sz_rec, m_sz_asc, m_chisq_rec, m_chisq_asc;
  Double_t m_de, m_mbc;
  Double_t m_f_bkg;
  Int_t m_elli;

  tree_sig->SetBranchAddress("exp",&m_exp);
  tree_sig->SetBranchAddress("chi2_vtx_d0",&m_chi2_vtx_d0);
  tree_sig->SetBranchAddress("costhBcms", &m_costhBcms);
  tree_sig->SetBranchAddress("bin_mc",&m_bin);
  tree_sig->SetBranchAddress("flv_mc",&m_flv);
  tree_sig->SetBranchAddress("tag_LH",&m_tag);
  tree_sig->SetBranchAddress("ntrk_sig",&m_ntrk_rec);
  tree_sig->SetBranchAddress("ntrk_asc",&m_ntrk_asc);
  tree_sig->SetBranchAddress("ndf_z_sig",&m_ndf_rec);
  tree_sig->SetBranchAddress("ndf_z_asc",&m_ndf_asc);
  tree_sig->SetBranchAddress("z_sig",&m_z_sig);
  tree_sig->SetBranchAddress("sz_sig",&m_sz_rec);
  tree_sig->SetBranchAddress("chisq_z_sig",&m_chisq_rec);
  tree_sig->SetBranchAddress("good_icpv",&m_good_icpv);
  tree_sig->SetBranchAddress("z_asc",&m_z_asc);
  tree_sig->SetBranchAddress("sz_asc",&m_sz_asc);
  tree_sig->SetBranchAddress("b0f",&m_b0f);
  tree_sig->SetBranchAddress("chisq_z_asc",&m_chisq_asc);
  tree_sig->SetBranchAddress("de",&m_de);
  tree_sig->SetBranchAddress("mbc",&m_mbc);

  toy_tree->Branch("exp",&m_exp,"exp/I");
  toy_tree->Branch("chi2_vtx_d0",&m_chi2_vtx_d0,"chi2_vtx_d0/D");
  toy_tree->Branch("costhBcms",&m_costhBcms,"costhBcms/D");
  toy_tree->Branch("bin_mc",&m_bin,"bin_mc/I");
  toy_tree->Branch("flv_mc",&m_flv,"flv_mc/I");
  toy_tree->Branch("de",&m_de,"de/D");
  toy_tree->Branch("mbc",&m_mbc,"mbc/D");
  toy_tree->Branch("good_icpv",&m_good_icpv,"good_icpv/I");
  toy_tree->Branch("tag_LH",&m_tag,"tag_LH/D");
  toy_tree->Branch("ntrk_rec",&m_ntrk_rec,"ntrk_rec/I");
  toy_tree->Branch("ntrk_asc",&m_ntrk_asc,"ntrk_asc/I");
  toy_tree->Branch("ndf_z_sig",&m_ndf_rec,"ndf_z_sig/I");
  toy_tree->Branch("ndf_asc",&m_ndf_asc,"ndf_asc/I");
  toy_tree->Branch("z_sig",&m_z_sig,"z_sig/D");
  toy_tree->Branch("z_asc",&m_z_asc,"z_asc/D");
  toy_tree->Branch("sz_sig",&m_sz_rec,"sz_sig/D");
  toy_tree->Branch("sz_asc",&m_sz_asc,"sz_asc/D");
  toy_tree->Branch("b0f",&m_b0f,"b0f/I");
  toy_tree->Branch("chisq_z_asc",&m_chisq_asc,"chisq_z_asc/D");
  toy_tree->Branch("chisq_z_sig",&m_chisq_rec,"chisq_z_sig/D");

  toy_tree->Branch("f_bkg",&m_f_bkg,"f_bkg/D");
  toy_tree->Branch("elli",&m_elli,"elli/I");

  const double cm2ps = 78.48566945838871754705;
  TRandom3* rndm = new TRandom3();
  rndm->SetSeed(0);
  cout << "Total events: " << NSigTot << endl;
  for(int i=0; i<NSigTot; i++){
    if(!(i%1000)) cout << "Event " << i << endl;
    tree_sig->GetEvent(i);
    if(!m_good_icpv) continue;
    if(m_b0f != 1 && m_b0f != 5 && m_b0f != 10) continue;
//    if(!IsInEllips(m_de,m_mbc,mbc_center.getVal(),de_center.getVal(),mbc_radius1.getVal(),de_radius1.getVal())) m_elli = 0;
//    else m_elli = 1;
    de.setVal(m_de); mbc.setVal(m_mbc);
    const int Bin_ind = bin_ind(m_bin);
    const int Flv_ind = flv_ind(m_flv);
    const double sig_val  = NSPredictedTagV[Flv_ind][Bin_ind]*pdf_sig.getVal(RooArgSet(de,mbc))/m_sigint;
    const double cmb_val  = NCPredictedTagV[Flv_ind][Bin_ind]*pdf_comb.getVal(RooArgSet(de,mbc))/m_combint;
    const double peak_val = NPPredictedTagV[Flv_ind][Bin_ind]*pdf_peak.getVal(RooArgSet(de,mbc))/m_peakint;
    m_f_bkg = (cmb_val + peak_val)/(cmb_val + peak_val + sig_val);
    if(m_f_bkg < 0 || m_f_bkg > 1){
      cout << m_f_bkg << " " << cmb_val << " " << peak_val << " " << sig_val << " " << m_de << " " << m_mbc << " " << m_bin << " " << m_flv << endl;
    }
    const double ElR2 = EllipsR2(m_de,m_mbc,mbc_center.getVal(),de_center.getVal(),mbc_radius1.getVal(),de_radius1.getVal());
    if(m_f_bkg>0.6) continue;
    if(ElR2 > 1) m_elli = 0;
    else         m_elli = 1;

    const bool sig_event = m_f_bkg < rndm->Rndm() ? true : false;

    RooDataSet* dtbkgset;
    RooArgSet* dtset;
    RooRealVar* DT;
    if(!sig_event){
      m_b0f = -1;
      const double SIGMA = 0.1*sqrt(m_sz_rec*m_sz_rec+m_sz_asc*m_sz_asc)*cm2ps;
//      cout << "Sigma: " << SIGMA << endl;
      f_bkg_delta.setVal(get_f_bkg_delta(_mode,m_exp,m_ndf_asc));
      mu_bkg_delta.setVal(get_mu_bkg_delta(_mode,m_exp));
      mu_bkg.setVal(get_mu_bkg(_mode,m_exp));
      sigma_mn_bkg.setVal(get_sigma_mn_bkg(_mode,m_exp,m_ndf_asc)*SIGMA);
      sigma_tl_bkg.setVal(get_sigma_tl_bkg(_mode,m_exp,m_ndf_asc)*SIGMA);
      tau_bkg.setVal(get_tau_bkg(_mode,m_exp,m_ndf_asc));
      f_bkg_tail.setVal(get_f_tl_bkg(_mode,m_exp,m_ndf_asc));
//      cout << f_bkg_delta.getVal() << " " << mu_bkg_delta.getVal() << " " << mu_bkg.getVal() << " " << sigma_mn_bkg.getVal() << " " << sigma_tl_bkg.getVal() << " " << tau_bkg.getVal() << " " << f_bkg_tail.getVal() << endl;
      dtbkgset = (RooDataSet*)bkg_dt_pdf.generate(dt,1);
      dtset    = (RooArgSet*)dtbkgset->get(0);
      DT       = (RooRealVar*)dtset->find("dt");
      m_z_sig = DT->getVal()*10./cm2ps;
      m_z_asc = 0;
      dtbkgset->Clear();
      dtset->Clear();
      DT->Clear();
    }
    toy_tree->Fill();
  }
  toy_tree->Write();
  toy_file->Close();
}
