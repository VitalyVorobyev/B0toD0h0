#include "cuts.h"
using namespace RooFit;

void Purity_2d_fit_etagg(bool data = false){
  TChain* tree = new TChain("TEvent");
  if(!data) tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1_full.root");
  else      tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_data.root");

  const int _mode = 2;
  const int _h0mode = 10;
  const int _b0f = -1;
//  gROOT->ProcessLine(".L pdfs/RooRhoDeltaEPdf.cxx+");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  b0f.defineType("fsr",10);
  b0f.defineType("bad_pi0",5);
//  b0f.defineType("peak1",3);
//  b0f.defineType("peak2",4);
//  b0f.defineType("peak3",11);
//  b0f.defineType("peak4",20);
  b0f.defineType("comb",-1);

  RooArgSet argset;

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  double deMin = -0.15;
  if(keysflag) deMin = -0.3;
  const double deMax = 0.3;
  const double elliscaleDe  = TMath::Sqrt(4./TMath::Pi());
  const double elliscaleMbc = TMath::Sqrt(4./TMath::Pi());

  RooRealVar mbc_center("mbc_center","mbc_center",0.5*(mbc_min+mbc_max),mbc_min,mbc_max); mbc_center.setConstant(kTRUE);
  RooRealVar mbc_center_eq("mbc_center_eq","mbc_center_eq",mr_argedge_3-0.5*(mbc_max-mbc_min)*elliscaleMbc,mbc_min,mbc_max); mbc_center_eq.setConstant(kTRUE);
  RooRealVar de_center("de_center","de_center",0.5*(de_min+de_max),de_min,de_max); de_center.setConstant(kTRUE);
  RooRealVar mbc_radius("mbc_radius","mbc_radius",0.5*(mbc_max-mbc_min)*elliscaleMbc,0,0.5*(mbcMax-mbcMin)); mbc_radius.setConstant(kTRUE);
  RooRealVar de_radius("de_radius","de_radius",0.5*(de_max-de_min)*elliscaleDe,0.,0.5*(deMax-deMin)); de_radius.setConstant(kTRUE);
  RooRealVar mbc_radius1("mbc_radius1","mbc_radius1",0.5*(mbc_max-mbc_min),0,0.5*(mbcMax-mbcMin)); mbc_radius1.setConstant(kTRUE);
  RooRealVar de_radius1("de_radius1","de_radius1",0.5*(de_max-de_min),0.,0.5*(deMax-deMin)); de_radius1.setConstant(kTRUE);

  cout << 0.5*(mbc_min+mbc_max) << " " << 0.5*(mbc_max-mbc_min) << endl;
  cout << 0.5*(de_min+de_max) << " " << 0.5*(de_max-de_min) << endl;

  mbc_center.Print();
  mbc_center_eq.Print();

  const double BDTG_MIN = bdtg_cut_etagg;
  const double BDTG_MAX = 1;
  RooCategory mode("mode","mode");
  mode.defineType("eta",2);
  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  argset.add(mode);
  argset.add(h0mode);

  RooRealVar mbc("mbc","M_{bc}",0.5*(mbc_min+mbc_max),mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);

  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_min,de_max);
  
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mh0("mh0","mh0",EtaMass-meta_cut,EtaMass+meta_cut,"GeV"); argset.add(mh0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,1.); argset.add(bdtg);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  if(!data) argset.add(b0f);
  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");
//  RooDataSet* ds0 = ds.reduce(RooArgSet(de,mbc));

  stringstream out;
  if(!data){
    out.str("");
    out << "de<" << de_max << " && de>" << de_min;
    out << " && mbc>" << mbc_min << " && mbc<" << mbc_max;
    Roo1DTable* sigtable = ds.table(b0f,out.str().c_str());
    sigtable->Print();
    sigtable->Print("v");

    Roo1DTable* fulltable = ds.table(b0f);
    fulltable->Print();
    fulltable->Print("v");
  }

//  RooDataHist* dh = ds0->binnedClone();

  ds.Print();

  ////////////////
  // Signal PDF //
  ////////////////
  ////////////
  // de pdf //
  ////////////
  if(!(_mode!=1 && _h0mode == 20 && _b0f == -1)){
  RooRealVar de0("de0","de0",get_de0(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG) de0.setConstant(kTRUE);
  RooRealVar s1("s1","s1",get_s1(_mode,_h0mode,_b0f),0.,0.5);       if(cSIG) s1.setConstant(kTRUE);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar deCBl("deCBl","deCBl",get_deCBl(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG) deCBl.setConstant(kTRUE);
  RooRealVar sCBl("sCBl","sCBl",get_sCBl(_mode,_h0mode,_b0f),0.,0.5);       if(cSIG) sCBl.setConstant(kTRUE);
  RooRealVar nl("nl","nl",get_nl(_mode,_h0mode,_b0f),0.,100.);              if(cSIG) nl.setConstant(kTRUE);
//  RooRealVar nl("nl","nl",1.,0.,100.); nl.setConstant(kTRUE);
  RooRealVar alphal("alphal","alphal",get_alphal(_mode,_h0mode,_b0f),-10.,10.); if(cSIG) alphal.setConstant(kTRUE);

  RooRealVar deCBr("deCBr","deCBr",get_deCBr(_mode,_h0mode,_b0f),-0.2,0.1);     if(cSIG) deCBr.setConstant(kTRUE);
  RooRealVar sCBr("sCBr","sCBr",get_sCBr(_mode,_h0mode,_b0f),0.,0.5);           if(cSIG) sCBr.setConstant(kTRUE);
  RooRealVar nr("nr","nr",get_nr(_mode,_h0mode,_b0f),0.,100.);                  if(cSIG) nr.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",get_alphar(_mode,_h0mode,_b0f),-10.,10.); if(cSIG) alphar.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",get_fCBl(_mode,_h0mode,_b0f),0.,1.); if(cSIG) fCBl.setConstant(kTRUE);
  RooRealVar fCBr("fCBr","fCBr",get_fCBr(_mode,_h0mode,_b0f),0.,1.); if(cSIG) fCBr.setConstant(kTRUE);

  RooAddPdf pdf_de_sig("pdf_de_sig","pdf_de_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  }

  if((_mode==3&&_h0mode == 20 || _mode==2) && _b0f == -1){
  RooRealVar de0_201("de0_201","de0_201",get_de0(_mode,_h0mode,1),-0.1,0.1); de0_201.setConstant(kTRUE);
  RooRealVar s1_201("s1_201","s1_201",get_s1(_mode,_h0mode,1),0.,0.5);       s1_201.setConstant(kTRUE);
  RooGaussian g1_201("g1_201","g1_201",de,de0_201,s1_201);

  RooRealVar deCBl_201("deCBl_201","deCBl_201",get_deCBl(_mode,_h0mode,1),-0.1,0.1);     deCBl_201.setConstant(kTRUE);
  RooRealVar sCBl_201("sCBl_201","sCBl_201",get_sCBl(_mode,_h0mode,1),0.,0.5);           sCBl_201.setConstant(kTRUE);
  RooRealVar nl_201("nl_201","nl_201",get_nl(_mode,_h0mode,1),0.,100.);                  nl_201.setConstant(kTRUE);
  RooRealVar alphal_201("alphal_201","alphal_201",get_alphal(_mode,_h0mode,1),-10.,10.); alphal_201.setConstant(kTRUE);
  RooRealVar deCBr_201("deCBr_201","deCBr_201",get_deCBr(_mode,_h0mode,1),-0.1,0.1);     deCBr_201.setConstant(kTRUE);
  RooRealVar sCBr_201("sCBr_201","sCBr_201",get_sCBr(_mode,_h0mode,1),0.,0.5);           sCBr_201.setConstant(kTRUE);
  RooRealVar nr_201("nr_201","nr_201",get_nr(_mode,_h0mode,1),0.,100.);                  nr_201.setConstant(kTRUE);
  RooRealVar alphar_201("alphar_201","alphar_201",get_alphar(_mode,_h0mode,1),-10.,10.); alphar_201.setConstant(kTRUE);
  RooCBShape CBl_201("CBl_201","CBl_201",de,deCBl_201,sCBl_201,alphal_201,nl_201);
  RooCBShape CBr_201("CBr_201","CBr_201",de,deCBr_201,sCBr_201,alphar_201,nr_201);

  RooRealVar fCBl_201("fCBl_201","fCBl_201",get_fCBl(_mode,_h0mode,1),0.,1.); fCBl_201.setConstant(kTRUE);
  RooRealVar fCBr_201("fCBr_201","fCBr_201",get_fCBr(_mode,_h0mode,1),0.,1.); fCBr_201.setConstant(kTRUE);

  RooAddPdf pdf_de_201("pdf_de_201","pdf_de_201",RooArgList(CBl_201,CBr_201,g1_201),RooArgSet(fCBl_201,fCBr_201));

  RooRealVar de0_205("de0_205","de0_205",get_de0(_mode,_h0mode,5),-0.1,0.1);  de0_205.setConstant(kTRUE);
  RooRealVar s1_205("s1_205","s1_205",get_s1(_mode,_h0mode,5),0.,0.5);        s1_205.setConstant(kTRUE);
  RooGaussian g1_205("g1_205","g1_205",de,de0_205,s1_205);

  RooRealVar deCBl_205("deCBl_205","deCBl_205",get_deCBl(_mode,_h0mode,5),-0.1,0.1);     deCBl_205.setConstant(kTRUE);
  RooRealVar sCBl_205("sCBl_205","sCBl_205",get_sCBl(_mode,_h0mode,5),0.,0.5);           sCBl_205.setConstant(kTRUE);
  RooRealVar nl_205("nl_205","nl_205",get_nl(_mode,_h0mode,5),0.,100.);                  nl_205.setConstant(kTRUE);
  RooRealVar alphal_205("alphal_205","alphal_205",get_alphal(_mode,_h0mode,5),-10.,10.); alphal_205.setConstant(kTRUE);
  RooCBShape CBl_205("CBl_205","CBl_205",de,deCBl_205,sCBl_205,alphal_205,nl_205);

  RooRealVar fCBl_205("fCBl_205","fCBl_205",get_fCBl(_mode,_h0mode,5),0.,1.); fCBl_205.setConstant(kTRUE);

  RooAddPdf pdf_de_205("pdf_de_205","pdf_de_205",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));

  RooRealVar f_201("f_201","f_201",get_f201(_mode,_h0mode),0.,1.); if(cSIG) f_201.setConstant(kTRUE);
  RooAddPdf pdf_de_sig("pdf_de_sig","pdf_de_sig",RooArgList(pdf_de_201,pdf_de_205),RooArgSet(f_201));
  }

  /////////////
  // mbc pdf //
  /////////////
  const bool cond_flag = false;//(b0f == 5);
  if(cond_flag){
    RooRealVar cond("cond","cond",mr_cond_3,-0.1.,1.);      if(cRHO) cond.setConstant(kTRUE);
    RooRealVar condr("condr","condr",mr_condr_3,-1000.,1.); if(cRHO) condr.setConstant(kTRUE);

    RooRealVar _sl("_sl","_sl",get_sl(_mode,_h0mode,_b0f),0.,0.5);     if(cSIG) _sl.setConstant(kTRUE);
    RooRealVar _sr("_sr","_sr",get_sr(_mode,_h0mode,_b0f),0.0005,0.5); if(cSIG) _sr.setConstant(kTRUE);
    RooFormulaVar sl("sl","sl","@0+@1*@2",RooArgList(_sl,cond,de));
    RooFormulaVar sr("sr","sr","@0+@1*@2",RooArgList(_sr,condr,de));

    RooRealVar _sll("_sll","_sll",get_sll(_mode,_h0mode,_b0f),0.,0.5);     if(cSIG) _sll.setConstant(kTRUE);
    RooRealVar _srr("_srr","_srr",get_srr(_mode,_h0mode,_b0f),0.0005,0.5); if(cSIG) _srr.setConstant(kTRUE);
    RooFormulaVar sll("sll","sll","@0+@1*@2",RooArgList(_sll,cond,de));
    RooFormulaVar srr("srr","srr","@0+@1*@2",RooArgList(_srr,condr,de));
  } else{
    RooRealVar sl("sl","sl",get_sl(_mode,_h0mode,_b0f),0.,0.5);         if(cSIG) sl.setConstant(kTRUE);
    RooRealVar sr("sr","sr",get_sr(_mode,_h0mode,_b0f),0.0005,0.5);     if(cSIG) sr.setConstant(kTRUE);
    RooRealVar sll("sll","sll",get_sll(_mode,_h0mode,_b0f),0.,0.5);     if(cSIG) sll.setConstant(kTRUE);
    RooRealVar srr("srr","srr",get_srr(_mode,_h0mode,_b0f),0.0005,0.5); if(cSIG) srr.setConstant(kTRUE);
  }

  RooRealVar mbc0("mbc0","mbc0",get_mbc0(_mode,_h0mode,_b0f),5.26,5.30); if(cSIG) mbc0.setConstant(kTRUE);
  RooBifurGauss bg("bg","bg",mbc,mbc0,sl,sr);
  RooRealVar mbc00("mbc00","mbc00",get_mbc00(_mode,_h0mode,_b0f),5.26,5.30); if(cSIG) mbc00.setConstant(kTRUE);
  if(_b0f != 5 || _mode == 3) RooBifurGauss bgg("bgg","bgg",mbc,mbc00,sll,srr);
  else RooBifurGauss bgg("bgg","bgg",mbc,mbc0,sll,srr);

  RooRealVar fmbc("fmbc","fmbc",get_fmbc(_mode,_h0mode,_b0f),0.,1.); if(cSIG) fmbc.setConstant(kTRUE);
  RooAddPdf pdf_mbc_sig("pdf_mbc_sig","pdf_mbc_sig",RooArgList(bg,bgg),RooArgSet(fmbc));

  /////////
  // pdf //
  /////////
  if(!cond_flag) RooProdPdf pdf_sig("pdf_sig","pdf_sig",RooArgList(pdf_de_sig,pdf_mbc_sig));
  else RooProdPdf pdf_sig("pdf_sig","pdf_sig",pdf_de_sig,Conditional(pdf_mbc,mbc));

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar c1("c1","c1",mc_c1,-10.,10.);// if(cComb || true) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",mc_c2,-10.,10.);// if(cComb || true) c2.setConstant(kTRUE);
  RooChebychev pdf_de_comb("pdf_de_comb","pdf_de_comb",de,RooArgSet(c1,c2));
//  RooRealVar c3("c3","c3",0,-100.,10.);// if(cComb || true) c2.setConstant(kTRUE);
//  RooExponential pdf_de_exp("pdf_de_exp","pdf_de_exp",de,c3);
//  RooRealVar fexp("fexp","fexp",0.1,0.,1.);
//  RooAddPdf pdf_de_comb("pdf_de_comb","pdf_de_comb",RooArgSet(pdf_de_exp,pdf_de_cheb),RooArgList(fexp));

  /////////////
  // mbc pdf //
  /////////////
  RooRealVar argpar("argpar","argus shape parameter",mc_argpar,-100.,-1.); if(cComb) argpar.setConstant(kTRUE);
  RooRealVar argedge("argedge","argedge",mc_argedge,5.285,5.292); argedge.setConstant(kTRUE);
  RooArgusBG pdf_mbc_comb("pdf_mbc_comb","Argus PDF",mbc,argedge,argpar);

  //////////////
  // pdf comb //
  //////////////
  RooProdPdf pdf_comb("pdf_comb","pdf_comb",RooArgList(pdf_de_comb,pdf_mbc_comb));
  
  /////////////////
  // BB comb PDF //
  /////////////////
  /////////////
  // mbc pdf //
  /////////////
  RooRealVar edge("edge","edge",5.29,5.28,5.30,"GeV");                      edge.setConstant(kTRUE);
  RooRealVar mbctau("mbctau","mbctau",m_mbctau_bbcomb_eta10,-100,0.,"GeV"); mbctau.setConstant(kTRUE);
  RooRealVar pow1("pow1","pow1",m_pow1_bbcomb_eta10,0.1,10);                pow1.setConstant(kTRUE);
  RooRealVar pow2("pow2","pow2",m_pow2_bbcomb_eta10,0.1,10);                pow2.setConstant(kTRUE);

  RooGenericPdf pdf_mbc_bb_comb("pdf_mbc_bb_comb","pow(@1-@0,@2)*exp(@3*pow(@1-@0,@4))",RooArgList(mbc,edge,pow1,mbctau,pow2));

  ////////////
  // de pdf //
  ////////////
//  RooRealVar d1("d1","d1",m_c1_bbcomb_eta10,-10.,10.); if(cComb && false) d1.setConstant(kTRUE);
//  RooRealVar d2("d2","d2",m_c2_bbcomb_eta10,-10.,10.); if(cComb && false) d2.setConstant(kTRUE);
//  RooChebychev pdf_de_bb_comb("pdf_de_bb_comb","pdf_de_bb_comb",de,RooArgSet(c1,c2));

  RooProdPdf pdf_bb_comb("pdf_bb_comb","pdf_bb_comb",RooArgList(pdf_mbc_bb_comb,pdf_de_comb));


  //////////////////
  // Complete PDF //
  //////////////////
  RooRealVar Nsig("Nsig","Nsig",200,10.,10000.);
  RooRealVar Ncmb("Ncmb","Ncmb",1000,10,100000);
  RooRealVar NBBcmb("NBBcmb","NBBcmb",300,10,100000.);
  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_comb,pdf_bb_comb),RooArgList(Nsig,Ncmb,NBBcmb));

  RooArgSet* params = pdf.getParameters(RooArgSet(de,mbc));
//  RooArgset* initParams = (RooArgSet*) params->snapshot();
  
  pdf.fitTo(ds,Verbose(),Timer(true));
  
  params->printLatex(OutputFile("PurityEtaGGFit.tex"));

   RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intRho  = pdf_bb_comb->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   const double nsig = intSig->getVal()*Nsig.getVal();
   const double nsig_err = intSig->getVal()*Nsig.getError();
   const double nsig_err_npq = TMath::Sqrt(nsig*(Nsig.getVal()-nsig)/Nsig.getVal());
   const double nsig_err_total = TMath::Sqrt(nsig_err*nsig_err+nsig_err_npq*nsig_err_npq);
   const double nrho = intRho->getVal()*NBBcmb.getVal();
   const double nrho_err = intRho->getVal()*NBBcmb.getError();
   const double nrho_err_npq = TMath::Sqrt(nrho*(NBBcmb.getVal()-nrho)/NBBcmb.getVal());
   const double nrho_err_total = TMath::Sqrt(nrho_err*nrho_err+nrho_err_npq*nrho_err_npq);
   const double ncmb = intCmb->getVal()*Ncmb.getVal();
   const double ncmb_err = intCmb->getVal()*Ncmb.getError();
   const double ncmb_err_npq = TMath::Sqrt(ncmb*(Ncmb.getVal()-ncmb)/Ncmb.getVal());
   const double ncmb_err_total = TMath::Sqrt(ncmb_err*ncmb_err+ncmb_err_npq*ncmb_err_npq);
   const double purity = nsig/(nsig+nrho+ncmb);
   const double purity_err = nsig_err_total/(nsig+nrho+ncmb);

   de.setRange("Ellips",de_min,de_max);
   RooFormulaVar mbclo("mbclo","(1-(@0-@3)/@4*(@0-@3)/@4) > 0 ? @1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4) : 0",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   RooFormulaVar mbchi("mbchi","(1-(@0-@3)/@4*(@0-@3)/@4) > 0 ? @1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4) : 0",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   mbc.setRange("Ellips",mbclo,mbchi);
   
   de.setRange("Elli",de_min,de_max);
   RooFormulaVar mbclo1("mbclo1","(1-(@0-@3)/@4*(@0-@3)/@4) > 0 ? @1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4) : 0",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   RooFormulaVar mbchi1("mbchi1","(1-(@0-@3)/@4*(@0-@3)/@4) > 0 ? @1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4) : 0",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   mbc.setRange("Elli",mbclo1,mbchi1);

   RooAbsReal* intSigEl = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intRhoEl = pdf_bb_comb->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intCmbEl = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   const double nsigEl = intSigEl->getVal()*Nsig.getVal();
   const double nsig_errEl = intSigEl->getVal()*Nsig.getError();
   const double nsig_errEl_npq = TMath::Sqrt(nsigEl*(Nsig.getVal()-nsigEl)/Nsig.getVal());
   const double nsig_errEl_total = TMath::Sqrt(nsig_errEl*nsig_errEl+nsig_errEl_npq*nsig_errEl_npq);
   const double nrhoEl = intRhoEl->getVal()*NBBcmb.getVal();
   const double nrho_errEl = intRhoEl->getVal()*NBBcmb.getError();
   const double nrho_errEl_npq = TMath::Sqrt(nrhoEl*(NBBcmb.getVal()-nrhoEl)/NBBcmb.getVal());
   const double nrho_errEl_total = TMath::Sqrt(nrho_errEl*nrho_errEl+nrho_errEl_npq*nrho_errEl_npq);
   const double ncmbEl = intCmbEl->getVal()*Ncmb.getVal();
   const double ncmb_errEl = intCmbEl->getVal()*Ncmb.getError();
   const double ncmb_errEl_npq = TMath::Sqrt(ncmbEl*(Ncmb.getVal()-ncmbEl)/Ncmb.getVal());
   const double ncmb_errEl_total = TMath::Sqrt(ncmb_errEl*ncmb_errEl+ncmb_errEl_npq*ncmb_errEl_npq);
   const double purityEl = nsigEl/(nsigEl+nrhoEl+ncmbEl);
   const double purity_errEl = nsig_errEl_total/(nsigEl+nrhoEl+ncmbEl);

   RooAbsReal* intSigEl1 = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   RooAbsReal* intRhoEl1 = pdf_bb_comb->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   RooAbsReal* intCmbEl1 = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   const double nsigEl1 = intSigEl1->getVal()*Nsig.getVal();
   const double nsig_errEl1 = intSigEl1->getVal()*Nsig.getError();
   const double nsig_errEl1_npq = TMath::Sqrt(nsigEl1*(Nsig.getVal()-nsigEl1)/Nsig.getVal());
   const double nsig_errEl1_total = TMath::Sqrt(nsig_errEl1*nsig_errEl1+nsig_errEl1_npq*nsig_errEl1_npq);
   const double nrhoEl1 = intRhoEl1->getVal()*NBBcmb.getVal();
   const double nrho_errEl1 = intRhoEl1->getVal()*NBBcmb.getError();
   const double nrho_errEl1_npq = TMath::Sqrt(nrhoEl1*(NBBcmb.getVal()-nrhoEl1)/NBBcmb.getVal());
   const double nrho_errEl1_total = TMath::Sqrt(nrho_errEl1*nrho_errEl1+nrho_errEl1_npq*nrho_errEl1_npq);
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
//  pdf.plotOn(deFrame,Components(pdf_bb_comb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(RooArgSet(pdf_comb,pdf_bb_comb)),LineStyle(kDashed),ProjectionRange("mbcSignal"));
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

  stringstream out1;
  TPaveText *pt = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out1.str("");
  out1 << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out1.str().c_str());
  out1.str("");
  if(!data) out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  else      out1 << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  pt->AddText(out1.str().c_str());
  out1.str("");
  if(!data) out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  else out1 << "Purity: " << std::fixed << std::setprecision(2) << purityEl*100. << " #pm " << purity_errEl*100;
  pt->AddText(out1.str().c_str());
  pt->Draw();

  TLine *de_line_RIGHT;
  if(!data) de_line_RIGHT = new TLine(de_max,0,de_max,120);
  else      de_line_RIGHT = new TLine(de_max,0,de_max,30);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  if(!data) de_line_LEFT = new TLine(de_min,0,de_min,120);
  else      de_line_LEFT = new TLine(de_min,0,de_min,30);
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
  pdf.plotOn(mbcFrame,Components(pdf_comb),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_sig),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf.plotOn(mbcFrame,Components(pdf_bb_comb),LineStyle(kDashed),ProjectionRange("deSignal"));
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

  TPaveText *ptmbc = new TPaveText(0.2,0.75,0.58,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out1.str("");
  out1 << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out1.str().c_str());
  out1.str("");
  if(!data) out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  else      out1 << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  ptmbc->AddText(out1.str().c_str());
  out1.str("");
  if(!data) out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  else out1 << "Purity: " << std::fixed << std::setprecision(2) << purityEl*100. << " #pm " << purity_errEl*100;
  ptmbc->AddText(out1.str().c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT;
  if(!data) mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,70);
  else      mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,40);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT;
  if(!data) mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,70);
  else      mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,40);
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
  elli1->SetLineColor(kBlue);
  elli1->SetLineWidth(2);
  TLine* l1 = new TLine(de_min,mbc_min,de_max,mbc_min);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->SetLineWidth(2);
  TLine* l2 = new TLine(de_min,mbc_max,de_max,mbc_max);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(1);
  l2->SetLineWidth(2);
  TLine* l3 = new TLine(de_min,mbc_min,de_min,mbc_max);
  l3->SetLineColor(kRed);
  l3->SetLineStyle(1);
  l3->SetLineWidth(2);
  TLine* l4 = new TLine(de_max,mbc_min,de_max,mbc_max);
  l4->SetLineColor(kRed);
  l4->SetLineStyle(1);
  l4->SetLineWidth(2);

  TCanvas* ellican = new TCanvas("ellican","ellican",400,400);
  ellican->cd();
  out.str("");
  out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX << " && de>-0.15 && de<0.20 && mbc>5.265";
  tree->Draw("mbc:de",out.str().c_str());
  elli->Draw(); elli1->Draw(); l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
  
  if(!data){
    TCanvas* sigcan = new TCanvas("sigcan","sigcan",400,400);
    sigcan->cd();
    out.str("");
    out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX << " && de>-0.15 && de<0.20 && mbc>5.265 && (b0f == 1 || b0f == 5 || b0f == 10)";
    tree->Draw("mbc:de",out.str().c_str());
    elli->Draw(); elli1->Draw(); l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
  
    TCanvas* backcan = new TCanvas("backcan","backcan",400,400);
    backcan->cd();
    out.str("");
    out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX << " && de>-0.15 && de<0.20 && mbc>5.265 && !(b0f == 1 || b0f == 5 || b0f == 10 || b0f == 0)";
    tree->Draw("mbc:de",out.str().c_str());
    elli->Draw(); elli1->Draw(); l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
    
    cout << "Rectangle:" << endl;
    out.str("");
    out << "de<" << de_max << " && de>" << de_min;
    out << " && mbc>" << mbc_min << " && mbc<" << mbc_max;
    Roo1DTable* recttable = ds.table(b0f,out.str().c_str());
    recttable->Print();
    recttable->Print("v");

    cout << "Ellips:" << endl;
    out.str("");
    out << "(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius.getVal() << "<1";
    cout << out.str() << endl;
    Roo1DTable* ellitable = ds.table(b0f,out.str().c_str());
    ellitable->Print();
    ellitable->Print("v");
     
    cout << "Elli:" << endl;
    out.str("");
    out << "(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "*(de-" << de_center.getVal() << ")/" << de_radius1.getVal() << "+(mbc-"<<mbc_center.getVal()<<")/" << mbc_radius1.getVal() << "*(mbc-" << mbc_center.getVal() << ")/" << mbc_radius1.getVal() << "<1";
    cout << out.str() << endl;
    Roo1DTable* ellitable1 = ds.table(b0f,out.str().c_str());
    ellitable1->Print();
    ellitable1->Print("v");

    Roo1DTable* fulltable = ds.table(b0f);
    fulltable->Print();
    fulltable->Print("v");
  }

  cout << "Rectangle:" << endl;
  cout << "Nsig    = " << nsig <<" +- " << nsig_err << " +- " << nsig_err_npq << " (" << nsig_err_total << ")" << endl;
  cout << "NBBcomb = " << nrho <<" +- " << nrho_err << " +- " << nrho_err_npq << " (" << nrho_err_total << ")" << endl;
  cout << "Ncmb    = " << ncmb <<" +- " << ncmb_err << " +- " << ncmb_err_npq << " (" << ncmb_err_total << ")" << endl;
  cout << "Pury    = " << purity << " +- " << purity_err << endl;

  cout << "Ellips:" << endl;
  cout << "Nsig    = " << nsigEl <<" +- " << nsig_errEl << " +- " << nsig_errEl_npq << " (" << nsig_errEl_total << ")" << endl;
  cout << "NBBcomb = " << nrhoEl <<" +- " << nrho_errEl << " +- " << nrho_errEl_npq << " (" << nrho_errEl_total << ")" << endl;
  cout << "Ncmb    = " << ncmbEl <<" +- " << ncmb_errEl << " +- " << ncmb_errEl_npq << " (" << ncmb_errEl_total << ")" << endl;
  cout << "Pury    = " << purityEl << " +- " << purity_errEl << endl;
  
  cout << "Elli:" << endl;
  cout << "Nsig    = " << nsigEl1 <<" +- " << nsig_errEl1 << " +- " << nsig_errEl1_npq << " (" << nsig_errEl1_total << ")" << endl;
  cout << "NBBcomb = " << nrhoEl1 <<" +- " << nrho_errEl1 << " +- " << nrho_errEl1_npq << " (" << nrho_errEl1_total << ")" << endl;
  cout << "Ncmb    = " << ncmbEl1 <<" +- " << ncmb_errEl1 << " +- " << ncmb_errEl1_npq << " (" << ncmb_errEl1_total << ")" << endl;
  cout << "Pury    = " << purityEl1 << " +- " << purity_errEl1 << endl;
}

