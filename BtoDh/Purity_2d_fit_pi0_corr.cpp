#include "cuts.h"
using namespace RooFit;

void Purity_2d_fit_pi0_corr(int type = 0){
  TChain* tree = new TChain("TEvent");
//  if(!type) tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_3-5.root");
//  else      tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_data.root");
  tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1_full.root");

  gROOT->ProcessLine(".L pdfs/RooRhoDeltaEPdf.cxx+");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
//  b0f.defineType("fsr",10);
//  b0f.defineType("bad_pi0",5);
//  b0f.defineType("rho",3);
//  b0f.defineType("comb",-1);

  RooArgSet argset;

  const double mbcMin = 5.20;
  const double mbcMax = mbc_max;
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

  RooRealVar mbc("mbc","M_{bc}",0.5*(mbc_min+mbc_max),mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);

  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_min,de_max);
  
  const double BDTG_MIN = bdtg_cut_pi0;
  const double BDTG_MAX = 1;
//   de.setRange("Ellips",de_min,de_max);
//   RooFormulaVar mbclo("mbclo","@1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.00001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
//   RooFormulaVar mbchi("mbchi","@1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.00001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
//   mbc.setRange("Ellips",mbclo,mbchi);
  
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",bdtg_cut_pi0,1.); argset.add(bdtg);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  if(!type) argset.add(b0f);
  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");
//  RooDataSet* ds0 = ds.reduce(RooArgSet(de,mbc));

  stringstream out;
  if(!type){
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
  RooRealVar de0("de0","de0",m_de0,-0.1,0.1); if(cSig) de0.setConstant(kTRUE);
  RooRealVar s1("s1","s1",m_s1,0.,0.5);       if(cSig) s1.setConstant(kTRUE);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar deCBl("deCBl","deCBl",m_deCBl,-0.1,0.1);     if(cSig) deCBl.setConstant(kTRUE);
  RooRealVar sCBl("sCBl","sCBl",m_sCBl,0.,0.5);           if(cSig) sCBl.setConstant(kTRUE);
  RooRealVar nl("nl","nl",m_nl,0.,100.);                  if(cSig) nl.setConstant(kTRUE);
  RooRealVar alphal("alphal","alphal",m_alphal,-10.,10.); if(cSig) alphal.setConstant(kTRUE);

  RooRealVar deCBr("deCBr","deCBr",m_deCBr,-0.1,0.1);     if(cSig) deCBr.setConstant(kTRUE);
  RooRealVar sCBr("sCBr","sCBr",m_sCBr,0.,0.5);           if(cSig) sCBr.setConstant(kTRUE);
  RooRealVar nr("nr","nr",m_nr,0.,100.);                  if(cSig) nr.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",m_alphar,-10.,10.); if(cSig) alphar.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",m_fCBl,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
  RooRealVar fCBr("fCBr","fCBr",m_fCBr,0.,1.); if(cSig) fCBr.setConstant(kTRUE);

  RooAddPdf pdf_de_sig("pdf_de_sig","pdf_de_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));

  RooConstVar DELO("DELO","DELO",-0.15);
  /////////////
  // mbc pdf //
  /////////////
  RooRealVar c1_mbc0("c1_mbc0","c1_mbc0",m_mbc0_c1_pi0,-0.1,0.); c1_mbc0.setConstant(kTRUE);
  RooRealVar c2_mbc0("c2_mbc0","c2_mbc0",m_mbc0_c2_pi0,0.,1.);   c2_mbc0.setConstant(kTRUE);
  RooRealVar c3_mbc0("c3_mbc0","c3_mbc0",m_mbc0_c3_pi0,0.,10.);  c3_mbc0.setConstant(kTRUE);
  RooRealVar mbc0("mbc0","mbc0",m_mbc0_c0_pi0,5.26,5.30); if(cSIG) mbc0.setConstant(kTRUE);
//  RooFormulaVar _mbc0("_mbc0","_mbc0","abs(@0+@1*@2+@1*@1*@3+@1*@1*@1*@4) > abs(@0+@5*@2+@5*@5*@3+@5*@5*@5*@4) ? (@0+@5*@2+@5*@5*@3+@5*@5*@5*@4) : @0+@1*@2+@1*@1*@3+@1*@1*@1*@4",RooArgList(mbc0,de,c1_mbc0,c2_mbc0,c3_mbc0,DELO));
  RooFormulaVar _mbc0("_mbc0","_mbc0","abs(@1*@2+@1*@1*@3+@1*@1*@1*@4) > abs(@5*@2+@5*@5*@3+@5*@5*@5*@4) ? @0 : @0+@1*@2+@1*@1*@3+@1*@1*@1*@4",RooArgList(mbc0,de,c1_mbc0,c2_mbc0,c3_mbc0,DELO));

  RooRealVar c1_alpha("c1_alpha","c1_alpha",m_alpha_c1_pi0,-2.,0.); c1_alpha.setConstant(kTRUE);
  RooRealVar c2_alpha("c2_alpha","c2_alpha",m_alpha_c2_pi0,10.,100.);c2_alpha.setConstant(kTRUE);
  RooRealVar c3_alpha("c3_alpha","c3_alpha",m_alpha_c3_pi0,10.,100.);c3_alpha.setConstant(kTRUE);
  RooRealVar alpha("alpha","alpha",m_alpha_c0_pi0,0.,0.1); if(cSIG) alpha.setConstant(kTRUE);
  RooFormulaVar _alpha("_alpha","_alpha","abs(@0+@1*@2+@1*@1*@3+@1*@1*@1*@4) > abs(@0+@5*@2+@5*@5*@3+@5*@5*@5*@4) ? (@0+@5*@2+@5*@5*@3+@5*@5*@5*@4) : (@0+@1*@2+@1*@1*@3+@1*@1*@1*@4)",RooArgList(alpha,de,c1_alpha,c2_alpha,c3_alpha,DELO));
  RooRealVar c1_width("c1_width","c1_width",m_width_c1_pi0,-2.,0.);  c1_width.setConstant(kTRUE);
  RooRealVar c2_width("c2_width","c2_width",m_width_c2_pi0,10.,100.);c2_width.setConstant(kTRUE);
  RooRealVar width("width","width",m_width_c0_pi0,0.,0.1); if(cSIG) width.setConstant(kTRUE);
  RooFormulaVar _width("_width","_width","@0+@1*@2+@1*@1*@3",RooArgList(width,de,c1_width,c2_width));

  RooNovosibirsk pdf_mbc_sig("pdf_mbc_sig","pdf_mbc_sig",mbc,_mbc0,_width,_alpha);

  /////////
  // pdf //
  /////////
  RooProdPdf pdf_sig("pdf_sig","pdf_sig",pdf_de_sig,Conditional(pdf_mbc_sig,mbc));

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar c1("c1","c1",mc_c1,-10.,10.); if(cComb) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",mc_c2,-10.,10.); if(cComb) c2.setConstant(kTRUE);
  //RooChebychev pdf_de_comb("pdf_de_comb","pdf_de_comb",de,RooArgSet(c1,c2));

  RooRealVar c20("c20","c20",m_c20_pi0,-10.,10.);c20.setConstant(kTRUE);
  RooRealVar c21("c21","c21",m_c21_pi0,-10.,10); c21.setConstant(kTRUE);
  RooRealVar c22("c22","c22",m_c22_pi0,-10.,10); c22.setConstant(kTRUE);
  RooRealVar c23("c23","c23",m_c23_pi0,-10.,10); c23.setConstant(kTRUE);
//  RooFormulaVar _c2("_c2","@0+@1*@4+@2*@4*@4+@3*@4*@4*@4",RooArgSet(c20,c21,c22,c23,mbc));
  RooRealVar c10("c10","c10",m_c10_pi0,-10.,10.);c10.setConstant(kTRUE);
  RooRealVar c11("c11","c11",m_c11_pi0,-10.,10); c11.setConstant(kTRUE);
  RooRealVar c12("c12","c12",m_c12_pi0,-10.,10); c12.setConstant(kTRUE);
//  RooFormulaVar _c1("_c1","@0+@1*@3+@2*@3*@3",RooArgSet(c10,c11,c12,mbc));


  RooRealVar c10("c10","c10",m_c10_eta10,-10.,10.);c10.setConstant(kTRUE);
  RooRealVar c11("c11","c11",m_c11_eta10,-10.,10); c11.setConstant(kTRUE);
  RooRealVar c12("c12","c12",m_c12_eta10,-10.,10); c12.setConstant(kTRUE);
  RooFormulaVar _c1("_c1","@0+@1*@3+@2*@3*@3",RooArgSet(c10,c11,c12,mbc));
  RooRealVar c20("c20","c20",m_c20_eta10,-10.,10.);c20.setConstant(kTRUE);
  RooRealVar c21("c21","c21",m_c21_eta10,-10.,10); c21.setConstant(kTRUE);
  RooFormulaVar _c2("_c2","@0+@1*@2",RooArgSet(c20,c21,mbc));

  RooChebychev pdf_de_comb("pdf_de_comb","pdf_de_comb",de,RooArgSet(c1,_c2));
  /////////////
  // mbc pdf //
  /////////////
  RooRealVar argpar("argpar","argus shape parameter",m_argpar_pi0,-100.,-1.); if(cComb) argpar.setConstant(kTRUE);
  RooRealVar argedge("argedge","argedge",m_argedge_pi0,5.285,5.292); argedge.setConstant(kTRUE);
  RooArgusBG pdf_mbc_comb("pdf_mbc_comb","Argus PDF",mbc,argedge,argpar);

  /////////
  // pdf //
  /////////
//  RooProdPdf pdf_comb("pdf_comb","pdf_comb",RooArgList(pdf_de_comb,pdf_mbc_comb));
  RooProdPdf pdf_comb("pdf_comb","pdf_comb",pdf_mbc_comb,Conditional(pdf_de_comb,de));
  
  /////////////
  // Rho PDF //
  /////////////
  if(keysflag){
//    TFile* file = TFile::Open("rho_pdf.root");
//    RooHistPdf* pdf_rho = (RooHistPdf*) file->Get("pdf");

    TFile *ifile1 = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-2.root");
    TTree *tree1 = (TTree*)ifile1->Get("TEvent");
    RooArgSet argset1;
    RooCategory b0f("b0f","b0f");
    b0f.defineType("rho",3);
    argset1.add(b0f);
    argset1.add(mbc);
    argset1.add(de);
    argset1.add(md);
    argset1.add(mk);
    argset1.add(mpi0);
    argset1.add(bdtgs);
    argset1.add(atckpi_max);
    RooDataSet ds1("ds1","ds1",tree1,argset1,"mbc>0||mbc<=0");
    RooNDKeysPdf* pdf_rho1 = new RooNDKeysPdf("pdf_rho1","pdf_rho1",RooArgList(de,mbc),ds1,"am",2);
    de.setBins(250);
    mbc.setBins(250);
    RooDataHist* hist = new RooDataHist("hist","hist",RooArgSet(de,mbc));
    pdf_rho1->fillDataHist(hist,new RooArgSet(de,mbc),1);
    RooHistPdf* pdf_rho = new RooHistPdf("pdf_rho","pdf_rho",RooArgSet(de,mbc),*hist);
    de.setBins(100);
    mbc.setBins(100);
  }
  else{
  ////////////
  // de pdf //
  ////////////
  if(de_rho_param == 0){
  RooRealVar exppar("exppar","exppar",mr_exppar,-40.,-25.);// if(cRho) exppar.setConstant(kTRUE);
  RooExponential pdf_de_rho("pdf_de_rho","pdf_de_rho",de,exppar);
  }
  
  RooRealVar de0r("de0r","de0r",mr_de0r,-0.2,0.12); if(cRho) de0r.setConstant(kTRUE);
  
  if(de_rho_param == 1){
   RooRealVar slopel("slopel","slopel",mr_slopel,-1000,-500.); if(cRho) slopel.setConstant(kTRUE);
   RooRealVar sloper("sloper","sloper",mr_sloper,-10000,0.); if(cRho) sloper.setConstant(kTRUE);
   RooRealVar steep("steep","steep",mr_steep,7.,9.); if(cRho) steep.setConstant(kTRUE);
   RooRealVar p5("p5","p5",mr_p5,0.01,1000.); if(cRho) p5.setConstant(kTRUE);
   RooRhoDeltaEPdf pdf_de_rho("pdf_de_rho","pdf_de_rho",de,de0r,slopel,sloper,steep,p5);
  }
  
  if(de_rho_param == -1){
   RooRealVar x0("x0","x0",mr_x0,-0.2,0.12); if(cRho) x0.setConstant(kTRUE);
   RooRealVar p1("p1","p1",mr_p1,-1000.,100.); if(cRho) p1.setConstant(kTRUE);
   RooRealVar p2("p2","p2",mr_p2,0.,100.); if(cRho) p2.setConstant(kTRUE);
   RooGenericPdf pdf_de_rho("pdf_de_rho","1+@0*@1-@2*TMath::Log(1+TMath::Exp(@2*(@0-@1)/@3))",RooArgSet(de,x0,p1,p2));
  }
  
  if(de_rho_param == 2){
     RooRealVar slopel("slopel","slopel",mr_slopel_2,-10,100.); if(cRho) slopel.setConstant(kTRUE);
     RooRealVar steep("steep","steep",mr_steep_2,0.,10000.); if(cRho) steep.setConstant(kTRUE);
     RooRealVar exppar("exppar","exppar",mr_exppar_2,1.,100.); if(cRho) exppar.setConstant(kTRUE);
     RooRealVar x0("x0","x0",mr_x0_2,-0.2,0.12); if(cRho) x0.setConstant(kTRUE);
     RooGenericPdf pdf_de_rho("pdf_de_rho","pdf_de_rho","1+@2*(@0-@1)+@3*TMath::Log(1+TMath::Exp(-@2*(@0-@1)+TMath::Exp(-@4*(@0-@1))))",RooArgSet(de,x0,slopel,steep,exppar));
   }
   
  /////////////
  // mbc pdf //
  /////////////
  if(mbc_rho_param == 0){
  RooRealVar rmbc0("rmbc0","rmbc0",mr_mbc0_0,5.26,5.30); if(cRho) rmbc0.setConstant(kTRUE);
  RooRealVar cond("cond","cond",mr_cond_0,-1000.,1.); if(cRho || true) cond.setConstant(kTRUE);
  RooRealVar condr("condr","condr",mr_condr_0,-1000.,1.); if(cRho || true) condr.setConstant(kTRUE);
  RooRealVar rsl("rsl","rsl",mr_sl_0,0.,0.5); if(cRho) rsl.setConstant(kTRUE);
  RooRealVar rsr("rsr","rsr",mr_sr_0,0.,0.5); if(cRho) rsr.setConstant(kTRUE);
  RooFormulaVar _sl("_sl","_sl","@0+@1*@2",RooArgList(rsl,cond,de));
  RooFormulaVar _sr("_sr","_sr","@0+@1*@2",RooArgList(rsr,condr,de));
  RooBifurGauss pdf_mbc_rho("pdf_mbc_rho","pdf_mbc_rho",mbc,rmbc0,_sl,_sr);
  }

  if(mbc_rho_param == 1){
   RooRealVar cond("cond","cond",mr_cond_1,-1000.,1.); if(cRho || true) cond.setConstant(kTRUE);
   RooRealVar rmbc0("rmbc0","rmbc0",mr_mbc0_1,5.26,5.30); if(cRho) rmbc0.setConstant(kTRUE);
   RooRealVar rsl("rsl","rsl",mr_sl_1,0.,0.5); if(cRho) rsl.setConstant(kTRUE);
   RooFormulaVar _rsl("_rsl","_rsl","@0+@1*(@2-@3)",RooArgList(rsl,cond,de,de0r));
   RooRealVar rsr("rsr","rsr",mr_sr_1,0.,0.5); if(cRho) rsr.setConstant(kTRUE);
   RooBifurGauss rbg("rbg","rbg",mbc,rmbc0,_rsl,rsr);

   RooRealVar rmbc00("rmbc00","rmbc00",mr_mbc00_1,5.26,5.30); if(cRho) rmbc00.setConstant(kTRUE);
   RooRealVar rsll("rsll","rsll",mr_sll_1,0.,0.5); if(cRho) rsll.setConstant(kTRUE);
   RooRealVar rsrr("rsrr","rsrr",mr_srr_1,0.,0.5); if(cRho) rsrr.setConstant(kTRUE);
   RooFormulaVar _rsll("_rsll","_rsll","@0+@1*(@2-@3)",RooArgList(rsll,cond,de,de0r));
   RooFormulaVar _rsrr("_rsrr","_rsrr","@0+@1*(@2-@3)",RooArgList(rsrr,cond,de,de0r));
   RooBifurGauss rbgg("rbgg","rbgg",mbc,rmbc00,_rsll,_rsrr);
 
   RooRealVar rfmbc("rfmbc","rfmbc",mr_fmbc_1,0.,1.); if(cRho) rfmbc.setConstant(kTRUE);
   RooAddPdf pdf_mbc_rho("pdf_mbc_rho","pdf_mbc_rho",RooArgList(rbg,rbgg),RooArgSet(rfmbc));
  }

  if(mbc_rho_param == 2){
   RooRealVar mbc1("mbc1","mbc1",mr_mbc1_2,5.26,5.30); if(cRho)
   mbc1.setConstant(kTRUE);
   RooRealVar ss("ss","ss",mr_ss_2,0.,0.5); if(cRho) ss.setConstant(kTRUE);
   RooRealVar alpha("alpha","alpha",mr_alpha_2,-3,0.); if(cRho) alpha.setConstant(kTRUE);
   RooRealVar n("n","n",mr_n_2,0.1,10.); n.setConstant(kTRUE);
   RooCBShape cb("cb","cb",mbc,mbc1,ss,alpha,n);

   RooRealVar rmbc0("rmbc0","rmbc0",mr_mbc0_2,5.26,5.30); if(cRho) rmbc0.setConstant(kTRUE);
   RooRealVar rsr("rsr","rsr",mr_sr_2,0.,0.5); if(cRho) rsr.setConstant(kTRUE);
   RooGaussian gauss("gauss","gauss",mbc,rmbc0,rsr);
   RooRealVar fcb("fcb","fcb",mr_fcb_2,0.,1.); if(cRho) fcb.setConstant(kTRUE);
   RooAddPdf pdf_mbc_rho("pdf_mbc_rho","pdf_mbc_rho",RooArgList(cb,gauss),RooArgSet(fcb));
  }
  
  if(mbc_rho_param == 3){
   RooRealVar rargpar("rargpar","argus shape parameter",mr_argpar_3,-100.,0.); if(cRho) rargpar.setConstant(kTRUE);
   RooRealVar rargedge("rargedge","rargedge",mr_argedge_3,5.285,5.292); if(cRho) rargedge.setConstant(kTRUE);
   RooArgusBG rargus("rargus","Argus PDF",mbc,rargedge,rargpar);

   RooRealVar rmbc0("rmbc0","rmbc0",mr_mbc0_3,5.26,5.30); if(cRho) rmbc0.setConstant(kTRUE);
   RooRealVar rcond("rcond","rcond",mr_cond_3,-0.1.,1.); if(cRho) rcond.setConstant(kTRUE);
   RooRealVar rcondr("rcondr","rcondr",mr_condr_3,-1000.,1.); if(cRho) rcondr.setConstant(kTRUE);
   RooRealVar rsl("rsl","rsl",mr_sl_3,0.,0.5); if(cRho) rsl.setConstant(kTRUE);
   RooRealVar rsr("rsr","rsr",mr_sr_3,0.,0.5); if(cRho) rsr.setConstant(kTRUE);
   RooFormulaVar _rsl("_rsl","_rsl","@0+@1*@2",RooArgList(rsl,rcond,de));
   RooFormulaVar _rsr("_rsr","_rsr","@0+@1*@2",RooArgList(rsr,rcondr,de));
   RooBifurGauss rbg("rbg","rbg",mbc,rmbc0,_rsl,_rsr);
   RooProdPdf pdf_mbc_rho("pdf_mbc_rho","pdf_mbc_rho",RooArgList(rbg,rargus));
  }

  /////////
  // pdf //
  /////////
  RooProdPdf *pdf_rho = new RooProdPdf("pdf_rho","pdf_rho",pdf_de_rho,Conditional(pdf_mbc_rho,mbc));
//  RooProdPdf pdf_rho("pdf_rho","pdf_rho",pdf_de_rho,pdf_mbc_rho);
  }

  //////////////////
  // Complete PDF //
  //////////////////
//  RooRealVar fsig("fsig","fsig",0.1,0.,1.);// fsig.setConstant(kTRUE);
//  RooRealVar frho("frho","frho",0.1,0.,1.);// frho.setConstant(kTRUE);
  RooRealVar Nsig("Nsig","Nsig",1150,0.,10000.);
//  Nsig.setVal(1); Nsig.setConstant(kTRUE);
  RooRealVar Nrho("Nrho","Nrho",9308,0,100000.);
//  Nrho.setVal(1); Nrho.setConstant(kTRUE);
  RooRealVar Ncmb("Ncmb","Ncmb",21288,0,100000);
//  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_rho,pdf_comb),RooArgList(fsig,frho));
  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,*pdf_rho,pdf_comb),RooArgList(Nsig,Nrho,Ncmb));

  RooArgSet* params = pdf.getParameters(RooArgSet(de,mbc));
//  RooArgset* initParams = (RooArgSet*) params->snapshot();
  
  pdf.fitTo(ds,Verbose(),Timer(true));
  
  params->printLatex(OutputFile("PurityFit.tex"));

   RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intRho  = pdf_rho->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Signal"));
   const double nsig = intSig->getVal()*Nsig.getVal();
   const double nsig_err = intSig->getVal()*Nsig.getError();
   const double nsig_err_npq = TMath::Sqrt(nsig*(Nsig.getVal()-nsig)/Nsig.getVal());
   const double nsig_err_total = TMath::Sqrt(nsig_err*nsig_err+nsig_err_npq*nsig_err_npq);
   const double nrho = intRho->getVal()*Nrho.getVal();
   const double nrho_err = intRho->getVal()*Nrho.getError();
   const double nrho_err_npq = TMath::Sqrt(nrho*(Nrho.getVal()-nrho)/Nrho.getVal());
   const double nrho_err_total = TMath::Sqrt(nrho_err*nrho_err+nrho_err_npq*nrho_err_npq);
   const double ncmb = intCmb->getVal()*Ncmb.getVal();
   const double ncmb_err = intCmb->getVal()*Ncmb.getError();
   const double ncmb_err_npq = TMath::Sqrt(ncmb*(Ncmb.getVal()-ncmb)/Ncmb.getVal());
   const double ncmb_err_total = TMath::Sqrt(ncmb_err*ncmb_err+ncmb_err_npq*ncmb_err_npq);
   const double purity = nsig/(nsig+nrho+ncmb);
   const double purity_err = nsig_err_total/(nsig+nrho+ncmb);

   de.setRange("Ellips",de_min,de_max);
   RooFormulaVar mbclo("mbclo","@1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.0000001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   RooFormulaVar mbchi("mbchi","@1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.0000001)",RooArgSet(de,mbc_center,mbc_radius,de_center,de_radius));
   mbc.setRange("Ellips",mbclo,mbchi);
   
   de.setRange("Elli",de_min,de_max);
   RooFormulaVar mbclo1("mbclo1","@1-@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   RooFormulaVar mbchi1("mbchi1","@1+@2*TMath::Sqrt(1-(@0-@3)/@4*(@0-@3)/@4+0.0000001)",RooArgSet(de,mbc_center,mbc_radius1,de_center,de_radius1));
   mbc.setRange("Elli",mbclo1,mbchi1);

   RooAbsReal* intSigEl = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intRhoEl = pdf_rho->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   RooAbsReal* intCmbEl = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Ellips"));
   const double nsigEl = intSigEl->getVal()*Nsig.getVal();
   const double nsig_errEl = intSigEl->getVal()*Nsig.getError();
   const double nsig_errEl_npq = TMath::Sqrt(nsigEl*(Nsig.getVal()-nsigEl)/Nsig.getVal());
   const double nsig_errEl_total = TMath::Sqrt(nsig_errEl*nsig_errEl+nsig_errEl_npq*nsig_errEl_npq);
   const double nrhoEl = intRhoEl->getVal()*Nrho.getVal();
   const double nrho_errEl = intRhoEl->getVal()*Nrho.getError();
   const double nrho_errEl_npq = TMath::Sqrt(nrhoEl*(Nrho.getVal()-nrhoEl)/Nrho.getVal());
   const double nrho_errEl_total = TMath::Sqrt(nrho_errEl*nrho_errEl+nrho_errEl_npq*nrho_errEl_npq);
   const double ncmbEl = intCmbEl->getVal()*Ncmb.getVal();
   const double ncmb_errEl = intCmbEl->getVal()*Ncmb.getError();
   const double ncmb_errEl_npq = TMath::Sqrt(ncmbEl*(Ncmb.getVal()-ncmbEl)/Ncmb.getVal());
   const double ncmb_errEl_total = TMath::Sqrt(ncmb_errEl*ncmb_errEl+ncmb_errEl_npq*ncmb_errEl_npq);
   const double purityEl = nsigEl/(nsigEl+nrhoEl+ncmbEl);
   const double purity_errEl = nsig_errEl_total/(nsigEl+nrhoEl+ncmbEl);

   RooAbsReal* intSigEl1 = pdf_sig.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   RooAbsReal* intRhoEl1 = pdf_rho->createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   RooAbsReal* intCmbEl1 = pdf_comb.createIntegral(RooArgSet(de,mbc),NormSet(RooArgSet(de,mbc)),Range("Elli"));
   const double nsigEl1 = intSigEl1->getVal()*Nsig.getVal();
   const double nsig_errEl1 = intSigEl1->getVal()*Nsig.getError();
   const double nsig_errEl1_npq = TMath::Sqrt(nsigEl1*(Nsig.getVal()-nsigEl1)/Nsig.getVal());
   const double nsig_errEl1_total = TMath::Sqrt(nsig_errEl1*nsig_errEl1+nsig_errEl1_npq*nsig_errEl1_npq);
   const double nrhoEl1 = intRhoEl1->getVal()*Nrho.getVal();
   const double nrho_errEl1 = intRhoEl1->getVal()*Nrho.getError();
   const double nrho_errEl1_npq = TMath::Sqrt(nrhoEl1*(Nrho.getVal()-nrhoEl1)/Nrho.getVal());
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
  pdf.plotOn(deFrame,Components(*pdf_rho),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf.plotOn(deFrame,Components(pdf_comb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
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
  if(!type) out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  else      out1 << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  pt->AddText(out1.str().c_str());
  out1.str("");
  if(!type) out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  else out1 << "Purity: " << std::fixed << std::setprecision(2) << purityEl*100. << " #pm " << purity_errEl*100;
  pt->AddText(out1.str().c_str());
  pt->Draw();

  TLine *de_line_RIGHT;
  if(!type) de_line_RIGHT = new TLine(de_max,0,de_max,120);
  else      de_line_RIGHT = new TLine(de_max,0,de_max,30);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  if(!type) de_line_LEFT = new TLine(de_min,0,de_min,120);
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
  pdf.plotOn(mbcFrame,Components(*pdf_rho),LineStyle(kDashed),ProjectionRange("deSignal"));
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
  if(!type) out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  else      out1 << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  ptmbc->AddText(out1.str().c_str());
  out1.str("");
  if(!type) out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  else out1 << "Purity: " << std::fixed << std::setprecision(2) << purityEl*100. << " #pm " << purity_errEl*100;
  ptmbc->AddText(out1.str().c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT;
  if(!type) mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,70);
  else      mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,40);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT;
  if(!type) mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,70);
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
  
  if(!type){
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
  cout << "Nsig = " << nsig <<" +- " << nsig_err << " +- " << nsig_err_npq << " (" << nsig_err_total << ")" << endl;
  cout << "Nrho = " << nrho <<" +- " << nrho_err << " +- " << nrho_err_npq << " (" << nrho_err_total << ")" << endl;
  cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << " +- " << ncmb_err_npq << " (" << ncmb_err_total << ")" << endl;
  cout << "Pury = " << purity << " +- " << purity_err << endl;

  cout << "Ellips:" << endl;
  cout << "Nsig = " << nsigEl <<" +- " << nsig_errEl << " +- " << nsig_errEl_npq << " (" << nsig_errEl_total << ")" << endl;
  cout << "Nrho = " << nrhoEl <<" +- " << nrho_errEl << " +- " << nrho_errEl_npq << " (" << nrho_errEl_total << ")" << endl;
  cout << "Ncmb = " << ncmbEl <<" +- " << ncmb_errEl << " +- " << ncmb_errEl_npq << " (" << ncmb_errEl_total << ")" << endl;
  cout << "Pury = " << purityEl << " +- " << purity_errEl << endl;
  
  cout << "Elli:" << endl;
  cout << "Nsig = " << nsigEl1 <<" +- " << nsig_errEl1 << " +- " << nsig_errEl1_npq << " (" << nsig_errEl1_total << ")" << endl;
  cout << "Nrho = " << nrhoEl1 <<" +- " << nrho_errEl1 << " +- " << nrho_errEl1_npq << " (" << nrho_errEl1_total << ")" << endl;
  cout << "Ncmb = " << ncmbEl1 <<" +- " << ncmb_errEl1 << " +- " << ncmb_errEl1_npq << " (" << ncmb_errEl1_total << ")" << endl;
  cout << "Pury = " << purityEl1 << " +- " << purity_errEl1 << endl;
}
