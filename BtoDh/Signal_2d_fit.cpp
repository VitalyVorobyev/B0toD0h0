#include "cuts.h"
using namespace RooFit;

void Signal_2d_fit(const int _mode = 1, const int _b0f=-1, const bool draw_bins = false){
  // mode 1 -> pi0
  // mode 2 -> eta
  // mode 3 -> omega
  const double cm2ps = 78.48566945838871754705;
  const bool projection_flag = true;
  const bool save_flag       = true;
  int mode,h0mode;
  double de_sig_min,de_sig_max;
  double BDTG_MIN = 0;
  double BDTG_MAX = 1;
  bool gg_flag = true;

  const bool cond_flag = false;//(b0f == 5);
  const bool nsk_mbc   = !cond_flag;
  string label;
  TFile *ifile;
  switch(_mode){
    case 1:
      mode = 1;
      h0mode = 10;
      de_sig_min = de_min;
      de_sig_max = de_max;
      BDTG_MIN = bdtg_cut_pi0;
//      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s5_full.root");
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
      label = string("#pi^{0}");
      break;
    case 2:
      mode = 2;
      h0mode = 10;
      de_sig_min = de_min;
      de_sig_max = de_max;
//      BDTG_MIN = bdtg_cut_etagg;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
      label = string("#eta#rightarrow#gamma#gamma");
      break;
    case 3:
      mode = 2;
      h0mode = 20;
      de_sig_min = de_min_etappp;
      de_sig_max = de_max_etappp;
//      BDTG_MIN = bdtg_cut_etappp;
      gg_flag = false;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
      label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
      break;
    case 4:
      mode = 3;
      h0mode = 20;
      de_sig_min = de_min_omega;
      de_sig_max = de_max_omega;
//      BDTG_MIN = bdtg_cut_omega;
      gg_flag = false;
      ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
      label = string("#omega");
      break;
    default:
      return;
  }

  bool remove_right_CB_flag  = false;
  if(_b0f == 5) remove_right_CB_flag = true;

  TTree *tree = (TTree*)ifile->Get("TEvent");
  RooArgSet argset;

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

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good",1);
  argset.add(good_icpv);

  RooCategory b0f("b0f","b0f");
  if(_b0f == 1 || _b0f == -1) b0f.defineType("signal",1);
  if(_b0f == 1 || _b0f == -1) b0f.defineType("fsr",10);
  if(_b0f == 5 || _b0f == -1){
    b0f.defineType("bad_pi0",5);
    b0f.defineType("bad_pi0_4",4);
  }
  argset.add(b0f);

  RooCategory m_mode("mode","mode");
  m_mode.defineType("sig",mode);
  argset.add(m_mode);

  RooCategory m_h0mode("h0mode","h0mode");
  m_h0mode.defineType("h0sig",h0mode);
  argset.add(m_h0mode);

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin  = -0.3;
  const double deMax  =  0.3;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_sig_min,de_sig_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_sig_min,de_sig_max);

  stringstream out;

  RooRealVar chi2_vtx_d0("chi2_vtx_d0","chi2_vtx_d0",0,50); argset.add(chi2_vtx_d0);

//  RooRealVar phi("phi","phi",0.,-1,1); phi.setConstant(kTRUE);
//  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
//  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  RooRealVar de0_201("de0_201","de0_201",get_de0(mode,h0mode,1),-0.1,0.1); if(cSIG) de0_201.setConstant(kTRUE);
  ////////////
  // de pdf //
  ////////////
  if(gg_flag){
    RooRealVar de0("de0","de0",get_de0(mode,h0mode,_b0f),-0.2,0.1); if(cSIG) de0.setConstant(kTRUE);
    RooRealVar s1("s1","s1",get_s1(mode,h0mode,_b0f),0.,0.5);       if(cSIG) s1.setConstant(kTRUE);
    RooGaussian g1("g1","g1",de,de0,s1);

    RooRealVar deCBl("deCBl","deCBl",get_deCBl(mode,h0mode,_b0f),-0.2,0.1); if(cSIG) deCBl.setConstant(kTRUE);
    RooRealVar sCBl("sCBl","sCBl",get_sCBl(mode,h0mode,_b0f),0.,0.5);       if(cSIG) sCBl.setConstant(kTRUE);
    RooRealVar nl("nl","nl",2.,0.,100.); nl.setConstant(kTRUE);
    RooRealVar alphal("alphal","alphal",get_alphal(mode,h0mode,_b0f),-10.,10.); if(cSIG) alphal.setConstant(kTRUE);

    RooRealVar deCBr("deCBr","deCBr",get_deCBr(mode,h0mode,_b0f),-0.2,0.1);     if(cSIG || remove_right_CB_flag) deCBr.setConstant(kTRUE);
    RooRealVar sCBr("sCBr","sCBr",get_sCBr(mode,h0mode,_b0f),0.,0.5);           if(cSIG || remove_right_CB_flag) sCBr.setConstant(kTRUE);
    RooRealVar nr("nr","nr",2,0.,100.); nr.setConstant(kTRUE);
    RooRealVar alphar("alphar","alphar",get_alphar(mode,h0mode,_b0f),-10.,10.); if(cSIG) alphar.setConstant(kTRUE);

    RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
    RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

    RooRealVar fCBl("fCBl","fCBl",get_fCBl(mode,h0mode,_b0f),0.,1.); if(cSIG) fCBl.setConstant(kTRUE);
    RooRealVar fCBr("fCBr","fCBr",get_fCBr(mode,h0mode,_b0f),0.,1.); if(cSIG || remove_right_CB_flag) fCBr.setConstant(kTRUE);
    if(remove_right_CB_flag) fCBr.setVal(0);

    RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  } else{
    RooRealVar s1_201("s1_201","s1_201",get_s1(mode,h0mode,1),0.,0.5); if(cSIG) s1_201.setConstant(kTRUE);
    RooGaussian g1_201("g1_201","g1_201",de,de0_201,s1_201);

    RooRealVar deCBl_201("deCBl_201","deCBl_201",get_deCBl(mode,h0mode,1),-0.1,0.1);     if(cSIG || _b0f == -1) deCBl_201.setConstant(kTRUE);
    RooRealVar sCBl_201("sCBl_201","sCBl_201",get_sCBl(mode,h0mode,1),0.,0.5);           if(cSIG || _b0f == -1) sCBl_201.setConstant(kTRUE);
    RooRealVar nl_201("nl_201","nl_201",2.,0.,100.); nl_201.setConstant(kTRUE);
    RooRealVar alphal_201("alphal_201","alphal_201",get_alphal(mode,h0mode,1),-10.,10.); if(cSIG || _b0f == -1) alphal_201.setConstant(kTRUE);
    RooRealVar deCBr_201("deCBr_201","deCBr_201",get_deCBr(mode,h0mode,1),-0.1,0.1);     if(cSIG || _b0f == -1) deCBr_201.setConstant(kTRUE);
    RooRealVar sCBr_201("sCBr_201","sCBr_201",get_sCBr(mode,h0mode,1),0.,0.5);           if(cSIG || _b0f == -1) sCBr_201.setConstant(kTRUE);
    RooRealVar nr_201("nr_201","nr_201",2.,0.,100.); nr_201.setConstant(kTRUE);
    RooRealVar alphar_201("alphar_201","alphar_201",get_alphar(mode,h0mode,1),-10.,10.); if(cSIG || _b0f == -1) alphar_201.setConstant(kTRUE);

    RooCBShape CBl_201("CBl_201","CBl_201",de,deCBl_201,sCBl_201,alphal_201,nl_201);
    RooCBShape CBr_201("CBr_201","CBr_201",de,deCBr_201,sCBr_201,alphar_201,nr_201);

    RooRealVar fCBl_201("fCBl_201","fCBl_201",get_fCBl(mode,h0mode,1),0.,1.); if(cSIG || _b0f == -1) fCBl_201.setConstant(kTRUE);
    RooRealVar fCBr_201("fCBr_201","fCBr_201",get_fCBr(mode,h0mode,1),0.,1.); if(cSIG || _b0f == -1) fCBr_201.setConstant(kTRUE);

    if(_b0f != 1) RooAddPdf pdf_de1("pdf_de1","pdf_de1",RooArgList(CBl_201,CBr_201,g1_201),RooArgSet(fCBl_201,fCBr_201));
    else          RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(CBl_201,CBr_201,g1_201),RooArgSet(fCBl_201,fCBr_201));

    RooRealVar de0_205("de0_205","de0_205",get_de0(mode,h0mode,5),-0.2,0.1);  if(cSIG || _b0f == -1) de0_205.setConstant(kTRUE);
    RooRealVar s1_205("s1_205","s1_205",get_s1(mode,h0mode,5),0.,0.5);        if(cSIG || _b0f == -1) s1_205.setConstant(kTRUE);
    RooGaussian g1_205("g1_205","g1_205",de,de0_205,s1_205);

    RooRealVar deCBl_205("deCBl_205","deCBl_205",get_deCBl(mode,h0mode,5),-0.1,0.1); if(cSIG || _b0f == -1) deCBl_205.setConstant(kTRUE);
    RooRealVar sCBl_205("sCBl_205","sCBl_205",get_sCBl(mode,h0mode,5),0.,0.5);       if(cSIG || _b0f == -1) sCBl_205.setConstant(kTRUE);
    RooRealVar nl_205("nl_205","nl_205",2,0.,100.); nl_205.setConstant(kTRUE);
    RooRealVar alphal_205("alphal_205","alphal_205",get_alphal(mode,h0mode,5),-10.,10.); if(cSIG || _b0f == -1) alphal_205.setConstant(kTRUE);
    RooCBShape CBl_205("CBl_205","CBl_205",de,deCBl_205,sCBl_205,alphal_205,nl_205);

    RooRealVar fCBl_205("fCBl_205","fCBl_205",get_fCBl(mode,h0mode,5),0.,1.); if(cSIG || _b0f == -1) fCBl_205.setConstant(kTRUE);

    if(_b0f != 5) RooAddPdf pdf_de5("pdf_de5","pdf_de5",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));
    else          RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));
  }

  /////////////
  // mbc pdf //
  /////////////
  const bool fflag = false;
  if(nsk_mbc){
    if(gg_flag){
      RooRealVar a_s("a_s","a_s",get_a_s(_mode)); a_s.setConstant(kTRUE);
      RooRealVar b_s("b_s","b_s",get_b_s(_mode)); b_s.setConstant(kTRUE);
      RooRealVar c_s("c_s","c_s",get_c_s(_mode),0.0015,0.0035);// c_s.setConstant(kTRUE);
      RooFormulaVar S("S","S","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s,b_s,a_s));

      RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);

      RooRealVar a_mbc0("a_mbc0","a_mbc0",get_a_mbc0(_mode)); a_mbc0.setConstant(kTRUE);
      RooRealVar b_mbc0("b_mbc0","b_mbc0",get_b_mbc0(_mode)); b_mbc0.setConstant(kTRUE);
      RooRealVar c_mbc0("c_mbc0","c_mbc0",get_c_mbc0(_mode),5.27,5.29);//    c_mbc0.setConstant(kTRUE);
      RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));
      RooNovosibirsk pdf_mbc("pdf_mbc","pdf_mbc",mbc,MBC0,S,alpha);
    } else{
      RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
      RooRealVar c0("c0","c0",get_c0(_mode)); c0.setConstant(kTRUE);
      RooRealVar c1("c1","c1",get_c1(_mode)); c1.setConstant(kTRUE);
      RooRealVar c2("c2","c2",get_c2(_mode)); c2.setConstant(kTRUE);
      RooRealVar mbc0("mbc0","mbc0",5.28,5.27,5.29);
      RooFormulaVar MBC("MBC","MBC","@0+@1*TMath::Erf((@2-@3))/@4",RooArgList(mbc0,c0,c1,de,c2));

      RooRealVar a_s1("a_s1","a_s1",get_a_s(_mode),0.15,0.45); a_s1.setConstant(kTRUE);
      RooRealVar b_s1("b_s1","b_s1",get_b_s(_mode),-0.05,0.05); b_s1.setConstant(kTRUE);
      RooRealVar c_s1("c_s1","c_s1",get_c_s(_mode),0.0015,0.0035);
      RooFormulaVar S1("S1","S1","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s1,b_s1,a_s1));
      if(_b0f != 1) RooNovosibirsk pdf_mbc1("pdf_mbc1","pdf_mbc1",mbc,MBC,S1,alpha);
      else          RooNovosibirsk pdf_mbc("pdf_mbc","pdf_mbc",mbc,MBC,S1,alpha);

      RooRealVar a_s5("a_s5","a_s5",get_a5_s(_mode)); a_s5.setConstant(kTRUE);
      RooRealVar b_s5("b_s5","b_s5",get_b5_s(_mode)); b_s5.setConstant(kTRUE);
      RooRealVar c_s5("c_s5","c_s5",get_c5_s(_mode),0.0015,0.0055);
      RooFormulaVar S5("S5","S5","@1+@2*@0+@3*@0*@0",RooArgList(de,c_s5,b_s5,a_s5));

      RooRealVar a_mbc0("a_mbc0","a_mbc0",get_a5_mbc0(_mode)); a_mbc0.setConstant(kTRUE);
      RooRealVar b_mbc0("b_mbc0","b_mbc0",get_b5_mbc0(_mode)); b_mbc0.setConstant(kTRUE);
      RooRealVar c_mbc0("c_mbc0","c_mbc0",get_c5_mbc0(_mode),5.27,5.29);
      RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));
      if(_b0f != 5) RooNovosibirsk pdf_mbc5("pdf_mbc5","pdf_mbc5",mbc,MBC0,S5,alpha);
      else          RooNovosibirsk pdf_mbc("pdf_mbc","pdf_mbc",mbc,MBC0,S5,alpha);
    }
  } else if(cond_flag){
    RooRealVar a_mbc0("a_mbc0","a_mbc0",0.1027);   a_mbc0.setConstant(kTRUE);
    RooRealVar b_mbc0("b_mbc0","b_mbc0",-0.01096); b_mbc0.setConstant(kTRUE);
    RooRealVar c_mbc0("c_mbc0","c_mbc0",5.28,5.27,5.29);    // c_mbc0.setConstant(kTRUE);
    RooFormulaVar MBC0("MBC0","MBC0","@1+@2*@0+@3*@0*@0",RooArgList(de,c_mbc0,b_mbc0,a_mbc0));

    RooRealVar a_sl("a_sl","a_sl", 0.2508);   a_sl.setConstant(kTRUE);
    RooRealVar b_sl("b_sl","b_sl",-0.004734); b_sl.setConstant(kTRUE);
    RooRealVar c_sl("c_sl","c_sl", 0.003002,0.002,0.004);// c_sl.setConstant(kTRUE);
    RooFormulaVar Sl("Sl","Sl","@1+@2*@0+@3*@0*@0",RooArgList(de,c_sl,b_sl,a_sl));

    RooRealVar a_sr("a_sr","a_sr",-0.01417);  a_sr.setConstant(kTRUE);
    RooRealVar b_sr("b_sr","b_sr", 0.002406); b_sr.setConstant(kTRUE);
    RooRealVar c_sr("c_sr","c_sr", 0.002413,0.0015,0.0035);// c_sr.setConstant(kTRUE);
    RooFormulaVar Sr("Sr","Sr","@1+@2*@0+@3*@0*@0",RooArgList(de,c_sr,b_sr,a_sr));
    RooBifurGauss pdf_mbc("pdf_mbc","pdf_mbc",mbc,MBC0,Sl,Sr);
  } else{
    RooRealVar sl("sl","sl",get_sl(mode,h0mode,_b0f),0.,0.5);         if(cSIG && fflag) sl.setConstant(kTRUE);
    RooRealVar sr("sr","sr",get_sr(mode,h0mode,_b0f),0.00000001,0.5); if(cSIG && fflag) sr.setConstant(kTRUE);
    RooRealVar sll("sll","sll",get_sll(mode,h0mode,_b0f),0.,0.5);     if(cSIG && fflag) sll.setConstant(kTRUE);
    RooRealVar srr("srr","srr",get_srr(mode,h0mode,_b0f),0.0005,0.5); if(cSIG && fflag) srr.setConstant(kTRUE);

    RooBifurGauss bg("bg","bg",mbc,mbc0,sl,sr);
    RooRealVar mbc00("mbc00","mbc00",get_mbc00(mode,h0mode,_b0f),5.26,5.30); if(cSIG) mbc00.setConstant(kTRUE);
    if(_b0f != 5 || mode == 3) RooBifurGauss bgg("bgg","bgg",mbc,mbc00,sll,srr);
    else                       RooBifurGauss bgg("bgg","bgg",mbc,mbc0,sll,srr);

    RooRealVar fmbc("fmbc","fmbc",get_fmbc(mode,h0mode,_b0f),0.,1.); if(cSIG && fflag) fmbc.setConstant(kTRUE);
    RooAddPdf pdf_mbc("pdf_mbc","pdf_mbc",RooArgList(bg,bgg),RooArgSet(fmbc));
  }

  /////////
  // pdf //
  /////////
  if(gg_flag || _b0f != -1){
    RooProdPdf pdf("pdf","pdf",pdf_de,Conditional(pdf_mbc,mbc));
  } else{
    RooRealVar f_201("f_201","f_201",get_f201(mode,h0mode),0.,1.);// if(cSIG) f_201.setConstant(kTRUE);
    RooProdPdf pdf1("pdf1","pdf1",pdf_de1,Conditional(pdf_mbc1,mbc));
    RooProdPdf pdf5("pdf5","pdf5",pdf_de5,Conditional(pdf_mbc5,mbc));
    RooAddPdf  pdf("pdf","pdf",RooArgList(pdf1,pdf5),RooArgSet(f_201));
  }

  pdf.fitTo(ds,Verbose(),Timer(true));
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  if(projection_flag){
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
    pdf.plotOn(deFrame,LineWidth(2),LineColor(kGreen));
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
    pdf.plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  } else{
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(deFrame,LineWidth(2));
  }
  ds.statOn(deFrame,Layout(0.55,0.98,0.9));

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
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  TPaveText *pt = new TPaveText(0.6,0.55,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(out.str().c_str());
  pt->AddText(label.c_str());
  pt->Draw();

  if(projection_flag){
  TLine *de_line_RIGHT = new TLine(de_sig_max,0,de_sig_max,250);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_sig_min,0,de_sig_min,250);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();
  }

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
  out.str("");
  out << "../Note/pics/de_sig_m" << _mode << "_b0f" << _b0f;
  out << ".eps";
  if(save_flag) cm->Print(out.str().c_str());

  // mbc //
  RooPlot* mbcFrame = mbc.frame();
  if(projection_flag){
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
    pdf.plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
    pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  } else{
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(mbcFrame,LineWidth(2));
  }

  ds.statOn(mbcFrame,Layout(0.2,0.68,0.9));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc.frame(Title("M_{bc} pull distribution"));
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

  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.3,0.55,0.68,0.7,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(label.c_str());
  ptmbc->Draw();

  if(projection_flag){
  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,800);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,800);
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();
  }

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
  out.str("");
  out << "../Note/pics/mbc_sig_m" << _mode << "_b0f" << _b0f;
  out << ".eps";
  if(save_flag) cmmbc->Print(out.str().c_str());

  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.3,0.3),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");
  cmmbc->Update();
  out.str("");
  out << "../Note/pics/2d_sig_mode" << mode << "_h0mode" << h0mode << "_b0f" << _b0f;
//  if(projection_flag) out << "_wproj";
//  if(deMin>-0.2) out << "015";
  out << ".eps";
  if(save_flag) hhc->Print(out.str().c_str());

  TLine* l1 = new TLine(de_sig_min,mbc_min,de_sig_max,mbc_min);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->SetLineWidth(2);
  TLine* l2 = new TLine(de_sig_min,mbc_max,de_sig_max,mbc_max);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(1);
  l2->SetLineWidth(2);
  TLine* l3 = new TLine(de_sig_min,mbc_min,de_sig_min,mbc_max);
  l3->SetLineColor(kRed);
  l3->SetLineStyle(1);
  l3->SetLineWidth(2);
  TLine* l4 = new TLine(de_sig_max,mbc_min,de_sig_max,mbc_max);
  l4->SetLineColor(kRed);
  l4->SetLineStyle(1);
  l4->SetLineWidth(2);

  TCanvas* ellican = new TCanvas("scatterplot","Scatter Plot Mbc dE",400,400);
  ellican->cd();
  out.str("");
  out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX;
  out << " && mbc>5.265 && de>" << deMin;
  out << " && mode == " << mode << " && h0mode == " << h0mode;
  if(_b0f == 1)      out << " && (b0f == 1 || b0f == 10)";
  else if(_b0f == 5) out << " && b0f != 0 && b0f != -1 && b0f != 1 && b0f != 10";
  else               out << " && b0f > 0";
  cout << "Cuts: " << out.str() << endl;
  tree->Draw("mbc:de",out.str().c_str());
  l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
  out.str("");
  out << "../Reports/pics/scatplot_sig_m" << _mode << "_b0f" << _b0f;
  if(deMin>-0.2) out << "_015" << ".png";
  if(save_flag) ellican->Print(out.str().c_str());

  cout << mode << " " << h0mode << " " << _b0f << endl;

  if(!draw_bins) return;
  for(int i=0; i<2; i++){
    for(int j=0; j<16; j++){
      int Flv = flv(i);
      if(Flv == 0) Flv = -1;
      const int Bin = bin(j);
      out.str("");
      out << "flv_mc == " << Flv << " && bin == " << Bin;
      RooDataSet* ds0 = (RooDataSet*)ds.reduce(RooArgSet(de,mbc),out.str().c_str());

      // de //
      RooPlot* deFrame0 = de.frame();
      ds0->plotOn(deFrame0,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
      pdf.plotOn(deFrame0,LineWidth(2),ProjectionRange("mbcSignal"));
      ds0->statOn(deFrame0,Layout(0.55,0.98,0.9));

      RooHist* hdepull0 = deFrame0->pullHist();
      out.str("");
      out << "#Delta E pull distribution, flv: " << Flv << ", bin: " << Bin;
      RooPlot* dePull0 = de.frame(Title(out.str().c_str()));
      dePull0->addPlotable(hdepull0,"P");
      dePull0->GetYaxis()->SetRangeUser(-5,5);

      out.str("");
      out << "#Delta E, Signal, flv: " << Flv << ", bin: " << Bin;
      TCanvas* cm0 = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
      cm0->cd();

      TPad *pad30 = new TPad("pad30","pad30",0.01,0.20,0.99,0.99);
      TPad *pad40 = new TPad("pad40","pad40",0.01,0.01,0.99,0.20);
      pad30->Draw();
      pad40->Draw();

      pad30->cd();
      pad30->SetLeftMargin(0.15);
      pad30->SetFillColor(0);

      deFrame0->GetXaxis()->SetTitleSize(0.05);
      deFrame0->GetXaxis()->SetTitleOffset(0.85);
      deFrame0->GetXaxis()->SetLabelSize(0.04);
      deFrame0->GetYaxis()->SetTitleOffset(1.6);
      deFrame0->Draw();

      stringstream out;
      out.str("");
      out << "#chi^{2}/n.d.f = " << deFrame0->chiSquare();
      TPaveText *pt0 = new TPaveText(0.6,0.55,0.98,0.7,"brNDC");
      pt0->SetFillColor(0);
      pt0->SetTextAlign(12);
      pt0->AddText(out.str().c_str());
      pt0->AddText(label.c_str());
      pt0->Draw();

      TLine *de_line_RIGHT0 = new TLine(de_sig_max,0,de_sig_max,250);
      de_line_RIGHT0->SetLineColor(kRed);
      de_line_RIGHT0->SetLineStyle(1);
      de_line_RIGHT0->SetLineWidth((Width_t)2.);
      de_line_RIGHT0->Draw();
      TLine *de_line_LEFT0 = new TLine(de_sig_min,0,de_sig_min,250);
      de_line_LEFT0->SetLineColor(kRed);
      de_line_LEFT0->SetLineStyle(1);
      de_line_LEFT0->SetLineWidth((Width_t)2.);
      de_line_LEFT0->Draw();

      pad40->cd(); pad40->SetLeftMargin(0.15); pad40->SetFillColor(0);
      dePull0->SetMarkerSize(0.05); dePull0->Draw();
      TLine *de_lineUP0 = new TLine(deMin,3,deMax,3);
      de_lineUP0->SetLineColor(kBlue);
      de_lineUP0->SetLineStyle(2);
      de_lineUP0->Draw();
      TLine *de_line0 = new TLine(deMin,0,deMax,0);
      de_line0->SetLineColor(kBlue);
      de_line0->SetLineStyle(1);
      de_line0->SetLineWidth((Width_t)2.);
      de_line0->Draw();
      TLine *de_lineDOWN0 = new TLine(deMin,-3,deMax,-3);
      de_lineDOWN0->SetLineColor(kBlue);
      de_lineDOWN0->SetLineStyle(2);
      de_lineDOWN0->Draw();

      cm0->Update();
      out.str("");
      out << "../Note/pics/de_sig_m" << _mode << "_b0f" << _b0f << "_flv" << Flv << "_bin" << Bin;
      out << ".eps";
      if(save_flag) cm0->Print(out.str().c_str());

      // mbc //
      RooPlot* mbcFrame0 = mbc.frame();
      ds0->plotOn(mbcFrame0,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
      pdf.plotOn(mbcFrame0,LineWidth(2),ProjectionRange("deSignal"));

      ds0->statOn(mbcFrame0,Layout(0.2,0.68,0.9));

      out.str("");
      out << "M_{bc} pull distribution, flv: " << Flv << ", bin: " << Bin;
      RooHist* hmbcpull0 = mbcFrame0->pullHist();
      RooPlot* mbcPull0 = mbc.frame(Title(out.str().c_str()));
      mbcPull0->addPlotable(hmbcpull0,"P");
      mbcPull0->GetYaxis()->SetRangeUser(-5,5);

      out.str("");
      out << "M_{bc}, Signal, flv: " << Flv << ", bin: " << Bin;
      TCanvas* cmmbc0 = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
      cmmbc0->cd();

      TPad *pad10 = new TPad("pad10","pad10",0.01,0.20,0.99,0.99);
      TPad *pad20 = new TPad("pad20","pad20",0.01,0.01,0.99,0.20);
      pad10->Draw();
      pad20->Draw();

      pad10->cd();
      pad10->SetLeftMargin(0.15);
      pad10->SetFillColor(0);

      mbcFrame0->GetXaxis()->SetTitleSize(0.05);
      mbcFrame0->GetXaxis()->SetTitleOffset(0.85);
      mbcFrame0->GetXaxis()->SetLabelSize(0.04);
      mbcFrame0->GetYaxis()->SetTitleOffset(1.6);
      mbcFrame0->Draw();

      out.str("");
      out << "#chi^{2}/n.d.f = " << mbcFrame0->chiSquare();
      TPaveText *ptmbc0 = new TPaveText(0.3,0.55,0.68,0.7,"brNDC");
      ptmbc0->SetFillColor(0);
      ptmbc0->SetTextAlign(12);
      ptmbc0->AddText(out.str().c_str());
      ptmbc0->AddText(label.c_str());
      ptmbc0->Draw();

      TLine *mbc_line_RIGHT0 = new TLine(mbc_max,0,mbc_max,800);
      mbc_line_RIGHT0->SetLineColor(kRed);
      mbc_line_RIGHT0->SetLineStyle(1);
      mbc_line_RIGHT0->SetLineWidth((Width_t)2.);
      mbc_line_RIGHT0->Draw();
      TLine *mbc_line_LEFT0 = new TLine(mbc_min,0,mbc_min,800);
      mbc_line_LEFT0->SetLineColor(kRed);
      mbc_line_LEFT0->SetLineStyle(1);
      mbc_line_LEFT0->SetLineWidth((Width_t)2.);
      mbc_line_LEFT0->Draw();

      pad20->cd();
      pad20->SetLeftMargin(0.15);
      pad20->SetFillColor(0);
      mbcPull0->SetMarkerSize(0.05);
      mbcPull0->Draw();
      TLine *mbc_lineUP0 = new TLine(mbcMin,3,mbcMax,3);
      mbc_lineUP0->SetLineColor(kBlue);
      mbc_lineUP0->SetLineStyle(2);
      mbc_lineUP0->Draw();
      TLine *mbc_line0 = new TLine(mbcMin,0,mbcMax,0);
      mbc_line0->SetLineColor(kBlue);
      mbc_line0->SetLineStyle(1);
      mbc_line0->SetLineWidth((Width_t)2.);
      mbc_line0->Draw();
      TLine *mbc_lineDOWN0 = new TLine(mbcMin,-3,mbcMax,-3);
      mbc_lineDOWN0->SetLineColor(kBlue);
      mbc_lineDOWN0->SetLineStyle(2);
      mbc_lineDOWN0->Draw();

      cmmbc0->Update();
      out.str("");
      out << "../Note/pics/mbc_sig_m" << _mode << "_b0f" << _b0f << "_flv" << Flv << "_bin" << Bin;
      out << ".eps";
      if(save_flag) cmmbc->Print(out.str().c_str());
    }
  }
}
