using namespace RooFit;

void zAscFit1(const int EXP = -1){
  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const double init_scale = 100.;
  const double z_min = -1;
  const double z_max =  1.5;

  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("badpi0",5);
  argset.add(b0f);

  RooCategory exp("exp","exp");
  if(EXP>0){
    exp.defineType("expNo",EXP);
  }
  if(EXP == -2 || EXP == -3){
    exp.defineType("7",7);
    exp.defineType("9",9);
    exp.defineType("11",11);
    exp.defineType("13",13);
    exp.defineType("15",15);
    exp.defineType("17",17);
    exp.defineType("21",21);
    exp.defineType("23",23);
    exp.defineType("25",25);
    exp.defineType("27",27);
  }
  if(EXP == -1 || EXP == -3){
    exp.defineType("31",31);
    exp.defineType("33",33);
    exp.defineType("35",35);
    exp.defineType("37",37);
    exp.defineType("39",39);
    exp.defineType("41",41);
    exp.defineType("43",43);
    exp.defineType("45",45);
    exp.defineType("47",47);
    exp.defineType("49",49);
    exp.defineType("51",51);
    exp.defineType("55",55);
    exp.defineType("61",61);
    exp.defineType("63",63);
    exp.defineType("65",65);
  }
  argset.add(exp); 

  RooCategory nptag("nptag","nptag");
  nptag.defineType("nonp",0);
  argset.add(nptag);

  RooRealVar z_asc("z_asc","z_asc",z_min,z_max,"mm"); argset.add(z_asc);
  RooRealVar sz_asc("sz_asc","dz_asc",0.01,0.15);     argset.add(sz_asc);
  RooRealVar chisq_z_asc("chisq_z_asc","chisq_z_asc",0.,200); argset.add(chisq_z_asc);

//  RooDataSet pids("pids","pids",pitree,argset);
  RooDataSet omegads("omegads","omegads",omegatree,argset,"(z_asc>0 || z_asc<0) && (sz_asc>0 || sz_asc<0)");
//  pids.Print();
  omegads.Print();

  RooRealVar tau("tau","tau",0.196,0.01,1,"mm"); tau.setConstant(kTRUE);
  RooRealVar dz01sigOm("dz01sigOm","dz01sigOm",0,3.*z_min,3.*z_max,"mm");// dz01sigOm.setConstant(kTRUE);
//  RooRealVar dz02sigOm("dz02sigOm","dz02sigOm",0,z_min,z_max,"mm"); dz02sigOm.setConstant(kTRUE);
  RooRealVar s01sigOm("s01sigOm","s01sigOm",1,0.1,5.);
  RooRealVar s02sigOm("s02sigOm","s02sigOm",4,0.1,27.);
  RooRealVar fsigOm("fsigOm","fsigOm",0.2,0.,1.);
  RooGaussModel g1sigOm("g1sigOm","g1sigOm",z_asc,dz01sigOm,s01sigOm,sz_asc);
  RooGaussModel g2sigOm("g2sigOm","g2sigOm",z_asc,dz01sigOm,s02sigOm,sz_asc);
  RooDecay dec1SigOm("dec1SigOm","dec1SigOm",z_asc,tau,g1sigOm,RooDecay::SingleSided);
  RooDecay dec2SigOm("dec2SigOm","dec2SigOm",z_asc,tau,g2sigOm,RooDecay::SingleSided);

  RooRealVar s_otlr("s_otlr","s_otlr",40.,1.,100.); s_otlr.setConstant(kTRUE);
  RooConstVar z0_otlr("z0_otlr","z0_otlr",0.);
  RooGaussian g_otlr("g_otlr","g_otlr",z_asc,z0_otlr,s_otlr);
  RooRealVar f_otlr("f_otlr","f_otlr",0.05,0.,1.);

  RooAddPdf pdfsigOm("pdfsigOm","pdfsigOm",RooArgList(g_otlr,dec2SigOm,dec1SigOm),RooArgSet(f_otlr,fsigOm));
  if(EXP!=55 && EXP>0){
    fsigOm.setVal(1.);
    fsigOm.setConstant(kTRUE);
    s02sigOm.setConstant(kTRUE);
  }
//  RooDecay pdfsigOm("pdfsigOm","pdpsigOm",z_sig,tau,g1sigOm,RooDecay::SingleSided);

//  pdfsigPi.fitTo(pids,ConditionalObservables(szSigPi));
  RooFitResult* fitRes = pdfsigOm.fitTo(omegads,Verbose(),ConditionalObservables(sz_asc));
//  fitRes->Print();

///////////
// Plots //
///////////
  RooPlot* zFrameSig = z_asc.frame();

  omegads.plotOn(zFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdfsigOm.plotOn(zFrameSig,LineColor(kRed),ProjWData(omegads));
  TCanvas* z_cm_sig = new TCanvas("z_{sig}, cm","z_{sig}, mm",700,500);
  z_cm_sig->cd();
  omegads.statOn(zFrameSig,Layout(0.6,0.98,0.9),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  zFrameSig->GetXaxis()->SetTitleOffset(0.75);
  zFrameSig->GetXaxis()->SetLabelSize(0.05);
  zFrameSig->GetXaxis()->SetTitleSize(0.06);
  zFrameSig->GetYaxis()->SetTitleOffset(1.6);
  zFrameSig->Draw();
  z_cm_sig->Update();
  stringstream out;
  out.str("");
  out << "zAscFit_Exp" << EXP << ".png";
  z_cm_sig->Print(out.str().c_str());
  out.str("");
  out << "zAscFit_Exp" << EXP << ".root";
  z_cm_sig->Print(out.str().c_str());

  return;
}

