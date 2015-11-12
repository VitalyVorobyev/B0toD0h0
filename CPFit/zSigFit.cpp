using namespace RooFit;

void zSigFit(const int EXP = -1){
  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const double init_scale = 100.;
  const double z_min = -1;
  const double z_max =  1.5;
  const bool true_params = true;

  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("badpi0",5);
  argset.add(b0f);

  RooCategory mode("mode","mode");
  mode.defineType("omega",3);
  argset.add(mode);

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

  RooRealVar z_sig("z_sig","z_sig",z_min,z_max,"mm"); argset.add(z_sig);
  RooRealVar sz_sig("sz_sig","dz_sig",0.01,0.15);     argset.add(sz_sig);
  RooRealVar xi_sig("xi_sig","xi_sig",0.,500);        argset.add(xi_sig);

  RooDataSet omegads("omegads","omegads",omegatree,argset,"(z_sig>0 || z_sig<0) && (sz_sig>0 || sz_sig<0)");
  omegads.Print();

  RooRealVar tau("tau","tau",0.196,0.01,1,"mm");                   // tau.setConstant(kTRUE);
  RooRealVar dzsigOm("dzsigOm","dzsigOm",0,3.*z_min,3.*z_max,"mm");// dzsigOm.setConstant(kTRUE);
  RooRealVar s0sigOm("s0sigOm","s0sigOm",1.,0.1,5.);
  RooRealVar s1sigOm("s1sigOm","s1sigOm",0.01,0.,27.);

  RooFormulaVar sigma("sigma","@0+@1*@2",RooArgSet(s0sigOm,s1sigOm,xi_sig));
  RooGaussModel gsigOm("gsigOm","gsigOm",z_sig,dzsigOm,sigma,sz_sig);
  RooDecay pdfsigOm("pdfsigOm","pdfsigOm",z_sig,tau,gsigOm,RooDecay::SingleSided);

//  RooRealVar s_otlr("s_otlr","s_otlr",40.,1.,100.); s_otlr.setConstant(kTRUE);
//  RooConstVar z0_otlr("z0_otlr","z0_otlr",0.);
//  RooGaussian g_otlr("g_otlr","g_otlr",z_sig,z0_otlr,s_otlr);
//  RooRealVar f_otlr("f_otlr","f_otlr",0.05,0.,1.);

  RooFitResult* fitRes = pdfsigOm.fitTo(omegads,Verbose(),ConditionalObservables(RooArgSet(xi_sig,sz_sig)));
//  fitRes->Print();

///////////
// Plots //
///////////
  RooPlot* zFrameSig = z_sig.frame();

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
  out << "zSigFit_Exp" << EXP << ".png";
  z_cm_sig->Print(out.str().c_str());
  out.str("");
  out << "zSigFit_Exp" << EXP << ".root";
  z_cm_sig->Print(out.str().c_str());

  return;
}

