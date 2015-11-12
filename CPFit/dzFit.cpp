using namespace RooFit;

void dzFit(const int _mode = 1, const int _h0mode = 10,const int EXP = -1,int fit_mode = 1){
  TFile *omegafile;
  switch(_mode){
  case 1:
    omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s4_full.root");
    break;
  case 2:
    omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
    break;
  case 3:
    omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
    break;
  }
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");
  if(_mode != 3) fit_mode = 1;

  string dzname,szname,xiname;
  switch(fit_mode){
  case 1:// Full fit
    dzname = string("dz_mc_sig");
    szname = string("sz_sig");
    xiname = string("chisq_z_sig");
    break;
  case 2:// pi+pi- fit
    dzname = string("dz_mc_sig_pipi");
    szname = string("sz_sig_pipi");
    xiname = string("xi_sig_pipi");
    break;
  case 3:// D0 fit
    dzname = string("dz_mc_sig_d0");
    szname = string("sz_sig_d0");
    xiname = string("chisq_sig_d0");
    break;
  }

  const double dz_min = -1;
  const double dz_max =  1;
  const bool true_params = true;

  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("badpi0",5);
  argset.add(b0f);

  RooCategory mode("mode","mode");
  mode.defineType("omega",_mode);
  argset.add(mode);

  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("h0mode",_h0mode);
  argset.add(h0mode);

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

  RooRealVar dz_mc_sig(dzname.c_str(),"#deltaz_{mc}^{sig}",dz_min,dz_max,"mm"); argset.add(dz_mc_sig);
  RooRealVar sz_sig(szname.c_str(),szname.c_str(),0.0001,0.5);   argset.add(sz_sig);
  RooRealVar xi_sig(xiname.c_str(),xiname.c_str(),0.,500);       argset.add(xi_sig);

//  cout << "Test 1" << endl;
  stringstream out;
  out << dzname << ">0 || 0>" << dzname << " && 0<" << szname << " && 0<" << xiname;
  RooDataSet omegads("omegads","omegads",omegatree,argset,out.str().c_str());
//  cout << "Test 2" << endl;
  omegads.Print();

  if(fit_mode != 3){
  RooRealVar dzsigOm("dzsigOm","dzsigOm",0,-0.2,0.2,"mm");// dzsigOm.setConstant(kTRUE);
  RooRealVar s0sigOm("s0sigOm","s0sigOm",1.76,0.01,5.);
  RooRealVar s1sigOm("s1sigOm","s1sigOm",0.038,0.,0.5);

  RooFormulaVar sigma("sigma","(@0+@1*@2)*@3",RooArgSet(s0sigOm,s1sigOm,xi_sig,sz_sig));
  RooGaussModel gsigOm("gsigOm","gsigOm",dz_mc_sig,dzsigOm,sigma);

  RooRealVar s_otlr("s_otlr","s_otlr",40.,0.01,100.);// s_otlr.setConstant(kTRUE);
  RooRealVar z0_otlr("z0_otlr","z0_otlr",0.,-0.5,0.5);
//  RooConstVar z0_otlr("z0_otlr","z0_otlr",0.);
//  RooGaussian g_otlr("g_otlr","g_otlr",dz_mc_sig,z0_otlr,s_otlr);
  RooGaussModel g_otlr("g_otlr","g_otlr",dz_mc_sig,z0_otlr,s_otlr);
  RooRealVar f_otlr("f_otlr","f_otlr",0.05,0.,1.);// f_otlr.setConstant(kTRUE);

  RooAddPdf pdf("pdf","pdf",RooArgList(g_otlr,gsigOm),RooArgSet(f_otlr));

  RooFitResult* fitRes = pdf.fitTo(omegads,Verbose(),ConditionalObservables(RooArgSet(xi_sig,sz_sig)));
  } else{
  RooRealVar dzsigOm("dzsigOm","dzsigOm",0,3.*dz_min,3.*dz_max,"mm"); dzsigOm.setConstant(kTRUE);
  RooRealVar sMainSigOm("sMainSigOm","sMainSigOm",1.76,0.01,5.);
  RooRealVar sTailSigOm("sTailSigOm","sTailSigOm",5,0.,27.);
  RooGaussModel gMainSigOm("gMainSigOm","gMainSigOm",dz_mc_sig,dzsigOm,sMainSigOm,sz_sig);
  RooGaussModel gTailSigOm("gTailSigOm","gTailSigOm",dz_mc_sig,dzsigOm,sTailSigOm,sz_sig);
  RooRealVar f_tail("f_tail","f_tail",0.1,0.,1.);
  RooAddPdf pdf("pdf","pdf",RooArgList(gTailSigOm,gMainSigOm),RooArgSet(f_tail));

  RooFitResult* fitRes = pdf.fitTo(omegads,Verbose(),ConditionalObservables(RooArgSet(sz_sig)));
  }

///////////
// Plots //
///////////
  RooPlot* dzFrameSig = dz_mc_sig.frame();
  omegads.plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
  if(fit_mode != 3) pdf.plotOn(dzFrameSig,Components(g_otlr),LineStyle(kDashed),LineColor(kBlue),ProjWData(omegads));
  else pdf.plotOn(dzFrameSig,Components(gTailSigOm),LineStyle(kDashed),LineColor(kBlue),ProjWData(omegads));
  pdf.plotOn(dzFrameSig,LineColor(kBlue),ProjWData(omegads));

  RooHist* hdzpull = dzFrameSig->pullHist();
  RooPlot* dzPull = dz_mc_sig.frame(Title("#deltaz pull distribution"));
  dzPull->addPlotable(hdzpull,"P");
  dzPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* dz_cm_sig = new TCanvas("#deltaz_{sig}, cm","#deltaz_{sig}, mm",700,800);
  dz_cm_sig->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd(); pad1->SetLeftMargin(0.15); pad1->SetFillColor(0);
  omegads.statOn(dzFrameSig,Label("#omega"),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetLogy();
  dzFrameSig->GetXaxis()->SetTitleOffset(0.75);
  dzFrameSig->GetXaxis()->SetLabelSize(0.05);
  dzFrameSig->GetXaxis()->SetTitleSize(0.06);
  dzFrameSig->GetYaxis()->SetTitleOffset(1.6);
  dzFrameSig->Draw();

  out.str("");
  TPaveText *pt = new TPaveText(0.63,0.65,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << dzFrameSig->chiSquare();
  pt->AddText(out.str().c_str());
  pt->Draw();

  pad2->cd(); pad2->SetLeftMargin(0.15); pad2->SetFillColor(0);
  dzPull->SetMarkerSize(0.05); dzPull->Draw();
  TLine *dz_lineUP = new TLine(dz_min,3,dz_max,3);
  dz_lineUP->SetLineColor(kBlue);
  dz_lineUP->SetLineStyle(2);
  dz_lineUP->Draw();
  TLine *dz_line = new TLine(dz_min,0,dz_max,0);
  dz_line->SetLineColor(kBlue);
  dz_line->SetLineStyle(1);
  dz_line->SetLineWidth((Width_t)2.);
  dz_line->Draw();
  TLine *dz_lineDOWN = new TLine(dz_min,-3,dz_max,-3);
  dz_lineDOWN->SetLineColor(kBlue);
  dz_lineDOWN->SetLineStyle(2);
  dz_lineDOWN->Draw();

  dz_cm_sig->Update();
//  dz_cm_sig->Print("dz_sig_fit.png");
//  dz_cm_sig->Print("dz_sig_fit.root");

}

