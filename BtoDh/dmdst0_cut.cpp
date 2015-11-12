using namespace RooFit;

void dmdst0_cut(const bool clean = false){
  TChain* tree = new TChain("TEvent","TEvent");
  tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");

  RooArgSet argset;

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good_icpv",1);
  argset.add(good_icpv);
  RooCategory mode("mode","mode");
  mode.defineType("pi0",10);
  mode.defineType("eta",20);
  argset.add(mode);
  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  h0mode.defineType("ppp",20);
  argset.add(h0mode);
  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  if(!clean) b0f.defineType("wrph",5);
  if(!clean) b0f.defineType("fsr",10);
  argset.add(b0f);

  RooRealVar chi2_mass_d0("chi2_mass_d0","chi2_mass_d0",0.,50.); argset.add(chi2_mass_d0);

  RooRealVar mbc("mbc","M_{bc}",5.24,5.291,"GeV");
  mbc.setRange("mbc_fit",5.28-0.01,5.28+0.01);
  mbc.setRange("mbc_plot",5.24,5.291);
  argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.15,0.3,"GeV");
  de.setRange("de_fit",-0.03,0.03);
  de.setRange("de_plot",-0.15,0.3);
  argset.add(de);
  RooRealVar dmdst0("dmdst0","m(D^{*0})-m(D^{0})",0.135,0.142+0.012,"GeV");
  dmdst0.setRange("dm_fit",0.142-0.003,0.142+0.003);
  dmdst0.setRange("dm_plot",0.135,0.142+0.012);
  argset.add(dmdst0);
  RooRealVar md("md_raw","m(D^{0})",1.85155,1.87833,"GeV");
  argset.add(md);

  // Pdf definition //
  RooRealVar mean("mean","mean",0.142,-0.1,10.);
  RooRealVar sl("sl","sl",0.1,0.,10.);
  RooRealVar sr("sr","sr",0.1,0.,10.);
  RooBifurGauss pdf_mbc("pdf_mbc","pdf_mbc",mbc,mean,sl,sr);
  RooBifurGauss pdf_de("pdf_de","pdf_de",de,mean,sl,sr);
  RooBifurGauss pdf_dm("pdf_dm","pdf_dm",dmdst0,mean,sl,sr);

  // * dm(D*) * //
  RooRealVar  s1("s1","s1",0.002,0.,0.1);
  RooGaussian g1("g1","g1",dmdst0,mean,s1);

  RooRealVar alphal("alphal","alphal",1., 0.,3.);// alphal.setConstant(kTRUE);
  RooRealVar nl("nl","nl",2,0.,100.); nl.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",-1.,-3.,0.);// alphar.setConstant(kTRUE);
  RooRealVar nr("nr","nr",2,0.,100.); nr.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",dmdst0,mean,sl,alphal,nl);
  RooCBShape CBr("CBr","CBr",dmdst0,mean,sr,alphar,nr);

  RooRealVar a0("a0","a0",0.,-10.,100.); a0.setConstant(kTRUE);
  RooRealVar a1("a1","a1",-1000.,-10000.,1.);
  RooRealVar edge("edge","edge",0.135,0.133,0.137); edge.setConstant(kTRUE);
  RooRealVar pow("pow","pow",0.3,0.1.,1.);          pow.setConstant(kTRUE);
  RooGenericPdf pdf_dm_5("pdf_dm_5","pdf_dm_5","pow(@0-@3,@1)*(1+@2*(@0-@3)+@4*(@0-@3)*(@0-@3)) > 0 ? pow(@0-@3,@1)*(1+@2*(@0-@3)+@4*(@0-@3)*(@0-@3)) : 0.0001",RooArgList(dmdst0,pow,a0,edge,a1));

  RooRealVar f5("f5","f5",0.3,0.,1.);
  RooRealVar fCBl("fCBl","fCBl",0.1,0.,1.);
  RooRealVar fCBr("fCBr","fCBr",0.1,0.,1.);
  RooAddPdf pdf_m_dm("pdf_m_dm","pdf_m_dm",RooArgSet(CBr,pdf_dm_5,g1),RooArgList(fCBr,f5));

  // * Delta E * //
  RooRealVar deCBl("deCBl","deCBl",0.,-0.2,0.1);
  RooRealVar deCBr("deCBr","deCBr",0.,-0.2,0.1);
  RooCBShape CBlde("CBlde","CBlde",de,deCBl,sl,alphal,nl);
//  RooCBShape CBrde("CBrde","CBrde",de,deCBr,sr,alphar,nr);
  RooGaussian g1de("g1de","g1de",de,mean,s1);

//  RooAddPdf pdf_m_de("pdf_m_de","pdf_m_de",RooArgList(CBlde,CBrde,g1de),RooArgSet(fCBl,fCBr));
  RooAddPdf pdf_m_de("pdf_m_de","pdf_m_de",RooArgList(CBlde,g1de),RooArgSet(fCBl));

  // * Mbc * //
  RooRealVar alpha_mbc("alpha_mbc","alpha_mbc",0.139,0.01,2.);// alpha.setConstant(kTRUE);
  RooNovosibirsk pdf_m_mbc("pdf_m_mbc","pdf_m_mbc",mbc,mean,s1,alpha_mbc);

  // ** D* pi0 ** //
  const string dm_m10h0m10_cut = string("abs(de)<0.1 && mbc>5.271 && mode == 10 && h0mode == 10");
  RooDataSet dm_m10h0m10_ds("dm_m10h0m10_ds","dm_m10h0m10_ds",tree,argset,dm_m10h0m10_cut.c_str());

  mean.setVal(1.42193e-01);
  sl.setVal(1.19097e-03);
  sr.setVal(1.15833e-03);

  pdf_dm.fitTo(dm_m10h0m10_ds,Timer(true),Range("dm_fit"));
  const double m_dm_m10h0m10_min = mean.getVal() - 3.*sl.getVal();
  const double m_dm_m10h0m10_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    a1.setVal(-5.75403e+02); a1.setConstant(kFALSE);
    alphar.setVal(-1.09037e+00);
    f5.setVal(2.86143e-01);  f5.setConstant(kFALSE);
    fCBr.setVal(3.13193e-01);
    mean.setVal(1.42166e-01);
    s1.setVal(6.79274e-04);
    sr.setVal(1.69249e-03);
  } else{
    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
    alphar.setVal(-1.66702e+00);
    fCBr.setVal(4.21899e-01);
    mean.setVal(1.42178e-01);
    s1.setVal(6.48126e-04);
    sr.setVal(1.53696e-03);
  }

  RooFitResult* r = pdf_m_dm.fitTo(dm_m10h0m10_ds,Timer(true),Range("dm_plot"),Save());
  dmdst0.setRange("dmdst0_sig_m10h0m10",m_dm_m10h0m10_min,m_dm_m10h0m10_max);
  const double int_dmdst0_m10h0m10_sig  = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dmdst0_sig_m10h0m10"))->getVal();
  const double int_dmdst0_m10h0m10_plot = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dm_plot"))->getVal();

  RooPlot* m_dmdst0_m10h0m10Frame = dmdst0.frame(Range("dm_plot"),Title("Signal range of m(D*^{0}) - m(D^{0}) for h^{0}=#pi^{0}"));
  dm_m10h0m10_ds.plotOn(m_dmdst0_m10h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_dm.plotOn(m_dmdst0_m10h0m10Frame,LineWidth(2));
  if(!clean) pdf_m_dm.plotOn(m_dmdst0_m10h0m10Frame,Components(pdf_dm_5),LineWidth(2));

  TCanvas* m_dmdst0_m10h0m10_cm = new TCanvas("m_dmdst0_m10h0m10_cm","m_dmdst0_m10h0m10_cm",600,400);
  m_dmdst0_m10h0m10_cm->cd();
  dm_m10h0m10_ds.statOn(m_dmdst0_m10h0m10Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_dmdst0_m10h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_dmdst0_m10h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_dmdst0_m10h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_dmdst0_m10h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_dmdst0_m10h0m10Frame->Draw();
  TLine *mdmdst0_m10h0m10lineLEFT = new TLine(m_dm_m10h0m10_min,0,m_dm_m10h0m10_min,320);
  mdmdst0_m10h0m10lineLEFT->SetLineColor(kRed);
  mdmdst0_m10h0m10lineLEFT->SetLineWidth(2);
  mdmdst0_m10h0m10lineLEFT->SetLineStyle(1);
  mdmdst0_m10h0m10lineLEFT->Draw();
  TLine *mdmdst0_m10h0m10lineRIGHT = new TLine(m_dm_m10h0m10_max,0,m_dm_m10h0m10_max,320);
  mdmdst0_m10h0m10lineRIGHT->SetLineColor(kRed);
  mdmdst0_m10h0m10lineRIGHT->SetLineWidth(2);
  mdmdst0_m10h0m10lineRIGHT->SetLineStyle(1);
  mdmdst0_m10h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_dmdst0_m10h0m10_sig/int_dmdst0_m10h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_dmdst0_m10h0m10_cm->Update();
  m_dmdst0_m10h0m10_cm->Print("pics/m_dmdst0_m10h0m10_cut.eps");
  m_dmdst0_m10h0m10_cm->Print("pics/m_dmdst0_m10h0m10_cut.root");
  // //

  // ** D* eta(->gg) ** //
  const string dm_m20h0m10_cut = string("abs(de)<0.1 && mbc>5.271 && mode == 20 && h0mode == 10");
  RooDataSet dm_m20h0m10_ds("dm_m20h0m10_ds","dm_m20h0m10_ds",tree,argset,dm_m20h0m10_cut.c_str());

  mean.setVal(1.42193e-01);
  sl.setVal(1.19097e-03);
  sr.setVal(1.15833e-03);

  pdf_dm.fitTo(dm_m20h0m10_ds,Timer(true),Range("dm_fit"));
  const double m_dm_m20h0m10_min = mean.getVal() - 3.*sl.getVal();
  const double m_dm_m20h0m10_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    a1.setVal(-5.75403e+02); a1.setConstant(kFALSE);
    alphar.setVal(-1.09037e+00);
    f5.setVal(2.86143e-01);  f5.setConstant(kFALSE);
    fCBr.setVal(3.13193e-01);
    mean.setVal(1.42166e-01);
    s1.setVal(6.79274e-04);
    sr.setVal(1.69249e-03);
  } else{
    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
    alphar.setVal(-1.66702e+00);
    fCBr.setVal(4.21899e-01);
    mean.setVal(1.42178e-01);
    s1.setVal(6.48126e-04);
    sr.setVal(1.53696e-03);
  }

  RooFitResult* r = pdf_m_dm.fitTo(dm_m20h0m10_ds,Timer(true),Range("dm_plot"),Save());
  dmdst0.setRange("dmdst0_sig_m20h0m10",m_dm_m20h0m10_min,m_dm_m20h0m10_max);
  const double int_dmdst0_m20h0m10_sig  = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dmdst0_sig_m20h0m10"))->getVal();
  const double int_dmdst0_m20h0m10_plot = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dm_plot"))->getVal();

  RooPlot* m_dmdst0_m20h0m10Frame = dmdst0.frame(Range("dm_plot"),Title("Signal range of m(D*^{0}) - m(D^{0}) for h^{0}=#eta(#rightarrow#gamma#gamma)"));
  dm_m20h0m10_ds.plotOn(m_dmdst0_m20h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_dm.plotOn(m_dmdst0_m20h0m10Frame,LineWidth(2));
  if(!clean) pdf_m_dm.plotOn(m_dmdst0_m20h0m10Frame,Components(pdf_dm_5),LineWidth(2));

  TCanvas* m_dmdst0_m20h0m10_cm = new TCanvas("m_dmdst0_m20h0m10_cm","m_dmdst0_m20h0m10_cm",600,400);
  m_dmdst0_m20h0m10_cm->cd();
  dm_m20h0m10_ds.statOn(m_dmdst0_m20h0m10Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_dmdst0_m20h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_dmdst0_m20h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_dmdst0_m20h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_dmdst0_m20h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_dmdst0_m20h0m10Frame->Draw();
  TLine *mdmdst0_m20h0m10lineLEFT = new TLine(m_dm_m20h0m10_min,0,m_dm_m20h0m10_min,320);
  mdmdst0_m20h0m10lineLEFT->SetLineColor(kRed);
  mdmdst0_m20h0m10lineLEFT->SetLineWidth(2);
  mdmdst0_m20h0m10lineLEFT->SetLineStyle(1);
  mdmdst0_m20h0m10lineLEFT->Draw();
  TLine *mdmdst0_m20h0m10lineRIGHT = new TLine(m_dm_m20h0m10_max,0,m_dm_m20h0m10_max,320);
  mdmdst0_m20h0m10lineRIGHT->SetLineColor(kRed);
  mdmdst0_m20h0m10lineRIGHT->SetLineWidth(2);
  mdmdst0_m20h0m10lineRIGHT->SetLineStyle(1);
  mdmdst0_m20h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_dmdst0_m20h0m10_sig/int_dmdst0_m20h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_dmdst0_m20h0m10_cm->Update();
  m_dmdst0_m20h0m10_cm->Print("pics/m_dmdst0_m20h0m10_cut.eps");
  m_dmdst0_m20h0m10_cm->Print("pics/m_dmdst0_m20h0m10_cut.root");
  // //

  // ** D* eta(->gg) ** //
  const string dm_m20h0m20_cut = string("abs(de)<0.1 && mbc>5.271 && mode == 20 && h0mode == 20");
  RooDataSet dm_m20h0m20_ds("dm_m20h0m20_ds","dm_m20h0m20_ds",tree,argset,dm_m20h0m20_cut.c_str());

  mean.setVal(1.42193e-01);
  sl.setVal(1.19097e-03);
  sr.setVal(1.15833e-03);

  pdf_dm.fitTo(dm_m20h0m20_ds,Timer(true),Range("dm_fit"));
  const double m_dm_m20h0m20_min = mean.getVal() - 3.*sl.getVal();
  const double m_dm_m20h0m20_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    a1.setVal(-5.75403e+02); a1.setConstant(kFALSE);
    alphar.setVal(-1.09037e+00);
    f5.setVal(2.86143e-01);  f5.setConstant(kFALSE);
    fCBr.setVal(3.13193e-01);
    mean.setVal(1.42166e-01);
    s1.setVal(6.79274e-04);
    sr.setVal(1.69249e-03);
  } else{
    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
    alphar.setVal(-1.66702e+00);
    fCBr.setVal(4.21899e-01);
    mean.setVal(1.42178e-01);
    s1.setVal(6.48126e-04);
    sr.setVal(1.53696e-03);
  }

  RooFitResult* r = pdf_m_dm.fitTo(dm_m20h0m20_ds,Timer(true),Range("dm_plot"),Save());
  dmdst0.setRange("dmdst0_sig_m20h0m20",m_dm_m20h0m20_min,m_dm_m20h0m20_max);
  const double int_dmdst0_m20h0m20_sig  = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dmdst0_sig_m20h0m20"))->getVal();
  const double int_dmdst0_m20h0m20_plot = pdf_m_dm.createIntegral(RooArgSet(dmdst0),NormSet(RooArgSet(dmdst0)),Range("dm_plot"))->getVal();

  RooPlot* m_dmdst0_m20h0m20Frame = dmdst0.frame(Range("dm_plot"),Title("Signal range of m(D*^{0}) - m(D^{0}) for h^{0}=#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  dm_m20h0m20_ds.plotOn(m_dmdst0_m20h0m20Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_dm.plotOn(m_dmdst0_m20h0m20Frame,LineWidth(2));
  if(!clean) pdf_m_dm.plotOn(m_dmdst0_m20h0m20Frame,Components(pdf_dm_5),LineWidth(2));

  TCanvas* m_dmdst0_m20h0m20_cm = new TCanvas("m_dmdst0_m20h0m20_cm","m_dmdst0_m20h0m20_cm",600,400);
  m_dmdst0_m20h0m20_cm->cd();
  dm_m20h0m20_ds.statOn(m_dmdst0_m20h0m20Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_dmdst0_m20h0m20Frame->GetXaxis()->SetTitleSize(0.05);
  m_dmdst0_m20h0m20Frame->GetXaxis()->SetTitleOffset(0.85);
  m_dmdst0_m20h0m20Frame->GetXaxis()->SetLabelSize(0.05);
  m_dmdst0_m20h0m20Frame->GetYaxis()->SetTitleOffset(1.6);
  m_dmdst0_m20h0m20Frame->Draw();
  TLine *mdmdst0_m20h0m20lineLEFT = new TLine(m_dm_m20h0m20_min,0,m_dm_m20h0m20_min,320);
  mdmdst0_m20h0m20lineLEFT->SetLineColor(kRed);
  mdmdst0_m20h0m20lineLEFT->SetLineWidth(2);
  mdmdst0_m20h0m20lineLEFT->SetLineStyle(1);
  mdmdst0_m20h0m20lineLEFT->Draw();
  TLine *mdmdst0_m20h0m20lineRIGHT = new TLine(m_dm_m20h0m20_max,0,m_dm_m20h0m20_max,320);
  mdmdst0_m20h0m20lineRIGHT->SetLineColor(kRed);
  mdmdst0_m20h0m20lineRIGHT->SetLineWidth(2);
  mdmdst0_m20h0m20lineRIGHT->SetLineStyle(1);
  mdmdst0_m20h0m20lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_dmdst0_m20h0m20_sig/int_dmdst0_m20h0m20_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_dmdst0_m20h0m20_cm->Update();
  m_dmdst0_m20h0m20_cm->Print("pics/m_dmdst0_m20h0m20_cut.eps");
  m_dmdst0_m20h0m20_cm->Print("pics/m_dmdst0_m20h0m20_cut.root");
  // //

  // ** DeltaE D* pi0 ** //
  const string de_m10h0m10_cut = string("mbc>5.271 && mode == 10 && h0mode == 10 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet de_m10h0m10_ds("de_m10h0m10_ds","de_m10h0m10_ds",tree,argset,de_m10h0m10_cut.c_str());

  mean.setVal(5.13614e-03);
  sl.setVal(5.91763e-02);
  sr.setVal(1.98861e-02);

  pdf_de.fitTo(de_m10h0m10_ds,Timer(true),Range("de_fit"));
  const double m_de_m10h0m10_min = mean.getVal() - 3.*sl.getVal() < -0.1 ? -0.1 : mean.getVal() - 3.*sl.getVal();
  const double m_de_m10h0m10_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    alphal.setVal(7.63119e-01);
    deCBl.setVal(-1.89424e-02);
    fCBl.setVal(8.80614e-01);
    mean.setVal(1.01210e-02);
    s1.setVal(1.40207e-02);
    sl.setVal(3.59638e-02);
  }// else{
//    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
//    alphar.setVal(-1.66702e+00);
//    fCBr.setVal(4.21899e-01);
//    mean.setVal(1.42178e-01);
//    s1.setVal(6.48126e-04);
//    sr.setVal(1.53696e-03);
//  }

  RooFitResult* r = pdf_m_de.fitTo(de_m10h0m10_ds,Timer(true),Range("de_plot"),Save());
  de.setRange("de_sig_m10h0m10",m_de_m10h0m10_min,m_de_m10h0m10_max);
  const double int_de_m10h0m10_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_sig_m10h0m10"))->getVal();
  const double int_de_m10h0m10_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* m_de_m10h0m10Frame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D*^{0}#pi^{0}"));
  de_m10h0m10_ds.plotOn(m_de_m10h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(m_de_m10h0m10Frame,LineWidth(2));
//  pdf_m_de.plotOn(m_de_m10h0m10Frame,Components(g1de),LineStyle(kDashed));
//  pdf_m_de.plotOn(m_de_m10h0m10Frame,Components(CBlde),LineStyle(kDashed));

  TCanvas* m_de_m10h0m10_cm = new TCanvas("m_de_m10h0m10_cm","m_de_m10h0m10_cm",600,400);
  m_de_m10h0m10_cm->cd();
  de_m10h0m10_ds.statOn(m_de_m10h0m10Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_de_m10h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_de_m10h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_de_m10h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_de_m10h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_de_m10h0m10Frame->Draw();
  TLine *mde_m10h0m10lineLEFT = new TLine(m_de_m10h0m10_min,0,m_de_m10h0m10_min,320);
  mde_m10h0m10lineLEFT->SetLineColor(kRed);
  mde_m10h0m10lineLEFT->SetLineWidth(2);
  mde_m10h0m10lineLEFT->SetLineStyle(1);
  mde_m10h0m10lineLEFT->Draw();
  TLine *mde_m10h0m10lineRIGHT = new TLine(m_de_m10h0m10_max,0,m_de_m10h0m10_max,320);
  mde_m10h0m10lineRIGHT->SetLineColor(kRed);
  mde_m10h0m10lineRIGHT->SetLineWidth(2);
  mde_m10h0m10lineRIGHT->SetLineStyle(1);
  mde_m10h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_m10h0m10_sig/int_de_m10h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_de_m10h0m10_cm->Update();
  m_de_m10h0m10_cm->Print("pics/m_de_m10h0m10_cut.eps");
  m_de_m10h0m10_cm->Print("pics/m_de_m10h0m10_cut.root");
  // //

  // ** DeltaE D* eta(->gg) ** //
  const string de_m20h0m10_cut = string("mbc>5.271 && mode == 20 && h0mode == 10 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet de_m20h0m10_ds("de_m20h0m10_ds","de_m20h0m10_ds",tree,argset,de_m20h0m10_cut.c_str());

  mean.setVal(5.13614e-03);
  sl.setVal(5.91763e-02);
  sr.setVal(1.98861e-02);

  pdf_de.fitTo(de_m20h0m10_ds,Timer(true),Range("de_fit"));
  const double m_de_m20h0m10_min = mean.getVal() - 3.*sl.getVal() < -0.1 ? -0.1 : mean.getVal() - 3.*sl.getVal();
  const double m_de_m20h0m10_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    alphal.setVal(1.37634e+00);
    deCBl.setVal(-2.95907e-02);
    fCBl.setVal(5.97526e-01);
    mean.setVal(6.62415e-04);
    s1.setVal(1.88803e-02);
    sl.setVal(3.93993e-02);
  }// else{
//    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
//    alphar.setVal(-1.66702e+00);
//    fCBr.setVal(4.21899e-01);
//    mean.setVal(1.42178e-01);
//    s1.setVal(6.48126e-04);
//    sr.setVal(1.53696e-03);
//  }

  RooFitResult* r = pdf_m_de.fitTo(de_m20h0m10_ds,Timer(true),Range("de_plot"),Save());
  de.setRange("de_sig_m20h0m10",m_de_m20h0m10_min,m_de_m20h0m10_max);
  const double int_de_m20h0m10_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_sig_m20h0m10"))->getVal();
  const double int_de_m20h0m10_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* m_de_m20h0m10Frame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D*^{0}#eta(#rightarrow#gamma#gamma)"));
  de_m20h0m10_ds.plotOn(m_de_m20h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(m_de_m20h0m10Frame,LineWidth(2));
//  pdf_m_de.plotOn(m_de_m20h0m10Frame,Components(g1de),LineStyle(kDashed));
//  pdf_m_de.plotOn(m_de_m20h0m10Frame,Components(CBlde),LineStyle(kDashed));

  TCanvas* m_de_m20h0m10_cm = new TCanvas("m_de_m20h0m10_cm","m_de_m20h0m10_cm",600,400);
  m_de_m20h0m10_cm->cd();
  de_m20h0m10_ds.statOn(m_de_m20h0m10Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_de_m20h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_de_m20h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_de_m20h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_de_m20h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_de_m20h0m10Frame->Draw();
  TLine *mde_m20h0m10lineLEFT = new TLine(m_de_m20h0m10_min,0,m_de_m20h0m10_min,320);
  mde_m20h0m10lineLEFT->SetLineColor(kRed);
  mde_m20h0m10lineLEFT->SetLineWidth(2);
  mde_m20h0m10lineLEFT->SetLineStyle(1);
  mde_m20h0m10lineLEFT->Draw();
  TLine *mde_m20h0m10lineRIGHT = new TLine(m_de_m20h0m10_max,0,m_de_m20h0m10_max,320);
  mde_m20h0m10lineRIGHT->SetLineColor(kRed);
  mde_m20h0m10lineRIGHT->SetLineWidth(2);
  mde_m20h0m10lineRIGHT->SetLineStyle(1);
  mde_m20h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_m20h0m10_sig/int_de_m20h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_de_m20h0m10_cm->Update();
  m_de_m20h0m10_cm->Print("pics/m_de_m20h0m10_cut.eps");
  m_de_m20h0m10_cm->Print("pics/m_de_m20h0m10_cut.root");
  // //

  // ** DeltaE D* eta(->ppp) ** //
  const string de_m20h0m20_cut = string("mbc>5.271 && mode == 20 && h0mode == 20 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet de_m20h0m20_ds("de_m20h0m20_ds","de_m20h0m20_ds",tree,argset,de_m20h0m20_cut.c_str());

  mean.setVal(5.13614e-03);
  sl.setVal(5.91763e-02);
  sr.setVal(1.98861e-02);

  pdf_de.fitTo(de_m20h0m20_ds,Timer(true),Range("de_fit"));
  const double m_de_m20h0m20_min = mean.getVal() - 3.*sl.getVal() < -0.1 ? -0.1 : mean.getVal() - 3.*sl.getVal();
  const double m_de_m20h0m20_max = mean.getVal() + 3.*sr.getVal();

  if(!clean){
    alphal.setVal(1.37634e+00);
    deCBl.setVal(-2.95907e-02);
    fCBl.setVal(5.97526e-01);
    mean.setVal(6.62415e-04);
    s1.setVal(1.88803e-02);
    sl.setVal(3.93993e-02);
  }// else{
//    f5.setVal(0); f5.setConstant(kTRUE); a1.setConstant(kTRUE);
//    alphar.setVal(-1.66702e+00);
//    fCBr.setVal(4.21899e-01);
//    mean.setVal(1.42178e-01);
//    s1.setVal(6.48126e-04);
//    sr.setVal(1.53696e-03);
//  }

  RooFitResult* r = pdf_m_de.fitTo(de_m20h0m20_ds,Timer(true),Range("de_plot"),Save());
  de.setRange("de_sig_m20h0m20",m_de_m20h0m20_min,m_de_m20h0m20_max);
  const double int_de_m20h0m20_sig  = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_sig_m20h0m20"))->getVal();
  const double int_de_m20h0m20_plot = pdf_m_de.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("de_plot"))->getVal();

  RooPlot* m_de_m20h0m20Frame = de.frame(Range("de_plot"),Title("Signal range of #DeltaE for D*^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  de_m20h0m20_ds.plotOn(m_de_m20h0m20Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_de.plotOn(m_de_m20h0m20Frame,LineWidth(2));
//  pdf_m_de.plotOn(m_de_m20h0m20Frame,Components(g1de),LineStyle(kDashed));
//  pdf_m_de.plotOn(m_de_m20h0m20Frame,Components(CBlde),LineStyle(kDashed));

  TCanvas* m_de_m20h0m20_cm = new TCanvas("m_de_m20h0m20_cm","m_de_m20h0m20_cm",600,400);
  m_de_m20h0m20_cm->cd();
  de_m20h0m20_ds.statOn(m_de_m20h0m20Frame,Layout(0.6,0.98,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_de_m20h0m20Frame->GetXaxis()->SetTitleSize(0.05);
  m_de_m20h0m20Frame->GetXaxis()->SetTitleOffset(0.85);
  m_de_m20h0m20Frame->GetXaxis()->SetLabelSize(0.05);
  m_de_m20h0m20Frame->GetYaxis()->SetTitleOffset(1.6);
  m_de_m20h0m20Frame->Draw();
  TLine *mde_m20h0m20lineLEFT = new TLine(m_de_m20h0m20_min,0,m_de_m20h0m20_min,320);
  mde_m20h0m20lineLEFT->SetLineColor(kRed);
  mde_m20h0m20lineLEFT->SetLineWidth(2);
  mde_m20h0m20lineLEFT->SetLineStyle(1);
  mde_m20h0m20lineLEFT->Draw();
  TLine *mde_m20h0m20lineRIGHT = new TLine(m_de_m20h0m20_max,0,m_de_m20h0m20_max,320);
  mde_m20h0m20lineRIGHT->SetLineColor(kRed);
  mde_m20h0m20lineRIGHT->SetLineWidth(2);
  mde_m20h0m20lineRIGHT->SetLineStyle(1);
  mde_m20h0m20lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_de_m20h0m20_sig/int_de_m20h0m20_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_de_m20h0m20_cm->Update();
  m_de_m20h0m20_cm->Print("pics/m_de_m20h0m20_cut.eps");
  m_de_m20h0m20_cm->Print("pics/m_de_m20h0m20_cut.root");
  // //

  // ** Mbc D* pi0 ** //
  const string mbc_m10h0m10_cut = string("de>-0.1 && de<0.0647985 && mode == 10 && h0mode == 10 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet mbc_m10h0m10_ds("mbc_m10h0m10_ds","mbc_m10h0m10_ds",tree,argset,mbc_m10h0m10_cut.c_str());

  mean.setVal(5.28013e+00);
  sl.setVal(3.47147e-03);
  sr.setVal(2.51821e-03);

  pdf_mbc.fitTo(mbc_m10h0m10_ds,Timer(true),Range("mbc_fit"));
  const double m_mbc_m10h0m10_min = mean.getVal() - 2.5*sl.getVal();
  const double m_mbc_m10h0m10_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(1.30133e-01);
  mean.setVal(5.27992e+00);
  s1.setVal(3.02195e-03);

  RooFitResult* r = pdf_m_mbc.fitTo(mbc_m10h0m10_ds,Timer(true),Range("mbc_plot"),Save());
  mbc.setRange("mbc_sig_m10h0m10",m_mbc_m10h0m10_min,m_mbc_m10h0m10_max);
  const double int_mbc_m10h0m10_sig  = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_m10h0m10"))->getVal();
  const double int_mbc_m10h0m10_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* m_mbc_m10h0m10Frame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D*^{0}#pi^{0}"));
  mbc_m10h0m10_ds.plotOn(m_mbc_m10h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(m_mbc_m10h0m10Frame,LineWidth(2));

  TCanvas* m_mbc_m10h0m10_cm = new TCanvas("m_mbc_m10h0m10_cm","m_mbc_m10h0m10_cm",600,400);
  m_mbc_m10h0m10_cm->cd();
  mbc_m10h0m10_ds.statOn(m_mbc_m10h0m10Frame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_mbc_m10h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_mbc_m10h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_mbc_m10h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_mbc_m10h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_mbc_m10h0m10Frame->Draw();
  TLine *mmbc_m10h0m10lineLEFT = new TLine(m_mbc_m10h0m10_min,0,m_mbc_m10h0m10_min,320);
  mmbc_m10h0m10lineLEFT->SetLineColor(kRed);
  mmbc_m10h0m10lineLEFT->SetLineWidth(2);
  mmbc_m10h0m10lineLEFT->SetLineStyle(1);
  mmbc_m10h0m10lineLEFT->Draw();
  TLine *mmbc_m10h0m10lineRIGHT = new TLine(m_mbc_m10h0m10_max,0,m_mbc_m10h0m10_max,320);
  mmbc_m10h0m10lineRIGHT->SetLineColor(kRed);
  mmbc_m10h0m10lineRIGHT->SetLineWidth(2);
  mmbc_m10h0m10lineRIGHT->SetLineStyle(1);
  mmbc_m10h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_m10h0m10_sig/int_mbc_m10h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_mbc_m10h0m10_cm->Update();
  m_mbc_m10h0m10_cm->Print("pics/m_mbc_m10h0m10_cut.eps");
  m_mbc_m10h0m10_cm->Print("pics/m_mbc_m10h0m10_cut.root");
  // //

  // ** Mbc D* eta(->gg) ** //
  const string mbc_m20h0m10_cut = string("de>-0.0853615 && de<0.0630641 && mode == 20 && h0mode == 10 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet mbc_m20h0m10_ds("mbc_m20h0m10_ds","mbc_m20h0m10_ds",tree,argset,mbc_m20h0m10_cut.c_str());

  mean.setVal(5.28013e+00);
  sl.setVal(3.47147e-03);
  sr.setVal(2.51821e-03);

  pdf_mbc.fitTo(mbc_m20h0m10_ds,Timer(true),Range("mbc_fit"));
  const double m_mbc_m20h0m10_min = mean.getVal() - 2.5*sl.getVal();
  const double m_mbc_m20h0m10_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(1.30133e-01);
  mean.setVal(5.27992e+00);
  s1.setVal(3.02195e-03);

  RooFitResult* r = pdf_m_mbc.fitTo(mbc_m20h0m10_ds,Timer(true),Range("mbc_plot"),Save());
  mbc.setRange("mbc_sig_m20h0m10",m_mbc_m20h0m10_min,m_mbc_m20h0m10_max);
  const double int_mbc_m20h0m10_sig  = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_m20h0m10"))->getVal();
  const double int_mbc_m20h0m10_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* m_mbc_m20h0m10Frame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D*^{0}#eta(#rightarrow#gamma#gamma)"));
  mbc_m20h0m10_ds.plotOn(m_mbc_m20h0m10Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(m_mbc_m20h0m10Frame,LineWidth(2));

  TCanvas* m_mbc_m20h0m10_cm = new TCanvas("m_mbc_m20h0m10_cm","m_mbc_m20h0m10_cm",600,400);
  m_mbc_m20h0m10_cm->cd();
  mbc_m20h0m10_ds.statOn(m_mbc_m20h0m10Frame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_mbc_m20h0m10Frame->GetXaxis()->SetTitleSize(0.05);
  m_mbc_m20h0m10Frame->GetXaxis()->SetTitleOffset(0.85);
  m_mbc_m20h0m10Frame->GetXaxis()->SetLabelSize(0.05);
  m_mbc_m20h0m10Frame->GetYaxis()->SetTitleOffset(1.6);
  m_mbc_m20h0m10Frame->Draw();
  TLine *mmbc_m20h0m10lineLEFT = new TLine(m_mbc_m20h0m10_min,0,m_mbc_m20h0m10_min,320);
  mmbc_m20h0m10lineLEFT->SetLineColor(kRed);
  mmbc_m20h0m10lineLEFT->SetLineWidth(2);
  mmbc_m20h0m10lineLEFT->SetLineStyle(1);
  mmbc_m20h0m10lineLEFT->Draw();
  TLine *mmbc_m20h0m10lineRIGHT = new TLine(m_mbc_m20h0m10_max,0,m_mbc_m20h0m10_max,320);
  mmbc_m20h0m10lineRIGHT->SetLineColor(kRed);
  mmbc_m20h0m10lineRIGHT->SetLineWidth(2);
  mmbc_m20h0m10lineRIGHT->SetLineStyle(1);
  mmbc_m20h0m10lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_m20h0m10_sig/int_mbc_m20h0m10_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_mbc_m20h0m10_cm->Update();
  m_mbc_m20h0m10_cm->Print("pics/m_mbc_m20h0m10_cut.eps");
  m_mbc_m20h0m10_cm->Print("pics/m_mbc_m20h0m10_cut.root");
  // //

  // ** Mbc D* eta(->ppp) ** //
  const string mbc_m20h0m20_cut = string("de>-0.0853615 && de<0.0630641 && mode == 20 && h0mode == 10 && dmdst0>0.13862 && dmdst0<0.145667");
  RooDataSet mbc_m20h0m20_ds("mbc_m20h0m20_ds","mbc_m20h0m20_ds",tree,argset,mbc_m20h0m20_cut.c_str());

  mean.setVal(5.28013e+00);
  sl.setVal(3.47147e-03);
  sr.setVal(2.51821e-03);

  pdf_mbc.fitTo(mbc_m20h0m20_ds,Timer(true),Range("mbc_fit"));
  const double m_mbc_m20h0m20_min = mean.getVal() - 2.5*sl.getVal();
  const double m_mbc_m20h0m20_max = mean.getVal() + 3.*sr.getVal();

  alpha_mbc.setVal(1.30133e-01);
  mean.setVal(5.27992e+00);
  s1.setVal(3.02195e-03);

  RooFitResult* r = pdf_m_mbc.fitTo(mbc_m20h0m20_ds,Timer(true),Range("mbc_plot"),Save());
  mbc.setRange("mbc_sig_m20h0m20",m_mbc_m20h0m20_min,m_mbc_m20h0m20_max);
  const double int_mbc_m20h0m20_sig  = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_sig_m20h0m20"))->getVal();
  const double int_mbc_m20h0m20_plot = pdf_m_mbc.createIntegral(RooArgSet(mbc),NormSet(RooArgSet(mbc)),Range("mbc_plot"))->getVal();

  RooPlot* m_mbc_m20h0m20Frame = mbc.frame(Range("mbc_plot"),Title("Signal range of M_{bc} for D*^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})"));
  mbc_m20h0m20_ds.plotOn(m_mbc_m20h0m20Frame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_m_mbc.plotOn(m_mbc_m20h0m20Frame,LineWidth(2));

  TCanvas* m_mbc_m20h0m20_cm = new TCanvas("m_mbc_m20h0m20_cm","m_mbc_m20h0m20_cm",600,400);
  m_mbc_m20h0m20_cm->cd();
  mbc_m20h0m20_ds.statOn(m_mbc_m20h0m20Frame,Layout(0.2,0.6,0.9));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  m_mbc_m20h0m20Frame->GetXaxis()->SetTitleSize(0.05);
  m_mbc_m20h0m20Frame->GetXaxis()->SetTitleOffset(0.85);
  m_mbc_m20h0m20Frame->GetXaxis()->SetLabelSize(0.05);
  m_mbc_m20h0m20Frame->GetYaxis()->SetTitleOffset(1.6);
  m_mbc_m20h0m20Frame->Draw();
  TLine *mmbc_m20h0m20lineLEFT = new TLine(m_mbc_m20h0m20_min,0,m_mbc_m20h0m20_min,320);
  mmbc_m20h0m20lineLEFT->SetLineColor(kRed);
  mmbc_m20h0m20lineLEFT->SetLineWidth(2);
  mmbc_m20h0m20lineLEFT->SetLineStyle(1);
  mmbc_m20h0m20lineLEFT->Draw();
  TLine *mmbc_m20h0m20lineRIGHT = new TLine(m_mbc_m20h0m20_max,0,m_mbc_m20h0m20_max,320);
  mmbc_m20h0m20lineRIGHT->SetLineColor(kRed);
  mmbc_m20h0m20lineRIGHT->SetLineWidth(2);
  mmbc_m20h0m20lineRIGHT->SetLineStyle(1);
  mmbc_m20h0m20lineRIGHT->Draw();

  TPaveText *pt = new TPaveText(0.2,0.6,0.4,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "Eff = " << std::setprecision(2) << int_mbc_m20h0m20_sig/int_mbc_m20h0m20_plot;
  pt->AddText(out.str().c_str());
  pt->Draw();

  m_mbc_m20h0m20_cm->Update();
  m_mbc_m20h0m20_cm->Print("pics/m_mbc_m20h0m20_cut.eps");
  m_mbc_m20h0m20_cm->Print("pics/m_mbc_m20h0m20_cut.root");
  // //

  cout << "dm(D*) limits:" << endl;
  cout << " pi0:        " << m_dm_m10h0m10_min << " " << m_dm_m10h0m10_max << endl;
  cout << " eta(->gg):  " << m_dm_m20h0m10_min << " " << m_dm_m20h0m10_max << endl;
  cout << " eta(->ppp): " << m_dm_m20h0m20_min << " " << m_dm_m20h0m20_max << endl;
  cout << "dE limits:" << endl;
  cout << " pi0:        " << m_de_m10h0m10_min << " " << m_de_m10h0m10_max << endl;
  cout << " eta(->gg):  " << m_de_m20h0m10_min << " " << m_de_m20h0m10_max << endl;
  cout << " eta(->ppp): " << m_de_m20h0m20_min << " " << m_de_m20h0m20_max << endl;
  cout << "Mbc limits:" << endl;
  cout << " pi0:        " << m_mbc_m10h0m10_min << " " << m_mbc_m10h0m10_max << endl;
  cout << " eta(->gg):  " << m_mbc_m20h0m10_min << " " << m_mbc_m20h0m10_max << endl;
  cout << " eta(->ppp): " << m_mbc_m20h0m20_min << " " << m_mbc_m20h0m20_max << endl;
}










