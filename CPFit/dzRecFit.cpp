using namespace RooFit;

TTree* GetGoodTree(TTree* omegatree,const int mode){
  Double_t m_dt,m_sz,m_h,m_mbc,m_de,m_bdtg;
  Int_t m_good_icpv,m_bin,m_exp,m_flv,m_ndf;
//  if(mode>2 && mode != 6){
    omegatree->SetBranchAddress("ndf_z_sig",&m_ndf);
    omegatree->SetBranchAddress("chisq_z_sig",&m_h);
//  }
  omegatree->SetBranchAddress("mbc",&m_mbc);
  omegatree->SetBranchAddress("de",&m_de);
  omegatree->SetBranchAddress("bdtg",&m_bdtg);
  omegatree->SetBranchAddress("bin_mc",&m_bin);
  omegatree->SetBranchAddress("exp",&m_exp);
  omegatree->SetBranchAddress("flv_mc",&m_flv);

  Double_t m_z_sig_mc,m_z_sig;

  omegatree->SetBranchAddress("z_sig_mc",&m_z_sig_mc);
  if(mode != 6){
    omegatree->SetBranchAddress("z_sig",&m_z_sig);
    omegatree->SetBranchAddress("sz_sig",&m_sz);
  } else{
    omegatree->SetBranchAddress("z_sig_d0",&m_z_sig);
    omegatree->SetBranchAddress("sz_sig_d0",&m_sz);
  }

  Int_t m_b0f, m_mode, m_h0mode;
  if(mode == 5 || mode == 6){
    omegatree->SetBranchAddress("bpf",&m_b0f);
    if(mode == 5){
      omegatree->SetBranchAddress("good_icpv_mlt",&m_good_icpv);
    } else{
      omegatree->SetBranchAddress("good_icpv_sgl",&m_good_icpv);
    }
  } else{
    omegatree->SetBranchAddress("b0f",&m_b0f);
    omegatree->SetBranchAddress("good_icpv",&m_good_icpv);
    omegatree->SetBranchAddress("mode",&m_mode);
    omegatree->SetBranchAddress("h0mode",&m_h0mode);
  }

  TTree* tree = new TTree("Event","Event");
  tree->Branch("dt",&m_dt,"dt/D");
  tree->Branch("sz",&m_sz,"sz/D");
  if(mode > 2 && mode != 6){
    tree->Branch("h",&m_h,"h/D");
  }
  if(mode != 5 && mode != 6){
    tree->Branch("mode",&m_mode,"mode/I");
    tree->Branch("h0mode",&m_h0mode,"h0mode/I");
  }
  tree->Branch("mbc",&m_mbc,"mbc/D");
  tree->Branch("de",&m_de,"de/D");
  tree->Branch("bdtg",&m_bdtg,"bdtg/D");

  tree->Branch("good_icpv",&m_good_icpv,"good_icpv/I");
  tree->Branch("b0f",&m_b0f,"b0f/I");
  tree->Branch("bin_mc",&m_bin,"bin_mc/I");
  tree->Branch("flv_mc",&m_flv,"flv_mc/I");
  tree->Branch("exp",&m_exp,"exp/I");

  const double mm2ps = 7.84857;
  const int NTot = omegatree->GetEntries();
  for(int i=0; i<NTot; i++){
    omegatree->GetEvent(i);
//    cout << m_dt << " " << m_sz << " " << m_h << " " << m_ndf << endl;
    if(m_h<0 || !m_bin) continue;
    if(mode>2 && m_ndf == 0 || mode<=2 && m_ndf>0) continue;
    if(!(m_mbc>0 || m_mbc<0)) continue;
    if(!(m_sz>0 || m_sz<0)) continue;
    m_dt = (m_z_sig-m_z_sig_mc)*mm2ps;
//    if(mode>4) m_dt *= 10;
    m_sz *= mm2ps;
    if(mode>2 && mode != 6) m_h /= m_ndf;
    tree->Fill();
  }
  return tree;
}

void dzRecFit(const int Mode = 4,const int SVD = 2, const double dt_max = 70){
    // mode 1 -> D0 pi0
    // mode 2 -> D0 eta -> gg
    // mode 3 -> D0 eta -> pi+pi-pi0
    // mode 4 -> D0 omega
    // mode 5 -> D0 pi+ (mlt)
    // mode 6 -> D0 pi+ (sgl)
  TChain *omegatree = new TChain("TEvent");
//  TFile *omegafile;
  string label;
  int _mode=3, _h0mode=20;
  switch(Mode){
  case 1:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
    label = string("#pi^{0}");
    _mode = 1; _h0mode = 10;
    break;
  case 2:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    label = string("#eta#rightarrow#gamma#gamma");
    _mode = 2; _h0mode = 10;
    break;
  case 3:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    _mode = 2; _h0mode = 20;
    break;
  case 4:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s6_full.root");
    label = string("#omega");
    _mode = 3; _h0mode = 20;
    break;
  case 5:
    omegatree->Add("/home/vitaly/B0toDh0/Bp2D0pi/FIL_bp2d0pip_sigmc.root");
    label = string("D^{0}#pi^{+} (multiple)");
    _mode = -1; _h0mode = -10;
    break;
  case 6:
    omegatree->Add("/home/vitaly/B0toDh0/Bp2D0pi/FIL_bp2d0pip_sigmc.root");
    label = string("D^{0}#pi^{+} (single)");
    _mode = -1; _h0mode = -10;
    break;
  }
  TTree* tree = GetGoodTree(omegatree,Mode);

  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("2",2);
  b0f.defineType("3",3);
  b0f.defineType("4",4);
  b0f.defineType("badpi0",5);
  b0f.defineType("10",10);
  argset.add(b0f);
  if(Mode<5){
    RooCategory mode("mode","mode");
    mode.defineType("omega",_mode);
    argset.add(mode);

    RooCategory h0mode("h0mode","h0mode");
    h0mode.defineType("h0mode",_h0mode);
    argset.add(h0mode);
  }

  RooCategory exp("exp","exp");
  if(SVD == 1){
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
  } else if(SVD == 2){
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

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good",1);
  argset.add(good_icpv);

//  const double dt_max =  10;
  const double dt_min = -dt_max;
  RooRealVar dt("dt","#deltat^{sig}",dt_min,dt_max,"ps"); argset.add(dt);
  RooRealVar sz("sz","sz",0.0,1);                         argset.add(sz);
  if(Mode>2 && Mode != 6) {RooRealVar h("h","h",0.,10000);argset.add(h);}
  RooRealVar mbc("mbc","mbc",5.272,5.286);                argset.add(mbc);
  RooRealVar de("de","de",-0.06,0.05);                    argset.add(de);
  RooRealVar bdtg("bdtg","bdtg",0.,1.);                   argset.add(bdtg);

  stringstream out;
  out.str("");
//  out << dzname << ">0 || 0>" << dzname << " && 0<" << szname << " && 0<" << xiname;
  RooDataSet ds("ds","ds",tree,argset,out.str().c_str());
  ds.Print();

//  RooConstVar MM2PS("mm2ps","mm2ps",mm2ps);
//  RooFormulaVar dt("dt","dt","@0*@1",RooArgSet(dz_mc_sig,MM2PS));

  const bool true_params = false;
  double Srec[2];
  double Smn_rec, Stl_rec, ftl_rec;
  if(SVD == 2){
    Srec[0] = +9.271430e-01; Srec[1] = +2.103700e-01;
    Smn_rec = +1.053040e+00; //* Smn_rec */
    Stl_rec = +4.320600e+00; //* Stl_rec */
    ftl_rec = +7.068970e-02; //* ftl_rec */
  } else{
    Srec[0] = +9.626350e-01; Srec[1] = +1.985560e-01;
    Smn_rec = +1.109750e+00; //* Smn_rec */
    Stl_rec = +1.000000e+00; //* Stl_rec */
    ftl_rec = +0.000000e+00; //* ftl_rec */
  }

  RooRealVar dz0sig("dz0sig","dz0sig",0,-1.27,1.27,"mm"); dz0sig.setConstant(kTRUE);
  RooRealVar sOtlr("sOtlr","sOtlr",3.063190e+01,0.,1000.); sOtlr.setConstant(kTRUE);
  RooRealVar fOtlr("fOtlr","fOtlr",8.843970e-05,0.,1.);    //fOtlr.setConstant(kTRUE);
//  RooRealVar fOtlr("fOtlr","fOtlr",0.,0.,1.);    fOtlr.setConstant(kTRUE);
  RooGaussModel Otlr("Otlr","Otlr",dt,dz0sig,sOtlr);
//  dz0sig.setConstant(kTRUE);

  if(Mode > 2 && Mode != 6){
//  RooRealVar dz1sigOm("dz1sigOm","dz1sigOm",0,-1.27,1.27,"mm");// dz1sigOm.setConstant(kTRUE);
  RooRealVar s0sig("s0sig","s0sig",Srec[0],0.0,10.5); if(true_params) s0sig.setConstant(kTRUE);
  RooRealVar s1sig("s1sig","s1sig",Srec[1],0.,10.);   if(true_params) s1sig.setConstant(kTRUE);
//  RooRealVar s2sig("s2sig","s2sig",0,0.0,10.5);       if(true_params) s2sig.setConstant(kTRUE);

  RooFormulaVar dz0("dz0","@0*@1",RooArgSet(dz0sig,h));

//  RooRealVar s2sigOm("s2sigOm","s2sigOm",1.5,1.,10.5);// if(true_params) s0sigOm.setConstant(kTRUE);
//  RooRealVar s3sigOm("s3sigOm","s3sigOm",1.5,0.,10.);    // if(true_params) s1sigOm.setConstant(kTRUE);

//  RooRealVar s4sigOm("s4sigOm","s4sigOm",2.5,1.5,10.5);// if(true_params) s0sigOm.setConstant(kTRUE);
//  RooRealVar s5sigOm("s5sigOm","s5sigOm",2.5,0.,10.);    // if(true_params) s1sigOm.setConstant(kTRUE);

//  RooFormulaVar mean("mean","@0+@1*@2",RooArgSet(dz0sigOm,dz1sigOm,xi_sig));
//  RooFormulaVar mean0("mean0","@0*@1",RooArgSet(dz0sigOm,xi_sig));
//  RooFormulaVar mean1("mean1","@0*@1",RooArgSet(dz1sigOm,xi_sig));


//  RooFormulaVar sigma("sigma","@4+(@0+@1*@2)*@3",RooArgSet(s0sig,s1sig,h,sz,s2sig));
  RooFormulaVar sigma("sigma","(@0+@1*@2)*@3",RooArgSet(s0sig,s1sig,h,sz));


  //  RooFormulaVar sigma_tail("sigma_tail","(@0+@1*@2)*@3",RooArgSet(s2sigOm,s3sigOm,xi_sig,sz_sig));
//  RooFormulaVar sigma_tail2("sigma_tail2","(@0+@1*@2)*@3",RooArgSet(s4sigOm,s5sigOm,xi_sig,sz_sig));
//  RooGaussModel pdf_tail("pdf_tail","pdf_tail",dz_mc_sig,mean1,sigma_tail);
//  RooGaussModel pdf_tail2("pdf_tail2","pdf_tail2",dz_mc_sig,mean,sigma_tail2);

//  RooRealVar f_tail("f_tail","f_tail",0.01,0.,1.);
//  RooRealVar f_tail2("f_tail2","f_tail2",0.01,0.,1.);
//  RooAddPdf pdf("pdf","pdf",RooArgList(Otlr,pdf_tail,pdf_tail2,pdf1),RooArgSet(fOtlr,f_tail,f_tail2));


  RooGaussian pdf("pdf","pdf",dt,dz0,sigma);
//  RooAddPdf pdf("pdf","pdf",RooArgList(Otlr,pdf1),RooArgSet(fOtlr));


  //  RooAddPdf pdf("pdf","pdf",RooArgList(Otlr,pdf_tail,pdf1),RooArgSet(fOtlr,f_tail));

//  RooFitResult* fitRes = pdf.fitTo(omegads,Verbose(),Minimizer("Minuit2","migrad"),ConditionalObservables(RooArgSet(xi_sig,sz_sig)));
  //RooFitResult* fitRes = pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(h,sz)));
  pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(h,sz)));
//  pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(h,sz)));
//  RooFitResult* fitRes = pdf.chi2fitTo(omegads,ConditionalObservables(RooArgSet(xi_sig,sz_sig)));
  } else{
//  RooRealVar dzsig("dzsig","dzsig",0,3.*dt_min,3.*dt_max,"mm");if(true_params) dzsig.setConstant(kTRUE);
  RooRealVar sMainSig("sMainSig","sMainSig",Smn_rec,0.01,5.);  if(true_params) sMainSig.setConstant(kTRUE);
  RooRealVar sTailSig("sTailSig","sTailSig",Stl_rec,0.,27.);   if(true_params) sTailSig.setConstant(kTRUE);
  RooGaussModel gMainSig("gMainSig","gMainSig",dt,dz0sig,sMainSig,sz);
  RooGaussModel gTailSig("gTailSig","gTailSig",dt,dz0sig,sTailSig,sz);
  RooRealVar f_tail("f_tail","f_tail",ftl_rec,0.,1.);          if(true_params) f_tail.setConstant(kTRUE);
  if(SVD == 1){
    RooAddPdf pdf("pdf","pdf",RooArgList(Otlr,gMainSig),RooArgSet(fOtlr));
  } else{
    RooAddPdf pdf("pdf","pdf",RooArgList(Otlr,gTailSig,gMainSig),RooArgSet(fOtlr,f_tail));
  }

  RooFitResult* fitRes = pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(sz)));
  }

  dt.setBins(250);
///////////
// Plots //
///////////
  RooPlot* dzFrameSig = dt.frame();
  ds.plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
//  pdf.plotOn(dzFrameSig,Components(Otlr),LineStyle(kDashed),LineColor(kBlue),ProjWData(ds));
//  pdf.plotOn(dzFrameSig,Components(RooArgSet(Otlr,pdf_tail)),LineStyle(kDashed),LineColor(kBlue),ProjWData(omegads));
//  pdf.plotOn(dzFrameSig,Components(RooArgSet(Otlr,pdf_tail,pdf_tail2)),LineStyle(kDashed),LineColor(kBlue),ProjWData(omegads));
//  if(mode < 2) pdf.plotOn(dzFrameSig,Components(gTailSigOm),LineStyle(kDashed),LineColor(kBlue),ProjWData(omegads));
  pdf.plotOn(dzFrameSig,LineColor(kBlue),ProjWData(ds));

  RooHist* hdzpull = dzFrameSig->pullHist();
  RooPlot* dzPull = dt.frame(Title("#deltat_{sig} pull distribution"));
  dzPull->addPlotable(hdzpull,"P");
  dzPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* dz_cm_sig = new TCanvas("#deltat_{sig} (ps)","#deltat_{sig} (ps)",700,800);
  dz_cm_sig->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd(); pad1->SetLeftMargin(0.15); pad1->SetFillColor(0);
  ds.statOn(dzFrameSig,Label(label.c_str()),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetGrid();
  pad1->SetLogy();
  dzFrameSig->GetXaxis()->SetTitleOffset(0.85);
  dzFrameSig->GetXaxis()->SetLabelSize(0.04);
  dzFrameSig->GetXaxis()->SetTitleSize(0.05);
  dzFrameSig->GetYaxis()->SetTitleOffset(1.6);
//  dzFrameSig->GetYaxis()->SetRangeUser(0.1,100000);
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
  TLine *dz_lineUP = new TLine(dt_min,3,dt_max,3);
  dz_lineUP->SetLineColor(kBlue);
  dz_lineUP->SetLineStyle(2);
  dz_lineUP->Draw();
  TLine *dz_line = new TLine(dt_min,0,dt_max,0);
  dz_line->SetLineColor(kBlue);
  dz_line->SetLineStyle(1);
  dz_line->SetLineWidth((Width_t)2.);
  dz_line->Draw();
  TLine *dz_lineDOWN = new TLine(dt_min,-3,dt_max,-3);
  dz_lineDOWN->SetLineColor(kBlue);
  dz_lineDOWN->SetLineStyle(2);
  dz_lineDOWN->Draw();

  dz_cm_sig->Update();
//  dz_cm_sig->Print("dz_sig_fit.png");
//  dz_cm_sig->Print("dz_sig_fit.root");

}

