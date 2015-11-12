using namespace RooFit;

TTree* GetGoodTree(TTree* omegatree,const bool single){
  Double_t m_dt,m_sz,m_h,m_mbc,m_de,m_bdtg;
  Int_t m_good_icpv,m_bin,m_exp,m_flv,m_ndf,m_nptag;
  omegatree->SetBranchAddress("nptag",&m_nptag);
  omegatree->SetBranchAddress("dz_mc_asc",&m_dt);
  omegatree->SetBranchAddress("sz_asc",&m_sz);
  if(!single){
    omegatree->SetBranchAddress("ndf_z_asc",&m_ndf);
    omegatree->SetBranchAddress("chisq_z_asc",&m_h);
  }
  omegatree->SetBranchAddress("mbc",&m_mbc);
  omegatree->SetBranchAddress("de",&m_de);
  omegatree->SetBranchAddress("bdtg",&m_bdtg);
  omegatree->SetBranchAddress("good_icpv",&m_good_icpv);
  omegatree->SetBranchAddress("bin_mc",&m_bin);
  omegatree->SetBranchAddress("exp",&m_exp);
  omegatree->SetBranchAddress("flv_mc",&m_flv);

  TTree* tree = new TTree("Event","Event");
  tree->Branch("nptag",&m_nptag,"nptag/I");
  tree->Branch("dt",&m_dt,"dt/D");
  tree->Branch("sz",&m_sz,"sz/D");
  if(!single){
    tree->Branch("h",&m_h,"h/D");
  }
  tree->Branch("mbc",&m_mbc,"mbc/D");
  tree->Branch("de",&m_de,"de/D");
  tree->Branch("bdtg",&m_bdtg,"bdtg/D");

  tree->Branch("good_icpv",&m_good_icpv,"good_icpv/I");
  tree->Branch("bin_mc",&m_bin,"bin_mc/I");
  tree->Branch("flv_mc",&m_flv,"flv_mc/I");
  tree->Branch("exp",&m_exp,"exp/I");

  const double mm2ps = 7.84857;
  const int NTot = omegatree->GetEntries();
  for(int i=0; i<NTot; i++){
    omegatree->GetEvent(i);
//    cout << m_dt << " " << m_sz << " " << m_h << " " << m_ndf << endl;
    if(m_nptag) continue;
    if(m_h<0 || !m_bin) continue;
    if(!single && m_ndf == 0 || single && m_ndf>0) continue;
    if(!(m_mbc>0 || m_mbc<0)) continue;
    if(!(m_sz>0 || m_sz<0)) continue;
    m_dt *= mm2ps;
    m_sz *= mm2ps;
    if(!single) m_h /= m_ndf;
    tree->Fill();
  }
  return tree;
}

void dzFitAsc(const int SVD = 2,const bool single){
//  TFile *omegafile;
  string label;
  if(single)   label = string("single, ");
  else         label = string("multi, ");
  if(SVD == 1) label += string("SVD1");
  else         label += string("SVD2");
  TChain *omegatree = new TChain("TEvent","TEvent");
//  switch(mode){
//  case 1:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
//    label = string("#pi^{0}");
//    break;
//  case 2:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
//    label = string("#eta#rightarrow#gamma#gamma");
//    break;
//  case 3:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
//    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
//    break;
//  case 4:
    omegatree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s6_full.root");
//    label = string("#omega");
//    break;
//  }
//  TTree *omegatree = (TTree*)omegafile->Get("TEvent");
//  omegatree->Print();
  TTree* tree = GetGoodTree(omegatree,single);
//  tree->Print();

  RooArgSet argset;

  RooCategory nptag("nptag","nptag");
  nptag.defineType("noNP",0);
  argset.add(nptag);

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
  } else{ // SVD2
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

  const double dt_max =  2;// ps
  const double dt_min = -dt_max;
  const bool true_params = false;
  RooRealVar dt("dt","#deltat^{tag}",dt_min,dt_max,"ps"); argset.add(dt);
  RooRealVar sz("sz","sz",0.00,1000);                     argset.add(sz);
  RooRealVar h("h","h",0.,10000);                         argset.add(h);

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  RooRealVar dtasc("dtasc","dtasc",0,dt_min,dt_max,"pm"); if(true_params) dtasc.setConstant(kTRUE);
  RooRealVar s_otlr("s_otlr","s_otlr",30.,0.1,100.);      s_otlr.setConstant(kTRUE);
  RooRealVar z0_otlr("z0_otlr","z0_otlr",0.,-0.1,0.1);    z0_otlr.setConstant(kTRUE);
  RooGaussian g_otlr("g_otlr","g_otlr",dt,z0_otlr,s_otlr);
  RooRealVar f_otlr("f_otlr","f_otlr",0.0,0.,1.);// f_otlr.setConstant(kTRUE);

  double Sasc[2];
  double Smn_asc, Stl_asc, ftl_asc;
  if(SVD == 2){
    Sasc[0] = +1.02836e+00; Sasc[1] = +1.96229e-01;
    Smn_asc = +1.053040e+00; //* Smn_asc */
    Stl_asc = +4.320600e+00; //* Stl_asc */
    ftl_asc = +7.068970e-02; //* ftl_asc */
  } else{
    Sasc[0] = +7.290850e-01; Sasc[1] = +1.719270e-01;
    Smn_asc = +1.16355e+00;  //* Smn_asc */
    Stl_asc = +1.000000e+00; //* Stl_asc */
    ftl_asc = +0.000000e+00; //* ftl_asc */
  }

  if(!single){
    RooRealVar s0asc("s0asc","s0asc",Sasc[0],0.1,5.); if(true_params) s0asc.setConstant(kTRUE);
    RooRealVar s1asc("s1asc","s1asc",Sasc[1],0.,27.); if(true_params) s1asc.setConstant(kTRUE);

    RooFormulaVar sigma("sigma","(@0+@1*@2)*@3",RooArgSet(s0asc,s1asc,h,sz));
    RooGaussModel pdf1("pdf1","pdf1",dt,dtasc,sigma);
    RooAddPdf pdf("pdf","pdf",RooArgSet(g_otlr,pdf1),RooArgList(f_otlr));

    RooFitResult* fitRes = pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(h,sz)));
  } else{
    RooRealVar s0asc("s0asc","s0asc",Smn_asc,0.1,5.);   if(true_params)  s0asc.setConstant(kTRUE);
    RooRealVar s1asc("s1asc","s1asc",Stl_asc,0.,27.);   if(true_params || SVD == 1) s1asc.setConstant(kTRUE);
    RooRealVar f_tail("f_tail","f_tail",ftl_asc,0.,1.); if(true_params || SVD == 1) f_tail.setConstant(kTRUE);

    RooFormulaVar sigma1("sigma1","sigma1","@0*@1",RooArgList(s0asc,sz));
    RooFormulaVar sigma2("sigma2","sigma2","@0*@1",RooArgList(s1asc,sz));
    RooGaussModel G1("G1","G1",dt,dtasc,sigma1);
    RooGaussModel G2("G2","G2",dt,dtasc,sigma2);

//    RooAddPdf gascOm("gascOm","gascOm",RooArgSet(G2,G1),RooArgList(f_tail));
    if(SVD == 1) RooAddPdf pdf("pdf","pdf",RooArgList(g_otlr,G1),RooArgSet(f_otlr));
    else         RooAddPdf pdf("pdf","pdf",RooArgList(g_otlr,G2,G1),RooArgSet(f_otlr,f_tail));

    RooFitResult* fitRes = pdf.fitTo(ds,Verbose(),ConditionalObservables(RooArgSet(sz)));
  }

///////////
// Plots //
///////////
  RooPlot* dzFrameSig = dt.frame();
  ds.plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
//  if(single){
    pdf.plotOn(dzFrameSig,Components(g_otlr),LineStyle(kDashed),LineColor(kBlue),ProjWData(ds));
//    if(SVD == 2){
//      pdf.plotOn(dzFrameSig,Components(G2),LineStyle(kDashed),LineColor(kBlue),ProjWData(ds));
//    }
//  }
  pdf.plotOn(dzFrameSig,LineColor(kBlue),ProjWData(ds));

  RooHist* hdzpull = dzFrameSig->pullHist();
  RooPlot* dzPull = dt.frame(Title("#deltat pull distribution"));
  dzPull->addPlotable(hdzpull,"P");
  dzPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* dz_cm_sig = new TCanvas("#deltat_{asc} (ps)","#deltat_{asc} (ps)",700,800);
  dz_cm_sig->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd(); pad1->SetLeftMargin(0.15); pad1->SetFillColor(0);
  ds.statOn(dzFrameSig,Label(label.c_str()),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
//  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->SetGrid();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
//  pad1->SetLogy();
  dzFrameSig->GetXaxis()->SetTitleOffset(0.75);
  dzFrameSig->GetXaxis()->SetLabelSize(0.05);
  dzFrameSig->GetXaxis()->SetTitleSize(0.06);
  dzFrameSig->GetYaxis()->SetTitleOffset(1.6);
  dzFrameSig->Draw();

  stringstream out;
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

