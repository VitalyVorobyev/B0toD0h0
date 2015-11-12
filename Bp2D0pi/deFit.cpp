#ifndef DEFIT_CPP
#define DEFIT_CPP

using namespace RooFit;

void deFit(const int mode = 0, const int svd = 2, const bool only_d0 = false){
  TChain* tree = new TChain("TEvent");

  const double bdtg_min = -0.44;
  const double deMin = -0.12;
  const double deMax = 0.3;
  const double mbcMin = 5.272;
  const double mbcMax = 5.287;
//  const double sz_sig_max = 0.2;
//  const double sz_asc_max = 0.2;
//  const double chisq_sig_max = only_d0 ? 10000 : 50;
//  const double chisq_asc_max = 50;
  const double atckpi_pi_max = 0.8;

  const double cm2ps = 78.48566945838871754705;

  const double de_min = -0.035;
  const double de_max =  0.035;

  const bool cComb = false;
  const bool cSig = true;
  const bool simple_peak = false;

  if(mode){
    tree->Add("FIL_b2dpi_charged_v2_0_10.root");
    tree->Add("FIL_b2dpi_charm_0_v2_10.root");
  } else{
    tree->Add("FIL_b2dpi_data_v2.root");
  }

  RooArgSet argset;

  RooCategory exp("exp","exp");
  if(svd == 1){
    exp.defineType("7",7);
    exp.defineType("9",9);
    exp.defineType("11",11);
    exp.defineType("13",13);
    exp.defineType("15",15);
    exp.defineType("17",17);
    exp.defineType("19",19);
    exp.defineType("21",21);
    exp.defineType("23",23);
    exp.defineType("25",25);
    exp.defineType("27",27);
  } else{
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

  RooCategory ndf_z_asc("ndf_z_asc","ndf_z_asc");
  ndf_z_asc.defineType("0",0);
  ndf_z_asc.defineType("2",2);
  ndf_z_asc.defineType("4",4);
  ndf_z_asc.defineType("6",6);
  ndf_z_asc.defineType("8",8);
  ndf_z_asc.defineType("10",10);
  ndf_z_asc.defineType("12",12);
  ndf_z_asc.defineType("14",14);
  ndf_z_asc.defineType("16",16);
  ndf_z_asc.defineType("18",18);
  argset.add(ndf_z_asc);

  RooCategory good_icpv_mlt("good_icpv_mlt","good_icpv_mlt");
  good_icpv_mlt.defineType("good",1);

  RooCategory good_icpv_sgl("good_icpv_sgl","good_icpv_sgl");
  good_icpv_sgl.defineType("good",1);

  if(only_d0){
    argset.add(good_icpv_sgl);
  } else{
    argset.add(good_icpv_mlt);
  }

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  if(!only_d0){ RooRealVar dz("dz","#Deltaz",-70./cm2ps,70./cm2ps,"cm"); argset.add(dz);}
  else        { RooRealVar dz("dz_d0","#Deltaz",-70./cm2ps,70./cm2ps,"cm"); argset.add(dz);}
  RooRealVar bdtg("bdtg","bdtg",bdtg_min,1.); argset.add(bdtg);
  RooRealVar atckpi_pi("atckpi_pi","atckpi_pi",0.,atckpi_pi_max); argset.add(atckpi_pi);
//  if(!only_d0){ RooRealVar sz_sig("sz_sig","sz_sig",0.,sz_sig_max,"mm"); argset.add(sz_sig);}
//  else        { RooRealVar sz_sig("sz_sig_d0","#sigma_{z}^{sig}",0.,sz_sig_max,"mm"); argset.add(sz_sig);}
//  RooRealVar sz_asc("sz_asc","sz_asc",0.,sz_asc_max,"mm"); argset.add(sz_asc);
//  if(!only_d0 || true){ RooRealVar chisq_z_sig("chisq_z_sig","chisq_z_sig",0.,chisq_sig_max); argset.add(chisq_z_sig);}
//  RooRealVar chisq_z_asc("chisq_z_asc","chisq_z_asc",0.,chisq_asc_max); argset.add(chisq_z_asc);

//  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0 && (ndf_z_asc == 0 || 1.*chisq_z_asc/(ndf_z_asc+0.001)<10)");
  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");

  RooRealVar de0DK("de0DK","de0DK",-0.049,-0.055,-0.040); de0DK.setConstant(kTRUE);
  RooRealVar sDK("sDK","sDK",0.017,0.013,0.016);          sDK.setConstant(kTRUE);
  RooGaussian gDK("gDK","gDK",de,de0DK,sDK);
  RooRealVar Nrho("NDK","NDK",50,0.,5000);// Nrho.setConstant(kTRUE);

  RooRealVar c1("c1","c1",-6.96922e-01,-10.,10.); if(cComb) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",1.72017e-01,-10.,10.); if(cComb) c2.setConstant(kTRUE);
  RooChebychev pdf_comb("pdf_comb","pdf_comb",de,RooArgSet(c1,c2));
  RooRealVar Ncmb("NComb","NComb",10000,0.,50000);

  RooRealVar de0("de0","de0",0.,-0.005,0.005);
  RooRealVar s("s","s",1.19644e-02,0.010,0.015);
  if(simple_peak){
    RooGaussian pdf_sig("pdf_sig","pdf_sig",de,de0,s);
  } else{
    RooGaussian g1("g1","g1",de,de0,s);

//    RooRealVar nl("nl","nl",4.93610e+00,0.,100.); if(cSig) nl.setConstant(kTRUE);
    RooRealVar nl("nl","nl",7.78037e+00,0.,100.); if(cSig) nl.setConstant(kTRUE);
    RooRealVar alphal("alphal","alphal",-1,-10.,10.); alphal.setConstant(kTRUE);

//    RooRealVar nr("nr","nr",4.91073e+00,0.,100.); if(cSig) nr.setConstant(kTRUE);
    RooRealVar nr("nr","nr",1.29892e+01,0.,100.); if(cSig) nr.setConstant(kTRUE);
    RooRealVar alphar("alphar","alphar",1,-10.,10.); alphar.setConstant(kTRUE);

    RooCBShape CBl("CBl","CBl",de,de0,s,alphal,nl);
    RooCBShape CBr("CBr","CBr",de,de0,s,alphar,nr);

//    RooRealVar fCBl("fCBl","fCBl",2.40571e-01,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
//    RooRealVar fCBr("fCBr","fCBr",2.20385e-01,0.,1.); if(cSig) fCBr.setConstant(kTRUE);
    RooRealVar fCBl("fCBl","fCBl",2.22046e-01,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
    RooRealVar fCBr("fCBr","fCBr",2.21964e-01,0.,1.); if(cSig) fCBr.setConstant(kTRUE);

    RooAddPdf pdf_sig("pdf_sig","pdf_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  }

  RooRealVar Nsig("NSig","NSig",20000,0.,25000);
  RooAddPdf pdf("pdf","pdf",RooArgSet(gDK,pdf_comb,pdf_sig),RooArgList(Nrho,Ncmb,Nsig));

  pdf.fitTo(ds,Verbose(),Timer(true));


  RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
  RooAbsReal* intRho  = gDK.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
  RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
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

  double sig_frac;
  double pdf_sig_val;
  double pdf_DK_val;
  double pdf_smooth_val;
  fstream ofile("de_sig_fraction.txt",fstream::out);
  for(int i=0; i<1000; i++){
    const double dde = 0.2/1000;
    de.setVal(-0.1+(i+0.5)*dde);
    pdf_sig_val = Nsig.getVal()*pdf_sig.getVal(de);
    pdf_DK_val = Nrho.getVal()*gDK.getVal(de);
    pdf_smooth_val = Ncmb.getVal()*pdf_comb.getVal(de);
    sig_frac = pdf_sig_val/(pdf_sig_val+pdf_DK_val+pdf_smooth_val);
    ofile << de.getVal() << " " << sig_frac << endl;
  }
  ofile.close();
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(gDK),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_comb),LineStyle(kDashed));
  pdf.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("Delta E","Delta E",600,700);
  cm->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);
  pad3->SetGrid();

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.04);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  const int height = svd == 2 ? 500 : 120;
  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,height);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,height);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  stringstream out1;
  TPaveText *pt = new TPaveText(0.4,0.65,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out1.str("");
  out1 << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out1.str().c_str());
  out1.str("");
  out1 << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err_total+0.5);
  pt->AddText(out1.str().c_str());
  out1.str("");
  out1 << "Purity: " << std::fixed << std::setprecision(2) << purity*100. << " #pm " << purity_err*100;
  pt->AddText(out1.str().c_str());
  pt->Draw();

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

//  cout << "Nsig = " << nsig <<" +- " << nsig_err << endl;
//  cout << "NDK  = " << nrho <<" +- " << nrho_err << endl;
//  cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << endl;
//  cout << "Purity = " << purity << " +- " << purity_err << endl;

  cout << "Nsig = " << nsig <<" +- " << nsig_err << " +- " << nsig_err_npq << " (" << nsig_err_total << ")" << endl;
  cout << "NDK  = " << nrho <<" +- " << nrho_err << " +- " << nrho_err_npq << " (" << nrho_err_total << ")" << endl;
  cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << " +- " << ncmb_err_npq << " (" << ncmb_err_total << ")" << endl;
  cout << "Pury = " << purity << " +- " << purity_err << endl;
}

#endif // DEFIT_CPP
