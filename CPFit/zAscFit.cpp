using namespace RooFit;

void zAscFit(void){
  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const double init_scale = 100.;
  const double z_min = -0.5;
  const double z_max =  0.5;

  RooArgSet argset;
  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig",1);
  b0f.defineType("badpi0",5);
  argset.add(b0f);

  RooRealVar z_asc("z_asc","z_asc",z_min,z_max,"mm"); argset.add(z_asc);
  RooRealVar sz_asc("sz_asc","dz_asc",0.,0.2); argset.add(sz_asc);
  RooRealVar chisq_z_asc("chisq_z_asc","chisq_z_asc",0.,200); argset.add(chisq_z_asc);

  RooDataSet omegads("omegads","omegads",omegatree,argset);
  omegads.Print();
  RooDataSet* szds = omegads.reduce(RooArgSet(sz_asc));
  RooDataHist* szdh = szds->binnedClone();
  RooHistPdf pdfsz("pdfsz","pdfsz",sz_asc,*szdh);

  RooRealVar tau("tau","tau",0.1,0.01,1,"mm");

  RooRealVar dz01ascOm("dz01ascOm","dz01ascOm",0,z_min,z_max,"mm"); dz01ascOm.setConstant(kTRUE);
//  RooRealVar dz02sigOm("dz02sigOm","dz02sigOm",0,z_min,z_max,"mm"); dz02sigOm.setConstant(kTRUE);
  RooRealVar s01ascOm("s01ascOm","s01ascOm",0.1*init_scale,0.,5.*init_scale,"mm");
//  RooRealVar s02sigOm("s02sigOm","s02sigOm",0.1*init_scale,0.,5.*init_scale,"mm");
//  RooFormulaVar szSigOm("szSigOm","100*sqrt(@0)",RooArgSet(sz_sig));
//  RooRealVar fsigOm("fsigOm","fsigOm",0.5,0.,1.);
  RooGaussModel g1ascOm("g1ascOm","g1ascOm",z_asc,dz01ascOm,s01ascOm,sz_asc);
//  RooGaussModel g2sigOm("g2sigOm","g2sigOm",z_sig,dz02sigOm,s02sigOm,sz_sig);
//  RooDecay dec1SigOm("dec1SigOm","dec1SigOm",z_sig,tau,g1sigOm,RooDecay::SingleSided);
//  RooDecay dec2SigOm("dec2SigOm","dec2SigOm",z_sig,tau,g2sigOm,RooDecay::SingleSided);
//  RooAddPdf pdfsigOm("pdfsigOm","pdfsigOm",RooArgList(g1sigOm,g2sigOm),RooArgSet(fsigOm));
  RooDecay pdfascOm("pdfascOm","pdpascOm",z_asc,tau,g1ascOm,RooDecay::SingleSided);

  RooProdPdf pdf("pdf","pdf",pdfsz,Conditional(pdfascOm,sz_asc));

  RooFitResult* fitRes = pdf.fitTo(omegads,Verbose());
//  fitRes->Print();

///////////
// Plots //
///////////
  RooPlot* zFrameAsc = z_asc.frame();
//  pids.plotOn(zFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
//  pdfsigPi.plotOn(zFrameSig,LineColor(kBlue));
  omegads.plotOn(zFrameAsc,DataError(RooAbsData::SumW2),MarkerSize(1));
//  pdf.plotOn(zFrameSig,LineColor(kRed),ProjWData(omegads));
  pdf.plotOn(zFrameAsc,LineColor(kRed));

  TCanvas* z_cm_asc = new TCanvas("z_{asc}, cm","z_{asc}, mm",700,500);
  z_cm_asc->cd();
//  pids.statOn(zFrameSig,Label("#pi^{0}"),Layout(0.6,0.98,0.9),What("MR"),Format("NE",AutoPrecision(1)));
  omegads.statOn(zFrameAsc,Layout(0.6,0.98,0.7),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  zFrameAsc->GetXaxis()->SetTitleOffset(0.75);
  zFrameAsc->GetXaxis()->SetLabelSize(0.05);
  zFrameAsc->GetXaxis()->SetTitleSize(0.06);
  zFrameAsc->GetYaxis()->SetTitleOffset(1.6);
  zFrameAsc->Draw();
  z_cm_asc->Update();
//  z_cm_sig->Print("dz_sig_fit.png");
//  z_cm_sig->Print("dz_sig_fit.root");
//  fitRes->Print();

}
