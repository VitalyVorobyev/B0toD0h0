using namespace RooFit;

void dzPullFit(const int EXP = -1){
  TFile *omegafile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
  TTree *omegatree = (TTree*)omegafile->Get("TEvent");

  const double pull_min = -10;
  const double pull_max =  10;
  const double z_min = -10;
  const double z_max =  10;

  RooArgSet argset;
  RooCategory exp("exp","exp");
  if(EXP>0){
    exp.defineType("expNo",EXP);
    argset.add(exp);
  }

  RooRealVar dz_pull_sig("dz_pull_sig","dz_pull_sig",pull_min,pull_max); argset.add(dz_pull_sig);
//  RooRealVar z_sig("z_sig","z_sig",z_min,z_max,"mm"); argset.add(z_sig);
  RooRealVar sz_sig("sz_sig","dz_sig",0.01,0.15);     argset.add(sz_sig);
  RooRealVar chisq_z_sig("chisq_z_sig","chisq_z_sig",0.,200); argset.add(chisq_z_sig);

  RooDataSet omegads("omegads","omegads",omegatree,argset,"(sz_sig>0 || sz_sig<0) && (dz_pull_sig>0 || dz_pull_sig<0)");
  omegads.Print();

  RooRealVar dz01sigOm("dz01sigOm","dz01sigOm",0,z_min,z_max,"mm");// dz01sigOm.setConstant(kTRUE);
  RooRealVar dz02sigOm("dz02sigOm","dz02sigOm",0,z_min,z_max,"mm");// dz01sigOm.setConstant(kTRUE);
  RooRealVar s01sigOm("s01sigOm","s01sigOm",2.7,0.1,5.);
  RooRealVar s02sigOm("s02sigOm","s02sigOm",6.,0.1,25.);
  RooRealVar fsigOm("fsigOm","fsigOm",0.9,0.,1.);
//  RooGaussModel g1sigOm("g1sigOm","g1sigOm",dz_pull_sig,dz01sigOm,s01sigOm,sz_sig);
//  RooGaussModel g2sigOm("g2sigOm","g2sigOm",dz_pull_sig,dz01sigOm,s02sigOm,sz_sig);
  RooGaussian g1sigOm("g1sigOm","g1sigOm",dz_pull_sig,dz01sigOm,s01sigOm);
  RooGaussian g2sigOm("g2sigOm","g2sigOm",dz_pull_sig,dz01sigOm,s02sigOm);
  RooAddPdf pdfsigOm("pdfsigOm","pdfsigOm",RooArgList(g1sigOm,g2sigOm),RooArgSet(fsigOm));

  RooFitResult* fitRes = pdfsigOm.fitTo(omegads,Verbose(),ConditionalObservables(sz_sig));

///////////
// Plots //
///////////
  RooPlot* dzFrameSig = dz_pull_sig.frame();
  omegads.plotOn(dzFrameSig,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdfsigOm.plotOn(dzFrameSig,LineColor(kRed),ProjWData(omegads));
  TCanvas* dz_cm_sig = new TCanvas("#deltaz_{sig}, cm","#deltaz_{sig}, mm",700,500);
  dz_cm_sig->cd();
  omegads.statOn(dzFrameSig,Layout(0.6,0.98,0.9),What("MR"),Label("#omega^{0}"),Format("NE",AutoPrecision(1)));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  dzFrameSig->GetXaxis()->SetTitleOffset(0.75);
  dzFrameSig->GetXaxis()->SetLabelSize(0.05);
  dzFrameSig->GetXaxis()->SetTitleSize(0.06);
  dzFrameSig->GetYaxis()->SetTitleOffset(1.6);
  dzFrameSig->Draw();
  dz_cm_sig->Update();
  stringstream out;
  out.str("");
  out << "zSigPullFit_Exp" << EXP << ".png";
  dz_cm_sig->Print(out.str().c_str());
  out.str("");
  out << "zSigPullFit_Exp" << EXP << ".root";
  dz_cm_sig->Print(out.str().c_str());

}

