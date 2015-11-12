using namespace RooFit;

void dzRooFit(const int _mode = 1, const int _h0mode = 10,const int EXP = -1){
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

  const double dz_min = -1.5;
  const double dz_max =  1.5;

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
//  argset.add(nptag);

  RooRealVar dz("dz","dz",dz_min,dz_max,"mm");    argset.add(dz);
  RooRealVar sz_asc("sz_asc","dz_asc",0.01,0.15); argset.add(sz_asc);
  RooRealVar sz_sig("sz_sig","dz_sig",0.01,0.15); argset.add(sz_sig);
  RooRealVar chisq_z_asc("chisq_z_asc","chisq_z_asc",0.,200); argset.add(chisq_z_asc);
  RooRealVar chisq_z_sig("chisq_z_sig","chisq_z_sig",0.,200); argset.add(chisq_z_sig);

//  RooDataSet pids("pids","pids",pitree,argset);
  RooDataSet omegads("omegads","omegads",omegatree,argset,"(dz>0 || dz<0) && sz_asc>0 && sz_sig>0 && chisq_z_asc>0 && chisq_z_sig>0");
//  pids.Print();
  omegads.Print();

  RooRealVar tau("tau","tau",0.196,0.01,1,"mm"); tau.setConstant(kTRUE);
  RooRealVar dz0Om("dz0Om","dz0Om",0,-3.,3.,"mm"); dz0Om.setConstant(kTRUE);
  RooRealVar s01sigOm("s01sigOm","s01sigOm",2,0.1,5.);
  RooRealVar s02sigOm("s02sigOm","s02sigOm",13,0.1,27.);
  RooRealVar fsigOm("fsigOm","fsigOm",0.9,0.,1.);

//  RooGaussModel g1Om("g1Om","g1Om",dz,dz0Om,s01Om,sz);
//  RooGaussModel g2Om("g2Om","g2Om",dz,dz0Om,s02Om,sz);
  RooDecay dec1SigOm("dec1SigOm","dec1SigOm",dz,tau,g1sigOm,RooDecay::DoubleSided);
  RooDecay dec2SigOm("dec2SigOm","dec2SigOm",dz,tau,g2sigOm,RooDecay::DoubleSided);
  RooAddPdf pdfsigOm("pdfsigOm","pdfsigOm",RooArgList(dec1SigOm,dec2SigOm),RooArgSet(fsigOm));
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
  out << "pics/zAscFit_Exp" << EXP << ".png";
  z_cm_sig->Print(out.str().c_str());
  out.str("");
  out << "pics/zAscFit_Exp" << EXP << ".root";
  z_cm_sig->Print(out.str().c_str());

  return;
}

