#ifndef DEFIT_CPP
#define DEFIT_CPP

using namespace RooFit;

const double bdtg_min =-0.44;
const double deMin    =-0.12;
const double deMax    = 0.3;
const double mbcMin   = 5.23;
const double mbcMax   = 5.29;
const double de_min   =-0.035;
const double de_max   = 0.035;
const double atckpi_pi_max = 0.8;
const double cm2ps = 78.48566945838871754705;
double f_bkg     = 0.10;
double f_bkg_err = 0.01;

RooRealVar*  de;
RooRealVar*  mbc;
RooCategory* bin;

TChain* tree;
TTree* good_tree;
TTree* sb_tree;

double K_sb[8],K_sb_err[8];
double Kb_sb[8], Kb_sb_err[8];
double K[8],K_err[8],K_err_npq[8],K_err_sb[8],K_err_bkg[8];
double Kb[8],Kb_err[8],Kb_err_npq[8],Kb_err_sb[8],Kb_err_bkg[8];

bool sig_flag = false;
bool sb_flag  = false;
bool tree_flag = false;

void SetSignalTree(const int svd, const bool single_track_vertex){
  good_tree = GetGoodTree(svd,single_track_vertex,false);
  sig_flag = true;
}

void SetSidebandTree(const int svd, const bool single_track_vertex){
  sb_tree = GetGoodTree(svd,single_track_vertex,true);
  sb_flag = true;
}

TTree* GetGoodTree(const int svd, const bool single_track_vertex, const bool sideband){
  stringstream out;
  out.str("");
  if(single_track_vertex) out << "good_icpv_sgl && abs(dz_d0)<" << 70./cm2ps;
  else                    out << "good_icpv_mlt && abs(dz)<" << 70./cm2ps;
  if(svd == 1) out << " && exp<30";
  if(svd == 2) out << " && exp>30";
  out << " && bdtg>" << bdtg_min;
  out << " && atckpi_pi<" << atckpi_pi_max;
  if(sideband) out << " && (mbc>5.23 && mbc<5.26 || mbc>5.26 && de>0.1)";
  else         out << " && mbc>5.27 && mbc<5.289";
  out << " && bin != 0";
  TTree* new_tree = tree->CopyTree(out.str().c_str());
  return new_tree;
}

RooDataSet* GetGoodDS(const int svd, const bool single_track_vertex, const bool sideband){
  SetSignalTree(svd,single_track_vertex);
  RooArgSet argset;

  mbc = new RooRealVar("mbc","M_{bc}",mbcMin,mbcMax,"GeV");
  de  = new RooRealVar("de","#DeltaE",deMin,deMax,"GeV");
   de->setRange("Signal",de_min,de_max);
  mbc->setRange("Sideband1",5.23,5.26);
  mbc->setRange("Sideband2",5.26,5.29);
   de->setRange("Sideband2",0.1,0.3);

  argset.add(*de);
  argset.add(*mbc);

  RooDataSet* ds = new RooDataSet("ds","ds",good_tree,argset,"mbc>0||mbc<=0");
  return ds;
}

void deFit(const int mode = 0, const int svd = 2, const bool only_d0 = false){
  tree = new TChain("TEvent");

  const bool cComb = false;
  const bool cSig = true;
  const bool simple_peak = false;

  if(mode){
    tree->Add("FIL_b2dpi_charged_v2_0_10.root");
    tree->Add("FIL_b2dpi_charm_0_v2_10.root");
  } else{
    tree->Add("FIL_b2dpi_data_v2.root");
  }
  tree_flag = true;

  RooDataSet* ds = GetGoodDS(svd,only_d0,false);

  RooRealVar de0DK("de0DK","de0DK",-0.049,-0.055,-0.040); de0DK.setConstant(kTRUE);
  RooRealVar sDK("sDK","sDK",0.017,0.013,0.016);          sDK.setConstant(kTRUE);
  RooGaussian gDK("gDK","gDK",*de,de0DK,sDK);
  RooRealVar Nrho("NDK","NDK",50,0.,5000);// Nrho.setConstant(kTRUE);

  RooRealVar c1("c1","c1",-6.96922e-01,-10.,10.); if(cComb) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",1.72017e-01,-10.,10.); if(cComb) c2.setConstant(kTRUE);
  RooChebychev pdf_comb("pdf_comb","pdf_comb",*de,RooArgSet(c1,c2));
  RooRealVar Ncmb("NComb","NComb",10000,0.,50000);

  RooRealVar de0("de0","de0",0.,-0.005,0.005);
  RooRealVar s("s","s",1.19644e-02,0.010,0.015);
  if(simple_peak){
    RooGaussian pdf_sig("pdf_sig","pdf_sig",*de,de0,s);
  } else{
    RooGaussian g1("g1","g1",*de,de0,s);

    RooRealVar nl("nl","nl",7.78037e+00,0.,100.); if(cSig) nl.setConstant(kTRUE);
    RooRealVar alphal("alphal","alphal",-1,-10.,10.); alphal.setConstant(kTRUE);

    RooRealVar nr("nr","nr",1.29892e+01,0.,100.); if(cSig) nr.setConstant(kTRUE);
    RooRealVar alphar("alphar","alphar",1,-10.,10.); alphar.setConstant(kTRUE);

    RooCBShape CBl("CBl","CBl",*de,de0,s,alphal,nl);
    RooCBShape CBr("CBr","CBr",*de,de0,s,alphar,nr);

    RooRealVar fCBl("fCBl","fCBl",2.22046e-01,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
    RooRealVar fCBr("fCBr","fCBr",2.21964e-01,0.,1.); if(cSig) fCBr.setConstant(kTRUE);

    RooAddPdf pdf_sig("pdf_sig","pdf_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  }

  RooRealVar Nsig("NSig","NSig",20000,0.,25000);
  RooAddPdf pdf("pdf","pdf",RooArgSet(gDK,pdf_comb,pdf_sig),RooArgList(Nrho,Ncmb,Nsig));

  pdf.fitTo(*ds,Verbose(),Timer(true));


  RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(*de),NormSet(RooArgSet(*de)),Range("Signal"));
  RooAbsReal* intRho  = gDK.createIntegral(RooArgSet(*de),NormSet(RooArgSet(*de)),Range("Signal"));
  RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(*de),NormSet(RooArgSet(*de)),Range("Signal"));
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
    pdf_sig_val = Nsig.getVal()*pdf_sig.getVal(*de);
    pdf_DK_val = Nrho.getVal()*gDK.getVal(*de);
    pdf_smooth_val = Ncmb.getVal()*pdf_comb.getVal(*de);
    sig_frac = pdf_sig_val/(pdf_sig_val+pdf_DK_val+pdf_smooth_val);
    ofile << de->getVal() << " " << sig_frac << endl;
  }
  ofile.close();
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de->frame();
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(gDK),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_comb),LineStyle(kDashed));
  pdf.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de->frame(Title("#Delta E pull distribution"));
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
  f_bkg = 1.-purity;
  f_bkg_err = purity_err;
  return;
}

int DrawDalitzPlots(void){
  tree = new TChain("TEvent");
  tree->Add("FIL_b2dpi_data_v2.root");

  SetSignalTree(0,false);
  SetSidebandTree(0,false);

  TCanvas* sig_dp = new TCanvas("sig_dp","sig_dp",600,600);
  sig_dp->cd();
  stringstream out;
  out.str("");
  out << "bin != 0 && de>" << de_min << " && de<" << de_max;
  good_tree->Draw("mm:mp",out.str().c_str());
//  for(int i=0; i<8; i++){
//    good_tree->SetMarkerColor(i+1);
//    out.str("");
//    out << "abs(bin) == " << i+1;
//    out << " && de>" << de_min << " && de<" << de_max;
//    good_tree->Draw("mm:mp",out.str().c_str(),"same");
//  }
  sig_dp->Print("sig_dp.root");

  TCanvas* sb_dp = new TCanvas("sb_dp","sb_dp",600,600);
  sb_dp->cd();
  out.str("");
  out << "bin != 0";// && de>" << de_min << " && de<" << de_max;
  sb_tree->Draw("mm:mp",out.str().c_str());
//  for(int i=0; i<8; i++){
//    sb_tree->SetMarkerColor(i+1);
//    out.str("");
//    out << "abs(bin) == " << i+1;
////    out << " && de>" << de_min << " && de<" << de_max;
//    sb_tree->Draw("mm:mp",out.str().c_str(),"same");
//  }
  sb_dp->Print("sb_dp.root");

  return 0;
}

void calc_sideband_K(void){
  if(!tree_flag){
    tree = new TChain("TEvent");
    tree->Add("FIL_b2dpi_data_v2.root");
    tree_flag = true;
  }

  int rawK[8];
  int rawKb[8];

  if(!sb_flag) SetSidebandTree(0,false);
  stringstream out;
  const int NTot = sb_tree->GetEntries();
  cout << NTot << " events" << endl;
  for(int i=0; i<8; i++){
    out.str("");
    out << "flv*bin == " << i+1;
    rawK[i] = sb_tree->Draw("mm:mp",out.str().c_str());
    out.str("");
    out << "flv*bin == " << -(i+1);
    rawKb[i] = sb_tree->Draw("mm:mp",out.str().c_str());
    cout << " bin " << i+1 << ": Nsb = " << rawK[i] << ", anti-Nsb = " << rawKb[i] << endl;
    K_sb[i]      = (double)rawK[i]/NTot;
    Kb_sb[i]     = (double)rawKb[i]/NTot;
    K_sb_err[i]  = sqrt(K_sb[i]*(1.-K_sb[i])/NTot);
    Kb_sb_err[i] = sqrt(Kb_sb[i]*(1.-Kb_sb[i])/NTot);
  }
  cout << "Sideband K values:" << endl;
  for(int i=0; i<8; i++){
    cout << " bin " << i+1 << ": K_sb = " << K_sb[i] << " +- " << K_sb_err[i];
    cout << ", Kb_sb = " << Kb_sb[i] << " +- " << Kb_sb_err[i] << endl;
  }
  return;
}

void calc_K(void){
  deFit(0,0);
  calc_sideband_K();

  int rawK[8];
  int rawKb[8];

  if(!sig_flag) SetSignalTree(0,false);
  stringstream out;
  out.str("");
  out << "de>" << de_min << " && de<" << de_max << " && bin != 0";
  const string pre_cut = out.str();
  const int NTot = good_tree->Draw("mp:mm",pre_cut.c_str());
  const int NSig = NTot*(1.-f_bkg);
  cout << NTot << " events" << endl;
  for(int i=0; i<8; i++){
    out.str("");
    out << pre_cut << " && -flv*bin == " << i+1;
    rawK[i] = good_tree->Draw("mm:mp",out.str().c_str());
    out.str("");
    out << pre_cut << " && -flv*bin == " << -(i+1);
    rawKb[i] = good_tree->Draw("mm:mp",out.str().c_str());
    cout << " bin " << i+1 << ": N = " << rawK[i] << ", anti-N = " << rawKb[i] << endl;
    K[i]      = (double)(rawK[i]-NTot*f_bkg*K_sb[i])/NSig;
    Kb[i]     = (double)(rawKb[i]-NTot*f_bkg*Kb_sb[i])/NSig;
    K_err_npq[i] = sqrt(K[i]*(1.-K[i])/NSig);
    K_err_sb[i]  = f_bkg*K_sb_err[i]/(1.-f_bkg);
    K_err_bkg[i] = ((double)rawK[i]/NTot-K_sb[i]*(1.+2.*f_bkg))*f_bkg_err;
    K_err[i] = sqrt(K_err_npq[i]*K_err_npq[i]+K_err_sb[i]*K_err_sb[i]+K_err_bkg[i]*K_err_bkg[i]);
    Kb_err_npq[i] = sqrt(Kb[i]*(1.-Kb[i])/NSig);
    Kb_err_sb[i]  = f_bkg*Kb_sb_err[i]/(1.-f_bkg);
    Kb_err_bkg[i] = ((double)rawKb[i]/NTot-Kb_sb[i]*(1.+2.*f_bkg))*f_bkg_err;
    Kb_err[i] = sqrt(Kb_err_npq[i]*Kb_err_npq[i]+Kb_err_sb[i]*Kb_err_sb[i]+Kb_err_bkg[i]*Kb_err_bkg[i]);
  }
  cout << "K values:" << endl;
  for(int i=0; i<8; i++){
    cout << " bin " << i+1 << ": K  = " << K[i] << " +- " << K_err_npq[i] << " +- " << K_err_sb[i] << " +- " << K_err_bkg[i] << " (" << K_err[i] << ")" << endl;
    cout << "        Kb = " << Kb[i] << " +- " << Kb_err_npq[i] << " +- " << Kb_err_sb[i] << " +- " << Kb_err_bkg[i] << " (" << Kb_err[i] << ")" << endl;
  }
  return;
}

#endif // DEFIT_CPP
