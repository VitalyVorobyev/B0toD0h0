#ifndef DEFIT_CPP
#define DEFIT_CPP

using namespace RooFit;

void deSigFit(const int svd = 2){
  TChain* tree = new TChain("TEvent");

  const double bdtg_min = -0.44;
  const double deMin = -0.1;
  const double deMax = 0.3;
  const double mbcMin = 5.272;
  const double mbcMax = 5.287;
  const double sz_sig_max = 0.2;
  const double sz_asc_max = 0.2;
  const double chisq_sig_max = 10;
  const double chisq_asc_max = 10;
  const double atckpi_pi_max = 10;

  const double de_min = -0.05;
  const double de_max =  0.05;

  const bool cComb = false;
  const bool cSig = false;
  tree->Add("FIL_b2dpi_charged_v2_0_10.root");

  RooArgSet argset;
  RooCategory bpf("bpf","bpf");
  bpf.defineType("signal",1);
  argset.add(bpf);

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
//  de.setRange("Signal",de_min,de_max);
//  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
//  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar bdtg("bdtg","bdtg",bdtg_min,1.); argset.add(bdtg);
  RooRealVar atckpi_pi("atckpi_pi","atckpi_pi",0.,atckpi_pi_max); argset.add(atckpi_pi);
  RooRealVar sz_sig("sz_sig","sz_sig",0.,sz_sig_max,"mm"); argset.add(sz_sig);
  RooRealVar sz_asc("sz_asc","sz_asc",0.,sz_asc_max,"mm"); argset.add(sz_asc);
  RooRealVar chisq_z_sig("chisq_z_sig","chisq_z_sig",0.,chisq_sig_max); argset.add(chisq_z_sig);
  RooRealVar chisq_z_asc("chisq_z_asc","chisq_z_asc",0.,chisq_asc_max); argset.add(chisq_z_asc);

  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");

  RooRealVar de0("de0","de0",1.04590e-03,-0.1,0.1); if(cSig) de0.setConstant(kTRUE);
  RooRealVar s1("s1","s1",9.53894e-03,0.,0.5); if(cSig) s1.setConstant(kTRUE);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar sCBl("sCBl","sCBl",1.45483e-02,0.,0.5); if(cSig) sCBl.setConstant(kTRUE);
  RooRealVar nl("nl","nl",5.63847e+00,0.,100.); if(cSig) nl.setConstant(kTRUE);
  RooRealVar alphal("alphal","alphal",-1,-10.,10.); alphal.setConstant(kTRUE);

  RooRealVar sCBr("sCBr","sCBr",1.57529e-02,0.,0.5); if(cSig) sCBr.setConstant(kTRUE);
  RooRealVar nr("nr","nr",4.79889e+00,0.,100.); if(cSig) nr.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",1,-10.,10.); alphar.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",de,de0,s1,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,de0,s1,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",2.03618e-01,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
  RooRealVar fCBr("fCBr","fCBr",1.42928e-01,0.,1.); if(cSig) fCBr.setConstant(kTRUE);

  RooAddPdf pdf("pdf","pdf",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
//    RooAddPdf pdf("pdf","pdf",RooArgList(CBl,CBr),RooArgSet(fCBl));

  pdf.fitTo(ds,Verbose(),Timer(true));

  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(CBl),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(CBr),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(g1),LineStyle(kDashed));
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

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.04);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,50);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,50);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

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
}

#endif // DEFIT_CPP

