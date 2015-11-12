#include "cuts.h"
using namespace RooFit;

void Signal_2d_fit_omega(const int _b0f=-1){
  // mode 1 -> pi0
  // mode 2 -> eta
  // mode 3 -> omega
  const bool projection_flag = true;
  const bool save_flag       = true;

  TFile *ifile;
  ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");

  const bool remove_left_CB_flag = false;
  bool remove_right_CB_flag  = false;
  if(_b0f == 5) remove_right_CB_flag = true;

  TTree *tree = (TTree*)ifile->Get("TEvent");
  RooArgSet argset;

  string cuts;
  RooCategory b0f("b0f","b0f");
  if(_b0f == 1 || _b0f == -1) b0f.defineType("signal",1);
  if(_b0f == 1 || _b0f == -1) b0f.defineType("fsr",10);
  if(_b0f == 5 || _b0f == -1){
    b0f.defineType("bad_pi0",5);
    b0f.defineType("bad_pi0_m11",-11);
    b0f.defineType("bad_pi0_4",4);
  }
  argset.add(b0f);

  RooCategory mode("mode","mode");
  mode.defineType("omega",3);
  argset.add(mode);

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin  = -0.15;
  const double deMax  =  0.30;

  double de_sig_min,de_sig_max;
  de_sig_min = de_min_omega;
  de_sig_max = de_max_omega;

  double BDTG_MIN = bdtg_cut_omega;
  double BDTG_MAX = 1;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_sig_min,de_sig_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_sig_min,de_sig_max);

  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  RooDataSet ds("ds","ds",tree,argset);
  RooDataSet* ds0 = ds.reduce(RooArgSet(de,mbc));
  ds.Print();

  const int _mode = 3;
  const int _h0mode = 20;
  RooRealVar mbc0("mbc0","mbc0",get_mbc0(_mode,_h0mode,_b0f),5.26,5.30);     if(cSIG && false) mbc0.setConstant(kTRUE);
  ////////////
  // de pdf //
  ////////////
  if(false){
  RooRealVar de0("de0","de0",get_de0(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG) de0.setConstant(kTRUE);
  RooRealVar s1("s1","s1",get_s1(_mode,_h0mode,_b0f),0.,0.5); if(cSIG) s1.setConstant(kTRUE);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar deCBl("deCBl","deCBl",get_deCBl(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG) deCBl.setConstant(kTRUE);
  RooRealVar sCBl("sCBl","sCBl",get_sCBl(_mode,_h0mode,_b0f),0.,0.5); if(cSIG) sCBl.setConstant(kTRUE);
  RooRealVar nl("nl","nl",get_nl(_mode,_h0mode,_b0f),0.,100.); if(cSIG) nl.setConstant(kTRUE);
//  RooRealVar nl("nl","nl",1.,0.,100.); nl.setConstant(kTRUE);
  RooRealVar alphal("alphal","alphal",get_alphal(_mode,_h0mode,_b0f),-10.,10.); if(cSIG) alphal.setConstant(kTRUE);

  RooRealVar deCBr("deCBr","deCBr",get_deCBr(_mode,_h0mode,_b0f),-0.2,0.1); if(cSIG || remove_right_CB_flag) deCBr.setConstant(kTRUE);
  RooRealVar sCBr("sCBr","sCBr",get_sCBr(_mode,_h0mode,_b0f),0.,0.5); if(cSIG || remove_right_CB_flag) sCBr.setConstant(kTRUE);
  RooRealVar nr("nr","nr",get_nr(_mode,_h0mode,_b0f),0.,100.); if(cSIG || remove_right_CB_flag) nr.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",get_alphar(_mode,_h0mode,_b0f),-10.,10.); if(cSIG || remove_right_CB_flag) alphar.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",get_fCBl(_mode,_h0mode,_b0f),0.,1.); if(cSIG) fCBl.setConstant(kTRUE);
  RooRealVar fCBr("fCBr","fCBr",get_fCBr(_mode,_h0mode,_b0f),0.,1.); if(cSIG || remove_right_CB_flag) fCBr.setConstant(kTRUE);
  if(remove_right_CB_flag) fCBr.setVal(0);

  RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));
  }

  if(true){
  RooRealVar de0_201("de0_201","de0_201",-1.17124e-02,-0.1,0.1);       if(_b0f != 1) de0_201.setConstant(kTRUE);
  RooRealVar s1_201("s1_201","s1_201",1.88306e-02,0.,0.5); if(_b0f != 1) s1_201.setConstant(kTRUE);

  RooGaussian g1_201("g1_201","g1_201",de,de0_201,s1_201);

  RooRealVar deCBl_201("deCBl_201","deCBl_201",2.19044e-04,-0.1,0.1);     if(_b0f != 1) deCBl_201.setConstant(kTRUE);
  RooRealVar sCBl_201("sCBl_201","sCBl_201",1.28877e-02,0.,0.5);           if(_b0f != 1) sCBl_201.setConstant(kTRUE);
  RooRealVar nl_201("nl_201","nl_201",3.00711e+00,0.,100.);                  if(_b0f != 1) nl_201.setConstant(kTRUE);
  RooRealVar alphal_201("alphal_201","alphal_201",1.12531e+00,-10.,10.); if(_b0f != 1) alphal_201.setConstant(kTRUE);
  RooRealVar deCBr_201("deCBr_201","deCBr_201",2.80491e-02,-0.1,0.1);     if(_b0f != 1) deCBr_201.setConstant(kTRUE);
  RooRealVar sCBr_201("sCBr_201","sCBr_201",1.22731e-02,0.,0.5);           if(_b0f != 1) sCBr_201.setConstant(kTRUE);
  RooRealVar nr_201("nr_201","nr_201",7.47296e+00,0.,100.);                  if(_b0f != 1) nr_201.setConstant(kTRUE);
  RooRealVar alphar_201("alphar_201","alphar_201",-6.50995e-01,-10.,10.); if(_b0f != 1) alphar_201.setConstant(kTRUE);

  RooCBShape CBl_201("CBl_201","CBl_201",de,deCBl_201,sCBl_201,alphal_201,nl_201);
  RooCBShape CBr_201("CBr_201","CBr_201",de,deCBr_201,sCBr_201,alphar_201,nr_201);

  RooRealVar fCBl_201("fCBl_201","fCBl_201",7.29526e-01,0.,1.); if(_b0f != 1) fCBl_201.setConstant(kTRUE);
  RooRealVar fCBr_201("fCBr_201","fCBr_201",1.30198e-01,0.,1.); if(_b0f != 1) fCBr_201.setConstant(kTRUE);

  RooAddPdf pdf_de_201("pdf_de_201","pdf_de_201",RooArgList(CBl_201,CBr_201,g1_201),RooArgSet(fCBl_201,fCBr_201));

  RooRealVar de0_205("de0_205","de0_205",-1.55385e-01,-0.15,0.1);  if(_b0f != 5) de0_205.setConstant(kTRUE);
  RooRealVar s1_205("s1_205","s1_205",1.37438e-01,0.,0.5);         if(_b0f != 5) s1_205.setConstant(kTRUE);
  RooGaussian g1_205("g1_205","g1_205",de,de0_205,s1_205);

  RooRealVar deCBl_205("deCBl_205","deCBl_205",-2.53735e-02,-0.1,0.1);   if(_b0f != 5) deCBl_205.setConstant(kTRUE);
  RooRealVar sCBl_205("sCBl_205","sCBl_205",3.89017e-02,0.,0.5);         if(_b0f != 5) sCBl_205.setConstant(kTRUE);
  RooRealVar nl_205("nl_205","nl_205",9.88285e+01,0.,100.);              if(_b0f != 5) nl_205.setConstant(kTRUE);
  RooRealVar alphal_205("alphal_205","alphal_205",3.69154e-01,-10.,10.); if(_b0f != 5) alphal_205.setConstant(kTRUE);
  RooCBShape CBl_205("CBl_205","CBl_205",de,deCBl_205,sCBl_205,alphal_205,nl_205);

  RooRealVar fCBl_205("fCBl_205","fCBl_205",7.32842e-01,0.,1.);          if(_b0f != 5) fCBl_205.setConstant(kTRUE);

  RooAddPdf pdf_de_205("pdf_de_205","pdf_de_205",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));
//  RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(CBl_205,g1_205),RooArgSet(fCBl_205));

  RooRealVar f_201("f_201","f_201",get_f201(_mode,_h0mode),0.,1.);// if(cSIG) f_201.setConstant(kTRUE);
  RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(pdf_de_201,pdf_de_205),RooArgSet(f_201));
//  RooGaussian pdf_de1("pdf_de1","pdf_de1",de,de0_201,s1_201);
//  RooGaussian pdf_de2("pdf_de2","pdf_de2",de,deCBl_201,sCBl_201);
//  RooAddPdf pdf_de("pdf_de","pdf_de",RooArgList(pdf_de1,pdf_de2),RooArgSet(f_201));
  }

  /////////////
  // mbc pdf //
  /////////////
  const bool fflag = false;
  const bool cond_flag = false;//(b0f == 5);
  RooRealVar c1_mbc0("c1_mbc0","c1_mbc0",-0.02723,-0.1,0.); c1_mbc0.setConstant(kTRUE);
  RooRealVar c2_mbc0("c2_mbc0","c2_mbc0",0.1456,0.,1.);     c2_mbc0.setConstant(kTRUE);
  RooRealVar c3_mbc0("c3_mbc0","c3_mbc0",7.212,0.,10.);     c3_mbc0.setConstant(kTRUE);
  RooRealVar mbc0("mbc0","mbc0",5.28,5.26,5.30);             mbc0.setConstant(kTRUE);
  RooFormulaVar _mbc0("_mbc0","_mbc0","@0+@1*@4+@2*@4*@4+@3*@4*@4*@4",RooArgList(mbc0,c1_mbc0,c2_mbc0,c3_mbc0,de));

  RooRealVar c1_alpha("c1_alpha","c1_alpha",-0.5706,-2.,0.);  c1_alpha.setConstant(kTRUE);
  RooRealVar c2_alpha("c2_alpha","c2_alpha",60.86,10.,100.);      c2_alpha.setConstant(kTRUE);
  RooRealVar alpha("alpha","alpha",0.1004,0.,0.1);           alpha.setConstant(kTRUE);
  RooFormulaVar _alpha("_alpha","_alpha","@0+@1*@3+@2*@3*@3",RooArgList(alpha,c1_alpha,c2_alpha,de));

  RooRealVar c1_width("c1_width","c1_width",0.0002038,-2.,2.);  c1_width.setConstant(kTRUE);
  RooRealVar c2_width("c2_width","c2_width",0.299,-10.,10);     c2_width.setConstant(kTRUE);
  RooRealVar width("width","width",0.002606,0.,0.1);            width.setConstant(kTRUE);
  RooFormulaVar _width("_width","_width","@0+@1*@3+@2*@3*@3",RooArgList(width,c1_width,c2_width,de));

  RooNovosibirsk pdf_mbc("pdf_mbc","pdf_mbc",mbc,_mbc0,_width,_alpha);
//  RooGaussian pdf_mbc("pdf_mbc","pdf_mbc",mbc,_mbc0,_width);

  /////////
  // pdf //
  /////////
  RooProdPdf pdf("pdf","pdf",pdf_de,Conditional(pdf_mbc,mbc));
//  RooProdPdf pdf("pdf","pdf",RooArgSet(pdf_mbc,pdf_de,de));

  pdf_de.fitTo(ds,Verbose(),Timer(true));

  de0_201.setConstant(kTRUE);
  s1_201.setConstant(kTRUE);
  deCBl_201.setConstant(kTRUE);
  sCBl_201.setConstant(kTRUE);
  nl_201.setConstant(kTRUE);
  alphal_201.setConstant(kTRUE);
  deCBr_201.setConstant(kTRUE);
  sCBr_201.setConstant(kTRUE);
  nr_201.setConstant(kTRUE);
  alphar_201.setConstant(kTRUE);
  fCBl_201.setConstant(kTRUE);
  fCBr_201.setConstant(kTRUE);

  de0_205.setConstant(kTRUE);
  s1_205.setConstant(kTRUE);
  deCBl_205.setConstant(kTRUE);
  sCBl_205.setConstant(kTRUE);
  nl_205.setConstant(kTRUE);
  alphal_205.setConstant(kTRUE);
  fCBl_205.setConstant(kTRUE);

  f_201.setConstant(kTRUE);

  pdf.fitTo(ds,Verbose(),Timer(true));
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  if(projection_flag){
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
    pdf.plotOn(deFrame,LineWidth(2),LineColor(kGreen));
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
    pdf.plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  } else{
    ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(deFrame,LineWidth(2));
  }
  ds.statOn(deFrame,Layout(0.55,0.98,0.9));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E, Signal","#Delta E, Signal",600,700);
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

  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  TPaveText *pt = new TPaveText(0.6,0.6,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(out.str().c_str());
  pt->Draw();

  if(projection_flag){
  TLine *de_line_RIGHT = new TLine(de_sig_max,0,de_sig_max,350);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_sig_min,0,de_sig_min,350);
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();
  }

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
  out.str("");
  out << "../Reports/pics/de_sig_omega" << "_b0f" << _b0f;
  if(projection_flag) out << "_wproj";
  if(deMin>-0.2) out << "015";
  out << ".png";
  if(save_flag) cm->Print(out.str().c_str());

  // mbc //
  RooPlot* mbcFrame = mbc.frame();
  if(projection_flag){
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
    pdf.plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
    pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  } else{
    ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf.plotOn(mbcFrame,LineWidth(2));
  }

  ds.statOn(mbcFrame,Layout(0.2,0.68,0.9));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc.frame(Title("#Delta E pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("M_{bc}, Signal","M_{bc}, Signal",600,700);
  cmmbc->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);

  mbcFrame->GetXaxis()->SetTitleSize(0.05);
  mbcFrame->GetXaxis()->SetTitleOffset(0.85);
  mbcFrame->GetXaxis()->SetLabelSize(0.04);
  mbcFrame->GetYaxis()->SetTitleOffset(1.6);
  mbcFrame->Draw();

  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.3,0.6,0.68,0.7,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  if(projection_flag){
  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,800);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,800);
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();
  }

  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  mbcPull->SetMarkerSize(0.05);
  mbcPull->Draw();
  TLine *mbc_lineUP = new TLine(mbcMin,3,mbcMax,3);
  mbc_lineUP->SetLineColor(kBlue);
  mbc_lineUP->SetLineStyle(2);
  mbc_lineUP->Draw();
  TLine *mbc_line = new TLine(mbcMin,0,mbcMax,0);
  mbc_line->SetLineColor(kBlue);
  mbc_line->SetLineStyle(1);
  mbc_line->SetLineWidth((Width_t)2.);
  mbc_line->Draw();
  TLine *mbc_lineDOWN = new TLine(mbcMin,-3,mbcMax,-3);
  mbc_lineDOWN->SetLineColor(kBlue);
  mbc_lineDOWN->SetLineStyle(2);
  mbc_lineDOWN->Draw();

  cmmbc->Update();
  out.str("");
  out << "../Reports/pics/mbc_sig_mode" << _mode << "_h0mode" << _h0mode << "_b0f" << _b0f;
  if(projection_flag) out << "_wproj";
  if(deMin>-0.2) out << "015";
  out << ".png";
  if(save_flag) cmmbc->Print(out.str().c_str());

  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.3,0.3),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");
  cmmbc->Update();
  out.str("");
  out << "../Reports/pics/2d_sig_mode" << _mode << "_h0mode" << _h0mode << "_b0f" << _b0f;
  if(projection_flag) out << "_wproj";
  if(deMin>-0.2) out << "015";
  out << ".png";
  if(save_flag) hhc->Print(out.str().c_str());

  TLine* l1 = new TLine(de_sig_min,mbc_min,de_sig_max,mbc_min);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->SetLineWidth(2);
  TLine* l2 = new TLine(de_sig_min,mbc_max,de_sig_max,mbc_max);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(1);
  l2->SetLineWidth(2);
  TLine* l3 = new TLine(de_sig_min,mbc_min,de_sig_min,mbc_max);
  l3->SetLineColor(kRed);
  l3->SetLineStyle(1);
  l3->SetLineWidth(2);
  TLine* l4 = new TLine(de_sig_max,mbc_min,de_sig_max,mbc_max);
  l4->SetLineColor(kRed);
  l4->SetLineStyle(1);
  l4->SetLineWidth(2);

  TCanvas* ellican = new TCanvas("scatterplot","Scatter Plot Mbc dE",400,400);
  ellican->cd();
  out.str("");
  out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX << " && mbc>5.265 && de>" << deMin;
  if(_mode != 1) out << " && mode == " << _mode << " && h0mode == " << _h0mode;
  if(_b0f == 1)      out << " && (b0f == 1 || b0f == 10)";
  else if(_b0f == 5) out << " && b0f != 0 && b0f != -1 && b0f != 1 && b0f != 10";
  else              out << " && b0f != 0 && b0f != -1";
//  cout << "Cuts: " << out.str() << endl;
  tree->Draw("mbc:de",out.str().c_str());
  l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
  out.str("");
  out << "../Reports/pics/scatplot_sig_mode" << _mode << "_h0mode" << _h0mode << "_b0f" << _b0f;
  if(deMin>-0.2) out << "_015" << ".png";
  if(save_flag) ellican->Print(out.str().c_str());

  cout << _mode << " " << _h0mode << " " << _b0f << endl;
}
