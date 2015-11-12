#include "mepdfcombinatorial.h"
#include "TMath.h"

MEPdfCombinatorial::MEPdfCombinatorial(RooRealVar* m_de, RooRealVar* m_mbc, const int mode, const int h0mode){
  cuts = new MyParams();
  de = m_de; mbc = m_mbc;
  m_mode = mode; m_h0mode = h0mode;
  InitParams(mode,h0mode);
  FixAll();
}

void MEPdfCombinatorial::InitParams(const int mode, const int h0mode){
  cout << "MEPdfCombinatorial::InitParams" << endl;
  cout << "de... ";
  c1_qq_cmb      = new RooRealVar("c1_qq_cmb","c1_qq_cmb",cuts->get_c1_qq_cmb(mode,h0mode),-10,50.);       de_param_vec.push_back(c1_qq_cmb);
  c2_qq_cmb      = new RooRealVar("c2_qq_cmb","c2_qq_cmb",cuts->get_c2_qq_cmb(mode,h0mode),-0.1,1);        de_param_vec.push_back(c2_qq_cmb);
  pdf_de_cmb_qq  = new RooChebychev("pdf_de_cmb_qq","pdf_de_cmb_qq",*de,RooArgSet(*c1_qq_cmb,*c2_qq_cmb));
  if(mode != 10){
//  a_c1_bb_cmb    = new RooRealVar("a_c1_bb_cmb","a_c1_bb_cmb",cuts->get_a_c1_bb_cmb(mode,h0mode),-10,50.); a_c1_bb_cmb->setConstant(kTRUE); de_param_vec.push_back(a_c1_bb_cmb);
//  b_c1_bb_cmb    = new RooRealVar("b_c1_bb_cmb","b_c1_bb_cmb",cuts->get_b_c1_bb_cmb(mode,h0mode),-50,0.);  b_c1_bb_cmb->setConstant(kTRUE); de_param_vec.push_back(b_c1_bb_cmb);
//  c1_bb_cmb      = new RooFormulaVar("c1_bb_cmb","@0+@1*@2",RooArgSet(*a_c1_bb_cmb,*b_c1_bb_cmb,*mbc));
//  c2_bb_cmb      = new RooRealVar("c2_bb_cmb","c2_bb_cmb",cuts->get_c2_bb_cmb(mode,h0mode),-10,1);        de_param_vec.push_back(c2_bb_cmb);
    c2_bb_cmb      = new RooRealVar("c2_bb_cmb","c2_bb_cmb",-6.5,-30,-0.1);        de_param_vec.push_back(c2_bb_cmb);
//  pdf_de_cmb_bb  = new RooChebychev("pdf_de_cmb_bb","pdf_de_cmb_bb",*de,RooArgSet(*c1_bb_cmb,*c2_bb_cmb));
    pdf_de_cmb_bb  = new RooExponential("pdf_de_cmb_bb","pdf_de_cmb_bb",*de,*c2_bb_cmb);
  } else{
//    const_bb_cmb   = new RooRealVar("const_bb_cmb","const_bb_cmb",0.5*TMath::Pi());       const_bb_cmb->setConstant(kTRUE); de_param_vec.push_back(const_bb_cmb);
//    de0_bb_cmb     = new RooRealVar("de0_bb_cmb","de0_bb_cmb",3.22957e-02,0.01,0.06);     de_param_vec.push_back(de0_bb_cmb);
//    steep_bb_cmb   = new RooRealVar("steep_bb_cmb","steep_bb_cmb",5.27767e-02,0.01,0.1); de_param_vec.push_back(steep_bb_cmb);
//    pdf_de_dst0_bb_cmb = new RooGenericPdf("pdf_de_bb_cmb","pdf_de_bb_cmb","@0-TMath::ATan((@1-@2)/@3)",RooArgList(*const_bb_cmb,*de,*de0_bb_cmb,*steep_bb_cmb));
//    pdf_de_dst0_bb_cmb = new RooGenericPdf("pdf_de_bb_cmb","pdf_de_bb_cmb","TMath::Erfc((@0-@1)/@2)",RooArgList(*de,*de0_bb_cmb,*steep_bb_cmb));
    c2_bb_cmb      = new RooRealVar("c2_bb_cmb","c2_bb_cmb",-6.5,-30,-0.1);        de_param_vec.push_back(c2_bb_cmb);
    pdf_de_cmb_bb  = new RooExponential("pdf_de_exp_cmb_bb","pdf_de_exp_cmb_bb",*de,*c2_bb_cmb);

    const_bb_cmb   = new RooRealVar("const_bb_cmb","const_bb_cmb",0.5*TMath::Pi()); const_bb_cmb->setConstant(kTRUE); de_param_vec.push_back(const_bb_cmb);
    de0_bb_cmb     = new RooRealVar("de0_bb_cmb","de0_bb_cmb",5.22833e-02,0.05,0.10);     de_param_vec.push_back(de0_bb_cmb);
    steep_bb_cmb   = new RooRealVar("steep_bb_cmb","steep_bb_cmb",4.02311e-02,0.01,0.03); de_param_vec.push_back(steep_bb_cmb);
    fatan          = new RooRealVar("f_atan","f_atan",0.7,0.,1.); de_param_vec.push_back(fatan);
    pdf_de_atan_bb_cmb = new RooGenericPdf("pdf_de_atan_bb_cmb","pdf_de_atan_bb_cmb","@0-TMath::ATan((@1-@2)/@3)",RooArgList(*const_bb_cmb,*de,*de0_bb_cmb,*steep_bb_cmb));
    pdf_de_dst0_bb_cmb = new RooAddPdf("pdf_de_cmb_bb","pdf_de_cmb_bb",RooArgSet(*pdf_de_atan_bb_cmb,*pdf_de_cmb_bb),RooArgList(*fatan));

//    de0_bb_cmb->setConstant(kTRUE); steep_bb_cmb->setConstant(kTRUE);
//    f_atan_cmb     = new RooRealVar("f_atan_cmb","f_atan_cmb",0.1,0.,1.); de_param_vec.push_back(f_atan_cmb);
//    pdf_de_dst0pi0_bb_cmb = new RooAddPdf("pdf_de_dst0pi0_bb_cmb","pdf_de_dst0pi0_bb_cmb",RooArgSet(*pdf_de_dst0_bb_cmb,*pdf_de_cmb_bb),RooArgList(*f_atan_cmb));
  }
  cout << "done." << endl;

  cout << "Mbc... ";
  argedge_cmb      = new RooRealVar("argedge_cmb","argedge_cmb",5.288941178,5.288,5.29); mbc_param_vec.push_back(argedge_cmb);
  argedge_cmb->setConstant(kTRUE);
  argpar_cmb_bb    = new RooRealVar("argpar_cmb_bb","argpar_cmb_bb",cuts->get_argpar_cmb_bb(mode,h0mode),-300,-1.); mbc_param_vec.push_back(argpar_cmb_bb);
//  argus_mbc_cmb_bb = new RooArgusBG("argus_mbc_cmb_bb","Argus PDF",*mbc,*argedge_cmb,*argpar_cmb_bb);
  pdf_mbc_cmb_bb = new RooArgusBG("pdf_mbc_cmb_bb","Argus PDF",*mbc,*argedge_cmb,*argpar_cmb_bb);

//  mbc0_cmb_bb      = new RooRealVar("mbc0_cmb_bb","mbc0_cmb_bb",cuts->get_mbc0_cmb_bb(mode,h0mode),5.25,5.29,"GeV"); mbc_param_vec.push_back(mbc0_cmb_bb);
//  s_mbc_cmb_bb     = new RooRealVar("s_mbc_cmb_bb","s_mbc_cmb_bb",cuts->get_s_mbc_cmb_bb(mode,h0mode),0.,0.1,"GeV"); mbc_param_vec.push_back(s_mbc_cmb_bb);
//  g_mbc_cmb_bb     = new RooGaussian("g_mbc_cmb_bb","g_mbc_cmb_bb",*mbc,*mbc0_cmb_bb,*s_mbc_cmb_bb);

//  fg_cmb_bb        = new RooRealVar("fg_cmb_bb","fg_cmb_bb",cuts->get_fg_cmb_bb(mode,h0mode),0.,1.); mbc_param_vec.push_back(fg_cmb_bb);
//  pdf_mbc_cmb_bb   = new RooAddPdf("pdf_mbc_cmb_bb","pdf_mbc_cmb_bb",RooArgList(*g_mbc_cmb_bb,*argus_mbc_cmb_bb),RooArgSet(*fg_cmb_bb));

  argpar_cmb_qq    = new RooRealVar("argpar_cmb_qq","argpar_cmb_qq",cuts->get_argpar_cmb_qq(mode,h0mode),-300,0.); mbc_param_vec.push_back(argpar_cmb_qq);
//  argpar_cmb_qq    = new RooRealVar("argpar_cmb_qq","argpar_cmb_qq",-10,-300,0.); mbc_param_vec.push_back(argpar_cmb_qq);
  pdf_mbc_cmb_qq   = new RooArgusBG("pdf_mbc_cmb_qq","pdf_mbc_cmb_qq",*mbc,*argedge_cmb,*argpar_cmb_qq);
  cout << "done." << endl;

  fbb_cmb          = new RooRealVar("fbb_cmb","fbb_cmb",0.13,0.,0.9);

//  if(!(mode==33 && h0mode == 20)){
//    fg_cmb_bb->setVal(0);
//    fg_cmb_bb->setConstant(kTRUE);
//    s_mbc_cmb_bb->setConstant(kTRUE);
//    mbc0_cmb_bb->setConstant(kTRUE);
//  }
  cout << "pdf.... ";
  pdf_cmb_qq     = new RooProdPdf("pdf_cmb_qq","pdf_cmb_qq",RooArgSet(*pdf_mbc_cmb_qq,*pdf_de_cmb_qq));
  if(mode != 10){
//    pdf_cmb_bb     = new RooProdPdf("pdf_cmb_bb","pdf_cmb_bb",*pdf_mbc_cmb_bb,Conditional(*pdf_de_cmb_bb,*de));
    pdf_cmb_bb     = new RooProdPdf("pdf_cmb_bb","pdf_cmb_bb",*pdf_mbc_cmb_bb,*pdf_de_cmb_bb);
  } else{
//    pdf_cmb_bb     = new RooProdPdf("pdf_cmb_bb","pdf_cmb_bb",*pdf_mbc_cmb_bb,Conditional(*pdf_de_dst0_bb_cmb,*de));
    pdf_cmb_bb     = new RooProdPdf("pdf_cmb_bb","pdf_cmb_bb",*pdf_mbc_cmb_bb,*pdf_de_dst0_bb_cmb);
  }
  cout << "done." << endl;
  pdf_cmb          = new RooAddPdf("pdf_cmb","pdf_cmb",RooArgSet(*pdf_cmb_bb,*pdf_cmb_qq),RooArgList(*fbb_cmb)); 
}

void MEPdfCombinatorial::ChangeParState(const int state_flag, const int component){
  //  component = 0 -> bb + qq
  //  component = 1 -> qq
  //  component = 2 -> bb
  ChangeMbcParState(state_flag,component);
  ChangeDeltaEParState(state_flag,component);
}

void MEPdfCombinatorial::ChangeMbcParState(const int state_flag, const int component){
  //  component = 0 -> bb + qq
  //  component = 1 -> qq
  //  component = 2 -> bb
  if(component != 1){
    argpar_cmb_bb->setConstant(state_flag);
    //if(m_mode==33 && m_h0mode == 20){
//      mbc0_cmb_bb->setConstant(state_flag);
//      s_mbc_cmb_bb->setConstant(state_flag);
//      fg_cmb_bb->setConstant(state_flag);
//    }
  }
  if(component != 2){
    argpar_cmb_qq->setConstant(state_flag);
  }
//  argedge_cmb->setConstant(state_flag);
}

void MEPdfCombinatorial::ChangeDeltaEParState(const int state_flag, const int component){
  //  component = 0 -> bb + qq
  //  component = 1 -> qq
  //  component = 2 -> bb
  if(component != 1){
//    a_c1_bb_cmb->setConstant(state_flag);
//    b_c1_bb_cmb->setConstant(state_flag);
    if(m_mode != 10){ c2_bb_cmb->setConstant(state_flag);
    } else{
      fatan->setConstant(state_flag);
      c2_bb_cmb->setConstant(state_flag);
      de0_bb_cmb->setConstant(state_flag);
      steep_bb_cmb->setConstant(state_flag);
    }
//        f_atan_cmb->setConstant(state_flag);
  }
  if(component != 2){
    c1_qq_cmb->setConstant(state_flag);
    c2_qq_cmb->setConstant(state_flag);
  }
}

int MEPdfCombinatorial::TryParameters(RooDataSet *ds){
  FixAll();
  pdf_cmb->Print();
//  ds->Print();
  pdf_cmb->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  return 0;
}

int MEPdfCombinatorial::FitParameters(RooDataSet *ds){
  FreeAll();
  FixQQ();
//  fg_cmb_bb->setConstant(kFALSE);
  pdf_cmb->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  FixAll();
  return 0;
}

int MEPdfCombinatorial::FitContParameters(RooDataSet *ds){
  FreeAll();
  FixBB();
//  argpar_cmb_qq->setConstant(kTRUE);
  //c1_qq_cmb->setConstant(kTRUE);
  //c2_qq_cmb->setConstant(kTRUE);
//  fg_cmb_bb->setVal(0); fg_cmb_bb->setConstant(kTRUE);
  pdf_cmb_qq->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds,true);
  PrintParameters();
  FixAll();
  return 0;
}

int MEPdfCombinatorial::TryContParameters(RooDataSet *ds){
//  fg_cmb_bb->setVal(0); fg_cmb_bb->setConstant(kTRUE);
  Draw(ds,true);
  PrintParameters();
  return 0;
}

void MEPdfCombinatorial::PrintParameters(void){
  cout << "Combinatorial PDF parameters for mode " << m_mode << ", h0mode " << m_h0mode << ":" << endl;
  const int NdePar = de_param_vec.size();
  cout << "Delta E (" << NdePar << " parameters):" << endl;
  for(int i=0; i<NdePar; i++){de_param_vec[i]->Print();}
  const int NmbcPar = mbc_param_vec.size();
  cout << "Mbc (" << NmbcPar << " parameters):" << endl;
  for(int i=0; i<NmbcPar; i++){mbc_param_vec[i]->Print();}
  cout << endl;
  fbb_cmb->Print();
  return;
}

void MEPdfCombinatorial::WriteParameters(const bool qqflag){
  stringstream out;
  out.str("");
  if(qqflag) out << "params/ContParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  else       out << "params/CombParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  cout << "Saving parameters in file " << out.str() << endl;
  ofstream ofile;
  ofile.open(out.str().c_str(),ofstream::out);
  const int NdePar = de_param_vec.size();
  ofile << "Delta E (" << NdePar << " parameters):" << endl;
  for(int i=0; i<NdePar; i++){ofile << de_param_vec[i]->format(4,"NE")->Data() << endl;}
  const int NmbcPar = mbc_param_vec.size();
  ofile << "Mbc (" << NmbcPar << " parameters):" << endl;
  for(int i=0; i<NmbcPar; i++){ofile << mbc_param_vec[i]->format(4,"NE")->Data() << endl;}
  return;
}

int MEPdfCombinatorial::GetParametersFromFile(void){
  stringstream out;
  out.str("");
  cout << "GPFF Cmb: mode: " << m_mode << ", h0mode: " << m_h0mode << endl;
  out << "params/CombParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  cout << "Getting parameters from file " << out.str() << endl;
  ifstream ifile;
  ifile.open(out.str().c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "Can't open file " << out.str() << endl;
    return -1;
  }
  string line,name;
  int npars;
  double val;
  char namech[15];
  getline(ifile,line);
  sscanf(line.c_str(),"Delta E (%d parameters):",&npars);
  for(int i=0; i<npars; i++){
    getline(ifile,line);
    cout << line << endl;
    sscanf(line.c_str(),"%s = %lf",namech,&val);
    name = string(namech);
//    if(name == string("a_c1_bb_cmb")){  a_c1_bb_cmb->setVal(val); continue;}
//    if(name == string("b_c1_bb_cmb")){  b_c1_bb_cmb->setVal(val); continue;}
    if(name == string("c2_bb_cmb")){    c2_bb_cmb->setVal(val); continue;}
    if(name == string("fatan")){        fatan->setVal(val); continue;}
//    if(name == string("f_atan_cmb")){   f_atan_cmb->setVal(val); continue;}

    if(name == string("c1_qq_cmb")){    c1_qq_cmb->setVal(val); continue;}
    if(name == string("c2_qq_cmb")){    c2_qq_cmb->setVal(val); continue;}

    if(name == string("const_bb_cmb")){ const_bb_cmb->setVal(val); continue;}
    if(name == string("de0_bb_cmb")){   de0_bb_cmb->setVal(val); continue;}
    if(name == string("steep_bb_cmb")){ steep_bb_cmb->setVal(val); continue;}

//    if(name == string("const_qq_cmb")){ const_qq_cmb->setVal(val); continue;}
//    if(name == string("de0_qq_cmb")){   de0_qq_cmb->setVal(val); continue;}
//    if(name == string("steep_qq_cmb")){ steep_qq_cmb->setVal(val); continue;}
    cout << "MEPdfCombinatorial::GetParametersFromFile: cannot find " << name << endl;
  }
  getline(ifile,line);
  sscanf(line.c_str(),"Mbc (%d parameters):",&npars);
  for(int i=0; i<npars; i++){
    getline(ifile,line);
    cout << line << endl;
    sscanf(line.c_str(),"%s = %lf",namech,&val);
    name = string(namech);
    if(name == string("argedge_cmb")){   argedge_cmb->setVal(val); continue;}
    if(name == string("argpar_cmb_bb")){ argpar_cmb_bb->setVal(val); continue;}
//    if(name == string("mbc0_cmb_bb")){   mbc0_cmb_bb->setVal(val); continue;}
//    if(name == string("s_mbc_cmb_bb")){  s_mbc_cmb_bb->setVal(val); continue;}
    if(name == string("argpar_cmb_qq")){ argpar_cmb_qq->setVal(val); continue;}
//    if(name == string("fg_cmb_bb")){     fg_cmb_bb->setVal(val); continue;}
    cout << "MEPdfCombinatorial::GetParametersFromFile: cannot find " << name << endl;
  }
  return 0;
}

void MEPdfCombinatorial::DrawDeltaE(RooDataSet* ds, const bool qqflag){
//  mbc->setRange("mbcSignal",cuts->get_mbc_min_h0(m_mode,m_h0mode),cuts->get_mbc_max_h0(m_mode,m_h0mode));
  RooPlot* deFrame = de->frame();
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  if(qqflag) pdf_cmb_qq->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  else{
    pdf_cmb->plotOn(deFrame,Components(*pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("mbcSignal"));
    pdf_cmb->plotOn(deFrame,Components(*pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
    pdf_cmb->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  }
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  if(qqflag) pdf_cmb_qq->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  else       pdf_cmb->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de->frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E, Combinatorics","#Delta E, Combinatorics",600,700);
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
  TPaveText *pt = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
  pt->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(cuts->get_de_max_h0(m_mode,m_h0mode),0,cuts->get_de_max_h0(m_mode,m_h0mode),de_line_size(qqflag));
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(cuts->get_de_min_h0(m_mode,m_h0mode),0,cuts->get_de_min_h0(m_mode,m_h0mode),de_line_size(qqflag));
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(cuts->get_de_fit_min(),3,cuts->get_de_fit_max(),3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(cuts->get_de_fit_min(),0,cuts->get_de_fit_max(),0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(cuts->get_de_fit_min(),-3,cuts->get_de_fit_max(),-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
  out.str("");
  if(qqflag) out << "pics/de_cont_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  else       out << "pics/de_comb_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cm->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());
  out.str("");
  if(qqflag) out << "pics/de_cont_m" << m_mode << "_h0m" << m_h0mode << ".root";
  else       out << "pics/de_comb_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cm->Print(out.str().c_str());
}

void MEPdfCombinatorial::DrawMbc(RooDataSet* ds,const bool qqflag){
//  de->setRange("deSignal",cuts->get_de_min_h0(m_mode,m_h0mode),cuts->get_de_max_h0(m_mode,m_h0mode));
  RooPlot* mbcFrame = mbc->frame();
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  if(qqflag) pdf_cmb_qq->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  else{
    pdf_cmb->plotOn(mbcFrame,Components(*pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("deSignal"));
    pdf_cmb->plotOn(mbcFrame,Components(*pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("deSignal"));
    pdf_cmb->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  }
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  if(qqflag) pdf_cmb_qq->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  else       pdf_cmb->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc->frame(Title("M_{bc} pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("M_{bc}, Combinatorics","M_{bc}, Combinatorics",600,700);
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

  stringstream out;
  TPaveText *ptmbc = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(cuts->get_mbc_max_h0(m_mode,m_h0mode),0,cuts->get_mbc_max_h0(m_mode,m_h0mode),mbc_line_size(qqflag));
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(cuts->get_mbc_min_h0(m_mode,m_h0mode),0,cuts->get_mbc_min_h0(m_mode,m_h0mode),mbc_line_size(qqflag));
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();

  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  mbcPull->SetMarkerSize(0.05);
  mbcPull->Draw();
  TLine *mbc_lineUP = new TLine(cuts->get_mbc_fit_min(),3,cuts->get_mbc_fit_max(),3);
  mbc_lineUP->SetLineColor(kBlue);
  mbc_lineUP->SetLineStyle(2);
  mbc_lineUP->Draw();
  TLine *mbc_line = new TLine(cuts->get_mbc_fit_min(),0,cuts->get_mbc_fit_max(),0);
  mbc_line->SetLineColor(kBlue);
  mbc_line->SetLineStyle(1);
  mbc_line->SetLineWidth((Width_t)2.);
  mbc_line->Draw();
  TLine *mbc_lineDOWN = new TLine(cuts->get_mbc_fit_min(),-3,cuts->get_mbc_fit_max(),-3);
  mbc_lineDOWN->SetLineColor(kBlue);
  mbc_lineDOWN->SetLineStyle(2);
  mbc_lineDOWN->Draw();

  cmmbc->Update();
  out.str("");
  if(qqflag) out << "pics/mbc_cont_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  else       out << "pics/mbc_comb_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cmmbc->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());
  out.str("");
  if(qqflag) out << "pics/mbc_cont_m" << m_mode << "_h0m" << m_h0mode << ".root";
  else       out << "pics/mbc_comb_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cmmbc->Print(out.str().c_str());
}

double MEPdfCombinatorial::de_line_size(const bool qqflag){
  switch(m_mode){
  case 1:  return qqflag ? 120 : 140;
  case 2:
    if(m_h0mode == 10) return qqflag ? 70 : 100;
    else               return qqflag ? 40 : 50;
  case 3:  return qqflag ? 140 : 200;
  case 5:  return qqflag ? 16  : 20;
  case 10: return qqflag ? 30  : 60;
  case 20: return qqflag ? 8   : 14;
  default:
    return -1;
  }
}

double MEPdfCombinatorial::mbc_line_size(const bool qqflag){
  switch(m_mode){
  case 1:  return qqflag ? 100 : 120;
  case 2:
    if(m_h0mode == 10) return qqflag ? 60 : 80;
    else               return qqflag ? 40 : 40;
  case 3:  return qqflag ? 120 : 160;
  case 5:  return qqflag ? 16  : 20;
  case 10: return qqflag ? 16  : 40;
  case 20: return qqflag ? 8   : 10;
  default:
    return -1;
  }
}
