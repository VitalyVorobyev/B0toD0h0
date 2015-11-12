#include "mepdfpartialb.h"

MEPdfPartialB::MEPdfPartialB(RooRealVar *m_de, RooRealVar *m_mbc, const int mode, const int h0mode){
  cout << "MEPdfPartialB constructor..." << endl;
  de = m_de; mbc = m_mbc;
  cuts = new MyParams();
  ggflag = h0mode == 10 ? true : false;
  m_mode = mode; m_h0mode = h0mode;
  InitParams(mode,h0mode);
  FixAll();
  cout << "MEPdfPartialB constructor done." << endl;
}

void MEPdfPartialB::InitParams(const int mode, const int h0mode){
  cout << "InitParams... ";
  de0_part           = new RooRealVar ("de0_part","de0_part",     cuts->get_de0_part(mode,h0mode),-0.2,0.12);    de_param_vec.push_back(de0_part);
  slopel_part        = new RooRealVar("slopel_part","slopel_part",cuts->get_slopel_part(mode,h0mode),-1.e5,0.);  de_param_vec.push_back(slopel_part);
  sloper_part        = new RooRealVar("sloper_part","sloper_part",cuts->get_sloper_part(mode,h0mode),-10000,0.); de_param_vec.push_back(sloper_part);
  steep_part         = new RooRealVar("steep_part","steep_part",  cuts->get_steep_part(mode,h0mode),0.,1000.);   de_param_vec.push_back(steep_part);
  p5_part            = new RooRealVar("p5_part","p5_part",        cuts->get_p5_part(mode,h0mode),0.01,1000.);    de_param_vec.push_back(p5_part);
  pdf_de_part        = new RooRhoDeltaEPdf("pdf_de_part","pdf_de_part",*de,*de0_part,*slopel_part,*sloper_part,*steep_part,*p5_part);

  de0_part->setConstant(kTRUE);
  slopel_part->setConstant(kTRUE);
  sloper_part->setConstant(kTRUE);
  steep_part->setConstant(kTRUE);
  p5_part->setConstant(kTRUE);

  if(ggflag){
    b_s_mbc_part     = new RooRealVar("b_s_mbc_part","b_s_mbc_part",cuts->get_b_s_mbc_part(mode,h0mode),-0.1,0.1); mbc_param_vec.push_back(b_s_mbc_part);
    k_s_mbc_part     = new RooRealVar("k_s_mbc_part","k_s_mbc_part",cuts->get_k_s_mbc_part(mode),-0.1,0.1); k_s_mbc_part->setConstant(kTRUE); mbc_param_vec.push_back(k_s_mbc_part);
    s_mbc_part       = new RooFormulaVar("s_mbc_part","s_mbc_part","@0+@1*@2",RooArgList(*b_s_mbc_part,*de,*k_s_mbc_part));
    alpha_part       = new RooRealVar("alpha_part","alpha_part",0.139,0.01,2.); mbc_param_vec.push_back(alpha_part);
    b_mbc0_part      = new RooRealVar("b_mbc0_part","b_mbc0_part",cuts->get_b_mbc0_part(mode,h0mode),5.25,5.29); mbc_param_vec.push_back(b_mbc0_part);
    k_mbc0_part      = new RooRealVar("k_mbc0_part","k_mbc0_part",cuts->get_k_mbc0_part(mode),-0.1,0.1); k_mbc0_part->setConstant(kTRUE); mbc_param_vec.push_back(k_mbc0_part);
    MBC0_part        = new RooFormulaVar("MBC0_part","MBC0_part","@0+@1*@2",RooArgList(*b_mbc0_part,*de,*k_mbc0_part));

    pdf_mbc_part_gg  = new RooNovosibirsk("pdf_mbc_part_gg","pdf_mbc_part_gg",*mbc,*MBC0_part,*s_mbc_part,*alpha_part);
    pdf_part         = new RooProdPdf("pdf_part","pdf_part",*pdf_de_part,Conditional(*pdf_mbc_part_gg,*mbc));
  } else{
    argedge_part_bb  = new RooRealVar("argedge_part_bb","argedge_part_bb",cuts->get_argedge(mode,h0mode),5.288,5.29); mbc_param_vec.push_back(argedge_part_bb);
    argpar_part_bb   = new RooRealVar("argpar_part_bb","argpar_part_bb",cuts->get_argpar_part_bb(mode,h0mode),-300,-10.); mbc_param_vec.push_back(argpar_part_bb);
    argus_mbc_part   = new RooArgusBG("argus_mbc_part","Argus PDF",*mbc,*argedge_part_bb,*argpar_part_bb);

    mbc0_part        = new RooRealVar("mbc0_part","mbc0_part",cuts->get_mbc0_part(mode,h0mode),5.25,5.291,"GeV"); mbc_param_vec.push_back(mbc0_part);
    s_mbc_part_ppp   = new RooRealVar("s_mbc_part_ppp","s_mbc_part_ppp",cuts->get_s_mbc_part(mode,h0mode),0.,0.1,"GeV"); mbc_param_vec.push_back(s_mbc_part_ppp);
    g_mbc_part       = new RooGaussian("g_mbc_part","g_mbc_part",*mbc,*mbc0_part,*s_mbc_part_ppp);
    fg_mbc_part      = new RooRealVar("fg_mbc_part","fg_mbc_part",cuts->get_fg_mbc_part(mode,h0mode),0.,1.); mbc_param_vec.push_back(fg_mbc_part);

    pdf_mbc_part_ppp = new RooAddPdf("pdf_mbc_part_ppp","pdf_mbc_part_ppp",RooArgList(*g_mbc_part,*argus_mbc_part),RooArgSet(*fg_mbc_part));
    pdf_part         = new RooProdPdf("pdf_part","pdf_part",*pdf_de_part,Conditional(*pdf_mbc_part_ppp,*mbc));
  }
  cout << "done." << endl;
}

void MEPdfPartialB::ChangeParState(const int state_flag){
//  de0_part->setConstant(state_flag);
//  slopel_part->setConstant(state_flag);
//  sloper_part->setConstant(state_flag);
//  steep_part->setConstant(state_flag);
//  p5_part->setConstant(state_flag);

  if(ggflag){
    b_s_mbc_part->setConstant(state_flag);
    k_s_mbc_part->setConstant(state_flag);
    alpha_part->setConstant(state_flag);
    b_mbc0_part->setConstant(state_flag);
    k_mbc0_part->setConstant(state_flag);
  } else{
    argedge_part_bb->setConstant(state_flag);
    argpar_part_bb->setConstant(state_flag);
    mbc0_part->setConstant(state_flag);
    s_mbc_part_ppp->setConstant(state_flag);
    fg_mbc_part->setConstant(state_flag);
  }
}

int MEPdfPartialB::TryParameters(RooDataSet *ds){
  FixAll();
  pdf_part->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  return 0;
}

int MEPdfPartialB::FitParameters(RooDataSet *ds,const bool etaggflag){
  FreeAll();
  if(etaggflag) k_s_mbc_part->setConstant(kTRUE);
  pdf_part->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  return 0;
}

void MEPdfPartialB::DrawDeltaE(RooDataSet* ds){
  mbc->setRange("mbcSignal",cuts->get_mbc_min_h0(m_mode,m_h0mode),cuts->get_mbc_max_h0(m_mode,m_h0mode));
  RooPlot* deFrame = de->frame();
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf_part->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  pdf_part->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

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

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.04);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  TPaveText *pt = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(out.str().c_str());
  pt->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(cuts->get_de_min_h0(m_mode,m_h0mode),0,cuts->get_de_min_h0(m_mode,m_h0mode),80);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(cuts->get_de_max_h0(m_mode,m_h0mode),0,cuts->get_de_max_h0(m_mode,m_h0mode),80);
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
  out << "pics/de_part_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cm->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());
  out.str("");
  out << "pics/de_part_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cm->Print(out.str().c_str());
}

void MEPdfPartialB::DrawMbc(RooDataSet* ds){
  de->setRange("deSignal",cuts->get_de_min_h0(m_mode,m_h0mode),cuts->get_de_max_h0(m_mode,m_h0mode));
  RooPlot* mbcFrame = mbc->frame();
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf_part->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf_part->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc->frame(Title("#Delta E pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("Mbc","Mbc",600,700);
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
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.3,0.75,0.68,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(cuts->get_mbc_max_h0(m_mode,m_h0mode),0,cuts->get_mbc_max_h0(m_mode,m_h0mode),50);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(cuts->get_mbc_min_h0(m_mode,m_h0mode),0,cuts->get_mbc_min_h0(m_mode,m_h0mode),50);
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
  out << "pics/mbc_part_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cmmbc->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());
  out.str("");
  out << "pics/mbc_part_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cmmbc->Print(out.str().c_str());
}

void MEPdfPartialB::PrintParameters(void){
  cout << "Partial B PDF parameters for mode " << m_mode << ", h0mode " << m_h0mode << ":" << endl;
  const int NdePar = de_param_vec.size();
  cout << "Delta E (" << NdePar << " parameters):" << endl;
  for(int i=0; i<NdePar; i++){de_param_vec[i]->Print();}
  const int NmbcPar = mbc_param_vec.size();
  cout << "Mbc (" << NmbcPar << " parameters):" << endl;
  for(int i=0; i<NmbcPar; i++){mbc_param_vec[i]->Print();}
  return;
}

void MEPdfPartialB::WriteParameters(void){
  stringstream out;
  out.str("");
  out << "params/PartParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
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

int MEPdfPartialB::GetParametersFromFile(void){
  stringstream out;
  out.str("");
  out << "params/PartParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
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
    if(name == string("de0_part")){    de0_part->setVal(val); continue;}
    if(name == string("slopel_part")){ slopel_part->setVal(val); continue;}
    if(name == string("sloper_part")){ sloper_part->setVal(val); continue;}
    if(name == string("steep_part")){  steep_part->setVal(val); continue;}
    if(name == string("p5_part")){     p5_part->setVal(val); continue;}
    cout << "MEPdfPartialB::GetParametersFromFile: can't find " << name << endl;
  }
  getline(ifile,line);
  sscanf(line.c_str(),"Mbc (%d parameters):",&npars);
  for(int i=0; i<npars; i++){
    getline(ifile,line);
    cout << line << endl;
    sscanf(line.c_str(),"%s = %lf",namech,&val);
    name = string(namech);
    if(name == string("b_s_mbc_part")){ b_s_mbc_part->setVal(val); continue;}
    if(name == string("k_s_mbc_part")){ k_s_mbc_part->setVal(val); continue;}
    if(name == string("alpha_part")){   alpha_part->setVal(val); continue;}
    if(name == string("b_mbc0_part")){  b_mbc0_part->setVal(val); continue;}
    if(name == string("k_mbc0_part")){  k_mbc0_part->setVal(val); continue;}

    if(name == string("argedge_part_bb")){ argedge_part_bb->setVal(val); continue;}
    if(name == string("argpar_part_bb")){  argpar_part_bb->setVal(val); continue;}
    if(name == string("mbc0_part")){       mbc0_part->setVal(val); continue;}
    if(name == string("s_mbc_part_ppp")){  s_mbc_part_ppp->setVal(val); continue;}
    if(name == string("fg_mbc_part")){     fg_mbc_part->setVal(val); continue;}
    cout << "MEPdfPartialB::GetParametersFromFile: can't find " << name << endl;
  }
  return 0;
}
