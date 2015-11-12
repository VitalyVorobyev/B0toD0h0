#include "cuts.h"
using namespace RooFit;

void Back_2d_fit(const int _mode = 1, const int type = 0, const bool draw_bins = false){
  // type 0 -> qq
  // type 1 -> bb
  // type 2 -> qq + bb
  const bool save_flag = true;
  TChain* tree = new TChain("TEvent");
  switch (type) {
  case 0:
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_cont_0-1.root");
    break;
  case 1:
//    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1.root");
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_bb_0-1.root");
    break;
  case 2:
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1.root");
    break;
  default:
    return;
    break;
  }

  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);

  RooCategory mode("mode","mode");
  RooCategory h0mode("h0mode","h0mode");

  string MODE,Mode;
  double mbcMin = 5.2;
  const double mbcMax = 5.288;
  const double deMin =-0.3;
  const double deMax = 0.3;
  double de_sig_min,de_sig_max;
  double BDTG_MIN = 0;
  double BDTG_MAX = 1;
  switch(_mode){
  case 1: BDTG_MIN = bdtg_cut_pi0;
          mode.defineType("pi0",1);
          h0mode.defineType("gg",10);
          MODE = string("pi0");
          Mode = string("#pi0");
          de_sig_min = de_min;
          de_sig_max = de_max;
          break;
  case 2: BDTG_MIN = bdtg_cut_etagg;
          mode.defineType("eta",2);
          h0mode.defineType("gg",10);
          MODE = string("etagg");
          Mode = string("#eta#rightarrow#gamma#gamma");
          de_sig_min = de_min;
          de_sig_max = de_max;
          break;
  case 3: BDTG_MIN = bdtg_cut_etappp;
          mode.defineType("eta",2);
          h0mode.defineType("ppp",20);
          MODE = string("etappp");
          de_sig_min = de_min_etappp;
          de_sig_max = de_max_etappp;
          Mode = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
          break;
  case 4: BDTG_MIN = bdtg_cut_omega;
          mode.defineType("omega",3);
          h0mode.defineType("ppp",20);
          MODE = string("omega");
          Mode = string("#omega");
          de_sig_min = de_min_omega;
          de_sig_max = de_max_omega;
          break;
  default:
          return;
  }

  RooArgSet argset;
  argset.add(mode);
  argset.add(h0mode);
  argset.add(b0f);

//  RooCategory flv("flv_mc","flv_mc");
//  flv.defineType("B0",1);
//  flv.defineType("anti-B0",-1);
//  argset.add(flv);

  RooCategory bin("bin","bin");
  bin.defineType("1",1); bin.defineType("-1",-1);
  bin.defineType("2",2); bin.defineType("-2",-2);
  bin.defineType("3",3); bin.defineType("-3",-3);
  bin.defineType("4",4); bin.defineType("-4",-4);
  bin.defineType("5",5); bin.defineType("-5",-5);
  bin.defineType("6",6); bin.defineType("-6",-6);
  bin.defineType("7",7); bin.defineType("-7",-7);
  bin.defineType("8",8); bin.defineType("-8",-8);
  argset.add(bin);
//  RooCategory good_icpv("good_icpv","good_icpv");
  
  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
//  de.setRange("Signal",de_min,de_max);
  de.setRange("Signal",de_sig_min,de_sig_max);
  de.setRange("mbcSignal",deMin,deMax);
//  de.setRange("deSignal",de_min,de_max);
  de.setRange("deSignal",de_sig_min,de_sig_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtg("bdtg","bdtg",BDTG_MIN,BDTG_MAX); argset.add(bdtg);
//  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut);// argset.add(atckpi_max);

  RooDataSet ds("ds","ds",tree,argset,"(mbc>0||mbc<=0) && (de>0||de<=0)");// to eliminate NANs

  ds.Print();

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  if(type != 0){
    RooRealVar c10("c10","c10",get_cmb_c10(_mode),-10,50.); if(type == 2) c10.setConstant(kTRUE);
    RooRealVar c11("c11","c11",get_cmb_c11(_mode),-50,0.);  if(type == 2) c11.setConstant(kTRUE);
    RooFormulaVar c1("c1","@0+@1*@2",RooArgSet(c10,c11,mbc));
    RooRealVar c2("c2","c2",get_cmb_c20(_mode),-0.1,1);     if(type == 2) c2.setConstant(kTRUE);
    if(type == 1) RooChebychev pdf_de_comb("pdf_de_comb","pdf_de_comb",de,RooArgSet(c1,c2));
    else          RooChebychev pdf_de_comb_bb("pdf_de_comb_bb","pdf_de_comb_bb",de,RooArgSet(c1,c2));
  }
  if(type != 1){
    RooRealVar C1("C1","C1",get_cmb_c1(_mode),-10,50.); if(type == 2) C1.setConstant(kTRUE);
    RooRealVar C2("C2","C2",get_cmb_c2(_mode),-0.1,1);  if(type == 2) C2.setConstant(kTRUE);
    if(type == 0) RooChebychev pdf_de_comb("pdf_de_comb","pdf_de_comb",de,RooArgSet(C1,C2));
    else          RooChebychev pdf_de_comb_qq("pdf_de_comb_qq","pdf_de_comb_qq",de,RooArgSet(C1,C2));
  }

  /////////////
  // mbc pdf //
  /////////////
  RooRealVar argedge("argedge","argedge",5.288,5.285,5.29); //argedge.setConstant(kTRUE);
  if(type != 0){
    RooRealVar argpar_bb("argpar_bb","argpar_bb",get_argpar_bb(_mode),-300,-10.); if(type == 2) argpar_bb.setConstant(kTRUE);
    RooArgusBG pdf_mbc_comb_ar("pdf_mbc_comb_ar","Argus PDF",mbc,argedge,argpar_bb);

    RooRealVar mbc0("mbc0","mbc0",get_mbc0_cmb_bb(_mode),5.25,5.29,"GeV");           if(type == 2) mbc0.setConstant(kTRUE);
    RooRealVar mbcWidth("mbcWidth","mbcWidth",get_mbcw_cmb_bb(_mode),0.,0.1,"GeV"); if(type == 2) mbcWidth.setConstant(kTRUE);
    RooGaussian mbcGaus("mbcGaus","mbcGaus",mbc,mbc0,mbcWidth);
    RooRealVar f_g("f_g","f_g",get_f_g_cmb_bb(_mode),0.,1.); if(type == 2) f_g.setConstant(kTRUE);

    if(type == 1) RooAddPdf pdf_mbc_comb("pdf_mbc_comb","pdf_mbc_comb",RooArgList(mbcGaus,pdf_mbc_comb_ar),RooArgSet(f_g));
    else          RooAddPdf pdf_mbc_comb_bb("pdf_mbc_comb_bb","pdf_mbc_comb_bb",RooArgList(mbcGaus,pdf_mbc_comb_ar),RooArgSet(f_g));
  }
  if(type != 1){
    RooRealVar argpar_qq("argpar_qq","argpar_qq",get_argpar_qq(_mode),-300,-10.); if(type == 2) argpar_qq.setConstant(kTRUE);
    if(type == 0) RooArgusBG pdf_mbc_comb("pdf_mbc_comb","Argus PDF",mbc,argedge,argpar_qq);
    else          RooArgusBG pdf_mbc_comb_qq("pdf_mbc_comb_qq","pdf_mbc_comb_qq",mbc,argedge,argpar_qq);
  }

  /////////
  // pdf //
  /////////
  if(type != 2){ RooProdPdf pdf("pdf","pdf",pdf_mbc_comb,Conditional(pdf_de_comb,de));}
  else{
    RooRealVar f_bb("f_bb","f_bb",0.2,0.,1.);
    RooProdPdf pdf_bb("pdf_bb","pdf_bb",pdf_mbc_comb_bb,Conditional(pdf_de_comb_bb,de));
    RooProdPdf pdf_qq("pdf_qq","pdf_qq",pdf_mbc_comb_qq,Conditional(pdf_de_comb_qq,de));
    RooAddPdf pdf("pdf","pdf",RooArgSet(pdf_bb,pdf_qq),RooArgList(f_bb));
  }

  RooFitResult* r = pdf.fitTo(ds,Verbose(),Timer(true));
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf.plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  pdf.plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
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
  out << "mode: " << Mode;
  pt->AddText(out.str().c_str());
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,120);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,120);
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
  out.str("");
  out << "../Reports/pics/de_2dfit_comb_" << MODE << ".png";
  cm->Print(out.str().c_str());
  out.str("");
  out << "../Reports/pics/de_2dfit_comb_" << MODE << ".root";
  cm->Print(out.str().c_str());
  
  // mbc //
  RooPlot* mbcFrame = mbc.frame();
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf.plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc.frame(Title("M_{bc} pull distribution"));
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

  TPaveText *ptmbc = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  out.str("");
  out << "mode: " << Mode;
  ptmbc->AddText(out.str().c_str());
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out.str().c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,70);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,70);
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();
  
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
  out << "../Reports/pics/mbc_2dfit_comb_" << MODE << ".png";
  cmmbc->Print(out.str().c_str());
  out.str("");
  out << "../Reports/pics/mbc_2dfit_comb_" << MODE << ".root";
  cmmbc->Print(out.str().c_str());
 
  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.15,0.1),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");
//  Chi2->Print();

  if(!draw_bins) return;
  for(int i=0; i<2; i++){
    for(int j=0; j<16; j++){
      int Flv = flv(i);
      if(Flv == 0) Flv = -1;
      const int Bin = bin(j);
      out.str("");
      out << "flv_mc == " << Flv << " && bin == " << Bin;
      RooDataSet* ds0 = (RooDataSet*)ds.reduce(RooArgSet(de,mbc),out.str().c_str());

      // de //
      RooPlot* deFrame0 = de.frame();
      ds0->plotOn(deFrame0,DataError(RooAbsData::SumW2),MarkerSize(1));//,CutRange("mbcSignal"));
      pdf.plotOn(deFrame0,LineWidth(2));//,ProjectionRange("mbcSignal"));
      ds0->statOn(deFrame0,Layout(0.55,0.98,0.9));

      RooHist* hdepull0 = deFrame0->pullHist();
      out.str("");
      out << "#Delta E pull distribution, flv: " << Flv << ", bin: " << Bin;
      RooPlot* dePull0 = de.frame(Title(out.str().c_str()));
      dePull0->addPlotable(hdepull0,"P");
      dePull0->GetYaxis()->SetRangeUser(-5,5);

      out.str("");
      out << "#Delta E, Combinatorics, flv: " << Flv << ", bin: " << Bin;
      TCanvas* cm0 = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
      cm0->cd();

      TPad *pad30 = new TPad("pad30","pad30",0.01,0.20,0.99,0.99);
      TPad *pad40 = new TPad("pad40","pad40",0.01,0.01,0.99,0.20);
      pad30->Draw();
      pad40->Draw();

      pad30->cd();
      pad30->SetLeftMargin(0.15);
      pad30->SetFillColor(0);

      deFrame0->GetXaxis()->SetTitleSize(0.05);
      deFrame0->GetXaxis()->SetTitleOffset(0.85);
      deFrame0->GetXaxis()->SetLabelSize(0.04);
      deFrame0->GetYaxis()->SetTitleOffset(1.6);
      deFrame0->Draw();

      stringstream out;
      out.str("");
      out << "#chi^{2}/n.d.f = " << deFrame0->chiSquare();
      TPaveText *pt0 = new TPaveText(0.6,0.55,0.98,0.7,"brNDC");
      pt0->SetFillColor(0);
      pt0->SetTextAlign(12);
      pt0->AddText(out.str().c_str());
      pt0->AddText("qq + BB comb");
      pt0->Draw();

      TLine *de_line_RIGHT0 = new TLine(de_sig_max,0,de_sig_max,250);
      de_line_RIGHT0->SetLineColor(kRed);
      de_line_RIGHT0->SetLineStyle(1);
      de_line_RIGHT0->SetLineWidth((Width_t)2.);
      de_line_RIGHT0->Draw();
      TLine *de_line_LEFT0 = new TLine(de_sig_min,0,de_sig_min,250);
      de_line_LEFT0->SetLineColor(kRed);
      de_line_LEFT0->SetLineStyle(1);
      de_line_LEFT0->SetLineWidth((Width_t)2.);
      de_line_LEFT0->Draw();

      pad40->cd(); pad40->SetLeftMargin(0.15); pad40->SetFillColor(0);
      dePull0->SetMarkerSize(0.05); dePull0->Draw();
      TLine *de_lineUP0 = new TLine(deMin,3,deMax,3);
      de_lineUP0->SetLineColor(kBlue);
      de_lineUP0->SetLineStyle(2);
      de_lineUP0->Draw();
      TLine *de_line0 = new TLine(deMin,0,deMax,0);
      de_line0->SetLineColor(kBlue);
      de_line0->SetLineStyle(1);
      de_line0->SetLineWidth((Width_t)2.);
      de_line0->Draw();
      TLine *de_lineDOWN0 = new TLine(deMin,-3,deMax,-3);
      de_lineDOWN0->SetLineColor(kBlue);
      de_lineDOWN0->SetLineStyle(2);
      de_lineDOWN0->Draw();

      cm0->Update();
      out.str("");
      out << "../Note/pics/de_comb_m" << _mode << "_flv" << Flv << "_bin" << Bin;
      out << ".eps";
      if(save_flag) cm0->Print(out.str().c_str());

      // mbc //
      RooPlot* mbcFrame0 = mbc.frame();
      ds0->plotOn(mbcFrame0,DataError(RooAbsData::SumW2),MarkerSize(1));//,CutRange("deSignal"));
      pdf.plotOn(mbcFrame0,LineWidth(2));//,ProjectionRange("deSignal"));

      ds0->statOn(mbcFrame0,Layout(0.2,0.68,0.9));

      out.str("");
      out << "M_{bc} pull distribution, flv: " << Flv << ", bin: " << Bin;
      RooHist* hmbcpull0 = mbcFrame0->pullHist();
      RooPlot* mbcPull0 = mbc.frame(Title(out.str().c_str()));
      mbcPull0->addPlotable(hmbcpull0,"P");
      mbcPull0->GetYaxis()->SetRangeUser(-5,5);

      out.str("");
      out << "M_{bc}, Combinatorics, flv: " << Flv << ", bin: " << Bin;
      TCanvas* cmmbc0 = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
      cmmbc0->cd();

      TPad *pad10 = new TPad("pad10","pad10",0.01,0.20,0.99,0.99);
      TPad *pad20 = new TPad("pad20","pad20",0.01,0.01,0.99,0.20);
      pad10->Draw();
      pad20->Draw();

      pad10->cd();
      pad10->SetLeftMargin(0.15);
      pad10->SetFillColor(0);

      mbcFrame0->GetXaxis()->SetTitleSize(0.05);
      mbcFrame0->GetXaxis()->SetTitleOffset(0.85);
      mbcFrame0->GetXaxis()->SetLabelSize(0.04);
      mbcFrame0->GetYaxis()->SetTitleOffset(1.6);
      mbcFrame0->Draw();

      out.str("");
      out << "#chi^{2}/n.d.f = " << mbcFrame0->chiSquare();
      TPaveText *ptmbc0 = new TPaveText(0.3,0.55,0.68,0.7,"brNDC");
      ptmbc0->SetFillColor(0);
      ptmbc0->SetTextAlign(12);
      ptmbc0->AddText(out.str().c_str());
      ptmbc0->AddText("qq + BB comb");
      ptmbc0->Draw();

      TLine *mbc_line_RIGHT0 = new TLine(mbc_max,0,mbc_max,800);
      mbc_line_RIGHT0->SetLineColor(kRed);
      mbc_line_RIGHT0->SetLineStyle(1);
      mbc_line_RIGHT0->SetLineWidth((Width_t)2.);
      mbc_line_RIGHT0->Draw();
      TLine *mbc_line_LEFT0 = new TLine(mbc_min,0,mbc_min,800);
      mbc_line_LEFT0->SetLineColor(kRed);
      mbc_line_LEFT0->SetLineStyle(1);
      mbc_line_LEFT0->SetLineWidth((Width_t)2.);
      mbc_line_LEFT0->Draw();

      pad20->cd();
      pad20->SetLeftMargin(0.15);
      pad20->SetFillColor(0);
      mbcPull0->SetMarkerSize(0.05);
      mbcPull0->Draw();
      TLine *mbc_lineUP0 = new TLine(mbcMin,3,mbcMax,3);
      mbc_lineUP0->SetLineColor(kBlue);
      mbc_lineUP0->SetLineStyle(2);
      mbc_lineUP0->Draw();
      TLine *mbc_line0 = new TLine(mbcMin,0,mbcMax,0);
      mbc_line0->SetLineColor(kBlue);
      mbc_line0->SetLineStyle(1);
      mbc_line0->SetLineWidth((Width_t)2.);
      mbc_line0->Draw();
      TLine *mbc_lineDOWN0 = new TLine(mbcMin,-3,mbcMax,-3);
      mbc_lineDOWN0->SetLineColor(kBlue);
      mbc_lineDOWN0->SetLineStyle(2);
      mbc_lineDOWN0->Draw();

      cmmbc0->Update();
      out.str("");
      out << "../Note/pics/mbc_comb_m" << _mode << "_flv" << Flv << "_bin" << Bin;
      out << ".eps";
      if(save_flag) cmmbc->Print(out.str().c_str());
    }
  }
}
