#include "cuts.h"
using namespace RooFit;

void Back_rho_2d_fit(const int _mode = 1){
  //TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-2.root");
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1_full.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");
  const int version = 0;

  RooCategory b0f("b0f","b0f");
  RooCategory mode("mode","mode");
  RooCategory h0mode("h0mode","h0mode");
  string label;
  bool gg_flag = true;
  switch (_mode) {
  case 1:
    mode.defineType("pi0",1);
    h0mode.defineType("gg",10);
    b0f.defineType("rho",3);
    label = string("#pi^{0}");
    break;
  case 2:
    mode.defineType("eta",2);
    h0mode.defineType("gg",10);
    b0f.defineType("rho2",2);
    b0f.defineType("rho3",3);
    b0f.defineType("rho4",4);
    b0f.defineType("rho11",11);
    label = string("#eta#rightarrow#gamma#gamma");
    break;
  case 3:
    mode.defineType("eta",2);
    h0mode.defineType("ppp",20);
    b0f.defineType("rho2",2);
    b0f.defineType("rho3",3);
    b0f.defineType("rho4",4);
    b0f.defineType("rho11",11);
    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    gg_flag = false;
    break;
  case 4:
    mode.defineType("omega",3);
    h0mode.defineType("ppp",20);
    b0f.defineType("rho2",2);
    b0f.defineType("rho3",3);
    b0f.defineType("rho4",4);
    b0f.defineType("rho11",11);
    label = string("#omega");
    gg_flag = false;
    break;
  default:
     return;
  }

  RooArgSet argset;
  argset.add(b0f); argset.add(mode); argset.add(h0mode);

  const double mbcMin = 5.20;
  const double mbcMax = 5.29;
  const double deMin = -0.15;
  const double deMax = 0.3;

  RooRealVar mbc("mbc","M_{bc}",mbcMin,mbcMax,"GeV"); argset.add(mbc);
  mbc.setRange("Signal",mbc_min,mbc_max);
  mbc.setRange("mbcSignal",mbc_min,mbc_max);
  mbc.setRange("deSignal",mbcMin,mbcMax);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  de.setRange("mbcSignal",deMin,deMax);
  de.setRange("deSignal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
//  RooRealVar bdtgs("bdtgs","bdtgs",bdtgs_cut,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");
  ds.Print();

  ////////////
  // de pdf //
  ////////////
  const bool deconst = true;
  RooRealVar de0r("de0r","de0r",get_de0r(_mode),-0.2,0.12); if(cRHO || deconst) de0r.setConstant(kTRUE);
  switch (de_rho_param) {
  case 0:
    RooRealVar exppar("exppar","exppar",mr_exppar,-100.,-1.); if(cRHO || deconst) exppar.setConstant(kTRUE);
    RooExponential pdf_de("pdf_de","pdf_de",de,exppar);
    break;
  case -1:
    RooRealVar x0("x0","x0",mr_x0,-0.2,0.12); if(cRHO || deconst) x0.setConstant(kTRUE);
    RooRealVar p1("p1","p1",mr_p1,-1000.,100.); if(cRHO || deconst) p1.setConstant(kTRUE);
    RooRealVar p2("p2","p2",mr_p2,0.,100.); if(cRHO || deconst) p2.setConstant(kTRUE);
    RooGenericPdf pdf_de("pdf_de","1+@0*@1-@2*TMath::Log(1+TMath::Exp(@2*(@0-@1)/@3))",RooArgSet(de,x0,p1,p2));
    break;
  case 1:
    RooRealVar slopel("slopel","slopel",get_slopel(_mode),-1.e5,0.);  if(cRHO || deconst) slopel.setConstant(kTRUE);
    RooRealVar sloper("sloper","sloper",get_sloper(_mode),-10000,0.); if(cRHO || deconst) sloper.setConstant(kTRUE);
    RooRealVar steep("steep","steep",get_steep(_mode),0.,1000.);      if(cRHO || deconst) steep.setConstant(kTRUE);
    RooRealVar p5("p5","p5",get_p5(_mode),0.01,1000.);                if(cRHO || deconst) p5.setConstant(kTRUE);
    RooGenericPdf pdf_de("pdf_de","1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) > 0 ? 1+(@0-@1)*@2+@4*TMath::Log(1+@5*TMath::Exp((@3-@2)*(@0-@1)/@4)) : 0.001",RooArgSet(de,de0r,slopel,sloper,steep,p5));
    break;
  case 2:
    RooRealVar slopel("slopel","slopel",get_slopel(_mode),-10,500.); if(cRHO || deconst) slopel.setConstant(kTRUE);
    RooRealVar steep("steep","steep",mr_steep_2,0.,50000.); if(cRHO || deconst) steep.setConstant(kTRUE);
    RooRealVar exppar("exppar","exppar",mr_exppar_2,1.,500.); if(cRHO || deconst) exppar.setConstant(kTRUE);
    RooRealVar x0("x0","x0",mr_x0_2,-0.2,0.12); if(cRHO || deconst) x0.setConstant(kTRUE);
    RooGenericPdf pdf_de("pdf_de","pdf_de","1+@2*(@0-@1)+@3*TMath::Log(1+TMath::Exp(-@2*(@0-@1)+TMath::Exp(-@4*(@0-@1))))",RooArgSet(de,x0,slopel,steep,exppar));
    pdf_de.Print();
    break;
  default:
     return;
  }

  /////////////
  // mbc pdf //
  /////////////
   switch (mbc_rho_param) {
   case 0:
     RooRealVar mbc0("mbc0","mbc0",mr_mbc0_0,5.26,5.30); if(cRHO) mbc0.setConstant(kTRUE);
     RooRealVar cond("cond","cond",mr_cond_0,-0.1.,1.); if(cRHO) cond.setConstant(kTRUE);
     RooRealVar condr("condr","condr",mr_condr_0,-1000.,1.); if(cRHO) condr.setConstant(kTRUE);
     RooRealVar sl("sl","sl",mr_sl_0,0.,0.5); if(cRHO) sl.setConstant(kTRUE);
     RooRealVar sr("sr","sr",mr_sr_0,0.,0.5); if(cRHO) sr.setConstant(kTRUE);
     RooFormulaVar _sl("_sl","_sl","@0+@1*@2",RooArgList(sl,cond,de));
     RooFormulaVar _sr("_sr","_sr","@0+@1*@2",RooArgList(sr,condr,de));
     RooBifurGauss pdf_mbc("pdf_mbc","pdf_mbc",mbc,mbc0,_sl,_sr);
     break;
   case 1:
     RooRealVar cond("cond","cond",mr_cond_1,-1.,1.); if(cRHO || true) cond.setConstant(kTRUE);
     RooRealVar mbc0("mbc0","mbc0",mr_mbc0_1,5.26,5.30); if(cRHO) mbc0.setConstant(kTRUE);
     RooRealVar sl("sl","sl",mr_sl_1,0.,0.5); if(cRHO) sl.setConstant(kTRUE);
     RooFormulaVar _sl("_sl","_sl","@0+@1*(@2-@3)",RooArgList(sl,cond,de,de0r));
     RooRealVar sr("sr","sr",mr_sr_1,0.,0.5); if(cRHO) sr.setConstant(kTRUE);
     RooBifurGauss bg("bg","bg",mbc,mbc0,_sl,sr);

     RooRealVar mbc00("mbc00","mbc00",mr_mbc00_1,5.26,5.30); if(cRHO) mbc00.setConstant(kTRUE);
     RooRealVar sll("sll","sll",mr_sll_1,0.,0.5); if(cRHO) sll.setConstant(kTRUE);
     RooRealVar srr("srr","srr",mr_srr_1,0.,0.5); if(cRHO) srr.setConstant(kTRUE);
     RooFormulaVar _sll("_sll","_sll","@0+@1*(@2-@3)",RooArgList(sll,cond,de,de0r));
     RooFormulaVar _srr("_srr","_srr","@0+@1*(@2-@3)",RooArgList(srr,cond,de,de0r));
     RooBifurGauss bgg("bgg","bgg",mbc,mbc00,_sll,_srr);

     RooRealVar fmbc("fmbc","fmbc",mr_fmbc_1,0.,1.); if(cRHO) fmbc.setConstant(kTRUE);
     RooAddPdf pdf_mbc("pdf_mbc","pdf_mbc",RooArgList(bg,bgg),RooArgSet(fmbc));
     break;
   case 2:
     RooRealVar mbc1("mbc1","mbc1",mr_mbc1_2,5.26,5.30); if(cRHO)
     mbc1.setConstant(kTRUE);
     RooRealVar ss("ss","ss",mr_ss_2,0.,0.5); if(cRHO) ss.setConstant(kTRUE);
     RooRealVar alpha("alpha","alpha",mr_alpha_2,-10,0.); if(cRHO) alpha.setConstant(kTRUE);
     RooRealVar n("n","n",mr_n_2,0.1,10.); if(cRHO || true) n.setConstant(kTRUE);
     RooCBShape cb("cb","cb",mbc,mbc1,ss,alpha,n);

     RooRealVar mbc0("mbc0","mbc0",mr_mbc0_2,5.26,5.30); if(cRHO) mbc0.setConstant(kTRUE);
     RooRealVar sr("sr","sr",mr_sr_2,0.,0.5); if(cRHO) sr.setConstant(kTRUE);
     RooGaussian gauss("gauss","gauss",mbc,mbc0,sr);
     RooRealVar fcb("fcb","fcb",mr_fcb_2,0.,1.); if(cRHO) fcb.setConstant(kTRUE);
     RooAddPdf pdf_mbc("pdf_mbc","pdf_mbc",RooArgList(cb,gauss),RooArgSet(fcb));
     break;
   case 3:
     RooRealVar argpar("argpar","argus shape parameter",mr_argpar_3,-100.,0.); if(cRHO) argpar.setConstant(kTRUE);
     RooRealVar argedge("argedge","argedge",mr_argedge_3,5.285,5.292);// argedge.setConstant(kTRUE);
     RooArgusBG argus("argus","Argus PDF",mbc,argedge,argpar);

     RooRealVar mbc0("mbc0","mbc0",mr_mbc0_3,5.26,5.30); if(cRHO) mbc0.setConstant(kTRUE);
     RooRealVar cond("cond","cond",mr_cond_3,-0.1.,1.); if(cRHO) cond.setConstant(kTRUE);
     RooRealVar condr("condr","condr",mr_condr_3,-1000.,1.); if(cRHO) condr.setConstant(kTRUE);
     RooRealVar sl("sl","sl",mr_sl_3,0.,0.5); if(cRHO) sl.setConstant(kTRUE);
     RooRealVar sr("sr","sr",mr_sr_3,0.,0.5); if(cRHO) sr.setConstant(kTRUE);
     RooFormulaVar _sl("_sl","_sl","@0+@1*@2",RooArgList(sl,cond,de));
     RooFormulaVar _sr("_sr","_sr","@0+@1*@2",RooArgList(sr,condr,de));
     RooBifurGauss bg("bg","bg",mbc,mbc0,_sl,_sr);
     RooProdPdf pdf_mbc("pdf_mbc","pdf_mbc",RooArgList(bg,argus));
   case 4:
     if(gg_flag){
       RooRealVar b_s("b_s","b_s",get_peak_b_s(_mode),-0.1,0.1);// a_s.setConstant(kTRUE);
       RooRealVar k_s("k_s","k_s",get_peak_k_s(_mode),-0.1,0.1); //b_s.setConstant(kTRUE);
       RooFormulaVar S("S","S","@0+@1*@2",RooArgList(b_s,de,k_s));
       RooRealVar alpha("alpha","alpha",0.139,0.01,2.); alpha.setConstant(kTRUE);
       RooRealVar b_mbc0("b_mbc0","b_mbc0",get_peak_b_mbc0(_mode),5.25,5.29);// a_mbc0.setConstant(kTRUE);
       RooRealVar k_mbc0("k_mbc0","k_mbc0",get_peak_k_mbc0(_mode),-0.1,0.1);// b_mbc0.setConstant(kTRUE);
       RooFormulaVar MBC0("MBC0","MBC0","@0+@1*@2",RooArgList(b_mbc0,de,k_mbc0));
       RooNovosibirsk pdf_mbc("pdf_mbc","pdf_mbc",mbc,MBC0,S,alpha);
     } else{
       RooRealVar argedge("argedge","argedge",5.288,5.285,5.29); //argedge.setConstant(kTRUE);
       RooRealVar argpar_bb("argpar_bb","argpar_bb",get_argpar_bb(_mode),-300,-10.); argpar_bb.setConstant(kTRUE);
       RooArgusBG pdf_mbc_comb_ar("pdf_mbc_comb_ar","Argus PDF",mbc,argedge,argpar_bb);

       RooRealVar mbc0("mbc0","mbc0",get_peak_b_mbc0(_mode),5.25,5.29,"GeV");// mbc0.setConstant(kTRUE);
       RooRealVar mbcWidth("mbcWidth","mbcWidth",get_peak_b_s(_mode),0.,0.1,"GeV");// mbcWidth.setConstant(kTRUE);
       RooGaussian mbcGaus("mbcGaus","mbcGaus",mbc,mbc0,mbcWidth);
       RooRealVar f_g("f_g","f_g",get_f_g_cmb_bb(_mode),0.,1.); //f_g.setConstant(kTRUE);

       RooAddPdf pdf_mbc("pdf_mbc_comb","pdf_mbc_comb",RooArgList(mbcGaus,pdf_mbc_comb_ar),RooArgSet(f_g));
     }
   }

  /////////
  // pdf //
  /////////
//  RooProdPdf pdf("pdf","pdf",RooArgList(pdf_de,pdf_mbc));
  RooProdPdf pdf("pdf","pdf",pdf_de,Conditional(pdf_mbc,mbc));
  pdf.fitTo(ds,Verbose(),Timer(true));

  if(version == 2){
    RooNDKeysPdf pdf("pdf","pdf",RooArgList(de,mbc),ds,RooNDKeysPdf::MirrorBoth,1);//"am");
  }

  /////////////
  //  Plots  //
  /////////////
  // Full region //
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
  pt->AddText(label.c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(de_max,0,de_max,80);
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(de_min,0,de_min,80);
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
  
  // mbc //
  RooPlot* mbcFrame = mbc.frame();
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  pdf.plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  ds.plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc.frame(Title("#Delta E pull distribution"));
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

  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.3,0.75,0.68,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(label.c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(mbc_max,0,mbc_max,50);
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(mbc_min,0,mbc_min,50);
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
  
  TH2D* hh_pdf = pdf.createHistogram("hh_data",de,Binning(50,-0.15,0.1),YVar(mbc,Binning(50,5.26,5.30)));
  hh_pdf->SetLineColor(kBlue);
  TCanvas* hhc = new TCanvas("hhc","hhc",600,600);
  hhc->cd();
  hh_pdf->Draw("SURF");

  cmmbc->Update();
}
