#include "cuts.h"
using namespace RooFit;

void Data_de_rho_fit(void){
  const string data_file("/home/vitaly/B0toDh0/TMVA/fil_b2dh_data.root");
  const string gen_file("/home/vitaly/B0toDh0/TMVA/fil_b2dh_gen.root");
  const string sig_file("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig_full.root");
  
  TFile *ifile = TFile::Open(data_file.c_str());
//  TFile *ifile = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/rooksfw/fil_b2dh_data.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  RooArgSet argset;
  RooRealVar mbc("mbc","mbc",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",-0.3,0.3,"GeV"); argset.add(de);
  de.setRange("signal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
//  RooRealVar bdt("bdt","bdt",bdt_cut,1.); argset.add(bdt);
  RooRealVar bdtgs("bdtgs","bdtgs",bdtgs_cut,1.); argset.add(bdtgs);
//  RooRealVar lh("lh","lh",0.85,1.); argset.add(lh);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  RooCategory exp("exp","exp");
  exp.defineType("7",7);
  exp.defineType("9",9);
  exp.defineType("11",11);
  exp.defineType("13",13);
  exp.defineType("15",15);
  exp.defineType("17",17);
  exp.defineType("19",19);
  exp.defineType("21",21);
  exp.defineType("23",23);
  exp.defineType("25",25);
  exp.defineType("27",27);
  exp.defineType("31",31);
  exp.defineType("33",33);
  exp.defineType("35",35);
  exp.defineType("37",37);
  exp.defineType("39",39);
  exp.defineType("41",41);
  if(!flag41){
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

  RooDataSet ds("ds","ds",tree,argset);
  ds.Print();

  RooRealVar shift("shift","shift",0,-0.1,0.1); shift.setConstant(kTRUE);
//  RooRealVar scale("scale","scale",1.,0.1,10.);
//  RooFormulaVar _de("_de","_de","@0*@1",RooArgList(de,scale));

  RooRealVar de0("de0","de0",1.47504e-02,-0.1,0.1);
  RooFormulaVar _de0("_de0","_de0","@0+@1",RooArgList(de0,shift));
  RooRealVar s1("s1","s1",2.31062e-02,0.,0.5);
  RooGaussian g1("g1","g1",de,_de0,s1);

  RooRealVar deCBl("deCBl","deCBl",-2.33551e-02,-0.1,0.1);
  RooFormulaVar _deCBl("_deCBl","_deCBl","@0+@1",RooArgList(deCBl,shift));
  RooRealVar sCBl("sCBl","sCBl",3.25778e-02,0.,0.5);
  RooRealVar nl("nl","nl",2.60160e+00,0.,100.);
  RooRealVar alphal("alphal","alphal",7.14010e-01,-10.,10.);

  RooRealVar deCBr("deCBr","deCBr",-2.90478e-02,-0.1,0.1);
  RooFormulaVar _deCBr("_deCBr","_deCBr","@0+@1",RooArgList(deCBr,shift));
  RooRealVar sCBr("sCBr","sCBr",5.08916e-02,0.,0.5);
  RooRealVar nr("nr","nr",9.34545e+00,0.,100.); 
  RooRealVar alphar("alphar","alphar",-9.83305e-01,-10.,10.);

  RooCBShape CBl("CBl","CBl",de,_deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,_deCBr,sCBr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",5.89364e-01,0.,1.); 
  RooRealVar fCBr("fCBr","fCBr",1.14349e-01,0.,1.);

  RooAddPdf pdf_sig("pdf_sig","pdf_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));

  RooRealVar c1("c1","c1",-4.27747e-01,-1.6,0.00001);
  RooRealVar c2("c2","c2",3.69914e-01,0.01,0.5);
  RooGenericPdf cheb("cheb","1+@1*@0+@2*(2*@0*@0-1)",RooArgSet(de,c1,c2));
//  RooChebychev cheb("cheb","cheb",de,RooArgSet(c1,c2));

  RooRealVar NSig("NSig","NSig",330,1.,1000.);
  RooRealVar NBack("NBack","NBack",1000,1.,4000.);

  stringstream out;
  if(fullfit){
    cout << "Signal fit:" << endl;
  }
  TFile *ifile_sig = TFile::Open(sig_file.c_str());
  TTree *tree_sig = (TTree*)ifile_sig->Get("TEvent");
  RooArgSet argset_sig = argset;
//  RooCategory b0f("b0f","b0f");
//  b0f.defineType("signal",1);
//  b0f.defineType("fsr",10);
//  b0f.defineType("bad_pi0",5);
////  argset_sig.add(argset);
//  argset_sig.add(b0f);
  RooDataSet ds_sig("ds_sig","ds_sig",tree_sig,argset_sig);
  if(fullfit){
    pdf_sig.fitTo(ds_sig,Timer(true));
  }

  if(!const_back_flag){
    cout << "Background fit:" << endl;
  }
  TFile *ifile_back = TFile::Open(gen_file.c_str());
  TTree *tree_back = (TTree*)ifile_back->Get("TEvent");
  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
  b0f.defineType("rho",3);
  b0f.defineType("unknown",0);
  RooArgSet argset_back(argset);
  argset_back.add(b0f);
  RooDataSet ds_back("ds_back","ds_back",tree_back,argset_back,"b0f != 3");
  RooDataSet ds_rho("ds_rho","ds_rho",tree_back,argset_back,"b0f == 3");
//  RooFormulaVar desh("desh","desh","@0+@1",RooArgList(de,shift));

  ds_back.Print();
  ds_rho.Print();

  RooKeysPdf kest("kest","kest",de,ds_rho,RooKeysPdf::MirrorBoth,1);
  de.setBins(1000);
  RooDataHist* hist = new RooDataHist("hist","hist",RooArgSet(de));
  hist->Print();
  RooArgSet* deset = new RooArgSet(de);
  deset->Print();
  kest.fillDataHist(hist,deset,1);
  hist->Print();
  RooHistPdf histpdf("histpdf","histpdf",RooArgSet(de),*hist);//????
  histpdf.Print();
  RooRealVar f("f","f",0.1,0.,1.);
  RooAddPdf pdf_back("pdf_back","pdf_back",RooArgList(cheb,kest),RooArgSet(f));
  de.setBins(60);
//  if(!const_back_flag){
    cheb.fitTo(ds_back,Range(-0.15,de_fit_max),Timer(true));
//  }
  out.str("");
  out << "de<" << de_max << "&&de>" << de_min;
  const int gen_back = (ds_back.sumEntries(out.str().c_str()) + ds_rho.sumEntries(out.str().c_str()))/6.;
  cout << "Gen background: " << gen_back << endl;

  de0.setConstant(kTRUE);
  s1.setConstant(kTRUE);
  deCBl.setConstant(kTRUE);
  sCBl.setConstant(kTRUE);
  nl.setConstant(kTRUE);
  alphal.setConstant(kTRUE);
  deCBr.setConstant(kTRUE);
  sCBr.setConstant(kTRUE);
  nr.setConstant(kTRUE);
  alphar.setConstant(kTRUE);
  fCBl.setConstant(kTRUE);
  fCBr.setConstant(kTRUE);
  shift.setConstant(kFALSE);

  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_back),RooArgSet(NSig,NBack));

  if(const_back_flag){
//    fsig.setConstant(kTRUE);
    c1.setConstant(kTRUE);
    c2.setConstant(kTRUE);
    shift.setConstant(kTRUE);
    pdf.fitTo(ds,Verbose(),Range(-0.15,de_fit_max),Timer(true));
  }
  else{
    c1.setConstant(kTRUE);
    c2.setConstant(kTRUE);
//    f.setConstant(kTRUE);
//    pdf.fitTo(ds,Verbose(),Range(de_fit_min,de_fit_max),Timer(true));
    pdf.fitTo(ds,Range(-0.15,de_fit_max),Timer(true));
  }

  RooAbsReal* intSig  = pdf_sig.createIntegral(de,NormSet(de),Range("signal"));
  RooAbsReal* intBack = pdf_back.createIntegral(de,NormSet(de),Range("signal"));
  const double nsig  = intSig->getVal()*NSig.getVal();
  const double nsig_err = intSig->getVal()*NSig.getError();
  const double nback = intBack->getVal()*NBack.getVal();
  const double nback_err = intBack->getVal()*NBack.getError();
  const double purity = nsig/(nsig+nback);
  cout << "Nsig  = " << nsig << " +- " << nsig_err << endl;
  cout << "Nback = " << nback << " +- " << nback_err << endl;
  cout << "Purity = " << purity << endl;
  /////////////
  //  Plots  //
  /////////////
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_back),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(kest),LineStyle(kDashed));
  pdf.plotOn(deFrame,LineWidth(2));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#DeltaE, Data","#DeltaE, Data",600,700);
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
  deFrame->GetXaxis()->SetLabelSize(0.05);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

//  TPaveText *pt = new TPaveText(0.5,0.6,0.98,0.9,"brNDC");
//  pt->SetFillColor(0);
//  pt->SetTextAlign(12);
//  out.str("");
//  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
//  pt->AddText(out.str().c_str());
//  out.str("");
//  out << "Signal region: (" << de_min << "," << de_max << ")";
//  pt->AddText(out.str().c_str());
//  out.str("");
//  out << "S: " << (int)(nsig+0.5) << " #pm " << (int)(nsig_err+0.5);
//  pt->AddText(out.str().c_str());
//  out.str("");
//  out << "B: " << (int)(nback+0.5) << " #pm " << (int)(nback_err+0.5) << " (" << gen_back << ")";
//  pt->AddText(out.str().c_str());
//  out.str("");
//  out << "Purity: " << purity*100.;
//  pt->AddText(out.str().c_str());
//  pt->Draw();
  
  TLine *de_lineLEFT = new TLine(de_min,0,de_min,40);
  de_lineLEFT->SetLineColor(kRed);
  de_lineLEFT->SetLineStyle(1);
  de_lineLEFT->Draw();
  
  TLine *de_lineRIGHT = new TLine(de_max,0,de_max,40);
  de_lineRIGHT->SetLineColor(kRed);
  de_lineRIGHT->SetLineStyle(1);
  de_lineRIGHT->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(-0.3,3,0.3,3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(-0.3,0,0.3,0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(-0.3,-3,0.3,-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
  
  ////////////
  // Signal //
  ////////////
  shift.setVal(0);
  RooPlot* deFrame_sig = de.frame();
  ds_sig.plotOn(deFrame_sig,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf_sig.plotOn(deFrame_sig,LineWidth(2));

  RooHist* hdepull_sig = deFrame_sig->pullHist();
  RooPlot* dePull_sig = de.frame(Title("#Delta E pull distribution"));
  dePull_sig->addPlotable(hdepull_sig,"P");
  dePull_sig->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm_sig = new TCanvas("#Delta E, Signal","#Delta E, Signal",600,700);
  cm_sig->cd();

  TPad *pad5 = new TPad("pad5","pad5",0.01,0.20,0.99,0.99);
  TPad *pad6 = new TPad("pad6","pad6",0.01,0.01,0.99,0.20);
  pad5->Draw();
  pad6->Draw();

  pad5->cd();
  pad5->SetLeftMargin(0.15);
  pad5->SetFillColor(0);

  deFrame_sig->GetXaxis()->SetTitleSize(0.05);
  deFrame_sig->GetXaxis()->SetTitleOffset(0.85);
  deFrame_sig->GetXaxis()->SetLabelSize(0.05);
  deFrame_sig->GetYaxis()->SetTitleOffset(1.6);
  deFrame_sig->Draw();

  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame_sig->chiSquare();
  TPaveText *pt_sig = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  pt_sig->SetFillColor(0);
  pt_sig->SetTextAlign(12);
  pt_sig->AddText(out.str().c_str());
  pt_sig->Draw();

  pad6->cd(); pad6->SetLeftMargin(0.15); pad6->SetFillColor(0);
  dePull_sig->SetMarkerSize(0.05); dePull_sig->Draw();
//  TLine *de_lineUP = new TLine(-0.3,3,0.3,3);
//  de_lineUP->SetLineColor(kBlue);
//  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
//  TLine *de_line = new TLine(-0.3,0,0.3,0);
//  de_line->SetLineColor(kBlue);
//  de_line->SetLineStyle(1);
//  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
//  TLine *de_lineDOWN = new TLine(-0.3,-3,0.3,-3);
//  de_lineDOWN->SetLineColor(kBlue);
//  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm_sig->Update();

  //////////////
  //  Gen MC  //
  //////////////
  RooPlot* deFrame_rho = de.frame();
  ds_rho.plotOn(deFrame_rho,DataError(RooAbsData::SumW2),MarkerSize(1));
  histpdf.plotOn(deFrame_rho,LineWidth(2));
  TCanvas* rho_cm = new TCanvas("rho_cm","rho_cm",600,400);
  rho_cm->cd();
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  deFrame_rho->GetXaxis()->SetTitleSize(0.05);
  deFrame_rho->GetXaxis()->SetTitleOffset(0.85);
  deFrame_rho->GetXaxis()->SetLabelSize(0.05);
  deFrame_rho->GetYaxis()->SetTitleOffset(1.6);
  deFrame_rho->Draw();
  
  RooPlot* deFrame_comb = de.frame();
  ds_back.plotOn(deFrame_comb,DataError(RooAbsData::SumW2),MarkerSize(1));
  cheb.plotOn(deFrame_comb,LineWidth(2));
  TCanvas* comb_cm = new TCanvas("comb_cm","comb_cm",600,400);
  comb_cm->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.00,0.99,0.99);
  pad2->Draw();
  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  deFrame_comb->GetXaxis()->SetTitleSize(0.05);
  deFrame_comb->GetXaxis()->SetTitleOffset(0.85);
  deFrame_comb->GetXaxis()->SetLabelSize(0.05);
  deFrame_comb->GetYaxis()->SetTitleOffset(1.6);
  deFrame_comb->Draw();

/*  RooHist* hdepull_back = deFrame_back->pullHist();
  RooPlot* dePull_back = de.frame(Title("#Delta E pull distribution"));
  dePull_back->addPlotable(hdepull,"P");
  dePull_back->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm_back = new TCanvas("#DeltaE, Genetic MC","#DeltaE, Genetic MC",600,700);
  cm_back->cd();

  TPad *pad7 = new TPad("pad7","pad7",0.01,0.20,0.99,0.99);
  TPad *pad8 = new TPad("pad8","pad8",0.01,0.01,0.99,0.20);
  pad7->Draw();
  pad8->Draw();

  pad7->cd();
  pad7->SetLeftMargin(0.15);
  pad7->SetFillColor(0);

  deFrame_back->GetXaxis()->SetTitleSize(0.05);
  deFrame_back->GetXaxis()->SetTitleOffset(0.85);
  deFrame_back->GetXaxis()->SetLabelSize(0.05);
  deFrame_back->GetYaxis()->SetTitleOffset(1.6);
  deFrame_back->Draw();
*/
  return;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame_back->chiSquare();
  TPaveText *pt_back = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
  pt_back->SetFillColor(0);
  pt_back->SetTextAlign(12);
  pt_back->AddText(out.str().c_str());
  pt_back->Draw();

  pad8->cd(); pad8->SetLeftMargin(0.15); pad8->SetFillColor(0);
  dePull_back->SetMarkerSize(0.05); dePull_back->Draw();
//  TLine *de_lineUP = new TLine(-0.3,3,0.3,3);
//  de_lineUP->SetLineColor(kBlue);
//  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
//  TLine *de_line = new TLine(-0.3,0,0.3,0);
//  de_line->SetLineColor(kBlue);
//  de_line->SetLineStyle(1);
//  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
//  TLine *de_lineDOWN = new TLine(-0.3,-3,0.3,-3);
//  de_lineDOWN->SetLineColor(kBlue);
//  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm_back->Update();
  
  ifile->Close();
  ifile_sig->Close();
  ifile_back->Close();
}
