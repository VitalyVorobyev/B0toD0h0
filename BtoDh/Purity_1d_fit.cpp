#include "cuts.h"
using namespace RooFit;

void Purity_1d_fit(int type = 0){
  TChain* tree = new TChain("TEvent");
  if(!type) tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1_full.root");
  else      tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_data.root");

  RooCategory b0f("b0f","b0f");
  b0f.defineType("signal",1);
  b0f.defineType("fsr",10);
  b0f.defineType("bad_pi0",5);
  b0f.defineType("rho",3);
  b0f.defineType("comb",-1);

  RooArgSet argset;

  const double deMin = -0.15;
  const double deMax = 0.3;

  RooRealVar mbc("mbc","M_{bc}",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",deMin,deMax,"GeV"); argset.add(de);
  de.setRange("Signal",de_min,de_max);
  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
  RooRealVar bdtgs("bdtgs","bdtgs",bdtgs_cut,1.); argset.add(bdtgs);
  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  if(!type) argset.add(b0f);

  RooDataSet ds("ds","ds",tree,argset,"mbc>0||mbc<=0");
//  RooDataSet* ds0 = ds.reduce(RooArgSet(de));
  
  stringstream out;
  if(!type){
    out.str("");
    out << "de<" << de_max << " && de>" << de_min;
    Roo1DTable* sigtable = ds.table(b0f,out.str().c_str());
    sigtable->Print();
    sigtable->Print("v");

    Roo1DTable* fulltable = ds.table(b0f);
    fulltable->Print();
    fulltable->Print("v");
  }

//  RooDataHist* dh = ds0->binnedClone();

//  ds0->Print();

  ////////////////
  // Signal PDF //
  ////////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar de0("de0","de0",m_de0,-0.1,0.1); if(cSig) de0.setConstant(kTRUE);
  RooRealVar s1("s1","s1",m_s1,0.,0.5); if(cSig) s1.setConstant(kTRUE);
  RooGaussian g1("g1","g1",de,de0,s1);

  RooRealVar deCBl("deCBl","deCBl",m_deCBl,-0.1,0.1); if(cSig) deCBl.setConstant(kTRUE);
  RooRealVar sCBl("sCBl","sCBl",m_sCBl,0.,0.5); if(cSig) sCBl.setConstant(kTRUE);
  RooRealVar nl("nl","nl",m_nl,0.,100.); if(cSig) nl.setConstant(kTRUE);
  RooRealVar alphal("alphal","alphal",m_alphal,-10.,10.); if(cSig) alphal.setConstant(kTRUE);

  RooRealVar deCBr("deCBr","deCBr",m_deCBr,-0.1,0.1); if(cSig) deCBr.setConstant(kTRUE);
  RooRealVar sCBr("sCBr","sCBr",m_sCBr,0.,0.5); if(cSig) sCBr.setConstant(kTRUE);
  RooRealVar nr("nr","nr",m_nr,0.,100.); if(cSig) nr.setConstant(kTRUE);
  RooRealVar alphar("alphar","alphar",m_alphar,-10.,10.); if(cSig) alphar.setConstant(kTRUE);

  RooCBShape CBl("CBl","CBl",de,deCBl,sCBl,alphal,nl);
  RooCBShape CBr("CBr","CBr",de,deCBr,sCBr,alphar,nr);

  RooRealVar fCBl("fCBl","fCBl",m_fCBl,0.,1.); if(cSig) fCBl.setConstant(kTRUE);
  RooRealVar fCBr("fCBr","fCBr",m_fCBr,0.,1.); if(cSig) fCBr.setConstant(kTRUE);

  RooAddPdf pdf_sig("pdf_sig","pdf_sig",RooArgList(CBl,CBr,g1),RooArgSet(fCBl,fCBr));

  //////////////
  // Comb PDF //
  //////////////
  ////////////
  // de pdf //
  ////////////
  RooRealVar c1("c1","c1",mc_c1_1d,-10.,10.); if(cComb) c1.setConstant(kTRUE);
  RooRealVar c2("c2","c2",mc_c2_1d,-10.,10.); if(cComb) c2.setConstant(kTRUE);
  RooChebychev pdf_comb("pdf_comb","pdf_comb",de,RooArgSet(c1,c2));

  /////////////
  // Rho PDF //
  /////////////
  ////////////
  // de pdf //
  ////////////
if(de_rho_param == 0){
  RooRealVar exppar("exppar","exppar",mr_exppar,-40.,-25.);// if(cRho) exppar.setConstant(kTRUE);
  RooExponential pdf_rho("pdf_rho","pdf_rho",de,exppar);
  }
  
  RooRealVar de0r("de0r","de0r",mr_de0r,-0.2,0.12); if(cRho) de0r.setConstant(kTRUE);
  
  if(de_rho_param == 1){
   RooRealVar slopel("slopel","slopel",mr_slopel,-1000,-500.); if(cRho) slopel.setConstant(kTRUE);
   RooRealVar sloper("sloper","sloper",mr_sloper,-10000,0.); if(cRho) sloper.setConstant(kTRUE);
   RooRealVar steep("steep","steep",mr_steep,7.,9.); if(cRho) steep.setConstant(kTRUE);
   RooRealVar p5("p5","p5",mr_p5,0.01,1000.); if(cRho) p5.setConstant(kTRUE);
   RooRhoDeltaEPdf pdf_rho("pdf_rho","pdf_rho",de,de0r,slopel,sloper,steep,p5);
  }
  
  if(de_rho_param == -1){
   RooRealVar x0("x0","x0",mr_x0_1d,-0.2,0.12); if(cRho) x0.setConstant(kTRUE);
   RooRealVar p1("p1","p1",mr_p1_1d,-1000.,100.); if(cRho) p1.setConstant(kTRUE);
   RooRealVar p2("p2","p2",mr_p2_1d,0.,100.); if(cRho) p2.setConstant(kTRUE);
   RooGenericPdf pdf_rho("pdf_rho","1+@0*@1-@2*TMath::Log(1+TMath::Exp(@2*(@0-@1)/@3))",RooArgSet(de,x0,p1,p2));
  }
  //////////////////
  // Complete PDF //
  //////////////////
  RooRealVar Nsig("Nsig","Nsig",700,100.,1500.);// fsig.setConstant(kTRUE);
  RooRealVar Nrho("Nrho","Nrho",400,100,1500.);// frho.setConstant(kTRUE);
  RooRealVar Ncmb("Ncmb","Ncmb",1000,100,100000);// frho.setConstant(kTRUE);
  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,pdf_rho,pdf_comb),RooArgList(Nsig,Nrho,Ncmb));

  pdf.fitTo(ds,Verbose(),Timer(true));

   RooAbsReal* intSig  = pdf_sig.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
   RooAbsReal* intRho  = pdf_rho.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
   RooAbsReal* intCmb  = pdf_comb.createIntegral(RooArgSet(de),NormSet(RooArgSet(de)),Range("Signal"));
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
   cout << "Nsig = " << nsig <<" +- " << nsig_err << endl;
   cout << "Nrho = " << nrho <<" +- " << nrho_err << endl;
   cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << endl;
   
  /////////////
  //  Plots  //
  /////////////
  // de //
  RooPlot* deFrame = de.frame();
  ds.plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
  pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_rho),LineStyle(kDashed));
  pdf.plotOn(deFrame,Components(pdf_comb),LineStyle(kDashed));
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

  stringstream out1;
  TPaveText *pt = new TPaveText(0.6,0.75,0.98,0.9,"brNDC");
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
  
  if(!type){
    out.str("");
    out << "de<" << de_max << " && de>" << de_min;
    Roo1DTable* sigtable = ds.table(b0f,out.str().c_str());
    sigtable->Print();
    sigtable->Print("v");
    
    Roo1DTable* fulltable = ds.table(b0f);
    fulltable->Print();
    fulltable->Print("v");
  }
  
  cout << "Nsig = " << nsig <<" +- " << nsig_err << " +- " << nsig_err_npq << " (" << nsig_err_total << ")" << endl;
  cout << "Nrho = " << nrho <<" +- " << nrho_err << " +- " << nrho_err_npq << " (" << nrho_err_total << ")" << endl;
  cout << "Ncmb = " << ncmb <<" +- " << ncmb_err << " +- " << ncmb_err_npq << " (" << ncmb_err_total << ")" << endl;
  cout << "Pury = " << purity << " +- " << purity_err << endl;
}
