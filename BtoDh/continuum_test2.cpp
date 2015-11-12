#include "cuts.h"

void continuum_test2(const int type = 11){
  TChain* gen_chain = new TChain("TEvent");
  gen_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed.root");
  TChain* sig_chain = new TChain("TEvent");
  TChain* bkg_chain = new TChain("TEvent");
  int nsig;
  int m_mode,m_h0mode;
  double eff_min = 0.4;
  double eff_max = 0.95;
  double signif_min = 0;
  double signif_max = 15;
  double rej_min = 0.8;
  double rej_max = 1.0;
//  double mh0_min, mh0_max;
  double de_min, de_max;
  double mbc_min, mbc_max;
  string label;
  switch(type){
  case 11:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcPi0_s7.root");
    m_mode = 1; m_h0mode = 10;
    label = string("#pi^{0}");
    signif_min = 13.; signif_max = 20;
    de_min  =-0.1;
    de_max  = 0.0921343;
    mbc_min = 5.27093;
    mbc_max = 5.28789;
    break;
  case 101:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    m_mode = 10; m_h0mode = 10;
    label = string("D^{*0}#pi^{0}");
    signif_min = 5.5; signif_max = 7.2;
    de_min  =-0.1;
    de_max  = 0.0647985;
    mbc_min = 5.27146;
    mbc_max = 5.28769;
    break;
  case 12:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s2.root");
    m_mode = 2; m_h0mode = 10;
    label = string("#eta#rightarrow#gamma#gamma");
    signif_min = 7.5; signif_max = 10;
    de_min  =-0.0994528;
    de_max  = 0.0758825;
    mbc_min = 5.27147;
    mbc_max = 5.28763;
    break;
  case 102:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    m_mode = 20; m_h0mode = 10;
    label = string("D^{*0}#eta(#rightarrow#gamma#gamma)");
    signif_min = 4.; signif_max = 5.;
    de_min  =-0.0853615;
    de_max  = 0.0630641;
    mbc_min = 5.27182;
    mbc_max = 5.28746;
    break;
  case 13:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s2.root");
    m_mode = 2; m_h0mode = 20;
    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    signif_min = 4.5; signif_max = 6.2;
    de_min  =-0.0534644;
    de_max  = 0.0453713;
    mbc_min = 5.27205;
    mbc_max = 5.28735;
    break;
  case 103:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    m_mode = 20; m_h0mode = 20;
    label = string("D^{*0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})");
    signif_min = 5; signif_max = 6.5;
    de_min  =-0.0567598;
    de_max  = 0.0459849;
    mbc_min = 5.27182;
    mbc_max = 5.28746;
    break;
  case 14:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcOMEGA_s5.root");
    m_mode = 3; m_h0mode = 20;
    label = string("#omega");
    signif_min = 11.5; signif_max = 15.0;
    de_min  =-0.0565915;
    de_max  = 0.0467748;
    mbc_min = 5.27198;
    mbc_max = 5.28743;
    break;
  case 15:
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETAP_s1.root");
    m_mode = 5; m_h0mode = 10;
    label = string("#eta'(#rightarrow#gamma#gamma)");
    signif_min = 3.; signif_max = 4;
    de_min  =-0.0873898;
    de_max  = 0.0642064;
    mbc_min = 5.27255;
    mbc_max = 5.28718;
    break;
  default:
      cout << "Wrong signal type " << type << endl;
      return;
  }

  bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_cont.root");
//  bkg_chain->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_charm_2_12.root");
//  bkg_chain->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_uds_2_12.root");

  RooArgSet argset;

  RooCategory mode("mode","mode");
  mode.defineType("mode",m_mode);

  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("h0mode",m_h0mode);

  argset.add(mode);
  argset.add(h0mode);

  RooCategory good_icpv("good_icpv","good_icpv");
  good_icpv.defineType("good_icpv",1);
  argset.add(good_icpv);

  RooRealVar mbc("mbc","mbc",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",de_min,de_max,"GeV"); argset.add(de);
  RooRealVar cos_thr("cos_thr","cos_thr",-1,1); argset.add(cos_thr);
  RooRealVar cos_hel("cos_hel","cos_hel",-1,1); if(type == 14) argset.add(cos_hel);
//  RooRealVar md("md_raw","md_raw",md_min,md_max,"GeV");     argset.add(md);
//  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mh0("mh0","mh0",mh0_min,mh0_max,"GeV"); argset.add(mh0);
//  RooRealVar mpi0("mpi0","mpi0",mpi0_min,mpi0_max,"GeV"); if(type != 12) argset.add(mpi0);
//  RooRealVar chi2_mass_d0("chi2_mass_d0","chi2_mass_d0",0.,50.); argset.add(chi2_mass_d0);

  const bool bdt_flag  = type<15;
  const bool bdtg_flag = false;
  const bool lh1_flag  = false;
  const bool lh0_flag  = type>=15;

  double bdt_best_cut  = 0;
  double bdtg_best_cut = 0;
  double lh1_best_cut  = 0;
  double lh0_best_cut  = 0;

  double bdt_best_sig  = 0;
  double bdtg_best_sig = 0;
  double lh1_best_sig  = 0;
  double lh0_best_sig  = 0;

  double bdt_best_eff  = 0;
  double bdtg_best_eff = 0;
  double lh1_best_eff  = 0;
  double lh0_best_eff  = 0;

  double bdt_best_sup  = 0;
  double bdtg_best_sup = 0;
  double lh1_best_sup  = 0;
  double lh0_best_sup  = 0;

  double bdt_best_pur  = 0;
  double bdtg_best_pur = 0;
  double lh1_best_pur  = 0;
  double lh0_best_pur  = 0;

  if(bdt_flag){  RooRealVar bdt("bdt","bdt",-1.,1.); argset.add(bdt);}
  if(bdtg_flag){ RooRealVar bdtg("bdtg","bdtg",-1.,1.); argset.add(bdtg);}
  if(lh1_flag){  RooRealVar lh1("lh1","lh1",0.,1.); argset.add(lh1);}
  if(lh0_flag){  RooRealVar lh0("lh0","lh0",0.,1.); argset.add(lh0);}
//  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);
//  RooRealVar cos_thr("cos_thr","cos_thr",-1.,1.); argset.add(cos_thr);

//  RooArgSet argset_gen(argset);

  RooCategory b0f("b0f","b0f");
  b0f.defineType("comb",-1);
  b0f.defineType("good",1);
  b0f.defineType("wrph",5);
  b0f.defineType("fsr",10);
  argset.add(b0f);

  RooDataSet ds_sig("ds_sig","ds_sig",sig_chain,argset,"b0f == 1 || b0f == 5 || b0f == 10");
  RooDataSet ds_con("ds_con","ds_con",bkg_chain,argset,"b0f == -1");


  RooDataSet ds_gen("ds_gen","ds_gen",gen_chain,argset,"b0f == 1 || b0f == 5 || b0f == 10");
  nsig = ds_gen.sumEntries()*0.25;
  cout << "Nsig: " << nsig << endl;
  ds_sig.Print();
  ds_con.Print();

  double sig_eff_bdt[100];
  double sig_eff_bdtg[100];
  double sig_eff_lh1[100];
  double sig_eff_lh0[100];

  double back_reg_bdt[100];
  double back_reg_bdtg[100];
  double back_reg_lh1[100];
  double back_reg_lh0[100];

  double signif_bdt[100];
  double signif_bdtg[100];
  double signif_lh1[100];
  double signif_lh0[100];

  const double gen = ds_con.sumEntries()*0.25;
  const double sig = ds_sig.sumEntries();

  stringstream out;
  double S,B;
  const int dots = 100;
  const double dbdt = 1./dots;
  cout << gen*4 << " " << sig << endl;
  for(int i=0; i<dots; i++){
    if(bdt_flag){
      out.str("");
      out << "bdt>" << i*dbdt;
      S = ds_sig.sumEntries(out.str().c_str());
      B = ds_con.sumEntries(out.str().c_str())*0.25;
      sig_eff_bdt[i]  = S/sig;
      back_reg_bdt[i] = (gen-B)/gen;
      S = nsig/sig*S;
//      B = nback/gen*B;
      signif_bdt[i] = S/sqrt(S+B+0.0001);
      if(bdt_best_sig<(signif_bdt[i]-0.07)){
        bdt_best_sig = signif_bdt[i];
        bdt_best_cut = i*dbdt;
        bdt_best_pur = S/(S+B);
        bdt_best_eff = sig_eff_bdt[i];
        bdt_best_sup = back_reg_bdt[i];
      }
    }

    if(lh1_flag){
      out.str("");
      out << "lh1>sqrt(1-" << i*dbdt*i*dbdt << ") && abs(cos_thr)<1.";
      if(type == 14) out << " && abs(cos_hel)>0.3";
      S = ds_sig.sumEntries(out.str().c_str());
      B = ds_con.sumEntries(out.str().c_str())*0.25;
      sig_eff_lh1[i]  = S/sig;
      back_reg_lh1[i] = (gen-B)/gen;
      S = nsig/sig*S;
//      B = nback/gen*B;
      signif_lh1[i] = S/sqrt(S+B+0.0001);
      if(lh1_best_sig<signif_lh1[i]){
        lh1_best_sig = signif_lh1[i];
        lh1_best_cut = sqrt(1-i*dbdt*i*dbdt);
        lh1_best_pur = S/(S+B);
        lh1_best_eff = sig_eff_lh1[i];
        lh1_best_sup = back_reg_lh1[i];
      }
    }

    if(lh0_flag){
      out.str("");
      out << "lh0>" << -1 +2*i*dbdt << " && abs(cos_thr)<1.";
      if(type == 14) out << " && abs(cos_hel)>0.3";
      S = ds_sig.sumEntries(out.str().c_str());
      B = ds_con.sumEntries(out.str().c_str())*0.25;
      sig_eff_lh0[i]  = S/sig;
      back_reg_lh0[i] = (gen-B)/gen;
      S = nsig/sig*S;
//      B = nback/gen*B;
      signif_lh0[i] = S/sqrt(S+B+0.0001);
      if(lh0_best_sig<signif_lh0[i]){
        lh0_best_sig = signif_lh0[i];
        lh0_best_cut = sqrt(1-i*dbdt*i*dbdt);
        lh0_best_pur = S/(S+B);
        lh0_best_eff = sig_eff_lh0[i];
        lh0_best_sup = back_reg_lh0[i];
      }
    }
    
    if(bdtg_flag){
      out.str("");
      out << "bdtg>sqrt(1-" << i*dbdt*i*dbdt << ")";
      S = ds_sig.sumEntries(out.str().c_str());
      B = ds_con.sumEntries(out.str().c_str())/0.949;
      sig_eff_bdtg[i]  = S/sig;
      back_reg_bdtg[i] = (gen-B)/gen;
      S = nsig/sig*S;
//      B = nback/gen*B;
      signif_bdtg[i] = S/sqrt(S+B+0.0001);
      if(bdtg_best_sig<signif_bdtg[i]){
        bdtg_best_sig = signif_bdtg[i];
        bdtg_best_cut = sqrt(1-i*dbdt*i*dbdt);
        bdtg_best_pur = S/(S+B);
        bdtg_best_eff = sig_eff_bdtg[i];
        bdtg_best_sup = back_reg_bdtg[i];
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1","Signal Efficiency vs. Background Rejection",600,600);
  TH1F* h1 = c1->Draw();
  c1->SetGrid();
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Background Rejection");
  TGraph* gr_bdt = new TGraph(dots,sig_eff_bdt,back_reg_bdt);
  gr_bdt->SetMarkerSize(1);
  gr_bdt->SetMarkerColor(kBlue);
  gr_bdt->SetMarkerStyle(21);
  gr_bdt->SetTitle("Background Rejection");
  gr_bdt->GetXaxis()->SetRangeUser(eff_min,eff_max);
  gr_bdt->GetXaxis()->SetTitle("Signal Efficiency");
  gr_bdt->GetXaxis()->SetTitleSize(0.06);
  gr_bdt->GetXaxis()->SetTitleOffset(0.75);
  gr_bdt->GetYaxis()->SetRangeUser(rej_min,rej_max);

  TGraph* gr_bdtg = new TGraph(dots,sig_eff_bdtg,back_reg_bdtg);
  gr_bdtg->SetMarkerSize(1);
  gr_bdtg->SetMarkerColor(kRed);
  gr_bdtg->SetMarkerStyle(21);
  
  TGraph* gr_ksfw1 = new TGraph(dots,sig_eff_lh1,back_reg_lh1);
  gr_ksfw1->SetMarkerSize(1);
  gr_ksfw1->SetMarkerColor(kRed);
  gr_ksfw1->SetMarkerStyle(22);

  TGraph* gr_ksfw0 = new TGraph(dots,sig_eff_lh0,back_reg_lh0);
  gr_ksfw0->SetMarkerSize(1.4);
  gr_ksfw0->SetMarkerColor(kBlue);
  gr_ksfw0->SetMarkerStyle(22);
  gr_ksfw0->SetTitle("Background Rejection");
  gr_ksfw0->GetXaxis()->SetRangeUser(eff_min,eff_max);
  gr_ksfw0->GetXaxis()->SetTitle("Signal Efficiency");
  gr_ksfw0->GetXaxis()->SetTitleSize(0.06);
  gr_ksfw0->GetXaxis()->SetTitleOffset(0.75);
  gr_ksfw0->GetYaxis()->SetRangeUser(rej_min,rej_max);

  if(bdt_flag)mg->Add(gr_bdt);
  if(bdtg_flag)mg->Add(gr_bdtg);
  if(lh1_flag)mg->Add(gr_ksfw1);
  if(lh0_flag)mg->Add(gr_ksfw0);

  if(bdt_flag && !bdtg_flag && !lh1_flag && !lh0_flag){
    gr_bdt->Draw("AP");
  } else if(lh0_flag && !bdtg_flag && !lh1_flag && !bdt_flag){
    gr_ksfw0->Draw("AP");
  } else{
    mg->Draw("AP");
  }

  if(type<15){
    TLine* cutline = new TLine(bdt_best_eff,rej_min,bdt_best_eff,rej_max);
    cutline->SetLineColor(kRed);
    cutline->SetLineWidth(2);
    cutline->Draw("AP");
  } else{
    TLine* cutline = new TLine(lh0_best_eff,rej_min,lh0_best_eff,rej_max);
    cutline->SetLineColor(kRed);
    cutline->SetLineWidth(2);
    cutline->Draw("AP");
  }

  TPaveText *pt = new TPaveText(0.7,0.75,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(label.c_str());
  out.str("");
  out << "Input S = " << nsig;
  pt->AddText(out.str().c_str());
  out.str("");
  out << "Input B = " << (int)gen;
  pt->AddText(out.str().c_str());
  pt->Draw();

  TPaveText *pt2 = new TPaveText(0.3,0.25,0.6,0.5,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  pt2->AddText("Best cut:");
  if(type<15){
    out.str("");
    out << setprecision(3);
    out << "Pur = " << bdt_best_pur*100;
    pt2->AddText(out.str().c_str());
    out.str("");
    out << "Eff = " << bdt_best_eff*100;
    pt2->AddText(out.str().c_str());
    out.str("");
    out << "Rej = " << bdt_best_sup*100;
    pt2->AddText(out.str().c_str());
    pt2->Draw();
  } else{
    out.str("");
    out << setprecision(3);
    out << "Pur = " << lh0_best_pur*100;
    pt2->AddText(out.str().c_str());
    out.str("");
    out << "Eff = " << lh0_best_eff*100;
    pt2->AddText(out.str().c_str());
    out.str("");
    out << "Rej = " << lh0_best_sup*100;
    pt2->AddText(out.str().c_str());
    pt2->Draw();
  }

  c1->Update();
  out.str("");
  out << "../Note/pics/eff-rej-m" << type << ".root";
  c1->Print(out.str().c_str());
  out.str("");
  out << "../Note/pics/eff-rej-m" << type << ".eps";
  c1->Print(out.str().c_str());

  TCanvas* c2 = new TCanvas("c2","Significance",600,600);
  c2->cd();
  c2->SetGrid();

  TMultiGraph *signif_mg = new TMultiGraph();
  signif_mg->SetTitle("Signal Significance");
  TGraph* signif_gr_bdt = new TGraph(dots,sig_eff_bdt,signif_bdt);
  signif_gr_bdt->SetMarkerSize(1);
  signif_gr_bdt->SetMarkerColor(kBlue);
  signif_gr_bdt->SetMarkerStyle(21);
  signif_gr_bdt->SetTitle("Signal Significance");
  signif_gr_bdt->GetXaxis()->SetRangeUser(eff_min,eff_max);
  signif_gr_bdt->GetXaxis()->SetTitle("Signal Efficiency");
  signif_gr_bdt->GetXaxis()->SetTitleSize(0.06);
  signif_gr_bdt->GetXaxis()->SetTitleOffset(0.75);
  signif_gr_bdt->GetYaxis()->SetRangeUser(signif_min,signif_max);

  TGraph* signif_gr_bdtg = new TGraph(dots,sig_eff_bdtg,signif_bdtg);
  signif_gr_bdtg->SetMarkerSize(1);
  signif_gr_bdtg->SetMarkerColor(kRed);
  signif_gr_bdtg->SetMarkerStyle(21);

  TGraph* signif_gr_ksfw1 = new TGraph(dots,sig_eff_lh1,signif_lh1);
  signif_gr_ksfw1->SetMarkerSize(1);
  signif_gr_ksfw1->SetMarkerColor(kRed);
  signif_gr_ksfw1->SetMarkerStyle(22);

  TGraph* signif_gr_ksfw0 = new TGraph(dots,sig_eff_lh0,signif_lh0);
  signif_gr_ksfw0->SetMarkerSize(1.4);
  signif_gr_ksfw0->SetMarkerColor(kBlue);
  signif_gr_ksfw0->SetMarkerStyle(22);
  signif_gr_ksfw0->SetTitle("Signal Significance");
  signif_gr_ksfw0->GetXaxis()->SetRangeUser(eff_min,eff_max);
  signif_gr_ksfw0->GetXaxis()->SetTitle("Signal Efficiency");
  signif_gr_ksfw0->GetXaxis()->SetTitleSize(0.06);
  signif_gr_ksfw0->GetXaxis()->SetTitleOffset(0.75);
  signif_gr_ksfw0->GetYaxis()->SetRangeUser(signif_min,signif_max);

  if(bdt_flag)  signif_mg->Add(signif_gr_bdt);
  if(bdtg_flag) signif_mg->Add(signif_gr_bdtg);
  if(lh1_flag)  signif_mg->Add(signif_gr_ksfw1);
  if(lh0_flag)  signif_mg->Add(signif_gr_ksfw0);

  if(bdt_flag && !bdtg_flag && !lh1_flag && !lh0_flag){
    signif_gr_bdt->Draw("AP");
  } else if(!bdt_flag && !bdtg_flag && !lh1_flag && lh0_flag){
    signif_gr_ksfw0->Draw("AP");
  } else{
    signif_mg->Draw("AP");
  }
  if(type<15){
    TLine* cutline2 = new TLine(bdt_best_eff,signif_min,bdt_best_eff,signif_max);
    cutline2->SetLineColor(kRed);
    cutline2->SetLineWidth(2);
    cutline2->Draw("AP");
  } else{
    TLine* cutline2 = new TLine(lh0_best_eff,signif_min,lh0_best_eff,signif_max);
    cutline2->SetLineColor(kRed);
    cutline2->SetLineWidth(2);
    cutline2->Draw("AP");
  }

  pt->Draw();
  pt2->Draw();

  c2->Update();
  out.str("");
  out << "../Note/pics/eff-fom-m" << type << ".root";
  c2->Print(out.str().c_str());
  out.str("");
  out << "../Note/pics/eff-fom-m" << type << ".eps";
  c2->Print(out.str().c_str());
  
  cout << "Summary:" << endl;
  
  if(bdtg_flag){
    cout << "bdtg decision: ";
    cout << "sig: "   << bdtg_best_sig;
    cout << ", cut: " << bdtg_best_cut;
    cout << ", pur: " << bdtg_best_pur;
    cout << ", eff: " << bdtg_best_eff;
    cout << ", sup: " << bdtg_best_sup << endl;
  }
  if(bdt_flag){
    cout << "bdt decision:  ";
    cout << "sig: "   << bdt_best_sig;
    cout << ", cut: " << bdt_best_cut;
    cout << ", pur: " << bdt_best_pur;
    cout << ", eff: " << bdt_best_eff;
    cout << ", sup: " << bdt_best_sup << endl;
  }
  if(lh1_flag){
    cout << "lh1 decision:  ";
    cout << "sig: "   << lh1_best_sig;
    cout << ", cut: " << lh1_best_cut;
    cout << ", pur: " << lh1_best_pur;
    cout << ", eff: " << lh1_best_eff;
    cout << ", sup: " << lh1_best_sup << endl;
  }
  if(lh0_flag){
    cout << "lh0 decision:  ";
    cout << "sig: "   << lh0_best_sig;
    cout << ", cut: " << lh0_best_cut;
    cout << ", pur: " << lh0_best_pur;
    cout << ", eff: " << lh0_best_eff;
    cout << ", sup: " << lh0_best_sup << endl;
  }

  return;
}
