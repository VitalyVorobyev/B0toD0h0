#include "cuts.h"

void continuum_test(const int type = 11,const double nsig = 500, const double nback = 3000){
  TFile *ifile_sig;
  switch(type){
    case 11:
      ifile_sig = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigPi0.root");
      break;
    case 111:
      ifile_sig = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sig.root");
      break;
    case 12:
      ifile_sig = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigEta.root");
      break;
    case 13:
      ifile_sig = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigOmega.root");
      break;
    default:
      cout << "Wrong signal type " << type << endl;
      return;
  }

  TFile *ifile_gen = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_cont.root");

  TTree *tree_sig = (TTree*)ifile_sig->Get("TEvent");
  TTree *tree_gen = (TTree*)ifile_gen->Get("TEvent");

  RooArgSet argset;
  RooRealVar mbc("mbc","mbc",mbc_min,mbc_max,"GeV"); argset.add(mbc);
  RooRealVar de("de","#DeltaE",de_min,de_max,"GeV"); argset.add(de);
//  RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
//  RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
//  RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);

  const bool bdt_flag      = false;
  const bool bdtg_flag     = false;
  const bool bdts_flag     = false;
  const bool bdtsl_flag    = false;
  const bool bdtlh_flag    = false;
  const bool bdtglh_flag   = false;
  const bool bdtlhs_flag   = false;
  const bool bdtglhs_flag  = false;
  const bool nn_flag       = false;

  const bool bdtgsl_flag   = true;
  const bool bdtgdes_flag  = false;
  const bool bdtgmbcs_flag = false;
  const bool bdtgfr_flag   = false;
  const bool bdtgs_flag    = false;
  const bool lh_flag       = true;
  const bool lh1_flag      = false;

  double lh_best_cut       = 0;
  double lh1_best_cut      = 0;
  double bdtgs_best_cut    = 0;
  double bdtgsl_best_cut   = 0;
  double bdtgmbcs_best_cut = 0;
  double bdtgdes_best_cut  = 0;
  double bdtgfr_best_cut   = 0;

  double lh_best_sig       = 0;
  double lh1_best_sig      = 0;
  double bdtgs_best_sig    = 0;
  double bdtgsl_best_sig   = 0;
  double bdtgmbcs_best_sig = 0;
  double bdtgdes_best_sig  = 0;
  double bdtgfr_best_sig   = 0;

  double lh_best_eff       = 0;
  double lh1_best_eff      = 0;
  double bdtgs_best_eff    = 0;
  double bdtgsl_best_eff   = 0;
  double bdtgmbcs_best_eff = 0;
  double bdtgdes_best_eff  = 0;
  double bdtgfr_best_eff   = 0;

  double lh_best_sup       = 0;
  double lh1_best_sup      = 0;
  double bdtgs_best_sup    = 0;
  double bdtgsl_best_sup   = 0;
  double bdtgmbcs_best_sup = 0;
  double bdtgdes_best_sup  = 0;
  double bdtgfr_best_sup   = 0;

  double lh_best_pur       = 0;
  double lh1_best_pur      = 0;
  double bdtgs_best_pur    = 0;
  double bdtgsl_best_pur   = 0;
  double bdtgmbcs_best_pur = 0;
  double bdtgdes_best_pur  = 0;
  double bdtgfr_best_pur   = 0;

  if(bdt_flag){    RooRealVar bdt("bdt","bdt",-1.,1.); argset.add(bdt);}
  if(bdtg_flag){   RooRealVar bdtg("bdtg","bdtg",-1.,1.); argset.add(bdtg);}
  if(bdts_flag){   RooRealVar bdts("bdts","bdts",-1.,1.); argset.add(bdts);}
  if(bdtgs_flag){  RooRealVar bdtgs("bdtgs","bdtgs",-1.,1.); argset.add(bdtgs);}
  if(bdtsl_flag){  RooRealVar bdtsl("bdtsl","bdtsl",-1.,1.); argset.add(bdtsl);}
  if(bdtgsl_flag){ RooRealVar bdtgsl("bdtgsl","bdtgsl",-1.,1.); argset.add(bdtgsl);}
  if(lh_flag){     RooRealVar lh("lh","lh",0.,1.); argset.add(lh);}
  if(lh1_flag){    RooRealVar lh1("lh1","lh1",0.,1.); argset.add(lh1);}
  if(bdtlh_flag){  RooRealVar bdtlh("bdtlh","bdtlh",-1.,1.); argset.add(bdtlh);}
  if(bdtglh_flag){ RooRealVar bdtglh("bdtglh","bdtglh",-1.,1.); argset.add(bdtglh);}
  if(bdtlhs_flag){ RooRealVar bdtlhs("bdtlhs","bdtlhs",-1.,1.); argset.add(bdtlhs);}
  if(bdtglhs_flag){RooRealVar bdtglhs("bdtglhs","bdtglhs",-1.,1.); argset.add(bdtglhs);}

  if(bdtgfr_flag){  RooRealVar bdtgfr("bdtgfr","bdtgfr",-1.,1.); argset.add(bdtgfr);}
  if(bdtgmbcs_flag){RooRealVar bdtgmbcs("bdtgmbcs","bdtgmbcs",-1.,1.); argset.add(bdtgmbcs);}
  if(bdtgdes_flag){ RooRealVar bdtgdes("bdtgdes","bdtgdes",-1.,1.); argset.add(bdtgdes);}
  if(nn_flag){      RooRealVar nn("nn","nn",-1.,1.);}

  RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

  RooArgSet argset_gen(argset);

  RooDataSet ds_sig("ds_sig","ds_sig",tree_sig,argset_gen);
  RooDataSet ds_gen("ds_gen","ds_gen",tree_gen,argset_gen);
  if(nn_flag){
    TFile *innfile_sig = TFile::Open("/home/vitaly/B0toDh0/Tuples/nn_sig.root");
    TFile *innfile_gen = TFile::Open("/home/vitaly/B0toDh0/Tuples/nn_cont.root");
    TTree *treenn_sig = (TTree*)innfile_sig->Get("TEvent");
    TTree *treenn_cont = (TTree*)innfile_gen->Get("TEvent");
    RooDataSet ds_nnsig("ds_nnsig","ds_nnsig",treenn_sig,RooArgSet(mbc,de,nn));
    RooDataSet ds_nncont("ds_nncont","ds_nncont",treenn_cont,RooArgSet(mbc,de,nn));
    ds_nnsig.Print();
    ds_nncont.Print();
  }

  ds_sig.Print();
  ds_gen.Print();

  double sig_eff_bdt[100];
  double sig_eff_bdtg[100];
  double sig_eff_bdts[100];
  double sig_eff_bdtgs[100];
  double sig_eff_bdtgmbcs[100];
  double sig_eff_bdtgdes[100];
  double sig_eff_bdtgfr[100];
  double sig_eff_bdtsl[100];
  double sig_eff_bdtgsl[100];
  double sig_eff_ksfw[100];
  double sig_eff_ksfw1[100];
  double sig_eff_nn[100];
  double back_reg_bdt[100];
  double back_reg_bdtg[100];
  double back_reg_bdts[100];
  double back_reg_bdtgmbcs[100];
  double back_reg_bdtgdes[100];
  double back_reg_bdtgfr[100];
  double back_reg_bdtgs[100];
  double back_reg_bdtsl[100];
  double back_reg_bdtgsl[100];
  double back_reg_ksfw[100];
  double back_reg_ksfw1[100];
  double back_reg_nn[100];
  double signif_bdt[100];
  double signif_bdtg[100];
  double signif_bdts[100];
  double signif_bdtgs[100];
  double signif_bdtgmbcs[100];
  double signif_bdtgdes[100];
  double signif_bdtgfr[100];
  double signif_bdtsl[100];
  double signif_bdtgsl[100];
  double signif_ksfw[100];
  double signif_ksfw1[100];
  double signif_nn[100];

  const double gen = ds_gen.sumEntries();
  const double sig = ds_sig.sumEntries();
//  const double gennn = ds_nncont.sumEntries();
//  const double signn = ds_nnsig.sumEntries();

  stringstream out;
  out.str("");
  out << "bdtgsl_" << nsig << "_" << nback << ".txt";
  ofstream ofile(out.str().c_str(),ofstream::out);
  double S,B;
  const int dots = 100;
  const double dbdt = 1./dots;
  const double dlh  = 2./dots;
  cout << gen << " " << sig << endl;// << " " << gennn << " " << signn << endl;
  for(int i=0; i<dots; i++){
    if(bdt_flag){
    out.str("");
    out << "bdt>" << i*dbdt;
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdt[i]  = S/sig;
    back_reg_bdt[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdt[i] = S/sqrt(S+B+0.0001);
    }

    if(lh_flag){
    out.str("");
    out << "lh>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_ksfw[i]  = S/sig;
    back_reg_ksfw[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_ksfw[i] = S/sqrt(S+B+0.0001);
    if(lh_best_sig<signif_ksfw[i]){
      lh_best_sig = signif_ksfw[i];
      lh_best_cut = sqrt(1-i*dbdt*i*dbdt);
      lh_best_pur = S/(S+B);
      lh_best_eff = sig_eff_ksfw[i];
      lh_best_sup = back_reg_ksfw[i];
    }
    }
    
    if(lh1_flag){
    out.str("");
    out << "lh1>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_ksfw1[i]  = S/sig;
    back_reg_ksfw1[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_ksfw1[i] = S/sqrt(S+B+0.0001);
    if(lh1_best_sig<signif_ksfw1[i]){
      lh1_best_sig = signif_ksfw1[i];
      lh1_best_cut = sqrt(1-i*dbdt*i*dbdt);
      lh1_best_pur = S/(S+B);
      lh1_best_eff = sig_eff_ksfw1[i];
      lh1_best_sup = back_reg_ksfw1[i];
    }
    }

    if(bdtg_flag){
    out.str("");
    out << "bdtg>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtg[i]  = S/sig;
    back_reg_bdtg[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtg[i] = S/sqrt(S+B+0.0001);
    }

    if(bdts_flag){
    out.str("");
    out << "bdts>" << i*dbdt;
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdts[i]  = S/sig;
    back_reg_bdts[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdts[i] = S/sqrt(S+B+0.0001);
    }

    if(bdtgs_flag){
    out.str("");
    out << "bdtgs>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtgs[i]  = S/sig;
    back_reg_bdtgs[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtgs[i] = S/sqrt(S+B+0.0001);
    if(bdtgs_best_sig<signif_bdtgs[i]){
      bdtgs_best_sig = signif_bdtgs[i];
      bdtgs_best_cut = sqrt(1-i*dbdt*i*dbdt);
      bdtgs_best_pur = S/(S+B);
      bdtgs_best_eff = sig_eff_bdtgs[i];
      bdtgs_best_sup = back_reg_bdtgs[i];
    }
    }

    if(bdtgmbcs_flag){
    out.str("");
    out << "bdtgmbcs>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtgmbcs[i]  = S/sig;
    back_reg_bdtgmbcs[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtgmbcs[i] = S/sqrt(S+B+0.0001);
    if(bdtgmbcs_best_sig<signif_bdtgmbcs[i]){
      bdtgmbcs_best_sig = signif_bdtgmbcs[i];
      bdtgmbcs_best_cut = sqrt(1-i*dbdt*i*dbdt);
      bdtgmbcs_best_pur = S/(S+B);
      bdtgmbcs_best_eff = sig_eff_bdtgmbcs[i];
      bdtgmbcs_best_sup = back_reg_bdtgmbcs[i];
    }
    }

    if(bdtgdes_flag){
    out.str("");
    out << "bdtgdes>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtgdes[i]  = S/sig;
    back_reg_bdtgdes[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtgdes[i] = S/sqrt(S+B+0.0001);
    if(bdtgdes_best_sig<signif_bdtgdes[i]){
      bdtgdes_best_sig = signif_bdtgdes[i];
      bdtgdes_best_cut = sqrt(1-i*dbdt*i*dbdt);
      bdtgdes_best_pur = S/(S+B);
      bdtgdes_best_eff = sig_eff_bdtgdes[i];
      bdtgdes_best_sup = back_reg_bdtgdes[i];
    }
    }

    if(bdtgfr_flag){
    out.str("");
    out << "bdtgfr>(-1+2*" << i*dbdt*i*dbdt << ")";
//    out << "bdtgfr>(0.0+0.009*" << i << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtgfr[i]  = S/sig;
    back_reg_bdtgfr[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtgfr[i] = S/sqrt(S+B+0.0001);
    if(bdtgfr_best_sig<signif_bdtgfr[i]){
      bdtgfr_best_sig = signif_bdtgfr[i];
      bdtgfr_best_cut = (-1+2*i*dbdt*i*dbdt);
      bdtgfr_best_pur = S/(S+B);
      bdtgfr_best_eff = sig_eff_bdtgfr[i];
      bdtgfr_best_sup = back_reg_bdtgfr[i];
    }
    }

    if(bdtsl_flag){
    out.str("");
    out << "bdtsl>" << i*dbdt;
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtsl[i]  = S/sig;
    back_reg_bdtsl[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtsl[i] = S/sqrt(S+B+0.0001);
    }

    if(bdtgsl_flag){
    out.str("");
    out << "bdtgsl>sqrt(1-" << i*dbdt*i*dbdt << ")";
    S = ds_sig.sumEntries(out.str().c_str());
    B = ds_gen.sumEntries(out.str().c_str());
    sig_eff_bdtgsl[i]  = S/sig;
    back_reg_bdtgsl[i] = (gen-B)/gen;
    S = nsig/sig*S;
    B = nback/gen*B;
    signif_bdtgsl[i] = S/sqrt(S+B+0.0001);
    if(bdtgsl_best_sig<signif_bdtgsl[i]){
      bdtgsl_best_sig = signif_bdtgsl[i];
      bdtgsl_best_cut = sqrt(1-i*dbdt*i*dbdt);
      bdtgsl_best_pur = S/(S+B);
      bdtgsl_best_eff = sig_eff_bdtgsl[i];
      bdtgsl_best_sup = back_reg_bdtgsl[i];
    }
//    ofile << sqrt(1-i*dbdt*i*dbdt) << " " << S << " " << B << " " << S/(S+B) << " " << signif_bdtgsl[i] << " " << sig_eff_bdtgsl[i] << " " << back_reg_bdtgsl[i] << endl;
    }

    if(nn_flag){
    out.str("");
    out << "nn>(-1+" << i*dlh << ")";
    S = ds_nnsig.sumEntries(out.str().c_str());
    B = ds_nncont.sumEntries(out.str().c_str());
    sig_eff_nn[i]  = S/signn;
    back_reg_nn[i] = (gennn-B)/gennn;
    S = nsig/signn*S;
    B = nback/gennn*B;
    signif_nn[i] = S/sqrt(S+B+0.0001);
    }
  }

  ofile.close();
  TCanvas* c1 = new TCanvas("c1","Signal Efficiency vs. Background Rejection",600,600);
  c1->SetGrid();
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Background Rejection");
//  mg->GetXaxis()->SetRangeUser(0.50,0.95);
//  mg->GetYaxis()->SetRangeUser(0.8,1.0);
  TGraph* gr_bdt = new TGraph(dots,sig_eff_bdt,back_reg_bdt);
  gr_bdt->SetMarkerSize(1);
  gr_bdt->SetMarkerColor(kBlue);
  gr_bdt->SetMarkerStyle(21);

  TGraph* gr_bdtg = new TGraph(dots,sig_eff_bdtg,back_reg_bdtg);
  gr_bdtg->SetMarkerSize(1);
  gr_bdtg->SetMarkerColor(kRed);
  gr_bdtg->SetMarkerStyle(21);

  TGraph* gr_bdts = new TGraph(dots,sig_eff_bdts,back_reg_bdts);
  gr_bdts->SetMarkerSize(1);
  gr_bdts->SetMarkerColor(kBlue);
  gr_bdts->SetMarkerStyle(20);

  TGraph* gr_bdtgs = new TGraph(dots,sig_eff_bdtgs,back_reg_bdtgs);
  gr_bdtgs->SetMarkerSize(1);
  gr_bdtgs->SetMarkerColor(kRed);
  gr_bdtgs->SetMarkerStyle(20);

  TGraph* gr_bdtgmbcs = new TGraph(dots,sig_eff_bdtgmbcs,back_reg_bdtgmbcs);
  gr_bdtgmbcs->SetMarkerSize(1);
  gr_bdtgmbcs->SetMarkerColor(kGreen);
  gr_bdtgmbcs->SetMarkerStyle(20);

  TGraph* gr_bdtgdes = new TGraph(dots,sig_eff_bdtgdes,back_reg_bdtgdes);
  gr_bdtgdes->SetMarkerSize(1);
  gr_bdtgdes->SetMarkerColor(kBlue);
  gr_bdtgdes->SetMarkerStyle(20);

  TGraph* gr_bdtgfr = new TGraph(dots,sig_eff_bdtgfr,back_reg_bdtgfr);
  gr_bdtgfr->SetMarkerSize(1);
  gr_bdtgfr->SetMarkerColor(kBlack);
  gr_bdtgfr->SetMarkerStyle(21);

  TGraph* gr_bdtsl = new TGraph(dots,sig_eff_bdtsl,back_reg_bdtsl);
  gr_bdtsl->SetMarkerSize(1);
  gr_bdtsl->SetMarkerColor(kBlack);
  gr_bdtsl->SetMarkerStyle(20);

  TGraph* gr_bdtgsl = new TGraph(dots,sig_eff_bdtgsl,back_reg_bdtgsl);
  gr_bdtgsl->SetMarkerSize(1);
  gr_bdtgsl->SetMarkerColor(kRed);
  gr_bdtgsl->SetMarkerStyle(22);

  TGraph* gr_ksfw = new TGraph(dots,sig_eff_ksfw,back_reg_ksfw);
  gr_ksfw->SetMarkerSize(1);
  gr_ksfw->SetMarkerColor(kBlue);
  gr_ksfw->SetMarkerStyle(22);
  
  TGraph* gr_ksfw1 = new TGraph(dots,sig_eff_ksfw1,back_reg_ksfw1);
  gr_ksfw1->SetMarkerSize(1);
  gr_ksfw1->SetMarkerColor(kYellow);
  gr_ksfw1->SetMarkerStyle(22);

  TGraph* gr_nn = new TGraph(dots,sig_eff_nn,back_reg_nn);
  gr_nn->SetMarkerSize(1);
  gr_nn->SetMarkerColor(kBlack);
  gr_nn->SetMarkerStyle(21);

  if(bdt_flag)mg->Add(gr_bdt);
  if(bdtg_flag)mg->Add(gr_bdtg);
  if(bdts_flag)mg->Add(gr_bdts);
  if(bdtgs_flag)mg->Add(gr_bdtgs);
  if(bdtgmbcs_flag)mg->Add(gr_bdtgmbcs);
  if(bdtgdes_flag)mg->Add(gr_bdtgdes);
  if(bdtgfr_flag)mg->Add(gr_bdtgfr);
  if(bdtgs_flag)mg->Add(gr_bdtgs);
  if(bdtsl_flag)mg->Add(gr_bdtsl);
  if(bdtgsl_flag)mg->Add(gr_bdtgsl);
  if(lh_flag)mg->Add(gr_ksfw);
  if(lh1_flag)mg->Add(gr_ksfw1);
  if(nn_flag)mg->Add(gr_nn);

  mg->Draw("AP");
  c1->Update();

  TCanvas* c2 = new TCanvas("c2","Significance",600,600);
  c2->cd();
  c2->SetGrid();

  TMultiGraph *signif_mg = new TMultiGraph();
  signif_mg->SetTitle("Signal Significance");
  TGraph* signif_gr_bdt = new TGraph(dots,sig_eff_bdt,signif_bdt);
  signif_gr_bdt->SetMarkerSize(1);
  signif_gr_bdt->SetMarkerColor(kBlue);
  signif_gr_bdt->SetMarkerStyle(21);

  TGraph* signif_gr_bdtg = new TGraph(dots,sig_eff_bdtg,signif_bdtg);
  signif_gr_bdtg->SetMarkerSize(1);
  signif_gr_bdtg->SetMarkerColor(kRed);
  signif_gr_bdtg->SetMarkerStyle(21);

  TGraph* signif_gr_bdts = new TGraph(dots,sig_eff_bdts,signif_bdts);
  signif_gr_bdts->SetMarkerSize(1);
  signif_gr_bdts->SetMarkerColor(kBlue);
  signif_gr_bdts->SetMarkerStyle(20);

  TGraph* signif_gr_bdtgs = new TGraph(dots,sig_eff_bdtgs,signif_bdtgs);
  signif_gr_bdtgs->SetMarkerSize(1);
  signif_gr_bdtgs->SetMarkerColor(kRed);
  signif_gr_bdtgs->SetMarkerStyle(20);
  
  TGraph* signif_gr_bdtgmbcs = new TGraph(dots,sig_eff_bdtgmbcs,signif_bdtgmbcs);
  signif_gr_bdtgmbcs->SetMarkerSize(1);
  signif_gr_bdtgmbcs->SetMarkerColor(kGreen);
  signif_gr_bdtgmbcs->SetMarkerStyle(20);
  
  TGraph* signif_gr_bdtgdes = new TGraph(dots,sig_eff_bdtgdes,signif_bdtgdes);
  signif_gr_bdtgdes->SetMarkerSize(1);
  signif_gr_bdtgdes->SetMarkerColor(kBlue);
  signif_gr_bdtgdes->SetMarkerStyle(20);
  
  TGraph* signif_gr_bdtgfr = new TGraph(dots,sig_eff_bdtgfr,signif_bdtgfr);
  signif_gr_bdtgfr->SetMarkerSize(1);
  signif_gr_bdtgfr->SetMarkerColor(kBlack);
  signif_gr_bdtgfr->SetMarkerStyle(21);

  TGraph* signif_gr_bdtsl = new TGraph(dots,sig_eff_bdtsl,signif_bdtsl);
  signif_gr_bdtsl->SetMarkerSize(1);
  signif_gr_bdtsl->SetMarkerColor(kBlack);
  signif_gr_bdtsl->SetMarkerStyle(20);

  TGraph* signif_gr_bdtgsl = new TGraph(dots,sig_eff_bdtgsl,signif_bdtgsl);
  signif_gr_bdtgsl->SetMarkerSize(1);
  signif_gr_bdtgsl->SetMarkerColor(kRed);
  signif_gr_bdtgsl->SetMarkerStyle(22);

  TGraph* signif_gr_ksfw = new TGraph(dots,sig_eff_ksfw,signif_ksfw);
  signif_gr_ksfw->SetMarkerSize(1);
  signif_gr_ksfw->SetMarkerColor(kBlue);
  signif_gr_ksfw->SetMarkerStyle(22);
  
  TGraph* signif_gr_ksfw1 = new TGraph(dots,sig_eff_ksfw1,signif_ksfw1);
  signif_gr_ksfw1->SetMarkerSize(1);
  signif_gr_ksfw1->SetMarkerColor(kYellow);
  signif_gr_ksfw1->SetMarkerStyle(22);

  TGraph* signif_gr_nn = new TGraph(dots,sig_eff_nn,signif_nn);
  signif_gr_nn->SetMarkerSize(1);
  signif_gr_nn->SetMarkerColor(kBlack);
  signif_gr_nn->SetMarkerStyle(21);

  if(bdt_flag)      signif_mg->Add(signif_gr_bdt);
  if(bdtg_flag)     signif_mg->Add(signif_gr_bdtg);
  if(bdts_flag)     signif_mg->Add(signif_gr_bdts);
  if(bdtgs_flag)    signif_mg->Add(signif_gr_bdtgs);
  if(bdtgmbcs_flag) signif_mg->Add(signif_gr_bdtgmbcs);
  if(bdtgdes_flag)  signif_mg->Add(signif_gr_bdtgdes);
  if(bdtgfr_flag)   signif_mg->Add(signif_gr_bdtgfr);
  if(bdtsl_flag)    signif_mg->Add(signif_gr_bdtsl);
  if(bdtgsl_flag)   signif_mg->Add(signif_gr_bdtgsl);
  if(lh_flag)       signif_mg->Add(signif_gr_ksfw);
  if(lh1_flag)      signif_mg->Add(signif_gr_ksfw1);
  if(nn_flag)       signif_mg->Add(signif_gr_nn);

//  signif_mg->GetXaxis()->SetRangeUser(0.50,0.95);
//  signif_mg->GetYaxis()->SetRangeUser(14.5,16.5);
  signif_mg->Draw("AP");
  c2->Update();
  
  cout << "Summary:" << endl;
  
  if(bdtgsl_flag){
    cout << "bdtgsl decision:   ";
    cout << "sig: " << bdtgsl_best_sig;
    cout << ", cut: "   << bdtgsl_best_cut;
    cout << ", pur: " << bdtgsl_best_pur;
    cout << ", eff: " << bdtgsl_best_eff;
    cout << ", sup: " << bdtgsl_best_sup << endl;
  }
  if(bdtgdes_flag){
    cout << "bdtgdes decision:  ";
    cout << "sig: " << bdtgdes_best_sig;
    cout << ", cut: "   << bdtgdes_best_cut;
    cout << ", pur: " << bdtgdes_best_pur;
    cout << ", eff: " << bdtgdes_best_eff;
    cout << ", sup: " << bdtgdes_best_sup << endl;
  }
  if(bdtgmbcs_flag){
    cout << "bdtgmbcs decision: ";
    cout << "sig: " << bdtgmbcs_best_sig;
    cout << ", cut: "   << bdtgmbcs_best_cut;
    cout << ", pur: " << bdtgmbcs_best_pur;
    cout << ", eff: " << bdtgmbcs_best_eff;
    cout << ", sup: " << bdtgmbcs_best_sup << endl;
  }
  if(bdtgfr_flag){
    cout << "bdtgdfr decision:  ";
    cout << "sig: " << bdtgfr_best_sig;
    cout << ", cut: "   << bdtgfr_best_cut;
    cout << ", pur: " << bdtgfr_best_pur;
    cout << ", eff: " << bdtgfr_best_eff;
    cout << ", sup: " << bdtgfr_best_sup << endl;
  }
  if(bdtgs_flag){
    cout << "bdtgs decision:    ";
    cout << "sig: " << bdtgs_best_sig;
    cout << ", cut: "   << bdtgs_best_cut;
    cout << ", pur: " << bdtgs_best_pur;
    cout << ", eff: " << bdtgs_best_eff;
    cout << ", sup: " << bdtgs_best_sup << endl;
  }
  if(lh_flag){
    cout << "lh decision:       ";
    cout << "sig: " << lh_best_sig;
    cout << ", cut: "   << lh_best_cut;
    cout << ", pur: " << lh_best_pur;
    cout << ", eff: " << lh_best_eff;
    cout << ", sup: " << lh_best_sup << endl;
  }
  if(lh1_flag){
    cout << "lh1 decision:      ";
    cout << "sig: " << lh1_best_sig;
    cout << ", cut: "   << lh1_best_cut;
    cout << ", pur: " << lh1_best_pur;
    cout << ", eff: " << lh1_best_eff;
    cout << ", sup: " << lh1_best_sup << endl;
  }
  return;
}
