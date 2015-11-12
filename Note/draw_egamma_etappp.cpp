
void draw_egamma_etappp(void){
  TChain* sig_tree = new TChain("TEvent","TEvent");
  sig_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s3_m2_h0m20.root");
  TChain* bkg_tree = new TChain("TEvent","TEvent");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_0_10_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_1_11_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_2_12_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_3_13_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_4_14_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_5_15_m2_h0m20.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_0_10_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_1_11_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_2_12_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_3_13_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_4_14_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_5_15_m2_h0m20.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_0_10_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_1_11_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_2_12_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_3_13_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_4_14_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_5_15_m2_h0m20.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_0_10_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_1_11_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_2_12_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_3_13_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_4_14_m2_h0m20.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_5_15_m2_h0m20.root");

  const int STot = sig_tree->GetEntries();
  const int BTot = bkg_tree->GetEntries();

  // Signal histograms //
  TH1D* h_th_gl_sig = new TH1D("h_th_gl_sig","h_th_gl_sig",60,-1.,1.);
  TH1D* h_th_gh_sig = new TH1D("h_th_gh_sig","h_th_gh_sig",60,-1.,1.);

  TH1D* h_e_gl_bar_sig = new TH1D("h_e_gl_bar_sig","h_e_gl_bar_sig",100,0.04,1);
  TH1D* h_e_gl_end_sig = new TH1D("h_e_gl_end_sig","h_e_gl_end_sig",100,0.04,1);
  TH1D* h_e_gh_bar_sig = new TH1D("h_e_gh_bar_sig","h_e_gh_bar_sig",100,0.04,2);
  TH1D* h_e_gh_end_sig = new TH1D("h_e_gh_end_sig","h_e_gh_end_sig",100,0.04,2);

  TH1D* h_ppi0_sig = new TH1D("h_ppi0_sig","h_ppi0_sig",100,0.,3.);
  TH1D* h_ppi1_sig = new TH1D("h_ppi1_sig","h_ppi1_sig",100,0.,3.);
  TH1D* h_ppi2_sig = new TH1D("h_ppi2_sig","h_ppi2_sig",100,0.,3.);
  TH1D* h_helpi0_sig  = new TH1D("h_helpi0_sig","h_helpi0_sig",60,0.,1.);

  TH1D* h_ptpi1_sig = new TH1D("h_ptpi1_sig","h_ptpi1_sig",100,0.,2.);
  TH1D* h_ptpi2_sig = new TH1D("h_ptpi2_sig","h_ptpi2_sig",100,0.,2.);

  // Background histograms //
  TH1D* h_th_gl_bkg = new TH1D("h_th_gl_bkg","h_th_gl_bkg",60,-1.,1.);
  TH1D* h_th_gh_bkg = new TH1D("h_th_gh_bkg","h_th_gh_bkg",60,-1.,1.);

  TH1D* h_e_gl_bar_bkg = new TH1D("h_e_gl_bar_bkg","h_e_gl_bar_bkg",100,0.04,1);
  TH1D* h_e_gl_end_bkg = new TH1D("h_e_gl_end_bkg","h_e_gl_end_bkg",100,0.04,1);
  TH1D* h_e_gh_bar_bkg = new TH1D("h_e_gh_bar_bkg","h_e_gh_bar_bkg",100,0.04,2);
  TH1D* h_e_gh_end_bkg = new TH1D("h_e_gh_end_bkg","h_e_gh_end_bkg",100,0.04,2);

  TH1D* h_ppi0_bkg = new TH1D("h_ppi0_bkg","h_ppi0_bkg",100,0.,3.);
  TH1D* h_ppi1_bkg = new TH1D("h_ppi1_bkg","h_ppi1_bkg",100,0.,3.);
  TH1D* h_ppi2_bkg = new TH1D("h_ppi2_bkg","h_ppi2_bkg",100,0.,3.);
  TH1D* h_helpi0_bkg  = new TH1D("h_helpi0_bkg","h_helpi0_bkg",60,0.,1.);

  TH1D* h_ptpi1_bkg = new TH1D("h_ptpi1_bkg","h_ptpi1_bkg",100,0.,2.);
  TH1D* h_ptpi2_bkg = new TH1D("h_ptpi2_bkg","h_ptpi2_bkg",100,0.,2.);

  const double edge = 0.7;
  int b0f,rndm_pi0;
  double e_g1,e_g2;
  double th_g1,th_g2;
  double de,mbc;
  double p_pi0,p_pi1,p_pi2;
  double hel_pi0;
  double bdt;
  double pt_pi1;
  double pt_pi2;

  sig_tree->SetBranchAddress("de",&de);
  sig_tree->SetBranchAddress("mbc",&mbc);
  sig_tree->SetBranchAddress("b0f",&b0f);
  sig_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  sig_tree->SetBranchAddress("e_g1",&e_g1);
  sig_tree->SetBranchAddress("e_g2",&e_g2);
  sig_tree->SetBranchAddress("th_g1",&th_g1);
  sig_tree->SetBranchAddress("th_g2",&th_g2);
  sig_tree->SetBranchAddress("p_pi0_h0",&p_pi0);
  sig_tree->SetBranchAddress("p_pip_h0",&p_pi1);
  sig_tree->SetBranchAddress("p_pim_h0",&p_pi2);
  sig_tree->SetBranchAddress("hel_pi0",&hel_pi0);
  sig_tree->SetBranchAddress("bdt",&bdt);
  sig_tree->SetBranchAddress("pt_pi1",&pt_pi1);
  sig_tree->SetBranchAddress("pt_pi2",&pt_pi2);

  bkg_tree->SetBranchAddress("de",&de);
  bkg_tree->SetBranchAddress("mbc",&mbc);
  bkg_tree->SetBranchAddress("b0f",&b0f);
  bkg_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  bkg_tree->SetBranchAddress("e_g1",&e_g1);
  bkg_tree->SetBranchAddress("e_g2",&e_g2);
  bkg_tree->SetBranchAddress("th_g1",&th_g1);
  bkg_tree->SetBranchAddress("th_g2",&th_g2);
  bkg_tree->SetBranchAddress("p_pi0_h0",&p_pi0);
  bkg_tree->SetBranchAddress("p_pip_h0",&p_pi1);
  bkg_tree->SetBranchAddress("p_pim_h0",&p_pi2);
  bkg_tree->SetBranchAddress("hel_pi0",&hel_pi0);
  bkg_tree->SetBranchAddress("bdt",&bdt);
  bkg_tree->SetBranchAddress("pt_pi1",&pt_pi1);
  bkg_tree->SetBranchAddress("pt_pi2",&pt_pi2);

  double el,eh;
  double thl,thh;
  cout << "Signal" << endl;
  for(int i=0; i<STot; i++){
    sig_tree->GetEvent(i);
    if(abs(de)>0.1 || mbc<5.271) continue;
    if(b0f != 1 && b0f != 5 && b0f != 10 && rndm_pi0 != 1) continue;
    thl = e_g1 > e_g2 ? th_g2 : th_g1;
    thh = e_g1 > e_g2 ? th_g1 : th_g2;
    el  = e_g1 > e_g2 ? e_g2  : e_g1;
    eh  = e_g1 > e_g2 ? e_g1  : e_g2;

    h_th_gl_sig->Fill(thl);
    h_th_gh_sig->Fill(thh);

    h_ppi0_sig->Fill(p_pi0);
    h_ppi1_sig->Fill(p_pi1);
    h_ppi2_sig->Fill(p_pi2);
    h_helpi0_sig->Fill(hel_pi0);

    h_ptpi1_sig->Fill(pt_pi1);
    h_ptpi2_sig->Fill(pt_pi2);

    if(!(i%10000)) cout << i << ": " << el << " " << eh << endl;

    if(fabs(thl)>edge) h_e_gl_end_sig->Fill(el);
    else               h_e_gl_bar_sig->Fill(el);

    if(fabs(thh)>edge) h_e_gh_end_sig->Fill(eh);
    else               h_e_gh_bar_sig->Fill(eh);
  }

  cout << "Background" << endl;
  for(int i=0; i<BTot; i++){
    bkg_tree->GetEvent(i);
    if(!(i%10000)) cout << i << ": " << el << " " << eh << endl;
    if(abs(de)>0.1 || mbc<5.271) continue;
    if(b0f == 1 || b0f == 5 || b0f == 10 || rndm_pi0 == 1) continue;
    thl = e_g1 > e_g2 ? th_g2 : th_g1;
    thh = e_g1 > e_g2 ? th_g1 : th_g2;
    el  = e_g1 > e_g2 ? e_g2  : e_g1;
    eh  = e_g1 > e_g2 ? e_g1  : e_g2;

    h_th_gl_bkg->Fill(thl);
    h_th_gh_bkg->Fill(thh);

    h_ppi0_bkg->Fill(p_pi0);
    h_ppi1_bkg->Fill(p_pi1);
    h_ppi2_bkg->Fill(p_pi2);
    h_helpi0_bkg->Fill(hel_pi0);

    h_ptpi1_bkg->Fill(pt_pi1);
    h_ptpi2_bkg->Fill(pt_pi2);

    if(fabs(thl)>edge) h_e_gl_end_bkg->Fill(el);
    else               h_e_gl_bar_bkg->Fill(el);

    if(fabs(thh)>edge) h_e_gh_end_bkg->Fill(eh);
    else               h_e_gh_bar_bkg->Fill(eh);
  }

  const bool draw_e   = false;
  const bool draw_bar = false;
  const bool draw_end = false;
  const bool draw_th  = false;
  const bool draw_mom = true;

  if(draw_mom){
  TCanvas* c_ptpi1_sig = new TCanvas("c_ptpi1_sig","c_ptpi1_sig",600,400);
  c_ptpi1_sig->cd();
  c_ptpi1_sig->SetGrid();
  h_ptpi1_sig->SetStats(0);
  h_ptpi1_sig->GetYaxis()->SetLabelSize(0.06);
  h_ptpi1_sig->GetXaxis()->SetLabelSize(0.06);
  h_ptpi1_sig->SetLineWidth(2);
  h_ptpi1_sig->SetTitle("p_{t} for #pi^{+} from #eta decay (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ptpi1_sig->GetXaxis()->SetTitle("p_{t} (GeV)");
  h_ptpi1_sig->GetXaxis()->SetTitleSize(0.06);
  h_ptpi1_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ptpi1_sig->Draw();
  c_ptpi1_sig->Print("pics/ptpi1_sig_etappp.eps");

  TCanvas* c_ptpi1_bkg = new TCanvas("c_ptpi1_bkg","c_ptpi1_bkg",600,400);
  c_ptpi1_bkg->cd();
  c_ptpi1_bkg->SetGrid();
  h_ptpi1_bkg->SetStats(0);
  h_ptpi1_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ptpi1_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ptpi1_bkg->SetLineWidth(2);
  h_ptpi1_bkg->SetTitle("p_{t} for #pi^{+} from #eta decay (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ptpi1_bkg->GetXaxis()->SetTitle("p_{t} (GeV)");
  h_ptpi1_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ptpi1_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ptpi1_bkg->Draw();
  c_ptpi1_bkg->Print("pics/ptpi1_bkg_etappp.eps");

  TCanvas* c_ptpi2_sig = new TCanvas("c_ptpi2_sig","c_ptpi2_sig",600,400);
  c_ptpi2_sig->cd();
  c_ptpi2_sig->SetGrid();
  h_ptpi2_sig->SetStats(0);
  h_ptpi2_sig->GetYaxis()->SetLabelSize(0.06);
  h_ptpi2_sig->GetXaxis()->SetLabelSize(0.06);
  h_ptpi2_sig->SetLineWidth(2);
  h_ptpi2_sig->SetTitle("p_{t} for #pi^{-} from #eta decay (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ptpi2_sig->GetXaxis()->SetTitle("p_{t} (GeV)");
  h_ptpi2_sig->GetXaxis()->SetTitleSize(0.06);
  h_ptpi2_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ptpi2_sig->Draw();
  c_ptpi2_sig->Print("pics/ptpi2_sig_etappp.eps");

  TCanvas* c_ptpi2_bkg = new TCanvas("c_ptpi2_bkg","c_ptpi2_bkg",600,400);
  c_ptpi2_bkg->cd();
  c_ptpi2_bkg->SetGrid();
  h_ptpi2_bkg->SetStats(0);
  h_ptpi2_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ptpi2_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ptpi2_bkg->SetLineWidth(2);
  h_ptpi2_bkg->SetTitle("p_{t} for #pi^{-} from #eta decay (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ptpi2_bkg->GetXaxis()->SetTitle("p_{t} (GeV)");
  h_ptpi2_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ptpi2_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ptpi2_bkg->Draw();
  c_ptpi2_bkg->Print("pics/ptpi2_bkg_etappp.eps");

  TCanvas* c_helpi0_sig = new TCanvas("c_helpi0_sig","c_helpi0_sig",600,400);
  c_helpi0_sig->cd();
  c_helpi0_sig->SetGrid();
  h_helpi0_sig->SetStats(0);
  h_helpi0_sig->GetYaxis()->SetLabelSize(0.06);
  h_helpi0_sig->GetXaxis()->SetLabelSize(0.06);
  h_helpi0_sig->SetLineWidth(2);
  h_helpi0_sig->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #pi^{0} from #eta decay (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_helpi0_sig->GetXaxis()->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L})");
  h_helpi0_sig->GetXaxis()->SetTitleSize(0.06);
  h_helpi0_sig->GetXaxis()->SetTitleOffset(0.75);
  h_helpi0_sig->Draw();
  c_helpi0_sig->Print("pics/helpi0_sig_etappp.eps");

  TCanvas* c_helpi0_bkg = new TCanvas("c_helpi0_bkg","c_helpi0_bkg",600,400);
  c_helpi0_bkg->cd();
  c_helpi0_bkg->SetGrid();
  h_helpi0_bkg->SetStats(0);
  h_helpi0_bkg->GetYaxis()->SetLabelSize(0.06);
  h_helpi0_bkg->GetXaxis()->SetLabelSize(0.06);
  h_helpi0_bkg->SetLineWidth(2);
  h_helpi0_bkg->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #pi^{0} from #eta decay (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_helpi0_bkg->GetXaxis()->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L})");
  h_helpi0_bkg->GetXaxis()->SetTitleSize(0.06);
  h_helpi0_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_helpi0_bkg->Draw();
  c_helpi0_bkg->Print("pics/helpi0_bkg_etappp.eps");

  TCanvas* c_ppi0_sig = new TCanvas("c_ppi0_sig","c_ppi0_sig",600,400);
  c_ppi0_sig->cd();
  c_ppi0_sig->SetGrid();
  h_ppi0_sig->SetStats(0);
  h_ppi0_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppi0_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppi0_sig->SetLineWidth(2);
  h_ppi0_sig->SetTitle("#pi^{0} momentum (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi0_sig->GetXaxis()->SetTitle("p(#pi^{0}) (GeV))");
  h_ppi0_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppi0_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppi0_sig->GetXaxis()->SetRangeUser(0.,1.);
  h_ppi0_sig->Draw();
  c_ppi0_sig->Print("pics/ppi0_sig_etappp.eps");

  TCanvas* c_ppi0_bkg = new TCanvas("c_ppi0_bkg","c_ppi0_bkg",600,400);
  c_ppi0_bkg->cd();
  c_ppi0_bkg->SetGrid();
  h_ppi0_bkg->SetStats(0);
  h_ppi0_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppi0_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppi0_bkg->SetLineWidth(2);
  h_ppi0_bkg->SetTitle("#pi^{0} momentum (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi0_bkg->GetXaxis()->SetTitle("p(#pi^{0}) (GeV))");
  h_ppi0_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppi0_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppi0_bkg->GetXaxis()->SetRangeUser(0.,1.);
  h_ppi0_bkg->Draw();
  c_ppi0_bkg->Print("pics/ppi0_bkg_etappp.eps");

  TCanvas* c_ppi1_sig = new TCanvas("c_ppi1_sig","c_ppi1_sig",600,400);
  c_ppi1_sig->cd();
  c_ppi1_sig->SetGrid();
  h_ppi1_sig->SetStats(0);
  h_ppi1_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppi1_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppi1_sig->SetLineWidth(2);
  h_ppi1_sig->SetTitle("#pi^{+} from #eta momentum (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi1_sig->GetXaxis()->SetTitle("p (GeV))");
  h_ppi1_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppi1_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppi1_sig->Draw();
  c_ppi1_sig->Print("pics/ppi1_sig_etappp.eps");

  TCanvas* c_ppi1_bkg = new TCanvas("c_ppi1_bkg","c_ppi1_bkg",600,400);
  c_ppi1_bkg->cd();
  c_ppi1_bkg->SetGrid();
  h_ppi1_bkg->SetStats(0);
  h_ppi1_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppi1_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppi1_bkg->SetLineWidth(2);
  h_ppi1_bkg->SetTitle("#pi^{+} from #eta momentum (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi1_bkg->GetXaxis()->SetTitle("p (GeV))");
  h_ppi1_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppi1_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppi1_bkg->Draw();
  c_ppi1_bkg->Print("pics/ppi1_bkg_etappp.eps");

  TCanvas* c_ppi2_sig = new TCanvas("c_ppi2_sig","c_ppi2_sig",600,400);
  c_ppi2_sig->cd();
  c_ppi2_sig->SetGrid();
  h_ppi2_sig->SetStats(0);
  h_ppi2_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppi2_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppi2_sig->SetLineWidth(2);
  h_ppi2_sig->SetTitle("#pi^{-} from #eta momentum (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi2_sig->GetXaxis()->SetTitle("p (GeV))");
  h_ppi2_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppi2_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppi2_sig->Draw();
  c_ppi2_sig->Print("pics/ppi2_sig_etappp.eps");

  TCanvas* c_ppi2_bkg = new TCanvas("c_ppi2_bkg","c_ppi2_bkg",600,400);
  c_ppi2_bkg->cd();
  c_ppi2_bkg->SetGrid();
  h_ppi2_bkg->SetStats(0);
  h_ppi2_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppi2_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppi2_bkg->SetLineWidth(2);
  h_ppi2_bkg->SetTitle("#pi^{-} from #eta momentum (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_ppi2_bkg->GetXaxis()->SetTitle("p (GeV))");
  h_ppi2_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppi2_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppi2_bkg->Draw();
  c_ppi2_bkg->Print("pics/ppi2_bkg_etappp.eps");
  }

  if(draw_th){
  TCanvas* c_th_gl_sig = new TCanvas("c_th_gl_sig","c_th_gl_sig",600,400);
  c_th_gl_sig->cd();
  c_th_gl_sig->cd();
  c_th_gl_sig->SetGrid();
  h_th_gl_sig->SetStats(0);
  h_th_gl_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_sig->SetLineWidth(2);
  h_th_gl_sig->SetTitle("cos(#theta) of low energy #gamma from #pi^{0} (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_th_gl_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_sig->Draw();
  c_th_gl_sig->Print("pics/th_gl_sig_etappp.eps");

  TCanvas* c_th_gh_sig = new TCanvas("c_th_gh_sig","c_th_gh_sig",600,400);
  c_th_gh_sig->cd();
  c_th_gh_sig->cd();
  c_th_gh_sig->SetGrid();
  h_th_gh_sig->SetStats(0);
  h_th_gh_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_sig->SetLineWidth(2);
  h_th_gh_sig->SetTitle("cos(#theta) of high energy #gamma from #pi^{0} (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_th_gh_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_sig->Draw();
  c_th_gh_sig->Print("pics/th_gh_sig_etappp.eps");
  }

  if(draw_e){
  if(draw_bar){
  TCanvas* c_e_gl_bar_sig  = new TCanvas("c_e_gl_bar_sig ","c_e_gl_bar_sig ",600,400);
  c_e_gl_bar_sig->cd();
  c_e_gl_bar_sig->SetGrid();
  h_e_gl_bar_sig->SetStats(0);
  h_e_gl_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_bar_sig->SetLineWidth(2);
  h_e_gl_bar_sig->SetTitle("Lower E_{#gamma} from #pi^{0} for barrel (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gl_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_sig->Draw();
  c_e_gl_bar_sig->Print("pics/e_gl_bar_sig_etappp.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_sig  = new TCanvas("c_e_gl_end_sig ","c_e_gl_end_sig ",600,400);
  c_e_gl_end_sig->cd();
  c_e_gl_end_sig->SetGrid();
  h_e_gl_end_sig->SetStats(0);
  h_e_gl_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->SetLineWidth(2);
  h_e_gl_end_sig->SetTitle("Lower E_{#gamma} from #pi^{0} for endcap (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gl_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_sig->Draw();
  c_e_gl_end_sig->Print("pics/e_gl_end_sig_etappp.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_sig  = new TCanvas("c_e_gh_bar_sig ","c_e_gh_bar_sig ",600,400);
  c_e_gh_bar_sig->cd();
  c_e_gh_bar_sig->SetGrid();
  h_e_gh_bar_sig->SetStats(0);
  h_e_gh_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->SetLineWidth(2);
  h_e_gh_bar_sig->SetTitle("Higher E_{#gamma} from #pi^{0} for barrel (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gh_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_sig->Draw();
  c_e_gh_bar_sig->Print("pics/e_gh_bar_sig_etappp.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_sig  = new TCanvas("c_e_gh_end_sig ","c_e_gh_end_sig ",600,400);
  c_e_gh_end_sig->cd();
  c_e_gh_end_sig->SetGrid();
  h_e_gh_end_sig->SetStats(0);
  h_e_gh_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->SetLineWidth(2);
  h_e_gh_end_sig->SetTitle("Higher E_{#gamma} from #pi^{0} for endcap (sig D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gh_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_sig->Draw();
  c_e_gh_end_sig->Print("pics/e_gh_end_sig_etappp.eps");
  }
  }

  if(draw_th){
  TCanvas* c_th_gl_bkg = new TCanvas("c_th_gl_bkg","c_th_gl_bkg",600,400);
  c_th_gl_bkg->cd();
  c_th_gl_bkg->SetGrid();
  h_th_gl_bkg->SetStats(0);
  h_th_gl_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_bkg->SetLineWidth(2);
  h_th_gl_bkg->SetTitle("cos(#theta) of low energy #gamma from #pi^{0} (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_th_gl_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_bkg->Draw();
  c_th_gl_bkg->Print("pics/th_gl_bkg_etappp.eps");

  TCanvas* c_th_gh_bkg = new TCanvas("c_th_gh_bkg","c_th_gh_bkg",600,400);
  c_th_gh_bkg->cd();
  c_th_gh_bkg->SetGrid();
  h_th_gh_bkg->SetStats(0);
  h_th_gh_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->SetLineWidth(2);
  h_th_gh_bkg->SetTitle("cos(#theta) of high energy #gamma from #pi^{0} (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_th_gh_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_bkg->Draw();
  c_th_gh_bkg->Print("pics/th_gh_bkg_etappp.eps");
  }

  if(draw_e){
  if(draw_bar){
  TCanvas* c_e_gl_bar_bkg  = new TCanvas("c_e_gl_bar_bkg ","c_e_gl_bar_bkg ",600,400);
  c_e_gl_bar_bkg->cd();
  c_e_gl_bar_bkg->SetGrid();
  h_e_gl_bar_bkg->SetStats(0);
  h_e_gl_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_bar_bkg->SetLineWidth(2);
  h_e_gl_bar_bkg->SetTitle("Lower E_{#gamma} from #pi^{0} for barrel (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gl_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_bkg->Draw();
  c_e_gl_bar_bkg->Print("pics/e_gl_bar_bkg_etappp.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_bkg  = new TCanvas("c_e_gl_end_bkg ","c_e_gl_end_bkg ",600,400);
  c_e_gl_end_bkg->cd();
  c_e_gl_end_bkg->SetGrid();
  h_e_gl_end_bkg->SetStats(0);
  h_e_gl_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->SetLineWidth(2);
  h_e_gl_end_bkg->SetTitle("Lower E_{#gamma} from #pi^{0} for endcap (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gl_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_bkg->Draw();
  c_e_gl_end_bkg->Print("pics/e_gl_end_bkg_etappp.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_bkg  = new TCanvas("c_e_gh_bar_bkg ","c_e_gh_bar_bkg ",600,400);
  c_e_gh_bar_bkg->cd();
  c_e_gh_bar_bkg->SetGrid();
  h_e_gh_bar_bkg->SetStats(0);
  h_e_gh_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->SetLineWidth(2);
  h_e_gh_bar_bkg->SetTitle("Higher E_{#gamma} from #pi^{0} for barrel (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gh_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_bkg->Draw();
  c_e_gh_bar_bkg->Print("pics/e_gh_bar_bkg_etappp.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_bkg  = new TCanvas("c_e_gh_end_bkg ","c_e_gh_end_bkg ",600,400);
  c_e_gh_end_bkg->cd();
  c_e_gh_end_bkg->SetGrid();
  h_e_gh_end_bkg->SetStats(0);
  h_e_gh_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->SetLineWidth(2);
  h_e_gh_end_bkg->SetTitle("Higher E_{#gamma} from #pi^{0} for endcap (bkg D^{0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0}))");
  h_e_gh_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_bkg->Draw();
  c_e_gh_end_bkg->Print("pics/e_gh_end_bkg_etappp.eps");
  }
  }

  return;
}
