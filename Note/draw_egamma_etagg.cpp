
void draw_egamma_etagg(void){
  TChain* sig_tree = new TChain("TEvent","TEvent");
  sig_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s3_m2_h0m10.root");
  TChain* bkg_tree = new TChain("TEvent","TEvent");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_0_10_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_1_11_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_2_12_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_3_13_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_4_14_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_5_15_m2_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_0_10_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_1_11_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_2_12_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_3_13_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_4_14_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_5_15_m2_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_0_10_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_1_11_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_2_12_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_3_13_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_4_14_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_5_15_m2_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_0_10_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_1_11_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_2_12_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_3_13_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_4_14_m2_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_5_15_m2_h0m10.root");

  const int STot = sig_tree->GetEntries();
  const int BTot = bkg_tree->GetEntries();

  // Signal histograms //
  TH1D* h_th_gl_sig = new TH1D("h_th_gl_sig","h_th_gl_sig",60,-1.,1.);
  TH1D* h_th_gh_sig = new TH1D("h_th_gh_sig","h_th_gh_sig",60,-1.,1.);

  TH1D* h_e_gl_bar_sig = new TH1D("h_e_gl_bar_sig","h_e_gl_bar_sig",100,0.08,1.8);
  TH1D* h_e_gl_end_sig = new TH1D("h_e_gl_end_sig","h_e_gl_end_sig",100,0.08,1.8);
  TH1D* h_e_gh_bar_sig = new TH1D("h_e_gh_bar_sig","h_e_gh_bar_sig",100,0.7,3.5);
  TH1D* h_e_gh_end_sig = new TH1D("h_e_gh_end_sig","h_e_gh_end_sig",100,0.7,3.5);

  TH1D* h_heleta_sig  = new TH1D("h_heleta_sig","h_heleta_sig",100,0.,1.);
  TH1D* h_peta_sig = new TH1D("h_peta_sig","h_peta_sig",100,1.4,4.);

  // Background histograms //
  TH1D* h_th_gl_bkg = new TH1D("h_th_gl_bkg","h_th_gl_bkg",60,-1.,1.);
  TH1D* h_th_gh_bkg = new TH1D("h_th_gh_bkg","h_th_gh_bkg",60,-1.,1.);

  TH1D* h_e_gl_bar_bkg = new TH1D("h_e_gl_bar_bkg","h_e_gl_bar_bkg",100,0.08,1.8);
  TH1D* h_e_gl_end_bkg = new TH1D("h_e_gl_end_bkg","h_e_gl_end_bkg",100,0.08,1.8);
  TH1D* h_e_gh_bar_bkg = new TH1D("h_e_gh_bar_bkg","h_e_gh_bar_bkg",100,0.7,3.5);
  TH1D* h_e_gh_end_bkg = new TH1D("h_e_gh_end_bkg","h_e_gh_end_bkg",100,0.7,3.5);

  TH1D* h_heleta_bkg  = new TH1D("h_heleta_sig","h_heleta_sig",100,0.,1.);
  TH1D* h_peta_bkg = new TH1D("h_peta_sig","h_peta_sig",100,1.4,4.);

  const double edge = 0.75;
  int b0f,rndm_pi0;
  double e_g1,e_g2;
  double th_g1,th_g2;
  double de,mbc;
  double hel_h0;
  double p_h0;

  sig_tree->SetBranchAddress("de",&de);
  sig_tree->SetBranchAddress("mbc",&mbc);
  sig_tree->SetBranchAddress("b0f",&b0f);
  sig_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  sig_tree->SetBranchAddress("e_g1",&e_g1);
  sig_tree->SetBranchAddress("e_g2",&e_g2);
  sig_tree->SetBranchAddress("th_g1",&th_g1);
  sig_tree->SetBranchAddress("th_g2",&th_g2);
  sig_tree->SetBranchAddress("hel_h0",&hel_h0);
  sig_tree->SetBranchAddress("p_h0",&p_h0);

  bkg_tree->SetBranchAddress("de",&de);
  bkg_tree->SetBranchAddress("mbc",&mbc);
  bkg_tree->SetBranchAddress("b0f",&b0f);
  bkg_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  bkg_tree->SetBranchAddress("e_g1",&e_g1);
  bkg_tree->SetBranchAddress("e_g2",&e_g2);
  bkg_tree->SetBranchAddress("th_g1",&th_g1);
  bkg_tree->SetBranchAddress("th_g2",&th_g2);
  bkg_tree->SetBranchAddress("hel_h0",&hel_h0);
  bkg_tree->SetBranchAddress("p_h0",&p_h0);

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

    h_heleta_sig->Fill(fabs(hel_h0));
    h_peta_sig->Fill(p_h0);

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

    h_heleta_bkg->Fill(fabs(hel_h0));
    h_peta_bkg->Fill(p_h0);

    if(fabs(thl)>edge) h_e_gl_end_bkg->Fill(el);
    else               h_e_gl_bar_bkg->Fill(el);

    if(fabs(thh)>edge) h_e_gh_end_bkg->Fill(eh);
    else               h_e_gh_bar_bkg->Fill(eh);
  }

  const bool draw_e   = true;
  const bool draw_bar = true;
  const bool draw_end = true;
  const bool draw_th  = false;
  const bool draw_mom = false;

  if(draw_mom){
  TCanvas* c_heleta_sig = new TCanvas("c_heleta_sig","c_heleta_sig",600,400);
  c_heleta_sig->cd();
  c_heleta_sig->SetGrid();
  h_heleta_sig->SetStats(0);
  h_heleta_sig->GetYaxis()->SetLabelSize(0.06);
  h_heleta_sig->GetXaxis()->SetLabelSize(0.06);
  h_heleta_sig->SetLineWidth(2);
  h_heleta_sig->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #eta (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_heleta_sig->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_heleta_sig->GetXaxis()->SetTitleSize(0.06);
  h_heleta_sig->GetXaxis()->SetTitleOffset(0.75);
  h_heleta_sig->Draw();
  c_heleta_sig->Print("pics/heleta_sig_etagg.eps");

  TCanvas* c_heleta_bkg = new TCanvas("c_heleta_bkg","c_heleta_bkg",600,400);
  c_heleta_bkg->cd();
  c_heleta_bkg->SetGrid();
  h_heleta_bkg->SetStats(0);
  h_heleta_bkg->GetYaxis()->SetLabelSize(0.06);
  h_heleta_bkg->GetXaxis()->SetLabelSize(0.06);
  h_heleta_bkg->SetLineWidth(2);
  h_heleta_bkg->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #eta (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_heleta_bkg->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_heleta_bkg->GetXaxis()->SetTitleSize(0.06);
  h_heleta_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_heleta_bkg->Draw();
  c_heleta_bkg->Print("pics/heleta_bkg_etagg.eps");

  TCanvas* c_peta_sig = new TCanvas("c_peta_sig","c_peta_sig",600,400);
  c_peta_sig->cd();
  c_peta_sig->SetGrid();
  h_peta_sig->SetStats(0);
  h_peta_sig->GetYaxis()->SetLabelSize(0.06);
  h_peta_sig->GetXaxis()->SetLabelSize(0.06);
  h_peta_sig->SetLineWidth(2);
  h_peta_sig->SetTitle("Lab frame momentum of #eta (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_peta_sig->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_peta_sig->GetXaxis()->SetTitleSize(0.06);
  h_peta_sig->GetXaxis()->SetTitleOffset(0.75);
  h_peta_sig->Draw();
  c_peta_sig->Print("pics/peta_sig_etagg.eps");

  TCanvas* c_peta_bkg = new TCanvas("c_peta_bkg","c_peta_bkg",600,400);
  c_peta_bkg->cd();
  c_peta_bkg->SetGrid();
  h_peta_bkg->SetStats(0);
  h_peta_bkg->GetYaxis()->SetLabelSize(0.06);
  h_peta_bkg->GetXaxis()->SetLabelSize(0.06);
  h_peta_bkg->SetLineWidth(2);
  h_peta_bkg->SetTitle("Lab frame momentum of #eta (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_peta_bkg->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_peta_bkg->GetXaxis()->SetTitleSize(0.06);
  h_peta_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_peta_bkg->Draw();
  c_peta_bkg->Print("pics/peta_bkg_etagg.eps");
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
  h_th_gl_sig->SetTitle("cos(#theta) of low energy #gamma from #eta (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_th_gl_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_sig->Draw();
  c_th_gl_sig->Print("pics/th_gl_sig_etagg.eps");

  TCanvas* c_th_gh_sig = new TCanvas("c_th_gh_sig","c_th_gh_sig",600,400);
  c_th_gh_sig->cd();
  c_th_gh_sig->cd();
  c_th_gh_sig->SetGrid();
  h_th_gh_sig->SetStats(0);
  h_th_gh_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_sig->SetLineWidth(2);
  h_th_gh_sig->SetTitle("cos(#theta) of high energy #gamma from #eta (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_th_gh_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_sig->Draw();
  c_th_gh_sig->Print("pics/th_gh_sig_etagg.eps");
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
  h_e_gl_bar_sig->SetTitle("Lower E_{#gamma} from #eta for barrel (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gl_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_sig->Draw();
  c_e_gl_bar_sig->Print("pics/e_gl_bar_sig_etagg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_sig  = new TCanvas("c_e_gl_end_sig ","c_e_gl_end_sig ",600,400);
  c_e_gl_end_sig->cd();
  c_e_gl_end_sig->SetGrid();
  h_e_gl_end_sig->SetStats(0);
  h_e_gl_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->SetLineWidth(2);
  h_e_gl_end_sig->SetTitle("Lower E_{#gamma} from #eta for endcap (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gl_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_sig->Draw();
  c_e_gl_end_sig->Print("pics/e_gl_end_sig_etagg.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_sig  = new TCanvas("c_e_gh_bar_sig ","c_e_gh_bar_sig ",600,400);
  c_e_gh_bar_sig->cd();
  c_e_gh_bar_sig->SetGrid();
  h_e_gh_bar_sig->SetStats(0);
  h_e_gh_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->SetLineWidth(2);
  h_e_gh_bar_sig->SetTitle("Higher E_{#gamma} from #eta for barrel (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gh_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_sig->Draw();
  c_e_gh_bar_sig->Print("pics/e_gh_bar_sig_etagg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_sig  = new TCanvas("c_e_gh_end_sig ","c_e_gh_end_sig ",600,400);
  c_e_gh_end_sig->cd();
  c_e_gh_end_sig->SetGrid();
  h_e_gh_end_sig->SetStats(0);
  h_e_gh_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->SetLineWidth(2);
  h_e_gh_end_sig->SetTitle("Higher E_{#gamma} from #eta for endcap (sig D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gh_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_sig->Draw();
  c_e_gh_end_sig->Print("pics/e_gh_end_sig_etagg.eps");
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
  h_th_gl_bkg->SetTitle("cos(#theta) of low energy #gamma from #eta (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_th_gl_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_bkg->Draw();
  c_th_gl_bkg->Print("pics/th_gl_bkg_etagg.eps");

  TCanvas* c_th_gh_bkg = new TCanvas("c_th_gh_bkg","c_th_gh_bkg",600,400);
  c_th_gh_bkg->cd();
  c_th_gh_bkg->SetGrid();
  h_th_gh_bkg->SetStats(0);
  h_th_gh_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->SetLineWidth(2);
  h_th_gh_bkg->SetTitle("cos(#theta) of high energy #gamma from #eta (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_th_gh_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_bkg->Draw();
  c_th_gh_bkg->Print("pics/th_gh_bkg_etagg.eps");
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
  h_e_gl_bar_bkg->SetTitle("Lower E_{#gamma} from #eta for barrel (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gl_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_bkg->Draw();
  c_e_gl_bar_bkg->Print("pics/e_gl_bar_bkg_etagg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_bkg  = new TCanvas("c_e_gl_end_bkg ","c_e_gl_end_bkg ",600,400);
  c_e_gl_end_bkg->cd();
  c_e_gl_end_bkg->SetGrid();
  h_e_gl_end_bkg->SetStats(0);
  h_e_gl_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->SetLineWidth(2);
  h_e_gl_end_bkg->SetTitle("Lower E_{#gamma} from #eta for endcap (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gl_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_bkg->Draw();
  c_e_gl_end_bkg->Print("pics/e_gl_end_bkg_etagg.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_bkg  = new TCanvas("c_e_gh_bar_bkg ","c_e_gh_bar_bkg ",600,400);
  c_e_gh_bar_bkg->cd();
  c_e_gh_bar_bkg->SetGrid();
  h_e_gh_bar_bkg->SetStats(0);
  h_e_gh_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->SetLineWidth(2);
  h_e_gh_bar_bkg->SetTitle("Higher E_{#gamma} from #eta for barrel (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gh_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_bkg->Draw();
  c_e_gh_bar_bkg->Print("pics/e_gh_bar_bkg_etagg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_bkg  = new TCanvas("c_e_gh_end_bkg ","c_e_gh_end_bkg ",600,400);
  c_e_gh_end_bkg->cd();
  c_e_gh_end_bkg->SetGrid();
  h_e_gh_end_bkg->SetStats(0);
  h_e_gh_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->SetLineWidth(2);
  h_e_gh_end_bkg->SetTitle("Higher E_{#gamma} from #eta for endcap (bkg D^{0}#eta(#rightarrow#gamma#gamma))");
  h_e_gh_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_bkg->Draw();
  c_e_gh_end_bkg->Print("pics/e_gh_end_bkg_etagg.eps");
  }
  }

  return;
}
