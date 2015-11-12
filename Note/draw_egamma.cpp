
void draw_egamma(void){
  TChain* sig_tree = new TChain("TEvent","TEvent");
  sig_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1_m10_h0m10.root");
  TChain* bkg_tree = new TChain("TEvent","TEvent");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_0_10_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_1_11_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_2_12_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_3_13_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_4_14_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_5_15_m10_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_0_10_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_1_11_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_2_12_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_3_13_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_4_14_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_5_15_m10_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_0_10_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_1_11_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_2_12_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_3_13_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_4_14_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_5_15_m10_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_0_10_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_1_11_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_2_12_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_3_13_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_4_14_m10_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_5_15_m10_h0m10.root");

  const int STot = sig_tree->GetEntries();
  const int BTot = bkg_tree->GetEntries();

  // Signal histograms //
  TH1D* h_th_gl_hard_sig = new TH1D("h_th_gl_hard_sig","h_th_gl_hard_sig",60,-1.,1.);
  TH1D* h_th_gh_hard_sig = new TH1D("h_th_gh_hard_sig","h_th_gh_hard_sig",60,-1.,1.);

  TH1D* h_th_gl_soft_sig = new TH1D("h_th_gl_soft_sig","h_th_gl_soft_sig",60,-1.,1.);
  TH1D* h_th_gh_soft_sig = new TH1D("h_th_gh_soft_sig","h_th_gh_soft_sig",60,-1.,1.);

  TH1D* h_e_gl_hard_bar_sig = new TH1D("h_e_gl_hard_bar_sig","h_e_gl_hard_bar_sig",60,0.04,2.0);
  TH1D* h_e_gl_hard_end_sig = new TH1D("h_e_gl_hard_end_sig","h_e_gl_hard_end_sig",60,0.04,2.0);
  TH1D* h_e_gh_hard_bar_sig = new TH1D("h_e_gh_hard_bar_sig","h_e_gh_hard_bar_sig",60,0.04,3.5);
  TH1D* h_e_gh_hard_end_sig = new TH1D("h_e_gh_hard_end_sig","h_e_gh_hard_end_sig",60,0.04,3.5);

  TH1D* h_e_gl_soft_bar_sig = new TH1D("h_e_gl_soft_bar_sig","h_e_gl_soft_bar_sig",60,0.04,0.20);
  TH1D* h_e_gl_soft_end_sig = new TH1D("h_e_gl_soft_end_sig","h_e_gl_soft_end_sig",60,0.04,0.20);
  TH1D* h_e_gh_soft_bar_sig = new TH1D("h_e_gh_soft_bar_sig","h_e_gh_soft_bar_sig",60,0.04,0.35);
  TH1D* h_e_gh_soft_end_sig = new TH1D("h_e_gh_soft_end_sig","h_e_gh_soft_end_sig",60,0.04,0.35);

  // Background histograms //
  TH1D* h_th_gl_hard_bkg = new TH1D("h_th_gl_hard_bkg","h_th_gl_hard_bkg",60,-1.,1.);
  TH1D* h_th_gh_hard_bkg = new TH1D("h_th_gh_hard_bkg","h_th_gh_hard_bkg",60,-1.,1.);

  TH1D* h_th_gl_soft_bkg = new TH1D("h_th_gl_soft_bkg","h_th_gl_soft_bkg",60,-1.,1.);
  TH1D* h_th_gh_soft_bkg = new TH1D("h_th_gh_soft_bkg","h_th_gh_soft_bkg",60,-1.,1.);

  TH1D* h_e_gl_hard_bar_bkg = new TH1D("h_e_gl_hard_bar_bkg","h_e_gl_hard_bar_bkg",60,0.04,2.0);
  TH1D* h_e_gl_hard_end_bkg = new TH1D("h_e_gl_hard_end_bkg","h_e_gl_hard_end_bkg",60,0.04,2.0);
  TH1D* h_e_gh_hard_bar_bkg = new TH1D("h_e_gh_hard_bar_bkg","h_e_gh_hard_bar_bkg",60,0.04,3.5);
  TH1D* h_e_gh_hard_end_bkg = new TH1D("h_e_gh_hard_end_bkg","h_e_gh_hard_end_bkg",60,0.04,3.5);

  TH1D* h_e_gl_soft_bar_bkg = new TH1D("h_e_gl_soft_bar_bkg","h_e_gl_soft_bar_bkg",60,0.04,0.20);
  TH1D* h_e_gl_soft_end_bkg = new TH1D("h_e_gl_soft_end_bkg","h_e_gl_soft_end_bkg",60,0.04,0.20);
  TH1D* h_e_gh_soft_bar_bkg = new TH1D("h_e_gh_soft_bar_bkg","h_e_gh_soft_bar_bkg",60,0.04,0.35);
  TH1D* h_e_gh_soft_end_bkg = new TH1D("h_e_gh_soft_end_bkg","h_e_gh_soft_end_bkg",60,0.04,0.35);

  const double edge = 0.70;
  int b0f,rndm_pi0;
  double e_g1,e_g2,e_g3,e_g4;
  double th_g1,th_g2,th_g3,th_g4;
  double de,mbc;

  sig_tree->SetBranchAddress("de",&de);
  sig_tree->SetBranchAddress("mbc",&mbc);
  sig_tree->SetBranchAddress("b0f",&b0f);
  sig_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  sig_tree->SetBranchAddress("e_g1",&e_g1);
  sig_tree->SetBranchAddress("e_g2",&e_g2);
  sig_tree->SetBranchAddress("e_g3",&e_g3);
  sig_tree->SetBranchAddress("e_g4",&e_g4);
  sig_tree->SetBranchAddress("th_g1",&th_g1);
  sig_tree->SetBranchAddress("th_g2",&th_g2);
  sig_tree->SetBranchAddress("th_g3",&th_g3);
  sig_tree->SetBranchAddress("th_g4",&th_g4);

  bkg_tree->SetBranchAddress("de",&de);
  bkg_tree->SetBranchAddress("mbc",&mbc);
  bkg_tree->SetBranchAddress("b0f",&b0f);
  bkg_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  bkg_tree->SetBranchAddress("e_g1",&e_g1);
  bkg_tree->SetBranchAddress("e_g2",&e_g2);
  bkg_tree->SetBranchAddress("e_g3",&e_g3);
  bkg_tree->SetBranchAddress("e_g4",&e_g4);
  bkg_tree->SetBranchAddress("th_g1",&th_g1);
  bkg_tree->SetBranchAddress("th_g2",&th_g2);
  bkg_tree->SetBranchAddress("th_g3",&th_g3);
  bkg_tree->SetBranchAddress("th_g4",&th_g4);

  double el_hard,eh_hard;
  double el_soft,eh_soft;
  double thl_hard,thh_hard;
  double thl_soft,thh_soft;
  cout << "Signal" << endl;
  for(int i=0; i<STot; i++){
    sig_tree->GetEvent(i);
    if(abs(de)>0.1 || mbc<5.271) continue;
    if(b0f != 1 && b0f != 5 && b0f != 10 && rndm_pi0 != 1) continue;
    thl_hard = e_g1 > e_g2 ? th_g2 : th_g1;
    thh_hard = e_g1 > e_g2 ? th_g1 : th_g2;
    el_hard  = e_g1 > e_g2 ? e_g2  : e_g1;
    eh_hard  = e_g1 > e_g2 ? e_g1  : e_g2;

    thl_soft = e_g3 > e_g4 ? th_g4 : th_g3;
    thh_soft = e_g3 > e_g4 ? th_g3 : th_g4;
    el_soft  = e_g3 > e_g4 ? e_g4  : e_g3;
    eh_soft  = e_g3 > e_g4 ? e_g3  : e_g4;

    h_th_gl_hard_sig->Fill(thl_hard);
    h_th_gh_hard_sig->Fill(thh_hard);

    h_th_gl_soft_sig->Fill(thl_soft);
    h_th_gh_soft_sig->Fill(thh_soft);

    if(!(i%1000)) cout << i << ": " << el_hard << " " << eh_hard << " " << el_soft << " " << eh_soft << endl;

    if(fabs(thl_hard)>edge) h_e_gl_hard_end_sig->Fill(el_hard);
    else              h_e_gl_hard_bar_sig->Fill(el_hard);

    if(fabs(thh_hard)>edge) h_e_gh_hard_end_sig->Fill(eh_hard);
    else              h_e_gh_hard_bar_sig->Fill(eh_hard);

    if(fabs(thl_soft)>edge) h_e_gl_soft_end_sig->Fill(el_soft);
    else              h_e_gl_soft_bar_sig->Fill(el_soft);

    if(fabs(thh_soft)>edge) h_e_gh_soft_end_sig->Fill(eh_soft);
    else              h_e_gh_soft_bar_sig->Fill(eh_soft);
  }

  cout << "Background" << endl;
  for(int i=0; i<BTot; i++){
    bkg_tree->GetEvent(i);
    if(!(i%1000)) cout << i << ": " << el_hard << " " << eh_hard << " " << el_soft << " " << eh_soft << endl;
    if(abs(de)>0.1 || mbc<5.271) continue;
    if(b0f == 1 || b0f == 5 || b0f == 10 || rndm_pi0 == 1) continue;
    thl_hard = e_g1 > e_g2 ? th_g2 : th_g1;
    thh_hard = e_g1 > e_g2 ? th_g1 : th_g2;
    el_hard  = e_g1 > e_g2 ? e_g2  : e_g1;
    eh_hard  = e_g1 > e_g2 ? e_g1  : e_g2;

    thl_soft = e_g3 > e_g4 ? th_g4 : th_g3;
    thh_soft = e_g3 > e_g4 ? th_g3 : th_g4;
    el_soft  = e_g3 > e_g4 ? e_g4  : e_g3;
    eh_soft  = e_g3 > e_g4 ? e_g3  : e_g4;

    h_th_gl_hard_bkg->Fill(thl_hard);
    h_th_gh_hard_bkg->Fill(thh_hard);

    h_th_gl_soft_bkg->Fill(thl_soft);
    h_th_gh_soft_bkg->Fill(thh_soft);

    if(fabs(thl_hard)>edge) h_e_gl_hard_end_bkg->Fill(el_hard);
    else              h_e_gl_hard_bar_bkg->Fill(el_hard);

    if(fabs(thh_hard)>edge) h_e_gh_hard_end_bkg->Fill(eh_hard);
    else              h_e_gh_hard_bar_bkg->Fill(eh_hard);

    if(fabs(thl_soft)>edge) h_e_gl_soft_end_bkg->Fill(el_soft);
    else              h_e_gl_soft_bar_bkg->Fill(el_soft);

    if(fabs(thh_soft)>edge) h_e_gh_soft_end_bkg->Fill(eh_soft);
    else              h_e_gh_soft_bar_bkg->Fill(eh_soft);
  }

  const bool draw_e   = true;
  const bool draw_bar = true;
  const bool draw_end = true;
  const bool draw_th  = true;

  if(draw_th){
  TCanvas* c_th_gl_hard_sig = new TCanvas("c_th_gl_hard_sig","c_th_gl_hard_sig",600,400);
  c_th_gl_hard_sig->cd();
  c_th_gl_hard_sig->cd();
  c_th_gl_hard_sig->SetGrid();
  h_th_gl_hard_sig->SetStats(0);
  h_th_gl_hard_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_hard_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_hard_sig->SetLineWidth(2);
  h_th_gl_hard_sig->SetTitle("cos(#theta) of low energy #gamma from hard #pi^{0} (sig D^{*0}#pi^{0})");
  h_th_gl_hard_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_hard_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_hard_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_hard_sig->Draw();
  c_th_gl_hard_sig->Print("pics/th_gl_hard_sig.eps");

  TCanvas* c_th_gh_hard_sig = new TCanvas("c_th_gh_hard_sig","c_th_gh_hard_sig",600,400);
  c_th_gh_hard_sig->cd();
  c_th_gh_hard_sig->cd();
  c_th_gh_hard_sig->SetGrid();
  h_th_gh_hard_sig->SetStats(0);
  h_th_gh_hard_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_hard_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_hard_sig->SetLineWidth(2);
  h_th_gh_hard_sig->SetTitle("cos(#theta) of high energy #gamma from hard #pi^{0} (sig D^{*0}#pi^{0})");
  h_th_gh_hard_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_hard_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_hard_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_hard_sig->Draw();
  c_th_gh_hard_sig->Print("pics/th_gh_hard_sig.eps");

  TCanvas* c_th_gl_soft_sig = new TCanvas("c_th_gl_soft_sig","c_th_gl_soft_sig",600,400);
  c_th_gl_soft_sig->cd();
  c_th_gl_soft_sig->cd();
  c_th_gl_soft_sig->SetGrid();
  h_th_gl_soft_sig->SetStats(0);
  h_th_gl_soft_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_soft_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_soft_sig->SetLineWidth(2);
  h_th_gl_soft_sig->SetTitle("cos(#theta) of low energy #gamma from soft #pi^{0} (sig D^{*0}#pi^{0})");
  h_th_gl_soft_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_soft_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_soft_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_soft_sig->Draw();
  c_th_gl_soft_sig->Print("pics/th_gl_soft_sig.eps");

  TCanvas* c_th_gh_soft_sig = new TCanvas("c_th_gh_soft_sig","c_th_gh_soft_sig",600,400);
  c_th_gh_soft_sig->cd();
  c_th_gh_soft_sig->cd();
  c_th_gh_soft_sig->SetGrid();
  h_th_gh_soft_sig->SetStats(0);
  h_th_gh_soft_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_soft_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_soft_sig->SetLineWidth(2);
  h_th_gh_soft_sig->SetTitle("cos(#theta) of high energy #gamma from soft #pi^{0} (sig D^{*0}#pi^{0})");
  h_th_gh_soft_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_soft_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_soft_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_soft_sig->Draw();
  c_th_gh_soft_sig->Print("pics/th_gh_soft_sig.eps");
  }

  if(draw_e){
  if(draw_bar){
  TCanvas* c_e_gl_hard_bar_sig  = new TCanvas("c_e_gl_hard_bar_sig ","c_e_gl_hard_bar_sig ",600,400);
  c_e_gl_hard_bar_sig->cd();
  c_e_gl_hard_bar_sig->SetGrid();
  h_e_gl_hard_bar_sig->SetStats(0);
  h_e_gl_hard_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_hard_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_hard_bar_sig->SetLineWidth(2);
  h_e_gl_hard_bar_sig->SetTitle("Lower E_{#gamma} from hard #pi^{0} for barrel (sig D^{*0}#pi^{0})");
  h_e_gl_hard_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_hard_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_hard_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_hard_bar_sig->Draw();
  c_e_gl_hard_bar_sig->Print("pics/e_gl_hard_bar_sig.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_hard_end_sig  = new TCanvas("c_e_gl_hard_end_sig ","c_e_gl_hard_end_sig ",600,400);
  c_e_gl_hard_end_sig->cd();
  c_e_gl_hard_end_sig->SetGrid();
  h_e_gl_hard_end_sig->SetStats(0);
  h_e_gl_hard_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_hard_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_hard_end_sig->SetLineWidth(2);
  h_e_gl_hard_end_sig->SetTitle("Lower E_{#gamma} from hard #pi^{0} for endcap (sig D^{*0}#pi^{0})");
  h_e_gl_hard_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_hard_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_hard_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_hard_end_sig->Draw();
  c_e_gl_hard_end_sig->Print("pics/e_gl_hard_end_sig.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_hard_bar_sig  = new TCanvas("c_e_gh_hard_bar_sig ","c_e_gh_hard_bar_sig ",600,400);
  c_e_gh_hard_bar_sig->cd();
  c_e_gh_hard_bar_sig->SetGrid();
  h_e_gh_hard_bar_sig->SetStats(0);
  h_e_gh_hard_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_hard_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_hard_bar_sig->SetLineWidth(2);
  h_e_gh_hard_bar_sig->SetTitle("Higher E_{#gamma} from hard #pi^{0} for barrel (sig D^{*0}#pi^{0})");
  h_e_gh_hard_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_hard_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_hard_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_hard_bar_sig->Draw();
  c_e_gh_hard_bar_sig->Print("pics/e_gh_hard_bar_sig.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_hard_end_sig  = new TCanvas("c_e_gh_hard_end_sig ","c_e_gh_hard_end_sig ",600,400);
  c_e_gh_hard_end_sig->cd();
  c_e_gh_hard_end_sig->SetGrid();
  h_e_gh_hard_end_sig->SetStats(0);
  h_e_gh_hard_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_hard_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_hard_end_sig->SetLineWidth(2);
  h_e_gh_hard_end_sig->SetTitle("Higher E_{#gamma} from hard #pi^{0} for endcap (sig D^{*0}#pi^{0})");
  h_e_gh_hard_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_hard_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_hard_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_hard_end_sig->Draw();
  c_e_gh_hard_end_sig->Print("pics/e_gh_hard_end_sig.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gl_soft_bar_sig  = new TCanvas("c_e_gl_soft_bar_sig ","c_e_gl_soft_bar_sig ",600,400);
  c_e_gl_soft_bar_sig->cd();
  c_e_gl_soft_bar_sig->SetGrid();
  h_e_gl_soft_bar_sig->SetStats(0);
  h_e_gl_soft_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_soft_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_soft_bar_sig->SetLineWidth(2);
  h_e_gl_soft_bar_sig->SetTitle("Lower E_{#gamma} from soft #pi^{0} for barrel (sig D^{*0}#pi^{0})");
  h_e_gl_soft_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_soft_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_soft_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_soft_bar_sig->Draw();
  c_e_gl_soft_bar_sig->Print("pics/e_gl_soft_bar_sig.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_soft_end_sig  = new TCanvas("c_e_gl_soft_end_sig ","c_e_gl_soft_end_sig ",600,400);
  c_e_gl_soft_end_sig->cd();
  c_e_gl_soft_end_sig->SetGrid();
  h_e_gl_soft_end_sig->SetStats(0);
  h_e_gl_soft_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_soft_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_soft_end_sig->SetLineWidth(2);
  h_e_gl_soft_end_sig->SetTitle("Lower E_{#gamma} from soft #pi^{0} for endcap (sig D^{*0}#pi^{0})");
  h_e_gl_soft_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_soft_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_soft_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_soft_end_sig->Draw();
  c_e_gl_soft_end_sig->Print("pics/e_gl_soft_end_sig.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_soft_bar_sig  = new TCanvas("c_e_gh_soft_bar_sig ","c_e_gh_soft_bar_sig ",600,400);
  c_e_gh_soft_bar_sig->cd();
  c_e_gh_soft_bar_sig->SetGrid();
  h_e_gh_soft_bar_sig->SetStats(0);
  h_e_gh_soft_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_soft_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_soft_bar_sig->SetLineWidth(2);
  h_e_gh_soft_bar_sig->SetTitle("Higher E_{#gamma} from soft #pi^{0} for barrel (sig D^{*0}#pi^{0})");
  h_e_gh_soft_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_soft_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_soft_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_soft_bar_sig->Draw();
  c_e_gh_soft_bar_sig->Print("pics/e_gh_soft_bar_sig.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_soft_end_sig  = new TCanvas("c_e_gh_soft_end_sig ","c_e_gh_soft_end_sig ",600,400);
  c_e_gh_soft_end_sig->cd();
  c_e_gh_soft_end_sig->SetGrid();
  h_e_gh_soft_end_sig->SetStats(0);
  h_e_gh_soft_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_soft_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_soft_end_sig->SetLineWidth(2);
  h_e_gh_soft_end_sig->SetTitle("Higher E_{#gamma} from soft #pi^{0} for endcap (sig D^{*0}#pi^{0})");
  h_e_gh_soft_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_soft_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_soft_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_soft_end_sig->Draw();
  c_e_gh_soft_end_sig->Print("pics/e_gh_soft_end_sig.eps");
  }
  }

  if(draw_th){
  TCanvas* c_th_gl_hard_bkg = new TCanvas("c_th_gl_hard_bkg","c_th_gl_hard_bkg",600,400);
  c_th_gl_hard_bkg->cd();
  c_th_gl_hard_bkg->SetGrid();
  h_th_gl_hard_bkg->SetStats(0);
  h_th_gl_hard_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_hard_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_hard_bkg->SetLineWidth(2);
  h_th_gl_hard_bkg->SetTitle("cos(#theta) of low energy #gamma from hard #pi^{0} (bkg D^{*0}#pi^{0})");
  h_th_gl_hard_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_hard_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_hard_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_hard_bkg->Draw();
  c_th_gl_hard_bkg->Print("pics/th_gl_hard_bkg.eps");

  TCanvas* c_th_gh_hard_bkg = new TCanvas("c_th_gh_hard_bkg","c_th_gh_hard_bkg",600,400);
  c_th_gh_hard_bkg->cd();
  c_th_gh_hard_bkg->SetGrid();
  h_th_gh_hard_bkg->SetStats(0);
  h_th_gh_hard_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_hard_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_hard_bkg->SetLineWidth(2);
  h_th_gh_hard_bkg->SetTitle("cos(#theta) of high energy #gamma from hard #pi^{0} (bkg D^{*0}#pi^{0})");
  h_th_gh_hard_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_hard_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_hard_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_hard_bkg->Draw();
  c_th_gh_hard_bkg->Print("pics/th_gh_hard_bkg.eps");

  TCanvas* c_th_gl_soft_bkg = new TCanvas("c_th_gl_soft_bkg","c_th_gl_soft_bkg",600,400);
  c_th_gl_soft_bkg->cd();
  c_th_gl_soft_bkg->SetGrid();
  h_th_gl_soft_bkg->SetStats(0);
  h_th_gl_soft_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_soft_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_soft_bkg->SetLineWidth(2);
  h_th_gl_soft_bkg->SetTitle("cos(#theta) of low energy #gamma from soft #pi^{0} (bkg D^{*0}#pi^{0})");
  h_th_gl_soft_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_soft_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_soft_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_soft_bkg->Draw();
  c_th_gl_soft_bkg->Print("pics/th_gl_soft_bkg.eps");

  TCanvas* c_th_gh_soft_bkg = new TCanvas("c_th_gh_soft_bkg","c_th_gh_soft_bkg",600,400);
  c_th_gh_soft_bkg->cd();
  c_th_gh_soft_bkg->SetGrid();
  h_th_gh_soft_bkg->SetStats(0);
  h_th_gh_soft_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_soft_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_soft_bkg->SetLineWidth(2);
  h_th_gh_soft_bkg->SetTitle("cos(#theta) of high energy #gamma from soft #pi^{0} (bkg D^{*0}#pi^{0})");
  h_th_gh_soft_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_soft_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_soft_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_soft_bkg->Draw();
  c_th_gh_soft_bkg->Print("pics/th_gh_soft_bkg.eps");
  }

  if(draw_e){
  if(draw_bar){
  TCanvas* c_e_gl_hard_bar_bkg  = new TCanvas("c_e_gl_hard_bar_bkg ","c_e_gl_hard_bar_bkg ",600,400);
  c_e_gl_hard_bar_bkg->cd();
  c_e_gl_hard_bar_bkg->SetGrid();
  h_e_gl_hard_bar_bkg->SetStats(0);
  h_e_gl_hard_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_hard_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_hard_bar_bkg->SetLineWidth(2);
  h_e_gl_hard_bar_bkg->SetTitle("Lower E_{#gamma} from hard #pi^{0} for barrel (bkg D^{*0}#pi^{0})");
  h_e_gl_hard_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_hard_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_hard_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_hard_bar_bkg->Draw();
  c_e_gl_hard_bar_bkg->Print("pics/e_gl_hard_bar_bkg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_hard_end_bkg  = new TCanvas("c_e_gl_hard_end_bkg ","c_e_gl_hard_end_bkg ",600,400);
  c_e_gl_hard_end_bkg->cd();
  c_e_gl_hard_end_bkg->SetGrid();
  h_e_gl_hard_end_bkg->SetStats(0);
  h_e_gl_hard_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_hard_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_hard_end_bkg->SetLineWidth(2);
  h_e_gl_hard_end_bkg->SetTitle("Lower E_{#gamma} from hard #pi^{0} for endcap (bkg D^{*0}#pi^{0})");
  h_e_gl_hard_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_hard_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_hard_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_hard_end_bkg->Draw();
  c_e_gl_hard_end_bkg->Print("pics/e_gl_hard_end_bkg.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_hard_bar_bkg  = new TCanvas("c_e_gh_hard_bar_bkg ","c_e_gh_hard_bar_bkg ",600,400);
  c_e_gh_hard_bar_bkg->cd();
  c_e_gh_hard_bar_bkg->SetGrid();
  h_e_gh_hard_bar_bkg->SetStats(0);
  h_e_gh_hard_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_hard_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_hard_bar_bkg->SetLineWidth(2);
  h_e_gh_hard_bar_bkg->SetTitle("Higher E_{#gamma} from hard #pi^{0} for barrel (bkg D^{*0}#pi^{0})");
  h_e_gh_hard_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_hard_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_hard_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_hard_bar_bkg->Draw();
  c_e_gh_hard_bar_bkg->Print("pics/e_gh_hard_bar_bkg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_hard_end_bkg  = new TCanvas("c_e_gh_hard_end_bkg ","c_e_gh_hard_end_bkg ",600,400);
  c_e_gh_hard_end_bkg->cd();
  c_e_gh_hard_end_bkg->SetGrid();
  h_e_gh_hard_end_bkg->SetStats(0);
  h_e_gh_hard_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_hard_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_hard_end_bkg->SetLineWidth(2);
  h_e_gh_hard_end_bkg->SetTitle("Higher E_{#gamma} from hard #pi^{0} for endcap (bkg D^{*0}#pi^{0})");
  h_e_gh_hard_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_hard_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_hard_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_hard_end_bkg->Draw();
  c_e_gh_hard_end_bkg->Print("pics/e_gh_hard_end_bkg.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gl_soft_bar_bkg  = new TCanvas("c_e_gl_soft_bar_bkg ","c_e_gl_soft_bar_bkg ",600,400);
  c_e_gl_soft_bar_bkg->cd();
  c_e_gl_soft_bar_bkg->SetGrid();
  h_e_gl_soft_bar_bkg->SetStats(0);
  h_e_gl_soft_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_soft_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_soft_bar_bkg->SetLineWidth(2);
  h_e_gl_soft_bar_bkg->SetTitle("Lower E_{#gamma} from soft #pi^{0} for barrel (bkg D^{*0}#pi^{0})");
  h_e_gl_soft_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_soft_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_soft_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_soft_bar_bkg->Draw();
  c_e_gl_soft_bar_bkg->Print("pics/e_gl_soft_bar_bkg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_soft_end_bkg  = new TCanvas("c_e_gl_soft_end_bkg ","c_e_gl_soft_end_bkg ",600,400);
  c_e_gl_soft_end_bkg->cd();
  c_e_gl_soft_end_bkg->SetGrid();
  h_e_gl_soft_end_bkg->SetStats(0);
  h_e_gl_soft_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_soft_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_soft_end_bkg->SetLineWidth(2);
  h_e_gl_soft_end_bkg->SetTitle("Lower E_{#gamma} from soft #pi^{0} for endcap (bkg D^{*0}#pi^{0})");
  h_e_gl_soft_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_soft_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_soft_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_soft_end_bkg->Draw();
  c_e_gl_soft_end_bkg->Print("pics/e_gl_soft_end_bkg.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_soft_bar_bkg  = new TCanvas("c_e_gh_soft_bar_bkg ","c_e_gh_soft_bar_bkg ",600,400);
  c_e_gh_soft_bar_bkg->cd();
  c_e_gh_soft_bar_bkg->SetGrid();
  h_e_gh_soft_bar_bkg->SetStats(0);
  h_e_gh_soft_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_soft_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_soft_bar_bkg->SetLineWidth(2);
  h_e_gh_soft_bar_bkg->SetTitle("Higher E_{#gamma} from soft #pi^{0} for barrel (bkg D^{*0}#pi^{0})");
  h_e_gh_soft_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_soft_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_soft_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_soft_bar_bkg->Draw();
  c_e_gh_soft_bar_bkg->Print("pics/e_gh_soft_bar_bkg.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_soft_end_bkg  = new TCanvas("c_e_gh_soft_end_bkg ","c_e_gh_soft_end_bkg ",600,400);
  c_e_gh_soft_end_bkg->cd();
  c_e_gh_soft_end_bkg->SetGrid();
  h_e_gh_soft_end_bkg->SetStats(0);
  h_e_gh_soft_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_soft_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_soft_end_bkg->SetLineWidth(2);
  h_e_gh_soft_end_bkg->SetTitle("Higher E_{#gamma} from soft #pi^{0} for endcap (bkg D^{*0}#pi^{0})");
  h_e_gh_soft_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_soft_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_soft_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_soft_end_bkg->Draw();
  c_e_gh_soft_end_bkg->Print("pics/e_gh_soft_end_bkg.eps");
  }
  }

  return;
}
