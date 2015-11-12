
void draw_egamma_pi0(const bool small_range = false){
  TChain* sig_tree = new TChain("TEvent","TEvent");
  sig_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcPi0_s8_m1_h0m10.root");
  TChain* bkg_tree = new TChain("TEvent","TEvent");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_0_10_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_1_11_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_2_12_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_3_13_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_4_14_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_5_15_m1_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_0_10_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_1_11_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_2_12_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_3_13_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_4_14_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charged_5_15_m1_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_0_10_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_1_11_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_2_12_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_3_13_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_4_14_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_mixed_5_15_m1_h0m10.root");

  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_0_10_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_1_11_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_2_12_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_3_13_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_4_14_m1_h0m10.root");
  bkg_tree->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_5_15_m1_h0m10.root");

  const int STot = sig_tree->GetEntries();
  const int BTot = bkg_tree->GetEntries();

  // Signal histograms //
  TH1D* h_th_gl_sig = new TH1D("h_th_gl_sig","h_th_gl_sig",60,-1.,1.);
  TH1D* h_th_gh_sig = new TH1D("h_th_gh_sig","h_th_gh_sig",60,-1.,1.);

  const double gl_limit = small_range ? 1. : 2.;
  const double gh_limit = small_range ? 1. : 3.5;
  TH1D* h_e_gl_bar_sig = new TH1D("h_e_gl_bar_sig","h_e_gl_bar_sig",100,0.04,gl_limit);
  TH1D* h_e_gl_end_sig = new TH1D("h_e_gl_end_sig","h_e_gl_end_sig",100,0.04,gl_limit);
  TH1D* h_e_gh_bar_sig = new TH1D("h_e_gh_bar_sig","h_e_gh_bar_sig",100,0.04,gh_limit);
  TH1D* h_e_gh_end_sig = new TH1D("h_e_gh_end_sig","h_e_gh_end_sig",100,0.04,gh_limit);

  const double ppi0_limit = small_range ? 2.4 : 4.;
  TH1D* h_helpi0_sig  = new TH1D("h_helpi0_sig","h_helpi0_sig",100,0.,1.);
  TH1D* h_ppi0_sig = new TH1D("h_ppi0_sig","h_ppi0_sig",100,1.4,ppi0_limit);

  const double ppi_limit = small_range ? 1. : 3.5;
  TH1D* h_ppip_sig = new TH1D("h_ppip_sig","h_ppip_sig",100,0.,ppi_limit);
  TH1D* h_ppim_sig = new TH1D("h_ppim_sig","h_ppim_sig",100,0.,ppi_limit);

  const double ptpi_limit = small_range ? 1. : 2.5;
  TH1D* h_ptpip_sig = new TH1D("h_ptpip_sig","h_ptpip_sig",100,0.,ptpi_limit);
  TH1D* h_ptpim_sig = new TH1D("h_ptpim_sig","h_ptpim_sig",100,0.,ptpi_limit);

  // Background histograms //
  TH1D* h_th_gl_bkg = new TH1D("h_th_gl_bkg","h_th_gl_bkg",60,-1.,1.);
  TH1D* h_th_gh_bkg = new TH1D("h_th_gh_bkg","h_th_gh_bkg",60,-1.,1.);

  TH1D* h_e_gl_bar_bkg = new TH1D("h_e_gl_bar_bkg","h_e_gl_bar_bkg",100,0.04,gl_limit);
  TH1D* h_e_gl_end_bkg = new TH1D("h_e_gl_end_bkg","h_e_gl_end_bkg",100,0.04,gl_limit);
  TH1D* h_e_gh_bar_bkg = new TH1D("h_e_gh_bar_bkg","h_e_gh_bar_bkg",100,0.04,gh_limit);
  TH1D* h_e_gh_end_bkg = new TH1D("h_e_gh_end_bkg","h_e_gh_end_bkg",100,0.04,gh_limit);

  TH1D* h_helpi0_bkg  = new TH1D("h_helpi0_bkg","h_helpi0_bkg",100,0.,1.);
  TH1D* h_ppi0_bkg = new TH1D("h_ppi0_bkg","h_ppi0_bkg",100,1.4,ppi0_limit);

  TH1D* h_ppip_bkg = new TH1D("h_ppip_bkg","h_ppip_bkg",100,0.,ppi_limit);
  TH1D* h_ppim_bkg = new TH1D("h_ppim_bkg","h_ppim_bkg",100,0.,ppi_limit);

  TH1D* h_ptpip_bkg = new TH1D("h_ptpip_bkg","h_ptpip_bkg",100,0.,ptpi_limit);
  TH1D* h_ptpim_bkg = new TH1D("h_ptpim_bkg","h_ptpim_bkg",100,0.,ptpi_limit);

  const double edge = 0.7;
  int b0f,rndm_pi0;
  double e_g1,e_g2;
  double th_g1,th_g2;
  double de,mbc;
  double pt_pip,pt_pim;
  double pz_pip,pz_pim;
  double p_h0;
  double hel_pi0;

  sig_tree->SetBranchAddress("de",&de);
  sig_tree->SetBranchAddress("mbc",&mbc);
  sig_tree->SetBranchAddress("b0f",&b0f);
  sig_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  sig_tree->SetBranchAddress("e_g1",&e_g1);
  sig_tree->SetBranchAddress("e_g2",&e_g2);
  sig_tree->SetBranchAddress("th_g1",&th_g1);
  sig_tree->SetBranchAddress("th_g2",&th_g2);
  sig_tree->SetBranchAddress("pt_pip",&pt_pip);
  sig_tree->SetBranchAddress("pt_pim",&pt_pim);
  sig_tree->SetBranchAddress("pz_pip",&pz_pip);
  sig_tree->SetBranchAddress("pz_pim",&pz_pim);
  sig_tree->SetBranchAddress("p_h0",&p_h0);
  sig_tree->SetBranchAddress("hel_h0",&hel_pi0);

  bkg_tree->SetBranchAddress("de",&de);
  bkg_tree->SetBranchAddress("mbc",&mbc);
  bkg_tree->SetBranchAddress("b0f",&b0f);
  bkg_tree->SetBranchAddress("rndm_pi0",&rndm_pi0);
  bkg_tree->SetBranchAddress("e_g1",&e_g1);
  bkg_tree->SetBranchAddress("e_g2",&e_g2);
  bkg_tree->SetBranchAddress("th_g1",&th_g1);
  bkg_tree->SetBranchAddress("th_g2",&th_g2);
  bkg_tree->SetBranchAddress("pt_pip",&pt_pip);
  bkg_tree->SetBranchAddress("pt_pim",&pt_pim);
  bkg_tree->SetBranchAddress("pz_pip",&pz_pip);
  bkg_tree->SetBranchAddress("pz_pim",&pz_pim);
  bkg_tree->SetBranchAddress("p_h0",&p_h0);
  bkg_tree->SetBranchAddress("hel_h0",&hel_pi0);

  double el,eh;
  double thl,thh;
  double p_pip,p_pim;
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

    p_pip = sqrt(pz_pip*pz_pip+pt_pip*pt_pip);
    p_pim = sqrt(pz_pim*pz_pim+pt_pim*pt_pim);

    h_helpi0_sig->Fill(fabs(hel_pi0));
    h_ppi0_sig->Fill(p_h0);

    h_ppip_sig->Fill(p_pip);
    h_ppim_sig->Fill(p_pim);

    h_ptpip_sig->Fill(pt_pip);
    h_ptpim_sig->Fill(pt_pim);

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

    p_pip = sqrt(pz_pip*pz_pip+pt_pip*pt_pip);
    p_pim = sqrt(pz_pim*pz_pim+pt_pim*pt_pim);

    h_helpi0_bkg->Fill(fabs(hel_pi0));
    h_ppi0_bkg->Fill(p_h0);

    h_ppip_bkg->Fill(p_pip);
    h_ppim_bkg->Fill(p_pim);

    h_ptpip_bkg->Fill(pt_pip);
    h_ptpim_bkg->Fill(pt_pim);

    if(fabs(thl)>edge) h_e_gl_end_bkg->Fill(el);
    else               h_e_gl_bar_bkg->Fill(el);

    if(fabs(thh)>edge) h_e_gh_end_bkg->Fill(eh);
    else               h_e_gh_bar_bkg->Fill(eh);
  }

  const bool draw_e   = false;
  const bool draw_bar = true;
  const bool draw_end = true;
  const bool draw_th  = false;
  const bool draw_mom = true;

  if(draw_mom){
  TCanvas* c_helpi0_sig = new TCanvas("c_helpi0_sig","c_helpi0_sig",600,400);
  c_helpi0_sig->cd();
  c_helpi0_sig->SetGrid();
  h_helpi0_sig->SetStats(0);
  h_helpi0_sig->GetYaxis()->SetLabelSize(0.06);
  h_helpi0_sig->GetXaxis()->SetLabelSize(0.06);
  h_helpi0_sig->SetLineWidth(2);
  h_helpi0_sig->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #pi^{0} (sig D^{0}#pi^{0})");
  h_helpi0_sig->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_helpi0_sig->GetXaxis()->SetTitleSize(0.06);
  h_helpi0_sig->GetXaxis()->SetTitleOffset(0.75);
  h_helpi0_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/helpi0_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/helpi0_sig_pi0_s.eps");

  TCanvas* c_helpi0_bkg = new TCanvas("c_helpi0_bkg","c_helpi0_bkg",600,400);
  c_helpi0_bkg->cd();
  c_helpi0_bkg->SetGrid();
  h_helpi0_bkg->SetStats(0);
  h_helpi0_bkg->GetYaxis()->SetLabelSize(0.06);
  h_helpi0_bkg->GetXaxis()->SetLabelSize(0.06);
  h_helpi0_bkg->SetLineWidth(2);
  h_helpi0_bkg->SetTitle("(E_{#gamma}^{H}-E_{#gamma}^{L})/(E_{#gamma}^{H}+E_{#gamma}^{L}) for #pi^{0} (bkg D^{0}#pi^{0})");
  h_helpi0_bkg->GetXaxis()->SetTitle("(Eh-El)/(Eh+El)");
  h_helpi0_bkg->GetXaxis()->SetTitleSize(0.06);
  h_helpi0_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_helpi0_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/helpi0_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/helpi0_bkg_pi0_s.eps");

  TCanvas* c_ppi0_sig = new TCanvas("c_ppi0_sig","c_ppi0_sig",600,400);
  c_ppi0_sig->cd();
  c_ppi0_sig->SetGrid();
  h_ppi0_sig->SetStats(0);
  h_ppi0_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppi0_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppi0_sig->SetLineWidth(2);
  h_ppi0_sig->SetTitle("Lab frame momentum of #pi^{0} (sig D^{0}#pi^{0})");
  h_ppi0_sig->GetXaxis()->SetTitle("p (GeV)");
  h_ppi0_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppi0_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppi0_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppi0_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppi0_sig_pi0_s.eps");

  TCanvas* c_ppi0_bkg = new TCanvas("c_ppi0_bkg","c_ppi0_bkg",600,400);
  c_ppi0_bkg->cd();
  c_ppi0_bkg->SetGrid();
  h_ppi0_bkg->SetStats(0);
  h_ppi0_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppi0_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppi0_bkg->SetLineWidth(2);
  h_ppi0_bkg->SetTitle("Lab frame momentum of #pi^{0} (bkg D^{0}#pi^{0})");
  h_ppi0_bkg->GetXaxis()->SetTitle("p (GeV)");
  h_ppi0_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppi0_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppi0_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppi0_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppi0_bkg_pi0_s.eps");

  TCanvas* c_ppip_sig = new TCanvas("c_ppip_sig","c_ppip_sig",600,400);
  c_ppip_sig->cd();
  c_ppip_sig->SetGrid();
  h_ppip_sig->SetStats(0);
  h_ppip_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppip_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppip_sig->SetLineWidth(2);
  h_ppip_sig->SetTitle("Lab frame momentum of #pi^{+} form D^{0} (sig D^{0}#pi^{0})");
  h_ppip_sig->GetXaxis()->SetTitle("p (GeV)");
  h_ppip_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppip_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppip_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppip_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppip_sig_pi0_s.eps");

  TCanvas* c_ppip_bkg = new TCanvas("c_ppip_bkg","c_ppip_bkg",600,400);
  c_ppip_bkg->cd();
  c_ppip_bkg->SetGrid();
  h_ppip_bkg->SetStats(0);
  h_ppip_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppip_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppip_bkg->SetLineWidth(2);
  h_ppip_bkg->SetTitle("Lab frame momentum of #pi^{+} form D^{0} (bkg D^{0}#pi^{0})");
  h_ppip_bkg->GetXaxis()->SetTitle("p (GeV)");
  h_ppip_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppip_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppip_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppip_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppip_bkg_pi0_s.eps");

  TCanvas* c_ppim_sig = new TCanvas("c_ppim_sig","c_ppim_sig",600,400);
  c_ppim_sig->cd();
  c_ppim_sig->SetGrid();
  h_ppim_sig->SetStats(0);
  h_ppim_sig->GetYaxis()->SetLabelSize(0.06);
  h_ppim_sig->GetXaxis()->SetLabelSize(0.06);
  h_ppim_sig->SetLineWidth(2);
  h_ppim_sig->SetTitle("Lab frame momentum of #pi^{-} form D^{0} (sig D^{0}#pi^{0})");
  h_ppim_sig->GetXaxis()->SetTitle("p (GeV)");
  h_ppim_sig->GetXaxis()->SetTitleSize(0.06);
  h_ppim_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ppim_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppim_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppim_sig_pi0_s.eps");

  TCanvas* c_ppim_bkg = new TCanvas("c_ppim_bkg","c_ppim_bkg",600,400);
  c_ppim_bkg->cd();
  c_ppim_bkg->SetGrid();
  h_ppim_bkg->SetStats(0);
  h_ppim_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ppim_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ppim_bkg->SetLineWidth(2);
  h_ppim_bkg->SetTitle("Lab frame momentum of #pi^{-} form D^{0} (bkg D^{0}#pi^{0})");
  h_ppim_bkg->GetXaxis()->SetTitle("p (GeV)");
  h_ppim_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ppim_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ppim_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ppim_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/ppim_bkg_pi0_s.eps");

  TCanvas* c_ptpip_sig = new TCanvas("c_ptpip_sig","c_ptpip_sig",600,400);
  c_ptpip_sig->cd();
  c_ptpip_sig->SetGrid();
  h_ptpip_sig->SetStats(0);
  h_ptpip_sig->GetYaxis()->SetLabelSize(0.06);
  h_ptpip_sig->GetXaxis()->SetLabelSize(0.06);
  h_ptpip_sig->SetLineWidth(2);
  h_ptpip_sig->SetTitle("Lab frame transverse momentum of #pi^{+} form D^{0} (sig D^{0}#pi^{0})");
  h_ptpip_sig->GetXaxis()->SetTitle("pt (GeV)");
  h_ptpip_sig->GetXaxis()->SetTitleSize(0.06);
  h_ptpip_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ptpip_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ptpip_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/ptpip_sig_pi0_s.eps");

  TCanvas* c_ptpip_bkg = new TCanvas("c_ptpip_bkg","c_ptpip_bkg",600,400);
  c_ptpip_bkg->cd();
  c_ptpip_bkg->SetGrid();
  h_ptpip_bkg->SetStats(0);
  h_ptpip_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ptpip_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ptpip_bkg->SetLineWidth(2);
  h_ptpip_bkg->SetTitle("Lab frame transverse momentum of #pi^{+} form D^{0} (bkg D^{0}#pi^{0})");
  h_ptpip_bkg->GetXaxis()->SetTitle("pt (GeV)");
  h_ptpip_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ptpip_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ptpip_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ptpip_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/ptpip_bkg_pi0_s.eps");

  TCanvas* c_ptpim_sig = new TCanvas("c_ptpim_sig","c_ptpim_sig",600,400);
  c_ptpim_sig->cd();
  c_ptpim_sig->SetGrid();
  h_ptpim_sig->SetStats(0);
  h_ptpim_sig->GetYaxis()->SetLabelSize(0.06);
  h_ptpim_sig->GetXaxis()->SetLabelSize(0.06);
  h_ptpim_sig->SetLineWidth(2);
  h_ptpim_sig->SetTitle("Lab frame transverse momentum of #pi^{-} form D^{0} (sig D^{0}#pi^{0})");
  h_ptpim_sig->GetXaxis()->SetTitle("pt (GeV)");
  h_ptpim_sig->GetXaxis()->SetTitleSize(0.06);
  h_ptpim_sig->GetXaxis()->SetTitleOffset(0.75);
  h_ptpim_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ptpim_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/ptpim_sig_pi0_s.eps");

  TCanvas* c_ptpim_bkg = new TCanvas("c_ptpim_bkg","c_ptpim_bkg",600,400);
  c_ptpim_bkg->cd();
  c_ptpim_bkg->SetGrid();
  h_ptpim_bkg->SetStats(0);
  h_ptpim_bkg->GetYaxis()->SetLabelSize(0.06);
  h_ptpim_bkg->GetXaxis()->SetLabelSize(0.06);
  h_ptpim_bkg->SetLineWidth(2);
  h_ptpim_bkg->SetTitle("Lab frame transverse momentum of #pi^{-} form D^{0} (bkg D^{0}#pi^{0})");
  h_ptpim_bkg->GetXaxis()->SetTitle("pt(GeV)");
  h_ptpim_bkg->GetXaxis()->SetTitleSize(0.06);
  h_ptpim_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_ptpim_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/ptpim_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/ptpim_bkg_pi0_s.eps");
  }

  if(draw_th){
  TCanvas* c_th_gl_sig = new TCanvas("c_th_gl_sig","c_th_gl_sig",600,400);
  c_th_gl_sig->cd();
  c_th_gl_sig->SetGrid();
  h_th_gl_sig->SetStats(0);
  h_th_gl_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gl_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gl_sig->SetLineWidth(2);
  h_th_gl_sig->SetTitle("cos(#theta) of low energy #gamma from #pi^{0} (sig D^{0}#pi^{0})");
  h_th_gl_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/th_gl_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/th_gl_sig_pi0_s.eps");

  TCanvas* c_th_gh_sig = new TCanvas("c_th_gh_sig","c_th_gh_sig",600,400);
  c_th_gh_sig->cd();
  c_th_gh_sig->SetGrid();
  h_th_gh_sig->SetStats(0);
  h_th_gh_sig->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_sig->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_sig->SetLineWidth(2);
  h_th_gh_sig->SetTitle("cos(#theta) of high energy #gamma from #pi^{0} (sig D^{0}#pi^{0})");
  h_th_gh_sig->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_sig->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_sig->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/th_gh_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/th_gh_sig_pi0_s.eps");
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
  h_e_gl_bar_sig->SetTitle("Lower E_{#gamma} from #pi^{0} for barrel (sig D^{0}#pi^{0})");
  h_e_gl_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gl_bar_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gl_bar_sig_pi0_s.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_sig  = new TCanvas("c_e_gl_end_sig ","c_e_gl_end_sig ",600,400);
  c_e_gl_end_sig->cd();
  c_e_gl_end_sig->SetGrid();
  h_e_gl_end_sig->SetStats(0);
  h_e_gl_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_sig->SetLineWidth(2);
  h_e_gl_end_sig->SetTitle("Lower E_{#gamma} from #pi^{0} for endcap (sig D^{0}#pi^{0})");
  h_e_gl_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gl_end_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gl_end_sig_pi0_s.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_sig  = new TCanvas("c_e_gh_bar_sig ","c_e_gh_bar_sig ",600,400);
  c_e_gh_bar_sig->cd();
  c_e_gh_bar_sig->SetGrid();
  h_e_gh_bar_sig->SetStats(0);
  h_e_gh_bar_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_sig->SetLineWidth(2);
  h_e_gh_bar_sig->SetTitle("Higher E_{#gamma} from #pi^{0} for barrel (sig D^{0}#pi^{0})");
  h_e_gh_bar_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gh_bar_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gh_bar_sig_pi0_s.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_sig  = new TCanvas("c_e_gh_end_sig ","c_e_gh_end_sig ",600,400);
  c_e_gh_end_sig->cd();
  c_e_gh_end_sig->SetGrid();
  h_e_gh_end_sig->SetStats(0);
  h_e_gh_end_sig->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_sig->SetLineWidth(2);
  h_e_gh_end_sig->SetTitle("Higher E_{#gamma} from #pi^{0} for endcap (sig D^{0}#pi^{0})");
  h_e_gh_end_sig->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_sig->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_sig->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_sig->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gh_end_sig_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gh_end_sig_pi0_s.eps");
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
  h_th_gl_bkg->SetTitle("cos(#theta) of low energy #gamma from #pi^{0} (bkg D^{0}#pi^{0})");
  h_th_gl_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gl_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gl_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gl_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/th_gl_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/th_gl_bkg_pi0_s.eps");

  TCanvas* c_th_gh_bkg = new TCanvas("c_th_gh_bkg","c_th_gh_bkg",600,400);
  c_th_gh_bkg->cd();
  c_th_gh_bkg->SetGrid();
  h_th_gh_bkg->SetStats(0);
  h_th_gh_bkg->GetYaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetLabelSize(0.06);
  h_th_gh_bkg->SetLineWidth(2);
  h_th_gh_bkg->SetTitle("cos(#theta) of high energy #gamma from #pi^{0} (bkg D^{0}#pi^{0})");
  h_th_gh_bkg->GetXaxis()->SetTitle("cos(#theta)");
  h_th_gh_bkg->GetXaxis()->SetTitleSize(0.06);
  h_th_gh_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_th_gh_bkg->Draw();
  c_th_gh_bkg->Print("pics/th_gh_bkg_pi0.eps");
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
  h_e_gl_bar_bkg->SetTitle("Lower E_{#gamma} from #pi^{0} for barrel (bkg D^{0}#pi^{0})");
  h_e_gl_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_bar_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gl_bar_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gl_bar_bkg_pi0_s.eps");
  }

  if(draw_end){
  TCanvas* c_e_gl_end_bkg  = new TCanvas("c_e_gl_end_bkg ","c_e_gl_end_bkg ",600,400);
  c_e_gl_end_bkg->cd();
  c_e_gl_end_bkg->SetGrid();
  h_e_gl_end_bkg->SetStats(0);
  h_e_gl_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gl_end_bkg->SetLineWidth(2);
  h_e_gl_end_bkg->SetTitle("Lower E_{#gamma} from #pi^{0} for endcap (bkg D^{0}#pi^{0})");
  h_e_gl_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gl_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gl_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gl_end_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gl_end_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gl_end_bkg_pi0_s.eps");
  }

  if(draw_bar){
  TCanvas* c_e_gh_bar_bkg  = new TCanvas("c_e_gh_bar_bkg ","c_e_gh_bar_bkg ",600,400);
  c_e_gh_bar_bkg->cd();
  c_e_gh_bar_bkg->SetGrid();
  h_e_gh_bar_bkg->SetStats(0);
  h_e_gh_bar_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_bar_bkg->SetLineWidth(2);
  h_e_gh_bar_bkg->SetTitle("Higher E_{#gamma} from #pi^{0} for barrel (bkg D^{0}#pi^{0})");
  h_e_gh_bar_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_bar_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_bar_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_bar_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("pics/e_gh_bar_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gh_bar_bkg_pi0_s.eps");
  }

  if(draw_end){
  TCanvas* c_e_gh_end_bkg  = new TCanvas("c_e_gh_end_bkg ","c_e_gh_end_bkg ",600,400);
  c_e_gh_end_bkg->cd();
  c_e_gh_end_bkg->SetGrid();
  h_e_gh_end_bkg->SetStats(0);
  h_e_gh_end_bkg->GetYaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetLabelSize(0.06);
  h_e_gh_end_bkg->SetLineWidth(2);
  h_e_gh_end_bkg->SetTitle("Higher E_{#gamma} from #pi^{0} for endcap (bkg D^{0}#pi^{0})");
  h_e_gh_end_bkg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h_e_gh_end_bkg->GetXaxis()->SetTitleSize(0.06);
  h_e_gh_end_bkg->GetXaxis()->SetTitleOffset(0.75);
  h_e_gh_end_bkg->Draw();
  if(!small_range) c_helpi0_sig->Print("picse_gh_end_bkg_pi0.eps");
  else             c_helpi0_sig->Print("pics/e_gh_end_bkg_pi0_s.eps");
  }
  }

  return;
}
