void draw_k(void){
  const int nbins = 8;
  const double bins[8]         = {1,2,3,4,5,6,7,8};
  const double bins_err[8]     = {0,0,0,0,0,0,0,0};
  const double model_kp[8]     = {16.88,11.85,9.62,7.42,9.07,3.08,10.54,7.79};
  const double model_km[8]     = { 8.80, 2.86,1.14,1.50,4.25,1.01, 2.32,1.78};

  const double pi0_kp[8]       = {16.90,11.95,9.60,7.54,9.35,3.02, 9.67,7.18};
  const double pi0_kp_err[8]   = { 0.21, 0.20,0.09,0.08,0.09,0.05, 0.09,0.08};
  const double pi0_km[8]       = { 9.43, 2.96,1.33,1.50,4.26,1.13, 2.29,1.90};
  const double pi0_km_err[8]   = { 0.20, 0.07,0.06,0.05,0.07,0.04, 0.06,0.05};

  const double omega_kp[8]     = {17.03,11.87,9.60,7.65,9.27,3.08,10.03,7.19};
  const double omega_kp_err[8] = { 0.24, 0.22,0.20,0.09,0.21,0.06, 0.21,0.09};
  const double omega_km[8]     = { 9.10, 2.98,1.16,1.50,4.40,1.07, 2.28,1.79};
  const double omega_km_err[8] = { 0.22, 0.08,0.07,0.07,0.09,0.05, 0.08,0.07};

  const double cleo_kp[8]      = {16.5, 12.4, 9.9, 7.1, 8.0, 3.0,  9.8, 7.7};
  const double cleo_kp_err[8]  = { 0.5,  0.4, 0.4, 0.3, 0.4, 0.2,  0.4, 0.4};
  const double cleo_km[8]      = { 8.8,  2.9, 1.6, 1.8, 4.0, 1.3,  3.2, 2.0};
  const double cleo_km_err[8]  = { 0.4,  0.2, 0.2, 0.2, 0.3, 0.2,  0.2, 0.2};

  const double bdpi_kp[8]      = {17.4, 12.5,10.5, 7.3, 9.4, 2.8, 10.2, 7.5};
  const double bdpi_kp_err[8]  = { 0.3,  0.3, 0.3, 0.2, 0.2, 0.2,  0.3, 0.2};
  const double bdpi_km[8]      = { 7.9,  2.4, 1.2, 1.7, 4.3, 1.2,  2.6, 1.3};
  const double bdpi_km_err[8]  = { 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2};

  const double cross_matrix[9][9] = {
           0.331,93.684, 2.017, 0.017, 0.017, 0.033, 0.099, 0.198, 3.605,
           0.501, 1.762,96.951, 0.596, 0.054, 0.014, 0.041, 0.041, 0.041,
           0.180, 0.060, 2.226,95.940, 1.353, 0.030, 0.030, 0.120, 0.060,
           0.277, 0.012, 0.035, 0.611,98.110, 0.876, 0.035, 0.023, 0.023,
           0.394, 0.069, 0.000, 0.034, 1.062,96.299, 2.022, 0.069, 0.051,
           0.426, 0.066, 0.016, 0.000, 0.033, 1.968,94.163, 3.132, 0.197,
           0.426, 0.178, 0.047, 0.000, 0.059, 0.083, 2.275,93.733, 3.199,
           0.203, 1.131, 0.044, 0.032, 0.006, 0.019, 0.064, 1.811,96.689};

  double avmc_kp[8];
  double avmc_kp_err[8];
  double avmc_km[8];
  double avmc_km_err[8];

  double effi_kp[8];
  double effi_kp_err[8];
  double effi_km[8];
  double effi_km_err[8];

  double diff_kp[8];
  double diff_kp_err[8];
  double diff_km[8];
  double diff_km_err[8];

  double raw_diff_kp[8];
  double raw_diff_kp_err[8];
  double raw_diff_km[8];
  double raw_diff_km_err[8];

  double raw_chisq = 0;
  double chisq = 0;
  cout << "Corrected K[i] from B+ -> D0 pi+:" << endl;
  for(int i=0; i<8; i++){
    raw_diff_kp[i] = bdpi_kp[i] - cleo_kp[i];
    raw_diff_kp_err[i] = sqrt(bdpi_kp_err[i]*bdpi_kp_err[i] + cleo_kp_err[i]*cleo_kp_err[i]);
    raw_diff_km[i] = bdpi_km[i] - cleo_km[i];
    raw_diff_km_err[i] = sqrt(bdpi_km_err[i]*bdpi_km_err[i] + cleo_km_err[i]*cleo_km_err[i]);
    raw_chisq += (raw_diff_kp[i]*raw_diff_kp[i])/(raw_diff_kp_err[i]*raw_diff_kp_err[i]);
    raw_chisq += (raw_diff_km[i]*raw_diff_km[i])/(raw_diff_km_err[i]*raw_diff_km_err[i]);

    avmc_kp[i]     = 0.5*(pi0_kp[i]+omega_kp[i]);
    avmc_kp_err[i] = 1./sqrt(1./pi0_kp_err[i]/pi0_kp_err[i]+1./omega_kp_err[i]/omega_kp_err[i]);
    avmc_km[i]     = 0.5*(pi0_km[i]+omega_km[i]);
    avmc_km_err[i] = 1./sqrt(1./pi0_km_err[i]/pi0_km_err[i]+1./omega_km_err[i]/omega_km_err[i]);

    effi_kp[i]     = avmc_kp[i]/model_kp[i];
    effi_kp_err[i] = avmc_kp_err[i]/model_kp[i];
    effi_km[i]     = avmc_km[i]/model_km[i];
    effi_km_err[i] = avmc_km_err[i]/model_km[i];

    diff_kp[i]     = bdpi_kp[i]*effi_kp[i] - cleo_kp[i];
    const double s1 = bdpi_kp[i]*effi_kp_err[i];
    const double s2 = bdpi_kp_err[i]*effi_kp[i];
    diff_kp_err[i] = sqrt(s1*s1+s2*s2+cleo_kp_err[i]*cleo_kp_err[i]);

    diff_km[i]     = bdpi_km[i]*effi_km[i] - cleo_km[i];
    const double s3 = bdpi_km[i]*effi_km_err[i];
    const double s4 = bdpi_km_err[i]*effi_km[i];
    diff_km_err[i] = sqrt(s3*s3+s4*s4+cleo_km_err[i]*cleo_km_err[i]);
    chisq += (diff_km[i]*diff_km[i])/(diff_km_err[i]*diff_km_err[i]);

    cout << "bin " << i+1;
    cout << ", K = "  << bdpi_kp[i]*effi_kp[i] << " +- " << sqrt(s1*s1+s2*s2);
    cout << ", Kb = " << bdpi_km[i]*effi_km[i] << " +- " << sqrt(s3*s3+s4*s4);
    cout << endl;
  }
  raw_chisq /= 16;
  chisq /= 16;
  stringstream out;

  TCanvas* c1 = new TCanvas("c1","c1",400,400);
  c1->cd();
  c1->SetGrid();
  TGraphErrors* effi_gr_p = new TGraphErrors(nbins,bins,effi_kp,bins_err,effi_kp_err);
  effi_gr_p->SetMarkerStyle(20);
  effi_gr_p->SetMarkerColor(kBlue);
  effi_gr_p->SetMarkerSize(1.3);
  TGraphErrors* effi_gr_m = new TGraphErrors(nbins,bins,effi_km,bins_err,effi_km_err);
  effi_gr_m->SetMarkerStyle(20);
  effi_gr_m->SetMarkerColor(kRed);
  effi_gr_m->SetMarkerSize(1.3);
  TMultiGraph* effi_mg = new TMultiGraph("effi_mg","Relative detection efficiency");
  effi_mg->Add(effi_gr_p);
  effi_mg->Add(effi_gr_m);
  effi_mg->Draw("ap");
  effi_mg->GetXaxis()->SetTitle("Dalitz bin");
  effi_mg->GetXaxis()->SetTitleSize(0.06);
  effi_mg->GetXaxis()->SetTitleOffset(0.75);
  effi_mg->GetXaxis()->SetLabelSize(0.06);
  effi_mg->GetYaxis()->SetRangeUser(0.85,1.15);
  effi_mg->GetYaxis()->SetLabelSize(0.05);
  effi_mg->Draw("ap");
  c1->Print("pics/k_effi.eps");

  TCanvas* c2 = new TCanvas("c2","c2",400,400);
  c2->cd();
  c2->SetGrid();
  TGraphErrors* diff_gr_p = new TGraphErrors(nbins,bins,diff_kp,bins_err,diff_kp_err);
  diff_gr_p->SetMarkerStyle(20);
  diff_gr_p->SetMarkerColor(kBlue);
  diff_gr_p->SetMarkerSize(1.3);
  TGraphErrors* diff_gr_m = new TGraphErrors(nbins,bins,diff_km,bins_err,diff_km_err);
  diff_gr_m->SetMarkerStyle(20);
  diff_gr_m->SetMarkerColor(kRed);
  diff_gr_m->SetMarkerSize(1.3);
  TMultiGraph* diff_mg = new TMultiGraph("diff_mg","K_{i} difference (10^{-2}) with CLEO measurement");
  diff_mg->Add(diff_gr_p);
  diff_mg->Add(diff_gr_m);
  diff_mg->Draw("ap");
  diff_mg->GetXaxis()->SetTitle("Dalitz bin");
  diff_mg->GetXaxis()->SetTitleSize(0.06);
  diff_mg->GetXaxis()->SetTitleOffset(0.75);
  diff_mg->GetXaxis()->SetLabelSize(0.06);
  diff_mg->GetYaxis()->SetRangeUser(-2,2);
  diff_mg->GetYaxis()->SetLabelSize(0.05);
  diff_mg->Draw("ap");

  TPaveText *pt1 = new TPaveText(0.25,0.18,0.75,0.30,"brNDC");
  pt1->SetFillColor(0);
  pt1->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt1->AddText(out.str().c_str());
  pt1->Draw();
  c2->Print("pics/k_diff.eps");

  TCanvas* c3 = new TCanvas("c3","c3",400,400);
  c3->cd();
  c3->SetGrid();
  TGraphErrors* raw_diff_gr_p = new TGraphErrors(nbins,bins,raw_diff_kp,bins_err,raw_diff_kp_err);
  raw_diff_gr_p->SetMarkerStyle(20);
  raw_diff_gr_p->SetMarkerColor(kBlue);
  raw_diff_gr_p->SetMarkerSize(1.3);
  TGraphErrors* raw_diff_gr_m = new TGraphErrors(nbins,bins,raw_diff_km,bins_err,raw_diff_km_err);
  raw_diff_gr_m->SetMarkerStyle(20);
  raw_diff_gr_m->SetMarkerColor(kRed);
  raw_diff_gr_m->SetMarkerSize(1.3);
  TMultiGraph* raw_diff_mg = new TMultiGraph("raw_diff_mg","K_{i} raw difference (10^{-2}) with CLEO measurement");
  raw_diff_mg->Add(raw_diff_gr_p);
  raw_diff_mg->Add(raw_diff_gr_m);
  raw_diff_mg->Draw("ap");
  raw_diff_mg->GetXaxis()->SetTitle("Dalitz bin");
  raw_diff_mg->GetXaxis()->SetTitleSize(0.06);
  raw_diff_mg->GetXaxis()->SetTitleOffset(0.75);
  raw_diff_mg->GetXaxis()->SetLabelSize(0.06);
  raw_diff_mg->GetYaxis()->SetRangeUser(-2,2);
  raw_diff_mg->GetYaxis()->SetLabelSize(0.05);
  raw_diff_mg->Draw("ap");

  TPaveText *pt2 = new TPaveText(0.25,0.18,0.76,0.30,"brNDC");
  pt2->SetFillColor(0);
  pt2->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << raw_chisq;
  pt2->AddText(out.str().c_str());
  pt2->Draw();
  c3->Print("pics/k_raw_diff.eps");

  return;
}
