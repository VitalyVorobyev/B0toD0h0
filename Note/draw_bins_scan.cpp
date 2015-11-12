void draw_bins_scan(void){
  // * Sampled offsets * //
  // sin cos sampled pi0
  const int NBins = 8;
//  char* FourModes[FourBins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega"};

  const double bins[NBins]        = {1,2,3,4,5,6,7,8};
  const double bins_err[NBins]    = {0,0,0,0,0,0,0,0};
  const double sin_pi0[NBins]     = {0.162463,-1.47185,0.107467,0.035329, 0.02022,-0.481105, 0.102757, 0.044575};
  const double sin_pi0_err[NBins] = {0.046884, 3.87315,0.057930,0.019464, 0.04978, 0.366296, 0.043770, 0.019559};
  const double cos_pi0[NBins]     = {0.146995, 0.04017,0.027093,2.305342,-0.05958,-0.024925,-0.011766,-0.233421};
  const double cos_pi0_err[NBins] = {0.058180, 0.02814,0.051949,3.3999,   0.04921, 0.040011, 0.056467, 0.477747};
  const double tau_pi0[NBins]     = {-0.0436,-0.02092,};
  const double tau_pi0_err[NBins] = { 0.0122, 0.01088,};

  TCanvas *c_sin_bins_scan = new TCanvas("c_sin_bins_scan","c_sin_bins_scan",400,400);
  c_sin_bins_scan->SetGrid();
  c_sin_bins_scan->Draw();
  TGraphErrors* gr_sin_pi0 = new TGraphErrors(NBins,bins,sin_pi0,bins_err,sin_pi0_err);
  gr_sin_pi0->SetTitle("sin offset");
  gr_sin_pi0->SetMarkerStyle(20);
  gr_sin_pi0->SetMarkerSize(1.);
  gr_sin_pi0->SetMarkerColor(kBlue);
  gr_sin_pi0->GetYaxis()->SetRangeUser(-0.2,0.2);
  gr_sin_pi0->GetYaxis()->SetLabelSize(0.06);
  gr_sin_pi0->GetXaxis()->SetLabelSize(0.06);
  gr_sin_pi0->Draw("ap");
  gr_sin_pi0->Fit("pol0");
  c_sin_bins_scan->Update();

  TCanvas *c_cos_bins_scan = new TCanvas("c_cos_bins_scan","c_cos_bins_scan",400,400);
  c_cos_bins_scan->SetGrid();
  c_cos_bins_scan->Draw();
  TGraphErrors* gr_cos_pi0 = new TGraphErrors(NBins,bins,cos_pi0,bins_err,cos_pi0_err);
  gr_cos_pi0->SetTitle("cos offset");
  gr_cos_pi0->SetMarkerStyle(20);
  gr_cos_pi0->SetMarkerSize(1.);
  gr_cos_pi0->SetMarkerColor(kBlue);
  gr_cos_pi0->GetYaxis()->SetRangeUser(-0.2,0.2);
  gr_cos_pi0->GetYaxis()->SetLabelSize(0.06);
  gr_cos_pi0->GetXaxis()->SetLabelSize(0.06);
  gr_cos_pi0->Draw("ap");
  gr_cos_pi0->Fit("pol0");
  c_cos_bins_scan->Update();
  return;
}
