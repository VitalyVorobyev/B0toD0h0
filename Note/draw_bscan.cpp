void draw_bscan(void){
  gStyle->SetOptFit(111);
  const double beta = 23.;
  const double true_sin = TMath::Sin(2.*beta/180.*TMath::Pi());
  const double true_cos = TMath::Cos(2.*beta/180.*TMath::Pi());

  const double sin_pi_raw[8]     = {5.37601,8.98176,8.1477,3.20178,4.95142,8.80586,-571.93,5.70282};
  const double sin_err_pi_raw[8] = {1.34583,3.30593,53.346,2.87926,1.30572,3.78047,729.357,3.14335};
  const double cos_pi_raw[8]     = {135.153,0.760521,2.66545,1.92436,-102.814,3.87988,8.63437,15.982};
  const double cos_err_pi_raw[8] = {165.604,3.45427,2.70724,3.5995,69.8559,3.47733,1.94763,4.3117};

  const double sin_omega_raw[8]     = {4.16782,2.98862,-116.055,-2.50392,1.63292,1.5513,428.066,9.56341};
  const double sin_err_omega_raw[8] = {1.55584,3.82523,62.6632,3.27219,1.48948,4.33989,653.912,3.60982};
  const double cos_omega_raw[8]     = {-225.285,-4.14024,0.674166,7.48107,-67.9957,1.02806,9.35539,3.9276};
  const double cos_err_omega_raw[8] = {190.382,3.99682,3.18223,4.09087,79.6059,3.99192,2.14351,4.95148};

  const double sin_pi[8]     = {5.02019,8.63824,7.00309,3.42892,5.19161,7.75691,-571.933,4.3004};
  const double sin_err_pi[8] = {1.3399,3.29287,52.4259,2.8871,1.3084,3.73242,731.448,3.07795};
  const double cos_pi[8]     = {137.66,0.462715,1.25226,2.15617,-103.01,2.77138,7.63314,14.2144};
  const double cos_err_pi[8] = {164.892,3.44062,2.66039,3.60931,69.9983,3.43312,1.92601,4.22202};

  const double sin_omega[8]     = {4.02945,2.45437,-116.547,-2.06525,1.53961,0.864418,428.066,8.50429};
  const double sin_err_omega[8] = {1.55337,3.80229,63.6026,3.2875,1.48811,4.31301,657.13,3.56651};
  const double cos_omega[8]     = {-225.762,-4.49556,1.81339,7.95109,-68.0433,0.497502,8.92347,2.94688};
  const double cos_err_omega[8] = {190.093,3.97285,3.23019,4.11002,79.5334,3.9672,2.1338,4.89207};

  const double bins[8]     = {1,2,3,4,5,6,7,8};
  const double bins_err[8] = {0,0,0,0,0,0,0,0};

  const double min_offset = -20.;
  const double max_offset =  20.;

  // sin pi0 //
  TCanvas *c_sin_bins_scan = new TCanvas("c_sin_bins_scan","c_sin_bins_scan",400,400);
  c_sin_bins_scan->SetGrid();
  c_sin_bins_scan->Draw();
  TGraphErrors* gr_sin_pi = new TGraphErrors(8,bins,sin_pi,bins_err,sin_err_pi);
  gr_sin_pi->SetMarkerStyle(20);
  gr_sin_pi->SetMarkerSize(1.3);
  gr_sin_pi->SetMarkerColor(kBlue);

  TGraphErrors* gr_sin_raw_pi = new TGraphErrors(8,bins,sin_pi_raw,bins_err,sin_err_pi_raw);
  gr_sin_raw_pi->SetMarkerStyle(21);
  gr_sin_raw_pi->SetMarkerSize(1.3);
  gr_sin_raw_pi->SetMarkerColor(kRed);

  TMultiGraph* mg_sin_pi = new TMultiGraph("mg_sin_pi","sin offset (10^{-2}), #pi^{0}");
  mg_sin_pi->Add(gr_sin_raw_pi);
  mg_sin_pi->Add(gr_sin_pi);
  mg_sin_pi->Draw("ap");

  mg_sin_pi->GetYaxis()->SetRangeUser(min_offset,max_offset);
  mg_sin_pi->GetYaxis()->SetLabelSize(0.06);
  mg_sin_pi->GetXaxis()->SetLabelSize(0.06);
  mg_sin_pi->GetXaxis()->SetTitle("Dalitz bin");
  mg_sin_pi->GetXaxis()->SetTitleSize(0.06);
  mg_sin_pi->GetXaxis()->SetTitleOffset(0.8);

  TFitResultPtr fit_res_pt_sin_pi = gr_sin_pi->Fit("pol0","s");
  TFitResult* fit_res_sin_pi = fit_res_pt_sin_pi.Get();
//  fit_res_sin_pi.Print();

  // cos pi0 //
  TCanvas *c_cos_bins_scan = new TCanvas("c_cos_bins_scan","c_cos_bins_scan",400,400);
  c_cos_bins_scan->SetGrid();
  c_cos_bins_scan->Draw();
  TGraphErrors* gr_cos_pi = new TGraphErrors(8,bins,cos_pi,bins_err,cos_err_pi);
  gr_cos_pi->SetMarkerStyle(20);
  gr_cos_pi->SetMarkerSize(1.3);
  gr_cos_pi->SetMarkerColor(kBlue);

  TGraphErrors* gr_cos_raw_pi = new TGraphErrors(8,bins,cos_pi_raw,bins_err,cos_err_pi_raw);
  gr_cos_raw_pi->SetMarkerStyle(21);
  gr_cos_raw_pi->SetMarkerSize(1.3);
  gr_cos_raw_pi->SetMarkerColor(kRed);

  TMultiGraph* mg_cos_pi = new TMultiGraph("mg_cos_pi","cos offset (10^{-2}), #pi^{0}");
  mg_cos_pi->Add(gr_cos_raw_pi);
  mg_cos_pi->Add(gr_cos_pi);
  mg_cos_pi->Draw("ap");

  mg_cos_pi->GetYaxis()->SetRangeUser(min_offset,max_offset);
  mg_cos_pi->GetYaxis()->SetLabelSize(0.06);
  mg_cos_pi->GetXaxis()->SetLabelSize(0.06);
  mg_cos_pi->GetXaxis()->SetTitle("Dalitz bin");
  mg_cos_pi->GetXaxis()->SetTitleSize(0.06);
  mg_cos_pi->GetXaxis()->SetTitleOffset(0.8);

  TFitResultPtr fit_res_pt_cos_pi = gr_cos_pi->Fit("pol0","s");
  TFitResult* fit_res_cos_pi = fit_res_pt_sin_pi.Get();
//  fit_res_sin_pi.Print();

  // sin omega //
  TCanvas *c_sin_bins_scan = new TCanvas("c_sin_bins_scan_omega","c_sin_bins_scan_omega",400,400);
  c_sin_bins_scan->SetGrid();
  c_sin_bins_scan->Draw();
  TGraphErrors* gr_sin_omega = new TGraphErrors(8,bins,sin_omega,bins_err,sin_err_omega);
  gr_sin_omega->SetMarkerStyle(20);
  gr_sin_omega->SetMarkerSize(1.3);
  gr_sin_omega->SetMarkerColor(kBlue);

  TGraphErrors* gr_sin_raw_omega = new TGraphErrors(8,bins,sin_omega_raw,bins_err,sin_err_omega_raw);
  gr_sin_raw_omega->SetMarkerStyle(21);
  gr_sin_raw_omega->SetMarkerSize(1.3);
  gr_sin_raw_omega->SetMarkerColor(kRed);

  TMultiGraph* mg_sin_omega = new TMultiGraph("mg_sin_omega","sin offset (10^{-2}), #omega");
  mg_sin_omega->Add(gr_sin_raw_omega);
  mg_sin_omega->Add(gr_sin_omega);
  mg_sin_omega->Draw("ap");

  mg_sin_omega->GetYaxis()->SetRangeUser(min_offset,max_offset);
  mg_sin_omega->GetYaxis()->SetLabelSize(0.06);
  mg_sin_omega->GetXaxis()->SetLabelSize(0.06);
  mg_sin_omega->GetXaxis()->SetTitle("Dalitz bin");
  mg_sin_omega->GetXaxis()->SetTitleSize(0.06);
  mg_sin_omega->GetXaxis()->SetTitleOffset(0.8);

  TFitResultPtr fit_res_pt_sin_omega = gr_sin_omega->Fit("pol0","s");
  TFitResult* fit_res_sin_omega = fit_res_pt_sin_omega.Get();
//  fit_res_sin_omega.Print();

  // cos omega //
  TCanvas *c_cos_bins_scan = new TCanvas("c_cos_bins_scan_omega","c_cos_bins_scan_omega",400,400);
  c_cos_bins_scan->SetGrid();
  c_cos_bins_scan->Draw();
  TGraphErrors* gr_cos_omega = new TGraphErrors(8,bins,cos_omega,bins_err,cos_err_omega);
  gr_cos_omega->SetMarkerStyle(20);
  gr_cos_omega->SetMarkerSize(1.3);
  gr_cos_omega->SetMarkerColor(kBlue);

  TGraphErrors* gr_cos_raw_omega = new TGraphErrors(8,bins,cos_omega_raw,bins_err,cos_err_omega_raw);
  gr_cos_raw_omega->SetMarkerStyle(21);
  gr_cos_raw_omega->SetMarkerSize(1.3);
  gr_cos_raw_omega->SetMarkerColor(kRed);

  TMultiGraph* mg_cos_omega = new TMultiGraph("mg_cos_omega","cos offset (10^{-2}), #omega");
  mg_cos_omega->Add(gr_cos_raw_omega);
  mg_cos_omega->Add(gr_cos_omega);
  mg_cos_omega->Draw("ap");

  mg_cos_omega->GetYaxis()->SetRangeUser(min_offset,max_offset);
  mg_cos_omega->GetYaxis()->SetLabelSize(0.06);
  mg_cos_omega->GetXaxis()->SetLabelSize(0.06);
  mg_cos_omega->GetXaxis()->SetTitle("Dalitz bin");
  mg_cos_omega->GetXaxis()->SetTitleSize(0.06);
  mg_cos_omega->GetXaxis()->SetTitleOffset(0.8);

  TFitResultPtr fit_res_pt_cos_omega = gr_cos_omega->Fit("pol0","s");
  TFitResult* fit_res_cos_omega = fit_res_pt_sin_omega.Get();
//  fit_res_sin_omega.Print();

}

