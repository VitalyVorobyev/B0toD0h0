void draw_gen(void){
  const int NBins = 8;

  const double bins[NBins]        = {1,2,3,4,5,6,7,8};
  const double bins_err[NBins]    = {0,0,0,0,0,0,0,0};

  const double C_true[NBins]     = {0.675798,0.431828,-0.036724,-0.642989,-0.935997,-0.615002,0.00341,0.571212};
  const double S_true[NBins]     = {-0.005376,0.413289,0.725041,0.514323,-0.017438,-0.668626,-0.814657,-0.416431};

  const double dC_pi0[NBins]     = {0.681018,0.447322,-0.0332043,-0.620657,-0.93471,-0.642679,-0.00706839,0.550723};
  const double dC_pi0_err[NBins] = {0.00909975,0.0139587,0.0186512,0.0175696,0.0116614,0.0209614,0.0148829,0.0174118};
  const double dS_pi0[NBins]     = {-0.00833695,0.413185,0.71594,0.527491,-0.0109006,-0.645878,-0.826182,-0.434764};
  const double dS_pi0_err[NBins] = {0.00942306,0.0144547,0.0193135,0.0181942,0.012074,0.0217063,0.0154116,0.0180302};

  const double dC_omega[NBins]     = {0.673161,0.432935,0.00515399,-0.579732,-0.907227,-0.604082,0.0477443,0.602641};
  const double dC_omega_err[NBins] = {0.0112583,0.0173618,0.023557,0.0214486,0.014249,0.0261389,0.0179881,0.0211782};
  const double dS_omega[NBins]     = {0.0119999,0.373264,0.737956,0.538924,-0.0089708,-0.666698,-0.842079,-0.425711};
  const double dS_omega_err[NBins] = {0.0116584,0.0179788,0.0243926,0.0222107,0.0147553,0.0270668,0.0186262,0.0219308};

  for(int i=0; i<8; i++){
    dC_pi0[i] -= C_true[i]; dS_pi0[i] -= S_true[i];
    dC_omega[i] -= C_true[i]; dS_omega[i] -= S_true[i];
  }

  TCanvas *c_scan = new TCanvas("c_scan","c_scan",400,400);
  c_scan->SetGrid();
  c_scan->Draw();
//  TGraphErrors* gr_c_pi0 = new TGraphErrors(NBins,bins,dC_pi0,bins_err,dC_pi0_err);
  TGraphErrors* gr_c_pi0 = new TGraphErrors(NBins,bins,dC_omega,bins_err,dC_omega_err);
  gr_c_pi0->SetTitle("C offset");
  gr_c_pi0->SetMarkerStyle(20);
  gr_c_pi0->SetMarkerSize(1.);
  gr_c_pi0->SetMarkerColor(kBlue);
  gr_c_pi0->GetYaxis()->SetRangeUser(-0.1,0.1);
  gr_c_pi0->GetYaxis()->SetLabelSize(0.06);
  gr_c_pi0->GetXaxis()->SetLabelSize(0.06);
  gr_c_pi0->Draw("ap");
  gr_c_pi0->Fit("pol0");
  c_scan->Update();

  TCanvas *s_scan = new TCanvas("s_scan","s_scan",400,400);
  s_scan->SetGrid();
  s_scan->Draw();
//  TGraphErrors* gr_s_pi0 = new TGraphErrors(NBins,bins,dS_pi0,bins_err,dS_pi0_err);
  TGraphErrors* gr_s_pi0 = new TGraphErrors(NBins,bins,dS_omega,bins_err,dS_omega_err);
  gr_s_pi0->SetTitle("S offset");
  gr_s_pi0->SetMarkerStyle(20);
  gr_s_pi0->SetMarkerSize(1.);
  gr_s_pi0->SetMarkerColor(kBlue);
  gr_s_pi0->GetYaxis()->SetRangeUser(-0.1,0.1);
  gr_s_pi0->GetYaxis()->SetLabelSize(0.06);
  gr_s_pi0->GetXaxis()->SetLabelSize(0.06);
  gr_s_pi0->Draw("ap");
  gr_s_pi0->Fit("pol0");
  s_scan->Update();
  return;
}

