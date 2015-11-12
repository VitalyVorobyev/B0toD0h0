void draw_smth(const double* vals,const double* errs, const int issin, const string& title, const string& lbl, const double& range){
  string streams[7] = {"S1","S2","S3","S4","S5","S6","All"};
  double strs[7]    = {1,2,3,4,5,6,7};
  double strerrs[7] = {0,0,0,0,0,0,0};

  double true_value = 0;
  if(issin == 1)      true_value = 0.719340;
  else if(issin == 0) true_value = 0.694658;

  double chisq1 = 0;
  const double sigmas = (vals[6]-true_value)/errs[6];
  const double chisq2 = sigmas*sigmas;//(vals[6]-true_value)*(vals[6]-true_value)/(errs[6]*errs[6]);
  for(int i=0; i<6; i++) chisq1 += (vals[i]-true_value)*(vals[i]-true_value)/(errs[i]*errs[i]);
//  for(int i=0; i<6; i++) chisq1 += (vals[i]-vals[6])*(vals[i]-vals[6])/(errs[i]*errs[i]);
  chisq1 /= 6.;

  TGraphErrors* gr = new TGraphErrors(7,strs,vals,strerrs,errs);

  TAxis *ay = gr->GetHistogram()->GetXaxis();
  Double_t y1 = ay->GetBinLowEdge(1);
  Double_t y2 = ay->GetBinUpEdge(ay->GetNbins()+1);
  gr->GetHistogram()->GetXaxis()->Set(7,y1,y2);
  for(Int_t k=0;k<7;k++) gr->GetHistogram()->GetXaxis()->SetBinLabel(k+1,streams[k].c_str());

  stringstream out;
  out.str("");
  out << "c_" << lbl;

  TCanvas* c1 = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
  c1->Draw();
  c1->cd();
  c1->SetGridy();
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerSize(1.3);
  gr->GetXaxis()->SetLabelSize(0.1);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetRangeUser(-range,range);
  gr->SetTitle(title.c_str());

  TMarker* mrk = new TMarker(7,vals[6],20);
  mrk->SetMarkerColor(kBlue);
  mrk->SetMarkerSize(1.5);

  TLine* line = new TLine(0.5,true_value,7.5,true_value);
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);

  TPaveText *pt = new TPaveText(0.60,0.13,0.9,0.33,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << std::setprecision(3);
  out << "#chi_{1}^{2}    = " << chisq2 << " (" << sigmas << "#sigma)";
  pt->AddText(out.str().c_str());
  out.str("");
  out << "#chi_{2}^{2}/6 = " << chisq1;
  pt->AddText(out.str().c_str());

  gr->Draw("AP");
  mrk->Draw();
  line->Draw();
  pt->Draw();

  out.str("");
  out << "pics/cpv_test_" << lbl << ".eps";
  c1->Print(out.str().c_str());

  return;
}

void draw_cpv_test(void){
    // pi0 //
  const double vals_sin_pi0_nb[7] = {0.509722,0.175482,0.513919,0.878458,-0.110159,0.503283,0.413961};
  const double errs_sin_pi0_nb[7] = {0.376131,0.326300,0.365814,0.309075, 0.369529,0.280405,0.136745};
  const double vals_cos_pi0_nb[7] = {0.668969,0.880117,0.050363,0.341849, 0.201591,1.452640,0.371254};
  const double errs_cos_pi0_nb[7] = {0.478672,0.452541,0.455987,0.460784, 0.473227,0.417083,0.192934};

//  draw_smth(vals_sin_pi0_nb,errs_sin_pi0_nb,1, string("Sin #pi^{0} (no Bkg)"),string("sin_pi0_nb"));
//  draw_smth(vals_cos_pi0_nb,errs_cos_pi0_nb,0,string("Cos #pi^{0} (no Bkg)"),string("cos_pi0_nb"));

  const double vals_sin_pi0[7] = {0.348882,-0.057525,0.712821,1.067600,0.148119,0.486814,0.418961};
  const double errs_sin_pi0[7] = {0.408618, 0.343347,0.390015,0.325847,0.363847,0.303462,0.142973};
  const double vals_cos_pi0[7] = {0.974837, 0.707821,0.235489,0.067891,0.384872,1.461050,0.347320};
  const double errs_cos_pi0[7] = {0.496554, 0.484875,0.488484,0.515366,0.517676,0.476137,0.207910};

//  draw_smth(vals_sin_pi0,errs_sin_pi0,1, string("Sin #pi^{0}"),string("sin_pi0"));
//  draw_smth(vals_cos_pi0,errs_cos_pi0,0,string("Cos #pi^{0}"),string("cos_pi0"));

  // omega //
  const double vals_sin_omega_nb[7] = {0.484669,0.940462,0.306539,0.352160,0.349865,0.780129,0.812037};
  const double errs_sin_omega_nb[7] = {0.388263,0.444313,0.415532,0.423731,0.435445,0.419841,0.165324};
  const double vals_cos_omega_nb[7] = {0.396300,0.899557,0.164807,2.104420,0.803494,2.517560,0.328612};
  const double errs_cos_omega_nb[7] = {0.586389,0.727503,0.741789,0.665011,0.614456,0.687360,0.271387};

//  draw_smth(vals_sin_omega_nb,errs_sin_omega_nb,1, string("Sin #omega (no Bkg)"),string("sin_omega_nb"));
//  draw_smth(vals_cos_omega_nb,errs_cos_omega_nb,0,string("Cos #omega (no Bkg)"),string("cos_omega_nb"));

  const double vals_sin_omega[7] = {0.487475,1.099620,0.127523,0.295108,0.221554,1.146650,0.872192};
  const double errs_sin_omega[7] = {0.421214,0.445156,0.458955,0.485310,0.480782,0.463100,0.191092};
  const double vals_cos_omega[7] = {0.290832,0.610225,0.162251,1.741410,0.801934,3.053750,0.198699};
  const double errs_cos_omega[7] = {0.644925,0.690908,0.786787,0.770146,0.671701,0.750426,0.293389};

//  draw_smth(vals_sin_omega,errs_sin_omega,1, string("Sin #omega"),string("sin_omega"));
//  draw_smth(vals_cos_omega,errs_cos_omega,0,string("Cos #omega"),string("cos_omega"));

  // pi0 eta omega //
//  const double vals_sin_full_nb[7] = {0.388864,0.591421,0.306592,0.657689,0.177763,0.517502,0.535898};
//  const double errs_sin_full_nb[7] = {0.237577,0.241454,0.226089,0.224243,0.254138,0.212273,0.094512};
//  const double vals_cos_full_nb[7] = {0.556817,0.723564,0.052946,0.957075,0.622058,1.745680,0.338003};
//  const double errs_cos_full_nb[7] = {0.307819,0.353687,0.342353,0.308988,0.342741,0.309331,0.138913};
  const double vals_sin_full_nb[7] = {1.001110,0.870265,0.890364,0.777768,0.811201,0.933496,0.817792};
  const double errs_sin_full_nb[7] = {0.235247,0.228112,0.235948,0.199889,0.221819,0.215416,0.091536};
  const double vals_cos_full_nb[7] = {0.822665,1.091760,0.534083,1.000000,0.978988,0.944387,0.896759};
  const double errs_cos_full_nb[7] = {0.281504,0.302339,0.323803,0.325072,0.327997,0.309553,0.124764};

//  draw_smth(vals_sin_full_nb,errs_sin_full_nb,1, string("Sin (no Bkg)"),string("sin_full_nb"));
//  draw_smth(vals_cos_full_nb,errs_cos_full_nb,0,string("Cos (no Bkg)"),string("cos_full_nb"));

//  const double vals_sin_full[7] = {0.411227,0.455597,0.293569,0.758373,0.049330,0.570909,0.520552};
//  const double errs_sin_full[7] = {0.255071,0.258172,0.250395,0.240265,0.266640,0.231753,0.102742};
//  const double vals_cos_full[7] = {0.567014,0.689669,0.168591,0.556845,0.704993,1.806570,0.281865};
//  const double errs_cos_full[7] = {0.344143,0.368504,0.369645,0.385445,0.375606,0.344406,0.151970};
  const double vals_sin_full[7] = {1.056660,0.758016,0.835610,0.846839,0.684047,0.948452,0.800653};
  const double errs_sin_full[7] = {0.259924,0.250115,0.254615,0.212623,0.227983,0.230069,0.098040};
  const double vals_cos_full[7] = {0.837837,1.041740,0.640939,0.652516,1.037320,1.014400,0.864006};
  const double errs_cos_full[7] = {0.322220,0.326810,0.346996,0.389921,0.366896,0.331662,0.138871};

  double bkg_effect_val_cos[7];
  double bkg_effect_err_cos[7];
  double bkg_effect_val_sin[7];
  double bkg_effect_err_sin[7];
  for(int i=0; i<7; i++){
    bkg_effect_val_sin[i] = vals_sin_full[i]-vals_sin_full_nb[i];
    bkg_effect_err_sin[i] = sqrt(errs_sin_full[i]*errs_sin_full[i]-errs_sin_full_nb[i]*errs_sin_full_nb[i]);
    bkg_effect_val_cos[i] = vals_cos_full[i]-vals_cos_full_nb[i];
    bkg_effect_err_cos[i] = sqrt(errs_cos_full[i]*errs_cos_full[i]-errs_cos_full_nb[i]*errs_cos_full_nb[i]);
  }
  draw_smth(bkg_effect_val_sin,bkg_effect_err_sin,2, string("Sin (background effect)"),string("sin_full_bkgeff"),0.5);
  draw_smth(bkg_effect_val_cos,bkg_effect_err_cos,2, string("Cos (background effect)"),string("cos_full_bkgeff"),0.5);

//  const double vals_sin_full_nuis_fbb   = 1.115510;
//  const double errs_sin_full_nuis_fbb   = 0.270419;
//  const double vals_cos_full_nuis_fbb   = 0.883140;
//  const double errs_cos_full_nuis_fbb   = 0.338869;

//  const double vals_sin_full_nuis_nsig  = 1.06748;
//  const double errs_sin_full_nuis_nsig  = 0.260928;
//  const double vals_cos_full_nuis_nsig  = 0.848884;
//  const double errs_cos_full_nuis_nsig  = 0.323833;

  const double vals_sin_full_nuis_ki  = 1.05648;
  const double errs_sin_full_nuis_ki  = 0.259957;
  const double vals_cos_full_nuis_ki  = 0.837352;
  const double errs_cos_full_nuis_ki  = 0.322435;

  const double vals_sin_full_nuis_csi = 1.08228;
  const double errs_sin_full_nuis_csi = 0.267542;
  const double vals_cos_full_nuis_csi = 0.940412;
  const double errs_cos_full_nuis_csi = 0.349435;

  const double vals_sin_full_nuis_tau = 1.0562;
  const double errs_sin_full_nuis_tau = 0.259533;
  const double vals_cos_full_nuis_tau = 0.835765;
  const double errs_cos_full_nuis_tau = 0.322043;

  const double vals_sin_full_nuis_tau_s1 = 0.759596;
  const double errs_sin_full_nuis_tau_s1 = 0.25061;
  const double vals_cos_full_nuis_tau_s1 = 1.04317;
  const double errs_cos_full_nuis_tau_s1 = 0.327018;

//  const double vals_sin_full_nuis_scale[3] = {1.058970,1.107590,1.07206};
//  const double errs_sin_full_nuis_scale[3] = {0.260193,0.266732,0.26180};
//  const double vals_cos_full_nuis_scale[3] = {0.838717,0.867927,0.84582};
//  const double errs_cos_full_nuis_scale[3] = {0.322405,0.326306,0.32334};

  draw_smth(vals_sin_full,errs_sin_full,1, string("Sin"),string("sin_full"),1.5);
  draw_smth(vals_cos_full,errs_cos_full,0,string("Cos"),string("cos_full"),1.5);

  return;
}
