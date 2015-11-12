void draw_each_mode_offset(void){
  const int Nbins = 4;
  const double offmin = -10;
  const double offmax =  10;

  string  cpvModes[Nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega"};
  string h0Modes[2] = {"Single","Multiple"};

//  double vals_sin_nobkg[Nbins] = { 0.6189,0.232785,0.617006, 0.6027,1.18236};
//  double errs_sin_nobkg[Nbins] = { 0.1548,0.268557,0.358842, 0.1882,0.88690};
//  double vals_cos_nobkg[Nbins] = {-0.1332,0.836078,0.505995,-0.2609,1.81729};
//  double errs_cos_nobkg[Nbins] = { 0.2275,0.395649,0.600699, 0.2944,0.99490};

//  double vals_sin[Nbins] =       { 0.595674,0.21658,0.407979, 0.59394,1.1692};
//  double errs_sin[Nbins] =       { 0.164558,0.30924,0.449164, 0.20332,1.0143};
//  double vals_cos[Nbins] =       {-0.181697,1.01185,0.718472,-0.40388,0.6837};
//  double errs_cos[Nbins] =       { 0.250834,0.48310,0.666282, 0.32315,0.9840};

  double vals_sin[Nbins] = {0.820785,0.941559,0.267048,0.808787};
  double errs_sin[Nbins] = {0.135001,0.284234,0.447967,0.181179};
  double vals_cos[Nbins] = {0.695985,1.275680,0.942595,0.959159};
  double errs_cos[Nbins] = {0.190945,0.387646,0.541324,0.257149};

  double vals_sin_nobkg[Nbins] = {0.859655,0.871896,0.555907,0.763584};
  double errs_sin_nobkg[Nbins] = {0.131422,0.233725,0.389491,0.168293};
  double vals_cos_nobkg[Nbins] = {0.737882,1.350460,0.730613,0.974639};
  double errs_cos_nobkg[Nbins] = {0.176188,0.306705,0.502757,0.230145};

  double vals_sin_diff[Nbins];
  double errs_sin_diff[Nbins];
  double vals_cos_diff[Nbins];
  double errs_cos_diff[Nbins];

  for(int i=0; i<Nbins; i++){
    vals_sin_diff[i] = (vals_sin[i] - vals_sin_nobkg[i])*10.;
    vals_cos_diff[i] = (vals_cos[i] - vals_cos_nobkg[i])*10.;
    errs_sin_diff[i] = 10.*sqrt(errs_sin[i]*errs_sin[i]-errs_sin_nobkg[i]*errs_sin_nobkg[i]);
    errs_cos_diff[i] = 10.*sqrt(errs_cos[i]*errs_cos[i]-errs_cos_nobkg[i]*errs_cos_nobkg[i]);
  }

//  double vals_sin_h0_nobkg[2] = {0.531301,0.658114};
//  double errs_sin_h0_nobkg[2] = {0.132553,0.162885};
//  double vals_cos_h0_nobkg[2] = {0.283538,0.478752};
//  double errs_cos_h0_nobkg[2] = {0.198171,0.254984};

//  double vals_sin_h0[2] = {0.515947,0.592382};
//  double errs_sin_h0[2] = {0.145799,0.182117};
//  double vals_cos_h0[2] = {0.348879,0.561866};
//  double errs_cos_h0[2] = {0.222419,0.278932};

  double vals_sin_h0_nobkg[2] = {0.870347,0.725135};
  double errs_sin_h0_nobkg[2] = {0.115232,0.153259};
  double vals_cos_h0_nobkg[2] = {0.871090,0.925734};
  double errs_cos_h0_nobkg[2] = {0.152677,0.211000};

  double vals_sin_h0[2] = {0.840389,0.720376};
  double errs_sin_h0[2] = {0.121666,0.165934};
  double vals_cos_h0[2] = {0.801738,0.952603};
  double errs_cos_h0[2] = {0.170384,0.233909};

  double vals_sin_diff_h0[2];
  double errs_sin_diff_h0[2];
  double vals_cos_diff_h0[2];
  double errs_cos_diff_h0[2];

  cout << "Background dt offsets:" << endl;
  for(int i=0; i<2; i++){
    vals_sin_diff_h0[i] = (vals_sin_h0[i] - vals_sin_h0_nobkg[i])*10.;
    vals_cos_diff_h0[i] = (vals_cos_h0[i] - vals_cos_h0_nobkg[i])*10.;
    errs_sin_diff_h0[i] = 10.*sqrt(errs_sin_h0[i]*errs_sin_h0[i]-errs_sin_h0_nobkg[i]*errs_sin_h0_nobkg[i]);
    errs_cos_diff_h0[i] = 10.*sqrt(errs_cos_h0[i]*errs_cos_h0[i]-errs_cos_h0_nobkg[i]*errs_cos_h0_nobkg[i]);
    cout << "  " << vals_sin_diff_h0[i] << " +- " << errs_sin_diff_h0[i];
    cout << ", " << vals_cos_diff_h0[i] << " +- " << errs_cos_diff_h0[i] << endl;
  }

  TCanvas *c1_h0_sin_diff = new TCanvas("c1_h0_sin_diff","c1_h0_sin_diff",400,400);
  c1_h0_sin_diff->SetGrid();
  TH1D *h_sin_h0 = new TH1D("sin_offset_small","Background effect on sin (10^{-1})",3,0,3);
  h_sin_h0->SetStats(0);
  h_sin_h0->SetFillColor(38);
  h_sin_h0->SetMarkerStyle(20);
  h_sin_h0->SetMarkerSize(1.3);
  h_sin_h0->SetMarkerColor(kBlue);
  for(int i=0; i<2; i++){
    h_sin_h0->Fill(h0Modes[i].c_str(),vals_sin_diff_h0[i]);
    h_sin_h0->SetBinError(i+1,errs_sin_diff_h0[i]);
  }
  h_sin_h0->LabelsDeflate();
  h_sin_h0->GetXaxis()->SetLabelSize(0.08);
  h_sin_h0->GetYaxis()->SetRangeUser(-1.,1.);
  h_sin_h0->GetYaxis()->SetLabelSize(0.06);
  h_sin_h0->Draw("e");
  c1_h0_sin_diff->Print("pics/sin_offset_h0.eps");

  TCanvas *c1_h0_cos_diff = new TCanvas("c1_h0_cos_diff","c1_h0_cos_diff",400,400);
  c1_h0_cos_diff->SetGrid();
  TH1D *h_cos_h0 = new TH1D("cos_offset_small","Background effect on cos (10^{-1})",3,0,3);
  h_cos_h0->SetStats(0);
  h_cos_h0->SetFillColor(38);
  h_cos_h0->SetMarkerStyle(20);
  h_cos_h0->SetMarkerSize(1.3);
  h_cos_h0->SetMarkerColor(kBlue);
  for(int i=0; i<2; i++){
    h_cos_h0->Fill(h0Modes[i].c_str(),vals_cos_diff_h0[i]);
    h_cos_h0->SetBinError(i+1,errs_cos_diff_h0[i]);
  }
  h_cos_h0->LabelsDeflate();
  h_cos_h0->GetXaxis()->SetLabelSize(0.08);
  h_cos_h0->GetYaxis()->SetRangeUser(-1.,1.);
  h_cos_h0->GetYaxis()->SetLabelSize(0.06);
  h_cos_h0->Draw("e");
  c1_h0_cos_diff->Print("pics/cos_offset_h0.eps");

  const double true_sin = 0.719340;
  const double true_cos = 0.694658;

  for(int i=0; i<Nbins; i++){
    vals_sin_nobkg[i] = (fabs(vals_sin_nobkg[i]) - true_sin)*10.;
    vals_cos_nobkg[i] = (fabs(vals_cos_nobkg[i]) - true_cos)*10.;
    errs_sin_nobkg[i] *= 10.; errs_cos_nobkg[i] *= 10.;
    cout << vals_sin_nobkg[i] << " " << vals_cos_nobkg[i] << endl;
  }

  for(int i=0; i<Nbins; i++){
    vals_sin[i] = (fabs(vals_sin[i]) - true_sin)*10.;
    vals_cos[i] = (fabs(vals_cos[i]) - true_cos)*10.;
    errs_sin[i] *= 10.; errs_cos[i] *= 10.;
    cout << vals_sin[i] << " " << vals_cos[i] << endl;
  }

  TCanvas *c_sin = new TCanvas("c_sin","c_sin",400,400);
  c_sin->SetGrid();
  TH1D *h_sin = new TH1D("h_sin","sin offset (10^{-1})",3,0,3);
  h_sin->SetStats(0);
  h_sin->SetFillColor(38);
  h_sin->SetMarkerStyle(20);
  h_sin->SetMarkerSize(1.3);
  h_sin->SetMarkerColor(kRed);
  for(int i=0; i<Nbins; i++){
    h_sin->Fill(cpvModes[i].c_str(),vals_sin[i]);
    h_sin->SetBinError(i+1,errs_sin[i]);
  }
  h_sin->LabelsDeflate();
  h_sin->GetXaxis()->SetLabelSize(0.08);
  h_sin->GetYaxis()->SetLabelSize(0.06);
  h_sin->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_sin_nobkg = new TH1D("h_sin_nobkg","sin offset (10^{-1})",3,0,3);
  h_sin_nobkg->SetStats(0);
  h_sin_nobkg->SetFillColor(38);
  h_sin_nobkg->SetMarkerStyle(20);
  h_sin_nobkg->SetMarkerSize(1.3);
  h_sin_nobkg->SetMarkerColor(kBlue);
  for(int i=0; i<Nbins; i++){
    h_sin_nobkg->Fill(cpvModes[i].c_str(),vals_sin_nobkg[i]);
    h_sin_nobkg->SetBinError(i+1,errs_sin_nobkg[i]);
  }

  h_sin->Draw("e");
  h_sin_nobkg->Draw("same");
  c_sin->Print("pics/sin_offset_modes.eps");

  TCanvas *c_sin_diff = new TCanvas("c_sin_diff","c_sin_diff",400,400);
  c_sin_diff->SetGrid();
  TH1D *h_sin_diff = new TH1D("h_sin_diff","Background effect on sin (10^{-1})",3,0,3);
  h_sin_diff->SetStats(0);
  h_sin_diff->SetFillColor(38);
  h_sin_diff->SetMarkerStyle(20);
  h_sin_diff->SetMarkerSize(1.3);
  h_sin_diff->SetMarkerColor(kRed);
  for(int i=0; i<Nbins; i++){
    h_sin_diff->Fill(cpvModes[i].c_str(),vals_sin_diff[i]);
    h_sin_diff->SetBinError(i+1,errs_sin_diff[i]);
  }
  h_sin_diff->LabelsDeflate();
  h_sin_diff->GetXaxis()->SetLabelSize(0.08);
  h_sin_diff->GetYaxis()->SetLabelSize(0.06);
  h_sin_diff->GetYaxis()->SetRangeUser(-5,5);
  h_sin_diff->Draw("e");
  c_sin_diff->Print("pics/sin_modes_diff.eps");

  TCanvas *c_cos_diff = new TCanvas("c_cos_diff","c_cos_diff",400,400);
  c_cos_diff->SetGrid();
  TH1D *h_cos_diff = new TH1D("h_cos_diff","Background effect on cos (10^{-1})",3,0,3);
  h_cos_diff->SetStats(0);
  h_cos_diff->SetFillColor(38);
  h_cos_diff->SetMarkerStyle(20);
  h_cos_diff->SetMarkerSize(1.3);
  h_cos_diff->SetMarkerColor(kRed);
  for(int i=0; i<Nbins; i++){
    h_cos_diff->Fill(cpvModes[i].c_str(),vals_cos_diff[i]);
    h_cos_diff->SetBinError(i+1,errs_cos_diff[i]);
  }
  h_cos_diff->LabelsDeflate();
  h_cos_diff->GetXaxis()->SetLabelSize(0.08);
  h_cos_diff->GetYaxis()->SetLabelSize(0.06);
  h_cos_diff->GetYaxis()->SetRangeUser(-5,5);
  h_cos_diff->Draw("e");
  c_cos_diff->Print("pics/cos_modes_diff.eps");

  TCanvas *c_cos = new TCanvas("c_cos","c_cos",400,400);
  c_cos->SetGrid();
  TH1D *h_cos = new TH1D("h_cos","cos offset (10^{-1})",3,0,3);
  h_cos->SetStats(0);
  h_cos->SetFillColor(38);
  h_cos->SetMarkerStyle(20);
  h_cos->SetMarkerSize(1.3);
  h_cos->SetMarkerColor(kRed);
  for(int i=0; i<Nbins; i++){
    h_cos->Fill(cpvModes[i].c_str(),vals_cos[i]);
    h_cos->SetBinError(i+1,errs_cos[i]);
  }
  h_cos->LabelsDeflate();
  h_cos->GetXaxis()->SetLabelSize(0.08);
  h_cos->GetYaxis()->SetLabelSize(0.06);
  h_cos->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_cos_nobkg = new TH1D("h_cos_nobkg","cos offset (10^{-1})",3,0,3);
  h_cos_nobkg->SetStats(0);
  h_cos_nobkg->SetFillColor(38);
  h_cos_nobkg->SetMarkerStyle(20);
  h_cos_nobkg->SetMarkerSize(1.3);
  h_cos_nobkg->SetMarkerColor(kBlue);
  for(int i=0; i<Nbins; i++){
    h_cos_nobkg->Fill(cpvModes[i].c_str(),vals_cos_nobkg[i]);
    h_cos_nobkg->SetBinError(i+1,errs_cos_nobkg[i]);
  }

  h_cos->Draw("e");
  h_cos_nobkg->Draw("same");

  c_cos->Print("pics/cos_offset_modes.eps");
}
