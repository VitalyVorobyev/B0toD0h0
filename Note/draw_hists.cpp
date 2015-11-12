void draw_hists(void){
  // * tau offset total * //
  const int nbins = 7;
  const double offmin = -10;
  const double offmax =  10;
  const int width = 400;
  char* Modes[nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#pi^{+}(D^{0})","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega","#rho","#pi^{+}"};
  double vals[nbins] = {-2.1,-1.3,-0.1,-1.1,-2.8,-3.3,-1.1};
  double errs[nbins] = { 0.4, 0.5, 0.5, 0.7, 0.4, 0.5, 0.5};
  
  TCanvas *c1 = new TCanvas("c1","c1",400,400);
  c1->SetGrid();
  TH1D *h_mean = new TH1D("tau_offset","#tau offset (10^{-2} ps)",3,0,3);
  h_mean->SetStats(0);
  h_mean->SetFillColor(38);
  h_mean->SetMarkerStyle(20);
  h_mean->SetMarkerSize(1.3);
  h_mean->SetMarkerColor(kBlue);
  h_mean->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean->Fill(Modes[i],vals[i]);
    h_mean->SetBinError(i+1,errs[i]);
  }
  h_mean->LabelsDeflate();
  h_mean->GetXaxis()->SetLabelSize(0.08);
  h_mean->GetYaxis()->SetLabelSize(0.06);
  h_mean->Draw("e");
  c1->Print("pics/tau_offset.eps");

  // * tau offset SVD1 single * //
  double vals_svd1_sgl[nbins] = {-1.4,-0.2,3.8,-6.1,-4.5,-6.6,-4.5};
  double errs_svd1_sgl[nbins] = { 2.3, 2.9,3.1, 4.4, 2.6, 3.2, 2.9};

  TCanvas *c1_sgl = new TCanvas("c1_sgl","c1_sgl",width,400);
  c1_sgl->SetGrid();
  TH1D *h_mean_svd1_sgl = new TH1D("tau_offset_svd1_sgl","#tau offset, SVD1 sgl (10^{-2} ps)",3,0,3);
  h_mean_svd1_sgl->SetStats(0);
  h_mean_svd1_sgl->SetFillColor(38);
  h_mean_svd1_sgl->SetMarkerStyle(20);
  h_mean_svd1_sgl->SetMarkerSize(1.3);
  h_mean_svd1_sgl->SetMarkerColor(kBlue);
  h_mean_svd1_sgl->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd1_sgl->Fill(Modes[i],vals_svd1_sgl[i]);
    h_mean_svd1_sgl->SetBinError(i+1,errs_svd1_sgl[i]);
  }
  h_mean_svd1_sgl->LabelsDeflate();
  h_mean_svd1_sgl->GetXaxis()->SetLabelSize(0.08);
  h_mean_svd1_sgl->GetYaxis()->SetLabelSize(0.06);
  h_mean_svd1_sgl->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean_svd1_sgl->Draw("e");
  c1_sgl->Print("pics/tau_offset_svd1_sgl.eps");
  
  // * tau offset SVD1 multiple * //
  double vals_svd1_mlt[nbins] = {4.8,7.5,5.6,2.0,1.1,-4.1,-1.7};
  double errs_svd1_mlt[nbins] = {1.2,1.5,1.7,2.3,1.4, 1.7, 1.5};
  
  TCanvas *c1_mlt = new TCanvas("c1_mlt","c1_mlt",width,400);
  c1_mlt->SetGrid();
  TH1D *h_mean_svd1_mlt = new TH1D("tau_offset_svd1_mlt","#tau offset, SVD1 mlt (10^{-2} ps)",3,0,3);
  h_mean_svd1_mlt->SetStats(0);
  h_mean_svd1_mlt->SetFillColor(38);
  h_mean_svd1_mlt->SetMarkerStyle(20);
  h_mean_svd1_mlt->SetMarkerSize(1.3);
  h_mean_svd1_mlt->SetMarkerColor(kBlue);
  h_mean_svd1_mlt->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd1_mlt->Fill(Modes[i],vals_svd1_mlt[i]);
    h_mean_svd1_mlt->SetBinError(i+1,errs_svd1_mlt[i]);
  }
  h_mean_svd1_mlt->LabelsDeflate();
  h_mean_svd1_mlt->GetXaxis()->SetLabelSize(0.08);
  h_mean_svd1_mlt->GetYaxis()->SetLabelSize(0.06);
  h_mean_svd1_mlt->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean_svd1_mlt->Draw("e");
  c1_mlt->Print("pics/tau_offset_svd1_mlt.eps");
  
  // * tau offset SVD2 single * //
  double vals_svd2_sgl[nbins] = {-5.7,-6.0,0.4,-3.7,-6.8,-8.5,0.3};
  double errs_svd2_sgl[nbins] = { 0.9, 1.2,1.2, 1.7, 1.0, 1.3,1.2};
  
  TCanvas *c2_sgl = new TCanvas("c2_sgl","c2_sgl",width,400);
  c2_sgl->SetGrid();
  TH1D *h_mean_svd2_sgl = new TH1D("tau_offset_svd2_sgl","#tau offset, SVD2 sgl (10^{-2} ps)",3,0,3);
  h_mean_svd2_sgl->SetStats(0);
  h_mean_svd2_sgl->SetFillColor(38);
  h_mean_svd2_sgl->SetMarkerStyle(20);
  h_mean_svd2_sgl->SetMarkerSize(1.3);
  h_mean_svd2_sgl->SetMarkerColor(kBlue);
  h_mean_svd2_sgl->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd2_sgl->Fill(Modes[i],vals_svd2_sgl[i]);
    h_mean_svd2_sgl->SetBinError(i+1,errs_svd2_sgl[i]);
  }
  h_mean_svd2_sgl->LabelsDeflate();
  h_mean_svd2_sgl->GetXaxis()->SetLabelSize(0.08);
  h_mean_svd2_sgl->GetYaxis()->SetLabelSize(0.06);
  h_mean_svd2_sgl->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean_svd2_sgl->Draw("e");
  c2_sgl->Print("pics/tau_offset_svd2_sgl.eps");
  
  // * tau offset SVD2 multiple * //
  double vals_svd2_mlt[nbins] = {-2.2,-1.5,-0.9,-0.8,-2.3,-2.5,-1.1};
  double errs_svd2_mlt[nbins] = { 0.4, 0.5, 0.6, 0.8, 0.5, 0.6, 0.6};
  
  TCanvas *c2_mlt = new TCanvas("c2_mlt","c2_mlt",width,400);
  c2_mlt->SetGrid();
  TH1D *h_mean_svd2_mlt = new TH1D("tau_offset_svd2_mlt","#tau offset, SVD2 mlt (10^{-2} ps)",3,0,3);
  h_mean_svd2_mlt->SetStats(0);
  h_mean_svd2_mlt->SetFillColor(38);
  h_mean_svd2_mlt->SetMarkerStyle(20);
  h_mean_svd2_mlt->SetMarkerSize(1.3);
  h_mean_svd2_mlt->SetMarkerColor(kBlue);
  h_mean_svd2_mlt->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd2_mlt->Fill(Modes[i],vals_svd2_mlt[i]);
    h_mean_svd2_mlt->SetBinError(i+1,errs_svd2_mlt[i]);
  }
  h_mean_svd2_mlt->LabelsDeflate();
  h_mean_svd2_mlt->GetXaxis()->SetLabelSize(0.08);
  h_mean_svd2_mlt->GetYaxis()->SetLabelSize(0.06);
  h_mean_svd2_mlt->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean_svd2_mlt->Draw("e");
  c2_mlt->Print("pics/tau_offset_svd2_mlt.eps");

//  // tau sampled pi0
//  char* fits[2] = {"f_{bkg}=0","f_{bkg}=0.4"};
//  TCanvas *c_tau_sampled_offset = new TCanvas("c_tau_sampled_offset","c_tau_sampled_offset",400,400);
//  double tau_smpl_offset_600[2]     = {-2.5,2.1};
//  double tau_smpl_offset_600_err[2] = { 3.4,5.4};
//  double tau_smpl_offset_300[2]     = {-2.6,2.2};
//  double tau_smpl_offset_300_err[2] = { 3.6,5.2};
//  c_tau_sampled_offset->SetGrid();
//  TH1D *h_tau_smpl_offset_600 = new TH1D("h_tau_smpl_offset_600","#tau offset (10^{-2} ps)",3,0,3);
//  h_tau_smpl_offset_600->SetStats(0);
//  h_tau_smpl_offset_600->SetFillColor(38);
//  h_tau_smpl_offset_600->SetMarkerStyle(20);
//  h_tau_smpl_offset_600->SetMarkerSize(1.3);
//  h_tau_smpl_offset_600->SetMarkerColor(kBlue);
//  h_tau_smpl_offset_600->SetBit(TH1::kCanRebin);
//  for(int i=0; i<2; i++){
//    h_tau_smpl_offset_600->Fill(fits[i],tau_smpl_offset_600[i]);
//    h_tau_smpl_offset_600->SetBinError(i+1,tau_smpl_offset_600_err[i]);
//  }
//  h_tau_smpl_offset_600->LabelsDeflate();
//  h_tau_smpl_offset_600->GetXaxis()->SetLabelSize(0.08);
//  h_tau_smpl_offset_600->GetYaxis()->SetLabelSize(0.06);
//  h_tau_smpl_offset_600->GetYaxis()->SetRangeUser(offmin,offmax);

//  TH1D *h_tau_smpl_offset_300 = new TH1D("h_tau_smpl_offset_300","#tau offset (10^{-2} ps)",3,0,3);
//  h_tau_smpl_offset_300->SetStats(0);
//  h_tau_smpl_offset_300->SetFillColor(38);
//  h_tau_smpl_offset_300->SetMarkerStyle(20);
//  h_tau_smpl_offset_300->SetMarkerSize(1.3);
//  h_tau_smpl_offset_300->SetMarkerColor(kRed);
//  h_tau_smpl_offset_300->SetBit(TH1::kCanRebin);
//  for(int i=0; i<2; i++){
//    h_tau_smpl_offset_300->Fill(fits[i],tau_smpl_offset_300[i]);
//    h_tau_smpl_offset_300->SetBinError(i+1,tau_smpl_offset_300_err[i]);
//  }

//  h_tau_smpl_offset_600->Draw("e");
//  h_tau_smpl_offset_300->Draw("same");
//  c_tau_sampled_offset->Print("pics/tau_sampled_offset.eps");

  // sin cos offset for full dataset
  const int Nbins = 5;
  char* Modes[Nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega","#rho"};
  double vals_sin_pt[Nbins] = {0.2,0.4,0.2,-0.8,-1.0};
  double errs_sin_pt[Nbins] = {0.7,0.9,1.3, 0.7, 0.9};
  double vals_sin[Nbins]    = {4.2,3.5,3.7,1.7,1.3};
  double errs_sin[Nbins]    = {1.1,1.5,2.2,1.3,1.6};
  double vals_cos_pt[Nbins] = {-5.0,-4.7,-6.0,-5.5,-4.3};
  double errs_cos_pt[Nbins] = { 0.9, 1.2, 1.7, 1.0, 1.3};
  double vals_cos[Nbins]    = {-4.3,-3.2,-0.2,-4.2,-3.1};
  double errs_cos[Nbins]    = { 1.6, 2.1, 3.0, 1.8, 2.2};

  TCanvas *c_s_fullds = new TCanvas("c_s_fullds","c_s_fullds",400,400);
  c_s_fullds->SetGrid();
  TH1D *h_sin_fullds_pt = new TH1D("h_sin_fullds_pt","sin offset (10^{-2})",3,0,3);
  h_sin_fullds_pt->SetStats(0);
  h_sin_fullds_pt->SetFillColor(38);
  h_sin_fullds_pt->SetMarkerStyle(20);
  h_sin_fullds_pt->SetMarkerSize(1.3);
  h_sin_fullds_pt->SetMarkerColor(kBlue);
  h_sin_fullds_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin_fullds_pt->Fill(Modes[i],vals_sin_pt[i]);
    h_sin_fullds_pt->SetBinError(i+1,errs_sin_pt[i]);
  }
  h_sin_fullds_pt->LabelsDeflate();
  h_sin_fullds_pt->GetXaxis()->SetLabelSize(0.08);
  h_sin_fullds_pt->GetYaxis()->SetLabelSize(0.06);
  h_sin_fullds_pt->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_sin_fullds = new TH1D("h_sin_fullds","sin offset (10^{-2})",3,0,3);
  h_sin_fullds->SetStats(0);
  h_sin_fullds->SetFillColor(38);
  h_sin_fullds->SetMarkerStyle(20);
  h_sin_fullds->SetMarkerSize(1.3);
  h_sin_fullds->SetMarkerColor(kRed);
//  h_sin_fullds->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin_fullds->Fill(Modes[i],vals_sin[i]);
    h_sin_fullds->SetBinError(i+1,errs_sin[i]);
  }

  h_sin_fullds_pt->Draw("e");
  h_sin_fullds->Draw("same");
  c_s_fullds->Print("pics/sin_full_ds_offset.eps");

  TCanvas *c_c_fullds = new TCanvas("c_c_fullds","c_c_fullds",400,400);
  c_c_fullds->SetGrid();
  TH1D *h_cos_fullds_pt = new TH1D("h_cos_fullds_pt","cos offset (10^{-2})",3,0,3);
  h_cos_fullds_pt->SetStats(0);
  h_cos_fullds_pt->SetFillColor(38);
  h_cos_fullds_pt->SetMarkerStyle(20);
  h_cos_fullds_pt->SetMarkerSize(1.3);
  h_cos_fullds_pt->SetMarkerColor(kBlue);
  h_cos_fullds_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos_fullds_pt->Fill(Modes[i],vals_cos_pt[i]);
    h_cos_fullds_pt->SetBinError(i+1,errs_cos_pt[i]);
  }
  h_cos_fullds_pt->LabelsDeflate();
  h_cos_fullds_pt->GetXaxis()->SetLabelSize(0.08);
  h_cos_fullds_pt->GetYaxis()->SetLabelSize(0.06);
  h_cos_fullds_pt->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_cos_fullds = new TH1D("h_cos_fullds","cos offset (10^{-2})",3,0,3);
  h_cos_fullds->SetStats(0);
  h_cos_fullds->SetFillColor(38);
  h_cos_fullds->SetMarkerStyle(20);
  h_cos_fullds->SetMarkerSize(1.3);
  h_cos_fullds->SetMarkerColor(kRed);
  h_cos_fullds->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos_fullds->Fill(Modes[i],vals_cos[i]);
    h_cos_fullds->SetBinError(i+1,errs_cos[i]);
  }

  h_cos_fullds_pt->Draw("e");
  h_cos_fullds->Draw("same");
  c_c_fullds->Print("pics/cos_full_ds_offset.eps");

//  // sin cos sampled pi0
//  string Fits[4] = {string("f_{b}=0, PT"),string("f_{b}=0"),string("f_{b}=0.4, PT"),string("f_{b}=0.4")};
//  TCanvas *c_sin_sampled_offset = new TCanvas("c_sin_sampled_offset","c_sin_sampled_offset",400,400);
//  double sin_smpl_offset_600[4]     = {0.3,5.2,10.4,13.0};
//  double sin_smpl_offset_600_err[4] = {0.7,1.2, 0.9, 1.5};
//  double sin_rms_smpl_offset_600[4] = {1.499,2.728,2.116,3.482};
//  double cos_smpl_offset_600[4]     = {-4.7,-3.5,1.7,2.0};
//  double cos_smpl_offset_600_err[4] = { 0.9, 1.7,1.3,2.0};
//  double cos_rms_smpl_offset_600[4] = {2.096,3.869,2.815,4.608};
//  double sin_smpl_offset_300[4]     = {0.6,6.2,13.1,15.0};
//  double sin_smpl_offset_300_err[4] = {0.7,1.2, 1.1, 1.5};
//  double sin_rms_smpl_offset_300[4] = {2.148,3.908,3.661,4.925};
//  double cos_smpl_offset_300[4]     = {-4.6,-1.4,3.2,3.8};
//  double cos_smpl_offset_300_err[4] = { 1.0, 1.8,1.4,2.1};
//  double cos_rms_smpl_offset_300[4] = {3.015,5.76,4.349,6.802};
//  // sin
//  c_sin_sampled_offset->SetGrid();
//  TH1D *h_sin_smpl_offset_600 = new TH1D("h_sin_smpl_offset_600","sin offset (10^{-2})",3,0,3);
//  h_sin_smpl_offset_600->SetStats(0);
//  h_sin_smpl_offset_600->SetFillColor(38);
//  h_sin_smpl_offset_600->SetMarkerStyle(20);
//  h_sin_smpl_offset_600->SetMarkerSize(1.3);
//  h_sin_smpl_offset_600->SetMarkerColor(kBlue);
//  h_sin_smpl_offset_600->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_sin_smpl_offset_600->Fill(Fits[i].c_str(),sin_smpl_offset_600[i]);
//    h_sin_smpl_offset_600->SetBinError(i+1,sin_smpl_offset_600_err[i]);
//  }
//  h_sin_smpl_offset_600->LabelsDeflate();
//  h_sin_smpl_offset_600->GetXaxis()->SetLabelSize(0.08);
//  h_sin_smpl_offset_600->GetYaxis()->SetLabelSize(0.06);
//  h_sin_smpl_offset_600->GetYaxis()->SetRangeUser(offmin-5,offmax+5);

//  TH1D *h_sin_smpl_offset_300 = new TH1D("h_sin_smpl_offset_300","sin offset (10^{-2})",3,0,3);
//  h_sin_smpl_offset_300->SetStats(0);
//  h_sin_smpl_offset_300->SetFillColor(38);
//  h_sin_smpl_offset_300->SetMarkerStyle(20);
//  h_sin_smpl_offset_300->SetMarkerSize(1.3);
//  h_sin_smpl_offset_300->SetMarkerColor(kBlack);
//  h_sin_smpl_offset_300->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_sin_smpl_offset_300->Fill(Fits[i].c_str(),sin_smpl_offset_300[i]);
//    h_sin_smpl_offset_300->SetBinError(i+1,sin_smpl_offset_300_err[i]);
//  }
//  h_sin_smpl_offset_600->Draw("e");
//  h_sin_smpl_offset_300->Draw("same");
//  c_sin_sampled_offset->Print("pics/sin_sampled_offset.eps");

//  // cos
//  TCanvas *c_cos_sampled_offset = new TCanvas("c_cos_sampled_offset","c_cos_sampled_offset",400,400);
//  c_cos_sampled_offset->SetGrid();
//  TH1D *h_cos_smpl_offset_600 = new TH1D("h_cos_smpl_offset_600","cos offset (10^{-2})",3,0,3);
//  h_cos_smpl_offset_600->SetStats(0);
//  h_cos_smpl_offset_600->SetFillColor(38);
//  h_cos_smpl_offset_600->SetMarkerStyle(20);
//  h_cos_smpl_offset_600->SetMarkerSize(1.3);
//  h_cos_smpl_offset_600->SetMarkerColor(kBlue);
//  h_cos_smpl_offset_600->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_cos_smpl_offset_600->Fill(Fits[i].c_str(),cos_smpl_offset_600[i]);
//    h_cos_smpl_offset_600->SetBinError(i+1,cos_smpl_offset_600_err[i]);
//  }
//  h_cos_smpl_offset_600->LabelsDeflate();
//  h_cos_smpl_offset_600->GetXaxis()->SetLabelSize(0.08);
//  h_cos_smpl_offset_600->GetYaxis()->SetLabelSize(0.06);
//  h_cos_smpl_offset_600->GetYaxis()->SetRangeUser(offmin,offmax);

//  TH1D *h_cos_smpl_offset_300 = new TH1D("h_cos_smpl_offset_300","cos offset (10^{-2})",3,0,3);
//  h_cos_smpl_offset_300->SetStats(0);
//  h_cos_smpl_offset_300->SetFillColor(38);
//  h_cos_smpl_offset_300->SetMarkerStyle(20);
//  h_cos_smpl_offset_300->SetMarkerSize(1.3);
//  h_cos_smpl_offset_300->SetMarkerColor(kBlack);
//  h_cos_smpl_offset_300->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_cos_smpl_offset_300->Fill(Fits[i].c_str(),cos_smpl_offset_300[i]);
//    h_cos_smpl_offset_300->SetBinError(i+1,cos_smpl_offset_300_err[i]);
//  }

//  h_cos_smpl_offset_600->Draw("e");
//  h_cos_smpl_offset_300->Draw("same");
//  c_cos_sampled_offset->Print("pics/cos_sampled_offset.eps");

//  // rms
//  TCanvas *c_rms_sampled_offset = new TCanvas("c_rms_sampled_offset","c_rms_sampled_offset",400,400);
//  c_rms_sampled_offset->SetGrid();
//  TH1D *h_rms_sin_smpl_offset_600 = new TH1D("h_rms_sin_smpl_offset_600","RMS (10^{-1})",3,0,3);
//  h_rms_sin_smpl_offset_600->SetStats(0);
//  h_rms_sin_smpl_offset_600->SetFillColor(38);
//  h_rms_sin_smpl_offset_600->SetMarkerStyle(21);
//  h_rms_sin_smpl_offset_600->SetMarkerSize(1.3);
//  h_rms_sin_smpl_offset_600->SetMarkerColor(kRed);
//  h_rms_sin_smpl_offset_600->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_rms_sin_smpl_offset_600->Fill(Fits[i].c_str(),sin_rms_smpl_offset_600[i]);
//    h_rms_sin_smpl_offset_600->SetBinError(i+1,0.1);
//  }
//  h_rms_sin_smpl_offset_600->LabelsDeflate();
//  h_rms_sin_smpl_offset_600->GetXaxis()->SetLabelSize(0.08);
//  h_rms_sin_smpl_offset_600->GetYaxis()->SetLabelSize(0.06);
//  h_rms_sin_smpl_offset_600->GetYaxis()->SetRangeUser(0,offmax);

//  TH1D *h_rms_sin_smpl_offset_300 = new TH1D("h_rms_sin_smpl_offset_300","RMS offset (10^{-1})",3,0,3);
//  h_rms_sin_smpl_offset_300->SetStats(0);
//  h_rms_sin_smpl_offset_300->SetFillColor(38);
//  h_rms_sin_smpl_offset_300->SetMarkerStyle(20);
//  h_rms_sin_smpl_offset_300->SetMarkerSize(1.3);
//  h_rms_sin_smpl_offset_300->SetMarkerColor(kRed);
//  h_rms_sin_smpl_offset_300->SetBit(TH1::kCanRebin);
//  h_rms_sin_smpl_offset_300->GetYaxis()->SetRangeUser(0,offmax);
//  for(int i=0; i<4; i++){
//    h_rms_sin_smpl_offset_300->Fill(Fits[i].c_str(),sin_rms_smpl_offset_300[i]);
//    h_rms_sin_smpl_offset_300->SetBinError(i+1,0.1);
//  }
//  TH1D *h_rms_cos_smpl_offset_600 = new TH1D("h_rms_cos_smpl_offset_600","RMS (10^{-1})",3,0,3);
//  h_rms_cos_smpl_offset_600->SetStats(0);
//  h_rms_cos_smpl_offset_600->SetFillColor(38);
//  h_rms_cos_smpl_offset_600->SetMarkerStyle(21);
//  h_rms_cos_smpl_offset_600->SetMarkerSize(1.3);
//  h_rms_cos_smpl_offset_600->SetMarkerColor(kBlue);
//  h_rms_cos_smpl_offset_600->SetBit(TH1::kCanRebin);
//  for(int i=0; i<4; i++){
//    h_rms_cos_smpl_offset_600->Fill(Fits[i].c_str(),cos_rms_smpl_offset_600[i]);
//    h_rms_cos_smpl_offset_600->SetBinError(i+1,0.1);
//  }
//  h_rms_cos_smpl_offset_600->LabelsDeflate();
//  h_rms_cos_smpl_offset_600->GetXaxis()->SetLabelSize(0.08);
//  h_rms_cos_smpl_offset_600->GetYaxis()->SetLabelSize(0.06);
//  h_rms_cos_smpl_offset_600->GetYaxis()->SetRangeUser(0,offmax);

//  TH1D *h_rms_cos_smpl_offset_300 = new TH1D("h_rms_cos_smpl_offset_300","RMS(cos) offset (10^{-1})",3,0,3);
//  h_rms_cos_smpl_offset_300->SetStats(0);
//  h_rms_cos_smpl_offset_300->SetFillColor(38);
//  h_rms_cos_smpl_offset_300->SetMarkerStyle(20);
//  h_rms_cos_smpl_offset_300->SetMarkerSize(1.3);
//  h_rms_cos_smpl_offset_300->SetMarkerColor(kBlue);
//  h_rms_cos_smpl_offset_300->SetBit(TH1::kCanRebin);
//  h_rms_cos_smpl_offset_300->GetYaxis()->SetRangeUser(0,offmax);
//  for(int i=0; i<4; i++){
//    h_rms_cos_smpl_offset_300->Fill(Fits[i].c_str(),cos_rms_smpl_offset_300[i]);
//    h_rms_cos_smpl_offset_300->SetBinError(i+1,0.1);
//  }

//  h_rms_cos_smpl_offset_600->Draw("e");
//  h_rms_cos_smpl_offset_300->Draw("same");
//  h_rms_sin_smpl_offset_600->Draw("same");
//  h_rms_sin_smpl_offset_300->Draw("same");
//  c_rms_sampled_offset->Print("pics/rms_sampled_offset.eps");

  // * Sampled offsets * //
  // sin cos sampled pi0
  const int FourBins = 4;
  char* FourModes[FourBins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega"};

  // sin
  TCanvas *c_sin_sampled_offset_600 = new TCanvas("c_sin_sampled_offset_600","c_sin_sampled_offset_600",400,400);
  double sin_smpl_600_pt_nb[FourBins]     = {0.29,0.49,0.43,-0.60};
  double sin_smpl_600_pt_nb_err[FourBins] = {0.67,0.89,1.17, 0.73};
  double sin_smpl_600_pt[FourBins]        = {10.40,15.86,10.45,6.75};
  double sin_smpl_600_pt_err[FourBins]    = { 0.94,1.86, 1.68,0.92};
  double sin_smpl_600_nb[FourBins]        = {5.20,4.24,4.09,2.44};
  double sin_smpl_600_nb_err[FourBins]    = {1.23,1.54,2.17,1.30};
  double sin_smpl_600[FourBins]           = {12.97,9.90,14.28,10.98};
  double sin_smpl_600_err[FourBins]       = { 1.55,1.97, 2.94,1.67};

  c_sin_sampled_offset_600->SetGrid();
  TH1D *h_sin_smpl_offset_600_pt_nb = new TH1D("h_sin_smpl_offset_600_pt_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_offset_600_pt_nb->SetStats(0);
  h_sin_smpl_offset_600_pt_nb->SetFillColor(38);
  h_sin_smpl_offset_600_pt_nb->SetMarkerStyle(20);
  h_sin_smpl_offset_600_pt_nb->SetMarkerSize(1.3);
  h_sin_smpl_offset_600_pt_nb->SetMarkerColor(kBlue);
  h_sin_smpl_offset_600_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_offset_600_pt_nb->Fill(FourModes[i],sin_smpl_600_pt_nb[i]);
    h_sin_smpl_offset_600_pt_nb->SetBinError(i+1,sin_smpl_600_pt_nb_err[i]);
  }
  h_sin_smpl_offset_600_pt_nb->LabelsDeflate();
  h_sin_smpl_offset_600_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_offset_600_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_offset_600_pt_nb->GetYaxis()->SetRangeUser(-15,16);

  TH1D *h_sin_smpl_offset_600_nb = new TH1D("h_sin_smpl_offset_600_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_offset_600_nb->SetStats(0);
  h_sin_smpl_offset_600_nb->SetFillColor(38);
  h_sin_smpl_offset_600_nb->SetMarkerStyle(20);
  h_sin_smpl_offset_600_nb->SetMarkerSize(1.3);
  h_sin_smpl_offset_600_nb->SetMarkerColor(kRed);
  h_sin_smpl_offset_600_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_offset_600_nb->Fill(FourModes[i],sin_smpl_600_nb[i]);
    h_sin_smpl_offset_600_nb->SetBinError(i+1,sin_smpl_600_nb_err[i]);
  }

  TH1D *h_sin_smpl_offset_600_pt = new TH1D("h_sin_smpl_offset_600_pt","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_offset_600_pt->SetStats(0);
  h_sin_smpl_offset_600_pt->SetFillColor(38);
  h_sin_smpl_offset_600_pt->SetMarkerStyle(20);
  h_sin_smpl_offset_600_pt->SetMarkerSize(1.3);
  h_sin_smpl_offset_600_pt->SetMarkerColor(6);
  h_sin_smpl_offset_600_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_offset_600_pt->Fill(FourModes[i],sin_smpl_600_pt[i]);
    h_sin_smpl_offset_600_pt->SetBinError(i+1,sin_smpl_600_pt_err[i]);
  }

  TH1D *h_sin_smpl_offset_600 = new TH1D("h_sin_smpl_offset_600","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_offset_600->SetStats(0);
  h_sin_smpl_offset_600->SetFillColor(38);
  h_sin_smpl_offset_600->SetMarkerStyle(20);
  h_sin_smpl_offset_600->SetMarkerSize(1.3);
  h_sin_smpl_offset_600->SetMarkerColor(kBlack);
  h_sin_smpl_offset_600->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_offset_600->Fill(FourModes[i],sin_smpl_600[i]);
    h_sin_smpl_offset_600->SetBinError(i+1,sin_smpl_600_err[i]);
  }

  h_sin_smpl_offset_600_pt_nb->Draw("e");
  h_sin_smpl_offset_600_nb->Draw("same");
  h_sin_smpl_offset_600_pt->Draw("same");
  h_sin_smpl_offset_600->Draw("same");

  c_sin_sampled_offset_600->Print("pics/sin_sampled_offset_600.eps");

  // cos
  TCanvas *c_cos_sampled_offset_600 = new TCanvas("c_cos_sampled_offset_600","c_cos_sampled_offset_600",400,400);
  double cos_smpl_600_pt_nb[FourBins]     = {-4.71,-4.80,-5.81,-5.25};
  double cos_smpl_600_pt_nb_err[FourBins] = { 0.94, 1.25, 1.75, 0.99};
  double cos_smpl_600_pt[FourBins]        = {1.72,5.49,5.70,1.24};
  double cos_smpl_600_pt_err[FourBins]    = {1.25,2.11,2.19,1.39};
  double cos_smpl_600_nb[FourBins]        = {-3.50,-2.56,0.77,-2.63};
  double cos_smpl_600_nb_err[FourBins]    = { 1.74, 2.35,2.94, 1.95};
  double cos_smpl_600[FourBins]           = {2.00,-3.42,6.17,3.24};
  double cos_smpl_600_err[FourBins]       = {2.05, 3.06,4.27,2.26};
  c_cos_sampled_offset_600->SetGrid();
  TH1D *h_cos_smpl_offset_600_pt_nb = new TH1D("h_cos_smpl_offset_600_pt_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_offset_600_pt_nb->SetStats(0);
  h_cos_smpl_offset_600_pt_nb->SetFillColor(38);
  h_cos_smpl_offset_600_pt_nb->SetMarkerStyle(20);
  h_cos_smpl_offset_600_pt_nb->SetMarkerSize(1.3);
  h_cos_smpl_offset_600_pt_nb->SetMarkerColor(kBlue);
  h_cos_smpl_offset_600_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_offset_600_pt_nb->Fill(FourModes[i],cos_smpl_600_pt_nb[i]);
    h_cos_smpl_offset_600_pt_nb->SetBinError(i+1,cos_smpl_600_pt_nb_err[i]);
  }
  h_cos_smpl_offset_600_pt_nb->LabelsDeflate();
  h_cos_smpl_offset_600_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_offset_600_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_offset_600_pt_nb->GetYaxis()->SetRangeUser(-15,15);

  TH1D *h_cos_smpl_offset_600_nb = new TH1D("h_cos_smpl_offset_600_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_offset_600_nb->SetStats(0);
  h_cos_smpl_offset_600_nb->SetFillColor(38);
  h_cos_smpl_offset_600_nb->SetMarkerStyle(20);
  h_cos_smpl_offset_600_nb->SetMarkerSize(1.3);
  h_cos_smpl_offset_600_nb->SetMarkerColor(kRed);
  h_cos_smpl_offset_600_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_offset_600_nb->Fill(FourModes[i],cos_smpl_600_nb[i]);
    h_cos_smpl_offset_600_nb->SetBinError(i+1,cos_smpl_600_nb_err[i]);
  }

  TH1D *h_cos_smpl_offset_600_pt = new TH1D("h_cos_smpl_offset_600_pt","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_offset_600_pt->SetStats(0);
  h_cos_smpl_offset_600_pt->SetFillColor(38);
  h_cos_smpl_offset_600_pt->SetMarkerStyle(20);
  h_cos_smpl_offset_600_pt->SetMarkerSize(1.3);
  h_cos_smpl_offset_600_pt->SetMarkerColor(6);
  h_cos_smpl_offset_600_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_offset_600_pt->Fill(FourModes[i],cos_smpl_600_pt[i]);
    h_cos_smpl_offset_600_pt->SetBinError(i+1,cos_smpl_600_pt_err[i]);
  }

  TH1D *h_cos_smpl_offset_600 = new TH1D("h_cos_smpl_offset_600","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_offset_600->SetStats(0);
  h_cos_smpl_offset_600->SetFillColor(38);
  h_cos_smpl_offset_600->SetMarkerStyle(20);
  h_cos_smpl_offset_600->SetMarkerSize(1.3);
  h_cos_smpl_offset_600->SetMarkerColor(kBlack);
  h_cos_smpl_offset_600->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_offset_600->Fill(FourModes[i],cos_smpl_600[i]);
    h_cos_smpl_offset_600->SetBinError(i+1,cos_smpl_600_err[i]);
  }

  h_cos_smpl_offset_600_pt_nb->Draw("e");
  h_cos_smpl_offset_600_nb->Draw("same");
  h_cos_smpl_offset_600_pt->Draw("same");
  h_cos_smpl_offset_600->Draw("same");

  c_cos_sampled_offset_600->Print("pics/cos_sampled_offset_600.eps");

  //tau
  TCanvas *c_tau_sampled_offset_600 = new TCanvas("c_tau_sampled_offset_600","c_tau_sampled_offset_600",400,400);
  double tau_smpl_600_nb[FourBins]     = {-2.5, -1.8, -1.4, -3.0};
  double tau_smpl_600_nb_err[FourBins] = {0.353,0.471,0.612,0.394};
  double tau_smpl_600[FourBins]        = {1.4,  0.8,  -1.5, -3.6};
  double tau_smpl_600_err[FourBins]    = {0.563,0.937,1.215,0.785};
  c_tau_sampled_offset_600->SetGrid();
  TH1D *h_tau_smpl_offset_600_nb = new TH1D("h_tau_smpl_offset_600_nb","tau offset (10^{-2} ps)",3,0,3);
  h_tau_smpl_offset_600_nb->SetStats(0);
  h_tau_smpl_offset_600_nb->SetFillColor(38);
  h_tau_smpl_offset_600_nb->SetMarkerStyle(20);
  h_tau_smpl_offset_600_nb->SetMarkerSize(1.3);
  h_tau_smpl_offset_600_nb->SetMarkerColor(kBlue);
  h_tau_smpl_offset_600_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_tau_smpl_offset_600_nb->Fill(FourModes[i],tau_smpl_600_nb[i]);
    h_tau_smpl_offset_600_nb->SetBinError(i+1,tau_smpl_600_nb_err[i]);
  }
  h_tau_smpl_offset_600_nb->LabelsDeflate();
  h_tau_smpl_offset_600_nb->GetXaxis()->SetLabelSize(0.08);
  h_tau_smpl_offset_600_nb->GetYaxis()->SetLabelSize(0.06);
  h_tau_smpl_offset_600_nb->GetYaxis()->SetRangeUser(-10,10);

  TH1D *h_tau_smpl_offset_600 = new TH1D("h_tau_smpl_offset_600","tau offset (10^{-2} ps)",3,0,3);
  h_tau_smpl_offset_600->SetStats(0);
  h_tau_smpl_offset_600->SetFillColor(38);
  h_tau_smpl_offset_600->SetMarkerStyle(20);
  h_tau_smpl_offset_600->SetMarkerSize(1.3);
  h_tau_smpl_offset_600->SetMarkerColor(kRed);
  h_tau_smpl_offset_600->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_tau_smpl_offset_600->Fill(FourModes[i],tau_smpl_600[i]);
    h_tau_smpl_offset_600->SetBinError(i+1,tau_smpl_600_err[i]);
  }

  h_tau_smpl_offset_600_nb->Draw("e");
  h_tau_smpl_offset_600->Draw("same");

  c_tau_sampled_offset_600->Print("pics/tau_sampled_offset_600.eps");

  return;
}
