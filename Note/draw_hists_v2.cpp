const bool draw_tau_svd1_sgl = false;
const bool draw_tau_svd1_mlt = false;
const bool draw_tau_svd2_sgl = false;
const bool draw_tau_svd2_mlt = false;
const bool draw_tau          = true;
const bool draw_tau_small    = true;

const bool draw_sin_cos      = false;

const double offmin_small = -50;
const double offmax_small =  50;

const double offmin = -10;
const double offmax =  10;

const int nbins = 7;
char* Modes[nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega","#eta`","D*#pi^{0}","D*#eta"};
void draw_hists_v2(void){
  // * tau offset GenMC * //
  if(draw_tau_small){


  const int width = 400;
  double vals_nobkg_small[nbins] = {-7.45,-3.96,14.28,-8.76,-2.6, 5.73, -2.35};
  double errs_nobkg_small[nbins] = { 5.11, 8.83,15.48, 5.57,27.27,12.06,20.0};

  double vals_small[nbins] = {-2.75,16.07,20.74,2.36,27.39,37.31,4.12};
  double errs_small[nbins] = { 5.46,10.36,17.28,6.26,25.33,14.04,22.63};

  TCanvas *c1 = new TCanvas("c1_small","c1_small",400,400);
  c1->SetGrid();
  TH1D *h_mean_small = new TH1D("tau_offset_small","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_small->SetStats(0);
  h_mean_small->SetFillColor(38);
  h_mean_small->SetMarkerStyle(20);
  h_mean_small->SetMarkerSize(1.3);
  h_mean_small->SetMarkerColor(kBlue);
  h_mean_small->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_small->Fill(Modes[i],vals_small[i]);
    h_mean_small->SetBinError(i+1,errs_small[i]);
  }
  TH1D *h_mean_nobkg_small = new TH1D("tau_offset_small_nobkg","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_nobkg_small->SetStats(0);
  h_mean_nobkg_small->SetFillColor(38);
  h_mean_nobkg_small->SetMarkerStyle(20);
  h_mean_nobkg_small->SetMarkerSize(1.3);
  h_mean_nobkg_small->SetMarkerColor(kRed);
  h_mean_nobkg_small->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_nobkg_small->Fill(Modes[i],vals_nobkg_small[i]);
    h_mean_nobkg_small->SetBinError(i+1,errs_nobkg_small[i]);
  }

  h_mean_small->LabelsDeflate();
  h_mean_small->GetXaxis()->SetLabelSize(0.08);
  h_mean_small->GetYaxis()->SetRangeUser(offmin_small,offmax_small);
  h_mean_small->GetYaxis()->SetLabelSize(0.06);
  h_mean_small->Draw("e");
  h_mean_nobkg_small->Draw("esame");
  c1->Print("pics/tau_offset_small.eps");
  }

  // * tau offset total * //
  if(draw_tau){
  const int nbins = 7;
  const double offmin = -10;
  const double offmax =  10;
  const int width = 400;

  double vals_nobkg[nbins] = {-2.713,-2.048,0.153,-0.843,-1.610,-2.869,-4.550};
  double errs_nobkg[nbins] = { 0.599, 0.746,0.790, 0.750, 2.860, 1.599, 2.273};

  double vals[nbins] = {-3.134,-2.573,-0.628,-1.691,-0.328,-1.819,-5.060};
  double errs[nbins] = { 0.631, 0.811, 0.828, 0.781, 3.100, 1.802, 2.416};

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

  TH1D *h_mean_nobkg = new TH1D("tau_offset_nobkg","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_nobkg->SetStats(0);
  h_mean_nobkg->SetFillColor(38);
  h_mean_nobkg->SetMarkerStyle(20);
  h_mean_nobkg->SetMarkerSize(1.3);
  h_mean_nobkg->SetMarkerColor(kRed);
  h_mean_nobkg->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_nobkg->Fill(Modes[i],vals_nobkg[i]);
    h_mean_nobkg->SetBinError(i+1,errs_nobkg[i]);
  }

  h_mean->LabelsDeflate();
  h_mean->GetXaxis()->SetLabelSize(0.08);
  h_mean->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean->GetYaxis()->SetLabelSize(0.06);
  h_mean->Draw("e");
  h_mean_nobkg->Draw("esame");
  c1->Print("pics/tau_offset.eps");
  }

  // * tau offset SVD1 single * //
  if(draw_tau_svd1_sgl){
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
  }
  
  // * tau offset SVD1 multiple * //
  if(draw_tau_svd1_mlt){
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
  }
  
  // * tau offset SVD2 single * //
  if(draw_tau_svd2_sgl){
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
  }
  
  // * tau offset SVD2 multiple * //
  if(draw_tau_svd2_mlt){
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
  }

  // sin cos offset for full dataset
  if(draw_sin_cos){
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
  }
  return;
}
