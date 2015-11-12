const bool draw_tau_svd1_sgl = false;
const bool draw_tau_svd1_mlt = false;
const bool draw_tau_svd2_sgl = false;
const bool draw_tau_svd2_mlt = false;
const bool draw_tau          = true;
const bool draw_tau_small    = true;
const bool draw_tau_h0       = true;

const bool draw_sin_cos       = true;
const bool draw_sin_cos_small = true;

const double offmin_small = -50;
const double offmax_small =  50;

const double offmin = -10;
const double offmax =  10;

const int nbins = 7;
string Modes[nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega","#eta`","D*#pi^{0}","D*#eta"};
//char* h0Modes[2] = {"h^{0}#rightarrow#gamma#gamma","h^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}"};
string h0Modes[2] = {"Single","Multiple"};

const int Nbins = 5;
string  cpvModes[Nbins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega","#eta'"};


void draw_hists_ww(void){
  // * h0 tau offset GenMC * //
  if(draw_tau_h0){
//  const int width = 400;
  const double vals_h0[2] = {1.209,-4.750};
  const double errs_h0[2] = {4.043, 5.079};

  const double vals_nobkg_h0[2] = {-1.524,-8.189};
  const double errs_nobkg_h0[2] = { 3.594, 4.555};

  const double vals_diff_h0[2] = {vals_h0[0]-vals_nobkg_h0[0],vals_h0[1]-vals_nobkg_h0[1]};
  const double errs_diff_h0[2] = {sqrt(errs_h0[0]*errs_h0[0]-errs_nobkg_h0[0]*errs_nobkg_h0[0]),sqrt(errs_h0[1]*errs_h0[1]-errs_nobkg_h0[1]*errs_nobkg_h0[1])};

  TCanvas *c1_h0_diff = new TCanvas("c1_h0_diff","c1_h0_diff",400,400);
  c1_h0_diff->SetGrid();
  TH1D *h_mean_h0 = new TH1D("tau_offset_small","#tau_{FIT}-#tau_{FIT}^{NoBkg} (10^{-2} ps)",3,0,3);
  h_mean_h0->SetStats(0);
  h_mean_h0->SetFillColor(38);
  h_mean_h0->SetMarkerStyle(20);
  h_mean_h0->SetMarkerSize(1.3);
  h_mean_h0->SetMarkerColor(kBlue);
//  h_mean_h0->SetBit(TH1::kCanRebin);
  for(int i=0; i<2; i++){
    h_mean_h0->Fill(h0Modes[i].c_str(),vals_diff_h0[i]);
    h_mean_h0->SetBinError(i+1,errs_diff_h0[i]);
  }
  h_mean_h0->LabelsDeflate();
  h_mean_h0->GetXaxis()->SetLabelSize(0.08);
  h_mean_h0->GetYaxis()->SetRangeUser(-10.,10.);
  h_mean_h0->GetYaxis()->SetLabelSize(0.06);
  h_mean_h0->Draw("e");
  c1_h0_diff->Print("pics/tau_offset_h0_ww.eps");
  }

  // * tau offset GenMC * //
  if(draw_tau_small){
  double vals_nobkg_small[nbins] = {-4.734,7.980,0.723,-10.76,11.469,1.887,-7.844};
  double errs_nobkg_small[nbins] = { 4.418,8.392,12.90, 4.965,27.270,10.62,17.528};

  double vals_small[nbins] = {-4.586,17.432, 4.333,-8.458,26.88,14.365,-16.35};
  double errs_small[nbins] = { 4.790, 9.828,14.272, 5.535,24.84,13.416, 20.55};

  double vals_small_fixed_cont[nbins] = {-5.322,14.375, 2.327,-9.891,27.606,15.182,-18.491};
  double errs_small_fixed_cont[nbins] = { 4.800, 9.812,14.256, 5.535,24.831,13.415, 20.5742};

  TCanvas *c1_small = new TCanvas("c1_small","c1_small",400,400);
  c1_small->SetGrid();
  TH1D *h_mean_small = new TH1D("tau_offset_small","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_small->SetStats(0);
  h_mean_small->SetFillColor(38);
  h_mean_small->SetMarkerStyle(20);
  h_mean_small->SetMarkerSize(1.3);
  h_mean_small->SetMarkerColor(kBlue);
//  h_mean_small->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_small->Fill(Modes[i].c_str(),vals_small[i]);
    h_mean_small->SetBinError(i+1,errs_small[i]);
  }
  TH1D *h_mean_nobkg_small = new TH1D("tau_offset_small_nobkg","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_nobkg_small->SetStats(0);
  h_mean_nobkg_small->SetFillColor(38);
  h_mean_nobkg_small->SetMarkerStyle(20);
  h_mean_nobkg_small->SetMarkerSize(1.3);
  h_mean_nobkg_small->SetMarkerColor(kRed);
//  h_mean_nobkg_small->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_nobkg_small->Fill(Modes[i].c_str(),vals_nobkg_small[i]);
    h_mean_nobkg_small->SetBinError(i+1,errs_nobkg_small[i]);
  }

  TH1D *h_mean_fixed_cont_small = new TH1D("tau_offset_small_fixed_cont","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_fixed_cont_small->SetStats(0);
  h_mean_fixed_cont_small->SetFillColor(38);
  h_mean_fixed_cont_small->SetMarkerStyle(20);
  h_mean_fixed_cont_small->SetMarkerSize(1.3);
  h_mean_fixed_cont_small->SetMarkerColor(kBlack);
//  h_mean_fixed_cont_small->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_fixed_cont_small->Fill(Modes[i].c_str(),vals_small_fixed_cont[i]);
    h_mean_fixed_cont_small->SetBinError(i+1,errs_small_fixed_cont[i]);
  }

  h_mean_small->LabelsDeflate();
  h_mean_small->GetXaxis()->SetLabelSize(0.08);
  h_mean_small->GetYaxis()->SetRangeUser(offmin_small,offmax_small);
  h_mean_small->GetYaxis()->SetLabelSize(0.06);
  h_mean_small->Draw("e");
  h_mean_nobkg_small->Draw("esame");
//  h_mean_fixed_cont_small->Draw("esame");
  c1_small->Print("pics/tau_offset_small_ww.eps");
  }

  // * tau offset total * //
  if(draw_tau){
  const int nbins = 7;
  const double offmin = -10;
  const double offmax =  10;

  double vals_nobkg_cut50[nbins] = {-2.865,-1.336,0.698,-1.101,-2.145,-2.159,-5.283};
  double errs_nobkg_cut50[nbins] = { 0.423, 0.534,0.761, 0.489, 2.878, 1.735, 2.440};

//  double vals_nobkg[nbins] = {-0.710,1.025,1.239,-0.458,-1.586,0.270,-3.814};
//  double errs_nobkg[nbins] = { 0.432,0.548,0.758, 0.487, 2.863,1.760, 2.488};

// With corrected chisq for multiple //
  double vals_nobkg[nbins] = {-0.710,1.025,-0.861,-0.458,-4.878,0.270,-3.814};
  double errs_nobkg[nbins] = { 0.432,0.548, 0.756, 0.487, 2.863,1.760, 2.488};

  TCanvas *c1 = new TCanvas("c1","c1",400,400);
  c1->SetGrid();
  TH1D *h_mean_nobkg = new TH1D("tau_offset_nobkg","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_nobkg->SetStats(0);
  h_mean_nobkg->SetFillColor(38);
  h_mean_nobkg->SetMarkerStyle(20);
  h_mean_nobkg->SetMarkerSize(1.3);
  h_mean_nobkg->SetMarkerColor(kRed);
//  h_mean_nobkg->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_nobkg->Fill(Modes[i].c_str(),vals_nobkg[i]);
    h_mean_nobkg->SetBinError(i+1,errs_nobkg[i]);
  }

  TH1D *h_mean_nobkg_cut50 = new TH1D("tau_offset_nobkg_cut50","#tau offset (10^{-2} ps)",3,0,3);
  h_mean_nobkg_cut50->SetStats(0);
  h_mean_nobkg_cut50->SetFillColor(38);
  h_mean_nobkg_cut50->SetMarkerStyle(20);
  h_mean_nobkg_cut50->SetMarkerSize(1.3);
  h_mean_nobkg_cut50->SetMarkerColor(kBlack);
//  h_mean_nobkg_cut50->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_nobkg_cut50->Fill(Modes[i].c_str(),vals_nobkg_cut50[i]);
    h_mean_nobkg_cut50->SetBinError(i+1,errs_nobkg_cut50[i]);
  }

  h_mean_nobkg->LabelsDeflate();
  h_mean_nobkg->GetXaxis()->SetLabelSize(0.08);
  h_mean_nobkg->GetYaxis()->SetRangeUser(offmin,offmax);
  h_mean_nobkg->GetYaxis()->SetLabelSize(0.06);
  h_mean_nobkg->Draw("e");
//  h_mean_nobkg_cut50->Draw("esame");
  c1->Print("pics/tau_offset_ww1.eps");
  }

  // * tau offset SVD1 single * //
  const int width = 400;
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
//  h_mean_svd1_sgl->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd1_sgl->Fill(Modes[i].c_str(),vals_svd1_sgl[i]);
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
//  h_mean_svd1_mlt->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd1_mlt->Fill(Modes[i].c_str(),vals_svd1_mlt[i]);
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
//  h_mean_svd2_sgl->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd2_sgl->Fill(Modes[i].c_str(),vals_svd2_sgl[i]);
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
//  h_mean_svd2_mlt->SetBit(TH1::kCanRebin);
  for(int i=0; i<nbins; i++){
    h_mean_svd2_mlt->Fill(Modes[i].c_str(),vals_svd2_mlt[i]);
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
  double vals_sin_tagcut[Nbins] = { 0.755402,0.742490,0.739773, 0.707418,0.719340};
  double errs_sin_tagcut[Nbins] = { 0.012988,0.016513,0.023106, 0.015077,0.085305};
  double vals_cos_tagcut[Nbins] = {-0.681964,0.735194,0.685695,-0.712986,0.718747};
  double errs_cos_tagcut[Nbins] = { 0.018976,0.023841,0.033229, 0.021542,0.117920};

  double vals_sin[Nbins] =        { 0.755427,0.742642,0.739840, 0.707484,0.719310};
  double errs_sin[Nbins] =        { 0.012989,0.016516,0.023109, 0.015077,0.085323};
  double vals_cos[Nbins] =        {-0.681922,0.735340,0.685587,-0.712974,0.718781};
  double errs_cos[Nbins] =        { 0.018978,0.023846,0.033235, 0.021543,0.117937};

  const double true_sin = 0.719340;
  const double true_cos = 0.694658;

  for(int i=0; i<Nbins; i++){
    vals_sin[i] = (fabs(vals_sin[i]) - true_sin)*100.;
    vals_cos[i] = (fabs(vals_cos[i]) - true_cos)*100.;
    errs_sin[i] *= 100.; errs_cos[i] *= 100.;
    cout << vals_sin[i] << " " << vals_cos[i] << endl;
  }

//  Before bin map shifting
//  double vals_sin_pt[Nbins] = { 0.759571,0.753903,0.735926, 0.734241,0.706593};
//  double errs_sin_pt[Nbins] = { 0.008120,0.010199,0.014332, 0.009301,0.055115};
//  double vals_cos_pt[Nbins] = {-0.721990,0.740594,0.685957,-0.717516,0.701805};
//  double errs_cos_pt[Nbins] = { 0.011673,0.014737,0.020601, 0.013286,0.075667};

// Without Ki correction
//  double vals_sin_pt[Nbins] = {0.7672690,0.753883,0.739913,0.743627,0.723939};
//  double errs_sin_pt[Nbins] = {0.0080795,0.010201,0.014308,0.009246,0.054549};
//  double vals_cos_pt[Nbins] = {0.7412460,0.760023,0.700753,0.735432,0.690649};
//  double errs_cos_pt[Nbins] = {0.0119106,0.015006,0.020987,0.013501,0.077553};

  double vals_sin_pt[Nbins] = {0.765481,0.751779,0.738390,0.74182,0.721523};
  double errs_sin_pt[Nbins] = {0.00806399,0.010181,0.014283,0.00922778,0.0544627};
  double vals_cos_pt[Nbins] = {0.731868,0.750344,0.699179,0.734377,0.690792};
  double errs_cos_pt[Nbins] = {0.0117726,0.0148286,0.020963,0.0134776,0.0774486};

  for(int i=0; i<Nbins; i++){
    vals_sin_pt[i] = (fabs(vals_sin_pt[i]) - true_sin)*100.;
    vals_cos_pt[i] = (fabs(vals_cos_pt[i]) - true_cos)*100.;
    errs_sin_pt[i] *= 100.; errs_cos_pt[i] *= 100.;
    cout << vals_sin_pt[i] << " " << vals_cos_pt[i] << endl;
  }

  double vals_sin_pb[Nbins] = { 0.761995,0.743690,0.746623, 0.707095,0.737619};
  double errs_sin_pb[Nbins] = { 0.012938,0.016492,0.023002, 0.015045,0.0849551};
  double vals_cos_pb[Nbins] = {-0.696160,0.729472,0.696994,-0.715748,0.693525};
  double errs_cos_pb[Nbins] = { 0.018940,0.023865,0.033245, 0.021557,0.118071};

  for(int i=0; i<Nbins; i++){
    vals_sin_pb[i] = (fabs(vals_sin_pb[i]) - true_sin)*100.;
    vals_cos_pb[i] = (fabs(vals_cos_pb[i]) - true_cos)*100.;
    errs_sin_pb[i] *= 100.; errs_cos_pb[i] *= 100.;
    cout << vals_sin_pb[i] << " " << vals_cos_pb[i] << endl;
  }

//  Before bin map shifting
//  double vals_sin_mc[Nbins] = { 0.7639310,0.754113,0.732412, 0.7341620,0.725278};
//  double errs_sin_mc[Nbins] = { 0.0081041,0.010190,0.014313, 0.0092915,0.054689};
//  double vals_cos_mc[Nbins] = {-0.7329740,0.741750,0.696646,-0.7213620,0.699351};
//  double errs_cos_mc[Nbins] = { 0.0116688,0.014702,0.020653, 0.0133058,0.075229};

// Without Ki correction
//  double vals_sin_mc[Nbins] = {0.773290,0.753599,0.737002,0.745775,0.725760};
//  double errs_sin_mc[Nbins] = {0.008060,0.010188,0.014297,0.009240,0.054544};
//  double vals_cos_mc[Nbins] = {0.752553,0.755853,0.712775,0.738804,0.684861};
//  double errs_cos_mc[Nbins] = {0.011900,0.015006,0.020922,0.013513,0.077339};

//  double vals_sin_mc[Nbins] = {0.77157,0.751474,0.735441,0.743957,0.723588};
//  double errs_sin_mc[Nbins] = {0.00804551,0.0101668,0.014271,0.00922251,0.0544614};
//  double vals_cos_mc[Nbins] = {0.742978,0.745983,0.711214,0.737736,0.685007};
//  double errs_cos_mc[Nbins] = {0.0117654,0.0148201,0.020899,0.0134883,0.0772523};

// With corrected chi2 value for multiple tracks modes //
  double vals_sin_mc[Nbins] =   {0.771570,0.751474,0.752171,0.752476,0.747389};
  double errs_sin_mc[Nbins] =   {0.008046,0.010167,0.014463,0.009341,0.055953};
  double vals_cos_mc[Nbins] =   {0.742978,0.745983,0.713963,0.746852,0.702061};
  double errs_cos_mc[Nbins] =   {0.011765,0.014820,0.021247,0.013642,0.077865};

  for(int i=0; i<Nbins; i++){
    vals_sin_mc[i] = (fabs(vals_sin_mc[i]) - true_sin)*100.;
    vals_cos_mc[i] = (fabs(vals_cos_mc[i]) - true_cos)*100.;
    errs_sin_mc[i] *= 100.; errs_cos_mc[i] *= 100.;
    cout << vals_sin_mc[i] << " " << vals_cos_mc[i] << endl;
  }

  TCanvas *c_sin = new TCanvas("c_sin","c_sin",400,400);
  c_sin->SetGrid();
  TH1D *h_sin = new TH1D("h_sin","sin offset (10^{-2})",3,0,3);
  h_sin->SetStats(0);
  h_sin->SetFillColor(38);
  h_sin->SetMarkerStyle(20);
  h_sin->SetMarkerSize(1.3);
  h_sin->SetMarkerColor(kRed);
//  h_sin->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin->Fill(cpvModes[i].c_str(),vals_sin[i]);
    h_sin->SetBinError(i+1,errs_sin[i]);
  }
  h_sin->LabelsDeflate();
  h_sin->GetXaxis()->SetLabelSize(0.08);
  h_sin->GetYaxis()->SetLabelSize(0.06);
  h_sin->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_sin_pt = new TH1D("h_sin_pt","sin offset (10^{-2})",3,0,3);
  h_sin_pt->SetStats(0);
  h_sin_pt->SetFillColor(38);
  h_sin_pt->SetMarkerStyle(20);
  h_sin_pt->SetMarkerSize(1.3);
  h_sin_pt->SetMarkerColor(kBlue);
//  h_sin_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin_pt->Fill(cpvModes[i].c_str(),vals_sin_pt[i]);
    h_sin_pt->SetBinError(i+1,errs_sin_pt[i]);
  }

  TH1D *h_sin_pb = new TH1D("h_sin_pb","sin offset (10^{-2})",3,0,3);
  h_sin_pb->SetStats(0);
  h_sin_pb->SetFillColor(38);
  h_sin_pb->SetMarkerStyle(21);
  h_sin_pb->SetMarkerSize(1.3);
  h_sin_pb->SetMarkerColor(kRed);
//  h_sin_pb->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin_pb->Fill(cpvModes[i].c_str(),vals_sin_pb[i]);
    h_sin_pb->SetBinError(i+1,errs_sin_pb[i]);
  }

  TH1D *h_sin_mc = new TH1D("h_sin_mc","sin offset (10^{-2})",3,0,3);
  h_sin_mc->SetStats(0);
  h_sin_mc->SetFillColor(38);
  h_sin_mc->SetMarkerStyle(21);
  h_sin_mc->SetMarkerSize(1.3);
  h_sin_mc->SetMarkerColor(kRed);
//  h_sin_mc->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_sin_mc->Fill(cpvModes[i].c_str(),vals_sin_mc[i]);
    h_sin_mc->SetBinError(i+1,errs_sin_mc[i]);
  }
  h_sin_mc->LabelsDeflate();
  h_sin_mc->GetXaxis()->SetLabelSize(0.08);
  h_sin_mc->GetYaxis()->SetLabelSize(0.06);
  h_sin_mc->GetYaxis()->SetRangeUser(offmin,offmax);

//  h_sin->Draw("e");
  h_sin_mc->Draw("e");
  h_sin_pt->Draw("same");
//  h_sin_pb->Draw("same");

  c_sin->Print("pics/sin_offset_nobkg_corrected.eps");

  TCanvas *c_cos = new TCanvas("c_cos","c_cos",400,400);
  c_cos->SetGrid();
  TH1D *h_cos = new TH1D("h_cos","cos offset (10^{-2})",3,0,3);
  h_cos->SetStats(0);
  h_cos->SetFillColor(38);
  h_cos->SetMarkerStyle(20);
  h_cos->SetMarkerSize(1.3);
  h_cos->SetMarkerColor(kRed);
//  h_cos->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos->Fill(cpvModes[i].c_str(),vals_cos[i]);
    h_cos->SetBinError(i+1,errs_cos[i]);
  }
  h_cos->LabelsDeflate();
  h_cos->GetXaxis()->SetLabelSize(0.08);
  h_cos->GetYaxis()->SetLabelSize(0.06);
  h_cos->GetYaxis()->SetRangeUser(offmin,offmax);

  TH1D *h_cos_pt = new TH1D("h_cos_pt","cos offset (10^{-2})",3,0,3);
  h_cos_pt->SetStats(0);
  h_cos_pt->SetFillColor(38);
  h_cos_pt->SetMarkerStyle(20);
  h_cos_pt->SetMarkerSize(1.3);
  h_cos_pt->SetMarkerColor(kBlue);
//  h_cos_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos_pt->Fill(cpvModes[i].c_str(),vals_cos_pt[i]);
    h_cos_pt->SetBinError(i+1,errs_cos_pt[i]);
  }

  TH1D *h_cos_pb = new TH1D("h_cos_pb","cos offset (10^{-2})",3,0,3);
  h_cos_pb->SetStats(0);
  h_cos_pb->SetFillColor(38);
  h_cos_pb->SetMarkerStyle(21);
  h_cos_pb->SetMarkerSize(1.3);
  h_cos_pb->SetMarkerColor(kRed);
//  h_cos_pb->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos_pb->Fill(cpvModes[i].c_str(),vals_cos_pb[i]);
    h_cos_pb->SetBinError(i+1,errs_cos_pb[i]);
  }

  TH1D *h_cos_mc = new TH1D("h_cos_mc","cos offset (10^{-2})",3,0,3);
  h_cos_mc->SetStats(0);
  h_cos_mc->SetFillColor(38);
  h_cos_mc->SetMarkerStyle(21);
  h_cos_mc->SetMarkerSize(1.3);
  h_cos_mc->SetMarkerColor(kRed);
//  h_cos_mc->SetBit(TH1::kCanRebin);
  for(int i=0; i<Nbins; i++){
    h_cos_mc->Fill(cpvModes[i].c_str(),vals_cos_mc[i]);
    h_cos_mc->SetBinError(i+1,errs_cos_mc[i]);
  }
  h_cos_mc->LabelsDeflate();
  h_cos_mc->GetXaxis()->SetLabelSize(0.08);
  h_cos_mc->GetYaxis()->SetLabelSize(0.06);
  h_cos_mc->GetYaxis()->SetRangeUser(offmin,offmax);

//  h_cos->Draw("e");
//  h_cos_pb->Draw("same");
  h_cos_mc->Draw("e");
  h_cos_pt->Draw("same");
  c_cos->Print("pics/cos_offset_nobkg_corrected.eps");

  double ave_cos_bin_offset = 0;
  double ave_sin_bin_offset = 0;
  double max_cos_bin_offset = 0;
  double max_sin_bin_offset = 0;
  for(int i=0; i<Nbins; i++){
    ave_cos_bin_offset += fabs(vals_cos_pt[i]-vals_cos_mc[i]);
    ave_sin_bin_offset += fabs(vals_sin_pt[i]-vals_sin_mc[i]);
    if(max_cos_bin_offset < fabs(vals_cos_pt[i]-vals_cos_mc[i])) max_cos_bin_offset = fabs(vals_cos_pt[i]-vals_cos_mc[i]);
    if(max_sin_bin_offset < fabs(vals_sin_pt[i]-vals_sin_mc[i])) max_sin_bin_offset = fabs(vals_sin_pt[i]-vals_sin_mc[i]);
  }
  ave_cos_bin_offset /= Nbins;
  ave_sin_bin_offset /= Nbins;

  double ave_cos_offset_mlt = vals_cos_mc[2]/errs_cos_mc[2] + vals_cos_mc[3]/errs_cos_mc[3] + vals_cos_mc[4]/errs_cos_mc[4];
  double ave_sin_offset_mlt = vals_sin_mc[2]/errs_sin_mc[2] + vals_sin_mc[3]/errs_sin_mc[3] + vals_sin_mc[4]/errs_sin_mc[4];
  double err_cos_offset_mlt = 1./errs_cos_mc[2] + 1./errs_cos_mc[3] + 1./errs_cos_mc[4];
  double err_sin_offset_mlt = 1./errs_sin_mc[2] + 1./errs_sin_mc[3] + 1./errs_sin_mc[4];
  ave_cos_offset_mlt /= err_cos_offset_mlt;
  ave_sin_offset_mlt /= err_sin_offset_mlt;
  err_cos_offset_mlt = 1./err_cos_offset_mlt;
  err_sin_offset_mlt = 1./err_sin_offset_mlt;

  double ave_cos_offset_sgl = vals_cos_mc[0]/errs_cos_mc[0] + vals_cos_mc[1]/errs_cos_mc[1];
  double ave_sin_offset_sgl = vals_sin_mc[0]/errs_sin_mc[0] + vals_sin_mc[1]/errs_sin_mc[1];
  double err_cos_offset_sgl = 1./errs_cos_mc[0] + 1./errs_cos_mc[0];
  double err_sin_offset_sgl = 1./errs_sin_mc[0] + 1./errs_sin_mc[0];
  ave_cos_offset_sgl /= err_cos_offset_sgl;
  ave_sin_offset_sgl /= err_sin_offset_sgl;
  err_cos_offset_sgl = 1./err_cos_offset_sgl;
  err_sin_offset_sgl = 1./err_sin_offset_sgl;

  cout << "Offsets summary:" << endl;
  cout << "** Signal resolution offset **" << endl;
  cout << "  cos: " << ave_cos_offset_sgl << " +- " << err_cos_offset_sgl << " - single, ";
  cout << ave_cos_offset_mlt << " +- " << err_cos_offset_mlt << " - multiple" << endl;
  cout << "  sin: " << ave_sin_offset_sgl << " +- " << err_sin_offset_sgl << " - single, ";
  cout << ave_sin_offset_mlt << " +- " << err_sin_offset_mlt << " - multiple" << endl;
  cout << "** Dalitz bin offset **" << endl;
  cout << "  cos: " << ave_cos_bin_offset << " - average, " << max_cos_bin_offset << " - max" << endl;
  cout << "  sin: " << ave_sin_bin_offset << " - average, " << max_sin_bin_offset << " - max" << endl;
  }
  return;
}
