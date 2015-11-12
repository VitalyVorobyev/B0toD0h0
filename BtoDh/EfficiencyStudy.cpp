void EfficiencyStudy(void){
  TFile *file = TFile::Open("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcOMEGA_s5.root");
  TTree *tree = (TTree*)file->Get("TEvent");

//  double Carr[8], Sarr[8];
  double Karr[8], Kbarr[8];
  double Karr_err[8], Kbarr_err[8];
  // CLEO measurements
//  const double Carr_CLEO[8] = { 0.365, 0.710, 0.481,0.008,-0.757,-0.884,-0.462, 0.106};
//  const double Sarr_CLEO[8] = {-0.179,-0.013,-0.147,0.938, 0.386,-0.162,-0.616,-1.063};
  double Karr_CLEO[8] = { 0.067, 0.093, 0.030,0.089, 0.079, 0.102, 0.123, 0.159};
  double Kbarr_CLEO[8]= { 0.016, 0.027, 0.013,0.040, 0.015, 0.012, 0.026, 0.109};

  const double Karr_CLEO_err[8] = { 0.004, 0.004, 0.002,0.004, 0.003, 0.004, 0.004, 0.004};
  const double Kbarr_CLEO_err[8]= { 0.002, 0.002, 0.001,0.003, 0.002, 0.002, 0.002, 0.004};
//  const double Carr_CLEO_stat_err[8] = {0.071,0.034,0.080,0.080,0.099,0.056,0.100,0.105};
//  const double Carr_CLEO_sist_err[8] = {0.078,0.038,0.087,0.087,0.065,0.054,0.082,0.100};
//  const double Sarr_CLEO_stat_err[8] = {0.166,0.097,0.177,0.120,0.208,0.130,0.188,0.174};
//  const double Sarr_CLEO_sist_err[8] = {0.048,0.031,0.107,0.047,0.067,0.041,0.052,0.066};

  // Belle model integrals
//  const double Carr_model[8] = {0.554166,-0.00454239,-0.608498,-0.934161,-0.565509,0.0800234,0.480612,0.68653};
//  const double Sarr_model[8] = {0.446559,0.825851,0.678567,-0.000976537,-0.572031,-0.734359,-0.372541,0.0279944};
  const double Karr_model[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
  const double Kbarr_model[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};

  for(int i=0; i<8; i++){
    Karr[i] = 0; Kbarr[i] = 0;
  }
  double mp,mm;
  int flv,bin,b0f;
  tree->SetBranchAddress("mp_mc",&mp);
  tree->SetBranchAddress("mm_mc",&mm);
  tree->SetBranchAddress("flv_mc",&flv);
  tree->SetBranchAddress("bin_mc",&bin);
  tree->SetBranchAddress("b0f",&b0f);

  const int NTot = tree->GetEntries();
  int DalitzBin;
  int Norm = 0;
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    if(bin == 0 || b0f != 1) continue;
    Norm++;
    DalitzBin = bin*flv;
    DalitzBin > 0 ? Kbarr[DalitzBin-1]++ : Karr[-DalitzBin-1]++;
  }
  cout << Norm << " good events in the tree" << endl;
  for(int i=0; i<8; i++){
    Karr_err[i] = sqrt(Karr[i])/Norm;
    Kbarr_err[i] = sqrt(Kbarr[i])/Norm;
    Karr[i] /= Norm; Kbarr[i] /= Norm;
    cout << "K[" << i+1 << "] = " << Karr[i] << " +- " << Karr_err[i];
    cout << ", Kb[" << i+1 << "] = " << Kbarr[i] << " +- " << Kbarr_err[i] << endl;
    Karr[i]      -= Karr_model[i]; Kbarr[i]      -= Kbarr_model[i];
    Karr_CLEO[i] -= Karr_model[i]; Kbarr_CLEO[i] -= Kbarr_model[i];
  }

  double bin_arr[8]     = {1,2,3,4,5,6,7,8};
  double bin_arr_err[8] = {0,0,0,0,0,0,0,0};

  TGraph* grK_model  = new TGraph(8,bin_arr,Karr_model);
  grK_model->SetMarkerStyle(20);
  grK_model->SetMarkerSize(1);
  grK_model->SetLineStyle(0);
  grK_model->SetMarkerColor(kBlue);

  TGraphErrors* grK_CLEO  = new TGraphErrors(8,bin_arr,Karr_CLEO,bin_arr_err,Karr_CLEO_err);
  grK_CLEO->SetMarkerStyle(20);
  grK_CLEO->SetMarkerSize(1);
  grK_CLEO->SetLineStyle(0);
  grK_CLEO->SetMarkerColor(kRed);

  TGraphErrors* grK = new TGraphErrors(8,bin_arr,Karr,bin_arr_err,Karr_err);
  grK->SetMarkerStyle(21);
  grK->SetMarkerSize(1);
  grK->SetLineStyle(0);
  grK->SetMarkerColor(kBlue);


  TMultiGraph* mgK = new TMultiGraph("mgK","mgK");
//  mgK->Add(grK_model);
  mgK->Add(grK_CLEO);
  mgK->Add(grK);

  TCanvas* c1 = new TCanvas("c1","c1",800,900);
  c1->cd();
  mgK->Draw("ap");
  c1->Update();

  TGraph* grKb_model = new TGraph(8,bin_arr,Kbarr_model);
  grKb_model->SetMarkerStyle(20);
  grKb_model->SetMarkerSize(1);
  grKb_model->SetLineStyle(0);
  grKb_model->SetMarkerColor(kBlue);

  TGraphErrors* grKb_CLEO = new TGraphErrors(8,bin_arr,Kbarr_CLEO,bin_arr_err,Kbarr_CLEO_err);
  grKb_CLEO->SetMarkerStyle(20);
  grKb_CLEO->SetMarkerSize(1);
  grKb_CLEO->SetLineStyle(0);
  grKb_CLEO->SetMarkerColor(kRed);

  TGraphErrors* grKb = new TGraphErrors(8,bin_arr,Kbarr,bin_arr_err,Kbarr_err);
  grKb->SetMarkerStyle(21);
  grKb->SetMarkerSize(1);
  grKb->SetLineStyle(0);
  grKb->SetMarkerColor(kBlue);

  TMultiGraph* mgKb = new TMultiGraph("mgKb","mgKb");
//  mgKb->Add(grKb_model);
  mgKb->Add(grKb_CLEO);
  mgKb->Add(grKb);

  TCanvas* c2 = new TCanvas("c2","c2",800,900);
  c2->cd();
  mgKb->Draw("ap");
  c2->Update();

  return;

}
