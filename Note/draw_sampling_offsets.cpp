void draw_sampling_offsets(void){
  // * Sampled offsets * //
  // sin cos sampled pi0
  const int FourBins = 4;
  char* FourModes[FourBins] = {"#pi^{0}","#eta#rightarrow#gamma#gamma","#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}","#omega"};

  // 5 //
  double sin_smpl_tagv_ebeb[FourBins]     = {4.1,5.7,8.9,6.1};
  double sin_smpl_tagv_ebeb_err[FourBins] = {1.5,2.5,4.7,2.0};
  double sin_smpl_tagv_ebeb_rms[FourBins] = {27,31,30,30};
  double cos_smpl_tagv_ebeb[FourBins]     = {4.0,0.3,-2.8,3.7};
  double cos_smpl_tagv_ebeb_err[FourBins] = {2.3,3.2,6.0,2.8};
  double cos_smpl_tagv_ebeb_rms[FourBins] = {43,41,38,42};

  // 4 //
  double sin_smpl_tagv[FourBins]     = {5.0,5.5,12.2,6.1};
  double sin_smpl_tagv_err[FourBins] = {1.5,2.4,5.1,1.9};
  double sin_smpl_tagv_rms[FourBins] = {27,39,30,28};
  double cos_smpl_tagv[FourBins]     = {4.6,-3.0,-0.4,4.4};
  double cos_smpl_tagv_err[FourBins] = {2.3,3.2,5.8,2.9};
  double cos_smpl_tagv_rms[FourBins] = {42,41,34,43};

  // 3 //
  double sin_smpl_tagv_nb[FourBins]     = {4.0,6.8,10.6,1.9};
  double sin_smpl_tagv_nb_err[FourBins] = {1.4,2.2,3.7,1.6};
  double sin_smpl_tagv_nb_rms[FourBins] = {26,28,25,25};
  double cos_smpl_tagv_nb[FourBins]     = {2.1,1.4,1.1,1.9};
  double cos_smpl_tagv_nb_err[FourBins] = {2.2,2.9,4.9,2.5};
  double cos_smpl_tagv_nb_rms[FourBins] = {41,37,33,37};

  // 2 //
  double sin_smpl_pt[FourBins]     = {0.9,1.9,2.5,2.8};
  double sin_smpl_pt_err[FourBins] = {0.9,1.3,2.7,1.1};
  double sin_smpl_pt_rms[FourBins] = {16,17,18,16};
  double cos_smpl_pt[FourBins]     = {1.6,-1.4,-5.9,0.5};
  double cos_smpl_pt_err[FourBins] = {1.3,1.8,3.0,1.5};
  double cos_smpl_pt_rms[FourBins] = {23,22,19,22};

  // 1 //
  double sin_smpl_pt_nb[FourBins]     = {0.0,1.4,2.4,-0.9};
  double sin_smpl_pt_nb_err[FourBins] = {0.8,1.2,2.1,0.9};
  double sin_smpl_pt_nb_rms[FourBins] = {15,15,14,13};
  double cos_smpl_pt_nb[FourBins]     = {0.0,-0.3,-3.1,-1.1};
  double cos_smpl_pt_nb_err[FourBins] = {1.3,1.6,2.2,1.3};
  double cos_smpl_pt_nb_rms[FourBins] = {22,20,15,19};

  // ** sin ** //
  const double SIN_MIN = -15;
  const double SIN_MAX =  15;
    TCanvas *c_sin_sampled_offset_600 = new TCanvas("c_sin_sampled_offset_600","c_sin_sampled_offset_600",400,400);
  c_sin_sampled_offset_600->SetGrid();
  // 5 //
  TH1D *h_sin_smpl_tagv_ebeb = new TH1D("h_sin_smpl_offset_600_pt_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_tagv_ebeb->SetStats(0);
  h_sin_smpl_tagv_ebeb->SetFillColor(38);
  h_sin_smpl_tagv_ebeb->SetMarkerStyle(21);
  h_sin_smpl_tagv_ebeb->SetMarkerSize(1.3);
  h_sin_smpl_tagv_ebeb->SetMarkerColor(kBlack);
  h_sin_smpl_tagv_ebeb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_tagv_ebeb->Fill(FourModes[i],sin_smpl_tagv_ebeb[i]);
    h_sin_smpl_tagv_ebeb->SetBinError(i+1,sin_smpl_tagv_ebeb_err[i]);
  }
  h_sin_smpl_tagv_ebeb->LabelsDeflate();
  h_sin_smpl_tagv_ebeb->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_tagv_ebeb->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_tagv_ebeb->GetYaxis()->SetRangeUser(SIN_MIN,SIN_MAX);

  // 4 //
  TH1D *h_sin_smpl_tagv = new TH1D("h_sin_smpl_offset_600_pt_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_tagv->SetStats(0);
  h_sin_smpl_tagv->SetFillColor(38);
  h_sin_smpl_tagv->SetMarkerStyle(21);
  h_sin_smpl_tagv->SetMarkerSize(1.3);
  h_sin_smpl_tagv->SetMarkerColor(kRed);
  h_sin_smpl_tagv->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_tagv->Fill(FourModes[i],sin_smpl_tagv[i]);
    h_sin_smpl_tagv->SetBinError(i+1,sin_smpl_tagv_err[i]);
  }
  h_sin_smpl_tagv->LabelsDeflate();
  h_sin_smpl_tagv->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_tagv->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_tagv->GetYaxis()->SetRangeUser(SIN_MIN,SIN_MAX);

  // 3 //
  TH1D *h_sin_smpl_tagv_nb = new TH1D("h_sin_smpl_offset_600_pt_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_tagv_nb->SetStats(0);
  h_sin_smpl_tagv_nb->SetFillColor(38);
  h_sin_smpl_tagv_nb->SetMarkerStyle(21);
  h_sin_smpl_tagv_nb->SetMarkerSize(1.3);
  h_sin_smpl_tagv_nb->SetMarkerColor(kBlue);
  h_sin_smpl_tagv_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_tagv_nb->Fill(FourModes[i],sin_smpl_tagv_nb[i]);
    h_sin_smpl_tagv_nb->SetBinError(i+1,sin_smpl_tagv_nb_err[i]);
  }
  h_sin_smpl_tagv_nb->LabelsDeflate();
  h_sin_smpl_tagv_nb->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_tagv_nb->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_tagv_nb->GetYaxis()->SetRangeUser(SIN_MIN,SIN_MAX);

  // 2 //
  TH1D *h_sin_smpl_pt = new TH1D("h_sin_smpl_offset_600_pt_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_pt->SetStats(0);
  h_sin_smpl_pt->SetFillColor(38);
  h_sin_smpl_pt->SetMarkerStyle(20);
  h_sin_smpl_pt->SetMarkerSize(1.3);
  h_sin_smpl_pt->SetMarkerColor(kRed);
  h_sin_smpl_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_pt->Fill(FourModes[i],sin_smpl_pt[i]);
    h_sin_smpl_pt->SetBinError(i+1,sin_smpl_pt_err[i]);
  }
  h_sin_smpl_pt->LabelsDeflate();
  h_sin_smpl_pt->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_pt->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_pt->GetYaxis()->SetRangeUser(SIN_MIN,SIN_MAX);

  // 1 //
  TH1D *h_sin_smpl_pt_nb = new TH1D("h_sin_smpl_offset_600_pt_nb_nb","sin offset (10^{-2})",3,0,3);
  h_sin_smpl_pt_nb->SetStats(0);
  h_sin_smpl_pt_nb->SetFillColor(38);
  h_sin_smpl_pt_nb->SetMarkerStyle(20);
  h_sin_smpl_pt_nb->SetMarkerSize(1.3);
  h_sin_smpl_pt_nb->SetMarkerColor(kBlue);
  h_sin_smpl_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_smpl_pt_nb->Fill(FourModes[i],sin_smpl_pt_nb[i]);
    h_sin_smpl_pt_nb->SetBinError(i+1,sin_smpl_pt_nb_err[i]);
  }
  h_sin_smpl_pt_nb->LabelsDeflate();
  h_sin_smpl_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_sin_smpl_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_sin_smpl_pt_nb->GetYaxis()->SetRangeUser(SIN_MIN,SIN_MAX);

  TH1D* sin_hists[5];
  sin_hists[0] = h_sin_smpl_pt_nb;
  sin_hists[1] = h_sin_smpl_pt;
  sin_hists[2] = h_sin_smpl_tagv_nb;
  sin_hists[3] = h_sin_smpl_tagv;
  sin_hists[4] = h_sin_smpl_tagv_ebeb;

  stringstream out;
  out.str("");
//  out << "pics/sin_sampled_offset_0.eps";
//  sin_hists[0]->Draw("e");
  out << "pics/sin_sampled_offset.eps";
  sin_hists[4]->Draw("e");
  c_sin_sampled_offset_600->Print(out.str().c_str());
//  for(int i=1; i<5; i++){
//    out.str("");
//    out << "pics/sin_sampled_offset_" << i << ".eps";
//    sin_hists[i]->Draw("same");
//    c_sin_sampled_offset_600->Print(out.str().c_str());
//  }

  // * RMS * //
  TCanvas *c_sin_sampled_rms_600 = new TCanvas("c_sin_sampled_rms_600","c_sin_sampled_rms_600",400,400);
  c_sin_sampled_rms_600->SetGrid();
  const double SIN_RMS_MIN = 0;
  const double SIN_RMS_MAX = 50;
  // 5 //
  TH1D *h_sin_rms_smpl_tagv_ebeb = new TH1D("h_sin_rms_smpl_offset_600_pt_nb","sin RMS (10^{-2})",3,0,3);
  h_sin_rms_smpl_tagv_ebeb->SetStats(0);
  h_sin_rms_smpl_tagv_ebeb->SetFillColor(38);
  h_sin_rms_smpl_tagv_ebeb->SetMarkerStyle(21);
  h_sin_rms_smpl_tagv_ebeb->SetMarkerSize(1.3);
  h_sin_rms_smpl_tagv_ebeb->SetMarkerColor(kBlack);
  h_sin_rms_smpl_tagv_ebeb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_rms_smpl_tagv_ebeb->Fill(FourModes[i],sin_smpl_tagv_ebeb_rms[i]);
    h_sin_rms_smpl_tagv_ebeb->SetBinError(i+1,0.1);
  }
  h_sin_rms_smpl_tagv_ebeb->LabelsDeflate();
  h_sin_rms_smpl_tagv_ebeb->GetXaxis()->SetLabelSize(0.08);
  h_sin_rms_smpl_tagv_ebeb->GetYaxis()->SetLabelSize(0.06);
  h_sin_rms_smpl_tagv_ebeb->GetYaxis()->SetRangeUser(SIN_RMS_MIN,SIN_RMS_MAX);

  // 4 //
  TH1D *h_sin_rms_smpl_tagv = new TH1D("h_sin_rms_smpl_offset_600_pt_nb","sin RMS (10^{-2})",3,0,3);
  h_sin_rms_smpl_tagv->SetStats(0);
  h_sin_rms_smpl_tagv->SetFillColor(38);
  h_sin_rms_smpl_tagv->SetMarkerStyle(21);
  h_sin_rms_smpl_tagv->SetMarkerSize(1.3);
  h_sin_rms_smpl_tagv->SetMarkerColor(kRed);
  h_sin_rms_smpl_tagv->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_rms_smpl_tagv->Fill(FourModes[i],sin_smpl_tagv_rms[i]);
    h_sin_rms_smpl_tagv->SetBinError(i+1,0.1);
  }
  h_sin_rms_smpl_tagv->LabelsDeflate();
  h_sin_rms_smpl_tagv->GetXaxis()->SetLabelSize(0.08);
  h_sin_rms_smpl_tagv->GetYaxis()->SetLabelSize(0.06);
  h_sin_rms_smpl_tagv->GetYaxis()->SetRangeUser(SIN_RMS_MIN,SIN_RMS_MAX);

  // 3 //
  TH1D *h_sin_rms_smpl_tagv_nb = new TH1D("h_sin_rms_smpl_offset_600_pt_nb","sin RMS (10^{-2})",3,0,3);
  h_sin_rms_smpl_tagv_nb->SetStats(0);
  h_sin_rms_smpl_tagv_nb->SetFillColor(38);
  h_sin_rms_smpl_tagv_nb->SetMarkerStyle(21);
  h_sin_rms_smpl_tagv_nb->SetMarkerSize(1.3);
  h_sin_rms_smpl_tagv_nb->SetMarkerColor(kBlue);
  h_sin_rms_smpl_tagv_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_rms_smpl_tagv_nb->Fill(FourModes[i],sin_smpl_tagv_nb_rms[i]);
    h_sin_rms_smpl_tagv_nb->SetBinError(i+1,0.1);
  }
  h_sin_rms_smpl_tagv_nb->LabelsDeflate();
  h_sin_rms_smpl_tagv_nb->GetXaxis()->SetLabelSize(0.08);
  h_sin_rms_smpl_tagv_nb->GetYaxis()->SetLabelSize(0.06);
  h_sin_rms_smpl_tagv_nb->GetYaxis()->SetRangeUser(SIN_RMS_MIN,SIN_RMS_MAX);

  // 2 //
  TH1D *h_sin_rms_smpl_pt = new TH1D("h_sin_rms_smpl_offset_600_pt_nb","sin RMS (10^{-2})",3,0,3);
  h_sin_rms_smpl_pt->SetStats(0);
  h_sin_rms_smpl_pt->SetFillColor(38);
  h_sin_rms_smpl_pt->SetMarkerStyle(20);
  h_sin_rms_smpl_pt->SetMarkerSize(1.3);
  h_sin_rms_smpl_pt->SetMarkerColor(kRed);
  h_sin_rms_smpl_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_rms_smpl_pt->Fill(FourModes[i],sin_smpl_pt_rms[i]);
    h_sin_rms_smpl_pt->SetBinError(i+1,0.1);
  }
  h_sin_rms_smpl_pt->LabelsDeflate();
  h_sin_rms_smpl_pt->GetXaxis()->SetLabelSize(0.08);
  h_sin_rms_smpl_pt->GetYaxis()->SetLabelSize(0.06);
  h_sin_rms_smpl_pt->GetYaxis()->SetRangeUser(SIN_RMS_MIN,SIN_RMS_MAX);

  // 1 //
  TH1D *h_sin_rms_smpl_pt_nb = new TH1D("h_sin_rms_smpl_offset_600_pt_nb_nb","sin RMS (10^{-2})",3,0,3);
  h_sin_rms_smpl_pt_nb->SetStats(0);
  h_sin_rms_smpl_pt_nb->SetFillColor(38);
  h_sin_rms_smpl_pt_nb->SetMarkerStyle(20);
  h_sin_rms_smpl_pt_nb->SetMarkerSize(1.3);
  h_sin_rms_smpl_pt_nb->SetMarkerColor(kBlue);
  h_sin_rms_smpl_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_sin_rms_smpl_pt_nb->Fill(FourModes[i],sin_smpl_pt_nb_rms[i]);
    h_sin_rms_smpl_pt_nb->SetBinError(i+1,0.1);
  }
  h_sin_rms_smpl_pt_nb->LabelsDeflate();
  h_sin_rms_smpl_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_sin_rms_smpl_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_sin_rms_smpl_pt_nb->GetYaxis()->SetRangeUser(SIN_RMS_MIN,SIN_RMS_MAX);

  TH1D* sin_rms_hists[5];
  sin_rms_hists[0] = h_sin_rms_smpl_pt_nb;
  sin_rms_hists[1] = h_sin_rms_smpl_pt;
  sin_rms_hists[2] = h_sin_rms_smpl_tagv_nb;
  sin_rms_hists[3] = h_sin_rms_smpl_tagv;
  sin_rms_hists[4] = h_sin_rms_smpl_tagv_ebeb;

  stringstream out;
  out.str("");
//  out << "pics/sin_sampled_rms_0.eps";
//  sin_rms_hists[0]->Draw("e");
  out << "pics/sin_sampled_rms.eps";
  sin_rms_hists[4]->Draw("e");
  c_sin_sampled_rms_600->Print(out.str().c_str());
//  for(int i=1; i<5; i++){
//    out.str("");
//    out << "pics/sin_sampled_rms_" << i << ".eps";
//    sin_rms_hists[i]->Draw("same");
//    c_sin_sampled_rms_600->Print(out.str().c_str());
//  }

  // ** cos ** //
  const double cos_MIN = -15;
  const double cos_MAX =  15;
    TCanvas *c_cos_sampled_offset_600 = new TCanvas("c_cos_sampled_offset_600","c_cos_sampled_offset_600",400,400);
  c_cos_sampled_offset_600->SetGrid();
  // 5 //
  TH1D *h_cos_smpl_tagv_ebeb = new TH1D("h_cos_smpl_offset_600_pt_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_tagv_ebeb->SetStats(0);
  h_cos_smpl_tagv_ebeb->SetFillColor(38);
  h_cos_smpl_tagv_ebeb->SetMarkerStyle(21);
  h_cos_smpl_tagv_ebeb->SetMarkerSize(1.3);
  h_cos_smpl_tagv_ebeb->SetMarkerColor(kBlack);
  h_cos_smpl_tagv_ebeb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_tagv_ebeb->Fill(FourModes[i],cos_smpl_tagv_ebeb[i]);
    h_cos_smpl_tagv_ebeb->SetBinError(i+1,cos_smpl_tagv_ebeb_err[i]);
  }
  h_cos_smpl_tagv_ebeb->LabelsDeflate();
  h_cos_smpl_tagv_ebeb->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_tagv_ebeb->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_tagv_ebeb->GetYaxis()->SetRangeUser(cos_MIN,cos_MAX);

  // 4 //
  TH1D *h_cos_smpl_tagv = new TH1D("h_cos_smpl_offset_600_pt_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_tagv->SetStats(0);
  h_cos_smpl_tagv->SetFillColor(38);
  h_cos_smpl_tagv->SetMarkerStyle(21);
  h_cos_smpl_tagv->SetMarkerSize(1.3);
  h_cos_smpl_tagv->SetMarkerColor(kRed);
  h_cos_smpl_tagv->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_tagv->Fill(FourModes[i],cos_smpl_tagv[i]);
    h_cos_smpl_tagv->SetBinError(i+1,cos_smpl_tagv_err[i]);
  }
  h_cos_smpl_tagv->LabelsDeflate();
  h_cos_smpl_tagv->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_tagv->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_tagv->GetYaxis()->SetRangeUser(cos_MIN,cos_MAX);

  // 3 //
  TH1D *h_cos_smpl_tagv_nb = new TH1D("h_cos_smpl_offset_600_pt_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_tagv_nb->SetStats(0);
  h_cos_smpl_tagv_nb->SetFillColor(38);
  h_cos_smpl_tagv_nb->SetMarkerStyle(21);
  h_cos_smpl_tagv_nb->SetMarkerSize(1.3);
  h_cos_smpl_tagv_nb->SetMarkerColor(kBlue);
  h_cos_smpl_tagv_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_tagv_nb->Fill(FourModes[i],cos_smpl_tagv_nb[i]);
    h_cos_smpl_tagv_nb->SetBinError(i+1,cos_smpl_tagv_nb_err[i]);
  }
  h_cos_smpl_tagv_nb->LabelsDeflate();
  h_cos_smpl_tagv_nb->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_tagv_nb->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_tagv_nb->GetYaxis()->SetRangeUser(cos_MIN,cos_MAX);

  // 2 //
  TH1D *h_cos_smpl_pt = new TH1D("h_cos_smpl_offset_600_pt_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_pt->SetStats(0);
  h_cos_smpl_pt->SetFillColor(38);
  h_cos_smpl_pt->SetMarkerStyle(20);
  h_cos_smpl_pt->SetMarkerSize(1.3);
  h_cos_smpl_pt->SetMarkerColor(kRed);
  h_cos_smpl_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_pt->Fill(FourModes[i],cos_smpl_pt[i]);
    h_cos_smpl_pt->SetBinError(i+1,cos_smpl_pt_err[i]);
  }
  h_cos_smpl_pt->LabelsDeflate();
  h_cos_smpl_pt->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_pt->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_pt->GetYaxis()->SetRangeUser(cos_MIN,cos_MAX);

  // 1 //
  TH1D *h_cos_smpl_pt_nb = new TH1D("h_cos_smpl_offset_600_pt_nb_nb","cos offset (10^{-2})",3,0,3);
  h_cos_smpl_pt_nb->SetStats(0);
  h_cos_smpl_pt_nb->SetFillColor(38);
  h_cos_smpl_pt_nb->SetMarkerStyle(20);
  h_cos_smpl_pt_nb->SetMarkerSize(1.3);
  h_cos_smpl_pt_nb->SetMarkerColor(kBlue);
  h_cos_smpl_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_smpl_pt_nb->Fill(FourModes[i],cos_smpl_pt_nb[i]);
    h_cos_smpl_pt_nb->SetBinError(i+1,cos_smpl_pt_nb_err[i]);
  }
  h_cos_smpl_pt_nb->LabelsDeflate();
  h_cos_smpl_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_cos_smpl_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_cos_smpl_pt_nb->GetYaxis()->SetRangeUser(cos_MIN,cos_MAX);

  TH1D* cos_hists[5];
  cos_hists[0] = h_cos_smpl_pt_nb;
  cos_hists[1] = h_cos_smpl_pt;
  cos_hists[2] = h_cos_smpl_tagv_nb;
  cos_hists[3] = h_cos_smpl_tagv;
  cos_hists[4] = h_cos_smpl_tagv_ebeb;

  stringstream out;
  out.str("");
//  out << "pics/cos_sampled_offset_0.eps";
//  cos_hists[0]->Draw("e");
  out << "pics/cos_sampled_offset.eps";
  cos_hists[4]->Draw("e");
  c_cos_sampled_offset_600->Print(out.str().c_str());
//  for(int i=1; i<5; i++){
//    out.str("");
//    out << "pics/cos_sampled_offset_" << i << ".eps";
//    cos_hists[i]->Draw("same");
//    c_cos_sampled_offset_600->Print(out.str().c_str());
//  }

  // * RMS * //
  TCanvas *c_cos_sampled_rms_600 = new TCanvas("c_cos_sampled_rms_600","c_cos_sampled_rms_600",400,400);
  c_cos_sampled_rms_600->SetGrid();
  const double cos_RMS_MIN = 0;
  const double cos_RMS_MAX = 50;
  // 5 //
  TH1D *h_cos_rms_smpl_tagv_ebeb = new TH1D("h_cos_rms_smpl_offset_600_pt_nb","cos RMS (10^{-2})",3,0,3);
  h_cos_rms_smpl_tagv_ebeb->SetStats(0);
  h_cos_rms_smpl_tagv_ebeb->SetFillColor(38);
  h_cos_rms_smpl_tagv_ebeb->SetMarkerStyle(21);
  h_cos_rms_smpl_tagv_ebeb->SetMarkerSize(1.3);
  h_cos_rms_smpl_tagv_ebeb->SetMarkerColor(kBlack);
  h_cos_rms_smpl_tagv_ebeb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_rms_smpl_tagv_ebeb->Fill(FourModes[i],cos_smpl_tagv_ebeb_rms[i]);
    h_cos_rms_smpl_tagv_ebeb->SetBinError(i+1,0.1);
  }
  h_cos_rms_smpl_tagv_ebeb->LabelsDeflate();
  h_cos_rms_smpl_tagv_ebeb->GetXaxis()->SetLabelSize(0.08);
  h_cos_rms_smpl_tagv_ebeb->GetYaxis()->SetLabelSize(0.06);
  h_cos_rms_smpl_tagv_ebeb->GetYaxis()->SetRangeUser(cos_RMS_MIN,cos_RMS_MAX);

  // 4 //
  TH1D *h_cos_rms_smpl_tagv = new TH1D("h_cos_rms_smpl_offset_600_pt_nb","cos RMS (10^{-2})",3,0,3);
  h_cos_rms_smpl_tagv->SetStats(0);
  h_cos_rms_smpl_tagv->SetFillColor(38);
  h_cos_rms_smpl_tagv->SetMarkerStyle(21);
  h_cos_rms_smpl_tagv->SetMarkerSize(1.3);
  h_cos_rms_smpl_tagv->SetMarkerColor(kRed);
  h_cos_rms_smpl_tagv->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_rms_smpl_tagv->Fill(FourModes[i],cos_smpl_tagv_rms[i]);
    h_cos_rms_smpl_tagv->SetBinError(i+1,0.1);
  }
  h_cos_rms_smpl_tagv->LabelsDeflate();
  h_cos_rms_smpl_tagv->GetXaxis()->SetLabelSize(0.08);
  h_cos_rms_smpl_tagv->GetYaxis()->SetLabelSize(0.06);
  h_cos_rms_smpl_tagv->GetYaxis()->SetRangeUser(cos_RMS_MIN,cos_RMS_MAX);

  // 3 //
  TH1D *h_cos_rms_smpl_tagv_nb = new TH1D("h_cos_rms_smpl_offset_600_pt_nb","cos RMS (10^{-2})",3,0,3);
  h_cos_rms_smpl_tagv_nb->SetStats(0);
  h_cos_rms_smpl_tagv_nb->SetFillColor(38);
  h_cos_rms_smpl_tagv_nb->SetMarkerStyle(21);
  h_cos_rms_smpl_tagv_nb->SetMarkerSize(1.3);
  h_cos_rms_smpl_tagv_nb->SetMarkerColor(kBlue);
  h_cos_rms_smpl_tagv_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_rms_smpl_tagv_nb->Fill(FourModes[i],cos_smpl_tagv_nb_rms[i]);
    h_cos_rms_smpl_tagv_nb->SetBinError(i+1,0.1);
  }
  h_cos_rms_smpl_tagv_nb->LabelsDeflate();
  h_cos_rms_smpl_tagv_nb->GetXaxis()->SetLabelSize(0.08);
  h_cos_rms_smpl_tagv_nb->GetYaxis()->SetLabelSize(0.06);
  h_cos_rms_smpl_tagv_nb->GetYaxis()->SetRangeUser(cos_RMS_MIN,cos_RMS_MAX);

  // 2 //
  TH1D *h_cos_rms_smpl_pt = new TH1D("h_cos_rms_smpl_offset_600_pt_nb","cos RMS (10^{-2})",3,0,3);
  h_cos_rms_smpl_pt->SetStats(0);
  h_cos_rms_smpl_pt->SetFillColor(38);
  h_cos_rms_smpl_pt->SetMarkerStyle(20);
  h_cos_rms_smpl_pt->SetMarkerSize(1.3);
  h_cos_rms_smpl_pt->SetMarkerColor(kRed);
  h_cos_rms_smpl_pt->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_rms_smpl_pt->Fill(FourModes[i],cos_smpl_pt_rms[i]);
    h_cos_rms_smpl_pt->SetBinError(i+1,0.1);
  }
  h_cos_rms_smpl_pt->LabelsDeflate();
  h_cos_rms_smpl_pt->GetXaxis()->SetLabelSize(0.08);
  h_cos_rms_smpl_pt->GetYaxis()->SetLabelSize(0.06);
  h_cos_rms_smpl_pt->GetYaxis()->SetRangeUser(cos_RMS_MIN,cos_RMS_MAX);

  // 1 //
  TH1D *h_cos_rms_smpl_pt_nb = new TH1D("h_cos_rms_smpl_offset_600_pt_nb_nb","cos RMS (10^{-2})",3,0,3);
  h_cos_rms_smpl_pt_nb->SetStats(0);
  h_cos_rms_smpl_pt_nb->SetFillColor(38);
  h_cos_rms_smpl_pt_nb->SetMarkerStyle(20);
  h_cos_rms_smpl_pt_nb->SetMarkerSize(1.3);
  h_cos_rms_smpl_pt_nb->SetMarkerColor(kBlue);
  h_cos_rms_smpl_pt_nb->SetBit(TH1::kCanRebin);
  for(int i=0; i<FourBins; i++){
    h_cos_rms_smpl_pt_nb->Fill(FourModes[i],cos_smpl_pt_nb_rms[i]);
    h_cos_rms_smpl_pt_nb->SetBinError(i+1,0.1);
  }
  h_cos_rms_smpl_pt_nb->LabelsDeflate();
  h_cos_rms_smpl_pt_nb->GetXaxis()->SetLabelSize(0.08);
  h_cos_rms_smpl_pt_nb->GetYaxis()->SetLabelSize(0.06);
  h_cos_rms_smpl_pt_nb->GetYaxis()->SetRangeUser(cos_RMS_MIN,cos_RMS_MAX);

  TH1D* cos_rms_hists[5];
  cos_rms_hists[0] = h_cos_rms_smpl_pt_nb;
  cos_rms_hists[1] = h_cos_rms_smpl_pt;
  cos_rms_hists[2] = h_cos_rms_smpl_tagv_nb;
  cos_rms_hists[3] = h_cos_rms_smpl_tagv;
  cos_rms_hists[4] = h_cos_rms_smpl_tagv_ebeb;

  stringstream out;
  out.str("");
//  out << "pics/cos_sampled_rms_0.eps";
//  cos_rms_hists[0]->Draw("e");
  out << "pics/cos_sampled_rms.eps";
  cos_rms_hists[4]->Draw("e");
  c_cos_sampled_rms_600->Print(out.str().c_str());
//  for(int i=1; i<5; i++){
//    out.str("");
//    out << "pics/cos_sampled_rms_" << i << ".eps";
//    cos_rms_hists[i]->Draw("same");
//    c_cos_sampled_rms_600->Print(out.str().c_str());
//  }

  return;
}
