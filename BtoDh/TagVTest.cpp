int TagBin(const double& tag){
  const double wbins[8] = {0.0,0.1,0.25,0.5,0.625,0.75,0.875,1.001};
  const double atag = fabs(tag);
//  cout << atag << endl;
  for(int i=1; i<8; i++){
    if(atag<wbins[i]) return i;
  }
  cout << "Wrong tag: " << tag << ", (" << atag << ")" << endl;
  return -1;
}

TTree* tree_svd1;
TTree* tree_svd2;

TTree* GetGoodTTree(TTree* tree, const int _mode, const int svd = 2){
  double demin, demax;
  double bdtg_cut=0;
  int mode, h0mode;
  const double mbcmin = 5.272;
  const double mbcmax = 5.289;
  switch(_mode){
  case 1:// D0 pi0
    demin = -0.08;
    demax =  0.07;
    mode = 1;
    h0mode = 10;
    bdtg_cut = 0.910;
    break;
  case 2:// D0 eta->gg
    demin = -0.08;
    demax =  0.07;
    mode = 2;
    h0mode = 10;
    bdtg_cut = 0.690;
    break;
  case 3:// D0 eta->pi+pi-pi0
    demin = -0.05;
    demax =  0.04;
    mode = 2;
    h0mode = 20;
    bdtg_cut = 0.500;
    break;
  case 4:// D0 omega->pi+pi-pi0
    demin = -0.06;
    demax =  0.05;
    mode = 3;
    h0mode = 20;
    bdtg_cut = 0.125;
    break;
  case 5:// D0 omega->pi+pi-pi0 (Ks rho model)
    demin = -0.06;
    demax =  0.05;
    mode = 3;
    h0mode = 20;
    bdtg_cut = 0.125;
    break;
  case 6:// D0 rho
    demin = -0.06;
    demax =  0.05;
    h0mode = 40;
    mode = 4;
    bdtg_cut = 0.;
    break;
  }

  stringstream out;
  out.str("");
  if(svd == 2)  out << "exp>30 && ";
  else          out << "exp<30 && ";
  out << "(b0f == 1 || b0f == 5 || b0f == 10) && ";

  out << "sz_sig>0 && ";
  out << "de<" << demax << " && " << "de>" << demin << " && ";
  out << "mbc<" << mbcmax << " && " << "mbc>" << mbcmin << " && ";
  out << "mode == " << mode << " && h0mode == " << h0mode << " && ";
  out << "chi2_vtx_d0<50 && ";
  out << "tag_LH_err<0.2 &&";
//  out << "bdtg>" << bdtg_cut << " && ";
//  if(sgl_asc) out << "ndf_z_asc==0 && ";
//  if(mlt_asc) out << "ndf_z_asc>0 && ";
  out << "good_icpv == 1 && ";
  out << "(z_sig-z_asc)*0.1*78.48>" << -70 << " && (z_sig-z_asc)*0.1*78.48<" << 70;

  return tree->CopyTree(out.str().c_str());
}

void TagVTest(const int mode = 1){
  TChain* tree = new TChain("TEvent");
  switch(mode){
  case 1:
    cout << "Mode: pi0" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
    break;
  case 2:
    cout << "Mode: eta -> gg" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
    break;
  case 3:
    cout << "Mode: eta -> pi+pi-pi0" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s2_full.root");
    break;
  case 4:
    cout << "Mode: omega" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s5_full.root");
    break;
  case 5:
    cout << "Mode: omega (rho Ks0)" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s3_full.root");
    break;
  case 6:
    cout << "Mode: rho" << endl;
    tree->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigRho_s1_full.root");
    break;
  default:
    cout << "Wrong mode" << endl;
    return -1;
  }

  Double_t tag;
  Int_t flv;
  tree->SetBranchAddress("tag_LH",&tag);
  tree->SetBranchAddress("flv_mc",&flv);
  const double w_mc_svd1[7]  = {0.5,0.420827,0.300296,0.219317,0.154636,0.0916131,0.0228891};
  const double w_mc_svd2[7]  = {0.5,0.412222,0.307838,0.212765,0.149933,0.0913264,0.0218754};
  const double dw_mc_svd1[7] = {0.0,0.0583019, 0.00573998,-0.0392635, 0.00474508,-0.0118737,-0.00585326};
  const double dw_mc_svd2[7] = {0.0,0.00408778,0.010326,  -0.00479522,0.00151989, 0.0143633, 0.00189979};

  tree_svd1 = GetGoodTTree(tree,mode,1);
  tree_svd2 = GetGoodTTree(tree,mode,2);
  const int NTotSVD1 = tree_svd1->GetEntries();
  const int NTotSVD2 = tree_svd2->GetEntries();
  cout << "SVD1: " << NTotSVD1 << ", SVD2: " << NTotSVD2 << endl;

  Double_t w_svd1_diff[7],w_svd1[7],w_svd2_diff[7],w_svd2[7],wb_svd1[7],wb_svd2[7],dw_svd1[7],dw_svd2[7];
  Double_t w_svd1_err[7],w_svd2_err[7],wb_svd1_err[7],wb_svd2_err[7],dw_svd1_err[7],dw_svd2_err[7];
  double t_svd1[7],t_svd2[7],tb_svd1[7],tb_svd2[7];
  for(int i=0; i<7; i++){
    w_svd1[i] = 0; wb_svd1[i] = 0; dw_svd1[i] = 0; t_svd1[i] = 0; tb_svd1[i] = 0;
    w_svd2[i] = 0; wb_svd2[i] = 0; dw_svd2[i] = 0; t_svd2[i] = 0; tb_svd2[i] = 0;
    w_svd1_diff[i] = 0; w_svd2_diff[i] = 0;
    w_svd1_err[i] = 0; wb_svd1_err[i] = 0; dw_svd1_err[i] = 0;
    w_svd2_err[i] = 0; wb_svd2_err[i] = 0; dw_svd2_err[i] = 0;
  }

  for(int i=0; i<NTotSVD1; i++){
    tree_svd1->GetEvent(i);
//    cout << tag << " " << flv << " " << TagBin(tag) << endl;
    if(!flv) continue;
    const int bin = TagBin(tag)-1;
    if(flv>0){
      if(tag*flv>0) w_svd1[bin] += 1;
      t_svd1[bin] += 1;
    } else{
      if(tag*flv>0) wb_svd1[bin] += 1;
      tb_svd1[bin] += 1;
    }
  }
  for(int i=0; i<NTotSVD2; i++){
    tree_svd2->GetEvent(i);
    if(!flv) continue;
    const int bin = TagBin(tag)-1;
    if(flv>0){
      if(tag*flv>0) w_svd2[bin] += 1;
      t_svd2[bin] += 1;
    } else{
      if(tag*flv>0) wb_svd2[bin] += 1;
      tb_svd2[bin] += 1;
    }
  }

  for(int i=0; i<7; i++){
    w_svd1[i] /= t_svd1[i]; wb_svd1[i] /= tb_svd1[i];
    w_svd2[i] /= t_svd2[i]; wb_svd2[i] /= tb_svd2[i];
    dw_svd1[i] = (w_svd1[i] - wb_svd1[i]);
    dw_svd2[i] = (w_svd2[i] - wb_svd2[i]);
    w_svd1_err[i]  = (1-0.5*0.6827)*sqrt(w_svd1[i]*(1-w_svd1[i])/t_svd1[i]);
    wb_svd1_err[i] = (1-0.5*0.6827)*sqrt(wb_svd1[i]*(1-wb_svd1[i])/tb_svd1[i]);
    w_svd2_err[i]  = (1-0.5*0.6827)*sqrt(w_svd2[i]*(1-w_svd2[i])/t_svd2[i]);
    wb_svd2_err[i] = (1-0.5*0.6827)*sqrt(wb_svd2[i]*(1-wb_svd2[i])/tb_svd2[i]);
    dw_svd1_err[i] = sqrt(w_svd1_err[i]*w_svd1_err[i]+wb_svd1_err[i]*wb_svd1_err[i]);
    dw_svd1_err[i] = sqrt(w_svd2_err[i]*w_svd2_err[i]+wb_svd2_err[i]*wb_svd2_err[i]);

    w_svd1[i] = 0.5*(w_svd1[i]+wb_svd1[i]);
    w_svd2[i] = 0.5*(w_svd2[i]+wb_svd2[i]);
    w_svd1_err[i] = 0.5*sqrt(w_svd1_err[i]*w_svd1_err[i]+wb_svd1_err[i]*wb_svd1_err[i]);
    w_svd2_err[i] = 0.5*sqrt(w_svd2_err[i]*w_svd2_err[i]+wb_svd2_err[i]*wb_svd2_err[i]);

    w_svd1_diff[i] = w_svd1[i] - w_mc_svd1[i];
    w_svd2_diff[i] = w_svd2[i] - w_mc_svd2[i];
  }

  cout << "SVD1:" << endl;
  for(int i=0; i<7; i++) cout << "  w[" << i+1 << "] = " << w_svd1[i] << " +- " << w_svd1_err[i] << " (" << w_mc_svd1[i] << ")" << endl;
  cout << endl;
  for(int i=0; i<7; i++) cout << " dw[" << i+1 << "] = " << dw_svd1[i] << " +- " << dw_svd1_err[i] << " (" << dw_mc_svd1[i] << ")" << endl;
  cout << endl << "SVD2:" << endl;
  for(int i=0; i<7; i++) cout << "  w[" << i+1 << "] = " << w_svd2[i] << " +- " << w_svd2_err[i] << " (" << w_mc_svd2[i] << ")" << endl;
  cout << endl;
  for(int i=0; i<7; i++) cout << " dw[" << i+1 << "] = " << dw_svd2[i] << " +- " << dw_svd2_err[i] << " (" << dw_mc_svd2[i] << ")" << endl;

  cout << "SVD1 anti-B0:" << endl;
  for(int i=0; i<7; i++) cout << "  w[" << i+1 << "] = " << wb_svd1[i] << " +- " << wb_svd1_err[i] << " (" << w_mc_svd1[i] << ")" << endl;
  cout << endl << "SVD2 anti-B0:" << endl;
  for(int i=0; i<7; i++) cout << "  w[" << i+1 << "] = " << wb_svd2[i] << " +- " << wb_svd2_err[i] << " (" << w_mc_svd2[i] << ")" << endl;

  const int NBins = 7;
  const double bins[7]     = {1,2,3,4,5,6,7};
  const double bins_err[7] = {0,0,0,0,0,0,0};
  TMultiGraph* mg1 = new TMultiGraph("mg1","w, SVD1");
  TGraph* g_w_mc1 = new TGraph(NBins,bins,w_mc_svd1);
  g_w_mc1->SetMarkerStyle(20);
  g_w_mc1->SetMarkerColor(kRed);
  g_w_mc1->SetMarkerSize(1.3);
  TGraphErrors* g_w1 = new TGraphErrors(NBins,bins,w_svd1,bins_err,w_svd1_err);
  g_w1->SetMarkerStyle(20);
  g_w1->SetMarkerColor(kBlue);
  g_w1->SetMarkerSize(1.3);
  mg1->Add(g_w_mc1);
  mg1->Add(g_w1);

  TGraphErrors* g_w1_diff = new TGraphErrors(NBins,bins,w_svd1_diff,bins_err,w_svd1_err);
  g_w1_diff->SetTitle("w difference, SVD1");
//  g_w1_diff->SetTitleSize(20);
  g_w1_diff->SetMarkerStyle(20);
  g_w1_diff->SetMarkerColor(kBlue);
  g_w1_diff->SetMarkerSize(1.3);

  TMultiGraph* mg2 = new TMultiGraph("mg2","w, SVD2");
  TGraph* g_w_mc2 = new TGraph(NBins,bins,w_mc_svd2);
  g_w_mc2->SetMarkerStyle(20);
  g_w_mc2->SetMarkerColor(kRed);
  g_w_mc2->SetMarkerSize(1.3);
  TGraphErrors* g_w2 = new TGraphErrors(NBins,bins,w_svd2,bins_err,w_svd2_err);
  g_w2->SetMarkerStyle(20);
  g_w2->SetMarkerColor(kBlue);
  g_w2->SetMarkerSize(1.3);
  mg2->Add(g_w_mc2);
  mg2->Add(g_w2);

  TGraphErrors* g_w2_diff = new TGraphErrors(NBins,bins,w_svd2_diff,bins_err,w_svd2_err);
  g_w2_diff->SetTitle("w difference, SVD2");
//  g_w2_diff->SetTitleSize(20);
  g_w2_diff->SetMarkerStyle(20);
  g_w2_diff->SetMarkerColor(kBlue);
  g_w2_diff->SetMarkerSize(1.3);

  TMultiGraph* mgd1 = new TMultiGraph("mgd1","#Deltaw, SVD1");
//  mgd1->SetTitleSize(20);
  TGraph* g_dw_mc1 = new TGraph(NBins,bins,dw_mc_svd1);
  g_dw_mc1->SetMarkerStyle(20);
  g_dw_mc1->SetMarkerColor(kRed);
  g_dw_mc1->SetMarkerSize(1.3);
  TGraphErrors* g_dw1 = new TGraphErrors(NBins,bins,dw_svd1,bins_err,dw_svd1_err);
  g_dw1->SetMarkerStyle(20);
  g_dw1->SetMarkerColor(kBlue);
  g_dw1->SetMarkerSize(1.3);
  mgd1->Add(g_dw_mc1);
  mgd1->Add(g_dw1);

  TMultiGraph* mgd2 = new TMultiGraph("mgd2","#Deltaw, SVD2");
  TGraph* g_dw_mc2 = new TGraph(NBins,bins,dw_mc_svd2);
  g_dw_mc2->SetMarkerStyle(20);
  g_dw_mc2->SetMarkerColor(kRed);
  g_dw_mc2->SetMarkerSize(1.3);
  TGraphErrors* g_dw2 = new TGraphErrors(NBins,bins,dw_svd2,bins_err,dw_svd2_err);
  g_dw2->SetMarkerStyle(20);
  g_dw2->SetMarkerColor(kBlue);
  g_dw2->SetMarkerSize(1.3);
  mgd2->Add(g_dw_mc2);
  mgd2->Add(g_dw2);

  TCanvas* c1 = new TCanvas("c1","c1",800,300);
  c1->Draw();
  c1->Divide(2,1);

  c1->cd(1);
  c1->Pad()->SetGrid();
  mg1->Draw("ap");
  mg1->GetXaxis()->SetLabelSize(0.06);
  mg1->GetYaxis()->SetLabelSize(0.06);
  gPad->Modified();

  c1->cd(2);
  c1->Pad()->SetGrid();
  mg2->Draw("ap");
  mg2->GetXaxis()->SetLabelSize(0.06);
  mg2->GetYaxis()->SetLabelSize(0.06);
  gPad->Modified();

  TCanvas* c2 = new TCanvas("c2","c2",800,300);
  c2->Draw();
  c2->Divide(2,1);
  c2->cd(1);
  c2->Pad()->SetGrid();
  g_w1_diff->Draw("ap");
  g_w1_diff->GetXaxis()->SetLabelSize(0.06);
  g_w1_diff->GetYaxis()->SetLabelSize(0.06);
  TLine* zeroline = new TLine(1,0,7,0);
  zeroline->SetLineWidth(2);
  zeroline->Draw();
  gPad->Modified();

  c2->cd(2);
  c2->Pad()->SetGrid();
  g_w2_diff->Draw("ap");
  g_w2_diff->GetXaxis()->SetLabelSize(0.06);
  g_w2_diff->GetYaxis()->SetLabelSize(0.06);
  zeroline->Draw();
  gPad->Modified();

  TCanvas* c3 = new TCanvas("c3","c3",800,300);
  c3->Draw();
  c3->Divide(2,1);
  c3->cd(1);
  c3->Pad()->SetGrid();
  mgd1->Draw("ap");
  mgd1->GetXaxis()->SetLabelSize(0.06);
  mgd1->GetYaxis()->SetLabelSize(0.06);
  gPad->Modified();

  c3->cd(2);
  c3->Pad()->SetGrid();
  mgd2->Draw("ap");
  mgd2->GetXaxis()->SetLabelSize(0.06);
  mgd2->GetYaxis()->SetLabelSize(0.06);
  gPad->Modified();

  c1->Update();
  stringstream out;
  out.str("");
  out << "../Note/pics/TagVTest1_m" << mode << ".eps";
  c1->Print(out.str().c_str());
  out.str("");
  out << "../Note/pics/TagVTest1_m" << mode << ".root";
  c1->Print(out.str().c_str());

  c2->Update();
  stringstream out;
  out.str("");
  out << "../Note/pics/TagVTest2_m" << mode << ".eps";
  c2->Print(out.str().c_str());
  out.str("");
  out << "../Note/pics/TagVTest2_m" << mode << ".root";
  c2->Print(out.str().c_str());

  c3->Update();
  stringstream out;
  out.str("");
  out << "../Note/pics/TagVTest3_m" << mode << ".eps";
  c3->Print(out.str().c_str());
  out.str("");
  out << "../Note/pics/TagVTest3_m" << mode << ".root";
  c3->Print(out.str().c_str());
  return;
}
