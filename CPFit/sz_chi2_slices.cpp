
void GetMoments(const int N, const int* nev, const double* sum, const double* sumsq, const double* sumcu, double* mean, double* mean_err, double* sig, double* sk, double* sk_err){
  for(int i=0; i<N; i++){
    if(nev[i]>1){
      mean[i] = sum[i]/nev[i];
      sig[i]  = sqrt((sumsq[i]/nev[i] - mean[i]*mean[i])*nev[i]/(nev[i]-1));
      mean_err[i] = sig[i]/sqrt(nev[i]);
      sk[i] = (sumcu[i]/nev[i] - 3.*mean[i]*sig[i]*sig[i] - mean[i]*mean[i]*mean[i])/(sig[i]*sig[i]*sig[i]);
      sk_err[i] = sqrt(6./nev[i]);
      cout << i+1  << ": " << nev[i] << " " << mean[i] << " +- " << mean_err[i] << ", " << sig[i] << " " << sk[i] << " +- " << sk_err[i] << endl;
    } else{
      mean[i] = 0; mean_err[i] = 0;
      sig[i] = 0;
      sk[i] = 0; sk_err[i] = 0;
    }
  }
  return;
}

void draw_slices(const string& label,const string& parname, const int svd, const int mode, const bool single, const int N,const double* par, const double* par_err, const double* mean, const double* mean_err, const double* sig, const double* sk, const double* sk_err){
  cout << label << ", " << parname << endl;
  stringstream out;
  out.str("");
  out << "c_" << parname << "_svd" << svd;
  if(mode){
    out << "_mode" << mode;
  } else{
    out << "_tag";
    if(single) out << "_sgl";
  }
  string parn;
  if(parname == string("s") && mode)    parn = string("#sigma^{sig} (ps)");
  if(parname == string("s") && !mode)   parn = string("#sigma^{asc} (ps)");
  if(parname == string("h") && mode)    parn = string("h^{sig}");
  if(parname == string("h") && !mode)   parn = string("h^{asc}");
  if(parname == string("hst") && mode)  parn = string("h#times#sigma^{sig}");
  if(parname == string("hst") && !mode) parn = string("h#times#sigma^{asc}");
//  TCanvas* c1 = new TCanvas(out.str().c_str(),out.str().c_str(),1200,400);
  TCanvas* c1 = new TCanvas(out.str().c_str(),out.str().c_str(),800,400);
  TGraphErrors* gr_mean = new TGraphErrors(N,par,mean,par_err,mean_err);
  string title = string("Mean, ") + label;
  gr_mean->SetTitle(title.c_str());
  gr_mean->SetMarkerStyle(21);
  gr_mean->SetMarkerColor(kBlue);

  TGraph* gr_sigma = new TGraph(N,par,sig);
  title = string("RMS, ") + label;
  gr_sigma->SetTitle(title.c_str());
  gr_sigma->SetMarkerStyle(21);
  gr_sigma->SetMarkerColor(kBlue);

//  TGraphErrors* gr_skew = new TGraphErrors(N,par,sk,par_err,sk_err);
//  gr_skew->SetTitle("Skewness");
//  gr_skew->SetMarkerStyle(21);
//  gr_skew->SetMarkerColor(kBlue);

//  c1->Divide(3,1);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetGrid();
  gr_mean->GetXaxis()->SetLabelSize(0.05);
  gr_mean->GetXaxis()->SetTitle(parn.c_str());
  gr_mean->GetXaxis()->SetTitleSize(0.06);
  gr_mean->GetXaxis()->SetTitleOffset(0.75);
  gr_mean->GetYaxis()->SetLabelSize(0.05);
  gr_mean->GetYaxis()->SetTitle("mean (ps)");
  gr_mean->GetYaxis()->SetTitleSize(0.05);
//  gr_mean->GetYaxis()->SetTitleOffset(0.85);
  gr_mean->GetYaxis()->SetRangeUser(-1.,1.);
  gr_mean->Draw("AP");

  c1->cd(2);
  gPad->SetGrid();
  gr_sigma->GetXaxis()->SetLabelSize(0.05);
  gr_sigma->GetXaxis()->SetTitle(parn.c_str());
  gr_sigma->GetXaxis()->SetTitleSize(0.06);
  gr_sigma->GetXaxis()->SetTitleOffset(0.75);
  gr_sigma->GetYaxis()->SetLabelSize(0.05);
  gr_sigma->GetYaxis()->SetTitle("RMS (ps)");
  gr_sigma->GetYaxis()->SetTitleSize(0.05);
//  gr_sigma->GetYaxis()->SetTitleOffset(0.85);
//  gr_sigma->GetYaxis()->SetRangeUser(-1.,1.);
  gr_sigma->Draw("AP");

//  c1->cd(3);
//  gPad->SetGrid();
//  gr_skew->GetXaxis()->SetLabelSize(0.06);
//  gr_skew->GetYaxis()->SetLabelSize(0.06);
//  gr_skew->Draw("AP");

  c1->Update();
  string str = out.str();
  out.str("");
  out << str << ".root";
  c1->Print(out.str().c_str());
//  out.str("");
//  out << str << ".eps";
//  c1->Print(out.str().c_str());
//  out.str("");
//  out << str << ".png";
//  c1->Print(out.str().c_str());
//  out.str("");
//  out << "display " << out.str() << " &";
//  system(out.str().c_str());
  return;
}

int sz_chi2_slices(const int mode,const int svd,const bool pipi=false){
  TFile *file;
  string label;
  switch(mode){
  case 1:
    file = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s7_full.root");
    label = string("#pi^{0}");
    break;
  case 2:
    file = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    label = string("#eta#rightarrow#gamma#gamma");
    break;
  case 3:
    file = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_s3_full.root");
    label = string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
    break;
  case 4:
    file = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s6_full.root");
    label = string("#omega");
    break;
  case 5:
    file = TFile::Open("/home/vitaly/B0toDh0/Bp2D0pi/FIL_bp2d0pip_sigmc.root");
    label = string("D^{-}#pi^{+}");
    break;
  case 6:
    file = TFile::Open("/home/vitaly/B0toDh0/Bp2D0pi/FIL_bp2d0pip_sigmc.root");
    label = string("D^{-}#pi^{+} (Only D^{0}})");
    break;
  case 7:
    file = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigRho_s1_full.root");
    label = string("#rho");
    break;
  }

  if(svd == 2)      label += string(", SVD2");
  else if(svd == 1) label += string(", SVD1");
  else return-1;

  TTree *tree = (TTree*)file->Get("TEvent");

//  const double bdtg_cut = 0.125;
  const double dt_max = 10;
  const double dt_min = -dt_max;

  const double mm2ps = 7.84857;
  const int NBins = 25;
  const double st_min = 0;
  const double st_sig_max     = mode > 2 ? 0.2*mm2ps : 0.5*mm2ps;
//  const double st_asc_sgl_max = 0.5*mm2ps;
  const double st_asc_sgl_max = 0.2*mm2ps;
  const double st_asc_mlt_max = 0.2*mm2ps;
  const double h_min = 0;
//  const double h_sig_max = 50;
//  const double h_asc_max = 50;
  const double h_sig_max = 10;
  const double h_asc_max = 10;
  const double dst_sig = (st_sig_max - st_min)/NBins;
  const double dst_mlt_asc = (st_asc_mlt_max - st_min)/NBins;
  const double dst_sgl_asc = (st_asc_sgl_max - st_min)/NBins;
  const double dh_sig = (h_sig_max - h_min)/NBins;
  const double dh_asc = (h_asc_max - h_min)/NBins;

  const double hst_min     = st_min*h_min;
  const double hst_sig_max = 10;//st_sig_max*h_sig_max;
  const double hst_asc_max = 5;//st_asc_mlt_max*h_asc_max;
  const double dhst_sig    = (hst_sig_max-hst_min)/NBins;
  const double dhst_asc    = (hst_asc_max-hst_min)/NBins;

  double dt_sig;
  double t_sig, t_sig_mc;
  double dt_asc;
  double t_asc, t_asc_mc;
  double st_sig;
  double st_asc;
  double hst_sig;
  double hst_asc;
  double h_sig;
  double h_asc;
  int ndf_sig;
  int ndf_asc;
  int nptag;
  int exp;
  int b0f;
  int good_icpv;
//  double dt;
  if(mode<5 || mode>6){
    tree->SetBranchAddress("b0f",&b0f);
    tree->SetBranchAddress("good_icpv",&good_icpv);
  } else{
    tree->SetBranchAddress("bpf",&b0f);
    if(mode == 5) tree->SetBranchAddress("good_icpv_mlt",&good_icpv);
    else          tree->SetBranchAddress("good_icpv_sgl",&good_icpv);
  }

  if(pipi){
    tree->SetBranchAddress("z_sig_pipi",&t_sig);
//    tree->SetBranchAddress("dz_pipi",&dt);
    tree->SetBranchAddress("sz_sig_pipi",&st_sig);
    tree->SetBranchAddress("chisq_sig_pipi",&h_sig);
  } else if(mode != 6){
    tree->SetBranchAddress("z_sig",&t_sig);
//    tree->SetBranchAddress("dz",&dt);
    tree->SetBranchAddress("sz_sig",&st_sig);
    tree->SetBranchAddress("chisq_z_sig",&h_sig);
  } else{//D0 pi+, only D0
    tree->SetBranchAddress("z_sig_d0",&t_sig);
//    tree->SetBranchAddress("dz_d0",&dt);
    tree->SetBranchAddress("sz_sig_d0",&st_sig);
    tree->SetBranchAddress("chisq_z_sig_d0",&h_sig);
  }

  tree->SetBranchAddress("z_sig_mc",&t_sig_mc);
  tree->SetBranchAddress("z_asc",&t_asc);
  tree->SetBranchAddress("z_asc_mc",&t_asc_mc);
  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("sz_asc",&st_asc);
  tree->SetBranchAddress("chisq_z_asc",&h_asc);
  tree->SetBranchAddress("ndf_z_sig",&ndf_sig);
  tree->SetBranchAddress("ndf_z_asc",&ndf_asc);
  tree->SetBranchAddress("nptag",&nptag);

//  double chisq_vec[12] = {0.25,0.75,1.25,1.75,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
//  double sz_vec[12] = {0.01,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.08,0.14,0.18};
  double h_err[NBins], s_err[NBins];
  double h_sig_vec[NBins];
  double h_asc_vec[NBins];
  double st_sig_vec[NBins];
  double st_asc_mlt_vec[NBins];
  double st_asc_sgl_vec[NBins];

  double hst_sig_vec[NBins];
  double hst_asc_vec[NBins];
  double sum_hst[NBins];
  double sumsq_hst[NBins];
  double sumcu_hst[NBins];
  int nev_hst[NBins];
  double sum_hst_asc[NBins];
  double sumsq_hst_asc[NBins];
  double sumcu_hst_asc[NBins];
  int nev_hst_asc[NBins];

  double sum_h[NBins];
  double sumsq_h[NBins];
  double sumcu_h[NBins];
  int nev_h[NBins];

  double sum_st[NBins];
  double sumsq_st[NBins];
  double sumcu_st[NBins];
  int nev_st[NBins];

  double sum_h_asc[NBins];
  double sumsq_h_asc[NBins];
  double sumcu_h_asc[NBins];
  int nev_h_asc[NBins];

  double sum_st_asc_mlt[NBins];
  double sumsq_st_asc_mlt[NBins];
  double sumcu_st_asc_mlt[NBins];
  int nev_st_asc_mlt[NBins];

  double sum_st_asc_sgl[NBins];
  double sumsq_st_asc_sgl[NBins];
  double sumcu_st_asc_sgl[NBins];
  int nev_st_asc_sgl[NBins];

  for(int i=0; i<NBins; i++){
    sum_h[i] = 0;
    sumsq_h[i] = 0;
    sumcu_h[i] = 0;
    nev_h[i] = 0;
    h_err[i] = 0;
    s_err[i] = 0;

    sum_st[i] = 0;
    sumsq_st[i] = 0;
    sumcu_st[i] = 0;
    nev_st[i] = 0;

    sum_h_asc[i] = 0;
    sumsq_h_asc[i] = 0;
    sumcu_h_asc[i] = 0;
    nev_h_asc[i] = 0;

    sum_st_asc_mlt[i] = 0;
    sumsq_st_asc_mlt[i] = 0;
    sumcu_st_asc_mlt[i] = 0;
    nev_st_asc_mlt[i] = 0;

    sum_st_asc_sgl[i] = 0;
    sumsq_st_asc_sgl[i] = 0;
    sumcu_st_asc_sgl[i] = 0;
    nev_st_asc_sgl[i] = 0;

    sum_hst[i] = 0;
    sumsq_hst[i] = 0;
    sumcu_hst[i] = 0;
    nev_hst[i] = 0;
    sum_hst_asc[i] = 0;
    sumsq_hst_asc[i] = 0;
    sumcu_hst_asc[i] = 0;
    nev_hst_asc[i] = 0;

    h_sig_vec[i] = h_min + (i+0.5)*dh_sig;
    h_asc_vec[i] = h_min + (i+0.5)*dh_asc;
    st_sig_vec[i] = st_min + (i+0.5)*dst_sig;
    st_asc_mlt_vec[i] = st_min + (i+0.5)*dst_mlt_asc;
    st_asc_sgl_vec[i] = st_min + (i+0.5)*dst_sgl_asc;
    hst_sig_vec[i] = hst_min + (i+0.5)*dhst_sig;
    hst_asc_vec[i] = hst_min + (i+0.5)*dhst_asc;
  }

  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    dt_sig = t_sig - t_sig_mc;
    dt_asc = t_asc - t_asc_mc;
    if(!good_icpv) continue;
    if(b0f<1)      continue;
    if(svd == 1 && exp>30 || svd == 2 && exp<30) continue;
    dt_sig *= mm2ps; dt_asc *= mm2ps;
    st_sig *= mm2ps; st_asc *= mm2ps;
    if(pipi) ndf_sig -= 2;
    if(abs(t_sig - t_asc)*mm2ps>70) continue;

    if(ndf_sig) h_sig /= ndf_sig;
    if(ndf_asc) h_asc /= ndf_asc;
    hst_sig = h_sig*st_sig;
    hst_asc = h_asc*st_asc;

    if(st_sig>0 && st_sig<st_sig_max){
      int bin_st = (int)((st_sig - st_min)/dst_sig);
      sum_st[bin_st]   += dt_sig;
      sumsq_st[bin_st] += dt_sig*dt_sig;
      sumcu_st[bin_st] += dt_sig*dt_sig*dt_sig;
      nev_st[bin_st]++;
    }

    if(h_sig>0 && h_sig<h_sig_max){
      int bin_h       = (int)((h_sig - h_min)/dh_sig);
      sum_h[bin_h]   += dt_sig;
      sumsq_h[bin_h] += dt_sig*dt_sig;
      sumcu_h[bin_h] += dt_sig*dt_sig*dt_sig;
      nev_h[bin_h]++;
    }

    if(hst_sig>0 && hst_sig<hst_sig_max){
      int bin_hst         = (int)((hst_sig - hst_min)/dhst_sig);
      sum_hst[bin_hst]   += dt_sig;
      sumsq_hst[bin_hst] += dt_sig*dt_sig;
      sumcu_hst[bin_hst] += dt_sig*dt_sig*dt_sig;
      nev_hst[bin_hst]++;
    }

    if(!nptag && ndf_asc && st_asc>0 && st_asc<st_asc_mlt_max){
      int bin_st                = (int)((st_asc - st_min)/dst_mlt_asc);
      sum_st_asc_mlt[bin_st]   += dt_asc;
      sumsq_st_asc_mlt[bin_st] += dt_asc*dt_asc;
      sumcu_st_asc_mlt[bin_st] += dt_asc*dt_asc*dt_asc;
      nev_st_asc_mlt[bin_st]++;
    }

    if(!nptag && !ndf_asc && st_asc>0 && st_asc<st_asc_sgl_max){
      int bin_st                = (int)((st_asc - st_min)/dst_sgl_asc);
      sum_st_asc_sgl[bin_st]   += dt_asc;
      sumsq_st_asc_sgl[bin_st] += dt_asc*dt_asc;
      sumcu_st_asc_sgl[bin_st] += dt_asc*dt_asc*dt_asc;
      nev_st_asc_sgl[bin_st]++;
    }

    if(!nptag && ndf_asc && h_asc>h_min && h_asc<h_asc_max){
      int bin_h = (int)((h_asc - h_min)/dh_asc);
      sum_h_asc[bin_h]   += dt_asc;
      sumsq_h_asc[bin_h] += dt_asc*dt_asc;
      sumcu_h_asc[bin_h] += dt_asc*dt_asc*dt_asc;
      nev_h_asc[bin_h]++;
    }

    if(!nptag && ndf_asc && hst_asc>hst_min && hst_asc<hst_asc_max){
      int bin_hst             = (int)((hst_asc - hst_min)/dhst_asc);
      sum_hst_asc[bin_hst]   += dt_asc;
      sumsq_hst_asc[bin_hst] += dt_asc*dt_asc;
      sumcu_hst_asc[bin_hst] += dt_asc*dt_asc*dt_asc;
      nev_hst_asc[bin_hst]++;
    }
  }

  double mean_vec[NBins];
  double sigma_vec[NBins];
  double mean_err_vec[NBins];
  double skew_vec[NBins];
  double skew_err_vec[NBins];

  GetMoments(NBins,nev_st,sum_st,sumsq_st,sumcu_st,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  draw_slices(label,string("s"),svd,mode,true,NBins,st_sig_vec,s_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  if(mode>2 && mode != 6){
    GetMoments(NBins,nev_h,sum_h,sumsq_h,sumcu_h,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
    draw_slices(label,string("h"),svd,mode,true,NBins,h_sig_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);

    GetMoments(NBins,nev_hst,sum_hst,sumsq_hst,sumcu_hst,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
    draw_slices(label,string("hst"),svd,mode,true,NBins,hst_sig_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  }

  GetMoments(NBins,nev_st_asc_mlt,sum_st_asc_mlt,sumsq_st_asc_mlt,sumcu_st_asc_mlt,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  draw_slices(label,string("s"),svd,0,false,NBins,st_asc_mlt_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);

  GetMoments(NBins,nev_st_asc_sgl,sum_st_asc_sgl,sumsq_st_asc_sgl,sumcu_st_asc_sgl,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  draw_slices(label,string("s"),svd,0,true,NBins,st_asc_sgl_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);

  GetMoments(NBins,nev_h_asc,sum_h_asc,sumsq_h_asc,sumcu_h_asc,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  draw_slices(label,string("h"),svd,0,false,NBins,h_asc_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);

  GetMoments(NBins,nev_hst_asc,sum_hst_asc,sumsq_hst_asc,sumcu_hst_asc,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);
  draw_slices(label,string("hst"),svd,0,false,NBins,hst_asc_vec,h_err,mean_vec,mean_err_vec,sigma_vec,skew_vec,skew_err_vec);

  return;
}
