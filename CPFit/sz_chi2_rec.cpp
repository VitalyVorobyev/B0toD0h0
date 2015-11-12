int sz_chi2_rec(const int svd = 2, const int fit_mode = 1){
  TFile* file = new TFile("../TMVA/FIL_b2dh_sigOmega_s2_full.root");
  TTree* tree = (TTree*)file->Get("TEvent");

  const double bdtg_cut = 0.125;
  const double dzmax = 10/7.848;
  const double dzmin = -dzmax;

  double dz_sig;
  double dz_asc;
  double sz_sig;
  double sz_asc;
  double chi2_sig;
  double chi2_asc;
  int exp;
  tree->SetBranchAddress("dz_mc_sig",&dz_sig);
  tree->SetBranchAddress("dz_mc_asc",&dz_asc);
  tree->SetBranchAddress("exp",&exp);

  tree->SetBranchAddress("sz_sig",&sz_rec);
  tree->SetBranchAddress("sz_asc",&sz_asc);
  tree->SetBranchAddress("chisq_z_sig",&chi2_rec);
  tree->SetBranchAddress("chisq_z_asc",&chi2_asc);
  tree->SetBranchAddress("ndf_z_asc",&ndf_z_asc);

  const int NBins = 12;
//  double chisq_vec[12] = {0.25,0.75,1.25,1.75,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
//  double sz_vec[12] = {0.01,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.08,0.14,0.18};
  double chisq_err[NBins];

  double sum_sig[NBins];
  double sumsq_sig[NBins];
  double sumcu_sig[NBins];
  int nev[NBins];

  double sum_sz_sig[NBins];
  double sumsq_sz_sig[NBins];
  double sumcu_sz_sig[NBins];
  int nev_sz_sigNBins;

  for(int i=0; i<NBins; i++){
    sum_sig[i] = 0;
    sumsq_sig[i] = 0;
    sumcu_sig[i] = 0;
    nev_sig[i] = 0;
    chisq_err[i] = 0;
    sum_sz_sig[i] = 0;
    sumsq_sz_sig[i] = 0;
    sumcu_sz_sig[i] = 0;
    nev_sz_sig[i] = 0;
  }

  const int NTot = tree->GetEntries();
  double dz;
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    if(svd == 1 && exp>30 || svd == 2 && exp<30) continue;
//    if(de<-0.06 || de>0.05) continue;
//    if(mbc<5.272 || mbc>5.286) continue;
    if(vt_err_rec>0.2 || vt_err_asc>0.2 || vt_chi2_rec>10 || vt_chi2_asc>10) continue;
//    if(bdtg<bdtg_cut) continue;
    if(abs(dz_sig)>dtmax) continue;

    for(int j=0; j<NBins; j++){
      if(j==0){
        if(vt_chi2_rec<0.5){
          sum[j] += dz;
          sumsq[j] += dz*dz;
          sumcu[j] += dz*dz*dz;
          nev[j]++;
//                cout << "chi2 bin " << j;
         }
                if(err<0.015){
                    sum_sz[j] += dz;
                    sumsq_sz[j] += dz*dz;
                    sumcu_sz[j] += dz*dz*dz;
                    nev_sz[j]++;
  //                cout << ", sz bin " << j << endl;
                }
            } else if(j==11){
                if(vt_chi2_rec>9){
                    sum[j] += dz;
                    sumsq[j] += dz*dz;
                    sumcu[j] += dz*dz*dz;
                    nev[j]++;
//                    cout << "chi2 bin " << j << endl;
                }
                if(err>0.16){
                    sum_sz[j] += dz;
                    sumsq_sz[j] += dz*dz;
                    sumcu_sz[j] += dz*dz*dz;
                    nev_sz[j]++;
//                    cout << ", sz bin " << j << endl;
                }
            } else{
                if(vt_chi2_rec<0.5*(chisq_vec[j]+chisq_vec[j+1]) && vt_chi2_rec>0.5*(chisq_vec[j-1]+chisq_vec[j])){
                    sum[j] += dz;
                    sumsq[j] += dz*dz;
                    sumcu[j] += dz*dz*dz;
                    nev[j]++;
//                    cout << "chi2 bin " << j << endl;
                }
                if(err<0.5*(sz_vec[j]+sz_vec[j+1]) && err>0.5*(sz_vec[j-1]+sz_vec[j])){
                    sum_sz[j] += dz;
                    sumsq_sz[j] += dz*dz;
                    sumcu_sz[j] += dz*dz*dz;
                    nev_sz[j]++;
//                    cout << ", sz bin " << j << endl;
                }
            }
//            if(j != 0 && j != 11) cout << 0.5*(sz_vec[j]+sz_vec[j+1]) << " - " << 0.5*(sz_vec[j-1]+sz_vec[j]) << endl;
        }
//        cout << "sz = " << vt_err_rec << endl;
    }

    double mean_vec[12];
    double sigma_vec[12];
    double mean_err_vec[12];
    double skew_vec[12];
    double skew_err_vec[12];
    cout << "Chi2 summary:" << endl;
    for(int i=0; i<12; i++){
        mean_vec[i] = sum[i]/nev[i];
        sigma_vec[i] = sqrt(sumsq[i]/nev[i] - mean_vec[i]*mean_vec[i]);
        mean_err_vec[i] = sigma_vec[i]/sqrt(nev[i]);
        skew_vec[i] = (sumcu[i]/nev[i] - 3.*mean_vec[i]*sigma_vec[i]*sigma_vec[i] - mean_vec[i]*mean_vec[i]*mean_vec[i])/(sigma_vec[i]*sigma_vec[i]*sigma_vec[i]);
        skew_err_vec[i] = sqrt(6./nev[i]);
        cout << " " << nev[i] << " " << chisq_vec[i] << " " << mean_vec[i] << " " << sigma_vec[i] << " " << skew_vec[i] << endl;
    }
    cout << endl;

    TCanvas* c = new TCanvas("chi2","chi2",1200,400);
    TGraphErrors* gr_mean = new TGraphErrors(12,chisq_vec,mean_vec,chisq_err,mean_err_vec);
    gr_mean->SetMarkerStyle(21);
    gr_mean->SetMarkerColor(kBlue);

    TGraph* gr_sigma = new TGraph(12,chisq_vec,sigma_vec);
    gr_sigma->SetMarkerStyle(21);
    gr_sigma->SetMarkerColor(kBlue);

    TGraphErrors* gr_skew = new TGraphErrors(12,chisq_vec,skew_vec,chisq_err,skew_err_vec);
    gr_skew->SetMarkerStyle(21);
    gr_skew->SetMarkerColor(kBlue);

    c->Divide(3,1);
    c->cd(1);
    gr_mean->Draw("AP");

    c->cd(2);
    gr_sigma->Draw("AP");

    c->cd(3);
    gr_skew->Draw("AP");

    c->Update();

    double mean_vec_sz[12];
    double sigma_vec_sz[12];
    double mean_err_vec_sz[12];
    double skew_vec_sz[12];
    double skew_err_vec_sz[12];
    cout << "Sz summary:" << endl;
    for(int i=0; i<12; i++){
        mean_vec_sz[i] = sum_sz[i]/nev_sz[i];
        sigma_vec_sz[i] = sqrt(sumsq_sz[i]/nev_sz[i] - mean_vec_sz[i]*mean_vec_sz[i]);
        mean_err_vec_sz[i] = sigma_vec_sz[i]/sqrt(nev_sz[i]);
        skew_vec_sz[i] = (sumcu_sz[i]/nev_sz[i] - 3.*mean_vec_sz[i]*sigma_vec_sz[i]*sigma_vec_sz[i] - mean_vec_sz[i]*mean_vec_sz[i]*mean_vec_sz[i])/(sigma_vec_sz[i]*sigma_vec_sz[i]*sigma_vec_sz[i]);
        skew_err_vec_sz[i] = sqrt(6./nev_sz[i]);
        cout << " " << nev_sz[i] << " " << sz_vec[i] << " " << mean_vec_sz[i] << " " << sigma_vec_sz[i] << " " << skew_vec_sz[i] << endl;
    }

    TCanvas* c = new TCanvas("sz","sz",1200,400);
    TGraphErrors* gr_mean_sz = new TGraphErrors(12,sz_vec,mean_vec_sz,chisq_err,mean_err_vec_sz);
    gr_mean_sz->SetMarkerStyle(21);
    gr_mean_sz->SetMarkerColor(kRed);

    TGraph* gr_sigma_sz = new TGraph(12,sz_vec,sigma_vec_sz);
    gr_sigma_sz->SetMarkerStyle(21);
    gr_sigma_sz->SetMarkerColor(kRed);

    TGraphErrors* gr_skew_sz = new TGraphErrors(12,sz_vec,skew_vec_sz,chisq_err,skew_err_vec_sz);
    gr_skew_sz->SetMarkerStyle(21);
    gr_skew_sz->SetMarkerColor(kRed);

    c->Divide(3,1);
    c->cd(1);
    gr_mean_sz->Draw("AP");

    c->cd(2);
    gr_sigma_sz->Draw("AP");

    c->cd(3);
    gr_skew_sz->Draw("AP");

    c->Update();
}
