void make_text(char* ifile,char* ofile){
  TFile file(ifile);
  TTree* tree = (TTree*)file.Get("TEvent");
  const int NTot = tree->GetEntries();
  cout << NTot << " event detected" << endl;
  
  Int_t exp, run, evt;
  Int_t b0f, d0f;
  Double_t bdtg;
  Int_t bin, mc_bin;
  Int_t mc_flv;
  Double_t tag_LH;
  Double_t vt_mc_rec, vt_mc_asc, mc_delta_t;
  Double_t vt_pos_rec, vt_err_rec;
  Int_t vt_ntrk_rec;
  Double_t vt_chi2_rec;
  Int_t vt_ndf_rec;
  Double_t vt_pos_asc, vt_err_asc;
  Int_t vt_ntrk_asc;
  Double_t vt_chi2_asc;
  Int_t vt_ndf_asc;
  Double_t deltae, mbc, costh;
  Int_t nptag;

  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evt);
  tree->SetBranchAddress("b0f",&b0f);
  tree->SetBranchAddress("d0f",&d0f);
  tree->SetBranchAddress("bdtgsl",&bdtg);
  tree->SetBranchAddress("bin",&bin);
//  tree->SetBranchAddress("bin_mc",&mc_bin);
//  tree->SetBranchAddress("mc_flv",&mc_flv);
  tree->SetBranchAddress("tag_LH",&tag_LH);
//  tree->SetBranchAddress("z_sig_mc",&vt_mc_rec);
//  tree->SetBranchAddress("z_asc_mc",&vt_mc_asc);
  tree->SetBranchAddress("z_sig",&vt_pos_rec);
  tree->SetBranchAddress("sz_sig",&vt_err_rec);
  tree->SetBranchAddress("ntrk_sig",&vt_ntrk_rec);
  tree->SetBranchAddress("chisq_z_sig",&vt_chi2_rec);
  tree->SetBranchAddress("ndf_z_sig",&vt_ndf_rec);
  tree->SetBranchAddress("z_asc",&vt_pos_asc);
  tree->SetBranchAddress("sz_asc",&vt_err_asc);
  tree->SetBranchAddress("ntrk_asc",&vt_ntrk_asc);
  tree->SetBranchAddress("chisq_z_asc",&vt_chi2_asc);
  tree->SetBranchAddress("ndf_z_asc",&vt_ndf_asc);
  tree->SetBranchAddress("de",&deltae);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("costhBcms",&costh);
  tree->SetBranchAddress("nptag",&nptag);

  ofstream output(ofile);
  for(int i=0; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);
    output << exp << " ";
    output << run << " ";
    output << evt << " ";
    output << b0f << " ";
    output << d0f << " ";
    output << bdtg << " ";
    output << bin << " ";
//    output << mc_bin << " ";
//    output << mc_flv << " ";
    output << tag_LH << " ";
//    output << vt_mc_rec << " ";
//    output << vt_mc_asc << " ";
    output << vt_pos_rec << " ";
    output << vt_err_rec << " ";
    output << vt_ntrk_rec << " ";
    output << vt_chi2_rec << " ";
    output << vt_ndf_rec << " ";
    output << vt_pos_asc << " ";
    output << vt_err_asc << " ";
    output << vt_ntrk_asc << " ";
    output << vt_chi2_asc << " ";
    output << vt_ndf_asc << " ";
    output << deltae << " ";
    output << mbc << " ";
    output << costh << " ";
    output << nptag;
    output << endl;
  }
  return;  
}

