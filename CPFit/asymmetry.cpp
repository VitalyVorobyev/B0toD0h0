void asymmetry(void){
  TChain* tree = new TChain("TEvent","TEvent");
  tree->Add("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m1_mh010.root");

  RooArgSet argset;
  RooRealVar z_sig("z_sig","z_sig",-15,15,"mm"); argset.add(z_sig);
  RooRealVar z_asc("z_asc","z_asc",-15,15,"mm"); argset.add(z_asc);
  RooConstant mm2ps("mm2ps","mm2ps",7.848566945838871754705);
  RooFormulaVar dt("dt","(@0-@1)*@2",RooArgList(z_sig,z_asc,mm2ps);

  RooRealVar sz_sig("sz_sig","sz_sig",0.,0.5,"mm"); argset.add(sz_sig);
  RooRealVar sz_asc("sz_asc","sz_asc",0.,0.5,"mm"); argset.add(sz_asc);

  RooRealVar chisq_sig("chisq_sig","chisq_sig",0.,1.e6); argset.add(chisq_sig);
  RooRealVar chisq_asc("chisq_sig","chisq_sig",0.,1.e6); argset.add(chisq_asc);

//  RooCategory ntrk_sig("ntrk_sig","ntrk_sig");

  const int ntrk_vec_size = 12;
  int ntrk_asc_vec[ntrk_vec_size] = {1,2,3,4,5,6,7,8,9,10,11,12};
  RooCategory ndf_asc("ndf_asc","ndf_asc");
  RooCategory ntrk_asc("ntrk_asc","ntrk_asc");
  stringstream out;
  for(int i=0; i<ntrk_vec_size; i++){
    out.str("");
    out << ntrk_asc_vec[i];
    ntrk_asc->defineType(out.str().c_str(),ntrk_asc_vec[i]);

    out.str("");
    out << 2*ntrk_asc_vec[i]-1;
    ndf_asc->defineType(out.str().c_str(),2*ntrk_asc_vec[i]-2);
  }

  
  
}
