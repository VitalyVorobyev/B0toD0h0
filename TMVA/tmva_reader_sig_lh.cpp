void tmva_reader_sig_lh(){
/////////////////
// TMVA Reader //
/////////////////
  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent:V");

  TFile *ifile = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/Tuples/b2dh_sigmc_s1.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");
//  ifile->Close();

  TFile ofile("fil_b2dh_sig_lh.root","RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");
  
  Double_t cos_b0,p_ks,p_pp,p_pm,p_pi0,chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx,cos_thr,thr_sig,thr_oth;
  Double_t lh;

  Double_t mbc,de,bdt,bdtg,mp,mm,dz,dt,atckpi_max,mpi0,mk,md;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err,dzerr;
  Int_t phsp,bin,exp,run,evtn;

  Float_t m_cos_b0,m_p_ks,m_p_pp,m_p_pm,m_p_pi0,m_chi2_ndf_D0,m_chi2_ndf_B0,m_chi2_tag_vtx,m_cos_thr,m_thr_sig,m_thr_oth;
  Float_t m_ks_dr,m_ks_dz,m_ks_dphi,m_ks_fl,m_tag_LH,m_tag_LH_err,m_dzerr;
  Float_t m_lh;

  Double_t mp_mc,mm_mc;
  Double_t t_mc,d0_t_mc;
  Double_t dt_mc,dz_mc;
  Int_t bin_mc;
  Int_t flv_mc,d0_flv_mc;
  Int_t mcflag,d0_mcflag;

  tree->SetBranchAddress("mp_mc",&mp_mc);
  tree->SetBranchAddress("mm_mc",&mm_mc);
  tree->SetBranchAddress("t_mc",&t_mc);
  tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
  tree->SetBranchAddress("dt_mc",&dt_mc);
  tree->SetBranchAddress("dz_mc",&dz_mc);
  tree->SetBranchAddress("bin_mc",&bin_mc);
  tree->SetBranchAddress("flv_mc",&flv_mc);
  tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
  tree->SetBranchAddress("mcflag",&mcflag);
  tree->SetBranchAddress("d0_mcflag",&d0_mcflag);

  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("bin",&bin);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("dt",&dt);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("phsp",&phsp);
  tree->SetBranchAddress("mpi0_raw",&mpi0);
  tree->SetBranchAddress("mks_raw",&mk);
  tree->SetBranchAddress("md0_raw",&md);

  tree->SetBranchAddress("cos_b0",&cos_b0);
  tree->SetBranchAddress("p_ks",&p_ks);
  tree->SetBranchAddress("p_pp",&p_pp);
  tree->SetBranchAddress("p_pm",&p_pm);
  tree->SetBranchAddress("p_pi0",&p_pi0);
  tree->SetBranchAddress("chi2_ndf_D0",&chi2_ndf_D0);
  tree->SetBranchAddress("chi2_ndf_B0",&chi2_ndf_B0);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);
  tree->SetBranchAddress("dzerr",&dzerr);
  
  tree->SetBranchAddress("ks_dr",&ks_dr);
  tree->SetBranchAddress("ks_dz",&ks_dz);
  tree->SetBranchAddress("ks_dphi",&ks_dphi);
  tree->SetBranchAddress("ks_fl",&ks_fl);

  tree->SetBranchAddress("lh",&lh);

  TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
  TEvent->Branch("mm_mc",&mm_mc,"mm_mc/D");
  TEvent->Branch("t_mc",&t_mc,"t_mc/D");
  TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
  TEvent->Branch("dt_mc",&dt_mc,"dt_mc/D");
  TEvent->Branch("dz_mc",&dz_mc,"dz_mc/D");
  TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
  TEvent->Branch("flv_mc",&flv_mc,"flv_mc/I");
  TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
  TEvent->Branch("mcflag",&mcflag,"mcflag/I");
  TEvent->Branch("d0_mcflag",&d0_mcflag,"d0_mcflag/I");
  
  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");
  
  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");
  TEvent->Branch("bdt",&bdt,"bdt/D");
  TEvent->Branch("bdtg",&bdtg,"bdtg/D");
  
  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");
  TEvent->Branch("bin",&bin,"bin/I");
  
  TEvent->Branch("dz",&dz,"dz/D");
  TEvent->Branch("dt",&dt,"dt/D");
  
  TEvent->Branch("tag_LH",&tag_LH,"tag_LH/D");
  TEvent->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");

  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("mpi0",&mpi0,"mpi0/D");
  TEvent->Branch("mk",&mk,"mk/D");
  TEvent->Branch("md",&md,"md/D");

  reader->AddVariable("abs(cos_b0)",&m_cos_b0);
//  reader->AddVariable("abs(cos_d0)");
//  reader->AddVariable("abs(pcm_d0-2.3)");
  reader->AddVariable("p_ks",&m_p_ks);
//  reader->AddVariable("p_pp",&m_p_pp);
//  reader->AddVariable("p_pm",&m_p_pm);
//  reader->AddVariable("abs(p_pi0-2.6)",&m_p_pi0);
//  reader->AddVariable("log(abs(r_pp))");
//  reader->AddVariable("log(abs(z_pp))");
//  reader->AddVariable("log(abs(r_pm))");
//  reader->AddVariable("log(abs(z_pm))");
  reader->AddVariable("log(chi2_ndf_D0)",&m_chi2_ndf_D0);
  reader->AddVariable("log(chi2_ndf_B0)",&m_chi2_ndf_B0);
  reader->AddVariable("log(chi2_tag_vtx)",&m_chi2_tag_vtx);
  reader->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader->AddVariable("thr_sig",&m_thr_sig);
  reader->AddVariable("thr_oth",&m_thr_oth);
  reader->AddVariable("log(tag_LH_err)",&m_tag_LH_err);
  reader->AddVariable("log(dzerr)",&m_dzerr);

//  reader->AddVariable("log(ks_dr)",&m_ks_dr);
//  reader->AddVariable("log(ks_dz)",&m_ks_dz);
//  reader->AddVariable("log(ks_dphi)",&m_ks_dphi);
//  reader->AddVariable("log(ks_fl)",&m_ks_fl);

  reader->AddVariable("lh",&m_lh);

  reader->BookMVA("BDT","weights/MVAnalysis_lh_BDT.weights.xml");
  reader->BookMVA("BDTG","weights/MVAnalysis_lh_BDTG.weights.xml");
  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);
//    if(!phsp) continue;
//    if(chi2_ndf_B0>1000) continue;

    m_cos_b0      = (float)cos_b0;
    m_p_ks        = (float)p_ks;
    m_p_pp        = (float)p_pp;
    m_p_pm        = (float)p_pm;
    m_p_pi0       = (float)p_pi0;
    m_chi2_ndf_D0 = (float)chi2_ndf_D0;
    m_chi2_ndf_B0 = (float)chi2_ndf_B0;
    m_cos_thr     = (float)cos_thr;
    m_thr_sig     = (float)thr_sig;
    m_thr_oth     = (float)thr_oth;
    m_dzerr       = (float)dzerr;
    m_ks_dr       = (float)ks_dr;
    m_ks_dz       = (float)ks_dz;
    m_ks_dphi     = (float)ks_dphi;
    m_ks_fl       = (float)ks_fl;
    m_lh          = (float)lh;

    bdt  = reader->EvaluateMVA("BDT");
    bdtg = reader->EvaluateMVA("BDTG");

    TEvent->Fill();
  }

  TEvent->Write();
  ofile.Write();
  ofile.Close();
}
