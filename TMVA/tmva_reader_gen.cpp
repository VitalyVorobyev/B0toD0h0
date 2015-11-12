void tmva_reader_gen(){
/////////////////
// TMVA Reader //
/////////////////
  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent:V");

  TFile *ifile = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/Tuples/b2dh_gen.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");
//  ifile->Close();

  TFile ofile("fil_b2dh_gen-2.root","RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  Double_t cos_b0,p_ks,p_pp,p_pm,p_pi0,chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx,cos_thr,thr_sig,thr_oth;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;

  Double_t mbc,de,bdt,bdtg,mp,mm,dz,dt,atckpi_max,mpi0,mk,md;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err,dzerr;
  Int_t phsp,bin,exp,run,evtn;

  Float_t m_cos_b0,m_p_ks,m_p_pp,m_p_pm,m_p_pi0,m_chi2_ndf_D0,m_chi2_ndf_B0,m_chi2_tag_vtx,m_cos_thr,m_thr_sig,m_thr_oth;
  Float_t m_ks_dr,m_ks_dz,m_ks_dphi,m_ks_fl,m_tag_LH,m_tag_LH_err,m_dzerr;
  Float_t m_k0mm2,m_k0et,m_k0hso00,m_k0hso02,m_k0hso04,m_k0hso10,m_k0hso12,m_k0hso14,m_k0hso20,m_k0hso22,m_k0hso24,m_k0hoo0,m_k0hoo1,m_k0hoo2,m_k0hoo3,m_k0hoo4;

  Int_t d0_chain[9];
  Int_t pi0_chain[9];
  Int_t pip_chain[9];
  Int_t pim_chain[9];
  Int_t k0s_chain[9];
  Int_t b0id,b0f;
  Int_t d0id,d0f;
  Int_t ksid,ksf;
  Int_t pi0id,pi0f;
  
  tree->SetBranchAddress("d0_chain",d0_chain);
  tree->SetBranchAddress("pi0_chain",pi0_chain);
  tree->SetBranchAddress("pip_chain",pip_chain);
  tree->SetBranchAddress("pim_chain",pim_chain);
  tree->SetBranchAddress("k0s_chain",k0s_chain);
  tree->SetBranchAddress("b0f",&b0f);
  tree->SetBranchAddress("b0id",&b0id);
  tree->SetBranchAddress("d0f",&d0f);
  tree->SetBranchAddress("d0id",&d0id);
  tree->SetBranchAddress("ksf",&ksf);
  tree->SetBranchAddress("ksid",&ksid);
  tree->SetBranchAddress("pi0f",&pi0f);
  tree->SetBranchAddress("pi0id",&pi0id);
  
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

  tree->SetBranchAddress("k0mm2",&k0mm2);
  tree->SetBranchAddress("k0et",&k0et);
  tree->SetBranchAddress("k0hso00",&k0hso00);
  tree->SetBranchAddress("k0hso02",&k0hso02);
  tree->SetBranchAddress("k0hso04",&k0hso04);
  tree->SetBranchAddress("k0hso10",&k0hso10);
  tree->SetBranchAddress("k0hso12",&k0hso12);
  tree->SetBranchAddress("k0hso14",&k0hso14);
  tree->SetBranchAddress("k0hso20",&k0hso20);
  tree->SetBranchAddress("k0hso22",&k0hso22);
  tree->SetBranchAddress("k0hso24",&k0hso24);
  tree->SetBranchAddress("k0hoo0",&k0hoo0);
  tree->SetBranchAddress("k0hoo1",&k0hoo1);
  tree->SetBranchAddress("k0hoo2",&k0hoo2);
  tree->SetBranchAddress("k0hoo3",&k0hoo3);
  tree->SetBranchAddress("k0hoo4",&k0hoo4);

  TEvent->Branch("pip_chain",pip_chain,"pipch0/I:pipch1/I:pipch2/I:pipch3/I:pipch4/I:pipch5/I:pipch6/I:pipch7/I:pipch8/I");
  TEvent->Branch("pim_chain",pim_chain,"pimch0/I:pimch1/I:pimch2/I:pimch3/I:pimch4/I:pimch5/I:pimch6/I:pimch7/I:pimch8/I");
  TEvent->Branch("k0s_chain",k0s_chain,"k0sch0/I:k0sch1/I:k0sch2/I:k0sch3/I:k0sch4/I:k0sch5/I:k0sch6/I:k0sch7/I:k0sch8/I");
  TEvent->Branch("pi0_chain",pi0_chain,"pi0ch0/I:pi0ch1/I:pi0ch2/I:pi0ch3/I:pi0ch4/I:pi0ch5/I:pi0ch6/I:pi0ch7/I:pi0ch8/I");
  TEvent->Branch("d0_chain",d0_chain,"d0ch0/I:d0ch1/I:d0ch2/I:d0ch3/I:d0ch4/I:d0ch5/I:d0ch6/I:d0ch7/I:d0ch8/I");

  TEvent->Branch("b0id",&b0id,"b0id/I");
  TEvent->Branch("b0f",&b0f,"b0f/I");
  TEvent->Branch("d0id",&d0id,"d0id/I");
  TEvent->Branch("d0f",&d0f,"d0f/I");
  TEvent->Branch("ksid",&ksid,"ksid/I");
  TEvent->Branch("ksf",&ksf,"ksf/I");
  TEvent->Branch("ksid",&ksid,"pi0id/I");
  TEvent->Branch("pi0f",&pi0f,"pi0f/I");
  
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

  reader->AddVariable("k0mm2",&m_k0mm2);
  reader->AddVariable("k0et",&m_k0et);
  reader->AddVariable("k0hso00",&m_k0hso00);
//  reader->AddVariable("k0hso01");
  reader->AddVariable("k0hso02",&m_k0hso02);
//  reader->AddVariable("k0hso03");
  reader->AddVariable("k0hso04",&m_k0hso04);
  reader->AddVariable("k0hso10",&m_k0hso10);
  reader->AddVariable("k0hso12",&m_k0hso12);
  reader->AddVariable("k0hso14",&m_k0hso14);
  reader->AddVariable("k0hso20",&m_k0hso20);
  reader->AddVariable("k0hso22",&m_k0hso22);
  reader->AddVariable("k0hso24",&m_k0hso24);
  reader->AddVariable("k0hoo0",&m_k0hoo0);
  reader->AddVariable("k0hoo1",&m_k0hoo1);
  reader->AddVariable("k0hoo2",&m_k0hoo2);
//  reader->AddVariable("k0hoo3",&m_k0hoo3);
  reader->AddVariable("k0hoo4",&m_k0hoo4);

  reader->BookMVA("BDT","weights/MVAnalysis_BDT.weights.xml");
  reader->BookMVA("BDTG","weights/MVAnalysis_BDTG.weights.xml");
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
    m_k0mm2       = (float)k0mm2;
    m_k0et        = (float)k0et;
    m_k0hso00     = (float)k0hso00;
    m_k0hso02     = (float)k0hso02;
    m_k0hso04     = (float)k0hso04;
    m_k0hso10     = (float)k0hso10;
    m_k0hso12     = (float)k0hso12;
    m_k0hso14     = (float)k0hso14;
    m_k0hso20     = (float)k0hso20;
    m_k0hso22     = (float)k0hso22;
    m_k0hso24     = (float)k0hso24;
    m_k0hoo0      = (float)k0hoo0;
    m_k0hoo1      = (float)k0hoo1;
    m_k0hoo2      = (float)k0hoo2;
    m_k0hoo3      = (float)k0hoo3;
    m_k0hoo4      = (float)k0hoo4;

    bdt  = reader->EvaluateMVA("BDT");
    bdtg = reader->EvaluateMVA("BDTG");

    TEvent->Fill();
  }
  
  TEvent->Write();
  ofile.Write();
  ofile.Close();
}
