#include "../rooksfw/ksfw.h"
#include "../rooksfw/rooksfw.h"
#include "../rooksfw/rooksfw.cc"

void tmva_reader(const int type = 1){
// type = 1 -> Signal
// type = 2 -> Continuum
// type = 0 -> Data
/////////////////
// TMVA Reader //
/////////////////
  rooksfw* ksfw  = new rooksfw("ksfw0", ksfw0_alpha, ksfw0_sigpdf, ksfw0_bkgpdf);
  TMVA::Reader* reader  = new TMVA::Reader("!Color:!Silent:V");
//  TMVA::Reader* reader1 = new TMVA::Reader("!Color:!Silent:V");

  string fname;
  if(type == 1)       fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmc.root");
  else if(type == 2 ) fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_cont.root");
  else if(type == 0)  fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_data.root");
  else{
      cout << "Wrong type " << type << endl;
      return;
  }
  stinrg treename;
  if(type) treename = string("TEventEx");
  else     treename = string("TEvent");

  TFile *ifile = TFile::Open(fname.c_str());
  TTree *tree = (TTree*)ifile->Get(treename.c_str());

  if(type == 1)       fname = string("fil_b2dh_sig.root");
  else if(type == 2) fname = string("fil_b2dh_cont.root");
  else if(type == 0)  fname = string("fil_b2dh_data.root");
  TFile ofile(fname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  double k0vars[17];

  Double_t cos_b0,p_ks,p_pp,p_pm,p_pi0,chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx,cos_thr,thr_sig,thr_oth;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;
  Double_t k1mm2,k1et,k1hso00,k1hso02,k1hso04,k1hso10,k1hso12,k1hso14,k1hso20,k1hso22,k1hso24,k1hoo0,k1hoo1,k1hoo2,k1hoo3,k1hoo4;
  Double_t pi0_chi2,egamma,ptgamma;
  Int_t ndf_tag_vtx;

  Double_t mbc,de,mp,mm,dz,dt,atckpi_max,mpi0,mk,md;
  Double_t bdt,bdtg,bdts,bdtgs,lh,bdtlh,bdtglh,bdtlhs,bdtglhs;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err,dzerr;
  Int_t phsp,bin,exp,run,evtn;

  Float_t m_cos_b0,m_p_ks,m_p_pp,m_p_pm,m_p_pi0,m_chi2_ndf_D0,m_chi2_ndf_B0,m_chi2_ndf_tag_vtx,m_cos_thr,m_thr_sig,m_thr_oth;
  Float_t m_ks_dr,m_ks_dz,m_ks_dphi,m_ks_fl,m_tag_LH_err,m_dzerr;
  Float_t m_pi0_chi2,m_egamma,m_ptgamma,m_lh;
  Float_t m_k1mm2,m_k1et,m_k1hso00,m_k1hso02,m_k1hso04,m_k1hso10,m_k1hso12,m_k1hso14,m_k1hso20,m_k1hso22,m_k1hso24,m_k1hoo0,m_k1hoo1,m_k1hoo2,m_k1hoo3,m_k1hoo4;

  Double_t mp_mc,mm_mc;
  Double_t t_mc,d0_t_mc;
  Double_t dt_mc,dz_mc;
  Int_t bin_mc;
  Int_t flv_mc,d0_flv_mc;
  Int_t mcflag,d0_mcflag;
  Int_t b0f,d0f;

  tree->SetBranchAddress("mp_mc",&mp_mc);
  tree->SetBranchAddress("mm_mc",&mm_mc);
  tree->SetBranchAddress("t_mc",&t_mc);
  tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
  tree->SetBranchAddress("dt_mc",&dt_mc);
  tree->SetBranchAddress("dz_mc",&dz_mc);
  tree->SetBranchAddress("bin_mc",&bin_mc);
  tree->SetBranchAddress("flv_mc",&flv_mc);
  tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
//  tree->SetBranchAddress("mcflag",&mcflag);
//  tree->SetBranchAddress("d0_mcflag",&d0_mcflag);
  tree->SetBranchAddress("b0f",&b0f);
  tree->SetBranchAddress("d0f",&d0f);

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
//  tree->SetBranchAddress("p_pp",&p_pp);
//  tree->SetBranchAddress("p_pm",&p_pm);
//  tree->SetBranchAddress("p_pi0",&p_pi0);
  tree->SetBranchAddress("chi2_ndf_D0",&chi2_ndf_D0);
  tree->SetBranchAddress("chi2_ndf_B0",&chi2_ndf_B0);
  tree->SetBranchAddress("chi2_tag_vtx",&chi2_tag_vtx);
  tree->SetBranchAddress("ndf_tag_vtx",&ndf_tag_vtx);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);
  tree->SetBranchAddress("dzerr",&dzerr);

  tree->SetBranchAddress("pi0_chi2",&pi0_chi2);
  tree->SetBranchAddress("egamma",&egamma);
//  tree->SetBranchAddress("ptgamma",&ptgamma);
//  tree->SetBranchAddress("lh",&m_lh);

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

  tree->SetBranchAddress("k1mm2",&k1mm2);
  tree->SetBranchAddress("k1et",&k1et);
  tree->SetBranchAddress("k1hso00",&k1hso00);
  tree->SetBranchAddress("k1hso02",&k1hso02);
  tree->SetBranchAddress("k1hso04",&k1hso04);
  tree->SetBranchAddress("k1hso10",&k1hso10);
  tree->SetBranchAddress("k1hso12",&k1hso12);
  tree->SetBranchAddress("k1hso14",&k1hso14);
  tree->SetBranchAddress("k1hso20",&k1hso20);
  tree->SetBranchAddress("k1hso22",&k1hso22);
  tree->SetBranchAddress("k1hso24",&k1hso24);
  tree->SetBranchAddress("k1hoo0",&k1hoo0);
  tree->SetBranchAddress("k1hoo1",&k1hoo1);
  tree->SetBranchAddress("k1hoo2",&k1hoo2);
  tree->SetBranchAddress("k1hoo3",&k1hoo3);
  tree->SetBranchAddress("k1hoo4",&k1hoo4);

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

  TEvent->Branch("b0f",&b0f,"b0f/I");
  TEvent->Branch("d0f",&d0f,"d0f/I");

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");
  TEvent->Branch("bdt",&bdt,"bdt/D");
  TEvent->Branch("bdtg",&bdtg,"bdtg/D");
  TEvent->Branch("bdts",&bdts,"bdts/D");
  TEvent->Branch("bdtgs",&bdtgs,"bdtgs/D");
  TEvent->Branch("bdtlh",&bdtlh,"bdtlh/D");
  TEvent->Branch("bdtglh",&bdtglh,"bdtglh/D");
  TEvent->Branch("bdtlhs",&bdtlhs,"bdtlhs/D");
  TEvent->Branch("bdtglhs",&bdtglhs,"bdtglhs/D");
  TEvent->Branch("lh",&lh,"lh/D");
  
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
//  reader->AddVariable("p_ks",&m_p_ks);
  reader->AddVariable("log(chi2_ndf_D0)",&m_chi2_ndf_D0);
//  reader->AddVariable("log(chi2_ndf_B0)",&m_chi2_ndf_B0);
//  reader->AddVariable("log(chi2_tag_vtx/ndf_tag_vtx)",&m_chi2_ndf_tag_vtx);
  reader->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader->AddVariable("abs(thr_sig-0.885)",&m_thr_sig);
  reader->AddVariable("thr_oth",&m_thr_oth);
  reader->AddVariable("log(tag_LH_err)",&m_tag_LH_err);
  reader->AddVariable("log(dzerr)",&m_dzerr);
  reader->AddVariable("log(pi0_chi2)",&m_pi0_chi2);
  reader->AddVariable("log(egamma)",&m_egamma);

  reader->AddVariable("k1mm2",&m_k1mm2);
  reader->AddVariable("k1et",&m_k1et);
  reader->AddVariable("k1hso00",&m_k1hso00);
  reader->AddVariable("k1hso02",&m_k1hso02);
  reader->AddVariable("k1hso04",&m_k1hso04);
  reader->AddVariable("k1hso10",&m_k1hso10);
  reader->AddVariable("k1hso12",&m_k1hso12);
  reader->AddVariable("k1hso14",&m_k1hso14);
  reader->AddVariable("k1hso20",&m_k1hso20);
  reader->AddVariable("k1hso22",&m_k1hso22);
  reader->AddVariable("k1hso24",&m_k1hso24);
  reader->AddVariable("k1hoo0",&m_k1hoo0);
  reader->AddVariable("k1hoo1",&m_k1hoo1);
  reader->AddVariable("k1hoo2",&m_k1hoo2);
  reader->AddVariable("k1hoo3",&m_k1hoo3);
  reader->AddVariable("k1hoo4",&m_k1hoo4);

  reader->BookMVA("BDT","weights/MVAnalysis_BDT.weights.xml");
  reader->BookMVA("BDTG","weights/MVAnalysis_BDTG.weights.xml");
  reader->BookMVA("BDTs","weights/MVAnalysis_sig_BDT.weights.xml");
  reader->BookMVA("BDTGs","weights/MVAnalysis_sig_BDTG.weights.xml");

  reader1->AddVariable("abs(cos_b0)",&m_cos_b0);
  reader1->AddVariable("p_ks",&m_p_ks);
  reader1->AddVariable("log(chi2_ndf_D0)",&m_chi2_ndf_D0);
  reader1->AddVariable("log(chi2_ndf_B0)",&m_chi2_ndf_B0);
  reader1->AddVariable("log(chi2_tag_vtx/ndf_tag_vtx)",&m_chi2_ndf_tag_vtx);
  reader1->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader1->AddVariable("abs(thr_sig-0.885)",&m_thr_sig);
  reader1->AddVariable("thr_oth",&m_thr_oth);
  reader1->AddVariable("log(tag_LH_err)",&m_tag_LH_err);
  reader1->AddVariable("log(dzerr)",&m_dzerr);
  reader1->AddVariable("log(pi0_chi2)",&m_pi0_chi2);
  reader1->AddVariable("log(egamma)",&m_egamma);
  reader1->AddVariable("lh",&m_lh);

  reader1->BookMVA("BDTlh","weights/MVAnalysis_lh_BDT.weights.xml");
  reader1->BookMVA("BDTGlh","weights/MVAnalysis_lh_BDTG.weights.xml");
  reader1->BookMVA("BDTlhs","weights/MVAnalysis_lh_sig_BDT.weights.xml");
  reader1->BookMVA("BDTGlhs","weights/MVAnalysis_lh_sig_BDTG.weights.xml");

  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);
//    if(!phsp) continue;
//    if(chi2_ndf_B0>1000) continue;

    m_cos_b0      = (float)TMath::Abs(cos_b0);
    m_p_ks        = (float)p_ks;
    m_p_pp        = (float)p_pp;
    m_p_pm        = (float)p_pm;
    m_p_pi0       = (float)p_pi0;
    m_chi2_ndf_D0 = (float)TMath::Log(chi2_ndf_D0);
    m_chi2_ndf_B0 = (float)TMath::Log(chi2_ndf_B0);
    m_chi2_ndf_tag_vtx = (float)TMath::Log(chi2_tag_vtx/ndf_tag_vtx);
    m_cos_thr     = (float)TMath::Abs(cos_thr);
    m_thr_sig     = (float)TMath::Abs(thr_sig-0.885);
    m_thr_oth     = (float)thr_oth;
    m_dzerr       = (float)TMath::Log(dzerr);
    m_ks_dr       = (float)ks_dr;
    m_ks_dz       = (float)ks_dz;
    m_ks_dphi     = (float)ks_dphi;
    m_ks_fl       = (float)ks_fl;
    m_pi0_chi2    = (float)TMath::Log(pi0_chi2);
    m_egamma      = (float)TMath::Log(egamma);
    m_ptgamma     = (float)ptgamma;
    m_k1mm2       = (float)k1mm2;
    m_k1et        = (float)k1et;
    m_k1hso00     = (float)k1hso00;
    m_k1hso02     = (float)k1hso02;
    m_k1hso04     = (float)k1hso04;
    m_k1hso10     = (float)k1hso10;
    m_k1hso12     = (float)k1hso12;
    m_k1hso14     = (float)k1hso14;
    m_k1hso20     = (float)k1hso20;
    m_k1hso22     = (float)k1hso22;
    m_k1hso24     = (float)k1hso24;
    m_k1hoo0      = (float)k1hoo0;
    m_k1hoo1      = (float)k1hoo1;
    m_k1hoo2      = (float)k1hoo2;
    m_k1hoo3      = (float)k1hoo3;
    m_k1hoo4      = (float)k1hoo4;
    m_lh          = 0;//(float)lh;

    k0vars[0]  = k0et;
    k0vars[1]  = k0hso00;
    k0vars[2]  = k0hso10;
    k0vars[3]  = k0hso20;
    k0vars[4]  = 0;//k0hso01;
    k0vars[5]  = k0hso02;
    k0vars[6]  = k0hso12;
    k0vars[7]  = k0hso22;
    k0vars[8]  = 0;//k0hso03;
    k0vars[9]  = k0hso04;
    k0vars[10] = k0hso14;
    k0vars[11] = k0hso24;
    k0vars[12] = k0hoo0;
    k0vars[13] = k0hoo1;
    k0vars[14] = k0hoo2;
    k0vars[15] = k0hoo3;
    k0vars[16] = k0hoo4;

    ksfw->input(k0mm2, k0vars);
    if(ksfw->ls() + ksfw->lb() > 0){
      lh = ksfw->ls()/(ksfw->ls()+ksfw->lb());
    } else{
      lh = -1;
    }

    bdt     = reader->EvaluateMVA("BDT");
    bdtg    = reader->EvaluateMVA("BDTG");
    bdts    = reader->EvaluateMVA("BDTs");
    bdtgs   = reader->EvaluateMVA("BDTGs");
    bdtlh   = reader1->EvaluateMVA("BDTlh");
    bdtglh  = reader1->EvaluateMVA("BDTGlh");
    bdtlhs  = reader1->EvaluateMVA("BDTlhs");
    bdtglhs = reader1->EvaluateMVA("BDTGlhs");

    TEvent->Fill();
  }

  TEvent->Write();
  ofile.Write();
  ofile.Close();
}

