#include "cuts.h"

void InD0Filter(const int type){
  TChain* tree = new TChain("TDEvent");
  string fname;
  switch(type){
  case 11://Signal MC pi0
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcPI0_s7.root");
//    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcPI0_s2.root");
    fname = string("fil_b2dh_sigmcPi0_s7.root");
    break;
  case 12://Signal MC eta
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcETA_DDALITZ_InclusiveD0_s1.root");
    fname = string("fil_b2dh_sigmcEta_DDALITZ_InclusiveD0_s1.root");
    break;
  case 13://Signal MC omega
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcOMEGA_s6.root");
    fname = string("fil_b2dh_sigmcOmega_s6.root");
    break;
  case 14://Signal MC rho
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcRHO_InclusiveD0_s1.root");
    fname = string("fil_b2dh_sigmcRho_InclusiveD0_s1.root");
    break;
  default:
    return;
  }
  Double_t mp,mm,z_sig,sz_sig,chisq_z_sig,p_d0,md,mk;
  Double_t mp_mc,mm_mc,d0_t_mc,t_sig_mc,z_sig_mc;
  Int_t exp,ntrk_sig,ndf_z_sig,bin;
  Int_t bin_mc,flv,mode,h0mode;
  Int_t good_icpv;

  tree->Print();

  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("ntrk_sig",&ntrk_sig);
  tree->SetBranchAddress("ndf_z_sig",&ndf_z_sig);
  tree->SetBranchAddress("bin",&bin);

  tree->SetBranchAddress("bin_mc",&bin_mc);
  tree->SetBranchAddress("flv",&flv);
  tree->SetBranchAddress("mode",&mode);
  tree->SetBranchAddress("h0mode",&h0mode);

  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("z_sig",&z_sig);
  tree->SetBranchAddress("sz_sig",&sz_sig);
  tree->SetBranchAddress("chisq_z_sig",&chisq_z_sig);
  tree->SetBranchAddress("p_d0",&p_d0);
  tree->SetBranchAddress("md0",&md);
  tree->SetBranchAddress("mks",&mk);

  tree->SetBranchAddress("mp_mc",&mp_mc);
  tree->SetBranchAddress("mm_mc",&mm_mc);
  tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
  tree->SetBranchAddress("t_sig_mc",&t_sig_mc);
  tree->SetBranchAddress("z_sig_mc",&z_sig_mc);

  TFile *ofile = new TFile(fname.c_str(),"RECREATE");
  fname = string("TEvent");

  TTree* TEvent = new TTree(fname.c_str(),fname.c_str());
  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("ntrk_sig",&ntrk_sig,"ntrk_sig/I");
  TEvent->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  TEvent->Branch("bin",&bin,"bin/I");

  TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
  TEvent->Branch("flv",&flv,"flv/I");
  TEvent->Branch("mode",&mode,"mode/I");
  TEvent->Branch("h0mode",&h0mode,"h0mode/I");

  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");
  TEvent->Branch("z_sig",&z_sig,"z_sig/D");
  TEvent->Branch("sz_sig",&sz_sig,"sz_sig/D");
  TEvent->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  TEvent->Branch("p_d0",&p_d0,"p_d0/D");
  TEvent->Branch("md0",&md,"md0/D");
  TEvent->Branch("mks",&mk,"mks/D");

  TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
  TEvent->Branch("mm_mc",&mm_mc,"mm_mc/D");
  TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
  TEvent->Branch("t_sig_mc",&t_sig_mc,"t_sig_mc/D");
  TEvent->Branch("z_sig_mc",&z_sig_mc,"z_sig_mc/D");

  TEvent->Branch("good_icpv",&good_icpv,"good_icpv/I");

  const int NTot = tree->GetEntries();
  for(int i=1; i<NTot; i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);

    if(md<(DMass-md_cut) || md>(DMass+md_cut)) continue;
    if(mk<(KMass-mk_cut) || mk>(KMass+mk_cut)) continue;

    z_sig *= 10;
    sz_sig = 10.*TMath::Sqrt(sz_sig);

    // * Standatd ICPV cuts * //
    if(sz_sig>0.5) good_icpv = 0;
    else           good_icpv = 1;
    // * ////////////////// * //

    TEvent->Fill();
  }
  TEvent->Print();
  TEvent->Write();
  ofile->Close();
  return;
}

