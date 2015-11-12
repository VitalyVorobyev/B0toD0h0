#include "ksfw.h"
#include "rooksfw.h"
#include "rooksfw.cc"
#include "../BtoDh/cuts.h"

using namespace std;

void ksfw_selection(const int type){
  rooksfw* ksfw  = new rooksfw("ksfw0", ksfw0_alpha, ksfw0_sigpdf, ksfw0_bkgpdf);
  TFile *ifile;
  if(type) ifile = TFile::Open("/home/vitaly/B0toDh0/Tuples/test_cont.root");
  else  ifile = TFile::Open("/home/vitaly/B0toDh0/Tuples/b2dh_sigmc_s1.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  string fname;
  if(type) fname = string("cont_w_lh.root");
  else fname = string("sig_w_lh.root");
  TFile ofile(fname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  Double_t cos_b0,p_ks,;
  Double_t chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx;
  Int_t ndf_tag_vtx;
  Double_t cos_thr,thr_sig,thr_oth,dzerr,tag_LH_err;
  Double_t pi0_chi2,egamma,ptgamma;

  Double_t mbc,de,mp,mm,dz,dt,atckpi_max,mpi0,mk,md;
  Double_t lh,bdt,bdtg;
  Int_t phsp;

  double k0mm2,k0vars[17];
  Int_t b0f,d0f;

  tree->SetBranchAddress("b0f",&b0f);
  tree->SetBranchAddress("d0f",&d0f);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("dt",&dt);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("phsp",&phsp);
  tree->SetBranchAddress("mpi0_raw",&mpi0);
  tree->SetBranchAddress("mks_raw",&mk);
  tree->SetBranchAddress("md0_raw",&md);

  tree->SetBranchAddress("cos_b0",&cos_b0);
  tree->SetBranchAddress("p_ks",&p_ks);
  tree->SetBranchAddress("chi2_ndf_D0",&chi2_ndf_D0);
  tree->SetBranchAddress("chi2_ndf_B0",&chi2_ndf_B0);
  tree->SetBranchAddress("chi2_tag_vtx",&chi2_tag_vtx);
  tree->SetBranchAddress("ndf_tag_vtx",&ndf_tag_vtx);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);
  tree->SetBranchAddress("dzerr",&dzerr);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("pi0_chi2",&pi0_chi2);
  tree->SetBranchAddress("egamma",&egamma);
  tree->SetBranchAddress("ptgamma",&ptgamma);

  tree->SetBranchAddress("k0mm2",   &k0mm2);
  tree->SetBranchAddress("k0et",    &k0vars[0]);
  tree->SetBranchAddress("k0hso00", &k0vars[1]);
  tree->SetBranchAddress("k0hso10", &k0vars[2]);
  tree->SetBranchAddress("k0hso20", &k0vars[3]);
  tree->SetBranchAddress("k0hso01", &k0vars[4]);
  tree->SetBranchAddress("k0hso02", &k0vars[5]);
  tree->SetBranchAddress("k0hso12", &k0vars[6]);
  tree->SetBranchAddress("k0hso22", &k0vars[7]);
  tree->SetBranchAddress("k0hso03", &k0vars[8]);
  tree->SetBranchAddress("k0hso04", &k0vars[9]);
  tree->SetBranchAddress("k0hso14", &k0vars[10]);
  tree->SetBranchAddress("k0hso24", &k0vars[11]);
  tree->SetBranchAddress("k0hoo0",  &k0vars[12]);
  tree->SetBranchAddress("k0hoo1",  &k0vars[13]);
  tree->SetBranchAddress("k0hoo2",  &k0vars[14]);
  tree->SetBranchAddress("k0hoo3",  &k0vars[15]);
  tree->SetBranchAddress("k0hoo4",  &k0vars[16]);

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");

  TEvent->Branch("d0f",&d0f,"d0f/I");
  TEvent->Branch("b0f",&b0f,"b0f/I");
  TEvent->Branch("phsp",&phsp,"phsp/I");

  TEvent->Branch("cos_b0",&cos_b0,"cos_b0/D");
  TEvent->Branch("cos_thr",&cos_thr,"cos_thr/D");
  TEvent->Branch("p_ks",&p_ks,"p_ks/D");
  TEvent->Branch("chi2_ndf_D0",&chi2_ndf_D0,"chi2_ndf_D0/D");
  TEvent->Branch("chi2_ndf_B0",&chi2_ndf_B0,"chi2_ndf_B0/D");
  TEvent->Branch("chi2_tag_vtx",&chi2_tag_vtx,"chi2_tag_vtx/D");
  TEvent->Branch("ndf_tag_vtx",&ndf_tag_vtx,"ndf_tag_vtx/I");
  TEvent->Branch("thr_sig",&thr_sig,"thr_sig/D");
  TEvent->Branch("thr_oth",&thr_oth,"thr_oth/D");
  TEvent->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");
  TEvent->Branch("dzerr",&dzerr,"dzerr/D");
  TEvent->Branch("pi0_chi2",&pi0_chi2,"pi0_chi2/D");
  TEvent->Branch("egamma",&egamma,"egamma/D");
  TEvent->Branch("ptgamma",&ptgamma,"ptgamma/D");

  TEvent->Branch("lh",&lh,"lh/D");

  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);
    
    if(!type){
      if(b0f != 1 && b0f != 5 && b0f != 10) continue;
    }
    else{
      if(!(b0f != 1 && b0f != 10 && !(d0f == 1 && (b0f == 5 || b0f == 4)))) continue;
    }
    if(md<(DMass-md_cut) || md>(DMass+md_cut)) continue;
    if(mpi0<(Pi0Mass-mpi0_cut) || mpi0>(Pi0Mass+mpi0_cut)) continue;
    if(mk<(KMass-mk_cut) || mk>(KMass+mk_cut)) continue;

    ksfw->input(k0mm2, k0vars);
    if(ksfw->ls() + ksfw->lb() > 0){
      lh = ksfw->ls() / (ksfw->ls() + ksfw->lb());
    } else{
      lh = -1;
    }
    TEvent->Fill();
  }

  TEvent->Write();
  ofile.Write();
  ofile.Close();
}

