#include "cuts.h"

using namespace std;

int min_decision(const vector<double> v){
  double vmin = v[0];
  double imin = 0;
  for(int i=1; i<v.size(); i++){
    if(v[i]<vmin){
      vmin = v[i];
      imin = i;
    }
  }
  return imin;
}

int max_decision(const vector<double> v){
  double vmax = v[0];
  double imax = 0;
  for(int i=1; i<v.size(); i++){
    if(v[i]>vmax){
      vmax = v[i];
      imax = i;
    }
  }
  return imax;
}

void multiplicity_data(void){
  TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/Tuples/b2dh_cont.root");
  TTree *tree = (TTree*)ifile->Get("TEvent");

  Double_t cos_b0,p_ks,p_pp,p_pm,p_pi0,chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx,cos_thr,thr_sig,thr_oth;
  Int_t ndf_tag_vtx,phsp,bin,exp,run,evtn,b0f;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;

  Double_t mbc,de,bdt,bdtg,mp,mm,dz,dt,atckpi_max,mpi0,mk,md;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err,dzerr;

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("md0_raw",&md);
  tree->SetBranchAddress("mpi0_raw",&mpi0);
  tree->SetBranchAddress("mks_raw",&mk);
  tree->SetBranchAddress("chi2_ndf_D0",&chi2_ndf_D0);
  tree->SetBranchAddress("chi2_ndf_B0",&chi2_ndf_B0);
  tree->SetBranchAddress("cos_thr",&cos_thr);

  double mult = 0,totsig = 0,mult_tot = 0;
  const int NTot = tree->GetEntries();
  vector<double> mbcv;
  int mbc_des,bchi2_des,dchi2_des,bdchi2_des,thr_des,rndm_des,true_des;
  int cur_evtn,EVTN;
  int cur_run,RUN;
  tree->GetEvent(0);
  mbcv.push_back(mbc);
  cur_evtn = evtn;
  for(int i=1; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);
    if(md<(DMass-md_cut) || md>(DMass+md_cut)) continue;
    if(mpi0<(Pi0Mass-mpi0_cut) || mpi0>(Pi0Mass+mpi0_cut)) continue;
    if(mk<(KMass-mk_cut) || mk>(KMass+mk_cut)) continue;
    if(cur_evtn == evtn && cur_run == run){ mbcv.push_back(mbc);}
    else{
      totsig += mbcv.size();
      mult++;
      mbcv.clear();

      mbcv.push_back(mbc);
      cur_evtn = evtn; cur_run = run;
    }
  }

  cout << "Totsig = " << totsig << endl;
  cout << "Mult   = " << mult << endl;
  if(mult) cout << "Multiplicity = " << (double)totsig/(double)mult << endl;
}
