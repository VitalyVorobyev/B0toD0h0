#include "cuts.h"

using namespace std;

int d0_des = 0;
int h0_des = 0;
int totl_des = 0;
int good_des = 0;
int n_sig_cand = 0;
int n_bkg_cand = 0;
int cand_struct1[10];
int cand_struct2[10];

Int_t good_icpv, b0f, mode, h0mode, run, evtn, exp;
Double_t mbc,de,bdt,mh0,mpi0,mk,md,chi2_mass_d0;

int m_mode,m_h0mode;
double mh0_min, mh0_max;
double de_min, de_max;
double mbc_min, mbc_max;
double True_mh0;
double bdt_min;

double GetGoodTree(TTree* oldtree, TTree* newtree,const int type,const bool test){
  oldtree->SetBranchAddress("exp",&exp);
  oldtree->SetBranchAddress("run",&run);
  oldtree->SetBranchAddress("evtn",&evtn);
  oldtree->SetBranchAddress("mbc",&mbc);
  oldtree->SetBranchAddress("de",&de);
  oldtree->SetBranchAddress("b0f",&b0f);
  oldtree->SetBranchAddress("chi2_mass_d0",&chi2_mass_d0);
  oldtree->SetBranchAddress("bdt",&bdt);
  oldtree->SetBranchAddress("good_icpv",&good_icpv);
  oldtree->SetBranchAddress("md_raw",&md);
  oldtree->SetBranchAddress("mk",&mk);
  oldtree->SetBranchAddress("mh0",&mh0);
  oldtree->SetBranchAddress("mpi0",&mpi0);
  oldtree->SetBranchAddress("mode",&mode);
  oldtree->SetBranchAddress("h0mode",&h0mode);

  newtree->Branch("exp",&exp,"exp/I");
  newtree->Branch("run",&run,"run/I");
  newtree->Branch("evtn",&evtn,"evtn/I");
  newtree->Branch("md",&md,"md/D");
  newtree->Branch("mh0",&mh0,"mh0/D");
  newtree->Branch("b0f",&b0f,"b0f/I");

  const int NTot = oldtree->GetEntries();
  cout << "Events before cuts: " << NTot << endl;
  for(int i=0; i<NTot; i++){
    oldtree->GetEvent(i);
    if(b0f == 0 || b0f<-1) continue;
    if(test){if(b0f != 1 && b0f != 5 && b0f != 10) continue;}
    if(!good_icpv) continue;
    if(mode != m_mode || h0mode != m_h0mode) continue;
    if(md<md_min || md>md_max) continue;
    if(mbc<mbc_min || mbc>mbc_max) continue;
    if(de<de_min || de>de_max) continue;
    if(mh0<mh0_min || mh0>mh0_max) continue;
    if(type != 2 && type != 5){ if(mpi0<mpi0_min || mpi0>mpi0_max) continue;}
    if(chi2_mass_d0>50) continue;
    newtree->Fill();
  }
  return newtree->GetEntries();
}

bool is_decision(const vector<int>& b0fvec){
  int nsig = 0;
  int nbkg = 0;
  for(int i=0; i<b0fvec.size(); i++){
    const int flag = b0fvec[i];
    if(flag == 0){ continue;}
    if(flag == 1 || flag == 5 || flag == 10){ nsig++;
    } else{ nbkg++;}
  }
  if(nsig && nbkg){
    totl_des++;
    n_sig_cand += nsig;
    n_bkg_cand += nbkg;
    if((nsig+nbkg)<10-2){
      cand_struct1[nsig+nbkg-2]++;
      cand_struct2[nbkg+1-2]++;
    } else if((1+nbkg-2)<10){
      cand_struct1[9]++;
      cand_struct2[nbkg+1-2]++;
    } else{
      cand_struct1[9]++;
      cand_struct2[9]++;
    }
    return true;
  }
  return false;
}

int my_decision(const vector<double>& d0mass,const vector<double>& h0mass, const vector<int>& b0fvec){
  if(d0mass.size() != h0mass.size()){
    cout << "My decision: wrong sizes " << d0mass.size() << ", " << h0mass.size() << endl;
    return 0;
  }
  if(d0mass.size()<2){
    cout << "My decision: size = " << d0mass.size() << endl;
    return -1;
  }

  vector<double> mdmins,h0mins;
  vector<int> indexes;
  double mdmin = d0mass[0];
  double imin = 0;
  mdmins.clear(); h0mins.clear(); indexes.clear();
  mdmins.push_back(d0mass[0]);
  h0mins.push_back(h0mass[0]);
  indexes.push_back(0);
  for(int i=1; i<d0mass.size(); i++){
    if(d0mass[i]<mdmin){
      mdmin = d0mass[i];
      imin = i;
      mdmins.clear(); h0mins.clear(); indexes.clear();
      mdmins.push_back(d0mass[i]);
      h0mins.push_back(h0mass[i]);
      indexes.push_back(i);
    } else if(abs(d0mass[i] - mdmin) < 0.0001){
      mdmins.push_back(d0mass[i]);
      h0mins.push_back(h0mass[i]);
      indexes.push_back(i);
    }
  }
  if(h0mins.size() == 1){
    d0_des++;
    if(b0fvec.size()){
      if(is_decision(b0fvec)){
        const int flag = b0fvec[indexes[0]];
//        cout << flag << " ";
        if(flag == 1 || flag == 5 || flag == 10) good_des++;
      }
    }
    return indexes[0];
  } else{
    h0_des++;
  }
  double h0min = h0mins[0];
  double iimin = 0;
  for(int i=1; i<h0mins.size(); i++){
    if(h0mins[i]<h0min){
      h0min = h0mins[i];
      iimin = i;
    }
  }
  if(b0fvec.size()){
    if(is_decision(b0fvec)){
      const int flag = b0fvec[indexes[iimin]];
//      cout << flag << " ";
      if(flag == 1 || flag == 5 || flag == 10) good_des++;
    }
  }
  return indexes[iimin];
}

void multiplicity2(const int type,const bool test = false){
    for(int i=0; i<10; i++){
      cand_struct1[i] = 0;
      cand_struct2[i] = 0;
    }
    TChain* tree0 = new TChain("TEvent");
    switch(type){
    case 1:// pi0
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcPi0_s7.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigPi0_s7.root");
      True_mh0 = Pi0Mass;
      m_mode = 1; m_h0mode = 10;
      mh0_min = mpi0_min;
      mh0_max = mpi0_max;
      de_min = de_min_pi0;
      de_max = de_max_pi0;
      mbc_min = mbc_min_pi0;
      mbc_max = mbc_max_pi0;
      bdt_min = bdt_cut_pi0;
      break;
    case 2:// eta->gg
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcEta_s2.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigEta_s2.root");
      True_mh0 = EtaMass;
      m_mode = 2; m_h0mode = 10;
      mh0_min = metagg_min;
      mh0_max = metagg_max;
      de_min = de_min_etagg;
      de_max = de_max_etagg;
      mbc_min = mbc_min_etagg;
      mbc_max = mbc_max_etagg;
      bdt_min = bdt_cut_etagg;
      break;
    case 3:// eta->ppp
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcEta_s2.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigEta_s2.root");
      True_mh0 = EtaMass;
      m_mode = 2; m_h0mode = 20;
      mh0_min = metappp_min;
      mh0_max = metappp_max;
      de_min = de_min_etappp;
      de_max = de_max_etappp;
      mbc_min = mbc_min_etappp;
      mbc_max = mbc_max_etappp;
      bdt_min = bdt_cut_etappp;
      break;
    case 4:// omega
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcOmega_s5.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigOmega_s5.root");
      True_mh0 = OmegaMass;
      m_mode = 3; m_h0mode = 20;
      mh0_min = momega_min;
      mh0_max = momega_max;
      de_min = de_min_omega;
      de_max = de_max_omega;
      mbc_min = mbc_min_omega;
      mbc_max = mbc_max_omega;
      bdt_min = bdt_cut_omega;
      break;
    case 5://eta', eta->gg
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcETAP_s1.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigETAP_s1.root");
      True_mh0 = EtaMass;
      m_mode = 5; m_h0mode = 10;
      mh0_min = metagg_min;
      mh0_max = metagg_max;
      de_min = de_min_etagg;
      de_max = de_max_etagg;
      mbc_min = mbc_min_etagg;
      mbc_max = mbc_max_etagg;
      bdt_min = bdt_cut_etagg;
      break;
    case 6://eta', eta->ppp
      if(!test) tree0->Add("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcETAP_s1.root");
      else      tree0->Add("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigETAP_s1.root");
      True_mh0 = EtaMass;
      m_mode = 5; m_h0mode = 20;
      mh0_min = metappp_min;
      mh0_max = metappp_max;
      de_min = de_min_etappp;
      de_max = de_max_etappp;
      mbc_min = mbc_min_etappp;
      mbc_max = mbc_max_etappp;
      bdt_min = bdt_cut_etappp;
      break;
    default:
      cout << "Wrong type " << type << endl;
      return;
    }

  TTree *tree = new TTree("TEvent","TEvent");
  GetGoodTree(tree0,tree,type,test);
  const int NTot = tree->GetEntries();
  cout << "Events after cuts:  " << NTot << endl;

  if(test) return;
  int nevents = 0;
  int nrecords = 0;

  int my_des;
  vector<double> mdv;
  vector<double> h0v;
  vector<int> bflags;
  vector<int> recnum;
  tree->GetEvent(0);
  mdv.push_back(TMath::Abs(md-DMass));
  h0v.push_back(TMath::Abs(mh0-True_mh0));
  recnum.push_back(0);
  bflags.push_back(b0f);
  int cur_evtn = evtn;
  int cur_run = run;
  for(int i=1; i<(NTot); i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);
    nrecords++;
    if(cur_evtn == evtn && cur_run == run){
      mdv.push_back(TMath::Abs(md-DMass));
      h0v.push_back(TMath::Abs(mh0-True_mh0));
      recnum.push_back(i);
      if(type) bflags.push_back(b0f);
    } else{
      cur_evtn = evtn; cur_run = run;
      nevents++;
      if(mdv.size() > 1){
        my_des = my_decision(mdv,h0v,bflags);
        tree->GetEvent(recnum[my_des]);
        tree->GetEvent(i);
      } else{
        tree->GetEvent(i-1);
        tree->GetEvent(i);
      }
      mdv.clear(); h0v.clear();
      recnum.clear(); bflags.clear();
      mdv.push_back(TMath::Abs(md-DMass));
      h0v.push_back(TMath::Abs(mh0-True_mh0));
      recnum.push_back(i);
      if(type) bflags.push_back(b0f);
    }
  }

  cout << "Multiplicity = " << nrecords << "/" << nevents << " = " << (double)nrecords/(double)nevents << " " << endl;
  cout << "  D0 decisions: " << d0_des << endl;
  cout << "  h0 decisions: " << h0_des << endl;
  if(type){
    if(totl_des){
      cout << "True mult: = " << (n_bkg_cand+n_sig_cand) << "/" << totl_des << " = " << (double)(n_bkg_cand+n_sig_cand)/totl_des << endl;
      cout << "Good decisions: " << good_des << "/" << totl_des << " = " << (double)good_des/totl_des << endl;
    }
    cout << "Multiplisity structure:" << endl;
    for(int i=0; i<10; i++){
      cout << "  " << i+1 << ": " << cand_struct1[i] << " (" << cand_struct2[i] << ")" << endl;
    }
  }

  return;
}
