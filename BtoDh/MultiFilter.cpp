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

int IsGoodICPV(const int ndf_z_sig, const double& sz_sig, const double& chisq_z_sig,const int ndf_z_asc, const double& sz_asc, const double& chisq_z_asc){
  // * Standatd ICPV cuts * //
  if(ndf_z_sig == 0 && sz_sig>0.5 || ndf_z_sig > 0 && (sz_sig>0.2 || chisq_z_sig/ndf_z_sig>50)){
    return 0;
  } else if(ndf_z_asc == 0 && sz_asc>0.5 || ndf_z_asc > 0 && (sz_asc>0.2 || chisq_z_asc/ndf_z_asc>50)){
    return 0;
  } else{
    return 1;
  }
  // * ////////////////// * //
}

int max_decision(const vector<double> v){
  if(v.size()<2){
    cout << "max_decision: size = " << v.size() << endl;
    return -1;
  }
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

int min_decision(const vector<double>& v){
//  if(v.size()<2){
//    cout << "min_decision: size = " << v.size() << endl;
//    return -1;
//  }
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

bool is_decision(const vector<int>& b0fvec){
  int nsig = 0;
  int nbkg = 0;
  for(int i=0; i<b0fvec.size(); i++){
    const int flag = b0fvec[i];
//    cout << flag << " ";
    if(flag == 0){ continue;}
    if(flag == 1 || flag == 5 || flag == 10){ nsig++;
    } else{ nbkg++;}
  }
//  cout << endl;
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

void MultiFilter(const int type){
  for(int i=0; i<10; i++){
    cand_struct1[i] = 0;
    cand_struct2[i] = 0;
  }
  string ofname;
  TFile *ifile;
  switch(type){
  case 0:// data
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_data.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_data.root");
    break;
  case 11:// pi0
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcPi0_s7.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigPi0_s7.root");
    break;
  case 12:// eta
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcEta_s2.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigEta_s2.root");
    break;
  case 13:// omega
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcOmega_s5.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_sigOmega_s5.root");
    break;
  case 14:// rho
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcRho_s1.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigRho_s1.root");
    break;
  case 15:// eta'
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_sigmcETAP_s1.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigETAP_s1.root");
    break;
  case 2:// gen
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_gen_0-1_full.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_gen_0-1.root");
    break;
  case 21:// charged
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_charged_2_12.root");
    ofname =     string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_charged_2_12.root");
    break;
  case 22:// mixed
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_mixed_2_12.root");
    ofname =     string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_mixed_2_12.root");
    break;
  case 23:// charm
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_charm_2_12.root");
    ofname =     string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_charm_2_12.root");
    break;
  case 24:// uds
    ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_uds_2_12.root");
    ofname =     string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_uds_2_12.root");
    break;
  case 3:// udsc
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil1_b2dh_cont_0-1.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL1_b2dh_cont_0-1.root");
    break;
  case 4:// B+-0
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_bb_0-1_full.root");
    ofname =      string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_bb_0-1.root");
    break;
  default:
    cout << "Wrong type " << type << endl;
    return;
  }

  TTree *tree = (TTree*)ifile->Get("TEvent");
  tree->Print();

  TFile* ofile = new TFile(ofname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  //Exp, run,evtn
  Int_t exp,run,evtn,good_icpv,good_icpv_d0;

  tree->SetBranchAddress("good_icpv",&good_icpv);
//  tree->SetBranchAddress("good_icpv_d0",&good_icpv_d0);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("exp",&exp);

  TEvent->Branch("good_icpv",&good_icpv,"good_icpv/I");
//  TEvent->Branch("good_icpv_d0",&good_icpv_d0,"good_icpv_d0/I");
  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");
  TEvent->Branch("exp",&exp,"exp/I");

  //Tatami variables
  Double_t chisq_z_sig,chisq_z_asc,cl_z_sig,cl_z_asc;
  Int_t ntrk_asc,ntrk_sig,ndf_z_sig,ndf_z_asc;

  tree->SetBranchAddress("ntrk_sig",&ntrk_sig);
  tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  tree->SetBranchAddress("ndf_z_sig",&ndf_z_sig);
  tree->SetBranchAddress("ndf_z_asc",&ndf_z_asc);
  tree->SetBranchAddress("chisq_z_sig",&chisq_z_sig);
  tree->SetBranchAddress("chisq_z_asc",&chisq_z_asc);
  tree->SetBranchAddress("cl_z_sig",&cl_z_sig);
  tree->SetBranchAddress("cl_z_asc",&cl_z_asc);

  TEvent->Branch("ntrk_sig",&ntrk_sig,"ntrk_sig/I");
  TEvent->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  TEvent->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  TEvent->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");
  TEvent->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  TEvent->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
  TEvent->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
  TEvent->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");

  Double_t px_pim,py_pim,pz_pim;
  Double_t px_pip,py_pip,pz_pip;
  Double_t px_ks,py_ks,pz_ks;

  tree->SetBranchAddress("px_pim",&px_pim);
  tree->SetBranchAddress("py_pim",&py_pim);
  tree->SetBranchAddress("pz_pim",&pz_pim);
  tree->SetBranchAddress("px_pip",&px_pip);
  tree->SetBranchAddress("py_pip",&py_pip);
  tree->SetBranchAddress("pz_pip",&pz_pip);
  tree->SetBranchAddress("px_ks",&px_ks);
  tree->SetBranchAddress("py_ks",&py_ks);
  tree->SetBranchAddress("pz_ks",&pz_ks);

  TEvent->Branch("px_pim",&px_pim,"px_pim/D");
  TEvent->Branch("py_pim",&py_pim,"py_pim/D");
  TEvent->Branch("pz_pim",&pz_pim,"pz_pim/D");
  TEvent->Branch("px_pip",&px_pip,"px_pip/D");
  TEvent->Branch("py_pip",&py_pip,"py_pip/D");
  TEvent->Branch("pz_pip",&pz_pip,"pz_pip/D");
  TEvent->Branch("px_ks",&px_ks,"px_ks/D");
  TEvent->Branch("py_ks",&py_ks,"py_ks/D");
  TEvent->Branch("pz_ks",&pz_ks,"pz_ks/D");

  //Continuum suppression
  Double_t bdtg,bdt,lh1,cos_thr,egamma,cos_hel;
  tree->SetBranchAddress("bdtg",&bdtg);
  TEvent->Branch("bdtg",&bdtg,"bdtg/D");

  tree->SetBranchAddress("bdt",&bdt);
  TEvent->Branch("bdt",&bdt,"bdt/D");

  tree->SetBranchAddress("lh1",&lh1);
  TEvent->Branch("lh1",&lh1,"lh1/D");

  tree->SetBranchAddress("cos_thr",&cos_thr);
  TEvent->Branch("cos_thr",&cos_thr,"cos_thr/D");

  tree->SetBranchAddress("cos_hel",&cos_hel);
  TEvent->Branch("cos_hel",&cos_hel,"cos_hel/D");

  tree->SetBranchAddress("egamma",&egamma);
  TEvent->Branch("egamma",&egamma,"egamma/D");

  //Tagging
  Double_t tag_LH,tag_LH_err;
  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);

  TEvent->Branch("tag_LH",&tag_LH,"tag_LH/D");
  TEvent->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");

  //Dalitz variables
  Int_t bin;
  Double_t mp,mm;
  tree->SetBranchAddress("bin",&bin);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);

  TEvent->Branch("bin",&bin,"bin/I");
  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");

  //MC
  //Double_t dz_pull_sig, dz_pull_asc;
  Double_t dz_mc_sig, dz_mc_asc, t_sig_mc, z_sig_mc, t_asc_mc, z_asc_mc;
  Double_t z_sig,z_asc,sz_sig,sz_asc;
  Int_t flv_mc, mode, h0mode;
  Int_t bin_mc;
  Double_t mp_mc,mm_mc,d0_t_mc,dt_mc,dz_mc;
  Double_t mp_raw,mm_raw;

  tree->SetBranchAddress("mode",&mode);
  tree->SetBranchAddress("h0mode",&h0mode);
  TEvent->Branch("mode",&mode,"mode/I");
  TEvent->Branch("h0mode",&h0mode,"h0mode/I");

  tree->SetBranchAddress("z_sig",&z_sig);
  TEvent->Branch("z_sig",&z_sig,"z_sig/D");
  tree->SetBranchAddress("z_asc",&z_asc);
  TEvent->Branch("z_asc",&z_asc,"z_asc/D");

  tree->SetBranchAddress("sz_sig",&sz_sig);
  TEvent->Branch("sz_sig",&sz_sig,"sz_sig/D");
  tree->SetBranchAddress("sz_asc",&sz_asc);
  TEvent->Branch("sz_asc",&sz_asc,"sz_asc/D");

//  Double_t dz_mc_sig_d0,z_sig_d0,sz_sig_d0,dz_pull_sig_d0;
//
//  if(type == 11 || type == 12 || type == 13 || type == 14){
//    tree->SetBranchAddress("dz_mc_sig_d0",&dz_mc_sig_d0);
//    tree->SetBranchAddress("z_sig_d0",&z_sig_d0);
//    tree->SetBranchAddress("sz_sig_d0",&sz_sig_d0);
//    tree->SetBranchAddress("dz_pull_sig_d0",&dz_pull_sig_d0);
//
//    TEvent->Branch("dz_mc_sig_d0",&dz_mc_sig_d0,"dz_mc_sig_d0/D");
//    TEvent->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
//    TEvent->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
//    TEvent->Branch("dz_pull_sig_d0",&dz_pull_sig_d0,"dz_pull_sig_d0/D");
//  }

  Int_t b0id,d0id,h0id,dst0id,dst0f,etapid,etapf,d0_flv_mc;

  if(type){
    Int_t b0f,d0f,h0f,pi0f;
    tree->SetBranchAddress("b0id",&b0id);
    tree->SetBranchAddress("b0f",&b0f);
    tree->SetBranchAddress("d0id",&d0id);
    tree->SetBranchAddress("d0f",&d0f);
    tree->SetBranchAddress("h0id",&h0id);
    tree->SetBranchAddress("h0f",&h0f);
    tree->SetBranchAddress("pi0f",&pi0f);
    tree->SetBranchAddress("dst0id",&dst0id);
    tree->SetBranchAddress("dst0f",&dst0f);
    tree->SetBranchAddress("etapid",&etapid);
    tree->SetBranchAddress("etapf",&etapf);

    TEvent->Branch("b0id",&b0id,"b0id/I");
    TEvent->Branch("b0f",&b0f,"b0f/I");
    TEvent->Branch("d0id",&d0id,"d0id/I");
    TEvent->Branch("d0f",&d0f,"d0f/I");
    TEvent->Branch("h0id",&h0id,"h0id/I");
    TEvent->Branch("h0f",&h0f,"h0f/I");
    TEvent->Branch("pi0f",&pi0f,"pi0f/I");
    TEvent->Branch("dst0id",&dst0id,"dst0id/I");
    TEvent->Branch("dst0f",&dst0f,"dst0f/I");
    TEvent->Branch("etapid",&etapid,"etapid/I");
    TEvent->Branch("etapf",&etapf,"etapf/I");

    if(type>10 && type<20){
      tree->SetBranchAddress("bin_mc",&bin_mc);
      tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
      tree->SetBranchAddress("mp_mc",&mp_mc);
      tree->SetBranchAddress("mm_mc",&mm_mc);
      tree->SetBranchAddress("mp_raw",&mp_raw);
      tree->SetBranchAddress("mm_raw",&mm_raw);
      tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
      tree->SetBranchAddress("dt_mc",&dt_mc);
      tree->SetBranchAddress("dz_mc",&dz_mc);

      tree->SetBranchAddress("dz_mc_sig",&dz_mc_sig);
      tree->SetBranchAddress("dz_mc_asc",&dz_mc_asc);
//      tree->SetBranchAddress("dz_pull_sig",&dz_pull_sig);
//      tree->SetBranchAddress("dz_pull_asc",&dz_pull_asc);
      tree->SetBranchAddress("flv_mc",&flv_mc);

      tree->SetBranchAddress("t_sig_mc",&t_sig_mc);
      tree->SetBranchAddress("z_sig_mc",&z_sig_mc);
      tree->SetBranchAddress("t_asc_mc",&t_asc_mc);
      tree->SetBranchAddress("z_asc_mc",&z_asc_mc);

      TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
      TEvent->Branch("flv_mc",&flv_mc,"flv_mc/I");
      TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
      TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
      TEvent->Branch("mm_mc",&mm_mc,"mm_mc/D");
      TEvent->Branch("mp_raw",&mp_raw,"mp_raw/D");
      TEvent->Branch("mm_raw",&mm_raw,"mm_raw/D");
      TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
      TEvent->Branch("dt_mc",&dt_mc,"dt_mc/D");
      TEvent->Branch("dz_mc",&dz_mc,"dz_mc/D");

      TEvent->Branch("dz_mc_sig",&dz_mc_sig,"dz_mc_sig/D");
      TEvent->Branch("dz_mc_asc",&dz_mc_asc,"dz_mc_asc/D");
//      TEvent->Branch("dz_pull_sig",&dz_pull_sig,"dz_pull_sig/D");
//      TEvent->Branch("dz_pull_asc",&dz_pull_asc,"dz_pull_asc/D");

      TEvent->Branch("t_sig_mc",&t_sig_mc,"t_sig_mc/D");
      TEvent->Branch("z_sig_mc",&z_sig_mc,"z_sig_mc/D");
      TEvent->Branch("t_asc_mc",&t_asc_mc,"t_asc_mc/D");
      TEvent->Branch("z_asc_mc",&z_asc_mc,"z_asc_mc/D");
    }
  }
  Int_t nptag;
  tree->SetBranchAddress("nptag",&nptag);
  TEvent->Branch("nptag",&nptag,"nptag/I");

  //Kinematic variables
  Double_t mbc,de,atckpi_max,mpi0,mh0,mk,costhBcms;
  Double_t md, md_raw, md_fit, mdpip, mdpim, p_d0;
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("mpi0",&mpi0);
  tree->SetBranchAddress("mh0",&mh0);
  tree->SetBranchAddress("mk",&mk);
  tree->SetBranchAddress("md_raw",&md_raw);
  tree->SetBranchAddress("md_fit",&md_fit);
  tree->SetBranchAddress("md",&md);
  tree->SetBranchAddress("mdpip",&mdpip);
  tree->SetBranchAddress("mdpim",&mdpim);
  tree->SetBranchAddress("p_d0",&p_d0);
  tree->SetBranchAddress("nptag",&nptag);
  tree->SetBranchAddress("costhBcms",&costhBcms);

  Double_t mdst0, metap, dmdst0, dmetap;
  tree->SetBranchAddress("mdst0",&mdst0);
  tree->SetBranchAddress("dmdst0",&dmdst0);
  tree->SetBranchAddress("metap",&metap);
  tree->SetBranchAddress("dmetap",&dmetap);

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");
  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("mpi0",&mpi0,"mpi0/D");
  TEvent->Branch("mh0",&mh0,"mh0/D");
  TEvent->Branch("mk",&mk,"mk/D");
  TEvent->Branch("md",&md,"md/D");
  TEvent->Branch("md_raw",&md_raw,"md_raw/D");
  TEvent->Branch("md_fit",&md_fit,"md_fit/D");
  TEvent->Branch("mdpip",&mdpip,"mdpip/D");
  TEvent->Branch("mdpim",&mdpim,"mdpim/D");
  TEvent->Branch("costhBcms",&costhBcms,"costhBcms/D");

  TEvent->Branch("mdst0",&mdst0,"mdst0/D");
  TEvent->Branch("dmdst0",&dmdst0,"dmdst0/D");
  TEvent->Branch("metap",&metap,"metap/D");
  TEvent->Branch("dmetap",&dmetap,"dmetap/D");

  //Time
  Double_t dz;
  tree->SetBranchAddress("dz",&dz);
  TEvent->Branch("dz",&dz,"dz/D");

  Double_t chi2_vtx_d0, chi2_mass_d0;
  tree->SetBranchAddress("chi2_vtx_d0",&chi2_vtx_d0);
  TEvent->Branch("chi2_vtx_d0",&chi2_vtx_d0,"chi2_vtx_d0/D");

  tree->SetBranchAddress("chi2_mass_d0",&chi2_mass_d0);
  TEvent->Branch("chi2_mass_d0",&chi2_mass_d0,"chi2_mass_d0/D");

  int nevents = 0;
  int nrecords = 0;

  const int NTot = tree->GetEntries();
  int my_des;
  vector<double> mdv;
  vector<double> h0v;
  vector<int> bflags;
  vector<int> recnum;
  tree->GetEvent(0);
  if(type == 2 || type == 3){
    bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
  }
  mdv.push_back(TMath::Abs(md-DMass));
  switch(mode){
  case 1:
    h0v.push_back(TMath::Abs(mh0-Pi0Mass));
    break;
  case 10:
    h0v.push_back(TMath::Abs(mh0-Pi0Mass));
    break;
  case 3:
    h0v.push_back(TMath::Abs(mh0-OmegaMass));
    break;
  case 30:
    h0v.push_back(TMath::Abs(mh0-OmegaMass));
    break;
  case 4:
    h0v.push_back(TMath::Abs(mh0-RhoMass));
    break;
  default:
    h0v.push_back(TMath::Abs(mh0-EtaMass));
    break;
  }
  recnum.push_back(0);
  if(type) bflags.push_back(b0f);
  int cur_evtn = evtn;
  int cur_run = run;
  for(int i=1; i<(NTot); i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);
//    if(type == 2 && (b0f == 1 || b0f == 5 || b0f == 10)) continue;
    if(type == 2 || type == 3){
      bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
    }
    nrecords++;
    if(cur_evtn == evtn && cur_run == run){
      mdv.push_back(TMath::Abs(md-DMass));
      switch(mode){
      case 1:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 10:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 3:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      case 30:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      case 4:
        h0v.push_back(TMath::Abs(mh0-RhoMass));
        break;
      default:
        h0v.push_back(TMath::Abs(mh0-EtaMass));
        break;
      }
      recnum.push_back(i);
      if(type) bflags.push_back(b0f);
    } else{
      cur_evtn = evtn; cur_run = run;
      nevents++;
      if(mdv.size() > 1){
          my_des = my_decision(mdv,h0v,bflags);
//        my_des = min_decision(mdv);
        tree->GetEvent(recnum[my_des]);
        if(type == 2 || type == 3){
          bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
        }
        TEvent->Fill();
        tree->GetEvent(i);
        if(type == 2 || type == 3){
          bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
        }
      } else{
        tree->GetEvent(i-1);
        if(type == 2 || type == 3){
          bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
        }
        TEvent->Fill();
        tree->GetEvent(i);
        if(type == 2 || type == 3){
          bin_mc = bin; flv_mc = tag_LH > 0 ? 1 : -1;
        }
      }
      mdv.clear(); h0v.clear();
      recnum.clear(); bflags.clear();
      mdv.push_back(TMath::Abs(md-DMass));
      switch(mode){
      case 1:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 10:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 3:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      case 30:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      case 4:
        h0v.push_back(TMath::Abs(mh0-RhoMass));
        break;
      default:
        h0v.push_back(TMath::Abs(mh0-EtaMass));
        break;
      }
      recnum.push_back(i);
      if(type) bflags.push_back(b0f);
    }
  }

  TEvent->Print();
  TEvent->Write();
  ofile->Close();

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
