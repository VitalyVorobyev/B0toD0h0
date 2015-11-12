#include "../BtoDh/cuts.h"

using namespace std;

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

int my_decision(const vector<double>& d0mass,const vector<double>& h0mass){
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
    } else if(d0mass[i] == mdmin){
      mdmins.push_back(d0mass[i]);
      h0mins.push_back(h0mass[i]);
      indexes.push_back(i);
    }
  }
//  cout << "h0mins.size() = " << h0mins.size() << endl;
  if(h0mins.size() == 1) return indexes[0];

  double h0min = h0mins[0];
  double iimin = 0;
  for(int i=1; i<h0mins.size(); i++){
    if(h0mins[i]<h0min){
      h0min = h0mins[i];
      iimin = i;
    }
  }
  return indexes[iimin];
}

void MultiFilter(const int type){
  string ofname;
  TFile *ifile;
  switch(type){
    case 0:// Data
      ifile  = TFile::Open("Fil_b2dpi_data_v2.root");
      ofname =      string("FIL_b2dpi_data_v2.root");
      break;
    case 1:// Charged
      ifile  = TFile::Open("Fil_b2dpi_charged_v2_0_10.root");
      ofname =      string("FIL_b2dpi_charged_v2_0_10.root");
      break;
    case 2:// Charm
      ifile  = TFile::Open("Fil_b2dpi_charm_v2_0_10.root");
      ofname =      string("FIL_b2dpi_charm_v2_0_10.root");
      break;
    case 3:// Charm
      ifile  = TFile::Open("Fil_bp2d0pip_sigmc.root");
      ofname =      string("FIL_bp2d0pip_sigmc.root");
      break;
    default:
      cout << "Wrong type " << type << endl;
      return;
  }

  TTree *tree = (TTree*)ifile->Get("TEvent");

  TFile* ofile = new TFile(ofname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  //Exp, run,evtn
  Int_t exp,run,evtn,flv;

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("flv",&flv);

  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");
  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("flv",&flv,"flv/I");

  //Tatami variables
  Double_t z_sig,z_asc,z_sig_d0,sz_sig,sz_sig_d0,sz_asc;
  Double_t chisq_z_sig,chisq_z_asc,cl_z_sig,cl_z_asc;
  Int_t ntrk_asc,ndf_z_sig,ndf_z_asc;
  Int_t good_icpv_mlt, good_icpv_sgl;

  tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  tree->SetBranchAddress("ndf_z_sig",&ndf_z_sig);
  tree->SetBranchAddress("ndf_z_asc",&ndf_z_asc);
  tree->SetBranchAddress("chisq_z_sig",&chisq_z_sig);
  tree->SetBranchAddress("chisq_z_asc",&chisq_z_asc);
  tree->SetBranchAddress("cl_z_sig",&cl_z_sig);
  tree->SetBranchAddress("cl_z_asc",&cl_z_asc);
  tree->SetBranchAddress("good_icpv_mlt",&good_icpv_mlt);
  tree->SetBranchAddress("good_icpv_sgl",&good_icpv_sgl);

  TEvent->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  TEvent->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  TEvent->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");
  TEvent->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  TEvent->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
  TEvent->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
  TEvent->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");
  TEvent->Branch("good_icpv_mlt",&good_icpv_mlt,"good_icpv_mlt/I");
  TEvent->Branch("good_icpv_sgl",&good_icpv_sgl,"good_icpv_sgl/I");

  //Continuum suppression
  Double_t bdtg;
  tree->SetBranchAddress("bdtg",&bdtg);

  TEvent->Branch("bdtg",&bdtg,"bdtg/D");

  //Dalitz variables
  Int_t bin;
  Double_t mp,mm;
  tree->SetBranchAddress("bin",&bin);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);

  TEvent->Branch("bin",&bin,"bin/I");
  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");

  tree->SetBranchAddress("z_sig",&z_sig);
  TEvent->Branch("z_sig",&z_sig,"z_sig/D");
  tree->SetBranchAddress("z_sig_d0",&z_sig_d0);
  TEvent->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
  tree->SetBranchAddress("z_asc",&z_asc);
  TEvent->Branch("z_asc",&z_asc,"z_asc/D");

  tree->SetBranchAddress("sz_sig",&sz_sig);
  TEvent->Branch("sz_sig",&sz_sig,"sz_sig/D");
  tree->SetBranchAddress("sz_sig_d0",&sz_sig_d0);
  TEvent->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
  tree->SetBranchAddress("sz_asc",&sz_asc);
  TEvent->Branch("sz_asc",&sz_asc,"sz_asc/D");

  Int_t nptag,bin_mc,flv_asc,d0_flv_mc;
  Double_t mp_mc,mm_mc,d0_t_mc,z_mc_sig,z_mc_asc;

  if(type == 3){
    tree->SetBranchAddress("mp_mc",&mp_mc);
    tree->SetBranchAddress("mm_mc",&mm_mc);
    tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
    tree->SetBranchAddress("z_sig_mc",&z_mc_sig);
    tree->SetBranchAddress("z_asc_mc",&z_mc_asc);
    tree->SetBranchAddress("bin_mc",&bin_mc);
    tree->SetBranchAddress("flv_mc",&flv_asc);
    tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);

    TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
    TEvent->Branch("mm_mc",&mm_mc,"mm_mc/D");
    TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
    TEvent->Branch("z_sig_mc",&z_mc_sig,"z_sig_mc/D");
    TEvent->Branch("z_asc_mc",&z_mc_asc,"z_asc_mc/D");
    TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
    TEvent->Branch("flv_mc",&flv_asc,"flv_asc/I");
    TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
  }

  if(type){
    Int_t bpf,d0f,pif;
    tree->SetBranchAddress("bpf",&bpf);
    tree->SetBranchAddress("d0f",&d0f);
    tree->SetBranchAddress("pif",&pif);
    tree->SetBranchAddress("nptag",&nptag);

    TEvent->Branch("bpf",&bpf,"bpf/I");
    TEvent->Branch("d0f",&d0f,"d0f/I");
    TEvent->Branch("pif",&pif,"pif/I");
    TEvent->Branch("nptag",&nptag,"nptag/I");
  }

  //Kinematic variables
  Double_t mbc,de,atckpi_max,atckpi_pi,mk,md,costhBcms;
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("atckpi_pi",&atckpi_pi);
  tree->SetBranchAddress("mk",&mk);
  tree->SetBranchAddress("md",&md);
  tree->SetBranchAddress("costhBcms",&costhBcms);

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");
  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("atckpi_pi",&atckpi_pi,"atckpi_pi/D");
  TEvent->Branch("mk",&mk,"mk/D");
  TEvent->Branch("md",&md,"md/D");
  TEvent->Branch("costhBcms",&costhBcms,"costhBcms/D");

  //Time
  Double_t dz,dz_d0;
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("dz_d0",&dz_d0);

  TEvent->Branch("dz",&dz,"dz/D");
  TEvent->Branch("dz_d0",&dz_d0,"dz_d0/D");

  int nevents = 0;
  int nrecords = 0;

  const int NTot = tree->GetEntries();
  int my_des;
  vector<double> mdv;
  vector<int> recnum;
  tree->GetEvent(0);
  mdv.push_back(TMath::Abs(md-DMass));

  recnum.push_back(0);
  int cur_evtn = evtn;
  int cur_run = run;
  for(int i=1; i<(NTot); i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);
    nrecords++;
    if(cur_evtn == evtn && cur_run == run){
      mdv.push_back(TMath::Abs(md-DMass));
      recnum.push_back(i);
    } else{
      cur_evtn = evtn; cur_run = run;
      nevents++;
      if(mdv.size() > 1){
        my_des = min_decision(mdv);
        tree->GetEvent(recnum[my_des]);
        TEvent->Fill();
        tree->GetEvent(i);
      }
      else{
        tree->GetEvent(i-1);
        TEvent->Fill();
        tree->GetEvent(i);
      }

      mdv.clear();
      recnum.clear();
      mdv.push_back(TMath::Abs(md-DMass));
      recnum.push_back(i);
    }
  }

  TEvent->Print();
  TEvent->Write();
  ofile->Close();

  cout << "Multiplicity = " << nrecords << "/" << nevents << " = " << (double)nrecords/(double)nevents << " " << endl;

  return;
}
