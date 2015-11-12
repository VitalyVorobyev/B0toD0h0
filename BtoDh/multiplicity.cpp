#include "cuts.h"

using namespace std;

int true_decision(const vector<int>& v,vector<int>& sigv,bool& true_mult_flag,bool& nosig_flag){
  true_mult_flag = false;
  nosig_flag = true;
  sigv.clear();

  for(int i=0; i<v.size(); i++){
//    if(v.size()>1) cout << v[i] << " ";
    if(v[i] == 1 || v[i] == 5 || v[i] == 10){
      nosig_flag = false;
      sigv.push_back(i);
      }
      else{
	true_mult_flag = true;
      }
  }
//  if(v.size()>1) cout << endl;
  return sigv.size();
}

int min_decision(const vector<double> v){
    if(v.size()<2){
        cout << "min_decision: size = " << v.size() << endl;
        return -1;
    }
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

void multiplicity(const int type){
  switch(type){
    case 11:
      TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigPi0_full.root");
      break;
    case 12:
      TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigEta_full.root");
      break;
    case 13:
      TFile *ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/fil_b2dh_sigOmega_full.root");
      break;
    default:
      cout << "Wrong type " << type << endl;
      return;
  }

  TTree *tree = (TTree*)ifile->Get("TEvent");
  TRandom3 rndm(0);

  int correct_dchi2 = 0;
  int correct_bchi2 = 0;
  int correct_bdchi2 = 0;
  int correct_thr = 0;
  int correct_mbc = 0;
  int correct_md = 0;
  int correct_my = 0;
  int correct_first = 0;
  int correct_random = 0;
  int correct_xi = 0;

  Double_t cos_b0,p_ks,p_pp,p_pm,p_pi0,chi2_ndf_D0,chi2_ndf_B0,chi2_tag_vtx,cos_thr,thr_sig,thr_oth;
  Int_t ndf_tag_vtx,phsp,bin,exp,run,evtn;
  Int_t b0f;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;

  Double_t mbc,de,bdtgs,mp,mm,dz,dt,atckpi_max,mh0,mpi0,mk,md;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err,dzerr;
  Double_t bdtgfr;
  Int_t mode;

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("b0f",&b0f);
  tree->SetBranchAddress("chi2_ndf_D0",&chi2_ndf_D0);
  tree->SetBranchAddress("chi2_ndf_B0",&chi2_ndf_B0);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("bdtgfr",&bdtgfr);

  tree->SetBranchAddress("md",&md);
  tree->SetBranchAddress("mk",&mk);
  tree->SetBranchAddress("mh0",&mh0);
  tree->SetBranchAddress("mpi0",&mpi0);
  tree->SetBranchAddress("mode",&mode);

  int nevents = 0;
  int nmult_events = 0;
  int nrecords = 0;
  int nback = 0;

  int nrecords_raw = 0;
  int multvector[15];
  for(int i=0; i<15; i++) multvector[i] = 0;

  const int NTot = tree->GetEntries();
  vector<double> mbcv,bchi2v,dchi2v,bdchi2v,costhrv,mdv,h0v;
  vector<int> b0fv;
  int mbc_des,bchi2_des,dchi2_des,bdchi2_des,thr_des,rndm_des,md_des;
  vector<int> sigv;
  bool true_mult_flag,nosig_flag;
  int nsig;
  int cur_evtn;
  int cur_run;
  tree->GetEvent(0);
  mbcv.push_back(mbc);
  bchi2v.push_back(chi2_ndf_B0);
  dchi2v.push_back(chi2_ndf_D0);
  bdchi2v.push_back(chi2_ndf_B0+chi2_ndf_D0);
  costhrv.push_back(TMath::Abs(cos_thr));
  b0fv.push_back(b0f);
  mdv.push_back(TMath::Abs(DMass-md));
  switch(mode){
  case 1:
    h0v.push_back(TMath::Abs(mh0-Pi0Mass));
    break;
  case 2:
    h0v.push_back(TMath::Abs(mh0-EtaMass));
    break;
  case 3:
    h0v.push_back(TMath::Abs(mh0-OmegaMass));
    break;
  }
  cout << TMath::Abs(DMass-md) << " " << md << endl;
  cur_evtn = evtn; cur_run = run;
  for(int i=1; i<(NTot); i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);
    if(md<md_min || md>md_max) continue;
    if(mpi0<(Pi0Mass-mpi0_cut) || mpi0>(Pi0Mass+mpi0_cut)) continue;
    if(mk<(KMass-mk_cut) || mk>(KMass+mk_cut)) continue;

    nrecords++;
    nrecords_raw++;
    //cout << b0f << "=";
    if(cur_evtn == evtn && cur_run == run){
      mbcv.push_back(mbc);
      bchi2v.push_back(chi2_ndf_B0);
      dchi2v.push_back(chi2_ndf_D0);
      bdchi2v.push_back(chi2_ndf_B0+chi2_ndf_D0);
      costhrv.push_back(TMath::Abs(cos_thr));
      b0fv.push_back(b0f);
      mdv.push_back(TMath::Abs(DMass-md));
      switch(mode){
      case 1:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 2:
        h0v.push_back(TMath::Abs(mh0-EtaMass));
        break;
      case 3:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      }
    }
    else{
      cur_evtn = evtn; cur_run = run;
//      cout << "___ " << b0fv[0] << " ____";
      true_decision(b0fv,sigv,true_mult_flag,nosig_flag);
//      if(b0fv.size() > 3) cout << " flags: " << true_mult_flag << " " << nosig_flag << endl;
      if(nosig_flag){
	nrecords -= b0fv.size(); nback++;
	nrecords_raw -= b0fv.size();
	mbcv.clear();
        bchi2v.clear();
        dchi2v.clear();
        bdchi2v.clear();
        costhrv.clear();
        b0fv.clear();
        mdv.clear();
        h0v.clear();

        mbcv.push_back(mbc);
        bchi2v.push_back(chi2_ndf_B0);
        dchi2v.push_back(chi2_ndf_D0);
        bdchi2v.push_back(chi2_ndf_B0+chi2_ndf_D0);
        costhrv.push_back(TMath::Abs(cos_thr));
        b0fv.push_back(b0f);
        mdv.push_back(TMath::Abs(DMass-md));
        switch(mode){
        case 1:
          h0v.push_back(TMath::Abs(mh0-Pi0Mass));
          break;
        case 2:
          h0v.push_back(TMath::Abs(mh0-EtaMass));
          break;
        case 3:
          h0v.push_back(TMath::Abs(mh0-OmegaMass));
          break;
        }
	continue;
      }
      nevents++;
      multvector[b0fv.size()-1]++;
      nrecords -= (sigv.size()-1);
      nsig = sigv.size();
      if(true_mult_flag){
	nmult_events++;
          mbc_des   = max_decision(mbcv);
          bchi2_des = min_decision(bchi2v);
          dchi2_des = min_decision(dchi2v);
          bdchi2_des= min_decision(bdchi2v);
          thr_des   = min_decision(costhrv);
          md_des    = min_decision(mdv);
          my_des    = my_decision(mdv,h0v);
          rndm_des  = (int)(mbcv.size()*rndm.Rndm());

          if(sigv[0] == 0) correct_first++;
	  if(b0fv.size()>3) cout << "Signal index: ";
          for(int j=0; j<nsig; j++){
	    if(b0fv.size()>3) cout << sigv[j] << " ";
              if(sigv[j] == mbc_des)   correct_mbc++;
              if(sigv[j] == bchi2_des) correct_bchi2++;
              if(sigv[j] == dchi2_des) correct_dchi2++;
              if(sigv[j] == bdchi2_des)correct_bdchi2++;
              if(sigv[j] == thr_des)   correct_thr++;
              if(sigv[j] == rndm_des)  correct_random++;
              if(sigv[j] == md_des)    correct_md++;
              if(sigv[j] == my_des)    correct_my++;
          }
          if(b0fv.size()>3) cout << endl;
      }
      mbcv.clear();
      bchi2v.clear();
      dchi2v.clear();
      bdchi2v.clear();
      costhrv.clear();
      b0fv.clear();
      mdv.clear();
      h0v.clear();

      mbcv.push_back(mbc);
      bchi2v.push_back(chi2_ndf_B0);
      dchi2v.push_back(chi2_ndf_D0);
      bdchi2v.push_back(chi2_ndf_B0+chi2_ndf_D0);
      costhrv.push_back(TMath::Abs(cos_thr));
      b0fv.push_back(b0f);
      mdv.push_back(TMath::Abs(DMass-md));
      switch(mode){
      case 1:
        h0v.push_back(TMath::Abs(mh0-Pi0Mass));
        break;
      case 2:
        h0v.push_back(TMath::Abs(mh0-EtaMass));
        break;
      case 3:
        h0v.push_back(TMath::Abs(mh0-OmegaMass));
        break;
      }
    }
  }

  if(nevents && nmult_events){
    cout << "First:    " << correct_first << ", " << 100.*correct_first/nmult_events << "\t"  << 100. - 100.*(nmult_events-correct_first)/nevents << endl;
    cout << "Random:   " << correct_random << ", " << 100.*correct_random/nmult_events << "\t" << 100. - 100.*(nmult_events-correct_random)/nevents << endl;
    cout << "Mbc:      " << correct_mbc << ", " << 100.*correct_mbc/nmult_events << "\t"    << 100. - 100.*(nmult_events-correct_mbc)/nevents << endl;
    cout << "mD:       " << correct_md << ", " << 100.*correct_md/nmult_events << "\t"     << 100. - 100.*(nmult_events-correct_md)/nevents << endl;
    cout << "My:       " << correct_my << ", " << 100.*correct_my/nmult_events << "\t"     << 100. - 100.*(nmult_events-correct_my)/nevents << endl;
//    cout << "B0 chi2:  " << correct_bchi2 << ", " << 100.*correct_bchi2/nmult_events << "\t"  << 100. - 100.*(nmult_events-correct_bchi2)/nevents << endl;
    cout << "D0 chi2:  " << correct_dchi2 << ", " << 100.*correct_dchi2/nmult_events << "\t"  << 100. - 100.*(nmult_events-correct_dchi2)/nevents << endl;
//    cout << "Sum chi2: " << correct_bdchi2 << ", " << 100.*correct_bdchi2/nmult_events << "\t" << 100. - 100.*(nmult_events-correct_bdchi2)/nevents << endl;
    cout << "Thrust:   " << correct_thr << ", " << 100.*correct_thr/nmult_events << "\t"    << 100. - 100.*(nmult_events-correct_thr)/nevents << endl;
    
    cout << "Multiplicity = " << nrecords << "/" << nevents << " = " << (double)nrecords/(double)nevents << endl;
    cout << "Background events: " << nback << endl;
    cout << "Raw multiplicity = " << nrecords_raw << "/" << nevents << " = " << (double)nrecords_raw/nevents << endl;
    cout << "Mult structure: ";
    for(int i=0; i<15; i++) cout << multvector[i] << " ";
    cout << endl;
  } else{
    cout << "nevents =      " << nevents << endl;
    cout << "nmult_events = " << nmult_events << endl;
  }

  return;
}
