#include "filter.h"
//#include "../rooksfw/ksfw_pi0.h"
//#include "../rooksfw/ksfw_etagg.h"
//#include "../rooksfw/ksfw_etappp.h"
//#include "../rooksfw/ksfw_omega.h"
//#include "../rooksfw/ksfw_etapgg.h"
//#include "../rooksfw/ksfw_etapppp.h"
//#include "../rooksfw/ksfw_dst0pi0.h"
//#include "../rooksfw/ksfw_dst0etagg.h"
//#include "../rooksfw/ksfw_dst0etappp.h"

#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>

using namespace std;

const bool make_substitution = true;

void switchAB(double& A, double& B){
  A += B;
  B = A - B;
  A -= B;
  return;
}

void GetShuffledVector(const int size,vector<int> &vec){
  vec.clear();
  for(int i=0; i<size; i++) vec.push_back(i);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  shuffle(vec.begin(),vec.end(),default_random_engine(seed));
  return;
}

filter::filter(const int type_in, const int type_out, const int str_num, const string& angle):
  d0_des(0),h0_des(0),totl_des(0),good_des(0),ambiguity(0),n_sig_cand(0),n_bkg_cand(0),sig_counter(0)
{
  Phsp = new PhaseSpace();
  cuts = new MyParams();
  M_SIGMC = false;
  M_DATA  = false;

  m_type_in  = type_in;
  m_type_out = type_out;
//  intree = new TChain("TEvent");
  prefix = string("/home/vitaly/B0toDh0/Tuples/");

  SetMVA();
  SetInput(str_num,angle);
  string toutstr;
  SetOutput(toutstr);
  line_out += toutstr + string(".root");

  cout << "InFile(s):" << endl;
  for(int i=0; i<line_in.size(); i++) cout << "  " << line_in[i] << endl;
  cout << "OutFile:" << endl;
  cout << "  " << line_out << endl;

//  outtree   = new TTree("rawTEvent","rawTEvent");
//  outfile      = new TFile(line_out.c_str(),"RECREATE");
//  mult_outtree = new TTree("TEvent","TEvent");

  for(int i=0; i<10; i++){
    cand_struct1[i] = 0;
    cand_struct2[i] = 0;
  }
}

void filter::MakeNTuples(void){
  Filter();
  MultiFilter();
//  outtree->Write();
//  mult_outtree->Write();
//  outfile->Close();
  return;
}

void filter::PrepareSigEvtVec(vector<ICPVEvent>& vec1,vector<ICPVEvent>& vec2,vector<int>& ind1,vector<int>& ind2){
  TChain insigchain("TEvent");
  const string fname = prefix + insinfile;
  insigchain.Add(fname.c_str());
  TTree*  insigtree_svd1 = insigchain.CopyTree("exp < 30 && (b0f == 1 || b0f == 5 || b0f == 10)");
  TTree*  insigtree_svd2 = insigchain.CopyTree("exp > 30 && (b0f == 1 || b0f == 5 || b0f == 10)");
  const int Nsvd1 = insigtree_svd1->GetEntries();
  const int Nsvd2 = insigtree_svd2->GetEntries();
  GetShuffledVector(Nsvd1,ind1);
  GetShuffledVector(Nsvd2,ind2);
  ICPVEvent::FillVectorWithTTree(vec1,insigtree_svd1,m_data_type,false);
  ICPVEvent::FillVectorWithTTree(vec2,insigtree_svd2,m_data_type,false);
  return;
}

void filter::Filter(void){
  ICPVEvent evt(m_data_type,false);
  const string rawfname = prefix + string("raw_") + line_out;
  cout << rawfname << endl;
  TFile rawfile;
  rawfile.Open(rawfname.c_str(),"RECREATE");
  TChain InTree("TEvent");
  for(int i=0; i<line_in.size(); i++){
    const string infname = prefix + line_in[i];
    cout << infname << endl;
    InTree.Add(infname.c_str());
  }

  evt.SetBrAddresses(&InTree);
  cout << "BranchAddresses set" << endl;

  TTree rawotree("rawTEvent","rawTEvent");
  evt.SetBranches(&rawotree);
  cout << "Branches set" << endl;

  vector<ICPVEvent> sigvec_svd1;
  vector<ICPVEvent> sigvec_svd2;
  vector<int> sigind_svd1;
  vector<int> sigind_svd2;
  if(make_substitution && m_type_in == 22) PrepareSigEvtVec(sigvec_svd1,sigvec_svd2,sigind_svd1,sigind_svd2);

  int svd1ind = 0;
  int svd2ind = 0;
  const int NTot = InTree.GetEntries();
  for(int i=1; i<NTot; i++){
    if(!(i%10000)){ cout << i << " events" << " " << evt.mbc << " " << evt.de << " " << evt.b0f << endl;}
    InTree.GetEvent(i);
    if(evt.mode != m_out_mode || evt.h0mode != m_out_h0mode) continue;
    if(make_substitution && m_type_in == 22 && (evt.b0f == 1 || evt.b0f == 5 || evt.b0f == 10)){
      if(evt.exp<30) evt = sigvec_svd1[sigind_svd1[svd1ind++]];//insigtree_svd1->GetEvent(sigindex_svd1[svd1ind++]);
      else           evt = sigvec_svd2[sigind_svd2[svd2ind++]];//insigtree_svd2->GetEvent(sigindex_svd2[svd2ind++]);
      cout << evt.mode << " " << evt.h0mode << " " << svd1ind << " " << svd2ind << endl;
    }

    if(!M_DATA && (evt.b0f == 0 || evt.b0f<-1)) continue;

    if(evt.mk<cuts->get_mk_min()       || evt.mk>cuts->get_mk_max())       continue;
    if(evt.de<cuts->get_de_fit_min()   || evt.de>cuts->get_de_fit_max())   continue;
    if(evt.mbc<cuts->get_mbc_fit_min() || evt.mbc>cuts->get_mbc_fit_max()) continue;
    if(evt.md_raw<cuts->get_md_min()   || evt.md_raw>cuts->get_md_max())   continue;
    evt.dmdst0 = evt.mdst0 - evt.md_raw;
    if(evt.mode>9 && (evt.dmdst0<cuts->get_dm_dst0_min() || evt.dmdst0>cuts->get_dm_dst0_max())) continue;
    if(evt.mode == 5 && (evt.dmetap<cuts->get_dm_etap_min(evt.h0mode) || evt.dmetap>cuts->get_dm_etap_max(evt.h0mode))) continue;
    if(evt.mh0<cuts->get_mh0_min(evt.mode,evt.h0mode) || evt.mh0>cuts->get_mh0_max(evt.mode,evt.h0mode)) continue;

    if( abs(evt.r_pip)>2 || abs(evt.r_pim)>2) continue;
    if( abs(evt.z_pip)>5 || abs(evt.z_pim)>2) continue;
    if((abs(evt.r_pi1)>2 || abs(evt.r_pi2)>2) && evt.r_pi1 != -99) continue;
    if((abs(evt.z_pi1)>5 || abs(evt.z_pi2)>5) && evt.z_pi1 != -99) continue;

    evt.flv = evt.tag_LH>0 ? -1 : 1;
    if(evt.e_g1 > evt.e_g2){
      switchAB(evt.e_g1,evt.e_g2);
      switchAB(evt.th_g1,evt.th_g2);
    }

    if(evt.e_g3 > evt.e_g4){
      switchAB(evt.e_g3,evt.e_g4);
      switchAB(evt.th_g3,evt.th_g4);
    }

    if(!M_DATA){
      evt.d0ch0 = evt.d0_chain[0];
      evt.d0ch1 = evt.d0_chain[1];
      evt.d0ch2 = evt.d0_chain[2];
      evt.d0ch3 = evt.d0_chain[3];

      evt.h0ch0 = evt.h0_chain[0];
      evt.h0ch1 = evt.h0_chain[1];
      if(!(evt.b0f == 1 || evt.b0f == 5 || evt.b0f == 10) && (evt.d0f == 1 || evt.d0f == 10) && (evt.h0f == 1 || evt.h0f == 3) && evt.b0f != 0 && evt.mbc>5.2 && abs(evt.d0ch1) == 423 && abs(evt.d0ch3) == 511 && (evt.mode == 10 || evt.mode == 20) && abs(evt.h0ch1) == 511 && (evt.h0ch0 == 111 || evt.h0ch0 == 221) && abs(evt.d0ch0) == 421){
        evt.rndm_pi0 = 1;
      } else{
        evt.rndm_pi0 = 0;
      }
    }

    if(!M_DATA){
      evt.bin_mc = Phsp->GetBin(evt.mp_mc,evt.mm_mc);
      if((evt.mp_mc>evt.mm_mc && evt.flv_mc*evt.bin_mc>0) || (evt.mp_mc<evt.mm_mc && evt.flv_mc*evt.bin_mc<0)) evt.bin_mc = -evt.bin_mc;
    }
    evt.bin = Phsp->GetBin(evt.mp,evt.mm);
    if((evt.mp>evt.mm && evt.tag_LH*evt.bin<0) || (evt.mp<evt.mm && evt.tag_LH*evt.bin>0)) evt.bin = -evt.bin;
    if(evt.bin) evt.phsp = 1;
    else        evt.phsp = 0;

// **     Bin sign conventions     ** //
//    bin*tag_LH>0 if mp>mm           //
//    bin_mc*flv_mc<0 if mp_mc>mm_mc  //
// ////////////////////////////////// //

    evt.z_sig *= 10; evt.z_asc *=10;// z_sig_d0 *= 10;
    evt.dz = evt.z_sig - evt.z_asc;
    if(M_SIGMC){
      evt.dz_mc_sig = evt.z_sig-evt.z_sig_mc;
      evt.dz_mc_asc = evt.z_asc-evt.z_asc_mc;
      evt.dt_mc     = evt.t_sig_mc - evt.t_asc_mc;
      evt.dz_mc     = evt.z_sig_mc - evt.z_asc_mc;
    }
    evt.sz_sig = 10.*TMath::Sqrt(evt.sz_sig);
    evt.sz_asc = 10.*TMath::Sqrt(evt.sz_asc);

    // * Standatd ICPV cuts * //
    evt.good_icpv = IsGoodICPV(evt.ndf_z_sig,evt.sz_sig,evt.chisq_z_sig,evt.ndf_z_asc,evt.sz_asc,evt.chisq_z_asc);
    if(!evt.good_icpv) continue;
    // * ////////////////// * //

    evt.lh1 = lh(evt.mode,evt.h0mode,evt.k1mm2,evt.k1vars,true);
    evt.lh0 = lh(evt.mode,evt.h0mode,evt.k0mm2,evt.k0vars,false);

    tmvaevt.Fill(evt);
//    m_lh0          = (float)lh0;//mode < 3 ? (float)lh1 : (float)(1./1.005-lh1*lh1);
//    m_costhBcms    = (float)abs(costhBcms);
//    m_chi2_mass_d0 = (float)log(chi2_mass_d0);
//    m_cos_thr      = (float)abs(cos_thr);
//    m_thr_sig      = (float)thr_sig;
//    m_h0_chi2      = (float)h0_chi2;
//    m_egamma       = (float)log(egamma);
//    m_p_pi0_h0     = (float)p_pi0_h0;
//    m_cos_hel      = (float)abs(cos_hel);

    evt.bdt = BDT(evt.mode,evt.h0mode);

    if(!M_DATA && !M_SIGMC){
      evt.bin_mc = evt.bin;
      evt.flv_mc = evt.tag_LH > 0 ? -1 : 1;
    }
//    outtree->Fill();
    rawotree.Fill();
  }
  cout << "After the first stage: " << rawotree.GetEntries() << " events" << endl;
  rawotree.Write();
  rawfile.Close();
  return;
}

void filter::MultiFilter(void){
  TChain intr("rawTEvent","rawTEvent");
  const string rawfname = prefix + string("raw_") + line_out;
  intr.Add(rawfname.c_str());
  cout << "... And still " << intr.GetEntries() << " events" << endl;
//  SetBranchAddresses(outtree,true);
//  SetBranches(mult_outtree);

  ICPVEvent evt(m_data_type,true);
  evt.SetBrAddresses(&intr);

  const string ofname = prefix + line_out;
  TFile ofile(ofname.c_str(),"RECREATE");
  TTree outr("TEvent","TEvent");
  evt.SetBranches(&outr);

  int nevents = 0;
  int nrecords = 0;

//  const int NTot = outtree->GetEntries();
  const int NTot =  intr.GetEntries();
  cout << "... and " << NTot << endl;
  int my_des;
  vector<double> mdv;
  vector<double> h0v;
  vector<double> dmv;
  vector<int> bflags;
  vector<int> modev;
  vector<int> h0modev;
  vector<int> recnum;
  bool dst0_flag = false;
  bool signal_flag = false;
//  outtree->GetEvent(0);
  intr.GetEvent(0);
  mdv.push_back(evt.md_raw);
  h0v.push_back(evt.mh0);
  modev.push_back(evt.mode); h0modev.push_back(evt.h0mode);
  if(evt.mode == 5)    dmv.push_back(evt.dmetap);
  else if(evt.mode>9){ dmv.push_back(evt.dmdst0); dst0_flag = true;}
  else                 dmv.push_back(0);
  recnum.push_back(0);
  if(!M_DATA){
    bflags.push_back(evt.b0f);
    if(evt.b0f == 1 || evt.b0f == 5 || evt.b0f == 10) signal_flag = true;
  } else bflags.push_back(0);
  int cur_evtn = evt.evtn;
  int cur_run  = evt.run;
  for(int i=1; i<NTot; i++){
    if(!(i%10000)){ cout << i << " events" << ", mbc = " << evt.mbc << endl;}
//    outtree->GetEvent(i);
    intr.GetEvent(i);
    nrecords++;
    if(cur_evtn == evt.evtn && cur_run == evt.run){
      mdv.push_back(evt.md_raw); h0v.push_back(evt.mh0);
      if(evt.mode == 5)    dmv.push_back(evt.dmetap);
      else if(evt.mode>9){ dmv.push_back(evt.dmdst0); dst0_flag = true;}
      else                 dmv.push_back(0);
      modev.push_back(evt.mode); h0modev.push_back(evt.h0mode);
      recnum.push_back(i);
      if(!M_DATA){
        bflags.push_back(evt.b0f);
        if(evt.b0f == 1 || evt.b0f == 5 || evt.b0f == 10) signal_flag = true;
      } else bflags.push_back(0);
    } else{
      cur_evtn = evt.evtn; cur_run = evt.run;
      nevents++;
      if(mdv.size() > 1){
        my_des = my_decision2(mdv,h0v,dmv,bflags,modev,h0modev);
//        outtree->GetEvent(recnum[my_des]);
        intr.GetEvent(recnum[my_des]);
      } else{
//        outtree->GetEvent(i-1);
        intr.GetEvent(i-1);
      }
//      mult_outtree->Fill();
//      outtree->GetEvent(i);
      outr.Fill();
      intr.GetEvent(i);
      if(signal_flag) sig_counter++;
      mdv.clear(); h0v.clear(); modev.clear(); h0modev.clear(); dmv.clear();
      recnum.clear(); bflags.clear(); dst0_flag = false; signal_flag = false;
      mdv.push_back(evt.md_raw); h0v.push_back(evt.mh0);
      if(evt.mode == 5)    dmv.push_back(evt.dmetap);
      else if(evt.mode>9){ dmv.push_back(evt.dmdst0); dst0_flag = true;}
      else                 dmv.push_back(0);
      modev.push_back(evt.mode); h0modev.push_back(evt.h0mode);
      recnum.push_back(i);
      if(!M_DATA){
        bflags.push_back(evt.b0f);
        if(evt.b0f == 1 || evt.b0f == 5 || evt.b0f == 10) signal_flag = true;
      } else bflags.push_back(0);
    }
  }

  cout << "Multiplicity = " << nrecords << "/" << nevents << " = " << (double)nrecords/(double)nevents << " " << endl;
  cout << "  D0 decisions: " << d0_des << endl;
  cout << "  h0 decisions: " << h0_des << endl;
  if(!M_DATA){
    if(totl_des){
      cout << "True mult: = " << (n_bkg_cand+n_sig_cand) << "/" << totl_des << " = " << (double)(n_bkg_cand+n_sig_cand)/totl_des << endl;
      cout << "Good decisions: " << good_des << "/" << totl_des << " = " << (double)good_des/totl_des << endl;
    }
    cout << "Multiplicity structure:" << endl;
    double aver_mult = 0;
    int sum_events = 0;
    for(int i=0; i<10; i++){
      aver_mult  += (i+2)*cand_struct2[i];
      sum_events += cand_struct2[i];
      cout << "  " << i+1 << ": " << cand_struct1[i] << " (" << cand_struct2[i] << ")" << endl;
    }
    cout << "Average multiplicity: " << aver_mult/sum_events << " with " << sum_events << " candidates and " << sig_counter << " signal events." << endl;
    const double mult = (double)(sig_counter + aver_mult - sum_events)/(double)(sig_counter);
    cout << "Multiplicity: " << mult << endl;
    cout << "Signal loss:  " << 100.*(totl_des-good_des)/sig_counter << endl;
  }

  outr.Write();
  ofile.Close();
  return;
}

void filter::SetInput(const int str_num, const string &angle){
  stringstream out;
  line_in.clear();
  switch(m_type_in){
  case 0://Data
    line_in.push_back(string("b2dh_data.root"));
    line_out = string("Fil_b2dh_data");//.root");
    M_DATA = true;
    break;
  case 11://Signal MC pi0
    if(angle != string("")){
      out.str("");
      out << string("Linearity/b2dh_sigmcPI0") << angle << string("_s1.root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Linearity/Fil_b2dh_sigmcPi0") << angle << string("_s1.root");
      line_out = out.str();
    } else if(str_num){
      out.str("");
      out << string("b2dh_sigmcPI0_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcPi0_s") << str_num;// << "_" << toutstr << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcPI0_s8.root"));
      line_out = string("Fil_b2dh_sigmcPi0_s8");// + toutstr + string(".root");
    }
    M_SIGMC = true;
    break;
  case 101://Signal MC D*0h0
    if(str_num){
      out.str("");
      out << string("b2dh_sigmcDST0_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcDST0_s") << str_num;// << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcDST0_s2.root"));
      line_out = string("Fil_b2dh_sigmcDST0_s2");//.root");
    }
    M_SIGMC = true;
    break;
  case 12://Signal MC eta
    if(str_num){
      out.str("");
      out << string("b2dh_sigmcETA_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcETA_s") << str_num;// << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcETA_s3.root"));
      line_out = string("Fil_b2dh_sigmcETA_s3");//.root");
    }
    M_SIGMC = true;
    break;
  case 13://Signal MC omega
    if(angle != string("")){
      out.str("");
      out << string("Linearity/b2dh_sigmcOMEGA") << angle << string("_s1.root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Linearity/Fil_b2dh_sigmcOMEGA") << angle << string("_s1.root");
      line_out = out.str();
    } else if(str_num){
      out.str("");
      out << string("b2dh_sigmcOMEGA_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcOMEGA_s") << str_num;// << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcOMEGA_s6.root"));
      line_out = string("Fil_b2dh_sigmcOMEGA_s6");//.root");
    }
    M_SIGMC = true;
    break;
  case 14://Signal MC rho
    if(str_num){
      out.str("");
      out << string("b2dh_sigmcRHO_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcRHO_s") << str_num;// << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcRHO_s1.root"));
      line_out = string("Fil_b2dh_sigmcRHO_s1");//.root");
    }
    M_SIGMC = true;
    break;
  case 15://Signal MC eta'
    if(str_num){
      out.str("");
      out << string("b2dh_sigmcETAP_s") << str_num << string(".root");
      line_in.push_back(out.str());
      out.str("");
      out << string("Fil_b2dh_sigmcETAP_s") << str_num;// << string(".root");
      line_out = out.str();
    } else{
      line_in.push_back(string("b2dh_sigmcETAP_s1.root"));
      line_out = string("Fil_b2dh_sigmcETAP_s1");//.root");
    }
    M_SIGMC = true;
    break;
  case 21://charged
    out.str("");
    out << string("GenMC/s") << str_num << string("/b2dh_charged_") << str_num+10 << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_charged_") << str_num << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("Fil_b2dh_charged_") << str_num << "_" << str_num+10;// << string(".root");
    line_out = out.str();
    break;
  case 22://mixed
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_mixed_") << str_num+10 << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_mixed_") << str_num << string(".root");
    line_in.push_back(out.str());
    out.str("");
    out << string("Fil_b2dh_mixed_") << str_num << "_" << str_num+10;// << string(".root");
    if(make_substitution) out << "_subs";
    line_out = out.str();
    break;
  case 23://charm
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_charm_") << str_num+10 << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_charm_") << str_num << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("Fil_b2dh_charm_") << str_num << "_" << str_num+10;// << string(".root");
    line_out = out.str();
    break;
  case 24://uds
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_uds_") << str_num+10 << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("GenMC/s") << str_num << string("/b2dh_uds_") << str_num << string(".root");
    line_in.push_back(out.str());
    out.str(""); out << string("Fil_b2dh_uds_") << str_num << "_" << str_num+10;// << string(".root");
    line_out = out.str();
    break;
  default:
    break;
  }
  line_out += "_v2";
  return;
}

void filter::SetOutput(string& toutstr){
  switch(m_type_out){
  case 1:
    toutstr = string("_m1_h0m10");
    m_out_mode = 1;
    m_out_h0mode = 10;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcPI0_s8.root");
    break;
  case 2:
    toutstr = string("_m2_h0m10");
    m_out_mode = 2;
    m_out_h0mode = 10;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcETA_s3.root");
    break;
  case 3:
    toutstr = string("_m2_h0m20");
    m_out_mode = 2;
    m_out_h0mode = 20;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcETA_s3.root");
    break;
  case 4:
    toutstr = string("_m3_h0m20");
    m_out_mode = 3;
    m_out_h0mode = 20;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcOMEGA_s6.root");
    break;
  case 5:
    toutstr = string("_m5_h0m10");
    m_out_mode = 5;
    m_out_h0mode = 10;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcETAP_s1.root");
    break;
  case 10:
    toutstr = string("_m10_h0m10");
    m_out_mode = 10;
    m_out_h0mode = 10;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcDST0_s2.root");
    break;
  case 20:
    toutstr = string("_m20_h0m10");
    m_out_mode = 20;
    m_out_h0mode = 10;
    if(make_substitution && m_type_in == 22) insinfile = string("b2dh_sigmcDST0_s2.root");
    break;
  default:
    cout << "Wrong out mode " << m_type_out << endl;
    return;
  }
}

void filter::SetMVA(void){
  ksfw1_pi0        = new rooksfw("ksfw1_pi0",   ksfw1_alpha_pi0,   ksfw1_sigpdf_pi0,   ksfw1_bkgpdf_pi0);
  ksfw1_etagg      = new rooksfw("ksfw1_etagg", ksfw1_alpha_etagg, ksfw1_sigpdf_etagg, ksfw1_bkgpdf_etagg);
  ksfw1_etappp     = new rooksfw("ksfw1_etappp",ksfw1_alpha_etappp,ksfw1_sigpdf_etappp,ksfw1_bkgpdf_etappp);
  ksfw1_omega      = new rooksfw("ksfw1_omega", ksfw1_alpha_omega, ksfw1_sigpdf_omega, ksfw1_bkgpdf_omega);

  ksfw0_pi0        = new rooksfw("ksfw0_pi0",   ksfw0_alpha_pi0,   ksfw0_sigpdf_pi0,   ksfw0_bkgpdf_pi0);
  ksfw0_etagg      = new rooksfw("ksfw0_etagg", ksfw0_alpha_etagg, ksfw0_sigpdf_etagg, ksfw0_bkgpdf_etagg);
  ksfw0_etappp     = new rooksfw("ksfw0_etappp",ksfw0_alpha_etappp,ksfw0_sigpdf_etappp,ksfw0_bkgpdf_etappp);
  ksfw0_omega      = new rooksfw("ksfw0_omega", ksfw0_alpha_omega, ksfw0_sigpdf_omega, ksfw0_bkgpdf_omega);

  ksfw1_dst0pi0    = new rooksfw("ksfw1_dst0pi0",   ksfw1_alpha_dst0pi0,   ksfw1_sigpdf_dst0pi0,   ksfw1_bkgpdf_dst0pi0);
  ksfw1_dst0etagg  = new rooksfw("ksfw1_dst0etagg", ksfw1_alpha_dst0etagg, ksfw1_sigpdf_dst0etagg, ksfw1_bkgpdf_dst0etagg);
  ksfw1_dst0etappp = new rooksfw("ksfw1_dst0etappp",ksfw1_alpha_dst0etappp,ksfw1_sigpdf_dst0etappp,ksfw1_bkgpdf_dst0etappp);
  ksfw1_etapgg     = new rooksfw("ksfw1_etapgg",    ksfw1_alpha_etapgg,    ksfw1_sigpdf_etapgg,    ksfw1_bkgpdf_etapgg);
  ksfw1_etapppp    = new rooksfw("ksfw1_etapppp",   ksfw1_alpha_etapppp,   ksfw1_sigpdf_etapppp,   ksfw1_bkgpdf_etapppp);

  ksfw0_dst0pi0    = new rooksfw("ksfw0_dst0pi0",   ksfw0_alpha_dst0pi0,   ksfw0_sigpdf_dst0pi0,   ksfw0_bkgpdf_dst0pi0);
  ksfw0_dst0etagg  = new rooksfw("ksfw0_dst0etagg", ksfw0_alpha_dst0etagg, ksfw0_sigpdf_dst0etagg, ksfw0_bkgpdf_dst0etagg);
  ksfw0_dst0etappp = new rooksfw("ksfw0_dst0etappp",ksfw0_alpha_dst0etappp,ksfw0_sigpdf_dst0etappp,ksfw0_bkgpdf_dst0etappp);
  ksfw0_etapgg     = new rooksfw("ksfw0_etapgg",    ksfw0_alpha_etapgg,    ksfw0_sigpdf_etapgg,    ksfw0_bkgpdf_etapgg);
  ksfw0_etapppp    = new rooksfw("ksfw0_etapppp",   ksfw0_alpha_etapppp,   ksfw0_sigpdf_etapppp,   ksfw0_bkgpdf_etapppp);

  if(m_type_out == 1 || m_type_out == 10){
    reader_pi0 = new TMVA::Reader("!Color:!Silent:V");
    reader_pi0->AddVariable("abs(costhBcms)",&tmvaevt.m_costhBcms);
    reader_pi0->AddVariable("log(chi2_mass_d0)",&tmvaevt.m_chi2_mass_d0);
    reader_pi0->AddVariable("abs(cos_thr)",&tmvaevt.m_cos_thr);
    reader_pi0->AddVariable("thr_sig",&tmvaevt.m_thr_sig);
    reader_pi0->AddVariable("log(egamma)",&tmvaevt.m_egamma);
    reader_pi0->AddVariable("lh0",&tmvaevt.m_lh0);
    reader_pi0->BookMVA("pi0","/home/vitaly/B0toDh0/TMVA/weights/MVA_softcut_ksfw_Dpi0_BDT.weights.xml");
  }

  if(m_type_out == 2 || m_type_out == 20 || m_type_out == 5){
    reader_gg  = new TMVA::Reader("!Color:!Silent:V");
    reader_gg->AddVariable("abs(costhBcms)",&tmvaevt.m_costhBcms);
    reader_gg->AddVariable("log(chi2_mass_d0)",&tmvaevt.m_chi2_mass_d0);
    reader_gg->AddVariable("abs(cos_thr)",&tmvaevt.m_cos_thr);
    reader_gg->AddVariable("thr_sig",&tmvaevt.m_thr_sig);
    reader_gg->AddVariable("log(h0_chi2)",&tmvaevt.m_h0_chi2);
    reader_gg->AddVariable("log(egamma)",&tmvaevt.m_egamma);
    reader_gg->AddVariable("lh0",&tmvaevt.m_lh0);
    reader_gg->BookMVA("gg","/home/vitaly/B0toDh0/TMVA/weights/MVA_softcut_ksfw_Detagg_BDT.weights.xml");
  }

  if(m_type_out == 3){
    reader_ppp = new TMVA::Reader("!Color:!Silent:V");
    reader_ppp->AddVariable("abs(costhBcms)",&tmvaevt.m_costhBcms);
    reader_ppp->AddVariable("log(chi2_mass_d0)",&tmvaevt.m_chi2_mass_d0);
    reader_ppp->AddVariable("abs(cos_thr)",&tmvaevt.m_cos_thr);
    reader_ppp->AddVariable("thr_sig",&tmvaevt.m_thr_sig);
    reader_ppp->AddVariable("log(egamma)",&tmvaevt.m_egamma);
    reader_ppp->AddVariable("p_pi0_h0",&tmvaevt.m_p_pi0_h0);
    reader_ppp->AddVariable("lh0",&tmvaevt.m_lh0);
    reader_ppp->BookMVA("etappp","/home/vitaly/B0toDh0/TMVA/weights/MVA_softcut_ksfw_Detappp_BDT.weights.xml");
  }

  if(m_type_out == 4){
    reader_omega = new TMVA::Reader("!Color:!Silent:V");
    reader_omega->AddVariable("abs(costhBcms)",&tmvaevt.m_costhBcms);
    reader_omega->AddVariable("log(chi2_mass_d0)",&tmvaevt.m_chi2_mass_d0);
    reader_omega->AddVariable("abs(cos_thr)",&tmvaevt.m_cos_thr);
    reader_omega->AddVariable("thr_sig",&tmvaevt.m_thr_sig);
    reader_omega->AddVariable("log(egamma)",&tmvaevt.m_egamma);
    reader_omega->AddVariable("p_pi0_h0",&tmvaevt.m_p_pi0_h0);
    reader_omega->AddVariable("abs(cos_hel)",&tmvaevt.m_cos_hel);
    reader_omega->AddVariable("lh0",&tmvaevt.m_lh0);
    reader_omega->BookMVA("omega","/home/vitaly/B0toDh0/TMVA/weights/MVA_softcut_ksfw_Domega_BDT.weights.xml");
  }
  return;
}

int filter::IsGoodICPV(const int ndf_z_sig, const double& sz_sig, const double& chisq_z_sig,const int ndf_z_asc, const double& sz_asc, const double& chisq_z_asc){
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

void filter::PrintVariants(const vector<double>& d0mass,const vector<double>& h0mass, const vector<int>& b0fvec, const vector<int>& modev){
  const int n = d0mass.size();
  cout << "PrintVariants. " << n << " candidates." << endl;
  for(int i=0; i<n; i++){
    cout << i+1 << ". " << modev[i] << " " << b0fvec[i] << " " << d0mass[i] << " " << h0mass[i] << endl;
  }
  return;
}

void filter::PrintVariants2(const vector<double>& d0mass,const vector<double>& h0mass, const vector<double>& dmvec, const vector<int>& b0fvec, const vector<int>& modev, const vector<int>& h0modev){
  const int n = d0mass.size();
  cout << "PrintVariants2. " << n << " candidates." << endl;
  for(int i=0; i<n; i++){
    cout << i+1 << ". " << modev[i] << " " << b0fvec[i] << " " << d0mass[i] << " " << h0mass[i] << " " << dmvec[i];
    cout << " " << chisq_cand(modev[i],h0modev[i],d0mass[i],h0mass[i],dmvec[i]) << endl;
  }
  return;
}

double filter::chisq_cand(const int mode, const int h0mode, const double& md0, const double& mh0, const double& dm){
  double dm0 = 0;
  double dmssq = 1;
  const double dmd0sq = (cuts->get_md()-md0)*(cuts->get_md()-md0);
  const double smd0sq = (cuts->get_md_max() - cuts->get_md_min())*(cuts->get_md_max() - cuts->get_md_min())/36;
//  double chisq_md0 = 12.*(cuts->get_md()-md0)/(cuts->get_md_max() - cuts->get_md_min()); chisq_md0 *= chisq_md0;
  const double dmh0sq = (cuts->get_mh0(mode) - mh0)*(cuts->get_mh0(mode) - mh0);
  const double smh0sq = (cuts->get_mh0_max(mode,h0mode) - cuts->get_mh0_min(mode,h0mode))*(cuts->get_mh0_max(mode,h0mode) - cuts->get_mh0_min(mode,h0mode))/36;
//  double chisq_mh0 = 6.*(cuts->get_mh0(mode) - mh0)/(cuts->get_mh0_max(mode,h0mode) - cuts->get_mh0_min(mode,h0mode)); chisq_mh0 *= chisq_mh0;
  if(mode == 5){
    dm0 = 0.410; dmssq =(cuts->get_dm_etagg_max()-cuts->get_dm_etagg_min())*(cuts->get_dm_etagg_max()-cuts->get_dm_etagg_min())/36;
//    double chisq_dm = 3.*( - dm)/(cuts->get_dm_etagg_max()-cuts->get_dm_etagg_min()); chisq_dm *= chisq_dm;
//    chisq = (chisq_md0+chisq_mh0+chisq_dm)/3.5;
  } else if(mode == 10 || mode == 20){
    dm0 = 0.142; dmssq = (cuts->get_dm_dst0_max()-cuts->get_dm_dst0_min())*(cuts->get_dm_dst0_max()-cuts->get_dm_dst0_min())/36;
//    double chisq_dm = 3.*(0.142 - dm)/(cuts->get_dm_dst0_max()-cuts->get_dm_dst0_min()); chisq_dm *= chisq_dm;
//    chisq = (chisq_md0+chisq_mh0+chisq_dm)/3.5;
  }
  const double ddmsq = dm0>0 ? (dm0-dm)*(dm0-dm) : 1;
//  cout << sqrt(smd0sq) << " " << sqrt(smh0sq) << " " << sqrt(dmssq) << endl;
//  cout << " " << sqrt(dmd0sq) << " " << sqrt(dmh0sq) << " " << sqrt(ddmsq) << endl;
//  cout << " D0 " << cuts->get_md() << " -> " << md0 << endl;
//  cout << " h0 " << cuts->get_mh0(mode) << " -> " << mh0 << endl;
//  cout << " dm " << dm0 << " -> " << dm << endl;
  const double chisq = dmd0sq/smd0sq + dmh0sq/smh0sq + ddmsq/dmssq;
//  cout << "Chi2 = " << chisq << endl;
//  const double chisq = ddmsq/dmssq;
  return chisq;
}

int filter::my_decision2(const vector<double>& d0mass,const vector<double>& h0mass, const vector<double>& dmass, const vector<int>& b0fvec, const vector<int>& modev, const vector<int>& h0modev){
//  cout << "mD0[0] = " << d0mass[0] << ", mh0[0] = " << h0mass[0] << ", dm[0] = " << dmass[0] << endl;
  if(d0mass.size() != h0mass.size()){
    cout << "My decision2: wrong sizes " << d0mass.size() << ", " << h0mass.size() << endl;
    return 0;
  }
  if(d0mass.size()<2){
    cout << "My decision2: size = " << d0mass.size() << endl;
    return -1;
  }

  bool printflag = true;
  vector<double> chisq_mins;
  vector<int> indexes;
  int imin = 0;
//  cout << "Cand flag[0]: " << b0fvec[0] << ", mode " << modev[0] << ", h0mode " << h0modev[0] << endl;
  double chisq_min = chisq_cand(modev[0],h0modev[0],d0mass[0],h0mass[0],dmass[0]);
  chisq_mins.clear(); indexes.clear();
  chisq_mins.push_back(chisq_min);
  indexes.push_back(imin);
  if(modev[0]==5 && (b0fvec[0] == 1 || b0fvec[0] == 5 || b0fvec[0] == 10) && printflag){
    PrintVariants2(d0mass,h0mass,dmass,b0fvec,modev,h0modev);
    printflag = false;
  }
  for(int i=1; i<d0mass.size(); i++){
    if(modev[i]==5 && (b0fvec[i] == 1 || b0fvec[i] == 5 || b0fvec[i] == 10) && printflag){
      PrintVariants2(d0mass,h0mass,dmass,b0fvec,modev,h0modev);
      printflag = false;
    }
//    cout << "Cand flag[i]: " << b0fvec[i] << ", mode " << modev[i] << ", h0mode " << h0modev[i] << endl;
    double chisq_cur = chisq_cand(modev[i],h0modev[i],d0mass[i],h0mass[i],dmass[i]);
    if((chisq_cur+0.0001)<chisq_min){
      chisq_min = chisq_cur;
      imin = i;
      chisq_mins.clear();
      indexes.clear();
      chisq_mins.push_back(chisq_min);
      indexes.push_back(i);
    } else if(abs(chisq_cur - chisq_min) < 0.0001){
      chisq_mins.push_back(chisq_cur);
      indexes.push_back(i);
    }
  }

  if(chisq_mins.size()>1) ambiguity++;
  if(is_decision(b0fvec)){
    const int flag = b0fvec[indexes[0]];
    if(flag == 1 || flag == 5 || flag == 10) good_des++;
  }

  return indexes[0];
}

int filter::my_decision(const vector<double>& d0mass,const vector<double>& h0mass, const vector<int>& b0fvec, const vector<int>& modev){
  if(d0mass.size() != h0mass.size()){
    cout << "My decision: wrong sizes " << d0mass.size() << ", " << h0mass.size() << endl;
    return 0;
  }
  if(d0mass.size()<2){
    cout << "My decision: size = " << d0mass.size() << endl;
    return -1;
  }

  bool printflag = true;
  vector<double> mdmins,h0mins;
  vector<int> indexes;
  double mdmin = TMath::Abs(d0mass[0]-cuts->get_md());
  int imin = 0;
  mdmins.clear(); h0mins.clear(); indexes.clear();
  mdmins.push_back(TMath::Abs(d0mass[0]-cuts->get_md()));
  h0mins.push_back(TMath::Abs(h0mass[0]-cuts->h0mass(modev[0])));
  indexes.push_back(0);
  if(modev[0]>9 && (b0fvec[0] == 1 || b0fvec[0] == 5 || b0fvec[0] == 10) && printflag){
    PrintVariants(d0mass,h0mass,b0fvec,modev);
    printflag = false;
  }
  for(int i=1; i<d0mass.size(); i++){
    if(modev[i]>9 && (b0fvec[i] == 1 || b0fvec[i] == 5 || b0fvec[i] == 10) && printflag){
      PrintVariants(d0mass,h0mass,b0fvec,modev);
      printflag = false;
    }
    if(TMath::Abs(d0mass[i]-cuts->get_md())<mdmin){
      mdmin = TMath::Abs(d0mass[i]-cuts->get_md());
      imin = i;
      mdmins.clear(); h0mins.clear(); indexes.clear();
      mdmins.push_back(TMath::Abs(d0mass[i]-cuts->get_md()));
      h0mins.push_back(TMath::Abs(h0mass[i]-cuts->h0mass(modev[i])));
      indexes.push_back(i);
    } else if(abs(TMath::Abs(d0mass[i]-cuts->get_md()) - mdmin) < 0.0001){
      mdmins.push_back(TMath::Abs(d0mass[i]-cuts->get_md()));
      h0mins.push_back(TMath::Abs(h0mass[i]-cuts->h0mass(modev[i])));
      indexes.push_back(i);
    }
  }
  if(h0mins.size() == 1){
    d0_des++;
    if(b0fvec.size()){
      if(is_decision(b0fvec)){
        const int flag = b0fvec[indexes[0]];
        if(flag == 1 || flag == 5 || flag == 10) good_des++;
      }
    }
    return indexes[0];
  } else{
    h0_des++;
  }
  double h0min = h0mins[0];
  int iimin = 0;
  for(int i=1; i<h0mins.size(); i++){
    if(h0mins[i]<h0min){
      h0min = h0mins[i];
      iimin = i;
    }
  }
  if(b0fvec.size()){
    if(is_decision(b0fvec)){
      const int flag = b0fvec[indexes[iimin]];
      if(flag == 1 || flag == 5 || flag == 10) good_des++;
    }
  }
  return indexes[iimin];
}

bool filter::is_decision(const vector<int>& b0fvec){
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

double filter::lh(const int mode, const int h0mode, const double& kmm2,const double* kvars,const bool fsf){
  switch(mode){
  case 1:
    if(fsf){
      ksfw1_pi0->input(kmm2, kvars);
      return ksfw1_pi0->ls()/(ksfw1_pi0->ls()+ksfw1_pi0->lb());
    } else{
      ksfw0_pi0->input(kmm2, kvars);
      return ksfw0_pi0->ls()/(ksfw0_pi0->ls()+ksfw0_pi0->lb());
    }
  case 10:
    if(fsf){
      ksfw1_dst0pi0->input(kmm2, kvars);
      return ksfw1_dst0pi0->ls()/(ksfw1_dst0pi0->ls()+ksfw1_dst0pi0->lb());
    } else{
      ksfw0_dst0pi0->input(kmm2, kvars);
      return ksfw0_dst0pi0->ls()/(ksfw0_dst0pi0->ls()+ksfw0_dst0pi0->lb());
    }
  case 2:
    if(h0mode == 10){
      if(fsf){
        ksfw1_etagg->input(kmm2, kvars);
        return ksfw1_etagg->ls()/(ksfw1_etagg->ls()+ksfw1_etagg->lb());
      } else{
        ksfw0_etagg->input(kmm2, kvars);
        return ksfw0_etagg->ls()/(ksfw0_etagg->ls()+ksfw0_etagg->lb());
      }
    } else{
      if(fsf){
        ksfw1_etappp->input(kmm2, kvars);
        return ksfw1_etappp->ls()/(ksfw1_etappp->ls()+ksfw1_etappp->lb());
      } else{
        ksfw0_etappp->input(kmm2, kvars);
        return ksfw0_etappp->ls()/(ksfw0_etappp->ls()+ksfw0_etappp->lb());
      }
    }
  case 20:
    if(h0mode == 10){
      if(fsf){
        ksfw1_dst0etagg->input(kmm2, kvars);
        return ksfw1_dst0etagg->ls()/(ksfw1_dst0etagg->ls()+ksfw1_dst0etagg->lb());
      } else{
        ksfw0_dst0etagg->input(kmm2, kvars);
        return ksfw0_dst0etagg->ls()/(ksfw0_dst0etagg->ls()+ksfw0_dst0etagg->lb());
      }
    } else{
      if(fsf){
        ksfw1_dst0etappp->input(kmm2, kvars);
        return ksfw1_dst0etappp->ls()/(ksfw1_dst0etappp->ls()+ksfw1_dst0etappp->lb());
      } else{
        ksfw0_dst0etappp->input(kmm2, kvars);
        return ksfw0_dst0etappp->ls()/(ksfw0_dst0etappp->ls()+ksfw0_dst0etappp->lb());
      }
    }
  case 3:
    if(fsf){
      ksfw1_omega->input(kmm2, kvars);
      return ksfw1_omega->ls()/(ksfw1_omega->ls()+ksfw1_omega->lb());
    } else{
      ksfw0_omega->input(kmm2, kvars);
      return ksfw0_omega->ls()/(ksfw0_omega->ls()+ksfw0_omega->lb());
    }
  case 5:
    if(h0mode == 10){
      if(fsf){
        ksfw1_etapgg->input(kmm2, kvars);
        return ksfw1_etapgg->ls()/(ksfw1_etapgg->ls()+ksfw1_etapgg->lb());
      } else{
        ksfw0_etapgg->input(kmm2, kvars);
        return ksfw0_etapgg->ls()/(ksfw0_etapgg->ls()+ksfw0_etapgg->lb());
      }
    } else{
      if(fsf){
        ksfw1_etapppp->input(kmm2, kvars);
        return ksfw1_etapppp->ls()/(ksfw1_etapppp->ls()+ksfw1_etapppp->lb());
      } else{
        ksfw0_etapppp->input(kmm2, kvars);
        return ksfw0_etapppp->ls()/(ksfw0_etapppp->ls()+ksfw0_etapppp->lb());
      }
    }
  default:
    return -2;
  }
}

double filter::BDT(const int mode, const int h0mode){
  switch(mode){
  case 1:
    return reader_pi0->EvaluateMVA("pi0");
    break;
  case 2:
    if(h0mode == 10) return reader_gg->EvaluateMVA("gg");
    else             return reader_ppp->EvaluateMVA("etappp");
    break;
  case 3:
    return reader_omega->EvaluateMVA("omega");
    break;
  case 4:
    return reader_ppp->EvaluateMVA("etappp");
    break;
  case 5://eta'
    if(h0mode == 10) return reader_gg->EvaluateMVA("gg");
    else             return reader_ppp->EvaluateMVA("etappp");
    break;
  case 10:// D* pi0
    return reader_pi0->EvaluateMVA("pi0");
    break;
  case 20://D* eta
    if(h0mode == 10) return reader_gg->EvaluateMVA("gg");
    else             return reader_ppp->EvaluateMVA("etappp");
    break;
  }
  return -5;
}
