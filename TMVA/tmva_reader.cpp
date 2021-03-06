#include "../rooksfw/ksfw_pi0.h"
#include "../rooksfw/ksfw_etagg.h"
#include "../rooksfw/ksfw_etappp.h"
#include "../rooksfw/ksfw_omega.h"
#include "../rooksfw/rooksfw.h"
#include "../rooksfw/rooksfw.cc"

using namespace std;
void tmva_reader(const int type = 1,const bool second_iter = false){
// type = 0 -> Data
/////////////////
// TMVA Reader //
/////////////////
//  rooksfw* ksfw0  = new rooksfw("ksfw0",ksfw0_alpha,ksfw0_sigpdf,ksfw0_bkgpdf);
  if(!second_iter){
    rooksfw* ksfw1_pi0    = new rooksfw("ksfw1_pi0",ksfw1_alpha_pi0,ksfw1_sigpdf_pi0,ksfw1_bkgpdf_pi0);
    rooksfw* ksfw1_etagg  = new rooksfw("ksfw1_etagg",ksfw1_alpha_etagg,ksfw1_sigpdf_etagg,ksfw1_bkgpdf_etagg);
    rooksfw* ksfw1_etappp = new rooksfw("ksfw1_etappp",ksfw1_alpha_etappp,ksfw1_sigpdf_etappp,ksfw1_bkgpdf_etappp);
    rooksfw* ksfw1_omega  = new rooksfw("ksfw1_omega",ksfw1_alpha_omega,ksfw1_sigpdf_omega,ksfw1_bkgpdf_omega);
  }
  TMVA::Reader* reader_pi0 = new TMVA::Reader("!Color:!Silent:V");
  TMVA::Reader* reader_gg  = new TMVA::Reader("!Color:!Silent:V");
  TMVA::Reader* reader_ppp = new TMVA::Reader("!Color:!Silent:V");

  vector<string> fnames;
  TString fname;
  string ofname;
  string prefix1,prefix2;
  if(!second_iter){
    prefix1 = string("/home/vitaly/B0toDh0/Tuples/fil_");
    prefix2 = string("/home/vitaly/B0toDh0/TMVA/Fil_");
  } else{
    prefix1 = string("/home/vitaly/B0toDh0/TMVA/Fil_");
    prefix2 = string("/home/vitaly/B0toDh0/TMVA/Fil1_");   
  }
  switch(type){
    case 11:
      const string fstr = prefix1 + string("b2dh_sigmcPi0_s7.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_sigmcPi0_s7");
      break;
    case 12:
      const string fstr = prefix1 + string("b2dh_sigmcEta_s2.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_sigmcEta_s2");
      break;
    case 13:
      const string fstr = prefix1 + string("b2dh_sigmcOmega_s5.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_sigmcOmega_s5");
      break;
    case 14:
      const string fstr = prefix1 + string("b2dh_sigmcOmega_s5.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_sigRho_s1");
      break;
    case 15:// eta'
      const string fstr = prefix1 + string("b2dh_sigmcETAP_s1.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_sigmcETAP_s1");
      break;
    case 2:
      fnames.push_back(new string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_cont_0-1.root"));
      ofname = string("fil_b2dh_cont_0-1");
      break;
    case 21://uds
      const string fstr = prefix1 + string("b2dh_uds_2_12.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_uds_2_12");
      break;
    case 22://charm
      const string fstr = prefix1 + string("b2dh_charm_2_12.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_charm_2_12");
      break;
    case 23://charged
      const string fstr = prefix1 + string("b2dh_charged_2_12.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_charged_2_12");
      break;
    case 24://mixed
      const string fstr = prefix1 + string("b2dh_mixed_2_12.root");
      fnames.push_back(fstr.c_str());
      ofname = prefix2 + string("b2dh_mixed_2_12");
      break;
    case 3:
      fnames.push_back("/home/vitaly/B0toDh0/Tuples/fil_b2dh_gen_0-1.root");
      ofname = string("fil_b2dh_gen_0-1_nosig");
      break;
    case 4:// b-bbar generic
      fname = TString("/home/vitaly/B0toDh0/Tuples/fil_bb_0-10.root");
      cout << fname << endl;
      fnames.push_back(fname);
      fname = TString("/home/vitaly/B0toDh0/Tuples/fil_bb_1-11.root");
      cout << fname << endl;
      fnames.push_back(fname);
      ofname = string("fil_b2dh_bb_0-1");
      break;
    case 0:
      fnames.push_back(new string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_data.root"));
      ofname = string("fil_b2dh_data");
      break;
    case 111:
      fnames.push_back(new string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmc.root"));
      ofname = string("fil_b2dh_sig");
      break;
    default:
      cout << "Wrong type " << type << endl;
      return;
  }

//  tmststr = new TString("TEvent");
//  treenames->pushback((TObject*)fname);

  ofname += string(".root");
  TFile ofile(ofname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

//  double k0vars[17];
  double k1vars[17];
  Double_t cos_thr,thr_sig,thr_oth,p_pi0_h0,p_pip_h0,p_pim_h0,cos_hel;
//  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;
  Double_t k1mm2,k1et,k1hso00,k1hso02,k1hso04,k1hso10,k1hso12,k1hso14,k1hso20,k1hso22,k1hso24,k1hoo0,k1hoo1,k1hoo2,k1hoo3,k1hoo4;
  Double_t egamma;
  Double_t p_d0,dz_mc_sig, dz_mc_asc;

  Double_t mbc,de,mp,mm,dz,atckpi_max,mpi0,mh0,mk,md;
  Double_t bdtg,bdt;
  Double_t lh,lh1;// RooKSFW
  Double_t tag_LH,tag_LH_err;
  Int_t phsp,bin,exp,run,evtn;
  Double_t p_h0;

  Float_t m_costhBcms,m_cos_thr,m_thr_sig,m_egamma,m_lh1;
  Float_t m_chi2_mass_d0,m_h0_chi2;
  Float_t m_p_d0,m_p_pi0_h0,m_p_pip_h0,m_p_pim_h0;
  Float_t m_k1mm2,m_k1et,m_k1hso00,m_k1hso02,m_k1hso04,m_k1hso10,m_k1hso12,m_k1hso14,m_k1hso20,m_k1hso22,m_k1hso24,m_k1hoo0,m_k1hoo1,m_k1hoo2,m_k1hoo3,m_k1hoo4;

  Double_t mp_mc,mm_mc;
  Double_t mp_raw,mm_raw;
  Double_t d0_t_mc;
  Double_t dt_mc,dz_mc;
  Int_t bin_mc;
  Int_t flv_mc,d0_flv_mc;
  Int_t b0f,d0f;
  Int_t nptag;
  Int_t good_icpv,good_icpv_d0;
//  Double_t dz_pull_sig,dz_pull_asc;

  Int_t mode,h0mode,h0f,pi0f;
  Double_t z_sig,z_asc,pi0_chi2;
  Double_t sz_sig,sz_asc;
  Int_t ntrk_sig,ntrk_asc,ndf_z_sig,ndf_z_asc;
  Double_t chisq_z_sig,chisq_z_asc,cl_z_sig,cl_z_asc,h0_chi2;
  Double_t costhB,costhBcms,Ecms;
  Double_t t_sig_mc,z_sig_mc,t_asc_mc,z_asc_mc;
//  Double_t dz_pull_sig_d0;
//  Double_t z_sig_d0;
//  Double_t sz_sig_d0;
//  Double_t dz_mc_sig_d0;
  Int_t b0id,d0id,h0id,dst0id,dst0f,etapid,etapf;

  reader_gg->AddVariable("abs(costhBcms)",&m_costhBcms);
  reader_gg->AddVariable("log(chi2_mass_d0)",&m_chi2_mass_d0);
  reader_gg->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader_gg->AddVariable("thr_sig",&m_thr_sig);
  reader_gg->AddVariable("log(h0_chi2)",&m_h0_chi2);
  reader_gg->AddVariable("log(egamma)",&m_egamma);

  reader_pi0->AddVariable("abs(costhBcms)",&m_costhBcms);
  reader_pi0->AddVariable("log(chi2_mass_d0)",&m_chi2_mass_d0);
  reader_pi0->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader_pi0->AddVariable("thr_sig",&m_thr_sig);
  reader_pi0->AddVariable("log(egamma)",&m_egamma);

  reader_ppp->AddVariable("abs(costhBcms)",&m_costhBcms);
  reader_ppp->AddVariable("log(chi2_mass_d0)",&m_chi2_mass_d0);
  reader_ppp->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader_ppp->AddVariable("thr_sig",&m_thr_sig);
  reader_ppp->AddVariable("log(egamma)",&m_egamma);
  reader_ppp->AddVariable("log(p_pi0_h0)",&m_p_pi0_h0);
  reader_ppp->AddVariable("log(p_pip_h0)",&m_p_pip_h0);
  reader_ppp->AddVariable("log(p_pim_h0)",&m_p_pim_h0);

  if(!second_iter){
    reader_gg->AddVariable("k1mm2",&m_k1mm2);
    reader_gg->AddVariable("k1et",&m_k1et);
    reader_gg->AddVariable("k1hso00",&m_k1hso00);
    reader_gg->AddVariable("k1hso02",&m_k1hso02);
    reader_gg->AddVariable("k1hso04",&m_k1hso04);
    reader_gg->AddVariable("k1hso10",&m_k1hso10);
    reader_gg->AddVariable("k1hso12",&m_k1hso12);
    reader_gg->AddVariable("k1hso14",&m_k1hso14);
    reader_gg->AddVariable("k1hso20",&m_k1hso20);
    reader_gg->AddVariable("k1hso22",&m_k1hso22);
    reader_gg->AddVariable("k1hso24",&m_k1hso24);
    reader_gg->AddVariable("k1hoo0",&m_k1hoo0);
    reader_gg->AddVariable("k1hoo1",&m_k1hoo1);
    reader_gg->AddVariable("k1hoo2",&m_k1hoo2);
    reader_gg->AddVariable("k1hoo3",&m_k1hoo3);
    reader_gg->AddVariable("k1hoo4",&m_k1hoo4);

    reader_pi0->AddVariable("k1mm2",&m_k1mm2);
    reader_pi0->AddVariable("k1et",&m_k1et);
    reader_pi0->AddVariable("k1hso00",&m_k1hso00);
    reader_pi0->AddVariable("k1hso02",&m_k1hso02);
    reader_pi0->AddVariable("k1hso04",&m_k1hso04);
    reader_pi0->AddVariable("k1hso10",&m_k1hso10);
    reader_pi0->AddVariable("k1hso12",&m_k1hso12);
    reader_pi0->AddVariable("k1hso14",&m_k1hso14);
    reader_pi0->AddVariable("k1hso20",&m_k1hso20);
    reader_pi0->AddVariable("k1hso22",&m_k1hso22);
    reader_pi0->AddVariable("k1hso24",&m_k1hso24);
    reader_pi0->AddVariable("k1hoo0",&m_k1hoo0);
    reader_pi0->AddVariable("k1hoo1",&m_k1hoo1);
    reader_pi0->AddVariable("k1hoo2",&m_k1hoo2);
    reader_pi0->AddVariable("k1hoo3",&m_k1hoo3);
    reader_pi0->AddVariable("k1hoo4",&m_k1hoo4);

    reader_ppp->AddVariable("k1mm2",&m_k1mm2);
    reader_ppp->AddVariable("k1et",&m_k1et);
    reader_ppp->AddVariable("k1hso00",&m_k1hso00);
    reader_ppp->AddVariable("k1hso02",&m_k1hso02);
    reader_ppp->AddVariable("k1hso04",&m_k1hso04);
    reader_ppp->AddVariable("k1hso10",&m_k1hso10);
    reader_ppp->AddVariable("k1hso12",&m_k1hso12);
    reader_ppp->AddVariable("k1hso14",&m_k1hso14);
    reader_ppp->AddVariable("k1hso20",&m_k1hso20);
    reader_ppp->AddVariable("k1hso22",&m_k1hso22);
    reader_ppp->AddVariable("k1hso24",&m_k1hso24);
    reader_ppp->AddVariable("k1hoo0",&m_k1hoo0);
    reader_ppp->AddVariable("k1hoo1",&m_k1hoo1);
    reader_ppp->AddVariable("k1hoo2",&m_k1hoo2);
    reader_ppp->AddVariable("k1hoo3",&m_k1hoo3);
    reader_ppp->AddVariable("k1hoo4",&m_k1hoo4);
  } else{
    reader_pi0->AddVariable("lh1",&m_lh1);
    reader_gg->AddVariable("lh1",&m_lh1);
    reader_ppp->AddVariable("lh1",&m_lh1);
  }

  if(!second_iter){
    reader_pi0->BookMVA("pi0","weights/MVA_softcut_Dpi0_BDTG.weights.xml");
    reader_gg->BookMVA("gg","weights/MVA_softcut_Detagg_BDTG.weights.xml");
    reader_ppp->BookMVA("ppp","weights/MVA_softcut_D_etappp_omega_BDTG.weights.xml");
  } else{
    reader_pi0->BookMVA("pi0","weights/MVA_softcut_ksfw_Dpi0_BDT.weights.xml");
    reader_gg->BookMVA("gg","weights/MVA_softcut_ksfw_Detagg_BDT.weights.xml");
    reader_ppp->BookMVA("etappp","weights/MVA_softcut_ksfw_Detappp_BDT.weights.xml");
    reader_ppp->BookMVA("omega","weights/MVA_softcut_ksfw_Domega_BDT.weights.xml");
  }
//  reader_gg->BookMVA("Dpi0","weights/MVA_softcut_Dpi0_BDTG.weights.xml");
//  reader_ppp->BookMVA("Detagg","weights/MVA_softcut_Detagg_BDTG.weights.xml");
//  reader->BookMVA("Detappp","weights/MVA_softcut_Detappp_BDTG.weights.xml");
//  reader->BookMVA("Domega","weights/MVA_softcut_Domega_BDTG.weights.xml");

//  const int NTrees = treenames->GetEntries();
  const int NFiles = fnames.size();
//  cout << "NTrees: " << NTrees << ", NFiles: " << NFiles << endl;
//  TListIter treeiter(treenames);
//  TString* tstr;
  for(int filenn = 0; filenn<NFiles; filenn++){
//    treeiter.Reset();
//    for(int treenn = 0; treenn<NTrees; treenn++){
//    cout << "Cycle: " << filenn << " " << treenn << endl;
//  tstr = (TString*)treeiter.Next();
//  cout << " " << fnames[filenn] << " " << *tstr << endl;
  TFile *ifile = TFile::Open(fnames[filenn].c_str());
  TTree *tree = (TTree*)ifile->Get("TEvent");
  if(type>10 && type<20){
//    tree->SetBranchAddress("dz_mc_sig_d0",&dz_mc_sig_d0);
//    tree->SetBranchAddress("z_sig_d0",&z_sig_d0);
//    tree->SetBranchAddress("sz_sig_d0",&sz_sig_d0);
//    tree->SetBranchAddress("dz_pull_sig_d0",&dz_pull_sig_d0);
    tree->SetBranchAddress("mp_mc",&mp_mc);
    tree->SetBranchAddress("mm_mc",&mm_mc);
    tree->SetBranchAddress("mp_raw",&mp_raw);
    tree->SetBranchAddress("mm_raw",&mm_raw);
    tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
    tree->SetBranchAddress("dt_mc",&dt_mc);
    tree->SetBranchAddress("dz_mc",&dz_mc);
    tree->SetBranchAddress("bin_mc",&bin_mc);
    tree->SetBranchAddress("flv_mc",&flv_mc);
    tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);

    tree->SetBranchAddress("t_sig_mc",&t_sig_mc);
    tree->SetBranchAddress("z_sig_mc",&z_sig_mc);
    tree->SetBranchAddress("t_asc_mc",&t_asc_mc);
    tree->SetBranchAddress("z_asc_mc",&z_asc_mc);
    tree->SetBranchAddress("dz_mc_sig",&dz_mc_sig);
    tree->SetBranchAddress("dz_mc_asc",&dz_mc_asc);

//    tree->SetBranchAddress("dz_pull_sig",&dz_pull_sig);
//    tree->SetBranchAddress("dz_pull_asc",&dz_pull_asc);
  }

  Double_t chi2_mass_d0,chi2_vtx_d0;
  tree->SetBranchAddress("chi2_mass_d0",&chi2_mass_d0);
  tree->SetBranchAddress("chi2_vtx_d0",&chi2_vtx_d0);

  tree->SetBranchAddress("good_icpv",&good_icpv);
//  tree->SetBranchAddress("good_icpv_d0",&good_icpv_d0);

  if(type){
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
  }

  Double_t md, md_raw, md_fit, mdpip, mdpim;
  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("bin",&bin);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("phsp",&phsp);
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

  tree->SetBranchAddress("mode",&mode);
  tree->SetBranchAddress("h0mode",&h0mode);
  tree->SetBranchAddress("h0f",&h0f);
  tree->SetBranchAddress("pi0f",&pi0f);

  tree->SetBranchAddress("z_sig",&z_sig);
  tree->SetBranchAddress("z_asc",&z_asc);

  tree->SetBranchAddress("sz_sig",&sz_sig);
  tree->SetBranchAddress("sz_asc",&sz_asc);

  tree->SetBranchAddress("ntrk_sig",&ntrk_sig);
  tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  tree->SetBranchAddress("ndf_z_sig",&ndf_z_sig);
  tree->SetBranchAddress("ndf_z_asc",&ndf_z_asc);

  tree->SetBranchAddress("chisq_z_sig",&chisq_z_sig);
  tree->SetBranchAddress("chisq_z_asc",&chisq_z_asc);
  tree->SetBranchAddress("cl_z_sig",&cl_z_sig);
  tree->SetBranchAddress("cl_z_asc",&cl_z_asc);

  tree->SetBranchAddress("costhBcms",&costhBcms);
  tree->SetBranchAddress("costhB",&costhB);

  tree->SetBranchAddress("p_h0",&p_h0);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("cos_hel",&cos_hel);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);

  tree->SetBranchAddress("h0_chi2",&h0_chi2);
  tree->SetBranchAddress("pi0_chi2",&pi0_chi2);
  tree->SetBranchAddress("egamma",&egamma);

  tree->SetBranchAddress("p_pi0_h0",&p_pi0_h0);
  tree->SetBranchAddress("p_pip_h0",&p_pip_h0);
  tree->SetBranchAddress("p_pim_h0",&p_pim_h0);

//  TEvent->Branch("lh",&lh,"lh/D");
  TEvent->Branch("lh1",&lh1,"lh1/D");

  Double_t mdst0, metap, dmdst0, dmetap;
  tree->SetBranchAddress("mdst0",&mdst0);
  tree->SetBranchAddress("dmdst0",&dmdst0);
  tree->SetBranchAddress("metap",&metap);
  tree->SetBranchAddress("dmetap",&dmetap);

  TEvent->Branch("mdst0",&mdst0,"mdst0/D");
  TEvent->Branch("dmdst0",&dmdst0,"dmdst0/D");
  TEvent->Branch("metap",&metap,"metap/D");
  TEvent->Branch("dmetap",&dmetap,"dmetap/D");

  TEvent->Branch("p_d0",&p_d0,"p_d0/D");
  TEvent->Branch("thr_oth",&thr_oth,"thr_oth/D");
  TEvent->Branch("nptag",&nptag,"nptag/I");
  TEvent->Branch("phsp",&phsp,"phsp/I");

  TEvent->Branch("egamma",&egamma,"egamma/D");
  TEvent->Branch("p_h0",&p_h0,"p_h0/D");
  TEvent->Branch("cos_thr",&cos_thr,"cos_thr/D");
  TEvent->Branch("cos_hel",&cos_hel,"cos_hel/D");
  TEvent->Branch("thr_sig",&thr_sig,"thr_sig/D");
  TEvent->Branch("p_pi0_h0",&p_pi0_h0,"p_pi0_h0/D");
  TEvent->Branch("p_pip_h0",&p_pip_h0,"p_pip_h0/D");
  TEvent->Branch("p_pim_h0",&p_pim_h0,"p_pim_h0/D");
//  tree->SetBranchAddress("k0mm2",&k0mm2);
//  tree->SetBranchAddress("k0et",&k0et);
//  tree->SetBranchAddress("k0hso00",&k0hso00);
//  tree->SetBranchAddress("k0hso02",&k0hso02);
//  tree->SetBranchAddress("k0hso04",&k0hso04);
//  tree->SetBranchAddress("k0hso10",&k0hso10);
//  tree->SetBranchAddress("k0hso12",&k0hso12);
//  tree->SetBranchAddress("k0hso14",&k0hso14);
//  tree->SetBranchAddress("k0hso20",&k0hso20);
//  tree->SetBranchAddress("k0hso22",&k0hso22);
//  tree->SetBranchAddress("k0hso24",&k0hso24);
//  tree->SetBranchAddress("k0hoo0",&k0hoo0);
//  tree->SetBranchAddress("k0hoo1",&k0hoo1);
//  tree->SetBranchAddress("k0hoo2",&k0hoo2);
//  tree->SetBranchAddress("k0hoo3",&k0hoo3);
//  tree->SetBranchAddress("k0hoo4",&k0hoo4);

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

  if(!second_iter){
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

  TEvent->Branch("k1mm2",&k1mm2,"k1mm2/D");
  TEvent->Branch("k1et",&k1et,"k1et/D");
  TEvent->Branch("k1hso00",&k1hso00,"k1hso00/D");
  TEvent->Branch("k1hso02",&k1hso02,"k1hso02/D");
  TEvent->Branch("k1hso04",&k1hso04,"k1hso04/D");
  TEvent->Branch("k1hso10",&k1hso10,"k1hso10/D");
  TEvent->Branch("k1hso12",&k1hso12,"k1hso12/D");
  TEvent->Branch("k1hso14",&k1hso14,"k1hso14/D");
  TEvent->Branch("k1hso20",&k1hso20,"k1hso20/D");
  TEvent->Branch("k1hso22",&k1hso22,"k1hso22/D");
  TEvent->Branch("k1hso24",&k1hso24,"k1hso24/D");
  TEvent->Branch("k1hoo0",&k1hoo0,"k1hoo0/D");
  TEvent->Branch("k1hoo1",&k1hoo1,"k1hoo1/D");
  TEvent->Branch("k1hoo2",&k1hoo2,"k1hoo2/D");
  TEvent->Branch("k1hoo3",&k1hoo3,"k1hoo3/D");
  TEvent->Branch("k1hoo4",&k1hoo4,"k1hoo4/D");
  } else{
    tree->SetBranchAddress("bdtg",&bdtg);
    tree->SetBranchAddress("lh1",&lh1);
    TEvent->Branch("bdt",&bdt,"bdt/D");
  }

  if(type>10 && type<20){
//    TEvent->Branch("dz_mc_sig_d0",&dz_mc_sig_d0,"dz_mc_sig_d0/D");
//    TEvent->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
//    TEvent->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
//    TEvent->Branch("dz_pull_sig_d0",&dz_pull_sig_d0,"dz_pull_sig_d0/D");
    TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
    TEvent->Branch("mm_mc",&mm_mc,"mm_mc/D");
    TEvent->Branch("mp_raw",&mp_raw,"mp_raw/D");
    TEvent->Branch("mm_raw",&mm_raw,"mm_raw/D");
    TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
    TEvent->Branch("dt_mc",&dt_mc,"dt_mc/D");
    TEvent->Branch("dz_mc",&dz_mc,"dz_mc/D");
    TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
    TEvent->Branch("flv_mc",&flv_mc,"flv_mc/I");
    TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");

    TEvent->Branch("t_sig_mc",&t_sig_mc,"t_sig_mc/D");
    TEvent->Branch("z_sig_mc",&z_sig_mc,"z_sig_mc/D");
    TEvent->Branch("t_asc_mc",&t_asc_mc,"t_asc_mc/D");
    TEvent->Branch("z_asc_mc",&z_asc_mc,"z_asc_mc/D");
    TEvent->Branch("dz_mc_sig",&dz_mc_sig,"dz_mc_sig/D");
    TEvent->Branch("dz_mc_asc",&dz_mc_asc,"dz_mc_asc/D");

//    TEvent->Branch("dz_pull_sig",&dz_pull_sig,"dz_pull_sig/D");
//    TEvent->Branch("dz_pull_asc",&dz_pull_asc,"dz_pull_asc/D");
  }

  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");
  TEvent->Branch("nptag",&nptag,"nptag/I");

  if(type){
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
  }

  TEvent->Branch("bdtg",&bdtg,"bdtg/D");
  TEvent->Branch("good_icpv",&good_icpv,"good_icpv/I");
//  TEvent->Branch("good_icpv_d0",&good_icpv_d0,"good_icpv_d0/I");

  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");
  TEvent->Branch("bin",&bin,"bin/I");

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");

  TEvent->Branch("chi2_mass_d0",&chi2_mass_d0,"chi2_mass_d0/D");
  TEvent->Branch("chi2_vtx_d0",&chi2_vtx_d0,"chi2_vtx_d0/D");

  TEvent->Branch("dz",&dz,"dz/D");

  TEvent->Branch("tag_LH",&tag_LH,"tag_LH/D");
  TEvent->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");

  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("mh0",&mh0,"mh0/D");
  TEvent->Branch("mpi0",&mpi0,"mpi0/D");
  TEvent->Branch("mk",&mk,"mk/D");
  TEvent->Branch("md",&md,"md/D");
  TEvent->Branch("md_raw",&md_raw,"md_raw/D");
  TEvent->Branch("md_fit",&md_fit,"md_fit/D");
  TEvent->Branch("mdpip",&mdpip,"mdpip/D");
  TEvent->Branch("mdpim",&mdpim,"mdpim/D");

  TEvent->Branch("mode",&mode,"mode/I");
  TEvent->Branch("h0mode",&h0mode,"h0mode/I");

  TEvent->Branch("z_sig",&z_sig,"z_sig/D");
  TEvent->Branch("z_asc",&z_asc,"z_asc/D");

  TEvent->Branch("sz_sig",&sz_sig,"sz_sig/D");
  TEvent->Branch("sz_asc",&sz_asc,"sz_asc/D");

  TEvent->Branch("ntrk_sig",&ntrk_sig,"ntrk_sig/I");
  TEvent->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  TEvent->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  TEvent->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");

  TEvent->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  TEvent->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
  TEvent->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
  TEvent->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");
  TEvent->Branch("h0_chi2",&h0_chi2,"h0_chi2/D");
  TEvent->Branch("pi0_chi2",&pi0_chi2,"pi0_chi2/D");

  TEvent->Branch("costhB",&costhB,"costhB/D");
  TEvent->Branch("costhBcms",&costhBcms,"costhBcms/D");
  TEvent->Branch("Ecms",&Ecms,"Ecm/D");

  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
//    if(type == 3 && (b0f == 1 || b0f == 5 || b0f == 10)) continue;
    if(chi2_mass_d0<0) continue;
    if(!(i%10000)) cout << i << " events" << endl;

    m_costhBcms      = (float)TMath::Abs(costhBcms);
//    m_chi2_vtx_d0 = (float)TMath::Log(chi2_vtx_d0);
    m_chi2_mass_d0= (float)TMath::Log(chi2_mass_d0);
    m_cos_thr     = (float)TMath::Abs(cos_thr);
    m_thr_sig     = (float)thr_sig;
    m_h0_chi2     = (float)TMath::Log(h0_chi2);
 //   m_thr_oth     = (float)thr_oth;
//    m_p_d0        = (float)p_d0;
    m_egamma      = (float)TMath::Log(egamma);
    if(p_pi0_h0>0) m_p_pi0_h0 = (float)TMath::Log(p_pi0_h0);
    if(p_pip_h0>0) m_p_pip_h0 = (float)TMath::Log(p_pip_h0);
    if(p_pim_h0>0) m_p_pim_h0 = (float)TMath::Log(p_pim_h0);
    if(!second_iter){
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

//    k0vars[0]  = k0et;
//    k0vars[1]  = k0hso00;
//    k0vars[2]  = k0hso10;
//    k0vars[3]  = k0hso20;
//    k0vars[4]  = 0;//k0hso01;
//    k0vars[5]  = k0hso02;
//    k0vars[6]  = k0hso12;
//    k0vars[7]  = k0hso22;
//    k0vars[8]  = 0;//k0hso03;
//    k0vars[9]  = k0hso04;
//    k0vars[10] = k0hso14;
//    k0vars[11] = k0hso24;
//    k0vars[12] = k0hoo0;
//    k0vars[13] = k0hoo1;
//    k0vars[14] = k0hoo2;
//    k0vars[15] = k0hoo3;
//    k0vars[16] = k0hoo4;

    k1vars[0]  = k1et;
    k1vars[1]  = k1hso00;
    k1vars[2]  = k1hso10;
    k1vars[3]  = k1hso20;
    k1vars[4]  = 0;//k1hso01;
    k1vars[5]  = k1hso02;
    k1vars[6]  = k1hso12;
    k1vars[7]  = k1hso22;
    k1vars[8]  = 0;//k1hso03;
    k1vars[9]  = k1hso04;
    k1vars[10] = k1hso14;
    k1vars[11] = k1hso24;
    k1vars[12] = k1hoo0;
    k1vars[13] = k1hoo1;
    k1vars[14] = k1hoo2;
    k1vars[15] = k1hoo3;
    k1vars[16] = k1hoo4;

    if(mode == 1){
      ksfw1_pi0->input(k1mm2, k1vars);
      if(ksfw1_pi0->ls() + ksfw1_pi0->lb() > 0){
        lh1 = ksfw1_pi0->ls()/(ksfw1_pi0->ls()+ksfw1_pi0->lb());
      } else{
        lh1 = -1;
      }
    } else if(mode == 2 && h0mode == 10){
      ksfw1_etagg->input(k1mm2, k1vars);
      if(ksfw1_etagg->ls() + ksfw1_etagg->lb() > 0){
        lh1 = ksfw1_etagg->ls()/(ksfw1_etagg->ls()+ksfw1_etagg->lb());
      } else{
        lh1 = -1;
      }
    } else if(mode == 2 && h0mode == 20){
      ksfw1_etappp->input(k1mm2, k1vars);
      if(ksfw1_etappp->ls() + ksfw1_etappp->lb() > 0){
        lh1 = ksfw1_etappp->ls()/(ksfw1_etappp->ls()+ksfw1_etappp->lb());
      } else{
        lh1 = -1;
      }
    } else{
      ksfw1_omega->input(k1mm2, k1vars);
      if(ksfw1_omega->ls() + ksfw1_omega->lb() > 0){
        lh1 = ksfw1_omega->ls()/(ksfw1_omega->ls()+ksfw1_omega->lb());
      } else{
        lh1 = -1;
      }
    }

    if(mode == 1)         bdtg = reader_pi0->EvaluateMVA("pi0");
    else if(h0mode == 20) bdtg = reader_ppp->EvaluateMVA("ppp");
    else                  bdtg = reader_gg->EvaluateMVA("gg");
    } else{// second iteration
      m_lh1 = (float)lh1;
      switch(mode){
      case 1:
        bdt = reader_pi0->EvaluateMVA("pi0");
        break;
      case 2:
        if(h0mode == 10) bdt = reader_gg->EvaluateMVA("gg");
        else             bdt = reader_ppp->EvaluateMVA("etappp");
        break;
      case 3:
        bdt = reader_ppp->EvaluateMVA("omega");
        break;
      case 4:
        bdt = reader_ppp->EvaluateMVA("omega");
        break;
      case 5://eta'
        if(h0mode == 10) bdt = reader_gg->EvaluateMVA("gg");
        else             bdt = reader_ppp->EvaluateMVA("etappp");
        break;
      }
    }
    TEvent->Fill();
  }
  }
  TEvent->Write();
  ofile.Write();
  ofile.Close();
}

