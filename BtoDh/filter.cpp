#include "cuts.h"

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

void filter(const int type){
  int M_MODE = 0;
  TChain* tree = new TChain("TEvent");
  string fname;
  switch(type){
  case 0://Data
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_data.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_data.root");
    break;
  case 11://Signal MC pi0
    M_MODE = 1;
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcPI0_s7.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcPi0_s7.root");
    break;
  case 12://Signal MC eta
    M_MODE = 2;
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcETA_s2.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcEta_s2.root");
    break;
  case 13://Signal MC omega
    M_MODE = 3;
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcOMEGA_s5.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcOmega_s5.root");
    break;
  case 14://Signal MC rho
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcRHO_s1.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcRho_s1.root");
    break;
  case 15://Signal MC eta'
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmcETAP_s1.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcETAP_s1.root");
    break;
  case 21://charged
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_2.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_12.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_charged_2_12.root");
    break;
  case 22://mixed
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_2.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_12.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_mixed_2_12.root");
    break;
  case 23://charm
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charm_12.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charm_2.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_charm_2_12.root");
    break;
  case 24://uds
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_uds_12.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_uds_2.root");
    fname = string("/home/vitaly/B0toDh0/Tuples/fil_b2dh_uds_2_12.root");
    break;
  case 2://Generic MC
//    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_gen.root");
//    stringstream out;
//    string path("/home/vitaly/B0toDh0/Tuples/");
    vector<string> types;
    types.push_back(string("charged"));
    types.push_back(string("mixed"));
    types.push_back(string("uds"));
    types.push_back(string("charm"));
    for(int i=0; i<types.size(); i++){
      for(int j=0; j<1; j++){
        out.str("");
        out << path << "b2dh_" << types[i] << "_" << j << ".root";
        tree->Add(out.str().c_str());
      }
    }
    fname = string("fil_b2dh_gen_0_2.root");
    break;
  case 3://Continuum
//    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_cont.root");
    stringstream out;
    string path("/home/vitaly/B0toDh0/Tuples/");
    vector<string> types;
    types.push_back(string("uds"));
    types.push_back(string("charm"));
    for(int i=0; i<types.size(); i++){
      for(int j=0; j<6; j++){
        out.str("");
	    out << path << "b2dh_" << types[i] << "_" << j << ".root";
	    tree->Add(out.str().c_str());
	  }
    }
    fname = string("fil_b2dh_cont.root");
    break;
  case 4://Custom
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_10.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_0.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_10.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_0.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_11.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_charged_1.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_11.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_mixed_1.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/fil_b2dh_charm_0_10.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/fil_b2dh_uds_0_10.root");
    fname = string("fil_b2dh_gen.root");
    break;
  default:
    return;
  }
  Double_t p_d0,p_h0,p_ks,p_pi0_h0,p_pip_h0,p_pim_h0,pi0_chi2,egamma,cos_thr,cos_hel,thr_sig,thr_oth;
  Int_t ndf_tag_vtx,phsp,bin,exp,run,evtn;
  Int_t b0f,d0f;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;
  Double_t k1mm2,k1et,k1hso00,k1hso01,k1hso02,k1hso03,k1hso04,k1hso10,k1hso12,k1hso14,k1hso20,k1hso22,k1hso24,k1hoo0,k1hoo1,k1hoo2,k1hoo3,k1hoo4;

  Double_t mbc,de,mp,mm,dz,atckpi_max,mpi0,mh0,mk;
  Double_t ks_dr,ks_dz,ks_dphi,ks_fl,tag_LH,tag_LH_err;
  Double_t mp_mc,mm_mc,d0_t_mc,dt_mc,dz_mc,bin_mc,flv_mc,d0_flv_mc;
//  Double_t z_sig_d0;
//  Double_t sz_sig_d0;
//  Double_t dz_mc_sig_d0;
  Double_t chi2_vtx_d0, chi2_mass_d0;
  Int_t good_icpv;//,good_icpv_d0;

  Int_t mode,h0mode,h0f,pi0f;
  Double_t z_sig,z_asc;
  Double_t sz_sig,sz_asc;
  Int_t ntrk_sig,ntrk_asc,ndf_z_sig,ndf_z_asc;
  Double_t chisq_z_sig,chisq_z_asc,h0_chi2;
  Double_t cl_z_sig,cl_z_asc;
  Double_t costhB,costhBcms,Ecms;
  Double_t t_sig_mc,z_sig_mc,t_asc_mc,z_asc_mc;
  Double_t dz_mc_sig, dz_mc_asc;
  Double_t z_upsilon;
//  Double_t dz_pull_sig,dz_pull_sig_d0,dz_pull_asc;
  Int_t nptag;

  tree->Print();

  TFile *ofile = new TFile(fname.c_str(),"RECREATE");
  TTree* TEvent = new TTree("TEvent","TEvent");

  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);

  tree->SetBranchAddress("p_d0",&p_d0);
  tree->SetBranchAddress("p_h0",&p_h0);
  tree->SetBranchAddress("p_ks",&p_ks);
  tree->SetBranchAddress("p_pi0_h0",&p_pi0_h0);
  tree->SetBranchAddress("p_pip_h0",&p_pip_h0);
  tree->SetBranchAddress("p_pim_h0",&p_pim_h0);
  tree->SetBranchAddress("egamma",&egamma);
  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("cos_hel",&cos_hel);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("ndf_tag_vtx",&ndf_tag_vtx);
  tree->SetBranchAddress("thr_oth",&thr_oth);
  tree->SetBranchAddress("phsp",&phsp);
  tree->SetBranchAddress("bin",&bin);

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

  Int_t b0id,d0id,h0id,dst0id,dst0f,etapid,etapf;

  tree->SetBranchAddress("chi2_vtx_d0",&chi2_vtx_d0);
  tree->SetBranchAddress("chi2_mass_d0",&chi2_mass_d0);

  TEvent->Branch("chi2_vtx_d0",&chi2_vtx_d0,"chi2_vtx_d0/D");
  TEvent->Branch("chi2_mass_d0",&chi2_mass_d0,"chi2_mass_d0/D");
  TEvent->Branch("good_icpv",&good_icpv,"good_icpv/I");

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

    Double_t mp_raw, mm_raw;
    if(type>10 && type<20 || type == 111){
      tree->SetBranchAddress("mp_mc",&mp_mc);
      tree->SetBranchAddress("mm_mc",&mm_mc);
      tree->SetBranchAddress("mp_raw",&mp_raw);
      tree->SetBranchAddress("mm_raw",&mm_raw);
      tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
      tree->SetBranchAddress("z_upsilon",&z_upsilon);
      tree->SetBranchAddress("bin_mc",&bin_mc);
      tree->SetBranchAddress("flv_mc",&flv_mc);
      tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);

      tree->SetBranchAddress("t_sig_mc",&t_sig_mc);
      tree->SetBranchAddress("z_sig_mc",&z_sig_mc);
      tree->SetBranchAddress("t_asc_mc",&t_asc_mc);
      tree->SetBranchAddress("z_asc_mc",&z_asc_mc);

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
//      if(type == 12 || type == 13 || type == 14){
//        tree->SetBranchAddress("z_sig_d0",&z_sig_d0);
//        tree->SetBranchAddress("sz_sig_d0",&sz_sig_d0);

//        TEvent->Branch("good_icpv_d0",&good_icpv_d0,"good_icpv_d0/I");
//        TEvent->Branch("dz_mc_sig_d0",&dz_mc_sig_d0,"dz_mc_sig_d0/D");
//        TEvent->Branch("dz_pull_sig_d0",&dz_pull_sig_d0,"dz_pull_sig_d0/D");
//        TEvent->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
//        TEvent->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
//      }
//      TEvent->Branch("dz_pull_asc",&dz_pull_asc,"dz_pull_asc/D");
      TEvent->Branch("bin_mc",&bin_mc,"bin_mc/I");
      TEvent->Branch("flv_mc",&flv_mc,"flv_mc/I");
      TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");

      TEvent->Branch("t_sig_mc",&t_sig_mc,"t_sig_mc/D");
      TEvent->Branch("z_sig_mc",&z_sig_mc,"z_sig_mc/D");
      TEvent->Branch("t_asc_mc",&t_asc_mc,"t_asc_mc/D");
      TEvent->Branch("z_asc_mc",&z_asc_mc,"z_asc_mc/D");

    }
  }

  tree->SetBranchAddress("nptag",&nptag);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("mh0_raw",&mh0);
  tree->SetBranchAddress("mpi0_raw",&mpi0);
  tree->SetBranchAddress("mks_raw",&mk);

  Double_t md, md_raw, md_fit, mdpip, mdpim;
  tree->SetBranchAddress("md0_raw",&md_raw);
  tree->SetBranchAddress("md0_fit",&md_fit);
  tree->SetBranchAddress("md0",&md);
  tree->SetBranchAddress("md0pip",&mdpip);
  tree->SetBranchAddress("md0pim",&mdpim);

  TEvent->Branch("md",&md,"md/D");
  TEvent->Branch("md_raw",&md_raw,"md_raw/D");
  TEvent->Branch("md_fit",&md_fit,"md_fit/D");
  TEvent->Branch("mdpip",&mdpip,"mdpip/D");
  TEvent->Branch("mdpim",&mdpim,"mdpim/D");

  Double_t mdst0, metap, dmdst0, dmetap;
  tree->SetBranchAddress("mdst0",&mdst0);
  tree->SetBranchAddress("dmdst0",&dmdst0);
  tree->SetBranchAddress("metap",&metap);
  tree->SetBranchAddress("dmetap",&dmetap);

  TEvent->Branch("mdst0",&mdst0,"mdst0/D");
  TEvent->Branch("dmdst0",&dmdst0,"dmdst0/D");
  TEvent->Branch("metap",&metap,"metap/D");
  TEvent->Branch("dmetap",&dmetap,"dmetap/D");

  tree->SetBranchAddress("mode",&mode);
  tree->SetBranchAddress("h0mode",&h0mode);

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

  tree->SetBranchAddress("pi0_chi2",&pi0_chi2);
  tree->SetBranchAddress("h0_chi2",&h0_chi2);

  TEvent->Branch("h0_chi2",&h0_chi2,"h0_chi2/D");
  TEvent->Branch("pi0_chi2",&pi0_chi2,"pi0_chi2/D");

  tree->SetBranchAddress("costhB",&costhB);
  tree->SetBranchAddress("costhBcms",&costhBcms);
  tree->SetBranchAddress("Ecms",&Ecms);

  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);

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
  tree->SetBranchAddress("k1hso01",&k1hso01);
  tree->SetBranchAddress("k1hso02",&k1hso02);
  tree->SetBranchAddress("k1hso03",&k1hso03);
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

  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");

  TEvent->Branch("p_d0",&p_d0,"p_d0/D");
  TEvent->Branch("p_h0",&p_h0,"p_h0/D");
  TEvent->Branch("p_ks",&p_ks,"p_ks/D");
  TEvent->Branch("p_pi0_h0",&p_pi0_h0,"p_pi0_h0/D");
  TEvent->Branch("p_pip_h0",&p_pip_h0,"p_pip_h0/D");
  TEvent->Branch("p_pim_h0",&p_pim_h0,"p_pim_h0/D");
  TEvent->Branch("egamma",&egamma,"egamma/D");
  TEvent->Branch("cos_thr",&cos_thr,"cos_thr/D");
  TEvent->Branch("cos_hel",&cos_hel,"cos_hel/D");
  TEvent->Branch("thr_sig",&thr_sig,"thr_sig/D");
  TEvent->Branch("thr_oth",&thr_oth,"thr_oth/D");
  TEvent->Branch("ndf_tag_vtx",&ndf_tag_vtx,"ndf_tag_vtx/I");
  TEvent->Branch("thr_oth",&thr_oth,"thr_oth/D");
  TEvent->Branch("phsp",&phsp,"phsp/I");
  TEvent->Branch("bin",&bin,"bin/I");

  TEvent->Branch("nptag",&nptag,"nptag/I");
  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");
  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");
  TEvent->Branch("dz",&dz,"dz/D");

  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("mh0",&mh0);
  TEvent->Branch("mpi0",&mpi0,"mpi0/D");
  TEvent->Branch("mk",&mk,"mk/D");

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

  TEvent->Branch("costhB",&costhB,"costhB/D");
  TEvent->Branch("costhBcms",&costhBcms,"costhBcms/D");
  TEvent->Branch("Ecms",&Ecms,"Ecm/D");

  TEvent->Branch("ks_dr",&ks_dr,"ks_dr/D");
  TEvent->Branch("ks_dz",&ks_dz,"ks_dz/D");
  TEvent->Branch("ks_dphi",&ks_dphi,"ks_dphi/D");
  TEvent->Branch("ks_fl",&ks_fl,"ks_fl/D");
  TEvent->Branch("tag_LH",&tag_LH,"tag_LH/D");
  TEvent->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");

  TEvent->Branch("k0mm2",&k0mm2,"k0mm2/D");
  TEvent->Branch("k0et",&k0et,"k0et/D");
  TEvent->Branch("k0hso00",&k0hso00,"k0hso00/D");
  TEvent->Branch("k0hso02",&k0hso02,"k0hso02/D");
  TEvent->Branch("k0hso04",&k0hso04,"k0hso04/D");
  TEvent->Branch("k0hso10",&k0hso10,"k0hso10/D");
  TEvent->Branch("k0hso12",&k0hso12,"k0hso12/D");
  TEvent->Branch("k0hso14",&k0hso14,"k0hso14/D");
  TEvent->Branch("k0hso20",&k0hso20,"k0hso20/D");
  TEvent->Branch("k0hso22",&k0hso22,"k0hso22/D");
  TEvent->Branch("k0hso24",&k0hso24,"k0hso24/D");
  TEvent->Branch("k0hoo0",&k0hoo0,"k0hoo0/D");
  TEvent->Branch("k0hoo1",&k0hoo1,"k0hoo1/D");
  TEvent->Branch("k0hoo2",&k0hoo2,"k0hoo2/D");
  TEvent->Branch("k0hoo3",&k0hoo3,"k0hoo3/D");
  TEvent->Branch("k0hoo4",&k0hoo4,"k0hoo4/D");

  TEvent->Branch("k1mm2",&k1mm2,"k1mm2/D");
  TEvent->Branch("k1et",&k1et,"k1et/D");
  TEvent->Branch("k1hso00",&k1hso00,"k1hso00/D");
  TEvent->Branch("k1hso01",&k1hso01,"k1hso01/D");
  TEvent->Branch("k1hso02",&k1hso02,"k1hso02/D");
  TEvent->Branch("k1hso03",&k1hso03,"k1hso03/D");
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

  const int NTot = tree->GetEntries();
  for(int i=1; i<NTot; i++){
    if(!(i%10000)){ cout << i << " events" << endl;}
    tree->GetEvent(i);
    if((b0f == 0 || b0f<-1) && type) continue;
    if(M_MODE){
      if(M_MODE != mode) continue;
    }
//    if(type == 11){
//      if((b0f != 1) && (b0f != 5) && (b0f != 10)) continue;
//    }
//    if(type == 11 || type == 12 || type == 13 || type == 14){
//      if(b0f < 1) continue;
//    }
//    if(type == 3){
//      if(!(b0f != 1 && b0f != 10 && !(d0f == 1 && (b0f == 5 || b0f == 4)))) continue;
//    }
//    if(type == 2){

//    }
//    if(md_raw<(DMass-md_cut) || md_raw>(DMass+md_cut)) continue;
//    if(!(mode == 2 && h0mode == 10) && mode != 14){
//      if(mpi0<(Pi0Mass-mpi0_cut) || mpi0>(Pi0Mass+mpi0_cut)) continue;
//    }
    if(mk<(KMass-mk_cut) || mk>(KMass+mk_cut)) continue;
//    if(chi2_vtx_d0>1000) continue;

    z_sig *= 10; z_asc *=10;// z_sig_d0 *= 10;
//    dz_mc_sig_d0   = z_sig_d0-z_sig_mc;
    dz_mc_sig      = z_sig-z_sig_mc;
    dz_mc_asc      = z_asc-z_asc_mc;
    dz = z_sig - z_asc;
    dt_mc = t_sig_mc - t_asc_mc;
    dz_mc = z_sig_mc - z_asc_mc;
    sz_sig      = 10.*TMath::Sqrt(sz_sig);
//    sz_sig_d0   = 10.*TMath::Sqrt(sz_sig_d0);
    sz_asc      = 10.*TMath::Sqrt(sz_asc);
//    dz_pull_sig      = dz_mc_sig/sz_sig;
//    dz_pull_sig_d0   = dz_mc_sig_d0/sz_sig_d0;
//    dz_pull_asc      = dz_mc_asc/sz_asc;

    // * Standatd ICPV cuts * //
    good_icpv    = IsGoodICPV(ndf_z_sig,sz_sig,chisq_z_sig,ndf_z_asc,sz_asc,chisq_z_asc);
//    good_icpv_d0 = IsGoodICPV(0,sz_sig_d0,1,ndf_z_asc,sz_asc,chisq_z_asc);
    // * ////////////////// * //

    TEvent->Fill();
  }
  TEvent->Print();
  TEvent->Write();
  ofile->Close();
  return;
}

