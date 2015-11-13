#include "icpvevent.h"

ICPVEvent::ICPVEvent(int type,const bool second_iter):
  m_type(type),m_second_iter(second_iter),
  exp(0),run(0),evtn(0),
  p_d0(0),p_h0(0),p_pi0_h0(0),p_pip_h0(0),p_pim_h0(0),egamma(0),
  cos_thr(-2),cos_hel(-2),thr_sig(-99),thr_oth(-99),
  phsp(-1),bin(0),
  hel_h0(0),hel_pi0(0),e_g1(0),e_g2(0),e_g3(0),e_g4(0),th_g1(-2),th_g2(-2),th_g3(-2),th_g4(-2),
  r_pip(-99),r_pim(-99),r_pi1(-99),r_pi2(-99),z_pip(-99),z_pim(-99),z_pi1(-99),z_pi2(-99),
  pt_pip(-99),pt_pim(-99),pt_pi1(-99),pt_pi2(-99),px_pim(-99),py_pim(-99),pz_pim(-99),
  px_pip(-99),py_pip(-99),pz_pip(-99),px_ks(-99),py_ks(-99),pz_ks(-99),
  chi2_vtx_d0(-1),chi2_mass_d0(-1),
  t_sig_mc(-99),z_sig_mc(-99),t_asc_mc(-99),z_asc_mc(-99),
  mp_raw(0),mm_raw(0),mp_mc(0),mm_mc(0),
  d0_t_mc(-99),dt_mc(-99),dz_mc(-99),dz_mc_sig(-99),dz_mc_asc(-99),
  b0f(0),d0f(0),h0f(0),pi0f(0),bin_mc(0),flv_mc(0),nptag(-1),d0ch0(0),d0ch1(0),d0ch2(0),d0ch3(0),h0ch0(0),h0ch1(0),
  rndm_pi0(-1),b0id(0),d0id(0),h0id(0),dst0id(0),dst0f(0),etapid(0),etapf(0),d0_flv_mc(0),
  mbc(0),de(-99),mp(0),mm(0),dz(-99),atckpi_max(-2),mpi0(0),mh0(0),mk(0),
  md(0),md_raw(0),md_fit(0),mdpip(0),mdpim(0),mdst0(0),metap(0),dmdst0(0),dmetap(0),
  mode(0),h0mode(0),
  z_sig(-99),z_asc(-99),sz_sig(-99),sz_asc(-99),
  ntrk_sig(-1),ntrk_asc(-1),ndf_z_sig(-1),ndf_z_asc(-1),
  chisq_z_sig(-1),chisq_z_asc(-1),cl_z_sig(-1),cl_z_asc(-1),
  h0_chi2(-1),pi0_chi2(-1),costhBcms(-2),tag_LH(-2),tag_LH_err(-1),flv(0),good_icpv(-1),lh0(-2),lh1(-2),bdt(-2)
{
  for(int i=0; i<8; i++){d0_chain[i] = 0; h0_chain[i] = 0;}
  k0mm2 = 0; k1mm2 = 0;
  for(int i=0; i<17; i++){k0vars[i] = 0; k1vars[i] = 0;}
}

void ICPVEvent::SetBrAddresses(TTree* tree){
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
  tree->SetBranchAddress("hel_h0",&hel_h0);
  tree->SetBranchAddress("hel_pi0",&hel_pi0);
  tree->SetBranchAddress("e_g1",&e_g1);
  tree->SetBranchAddress("e_g2",&e_g2);
  tree->SetBranchAddress("e_g3",&e_g3);
  tree->SetBranchAddress("e_g4",&e_g4);
  tree->SetBranchAddress("th_g1",&th_g1);
  tree->SetBranchAddress("th_g2",&th_g2);
  tree->SetBranchAddress("th_g3",&th_g3);
  tree->SetBranchAddress("th_g4",&th_g4);
  tree->SetBranchAddress("r_pip",&r_pip);
  tree->SetBranchAddress("r_pim",&r_pim);
  tree->SetBranchAddress("r_pi1",&r_pi1);
  tree->SetBranchAddress("r_pi2",&r_pi2);
  tree->SetBranchAddress("z_pip",&z_pip);
  tree->SetBranchAddress("z_pim",&z_pim);
  tree->SetBranchAddress("z_pi1",&z_pi1);
  tree->SetBranchAddress("z_pi2",&z_pi2);
  tree->SetBranchAddress("pt_pip",&pt_pip);
  tree->SetBranchAddress("pt_pim",&pt_pim);
  tree->SetBranchAddress("pt_pi1",&pt_pi1);
  tree->SetBranchAddress("pt_pi2",&pt_pi2);
  tree->SetBranchAddress("px_pim",&px_pim);
  tree->SetBranchAddress("py_pim",&py_pim);
  tree->SetBranchAddress("pz_pim",&pz_pim);
  tree->SetBranchAddress("px_pip",&px_pip);
  tree->SetBranchAddress("py_pip",&py_pip);
  tree->SetBranchAddress("pz_pip",&pz_pip);
  tree->SetBranchAddress("px_ks",&px_ks);
  tree->SetBranchAddress("py_ks",&py_ks);
  tree->SetBranchAddress("pz_ks",&pz_ks);
  tree->SetBranchAddress("chi2_vtx_d0",&chi2_vtx_d0);
  tree->SetBranchAddress("chi2_mass_d0",&chi2_mass_d0);
  if(m_type){
    tree->SetBranchAddress("nptag",&nptag);
    if(!m_second_iter){
      tree->SetBranchAddress("d0_chain",&d0_chain);
      tree->SetBranchAddress("h0_chain",&h0_chain);
    }
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
    if(m_type == 1){
      tree->SetBranchAddress("mp_mc",&mp_mc);
      tree->SetBranchAddress("mm_mc",&mm_mc);
      tree->SetBranchAddress("mp_raw",&mp_raw);
      tree->SetBranchAddress("mm_raw",&mm_raw);
      tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
      tree->SetBranchAddress("flv_mc",&flv_mc);
      tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
      tree->SetBranchAddress("t_sig_mc",&t_sig_mc);
      tree->SetBranchAddress("z_sig_mc",&z_sig_mc);
      tree->SetBranchAddress("t_asc_mc",&t_asc_mc);
      tree->SetBranchAddress("z_asc_mc",&z_asc_mc);
    }
  }
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  if(!m_second_iter){
    tree->SetBranchAddress("mks_raw",&mk);
    tree->SetBranchAddress("md0pip",&mdpip);
    tree->SetBranchAddress("md0pim",&mdpim);
    tree->SetBranchAddress("md0_raw",&md_raw);
    tree->SetBranchAddress("md0",&md);
    tree->SetBranchAddress("md0_fit",&md_fit);
    tree->SetBranchAddress("mh0_raw",&mh0);
    tree->SetBranchAddress("mpi0_raw",&mpi0);
  } else{
    tree->SetBranchAddress("mk",&mk);
    tree->SetBranchAddress("mdpip",&mdpip);
    tree->SetBranchAddress("mdpim",&mdpim);
    tree->SetBranchAddress("md_raw",&md_raw);
    tree->SetBranchAddress("md",&md);
    tree->SetBranchAddress("md_fit",&md_fit);
    tree->SetBranchAddress("mh0",&mh0);
    tree->SetBranchAddress("mpi0",&mpi0);
  }
  tree->SetBranchAddress("mdst0",&mdst0);
  tree->SetBranchAddress("metap",&metap);
  tree->SetBranchAddress("dmetap",&dmetap);
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
  tree->SetBranchAddress("costhBcms",&costhBcms);
  tree->SetBranchAddress("tag_LH",&tag_LH);
  tree->SetBranchAddress("tag_LH_err",&tag_LH_err);
  tree->SetBranchAddress("k0mm2",  &k0mm2);
  tree->SetBranchAddress("k0et",   &k0vars[0]);
  tree->SetBranchAddress("k0hso00",&k0vars[1]);
  tree->SetBranchAddress("k0hso10",&k0vars[2]);
  tree->SetBranchAddress("k0hso20",&k0vars[3]);
  tree->SetBranchAddress("k0hso02",&k0vars[5]);
  tree->SetBranchAddress("k0hso12",&k0vars[6]);
  tree->SetBranchAddress("k0hso22",&k0vars[7]);
  tree->SetBranchAddress("k0hso04",&k0vars[9]);
  tree->SetBranchAddress("k0hso14",&k0vars[10]);
  tree->SetBranchAddress("k0hso24",&k0vars[11]);
  tree->SetBranchAddress("k0hoo0", &k0vars[12]);
  tree->SetBranchAddress("k0hoo1", &k0vars[13]);
  tree->SetBranchAddress("k0hoo2", &k0vars[14]);
  tree->SetBranchAddress("k0hoo3", &k0vars[15]);
  tree->SetBranchAddress("k0hoo4", &k0vars[16]);
  tree->SetBranchAddress("k1mm2",  &k1mm2);
  tree->SetBranchAddress("k1et",   &k1vars[0]);
  tree->SetBranchAddress("k1hso00",&k1vars[1]);
  tree->SetBranchAddress("k1hso10",&k1vars[2]);
  tree->SetBranchAddress("k1hso20",&k1vars[3]);
  tree->SetBranchAddress("k1hso02",&k1vars[5]);
  tree->SetBranchAddress("k1hso12",&k1vars[6]);
  tree->SetBranchAddress("k1hso22",&k1vars[7]);
  tree->SetBranchAddress("k1hso04",&k1vars[9]);
  tree->SetBranchAddress("k1hso14",&k1vars[10]);
  tree->SetBranchAddress("k1hso24",&k1vars[11]);
  tree->SetBranchAddress("k1hoo0", &k1vars[12]);
  tree->SetBranchAddress("k1hoo1", &k1vars[13]);
  tree->SetBranchAddress("k1hoo2", &k1vars[14]);
  tree->SetBranchAddress("k1hoo3", &k1vars[15]);
  tree->SetBranchAddress("k1hoo4", &k1vars[16]);
  if(m_second_iter){
    tree->SetBranchAddress("flv",&flv);
    tree->SetBranchAddress("good_icpv",&good_icpv);
    tree->SetBranchAddress("lh0",&lh0);
    tree->SetBranchAddress("lh1",&lh1);
    tree->SetBranchAddress("bdt",&bdt);
  }
  return;
}

void ICPVEvent::SetBranches(TTree* tree){
  tree->Branch("exp",&exp,"exp/I");
  tree->Branch("run",&run,"run/I");
  tree->Branch("evtn",&evtn,"evtn/I");
  tree->Branch("p_d0",&p_d0,"p_d0/D");
  tree->Branch("p_h0",&p_h0,"p_h0/D");
  tree->Branch("p_ks",&p_ks,"p_ks/D");
  tree->Branch("p_pi0_h0",&p_pi0_h0,"p_pi0_h0/D");
  tree->Branch("p_pip_h0",&p_pip_h0,"p_pip_h0/D");
  tree->Branch("p_pim_h0",&p_pim_h0,"p_pim_h0/D");
  tree->Branch("egamma",&egamma,"egamma/D");
  tree->Branch("cos_thr",&cos_thr,"cos_thr/D");
  tree->Branch("cos_hel",&cos_hel,"cos_hel/D");
  tree->Branch("thr_sig",&thr_sig,"thr_sig/D");
  tree->Branch("thr_oth",&thr_oth,"thr_oth/D");
  tree->Branch("phsp",&phsp,"phsp/I");
  tree->Branch("bin",&bin,"bin/I");
  tree->Branch("hel_h0",&hel_h0,"hel_h0/D");
  tree->Branch("hel_pi0",&hel_pi0,"hel_pi0/D");
  tree->Branch("e_g1",&e_g1,"e_g1/D");
  tree->Branch("e_g2",&e_g2,"e_g2/D");
  tree->Branch("e_g3",&e_g3,"e_g3/D");
  tree->Branch("e_g4",&e_g4,"e_g4/D");
  tree->Branch("th_g1",&th_g1,"th_g1/D");
  tree->Branch("th_g2",&th_g2,"th_g2/D");
  tree->Branch("th_g3",&th_g3,"th_g3/D");
  tree->Branch("th_g4",&th_g4,"th_g4/D");
  tree->Branch("r_pip",&r_pip,"r_pip/D");
  tree->Branch("r_pim",&r_pim,"r_pim/D");
  tree->Branch("r_pi1",&r_pi1,"r_pi1/D");
  tree->Branch("r_pi2",&r_pi2,"r_pi2/D");
  tree->Branch("z_pip",&z_pip,"z_pip/D");
  tree->Branch("z_pim",&z_pim,"z_pim/D");
  tree->Branch("z_pi1",&z_pi1,"z_pi1/D");
  tree->Branch("z_pi2",&z_pi2,"z_pi2/D");
  tree->Branch("pt_pip",&pt_pip,"pt_pip/D");
  tree->Branch("pt_pim",&pt_pim,"pt_pim/D");
  tree->Branch("pt_pi1",&pt_pi1,"pt_pi1/D");
  tree->Branch("pt_pi2",&pt_pi2,"pt_pi2/D");
  tree->Branch("px_pim",&px_pim,"px_pim/D");
  tree->Branch("py_pim",&py_pim,"py_pim/D");
  tree->Branch("pz_pim",&pz_pim,"pz_pim/D");
  tree->Branch("px_pip",&px_pip,"px_pip/D");
  tree->Branch("py_pip",&py_pip,"py_pip/D");
  tree->Branch("pz_pip",&pz_pip,"pz_pip/D");
  tree->Branch("px_ks",&px_ks,"px_ks/D");
  tree->Branch("py_ks",&py_ks,"py_ks/D");
  tree->Branch("pz_ks",&pz_ks,"pz_ks/D");
  tree->Branch("chi2_vtx_d0",&chi2_vtx_d0,"chi2_vtx_d0/D");
  tree->Branch("chi2_mass_d0",&chi2_mass_d0,"chi2_mass_d0/D");
  if(m_type){
    tree->Branch("nptag",&nptag,"nptag/I");
    tree->Branch("bin_mc",&bin_mc,"bin_mc/I");
    tree->Branch("flv_mc",&flv_mc,"flv_mc/I");
    tree->Branch("d0ch0",&d0ch0,"d0ch0/I");
    tree->Branch("d0ch1",&d0ch1,"d0ch1/I");
    tree->Branch("d0ch2",&d0ch2,"d0ch2/I");
    tree->Branch("d0ch3",&d0ch3,"d0ch3/I");
    tree->Branch("h0ch0",&h0ch0,"h0ch0/I");
    tree->Branch("h0ch1",&h0ch1,"h0ch1/I");
    tree->Branch("rndm_pi0",&rndm_pi0,"rndm_pi0/I");
    tree->Branch("b0id",&b0id,"b0id/I");
    tree->Branch("b0f",&b0f,"b0f/I");
    tree->Branch("d0id",&d0id,"d0id/I");
    tree->Branch("d0f",&d0f,"d0f/I");
    tree->Branch("h0id",&h0id,"h0id/I");
    tree->Branch("h0f",&h0f,"h0f/I");
    tree->Branch("pi0f",&pi0f,"pi0f/I");
    tree->Branch("dst0id",&dst0id,"dst0id/I");
    tree->Branch("dst0f",&dst0f,"dst0f/I");
    tree->Branch("etapid",&etapid,"etapid/I");
    tree->Branch("etapf",&etapf,"etapf/I");
    if(m_type == 1){
      tree->Branch("mp_mc",&mp_mc,"mp_mc/D");
      tree->Branch("mm_mc",&mm_mc,"mm_mc/D");
      tree->Branch("mp_raw",&mp_raw,"mp_raw/D");
      tree->Branch("mm_raw",&mm_raw,"mm_raw/D");
      tree->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
      tree->Branch("dt_mc",&dt_mc,"dt_mc/D");
      tree->Branch("dz_mc",&dz_mc,"dz_mc/D");
      tree->Branch("dz_mc_sig",&dz_mc_sig,"dz_mc_sig/D");
      tree->Branch("dz_mc_asc",&dz_mc_asc,"dz_mc_asc/D");
      tree->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
      tree->Branch("t_sig_mc",&t_sig_mc,"t_sig_mc/D");
      tree->Branch("z_sig_mc",&z_sig_mc,"z_sig_mc/D");
      tree->Branch("t_asc_mc",&t_asc_mc,"t_asc_mc/D");
      tree->Branch("z_asc_mc",&z_asc_mc,"z_asc_mc/D");
    }
  }
  tree->Branch("mbc",&mbc,"mbc/D");
  tree->Branch("de",&de,"de/D");
  tree->Branch("mp",&mp,"mp/D");
  tree->Branch("mm",&mm,"mm/D");
  tree->Branch("dz",&dz,"dz/D");
  tree->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  tree->Branch("mh0",&mh0);
  tree->Branch("mpi0",&mpi0,"mpi0/D");
  tree->Branch("mk",&mk,"mk/D");
  tree->Branch("md",&md,"md/D");
  tree->Branch("md_raw",&md_raw,"md_raw/D");
  tree->Branch("md_fit",&md_fit,"md_fit/D");
  tree->Branch("mdpip",&mdpip,"mdpip/D");
  tree->Branch("mdpim",&mdpim,"mdpim/D");
  tree->Branch("mdst0",&mdst0,"mdst0/D");
  tree->Branch("dmdst0",&dmdst0,"dmdst0/D");
  tree->Branch("metap",&metap,"metap/D");
  tree->Branch("dmetap",&dmetap,"dmetap/D");
  tree->Branch("mode",&mode,"mode/I");
  tree->Branch("h0mode",&h0mode,"h0mode/I");
  tree->Branch("z_sig",&z_sig,"z_sig/D");
  tree->Branch("z_asc",&z_asc,"z_asc/D");
  tree->Branch("sz_sig",&sz_sig,"sz_sig/D");
  tree->Branch("sz_asc",&sz_asc,"sz_asc/D");
  tree->Branch("ntrk_sig",&ntrk_sig,"ntrk_sig/I");
  tree->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  tree->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  tree->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");
  tree->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  tree->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
  tree->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
  tree->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");
  tree->Branch("h0_chi2",&h0_chi2,"h0_chi2/D");
  tree->Branch("pi0_chi2",&pi0_chi2,"pi0_chi2/D");
  tree->Branch("costhBcms",&costhBcms,"costhBcms/D");
  tree->Branch("flv",&flv,"flv/I");
  tree->Branch("tag_LH",&tag_LH,"tag_LH/D");
  tree->Branch("tag_LH_err",&tag_LH_err,"tag_LH_err/D");
  tree->Branch("k0mm2",  &k0mm2,    "k0mm2/D");
  tree->Branch("k0et",   &k0vars[0],"k0et/D");
  tree->Branch("k0hso00",&k0vars[1],"k0hso00/D");
  tree->Branch("k0hso10",&k0vars[2],"k0hso10/D");
  tree->Branch("k0hso20",&k0vars[3],"k0hso20/D");
  tree->Branch("k0hso02",&k0vars[5],"k0hso02/D");
  tree->Branch("k0hso12",&k0vars[6],"k0hso12/D");
  tree->Branch("k0hso22",&k0vars[7],"k0hso22/D");
  tree->Branch("k0hso04",&k0vars[9],"k0hso04/D");
  tree->Branch("k0hso14",&k0vars[10],"k0hso14/D");
  tree->Branch("k0hso24",&k0vars[11],"k0hso24/D");
  tree->Branch("k0hoo0", &k0vars[12],"k0hoo0/D");
  tree->Branch("k0hoo1", &k0vars[13],"k0hoo1/D");
  tree->Branch("k0hoo2", &k0vars[14],"k0hoo2/D");
  tree->Branch("k0hoo3", &k0vars[15],"k0hoo3/D");
  tree->Branch("k0hoo4", &k0vars[16],"k0hoo4/D");
  tree->Branch("k1mm2",  &k1mm2,"k1mm2/D");
  tree->Branch("k1et",   &k1vars[0],"k1et/D");
  tree->Branch("k1hso00",&k1vars[1],"k1hso00/D");
  tree->Branch("k1hso10",&k1vars[2],"k1hso10/D");
  tree->Branch("k1hso20",&k1vars[3],"k1hso20/D");
  tree->Branch("k1hso02",&k1vars[5],"k1hso02/D");
  tree->Branch("k1hso12",&k1vars[6],"k1hso12/D");
  tree->Branch("k1hso22",&k1vars[7],"k1hso22/D");
  tree->Branch("k1hso04",&k1vars[9],"k1hso04/D");
  tree->Branch("k1hso14",&k1vars[10],"k1hso14/D");
  tree->Branch("k1hso24",&k1vars[11],"k1hso24/D");
  tree->Branch("k1hoo0", &k1vars[12],"k1hoo0/D");
  tree->Branch("k1hoo1", &k1vars[13],"k1hoo1/D");
  tree->Branch("k1hoo2", &k1vars[14],"k1hoo2/D");
  tree->Branch("k1hoo3", &k1vars[15],"k1hoo3/D");
  tree->Branch("k1hoo4", &k1vars[16],"k1hoo4/D");
  tree->Branch("good_icpv",&good_icpv,"good_icpv/I");
  tree->Branch("lh0",&lh0,"lh0/D");
  tree->Branch("lh1",&lh1,"lh1/D");
  tree->Branch("bdt",&bdt,"bdt/D");
  return;
}

int ICPVEvent::FillVectorWithTTree(std::vector<ICPVEvent>& vec, TTree* tree,int type,const bool second_iter,const int svd){
  vec.clear();
  const int NTot = tree->GetEntries();
  ICPVEvent evt(type,second_iter);
  evt.SetBrAddresses(tree);
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    if(svd == 1 && evt.exp>30) continue;
    if(svd == 2 && evt.exp<30) continue;
    vec.push_back(evt);
  }
  return vec.size();
}

ICPVEvent& ICPVEvent::operator=(const ICPVEvent& oth){
  this->exp = oth.exp;
  this->run = oth.run;
  this->evtn = oth.evtn;
  this->p_d0 = oth.p_d0;
  this->p_h0 = oth.p_h0;
  this->p_pi0_h0 = oth.p_pi0_h0;
  this->p_pip_h0 = oth.p_pip_h0;
  this->p_pim_h0 = oth.p_pim_h0;
  this->egamma   = oth.egamma;
  this->cos_thr  = oth.cos_thr;
  this->cos_hel  = oth.cos_hel;
  this->thr_sig  = oth.thr_sig;
  this->thr_oth  = oth.thr_oth;
  this->phsp     = oth.phsp;
  this->bin      = oth.bin;
  this->hel_h0   = oth.hel_h0;
  this->hel_pi0  = oth.hel_pi0;
  this->e_g1   = oth.e_g1;
  this->e_g2   = oth.e_g2;
  this->e_g3   = oth.e_g3;
  this->e_g4   = oth.e_g4;
  this->th_g1   = oth.th_g1;
  this->th_g2   = oth.th_g2;
  this->th_g3   = oth.th_g3;
  this->th_g4   = oth.th_g4;
  this->r_pip   = oth.r_pip;
  this->r_pim   = oth.r_pim;
  this->r_pi1   = oth.r_pi1;
  this->r_pi2   = oth.r_pi2;
  this->z_pip   = oth.z_pip;
  this->z_pim   = oth.z_pim;
  this->z_pi1   = oth.z_pi1;
  this->z_pi2   = oth.z_pi2;
  this->pt_pip   = oth.pt_pip;
  this->pt_pim   = oth.pt_pim;
  this->pt_pi1   = oth.pt_pi1;
  this->pt_pi2   = oth.pt_pi2;
  this->px_pim   = oth.px_pim;
  this->py_pim   = oth.py_pim;
  this->pz_pim   = oth.pz_pim;
  this->px_pip   = oth.px_pip;
  this->py_pip   = oth.py_pip;
  this->pz_pip   = oth.pz_pip;
  this->px_ks   = oth.px_ks;
  this->py_ks   = oth.py_ks;
  this->pz_ks   = oth.pz_ks;
  this->chi2_vtx_d0   = oth.chi2_vtx_d0;
  this->chi2_mass_d0   = oth.chi2_mass_d0;
  this->t_sig_mc   = oth.t_sig_mc;
  this->z_sig_mc   = oth.z_sig_mc;
  this->t_asc_mc   = oth.t_asc_mc;
  this->z_asc_mc   = oth.z_asc_mc;
  this->mp_raw   = oth.mp_raw;
  this->mm_raw   = oth.mm_raw;
  this->mp_mc   = oth.mp_mc;
  this->mm_mc   = oth.mm_mc;
  this->d0_t_mc   = oth.d0_t_mc;
  this->dt_mc   = oth.dt_mc;
  this->dz_mc   = oth.dz_mc;
  this->dz_mc_sig   = oth.dz_mc_sig;
  this->dz_mc_asc   = oth.dz_mc_asc;
  this->b0f   = oth.b0f;
  this->d0f   = oth.d0f;
  this->h0f   = oth.h0f;
  this->pi0f   = oth.pi0f;
  this->bin_mc   = oth.bin_mc;
  this->flv_mc   = oth.flv_mc;
  this->nptag   = oth.nptag;
  this->d0ch0   = oth.d0ch0;
  this->d0ch1   = oth.d0ch1;
  this->d0ch2   = oth.d0ch2;
  this->d0ch3   = oth.d0ch3;
  this->h0ch0   = oth.h0ch0;
  this->h0ch1   = oth.h0ch1;
  this->rndm_pi0   = oth.rndm_pi0;
  this->b0id   = oth.b0id;
  this->d0id   = oth.d0id;
  this->h0id   = oth.h0id;
  this->dst0id   = oth.dst0id;
  this->dst0f   = oth.dst0f;
  this->etapid   = oth.etapid;
  this->etapf   = oth.etapf;
  this->d0_flv_mc   = oth.d0_flv_mc;
  this->mbc   = oth.mbc;
  this->de   = oth.de;
  this->mp   = oth.mp;
  this->mm   = oth.mm;
  this->dz   = oth.dz;
  this->atckpi_max   = oth.atckpi_max;
  this->mpi0   = oth.mpi0;
  this->mh0   = oth.mh0;
  this->mk   = oth.mk;
  this->md   = oth.md;
  this->md_raw   = oth.md_raw;
  this->md_fit   = oth.md_fit;
  this->mdpip   = oth.mdpip;
  this->mdpim   = oth.mdpim;
  this->mdst0   = oth.mdst0;
  this->metap   = oth.metap;
  this->dmdst0   = oth.dmdst0;
  this->dmetap   = oth.dmetap;
  this->mode   = oth.mode;
  this->h0mode   = oth.h0mode;
  this->z_sig   = oth.z_sig;
  this->z_asc   = oth.z_asc;
  this->sz_sig   = oth.sz_sig;
  this->sz_asc   = oth.sz_asc;
  this->ntrk_sig   = oth.ntrk_sig;
  this->ntrk_asc   = oth.ntrk_asc;
  this->ndf_z_sig   = oth.ndf_z_sig;
  this->ndf_z_asc   = oth.ndf_z_asc;
  this->chisq_z_sig   = oth.chisq_z_sig;
  this->chisq_z_asc   = oth.chisq_z_asc;
  this->cl_z_sig   = oth.cl_z_sig;
  this->cl_z_asc   = oth.cl_z_asc;
  this->h0_chi2   = oth.h0_chi2;
  this->pi0_chi2   = oth.pi0_chi2;
  this->costhBcms   = oth.costhBcms;
  this->tag_LH   = oth.tag_LH;
  this->tag_LH_err   = oth.tag_LH_err;
  this->flv   = oth.flv;
  this->good_icpv   = oth.good_icpv;
  this->lh0   = oth.lh0;
  this->lh1   = oth.lh1;
  this->bdt   = oth.bdt;
  for(int i=0; i<8; i++){
    this->d0_chain[i] = oth.d0_chain[i];
    this->h0_chain[i] = oth.h0_chain[i];
  }
  this->k0mm2   = oth.k0mm2;
  this->k1mm2   = oth.k1mm2;
  for(int i=0; i<17; i++){
    this->k0vars[i] = oth.k0vars[i];
    this->k1vars[i] = oth.k1vars[i];
  }
  return *this;
}

void TMVAEvent::Fill(const ICPVEvent& evt){
  m_lh0          = (float)evt.lh0;//mode < 3 ? (float)lh1 : (float)(1./1.005-lh1*lh1);
  m_costhBcms    = (float)abs(evt.costhBcms);
  m_chi2_mass_d0 = (float)log(evt.chi2_mass_d0);
  m_cos_thr      = (float)abs(evt.cos_thr);
  m_thr_sig      = (float)evt.thr_sig;
  m_h0_chi2      = (float)evt.h0_chi2;
  m_egamma       = (float)log(evt.egamma);
  m_p_pi0_h0     = (float)evt.p_pi0_h0;
  m_cos_hel      = (float)abs(evt.cos_hel);
}
