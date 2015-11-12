#ifndef ICPVEVENT_H
#define ICPVEVENT_H

#include <vector>
#include "TTree.h"

class ICPVEvent{
public:
  ICPVEvent(int type,const bool second_iter);
  // type 0 -> Data
  // type 1 -> Sig MC
  // type 2 -> Gen MC
  Int_t exp,run,evtn;
  double p_d0,p_h0,p_ks,p_pi0_h0,p_pip_h0,p_pim_h0,egamma,cos_thr,cos_hel,thr_sig,thr_oth;
  Int_t phsp,bin;
  double hel_h0,hel_pi0;
  double e_g1,e_g2,e_g3,e_g4;
  double th_g1,th_g2,th_g3,th_g4;
  double r_pip, r_pim, r_pi1, r_pi2;
  double z_pip, z_pim, z_pi1, z_pi2;
  double pt_pip, pt_pim, pt_pi1, pt_pi2;
  double px_pim,py_pim,pz_pim;
  double px_pip,py_pip,pz_pip;
  double px_ks,py_ks,pz_ks;
  double chi2_vtx_d0, chi2_mass_d0;
  double t_sig_mc,z_sig_mc,t_asc_mc,z_asc_mc;
  double mp_raw, mm_raw,mp_mc,mm_mc,d0_t_mc,dt_mc,dz_mc;
  double dz_mc_sig, dz_mc_asc;
  Int_t b0f,d0f,h0f,pi0f;
  Int_t bin_mc,flv_mc;
  Int_t nptag;
  int d0ch0,d0ch1,d0ch2,d0ch3;
  int h0ch0,h0ch1;
  int d0_chain[8];
  int h0_chain[8];
  int rndm_pi0;
  Int_t b0id,d0id,h0id,dst0id,dst0f,etapid,etapf;
  Int_t d0_flv_mc;
  double mbc,de,mp,mm,dz,atckpi_max,mpi0,mh0,mk;
  double md, md_raw, md_fit, mdpip, mdpim;
  double mdst0, metap, dmdst0, dmetap;
  Int_t mode,h0mode;
  double z_sig,z_asc;
  double sz_sig,sz_asc;
  Int_t ntrk_sig,ntrk_asc,ndf_z_sig,ndf_z_asc;
  double chisq_z_sig,chisq_z_asc,cl_z_sig,cl_z_asc;
  double h0_chi2,pi0_chi2;
  double costhBcms;
  double tag_LH,tag_LH_err;
  int flv;
  double k0mm2;
  double k0vars[17];
  double k1mm2;
  double k1vars[17];
  Int_t good_icpv;
  double lh0,lh1,bdt;

//  ICPVEvent& operator=(ICPVEvent& othevt);
  void SetBrAddresses(TTree* tree);
  void SetBranches(TTree* tree);
  static int FillVectorWithTTree(std::vector<ICPVEvent>& vec, TTree* tree, int type, const bool second_iter, const int svd = 0);

private:
  int  m_type;
  bool m_second_iter;
};

class TMVAEvent{
  public:
  TMVAEvent(void);
  Float_t m_costhBcms, m_chi2_mass_d0, m_cos_thr, m_thr_sig, m_h0_chi2, m_egamma, m_cos_hel;
  Float_t m_p_pi0_h0, m_p_pip_h0, m_p_pim_h0, m_lh0;

  void Fill(const ICPVEvent& evt);
}

#endif // ICPVEVENT_H
