#ifndef MYPARAMS_H
#define MYPARAMS_H

#include <string>

using namespace std;

class MyParams{
public:
  MyParams(void);

  int bin(const int j);
  int flv_ind(const int flv);
  int bin_ind(const int bin);
  int flv(const int k);
  int q(const double& tag_LH);
  double K(const int bin);
  //inline double N(const int bin, const int flv);
  double N(const int bin, const int flv, const double& wt);
  int get_wbin(const double& tag);
  double get_wtag_prob(const double& wt, const int exp,const bool data);

  string GetLabel(const int mode,const int h0mode);
  double get_argedge(const int mode,const int h0mode);
  double f_p_f_bbc(const int mode,const int h0mode);

  double get_de_fit_min(void)  const {return de_fit_min;}
  double get_de_fit_max(void)  const {return de_fit_max;}
  double get_mbc_fit_min(void) const {return mbc_fit_min;}
  double get_mbc_fit_max(void) const {return mbc_fit_max;}
  double get_mbc_SB_min(void)  const {return mbc_fit_min;}
  double get_mbc_SB_max(void)  const {return mbc_side;}

  double get_de_min_h0(const int mode,const int h0mode) const;
  double get_de_max_h0(const int mode,const int h0mode) const;
  double get_mbc_min_h0(const int mode,const int h0mode) const;
  double get_mbc_max_h0(const int mode,const int h0mode) const;
  double get_mh0(const int mode) const;
  double get_mh0_min(const int mode,const int h0mode) const;
  double get_mh0_max(const int mode,const int h0mode) const;
  double get_dm_etap_min(const int h0mode) const;
  double get_dm_etap_max(const int h0mode) const;

  double get_md(void) const {return DMass;}
  double get_md_min(void) const {return md_min;}
  double get_md_max(void) const {return md_max;}
  double get_mk(void) const {return KMass;}
  double get_mk_min(void) const {return mk_min;}
  double get_mk_max(void) const {return mk_max;}
  double get_mpi0(void) const {return Pi0Mass;}
  double get_mpi0_min(void) const {return mpi0_min;}
  double get_mpi0_max(void) const {return mpi0_max;}
  double get_meta(void) const {return EtaMass;}
  double get_metagg_min(void) const {return metagg_min;}
  double get_metagg_max(void) const {return metagg_max;}
  double get_metappp_min(void) const {return metappp_min;}
  double get_metappp_max(void) const {return metappp_max;}
  double get_dm_etagg_min(void) const {return dmetapgg_min;}
  double get_dm_etagg_max(void) const {return dmetapgg_max;}
  double get_dm_etappp_min(void) const {return dmetapppp_min;}
  double get_dm_etappp_max(void) const {return dmetapppp_max;}
  double get_momega(void) const {return OmegaMass;}
  double get_momega_min(void) const {return momega_min;}
  double get_momega_max(void) const {return momega_max;}
  double get_dm_dst0_min(void) const {return dmdst0_min;}
  double get_dm_dst0_max(void) const {return dmdst0_max;}

  double get_atckpi_cut(void) const {return atckpi_cut;}
  double bdt_cut(const int mode,const int h0mode=10) const;
  double lh0_cut(const int mode,const int h0mode=10) const;

  double get_de0_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_s_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_de0_CBl_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_de0_CBr_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_s_CBl_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_s_CBr_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_alphal_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_alphar_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_f_CBl_de_sig(const int mode,const int h0mode, const int b0f=1) const;
  double get_f_CBr_de_sig(const int mode,const int h0mode, const int b0f=1) const;

  double get_a_s_mbc_sig(const int mode,const int h0mode) const;
  double get_b_s_mbc_sig(const int mode,const int h0mode) const;
  double get_c_s_mbc_sig(const int mode,const int h0mode) const;

  double get_a_mbc0_sig(const int mode, const int h0mode) const;
  double get_b_mbc0_sig(const int mode, const int h0mode) const;
  double get_c_mbc0_sig(const int mode, const int h0mode) const;

  double get_a_mbc0_sig_tail(const int mode) const;
  double get_b_mbc0_sig_tail(const int mode) const;
  double get_c_mbc0_sig_tail(const int mode) const;

  double get_a_s_mbc_sig_tail(const int mode) const;
  double get_b_s_mbc_sig_tail(const int mode) const;
  double get_c_s_mbc_sig_tail(const int mode) const;

  // Combinatorial //
  double get_a_c1_bb_cmb(const int mode,const int h0mode);
  double get_b_c1_bb_cmb(const int mode,const int h0mode);
  double get_c2_bb_cmb(const int mode,const int h0mode);
  double get_c1_qq_cmb(const int mode,const int h0mode);
  double get_c2_qq_cmb(const int mode,const int h0mode);
  double get_argpar_cmb_bb(const int mode,const int h0mode);
  double get_mbc0_cmb_bb(const int mode,const int h0mode);
  double get_s_mbc_cmb_bb(const int mode,const int h0mode);
  double get_fg_cmb_bb(const int mode,const int h0mode);
  double get_argpar_cmb_qq(const int mode,const int h0mode);

  // Partial //
  double get_de0_part(const int mode,const int h0mode);
  double get_slopel_part(const int mode,const int h0mode);
  double get_sloper_part(const int mode,const int h0mode);
  double get_steep_part(const int mode,const int h0mode);
  double get_p5_part(const int mode,const int h0mode);
  double get_b_s_mbc_part(const int mode,const int h0mode);
  double get_k_s_mbc_part(const int mode);
  double get_b_mbc0_part(const int mode,const int h0mode);
  double get_k_mbc0_part(const int mode);
  double get_argpar_part_bb(const int mode,const int h0mode);
  double get_mbc0_part(const int mode,const int h0mode);
  double get_s_mbc_part(const int mode,const int h0mode);
  double get_fg_mbc_part(const int mode,const int h0mode);

  bool IsInEllips(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde);
  double EllipsR2(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde);

  double h0mass(const int mode);

  double get_f_cont_in_bkg_sig(const int mode, const int h0mode);
  double get_f_cont_in_bkg_sideband(const int mode, const int h0mode);

private:
  double de_fit_min,  de_fit_max;
  double mbc_fit_min, mbc_fit_max, mbc_side;

  double de_min_pi0,     de_max_pi0;
  double de_min_etagg,   de_max_etagg;
  double de_min_etappp,  de_max_etappp;
  double de_min_etapgg,  de_max_etapgg;
  double de_min_etapppp, de_max_etapppp;
  double de_min_omega,   de_max_omega;
  double de_min_dst0pi0, de_max_dst0pi0;
  double de_min_dst0etagg,  de_max_dst0etagg;
  double de_min_dst0etappp, de_max_dst0etappp;

  double mbc_min_pi0,     mbc_max_pi0;
  double mbc_min_etagg,   mbc_max_etagg;
  double mbc_min_etappp,  mbc_max_etappp;
  double mbc_min_etapgg,  mbc_max_etapgg;
  double mbc_min_etapppp, mbc_max_etapppp;
  double mbc_min_omega,   mbc_max_omega;
  double mbc_min_dst0pi0, mbc_max_dst0pi0;
  double mbc_min_dst0etagg,  mbc_max_dst0etagg;
  double mbc_min_dst0etappp, mbc_max_dst0etappp;

  double DMass, md_min, md_max;
  double KMass, mk_min, mk_max;
  double Pi0Mass, mpi0_min, mpi0_max;
  double EtaMass, metagg_min, metagg_max, metappp_min, metappp_max;
  double dmetapgg_min, dmetapgg_max, dmetapppp_min, dmetapppp_max;
  double OmegaMass, momega_min, momega_max;
  double dmdst0_min, dmdst0_max;

  double atckpi_cut;

  double bdt_cut_pi0, bdt_cut_etagg, bdt_cut_etappp, bdt_cut_omega;
  double lh0_cut_etap, lh0_cut_dst0pi0, lh0_cut_dst0eta;

  // * dE-Mbc Signal * //
  // pi0 dE-Mbc shape //
  // mbc
  double m_a_mbc0_pi1m;
  double m_b_mbc0_pi1m;
  double m_c_mbc0_pi1m;
  double m_a_s_pi1m;
  double m_b_s_pi1m;
  double m_c_s_pi1m;
  double m_alpha_nks_pi1m;

  // de
  double m_nl_pi1m;
  double m_alphal_pi1m;
  double m_nr_pi1m;
  double m_alphar_pi1m;
  double m_de0_pi1m;
  double m_deCBl_pi1m;
  double m_deCBr_pi1m;
  double m_fCBl_pi1m;
  double m_fCBr_pi1m;
  double m_s1_pi1m;
  double m_sCBl_pi1m;
  double m_sCBr_pi1m;

  // eta -> gg dE-Mbc shape //
  // mbc
  double m_a_mbc0_eta10m1;
  double m_b_mbc0_eta10m1;
  double m_c_mbc0_eta10m1;
  double m_a_s_eta10m1;
  double m_b_s_eta10m1;
  double m_c_s_eta10m1;
  double m_alpha_nks_eta10m1;

  // de
  double m_alphal_eta10m1;
  double m_alphar_eta10m1;
  double m_de0_eta10m1;
  double m_deCBl_eta10m1;
  double m_deCBr_eta10m1;
  double m_fCBl_eta10m1;
  double m_fCBr_eta10m1;
  double m_nl_eta10m1;
  double m_nr_eta10m1;
  double m_s1_eta10m1;
  double m_sCBl_eta10m1;
  double m_sCBr_eta10m1;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // mbc
  double m_c0_eta201;
  double m_c1_eta201;
  double m_c2_eta201;
  double m_a_s_eta201;
  double m_b_s_eta201;
  double m_c_s_eta201;
  double m_alpha_nks_eta201;

  double m_a_mbc0_eta205;
  double m_b_mbc0_eta205;
  double m_c_mbc0_eta205;
  double m_a_s_eta205;
  double m_b_s_eta205;
  double m_c_s_eta205;
  double m_alpha_nks_eta205;

  // de
  double m_alphal_eta201;
  double m_alphar_eta201;
  double m_de0_eta201;
  double m_deCBl_eta201;
  double m_deCBr_eta201;
  double m_fCBl_eta201;
  double m_fCBr_eta201;
  double m_nl_eta201;
  double m_nr_eta201;
  double m_s1_eta201;
  double m_sCBl_eta201;
  double m_sCBr_eta201;

  double m_alphal_eta205;
  double m_alphar_eta205;
  double m_de0_eta205;
  double m_deCBl_eta205;
  double m_deCBr_eta205;
  double m_fCBl_eta205;
  double m_fCBr_eta205;
  double m_nl_eta205;
  double m_nr_eta205;
  double m_s1_eta205;
  double m_sCBl_eta205;
  double m_sCBr_eta205;

  // omega dE-Mbc shape //
  // mbc
  double m_mbc0_omega1;
  double m_c0_omega1;
  double m_c1_omega1;
  double m_c2_omega1;
  double m_a_s_omega1;
  double m_b_s_omega1;
  double m_c_s_omega1;
  double m_alpha_nks_omega1;

  double m_a_mbc0_omega5;
  double m_b_mbc0_omega5;
  double m_c_mbc0_omega5;
  double m_a_s_omega5;
  double m_b_s_omega5;
  double m_c_s_omega5;
  double m_alpha_nks_omega5;

  // de
  double m_alphal_omega201;
  double m_alphar_omega201;
  double m_de0_omega201;
  double m_deCBl_omega201;
  double m_deCBr_omega201;
  double m_fCBl_omega201;
  double m_fCBr_omega201;
  double m_nl_omega201;
  double m_nr_omega201;
  double m_s1_omega201;
  double m_sCBl_omega201;
  double m_sCBr_omega201;

  double m_alphal_omega205;
  double m_alphar_omega205;
  double m_de0_omega205;
  double m_deCBl_omega205;
  double m_deCBr_omega205;
  double m_fCBl_omega205;
  double m_fCBr_omega205;
  double m_nl_omega205;
  double m_nr_omega205;
  double m_s1_omega205;
  double m_sCBl_omega205;
  double m_sCBr_omega205;

  // * dE-Mbc Combinatorial * //
  // pi0 dE-Mbc shape //
  double m_cmb_argpar_qq_pi0;
  double m_cmb_c1_pi0;
  double m_cmb_c2_pi0;
  double m_cmb_argpar_bb_pi0;
  double m_cmb_mbc0_bb_pi0;
  double m_cmb_mbcw_bb_pi0;
  double m_cmb_f_g_bb_pi0;
  double m_cmb_c10_pi0;
  double m_cmb_c11_pi0;
  double m_cmb_c20_pi0;
  double m_cmb_fbb_pi0;

  // eta -> gg dE-Mbc shape //
  // gg QQ combinatorial
  double m_cmb_argpar_qq_eta10;
  double m_cmb_c1_eta10;
  double m_cmb_c2_eta10;
  // gg bb  combinatorial
  double m_cmb_argpar_bb_eta10;
  double m_cmb_mbc0_bb_eta10;
  double m_cmb_mbcw_bb_eta10;
  double m_cmb_f_g_bb_eta10;
  double m_cmb_c10_eta10;
  double m_cmb_c11_eta10;
  double m_cmb_c20_eta10;
  double m_cmb_fbb_eta10;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // ppp QQ combinatorial
  double m_cmb_argpar_qq_eta20;
  double m_cmb_c1_eta20;
  double m_cmb_c2_eta20;
  // ppp bb  combinatorial
  double m_cmb_argpar_bb_eta20;
  double m_cmb_mbc0_bb_eta20;
  double m_cmb_mbcw_bb_eta20;
  double m_cmb_f_g_bb_eta20;
  double m_cmb_c10_eta20;
  double m_cmb_c11_eta20;
  double m_cmb_c20_eta20;
  double m_cmb_fbb_eta20;

  // omega dE-Mbc shape //
  // QQ combinatorial
  double m_cmb_argpar_qq_omega;
  double m_cmb_c1_omega;
  double m_cmb_c2_omega;
  // bb combinatorial
  double m_cmb_argpar_bb_omega;
  double m_cmb_mbc0_bb_omega;
  double m_cmb_mbcw_bb_omega;
  double m_cmb_f_g_bb_omega;
  double m_cmb_c10_omega;
  double m_cmb_c11_omega;
  double m_cmb_c20_omega;

  // * dE-Mbc Partial BB * //
  // pi0 dE-Mbc shape //
  // Mbc
  double m_peak_b_mbc0_pi0;
  double m_peak_k_mbc0_pi0;
  double m_peak_b_s_pi0;
  double m_peak_k_s_pi0;

  // dE
  double m_peak_de0r_pi0;
  double m_peak_slopel_pi0;
  double m_peak_sloper_pi0;
  double m_peak_steep_pi0;
  double m_peak_p5_pi0;

  // eta -> gg dE-Mbc shape //
  // Mbc
  double m_peak_b_mbc0_eta10;
  double m_peak_k_mbc0_eta10;
  double m_peak_b_s_eta10;
  double m_peak_k_s_eta10;
  // dE
  double m_peak_de0r_eta10;
  double m_peak_slopel_eta10;
  double m_peak_sloper_eta10;
  double m_peak_steep_eta10;
  double m_peak_p5_eta10;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // Mbc
  double m_peak_mbc0_eta20;
  double m_peak_s_eta20;
  double m_peak_fg_eta20;
  // dE
  double m_peak_de0r_eta20;
  double m_peak_slopel_eta20;
  double m_peak_sloper_eta20;
  double m_peak_steep_eta20;
  double m_peak_p5_eta20;

  // omega dE-Mbc shape //
  // Mbc
  double m_peak_mbc0_omega;
  double m_peak_s_omega;
  double m_peak_fg_omega;
  // dE
  double m_peak_de0r_omega;
  double m_peak_slopel_omega;
  double m_peak_sloper_omega;
  double m_peak_steep_omega;
  double m_peak_p5_omega;
};

#endif // MYPARAMS_H
