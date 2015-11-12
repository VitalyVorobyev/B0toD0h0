#ifndef RKRDETRNPRDF_H
#define RKRDETRNPRDF_H

#include "cnvl.h"
#include "conv_coef.h"
#include "ResConst.h"

class RkRdetRnpPdf : public cnvl, public conv_coef, public ResConst{
public:
  void SetTag(const double& x);
  void SetFlv(const int x);
  RkRdetRnpPdf(const int mc, const int svd, const bool charged = false, const int error = 0);

  double PdfAB(const double& dt,const int ntrk_rec,const double& sz_rec,const double& chisq_z_rec,const int ndf_z_rec,const int ntrk_asc,const double& sz_asc,const double& chisq_z_asc,const int ndf_z_asc, const bool otlr = true, const bool no_interf = false);
  double Pdf(const double& dt,const int ntrk_rec,const double& sz_rec,const double& chisq_z_rec,const int ndf_z_rec,const int ntrk_asc,const double& sz_asc,const double& chisq_z_asc,const int ndf_z_asc, const bool otlr = true, const bool no_interf = false);
  double NoNPPdf(const double& dt, const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc, const bool otlr = true, const bool no_interf = false);
  double PdfRrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec);
  double PdfRascRnp(const double& x, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc);

  double Rrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec);
  double norm_Rrec(const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec);
  double RascRnp(const double& x, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc);
  double norm_RascRnp(const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc);

  double EfRkRdetRnp_fullrec(const double& x,
                             const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);

  double norm_EfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);

  double MfRkRdetRnp_fullrec(const double& x,
                             const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);

  double norm_MfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);

  double AfRkRdetRnp_fullrec(const double& x,
                             const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);

  double norm_AfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                             const double& chisq_z_rec, const int ndf_z_rec,
                             const int ntrk_asc, const double& sz_asc,
                             const double& chisq_z_asc, const int ndf_z_asc);
// * RkRdet *
  double EfRkRdet_fullrec(const double& x,
                          const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

  double AfRkRdet_fullrec(const double& x,
                          const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

  double MfRkRdet_fullrec(const double& x,
                          const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

  double norm_EfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

  double norm_AfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

  double norm_MfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,
                          const double& chisq_z_rec, const int ndf_z_rec,
                          const int ntrk_asc, const double& sz_asc,
                          const double& chisq_z_asc, const int ndf_z_asc);

private:
  double EfRkRdetRnp_full_sup(const double& x,
                              const double& fact_n, const double& ntau_n,
                              const double& fact_p, const double& ntau_p,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double norm_EfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n,
                              const double& fact_p, const double& ntau_p,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double MfRkRdetRnp_full_sup(const double& x,
                              const double& fact_n, const double& ntau_n, const double& ndm_n,
                              const double& fact_p, const double& ntau_p, const double& ndm_p,
                              const double& ndmtau, const double& cktau,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double norm_MfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n,
                              const double& fact_p, const double& ntau_p, const double& ndm_p,
                              const double& ndmtau, const double& cktau,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double AfRkRdetRnp_full_sup(const double& x,
                              const double& fact_n, const double& ntau_n, const double& ndm_n,
                              const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double norm_AfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n,
                              const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau,
                              const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                              const double& mu_det, const double& sigma_det);

  double EfRkRdet_full_sup(const double& x,const double& fact_n, const double& ntau_n, const double& fact_p, const double& ntau_p,const double& mu_det, const double& sigma_det);
  double AfRkRdet_full_sup(const double& x, const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& mu, const double& sigma);
  double MfRkRdet_full_sup(const double& x,const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p,const double& ndmtau, const double& cktau,const double& mu, const double& sigma);
  double norm_EfRkRdet_full_sup(const double& fact_n, const double& ntau_n,const double& fact_p, const double& ntau_p, const double& mu_det, const double& sigma_det);
  double norm_AfRkRdet_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& mu, const double& sigma);
  double norm_MfRkRdet_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& cktau, const double& mu, const double& sigma);

  double sum_sigma(const double& s1, const double& s2);
  void constraint(double& x, const double& ll, const double& ul);
  void constraint(double& x, const double& ll);

  void Rasc_param(const int ntrk_asc,const double& xi_asc,const double& st_asc,double& ftail_asc,double& mu_main_asc, double& Smain_asc,double& mu_tail_asc, double& Stail_asc);
  void Rrec_param(const int ntrk_rec,const double& xi_rec,const double& st_rec,double& ftail_rec,double& mu_main_rec, double& Smain_rec,double& mu_tail_rec, double& Stail_rec);
  void Rnp_param(const int ntrk_asc,const double& Smain_asc,const double& Stail_asc,double& fd, double& fp,double& tau_np_p, double& tau_np_n,double& tau_np_p_tl, double& tau_np_n_tl);
  void Rnp_param(const int ntrk_asc, const double& xi_asc, const double& st_asc,const double& Smain_asc, const double& Stail_asc,double& fd, double& fp,double& tau_np_p, double& tau_np_n,double& tau_np_p_tl, double& tau_np_n_tl);
  void Rnp_param_03(const int ntrk_asc,const double& xi_asc, const double& st_asc,const double& Smain_asc, const double& Stail_asc,double& fd, double& fp,double& tau_np_p, double& tau_np_n,double& tau_np_p_tl, double& tau_np_n_tl);
  void Rnp_param_10(const int ntrk_asc,const double& xi_asc, const double& st_asc,double& fd, double& fp,double& tau_np_p, double& tau_np_n,double& tau_np_p_tl, double& tau_np_n_tl);

  void swap_rnp_param(double& fp,double& tau_np_p, double& tau_np_n,double& tau_np_p_tl, double& tau_np_n_tl);
  void calc_vtxparam_asc(const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc,double& xi_asc, double& st_asc);
  void calc_vtxparam_rec(const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec,double& xi_rec, double& st_rec);

  double AddOutlier(const double& x,const double& Lin, const int ntrk_rec, const int ntrk_asc, const double& nLi = 1.0, const double& alpha = 1.0);
  double Add_Outlier(const double& x, const double Lin,const int ntrk_rec, const int ntrk_asc, const double& nLi, const double& alpha);
private:
  int flavor;
  int keeptagl;
  int m_svd;
  double cexp, amix;

  double cm2ps;

};

#endif // RKRDETRNPRDF_H
