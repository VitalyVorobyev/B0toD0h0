double EfRkRdetRnp_fullrec(const double x, const int flavor,  const double tau,
                           const double ak, const double ck,
                           const int ntrk_rec, const double sz_rec,
                           const double chisq_z_rec, const int ndf_z_rec,
                           const int ntrk_asc, const double sz_asc,
                           const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
                           const dtres_param_t * const param){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
             &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
             &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
            &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double r_ckak = ck/ak;
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_full_sup(x,
                                               fact_n, ntau_n, fact_p, ntau_p,
                                               fd, fp, tau_np_p, tau_np_n,
                                               mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_full_sup(x,
                                                 fact_n, ntau_n, fact_p, ntau_p,
                                                 fd, fp, tau_np_p_tl, tau_np_n_tl,
                                                 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_full_sup(x,
                                                   fact_n, ntau_n, fact_p, ntau_p,
                                                   fd, fp, tau_np_p, tau_np_n,
                                                   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_full_sup(x,
                                                   fact_n, ntau_n, fact_p, ntau_p,
                                                   fd, fp, tau_np_p_tl, tau_np_n_tl,
                                                   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                         +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdetRnp_full_sup(x,
                                                 fact_n, ntau_n, fact_p, ntau_p,
                                                 fd, fp, tau_np_p, tau_np_n,
                                                 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline void calc_vtxparam_asc(const int ntrk_asc, const double sz_asc,
                              const double chisq_z_asc, const int ndf_z_asc,
                              double*  xi_asc, double* st_asc){
  const double& inv_bgc = dt_resol_global::inv_bgc;
  *xi_asc = (ndf_z_asc > 0 ? chisq_z_asc/ndf_z_asc : 1);
  *st_asc = sz_asc*inv_bgc;
  return;
}

inline void calc_vtxparam_rec(const int ntrk_rec, const double sz_rec,
                              const double chisq_z_rec, const int ndf_z_rec,
                              double*  xi_rec, double* st_rec){
  const double& inv_bgc = dt_resol_global::inv_bgc;;
  *xi_rec = (ndf_z_rec > 0 ? chisq_z_rec/ndf_z_rec : 1);
  *st_rec = sz_rec*inv_bgc;
  return;
}

inline double norm_EfRkRdetRnp_full_sup(const double ll, const double ul,
                                        const double fact_n, const double ntau_n,
                                        const double fact_p, const double ntau_p,
                                        const double fd, const double fp, const double tau_np_p, const double tau_np_n,
                                        const double mu_det, const double sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n*fd;
  double fEp = fact_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  add_EnEn_coef(&fEn, &fEn_np, &fxEn, ntau_n, tau_np_n, fact_n*nfn);
  add_EnEp_coef(&fEn, &fEp_np,        ntau_n, tau_np_p, fact_n*nfp);
  add_EpEn_coef(&fEp, &fEn_np,        ntau_p, tau_np_n, fact_p*nfn);
  add_EpEp_coef(&fEp, &fEp_np, &fxEp, ntau_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(ll, ul, tau_np_n, mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_det, sigma_det);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  return Li;
}

