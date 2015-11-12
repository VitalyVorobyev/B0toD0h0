double EfRdet(const double x, const double tau,
              const int ntrk_rec, const double sz_rec,
              const double chisq_z_rec, const int ndf_z_rec,
              const int ntrk_asc, const double sz_asc,
              const double chisq_z_asc, const int ndf_z_asc,
              const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
             &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
             &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = Ef_conv_gauss(x, tau,  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = Ef_conv_gauss(x, tau,  mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = Ef_conv_gauss(x, tau,  mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = Ef_conv_gauss(x, tau,  mu_tt, sigma_tt);
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
    const double Li_tm    = Ef_conv_gauss(x, tau,  mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

void Rasc_param(const int ntrk_asc,
                const double xi_asc, const double st_asc,
                const dtres_param_t * const param,
                double* ftail_asc,
                double* mu_main_asc, double* Smain_asc,
                double* mu_tail_asc, double* Stail_asc){
  if(ntrk_asc>1){ /*  multiple track vertex  */
    *ftail_asc = (param->ftl_asc_mlt)[0] + (param->ftl_asc_mlt)[1] * xi_asc;
    *Smain_asc = ((param->Sasc)[0] + (param->Sasc)[1] * xi_asc ) * st_asc;
    *Stail_asc = param->Stl_asc_mlt*(*Smain_asc);
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    *ftail_asc = param->ftl_asc;
    *Smain_asc = param->Smn_asc * st_asc;
    *Stail_asc = param->Stl_asc * st_asc;
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
    constraint(Stail_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  *mu_main_asc = 0.0;
  *mu_tail_asc = 0.0;
  return;
}
