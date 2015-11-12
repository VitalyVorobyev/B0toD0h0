double Rrec(const double x,
            const int ntrk_rec, const double sz_rec,
            const double chisq_z_rec, const int ndf_z_rec,
            const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
             &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  const double Li_mn = gaussian(x, mu_main_rec, Smain_rec);
  if(ftail_rec>0.0){ /* Mainly single track case */
    const double Li_tl = gaussian(x, mu_tail_rec, Stail_rec);
    const double Li    = (1.0-ftail_rec)*Li_mn + ftail_rec*Li_tl;
    return Li;
  }
  return Li_mn;
}

inline void calc_vtxparam_rec(const int ntrk_rec, const double sz_rec,
                              const double chisq_z_rec, const int ndf_z_rec,
                              double*  xi_rec, double* st_rec){
  const double& inv_bgc = dt_resol_global::inv_bgc;;
  *xi_rec = (ndf_z_rec > 0 ? chisq_z_rec/ndf_z_rec : 1);
  *st_rec = sz_rec*inv_bgc;
  return;
}

void Rrec_param(const int ntrk_rec,
                const double xi_rec, const double st_rec,
                const dtres_param_t * const param,
                double* ftail_rec,
                double* mu_main_rec, double* Smain_rec,
                double* mu_tail_rec, double* Stail_rec){
  if(ntrk_rec>1){ /*  multiple track vertex  */
    *ftail_rec = (param->ftl_rec_mlt)[0] + (param->ftl_rec_mlt)[1] * xi_rec;
    *Smain_rec = ((param->Srec)[0] + (param->Srec)[1] * xi_rec ) * st_rec;
    *Stail_rec = param->Stl_rec_mlt*(*Smain_rec);
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    *ftail_rec = param->ftl_rec;
    *Smain_rec = param->Smn_rec * st_rec;
    *Stail_rec = param->Stl_rec * st_rec;
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
    constraint(Stail_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  *mu_main_rec = 0.0;
  *mu_tail_rec = 0.0;
  return;
}
