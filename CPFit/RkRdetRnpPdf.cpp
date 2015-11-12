#include "RkRdetRnpPdf.h"

#ifndef DTRES_EXTERAM_THRE
#define DTRES_EXTERAM_THRE -FLT_MAX
#endif

// *** Omega *** //
// Omega SVD1:
//  w[1] = 0.480056 +- 0.00453218 (0.5)
//  w[2] = 0.417534 +- 0.00511032 (0.420827)
//  w[3] = 0.292795 +- 0.00465175 (0.300296)
//  w[4] = 0.215926 +- 0.00524856 (0.219317)
//  w[5] = 0.149944 +- 0.00441952 (0.154636)
//  w[6] = 0.103546 +- 0.0038593 (0.0916131)
//  w[7] = 0.0185879 +- 0.00135597 (0.0228891)
//dw[1] = -0.0552954 +- 0.0033853 (0)
//dw[2] = 0.0486283 +- 0.0039842 (0.0583019)
//dw[3] = 0.00770271 +- 0.00331092 (0.00573998)
//dw[4] = -0.00286304 +- 0.00355698 (-0.0392635)
//dw[5] = -0.00465035 +- 0.00346327 (0.00474508)
//dw[6] = 0.00589101 +- 0.00275603 (-0.0118737)
//dw[7] = -0.00420355 +- 0.00109119 (-0.00585326)

// Omega SVD2:
//  w[1] = 0.469111 +- 0.00169265 (0.5)
//  w[2] = 0.411263 +- 0.0019921 (0.412222)
//  w[3] = 0.304424 +- 0.00165546 (0.307838)
//  w[4] = 0.209285 +- 0.00177849 (0.212765)
//  w[5] = 0.149856 +- 0.00173164 (0.149933)
//  w[6] = 0.0860236 +- 0.00137802 (0.0913264)
//  w[7] = 0.0233198 +- 0.000545593 (0.0218754)
//dw[1] = 0.0244654 +- 0 (0)
//dw[2] = 0.0125918 +- 0 (0.00408778)
//dw[3] = 0.0138589 +- 0 (0.010326)
//dw[4] = -0.0166051 +- 0 (-0.00479522)
//dw[5] = -0.00348628 +- 0 (0.00151989)
//dw[6] = 0.00554697 +- 0 (0.0143633)
//dw[7] = -0.00263454 +- 0 (0.00189979)

//const double w_mc_svd1[7]    = {0.480056,0.417534,0.292795,0.215926,0.149944,0.103546 ,0.0185879};
//const double w_mc_svd2[7]    = {0.469111,0.411263,0.304424,0.209285,0.149856,0.0860236,0.0233198};
//const double dw_mc_svd1[7]   = {-0.0552954,0.0486283,0.00770271,-0.00286304,-0.00465035,0.00589101,-0.00420355};
//const double dw_mc_svd2[7]   = { 0.0244654,0.0125918,0.0138589,-0.0166051,-0.00348628,0.00554697,-0.00263454};
// ***  *** //

// *** pi0 *** //
//SVD1:
//  w[1] = 0.474916 +- 0.00364659 (0.5)
//  w[2] = 0.408793 +- 0.00409292 (0.420827)
//  w[3] = 0.305029 +- 0.00377445 (0.300296)
//  w[4] = 0.222252 +- 0.00428229 (0.219317)
//  w[5] = 0.155031 +- 0.00353754 (0.154636)
//  w[6] = 0.080184 +- 0.00279251 (0.0916131)
//  w[7] = 0.021744 +- 0.00116426 (0.0228891)

// dw[1] = 0.0409622 +- 0.00290129 (0)
// dw[2] = -0.0542263 +- 0.00343594 (0.0583019)
// dw[3] = 0.00754399 +- 0.00284783 (0.00573998)
// dw[4] = 0.0245712 +- 0.00305279 (-0.0392635)
// dw[5] = -0.0138642 +- 0.00297483 (0.00474508)
// dw[6] = 0.00111057 +- 0.00240626 (-0.0118737)
// dw[7] = 0.00442319 +- 0.000897826 (-0.00585326)

//SVD2:
//  w[1] = 0.4712 +- 0.00145065 (0.5)
//  w[2] = 0.410327 +- 0.00171797 (0.412222)
//  w[3] = 0.303792 +- 0.00142392 (0.307838)
//  w[4] = 0.209676 +- 0.00152639 (0.212765)
//  w[5] = 0.149952 +- 0.00148742 (0.149933)
//  w[6] = 0.0899141 +- 0.00120313 (0.0913264)
//  w[7] = 0.0214378 +- 0.000448913 (0.0218754)

// dw[1] = -0.028217 +- 0 (0)
// dw[2] = -0.00790949 +- 0 (0.00408778)
// dw[3] = -0.0078825 +- 0 (0.010326)
// dw[4] = 0.0163752 +- 0 (-0.00479522)
// dw[5] = 0.00584695 +- 0 (0.00151989)
// dw[6] = -0.0037633 +- 0 (0.0143633)
// dw[7] = -0.00361328 +- 0 (0.00189979)
// ***  *** //

const double wbins[8]        = {0.0,0.1,0.25,0.5,0.625,0.75,0.875,1.01};
const double w_mc_svd1[7]    = {0.5,0.420827,0.300296,0.219317,0.154636,0.0916131,0.0228891};
const double w_data_svd1[7]  = {0.5,0.418852,0.329879,0.233898,0.170608,0.099791, 0.0228501};
const double w_mc_svd2[7]    = {0.5,0.412222,0.307838,0.212765,0.149933,0.0913264,0.0218754};
const double w_data_svd2[7]  = {0.5,0.418826,0.319303,0.222948,0.163191,0.104085, 0.0251454};

const double dw_mc_svd1[7]   = {0.0, 0.0583019, 0.00573998,-0.0392635, 0.00474508,-0.0118737, -0.00585326};
const double dw_data_svd1[7] = {0.0, 0.0569661, 0.0126192, -0.0147724,-0.000550289,0.00887704, 0.00465683};
const double dw_mc_svd2[7]   = {0.0, 0.00408778,0.010326,  -0.00479522,0.00151989, 0.0143633,  0.00189979};
const double dw_data_svd2[7] = {0.0,-0.00877001,0.0103515, -0.0109253,-0.0186365,  0.00168037,-0.0036441};

const double w_data_svd1_posi[7]   = {0.,7.235697e-03,7.129388e-03,7.417778e-03,6.885875e-03,6.761047e-03,4.336734e-03};
const double w_data_svd1_nega[7]   = {0.,6.001569e-03,6.430566e-03,7.693083e-03,6.416449e-03,8.807757e-03,4.587614e-03};
const double w_data_svd2_posi[7]   = {0.,4.152612e-03,3.243236e-03,3.721417e-03,3.315138e-03,3.180302e-03,2.175087e-03};
const double w_data_svd2_nega[7]   = {0.,3.577812e-03,2.803811e-03,3.486607e-03,4.241595e-03,3.696399e-03,3.077622e-03};

const double dw_data_svd1_posi[7]  = {0.,3.961548e-03,3.543129e-03,4.129422e-03,4.169570e-03,3.998982e-03,2.433324e-03};
const double dw_data_svd1_nega[7]  = {0.,3.927049e-03,3.698619e-03,4.179366e-03,4.602366e-03,3.914627e-03,2.360543e-03};

double w[7],dw[7];

void RkRdetRnpPdf::SetFlv(const int x){
  flv = x;
  amix = flv;
  cexp = 1;
}

void RkRdetRnpPdf::SetTag(const double& x){
  if(abs(x)>1.0001){
    cout << x << " --- tag>1 !!!" << endl;
    return;
  }
  flv = x>0 ? 1 : -1;
  int i=0;
  while(wbins[i+1] < abs(x)) i++;
  amix = flv*(1-2*w[i]);
  cexp = 1;//-flv*dw[i];
//  cout << "SetTag:  tag = " << x << ", bin = " << i << ", w = " << w[i] << ", dw = " << dw[i] << ", flv = " << flv << ", amix = " << amix << ", cexp = " << cexp << endl;
  return;
}

RkRdetRnpPdf::RkRdetRnpPdf(const int mc, const int svd, const bool charged, const int error): cnvl(), ResConst(mc,svd,charged,error) {
  cm2ps = 78.48566945838871754705;
  flavor = 0;
  keeptagl = 0;
  m_svd = svd;
  cexp = 1; amix = 1;
  if(mc == 1){
    if(svd == 1){
      for(int i=0; i<7; i++){w[i] = w_mc_svd1[i]; dw[i] = dw_mc_svd1[i];}
    }
    else{
      for(int i=0; i<7; i++){w[i] = w_mc_svd2[i]; dw[i] = dw_mc_svd2[i];}
    }
  } else{
    if(svd == 1){
      for(int i=0; i<7; i++){w[i] = w_data_svd1[i]; dw[i] = dw_data_svd1[i];}
    }
    else{
      for(int i=0; i<7; i++){w[i] = w_data_svd2[i]; dw[i] = dw_data_svd2[i];}
    }
  }
}

double RkRdetRnpPdf::sum_sigma(const double& s1, const double& s2){
  return sqrt(s1*s1+s2*s2);
}

void RkRdetRnpPdf::constraint(double& x, const double& ll, const double& ul){
  if(x<ll){
    x = ll;
  } else if(x>ul){
    x = ul;
  }
  return;
}

void RkRdetRnpPdf::constraint(double& x, const double& ll){
  if(x<ll) x = ll;
  return;
}

void RkRdetRnpPdf::Rasc_param(const int ntrk_asc,
                const double& xi_asc, const double& st_asc,
                double& ftail_asc,
                double& mu_main_asc, double& Smain_asc,
                double& mu_tail_asc, double& Stail_asc){
  if(ntrk_asc>1){ /*  multiple track vertex  */
    ftail_asc = ftl_asc_mlt[0] + ftl_asc_mlt[1]*xi_asc;
    Smain_asc = (Sasc[0] + Sasc[1]*xi_asc)*st_asc;
    Stail_asc = Stl_asc_mlt*Smain_asc;
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    ftail_asc = ftl_asc;
    Smain_asc = Smn_asc*st_asc;
    Stail_asc = Stl_asc*st_asc;
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
    constraint(Stail_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  mu_main_asc = 0.0;
  mu_tail_asc = 0.0;
  return;
}

void RkRdetRnpPdf::Rrec_param(const int ntrk_rec, const double& xi_rec, const double& st_rec,
                double& ftail_rec,
                double& mu_main_rec, double& Smain_rec,
                double& mu_tail_rec, double& Stail_rec){
  if(ntrk_rec>1){ /*  multiple track vertex  */
    ftail_rec = ftl_rec_mlt[0] + ftl_rec_mlt[1]*xi_rec;
    Smain_rec = (Srec[0] + Srec[1]*xi_rec)*st_rec;
    Stail_rec = Stl_rec_mlt*Smain_rec;
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    ftail_rec = ftl_rec;
    Smain_rec = Smn_rec*st_rec;
    Stail_rec = Stl_rec*st_rec;
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
    constraint(Stail_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  mu_main_rec = 0.0;
  mu_tail_rec = 0.0;
  return;
}

void RkRdetRnpPdf::Rnp_param(const int ntrk_asc,
               const double& Smain_asc, const double& Stail_asc,
               double& fd, double& fp,
               double& tau_np_p, double& tau_np_n,
               double& tau_np_p_tl, double& tau_np_n_tl){
//  const int flvidx = (flavor==0) ? 0 : 1;
//  const int keeptlidx = (keeptagl==0)? 0 : 1;
  if(ntrk_asc>1){ /*  multiple track vertex  */
    fd = fd_np_mlt[keeptagl];
    fp = fp_np_mlt;
    tau_np_p = tau_np_p_mlt[0] + tau_np_p_mlt[1]*Smain_asc;
    tau_np_n = tau_np_n_mlt[0] + tau_np_n_mlt[1]*Smain_asc;
    tau_np_p_tl = 0.0;
    tau_np_n_tl = 0.0;
  }else{ /*  single track vertex */
    fd = fd_np_sgl[keeptagl];
    fp = fp_np_sgl;
    tau_np_p    = tau_np_p_sgl[0] + tau_np_p_sgl[1]*Smain_asc;
    tau_np_n    = tau_np_n_sgl[0] + tau_np_n_sgl[1]*Smain_asc;
    tau_np_p_tl = tau_np_p_sgl[0] + tau_np_p_sgl[1]*Stail_asc;
    tau_np_n_tl = tau_np_n_sgl[0] + tau_np_n_sgl[1]*Stail_asc;
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(tau_np_p_tl, __POSTIVEPARAM_MIN__);
    constraint(tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp,0.,1.);
  constraint(fd,0.,1.);
  return;
}

void RkRdetRnpPdf::Rnp_param(const int ntrk_asc,
                      const double& xi_asc, const double& st_asc,
                      const double& Smain_asc, const double& Stail_asc,
                      double& fd, double& fp,
                      double& tau_np_p, double& tau_np_n,
                      double& tau_np_p_tl, double& tau_np_n_tl){
  //if((param->fn_np_mlt)[0]>DTRES_EXTERAM_THRE){
  if(fn_np_mlt >= -0.00001) Rnp_param_10(ntrk_asc, xi_asc, st_asc,                        fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  else                         Rnp_param_03(ntrk_asc, xi_asc, st_asc, Smain_asc,  Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  return;
}

void RkRdetRnpPdf::Rnp_param_10(const int ntrk_asc,
                  const double& xi_asc, const double& st_asc,
                  double& fd, double& fp,
                  double& tau_np_p, double& tau_np_n,
                  double& tau_np_p_tl, double& tau_np_n_tl){
  if(ntrk_asc>1){ /*  multiple track vertex  */
    double st_asc_forf = (st_asc > rnp_kink_st ? rnp_kink_st : st_asc);
    double xi_asc_forf = (xi_asc > rnp_kink_xi ? rnp_kink_xi : xi_asc);

    double f_neg = fn_np_mlt;
    double f_delta = fd_np_mlt[keeptagl] + fd_np_st_mlt*st_asc_forf + fd_np_xi_mlt*xi_asc_forf + fd_np_stxi_mlt*st_asc_forf*xi_asc_forf;

    constraint(f_delta, 0.0, 1.0);
    constraint(f_neg, 0.0, 1.0);

    fd = (1. - f_neg) * f_delta;
    fp = (fd < 1. ? (1. - f_neg)*(1. - f_delta)/(1. - fd) : 0.5);
    tau_np_p = Snp_global*(tau_np_p_mlt[0] + tau_np_p_mlt[1]*st_asc + tau_np_p_xi_mlt*xi_asc + tau_np_p_stxi_mlt*st_asc*xi_asc);
    tau_np_n = Snp_global*(tau_np_n_mlt[0] + tau_np_n_mlt[1]*st_asc + tau_np_n_xi_mlt*xi_asc + tau_np_n_stxi_mlt*st_asc*xi_asc);
  }else{ /*  single track vertex */
    fd = fd_np_sgl[keeptagl];
    fp = fp_np_sgl;
    tau_np_p = Snp_global*(tau_np_p_sgl[0] + tau_np_p_sgl[1]*st_asc);
    tau_np_n = Snp_global*(tau_np_n_sgl[0] + tau_np_n_sgl[1]*st_asc);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp, 0.0, 1.0);
  tau_np_p_tl = tau_np_p;
  tau_np_n_tl = tau_np_n;
  return;
}

void RkRdetRnpPdf::Rnp_param_03(const int ntrk_asc,
                  const double& xi_asc, const double& st_asc,
                  const double& Smain_asc, const double& Stail_asc,
                  double& fd, double& fp,
                  double& tau_np_p, double& tau_np_n,
                  double& tau_np_p_tl, double& tau_np_n_tl){
  if(ntrk_asc>1){ /*  multiple track vertex  */
    const double Smain_asc_in = (Snp <= DTRES_EXTERAM_THRE) ? Smain_asc : (1.0+Snp*xi_asc)*st_asc;
    fd = fd_np_mlt[keeptagl];
    fp = fp_np_mlt;
    tau_np_p = Snp_global*(tau_np_p_mlt[0] + tau_np_p_mlt[1]*Smain_asc_in);
    tau_np_n = Snp_global*(tau_np_n_mlt[0] + tau_np_n_mlt[1]*Smain_asc_in);
    const double Stail_asc_in = (Snp <= DTRES_EXTERAM_THRE) ? Stail_asc : (1.0+Snp*xi_asc)*st_asc;
    tau_np_p_tl = Snp_global*(tau_np_p_mlt[0] + tau_np_p_mlt[1]*Stail_asc_in);
    tau_np_n_tl = Snp_global*(tau_np_n_mlt[0] + tau_np_n_mlt[1]*Stail_asc_in);
  }else{ /*  single track vertex */
    const double Smain_asc_in = (Snp <= DTRES_EXTERAM_THRE) ? Smain_asc : st_asc;
    const double Stail_asc_in = (Snp <= DTRES_EXTERAM_THRE) ? Stail_asc : st_asc;
    fd = fd_np_sgl[keeptagl];
    fp = fp_np_sgl;
    tau_np_p    = Snp_global*(tau_np_p_sgl[0] + tau_np_p_sgl[1]*Smain_asc_in);
    tau_np_n    = Snp_global*(tau_np_n_sgl[0] + tau_np_n_sgl[1]*Smain_asc_in);
    tau_np_p_tl = Snp_global*(tau_np_p_sgl[0] + tau_np_p_sgl[1]*Stail_asc_in);
    tau_np_n_tl = Snp_global*(tau_np_n_sgl[0] + tau_np_n_sgl[1]*Stail_asc_in);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_p_tl, __POSTIVEPARAM_MIN__);
  constraint(tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp,0.0,1.0);
  constraint(fd,0.0,1.0);
  return;
}

void RkRdetRnpPdf::swap_rnp_param(double& fp,
                    double& tau_np_p, double& tau_np_n,
                    double& tau_np_p_tl, double& tau_np_n_tl){
  fp = 1.0 - fp;
  const double taunp_n_tmp    = tau_np_n;
  tau_np_n = tau_np_p;
  tau_np_p = taunp_n_tmp;
  const double taunp_n_tl_tmp = tau_np_n_tl;
  tau_np_n_tl = tau_np_p_tl;
  tau_np_p_tl = taunp_n_tl_tmp;
  return;
}

void RkRdetRnpPdf::calc_vtxparam_asc(const int ntrk_asc, const double& sz_asc,
                              const double& chisq_z_asc, const int ndf_z_asc,
                              double& xi_asc, double& st_asc){
  xi_asc = (ndf_z_asc > 0 ? chisq_z_asc/ndf_z_asc : 1);
  st_asc = sz_asc*cm2ps;
  return;
}

void RkRdetRnpPdf::calc_vtxparam_rec(const int ntrk_rec, const double& sz_rec,
                              const double& chisq_z_rec, const int ndf_z_rec,
                              double& xi_rec, double& st_rec){
  xi_rec = (ndf_z_rec > 0 ? chisq_z_rec/ndf_z_rec : 1);
  st_rec = sz_rec*cm2ps;
  return;
}

double RkRdetRnpPdf::EfRkRdetRnp_fullrec(const double& x,
                           const int ntrk_rec, const double& sz_rec,
                           const double& chisq_z_rec, const int ndf_z_rec,
                           const int ntrk_asc, const double& sz_asc,
                           const double& chisq_z_asc, const int ndf_z_asc){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec,
             ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc,
             ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc,
            fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

//  cout << "tau = " << m_tau << " ";
  const double r_ckak = ck/ak;
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_full_sup(x,fact_n, ntau_n, fact_p, ntau_p,fd, fp, tau_np_p, tau_np_n,mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_full_sup(x,fact_n, ntau_n, fact_p, ntau_p,fd, fp, tau_np_p_tl, tau_np_n_tl,mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_full_sup(x,fact_n, ntau_n, fact_p, ntau_p,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_full_sup(x,fact_n, ntau_n, fact_p, ntau_p,fd,fp,tau_np_p_tl,tau_np_n_tl,mu_tt,sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                         +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    } else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  } else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdetRnp_full_sup(x,fact_n, ntau_n, fact_p, ntau_p,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::EfRkRdetRnp_full_sup(const double& x,
                                   const double& fact_n, const double& ntau_n,
                                   const double& fact_p, const double& ntau_p,
                                   const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                   const double& mu_det, const double& sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n*fd;
  double fEp = fact_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  add_EnEn_coef(fEn, fEn_np, fxEn,  ntau_n, tau_np_n, fact_n*nfn);
  add_EnEp_coef(fEn, fEp_np,        ntau_n, tau_np_p, fact_n*nfp);
  add_EpEn_coef(fEp, fEn_np,        ntau_p, tau_np_n, fact_p*nfn);
  add_EpEp_coef(fEp, fEp_np, fxEp,  ntau_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  En_conv_gauss(x, ntau_n,   mu_det, sigma_det);
  if(fEp!=0.0)    Li += fEp    *  Ep_conv_gauss(x, ntau_p,   mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(x, tau_np_n, mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(x, tau_np_p, mu_det, sigma_det);
  if(fxEn!=0.0)   Li += fxEn   * xEn_conv_gauss(x, ntau_n,   mu_det, sigma_det);
  if(fxEp!=0.0)   Li += fxEp   * xEp_conv_gauss(x, ntau_p,   mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::AfRkRdetRnp_fullrec(const double& x,
                           const int ntrk_rec, const double& sz_rec,
                           const double& chisq_z_rec, const int ndf_z_rec,
                           const int ntrk_asc, const double& sz_asc,
                           const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,fd, fp, tau_np_p, tau_np_n,mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,fd, fp, tau_np_p_tl, tau_np_n_tl,mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,fd, fp, tau_np_p_tl, tau_np_n_tl,mu_tt, sigma_tt);
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
    const double Li_tm    = AfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::AfRkRdetRnp_full_sup(const double& x,
                                   const double& fact_n, const double& ntau_n, const double& ndm_n,
                                   const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau,
                                   const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                   const double& mu_det, const double& sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_n*ndmtau;
  const double w_mn_p = -fact_p*ndmtau;
  double fAn = fact_n*fd;
  double fAp = fact_p*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
  add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np * En_conv_gauss(x, tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np * Ep_conv_gauss(x, tau_np_p,      mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::MfRkRdetRnp_fullrec(const double& x,
                           const int ntrk_rec, const double& sz_rec,
                           const double& chisq_z_rec, const int ndf_z_rec,
                           const int ntrk_asc, const double& sz_asc,
                           const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*m_tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec,
             ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc,
             ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc,
            fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,fd, fp, tau_np_p, tau_np_n,mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,fd, fp, tau_np_p_tl, tau_np_n_tl,mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,fd, fp, tau_np_p_tl, tau_np_n_tl,mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                           +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec,Smain_asc);
    const double Li_tm    = MfRkRdetRnp_full_sup(x,fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::MfRkRdetRnp_full_sup(const double& x,
                                   const double& fact_n, const double& ntau_n, const double& ndm_n,
                                   const double& fact_p, const double& ntau_p, const double& ndm_p,
                                   const double& ndmtau, const double& cktau,
                                   const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                   const double& mu_det, const double& sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_n*ndmtau;
  const double w_an_p = fact_p*ndmtau;
  double fMn = fact_n*fd;
  double fMp = fact_p*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);

  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  An_conv_gauss(x,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    *  Ap_conv_gauss(x,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    *  Mn_conv_gauss(x,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    *  Mp_conv_gauss(x,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_k!=0.0)  Li += fEn_k  *  En_conv_gauss(x, -cktau,         mu_det, sigma_det);
  if(fEp_k!=0.0)  Li += fEp_k  *  Ep_conv_gauss(x,  cktau,         mu_det, sigma_det);
  if(fxEn_k!=0.0) Li += fxEn_k * xEn_conv_gauss(x, -cktau,         mu_det, sigma_det);
  if(fxEp_k!=0.0) Li += fxEp_k * xEp_conv_gauss(x,  cktau,         mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(x,  tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(x,  tau_np_p,      mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                                const double& chisq_z_rec, const int ndf_z_rec,
                                const int ntrk_asc, const double& sz_asc,
                                const double& chisq_z_asc, const int ndf_z_asc){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);
  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double r_ckak = ck/ak;
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdetRnp_full_sup(fact_n, ntau_n, fact_p, ntau_p,
                                                    fd, fp, tau_np_p, tau_np_n, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdetRnp_full_sup(fact_n, ntau_n, fact_p, ntau_p,
                                                      fd, fp, tau_np_p_tl, tau_np_n_tl,
                                                      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
        const double mu_tm    = mu_tail_rec+mu_main_asc;
        const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
        const double Li_tm    = norm_EfRkRdetRnp_full_sup(fact_n, ntau_n, fact_p, ntau_p,
                                       fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
        const double mu_tt    = mu_tail_rec+mu_tail_asc;
        const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
        const double Li_tt    = norm_EfRkRdetRnp_full_sup(fact_n, ntau_n, fact_p, ntau_p,
                                       fd, fp, tau_np_p_tl, tau_np_n_tl,mu_tt, sigma_tt);
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
    const double Li_tm    = norm_EfRkRdetRnp_full_sup(fact_n, ntau_n, fact_p, ntau_p,
                                 fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_EfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n,
                                        const double& fact_p, const double& ntau_p,
                                        const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                        const double& mu_det, const double& sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n*fd;
  double fEp = fact_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  add_EnEn_coef(fEn,fEn_np,fxEn,ntau_n, tau_np_n, fact_n*nfn);
  add_EnEp_coef(fEn,fEp_np,     ntau_n, tau_np_p, fact_n*nfp);
  add_EpEn_coef(fEp,fEn_np,     ntau_p, tau_np_n, fact_p*nfn);
  add_EpEp_coef(fEp,fEp_np,fxEp,ntau_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(ll, ul, tau_np_n, mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_det, sigma_det);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                                const double& chisq_z_rec, const int ndf_z_rec,
                                const int ntrk_asc, const double& sz_asc,
                                const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*m_tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
                                                    fd, fp, tau_np_p, tau_np_n, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
                                                      fd, fp, tau_np_p_tl, tau_np_n_tl,mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
                                                          fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
                                                        fd, fp, tau_np_p_tl, tau_np_n_tl, mu_tt, sigma_tt);
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
    const double Li_tm    = norm_MfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
                                                      fd, fp, tau_np_p, tau_np_n,mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_MfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n,
                                        const double& fact_p, const double& ntau_p, const double& ndm_p,
                                        const double& ndmtau, const double& cktau,
                                        const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                        const double& mu_det, const double& sigma_det){

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_n*ndmtau;
  const double w_an_p = fact_p*ndmtau;
  double fMn = fact_n*fd;
  double fMp = fact_p*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  norm_An_conv_gauss(ll, ul,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    *  norm_Ap_conv_gauss(ll, ul,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    *  norm_Mn_conv_gauss(ll, ul,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    *  norm_Mp_conv_gauss(ll, ul,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_k!=0.0)  Li += fEn_k  *  norm_En_conv_gauss(ll, ul, -cktau,         mu_det, sigma_det);
  if(fEp_k!=0.0)  Li += fEp_k  *  norm_Ep_conv_gauss(ll, ul,  cktau,         mu_det, sigma_det);
  if(fxEn_k!=0.0) Li += fxEn_k * norm_xEn_conv_gauss(ll, ul, -cktau,         mu_det, sigma_det);
  if(fxEp_k!=0.0) Li += fxEp_k * norm_xEp_conv_gauss(ll, ul,  cktau,         mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(ll, ul,  tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(ll, ul,  tau_np_p,      mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_fullrec(const int ntrk_rec, const double& sz_rec,
                                const double& chisq_z_rec, const int ndf_z_rec,
                                const int ntrk_asc, const double& sz_asc,
                                const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc,
            fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  swap_rnp_param(fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
                                                    fd, fp, tau_np_p, tau_np_n, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
                                                      fd, fp, tau_np_p_tl, tau_np_n_tl, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
                                                        fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
                                                        fd, fp, tau_np_p_tl, tau_np_n_tl, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                         +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      return (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_AfRkRdetRnp_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
                                                      fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_AfRkRdetRnp_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n,
                                        const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau,
                                        const double& fd, const double& fp, const double& tau_np_p, const double& tau_np_n,
                                        const double& mu_det, const double& sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_n*ndmtau;
  const double w_mn_p = -fact_p*ndmtau;
  double fAn = fact_n*fd;
  double fAp = fact_p*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(fAn, fMn, fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_AnEp_coef(fAn, fMn, fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_MnEn_coef(fMn, fAn, fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  add_MnEp_coef(fMn, fAn, fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  add_ApEn_coef(fAp, fMp, fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_ApEp_coef(fAp, fMp, fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
  add_MpEn_coef(fMp, fAp, fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  add_MpEp_coef(fMp, fAp, fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np * norm_En_conv_gauss(ll, ul, tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np * norm_Ep_conv_gauss(ll, ul, tau_np_p,      mu_det, sigma_det);
  return Li;
}

double RkRdetRnpPdf::EfRkRdet_fullrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  const double r_ckak = ck/ak;
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdet_full_sup(x, fact_n, ntau_n, fact_p, ntau_p,mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdet_full_sup(x, fact_n, ntau_n, fact_p, ntau_p,mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdet_full_sup(x, fact_n, ntau_n, fact_p, ntau_p, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdet_full_sup(x, fact_n, ntau_n, fact_p, ntau_p, mu_tt, sigma_tt);
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
    const double Li_tm    = EfRkRdet_full_sup(x, fact_n, ntau_n, fact_p, ntau_p, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::EfRkRdet_full_sup(const double& x,const double& fact_n, const double& ntau_n,
                                const double& fact_p, const double& ntau_p,const double& mu_det, const double& sigma_det){
  return fact_n*En_conv_gauss(x, ntau_n, mu_det, sigma_det) + fact_p*Ep_conv_gauss(x, ntau_p, mu_det, sigma_det);
}

double RkRdetRnpPdf::AfRkRdet_fullrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tt, sigma_tt);
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
    const double Li_tm    = AfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::AfRkRdet_full_sup(const double& x, const double& fact_n, const double& ntau_n, const double& ndm_n,
                                const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& mu, const double& sigma){
  return fact_n*(An_conv_gauss(x, ntau_n, ndm_n, mu, sigma)-ndmtau*Mn_conv_gauss(x, ntau_n, ndm_n, mu, sigma))
       + fact_p*(Ap_conv_gauss(x, ntau_p, ndm_p, mu, sigma)-ndmtau*Mp_conv_gauss(x, ntau_p, ndm_p, mu, sigma));
}

double RkRdetRnpPdf::MfRkRdet_fullrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*m_tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tt, sigma_tt);
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
    const double Li_tm    = MfRkRdet_full_sup(x, fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::MfRkRdet_full_sup(const double& x,const double& fact_n, const double& ntau_n, const double& ndm_n,
                                const double& fact_p, const double& ntau_p, const double& ndm_p,const double& ndmtau, const double& cktau,const double& mu, const double& sigma){
  return fact_n*(Mn_conv_gauss(x, ntau_n, ndm_n, mu, sigma)+ndmtau*An_conv_gauss(x, ntau_n, ndm_n, mu, sigma))
       + fact_p*(Mp_conv_gauss(x, ntau_p, ndm_p, mu, sigma)+ndmtau*Ap_conv_gauss(x, ntau_p, ndm_p, mu, sigma));
}

double RkRdetRnpPdf::norm_EfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec,mu_main_rec,Smain_rec,mu_tail_rec,Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc,mu_main_asc,Smain_asc,mu_tail_asc,Stail_asc);

  const double r_ckak = ck/ak;
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdet_full_sup(fact_n, ntau_n, fact_p, ntau_p,mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdet_full_sup(fact_n, ntau_n, fact_p, ntau_p, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdet_full_sup(fact_n, ntau_n, fact_p, ntau_p, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdet_full_sup(fact_n, ntau_n, fact_p, ntau_p, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
                        +(1.0-ftail_rec)*ftail_asc*Li_mt +       ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  } else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdet_full_sup(fact_n, ntau_n, fact_p, ntau_p,mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_EfRkRdet_full_sup(const double& fact_n, const double& ntau_n,const double& fact_p, const double& ntau_p,const double& mu_det, const double& sigma_det){
  return fact_n*norm_En_conv_gauss(ll, ul, ntau_n, mu_det, sigma_det) + fact_p*norm_Ep_conv_gauss(ll, ul, ntau_p, mu_det, sigma_det);
}

double RkRdetRnpPdf::norm_AfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tt, sigma_tt);
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
    const double Li_tm    = norm_AfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_AfRkRdet_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& mu, const double& sigma){
  return fact_n*(norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma)-ndmtau*norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma))
       + fact_p*(norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma)-ndmtau*norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma));
}

double RkRdetRnpPdf::norm_MfRkRdet_fullrec(const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec,
                        const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*m_dm*m_tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = m_dm/(ak-ck);
  const double ndm_p  = m_dm/(ak+ck);
  const double ntau_n = m_tau*(ak-ck);
  const double ntau_p = m_tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*m_tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tt, sigma_tt);
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
    const double Li_tm    = norm_MfRkRdet_full_sup(fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double RkRdetRnpPdf::norm_MfRkRdet_full_sup(const double& fact_n, const double& ntau_n, const double& ndm_n, const double& fact_p, const double& ntau_p, const double& ndm_p, const double& ndmtau, const double& cktau, const double& mu, const double& sigma){
  return fact_n*(norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma)+ndmtau*norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma))
       + fact_p*(norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma)+ndmtau*norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma));
}

double RkRdetRnpPdf::RascRnp(const double& x, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc){
  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double Li_md = gaussian(x, mu_main_asc, Smain_asc);
  const double Li_mp = Ep_conv_gauss(x, tau_np_p, mu_main_asc, Smain_asc);
  const double Li_mn = En_conv_gauss(x, tau_np_n, mu_main_asc, Smain_asc);
  const double Li_me = fp*Li_mp + (1.0 - fp) * Li_mn;
  const double Li_mt = fd*Li_md + (1.0 - fd) * Li_me;
  if(ftail_asc==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = gaussian(x, mu_tail_asc, Stail_asc);
  const double Li_tp = Ep_conv_gauss(x, tau_np_p_tl, mu_tail_asc, Stail_asc);
  const double Li_tn = En_conv_gauss(x, tau_np_n_tl, mu_tail_asc, Stail_asc);
  const double Li_te = fp*Li_tp + (1.0 - fp) * Li_tn;
  const double Li_tt = fd*Li_td + (1.0 - fd) * Li_te;
  double Li  = (1.0 - ftail_asc) * Li_mt + ftail_asc * Li_tt;
  return Li;
}

double RkRdetRnpPdf::norm_RascRnp(const int ntrk_asc, const double& sz_asc,const double& chisq_z_asc, const int ndf_z_asc){
  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, xi_asc, st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(ntrk_asc, xi_asc, st_asc, Smain_asc, Stail_asc, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);

  const double Li_md = norm_gaussian(ll, ul, mu_main_asc, Smain_asc);
  const double Li_mp = norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_main_asc, Smain_asc);
  const double Li_mn = norm_En_conv_gauss(ll, ul, tau_np_n, mu_main_asc, Smain_asc);
  const double Li_me = fp*Li_mp + (1.0 - fp) * Li_mn;
  const double Li_mt = fd*Li_md + (1.0 - fd) * Li_me;
  if(ftail_asc==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = norm_gaussian(ll, ul, mu_tail_asc, Stail_asc);
  const double Li_tp = norm_Ep_conv_gauss(ll, ul, tau_np_p_tl, mu_tail_asc, Stail_asc);
  const double Li_tn = norm_En_conv_gauss(ll, ul, tau_np_n_tl, mu_tail_asc, Stail_asc);
  const double Li_te = fp*Li_tp + (1.0 - fp) * Li_tn;
  const double Li_tt = fd*Li_td + (1.0 - fd) * Li_te;
  double Li  = (1.0 - ftail_asc) * Li_mt + ftail_asc * Li_tt;
  return Li;
}

double RkRdetRnpPdf::Rrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);
  const double Li_mn = gaussian(x, mu_main_rec, Smain_rec);
  if(ftail_rec>0.0){ /* Mainly single track case */
    const double Li_tl = gaussian(x, mu_tail_rec, Stail_rec);
    const double Li    = (1.0-ftail_rec)*Li_mn + ftail_rec*Li_tl;
    return Li;
  }
  return Li_mn;
}

double RkRdetRnpPdf::norm_Rrec(const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, xi_rec, st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec);
  const double Li_mn = norm_gaussian(ll, ul, mu_main_rec, Smain_rec);
  if(ftail_rec>0.0){ /* Mainly multiple track case */
    const double Li_tl = norm_gaussian(ll, ul, mu_tail_rec, Stail_rec);
    const double Li    = (1.0-ftail_rec)*Li_mn + ftail_rec*Li_tl;
    return Li;
  }
  return Li_mn;
}

double RkRdetRnpPdf::PdfRrec(const double& x,const int ntrk_rec, const double& sz_rec,const double& chisq_z_rec, const int ndf_z_rec){
  const double pdf = Rrec(x,ntrk_rec,sz_rec,chisq_z_rec,ndf_z_rec);
  const double norm_pdf = norm_Rrec(ntrk_rec,sz_rec,chisq_z_rec,ndf_z_rec);
  return AddOutlier(x,pdf,ntrk_rec,2,norm_pdf);
}

double RkRdetRnpPdf::PdfRascRnp(const double& x, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc){
  const double pdf = RascRnp(x,ntrk_asc,sz_asc,chisq_z_asc,ndf_z_asc);
  const double norm_pdf = norm_RascRnp(ntrk_asc,sz_asc,chisq_z_asc,ndf_z_asc);
  return AddOutlier(x,pdf,2,ntrk_asc,norm_pdf);
}

double RkRdetRnpPdf::PdfAB(const double& dt,const int ntrk_rec,const double& sz_rec,const double& chisq_z_rec,const int ndf_z_rec,const int ntrk_asc,const double& sz_asc,const double& chisq_z_asc,const int ndf_z_asc, const bool otlr, const bool no_interf){
  const double Ef = EfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  const double norm_Ef = norm_EfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  if(!no_interf){
    const double Mf = MfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double Af = AfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double norm_Mf = norm_MfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  //  const double norm_Af = norm_AfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  //  cout << "ak = " << ak << ", ck = " << ck << endl;
    const double pdf = Ef*cexp - 0.5*amix*flv/m_tau*(A*Mf + B*Af);
    const double pdf_norm = norm_Ef*cexp - 0.5*amix*flv/m_tau*A*norm_Mf;// + 0.5/m_tau*B*norm_Af;
    if(false && (pdf<0 || pdf_norm<=0)){
      cout << "PdfAB. pdf: " << pdf << ", norm: " << pdf_norm;
      cout << ", cexp: " << cexp << ", amix: " << amix << ", flv: " << flv;
      cout << endl;
  //    cout << "       Ef  = " << Ef      << ", Mf  = " << Mf      << ", Af  = " << Af      << endl;
  //    cout << "       Efn = " << norm_Ef << ", Mfn = " << norm_Mf << ", Afn = " << norm_Af << endl;
      return -fabs(pdf/pdf_norm);
    }
    if(otlr) return AddOutlier(dt,pdf,ntrk_rec,ntrk_asc,pdf_norm);
    else     return pdf/pdf_norm;
  } else{
    if(otlr) return AddOutlier(dt,Ef,ntrk_rec,ntrk_asc,norm_Ef);
    else     return Ef/norm_Ef;
  }
  return -999;
}

double RkRdetRnpPdf::Pdf(const double& dt, const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc, const bool otlr, const bool no_interf){
//  cout << "Pdf. flv: " << flv << ", amix: " << amix << ", cexp: " << cexp << endl;
  double m_pdf = 0l;
  const double Ef = EfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  const double norm_Ef = norm_EfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  if(!no_interf){
    const double Mf = MfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double Af = AfRkRdetRnp_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double norm_Mf = norm_MfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
//    const double norm_Af = norm_AfRkRdetRnp_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double pdf = (K+Kb)*Ef*cexp - 0.5*amix/m_tau*((K-Kb)*Mf + 2.*xi*sqrt(K*Kb)*Af*(C*sin2beta+S*cos2beta));
    const double pdf_norm = (K+Kb)*norm_Ef*cexp - 0.5*amix/m_tau*(K-Kb)*norm_Mf;// + 2.*flv*xi*sqrt(K*Kb)*norm_Af*(C*sin2beta+S*cos2beta);
    if(pdf<=0.00000 || pdf_norm<=0.00000){
      cout << "Pdf1: pdf = " << pdf << ", norm = " << pdf_norm << ", cexp: " << cexp << ", amix: " << amix << endl;
      return pdf/pdf_norm;
    }
    if(otlr) m_pdf = AddOutlier(dt,pdf,ntrk_rec,ntrk_asc,pdf_norm);
    else     m_pdf = pdf/pdf_norm;
//    cout << "RkRdetRnpPdf::Pdf " << m_pdf << " " << pdf << " " << pdf_norm << endl;
    return m_pdf;
  } else{
    if(Ef<=0.0000 || norm_Ef<=0.00000){
      cout << "Pdf: pdf = " << Ef << ", norm = " << norm_Ef << endl;
      return Ef/norm_Ef;
    }
    if(otlr) m_pdf = AddOutlier(dt,Ef,ntrk_rec,ntrk_asc,norm_Ef);
    else     m_pdf = Ef/norm_Ef;
//    cout << "RkRdetRnpPdf::Pdf " << m_pdf << " " << Ef << " " << norm_Ef << endl;
    return m_pdf;
  }
  return -999;
}

double RkRdetRnpPdf::NoNPPdf(const double& dt, const int ntrk_rec, const double& sz_rec, const double& chisq_z_rec, const int ndf_z_rec, const int ntrk_asc, const double& sz_asc, const double& chisq_z_asc, const int ndf_z_asc, const bool otlr, const bool no_interf){
  const double Ef = EfRkRdet_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  const double norm_Ef = norm_EfRkRdet_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
  if(!no_interf){
    const double Mf = MfRkRdet_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double Af = AfRkRdet_fullrec(dt, ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double norm_Mf = norm_MfRkRdet_fullrec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc);
    const double pdf = (K+Kb)*Ef - 0.5/m_tau*flv*(K-Kb)*Mf + 0.5/m_tau*2.*flv*xi*sqrt(K*Kb)*Af*(C*sin2beta+S*cos2beta);
    const double pdf_norm = (K+Kb)*norm_Ef - 0.5/m_tau*flv*(K-Kb)*norm_Mf;
    if(pdf<0 || pdf_norm<=0){
      cout << "Pdf: pdf = " << pdf << ", norm = " << pdf_norm << endl;
      return 0;
    }
    if(otlr) return AddOutlier(dt,pdf,ntrk_rec,ntrk_asc,pdf_norm);
    else     return pdf/pdf_norm;
  } else{
    if(otlr) return AddOutlier(dt,Ef,ntrk_rec,ntrk_asc,norm_Ef);
    else     return Ef/norm_Ef;
  }
  return -999;
}

double RkRdetRnpPdf::AddOutlier(const double& x,const double& Lin,
                         const int ntrk_rec, const int ntrk_asc,
                         const double& nLi,const double& alpha){
  return Add_Outlier(x,Lin,ntrk_rec,ntrk_asc,nLi,alpha);
}

double RkRdetRnpPdf::Add_Outlier(const double& x, const double Lin,
                   const int ntrk_rec, const int ntrk_asc,
                   const double& nLi, const double& alpha){
  const double fol = ((m_svd == 2 || ntrk_rec>1) && ntrk_asc>1) ? f_ol_mul : f_ol_sgl;
//  const double fol = (ntrk_rec>1 && ntrk_asc>1) ? f_ol_mul : f_ol_sgl;
  const double m = 0.0;
  const double Lol = gaussian(x,m,sigma_ol);
  const double nLol = norm_gaussian(ll,ul,m,sigma_ol);
  const double Li = (1.0-fol)*Lin/nLi + fol*alpha*Lol/nLol;
//  cout << "Otlr: " << Li << " " << fol << " " << Lin << " " << nLi << " " << alpha << " " << Lol << " " << nLol << endl;
  return Li;
}
