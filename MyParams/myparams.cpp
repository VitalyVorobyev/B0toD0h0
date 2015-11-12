#include "myparams.h"
#include <iomanip>
#include "TMath.h"
#include "math.h"

//const double Karr[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
//const double Kbarr[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};
const double Karr[8] = { 0.168795,0.118502, 0.096164, 0.074168, 0.090667, 0.030814, 0.105442, 0.077863};
const double Kbarr[8]= { 0.088009,0.028629, 0.012360, 0.014960, 0.042507, 0.010055, 0.023245, 0.017819};

const double m_btau = 1.534;//1.520;//0.460/sol*cm2ps;//0.507;// ps
const double m_dm   = 0.510;//0.452;//0.510/sol;//0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
const double m_xi   = 1+m_btau*m_dm*m_btau*m_dm;

const double wbins[8]        = {0.0,0.1,0.25,0.5,0.625,0.75,0.875,1.01};
const double w_mc_svd1[7]    = {0.5,0.420827,0.300296,0.219317,0.154636,0.0916131,0.0228891};
const double w_data_svd1[7]  = {0.5,0.418852,0.329879,0.233898,0.170608,0.099791, 0.0228501};
const double w_mc_svd2[7]    = {0.5,0.412222,0.307838,0.212765,0.149933,0.0913264,0.0218754};
const double w_data_svd2[7]  = {0.5,0.418826,0.319303,0.222948,0.163191,0.104085, 0.0251454};

MyParams::MyParams(void){
  de_fit_min = -0.15;
  de_fit_max = 0.3;

  mbc_fit_min = 5.20;
  mbc_fit_max = 5.289;

//  de_min_pi0 = -0.128342;
  de_min_pi0 =       -0.1;
  de_max_pi0 =        0.0921343;

  de_min_dst0pi0 =   -0.1;
  de_max_dst0pi0 =    0.0647985;

  de_min_etagg =     -0.0994528;
  de_max_etagg =      0.0758825;

  de_min_etapgg =    -0.0873898;
  de_max_etapgg =     0.0642064;

  de_min_dst0etagg = -0.0853615;
  de_max_dst0etagg =  0.0630641;

  de_min_etappp =    -0.0534644;
  de_max_etappp =     0.0453713;

  de_min_etapppp =   -0.0394845;
  de_max_etapppp =    0.0357815;

  de_min_dst0etappp =-0.0567598;
  de_max_dst0etappp = 0.0459849;

  de_min_omega = -0.0565915;
  de_max_omega = 0.0467748;

  mbc_min_pi0 = 5.27093;
  mbc_max_pi0 = 5.28789;

  mbc_min_dst0pi0 = 5.27146;
  mbc_max_dst0pi0 = 5.28769;

  mbc_min_etagg = 5.27147;
  mbc_max_etagg = 5.28763;

  mbc_min_etapgg = 5.27255;
  mbc_max_etapgg = 5.28718;

  mbc_min_dst0etagg = 5.27182;
  mbc_max_dst0etagg = 5.28746;

  mbc_min_etappp = 5.27205;
  mbc_max_etappp = 5.28735;

  mbc_min_etapppp = 5.27288;
  mbc_max_etapppp = 5.28702;

  mbc_min_dst0etappp = 5.27182;
  mbc_max_dst0etappp = 5.28746;

  mbc_min_omega = 5.27198;
  mbc_max_omega = 5.28743;

  mbc_side = 5.25;

  KMass = 0.4975;
  mk_min = KMass - 0.009;
  mk_max = KMass + 0.009;

  Pi0Mass = 0.135;
  mpi0_min = 0.115738;
  mpi0_max = 0.153714;

  EtaMass = 0.5478;
  metagg_min = 0.517603;
  metagg_max = 0.573654;
  metappp_min = 0.537626;
  metappp_max = 0.557443;

  DMass = 1.865;
  md_min = 1.85155;
  md_max = 1.87833;

  dmetapgg_min = 0.401668;
  dmetapgg_max = 0.417717;

  dmetapppp_min = 0.402493;
  dmetapppp_max = 0.417109;

  OmegaMass = 0.781;
  momega_min = 0.760430;
  momega_max = 0.803869;

  dmdst0_min = 0.1386;
  dmdst0_max = 0.1457;

  atckpi_cut = 1.;

  bdt_cut_pi0    = 0.27;
  bdt_cut_etagg  = 0.15;
  bdt_cut_etappp = 0.18;
  bdt_cut_omega  = 0.18;

  lh0_cut_etap    = 0.56;
  lh0_cut_dst0pi0 = 0.60;
  lh0_cut_dst0eta = 0.68;

  // pi0 dE-Mbc shape //
  m_a_mbc0_pi1m =-0.07404;
  m_b_mbc0_pi1m =-0.005278;
  m_c_mbc0_pi1m = 5.27994e+00;
  m_a_s_pi1m    = 0.1294;
  m_b_s_pi1m    =-0.0009067;
  m_c_s_pi1m    = 0.002836;
  m_alpha_nks_pi1m = 0.139;
  //de
  m_nl_pi1m     = 2;
  m_alphal_pi1m = 7.21400e-01;
  m_nr_pi1m     = 2;
  m_alphar_pi1m =-1.93372e+00;
  m_de0_pi1m    = 2.06276e-02;
  m_deCBl_pi1m  =-1.12822e-02;
  m_deCBr_pi1m  =-1.61835e-02;
  m_fCBl_pi1m   = 6.17109e-01;
  m_fCBr_pi1m   = 2.37276e-01;
  m_s1_pi1m     = 1.95009e-02;
  m_sCBl_pi1m   = 2.82849e-02;
  m_sCBr_pi1m   = 4.27106e-02;

  // eta -> gg dE-Mbc shape //
  // mbc
  m_a_mbc0_eta10m1 =-0.06697;
  m_b_mbc0_eta10m1 =-0.007232;
  m_c_mbc0_eta10m1 = 5.28;
  m_a_s_eta10m1 = 0.1437;
  m_b_s_eta10m1 =-0.000311;
  m_c_s_eta10m1 = 0.002828;
  m_alpha_nks_eta10m1 = 0.139;

  // de
  m_alphal_eta10m1 = 1.76432e+00;
  m_alphar_eta10m1 =-2.28305e+00;
  m_de0_eta10m1    = 3.93738e-03;
  m_deCBl_eta10m1  =-2.32713e-02;
  m_deCBr_eta10m1  =-5.13742e-02;
  m_fCBl_eta10m1   = 4.33081e-01;
  m_fCBr_eta10m1   = 1.61207e-01;
  m_nl_eta10m1     = 2;
  m_nr_eta10m1     = 2;
  m_s1_eta10m1     = 2.10044e-02;
  m_sCBl_eta10m1   = 3.60619e-02;
  m_sCBr_eta10m1   = 7.09575e-02;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // mbc
  m_c0_eta201        = 0.0006327;
  m_c1_eta201        = 0.006709;
  m_c2_eta201        = 0.01999;
  m_a_s_eta201       = 0.1868;
  m_b_s_eta201       = 0.001764;
  m_c_s_eta201       = 0.002715;
  m_alpha_nks_eta201 = 0.139;

  m_a_mbc0_eta205    =-0.06013;
  m_b_mbc0_eta205    =-0.0004072;
  m_c_mbc0_eta205    = 5.28;
  m_a_s_eta205       = 0.1223;
  m_b_s_eta205       =-0.003813;
  m_c_s_eta205       = 0.003251;
  m_alpha_nks_eta205 = 0.139;

  // de
  m_alphal_eta201 = 1.34806e+00;
  m_alphar_eta201 =-2.00612e+00;
  m_de0_eta201    = 8.71661e-04;
  m_deCBl_eta201  =-1.19361e-02;
  m_deCBr_eta201  =-4.14426e-03;
  m_fCBl_eta201   = 6.72362e-02;
  m_fCBr_eta201   = 4.97299e-01;
  m_nl_eta201     = 2;
  m_nr_eta201     = 2;
  m_s1_eta201     = 1.16465e-02;
  m_sCBl_eta201   = 4.93207e-02;
  m_sCBr_eta201   = 2.15750e-02;

  m_alphal_eta205 = 7.11846e-01;
  m_alphar_eta205 = 7.11846e-01;// fixed
  m_de0_eta205    =-1.71891e-01;
  m_deCBl_eta205  =-2.51194e-02;
  m_deCBr_eta205  = 0;// fixed
  m_fCBl_eta205   = 6.22736e-01;
  m_fCBr_eta205   = 0;// fixed
  m_nl_eta205     = 2;
  m_nr_eta205     = 2;
  m_s1_eta205     = 1.45650e-01;
  m_sCBl_eta205   = 3.46718e-02;
  m_sCBr_eta205   = 0.01;// fixed

  // omega dE-Mbc shape //
  // mbc
  m_mbc0_omega1     = 5.28;
  m_c0_omega1       = 0.0006307;
  m_c1_omega1       = 0.007338;
  m_c2_omega1       = 0.02292;
  m_a_s_omega1      = 0.1595;
  m_b_s_omega1      = 0.000648;
  m_c_s_omega1      = 0.00216002;
  m_alpha_nks_omega1= 0.139;

  m_a_mbc0_omega5   =-0.06405;
  m_b_mbc0_omega5   =-0.0005805;
  m_c_mbc0_omega5   = 5.28;
  m_a_s_omega5      = 0.1315;
  m_b_s_omega5      =-0.00421;
  m_c_s_omega5      = 0.003216;
  m_alpha_nks_omega5= 0.139;

  // de
  m_alphal_omega201 = 3.00093e+00;
  m_alphar_omega201 =-1.80433e+00;
  m_de0_omega201    = 2.95816e-04;
  m_deCBl_omega201  =-3.78618e-04;
  m_deCBr_omega201  = 1.49255e-03;
  m_fCBl_omega201   = 5.33498e-01;
  m_fCBr_omega201   = 2.21978e-02;
  m_nl_omega201     = 2;// fixed
  m_nr_omega201     = 2;// fixed
  m_s1_omega201     = 2.19207e-02;
  m_sCBl_omega201   = 1.18658e-02;
  m_sCBr_omega201   = 4.95480e-02;

  m_alphal_omega205 = 1.09220e+00;
  m_alphar_omega205 =-1.09220e+00;// fixed
  m_de0_omega205    =-1.39502e-01;
  m_deCBl_omega205  =-2.66137e-02;
  m_deCBr_omega205  = 0;// fixed
  m_fCBl_omega205   = 8.26807e-01;
  m_fCBr_omega205   = 0;// fixed
  m_nl_omega205     = 2;// fixed
  m_nr_omega205     = 2;// fixed
  m_s1_omega205     = 1.31967e-01;
  m_sCBl_omega205   = 3.96605e-02;
  m_sCBr_omega205   = 0.01;// fixed

  // * dE-Mbc Combinatorial * //
  // pi0 dE-Mbc shape //
//  m_cmb_argedge_qq_pi0 =5.28894e+00;
  m_cmb_argpar_qq_pi0 =-2.51995e+01;
  m_cmb_c1_pi0        =-2.99258e-01;
  m_cmb_c2_pi0        =-2.18132e-02;
  m_cmb_argpar_bb_pi0 =-2.82125e+01;
  m_cmb_mbc0_bb_pi0   = 5.27771e+00;
  m_cmb_mbcw_bb_pi0   = 1.93579e-02;
  m_cmb_f_g_bb_pi0    = 1.63495e-01;
  m_cmb_c10_pi0       = 1.85215e+01;
  m_cmb_c11_pi0       =-3.74672e+00;
  m_cmb_c20_pi0       = 2.89562e-01;
  m_cmb_fbb_pi0       = 3.14076e-01;

  // eta -> gg dE-Mbc shape //
  // gg QQ combinatorial
  m_cmb_argpar_qq_eta10 =-2.13770e+01;
  m_cmb_c1_eta10        =-3.13403e-01;
  m_cmb_c2_eta10        =-2.04026e-02;
  // gg bb  combinatorial
  m_cmb_argpar_bb_eta10 =-5.41534e+01;
  m_cmb_mbc0_bb_eta10   = 5.27912e+00;
  m_cmb_mbcw_bb_eta10   = 1.19674e-02;
  m_cmb_f_g_bb_eta10    = 5.17990e-01;
  m_cmb_c10_eta10       = 1.84714e+01;
  m_cmb_c11_eta10       =-3.74995e+00;
  m_cmb_c20_eta10       = 3.90586e-01;
  m_cmb_fbb_eta10       = 3.14566e-01;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // ppp QQ combinatorial
  m_cmb_argpar_qq_eta20 =-1.18837e+01;
  m_cmb_c1_eta20        =-4.47754e-01;
  m_cmb_c2_eta20        = 3.01497e-02;
  // ppp bb  combinatorial
  m_cmb_argpar_bb_eta20 =-4.57659e+01;
  m_cmb_mbc0_bb_eta20   = 5.27928e+00;
  m_cmb_mbcw_bb_eta20   = 5.52250e-03;
  m_cmb_f_g_bb_eta20    = 1.08317e-01;
  m_cmb_c10_eta20       = 1.43302e+01;
  m_cmb_c11_eta20       =-2.89865e+00;
  m_cmb_c20_eta20       = 1.76120e-01;
  m_cmb_fbb_eta20       = 2.09175e-01;

  // omega dE-Mbc shape //
  // QQ combinatorial
  m_cmb_argpar_qq_omega =-1.35393e+01;
  m_cmb_c1_omega        =-2.42659e-01;
  m_cmb_c2_omega        =-6.59673e-02;
  // bb combinatorial
  m_cmb_argpar_bb_omega =-5.53645e+01;
  m_cmb_mbc0_bb_omega   = 5.27700e+00;
  m_cmb_mbcw_bb_omega   = 1.57343e-02;
  m_cmb_f_g_bb_omega    = 0.69;
  m_cmb_c10_omega       = 1.08011e+01;
  m_cmb_c11_omega       =-2.18718e+00;
  m_cmb_c20_omega       = 8.78132e-02;

  // * dE-Mbc Partial BB * //
  // pi0 dE-Mbc shape //
  // Mbc
  m_peak_b_mbc0_pi0 = 5.27601e+00;
  m_peak_k_mbc0_pi0 =-2.12722e-02;
  m_peak_b_s_pi0    = 1.23366e-02;
  m_peak_k_s_pi0    = 5.84569e-02;
  // dE
  m_peak_de0r_pi0   =-1.11131e-01;
  m_peak_slopel_pi0 =-7.03669e+02;
  m_peak_sloper_pi0 =-1.52863e+00;
  m_peak_steep_pi0  = 9.52444e+00;
  m_peak_p5_pi0     = 9.39841e-01;

  // eta -> gg dE-Mbc shape //
  // Mbc
  m_peak_b_mbc0_eta10 = 5.28268e+00;
  m_peak_k_mbc0_eta10 = 4.61288e-02;
  m_peak_b_s_eta10    = 1.53838e-02;
  m_peak_k_s_eta10    = 6.67815e-02;
  // dE
  m_peak_de0r_eta10   =-1.11131e-01;
  m_peak_slopel_eta10 =-7.03669e+02;
  m_peak_sloper_eta10 =-1.52863e+00;
  m_peak_steep_eta10  = 9.52444e+00;
  m_peak_p5_eta10     = 9.39841e-01;

  // eta -> pi+pi-pi0 dE-Mbc shape //
  // Mbc
  m_peak_mbc0_eta20   = 5.27832e+00;
  m_peak_s_eta20      = 3.92662e-03;
  m_peak_fg_eta20     = 3.19374e-01;
  // dE
  m_peak_de0r_eta20   =-1.10127e-01;
  m_peak_slopel_eta20 =-1.83356e+02;
  m_peak_sloper_eta20 =-1.23245e+00;
  m_peak_steep_eta20  = 8.13795e+00;
  m_peak_p5_eta20     = 9.47792e-01;

  // omega dE-Mbc shape //
  // Mbc
  m_peak_mbc0_omega   = 5.28531e+00;
  m_peak_s_omega      = 6.39364e-03;
  m_peak_fg_omega     = 6.28950e-01;
  // dE
  m_peak_de0r_omega   =-1.10127e-01;
  m_peak_slopel_omega =-1.83356e+02;
  m_peak_sloper_omega =-1.23245e+00;
  m_peak_steep_omega  = 8.13795e+00;
  m_peak_p5_omega     = 9.47792e-01;
}

double MyParams::bdt_cut(const int mode,const int h0mode) const{
  switch(mode){
  case 1:  return bdt_cut_pi0;
  case 2:
    if(h0mode == 10) return bdt_cut_etagg;
    else             return bdt_cut_etappp;
  case 3:  return bdt_cut_omega;
  case 4:  return bdt_cut_omega;// rho
  default:
    return -1;
  }
  return -1;
}

double MyParams::lh0_cut(const int mode,const int h0mode) const{
  switch(mode){
  case 10: return lh0_cut_dst0pi0;
  case 20: return lh0_cut_dst0eta;
  case 5:  return lh0_cut_etap;
  default:
    return 0;
  }
  return 0;
}

double MyParams::get_de0_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_de0_pi1m;
  case 10: return 1.01146e-02;
  case 2:
    if(h0mode == 10)  return m_de0_eta10m1;
    else if(b0f == 1) return m_de0_eta201;
    else              return m_de0_eta205;
  case 20:
    if(h0mode == 10)  return 6.62297e-04;
//    else if(b0f == 1) return m_de0_eta201;
    else              return -7.59639e-04;
  case 3:
    if(b0f == 1) return m_de0_omega201;
    else         return m_de0_omega205;
  case 5:
    if(h0mode == 10)  return m_de0_eta10m1;
    else if(b0f == 1) return m_de0_eta201;
    else              return m_de0_eta205;
  }
  return -1;
}

double MyParams::get_s_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_s1_pi1m;
  case 10: return 1.40192e-02;
  case 2:
    if(h0mode == 10)  return m_s1_eta10m1;
    else if(b0f == 1) return m_s1_eta201;
    else              return m_s1_eta205;
  case 20:
    if(h0mode == 10)  return 1.88802e-02;
//    else if(b0f == 1) return m_s1_eta201;
    else              return 1.42191e-02;
  case 3:
    if(b0f == 1) return m_s1_omega201;
    else         return m_s1_omega205;
  case 5:
    if(h0mode == 10)  return m_s1_eta10m1;
    else if(b0f == 1) return m_s1_eta201;
    else              return m_s1_eta205;
  }
  return -1;
}

double MyParams::get_de0_CBl_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_deCBl_pi1m;
  case 10: return -1.89412e-02;
  case 2:
    if(h0mode == 10)  return m_deCBl_eta10m1;
    else if(b0f == 1) return m_deCBl_eta201;
    else              return m_deCBl_eta205;
  case 20:
    if(h0mode == 10)  return -2.95909e-02;
//    else if(b0f == 1) return m_deCBl_eta201;
    else              return -2.25347e-02;
  case 3:
    if(b0f == 1) return m_deCBl_omega201;
    else         return m_deCBl_omega205;
  case 5:
    if(h0mode == 10)  return m_deCBl_eta10m1;
    else if(b0f == 1) return m_deCBl_eta201;
    else              return m_deCBl_eta205;
  }
  return -1;
}

double MyParams::get_s_CBl_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_sCBl_pi1m;
  case 10: return 3.59628e-02;
  case 2:
    if(h0mode == 10)  return m_sCBl_eta10m1;
    else if(b0f == 1) return m_sCBl_eta201;
    else              return m_sCBl_eta205;
  case 20:
    if(h0mode == 10)  return 3.93993e-02;
//    else if(b0f == 1) return m_sCBl_eta201;
    else              return 3.69556e-02;
  case 3:
    if(b0f == 1) return m_sCBl_omega201;
    else         return m_sCBl_omega205;
  case 5:
    if(h0mode == 10)  return m_sCBl_eta10m1;
    else if(b0f == 1) return m_sCBl_eta201;
    else              return m_sCBl_eta205;
  }
  return -1;
}

double MyParams::get_alphal_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_alphal_pi1m;
  case 10: return 7.63149e-01;
  case 2:
    if(h0mode == 10)  return m_alphal_eta10m1;
    else if(b0f == 1) return m_alphal_eta201;
    else              return m_alphal_eta205;
  case 20:
    if(h0mode == 10)  return 1.37630e+00;
//    else if(b0f == 1) return m_alphal_eta201;
    else              return 1.19340e+00;
  case 3:
    if(b0f == 1) return m_alphal_omega201;
    else         return m_alphal_omega205;
  case 5:
    if(h0mode == 10)  return m_alphal_eta10m1;
    else if(b0f == 1) return m_alphal_eta201;
    else              return m_alphal_eta205;
  }
  return -1;
}

double MyParams::get_de0_CBr_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1: return m_deCBr_pi1m;
  case 10: return m_deCBr_pi1m;
  case 2:
    if(h0mode == 10)  return m_deCBr_eta10m1;
    else if(b0f == 1) return m_deCBr_eta201;
    else              return m_deCBr_eta205;
  case 20:
    if(h0mode == 10)  return m_deCBr_eta10m1;
    else if(b0f == 1) return m_deCBr_eta201;
    else              return m_deCBr_eta205;
  case 3:
    if(b0f == 1) return m_deCBr_omega201;
    else         return m_deCBr_omega205;
  case 5:
    if(h0mode == 10)  return m_deCBr_eta10m1;
    else if(b0f == 1) return m_deCBr_eta201;
    else              return m_deCBr_eta205;
  }
  return -1;
}

double MyParams::get_s_CBr_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_sCBr_pi1m;
  case 10: return m_sCBr_pi1m;
  case 2:
    if(h0mode == 10)  return m_sCBr_eta10m1;
    else if(b0f == 1) return m_sCBr_eta201;
    else              return m_sCBr_eta205;
  case 20:
    if(h0mode == 10)  return m_sCBr_eta10m1;
    else if(b0f == 1) return m_sCBr_eta201;
    else              return m_sCBr_eta205;
  case 3:
    if(b0f == 1) return m_sCBr_omega201;
    else         return m_sCBr_omega205;
  case 5:
    if(h0mode == 10)  return m_sCBr_eta10m1;
    else if(b0f == 1) return m_sCBr_eta201;
    else              return m_sCBr_eta205;
  }
  return -1;
}

double MyParams::get_alphar_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_alphar_pi1m;
  case 10: return m_alphar_pi1m;
  case 2:
    if(h0mode == 10)  return m_alphar_eta10m1;
    else if(b0f == 1) return m_alphar_eta201;
    else              return m_alphar_eta205;
  case 20:
    if(h0mode == 10)  return m_alphar_eta10m1;
    else if(b0f == 1) return m_alphar_eta201;
    else              return m_alphar_eta205;
  case 3:
    if(b0f == 1) return m_alphar_omega201;
    else         return m_alphar_omega205;
  case 5:
    if(h0mode == 10)  return m_alphar_eta10m1;
    else if(b0f == 1) return m_alphar_eta201;
    else              return m_alphar_eta205;
  }
  return -1;
}

double MyParams::get_f_CBl_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_fCBl_pi1m;
  case 10: return 8.80589e-01;
  case 2:
    if(h0mode == 10)  return m_fCBl_eta10m1;
    else if(b0f == 1) return m_fCBl_eta201;
    else              return m_fCBl_eta205;
  case 20:
    if(h0mode == 10)  return 5.97526e-01;
//    else if(b0f == 1) return m_fCBl_eta201;
    else              return 4.12468e-01;
  case 3:
    if(b0f == 1) return m_fCBl_omega201;
    else         return m_fCBl_omega205;
  case 5:
    if(h0mode == 10)  return m_fCBl_eta10m1;
    else if(b0f == 1) return m_fCBl_eta201;
    else              return m_fCBl_eta205;
  }
  return -1;
}

double MyParams::get_f_CBr_de_sig(const int mode,const int h0mode, const int b0f) const{
  switch(mode){
  case 1:  return m_fCBr_pi1m;
  case 10: return 0;
  case 2:
    if(h0mode == 10)  return m_fCBr_eta10m1;
    else if(b0f == 1) return m_fCBr_eta201;
    else              return m_fCBr_eta205;
  case 20:
    if(h0mode == 10)  return 0;
//    else if(b0f == 1) return m_fCBr_eta201;
    else              return 0;
  case 3:
    if(b0f == 1) return m_fCBr_omega201;
    else         return m_fCBr_omega205;
  case 5:
    if(h0mode == 10)  return m_fCBr_eta10m1;
    else if(b0f == 1) return m_fCBr_eta201;
    else              return m_fCBr_eta205;
  }
  return -1;
}

double MyParams::get_a_s_mbc_sig(const int mode,const int h0mode) const{
  switch(mode){
  case 1:  return m_a_s_pi1m;
  case 10: return m_a_s_pi1m;
  case 2:
    if(h0mode == 10) return m_a_s_eta10m1;
    else             return m_a_s_eta201;
  case 20:
    if(h0mode == 10) return m_a_s_eta10m1;
    else             return m_a_s_eta201;
  case 3: return m_a_s_omega1;
  case 5:
    if(h0mode == 10) return m_a_s_eta10m1;
    else             return m_a_s_eta201;
  }
  return 0;
}

double MyParams::get_b_s_mbc_sig(const int mode,const int h0mode) const{
  switch(mode){
  case 1: return m_b_s_pi1m;
  case 10: return m_b_s_pi1m;
  case 2:
    if(h0mode == 10) return m_b_s_eta10m1;
    else             return m_b_s_eta201;
  case 20:
    if(h0mode == 10) return m_b_s_eta10m1;
    else             return m_b_s_eta201;
  case 3: return m_b_s_omega1;
  }
  return 0;
}

double MyParams::get_c_s_mbc_sig(const int mode,const int h0mode) const{
  switch(mode){
  case 1: return m_c_s_pi1m;
  case 10: return m_c_s_pi1m;
  case 2:
    if(h0mode == 10) return m_c_s_eta10m1;
    else             return m_c_s_eta201;
  case 20:
    if(h0mode == 10) return m_c_s_eta10m1;
    else             return m_c_s_eta201;
  case 3: return m_c_s_omega1;
  }
  return 0;
}

double MyParams::get_a_mbc0_sig(const int mode,const int h0mode) const{
  switch(mode){
  case 1:  return m_a_mbc0_pi1m;
  case 10: return m_a_mbc0_pi1m;
  case 2:
    if(h0mode == 10) return m_a_mbc0_eta10m1;
    else             return 0.0006327;
  case 20: return m_a_mbc0_eta10m1;
  case 3:  return 0.0006307;
  }
  return 0;
}

double MyParams::get_b_mbc0_sig(const int mode, const int h0mode) const{
  switch(mode){
  case 1:  return m_b_mbc0_pi1m;
  case 10: return m_b_mbc0_pi1m;
  case 2:
    if(h0mode == 10) return m_b_mbc0_eta10m1;
    else             return 0.006709;
  case 20: return m_b_mbc0_eta10m1;
  case 3:  return 0.007338;
  }
  return 0;
}

double MyParams::get_c_mbc0_sig(const int mode, const int h0mode) const{
  switch(mode){
  case 1: return m_c_mbc0_pi1m;
  case 10: return m_c_mbc0_pi1m;
  case 2:
    if(h0mode == 10) return m_c_mbc0_eta10m1;
    else             return 0.01999;
  case 20: return m_c_mbc0_eta10m1;
  case 3:  return 0.02292;
  }
  return 0;
}

double MyParams::get_a_mbc0_sig_tail(const int mode) const{
  switch(mode){
  case 2: return m_a_mbc0_eta205;
  case 20: return m_a_mbc0_eta205;
  case 3: return m_a_mbc0_omega5;
  }
  return 0;
}

double MyParams::get_b_mbc0_sig_tail(const int mode) const{
  switch(mode){
  case 2: return m_b_mbc0_eta205;
  case 20: return m_b_mbc0_eta205;
  case 3: return m_b_mbc0_omega5;
  }
  return 0;
}

double MyParams::get_c_mbc0_sig_tail(const int mode) const{
  switch(mode){
  case 2: return m_c_mbc0_eta205;
  case 20: return m_c_mbc0_eta205;
  case 3: return m_c_mbc0_omega5;
  }
  return 0;
}

double MyParams::get_a_s_mbc_sig_tail(const int mode) const{
  switch(mode){
  case 2: return m_a_s_eta205;
  case 20: return m_a_s_eta205;
  case 3: return m_a_s_omega5;
  }
  return 0;
}

double MyParams::get_b_s_mbc_sig_tail(const int mode) const{
  switch(mode){
  case 2:  return m_b_s_eta205;
  case 20: return m_b_s_eta205;
  case 3:  return m_b_s_omega5;
  }
  return 0;
}

double MyParams::get_c_s_mbc_sig_tail(const int mode) const{
  switch(mode){
  case 2:  return m_c_s_eta205;
  case 20: return m_c_s_eta205;
  case 3:  return m_c_s_omega5;
  }
  return 0;
}

double MyParams::get_de_min_h0(const int mode,const int h0mode) const{
  switch(mode){
  case 1:            return de_min_pi0;
  case 10:           return de_min_dst0pi0;
  case 2:
    if(h0mode == 10) return de_min_etagg;
    else             return de_min_etappp;
  case 20:
    if(h0mode == 10) return de_min_dst0etagg;
    else             return de_min_dst0etappp;
  case 3:            return de_min_omega;
  case 5:
    if(h0mode == 10) return de_min_etapgg;
    else             return de_min_etapppp;
  }
  return -1;
}

double MyParams::get_de_max_h0(const int mode,const int h0mode) const{
  switch(mode){
  case 1:            return de_max_pi0;
  case 10:           return de_max_dst0pi0;
  case 2:
    if(h0mode == 10) return de_max_etagg;
    else             return de_max_etappp;
  case 20:
    if(h0mode == 10) return de_max_dst0etagg;
    else             return de_max_dst0etappp;
  case 3:            return de_max_omega;
  case 5:
    if(h0mode == 10) return de_max_etapgg;
    else             return de_max_etapppp;
  }
  return -1;
}

double MyParams::get_mbc_min_h0(const int mode,const int h0mode) const{
  switch(mode){
  case 1:            return mbc_min_pi0;
  case 10:           return mbc_min_dst0pi0;
  case 2:
    if(h0mode == 10) return mbc_min_etagg;
    else             return mbc_min_etappp;
  case 20:
    if(h0mode == 10) return mbc_min_dst0etagg;
    else             return mbc_min_dst0etappp;
  case 3:            return mbc_min_omega;
  case 5:
    if(h0mode == 10) return mbc_min_etapgg;
    else             return mbc_min_etapppp;
  }
  return -1;
}

double MyParams::get_mbc_max_h0(const int mode,const int h0mode) const{
  switch(mode){
  case 1:            return mbc_max_pi0;
  case 10:           return mbc_max_dst0pi0;
  case 2:
    if(h0mode == 10) return mbc_max_etagg;
    else             return mbc_max_etappp;
  case 20:
    if(h0mode == 10) return mbc_max_dst0etagg;
    else             return mbc_max_dst0etappp;
  case 3:            return mbc_max_omega;
  case 5:
    if(h0mode == 10) return mbc_max_etapgg;
    else             return mbc_max_etapppp;
  }
  return -1;
}

double MyParams::get_mh0(const int mode) const{
  switch(mode){
  case 1:  return Pi0Mass;
  case 10: return Pi0Mass;
  case 3: return OmegaMass;
  case 30: return OmegaMass;
  default: return EtaMass;
  }
  return -1;
}

double MyParams::get_mh0_min(const int mode,const int h0mode) const{
  switch(mode){
  case 1:  return mpi0_min;
  case 10: return mpi0_min;
  case 2:
    if(h0mode == 10) return metagg_min;
    else             return metappp_min;
  case 20:
    if(h0mode == 10) return metagg_min;
    else             return metappp_min;
  case 3: return momega_min;
  case 5: return metagg_min;
  }
  return -1;
}

double MyParams::get_mh0_max(const int mode,const int h0mode) const{
  switch(mode){
  case 1:  return mpi0_max;
  case 10: return mpi0_max;
  case 2:
    if(h0mode == 10) return metagg_max;
    else             return metappp_max;
  case 20:
    if(h0mode == 10) return metagg_max;
    else             return metappp_max;
  case 3: return momega_max;
  case 5: return metagg_max;
  }
  return -1;
}

double MyParams::get_dm_etap_min(const int h0mode) const{
  if(h0mode == 10) return dmetapgg_min;
  return dmetapppp_min;
}

double MyParams::get_dm_etap_max(const int h0mode) const{
  if(h0mode == 10) return dmetapgg_max;
  return dmetapppp_max;
}

double MyParams::get_a_c1_bb_cmb(const int mode,const int h0mode){
  switch(mode){
  case 1: return m_cmb_c10_pi0;
  case 10: return m_cmb_c10_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_c10_eta10;
    else             return m_cmb_c10_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_c10_eta10;
    else             return m_cmb_c10_eta20;
  case 3: return m_cmb_c10_omega;
  case 5: return m_cmb_c10_eta10;
  }
  return 0;
}

double MyParams::get_b_c1_bb_cmb(const int mode,const int h0mode){
  switch(mode){
  case 1: return m_cmb_c11_pi0;
  case 10: return m_cmb_c11_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_c11_eta10;
    else             return m_cmb_c11_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_c11_eta10;
    else             return m_cmb_c11_eta20;
  case 3: return m_cmb_c11_omega;
  case 5: return m_cmb_c11_eta10;
  }
  return 0;
}

double MyParams::get_c2_bb_cmb(const int mode,const int h0mode){
  switch(mode){
  case 1: return m_cmb_c20_pi0;
  case 10: return m_cmb_c20_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_c20_eta10;
    else             return m_cmb_c20_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_c20_eta10;
    else             return m_cmb_c20_eta20;
  case 3: return m_cmb_c20_omega;
  case 5: return m_cmb_c20_eta10;
  }
  return 0;
}

double MyParams::get_c1_qq_cmb(const int mode,const int h0mode){
  switch(mode){
  case 1: return m_cmb_c1_pi0;
  case 10: return m_cmb_c1_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_c1_eta10;
    else             return m_cmb_c1_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_c1_eta10;
    else             return m_cmb_c1_eta20;
  case 3: return m_cmb_c1_omega;
  case 5: return m_cmb_c1_eta10;
  }
  return 0;
}

double MyParams::get_c2_qq_cmb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_c2_pi0;
  case 10: return m_cmb_c2_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_c2_eta10;
    else             return m_cmb_c2_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_c2_eta10;
    else             return m_cmb_c2_eta20;
  case 3: return m_cmb_c2_omega;
  case 5: return m_cmb_c2_eta10;
  }
  return 0;
}

double MyParams::get_argpar_cmb_bb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_argpar_bb_pi0;
  case 10: return -129.5;
  case 2:
    if(h0mode == 10) return m_cmb_argpar_bb_eta10;
    else             return m_cmb_argpar_bb_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_argpar_bb_eta10;
    else             return m_cmb_argpar_bb_eta20;
  case 3: return m_cmb_argpar_bb_omega;
  case 5: return m_cmb_argpar_bb_eta10;
  }
  return 0;
}

double MyParams::get_mbc0_cmb_bb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_mbc0_bb_pi0;
  case 10: return m_cmb_mbc0_bb_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_mbc0_bb_eta10;
    else             return m_cmb_mbc0_bb_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_mbc0_bb_eta10;
    else             return m_cmb_mbc0_bb_eta20;
  case 3: return m_cmb_mbc0_bb_omega;
  case 5: return m_cmb_mbc0_bb_eta10;
  }
  return 0;
}

double MyParams::get_s_mbc_cmb_bb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_mbcw_bb_pi0;
  case 10: return m_cmb_mbcw_bb_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_mbcw_bb_eta10;
    else             return m_cmb_mbcw_bb_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_mbcw_bb_eta10;
    else             return m_cmb_mbcw_bb_eta20;
  case 3: return m_cmb_mbcw_bb_omega;
  case 5: return m_cmb_mbcw_bb_eta10;
  }
  return 0;
}

double MyParams::get_fg_cmb_bb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_f_g_bb_pi0;
  case 10: return 0;
  case 2:
    if(h0mode == 10) return m_cmb_f_g_bb_eta10;
    else             return m_cmb_f_g_bb_eta20;
  case 20:
    if(h0mode == 10) return 0;
    else             return 0;
  case 3: return m_cmb_f_g_bb_omega;
  case 5: return m_cmb_f_g_bb_eta10;
  }
  return 0;
}

double MyParams::get_argpar_cmb_qq(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_argpar_qq_pi0;
  case 10: return m_cmb_argpar_qq_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_argpar_qq_eta10;
    else             return m_cmb_argpar_qq_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_argpar_qq_eta10;
    else             return m_cmb_argpar_qq_eta20;
  case 3: return m_cmb_argpar_qq_omega;
  case 5: return m_cmb_argpar_qq_eta10;
  }
  return 0;
}

double MyParams::get_de0_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_de0r_pi0;
  case 10: return m_peak_de0r_pi0;
  case 2:
    if(h0mode == 10) return m_peak_de0r_eta10;
    else             return m_peak_de0r_eta20;
  case 20:
    if(h0mode == 10) return m_peak_de0r_eta10;
    else             return m_peak_de0r_eta20;
  case 3: return m_peak_de0r_omega;
  case 5: return m_peak_de0r_eta10;
  }
  return 0;
}

double MyParams::get_slopel_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_slopel_pi0;
  case 10: return m_peak_slopel_pi0;
  case 2:
    if(h0mode == 10) return m_peak_slopel_eta10;
    else             return m_peak_slopel_eta20;
  case 20:
    if(h0mode == 10) return m_peak_slopel_eta10;
    else             return m_peak_slopel_eta20;
  case 3: return m_peak_slopel_omega;
  case 5: return m_peak_slopel_eta10;
  }
  return 0;
}

double MyParams::get_sloper_part(const int mode,const int h0mode){
  switch(mode){
  case 1: return m_peak_sloper_pi0;
  case 2:
    if(h0mode == 10) return m_peak_sloper_eta10;
    else             return m_peak_sloper_eta20;
  case 3: return m_peak_sloper_omega;
  case 5: return m_peak_sloper_eta10;
  }
  return 0;
}

double MyParams::get_steep_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_steep_pi0;
  case 10: return m_peak_steep_pi0;
  case 2:
    if(h0mode == 10) return m_peak_steep_eta10;
    else             return m_peak_steep_eta20;
  case 20:
    if(h0mode == 10) return m_peak_steep_eta10;
    else             return m_peak_steep_eta20;
  case 3: return m_peak_steep_omega;
  case 5: return m_peak_steep_eta10;
  }
  return 0;
}

double MyParams::get_p5_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_p5_pi0;
  case 10: return m_peak_p5_pi0;
  case 2:
    if(h0mode == 10) return m_peak_p5_eta10;
    else             return m_peak_p5_eta20;
  case 20:
    if(h0mode == 10) return m_peak_p5_eta10;
    else             return m_peak_p5_eta20;
  case 3: return m_peak_p5_omega;
  case 5: return m_peak_p5_eta10;
  }
  return 0;
}

double MyParams::get_b_s_mbc_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_b_s_pi0;
  case 10: return m_peak_b_s_pi0;
  case 2:
    if(h0mode == 10) return m_peak_b_s_eta10;
    else             return m_peak_s_eta20;
  case 20:
    if(h0mode == 10) return m_peak_b_s_eta10;
    else             return m_peak_s_eta20;
  case 3: return m_peak_s_omega;
  case 5: return m_peak_b_s_eta10;
  }
  return 0;
}

double MyParams::get_k_s_mbc_part(const int mode){
  switch(mode){
  case 1:  return m_peak_k_s_pi0;
  case 10: return m_peak_k_s_pi0;
  case 2:  return m_peak_k_s_eta10;
  case 20: return m_peak_k_s_eta10;
  case 5:  return m_peak_k_s_eta10;
  }
  return 0;
}

double MyParams::get_b_mbc0_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_b_mbc0_pi0;
  case 10: return m_peak_b_mbc0_pi0;
  case 2:
    if(h0mode == 10) return m_peak_b_mbc0_eta10;
    else             return m_peak_mbc0_eta20;
  case 20:
    if(h0mode == 10) return m_peak_b_mbc0_eta10;
    else             return m_peak_mbc0_eta20;
  case 3: return m_peak_mbc0_omega;
  case 5: return m_peak_b_mbc0_eta10;
  }
  return 0;
}

double MyParams::get_k_mbc0_part(const int mode){
  switch(mode){
  case 1:  return m_peak_k_mbc0_pi0;
  case 10: return m_peak_k_mbc0_pi0;
  case 2:  return m_peak_k_mbc0_eta10;
  case 20: return m_peak_k_mbc0_eta10;
  case 5:  return m_peak_k_mbc0_eta10;
  }
  return 0;
}

double MyParams::get_argedge(const int mode,const int h0mode){
  switch(mode){
  case 1:  return 5.28814e+00;
  case 10: return 5.28814e+00;
  case 2:
    if(h0mode == 10) return 5.28814e+00;
    else             return 5.28814e+00;
  case 20:
    if(h0mode == 10) return 5.28814e+00;
    else             return 5.28814e+00;
  case 3: return 5.28814e+00;
  case 5: return 5.28814e+00;
  }
  return 0;
}

double MyParams::get_argpar_part_bb(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_argpar_bb_pi0;
  case 10: return m_cmb_argpar_bb_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_argpar_bb_eta10;
    else             return m_cmb_argpar_bb_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_argpar_bb_eta10;
    else             return m_cmb_argpar_bb_eta20;
  case 3: return m_cmb_argpar_bb_omega;
  case 5: return m_cmb_argpar_bb_eta10;
  }
  return 0;
}

double MyParams::get_mbc0_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_b_mbc0_pi0;
  case 10: return m_peak_b_mbc0_pi0;
  case 2:
    if(h0mode == 10) return m_peak_b_mbc0_eta10;
    else             return m_peak_mbc0_eta20;
  case 20:
    if(h0mode == 10) return m_peak_b_mbc0_eta10;
    else             return m_peak_mbc0_eta20;
  case 3: return m_peak_mbc0_omega;
  case 5: return m_peak_b_mbc0_eta10;
  }
  return 0;
}

double MyParams::get_s_mbc_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_peak_b_s_pi0;
  case 10: return m_peak_b_s_pi0;
  case 2:
    if(h0mode == 10) return m_peak_b_s_eta10;
    else             return m_peak_s_eta20;
  case 20:
    if(h0mode == 10) return m_peak_b_s_eta10;
    else             return m_peak_s_eta20;
  case 3: return m_peak_s_omega;
  case 5: return m_peak_b_s_eta10;
  }
  return 0;
}

double MyParams::get_fg_mbc_part(const int mode,const int h0mode){
  switch(mode){
  case 1:  return m_cmb_f_g_bb_pi0;
  case 10: return m_cmb_f_g_bb_pi0;
  case 2:
    if(h0mode == 10) return m_cmb_f_g_bb_eta10;
    else             return m_cmb_f_g_bb_eta20;
  case 20:
    if(h0mode == 10) return m_cmb_f_g_bb_eta10;
    else             return m_cmb_f_g_bb_eta20;
  case 3: return m_cmb_f_g_bb_omega;
  case 5: return m_cmb_f_g_bb_eta10;
  }
  return 0;
}

double MyParams::f_p_f_bbc(const int mode,const int h0mode){
  switch(mode){
  case 1:  return 4.67872134323538913e-01;//0.0051;
  case 10: return 2.78008298755186734e-01;//5.39867109634551492e-02;//0.0051;
  case 2:
    if(h0mode == 10) return 5.39867109634551492e-02;//0.0051;
    else             return 1.03274559193954660e-01;//0.0081;
  case 20:
    if(h0mode == 10) return 7.06521739130434728e-02;//0.0051;
    else             return 1.03274559193954660e-01;//0.0081;
  case 3: return 4.23086329175837647e-02;//0.0031;
  case 5: return 8.12720848056537049e-02;//0.0051;
  }
  return 0;
}

string MyParams::GetLabel(const int mode,const int h0mode){
  switch(mode){
  case 1:
    return string("#pi^{0}");
  case 10:
    return string("D^{*0}#pi^{0}");
  case 2:
    if(h0mode == 10) return string("#eta#rightarrow#gamma#gamma");
    else             return string("#eta#rightarrow#pi^{+}#pi^{-}#pi^{0}");
  case 20:
    if(h0mode == 10) return string("D^{*0}#eta(#rightarrow#gamma#gamma)");
    else             return string("D^{*0}#eta(#rightarrow#pi^{+}#pi^{-}#pi^{0})");
  case 3:
    return string("#omega");
  case 5:
    return string("#eta`");
  default:
    break;
  }
  return string("");
}

int MyParams::bin(const int j){
  if(j<8) return j-8;
  else    return j-7;
}
int MyParams::flv_ind(const int flv){
  if(flv == 1) return 0;
  else         return 1;
}
int MyParams::bin_ind(const int bin){
  if(bin<0) return bin+8;
  else      return bin+7;
}
int MyParams::flv(const int k){
  if(k) return -1;
  else  return  1;
}

double MyParams::K(const int bin){
  return bin>0 ? Karr[abs(bin)-1] : Kbarr[abs(bin)-1];
}

//inline double MyParams::N(const int bin, const int flv){
//  return ((K(bin) + K(-bin))*m_xi + flv*(K(bin) - K(-bin)))/(2*m_xi);
//}

double MyParams::N(const int bin, const int flv, const double& wt){
  return ((K(bin) + K(-bin))*m_xi + (1-2*wt)*flv*(K(bin) - K(-bin)))/(2*m_xi);
}

int MyParams::get_wbin(const double& tag){
  const double atag = TMath::Abs(tag);
  for(int i=0; i<6; i++){ if(atag<wbins[i+1]){ return i;}}
  return 6;
}

double MyParams::get_wtag_prob(const double& wt, const int exp,const bool data){
  int wbin = get_wbin(wt);
  double wprob;
  if(data){wprob = exp > 30 ? w_data_svd2[wbin] : w_data_svd1[wbin];}
  else    {wprob = exp > 30 ? w_mc_svd2[wbin]   : w_mc_svd1[wbin];}
  return wprob;
}

bool MyParams::IsInEllips(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde){
  const double xi_de = (de-de0)/Rde;
  const double xi_mbc = (mbc-mbc0)/Rmbc;
  if(xi_de*xi_de + xi_mbc*xi_mbc < 1) return true;
  else return false;
}

double MyParams::EllipsR2(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde){
  const double xi_de = (de-de0)/Rde;
  const double xi_mbc = (mbc-mbc0)/Rmbc;
  return xi_de*xi_de + xi_mbc*xi_mbc;
}

int MyParams::q(const double& tag_LH){
  if(tag_LH<0) return 1;
   return -1;
}

double MyParams::h0mass(const int mode){
  switch(mode){
  case 1:  return Pi0Mass;
  case 10: return Pi0Mass;
  case 3:  return OmegaMass;
  case 30: return OmegaMass;
//  case 4:  return RhoMass;
  default: return EtaMass;
  }
  return 0;
}

double MyParams::get_f_cont_in_bkg_sig(const int mode, const int h0mode){
  switch(mode){
  case 1:  return 439./(439+226+152);
  case 2:
    if(h0mode == 10) return 189./(189+223+9);
    else             return 41./(41+19+2);
  case 3:  return 262./(262+308+11);
  case 5:  return 35./(35+31+3);
  case 10: return 87./(87+327+39);
  case 20: return 7./(7+48);
  }
  return 0;
}
double MyParams::get_f_cont_in_bkg_sideband(const int mode, const int h0mode){
  switch(mode){
  case 1:  return 3955./(3955+1209+21);
  case 2:
    if(h0mode == 10) return 2411./(2411+353+0);
    else             return 1102./(1102+142+5);
  case 3:  return 5548./(5548+1842+12);
  case 5:  return 446./(446+115);
  case 10: return 901./(901+405+6);
  case 20: return 165./(165+88);
  }
  return 0;
}
