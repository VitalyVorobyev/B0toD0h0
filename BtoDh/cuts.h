#ifndef CUTS_H
#define CUTS_H

#include "shapes_pi0.h"
#include "shapes_eta.h"
#include "shapes_omega.h"

int bin(const int j){
  if(j<8) return j-8;
  else    return j-7;
}
int flv_ind(const int flv){
  if(flv == 1) return 0;
  else         return 1;
}
int bin_ind(const int bin){
  if(bin<0) return bin+8;
  else      return bin+7;
}
int flv(const int k){
  if(k) return -1;
  else  return  1;
}

const double Karr[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
const double Kbarr[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};
inline double K(const int bin){
  return bin>0 ? Karr[abs(bin)-1] : Kbarr[abs(bin)-1];
}
const double m_btau     = 1.534;//1.520;//0.460/sol*cm2ps;//0.507;// ps
const double m_dm       = 0.510;//0.452;//0.510/sol;//0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
const double m_xi = 1+m_btau*m_dm*m_btau*m_dm;

inline double N(const int bin, const int flv){
  return ((K(bin) + K(-bin))*m_xi + flv*(K(bin) - K(-bin)))/(2*m_xi);
}

const double wbins[8]        = {0.0,0.1,0.25,0.5,0.625,0.75,0.875,1.01};
const double w_mc_svd1[7]    = {0.5,0.420827,0.300296,0.219317,0.154636,0.0916131,0.0228891};
//const double w_data_svd1[7]  = {0.5,0.418852,0.329879,0.233898,0.170608,0.099791, 0.0228501};
const double w_mc_svd2[7]    = {0.5,0.412222,0.307838,0.212765,0.149933,0.0913264,0.0218754};
//const double w_data_svd2[7]  = {0.5,0.418826,0.319303,0.222948,0.163191,0.104085, 0.0251454};

int get_wbin(const double& tag){
  const double atag = TMath::Abs(tag);
//  cout << atag << " ";
  for(int i=0; i<6; i++){
    if(atag<wbins[i+1]){
//      cout << "return bin " << i << " ";
      return i;
    }
  }
//  cout << "return bin " << i << " ";
  return 6;
}

double get_wtag_prob(const double& wt, const int exp){
  int wbin = get_wbin(wt);
//  cout << wt << " " << wbin << " " << w_mc_svd2[wbin] << " " << w_mc_svd1[wbin] << " ";
  return exp > 30 ? w_mc_svd2[wbin] : w_mc_svd1[wbin];
}

double N(const int bin, const int flv, const double& wt){
  return ((K(bin) + K(-bin))*m_xi + (1-2*wt)*flv*(K(bin) - K(-bin)))/(2*m_xi);
}

const bool fullfit = false;

const bool const_back_flag = false;
const bool flag41 = false;

const double de_fit_min = -0.15;
const double de_fit_max = 0.3;

const double mbc_fit_min = 5.20;
const double mbc_fit_max = 5.29;

const double de_min_pi0 = -0.128342;
const double de_max_pi0 = 0.0921343;

const double de_min_etagg = -0.0994528;
const double de_max_etagg = 0.0758825;

const double de_min_etappp = -0.0534644;
const double de_max_etappp = 0.0453713;

const double de_min_etapgg = -0.0873898;
const double de_max_etapgg =  0.0642064;

const double de_min_etapppp = -0.0394845;
const double de_max_etapppp =  0.0357815;

const double de_min_omega = -0.0565915;
const double de_max_omega = 0.0467748;

//const double mbc_min = 5.272;
//const double mbc_max = 5.286;

const double mbc_min_pi0 = 5.27093;
const double mbc_max_pi0 = 5.28789;

const double mbc_min_etagg = 5.27147;
const double mbc_max_etagg = 5.28763;

const double mbc_min_etappp = 5.27205;
const double mbc_max_etappp = 5.28735;

const double mbc_min_etapgg = 5.27255;
const double mbc_max_etapgg = 5.28718;

const double mbc_min_etapppp = 5.27288;
const double mbc_max_etapppp = 5.28702;

const double mbc_min_omega = 5.27198;
const double mbc_max_omega = 5.28743;

const double mbc_side = 5.25;

//const double mbc_min_omega = 5.273;
//const double mbc_max_omega = 5.287;

const double DMass = 1.865;
const double md_cut = 0.013;

const double md_min = 1.85155;
const double md_max = 1.87833;

const double KMass = 0.4975;
const double mk_cut = 0.009;

const double Pi0Mass = 0.135;
const double mpi0_cut = 0.012;

const double mpi0_min = 0.115738;
const double mpi0_max = 0.153714;

const double EtaMass = 0.5478;
const double EtaGGMass = 0.541;
const double metagg_cut = 0.035;
const double metappp_cut = 0.015;

const double metagg_min = 0.517603;
const double metagg_max = 0.573654;

const double metappp_min = 0.537626;
const double metappp_max = 0.557443;

const double dmetapgg_min = 0.401668;
const double dmetapgg_max = 0.417717;

const double dmetapppp_min = 0.402493;
const double dmetapppp_max = 0.417109;

const double OmegaMass = 0.781;
const double momega_cut = 0.03;

const double momega_min = 0.760430;
const double momega_max = 0.803869;

const double RhoMass = 0.7745;
const double mrho_cut = 0.150;

const double dmdst0_min = 0.149 - 0.007;
const double dmdst0_max = 0.149 + 0.007;

double get_mode(const int type){
  switch(type){
  case 1:  return 1; // pi0
  case 10: return 10;// D*0 pi0
  case 2:  return 2; // eta -> gg
  case 20: return 20;// D*0 eta
  case 3:  return 2; // eta -> pi+pi-pi0
  case 4:  return 3; // omega
  case 5:  return 5; // eta'
  }
  return 0;
}

double get_h0mode(const int type){
  switch(type){
  case 1:  return 10;// pi0
  case 10: return 10;// D*0 pi0
  case 2:  return 10;// eta -> gg
  case 20: return 10;// D*0 eta
  case 3:  return 20;// eta -> pi+pi-pi0
  case 4:  return 20;// omega
  case 5:  return 10;// eta'
  }
  return 0;
}

double h0mass(const int mode){
  switch(mode){
  case 1:  return Pi0Mass;
  case 10: return Pi0Mass;
  case 3:  return OmegaMass;
  case 30: return OmegaMass;
  case 4:  return RhoMass;
  default: return EtaMass;
  }
  return 0;
}

int get_mode_number(const int type){
  switch(type){
  case 1:  return 1;// pi0
  case 10: return 1;// D*0 pi0
  case 2:  return 2;// eta -> gg
  case 20: return 2;// D*0 eta
  case 3:  return 3;// eta -> pi+pi-pi0
  case 4:  return 4;// omega
  case 5:  return 5;// eta'
  }
}

double h0_min_mass(const int mode,const int h0mode=10){
  switch(mode){
  case 1:  return mpi0_min;
  case 10: return mpi0_min;
  case 3:  return momega_min;
  case 30: return momega_min;
  case 4:  return RhoMass-mrho_cut;
  default:
    if(h0mode == 10) return metagg_min;
    else             return metappp_min;
  }
  return 0;
}

double h0_max_mass(const int mode,const int h0mode=10){
  switch(mode){
  case 1:  return mpi0_max;
  case 10: return mpi0_max;
  case 3:  return momega_max;
  case 30: return momega_max;
  case 4:  return RhoMass+mrho_cut;
  default:
    if(h0mode == 10) return metagg_max;
    else             return metappp_max;
  }
  return 0;
}

const double atckpi_cut = 1.;

double bdt_cut_pi0    = 0.26;
double bdt_cut_etagg  = 0.21;
double bdt_cut_etappp = 0.24;
double bdt_cut_omega  = 0.25;
//double bdtgs_cut = 0.916515;//0.98;
//double bdtgfr_cut = 0.2168;
//double bdtgsl_cut = 0.821645;
//double bdtgmbcs_cut = 0.759934;
//double bdtgdes_cut = 0.6;
double lh_cut = 0.814616;
double lh1_cut = 0.742361;
const double bdtg_cut_pi0    = 0.910;
const double bdtg_cut_etagg  = 0.690;
const double bdtg_cut_etappp = 0.500;
//const double bdtg_cut_omega  = 0.125;
const double bdtg_cut_omega  = 0.8;

double bdt_cut(const int mode,const int h0mode=10){
  switch(mode){
  case 1:  return bdt_cut_pi0;
  case 10: return bdt_cut_pi0;
  case 3:  return bdt_cut_omega;
  case 30: return bdt_cut_omega;
  case 4:  return bdt_cut_omega;// rho
  default:
    if(h0mode == 10) return bdt_cut_etagg;
    else             return bdt_cut_etappp;
  }
  return 0;
}

// Signal shape
const bool cSig = true;
const bool cSIG = true;
// mbc
const double m_mbc0 = 5.28291e+00;
const double m_mbc00 = 5.27979e+00;
const double m_sl = 6.28487e-03;
const double m_sll = 2.96898e-03;
const double m_sr = 1.92172e-03;
const double m_srr = 2.17430e-03;
const double m_fmbc = 2.80705e-01;

// de
const double m_alphal = 3.85675e-01;
const double m_alphar = -1.56269e+00;
const double m_de0 = 1.51293e-02;
const double m_deCBl = -2.58588e-02;
const double m_deCBr = -1.22991e-02;
const double m_fCBl = 4.62725e-01;
const double m_fCBr = 2.11382e-01;
const double m_nl = 1.71068e+01;
const double m_nr = 3.79904e+00;
const double m_s1 = 2.24893e-02;
const double m_sCBl = 2.71645e-02;
const double m_sCBr = 3.97215e-02;

// Rho shape
const bool keysflag = false;
const bool cRho = true;
const bool cRHO = false;
// mbc
const int mbc_rho_param = 4;
const double mr_mbc0_0 = 5.28177e+00;
const double mr_sl_0 = 1.88012e-02;
const double mr_sr_0 = 3.21321e-03;
const double mr_cond_0 = 8.97366e-02;
const double mr_condr_0 = 8.40921e-03;

const double mr_mbc0_1 = 5.28064e+00;
const double mr_sl_1 = 2.76386e-02;
const double mr_sr_1 = 1.04874e-06;
const double mr_cond_1 = 0.006;
const double mr_mbc00_1 = 5.27999e+00;
const double mr_sll_1 = 6.24356e-03;
const double mr_srr_1 = 4.45965e-03;
const double mr_fmbc_1 = 5.04927e-02;

const double mr_mbc0_2 = 5.27878e+00;
const double mr_sr_2 = 3.80005e-03;
const double mr_mbc1_2 = 5.26645e+00;
const double mr_fcb_2 = 6.11183e-02;
const double mr_alpha_2 = -3.74609e+00;
const double mr_ss_2 = 1.15691e-02;
const double mr_n_2 = 0.48;

const double mr_argpar_3 = -9.53674e-05;
const double mr_argedge_3 = 5.28616e+00;
const double mr_mbc0_3 = 5.28440e+00;
const double mr_cond_3 = 3.97230e-02;
const double mr_condr_3 = 3.53054e-02;
const double mr_sl_3 = 1.21988e-02;
const double mr_sr_3 = 3.16735e-03;
//de
const int de_rho_param = 1;

const double mr_exppar = -3.37184e+01;

const double mr_de0r = -1.11131e-01;
const double mr_slopel = -7.03669e+02;
const double mr_sloper = -1.52863e+00;
const double mr_steep = 9.52444e+00;
const double mr_p5 = 9.39841e-01;

const double mr_x0 = -1.27853e-01;
const double mr_p1 = -2.80709e+02;
const double mr_p2 = 5.31945e+00;

const double mr_x0_1d = -1.21409e-01;
const double mr_p1_1d = -4.26629e+02;
const double mr_p2_1d = 7.27530e+00;

const double mr_x0_2 = -5.48334e-02;
const double mr_slopel_2 = 1.29144e+01;
const double mr_steep_2 = 2.47426e+03;
const double mr_exppar_2 = 4.38643e+01;

// Comb shape
const bool cComb = true;
const bool cPeak = true;
// mbc
//const double mc_argpar = -2.37807e+01;
//const double mc_argedge = 5.28941e+00;
const double mc_argpar = -2.04306e+01;
const double mc_argedge = 5.28941e+00;

// de
//const double mc_c1 = -4.57562e-01;
//const double mc_c2 = 6.94555e-02;
const double mc_c1 = -4.55023e-01;
const double mc_c2 = 6.71316e-02;

const double mc_c1_1d = -5.39318e-01;
const double mc_c2_1d = 9.53054e-02;

// RooKSFW pre cuts
const double k0mm2_min = -20.5;
const double k0mm2_max = 22.;

const double k0et_min = 2.8;
const double k0et_max = 12.;

const double k0hso00_min = 0.;
const double k0hso00_max = 1.3;

const double k0hso02_min = -0.5;
const double k0hso02_max = 0.75;

//const double k0hso4_min = -0.35;
//const double k0hso4_max = 0.6;

const double k0hso10_min = 0;
const double k0hso10_max = 1.;

const double k0hso12_min = -0.35;
const double k0hso12_max = 0.6;

const double k0hso14_min = -0.3;
const double k0hso14_max = 0.45;

const double k0hso20_min = 0.;
const double k0hso20_max = 0.75;

const double k0hso22_min = -0.3;
const double k0hso22_max = 0.6;

const double k0hso24_min = -0.3;
const double k0hso24_max = 0.5;

const double k0hoo0_min = 0.;
const double k0hoo0_max = 0.5;

const double k0hoo1_min = -0.05;
const double k0hoo1_max = 0.07;

const double k0hoo2_min = -0.03;
const double k0hoo2_max = 0.02;

const double k0hoo3_min = -0.04;
const double k0hoo3_max = 0.05;

const double k0hoo4_min = -0.035;
const double k0hoo4_max = 0.150;

const double k1mm2_min = -20.5;
const double k1mm2_max = 22.;

const double k1et_min = 3.7;
const double k1et_max = 13.;

const double k1hso00_min = 0.;
const double k1hso00_max = 3.8;

const double k1hso01_min = -0.8;
const double k1hso01_max = 0.8;

const double k1hso02_min = -1.;
const double k1hso02_max = 1.6;

const double k0hso4_min = -0.75;
const double k0hso4_max = 1.;

const double k1hso10_min = 0;
const double k1hso10_max = 3.;

const double k1hso12_min = -0.75;
const double k1hso12_max = 1.;

const double k1hso14_min = -0.5;
const double k1hso14_max = 0.7;

const double k1hso20_min = 0.;
const double k1hso20_max = 2.3;

const double k1hso22_min = -0.7;
const double k1hso22_max = 1.2;

const double k1hso24_min = -0.6;
const double k1hso24_max = 0.9;

const double k1hoo0_min = 0.;
const double k1hoo0_max = 0.5;

const double k1hoo1_min = -0.05;
const double k1hoo1_max = 0.07;

const double k1hoo2_min = -0.03;
const double k1hoo2_max = 0.02;

const double k1hoo3_min = -0.04;
const double k1hoo3_max = 0.05;

const double k1hoo4_min = -0.035;
const double k1hoo4_max = 0.150;

// de
double get_alphal(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_alphal_pi1;
    else if(b0f == 5) return m_alphal_pi5;
    else              return m_alphal_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_alphal_eta10;
      else if(b0f == 5) return m_alphal_eta105;
      else              return m_alphal_eta10m1;
    } else{
      if(b0f == 1)      return m_alphal_eta201;
      else if(b0f == 5) return m_alphal_eta205;
      else              return m_alphal_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_alphal_omega201;
    else if(b0f == 5) return m_alphal_omega205;
    else              return m_alphal_omega20m1;
  }
  return -1;
}
double get_alphar(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_alphar_pi1;
    else if(b0f == 5) return m_alphar_pi5;
    else              return m_alphar_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_alphar_eta10;
      else if(b0f == 5) return m_alphar_eta105;
      else              return m_alphar_eta10m1;
    } else{
      if(b0f == 1)      return m_alphar_eta201;
      else if(b0f == 5) return m_alphar_eta205;
      else              return m_alphar_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_alphar_omega201;
    else if(b0f == 5) return m_alphar_omega205;
    else              return m_alphar_omega20m1;
  }
  return -1;
}
double get_de0(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_de0_pi1;
    else if(b0f == 5) return m_de0_pi5;
    else              return m_de0_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_de0_eta10;
      else if(b0f == 5) return m_de0_eta105;
      else              return m_de0_eta10m1;
    } else{
      if(b0f == 1)      return m_de0_eta201;
      else if(b0f == 5) return m_de0_eta205;
      else              return m_de0_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_de0_omega201;
    else if(b0f == 5) return m_de0_omega205;
    else              return m_de0_omega20m1;
  }
  return -1;
}
double get_deCBl(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_deCBl_pi1;
    else if(b0f == 5) return m_deCBl_pi5;
    else              return m_deCBl_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_deCBl_eta10;
      else if(b0f == 5) return m_deCBl_eta105;
      else              return m_deCBl_eta10m1;
    } else{
      if(b0f == 1)      return m_deCBl_eta201;
      else if(b0f == 5) return m_deCBl_eta205;
      else              return m_deCBl_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_deCBl_omega201;
    else if(b0f == 5) return m_deCBl_omega205;
    else              return m_deCBl_omega20m1;
  }
  return -1;
}
double get_deCBr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_deCBr_pi1;
    else if(b0f == 5) return m_deCBr_pi5;
    else              return m_deCBr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_deCBr_eta10;
      else if(b0f == 5) return m_deCBr_eta105;
      else              return m_deCBr_eta10m1;
    } else{
      if(b0f == 1)      return m_deCBr_eta201;
      else if(b0f == 5) return m_deCBr_eta205;
      else              return m_deCBr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_deCBr_omega201;
    else if(b0f == 5) return m_deCBr_omega205;
    else              return m_deCBr_omega20m1;
  }
  return -1;
}
double get_fCBl(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_fCBl_pi1;
    else if(b0f == 5) return m_fCBl_pi5;
    else              return m_fCBl_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_fCBl_eta10;
      else if(b0f == 5) return m_fCBl_eta105;
      else              return m_fCBl_eta10m1;
    } else{
      if(b0f == 1)      return m_fCBl_eta201;
      else if(b0f == 5) return m_fCBl_eta205;
      else              return m_fCBl_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_fCBl_omega201;
    else if(b0f == 5) return m_fCBl_omega205;
    else              return m_fCBl_omega20m1;
  }
  return -1;
}
double get_fCBr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_fCBr_pi1;
    else if(b0f == 5) return m_fCBr_pi5;
    else              return m_fCBr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_fCBr_eta10;
      else if(b0f == 5) return m_fCBr_eta105;
      else              return m_fCBr_eta10m1;
    } else{
      if(b0f == 1)      return m_fCBr_eta201;
      else if(b0f == 5) return m_fCBr_eta205;
      else              return m_fCBr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_fCBr_omega201;
    else if(b0f == 5) return m_fCBr_omega205;
    else              return m_fCBr_omega20m1;
  }
  return -1;
}
double get_nl(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_nl_pi1;
    else if(b0f == 5) return m_nl_pi5;
    else              return m_nl_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_nl_eta10;
      else if(b0f == 5) return m_nl_eta105;
      else              return m_nl_eta10m1;
    } else{
      if(b0f == 1)      return m_nl_eta201;
      else if(b0f == 5) return m_nl_eta205;
      else              return m_nl_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_nl_omega201;
    else if(b0f == 5) return m_nl_omega205;
    else              return m_nl_omega20m1;
  }
  return -1;
}
double get_nr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_nr_pi1;
    else if(b0f == 5) return m_nr_pi5;
    else              return m_nr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_nr_eta10;
      else if(b0f == 5) return m_nr_eta105;
      else              return m_nr_eta10m1;
    } else{
      if(b0f == 1)      return m_nr_eta201;
      else if(b0f == 5) return m_nr_eta205;
      else              return m_nr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_nr_omega201;
    else if(b0f == 5) return m_nr_omega205;
    else              return m_nr_omega20m1;
  }
  return -1;
}
double get_s1(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_s1_pi1;
    else if(b0f == 5) return m_s1_pi5;
    else              return m_s1_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_s1_eta10;
      else if(b0f == 5) return m_s1_eta105;
      else              return m_s1_eta10m1;
    } else{
      if(b0f == 1)      return m_s1_eta201;
      else if(b0f == 5) return m_s1_eta205;
      else              return m_s1_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_s1_omega201;
    else if(b0f == 5) return m_s1_omega205;
    else              return m_s1_omega20m1;
  }
  return -1;
}
double get_sCBl(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_sCBl_pi1;
    else if(b0f == 5) return m_sCBl_pi5;
    else              return m_sCBl_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_sCBl_eta10;
      else if(b0f == 5) return m_sCBl_eta105;
      else              return m_sCBl_eta10m1;
    } else{
      if(b0f == 1)      return m_sCBl_eta201;
      else if(b0f == 5) return m_sCBl_eta205;
      else              return m_sCBl_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_sCBl_omega201;
    else if(b0f == 5) return m_sCBl_omega205;
    else              return m_sCBl_omega20m1;
  }
  return -1;
}
double get_sCBr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_sCBr_pi1;
    else if(b0f == 5) return m_sCBr_pi5;
    else return m_sCBr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_sCBr_eta10;
      else if(b0f == 5) return m_sCBr_eta105;
      else              return m_sCBr_eta10m1;
    } else{
      if(b0f == 1)      return m_sCBr_eta201;
      else if(b0f == 5) return m_sCBr_eta205;
      else              return m_sCBr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_sCBr_omega201;
    else if(b0f == 5) return m_sCBr_omega205;
    else              return m_sCBr_omega20m1;
  }
  return -1;
}
double get_mbc0(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_mbc0_pi1;
    else if(b0f == 5) return m_mbc0_pi5;
    else              return m_mbc0_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_mbc0_eta10;
      else if(b0f == 5) return m_mbc0_eta105;
      else              return m_mbc0_eta10m1;
    } else{
      if(b0f == 1)      return m_mbc0_eta201;
      else if(b0f == 5) return m_mbc0_eta205;
      else              return m_mbc0_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_mbc0_omega201;
    else if(b0f == 5) return m_mbc0_omega205;
    else              return m_mbc0_omega20m1;
  }
  return -1;
}
double get_mbc00(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_mbc00_pi1;
    else if(b0f == 5) return m_mbc00_pi5;
    else              return m_mbc00_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_mbc00_eta10;
      else if(b0f == 5) return m_mbc00_eta105;
      else              return m_mbc00_eta10m1;
    } else{
      if(b0f == 1)      return m_mbc00_eta201;
      else if(b0f == 5) return m_mbc00_eta205;
      else              return m_mbc00_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_mbc00_omega201;
    else if(b0f == 5) return m_mbc00_omega205;
    else              return m_mbc00_omega20m1;
  }
  return -1;
}
double get_sl(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_sl_pi1;
    else if(b0f == 5) return m_sl_pi5;
    else              return m_sl_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_sl_eta10;
      else if(b0f == 5) return m_sl_eta105;
      else              return m_sl_eta10m1;
    } else{
      if(b0f == 1)      return m_sl_eta201;
      else if(b0f == 5) return m_sl_eta205;
      else              return m_sl_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_sl_omega201;
    else if(b0f == 5) return m_sl_omega205;
    else              return m_sl_omega20m1;
  }
  return -1;
}
double get_sll(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_sll_pi1;
    else if(b0f == 5) return m_sll_pi5;
    else              return m_sll_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_sll_eta10;
      else if(b0f == 5) return m_sll_eta105;
      else              return m_sll_eta10m1;
    } else{
      if(b0f == 1)      return m_sll_eta201;
      else if(b0f == 5) return m_sll_eta205;
      else              return m_sll_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_sll_omega201;
    else if(b0f == 5) return m_sll_omega205;
    else              return m_sll_omega20m1;
  }
  return -1;
}
double get_sr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_sr_pi1;
    else if(b0f == 5) return m_sr_pi5;
    else              return m_sr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_sr_eta10;
      else if(b0f == 5) return m_sr_eta105;
      else              return m_sr_eta10m1;
    } else{
      if(b0f == 1)      return m_sr_eta201;
      else if(b0f == 5) return m_sr_eta205;
      else              return m_sr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_sr_omega201;
    else if(b0f == 5) return m_sr_omega205;
    else              return m_sr_omega20m1;
  }
  return -1;
}
double get_srr(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_srr_pi1;
    else if(b0f == 5) return m_srr_pi5;
    else              return m_srr_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_srr_eta10;
      else if(b0f == 5) return m_srr_eta105;
      else              return m_srr_eta10m1;
    } else{
      if(b0f == 1)      return m_srr_eta201;
      else if(b0f == 5) return m_srr_eta205;
      else              return m_srr_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_srr_omega201;
    else if(b0f == 5) return m_srr_omega205;
    else              return m_srr_omega20m1;
  }
  return -1;
}
double get_fmbc(const int mode,const int h0mode, const int b0f=1){
  if(mode == 1){
    if(b0f == 1)      return m_fmbc_pi1;
    else if(b0f == 5) return m_fmbc_pi5;
    else              return m_fmbc_pi1m;
  }
  if(mode == 2){
    if(h0mode == 10){
      if(b0f == 1)      return m_fmbc_eta10;
      else if(b0f == 5) return m_fmbc_eta105;
      else              return m_fmbc_eta10m1;
    } else{
      if(b0f == 1)      return m_fmbc_eta201;
      else if(b0f == 5) return m_fmbc_eta205;
      else              return m_fmbc_eta20m1;
    }
  }
  if(mode == 3){
    if(b0f == 1)      return m_fmbc_omega201;
    else if(b0f == 5) return m_fmbc_omega205;
    else              return m_fmbc_omega20m1;
  }
  return -1;
}
double get_f201(const int mode,const int h0mode){
  if(mode == 2) return m_f201_eta20;
  else          return m_f201_omega;
}

double get_a_mbc0(const int mode){
  switch(mode){
  case 1: return m_a_mbc0_pi1m;
  case 2: return m_a_mbc0_eta10m1;
  }
  return 0;
}

double get_b_mbc0(const int mode){
  switch(mode){
  case 1: return m_b_mbc0_pi1m;
  case 2: return m_b_mbc0_eta10m1;
  }
  return 0;
}

double get_c_mbc0(const int mode){
  switch(mode){
  case 1: return m_c_mbc0_pi1m;
  case 2: return m_c_mbc0_eta10m1;
  }
  return 0;
}

double get_a5_mbc0(const int mode){
  switch(mode){
  case 3: return m_a_mbc0_eta205;
  case 4: return m_a_mbc0_omega5;
  }
  return 0;
}

double get_b5_mbc0(const int mode){
  switch(mode){
  case 3: return m_b_mbc0_eta205;
  case 4: return m_b_mbc0_omega5;
  }
  return 0;
}

double get_c5_mbc0(const int mode){
  switch(mode){
  case 3: return m_c_mbc0_eta205;
  case 4: return m_c_mbc0_omega5;
  }
  return 0;
}

double get_a_s(const int mode){
  switch(mode){
  case 1: return m_a_s_pi1m;
  case 2: return m_a_s_eta10m1;
  case 3: return m_a_s_eta201;
  case 4: return m_a_s_omega1;
  }
  return 0;
}

double get_b_s(const int mode){
  switch(mode){
  case 1: return m_b_s_pi1m;
  case 2: return m_b_s_eta10m1;
  case 3: return m_b_s_eta201;
  case 4: return m_b_s_omega1;
  }
  return 0;
}

double get_c_s(const int mode){
  switch(mode){
  case 1: return m_c_s_pi1m;
  case 2: return m_c_s_eta10m1;
  case 3: return m_c_s_eta201;
  case 4: return m_c_s_omega1;
  }
}

double get_a5_s(const int mode){
  switch(mode){
  case 3: return m_a_s_eta205;
  case 4: return m_a_s_omega5;
  }
  return 0;
}

double get_b5_s(const int mode){
  switch(mode){
  case 3: return m_b_s_eta205;
  case 4: return m_b_s_omega5;
  }
  return 0;
}

double get_c5_s(const int mode){
  switch(mode){
  case 3: return m_c_s_eta205;
  case 4: return m_c_s_omega5;
  }
  return 0;
}

double get_c0(const int mode){
  switch(mode){
  case 3: return m_c0_eta201;
  case 4: return m_c0_omega1;
  }
  return 0;
}

double get_c1(const int mode){
  switch(mode){
  case 3: return m_c1_eta201;
  case 4: return m_c1_omega1;
  }
  return 0;
}

double get_c2(const int mode){
  switch(mode){
  case 3: return m_c2_eta201;
  case 4: return m_c2_omega1;
  }
  return 0;
}

double get_cmb_c10(const int mode){
  switch(mode){
  case 1: return m_cmb_c10_pi0;
  case 2: return m_cmb_c10_eta10;
  case 3: return m_cmb_c10_eta20;
  case 4: return m_cmb_c10_omega;
  }
  return 0;
}

double get_cmb_c11(const int mode){
  switch(mode){
  case 1: return m_cmb_c11_pi0;
  case 2: return m_cmb_c11_eta10;
  case 3: return m_cmb_c11_eta20;
  case 4: return m_cmb_c11_omega;
  }
  return 0;
}

double get_cmb_c20(const int mode){
  switch(mode){
  case 1: return m_cmb_c20_pi0;
  case 2: return m_cmb_c20_eta10;
  case 3: return m_cmb_c20_eta20;
  case 4: return m_cmb_c20_omega;
  }
  return 0;
}

double get_cmb_c1(const int mode){
  switch(mode){
  case 1: return m_cmb_c1_pi0;
  case 2: return m_cmb_c1_eta10;
  case 3: return m_cmb_c1_eta20;
  case 4: return m_cmb_c1_omega;
  }
  return 0;
}

double get_cmb_c2(const int mode){
  switch(mode){
  case 1: return m_cmb_c2_pi0;
  case 2: return m_cmb_c2_eta10;
  case 3: return m_cmb_c2_eta20;
  case 4: return m_cmb_c2_omega;
  }
  return 0;
}

double get_argpar_bb(const int mode){
  switch(mode){
  case 1: return m_cmb_argpar_bb_pi0;
  case 2: return m_cmb_argpar_bb_eta10;
  case 3: return m_cmb_argpar_bb_eta20;
  case 4: return m_cmb_argpar_bb_omega;
  }
  return 0;
}

double get_argpar_qq(const int mode){
  switch(mode){
  case 1: return m_cmb_argpar_qq_pi0;
  case 2: return m_cmb_argpar_qq_eta10;
  case 3: return m_cmb_argpar_qq_eta20;
  case 4: return m_cmb_argpar_qq_omega;
  }
  return 0;
}

double get_mbc0_cmb_bb(const int mode){
  switch(mode){
  case 1: return m_cmb_mbc0_bb_pi0;
  case 2: return m_cmb_mbc0_bb_eta10;
  case 3: return m_cmb_mbc0_bb_eta20;
  case 4: return m_cmb_mbc0_bb_omega;
  }
  return 0;
}

double get_mbcw_cmb_bb(const int mode){
  switch(mode){
  case 1: return m_cmb_mbcw_bb_pi0;
  case 2: return m_cmb_mbcw_bb_eta10;
  case 3: return m_cmb_mbcw_bb_eta20;
  case 4: return m_cmb_mbcw_bb_omega;
  }
  return 0;
}

double get_f_g_cmb_bb(const int mode){
  switch(mode){
  case 1: return m_cmb_f_g_bb_pi0;
  case 2: return m_cmb_f_g_bb_eta10;
  case 3: return m_cmb_f_g_bb_eta20;
  case 4: return m_cmb_f_g_bb_omega;
  }
  return 0;
}

double get_de0r(const int mode){
  switch (mode) {
  case 1: return m_peak_de0r_pi0;
  case 2: return m_peak_de0r_eta10;
  case 3: return m_peak_de0r_eta20;
  case 4: return m_peak_de0r_omega;
  }
  return 0;
}

double get_slopel(const int mode){
  switch (mode) {
  case 1: return m_peak_slopel_pi0;
  case 2: return m_peak_slopel_eta10;
  case 3: return m_peak_slopel_eta20;
  case 4: return m_peak_slopel_omega;
  }
  return 0;
}

double get_sloper(const int mode){
  switch (mode) {
  case 1: return m_peak_sloper_pi0;
  case 2: return m_peak_sloper_eta10;
  case 3: return m_peak_sloper_eta20;
  case 4: return m_peak_sloper_omega;
  }
  return 0;
}

double get_steep(const int mode){
  switch (mode) {
  case 1: return m_peak_steep_pi0;
  case 2: return m_peak_steep_eta10;
  case 3: return m_peak_steep_eta20;
  case 4: return m_peak_steep_omega;
  }
  return 0;
}

double get_p5(const int mode){
  switch (mode) {
  case 1: return m_peak_p5_pi0;
  case 2: return m_peak_p5_eta10;
  case 3: return m_peak_p5_eta20;
  case 4: return m_peak_p5_omega;
  }
  return 0;
}

double get_peak_k_mbc0(const int mode){
  switch(mode){
  case 1: return m_peak_k_mbc0_pi0;
  case 2: return m_peak_k_mbc0_eta10;
  }
  return 0;
}

double get_peak_b_mbc0(const int mode){
  switch(mode){
  case 1: return m_peak_b_mbc0_pi0;
  case 2: return m_peak_b_mbc0_eta10;
  case 3: return m_peak_mbc0_eta20;
  case 4: return m_peak_mbc0_omega;
  }
  return 0;
}

double get_peak_k_s(const int mode){
  switch(mode){
  case 1: return m_peak_k_s_pi0;
  case 2: return m_peak_k_s_eta10;
  }
  return 0;
}

double get_peak_b_s(const int mode){
  switch(mode){
  case 1: return m_peak_b_s_pi0;
  case 2: return m_peak_b_s_eta10;
  case 3: return m_peak_s_eta20;
  case 4: return m_peak_s_omega;
  }
  return 0;
}

double get_tau_bkg(const int mode, const int exp, const int ndf){
  switch(mode){
  case 1: return exp > 30 ? 0.59042 : 0.79461;
  case 2: return 1.20152;
  case 3: return 1.20152;
  case 4: return exp > 30 ? 0.72333 : 1.07632;
  }
  return 0;
}

double get_mu_bkg_delta(const int mode, const int exp){
  switch(mode){
  case 1: return exp > 30 ? -0.0148413 : -0.0277254;
  case 2: return exp > 30 ? -0.0194416 : -0.0111621;
  case 3: return exp > 30 ?  0.0259721 : -0.0356215;
  case 4: return exp > 30 ?  0.0028196 : -0.00802693;
  }
  return 0;
}

double get_mu_bkg(const int mode, const int exp){
  switch(mode){
  case 1: return exp > 30 ? -0.0170737 : -0.018311;
  case 2: return exp > 30 ? -0.0812386 : -0.400879;
  case 3: return exp > 30 ?  0.0024545 :  0.0582829;
  case 4: return exp > 30 ?  0.0543665 :  0.0886096;
  }
  return 0;
}

double get_s_otlr_bkg(const int mode, const int exp){
  switch(mode){
  case 1: return exp < 30 ? 41.2166 : 25.8435;
  case 2: return exp < 30 ? 32.7912 : 31.5617;
  case 3: return exp < 30 ? 22.8887 : 17.979;
  case 4: return exp < 30 ? 30.1879 : 23.3021;
  }
  return 0;
}

double get_f_otlr_bkg(const int mode, const int exp){
  switch(mode){
  case 1: return exp < 30 ? 0.0122323 : 0.0174537;
  case 2: return exp < 30 ? 0.0252212 : 0.0125753;
  case 3: return exp < 30 ? 0.0088723 : 0.0135348;
  case 4: return exp < 30 ? 0.0048445 : 0.00579139;
  }
  return 0;
}

double get_f_bkg_delta(const int mode, const int exp, const int ndf){
  switch(mode){
  case 1:
    if(ndf>0) {return exp > 30 ? 0.654693 : 0.676119;}
    else      {return exp > 30 ? 0.567666 : 0.651968;}
    break;
  case 2:
    if(ndf>0) {return exp > 30 ? 0.896066 : 0.999606;}
    else      {return exp > 30 ? 0.770724 : 0.86151;}
    break;
  case 3:
    if(ndf>0) {return exp > 30 ? 0.828073 : 0.748569;}
    else      {return exp > 30 ? 0.448808 : 0.508574;}
    break;
  case 4:
    if(ndf>0) {return exp > 30 ? 0.514367 : 0.621347;}
    else      {return exp > 30 ? 0.370363 : 0.487717;}
    break;
  }
  return 0;
}

double get_sigma_mn_bkg(const int mode, const int exp, const int ndf){
  switch(mode){
  case 1:
    if(ndf>0) {return exp > 30 ? 1.18994 : 1.21738;}
    else      {return exp > 30 ? 1.12994 : 1.16497;}
    break;
  case 2:
    if(ndf>0) {return exp > 30 ? 1.26948 : 1.22732;}
    else      {return exp > 30 ? 1.14049 : 1.15299;}
    break;
  case 3:
    if(ndf>0) {return exp > 30 ? 1.18743 : 1.10097;}
    else      {return exp > 30 ? 0.12624 : 0.110032;}
    break;
  case 4:
    if(ndf>0) {return exp > 30 ? 1.26818 : 1.27635;}
    else      {return exp > 30 ? 1.08310 : 1.15065;}
    break;
  }
  return 0;
}

double get_sigma_tl_bkg(const int mode, const int exp, const int ndf){
  switch(mode){
  case 1:
    if(ndf>0) {return exp > 30 ? 3.54445 : 4.21694;}
    else      {return exp > 30 ? 3.94822 : 5.52881;}
    break;
  case 2:
    if(ndf>0) {return exp > 30 ? 3.36794 : 2.58338;}
    else      {return exp > 30 ? 5.47538 : 3.55447;}
    break;
  case 3:
    if(ndf>0) {return exp > 30 ? 2.30686 : 2.56057;}
    else      {return exp > 30 ? 8.38917 : 9.95;}
    break;
  case 4:
    if(ndf>0) {return exp > 30 ? 4.47635 : 4.79342;}
    else      {return exp > 30 ? 6.81437 : 15;}
    break;
  }
  return 0;
}

double get_f_tl_bkg(const int mode, const int exp, const int ndf){
  switch(mode){
  case 1:
    if(ndf>0) {return exp > 30 ? 0.149367 : 0.105607;}
    else      {return exp > 30 ? 0.126601 : 0.109805;}
    break;
  case 2:
    if(ndf>0) {return exp > 30 ? 0.141492 : 0.266935;}
    else      {return exp > 30 ? 0.089761 : 0.110071;}
    break;
  case 3:
    if(ndf>0) {return exp > 30 ? 0.372017 : 0.318554;}
    else      {return exp > 30 ? 0 : 0.070058;}
    break;
  case 4:
    if(ndf>0) {return exp > 30 ? 0.090121 : 0.0590122;}
    else      {return exp > 30 ? 0.075740 : 0.0471682;}
    break;
  }
  return 0;
}

#endif // CUTS_H
