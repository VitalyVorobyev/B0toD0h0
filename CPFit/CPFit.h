#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnStrategy.h"

#include "TPaveStats.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <utility>
#include <fstream>
#include <cstdlib>

#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>

#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TStyle.h"

//#include "Math/SMatrix.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"

#include "RkRdetRnpPdf.h"
#include "RbkgPdf.h"

#include "../DalitzModelStudy/dalitzmodel.h"
#include "../MyParams/myparams.h"

using namespace std;
//typedef ROOT::Math::SMatrix<double,8> SymD8Matrix;
//typedef ROOT::Math::SMatrix<double,16> SymD16Matrix;

const int toy_composition[4] = {414, 124, 40, 246};
int NToySamples = 0;
TTree* tree_ww_pi0;
TTree* tree_ww_etagg;
TTree* tree_ww_etappp;
TTree* tree_ww_omega;
TChain *m_tree_pi0;
TChain *m_tree_etagg;
TChain *m_tree_etappp;
TChain *m_tree_omega;
vector<int> index_pi0;
vector<int> index_etagg;
vector<int> index_etappp;
vector<int> index_omega;

void GetShuffledVector(const int size,vector<int> &vec){
  vec.clear();
  for(int i=0; i<size; i++) vec.push_back(i);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  shuffle(vec.begin(),vec.end(),default_random_engine(seed));
  return;
}

double Carr[8], Sarr[8], Karr[8], Kbarr[8];

inline double C(const int bin){ return Carr[abs(bin)-1];}
inline double S(const int bin){ return bin>0 ? Sarr[abs(bin)-1] : -Sarr[abs(bin)-1];}
inline double K(const int bin){ return bin>0 ? Karr[abs(bin)-1] : Kbarr[abs(bin)-1];}

const int nevbins[7]             = {48443,33788,42434,28549,24569,24584,48516};
const double w_data_svd1_posi[7] = {0.,7.235697e-03,7.129388e-03,7.417778e-03,6.885875e-03,6.761047e-03,4.336734e-03};
const double w_data_svd1_nega[7] = {0.,6.001569e-03,6.430566e-03,7.693083e-03,6.416449e-03,8.807757e-03,4.587614e-03};
const double w_data_svd2_posi[7] = {0.,4.152612e-03,3.243236e-03,3.721417e-03,3.315138e-03,3.180302e-03,2.175087e-03};
const double w_data_svd2_nega[7] = {0.,3.577812e-03,2.803811e-03,3.486607e-03,4.241595e-03,3.696399e-03,3.077622e-03};

double err_wt_svd1 = 0.5;
double err_wt_svd2 = 0.5;

vector< vector<double> > NsigVec;
vector< vector<double> > NsigErrs;
vector< vector<double> > NbsigVec;
vector< vector<double> > NbsigErrs;
vector<double> fbbVec;
vector<double> fbbErrs;
vector<double> fprtVec;
vector<double> fprtErrs;
double aver_wtag = 0.5;

void init_wt_err(void){
  int NTot = 0;
  err_wt_svd1 = 0;
  err_wt_svd2 = 0;
  for(int i=0; i<7; i++){
    NTot += nevbins[i];
    err_wt_svd1 += nevbins[i]*w_data_svd1_posi[i];
    err_wt_svd2 += nevbins[i]*w_data_svd2_posi[i];
  }
  err_wt_svd1 /= NTot;
  err_wt_svd2 /= NTot;

  cout << "Wrong tagging errors are initialized:" << endl;
  cout << "  SVD1: " << err_wt_svd1 << endl;
  cout << "  SVD2: " << err_wt_svd2 << endl;
  return;
}

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

typedef TMatrixTSym<double> TMatrixDSym;
typedef TMatrixT<double> TMatrixD;

const double Carr_CLEO[8] = { 0.710, 0.365,0.008,-0.757,-0.884,-0.462, 0.106, 0.481};
const double Sarr_CLEO[8] = {-0.013,-0.179,0.938, 0.386,-0.162,-0.616,-1.063,-0.147};
const double Karr_CLEO[8] = { 0.159, 0.123,0.102, 0.079, 0.089, 0.030, 0.093, 0.067};
const double Kbarr_CLEO[8]= { 0.109, 0.026,0.012, 0.015, 0.040, 0.013, 0.027, 0.016};

// Numerical integral //
const double Carr_model[8] = { 0.675798,0.431828,-0.036724,-0.642989,-0.935997,-0.615002, 0.003410, 0.571212};
const double Sarr_model[8] = {-0.005376,0.413289, 0.725041, 0.514323,-0.017438,-0.668626,-0.814657,-0.416431};
const double Karr_model[8] = { 0.168795,0.118502, 0.096164, 0.074168, 0.090667, 0.030814, 0.105442, 0.077863};
const double Kbarr_model[8]= { 0.088009,0.028629, 0.012360, 0.014960, 0.042507, 0.010055, 0.023245, 0.017819};

// Corrected with gen fit //
//const double Carr_model[8]   = { 0.681017,0.447322,-0.033204,-0.620657,-0.934699,-0.642683,-0.007068, 0.550722};
//const double Sarr_model[8]   = {-0.008337,0.413185, 0.715935, 0.527492,-0.010901,-0.645881,-0.826176,-0.434764};
//const double Karr_model[8] = {};
//const double Kbarr_model[8]= {};

const double Karr_omega[8] = {0.170277,0.118646,0.0959634,0.0765277,0.0926909,0.0308337,0.100305,0.0719257};
const double Kbarr_omega[8]= {0.0909526,0.0297611,0.0115792,0.015016,0.0440408,0.0107163,0.0228221,0.0179427};

const double Karr_pi0[8] = {0.168992,0.119478,0.0960086,0.0753697,0.0935319,0.0301996,0.0967017,0.0718133};
const double Kbarr_pi0[8]= {0.0942696,0.0295868,0.0132924,0.0149604,0.0425527,0.0112847,0.0229337,0.0190245};

const double Karr_err[8] = {0.003238,0.002848,0.002648,0.002231,0.00250,0.001425,0.002588,0.002251};
const double Kbarr_err[8]= {0.002424,0.001402,0.001001,0.001161,0.00178,0.000963,0.001438,0.001054};

const double CS_stat_err_arr[256] = {0, -3, -2,  5, 10,  0,  1, -3,  0,  0,  0,  0,  0,  0,  0,  0,
                                     0,  0,  2, -3,  8,  1,  1,  0,  0, -3,  0,  0, -1,  0,  0,  0,
                                     0,  0,  0, -2,  1, -1, -1,  2,  1,  0,  3, -1,  0,  2, -2, -5,
                                     0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  5,  0,  0,  0,  0,
                                     0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  3,  0,  0,  0,
                                     0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  0, -1,  0,  1,  0,  2,
                                     0,  0,  0,  0,  0,  0,  0,  0,  2,  1,  1, -5,  0,  1, -2,  1,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  9,  3,-18,  8,  5, 23,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3, -7,  9,-10, -5,-13,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 24,  1,  3, 57,  9,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  7, 33, 12,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -5,  0, -6,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  8,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

const double CS_syst_err_arr[256] = {0, 87, 93, 72, 85, 87, 85, 92, 31, 21, 29,  0, 20,  0, -4, -6,
                                     0,  0, 89, 77, 87, 89, 88, 89, 20,  0, 20,-10,  1, 11,-13,  8,
                                     0,  0,  0, 73, 88, 90, 89, 94, 29, 17, 29, -3, 13,  6, -6, -1,
                                     0,  0,  0,  0, 84, 76, 78, 70,  2,-20,  6,-17,-13, 20,-22, 18,
                                     0,  0,  0,  0,  0, 86, 86, 84, 10, -4, 14,-12,  1, 12,-18,  7,
                                     0,  0,  0,  0,  0,  0, 88, 88, 15,  0, 15,-10,  2, 14,-15,  4,
                                     0,  0,  0,  0,  0,  0,  0, 84, 15, -8, 11,-15, -7, 20,-18, 13,
                                     0,  0,  0,  0,  0,  0,  0,  0, 30, 23, 32,  0, 22,  0, -2, -5,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0, 43, 42, 20, 19,  0, 32, -6,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 44, 37, 65,-37, 44,-52,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 21, 43,  0, 51,  0,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 31, -9, 38,-13,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-41, 33,-42,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-19, 44,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-11,
                                     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

TMatrixDSym CS_covMtx(16);
TMatrixDSym CS_syst_covMtx(16);
TMatrixDSym CS_inv_covMtx(16);

const double C_stat_err[8] = {3.4,7.1,10.5,10.0,5.6,9.9,8.0,8.0};
const double S_stat_err[8] = {9.7,16.6,17.4,18.8,13.0,20.8,12.0,17.7};

const double C_syst_err[8] = {3.8,7.8,10.0,8.2,5.4,6.5,8.7,7.0};
const double S_syst_err[8] = {3.1,4.8,6.6,5.2,4.1,6.7,4.7,10.7};

double K0[8],Kb0[8],C0[8],S0[8];

void calc_covariance(void){
  cout << "Calculating Covariance and Concentration matrices..." << endl;
  CS_covMtx.SetMatrixArray(CS_stat_err_arr,"F");
  CS_syst_covMtx.SetMatrixArray(CS_syst_err_arr,"F");
  double stat_err1, stat_err2;
  double syst_err1, syst_err2;
  for(int i=0; i<16; i++){
    for(int j=0; j<=i; j++){
      if(j == i){
        if(i<8){// C
          CS_covMtx[i][i] = 0.0001*C_stat_err[i]*C_stat_err[i];
          CS_syst_covMtx[i][i] = 0.0001*C_syst_err[i]*C_syst_err[i];
        } else{// S
          CS_covMtx[i][i] = 0.0001*S_stat_err[i-8]*S_stat_err[i-8];
          CS_syst_covMtx[i][i] = 0.0001*S_syst_err[i-8]*S_syst_err[i-8];
        }
      } else{// j != i
        if(i<8){// C x C
          stat_err1 = C_stat_err[i]; stat_err2 = C_stat_err[j];
          syst_err1 = C_syst_err[i]; syst_err2 = C_syst_err[j];
        } else if(j<8){// S x C
          stat_err1 = S_stat_err[i-8]; stat_err2 = C_stat_err[j];
          syst_err1 = S_syst_err[i-8]; syst_err2 = C_syst_err[j];
        } else{// S x S
          stat_err1 = S_stat_err[i-8]; stat_err2 = S_stat_err[j-8];
          syst_err1 = S_syst_err[i-8]; syst_err2 = S_syst_err[j-8];
        }
        CS_covMtx[i][j] *= 0.000001*stat_err1*stat_err2;
        CS_syst_covMtx[i][j] *= 0.000001*syst_err1*syst_err2;
        CS_covMtx[j][i] = CS_covMtx[i][j];
        CS_syst_covMtx[j][i] = CS_syst_covMtx[i][j];
      }
    }
  }
  ofstream stat_cov_file("CS_stat_cov_matrix.txt",ofstream::out);
  ofstream syst_cov_file("CS_syst_cov_matrix.txt",ofstream::out);
  ofstream full_con_file("CS_full_con_matrix.txt",ofstream::out);
  stat_cov_file << "Matrix of statistical covariance for S and C" << endl;
  syst_cov_file << "Matrix of systematic covariance for S and C" << endl;
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      stat_cov_file << CS_covMtx[i][j] << " ";
      syst_cov_file << CS_syst_covMtx[i][j] << " ";
    }
    stat_cov_file << endl;
    syst_cov_file << endl;
  }
  CS_covMtx += CS_syst_covMtx;
  CS_inv_covMtx = CS_covMtx.Invert();

  full_con_file << "Concentration matrix for S and C" << endl;
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      full_con_file << CS_covMtx[i][j] << " ";
    }
    full_con_file << endl;
  }


  stat_cov_file.close();
  syst_cov_file.close();
  full_con_file.close();

  cout << "Done." << endl;
  return;
}

#define SIGNALONLY 1
#define BACKGROUND 2
#define ABFIT     -1

bool nega_pdf_flag = false;

RkRdetRnpPdf* m_pdf_svd2,      *m_pdf_svd1;
RbkgPdf*      m_pdf_back_svd2, *m_pdf_back_svd1;

RbkgPdf*      m_pdf_back_svd2_bb_gg, *m_pdf_back_svd1_bb_gg;
RbkgPdf*      m_pdf_back_svd2_cont_gg, *m_pdf_back_svd1_cont_gg;
RbkgPdf*      m_pdf_back_svd2_bb_ppp, *m_pdf_back_svd1_bb_ppp;
RbkgPdf*      m_pdf_back_svd2_cont_ppp, *m_pdf_back_svd1_cont_ppp;
TChain* m_tree;
TChain* m_test_tree;
TTree*  m_good_tree;
TTree*  m_good_sb_tree;
TTree*  m_sig_tree,*m_bkg_tree;
TChain* m_toy_bkg_tree;
TTree*  m_ww_tree;
TTree*  m_ww_sb_tree;

Int_t m_exp,m_bin,m_flv;//,m_good_icpv;
Int_t m_bin_mc, m_flv_mc;
Int_t m_b0f;//,m_d0f,m_h0f;
Double_t m_dt,m_t_sig,m_t_asc;
Double_t m_t_sig_mc, m_t_asc_mc;
Double_t m_A,m_B;
Double_t m_costhBcms;
Double_t m_f_cont_in_comb,m_f_cont_in_comb_bin_mc,m_f_cont_in_comb_flv_mc,m_f_cont_in_comb_mc;
Double_t m_f_cont,m_f_cont_bin_mc,m_f_cont_flv_mc,m_f_cont_mc;
Double_t m_f_bkg,m_f_bkg_bin_mc,m_f_bkg_flv_mc,m_f_bkg_mc;
double m_f_bkg_err,m_f_cont_in_comb_err;
Int_t m_sigarea;
//Double_t m_chi2_vtx_d0;
int m_svd = 0;
int m_mode = 4;
int mm_mode, mm_h0mode;
int m_fitbin = 0;
int m_fitflv = 0;
int m_type_flag = 0;// 1 -> cont, 2 -> BB
//double m_f_cont = 0.7757;//0.4293;

double m_S_main_mlt_cont_svd1_gg, m_S_main_sgl_cont_svd1_gg;
double m_S_main_mlt_cont_svd2_gg, m_S_main_sgl_cont_svd2_gg;
double m_S_main_mlt_bb_svd1_gg,   m_S_main_sgl_bb_svd1_gg;
double m_S_main_mlt_bb_svd2_gg,   m_S_main_sgl_bb_svd2_gg;

double m_S_main_mlt_cont_svd1_ppp, m_S_main_sgl_cont_svd1_ppp;
double m_S_main_mlt_cont_svd2_ppp, m_S_main_sgl_cont_svd2_ppp;
double m_S_main_mlt_bb_svd1_ppp,   m_S_main_sgl_bb_svd1_ppp;
double m_S_main_mlt_bb_svd2_ppp,   m_S_main_sgl_bb_svd2_ppp;

// * Event-by-event parameters
int m_ntrk_rec;
int m_ntrk_asc;
int m_ndf_rec;
int m_ndf_asc;
double m_sz_rec;
double m_sz_asc;
double m_chisq_rec;
double m_chisq_asc;
double m_mp;
double m_mm;
double m_mp_mc;
double m_mm_mc;

double m_Nsig;
double m_Nsig_err;
int m_Ntot;
double m_psig;
double m_pcnt;
double m_pprt;
double m_pcmb;
double m_fbb;
double m_fbb_err;
double m_fprt;
double m_fprt_err;
double m_wtag;
// *

// Nuisance parameters
double m_wr_tag_offset_svd1 = 0;
double m_wr_tag_offset_svd2 = 0;
vector<double> m_f_bb_offset;
vector<double> m_f_prt_offset;
vector< vector<double> > m_Nsig_offset;
vector< vector<double> > m_Nbsig_offset;
vector<int> modes_set;

using namespace std;
using namespace ROOT;
using namespace Minuit2;

int Mode(const int mode, const int h0mode){
  if(mode == 1 && h0mode == 10)  return 1;
  if(mode == 2 && h0mode == 10)  return 2;
  if(mode == 2 && h0mode == 20)  return 3;
  if(mode == 3 && h0mode == 20)  return 4;
  if(mode == 5 && h0mode == 10)  return 5;
  if(mode == 10 && h0mode == 10) return 10;
  if(mode == 20 && h0mode == 10) return 20;
  return 0;
}

string get_label(const int mode){
  string label;
  switch(mode) {
  case 0:
    label = string(", Data");
    break;
  case 1:
    label = string(", #pi^{0}");
    break;
  case 10:
    label = string(", D^{*}#pi^{0}");
    break;
  case 2:
    label = string(", #eta#rightarrow#gamma#gamma");
    break;
  case 20:
    label = string(", D^{*}#eta");
    break;
  case 3:
    label = string(", #eta#rightarrow#pi^{+}#pi^{+}#pi^{0}");
    break;
  case 4:
    label = string(", #omega");
    break;
  case 5:
    label = string(", #eta`");
//    label = string(", #omega (#rho^{0}(770)K_{S}^{0})");
    break;
  case 6:
    label = string(", #rho");
    break;
  default:
    label = string("");
    break;
  }
  return label;
}

const double sol = 2.99792458;
bool draw_plots = true;
const double cm2ps  = 78.48566945838871754705;
const double m_btau = 1.534;//1.520;//0.460/sol*cm2ps;//0.507;// ps
double mm_btau = m_btau;
const double m_btau_err = 0.05;
const double m_dm   = 0.510;//0.452;//0.510/sol;//0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
const double m_dm_err = 0.003;
double mm_dm = m_dm;
double beta = 67.;
double m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
double m_cos2beta = TMath::Cos(2.*beta/180.*TMath::Pi());
double mm_sin2beta = m_sin2beta;
double mm_cos2beta = m_cos2beta;
int xi;
const int NDots = 250;
const int NBins = NDots;
double dtmax = 70;
double dtmin = -dtmax;

int m_NSig    = 600;
double m_fSig = 0.6;
int m_NBkg    = m_NSig/m_fSig-m_NSig;
int m_nsets   = 1000000;
double m_scale1 = 1.;
double m_scale2 = 1.;
double m_scale1_err = 100.;
double m_scale2_err = 100.;

//const double dzmax =  dtmax/78.4857;
//const double dzmin = -dzmax;

bool pipi_fit      = false;
bool d0_fit        = false && !pipi_fit;
const bool draw_b_bbar = false;
const bool make_k_fit  = false;
bool make_Rrec_fit = false;
bool make_Rasc_fit = false;
bool make_Rnp_fit  = false;
bool sideband_fit  = false;
bool add_otlr      = false;
bool no_interf     = false || make_Rnp_fit;
bool no_np         = false && !make_Rnp_fit;
bool fix_btau      = true  || make_k_fit || make_Rnp_fit;
bool fix_dm        = true  || make_k_fit || make_Rnp_fit || no_interf;
bool fix_sin2beta  = false || make_k_fit || make_Rnp_fit || no_interf;
bool fix_cos2beta  = false || make_k_fit || make_Rnp_fit || no_interf;
bool fix_A         = false;
bool fix_B         = false;
bool m_cleo        = false;
bool sgl_asc       = false;
bool mlt_asc       = false;
bool include_s0    = false;
bool include_dt0   = false;
bool no_bkg        = false;
bool sigmc         = false;
bool perftag       = false;
bool perfbin       = false;
bool full_ds_fit   = false;
bool fix_pdf       = false;
bool ABfit         = false;
bool FullFit       = true;
bool make_bins_scan= false;
bool toybkg        = false;
bool m_ebeb        = false;
bool m_data        = false;
bool m_norm_test_flag = false;
int m_Neve = 0;
bool m_ww  = false;
bool m_gg  = false;
bool m_ppp = false;
bool m_calc_K = false;
bool m_calc_CS = false;
bool m_nuisance = false;
bool m_genfit = false;
bool m_line_test = false;
bool m_toyfit = false;

double m_tag;
double m_s0;
int m_NSigTot  = 0;
int m_NBkgTot  = 0;
int m_NSBTot   = 0;
int m_NGoodTot = 0;
double m_fbkg_tot = 0;
int m_sig_map[2][16];
int m_bkg_map[2][16];
double m_fbkg_map[2][16];
vector<double> m_fb_vec;

double m_f_ol_sgl_svd1 = 0;
double m_f_ol_mlt_svd1 = 0;
double m_s_ol_svd1     = 0;
double m_f_ol_sgl_svd2 = 0;
double m_f_ol_mlt_svd2 = 0;
double m_s_ol_svd2     = 0;
int m_toysize = 0;
int m_ns = 1;
int m_cs = 0;

double sum_sigma(const double& s1, const double& s2){
  return sqrt(s1*s1+s2*s2);
}

typedef struct ICPVEVT{
  int exp;
  int flv;
  double tag;
  int bin;
  double costhBcms;

  int nrtk_sig;
  int ndf_sig;
  double chisq_sig;
  double t_sig;
  double sz_sig;

  int nrtk_asc;
  int m_ndf_asc;
  double m_chisq_asc;
  double t_asc;
  double m_sz_asc;
} ICPVEvt;

vector<ICPVEvt> m_sig_vec;
vector<ICPVEvt> m_bkg_vec;

void PrintEvent(void){
  cout << " " << m_ntrk_rec  << " " << m_sz_rec;
  cout << " " << m_chisq_rec << " " << m_ndf_rec;
  cout << " " << m_ntrk_asc  << " " << m_sz_asc;
  cout << " " << m_chisq_asc << " " << m_ndf_asc << endl;
}

void calc_fbkg_fcnt(void){
  if(no_bkg){
    m_f_bkg = 0;
    m_f_cont_in_comb = 0.5;
    return;
  }
  MyParams cuts;
  int mode_ind;
  const double mode = Mode(mm_mode,mm_h0mode);
  for(int j=0; j<modes_set.size(); j++){
    if(modes_set[j] == mode){ mode_ind = j; break;}
  }
  const int bin_ind = cuts.bin_ind(m_bin);
  const double Nsig_offset = m_flv>0 ? m_Nsig_offset[mode_ind][bin_ind] : m_Nbsig_offset[mode_ind][bin_ind];// Nsig offset
  const double wr_tag_offset = m_exp > 30 ? m_wr_tag_offset_svd2 : m_wr_tag_offset_svd1;
  const double Nsig     = 0.5*m_Nsig*cuts.N(m_bin,m_flv,m_wtag+wr_tag_offset)+Nsig_offset;// wr tag offset
  const double Nbkg     = m_Ntot - Nsig > 0 ? m_Ntot - Nsig : 0;
  const double fprt_bb  = m_fprt+m_f_prt_offset[mode_ind];// fprt offset
  const double fbb      = m_fbb+m_f_bb_offset[mode_ind];// fbb offset
  const double fprt     = fbb*fprt_bb;
  const double Nprt     = fprt/(1.+fprt)*Nbkg;
  const double Ncnt     = (1.-fbb)/(1+fprt)*Nbkg;
  const double Ncmb     = 1./(1.+fprt)*Nbkg;

  const double Vsig = m_psig*Nsig;
  const double Vcnt = m_pcnt*Ncnt;
  const double Vprt = m_pprt*Nprt;
  const double Vcmb = m_pcmb*Ncmb;

  m_f_bkg = (Vcmb+Vprt+Vsig)>0 ? (Vcmb+Vprt)/(Vcmb+Vprt+Vsig) : 0;
  m_f_cont_in_comb = Vcmb>0 ? Vcnt/Vcmb : 0;
  return;
}


inline double A(const int flv,const int bin){ return flv*(K(bin)-K(-bin))/(K(bin)+K(-bin));}
inline double B(const int flv,const int bin){ return 2.*flv*xi*sqrt(K(bin)*K(-bin))/(K(bin)+K(-bin))*(C(bin)*m_sin2beta+S(bin)*m_cos2beta);}

// CLEO measurements
//const double Carr_CLEO[8] = { 0.365, 0.710, 0.481,0.008,-0.757,-0.884,-0.462, 0.106};
//const double Sarr_CLEO[8] = {-0.179,-0.013,-0.147,0.938, 0.386,-0.162,-0.616,-1.063};
//const double Carr_CLEO[8] = {  0.481, 0.106,-0.462,-0.884,-0.757, 0.008, 0.365, 0.710};
//const double Sarr_CLEO[8] = { -0.147,-1.063,-0.616,-0.162, 0.386, 0.938,-0.179,-0.013};
//const double Karr_CLEO[8] = {  0.067, 0.093, 0.030, 0.089, 0.079, 0.102, 0.123, 0.159};
//const double Kbarr_CLEO[8]= {  0.016, 0.027, 0.013, 0.040, 0.015, 0.012, 0.026, 0.109};

//Model integrals
//const double Carr_model[8] = {0.594511, -0.329868,-0.61714, -0.758994,-0.368921,  0.133975, 0.484073, 0.742303};
//const double Sarr_model[8] = {-0.357055,-0.680259,-0.277483,-0.259506, 0.440863,  0.623397, 0.358263,-0.0558073};
//const double Karr_model[8] = {0.0665352,0.0835056, 0.0360653,0.0961057,0.0722859, 0.0980402,0.123951, 0.198887};
//const double Kbarr_model[8]= {0.0304473,0.00571654,0.0022683,0.0264409,0.00969127,0.0105694,0.0311109,0.10838};

//const double Carr_model[8] = { 0.554166,-0.00454239,-0.608498,-0.934161,  -0.565509,0.0800234,0.480612, 0.68653};
//const double Sarr_model[8] = {-0.446559,-0.825851,  -0.678567, 0.000976537,0.572031,0.734359, 0.372541,-0.0279944};
//const double Karr_model[8] = {0.0721528,0.0963476,0.0314874,0.0973232,0.079597,0.0991814,0.120779,0.158244};
//const double Kbarr_model[8]= {0.0188839,0.0241389,0.0105323,0.045323,0.0152623,0.0124246,0.0286007,0.0897215};

//const double Carr_model[8] = {0.586785,0.00996666,-0.615205,-0.933178,-0.576292,0.0177753,0.448619,0.681499};
//const double Sarr_model[8] = {-0.411777,-0.819878,-0.629422,0.0157803,0.566233,0.736251,0.419118,0.00536454};
//const double Karr_model[8] = {0.0823076,0.101649,0.0394897,0.0949997,0.0792316,0.0893498,0.114192,0.160569};
//const double Kbarr_model[8]= {0.0186703,0.0234387,0.0114402,0.0438988,0.0139401,0.0117233,0.0238303,0.0912706};

//const double Carr_model[8] = {0.681499,0.448619,0.017775,-0.576292,-0.933178,-0.615205, 0.009967, 0.586785};
//const double Sarr_model[8] = {0.005365,0.419118,0.736251, 0.566233, 0.015780,-0.629422,-0.819878,-0.411777};
//const double Karr_model[8] = {0.160569,0.114192,0.089350, 0.079232, 0.095000, 0.039490, 0.101649, 0.082308};
//const double Kbarr_model[8]= {0.091271,0.023830,0.011723, 0.013940, 0.043899, 0.011440, 0.023439, 0.018670};

//pi0 K parameters:
// bin 1, K = 16.8992 +- 0.114457, Kb = 9.42696 +- 0.0999748
// bin 2, K = 11.9478 +- 0.0952888, Kb = 2.95868 +- 0.0690292
// bin 3, K = 9.60086 +- 0.0851777, Kb = 1.32924 +- 0.0556704
// bin 4, K = 7.53697 +- 0.0767571, Kb = 1.49604 +- 0.052837
// bin 5, K = 9.35319 +- 0.0874837, Kb = 4.25527 +- 0.0718956
// bin 6, K = 3.01996 +- 0.0506911, Kb = 1.12847 +- 0.0393301
// bin 7, K = 9.67017 +- 0.0865278, Kb = 2.29337 +- 0.0618354
// bin 8, K = 7.18133 +- 0.0756784, Kb = 1.90245 +- 0.0549534

//Omega K parameters:
// bin 1, K = 17.0277 +- 0.140264, Kb = 9.09526 +- 0.121309
// bin 2, K = 11.8646 +- 0.116381, Kb = 2.97611 +- 0.0845015
// bin 3, K = 9.59634 +- 0.104054, Kb = 1.15792 +- 0.0667511
// bin 4, K = 7.65277 +- 0.0946452, Kb = 1.5016 +- 0.0650375
// bin 5, K = 9.26909 +- 0.10695, Kb = 4.40408 +- 0.0887972
// bin 6, K = 3.08337 +- 0.062504, Kb = 1.07163 +- 0.0477013
// bin 7, K = 10.0305 +- 0.107617, Kb = 2.28221 +- 0.076375
// bin 8, K = 7.19257 +- 0.0925689, Kb = 1.79427 +- 0.0664242


// rho Ks
const double Carr_rhoKs[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
const double Sarr_rhoKs[8] = {0,0,0,0,0,0,0,0};
const double Karr_rhoKs[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_rhoKs[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

// no interfer
const double Carr_noint[8] = {0,0,0,0,0,0,0,0};
const double Sarr_noint[8] = {0,0,0,0,0,0,0,0};
const double Karr_noint[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_noint[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

int xil(const int mode){
  if(!m_data) return 1;
  switch (mode) {
  case 1:  return  1; // pi0
//  case 10: return -1; // D*0 pi0
  case 2:  return  1; // eta->gg
//  case 20: return -1; // D*0 eta
  case 3:  return 1; // eta->ppp
//  case 4:  return -1; // omega
//  case 5:  return -1; // eta'
  }
  return -1;
}

int Xil(const int mode, const int h0mode){
  if(!m_data) return 1;
//  return 1;
  switch (mode) {
  case 1:  return  1; // pi0
//  case 10: return -1; // D*0 pi0
  case 2:  return  1; // eta
//  case 20: return -1; // D*0 eta
//  case 3:  return -1; // omega
//  case 5:  return -1; // eta'
  }
  return -1;
}

int Beta(const int mode){
  return 23;
  switch (mode) {
//  case 1:  return 67.; // pi0
//  case 10: return 67; // D*0 pi0
  case 2:  return 23; // eta->gg
//  case 20: return 67; // D*0 eta
  case 3:  return 23; // eta->ppp
//  case 4:  return 67; // omega
  case 5:  return 23; // eta'
  }
  return 67;
}

int Beta(const int mode, const int h0mode){
  return 23;
  switch (mode) {
//  case 1:  return 67.; // pi0
//  case 10: return 67; // D*0 pi0
  case 2:  return 23; // eta->gg
//  case 20: return 67; // D*0 eta
//  case 4:  return 67; // omega
  case 5:  return 23; // eta'
  }
  return 67;
}

int ppp_flag(const int mode){
  switch (mode) {
  case 1:  return 2; // pi0
  case 10: return 2; // D*0 pi0
  case 2:  return 2; // eta->gg
  case 20: return 2; // D*0 eta
  case 3:  return 1; // eta->ppp
  case 4:  return 1; // omega
  case 5:  return 1; // eta'
  }
  return 1;
}

int Mode(const int mode){
  switch (mode) {
  case 1:  return 1; // pi0
  case 10: return 10;// D*0 pi0
  case 20: return 20;// D*0 eta
  case 4:  return 3; // omega
  case 5:  return 5; // eta'
  case 77: return 77; // eta'
  default: return 2; // eta
  }
}

int h0Mode(const int mode){
  switch (mode) {
  case 1:  return 10; // pi0
  case 10: return 10; // D*0 pi0
  case 2:  return 10; // eta->gg
  case 20: return 10; // D*0 eta
  case 3:  return 20; // eta->ppp
  case 4:  return 20; // omega
  case 5:  return 10; // eta'
  case 77: return 77; // eta'
  }
}

void DrawSigPredictions(void){
  MyParams cuts;
  const double bins[16]     = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
  const double bins_err[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double true_sig[2][16];
  stringstream out;
  for(int i=0; i<NsigVec.size(); i++){
    if(!m_data){
      for(int k=0; k<2; k++){
        for(int j=0; j<16; j++){
          out.str("");
//          out << "flv == " << cuts.flv(k);
          out << "flv_mc == " << cuts.flv(k);
//          out << " && bin == " << cuts.bin(j);
          out << " && bin_mc == " << cuts.bin(j);
          out << " && mode == " << Mode(modes_set[i]) << " && h0mode == " << h0Mode(modes_set[i]);
          out << " && (b0f == 1 || b0f == 5 || b0f == 10)";
          true_sig[k][j] = m_ww_tree->Draw("de",out.str().c_str());
//          cout << k << " " << j << " " << true_sig[k][j] << endl;
        }
      }
    }
    TCanvas* c1 = new TCanvas("c1","c1",800,400);
    c1->Draw();
    TPad* pad1 = new TPad("pad1","pad1",0.01,0.0,0.49,0.99);
    pad1->Draw();
    pad1->SetGrid();
    TPad* pad2 = new TPad("pad2","pad2",0.51,0.0,0.99,0.99);
    pad2->Draw();
    pad2->SetGrid();
    pad1->cd();
    TGraphErrors* gr = new TGraphErrors(16,bins,&NsigVec[i][0],bins_err,&NsigErrs[i][0]);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerColor(kBlue);
    if(m_data){
      gr->Draw("ap");
    } else{
      TMultiGraph* mg = new TMultiGraph("mg","mg");
      TGraph* gr0 = new TGraph(16,bins,true_sig[0]);
      gr0->SetMarkerStyle(21);
      gr0->SetMarkerSize(1.2);
      gr0->SetMarkerColor(kRed);

      mg->Add(gr);
      mg->Add(gr0);
      mg->Draw("ap");
    }
    pad2->cd();
    TGraphErrors* grb = new TGraphErrors(16,bins,&NbsigVec[i][0],bins_err,&NbsigErrs[i][0]);
    grb->SetMarkerStyle(20);
    grb->SetMarkerSize(1.2);
    grb->SetMarkerColor(kBlue);
    if(m_data){
      grb->Draw("ap");
    } else{
      TMultiGraph* mgb = new TMultiGraph("mgb","mgb");
      TGraph* grb0 = new TGraph(16,bins,true_sig[1]);
      grb0->SetMarkerStyle(21);
      grb0->SetMarkerSize(1.2);
      grb0->SetMarkerColor(kRed);

      mgb->Add(grb);
      mgb->Add(grb0);
      mgb->Draw("ap");
    }
    c1->Print("test.eps");
    system("evince test.eps &");
  }
}

void GetEvent(const int i,const int SetNum=0,const int type=0){
// type 0 -> Signal
// type 1 -> Background
// type 2 -> Sideband
  if(m_ebeb){
    if(!SetNum) m_good_tree->GetEvent(i);
    else{
      const int j = (m_NBkg+m_NSig)*SetNum + i;
      m_good_tree->GetEvent(j);
    }
  } else if(type == 2){
    m_good_sb_tree->GetEvent(i);
//    cout << m_sz_rec << " " << m_sz_asc << endl;
  }
  else if(type == 1){
    m_bkg_tree->GetEvent(i);
//    cout << m_sz_rec << " " << m_sz_asc << " " << m_t_sig << " " << m_t_asc << " " << m_ndf_asc << endl;
  }
  else if(!m_mode){ m_good_tree->GetEvent(i);}
  else if(type == 0){
    if(SetNum){
      const int cycle = m_NSig*SetNum/m_NSigTot;
      int j = (SetNum-1)*m_NSig + i - cycle*m_NSigTot;
      m_sig_tree->GetEvent(j);
    } else{
      if(i<m_NSigTot) m_sig_tree->GetEvent(i);
      else{
        cout << "GetEvent: i = " << i << " > NSigTot = " << m_NSigTot << endl;
        return;
      }
    }
  } else{
    if(SetNum){
      const int cycle = m_NBkg*SetNum/m_NBkgTot;
      int j = (SetNum-1)*m_NBkg + i - cycle*m_NBkgTot;
        m_bkg_tree->GetEvent(j);
    } else{
      if((i<m_NBkgTot)){
        m_bkg_tree->GetEvent(i);
      } else{
        cout << "GetEvent: i = " << i << " > NBkgTot = " << m_NBkgTot << endl;
        return;
      }
    }
  }
  m_dt = (m_t_sig - m_t_asc)*0.1*cm2ps;
  if(pipi_fit){
    m_ntrk_rec = 2;
    m_ndf_rec = 2;
  }
  if(d0_fit){
    m_ntrk_rec = 1;
    m_ndf_rec = 0;
  }
  m_sz_rec /= 10.;
  m_sz_asc /= 10.;

  if(!perftag){ m_flv = m_tag > 0 ? -1 : 1;}
  if(no_bkg) m_f_bkg = 0;
  return;
}

void GetEventWW(const int i,const int mode = 0){
  if(m_toyfit){
    switch (mode) {
    case 0:
      tree_ww_pi0->GetEvent(index_pi0[i]);
      break;
    case 1:
      tree_ww_etagg->GetEvent(index_etagg[i]);
      break;
    case 2:
      tree_ww_etappp->GetEvent(index_etappp[i]);
      break;
    case 3:
      tree_ww_omega->GetEvent(index_omega[i]);
      break;
    default:
      cout << "GetEventWW: Wrong mode " << mode << endl;
      return;
    }
  } else{
    m_ww_tree->GetEvent(i);
  }
  m_dt = (m_t_sig - m_t_asc)*0.1*cm2ps;
  m_sz_rec /= 10.;
  m_sz_asc /= 10.;
//  m_bin = -m_bin;
//  m_bin_mc = -m_bin_mc;
  if(no_bkg) m_f_bkg = 0;
  return;
}

void GetEventGen(const int i){
  m_ww_tree->GetEvent(i);
  m_dt = (m_t_sig_mc - m_t_asc_mc)*0.1*cm2ps;
  return;
}

void GetSBEventWW(const int i){
  m_ww_sb_tree->GetEvent(i);
  m_dt = (m_t_sig - m_t_asc)*0.1*cm2ps;
  m_sz_rec /= 10.;
  m_sz_asc /= 10.;

//  if(!perftag){ m_flv = m_tag > 0 ? -1 : 1;}
  return;
}

//void SetPDFParams(const vector<double>& par){
//  if(m_svd != 2){
//    m_pdf_svd1->SetTauDm(par[0],par[1]);
//    m_pdf_svd1->SetSinCos(par[2],par[3]);
//    m_pdf_svd1->Set_f_ol_sgl(par[4]);
//    m_pdf_svd1->Set_f_ol_mlt(par[5]);
//    m_pdf_svd1->Set_sigma_ol(par[6]);
//  }
//  if(m_svd != 1){
//    m_pdf_svd2->SetTauDm(par[0],par[1]);
//    m_pdf_svd2->SetSinCos(par[2],par[3]);
//    m_pdf_svd2->Set_f_ol_sgl(par[7]);
//    m_pdf_svd2->Set_f_ol_mlt(par[8]);
//    m_pdf_svd2->Set_sigma_ol(par[9]);
//  }

////  m_pdf->Set_Srec(par[i],par[i+1]); i += 2;
////  m_pdf->Set_Smn_rec(par[i++]);
////  m_pdf->Set_Stl_rec(par[i++]);
////  m_pdf->Set_ftl_rec(par[i++]);

////  m_pdf->Set_Sasc(par[i],par[i+1]); i += 2;
////  m_pdf->Set_Smn_asc(par[i++]);
////  m_pdf->Set_Stl_asc(par[i++]);
////  m_pdf->Set_ftl_asc(par[i++]);

////  m_pdf->Set_fd_np_sgl(par[i],par[i+1]); i+=2;
////  m_pdf->Set_fp_np_sgl(par[i++]);
////  m_pdf->Set_tau_np_p_sgl(par[i],par[i+1]); i += 2;
////  m_pdf->Set_tau_np_n_sgl(par[i],par[i+1]); i += 2;

////  m_pdf->Set_fd_np_mlt(par[i],par[i+1]); i += 2;
////  m_pdf->Set_fd_np_st_mlt(par[i++]);
////  m_pdf->Set_fd_np_xi_mlt(par[i++]);
////  m_pdf->Set_fd_np_stxi_mlt(par[i++]);
////  m_pdf->Set_fp_np_mlt(par[i++]);
////  m_pdf->Set_fn_np_mlt(par[i++]);
////  m_pdf->Set_tau_np_p_mlt(par[i],par[i+1]); i += 2;
////  m_pdf->Set_tau_np_p_xi_mlt(par[i++]);
////  m_pdf->Set_tau_np_p_stxi_mlt(par[i++]);
////  m_pdf->Set_tau_np_n_mlt(par[i],par[i+1]); i += 2;
////  m_pdf->Set_tau_np_n_xi_mlt(par[i++]);
////  m_pdf->Set_tau_np_n_stxi_mlt(par[i++]);

//  return;
//}

void SetPDFParamsWW(const vector<double>& par){
  if(m_svd != 2){
    m_pdf_svd1->SetTauDm(par[0],0);
    m_pdf_svd1->Set_f_ol_sgl(par[1]);
    m_pdf_svd1->Set_f_ol_mlt(par[2]);
    m_pdf_svd1->Set_sigma_ol(par[3]);
  }
  if(m_svd != 1){
    m_pdf_svd2->SetTauDm(par[0],0);
    m_pdf_svd2->Set_f_ol_sgl(par[4]);
    m_pdf_svd2->Set_f_ol_mlt(par[5]);
    m_pdf_svd2->Set_sigma_ol(par[6]);
  }
  return;
}

void SetSBPDFParamsWW(const vector<double>& par){
  if(m_svd != 2){
      m_pdf_back_svd1_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd1_gg*par[0]);
    m_pdf_back_svd1_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd1_gg*par[0]);
      m_pdf_back_svd1_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd1_gg*par[0]);
    m_pdf_back_svd1_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd1_gg*par[0]);
      m_pdf_back_svd1_bb_gg->Set_f_otlr(par[2]);
    m_pdf_back_svd1_cont_gg->Set_f_otlr(par[3]);

      m_pdf_back_svd1_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd1_gg*par[0]);
    m_pdf_back_svd1_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd1_gg*par[0]);
      m_pdf_back_svd1_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd1_gg*par[0]);
    m_pdf_back_svd1_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd1_gg*par[0]);
      m_pdf_back_svd1_bb_ppp->Set_f_otlr(par[2]);
    m_pdf_back_svd1_cont_ppp->Set_f_otlr(par[3]);
  }

  if(m_svd != 1){
      m_pdf_back_svd2_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd2_gg*par[1]);
    m_pdf_back_svd2_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd2_gg*par[1]);
      m_pdf_back_svd2_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd2_gg*par[1]);
    m_pdf_back_svd2_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd2_gg*par[1]);
      m_pdf_back_svd2_bb_gg->Set_f_otlr(par[4]);
    m_pdf_back_svd2_cont_gg->Set_f_otlr(par[5]);

      m_pdf_back_svd2_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd2_gg*par[1]);
    m_pdf_back_svd2_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd2_gg*par[1]);
      m_pdf_back_svd2_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd2_gg*par[1]);
    m_pdf_back_svd2_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd2_gg*par[1]);
      m_pdf_back_svd2_bb_ppp->Set_f_otlr(par[4]);
    m_pdf_back_svd2_cont_ppp->Set_f_otlr(par[5]);
  }

  return;
}

void SetPDFParamsCPVWW(const vector<double>& par){
  if(m_svd != 2){
    m_pdf_svd1->SetSinCos(par[0],par[1]);
    m_pdf_svd1->Set_f_ol_sgl(par[2]);
    m_pdf_svd1->Set_f_ol_mlt(par[3]);
    m_pdf_svd1->Set_sigma_ol(par[4]);
  }
  if(m_svd != 1){
    m_pdf_svd2->SetSinCos(par[0],par[1]);
    m_pdf_svd2->Set_f_ol_sgl(par[5]);
    m_pdf_svd2->Set_f_ol_mlt(par[6]);
    m_pdf_svd2->Set_sigma_ol(par[7]);
  }
  return;
}

void SetPDFParamsExt(const vector<double>& par){
  if(m_svd != 2){
    // lifetime + cpv
    m_pdf_svd1->SetTauDm(par[0],par[1]);
    m_pdf_svd1->SetSinCos(par[2],par[3]);

    // Outlier
    m_pdf_svd1->Set_f_ol_sgl(par[4]);
    m_pdf_svd1->Set_f_ol_mlt(par[5]);
//    m_pdf_svd1->Set_sigma_ol(par[6]);

    // Background Scale
    if(!no_bkg){
        m_pdf_back_svd1_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd1_gg*par[10]);
      m_pdf_back_svd1_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd1_gg*par[10]);
        m_pdf_back_svd1_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd1_gg*par[10]);
      m_pdf_back_svd1_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd1_gg*par[10]);

        m_pdf_back_svd1_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd1_ppp*par[10]);
      m_pdf_back_svd1_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd1_ppp*par[10]);
        m_pdf_back_svd1_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd1_ppp*par[10]);
      m_pdf_back_svd1_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd1_ppp*par[10]);
    }
  }
  if(m_svd != 1){
    // lifetime + cpv
    m_pdf_svd2->SetTauDm(par[0],par[1]);
    m_pdf_svd2->SetSinCos(par[2],par[3]);

    // Outlier
    m_pdf_svd2->Set_f_ol_sgl(par[7]);
    m_pdf_svd2->Set_f_ol_mlt(par[8]);
//    m_pdf_svd1->Set_sigma_ol(par[9]);

    // Background Scale
    if(!no_bkg){
        m_pdf_back_svd2_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd2_gg*par[11]);
      m_pdf_back_svd2_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd2_gg*par[11]);
        m_pdf_back_svd2_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd2_gg*par[11]);
      m_pdf_back_svd2_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd2_gg*par[11]);

        m_pdf_back_svd2_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd2_ppp*par[11]);
      m_pdf_back_svd2_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd2_ppp*par[11]);
        m_pdf_back_svd2_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd2_ppp*par[11]);
      m_pdf_back_svd2_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd2_ppp*par[11]);
    }
  }
  for(int i=0; i<8; i++){
    Karr[i]  = par[12+i];
    Kbarr[i] = par[20+i];
    Carr[i]  = par[28+i];
    Sarr[i]  = par[36+i];
  }

  m_wr_tag_offset_svd1 = par[44];
  m_wr_tag_offset_svd2 = par[45];

  for(int j=0; j<modes_set.size(); j++){
    m_f_bb_offset[j]  = par[46+34*j];
    m_f_prt_offset[j] = par[47+34*j];
    for(int i=0; i<16; i++){
      m_Nsig_offset[j][i]  = par[48+34*j+i];
      m_Nbsig_offset[j][i] = par[48+34*j+i+16];
    }
  }

  return;
}

void SetSBTestPdfParams(const vector<double>& par){
  int i=0;
  if(m_svd != 2){
    m_pdf_back_svd1->Set_f_delta_sgl(par[i++]);
    m_pdf_back_svd1->Set_f_delta_mlt(par[i++]);
    m_pdf_back_svd1->Set_f_otlr(par[i++]);
  } else{ i += 3;}
  if(m_svd != 1){
    m_pdf_back_svd2->Set_f_delta_sgl(par[i++]);
    m_pdf_back_svd2->Set_f_delta_mlt(par[i++]);
    m_pdf_back_svd2->Set_f_otlr(par[i++]);
  }
  return;
}

void SetBkgPdfParams(const vector<double>& par){
  int i=0;
  if(m_svd != 2){
    if(m_type_flag){
      m_pdf_back_svd1->SetTau(par[i++]);
      m_pdf_back_svd1->Set_mu(par[i++]);
      m_pdf_back_svd1->Set_mu_delta(par[i++]);

      m_pdf_back_svd1->Set_f_delta_mlt(par[i++]);
      m_pdf_back_svd1->Set_f_tail_mlt(par[i++]);
      m_pdf_back_svd1->Set_S_main_mlt(par[i++]);
      m_pdf_back_svd1->Set_S_tail_mlt(par[i++]);

      m_pdf_back_svd1->Set_f_delta_sgl(par[i++]);
      m_pdf_back_svd1->Set_f_tail_sgl(par[i++]);
      m_pdf_back_svd1->Set_S_main_sgl(par[i++]);
      m_pdf_back_svd1->Set_S_tail_sgl(par[i++]);

      m_pdf_back_svd1->Set_f_otlr(par[i++]);
      m_pdf_back_svd1->Set_s_otlr(par[i++]);
    } else{
        m_pdf_back_svd1_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd1_gg*par[26]);
      m_pdf_back_svd1_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd1_gg*par[26]);
        m_pdf_back_svd1_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd1_gg*par[26]);
      m_pdf_back_svd1_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd1_gg*par[26]);
        m_pdf_back_svd1_bb_gg->Set_f_otlr(par[11]);
      m_pdf_back_svd1_cont_gg->Set_f_otlr(par[11]);

        m_pdf_back_svd1_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd1_ppp*par[26]);
      m_pdf_back_svd1_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd1_ppp*par[26]);
        m_pdf_back_svd1_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd1_ppp*par[26]);
      m_pdf_back_svd1_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd1_ppp*par[26]);
        m_pdf_back_svd1_bb_ppp->Set_f_otlr(par[11]);
      m_pdf_back_svd1_cont_ppp->Set_f_otlr(par[11]);
    }
  } else{
    i += 13;
  }
  if(m_svd != 1){
    if(m_type_flag){
      m_pdf_back_svd2->SetTau(par[0]);
      m_pdf_back_svd2->Set_mu(par[i++]);
      m_pdf_back_svd2->Set_mu_delta(par[i++]);

      m_pdf_back_svd2->Set_f_delta_mlt(par[i++]);
      m_pdf_back_svd2->Set_f_tail_mlt(par[i++]);
      m_pdf_back_svd2->Set_S_main_mlt(par[i++]);
      m_pdf_back_svd2->Set_S_tail_mlt(par[i++]);

      m_pdf_back_svd2->Set_f_delta_sgl(par[i++]);
      m_pdf_back_svd2->Set_f_tail_sgl(par[i++]);
      m_pdf_back_svd2->Set_S_main_sgl(par[i++]);
      m_pdf_back_svd2->Set_S_tail_sgl(par[i++]);

      m_pdf_back_svd2->Set_f_otlr(par[i++]);
      m_pdf_back_svd2->Set_s_otlr(par[i++]);
    } else{
        m_pdf_back_svd2_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd2_gg*par[26]);
      m_pdf_back_svd2_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd2_gg*par[26]);
        m_pdf_back_svd2_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd2_gg*par[26]);
      m_pdf_back_svd2_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd2_gg*par[26]);
        m_pdf_back_svd2_bb_gg->Set_f_otlr(par[23]);
      m_pdf_back_svd2_cont_gg->Set_f_otlr(par[23]);

        m_pdf_back_svd2_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd2_ppp*par[26]);
      m_pdf_back_svd2_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd2_ppp*par[26]);
        m_pdf_back_svd2_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd2_ppp*par[26]);
      m_pdf_back_svd2_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd2_ppp*par[26]);
        m_pdf_back_svd2_bb_ppp->Set_f_otlr(par[23]);
      m_pdf_back_svd2_cont_ppp->Set_f_otlr(par[23]);
    }
  } else{
    i += 12;
  }
  if(!m_type_flag){
    m_f_cont = par[25];
//    m_scale  = par[26];
  }
  return;
}

void GetMoments(const vector<double> v,double& mean, double& rms, double& merr){
  const int N = v.size();
  double sum = 0;
  double sumsq = 0;
  for(int i=0; i<N; i++){
    sum   += v[i];
    sumsq += v[i]*v[i];
  }
  mean = sum/N;
  rms  = (sumsq/N - mean*mean)*N/(N-1);
  merr = rms/sqrt(N);
  return;
}

void AnalyseDataSet1(const int SetNum = 0){
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      m_sig_map[k][j] = 0;
      m_bkg_map[k][j] = 0;
      m_fbkg_map[k][j] = 0;
    }
  }
  int NTot = SetNum ? m_NSig+m_NBkg : m_NSigTot+m_NBkgTot;
  const int NSig = SetNum ? m_NSig        : m_NSigTot;
  if(no_bkg) NTot = NSig;
  int BadEvents = 0;
  for(int i=0; i<NTot; i++){
    i<NSig ? GetEvent(i,SetNum,0) : GetEvent(i-NSig,SetNum,1);
    if(!perftag) m_flv = m_tag > 0 ? 1 : -1;
    i<NSig ? m_sig_map[flv_ind(m_flv)][bin_ind(m_bin)]++ : m_bkg_map[flv_ind(m_flv)][bin_ind(m_bin)]++;
  }
  if(no_bkg){
    m_fbkg_tot = 0;
    return;
  }
//  for(int i=0; i<m_NSigTot; i++){
//    sig_tree->GetEvent(i);
//    if(!perftag) m_flv = m_tag > 0 ? 1 : -1;
//    m_sig_map[flv_ind(m_flv)][bin_ind(m_bin)]++;
//  }
//  for(int i=0; i<m_NBkgTot; i++){
//    bkg_tree->GetEvent(i);
//    if(!perftag) m_flv = m_tag > 0 ? 1 : -1;
//    m_bkg_map[flv_ind(m_flv)][bin_ind(m_bin)]++;
//  }

  m_fbkg_tot = (double)(NTot-NSig)/((double)NTot);
  for(int j=0; j<16; j++){
    for(int k=0; k<2; k++){
      if(m_bkg_map[k][j]+m_sig_map[k][j] != 0){
        m_fbkg_map[k][j] = (double)m_bkg_map[k][j]/((double)m_bkg_map[k][j]+(double)m_sig_map[k][j]);
        if(true || !SetNum) cout << "fbkg[" << flv(k) << "," << bin(j) << "]: " << m_fbkg_map[k][j] << ", ";
      }
    }
    if(true || !SetNum) {if(m_bkg_map[0][j]+m_sig_map[0][j] + m_bkg_map[1][j]+m_sig_map[1][j] != 0) cout << endl;}
  }
  if(!SetNum){
    cout << "Data samples analized: " << endl;
    cout << "  NSig: " << m_NSigTot << endl;
    cout << "  NBkg: " << m_NBkgTot << endl;
    cout << "  NBad: " << BadEvents << endl;
    cout << "  fBkg: " << m_fbkg_tot << endl;
  }
  return;
}

void AnalyseDataSet2(const int SetNum = 0){
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      m_sig_map[k][j] = 0;
      m_bkg_map[k][j] = 0;
      m_fbkg_map[k][j] = 0;
    }
  }
  int NTot = SetNum ? m_NSig + m_NBkg : m_NGoodTot;
  int NSig = 0;
  for(int i=0; i<NTot; i++){
    GetEvent(i,SetNum,-1);
    if(!perftag) m_flv = m_tag > 0 ? -1 : 1;
    if(m_b0f != -1){
      m_sig_map[flv_ind(m_flv)][bin_ind(m_bin)]++;
      NSig++;
    } else{
      m_bkg_map[flv_ind(m_flv)][bin_ind(m_bin)]++;
    }
  }
  m_fbkg_tot = (double)(NTot-NSig)/((double)NTot);
  for(int j=0; j<16; j++){
    for(int k=0; k<2; k++){
      if(m_bkg_map[k][j]+m_sig_map[k][j] != 0){
        m_fbkg_map[k][j] = (double)m_bkg_map[k][j]/((double)m_bkg_map[k][j]+(double)m_sig_map[k][j]);
        if(true || !SetNum) cout << "fbkg[" << flv(k) << "," << bin(j) << "]: " << m_fbkg_map[k][j] << ", ";
      }
    }
    if(true || !SetNum) {if(m_bkg_map[0][j]+m_sig_map[0][j] + m_bkg_map[1][j]+m_sig_map[1][j] != 0) cout << endl;}
  }
  m_NSigTot = NSig;
  m_NBkgTot = 0;
  if(!SetNum){
    cout << "Data samples analized: " << endl;
    cout << "  NSig: " << NSig << endl;
    cout << "  NBkg: " << NTot-NSig << endl;
    cout << "  NBad: " << 0 << endl;
    cout << "  fBkg: " << m_fbkg_tot << endl;
  }
  return;
}

double PdfGen(const double& dt){
  const double tau = mm_btau;
  const double dm  = mm_dm;
  const double k   = K( m_bin_mc);
  const double kb  = K(-m_bin_mc);
  const double s   = S(m_bin_mc);
  const double c   = C(m_bin_mc);
  const int flv    = -m_flv_mc;
  const int xil    = 1;
  const double sin2beta = mm_sin2beta;
  const double cos2beta = mm_cos2beta;

  double pdf = k + kb;
  pdf += flv*(k-kb)*cos(dm*dt);
  pdf += 2*flv*xil*sqrt(k*kb)*sin(dm*dt)*(c*sin2beta+s*cos2beta);
  pdf *= exp(-abs(dt)/tau);
  const double norm = 2*tau*(k+kb+flv*(k-kb)/(1.+tau*dm*tau*dm));
//  const double norm = 0.5/tau*(k+kb);
  if(std::isnan(pdf) || std::isnan(norm) || pdf<0 || norm<0){
    cout << "!!! Bad gen PDF: " << pdf << ", norm: " << norm << endl;
    return 0;
  }
  return pdf/norm;
}

double Pdf(const double& dt, const int mode = 0){
  // mode  0 -> defalt
  // mode  1 -> only signal
  // mode  2 -> only background (sideband)
  // mode -1 -> AB fit
  double pdf = 0;
  const int svd = m_exp<30 ? 1 : 2;
  const double flv = perftag ? m_flv_mc :-m_tag;
  const int    bin = perfbin ? m_bin_mc : m_bin;
  const double s = sum_sigma(m_sz_rec,m_sz_asc);
  double f_bkg = 0, f_cont = 0.5;
  if(!no_bkg){
    if(!perfbin && !perftag){ f_bkg = m_f_bkg;        f_cont = m_f_cont_in_comb;}
    if(perfbin && !perftag){  f_bkg = m_f_bkg_bin_mc; f_cont = m_f_cont_in_comb_bin_mc;}
    if(!perfbin && perftag){  f_bkg = m_f_bkg_flv_mc; f_cont = m_f_cont_in_comb_flv_mc;}
    if(perfbin && perftag){   f_bkg = m_f_bkg_mc;     f_cont = m_f_cont_in_comb_mc;}
  }
  RbkgPdf* bkg_pdf_bb;
  RbkgPdf* bkg_pdf_cont;
  if(ppp_flag(Mode(mm_mode,mm_h0mode))){
      if(m_exp>30){
        bkg_pdf_bb = m_pdf_back_svd2_bb_ppp;
        bkg_pdf_cont = m_pdf_back_svd2_cont_ppp;
      } else{
        bkg_pdf_bb = m_pdf_back_svd1_bb_ppp;
        bkg_pdf_cont = m_pdf_back_svd1_cont_ppp;
      }
  } else{
    if(m_exp>30){
      bkg_pdf_bb = m_pdf_back_svd2_bb_gg;
      bkg_pdf_cont = m_pdf_back_svd2_cont_gg;
    } else{
      bkg_pdf_bb = m_pdf_back_svd1_bb_gg;
      bkg_pdf_cont = m_pdf_back_svd1_cont_gg;
    }
  }

//  cout << "Pdf mode " << mode << endl;
  if(svd == 1){
    if(mode != 2){
      m_pdf_svd1->SetAkCk(m_costhBcms,0.5*10.58);
      m_pdf_svd1->SetTag(flv);
      if(mode != -1){
        m_pdf_svd1->SetKKCS(K(bin),K(-bin),C(bin),S(bin));
        pdf = (1-f_bkg)*m_pdf_svd1->Pdf(dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc,true,no_interf);
        if(pdf<=0){
          cout << "SigPdf: SVD1 " << m_pdf_svd1->Pdf(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc) << endl;
          cout << " (" << m_dt << "," << m_ntrk_rec << "," << m_sz_rec << "," << m_chisq_rec << "," << m_ndf_rec << "," << m_ntrk_asc << "," << m_sz_asc << "," << m_chisq_asc << "," << m_ndf_asc << ")" << endl;
          cout << " (" << K(bin) << "," << K(-bin) << "," << C(bin) << "," << S(bin) << ")" << endl;
        } else{
//          cout << "Good sig Pdf SVD1: " << pdf << endl;
        }
      } else{
        pdf = (1-f_bkg)*m_pdf_svd1->PdfAB(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc);
        if(pdf<=0){
          cout << "SigPdf: AB SVD1 " << m_pdf_svd1->PdfAB(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc) << endl;
          cout << " (" << m_dt << "," << m_ntrk_rec << "," << m_sz_rec << "," << m_chisq_rec << "," << m_ndf_rec << "," << m_ntrk_asc << "," << m_sz_asc << "," << m_chisq_asc << "," << m_ndf_asc << ")" << endl;
        }
      }
    }
    if(!no_bkg){
      if(m_type_flag) pdf += f_bkg*m_pdf_back_svd1->Pdf(dt,s,m_ndf_asc);
      else            pdf += f_bkg*((1.-f_cont)*bkg_pdf_bb->Pdf(dt,s,m_ndf_asc)+f_cont*bkg_pdf_cont->Pdf(dt,s,m_ndf_asc));
    }
  } else{
    if(mode != 2){
      m_pdf_svd2->SetAkCk(m_costhBcms,0.5*10.58);
      m_pdf_svd2->SetTag(flv);
      if(mode != -1){
        m_pdf_svd2->SetKKCS(K(bin),K(-bin),C(bin),S(bin));
        pdf = (1-f_bkg)*m_pdf_svd2->Pdf(dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc,true,no_interf);
        if(pdf<=0){
          cout << "SigPdf: SVD2 " << m_pdf_svd1->Pdf(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc) << endl;
          cout << " (" << m_dt << "," << m_ntrk_rec << "," << m_sz_rec << "," << m_chisq_rec << "," << m_ndf_rec << "," << m_ntrk_asc << "," << m_sz_asc << "," << m_chisq_asc << "," << m_ndf_asc << ")" << endl;
          cout << " (" << K(bin) << "," << K(-bin) << "," << C(bin) << "," << S(bin) << ")" << endl;
        } else{
//          cout << "Good sig Pdf SVD2: " << pdf << endl;
        }
      } else{
        pdf = (1-f_bkg)*m_pdf_svd2->PdfAB(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc);
        if(pdf<=0){
          cout << "SigPdf: AB SVD1 " << m_pdf_svd1->PdfAB(m_dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc) << endl;
          cout << " (" << m_dt << "," << m_ntrk_rec << "," << m_sz_rec << "," << m_chisq_rec << "," << m_ndf_rec << "," << m_ntrk_asc << "," << m_sz_asc << "," << m_chisq_asc << "," << m_ndf_asc << ")" << endl;
        }
      }
    }
    if(!no_bkg){
      if(m_type_flag) pdf += f_bkg*m_pdf_back_svd2->Pdf(dt,s,m_ndf_asc);
      else            pdf += f_bkg*((1.-m_f_cont_in_comb)*bkg_pdf_bb->Pdf(dt,s,m_ndf_asc)+m_f_cont_in_comb*bkg_pdf_cont->Pdf(dt,s,m_ndf_asc));
    }
  }
  if(std::isnan(pdf) || pdf<= 0){
    cout << "pdf: " << pdf << ", dt = " << dt << ", fbkg = " << f_bkg << ", fcnt = " << m_f_cont_in_comb << endl;
    cout << " sz_sig = " << m_sz_rec << ", sz_asc = " << m_sz_asc << endl;
    nega_pdf_flag = true;
  }
  return pdf;
}

int GetGoodTTree(TTree* tree, const int bin = 0, const int flv = 0, const bool SB = false){
  stringstream out;
  out.str("");
  out << "1";
  if(m_svd == 2)            out << " && exp>30";
  else if(m_svd == 1)       out << " && exp<30";
  if(bin &&  ABfit)         out << " && bin == " << bin;
  if(bin && !ABfit)         out << " && abs(bin) == " << bin;
  if(flv && perftag){       out << " && flv == " << flv;}
  if(flv ==  1 && !perftag) out << " && tag_LH < 0";
  if(flv == -1 && !perftag) out << " && tag_LH > 0";
//  out << " && mode == "   << Mode(m_mode);
//  out << " && h0mode == " << h0Mode(m_mode);
  if(no_bkg) out << " && sigarea";

  if(sgl_asc) out << " && ndf_asc==0";
  if(mlt_asc) out << " && ndf_asc>0";
  if(no_np)   out << " && !nptag";

  if(m_mode && !SB && !m_ebeb){
    const string sig_cut = out.str() + string(" && (b0f == 1 || b0f == 5 || b0f == 10)");
    m_sig_tree = tree->CopyTree(sig_cut.c_str());
    m_NSigTot = m_sig_tree->GetEntries();
    const string bkg_cut = out.str() + string(" && (b0f == -1 || b0f>0) && b0f != 1 && b0f != 5 && b0f != 10");
    m_bkg_tree = tree->CopyTree(bkg_cut.c_str());
    m_NBkgTot = m_bkg_tree->GetEntries();
    AnalyseDataSet1(0);
  } else if(SB){
    out << " && b0f != 1 && b0f != 5 && b0f != 10";
    m_good_sb_tree = tree->CopyTree(out.str().c_str());
    m_NSBTot = m_good_sb_tree->GetEntries();
    if(!m_data) m_bkg_tree = m_test_tree->CopyTree(out.str().c_str());
  } else{
    m_good_tree = tree->CopyTree(out.str().c_str());
    m_NGoodTot = m_good_tree->GetEntries();
    if(m_ebeb) AnalyseDataSet2();
  }
  return 0;
}

int GetGoodWWTTree(TChain* tree, TTree* outtree){
  stringstream out;
  out.str("");
  out << "b0f!=0";
  out << " && sigarea";
  out << " && sz_sig>0.0001 && sz_asc>0.0001";
  if(m_fitbin && perfbin && ABfit) out << " && bin_mc == " << m_fitbin;
  if(m_fitbin &&!perfbin && ABfit) out << " && bin == "    << m_fitbin;
  if(m_fitbin && perfbin &&!ABfit) out << " && abs(bin_mc) == " << m_fitbin;
  if(m_fitbin &&!perfbin &&!ABfit) out << " && abs(bin) == "    << m_fitbin;
  if(m_fitflv && perftag && ABfit) out << " && flv_mc == " << m_fitflv;
  if(m_fitflv &&!perftag && ABfit) out << " && flv == "    << m_fitflv;
  if(m_fitflv && perftag &&!ABfit) out << " && flv_mc == " << m_fitflv;
  if(m_fitflv &&!perftag &&!ABfit) out << " && flv == "    << m_fitflv;
  if(no_bkg || sigmc) out << " && (b0f == 1 || b0f == 5 || b0f == 10)";
  if(!perftag && !no_interf && !m_genfit) out << " && abs(tag_LH)>0.1";
  if(m_svd == 2) out << " && exp>30";
  if(m_svd == 1) out << " && exp<30";
  cout << out.str() << endl;
  cout << tree->GetEntries() << endl;
  outtree = tree->CopyTree(out.str().c_str());

//  if(!no_bkg){
//    out.str("");
//    out << "mbc>5.23 && mbc<5.26 || mbc>5.25 && de>0.12";
//    m_ww_sb_tree = tree->CopyTree(out.str().c_str());
//  }

  return 0;
}

//class pdfFcn : public FCNBase{
//public:
//  pdfFcn(const int Set=0){
//    SetNum = Set;
//    theErrorDef = 1;
//    if(!m_mode){           NTot = m_NGoodTot;          NSig = NTot;}
//    else if(!full_ds_fit){ NTot = m_NSig+m_NBkg;       NSig = m_NSig;}
//    else{                  NTot = m_NSigTot+m_NBkgTot; NSig = m_NSigTot;}
//    if(no_bkg) NTot = NSig;
//    if(m_ebeb) NSig = NTot;
//    cout << NTot << " events to process, fbkg = " << m_fbkg_tot << ", Set = " << Set << endl;
//  }
//  ~pdfFcn() {}

//  double Up() const {return theErrorDef;}
//  double operator()(const vector<double>& par) const {
//    double loglh = 0;
//    double pdf;
//    SetPDFParams(par);
//    for(int i=0; i<NTot; i++){
//      i<NSig ? GetEvent(i,SetNum,0) : GetEvent(i-NSig,SetNum,1);
//      pdf = Pdf(m_dt);
//      if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
//      else{
//        if(!std::isnan(pdf)) loglh +=-10000*pdf;
//        else            loglh += 10000;
//        cout << "pdf = " << pdf; PrintEvent();
//      }
//    }
//    if(full_ds_fit){
//      cout << "loglh: " << loglh << ", tau: " << par[0] << ", fbkg: " << m_f_bkg;
//      if(!no_interf) cout << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3];
//      cout << endl;
//    }
//    return loglh;
//  }
//private:
//  double theErrorDef;
//  int NTot,NSig;
//  int SetNum;
//};

void SetPDFParamsGen(const vector<double>& par){
  mm_btau = par[0];
  mm_dm   = par[1];
  mm_sin2beta = par[2];
  mm_cos2beta = par[3];

  for(int i=0; i<8; i++){
    Karr[i]  = par[4+i];
    Kbarr[i] = par[4+8+i];
    Carr[i]  = par[4+8+8+i];
    Sarr[i]  = par[4+8+8+8+i];
  }
  return;
}

class pdfFcnGen : public FCNBase{
public:
  pdfFcnGen(void){
    theErrorDef = 1;
    NTot = m_ww_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnGen() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double loglh = 0;
    double pdf;
    SetPDFParamsGen(par);
    for(int i=0; i<NTot; i++){
      GetEventGen(i);
      pdf = PdfGen(m_dt);
      if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
      else{
        if(!std::isnan(pdf)) loglh +=-10000*pdf;
        else            loglh += 10000;
        cout << "pdf = " << pdf; PrintEvent();
      }
    }
    cout << "loglh: " << loglh << ", tau: " << par[0] << ", dm = " << par[1];
    cout << ", sin: " << par[2] << ", cos: " << par[3];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

class pdfFcnWW : public FCNBase{
public:
  pdfFcnWW(const int smpl){
    theErrorDef = 1;
    NIni = 0;
    NTot = m_toyfit ? 0 : m_ww_tree->GetEntries();
    m_smpl = smpl;
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnWW() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double loglh = 0;
    double pdf;
    SetPDFParamsWW(par);
    if(m_toyfit){
      for(int i=0; i<4; i++){
        for(int j=0; j<toy_composition[i]; j++){
          GetEventWW(m_smpl*toy_composition[i] + j,i);
          pdf = Pdf(m_dt,1);
          if(!std::isnan(pdf) && pdf>0){
            loglh += -2*TMath::Log(pdf);
          } else{
            if(!std::isnan(pdf)) loglh +=-10000*pdf;
            else                 loglh += 10000;
            cout << "pdf = " << pdf; PrintEvent();
          }
        }
      }
    } else{
      for(int i=NIni; i<NIni+NTot; i++){
        GetEventWW(i);
        pdf = no_bkg ? Pdf(m_dt,1) : Pdf(m_dt);
        if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
        else{
          if(!std::isnan(pdf)) loglh +=-10000*pdf;
          else            loglh += 10000;
          cout << "pdf = " << pdf; PrintEvent();
        }
      }
    }
    cout << "loglh: " << loglh << ", tau: " << par[0];
//    cout << ", scale1: " << par[1] << ", scale2: " << par[2];
//    cout << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot,NIni;
  int m_smpl;
};

class pdfFcnSBWW : public FCNBase{
public:
  pdfFcnSBWW(void){
    theErrorDef = 1;
    NTot = m_ww_sb_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnSBWW() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double loglh = 0;
    double pdf;
    SetSBPDFParamsWW(par);
    for(int i=0; i<NTot; i++){
      GetSBEventWW(i);
      pdf = Pdf(m_dt,2);
      if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
      else{
        if(!std::isnan(pdf)) loglh +=-10000*pdf;
        else            loglh += 10000;
        cout << "pdf = " << pdf; PrintEvent();
      }
    }
    cout << "loglh: " << loglh << ", scale1: " << par[0];
    cout << ", scale2: " << par[1];
//    cout << ", dm = " << par[1] << ", sin: " << par[2] << ", cos: " << par[3];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
};

class pdfFcnCPVWW : public FCNBase{
public:
  pdfFcnCPVWW(const int smpl = 0){
    theErrorDef = 1;
    NIni = 0;
    NTot = m_toyfit ? 0 : m_ww_tree->GetEntries();
    m_smpl = smpl;
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnCPVWW() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double loglh = 0;
    double pdf;
    SetPDFParamsCPVWW(par);
    if(m_toyfit){
      for(int i=0; i<4; i++){
        for(int j=0; j<toy_composition[i]; j++){
          GetEventWW(m_smpl*toy_composition[i] + j,i);
          if(m_svd != 2) m_pdf_svd1->SetXi(Xil(mm_mode,mm_h0mode));
          if(m_svd != 1) m_pdf_svd2->SetXi(Xil(mm_mode,mm_h0mode));
          pdf = Pdf(m_dt,1);
          if(!std::isnan(pdf) && pdf>0){
            loglh += -2*TMath::Log(pdf);
          } else{
            if(!std::isnan(pdf)) loglh +=-10000*pdf;
            else                 loglh += 10000;
            cout << "pdf = " << pdf; PrintEvent();
          }
        }
      }
    } else{
      for(int i=NIni; i<NIni+NTot; i++){
        GetEventWW(i);
        if(m_svd != 2){
          m_pdf_svd1->SetXi(Xil(mm_mode,mm_h0mode));
          if(false && !m_line_test && Beta(mm_mode,mm_h0mode) != 23) m_pdf_svd1->SetSinCos(par[0],-par[1]);
          else                              m_pdf_svd1->SetSinCos(par[0],par[1]);
        }
        if(m_svd != 1){
          m_pdf_svd2->SetXi(Xil(mm_mode,mm_h0mode));
          if(false && !m_line_test && Beta(mm_mode,mm_h0mode) != 23) m_pdf_svd2->SetSinCos(par[0],-par[1]);
          else                              m_pdf_svd2->SetSinCos(par[0],par[1]);
        }
        pdf = no_bkg ? Pdf(m_dt,1) : Pdf(m_dt);
        if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
        else{
          if(!std::isnan(pdf)) loglh +=-10000*pdf;
          else            loglh += 10000;
          cout << "pdf = " << pdf; PrintEvent();
        }
      }
    }
    cout << "loglh: " << loglh << ", sin: " << par[0] << ", cos: " << par[1];
//    cout << ", scale1: " << par[2] << ", scale2: " << par[3];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot,NIni;
  int m_smpl;
};

class pdfFcnExt : public FCNBase{
public:
  pdfFcnExt(const int Nini = 0, const int NumEve = 0){
    theErrorDef = 1;
    if(!Nini && !NumEve){
      NIni = 0;
      NTot = m_ww_tree->GetEntries();
    } else{
      NIni = Nini;
      NTot = NumEve;
    }
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnExt() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    double pdf;
//    cout << "Size: " << par.size() << endl;
    SetPDFParamsExt(par);
    double Nuisance = 0;
    Nuisance += pow((par[0]-m_btau)/m_btau_err,2);
    Nuisance += pow((par[1]-m_dm)/m_dm_err,2);
    Nuisance += pow((par[10]-m_scale1)/(3.*m_scale1_err),2);
    Nuisance += pow((par[11]-m_scale2)/(3.*m_scale2_err),2);
    // ** dE-Mbc fit uncertainties ** //
    for(int j=0; j<modes_set.size(); j++){
      if(fbbErrs[j]>0)  Nuisance += pow(par[46+34*j]/fbbErrs[j],2);
      if(fprtErrs[j]>0) Nuisance += pow(par[47+34*j]/fprtErrs[j],2);
      for(int i=0; i<16; i++){
        if(NsigErrs[j][i]>0)  Nuisance += pow(par[48+34*j+i]/NsigErrs[j][i],2);
        else cout << "NsigErr " << j << " " << i << " = " << NsigErrs[j][i] << endl;
        if(NbsigErrs[j][i]>0) Nuisance += pow(par[48+34*j+i+16]/NbsigErrs[j][i],2);
        else cout << "NbsigErr " << j << " " << i << " = " << NbsigErrs[j][i] << endl;
      }
    }
    Nuisance += pow(par[44]/err_wt_svd1,2);
    Nuisance += pow(par[45]/err_wt_svd2,2);
    // ** ** //
    for(int i=0; i<8; i++){
      Nuisance += pow((par[12+i]-K0[i])/Karr_err[i],2) + pow((par[20+i]-Kb0[i])/Kbarr_err[i],2);
//      Nuisance += 10000.*(Carr[i]-C0[i])*(Carr[i]-C0[i])/(C_stat_err[i]*C_stat_err[i]);//+C_syst_err[i]*C_syst_err[i]);
//      Nuisance += 10000.*(Sarr[i]-S0[i])*(Sarr[i]-S0[i])/(S_stat_err[i]*S_stat_err[i]);//+S_syst_err[i]*S_syst_err[i]);
      Nuisance += (Carr[i]-C0[i])*(Carr[i]-C0[i])*CS_inv_covMtx[i][i];
      Nuisance += (Sarr[i]-S0[i])*(Sarr[i]-S0[i])*CS_inv_covMtx[i+8][i+8];
      for(int j=0; j<8; j++){ if(i!=j){
        Nuisance += 0.5*(Carr[i]-C0[i])*(Carr[j]-C0[j])*CS_inv_covMtx[i][j];
        Nuisance += 0.5*(Sarr[i]-S0[i])*(Sarr[j]-S0[j])*CS_inv_covMtx[i+8][j+8];
        Nuisance += (Sarr[i]-S0[i])*(Carr[j]-C0[j])*CS_inv_covMtx[i+8][j];
//        cout << Nuisance << " " << Sarr[i]-S0[i] << " " << CS_inv_covMtx[i+8][j] << endl;
      }}
    }

    double loglh = Nuisance;

    for(int i=NIni; i<NIni+NTot; i++){
      GetEventWW(i);
      calc_fbkg_fcnt();

      if(m_svd != 2){
        m_pdf_svd1->SetXi(Xil(mm_mode,mm_h0mode));
        if(Beta(mm_mode,mm_h0mode) != 23) m_pdf_svd1->SetSinCos(par[2],-par[3]);
        else                              m_pdf_svd1->SetSinCos(par[2],par[3]);
      }
      if(m_svd != 1){
        m_pdf_svd2->SetXi(Xil(mm_mode,mm_h0mode));
        if(Beta(mm_mode,mm_h0mode) != 23) m_pdf_svd2->SetSinCos(par[2],-par[3]);
        else                              m_pdf_svd2->SetSinCos(par[2],par[3]);
      }
      pdf = no_bkg ? Pdf(m_dt,1) : Pdf(m_dt);
      if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
      else{
        if(!std::isnan(pdf)) loglh +=-10000*pdf;
        else            loglh += 10000;
        cout << "pdf = " << pdf; PrintEvent();
        continue;
      }

    }
    cout << "loglh: " << loglh << ", Nuis: " << Nuisance;
    cout << ", tau: " << par[0] << ", dm:" << par[1];
    cout << ", sin: " << par[2] << ", cos: " << par[3];
    cout << ", C1: " << par[28] << ", S1: " << par[36];
    cout << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot,NIni;
};

class pdfFcnSingle : public FCNBase{
public:
  pdfFcnSingle(void){
    theErrorDef = 1;
    NTot = m_ww_tree->GetEntries();
    cout << NTot << " events to process";//, flv: " << Flv << ", bin: " << Bin << endl;
    cout << endl;
  }
  ~pdfFcnSingle() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    const double& btau = par[0];
    const double& dm   = par[1];
    const double& A    = par[2];
    const double& B    = par[3];
    double loglh = 0;
    double pdf;
    m_pdf_svd1->SetTauDm(btau,dm);
    m_pdf_svd1->SetAB(A,B);
    m_pdf_svd2->SetTauDm(btau,dm);
    m_pdf_svd2->SetAB(A,B);
    for(int i=0; i<NTot; i++){
      GetEventWW(i);
      pdf = Pdf(m_dt,ABFIT);
      if(!std::isnan(pdf) && pdf>0){ loglh += -2*TMath::Log(pdf);}
//      else{ cout << "pdfFcnSingle: pdf = " << pdf << endl;}
    }
    cout << "loglh: " << loglh << ", tau: " << par[0] << ", dm = " << par[1] << ", A: " << par[2] << ", B: " << par[3] << endl;
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;//,NSig;
//  int SetNum;
//  int Flv,Bin;
};

class pdfFcnBkg : public FCNBase{
public:
  pdfFcnBkg(const int type = 2){
    // type 2 -> Sideband
    // type 1 -> SB test
    m_type = type;
    theErrorDef = 1;
    if(m_type == 2) NTot = m_NSBTot;
    else            NTot = m_bkg_tree->GetEntries();
    cout << NTot << " events to process" << endl;
  }
  ~pdfFcnBkg() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    if(m_type == 2) SetBkgPdfParams(par);
    else            SetSBTestPdfParams(par);
    double loglh = 0;
    double pdf;
    for(int i=0; i<NTot; i++){
      GetEvent(i,0,m_type);
      pdf = Pdf(m_dt,BACKGROUND);
//      if(m_type == 1) cout << m_dt << " " << pdf << endl;
      if(!std::isnan(pdf) && pdf>0) loglh += -2*TMath::Log(pdf);
//      else cout << "pdf is nan or zero " << pdf << " " << sum_sigma(m_sz_rec,m_sz_asc) << endl;
    }
    cout << "lh: " << loglh;
    if(m_type == 2){
      if(m_type_flag){
        for(int i=0; i<10; i++) cout << " " << par[i];
        cout << endl << "           ";
        for(int i=11; i<20; i++) cout << " " << par[i];
        cout << endl;
      } else{
        cout << " " << par[11];
        cout << " " << par[23];
        cout << " " << par[25] << " " << par[26] << endl;
      }
    } else{
      for(int i=0; i<6; i++) cout << " " << par[i];
      cout << endl;
    }
    return loglh;
  }
private:
  double theErrorDef;
  int NTot;
  int m_type;
};

int draw_fit_results(const double vals[][2],const double errs[][2], const double* Aref, const double* Bref){
  const int NFits = 16;
  double Fits[NFits],FitsErr[NFits];
  double A[NFits],AErr[NFits], dA[NFits];
  double B[NFits],BErr[NFits], dB[NFits];
  double ARef[NFits],BRef[NFits];

  for(int i=0; i<NFits; i++){
    if(NFits == 16) { Fits[i] = bin(i); FitsErr[i] = 0;}
    else            { Fits[i] = i+1; FitsErr[i] = 0;}
    A[i]   = vals[i][0];     AErr[i] = errs[i][0];
    B[i]   = vals[i][1];     BErr[i] = errs[i][1];
    ARef[i]= Aref[i];        BRef[i] = Bref[i];
    dA[i]  = A[i] - ARef[i]; dB[i]   = B[i] - BRef[i];
  }
  stringstream out;
  out.str("");
  string label(""),rootname,pngname;
  if(m_svd == 2) label += string("SVD2");
  if(m_svd == 1) label += string("SVD1");
  label += get_label(m_mode);
  if(m_flv == -1) label += string(", anti-B^{0}");
  if(sgl_asc) label += string(", sgl");
  if(mlt_asc) label += string(", mlt");
  if(no_bkg)  label += string(", no Bkg");
  if(perftag) label += string(", perf tag");

  out.str("");
  out << "A" << label;
  TMultiGraph* mgA = new TMultiGraph("mgA",out.str().c_str());
  out.str("");
  out << "B" << label;
  TMultiGraph* mgB = new TMultiGraph("mgB",out.str().c_str());
  TGraphErrors* grA = new TGraphErrors(NFits,Fits,A,FitsErr,AErr);
  grA->SetMarkerStyle(20);
  grA->SetMarkerSize(1.);
  grA->SetMarkerColor(kBlue);
  TGraph* grAref = new TGraphErrors(NFits,Fits,ARef);
  grAref->SetMarkerStyle(20);
  grAref->SetMarkerSize(1.);
  grAref->SetMarkerColor(kRed);

  TGraphErrors* grdA = new TGraphErrors(NFits,Fits,dA,FitsErr,AErr);
  grdA->SetMarkerStyle(20);
  grdA->SetMarkerSize(1.);
  grdA->SetMarkerColor(kBlue);

  TGraphErrors* grB = new TGraphErrors(NFits,Fits,B,FitsErr,BErr);
  grB->SetMarkerStyle(20);
  grB->SetMarkerSize(1.);
  grB->SetMarkerColor(kBlue);
  mgA->Add(grA);
  mgA->Add(grAref);
  TGraph* grBref = new TGraphErrors(NFits,Fits,BRef);
  grBref->SetMarkerStyle(20);
  grBref->SetMarkerSize(1.);
  grBref->SetMarkerColor(kRed);
  mgB->Add(grB);
  mgB->Add(grBref);

  TGraphErrors* grdB = new TGraphErrors(NFits,Fits,dB,FitsErr,BErr);
  grdB->SetMarkerStyle(20);
  grdB->SetMarkerSize(1.);
  grdB->SetMarkerColor(kBlue);

  TCanvas* c1 = new TCanvas("c1","c1",800,500);
  c1->cd();
  c1->Draw();
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.21,0.49,0.99);
  pad1->Draw();
  pad1->SetGrid();
  TPad* pad2 = new TPad("pad2","pad2",0.51,0.21,0.99,0.99);
  pad2->Draw();
  pad2->SetGrid();
  TPad* pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.19);
  pad3->Draw();
  pad3->SetGrid();
  TPad* pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.19);
  pad4->Draw();
  pad4->SetGrid();

  pad1->cd();
  grAref->GetXaxis()->SetTitle("Dalitz bin");
  grAref->GetXaxis()->SetTitleSize(0.05);
  grAref->GetXaxis()->SetTitleOffset(0.9);
  grAref->GetYaxis()->SetTitleSize(0.05);
  mgA->Draw("ap");
  pad2->cd();
  grBref->GetXaxis()->SetTitle("Dalitz bin");
  grBref->GetXaxis()->SetTitleSize(0.05);
  grBref->GetXaxis()->SetTitleOffset(0.9);
  grBref->GetYaxis()->SetTitleSize(0.05);
  mgB->Draw("ap");

  TLine* zeroline = new TLine(-8,0,8,0);
  zeroline->SetLineColor(kRed);
  zeroline->SetLineWidth(2);

  pad3->cd();
  grdA->GetYaxis()->SetRangeUser(-0.3,0.3);
  grdA->GetXaxis()->SetTitle("Dalitz bin");
  grdA->GetXaxis()->SetTitleSize(0.05);
  grdA->GetXaxis()->SetTitleOffset(0.9);
  grdA->GetYaxis()->SetTitleSize(0.05);
  grdA->Draw("ap");
  zeroline->Draw();

  pad4->cd();
  grdB->GetYaxis()->SetRangeUser(-0.3,0.3);
  grdB->GetXaxis()->SetTitle("Dalitz bin");
  grdB->GetXaxis()->SetTitleSize(0.05);
  grdB->GetXaxis()->SetTitleOffset(0.9);
  grdB->GetYaxis()->SetTitleSize(0.05);
  grdB->Draw("ap");
  zeroline->Draw();

  c1->Update();
  out.str("");
  out << "pics/ABFit" << m_mode;
  if(mlt_asc)    out << "_mlt";
  if(sgl_asc)    out << "_sgl";
  if(m_svd == 2) out << "_svd2";
  if(m_svd == 1) out << "_svd1";
  if(no_bkg)     out << "_nobkg";
  if(perftag)    out << "_perftag";
  if(m_flv == -1)out << "_b0b";
  rootname = out.str() + string(".root");

  c1->Print(rootname.c_str());
  pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  return 0;
}

void Draw_BBbar_CPV(const int SetNum = 0){
  int NTot,NSig;
  if(!m_mode)         { NTot = m_NGoodTot;          NSig = NTot;}
  else if(full_ds_fit){ NTot = m_NSigTot+m_NBkgTot; NSig = m_NSigTot;}
  else{                 NTot = m_NSig+m_NBkg;       NSig = m_NSig;}

  const double ddt = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  int Nev[2] = {0,0};
  for(int j=0; j<NDots; j++){
    dt_arr[j] = dtmin+(j+0.5)*ddt;
    if(include_dt0) dt_arr[j] += -0.0224885;
  }
  TH1I* dh[2];
  for(int i=0; i<2; i++){
    dh[i] = !i ? new TH1I("dh","dh",NBins,dtmin,dtmax) : new TH1I("dhb","dhb",NBins,dtmin,dtmax);
    dh[i]->SetMarkerStyle(21);
    !i ? dh[i]->SetMarkerColor(kBlue) : dh[i]->SetMarkerColor(kRed);
    dh[i]->SetMarkerSize(1.1);
  }
  double pdf_arr[2][NDots];
  double norm[2];

  for(int k=0; k<2; k++){
    norm[k] = 0;
    for(int j=0; j<NDots; j++){
      pdf_arr[k][j] = 0;
    }
  }
  for(int l=0; l<NTot; l++){
    l<NSig ? GetEvent(l,SetNum,0) : GetEvent(l-NSig,SetNum,1);
    dh[flv_ind(m_flv)]->Fill(m_dt);
    Nev[flv_ind(m_flv)]++;
    for(int k=0; k<2; k++){
      for(int j=0; j<NDots; j++){
        pdf_arr[k][j] += Pdf(dt_arr[j]);
      }
    }
  }
  for(int k=0; k<2; k++){
    for(int j=0; j<NDots; j++){
      pdf_arr[k][j] *= ddtb;
      norm[k] += pdf_arr[k][j];
    }
    cout << "Nev[" << k << "] = " << Nev[k] << ", ";
    cout << "norm[" << k << "] = " << norm[k] << endl;
  }

  TGraph* gr = new TGraph(NDots,dt_arr,pdf_arr[0]);
  gr->SetMarkerStyle(kDot);
  gr->SetMarkerSize(1.5);
  gr->SetLineWidth(2);

  TGraph* grb = new TGraph(NDots,dt_arr,pdf_arr[1]);
  grb->SetMarkerStyle(kDot);
  grb->SetMarkerSize(1.5);
  grb->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  dh[0]->Draw("e");
  dh[1]->Draw("e,same");
  gr->Draw("same");
  grb->Draw("same");

  c1->Print("pics/fit_no_int.png");
  c1->Print("pics/fit_no_int.root");

  system("display pics/fit_no_int.png &");
  return;
}

void Draw_NoTag(const int NFreePar,const int SetNum = 0){
  int NTot,NSig;
  if(!m_mode)         { NTot = m_NGoodTot;          NSig = NTot;}
  else if(full_ds_fit){ NTot = m_NSigTot+m_NBkgTot; NSig = m_NSigTot;}
  else{                 NTot = m_NSig+m_NBkg;       NSig = m_NSig;}
  if(no_bkg) NTot = NSig;
  const double ddt = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];

  for(int j=0; j<NDots; j++){
    dt_arr[j] = dtmin+(j+0.5)*ddt;
    if(include_dt0) dt_arr[j] += -0.0224885;
  }
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }

  double Norm = 0;
  double pdf_array[NDots];

  stringstream out;
  out.str("");
  out << "Lifetime fit";
  if(m_svd == 2)      out << ", SVD2";
  else if(m_svd == 1) out << ", SVD1";
  out << get_label(m_mode);
  if(sgl_asc) out << " (sgl)";
  if(mlt_asc) out << " (mlt)";

  gStyle->SetStatW(0.4); gStyle->SetStatH(0.3);
  TH1I* DH = new TH1I("DH",out.str().c_str(),NBins,dtmin,dtmax);
  DH->SetMarkerStyle(20);
  DH->SetMarkerColor(kBlack);
  DH->SetLineColor(kBlack);
  DH->SetMarkerSize(1.1);

  for(int i=0; i<NDots; i++) pdf_array[i] = 0;
  for(int j=0; j<NTot; j++){
    j<NSig ? GetEvent(j,SetNum,0) : GetEvent(j-NSig,SetNum,1);
    DH->Fill(m_dt);
    for(int i=0; i<NDots; i++){
      double pdfval;
      pdfval = Pdf(dt_arr[i]);
      if(!std::isnan(pdfval)){ pdf_array[i] += pdfval;}
      else{
        cout << "pdf is nan: " << dt_arr[i] << " " << m_ntrk_rec << " " << m_sz_rec << " " << m_chisq_rec << " " << m_ndf_rec << " " << m_ntrk_asc << " " << m_sz_asc << " " << m_chisq_asc << " " << m_ndf_asc << endl;
      }
    }
  }

  double pull_array[NBins];
  const int nBD = NDots/NBins;
  double chisq = 0;
  int NDF = 0;
  for(int i=0; i<NDots; i++){
    pdf_array[i] *= ddtb;
    Norm += pdf_array[i];
    if(!(i%nBD)){
      const int bin = i/nBD;
      const int bin_content = DH->GetBinContent(bin+1);
//      cout << "bin " << bin << " " << bin_content << " " << pdf_array[i] << endl;
      pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
      if(bin_content>5){
        chisq += pull_array[bin]*pull_array[bin]; NDF++;
      }
    }
  }
  cout << "Norm = " << Norm/nBD << ", Nevents = " << NTot << endl;
  chisq /= NDF;
  cout << "chi2/n.d.f. = " << chisq << endl;

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetLineWidth(2);
  GR->SetMarkerColor(kBlue);
  GR->SetLineColor(kBlue);

  TGraphErrors* GRP = new TGraphErrors(NBins,dt_barr,pull_array,dt_barr_err,dt_pull_err);
  GRP->SetMarkerStyle(20);
  GRP->SetMarkerSize(1.2);
  GRP->SetMarkerColor(kBlack);
  GRP->SetLineColor(kBlack);

  TCanvas* c1 = new TCanvas("c1","c1",600,800);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.19);
  pad1->Draw(); pad2->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  DH->GetXaxis()->SetTitle("#Deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.06);
  DH->GetXaxis()->SetTitleOffset(0.75);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->Draw("e1");

  TPaveText *pt = new TPaveText(0.60,0.53,0.98,0.6,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt->AddText(out.str().c_str());
  pt->Draw();

  GR->Draw("same");

  pad2->cd();
  GRP->GetYaxis()->SetRangeUser(-5,5);
  GRP->GetXaxis()->SetRangeUser(dtmin,dtmax);
  GRP->Draw("ap");
  TLine* zeroline = new TLine(dtmin,0,dtmax,0);
  zeroline->SetLineColor(kBlue);
  zeroline->SetLineWidth(2);
  zeroline->Draw();

  TLine* plusline = new TLine(dtmin,3,dtmax,3);
  plusline->SetLineColor(kBlue);
  plusline->SetLineWidth(2);
  plusline->SetLineStyle(kDashed);
  plusline->Draw();

  TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
  minusline->SetLineColor(kBlue);
  minusline->SetLineWidth(2);
  minusline->SetLineStyle(kDashed);
  minusline->Draw();

  out.str("");
  out << "pics/lifetime_m" << m_mode;
  if(no_bkg)          out << "_nobkg";
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  else if(m_svd == 1) out << "_svd1";
  if(!NFreePar)       out << "_def";
  string rootname = out.str() + string(".root");
  c1->Print(rootname.c_str());
  string epsname = out.str() + string(".eps");
  c1->Print(epsname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());
  return;
}

void DrawSBTest(const int NFreePar,const bool fitted = false){
  const int NTot = m_bkg_tree->GetEntries();
  const double ddt  = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];

  cout << "m_f_bkg: " << m_f_bkg << endl;

  for(int j=0; j<NDots; j++){
    dt_arr[j] = dtmin+(j+0.5)*ddt;
  }
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }
  double Norm = 0;
  double pdf_array[NDots];
  for(int i=0; i<NDots; i++) pdf_array[i] = 0;
  const int nBD = NDots/NBins;

  stringstream out;
  out.str("");
  string label("");
  if(m_svd == 2)      label += string("SVD2");
  else if(m_svd == 1) label += string("SVD1");
  label += get_label(m_mode);
  if(sgl_asc) label += string(" (sgl)");
  if(mlt_asc) label += string(" (mlt)");

  string hname = string("SB test") + label;
  TH1I* DH = new TH1I("DH",hname.c_str(),NBins,dtmin,dtmax);
  DH->SetMarkerStyle(20);
  DH->SetMarkerColor(kBlack);
  DH->SetLineColor(kBlack);
  DH->SetMarkerSize(1.1);

  for(int j=0; j<NTot; j++){
    GetEvent(j,0,1);
    DH->Fill(m_dt);
    for(int i=0; i<NDots; i++){
      double pdfval = Pdf(dt_arr[i],BACKGROUND);
      if(pdfval<0) cout << "pdf < 0: " << pdfval << endl;
      if(!std::isnan(pdfval)){ pdf_array[i] += pdfval;}
      else{
        cout << "pdf is nan: " << dt_arr[i] << endl;
        cout << "  " << m_ntrk_rec << " " << m_sz_rec << " " << m_chisq_rec << " " << m_ndf_rec << endl;
        cout << "  " << m_ntrk_asc << " " << m_sz_asc << " " << m_chisq_asc << " " << m_ndf_asc << endl;
      }
    }
  }

  double pull_array[NBins];
  double chisq = 0;
  int NDF = 0;
  for(int i=0; i<NDots; i++){
    pdf_array[i] *= ddtb;
    Norm += pdf_array[i];
  }
  for(int i=0; i<NDots; i++){
    pdf_array[i] /= Norm/NTot*nBD;
    if(!(i%nBD)){
      const int bin = i/nBD;
      const int bin_content = DH->GetBinContent(bin+1);
      if(bin_content){
        pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
        if(bin_content>5){
          chisq += pull_array[bin]*pull_array[bin];
          NDF++;
        }
      }
    }
  }

//  for(int i=0; i<NDots; i++){
//    pdf_array[i] /= Norm/NTot*nBD;
//    cout << pdf_array[i] << endl;
//  }
  cout << "Norm = " << Norm/nBD << ", Nevents = " << NTot << endl;
  chisq /= NDF;//NBins;
  cout << "chi2/n.d.f. = " << chisq << endl;

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetLineWidth(2);
  GR->SetMarkerColor(kBlue);
  GR->SetLineColor(kBlue);

  TGraphErrors* GRP = new TGraphErrors(NBins,dt_barr,pull_array,dt_barr_err,dt_pull_err);
  GRP->SetMarkerStyle(20);
  GRP->SetMarkerSize(1.2);
  GRP->SetMarkerColor(kBlack);
  GRP->SetLineColor(kBlack);

  TCanvas* c1 = new TCanvas("c1","c1",600,800);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.19);
  pad1->Draw(); pad2->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  DH->GetXaxis()->SetTitle("#Deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.06);
  DH->GetXaxis()->SetTitleOffset(0.75);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->Draw("e1");

  TPaveText *pt = new TPaveText(0.60,0.53,0.98,0.6,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt->AddText(out.str().c_str());
  pt->Draw();

  GR->Draw("same");

  pad2->cd();
  GRP->GetYaxis()->SetRangeUser(-5,5);
  GRP->GetXaxis()->SetRangeUser(dtmin,dtmax);
  GRP->Draw("AP");
  TLine* zeroline = new TLine(dtmin,0,dtmax,0);
  zeroline->SetLineColor(kBlue);
  zeroline->SetLineWidth(2);
  zeroline->Draw();

  TLine* plusline = new TLine(dtmin,3,dtmax,3);
  plusline->SetLineColor(kBlue);
  plusline->SetLineWidth(1);
  plusline->SetLineStyle(kDashed);
  plusline->Draw();

  TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
  minusline->SetLineColor(kBlue);
  minusline->SetLineWidth(1);
  minusline->SetLineStyle(kDashed);
  minusline->Draw();

  out.str("");
  out << "pics/sb_test_m" << m_mode;
  if(m_type_flag == 1)      out << "_cont";
  else if(m_type_flag == 2) out << "_BB";
  if(m_gg) out << "_gg";
  if(m_ppp) out << "_ppp";
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  else if(m_svd == 1) out << "_svd1";
  if(!NFreePar)       out << "_def";
  if(fitted)          out << "_fitted";
  string rootname =   out.str() + string(".root");
  c1->Print(rootname.c_str());

  string epsname = out.str() + string(".eps");
  c1->Print(epsname.c_str());

  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  return;
}

void DrawSideband(const int NFreePar){
  cout << "DrawSideband..." << endl;
//  int NTot;//,NSig;
//  if(!m_mode)         { NTot = m_NGoodTot;}// NSig = NTot;}
//  else if(full_ds_fit){ NTot = m_NBkgTot;}//  NSig = m_NSigTot;}
//  else{                 NTot = m_NBkg;}//     NSig = m_NSig;}
//  if(no_bkg) NTot = 0;//NSig;
//  const int NTot = m_good_sb_tree->GetEntries();
  const int NTotSB  = m_NSBTot;
  const double ddt  = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];

  for(int j=0; j<NDots; j++){
    dt_arr[j] = dtmin+(j+0.5)*ddt;
  }
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }

  double pdf_array[NDots];
  stringstream out;
  out.str("");
  string label("");
  if(m_svd == 2)      label += string("SVD2");
  else if(m_svd == 1) label += string("SVD1");
  label += get_label(m_mode);
  if(sgl_asc) label += string(" (sgl)");
  if(mlt_asc) label += string(" (mlt)");

  cout << label << endl;

  const int nBD = NDots/NBins;

  gStyle->SetStatW(0.4); gStyle->SetStatH(0.3);
  string hname = string("Sideband fit") + label;
  TH1I* DH = new TH1I("DH",hname.c_str(),NBins,dtmin,dtmax);
  DH->SetMarkerStyle(20);
  DH->SetMarkerColor(kBlack);
  DH->SetLineColor(kBlack);
  DH->SetMarkerSize(1.1);

  for(int i=0; i<NDots; i++) pdf_array[i] = 0;
  for(int j=0; j<NTotSB; j++){
    GetEvent(j,0,2);
    DH->Fill(m_dt);
    for(int i=0; i<NDots; i++){
      double pdfval = Pdf(dt_arr[i],BACKGROUND);
      if(!std::isnan(pdfval)){ pdf_array[i] += pdfval;}
      else{
        cout << "pdf is nan: " << dt_arr[i] << endl;
        cout << "  " << m_ntrk_rec << " " << m_sz_rec << " " << m_chisq_rec << " " << m_ndf_rec << endl;
        cout << "  " << m_ntrk_asc << " " << m_sz_asc << " " << m_chisq_asc << " " << m_ndf_asc << endl;
      }
    }
  }
  double pull_array[NBins];
  double chisq = 0;
  int NDF = 0;
  double Norm = 0;
  for(int i=0; i<NDots; i++){
    pdf_array[i] *= ddtb;
    Norm += pdf_array[i];
  }
  for(int i=0; i<NDots; i++){
    pdf_array[i] /= Norm/NTotSB*nBD;
    if(!(i%nBD)){
      const int bin = i/nBD;
      const int bin_content = DH->GetBinContent(bin+1);
      if(bin_content){
        pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
        if(bin_content>5){ 
          chisq += pull_array[bin]*pull_array[bin];
          NDF++;
        }
      }
    }
  }
  cout << "Norm = " << Norm/nBD << ", Nevents = " << NTotSB << endl;
  chisq /= NDF;
  cout << "chi2/n.d.f. = " << chisq << endl;

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetLineWidth(2);
  GR->SetMarkerColor(kBlue);
  GR->SetLineColor(kBlue);

  TGraphErrors* GRP = new TGraphErrors(NBins,dt_barr,pull_array,dt_barr_err,dt_pull_err);
  GRP->SetMarkerStyle(20);
  GRP->SetMarkerSize(1.2);
  GRP->SetMarkerColor(kBlack);
  GRP->SetLineColor(kBlack);

  TCanvas* c1 = new TCanvas("c1","c1",600,800);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.19);
  pad1->Draw(); pad2->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  DH->GetXaxis()->SetTitle("#Deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.06);
  DH->GetXaxis()->SetTitleOffset(0.75);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->Draw("e1");

  TPaveText *pt = new TPaveText(0.60,0.53,0.98,0.6,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt->AddText(out.str().c_str());
  pt->Draw();

  GR->Draw("same");

  pad2->cd();
  GRP->GetYaxis()->SetRangeUser(-5,5);
  GRP->GetXaxis()->SetRangeUser(dtmin,dtmax);
  GRP->Draw("AP");
  TLine* zeroline = new TLine(dtmin,0,dtmax,0);
  zeroline->SetLineColor(kBlue);
  zeroline->SetLineWidth(2);
  zeroline->Draw();

  TLine* plusline = new TLine(dtmin,3,dtmax,3);
  plusline->SetLineColor(kBlue);
  plusline->SetLineWidth(1);
  plusline->SetLineStyle(kDashed);
  plusline->Draw();

  TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
  minusline->SetLineColor(kBlue);
  minusline->SetLineWidth(1);
  minusline->SetLineStyle(kDashed);
  minusline->Draw();

  out.str("");
  out << "pics/sideband_m" << m_mode;
  if(m_type_flag == 1)      out << "_cont";
  else if(m_type_flag == 2) out << "_BB";
  if(m_gg) out << "_gg";
  if(m_ppp) out << "_ppp";
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  else if(m_svd == 1) out << "_svd1";
  if(!NFreePar)       out << "_def";
  string rootname =   out.str() + string(".root");
  c1->Print(rootname.c_str());

  string epsname = out.str() + string(".eps");
  c1->Print(epsname.c_str());

  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  if(m_mode) DrawSBTest(NFreePar);
  return;
}

void Draw_All(const int NFreePar,const int SetNum = 0){
  stringstream out;
  int NTot,NSig;
  if(!m_mode)         { NTot = m_NGoodTot;          NSig = NTot;}
  else if(full_ds_fit){ NTot = m_NSigTot+m_NBkgTot; NSig = m_NSigTot;}
  else{                 NTot = m_NSig+m_NBkg;       NSig = m_NSig;}
  if(no_bkg) NTot = NSig;
  if(m_ebeb) NSig = NTot;
  const double ddt  = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  const int nBD = NDots/NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];
  for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }

  int Nev[2][16];
  TH1I* dh[2][16];
  double pdf_arr[2][16][NDots];
  double norm[2][16];

  int NEveCounter = 0;

  // Initialization
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      Nev[k][i]  = 0;
      norm[k][i] = 0;
      out.str("");
      out << "dh" << bin(i) << "_" << k;
      dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dtmin,dtmax);
      dh[k][i]->SetMarkerStyle(21);
      dh[k][i]->SetMarkerSize(1.1);
      !k ? dh[k][i]->SetMarkerColor(kBlue) : dh[k][i]->SetMarkerColor(kRed);
      for(int j=0; j<NDots; j++){
        pdf_arr[k][i][j] = 0;
      }
    }
  }
  // Filling
  int flv_index, bin_index;
  for(int l=0; l<NTot; l++){
    l<NSig ? GetEvent(l,SetNum,0) : GetEvent(l-NSig,SetNum,1);
    flv_index = flv_ind(m_flv);
    bin_index = bin_ind(m_bin);
    dh[flv_index][bin_index]->Fill(m_dt);
    Nev[flv_index][bin_index]++;
    NEveCounter++;
    for(int j=0; j<NDots; j++){
      double& dt = dt_arr[j];
      double pdfval = Pdf(dt);
      if(!std::isnan(pdfval)){ pdf_arr[flv_index][bin_index][j] += pdfval;}
      else{
        cout << "pdf is nan: dt = " << dt << ", flv: " << m_flv << ", bin: " << m_bin << endl;
        cout << "  m_ntrk_rec: " << m_ntrk_rec << ", m_sz_rec: " << m_sz_rec << ", m_chisq_rec: " << m_chisq_rec << ", m_ndf_rec: " << m_ndf_rec << endl;
        cout << "  m_ntrk_asc: " << m_ntrk_asc << ", m_sz_asc: " << m_sz_asc << ", m_chisq_asc: " << m_chisq_asc << ", m_ndf_asc: " << m_ndf_asc << endl;
      }
    }
  }
  cout << "NEveCounter = " << NEveCounter << endl;

  double pull_array[2][16][NBins];
  double chisq[2][16];

  cout << "Set norm" << endl;
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      chisq[k][i] = 0;
      for(int j=0; j<NDots; j++){
        pdf_arr[k][i][j] *= ddtb;
        norm[k][i] += pdf_arr[k][i][j];
        if(!(j%nBD)){
          const int bin = j/nBD;
          const int bin_content = dh[k][i]->GetBinContent(bin+1);
          if(bin_content){
            pull_array[k][i][bin] = (bin_content-pdf_arr[k][i][j])/sqrt(bin_content);
            chisq[k][i] += pull_array[k][i][bin]*pull_array[k][i][bin];
          }
        }
      }
    }
  }

  TGraph* gr[2][16];
  TGraphErrors* grp[2][16];
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      gr[k][i] = new TGraph(NDots,dt_arr,pdf_arr[k][i]);
      gr[k][i]->SetMarkerStyle(kDot);
      gr[k][i]->SetMarkerSize(1.5);
      gr[k][i]->SetLineWidth(2);
      k == 0 ? gr[k][i]->SetMarkerColor(kBlue) : gr[k][i]->SetMarkerColor(kRed);

      grp[k][i] = new TGraphErrors(NBins,dt_barr,pull_array[k][i],dt_barr_err,dt_pull_err);
      grp[k][i]->SetMarkerStyle(21);
      grp[k][i]->SetMarkerSize(1.2);
      k == 0 ? grp[k][i]->SetMarkerColor(kBlue) : grp[k][i]->SetMarkerColor(kRed);
      grp[k][i]->SetLineColor(0);
    }
  }

  cout << "Draw" << endl;
  TCanvas* c1 = new TCanvas("c1","c1",2600,1000);
  TPad *pad[8];
  TPad *padb[8];
  TPad *padp[8];
  TPad *padpb[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad" << i;
    pad[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.61,0.125*(i+1)-0.01,0.99);
    pad[i]->Draw();
    out.str("");
    out << "padb" << i;
    padb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.11,0.125*(i+1)-0.01,0.49);
    padb[i]->Draw();

    out.str("");
    out << "padp" << i;
    padp[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.59);
    padp[i]->Draw();
    out.str("");
    out << "padpb" << i;
    padpb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.09);
    padpb[i]->Draw();
  }
  int k = 0;
  for(int i=0; i<8; i++){
    pad[i]->cd();
    pad[i]->SetGrid();
    pad[i]->SetLogy();
    dh[k][i]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i]->Draw("e");
    pad[i]->Update();
    TPaveStats *st = (TPaveStats*)dh[k][i]->FindObject("stats");
    st->SetOptStat(1110);
    st->PaintPave(0.55,0.75,0.99,0.9);
    pad[i]->Update();
    gr[k][i]->Draw("same");

    padb[i]->cd();
    padb[i]->SetGrid();
    padb[i]->SetLogy();
    dh[k][i+8]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i+8]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i+8]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i+8]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i+8]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i+8]->Draw("e");
    padb[i]->Update();
    st = (TPaveStats*)dh[k][i+8]->FindObject("stats");
    st->SetOptStat(1110);
    st->PaintPave(0.6,0.7,0.99,0.9);
    padb[i]->Update();
    gr[k][i+8]->Draw("same");

    padp[i]->cd();
    grp[k][i]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i]->Draw();
    TLine* zeroline = new TLine(dtmin,0,dtmax,0);
    zeroline->SetLineWidth(2);
    zeroline->Draw("AP");

    TLine* plusline = new TLine(dtmin,3,dtmax,3);
    plusline->SetLineWidth(1);
    plusline->SetLineStyle(kDashed);
    plusline->Draw();

    TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
    minusline->SetLineWidth(1);
    minusline->SetLineStyle(kDashed);
    minusline->Draw();

    padpb[i]->cd();
    grp[k][i+8]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i+8]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i+8]->Draw();
    TLine* zerolineb = new TLine(dtmin,0,dtmax,0);
    zerolineb->SetLineWidth(2);
    zerolineb->Draw("AP");

    TLine* pluslineb = new TLine(dtmin,3,dtmax,3);
    pluslineb->SetLineWidth(1);
    pluslineb->SetLineStyle(kDashed);
    pluslineb->Draw();

    TLine* minuslineb = new TLine(dtmin,-3,dtmax,-3);
    minuslineb->SetLineWidth(1);
    minuslineb->SetLineStyle(kDashed);
    minuslineb->Draw();
  }

  out.str("");
  out << "pics/full_cpfit_m" << m_mode;
  if(mlt_asc)         out << "_mlt_";
  if(sgl_asc)         out << "_sgl_";
  if(m_svd == 2)      out << "svd2";
  else if(m_svd == 1) out << "svd1";
  if(!NFreePar)       out << "_def";
  if(no_bkg)          out << "_nobkg";
  if(perftag)         out << "_perftag";
  if(!full_ds_fit){
    out << "_nsig_" << m_NSig;
    out << "_fsig_" << m_fSig;
  }
  string rootname =   out.str() + string(".root");

  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  TCanvas* c2 = new TCanvas("c2","c2",2600,1000);
  c2->cd();
  TPad *pad1[8];
  TPad *padb1[8];
  TPad *padp1[8];
  TPad *padpb1[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad1" << i;
    pad1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.61,0.125*(i+1)-0.01,0.99);
    pad1[i]->Draw();
    out.str("");
    out << "padb1" << i;
    padb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.11,0.125*(i+1)-0.01,0.49);
    padb1[i]->Draw();

    out.str("");
    out << "padp1" << i;
    padp1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.59);
    padp1[i]->Draw();
    out.str("");
    out << "padpb1" << i;
    padpb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.09);
    padpb1[i]->Draw();
  }
  k = 1;
  for(int i=0; i<8; i++){
    pad1[i]->cd();
//    gStyle->SetOptStat("RMe");
    pad1[i]->SetGrid();
    pad1[i]->SetLogy();
    dh[k][i]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i]->Draw("e");
    pad1[i]->Update();
    TPaveStats *st = (TPaveStats*)dh[k][i]->FindObject("stats");
    st->SetOptStat(2210);
    st->PaintPave(0.60,0.65,0.99,0.9);
    pad1[i]->Update();
    gr[k][i]->Draw("same");

    padb1[i]->cd();
    padb1[i]->SetGrid();
    padb1[i]->SetLogy();
    dh[k][i+8]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i+8]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i+8]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i+8]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i+8]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i+8]->Draw("e");
    padb1[i]->Update();
    st = (TPaveStats*)dh[k][i+8]->FindObject("stats");
    st->SetOptStat(2210);
    st->PaintPave(0.65,0.65,0.99,0.9);
    padb1[i]->Update();
    gr[k][i+8]->Draw("same");

    padp1[i]->cd();
    grp[k][i]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i]->Draw();
    TLine* zeroline1 = new TLine(dtmin,0,dtmax,0);
    zeroline1->SetLineWidth(2);
    zeroline1->Draw("AP");

    TLine* plusline = new TLine(dtmin,3,dtmax,3);
    plusline->SetLineWidth(1);
    plusline->SetLineStyle(kDashed);
    plusline->Draw();

    TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
    minusline->SetLineWidth(1);
    minusline->SetLineStyle(kDashed);
    minusline->Draw();

    padpb1[i]->cd();
    grp[k][i+8]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i+8]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i+8]->Draw();
    TLine* zerolineb = new TLine(dtmin,0,dtmax,0);
    zerolineb->SetLineWidth(2);
    zerolineb->Draw("AP");

    TLine* pluslineb = new TLine(dtmin,3,dtmax,3);
    pluslineb->SetLineWidth(1);
    pluslineb->SetLineStyle(kDashed);
    pluslineb->Draw();

    TLine* minuslineb = new TLine(dtmin,-3,dtmax,-3);
    minuslineb->SetLineWidth(1);
    minuslineb->SetLineStyle(kDashed);
    minuslineb->Draw();
  }

  out.str("");
  out << "pics/full_cpfit_m" << m_mode;
  if(mlt_asc)         out << "_mlt_";
  if(sgl_asc)         out << "_sgl_";
  if(m_svd == 2)      out << "svd2";
  else if(m_svd == 1) out << "svd1";
  if(!NFreePar)       out << "_def";
  out << "_b0b";
  rootname = out.str() + string(".root");

  c2->Print(rootname.c_str());
  pngname = out.str() + string(".png");
  c2->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  cout << "NEveCounter = " << NEveCounter << endl;
  return;
}

void draw_single_fit(const int SetNum=0){
  TH1I dh("dh","dh",NBins,dtmin,dtmax);
  int NTot,NSig;
  if(!m_mode)         { NTot = m_NGoodTot;          NSig = NTot;}
  else if(full_ds_fit){ NTot = m_NSigTot+m_NBkgTot; NSig = m_NSigTot;}
  else{                 NTot = m_NSig+m_NBkg;       NSig = m_NSig;}
  if(no_bkg) NTot = NSig;
  if(m_ebeb) NSig = NTot;
  for(int i=0; i<NTot; i++){
    i<NSig ? GetEvent(i,SetNum,0) : GetEvent(i-NSig,SetNum,0);
    dh.Fill(m_dt);
    dh.SetMarkerStyle(20);
    dh.SetMarkerColor(kBlue);
    dh.SetMarkerSize(1.2);
  }
  double pdf_arr[NDots],dt_arr[NDots];
  double norm = 0;

  const double ddt = (dtmax-dtmin)/(double)NDots;
  const double Nddt = ddt*NDots/NBins;

  double pdf = 0;
  for(int i=0; i<NDots; i++){
    pdf_arr[i] = 0;
    dt_arr[i] = dtmin+(i+0.5)*ddt;
    const double dt = dt_arr[i];//cm2ps;
    for(int j=0; j<NTot; j++){
      j<NSig ? GetEvent(j,SetNum,0) : GetEvent(j-NSig,SetNum,0);
      pdf = Pdf(dt,0);
      if(!std::isnan(pdf)) pdf_arr[i] += pdf;
    }
    norm += pdf_arr[i];
  }
  for(int i=0; i<NDots; i++) pdf_arr[i] *= NTot/norm*Nddt;
  norm *= ddt/NTot;

  cout << "norm = " << norm << endl;

  TGraph* gr = new TGraph(NDots,dt_arr,pdf_arr);
  gr->SetLineWidth(2.);

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  c1->SetLogy();
  dh.Draw("e");
  gr->Draw("same");

  c1->Print("fit.png");
  c1->Print("fit.root");
  system("display fit.png &");
  return;
}

void DrawToyLifetime(vector<double> val_vec,vector<double> err_vec){
  cout << "Drawing toy lifetime" << endl;
  stringstream out;
  out.str("");
  out << "pics/toy_tau.root";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname1 = out.str();
  out.str("");
  out << "pics/toy_tau_pull";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname2 = out.str();
//  const string label = get_label(m_mode);

  TCanvas* c1  = new TCanvas("c1","c1",400,400);
  TCanvas* c2  = new TCanvas("c2","c2",400,400);
  const int NFits = val_vec.size();

  string hpname = string("#tau pull");// + label;
  TH1D *hp_tau = new TH1D("pullhist",hpname.c_str(),100,-5.,5.);
  hp_tau->SetMarkerStyle(21);
  hp_tau->SetMarkerColor(kBlue);
  hp_tau->SetMarkerSize(1.1);
  hp_tau->GetXaxis()->SetTitle("#tau pull");
  hp_tau->GetXaxis()->SetTitleSize(0.06);
  hp_tau->GetXaxis()->SetLabelSize(0.05);
  hp_tau->GetXaxis()->SetTitleOffset(0.75);
  hp_tau->GetYaxis()->SetLabelSize(0.05);
  hp_tau->GetYaxis()->SetTitleSize(0.05);

  hpname = string("#tau (ps)");// + label;
  TH1D *h_tau  = new TH1D("tauhist",hpname.c_str(),100,1.,2.);
  h_tau->SetMarkerStyle(21);
  h_tau->SetMarkerColor(kBlue);
  h_tau->SetMarkerSize(1.1);
  h_tau->GetXaxis()->SetTitle("#tau (ps)");
  h_tau->GetXaxis()->SetTitleSize(0.06);
  h_tau->GetXaxis()->SetLabelSize(0.05);
  h_tau->GetXaxis()->SetTitleOffset(0.75);
  h_tau->GetYaxis()->SetLabelSize(0.05);

  for(int i=0; i<NFits; i++){
    h_tau->Fill(val_vec[i]);
    hp_tau->Fill((val_vec[i] - m_btau)/err_vec[i]);
  }

  c1->cd();
  c1->SetGrid();
  h_tau->Draw();
  TLine* tauline = new TLine(m_btau,0,m_btau,0.07*h_tau->GetEntries());
  tauline->SetLineColor(kRed);
  tauline->SetLineWidth(2);
  tauline->SetLineStyle(9);
  tauline->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)h_tau->FindObject("stats");
  st->SetOptStat(220002220);
  st->PaintPave(0.55,0.75,0.99,0.9);
  c1->Update();

  c2->cd();
  c2->SetGrid();
  hp_tau->Draw();
  hp_tau->Fit("gaus");
  c2->Update();

  out.str("");
  out << fname1 << ".root";
  c1->Print(out.str().c_str());

  out.str("");
  out << fname1 << ".eps";
  c1->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname1 << ".eps &";
  system(out.str().c_str());

  out.str("");
  out << fname2 << ".root";
  c2->Print(out.str().c_str());

  out.str("");
  out << fname2 << ".eps";
  c2->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname2 << ".eps &";
  system(out.str().c_str());
  return;
}

void DrawToyCPV(vector<double> *val_vec,vector<double> *err_vec){
  cout << "Drawing toy CPV fit" << endl;
  stringstream out;
  out.str("");
  out << "pics/toy_sin";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname_sin1 = out.str();
  out.str("");
  out << "pics/toy_sin_pull";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname_sin2 = out.str();

  out.str("");
  out << "pics/toy_cos";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname_cos1 = out.str();
  out.str("");
  out << "pics/toy_cos_pull";// << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
//  out << "_n" << m_toysize;
//  if(perfbin) out << "_pb";
//  if(perftag) out << "_pt";
  const string fname_cos2 = out.str();


  const string label = get_label(m_mode);

  TCanvas* c1_sin  = new TCanvas("c1_sin","c1_sin",400,400);
  TCanvas* c2_sin  = new TCanvas("c2_sin","c2_sin",400,400);
  TCanvas* c1_cos  = new TCanvas("c1_cos","c1_cos",400,400);
  TCanvas* c2_cos  = new TCanvas("c2_cos","c2_cos",400,400);
  const int NFits = val_vec[0].size();

  string hpname = string("sin(2#phi_{1}) pull");// + label;
  TH1D *hp_sin = new TH1D("pullsinhist",hpname.c_str(),100,-5.,5.);
  hp_sin->SetMarkerStyle(21);
  hp_sin->SetMarkerColor(kBlue);
  hp_sin->SetMarkerSize(1.1);
  hp_sin->GetXaxis()->SetTitle("sin(2#phi_{1}) pull");
  hp_sin->GetXaxis()->SetTitleSize(0.06);
  hp_sin->GetXaxis()->SetLabelSize(0.05);
  hp_sin->GetXaxis()->SetTitleOffset(0.75);
  hp_sin->GetYaxis()->SetLabelSize(0.05);
  hp_sin->GetYaxis()->SetTitleSize(0.05);

  hpname = string("sin(2#phi_{1})");// + label;
  TH1D *h_sin  = new TH1D("sinhist",hpname.c_str(),100,-2.,2.);
  h_sin->SetMarkerStyle(21);
  h_sin->SetMarkerColor(kBlue);
  h_sin->SetMarkerSize(1.1);
  h_sin->GetXaxis()->SetTitle("sin(2#phi_{1})");
  h_sin->GetXaxis()->SetTitleSize(0.06);
  h_sin->GetXaxis()->SetLabelSize(0.05);
  h_sin->GetXaxis()->SetTitleOffset(0.75);
  h_sin->GetYaxis()->SetLabelSize(0.05);

  hpname = string("cos(2#phi_{1}) pull");// + label;
  TH1D* hp_cos = new TH1D("pullcoshist",hpname.c_str(),100,-5.,5.);
  hp_cos->SetMarkerStyle(21);
  hp_cos->SetMarkerColor(kBlue);
  hp_cos->SetMarkerSize(1.1);
  hp_cos->GetXaxis()->SetTitle("cos(2#phi_{1}) pull");
  hp_cos->GetXaxis()->SetTitleSize(0.06);
  hp_cos->GetXaxis()->SetLabelSize(0.05);
  hp_cos->GetXaxis()->SetTitleOffset(0.75);
  hp_cos->GetYaxis()->SetLabelSize(0.05);
  hp_cos->GetYaxis()->SetTitleSize(0.05);

  hpname = string("cos(2#phi_{1})");// + label;
  TH1D *h_cos  = new TH1D("coshist",hpname.c_str(),100,-2.,2.);
  h_cos->SetMarkerStyle(21);
  h_cos->SetMarkerColor(kBlue);
  h_cos->SetMarkerSize(1.1);
  h_cos->GetXaxis()->SetTitle("cos(2#phi_{1})");
  h_cos->GetXaxis()->SetTitleSize(0.06);
  h_cos->GetXaxis()->SetLabelSize(0.05);
  h_cos->GetXaxis()->SetTitleOffset(0.75);
  h_cos->GetYaxis()->SetLabelSize(0.05);

  for(int i=0; i<NFits; i++){
    h_sin->Fill(val_vec[0][i]);
    hp_sin->Fill((val_vec[0][i] - m_sin2beta)/err_vec[0][i]);
    h_cos->Fill(val_vec[1][i]);
    hp_cos->Fill((val_vec[1][i] - m_cos2beta)/err_vec[1][i]);
  }

  c1_sin->cd();
  c1_sin->SetGrid();
  h_sin->Draw();
  TLine* sinline = new TLine(m_sin2beta,0,m_sin2beta,0.07*h_sin->GetEntries());
  sinline->SetLineColor(kRed);
  sinline->SetLineWidth(2);
  sinline->SetLineStyle(9);
  sinline->Draw();
  c1_sin->Update();
  TPaveStats *st = (TPaveStats*)h_sin->FindObject("stats");
  st->SetOptStat(220002220);
  st->PaintPave(0.55,0.75,0.99,0.9);
  c1_sin->Update();

  c2_sin->cd();
  c2_sin->SetGrid();
  hp_sin->Draw();
  hp_sin->Fit("gaus");
  c2_sin->Update();

  out.str("");
  out << fname_sin1 << ".root";
  c1_sin->Print(out.str().c_str());

  out.str("");
  out << fname_sin1 << ".eps";
  c1_sin->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname_sin1 << ".eps &";
  system(out.str().c_str());

  out.str("");
  out << fname_sin2 << ".root";
  c2_sin->Print(out.str().c_str());

  out.str("");
  out << fname_sin2 << ".eps";
  c2_sin->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname_sin2 << ".eps &";
  system(out.str().c_str());

  c1_cos->cd();
  c1_cos->SetGrid();
  h_cos->Draw();
  TLine* cosline = new TLine(m_cos2beta,0,m_cos2beta,0.07*h_cos->GetEntries());
  cosline->SetLineColor(kRed);
  cosline->SetLineWidth(2);
  cosline->SetLineStyle(9);
  cosline->Draw();
  c1_cos->Update();
  st = (TPaveStats*)h_cos->FindObject("stats");
  st->SetOptStat(220002220);
  st->PaintPave(0.55,0.75,0.99,0.9);
  c1_cos->Update();

  c2_cos->cd();
  c2_cos->SetGrid();
  hp_cos->Draw();
  hp_cos->Fit("gaus");
  c2_cos->Update();

  out.str("");
  out << fname_cos1 << ".root";
  c1_cos->Print(out.str().c_str());

  out.str("");
  out << fname_cos1 << ".eps";
  c1_cos->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname_cos1 << ".eps &";
  system(out.str().c_str());

  out.str("");
  out << fname_cos2 << ".root";
  c2_cos->Print(out.str().c_str());

  out.str("");
  out << fname_cos2 << ".eps";
  c2_cos->Print(out.str().c_str());
  out.str("");
  out << "evince " << fname_cos2 << ".eps &";
  system(out.str().c_str());
  return;
}

void DrawFullSampling(vector<double>* val_vec,vector<double>* err_vec){
  cout << "DrawFullSampling" << endl;

  stringstream out;
  out.str("");
  out << "hists/";
  if(no_interf) out << "tau";
  else          out << "sin_cos";
  out << "_sampling_m" << m_mode;
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  else if(m_svd == 1) out << "_svd1";
  if(no_bkg)          out << "_nobkg";
  if(perftag)         out << "_perftag";
  if(m_ebeb)          out << "_ebeb";
  if(!full_ds_fit){
    out << "_nsig_" << m_NSig;
    out << "_fsig_" << m_fSig;
  }
  out << "_v3";
  string filename = out.str() + string(".root");
  TFile f(filename.c_str(),"recreate");

  double mean_sin, mean_cos, mean_tau;
  double err_sin, err_cos, err_tau;
  double merr_sin, merr_cos, merr_tau;
  for(int i=0; i<16; i++){
    if(no_interf) GetMoments(val_vec[0],mean_tau,err_tau,merr_tau);
    else{
      GetMoments(val_vec[1],mean_sin,err_sin,merr_sin);
      GetMoments(val_vec[2],mean_cos,err_cos,merr_cos);
    }
  }

  string label("");
  label += get_label(m_mode);
  if(m_svd == 2) label += string(", SVD2");
  if(m_svd == 1) label += string(", SVD1");
  if(no_bkg)     label += string(", no Bkg");
  if(perftag)    label += string(", perf. tag");
  if(m_fitbin){
    out.str("");
    out << m_fitbin;
    label += string(", bin ") + out.str();
  }

  TCanvas* c1;
  TCanvas* c2;
  TH1D *h_sin, *hp_sin, *h_cos, *hp_cos, *h_tau, *hp_tau;
  const int NFits = val_vec[0].size();
  if(no_interf){
    string hpname = string("#tau pull") + label;
    hp_tau = new TH1D("#tau pull",hpname.c_str(),100,-10.,10.);
    hp_tau->SetMarkerStyle(21);
    hp_tau->SetMarkerColor(kBlue);
    hp_tau->SetMarkerSize(1.1);
    hp_tau->GetXaxis()->SetTitle("#tau pull");
    hp_tau->GetXaxis()->SetTitleSize(0.06);
    hp_tau->GetXaxis()->SetLabelSize(0.05);
    hp_tau->GetXaxis()->SetTitleOffset(0.75);
    hp_tau->GetYaxis()->SetLabelSize(0.05);
    hp_tau->GetYaxis()->SetTitleSize(0.05);
    hpname = string("#tau (ps)") + label;
    double taumin = 0.5;//m_NSig > 300 ? 1 : 0.5;
    double taumax = 3;//m_NSig > 300 ? 2 : 3;
    h_tau = new TH1D("#tau",hpname.c_str(),100,taumin,taumax);
    h_tau->SetMarkerStyle(21);
    h_tau->SetMarkerColor(kBlue);
    h_tau->SetMarkerSize(1.1);
    h_tau->GetXaxis()->SetTitle("#tau (ps)");
    h_tau->GetXaxis()->SetTitleSize(0.06);
    h_tau->GetXaxis()->SetLabelSize(0.05);
    h_tau->GetXaxis()->SetTitleOffset(0.75);
    h_tau->GetYaxis()->SetLabelSize(0.05);
    for(int i=0; i<NFits; i++){
      h_tau->Fill(val_vec[0][i]);
      hp_tau->Fill((val_vec[0][i] - m_btau)/err_vec[0][i]);
    }
    c1 = new TCanvas("c1","c1",400,400);
    c1->cd();
    c1->SetGrid();
    h_tau->Draw();
    TLine* tauline = new TLine(m_btau,0,m_btau,0.07*h_tau->GetEntries());
    tauline->SetLineColor(kRed);
    tauline->SetLineWidth(2);
    tauline->SetLineStyle(9);
    tauline->Draw();
    c1->Update();
    TPaveStats *st = (TPaveStats*)h_tau->FindObject("stats");
    st->SetOptStat(220002220);
    st->PaintPave(0.55,0.75,0.99,0.9);
    c1->Update();
    c2 = new TCanvas("c2","c2",400,400);
    c2->cd();
    c2->SetGrid();
    hp_tau->Draw();
    hp_tau->Fit("gaus");
    c2->Update();

    if(!fix_btau){
      h_tau->Write();
      hp_tau->Write();
    }
  } else{
    string hpname;
    if(!fix_sin2beta){
      const double sin_max = 3.;
      hpname = string("sin(2#varphi_{1})") + label;
      h_sin = new TH1D("sin",hpname.c_str(),100,-sin_max,sin_max);
      h_sin->SetMarkerStyle(21);
      h_sin->SetMarkerColor(kBlue);
      h_sin->SetMarkerSize(1.1);
      h_sin->GetXaxis()->SetTitle("sin(2#varphi_{1})");
      h_sin->GetXaxis()->SetTitleSize(0.06);
      h_sin->GetXaxis()->SetLabelSize(0.05);
      h_sin->GetXaxis()->SetTitleOffset(0.75);
      h_sin->GetYaxis()->SetLabelSize(0.05);
      h_sin->GetYaxis()->SetTitleSize(0.05);
    }
    if(!fix_cos2beta){
      const double cos_max = 3.;
      hpname = string("cos(2#varphi_{1})") + label;
      h_cos = new TH1D("cos",hpname.c_str(),100,-cos_max,cos_max);
      h_cos->SetMarkerStyle(21);
      h_cos->SetMarkerColor(kBlue);
      h_cos->SetMarkerSize(1.1);
      h_cos->GetXaxis()->SetTitle("cos(2#varphi_{1})");
      h_cos->GetXaxis()->SetTitleSize(0.06);
      h_cos->GetXaxis()->SetLabelSize(0.05);
      h_cos->GetXaxis()->SetTitleOffset(0.75);
      h_cos->GetYaxis()->SetLabelSize(0.05);
      h_cos->GetYaxis()->SetTitleSize(0.05);
    }
    if(!fix_sin2beta){
      hpname = string("sin(2#varphi_{1}) pull") + label;
      hp_sin = new TH1D("sin_pull",hpname.c_str(),100,-10.,10.);
      hp_sin->SetMarkerStyle(21);
      hp_sin->SetMarkerColor(kBlue);
      hp_sin->SetMarkerSize(1.1);
      hp_sin->GetXaxis()->SetTitle("sin(2#varphi_{1}) pull");
      hp_sin->GetXaxis()->SetTitleSize(0.06);
      hp_sin->GetXaxis()->SetLabelSize(0.05);
      hp_sin->GetXaxis()->SetTitleOffset(0.75);
      hp_sin->GetYaxis()->SetLabelSize(0.05);
      hp_sin->GetYaxis()->SetTitleSize(0.05);
    }
    if(!fix_cos2beta){
      hpname = string("cos(2#varphi_{1}) pull") + label;
      hp_cos = new TH1D("cos_pull",hpname.c_str(),100,-10.,10.);
      hp_cos->SetMarkerStyle(21);
      hp_cos->SetMarkerColor(kBlue);
      hp_cos->SetMarkerSize(1.1);
      hp_cos->GetXaxis()->SetTitle("cos(2#varphi_{1}) pull");
      hp_cos->GetXaxis()->SetTitleSize(0.06);
      hp_cos->GetXaxis()->SetLabelSize(0.05);
      hp_cos->GetXaxis()->SetTitleOffset(0.75);
      hp_cos->GetYaxis()->SetLabelSize(0.05);
      hp_cos->GetYaxis()->SetTitleSize(0.05);
    }
    for(int i=0; i<NFits; i++){
      if(!fix_sin2beta) h_sin->Fill(val_vec[1][i]);
      if(!fix_cos2beta) h_cos->Fill(val_vec[2][i]);
      if(!fix_sin2beta) hp_sin->Fill((val_vec[1][i] - m_sin2beta)/err_vec[1][i]);
      if(!fix_cos2beta) hp_cos->Fill((val_vec[2][i] - m_cos2beta)/err_vec[2][i]);
    }
    if(!fix_sin2beta && !fix_cos2beta){
      c1 = new TCanvas("c1","c1",800,400);
      c1->Divide(2,1);
      c1->cd(1);
      c1->Pad()->SetGrid();
      h_sin->Draw();
      TLine* sinline = new TLine(m_sin2beta,0,m_sin2beta,0.07*h_sin->GetEntries());
      sinline->SetLineColor(kRed);
      sinline->SetLineWidth(2);
      sinline->SetLineStyle(9);
      sinline->Draw();
      c1->Update();
      TPaveStats *st_sin = (TPaveStats*)h_sin->FindObject("stats");
      st_sin->SetOptStat(220002220);
      st_sin->PaintPave(0.55,0.75,0.99,0.9);
      c1->cd(2);
      c1->Pad()->SetGrid();
      h_cos->Draw();
      TLine* cosline = new TLine(m_cos2beta,0,m_cos2beta,0.07*h_cos->GetEntries());
      cosline->SetLineColor(kRed);
      cosline->SetLineWidth(2);
      cosline->SetLineStyle(9);
      cosline->Draw();
      c1->Update();
      TPaveStats *st_cos = (TPaveStats*)h_cos->FindObject("stats");
      st_cos->SetOptStat(220002220);
      st_cos->PaintPave(0.55,0.75,0.99,0.9);
      c1->Update();
      c2 = new TCanvas("c2","c2",800,400);
      c2->Divide(2,1);
      c2->cd(1);
      c2->Pad()->SetGrid();
      hp_sin->Draw();
      hp_sin->Fit("gaus");
      c2->cd(2);
      c2->Pad()->SetGrid();
      hp_cos->Draw();
      hp_cos->Fit("gaus");
      c2->Update();
    } else if(!fix_cos2beta){
      c1 = new TCanvas("c1","c1",400,400);
      c1->cd();
      c1->SetGrid();
      h_cos->Draw();
      TLine* cosline = new TLine(m_cos2beta,0,m_cos2beta,0.07*h_cos->GetEntries());
      cosline->SetLineColor(kRed);
      cosline->SetLineWidth(2);
      cosline->SetLineStyle(9);
      cosline->Draw();
      c1->Update();
      TPaveStats *st_cos = (TPaveStats*)h_cos->FindObject("stats");
      st_cos->SetOptStat(220002220);
      st_cos->PaintPave(0.55,0.75,0.99,0.9);
      c1->Update();
      c2 = new TCanvas("c2","c2",400,400);
      c2->cd();
      c2->SetGrid();
      hp_cos->Draw();
      hp_cos->Fit("gaus");
      c2->Update();
    } else if(!fix_sin2beta){
      c1 = new TCanvas("c1","c1",400,400);
      c1->cd();
      c1->SetGrid();
      h_sin->Draw();
      TLine* sinline = new TLine(m_sin2beta,0,m_sin2beta,0.07*h_sin->GetEntries());
      sinline->SetLineColor(kRed);
      sinline->SetLineWidth(2);
      sinline->SetLineStyle(9);
      sinline->Draw();
      c1->Update();
      c2 = new TCanvas("c2","c2",400,400);
      c2->cd();
      c2->SetGrid();
      hp_sin->Draw();
      hp_sin->Fit("gaus");
      c2->Update();
    } else{
      cout << "Nothing to draw.." << endl;
      return;
    }

    if(!fix_cos2beta){
      h_cos->Write();
      hp_cos->Write();
    }
    if(!fix_sin2beta){
      h_sin->Write();
      hp_sin->Write();
    }
  }

  cout << "Almost done.." << endl;
  out.str("");
  out << "pics/";
  if(no_interf) out << "tau";
  else          out << "sin_cos";
  if(fix_btau && no_interf)  out << "_fixtau";
  if(fix_sin2beta && !no_interf)  out << "_fixsin";
  if(fix_cos2beta && !no_interf)  out << "_fixcos";
  out << "_sampling_m" << m_mode;
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  else if(m_svd == 1) out << "_svd1";
  if(no_bkg)          out << "_nobkg";
  if(perftag)         out << "_perftag";
  if(m_fitbin)        out << "_bin_" << m_fitbin;
  if(m_ebeb)          out << "_ebeb";
  out << "_v3";

  if(!full_ds_fit){
    out << "_nsig_" << m_NSig;
    out << "_fsig_" << m_fSig;
  }
  string rootname = out.str() + string(".root");

  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  string fn = string("display ") + pngname + string(" &");
  system(fn.c_str());

  out << "_pull";
  rootname = out.str() + string(".root");
  c2->Print(rootname.c_str());
  pngname = out.str() + string(".png");
  c2->Print(pngname.c_str());
  fn = string("display ") + pngname + string(" &");
  system(fn.c_str());

  f.Write();
  f.Close();
  return;
}

void DrawABSampling(vector<double> val_vec[][2],vector<double> err_vec[][2], const double* Aref, const double* Bref){
  double mean[16][2], rms[16][2], merr[16][2];
  for(int i=0; i<16; i++){
    GetMoments(val_vec[i][0],mean[i][0],rms[i][0],merr[i][0]);
    GetMoments(val_vec[i][1],mean[i][1],rms[i][1],merr[i][1]);
  }
  draw_fit_results(mean,merr,Aref,Bref);

  const int NFits = val_vec[0][0].size();
  TH1D* h_A_posi = new TH1D("A_posi","A_posi",100,-10.,10.);
  h_A_posi->SetMarkerStyle(21);
  h_A_posi->SetMarkerColor(kBlue);
  h_A_posi->SetMarkerSize(1.1);
  TH1D* h_A_nega = new TH1D("A_nega","A_nega",100,-10.,10.);
  h_A_nega->SetMarkerStyle(21);
  h_A_nega->SetMarkerColor(kBlue);
  h_A_nega->SetMarkerSize(1.1);
  TH1D* h_B_posi = new TH1D("B_posi","B_posi",100,-10.,10.);
  h_B_posi->SetMarkerStyle(21);
  h_B_posi->SetMarkerColor(kBlue);
  h_B_posi->SetMarkerSize(1.1);
  TH1D* h_B_nega = new TH1D("B_nega","B_nega",100,-10.,10.);
  h_B_nega->SetMarkerStyle(21);
  h_B_nega->SetMarkerColor(kBlue);
  h_B_nega->SetMarkerSize(1.1);
  for(int i=0; i<NFits; i++){
    for(int j=0; j<16; j++){
      if(j<8){
        h_A_nega->Fill((val_vec[j][0][i] - Aref[j])/err_vec[j][0][i]);
        h_B_nega->Fill((val_vec[j][1][i] - Bref[j])/err_vec[j][1][i]);
      } else{
        h_A_posi->Fill((val_vec[j][0][i] - Aref[j])/err_vec[j][0][i]);
        h_B_posi->Fill((val_vec[j][1][i] - Bref[j])/err_vec[j][1][i]);
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  h_A_posi->Draw();
  c1->cd(2);
  h_A_nega->Draw();
  c1->cd(3);
  h_B_posi->Draw();
  c1->cd(4);
  h_B_nega->Draw();
  c1->Update();
  c1->Draw();

  stringstream out;
  out.str("");
  out << "pics/AB_sampling_m" << m_mode;
  if(mlt_asc)         out << "_mlt";
  if(sgl_asc)         out << "_sgl";
  if(m_svd == 2)      out << "_svd2";
  if(m_svd == 1)      out << "_svd1";
  if(m_ebeb)          out << "_ebeb";
//  if(!NFreePar)       out << "_def";
  string rootname = out.str() + string(".root");

  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());
  return;
}

void norm_test(const int Neve){
  cout << "Normalization test:" << endl;
  m_pdf_svd2->SetTauDm(m_btau,m_dm);
  m_pdf_svd2->SetFlvXi(1,xi);

  const int NRanges = 4;
  const int NBinnings = 7;
  double tRanges[NRanges] = {2,10,15,70};
  int Binnings[NBinnings] = {1000,500,250,200,100,50,25};
  for(int k=0; k<NRanges; k++){
    m_pdf_svd2->SetRange(tRanges[k]);
    double max = tRanges[k];
    double min = -max;
    cout << "*** Range (" << min << "," << max << ") ***" << endl;
    double ddt = (max-min)/1000.;
    double dt;
    double NormsRkRdetRnp[NBinnings];
    double NormsRascRnp[NBinnings];
    double valRkRdetRnp, valRascRnp;
    for(int i=0; i<Neve; i++){
      for(int l=0; l<NBinnings; l++){
        NormsRkRdetRnp[l] = 0;
        NormsRascRnp[l]   = 0;
      }
      cout << "--- Event " << i+1 << " ---" <<endl;
      GetEvent(i,0,1);
      m_pdf_svd2->SetAkCk(m_costhBcms,0.5*10.580);
      for(int j=0; j<Binnings[0]; j++){
        dt = min + (j+0.5)*ddt;
        valRkRdetRnp = m_pdf_svd2->Pdf(dt,m_ntrk_rec,m_sz_rec,m_chisq_rec,m_ndf_rec,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc,true,true);
        if(std::isnan(valRkRdetRnp)){
          cout << "RkRdetRnp Pdf is nan: dt = " << dt << endl;
          valRkRdetRnp = 0;
        }
        valRascRnp = m_pdf_svd2->PdfRascRnp(dt,m_ntrk_asc,m_sz_asc,m_chisq_asc,m_ndf_asc);
        if(std::isnan(valRascRnp)){
          cout << "RascRnp Pdf is nan: dt = " << dt << endl;
          valRascRnp = 0;
        }
        NormsRkRdetRnp[0] += valRkRdetRnp*ddt;
        NormsRascRnp[0]   += valRascRnp*ddt;
        if(!(j%2)){
          NormsRkRdetRnp[1] += valRkRdetRnp*2*ddt;
          NormsRascRnp[1]   += valRascRnp*2*ddt;
        }
        if(!(j%4)){
          NormsRkRdetRnp[2] += valRkRdetRnp*4*ddt;
          NormsRascRnp[2]   += valRascRnp*4*ddt;
        }
        if(!(j%5)){
          NormsRkRdetRnp[3] += valRkRdetRnp*5*ddt;
          NormsRascRnp[3]   += valRascRnp*5*ddt;
        }
        if(!(j%10)){
          NormsRkRdetRnp[4] += valRkRdetRnp*10*ddt;
          NormsRascRnp[4]   += valRascRnp*10*ddt;
        }
        if(!(j%20)){
          NormsRkRdetRnp[5] += valRkRdetRnp*20*ddt;
          NormsRascRnp[5]   += valRascRnp*20*ddt;
        }
        if(!(j%40)){
          NormsRkRdetRnp[6] += valRkRdetRnp*40*ddt;
          NormsRascRnp[6]   += valRascRnp*40*ddt;
//          cout << "dt = " << dt
        }
      }
      for(int j=NBinnings-1; j>=0; j--){
        cout << "Nbins = " << Binnings[j] << ", NormRkRdetRnp = " << NormsRkRdetRnp[j] << ", NormRascRnp = " << NormsRascRnp[j] << endl;
      }
      cout << endl;
    }
  }
  return;
}

void DrawBinsScanCPV(const double val_vec[2][8], const double err_vec[2][8]){
  const int NBins = 8;
  const double bins[NBins]      = {1,2,3,4,5,6,7,8};
  const double bins_err[NBins]  = {0,0,0,0,0,0,0,0};

  gStyle->SetOptFit(111);
  stringstream out;
  out.str("");
  out << "pics/sin_bscan_m" << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
  if(perfbin) out << "_pb";
  if(perftag) out << "_pt";
  if(m_genfit) out << "_genfit";
  const string sin_fname = out.str();
  out.str("");
  out << "pics/cos_bscan_m" << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
  if(perfbin) out << "_pb";
  if(perftag) out << "_pt";
  if(m_genfit) out << "_genfit";
  const string cos_fname = out.str();

  TCanvas *c_sin_bins_scan = new TCanvas("c_sin_bins_scan","c_sin_bins_scan",400,400);
  c_sin_bins_scan->SetGrid();
  c_sin_bins_scan->Draw();
  TGraphErrors* gr_sin = new TGraphErrors(NBins,bins,val_vec[0],bins_err,err_vec[0]);
  out.str("");
  out << "sin offset (10^{-2})" << get_label(m_mode);
  gr_sin->SetTitle(out.str().c_str());
  gr_sin->SetMarkerStyle(20);
  gr_sin->SetMarkerSize(1.);
  gr_sin->SetMarkerColor(kBlue);
  gr_sin->GetYaxis()->SetRangeUser(-300.,300.);
  gr_sin->GetYaxis()->SetLabelSize(0.06);
  gr_sin->GetXaxis()->SetLabelSize(0.06);
  gr_sin->GetXaxis()->SetTitle("Dalitz bin");
  gr_sin->GetXaxis()->SetTitleSize(0.06);
  gr_sin->GetXaxis()->SetTitleOffset(0.8);
  gr_sin->Draw("ap");
  gr_sin->Fit("pol0");
  c_sin_bins_scan->Update();
  out.str("");
  out << sin_fname << ".root";
  c_sin_bins_scan->Print(out.str().c_str());
  out.str("");
  out << sin_fname << ".eps";
  c_sin_bins_scan->Print(out.str().c_str());
  out.str("");
  out << "evince " << sin_fname << ".eps &";
  system(out.str().c_str());

  TCanvas *c_cos_bins_scan = new TCanvas("c_cos_bins_scan","c_cos_bins_scan",400,400);
  c_cos_bins_scan->SetGrid();
  c_cos_bins_scan->Draw();
  TGraphErrors* gr_cos = new TGraphErrors(NBins,bins,val_vec[1],bins_err,err_vec[1]);
  out.str("");
  out << "cos offset (10^{-2})" << get_label(m_mode);
  gr_cos->SetTitle(out.str().c_str());
  gr_cos->SetMarkerStyle(20);
  gr_cos->SetMarkerSize(1.);
  gr_cos->SetMarkerColor(kBlue);
  gr_cos->GetYaxis()->SetRangeUser(-300.,300.);
  gr_cos->GetYaxis()->SetLabelSize(0.06);
  gr_cos->GetXaxis()->SetLabelSize(0.06);
  gr_cos->GetXaxis()->SetTitle("Dalitz bin");
  gr_cos->GetXaxis()->SetTitleSize(0.06);
  gr_cos->GetXaxis()->SetTitleOffset(0.8);
  gr_cos->Draw("ap");
  gr_cos->Fit("pol0");
  c_cos_bins_scan->Update();
  out.str("");
  out << cos_fname << ".root";
  c_cos_bins_scan->Print(out.str().c_str());
  out.str("");
  out << cos_fname << ".eps";
  c_cos_bins_scan->Print(out.str().c_str());
  out.str("");
  out << "evince " << cos_fname << ".eps &";
  system(out.str().c_str());
  return;
}

void DrawBinsScanLifetime(const double val_vec[8], const double err_vec[8]){
  const int NBins = 8;
  const double bins[NBins]      = {1,2,3,4,5,6,7,8};
  const double bins_err[NBins]  = {0,0,0,0,0,0,0,0};

  gStyle->SetOptFit(111);
  stringstream out;
  out.str("");
  out << "pics/tau_bscan_m" << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
  if(m_genfit) out << "_genfit";
  const string tau_fname = out.str();

  TCanvas *c_tau_bins_scan = new TCanvas("c_tau_bins_scan","c_tau_bins_scan",400,400);
  c_tau_bins_scan->SetGrid();
  c_tau_bins_scan->Draw();
  TGraphErrors* gr_tau = new TGraphErrors(NBins,bins,val_vec,bins_err,err_vec);
  out.str("");
  out << "#tau offset (10^{-2} ps)" << get_label(m_mode);
  gr_tau->SetTitle(out.str().c_str());
  gr_tau->SetMarkerStyle(20);
  gr_tau->SetMarkerSize(1.);
  gr_tau->SetMarkerColor(kBlue);
  gr_tau->GetYaxis()->SetRangeUser(-10.,10.);
  gr_tau->GetYaxis()->SetLabelSize(0.06);
  gr_tau->GetXaxis()->SetLabelSize(0.06);
  gr_tau->GetXaxis()->SetTitle("Dalitz bin");
  gr_tau->GetXaxis()->SetTitleSize(0.06);
  gr_tau->GetXaxis()->SetTitleOffset(0.8);
  gr_tau->Draw("ap");
  gr_tau->Fit("pol0");
  c_tau_bins_scan->Update();
  out.str("");
  out << tau_fname << ".root";
  c_tau_bins_scan->Print(out.str().c_str());
  out.str("");
  out << tau_fname << ".eps";
  c_tau_bins_scan->Print(out.str().c_str());
  out.str("");
  out << "evince " << tau_fname << ".eps &";
  system(out.str().c_str());

  return;
}

double GetBkgMajorant(TTree* SigTree, TRandom3& rndm, const int NTries = 1000){
  double maj = 0;
  double pdf_val = 0;
  double dt;
  Double_t sz_sig,sz_asc;
  SigTree->SetBranchAddress("sz_sig",&sz_sig);
  SigTree->SetBranchAddress("sz_asc",&sz_asc);
  for(int i=0; i<NTries; i++){
    SigTree->GetEvent(i);
    double sigma = 0.1*sum_sigma(sz_sig,sz_asc);
    dt = dtmin + (dtmax-dtmin)*rndm.Rndm();
    pdf_val = m_pdf_back_svd2->Pdf(dt,sigma,1);
    if(pdf_val>maj) maj = pdf_val;
  }
  cout << "Majorant is " << maj << endl;
  return maj;
}

int GenerateToyBackground(TTree* SigTree,const int _mode){
  const int NTot = SigTree->GetEntries();

  TRandom3 rndm;
  rndm.SetSeed(0);
  double Maj = 1.2;//GetBkgMajorant(SigTree,pdf,rndm,1000);

  Double_t sz_sig;
  Double_t sz_asc;
  Double_t chisq_asc;
  Double_t chisq_sig;
  Int_t exp;
  Int_t ndf_rec;
  Int_t ntrk_rec;
  Int_t ndf_asc;
  Int_t ntrk_asc;
  Double_t tag;
  Int_t Bin;
  Int_t flv;
  Double_t costhBcms;
  Double_t chi2_vtx_d0;
  Double_t de;
  Double_t mbc;
  Int_t good_icpv;

  SigTree->SetBranchAddress("sz_sig",&sz_sig);
  SigTree->SetBranchAddress("sz_asc",&sz_asc);
  SigTree->SetBranchAddress("chisq_z_sig",&chisq_sig);
  SigTree->SetBranchAddress("chisq_z_asc",&chisq_asc);
  SigTree->SetBranchAddress("good_icpv",&good_icpv);
  SigTree->SetBranchAddress("exp",&exp);
  SigTree->SetBranchAddress("ndf_z_sig",&ndf_rec);
  SigTree->SetBranchAddress("ndf_z_asc",&ndf_asc);
  SigTree->SetBranchAddress("ntrk_sig",&ntrk_rec);
  SigTree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  SigTree->SetBranchAddress("de",&de);
  SigTree->SetBranchAddress("mbc",&mbc);
//  if(m_ebeb) SigTree->SetBranchAddress("f_bkg",&m_f_bkg);
  SigTree->SetBranchAddress("costhBcms",&costhBcms);
  SigTree->SetBranchAddress("bin_mc",&Bin);
  SigTree->SetBranchAddress("flv_mc",&flv);
  SigTree->SetBranchAddress("tag_LH",&tag);
  SigTree->SetBranchAddress("chi2_vtx_d0",&chi2_vtx_d0);

  Double_t z_asc = 0;
  Double_t z_sig;
  Int_t b0f = -1;
  Int_t h0f = 1;

  TFile* file = new TFile("toy_bkg_tree.root","RECREATE");

  int mode = 1;
  int h0mode = 10;
  switch(_mode){
  case 1:// D0 pi0
    mode = 1; h0mode = 10; break;
  case 2:// D0 eta->gg
    mode = 2; h0mode = 10; break;
  case 3:// D0 eta->pi+pi-pi0
    mode = 2; h0mode = 20; break;
  case 4:// D0 omega->pi+pi-pi0
    mode = 3; h0mode = 20; break;
  case 5:// D0 omega->pi+pi-pi0 (Ks rho model)
    mode = 3; h0mode = 20; break;
  case 6:// D0 rho
    h0mode = 40; mode = 4; break;
  }

  TTree* tree = new TTree("TEvent","TEvent");
  tree->Branch("exp",&exp,"exp/I");
  tree->Branch("mode",&mode,"mode/I");
  tree->Branch("h0mode",&h0mode,"h0mode/I");
  tree->Branch("chi2_vtx_d0",&chi2_vtx_d0,"chi2_vtx_d0/D");

  tree->Branch("de",&de,"de/D");
  tree->Branch("mbc",&mbc,"mbc/D");
  tree->Branch("costhBcms",&costhBcms,"costhBcms/D");
  tree->Branch("bin_mc",&Bin,"bin_mc/I");
  tree->Branch("flv_mc",&flv,"flv_mc/I");
  tree->Branch("tag_LH",&tag,"tag_LH/D");
  tree->Branch("ntrk_sig",&ntrk_rec,"ntrk_sig/I");
  tree->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  tree->Branch("ndf_z_sig",&ndf_rec,"ndf_z_sig/I");
  tree->Branch("ndf_z_asc",&ndf_asc,"ndf_z_asc/I");
  tree->Branch("z_sig",&z_sig,"z_sig/D");
  tree->Branch("sz_sig",&sz_sig,"sz_sig/D");
  tree->Branch("chisq_z_sig",&chisq_sig,"chisq_z_sig/D");
  tree->Branch("good_icpv",&good_icpv,"good_icpv/I");
  tree->Branch("z_sig_pipi",&z_sig,"z_sig_pipi/D");
  tree->Branch("sz_sig_pipi",&sz_sig,"sz_sig_pipi/D");
  tree->Branch("chisq_sig_pipi",&chisq_sig,"chisq_sig_pipi/D");
  tree->Branch("z_sig_d0",&z_sig,"z_sig_d0/D");
  tree->Branch("sz_sig_d0",&sz_sig,"sz_sig_d0/D");
  tree->Branch("chisq_sig_d0",&chisq_sig,"chisq_sig_d0/D");
  tree->Branch("z_asc",&z_asc,"z_asc/D");
  tree->Branch("sz_asc",&sz_asc,"sz_asc/D");
  tree->Branch("chisq_z_asc",&chisq_asc,"chisq_z_asc/D");
  tree->Branch("b0f",&b0f,"b0f/I");
  tree->Branch("h0f",&h0f,"h0f/I");

  double xi,dt;
  bool flag = false;
  for(int i=0; i<NTot; i++){
    SigTree->GetEvent(i);
    if(!good_icpv) continue;
    double sigma = 0.1*sum_sigma(sz_sig,sz_asc);
    while(!flag){
      dt = dtmin + (dtmax-dtmin)*rndm.Rndm();
      xi = Maj*rndm.Rndm();
      double pdf_val = exp>30 ?  m_pdf_back_svd2->Pdf(dt,sigma,ndf_asc) : m_pdf_back_svd1->Pdf(dt,sigma,ndf_asc);
      if(pdf_val>Maj){
        cout << "Maj update: " << Maj << " -> " << pdf_val << endl;
        Maj = pdf_val;
      }
      if(xi < pdf_val){
        z_sig = dt/(0.1*cm2ps);
//        flv = (rndm.Rndm() > 0.5) ? 1 : -1;
//        tag = 0.9*flv;
//        Bin = bin((int)(16*rndm.Rndm()));
        if(!(i%10000)) cout << "event " << i << " " << z_sig << " " << flv << " " << tag << " " << Bin << endl;
        tree->Fill();
        flag = true;
      }
    }
    flag = false;
  }
  tree->Write();
  file->Write();
  file->Close();
//  tree->Show(0);
//  cout << "Toy Bkg Tree with " << tree->GetEntries() << " event is generated" << endl;
  return 0;
}

double get_f_cont_in_bkg_sig(const int mode, const int h0mode){
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
double get_f_cont_in_bkg_sideband(const int mode, const int h0mode){
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

string GenFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/mixtree_m" << Mode(mode) << "_mh0" << h0Mode(mode) << ".root";
  return out.str();
}

string SigMCFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/sigmc_cpv_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  out << ".root";
  return out.str();
}

string DataFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/data_cpv_tree_m" << Mode(mode) << "_mh0" << h0Mode(mode) << ".root";
  return out.str();
}

string SidebandFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/data_sideband_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  out << ".root";
  return out.str();
}

string GenSidebandFile(const int mode,const int flag = 0){
    // flag == 0 -> all events
    // flag == 1 -> continuum
    // flag == 2 -> BB
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/genmc_sideband_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  if(flag == 1)      out << "_cont";
  else if(flag == 2) out << "_BB";
  out << ".root";
  return out.str();
}

string GenTestFile(const int mode, const int flag = 0){
    // flag == 0 -> all events
    // flag == 1 -> continuum
    // flag == 2 -> BB
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/genmc_cpv_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  if(flag == 1)      out << "_cont";
  else if(flag == 2) out << "_BB";
  out << ".root";
  return out.str();
}

string GenWWFile(const int mode, const int flag,const int nstr, const int cstr){
    // flag == 0 -> all events
    // flag == 1 -> continuum
    // flag == 2 -> BB
  stringstream out;
  out.str("");
//  out << "/home/vitaly/B0toDh0/PurityFit/data/genmc_cpv_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  out << "/home/vitaly/B0toDh0/PurityFit/data/mixcpvtree_m" << Mode(mode) << "_mh0" << h0Mode(mode);
  if(flag == 1)      out << "_cont";
  else if(flag == 2) out << "_BB";
  out << "_ns" << nstr;
  out << "_cs" << cstr;
//  out << "_ww.root";
//  out << "_smpl1";
  out << ".root";
  cout << "GenWWFile: " << out.str() << endl;
  return out.str();
}

string SigLineFile(const int mode, const string& angle){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/sigmc_cpv_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  out << "_line" << angle << ".root";
  return out.str();
}

string OneStreamGenWWFile(const int mode){
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/mixcpvtree_m" << Mode(mode) << "_mh0" << h0Mode(mode);
  out << "_one_stream.root";
  return out.str();
}

int WriteBkgParameters(const MnUserParameterState& pstate,const int svd){
  stringstream out;
  out.str("");
  out << "params/BackParams_m" << Mode(m_mode) << "_mh0" << h0Mode(m_mode) << "_svd" << svd;
  if(!m_data) out << "_mc";
  if(m_type_flag == 1)      out << "_cont";
  else if(m_type_flag == 2) out << "_BB";
  if(m_gg) out << "_gg";
  if(m_ppp) out << "_ppp";
  out << ".txt";
  cout << "Saving parameters in file " << out.str() << endl;
  ofstream ofile;
  ofile.open(out.str().c_str(),ofstream::out);
  ofile << pstate.Name(0) << " = " << pstate.Value(0) << " +- " << pstate.Error(0) << endl;
  const int imin = svd == 1 ? 1  : 13;
  const int imax = svd == 1 ? 13 : 25;
  for(int i=imin; i<imax; i++){
    ofile << pstate.Name(i) << " = " << pstate.Value(i) << " +- " << pstate.Error(i) << endl;
  }
  ofile.close();
  return 0;
}

int read_com_line_params(const int argn, const char** argv){
  string str;
  for(int i=1; i<argn; i++){
    str = string(argv[i]);
    if(string("-b") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_fitbin)) return -2;
      if(abs(m_fitbin)>8) return -3;
      continue;
    }
    if(string("svd") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_svd)) return -2;
      if(m_svd != 1 && m_svd != 2) return -3;
      cout << "SVD " << m_svd << " setted" << endl;
      continue;
    }
    if(string("-m") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_mode)) return -2;
//      if(m_mode < 1 || m_mode > 6) return -3;
    }
    if(string("-f") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_fitflv)) return -2;
      if(m_fitflv != 1 && m_fitflv != -1) return -3;
    }
    if(string("nsig") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_NSig)) return -2;
      if(!(m_NSig>0)) return -3;
      m_NBkg = m_NSig/m_fSig-m_NSig;
    }
    if(string("toy") == str){
//      if(++i == argn) return -1;
//      if(1 != sscanf(argv[i],"%d",&m_toysize)) return -2;
//      if(!(m_toysize>0)) return -3;
      no_bkg = true; sigmc = true;
      m_toyfit = true;
    }
    if(string("fsig") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%lf",&m_fSig)) return -2;
      if(m_fSig<0 || m_fSig>1) return -3;
      m_NBkg = m_NSig/m_fSig-m_NSig;
    }
    if(string("nsets") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_nsets)) return -2;
    }
    if(string("normtest") == str){
      if(++i == argn) return -1;
      m_norm_test_flag = true;
      if(1 != sscanf(argv[i],"%d",&m_Neve)) return -2;
      if(m_Neve<1) return -3;
    }
    if(string("abfit") == str){
      ABfit = true; FullFit = false;
      if(!m_fitflv) m_fitflv = 1;
    }
    if(string("full") == str){
      ABfit = false; FullFit = true;
    }
    if(string("toybkg") == str){
      toybkg = true;
    }
    if(string("Rsig") == str){
      make_Rrec_fit = true;
    }
    if(string("Rtag") == str){
      make_Rasc_fit = true;
    }
    if(string("Rnp") == str){
      make_Rnp_fit = true;
    }
    if(string("pipi") == str){
      pipi_fit = true;
    }
    if(string("d0") == str){
      d0_fit = true;
    }
    if(string("noint") == str){
      no_interf = true;
      fix_btau = false;
      fix_cos2beta = true;
      fix_sin2beta = true;
    }
    if(string("nonp") == str){
      no_np = true;
    }
    if(string("bscan") == str){
      make_bins_scan = true;
      full_ds_fit = true;
    }
    if(string("fixtau") == str){
      fix_btau = true;
    }
    if(string("freetau") == str){
      fix_btau = false;
    }
    if(string("fixdm") == str){
      fix_dm = true;
    }
    if(string("freedm") == str){
      fix_dm = false;
    }
    if(string("fixsin") == str){
      fix_sin2beta = true;
    }
    if(string("freesin") == str){
      fix_sin2beta = false;
    }
    if(string("fixcos") == str){
      fix_cos2beta = true;
    }
    if(string("freecos") == str){
      fix_cos2beta = false;
    }
    if(string("fixA") == str){
      fix_A = true;
    }
    if(string("fixB") == str){
      fix_B = true;
    }
    if(string("fix_pdf") == str){
      fix_pdf = true;
    }
    if(string("cleo") == str){
      m_cleo = true;
    }
    if(string("sgl") == str){
      sgl_asc = true;
    }
    if(string("mlt") == str){
      mlt_asc = true;
    }
    if(string("s0") == str){
      include_s0 = true;
    }
    if(string("dt0") == str){
      include_dt0 = true;
    }
    if(string("nobkg") == str){
      no_bkg = true;
    }
    if(string("fullds") == str){
      full_ds_fit = true;
    }
    if(string("perftag") == str){
      perftag = true;
    }
    if(string("perfbin") == str){
      perfbin = true;
    }
    if(string("dtmax") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%lf",&dtmax)) return -2;
      dtmax = fabs(dtmax);
      dtmin = -dtmax;
    }
    if(string("sideband") == str){
      sideband_fit = true;
    }
    if(string("Olr") == str){
      add_otlr = true;
    }
    if(string("ebeb") == str){
      m_ebeb = true;
    }
    if(string("data") == str){
      m_data = true;
    }
    if(string("cont") == str){
      m_type_flag = 1;
    }
    if(string("BB") == str){
      m_type_flag = 2;
    }
    if(string("ww") == str){
      m_ww = true;
    }
    if(string("gg") == str){
      m_gg = true;
    }
    if(string("ppp") == str){
      m_ppp = true;
    }
    if(string("sigmc") == str){
      sigmc = true; no_bkg = true;
    }
    if(string("calcK") == str){
      m_calc_K = true;
      sigmc = true;
      no_bkg = true;
    }
    if(string("calcCS") == str){
      m_calc_CS = true;
      sigmc = true;
      no_bkg = true;
    }
    if(string("nuis") == str){
      m_nuisance = true;
    }
    if(string("genfit") == str){
      cout << "genfit" << endl;
      m_genfit = true;
      no_bkg = true;
      sigmc = true;
    }
    if(string("line") == str){
      m_line_test = true;
    }
    if(string("ns") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_ns)) return -2;
    }
    if(string("cs") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_cs)) return -2;
    }
  }
  return 0;
}

void DrawCPV(const int NFreePar){
  stringstream out;
  int NTot = m_ww_tree->GetEntries();
  const double ddt  = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  const int nBD = NDots/NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];
  for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }

  int Nev[2][16];
  TH1I* dh[2][16];
  double pdf_arr[2][16][NDots];
  double norm[2][16];

  int NEveCounter = 0;

  // Initialization
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      Nev[k][i]  = 0;
      norm[k][i] = 0;
      out.str("");
      out << "dh" << bin(i) << "_" << k;
      dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dtmin,dtmax);
      dh[k][i]->SetMarkerStyle(21);
      dh[k][i]->SetMarkerSize(1.1);
      !k ? dh[k][i]->SetMarkerColor(kBlue) : dh[k][i]->SetMarkerColor(kRed);
      for(int j=0; j<NDots; j++) pdf_arr[k][i][j] = 0;
    }
  }
  // Filling
  int flv_index, bin_index;
  for(int l=0; l<NTot; l++){
    GetEventWW(l);
    flv_index = perftag ? flv_ind(m_flv_mc) : flv_ind(m_flv);
    bin_index = perfbin ? bin_ind(m_bin_mc) : bin_ind(m_bin);
    dh[flv_index][bin_index]->Fill(m_dt);
    Nev[flv_index][bin_index]++;
    NEveCounter++;
    for(int j=0; j<NDots; j++){
      double& dt = dt_arr[j];
      double pdfval = Pdf(dt);
      if(!std::isnan(pdfval)){ pdf_arr[flv_index][bin_index][j] += pdfval;}
      else{
        cout << "pdf is nan: dt = " << dt << ", flv: " << m_flv << ", bin: " << m_bin << endl;
        cout << "  m_ntrk_rec: " << m_ntrk_rec << ", m_sz_rec: " << m_sz_rec << ", m_chisq_rec: " << m_chisq_rec << ", m_ndf_rec: " << m_ndf_rec << endl;
        cout << "  m_ntrk_asc: " << m_ntrk_asc << ", m_sz_asc: " << m_sz_asc << ", m_chisq_asc: " << m_chisq_asc << ", m_ndf_asc: " << m_ndf_asc << endl;
      }
    }
  }
  cout << "NEveCounter = " << NEveCounter << endl;

  double pull_array[2][16][NBins];
  double chisq[2][16];

  cout << "Set norm" << endl;
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      chisq[k][i] = 0;
      for(int j=0; j<NDots; j++){
        pdf_arr[k][i][j] *= ddtb;
        norm[k][i] += pdf_arr[k][i][j];
        if(!(j%nBD)){
          const int bin = j/nBD;
          const int bin_content = dh[k][i]->GetBinContent(bin+1);
          if(bin_content){
            pull_array[k][i][bin] = (bin_content-pdf_arr[k][i][j])/sqrt(bin_content);
            chisq[k][i] += pull_array[k][i][bin]*pull_array[k][i][bin];
          }
        }
      }
    }
  }

  TGraph* gr[2][16];
  TGraphErrors* grp[2][16];
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      gr[k][i] = new TGraph(NDots,dt_arr,pdf_arr[k][i]);
      gr[k][i]->SetMarkerStyle(kDot);
      gr[k][i]->SetMarkerSize(1.5);
      gr[k][i]->SetLineWidth(2);
      k == 0 ? gr[k][i]->SetMarkerColor(kBlue) : gr[k][i]->SetMarkerColor(kRed);

      grp[k][i] = new TGraphErrors(NBins,dt_barr,pull_array[k][i],dt_barr_err,dt_pull_err);
      grp[k][i]->SetMarkerStyle(21);
      grp[k][i]->SetMarkerSize(1.2);
      k == 0 ? grp[k][i]->SetMarkerColor(kBlue) : grp[k][i]->SetMarkerColor(kRed);
      grp[k][i]->SetLineColor(0);
    }
  }

  cout << "Draw" << endl;
  TCanvas* c1 = new TCanvas("c1","c1",2600,1000);
  TPad *pad[8];
  TPad *padb[8];
  TPad *padp[8];
  TPad *padpb[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad" << i;
    pad[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.61,0.125*(i+1)-0.01,0.99);
    pad[i]->Draw();
    out.str("");
    out << "padb" << i;
    padb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.11,0.125*(i+1)-0.01,0.49);
    padb[i]->Draw();

    out.str("");
    out << "padp" << i;
    padp[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.59);
    padp[i]->Draw();
    out.str("");
    out << "padpb" << i;
    padpb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.09);
    padpb[i]->Draw();
  }
  int k = 0;
  for(int i=0; i<8; i++){
    pad[i]->cd();
    pad[i]->SetGrid();
    pad[i]->SetLogy();
    dh[k][i]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i]->Draw("e");
    pad[i]->Update();
    TPaveStats *st = (TPaveStats*)dh[k][i]->FindObject("stats");
    st->SetOptStat(1110);
    st->PaintPave(0.55,0.75,0.99,0.9);
    pad[i]->Update();
    gr[k][i]->Draw("same");

    padb[i]->cd();
    padb[i]->SetGrid();
    padb[i]->SetLogy();
    dh[k][i+8]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i+8]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i+8]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i+8]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i+8]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i+8]->Draw("e");
    padb[i]->Update();
    st = (TPaveStats*)dh[k][i+8]->FindObject("stats");
    st->SetOptStat(1110);
    st->PaintPave(0.6,0.7,0.99,0.9);
    padb[i]->Update();
    gr[k][i+8]->Draw("same");

    padp[i]->cd();
    grp[k][i]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i]->Draw();
    TLine* zeroline = new TLine(dtmin,0,dtmax,0);
    zeroline->SetLineWidth(2);
    zeroline->Draw("AP");

    TLine* plusline = new TLine(dtmin,3,dtmax,3);
    plusline->SetLineWidth(1);
    plusline->SetLineStyle(kDashed);
    plusline->Draw();

    TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
    minusline->SetLineWidth(1);
    minusline->SetLineStyle(kDashed);
    minusline->Draw();

    padpb[i]->cd();
    grp[k][i+8]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i+8]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i+8]->Draw();
    TLine* zerolineb = new TLine(dtmin,0,dtmax,0);
    zerolineb->SetLineWidth(2);
    zerolineb->Draw("AP");

    TLine* pluslineb = new TLine(dtmin,3,dtmax,3);
    pluslineb->SetLineWidth(1);
    pluslineb->SetLineStyle(kDashed);
    pluslineb->Draw();

    TLine* minuslineb = new TLine(dtmin,-3,dtmax,-3);
    minuslineb->SetLineWidth(1);
    minuslineb->SetLineStyle(kDashed);
    minuslineb->Draw();
  }

  out.str("");
  out << "pics/full_cpfit_m";
  if(!m_gg && !m_ppp) out << m_mode;
  else if(m_gg) out << "gg";
  else out << "ppp";
  if(mlt_asc)         out << "_mlt_";
  if(sgl_asc)         out << "_sgl_";
  if(m_svd == 2)      out << "svd2";
  else if(m_svd == 1) out << "svd1";
  if(!NFreePar)       out << "_def";
  if(no_bkg)          out << "_nobkg";
  if(perftag)         out << "_perftag";
  if(!full_ds_fit){
    out << "_nsig_" << m_NSig;
    out << "_fsig_" << m_fSig;
  }
  string rootname =   out.str() + string(".root");

  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  TCanvas* c2 = new TCanvas("c2","c2",2600,1000);
  c2->cd();
  TPad *pad1[8];
  TPad *padb1[8];
  TPad *padp1[8];
  TPad *padpb1[8];
  for(int i=0; i<8; i++){
    out.str("");
    out << "pad1" << i;
    pad1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.61,0.125*(i+1)-0.01,0.99);
    pad1[i]->Draw();
    out.str("");
    out << "padb1" << i;
    padb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.11,0.125*(i+1)-0.01,0.49);
    padb1[i]->Draw();

    out.str("");
    out << "padp1" << i;
    padp1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.59);
    padp1[i]->Draw();
    out.str("");
    out << "padpb1" << i;
    padpb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.09);
    padpb1[i]->Draw();
  }
  k = 1;
  for(int i=0; i<8; i++){
    pad1[i]->cd();
//    gStyle->SetOptStat("RMe");
    pad1[i]->SetGrid();
    pad1[i]->SetLogy();
    dh[k][i]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i]->Draw("e");
    pad1[i]->Update();
    TPaveStats *st = (TPaveStats*)dh[k][i]->FindObject("stats");
    st->SetOptStat(2210);
    st->PaintPave(0.60,0.65,0.99,0.9);
    pad1[i]->Update();
    gr[k][i]->Draw("same");

    padb1[i]->cd();
    padb1[i]->SetGrid();
    padb1[i]->SetLogy();
    dh[k][i+8]->GetXaxis()->SetTitle("#Deltat (ps)");
    dh[k][i+8]->GetXaxis()->SetTitleSize(0.06);
    dh[k][i+8]->GetXaxis()->SetTitleOffset(0.75);
    dh[k][i+8]->GetXaxis()->SetLabelSize(0.05);
    dh[k][i+8]->GetYaxis()->SetLabelSize(0.05);
    dh[k][i+8]->Draw("e");
    padb1[i]->Update();
    st = (TPaveStats*)dh[k][i+8]->FindObject("stats");
    st->SetOptStat(2210);
    st->PaintPave(0.65,0.65,0.99,0.9);
    padb1[i]->Update();
    gr[k][i+8]->Draw("same");

    padp1[i]->cd();
    grp[k][i]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i]->Draw();
    TLine* zeroline1 = new TLine(dtmin,0,dtmax,0);
    zeroline1->SetLineWidth(2);
    zeroline1->Draw("AP");

    TLine* plusline = new TLine(dtmin,3,dtmax,3);
    plusline->SetLineWidth(1);
    plusline->SetLineStyle(kDashed);
    plusline->Draw();

    TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
    minusline->SetLineWidth(1);
    minusline->SetLineStyle(kDashed);
    minusline->Draw();

    padpb1[i]->cd();
    grp[k][i+8]->GetYaxis()->SetRangeUser(-5,5);
    grp[k][i+8]->GetXaxis()->SetRangeUser(dtmin,dtmax);
    grp[k][i+8]->Draw();
    TLine* zerolineb = new TLine(dtmin,0,dtmax,0);
    zerolineb->SetLineWidth(2);
    zerolineb->Draw("AP");

    TLine* pluslineb = new TLine(dtmin,3,dtmax,3);
    pluslineb->SetLineWidth(1);
    pluslineb->SetLineStyle(kDashed);
    pluslineb->Draw();

    TLine* minuslineb = new TLine(dtmin,-3,dtmax,-3);
    minuslineb->SetLineWidth(1);
    minuslineb->SetLineStyle(kDashed);
    minuslineb->Draw();
  }

  out.str("");
  out << "pics/full_cpfit_m" << m_mode;
  if(mlt_asc)         out << "_mlt_";
  if(sgl_asc)         out << "_sgl_";
  if(m_svd == 2)      out << "svd2";
  else if(m_svd == 1) out << "svd1";
  if(!NFreePar)       out << "_def";
  out << "_b0b";
  rootname = out.str() + string(".root");

  c2->Print(rootname.c_str());
  pngname = out.str() + string(".png");
  c2->Print(pngname.c_str());

  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());

  cout << "NEveCounter = " << NEveCounter << endl;
  return;
}

void SetBranchAddresses(TTree* tree){
  tree->SetBranchAddress("exp",&m_exp);
  tree->SetBranchAddress("bin",&m_bin);
  tree->SetBranchAddress("flv",&m_flv);
  tree->SetBranchAddress("mode",&mm_mode);
  tree->SetBranchAddress("h0mode",&mm_h0mode);
  tree->SetBranchAddress("ndf_sig",&m_ndf_rec);
  tree->SetBranchAddress("ndf_asc",&m_ndf_asc);
  tree->SetBranchAddress("ntrk_sig",&m_ntrk_rec);
  tree->SetBranchAddress("ntrk_asc",&m_ntrk_asc);
  tree->SetBranchAddress("tag_LH",&m_tag);
  tree->SetBranchAddress("z_sig",&m_t_sig);
  tree->SetBranchAddress("sz_sig",&m_sz_rec);
  tree->SetBranchAddress("chisq_sig",&m_chisq_rec);
  tree->SetBranchAddress("z_asc",&m_t_asc);
  tree->SetBranchAddress("sz_asc",&m_sz_asc);
  tree->SetBranchAddress("chisq_asc",&m_chisq_asc);
  tree->SetBranchAddress("costhBcms",&m_costhBcms);
  tree->SetBranchAddress("sigarea",&m_sigarea);
  if(!m_data){
    tree->SetBranchAddress("bin_mc",&m_bin_mc);
    tree->SetBranchAddress("flv_mc",&m_flv_mc);
    tree->SetBranchAddress("b0f",&m_b0f);
  }
  if(sigmc){
    tree->SetBranchAddress("z_sig_mc",&m_t_sig_mc);
    tree->SetBranchAddress("z_asc_mc",&m_t_asc_mc);
  }
  if(!no_bkg){
    tree->SetBranchAddress("f_cont_in_comb",&m_f_cont_in_comb);
    tree->SetBranchAddress("f_cont",&m_f_cont);
    tree->SetBranchAddress("f_bkg",&m_f_bkg);
    if(!m_data){
      tree->SetBranchAddress("f_cont_in_comb_mc",&m_f_cont_in_comb_mc);
      tree->SetBranchAddress("f_cont_mc",&m_f_cont_mc);
      tree->SetBranchAddress("f_bkg_mc",&m_f_bkg_mc);

      tree->SetBranchAddress("f_cont_in_comb_bin_mc",&m_f_cont_in_comb_bin_mc);
      tree->SetBranchAddress("f_cont_bin_mc",&m_f_cont_bin_mc);
      tree->SetBranchAddress("f_bkg_bin_mc",&m_f_bkg_bin_mc);

      tree->SetBranchAddress("f_cont_in_comb_flv_mc",&m_f_cont_in_comb_flv_mc);
      tree->SetBranchAddress("f_cont_flv_mc",&m_f_cont_flv_mc);
      tree->SetBranchAddress("f_bkg_flv_mc",&m_f_bkg_flv_mc);
    }
  }
  tree->SetBranchAddress("mp",&m_mp);
  tree->SetBranchAddress("mm",&m_mm);
  tree->SetBranchAddress("mp_mc",&m_mp_mc);
  tree->SetBranchAddress("mm_mc",&m_mm_mc);
  if(!sigmc && !m_data){
    tree->SetBranchAddress("Nsig",&m_Nsig);
    tree->SetBranchAddress("Nsig_err",&m_Nsig_err);
    tree->SetBranchAddress("Ntot",&m_Ntot);
    tree->SetBranchAddress("psig",&m_psig);
    tree->SetBranchAddress("pcnt",&m_pcnt);
    tree->SetBranchAddress("pprt",&m_pprt);
    tree->SetBranchAddress("pcmb",&m_pcmb);
    tree->SetBranchAddress("fbb",&m_fbb);
    tree->SetBranchAddress("fbb_err",&m_fbb_err);
    tree->SetBranchAddress("fbb",&m_fprt);
    tree->SetBranchAddress("fbb_err",&m_fprt_err);
    tree->SetBranchAddress("wrtag",&m_wtag);
  }
  return;
}

int init(const int _mode){
  cout << "init: mode " << _mode << endl;

  xi = xil(_mode);
  beta = Beta(_mode);
  m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
  m_cos2beta = TMath::Cos(2.*beta/180.*TMath::Pi());

  if(m_svd != 2 && !m_norm_test_flag){
    m_pdf_svd1 = new RkRdetRnpPdf(1,1);
    m_pdf_svd1->SetTauDm(m_btau,m_dm);
    m_f_ol_sgl_svd1 = m_pdf_svd1->Get_f_ol_sgl();
    m_f_ol_mlt_svd1 = m_pdf_svd1->Get_f_ol_mlt();
    m_s_ol_svd1     = m_pdf_svd1->Get_sigma_ol();
    cout << "SVD1 Outlier " << m_f_ol_sgl_svd1 << " " << m_f_ol_mlt_svd1 << " " << m_s_ol_svd1 << endl;
    const int MODE = toybkg ? 100+m_mode : m_mode;
    if(!no_bkg){
      if(m_type_flag){
      m_pdf_back_svd1 = m_mode ? new RbkgPdf(MODE,1,true,false) : new RbkgPdf(MODE,1,false,false);
      m_pdf_back_svd1->GetParametersFromFile(Mode(m_mode),h0Mode(m_mode),1,!m_data,0);
      } else{
        cout << "Initializing background PDFs..." << endl;
//        m_pdf_back_svd1_cont = m_mode ? new RbkgPdf(MODE,1,true,false) : new RbkgPdf(MODE,1,false,false);
//        m_pdf_back_svd1_cont->GetParametersFromFile(77,77,1,!m_data,1,ppp_flag(m_mode));
//        m_pdf_back_svd1_bb = m_mode ? new RbkgPdf(MODE,1,true,false) : new RbkgPdf(MODE,1,false,false);
//        m_pdf_back_svd1_bb->GetParametersFromFile(77,77,1,!m_data,2,ppp_flag(m_mode));
//        m_S_main_mlt_cont_svd1 = m_pdf_back_svd1_cont->Get_S_main_mlt();
//        m_S_main_sgl_cont_svd1 = m_pdf_back_svd1_cont->Get_S_main_sgl();
//        m_S_main_mlt_bb_svd1   = m_pdf_back_svd1_bb->Get_S_main_mlt();
//        m_S_main_sgl_bb_svd1   = m_pdf_back_svd1_bb->Get_S_main_sgl();

        m_pdf_back_svd1_cont_gg = new RbkgPdf(MODE,1,true,false);
          m_pdf_back_svd1_bb_gg = new RbkgPdf(MODE,1,true,false);
        m_pdf_back_svd1_cont_gg->GetParametersFromFile(77,77,1,true,1,false);
          m_pdf_back_svd1_bb_gg->GetParametersFromFile(77,77,1,true,2,false);

        m_pdf_back_svd1_cont_ppp = new RbkgPdf(MODE,1,true,false);
          m_pdf_back_svd1_bb_ppp = new RbkgPdf(MODE,1,true,false);
        m_pdf_back_svd1_cont_ppp->GetParametersFromFile(77,77,1,true,1,true);
          m_pdf_back_svd1_bb_ppp->GetParametersFromFile(77,77,1,true,2,true);

        m_S_main_mlt_cont_svd1_gg = m_pdf_back_svd1_cont_gg->Get_S_main_mlt();
        m_S_main_sgl_cont_svd1_gg = m_pdf_back_svd1_cont_gg->Get_S_main_sgl();
          m_S_main_mlt_bb_svd1_gg =   m_pdf_back_svd1_bb_gg->Get_S_main_mlt();
          m_S_main_sgl_bb_svd1_gg =   m_pdf_back_svd1_bb_gg->Get_S_main_sgl();

        m_S_main_mlt_cont_svd1_ppp = m_pdf_back_svd1_cont_ppp->Get_S_main_mlt();
        m_S_main_sgl_cont_svd1_ppp = m_pdf_back_svd1_cont_ppp->Get_S_main_sgl();
          m_S_main_mlt_bb_svd1_ppp =   m_pdf_back_svd1_bb_ppp->Get_S_main_mlt();
          m_S_main_sgl_bb_svd1_ppp =   m_pdf_back_svd1_bb_ppp->Get_S_main_sgl();
      }
    }
    m_pdf_svd1->SetXi(xi);
  }

  if(m_svd != 1 || m_norm_test_flag){
    m_pdf_svd2 = new RkRdetRnpPdf(1,2);
    m_pdf_svd2->SetTauDm(m_btau,m_dm);
    m_f_ol_sgl_svd2 = m_pdf_svd2->Get_f_ol_sgl();
    m_f_ol_mlt_svd2 = m_pdf_svd2->Get_f_ol_mlt();
    m_s_ol_svd2     = m_pdf_svd2->Get_sigma_ol();
    cout << "SVD2 Outlier " << m_f_ol_sgl_svd2 << " " << m_f_ol_mlt_svd2 << " " << m_s_ol_svd2 << endl;
    const int MODE = toybkg ? 100+m_mode : m_mode;
    cout << "MODE: " << MODE << endl;
    if(!no_bkg){
    if(m_type_flag){
      m_pdf_back_svd2 = m_mode ? new RbkgPdf(MODE,2,true,false) : new RbkgPdf(MODE,2,false,false);
      m_pdf_back_svd2->GetParametersFromFile(Mode(m_mode),h0Mode(m_mode),2,!m_data,0);
    } else{
//      m_pdf_back_svd2_cont = m_mode ? new RbkgPdf(MODE,2,true,false) : new RbkgPdf(MODE,2,false,false);
//      m_pdf_back_svd2_cont->GetParametersFromFile(77,77,2,!m_data,1,ppp_flag(m_mode));
//      m_pdf_back_svd2_bb = m_mode ? new RbkgPdf(MODE,2,true,false) : new RbkgPdf(MODE,2,false,false);
//      m_pdf_back_svd2_bb->GetParametersFromFile(77,77,2,!m_data,2,ppp_flag(m_mode));
//      m_S_main_mlt_cont_svd2 = m_pdf_back_svd2_cont->Get_S_main_mlt();
//      m_S_main_sgl_cont_svd2 = m_pdf_back_svd2_cont->Get_S_main_sgl();
//      m_S_main_mlt_bb_svd2   = m_pdf_back_svd2_bb->Get_S_main_mlt();
//      m_S_main_sgl_bb_svd2   = m_pdf_back_svd2_bb->Get_S_main_sgl();
      m_pdf_back_svd2_cont_gg = new RbkgPdf(MODE,2,true,false);
        m_pdf_back_svd2_bb_gg = new RbkgPdf(MODE,2,true,false);
      m_pdf_back_svd2_cont_gg->GetParametersFromFile(77,77,2,true,1,false);
        m_pdf_back_svd2_bb_gg->GetParametersFromFile(77,77,2,true,2,false);

      m_pdf_back_svd2_cont_ppp = new RbkgPdf(MODE,2,true,false);
        m_pdf_back_svd2_bb_ppp = new RbkgPdf(MODE,2,true,false);
      m_pdf_back_svd2_cont_ppp->GetParametersFromFile(77,77,2,true,1,true);
        m_pdf_back_svd2_bb_ppp->GetParametersFromFile(77,77,2,true,2,true);

      m_S_main_mlt_cont_svd2_gg = m_pdf_back_svd2_cont_gg->Get_S_main_mlt();
      m_S_main_sgl_cont_svd2_gg = m_pdf_back_svd2_cont_gg->Get_S_main_sgl();
        m_S_main_mlt_bb_svd2_gg =   m_pdf_back_svd2_bb_gg->Get_S_main_mlt();
        m_S_main_sgl_bb_svd2_gg =   m_pdf_back_svd2_bb_gg->Get_S_main_sgl();

      m_S_main_mlt_cont_svd2_ppp = m_pdf_back_svd2_cont_ppp->Get_S_main_mlt();
      m_S_main_sgl_cont_svd2_ppp = m_pdf_back_svd2_cont_ppp->Get_S_main_sgl();
        m_S_main_mlt_bb_svd2_ppp =   m_pdf_back_svd2_bb_ppp->Get_S_main_mlt();
        m_S_main_sgl_bb_svd2_ppp =   m_pdf_back_svd2_bb_ppp->Get_S_main_sgl();
      }
    }
    m_pdf_svd2->SetXi(xi);
  }

  if(m_svd != 2){
    m_pdf_svd1->SetRange(dtmax);
    if(!no_bkg){
    if(m_type_flag){
      m_pdf_back_svd1->SetRange(dtmax);
      m_pdf_back_svd1->Set_f_otlr(0);
    } else{
//      m_pdf_back_svd1_cont->SetRange(dtmax);
//      m_pdf_back_svd1_cont->Set_f_otlr(0);
//      m_pdf_back_svd1_bb->SetRange(dtmax);
//      m_pdf_back_svd1_bb->Set_f_otlr(0);
      m_pdf_back_svd1_cont_gg->SetRange(dtmax);
      m_pdf_back_svd1_cont_gg->Set_f_otlr(0);
      m_pdf_back_svd1_bb_gg->SetRange(dtmax);
      m_pdf_back_svd1_bb_gg->Set_f_otlr(0);

      m_pdf_back_svd1_cont_ppp->SetRange(dtmax);
      m_pdf_back_svd1_cont_ppp->Set_f_otlr(0);
      m_pdf_back_svd1_bb_ppp->SetRange(dtmax);
      m_pdf_back_svd1_bb_ppp->Set_f_otlr(0);
    }
    }
  }
  if(m_svd != 1){
    m_pdf_svd2->SetRange(dtmax);
    if(!no_bkg){
    if(m_type_flag){
      m_pdf_back_svd2->SetRange(dtmax);
      m_pdf_back_svd2->Set_f_otlr(0);
    } else{
//      m_pdf_back_svd2_cont->SetRange(dtmax);
//      m_pdf_back_svd2_cont->Set_f_otlr(0);
//      m_pdf_back_svd2_bb->SetRange(dtmax);
//      m_pdf_back_svd2_bb->Set_f_otlr(0);
      m_pdf_back_svd2_cont_gg->SetRange(dtmax);
      m_pdf_back_svd2_cont_gg->Set_f_otlr(0);
      m_pdf_back_svd2_bb_gg->SetRange(dtmax);
      m_pdf_back_svd2_bb_gg->Set_f_otlr(0);

      m_pdf_back_svd2_cont_ppp->SetRange(dtmax);
      m_pdf_back_svd2_cont_ppp->Set_f_otlr(0);
      m_pdf_back_svd2_bb_ppp->SetRange(dtmax);
      m_pdf_back_svd2_bb_ppp->Set_f_otlr(0);
    }
    }
  }
  if(no_interf){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_noint[i];
      Sarr[i] = Sarr_noint[i];
      Karr[i] = Karr_noint[i];
      Kbarr[i] = Kbarr_noint[i];
    }
  //} else if(_mode == 5){
   // for(int i=0; i<8; i++){
   //   Carr[i] = Carr_rhoKs[i];
   //   Sarr[i] = Sarr_rhoKs[i];
   //   Karr[i] = Karr_rhoKs[i];
   //   Kbarr[i] = Kbarr_rhoKs[i];
   // }
  } else if(!m_cleo){
    for(int i=0; i<8; i++){
      Carr[i] = Carr_model[i];
      Sarr[i] = Sarr_model[i];
      // Karr[i] = Karr_model[i];
      // Kbarr[i] = Kbarr_model[i];
      if(m_mode != 3 && m_mode != 4 && m_mode != 5){
        Karr[i] = Karr_pi0[i];
        Kbarr[i] = Kbarr_pi0[i];
      } else{
        Karr[i] = Karr_omega[i];
        Kbarr[i] = Kbarr_omega[i];
      }
    }
  } else{
    for(int i=0; i<8; i++){
      Carr[i] = Carr_CLEO[i];
      Sarr[i] = Sarr_CLEO[i];
      Karr[i] = Karr_CLEO[i];
      Kbarr[i] = Kbarr_CLEO[i];
    }
  }
  cout << "Initialized arrays:" << endl;
  for(int i=0;i<8;i++){
    cout << i+1 << " " << Karr[i] << " " << Kbarr[i] << " " << Carr[i] << " " << Sarr[i] << endl;
  }
  for(int i=0; i<8; i++){
    K0[i] = Karr[i];
    Kb0[i] = Kbarr[i];
    C0[i] = Carr[i];
    S0[i] = Sarr[i];
  }
  if(!m_toyfit) SetBranchAddresses(m_tree);

  return 0;
}

void SetBkgScales(const double& scale1, const double& scale2){
  if(!no_bkg && m_svd != 2){
      m_pdf_back_svd1_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd1_gg*scale1);
    m_pdf_back_svd1_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd1_gg*scale1);
      m_pdf_back_svd1_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd1_gg*scale1);
    m_pdf_back_svd1_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd1_gg*scale1);

      m_pdf_back_svd1_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd1_ppp*scale1);
    m_pdf_back_svd1_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd1_ppp*scale1);
      m_pdf_back_svd1_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd1_ppp*scale1);
    m_pdf_back_svd1_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd1_ppp*scale1);
  }
  if(!no_bkg && m_svd != 1){
      m_pdf_back_svd2_bb_gg->Set_S_main_mlt(m_S_main_mlt_bb_svd2_gg*scale2);
    m_pdf_back_svd2_cont_gg->Set_S_main_mlt(m_S_main_mlt_cont_svd2_gg*scale2);
      m_pdf_back_svd2_bb_gg->Set_S_main_sgl(m_S_main_sgl_cont_svd2_gg*scale2);
    m_pdf_back_svd2_cont_gg->Set_S_main_sgl(m_S_main_sgl_bb_svd2_gg*scale2);

      m_pdf_back_svd2_bb_ppp->Set_S_main_mlt(m_S_main_mlt_bb_svd2_ppp*scale2);
    m_pdf_back_svd2_cont_ppp->Set_S_main_mlt(m_S_main_mlt_cont_svd2_ppp*scale2);
      m_pdf_back_svd2_bb_ppp->Set_S_main_sgl(m_S_main_sgl_cont_svd2_ppp*scale2);
    m_pdf_back_svd2_cont_ppp->Set_S_main_sgl(m_S_main_sgl_bb_svd2_ppp*scale2);
  }
  return;
}
