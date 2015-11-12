#ifndef FILTER_H
#define FILTER_H

#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
//#include "Int_t.h"
//#include "Double_t.h"

#include "TMVA/Reader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include "../rooksfw/rooksfw.h"
#include "../MyParams/myparams.h"
#include "../DalitzModelStudy/phasespace.h"
#include "icpvevent.h"

class filter{
public:
  filter(const int type_in,const int type_out, const int str_num=0, const std::string& angle = std::string(""));
  void MakeNTuples(void);
private:
  void SetMVA(void);
  void SetInput(const int str_num, const std::string& angle);
  void SetOutput(std::string& toutstr);
  double chisq_cand(const int mode, const int h0mode, const double& md0, const double& mh0, const double& dm = 0);
  int my_decision2(const std::vector<double>& d0mass,const std::vector<double>& h0mass, const std::vector<double>& dmass, const std::vector<int>& b0fvec, const std::vector<int>& modev, const std::vector<int>& h0modev);
  void PrintVariants2(const std::vector<double>& d0mass,const std::vector<double>& h0mass, const std::vector<double>& dmvec, const std::vector<int>& b0fvec, const std::vector<int>& modev, const std::vector<int>& h0modev);
  void SetBranchAddresses(TTree* tree, const bool second_iter);
  void SetBranches(TTree* tree);
  std::string prefix;
  std::string line_out;
  std::string line_in;
  TChain* intree;
//  TChain* insigchain;
//  TTree* insigtree_svd1;
//  TTree* insigtree_svd2;
  TFile* outfile;
  TTree* outtree;
  TTree* mult_outtree;
  std::vector<int> sigindex_svd1;
  std::vector<int> sigindex_svd2;
  std::string insinfile;
  rooksfw* ksfw1_pi0;
  rooksfw* ksfw1_etagg;
  rooksfw* ksfw1_etappp;
  rooksfw* ksfw1_omega;
  rooksfw* ksfw0_pi0;
  rooksfw* ksfw0_etagg;
  rooksfw* ksfw0_etappp;
  rooksfw* ksfw0_omega;

  rooksfw* ksfw0_dst0pi0;
  rooksfw* ksfw0_dst0etagg;
  rooksfw* ksfw0_dst0etappp;
  rooksfw* ksfw0_etapgg;
  rooksfw* ksfw0_etapppp;

  rooksfw* ksfw1_dst0pi0;
  rooksfw* ksfw1_dst0etagg;
  rooksfw* ksfw1_dst0etappp;
  rooksfw* ksfw1_etapgg;
  rooksfw* ksfw1_etapppp;

  double lh(const int mode, const int h0mode, const double& kmm2,const double* kvars,const bool fsf = false);
  double BDT(const int mode, const int h0mode);

  MyParams* cuts;

  TMVA::Reader* reader_pi0;
  TMVA::Reader* reader_gg;
  TMVA::Reader* reader_ppp;
  TMVA::Reader* reader_omega;
  TMVAEvent tmvaevt;
//  Float_t m_costhBcms, m_chi2_mass_d0, m_cos_thr, m_thr_sig, m_h0_chi2, m_egamma, m_cos_hel;
//  Float_t m_p_pi0_h0, m_p_pip_h0, m_p_pim_h0, m_lh0;

  void Filter(void);
  void MultiFilter(void);
  int IsGoodICPV(const int ndf_z_sig, const double& sz_sig, const double& chisq_z_sig,const int ndf_z_asc, const double& sz_asc, const double& chisq_z_asc);
  int my_decision(const std::vector<double>& d0mass,const std::vector<double>& h0mass, const std::vector<int>& b0fvec, const std::vector<int>& modev);
  bool is_decision(const std::vector<int>& b0fvec);
  void PrintVariants(const std::vector<double>& d0mass,const std::vector<double>& h0mass, const std::vector<int>& b0fvec, const std::vector<int>& modev);
  bool M_SIGMC,M_DATA;

  int d0_des;
  int h0_des;
  int totl_des;
  int good_des;
  int ambiguity;
  int n_sig_cand;
  int n_bkg_cand;
  int sig_counter;
  int cand_struct1[10];
  int cand_struct2[10];

  PhaseSpace* Phsp;
  int m_data_type;
  int m_type_in, m_type_out;
  int m_out_mode, m_out_h0mode;

  ICPVEvent evt;
};

#endif // FILTER_H
