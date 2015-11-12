#ifndef RESCONST_H
#define RESCONST_H

#include <iostream>

using namespace std;

class ResConst{
public:
  ResConst(const int mc, const int svd, const bool charged = false, const int error = 0);

  void Set_ftl_rec_mlt(const double& x1, const double& x2){ftl_rec_mlt[0] = x1; ftl_rec_mlt[1] = x2; return;}
  void Set_Srec(const double& x1, const double& x2){Srec[0] = x1; Srec[1] = x2; return;}
  void Set_Stl_rec_mlt(const double& x){ Stl_rec_mlt = x; return;}
  void Set_ftl_rec(const double& x){ ftl_rec = x; return;}
  void Set_Smn_rec(const double& x){ Smn_rec = x; return;}
  void Set_Stl_rec(const double& x){ Stl_rec = x; return;}

//  double Get_flt_rec_mlt(const int i) const {return (i>=0 && i<2) ? flt_rec_mlt[i] : -1;}
  double Get_Srec(const int i) const {return Srec[i];}
//  double Get_Stl_rec_mlt(void){return Stl_rec_mlt;}
  double Get_ftl_rec(void) const {return ftl_rec;}
  double Get_Smn_rec(void) const {return Smn_rec;}
  double Get_Stl_rec(void) const {return Stl_rec;}

  void Set_ftl_asc_mlt(const double& x1, const double& x2){ftl_asc_mlt[0] = x1; ftl_asc_mlt[1] = x2; return;}
  void Set_Sasc(const double& x1, const double& x2){Sasc[0] = x1; Sasc[1] = x2; return;}
  void Set_Stl_asc_mlt(const double& x){ Stl_asc_mlt = x; return;}
  void Set_ftl_asc(const double& x){ ftl_asc = x; return;}
  void Set_Smn_asc(const double& x){ Smn_asc = x; return;}
  void Set_Stl_asc(const double& x){ Stl_asc = x; return;}
  // * NP tracks *
  void Set_fd_np_mlt(const double& x1, const double& x2){fd_np_mlt[0] = x1; fd_np_mlt[1] = x2; return;}
  void Set_fp_np_mlt(const double& x){fp_np_mlt = x; return;}
  void Set_tau_np_p_mlt(const double& x1, const double& x2){tau_np_p_mlt[0] = x1; tau_np_p_mlt[1] = x2; return;}
  void Set_tau_np_n_mlt(const double& x1, const double& x2){tau_np_n_mlt[0] = x1; tau_np_n_mlt[1] = x2; return;}
  void Set_fd_np_sgl(const double& x1, const double& x2){fd_np_sgl[0] = x1; fd_np_sgl[1] = x2; return;}
  void Set_fp_np_sgl(const double& x){fp_np_sgl = x; return;}
  void Set_tau_np_p_sgl(const double& x1, const double& x2){tau_np_p_sgl[0] = x1; tau_np_p_sgl[1] = x2; return;}
  void Set_tau_np_n_sgl(const double& x1, const double& x2){tau_np_n_sgl[0] = x1; tau_np_n_sgl[1]; return;}
  void Set_fn_np_mlt(const double& x){fn_np_mlt = x; return;}
  void Set_rnp_kink_st(const double& x){rnp_kink_st = x; return;}
  void Set_rnp_kink_xi(const double& x){rnp_kink_xi = x; return;}
  void Set_fd_np_st_mlt(const double& x){fd_np_st_mlt = x; return;}
  void Set_fd_np_xi_mlt(const double& x){fd_np_xi_mlt = x; return;}
  void Set_fd_np_stxi_mlt(const double& x){fd_np_stxi_mlt = x; return;}
  void Set_Snp_global(const double& x){Snp_global = x; return;}
  void Set_tau_np_p_xi_mlt(const double& x){tau_np_p_xi_mlt = x; return;}
  void Set_tau_np_p_stxi_mlt(const double& x){tau_np_p_stxi_mlt = x; return;}
  void Set_tau_np_n_xi_mlt(const double& x){tau_np_n_xi_mlt = x; return;}
  void Set_tau_np_n_stxi_mlt(const double& x){tau_np_n_stxi_mlt = x; return;}
  void Set_Snp(const double& x){Snp = x; return;}

  void Set_f_ol_sgl(const double& x) {f_ol_sgl = x; return;}
  void Set_f_ol_mlt(const double& x) {f_ol_mul = x; return;}
  void Set_sigma_ol(const double& x) {sigma_ol = x; return;}


  double Get_ftl_asc_mlt(const int i) const {return ftl_asc_mlt[i];}
  double Get_Sasc(const int i) const {return Sasc[i];}
  double Get_Stl_asc_mlt(void) const {return Stl_asc_mlt;}
  double Get_ftl_asc(void) const {return ftl_asc;}
  double Get_Smn_asc(void) const {return Smn_asc;}
  double Get_Stl_asc(void) const {return Stl_asc ;}

  // * NP tracks *
  double Get_fd_np_mlt(const int i) const {return fd_np_mlt[i];}
  double Get_fp_np_mlt(void) const {return fp_np_mlt;}
  double Get_tau_np_p_mlt(const int i) const {return tau_np_p_mlt[i];}
  double Get_tau_np_n_mlt(const int i) const {return tau_np_n_mlt[i];}
  double Get_fd_np_sgl(const int i) const {return fd_np_sgl[i];}
  double Get_fp_np_sgl(void) const {return fp_np_sgl;}
  double Get_tau_np_p_sgl(const int i) const {return tau_np_p_sgl[i];}
  double Get_tau_np_n_sgl(const int i) const {return tau_np_n_sgl[i];}
  double Get_fn_np_mlt(void) const {return fn_np_mlt;}
  double Get_rnp_kink_st(void) const {return rnp_kink_st;}
  double Get_rnp_kink_xi(void) const {return rnp_kink_xi;}
  double Get_fd_np_st_mlt(void) const {return fd_np_st_mlt;}
  double Get_fd_np_xi_mlt(void) const {return fd_np_xi_mlt;}
  double Get_fd_np_stxi_mlt(void) const {return fd_np_stxi_mlt;}
  double Get_Snp_global(void) const {return Snp_global;}
  double Get_tau_np_p_xi_mlt(void) const {return tau_np_p_xi_mlt;}
  double Get_tau_np_p_stxi_mlt(void) const {return tau_np_p_stxi_mlt;}
  double Get_tau_np_n_xi_mlt(void) const {return tau_np_n_xi_mlt;}
  double Get_tau_np_n_stxi_mlt(void) const {return tau_np_n_stxi_mlt;}
  double Get_Snp(void) const {return Snp;}

  double Get_f_ol_sgl(void) const {return f_ol_sgl;}
  double Get_f_ol_mlt(void) const {return f_ol_mul;}
  double Get_sigma_ol(void) const {return sigma_ol;}

protected:
  // ** Resolution parameters **
  // * Rec side *
  double ftl_rec_mlt[2];
  double Srec[2];
  double Stl_rec_mlt;
  double ftl_rec;
  double Smn_rec;
  double Stl_rec;

  // * Tag side *
  double ftl_asc_mlt[2];
  double Sasc[2];
  double Stl_asc_mlt;
  double ftl_asc;
  double Smn_asc;
  double Stl_asc;

  // * NP tracks *
  double fd_np_mlt[2];
  double fp_np_mlt;
  double tau_np_p_mlt[2];
  double tau_np_n_mlt[2];
  double fd_np_sgl[2];
  double fp_np_sgl;
  double tau_np_p_sgl[2];
  double tau_np_n_sgl[2];
  double fn_np_mlt;
  double rnp_kink_st;
  double rnp_kink_xi;
  double fd_np_st_mlt;
  double fd_np_xi_mlt;
  double fd_np_stxi_mlt;
  double Snp_global;
  double tau_np_p_xi_mlt;
  double tau_np_p_stxi_mlt;
  double tau_np_n_xi_mlt;
  double tau_np_n_stxi_mlt;
  double Snp;

  // * Outlayer *
  double sigma_ol;
  double f_ol_sgl;
  double f_ol_mul;

};

#endif
