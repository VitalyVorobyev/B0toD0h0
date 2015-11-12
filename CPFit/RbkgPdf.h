#ifndef RbkgPdf_H
#define RbkgPdf_H

#include "cnvl.h"

class RbkgPdf : public cnvl{
public:
  RbkgPdf(const int mode, const int svd = 2, const bool mc = false, const bool tune = false);
  double Pdf(const double& x, const double& s, const int ndf);
  void SetTau(const double& x){tau = x; return;}

  void Set_f_tail_mlt(const double& x){f_tail_mlt = x; return;}
  void Set_S_main_mlt(const double& x){S_main_mlt = x; return;}
  void Set_S_tail_mlt(const double& x){S_tail_mlt = x; return;}
  void Set_f_tail_sgl(const double& x){f_tail_sgl = x; return;}
  void Set_S_main_sgl(const double& x){S_main_sgl = x; return;}
  void Set_S_tail_sgl(const double& x){S_tail_sgl = x; return;}

  void Set_f_delta_mlt(const double& x){f_delta_mlt = x; return;}
  void Set_f_delta_sgl(const double& x){f_delta_sgl = x; return;}
  void Set_mu_delta(const double& x){mu_delta = x; return;}
  void Set_mu(const double& x){mu = x; return;}

  void Set_f_otlr(const double& x){f_otlr = x; return;}
  void Set_s_otlr(const double& x){s_otlr = x; return;}

  double GetTau(void) const {return tau;}
  double Get_f_tail_mlt(void) const {return f_tail_mlt;}
  double Get_S_main_mlt(void){return S_main_mlt;}
  double Get_S_tail_mlt(void){return S_tail_mlt;}
  double Get_f_tail_sgl(void){return f_tail_sgl;}
  double Get_S_main_sgl(void){return S_main_sgl;}
  double Get_S_tail_sgl(void){return S_tail_sgl;}
  double Get_f_delta_mlt(void){return f_delta_mlt;}
  double Get_f_delta_sgl(void){return f_delta_sgl;}
  double Get_mu_delta(void){return mu_delta;}
  double Get_mu(void){return mu;}

  double Get_f_otlr(void){return f_otlr;}
  double Get_s_otlr(void){return s_otlr;}

  int GetParametersFromFile(const int mode, const int h0mode, const int svd, const bool mc, const int type_flag,const int ppp_flag = 0);

private:
  double AddOutlier(const double& x, const double Lin,const double& nLi);

  //Res params
  double f_tail_mlt;
  double S_main_mlt;
  double S_tail_mlt;

  double f_tail_sgl;
  double S_main_sgl;
  double S_tail_sgl;

  //p.d.f. params
  double f_delta_mlt;
  double mu_delta;
  double mu;
  double tau;

  double f_delta_sgl;

  //Outlier
  double s_otlr;
  double f_otlr;

  double cm2ps;
};

#endif // RbkgPdf_H
