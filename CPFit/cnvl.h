#ifndef CNVL_H
#define CNVL_H

#include <iostream>
#include <cmath>
#include <float.h>

#include "RooMath.h"
#include <complex>

using namespace std;
using namespace ROOT;

typedef double (*nXi_t)(const double& t, const double& xd);

class cnvl
{
public:
  cnvl(): ll(-10), ul(10), ak(1), ck(0){
    inv_sqrt2    = 0.707106781186547461715008466853760182857513427734375;
    inv_sqrt_pi  = 0.56418958354775627928034964497783221304416656494140625;
    sqrt_pi      = 1.772453850905515881919427556567825376987457275390625;
    inv_sqrt_2pi = 0.398942280401432702863218082711682654917240142822265625;
    sqrt2        = 1.4142135623730951454746218587388284504413604736328125;

    beta = 0.39114064034485169394; // 0.425/sqrt(1.0+0.425^2)
    mbzero = 5.2794;

    roomath = new RooMath();
  }
  ~cnvl() {delete roomath;}

  // ** Interface parameters for inherits classes **

  int SetRange(const double& dtmax);
  int SetAkCk(const double& costh, const double& ecm);
  int SetTauDm(const double& _tau, const double& _dm);
  int SetKKCS(const double& _K, const double& _Kb, const double& _C, const double& _S);
  int SetAB(const double& _A, const double& _B);
  int SetFlvXi(const int _flv, const int _xi);
  int SetXi(const int _xi);
  int SetSinCos(const double& sin, const double& cos);

  double GetSin2Beta(void) const {return sin2beta;}
  double GetCos2Beta(void) const {return cos2beta;}

protected:
  double recexp(const double& re, const double& im);
  double imcexp(const double& re, const double& im);
  double rewerf(const double& re, const double& im);
  double imwerf(const double& re, const double& im);

  double Ep(const double& t, const double& tau);
  double En(const double& t, const double& tau);
  double Ef(const double& t, const double& tau);
  double Enp(const double& t, const double& tau_n, const double& tau_p);
  double xEp(const double& t, const double& tau);
  double xEn(const double& t, const double& tau);
  double xEf(const double& t, const double& tau);
  static double nMp(const double& t, const double& xd);
  double Mp(const double& t, const double& tau, const double& dm);
  static double nMn(const double& t, const double& xd);
  double Mn(const double& t, const double& tau, const double& dm);
  double nMf(const double& t, const double& xd);
  double Mf(const double& t, const double& tau, const double& dm);
  static double nAp(const double& t, const double& xd);
  double Ap(const double& t, const double& tau, const double& dm);
  static double nAn(const double& t, const double& xd);
  double An(const double& t, const double& tau, const double& dm);
  double nAf(const double& t, const double& xd);
  double Af(const double& t, const double& tau, const double& dm);

  double norm_nEp(const double & _ll, const double & _ul, const double& o = 0);
  double norm_nEn(const double & _ll, const double & _ul, const double& o = 0);
  double norm_nEf(const double & _ll, const double & _ul, const double& o = 0);
  double norm_Ep(const double& _ll, const double& _ul, const double& tau, const double& o = 0);
  double norm_En(const double& _ll, const double& _ul, const double& tau, const double& o = 0);
  double norm_Ef(const double& _ll, const double& _ul, const double& tau, const double& o = 0);
  double norm_Ap(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0);
  double norm_Ax_sup(const double& x, const double& tau, const double& dm);
  double norm_Mp(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0);
  double norm_Mn(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0);
  double norm_Mx_sup(const double& x, const double& tau, const double& dm);
  double norm_An(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o = 0);

  double nEp_conv_gauss(const double& t, const double& m, const double& s);
  double nEn_conv_gauss(const double& t, const double& m, const double& s);
  double Ep_conv_gauss(const double& t, const double& tau, const double& m, const double& s);
  double En_conv_gauss(const double& t, const double& tau, const double& m, const double& s);
  double Ef_conv_gauss(const double& t, const double& tau, const double& m, const double& s);
  double Enp_conv_gauss(const double& t, const double& tau_n, const double& tau_p, const double& m, const double& s);
  double approx_exp2erfc(const double& x);

  double nMp_conv_gauss(const double& t, const double& xd, const double& m, const double& s);
  double nMn_conv_gauss(const double& t, const double& xd, const double& m, const double& s);
  double nAp_conv_gauss(const double& t, const double& xd, const double& m, const double& s);
  double nAn_conv_gauss(const double& t, const double& xd, const double& m, const double& s);
  double Mp_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);
  double Mn_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);
  double Mf_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);
  double Ap_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);
  double An_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);
  double Af_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s);

  double _IM(const double& x, const double& m, const double& s, const double& beta, const double& gamma, const double& b);
  double _IA(const double& x, const double& m, const double& s, const double& beta, const double& gamma, const double& b);

  double gaussian(const double& x, const double& m, const double& s);
  double norm_gaussian_w_cutoff(const double& cutoff, const double& m, const double& s);
  double norm_gaussian(const double& ll, const double& ul, const double& m, const double& s);
  double DiracDelta(const double& x);

  double xXi_conv_gauss_by_int(nXi_t p_func, const double& t, const double& xd, const double& mu, const double& sigma);
  double xEp_conv_gauss(const double& t, const double& tau, const double& m, const double& s);
  double xEn_conv_gauss(const double& t, const double& tau, const double& m, const double& s);
  double xEf_conv_gauss(const double& t, const double& tau, const double& m, const double& s);

  double norm_neg_sup(const double& m, const double& s);
  double norm_nEp_conv_gauss_sub(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.);
  double norm_nEn_conv_gauss_sub(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.);
  double norm_nEp_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.);
  double norm_nEn_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.);
  double norm_Ep_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);
  double norm_En_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);
  double norm_nEf_conv_gauss(const double& _ll, const double& _ul,const double& m, const double& s, const double& o = 0.);
  double norm_Ef_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);
  double norm_Enp_conv_gauss(const double& _ll, const double& _ul,const double& tau_n, const double& tau_p,const double& m, const double& s, const double& o = 0.);
  double norm_nag_sup(const double& m, const double& s, const double& xd);
  double norm_nmg_sup(const double& m, const double& s, const double& xd);
  double norm_nAn_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_nAp_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_nAf_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_nMn_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_nMp_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_nMf_conv_gauss(const double& ll, const double& ul, const double& xd, const double& m, const double& s, const double& o = 0.);
  double norm_An_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);
  double norm_Ap_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);
  double norm_Af_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);
  double norm_Mn_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);
  double norm_Mp_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);
  double norm_Mf_conv_gauss(const double& ll, const double& ul,  const double& tau, const double& dm, const double& m, const double& s, const double& o = 0.);

  double int_polyexp2(const double& _ll, const double& _ul,const double& alpha, const double& beta, const double& gamma,const double& a);
  double int_polyexp_erfc(const double& _ll, const double& _ul,const double& alpha, const double& beta, const double& gamma,const double& a);
  double norm_xEp_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);
  double norm_xEn_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);
  double norm_xEf_conv_gauss(const double& _ll, const double& _ul, const double& tau,const double& m, const double& s, const double& o = 0.);

  double norm_xEp(const double& _ll, const double& _ul, const double& tau, const double& o = 0.);
  double norm_xEn(const double& _ll, const double& _ul, const double& tau, const double& o = 0.);
  double norm_xEf(const double& _ll, const double& _ul, const double& tau, const double& o = 0.);

protected:
  // ** Interface parameters for inherits classes **
  double ll, ul;
  double ak,ck;
  double beta;
  double mbzero;
  double K,Kb,C,S;
  double m_tau,m_dm;
  double A,B;
  int flv, xi;
  double sin2beta;
  double cos2beta;

  // ** Constants **
  double inv_sqrt2;
  double inv_sqrt_pi;
  double sqrt_pi;
  double inv_sqrt_2pi;
  double sqrt2;
  RooMath* roomath;
};

#endif // CNVL_H
