#include "cnvl.h"

int cnvl::SetRange(const double& dtmax){
  ul = dtmax; ll = -ul;
  return 0;
}

int cnvl::SetAkCk(const double& costh, const double& ecm){
  ak = ecm/mbzero;
  const double pcm = sqrt(ecm*ecm-mbzero*mbzero);
  ck = pcm*costh/(beta*mbzero);
  if(ak==0.0||ak==ck){
    cout << "SetAkCk: Invalid ak = " << ak << ", ck = " << ck << ", where they should be ak!=0.0&&ak!=ck." << endl;
    return -1;
  }// else{
//    cout << "SetAkCk: ak = " << ak << ", ck = " << ck << endl;
//  }
  return 0;
}

int cnvl::SetTauDm(const double& _tau, const double& _dm){
  m_tau = _tau; m_dm = _dm;
  return 0;
}

int cnvl::SetKKCS(const double& _K, const double& _Kb, const double& _C, const double& _S){
  K = _K; Kb = _Kb; C = _C; S = _S;
  return 0;
}
int cnvl::SetAB(const double& _A, const double& _B){
  A = _A; B = _B;
  return 0;
}
int cnvl::SetFlvXi(const int _flv, const int _xi){
  flv = _flv, xi = _xi;
  return 0;
}
int cnvl::SetXi(const int _xi){
  xi = _xi;
  return 0;
}
int cnvl::SetSinCos(const double& sin, const double& cos){
  sin2beta = sin; cos2beta = cos;
  return 0;
}

double cnvl::xXi_conv_gauss_by_int(nXi_t p_func,
                                    const double& t, const double& xd,
                                    const double& mu, const double& sigma){
  double area=0.0;
  if(sigma==0.0) return (*p_func)(t,xd);
  /* integral (-20sigma -- +20 sigma)*/
  const double t_ini = -20.0*sigma+mu;
  const double t_end = +20.0*sigma+mu;
  const int ndiv = 400;
  const double t_del = (t_end-t_ini)/ndiv;

  for(int i=0;i<ndiv;++i){
    double t_run = t_ini+t_del*i;
    const double err = gaussian(t_run, t-mu, sigma);
    area += err*(*p_func)(t_run, xd);
  }
  area -= 0.5*gaussian(t_ini, t-mu, sigma)*(*p_func)(t_ini, xd);
  area -= 0.5*gaussian(t_end, t-mu, sigma)*(*p_func)(t_end, xd);
  area *= t_del;
  return area;
}

double cnvl::recexp(const double& re, const double& im){
  complex<double> z(re,im);
  return exp(z).real();
}

double cnvl::imcexp(const double& re, const double& im){
  complex<double> z(re,im);
  return exp(z).imag();
}

double cnvl::rewerf(const double& re, const double& im){
  complex<double> z(re,im);
//  return roomath->erf(z).real();
  return roomath->faddeeva(z).real();
}

double cnvl::imwerf(const double& re, const double& im){
  complex<double> z(re,im);
//  return roomath->erf(z).imag();
  return roomath->faddeeva(z).imag();
}

double cnvl::Ep(const double& t, const double& tau){
  if(t<0.0) return 0.0;
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  return inv_atau*exp(-at*inv_atau);
}

double cnvl::En(const double& t, const double& tau){
  if(t>=0.0) return 0.0;
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  return inv_atau*exp(-at*inv_atau);
}

double cnvl::Ef(const double& t, const double& tau){
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  return 0.5*inv_atau*exp(-at*inv_atau);
}

double cnvl::Enp(const double& t, const double& tau_n, const double& tau_p){
  const double at   = fabs(t);
  const double norm = 1.0/(fabs(tau_n)+fabs(tau_p));
  double inv_atau = (t>=0) ? 1.0/fabs(tau_p) : 1.0/fabs(tau_n) ;
  return norm*exp(-at*inv_atau);
}

double cnvl::xEp(const double& t, const double& tau){
  const double inv_atau = 1.0/fabs(tau);
  return t*inv_atau*Ep(t, tau);
}

double cnvl::xEn(const double& t, const double& tau){
  const double inv_atau = 1.0/fabs(tau);
  return -t*inv_atau*En(t, tau);
}

double cnvl::xEf(const double& t, const double& tau){
  const double inv_atau = 1.0/fabs(tau);
  return fabs(t)*inv_atau*Ef(t, tau);
}

double cnvl::nMp(const double& t, const double& xd){
  const double at = fabs(t);
  if(t<0.0) return 0.0;
  return exp(-at)*cos(xd*t);
}

double cnvl::Mp(const double& t, const double& tau, const double& dm){
  if (t<0.0) return 0.0;
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMp(nt,xd);
}

double cnvl::nMn(const double& t, const double& xd){
  if (t>=0.0) return 0.0;
  const double at = fabs(t);
  return exp(-at)*cos(xd*t);
}

double cnvl::Mn(const double& t, const double& tau, const double& dm){
  if (t>=0.0) return 0.0;
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMn(nt,xd);
}

double cnvl::nMf(const double& t, const double& xd){
  const double at = fabs(t);
  return exp(-at)*cos(xd*t);
}

double cnvl::Mf(const double& t, const double& tau, const double& dm){
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMf(nt, xd);
}

double cnvl::nAp(const double& t, const double& xd){
  if(t<0.0) return 0.0;
  const double at = fabs(t);
  return exp(-at)*sin(xd*t);
}

double cnvl::Ap(const double& t, const double& tau, const double& dm){
  if(t<0.0) return 0.0;
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAp(nt, xd);
}

double cnvl::nAn(const double& t, const double& xd){
  if(t>=0.0) return 0.0;
  const double at = fabs(t);
  return exp(-at)*sin(xd*t);
}

double cnvl::An(const double& t, const double& tau, const double& dm){
  if(t>=0.0) return 0.0;
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAn(nt, xd);
}

double cnvl::nAf(const double& t, const double& xd){
  const double at   = fabs(t);
  return exp(-at)*sin(xd*t);
}

double cnvl::Af(const double& t, const double& tau, const double& dm){
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAf(nt, xd);
}

double cnvl::norm_Ap(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o){
  double nll = (_ll<_ul) ? _ll : _ul;
  double nul = (_ll<_ul) ? _ul : _ll;
  nll -= o;
  nul -= o;
  if(nul<=0.0) return 0.0;
  double f = norm_Ax_sup(nul, tau, dm);
  f -= norm_Ax_sup((nll>0.0) ? nll : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double cnvl::norm_Ax_sup(const double& x, const double& tau, const double& dm){
  const double dmt   = dm*x;
  const double dmtau = dm*tau;
  double f = -sin(dmt)-dmtau*cos(dmt);
  f *= tau*exp(-x/tau)/(1+dmtau*dmtau);
  return f;
}

double cnvl::norm_Mp(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o){
  double nll = (_ll<_ul) ? _ll : _ul;
  double nul = (_ll<_ul) ? _ul : _ll;
  nll -= o;
  nul -= o;
  if(nul <= 0.0) return 0.0;
  double f = norm_Mx_sup(nul, tau, dm);
  f -= norm_Mx_sup((nll>0.0) ? nll : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double cnvl::norm_Mn(const double& _ll, const double& _ul, const double& tau, const double& dm, const double&o){
  double nll = (_ll<_ul) ? _ll : _ul;
  double nul = (_ll<_ul) ? ul : _ll;
  nll -= o;
  nul -= o;
  if(nll>0.0) return 0.0;
  double f = norm_Mx_sup(-nll, tau, dm);
  f -= norm_Mx_sup((nul<0.0) ? -nul : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double cnvl::norm_Mx_sup(const double& x, const double& tau, const double& dm){
  const double dmt   = dm*x;
  const double dmtau = dm*tau;
  double f = -cos(dmt)-dmtau*sin(dmt);
  f *= tau*exp(-x/tau)/(1+dmtau*dmtau);
  return f;
}

double cnvl::norm_An(const double& _ll, const double& _ul, const double& tau, const double& dm, const double& o){
  double nll = (_ll<_ul) ? _ll : _ul;
  double nul = (_ll<_ul) ? _ul : _ll;
  nll -= o;
  nul -= o;
  if(nll>0.0) return 0.0;
  double f = norm_Ax_sup((nul<0.0) ? -nul : 0.0, tau, dm);
  f -= norm_Ax_sup(-nll, tau, dm);
  return (ll<ul) ? f : -f;
}

double cnvl::norm_nEp(const double& _ll, const double& _ul, const double& o){
  const double nul = (_ul-o>=0.0) ? (_ul-o) : 0.0;
  const double nll = (_ll-o>=0.0) ? (_ll-o) : 0.0;
  if(nul==nll) return 0;
  double f = exp(-nll)-exp(-nul);
  return f;
}

double cnvl::norm_nEn(const double& _ll, const double& _ul, const double& o){
  const double nul = (_ul-o<0.0) ? (_ul-o) : 0.0;
  const double nll = (_ll-o<0.0) ? (_ll-o) : 0.0;
  if(nul==nll) return 0;
  double f = exp(nul)-exp(nll);
  return f;
}

double cnvl::norm_nEf(const double& _ll, const double& _ul, const double& o){
  return 0.5*(norm_nEn(_ll,_ul,o)+norm_nEp(_ll,_ul,o));
}

double cnvl::norm_Ep(const double& _ll, const double& _ul, const double& tau, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEp(nll,nul,no);
}

double cnvl::norm_En(const double& _ll, const double& _ul, const double& tau, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEn(nll,nul,no);
}

double cnvl::norm_Ef(const double& _ll, const double& _ul, const double& tau, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEf(nll,nul,no);
}

double cnvl::nEp_conv_gauss(const double& t, const double& m, const double& s){
  static const double Tc = DBL_MAX_10_EXP*log(10.0);
  if(s == 0.0) return Ep(t-m, 1.0);
  double inv_s = 1.0/fabs(s);
  double dt = -t+m;
  double ex = 0.5*s*s+dt;
  if(ex<Tc) return 0.5*exp(ex)*erfc(inv_sqrt2*inv_s*(s*s+dt));
  const double gf = exp(-0.5*dt*dt*inv_s*inv_s);
  const double x  = inv_sqrt2*(s+dt*inv_s);
  return 0.5*approx_exp2erfc(x)*gf;
}

double cnvl::nEn_conv_gauss(const double& t, const double& m, const double& s){
  if(s==0.0) return En(t-m, 1.0);
  static const double Tc = DBL_MAX_10_EXP*log(10.0);
  double inv_s = 1.0/fabs(s);
  double dt = t-m;
  double ex = 0.5*s*s+dt;
  if(ex<Tc) return 0.5*exp(ex)*erfc(inv_sqrt2*inv_s*(s*s+dt));
  const double gf = exp(-0.5*dt*dt*inv_s*inv_s);
  const double x  = inv_sqrt2*(s+dt*inv_s);
  return 0.5*approx_exp2erfc(x)*gf;
}

double cnvl::Ep_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double f = inv_atau*nEp_conv_gauss(nt, nm, ns);
  return f;
}

double cnvl::En_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double f = inv_atau*nEn_conv_gauss(nt, nm, ns);
  return f;
}

double cnvl::Ef_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  double f = nEn_conv_gauss(nt, nm, ns)+nEp_conv_gauss(nt, nm, ns);
  f *= 0.5*inv_atau;
  return f;
}

double cnvl::Enp_conv_gauss(const double& t, const double& tau_n, const double& tau_p, const double& m, const double& s){
  const double inv_atau_n = 1.0/fabs(tau_n);
  const double nt_n = t*inv_atau_n;
  const double nm_n = m*inv_atau_n;
  const double ns_n = fabs(s)*inv_atau_n;

  const double inv_atau_p = 1.0/fabs(tau_p);
  const double nt_p = t*inv_atau_p;
  const double nm_p = m*inv_atau_p;
  const double ns_p = fabs(s)*inv_atau_p;

  double f = nEn_conv_gauss(nt_n, nm_n, ns_n)+nEp_conv_gauss(nt_p, nm_p, ns_p);
  f *= 1.0/(fabs(tau_n)+fabs(tau_p));
  return f;
}

double cnvl::approx_exp2erfc(const double& x){
  const double inv_x = 1./x;
  return inv_sqrt_pi*(inv_x-0.5*inv_x*inv_x*inv_x);
}

double cnvl::nMp_conv_gauss(const double& t, const double& xd, const double& m, const double& s){
  if(!finite(s)) return 0.0;
  if(s==0.0)     return Mp(t-m, 1.0, xd);
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)) return Mp(t-m, 1.0, xd);
  const double inv_s2 = 1.0/(s*s);
  double f = _IM(t, m, s, 0.5*inv_s2, 1-(t-m)*inv_s2, xd);
  if(!finite(f)) return xXi_conv_gauss_by_int(nMp,t,xd,m,s);
  return f;
}

double cnvl::nMn_conv_gauss(const double& t, const double& xd, const double& m, const double& s){
  if(!finite(s)) return 0.0;
  if(s==0.0)     return Mn(t-m, 1.0, xd);
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)) return Mn(t-m, 1.0, xd);
  const double inv_s2 = 1.0/(s*s);
  double f = _IM(t, m, s, 0.5*inv_s2, 1+(t-m)*inv_s2, xd);
  if(!finite(f)) return xXi_conv_gauss_by_int(nMn,t,xd,m,s);
  return f;
}

double cnvl::nAp_conv_gauss(const double& t, const double& xd, const double& m, const double& s){
  if(!finite(s)) return 0.0;
  if(s==0.0) return Ap(t-m, 1.0, xd);
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)) return Ap(t-m, 1.0, xd);
  const double inv_s2 = 1.0/s/s;
  double f = _IA(t, m, s, 0.5*inv_s2, 1-(t-m)*inv_s2, xd);
  if(!finite(f)) return xXi_conv_gauss_by_int(nAp,t,xd,m,s);
  return f;
}

double cnvl::nAn_conv_gauss(const double& t, const double& xd, const double& m, const double& s){
  if(!finite(s)) return 0.0;
  if(s==0.0) return An(t-m, 1.0, xd);
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)) return An(t-m, 1.0, xd);
  const double inv_s2 = 1.0/s/s;
  double f = -1.0*_IA(t, m, s, 0.5*inv_s2, 1+(t-m)*inv_s2, xd);
  if(!finite(f)) return xXi_conv_gauss_by_int(nAn, t, xd, m, s);
  return f;
}

double cnvl::Mp_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMp_conv_gauss(nt, xd, nm,ns);
  return f;
}

double cnvl::Mn_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMn_conv_gauss(nt,xd,nm,ns);
  return f;
}

double cnvl::Mf_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMn_conv_gauss(nt,xd,nm,ns)+nMp_conv_gauss(nt,xd,nm,ns);
//    f *= 0.5;
  return f;
}

double cnvl::Ap_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAp_conv_gauss(nt,xd,nm,ns);
  return f;
}

double cnvl::An_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAn_conv_gauss(nt,xd,nm,ns);
  return f;
}

double cnvl::Af_conv_gauss(const double& t, const double& tau, const double& dm, const double& m, const double& s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAn_conv_gauss(nt,xd,nm,ns)+nAp_conv_gauss(nt,xd,nm,ns);
//    f *= 0.5;
  return f;
}

double cnvl::_IM(const double& x, const double& m, const double& s,
                  const double& beta, const double& gamma, const double& b){
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA==0)
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = rewerf(inv_2sqrtbeta*b, inv_2sqrtbeta*gamma);
  f *= gaussian(x, m, s)*sqrt_pi*inv_2sqrtbeta;
  return f;
#else /* USE_CERNLIB_WERF */
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = rewerf(inv_2sqrtbeta*b, inv_2sqrtbeta*fabs(gamma));
  f *= gaussian(x,m,s)*sqrt_pi*inv_2sqrtbeta;
  if(/* gamma */ inv_2sqrtbeta*gamma >=0){
    return f;
  }
  f *= -1.0;
  const double inv_beta = 1.0/beta;
  const double dx = x-m;
  const double inv_s = 1.0/fabs(s);
  const double rex = 0.25*(gamma*gamma-b*b)*inv_beta-0.5*dx*dx*inv_s*inv_s;
  double imx = 0.5*fabs(gamma)*b*inv_beta;
  f += inv_sqrt_2pi*inv_s*sqrt(M_PI*inv_beta)*recexp(rex, imx);
  return f;
#endif
}

double cnvl::_IA(const double& x, const double& m, const double& s,
                  const double& beta, const double& gamma, const double& b){
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA==0)
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = imwerf(inv_2sqrtbeta*b, inv_2sqrtbeta*gamma);
  f *= gaussian(x,m,s)*sqrt_pi*inv_2sqrtbeta;
  return f;
#else /* USE_CERNLIB_WERF */
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = imwerf(inv_2sqrtbeta*b, inv_2sqrtbeta*fabs(gamma));
  f *= gaussian(x,m,s)*sqrt_pi*inv_2sqrtbeta;
  if( /*  gamma */  inv_2sqrtbeta*gamma  >= 0){
    return f;
  }
  const double inv_beta = 1.0/beta;
  const double dx = x-m;
  const double inv_s = 1.0/fabs(s);
  const double rex = 0.25*(gamma*gamma-b*b)*inv_beta-0.5*dx*dx*inv_s*inv_s;
  const double imx = 0.5*fabs(gamma)*b*inv_beta;
  f += inv_sqrt_2pi*inv_s*sqrt(M_PI*inv_beta)*imcexp(rex,imx);
  return f;
#endif /* USE_CERNLIB_WERF */
}

double cnvl::gaussian(const double& x, const double& m, const double& s){
  if(!finite(s)) return 0;
  if(s==0.0) return DiracDelta(x-m);
  double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)) return DiracDelta(x-m);
  const double dx = x-m;
  return inv_sqrt_2pi*inv_s*exp(-0.5*dx*dx*inv_s*inv_s);
}

double cnvl::norm_gaussian_w_cutoff(const double& cutoff, const double& m, const double& s){
  double a_s = fabs(s);
  if(s==0.0) return 1;
  const double inv_sqrt2 = 0.707106781;
  double inv_s = 1.0/a_s;
  double x1 = (cutoff+m)*inv_sqrt2*inv_s;
  double x2 = (cutoff-m)*inv_sqrt2*inv_s;
  return 1.0-0.5*erfc(x1)-0.5*erfc(x2);
}

double cnvl::norm_gaussian(const double& ll, const double& ul, const double& m, const double& s){
  double a_s = fabs(s);
  if(s==0.0) return 1;
  const double inv_sqrt2 = 0.707106781;
  double inv_s = 1.0/a_s;
  double x1 = (-ll+m)*inv_sqrt2*inv_s;
  double x2 = (ul-m)*inv_sqrt2*inv_s;
  return 1.0-0.5*erfc(x1)-0.5*erfc(x2);
}

double cnvl::DiracDelta(const double& x){
  if(x==0.0) return FLT_MAX;
  return 0.0;
}

double cnvl::xEn_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  if(s==0.0) return xEn(t-m, tau);
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  return (-ns2-nt+nm)*En_conv_gauss(t, tau, m, s) + ns2*gaussian(t, m, s);
}

double cnvl::xEf_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  if(s==0.0) return xEn(t-m, tau);
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  return 0.5*(-ns2+nt-nm)*Ep_conv_gauss(t, tau, m, s) +
         0.5*(-ns2-nt+nm)*En_conv_gauss(t, tau, m, s) +
         ns2*gaussian(t, m, s);
}

double cnvl::xEp_conv_gauss(const double& t, const double& tau, const double& m, const double& s){
  if(s==0.0) return xEp(t-m, tau);
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  return (-ns2+nt-nm)*Ep_conv_gauss(t,tau,m,s) + ns2*gaussian(t,m,s);
}

double cnvl::norm_neg_sup(const double& m, const double& s){
  const double Tc = DBL_MAX_10_EXP*log(10.0);

  const double as = fabs(s);
  const double ex = 0.5*as*as-m;
  const double x = inv_sqrt2*(as-m/as);
  if(ex<Tc){
    double f = exp(ex);
    double g = erfc(x);
    if(g>0.0&&finite(f)){
      return f*g;
    }
  }
  double f = exp(-0.5*m*m/(s*s));
  double g = approx_exp2erfc(x);
  return f*g;
}

double cnvl::norm_nEp_conv_gauss_sub(const double& _ll, const double& _ul,
                               const double& m, const double& s, const double& o){
  const double f1 = erfc(inv_sqrt2*(-ll+o+m)/s) - norm_neg_sup(ll-o-m,s);
  const double f2 = erfc(inv_sqrt2*(ul-o-m)/s)  + norm_neg_sup(ul-o-m,s);
  return 0.5*(f1+f2);
}

double cnvl::norm_nEn_conv_gauss_sub(const double& _ll, const double& _ul,
                                      const double& m, const double& s, const double& o){
  const double f1 = erfc(inv_sqrt2*(-_ll+o+m)/s)  + norm_neg_sup(-_ll+o+m,s);
  const double f2 = erfc(inv_sqrt2*(_ul-o-m)/s) - norm_neg_sup(-_ul+o+m,s);
  return 0.5*(f1+f2);
}

double cnvl::norm_nEp_conv_gauss(const double& _ll, const double& _ul,
                           const double& m, const double& s, const double& o){
  return 1.0-norm_nEp_conv_gauss_sub(_ll,_ul,m,s,o);
}

double cnvl::norm_nEn_conv_gauss(const double& _ll, const double& _ul,
                           const double& m, const double& s, const double& o){
  return 1.0-norm_nEn_conv_gauss_sub(_ll,_ul,m,s,o);
}

double cnvl::norm_Ep_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  return 1.0-norm_nEp_conv_gauss_sub(nll, nul, nm, ns, no);
}

double cnvl::norm_En_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  return 1.0-norm_nEn_conv_gauss_sub(nll, nul, nm, ns, no);
}

double cnvl::norm_nEf_conv_gauss(const double& _ll, const double& _ul,
                           const double& m, const double& s, const double& o){
  return 0.5*(norm_nEn_conv_gauss(_ll,_ul,m,s,o)+norm_nEp_conv_gauss(_ll,_ul,m,s,o));
}

double cnvl::norm_Ef_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = _ll*inv_atau;
  const double nul = _ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = (o==0.0) ? 0.0 : o*inv_atau;
  return 1.0-0.5*(norm_nEp_conv_gauss_sub(nll, nul, nm, ns, no)
                 +norm_nEn_conv_gauss_sub(nll, nul, nm, ns, no));
}

double cnvl::norm_Enp_conv_gauss(const double& _ll, const double& _ul,
                           const double& tau_n, const double& tau_p,
                           const double& m, const double& s, const double& o){
  const double inv_atau_n = 1.0/fabs(tau_n);
  const double nll_n = _ll*inv_atau_n;
  const double nul_n = _ul*inv_atau_n;
  const double nm_n  = m*inv_atau_n;
  const double ns_n  = s*inv_atau_n;
  const double no_n  = o*inv_atau_n;

  const double inv_atau_p = 1.0/fabs(tau_p);
  const double nll_p = _ll*inv_atau_p;
  const double nul_p = _ul*inv_atau_p;
  const double nm_p  = m*inv_atau_p;
  const double ns_p  = s*inv_atau_p;
  const double no_p  = o*inv_atau_p;

  const double fn = fabs(tau_n)/(fabs(tau_n)+fabs(tau_p));
  const double fp = fabs(tau_p)/(fabs(tau_n)+fabs(tau_p));
  double f = fn*norm_nEn_conv_gauss_sub(nll_n, nul_n, nm_n, ns_n, no_n);
  f       += fp*norm_nEp_conv_gauss_sub(nll_p, nul_p, nm_p, ns_p, no_p);
  return 1.0-f;
}

double cnvl::norm_nag_sup(const double& m, const double& s, const double& xd){
  const double as = fabs(s);
  const double inv_as = 1.0/as;

  double rex = inv_sqrt2*xd*as;
  double imx = inv_sqrt2*(as+m*inv_as);
  double f=0;
  if(imx>=0){
    double f = imwerf(rex, imx)+xd*rewerf(rex, imx);
    f *= -exp(-0.5*m*m*inv_as*inv_as);
    return f;
  }

  double rexx = -rex*rex+imx*imx-0.5*m*m*inv_as*inv_as;
  double imxx = -2.0*rex*imx;

  f = exp(-0.5*m*m*inv_as*inv_as)*(imwerf(-rex, -imx)+xd*rewerf(-rex, -imx));
  f += -2.0*(imcexp(rexx, imxx)+xd*recexp(rexx, imxx));
  return f;
}

double cnvl::int_polyexp2(const double& _ll, const double& _ul,
                    const double& alpha, const double& beta, const double& gamma,
                    const double& a){
  const double inv_2a = 0.5/a;
  const double sqrt_a = sqrt(a);

  double f = -inv_2a*(alpha*_ul+beta)*exp(-a*_ul*_ul);
  f       +=  inv_2a*(alpha*_ll+beta)*exp(-a*_ll*_ll);
  f += (alpha*inv_2a+gamma)*sqrt_pi/sqrt_a*(1.0-0.5*erfc(sqrt_a*_ul)-0.5*erfc(sqrt_a*_ll));
  return f;
}

double cnvl::int_polyexp_erfc(const double& _ll, const double& _ul,
                        const double& alpha, const double& beta, const double& gamma,
                        const double& a){
  const double inv_a = 1.0/a;
  double f = 0.0;

  const double a1 = -2.0*alpha*inv_a + beta;
  const double a0 = 2.0*alpha*alpha*inv_a*inv_a + beta*inv_a + a*gamma;

  f +=  inv_a*exp(a*_ul)*(alpha*_ul*_ul+a1*_ul+a0)*erfc(_ul);
  f += -inv_a*exp(a*_ll)*(alpha*_ll*_ll+a1*_ll+a0)*erfc(_ll);

  const double aa1 = a1 + a*alpha;
  const double aa0 = a0 + 0.5*a*a1 + 0.25*a*a*alpha;

  f += exp(0.25*a*a)*int_polyexp2(_ll-0.5*a, _ul-0.5*a, alpha, aa1, aa0, 1.0);
  return f;
}

double cnvl::norm_xEp_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                           const double& m, const double& s, const double& o){
  if(s==0.0) return norm_xEp(_ll-m,_ul-m,tau,o);
  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double nll = (-(_ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = (-(_ul-m-o)-ns*s)*inv_sqrt2*inv_s;

  const double f
    = ns2*exp(1.5*ns2)*int_polyexp_erfc(nll,nul,0.0,1.0,0.0,ns*sqrt2)
    + ns2*norm_gaussian(_ll,_ul,m+o,s);
  return f;
}

double cnvl::norm_xEn_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                           const double& m, const double& s, const double& o){
  if(s==0.0) return norm_xEn(_ll-m,_ul-m,tau,o);
  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double nll = ((_ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = ((_ul-m-o)-ns*s)*inv_sqrt2*inv_s;

  const double f
    = ns2*exp(1.5*ns2)*int_polyexp_erfc(nll, nul, 0.0, 1.0, 0.0, ns*sqrt2)
    + ns2*norm_gaussian(_ll,_ul,m+o,s);
  return f;
}

double cnvl::norm_xEf_conv_gauss(const double& _ll, const double& _ul, const double& tau,
                           const double& m, const double& s, const double& o){
  return 0.5*(norm_xEp_conv_gauss(_ll,_ul, tau, m, s, o)
              +norm_xEn_conv_gauss(_ll,_ul, tau, m, s, o));
}

double cnvl::norm_xEp(const double& _ll, const double& _ul, const double& tau, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (_ul-o>=0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (_ll-o>=0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll) return 0;
  const double vu = exp(-nul)*(-1.0-nul);
  const double vl = exp(-nll)*(-1.0-nll);
  return (vu-vl);
}

double cnvl::norm_xEn(const double& _ll, const double& _ul, const double& tau, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (_ul-o<0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (_ll-o<0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll) return 0;
  const double vu = exp(nul)*(1.0-nul);
  const double vl = exp(nll)*(1.0-nll);
  return (vu-vl);
}

double cnvl::norm_xEf(const double& _ll, const double& _ul, const double& tau, const double& o){
  return 0.5*(norm_xEn(_ll,_ul,tau,o)+norm_xEp(_ll,_ul,tau,o));
}

double cnvl::norm_nmg_sup(const double& m, const double& s, const double& xd){
  const double as = fabs(s);
  const double inv_as = 1.0/as;

  double rex = inv_sqrt2*xd*as;
  double imx = inv_sqrt2*(as+m*inv_as);
  double f=0;
  if(imx>=0){
    f = rewerf(rex, imx)-xd*imwerf(rex, imx);
    f *= -exp(-0.5*m*m*inv_as*inv_as);
    return f;
  }
  double rexx = -rex*rex+imx*imx-0.5*m*m*inv_as*inv_as;
  double imxx = -2.0*rex*imx;

  f = exp(-0.5*m*m*inv_as*inv_as)*(rewerf(-rex, -imx)-xd*imwerf(-rex, -imx));
  f += -2.0*(recexp(rexx, imxx)-xd*imcexp(rexx, imxx));
  return f;
}

double cnvl::norm_nAn_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  const double inv_as = 1.0/fabs(s);

  double f = -xd;
  f += 0.5*(xd*erfc(inv_sqrt2*inv_as*(-ll+o+m))-norm_nag_sup(ll-o-m, s, xd));
  f += 0.5*(xd*erfc(inv_sqrt2*inv_as*(ul-o-m)) +norm_nag_sup(ul-o-m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double cnvl::norm_nAp_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  const double inv_as = 1.0/fabs(s);
  double f = xd;
  f += -0.5*(xd*erfc(inv_sqrt2*inv_as*(-ll+o+m))+norm_nag_sup(-ll+o+m, s, xd));
  f += -0.5*(xd*erfc(inv_sqrt2*inv_as*(ul-o-m)) -norm_nag_sup(-ul+o+m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double cnvl::norm_nAf_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  return norm_nAn_conv_gauss(ll,ul,xd,m,s,o) + norm_nAp_conv_gauss(ll,ul,xd,m,s,o);
}

double cnvl::norm_nMn_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  const double inv_as = 1.0/fabs(s);
  double f = 1.0;
  f += -0.5*(erfc(inv_sqrt2*inv_as*(-ll+o+m))-norm_nmg_sup(ll-o-m, s, xd));
  f += -0.5*(erfc(inv_sqrt2*inv_as*(ul-o-m)) +norm_nmg_sup(ul-o-m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double cnvl::norm_nMp_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  const double inv_as = 1.0/fabs(s);
  double f = 1.0;
  f += -0.5*(erfc(inv_sqrt2*inv_as*(-ll+o+m))+norm_nmg_sup(-ll+o+m, s, xd));
  f += -0.5*(erfc(inv_sqrt2*inv_as*(ul-o-m)) -norm_nmg_sup(-ul+o+m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double cnvl::norm_nMf_conv_gauss(const double& ll, const double& ul, const double& xd,
                           const double& m, const double& s, const double& o){
  return norm_nMn_conv_gauss(ll,ul,xd,m,s,o) + norm_nMp_conv_gauss(ll,ul,xd,m,s,o);
}

double cnvl::norm_An_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAn_conv_gauss(nll, nul, xd, nm, ns, no);
}

double cnvl::norm_Ap_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAp_conv_gauss(nll, nul, xd, nm, ns, no);
}

double cnvl::norm_Af_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAf_conv_gauss(nll, nul, xd, nm, ns, no);
}

double cnvl::norm_Mn_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMn_conv_gauss(nll,nul,xd,nm,ns,no);
}

double cnvl::norm_Mp_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMp_conv_gauss(nll,nul,xd,nm,ns,no);
}

double cnvl::norm_Mf_conv_gauss(const double& ll, const double& ul,
                          const double& tau, const double& dm,
                          const double& m, const double& s, const double& o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMf_conv_gauss(nll,nul,xd,nm,ns,no);
}

