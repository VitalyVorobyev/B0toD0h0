#include "Rk.h"

#ifndef RK_CPP
#define RK_CPP

double RkPdf::norm_EfRk(const double& tau){
  const double r_ckak = ck/ak;
  return 0.5*((1-r_ckak)*norm_En(ll,0,tau*(ak-ck)) + (1+r_ckak)*norm_Ep(0,ul,tau*(ak+ck)));
}

double RkPdf::norm_AfRk(const double& tau, const double& dm){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ntau_n = tau*(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_p = tau*(ak+ck);
  return inv_ak*fact*(norm_An(ll,0,ntau_n,ndm_n) - ndmtau*norm_Mn(ll,0,ntau_n,ndm_n) + norm_Ap(0,ul,ntau_p,ndm_p) - ndmtau*norm_Mp(0,ul,ntau_p,ndm_p));
}

double RkPdf::norm_MfRk(const double& tau, const double& dm){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ntau_n = tau*(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_p = tau*(ak+ck);
  return inv_ak*fact*(norm_Mn(ll,0,ntau_n,ndm_n) + ndmtau*norm_An(ll,0,ntau_n,ndm_n) + norm_Mp(0,ul,ntau_p,ndm_p) + ndmtau*norm_Ap(0,ul,ntau_p,ndm_p));
}

double RkPdf::EfRk(const double& x, const double& tau){
  const double r_ckak = ck/ak;
  return (x<0.0 ? 0.5*(1-r_ckak)*En(x,tau*(ak-ck)) : 0.5*(1+r_ckak)*Ep(x, tau*(ak+ck)));
}

double RkPdf::AfRk(const double& x, const double& tau, const double& dm){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;
  if(x<0.0){
    const double ndm  = dm/(ak-ck);
    const double ntau = tau*(ak-ck);
    Li = inv_ak*fact*(An(x,ntau,ndm) - ndmtau*Mn(x,ntau,ndm));
  }else{
    const double ndm  = dm/(ak+ck);
    const double ntau = tau*(ak+ck);
    Li = inv_ak*fact*(Ap(x,ntau,ndm) - ndmtau*Mp(x,ntau,ndm));
  }
  return Li;
}

double RkPdf::MfRk(const double& x, const double& tau, const double& dm){
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;
  if(x<0.0){
    const double ndm  = dm/(ak-ck);
    const double ntau = tau*(ak-ck);
    Li = inv_ak*fact*(Mn(x,ntau,ndm)+ndmtau*An(x,ntau,ndm));
  }else{
    const double ndm  = dm/(ak+ck);
    const double ntau = tau*(ak+ck);
    Li = inv_ak*fact*(Mp(x,ntau,ndm)+ndmtau*Ap(x,ntau,ndm));
  }
  return Li;
}

double RkPdf::PdfAB(const double& dt, const bool no_interf){
  double pdf      = EfRk(dt,m_tau);
  double pdf_norm = norm_EfRk(m_tau);
  if(!no_interf){
    pdf      += - 0.5/m_tau*A*MfRk(dt,m_tau,m_dm)   + 0.5/m_tau*B*AfRk(dt,m_tau,m_dm);
    pdf_norm += - 0.5/m_tau*A*norm_MfRk(m_tau,m_dm) + 0.5/m_tau*B*norm_AfRk(m_tau,m_dm);
  }
  if(pdf<0 || pdf_norm<=0){
    cout << "PdfAB: pdf = " << pdf << ", norm = " << pdf_norm << endl;
    return 0;
  }
  return pdf/pdf_norm;
}

double RkPdf::Pdf(const double& dt, const bool no_interf){
  double pdf = (K+Kb)*EfRk(dt,m_tau);
  double pdf_norm = (K+Kb)*norm_EfRk(m_tau);
  if(!no_interf){
    pdf      += - 0.5/m_tau*flv*(K-Kb)*MfRk(dt,m_tau,m_dm) + 0.5/m_tau*2.*flv*xi*sqrt(K*Kb)*AfRk(dt,m_tau,m_dm)*(C*sin2beta+S*cos2beta);
    pdf_norm += - 0.5/m_tau*flv*(K-Kb)*norm_MfRk(m_tau,m_dm);// + 2.*flv*xi*sqrt(K*Kb)*norm_AfRk(m_tau,m_dm)*(C*sin2beta+S*cos2beta);
  }
  if(pdf<0 || pdf_norm<=0){
    cout << "Pdf: pdf = " << pdf << ", norm = " << pdf_norm << endl;
    return 0;
  }
  return pdf/pdf_norm;
}

#endif // RK_CPP
