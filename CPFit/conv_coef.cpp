#include "conv_coef.h"

void conv_coef::add_EpEn_coef(double& fEp1, double& fEn2,
                          const double& tau1, const double& tau2,
                          const double& weight){
  const double inv_tausum = 1.0/(tau1+tau2);
  fEp1 += weight*tau1*inv_tausum;
  fEn2 += weight*tau2*inv_tausum;
  return;
}

void conv_coef::add_EnEp_coef(double& fEn1, double& fEp2,
                          const double& tau1, const double& tau2,
                          const double& weight){
  const double inv_tausum = 1.0/(tau1+tau2);
  fEn1 += weight*tau1*inv_tausum;
  fEp2 += weight*tau2*inv_tausum;
  return;
}

void conv_coef::add_EnEn_coef(double& fEn1, double& fEn2, double& fxEn1,
                          const double& tau1, const double& tau2,
                          const double& weight){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */){
    fxEn1 += weight;
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    fEn1 +=  weight*tau1*inv_tausub;
    fEn2 += -weight*tau2*inv_tausub;
  }
  return;
}

void conv_coef::add_EpEp_coef(double& fEp1, double& fEp2, double& fxEp1,
                          const double& tau1, const double& tau2,
                          const double& weight){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */){
    fxEp1 += weight;
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    fEp1 +=  weight*tau1*inv_tausub;
    fEp2 += -weight*tau2*inv_tausub;
  }
  return;
}

void conv_coef::add_ApEp_coef(double& fAp1, double& fMp1, double& fEp2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fAp1 += -weight*A*inv_at2*inv_tau;
  fMp1 += -weight*A*inv_at2*dm;
  fEp2 +=  weight*A*dm;
  return;
}

void conv_coef::add_AnEn_coef(double& fAn1, double& fMn1, double& fEn2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fAn1 += -weight*A*inv_at2*inv_tau;
  fMn1 += +weight*A*inv_at2*dm;
  fEn2 += -weight*A*dm;
  return;
}

void conv_coef::add_ApEn_coef(double& fAp1, double& fMp1, double& fEn2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fAp1 += weight*A*inv_at2*inv_tau;
  fMp1 += weight*A*inv_at2*dm;
  fEn2 += weight*A*dm;
  return;
}

void conv_coef::add_AnEp_coef(double& fAn1, double& fMn1, double& fEp2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fAn1 +=  weight*A*inv_at2*inv_tau;
  fMn1 += -weight*A*inv_at2*dm;
  fEp2 += -weight*A*dm;
  return;
}

void conv_coef::add_MpEp_coef(double& fMp1, double& fAp1, double& fEp2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fMp1 += -weight*A*inv_at2*inv_tau;
  fAp1 +=  weight*A*inv_at2*dm;
  fEp2 +=  weight*A*inv_tau;
  return;
}

void conv_coef::add_MnEn_coef(double& fMn1, double& fAn1, double& fEn2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fMn1 += -weight*A*inv_at2*inv_tau;
  fAn1 += -weight*A*inv_at2*dm;
  fEn2 +=  weight*A*inv_tau;
  return;
}

void conv_coef::add_MpEn_coef(double& fMp1, double& fAp1, double& fEn2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fMp1 +=  weight*A*inv_at2*inv_tau;
  fAp1 += -weight*A*inv_at2*dm;
  fEn2 +=  weight*A*inv_tau;
  return;
}

void conv_coef::add_MnEp_coef(double& fMn1, double& fAn1, double& fEp2,
                          const double& tau1, const double& dm,
                          const double& tau2, const double& weight){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  fMn1 +=  weight*A*inv_at2*inv_tau;
  fAn1 +=  weight*A*inv_at2*dm;
  fEp2 +=  weight*A*inv_tau;
  return;
}
