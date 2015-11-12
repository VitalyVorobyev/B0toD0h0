#ifndef CONV_COEF_H
#define CONV_COEF_H

#include <iostream>
#include <cmath>

using namespace std;

class conv_coef{
protected:
  void add_EpEn_coef(double& fEp1, double& fEn2,                const double& tau1, const double& tau2, const double& weight = 1);
  void add_EnEp_coef(double& fEn1, double& fEp2,                const double& tau1, const double& tau2, const double& weight = 1);
  void add_EnEn_coef(double& fEn1, double& fEn2, double& fxEn1, const double& tau1, const double& tau2, const double& weight = 1);
  void add_EpEp_coef(double& fEp1, double& fEp2, double& fxEp1, const double& tau1, const double& tau2, const double& weight = 1);

  void add_ApEp_coef(double& fAp1, double& fMp1, double& fEp2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_AnEn_coef(double& fAn1, double& fMn1, double& fEn2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_ApEn_coef(double& fAp1, double& fMp1, double& fEn2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_AnEp_coef(double& fAn1, double& fMn1, double& fEp2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);

  void add_MpEp_coef(double& fMp1, double& fAp1, double& fEp2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_MnEn_coef(double& fMn1, double& fAn1, double& fEn2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_MpEn_coef(double& fMp1, double& fAp1, double& fEn2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
  void add_MnEp_coef(double& fMn1, double& fAn1, double& fEp2, const double& tau1, const double& dm, const double& tau2, const double& weight = 1.0);
};

#endif // CONV_COEF_H
