#ifndef RK_H
#define RK_H

#include "cnvl.h"

class RkPdf: public cnvl{
public:
  RkPdf(): cnvl() {}
//  ~RkPdf(){}
  double PdfAB(const double& dt,const bool no_interf = false);
  double Pdf(const double& dt,const bool no_interf = false);

private:
  double norm_EfRk(const double& tau);
  double norm_AfRk(const double& tau, const double& dm);
  double norm_MfRk(const double& tau, const double& dm);
  double EfRk(const double& x, const double& tau);
  double AfRk(const double& x, const double& tau, const double& dm);
  double MfRk(const double& x, const double& tau, const double& dm);

//  double ll,ul;
//  double ak,ck;
//  double beta;
//  double mbzero;
//  double K,Kb,C,S;
//  double m_tau,m_dm;
//  double A,B;
//  int flv, xi;
//  double sin2beta;
//  double cos2beta;
//  double mbplus;
};

#endif // RK_H
