#ifndef PHASESPACE_H
#define PHASESPACE_H

#include "TH2F.h"
#include "TFile.h"

#include <math.h>
#include <iomanip>
using namespace std;

class PhaseSpace{
public:
  PhaseSpace(void);
  bool IsInPlot(const double& mp,const double& mm);
  inline double Mpp(const double& mp,const double& mm);

  int GetBin(const double&, const double&);

  double m_pi() const {return M_pi;}
  double m_Ks() const {return M_Ks;}
  double m_D() const {return M_D;}
  double m_min() const {return M_min;}
  double m_max() const {return M_max;}
  double mm_max(const double& mp);
  double mm_min(const double& mp);
private:
  double M_D,M_pi,M_Ks;
  double M_D_sq,M_pi_sq,M_Ks_sq;
  double M_min,M_max;

  inline double E_Ks_star(const double& mp);
  inline double E_pi_star(const double& mp);

  TH2F *dkpp_bin_h;
};

#endif // PHASESPACE_H
