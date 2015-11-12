#ifndef DALITZMODEL_H
#define DALITZMODEL_H

#include "EvtResonance2.h"
#include "phasespace.h"

#include <math.h>
#include <vector>

class DalitzModel{
public:
  DalitzModel();
  EvtComplex Amp(const double& mp, const double& mm);
  EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  double P(const double& mp, const double& mm);
  double Arg(const double& mp, const double& mm);
  double delta(const double& mp, const double& mm);
  void PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta);
  int GetBin(const double& mp, const double& mm);
private:
  EvtComplex amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  void Get4Vs(const double& mp,const double& mm, EvtVector4R& pd, EvtVector4R& pks,EvtVector4R& ppip,EvtVector4R& ppim);
  double del_min, del_max;
  PhaseSpace* phsp;
};

#endif // DALITZMODEL_H
