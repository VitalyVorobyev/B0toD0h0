#ifndef EVTRESONANCE2_H
#define EVTRESONANCE2_H

#include "EvtVector4R.h"
#include "EvtComplex.h"

class EvtResonance2{
public:
  EvtResonance2& operator = (const EvtResonance2 &);
  //constructor with all information about the resonance
  EvtResonance2(const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                const EvtVector4R& p4_d2,
                double ampl = 0.0, double theta = 0.0, double gamma = 0.0,
                double bwm = 0.0, int spin = 0);
  //destructor
  virtual ~EvtResonance2();

  //accessors
    //return 4-momenta of the particles involved
  inline const EvtVector4R& p4_p() { return _p4_p; }
  inline const EvtVector4R& p4_d1() { return _p4_d1; }
  inline const EvtVector4R& p4_d2() { return _p4_d2; }

  //return amplitude
  inline double amplitude() { return _ampl; }

  //return theta
  inline double theta() { return _theta; }

  //return gamma
  inline double gamma() { return _gamma; }

  //return bwm
  inline double bwm() { return _bwm; }

  //return spin
  inline int spin() { return _spin; }

//functions
  //calculate amplitude for this resonance
  EvtComplex resAmpl();

private:
  EvtVector4R _p4_p, _p4_d1, _p4_d2;
  double _ampl, _theta, _gamma, _bwm;
  int _spin;
};

#endif // EVTRESONANCE2_H
