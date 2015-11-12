#ifndef EVTKINE_H
#define EVTKINE_H

#include "EvtVector4R.h"

double EvtDecayAngle(const EvtVector4R&, const EvtVector4R&,
                     const EvtVector4R&);

double EvtDecayAngleChi(const EvtVector4R&, const EvtVector4R&,
                        const EvtVector4R&, const EvtVector4R&,
                        const EvtVector4R& );

//
// This routine calculates the cosine of the angle between
// the normal of the decay plane and the flight direction of particle q
// in the parent frame.
//
double EvtDecayPlaneNormalAngle(const EvtVector4R& p,const EvtVector4R& q,
                          const EvtVector4R& d1,const EvtVector4R& d2);

#endif // EVTKINE_H
