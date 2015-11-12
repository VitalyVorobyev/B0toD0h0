#include "EvtComplex.h"

std::ostream& operator<<(std::ostream& s, const EvtComplex& c){
  s << "(" << c._rpart << "," << c._ipart << ")";
  return s;
}

EvtComplex& EvtComplex::operator*=(EvtComplex c){
  double r=_rpart*c._rpart-_ipart*c._ipart;
  double i=_rpart*c._ipart+_ipart*c._rpart;
  _rpart=r;
  _ipart=i;
  return *this;
}

EvtComplex& EvtComplex::operator/=(EvtComplex c){
  double inv=1.0/(c._rpart*c._rpart+c._ipart*c._ipart);
  double r=inv*(_rpart*c._rpart+_ipart*c._ipart);
  double i=inv*(_ipart*c._rpart-_rpart*c._ipart);
  _rpart=r;
  _ipart=i;
  return *this;
}

EvtComplex EvtComplex::conj(void) const{
  return EvtComplex(this->_rpart,-this->_ipart);
}

