#ifndef FIT3D_H
#define FIT3D_H

#include "thefcn.h"

class Fit3D{
public:
  Fit3D(TTree *tree);
  double FixPar(const int parnum);
  double FixPar(const char* parname);
  double ReleasePar(const int parnum);
  double ReleasePar(const char* parname);

  FunctionMinimum Fit();

private:
  TheFcn* fcn;
  MnUserParameters* upar;
  MnMigrad* migrad;

  TTree* m_tree;
};

#endif // FIT3D_H
