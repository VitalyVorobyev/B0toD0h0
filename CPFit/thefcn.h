#ifndef THEFCN_H
#define THEFCN_H

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/FunctionMinimum.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include <cstdlib>

#include "TTree.h"

using namespace std;

class TheFcn : public FCNBase{
public:
  TheFcn(TTree* tree);
  ~pdfFcn() {}
  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const;
  double Pdf(const vector<double>& par);
private:
  TTree* m_tree;
  double theErrorDef;
  double phi, c1, c2;
};

#endif // THEFCN_H
