#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/FunctionMinimum.h"

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"

#include <math.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace ROOT;
using namespace Minuit2;

class pdfFcn : public FCNBase{
public:
  pdfFcn(TTree* tree){
    m_tree = tree;
    m_NTot = m_tree->GetEntries();
    m_tree->SetBranchAddress("mh0",&m_meta);
    m_tree->SetBranchAddress("b0f",&m_b0f);
    m_tree->SetBranchAddress("mode",&m_mode);
    theErrorDef = 1.;
  }
  ~pdfFcn() {}

  double Up() const {return theErrorDef;}
  double operator()(const vector<double>& par) const {
    assert(par.size() == 3);
    const double norm = par[0];
    const double mean = par[1];
    const double widt = par[2];
    double loglh = 0;
    for(int i=0; i<1000; i++){
      m_tree->GetEvent(i);
      if(m_b0f != 1 || m_mode != 3 || isnan(m_meta)) continue;
      const double pdf = norm*1./(TMath::Sqrt(2*TMath::Pi())*widt)*TMath::Exp(-(m_meta-mean)*(m_meta-mean)/(2.*widt*widt));
      cout << pdf << endl;
      loglh += -2*TMath::Log(pdf);
    }
    cout << "loglh: " << loglh << " " << norm << " " << mean << " " << widt << endl;
    return loglh;
  }
private:
  TTree* m_tree;
  double m_meta;
  int m_b0f;
  int m_mode;
  int m_NTot;

  double theErrorDef;
};

