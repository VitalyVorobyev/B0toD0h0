#ifndef __RESPDF_H__
#define __RESPDF_H__

#include <cassert>
#include <iostream>

#include "TTree.h"
#include "TMath.h"
#include "Minuit2/FCNBase.h"

using namespace ROOT;
using namespace Minuit2;
using namespace std;

const double DZMAX = -0.5;
const double DZMIN =  0.5;
const double SZMAX =  0.4;

class dzSigFcn : public FCNBase{
public:
  dzSigFcn(TTree* tree): m_errdef(1){
    m_tree = tree;
    m_NTot = m_tree->GetEntries();
    m_tree->SetBranchAddress("sz_sig",&m_sz_sig);
    m_tree->SetBranchAddress("dz_mc_sig1",&m_dz_sig);
  }
  ~dzSigFcn() {}

  virtual double Up() const {return m_errdef;}
  virtual double operator() (const vector<double> &) const;
private:
  TTree* m_tree;
  double m_errdef;
  int m_NTot;
  double m_sz_sig;
  double m_dz_sig;
};

int calc_dz_sig_pdf(const vector<double>& par, TTree* tree, vector<double>& dz_vec, vector<double>& pdf_vec, const double& dzmin = -0.5, const double& dzmax = 0.5, const int ndots = 100);

#endif

