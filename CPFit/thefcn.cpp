#include "thefcn.h"
#include <cmath>

TheFcn::TheFcn(TTree* tree){
  theErrorDef = 1;
  m_tree = tree;
  NTot = m_tree->GetEntrie();
  m_tree->SetBranchAddress("CosThetaPhi",&c1);
  m_tree->SetBranchAddress("CosThetaKst",&c2);
  m_tree->SetBranchAddress("dPhi",       &phi);
  cout << "Tree with " << NTot << " events will be fitted." << endl;
}

double TheFcn::Pdf(const vector<double>& par){
  const double c1sq = c1*c1;
  const double c2sq = c2*c2;
  const double s1sq = 1.-c1sq;
  const double s2sq = 1.-c2sq;
  const double sd1  = 2*c1*sqrt(s1sq);
  const double sd2  = 2*c2*sqrt(s2sq);

  double pdf = 0;
  pdf += par[1]*c1sq*c2sq;
  pdf += par[2]*0.25*s1sq*s2sq;
  pdf += 0.25*s1sq*s2sq*(par[3]*cos(2.*phi)-par[4]*2*sin(2.*phi));
  pdf += 1./(2.*sqrt(2.))*sd1*sd2*(par[5]*cos(phi)-par[6]*sin(phi));
  return par[0]*pdf;
}

TheFcn::operator()(const vector<double>& par) const {
  double loglh;
  double pdf;
  for(int i=0; i<NTot; i++){
    m_tree->GetEvent(i);
    pdf = Pdf(par);
    if(!isnan(pdf) && pdf>0) loglh += -2.*TMath::Log(pdf);
    else cout << "Bad pdf: " << pdf << end;
  }
  return loglh;
}
