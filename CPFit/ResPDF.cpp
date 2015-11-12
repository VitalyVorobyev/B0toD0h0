#include "ResPDF.h"

double logPoisson(const int n, const double& mu){
  return -2*(n*log(mu)-mu-TMath::LnGamma(n+1));
}

double logGaus(const double& x, const double& x0, const double& s){// -2ln(gauss)
//  cout << TMath::Sqrt(2.*TMath::Pi())*s << " " << 1./(TMath::Sqrt(2.*TMath::Pi())*s)*TMath::Exp(-0.5*(x-x0)*(x-x0)/(s*s)) << endl;
  return (x-x0)*(x-x0)/(s*s) + 2*TMath::Log(TMath::Sqrt(2.*TMath::Pi())*s);
}

double dz_sig_pdf(const double& dz, const double& m1, const double& s1, const double& m2, const double& s2, const double& fg1){
  return fg1*logGaus(dz,m1,s1) + (1-fg1)*logGaus(dz,m2,s2);
}

int calc_dz_sig_pdf(const vector<double>& par, TTree* tree, vector<double>& dz_vec, vector<double>& pdf_vec, const double& dzmin, const double& dzmax, const int ndots){
  assert(par.size() == 5);
  dz_vec.clear(); pdf_vec.clear();
  const int NTot = tree->GetEntries();
  double sz_sig,dz_sig;
  tree->SetBranchAddress("sz_sig",&sz_sig);
  tree->SetBranchAddress("dz_mc_sig1",&dz_sig);
  cout << "m1 = " << par[0] << endl;
  cout << "s1 = " << par[1] << endl;
  cout << "m2 = " << par[2] << endl;
  cout << "s2 = " << par[3] << endl;
  cout << "fr = " << par[4] << endl;

  const double ddz = (dzmax-dzmin)/((double)ndots);

  double dzvec[ndots], pdfvec[ndots];
  for(int i=0; i<ndots; i++){
    dzvec[i] = dzmin + (i+0.5)*ddz;
    pdfvec[i] = 0;
  }

  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    if(abs(dz_sig)>0.5 || isnan(dz_sig) || isnan(sz_sig)) continue;
    for(int j=0; j<ndots; j++){
//      pdfvec[j] += dz_sig_pdf(dzvec[j],par[0],par[1]*sz_sig,par[2],par[3]*sz_sig,par[4]);
      pdfvec[j] += par[4]*TMath::Gaus(dzvec[j],par[0],par[1]*sz_sig) + (1-par[4])*TMath::Gaus(dzvec[j],par[2],par[3]*sz_sig);
    }
  }
  double norm = 0;
  for(int i=0; i<ndots; i++) norm += pdfvec[i];

  norm *= (dzmax-dzmin);
  for(int i=0; i<ndots; i++){
    dz_vec.push_back(dzvec[i]);
    pdf_vec.push_back(pdfvec[i]/norm);
  }

  return 0;
}

double dzSigFcn::operator()(const vector<double>& par) const{
  assert(par.size() == 5);

  const double m1 = par[0];
  const double s1 = par[1];
  const double m2 = par[2];
  const double s2 = par[3];
  const double fg1 = par[4];

  double pdf = 0 ;
  for(int i=0; i<100; i++){
    m_tree->GetEvent(i);
    if(abs(m_dz_sig)>0.5 || isnan(m_dz_sig) || isnan(m_sz_sig)) continue;
    pdf += dz_sig_pdf(m_dz_sig,m1,s1*m_sz_sig,m2,s2*m_sz_sig,fg1);
  }
  cout << "lh: " << pdf << ", s1: " << s1 << ", s2: " << s2 << ", f: " << fg1 << endl;
  return pdf;
}

