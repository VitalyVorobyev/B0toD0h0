#include "phasespace.h"

PhaseSpace::PhaseSpace(void){
//  TFile* binfile = TFile::Open("/home/vitaly/B0toDh0/binning/dkpp_belle_ddd.root");
  TFile* binfile = TFile::Open("/home/vitaly/B0toDh0/binning/cleo_binning/binning/dkpp_belle_eqdel.root");
//  TFile* binfile = TFile::Open("/home/vitaly/B0toDh0/binning/cleo_binning/binning/dkpp_babar_mod.root");
//  TFile* binfile = TFile::Open("/home/vitaly/B0toDh0/binning/cleo_binning/binning/dkpp_babar_eqdel.root");
  dkpp_bin_h = (TH2F *)binfile->Get("dkpp_bin_h");
  dkpp_bin_h = (TH2F *)dkpp_bin_h->Clone();

  M_D = 1.86484;
//  M_D = 1.864;
  M_D_sq = M_D*M_D;

  M_Ks = 0.497614;
//  M_Ks = 0.50728;//0.505;
  M_Ks_sq = M_Ks*M_Ks;

  M_pi = 0.13957018;
  M_pi_sq = M_pi*M_pi;

  M_min = (M_pi+M_Ks)*(M_pi+M_Ks);
  M_max = (M_D-M_pi)*(M_D-M_pi);
}

inline double PhaseSpace::Mpp(const double& mp,const double& mm){return M_D_sq+2*M_pi_sq+M_Ks_sq-mp-mm;}
inline double PhaseSpace::E_Ks_star(const double& mp){return (mp-M_pi_sq+M_Ks_sq)*0.5/sqrt(mp);}
inline double PhaseSpace::E_pi_star(const double& mp){return (M_D_sq-mp-M_pi_sq)*0.5/sqrt(mp);}

double PhaseSpace::mm_max(const double& mp){
  if(mp > M_max || mp < M_min) return -1;
  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
  const double sum = EKs+Epi;
  const double PKs2 = EKs*EKs-M_Ks_sq;
  const double Ppi2 = Epi*Epi-M_pi_sq;
  const double dif = sqrt(PKs2)-sqrt(Ppi2);
  return (sum-dif)*(sum+dif);
}
double PhaseSpace::mm_min(const double& mp){
  if(mp > M_max || mp < M_min) return -1;
  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
  const double sum = EKs+Epi;
  const double PKs2 = EKs*EKs-M_Ks_sq;
  const double Ppi2 = Epi*Epi-M_pi_sq;
  const double sumP = sqrt(PKs2)+sqrt(Ppi2);
  return (sum-sumP)*(sum+sumP);
}

bool PhaseSpace::IsInPlot(const double& mp,const double& mm){
  if(mm >= mm_min(mp) && mm <= mm_max(mp)) return true;
  else return false;
}
int PhaseSpace::GetBin(const double& mp,const double& mm){
  int gbin = dkpp_bin_h->FindBin(mp+0.00527,mm+0.00527);
//  int gbin = dkpp_bin_h->FindBin(mp,mm);
  return (int)floor(dkpp_bin_h->GetBinContent(gbin)+0.01);
}
//Good: 1630
//Bad:  890

//Good: 1204
//Bad:  2861

