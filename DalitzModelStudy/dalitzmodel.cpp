#include "dalitzmodel.h"

//const int BINS[8] = {8,7,6,5,4,3,2,1};
const int BINS[8] = {8,1,2,3,4,5,6,7};

DalitzModel::DalitzModel(){
  del_min = -M_PI/8.; del_max = 15.*M_PI/8.;
  phsp = new PhaseSpace();
}

int DalitzModel::GetBin(const double& mp, const double& mm){
  if(!phsp->IsInPlot(mp,mm)) return 0;
  const double delt = mp>mm ? delta(mp,mm) : delta(mm,mp);
  for(int i=1; i<=8; i++){
    if(M_PI*(i-1.5)/4 < delt && delt < M_PI*(i-0.5)/4){
      if(mp>mm) return  BINS[i-1];
      else      return -BINS[i-1];
    }
  }
  cout << "delta  = " << delt << ", mp: " << mp << ", mm: " << mm << endl;
  return 0;
}

EvtComplex DalitzModel::Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  return amp_Belle2010(p4_p,moms1,moms2,moms3);
}

EvtComplex DalitzModel::amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  // ** A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010 **
  vector<EvtResonance2> m_res_v;
//  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 1.638, 133.2, 0.0508, .89166, 1));//K*(892)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 1.638, 133.2, 0.0484, .8937, 1));//K*(892)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 2.210, 358.9, 0.294 , 1.412 , 0));//K0*(1430)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 0.890, 314.8, 0.0985, 1.4256, 2));//K2*(1430)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 0.880,  82.0, 0.322 , 1.717 , 1));//K*(1680)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms2, 0.650, 120.0, 0.232 , 1.414 , 1));//K*(1410)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms3, 0.149, 325.4, 0.0508, .89166, 1));//DCS K*(892)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms3, 0.360,  87.0, 0.294 , 1.412 , 0));//DCS K0*(1430)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms3, 0.230, 275.0, 0.0985, 1.4256, 2));//DCS K2*(1430)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms3, 2.100, 130.0, 0.322 , 1.717 , 1));//DCS K*(1680)
  m_res_v.push_back(EvtResonance2(p4_p,moms1,moms3, 0.420, 253.0, 0.232 , 1.414 , 1));//DCS K*(1410)
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.000,   0.0, 0.1490, 0.7717, 1));//Rho
//  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.000,   0.0, 0.1360, 0.7717, 1));//Rho
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, .0343, 112.0, .00849, .78265, 1));//Omega
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 0.385, 207.3, 0.05  , 0.977,  0));//f0(980)
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.250,  69.0, 0.272 , 1.31 ,  0));//f0(1370)
//  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.56,  110.0, 0.272 , 1.31 ,  0));//f0(1370)
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.440, 342.9, 0.1851, 1.2754, 2));//f2(1270)
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 0.490,  64.0, 0.400 , 1.465,  1));//Rho(1450)
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 1.560, 214.0, 0.453 , 0.522,  0));//sigma1
  m_res_v.push_back(EvtResonance2(p4_p,moms3,moms2, 0.200, 212.0, 0.088 , 1.033,  0));//sigma2
  EvtComplex amp(-2.537,0.923);
  for(int i=0; i<18; i++){ amp += m_res_v[i].resAmpl();}
  return amp;
}

EvtComplex DalitzModel::Amp(const double& mp, const double& mm){
  EvtVector4R pd;
  EvtVector4R pks;
  EvtVector4R ppip;
  EvtVector4R ppim;
  Get4Vs(mp,mm,pd,pks,ppip,ppim);
  return Amp(pd,pks,ppip,ppim);
}

double DalitzModel::P(const double& mp,const double& mm){
  return abs2(Amp(mp,mm));
}

double DalitzModel::Arg(const double& mp,const double& mm){
  return arg(Amp(mp,mm));
}

double DalitzModel::delta(const double& mp,const double& mm){
  EvtVector4R pd;
  EvtVector4R pks;
  EvtVector4R ppip;
  EvtVector4R ppim;
  Get4Vs(mp,mm,pd,pks,ppip,ppim);
  double del = - (arg(Amp(pd,pks,ppip,ppim)) - arg(Amp(pd,pks,ppim,ppip)));
  if(del<del_min) return del + 2.*M_PI;
  if(del>del_max) return del - 2.*M_PI;
  return del;
}

void DalitzModel::PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta){
  EvtVector4R pd;
  EvtVector4R pks;
  EvtVector4R ppip;
  EvtVector4R ppim;
//  if(mp)
  Get4Vs(mp,mm,pd,pks,ppip,ppim);
  EvtComplex A  = Amp(pd,pks,ppip,ppim);
  EvtComplex Ab = Amp(pd,pks,ppim,ppip);
  delta = arg(A) - arg(Ab);
  if(delta<del_min) delta += 2.*M_PI;
  if(delta>del_max) delta -= 2.*M_PI;
  P = abs2(A); Pbar = abs2(Ab);
  return;
}

void DalitzModel::Get4Vs(const double& mp,const double& mm, EvtVector4R& pd, EvtVector4R& pks,EvtVector4R& ppip,EvtVector4R& ppim){
  const double M_D = phsp->m_D();
  const double M_D_sq = M_D*M_D;
  const double M_pi = phsp->m_pi();
  const double M_pi_sq = M_pi*M_pi;
  const double M_Ks = phsp->m_Ks();
  const double M_Ks_sq = M_Ks*M_Ks;

  double eks,pxks = 0 ,pyks = 0,pzks;
  double epip,pxpip,pypip = 0,pzpip;
  double epim,pxpim,pypim = 0,pzpim;
  double ed = M_D,pxd = 0,pyd = 0,pzd = 0;
  double mm_test, mp_test;

  eks   = (mp+mm-2*M_pi_sq)/(2*M_D);
  epip  = (M_D_sq+M_pi_sq-mm)/(2*M_D);
  epim  = (M_D_sq+M_pi_sq-mp)/(2*M_D);
  pzks  = sqrt(eks*eks - M_Ks_sq);
  pzpip = (M_pi_sq+M_Ks_sq+2*eks*epip-mp)/(2*pzks);
  pzpim = (M_pi_sq+M_Ks_sq+2*eks*epim-mm)/(2*pzks);
  pxpip = sqrt(epip*epip - pzpip*pzpip - M_pi_sq);
  pxpim = -pxpip;

  mm_test = (eks+epim)*(eks+epim)-(pxks+pxpim)*(pxks+pxpim)-(pyks+pypim)*(pyks+pypim)-(pzks+pzpim)*(pzks+pzpim);
  mp_test = (eks+epip)*(eks+epip)-(pxks+pxpip)*(pxks+pxpip)-(pyks+pypip)*(pyks+pypip)-(pzks+pzpip)*(pzks+pzpip);

  if(abs(mm-mm_test)>0.0001 || abs(mp-mp_test)>0.0001){
    const double mpip_test = sqrt(epip*epip-pxpip*pxpip-pypip*pypip-pzpip*pzpip);
    const double mpim_test = sqrt(epim*epim-pxpim*pxpim-pypim*pypim-pzpim*pzpim);
    const double mks_test = sqrt(eks*eks-pxks*pxks-pyks*pyks-pzks*pzks);
    cout << "Wrong (mm,mp): (" << mm << "," << mp << ") -> (" << mm_test << "," << mp_test << "):" << endl;
    cout << " pi+: (" << epip << "," << pxpip << "," << pypip << "," << pzpip << ") -> " << mpip_test << endl;
    cout << " pi-: (" << epim << "," << pxpim << "," << pypim << "," << pzpim << ") -> " << mpim_test << endl;
    cout << " k0s: (" << eks  << "," << pxks  << "," << pyks  << "," << pzks  << ") -> " << mks_test  << endl;
  }
  pd   = EvtVector4R(ed,pxd,pyd,pzd);
  pks  = EvtVector4R(eks,pxks,pyks,pzks);
  ppip = EvtVector4R(epip,pxpip,pypip,pzpip);
  ppim = EvtVector4R(epim,pxpim,pypim,pzpim);
  return;
}
