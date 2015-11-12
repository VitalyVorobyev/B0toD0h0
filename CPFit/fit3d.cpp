#include "fit3d.h"

Fit3D::Fit3D(TTree* tree){
  m_tree = tree;
  fcn = new TheFcn(m_tree);
  upar = new MnUserParameters();
  upar->Add("norm",m_tree->GetEntries());
  upar->Add("a1",0.48,0.1,0.,1.);
  upar->Add("a2",0.52,0.1,0.,1.);
  upar->Add("a3",0.00,0.1,-1.,1.);
  upar->Add("a4",0.00,0.1,-1.,1.);
  upar->Add("a5",0.00,0.1,-1.,1.);
  upar->Add("a6",0.00,0.1,-1.,1.);

  migrad = new MnMigrad(*fcn,*upar);
  migrad->Fix("a1");
}

double Fit3D::FixPar(const int parnum){
  migrad->Fix(parnum);
  return upar->Value(parnum);
}

double Fit3D::ReleasePar(const int parnum){
  migrad->Release(parnum);
  return migrad->Value(parnum);
}

double Fit3D::FixPar(const char* parname){
  migrad->Fix(parname);
  return upar->Value(parname);
}

double Fit3D::ReleasePar(const char* parname){
  migrad->Release(parname);
  return migrad->Value(parname);
}

FunctionMinimum Fit3D::Fit(){
  return migrad();
}
