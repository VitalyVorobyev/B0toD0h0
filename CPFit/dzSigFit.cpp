#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
//#include "Minuit2/MinimumPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnPrint.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "ResPDF.h"

using namespace ROOT;
using namespace Minuit2;
using namespace std;

int main(void){

  TFile file("../Tuples/fil_b2dh_sigmcOmega_s1new.root");
  TTree* tree = (TTree*)file.Get("TEventTr");

  MnUserParameters upar;
  upar.Add("mean1",0.,0.1);
  upar.Add("scale1",1.92,1);
  upar.Add("mean2",0.,0.1);
  upar.Add("scale2",10.,5);
  upar.Add("frac",1.,0.1);

  upar.Fix("mean1");
  upar.Fix("mean2");
  upar.Fix("frac");
  upar.Fix("scale2");
  upar.Fix("scale1");

  upar.SetLimits("frac",0.,1.);
  upar.SetLimits("scale1",0.1,10.);
  upar.SetLimits("scale2",3.1,100.);
  upar.SetLimits("mean2",-0.1,0.1);

  dzSigFcn theFCN(tree);

  MnMigrad migrad(theFCN,upar);

  FunctionMinimum min = migrad();
  MnUserParameterState parstate = min.UserState();

  cout << "Minimization status: " << min.IsValid() << endl;
  cout << "m1 = " << parstate.Value(0) << " +- " << parstate.Error(0) << endl;
  cout << "s1 = " << parstate.Value(1) << " +- " << parstate.Error(1) << endl;
  cout << "m2 = " << parstate.Value(2) << " +- " << parstate.Error(2) << endl;
  cout << "s2 = " << parstate.Value(3) << " +- " << parstate.Error(3) << endl;
  cout << "f = "  << parstate.Value(4) << " +- " << parstate.Error(4) << endl;

  vector<double> dzvec, pdfvec, parvec;
  for(int i=0; i<5; i++) parvec.push_back(parstate.Value(i));
  calc_dz_sig_pdf(parvec,tree,dzvec,pdfvec);

  ofstream ofile;
  ofile.open("dz_sig_rf.txt");
  for(int i=0; i<dzvec.size(); i++){
    ofile << dzvec[i] << " " << pdfvec[i] << endl;
  }
  ofile.close();

  return 0;
}

