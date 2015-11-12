#include "fit3d.h"

int main(int argc, char** argv){
  TChain* tree = new TChain("tree");
  tree->Add("/home/vitaly/Документы/Vera/work/b0tophiK.root");
  Fit3D fitter(tree);

  FunctionMinimum min = fitter.Fit();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
    return -1;
  }
  MnUserParameterState pstate = min.UserState();
  cout << "Fit results:" << endl;
  for(int i=0; i<7; i++){
    cout << pstate.Value(i) << " = " << pstate.Value(i) << " +- " << pstate.Error(i) << endl;
  }
}
