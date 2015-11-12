#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test_cut(){
  const string outfileName("TMVA.root");
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis",outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *inputSig = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/rooksfw/fil_b2dh_sig.root");
  TFile *inputBack = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/rooksfw/fil_b2dh_cont.root");

  TTree* SigTree = (TTree*)inputSig->Get("TEvent");
  TTree* BackTree = (TTree*)inputBack->Get("TEvent");
    
  factory->AddVariable("lh");
//  factory->AddVariable("mbc");

  factory->AddSignalTree(SigTree,1.0);
  factory->AddBackgroundTree(BackTree,1.0);

  string Common_precuts("");//phsp && mbc>5.26 && de<0.1 && de>-0.2 && abs(mpi0_raw-0.135)<0.03 && abs(md0_raw-1.865)<0.03 && chi2_ndf_B0<100");
  string sig_cuts = Common_precuts;// + string(" && ") + string("(mcflag == 1 || mcflag == 10)");
  string back_cuts = Common_precuts;// + string(" && ") + string("b0f != 1 && b0f != 10");
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());

  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
//  factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
  factory->BookMethod( TMVA::Types::kCuts, "Cuts");

//  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggedFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrenght=50:NNodesMax=5");
//  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");
//  factory->BookMethod(TMVA::Types::kBDT,"BDT","!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVAGui( outfileName.c_str() );
}
