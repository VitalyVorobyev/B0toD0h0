#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test2(const int sigreg = true){
  const string outfileName("TMVA.root");
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");

  string name;
  if(!sigreg) name = string("MVAnalysis_lh");
  else name = string("MVAnalysis_lh_sig");
  TMVA::Factory *factory = new TMVA::Factory(name.c_str(),outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *inputSig  = TFile::Open("/home/vitaly/B0toDh0/rooksfw/sig_w_lh.root");
  TFile *inputBack = TFile::Open("/home/vitaly/B0toDh0/rooksfw/cont_w_lh.root");

  TTree* insigtree = (TTree*)inputSig->Get("TEvent");
  TTree* inbacktree = (TTree*)inputBack->Get("TEvent");

  factory->AddVariable("abs(cos_b0)");
  factory->AddVariable("p_ks");
  factory->AddVariable("log(chi2_ndf_D0)");
  factory->AddVariable("log(chi2_ndf_B0)");
  factory->AddVariable("log(chi2_tag_vtx/ndf_tag_vtx)");
  factory->AddVariable("abs(cos_thr)");
  factory->AddVariable("abs(thr_sig-0.885)");
  factory->AddVariable("thr_oth");
  factory->AddVariable("log(tag_LH_err)");
  factory->AddVariable("log(dzerr)");
  factory->AddVariable("log(pi0_chi2)");
  factory->AddVariable("log(egamma)");
//  factory->AddVariable("log(ptgamma)");
  
  factory->AddVariable("lh");

  factory->AddSignalTree(insigtree,1.0);
  factory->AddBackgroundTree(inbacktree,1.0);

  string Common_precuts("chi2_ndf_B0<1000");
  if(sigreg) Common_precuts += string(" && mbc>5.271 && mbc<5.289 && de<0.08 && de>-0.1");
  string sig_cuts = Common_precuts;
  string back_cuts = Common_precuts;
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());
  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

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
