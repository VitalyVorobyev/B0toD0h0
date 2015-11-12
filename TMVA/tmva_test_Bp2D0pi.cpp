#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test_Bp2D0pi(void){
  const string outfileName("TMVA.root");
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");
  string name("MVA_Bp2D0pi");

  TFile *inputSig  = TFile::Open("/home/vitaly/B0toDh0/Bp2D0pi/fil_b2dpi_charged_v2_0_10.root");;
  TFile *inputBack = TFile::Open("/home/vitaly/B0toDh0/Bp2D0pi/fil_b2dpi_charm_0_10_v2.root");
  TMVA::Factory *factory = new TMVA::Factory(name.c_str(),outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TTree* insigtree   = (TTree*)inputSig->Get("TEvent");
  TTree* inbacktree  = (TTree*)inputBack->Get("TEvent");

  factory->AddVariable("abs(costhBcms)");
  factory->AddVariable("log(chisq_vtx_d0)");
  factory->AddVariable("log(chisq_mass_d0)");
  factory->AddVariable("abs(cos_thr)");
  factory->AddVariable("thr_sig");
  factory->AddVariable("thr_oth");
  factory->AddVariable("p_pi");
  factory->AddVariable("k1mm2");
  factory->AddVariable("k1et");
  factory->AddVariable("k1hso00");
//  factory->AddVariable("k1hso01");
  factory->AddVariable("k1hso02");
//  factory->AddVariable("k1hso03");
  factory->AddVariable("k1hso04");
  factory->AddVariable("k1hso10");
  factory->AddVariable("k1hso12");
  factory->AddVariable("k1hso14");
  factory->AddVariable("k1hso20");
  factory->AddVariable("k1hso22");
  factory->AddVariable("k1hso24");
  factory->AddVariable("k1hoo0");
  factory->AddVariable("k1hoo1");
  factory->AddVariable("k1hoo2");
  factory->AddVariable("k1hoo3");
  factory->AddVariable("k1hoo4");

  factory->AddSignalTree(insigtree,1.0);
  factory->AddBackgroundTree(inbacktree,1.0);

  string Common_precuts("mbc>0 && de<0.1 && de>-0.1 && mbc>5.25 && mbc<5.289 && chisq_z_sig<1000 && chisq_z_sig>0 && chisq_z_asc<1000 && chisq_z_asc>0");

  string sig_cuts = Common_precuts;
  string back_cuts = Common_precuts;
  sig_cuts += string(" && bpf == 1");
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());

  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=800:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");

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
