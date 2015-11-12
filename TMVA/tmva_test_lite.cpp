#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test_lite(const int sigreg = true){
  const string outfileName("TMVA.root");
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");
  string name;
  if(!sigreg) name = string("MVAnalysisLite");
  else name = string("MVAnalysis_sig_Lite");
  TMVA::Factory *factory = new TMVA::Factory(name.c_str(),outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmc.root");
  TFile *inputBack = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_cont.root");

//  TFile *tfile = new TFile("","RECREATE");

  TTree* insigtree   = (TTree*)inputSig->Get("TEventTr");
  TTree* inbacktree  = (TTree*)inputBack->Get("TEventTr");

  factory->AddVariable("abs(cos_b0)");
//  factory->AddVariable("p_ks");
  factory->AddVariable("log(chi2_ndf_D0)");
//  factory->AddVariable("log(chi2_ndf_B0)");
//  factory->AddVariable("log(chi2_tag_vtx/ndf_tag_vtx)");
  factory->AddVariable("abs(cos_thr)");
  factory->AddVariable("thr_sig");
  factory->AddVariable("thr_oth");
  factory->AddVariable("log(tag_LH_err)");
//  factory->AddVariable("log(dzerr)");
//  factory->AddVariable("log(pi0_chi2)");
//  factory->AddVariable("log(egamma)");
//  factory->AddVariable("log(ptgamma)");
  
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

//  factory->AddVariable("k0mm2");
//  factory->AddVariable("k0et");
//  factory->AddVariable("k0hso00");
//  factory->AddVariable("k0hso02");
//  factory->AddVariable("k0hso04");
//  factory->AddVariable("k0hso10");
//  factory->AddVariable("k0hso12");
//  factory->AddVariable("k0hso14");
//  factory->AddVariable("k0hso20");
//  factory->AddVariable("k0hso22");
//  factory->AddVariable("k0hso24");
//  factory->AddVariable("k0hoo0");
//  factory->AddVariable("k0hoo1");
//  factory->AddVariable("k0hoo2");
//  factory->AddVariable("k0hoo4");

  factory->AddSignalTree(insigtree,1.0);
  factory->AddBackgroundTree(inbacktree,1.0);

//  string Common_precuts("mbc>5.271 && mbc<5.289 && de<0.08 && de>-0.1 && chi2_ndf_B0<1000 && abs(mks_raw-0.4975)<0.009 && abs(md0_raw-1.865)<0.015 && abs(mpi0_raw-0.135)<0.012");
  string Common_precuts("chi2_ndf_B0<1000 && de<0.1 && de>-0.1 && mbc>5.25 && mbc<5.289");
  if(sigreg) Common_precuts += string(" && mbc>5.271 && mbc<5.289 && de<0.08 && de>-0.1");
  string sig_cuts = Common_precuts;// + string("");
  string back_cuts = Common_precuts;// + string(" && abs(md0_raw-1.865)<0.030 && abs(mpi0_raw-0.135)<0.030");
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());

  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
//  factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

//  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggedFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrenght=50:NNodesMax=5");
  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=800:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");
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
