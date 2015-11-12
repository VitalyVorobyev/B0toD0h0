#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test_ksfw(){
  const string outfileName("TMVA.root");
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");
  
  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis_ksfw",outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *inputSig = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/Tuples/b2dh_sigmc_s1.root");
  TFile *inputBack = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/Tuples/b2dh_charm.root");
  TFile *inputBack1 = TFile::Open("/home/vitaly/Belle_analysis/B0toDh0_Belle/Tuples/b2dh_uds.root");

//  TFile *tfile = new TFile("","RECREATE");

  TTree* insigtree = (TTree*)inputSig->Get("TEvent");
  TTree* inbacktree = (TTree*)inputBack->Get("TEvent");
  TTree* inbacktree1 = (TTree*)inputBack->Get("TEvent");
  
  TTree* SigTree  = insigtree.CopyTree("mcflag == 1 || mcflag == 10 || mcflag == 5");
  TTree* BackTree = inbacktree.CopyTree("b0f != 1 && b0f != 10 && !(d0f == 1 && (b0f == 5 || b0f == 4))");
  TTree* BackTree1 = inbacktree1.CopyTree("b0f != 1 && b0f != 10 && !(d0f == 1 && (b0f == 5 || b0f == 4))");
  
//  factory->AddVariable("abs(cos_b0)");
//  factory->AddVariable("abs(cos_d0)");
//  factory->AddVariable("abs(pcm_d0-2.3)");
//  factory->AddVariable("p_ks");
//  factory->AddVariable("p_pp");
//  factory->AddVariable("p_pm");
//  factory->AddVariable("abs(p_pi0-2.6)");
//  factory->AddVariable("log(abs(r_pp))");
//  factory->AddVariable("log(abs(z_pp))");
//  factory->AddVariable("log(abs(r_pm))");
//  factory->AddVariable("log(abs(z_pm))");
//  factory->AddVariable("log(chi2_ndf_D0)");
//  factory->AddVariable("log(chi2_ndf_B0)");
//  factory->AddVariable("log(chi2_tag_vtx)");
//  factory->AddVariable("abs(cos_thr)");
//  factory->AddVariable("thr_sig");
//  factory->AddVariable("thr_oth");
//  factory->AddVariable("log(tag_LH_err)");
//  factory->AddVariable("log(dzerr)");
  
//  factory->AddVariable("log(ks_dr)");
//  factory->AddVariable("log(ks_dz)");
//  factory->AddVariable("log(ks_dphi)");
//  factory->AddVariable("log(ks_fl)");

  factory->AddVariable("k0mm2");
  factory->AddVariable("k0et");
  factory->AddVariable("k0hso00");
//  factory->AddVariable("k0hso01");
  factory->AddVariable("k0hso02");
//  factory->AddVariable("k0hso03");
  factory->AddVariable("k0hso04");
  factory->AddVariable("k0hso10");
  factory->AddVariable("k0hso12");
  factory->AddVariable("k0hso14");
  factory->AddVariable("k0hso20");
  factory->AddVariable("k0hso22");
  factory->AddVariable("k0hso24");
  factory->AddVariable("k0hoo0");
  factory->AddVariable("k0hoo1");
  factory->AddVariable("k0hoo2");
  factory->AddVariable("k0hoo3");
  factory->AddVariable("k0hoo4");

  factory->AddSignalTree(SigTree,1.0);
  factory->AddBackgroundTree(BackTree,1.0);
  factory->AddBackgroundTree(BackTree1,1.0);
  
  string Common_precuts("mbc>5.271 && mbc<5.289 && de<0.08 && de>-0.1 && abs(mpi0_raw-0.135)<0.012 && abs(md0_raw-1.865)<0.015 && abs(mks_raw-0.4975)<0.009 && chi2_ndf_B0<1000");
  string sig_cuts = Common_precuts;// + string(" && ") + string("(mcflag == 1 || mcflag == 10)");
  string back_cuts = Common_precuts;// + string(" && ") + string("b0f != 1 && b0f != 10");
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());

  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
//  factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

//  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggedFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrenght=50:NNodesMax=5");
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
