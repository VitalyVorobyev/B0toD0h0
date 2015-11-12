#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"

using namespace std;

void tmva_test(const int decay = 1, const bool kfsw_flag = true, const int mode = 1){
  // mode = 0 -> Full region
  // mode = 1 -> Soft cut
  // mode = 2 -> Signal region
  // mode = 3 -> mbc signal
  // mode = 4 -> de  signal

  // decay = 1 -> D0 pi0
  // decay = 2 -> D0 eta->gg
  // decay = 3 -> D0 eta->pi+pi-p0
  // decay = 4 -> D0 omega
  string outfileName;
  string name;
  string mode_cut;
  switch(mode){
    case 0:
      name = string("MVA_fullrange");
      break;
    case 1:
      name = string("MVA_softcut");
      break;
    case 2:
      name = string("MVA_signalreg");
      break;
    case 3:
      name = string("MVA_mbcsig");
      break;
    case 4:
      name = string("MVA_desig");
      break;
    case 5:
      name = string("MVA_fullrange_signal");
      break;
    default:
      cout << "Unknown mode " << mode << endl;
      return;
  }
//  TFile *inputSig;
//  TFile *inputSig2;
//  TFile *inputBack = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_cont_0-1.root");

//  TFile *inputBack_charm;
//  TFile *inputBack_uds;
  TChain* sig_chain = new TChain("TEvent","TEvent");
  TChain* bkg_chain = new TChain("TEvent","TEvent");

  if(kfsw_flag){
//    inputBack_charm = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_charm_2_12.root");
//    inputBack_uds = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_uds_2_12.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_3_13.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_3_13.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_4_14.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_4_14.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_5_15.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_5_15.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_charm_0_0.root");
    bkg_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_uds_0_0.root");
  } else{
    inputBack_charm = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_charm_2_12.root");
    inputBack_uds = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_uds_2_12.root");
  }
//  name += string("_Domega");
//  inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcOmega_s5.root");
//  outfileName = string("TMVA_omega.root");

  switch(decay){
  case 1:
    if(kfsw_flag){
      name += string("_ksfw_Dpi0");
      sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcPi0_s7.root");
//      inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigPi0_s7.root");
      outfileName = string("TMVA_ksfw_pi0.root");
    } else{
    name += string("_Dpi0");
    inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcPi0_s7.root");
    outfileName = string("TMVA_pi0.root");
    }
    break;
  case 10:// D*0 pi0
    name += string("_ksfw_Dst0pi0");
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    outfileName = string("TMVA_ksfw_dst0pi0.root");
    break;
  case 2:
    if(kfsw_flag){
      name += string("_ksfw_Detagg");
      sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s2.root");
      //inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigEta_s2.root");
      outfileName = string("TMVA_ksfw_etagg.root");
    } else{
    name += string("_Detagg");
    inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcEta_s2.root");
    outfileName = string("TMVA_etagg.root");
    }
    break;
  case 11:// D*0 eta(->gg)
    name += string("_ksfw_Dst0etagg");
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    outfileName = string("TMVA_ksfw_dst0etagg.root");
    break;
  case 12:// D*0 eta(->ppp)
    name += string("_ksfw_Dst0etappp");
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1.root");
    outfileName = string("TMVA_ksfw_dst0etappp.root");
    break;
  case 20:// D0 eta'(->gg)
    name += string("_ksfw_etapgg");
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETAP_s1.root");
    outfileName = string("TMVA_ksfw_etapgg.root");
    break;
  case 21:// D0 eta'(->ppp)
    name += string("_ksfw_etapppp");
    sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETAP_s1.root");
    outfileName = string("TMVA_ksfw_etapppp.root");
    break;
  case 3:
    if(kfsw_flag){
      name += string("_ksfw_Detappp");
      sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s2.root");
      //inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigEta_s2.root");
      outfileName = string("TMVA_ksfw_etappp.root");
    } else{
      name += string("_Detappp");
      inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcEta_s2.root");
      outfileName = string("TMVA_etappp.root");
    }
    break;
  case 4:
    if(kfsw_flag){
      name += string("_ksfw_Domega");
      sig_chain->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcOMEGA_s5.root");
      //inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigOmega_s5.root");
      outfileName = string("TMVA_ksfw_omega.root");
    } else{
    name += string("_Domega");
    inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcOmega_s5.root");
    outfileName = string("TMVA_omega.root");
    }
    break;
  case 5:// pi0 + etagg
    if(kfsw_flag){
      inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigPi0_s7.root");
      inputSig2 = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigEta_s2.root");
      outfileName = string("TMVA_ksfw_pi0_etagg.root");
      name += string("_ksfw_D_pi0_etagg");
    } else{
      inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcPi0_s7.root");
      inputSig2 = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcEta_s2.root");
      outfileName = string("TMVA_pi0_etagg.root");
      name += string("_D_pi0_etagg");
    }
    break;
  case 6:// etappp + omega
    if(kfsw_flag){
      name += string("_ksfw_D_etappp_omega");
      inputSig = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigEta_s2.root");
      inputSig2 = TFile::Open("/home/vitaly/B0toDh0/TMVA/Fil_b2dh_sigOmega_s5.root");
      outfileName = string("TMVA_ksfw_etappp_omega.root");
    } else{
      name += string("_D_etappp_omega");
      inputSig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcEta_s2.root");
      inputSig2 = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmcOmega_s5.root");
      outfileName = string("TMVA_etappp_omega.root");
    }
    break;
  default:
    cout << "Unknown decay " << decay << endl;
    return;
  }
  TFile* outputFile = TFile::Open(outfileName.c_str(),"RECREATE");

  TMVA::Factory *factory = new TMVA::Factory(name.c_str(),outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  if(!kfsw_flag){
    TTree* insigtree   = (TTree*)inputSig->Get("TEvent");
    if(decay == 5 || decay == 6){
      TTree* insigtree2   = (TTree*)inputSig2->Get("TEvent");
    }
    TTree* inbacktree_charm   = (TTree*)inputBack_charm->Get("TEvent");
    TTree* inbacktree_uds     = (TTree*)inputBack_uds->Get("TEvent");
  }

//  insigtree->Print();

  factory->AddVariable("abs(costhBcms)");
  factory->AddVariable("log(chi2_mass_d0)");
  factory->AddVariable("abs(cos_thr)");
  factory->AddVariable("thr_sig");
//  factory->AddVariable("thr_oth");
//  factory->AddVariable("p_h0");
  if(decay == 2 || decay == 11 || decay == 20) factory->AddVariable("log(h0_chi2)");
//  else           factory->AddVariable("log(pi0_chi2)");
//  factory->AddVariable("log(tag_LH_err)");
  factory->AddVariable("log(egamma)");
  if(decay>2 && decay != 5 && decay != 11 && decay != 20 && decay != 10){
    factory->AddVariable("p_pi0_h0");
//    factory->AddVariable("p_pip_h0");
//    factory->AddVariable("p_pim_h0");
  }
  if(decay==4) factory->AddVariable("abs(cos_hel)");

  if(kfsw_flag){
//    if(decay > 2) factory->AddVariable("1./(1.005 - lh1*lh1)");
//    else factory->AddVariable("lh1");
    factory->AddVariable("lh0");
  } else{
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
  }

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

  if(kfsw_flag){
    factory->AddSignalTree(sig_chain,1.0);
    factory->AddBackgroundTree(bkg_chain,1.0);
  } else{
    factory->AddBackgroundTree(inbacktree_charm,1.0);
    factory->AddBackgroundTree(inbacktree_uds,1.0);

    factory->AddSignalTree(insigtree,1.0);
    if(decay == 5 || decay == 6) factory->AddSignalTree(insigtree2,1.0);
  }

//  string Common_precuts("mbc>5.271 && mbc<5.289 && de<0.08 && de>-0.1 && chi2_mass_d0>0 && abs(mks_raw-0.4975)<0.009 && abs(md0_raw-1.865)<0.015 && abs(mpi0_raw-0.135)<0.012");
  string Common_precuts = string("mbc>0");// + mode_cut;
  Common_precuts += string(" && egamma>0.05 && chi2_vtx_d0<50 && abs(md_raw-1.865)<0.015 && good_icpv==1 && chi2_vtx_d0>0 && chi2_mass_d0>0 && log(chi2_mass_d0)<50");
  if(decay == 2) Common_precuts += string(" && egamma>0.08 && log(h0_chi2)>-10");
  switch(mode){
    case 0:// Full range
      break;
    case 1:// Soft cut
      Common_precuts += string(" && de<0.1 && de>-0.1 && mbc>5.271  && mbc<5.289");
      break;
    case 2:// Signal region
      Common_precuts += string(" && de<0.08 && de>-0.1 && mbc>5.271 && mbc<5.289");
      break;
    case 3:// Mbc signal
      Common_precuts += string(" && mbc>5.271 && mbc<5.289");
      break;
    case 4:// dE signal
      Common_precuts += string(" && de<0.08 && de>-0.1");
      break;
  }
  switch(decay){
  case 1:// pi0
    Common_precuts += string(" && mode == 1");
    break;
  case 10:// D*0 pi0
    Common_precuts += string(" && mode == 10");
    break;
  case 2:// eta->gg
    Common_precuts += string(" && mode == 2 && h0mode == 10");
    break;
  case 11:// D*0 eta(->gg)
    Common_precuts += string(" && mode == 20 && h0mode == 10");
    break;
  case 12:// D*0 eta(->ppp)
    Common_precuts += string(" && mode == 20 && h0mode == 20");
    break;
  case 20:// D0 eta'(->gg)
    Common_precuts += string(" && mode == 5 && h0mode == 10");
    break;
  case 21:// D0 eta'(->ppp)
    Common_precuts += string(" && mode == 5 && h0mode == 20");
    break;
  case 3:// eta->ppp0
    Common_precuts += string(" && mode == 2 && h0mode == 20");
    break;
  case 4:// omega
    Common_precuts += string(" && mode == 3");
    break;
  case 5:// pi0 + etagg
    Common_precuts += string(" && h0mode == 10");
    break;
  case 6:// etappp + omega
    Common_precuts += string(" && h0mode == 20");
    break;
  }
  string sig_cuts = Common_precuts;
  string back_cuts = Common_precuts;
  if(mode == 5) sig_cuts += string("&& (b0f == 1 || b0f == 5 || b0f == 10) && de<0.08 && de>-0.1 && mbc>5.271 && mbc<5.289");
  cout << "sig cuts: " << sig_cuts << endl;
  cout << "bkg cuts: " << back_cuts << endl;
  TCut cutsig(sig_cuts.c_str());
  TCut cutback(back_cuts.c_str());

  factory->PrepareTrainingAndTestTree(cutsig,cutback,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
//  factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
//  factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggedFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrenght=50:NNodesMax=5");

  if(!kfsw_flag){
    factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=800:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");
  } else{
    factory->BookMethod(TMVA::Types::kBDT,"BDT","!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");
//    factory->BookMethod(TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=800:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:NNodesMax=5");
  }

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
