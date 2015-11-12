#include "stdio.h"

#include "purityfit.h"

int main(int argc, char** argv){
  int mode = 1;
  int nstr = 1;
  int cstr = 1;
//  char data[10];
  if(argc == 2) sscanf(argv[1],"%d",&mode);
//  if(argc == 3){
//    sscanf(argv[1],"%d",&mode);
//    sscanf(argv[2],"%s",data);
//  }
  if(argc == 4){
    sscanf(argv[1],"%d",&mode);
    sscanf(argv[2],"%d",&nstr);
    sscanf(argv[3],"%d",&cstr);
  }
//  PurityFit pfit(mode);
//  pfit.MakeFit(true);
//  return 0;

  // mode 1  -> pi0
  // mode 10 -> D*0pi0
  // mode 2  -> eta -> gg
  // mode 20 -> D*0 eta(->gg)
  // mode 21 -> D*0 eta(->ppp)
  // mode 3  -> eta -> pi+pi-pi0
  // mode 4  -> omega
  // mode 5  -> eta'(->gg)

//  const bool line_flag = true;

//  vector<string> angle_vec;
//  angle_vec.push_back(string("0"));
//  angle_vec.push_back(string("15"));
//  angle_vec.push_back(string("22_5"));
//  angle_vec.push_back(string("30"));
//  angle_vec.push_back(string("45"));
//  angle_vec.push_back(string("60"));
//  angle_vec.push_back(string("67_5"));
//  angle_vec.push_back(string("75"));
//  angle_vec.push_back(string("90"));
//  angle_vec.push_back(string("105"));
//  angle_vec.push_back(string("112_5"));
//  angle_vec.push_back(string("120"));
//  angle_vec.push_back(string("135"));
//  angle_vec.push_back(string("150"));
//  angle_vec.push_back(string("157_5"));
//  angle_vec.push_back(string("165"));
//  for(int i=0; i<angle_vec.size(); i++){
//    PurityFit pfit(mode);
//    pfit.SaveSigCPVTree(true,angle_vec[i]);
//  }
//  PurityFit pfit(mode);
//////  pfit.SaveSigLineCPVTrees();
////  return 0;

//  pfit.SaveSigCPVTree();

//  return 0;
//  pfit.ReadAllParamsFromFile();
////  pfit.SaveSidebandTree();
////  pfit.CheckSigPdf();
////  pfit.MakeFit();

//  PurityFit pfit(mode);
//  pfit.SaveSigCPVTree();
//  return 0;

//  pfit.FitSigPdf();

//  return 0;

//  vector<int> streams;
//  streams.push_back(0);
//  streams.push_back(1);
//  streams.push_back(2);
//  streams.push_back(3);
////  streams.push_back(4);
////  streams.push_back(5);

////////  pfit.CheckPartPdf(streams);
//  pfit.FitPartPdf(streams);
//  pfit.FitCombPdf(streams);
//  pfit.MakeGenMCPutiryFit(streams,0,1,0,0,true);// No mixing

//return 0;

  PurityFit pfit2(mode);
  vector<int> streams2;
//  if(nstr == 6){
//    streams2.push_back(cstr);
//  } else{
    streams2.push_back(0);
    streams2.push_back(1);
    streams2.push_back(2);
    streams2.push_back(3);
    streams2.push_back(4);
    streams2.push_back(5);
//  }
// //  pfit2.ZeroSigTest(streams2,0,true);
//  pfit2.MakeGenMCPutiryFit(streams2,0,1,0,0,true);// No mixing
// //  pfit2.MakeGenMCPutiryFit(streams2,0,2,0,0,true,nstr,cstr);// Mixing
  pfit2.MakeGenMCPutiryFit2(streams2,true,true);// Mixing

//  return 0;

////  pfit.MakeFit(true); return 0;

//  pfit.FitPartPdf(streams);
//  pfit.FitCombPdf(streams);
////  pfit.ZeroSigTest(streams);
////  pfit.MakeGenMCPutiryFit(streams,0,1,0,0,true);// No mixing
//  pfit.MakeGenMCPutiryFit(streams,0,2,0,0,true);// Mixing
////  pfit.MakeGenMCPutiryFit(streams,0,1,0,1,true);// No mixing cont
////  pfit.MakeGenMCPutiryFit(streams,0,1,0,2,true);// No mixing BB
////  pfit.MakeGenMCPutiryFit(streams,0,2,0,0,false);// Mixing
////  pfit.SaveSidebandTree(streams,0,0);// qq + BB
////  pfit.SaveSidebandTree(streams,0,1);// qq
////  pfit.SaveSidebandTree(streams,0,2);// BB

////  RooDataSet* ds = pfit.GetBBBackMCDataSet(streams);
////  pfit.DefineElliRange();
////  pfit.CountTrueNumbers(ds);
////  pfit.PrintTrueNumbers();
////
////  pfit.CheckAllPdfs(streams);
////  pfit.CheckCombPdf(streams);
////  pfit.FitContPdf(streams);
////  cout << "FitCombPdf" << endl;
////  pfit.FitCombPdf(streams);
////  pfit.CheckPartPdf(streams);
////  pfit.FitPartPdf(streams);

//  return 0;
}
