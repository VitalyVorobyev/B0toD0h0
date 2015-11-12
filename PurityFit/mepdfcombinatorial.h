#ifndef MEPDFCOMBINATORIAL_H
#define MEPDFCOMBINATORIAL_H

#include "../MyParams/myparams.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooArgusBG.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "RooDataSet.h"
#include "RooHist.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace RooFit;

class MEPdfCombinatorial{
public:
  MEPdfCombinatorial(RooRealVar* m_de, RooRealVar* m_mbc, const int mode, const int h0mode);

  RooAddPdf*  GetPdf(void) {return pdf_cmb;}
  RooRealVar* Get_fbb(void) {return fbb_cmb;}
  RooProdPdf* GetPdfBB(void) {return pdf_cmb_bb;}
  RooProdPdf* GetPdfQQ(void) {return pdf_cmb_qq;}

  int TryParameters(RooDataSet *ds);
  int FitParameters(RooDataSet *ds);
  int FitContParameters(RooDataSet *ds);
  int TryContParameters(RooDataSet *ds);

  void FixAll(void) {ChangeParState(kTRUE,0);}
  void FixBB(void)  {ChangeParState(kTRUE,2);}
  void FixQQ(void)  {ChangeParState(kTRUE,1);}
  void FreeAll(void){ChangeParState(kFALSE,0);}
  void FreeBB(void) {ChangeParState(kFALSE,2);}
  void FreeQQ(void) {ChangeParState(kFALSE,1);}

  void Release_de_QQ(void){
    c1_qq_cmb->setConstant(kFALSE);
    c2_qq_cmb->setConstant(kFALSE);
  }
  void Fix_de_QQ(void){
    c1_qq_cmb->setConstant(kTRUE);
    c2_qq_cmb->setConstant(kTRUE);
  }

  void DrawDeltaE(RooDataSet* ds, const bool qqflag);
  void DrawMbc(RooDataSet* ds, const bool qqflag);
  void Draw(RooDataSet* ds, const bool qqflag = false){
    DrawDeltaE(ds,qqflag);
    DrawMbc(ds,qqflag);
    return;
  }

  void PrintParameters(void);
  void WriteParameters(const bool qqflag = false);
  int GetParametersFromFile(void);

private:
  void InitParams(const int mode, const int h0mode);
  void ChangeParState(const int state_flag, const int component = 0);
  void ChangeMbcParState(const int state_flag, const int component = 0);
  void ChangeDeltaEParState(const int state_flag, const int component = 0 );

  double de_line_size(const bool qqflag);
  double mbc_line_size(const bool qqflag);
  MyParams* cuts;
  vector<RooRealVar* > de_param_vec;
  vector<RooRealVar* > mbc_param_vec;

  RooRealVar* de;
  RooRealVar* mbc;
  int m_mode, m_h0mode;

  // de
  // BB
  //RooRealVar* a_c1_bb_cmb;
  //RooRealVar* b_c1_bb_cmb;
  //RooFormulaVar* c1_bb_cmb;
  RooRealVar* c2_bb_cmb;
//  RooChebychev* pdf_de_cmb_bb;
  RooExponential* pdf_de_cmb_bb;

  RooRealVar* const_bb_cmb;
  RooRealVar* de0_bb_cmb;
  RooRealVar* steep_bb_cmb;
//  RooRealVar* f_atan_cmb;
  RooAddPdf*  pdf_de_dst0pi0_bb_cmb;
  RooGenericPdf* pdf_de_atan_bb_cmb;
  RooRealVar* fatan;
  RooAddPdf* pdf_de_dst0_bb_cmb;

  // qq
  RooRealVar* c1_qq_cmb;
  RooRealVar* c2_qq_cmb;
  RooChebychev* pdf_de_cmb_qq;

  //RooRealVar* const_qq_cmb;
  //RooRealVar* de0_qq_cmb;
  //RooRealVar* steep_qq_cmb;
  //RooGenericPdf* pdf_de_dst0_qq_cmb;

  // mbc
  RooRealVar* argedge_cmb;
  RooRealVar* argpar_cmb_bb;
  RooArgusBG* argus_mbc_cmb_bb;

//  RooRealVar* mbc0_cmb_bb;
//  RooRealVar* s_mbc_cmb_bb;
//  RooGaussian* g_mbc_cmb_bb;

//  RooRealVar* fg_cmb_bb;
//  RooAddPdf* pdf_mbc_cmb_bb;
  RooArgusBG* pdf_mbc_cmb_bb;

  RooRealVar* argpar_cmb_qq;
  RooArgusBG* pdf_mbc_cmb_qq;

  RooRealVar* fbb_cmb;
  RooProdPdf* pdf_cmb_bb;
  RooProdPdf* pdf_cmb_qq;
  RooAddPdf* pdf_cmb;
};

#endif // MEPDFCOMBINATORIAL_H
