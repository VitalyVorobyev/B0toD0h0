#ifndef MEPDFPARTIALB_H
#define MEPDFPARTIALB_H

#include <iostream>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooNovosibirsk.h"
#include "RooArgusBG.h"
#include "RooGaussian.h"
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

#include "../MyParams/myparams.h"
#include "../BtoDh/pdfs/RooRhoDeltaEPdf.h"

using namespace std;
using namespace RooFit;

class MEPdfPartialB{
public:
  MEPdfPartialB(RooRealVar* m_de, RooRealVar* m_mbc, const int mode, const int h0mode);

  RooProdPdf* GetPdf(void) {return pdf_part;}

  int TryParameters(RooDataSet *ds);
  int FitParameters(RooDataSet *ds, const bool etaggflag = false);

  void FixAll(void) {ChangeParState(kTRUE);}
  void FreeAll(void){ChangeParState(kFALSE);}

  void DrawDeltaE(RooDataSet* ds);
  void DrawMbc(RooDataSet* ds);
  void Draw(RooDataSet* ds){
    DrawDeltaE(ds);
    DrawMbc(ds);
    return;
  }

  void PrintParameters(void);
  void WriteParameters(void);
  int GetParametersFromFile(void);

  void Release_de0(void)     {de0_part->setConstant(kFALSE);}
  void Fix_de0(void)         {de0_part->setConstant(kTRUE);}
  void Release_k_s_mbc(void) {k_s_mbc_part->setConstant(kFALSE);}
  void Fix_k_s_mbc(void)     {k_s_mbc_part->setConstant(kTRUE);}
  void Release_alpha(void)   {alpha_part->setConstant(kFALSE);}
  void Fix_alpha(void)       {alpha_part->setConstant(kTRUE);}

private:
  void InitParams(const int mode, const int h0mode);
  void ChangeParState(const int state_flag);

  MyParams* cuts;
  vector<RooRealVar* > de_param_vec;
  vector<RooRealVar* > mbc_param_vec;

  RooRealVar* de;
  RooRealVar* mbc;
  bool ggflag;
  int m_mode, m_h0mode;

  RooRealVar* de0_part;
  RooRealVar* slopel_part;
  RooRealVar* sloper_part;
  RooRealVar* steep_part;
  RooRealVar* p5_part;
  RooRhoDeltaEPdf* pdf_de_part;

  // gg
  RooRealVar* b_s_mbc_part;
  RooRealVar* k_s_mbc_part;
  RooFormulaVar* s_mbc_part;
  RooRealVar* alpha_part;
  RooRealVar* b_mbc0_part;
  RooRealVar* k_mbc0_part;
  RooFormulaVar* MBC0_part;
  RooNovosibirsk* pdf_mbc_part_gg;

  // pi+pi-pi0
  RooRealVar* argedge_part_bb;
  RooRealVar* argpar_part_bb;
  RooArgusBG* argus_mbc_part;

  RooRealVar* mbc0_part;
  RooRealVar* s_mbc_part_ppp;
  RooGaussian* g_mbc_part;
  RooRealVar* fg_mbc_part;

  RooAddPdf* pdf_mbc_part_ppp;

  //////////////
  // part pdf //
  //////////////
  RooProdPdf* pdf_part;
};

#endif // MEPDFPARTIALB_H
