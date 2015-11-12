#ifndef MEPDFSIGNAL_H
#define MEPDFSIGNAL_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooNovosibirsk.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "TTree.h"
#include "RooPrintable.h"
#include "RooDataSet.h"

#include "RooPlot.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooHist.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TAxis.h"

#include "../MyParams/myparams.h"

using namespace std;
using namespace RooFit;

//using namespace Root;

class MEPdfSignal{
public:
  MEPdfSignal(RooRealVar* m_de, RooRealVar* m_mbc, const int mode, const int h0mode, const int b0f = -1);

  RooAbsPdf* GetPdf(void){if(ggflag){ return (RooAbsPdf*)pdf_sig_gg;}
                          else        return (RooAbsPdf*)pdf_sig_ppp;}
  int TryParameters(RooDataSet* ds);
  int FitParameters(RooDataSet *ds);
  int FitTailParameters(RooDataSet *ds);
  int FitPeakParameters(RooDataSet *ds);
  int FitMbcParameters(RooDataSet *ds);

  void FixAll(void)  {ChangeParState(kTRUE,false);}
  void FixTail(void) {ChangeParState(kTRUE,true);}
  void Fix(void)     {ChangeParState(kTRUE,false);}
  void FreeAll(void) {ChangeParState(kFALSE,false);}
  void FreeMbc(void) {ChangeMbcParState(kFALSE,false);}

  void Free_de0(void) {de_offset->setConstant(kFALSE);}
  void Fix_de0(void)  {de_offset->setConstant(kTRUE);}
  void Free_mbc0(void){mbc_offset->setConstant(kFALSE);}
  void Fix_mbc0(void) {mbc_offset->setConstant(kTRUE);}

  void DrawDeltaE(RooDataSet* ds, const int tail);
  void DrawMbc(RooDataSet* ds, const int tail);
  void Draw(RooDataSet* ds, const int tail = 0){
    DrawDeltaE(ds,tail);
    DrawMbc(ds,tail);
    return;
  }

  void PrintParameters(void);
  void WriteParameters(void);
  int GetParametersFromFile(void);

private:
  void InitParams(const int b0f);
  void ChangeParState(const int state_flag, const bool tail_flag);
  void ChangeMbcParState(const int state_flag, const bool tail_flag);
  void ChangeDeltaEParState(const int state_flag, const bool tail_flag);
  double de_line_size(void);
  double mbc_line_size(void);
  MyParams* cuts;

  vector<RooRealVar* > de_param_vec;
  vector<RooRealVar* > mbc_param_vec;

  bool ggflag;
  int m_mode, m_h0mode;
  RooRealVar* de;
  RooRealVar* mbc;

  // Stuff for PDFs //
  // Delta E signal //
  RooRealVar* de_offset;
  RooFormulaVar* de0_sig_offset;
  RooFormulaVar* de0_CBl_sig_offset;
  RooFormulaVar* de0_CBr_sig_offset;
  RooFormulaVar* de0_sig_tail_offset;
  RooFormulaVar* de0_CBl_sig_tail_offset;

  RooRealVar* de0_sig;
  RooRealVar* s_de_sig;
  RooGaussian* g_de_sig;

  RooRealVar* de0_CBl_sig;
  RooRealVar* s_CBl_de_sig;
  RooRealVar* alphal_de_sig;
  RooRealVar* nl_de_sig;

  RooRealVar* de0_CBr_sig;
  RooRealVar* s_CBr_de_sig;
  RooRealVar* alphar_de_sig;
  RooRealVar* nr_de_sig;

  RooCBShape* CBl_de_sig;
  RooCBShape* CBr_de_sig;

  RooRealVar* f_CBl_de_sig;
  RooRealVar* f_CBr_de_sig;

  RooAddPdf* pdf_de_sig;
  // h0 -> pi+pi-pi0
  RooRealVar*  de0_sig_tail;
  RooRealVar*  s_de_sig_tail;
  RooGaussian* g_de_sig_tail;

  RooRealVar* de0_CBl_sig_tail;
  RooRealVar* s_CBl_de_sig_tail;
  RooRealVar* nl_de_sig_tail;
  RooRealVar* alphal_de_sig_tail;
  RooCBShape* CBl_de_sig_tail;

  RooRealVar* fCBl_de_sig_tail;

  RooAddPdf* pdf_de_sig_tail;
  // * //
  // Mbc signal //
  // h0 -> gg
  RooRealVar* mbc_offset;

  RooRealVar* a_s_mbc_sig;
  RooRealVar* b_s_mbc_sig;
  RooRealVar* c_s_mbc_sig;
  RooFormulaVar* S_mbc_sig;

  RooRealVar* alpha_mbc_sig;

  RooRealVar* a_mbc0_sig;
  RooRealVar* b_mbc0_sig;
  RooRealVar* c_mbc0_sig;
  RooRealVar* d_mbc0_sig;
  RooFormulaVar* MBC0_sig;
  RooNovosibirsk* pdf_mbc_sig;
  // h0 -> pi+pi-pi0
  RooRealVar* a_s_mbc_sig_tail;
  RooRealVar* b_s_mbc_sig_tail;
  RooRealVar* c_s_mbc_sig_tail;
  RooFormulaVar* S_mbc_sig_tail;

  RooRealVar* a_mbc0_sig_tail;
  RooRealVar* b_mbc0_sig_tail;
  RooRealVar* c_mbc0_sig_tail;
  RooFormulaVar* MBC0_sig_tail;
  RooNovosibirsk* pdf_mbc_sig_tail;
  // //

  RooProdPdf* pdf_peak;
  RooProdPdf* pdf_tail;
  RooRealVar* f_tail;
  RooAddPdf* pdf_sig_ppp;
  RooProdPdf* pdf_sig_gg;
};

#endif // MEPDFSIGNAL_H
