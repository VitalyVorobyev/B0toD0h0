#ifndef PURITYFIT_H
#define PURITYFIT_H

#include "mepdfsignal.h"
#include "mepdfcombinatorial.h"
#include "mepdfpartialb.h"
#include "../MyParams/myparams.h"
//#include "icpvevent.h"

#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>

#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "Roo1DTable.h"
#include "TMath.h"
#include "RooConstVar.h"

const bool m_release_params = true;

typedef struct{
  int exp;
  //int run;
  //int evtn;

  int flv;
  int bin;

  int flv_mc;
  int bin_mc;

  int mode;
  int h0mode;
  int b0f;

  double z_sig;
  double z_asc;
  double sz_sig;
  double sz_asc;
  double chisq_sig;
  double chisq_asc;
  int ntrk_sig;
  int ntrk_asc;
  int ndf_sig;
  int ndf_asc;

  double z_sig_mc;
  double z_asc_mc;

  double f_bkg;
  double f_cont;
  double f_cont_in_comb;
  int sigarea;

  double f_bkg_mc;
  double f_cont_mc;
  double f_cont_in_comb_mc;

  double f_bkg_bin_mc;
  double f_cont_bin_mc;
  double f_cont_in_comb_bin_mc;

  double f_bkg_flv_mc;
  double f_cont_flv_mc;
  double f_cont_in_comb_flv_mc;

  double de;
  double mbc;
  double tag_LH;
  double costhBcms;

  double mp;
  double mm;
  double mp_mc;
  double mm_mc;

  double psig;
  double pcnt;
  double pprt;
  double pcmb;

  double fbb;
  double fbb_err;
  double fprt;
  double fprt_err;
  double Nsig;
  double Nsig_err;

  int NTot;// total number of events in Dalitz-flv bin
  double wrtag;// avetage wrtag for Dalitz bin
  int NTot_bin_mc;
  double wrtag_bin_mc;
  int NTot_flv_mc;
  double wrtag_flv_mc;
  int NTot_mc;
  double wrtag_mc;
} ICPVEv;

class PurityFit{
public:
  PurityFit(const int type, const int _b0f = -1);

  void FitSigPdf(void);
  void FitMbcSigPdf(void);
  void CheckSigPdf(void);

  void FitCombPdf(const vector<int> streams);
  void CheckCombPdf(const vector<int> streams);

  void FitContPdf(const vector<int> streams);
  void CheckContPdf(const vector<int> streams);

  void FitPartPdf(const vector<int> streams);
  void CheckPartPdf(const vector<int> streams);

  void MakeGenMCPutiryFit(const vector<int> streams, const bool fixshape, const int svd, const int icpv_flag, const int bb_or_cnt_flag, const int sig_bkg_flag, const int nstr=1, const int cstr=0);
  void MakeGenMCPutiryFit2(const vector<int> streams, const bool singlefbb, const bool fixshape, const int svd = 0, const int nstr=1, const int cstr=0);
  void ZeroSigTest(const vector<int> streams, const bool fixshape,const int svd = 0){
    MakeGenMCPutiryFit(streams,fixshape,svd,1,2,0);
  }
  void ZeroSigTest2(const vector<int> streams, const bool singlefbb, const bool fixshape,const int svd = 0){
    MakeGenMCPutiryFit2(streams,singlefbb,fixshape,svd,2,0);
  }
  void ZeroBkgTest(const vector<int> streams, const bool fixshape,const int svd,const bool wide_window){
    MakeGenMCPutiryFit(streams,fixshape,svd,1,1,0,wide_window);
  }
  void ZeroBkgTest2(const vector<int> streams, const bool singlefbb, const bool fixshape,const int svd = 0){
    MakeGenMCPutiryFit2(streams,singlefbb,fixshape,svd,1,0);
  }
  void MakeFit(const bool fixshape);
  void MakeFit2(const bool singlefbb, const bool fixshape);

  void ReadSigParamsFromFile(void){  pdf_gen_sig->GetParametersFromFile();}
  void ReadCmbParamsFromFile(void){  pdf_gen_cmb->GetParametersFromFile();}
  void ReadPartParamsFromFile(void){ pdf_gen_part->GetParametersFromFile();}

  void SaveSidebandTree(const int svd = 0);
  void SaveSidebandTree(const vector<int> streams, const int svd = 0, const int bb_or_cnt_flag = 0);
  void SaveSigLineCPVTrees(void);
  void SaveSigCPVTree(const bool line_flag = false,const string& angle = string(""));

  void CheckAllPdfs(const vector<int> streams){
    CheckCombPdf(streams);
    CheckContPdf(streams);
    CheckPartPdf(streams);
  }

  void ReadAllParamsFromFile(void){
    cout << "Read all parameters from files" << endl;
    ReadCmbParamsFromFile();
    ReadSigParamsFromFile();
    ReadPartParamsFromFile();
  }

  RooDataSet* GetGenMCDataSet(const vector<int> streams);
  RooDataSet* GetBBBackMCDataSet(const vector<int> streams);
  RooDataSet* GetGenMCMbcSidebandSet(const vector<int> streams);
  void DefineElliRange(void);
  void CountTrueNumbers(RooDataSet* ds, RooDataSet *cntds);
  void PrintTrueNumbers(void);
  void WriteTrueNumbers(void);

private:
  bool m_fixshape;
  bool m_singlefbb;
  int m_bb_or_cnt_flag;
  int m_sig_bkg_flag;
  int m_svd;
  bool m_mcflag;
  int m_nstr;
  int m_cstr;

  void GetShuffledVector(const int size,vector<int> &vec);
  void WriteBinFlvMap(RooDataSet* ds,const char* fname);
  void FillEvtVec(vector<ICPVEv> &ev_vec,RooDataSet* ds);
  void FillFlvBinsArrays(vector<ICPVEv> (&ev_arr)[2][16], RooDataSet *ds);
  double GetFBkg(int& evnum, double &de, double &mbc, const int flv_ind, const int bin_ind);
  int CalcWTagAndNEv(RooDataSet* ds, double* wtag_arr, int (&nev_arr)[2][16], const bool mcbin, const bool mcflv);
  int MakePredictions(const double* wtag_arr,const int (&nev_arr)[2][16],double (&sig_arr)[2][16],double (&sig_err_arr)[2][16],double (&cmb_arr)[2][16],double (&cnt_arr)[2][16],double (&prt_arr)[2][16],const char* fname);
  int MakePredictions2(const int (&nev_arr)[2][16],double (&sig_arr)[2][16],double (&sig_err_arr)[2][16],double (&cmb_arr)[2][16],double (&cnt_arr)[2][16],double (&prt_arr)[2][16],const char* fname);

  double de_line_size(void);
  double mbc_line_size(void);
  void GenerateDS(const int NTot = 1e6);
  void FillEvent(const RooArgSet *aset, ICPVEv& ev,const int mctype);
  void CopyEvent(const ICPVEv& ev_from, ICPVEv& ev_to);
  void CopyTimeEvent(const ICPVEv& ev_from, ICPVEv& ev_to);
//  void SetArgSet(RooArgSet* aset, const bool include_bdt = true, const int svd = 0, const int mcflag);
  RooArgSet* GetArgSet(const int mctype,const bool include_bdt = true);

  void Fit(RooDataSet* ds);
  void Fit2(RooDataSet* ds);
  void CalcSigIntegrals(void);
  void CalcSigIntegrals2(void);
  void CalcSB1Integrals(void);
  void CalcSB2Integrals(void);
  void PrintIntegrals(void);
  void WriteIntegrals(void);
  void WriteIntegrals2(void);
  void Draw(RooAddPdf& pdf,RooDataSet* ds){
    DrawDeltaE(&pdf,ds);
    DrawMbc(&pdf,ds);
  }
  void DrawData(RooAddPdf& pdf,RooDataSet* ds){
    DrawDeltaE(&pdf,ds);
    DrawMbc(&pdf,ds);
  }
  void DrawDeltaE(RooAddPdf* pdf, RooDataSet* ds);
  void DrawDeltaE2(RooDataSet* ds);
  void DrawMbc(RooAddPdf* pdf, RooDataSet* ds);
  void CalcBinsFractions(RooDataSet* ds);
  void MixCPVTree(RooDataSet *ds);
  void MixCPVTree2(RooDataSet *ds);
  bool IsGoodNumbers(const ICPVEv& ev);
  void SaveCPVTree(RooDataSet* ds);
  TTree* GetCPVTree(ICPVEv& ev);

  RooDataSet* GetSigMCDataSet(void);
  RooDataSet* GetSigMCLineDataSet(const string& angle = string(""));
  RooDataSet* GetContMCDataSet(const vector<int> streams);
  RooDataSet* GetCombMCDataSet(const vector<int> streams, const int svd = 0);
  RooDataSet* GetRawBackMCDataSet(const vector<int> streams, const int svd = 0);
  RooDataSet* GetPartMCDataSet(const vector<int> streams);
  RooDataSet* GetDataSet(const int svd = 0);
  RooDataSet* GetMbcSidebandSet(const int svd = 0);
  RooDataSet* GetBBMCDataSet(const vector<int> streams);
  string GetCut(void);

  RooArgSet* argset;
  RooDataSet* sig_ds_svd1;
  RooDataSet* sig_ds_svd2;
  RooDataSet* raw_back_ds_svd1;
  RooDataSet* raw_back_ds_svd2;

  MyParams* cuts;
  // Variables for selections //
  // discrete variables //
  RooCategory* rndm_pi0;
  RooCategory* b0f;
  RooCategory* d0f;
  RooCategory* h0f;
  RooCategory* mode;
  RooCategory* h0mode;
  RooCategory* flv;
  RooCategory* bin;
  RooCategory* flv_mc;
  RooCategory* bin_mc;
  RooCategory* good_icpv;
//  RooSuperCategory* binflv;

  RooCategory* ndf_asc;
  RooCategory* ndf_rec;
  RooCategory* ntrk_rec;
  RooCategory* ntrk_asc;

  RooCategory* exp;
  RooCategory* exp_svd1;
  RooCategory* exp_svd2;
//  RooCategory* run;
//  RooCategory* evtn;

  RooRealVar* mp;
  RooRealVar* mm;
  RooRealVar* mp_mc;
  RooRealVar* mm_mc;

  RooRealVar* z_rec;
  RooRealVar* sz_rec;
  RooRealVar* chisq_rec;
  RooRealVar* z_rec_mc;

  RooRealVar* z_asc;
  RooRealVar* sz_asc;
  RooRealVar* chisq_asc;
  RooRealVar* z_asc_mc;

  RooRealVar* costhBcms;

  RooRealVar* de;
  RooRealVar* mbc;

  RooRealVar* mbc_center;
  RooRealVar* de_center;
  RooRealVar* mbc_radius;
  RooRealVar* de_radius;

  RooRealVar* md;
  RooRealVar* mk;
  RooRealVar* mh0;
  RooRealVar* mpi0;
  RooRealVar* dmetap;
  RooRealVar* dmdts0;

  RooRealVar* chi2_vtx_d0;

  RooRealVar* bdt;
  RooRealVar* lh0;
  RooRealVar* tag_LH;

  RooRealVar* e_g1;
//  RooRealVar* e_g2;
//  RooRealVar* e_g3;
//  RooRealVar* e_g4;

//  RooRealVar* th_g1;
//  RooRealVar* th_g2;
//  RooRealVar* th_g3;
//  RooRealVar* th_g4;

  RooRealVar* pt_pip;
  RooRealVar* pt_pim;
  RooRealVar* pt_pi1;
  RooRealVar* pt_pi2;
  RooRealVar* p_pi0_h0;
  RooRealVar* cos_hel;

  RooRealVar* r_pip;
  RooRealVar* r_pim;
  RooRealVar* z_pip;
  RooRealVar* z_pim;

  RooRealVar* r_pi1;
  RooRealVar* r_pi2;
  RooRealVar* z_pi1;
  RooRealVar* z_pi2;

  string afretcuts;
  // //

  MEPdfSignal*        pdf_gen_sig;
  MEPdfCombinatorial* pdf_gen_cmb;
  MEPdfPartialB*      pdf_gen_part;
  RooAddPdf* pdf;
  RooRealVar* Nsig;
  RooRealVar* Ncmb;
  RooRealVar* fbb;
  RooRealVar* f_p_f_bbc;
  RooFormulaVar* Npart;


  // Fit v.2 //
  RooConstVar*      wrtagI[2][16];
  RooConstVar*      KI[2][16];
  RooConstVar*      xi;
  RooSuperCategory* binflv;
  RooSimultaneous*  simpdf;
  RooFormulaVar*    NsigI[2][16];
  RooFormulaVar*    NprtI[2][16];
  RooRealVar*       NcmbI[2][16];
//  RooFormulaVar*    f_p_f_bbcI[2][16];
  RooRealVar*       fbbI[2][16];
  RooAddPdf*        pdfI[2][16];
  RooAbsPdf*  m_pdf_sig;
  RooProdPdf* m_pdf_comb_bb;
  RooProdPdf* m_pdf_comb_qq;
  RooProdPdf* m_pdf_part;
  RooAddPdf*  m_pdf_comb[2][16];
  // ////// //

  double nsigEl,  nsig_errEl_total;
  double ncmbEl,  ncmb_errEl_total;
  double ncntEl,  ncnt_errEl_total;
  double npartEl, npart_errEl_total;
  double purityEl,purity_errEl;

  double nsigElI[2][16], nsig_errEl_totalI[2][16];
  double ncmbElI[2][16], ncmb_errEl_totalI[2][16];
  double ncntElI[2][16], ncnt_errEl_totalI[2][16];
  double nprtElI[2][16], nprt_errEl_totalI[2][16];
  double purElI[2][16],  pur_errElI[2][16];

  double nsig_sb1,  nsig_err_sb1_total;
  double ncmb_sb1,  ncmb_err_sb1_total;
  double ncnt_sb1,  ncnt_err_sb1_total;
  double npart_sb1, npart_err_sb1_total;
  double purity_sb1,purity_err_sb1;

  double nsig_sb2,  nsig_err_sb2_total;
  double ncmb_sb2,  ncmb_err_sb2_total;
  double ncnt_sb2,  ncnt_err_sb2_total;
  double npart_sb2, npart_err_sb2_total;
  double purity_sb2,purity_err_sb2;

  double sigint,     cmbint,     partint,     cntint;
  double sigint_sb1, cmbint_sb1, partint_sb1, cntint_sb1;
  double sigint_sb2, cmbint_sb2, partint_sb2, cntint_sb2;

  double NSigPredicted[2][16];
  double NSigPredicted_err[2][16];
  int EventsMap[2][16];
  double WrTagMap[16];
  double NCmbPredicted[2][16];
  double NCntPredicted[2][16];
  double NPrtPredicted[2][16];

  double NSigPredicted_mc[2][16];
  double NSigPredicted_err_mc[2][16];
  int EventsMap_mc[2][16];
  double WrTagMap_mc[16];
  double NCmbPredicted_mc[2][16];
  double NCntPredicted_mc[2][16];
  double NPrtPredicted_mc[2][16];

  double NSigPredicted_bin_mc[2][16];
  double NSigPredicted_err_bin_mc[2][16];
  int EventsMap_bin_mc[2][16];
  double WrTagMap_bin_mc[16];
  double NCmbPredicted_bin_mc[2][16];
  double NCntPredicted_bin_mc[2][16];
  double NPrtPredicted_bin_mc[2][16];

  double NSigPredicted_flv_mc[2][16];
  double NSigPredicted_err_flv_mc[2][16];
  int EventsMap_flv_mc[2][16];
  double WrTagMap_flv_mc[16];
  double NCmbPredicted_flv_mc[2][16];
  double NCntPredicted_flv_mc[2][16];
  double NPrtPredicted_flv_mc[2][16];

  RooDataSet* de_mbc_ds;

  int NSIGNAL, NCOMB, NCONT, NPART;
  int NSIGNAL_ELLI, NCOMB_ELLI, NCONT_ELLI, NPART_ELLI;
  double PURITY, PURITY_SB1, PURITY_SB2;
  int NSIGNAL_SB1, NCOMB_SB1, NCONT_SB1, NPART_SB1;
  int NSIGNAL_SB2, NCOMB_SB2, NCONT_SB2, NPART_SB2;
  int NSIGARR[2][16], NSIGARR_ELLI[2][16];
  int NCMBARR[2][16], NCMBARR_ELLI[2][16];
  int NCNTARR[2][16], NCNTARR_ELLI[2][16];
  int NPRTARR[2][16], NPRTARR_ELLI[2][16];

  int m_mode,m_h0mode;

  double mbc_min;
  double mbc_max;
  double de_min;
  double de_max;

  int get_mode(const int type) const;
  int get_h0mode(const int type) const;
};

#endif // PURITYFIT_H
