#ifndef ICPVEVENT_H
#define ICPVEVENT_H

class ICPVEvent{
public:
  ICPVEvent(char* fname);

private:
  int exp;
  //int run;
  //int evtn;

  int flv;
  int bin;

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
  double de;
  double mbc;
  double tag_LH;
  double costhBcms;

  TTree* tree;
  TFile* file;
};

#endif // ICPVEVENT_H
