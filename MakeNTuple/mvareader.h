#ifndef MVAREADER_H
#define MVAREADER_H

#include "icpvevent.h"

class MVAReader
{
public:
    MVAReader(const int mode);
    double BDT(const TMVAEvent& evt);
    double lh(const ICPVEvent& evt);

private:
  rooksfw* ksfw1_pi0;
  rooksfw* ksfw1_etagg;
  rooksfw* ksfw1_etappp;
  rooksfw* ksfw1_omega;
  rooksfw* ksfw0_pi0;
  rooksfw* ksfw0_etagg;
  rooksfw* ksfw0_etappp;
  rooksfw* ksfw0_omega;

  rooksfw* ksfw0_dst0pi0;
  rooksfw* ksfw0_dst0etagg;
  rooksfw* ksfw0_dst0etappp;
  rooksfw* ksfw0_etapgg;
  rooksfw* ksfw0_etapppp;

  rooksfw* ksfw1_dst0pi0;
  rooksfw* ksfw1_dst0etagg;
  rooksfw* ksfw1_dst0etappp;
  rooksfw* ksfw1_etapgg;
  rooksfw* ksfw1_etapppp;

  TMVA::Reader* reader_pi0;
  TMVA::Reader* reader_gg;
  TMVA::Reader* reader_ppp;
  TMVA::Reader* reader_omega;
};

#endif // MVAREADER_H
