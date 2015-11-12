#include "cuts.h"
#include "ToyMCParams.h"

void dtfit_sigmc(){
    TChain* tree = new TChain("TEvent");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmc_s1.root");
    tree->Add("/home/vitaly/B0toDh0/Tuples/b2dh_sigmc_s2.root");

    RooArgSet argset;
    RooRealVar mbc("mbc","mbc",mbc_min,mbc_max,"GeV"); argset.add(mbc);
    RooRealVar de("de","#DeltaE",-0.15,0.3,"GeV"); argset.add(de);
    de.setRange("signal",de_min,de_max);
    RooRealVar md("md","md",DMass-md_cut,DMass+md_cut,"GeV"); argset.add(md);
    RooRealVar mk("mk","mk",KMass-mk_cut,KMass+mk_cut,"GeV"); argset.add(mk);
    RooRealVar mpi0("mpi0","mpi0",Pi0Mass-mpi0_cut,Pi0Mass+mpi0_cut,"GeV"); argset.add(mpi0);
    RooRealVar bdtgs("bdtgs","bdtgs",0.98,1.); argset.add(bdtgs);
    RooRealVar atckpi_max("atckpi_max","atckpi_max",0.,atckpi_cut); argset.add(atckpi_max);

    RooDataSet ds("ds","ds",tree,argset);
    ds.Print();
}
