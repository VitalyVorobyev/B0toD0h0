#include "cuts.h"

void time_fit(void){
    TFile *ifile_sig = TFile::Open("/home/vitaly/B0toDh0/Tuples/fil_b2dh_sigmc.root");
    TTree *tree = (TTree*)ifile_sig->Get("TEvent");

    RooArgSet argset;
    RooRealVar dz("dz","dz",-0.2,0.2,"mm"); argset.add(dz); argset.add(dz);
    RooRealVar dzerr("dzerr","dzerr",0.000,0.002,"mm"); argset.add(dzerr); argset.add(dzerr);

    RooCategory flv_mc("flv_mc","flv_mc");
    flv_mc.AddCategory("B0",1);
    flv_mc.AddCategory("B0B",-1);
    argset.add(flv_mc);

    RooCategory bin("bin","bin");
    bin.AddCategory("1",1);
    bin.AddCategory("-1",-1);
    argset.add(bin);

    RooDataSet ds("ds","ds",tree,argset);

    //Resolution function
    RooRealVar bias("bias","bias",0,-0.01,0.01); bias.setConstant(kTRUE);
    RooRealVar s1("s1","s1",0.001,0.00001,0.1);
    RooGaussModel g1("g1","g1",dz,bias,s1);

    //B0 decay
    RooRealVar tau("tau","tau",0.001,0.00001,0.1,"cm");
    RooRealVar dgamma("dgamma","dgamma",0.); dgamma.setConstant(kTRUE);
    RooRealVar f0("f0","f0",0.); f0.setConstant(kTRUE);
    RooRealVar f1("f1","f1",0.); f1.setConstant(kTRUE);

    RooRealVar C("C","C",0.); C.setConstant(kTRUE);
    RooRealVar S("S","S",0.); S.setConstant(kTRUE);

    RooRealVar K("K","K",0.);

    RooRealVar sin2beta("sin2beta","sin2beta",0.667,-1.,1.);
    RooRormulaVar cos2beta("cos2beta","cos2beta","1-@0*@0",RooArgSet(sin2beta));

    RooFormulaVar f2("f0","f0",0.); f0.setConstant(kTRUE);
    RooFormulaVar f3("f1","f1",0.); f1.setConstant(kTRUE);

    RooRealVar dm("dm","dm",0,0.,0.1,"GeV");

    RooBDecay bdecay("bdecay","bdecay",dz,tau,dgamma,f0,f1, RooAbsReal& f2, RooAbsReal& f3, RooAbsReal& dm,g1,RooBDecay::DoubleSided)
}
