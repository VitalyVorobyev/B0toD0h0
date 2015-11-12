//#include "ToyMCParams.h"
#include <unistd.h>

using namespace RooFit;
using namespace std;

const double cm2ps = 78.4857;
const double deg2rad = TMath::Pi()/180.;
const double beta = 67.;
const double _btau = 0.0454193;//1.533*2.99792458/cm2ps;// ps
const double _dm   = 0.507/2.99792458/2.99792458*cm2ps;// ps^{-1}
//const double _btau     = 1.533/cm2ps;// ps
//const double _dm       = 0.507/2.99792458*cm2ps;// ps^{-1}//0.507
const double _sin2beta = TMath::Sin(2.*beta*deg2rad);
const double _cos2beta = TMath::Cos(2.*beta*deg2rad);
int xi;

const bool fix_tau      = false;
const bool fix_dm       = false;
const bool fix_sin2beta = false;
const bool fix_cos2beta = false;
const bool cK           = true;

const double dtmax = 30.;
const double dtmin = -dtmax;
const double dzmax = dtmax/cm2ps;
const double dzmin = -dzmax;

RooRealVar* dz   = new RooRealVar("dt_mc","#Deltaz",dzmin,dzmax,"mm");
RooCategory* tag = new RooCategory("flv_mc","flv_mc");
RooCategory* bin = new RooCategory("bin_mc","bin_mc");
RooDataSet* d;

//int bins[]

class FR{
public:
  FR(): sin2beta(0), sin2betaErr(0), cos2beta(0), cos2betaErr(0), tau(0), tauErr(0), dm(0), dmErr(0), cS(false), cC(false), cTau(false), cDm(false) {}
  Print(){
    if(!cS)   cout << "sin2beta = " << sin2beta << " +- " << sin2betaErr << endl;
    if(!cC)   cout << "cos2beta = " << cos2beta << " +- " << cos2betaErr << endl;
    if(!cTau) cout << "tau      = " << tau      << " +- " << tauErr      << endl;
    if(!cDm)  cout << "dm       = " << dm       << " +- " << dmErr       << endl;
  }
  void SetSin(const double& x){sin2beta = x;}
  void SetSinErr(const double& x){sin2betaErr = x;}
  void SetSin(const double& x,const double& y){sin2beta = x; sin2betaErr = y;}

  void SetCos(const double& x){cos2beta = x;}
  void SetCos(const double& x,const double& y){cos2beta = x; cos2betaErr = y;}
  void SetCosErr(const double& x){cos2betaErr = x;}

  void SetTau(const double& x){tau = x;}
  void SetTau(const double& x,const double& y){tau = x; tauErr = y;}
  void SetTauErr(const double& x){tauErr = x;}

  void SetDm(const double& x){dm = x;}
  void SetDm(const double& x,const double& y){dm = x; dmErr = y;}
  void SetDmErr(const double& x){dmErr = x;}

  double Sin(void)    const {return sin2beta;}
  double SinErr(void) const {return sin2betaErr;}
  double Cos(void)    const {return cos2beta;}
  double CosErr(void) const {return cos2betaErr;}
  double Tau(void)    const {return tau;}
  double TauErr(void) const {return tauErr;}
  double dm(void)     const {return dm;}
  double dmErr(void)  const {return dmErr;}

  FR operator=(const FR& x){
    this->sin2beta    = x.Sin();
    this->sin2betaErr = x.SinErr();
    this->cos2beta    = x.Cos();
    this->cos2betaErr = x.CosErr();
    this->tau         = x.Tau();
    this->tauErr      = x.TauErr();
    this->dm          = x.dm();
    this->dmErr       = x.dmErr();
    return *this;
  }

private:
  double sin2beta, sin2betaErr;
  double cos2beta, cos2betaErr;
  double tau, tauErr;
  double dm, dmErr;
  bool cS,cC,cTau,cDm;
};

double Carr[8],Sarr[8],Karr[8],Kbarr[8];

//Model integrals
const double Carr_model[8] = { 0.594511, -0.329868, -0.61714, -0.758994,-0.368921,  0.133975, 0.484073, 0.742303};
const double Sarr_model[8] = {-0.357055, -0.680259, -0.277483,-0.259506, 0.440863,  0.623397, 0.358263,-0.0558073};
const double Karr_model[8] = { 0.0665352, 0.0835056, 0.0360653,0.0961057,0.0722859, 0.0980402,0.123951, 0.198887};
const double Kbarr_model[8]= { 0.0304473, 0.00571654,0.0022683,0.0264409,0.00969127,0.0105694,0.0311109,0.10838};

// rho Ks
const double Carr_rhoKs[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
const double Sarr_rhoKs[8] = {0,0,0,0,0,0,0,0};
const double Karr_rhoKs[8] = {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};
const double Kbarr_rhoKs[8]= {0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625};

void init_arrs(const int mode){
  if(mode == 5){
     for(int i=0; i<8; i++){
       Carr[i]  = Carr_rhoKs[i];
       Sarr[i]  = Sarr_rhoKs[i];
       Karr[i]  = Karr_rhoKs[i];
       Kbarr[i] = Kbarr_rhoKs[i];
     }
  } else{
      for(int i=0; i<8; i++){
        Carr[i]  = Carr_model[i];
        Sarr[i]  = Sarr_model[i];
        Karr[i]  = Karr_model[i];
        Kbarr[i] = Kbarr_model[i];
      }
  }
}

FR RooToyMCFit(const int mode){
  stringstream out;
  FR FitRes;

  const bool dump = true;
  RooRealVar tau("tau","tau",_btau,0.,10.*_btau,"mm");
  RooRealVar dm("dm","dm",_dm,0.,5.*_dm,"mm^{-1}");
  RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta,-3.,3.);
  RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta,-3.,3.);

  if(fix_tau)      tau.setConstant(kTRUE);
  if(fix_dm)       dm.setConstant(kTRUE);
  if(fix_sin2beta) sin2beta.setConstant(kTRUE);
  if(fix_cos2beta) cos2beta.setConstant(kTRUE);

//  RooRealVar dz("dz_mc","#Deltaz",-10./cm2ps,10./cm2ps,"mm");
  RooRealVar avgMistag("avgMistag","avgMistag",0); avgMistag.setConstant(kTRUE);
  RooRealVar delMistag("delMistag","delMistag",0); delMistag.setConstant(kTRUE);
  RooRealVar mu("mu","mu",0);                      mu.setConstant(kTRUE);

  RooRealVar moment("moment","moment",0.);  moment.setConstant(kTRUE);
  RooRealVar parity("parity","parity",xi);  parity.setConstant(kTRUE);

  RooRealVar* K[8];
  RooAbsReal* Kb[8];
  RooArgList Kset;
  RooRealVar* C[8];
  RooRealVar* S[8];
  RooFormulaVar* a1[2][8];
  RooFormulaVar* b1[2][8];
  RooFormulaVar* a2[2][8];
  RooFormulaVar* b2[2][8];

  for(int i=0; i<8; i++){
    out.str("");
    out << "K" << i+1;
    K[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),Karr[i],0.,1.);
    Kset.add(*K[i]);
    if(cK) K[i]->setConstant(kTRUE);
    if(i!=7){
      out.str("");
      out << "Kb" << i+1;
      Kb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),Kbarr[i],0.,1.);
      Kset.add(*Kb[i]);
      if(cK) ((RooRealVar*)Kb[i])->setConstant(kTRUE);
    }
    out.str("");
    out << "C" << i+1;
    C[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),Carr[i]);
    C[i]->setConstant(kTRUE);
    out.str("");
    out << "S" << i+1;
    S[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),Sarr[i]);
    S[i]->setConstant(kTRUE);
  }
  Kb[7] = new RooFormulaVar("K8b","K8b","1-@0-@1-@2-@3-@4-@5-@6-@7-@8-@9-@10-@11-@12-@13-@14",Kset);

//  const bool psioff = false;
  for(int i=0; i<8; i++){
    out.str("");
    out << "a10_" << i+1;
//    if(!psioff) a1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1*@2*@2)/(@0+@1*@2*@2)",RooArgList(*K[i],*Kb[i],*psi));
//    else
    a1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
    if(dump) a1[0][i]->Print();
    out.str("");
    out << "a11_" << i+1;
//    if(!psioff) a1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0*@2*@2-@1)/(@0*@2*@2+@1)",RooArgList(*K[i],*Kb[i],*psi));
//    else
    a1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
    if(dump) a1[1][i]->Print();
    out.str("");
    out << "a20_" << i+1;
//    if(!psioff) a2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0*@2*@2-@1)/(@0*@2*@2+@1)",RooArgList(*K[i],*Kb[i],*psi));
//    else
    a2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
    if(dump) a2[0][i]->Print();
    out.str("");
    out << "a21_" << i+1;
//    if(!psioff) a2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1*@2*@2)/(@0+@1*@2*@2)",RooArgList(*K[i],*Kb[i],*psi));
//    else
    a2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
    if(dump) a2[1][i]->Print();

    out.str("");
    out << "b10_" << i+1;
//    if(!psioff) b1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1*@6*@6)/(@0+@1*@6*@6)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],*sin2beta,*cos2beta,*psi));
//    else
    b1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*@6*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,parity));
    if(dump) b1[0][i]->Print();
    out.str("");
    out << "b11_" << i+1;
//    if(!psioff) b1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1*@6*@6)/(@0*@6*@6+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],*sin2beta,*cos2beta,*psi));
//    else
    b1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*@6*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,parity));
    if(dump) b1[1][i]->Print();
    out.str("");
    out << "b20_" << i+1;
//    if(!psioff) b2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1*@6*@6)/(@0*@6*@6+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],*sin2beta,*cos2beta,*psi));
//    else
    b2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*@6*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,parity));
    if(dump) b2[0][i]->Print();
    out.str("");
    out << "b21_" << i+1;
//    if(!psioff) b2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1*@6*@6)/(@0+@1*@6*@6)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],*sin2beta,*cos2beta,*psi));
//    else
    b2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*@6*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,parity));
    if(dump) b2[1][i]->Print();

    if(dump){
      cout << "bin = " << i+1 << endl;
      cout << "  a11 = " << a1[0][i]->getVal() << " b11 = " << b1[0][i]->getVal() << endl;
      cout << "  a12 = " << a1[1][i]->getVal() << " b12 = " << b1[1][i]->getVal() << endl;
      cout << "  a21 = " << a2[0][i]->getVal() << " b21 = " << b2[0][i]->getVal() << endl;
      cout << "  a22 = " << a2[1][i]->getVal() << " b22 = " << b2[1][i]->getVal() << endl;
    }
  }

  RooRealVar* dgamma = new RooRealVar("dgamma","dgamma",0.); dgamma->setConstant(kTRUE);
  RooRealVar* f0 = new RooRealVar("f0","f0",1.);             f0->setConstant(kTRUE);
  RooRealVar* f1 = new RooRealVar("f1","f1",0.);             f1->setConstant(kTRUE);

//  RooCategory tag("flv_mc","flv_mc");
//  tag.defineType("B0",1);
//  tag.defineType("anti-B0",-1);

//  RooCategory bin("bin_mc","bin_mc");
//  bin.defineType("1",1); bin.defineType("-1",-1);
//  bin.defineType("2",2); bin.defineType("-2",-2);
//  bin.defineType("3",3); bin.defineType("-3",-3);
//  bin.defineType("4",4); bin.defineType("-4",-4);
//  bin.defineType("5",5); bin.defineType("-5",-5);
//  bin.defineType("6",6); bin.defineType("-6",-6);
//  bin.defineType("7",7); bin.defineType("-7",-7);
//  bin.defineType("8",8); bin.defineType("-8",-8);

  RooSuperCategory bintag("bintag","bintag",RooArgSet(*bin,*tag));

//  RooRealVar mean("mean","mean",0.,"cm"); mean.setConstant(kTRUE);
//  RooRealVar sigma("sigma","sigma",_btau/cm2ps*0.1,"cm"); sigma.setConstant(kTRUE);
//  RooGaussModel rf("rf","rf",*dz,mean,sigma);
  RooTruthModel deltf("deltf","deltf",*dz);
//  RooGaussian rfpdf("rfpdf","rfpdf",dt,mean,sigma);

//  cout << "Preparing PDFs..." << endl;
  RooBDecay* sigpdfs1[8];
  RooBDecay* sigpdfs1b[8];
  RooBDecay* sigpdfs2[8];
  RooBDecay* sigpdfs2b[8];
//  RooAddPdf* PDFs1[8];
//  RooAddPdf* PDFs1b[8];
//  RooAddPdf* PDFs2[8];
//  RooAddPdf* PDFs2b[8];
//  RooAddPdf* pdfs1[8];
//  RooAddPdf* pdfs1b[8];
//  RooAddPdf* pdfs2[8];
//  RooAddPdf* pdfs2b[8];
  RooSimultaneous pdf("pdf","pdf",bintag);
//  RooRealVar fsig("fsig","fsig",PUR,0.,1.); if(cFS) fsig.setConstant(kTRUE);
  for(int j=0; j<8; j++){
    out.str("");
    out << "sigpdf1" << j+1;
    sigpdfs1[j]  = new RooBDecay(out.str().c_str(),out.str().c_str(),*dz,tau,*dgamma,*f0,*f1,*a1[0][j],*b1[0][j],dm,deltf,RooBDecay::DoubleSided);
    out.str("");
    out << "sigpdf1b" << j+1;
    sigpdfs1b[j] = new RooBDecay(out.str().c_str(),out.str().c_str(),*dz,tau,*dgamma,*f0,*f1,*a1[1][j],*b1[1][j],dm,deltf,RooBDecay::DoubleSided);
    out.str("");
    out << "sigpdf2" << j+1;
    sigpdfs2[j]  = new RooBDecay(out.str().c_str(),out.str().c_str(),*dz,tau,*dgamma,*f0,*f1,*a2[0][j],*b2[0][j],dm,deltf,RooBDecay::DoubleSided);
    out.str("");
    out << "sigpdf2b" << j+1;
    sigpdfs2b[j] = new RooBDecay(out.str().c_str(),out.str().c_str(),*dz,tau,*dgamma,*f0,*f1,*a2[1][j],*b2[1][j],dm,deltf,RooBDecay::DoubleSided);

//    out.str("");
//    out << "PDF1" << j+1;
//    PDFs1[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs1[j],rfpdf),RooArgList(fsig));
//    out.str("");
//    out << "PDF1b" << j+1;
//    PDFs1b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs1b[j],rfpdf),RooArgList(fsig));
//    out.str("");
//    out << "PDF2" << j+1;
//    PDFs2[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs2[j],rfpdf),RooArgList(fsig));
//    out.str("");
//    out << "PDF2b" << j+1;
//    PDFs2b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs2b[j],rfpdf),RooArgList(fsig));

    //Adding mistaging
//    out.str("");
//    out << "pdf1" << j+1;
//    pdfs1[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2[j],*PDFs1[j]),RooArgList(avgMistag));
//    out.str("");
//    out << "pdf1b" << j+1;
//    pdfs1b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2b[j],*PDFs1b[j]),RooArgList(avgMistag));
//    out.str("");
//    out << "pdf2" << j+1;
//    pdfs2[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1[j],*PDFs2[j]),RooArgList(avgMistag));
//    out.str("");
//    out << "pdf2b" << j+1;
//    pdfs2b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1b[j],*PDFs2b[j]),RooArgList(avgMistag));

    out.str("");
    out << "{" << j+1 << ";B0}";
    pdf.addPdf(*sigpdfs1[j],out.str().c_str());
    out.str("");
    out << "{" << -(j+1) << ";B0}";
    pdf.addPdf(*sigpdfs1b[j],out.str().c_str());
    out.str("");
    out << "{" << j+1 << ";anti-B0}";
    pdf.addPdf(*sigpdfs2[j],out.str().c_str());
    out.str("");
    out << "{" << -(j+1) << ";anti-B0}";
    pdf.addPdf(*sigpdfs2b[j],out.str().c_str());
  }

//  RooRealVar mean("mean","mean",0.,"ps");
//  RooRealVar sigma("sigma","sigma",0.002,0.,1.5,"cm");
//  RooGaussModel pdf1("pdf1","pdf1",*dz,mean,sigma);

  pdf.fitTo(*d,SumW2Error(kTRUE),PrintLevel(1),Verbose(kTRUE));

  FitRes.SetCos(cos2beta.getVal(),cos2beta.getError());
  FitRes.SetSin(sin2beta.getVal(),sin2beta.getError());
  FitRes.SetTau(tau.getVal(),tau.getError());
  FitRes.SetDm(dm.getVal(),dm.getError());
//  RooDataSet* reducedDS;

  if(mode == 5){
    RooDataSet* reducedDS = (RooDataSet*)d->reduce("flv_mc == 1");
    out.str("");
    out << "dz_B0";
    draw_plots(reducedDS,sigpdfs1[3],out.str());
    RooDataSet* reducedDSb = (RooDataSet*)d->reduce("flv_mc == -1");
    out.str("");
    out << "dz_B0b";
    draw_plots(reducedDSb,sigpdfs2[3],out.str());
  } else{
    for(int i=1; i<8; i++){
      out.str("");
      out << "flv_mc == 1 && bin_mc == " << i+1;
      RooDataSet* reducedDS1 = (RooDataSet*)d->reduce(out.str().c_str());
      reducedDS1->Print();
      out.str("");
      out << "dz_B0_bin" << i+1;
      draw_plots(reducedDS1,sigpdfs1[i],out.str());

      out.str("");
      out << "flv_mc == 1 && bin_mc == " << -(i+1);
      RooDataSet* reducedDS2 = (RooDataSet*)d->reduce(out.str().c_str());
      reducedDS2->Print();
      out.str("");
      out << "dz_B0_bin" << -(i+1);
      draw_plots(reducedDS2,sigpdfs1b[i],out.str());

      out.str("");
      out << "flv_mc == -1 && bin_mc == " << i+1;
      RooDataSet* reducedDS3 = (RooDataSet*)d->reduce(out.str().c_str());
      reducedDS3->Print();
      out.str("");
      out << "dz_B0b_bin" << i+1;
      draw_plots(reducedDS3,sigpdfs2[i],out.str());

      out.str("");
      out << "flv_mc == -1 && bin_mc == " << -(i+1);
      RooDataSet* reducedDS4 = (RooDataSet*)d->reduce(out.str().c_str());
      reducedDS4->Print();
      out.str("");
      out << "dz_B0b_bin" << -(i+1);
      draw_plots(reducedDS4,sigpdfs2b[i],out.str());
    }
  }

  return FitRes;
}

int draw_plots(RooDataSet* ds, RooAbsPdf* pdf,string& fname){
    /////////////
    //  Plots  //
    /////////////
    stringstream out;
    // dz //
    RooPlot* dzFrame = dz->frame();
    ds->plotOn(dzFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    pdf->plotOn(dzFrame,LineWidth(2));

    RooHist* hdepull = dzFrame->pullHist();
    out.str("");
    out << fname << "_pull";
    RooPlot* dzPull = dz->frame(Title(out.str().c_str()));
    dzPull->addPlotable(hdepull,"P");
    dzPull->GetYaxis()->SetRangeUser(-5,5);

    out.str("");
    out << fname << "_canv";
    TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
    cm->cd();

    out.str("");
    out << fname << "_pad3";
    TPad *pad3 = new TPad(out.str().c_str(),out.str().c_str(),0.01,0.20,0.99,0.99);
    out.str("");
    out << fname << "_pad4";
    TPad *pad4 = new TPad(out.str().c_str(),out.str().c_str(),0.01,0.01,0.99,0.20);
    pad3->Draw();
    pad4->Draw();

    pad3->cd();
    pad3->SetLeftMargin(0.15);
    pad3->SetFillColor(0);

    dzFrame->GetXaxis()->SetTitleSize(0.05);
    dzFrame->GetXaxis()->SetTitleOffset(0.85);
    dzFrame->GetXaxis()->SetLabelSize(0.04);
    dzFrame->GetYaxis()->SetTitleOffset(1.6);
    dzFrame->Draw();

    TPaveText *pt = new TPaveText(0.6,0.8,0.98,0.9,"brNDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    out.str("");
    out << "#chi^{2}/n.d.f = " << dzFrame->chiSquare();
    pt->AddText(out.str().c_str());
    pt->Draw();

    pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
    dzPull->SetMarkerSize(0.05); dzPull->Draw();
    TLine *dz_lineUP = new TLine(dzmin,3,dzmax,3);
    dz_lineUP->SetLineColor(kBlue);
    dz_lineUP->SetLineStyle(2);
    dz_lineUP->Draw();
    TLine *dz_line = new TLine(dzmin,0,dzmax,0);
    dz_line->SetLineColor(kBlue);
    dz_line->SetLineStyle(1);
    dz_line->SetLineWidth((Width_t)2.);
    dz_line->Draw();
    TLine *dz_lineDOWN = new TLine(dzmin,-3,dzmax,-3);
    dz_lineDOWN->SetLineColor(kBlue);
    dz_lineDOWN->SetLineStyle(2);
    dz_lineDOWN->Draw();

    cm->Update();
    out.str("");
    out << "pics/" << fname << ".root";
    cm->Print(out.str().c_str());

    out.str("");
    out << "pics/" << fname << ".png";
    cm->Print(out.str().c_str());

//    sleep(0.1);
    string cmd = string("display ") + out.str() + string(" &");
    system(cmd.c_str());
//    sleep(0.1);

    delete cm;

    return 0;
}

int RooCPGenFit(const int _mode = 5){
    string ifilename;
    switch(_mode){
    case 1:
      cout << "Mode: pi0" << endl;
      ifilename = string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s5_full.root");
      xi = 1;
      break;
    case 2:
      cout << "Mode: eta -> gg" << endl;
      ifilename = string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
      xi = 1;
      break;
    case 3:
      cout << "Mode: eta -> pi+pi-pi0" << endl;
      ifilename = string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
      xi = 1;
      break;
    case 4:
      cout << "Mode: omega" << endl;
      ifilename = string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s2_full.root");
      xi = -1;
      break;
    case 5:
      cout << "Mode: omega (rho Ks0)" << endl;
      ifilename = string("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s3_full.root");
      xi = -1;
      break;
    default:
      cout << "Wrong mode" << endl;
      return -1;
    }
    ifile = TFile::Open(ifilename.c_str());
//    stringstream out;

    if(!ifile->IsOpen()){
      cout << "Can't open file " << ifile << endl;
      return -1;
    }
    TTree* tree = (TTree*)ifile->Get("TEvent");
    if(tree == NULL){
      cout << "Can't find TTree TEvent in the file " << ifile << endl;
      return -2;
    }
    init_arrs(_mode);

    tag->defineType("B0",1);
    tag->defineType("anti-B0",-1);

    bin->defineType("1",1); bin->defineType("-1",-1);
    bin->defineType("2",2); bin->defineType("-2",-2);
    bin->defineType("3",3); bin->defineType("-3",-3);
    bin->defineType("4",4); bin->defineType("-4",-4);
    bin->defineType("5",5); bin->defineType("-5",-5);
    bin->defineType("6",6); bin->defineType("-6",-6);
    bin->defineType("7",7); bin->defineType("-7",-7);
    bin->defineType("8",8); bin->defineType("-8",-8);

    d = new RooDataSet("data","data",tree,RooArgSet(*dz,*bin,*tag),"(dt_mc>0 || dt_mc<0)");
//    d = new RooDataSet("data","data",tree,RooArgSet(*dz),"(dz_mc>0 || dz_mc<0)");
//    d = (RooDataSet*)d->reduce("flv_mc == -1");
    d->Print();

//    RooRealVar mean("mean","mean",0.0,-0.1,0.1,"cm");
//    RooRealVar sigma("sigma","sigma",0.02,0.,1.5,"cm");
//    RooGaussian pdf1("pdf1","pdf1",*dz,mean,sigma);

//    pdf1.fitTo(*d);
    FR FitRes = RooToyMCFit(_mode);
    FitRes.Print();

    ifile->Close();

    return 0;
}

