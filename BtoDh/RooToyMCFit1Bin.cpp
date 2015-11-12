#include "ToyMCParams.h"

using namespace RooFit;

void RooToyMCFit1Bin(const int _bin){
    if(_bin<1 || _bin>8){
        cout << "Bin number should be between 1 and 8" << endl;
        return;
    }
    stringstream out;
    out.str("");
    out << "toyMC_" << _sigma_over_tau << "_" << _purity << "_" << mistag_rate << ".root";
    TFile* file = TFile::Open(out.str().c_str());
    TTree* tree = (TTree*)file->Get("ToyTree");
    cout << "The tree has been readed from the file " << out.str().c_str() << endl;

    RooRealVar tau("tau","tau",_tau,"ps"); tau.setConstant(kTRUE);
    RooRealVar dm("dm","dm",_dm,"ps^{-1}"); dm.setConstant(kTRUE);
    RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta,-5.,5.); if(constBeta) sin2beta.setConstant(kTRUE);
    RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta,-5.,5.); if(constBeta) cos2beta.setConstant(kTRUE);
    RooRealVar dt("dt","#Deltat",-5.,5.,"ps");
    RooRealVar avgMisgat("avgMisgat","avgMisgat",mistag_rate,0.0,0.5); if(constMistag) avgMisgat.setConstant(kTRUE);
    RooRealVar delMisgat("delMisgat","delMisgat",0); delMisgat.setConstant(kTRUE);
    RooRealVar mu("mu","mu",0); mu.setConstant(kTRUE);
    RooRealVar moment("moment","moment",0.);  moment.setConstant(kTRUE);
    RooRealVar parity("parity","parity",-1.); parity.setConstant(kTRUE);

    cout << "Preparing coefficients..." << endl;
    RooRealVar*    K  = new RooRealVar("K","K",K8[_bin-1],0.,1.); if(constK) K->setConstant(kTRUE);
    RooFormulaVar* Kb = new RooFormulaVar("Kb","Kb","1-@0",RooArgList(*K));

    RooRealVar* C = new RooRealVar("C","C",_C[_bin-1]); C->setConstant(kTRUE);
    RooRealVar* S = new RooRealVar("S","S",_S[_bin-1]); S->setConstant(kTRUE);

    RooFormulaVar* a1  = new RooFormulaVar("a1","a1","-(@0-@1)/(@0+@1)",RooArgList(*K,*Kb));
    RooFormulaVar* a1b = new RooFormulaVar("a1b","a1b","(@0-@1)/(@0+@1)",RooArgList(*K,*Kb));
    RooFormulaVar* a2  = new RooFormulaVar("a2","a2","(@0-@1)/(@0+@1)",RooArgList(*K,*Kb));
    RooFormulaVar* a2b = new RooFormulaVar("a2b","a2b","-(@0-@1)/(@0+@1)",RooArgList(*K,*Kb));

    RooFormulaVar* b1  = new RooFormulaVar("b1","b1","2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1);",RooArgList(*K,*Kb,*C,*S,sin2beta,cos2beta));
    RooFormulaVar* b1b = new RooFormulaVar("b1b","b1b","2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1);",RooArgList(*K,*Kb,*C,*S,sin2beta,cos2beta));
    RooFormulaVar* b2  = new RooFormulaVar("b2","b2","-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1);",RooArgList(*K,*Kb,*C,*S,sin2beta,cos2beta));
    RooFormulaVar* b2b = new RooFormulaVar("b2b","b2b","-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1);",RooArgList(*K,*Kb,*C,*S,sin2beta,cos2beta));

    RooRealVar* dgamma = new RooRealVar("dgamma","dgamma",0.); dgamma->setConstant(kTRUE);
    RooRealVar* f0 = new RooRealVar("f0","f0",1.); f0->setConstant(kTRUE);
    RooRealVar* f1 = new RooRealVar("f1","f1",0.); f1->setConstant(kTRUE);

    RooCategory tag("tag","tag");
    tag.defineType("B0",1);
    tag.defineType("anti-B0",-1);

    RooCategory bin("bin","bin");
    bin.defineType("bin",_bin);
    bin.defineType("binb",-_bin);

    RooSuperCategory bintag("bintag","bintag",RooArgSet(bin,tag));

    RooDataSet d("data","data",tree,RooArgSet(dt,bin,tag));
    cout << "DataSet is ready." << endl;
    d.Print();

    RooRealVar mean("mean","mean",0.,"ps"); mean.setConstant(kTRUE);
    RooRealVar sigma("sigma","sigma",_sigma_over_tau*_tau,0.,_tau,"ps"); if(constSigma) sigma.setConstant(kTRUE);
    RooGaussModel rf("rf","rf",dt,mean,sigma);
//    RooTruthModel rf("rf","rf",dt);
    RooGaussian rfpdf("rfpdf","rfpdf",dt,mean,sigma);

    cout << "Preparing PDFs..." << endl;
    RooBDecay* sigpdf1  = new RooBDecay("sigpdf1","sigpdf1",dt,tau,*dgamma,*f0,*f1,*a1,*b1,dm,rf,RooBDecay::DoubleSided);
    RooBDecay* sigpdf2  = new RooBDecay("sigpdf2","sigpdf2",dt,tau,*dgamma,*f0,*f1,*a2,*b2,dm,rf,RooBDecay::DoubleSided);
    RooBDecay* sigpdf1b  = new RooBDecay("sigpdf1b","sigpdf1b",dt,tau,*dgamma,*f0,*f1,*a1b,*b1b,dm,rf,RooBDecay::DoubleSided);
    RooBDecay* sigpdf2b  = new RooBDecay("sigpdf2b","sigpdf2b",dt,tau,*dgamma,*f0,*f1,*a2b,*b2b,dm,rf,RooBDecay::DoubleSided);

    RooRealVar fsig("fsig","fsigs",_purity,0.,1.);  if(constFSig) fsig.setConstant(kTRUE);
    RooAddPdf* PDF1 = new RooAddPdf("PDF1","PDF1",RooArgList(*sigpdf1,rfpdf),RooArgList(fsig));
    RooAddPdf* PDF2 = new RooAddPdf("PDF2","PDF2",RooArgList(*sigpdf2,rfpdf),RooArgList(fsig));
    RooAddPdf* PDF1b= new RooAddPdf("PDF1b","PDF1b",RooArgList(*sigpdf1b,rfpdf),RooArgList(fsig));
    RooAddPdf* PDF2b= new RooAddPdf("PDF2b","PDF2b",RooArgList(*sigpdf2b,rfpdf),RooArgList(fsig));

    //Adding mistaging
    RooAddPdf* pdf1 = new RooAddPdf("pdf1","pdf1",RooArgList(*PDF2,*PDF1),RooArgList(avgMisgat));
    RooAddPdf* pdf2 = new RooAddPdf("pdf2","pdf2",RooArgList(*PDF1,*PDF2),RooArgList(avgMisgat));
    RooAddPdf* pdf1b= new RooAddPdf("pdf1b","pdf1b",RooArgList(*PDF2b,*PDF1b),RooArgList(avgMisgat));
    RooAddPdf* pdf2b= new RooAddPdf("pdf2b","pdf2b",RooArgList(*PDF1b,*PDF2b),RooArgList(avgMisgat));

    RooSimultaneous pdf("pdf","pdf",bintag);
    pdf.addPdf(*pdf1,"{bin;B0}");
    pdf.addPdf(*pdf2,"{bin;anti-B0}");
    pdf.addPdf(*pdf1b,"{binb;B0}");
    pdf.addPdf(*pdf2b,"{binb;anti-B0}");

    cout << "Fitting..." << endl;
    pdf.fitTo(d,Verbose(),Timer());

    cout << "Drawing plots." << endl;
    // Plus bin
    RooPlot* dtFrame = dt.frame();

    // B0
    out.str("");
    out << "tag == 1 && bin == " << _bin;
    cout << out.str() << endl;
    RooDataSet* ds = d.reduce(out.str().c_str());
    ds->Print();
    ds->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kBlue));
    bintag = "{bin;B0}";
    pdf.plotOn(dtFrame,ProjWData(RooArgSet(),*ds),Slice(bintag),LineColor(kBlue));
    double chi2 = dtFrame->chiSquare();

    // anti-B0
    out.str("");
    out << "tag == -1 && bin == " << _bin;
    cout << out.str() << endl;
    RooDataSet* dsb = d.reduce(out.str().c_str());
    dsb->Print();
    dsb->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kRed));
    bintag = "{bin;anti-B0}";
    pdf.plotOn(dtFrame,ProjWData(RooArgSet(),*dsb),Slice(bintag),LineColor(kRed));
    double chi2b = dtFrame->chiSquare();

    // Canvas
    out.str("");
    out << "#Delta t, toy MC, bin == " << _bin;
    TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
    cm->cd();
    dtFrame->GetXaxis()->SetTitleSize(0.05);
    dtFrame->GetXaxis()->SetTitleOffset(0.85);
    dtFrame->GetXaxis()->SetLabelSize(0.05);
    dtFrame->GetYaxis()->SetTitleOffset(1.6);

    TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.99,"brNDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    out.str("");
    out << "bin = " << _bin;
    pt->AddText(out.str().c_str());
    out.str("");
    out << "#chi^{2}(B^{0}) = " << chi2;
    pt->AddText(out.str().c_str());
    out.str("");
    out << "#chi^{2}(#barB^{0}) = " << chi2b;
    pt->AddText(out.str().c_str());
    out.str("");
    out << "#sigma/#tau = " << sigma.getVal()/tau.getVal();
    pt->AddText(out.str().c_str());
    out.str("");
    out << "purity = " << _purity;
    pt->AddText(out.str().c_str());
    out.str("");
    out << "mistag = " << mistag_rate;
    pt->AddText(out.str().c_str());

    dtFrame->Draw();
    pt->Draw();

    // Minus bin
    RooPlot* dtFrameB = dt.frame();

    // B0
    out.str("");
    out << "tag == 1 && bin == " << -_bin;
    cout << out.str() << endl;
    RooDataSet* ds2 = d.reduce(out.str().c_str());
    ds2->Print();
    ds2->plotOn(dtFrameB,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kBlue));
    bintag = "{binb;B0}";
    pdf.plotOn(dtFrameB,ProjWData(RooArgSet(),*ds2),Slice(bintag),LineColor(kBlue));
    chi2 = dtFrameB->chiSquare();

    //anti-B0
    out.str("");
    out << "tag == -1 && bin == " << -_bin;
    cout << out.str() << endl;
    RooDataSet* dsb2 = d.reduce(out.str().c_str());
    dsb2->Print();
    dsb2->plotOn(dtFrameB,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kRed));
    bintag = "{binb;anti-B0}";
    pdf.plotOn(dtFrameB,ProjWData(RooArgSet(),*dsb2),Slice(bintag),LineColor(kRed));
    chi2b = dtFrameB->chiSquare();

    //Canvas
    out.str("");
    out << "#Delta t, toy MC, bin == " << -_bin;
    TCanvas* cm2 = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
    cm2->cd();
    dtFrameB->GetXaxis()->SetTitleSize(0.05);
    dtFrameB->GetXaxis()->SetTitleOffset(0.85);
    dtFrameB->GetXaxis()->SetLabelSize(0.05);
    dtFrameB->GetYaxis()->SetTitleOffset(1.6);

    TPaveText *ptB = new TPaveText(0.7,0.65,0.98,0.99,"brNDC");
    ptB->SetFillColor(0);
    ptB->SetTextAlign(12);
    out.str("");
    out << "bin = " << -_bin;
    ptB->AddText(out.str().c_str());
    out.str("");
    out << "#chi^{2}(B^{0}) = " << chi2;
    ptB->AddText(out.str().c_str());
    out.str("");
    out << "#chi^{2}(#barB^{0}) = " << chi2b;
    ptB->AddText(out.str().c_str());
    out.str("");
    out << "#sigma/#tau = " << sigma.getVal()/tau.getVal();
    ptB->AddText(out.str().c_str());
    out.str("");
    out << "purity = " << _purity;
    ptB->AddText(out.str().c_str());
    out.str("");
    out << "mistag = " << mistag_rate;
    ptB->AddText(out.str().c_str());

    dtFrameB->Draw();
    ptB->Draw();

    return;
}
