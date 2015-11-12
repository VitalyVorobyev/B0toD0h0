#include "ToyMCParams.h"

void toyMCGen(const int nsig, const float wtag, const float purity){
    TFile* file = new TFile("toyMC.root","RECREATE");

//    gROOT->ProcessLine(".L ../PDFs/Faddeeva.hh+");
    gROOT->LoadMacro("../PDFs/libFadeeva.so");
    gROOT->ProcessLine(".L ../PDFs/RooBtoDhFadeeva.cxx++");

    RooRealVar tau("tau","tau",_tau,"ps"); tau.setConstant(kTRUE);
    RooRealVar dm("dm","dm",_dm,"ps^{-1}"); dm.setConstant(kTRUE);
    RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta); sin2beta.setConstant(kTRUE);
    RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta); cos2beta.setConstant(kTRUE);
    RooRealVar dt("dt","dt",-3.,3.,"ps");

    RooRealVar moment("moment","moment",-1.); moment.setConstant(kTRUE);
    RooRealVar parity("parity","parity",0.); parity.setConstant(kTRUE);

    RooRealVar* K[8];
    RooRealVar* Kb[8];
    RooRealVar* C[8];
    RooRealVar* S[8];
    RooRealVar* Sb[8];
    stringstream out;
    for(int i=0; i<8; i++){
        out.str("");
        out << "K" << i;
        K[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_K[i]); K[i]->setConstant(kTRUE);
        out.str("");
        out << "Kb" << i;
        Kb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_Kb[i]); Kb[i]->setConstant(kTRUE);
        out.str("");
        out << "C" << i;
        C[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_C[i]); C[i]->setConstant(kTRUE);
        out.str("");
        out << "S" << i;
        S[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_S[i]); S[i]->setConstant(kTRUE);
        out.str("");
        out << "Sb" << i;
        Sb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),-_S[i]); Sb[i]->setConstant(kTRUE);
    }
    cout << "Hello world!" << endl;

    RooCategory tag("tag","tag");
    tag.defineType("B0",1);
    tag.defineType("anti-B0",-1);

    RooCategory bin("bin","bin");
    bin.defineType("1",1);
    bin.defineType("2",2);
    bin.defineType("3",3);
    bin.defineType("4",4);
    bin.defineType("5",5);
    bin.defineType("6",6);
    bin.defineType("7",7);
    bin.defineType("8",8);
    bin.defineType("-1",-1);
    bin.defineType("-2",-2);
    bin.defineType("-3",-3);
    bin.defineType("-4",-4);
    bin.defineType("-5",-5);
    bin.defineType("-6",-6);
    bin.defineType("-7",-7);
    bin.defineType("-8",-8);

    RooRealVar mean("mean","mean",0.,"ps"); mean.setConstant(kTRUE);
    RooRealVar sigma("sigma","sigma",_tau,"ps"); sigma.setConstant(kTRUE);
//    RooGaussModel rf("rf","rf",dt,mean,sigma);
//    RooTruthModel rf("rf","rf",dt);

    cout << "Raw PDFs..." << endl;

    RooBtoDhFadeeva pdfs[16];
    RooBtoDhFadeeva pdfsb[16];
    for(int i=0; i<8; i++){
        out.str("");
        out << "B0 pdf " << i+1;
        pdfs[i]   = RooBtoDhFadeeva(out.str().c_str(),out.str().c_str(),dt,tau,dm,sin2beta,cos2beta,*K[i],*Kb[i],*C[i],*S[i],moment,parity);
        out.str("");
        out << "B0 pdf " << -(i+1);
        pdfs[i+8] = RooBtoDhFadeeva(out.str().c_str(),out.str().c_str(),dt,tau,dm,sin2beta,cos2beta,*Kb[i],*K[i],*C[i],*Sb[i],moment,parity);
        out.str("");
        out << "anti-B0 pdf " << i+1;
        pdfsb[i]  = RooBtoDhFadeeva(out.str().c_str(),out.str().c_str(),dt,tau,dm,sin2beta,cos2beta,*K[i],*Kb[i],*C[i],*Sb[i],moment,parity);
        out.str("");
        out << "anti-B0 pdf " << -(i+1);
        pdfsb[i+8]= RooBtoDhFadeeva(out.str().c_str(),out.str().c_str(),dt,tau,dm,sin2beta,cos2beta,*Kb[i],*K[i],*C[i],*S[i],moment,parity);
    }

    cout << "Numerical convolution..." << endl;

//    RooFFTConvPdf PDFs[16];
//    RooFFTConvPdf PDFsb[16];
//    for(int i=0; i<8; i++){
//        out.str("");
//        out << "B0 PDF " << i+1;
//        PDFs[i]   = RooFFTConvPdf(out.str().c_str(),out.str().c_str(),dt,pdfs[i],rf);
//        out.str("");
//        out << "B0 PDF " << -(i+1);
//        PDFs[i+8] = RooFFTConvPdf(out.str().c_str(),out.str().c_str(),dt,pdfs[i+8],rf);
//        out.str("");
//        out << "anti-B0 PDF " << i+1;
//        PDFsb[i]  = RooFFTConvPdf(out.str().c_str(),out.str().c_str(),dt,pdfsb[i],rf);
//        out.str("");
//        out << "anti-B0 PDF " << -(i+1);
//        PDFsb[i+8]= RooFFTConvPdf(out.str().c_str(),out.str().c_str(),dt,pdfsb[i+8],rf);
//    }

    cout << "Filling dataset..." << endl;

    RooDataSet d("data","data",RooArgSet(dt,bin,tag));
    RooDataSet *data;
    int N;
    for(int i=0; i<8; i++){
        out.str("");
        out << i+1;
        bin.setLabel(out.str().c_str());
        // Signal with true tag
        N = (int)(nsig*(1.-wtag)*_K[i]*0.5);
        // B0
        data = PDFs[i].generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = PDFsb[i].generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // Signal with wrong tag
        N = (int)(nsig*wtag*_K[i]*0.5);
        // B0
        data = PDFsb[i].generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = PDFs[i].generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        //Background
        N = (int)(nsig/purity*_K[i]*0.5);
        // B0
        data = rf.generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = rf.generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }
        //////////////
        // anti-bin //
        //////////////
        out.str("");
        out << -(i+1);
        bin.setLabel(out.str().c_str());
        // Signal with true tag
        N = (int)(nsig*(1.-wtag)*_Kb[i]*0.5);
        // B0
        data = PDFs[i+8].generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = PDFsb[i+8].generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // Signal with wrong tag
        N = (int)(nsig*wtag*_Kb[i]*0.5);
        // B0
        data = PDFsb[i+8].generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = PDFs[i+8].generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        //Background
        N = (int)(nsig/purity*_Kb[i]*0.5);
        // B0
        data = rf.generate(RooArgSet(dt),N);
        tag.setLabel("B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }

        // anti-B0
        data = rf.generate(RooArgSet(dt),N);
        tag.setLabel("anti-B0");
        for(int j=0; j<N; j++){
            dt = data->get(j)->find(dt.GetName())->getVal();
            d.add(RooArgSet(dt,bin,tag));
        }
    }

    d.Print();
    d.Write();
    file->Close();

    return;
}
