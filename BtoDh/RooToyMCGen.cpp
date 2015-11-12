#include "ToyMCParams.h"

using namespace RooFit;

void RooToyMCGen(const int nsig){
    stringstream out;
    out.str("");
    out << "toyMC_" << _sigma_over_tau << "_" << _purity << "_" << mistag_rate << ".root";
    TFile* file = new TFile(out.str().c_str(),"RECREATE");
    TTree* tree = new TTree("ToyTree","ToyTree");
    double _dt;
    int _bin,_tag;
    tree->Branch("dt",&_dt,"dt/D");
    tree->Branch("bin",&_bin,"bin/I");
    tree->Branch("tag",&_tag,"tag/I");
    TRandom3 rndm(0);

    RooRealVar tau("tau","tau",_tau,"ps"); tau.setConstant(kTRUE);
    RooRealVar dm("dm","dm",_dm,"ps^{-1}"); dm.setConstant(kTRUE);
    RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta); sin2beta.setConstant(kTRUE);
    RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta); cos2beta.setConstant(kTRUE);
    RooRealVar dt("dt","#Deltat",-5.,5.,"ps");
    RooRealVar avgMisgat("avgMisgat","avgMisgat",mistag_rate); avgMisgat.setConstant(kTRUE);
    RooRealVar delMisgat("delMisgat","delMisgat",0); delMisgat.setConstant(kTRUE);
    RooRealVar mu("mu","mu",0); mu.setConstant(kTRUE);

    RooRealVar moment("moment","moment",0.); moment.setConstant(kTRUE);
    RooRealVar parity("parity","parity",-1.); parity.setConstant(kTRUE);

    RooRealVar a1("a1","a1",-1000.,1000.);
    RooRealVar b1("b1","b1",-1000.,1000.);
    RooFormulaVar a2("a2","a2","-@0",RooArgList(a1));
    RooFormulaVar b2("b2","b2","-@0",RooArgList(b1));

    RooRealVar* K[8];
    RooRealVar* Kb[8];
    RooRealVar* C[8];
    RooRealVar* S[8];
    for(int i=0; i<8; i++){
        out.str("");
        out << "K" << i;
        K[i]  = new RooRealVar(out.str().c_str(),out.str().c_str(),_K[i]);  K[i]->setConstant(kTRUE);
        out.str("");
        out << "Kb" << i;
        Kb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_Kb[i]); Kb[i]->setConstant(kTRUE);
        out.str("");
        out << "C" << i;
        C[i]  = new RooRealVar(out.str().c_str(),out.str().c_str(),_C[i]);  C[i]->setConstant(kTRUE);
        out.str("");
        out << "S" << i;
        S[i]  = new RooRealVar(out.str().c_str(),out.str().c_str(),_S[i]);  S[i]->setConstant(kTRUE);
    }

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
    RooRealVar sigma("sigma","sigma",_sigma_over_tau*_tau,"ps"); sigma.setConstant(kTRUE);
    RooGaussModel rf("rf","rf",dt,mean,sigma);

    RooRealVar dgamma("dgamma","dgamma",0.); dgamma.setConstant(kTRUE);
    RooRealVar f0("f0","f0",1.); f0.setConstant(kTRUE);
    RooRealVar f1("f1","f1",0.); f1.setConstant(kTRUE);

//    RooBCPGenDecay pdf("pdf","pdf",dt,tag,tau,dm,avgMisgat,a,b,delMisgat,mu,rf,RooBCPGenDecay::DoubleSided);
    RooBDecay pdf1("sigpdf1","sigpdf1",dt,tau,dgamma,f0,f1,a1,b1,dm,rf,RooBDecay::DoubleSided);
    RooBDecay pdf2("sigpdf2","sigpdf2",dt,tau,dgamma,f0,f1,a2,b2,dm,rf,RooBDecay::DoubleSided);
    RooAddPdf pdf("pdf","pdf",RooArgList(pdf2,pdf1),RooArgList(avgMisgat));
    RooRandom::randomGenerator()->SetSeed(000);

    RooDataSet d("data","data",RooArgSet(dt,bin,tag));
    RooDataSet *data;
    int N;
    for(int i=0; i<8; i++){
        out.str("");
        out << i+1;
        bin.setLabel(out.str().c_str());
        _bin = i+1;
        // B0
        if(PoissonFlag) N = rndm.Poisson(nsig*_K[i]*0.5);
        else            N = (int)(nsig*_K[i]*0.5);
        a1 = -(_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
        b1 = 2.*(_C[i]*_sin2beta+_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
//        cout << "K = " << _K[i] << "Kb = " << _Kb[i] << endl;
//        cout << "C = " << _C[i] << "S = " << _S[i] << endl;
//        cout << "sin = " << _sin2beta << " cos = " << _cos2beta << endl;
//        cout << "bin = " << i+1 << " a = " << a.getVal() << " b = " << b.getVal() << " N = " << N << endl;

        tag.setLabel("B0"); _tag = 1;
        data = pdf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        // anti-B0
        a1 = (_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
        b1 = -2.*(_C[i]*_sin2beta+_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
        tag.setLabel("anti-B0"); _tag = -1;
        if(PoissonFlag) N = rndm.Poisson(nsig*_K[i]*0.5);
        data = pdf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        //Background
        if(PoissonFlag) N = rndm.Poisson((nsig/_purity-nsig)*_K[i]*0.5);
        else            N = (int)((nsig/_purity-nsig)*_K[i]*0.5);
        // B0
        tag.setLabel("B0"); _tag = 1;
        data = rf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        // anti-B0
        tag.setLabel("anti-B0"); _tag = -1;
        if(PoissonFlag) N = rndm.Poisson((nsig/_purity-nsig)*_K[i]*0.5);
        data = rf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }
        //////////////
        // anti-bin //
        //////////////
        out.str("");
        out << -(i+1);
        bin.setLabel(out.str().c_str());
        _bin = -(i+1);
        if(PoissonFlag) N = rndm.Poisson(nsig*_Kb[i]*0.5);
        else            N = (int)(nsig*_Kb[i]*0.5);
        // B0
        a1 = (_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
        b1 = 2.*(_C[i]*_sin2beta-_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
//        cout << "bin = " << -(i+1) << " a = " << a.getVal() << " b = " << b.getVal() << " N = " << N << endl;
        tag.setLabel("B0"); _tag = 1;
        data = pdf.generate(RooArgSet(dt),N);
//        cout << "N = " << N << ", a = " << a << ", b = " << b << endl;
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        // anti-B0
        a1 = -(_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
        b1 = -2.*(_C[i]*_sin2beta-_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
        tag.setLabel("anti-B0"); _tag = -1;
        if(PoissonFlag) N = rndm.Poisson(nsig*_Kb[i]*0.5);
        data = pdf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        //Background
        if(PoissonFlag) N = rndm.Poisson((nsig/_purity-nsig)*_Kb[i]*0.5);
        else            N = (int)((nsig/_purity-nsig)*_Kb[i]*0.5);
        // B0
        tag.setLabel("B0"); _tag = 1;
        data = rf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }

        // anti-B0
        tag.setLabel("anti-B0"); _tag = -1;
        if(PoissonFlag) N = rndm.Poisson((nsig/_purity-nsig)*_Kb[i]*0.5);
        data = rf.generate(RooArgSet(dt),N);
        for(int j=0; j<N; j++){
            dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal(); _dt = dt.getVal();
            d.add(RooArgSet(dt,bin,tag)); tree->Fill();
        }
    }

    tree->Print();
    tree->Write();
    file->Close();
    dt.setBins(50);

    /////////////
    //  Plots  //
    /////////////

    for(int i=-8; i<=8; i++){if(i){
            out.str("");
            out << "tag == 1 && bin == " << i;
            RooDataSet* ds1 = (RooDataSet*)d.reduce(SelectVars(dt),Cut(out.str().c_str()));
            out.str("");
            out << "tag == -1 && bin == " << i;
            RooDataSet* ds2 = (RooDataSet*)d.reduce(SelectVars(dt),Cut(out.str().c_str()));
            RooPlot* deFrame = dt.frame();
            ds1->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kBlue));
            ds2->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kRed));
            out.str("");
            out << "#Delta t, toy MC, bin " << i;
            TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
            cm->cd();

            deFrame->GetXaxis()->SetTitleSize(0.05);
            deFrame->GetXaxis()->SetTitleOffset(0.85);
            deFrame->GetXaxis()->SetLabelSize(0.05);
            deFrame->GetYaxis()->SetTitleOffset(1.6);
            deFrame->Draw();
        }
    }
    RooDataSet* ds1 = (RooDataSet*)d.reduce(SelectVars(dt),Cut("tag == 1"));
    RooDataSet* ds2 = (RooDataSet*)d.reduce(SelectVars(dt),Cut("tag == -1"));

    RooPlot* deFrame = dt.frame();
    ds1->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kBlue));
    ds2->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kRed));

    TCanvas* cm = new TCanvas("#Delta t, toy MC","#Delta t, toy MC",600,400);
    cm->cd();

    deFrame->GetXaxis()->SetTitleSize(0.05);
    deFrame->GetXaxis()->SetTitleOffset(0.85);
    deFrame->GetXaxis()->SetLabelSize(0.05);
    deFrame->GetYaxis()->SetTitleOffset(1.6);
    deFrame->Draw();

    return;
}
