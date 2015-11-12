#include "ToyMCParams.h"

using namespace RooFit;

void ToyMCFitter(const int nfits = 10000){
    ifstream ifile(fnamefit,ifstream::in);
    if(!ifile.is_open()){
        cout << "Can't open file " << fnamefit << endl;
        return;
    }
    string line;
    double nsig,pur,signif,cut;

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

    RooRealVar cosPull("cosPull","cosPull",-5.,5.);
    RooRealVar sinPull("sinPull","sinPull",-5.,5.);
    RooRealVar sigmaPull("sigmaPull","sigmaPull",-5.,5.);
    RooRealVar mistagPull("mistagPull","mistagPull",-5.,5.);

    RooRealVar tau("tau","tau",_tau,"ps"); tau.setConstant(kTRUE);
    RooRealVar dm("dm","dm",_dm,"ps^{-1}"); dm.setConstant(kTRUE);
    RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta,-3.,3.); if(constBeta) sin2beta.setConstant(kTRUE);
    RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta,-3.,3.); if(constBeta) cos2beta.setConstant(kTRUE);
    RooRealVar dt("dt","#Deltat",-5.,5.,"ps");
    RooRealVar avgMistag("avgMistag","avgMistag",mistag_rate,0.0,0.5); if(constMistag) avgMistag.setConstant(kTRUE);
    RooRealVar delMistag("delMistag","delMistag",0); delMistag.setConstant(kTRUE);
    RooRealVar mu("mu","mu",0); mu.setConstant(kTRUE);

    RooRealVar moment("moment","moment",0.); moment.setConstant(kTRUE);
    RooRealVar parity("parity","parity",-1.); parity.setConstant(kTRUE);

    RooRealVar* K[8];
    RooRealVar* Kb[7];
    RooFormulaVar* Kb8;
    RooArgList Kset;
    RooRealVar* C[8];
    RooRealVar* S[8];
    RooFormulaVar* a1[2][8];
    RooFormulaVar* b1[2][8];
    RooFormulaVar* a2[2][8];
    RooFormulaVar* b2[2][8];

    stringstream out;
    out.str("");
    for(int i=0; i<8; i++){
        out.str("");
        out << "K" << i+1;
        K[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_K[i],0.,1.); Kset.add(*K[i]); if(constK) K[i]->setConstant(kTRUE);
        if(i!=7){
          out.str("");
          out << "Kb" << i+1;
          Kb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_Kb[i],0.,1.); Kset.add(*Kb[i]); if(constK) Kb[i]->setConstant(kTRUE);
        }
        out.str("");
        out << "C" << i+1;
        C[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_C[i]); C[i]->setConstant(kTRUE);
        out.str("");
        out << "S" << i+1;
        S[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_S[i]); S[i]->setConstant(kTRUE);
    }
    Kb8 = new RooFormulaVar("K8b","K8b","1-@0-@1-@2-@3-@4-@5-@6-@7-@8-@9-@10-@11-@12-@13-@14",Kset);

    for(int i=0; i<7; i++){
        out.str("");
        out << "a10_" << i+1;
        a1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
        out.str("");
        out << "a11_" << i+1;
        a1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
        out.str("");
        out << "a20_" << i+1;
        a2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));
        out.str("");
        out << "a21_" << i+1;
        a2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1)/(@0+@1)",RooArgList(*K[i],*Kb[i]));

        out.str("");
        out << "b10_" << i+1;
        b1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1);",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta));
        out.str("");
        out << "b11_" << i+1;
        b1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta));
        out.str("");
        out << "b20_" << i+1;
        b2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta));
        out.str("");
        out << "b21_" << i+1;
        b2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta));
    }
    a1[0][7] = new RooFormulaVar("a10_8","a10_8","-(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a1[1][7] = new RooFormulaVar("a11_8","a11_8","(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a2[0][7] = new RooFormulaVar("a20_8","a20_8","(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a2[1][7] = new RooFormulaVar("a21_8","a21_8","-(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    b1[0][7] = new RooFormulaVar("b10_8","b10_8","2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b1[1][7] = new RooFormulaVar("b11_8","b11_8","2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b2[0][7] = new RooFormulaVar("b20_8","b20_8","-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b2[1][7] = new RooFormulaVar("b21_8","b21_8","-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));

    RooRealVar* dgamma = new RooRealVar("dgamma","dgamma",0.); dgamma->setConstant(kTRUE);
    RooRealVar* f0 = new RooRealVar("f0","f0",1.); f0->setConstant(kTRUE);
    RooRealVar* f1 = new RooRealVar("f1","f1",0.); f1->setConstant(kTRUE);

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

    RooSuperCategory bintag("bintag","bintag",RooArgSet(bin,tag));

    RooRealVar mean("mean","mean",0.,"ps"); mean.setConstant(kTRUE);
    RooRealVar sigma("sigma","sigma",_sigma_over_tau*_tau,0.,_tau,"ps"); if(constSigma) sigma.setConstant(kTRUE);
    RooGaussModel rf("rf","rf",dt,mean,sigma);
    RooGaussian rfpdf("rfpdf","rfpdf",dt,mean,sigma);

    cout << "Preparing PDFs..." << endl;
//    RooRealVar* fsigs1[8];
//    RooRealVar* fsigs1b[8];
//    RooRealVar* fsigs2[8];
//    RooRealVar* fsigs2b[8];
    RooRealVar fsig("fsig","fsig",_purity,0.,1.); if(constFSig) fsig.setConstant(kTRUE);
    RooBDecay* sigpdfs1[8];
    RooBDecay* sigpdfs1b[8];
    RooBDecay* sigpdfs2[8];
    RooBDecay* sigpdfs2b[8];
    RooAddPdf* PDFs1[8];
    RooAddPdf* PDFs1b[8];
    RooAddPdf* PDFs2[8];
    RooAddPdf* PDFs2b[8];
    RooAddPdf* pdfs1[8];
    RooAddPdf* pdfs1b[8];
    RooAddPdf* pdfs2[8];
    RooAddPdf* pdfs2b[8];
    RooSimultaneous pdf("pdf","pdf",bintag);
    for(int j=0; j<8; j++){
//        out.str("");
//        out << "fsig1" << j+1;
//        fsigs1[j] = new RooRealVar(out.str().c_str(),out.str().c_str(),_purity,0.,1.);  if(constFSig) fsigs1[j]->setConstant(kTRUE);
//        out.str("");
//        out << "fsig1b" << j+1;
//        fsigs1b[j] = new RooRealVar(out.str().c_str(),out.str().c_str(),_purity,0.,1.); if(constFSig) fsigs1b[j]->setConstant(kTRUE);
//        out.str("");
//        out << "fsig2" << j+1;
//        fsigs2[j] = new RooRealVar(out.str().c_str(),out.str().c_str(),_purity,0.,1.);  if(constFSig) fsigs2[j]->setConstant(kTRUE);
//        out.str("");
//        out << "fsig2b" << j+1;
//        fsigs2b[j] = new RooRealVar(out.str().c_str(),out.str().c_str(),_purity,0.,1.); if(constFSig) fsigs2b[j]->setConstant(kTRUE);

        out.str("");
        out << "sigpdf1" << j+1;
        sigpdfs1[j]  = new RooBDecay(out.str().c_str(),out.str().c_str(),dt,tau,*dgamma,*f0,*f1,*a1[0][j],*b1[0][j],dm,rf,RooBDecay::DoubleSided);
        out.str("");
        out << "sigpdf1b" << j+1;
        sigpdfs1b[j] = new RooBDecay(out.str().c_str(),out.str().c_str(),dt,tau,*dgamma,*f0,*f1,*a1[1][j],*b1[1][j],dm,rf,RooBDecay::DoubleSided);
        out.str("");
        out << "sigpdf2" << j+1;
        sigpdfs2[j]  = new RooBDecay(out.str().c_str(),out.str().c_str(),dt,tau,*dgamma,*f0,*f1,*a2[0][j],*b2[0][j],dm,rf,RooBDecay::DoubleSided);
        out.str("");
        out << "sigpdf2b" << j+1;
        sigpdfs2b[j] = new RooBDecay(out.str().c_str(),out.str().c_str(),dt,tau,*dgamma,*f0,*f1,*a2[1][j],*b2[1][j],dm,rf,RooBDecay::DoubleSided);

        out.str("");
        out << "PDF1" << j+1;
        PDFs1[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs1[j],rfpdf),RooArgList(fsig));
        out.str("");
        out << "PDF1b" << j+1;
        PDFs1b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs1b[j],rfpdf),RooArgList(fsig));
        out.str("");
        out << "PDF2" << j+1;
        PDFs2[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs2[j],rfpdf),RooArgList(fsig));
        out.str("");
        out << "PDF2b" << j+1;
        PDFs2b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*sigpdfs2b[j],rfpdf),RooArgList(fsig));

        //Adding mistaging
        out.str("");
        out << "pdf1" << j+1;
        pdfs1[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2[j],*PDFs1[j]),RooArgList(avgMistag));
        out.str("");
        out << "pdf1b" << j+1;
        pdfs1b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2b[j],*PDFs1b[j]),RooArgList(avgMistag));
        out.str("");
        out << "pdf2" << j+1;
        pdfs2[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1[j],*PDFs2[j]),RooArgList(avgMistag));
        out.str("");
        out << "pdf2b" << j+1;
        pdfs2b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1b[j],*PDFs2b[j]),RooArgList(avgMistag));

        out.str("");
        out << "{" << j+1 << ";B0}";
        pdf.addPdf(*pdfs1[j],out.str().c_str());
        out.str("");
        out << "{" << -(j+1) << ";B0}";
        pdf.addPdf(*pdfs1b[j],out.str().c_str());
        out.str("");
        out << "{" << j+1 << ";anti-B0}";
        pdf.addPdf(*pdfs2[j],out.str().c_str());
        out.str("");
        out << "{" << -(j+1) << ";anti-B0}";
        pdf.addPdf(*pdfs2b[j],out.str().c_str());
    }

    TFile* resfile;
    TTree* restree;

    RooArgSet ods_ds(cos2beta,sin2beta,sigma,avgMistag,cosPull,sinPull,sigmaPull,mistagPull);
    while(!ifile.eof()){
        getline(ifile,line);
        const int flag = sscanf(line.c_str(),"%lf %lf %lf %lf",&nsig,&pur,&signif,&cut);
        if(flag != 4) continue;
        cout << line << endl;
        int NSIG = (int)nsig;
        out.str("");
        out << "ToyData/toyMC_" << NSIG << "_" << pur << "_" << mistag_rate << "_" << _sigma_over_tau << "_" << PoissonFlag << ".root";
        TFile* file = TFile::Open(out.str().c_str());
        if(!file->IsOpen()){
          cout << "Can't open file " << out.str() << endl;
          continue;
        }
        TTree* tree = (TTree*)file->Get("ToyTree0");
        if(tree == NULL){
          cout << "The first tree in file " << out.str() << " is empty." << endl;
          continue;
        }

        out.str("");
        out << "ToyFitRes_" << NSIG << "_" << pur << "_" << mistag_rate << "_" << _sigma_over_tau << "_" << PoissonFlag << "fit_";
        if(constK)      out << "1";
        else            out << "0";
        if(constMistag) out << "1";
        else            out << "0";
        if(constFSig)   out << "1";
        else            out << "0";
        if(constSigma)  out << "1";
        else            out << "0";
        if(constBeta)   out << "1";
        else            out << "0";
        out << ".root";
        resfile = new TFile(out.str().c_str(),"RECREATE");
        restree = new TTree("ToyFitTree","ToyFitTree");
        double m_cos2beta,m_sin2beta,m_sigma,m_mistag,m_purity;
        double m_cos2betaErr,m_sin2betaErr,m_sigmaErr,m_mistagErr,m_purityErr;
        double m_cos2betaPull,m_sin2betaPull,m_sigmaPull,m_mistagPull,m_purityPull;
        restree->Branch("pur",&pur,"pur/D");

        restree->Branch("cos2beta",&m_cos2beta,"cos2beta/D");
        restree->Branch("sin2beta",&m_sin2beta,"sin2beta/D");
        restree->Branch("sigma",&m_sigma,"sigma/D");
        restree->Branch("mistag",&m_mistag,"mistag/D");
        restree->Branch("purity",&m_purity,"purity/D");

        restree->Branch("cos2betaErr",&m_cos2betaErr,"cos2betaErr/D");
        restree->Branch("sin2betaErr",&m_sin2betaErr,"sin2betaErr/D");
        restree->Branch("sigmaErr",&m_sigmaErr,"sigmaErrErr/D");
        restree->Branch("mistagErr",&m_mistagErr,"mistagErr/D");
        restree->Branch("purityErr",&m_purityErr,"purityErr/D");

        restree->Branch("cos2betaPull",&m_cos2betaPull,"cos2betaPull/D");
        restree->Branch("sin2betaPull",&m_sin2betaPull,"sin2betaPull/D");
        restree->Branch("sigmaPull",&m_sigmaPull,"sigmaPull/D");
        restree->Branch("mistagPull",&m_mistagPull,"mistagPull/D");
        restree->Branch("purityPull",&m_purityPull,"purityPull/D");

        int jj = 0;
        RooDataSet ods("ods","ods",ods_ds);
        while(tree != NULL && jj<nfits){
            RooDataSet d("data","data",tree,RooArgSet(dt,bin,tag));
            d.Print();

//            tau.setValue(_tau);
            sin2beta.setVal(_sin2beta);
            cos2beta.setVal(_cos2beta);
            sigma.setVal(_tau*_sigma_over_tau);
            fsig.setVal(pur);
            avgMistag.setVal(mistag_rate);
            for(int i=0; i<8; i++){
                K[i]->setVal(_K[i]);
                if(i != 7)Kb[i]->setVal(_Kb[i]);
            }

            pdf.fitTo(d,PrintLevel(-1),NumCPU(1));

            m_cos2beta = cos2beta.getVal();  m_cos2betaErr = cos2beta.getError();  m_cos2betaPull = m_cos2betaErr != 0 ? (m_cos2beta-_cos2beta)/m_cos2betaErr     : 0;
            m_sin2beta = sin2beta.getVal();  m_sin2betaErr = sin2beta.getError();  m_sin2betaPull = m_sin2betaErr != 0 ? (m_sin2beta-_sin2beta)/m_sin2betaErr     : 0;
            m_sigma    = sigma.getVal();     m_sigmaErr    = sigma.getError();     m_sigmaPull    = m_sigmaErr != 0 ? (m_sigma-_tau*_sigma_over_tau)/(m_sigmaErr) : 0;
            m_mistag   = avgMistag.getVal(); m_mistagErr   = avgMistag.getError(); m_mistagPull   = m_mistagErr != 0 ? (m_mistag-mistag_rate)/(m_mistagErr)       : 0;
            m_purity   = fsig.getVal();      m_purityErr   = fsig.getError();      m_purityPull   = m_purityErr != 0 ? (m_purity-_purity)/(m_purityErr)            : 0;
            restree->Fill();

            cosPull    = m_cos2betaPull;
            sinPull    = m_sin2betaPull;
            sigmaPull  = m_sigmaPull;
            mistagPull = m_mistagPull;
            ods.add(ods_ds);

            jj++;
            if(!(jj%100) && jj) cout << jj << " fits performed" << endl;
//            delete tree;
            out.str("");
            out << "ToyTree" << jj;
            tree = (TTree*)file->Get(out.str().c_str());
//            tree->Print();
//            cout << "Entries: " << tree->GetEntries() << endl;
        }
        file->Close();

        RooPlot* cosFrame = cos2beta.frame();
        ods.plotOn(cosFrame,MarkerSize(1));
        ods.statOn(cosFrame);
        out.str("");
        out << "cos2#beta, pur = " << pur;
        TCanvas* ccos = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
        ccos->cd();
        cosFrame->GetXaxis()->SetTitleSize(0.05);
        cosFrame->GetXaxis()->SetTitleOffset(0.85);
        cosFrame->GetXaxis()->SetLabelSize(0.05);
        cosFrame->GetYaxis()->SetTitleOffset(1.6);
        cosFrame->Draw();

        RooPlot* sinFrame = sin2beta.frame();
        ods.plotOn(sinFrame,MarkerSize(1));
        ods.statOn(sinFrame);
        out << "sin2#beta, pur = " << pur;
        TCanvas* csin = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
        csin->cd();
        sinFrame->GetXaxis()->SetTitleSize(0.05);
        sinFrame->GetXaxis()->SetTitleOffset(0.85);
        sinFrame->GetXaxis()->SetLabelSize(0.05);
        sinFrame->GetYaxis()->SetTitleOffset(1.6);
        sinFrame->Draw();

        resfile->cd();
        restree->Print();
        restree->Write();
        resfile->Close();
    }
    ifile.close();

    return;
}

