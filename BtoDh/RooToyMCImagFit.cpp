#include "ToyMCParams.h"

using namespace RooFit;

void RooToyMCImagFit(){
    stringstream out;
   out.str("");
   out << "toyMC_" << _sigma_over_tau << "_" << _purity << "_" << mistag_rate << ".root";
   TFile* file = TFile::Open(out.str().c_str());
   TTree* tree = (TTree*)file->Get("ToyTree");

//    TFile* file = TFile::Open("ToyData/toyMC_10000_0.9_0.3_0.8_0.root");
    TTree* tree = (TTree*)file->Get("ToyTree");

    RooRealVar tau("tau","tau",_tau,"ps"); tau.setConstant(kTRUE);
    RooRealVar dm("dm","dm",_dm,"ps^{-1}"); dm.setConstant(kTRUE);

    RooRealVar phi("phi","phi",_phi,0.,2*TMath::Pi()); if(constBeta) phi.setConstant(kTRUE);
    RooRealVar psi("psi","psi",1.,0.1,10.); if(constPsi) psi.setConstant(kTRUE);

    RooFormulaVar sin2beta("sin2beta","sin2beta","TMath::Sin(@0)",RooArgSet(phi));
    RooFormulaVar cos2beta("cos2beta","cos2beta","TMath::Cos(@0)",RooArgSet(phi));

//    RooRealVar sin2beta("sin2beta","sin2beta",_sin2beta,-3.,3.); if(constBeta) sin2beta.setConstant(kTRUE);
//    RooRealVar cos2beta("cos2beta","cos2beta",_cos2beta,-3.,3.); if(constBeta) cos2beta.setConstant(kTRUE);
    RooRealVar dt("dt","#Deltat",-5.,5.,"ps");
    RooRealVar avgMisgat("avgMisgat","avgMisgat",mistag_rate,0.0,0.5); if(constMistag) avgMisgat.setConstant(kTRUE);
    RooRealVar delMisgat("delMisgat","delMisgat",0); delMisgat.setConstant(kTRUE);
    RooRealVar mu("mu","mu",0); mu.setConstant(kTRUE);

    tau.Print();
    dm.Print();
    sin2beta.Print();
    cos2beta.Print();
    avgMisgat.Print();

    RooRealVar moment("moment","moment",0.); moment.setConstant(kTRUE);
    RooRealVar parity("parity","parity",-1.); parity.setConstant(kTRUE);

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
        K[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_K[i],0.,1.); Kset.add(*K[i]); if(constK) K[i]->setConstant(kTRUE);
        K[i]->Print();
        if(i!=7){
          out.str("");
          out << "Kb" << i+1;
          Kb[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_Kb[i],0.,1.); Kset.add(*Kb[i]); if(constK)((RooRealVar*)Kb[i])->setConstant(kTRUE);
          Kb[i]->Print();
        }
        out.str("");
        out << "C" << i+1;
        C[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_C[i]); C[i]->setConstant(kTRUE);
        C[i]->Print();
        out.str("");
        out << "S" << i+1;
        S[i] = new RooRealVar(out.str().c_str(),out.str().c_str(),_S[i]); S[i]->setConstant(kTRUE);
        S[i]->Print();
    }
    Kb[7] = new RooFormulaVar("K8b","K8b","1-@0-@1-@2-@3-@4-@5-@6-@7-@8-@9-@10-@11-@12-@13-@14",Kset);

    for(int i=0; i<8; i++){
        out.str("");
        out << "a10_" << i+1;
        a1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1*@2)/(@0+@1*@2)",RooArgList(*K[i],*Kb[i],psi));
        a1[0][i]->Print();
        out.str("");
        out << "a11_" << i+1;
        a1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0*@2-@1)/(@0*@2+@1)",RooArgList(*K[i],*Kb[i],psi));
        a1[1][i]->Print();
        out.str("");
        out << "a20_" << i+1;
        a2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"(@0/@2-@1)/(@0/@2+@1)",RooArgList(*K[i],*Kb[i],psi));
        a2[0][i]->Print();
        out.str("");
        out << "a21_" << i+1;
        a2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-(@0-@1/@2)/(@0+@1/@2)",RooArgList(*K[i],*Kb[i],psi));
        a2[1][i]->Print();

        out.str("");
        out << "b10_" << i+1;
        b1[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1*@6)/(@0+@1*@6)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,psi));
        b1[0][i]->Print();
        out.str("");
        out << "b11_" << i+1;
        b1[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1*@6)/(@0*@6+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,psi));
        b1[1][i]->Print();
        out.str("");
        out << "b20_" << i+1;
        b2[0][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1/@6)/(@0/@6+@1)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,psi));
        b2[0][i]->Print();
        out.str("");
        out << "b21_" << i+1;
        b2[1][i] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1/@6)/(@0+@1/@6)",RooArgList(*K[i],*Kb[i],*C[i],*S[i],sin2beta,cos2beta,psi));
        b2[1][i]->Print();

        cout << "bin = " << i+1 << endl;
        cout << "  a11 = " << a1[0][i]->getVal() << " b11 = " << b1[0][i]->getVal() << endl;
        cout << "  a12 = " << a1[1][i]->getVal() << " b12 = " << b1[1][i]->getVal() << endl;
        cout << "  a21 = " << a2[0][i]->getVal() << " b21 = " << b2[0][i]->getVal() << endl;
        cout << "  a22 = " << a2[1][i]->getVal() << " b22 = " << b2[1][i]->getVal() << endl;
    }
/*    a1[0][7] = new RooFormulaVar("a10_8","a10_8","-(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a1[1][7] = new RooFormulaVar("a11_8","a11_8","(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a2[0][7] = new RooFormulaVar("a20_8","a20_8","(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    a2[1][7] = new RooFormulaVar("a21_8","a21_8","-(@0-@1)/(@0+@1)",RooArgList(*K[7],*Kb8));
    b1[0][7] = new RooFormulaVar("b10_8","b10_8","2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b1[1][7] = new RooFormulaVar("b11_8","b11_8","2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b2[0][7] = new RooFormulaVar("b20_8","b20_8","-2.*(@2*@4+@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));
    b2[1][7] = new RooFormulaVar("b21_8","b21_8","-2.*(@2*@4-@3*@5)*TMath::Sqrt(@0*@1)/(@0+@1)",RooArgList(*K[7],*Kb8,*C[7],*S[7],sin2beta,cos2beta));

    a1[0][7]->Print();
    a1[1][7]->Print();
    a2[0][7]->Print();
    a2[1][7]->Print();
    b1[0][7]->Print();
    b1[1][7]->Print();
    b2[0][7]->Print();
    b2[1][7]->Print();
*/
    RooRealVar* dgamma = new RooRealVar("dgamma","dgamma",0.); dgamma->setConstant(kTRUE);
    RooRealVar* f0 = new RooRealVar("f0","f0",1.);             f0->setConstant(kTRUE);
    RooRealVar* f1 = new RooRealVar("f1","f1",0.);             f1->setConstant(kTRUE);

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

    RooDataSet d("data","data",tree,RooArgSet(dt,bin,tag));
    cout << "DataSet is ready." << endl;
    d.Print();

    RooRealVar mean("mean","mean",0.,"ps"); mean.setConstant(kTRUE);
    RooRealVar sigma("sigma","sigma",_sigma_over_tau*_tau,0.,1.5*_tau,"ps"); if(constSigma) sigma.setConstant(kTRUE);
    RooGaussModel rf("rf","rf",dt,mean,sigma);
//    RooTruthModel rf("rf","rf",dt);
    RooGaussian rfpdf("rfpdf","rfpdf",dt,mean,sigma);

    cout << "Preparing PDFs..." << endl;
//    RooRealVar* fsigs1[8];
//    RooRealVar* fsigs1b[8];
//    RooRealVar* fsigs2[8];
//    RooRealVar* fsigs2b[8];
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
    RooRealVar fsig("fsig","fsig",_purity,0.,1.); if(constFSig) fsig.setConstant(kTRUE);
    fsig.Print();
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
        pdfs1[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2[j],*PDFs1[j]),RooArgList(avgMisgat));
        out.str("");
        out << "pdf1b" << j+1;
        pdfs1b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs2b[j],*PDFs1b[j]),RooArgList(avgMisgat));
        out.str("");
        out << "pdf2" << j+1;
        pdfs2[j]  = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1[j],*PDFs2[j]),RooArgList(avgMisgat));
        out.str("");
        out << "pdf2b" << j+1;
        pdfs2b[j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*PDFs1b[j],*PDFs2b[j]),RooArgList(avgMisgat));

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

    cout << "Fitting..." << endl;
    pdf.fitTo(d,Verbose(),Timer());
    dt.setBins(50);

    cout << "Drawing plots." << endl;
    for(int i=0; i<8; i++){
        RooPlot* dtFrame = dt.frame();
        out.str("");
        out << "tag == 1 && bin == " << i+1;
        RooDataSet* ds = d.reduce(out.str().c_str());
        ds->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),Cut(out.str().c_str()),MarkerColor(kBlue));

        out.str("");
        out << "{" << j+1 << ";B0}";
        bintag = out.str().c_str();
        pdf.plotOn(dtFrame,ProjWData(*ds),Slice(bintag),LineColor(kBlue));
        double chi2 = dtFrame->chiSquare();

        out.str("");
        out << "tag == -1 && bin == " << i+1;
        RooDataSet* dsb = d.reduce(out.str().c_str());
        dsb->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),Cut(out.str().c_str()),MarkerColor(kRed));

        out.str("");
        out << "{" << j+1 << ";anti-B0}";
        bintag = out.str().c_str();
        pdf.plotOn(dtFrame,ProjWData(*dsb),Slice(bintag),LineColor(kRed));
        double chi2b = dtFrame->chiSquare();

        out.str("");
        out << "#Delta t, toy MC, bin == " << i+1;
        TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
        cm->cd();
        dtFrame->GetXaxis()->SetTitleSize(0.05);
        dtFrame->GetXaxis()->SetTitleOffset(0.85);
        dtFrame->GetXaxis()->SetLabelSize(0.05);
        dtFrame->GetYaxis()->SetTitleOffset(1.6);

        TPaveText *ptB = new TPaveText(0.7,0.65,0.98,0.99,"brNDC");
        ptB->SetFillColor(0);
        ptB->SetTextAlign(12);
        out.str("");
        out << "bin = " << (i+1);
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#chi^{2}(B^{0}) = " << chi2;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#chi^{2}(#barB^{0}) = " << chi2b;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#sigma/#tau = " << _sigma_over_tau;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "purity = " << _purity;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "mistag = " << mistag_rate;
        ptB->AddText(out.str().c_str());

        dtFrame->Draw();
        ptB->Draw();

//        out.str("");
//        out << "../Reports/cpfitfigs/dt_bin" << (i+1) << ".root";
//        cm->Print(out.str().c_str());
//        out.str("");
//        out << "../Reports/cpfitfigs/dt_bin" << (i+1) << ".png";
//        cm->Print(out.str().c_str());
    }

    for(int i=0; i<8; i++){
        RooPlot* dtFrame = dt.frame();
        out.str("");
        out << "tag == 1 && bin == " << -(i+1);
        RooDataSet* ds2 = (RooDataSet*)d.reduce(out.str().c_str());
        ds2->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),Cut(out.str().c_str()),MarkerColor(kBlue));

        out.str("");
        out << "{" << -(j+1) << ";B0}";
        bintag = out.str().c_str();
        pdf.plotOn(dtFrame,ProjWData(*ds2),Slice(bintag),LineColor(kBlue));
        double chi2 = dtFrame->chiSquare();

        out.str("");
        out << "tag == -1 && bin == " << -(i+1);
        RooDataSet* ds2b = (RooDataSet*)d.reduce(out.str().c_str());
        ds2b->plotOn(dtFrame,DataError(RooAbsData::SumW2),MarkerSize(1),Cut(out.str().c_str()),MarkerColor(kRed));

        out.str("");
        out << "{" << -(j+1) << ";anti-B0}";
        bin = out.str().c_str();
        pdf.plotOn(dtFrame,ProjWData(*ds2b),Slice(bintag),LineColor(kRed));
        double chi2b = dtFrame->chiSquare();

        out.str("");
        out << "#Delta t, toy MC, bin == " << -(i+1);
        TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,400);
        cm->cd();
        dtFrame->GetXaxis()->SetTitleSize(0.05);
        dtFrame->GetXaxis()->SetTitleOffset(0.85);
        dtFrame->GetXaxis()->SetLabelSize(0.05);
        dtFrame->GetYaxis()->SetTitleOffset(1.6);

        TPaveText *ptB = new TPaveText(0.7,0.65,0.98,0.99,"brNDC");
        ptB->SetFillColor(0);
        ptB->SetTextAlign(12);
        out.str("");
        out << "bin = " << -(i+1);
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#chi^{2}(B^{0}) = " << chi2;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#chi^{2}(#barB^{0}) = " << chi2b;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "#sigma/#tau = " << _sigma_over_tau;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "purity = " << _purity;
        ptB->AddText(out.str().c_str());
        out.str("");
        out << "mistag = " << mistag_rate;
        ptB->AddText(out.str().c_str());

        dtFrame->Draw();
        ptB->Draw();

//        out.str("");
//        out << "../Reports/cpfitfigs/dt_bin" << -(i+1) << ".root";
//        cm->Print(out.str().c_str());
//        out.str("");
//        out << "../Reports/cpfitfigs/dt_bin" << -(i+1) << ".png";
//        cm->Print(out.str().c_str());
    }

    return;
}

