#include "ToyMCParams.h"
using namespace RooFit;

void RooToyFitAnalysis(const char* input){
    TFile* file = TFile::Open(input);
    TTree* tree = (TTree*)file->Get("ToyFitTree");
    tree->Print();
    double SINErr,COSErr,PURErr,MISErr,SIGErr;
    tree->SetBranchAddress("sin2betaErr",&SINErr);
    tree->SetBranchAddress("cos2betaErr",&COSErr);
    tree->SetBranchAddress("purityErr",&PURErr);
    tree->SetBranchAddress("mistagErr",&MISErr);
    tree->SetBranchAddress("sigmaErr",&SIGErr);
    tree->GetEvent(0);
    const bool run_sin = SINErr != 0 ? true : false;
    const bool run_cos = COSErr != 0 ? true : false;
    const bool run_pur = PURErr != 0 ? true : false;
    const bool run_mis = MISErr != 0 ? true : false;
    const bool run_sig = SIGErr != 0 ? true : false;

    const double dsin = 1.5;

    RooArgSet argset;
    RooRealVar sin2beta("sin2beta","sin(2#beta)",_sin2beta-dsin,_sin2beta+dsin); argset.add(sin2beta);
    RooRealVar cos2beta("cos2beta","cos(2#beta)",_cos2beta-dsin,_cos2beta+dsin); argset.add(cos2beta);
    RooRealVar purity("purity","Purity",0.,1.);                                  argset.add(purity);
    RooRealVar mistag("mistag","Miss tag rate",0.,0.5);                          argset.add(mistag);
    RooRealVar sigma("sigma","Sigma",0.,_sigma_over_tau*_tau+dsin);              argset.add(sigma);

    RooRealVar sin2betaPull("sin2betaPull","sin(2#beta) Pull",-5.,5.); argset.add(sin2betaPull);
    RooRealVar cos2betaPull("cos2betaPull","cos(2#beta) Pull",-5.,5.); argset.add(cos2betaPull);
    RooRealVar purityPull("purityPull","Purity Pull",-5.,5.);          argset.add(purityPull);
    RooRealVar mistagPull("mistagPull","Miss tag rate Pull",-5.,5.);   argset.add(mistagPull);
    RooRealVar sigmaPull("sigmaPull","Sigma Pull",-5.,5.);             argset.add(sigmaPull);

    RooDataSet ds("ds","ds",tree,argset);
    ds.Print();

    /////////
    // sin //
    /////////
    if(run_sin){
    RooRealVar msin("msin","msin",_sin2beta,_sin2beta-0.3,_sin2beta+0.3);
    RooRealVar ssin("ssin","ssin",0.1,0.01,0.5);
    RooGaussian gsin("gsin","gsin",sin2beta,msin,ssin);
    gsin.fitTo(ds);

    RooPlot* sinFrame = sin2beta.frame();
    ds.plotOn(sinFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    gsin.plotOn(sinFrame);
    ds.statOn(sinFrame,Layout(0.6,0.91,0.91));

    TCanvas* csin = new TCanvas("sin(2#beta)","sin(2#beta)",600,400);
    csin->cd();

    sinFrame->GetXaxis()->SetTitleSize(0.05);
    sinFrame->GetXaxis()->SetTitleOffset(0.85);
    sinFrame->GetXaxis()->SetLabelSize(0.05);
    sinFrame->GetYaxis()->SetTitleOffset(1.6);

    TLine *sinline = new TLine(_sin2beta,0,_sin2beta,40);
    sinline->SetLineColor(kRed);
    sinline->SetLineStyle(1);
    sinline->SetLineWidth(2);
    sinFrame->Draw();
    sinline->Draw();
    csin->Print("../Reports/cpfitfigs/sin.png");
    csin->Print("../Reports/cpfitfigs/sin.root");

    RooRealVar msinp("msinp","msinp",0,-5.,5.);
    RooRealVar ssinp("ssinp","ssinp",1.,0.3,3.);
    RooGaussian gsinp("gsinp","gsinp",sin2betaPull,msinp,ssinp);
    gsinp.fitTo(ds);

    RooPlot* sinFramePull = sin2betaPull.frame();
    ds.plotOn(sinFramePull,DataError(RooAbsData::SumW2),MarkerSize(1));
    gsinp.plotOn(sinFramePull);
    ds.statOn(sinFramePull,Layout(0.6,0.91,0.91));

    TCanvas* csinp = new TCanvas("sin(2#beta) pull","sin(2#beta) pull",600,400);
    csinp->cd();

    sinFramePull->GetXaxis()->SetTitleSize(0.05);
    sinFramePull->GetXaxis()->SetTitleOffset(0.85);
    sinFramePull->GetXaxis()->SetLabelSize(0.05);
    sinFramePull->GetYaxis()->SetTitleOffset(1.6);
    sinFramePull->Draw();
    csinp->Print("../Reports/cpfitfigs/sinpull.png");
    csinp->Print("../Reports/cpfitfigs/sinpull.root");
    }


    /////////
    // cos //
    /////////
    if(run_cos){
    RooRealVar mcos("mcos","mcos",_cos2beta,_cos2beta-0.3,_cos2beta+0.3);
    RooRealVar scos("scos","scos",0.1,0.01,0.5);
    RooGaussian gcos("gcos","gcos",cos2beta,mcos,scos);
    gcos.fitTo(ds);

    RooPlot* cosFrame = cos2beta.frame();
    ds.plotOn(cosFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    gcos.plotOn(cosFrame);
    ds.statOn(cosFrame,Layout(0.6,0.91,0.91));

    TCanvas* ccos = new TCanvas("cos(2#beta)","cos(2#beta)",600,400);
    ccos->cd();

    cosFrame->GetXaxis()->SetTitleSize(0.05);
    cosFrame->GetXaxis()->SetTitleOffset(0.85);
    cosFrame->GetXaxis()->SetLabelSize(0.05);
    cosFrame->GetYaxis()->SetTitleOffset(1.6);

    TLine *cosline = new TLine(_cos2beta,0,_cos2beta,40);
    cosline->SetLineColor(kRed);
    cosline->SetLineStyle(1);
    cosline->SetLineWidth(2);
    cosFrame->Draw();
    cosline->Draw();
    ccos->Print("../Reports/cpfitfigs/cos.png");
    ccos->Print("../Reports/cpfitfigs/cos.root");

    RooRealVar mcosp("mcosp","mcosp",0,-5.,5.);
    RooRealVar scosp("scosp","scosp",1.,0.3,3.);
    RooGaussian gcosp("gcosp","gcosp",cos2betaPull,mcosp,scosp);
    gcosp.fitTo(ds);

    RooPlot* cosFramePull = cos2betaPull.frame();
    ds.plotOn(cosFramePull,DataError(RooAbsData::SumW2),MarkerSize(1));
    gcosp.plotOn(cosFramePull);
    ds.statOn(cosFramePull,Layout(0.6,0.91,0.91));

    TCanvas* ccosp = new TCanvas("cos(2#beta) pull","cos(2#beta) pull",600,400);
    ccosp->cd();

    cosFramePull->GetXaxis()->SetTitleSize(0.05);
    cosFramePull->GetXaxis()->SetTitleOffset(0.85);
    cosFramePull->GetXaxis()->SetLabelSize(0.05);
    cosFramePull->GetYaxis()->SetTitleOffset(1.6);
    cosFramePull->Draw();
    ccosp->Print("../Reports/cpfitfigs/cospull.png");
    ccosp->Print("../Reports/cpfitfigs/cospull.root");
    }

    ////////////
    // purity //
    ////////////
    if(run_pur){
    RooRealVar mpur("mpur","mpur",_purity,_purity-0.3,_purity+0.3);
    RooRealVar spur("spur","spur",0.1,0.01,0.5);
    RooGaussian gpur("gpur","gpur",purity,mpur,spur);
    gpur.fitTo(ds);

    RooPlot* purFrame = purity.frame();
    ds.plotOn(purFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    gpur.plotOn(purFrame);
    ds.statOn(purFrame,Layout(0.6,0.91,0.91));

    TCanvas* cpur = new TCanvas("Purity","Purity",600,400);
    cpur->cd();

    purFrame->GetXaxis()->SetTitleSize(0.05);
    purFrame->GetXaxis()->SetTitleOffset(0.85);
    purFrame->GetXaxis()->SetLabelSize(0.05);
    purFrame->GetYaxis()->SetTitleOffset(1.6);

    TLine *purline = new TLine(_purity,0,_purity,40);
    purline->SetLineColor(kRed);
    purline->SetLineStyle(1);
    purline->SetLineWidth(2);
    purFrame->Draw();
    purline->Draw();
    cpur->Print("../Reports/cpfitfigs/pur.png");
    cpur->Print("../Reports/cpfitfigs/pur.root");

    RooRealVar mpurp("mpurp","mpurp",0,-5.,5.);
    RooRealVar spurp("spurp","spurp",1.,0.3,3.);
    RooGaussian gpurp("gpurp","gpurp",purityPull,mpurp,spurp);
    gpurp.fitTo(ds);

    RooPlot* purFramePull = purityPull.frame();
    ds.plotOn(purFramePull,DataError(RooAbsData::SumW2),MarkerSize(1));
    gpurp.plotOn(purFramePull);
    ds.statOn(purFramePull,Layout(0.6,0.91,0.91));

    TCanvas* cpurp = new TCanvas("purity pull","purity pull",600,400);
    cpurp->cd();

    purFramePull->GetXaxis()->SetTitleSize(0.05);
    purFramePull->GetXaxis()->SetTitleOffset(0.85);
    purFramePull->GetXaxis()->SetLabelSize(0.05);
    purFramePull->GetYaxis()->SetTitleOffset(1.6);
    purFramePull->Draw();
    cpurp->Print("../Reports/cpfitfigs/purpull.png");
    cpurp->Print("../Reports/cpfitfigs/purpull.root");
    }

    ////////////
    // mistag //
    ////////////
    if(run_mis){
    RooRealVar mmis("mmis","mmis",mistag_rate,0.,0.5);
    RooRealVar smis("smis","smis",0.1,0.01,0.5);
    RooGaussian gmis("gmis","gmis",mistag,mmis,smis);
    gmis.fitTo(ds);

    RooPlot* misFrame = mistag.frame();
    ds.plotOn(misFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    gmis.plotOn(misFrame);
    ds.statOn(misFrame,Layout(0.6,0.91,0.91));

    TCanvas* cmis = new TCanvas("Mistag","Mistag",600,400);
    cmis->cd();

    misFrame->GetXaxis()->SetTitleSize(0.05);
    misFrame->GetXaxis()->SetTitleOffset(0.85);
    misFrame->GetXaxis()->SetLabelSize(0.05);
    misFrame->GetYaxis()->SetTitleOffset(1.6);

    TLine *misline = new TLine(mistag_rate,0,mistag_rate,40);
    misline->SetLineColor(kRed);
    misline->SetLineStyle(1);
    misline->SetLineWidth(2);
    misFrame->Draw();
    misline->Draw();
    cmis->Print("../Reports/cpfitfigs/mis.png");
    cmis->Print("../Reports/cpfitfigs/mis.root");

    RooRealVar mmisp("mmisp","mmisp",0,-5.,5.);
    RooRealVar smisp("smisp","smisp",1.,0.3,3.);
    RooGaussian gmisp("gmisp","gmisp",mistagPull,mmisp,smisp);
    gmisp.fitTo(ds);

    RooPlot* misFramePull = mistagPull.frame();
    ds.plotOn(misFramePull,DataError(RooAbsData::SumW2),MarkerSize(1));
    gmisp.plotOn(misFramePull);
    ds.statOn(misFramePull,Layout(0.6,0.91,0.91));

    TCanvas* cmisp = new TCanvas("mistag pull","mistag pull",600,400);
    cmisp->cd();

    misFramePull->GetXaxis()->SetTitleSize(0.05);
    misFramePull->GetXaxis()->SetTitleOffset(0.85);
    misFramePull->GetXaxis()->SetLabelSize(0.05);
    misFramePull->GetYaxis()->SetTitleOffset(1.6);
    misFramePull->Draw();
    cmisp->Print("../Reports/cpfitfigs/mispull.png");
    cmisp->Print("../Reports/cpfitfigs/mispull.root");
    }

    ///////////
    // sigma //
    ///////////
    if(run_sig){
    RooRealVar msig("msig","msig",_sigma_over_tau*_tau,0.,_tau);
    RooRealVar ssig("ssig","ssig",0.1,0.01,0.5);
    RooGaussian gsig("gsig","gsig",sigma,msig,ssig);
    gsig.fitTo(ds);

    RooPlot* sigFrame = sigma.frame();
    ds.plotOn(sigFrame,DataError(RooAbsData::SumW2),MarkerSize(1));
    gsig.plotOn(sigFrame);
    ds.statOn(sigFrame,Layout(0.6,0.91,0.91));

    TCanvas* csig = new TCanvas("Sigma","Sigma",600,400);
    csig->cd();

    sigFrame->GetXaxis()->SetTitleSize(0.05);
    sigFrame->GetXaxis()->SetTitleOffset(0.85);
    sigFrame->GetXaxis()->SetLabelSize(0.05);
    sigFrame->GetYaxis()->SetTitleOffset(1.6);

    TLine *sigline = new TLine(_sigma_over_tau*_tau,0,_sigma_over_tau*_tau,70);
    sigline->SetLineColor(kRed);
    sigline->SetLineStyle(1);
    sigline->SetLineWidth(2);
    sigFrame->Draw();
    sigline->Draw();
    csig->Print("../Reports/cpfitfigs/sig.png");
    csig->Print("../Reports/cpfitfigs/sig.root");

    RooRealVar msigp("msigp","msigp",0,-5.,5.);
    RooRealVar ssigp("ssigp","ssigp",1.,0.3,3.);
    RooGaussian gsigp("gsigp","gsigp",sigmaPull,msigp,ssigp);
    gsigp.fitTo(ds);

    RooPlot* sigFramePull = sigmaPull.frame();
    ds.plotOn(sigFramePull,DataError(RooAbsData::SumW2),MarkerSize(1));
    gsigp.plotOn(sigFramePull);
    ds.statOn(sigFramePull,Layout(0.6,0.91,0.91));

    TCanvas* csigp = new TCanvas("sigma pull","sigma pull",600,400);
    csigp->cd();

    sigFramePull->GetXaxis()->SetTitleSize(0.05);
    sigFramePull->GetXaxis()->SetTitleOffset(0.85);
    sigFramePull->GetXaxis()->SetLabelSize(0.05);
    sigFramePull->GetYaxis()->SetTitleOffset(1.6);
    sigFramePull->Draw();
    csigp->Print("../Reports/cpfitfigs/sigpull.png");
    csigp->Print("../Reports/cpfitfigs/sigpull.root");
    }

    return;
}
