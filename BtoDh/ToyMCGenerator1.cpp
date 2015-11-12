#include "ToyMCParams.h"

TTree* GenerateTree(const string& treename,const int nsig,const double& pur){
    TRandom3 rndm(0);
    int N;
    stringstream out;
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

    RooRealVar mean("mean","mean",0.,"ps"); mean.setConstant(kTRUE);
    RooRealVar sigma("sigma","sigma",_sigma_over_tau*_tau,"ps"); sigma.setConstant(kTRUE);
    RooGaussModel rf("rf","rf",dt,mean,sigma);

    RooRealVar dgamma("dgamma","dgamma",0.); dgamma.setConstant(kTRUE);
    RooRealVar f0("f0","f0",1.); f0.setConstant(kTRUE);
    RooRealVar f1("f1","f1",0.); f1.setConstant(kTRUE);

    RooBDecay pdf1("sigpdf1","sigpdf1",dt,tau,dgamma,f0,f1,a1,b1,dm,rf,RooBDecay::DoubleSided);
    RooBDecay pdf2("sigpdf2","sigpdf2",dt,tau,dgamma,f0,f1,a2,b2,dm,rf,RooBDecay::DoubleSided);
    RooAddPdf pdf("pdf","pdf",RooArgList(pdf2,pdf1),RooArgList(avgMisgat));

    double _dt;
    int _bin,_tag;
    int flag;
    TTree* tree = new TTree(treename.c_str(),treename.c_str());
    tree->Branch("dt",&_dt,"dt/D");
    tree->Branch("bin",&_bin,"bin/I");
    tree->Branch("tag",&_tag,"tag/I");
    RooDataSet* data;

    for(int i=0; i<8; i++){
      _bin = i+1;
      // B0
      if(PoissonFlag) N = rndm.Poisson(nsig*_K[i]*0.5);
      else            N = (int)(nsig*_K[i]*0.5);
      a1 = -(_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
      b1 = 2.*(_C[i]*_sin2beta+_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
      _tag = 1;
      data = pdf.generate(RooArgSet(dt),N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      // anti-B0
      a1 = (_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
      b1 = -2.*(_C[i]*_sin2beta+_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
      _tag = -1;
      if(PoissonFlag) N = rndm.Poisson(nsig*_K[i]*0.5);
      data = pdf.generate(RooArgSet(dt),N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      //Background
      if(PoissonFlag) N = rndm.Poisson((nsig/pur-nsig)*_K[i]*0.5);
      else            N = (int)((nsig/pur-nsig)*_K[i]*0.5);
      // B0
      _tag = 1;
      data = rf.generate(dt,N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      // anti-B0
      _tag = -1;
      if(PoissonFlag) N = rndm.Poisson((nsig/pur-nsig)*_K[i]*0.5);
      data = rf.generate(dt,N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      //////////////
      // anti-bin //
      //////////////
      _bin = -(i+1);
      if(PoissonFlag) N = rndm.Poisson(nsig*_Kb[i]*0.5);
      else            N = (int)(nsig*_Kb[i]*0.5);
      // B0
      a1 = (_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
      b1 = 2.*(_C[i]*_sin2beta-_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
      _tag = 1;
      data = pdf.generate(RooArgSet(dt),N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      // anti-B0
      a1 = -(_K[i]-_Kb[i])/(_K[i]+_Kb[i]);
      b1 = -2.*(_C[i]*_sin2beta-_S[i]*_cos2beta)*TMath::Sqrt(_K[i]*_Kb[i])/(_K[i]+_Kb[i]);
      _tag = -1;
      if(PoissonFlag) N = rndm.Poisson(nsig*_Kb[i]*0.5);
      data = pdf.generate(RooArgSet(dt),N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      //Background
      if(PoissonFlag) N = rndm.Poisson((nsig/pur-nsig)*_Kb[i]*0.5);
      else            N = (int)((nsig/pur-nsig)*_Kb[i]*0.5);
      // B0
      _tag = 1;
      data = rf.generate(dt,N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }

      // anti-B0
      _tag = -1;
      if(PoissonFlag) N = rndm.Poisson((nsig/pur-nsig)*_Kb[i]*0.5);
      data = rf.generate(dt,N);
      for(int j=0; j<N; j++){
	_dt = ((RooRealVar*)data->get(j)->find(dt.GetName()))->getVal();
	tree->Fill();
      }
    }

    return tree;
}

void ToyMCGenerator1(const int nsamples = 1000){
    ifstream ifile(fname,ifstream::in);
    if(!ifile.is_open()){
        cout << "Can't open file " << fname << endl;
        return;
    }
    string line;
    stringstream out;
    TTree* tree;
    double nsig,pur,signif,cut;
    while(!ifile.eof()){
        getline(ifile,line);
        int flag = sscanf(line.c_str(),"%lf %lf %lf %lf",&nsig,&pur,&signif,&cut);
        if(flag != 4) continue;
        cout << line << endl;
        const int NSIG = (int)nsig;
        out.str("");
        out << "ToyData/toyMC_" << NSIG << "_" << pur << "_" << mistag_rate << "_" << _sigma_over_tau << "_" << PoissonFlag << ".root";
        TFile* file = new TFile(out.str().c_str(),"RECREATE");
	for(int jj=0; jj<nsamples;jj++){
            if(!(jj%100) && jj) cout << jj << " trees generated" << endl;
            out.str("");
            out << "ToyTree" << jj;
	    string treename = out.str();
	    TTree* tree = GenerateTree(treename,NSIG,pur);
	    tree->Write();
	}
	file->Close();
    }
    ifile.close();
    return;
}
