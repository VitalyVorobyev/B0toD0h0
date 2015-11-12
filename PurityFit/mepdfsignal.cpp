#include "mepdfsignal.h"

MEPdfSignal::MEPdfSignal(RooRealVar* m_de, RooRealVar* m_mbc, const int mode, const int h0mode, const int b0f){
  cuts = new MyParams();
  de = m_de; mbc = m_mbc;
  m_mode = mode; m_h0mode = h0mode;
  cout << "Sig constructor. mode: " << m_mode << ", h0mode: " << m_h0mode << endl;
  ggflag = m_h0mode == 10 || m_mode>9? true : false;

  if(!ggflag && b0f == -1) b0f == 1;
  InitParams(b0f);
  Fix();
}

void MEPdfSignal::InitParams(const int b0f){
  // * Delta E * //
  cout << "InitParams... ";
  de_offset = new RooRealVar("de_offset","de_offset",0,-0.015,0.015); de_offset->setConstant(kTRUE);

  de0_sig        = new RooRealVar("de0_sig","de0_sig",cuts->get_de0_sig(m_mode,m_h0mode,b0f),-0.2,0.1);
  de0_sig_offset = new RooFormulaVar("de0_sig_offset","de0_sig_offset","@0+@1",RooArgList(*de0_sig,*de_offset));
  de_param_vec.push_back(de0_sig);
  s_de_sig      = new RooRealVar("s_de_sig","s_de_sig",cuts->get_s_de_sig(m_mode,m_h0mode,b0f),0.01,0.5);
  de_param_vec.push_back(s_de_sig);
//  g_de_sig      = new RooGaussian("g_de_sig","g_de_sig",*de,*de0_sig,*s_de_sig);
  g_de_sig      = new RooGaussian("g_de_sig","g_de_sig",*de,*de0_sig_offset,*s_de_sig);

  de0_CBl_sig   = new RooRealVar("de0_CBl_sig","de0_CBl_sig",cuts->get_de0_CBl_sig(m_mode,m_h0mode,b0f),-0.2,0.1); de_param_vec.push_back(de0_CBl_sig);
  de0_CBl_sig_offset = new RooFormulaVar("de0_CBl_sig_offset","de0_CBl_sig_offset","@0+@1",RooArgList(*de0_CBl_sig,*de_offset));
  s_CBl_de_sig  = new RooRealVar("s_CBl_de_sig","s_CBl_de_sig",cuts->get_s_CBl_de_sig(m_mode,m_h0mode,b0f),0.01,0.5); de_param_vec.push_back(s_CBl_de_sig);
  alphal_de_sig = new RooRealVar("alphal_de_sig","alphal_de_sig",cuts->get_alphal_de_sig(m_mode,m_h0mode,b0f),-10,10.); de_param_vec.push_back(alphal_de_sig);
  nl_de_sig     = new RooRealVar("nl_de_sig","nl_de_sig",2.,0.,100.); nl_de_sig->setConstant(kTRUE); de_param_vec.push_back(nl_de_sig);
//  CBl_de_sig    = new RooCBShape("CBl_de_sig","CBl_de_sig",*de,*de0_CBl_sig,*s_CBl_de_sig,*alphal_de_sig,*nl_de_sig);
  CBl_de_sig    = new RooCBShape("CBl_de_sig","CBl_de_sig",*de,*de0_CBl_sig_offset,*s_CBl_de_sig,*alphal_de_sig,*nl_de_sig);

  if(m_mode<20){
    de0_CBr_sig   = new RooRealVar("de0_CBr_sig","de0_CBr_sig",cuts->get_de0_CBr_sig(m_mode,m_h0mode,b0f),-0.2,0.1); de_param_vec.push_back(de0_CBr_sig);
    de0_CBr_sig_offset = new RooFormulaVar("de0_CBr_sig_offset","de0_CBr_sig_offset","@0+@1",RooArgList(*de0_CBr_sig,*de_offset));
    s_CBr_de_sig  = new RooRealVar("s_CBr_de_sig","s_CBr_de_sig",cuts->get_s_CBr_de_sig(m_mode,m_h0mode,b0f),0.01,0.5); de_param_vec.push_back(s_CBr_de_sig);
    alphar_de_sig = new RooRealVar("alphar_de_sig","alphar_de_sig",cuts->get_alphar_de_sig(m_mode,m_h0mode,b0f),-10.,10.); de_param_vec.push_back(alphar_de_sig);
    nr_de_sig     = new RooRealVar("nr_de_sig","nr_de_sig",2,0.,100.); nr_de_sig->setConstant(kTRUE); de_param_vec.push_back(nr_de_sig);
//    CBr_de_sig    = new RooCBShape("CBr_de_sig","CBr_de_sig",*de,*de0_CBr_sig,*s_CBr_de_sig,*alphar_de_sig,*nr_de_sig);
    CBr_de_sig    = new RooCBShape("CBr_de_sig","CBr_de_sig",*de,*de0_CBr_sig_offset,*s_CBr_de_sig,*alphar_de_sig,*nr_de_sig);

    alphar_de_sig->setConstant(kTRUE);
    s_CBr_de_sig->setConstant(kTRUE);
    de0_CBr_sig->setConstant(kTRUE);
  }

  alphal_de_sig->setConstant(kTRUE);
  s_CBl_de_sig->setConstant(kTRUE);
  de0_CBl_sig->setConstant(kTRUE);

  s_de_sig->setConstant(kTRUE);

  f_CBl_de_sig  = new RooRealVar("f_CBl_de_sig","f_CBl_de_sig",cuts->get_f_CBl_de_sig(m_mode,m_h0mode,b0f),0.,1.); de_param_vec.push_back(f_CBl_de_sig);
  f_CBr_de_sig  = new RooRealVar("f_CBr_de_sig","f_CBr_de_sig",cuts->get_f_CBr_de_sig(m_mode,m_h0mode,b0f),0.,1.); de_param_vec.push_back(f_CBr_de_sig);
  f_CBl_de_sig->setConstant(kTRUE);
  f_CBr_de_sig->setConstant(kTRUE);

  if(m_mode<20){
    pdf_de_sig  = new RooAddPdf("pdf_de_sig","pdf_de_sig",RooArgList(*CBl_de_sig,*CBr_de_sig,*g_de_sig),RooArgSet(*f_CBl_de_sig,*f_CBr_de_sig));
  } else{
    pdf_de_sig  = new RooAddPdf("pdf_de_sig","pdf_de_sig",RooArgList(*CBl_de_sig,*g_de_sig),RooArgSet(*f_CBl_de_sig));
  }

  // * Mbc * //
  mbc_offset = new RooRealVar("mbc_offset","mbc_offset",0,-0.005,0.005); mbc_offset->setConstant(kTRUE);

  a_s_mbc_sig   = new RooRealVar("a_s_mbc_sig","a_s_mbc_sig",cuts->get_a_s_mbc_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(a_s_mbc_sig); a_s_mbc_sig->setConstant(kTRUE);
  b_s_mbc_sig   = new RooRealVar("b_s_mbc_sig","b_s_mbc_sig",cuts->get_b_s_mbc_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(b_s_mbc_sig); b_s_mbc_sig->setConstant(kTRUE);
  c_s_mbc_sig   = new RooRealVar("c_s_mbc_sig","c_s_mbc_sig",cuts->get_c_s_mbc_sig(m_mode,m_h0mode),0.0015,0.0035); mbc_param_vec.push_back(c_s_mbc_sig);

  S_mbc_sig     = new RooFormulaVar("S_mbc_sig","S_mbc_sig","@1+@2*@0+@3*@0*@0",RooArgList(*de,*c_s_mbc_sig,*b_s_mbc_sig,*a_s_mbc_sig));

  if(ggflag){
    a_mbc0_sig    = new RooRealVar("a_mbc0_sig","a_mbc0_sig",cuts->get_a_mbc0_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(a_mbc0_sig); a_mbc0_sig->setConstant(kTRUE);
    b_mbc0_sig    = new RooRealVar("b_mbc0_sig","b_mbc0_sig",cuts->get_b_mbc0_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(b_mbc0_sig); b_mbc0_sig->setConstant(kTRUE);
    c_mbc0_sig    = new RooRealVar("c_mbc0_sig","c_mbc0_sig",cuts->get_c_mbc0_sig(m_mode,m_h0mode),5.277,5.285); mbc_param_vec.push_back(c_mbc0_sig);
    MBC0_sig      = new RooFormulaVar("MBC0_sig","MBC0_sig","@1+@2*@0+@3*@0*@0+@4",RooArgList(*de,*c_mbc0_sig,*b_mbc0_sig,*a_mbc0_sig,*mbc_offset));
  } else{
    a_mbc0_sig    = new RooRealVar("a_mbc0_sig","a_mbc0_sig",cuts->get_a_mbc0_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(a_mbc0_sig); a_mbc0_sig->setConstant(kTRUE);
    b_mbc0_sig    = new RooRealVar("b_mbc0_sig","b_mbc0_sig",cuts->get_b_mbc0_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(b_mbc0_sig); b_mbc0_sig->setConstant(kTRUE);
    c_mbc0_sig    = new RooRealVar("c_mbc0_sig","c_mbc0_sig",5.284,5.277,5.285); mbc_param_vec.push_back(c_mbc0_sig);
    d_mbc0_sig    = new RooRealVar("d_mbc0_sig","d_mbc0_sig",cuts->get_c_mbc0_sig(m_mode,m_h0mode)); mbc_param_vec.push_back(d_mbc0_sig); d_mbc0_sig->setConstant(kTRUE);
    MBC0_sig      = new RooFormulaVar("MBC0_sig","MBC0_sig","@0+@1*TMath::Erf((@2-@3))/@4+@5",RooArgList(*c_mbc0_sig,*a_mbc0_sig,*b_mbc0_sig,*de,*d_mbc0_sig,*mbc_offset));
  }

  alpha_mbc_sig = new RooRealVar("alpha_mbc_sig","alpha_mbc_sig",0.139,0.01,2.); mbc_param_vec.push_back(alpha_mbc_sig);

  pdf_mbc_sig   = new RooNovosibirsk("pdf_mbc_sig","pdf_mbc_sig",*mbc,*MBC0_sig,*S_mbc_sig,*alpha_mbc_sig);

  if(!ggflag){
    // * Delta E * //
    de0_sig_tail       = new RooRealVar("de0_sig_tail","de0_sig_tail",cuts->get_de0_sig(m_mode,m_h0mode,5),-0.2,0.1); de_param_vec.push_back(de0_sig_tail);
    de0_sig_tail_offset = new RooFormulaVar("de0_sig_tail_offset","de0_sig_tail_offset","@0+@1",RooArgList(*de0_sig_tail,*de_offset));
    s_de_sig_tail      = new RooRealVar("s_de_sig_tail","s_de_sig_tail",cuts->get_s_de_sig(m_mode,m_h0mode,5),0.,0.5); de_param_vec.push_back(s_de_sig_tail);
//    g_de_sig_tail      = new RooGaussian("g_de_sig_tail","g_de_sig_tail",*de,*de0_sig_tail,*s_de_sig_tail);
    g_de_sig_tail      = new RooGaussian("g_de_sig_tail","g_de_sig_tail",*de,*de0_sig_tail_offset,*s_de_sig_tail);

    de0_CBl_sig_tail   = new RooRealVar("de0_CBl_sig_tail","de0_CBl_sig_tail",cuts->get_de0_CBl_sig(m_mode,m_h0mode,5),-0.1,0.1); de_param_vec.push_back(de0_CBl_sig_tail);
    de0_CBl_sig_tail_offset = new RooFormulaVar("de0_CBl_sig_tail_offset","de0_CBl_sig_tail_offset","@0+@1",RooArgList(*de0_CBl_sig_tail,*de_offset));
    s_CBl_de_sig_tail  = new RooRealVar("s_CBl_de_sig_tail","s_CBl_de_sig_tail",cuts->get_s_CBl_de_sig(m_mode,m_h0mode,5),0.,0.5); de_param_vec.push_back(s_CBl_de_sig_tail);
    alphal_de_sig_tail = new RooRealVar("alphal_de_sig_tail","alphal_de_sig_tail",cuts->get_alphal_de_sig(m_mode,m_h0mode,5),-10.,10.); de_param_vec.push_back(alphal_de_sig_tail);
    nl_de_sig_tail     = new RooRealVar("nl_de_sig_tail","nl_de_sig_tail",2.,0.,100.); nl_de_sig_tail->setConstant(kTRUE); de_param_vec.push_back(nl_de_sig_tail);
//    CBl_de_sig_tail    = new RooCBShape("CBl_de_sig_tail","CBl_de_sig_tail",*de,*de0_CBl_sig_tail,*s_CBl_de_sig_tail,*alphal_de_sig_tail,*nl_de_sig_tail);
    CBl_de_sig_tail    = new RooCBShape("CBl_de_sig_tail","CBl_de_sig_tail",*de,*de0_CBl_sig_tail_offset,*s_CBl_de_sig_tail,*alphal_de_sig_tail,*nl_de_sig_tail);

    fCBl_de_sig_tail   = new RooRealVar("fCBl_de_sig_tail","fCBl_de_sig_tail",cuts->get_f_CBl_de_sig(m_mode,m_h0mode,5),0.,1.); de_param_vec.push_back(fCBl_de_sig_tail);

    pdf_de_sig_tail    = new RooAddPdf("pdf_de_sig_tail","pdf_de_sig_tail",RooArgList(*CBl_de_sig_tail,*g_de_sig_tail),RooArgSet(*fCBl_de_sig_tail));

    // * Mbc * //
    a_mbc0_sig_tail  = new RooRealVar("a_mbc0_sig_tail","a_mbc0_sig_tail",cuts->get_a_mbc0_sig_tail(m_mode)); mbc_param_vec.push_back(a_mbc0_sig_tail);
    a_mbc0_sig_tail->setConstant(kTRUE);
    b_mbc0_sig_tail  = new RooRealVar("b_mbc0_sig_tail","b_mbc0_sig_tail",cuts->get_b_mbc0_sig_tail(m_mode)); mbc_param_vec.push_back(b_mbc0_sig_tail);
    b_mbc0_sig_tail->setConstant(kTRUE);
    c_mbc0_sig_tail  = new RooRealVar("c_mbc0_sig_tail","c_mbc0_sig_tail",cuts->get_c_mbc0_sig_tail(m_mode),5.27,5.29); mbc_param_vec.push_back(c_mbc0_sig_tail);
    MBC0_sig_tail    = new RooFormulaVar("MBC0_sig_tail","MBC0_sig_tail","@1+@2*@0+@3*@0*@0+@4",RooArgList(*de,*c_mbc0_sig_tail,*b_mbc0_sig_tail,*a_mbc0_sig_tail,*mbc_offset));

    a_s_mbc_sig_tail = new RooRealVar("a_s_mbc_sig_tail","a_s_mbc_sig_tail",cuts->get_a_s_mbc_sig_tail(m_mode)); mbc_param_vec.push_back(a_s_mbc_sig_tail);
    a_s_mbc_sig_tail->setConstant(kTRUE);
    b_s_mbc_sig_tail = new RooRealVar("b_s_mbc_sig_tail","b_s_mbc_sig_tail",cuts->get_b_s_mbc_sig_tail(m_mode)); mbc_param_vec.push_back(b_s_mbc_sig_tail);
    b_s_mbc_sig_tail->setConstant(kTRUE);
    c_s_mbc_sig_tail = new RooRealVar("c_s_mbc_sig_tail","c_s_mbc_sig_tail",cuts->get_c_s_mbc_sig_tail(m_mode),0.0015,0.0055); mbc_param_vec.push_back(c_s_mbc_sig_tail);
    S_mbc_sig_tail   = new RooFormulaVar("S_mbc_sig_tail","S_mbc_sig_tail","@1+@2*@0+@3*@0*@0",RooArgList(*de,*c_s_mbc_sig_tail,*b_s_mbc_sig_tail,*a_s_mbc_sig_tail));

    pdf_mbc_sig_tail = new RooNovosibirsk("pdf_mbc_sig_tail","pdf_mbc_sig_tail",*mbc,*MBC0_sig_tail,*S_mbc_sig_tail,*alpha_mbc_sig);

    pdf_peak = new RooProdPdf("pdf_peak","pdf_peak",*pdf_de_sig,Conditional(*pdf_mbc_sig,*mbc));
    pdf_tail = new RooProdPdf("pdf_tail","pdf_tail",*pdf_de_sig_tail,Conditional(*pdf_mbc_sig_tail,*mbc));

    f_tail      = new RooRealVar("f_tail","f_tail",0.3,0.,1.); de_param_vec.push_back(f_tail);
    pdf_sig_ppp = new RooAddPdf("pdf_sig_ppp","pdf_sig_ppp",RooArgSet(*pdf_tail,*pdf_peak),RooArgList(*f_tail));
  } else{
    pdf_sig_gg = new RooProdPdf("pdf_sig","pdf_sig",*pdf_de_sig,Conditional(*pdf_mbc_sig,*mbc));
  }
  cout << "done." << endl;
}

int MEPdfSignal::TryParameters(RooDataSet *ds){
  FixAll();
  if(ggflag) pdf_sig_gg->fitTo(*ds,Verbose(),Timer(true));
  else       pdf_sig_ppp->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  return 0;
}

int MEPdfSignal::FitParameters(RooDataSet *ds){
  FreeAll();
  if(ggflag) pdf_sig_gg->fitTo(*ds,Verbose(),Timer(true));
  else{
    FitTailParameters(ds);
    FitPeakParameters(ds);
    FixAll();
    de0_sig->setConstant(kFALSE);
    s_de_sig->setConstant(kFALSE);
    s_de_sig_tail->setConstant(kFALSE);

    c_mbc0_sig->setConstant(kFALSE);
    c_s_mbc_sig->setConstant(kFALSE);
    c_mbc0_sig_tail->setConstant(kFALSE);
    c_s_mbc_sig_tail->setConstant(kFALSE);

    pdf_sig_ppp->fitTo(*ds,Verbose(),Timer(true));
  }
  Draw(ds);
  PrintParameters();
  return 0;
}

int MEPdfSignal::FitTailParameters(RooDataSet *ds){
  FreeAll();
  RooDataSet* ds0 = (RooDataSet*)ds->reduce(RooArgSet(*mbc,*de),"b0f==5 || b0f == 10");
  pdf_tail->fitTo(*ds0,Verbose(),Timer(true));
  Draw(ds0,true);
  PrintParameters();
  return 0;
}

int MEPdfSignal::FitPeakParameters(RooDataSet *ds){
  FreeAll();
  RooDataSet* ds0 = (RooDataSet*)ds->reduce(RooArgSet(*mbc,*de),"b0f==1");
  pdf_peak->fitTo(*ds0,Verbose(),Timer(true));
  Draw(ds0,2);
  PrintParameters();
  return 0;
}

int MEPdfSignal::FitMbcParameters(RooDataSet *ds){
  Fix();
  FreeMbc();
  if(ggflag) pdf_sig_gg->fitTo(*ds,Verbose(),Timer(true));
  else       pdf_sig_ppp->fitTo(*ds,Verbose(),Timer(true));
  Draw(ds);
  PrintParameters();
  return 0;
}

void MEPdfSignal::PrintParameters(void){
  cout << "Signal PDF parameters for mode " << m_mode << ", h0mode " << m_h0mode << ":" << endl;
  const int NdePar = de_param_vec.size();
  cout << "Delta E (" << NdePar << " parameters):" << endl;
  for(int i=0; i<NdePar; i++){
    de_param_vec[i]->Print();
  }
  const int NmbcPar = mbc_param_vec.size();
  cout << "Mbc (" << NmbcPar << " parameters):" << endl;
  for(int i=0; i<NmbcPar; i++){
    mbc_param_vec[i]->Print();
  }
  return;
}

void MEPdfSignal::WriteParameters(void){
  stringstream out;
  out.str("");
  out << "params/SigParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  cout << "Saving parameters in file " << out.str() << endl;
  ofstream ofile;
  ofile.open(out.str().c_str(),ofstream::out);
  const int NdePar = de_param_vec.size();
  ofile << "Delta E (" << NdePar << " parameters):" << endl;
  for(int i=0; i<NdePar; i++){
    ofile << de_param_vec[i]->format(4,"NE")->Data() << endl;
  }
  const int NmbcPar = mbc_param_vec.size();
  ofile << "Mbc (" << NmbcPar << " parameters):" << endl;
  for(int i=0; i<NmbcPar; i++){
    ofile << mbc_param_vec[i]->format(4,"NE")->Data() << endl;
  }
  ofile.close();
  return;
}

int MEPdfSignal::GetParametersFromFile(void){
  stringstream out;
  cout << "GPFF Sig: mode: " << m_mode << ", h0mode: " << m_h0mode << endl;
  out.str("");
  out << "params/SigParams_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  cout << "Loading parameters in file " << out.str() << endl;
  ifstream ifile;
  ifile.open(out.str().c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "Can't open file " << out.str() << endl;
    return -1;
  }

  string line,name;
  int npars;
  double val;
  char namech[15];
  getline(ifile,line);
  sscanf(line.c_str(),"Delta E (%d parameters):",&npars);
  for(int i=0; i<npars; i++){
    getline(ifile,line);
    cout << line << endl;
    sscanf(line.c_str(),"%s = %lf",namech,&val);
    name = string(namech);
    cout << name << " " << val << endl;
    if(name == string("de0_sig")){       de0_sig->setVal(val); continue;}
    if(name == string("s_de_sig")){      s_de_sig->setVal(val); continue;}
    if(name == string("de0_CBl_sig")){   de0_CBl_sig->setVal(val); continue;}
    if(name == string("s_CBl_de_sig")){  s_CBl_de_sig->setVal(val); continue;}
    if(name == string("alphal_de_sig")){ alphal_de_sig->setVal(val); continue;}
    if(name == string("nl_de_sig")){     nl_de_sig->setVal(val); continue;}
    if(name == string("de0_CBr_sig")){   de0_CBr_sig->setVal(val); continue;}
    if(name == string("s_CBr_de_sig")){  s_CBr_de_sig->setVal(val); continue;}
    if(name == string("alphar_de_sig")){ alphar_de_sig->setVal(val); continue;}
    if(name == string("nr_de_sig")){     nr_de_sig->setVal(val); continue;}
    if(name == string("f_CBl_de_sig")){  f_CBl_de_sig->setVal(val); continue;}
    if(name == string("f_CBr_de_sig")){  f_CBr_de_sig->setVal(val); continue;}

    if(name == string("de0_sig_tail")){       de0_sig_tail->setVal(val); continue;}
    if(name == string("s_de_sig_tail")){      s_de_sig_tail->setVal(val); continue;}
    if(name == string("de0_CBl_sig_tail")){   de0_CBl_sig_tail->setVal(val); continue;}
    if(name == string("s_CBl_de_sig_tail")){  s_CBl_de_sig_tail->setVal(val); continue;}
    if(name == string("alphal_de_sig_tail")){ alphal_de_sig_tail->setVal(val); continue;}
    if(name == string("nl_de_sig_tail")){     nl_de_sig_tail->setVal(val); continue;}
    if(name == string("fCBl_de_sig_tail")){   fCBl_de_sig_tail->setVal(val); continue;}
    if(name == string("f_tail")){             f_tail->setVal(val); continue;}
    cout << "MEPdfSignal::GetParametersFromFile: can't find " << name << endl;
  }
  getline(ifile,line);
  sscanf(line.c_str(),"Mbc (%d parameters):",&npars);
  for(int i=0; i<npars; i++){
    getline(ifile,line);
    cout << line << endl;
    sscanf(line.c_str(),"%s = %lf",namech,&val);
    name = string(namech);
    if(name == string("a_s_mbc_sig")){   a_s_mbc_sig->setVal(val);   continue;}
    if(name == string("b_s_mbc_sig")){   b_s_mbc_sig->setVal(val);   continue;}
    if(name == string("c_s_mbc_sig")){   c_s_mbc_sig->setVal(val);   continue;}
    if(name == string("a_mbc0_sig")){    a_mbc0_sig->setVal(val);    continue;}
    if(name == string("b_mbc0_sig")){    b_mbc0_sig->setVal(val);    continue;}
    if(name == string("c_mbc0_sig")){    c_mbc0_sig->setVal(val);    continue;}
    if(name == string("d_mbc0_sig")){    d_mbc0_sig->setVal(val);    continue;}
    if(name == string("alpha_mbc_sig")){ alpha_mbc_sig->setVal(val); continue;}

    if(name == string("a_mbc0_sig_tail")){  a_mbc0_sig_tail->setVal(val);  continue;}
    if(name == string("b_mbc0_sig_tail")){  b_mbc0_sig_tail->setVal(val);  continue;}
    if(name == string("c_mbc0_sig_tail")){  c_mbc0_sig_tail->setVal(val);  continue;}
    if(name == string("a_s_mbc_sig_tail")){ a_s_mbc_sig_tail->setVal(val); continue;}
    if(name == string("b_s_mbc_sig_tail")){ b_s_mbc_sig_tail->setVal(val); continue;}
    if(name == string("c_s_mbc_sig_tail")){ c_s_mbc_sig_tail->setVal(val); continue;}
    cout << "MEPdfSignal::GetParametersFromFile: can't find " << name << endl;
  }
  cout << "Done." << endl;
  return 0;
}

void MEPdfSignal::ChangeParState(const int state_flag, const bool tail_flag){
  ChangeDeltaEParState(state_flag,tail_flag);
  ChangeMbcParState(state_flag,tail_flag);
  return;
}

void MEPdfSignal::ChangeDeltaEParState(const int state_flag, const bool tail_flag){
  cout << "Changing dE par state " << state_flag << " " << tail_flag << "... ";
  if(!tail_flag){
    de0_sig->setConstant(state_flag);
    s_de_sig->setConstant(state_flag);

    f_CBl_de_sig->setConstant(state_flag);
    de0_CBl_sig->setConstant(state_flag);
    s_CBl_de_sig->setConstant(state_flag);
    alphal_de_sig->setConstant(state_flag);

    if(m_mode<20){
      f_CBr_de_sig->setConstant(state_flag);
      de0_CBr_sig->setConstant(state_flag);
      s_CBr_de_sig->setConstant(state_flag);
      alphar_de_sig->setConstant(state_flag);
    }
  }
  if(!ggflag){// h0 -> pi+pi-pi0
    de0_sig_tail->setConstant(state_flag);
    s_de_sig_tail->setConstant(state_flag);
    de0_CBl_sig_tail->setConstant(state_flag);
    s_CBl_de_sig_tail->setConstant(state_flag);
    alphal_de_sig_tail->setConstant(state_flag);
    fCBl_de_sig_tail->setConstant(state_flag);
  }
  cout << "done." << endl;
}

void MEPdfSignal::ChangeMbcParState(const int state_flag, const bool tail_flag){
  cout << "Changing Mbc par state " << state_flag << " " << tail_flag << "... ";
  if(!tail_flag){
    c_s_mbc_sig->setConstant(state_flag);
    alpha_mbc_sig->setConstant(state_flag);
    c_mbc0_sig->setConstant(state_flag);
  }
  if(!ggflag){// h0 -> pi+pi-pi0
    c_s_mbc_sig_tail->setConstant(state_flag);
    c_mbc0_sig_tail->setConstant(state_flag);
  }
  cout << "done." << endl;
}

void MEPdfSignal::DrawDeltaE(RooDataSet* ds,const int tail){
  mbc->setRange("mbcSignal",cuts->get_mbc_min_h0(m_mode,m_h0mode),cuts->get_mbc_max_h0(m_mode,m_h0mode));
  RooPlot* deFrame = de->frame();
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  if(tail == 1){
    pdf_tail->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  } else if(tail == 2){
    pdf_peak->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  } else{
    if(ggflag) pdf_sig_gg->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
    else      pdf_sig_ppp->plotOn(deFrame,LineWidth(2),LineColor(kGreen));
  }
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
  if(tail == 1){
    pdf_tail->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  } else if(tail == 2){
    pdf_peak->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  } else{
    if(ggflag) pdf_sig_gg->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
    else      pdf_sig_ppp->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
  }

  ds->statOn(deFrame,Layout(0.55,0.98,0.9));
  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull = de->frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E, Signal","#Delta E, Signal",600,700);
  cm->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.04);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  TPaveText *pt = new TPaveText(0.6,0.55,0.98,0.7,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(out.str().c_str());
  pt->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  pt->Draw();

  TLine *de_line_RIGHT = new TLine(cuts->get_de_max_h0(m_mode,m_h0mode),0,cuts->get_de_max_h0(m_mode,m_h0mode),de_line_size());
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT = new TLine(cuts->get_de_min_h0(m_mode,m_h0mode),0,cuts->get_de_min_h0(m_mode,m_h0mode),de_line_size());
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(cuts->get_de_fit_min(),3,cuts->get_de_fit_max(),3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(cuts->get_de_fit_min(),0,cuts->get_de_fit_max(),0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(cuts->get_de_fit_min(),-3,cuts->get_de_fit_max(),-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();
  out.str("");
  if(tail == 1)      out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << "_tail.eps";
  else if(tail == 2) out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << "_peak.eps";
  else               out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cm->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());

  out.str("");
  if(tail == 1)      out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << "_tail.root";
  else if(tail == 2) out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << "_peak.root";
  else               out << "pics/de_sig_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cm->Print(out.str().c_str());

  return;
}

void MEPdfSignal::DrawMbc(RooDataSet* ds, const int tail){
  de->setRange("deSignal",cuts->get_de_min_h0(m_mode,m_h0mode),cuts->get_de_max_h0(m_mode,m_h0mode));
  RooPlot* mbcFrame = mbc->frame();
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),MarkerColor(kGreen));
  if(tail == 1){
    pdf_tail->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  } else if(tail == 2){
    pdf_peak->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  } else{
    if(ggflag) pdf_sig_gg->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
    else      pdf_sig_ppp->plotOn(mbcFrame,LineWidth(2),LineColor(kGreen));
  }
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  if(tail == 1){
    pdf_tail->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  } else if(tail == 2){
    pdf_peak->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  }else {
    if(ggflag) pdf_sig_gg->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
    else      pdf_sig_ppp->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));
  }
  ds->statOn(mbcFrame,Layout(0.2,0.68,0.9));
  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull = mbc->frame(Title("M_{bc} pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("M_{bc}, Signal","M_{bc}, Signal",600,700);
  cmmbc->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);

  mbcFrame->GetXaxis()->SetTitleSize(0.05);
  mbcFrame->GetXaxis()->SetTitleOffset(0.85);
  mbcFrame->GetXaxis()->SetLabelSize(0.04);
  mbcFrame->GetYaxis()->SetTitleOffset(1.6);
  mbcFrame->Draw();

  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  TPaveText *ptmbc = new TPaveText(0.3,0.55,0.68,0.7,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT = new TLine(cuts->get_mbc_max_h0(m_mode,m_h0mode),0,cuts->get_mbc_max_h0(m_mode,m_h0mode),mbc_line_size());
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT = new TLine(cuts->get_mbc_min_h0(m_mode,m_h0mode),0,cuts->get_mbc_min_h0(m_mode,m_h0mode),mbc_line_size());
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();

  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  mbcPull->SetMarkerSize(0.05);
  mbcPull->Draw();
  TLine *mbc_lineUP = new TLine(cuts->get_mbc_fit_min(),3,cuts->get_mbc_fit_max(),3);
  mbc_lineUP->SetLineColor(kBlue);
  mbc_lineUP->SetLineStyle(2);
  mbc_lineUP->Draw();
  TLine *mbc_line = new TLine(cuts->get_mbc_fit_min(),0,cuts->get_mbc_fit_max(),0);
  mbc_line->SetLineColor(kBlue);
  mbc_line->SetLineStyle(1);
  mbc_line->SetLineWidth((Width_t)2.);
  mbc_line->Draw();
  TLine *mbc_lineDOWN = new TLine(cuts->get_mbc_fit_min(),-3,cuts->get_mbc_fit_max(),-3);
  mbc_lineDOWN->SetLineColor(kBlue);
  mbc_lineDOWN->SetLineStyle(2);
  mbc_lineDOWN->Draw();

  cmmbc->Update();
  out.str("");
  if(tail == 1)      out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << "_tail.eps";
  else if(tail == 2) out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << "_peak.eps";
  else               out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << ".eps";
  cmmbc->Print(out.str().c_str());
  string line = string("evince ") + out.str() + string(" &");
  system(line.c_str());

  out.str("");
  if(tail == 1)      out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << "_tail.root";
  else if(tail == 2) out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << "_peak.root";
  else               out << "pics/mbc_sig_m" << m_mode << "_h0m" << m_h0mode << ".root";
  cmmbc->Print(out.str().c_str());

  return;
}

double MEPdfSignal::de_line_size(void){
  switch(m_mode){
  case 1: return 8000;
  case 2:
    if(m_h0mode == 10) return 5000;
    else               return 3000;
  case 3:  return 8000;
  case 5:  return 300;
  case 10: return 600;
  case 20: return 400;
  default:
    return -1;
  }
}

double MEPdfSignal::mbc_line_size(void){
  switch(m_mode){
  case 1:  return 20000;
  case 2:
    if(m_h0mode == 10) return 12000;
    else               return 8000;
  case 3:  return 20000;
  case 5:  return 500;
  case 10: return 1400;
  case 20: return 800;
  default:
    return -1;
  }
}

//void MEPdfSignal::DrawScatterPlot(RooDataSet *ds){
//  TCanvas* ellican = new TCanvas("scatterplot","Scatter Plot Mbc dE",400,400);
//  ellican->cd();
//  out.str("");
//  out << "bdtg>" << BDTG_MIN << " && bdtg<" << BDTG_MAX;
//  out << " && mbc>5.265 && de>" << deMin;
//  out << " && mode == " << mode << " && h0mode == " << h0mode;
//  if(_b0f == 1)      out << " && (b0f == 1 || b0f == 10)";
//  else if(_b0f == 5) out << " && b0f != 0 && b0f != -1 && b0f != 1 && b0f != 10";
//  else               out << " && b0f > 0";
//  cout << "Cuts: " << out.str() << endl;
//  tree->Draw("mbc:de",out.str().c_str());
//  out.str("");
//  out << "../Reports/pics/scatplot_sig_m" << _mode << "_b0f" << _b0f;
//  if(deMin>-0.2) out << "_015" << ".png";
//  if(save_flag) ellican->Print(out.str().c_str());
//}
