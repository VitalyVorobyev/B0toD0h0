#include "CPFit.h"

int init_modes(void){
  cout << "Init modes..." << endl;
  MyParams cuts;
  const int NTot = m_ww_tree->GetEntries();
  modes_set.clear();
  bool flag;
  stringstream out;
  for(int i=0; i<NTot; i++){
    flag = false;
    m_ww_tree->GetEvent(i);
    for(int j=0; j<modes_set.size(); j++){
      if(modes_set[j] == Mode(mm_mode,mm_h0mode)){ flag = true; break;}
    }
    if(flag) continue;
    modes_set.push_back(Mode(mm_mode,mm_h0mode));
    if(!no_bkg && !sigmc){
      fbbVec.push_back(m_fbb);
      fbbErrs.push_back(m_fbb_err);
      fprtVec.push_back(m_fprt);
      fprtErrs.push_back(m_fprt_err);
    } else{
      fbbVec.push_back(0);
      fbbErrs.push_back(0);
      fprtVec.push_back(0);
      fprtErrs.push_back(0);
    }

    m_f_bb_offset.push_back(0);
    m_f_prt_offset.push_back(0);
    vector<double> NsigOff1,NsigOff2;
    m_Nsig_offset.push_back(NsigOff1);
    m_Nbsig_offset.push_back(NsigOff2);

    out.str("");
    out << "mode == " << mm_mode << " && h0mode == " << mm_h0mode;
    vector<double> Nsigs1,Nsigs2;
    vector<double> SigErrs1,SigErrs2;
    NsigVec.push_back(Nsigs1);
    NbsigVec.push_back(Nsigs2);
    NsigErrs.push_back(SigErrs1);
    NbsigErrs.push_back(SigErrs2);
    if(sigmc || no_bkg){
      m_Nsig = m_ww_tree->GetEntries(out.str().c_str());
      m_Nsig_err = 0;
    }
    for(int j=0; j<16; j++){
      m_Nsig_offset[modes_set.size()-1].push_back(0);
      m_Nbsig_offset[modes_set.size()-1].push_back(0);
      const int bin = cuts.bin(j);
      const double Ni  = 0.5*cuts.N(bin,-1,0.338);
      const double Nib = 0.5*cuts.N(bin, 1,0.338);
      NsigVec[modes_set.size()-1].push_back(m_Nsig*Ni);
      NbsigVec[modes_set.size()-1].push_back(m_Nsig*Nib);
      const double error1 = m_Nsig_err*Ni;
      const double error2sq = m_Nsig*Ni*(1.-Ni);
      const double error = sqrt(error1*error1+error2sq);
      NsigErrs[modes_set.size()-1].push_back(error);
      const double error1b = m_Nsig_err*Nib;
      const double error2sqb = m_Nsig*Nib*(1.-Nib);
      const double errorb = sqrt(error1b*error1b+error2sqb);
      NbsigErrs[modes_set.size()-1].push_back(errorb);
    }
  }
  cout << "Modes are initialized: " << endl;
  for(int i=0; i<modes_set.size(); i++){
    cout << "Mode " << modes_set[i] << ":" << endl;
    for(int j=0; j<16; j++){
      cout << " Bin " << cuts.bin(j) << ": Nsig = " << NsigVec[i][j] << " +- " << NsigErrs[i][j];
      cout << ", Nbsig = " << NbsigVec[i][j] << " +- " << NbsigErrs[i][j];
      cout << endl;
    }
  }
  cout << endl;
// DrawSigPredictions();
  return modes_set.size();
}

void nuisance_cpv_fit(MnUserParameterState& pstate){
  init_modes();
  init_wt_err();
  calc_covariance();
  int NPar = 0;
  MnUserParameters upar;

  m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
  m_cos2beta = fabs(TMath::Cos(2.*beta/180.*TMath::Pi()));
  upar.Add(string("btau"),     m_btau,0.01,1.,2.);     NPar++;// 0
  upar.Add(string("dm"),       m_dm,0.01,0.500,0.520); NPar++;// 1
  upar.Add(string("sin2beta"), m_sin2beta,0.1,-5.,5.); NPar++;// 2
  upar.Add(string("cos2beta"), m_cos2beta,0.1,-5.,5.); NPar++;// 3

  upar.Add(string("f_ol_sgl_svd1"), m_f_ol_sgl_svd1,0.1,0.,1.); NPar++;// 4
  upar.Add(string("f_ol_mlt_svd1"), m_f_ol_mlt_svd1,0.1,0.,1.); NPar++;// 5
  upar.Add(string("s_ol_svd1"),     m_s_ol_svd1,1,15.,100.);    NPar++;// 6
  upar.Add(string("f_ol_sgl_svd2"), m_f_ol_sgl_svd2,0.1,0.,1.); NPar++;// 7
  upar.Add(string("f_ol_mlt_svd2"), m_f_ol_mlt_svd2,0.1,0.,1.); NPar++;// 8
  upar.Add(string("s_ol_svd2"),     m_s_ol_svd2,1,15.,100.);    NPar++;// 9

  upar.Add(string("scale1"), m_scale1,0.05,0.5,2.); NPar++;// 10
  upar.Add(string("scale2"), m_scale2,0.05,0.5,2.); NPar++;// 11

  upar.Add(string("K1"),K(1),0.01,0.,1.); NPar++;// 12
  upar.Add(string("K2"),K(2),0.01,0.,1.); NPar++;// 13
  upar.Add(string("K3"),K(3),0.01,0.,1.); NPar++;// 14
  upar.Add(string("K4"),K(4),0.01,0.,1.); NPar++;// 15
  upar.Add(string("K5"),K(5),0.01,0.,1.); NPar++;// 16
  upar.Add(string("K6"),K(6),0.01,0.,1.); NPar++;// 17
  upar.Add(string("K7"),K(7),0.01,0.,1.); NPar++;// 18
  upar.Add(string("K8"),K(8),0.01,0.,1.); NPar++;// 19

  upar.Add(string("Kb1"),K(-1),0.01,0.,1.); NPar++;// 20
  upar.Add(string("Kb2"),K(-2),0.01,0.,1.); NPar++;// 21
  upar.Add(string("Kb3"),K(-3),0.01,0.,1.); NPar++;// 22
  upar.Add(string("Kb4"),K(-4),0.01,0.,1.); NPar++;// 23
  upar.Add(string("Kb5"),K(-5),0.01,0.,1.); NPar++;// 24
  upar.Add(string("Kb6"),K(-6),0.01,0.,1.); NPar++;// 25
  upar.Add(string("Kb7"),K(-7),0.01,0.,1.); NPar++;// 26
  upar.Add(string("Kb8"),K(-8),0.01,0.,1.); NPar++;// 27

  upar.Add(string("C1"),C(1),0.1,-3.,3.); NPar++;// 28
  upar.Add(string("C2"),C(2),0.1,-3.,3.); NPar++;// 29
  upar.Add(string("C3"),C(3),0.1,-3.,3.); NPar++;// 30
  upar.Add(string("C4"),C(4),0.1,-3.,3.); NPar++;// 31
  upar.Add(string("C5"),C(5),0.1,-3.,3.); NPar++;// 32
  upar.Add(string("C6"),C(6),0.1,-3.,3.); NPar++;// 33
  upar.Add(string("C7"),C(7),0.1,-3.,3.); NPar++;// 34
  upar.Add(string("C8"),C(8),0.1,-3.,3.); NPar++;// 35

  upar.Add(string("S1"),S(1),0.1,-3.,3.); NPar++;// 36
  upar.Add(string("S2"),S(2),0.1,-3.,3.); NPar++;// 37
  upar.Add(string("S3"),S(3),0.1,-3.,3.); NPar++;// 38
  upar.Add(string("S4"),S(4),0.1,-3.,3.); NPar++;// 39
  upar.Add(string("S5"),S(5),0.1,-3.,3.); NPar++;// 40
  upar.Add(string("S6"),S(6),0.1,-3.,3.); NPar++;// 41
  upar.Add(string("S7"),S(7),0.1,-3.,3.); NPar++;// 42
  upar.Add(string("S8"),S(8),0.1,-3.,3.); NPar++;// 43

  upar.Add(string("wrtag1"),0.003,0.01,-0.1,0.1); NPar++;// 44
  upar.Add(string("wrtag2"),0.003,0.01,-0.1,0.1); NPar++;// 45

  stringstream out;
  cout << "Modes size: " << modes_set.size() << endl;
  cout << "Npar: " << NPar << endl;
  for(int j=0; j<modes_set.size(); j++){
    out.str("");
    out << "fbb_m" << modes_set[j];
    upar.Add(out.str(),0,0.05,-fbbVec[j],1.-fbbVec[j]); NPar++;// 46+34*j
    out.str("");
    out << "fprt_m" << modes_set[j];
    upar.Add(out.str(),0,0.05,-fprtVec[j],1.-fprtVec[j]); NPar++;// 47+34*j
    for(int i=0; i<16; i++){
      out.str("");
      out << "Nsig_m" << modes_set[j] << "_b" << i;
      upar.Add(out.str(),0,sqrt(NsigVec[j][i]),-NsigVec[j][i],3.*NsigVec[j][i]); NPar++;// 48+34*j+i
      out.str("");
      out << "Nbsig_m" << modes_set[j] << "_b" << i;
      upar.Add(out.str(),0,sqrt(NbsigVec[j][i]),-NbsigVec[j][i],3.*NbsigVec[j][i]); NPar++;// 48+34*j+16+i
    }
  }
  cout << "Npar: " << NPar << endl;

  pdfFcnExt* theFCN = new pdfFcnExt();
  MnMigrad migrad(*theFCN,upar);

//  if(no_bkg){
    migrad.Fix("scale1");
    migrad.Fix("scale2");
//  }

//  migrad.Fix("f_ol_sgl_svd1");
//  migrad.Fix("f_ol_mlt_svd1");
//  migrad.Fix("f_ol_sgl_svd2");
//  migrad.Fix("f_ol_mlt_svd2");

  migrad.Fix("wrtag1");
  migrad.Fix("wrtag2");
  for(int j=0; j<modes_set.size(); j++){
    out.str("");
    out << "fbb_m" << modes_set[j];
    migrad.Fix(out.str().c_str());
    out.str("");
    out << "fprt_m" << modes_set[j];
    migrad.Fix(out.str().c_str());
    for(int i=0; i<16; i++){
      out.str("");
      out << "Nsig_m" << modes_set[j] << "_b" << i;
      migrad.Fix(out.str().c_str());
      out.str("");
      out << "Nbsig_m" << modes_set[j] << "_b" << i;
      migrad.Fix(out.str().c_str());
    }
  }

  migrad.Fix("C1");
  migrad.Fix("C2");
  migrad.Fix("C3");
  migrad.Fix("C4");
  migrad.Fix("C5");
  migrad.Fix("C6");
  migrad.Fix("C7");
  migrad.Fix("C8");

  migrad.Fix("S1");
  migrad.Fix("S2");
  migrad.Fix("S3");
  migrad.Fix("S4");
  migrad.Fix("S5");
  migrad.Fix("S6");
  migrad.Fix("S7");
  migrad.Fix("S8");

  migrad.Fix("dm");
  migrad.Fix("btau");

  migrad.Fix("K1");
  migrad.Fix("K2");
  migrad.Fix("K3");
  migrad.Fix("K4");
  migrad.Fix("K5");
  migrad.Fix("K6");
  migrad.Fix("K7");
  migrad.Fix("K8");

  migrad.Fix("Kb1");
  migrad.Fix("Kb2");
  migrad.Fix("Kb3");
  migrad.Fix("Kb4");
  migrad.Fix("Kb5");
  migrad.Fix("Kb6");
  migrad.Fix("Kb7");
  migrad.Fix("Kb8");

  if(m_mode == 1 || m_mode == 2 || m_mode == 10 || m_mode == 20 || m_gg) migrad.Fix("f_ol_mlt_svd1");
  migrad.Fix("s_ol_svd1");
  migrad.Fix("s_ol_svd2");

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
//    return -1;
  }

  pstate = min.UserState();
  const double Fval = min.Fval();
  const double Edm = min.Edm();
  const int NFcn = min.NFcn();

  if(false){
  MnContours cont(*theFCN,min);
  vector< pair<double,double> > Cnt = cont(2,3);
  TCanvas* c1 = new TCanvas("c1","c1",400,400);
  c1->Draw();
  cout << "Size: " << Cnt.size() << endl;
  double sin_con[20], cos_con[20];
  for(int i=0; i<20; i++){
    sin_con[i] = Cnt[i].first;
    cos_con[i] = Cnt[i].second;
    cout << sin_con[i] <<  " " << cos_con[i] << endl;
  }
  double sin_fit[1] = {pstate.Value(2)};
  double cos_fit[1] = {pstate.Value(3)};

  double sin_true[1] = {upar.Value(2)};
  double cos_true[1] = {upar.Value(3)};

  TMultiGraph* mg = new TMultiGraph("mg","mg");
  TGraph* gr = new TGraph(1,sin_true,cos_true);
  TGraph* gr0 = new TGraph(1,sin_fit,cos_fit);
  TGraph* gr1 = new TGraph(20,sin_con,cos_con);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1.3);
  gr->SetMarkerColor(kRed);
  gr0->SetMarkerStyle(20);
  gr0->SetMarkerSize(1.3);
  gr0->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.1);
  gr1->SetMarkerColor(kBlue);

  mg->Add(gr1);
  mg->Add(gr0);
  mg->Add(gr);
  mg->Draw("ap");
  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetYaxis()->SetRangeUser(0.,1.);
  c1->Update();
  c1->Print("nuis_cont.root");
  c1->Print("nuis_cont.eps");
  system("evince nuis_cont.eps &");
  }

  cout << "Minimization completed!" << endl;
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }
  cout << "LH: " << Fval << ", EDM: " << Edm << ", NDF: " << NFcn << endl;

  return;
}

void fit_gen(MnUserParameterState& pstate){
  int NPar = 0;
  MnUserParameters upar;

  upar.Add(string("btau"),     m_btau,0.01,1.,2.);    NPar++;
  upar.Add(string("dm"),       m_dm,0.01,0.500,0.520); NPar++;
  upar.Add(string("sin2beta"), m_sin2beta,0.1,-5.,5.); NPar++;
  upar.Add(string("cos2beta"), m_cos2beta,0.1,-5.,5.); NPar++;

  upar.Add(string("K1"),K(1),0.01,0.,1.); NPar++;
  upar.Add(string("K2"),K(2),0.01,0.,1.); NPar++;
  upar.Add(string("K3"),K(3),0.01,0.,1.); NPar++;
  upar.Add(string("K4"),K(4),0.01,0.,1.); NPar++;
  upar.Add(string("K5"),K(5),0.01,0.,1.); NPar++;
  upar.Add(string("K6"),K(6),0.01,0.,1.); NPar++;
  upar.Add(string("K7"),K(7),0.01,0.,1.); NPar++;
  upar.Add(string("K8"),K(8),0.01,0.,1.); NPar++;

  upar.Add(string("Kb1"),K(-1),0.01,0.,1.); NPar++;
  upar.Add(string("Kb2"),K(-2),0.01,0.,1.); NPar++;
  upar.Add(string("Kb3"),K(-3),0.01,0.,1.); NPar++;
  upar.Add(string("Kb4"),K(-4),0.01,0.,1.); NPar++;
  upar.Add(string("Kb5"),K(-5),0.01,0.,1.); NPar++;
  upar.Add(string("Kb6"),K(-6),0.01,0.,1.); NPar++;
  upar.Add(string("Kb7"),K(-7),0.01,0.,1.); NPar++;
  upar.Add(string("Kb8"),K(-8),0.01,0.,1.); NPar++;


  upar.Add(string("C1"),C(1),0.1,-1.,1.); NPar++;
  upar.Add(string("C2"),C(2),0.1,-1.,1.); NPar++;
  upar.Add(string("C3"),C(3),0.1,-1.,1.); NPar++;
  upar.Add(string("C4"),C(4),0.1,-1.,1.); NPar++;
  upar.Add(string("C5"),C(5),0.1,-1.,1.); NPar++;
  upar.Add(string("C6"),C(6),0.1,-1.,1.); NPar++;
  upar.Add(string("C7"),C(7),0.1,-1.,1.); NPar++;
  upar.Add(string("C8"),C(8),0.1,-1.,1.); NPar++;

  upar.Add(string("S1"),S(1),0.1,-1.,1.); NPar++;
  upar.Add(string("S2"),S(2),0.1,-1.,1.); NPar++;
  upar.Add(string("S3"),S(3),0.1,-1.,1.); NPar++;
  upar.Add(string("S4"),S(4),0.1,-1.,1.); NPar++;
  upar.Add(string("S5"),S(5),0.1,-1.,1.); NPar++;
  upar.Add(string("S6"),S(6),0.1,-1.,1.); NPar++;
  upar.Add(string("S7"),S(7),0.1,-1.,1.); NPar++;
  upar.Add(string("S8"),S(8),0.1,-1.,1.); NPar++;

  pdfFcnGen* theFCN = new pdfFcnGen();
  MnMigrad migrad(*theFCN,upar);

  migrad.Fix("dm");
//  if(!make_bins_scan){
//    migrad.Fix("btau");
//    migrad.Fix("sin2beta");
//    migrad.Fix("cos2beta");
//  } else{
    if(!no_interf){ migrad.Fix("btau");
    } else{
      migrad.Fix("sin2beta");
      migrad.Fix("cos2beta");
    }
//  }

  migrad.Fix("K1");
  migrad.Fix("K2");
  migrad.Fix("K3");
  migrad.Fix("K4");
  migrad.Fix("K5");
  migrad.Fix("K6");
  migrad.Fix("K7");
  migrad.Fix("K8");

  migrad.Fix("Kb1");
  migrad.Fix("Kb2");
  migrad.Fix("Kb3");
  migrad.Fix("Kb4");
  migrad.Fix("Kb5");
  migrad.Fix("Kb6");
  migrad.Fix("Kb7");
  migrad.Fix("Kb8");

//  if(make_bins_scan){
    migrad.Fix("C1");
    migrad.Fix("C2");
    migrad.Fix("C3");
    migrad.Fix("C4");
    migrad.Fix("C5");
    migrad.Fix("C6");
    migrad.Fix("C7");
    migrad.Fix("C8");

    migrad.Fix("S1");
    migrad.Fix("S2");
    migrad.Fix("S3");
    migrad.Fix("S4");
    migrad.Fix("S5");
    migrad.Fix("S6");
    migrad.Fix("S7");
    migrad.Fix("S8");
//  }

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
//    return -1;
  }

  pstate = min.UserState();
  //min.UserParameters();
  //min.UserCovariance();

  cout << "Minimization completed!" << endl;
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

//  MnContours cont(*theFCN,min,2);
//  vector<pair<double,double> > Cnt = cont(2,3);

  return;
}

int calc_CS(void){
  const int NTot = m_ww_tree->GetEntries();
  const double xi = m_btau*m_dm;
  const double alpha = 1.+xi*xi;
  const double prob_not_change_flavor = (alpha+1.)/(2.*alpha);
  const double C1 = 0.5*(1.+alpha);
  const double C2 = 0.5*(1.-alpha);
  DalitzModel model;
  double intC[8];
  double intS[8];
  double intP[8];
  double intPb[8];
  double intK[8];
  double intKb[8];
  int nevBin[8];
  for(int i=0; i<8; i++){
    intC[i] = 0;
    intS[i] = 0;
    intP[i] = 0;
    intPb[i] = 0;
    nevBin[i] = 0;
  }
  double p,pbar,delta;
  double ppbar_sqrt;
  int index;
  for(int i=0; i<NTot; i++){
    GetEventWW(i);
    if(m_flv_mc>0) continue;
    if(m_flv_mc*m_bin_mc>0) model.PPbarDelta(m_mp_mc,m_mm_mc,p,pbar,delta);
    else                    model.PPbarDelta(m_mm_mc,m_mp_mc,p,pbar,delta);
    if(!(i%10000)) cout << i << " event: " << m_mp_mc << " " << m_mm_mc << " " << p << " " << pbar << " " << delta << endl;
    if(std::isnan(p) || std::isnan(pbar) || std::isnan(delta)){
      cout << i << " - bad event: " << m_mp_mc << " " << m_mm_mc << " " << p << " " << pbar << " " << delta << endl;
      continue;
    }
    index = abs(m_bin_mc)-1;
    if(index<0) continue;
    ppbar_sqrt = sqrt(p*pbar);
    intC[index]  += 0.0001*ppbar_sqrt*TMath::Cos(delta);
    intS[index]  += 0.0001*ppbar_sqrt*TMath::Sin(delta);
    intP[index]  += 0.0001*p;
    intPb[index] += 0.0001*pbar;
    nevBin[index]++;
  }
  cout << "Phases with efficiency correction:" << endl;
  double norm = 0;
  cout << "C1 " << C1 << ", C2 " << C2 << endl;
  for(int i=0; i<8; i++){
//    intK[i]  = intP[i]*C1  + intPb[i]*C2;
//    intKb[i] = intPb[i]*C1 + intP[i]*C2;
    intK[i] = (intP[i]*prob_not_change_flavor-intPb[i]*(1.-prob_not_change_flavor))/(2.*prob_not_change_flavor-1.);
    intKb[i] = (intPb[i]*prob_not_change_flavor-intP[i]*(1.-prob_not_change_flavor))/(2.*prob_not_change_flavor-1.);
//    cout << intP[i] << " " << intPb[i] << " " << intK[i] << " " << intKb[i] << endl;
    norm += intK[i]+intKb[i];

    intC[i] /= sqrt(intK[i]*intKb[i]);
    intS[i] /= sqrt(intK[i]*intKb[i])*(2.*prob_not_change_flavor-1.);
//    intK[i] /= nevBin[i]*0.0001;
//    intKb[i]/= nevBin[i]*0.0001;
    cout << "bin " << i+1 << ", C = " << intC[i] << ", S = " << intS[i] << endl;
  }
  for(int i=0; i<8; i++){
    cout << "bin " << i+1 << ", K = " << intK[i]/norm << ", Kb = " << intKb[i]/norm << endl;
  }
  return 0;
}

int calc_K(void){
  int rawK[8];
  int rawKb[8];
  double newN[8], newN_err[8];
  double newNb[8],newNb_err[8];
  double newK[8], newK_err[8];
  double newKb[8],newKb_err[8];
  const double alpha = 1.+m_btau*m_btau*m_dm*m_dm;
  const double C1 = 0.5*(1.+alpha);
  const double C2 = 0.5*(1.-alpha);

  stringstream out;
  cout << "Calculating N with mode " << get_label(m_mode) << endl;
  const int NTot = m_ww_tree->GetEntries();
  for(int i=0; i<8; i++){
    out.str("");
    out << "-bin_mc*flv_mc == " << i+1;
    rawK[i] = m_ww_tree->Draw("bin",out.str().c_str());
    out.str("");
    out << "-bin_mc*flv_mc == " << -(i+1);
    rawKb[i] = m_ww_tree->Draw("bin",out.str().c_str());
    newN[i]  = (double)rawK[i]/NTot;
    newNb[i] = (double)rawKb[i]/NTot;
    newN_err[i]  = sqrt(newN[i]*(1.-newN[i])/NTot);
    newNb_err[i] = sqrt(newNb[i]*(1.-newNb[i])/NTot);
    cout << " bin " << i+1 << ", N = " << newN[i] << " +- " << newN_err[i];
    cout << ", Nb = " << newNb[i] << " +- " << newNb_err[i] << endl;

    newK[i]  = newN[i]*C1  + newNb[i]*C2;
    newKb[i] = newNb[i]*C1 + newN[i]*C2;
    newK_err[i]  = sqrt(newN_err[i]*C1*newN_err[i]*C1 + newNb_err[i]*C2*newNb_err[i]*C2);
    newKb_err[i] = sqrt(newNb_err[i]*C1*newNb_err[i]*C1 + newN_err[i]*C2*newN_err[i]*C2);
  }
  cout << "K parameters:" << endl;
  for(int i=0; i<8; i++){
    cout << " bin " << i+1 << ", K = " << newK[i]*100. << " +- " << newK_err[i]*100.;
    cout << ", Kb = " << newKb[i]*100. << " +- " << newKb_err[i]*100. << endl;
  }

  double KK[8],KK_err[8],KbKb[8],KbKb_err[8];
  double KKb[8],KKb_err[8];
  double bins[8],bins_err[8];

  cout << "Efficiency map:" << endl;
  for(int i=0; i<8; i++){
    KK[i] = newK[i]/Karr_model[i];
    KbKb[i] = newKb[i]/Kbarr_model[i];
    KK_err[i] = newK_err[i]/Karr_model[i];
    KbKb_err[i] = newKb_err[i]/Kbarr_model[i];

    KKb[i] = KK[i]/KbKb[i];
    KKb_err[i] = sqrt(KK_err[i]*KK_err[i]/(KK[i]*KK[i])+KKb[i]*KKb[i]*KbKb_err[i]*KbKb_err[i]/(KbKb[i]*KbKb[i]));

    bins[i] = i+1; bins_err[i] = 0;
    cout << " bin " << i+1;
    cout << ", KK   = " << KK[i] << " + - " << KK_err[i];
    cout << ", KbKb = " << KbKb[i] << " + - " << KbKb_err[i];
    cout << ", KKb  = " << KKb[i] << " + - " << KKb_err[i];
    cout << endl;
  }

  TCanvas* c1 = new TCanvas("c1","c1",400,400);
  c1->cd();
  TGraphErrors* gr1 = new TGraphErrors(8,bins,KKb,bins_err,KKb_err);
  out.str("");
  out << "Binned Dalitz plot efficiency asymmetry" << get_label(m_mode);
  gr1->SetTitle(out.str().c_str());
  c1->SetGrid();
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerColor(kBlue);
  gr1->GetXaxis()->SetTitle("Dalitz bin");
  gr1->GetXaxis()->SetTitleSize(0.06);
  gr1->GetXaxis()->SetLabelSize(0.06);
  gr1->GetYaxis()->SetLabelSize(0.06);
  gr1->GetXaxis()->SetTitleOffset(0.8);
  gr1->GetYaxis()->SetRangeUser(0.8,1.2);
  gr1->Draw("ap");
  out.str("");
  out << "pics/asym_m" << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
  const string str1 = out.str();
  out.str("");
  out << str1 << ".root";
  c1->Print(out.str().c_str());
  out.str("");
  out << str1 << ".eps";
  c1->Print(out.str().c_str());
  out.str("");
  out << "evince " << str1 << ".eps &";
  system(out.str().c_str());

  TCanvas* c2 = new TCanvas("c2","c2",400,400);
  c2->cd();
  c2->SetGrid();
  out.str("");
  out << "K_{i}/K_{i}^{model}" << get_label(m_mode);
  TMultiGraph* mg = new TMultiGraph("mg",out.str().c_str());

  TGraphErrors* gr2 = new TGraphErrors(8,bins,KK,bins_err,KK_err);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerColor(kBlue);

  TGraphErrors* gr3 = new TGraphErrors(8,bins,KbKb,bins_err,KbKb_err);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerColor(kRed);

  mg->Add(gr2);
  mg->Add(gr3);

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("Dalitz bin");
  mg->GetXaxis()->SetTitleSize(0.06);
  mg->GetXaxis()->SetLabelSize(0.06);
  mg->GetYaxis()->SetLabelSize(0.06);
  mg->GetXaxis()->SetTitleOffset(0.8);
  mg->GetYaxis()->SetRangeUser(0.8,1.2);
  gPad->Modified();
  out.str("");
  out << "pics/effit_m" << Mode(m_mode) << "_h0m" << h0Mode(m_mode);
  const string str2 = out.str();
  out.str("");
  out << str2 << ".root";
  c2->Print(out.str().c_str());
  out.str("");
  out << str2 << ".eps";
  c2->Print(out.str().c_str());
  out.str("");
  out << "evince " << str2 << ".eps &";
  system(out.str().c_str());

  return 0;
}

//int make_single_fit(MnUserParameterState& pstate, const int flv, const int bin, const int SetNum = 0, const bool plot_flag = false){
//  int NPar = 0;
//  MnUserParameters upar;
//  upar.Add(string("btau"),m_btau,0.2*m_btau,0.,10*m_btau); NPar++;
//  upar.Add(string("dm"),  m_dm,  0.1*m_dm,  0.,5.*m_dm);   NPar++;
//  upar.Add(string("A"),   m_A,   0.1,      -3.,3.);        NPar++;
//  upar.Add(string("B"),   m_B,   0.1,      -3.,3.);        NPar++;

////  if(m_NBkgTot && m_svd != 2) m_pdf_back_svd1->Set_f_otlr(0);
////  if(m_NBkgTot && m_svd != 1) m_pdf_back_svd2->Set_f_otlr(0);

//  pdfFcnSingle* theFCN = new pdfFcnSingle();
//  MnMigrad migrad(*theFCN,upar);
//  if(fix_btau || true)   migrad.Fix("btau");
//  if(fix_dm   || true)   migrad.Fix("dm");
//  if(fix_A)      migrad.Fix("A");
//  if(fix_B)      migrad.Fix("B");
//  FunctionMinimum min = migrad();
//  if(!min.IsValid()) return -1;
//  pstate = min.UserState();
//  for(int i=0; i<NPar; i++){
//    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
//    upar.SetValue(i,pstate.Value(i));
//    cout << pstate.Value(i) << " +- " << pstate.Error(i) << endl;
//  }

//  if(plot_flag){
//    const double& btau = upar.Value(0);
//    const double& dm   = upar.Value(1);
//    const double& A    = upar.Value(2);
//    const double& B    = upar.Value(3);
//    cout << "Set parameters: " << btau << ", " << dm << ", " << A << ", " << B << endl;
//    if(m_svd != 2){
//      m_pdf_svd1->SetTauDm(btau,dm);
//      m_pdf_svd1->SetAB(A,B);
//    }
//    if(m_svd != 1){
//      m_pdf_svd2->SetTauDm(btau,dm);
//      m_pdf_svd2->SetAB(A,B);
//    }
//    draw_single_fit(SetNum);
//  }
//  return 0;
//}

int lifetime_wide_window(MnUserParameterState& pstate, const int smpl = 0){
  int NPar = 0;
  MnUserParameters upar;
  m_f_cont_in_comb = get_f_cont_in_bkg_sig(Mode(m_mode),h0Mode(m_mode));
  upar.Add(string("btau"),m_btau,0.2*m_btau,0, 10*m_btau);      NPar++;
  upar.Add(string("f_ol_sgl_svd1"), m_f_ol_sgl_svd1,0.1,0.,1.); NPar++;
  upar.Add(string("f_ol_mlt_svd1"), m_f_ol_mlt_svd1,0.1,0.,1.); NPar++;
  upar.Add(string("s_ol_svd1"),     m_s_ol_svd1,1,15.,100.);    NPar++;
  upar.Add(string("f_ol_sgl_svd2"), m_f_ol_sgl_svd2,0.1,0.,1.); NPar++;
  upar.Add(string("f_ol_mlt_svd2"), m_f_ol_mlt_svd2,0.1,0.,1.); NPar++;
  upar.Add(string("s_ol_svd2"),     m_s_ol_svd2,1,15.,100.);    NPar++;

  SetBkgScales(m_scale1,m_scale2);

  int NFreePar = NPar;
  pdfFcnWW* theFCN = new pdfFcnWW(smpl);
  MnMigrad migrad(*theFCN,upar);

  if(m_mode == 1 || m_mode == 2 || m_mode == 10 || m_mode == 20 || m_gg) migrad.Fix("f_ol_mlt_svd1"); NFreePar--;
  migrad.Fix("s_ol_svd1"); NFreePar--;
  migrad.Fix("s_ol_svd2"); NFreePar--;

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
  }
//  if(nega_pdf_flag){
//    cout << "nega_pdf_flag!" << endl;
//    nega_pdf_flag = false;
//    return -1;
//  }
  pstate = min.UserState();
  if(NFreePar) cout << "Minimization completed!" << endl;
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  if(NFreePar) cout << "Initializing PDF with optimized parameters." << endl;
  SetPDFParamsWW(upar.Params());

//  cout << "Drawing..." << endl;
//  Draw_NoTag(NFreePar,SetNum);

  return 0;
}

int sideband_wide_window(MnUserParameterState& pstate){
  int NPar = 0;
  MnUserParameters upar;
  m_f_cont_in_comb = get_f_cont_in_bkg_sideband(Mode(m_mode),h0Mode(m_mode));
  upar.Add(string("scale_svd1"),    1,0.2,0.8,1.2);    NPar++;
  upar.Add(string("scale_svd2"),    1,0.2,0.8,2.2);    NPar++;
  upar.Add(string("f_ol_sgl_svd1"), m_f_ol_sgl_svd1,0.1,0.,1.); NPar++;
  upar.Add(string("f_ol_mlt_svd1"), m_f_ol_mlt_svd1,0.1,0.,1.); NPar++;
//  upar.Add(string("s_ol_svd1"),     m_s_ol_svd1,1,15.,100.);    NPar++;
  upar.Add(string("f_ol_sgl_svd2"), m_f_ol_sgl_svd2,0.1,0.,1.); NPar++;
  upar.Add(string("f_ol_mlt_svd2"), m_f_ol_mlt_svd2,0.1,0.,1.); NPar++;
//  upar.Add(string("s_ol_svd2"),     m_s_ol_svd2,1,15.,100.);    NPar++;

  int NFreePar = NPar;
  pdfFcnSBWW* theFCN = new pdfFcnSBWW();
  MnMigrad migrad(*theFCN,upar);
  if(m_mode == 1 || m_mode == 2 || m_mode == 10 || m_mode == 20) migrad.Fix("f_ol_mlt_svd1"); NFreePar--;
//  migrad.Fix("s_ol_svd1"); NFreePar--;
//  migrad.Fix("s_ol_svd2"); NFreePar--;

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
//    return -1;
  }

  pstate = min.UserState();
  if(NFreePar) cout << "Minimization completed!" << endl;
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  if(NFreePar) cout << "Initializing PDF with optimized parameters." << endl;
  SetSBPDFParamsWW(upar.Params());

  return 0;
}

int cpv_wide_window(MnUserParameterState& pstate, const int smpl=0){
  int NPar = 0;
  MnUserParameters upar;
//  m_f_cont_in_comb = get_f_cont_in_bkg_sig(Mode(m_mode),h0Mode(m_mode));
  m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
  m_cos2beta = fabs(TMath::Cos(2.*beta/180.*TMath::Pi()));
  upar.Add(string("sin2beta"),      m_sin2beta,0.1,-5.,5.);      NPar++;
  upar.Add(string("cos2beta"),      m_cos2beta,0.1,-5.,5.);      NPar++;
  upar.Add(string("f_ol_sgl_svd1"), m_f_ol_sgl_svd1,0.1,0.,0.1); NPar++;
  upar.Add(string("f_ol_mlt_svd1"), m_f_ol_mlt_svd1,0.1,0.,0.1); NPar++;
  upar.Add(string("s_ol_svd1"),     m_s_ol_svd1,1,15.,100.);     NPar++;
  upar.Add(string("f_ol_sgl_svd2"), m_f_ol_sgl_svd2,0.1,0.,0.1); NPar++;
  upar.Add(string("f_ol_mlt_svd2"), m_f_ol_mlt_svd2,0.1,0.,0.1); NPar++;
  upar.Add(string("s_ol_svd2"),     m_s_ol_svd2,1,15.,100.);     NPar++;

  SetBkgScales(m_scale1,m_scale2);

  int NFreePar = NPar;
  pdfFcnCPVWW* theFCN = new pdfFcnCPVWW(smpl);
  MnMigrad migrad(*theFCN,upar);

  if(m_mode == 1 || m_mode == 2 || m_mode == 10 || m_mode == 20 || m_gg){
    migrad.Fix("f_ol_mlt_svd1"); NFreePar--;
    //migrad.Fix("f_ol_mlt_svd2"); NFreePar--;
  } else{
//    migrad.Fix("f_ol_sgl_svd1"); NFreePar--;
  }
  migrad.Fix("s_ol_svd1"); NFreePar--;
  migrad.Fix("s_ol_svd2"); NFreePar--;

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
//    return -1;
  }

  pstate = min.UserState();
  if(NFreePar) cout << "Minimization completed!" << endl;
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  if(NFreePar) cout << "Initializing PDF with optimized parameters." << endl;
  SetPDFParamsCPVWW(upar.Params());

  cout << "Drawing..." << endl;
  if(!m_toyfit && !make_bins_scan && !m_nuisance && !m_line_test) DrawCPV(NFreePar);

  return 0;
}

//int make_full_fit(MnUserParameterState& pstate, const int SetNum = 0,const bool draw_flag = true){
//  int NPar = 0;
////  make_sideband_fit();
//  m_f_cont = 0.9;
//  m_sin2beta = TMath::Sin(2.*beta/180.*TMath::Pi());
//  m_cos2beta = TMath::Cos(2.*beta/180.*TMath::Pi());
//  MnUserParameters upar;
//  upar.Add(string("btau"),    m_btau,    0.2*m_btau,0, 10*m_btau); NPar++;
//  upar.Add(string("dm"),      m_dm,      0.1*m_dm,  0.,5.*m_dm);   NPar++;
//  upar.Add(string("sin2beta"),m_sin2beta,0.1,      -5.,5.);        NPar++;
//  upar.Add(string("cos2beta"),m_cos2beta,0.1,      -5.,5.);        NPar++;

//  upar.Add(string("f_ol_sgl_svd1"), m_f_ol_sgl_svd1,0.1,0.,1.); NPar++;
//  upar.Add(string("f_ol_mlt_svd1"), m_f_ol_mlt_svd1,0.1,0.,1.); NPar++;
//  upar.Add(string("s_ol_svd1"),     m_s_ol_svd1,1,15.,100.);    NPar++;
//  upar.Add(string("f_ol_sgl_svd2"), m_f_ol_sgl_svd2,0.1,0.,1.); NPar++;
//  upar.Add(string("f_ol_mlt_svd2"), m_f_ol_mlt_svd2,0.1,0.,1.); NPar++;
//  upar.Add(string("s_ol_svd2"),     m_s_ol_svd2,1,15.,100.);    NPar++;

//  int NFreePar = NPar;
//  pdfFcn* theFCN = new pdfFcn(SetNum);
//  MnMigrad migrad(*theFCN,upar);

//  migrad.Fix("s_ol_svd1"); NFreePar--;
//  migrad.Fix("s_ol_svd2"); NFreePar--;

//  if(m_svd == 2 || SetNum){
//    migrad.Fix("f_ol_mlt_svd1"); NFreePar--;
//    migrad.Fix("f_ol_sgl_svd1"); NFreePar--;
//  } else{
//    if(m_mode<3 || m_mode>9 || sgl_asc){ migrad.Fix("f_ol_mlt_svd1"); NFreePar--;}
//    else{                    migrad.Fix("f_ol_sgl_svd1"); NFreePar--;}
//  }
//  if(m_svd == 1 || SetNum){
//    migrad.Fix("f_ol_sgl_svd2"); NFreePar--;
//    migrad.Fix("f_ol_mlt_svd2"); NFreePar--;
//  } else{
//    if(sgl_asc){ migrad.Fix("f_ol_mlt_svd2"); NFreePar--;}
//    if(mlt_asc){ migrad.Fix("f_ol_sgl_svd2"); NFreePar--;}
//  }

//  if(fix_btau){                  migrad.Fix("btau");     NFreePar--;}
//  if(fix_dm       || no_interf){ migrad.Fix("dm");       NFreePar--;}
//  if(fix_sin2beta || no_interf){ migrad.Fix("sin2beta"); NFreePar--;}
//  if(fix_cos2beta || no_interf){ migrad.Fix("cos2beta"); NFreePar--;}

//  if(NFreePar) cout << "Starting minimization with " << NFreePar << " free parameters" << endl;
//  FunctionMinimum min = migrad();
//  if(!min.IsValid()){
//    cout << "Fit is not valid" << endl;
//    return -1;
//  }
////  if(nega_pdf_flag){
////    cout << "nega_pdf_flag!" << endl;
////    nega_pdf_flag = false;
////    return -1;
////  }
//  pstate = min.UserState();
//  if(NFreePar) cout << "Minimization completed!" << endl;
//  for(int i=0; i<NPar; i++){
//    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
//    upar.SetValue(i,pstate.Value(i));
//    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
//  }

//  if(NFreePar) cout << "Initializing PDF with optimized parameters." << endl;
//  SetPDFParams(upar.Params());

//  if(draw_flag){
//    cout << "Drawing..." << endl;
//    if(draw_b_bbar){
//      Draw_BBbar_CPV(SetNum);
//    } else if(no_interf){
//      Draw_NoTag(NFreePar,SetNum);
//    } else{
//      Draw_All(NFreePar,SetNum);
//    }
//  }
//  return 0;
//}

int bins_scan_wide_window(void){
  m_fitflv = 0;
  double val_vec[2][8];// [0] -> sin
  double err_vec[2][8];// [1] -> cos
  for(int i=0; i<8; i++){
    for(int j=0; j<2; j++){
      val_vec[j][i] = 0;
      err_vec[j][i] = 0;
    }
  }
  const int sin_index = m_genfit ? 2 : 0;
  const int cos_index = m_genfit ? 3 : 1;
  for(int i=0; i<8; i++){
    m_fitbin = i+1;
    GetGoodWWTTree(m_tree,m_ww_tree);
    MnUserParameterState pstate;
    if(!m_genfit) cpv_wide_window(pstate);
    else          fit_gen(pstate);

    val_vec[0][i] = (pstate.Params().at(sin_index) - m_sin2beta)*100.;// sin;
    err_vec[0][i] = pstate.Errors().at(sin_index)*100.;//
    val_vec[1][i] = (pstate.Params().at(cos_index) - m_cos2beta)*100.;// cos;
    err_vec[1][i] = pstate.Errors().at(cos_index)*100.;//
  }
  for(int i=0; i<8; i++){
    cout         << val_vec[0][i] << " +- " << err_vec[0][i];
    cout << ", " << val_vec[1][i] << " +- " << err_vec[1][i] << endl;
  }
  DrawBinsScanCPV(val_vec,err_vec);
  return 0;
}

int bins_scan_wide_lifetime(void){
  m_fitflv = 0;
  double val_vec[8];
  double err_vec[8];
  for(int i=0; i<8; i++){
    val_vec[i] = 0;
    err_vec[i] = 0;
  }
  for(int i=0; i<8; i++){
    m_fitbin = i+1;
    GetGoodWWTTree(m_tree,m_ww_tree);
    MnUserParameterState pstate;
    if(!m_genfit) lifetime_wide_window(pstate);
    else          fit_gen(pstate);
    val_vec[i] = (pstate.Params().at(0) - m_btau)*100.;
    err_vec[i] = pstate.Errors().at(0)*100.;
  }
  for(int i=0; i<8; i++){
    cout << "Bin " << i+1 << ": tau = " << val_vec[i] << " +- " << err_vec[i];
  }
  DrawBinsScanLifetime(val_vec,err_vec);
  return 0;
}

//int bins_scan(void){
//  m_fitflv = 0;
//  double val_vec[3][8];
//  double err_vec[3][8];
//  for(int i=0; i<8; i++){
//    for(int j=0; j<3; j++){
//      val_vec[j][i] = 0;
//      err_vec[j][i] = 0;
//    }
//  }
//  MnUserParameterState pstate;
//  for(int i=0; i<8; i++){
//    m_fitbin = i+1;
//    GetGoodTTree(m_tree,m_mode,m_fitbin,m_fitflv);
//    if(make_full_fit(pstate,0,false)) continue;
//    val_vec[0][i] = pstate.Params().at(0) - m_btau;// tau;
//    err_vec[0][i] = pstate.Errors().at(0);//
//    val_vec[1][i] = pstate.Params().at(2) - m_sin2beta;// sin;
//    err_vec[1][i] = pstate.Errors().at(2);//
//    val_vec[2][i] = pstate.Params().at(3) - m_cos2beta;// cos;
//    err_vec[2][i] = pstate.Errors().at(3);//
//  }
////  DrawBinsScan(val_vec,err_vec);
//  return 0;
//}

int make_sb_test_fit(void){
  m_f_cont = 0.4293;
  cout << "Sideband Test Fit" << endl;
  int NPar = 0;
  double f_delta_mlt_svd1, f_delta_mlt_svd2;
  double f_delta_sgl_svd1, f_delta_sgl_svd2;
  double f_otlr_svd1, f_otlr_svd2;
  MnUserParameters upar;
  if(m_svd != 2){
    f_delta_mlt_svd1 = m_pdf_back_svd1->Get_f_delta_mlt();
    f_delta_sgl_svd1 = m_pdf_back_svd1->Get_f_delta_sgl();
    f_otlr_svd1 = m_pdf_back_svd1->Get_f_otlr();
    cout << "f_delta_mlt_svd1: " << f_delta_mlt_svd1 << endl;
    cout << "f_delta_sgl_svd1: " << f_delta_sgl_svd1 << endl;
    cout << "f_otlr_svd1:      " << f_otlr_svd1 << endl;
  }
  if(m_svd != 1){
    f_delta_mlt_svd2 = m_pdf_back_svd2->Get_f_delta_mlt();
    f_delta_sgl_svd2 = m_pdf_back_svd2->Get_f_delta_sgl();
    f_otlr_svd2 = m_pdf_back_svd2->Get_f_otlr();
    cout << "f_delta_mlt_svd2: " << f_delta_mlt_svd2 << endl;
    cout << "f_delta_sgl_svd2: " << f_delta_sgl_svd2 << endl;
    cout << "f_otlr_svd2:      " << f_otlr_svd2 << endl;
  }
  upar.Add(string("f_delta_sgl_svd1"), f_delta_sgl_svd1, 0.2, 0.,1.); NPar++;
  upar.Add(string("f_delta_mlt_svd1"), f_delta_mlt_svd1, 0.2, 0.,1.); NPar++;
  upar.Add(string("f_otlr_svd1"),      f_otlr_svd1,      0.2, 0.,1.); NPar++;

  upar.Add(string("f_delta_sgl_svd2"), f_delta_sgl_svd2, 0.2, 0.,1.); NPar++;
  upar.Add(string("f_delta_mlt_svd2"), f_delta_mlt_svd2, 0.2, 0.,1.); NPar++;
  upar.Add(string("f_otlr_svd2"),      f_otlr_svd2,      0.2, 0.,1.); NPar++;

  pdfFcnBkg* theFCN = new pdfFcnBkg(1);
  MnMigrad migrad(*theFCN,upar);

  int NFreePar = NPar;
  if(!add_otlr){
    migrad.Fix("f_otlr_svd1"); NFreePar--;
    migrad.Fix("f_otlr_svd2"); NFreePar--;
  }
//  if(sgl_asc){
    migrad.Fix("f_delta_mlt_svd1"); NFreePar--;
    migrad.Fix("f_delta_mlt_svd2"); NFreePar--;
//  }
//  if(mlt_asc){
    migrad.Fix("f_delta_sgl_svd1"); NFreePar--;
    migrad.Fix("f_delta_sgl_svd2"); NFreePar--;
//  }
  if(m_svd == 1){
    if(add_otlr){ migrad.Fix("f_otlr_svd2"); NFreePar--;}
    if(!sgl_asc){ migrad.Fix("f_delta_mlt_svd2"); NFreePar--;}
    if(!mlt_asc){ migrad.Fix("f_delta_sgl_svd2"); NFreePar--;}
  }
  if(m_svd == 2){
    if(add_otlr){ migrad.Fix("f_otlr_svd1"); NFreePar--;}
    if(!sgl_asc){ migrad.Fix("f_delta_mlt_svd1"); NFreePar--;}
    if(!mlt_asc){ migrad.Fix("f_delta_sgl_svd1"); NFreePar--;}
  }

  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }
  SetSBTestPdfParams(upar.Params());
  DrawSBTest(NFreePar,true);

  return 0;
}

int make_sideband_fit(void){
  int NPar = 0;
  cout << "Init params" << endl;
  double tau_svd1, mu_svd1, mu_delta_svd1, f_delta_mlt_svd1, f_tail_mlt_svd1, S_main_mlt_svd1, S_tail_mlt_svd1, f_delta_sgl_svd1, f_tail_sgl_svd1, S_main_sgl_svd1, S_tail_sgl_svd1, f_otlr_svd1, s_otlr_svd1;
  double tau_svd2, mu_svd2, mu_delta_svd2, f_delta_mlt_svd2, f_tail_mlt_svd2, S_main_mlt_svd2, S_tail_mlt_svd2, f_delta_sgl_svd2, f_tail_sgl_svd2, S_main_sgl_svd2, S_tail_sgl_svd2, f_otlr_svd2, s_otlr_svd2;
  MnUserParameters upar;
  if(m_svd != 2){
    tau_svd1         = m_pdf_back_svd1->GetTau();
    mu_svd1          = m_pdf_back_svd1->Get_mu();
    mu_delta_svd1    = m_pdf_back_svd1->Get_mu_delta();
    f_delta_mlt_svd1 = m_pdf_back_svd1->Get_f_delta_mlt();
    f_tail_mlt_svd1  = m_pdf_back_svd1->Get_f_tail_mlt();
    S_main_mlt_svd1  = m_pdf_back_svd1->Get_S_main_mlt();
    S_tail_mlt_svd1  = m_pdf_back_svd1->Get_S_tail_mlt();
    f_delta_sgl_svd1 = m_pdf_back_svd1->Get_f_delta_sgl();
    f_tail_sgl_svd1  = m_pdf_back_svd1->Get_f_tail_sgl();
    S_main_sgl_svd1  = m_pdf_back_svd1->Get_S_main_sgl();
    S_tail_sgl_svd1  = m_pdf_back_svd1->Get_S_tail_sgl();
    f_otlr_svd1      = m_pdf_back_svd1->Get_f_otlr();
    s_otlr_svd1      = m_pdf_back_svd1->Get_s_otlr();
  } else{
    tau_svd1         = m_pdf_back_svd2->GetTau();
    mu_svd1          = m_pdf_back_svd2->Get_mu();
    mu_delta_svd1    = m_pdf_back_svd2->Get_mu_delta();
    f_delta_mlt_svd1 = m_pdf_back_svd2->Get_f_delta_mlt();
    f_tail_mlt_svd1  = m_pdf_back_svd2->Get_f_tail_mlt();
    S_main_mlt_svd1  = m_pdf_back_svd2->Get_S_main_mlt();
    S_tail_mlt_svd1  = m_pdf_back_svd2->Get_S_tail_mlt();
    f_delta_sgl_svd1 = m_pdf_back_svd2->Get_f_delta_sgl();
    f_tail_sgl_svd1  = m_pdf_back_svd2->Get_f_tail_sgl();
    S_main_sgl_svd1  = m_pdf_back_svd2->Get_S_main_sgl();
    S_tail_sgl_svd1  = m_pdf_back_svd2->Get_S_tail_sgl();
    f_otlr_svd1      = m_pdf_back_svd2->Get_f_otlr();
    s_otlr_svd1      = m_pdf_back_svd2->Get_s_otlr();
  }
  upar.Add(string("tau"),             tau_svd1,         0.2, 0.01,10); NPar++;
//  upar.Add(string("tau_svd1"),         tau_svd1,         0.2, 0.01,10); NPar++;
  upar.Add(string("mu_svd1"),          mu_svd1,          0.2,-0.2,0.2); NPar++;
  upar.Add(string("mu_delta_svd1"),    mu_delta_svd1,    0.2,-0.2,0.2); NPar++;
  upar.Add(string("f_delta_mlt_svd1"), f_delta_mlt_svd1, 0.2, 0.,1.);  NPar++;
  upar.Add(string("f_tail_mlt_svd1"),  f_tail_mlt_svd1,  0.2, 0.,1.);  NPar++;
  upar.Add(string("S_main_mlt_svd1"),  S_main_mlt_svd1,  0.2, 0.5,2.); NPar++;
  upar.Add(string("S_tail_mlt_svd1"),  S_tail_mlt_svd1,  0.2, 2. ,8.); NPar++;
  upar.Add(string("f_delta_sgl_svd1"), f_delta_sgl_svd1, 0.2, 0.,1.);  NPar++;
  upar.Add(string("f_tail_sgl_svd1"),  f_tail_sgl_svd1,  0.2, 0.,1.);  NPar++;
  upar.Add(string("S_main_sgl_svd1"),  S_main_sgl_svd1,  0.2, 0.5 ,2.); NPar++;
  upar.Add(string("S_tail_sgl_svd1"),  S_tail_sgl_svd1,  0.2, 2. ,8.); NPar++;
  upar.Add(string("f_otlr_svd1"),      f_otlr_svd1,      0.2, 0.,1.);   NPar++;
  upar.Add(string("s_otlr_svd1"),      s_otlr_svd1,      0.2,15.,70.);  NPar++;

  if(m_svd != 1){
    tau_svd2 =         m_pdf_back_svd2->GetTau();
    mu_svd2  =         m_pdf_back_svd2->Get_mu();
    mu_delta_svd2 =    m_pdf_back_svd2->Get_mu_delta();
    f_delta_mlt_svd2 = m_pdf_back_svd2->Get_f_delta_mlt();
    f_tail_mlt_svd2 =  m_pdf_back_svd2->Get_f_tail_mlt();
    S_main_mlt_svd2 =  m_pdf_back_svd2->Get_S_main_mlt();
    S_tail_mlt_svd2 =  m_pdf_back_svd2->Get_S_tail_mlt();
    f_delta_sgl_svd2 = m_pdf_back_svd2->Get_f_delta_sgl();
    f_tail_sgl_svd2  = m_pdf_back_svd2->Get_f_tail_sgl();
    S_main_sgl_svd2 =  m_pdf_back_svd2->Get_S_main_sgl();
    S_tail_sgl_svd2 =  m_pdf_back_svd2->Get_S_tail_sgl();
    f_otlr_svd2 =      m_pdf_back_svd2->Get_f_otlr();
    s_otlr_svd2 =      m_pdf_back_svd2->Get_s_otlr();
  } else{
    tau_svd2 =         m_pdf_back_svd1->GetTau();
    mu_svd2  =         m_pdf_back_svd1->Get_mu();
    mu_delta_svd2 =    m_pdf_back_svd1->Get_mu_delta();
    f_delta_mlt_svd2 = m_pdf_back_svd1->Get_f_delta_mlt();
    f_tail_mlt_svd2 =  m_pdf_back_svd1->Get_f_tail_mlt();
    S_main_mlt_svd2 =  m_pdf_back_svd1->Get_S_main_mlt();
    S_tail_mlt_svd2 =  m_pdf_back_svd1->Get_S_tail_mlt();
    f_delta_sgl_svd2 = m_pdf_back_svd1->Get_f_delta_sgl();
    f_tail_sgl_svd2 =  m_pdf_back_svd1->Get_f_tail_sgl();
    S_main_sgl_svd2 =  m_pdf_back_svd1->Get_S_main_sgl();
    S_tail_sgl_svd2 =  m_pdf_back_svd1->Get_S_tail_sgl();
    f_otlr_svd2 =      m_pdf_back_svd1->Get_f_otlr();
    s_otlr_svd2 =      m_pdf_back_svd1->Get_s_otlr();
  }
//  upar.Add(string("tau_svd2"),         tau_svd2,         0.2, 0.1,10);  NPar++;
  upar.Add(string("mu_svd2"),          mu_svd2,          0.2,-0.2,0.2);NPar++;
  upar.Add(string("mu_delta_svd2"),    mu_delta_svd2,    0.2,-0.2,0.2);NPar++;
  upar.Add(string("f_delta_mlt_svd2"), f_delta_mlt_svd2, 0.2, 0.,1.);  NPar++;
  upar.Add(string("f_tail_mlt_svd2"),  f_tail_mlt_svd2,  0.2, 0.,1.);  NPar++;
  upar.Add(string("S_main_mlt_svd2"),  S_main_mlt_svd2,  0.2, 0.5,2.); NPar++;
  upar.Add(string("S_tail_mlt_svd2"),  S_tail_mlt_svd2,  0.2, 2. ,8.); NPar++;
  upar.Add(string("f_delta_sgl_svd2"), f_delta_sgl_svd2, 0.2, 0.,1.);  NPar++;
  upar.Add(string("f_tail_sgl_svd2"),  f_tail_sgl_svd2,  0.2, 0.,1.);  NPar++;
  upar.Add(string("S_main_sgl_svd2"),  S_main_sgl_svd2,  0.2, 0.6,2.); NPar++;
  upar.Add(string("S_tail_sgl_svd2"),  S_tail_sgl_svd2,  0.2, 2. ,8.); NPar++;
  upar.Add(string("f_otlr_svd2"),      f_otlr_svd2,      0.2, 0.,1.);  NPar++;
  upar.Add(string("s_otlr_svd2"),      s_otlr_svd2,      0.2,10.,70.); NPar++;

  upar.Add(string("f_cont"),           m_f_cont,         0.2,0.,1.);   NPar++;
  upar.Add(string("scale"),            1,                0.2,0.8,1.2); NPar++;

  cout << "Init pdf" << endl;
  pdfFcnBkg* theFCN = new pdfFcnBkg();
  MnMigrad migrad(*theFCN,upar);

//  migrad.Fix("tau_svd1");
//  migrad.Fix("mu_svd1");
//  migrad.SetValue("f_delta_mlt_svd1",1);
//  migrad.SetValue("f_delta_sgl_svd1",1);
//  migrad.Fix("f_delta_mlt_svd1");
//  migrad.Fix("f_delta_sgl_svd1");

  int NFreePar = NPar;
  if(!m_type_flag){
    NFreePar = 2;
    migrad.Fix("f_cont");

    migrad.Fix("tau");
    migrad.Fix("mu_svd1");
    migrad.Fix("mu_delta_svd1");
    migrad.Fix("f_delta_mlt_svd1");
    migrad.Fix("f_tail_mlt_svd1");
    migrad.Fix("S_main_mlt_svd1");
    migrad.Fix("S_tail_mlt_svd1");
    migrad.Fix("f_delta_sgl_svd1");
    migrad.Fix("f_tail_sgl_svd1");
    migrad.Fix("S_main_sgl_svd1");
    migrad.Fix("S_tail_sgl_svd1");
//    migrad.Fix("f_otlr_svd1");
    migrad.Fix("s_otlr_svd1");

    migrad.Fix("tau");
    migrad.Fix("mu_svd2");
    migrad.Fix("mu_delta_svd2");
    migrad.Fix("f_delta_mlt_svd2");
    migrad.Fix("f_tail_mlt_svd2");
    migrad.Fix("S_main_mlt_svd2");
    migrad.Fix("S_tail_mlt_svd2");
    migrad.Fix("f_delta_sgl_svd2");
    migrad.Fix("f_tail_sgl_svd2");
    migrad.Fix("S_main_sgl_svd2");
    migrad.Fix("S_tail_sgl_svd2");
//    migrad.Fix("f_otlr_svd2");
    migrad.Fix("s_otlr_svd2");
  } else{
    migrad.Fix("f_cont");
    migrad.Fix("scale");
//  migrad.Fix("s_otlr_svd1"); NFreePar--;
//  migrad.Fix("s_otlr_svd2"); NFreePar--;
  if(!add_otlr){
    migrad.Fix("f_otlr_svd1"); NFreePar--;
    migrad.Fix("f_otlr_svd2"); NFreePar--;
  }
  if(sgl_asc){
    migrad.Fix("f_delta_mlt_svd1"); NFreePar--;
    migrad.Fix("f_tail_mlt_svd1"); NFreePar--;
    migrad.Fix("S_main_mlt_svd1"); NFreePar--;
    migrad.Fix("S_tail_mlt_svd1"); NFreePar--;

    migrad.Fix("f_delta_mlt_svd2"); NFreePar--;
    migrad.Fix("f_tail_mlt_svd2"); NFreePar--;
    migrad.Fix("S_main_mlt_svd2"); NFreePar--;
    migrad.Fix("S_tail_mlt_svd2"); NFreePar--;
  }
  if(mlt_asc){
    migrad.Fix("f_delta_sgl_svd1"); NFreePar--;
    migrad.Fix("f_tail_sgl_svd1"); NFreePar--;
    migrad.Fix("S_main_sgl_svd1"); NFreePar--;
    migrad.Fix("S_tail_sgl_svd1"); NFreePar--;

    migrad.Fix("f_delta_sgl_svd2"); NFreePar--;
    migrad.Fix("f_tail_sgl_svd2"); NFreePar--;
    migrad.Fix("S_main_sgl_svd2"); NFreePar--;
    migrad.Fix("S_tail_sgl_svd2"); NFreePar--;
  }
  }

  cout << "Starting minimization, " << NFreePar << " free parameters" << endl;
  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  cout << "after migrad" << endl;

  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  if(m_svd != 2 && m_type_flag) WriteBkgParameters(pstate,1);
  if(m_svd != 1 && m_type_flag) WriteBkgParameters(pstate,2);

  SetBkgPdfParams(upar.Params());
  DrawSideband(NFreePar);

  make_sb_test_fit();
  return 0;
}

//int AB_fit(double vals[][2], double errs[][2],const int SetNum = 0, const bool plot_flag = true){
//  // AB fit from here
//  cout << "AB fit";
//  if(SetNum) cout << ", SetNum " << SetNum;
//  cout << endl;

//  m_cos2beta = TMath::Cos(2.*Beta(Mode(m_mode),h0Mode(m_mode))/180.*TMath::Pi());

////  vector<double> vals_vec[16], errs_vec[16];
//  double Aref[16],Bref[16];
//  for(int i=0; i<16; i++){
//    m_fitbin = bin(i);
//    Aref[i] = A(m_fitflv,m_fitbin);
//    Bref[i] = B(m_fitflv,m_fitbin);
//    m_A = Aref[i]; m_B = Bref[i];
//    cout << "A = " << m_A << ", B = " << m_B << endl;
//    GetGoodWWTTree(m_tree,m_ww_tree);
//    MnUserParameterState pstate;
//    if(make_single_fit(pstate,m_fitflv,m_fitbin,SetNum)){
//      vals[i][0] = -9;
//      errs[i][0] = -9;
//      vals[i][1] = -9;
//      errs[i][1] = -9;
//    } else{
//      vals[i][0] = pstate.Params().at(2);
//      errs[i][0] = pstate.Errors().at(2);
//      vals[i][1] = pstate.Params().at(3);
//      errs[i][1] = pstate.Errors().at(3);
//    }
//  }

////  if(!SetNum){
//  for(int i=0; i<16; i++){
//    cout << "Fit procedure " << i+1 << ":" << endl;
//    cout << "  Aref = " << Aref[i] << ", Bref = " << Bref[i] << endl;
//    for(int j=0; j<2; j++) cout << "  " << vals[i][j] << " +- " << errs[i][j] << endl;
//  }
////  }

//  if(plot_flag) draw_fit_results(vals,errs,Aref,Bref);
//}

//void AB_sampling_fit(void){
//  vector<double> val_vec[16][2];
//  vector<double> err_vec[16][2];
//  const int NSigSets = (m_NSigTot / m_NSig) < m_nsets ? (m_NSigTot / m_NSig) : m_nsets;
//  const int NBkgSets =  m_NBkgTot / m_NBkg;
//  cout <<   "NSigSets: " << NSigSets;
//  cout << ", NBkgSets: " << NBkgSets << endl;
//  bool draw_flag = false;
//  double vals[16][2],errs[16][2];
//  for(int i=1; i<=NSigSets; i++){
////    if(i == NSigSets) draw_flag = true;
//    AB_fit(vals,errs,i,draw_flag);
//    for(int j=0; j<16; j++){
//      for(int k=0; k<2; k++){
//        val_vec[j][k].push_back(vals[j][k]);
//        err_vec[j][k].push_back(errs[j][k]);
//      }
//    }
//  }
//  double Aref[16], Bref[16];
//  for(int i=0; i<16; i++){
//    int m_bin = bin(i);
//    Aref[i] = A(m_fitflv,m_bin);
//    Bref[i] = B(m_fitflv,m_bin);
//  }
//  DrawABSampling(val_vec,err_vec,Aref,Bref);
//  return;
//}

void pseudo_toy_lifetime(const int Nsmpl){
  vector<double> val_vec;
  vector<double> err_vec;
  MnUserParameterState pstate;
  for(int i=0; i<Nsmpl; i++){
    lifetime_wide_window(pstate,i);
    val_vec.push_back(pstate.Params().at(0));
    err_vec.push_back(pstate.Errors().at(0));
  }
  DrawToyLifetime(val_vec,err_vec);
  return;
}

void pseudo_toy_cpv(const int Nsmpl){
  vector<double> val_vec[2];
  vector<double> err_vec[2];
  MnUserParameterState pstate;
  for(int i=0; i<Nsmpl; i++){
    cpv_wide_window(pstate,i);
    val_vec[0].push_back(pstate.Params().at(0));
    err_vec[0].push_back(pstate.Errors().at(0));
    val_vec[1].push_back(pstate.Params().at(1));
    err_vec[1].push_back(pstate.Errors().at(1));
  }
  DrawToyCPV(val_vec,err_vec);
  return;
}

//void full_sampling_fit(void){
//  cout << "full_sampling_fit" << endl;
//  vector<double> val_vec[3];
//  vector<double> err_vec[3];
//  const int NSigSets = (m_NSigTot / m_NSig) < m_nsets ? (m_NSigTot / m_NSig) : m_nsets;
//  const int NBkgSets = m_NBkgTot / m_NBkg;
//  const int NGoodSets = m_NGoodTot / (m_NSig + m_NBkg) - 1;
//  cout <<   "NSigSets: " << NSigSets;
//  cout << ", NBkgSets: " << NBkgSets << endl;
//  cout << " NGoodSets: " << NGoodSets << endl;
//  bool draw_flag = false;
//  MnUserParameterState pstate;
//  const int NSets = m_ebeb ? NGoodSets : NSigSets;
//  for(int i=1; i<=NSets; i++){
////    if(i == NSigSets) draw_flag = true;
//    m_ebeb ? AnalyseDataSet2(i) : AnalyseDataSet1(i);// Setting fbkg(bin)
//    if(make_full_fit(pstate,i,draw_flag)) continue;
//    val_vec[0].push_back(pstate.Params().at(0));//tau
//    err_vec[0].push_back(pstate.Errors().at(0));
//    val_vec[1].push_back(pstate.Params().at(2));//sin
//    err_vec[1].push_back(pstate.Errors().at(2));
//    val_vec[2].push_back(pstate.Params().at(3));//cos
//    err_vec[2].push_back(pstate.Errors().at(3));
//  }
//  DrawFullSampling(val_vec,err_vec);
//  return;
//}

void LineTest(void){
  vector<string> angle_vec;
  angle_vec.push_back(string("0"));
  angle_vec.push_back(string("15"));
  angle_vec.push_back(string("22_5"));
  angle_vec.push_back(string("30"));
  angle_vec.push_back(string("45"));
  angle_vec.push_back(string("60"));
  angle_vec.push_back(string("67_5"));
  angle_vec.push_back(string("75"));
  angle_vec.push_back(string("90"));
  angle_vec.push_back(string("105"));
  angle_vec.push_back(string("112_5"));
  angle_vec.push_back(string("120"));
  angle_vec.push_back(string("135"));
  angle_vec.push_back(string("150"));
  angle_vec.push_back(string("157_5"));
  angle_vec.push_back(string("165"));
  const double deg_to_rad = 2.*TMath::Pi()/180.;
  const int VSize = angle_vec.size();
  double zero_vec[VSize];
  double sin_val_vec[VSize] = {TMath::Sin(0),TMath::Sin(15*deg_to_rad),TMath::Sin(22.5*deg_to_rad),TMath::Sin(30*deg_to_rad),
                               TMath::Sin(45*deg_to_rad),TMath::Sin(60*deg_to_rad),TMath::Sin(67.5*deg_to_rad),
                               TMath::Sin(75*deg_to_rad),TMath::Sin(90*deg_to_rad),TMath::Sin(105*deg_to_rad),
                               TMath::Sin(112.5*deg_to_rad),TMath::Sin(120*deg_to_rad),TMath::Sin(135*deg_to_rad),
                               TMath::Sin(150*deg_to_rad),TMath::Sin(157.5*deg_to_rad),TMath::Sin(165*deg_to_rad)};
  double cos_val_vec[VSize] = {TMath::Cos(0),TMath::Cos(15*deg_to_rad),TMath::Cos(22.5*deg_to_rad),TMath::Cos(30*deg_to_rad),
                               TMath::Cos(45*deg_to_rad),TMath::Cos(60*deg_to_rad),TMath::Cos(67.5*deg_to_rad),
                               TMath::Cos(75*deg_to_rad),TMath::Cos(90*deg_to_rad),TMath::Cos(105*deg_to_rad),
                               TMath::Cos(112.5*deg_to_rad),TMath::Cos(120*deg_to_rad),TMath::Cos(135*deg_to_rad),
                               TMath::Cos(150*deg_to_rad),TMath::Cos(157.5*deg_to_rad),TMath::Cos(165*deg_to_rad)};
  double sin_fit_vec[VSize];
  double cos_fit_vec[VSize];
  double sin_err_vec[VSize];
  double cos_err_vec[VSize];

  for(int i=0; i<VSize; i++){
    cout << "sin "  << sin_val_vec[i];
    cout << ", cos " << cos_val_vec[i];
    cout << endl;
  }

  for(int i=0; i<VSize; i++){
    zero_vec[i] = 0;
    const string fname = SigLineFile(m_mode,angle_vec[i]);
    TChain* chain = new TChain("TEvent");
    chain->Add(fname.c_str());
    SetBranchAddresses(chain);
    GetGoodWWTTree(chain,m_ww_tree);
    MnUserParameterState pstate;
    cpv_wide_window(pstate);
    sin_fit_vec[i] = pstate.Value(0);
    sin_err_vec[i] = pstate.Error(0);
    cos_fit_vec[i] = pstate.Value(1);
    cos_err_vec[i] = pstate.Error(1);
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,400);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.99);
  pad1->Draw(); pad2->Draw();

  gStyle->SetOptFit(111);
  TGraphErrors* gr_sin = new TGraphErrors(VSize,sin_val_vec,sin_fit_vec,zero_vec,sin_err_vec);
  if(m_mode == 4) gr_sin->SetTitle("Sin linearity (#omega)");
  else            gr_sin->SetTitle("Sin linearity (#pi^{0})");
  gr_sin->SetMarkerStyle(20);
  gr_sin->SetMarkerSize(1.2);
  gr_sin->SetMarkerColor(kBlue);
  gr_sin->GetXaxis()->SetTitle("Generated");
  gr_sin->GetXaxis()->SetTitleSize(0.06);
  gr_sin->GetXaxis()->SetTitleOffset(0.7);
  gr_sin->GetXaxis()->SetLabelSize(0.05);

  gr_sin->GetYaxis()->SetTitle("Fitted");
  gr_sin->GetYaxis()->SetTitleSize(0.06);
  gr_sin->GetYaxis()->SetTitleOffset(0.7);
  gr_sin->GetYaxis()->SetLabelSize(0.05);
  TGraphErrors* gr_cos = new TGraphErrors(VSize,cos_val_vec,cos_fit_vec,zero_vec,cos_err_vec);
  if(m_mode == 4) gr_cos->SetTitle("Cos linearity (#omega)");
  else            gr_cos->SetTitle("Cos linearity (#pi^{0})");
  gr_cos->SetMarkerStyle(20);
  gr_cos->SetMarkerSize(1.2);
  gr_cos->SetMarkerColor(kBlue);
  gr_cos->GetXaxis()->SetTitle("Generated");
  gr_cos->GetXaxis()->SetTitleSize(0.06);
  gr_cos->GetXaxis()->SetTitleOffset(0.7);
  gr_cos->GetXaxis()->SetLabelSize(0.05);

  gr_cos->GetYaxis()->SetTitle("Fitted");
  gr_cos->GetYaxis()->SetTitleSize(0.06);
  gr_cos->GetYaxis()->SetTitleOffset(0.7);
  gr_cos->GetYaxis()->SetLabelSize(0.05);

  pad1->cd();
  gr_sin->Draw("ap");
  gr_sin->Fit("pol1");

  pad2->cd();
  gr_cos->Draw("ap");
  gr_cos->Fit("pol1");

  if(m_mode == 4){
    c1->Print("linearity_omega.eps");
    system("linearity_omega.eps &");
    c1->Print("linearity_omega.root");
  } else{
    c1->Print("linearity_pi0.eps");
    system("linearity_pi0.eps &");
    c1->Print("linearity_pi0.root");
  }

  for(int i=0; i<VSize; i++){
    cout << "sin "  << sin_val_vec[i] << " -> " << sin_fit_vec[i] << " +- " << sin_err_vec[i];
    cout << ", cos " << cos_val_vec[i] << " -> " << cos_fit_vec[i] << " +- " << cos_err_vec[i];
    cout << endl;
  }
  return;
}

int main(const int argn, const char** argv){
  read_com_line_params(argn,argv);
  const int nstr = m_ns;
  const int cstr = m_cs;

  m_tree = new TChain("TEvent");

  if(m_toyfit){
    m_tree_pi0 = new TChain("TEvent");
    m_tree_pi0->Add(SigMCFile(1).c_str());
    cout << m_tree_pi0->GetEntries() << endl;
//    GetGoodWWTTree(m_tree_pi0,tree_ww_pi0);
    tree_ww_pi0 = m_tree_pi0->CopyTree("sz_sig>0.0001 && sz_asc>0.0001 && sigarea");
    const int pi0_evt = tree_ww_pi0->GetEntries();
    const int toy_samples_pi0 = tree_ww_pi0->GetEntries()/toy_composition[0];
    NToySamples = toy_samples_pi0;
    cout << toy_samples_pi0 << " samples available for pi0" << endl;
    SetBranchAddresses(tree_ww_pi0);
    GetShuffledVector(pi0_evt,index_pi0);

    m_tree_etagg = new TChain("TEvent");
    m_tree_etagg->Add(SigMCFile(2).c_str());
    cout << m_tree_etagg->GetEntries() << endl;
//    GetGoodWWTTree(m_tree_etagg,tree_ww_etagg);
    tree_ww_etagg = m_tree_etagg->CopyTree("sz_sig>0.0001 && sz_asc>0.0001 && sigarea");
    const int etagg_evt = tree_ww_etagg->GetEntries();
    const int toy_samples_etagg = tree_ww_etagg->GetEntries()/toy_composition[1];
    NToySamples = NToySamples < toy_samples_etagg ? NToySamples : toy_samples_etagg;
    cout << toy_samples_etagg << " samples available for eta->gg" << endl;
    SetBranchAddresses(tree_ww_etagg);
    GetShuffledVector(etagg_evt,index_etagg);

    m_tree_etappp = new TChain("TEvent");
    m_tree_etappp->Add(SigMCFile(3).c_str());
//    GetGoodWWTTree(m_tree_etappp,tree_ww_etappp);
    tree_ww_etappp = m_tree_etappp->CopyTree("sz_sig>0.0001 && sz_asc>0.0001 && sigarea");
    const int etappp_evt = tree_ww_etappp->GetEntries();
    const int toy_samples_etappp = tree_ww_etappp->GetEntries()/toy_composition[2];
    NToySamples = NToySamples < toy_samples_etappp ? NToySamples : toy_samples_etappp;
    cout << toy_samples_etappp << " samples available for eta->pi+pi-pi0" << endl;
    SetBranchAddresses(tree_ww_etappp);
    GetShuffledVector(etappp_evt,index_etappp);

    m_tree_omega = new TChain("TEvent");
    m_tree_omega->Add(SigMCFile(4).c_str());
//    GetGoodWWTTree(m_tree_omega,tree_ww_omega);
    tree_ww_omega = m_tree_omega->CopyTree("sz_sig>0.0001 && sz_asc>0.0001 && sigarea");
    const int omega_evt = tree_ww_omega->GetEntries();
    const int toy_samples_omega = tree_ww_omega->GetEntries()/toy_composition[3];
    NToySamples = NToySamples < toy_samples_omega ? NToySamples : toy_samples_omega;
    cout << toy_samples_omega << " samples available for omega" << endl;
    SetBranchAddresses(tree_ww_omega);
    GetShuffledVector(omega_evt,index_omega);

    init(4);
    cout << NToySamples << " samples accepted" << endl;
    if(no_interf) pseudo_toy_lifetime(NToySamples);
    else          pseudo_toy_cpv(NToySamples);

    return 0;
  } else if(m_calc_CS){
    cout << "calc_CS" << endl;
//    m_tree->Add(SigMCFile(1).c_str());
//    m_tree->Add(SigMCFile(2).c_str());
//    m_tree->Add(SigMCFile(3).c_str());
//    m_tree->Add(SigMCFile(4).c_str());
    m_tree->Add(SigMCFile(5).c_str());
//    m_tree->Add(SigMCFile(10).c_str());
//    m_tree->Add(SigMCFile(20).c_str());
  } else if(m_nuisance){
    cout << "GenMC!" << endl;
//     m_tree->Add(GenWWFile(m_mode,0,nstr,cstr).c_str());
    m_tree->Add(GenWWFile(1,0,nstr,cstr).c_str());
    m_tree->Add(GenWWFile(2,0,nstr,cstr).c_str());
    m_tree->Add(GenWWFile(3,0,nstr,cstr).c_str());
    m_tree->Add(GenWWFile(4,0,nstr,cstr).c_str());
    //    m_tree->Add(GenWWFile(5,0,nstr,cstr).c_str());
//        m_tree->Add(SigMCFile(1).c_str());
//    m_tree->Add(SigLineFile(1,string("22_5")).c_str());
//    m_tree->Add(SigLineFile(4,string("22_5")).c_str());
  } else if(m_data){
    cout << "Data!" << endl;
    m_tree->Add(DataFile(m_mode).c_str());
  } else if(sigmc){
    cout << "SigMC!" << endl;
    if(m_genfit){
      m_tree->Add(SigMCFile(1).c_str());
      m_tree->Add(SigMCFile(4).c_str());
    } else{
      m_tree->Add(SigMCFile(m_mode).c_str());
    }
  } else{
    if(m_gg){
      cout << "GG!" << endl;
      m_tree->Add(GenWWFile(1,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(2,0,nstr,cstr).c_str());
      if(no_interf){
        m_tree->Add(GenWWFile(10,0,nstr,cstr).c_str());
        m_tree->Add(GenWWFile(20,0,nstr,cstr).c_str());
      }
    }
    if(m_ppp){
        cout << "PPP!" << endl;
      m_tree->Add(GenWWFile(3,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(4,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(5,0,nstr,cstr).c_str());
    }
    if(!m_gg && !m_ppp){
      cout << "GenWW!" << endl;
//      m_tree->Add(GenWWFile(m_mode,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(1,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(2,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(3,0,nstr,cstr).c_str());
      m_tree->Add(GenWWFile(4,0,nstr,cstr).c_str());
//      m_tree->Add(GenWWFile(1,0,nstr,cstr).c_str());
//      m_tree->Add(GenWWFile(2,0,nstr,cstr).c_str());
//      m_tree->Add(GenWWFile(3,0,nstr,cstr).c_str());
//      m_tree->Add(GenWWFile(4,0,nstr,cstr).c_str());
//      m_tree->Add(GenWWFile(5,0,nstr,cstr).c_str());
    }
  }

  init(m_mode);
  cout << "Initialization complete" << endl;

  if(m_line_test){
    LineTest();
    return 0;
  }

  if(make_bins_scan){
    if(!no_interf) bins_scan_wide_window();
    else           bins_scan_wide_lifetime();
    return 0;
  }

//  if(ABfit){
//    double vals[16][2],errs[16][2];
//    AB_fit(vals,errs);
//    return 0;
//  }

  GetGoodWWTTree(m_tree,m_ww_tree);
  cout << "Good tree is ready" << endl;
  cout << m_tree->GetEntries() << " events total" << endl;
  cout << m_ww_tree->GetEntries() << " good events" << endl;
  if(!no_bkg) cout << m_ww_sb_tree->GetEntries() << " events in sideband" << endl;

  if(m_genfit){
    MnUserParameterState pstate;
    fit_gen(pstate);
    return 0;
  }

  if(m_calc_K){
    calc_K();
    return 0;
  }

  if(m_calc_CS){
    calc_CS();
    return 0;
  }

//  if(m_toysize){
//    if(no_interf) pseudo_toy_lifetime();
//    else          pseudo_toy_cpv();
//    return 0;
//  }

  if(!no_bkg){
    MnUserParameterState sbpstate;
    sideband_wide_window(sbpstate);
    m_scale1 = sbpstate.Value(0);
    m_scale2 = sbpstate.Value(1);
    m_scale1_err = sbpstate.Error(0);
    m_scale2_err = sbpstate.Error(1);
  }
  if(no_interf){// lifetime fit
    MnUserParameterState pstate;
    lifetime_wide_window(pstate);
  } else{// CPV fit
    MnUserParameterState pstate;
    if(m_nuisance) nuisance_cpv_fit(pstate);
    else cpv_wide_window(pstate);
  }
  return 0;

//  GetGoodTTree(m_tree,m_fitbin,m_fitflv,sideband_fit);

//  if(m_norm_test_flag){
//    norm_test(m_Neve);
//    return 0;
//  }

//  if(sideband_fit){
//                cout << "NEntries: Sideband:      " << m_good_sb_tree->GetEntries() << endl;
//    if(!m_data) cout << "NEntries: Signal region: " << m_test_tree->GetEntries() << endl;
////    if(m_mode) cout << ", sig: " << m_sig_tree->GetEntries() << ", bkg: " <<  m_bkg_tree->GetEntries() << endl;
////    else cout << ", good: " << m_good_tree->GetEntries() << endl;
//    make_sideband_fit();
//    return 0;
//  }

////  if(m_mode && !full_ds_fit){
////    if(FullFit){ full_sampling_fit(); return 0;}
//////    else{ AB_sampling_fit(); return 0;}
////  }

//  if(FullFit){
//    MnUserParameterState pstate;
//    make_full_fit(pstate);
//    return 0;
//  } else{
//    double vals[16][2],errs[16][2];
//    AB_fit(vals,errs);
//    return 0;
//  }

//  return 0;
}
