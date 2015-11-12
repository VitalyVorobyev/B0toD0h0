#include "BpTauFit.h"

int make_fit(void){
  int NPar = 0;
  MnUserParameters upar;
  upar.Add(string("btau"),          m_btau,0.2,1.,2.);                        NPar++;
  upar.Add(string("f_ol_sgl_svd1"), m_pdf_sig_svd1->Get_f_ol_sgl(),0.1,0.,1.);NPar++;
  upar.Add(string("f_ol_mlt_svd1"), m_pdf_sig_svd1->Get_f_ol_mlt(),0.1,0.,1.);NPar++;
  upar.Add(string("f_ol_sgl_svd2"), m_pdf_sig_svd2->Get_f_ol_sgl(),0.1,0.,1.);NPar++;
  upar.Add(string("f_ol_mlt_svd2"), m_pdf_sig_svd2->Get_f_ol_mlt(),0.1,0.,1.);NPar++;

  int NFreePar = NPar;
  pdfFcn* theFCN = new pdfFcn();
  MnMigrad migrad(*theFCN,upar);

  if(m_svd == 2){
    migrad.Fix("f_ol_sgl_svd1"); NFreePar--;
    migrad.Fix("f_ol_mlt_svd1"); NFreePar--;
  } else{
    if(single_d0 || sgl_asc){ migrad.Fix("f_ol_mlt_svd1"); NFreePar--;}
    else{                     migrad.Fix("f_ol_sgl_svd1"); NFreePar--;}
  }
  if(m_svd == 1){
    migrad.Fix("f_ol_sgl_svd2"); NFreePar--;
    migrad.Fix("f_ol_mlt_svd2"); NFreePar--;
  } else{
    if(sgl_asc){ migrad.Fix("f_ol_mlt_svd2"); NFreePar--;}
    if(mlt_asc){ migrad.Fix("f_ol_sgl_svd2"); NFreePar--;}
  }

  FunctionMinimum min = migrad();
  if(!min.IsValid()){
    cout << "Fit is not valid" << endl;
    return min.IsValid();
  }
  MnUserParameterState pstate = min.UserState();
  cout << "after migrad" << endl;

  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  SetPDFParams(upar.Params());
  Draw_NoTag(NFreePar);
  return min.IsValid();
}

int make_sideband_fit(void){
  int NPar = 0;
  MnUserParameters upar;
  upar.Add(string("tau_svd1"),         m_pdf_back_svd1->GetTau(),          0.2, 0.1  ,10);  NPar++;
  upar.Add(string("mu_svd1"),          m_pdf_back_svd1->Get_mu(),          0.2,-10.,10.); NPar++;
  upar.Add(string("mu_delta_svd1"),    m_pdf_back_svd1->Get_mu_delta(),    0.2,-10.,10.); NPar++;
  upar.Add(string("f_delta_mlt_svd1"), m_pdf_back_svd1->Get_f_delta_mlt(), 0.2, 0. ,1.);  NPar++;
  upar.Add(string("f_tail_mlt_svd1"),  m_pdf_back_svd1->Get_f_tail_mlt(),  0.2, 0. ,1.);  NPar++;
  upar.Add(string("S_main_mlt_svd1"),  m_pdf_back_svd1->Get_S_main_mlt(),  0.2, 0. ,10.); NPar++;
  upar.Add(string("S_tail_mlt_svd1"),  m_pdf_back_svd1->Get_S_tail_mlt(),  0.2, 0  ,1000.); NPar++;
  upar.Add(string("f_delta_sgl_svd1"), m_pdf_back_svd1->Get_f_delta_sgl(), 0.2, 0. ,1.);  NPar++;
  upar.Add(string("f_tail_sgl_svd1"),  m_pdf_back_svd1->Get_f_tail_sgl(),  0.2, 0. ,1.);  NPar++;
  upar.Add(string("S_main_sgl_svd1"),  m_pdf_back_svd1->Get_S_main_sgl(),  0.2, 0  ,10.); NPar++;
  upar.Add(string("S_tail_sgl_svd1"),  m_pdf_back_svd1->Get_S_tail_sgl(),  0.2, 0  ,100.); NPar++;

  upar.Add(string("tau_svd2"),         m_pdf_back_svd2->GetTau(),          0.2, 0.1  ,10);  NPar++;
  upar.Add(string("mu_svd2"),          m_pdf_back_svd2->Get_mu(),          0.2,-10.,10.); NPar++;
  upar.Add(string("mu_delta_svd2"),    m_pdf_back_svd2->Get_mu_delta(),    0.2,-10.,10.); NPar++;
  upar.Add(string("f_delta_mlt_svd2"), m_pdf_back_svd2->Get_f_delta_mlt(), 0.2, 0. ,1.);  NPar++;
  upar.Add(string("f_tail_mlt_svd2"),  m_pdf_back_svd2->Get_f_tail_mlt(),  0.2, 0. ,1.);  NPar++;
  upar.Add(string("S_main_mlt_svd2"),  m_pdf_back_svd2->Get_S_main_mlt(),  0.2, 0. ,10.); NPar++;
  upar.Add(string("S_tail_mlt_svd2"),  m_pdf_back_svd2->Get_S_tail_mlt(),  0.2, 0  ,1000.); NPar++;
  upar.Add(string("f_delta_sgl_svd2"), m_pdf_back_svd2->Get_f_delta_sgl(), 0.2, 0. ,1.);  NPar++;
  upar.Add(string("f_tail_sgl_svd2"),  m_pdf_back_svd2->Get_f_tail_sgl(),  0.2, 0. ,1.);  NPar++;
  upar.Add(string("S_main_sgl_svd2"),  m_pdf_back_svd2->Get_S_main_sgl(),  0.2, 0  ,10.); NPar++;
  upar.Add(string("S_tail_sgl_svd2"),  m_pdf_back_svd2->Get_S_tail_sgl(),  0.2, 0  ,100.); NPar++;

  pdfFcnBkg* theFCN = new pdfFcnBkg();
  MnMigrad migrad(*theFCN,upar);

  migrad.Fix("tau_svd1");
  migrad.Fix("tau_svd2");

  if(m_svd == 1){
    migrad.Fix("f_delta_mlt_svd2");
    migrad.Fix("f_tail_mlt_svd2");
    migrad.Fix("S_main_mlt_svd2");
    migrad.Fix("S_tail_mlt_svd2");
    migrad.Fix("f_delta_sgl_svd2");
    migrad.Fix("f_tail_sgl_svd2");
    migrad.Fix("S_main_sgl_svd2");
    migrad.Fix("S_tail_sgl_svd2");
  }
  if(m_svd == 2){
    migrad.Fix("f_delta_mlt_svd1");
    migrad.Fix("f_tail_mlt_svd1");
    migrad.Fix("S_main_mlt_svd1");
    migrad.Fix("S_tail_mlt_svd1");
    migrad.Fix("f_delta_sgl_svd1");
    migrad.Fix("f_tail_sgl_svd1");
    migrad.Fix("S_main_sgl_svd1");
    migrad.Fix("S_tail_sgl_svd1");
  }
  if(sgl_asc){
    migrad.Fix("f_delta_mlt_svd1");
    migrad.Fix("f_tail_mlt_svd1");
    migrad.Fix("S_main_mlt_svd1");
    migrad.Fix("S_tail_mlt_svd1");

    migrad.Fix("f_delta_mlt_svd2");
    migrad.Fix("f_tail_mlt_svd2");
    migrad.Fix("S_main_mlt_svd2");
    migrad.Fix("S_tail_mlt_svd2");
  }
  if(mlt_asc){
    migrad.Fix("f_delta_sgl_svd1");
    migrad.Fix("f_tail_sgl_svd1");
    migrad.Fix("S_main_sgl_svd1");
    migrad.Fix("S_tail_sgl_svd1");

    migrad.Fix("f_delta_sgl_svd2");
    migrad.Fix("f_tail_sgl_svd2");
    migrad.Fix("S_main_sgl_svd2");
    migrad.Fix("S_tail_sgl_svd2");
  }

  cout << "Starting minimization" << endl;
  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  cout << "after migrad" << endl;

  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  const int NTot = m_good_tree->GetEntries();

  const double ddt = (dtmax-dtmin)/(double)NDots;
  double dt_arr[NDots];

  cout << "Init hists and fill graphs" << endl;
  for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;

  SetBkgPdfParams(upar.Params());

  double Norm = 0;
  double pdf_array[NDots];
  for(int i=0; i<NDots; i++){
    pdf_array[i] = 0;
    const double& dt = dt_arr[i];
    for(int j=0; j<NTot; j++){
      GetEvent(j);
      if(m_exp>30) pdf_array[i] += m_pdf_back_svd2->Pdf(dt,sum_sigma(sz_rec,sz_asc),ndf_asc);
      else         pdf_array[i] += m_pdf_back_svd1->Pdf(dt,sum_sigma(sz_rec,sz_asc),ndf_asc);
      Norm += pdf_array[i];
    }
    pdf_array[i] *= ddt*NDots/(double)NBins;
  }
  Norm *= ddt/NTot;
  cout << "Norm = " << Norm << "Nevents = " << NTot << endl;
  TH1I* DH = new TH1I("DH","DH",NBins,dtmin,dtmax);
  for(int i=0; i<NTot; i++){
    GetEvent(i);
    DH->Fill(m_dt);
  }
  DH->SetMarkerStyle(20);
  DH->SetMarkerColor(kBlue);
  DH->SetMarkerSize(1.2);

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetMarkerColor(kBlue);
  GR->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  c1->SetGrid();
  DH->GetXaxis()->SetTitle("#Deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.05);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->GetXaxis()->SetTitleOffset(0.9);
  DH->Draw("e");
  GR->Draw("same");
  c1->Print("fit_sideband.root");
  c1->Print("fit_sideband.png");

  system("display fit_sideband.png &");

  return 0;
}

int main(const int argn, const char** argv){
  bool sideband = false;
  for(int i=1; i<argn; i++){
    if(string(argv[i]) == string("sideband")) sideband = true;
    if(string(argv[i]) == string("generic"))  generic = true;
    if(string(argv[i]) == string("signal"))   signal = true;
    if(string(argv[i]) == string("svd1"))     m_svd = 1;
    if(string(argv[i]) == string("svd2"))     m_svd = 2;
    if(string(argv[i]) == string("d0"))       single_d0 = true;
    if(string(argv[i]) == string("mlt"))      mlt_asc = true;
    if(string(argv[i]) == string("sgl"))      sgl_asc = true;
  }

  m_tree = new TChain("TEvent");
//  TFile *ifile;
  if(signal){
    m_tree->Add("/home/vitaly/B0toDh0/Bp2D0pi/FIL_bp2d0pip_sigmc.root");
  }
  else if(generic){
    m_tree->Add("/home/vitaly/B0toDh0/Bp2D0pi/FIL_b2dpi_gen_v2_0_10.root");
  } else{
    m_tree->Add("/home/vitaly/B0toDh0/Bp2D0pi/FIL_b2dpi_data_v2.root");
  }

  const bool mcflag = generic || signal;
  if(single_d0){
    m_pdf_back_svd1 = new RbkgPdf(66,1,mcflag);
    m_pdf_back_svd2 = new RbkgPdf(66,2,mcflag);
  } else{
    m_pdf_back_svd1 = new RbkgPdf(55,1,mcflag);
    m_pdf_back_svd2 = new RbkgPdf(55,2,mcflag);
  }
  m_pdf_sig_svd1 = new RkRdetRnpPdf(mcflag,1,true);
  m_pdf_sig_svd2 = new RkRdetRnpPdf(mcflag,2,true);
  init();
  m_tree->SetBranchAddress("exp",&m_exp);
  if(single_d0){
    m_tree->SetBranchAddress("dz_d0", &m_dt);
    m_tree->SetBranchAddress("sz_sig_d0",&sz_rec);
    m_tree->SetBranchAddress("good_icpv_sgl", &good_icpv_sgl);
  } else{
    m_tree->SetBranchAddress("dz", &m_dt);
    m_tree->SetBranchAddress("sz_sig",&sz_rec);
    m_tree->SetBranchAddress("good_icpv_mlt", &good_icpv_mlt);
  }
  m_tree->SetBranchAddress("ndf_z_asc",&ndf_asc);
  m_tree->SetBranchAddress("ndf_z_sig",&ndf_rec);
  m_tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  m_tree->SetBranchAddress("sz_asc",&sz_asc);
  m_tree->SetBranchAddress("chisq_z_sig",&chisq_rec);
  m_tree->SetBranchAddress("chisq_z_asc",&chisq_asc);
  m_tree->SetBranchAddress("costhBcms",&m_costhBcms);
  m_good_tree = GetGoodTTree(m_tree,sideband);

  if(sideband){
    make_sideband_fit();
  } else{
    make_fit();
  }

  return 0;
}
