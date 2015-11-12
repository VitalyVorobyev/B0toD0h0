#include "PseudoToyFit.h"

MnUserParameterState make_single_fit(){
  int NPar = 0;
  MnUserParameters upar;
  upar.Add(string("btau"),m_btau,0.2*m_btau,0.,10*m_btau); NPar++;
  upar.Add(string("dm"),  m_dm,  0.1*m_dm,  0.,5.*m_dm);   NPar++;
  upar.Add(string("A"),   m_A,   0.1,      -3.,3.);        NPar++;
  upar.Add(string("B"),   m_B,   0.1,      -3.,3.);        NPar++;

  pdfFcnSingle* theFCN = new pdfFcnSingle();
  MnMigrad migrad(*theFCN,upar);
  migrad.Fix("btau");
  migrad.Fix("dm");
//  migrad.Fix("A");
//  migrad.Fix("B");
  FunctionMinimum min = migrad();
  MnUserParameterState pstate = min.UserState();
  for(int i=0; i<NPar; i++){
    cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
    upar.SetValue(i,pstate.Value(i));
    cout << pstate.Value(i) << " +- " << pstate.Error(i) << endl;
  }

  const double& btau = upar.Value(0);
  const double& dm   = upar.Value(1);
  const double& A    = upar.Value(2);
  const double& B    = upar.Value(3);

  TH1I dh("dh","dh",NBins,dtmin,dtmax);
  const int NTot = m_good_tree->GetEntries();
  for(int i=0; i<NTot; i++){
    GetEvent(i);
    dh.Fill(m_dt);
  }
  double pdf_arr[NDots],dt_arr[NDots];
  double norm = 0;

  const double ddt = (dtmax-dtmin)/(double)NDots;
  const double Nddt = ddt*NDots/NBins;

  double pdf = 0;
  cout << "Set parameters: " << btau << ", " << dm << ", " << A << ", " << B << endl;
  m_pdf->SetTauDm(btau,dm);
  m_pdf->SetAB(A,B);
  for(int i=0; i<NDots; i++){
    pdf_arr[i] = 0;
    dt_arr[i] = dtmin+(i+0.5)*ddt;
    const double dt = dt_arr[i];//cm2ps;
    for(int j=0; j<NTot; j++){
      GetEvent(j);
      m_pdf->SetAkCk(m_costhBcms,0.5*10.58);
      pdf = m_pdf->PdfAB(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
      pdf_arr[i] += pdf;
    }
    pdf_arr[i] *= Nddt;//norm*NTot;
    norm       += pdf_arr[i];
  }
  for(int i=0; i<NDots; i++) pdf_arr[i] *= NTot/norm;
  norm *= Nddt;

  cout << "norm = " << norm << endl;

  TGraph* gr = new TGraph(NDots,dt_arr,pdf_arr);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.);
  gr->SetMarkerColor(kRed);

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
//  c1->SetLogy();
  dh.Draw("e");
  gr->Draw("same");

  if(draw_plots){
    c1->Print("fit.png");
    c1->Print("fit.root");
    system("display fit.png &");
  }

  return pstate;
}

MnUserParameterState make_full_fit(void){
    stringstream out;
    int NPar = 0;
    MnUserParameters upar;
    upar.Add(string("btau"),    m_btau,    0.2*m_btau,0, 10*m_btau); NPar++;
    upar.Add(string("dm"),      m_dm,      0.1*m_dm,  0.,5.*m_dm);   NPar++;
    upar.Add(string("sin2beta"),m_sin2beta,0.1,      -3.,3.);        NPar++;
    upar.Add(string("cos2beta"),m_cos2beta,0.1,      -3.,3.);        NPar++;

    pdfFcn* theFCN = new pdfFcn();
    MnMigrad migrad(*theFCN,upar);

    if(fix_btau)                  migrad.Fix("btau");
    if(fix_dm       || no_interf) migrad.Fix("dm");
    if(fix_sin2beta || no_interf) migrad.Fix("sin2beta");
    if(fix_cos2beta || no_interf) migrad.Fix("cos2beta");

    cout << "Starting minimization" << endl;
    FunctionMinimum min = migrad();
    MnUserParameterState pstate = min.UserState();
    cout << "after migrad" << endl;

    for(int i=0; i<NPar; i++){
      cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
      upar.SetValue(i,pstate.Value(i));
      cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
    }

    return pstate;

    SetPDFParams(upar.Params());
    const int NTot = m_good_tree->GetEntries();
    const double ddt = (dtmax-dtmin)/(double)NDots;
    double dt_arr[NDots];

    cout << "Init hists and fill graphs" << endl;
    for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;
    if(no_interf){
      double Norm = 0;
      double pdf_array[NDots];
      for(int i=0; i<NDots; i++){
        pdf_array[i] = 0;
        const double& dt = dt_arr[i];
        for(int j=SampleNumber*EventsInASample; j<(SampleNumber+1)*EventsInASample; j++){
          GetEvent(j);
          m_pdf->SetAkCk(m_costhBcms,0.5*10.58);
          if(!no_np){
            pdf_array[i] += m_pdf->Pdf(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,false);
          } else{
            pdf_array[i] += m_pdf->NoNPPdf(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc,true,false);
          }
          Norm += pdf_array[i];
        }
        pdf_array[i] *= ddt*NDots/(double)NBins;
      }
      Norm *= ddt/EventsInASample;
      cout << "Norm = " << Norm << "Nevents = " << EventsInASample << endl;
      TH1I* DH = new TH1I("DH","DH",NBins,dtmin,dtmax);
      for(int i=SampleNumber*EventsInASample; i<(SampleNumber+1)*EventsInASample; i++){
        GetEvent(i);
        DH->Fill(m_dt);
      }
      TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
      GR->SetMarkerStyle(kDot);
      GR->SetMarkerSize(1.5);
      GR->SetMarkerColor(kBlue);

      TCanvas* c1 = new TCanvas("c1","c1",800,600);
      c1->cd();
      DH->Draw("e");
      GR->Draw("same");
      c1->Print("fit_noint.root");
      c1->Print("fit_noint.png");

      system("display fit_noint.png &");
      return pstate;
    }
    else{
    int Nev[2][16];
    TH1I* dh[2][16];
    double pdf_arr[2][16][NDots];
    double norm[2][16];
    double norm1[2][16];
    for(int k=0; k<2; k++){
      cout << "flv = " << flv(k) << endl;
      m_pdf->SetFlvXi(flv(k),xi);
      for(int i=0; i<16; i++){
        Nev[k][i]    = 0;
        norm[k][i]   = 0;
        norm1[k][i]  = 0;
        out.str("");
        out << "dh" << bin(i) << "_" << k;
        dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dtmin,dtmax);
        m_bin = bin(i);
        cout << "  bin = " << m_bin << endl;
        m_pdf->SetKKCS(K(m_bin),K(-m_bin),C(m_bin),S(m_bin));
        for(int j=0; j<NDots; j++){
          const double& dt = dt_arr[j];
//          cout << "    dt = " << dt << endl;
          pdf_arr[k][i][j] = 0;
          for(int l=SampleNumber*EventsInASample; l<(SampleNumber+1)*EventsInASample; l++){
            GetEvent(l);
            m_pdf->SetAkCk(m_costhBcms,0.5*10.58);
            if(!no_np){
              pdf_arr[k][i][j] += m_pdf->Pdf(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
            } else{
              pdf_arr[k][i][j] += m_pdf->NoNPPdf(dt,ntrk_rec,sz_rec,chisq_rec,ndf_rec,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
            }
            norm[k][i] += pdf_arr[k][i][j];
          }
          pdf_arr[k][i][j] *= ddt/EventsInASample;
        }
        norm[k][i] *= ddt/EventsInASample;
      }
    }

    int NEveCounter = 0;
    for(int i=SampleNumber*EventsInASample; i<(SampleNumber+1)*EventsInASample; i++){
      GetEvent(i);
      const int k = flv_ind(m_flv);
      const int j = bin_ind(m_bin);
      dh[k][j]->Fill(m_dt);
      Nev[k][j]++;
      NEveCounter++;
    }
    cout << "NEveCounter = " << NEveCounter << endl;

    cout << "Set norm" << endl;
    for(int k=0; k<2; k++){
      for(int i=0; i<16; i++){
        for(int j=0; j<NDots; j++){
          pdf_arr[k][i][j]  *= Nev[k][i]*NDots/(double)NBins;
          norm1[k][i]       += pdf_arr[k][i][j];
        }
      }
    }

    TGraph* gr[2][16];
    for(int k=0; k<2; k++){
      for(int i=0; i<16; i++){
        gr[k][i]  = new TGraph(NDots,dt_arr,pdf_arr[k][i]);
        gr[k][i]->SetMarkerStyle(kDot);
        gr[k][i]->SetMarkerSize(1.5);
        k == 0 ? gr[k][i]->SetMarkerColor(kRed) : gr[k][i]->SetMarkerColor(kRed);
      }
    }

    cout << "Draw" << endl;
    TCanvas* c1 = new TCanvas("c1","c1",2600,800);
    TPad *pad[8];
    TPad *padb[8];
    for(int i=0; i<8; i++){
      out.str("");
      out << "pad" << i;
      pad[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.99);
      pad[i]->Draw();
      out.str("");
      out << "padb" << i;
      padb[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.49);
      padb[i]->Draw();
    }
    int k = 0;
    for(int i=0; i<8; i++){
      pad[i]->cd();
      dh[k][i]->Draw("e");
      gr[k][i]->Draw("same");

      padb[i]->cd();
      dh[k][i+8]->Draw("e");
      gr[k][i+8]->Draw("same");
    }
    c1->Print("gen_fit.root");
    c1->Print("gen_fit.png");

    system("display gen_fit.png &");

    TCanvas* c2 = new TCanvas("c2","c2",2600,800);
    c2->cd();
    TPad *pad1[8];
    TPad *padb1[8];
    for(int i=0; i<8; i++){
      out.str("");
      out << "pad1" << i;
      pad1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.51,0.125*(i+1)-0.01,0.99);
      pad1[i]->Draw();
      out.str("");
      out << "padb1" << i;
      padb1[i] = new TPad(out.str().c_str(),out.str().c_str(),0.125*i+0.01,0.01,0.125*(i+1)-0.01,0.49);
      padb1[i]->Draw();
    }
    k = 1;
    for(int i=0; i<8; i++){
      pad1[i]->cd();
      dh[k][i]->Draw("e");
      gr[k][i]->Draw("same");

      padb1[i]->cd();
      dh[k][i+8]->Draw("e");
      gr[k][i+8]->Draw("same");
    }
    out.str("");
    out << "ptoy_fit" << SampleNumber << ".root";
    c2->Print(out.str().c_str());
    out.str("");
    out << "ptoy_fit" << SampleNumber << ".png";
    c2->Print(out.str().c_str());

    out.str("");
    out << "display ptoy_fit" << SampleNumber << ".png &";
    system(out.str().c_str());
    cout << "NEveCounter = " << NEveCounter << endl;
    }
    return pstate;
}

int main(const int argn, const char** argv){
  TFile *ifile;
  int _mode = 4;
  int _binfit = 1;
  int _flv = 1;
  if(argn == 2) sscanf(argv[1],"%d",&_mode);
  if(argn == 3){
    sscanf(argv[1],"%d",&_mode);
    sscanf(argv[2],"%d",&_binfit);
    if(_binfit) return -2;
  }
  if(argn == 4){
    sscanf(argv[1],"%d",&_mode);
    sscanf(argv[2],"%d",&_binfit);
    sscanf(argv[3],"%d",&_flv);
    if(abs(_flv) != 1) return -3;
  }
  switch(_mode){
  case 1:
    cout << "Mode: pi0" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigPi0_s5_full.root");
    xi = 1;
    break;
  case 2:
    cout << "Mode: eta -> gg" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
    xi = 1;
    break;
  case 3:
    cout << "Mode: eta -> pi+pi-pi0" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigEta_full.root");
    xi = 1;
    break;
  case 4:
    cout << "Mode: omega" << endl;
    if(no_interf){
      ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s4_full.root");
    } else{
      ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s4_full.root");
    }
    xi = -1;
    break;
  case 5:
    cout << "Mode: omega (rho Ks0)" << endl;
    ifile  = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s3_full.root");
    xi = -1;
    break;
  default:
    cout << "Wrong mode" << endl;
    return -1;
  }
  m_pdf = new RkRdetRnpPdf(1,m_svd);
  init(_mode);
  m_tree = (TTree*)ifile->Get("TEvent");
  m_tree->SetBranchAddress("exp",&m_exp);
  m_tree->SetBranchAddress("dz", &m_dt);
  m_tree->SetBranchAddress("costhBcms", &m_costhBcms);
  m_tree->SetBranchAddress("bin_mc",&m_bin);
  m_tree->SetBranchAddress("flv_mc",&m_flv);
  m_tree->SetBranchAddress("ntrk_sig",&ntrk_rec);
  m_tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  m_tree->SetBranchAddress("ndf_z_sig",&ndf_rec);
  m_tree->SetBranchAddress("ndf_z_asc",&ndf_asc);
  m_tree->SetBranchAddress("sz_sig",&sz_rec);
  m_tree->SetBranchAddress("sz_asc",&sz_asc);
  m_tree->SetBranchAddress("chisq_z_sig",&chisq_rec);
  m_tree->SetBranchAddress("chisq_z_asc",&chisq_asc);
  m_good_tree = GetGoodTTree(m_tree,_mode);

  const int NTot = m_good_tree->GetEntries();
  NSamples = NTot / EventsInASample;
  cout << NSamples << " to be considered" << endl;
  int j,k;

  vector< vector<double> > vals_vec, errs_vec;
  vector<double> Aref,Bref;

  if(argn == 4){
    if(abs(_binfit)>8) return -1;
    if(_binfit){
      m_good_tree = GetGoodTTree(m_tree,_mode,_binfit,_flv);
//      m_good_tree->Print();
      m_flv = _flv; m_bin = _binfit;
      Aref.push_back(A(m_flv,m_bin));
      Bref.push_back(B(m_flv,m_bin));
      m_A = Aref[0]; m_B = Bref[0];
      MnUserParameterState upar = make_single_fit();
      vals_vec.push_back(upar.Params());
      errs_vec.push_back(upar.Errors());
    } else{
      draw_plots = false;
      m_flv = _flv;
      for(int i=0; i<16; i++){
        m_bin = bin(i);
        Aref.push_back(A(m_flv,m_bin));
        Bref.push_back(B(m_flv,m_bin));
        m_A = Aref[i]; m_B = Bref[i];
        m_good_tree = GetGoodTTree(m_tree,_mode,m_bin,m_flv);
//        m_good_tree->Print();
        MnUserParameterState upar = make_single_fit();
        vals_vec.push_back(upar.Params());
        errs_vec.push_back(upar.Errors());
      }
    }
  } else if(argn == 3){
    draw_plots = false;
    for(k=0; k<2; k++){
      for(j=0; j<16; j++){
        m_flv = flv(k);
        m_bin = bin(j);
        Aref.push_back(A(m_flv,m_bin));
        Bref.push_back(B(m_flv,m_bin));
        m_A = Aref[16*k+j]; m_B = Bref[16*k+j];
        m_good_tree = GetGoodTTree(m_tree,_mode,m_bin,m_flv);
//        m_good_tree->Print();
        MnUserParameterState upar = make_single_fit();
        vals_vec.push_back(upar.Params());
        errs_vec.push_back(upar.Errors());
      }
    }
  } else {
    m_good_tree = GetGoodTTree(m_tree,_mode);
    for(int i=0; i<NSamples; i++){
      MnUserParameterState upar = make_full_fit();
      vals_vec.push_back(upar.Params());
      errs_vec.push_back(upar.Errors());
      SampleNumber++;
    }
    draw_fit_results2(vals_vec,errs_vec);
    return 0;
  }

  ofstream ofile;
  ofile.open("fit.txt",std::ofstream::out);
  ofile << "argn = " << argn << ", bin = " << _binfit << ", flv = " << _flv << endl;
  for(int i=0; i<vals_vec.size(); i++){
    cout << "Fit procedure " << i+1 << ":" << endl;
    cout << "  Aref = " << Aref[i] << ", Bref = " << Bref[i] << endl;
    for(int j=0; j<vals_vec[i].size(); j++){
      cout << "  " << vals_vec[i][j] << " +- " << errs_vec[i][j] << endl;
      ofile << vals_vec[i][j] << " +- " << errs_vec[i][j] << " ";
    }
    ofile << endl;
  }

  ofile.close();
  draw_fit_results(vals_vec,errs_vec,Aref,Bref);
  return 0;
}

