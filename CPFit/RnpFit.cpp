#include "RnpFit.h"
#include "TStyle.h"

MnUserParameterState make_single_fit(void){
  MnUserParameters upar;
  return upar;
}

int make_full_fit(void){
//  stringstream out;
  int NPar = 0;
  MnUserParameters upar;
  upar.Add(string("f_ol_sgl"),   m_pdf->Get_f_ol_sgl(),0.1,0.,1.);       NPar++;
  upar.Add(string("f_ol_mlt"),   m_pdf->Get_f_ol_mlt(),0.1,0.,1.);       NPar++;
  upar.Add(string("f_sigma_ol"), m_pdf->Get_sigma_ol(),0.1,20.,150.);  NPar++;

  // * multiple Rdet * //
  upar.Add(string("Sasc0"), m_pdf->Get_Sasc(0),0.1,0.,10.);  NPar++;
  upar.Add(string("Sasc1"), m_pdf->Get_Sasc(1),0.1,0.,10.);  NPar++;

  // * single Rdet * //
  upar.Add(string("Smn_asc"), m_pdf->Get_Smn_asc(),0.1,0.,10.);  NPar++;
  upar.Add(string("Stl_asc"), m_pdf->Get_Stl_asc(),0.1,0.,10.);  NPar++;
  upar.Add(string("ftl_asc"), m_pdf->Get_ftl_asc(),0.1,0.,1.);   NPar++;

  // * single Rnp * //
  upar.Add(string("fd_np_sgl0"),m_pdf->Get_fd_np_sgl(0),0.1,-1.,1.); NPar++;
  upar.Add(string("fd_np_sgl1"),m_pdf->Get_fd_np_sgl(1),0.1,-1.,1.); NPar++;
  upar.Add(string("fp_np_sgl"),m_pdf->Get_fp_np_sgl(),0.1,-1.,1.); NPar++;
  upar.Add(string("tau_np_p_sgl0"),m_pdf->Get_tau_np_p_sgl(0),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_p_sgl1"),m_pdf->Get_tau_np_p_sgl(1),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_sgl0"),m_pdf->Get_tau_np_n_sgl(0),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_sgl1"),m_pdf->Get_tau_np_n_sgl(1),0.1,-10.,10.); NPar++;

  // * multiple Rnp */
  upar.Add(string("fd_np_mlt0"),m_pdf->Get_fd_np_mlt(0),0.1,-10.,10.); NPar++;
  upar.Add(string("fd_np_mlt1"),m_pdf->Get_fd_np_mlt(1),0.1,-1.,1.); NPar++;
  upar.Add(string("fd_np_st_mlt"),m_pdf->Get_fd_np_st_mlt(),0.1,-10.,10.); NPar++;
  upar.Add(string("fd_np_xi_mlt"),m_pdf->Get_fd_np_xi_mlt(),0.1,-1.,1.); NPar++;
  upar.Add(string("fd_np_stxi_mlt"),m_pdf->Get_fd_np_stxi_mlt(),0.1,-1.,1.); NPar++;
  upar.Add(string("fp_np_mlt"),m_pdf->Get_fp_np_mlt(),0.1,-1.,1.); NPar++;
  upar.Add(string("fn_np_mlt"),m_pdf->Get_fn_np_mlt(),0.1,-1.,1.); NPar++;
  upar.Add(string("tau_np_p_mlt0"),m_pdf->Get_tau_np_p_mlt(0),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_p_mlt1"),m_pdf->Get_tau_np_p_mlt(1),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_p_xi_mlt"),m_pdf->Get_tau_np_p_xi_mlt(),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_p_stxi_mlt"),m_pdf->Get_tau_np_p_stxi_mlt(),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_mlt0"),m_pdf->Get_tau_np_n_mlt(0),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_mlt1"),m_pdf->Get_tau_np_n_mlt(1),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_xi_mlt"),m_pdf->Get_tau_np_n_xi_mlt(),0.1,-10.,10.); NPar++;
  upar.Add(string("tau_np_n_stxi_mlt"),m_pdf->Get_tau_np_n_stxi_mlt(),0.1,-10.,10.); NPar++;
  upar.Add(string("s0"),0.001,0.1,0.,0.5);    NPar++;

  cout << NPar << " parameters" << endl;
  pdfFcn* theFCN = new pdfFcn();
  MnMigrad migrad(*theFCN,upar);

  int NFreePar = NPar;

  if(!include_s0) { migrad.Fix("s0"); NFreePar--;}

  migrad.Fix("fd_np_mlt1"); NFreePar--;
  migrad.Fix("fp_np_mlt"); NFreePar--; // == 0
  migrad.Fix("tau_np_n_sgl0"); NFreePar--;
  migrad.Fix("tau_np_n_sgl1"); NFreePar--;// == 0

  migrad.Fix("fd_np_sgl1"); NFreePar--;
  migrad.Fix("fd_np_sgl0"); NFreePar--;
  migrad.Fix("tau_np_p_sgl1"); NFreePar--;// == 0
  migrad.Fix("fp_np_sgl"); NFreePar--;
  if(!(make_Rnp_fit && make_sgl_fit)){
    migrad.Fix("tau_np_p_sgl0"); NFreePar--;
  }
  if(!(make_Rnp_fit && make_mlt_fit)){
    migrad.Fix("fd_np_mlt0"); NFreePar--;
//    migrad.Fix("fd_np_mlt1"); NFreePar--;
    migrad.Fix("fd_np_st_mlt"); NFreePar--;
    migrad.Fix("fd_np_xi_mlt"); NFreePar--;
    migrad.Fix("fd_np_stxi_mlt"); NFreePar--;
    migrad.Fix("fn_np_mlt"); NFreePar--;
    migrad.Fix("tau_np_p_mlt0"); NFreePar--;
    migrad.Fix("tau_np_p_mlt1"); NFreePar--;
    migrad.Fix("tau_np_p_xi_mlt"); NFreePar--;
    migrad.Fix("tau_np_p_stxi_mlt"); NFreePar--;
    migrad.Fix("tau_np_n_mlt0"); NFreePar--;
    migrad.Fix("tau_np_n_mlt1"); NFreePar--;
    migrad.Fix("tau_np_n_xi_mlt"); NFreePar--;
    migrad.Fix("tau_np_n_stxi_mlt"); NFreePar--;
  }
  if(!(make_Rdet_fit && make_sgl_fit)){
    migrad.Fix("Smn_asc"); NFreePar--;
    migrad.Fix("Stl_asc"); NFreePar--;
    migrad.Fix("ftl_asc"); NFreePar--;
  }
  if(!(make_Rdet_fit && make_mlt_fit)){
    migrad.Fix("Sasc0"); NFreePar--;
    migrad.Fix("Sasc1"); NFreePar--;
  }

  if(!(make_Olr_fit && make_mlt_fit)){
    migrad.Fix("f_ol_mlt"); NFreePar--;
  }
  if(!(make_Olr_fit && make_sgl_fit)){
    migrad.Fix("f_ol_sgl"); NFreePar--;
  }
  if(m_svd == 1){
    migrad.Fix("Stl_asc"); NFreePar--;
    migrad.Fix("ftl_asc"); NFreePar--;
  }
//  if(!make_Olr_fit){
    migrad.Fix("f_sigma_ol"); NFreePar--;
//  }

  cout <<  NFreePar << " free parameters" << endl;

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
  string str;
  if(m_svd == 2 && make_mlt_fit) str = string("Tag side resolution (mlt, SVD2");
  if(m_svd == 2 &&!make_mlt_fit) str = string("Tag side resolution (sgl, SVD2");
  if(m_svd == 1 && make_mlt_fit) str = string("Tag side resolution (mlt, SVD1");
  if(m_svd == 1 &&!make_mlt_fit) str = string("Tag side resolution (sgl, SVD1");
  switch(m_mode){
  case 1:
    str += string(", #pi^{0})");
    break;
  case 2:
    str += string(", #eta#rightarrow#gamma#gamma)");
    break;
  case 3:
    str += string(", #eta#rightarrow#pi^{+}#pi^{-}#pi^{0})");
    break;
  case 4:
    str += string(", #omega)");
    break;
  case 6:
    str += string(", D^{0}#pi^{+})");
    break;
  default:
    str += string(")");
    break;
  }
  gStyle->SetStatW(0.4); gStyle->SetStatH(0.3);
  TH1I* DH = new TH1I("DH",str.c_str(),NBins,dtmin,dtmax);

  DH->SetMarkerStyle(20);
  DH->SetMarkerSize(1.2);
  DH->SetMarkerColor(kBlack);
  DH->SetLineColor(kBlack);
  for(int i=0; i<NTot; i++){
    GetEvent(i);
    DH->Fill(m_dt);
  }

  SetPDFParams(upar.Params());
  const double ddt  = (dtmax-dtmin)/(double)NDots;
  const double ddtb = (dtmax-dtmin)/(double)NBins;
  double dt_arr[NDots];
  double dt_barr[NBins];
  double dt_barr_err[NBins];
  double dt_pull_err[NBins];

  cout << "Init hists and fill graphs" << endl;
  for(int j=0; j<NDots; j++) dt_arr[j] = dtmin+(j+0.5)*ddt;
  for(int j=0; j<NBins; j++){
    dt_barr_err[j] = 0;
    dt_pull_err[j] = 1;
    dt_barr[j] = dtmin+(j+0.5)*ddtb;
  }
  double Norm = 0;
  double NormB = 0;
  double pdf_array[NDots];
  double pull_array[NBins];
  const int nBD = NDots/NBins;
  double chisq = 0;
  for(int i=0; i<NDots; i++){
    pdf_array[i] = 0;
    const double& dt = dt_arr[i];
    for(int j=0; j<NTot; j++){
      GetEvent(j);
      if(include_s0) { sz_asc += upar.Value(NPar-1);}
      double pdfval = m_pdf->PdfRascRnp(dt,ntrk_asc,sz_asc,chisq_asc,ndf_asc);
      if(!std::isnan(pdfval)){
        pdf_array[i] += pdfval;
        Norm += pdf_array[i];
      }
    }
    pdf_array[i] *= ddtb;
    if(!(i%nBD)){
      const int bin = i/nBD;
      const int bin_content = DH->GetBinContent(bin+1);
//      cout << "bin " << bin << " " << bin_content << " " << pdf_array[i] << endl;
      if(bin_content){
        pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
        chisq += pull_array[bin]*pull_array[bin];
        NormB += pdf_array[i];
      }
    }
  }
  chisq /= NBins;
  cout << "chi2/n.d.f. = " << chisq << endl;
  //Norm *= ddt/NTot;
  cout << "Norm = " << Norm*ddt/NTot << ", Nevents = " << NTot << endl;

  TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
  GR->SetMarkerStyle(kDot);
  GR->SetMarkerSize(1.5);
  GR->SetMarkerColor(kBlue);
  GR->SetLineColor(kBlue);
  GR->SetLineWidth(2);

  TGraphErrors* GRP = new TGraphErrors(NBins,dt_barr,pull_array,dt_barr_err,dt_pull_err);
  GRP->SetMarkerStyle(20);
  GRP->SetMarkerSize(1.2);
  GRP->SetMarkerColor(kBlack);
  GRP->SetLineColor(kBlack);

  TCanvas* c1 = new TCanvas("c1","c1",600,800);
  TPad* pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.19);
  pad1->Draw(); pad2->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLogy();
  DH->GetXaxis()->SetTitle("#deltat (ps)");
  DH->GetXaxis()->SetTitleSize(0.06);
  DH->GetXaxis()->SetTitleOffset(0.75);
  DH->GetXaxis()->SetLabelSize(0.05);
  DH->GetYaxis()->SetLabelSize(0.05);
  DH->Draw("e1");

  stringstream out;
  TPaveText *pt = new TPaveText(0.60,0.53,0.98,0.6,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << chisq;
  pt->AddText(out.str().c_str());
  pt->Draw();

  GR->Draw("same");

  pad2->cd();
//  pad2->SetGrid();
  GRP->GetYaxis()->SetRangeUser(-5,5);
  GRP->GetXaxis()->SetRangeUser(dtmin,dtmax);
  GRP->Draw("ap");
  TLine* zeroline = new TLine(dtmin,0,dtmax,0);
  zeroline->SetLineColor(kBlue);
  zeroline->SetLineWidth(2);
  zeroline->Draw();

  TLine* plusline = new TLine(dtmin,3,dtmax,3);
  plusline->SetLineColor(kBlue);
  plusline->SetLineWidth(1);
  plusline->SetLineStyle(kDashed);
  plusline->Draw();

  TLine* minusline = new TLine(dtmin,-3,dtmax,-3);
  minusline->SetLineColor(kBlue);
  minusline->SetLineWidth(1);
  minusline->SetLineStyle(kDashed);
  minusline->Draw();

  out.str("");
  out << "res_tag_m" << m_mode;
  if(make_mlt_fit) out << "_mlt_";
  else             out << "_sgl_";
  if(m_svd == 2)   out << "svd2";
  else             out << "svd1";
  if(!NFreePar)    out << "_def";

  string rootname = out.str() + string(".root");
  c1->Print(rootname.c_str());
  string pngname = out.str() + string(".eps");
  c1->Print(pngname.c_str());

  pngname = out.str() + string(".png");
  c1->Print(pngname.c_str());
  out.str("");
  out << "display " << pngname << " &";
  system(out.str().c_str());
  return 0;
}

int main(const int argn, const char** argv){
//  TFile *ifile;
  string str;
  for(int i=1; i<argn; i++){
    str = string(argv[i]);
    if(string("-m") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_mode)) return -2;
      if(m_mode < 1 || m_mode > 6) return -3;
    }
    if(string("svd") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%d",&m_svd)) return -2;
      if(m_svd != 1 && m_svd != 2) return -3;
    }
    if(string("Rdet") == str){
      make_Rdet_fit = true;
    }
    if(string("Rnp") == str){
      make_Rnp_fit = true;
    }
    if(string("Olr") == str){
      make_Olr_fit = true;
    }
    if(string("sgl") == str){
      make_sgl_fit = true;
      make_mlt_fit = false;
    }
    if(string("mlt") == str){
      make_mlt_fit = true;
      make_sgl_fit = false;
    }
    if(string("dtmax") == str){
      if(++i == argn) return -1;
      if(1 != sscanf(argv[i],"%lf",&dtmax)) return -2;
      dtmax = fabs(dtmax);
      dtmin = -dtmax;
    }
    if(string("s0") == str){
      include_s0 = true;
    }
  }

  m_tree = new TChain("TEvent","TEvent");
  m_tree->Add(GenFile(m_mode).c_str());

//  m_pdf = m_mode<5 ? new RkRdetRnpPdf(1,m_svd) : new RkRdetRnpPdf(1,m_svd,true);
  m_pdf = new RkRdetRnpPdf(1,m_svd);
  init(m_mode);
  m_tree->SetBranchAddress("exp",&m_exp);
  m_tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  m_tree->SetBranchAddress("ndf_asc",&ndf_asc);
  m_tree->SetBranchAddress("sz_asc",&sz_asc);
  m_tree->SetBranchAddress("chisq_asc",&chisq_asc);
  m_tree->SetBranchAddress("z_asc",&m_z_asc);
  m_tree->SetBranchAddress("z_asc_mc",&m_z_asc_mc);
  m_tree->SetBranchAddress("b0f",&m_b0f);
  m_good_tree = GetGoodTTree(m_tree);

  const int NTot = m_good_tree->GetEntries();
//  m_good_tree->Print();

  make_full_fit();

  return 0;
}
