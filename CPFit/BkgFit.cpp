#include "BkgFit.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPaveText.h"

MnUserParameterState make_single_fit(void){
  MnUserParameters upar;
  return upar;
}

int make_full_fit(void){
    stringstream out;
    int NPar = 0;
    MnUserParameters upar;
    upar.Add(string("tau"),         m_pdf->GetTau(),          0.2, 0  ,10);  NPar++;
    upar.Add(string("mu"),          m_pdf->Get_mu(),          0.2,-10.,10.); NPar++;
    upar.Add(string("mu_delta"),    m_pdf->Get_mu_delta(),    0.2,-10.,10.); NPar++;
    upar.Add(string("f_delta_mlt"), m_pdf->Get_f_delta_mlt(), 0.2, 0. ,1.);  NPar++;
    upar.Add(string("f_tail_mlt"),  m_pdf->Get_f_tail_mlt(),  0.2, 0. ,1.);  NPar++;
    upar.Add(string("S_main_mlt"),  m_pdf->Get_S_main_mlt(),  0.2, 0. ,10.); NPar++;
    upar.Add(string("S_tail_mlt"),  m_pdf->Get_S_tail_mlt(),  0.2, 0  ,10.); NPar++;
    upar.Add(string("f_delta_sgl"), m_pdf->Get_f_delta_sgl(), 0.2, 0. ,1.);  NPar++;
    upar.Add(string("f_tail_sgl"),  m_pdf->Get_f_tail_sgl(),  0.2, 0. ,1.);  NPar++;
    upar.Add(string("S_main_sgl"),  m_pdf->Get_S_main_sgl(),  0.2, 0  ,10.); NPar++;
    upar.Add(string("S_tail_sgl"),  m_pdf->Get_S_tail_sgl(),  0.2, 0  ,10.); NPar++;

    pdfFcn* theFCN = new pdfFcn();
    MnMigrad migrad(*theFCN,upar);

    if(fix_mlt){
      if(!fix_all) migrad.Fix("f_delta_mlt");
      migrad.Fix("f_tail_mlt");
      migrad.Fix("S_main_mlt");
      migrad.Fix("S_tail_mlt");
    }
    if(fix_sgl){
      if(!fix_all) migrad.Fix("f_delta_sgl");
      migrad.Fix("f_tail_sgl");
      migrad.Fix("S_main_sgl");
      migrad.Fix("S_tail_sgl");
    }
    if(fix_all){
      migrad.Fix("tau");
      migrad.Fix("mu");
      migrad.Fix("mu_delta");
    }
    if(fix_delta){
      migrad.Fix("mu_delta");
//      migrad.Fix("f_delta_mlt");
//      migrad.Fix("f_delta_sgl");
    }
    if(fix_lt){
      migrad.Fix("tau");
      migrad.Fix("mu");
      migrad.Fix("f_tail_mlt");
      migrad.Fix("S_main_mlt");
      migrad.Fix("S_tail_mlt");
      migrad.Fix("f_tail_sgl");
      migrad.Fix("S_main_sgl");
      migrad.Fix("S_tail_sgl");
    }
    if(fix_main){
      migrad.Fix("S_main_mlt");
      migrad.Fix("S_main_sgl");
    }
    if(fix_tail){
      migrad.Fix("S_tail_mlt");
      migrad.Fix("S_tail_sgl");
    }

    cout << "Starting minimization" << endl;
    FunctionMinimum min = migrad();
    MnUserParameterState pstate = min.UserState();
    cout << "after migrad" << endl;

//    vector<double> pars_vec;
    for(int i=0; i<NPar; i++){
      cout << upar.Name(i) << ": " << upar.Value(i) << " -> ";
      upar.SetValue(i,pstate.Value(i));
      cout << upar.Value(i) << " +- " << pstate.Error(i) << endl;
//      pars_vec.push_back(upar.Value(i));
    }

    const int NTot = m_good_tree->GetEntries();

    const double ddt = (dtmax-dtmin)/(double)NDots;
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

    SetPdfParams(upar.Params());

    if(!separate_plots){
      gStyle->SetStatW(0.4); gStyle->SetStatH(0.3);
      stringstream out;
      out.str("");
      out << "Background #Deltat fit";
      if(only_mlt && m_svd == 2)      out << " (mlt,svd2)";
      else if(only_mlt && m_svd == 1) out << " (mlt,svd1)";
      else if(only_sgl && m_svd == 2) out << " (sgl,svd2)";
      else if(only_sgl && m_svd == 1) out << " (sgl,svd1)";
      else if(m_svd == 2)             out << " (svd2)";
      else                            out << " (svd1)";
      TH1I* DH = new TH1I("DH",out.str().c_str(),NBins,dtmin,dtmax);
      DH->SetMarkerStyle(20);
      DH->SetMarkerColor(kBlack);
      DH->SetLineColor(kBlack);
      DH->SetMarkerSize(1.1);

      for(int i=0; i<NTot; i++){
        GetEvent(i);
        DH->Fill(m_dt);
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
          double pdfval = m_pdf->Pdf(dt,sum_sigma(sz_rec,sz_asc),ndf_asc);
          if(!std::isnan(pdfval)){
            pdf_array[i] += pdfval;
            Norm += pdf_array[i];
          }
        }
        pdf_array[i] *= ddtb;
        if(!(i%nBD)){
          const int bin = i/nBD;
          const int bin_content = DH->GetBinContent(bin+1);
//          cout << "bin " << bin << " " << bin_content << " " << pdf_array[i] << endl;
          if(bin_content){
            pull_array[bin] = (bin_content - pdf_array[i])/sqrt(bin_content);
            chisq += pull_array[bin]*pull_array[bin];
            NormB += pdf_array[i];
          }
        }
      }
      chisq /= NBins;
      cout << "chi2/n.d.f. = " << chisq << endl;
//      Norm *= ddt/NTot;
//      cout << "Norm = " << Norm << "Nevents = " << NTot << endl;

      TGraph* GR = new TGraph(NDots,dt_arr,pdf_array);
      GR->SetMarkerStyle(kDot);
      GR->SetMarkerSize(1.5);
      GR->SetLineWidth(2);
      GR->SetMarkerColor(kBlue);
      GR->SetLineColor(kBlue);

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
      DH->GetXaxis()->SetTitle("#Deltat (ps)");
      DH->GetXaxis()->SetTitleSize(0.06);
      DH->GetXaxis()->SetTitleOffset(0.75);
      DH->GetXaxis()->SetLabelSize(0.05);
      DH->GetYaxis()->SetLabelSize(0.05);
      DH->Draw("e1");

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
      out << "fit_bkg";
      if(only_mlt)   out << "_mlt";
      if(only_sgl)   out << "_sgl";
      if(m_svd == 2) out << "_svd2";
      else           out << "_svd1";

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
    else{
    int Nev[2][16];
    TH1I* dh[2][16];
    double pdf_arr[2][16][NDots];
    double norm[2][16];
    double norm1[2][16];
    for(int k=0; k<2; k++){
      cout << "flv = " << flv(k) << endl;
      for(int i=0; i<16; i++){
        Nev[k][i]    = 0;
        norm[k][i]   = 0;
        norm1[k][i]  = 0;
        out.str("");
        out << "dh" << bin(i) << "_" << k;
        dh[k][i]  = new TH1I(out.str().c_str(),out.str().c_str(),NBins,dtmin,dtmax);
        m_bin = bin(i);
        cout << "  bin = " << m_bin << endl;
        for(int j=0; j<NDots; j++){
          const double& dt = dt_arr[j];
          pdf_arr[k][i][j] = 0;
          for(int l=0; l<NTot; l++){
            GetEvent(l);
            pdf_arr[k][i][j] += m_pdf->Pdf(dt,sum_sigma(sz_rec,sz_asc),ndf_asc);
            norm[k][i] += pdf_arr[k][i][j];
          }
          pdf_arr[k][i][j] *= ddt/NTot;
        }
        norm[k][i] *= ddt/NTot;
      }
    }

    int NEveCounter = 0;
    for(int i=0; i<NTot; i++){
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
    c2->Print("genb_fit.root");
    c2->Print("genb_fit.png");

    system("display genb_fit.png &");

    cout << "NEveCounter = " << NEveCounter << endl;
    }
    return 0;
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
  ifile  = TFile::Open("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m1_mh010.root");

  m_pdf = new RbkgPdf(_mode,m_svd);
  init(_mode);
  m_tree = (TTree*)ifile->Get("TEvent");
  m_tree->SetBranchAddress("exp",&m_exp);
  m_tree->SetBranchAddress("z_sig", &m_z_sig);
  m_tree->SetBranchAddress("z_asc", &m_z_asc);
  m_tree->SetBranchAddress("bin",&m_bin);
  m_tree->SetBranchAddress("tag_LH",&m_tag_LH);
  m_tree->SetBranchAddress("ndf_asc",&ndf_asc);
  m_tree->SetBranchAddress("sz_sig",&sz_rec);
  m_tree->SetBranchAddress("sz_asc",&sz_asc);
  m_good_tree = GetGoodTTree(m_tree);

  int j,k;

  vector< vector<double> > vals_vec, errs_vec;
  vector<double> Aref,Bref;

  if(argn == 4){
    if(abs(_binfit)>8) return -1;
    if(_binfit){
      m_good_tree = GetGoodTTree(m_tree,_binfit,_flv);
      MnUserParameterState upar = make_single_fit();
      vals_vec.push_back(upar.Params());
      errs_vec.push_back(upar.Errors());
    } else{
      draw_plots = false;
      m_flv = _flv;
      for(int i=0; i<16; i++){
        m_bin = bin(i);
        m_good_tree = GetGoodTTree(m_tree,m_bin,m_flv);
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
        m_good_tree = GetGoodTTree(m_tree,m_bin,m_flv);
        MnUserParameterState upar = make_single_fit();
        vals_vec.push_back(upar.Params());
        errs_vec.push_back(upar.Errors());
      }
    }
  } else {
    m_good_tree = GetGoodTTree(m_tree);
    make_full_fit();
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
