#ifndef MYTOOLS_H
#define MYTOOLS_H

using namespace RooFit;

bool IsInEllips(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde){
  const double xi_de = (de-de0)/Rde;
  const double xi_mbc = (mbc-mbc0)/Rmbc;
  if(xi_de*xi_de + xi_mbc*xi_mbc < 1) return true;
  else return false;
//  const double mbc_max_offset = Rmbc*sqrt(TMath::Abs(1-xi*xi));
//  if(TMath::Abs(mbc-mbc0) > mbc_max_offset) return false;
//  else return true;
}

double EllipsR2(const double& de, const double& mbc, const double& mbc0, const double& de0, const double& Rmbc, const double& Rde){
  const double xi_de = (de-de0)/Rde;
  const double xi_mbc = (mbc-mbc0)/Rmbc;
  return xi_de*xi_de + xi_mbc*xi_mbc;
}

int GetGoodTTree(TTree* tree, TTree* outtree, const int _mode){
  int mode, h0mode;
  switch(_mode){
  case 1: mode = 1; h0mode = 10; break; // D0 pi0
  case 2: mode = 2; h0mode = 10; break; // D0 eta->gg
  case 3: mode = 2; h0mode = 20; break; // D0 eta->pi+pi-pi0
  case 4: mode = 3; h0mode = 20; break; // D0 omega->pi+pi-pi0
  case 5: mode = 3; h0mode = 20; break; // D0 omega->pi+pi-pi0 (Ks rho model)
  case 6: mode = 4; h0mode = 40; break; // D0 rho
  default: return -1;
  }

  stringstream out;
  out.str("");
  out << "sz_sig>0 && ";
  out << "mode == " << mode << " && h0mode == " << h0mode << " && ";
  out << "chi2_vtx_d0<50 && ";
  out << "good_icpv == 1 && ";
  out << "(z_sig-z_asc)*0.1*78.48>" << dtmin << " && (z_sig-z_asc)*0.1*78.48<" << dtmax;
  const string sig_cut = out.str() + string(" && (b0f == 1 || b0f == 5 || b0f == 10)");
  outtree = tree->CopyTree(sig_cut.c_str());

  return 0;
}

#endif // MYTOOLS_H

//      /////////////
//      //  Plots  //
//      /////////////
//      // de //
//      RooPlot* deFrame = de.frame();
//      dataset[i][j]->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
//      new_pdf.plotOn(deFrame,Components(pdf_sig),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      new_pdf.plotOn(deFrame,Components(RooArgSet(pdf_peak,pdf_cmb_bb,pdf_cmb_qq)),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      new_pdf.plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

//      RooHist* hdepull = deFrame->pullHist();
//      RooPlot* dePull = de.frame(Title("#Delta E pull distribution"));
//      dePull->addPlotable(hdepull,"P");
//      dePull->GetYaxis()->SetRangeUser(-5,5);

//      out.str("");
//      out << "#Delta E, Signal " << Flv << Bin;
//      TCanvas* cm = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
//      cm->cd();

//      TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
//      TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
//      pad3->Draw();
//      pad4->Draw();

//      pad3->cd();
//      pad3->SetLeftMargin(0.15);
//      pad3->SetFillColor(0);

//      deFrame->GetXaxis()->SetTitleSize(0.05);
//      deFrame->GetXaxis()->SetTitleOffset(0.85);
//      deFrame->GetXaxis()->SetLabelSize(0.04);
//      deFrame->GetYaxis()->SetTitleOffset(1.6);
//      deFrame->Draw();

//      const int M_NSIGNAL_ELLI = m_ellitable1->get("signal") + m_ellitable1->get("fsr") + m_ellitable1->get("bad_pi0");

//      stringstream out1;
//      TPaveText *pt = new TPaveText(0.5,0.6,0.98,0.9,"brNDC");
//      pt->SetFillColor(0);
//      pt->SetTextAlign(12);
//      out1.str("");
//      out1 << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
//      pt->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << "S: " << (int)(m_nsigEl1+0.5) << " #pm " << (int)(m_nsig_errEl1_total+0.5) << " (" << M_NSIGNAL_ELLI << ")";
//      pt->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << "P: " << std::fixed << std::setprecision(2) << m_purityEl1*100. << " #pm " << m_purity_errEl1*100;
//      pt->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << label.c_str() << " (Bin " << Bin << ", flv " << Flv << ")";
//      pt->AddText(out1.str().c_str());
//      pt->Draw();

//      TLine *de_line_RIGHT;
//      de_line_RIGHT = new TLine(dE_max,0,dE_max,120);
//      de_line_RIGHT->SetLineColor(kRed);
//      de_line_RIGHT->SetLineStyle(1);
//      de_line_RIGHT->SetLineWidth((Width_t)2.);
//      de_line_RIGHT->Draw();
//      TLine *de_line_LEFT;
//      de_line_LEFT = new TLine(dE_min,0,dE_min,120);
//      de_line_LEFT->SetLineColor(kRed);
//      de_line_LEFT->SetLineStyle(1);
//      de_line_LEFT->SetLineWidth((Width_t)2.);
//      de_line_LEFT->Draw();

//      pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
//      dePull->SetMarkerSize(0.05); dePull->Draw();
//      TLine *de_lineUP = new TLine(deMin,3,deMax,3);
//      de_lineUP->SetLineColor(kBlue);
//      de_lineUP->SetLineStyle(2);
//      de_lineUP->Draw();
//      TLine *de_line = new TLine(deMin,0,deMax,0);
//      de_line->SetLineColor(kBlue);
//      de_line->SetLineStyle(1);
//      de_line->SetLineWidth((Width_t)2.);
//      de_line->Draw();
//      TLine *de_lineDOWN = new TLine(deMin,-3,deMax,-3);
//      de_lineDOWN->SetLineColor(kBlue);
//      de_lineDOWN->SetLineStyle(2);
//      de_lineDOWN->Draw();

//      cm->Update();

//      // mbc //
//      RooPlot* mbcFrame = mbc.frame();
//      dataset[i][j]->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
//      new_pdf.plotOn(mbcFrame,Components(RooArgSet(pdf_peak,pdf_cmb_bb,pdf_cmb_qq)),LineStyle(kDashed),ProjectionRange("deSignal"));
//      new_pdf.plotOn(mbcFrame,Components(pdf_sig),LineStyle(kDashed),ProjectionRange("deSignal"));
//      new_pdf.plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

//      RooHist* hmbcpull = mbcFrame->pullHist();
//      RooPlot* mbcPull  = mbc.frame(Title("#Delta E pull distribution"));
//      mbcPull->addPlotable(hmbcpull,"P");
//      mbcPull->GetYaxis()->SetRangeUser(-5,5);

//      out.str("");
//      out << "M_{bc}, Signal" << Flv << Bin;
//      TCanvas* cmmbc = new TCanvas(out.str().c_str(),out.str().c_str(),600,700);
//      cmmbc->cd();

//      TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
//      TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
//      pad1->Draw();
//      pad2->Draw();

//      pad1->cd();
//      pad1->SetLeftMargin(0.15);
//      pad1->SetFillColor(0);

//      mbcFrame->GetXaxis()->SetTitleSize(0.05);
//      mbcFrame->GetXaxis()->SetTitleOffset(0.85);
//      mbcFrame->GetXaxis()->SetLabelSize(0.04);
//      mbcFrame->GetYaxis()->SetTitleOffset(1.6);
//      mbcFrame->Draw();

//      TPaveText *ptmbc = new TPaveText(0.2,0.6,0.7,0.9,"brNDC");
//      ptmbc->SetFillColor(0);
//      ptmbc->SetTextAlign(12);
//      out1.str("");
//      out1 << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
//      ptmbc->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << "S: " << (int)(m_nsigEl1+0.5) << " #pm " << (int)(m_nsig_errEl1_total+0.5) << " (" << M_NSIGNAL_ELLI << ")";
//      ptmbc->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << "P: " << std::fixed << std::setprecision(2) << m_purityEl1*100. << " #pm " << m_purity_errEl1*100;
//      ptmbc->AddText(out1.str().c_str());
//      out1.str("");
//      out1 << label.c_str() << " (Bin " << Bin << ", flv " << Flv << ")";
//      ptmbc->AddText(out1.str().c_str());
//      ptmbc->Draw();

//      TLine *mbc_line_RIGHT;
//      mbc_line_RIGHT = new TLine(Mbc_max,0,Mbc_max,40);
//      mbc_line_RIGHT->SetLineColor(kRed);
//      mbc_line_RIGHT->SetLineStyle(1);
//      mbc_line_RIGHT->SetLineWidth((Width_t)2.);
//      mbc_line_RIGHT->Draw();
//      TLine *mbc_line_LEFT;
//      mbc_line_LEFT = new TLine(Mbc_min,0,Mbc_min,40);
//      mbc_line_LEFT->SetLineColor(kRed);
//      mbc_line_LEFT->SetLineStyle(1);
//      mbc_line_LEFT->SetLineWidth((Width_t)2.);
//      mbc_line_LEFT->Draw();

//      pad2->cd();
//      pad2->SetLeftMargin(0.15);
//      pad2->SetFillColor(0);
//      mbcPull->SetMarkerSize(0.05);
//      mbcPull->Draw();
//      TLine *mbc_lineUP = new TLine(mbcMin,3,mbcMax,3);
//      mbc_lineUP->SetLineColor(kBlue);
//      mbc_lineUP->SetLineStyle(2);
//      mbc_lineUP->Draw();
//      TLine *mbc_line = new TLine(mbcMin,0,mbcMax,0);
//      mbc_line->SetLineColor(kBlue);
//      mbc_line->SetLineStyle(1);
//      mbc_line->SetLineWidth((Width_t)2.);
//      mbc_line->Draw();
//      TLine *mbc_lineDOWN = new TLine(mbcMin,-3,mbcMax,-3);
//      mbc_lineDOWN->SetLineColor(kBlue);
//      mbc_lineDOWN->SetLineStyle(2);
//      mbc_lineDOWN->Draw();

//      cmmbc->Update();
