void drawZShiftsOmega(){
  TCanvas* c1 = new TCanvas("c1","c1",600,600);

  TMultiGraph* mgr = new TMultiGraph();
  mgr->SetTitle("Omega: dz_{0}");

  TGraphErrors* gr_dz0_sig = new TGraphErrors("ZShiftSigOmega.txt","%lg %*lg %lg %lg");
  gr_dz0_sig->SetMarkerStyle(21);
  gr_dz0_sig->SetMarkerColor(kBlue);
  mgr->Add(gr_dz0_sig);

  TGraphErrors* gr_dz0_asc = new TGraphErrors("ZShiftAscOmega.txt","%lg %*lg %lg %lg");
  gr_dz0_asc->SetMarkerStyle(21);
  gr_dz0_asc->SetMarkerColor(kRed);
  mgr->Add(gr_dz0_asc);

  mgr->Draw("ASP");

  c1->Update();
  c1->Print("c1.png");
  c1->Print("c1.root");

  return;
}
