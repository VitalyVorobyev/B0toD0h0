void drawPullParamsOmega(){
  TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
  c1->Divide(2,2);

  c1->cd(1);
  TMultiGraph* mgr_dz0 = new TMultiGraph();
  mgr_dz0->SetTitle("Omega: dz_{0}(pull)");

  TGraphErrors* gr_dz0_sig = new TGraphErrors("ResolutionZSigPullOmega.txt","%lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gr_dz0_sig->SetMarkerStyle(21);
  gr_dz0_sig->SetMarkerColor(kBlue);
  mgr_dz0->Add(gr_dz0_sig);

  TGraphErrors* gr_dz0_asc = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gr_dz0_asc->SetMarkerStyle(21);
  gr_dz0_asc->SetMarkerColor(kRed);
  mgr_dz0->Add(gr_dz0_asc);

  mgr_dz0->Draw("ASP");

  c1->cd(2);
  TMultiGraph* mgr_f = new TMultiGraph();
  mgr_f->SetTitle("Omega: f(s1)(pull)");

  TGraphErrors* gr_f_sig = new TGraphErrors("ResolutionZSigPullOmega.txt","%lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg");
  gr_f_sig->SetMarkerStyle(21);
  gr_f_sig->SetMarkerColor(kBlue);
  mgr_f->Add(gr_f_sig);

  TGraphErrors* gr_f_asc = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg");
  gr_f_asc->SetMarkerStyle(21);
  gr_f_asc->SetMarkerColor(kRed);
  mgr_f->Add(gr_f_asc);

  mgr_f->Draw("ASP");

  c1->cd(3);
  TMultiGraph* mgr_s1 = new TMultiGraph();
  mgr_s1->SetTitle("Omega: s1(pull)");

  TGraphErrors* gr_s1_sig = new TGraphErrors("ResolutionZSigPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg");
  gr_s1_sig->SetMarkerStyle(21);
  gr_s1_sig->SetMarkerColor(kBlue);
  mgr_s1->Add(gr_s1_sig);

  TGraphErrors* gr_s1_asc = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg");
  gr_s1_asc->SetMarkerStyle(21);
  gr_s1_asc->SetMarkerColor(kRed);
  mgr_s1->Add(gr_s1_asc);

  mgr_s1->Draw("ASP");

  c1->cd(4);
  TMultiGraph* mgr_s2 = new TMultiGraph();
  mgr_s2->SetTitle("Omega: s2(pull)");

  TGraphErrors* gr_s2_sig = new TGraphErrors("ResolutionZSigPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
  gr_s2_sig->SetMarkerStyle(21);
  gr_s2_sig->SetMarkerColor(kBlue);
  mgr_s2->Add(gr_s2_sig);

  TGraphErrors* gr_s2_asc = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
  gr_s2_asc->SetMarkerStyle(21);
  gr_s2_asc->SetMarkerColor(kRed);
  mgr_s2->Add(gr_s2_asc);

  mgr_s2->Draw("ASP");


  c1->Update();
  c1->Print("c1.png");
  c1->Print("c1.root");

  return;
}
