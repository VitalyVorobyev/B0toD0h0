void drawPullParamsAscOmega(){
  TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
  c1->Divide(2,2);

  c1->cd(1);
  TGraphErrors* gr_dz0 = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gr_dz0->SetTitle("Omega asc: dz_{0}(pull)");
  gr_dz0->SetMarkerStyle(21);
  gr_dz0->SetMarkerColor(kBlue);
  gr_dz0->Draw("ASP");

  c1->cd(2);
  TGraphErrors* gr_f = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg");
  gr_f->SetTitle("Omega asc: f_{s1}(pull)");
  gr_f->SetMarkerStyle(21);
  gr_f->SetMarkerColor(kBlue);
  gr_f->Draw("ASP");

  c1->cd(3);
  TGraphErrors* gr_s1 = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg");
  gr_s1->SetTitle("Omega asc: s1(pull)");
  gr_s1->SetMarkerStyle(21);
  gr_s1->SetMarkerColor(kBlue);
  gr_s1->Draw("ASP");

  c1->cd(4);
  TGraphErrors* gr_s2 = new TGraphErrors("ResolutionZAscPullOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
  gr_s2->SetTitle("Omega asc: s2(pull)");
  gr_s2->SetMarkerStyle(21);
  gr_s2->SetMarkerColor(kBlue);
  gr_s2->Draw("ASP");

  c1->Update();
  c1->Print("c1.png");
  c1->Print("c1.root");

  return;
}
