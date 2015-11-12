void drawParamsOmega(){
  TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
  c1->Divide(2,2);

  c1->cd(1);
  TGraphErrors* gr_dz0 = new TGraphErrors("ResolutionZSigOmega.txt","%lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  gr_dz0->SetTitle("Omega sig: dz_{0}");
  gr_dz0->SetMarkerStyle(21);
  gr_dz0->SetMarkerColor(kBlue);
  gr_dz0->Draw("ASP");

  c1->cd(2);
  TGraphErrors* gr_f = new TGraphErrors("ResolutionZSigOmega.txt","%lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg");
  gr_f->SetTitle("Omega sig: f_{s1}");
  gr_f->SetMarkerStyle(21);
  gr_f->SetMarkerColor(kBlue);
  gr_f->Draw("ASP");

  c1->cd(3);
  TGraphErrors* gr_s1 = new TGraphErrors("ResolutionZSigOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg");
  gr_s1->SetTitle("Omega sig: s1");
  gr_s1->SetMarkerStyle(21);
  gr_s1->SetMarkerColor(kBlue);
  gr_s1->Draw("ASP");

  c1->cd(4);
  TGraphErrors* gr_s2 = new TGraphErrors("ResolutionZSigOmega.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
  gr_s2->SetTitle("Omega sig: s2");
  gr_s2->SetMarkerStyle(21);
  gr_s2->SetMarkerColor(kBlue);
  gr_s2->Draw("ASP");

  c1->Update();
  c1->Print("c1.png");
  c1->Print("c1.root");

  return;
}
