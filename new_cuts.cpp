string GetElliCut(const double& demin,const double& demax,const double& mbcmin,const double& mbcmax){
  stringstream out;
  const double de0    = 0.5*(demax+demin);
  const double mbc0   = 0.5*(mbcmax+mbcmin);
  const double dersq  = 0.25*(demax-demin)*(demax-demin);
  const double mbcrsq = 0.25*(mbcmax-mbcmin)*(mbcmax-mbcmin);
  out << "((" << de0 << "-de)*(" << de0 << "-de)/" << dersq;
  out << "+(" << mbc0 << "-mbc)*(" << mbc0<<"-mbc)/" << mbcrsq;
  out << ")<1";
  return out.str();
}

const string prefix("Tuples/");

void new_cuts_pi0(void){
  stringstream out;
  const string ellicut = GetElliCut(-0.1,0.0921343,5.27093,5.28789);
  const string bdtcut("bdt>0.27");
//  const string bdtcut("bdt>-1");
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10)");
  const string bkgcut("!(b0f == 1 || b0f == 5 || b0f == 10)");

  TChain* tree = new TChain("TEvent","TEvent");
  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
  int streams[6] = {0,1,2,3,4,5};
  for(int i=0; i<4; i++){
    for(int j=0; j<6; j++){
      out.str("");
      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m1_h0m10.root";
      tree->Add(out.str().c_str());
    }
  }
  double de,mbc,bdt,pt_pip,pt_pim,e_g1;
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("bdt",&bdt);
  tree->SetBranchAddress("pt_pim",&pt_pim);
  tree->SetBranchAddress("pt_pip",&pt_pip);
  tree->SetBranchAddress("e_g1",&e_g1);

  TGraph* gr_fom_g1 = new TGraph(20);
  double S,B,FOM;
  double eg1cut;
  for(int i=0; i<20; i++){
    eg1cut = 0.04+i*0.005;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "e_g1>" << eg1cut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "e_g1>" << eg1cut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_g1->SetPoint(i,eg1cut,FOM);
  }

  gr_fom_g1->SetMarkerColor(kBlue);
  gr_fom_g1->SetMarkerStyle(20);
  gr_fom_g1->SetMarkerSize(1.3);

  TCanvas* c1 = new TCanvas("c1","c1",600,400);
  c1->cd();
  c1->SetGrid();
  c1->Draw();
  gr_fom_g1->Draw("ap");

  cout << endl;

  TGraph* gr_fom_pt = new TGraph(20);
  double ptcut;
  for(int i=0; i<20; i++){
    ptcut = 0.04+i*0.005;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& pt_pip>" << ptcut << "&& pt_pim>" << ptcut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& pt_pip>" << ptcut << "&& pt_pim>" << ptcut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_pt->SetPoint(i,ptcut,FOM);
  }

  gr_fom_pt->SetMarkerColor(kBlue);
  gr_fom_pt->SetMarkerStyle(20);
  gr_fom_pt->SetMarkerSize(1.3);

  TCanvas* c2 = new TCanvas("c2","c2",600,400);
  c2->cd();
  c2->SetGrid();
  c2->Draw();
  gr_fom_pt->Draw("ap");
  return;
}

void new_cuts_dstpi0(void){
  stringstream out;
  const string ellicut = GetElliCut(-0.1,0.0921343,5.27093,5.28789);
  const string bdtcut("lh0>0.60");
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10)");
  const string bkgcut("!(b0f == 1 || b0f == 5 || b0f == 10)");

  TChain* tree = new TChain("TEvent","TEvent");
  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
  int streams[6] = {0,1,2,3,4,5};
  for(int i=0; i<4; i++){
    for(int j=0; j<6; j++){
      out.str("");
      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m10_h0m10.root";
      tree->Add(out.str().c_str());
    }
  }
  double de,mbc,lh0,pt_pip,pt_pim,e_g1,e_g3;
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("lh0",&lh0);
  tree->SetBranchAddress("pt_pim",&pt_pim);
  tree->SetBranchAddress("pt_pip",&pt_pip);
  tree->SetBranchAddress("e_g1",&e_g1);
  tree->SetBranchAddress("e_g3",&e_g3);

  double S,B,FOM;
  double eg1cut;
//  TGraph* gr_fom_g1 = new TGraph(20);
//  for(int i=0; i<20; i++){
//    eg1cut = 0.04+i*0.005;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "e_g1>" << eg1cut;
//    S = tree->GetEntries(out.str().c_str())/6;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "e_g1>" << eg1cut;
//    B = tree->GetEntries(out.str().c_str())/6;
//    FOM = S/sqrt(S+B);
//    cout << S << " " << B << " " << FOM << endl;
//    gr_fom_g1->SetPoint(i,eg1cut,FOM);
//  }

//  gr_fom_g1->SetMarkerColor(kBlue);
//  gr_fom_g1->SetMarkerStyle(20);
//  gr_fom_g1->SetMarkerSize(1.3);

//  TCanvas* c1 = new TCanvas("c1","c1",600,400);
//  c1->cd();
//  c1->SetGrid();
//  c1->Draw();
//  gr_fom_g1->Draw("ap");

//  cout << endl;

  TGraph* gr_fom_g3 = new TGraph(30);
  for(int i=0; i<30; i++){
    eg1cut = 0.04+i*0.002;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "e_g3>" << eg1cut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "e_g3>" << eg1cut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_g3->SetPoint(i,eg1cut,FOM);
  }

  gr_fom_g3->SetMarkerColor(kBlue);
  gr_fom_g3->SetMarkerStyle(20);
  gr_fom_g3->SetMarkerSize(1.3);

  TCanvas* c3 = new TCanvas("c3","c3",600,400);
  c3->cd();
  c3->SetGrid();
  c3->Draw();
  gr_fom_g3->Draw("ap");

  cout << endl;

  return;
}

void new_cuts_etagg(void){
  stringstream out;
  const string ellicut = GetElliCut(-0.0994528,0.0758825,5.27147,5.28763);
  const string bdtcut("bdt>0.15");
//  const string bdtcut("bdt>-1");
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10)");
  const string bkgcut("!(b0f == 1 || b0f == 5 || b0f == 10)");

  TChain* tree = new TChain("TEvent","TEvent");
  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
  int streams[6] = {0,1,2,3,4,5};
  for(int i=0; i<4; i++){
    for(int j=0; j<6; j++){
      out.str("");
      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m2_h0m10.root";
      tree->Add(out.str().c_str());
    }
  }
  double de,mbc,bdt,pt_pip,pt_pim,e_g1;
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("bdt",&bdt);
  tree->SetBranchAddress("pt_pim",&pt_pim);
  tree->SetBranchAddress("pt_pip",&pt_pip);
  tree->SetBranchAddress("e_g1",&e_g1);

  TGraph* gr_fom_g1 = new TGraph(30);
  double S,B,FOM;
  double eg1cut;
  for(int i=0; i<30; i++){
    eg1cut = 0.08+i*0.01;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "e_g1>" << eg1cut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "e_g1>" << eg1cut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_g1->SetPoint(i,eg1cut,FOM);
  }

  gr_fom_g1->SetMarkerColor(kBlue);
  gr_fom_g1->SetMarkerStyle(20);
  gr_fom_g1->SetMarkerSize(1.3);

  TCanvas* c1 = new TCanvas("c1","c1",600,400);
  c1->cd();
  c1->SetGrid();
  c1->Draw();
  gr_fom_g1->Draw("ap");

//  cout << endl;

//  TGraph* gr_fom_pt = new TGraph(20);
//  double ptcut;
//  for(int i=0; i<20; i++){
//    ptcut = 0.04+i*0.005;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& pt_pip>" << ptcut << "&& pt_pim>" << ptcut;
//    S = tree->GetEntries(out.str().c_str())/6;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& pt_pip>" << ptcut << "&& pt_pim>" << ptcut;
//    B = tree->GetEntries(out.str().c_str())/6;
//    FOM = S/sqrt(S+B);
//    cout << S << " " << B << " " << FOM << endl;
//    gr_fom_pt->SetPoint(i,ptcut,FOM);
//  }

//  gr_fom_pt->SetMarkerColor(kBlue);
//  gr_fom_pt->SetMarkerStyle(20);
//  gr_fom_pt->SetMarkerSize(1.3);

//  TCanvas* c2 = new TCanvas("c2","c2",600,400);
//  c2->cd();
//  c2->SetGrid();
//  c2->Draw();
//  gr_fom_pt->Draw("ap");
  return;
}

void new_cuts_etappp(void){
  stringstream out;
  const string ellicut = GetElliCut(-0.0534644,0.0453713,5.27205,5.28735);
  const string bdtcut("bdt>0.18");
//  const string bdtcut("bdt>-1");
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10)");
  const string bkgcut("!(b0f == 1 || b0f == 5 || b0f == 10)");

  TChain* tree = new TChain("TEvent","TEvent");
  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
  int streams[6] = {0,1,2,3,4,5};
  for(int i=0; i<4; i++){
    for(int j=0; j<6; j++){
      out.str("");
      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m2_h0m20.root";
      tree->Add(out.str().c_str());
    }
  }
  double de,mbc,bdt,pt_pip,pt_pim,e_g1,p_pi0;
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("bdt",&bdt);
  tree->SetBranchAddress("pt_pi1",&pt_pim);
  tree->SetBranchAddress("pt_pi2",&pt_pip);
  tree->SetBranchAddress("p_pi0_h0",&p_pi0);
  tree->SetBranchAddress("e_g1",&e_g1);

  TGraph* gr_fom_g1 = new TGraph(30);
  double S,B,FOM;
  double eg1cut;
  for(int i=0; i<30; i++){
    eg1cut = 0.04+i*0.005;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "e_g1>" << eg1cut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "e_g1>" << eg1cut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_g1->SetPoint(i,eg1cut,FOM);
  }

  gr_fom_g1->SetMarkerColor(kBlue);
  gr_fom_g1->SetMarkerStyle(20);
  gr_fom_g1->SetMarkerSize(1.3);

  TCanvas* c1 = new TCanvas("c1","c1",600,400);
  c1->cd();
  c1->SetGrid();
  c1->Draw();
  gr_fom_g1->Draw("ap");

  cout << endl;

  TGraph* gr_fom_ppi0 = new TGraph(30);
  double ppi0cut;
  for(int i=0; i<30; i++){
    ppi0cut = 0.0+i*0.02;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "p_pi0_h0>" << ppi0cut;
//    cout << out.str() << endl;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "p_pi0_h0>" << ppi0cut;
//    cout << out.str() << endl;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_ppi0->SetPoint(i,ppi0cut,FOM);
  }

  gr_fom_ppi0->SetMarkerColor(kBlue);
  gr_fom_ppi0->SetMarkerStyle(20);
  gr_fom_ppi0->SetMarkerSize(1.3);

  TCanvas* c3 = new TCanvas("c3","c3",600,400);
  c3->cd();
  c3->SetGrid();
  c3->Draw();
  gr_fom_ppi0->Draw("ap");

  cout << endl;

  TGraph* gr_fom_pt = new TGraph(30);
  double ptcut;
  for(int i=0; i<30; i++){
    ptcut = 0.04+i*0.01;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& pt_pi1>" << ptcut << "&& pt_pi2>" << ptcut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& pt_pi1>" << ptcut << "&& pt_pi2>" << ptcut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_pt->SetPoint(i,ptcut,FOM);
  }

  gr_fom_pt->SetMarkerColor(kBlue);
  gr_fom_pt->SetMarkerStyle(20);
  gr_fom_pt->SetMarkerSize(1.3);

  TCanvas* c2 = new TCanvas("c2","c2",600,400);
  c2->cd();
  c2->SetGrid();
  c2->Draw();
  gr_fom_pt->Draw("ap");
  return;
}

void new_cuts_omega(void){
  stringstream out;
  const string ellicut = GetElliCut(-0.0565915,0.0467748,5.27198,5.28743);
  const string bdtcut("bdt>0.18");
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10)");
  const string bkgcut("!(b0f == 1 || b0f == 5 || b0f == 10)");

  TChain* tree = new TChain("TEvent","TEvent");
  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
  int streams[6] = {0,1,2,3,4,5};
  for(int i=0; i<4; i++){
    for(int j=0; j<6; j++){
      out.str("");
      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m3_h0m20.root";
      tree->Add(out.str().c_str());
    }
  }
  double de,mbc,bdt,pt_pip,pt_pim,e_g1,p_pi0,cos_hel,th_g1;
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("bdt",&bdt);
  tree->SetBranchAddress("pt_pi1",&pt_pim);
  tree->SetBranchAddress("pt_pi2",&pt_pip);
  tree->SetBranchAddress("p_pi0_h0",&p_pi0);
  tree->SetBranchAddress("e_g1",&e_g1);
  tree->SetBranchAddress("cos_hel",&cos_hel);
  tree->SetBranchAddress("th_g1",&th_g1);

  double S,B,FOM;
  double eg1cut;
//  TGraph* gr_fom_g1_bar = new TGraph(30);
//  for(int i=0; i<30; i++){
//    eg1cut = 0.04+i*0.005;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& abs(th_g1)<0.7 && e_g1>" << eg1cut;
//    S = tree->GetEntries(out.str().c_str())/6;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& abs(th_g1)<0.7 && e_g1>" << eg1cut;
//    B = tree->GetEntries(out.str().c_str())/6;
//    FOM = S/sqrt(S+B);
//    cout << S << " " << B << " " << FOM << endl;
//    gr_fom_g1_bar->SetPoint(i,eg1cut,FOM);
//  }

//  gr_fom_g1_bar->SetMarkerColor(kBlue);
//  gr_fom_g1_bar->SetMarkerStyle(20);
//  gr_fom_g1_bar->SetMarkerSize(1.3);

//  TCanvas* c1 = new TCanvas("c1bar","c1bar",600,400);
//  c1->cd();
//  c1->SetGrid();
//  c1->Draw();
//  gr_fom_g1_bar->Draw("ap");

//  cout << endl;

//  TGraph* gr_fom_g1_end = new TGraph(30);
//  for(int i=0; i<30; i++){
//    eg1cut = 0.04+i*0.005;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& abs(th_g1)>0.7 && e_g1>" << eg1cut;
//    S = tree->GetEntries(out.str().c_str())/6;
//    out.str("");
//    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& abs(th_g1)>0.7 && e_g1>" << eg1cut;
//    B = tree->GetEntries(out.str().c_str())/6;
//    FOM = S/sqrt(S+B);
//    cout << S << " " << B << " " << FOM << endl;
//    gr_fom_g1_end->SetPoint(i,eg1cut,FOM);
//  }

//  gr_fom_g1_end->SetMarkerColor(kBlue);
//  gr_fom_g1_end->SetMarkerStyle(20);
//  gr_fom_g1_end->SetMarkerSize(1.3);

//  TCanvas* c4 = new TCanvas("c4end","c4end",600,400);
//  c4->cd();
//  c4->SetGrid();
//  c4->Draw();
//  gr_fom_g1_end->Draw("ap");

//  cout << endl;

  TGraph* gr_fom_ppi0 = new TGraph(30);
  double ppi0cut;
  for(int i=0; i<30; i++){
    ppi0cut = 0.0+i*0.02;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&&" << "p_pi0_h0>" << ppi0cut;
//    cout << out.str() << endl;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&&" << "p_pi0_h0>" << ppi0cut;
//    cout << out.str() << endl;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_ppi0->SetPoint(i,ppi0cut,FOM);
  }

  gr_fom_ppi0->SetMarkerColor(kBlue);
  gr_fom_ppi0->SetMarkerStyle(20);
  gr_fom_ppi0->SetMarkerSize(1.3);

  TCanvas* c3 = new TCanvas("c3","c3",600,400);
  c3->cd();
  c3->SetGrid();
  c3->Draw();
  gr_fom_ppi0->Draw("ap");

  cout << endl;

  TGraph* gr_fom_pt = new TGraph(30);
  double ptcut;
  for(int i=0; i<30; i++){
    ptcut = 0.04+i*0.01;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& pt_pi1>" << ptcut << "&& pt_pi2>" << ptcut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& pt_pi1>" << ptcut << "&& pt_pi2>" << ptcut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_pt->SetPoint(i,ptcut,FOM);
  }

  gr_fom_pt->SetMarkerColor(kBlue);
  gr_fom_pt->SetMarkerStyle(20);
  gr_fom_pt->SetMarkerSize(1.3);

  TCanvas* c2 = new TCanvas("c2","c2",600,400);
  c2->cd();
  c2->SetGrid();
  c2->Draw();
  gr_fom_pt->Draw("ap");

  TGraph* gr_fom_hel = new TGraph(30);
  double helcut;
  for(int i=0; i<30; i++){
    helcut = 0.1+i*0.01;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << sigcut << "&& abs(cos_hel)>" << helcut;
    S = tree->GetEntries(out.str().c_str())/6;
    out.str("");
    out << ellicut << "&&" << bdtcut << "&&" << bkgcut << "&& abs(cos_hel)>" << helcut;
    B = tree->GetEntries(out.str().c_str())/6;
    FOM = S/sqrt(S+B);
    cout << S << " " << B << " " << FOM << endl;
    gr_fom_hel->SetPoint(i,helcut,FOM);
  }

  gr_fom_hel->SetMarkerColor(kBlue);
  gr_fom_hel->SetMarkerStyle(20);
  gr_fom_hel->SetMarkerSize(1.3);

  TCanvas* c5 = new TCanvas("c5","c5",600,400);
  c5->cd();
  c5->SetGrid();
  c5->Draw();
  gr_fom_hel->Draw("ap");
  return;
}
