int Mode(const int mode){
  switch (mode) {
  case 1:  return 1; // pi0
  case 10: return 10;// D*0 pi0
  case 20: return 20;// D*0 eta
  case 4:  return 3; // omega
  case 5:  return 5; // eta'
  case 77: return 77; // eta'
  default: return 2; // eta
  }
}

int h0Mode(const int mode){
  switch (mode) {
  case 1:  return 10; // pi0
  case 10: return 10; // D*0 pi0
  case 2:  return 10; // eta->gg
  case 20: return 10; // D*0 eta
  case 3:  return 20; // eta->ppp
  case 4:  return 20; // omega
  case 5:  return 10; // eta'
  case 77: return 77; // eta'
  }
}

string GenWWFile(const int mode, const int flag = 0){
    // flag == 0 -> all events
    // flag == 1 -> continuum
    // flag == 2 -> BB
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/PurityFit/data/genmc_cpv_tree_m" << Mode(mode) << "_hm" << h0Mode(mode);
  if(flag == 1)      out << "_cont";
  else if(flag == 2) out << "_BB";
  out << "_ww.root";
  return out.str();
}

void draw_de_mbc_regions(const int mode){
  TChain* tree = new TChain("TEvent");
  tree->Add(GenWWFile(mode).c_str());

  TCanvas* ellican = new TCanvas("ellican","ellican",800,600);
  ellican->cd();
  stringstream out;
  out.str("mbc>5.2");
  tree->Draw("de:mbc",out.str().c_str());
  tree->SetMarkerStyle(6);

  tree->SetMarkerColor(kBlue);
  out.str("");
  out << "sigarea";
  tree->Draw("de:mbc",out.str().c_str(),"same");

  tree->SetMarkerColor(kRed);
  out.str("");
  out << "(mbc>5.23 && mbc<5.26 || mbc>5.25 && de>0.12)";
  tree->Draw("de:mbc",out.str().c_str(),"same");
  ellican->Draw();
}

