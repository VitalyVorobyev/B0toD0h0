#include "cuts.h"

void draw_f_bkg(const int _mode){
  TChain* tree = new TChain("TEvent","TEvent");
  switch(_mode){
  case 1: tree->Add("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m1_mh010.root"); break;
  case 2: tree->Add("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m2_mh010.root"); break;
  case 3: tree->Add("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m2_mh020.root"); break;
  case 4: tree->Add("/home/vitaly/B0toDh0/PurityFit/data/mixtree_m3_mh020.root"); break;
  }

  TCanvas* ellican = new TCanvas("ellican","ellican",1600,800);
  ellican->Divide(8,4);
  ellican->Draw();
  stringstream out;
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      const int Bin = bin(j);
      const int padnum = k*16+j;
      ellican->cd(padnum+1);
      out.str("");
      if(k == 0) out << "bin == " << Bin << " && flv==1";// << " && elli";
      else       out << "bin == " << Bin << " && flv==-1";// << " && elli";
      const string BF = out.str();
      tree->Draw("mbc:de",BF.c_str());
      for(int i=0; i<10; i++){
        out.str("");
        out << BF << " && f_bkg>" << i*0.1 << " && f_bkg<" << (i+1)*0.1;// << " && elli";
        tree->SetMarkerStyle(8);
        tree->SetMarkerSize(0.2);
        tree->SetMarkerColor(i+1);
        tree->Draw("mbc:de",out.str().c_str(),"same");
//        ellican->Update();
      }
    }
  }
  ellican->cd();
  ellican->Print("pics/f_bkg_m2.pdf");
  return;
}
