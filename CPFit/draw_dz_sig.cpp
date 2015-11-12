void draw_dz_sig(void){
  TCanvas* c1 = new TCanvas("c1","c1",700,500);
  c1->cd();

  TGraph* grsig = new TGraph("dz_sig_rf.txt","%lg %lg");
  grsig->SetMarkerStyle(21);
  grsig->SetMarkerColor(kBlue);
//  grsig->Print();

  TFile *file = TFile::Open("../Tuples/fil_b2dh_sigmcOmega_s1new.root");
  TTree* tree = (TTree*)file->Get("TEventTr");
//  RooArgSet argset;
//  RooRealVar dz_mc_sig1("dz_mc_sig1","dz_mc_sig1",-0.5,0.5,"mm"); argset.add(dz_mc_sig1);
//  RooDataSet ds("ds","ds",tree,argset);
//  RooDataHist* dh = ds.binnedclone();
//  RooHistPdf* hpdf = new RooHistPdf("hpdf","hpdf",RooArgSet(dz_mc_sig1),*dh);
  tree->Draw("dz_mc_sig1","abs(dz_mc_sig1)<0.5 && sz_sig<0.4");
  TH1* hist = tree->GetHistogram();
  const Double_t scale = 1./hist->Integral();
//  const Double_t scale = hist->GetXaxis()->GetBinWidth(1)/hist->GetIntegral();
  hist->Scale(scale);

  ///////////
  // Plots //
  ///////////
  hist->Draw();
  grsig->Draw("LP");

  c1->Update();
  c1->Print("c1.png");
  c1->Print("c1.root");
 
  return;
}
