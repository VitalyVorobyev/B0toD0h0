void dt_graph(const int npflag = 0){
  TFile file("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s1_full.root");
  TTree* tree = (TTree*)file.Get("TEvent");
  const int NTot = tree->GetEntries();
  int NCount = 0;
  int NCountNoNPTAG = 0;
  const double NtotB = 1./(double)NTot;
  const int NBins = 100;
  const double dtmin = -1.5;
  const double dtmax =  1.5;
  const double ddt = (dtmax-dtmin)/NBins;
  double time[NBins];
  double dthist[NBins],dthistNoNPTAG[NBins];
  double zsig,zasc;
  int nptag;
  tree->SetBranchAddress("z_sig",&zsig);
  tree->SetBranchAddress("z_asc",&zasc);
  tree->SetBranchAddress("nptag",&nptag);
  for(int j=0; j<NBins; j++){
    time[j] = -(dtmin + (j+0.5)*ddt);
    dthist[j] = 0;
    dthistNoNPTAG[j] = 0;
  }
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    NCount++;
    int bin = ((zsig-zasc)-dtmin)/ddt;
    if(bin<0 || bin>99) continue;
    dthist[bin] += 1;//NtotB;
    if(!nptag){
      NCountNoNPTAG++;
      dthistNoNPTAG[bin] += 1;
    }
  }
  for(int i=0; i<NBins; i++){
    dthist[i] /= NCount;
    dthistNoNPTAG[i] /= NCountNoNPTAG;
  }
  TGraph* grdata = new TGraph(NBins,time,dthist);
  grdata->SetMarkerStyle(21);
  grdata->SetMarkerColor(kBlue);
  grdata->Print();

  TGraph* grdataNoNPTAG = new TGraph(NBins,time,dthistNoNPTAG);
  grdataNoNPTAG->SetMarkerStyle(21);
  grdataNoNPTAG->SetMarkerColor(kRed);

  char* fname;
  if(npflag) fname = "dt_graph.txt";
  else       fname = "dt_graph_np.txt";
  TGraph* grsig = new TGraph(fname,"%lg %lg %*lg %*lg %*lg %*lg %*lg %*lg");
  grsig->SetMarkerStyle(21);
  grsig->SetMarkerColor(kBlack);
//  grsig->Print();
  TGraph* grb0 = new TGraph(fname,"%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg");
//  grb0->Print();
  grb0->SetMarkerStyle(21);
  grb0->SetMarkerColor(kBlue);
  TGraph* grb0b = new TGraph(fname,"%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
//  grb0b->Print();
  grb0b->SetMarkerStyle(21);
  grb0b->SetMarkerColor(kRed);

  TGraph* grEn = new TGraph(fname,"%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
//  grb0b->Print();
  grEn->SetMarkerStyle(21);
  grEn->SetMarkerColor(kRed);

  TGraph* grxEn = new TGraph(fname,"%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
//  grb0b->Print();
  grxEn->SetMarkerStyle(20);
  grxEn->SetMarkerColor(kRed);

  TGraph* grnp = new TGraph(fname,"%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
//  grb0b->Print();
  grnp->SetMarkerStyle(20);
  grnp->SetMarkerColor(kGreen);

  TGraph* grdet = new TGraph(fname,"%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
//  grdet->Print();
  grdet->SetMarkerStyle(20);
  grdet->SetMarkerColor(kRed);

  TCanvas* c1 = new TCanvas("delta t distr","delta t distr",700,500);
  c1->cd();
  c1->SetGrid();
  TMultiGraph* mg = new TMultiGraph();
//  mg->Add(grb0);
//  mg->Add(grb0b);
  mg->Add(grsig);
  if(!npflag)mg->Add(grdata);
  if(npflag) mg->Add(grdataNoNPTAG);
//  mg->Add(grEn);
//  mg->Add(grxEn);
//  mg->Add(grnp);
//  mg->Add(grdet);

  mg->Draw("AP");

  c1->Update();

  return;
}
