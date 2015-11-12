
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

   const int numchains=3;

   double chain[numchains];
   double chainerr[numchains];

   for(int i=0; i<numchains; i++){
      chain[i]=i;
      chainerr[i]=0.5;
   }

   Double_t times[numchains] = {
      0.32528,
      2.48328,
      0.3505
   };

   Double_t l2timerms[numchains] = {
      0.305687,
      12.0619,
      0.189249
   };

   std::string names[numchains] = {
      "sequence_L2_muon_standalone_mu6l",
      "sequence_L2_mu6l",
      "sequence_L2_te650"
   };

   TGraphErrors *l2timevschain = new TGraphErrors(numchains,chain,times,chainerr,l2timerms);
   TAxis *ax = l2timevschain->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   l2timevschain->GetHistogram()->GetXaxis()->Set(3,x1,x2);

   for(Int_t k=0;k<numchains;k++){
      l2timevschain->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   l2timevschain->SetMarkerStyle(21);
   l2timevschain->Draw("AP");
}
