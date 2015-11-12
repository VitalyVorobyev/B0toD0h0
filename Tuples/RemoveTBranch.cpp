void RemoveTBranch(char* file, char* trname, char* brname){
  TFile f(file,"update");
  TTree *T = (TTree*)f.Get(trname);
  TBranch *b = T->GetBranch(brname);
  T->GetListOfBranches()->Remove(b);
  T->Write();
  T->Print();
  return;
}
