void RemoveTBranch(char* file, char* trname, char* brname){
  TFile f(file,"update");
  TTree *T = (TTree*)f.Get(trname);
  TBranch *b = T->GetBranch(brname);
  T->GetListOfBranches()->Remove(b);
  TLeaf* l = T->GetLeaf(brname);
  T->GetListOfLeaves()->Remove(l);
  T->Write();
  T->Print();
  return;
}
