const int expv[27] = {7,9,11,13,15,17,19,21,23,25,27,31,33,35,37,39,41,43,45,47,49,51,53,55,61,63,65};
const int expv_svd1[11] = {7,9,11,13,15,17,19,21,23,25,27};
const int expv_svd2[16] = {31,33,35,37,39,41,43,45,47,49,51,53,55,61,63,65};

int exp(const int bin,const int SVD){
  switch(SVD){
  case 0:
    if(bin<0 || bin>26) return 0;
    return expv[bin];
  case 1:
    if(bin<0 || bin>10) return 0;
    return expv_svd1[bin];
  case 2:
    if(bin<0 || bin>15) return 0;
    return expv_svd2[bin];
  }
  return -1;
}

int binexp(const int exp,const int SVD){
  switch(SVD){
  case 0:
    for(int i=0; i<27; i++){if(exp == expv[i]) return i;}
    break;
  case 1:
    for(int i=0; i<11; i++){if(exp == expv_svd1[i]) return i;}
    break;
  case 2: 
    for(int i=0; i<16; i++){if(exp == expv_svd2[i]) return i;}
    break;
  }
  return -1;
}

int lumi_test(char* file,const int SVD = 0){
  int exp_num;
  double total_on_res;
  switch(SVD){
  case 0:
    exp_num = 27;
    total_on_res = 715000;
    break;
  case 1:
    exp_num = 11;
    total_on_res = 140000;
    break;
  case 2:
    exp_num = 16;
    total_on_res = 571000;
    break;
  default:
    return -1;
  }
  vector<double> lum_v[exp_num];
  vector<int> use_v[exp_num];

  double total_lum = 0;
  double used_lum = 0;

  stringstream out;
  string line;
  double lum;
  int run;
  for(int i=0; i<exp_num; i++){
    cout << "Exp " << exp(i,SVD) << ":" << endl;
    out.str("");
    out << "/home/vitaly/B0toDh0/lumi/lumexp" << exp(i,SVD) << ".txt";
    ifstream lumfile(out.str().c_str());
    if(!lumfile.is_open()){
      cout << "Can't open file " << out.str() << endl;
      continue;
    }
    while(getline(lumfile,line)){
      istringstream iss(line);
      iss >> run; iss >> lum;
//      cout << "  run: " << run << ", lum: " << lum << endl;
      lum_v[i].push_back(lum);
      use_v[i].push_back(0);
      total_lum += lum;
    }
  }

  TChain* tree = new TChain("TEvent");
  tree->Add(file);
  const int NTot = tree->GetEntries();
  int exp;
  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  for(int i=0; i<NTot; i++){
    tree->GetEvent(i);
    const int bin = binexp(exp,SVD);
    if(bin<0){
      cout << "Negative bin!!!" << endl;
      continue;
    }
    if(!use_v[bin][run]){
      use_v[bin][run] = 1;
      used_lum += lum_v[bin][run];
    }
  }
//  cout << "Test completed:" << endl;
//  cout << "Missing runs:" << endl;
//  for(int i=0; i<exp_num; i++){
//    cout << "Exp " << exp(i,SVD) << ":" << endl;
//    for(int j=0; j<lum_v[i].size(); j++){
//      if(!use_v[i][j] && lum_v[i][j]>0) cout << " " << j;
//    }
//    cout << endl;
//  }
//  cout << "Used runs:" << endl;
//  for(int i=0; i<exp_num; i++){
//    cout << "Exp " << exp(i,SVD) << ":" << endl;
//    for(int j=0; j<lum_v[i].size(); j++){
//      if(use_v[i][j] && lum_v[i][j]>0) cout << " " << j;
//    }
//    cout << endl;
//  }
//  cout << "Used lum / Total lum = " << used_lum << " / " << total_lum << " = " << used_lum/total_lum << endl;
  cout << "Used lum / Total lum = " << used_lum << " / " << total_on_res << " = " << used_lum/total_on_res << endl;
  return 0;
}

