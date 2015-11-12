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

const string prefix("../Tuples/");

void CrossFeed(void){
  stringstream out;
  double FeedTable[7][7];
//  const string ellicut = GetElliCut(-0.1,0.0921343,5.27093,5.28789);
  const string sigcut("(b0f == 1 || b0f == 5 || b0f == 10 || rndm_pi0)");
  const string crfcut("!rndm_pi0 && d0f == 1 && b0f == -1 && abs(d0ch0) == 421 && (abs(d0ch1) == 511 || abs(d0ch1) == 423)");

  string strmodes[7] = {"Pi0_s8_m","ETA_s3_m","ETA_s3_m","OMEGA_s6_m","ETAP_s1_m","DST0_s1_m","DST0_s1_m"};
  int modes[7]       = {1,2,2,3,5,10,20};
  int h0modes[7]     = {10,10,20,20,10,10,10};
//  double bdt_cuts[7] = {0.27,0.15,0.18,0.18,0.56,0.60,0.68};
  double de_mins[7]  = {-0.1,-0.0994528,-0.0534644,-0.0565915,-0.0873898,-0.1,-0.0853615};
  double de_maxs[7]  = {0.0921343,0.0758825,0.0453713,0.0467748,0.0642064,0.0647985,0.0630641};
  double mbc_mins[7] = {5.27093,5.27147,5.27205,5.27198,5.27288,5.27146,5.27182};
  double mbc_maxs[7] = {5.28789,5.28763,5.28735,5.28743,5.28702,5.28769,5.28746};

  const string ptcut("pt_pip>0.05 && pt_pim>0.05");
  const string ptcut2("pt_pi1>0.1 && pt_pi2>0.1");
  const string helcut("abs(cos_hel)>0.2");
  const string ppi0cut("p_pi0_h0>0.2");
  const string addstr(" && ");
  string modecuts[7];
  modecuts[0] = ptcut + string(" && bdt>0.27");
  modecuts[1] = ptcut + string(" && bdt>0.15 && e_g1>0.15");;
  modecuts[2] = ptcut + addstr + ptcut2 + addstr + ppi0cut + string(" && bdt>0.18");
  modecuts[3] = ptcut + addstr + ptcut2 + addstr + ppi0cut + addstr + helcut + string(" && bdt>0.18");
  modecuts[4] = ptcut + string(" && lh0>0.56");
  modecuts[5] = ptcut + string(" && lh0>0.60");
  modecuts[6] = ptcut + string(" && lh0>0.68");

//  TChain* tree = new TChain("TEvent","TEvent");
//  string types[4] = {string("uds"),string("charm"),string("charged"),string("mixed")};
//  int streams[6] = {0,1,2,3,4,5};
//  for(int i=0; i<4; i++){
//    for(int j=0; j<6; j++){
//      out.str("");
//      out << prefix << "Fil_b2dh_" << types[i] << "_" << streams[j] << "_1" << streams[j] << "_m1_h0m10.root";
//      tree->Add(out.str().c_str());
//    }
//  }
//  double de,mbc,bdt,pt_pip,pt_pim,cos_hel,p_pi0_h0,pt_pi1,pt_pi2,lh0;
//  tree->SetBranchAddress("de",&de);
//  tree->SetBranchAddress("mbc",&mbc);
//  tree->SetBranchAddress("bdt",&bdt);
//  tree->SetBranchAddress("lh0",&lh0);
//  tree->SetBranchAddress("pt_pim",&pt_pim);
//  tree->SetBranchAddress("pt_pip",&pt_pip);
//  tree->SetBranchAddress("pt_pi1",&pt_pi1);
//  tree->SetBranchAddress("pt_pi2",&pt_pi2);
//  tree->SetBranchAddress("p_pi0_h0",&p_pi0_h0);
//  tree->SetBranchAddress("cos_hel",&cos_hel);

  TChain* tree[7][7];
  for(int i=0; i<7; i++){
    for(int j=0; j<7; j++){
//      cout << i << " " << j << endl;
      out.str("");
      out << prefix << "Fil_b2dh_sigmc" << strmodes[i];
      out << modes[j] << "_h0m" << h0modes[j] << ".root";
      tree[i][j] = new TChain("TEvent","TEvent");
      tree[i][j]->Add(out.str().c_str());
      out.str("");
      if(i == j){
        out << modecuts[j] << " && " << sigcut << " && " << GetElliCut(de_mins[j],de_maxs[j],mbc_mins[j],mbc_maxs[j]);
      } else{
        out << modecuts[j] << " && " << crfcut << " && " << GetElliCut(de_mins[j],de_maxs[j],mbc_mins[j],mbc_maxs[j]);
      }
//      cout << out.str() << endl;
      FeedTable[i][j] = tree[i][j]->GetEntries(out.str().c_str());
      cout << FeedTable[i][j] << " ";
    }
    cout << endl;
  }
  return;
}
