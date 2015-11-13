void DrawComparison(const int N, const int ref[2][16], double arrs[N][2][16], double errs[N][2][16], const string& label){
  TMultiGraph* mg[2];
  const int Nb = 16;
  double bvals[Nb] = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
  double berrs[Nb] = { 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0};

  int colors[N], shapes[N];
  for(int i=0; i<N; i++){
    colors[i] = kBlue+10*i;
    shapes[i] = 20+i;
    for(int k=0; k<2; k++){
      for(int j=0; j<16; j++){
        arrs[i][k][j] -= ref[k][15-j];
      }
    }
  }

  stringstream out;
  out.str("");
  out << "c1_" << label;
  TCanvas* c1 = new TCanvas(out.str().c_str(),out.str().c_str(),800,800);
  c1->Draw();
  TPad* pads[2];
  TGraphErrors* grphs[2][N];
  for(int k=0; k<2; k++){
    c1->cd();
    mg[k] = new TMultiGraph();
    out.str(""); out << "pad" << k+1;
    pads[k] = new TPad(out.str().c_str(),out.str().c_str(),0.01,0.5*k+0.01,0.99,0.5*(k+1)-0.01);
    pads[k]->Draw();
    pads[k]->cd();
    for(int i=0; i<N; i++){
      grphs[k][i] = new TGraphErrors(Nb,bvals,arrs[i][k],berrs,errs[i][k]);
      grphs[k][i]->SetMarkerSize(1.2);
      grphs[k][i]->SetMarkerColor(colors[i]);
      grphs[k][i]->SetMarkerStyle(shapes[i]);
      mg[k]->Add(grphs[k][i]);
    }
    mg[k]->Draw("ap");
  }
}

void AnalysePurity(void){
  const string prefix("../PurityFit/params/");
  string fname[5];
  fname[0] = prefix + string("Integrals2_m1_mh010_mc.txt");
  fname[1] = prefix + string("Integrals2_m1_mh010_singlefbb_mc.txt");
  fname[2] = prefix + string("Integrals2_m1_mh010_fixedshape_mc.txt");
  fname[3] = prefix + string("Integrals2_m1_mh010_fixedshape_singlefbb_mc.txt");
  fname[4] = prefix + string("Integrals_m1_mh010_mc.txt");
  string rfname = prefix + string("TrueNumbers_m1_mh010.txt");

  int bin, flv;
  double nsig[5][2][16], nsig_err[5][2][16];
  double ncmb[5][2][16], ncmb_err[5][2][16];
  double ncnt[5][2][16], ncnt_err[5][2][16];
  double nprt[5][2][16], nprt_err[5][2][16];
  double pur[5][2][16], pur_err[5][2][16];

  int rnsig[2][16];
  int rncmb[2][16];
  int rncnt[2][16];
  int rnprt[2][16];
  double rnpur[2][16];

  string line;
  ifstream rfile(rfname.c_str(),std::ifstream::in);
  for(int i=0; i<46; i++) getline(rfile,line);
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      getline(rfile,line);
      sscanf(line.c_str(),"bin %d: %d %d %d %d",&bin,&rnsig[k][j],&rncmb[k][j],&rncnt[k][j],&rnprt[k][j]);
    }
    getline(rfile,line);
  }

  ifstream ifiles[5];
  for(int i=0; i<4; i++){
    ifiles[i] = ifstream(fname[i].c_str(),std::ifstream::in);
    getline(ifiles[i],line);
    for(int k=0; k<2; k++){
      for(int j=0; j<16; j++){
        getline(ifiles[i],line);
        sscanf(line.c_str(),"%d %d %lf +- %lf %lf +- %lf %lf +- %lf %lf +- %lf %lf +- %lf",&flv,&bin,&nsig[i][k][j],&nsig_err[i][k][j],&ncmb[i][k][j],&ncmb_err[i][k][j],&ncnt[i][k][j],&ncnt_err[i][k][j],&nprt[i][k][j],&nprt_err[i][k][j],&pur[i][k][j],&pur_err[i][k][j]);
      }
    }
  }
  DrawComparison(4,rnsig,nsig,nsig_err,string("sig"));
  DrawComparison(4,rncnt,ncnt,ncnt_err,string("cnt"));
  return;
}
 
