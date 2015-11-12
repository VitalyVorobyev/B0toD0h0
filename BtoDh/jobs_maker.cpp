//#include <vector>
//#include <fstream>

using namespace std;

const int NEXP = 27;
const int expv[NEXP] = {7,9,11,13,15,17,19,21,23,25,27,31,33,35,37,39,41,43,45,47,49,51,53,55,61,63,65};

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

void jobs_maker(const int mode = 0,const int SVD = 0, const int STREAM = 0){
  stringstream out;
  string line;
  ofstream *ofile;
//  vector< ofstream > ofilev;
//  vector<string> gens;
  const bool ntflag = true;
  ofstream* ofilev[4];
  string gens[4];
  if(!mode){
    if(!ntflag) ofile = new ofstream("scripts_data_skim.csh");
    else        ofile = new ofstream("scripts_data_nt.csh");
  } else if(mode == 1){// skim
    ofilev[0] = new ofstream("scripts_gen_uds_skim_2.csh");
    ofilev[1] = new ofstream("scripts_gen_charm_skim_2.csh");
    ofilev[2] = new ofstream("scripts_gen_charged_skim_2.csh");
    ofilev[3] = new ofstream("scripts_gen_mixed_skim_2.csh");

    gens[0] = string("uds");
    gens[1] = string("charm");
    gens[2] = string("charged");
    gens[3] = string("mixed");
  } else{// NTuple
    ofilev[0] = new ofstream("scripts_gen_uds_ntuple_2.csh");
    ofilev[1] = new ofstream("scripts_gen_charm_ntuple_2.csh");
    ofilev[2] = new ofstream("scripts_gen_charged_ntuple_2.csh");
    ofilev[3] = new ofstream("scripts_gen_mixed_ntuple_2.csh");

    gens[0] = string("uds");
    gens[1] = string("charm");
    gens[2] = string("charged");
    gens[3] = string("mixed");
  }

  const double max_lum = !mode ? 1200 : 2000;
  double cur_lum = 0;
  double lum;
  int min_run = 0;
  int max_run = 0;
  int cur_run;
  for(int i=0; i<NEXP; i++){
    min_run = 0; max_run = 0; cur_lum = 0;
    out.str("");
    out << "/home/vitaly/B0toDh0/lumi/lumexp" << exp(i,SVD) << ".txt";
    ifstream lumfile(out.str().c_str());
    if(!lumfile.is_open()){
      cout << "Can't open file " << out.str() << endl;
      continue;
    }
    while(getline(lumfile,line)){
      istringstream iss(line);
      iss >> cur_run; iss >> lum;
      if(cur_lum+lum < max_lum){
        cur_lum += lum;
      } else{
        min_run = max_run+1; max_run = cur_run-1;
        switch(mode){
        case 0:// data
          if(!ntflag){
            out.str("");
            out << "./DataSkim.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " > " << "scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh";
            *ofile << out.str() << endl;
            *ofile << "chmod 755 " << "./scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
            *ofile << "bsub -q b_a " << "./scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
          } else{
            out.str("");
            out << "./DataNT.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " > " << "scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh";
            *ofile << out.str() << endl;
            *ofile << "chmod 755 " << "./scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
            *ofile << "bsub -q l " << "./scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
          }
          break;
        case 1:// generic
          int stream = exp(i,SVD) > 30 ? STREAM : STREAM + 10;
          for(int j=0; j<4; j++){
            out.str("");
            out << "./GenSkim.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " " << stream << " " << gens[j] << " > scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh";
            *ofilev[j] << out.str() << endl;
            *ofilev[j] << "chmod 755 ./scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
            *ofilev[j] << "bsub -q b_a ./scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
          }
          break;
        case 2:// nt
          int stream = exp(i,SVD) > 30 ? STREAM : STREAM + 10;
          for(int j=0; j<4; j++){
            out.str("");
            out << "./GenNT.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " " << stream << " " << gens[j] << " > scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh";
            *ofilev[j] << out.str() << endl;
            *ofilev[j] << "chmod 755 ./scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
            *ofilev[j] << "bsub -q l ./scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
          }
          break;
        }
        cout << exp(i,SVD) << " " << min_run << " " << max_run << " " << cur_lum << endl;
        cur_lum = 0; cur_lum += lum;
      }
    }
    min_run = max_run+1; max_run = cur_run-1;
    switch(mode){
    case 0:// data
      if(!ntflag){
        out.str("");
        out << "./DataSkim.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " > " << "scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh";
        *ofile << out.str() << endl;
        *ofile << "chmod 755 " << "./scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
        *ofile << "bsub -q b_a " << "./scripts/Data_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
      } else{
        out.str("");
        out << "./DataNT.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " > " << "scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh";
        *ofile << out.str() << endl;
        *ofile << "chmod 755 " << "./scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
        *ofile << "bsub -q l " << "./scripts/DataNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << ".csh" << endl;
      }
      break;
    case 1:
      int stream = exp(i,SVD) > 30 ? STREAM : STREAM + 10;
      for(int j=0; j<4; j++){
        out.str("");
        out << "./GenSkim.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " " << stream << " " << gens[j] << " > scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh";
        *ofilev[j] << out.str() << endl;
        *ofilev[j] << "chmod 755 ./scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
        *ofilev[j] << "bsub -q b_a ./scripts/Gen_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
      }
      break;
    case 2:
      int stream = exp(i,SVD) > 30 ? STREAM : STREAM + 10;
      for(int j=0; j<4; j++){
        out.str("");
        out << "./GenNT.csh " << exp(i,SVD) << " " << min_run << " " << max_run << " " << stream << " " << gens[j] << " > scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh";
        *ofilev[j] << out.str() << endl;
        *ofilev[j] << "chmod 755 ./scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
        *ofilev[j] << "bsub -q l ./scripts/GenNT_" << exp(i,SVD) << "_" << min_run << "_" << max_run << "_" << stream << "_" << gens[j] << ".csh" << endl;
      }
    }
    cout << exp(i,SVD) << " " << min_run << " " << max_run << " " << cur_lum << endl;
    cur_lum = 0; cur_lum += lum;
  }
  return;
}
