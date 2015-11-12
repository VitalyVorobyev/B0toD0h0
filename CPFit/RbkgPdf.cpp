#include "RbkgPdf.h"
#include <fstream>

int RbkgPdf::GetParametersFromFile(const int mode, const int h0mode, const int svd, const bool mc, const int type_flag,const int ppp_flag){
  stringstream out;
  out.str("");
  out << "params/BackParams_m" << mode << "_mh0" << h0mode << "_svd" << svd;
  if(mc) out << "_mc";
  if(type_flag == 1)      out << "_cont";
  else if(type_flag == 2) out << "_BB";
  if(ppp_flag == 1) out << "_ppp";
  else if(ppp_flag == 2) out << "_gg";
  out << ".txt";
  ifstream ifile;
  ifile.open(out.str().c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "Can't open file " << out.str() << endl;
    return -1;
  } else{
    cout << "Getting background description from file " << out.str() << endl;
  }
  string postfix;
  if(svd == 1) postfix = string("_svd1");
  else         postfix = string("_svd2");
  string line,name;
  double val,err;
  char namech[15];
  int counter = 0;
  for(int i=0; i<13; i++){
    getline(ifile,line);
    sscanf(line.c_str(),"%s = %lf +- %lf",namech,&val,&err);
    name = string(namech);
    cout << name << " " << val << " +- " << err << endl;
    if(name == string("tau"))                { tau = val;         counter++; continue;}
    if(name == string("mu")+postfix)         { mu = val;          counter++; continue;}
    if(name == string("mu_delta")+postfix)   { mu_delta = val;    counter++; continue;}
    if(name == string("f_delta_mlt")+postfix){ f_delta_mlt = val; counter++; continue;}
    if(name == string("f_tail_mlt")+postfix) { f_tail_mlt = val;  counter++; continue;}
    if(name == string("S_main_mlt")+postfix) { S_main_mlt = val;  counter++; continue;}
    if(name == string("S_tail_mlt")+postfix) { S_tail_mlt = val;  counter++; continue;}
    if(name == string("f_delta_sgl")+postfix){ f_delta_sgl = val; counter++; continue;}
    if(name == string("f_tail_sgl")+postfix) { f_tail_sgl = val;  counter++; continue;}
    if(name == string("S_main_sgl")+postfix) { S_main_sgl = val;  counter++; continue;}
    if(name == string("S_tail_sgl")+postfix) { S_tail_sgl = val;  counter++; continue;}
    if(name == string("f_otlr")+postfix)     { f_otlr = val;      counter++; continue;}
    if(name == string("s_otlr")+postfix)     { s_otlr = val;      counter++; continue;}
  }
  cout << "GetParametersFromFile: " << counter << endl;
  return counter;
}

RbkgPdf::RbkgPdf(const int mode, const int svd, const bool mc, const bool tune): cnvl() {
  cm2ps = 78.48566945838871754705;
  cout << "RbkgPdf mode " << mode << endl;
  f_otlr = 0.;
  s_otlr = 30.;
  tau = 1.20152;
  switch(mode){
  case 1:// pi0
    cout << "Init RbkgPdf for pi0" << endl;
    if(svd == 2){
      tau         = 0.780211;// +- 0.0901569
      mu          =-0.166954;// +- 0.0534155
      mu_delta    = 0.0285984;// +- 0.0248846
      f_delta_mlt = tune ? 0.166906 : 0.572719;// +- 0.0640259
      f_tail_mlt  = 0.0738474;// +- 0.0164615
      S_main_mlt  = 1.11188;// +- 0.0354149
      S_tail_mlt  = 4.26611;// +- 0.475523
      f_delta_sgl = tune ? 0.0234418 : 0.479424;// +- 0.101054
      f_tail_sgl  = 0.0558646;// +- 0.0116118
      S_main_sgl  = 1.11764;// +- 0.0460215
      S_tail_sgl  = 11.8508;// +- 2.91068
      f_otlr      = 0;//tune ? 0.00779259 : 0.0047021;// +- 0.00190822
      s_otlr      = 25.8435;
    } else{
      tau         = 0.79461;
      mu          =-0.018311;
      mu_delta    =-0.00896778;// +- 0.0385935
      f_delta_mlt = tune ? 0 : 1;// +- 0.0932838
      f_tail_mlt  = 0.0982681;// +- 0.0322716
      S_main_mlt  = 1.24241;// +- 0.0571264
      S_tail_mlt  = 3.33009;// +- 0.544847
      f_delta_sgl = tune ? 0 : 1;// +- 0.174194
      f_tail_sgl  = 0.107781;// +- 0.050587
      S_main_sgl  = 1.24782;// +- 0.0859298
      S_tail_sgl  = 3.36498;// +- 1.14478
      f_otlr      = 0;//tune ? 0.0182992 : 0.0189294;// +- 0.0071703
      s_otlr      = 41.2166;
    }
    break;
  case 2:// eta -> gg
    if(svd == 2){
      tau         = 2.30926;// +- 0.315209
      mu          =-0.598281;// +- 0.304824
      mu_delta    =-0.0871882;// +- 0.0175144
      f_delta_mlt = 0.214413;// +- 0.0357144
      f_tail_mlt  = 0.141492;// +- 0.0137353
      S_main_mlt  = 1.26948;// +- 0.0183759
      S_tail_mlt  = 3.36794;// +- 0.139247
      f_delta_sgl = 0.770724;// +- 0.036027
      f_tail_sgl  = 0.0897614;// +- 0.00979697
      S_main_sgl  = 1.14049;// +- 0.0251209
      S_tail_sgl  = 5.47538;// +- 0.419777
      f_otlr      = 0.0125753;// +- 0.00123563
      s_otlr      = 31.5617;// +- 2.49151
    } else{
      tau         = 0.79461;
      mu          =-0.400879;
      mu_delta    =-0.0111621;// +- 0.0234151;
      f_delta_mlt = 0.999606;// +- 0.00390683;
      f_tail_mlt  = 0.266935;// +- 0.0445861;
      S_main_mlt  = 1.22732;// +- 0.0566526;
      S_tail_mlt  = 2.58338;// +- 0.159928;
      f_delta_sgl = 0.86151;// +- 0.14076
      f_tail_sgl  = 0.110071;// +- 0.0458817
      S_main_sgl  = 1.15299;// +- 0.0660206
      S_tail_sgl  = 3.55447;// +- 0.825178
      f_otlr      = 0.0252212;// +- 0.00390868;
      s_otlr      = 32.7912;// +- 4.63905
    }
    break;
  case 3:// eta -> pi+ pi- pi0
    if(svd == 2){
      mu          = 0.00245452;// +- 0.0390824
      mu_delta    = 0.0259721;// +- 0.00823168
      f_delta_mlt = 0.828073;// +- 0.0177527
      f_tail_mlt  = 0.372017;// +- 0.0384497
      S_main_mlt  = 1.18743;// +- 0.0452916
      S_tail_mlt  = 2.30686;// +- 0.0731427
      f_delta_sgl = 0.448808;// +- 0.0275643
      f_tail_sgl  = 0;// +- 0.0183342
      S_main_sgl  = 0.126242;// +- 0.034274
      S_tail_sgl  = 8.38917;// +- 2.12927
      f_otlr      = 0.0135348;// +- 0.00122732
      s_otlr      = 17.979;// +- 0.982226
    } else{
      mu          = 0.0582829;// +- 0.131984
      mu_delta    =-0.0356215;// +- 0.0382544
      f_delta_mlt = 0.748569;// +- 0.0766957
      f_tail_mlt  = 0.318554;// +- 0.0955359
      S_main_mlt  = 1.10097;// +- 0.115359
      S_tail_mlt  = 2.56057;// +- 0.223232
      f_delta_sgl = 0.508574;// +- 0.102941
      f_tail_sgl  = 0.070058;// +- 0.0152614
      S_main_sgl  = 0.110032;// +- 0.00997487
      S_tail_sgl  = 9.95;// +- 1.89514
      f_otlr      = 0.00887226;// +- 0.00306883
      s_otlr      = 22.8887;// +- 5.31512
    }
    break;
  case 4:// omega
    if(svd == 2){
      tau         = 0.72333;// +- 0.0148013
      mu          = 0.0543665;// +- 0.00751679
      mu_delta    = 0.00281958;// +- 0.00405444
      f_delta_mlt = 0.514367;// +- 0.013092
      f_tail_mlt  = 0.090121;// +- 0.00451183
      S_main_mlt  = 1.26818;// +- 0.010527
      S_tail_mlt  = 4.47635;// +- 0.106123
      f_delta_sgl = 0.370363;// +- 0.0149424
      f_tail_sgl  = 0.0757399;// +- 0.00342845
      S_main_sgl  = 1.0831;// +- 0.0103165
      S_tail_sgl  = 6.81437;// +- 0.248745
      f_otlr      = 0.00579139;// +- 0.000355526
      s_otlr      = 23.3021;// +- 0.845739
    } else{
      tau         = 1.07632;// +- 0.0626816
      mu          = 0.0886096;// +- 0.0311117
      mu_delta    =-0.00802693;// +- 0.0126699
      f_delta_mlt = 0.621347;// +- 0.0292674
      f_tail_mlt  = 0.0590122;// +- 0.00920791
      S_main_mlt  = 1.27635;// +- 0.0259536
      S_tail_mlt  = 4.79342;// +- 0.199744
      f_delta_sgl = 0.487717;// +- 0.0455091
      f_tail_sgl  = 0.0471682;// +- 0.00434327
      S_main_sgl  = 1.15065;// +- 0.0267731
      S_tail_sgl  = 15;// +- 0.236036
      f_otlr      = 0.00484446;// +- 0.000891845
      s_otlr      = 30.1879;// +- 3.90786
    }
    break;
  case 55:// B+ -> D0 pi+ sideband
    f_delta_sgl = 1.;
    f_delta_mlt = 1.;
    mu          =-0.174741;// is not used
    tau         = 1.24768; // is not used
    if(svd == 2){
      if(!mc){// SVD2, Data
        f_delta_sgl = 0.553644;// +- 0.0722156
        f_delta_mlt = 0.731629;// +- 0.0526017
        mu          =-0.174741;// +- 0.0745981
        tau         = 1.24768;// +- 0.105606
        f_tail_mlt  = 0.232202;// +- 0.0739445
        S_main_mlt  = 1.27208;// +- 0.060037
        S_tail_mlt  = 2.21991;// +- 0.195238
        f_tail_sgl  = 0.186818;// +- 0.177549
        S_main_sgl  = 1.06507;// +- 0.137236
        S_tail_sgl  = 2.25845;// +- 0.72261
        mu_delta    = 0.00362535;// +- 0.0178782
      } else{// SVD2, Generic MC
        f_tail_mlt = 0.141636;// +- 0.00940239
        S_main_mlt = 1.49728;// +- 0.0211451
        S_tail_mlt = 3.3682;// +- 0.0841069
        f_tail_sgl = 0.839901;// +- 0.0156033
        S_main_sgl = 5.18187;// +- 0.311361
        S_tail_sgl = 0.251266;// +- 0.0132309
        mu_delta   =-0.0210223;// +- 0.00642974
      }
    } else{//SVD1 Data
        f_delta_sgl = 0.56811;// +- 0.180594
        f_delta_mlt = 0.651789;// +- 0.119589
        mu          =-0.287062;// +- 0.207576
        tau         = 1.24768;// fixed
        f_tail_mlt  = 0.0358332;// +- 0.0342834
        S_main_mlt  = 1.29706;// +- 0.115975
        S_tail_mlt  = 5.26935;// +- 2.00908
        f_tail_sgl  = 0.0768885;// +- 0.0637036
        S_main_sgl  = 1.08073;// +- 0.127863
        S_tail_sgl  = 4.1387;// +- 2.10512
        mu_delta    = -0.04701;// +- 0.0636687
    }
    break;
  case 66:// B+ single D0
    if(svd == 2){
      f_delta_sgl = 0.708899;// +- 0.0709401
      f_delta_mlt = 0.739085;// +- 0.0359394
      mu          = -0.248165;// +- 0.102662
      tau         = 1.24768;// (fixed)
      f_tail_mlt  = 0.205906;// +- 0.0440248
      S_main_mlt  = 1.16552;// +- 0.0428438
      S_tail_mlt  = 2.38594;// +- 0.163264
      f_tail_sgl  = 0.41778;// +- 0.120553
      S_main_sgl  = 0.906777;// +- 0.113909
      S_tail_sgl  = 2.27074;// +- 0.175665
      mu_delta    =-0.0086097;// +- 0.0226586
    } else{
      f_delta_sgl = 0.351872;// +- 0.242127
      f_delta_mlt = 0.817148;// +- 0.161986
      mu          =-0.14194;// +- 0.220547
      tau         = 1.24768;// (fixed)
      f_tail_mlt  = 0.11197;// +- 0.0555087
      S_main_mlt  = 1.24885;// +- 0.0997688
      S_tail_mlt  = 4.10095;// +- 0.87803
      f_tail_sgl  = 0.125882;// +- 0.0935304
      S_main_sgl  = 0.932641;// +- 0.150454
      S_tail_sgl  = 3.57926;// +- 1.24718
      mu_delta    =-0.0704859;// +- 0.0634029
    }
    break;
  case 101:// pi0 toybkg
    if(svd == 2){
      tau         = 0.67958;
      mu          = -0.011487;
      mu_delta    = 0.0573504;
      f_delta_mlt = 0.0127504;
      f_tail_mlt  = 0.114427;
      S_main_mlt  = 0.184271;
      S_tail_mlt  = 15;// +- 0.383029
      f_delta_sgl = 0.0102088;
      f_tail_sgl  = 0.217907;
      S_main_sgl  = 0.153081;
      S_tail_sgl  = 15;// +- 1.6698
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    } else{
      tau         = 0.924677;
      mu          = 0.013262;
      mu_delta    = -0.0735636;
      f_delta_mlt = 0.0449784;
      f_tail_mlt  = 0.0279165;
      S_main_mlt  = 0.357818;
      S_tail_mlt  = 15;// +- 0.383029
      f_delta_sgl = 0.0483864;
      f_tail_sgl  = 0.0737662;
      S_main_sgl  = 0.246649;
      S_tail_sgl  = 15;// +- 1.6698
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    }
    break;
  case 102:// etagg toybkg
    if(svd == 2){
      tau         =1.30318;
      mu          =-0.0144213;
      mu_delta    =-0.0480987;
      f_delta_mlt =0.000106607;
      f_tail_mlt  =0.547466;
      S_main_mlt  =0.00389399;
      S_tail_mlt  =0.0439527;
      f_delta_sgl =0.000540051;
      f_tail_sgl  =0.139334;
      S_main_sgl  =0.160383;
      S_tail_sgl  =0.158175;
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    } else{
      tau         = 1.29233;
      mu          =-0.0215312;
      mu_delta    = 0.0503878;
      f_delta_mlt = 0.00239699;
      f_tail_mlt  = 0.479131;
      S_main_mlt  = 0.108825;
      S_tail_mlt  = 1.26172;
      f_delta_sgl = 0.00353238;
      f_tail_sgl  = 0.0824128;
      S_main_sgl  = 0.025976;
      S_tail_sgl  = 5.5189;
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    }
    break;
  case 103:// etappp toybkg
    if(svd == 2){
      tau         = 1.22802;
      mu          =-0.00620251;
      mu_delta    = 0.00797301;
      f_delta_mlt = 0.00591479;
      f_tail_mlt  = 0.215427;
      S_main_mlt  = 0.190978;
      S_tail_mlt  = 4.29344;
      f_delta_sgl = 0.0108685;
      f_tail_sgl  = 0.148488;
      S_main_sgl  = 0.0581581;
      S_tail_sgl  = 6.93318;
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    } else{
      tau         = 1.35505;
      mu          = 0.0298821;
      mu_delta    =-0.0945215;
      f_delta_mlt = 7.20076e-09;
      f_tail_mlt  = 0.0228605;
      S_main_mlt  = 0.160168;
      S_tail_mlt  = 1.26804;
      f_delta_sgl = 0.0427675;
      f_tail_sgl  = 0.999957;
      S_main_sgl  = 0.0975763;
      S_tail_sgl  = 1.85798;
      f_otlr      = 0;
      s_otlr      = 41.2166;// +- 9.86021
    }
    break;
  case 104:// omega toybkg
    if(svd == 2){
      tau         = 1.22826;
      mu          =-0.00382201;
      mu_delta    = -0.0282578;
      f_delta_mlt = 0.0038321;
      f_tail_mlt  = 0.219835;
      S_main_mlt  = 0.161302;
      S_tail_mlt  = 4.88746;
      f_delta_sgl = 0.000517077;
      f_tail_sgl  = 1.73052e-06;
      S_main_sgl  = 1.74157e-05;
      S_tail_sgl  = 0.00848174;
      f_otlr      = 0;
      s_otlr      = 25.8435;// +- 1.23107
    } else{
      tau         = 1.35509;
      mu          = 0.0298602;
      mu_delta    =-0.0947775;
      f_delta_mlt = 8.04556e-08;
      f_tail_mlt  = 0.000504229;
      S_main_mlt  = 0.160541;
      S_tail_mlt  = 0.090424;
      f_delta_sgl = 0.0427621;
      f_tail_sgl  = 5.52021e-08;
      S_main_sgl  = 0.180948;
      S_tail_sgl  = 4.80727;
      f_otlr      = 0;
      s_otlr      = 41.2166;// +- 9.86021
    }
    break;
  default:
    mu_delta    =-0.000501673;
    mu          =-0.0129797;
    tau         =1.20152;
    f_tail_mlt = 0.318409;
    S_main_mlt = 1.14702;
    S_tail_mlt = 1.94894;
    f_tail_sgl = 0.0410285;
    S_main_sgl = 1.11925;
    S_tail_sgl = 0.147687;
    f_delta_mlt = 0.805273;
    f_delta_sgl = 0.58846;
    break;
  }
  return;
}

double RbkgPdf::Pdf(const double &x, const double &s, const int ndf){
  const double sigma_main = ndf ? s*S_main_mlt*cm2ps    : s*S_main_sgl*cm2ps;
  const double sigma_tail = ndf ? sigma_main*S_tail_mlt : sigma_main*S_tail_sgl;
  const double f_tail     = ndf ? f_tail_mlt            : f_tail_sgl;
  const double f_delta    = ndf ? f_delta_mlt           : f_delta_sgl;

  const double pdf_l = (1-f_tail)*Enp_conv_gauss(x,tau,tau,mu,sigma_main)+f_tail*Enp_conv_gauss(x,tau,tau,mu,sigma_tail);
  const double int_pdf_l = (1-f_tail)*norm_Enp_conv_gauss(ll,ul,tau,tau,mu,sigma_main)+f_tail*norm_Enp_conv_gauss(ll,ul,tau,tau,mu,sigma_tail);

  const double pdf_d = (1-f_tail)*gaussian(x,mu_delta,sigma_main)+f_tail*gaussian(x,mu_delta,sigma_tail);
  const double int_pdf_d = (1-f_tail)*norm_gaussian(ll,ul,mu_delta,sigma_main)+f_tail*norm_gaussian(ll,ul,mu_delta,sigma_tail);

  const double pdf = f_delta*pdf_d + (1-f_delta)*pdf_l;
  const double int_pdf = f_delta*int_pdf_d + (1-f_delta)*int_pdf_l;
  if(pdf>0 && int_pdf>0){
    if(f_otlr>0.0001){ return AddOutlier(x,pdf,int_pdf);}
    else{              return pdf/int_pdf;}
  } else{
//    cout << "RbkgPdf::Pdf: " << pdf << ", norm = " << int_pdf << endl;
    return 0;
  }
}

double RbkgPdf::AddOutlier(const double& x, const double Lin,const double& nLi){
  const double Lol = gaussian(x,0.,s_otlr);
  const double nLol = norm_gaussian(ll,ul,0.,s_otlr);
  const double Li = (1.0-f_otlr)*Lin/nLi + f_otlr*Lol/nLol;
  return Li;
}
