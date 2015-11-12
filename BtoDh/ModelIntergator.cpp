using namespace std;

const double M_D = 1.86483;
const double M_D_sq = M_D*M_D;
const double M_Ks = 0.497614;
const double M_Ks_sq = M_Ks*M_Ks;
const double M_pi = 0.13957018;
const double M_pi_sq = M_pi*M_pi;
const double M_min = (M_pi+M_Ks)*(M_pi+M_Ks);
const double M_max = (M_D-M_pi)*(M_D-M_pi);

TH2F* hist;

int GetBin(const double& mp,const double& mm){
  int gbin = hist->FindBin(mp,mm);
  int dpbin = (int)floor(hist->GetBinContent(gbin)+0.01);

  return dpbin;
}

double CalcPhase(const double& Re, const double& Im){
  if(Re>0) return atan(Im/Re);
  if(Re<0 && Im>0) return atan(Im/Re) + TMath::Pi();
  if(Re<0 && Im<0) return atan(Im/Re) - TMath::Pi();
  if(Im == 0 && Re != 0) return TMath::Pi();
  if(Re == 0 && Im>0) return 0.5*TMath::Pi();
  if(Re == 0 && Im<0) return -0.5*TMath::Pi();
  return -999;
}

int CalcDelta(const double& ampRe, const double& ampIm, const double& ampbRe, const double& ampbIm, double& ampMod, double& ampbMod, double& phase){
  ampMod = ampRe*ampRe + ampIm*ampIm;
  ampbMod = ampbRe*ampbRe + ampbIm*ampbIm;
  double delt1 = CalcPhase(ampRe,ampIm);
  double delt2 = CalcPhase(ampbRe,ampbIm);
  phase = delt2 - delt1;
  return 0;
}

double E_Ks_star(const double& mp){return (mp-M_pi_sq+M_Ks_sq)*0.5/sqrt(mp);}
double E_pi_star(const double& mp){return (M_D_sq-mp-M_pi_sq)*0.5/sqrt(mp);}

double mm_max(const double& mp){
  if(mp > M_max || mp < M_min) return -1;
  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
  const double sum = EKs+Epi;
  const double PKs2 = EKs*EKs-M_Ks_sq;
  const double Ppi2 = Epi*Epi-M_pi_sq;
  const double dif = sqrt(PKs2)-sqrt(Ppi2);
  return (sum-dif)*(sum+dif);
}

double mm_min(const double& mp){
  if(mp > M_max || mp < M_min) return -1;
  const double EKs = E_Ks_star(mp), Epi = E_pi_star(mp);
  const double sum = EKs+Epi;
  const double PKs2 = EKs*EKs-M_Ks_sq;
  const double Ppi2 = Epi*Epi-M_pi_sq;
  const double sumP = sqrt(PKs2)+sqrt(Ppi2);
  return (sum-sumP)*(sum+sumP);
}

bool IsInPlot(const double& mp,const double& mm){
  if(mm >= mm_min(mp) && mm <= mm_max(mp)) return true;
  else return false;
}

void changeAB(double& A, double& B){
    A = A + B;
    B = A - B;
    A = A - B;
    return;
}

void ModelIntergator(void){
    TFile* binning = TFile::Open("../Tuples/dkpp_belle_ddd.root");
    hist = (TH2F*)binning->Get("dkpp_bin_h");
    ifstream table("../Tuples/amp_table.txt");
    string line;
    double mp,mm,ampRe,ampIm,ampbRe,ampbIm;
    double intC[8], intS[8], intK[8], intKb[8];
    double delta, abssq_amp, abssq_ampb;
    double Majorant = 0;
    int bin;
    for(int i=0; i<8; i++){
        intC[i] = 0;
        intS[i] = 0;
        intK[i] = 0;
        intKb[i] = 0;
    }
    while(!table.eof()){
        getline(table,line);
        sscanf(line.c_str(),"%lf %lf: amp(%lf,%lf), ampb(%lf,%lf)",&mp,&mm,&ampRe,&ampIm,&ampbRe,&ampbIm);
        if(!IsInPlot(mp,mm)) continue;
        bin = GetBin(mp,mm);
        CalcDelta(ampRe,ampIm,ampbRe,ampbIm,abssq_amp,abssq_ampb,delta);
        if(abssq_amp+abssq_ampb > Majorant) Majorant = abssq_amp+abssq_ampb;
        if(!bin || abs(bin)>8) continue;
        if(bin < 0){
            changeAB(abssq_amp,abssq_ampb);
            delta = -delta;
        }
        bin = abs(bin) - 1;
        intC[bin]  += sqrt(abssq_amp*abssq_ampb)*TMath::Cos(delta);
        intS[bin]  += sqrt(abssq_amp*abssq_ampb)*TMath::Sin(delta);
        intK[bin]  += abssq_amp;
        intKb[bin] += abssq_ampb;
    }

    double norm = 0;
    for(int i=0; i<8; i++){
        intC[i] /= sqrt(intK[i]*intKb[i]);
        intS[i] /= sqrt(intK[i]*intKb[i]);
        norm += intK[i] + intKb[i];
    }

    cout << "Norm     = " << norm << endl;
    cout << "Majorant = " << Majorant << endl;
    for(int i=0; i<8; i++){
        intK[i]  /= norm;
        intKb[i] /= norm;
        cout << i+1 << ": C = " << intC[i] << ", S = " << intS[i] << ", K = " << intK[i] << ", Kbar = " << intKb[i] << ", Q = " << intC[i]*intC[i] +intS[i]*intS[i] << endl;
    }
    return;
}
