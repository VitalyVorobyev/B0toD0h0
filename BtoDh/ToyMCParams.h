const double _tau = 1.519;// +- 0.007 ps
const double _dm  = 0.507;// +- 0.004 ps^{-1}
const double _sin2beta = 0.682;// +- 0.019
const double _cos2beta = TMath::Sqrt(1.-_sin2beta*_sin2beta);
const double _sigma_over_tau = 0.8;
const double _phi = 0.75049;

const double mistag_rate = 0.3;
const double _purity = 0.6;

const bool constK      = true;
const bool constMistag = true;
const bool constFSig   = true;
const bool constSigma  = true;
const bool constBeta   = false;
const bool constPsi    = false;

const bool softlimit = true;

const bool PoissonFlag = false;
const char* fname = "/home/vitaly/B0toDh0/BtoDh/bdtg_1000_5000.txt";
//const char* fname = "/home/vitaly/B0toDh0/BtoDh/bdtg_500_3000.txt";
//const char* fnamefit = "/home/vitaly/B0toDh0/BtoDh/bdtg_500_3000.txt";
const char* fnamefit = "/home/vitaly/B0toDh0/BtoDh/bdtg_1000_5000.txt";

const float _C[8] = {
0.365,
0.710,
0.481,
0.008,
-0.757,
-0.884,
-0.462,
0.106
};

const float _S[8] = {
-0.179,
-0.013,
-0.147,
0.938,
0.386,
-0.162,
-0.616,
-1.063
};

const float K8[8] = {
0.067/(0.016+0.067),
0.093/(0.027+0.093),
0.030/(0.013+0.030),
0.089/(0.040+0.089),
0.079/(0.015+0.079),
0.102/(0.012+0.102),
0.123/(0.026+0.123),
0.159/(0.109+0.159)
};

const float _K[8] = {
0.067,
0.093,
0.030,
0.089,
0.079,
0.102,
0.123,
0.159
};

const float _Kb[8] = { //{0,0,0,0,0,0,0,0};
0.016,
0.027,
0.013,
0.040,
0.015,
0.012,
0.026,
0.109
};

