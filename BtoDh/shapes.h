// Combinatorics 2d fit
FCN=-42519.4 FROM HESSE     STATUS=OK             23 CALLS         122 TOTAL
                     EDM=0.000101476    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  argedge      5.28941e+00   2.16695e-05   7.66376e-05   8.67961e-01
   2  argpar      -2.37807e+01   9.77140e-01   6.01399e-04  -3.87711e-01
   3  c1          -4.57562e-01   1.42523e-02   4.05230e-05  -4.57722e-02
   4  c2           6.94555e-02   1.38983e-02   3.94823e-05   6.94561e-03

// Rho 2d fit
FCN=-5134.26 FROM HESSE     STATUS=OK             61 CALLS         409 TOTAL
                     EDM=0.000191754    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  cond         6.46863e-02   7.87968e-03   3.59929e-07   1.50965e+00
   2  de0r        -1.25915e-01   1.38529e-02   2.94686e-04  -5.66837e-01
   3  mbc0         5.28108e+00   3.55011e-04   8.04218e-05   5.42222e-02
   4  sl           1.52139e-02   9.99278e-04   2.23692e-05  -1.22013e+00
   5  slopel      -8.36595e+01   8.83161e+02   2.22665e-04   7.28932e-01
   6  sloper       1.39347e+02   9.29988e+02   2.63707e-04  -8.05685e-01
   7  sr           1.14841e-02   1.20591e-03   1.45184e-05  -1.26652e+00
   8  steep        4.15894e+00   3.07048e+00   9.22306e-05  -1.16005e+00

// Signal
FCN=-748601 FROM HESSE     STATUS=OK            250 CALLS        6467 TOTAL
                     EDM=0.00839972    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  alphal       3.85675e-01   1.34143e-01   5.09654e-05   3.85770e-02
   2  alphar      -1.56269e+00   1.88367e-01   4.27062e-05  -1.56912e-01
   3  de0          1.51293e-02   2.26574e-03   2.03775e-04   1.51877e-01
   4  deCBl       -2.58588e-02   4.89473e-03   3.65526e-04  -2.61560e-01
   5  deCBr       -1.22991e-02   5.02117e-03   6.52075e-04  -1.23303e-01
   6  fCBl         4.62725e-01   1.08468e-01   4.58604e-04  -7.46192e-02
   7  fCBr         2.11382e-01   6.59724e-02   8.34192e-04  -6.15341e-01
   8  fmbc         2.80705e-01   2.24472e-02   1.87540e-04  -4.54030e-01
   9  mbc0         5.28291e+00   1.71261e-04   2.67240e-05   1.46044e-01
  10  mbc00        5.27979e+00   4.62653e-05   1.32308e-05  -1.06595e-02
  11  nl           1.71068e+01   5.79516e+01   3.07191e-03  -7.17980e-01
  12  nr           3.79904e+00   7.79509e-01   2.00482e-04  -1.17846e+00
  13  s1           2.24893e-02   6.95252e-04   3.09875e-05  -1.14339e+00
  14  sCBl         2.71645e-02   4.38709e-03   2.67755e-04  -1.10030e+00
  15  sCBr         3.97215e-02   2.61255e-03   3.50630e-04  -9.99340e-01
  16  sl           6.28487e-03   1.49422e-04   1.27334e-05  -1.34609e+00
  17  sll          2.96898e-03   4.70941e-05   7.27391e-06  -1.41653e+00
  18  sr           1.92172e-03   5.99257e-05   1.10586e-05  -1.44673e+00
  19  srr          2.17430e-03   4.88485e-05   7.93512e-06  -1.43881e+00

  
// Generic MC
  Table b0f : ds(de<0.07 && de>-0.08 && mbc>5.272 && mbc<5.288)
  +---------+-----+
  |  signal | 776 |
  |     fsr |   8 |
  | bad_pi0 | 134 |
  |     rho |  34 |
  |    comb | 614 |
  +---------+-----+

  Table b0f : ds
  +---------+-------+
  |  signal |   863 |
  |     fsr |     8 |
  | bad_pi0 |   174 |
  |     rho |   898 |
  |    comb | 1.317e+04 |
  +---------+-------+
  
// 2d fit
FCN=-181050 FROM HESSE     STATUS=OK             50 CALLS         206 TOTAL
                     EDM=7.05268e-05    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  Ncomb        1.30854e+04   1.19380e+02   4.25875e-05  -8.56954e-01
   2  Nrho         9.49903e+02   4.07573e+01   8.41944e-04  -9.76071e-01
   3  Nsig         1.07986e+03   4.37099e+01   8.46282e-04  -9.30716e-01
   4  c1          -4.31413e-01   1.47290e-02   1.70660e-05  -4.31547e-02
   5  c2           4.40049e-02   1.45117e-02   1.67795e-05   4.40050e-03
   6  de0          1.54475e-02   2.57058e-03   3.02690e-04   1.55096e-01
   7  mbc0         5.28268e+00   4.71442e-04   2.77499e-04   1.34351e-01

// Rectangle
  +---------+-----+
  |  signal | 774 |
  |     fsr |   8 |
  | bad_pi0 | 133 |
  |     rho |  30 |
  |    comb | 573 |
  +---------+-----+

// Ellips
  +---------+-----+
  |  signal | 777 |
  |     fsr |   8 |
  | bad_pi0 | 133 |
  |     rho |  30 |
  |    comb | 549 |
  +---------+-----+

// Elli
  +---------+-----+
  |  signal | 750 |
  |     fsr |   8 |
  | bad_pi0 | 121 |
  |     rho |  21 |
  |    comb | 417 |
  +---------+-----+

// Full
  Table b0f : ds
  +---------+-------+
  |  signal |   863 |
  |     fsr |     8 |
  | bad_pi0 |   174 |
  |     rho |   898 |
  |    comb | 1.317e+04 |
  +---------+-------+
   
// 1d fit
  Table b0f : ds(de<0.07 && de>-0.08)
  +---------+-----+
  |  signal | 776 |
  |     fsr |   8 |
  | bad_pi0 | 134 |
  |     rho |  34 |
  |    comb | 614 |
  +---------+-----+

  Table b0f : ds
  +---------+------+
  |  signal |  845 |
  |     fsr |    8 |
  | bad_pi0 |  164 |
  |     rho |  803 |
  |    comb | 1588 |
  +---------+------+
  
Nsig = 911.6
Nrho = 47.14

FCN=-27998.4 FROM HESSE     STATUS=OK             31 CALLS         179 TOTAL
                     EDM=5.51465e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  Ncomb        1.50415e+03   9.04351e+01   9.73172e-05  -1.33312e+00
   2  Nrho         8.38113e+02   5.18867e+01   2.35546e-04   5.44745e-02
   3  Nsig         1.04981e+03   6.29243e+01   3.05881e-04   3.64914e-01
   4  c1          -4.67115e-01   7.60315e-02   2.20417e-05  -4.67285e-02
   5  exppar      -4.10853e+01   3.03064e+00   2.28651e-04   1.91362e-01

   // 1d fit
FCN=-28002 FROM HESSE     STATUS=OK             40 CALLS         232 TOTAL
                     EDM=4.176e-05    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  Ncomb        1.34720e+03   1.11176e+02   9.95935e-05  -1.34686e+00
   2  Nrho         9.21900e+02   5.84328e+01   2.45482e-04   1.75035e-01
   3  Nsig         1.12294e+03   7.34255e+01   3.24053e-04   4.79508e-01
   4  c1          -3.61248e-01   1.21386e-01   2.48322e-05  -3.61327e-02
   5  c2          -1.08085e-01   9.37308e-02   2.64929e-05  -1.08087e-02
   6  de0          1.71219e-02   2.64876e-03   5.88002e-04   1.72067e-01
   
Nsig = 974.4 +- 63.71
Nrho = 47.35 +- 3.001
Ncmb = 520.7 +- 42.97


FCN=-27998.6 FROM HESSE     STATUS=OK             16 CALLS          93 TOTAL
                     EDM=9.08594e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  Ncomb        1.45403e+03   6.02648e+01   1.98426e-05  -1.33743e+00
   2  Nrho         8.55362e+02   3.83751e+01   2.41121e-04   7.91714e-02
   3  Nsig         1.08260e+03   5.03166e+01   3.12935e-04   4.15565e-01
Nsig = 940.1 +- 43.69
Nrho = 43.93 +- 1.971
Ncmb = 553.7 +- 22.95


// Data fit
FCN=-57332.6 FROM HESSE     STATUS=OK             50 CALLS         207 TOTAL
                     EDM=7.19383e-05    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  Ncomb        4.61591e+03   7.09012e+01   1.24074e-04  -1.18620e+00
   2  Nrho         2.52074e+02   2.22824e+01   1.17553e-04  -1.32228e+00
   3  Nsig         3.91828e+02   2.63296e+01   1.01093e-04  -1.22571e+00
   4  c1          -4.12708e-01   2.48884e-02   1.62144e-05  -4.12825e-02
   5  c2           5.90759e-02   2.42794e-02   7.88434e-05   5.90762e-03
   6  de0          2.16901e-02   3.96157e-03   2.65119e-04   2.18638e-01
   7  mbc0         5.28232e+00   8.64084e-04   2.92586e-04   1.16321e-01
   
Nsig = 328.9 +- 22.1
Nrho = 7.434 +- 0.6571
Ncmb = 196.6 +- 3.019
Pury = 0.6172 +- 0.04147
Ellips:
Nsig = 327.3 +- 21.99
Nrho = 6.262 +- 0.5536
Ncmb = 187.6 +- 2.881
Pury = 0.6281 +- 0.04221
Elli:
Nsig = 314.1 +- 21.1
Nrho = 4.701 +- 0.4155
Ncmb = 155.3 +- 2.386
Pury = 0.6624 +- 0.04451

// Sidebands
// SB comb
FCN=-14551.9 FROM MIGRAD    STATUS=CONVERGED      39 CALLS          40 TOTAL
                     EDM=6.7464e-06    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  c1          -5.50812e-01   9.80517e-03   8.14900e-05   1.85944e-01
   2  c2           9.64572e-02   9.55882e-03   7.92930e-05  -2.69234e+00

// SR comb
FCN=-3119.53 FROM HESSE     STATUS=OK             10 CALLS          43 TOTAL
                     EDM=1.06421e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  c1          -7.41045e-01   2.02999e-02   1.56402e-05  -7.41724e-02
   2  c2           1.60880e-01   2.04495e-02   1.56963e-05   1.60887e-02

// SB data
FCN=-2429.7 FROM HESSE     STATUS=OK             10 CALLS          40 TOTAL
                     EDM=3.46393e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  c1          -4.90309e-01   2.42423e-02   1.64238e-05  -4.90505e-02
   2  c2           4.28379e-02   2.38680e-02   1.61549e-05   4.28381e-03
 
 
// mbc_rho_param = 2
Rectangle:
Nsig = 900.5 +- 36.82 +- 10.95 (38.41)
Nrho = 41.98 +- 2.05 +- 1.28 (2.417)
Ncmb = 549.6 +- 5.183 +- 4.715 (7.006)
Pury = 0.6035 +- 0.02574
Ellips:
Nsig = 897.4 +- 36.69 +- 11.01 (38.31)
Nrho = 38.13 +- 1.862 +- 1.165 (2.197)
Ncmb = 524.5 +- 4.946 +- 4.504 (6.69)
Pury = 0.6146 +- 0.02624
Elli:
Nsig = 860.4 +- 35.18 +- 11.64 (37.05)
Nrho = 31.05 +- 1.517 +- 0.9522 (1.791)
Ncmb = 434.9 +- 4.101 +- 3.748 (5.556)
Pury = 0.6487 +- 0.02793

// mbc_rho_param = 2
Rectangle:
Nsig = 910.6 +- 36.9 +- 11 (38.5)
Nrho = 38.1 +- 1.612 +- 1.187 (2.002)
Ncmb = 550.7 +- 5.034 +- 4.719 (6.9)
Pury = 0.6073 +- 0.02568
Ellips:
Nsig = 907.4 +- 36.77 +- 11.06 (38.4)
Nrho = 33.29 +- 1.408 +- 1.04 (1.751)
Ncmb = 525.5 +- 4.804 +- 4.508 (6.587)
Pury = 0.6189 +- 0.02619
Elli:
Nsig = 870.8 +- 35.29 +- 11.68 (37.17)
Nrho = 26.03 +- 1.101 +- 0.8162 (1.371)
Ncmb = 435.8 +- 3.983 +- 3.751 (5.471)
Pury = 0.6535 +- 0.02789

// mbc_rho_param = 0
Rectangle:
Nsig = 928 +- 37.17 +- 11.14 (38.8)
Nrho = 33.29 +- 1.447 +- 1.069 (1.799)
Ncmb = 545.8 +- 4.98 +- 4.672 (6.829)
Pury = 0.6158 +- 0.02575
Ellips:
Nsig = 924.8 +- 37.04 +- 11.2 (38.69)
Nrho = 28.04 +- 1.219 +- 0.9029 (1.517)
Ncmb = 520.8 +- 4.752 +- 4.462 (6.519)
Pury = 0.6276 +- 0.02626
Elli:
Nsig = 885.9 +- 35.48 +- 11.84 (37.41)
Nrho = 21.16 +- 0.92 +- 0.6839 (1.146)
Ncmb = 432.1 +- 3.943 +- 3.715 (5.417)
Pury = 0.6616 +- 0.02793