#include "ResConst.h"

ResConst::ResConst(const int mc, const int svd, const bool charged, const int error){
  rnp_kink_xi = 3.5;  /* for Rnp from 2010*/
  rnp_kink_st = 0.75; /* for Rnp from 2010*/
  const bool my_Rdet_params = false;
  const bool my_Rdet_asc_params = false;
  const bool my_Rnp_params = false;

  if(mc && svd == 2 && !error){
  //*** Rdet ***/
  //* multiple track vertex */
  if(!my_Rdet_params){ Srec[0]      = +9.271430e-01; Srec[1]        = +2.103700e-01;} //* Srec[2] */
  else               { Srec[0]      = +0.715815;     Srec[1]        = +0.152736; } //* Srec[2] */
  ftl_rec_mlt[0]    = +0.000000e+00; ftl_rec_mlt[1] = +0.000000e+00; //* ftl_rec_mlt[2] */
  Stl_rec_mlt       = +1.000000e+00; //* Stl_rec_mlt */
  if(!my_Rdet_asc_params){ Sasc[0]           = +8.210910e-01; Sasc[1]        = +1.408000e-01;} //* Sasc[2] */
//  else                   { Sasc[0]           = +1.02836e+00;  Sasc[1]        = +1.96229e-01; } //* Sasc[2] */
  else                   { Sasc[0]           = +0.866;  Sasc[1]       = +0.0992; } //* Sasc[2] */
  Snp               = +0.000000e+00; Snp_global     = +1.000000e+00; //* Snp, Snp_global */
  ftl_asc_mlt[0]    = +0.000000e+00; ftl_asc_mlt[1] = +0.000000e+00; //* ftl_asc_mlt[2] */
  Stl_asc_mlt       = +1.000000e+00; //* Stl_asc_mlt */
  //* single track vertex */
  if(!my_Rdet_params){
    Smn_rec           = +1.053040e+00; //* Smn_rec */
    Stl_rec           = +4.320600e+00; //* Stl_rec */
    ftl_rec           = +0.0535341; //* ftl_rec */
  } else{
    Smn_rec           = +1.168; //* Smn_rec */
    Stl_rec           = +4.38; //* Stl_rec */
    ftl_rec           = +0.2015; //* ftl_rec */
  }
  if(!my_Rdet_asc_params){
    Smn_asc           = +1.053040e+00; //* Smn_asc */
    Stl_asc           = +4.320600e+00; //* Stl_asc */
    ftl_asc           = +7.068970e-02; //* ftl_asc */
  } else{
    Smn_asc           = +1.027; //* Smn_asc */
    Stl_asc           = +3.30; //* Stl_asc */
    ftl_asc           = +0.08; //* ftl_asc */
  }
  //*** Rnp ***/
  //* single track vertex */
  if(!my_Rnp_params){
    fd_np_sgl[0]      = +7.737830e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.013430e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.625960e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +9.181350e-01; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  } else{
    fd_np_sgl[0]      = +7.737830e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = 0.798378; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.000000e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = 0.000150321; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  }
  //* multiple track vertex */
  if(!my_Rnp_params){
    fd_np_mlt[0]      = +5.600550e-01; fd_np_mlt[1]    = +7.507230e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +1.569090e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -2.021350e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.324180e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.224260e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = +3.871670e-02; tau_np_p_mlt[1] = +7.653070e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.632670e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +3.214540e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +8.292430e-02; tau_np_n_mlt[1] = +5.342620e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -3.013640e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.899160e-01; //* tau_np_n_mlt_stxi[2] */
  } else{
    fd_np_mlt[0]      = 1.; fd_np_mlt[1]    = 0.292532; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = -0.187605; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -0.31194; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = 0.285766; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = 0; //* fp_np_mlt[2] */
    fn_np_mlt         = 0.178399; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = 0.136984; tau_np_p_mlt[1] = 0.713229; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -0.0388166; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = 0.258934; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = 0.319613; tau_np_n_mlt[1] = 0.181811; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -0.0343582; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = 0.126289; //* tau_np_n_mlt_stxi[2] */
  }
  // Outlayer MC SVD2
  sigma_ol          = +3.063190e+01; //* sigma_ol */
//  f_ol_sgl          = +2.231740e-02; //* f_ol_sgl */
//  f_ol_mul          = +8.843970e-05; //* f_ol_mul */
  f_ol_sgl          = +0.0171085;// +- 0.000693955 //* f_ol_sgl */
  f_ol_mul          = +0.000655638;// 8.60511e-05 //* f_ol_mul */
  }
  if(mc && svd == 1 && !error){
  //*** Rdet ***/
  //* multiple track vertex */
  if(!my_Rdet_params){ Srec[0]           = +9.626350e-01; Srec[1]        = +1.985560e-01;} //* Srec[2] */
  else               { Srec[0]           = +0.395;  Srec[1]       = +0.091; } //* Srec[2] */
  ftl_rec_mlt[0]    = +0.000000e+00; ftl_rec_mlt[1] = +0.000000e+00; //* ftl_rec_mlt[2] */
  Stl_rec_mlt       = +1.000000e+00; //* Stl_rec_mlt */
  if(!my_Rdet_asc_params){ Sasc[0]           = +7.290850e-01; Sasc[1]        = +1.719270e-01;} //* Sasc[2] */
  else               { Sasc[0]           = +0.877;  Sasc[1]        = +0.087;  }//* Sasc[2] */
  Snp               = +0.000000e+00; Snp_global     = +1.000000e+00; //* Snp, Snp_global */
  ftl_asc_mlt[0]    = +0.000000e+00; ftl_asc_mlt[1] = +0.000000e+00; //* ftl_asc_mlt[2] */
  Stl_asc_mlt       = +1.000000e+00; //* Stl_asc_mlt */
  //* single track vertex */
  if(!my_Rdet_params){
    Smn_rec           = +1.109750e+00; //* Smn_rec */
    Stl_rec           = +1.000000e+00; //* Stl_rec */
    ftl_rec           = +0.000000e+00; //* ftl_rec */
  } else{
    Smn_rec           = +1.210; //* Smn_rec */
    Stl_rec           = +4.20; //* Stl_rec */
    ftl_rec           = +0.1810; //* ftl_rec */
  }
//  Smn_asc           = +1.109750e+00; //* Smn_asc */
  if(!my_Rdet_asc_params){
    Smn_asc           = +1.16355e+00; //* Smn_asc */
  } else{
    Smn_asc           = +1.085; //* Smn_asc */
  }
  Stl_asc           = +1.000000e+00; //* Stl_asc */
  ftl_asc           = +0.000000e+00; //* ftl_asc */

  //*** Rnp ***/
  //* single track vertex */
  if(!my_Rnp_params){
    fd_np_sgl[0]      = +7.817180e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.186460e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.847670e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +2.041140e+00; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  } else{
    fd_np_sgl[0]      = +7.817180e-01; fd_np_sgl[1]    = 0.8496; //* fd_np_sgl[2][2] */
    fp_np_sgl         = 0.798855; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = 1.1171; tau_np_p_sgl[1] = 0.981453; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = 2.35317; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  }
  //* multiple track vertex */
  if(!my_Rnp_params){
    fd_np_mlt[0]      = +4.664410e-01; fd_np_mlt[1]    = +6.371510e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +2.706290e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -2.204070e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.228050e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.232800e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = -5.202290e-03; tau_np_p_mlt[1] = +7.168080e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.966400e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.514870e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +4.517970e-02; tau_np_n_mlt[1] = +5.151980e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -7.831950e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +4.304680e-01; //* tau_np_n_mlt_stxi[2] */
  } else{
    fd_np_mlt[0]      = +4.664410e-01; fd_np_mlt[1]    = +6.371510e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +2.706290e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -2.204070e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.228050e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.232800e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = -5.202290e-03; tau_np_p_mlt[1] = +7.168080e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.966400e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.514870e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +4.517970e-02; tau_np_n_mlt[1] = +5.151980e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -7.831950e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +4.304680e-01; //* tau_np_n_mlt_stxi[2] */
  }
  // Outlayer MC svd1
  sigma_ol          = +3.319470e+01; //* sigma_ol */
//  f_ol_sgl          = +3.784100e-02; //* f_ol_sgl */
  f_ol_sgl          = +0.0117525;// +- 0.00073445 //* f_ol_sgl */
  f_ol_mul          = +2.146080e-04; //* f_ol_mul */
  }
  if(!mc && svd == 1 && !error){
  //*** Rdet ***/
  //* multiple track vertex */
  Srec[0]           = +7.046620e-01; Srec[1]        = +2.120840e-01; //* Srec[2] */
  ftl_rec_mlt[0]    = +0.000000e+00; ftl_rec_mlt[1] = +0.000000e+00; //* ftl_rec_mlt[2] */
  Stl_rec_mlt       = +1.000000e+00; //* Stl_rec_mlt */
  Sasc[0]           = +4.834940e-01; Sasc[1]        = +2.366200e-01; //* Sasc[2] */
  Snp               = +0.000000e+00; Snp_global     = +1.057080e+00; //* Snp, Snp_global */
  ftl_asc_mlt[0]    = +0.000000e+00; ftl_asc_mlt[1] = +0.000000e+00; //* ftl_asc_mlt[2] */
  Stl_asc_mlt       = +1.000000e+00; //* Stl_asc_mlt */
  //* single track vertex */
  Smn_rec           = +9.798010e-01; //* Smn_rec */
  Stl_rec           = +1.000000e+00; //* Stl_rec */
  ftl_rec           = +0.000000e+00; //* ftl_rec */
  Smn_asc           = +9.798010e-01; //* Smn_asc */
  Stl_asc           = +1.000000e+00; //* Stl_asc */
  ftl_asc           = +0.000000e+00; //* ftl_asc */
  //*** Rnp ***/
  //* single track vertex */
  fd_np_sgl[0]      = +7.817180e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
  fp_np_sgl         = +8.186460e-01; //* fp_np_sgl[2] */
  tau_np_p_sgl[0]   = +1.847670e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
  tau_np_n_sgl[0]   = +2.041140e+00; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  //* multiple track vertex */
  fd_np_mlt[0]      = +4.664410e-01; fd_np_mlt[1]    = +6.371510e-01; //* fd_np_mlt[2][2] */
  fd_np_st_mlt      = +2.706290e-01; //* fd_np_mlt_st[2] */
  fd_np_xi_mlt      = -2.204070e-01; //* fd_np_mlt_xi[2] */
  fd_np_stxi_mlt    = +2.228050e-01; //* fd_np_mlt_stxi[2] */
  fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
  fn_np_mlt         = +1.232800e-01; //* fn_np_mlt[2] */
  tau_np_p_mlt[0]   = -5.202290e-03; tau_np_p_mlt[1] = +7.168080e-01; //* tau_np_p_mlt[2][2] */
  tau_np_p_xi_mlt   = -2.966400e-02; //* tau_np_p_mlt_xi[2] */
  tau_np_p_stxi_mlt = +2.514870e-01; //* tau_np_p_mlt_stxi[2] */
  tau_np_n_mlt[0]   = +4.517970e-02; tau_np_n_mlt[1] = +5.151980e-01; //* tau_np_n_mlt[2][2] */
  tau_np_n_xi_mlt   = -7.831950e-02; //* tau_np_n_mlt_xi[2] */
  tau_np_n_stxi_mlt = +4.304680e-01; //* tau_np_n_mlt_stxi[2] */
  // Outlayer
  sigma_ol          = +4.369930e+01; //* sigma_ol */
  f_ol_sgl          = +3.700380e-02; //* f_ol_sgl */
  f_ol_mul          = +1.141840e-04; //* f_ol_mul */
  }
  if(!mc && svd == 2 && error == -1){// svd2 data nega error
  //*** Rdet ***/
  //* multiple track vertex */
    Srec[0]        = -1.478871e-01; Srec[1]        = -5.859953e-02; //* Srec[2] */
    ftl_rec_mlt[0] = +0.000000e+00; ftl_rec_mlt[1] = +0.000000e+00; //* ftl_rec_mlt[2] */
    Stl_rec_mlt    = +0.000000e+00; //* Stl_rec_mlt */
    Sasc[0]        = -7.211574e-02; Sasc[1]        = -5.239444e-02; //* Sasc[2] */
    Snp            = +0.000000e+00; Snp_global     = -1.643656e-01; //* Snp, Snp_global */
    ftl_asc_mlt[0] = +0.000000e+00; ftl_asc_mlt[1] = +0.000000e+00; //* ftl_asc_mlt[2] */
    Stl_asc_mlt    = +0.000000e+00; //* Stl_asc_mlt */
  //* single track vertex */
    Smn_rec        = -3.897418e-02; //* Smn_rec */
    Stl_rec        = -3.911839e-01; //* Stl_rec */
    ftl_rec        = -3.946462e-02; //* ftl_rec */
    Smn_asc        = -3.897418e-02; //* Smn_asc */
    Stl_asc        = -3.911839e-01; //* Stl_asc */
    ftl_asc        = -3.946462e-02; //* ftl_asc */
  //*** Rnp ***/
  //* single track vertex */
    fd_np_sgl[0]      = -8.599940e-03; fd_np_sgl[1] = +0.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = -1.769460e-02; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = -4.888840e-02; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = -7.819400e-02; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  //* multiple track vertex */
    fd_np_mlt[0]      = -1.248950e-02; fd_np_mlt[1] = -1.255940e-02; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = -2.910630e-02; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -4.845690e-03; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = -1.150930e-02; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = -2.553030e-03; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = -5.430320e-03; tau_np_p_mlt[1] = -2.052420e-02; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -1.737370e-03; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = -7.070600e-03; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = -9.286360e-03; tau_np_n_mlt[1] = -3.092940e-02; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -4.304750e-03; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = -1.567680e-02; //* tau_np_n_mlt_stxi[2] */
  // Outlayer
    sigma_ol          = 9.217143e+00; //* sigma_ol */
    f_ol_sgl          = -4.659891e-03; //* f_ol_sgl */
    f_ol_mul          = -7.195453e-05; //* f_ol_mul */
  }
  if(!mc && svd == 2 && !error){// svd2 data central values
  //*** Rdet ***/
  //* multiple track vertex */
    Srec[0]        = +8.075310e-01; Srec[1]        = +2.326080e-01; //* Srec[2] */
    ftl_rec_mlt[0] = +0.000000e+00; ftl_rec_mlt[1] = +0.000000e+00; //* ftl_rec_mlt[2] */
    Stl_rec_mlt    = +1.000000e+00; //* Stl_rec_mlt */
    Sasc[0]        = +6.429470e-01; Sasc[1]        = +2.290660e-01; //* Sasc[2] */
    Snp            = +0.000000e+00; Snp_global     = +1.014100e+00; //* Snp, Snp_global */
    ftl_asc_mlt[0] = +0.000000e+00; ftl_asc_mlt[1] = +0.000000e+00; //* ftl_asc_mlt[2] */
    Stl_asc_mlt    = +1.000000e+00; //* Stl_asc_mlt */
  //* single track vertex */
    Smn_rec        = +1.014710e+00; //* Smn_rec */
    Stl_rec        = +3.662940e+00; //* Stl_rec */
    ftl_rec        = +1.110420e-01; //* ftl_rec */
    Smn_asc        = +1.014710e+00; //* Smn_asc */
    Stl_asc        = +3.662940e+00; //* Stl_asc */
    ftl_asc        = +1.110420e-01; //* ftl_asc */
  //*** Rnp ***/
  //* single track vertex */
    fd_np_sgl[0]      = +7.737830e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.013430e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.625960e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +9.181350e-01; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
  //* multiple track vertex */
    fd_np_mlt[0]      = +5.600550e-01; fd_np_mlt[1]    = +7.507230e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +1.569090e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -2.021350e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.324180e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.224260e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = +3.871670e-02; tau_np_p_mlt[1] = +7.653070e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.632670e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +3.214540e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +8.292430e-02; tau_np_n_mlt[1] = +5.342620e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -3.013640e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.899160e-01; //* tau_np_n_mlt_stxi[2] */
  // Outlayer
    sigma_ol          = +3.352530e+01; //* sigma_ol */
    f_ol_sgl          = +2.730960e-02; //* f_ol_sgl */
    f_ol_mul          = +1.529520e-04; //* f_ol_mul */
  }

  // /////////////////// //
  //       Charged       //
  // /////////////////// //

  if(!mc && svd == 2 && !error && charged){
    //*** Rnp ***/
    //* single track vertex */
    fd_np_sgl[0]      = +8.062880e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.274010e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +9.860720e-01; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +4.316230e-01; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
    //* multiple track vertex */
    fd_np_mlt[0]      = +5.345170e-01; fd_np_mlt[1]    = +7.087410e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +1.870940e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -1.879030e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.255560e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.210540e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = +2.442560e-02; tau_np_p_mlt[1] = +7.411070e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -1.976910e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.753010e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +6.553030e-02; tau_np_n_mlt[1] = +5.239120e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -1.906810e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.309140e-01; //* tau_np_n_mlt_stxi[2] */
  }
  if(!mc && svd == 1 && !error && charged){
    //*** Rnp ***/
    //* single track vertex */
    fd_np_sgl[0]      = +8.295620e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.495510e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.417050e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +1.982780e+00; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
    //* multiple track vertex */
    fd_np_mlt[0]      = +4.053710e-01; fd_np_mlt[1]    = +6.110950e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +3.052200e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -1.830920e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +1.983780e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.371750e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = -6.925350e-03; tau_np_p_mlt[1] = +6.471710e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.786130e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.420110e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = -3.982540e-02; tau_np_n_mlt[1] = +6.046470e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -3.812870e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.372510e-01; //* tau_np_n_mlt_stxi[2] */
  }
  if(mc && svd == 2 && !error && charged){
    if(my_Rdet_params){
      Srec[0] = +1.053;  Srec[1] = +0.213; //* Srec[2] */
      Smn_rec = +1.070; //* Smn_rec */
      Stl_rec = +3.28; //* Stl_rec */
      ftl_rec = +0.1344; //* ftl_rec */
    }
    //*** Rnp ***/
    //* single track vertex */
    fd_np_sgl[0]      = +8.062880e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.274010e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +9.860720e-01; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +4.316230e-01; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
    //* multiple track vertex */
    fd_np_mlt[0]      = +5.345170e-01; fd_np_mlt[1]    = +7.087410e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +1.870940e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -1.879030e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +2.255560e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.210540e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = +2.442560e-02; tau_np_p_mlt[1] = +7.411070e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -1.976910e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.753010e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = +6.553030e-02; tau_np_n_mlt[1] = +5.239120e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -1.906810e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.309140e-01; //* tau_np_n_mlt_stxi[2] */
  }
  if(mc && svd == 1 && !error && charged){
    if(my_Rdet_params){
      Srec[0] = +1.020;  Srec[1] = +0.159; //* Srec[2] */
      Smn_rec = +1.330; //* Smn_rec */
    }
    //*** Rnp ***/
    //* single track vertex */
    fd_np_sgl[0]      = +8.295620e-01; fd_np_sgl[1]    = +1.000000e+00; //* fd_np_sgl[2][2] */
    fp_np_sgl         = +8.495510e-01; //* fp_np_sgl[2] */
    tau_np_p_sgl[0]   = +1.417050e+00; tau_np_p_sgl[1] = +0.000000e+00; //* tau_np_p_sgl[2][2] */
    tau_np_n_sgl[0]   = +1.982780e+00; tau_np_n_sgl[1] = +0.000000e+00; //* tau_np_n_sgl[2][2] */
    //* multiple track vertex */
    fd_np_mlt[0]      = +4.053710e-01; fd_np_mlt[1]    = +6.110950e-01; //* fd_np_mlt[2][2] */
    fd_np_st_mlt      = +3.052200e-01; //* fd_np_mlt_st[2] */
    fd_np_xi_mlt      = -1.830920e-01; //* fd_np_mlt_xi[2] */
    fd_np_stxi_mlt    = +1.983780e-01; //* fd_np_mlt_stxi[2] */
    fp_np_mlt         = +0.000000e+00; //* fp_np_mlt[2] */
    fn_np_mlt         = +1.371750e-01; //* fn_np_mlt[2] */
    tau_np_p_mlt[0]   = -6.925350e-03; tau_np_p_mlt[1] = +6.471710e-01; //* tau_np_p_mlt[2][2] */
    tau_np_p_xi_mlt   = -2.786130e-02; //* tau_np_p_mlt_xi[2] */
    tau_np_p_stxi_mlt = +2.420110e-01; //* tau_np_p_mlt_stxi[2] */
    tau_np_n_mlt[0]   = -3.982540e-02; tau_np_n_mlt[1] = +6.046470e-01; //* tau_np_n_mlt[2][2] */
    tau_np_n_xi_mlt   = -3.812870e-02; //* tau_np_n_mlt_xi[2] */
    tau_np_n_stxi_mlt = +3.372510e-01; //* tau_np_n_mlt_stxi[2] */
  }

  return;
}

