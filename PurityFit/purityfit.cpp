#include "purityfit.h"

#include "TChain.h"
#include "TRandom3.h"
#include "TFile.h"

#include "Roo1DTable.h"

#include <iomanip>

void PurityFit::GetShuffledVector(const int size,vector<int> &vec){
  vec.clear();
  for(int i=0; i<size; i++) vec.push_back(i);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  shuffle(vec.begin(),vec.end(),default_random_engine(seed));
  return;
}

PurityFit::PurityFit(const int type, const int _b0f){
  m_fixshape  = false;
  m_singlefbb = false;
  m_bb_or_cnt_flag = 0;
  m_sig_bkg_flag   = 0;
  m_mcflag = true;
  m_nstr = 1;
  m_cstr = 0;

  stringstream out;
//  const double cm2ps = 78.48566945838871754705;
  cuts = new MyParams();
  m_mode = get_mode(type); m_h0mode = get_h0mode(type);
  cout << "mode: " << m_mode << ", h0mode: " << m_h0mode << endl;
  afretcuts = GetCut();

  const int NExp = 26;
  const int NExp_svd1 = 11;
  const int NExp_svd2 = 15;
  int exp_vec[NExp]           = {7,9,11,13,15,17,19,21,23,25,27,31,33,35,37,39,41,43,45,47,49,51,55,61,63,65};
  int exp_vec_svd1[NExp_svd1] = {7,9,11,13,15,17,19,21,23,25,27};
  int exp_vec_svd2[NExp_svd2] =                                {31,33,35,37,39,41,43,45,47,49,51,55,61,63,65};
  exp = new RooCategory("exp","exp");
  for(int i=0; i<NExp; i++){
    out.str(""); out << exp_vec[i];
    exp->defineType(out.str().c_str(),exp_vec[i]);
  }
  exp_svd1 = new RooCategory("exp","exp_svd1");
  for(int i=0; i<NExp_svd1; i++){
    out.str(""); out << exp_vec_svd1[i];
    exp_svd1->defineType(out.str().c_str(),exp_vec_svd1[i]);
  }
  exp_svd2 = new RooCategory("exp","exp_svd2");
  for(int i=0; i<NExp_svd2; i++){
    out.str(""); out << exp_vec_svd2[i];
    exp_svd2->defineType(out.str().c_str(),exp_vec_svd2[i]);
  }

  pt_pip = new RooRealVar("pt_pip","pt_pip",0.0,6.);
  pt_pim = new RooRealVar("pt_pim","pt_pim",0.0,6.);

  pt_pi1 = new RooRealVar("pt_pi1","pt_pi1",0.0,6.);
  pt_pi2 = new RooRealVar("pt_pi2","pt_pi2",0.0,6.);

  p_pi0_h0 = new RooRealVar("p_pi0_h0","p_pi0_h0",0.0,6.);
  cos_hel = new RooRealVar("cos_hel","cos_hel",-1.,1.);

  r_pip = new RooRealVar("r_pip","r_pip",-2.,2.);
  r_pim = new RooRealVar("r_pim","r_pim",-2.,2.);
  z_pip = new RooRealVar("z_pip","z_pip",-5.,5.);
  z_pim = new RooRealVar("z_pim","z_pim",-5.,5.);

  r_pi1 = new RooRealVar("r_pi1","r_pi1",-2.,2.);
  r_pi2 = new RooRealVar("r_pi2","r_pi2",-2.,2.);
  z_pi1 = new RooRealVar("z_pi1","z_pi1",-5.,5.);
  z_pi2 = new RooRealVar("z_pi2","z_pi2",-5.,5.);
  
  e_g1 = new RooRealVar("e_g1","e_g1",0.04,6.);
//  e_g2 = new RooRealVar("e_g2","e_g2",0.04,4.);
//  e_g3 = new RooRealVar("e_g3","e_g3",0.04,4.);
//  e_g4 = new RooRealVar("e_g4","e_g4",0.04,4.);

//  th_g1 = new RooRealVar("th_g1","th_g1",-1.,1.);
//  th_g2 = new RooRealVar("th_g2","th_g2",-1.,1.);
//  th_g3 = new RooRealVar("th_g3","th_g3",-1.,1.);
//  th_g4 = new RooRealVar("th_g4","th_g4",-1.,1.);

  b0f = new RooCategory("b0f","b0f");
  b0f->defineType("signal",1);
  b0f->defineType("fsr",10);
  b0f->defineType("bad_pi0",5);
  b0f->defineType("rho2",2);
  b0f->defineType("rho3",3);
  b0f->defineType("rho4",4);
  b0f->defineType("rho6",6);
  b0f->defineType("rho11",11);
  b0f->defineType("comb",-1);

  rndm_pi0 = new RooCategory("rndm_pi0","rndm_pi0");
  rndm_pi0->defineType("no",0);
  rndm_pi0->defineType("yes",1);

  d0f = new RooCategory("d0f","d0f");
  d0f->defineType("signal",1);
  d0f->defineType("fsr",10);
//  d0f->defineType("bad_pi0",5);
//  d0f->defineType("rho2",2);
  d0f->defineType("rho3",3);
//  d0f->defineType("rho4",4);
//  d0f->defineType("rho6",6);
  d0f->defineType("rho11",11);
  d0f->defineType("comb",-1);

  h0f = new RooCategory("h0f","h0f");
  h0f->defineType("signal",1);
//  h0f->defineType("fsr",10);
//  d0f->defineType("bad_pi0",5);
//  d0f->defineType("rho2",2);
  h0f->defineType("rho3",3);
//  d0f->defineType("rho4",4);
//  d0f->defineType("rho6",6);
  h0f->defineType("rho11",11);
  h0f->defineType("comb",-1);

  mode   = new RooCategory("mode","mode");
  if(m_mode != 100 && m_mode != 200){ mode->defineType("mode",m_mode);}
  else if(m_mode == 100){
    mode->defineType("pi0",1);
    mode->defineType("eta",2);
    mode->defineType("Dst_pi0",10);
    mode->defineType("Dst_eta",20);
  }
  h0mode = new RooCategory("h0mode","h0mode");
  h0mode->defineType("h0mode",m_h0mode);

  flv = new RooCategory("flv","flv");
  flv->defineType("B0",1);
  flv->defineType("anti-B0",-1);

  flv_mc = new RooCategory("flv_mc","flv_mc");
  flv_mc->defineType("B0",1);
  flv_mc->defineType("anti-B0",-1);

  good_icpv = new RooCategory("good_icpv","good_icpv");
  good_icpv->defineType("good",1);

  const int NBins = 16;
  int bins_vec[NBins] = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
  bin = new RooCategory("bin","bin");
  bin_mc = new RooCategory("bin_mc","bin_mc");
  for(int i=0; i<NBins; i++){
    out.str(""); out << bins_vec[i];
    bin->defineType(out.str().c_str(),bins_vec[i]);
    bin_mc->defineType(out.str().c_str(),bins_vec[i]);
  }

  binflv = new RooSuperCategory("binflv","binflv",RooArgSet(*bin,*flv));
  binflv->Print("t");

  mbc_min = cuts->get_mbc_min_h0(m_mode,m_h0mode);
  mbc_max = cuts->get_mbc_max_h0(m_mode,m_h0mode);
  de_min  = cuts->get_de_min_h0(m_mode,m_h0mode);
  de_max  = cuts->get_de_max_h0(m_mode,m_h0mode);

  de = new RooRealVar("de","#DeltaE",cuts->get_de_fit_min(),cuts->get_de_fit_max(),"GeV");
  de->setRange("deSignal",de_min,de_max);
  de->setBins(50);
//  de->setRange("Sideband1",-0.15,0.3);
  de->setRange("Sideband2",0.12,0.3);

//  mbc = new RooRealVar("mbc","M_{bc}",cuts->get_mbc_fit_min(),cuts->get_mbc_fit_max(),"GeV");
  mbc = new RooRealVar("mbc","M_{bc}",cuts->get_mbc_fit_min(),5.2889,"GeV");
  mbc->setBins(50);
  mbc->setRange("mbcSignal",mbc_min,mbc_max);
  mbc->setRange("Sideband1",5.23,5.26);
  mbc->setRange("Sideband2",5.25,5.2889);

  cout << "Init Signal PDF maker..." << endl;
  pdf_gen_sig  = new MEPdfSignal(de,mbc,m_mode,m_h0mode,_b0f);
  pdf_gen_sig->GetParametersFromFile();
  cout << "Init Comb PDF maker..." << endl;
  pdf_gen_cmb  = new MEPdfCombinatorial(de,mbc,m_mode,m_h0mode);
  pdf_gen_cmb->GetParametersFromFile();
  cout << "Init Partial PDF maker..." << endl;
  pdf_gen_part = new MEPdfPartialB(de,mbc,m_mode,m_h0mode);
  pdf_gen_part->GetParametersFromFile();

//  if(m_release_params){
//    pdf_gen_sig->Free_de0();
//    pdf_gen_sig->Free_mbc0();
//    pdf_gen_cmb->Release_de_QQ();
//    if(m_mode == 1 || m_mode == 10) pdf_gen_part->Release_de0();
//  }

  md     = new RooRealVar("md_raw","M(D^{0}})",cuts->get_md_min(),cuts->get_md_max(),"GeV");
  mk     = new RooRealVar("mk","m(K_{S}^{0}})",cuts->get_mk_min(),cuts->get_mk_max(),"GeV");
  mh0    = new RooRealVar("mh0","m(h^{0}})",cuts->get_mh0_min(m_mode,m_h0mode),cuts->get_mh0_max(m_mode,m_h0mode),"GeV");
  mpi0   = new RooRealVar("mpi0","m(#pi^{0}})",cuts->get_mpi0_min(),cuts->get_mpi0_max(),"GeV");
  dmetap = new RooRealVar("dmetap","#DeltaM(#eta`})",cuts->get_dm_etap_min(m_h0mode),cuts->get_dm_etap_max(m_h0mode),"GeV");
  dmdts0 = new RooRealVar("dmdst0","#DeltaM(D^{*}})",cuts->get_dm_dst0_min(),cuts->get_dm_dst0_max(),"GeV");
  chi2_vtx_d0 = new RooRealVar("chi2_vtx_d0","#Chi^{2}/n.d.f. of D^{0} vertex fit",0.,500.);
  bdt    = new RooRealVar("bdt","BDT",cuts->bdt_cut(m_mode,m_h0mode),1.);
  lh0    = new RooRealVar("lh0","lh0",cuts->lh0_cut(m_mode,m_h0mode),1.);
  tag_LH = new RooRealVar("tag_LH","tag_LH",-1.,1.);
  costhBcms = new RooRealVar("costhBcms","costhBcms",-1.,1.);

  mp    = new RooRealVar("mp","mp",0.,3.,"GeV^{2}");
  mm    = new RooRealVar("mm","mm",0.,3.,"GeV^{2}");
  mp_mc = new RooRealVar("mp_mc","mp_mc",0.,3.,"GeV^{2}");
  mm_mc = new RooRealVar("mm_mc","mm_mc",0.,3.,"GeV^{2}");

  const int ntrk_vec_size = 12;
  int ntrk_asc_vec[ntrk_vec_size] = {1,2,3,4,5,6,7,8,9,10,11,12};
  ndf_asc  = new RooCategory("ndf_z_asc","ndf_z_asc");
  ntrk_asc = new RooCategory("ntrk_asc","ntrk_asc");
  for(int i=0; i<ntrk_vec_size; i++){
    out.str("");
    out << ntrk_asc_vec[i];
    ntrk_asc->defineType(out.str().c_str(),ntrk_asc_vec[i]);

    out.str("");
    out << 2*ntrk_asc_vec[i]-1;
    ndf_asc->defineType(out.str().c_str(),2*ntrk_asc_vec[i]-2);
  }

//  exp  = new RooCategory("exp","exp"); argset->add(*exp);
//  run  = new RooCategory("run","run");
//  evtn = new RooCategory("evtn","evtn");

  ndf_rec = new RooCategory("ndf_z_sig","ndf_z_sig");
  ndf_rec->defineType("0",0);
  ndf_rec->defineType("4",4);
  ndf_rec->defineType("8",8);

  ntrk_rec = new RooCategory("ntrk_sig","ntrk_sig");
  ntrk_rec->defineType("1",1);
  ntrk_rec->defineType("3",3);
  ntrk_rec->defineType("5",5);

  z_rec     = new RooRealVar("z_sig","z_sig",-100,100);
  sz_rec    = new RooRealVar("sz_sig","sz_sig",0.,100);
  chisq_rec = new RooRealVar("chisq_z_sig","chisq_z_sig",0.,1e6);
  z_asc     = new RooRealVar("z_asc","z_asc",-100,100);
  sz_asc    = new RooRealVar("sz_asc","sz_asc",0,100);
  chisq_asc = new RooRealVar("chisq_z_asc","chisq_z_asc",0.,1e6);

  z_rec_mc  = new RooRealVar("z_sig_mc","z_sig_mc",-100,100);
  z_asc_mc  = new RooRealVar("z_asc_mc","z_asc_mc",-100,100);

//  argset = new RooArgSet();
//  SetArgSet(argset);
  cout << "PurityFit constructor finished..." << endl;
}

int PurityFit::get_mode(const int type) const{
  switch(type){
  case 1:  return 1; // pi0
  case 10: return 10;// D*0 pi0
  case 2:  return 2; // eta -> gg
  case 20: return 20;// D*0 eta(->gg)
  case 21: return 20;// D*0 eta(->ppp)
  case 3:  return 2; // eta -> pi+pi-pi0
  case 4:  return 3; // omega
  case 5:  return 5; // eta'
  case 100: return 100;// single
  case 200: return 200;// multiple
  }
  return 0;
}

int PurityFit::get_h0mode(const int type) const{
  switch(type){
  case 1:  return 10;// pi0
  case 10: return 10;// D*0 pi0
  case 2:  return 10;// eta -> gg
  case 20: return 10;// D*0 eta(->gg)
  case 21: return 20;// D*0 eta(->ppp)
  case 3:  return 20;// eta -> pi+pi-pi0
  case 4:  return 20;// omega
  case 5:  return 10;// eta'
  case 100: return 10;// single
  case 200: return 20;// multiple
  }
  return 0;
}

string PurityFit::GetCut(void){
  string cuts("pt_pip>0.05 && pt_pim>0.05 && abs(z_pip)<5 && abs(z_pim)<5 && abs(r_pip)<2 && abs(r_pim)<2");
  switch(m_mode){
  case 1:
//    cuts += string(" && (e_g1>0.06 || abs(th_g1)>0.7)");
    break;
  case 2:
    if(m_h0mode == 10){
//      cuts += string(" && e_g1>0.15");
      break;
    } else{
      cuts += string(" && pt_pi1>0.10 && pt_pi2>0.10 && p_pi0_h0>0.2 && abs(z_pi1)<5 && abs(z_pi2)<5 && abs(r_pi1)<2 && abs(r_pi2)<2");
      break;
    }
  case 3:
    cuts += string(" && abs(cos_hel)>0.2 && pt_pi1>0.10 && pt_pi2>0.10 && p_pi0_h0>0.2 && abs(z_pi1)<5 && abs(z_pi2)<5 && abs(r_pi1)<2 && abs(r_pi2)<2");
    break;
  }
  return cuts;
}

RooDataSet* PurityFit::GetSigMCDataSet(void){
  TChain* treeSig = new TChain("TEvent","TEventSig");
  switch(m_mode){
  case 1:
    treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcPi0_s8_m1_h0m10.root");
    break;
  case 10:
    treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1_m10_h0m10.root");
    break;
  case 2:
    if(m_h0mode == 10) treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s3_m2_h0m10.root");
    else               treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETA_s3_m2_h0m20.root");
    break;
  case 20:
    treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcDST0_s1_m20_h0m10.root");
    break;
  case 3:
    treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcOMEGA_s6_m3_h0m20.root");
    break;
  case 5:
    treeSig->Add("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_sigmcETAP_s1_m5_h0m10.root");
    break;
  default:
    break;
  }
  cout << "Events in the TTree: " << treeSig->GetEntries() << endl;
  stringstream out;
  out.str("");
  out << afretcuts;
  out << " && (b0f == 1 || b0f == 5 || b0f == 10 || rndm_pi0 == 1)";
  out << "&& abs(z_sig-z_asc)<70./7.848";
  cout << out.str() << endl;
  string title;
  RooArgSet* aset = GetArgSet(1);
  switch(m_svd){
  case 2:
    title = string("ds_sig_svd2");
    break;
  case 1:
    title = string("ds_sig_svd1");
    break;
  default:
    title = string("ds_sig");
    break;
  }
  RooDataSet* dsSig = new RooDataSet(title.c_str(),title.c_str(),treeSig,*aset,out.str().c_str());
  delete treeSig;
  return dsSig;
}

RooDataSet* PurityFit::GetSigMCLineDataSet(const string& angle){
  TChain* treeSig = new TChain("TEvent","TEventSig");
  const string prefix("/home/vitaly/B0toDh0/Tuples/Linearity/Fil_b2dh_sigmc");
  stringstream out;
  string infile;
  switch(m_mode){
  case 1:
    out.str("");
    out << prefix << "Pi0" << angle << "_s1.root";
    infile = out.str();
    break;
  case 10:
    out.str("");
    out << prefix << "DST0" << angle << "_s1.root";
    infile = out.str();
    break;
  case 2:
    out.str("");
    out << prefix << "ETA" << angle << "_s1.root";
    infile = out.str();
    break;
  case 20:
    out.str("");
    out << prefix << "DST0" << angle << "_s1.root";
    infile = out.str();
    break;
  case 3:
    out.str("");
    out << prefix << "OMEGA" << angle << "_s1.root";
    infile = out.str();
    cout << infile << endl;
    break;
  case 5:
    out.str("");
    out << prefix << "ETAP" << angle << "_s1.root";
    infile = out.str();
    break;
  default:
    break;
  }
  treeSig->Add(infile.c_str());
  cout << "Events in the TTree: " << treeSig->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && (b0f == 1 || b0f == 5 || b0f == 10)";
  out << "&& abs(z_sig-z_asc)<70./7.848";
  cout << out.str() << endl;
  string title;
  RooArgSet* aset = GetArgSet(1);
  switch(m_svd){
  case 2:
    title = string("ds_sig_svd2");
    break;
  case 1:
    title = string("ds_sig_svd1");
    break;
  default:
    title = string("ds_sig");
    break;
  }
  RooDataSet* dsSig = new RooDataSet(title.c_str(),title.c_str(),treeSig,*aset,out.str().c_str());
  delete treeSig;
  return dsSig;
}

RooDataSet* PurityFit::GetGenMCDataSet(const vector<int> streams){
//  sig_bkg_flag 0 -> normal mode
//  sig_bkg_flag 1 -> signal only
//  sig_bkg_flag 2 -> background only
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    if(m_bb_or_cnt_flag != 2){
      line = prefix + string("uds")     + out.str();
      tree->Add(line.c_str());
      line = prefix + string("charm")   + out.str();
      tree->Add(line.c_str());
    }
    if(m_bb_or_cnt_flag != 1){
      line = prefix + string("mixed")   + out.str();
      tree->Add(line.c_str());
      line = prefix + string("charged") + out.str();
      tree->Add(line.c_str());
    }
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && b0f != 0 && b0f>=-1 && abs(z_sig-z_asc)<70./7.848";
  if(m_sig_bkg_flag == 1) out << " && (b0f == 1 || b0f == 5 || b0f == 10 || rndm_pi0 == 1)";
  if(m_sig_bkg_flag == 2) out << " && !(b0f == 1 || b0f == 5 || b0f == 10) && rndm_pi0 == 0";
  string title;
  RooArgSet* aset = GetArgSet(2);
  aset->Print();
  switch(m_svd){
  case 2:
    title = string("ds_genmc_svd2");
    break;
  case 1:
    title = string("ds_genmc_svd1");
    break;
  default:
    title = string("ds_genmc");
    break;
  }
  RooDataSet* dsGen = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsGen;
}

RooDataSet* PurityFit::GetDataSet(const int svd){
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/Tuples/Fil_b2dh_data_m" << m_mode << "_h0m" << m_h0mode << ".root";
  tree->Add(out.str().c_str());
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && abs(z_sig-z_asc)<70./7.848";
  string title;
  RooArgSet* aset = GetArgSet(0);
  aset->Print();
  switch(svd){
  case 2:
    title = string("ds_data_svd2");
    break;
  case 1:
    title = string("ds_data_svd1");
    break;
  default:
    title = string("ds_data");
    break;
  }
  RooDataSet* dsGen = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsGen;
}

void PurityFit::SaveSidebandTree(const int svd){
  RooDataSet* ds = GetMbcSidebandSet(svd);
  const int NTot = ds->sumEntries();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  ICPVEv ev; ev.f_bkg = 1; ev.sigarea = 0;
  stringstream out;
  out.str("");
  out << "data/";
  out << "data_sideband_tree_m" << m_mode << "_hm" << m_h0mode << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  double cmb_val, cnt_val, prt_val;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,2);
    de->setVal(ev.de); mbc->setVal(ev.mbc);
    cmb_val = pdf_comb->getVal(RooArgSet(*de,*mbc));
    cnt_val = pdf_cont->getVal(RooArgSet(*de,*mbc));
    prt_val = pdf_part->getVal(RooArgSet(*de,*mbc));
    ev.f_cont_in_comb = f_cont_in_comb*cnt_val/cmb_val;
    ev.f_cont         = ev.f_cont_in_comb*cmb_val/(cmb_val+prt_val);
    tree->Fill();
  }
  tree->Write();
  file->Close();
}

void PurityFit::SaveSidebandTree(const vector<int> streams, const int svd, const int bb_or_cnt_flag){
  m_svd = svd;
  m_bb_or_cnt_flag = bb_or_cnt_flag;
  RooDataSet* ds = GetGenMCMbcSidebandSet(streams);

  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  const int NTot = ds->sumEntries();
  ICPVEv ev;
  ev.f_bkg = 1;
  ev.sigarea = 0;
  stringstream out;
  out.str("");
  out << "data/";
  out << "genmc_sideband_tree_m" << m_mode << "_h0m" << m_h0mode;
  if(m_bb_or_cnt_flag == 1)      out << "_cont";
  else if(m_bb_or_cnt_flag == 2) out << "_BB";
  out << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  double cmb_val, cnt_val, prt_val;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,2);
    de->setVal(ev.de); mbc->setVal(ev.mbc);
    cmb_val = pdf_comb->getVal(RooArgSet(*de,*mbc));
    cnt_val = pdf_cont->getVal(RooArgSet(*de,*mbc));
    prt_val = pdf_part->getVal(RooArgSet(*de,*mbc));
    ev.f_cont_in_comb = f_cont_in_comb*cnt_val/cmb_val;
    ev.f_cont         = ev.f_cont_in_comb*cmb_val/(cmb_val+prt_val);
    tree->Fill();
  }
  tree->Write();
  file->Close();
}

RooDataSet* PurityFit::GetMbcSidebandSet(const int svd){
  m_svd = svd;
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  out.str("");
  out << "/home/vitaly/B0toDh0/Tuples/Fil_b2dh_data_m" << m_mode << "_h0m" << m_h0mode << ".root";
  tree->Add(out.str().c_str());
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  string title;
  RooArgSet* aset = GetArgSet(2);
  aset->Print();
  switch(m_svd){
  case 2:
    title = string("ds_sideband_svd2");
    break;
  case 1:
    title = string("ds_sideband_svd1");
    break;
  default:
    title = string("ds_sideband");
    break;
  }
  out.str("");
  out << "abs(z_sig-z_asc)<70./7.848";
//  out << " && mbc<5.25";
  if(m_svd == 1) out << " && exp<30";
  if(m_svd == 2) out << " && exp>30";
  out << " && (mbc>5.23 && mbc<5.25 || mbc>5.25 && de>0.1)";
  RooDataSet* dsGen = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsGen;
}

RooDataSet* PurityFit::GetGenMCMbcSidebandSet(const vector<int> streams){
  // flag == 0 -> all events
  // flag == 1 -> continuum
  // flag == 2 -> BBbar
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    if(m_bb_or_cnt_flag != 2){
      line = prefix + string("uds")     + out.str();
      tree->Add(line.c_str());
      line = prefix + string("charm")   + out.str();
      tree->Add(line.c_str());
    }
    if(m_bb_or_cnt_flag != 1){
      line = prefix + string("mixed")   + out.str();
      tree->Add(line.c_str());
      line = prefix + string("charged") + out.str();
      tree->Add(line.c_str());
    }
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << "b0f != 0 && b0f>=-1 && abs(z_sig-z_asc)<70./7.848";
//  out << " && mbc<5.25";
  out << " && (mbc>5.23 && mbc<5.25 || mbc>5.25 && de>0.1)";
  cout << out.str() << endl;
  string title;
  RooArgSet* aset = GetArgSet(2);
  aset->Print();
  switch(m_svd){
  case 2:
    title = string("ds_genmc_sideband_svd2");
    break;
  case 1:
    title = string("ds_genmc_sideband_svd1");
    break;
  default:
    title = string("ds_genmc_sideband");
    break;
  }
  RooDataSet* dsGen = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsGen;
}

RooDataSet* PurityFit::GetCombMCDataSet(const vector<int> streams, const int svd){
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("uds")     + out.str();
    tree->Add(line.c_str());
    line = prefix + string("charm")   + out.str();
    tree->Add(line.c_str());
    line = prefix+ string("mixed")   + out.str();
    tree->Add(line.c_str());
    line = prefix + string("charged") + out.str();
    tree->Add(line.c_str());
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && b0f==-1 && rndm_pi0 == 0";
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(svd){
  case 2:
    title = string("ds_comb_svd2");
    break;
  case 1:
    title = string("ds_comb_svd1");
    break;
  default:
    title = string("ds_comb");
    break;
  }
  RooDataSet* dsComb = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsComb;
}

RooDataSet* PurityFit::GetRawBackMCDataSet(const vector<int> streams,const int svd){
  TChain* treeRawBack = new TChain("TEvent","TEventRawBack");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("uds")     + out.str();
    treeRawBack->Add(line.c_str());
    line = prefix + string("charm")   + out.str();
    treeRawBack->Add(line.c_str());
    line = prefix + string("mixed")   + out.str();
    treeRawBack->Add(line.c_str());
    line = prefix + string("charged") + out.str();
    treeRawBack->Add(line.c_str());
  }
  cout << "Events in the TTree: " << treeRawBack->GetEntries() << endl;
  out.str("");
  out << "b0f!=1 && b0f!=5 && b0f!=10";
  out << " && abs(z_sig-z_asc)<70./7.848";
  cout << out.str() << endl;
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(svd){
  case 2:
    title = string("ds_rawback_svd2");
    break;
  case 1:
    title = string("ds_rawback_svd1");
    break;
  default:
    title = string("ds_rawback");
    break;
  }
  RooDataSet* dsRawBack = new RooDataSet(title.c_str(),title.c_str(),treeRawBack,*aset,out.str().c_str());
  delete treeRawBack;
  delete aset;
  return dsRawBack;
}

RooDataSet* PurityFit::GetContMCDataSet(const vector<int> streams){
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("uds")     + out.str();
    tree->Add(line.c_str());
    line = prefix + string("charm")   + out.str();
    tree->Add(line.c_str());
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && b0f != 0 && b0f>=-1 && abs(z_sig-z_asc)<70./7.848";
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(m_svd){
  case 2:
    title = string("ds_cont_svd2");
    break;
  case 1:
    title = string("ds_cont_svd1");
    break;
  default:
    title = string("ds_cont");
    break;
  }
  RooDataSet* dsCont = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsCont;
}

RooDataSet* PurityFit::GetBBMCDataSet(const vector<int> streams){
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("charged") + out.str();
    tree->Add(line.c_str());
    line = prefix + string("mixed")   + out.str();
    tree->Add(line.c_str());
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << "b0f != 0 && b0f>=-1 && b0f != 1 && b0f != 5 && b0f != 10 && abs(z_sig-z_asc)<70./7.848";
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(m_svd){
  case 2:
    title = string("ds_bb_svd2");
    break;
  case 1:
    title = string("ds_bb_svd1");
    break;
  default:
    title = string("ds_bb");
    break;
  }
  RooDataSet* dsCont = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsCont;
}

RooDataSet* PurityFit::GetPartMCDataSet(const vector<int> streams){
  TChain* tree = new TChain("TEvent","TEvent");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("mixed")   + out.str();
    tree->Add(line.c_str());
    line = prefix + string("charged") + out.str();
    tree->Add(line.c_str());
  }
  cout << "Events in the TTree: " << tree->GetEntries() << endl;
  out.str("");
  out << afretcuts;
  out << " && (b0f == 2 || b0f == 3 || b0f == 4 || b0f == 6 || b0f == 11 && abs(z_sig-z_asc)<70./7.848)";
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(m_svd){
  case 2:
    title = string("ds_part_svd2");
    break;
  case 1:
    title = string("ds_part_svd1");
    break;
  default:
    title = string("ds_part");
    break;
  }
  RooDataSet* dsPart = new RooDataSet(title.c_str(),title.c_str(),tree,*aset,out.str().c_str());
  delete tree;
  delete aset;
  return dsPart;
}

RooDataSet* PurityFit::GetBBBackMCDataSet(const vector<int> streams){
  TChain* treeBack = new TChain("TEvent","TEventBBBack");
  stringstream out;
  for(int i=0; i<streams.size(); i++){
    out.str("");
    out << "_" << streams[i] << "_1" << streams[i] << "_m" << m_mode << "_h0m" << m_h0mode << ".root";
    const string prefix("/home/vitaly/B0toDh0/Tuples/Fil_b2dh_");
    string line;
    line = prefix + string("mixed")   + out.str();
    treeBack->Add(line.c_str());
    line = prefix + string("charged") + out.str();
    treeBack->Add(line.c_str());
  }
  cout << "Events in the TTree: " << treeBack->GetEntries() << endl;
  out.str("");
  out << "b0f!=1 && b0f!=5 && b0f!=10";
  out << " && abs(z_sig-z_asc)<70./7.848";
  cout << out.str() << endl;
  string title;
  RooArgSet* aset = GetArgSet(2);
  switch(m_svd){
  case 2:
    title = string("ds_bbback_svd2");
    break;
  case 1:
    title = string("ds_bbback_svd1");
    break;
  default:
    title = string("ds_bbback");
    break;
  }
  RooDataSet* dsBack = new RooDataSet(title.c_str(),title.c_str(),treeBack,*aset,out.str().c_str());
  delete treeBack;
  delete aset;
  return dsBack;
}

void PurityFit::FitSigPdf(void){
  RooDataSet* ds = GetSigMCDataSet();
  ds->Print();
  pdf_gen_sig->FitParameters(ds);
  pdf_gen_sig->WriteParameters();
  return;
}

void PurityFit::FitMbcSigPdf(void){
  RooDataSet* ds = GetSigMCDataSet();
  ds->Print();
  pdf_gen_sig->FitMbcParameters(ds);
  pdf_gen_sig->WriteParameters();
  return;
}

void PurityFit::CheckSigPdf(void){
  RooDataSet* ds = GetSigMCDataSet();
  ds->Print();
  pdf_gen_sig->TryParameters(ds);
  return;
}

void PurityFit::FitCombPdf(const vector<int> streams){
  FitContPdf(streams);
  RooDataSet* ds = GetCombMCDataSet(streams);
  stringstream out;
  out.str("");
  out << "params/BinMap_comb_m" << m_mode << "_h0m" << m_h0mode << ".txt";
  WriteBinFlvMap(ds,out.str().c_str());
  pdf_gen_cmb->FitParameters(ds);
  pdf_gen_cmb->WriteParameters();
  return;
}

void PurityFit::CheckCombPdf(const vector<int> streams){
  RooDataSet* ds = GetCombMCDataSet(streams);
  pdf_gen_cmb->TryParameters(ds);
  return;
}

void PurityFit::FitContPdf(const vector<int> streams){
  RooDataSet* ds = GetContMCDataSet(streams);
  stringstream out;
  out.str("");
  out << "params/BinMap_cont_m" << m_mode << "_h0m" << m_h0mode << ".txt";
  WriteBinFlvMap(ds,out.str().c_str());
  pdf_gen_cmb->FitContParameters(ds);
  pdf_gen_cmb->WriteParameters(true);
  return;
}

void PurityFit::CheckContPdf(const vector<int> streams){
  RooDataSet* ds = GetContMCDataSet(streams);
  pdf_gen_cmb->TryContParameters(ds);
  return;
}

void PurityFit::FitPartPdf(const vector<int> streams){
  RooDataSet* ds = GetPartMCDataSet(streams);
  stringstream out;
  out.str("");
  out << "params/BinMap_part_m" << m_mode << "_h0m" << m_h0mode << ".txt";
  WriteBinFlvMap(ds,out.str().c_str());
  bool etaggflag = false;
  if(m_mode == 2 && m_h0mode == 10) etaggflag = true;

  pdf_gen_part->FitParameters(ds,etaggflag);
  pdf_gen_part->WriteParameters();
  return;
}

void PurityFit::CheckPartPdf(const vector<int> streams){
  RooDataSet* ds = GetPartMCDataSet(streams);
  pdf_gen_part->TryParameters(ds);
  return;
}

void PurityFit::MakeGenMCPutiryFit(const vector<int> streams,const bool fixshape, const int svd, const int icpv_flag, const int bb_or_cnt_flag, const int sig_bkg_flag,const int nstr,const int cstr){
  m_svd = svd;
  m_bb_or_cnt_flag = bb_or_cnt_flag;
  m_sig_bkg_flag = sig_bkg_flag;
  m_fixshape = fixshape;
  m_svd = svd;
  m_cstr = cstr;
  m_nstr = nstr;
  m_mcflag = true;

  RooDataSet* ds    = GetGenMCDataSet(streams);
  RooDataSet* cntds = GetContMCDataSet(streams);
  stringstream out;
  out.str("");
  out << "params/BinMap_genmc_m" << m_mode << "_h0m" << m_h0mode << ".txt";
  WriteBinFlvMap(ds,out.str().c_str());
  Fit(ds);

  if(icpv_flag>1){
    m_svd = 1;
    sig_ds_svd1 = GetSigMCDataSet();
    m_svd = 2;
    sig_ds_svd2 = GetSigMCDataSet();
    m_svd = svd;
  }

  DefineElliRange();
  CalcSigIntegrals();
  CalcSB1Integrals();
  CalcSB2Integrals();
  CountTrueNumbers(ds,cntds);
  Draw(*pdf,ds);

  PrintIntegrals();
  PrintTrueNumbers();
  WriteIntegrals();
  WriteTrueNumbers();

  CalcBinsFractions(ds);

  switch(icpv_flag){
  case 1:// Save GenMC tree
    SaveCPVTree(ds);
    break;
  case 2:// Save mixed tree
    MixCPVTree(ds);
    break;
  case 3:// Save large mixed tree
//    MixBigCPVTree(ds,streams);
    break;
  default:
    break;
  }
  return;
}

void PurityFit::MakeGenMCPutiryFit2(const vector<int> streams, const bool singlefbb, const bool fixshape, const int svd, const int nstr, const int cstr){
  m_svd = svd;
  m_fixshape = fixshape;
  m_singlefbb = singlefbb;
  m_mcflag = true;
  m_nstr = nstr;
  m_cstr = cstr;
  cout << "MakeGenMCPutiryFit2..." << endl;
  RooDataSet* ds = GetGenMCDataSet(streams);
  RooDataSet* cntds = GetContMCDataSet(streams);
  cout << "  dataset..." << endl;
//  CalcBinsFractions(ds,true);
  CalcWTagAndNEv(ds,WrTagMap,EventsMap,false,false);
  cout << "  wrongtag..." << endl;

  stringstream out;
  out.str("");
  out << "params/BinMap2_genmc_m" << m_mode << "_h0m" << m_h0mode << ".txt";
  WriteBinFlvMap(ds,out.str().c_str());

  cout << "  binflvmap..." << endl;
  Fit2(ds);

  out.str("");
  out << "params/Predictions2_m" << m_mode << "_h0m" << m_h0mode;
  out << "_mc";
  if(m_singlefbb) out << "_sglfbb";
  if(m_fixshape)  out << "_fixedshape";
  out << ".txt";
  MakePredictions2(EventsMap,NSigPredicted,NSigPredicted_err,NCmbPredicted,NCntPredicted,NPrtPredicted,out.str().c_str());

  m_svd = 1;
  sig_ds_svd1 = GetSigMCDataSet();
  m_svd = 2;
  sig_ds_svd2 = GetSigMCDataSet();
  m_svd = svd;

  DefineElliRange();
  CalcSigIntegrals2();
//  CalcSB1Integrals2();
//  CalcSB2Integrals2();
  CountTrueNumbers(ds,cntds);

//  PrintIntegrals2();
//  PrintTrueNumbers2();
  WriteIntegrals2();
  WriteTrueNumbers();

  DrawDeltaE2(ds);
  MixCPVTree2(ds);
  return;
}

void PurityFit::MakeFit(const bool fixshape){
  m_fixshape = fixshape;
  m_mcflag = false;
  RooDataSet* ds = GetDataSet();
  pdf_gen_sig->PrintParameters();
  pdf_gen_cmb->PrintParameters();
  pdf_gen_part->PrintParameters();
  Fit(ds);

  DefineElliRange();
  CalcSigIntegrals();
  DrawData(*pdf,ds);

  PrintIntegrals();

  CalcBinsFractions(ds);
  SaveCPVTree(ds);

  return;
}

void PurityFit::MakeFit2(const bool singlefbb,const bool fixshape){
  m_fixshape  = fixshape;
  m_singlefbb = singlefbb;
  m_mcflag = false;
  RooDataSet* ds = GetDataSet();
//  CalcBinsFractions(ds);
  Fit2(ds);

  DefineElliRange();
  CalcSigIntegrals();
  DrawData(*pdf,ds);
  PrintIntegrals();

//  SaveCPVTree(ds);
  return;
}

void PurityFit::Fit(RooDataSet* ds){
  cout << "Fit!" << endl;
  RooAbsPdf* pdf_sig   = pdf_gen_sig->GetPdf();
  RooAddPdf* pdf_comb  = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  pdf_sig->Print("t");
  pdf_comb->Print("t");
  pdf_part->Print("t");
  pdf_gen_sig->FixAll();
  pdf_gen_cmb->FixAll();
  pdf_gen_part->FixAll();

  if(!m_fixshape){
    if(true || m_mode != 10){
      pdf_gen_sig->Free_de0();
      pdf_gen_sig->Free_mbc0();
      pdf_gen_cmb->Release_de_QQ();
    }
    if(m_mode == 1 || m_mode == 10) pdf_gen_part->Release_de0();
  }

  Nsig = new RooRealVar("Nsig","Nsig",1150,-20,10000.);
  Ncmb = new RooRealVar("Ncmb","Ncmb",2288,0.,100000);
  fbb = pdf_gen_cmb->Get_fbb();
  f_p_f_bbc = new RooRealVar("f_p_f_bbc","f_p_f_bbc",cuts->f_p_f_bbc(m_mode,m_h0mode),0,10);
  if(m_mode != 1 && m_mode != 10){ f_p_f_bbc->setConstant(kTRUE);
  } else{
    f_p_f_bbc->setVal(0.5);
    f_p_f_bbc->setConstant(kFALSE);
  }
  Npart = new RooFormulaVar("Npart","Npart","@0*@1*@2",RooArgList(*Ncmb,*fbb,*f_p_f_bbc));

  pdf = new RooAddPdf("pdf","pdf",RooArgList(*pdf_sig,*pdf_part,*pdf_comb),RooArgList(*Nsig,*Npart,*Ncmb));

  Npart->Print("t");
  pdf->fitTo(*ds,Verbose(),Timer(true));

  return;
}

void PurityFit::Fit2(RooDataSet* ds){
  const double m_btau = 1.534;
  const double m_dm   = 0.510;
  const double m_xi   = 1+m_btau*m_dm*m_btau*m_dm;
  cout << "Fit2!" << endl;
  m_pdf_sig     = pdf_gen_sig->GetPdf();
  m_pdf_comb_bb = pdf_gen_cmb->GetPdfBB();
  m_pdf_comb_qq = pdf_gen_cmb->GetPdfQQ();
  m_pdf_part    = pdf_gen_part->GetPdf();
  f_p_f_bbc = new RooRealVar("f_p_f_bbc","f_p_f_bbc",cuts->f_p_f_bbc(m_mode,m_h0mode),0,10);
  if(m_mode != 1 && m_mode != 10){ f_p_f_bbc->setConstant(kTRUE);
  } else{
    f_p_f_bbc->setVal(0.5);
    f_p_f_bbc->setConstant(kFALSE);
  }
  fbb = new RooRealVar("fbb","fbb",0.5,0.,1.);
  if(!m_fixshape){
    if(true || m_mode != 10){
      pdf_gen_sig->Free_de0();
      pdf_gen_sig->Free_mbc0();
      pdf_gen_cmb->Release_de_QQ();
    }
    if(m_mode == 1 || m_mode == 10) pdf_gen_part->Release_de0();
  }

//  binflv = new RooSuperCategory("binflv","binflv",RooArgSet(*bin,*flv));
  Nsig = new RooRealVar("Nsig","Nsig",1150,-20,10000.);
  simpdf = new RooSimultaneous("simpdf","simpdf",*binflv);
  binflv->Print("t");
  stringstream out;
  const string flavorls[2] = {"anti-B0","B0"};
  cout << "Loop over bin x flv" << endl;
  for(int k=0; k<2; k++){
    for(int j=0; j<8;j++){
      const int flv = cuts->flv(k);
      const int bin = cuts->bin(j);
      const int bbind = cuts->bin_ind(-bin);
      cout << "  constants..." << endl;
      out.str(""); out << "wrtag_"<< k << "_" << j;
      wrtagI[k][j]   = new RooConstVar(out.str().c_str(),out.str().c_str(),WrTagMap[j]);
      out.str(""); out << "wrtag_"<< k << "_" << bbind;
      wrtagI[k][bbind] = new RooConstVar(out.str().c_str(),out.str().c_str(),WrTagMap[bbind]);
      out.str(""); out << "k_"<< k << "_" << j;
      KI[k][j]   = new RooConstVar(out.str().c_str(),out.str().c_str(),cuts->K( bin));
      out.str(""); out << "k_"<< k << "_" << bbind;
      KI[k][bbind] = new RooConstVar(out.str().c_str(),out.str().c_str(),cuts->K(-bin));
      out.str(""); out << "Nsig_"<< k << "_" << j;
      const string eq1("0.25*@5*((@0+@1)*@2-(1.-2.*@3)*@4*(@0-@1))/@2");
      RooConstVar* xi = new RooConstVar("xi","xi",m_xi);
      RooConstVar* curflv  = new RooConstVar("curflv","curflv",flv);
      cout << "  Nsig..." << endl;
      NsigI[k][j]   = new RooFormulaVar(out.str().c_str(),out.str().c_str(),eq1.c_str(),RooArgList(*KI[k][j],*KI[k][bbind],*xi,*wrtagI[k][j],*curflv,*Nsig));
      out.str(""); out << "Nsig_"<< k << "_" << bbind;
      NsigI[k][bbind] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),eq1.c_str(),RooArgList(*KI[k][bbind],*KI[k][j],*xi,*wrtagI[k][bbind],*curflv,*Nsig));
      cout << "  Ncmb..." << endl;
      out.str(""); out << "Ncmb_"<< k << "_" << j;
      NcmbI[k][j]   = new RooRealVar(out.str().c_str(),out.str().c_str(),100,0.,1e4);
      out.str(""); out << "Ncmb_"<< k << "_" << bbind;
      NcmbI[k][bbind] = new RooRealVar(out.str().c_str(),out.str().c_str(),100,0.,1e4);
      out.str(""); out << "fbb_"<< k << "_" << j;
      if(!m_singlefbb) fbbI[k][j]   = new RooRealVar(out.str().c_str(),out.str().c_str(),0.5,0.,1.);
      else           fbbI[k][j] = fbb;
      out.str(""); out << "fbb_"<< k << "_" << bbind;
      if(!m_singlefbb) fbbI[k][bbind] = new RooRealVar(out.str().c_str(),out.str().c_str(),0.5,0.,1.);
      else           fbbI[k][bbind] = fbb;
      out.str(""); out << "pdf_comb_" << k << "_" << j;
      m_pdf_comb[k][j] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*m_pdf_comb_bb,*m_pdf_comb_qq),RooArgList(*fbbI[k][j]));
      out.str(""); out << "pdf_comb_" << k << "_" << bbind;
      m_pdf_comb[k][bbind] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*m_pdf_comb_bb,*m_pdf_comb_qq),RooArgList(*fbbI[k][bbind]));
      out.str(""); out << "Nprt_" << k << "_" << j;
      cout << "  Nprt..." << endl;
      NprtI[k][j]   = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"@0*@1*@2",RooArgList(*NcmbI[k][j],*fbbI[k][j],*f_p_f_bbc));
      out.str(""); out << "Nprt_"<< k << "_" << bbind;
      NprtI[k][bbind] = new RooFormulaVar(out.str().c_str(),out.str().c_str(),"@0*@1*@2",RooArgList(*NcmbI[k][bbind],*fbbI[k][bbind],*f_p_f_bbc));
      out.str(""); out << "pdf_"<< k << "_" << j;
      cout << "  pdfs..." << endl;
      pdfI[k][j]   = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*m_pdf_sig,*m_pdf_part,*m_pdf_comb[k][j]),RooArgList(*NsigI[k][j],*NprtI[k][j],*NcmbI[k][j]));
      out.str(""); out << "pdf_"<< k << "_" << bbind;
      pdfI[k][bbind] = new RooAddPdf(out.str().c_str(),out.str().c_str(),RooArgList(*m_pdf_sig,*m_pdf_part,*m_pdf_comb[k][bbind]),RooArgList(*NsigI[k][bbind],*NprtI[k][bbind],*NcmbI[k][bbind]));

      cout << "  add pdfs to SimPdf" << endl;
      out.str(""); out << "{" <<  bin << ";" << flavorls[k] << "}";
      simpdf->addPdf(*pdfI[k][j],out.str().c_str());
      out.str(""); out << "{" << -bin << ";" << flavorls[k] << "}";
      simpdf->addPdf(*pdfI[k][bbind],out.str().c_str());
    }
  }
  simpdf->fitTo(*ds,Verbose(),Timer(true));
}

void PurityFit::GenerateDS(const int NTot){
  de_mbc_ds = const_cast< RooDataSet* > (pdf->generate(RooArgSet(*de,*mbc),NTot));
  de_mbc_ds->Print();
}

void PurityFit::CalcSigIntegrals(void){
  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();

  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  RooAbsReal* intSigElli  = pdf_sig->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  sigint = intSigElli->getVal();
  nsigEl = sigint*Nsig->getVal();
  const double nsig_errEl = sigint*Nsig->getError();
  const double nsig_errEl_npq = TMath::Sqrt(fabs(nsigEl*(Nsig->getVal()-nsigEl)/Nsig->getVal()));
  nsig_errEl_total = TMath::Sqrt(fabs(nsig_errEl*nsig_errEl+nsig_errEl_npq*nsig_errEl_npq));

  RooAbsReal* intPartEl = pdf_part->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  partint = intPartEl->getVal();
  const double npart = f_p_f_bbc->getVal()*Ncmb->getVal()*fbb->getVal();
  npartEl = partint*npart;
  const double npart_errEl     = partint*fbb->getError()*Ncmb->getVal()*cuts->f_p_f_bbc(m_mode,m_h0mode);
  const double npart_errEl_npq = TMath::Sqrt(fabs(npartEl*(npart-npartEl)/npart));
  npart_errEl_total = TMath::Sqrt(fabs(npart_errEl*npart_errEl+npart_errEl_npq*npart_errEl_npq));

  RooAbsReal* intCombEl = pdf_comb->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  cmbint = intCombEl->getVal();
  ncmbEl = cmbint*Ncmb->getVal();
  const double ncmb_errEl     = cmbint*Ncmb->getError();
  const double ncmb_errEl_npq = TMath::Sqrt(ncmbEl*(Ncmb->getVal()-ncmbEl)/Ncmb->getVal());
  ncmb_errEl_total = TMath::Sqrt(ncmb_errEl*ncmb_errEl+ncmb_errEl_npq*ncmb_errEl_npq);

  RooAbsReal* intContEl = pdf_cont->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  const double Ncnt = Ncmb->getVal()*f_cont_in_comb;
  cntint = intContEl->getVal();
  ncntEl = cntint*Ncmb->getVal()*f_cont_in_comb;
  const double ncnt_errEl = ncntEl*Ncmb->getError();
  const double ncnt_errEl_npq = TMath::Sqrt(ncntEl*(Ncnt-ncntEl)/Ncnt);
  ncnt_errEl_total = TMath::Sqrt(ncnt_errEl*ncnt_errEl+ncnt_errEl_npq*ncnt_errEl_npq);

  purityEl = nsigEl/(nsigEl+npartEl+ncmbEl)*100;
  purity_errEl = nsig_errEl_total/(nsigEl+npartEl+ncmbEl)*100;
  return;
}

void PurityFit::CalcSigIntegrals2(void){
  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();

  RooAbsReal* intSigElli = pdf_sig->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  RooAbsReal* intPrtElli = pdf_part->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  RooAbsReal* intCmbElli = pdf_comb->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  RooAbsReal* intCntElli = pdf_cont->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Elli"));
  sigint  = intSigElli->getVal();
  partint = intPrtElli->getVal();
  cmbint  = intCmbElli->getVal();
  cntint  = intCntElli->getVal();

  const double SigRelErr = Nsig->getError()/Nsig->getVal();
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      nsigElI[k][j] = sigint*NsigI[k][j]->getVal();
      const double nsig_errEl = nsigElI[k][j]*SigRelErr;
      const double nsig_errEl_npq = TMath::Sqrt(fabs(nsigElI[k][j]*(Nsig->getVal()-nsigElI[k][j])/Nsig->getVal()));
      nsig_errEl_totalI[k][j] = TMath::Sqrt(fabs(nsig_errEl*nsig_errEl+nsig_errEl_npq*nsig_errEl_npq));

      ncmbElI[k][j] = cmbint*NcmbI[k][j]->getVal();
      const double ncmb_errEl     = cmbint*NcmbI[k][j]->getError();
      const double ncmb_errEl_npq = TMath::Sqrt(ncmbElI[k][j]*(NcmbI[k][j]->getVal()-ncmbElI[k][j])/NcmbI[k][j]->getVal());
      ncmb_errEl_totalI[k][j]  = TMath::Sqrt(ncmb_errEl*ncmb_errEl+ncmb_errEl_npq*ncmb_errEl_npq);

      const double fcnt = 1.-fbbI[k][j]->getVal();
      ncntElI[k][j] = fcnt*cntint*NcmbI[k][j]->getVal();
      const double ncnt_full  = fcnt*NcmbI[k][j]->getVal();
      const double ncnt_errEl = fcnt*cntint*NcmbI[k][j]->getError();
      const double ncnt_errEl_npq = TMath::Sqrt(ncntElI[k][j]*(ncnt_full-ncntElI[k][j])/ncnt_full);
      ncnt_errEl_totalI[k][j] = TMath::Sqrt(ncnt_errEl*ncnt_errEl+ncnt_errEl_npq*ncnt_errEl_npq);

      const double fprt = fbbI[k][j]->getVal()*f_p_f_bbc->getVal();
      nprtElI[k][j] = fprt*partint*NprtI[k][j]->getVal();
      const double nprt_full  = fprt*NcmbI[k][j]->getVal();
      const double nprt_errEl = fcnt*partint*NcmbI[k][j]->getError();
      const double nprt_errEl_npq = TMath::Sqrt(nprtElI[k][j]*(nprt_full-nprtElI[k][j])/nprt_full);
      nprt_errEl_totalI[k][j] = TMath::Sqrt(nprt_errEl*nprt_errEl+nprt_errEl_npq*nprt_errEl_npq);

      purElI[k][j] = nsigElI[k][j]/(nsigElI[k][j]+nprtElI[k][j]+ncmbElI[k][j])*100;
      pur_errElI[k][j] = nsig_errEl_totalI[k][j]/(nsigElI[k][j]+nprtElI[k][j]+ncmbElI[k][j])*100;
    }
  }

  return;
}

void PurityFit::CalcSB1Integrals(void){
  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  RooAbsReal* intSig_sb1  = pdf_sig->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband1"));
  sigint_sb1 = intSig_sb1->getVal();
  nsig_sb1 = sigint_sb1*Nsig->getVal();
  const double nsig_err_sb1 = sigint_sb1*Nsig->getError();
  const double nsig_err_sb1_npq = TMath::Sqrt(fabs(nsig_sb1*(Nsig->getVal()-nsig_sb1)/Nsig->getVal()));
  nsig_err_sb1_total = TMath::Sqrt(fabs(nsig_err_sb1*nsig_err_sb1+nsig_err_sb1_npq*nsig_err_sb1_npq));

  RooAbsReal* intPart_sb1 = pdf_part->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband1"));
  partint_sb1 = intPart_sb1->getVal();
  const double npart = f_p_f_bbc->getVal()*Ncmb->getVal()*fbb->getVal();
  npart_sb1 = partint_sb1*npart;
  const double npart_err_sb1     = partint*fbb->getError()*Ncmb->getVal()*cuts->f_p_f_bbc(m_mode,m_h0mode);
  const double npart_err_sb1_npq = TMath::Sqrt(fabs(npart_sb1*(npart-npart_sb1)/npart));
  npart_err_sb1_total = TMath::Sqrt(fabs(npart_err_sb1*npart_err_sb1+npart_err_sb1_npq*npart_err_sb1_npq));

  RooAbsReal* intComb_sb1 = pdf_comb->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband1"));
  cmbint_sb1 = intComb_sb1->getVal();
  ncmb_sb1 = cmbint_sb1*Ncmb->getVal();
  const double ncmb_err_sb1     = cmbint*Ncmb->getError();
  const double ncmb_err_sb1_npq = TMath::Sqrt(ncmb_sb1*(Ncmb->getVal()-ncmb_sb1)/Ncmb->getVal());
  ncmb_err_sb1_total = TMath::Sqrt(ncmb_err_sb1*ncmb_err_sb1+ncmb_err_sb1_npq*ncmb_err_sb1_npq);

  RooAbsReal* intCont_sb1 = pdf_cont->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband1"));
  const double Ncnt_sb1 = Ncmb->getVal()*f_cont_in_comb;
  cntint_sb1 = intCont_sb1->getVal();
  ncnt_sb1 = cntint_sb1*Ncmb->getVal()*f_cont_in_comb;
  const double ncnt_err_sb1 = ncnt_sb1*Ncmb->getError();
  const double ncnt_err_sb1_npq = TMath::Sqrt(ncnt_sb1*(Ncnt_sb1-ncnt_sb1)/Ncnt_sb1);
  ncnt_err_sb1_total = TMath::Sqrt(ncnt_err_sb1*ncnt_err_sb1+ncnt_err_sb1_npq*ncnt_err_sb1_npq);

  purity_sb1 = nsig_sb1/(nsig_sb1+npart_sb1+ncmb_sb1)*100;
  purity_err_sb1 = nsig_err_sb1_total/(nsig_sb1+npart_sb1+ncmb_sb1)*100;
  return;
}

void PurityFit::CalcSB2Integrals(void){
  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  RooAbsReal* intSig_sb2  = pdf_sig->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband2"));
  sigint_sb2 = intSig_sb2->getVal();
  nsig_sb2 = sigint_sb2*Nsig->getVal();
  const double nsig_err_sb2 = sigint_sb2*Nsig->getError();
  const double nsig_err_sb2_npq = TMath::Sqrt(fabs(nsig_sb2*(Nsig->getVal()-nsig_sb2)/Nsig->getVal()));
  nsig_err_sb2_total = TMath::Sqrt(fabs(nsig_err_sb2*nsig_err_sb2+nsig_err_sb2_npq*nsig_err_sb2_npq));

  RooAbsReal* intPart_sb2 = pdf_part->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband2"));
  partint_sb2 = intPart_sb2->getVal();
  const double npart = f_p_f_bbc->getVal()*Ncmb->getVal()*fbb->getVal();
  npart_sb2 = partint_sb2*npart;
  const double npart_err_sb2     = partint*fbb->getError()*Ncmb->getVal()*cuts->f_p_f_bbc(m_mode,m_h0mode);
  const double npart_err_sb2_npq = TMath::Sqrt(fabs(npart_sb2*(npart-npart_sb2)/npart));
  npart_err_sb2_total = TMath::Sqrt(fabs(npart_err_sb2*npart_err_sb2+npart_err_sb2_npq*npart_err_sb2_npq));

  RooAbsReal* intComb_sb2 = pdf_comb->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband2"));
  cmbint_sb2 = intComb_sb2->getVal();
  ncmb_sb2 = cmbint_sb2*Ncmb->getVal();
  const double ncmb_err_sb2     = cmbint*Ncmb->getError();
  const double ncmb_err_sb2_npq = TMath::Sqrt(ncmb_sb2*(Ncmb->getVal()-ncmb_sb2)/Ncmb->getVal());
  ncmb_err_sb2_total = TMath::Sqrt(ncmb_err_sb2*ncmb_err_sb2+ncmb_err_sb2_npq*ncmb_err_sb2_npq);

  RooAbsReal* intCont_sb2 = pdf_cont->createIntegral(RooArgSet(*de,*mbc),NormSet(RooArgSet(*de,*mbc)),Range("Sideband2"));
  const double Ncnt_sb2 = Ncmb->getVal()*f_cont_in_comb;
  cntint_sb2 = intCont_sb2->getVal();
  ncnt_sb2 = cntint_sb2*Ncmb->getVal()*f_cont_in_comb;
  const double ncnt_err_sb2 = ncnt_sb2*Ncmb->getError();
  const double ncnt_err_sb2_npq = TMath::Sqrt(ncnt_sb2*(Ncnt_sb2-ncnt_sb2)/Ncnt_sb2);
  ncnt_err_sb2_total = TMath::Sqrt(ncnt_err_sb2*ncnt_err_sb2+ncnt_err_sb2_npq*ncnt_err_sb2_npq);

  purity_sb2 = nsig_sb2/(nsig_sb2+npart_sb2+ncmb_sb2)*100;
  purity_err_sb2 = nsig_err_sb2_total/(nsig_sb2+npart_sb2+ncmb_sb2)*100;
  return;
}

void PurityFit::PrintIntegrals(void){
  cout << " *** Fit results for signal region ***" << endl;
  cout << "*               *  Signal range  *" << endl;
  cout << "1. Signal:        " << nsigEl  << " +- " << nsig_errEl_total << endl;
  cout << "2. Combinatorial: " << ncmbEl  << " +- " << ncmb_errEl_total << endl;
  cout << "3. Partial:       " << npartEl << " +- " << npart_errEl_total << endl;
  cout << "Purity: " << purityEl << " +- " << purity_errEl << endl;
  return;
}

void PurityFit::WriteIntegrals(void){
  ofstream ofile;
  stringstream out;
  out.str("");
  out << "params/Integrals_m" << m_mode << "_mh0" << m_h0mode;
  if(m_fixshape)  out << "_fixedshape";
  if(m_singlefbb) out << "_singlefbb";
  if(m_mcflag)    out << "_mc";
  out << ".txt";
  ofile.open(out.str().c_str(),ofstream::out);
  ofile << "Fit results for fit region" << endl;
  ofile << "1. Signal:        " << Nsig->getVal()  << " +- " << Nsig->getError() << endl;
  ofile << "2. Combinatorial: " << Ncmb->getVal()  << " +- " << Ncmb->getError() << endl;
  const double ncnt = Ncmb->getVal()*(1.-fbb->getVal());
  const double ncnt_rerr = ncnt*(Ncmb->getError()/Ncmb->getVal() + fbb->getError()/(1.-fbb->getVal()));
  cout << Ncmb->getError() << " " << Ncmb->getVal() << " " << fbb->getError() << " " << 1.-fbb->getVal() << endl;
  ofile << "2.1 Continuum:    " << ncnt  << " +- " << ncnt_rerr << endl;
  ofile << "3. Partial:       " << Npart->getVal() << " +- " << 0 << endl;
  const double purity = Nsig->getVal()/(Nsig->getVal()+Ncmb->getVal()+Npart->getVal())*100.;
  ofile << "Purity: " << purity << " +- " << 0 << endl;
  ofile << endl;
  ofile << "Fit results for signal region" << endl;
  ofile << "1. Signal:        " << nsigEl  << " +- " << nsig_errEl_total << endl;
  ofile << "2. Combinatorial: " << ncmbEl  << " +- " << ncmb_errEl_total << endl;
  ofile << "2.1 Continuum:    " << ncntEl  << " +- " << ncnt_errEl_total << endl;
  ofile << "3. Partial:       " << npartEl << " +- " << npart_errEl_total << endl;
  ofile << "Purity: " << purityEl << " +- " << purity_errEl << endl;
  ofile << endl;
  ofile << "Fit results for sideband 1" << endl;
  ofile << "1. Signal:        " << nsig_sb1  << " +- " << nsig_err_sb1_total << endl;
  ofile << "2. Combinatorial: " << ncmb_sb1  << " +- " << ncmb_err_sb1_total << endl;
  ofile << "2.1 Continuum:    " << ncnt_sb1  << " +- " << ncnt_err_sb1_total << endl;
  ofile << "3. Partial:       " << npart_sb1 << " +- " << npart_err_sb1_total << endl;
  ofile << "Purity: " << purity_sb1 << " +- " << purity_err_sb1 << endl;
  ofile << endl;
  ofile << "Fit results for sideband 2" << endl;
  ofile << "1. Signal:        " << nsig_sb2  << " +- " << nsig_err_sb2_total << endl;
  ofile << "2. Combinatorial: " << ncmb_sb2  << " +- " << ncmb_err_sb2_total << endl;
  ofile << "2.1 Continuum:    " << ncnt_sb2  << " +- " << ncnt_err_sb2_total << endl;
  ofile << "3. Partial:       " << npart_sb2 << " +- " << npart_err_sb2_total << endl;
  ofile << "Purity: " << purity_sb2 << " +- " << purity_err_sb2 << endl;

  ofile.close();
  return;
}

void PurityFit::WriteIntegrals2(void){
  ofstream ofile;
  stringstream out;
  out.str("");
  out << "params/Integrals2_m" << m_mode << "_mh0" << m_h0mode;
  if(m_fixshape)  out << "_fixedshape";
  if(m_singlefbb) out << "_singlefbb";
  if(m_mcflag)    out << "_mc";
  out << ".txt";
  ofile.open(out.str().c_str(),ofstream::out);
//  ofile << "Fit results for fit region" << endl;
//  ofile << "1. Signal:        " << Nsig->getVal()  << " +- " << Nsig->getError() << endl;
//  ofile << "2. Combinatorial: " << Ncmb->getVal()  << " +- " << Ncmb->getError() << endl;
//  const double ncnt = Ncmb->getVal()*(1.-fbb->getVal());
//  const double ncnt_rerr = ncnt*(Ncmb->getError()/Ncmb->getVal() + fbb->getError()/(1.-fbb->getVal()));
//  cout << Ncmb->getError() << " " << Ncmb->getVal() << " " << fbb->getError() << " " << 1.-fbb->getVal() << endl;
//  ofile << "2.1 Continuum:    " << ncnt  << " +- " << ncnt_rerr << endl;
//  ofile << "3. Partial:       " << Npart->getVal() << " +- " << 0 << endl;
//  const double purity = Nsig->getVal()/(Nsig->getVal()+Ncmb->getVal()+Npart->getVal())*100.;
//  ofile << "Purity: " << purity << " +- " << 0 << endl;
//  ofile << endl;
  ofile << "Fit results for signal region" << endl;
  for(int k=0; k<2; k++){
    const int flvav = cuts->flv(k);
    for(int j=0; j<16; j++){
      const int binnum = cuts->bin(j);
      ofile << flvav << " " << binnum;
      ofile << " " << nsigElI[k][j] << " +- " << nsig_errEl_totalI[k][j];
      ofile << " " << ncmbElI[k][j] << " +- " << ncmb_errEl_totalI[k][j];
      ofile << " " << ncntElI[k][j] << " +- " << ncnt_errEl_totalI[k][j];
      ofile << " " << nprtElI[k][j] << " +- " << nprt_errEl_totalI[k][j];
      ofile << " " << purElI[k][j]  << " +- " << pur_errElI[k][j];
      ofile << endl;
    }
  }
//  ofile << "1. Signal:        " << nsigEl  << " +- " << nsig_errEl_total << endl;
//  ofile << "2. Combinatorial: " << ncmbEl  << " +- " << ncmb_errEl_total << endl;
//  ofile << "2.1 Continuum:    " << ncntEl  << " +- " << ncnt_errEl_total << endl;
//  ofile << "3. Partial:       " << npartEl << " +- " << npart_errEl_total << endl;
//  ofile << "Purity: " << purityEl << " +- " << purity_errEl << endl;
//  ofile << endl;
//  ofile << "Fit results for sideband 1" << endl;
//  ofile << "1. Signal:        " << nsig_sb1  << " +- " << nsig_err_sb1_total << endl;
//  ofile << "2. Combinatorial: " << ncmb_sb1  << " +- " << ncmb_err_sb1_total << endl;
//  ofile << "2.1 Continuum:    " << ncnt_sb1  << " +- " << ncnt_err_sb1_total << endl;
//  ofile << "3. Partial:       " << npart_sb1 << " +- " << npart_err_sb1_total << endl;
//  ofile << "Purity: " << purity_sb1 << " +- " << purity_err_sb1 << endl;
//  ofile << endl;
//  ofile << "Fit results for sideband 2" << endl;
//  ofile << "1. Signal:        " << nsig_sb2  << " +- " << nsig_err_sb2_total << endl;
//  ofile << "2. Combinatorial: " << ncmb_sb2  << " +- " << ncmb_err_sb2_total << endl;
//  ofile << "2.1 Continuum:    " << ncnt_sb2  << " +- " << ncnt_err_sb2_total << endl;
//  ofile << "3. Partial:       " << npart_sb2 << " +- " << npart_err_sb2_total << endl;
//  ofile << "Purity: " << purity_sb2 << " +- " << purity_err_sb2 << endl;

  ofile.close();
  return;
}

void PurityFit::CountTrueNumbers(RooDataSet* ds,RooDataSet* cntds){
  stringstream out;
  out.str("");
  out << "(de-" << de_center->getVal() << ")/" << de_radius->getVal() << "*(de-" << de_center->getVal() << ")/" << de_radius->getVal() << "+(mbc-"<<mbc_center->getVal()<<")/" << mbc_radius->getVal() << "*(mbc-" << mbc_center->getVal() << ")/" << mbc_radius->getVal() << "<1";
  const string ellicut = out.str();
  const string sb1cut("mbc>5.23 && mbc<5.25");
  const string sb2cut("mbc>5.25 && de>0.1");

  int rndm_pi0      = 0;
  int rndm_pi0_elli = 0;
  int rndm_pi0_sb1  = 0;
  int rndm_pi0_sb2  = 0;

  const string rndm_pi0_cut("rndm_pi0 == 1");
  const string strand(" && ");
  const string rndm_pi0_cut_elli = rndm_pi0_cut + strand + ellicut;
  const string rndm_pi0_cut_sb1  = rndm_pi0_cut + strand + sb1cut;
  const string rndm_pi0_cut_sb2  = rndm_pi0_cut + strand + sb2cut;

  rndm_pi0      = ds->sumEntries(rndm_pi0_cut.c_str());
  rndm_pi0_elli = ds->sumEntries(rndm_pi0_cut_elli.c_str());
  rndm_pi0_sb1  = ds->sumEntries(rndm_pi0_cut_sb1.c_str());
  rndm_pi0_sb2  = ds->sumEntries(rndm_pi0_cut_sb2.c_str());
  cout << "Random soft pi0:" << endl;
  cout << "  Fit area: " << rndm_pi0      << endl;
  cout << "  Sig area: " << rndm_pi0_elli << endl;
  cout << "  Sb1 area: " << rndm_pi0_sb1  << endl;
  cout << "  Sb2 area: " << rndm_pi0_sb2  << endl;

  Roo1DTable* fulltable = ds->table(*b0f);
  Roo1DTable* cnt_fulltable = cntds->table(*b0f);
  NSIGNAL = fulltable->get("signal") + fulltable->get("fsr") + fulltable->get("bad_pi0") + rndm_pi0;
  NCOMB   = fulltable->get("comb") - rndm_pi0;
  NCONT   = cnt_fulltable->get("comb");
  NPART   = fulltable->get("rho2") + fulltable->get("rho3") + fulltable->get("rho4") + fulltable->get("rho11");

  Roo1DTable* ellitable = ds->table(*b0f,ellicut.c_str());
  Roo1DTable* cnt_ellitable = cntds->table(*b0f,ellicut.c_str());
  NSIGNAL_ELLI = ellitable->get("signal") + ellitable->get("fsr") + ellitable->get("bad_pi0") + rndm_pi0_elli;
  NCOMB_ELLI   = ellitable->get("comb") - rndm_pi0_elli;
  NCONT_ELLI   = cnt_ellitable->get("comb");
  NPART_ELLI   = ellitable->get("rho2") + ellitable->get("rho3") + ellitable->get("rho4") + ellitable->get("rho11");
  PURITY = (double)NSIGNAL_ELLI/(NSIGNAL_ELLI+NCOMB_ELLI+NPART_ELLI)*100;

  Roo1DTable* sb1_table = ds->table(*b0f,sb1cut.c_str());
  Roo1DTable* cnt_sb1_table = cntds->table(*b0f,sb1cut.c_str());
  NSIGNAL_SB1 = sb1_table->get("signal") + sb1_table->get("fsr") + sb1_table->get("bad_pi0") + rndm_pi0_sb1;
  NCOMB_SB1   = sb1_table->get("comb") - rndm_pi0_sb1;
  NCONT_SB1   = cnt_sb1_table->get("comb");
  NPART_SB1   = sb1_table->get("rho2") + sb1_table->get("rho3") + sb1_table->get("rho4") + sb1_table->get("rho11");
  PURITY_SB1  = (double)NSIGNAL_SB1/(NSIGNAL_SB1+NCOMB_SB1+NPART_SB1);

  Roo1DTable* sb2_table = ds->table(*b0f,sb2cut.c_str());
  Roo1DTable* cnt_sb2_table = cntds->table(*b0f,sb2cut.c_str());
  NSIGNAL_SB2 = sb2_table->get("signal") + sb2_table->get("fsr") + sb2_table->get("bad_pi0") + rndm_pi0_sb2;
  NCOMB_SB2   = sb2_table->get("comb") - rndm_pi0_sb2;
  NCONT_SB2   = cnt_sb2_table->get("comb");
  NPART_SB2   = sb2_table->get("rho2") + sb2_table->get("rho3") + sb2_table->get("rho4") + sb2_table->get("rho11");
  PURITY_SB2  = (double)NSIGNAL_SB2/(NSIGNAL_SB2+NCOMB_SB2+NPART_SB2);

  RooArgSet b0fbinflv;
  b0fbinflv.add(*b0f);
  b0fbinflv.add(*bin);
  b0fbinflv.add(*flv);
  Roo1DTable* fulltable_binflv          = ds->table(b0fbinflv);
  Roo1DTable* cnt_fulltable_binflv      = cntds->table(b0fbinflv);
  Roo1DTable* fulltable_binflv_elli     = ds->table(b0fbinflv,ellicut.c_str());
  Roo1DTable* cnt_fulltable_binflv_elli = cntds->table(b0fbinflv,ellicut.c_str());

  const string flavors[2] = {"anti-B0","B0"};
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      out.str(""); out << "{signal;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NSIGARR[k][j]      = fulltable_binflv->get(out.str().c_str());
      NSIGARR_ELLI[k][j] = fulltable_binflv_elli->get(out.str().c_str());
      out.str(""); out << "{comb;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NCMBARR[k][j]      = fulltable_binflv->get(out.str().c_str());
      NCMBARR_ELLI[k][j] = fulltable_binflv_elli->get(out.str().c_str());
      NCNTARR[k][j]      = cnt_fulltable_binflv->get(out.str().c_str());
      NCNTARR_ELLI[k][j] = cnt_fulltable_binflv_elli->get(out.str().c_str());
      out.str(""); out << "{rho2;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NPRTARR[k][j]      = fulltable_binflv->get(out.str().c_str());
      NPRTARR_ELLI[k][j] = fulltable_binflv_elli->get(out.str().c_str());
      out.str(""); out << "{rho3;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NPRTARR[k][j]     += fulltable_binflv->get(out.str().c_str());
      NPRTARR_ELLI[k][j]+= fulltable_binflv_elli->get(out.str().c_str());
      out.str(""); out << "{rho4;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NPRTARR[k][j]     += fulltable_binflv->get(out.str().c_str());
      NPRTARR_ELLI[k][j]+= fulltable_binflv_elli->get(out.str().c_str());
      out.str(""); out << "{rho11;" << cuts->bin(j) << ";" << flavors[k] << "}";
      NPRTARR[k][j]     += fulltable_binflv->get(out.str().c_str());
      NPRTARR_ELLI[k][j]+= fulltable_binflv_elli->get(out.str().c_str());
    }
  }
  cout << "True numbers recorded!" << endl;
  return;
}

void PurityFit::PrintTrueNumbers(void){
  cout << " *** Data Set Structure ***" << endl;
  cout << "*               *  Signal range  *  Full range  *  Fraction " << endl;
  cout << "1. Signal:        " << NSIGNAL_ELLI << "             " << NSIGNAL << "            " << (double)NSIGNAL_ELLI/NSIGNAL << endl;
  cout << "2. Combinatorial: " <<   NCOMB_ELLI << "             " <<   NCOMB << "            " << (double)NCOMB_ELLI/NCOMB     << endl;
  cout << "2.1 Continuum:    " <<   NCONT_ELLI << "             " <<   NCONT << "            " << (double)NCONT_ELLI/NCONT     << endl;
  cout << "3. Partial:       " <<   NPART_ELLI << "             " <<   NPART << "            " << (double)NPART_ELLI/NPART     << endl;
  cout << "Purity: " << PURITY << endl;
  cout << endl;
  cout << "For each bin and flavor:" << endl;
  cout << "  Full range:" << endl;
  const string flavors[2] = {"anti-B0","B0"};
  for(int k=0; k<2; k++){
    cout << "   " << flavors[k] << ":" << endl;
    for(int j=0; j<16; j++){
      cout << "bin " << cuts->bin(j) << ": " << NSIGARR[k][j];
      cout << " " << NCMBARR[k][j];
      cout << " " << NCNTARR[k][j];
      cout << " " << NPRTARR[k][j];
      cout << endl;
    }
  }
  cout << "  Signal range:" << endl;
  for(int k=0; k<2; k++){
    cout << "   " << flavors[k] << ":" << endl;
    for(int j=0; j<16; j++){
      cout << "bin " << cuts->bin(j) << ": " << NSIGARR_ELLI[k][j];
      cout << " " << NCMBARR_ELLI[k][j];
      cout << " " << NCNTARR_ELLI[k][j];
      cout << " " << NPRTARR_ELLI[k][j];
      cout << endl;
    }
  }
  return;
}

void PurityFit::WriteTrueNumbers(void){
  ofstream ofile;
  stringstream out;
  out.str("");
  out << "params/TrueNumbers_m" << m_mode << "_mh0" << m_h0mode << ".txt";
  ofile.open(out.str().c_str(),ofstream::out);

  ofile << "Data Set Structure. Mode " << m_mode << ", h0 mode " << m_h0mode << endl;
  ofile << "*               * Total * Sig * Frac * SB1 * FrSB1 * SB2 * FrSB2 *" << endl;
  ofile << "1. Signal:        " << NSIGNAL << " " << NSIGNAL_ELLI << " " << (double)NSIGNAL_ELLI/NSIGNAL;
  ofile << " " << NSIGNAL_SB1 << " " << (double)NSIGNAL_SB1/NSIGNAL;
  ofile << " " << NSIGNAL_SB2 << " " << (double)NSIGNAL_SB2/NSIGNAL << endl;
  ofile << "2. Combinatorial: " <<   NCOMB << " " <<   NCOMB_ELLI << " " << (double)NCOMB_ELLI/NCOMB;
  ofile << " " << NCOMB_SB1 << " " << (double)NCOMB_SB1/NCOMB;
  ofile << " " << NCOMB_SB2 << " " << (double)NCOMB_SB2/NCOMB << endl;
  ofile << "2.1 Continuum:    " <<   NCONT << " " <<   NCONT_ELLI << " " << (double)NCONT_ELLI/NCONT;
  ofile << " " << NCONT_SB1 << " " << (double)NCONT_SB1/NCONT;
  ofile << " " << NCONT_SB2 << " " << (double)NCONT_SB2/NCONT << endl;
  ofile << "3. Partial:       " <<   NPART << " " <<   NPART_ELLI << " " << (double)NPART_ELLI/NPART;
  ofile << " " << NPART_SB1 << " " << (double)NPART_SB1/NPART;
  ofile << " " << NPART_SB2 << " " << (double)NPART_SB2/NPART << endl;
  ofile << "Purity: " << PURITY << " " << PURITY_SB1 << " " << PURITY_SB2 << endl;
  ofile << endl;
  ofile << "For each bin and flavor:" << endl;
  ofile << "  Full range:" << endl;
  const string flavors[2] = {"anti-B0","B0"};
  for(int k=0; k<2; k++){
    ofile << "   " << flavors[k] << ":" << endl;
    for(int j=0; j<16; j++){
      ofile << "bin " << cuts->bin(j) << ": " << NSIGARR[k][j];
      ofile << " " << NCMBARR[k][j];
      ofile << " " << NCNTARR[k][j];
      ofile << " " << NPRTARR[k][j];
      ofile << endl;
    }
  }
  ofile << "  Signal range:" << endl;
  for(int k=0; k<2; k++){
    ofile << "   " << flavors[k] << ":" << endl;
    for(int j=0; j<16; j++){
      ofile << "bin " << cuts->bin(j) << ": " << NSIGARR_ELLI[k][j];
      ofile << " " << NCMBARR_ELLI[k][j];
      ofile << " " << NCNTARR_ELLI[k][j];
      ofile << " " << NPRTARR_ELLI[k][j];
      ofile << endl;
    }
  }

  ofile.close();
  return;
}


void PurityFit::DefineElliRange(void){
 mbc_center = new RooRealVar("mbc_center","mbc_center",0.5*(mbc_min+mbc_max)); mbc_center->setConstant(kTRUE);
 mbc_radius = new RooRealVar("mbc_radius","mbc_radius",0.5*(mbc_max-mbc_min)); mbc_radius->setConstant(kTRUE);
 de_center  = new RooRealVar("de_center","de_center",0.5*(de_min+de_max));     de_center->setConstant(kTRUE);
 de_radius  = new RooRealVar("de_radius","de_radius",0.5*(de_max-de_min));     de_radius->setConstant(kTRUE);

 de->setRange("Elli",de_min,de_max);
 RooFormulaVar* mbclo = new RooFormulaVar("mbclo","@1-@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(*de,*mbc_center,*mbc_radius,*de_center,*de_radius));
 RooFormulaVar* mbchi = new RooFormulaVar("mbchi","@1+@2*TMath::Sqrt(TMath::Abs(1-(@0-@3)/@4*(@0-@3)/@4)+0.0000001)",RooArgSet(*de,*mbc_center,*mbc_radius,*de_center,*de_radius));
 mbc->setRange("Elli",*mbclo,*mbchi);
}

void PurityFit::DrawDeltaE(RooAddPdf* pdf,RooDataSet* ds){
  RooPlot* deFrame = de->frame();
  RooProdPdf* pdf_cmb_bb = pdf_gen_cmb->GetPdfBB();
  RooProdPdf* pdf_cmb_qq = pdf_gen_cmb->GetPdfQQ();
  RooAbsPdf*  pdf_sig    = pdf_gen_sig->GetPdf();
  RooProdPdf* pdf_part   = pdf_gen_part->GetPdf();
  ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
//  if(m_mcflag){
//    ds->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"));
//  }
  pdf->plotOn(deFrame,Components(*pdf_sig),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf->plotOn(deFrame,Components(*pdf_part),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf->plotOn(deFrame,Components(*pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf->plotOn(deFrame,Components(*pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("mbcSignal"));
  pdf->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));

  RooHist* hdepull = deFrame->pullHist();
  RooPlot* dePull  = de->frame(Title("#Delta E pull distribution"));
  dePull->addPlotable(hdepull,"P");
  dePull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cm = new TCanvas("#Delta E, Signal","#Delta E, Signal",600,700);
  cm->cd();

  TPad *pad3 = new TPad("pad3","pad3",0.01,0.20,0.99,0.99);
  TPad *pad4 = new TPad("pad4","pad4",0.01,0.01,0.99,0.20);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  pad3->SetLeftMargin(0.15);
  pad3->SetFillColor(0);

  deFrame->GetXaxis()->SetTitleSize(0.05);
  deFrame->GetXaxis()->SetTitleOffset(0.85);
  deFrame->GetXaxis()->SetLabelSize(0.04);
  deFrame->GetYaxis()->SetTitleOffset(1.6);
  deFrame->Draw();

  stringstream out;

  TPaveText *pt = new TPaveText(0.5,0.6,0.98,0.9,"brNDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  out.str("");
  out << "#chi^{2}/n.d.f = " << deFrame->chiSquare();
  pt->AddText(out.str().c_str());
  out.str("");
  out << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  if(m_mcflag) out << " (" << NSIGNAL_ELLI << ")";
  pt->AddText(out.str().c_str());
  out.str("");
  out << "P: " << std::fixed << std::setprecision(2) << purityEl << " #pm " << purity_errEl;
  if(m_mcflag) out << " (" << PURITY << ")";
  pt->AddText(out.str().c_str());
  pt->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  pt->Draw();

  TLine *de_line_RIGHT;
  de_line_RIGHT = new TLine(cuts->get_de_max_h0(m_mode,m_h0mode),0,cuts->get_de_max_h0(m_mode,m_h0mode),de_line_size());
  de_line_RIGHT->SetLineColor(kRed);
  de_line_RIGHT->SetLineStyle(1);
  de_line_RIGHT->SetLineWidth((Width_t)2.);
  de_line_RIGHT->Draw();
  TLine *de_line_LEFT;
  de_line_LEFT = new TLine(cuts->get_de_min_h0(m_mode,m_h0mode),0,cuts->get_de_min_h0(m_mode,m_h0mode),de_line_size());
  de_line_LEFT->SetLineColor(kRed);
  de_line_LEFT->SetLineStyle(1);
  de_line_LEFT->SetLineWidth((Width_t)2.);
  de_line_LEFT->Draw();

  pad4->cd(); pad4->SetLeftMargin(0.15); pad4->SetFillColor(0);
  dePull->SetMarkerSize(0.05); dePull->Draw();
  TLine *de_lineUP = new TLine(cuts->get_de_fit_min(),3,cuts->get_de_fit_max(),3);
  de_lineUP->SetLineColor(kBlue);
  de_lineUP->SetLineStyle(2);
  de_lineUP->Draw();
  TLine *de_line = new TLine(cuts->get_de_fit_min(),0,cuts->get_de_fit_max(),0);
  de_line->SetLineColor(kBlue);
  de_line->SetLineStyle(1);
  de_line->SetLineWidth((Width_t)2.);
  de_line->Draw();
  TLine *de_lineDOWN = new TLine(cuts->get_de_fit_min(),-3,cuts->get_de_fit_max(),-3);
  de_lineDOWN->SetLineColor(kBlue);
  de_lineDOWN->SetLineStyle(2);
  de_lineDOWN->Draw();

  cm->Update();

  out.str("");
  if(m_mcflag) out << "pics/de_purity";
  else         out << "pics/de_data";
  if(m_sig_bkg_flag == 1) out << "_zero_bkg";
  if(m_sig_bkg_flag == 2) out << "_zero_sig";
  out << "_m" << m_mode << "_h0m" << m_h0mode;
  if(m_bb_or_cnt_flag == 1)      out << "_cont";
  else if(m_bb_or_cnt_flag == 2) out << "_BB";
  out << "_ns" << m_nstr << "_cs" << m_cstr;
  const string name = out.str();
  const string epsname = name + string(".eps");
  cm->Print(epsname.c_str());
  out.str(""); out << "evince " << epsname << " &";
  system(out.str().c_str());
  const string rootname = name + string(".root");
  cm->Print(rootname.c_str());
}

void PurityFit::DrawDeltaE2(RooDataSet* ds){
//  RooProdPdf* pdf_cmb_bb = pdf_gen_cmb->GetPdfBB();
//  RooProdPdf* pdf_cmb_qq = pdf_gen_cmb->GetPdfQQ();
//  RooAbsPdf*  pdf_sig    = pdf_gen_sig->GetPdf();
//  RooProdPdf* pdf_part   = pdf_gen_part->GetPdf();

  stringstream out;
  TCanvas* cvec[2];
  cvec[0] = new TCanvas("c1","c1",1600,600);// Canvas for B0
  cvec[0]->Divide(8,2);
  cvec[1] = new TCanvas("c2","c2",1600,600);// Canvas for anti-B0
  cvec[1]->Divide(8,2);
  for(int k=0; k<2; k++){
    for(int j=0; j<16; j++){
      RooPlot* deFrame = de->frame(50);
      out.str(""); out << "flv == " << cuts->flv(k) << " && bin == " << cuts->bin(j);
      RooDataSet* ds0 = (RooDataSet*)ds->reduce(Cut(out.str().c_str()));
      ds0->plotOn(deFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("mbcSignal"),Cut(out.str().c_str()));
//      pdfI[k][j]->plotOn(deFrame,Components(*m_pdf_sig),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      pdfI[k][j]->plotOn(deFrame,Components(*m_pdf_part),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      pdfI[k][j]->plotOn(deFrame,Components(*m_pdf_comb_bb),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      pdfI[k][j]->plotOn(deFrame,Components(*m_pdf_comb_qq),LineStyle(kDashed),ProjectionRange("mbcSignal"));
//      pdfI[k][j]->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"));
      simpdf->plotOn(deFrame,LineWidth(2),ProjectionRange("mbcSignal"),ProjWData(*binflv,*ds0));
      cvec[k]->cd(j+1);
      out.str(""); out << "pad_" << k << "_" << j;
      TPad *pad = new TPad(out.str().c_str(),out.str().c_str(),0.01,0.20,0.99,0.99);
      pad->Draw();
      pad->cd();
      pad->SetLeftMargin(0.15);
      pad->SetFillColor(0);

      deFrame->GetXaxis()->SetTitleSize(0.05);
      deFrame->GetXaxis()->SetTitleOffset(0.85);
      deFrame->GetXaxis()->SetLabelSize(0.04);
      deFrame->GetYaxis()->SetTitleOffset(1.6);
      deFrame->Draw();
    }
    out.str("");
    if(m_mcflag) out << "pics/de_purity2";
    else       out << "pics/de_data2";
    out << "_m" << m_mode << "_h0m" << m_h0mode;
    out << "_" << k;
    const string name = out.str();
    const string epsname = name + string(".eps");
    cvec[k]->Print(epsname.c_str());
    out.str(""); out << "evince " << epsname << " &";
    system(out.str().c_str());
    const string rootname = name + string(".root");
    cvec[k]->Print(rootname.c_str());
  }
  return;
}

void PurityFit::DrawMbc(RooAddPdf* pdf,RooDataSet* ds){
  RooPlot* mbcFrame = mbc->frame();
  RooProdPdf* pdf_cmb_bb = pdf_gen_cmb->GetPdfBB();
  RooProdPdf* pdf_cmb_qq = pdf_gen_cmb->GetPdfQQ();
  RooAbsPdf*  pdf_sig    = pdf_gen_sig->GetPdf();
  RooProdPdf* pdf_part   = pdf_gen_part->GetPdf();
  ds->plotOn(mbcFrame,DataError(RooAbsData::SumW2),MarkerSize(1),CutRange("deSignal"));
  pdf->plotOn(mbcFrame,Components(*pdf_cmb_bb),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf->plotOn(mbcFrame,Components(*pdf_cmb_qq),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf->plotOn(mbcFrame,Components(*pdf_sig),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf->plotOn(mbcFrame,Components(*pdf_part),LineStyle(kDashed),ProjectionRange("deSignal"));
  pdf->plotOn(mbcFrame,LineWidth(2),ProjectionRange("deSignal"));

  RooHist* hmbcpull = mbcFrame->pullHist();
  RooPlot* mbcPull  = mbc->frame(Title("#Delta E pull distribution"));
  mbcPull->addPlotable(hmbcpull,"P");
  mbcPull->GetYaxis()->SetRangeUser(-5,5);

  TCanvas* cmmbc = new TCanvas("M_{bc}, Signal","M_{bc}, Signal",600,700);
  cmmbc->cd();

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.20);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);

  mbcFrame->GetXaxis()->SetTitleSize(0.05);
  mbcFrame->GetXaxis()->SetTitleOffset(0.85);
  mbcFrame->GetXaxis()->SetLabelSize(0.04);
  mbcFrame->GetYaxis()->SetTitleOffset(1.6);
  mbcFrame->Draw();

  TPaveText *ptmbc = new TPaveText(0.2,0.6,0.7,0.9,"brNDC");
  ptmbc->SetFillColor(0);
  ptmbc->SetTextAlign(12);
  stringstream out;
  out.str("");
  out << "#chi^{2}/n.d.f = " << mbcFrame->chiSquare();
  ptmbc->AddText(out.str().c_str());
  out.str("");
  out << "S: " << (int)(nsigEl+0.5) << " #pm " << (int)(nsig_errEl_total+0.5);
  if(m_mcflag) out << " (" << NSIGNAL_ELLI << ")";
  ptmbc->AddText(out.str().c_str());
  out.str("");
  out << "P: " << std::fixed << std::setprecision(2) << purityEl << " #pm " << purity_errEl;
  if(m_mcflag) out << " (" << PURITY << ")";
  ptmbc->AddText(out.str().c_str());
  ptmbc->AddText(cuts->GetLabel(m_mode,m_h0mode).c_str());
  ptmbc->Draw();

  TLine *mbc_line_RIGHT;
  mbc_line_RIGHT = new TLine(cuts->get_mbc_max_h0(m_mode,m_h0mode),0,cuts->get_mbc_max_h0(m_mode,m_h0mode),mbc_line_size());
  mbc_line_RIGHT->SetLineColor(kRed);
  mbc_line_RIGHT->SetLineStyle(1);
  mbc_line_RIGHT->SetLineWidth((Width_t)2.);
  mbc_line_RIGHT->Draw();
  TLine *mbc_line_LEFT;
  mbc_line_LEFT = new TLine(cuts->get_mbc_min_h0(m_mode,m_h0mode),0,cuts->get_mbc_min_h0(m_mode,m_h0mode),mbc_line_size());
  mbc_line_LEFT->SetLineColor(kRed);
  mbc_line_LEFT->SetLineStyle(1);
  mbc_line_LEFT->SetLineWidth((Width_t)2.);
  mbc_line_LEFT->Draw();

  pad2->cd();
  pad2->SetLeftMargin(0.15);
  pad2->SetFillColor(0);
  mbcPull->SetMarkerSize(0.05);
  mbcPull->Draw();
  TLine *mbc_lineUP = new TLine(cuts->get_mbc_fit_min(),3,cuts->get_mbc_fit_max(),3);
  mbc_lineUP->SetLineColor(kBlue);
  mbc_lineUP->SetLineStyle(2);
  mbc_lineUP->Draw();
  TLine *mbc_line = new TLine(cuts->get_mbc_fit_min(),0,cuts->get_mbc_fit_max(),0);
  mbc_line->SetLineColor(kBlue);
  mbc_line->SetLineStyle(1);
  mbc_line->SetLineWidth((Width_t)2.);
  mbc_line->Draw();
  TLine *mbc_lineDOWN = new TLine(cuts->get_mbc_fit_min(),-3,cuts->get_mbc_fit_max(),-3);
  mbc_lineDOWN->SetLineColor(kBlue);
  mbc_lineDOWN->SetLineStyle(2);
  mbc_lineDOWN->Draw();
  cmmbc->Update();

  out.str("");
  if(m_mcflag) out << "pics/mbc_purity";
  else         out << "pics/mbc_data";
  if(m_sig_bkg_flag == 1) out << "_zero_bkg";
  if(m_sig_bkg_flag == 2) out << "_zero_sig";
  out << "_m" << m_mode << "_h0m" << m_h0mode;
  if(m_bb_or_cnt_flag == 1)      out << "_cont";
  else if(m_bb_or_cnt_flag == 2) out << "_BB";
  out << "_ns" << m_nstr << "_cs" << m_cstr;
  const string name = out.str();
  const string epsname = name + string(".eps");
  cmmbc->Print(epsname.c_str());
  out.str(""); out << "evince " << epsname << " &";
  system(out.str().c_str());
  const string rootname = name + string(".root");
  cmmbc->Print(rootname.c_str());
}

int PurityFit::CalcWTagAndNEv(RooDataSet* ds,double* wtag_arr, int (&nev_arr)[2][16], const bool mcbin, const bool mcflv){
  double m_tag_LH,m_t_asc,m_t_rec;
  int m_exp,m_bin,flv;
  const int NTot = ds->sumEntries();
  for(int i=0; i<16; i++){
    wtag_arr[i] = 0;
    nev_arr[0][i] = 0;
    nev_arr[1][i] = 0;
  }
  // ** Average wrong tagging probability calculation ** //
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    m_t_asc = aset->getRealValue("z_asc");
    m_t_rec = aset->getRealValue("z_sig");
    if(abs(m_t_rec-m_t_asc)>70./7.848) continue;

    m_bin = mcbin ? aset->getCatIndex("bin_mc") : aset->getCatIndex("bin");
    if(!mcflv){
      m_tag_LH = aset->getRealValue("tag_LH");
      m_exp    = aset->getCatIndex("exp");
      flv = cuts->q(m_tag_LH);
    } else{
      flv = aset->getCatIndex("flv_mc");
    }

    const int bin_ind = cuts->bin_ind(m_bin);
    if(!mcflv)wtag_arr[bin_ind] += cuts->get_wtag_prob(m_tag_LH,m_exp,!m_mcflag);

    const int flv_ind = cuts->flv_ind(flv);
    nev_arr[flv_ind][bin_ind]++;
  }
  cout << "Data set in bins:" << endl;
  for(int i=0; i<16; i++){
    cout << "bin " << cuts->bin(i);
    cout << " N: " << nev_arr[0][i];
    cout << " Nbar: " << nev_arr[1][i];
    cout << endl;
    if(!mcflv) wtag_arr[i] /= (nev_arr[0][i] + nev_arr[1][i]);
  }
  if(!mcflv){
    cout << "Average wrong tag probability: " << endl;
    for(int i=0; i<8; i++){
      cout <<   i+1  << " bin :" << wtag_arr[cuts->bin_ind(i+1)] << ", ";
      cout << -(i+1) << " bin :" << wtag_arr[cuts->bin_ind(-(i+1))] << endl;
    }
  }
  // * wrong tagging probability is calculated * //
  return 0;
}

void PurityFit::WriteBinFlvMap(RooDataSet* ds,const char* fname){
  cout << "Saving Bin Map to file " << fname << endl;
  ofstream ofile;
  ofile.open(fname);
  Roo1DTable* map = ds->table(RooArgSet(*flv,*bin));
  map->Print();
  stringstream out;
  ofile << "Bin Flavor Map. Mode " << m_mode << ", h0 mode " << m_h0mode << endl;
  for(int k=0; k<2; k++){
    ofile << "Flv " << cuts->flv(k) << endl;
    for(int i=0; i<16; i++){
      out.str("");
      out << (cuts->flv(k)>0 ? "tag_LH<0" : "tag_LH>0");
      out << " && bin == " << cuts->bin(i);
      ofile << " bin " << cuts->bin(i) << ", value " << ds->sumEntries(out.str().c_str()) << endl;
    }
  }
  ofile.close();
  return;
}

int PurityFit::MakePredictions(const double* wtag_arr,const int (&nev_arr)[2][16],double (&sig_arr)[2][16],double (&sig_err_arr)[2][16],double (&cmb_arr)[2][16],double (&cnt_arr)[2][16],double (&prt_arr)[2][16],const char* fname){
  const bool write_flag = string(fname) == string("") ? false : true;
  cout << "Making predictions..." << endl;
  const double f_bb        = fbb->getVal();
  const double denominator = 1.+f_p_f_bbc->getVal()*f_bb;
  const double part_coeff  = (f_p_f_bbc->getVal()*f_bb)/denominator;
  const double comb_coeff  = 1./denominator;

  int m_bin,m_flv;
  double bkg_arr[2][16];
  for(int i=0; i<16; i++){
    sig_arr[0][i] = 0; sig_arr[1][i] = 0;
    cmb_arr[0][i] = 0; cmb_arr[1][i] = 0;
    cnt_arr[0][i] = 0; cnt_arr[1][i] = 0;
    prt_arr[0][i] = 0; prt_arr[1][i] = 0;
    bkg_arr[0][i] = 0; bkg_arr[1][i] = 0;
  }
  for(int i=0; i<16; i++){
    for(int j=0; j<2; j++){
      m_bin = cuts->bin(i);
      m_flv = cuts->flv(j);
      sig_arr[j][i] = 0.5*cuts->N(m_bin,m_flv,wtag_arr[i])*Nsig->getVal();
      const double F = cuts->N(m_bin,m_flv,wtag_arr[i]);
      const double err1   = 0.5*F*Nsig->getError();
      const double err2sq = 0.5*F*(1.-F)*Nsig->getVal();
      sig_err_arr[j][i] = sqrt(err1*err1+err2sq);
      bkg_arr[j][i] = nev_arr[j][i] - sig_arr[j][i];
      if(bkg_arr[j][i]<0) bkg_arr[j][i] = 0;
      prt_arr[j][i] = bkg_arr[j][i]*part_coeff;
      cmb_arr[j][i] = bkg_arr[j][i]*comb_coeff;
      cnt_arr[j][i] = cmb_arr[j][i]*(1.-f_bb);
      cout << "bin " << m_bin << " flv " << m_flv;
      cout << " Tot: " << nev_arr[j][i];
      cout << " Sig: " << sig_arr[j][i];
      cout << " Cmb: " << cmb_arr[j][i];
      cout << " Cnt: " << cnt_arr[j][i];
      cout << " Prt: " << prt_arr[j][i];
      cout << " Pur: " << 100.*sig_arr[j][i]/(sig_arr[j][i]+cmb_arr[j][i]+prt_arr[j][i]);
      cout << endl;
    }
  }
  if(write_flag){
    ofstream ofile;
    ofile.open(fname);
    ofile << "Predictions for mode " << m_mode << ", h0 mode " << m_h0mode << endl;
    ofile << "NSig = " << Nsig->getVal() << " +- " << Nsig->getError() << endl;
    ofile << "fbb  = " << fbb->getVal()  << " +- " << fbb->getError()  << endl;
    ofile << "fprt = " << f_p_f_bbc->getVal() << " +- " << f_p_f_bbc->getError() << endl;
    for(int i=0; i<16; i++){
      for(int j=0; j<2; j++){
        m_bin = cuts->bin(i);
        m_flv = cuts->flv(j);
        ofile << "bin " << m_bin << " flv " << m_flv;
        ofile << " Tot: " << nev_arr[j][i];
        ofile << " Sig: " << sig_arr[j][i];
        ofile << " Cmb: " << cmb_arr[j][i];
        ofile << " Cnt: " << cnt_arr[j][i];
        ofile << " Prt: " << prt_arr[j][i];
        ofile << " Pur: " << 100.*sig_arr[j][i]/(sig_arr[j][i]+cmb_arr[j][i]+prt_arr[j][i]);
        ofile << endl;
      }
    }
    ofile.close();
  }

  return 0;
}

int PurityFit::MakePredictions2(const int (&nev_arr)[2][16],double (&sig_arr)[2][16],double (&sig_err_arr)[2][16],double (&cmb_arr)[2][16],double (&cnt_arr)[2][16],double (&prt_arr)[2][16],const char* fname){
  const bool write_flag = string(fname) == string("") ? false : true;
  cout << "Making predictions 2..." << endl;
  double cnt_arr_err[2][16];
  double prt_arr_err[2][16];
  for(int i=0; i<16; i++){
    sig_arr[0][i] = 0; sig_arr[1][i] = 0;
    cmb_arr[0][i] = 0; cmb_arr[1][i] = 0;
    cnt_arr[0][i] = 0; cnt_arr[1][i] = 0;
    prt_arr[0][i] = 0; prt_arr[1][i] = 0;
  }
  for(int i=0; i<16; i++){
    for(int j=0; j<2; j++){
      sig_arr[j][i] = NsigI[j][i]->getVal();
      const double F = sig_arr[j][i]/Nsig->getVal();
      const double err1   = 0.5*F*Nsig->getError();
      const double err2sq = 0.5*(1.-F)*NsigI[j][i]->getVal();
      sig_err_arr[j][i] = sqrt(err1*err1+err2sq);

      cmb_arr[j][i] = NcmbI[j][i]->getVal();

      cnt_arr[j][i] = NcmbI[j][i]->getVal()*(1.-fbbI[j][i]->getVal());
      const double cnt_err1 = (1.-fbbI[j][i]->getVal())*NcmbI[j][i]->getError();
      const double cnt_err2 = NcmbI[j][i]->getVal()*fbbI[j][i]->getError();
      cnt_arr_err[j][i] = sqrt(cnt_err1*cnt_err1+cnt_err2*cnt_err2);

      prt_arr[j][i] = NcmbI[j][i]->getVal()*fbbI[j][i]->getVal()*f_p_f_bbc->getVal();
      const double prt_err1 = fbbI[j][i]->getVal()*f_p_f_bbc->getVal()*NcmbI[j][i]->getError();
      const double prt_err2 = fbbI[j][i]->getVal()*f_p_f_bbc->getError()*NcmbI[j][i]->getVal();
      const double prt_err3 = fbbI[j][i]->getError()*f_p_f_bbc->getVal()*NcmbI[j][i]->getVal();
      prt_arr_err[j][i] = sqrt(prt_err1*prt_err1+prt_err2*prt_err2+prt_err3*prt_err3);
      cout << "bin " << cuts->bin(i) << " flv " << cuts->flv(j);
      cout << " Tot: " << nev_arr[j][i];
      cout << " Sig: " << sig_arr[j][i];
      cout << " Cmb: " << cmb_arr[j][i];
      cout << " Cnt: " << cnt_arr[j][i];
      cout << " Prt: " << prt_arr[j][i];
      cout << " Pur: " << 100.*sig_arr[j][i]/(sig_arr[j][i]+cmb_arr[j][i]+prt_arr[j][i]);
      cout << endl;
    }
  }
  if(write_flag){
    ofstream ofile;
    ofile.open(fname);
    ofile << "Predictions for mode " << m_mode << ", h0 mode " << m_h0mode << endl;
    ofile << "NSig = " << Nsig->getVal() << " +- " << Nsig->getError() << endl;
    ofile << "fbb  = " << fbb->getVal()  << " +- " << fbb->getError()  << endl;
    ofile << "fprt = " << f_p_f_bbc->getVal() << " +- " << f_p_f_bbc->getError() << endl;
    for(int i=0; i<16; i++){
      for(int j=0; j<2; j++){
        ofile << "bin " << cuts->bin(i) << " flv " << cuts->flv(j);
        ofile << " Tot: " << nev_arr[j][i];
        ofile << " Sig: " << sig_arr[j][i] << " +- " << sig_err_arr[j][i];
        ofile << " Cmb: " << cmb_arr[j][i] << " +- " << NcmbI[j][i]->getError();
        ofile << " Cnt: " << cnt_arr[j][i] << " +- " << cnt_arr_err[j][i];
        ofile << " Prt: " << prt_arr[j][i] << " +- " << prt_arr_err[j][i];
        ofile << " Pur: " << 100.*sig_arr[j][i]/(sig_arr[j][i]+cmb_arr[j][i]+prt_arr[j][i]);
        ofile << endl;
      }
    }
    ofile.close();
  }

  return 0;
}

void PurityFit::CalcBinsFractions(RooDataSet* ds){
  const int NTot = ds->sumEntries();
  cout << "Calculating wrong tag probability in bins with " << NTot << " events...  " << endl;

  stringstream out;
  if(m_mcflag){
    out.str("");
    out << "params/Predictions_m" << m_mode << "_h0m" << m_h0mode;
    if(m_fixshape) out << "_fixshape";
    out << "_bin_mc.txt";
    CalcWTagAndNEv(ds,WrTagMap_bin_mc,EventsMap_bin_mc,true,false);
    MakePredictions(WrTagMap_bin_mc,EventsMap_bin_mc,NSigPredicted_bin_mc,NSigPredicted_err_bin_mc,NCmbPredicted_bin_mc,NCntPredicted_bin_mc,NPrtPredicted_bin_mc,out.str().c_str());

    out.str("");
    out << "params/Predictions_m" << m_mode << "_h0m" << m_h0mode;
    if(m_fixshape) out << "_fixshape";
    out << "_flv_mc.txt";
    CalcWTagAndNEv(ds,WrTagMap_flv_mc,EventsMap_flv_mc,false,true);
    MakePredictions(WrTagMap_flv_mc,EventsMap_flv_mc,NSigPredicted_flv_mc,NSigPredicted_err_flv_mc,NCmbPredicted_flv_mc,NCntPredicted_flv_mc,NPrtPredicted_flv_mc,out.str().c_str());

    out.str("");
    out << "params/Predictions_m" << m_mode << "_h0m" << m_h0mode;
    if(m_fixshape) out << "_fixshape";
    out << "_bin_flv_mc.txt";
    CalcWTagAndNEv(ds,WrTagMap_mc,EventsMap_mc,true,true);
    MakePredictions(WrTagMap_mc,EventsMap_mc,NSigPredicted_mc,NSigPredicted_err_mc,NCmbPredicted_mc,NCntPredicted_mc,NPrtPredicted_mc,out.str().c_str());
  }

  CalcWTagAndNEv(ds,WrTagMap,EventsMap,false,false);
  out.str("");
  out << "params/Predictions_m" << m_mode << "_h0m" << m_h0mode;
  if(m_fixshape) out << "_fixshape";
  if(m_mcflag) out << "_mc";
  out << ".txt";
  MakePredictions(WrTagMap,EventsMap,NSigPredicted,NSigPredicted_err,NCmbPredicted,NCntPredicted,NPrtPredicted,out.str().c_str());

  return;
}

void PurityFit::SaveCPVTree(RooDataSet* ds){
  const int NTot = ds->sumEntries();
  const double mbc0 = mbc_center->getVal();
  const double de0  = de_center->getVal();
  const double mbcr = mbc_radius->getVal();
  const double der  = de_radius->getVal();

  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  ICPVEv ev;
  stringstream out;
  out.str("");
  out << "data/";
  if(!m_mcflag) out << "data_cpv_tree";
  else        out << "genmc_cpv_tree";
  if(m_sig_bkg_flag == 1) out << "_zero_sig";
  out << "_m" << m_mode << "_hm" << m_h0mode;
  if(m_bb_or_cnt_flag == 1)      out << "_cont";
  else if(m_bb_or_cnt_flag == 2) out << "_BB";
  if(true) out << "_ww";
  out << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  double sig_val, cmb_val, part_val, f_cmb;

  int zerobin = 0;
  const int mctype = m_mcflag ? 2 : 0;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,mctype);
    ev.sigarea = cuts->EllipsR2(ev.de,ev.mbc,mbc0,de0,mbcr,der)>1 ? 0 : 1;
    if(!true && !ev.sigarea) continue;
    if(!ev.bin){
      zerobin++;
      continue;
    }
    de->setVal(ev.de); mbc->setVal(ev.mbc);
    const int bin_ind = cuts->bin_ind(ev.bin);
    const int flv_ind = cuts->flv_ind(cuts->q(ev.tag_LH));
    ev.fbb      = pdf_gen_cmb->Get_fbb()->getVal();
    ev.fbb_err  = pdf_gen_cmb->Get_fbb()->getError();
    ev.Nsig     = Nsig->getVal();//NSigPredicted[flv_ind][bin_ind];
    ev.Nsig_err = Nsig->getError();//NSigPredicted_err[flv_ind][bin_ind];
    ev.psig     = pdf_sig->getVal(RooArgSet(*de,*mbc));
    ev.pcmb     = pdf_comb->getVal(RooArgSet(*de,*mbc));
    ev.pprt     = pdf_part->getVal(RooArgSet(*de,*mbc));
    ev.pcnt     = pdf_cont->getVal(RooArgSet(*de,*mbc));
    ev.fprt     = f_p_f_bbc->getVal();
    ev.fprt_err = f_p_f_bbc->getError();

    ev.wrtag    = WrTagMap[bin_ind];
    ev.NTot     = EventsMap[flv_ind][bin_ind];
    ev.wrtag_bin_mc = WrTagMap_bin_mc[bin_ind];
    ev.NTot_bin_mc  = EventsMap_bin_mc[flv_ind][bin_ind];
    ev.wrtag_flv_mc = WrTagMap_flv_mc[bin_ind];
    ev.NTot_flv_mc  = EventsMap_flv_mc[flv_ind][bin_ind];
    ev.wrtag_mc = WrTagMap_mc[bin_ind];
    ev.NTot_mc  = EventsMap_mc[flv_ind][bin_ind];

    sig_val  = NSigPredicted[flv_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted[flv_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted[flv_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont         = ev.f_cont_in_comb*f_cmb;

    if(m_mcflag){
      const int bin_mc_ind = cuts->bin_ind(ev.bin_mc);
      const int flv_mc_ind = cuts->flv_ind(ev.flv_mc);
      sig_val  = NSigPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
      cmb_val  = NCmbPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
      part_val = NPrtPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

      ev.f_bkg_bin_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
      f_cmb = cmb_val/(cmb_val + part_val + sig_val);
      ev.f_cont_in_comb_bin_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_bin_mc[flv_ind][bin_ind]/(cmb_val+part_val);
      ev.f_cont_bin_mc         = ev.f_cont_in_comb_bin_mc*f_cmb;

      sig_val  = NSigPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
      cmb_val  = NCmbPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
      part_val = NPrtPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

      ev.f_bkg_flv_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
      f_cmb = cmb_val/(cmb_val + part_val + sig_val);
      ev.f_cont_in_comb_flv_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_flv_mc[flv_ind][bin_ind]/(cmb_val+part_val);
      ev.f_cont_flv_mc         = ev.f_cont_in_comb_flv_mc*f_cmb;

      sig_val  = NSigPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
      cmb_val  = NCmbPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
      part_val = NPrtPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

      ev.f_bkg_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
      f_cmb = cmb_val/(cmb_val + part_val + sig_val);
      ev.f_cont_in_comb_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_mc[flv_ind][bin_ind]/(cmb_val+part_val);
      ev.f_cont_mc         = ev.f_cont_in_comb_mc*f_cmb;
    }

    tree->Fill();
  }
  cout << "Writing a tree... ( zerobin: " << zerobin << ")" << endl;

  tree->Write();
  file->Close();
}

void PurityFit::SaveSigLineCPVTrees(void){
  vector<string> angle_vec;
  angle_vec.push_back(string("0"));
  angle_vec.push_back(string("15"));
  angle_vec.push_back(string("22_5"));
  angle_vec.push_back(string("30"));
  angle_vec.push_back(string("45"));
  angle_vec.push_back(string("60"));
  angle_vec.push_back(string("67_5"));
  angle_vec.push_back(string("75"));
  angle_vec.push_back(string("90"));
  angle_vec.push_back(string("105"));
  angle_vec.push_back(string("112_5"));
  angle_vec.push_back(string("120"));
  angle_vec.push_back(string("135"));
  angle_vec.push_back(string("150"));
  angle_vec.push_back(string("157_5"));
  angle_vec.push_back(string("165"));
  for(int i=0; i<angle_vec.size(); i++){
    SaveSigCPVTree(true,angle_vec[i]);
  }
  return;
}

void PurityFit::SaveSigCPVTree(const bool line_flag,const string& angle){
  RooDataSet* ds = line_flag ? GetSigMCLineDataSet(angle) : GetSigMCDataSet();
  const int NTot = ds->sumEntries();
  DefineElliRange();
  const double mbc0 = mbc_center->getVal();
  const double de0  = de_center->getVal();
  const double mbcr = mbc_radius->getVal();
  const double der  = de_radius->getVal();

  ICPVEv ev;
  stringstream out;
  out.str("");
  out << "data/";
  out << "sigmc_cpv_tree";
  out << "_m" << m_mode << "_hm" << m_h0mode;
  if(line_flag) out << "_line" << angle;
  out << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  int zerobin = 0;
  ev.f_bkg = 0;
  ev.f_cont_in_comb = 0;
  ev.f_cont= 0;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,1);
    ev.sigarea = cuts->EllipsR2(ev.de,ev.mbc,mbc0,de0,mbcr,der)>1 ? 0 : 1;
    if(!ev.bin){
      zerobin++;
      continue;
    }
    tree->Fill();
  }
  cout << "Writing a tree... ( zerobin: " << zerobin << ")" << endl;

  tree->Write();
  file->Close();
}

//void PurityFit::GenerateCPVTree(void){
//  cout << "*** Preparing simulated CPV tree (SigMC + Simulated background) ***" << endl;
//  const int  NSigSvd1 = sig_ds_svd1->sumEntries();
//  const int  NSigSvd2 = sig_ds_svd2->sumEntries();

//  const double mbc0 = mbc_center->getVal();
//  const double de0  = de_center->getVal();
//  const double mbcr = mbc_radius->getVal();
//  const double der  = de_radius->getVal();

//  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
//  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
//  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();

//  ICPVEv ev;
//  TFile* file = new TFile("mixtree.root","RECREATE");
//  cout << "Initiating tree..." << endl;
//  TTree* tree = GetCPVTree(ev);
//}

void PurityFit::FillFlvBinsArrays(vector<ICPVEv> (&ev_arr)[2][16],RooDataSet* ds){
  const int  NTot = ds->sumEntries();
  ICPVEv ev;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,1);
    if(abs(ev.flv)>1 || abs(ev.bin)>8) continue;
    ev_arr[cuts->flv_ind(ev.flv)][cuts->bin_ind(ev.bin)].push_back(ev);
  }
  return;
}

void PurityFit::FillEvtVec(vector<ICPVEv> &ev_vec,RooDataSet* ds){
  const int  NTot = ds->sumEntries();
  ICPVEv ev;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,1);
    ev_vec.push_back(ev);
  }
  return;
}

double PurityFit::GetFBkg(int& evnum, double& m_de, double& m_mbc, const int flv_ind, const int bin_ind){
  const double mbc0 = mbc_center->getVal();
  const double de0  = de_center->getVal();
  const double mbcr = mbc_radius->getVal();
  const double der  = de_radius->getVal();

  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  do{
    if(evnum>=de_mbc_ds->sumEntries()) return -1;
    const RooArgSet* aset  = de_mbc_ds->get(evnum);
    m_de  = aset->getRealValue("de");
    m_mbc = aset->getRealValue("mbc");
    evnum++;
  } while(cuts->EllipsR2(m_de,m_mbc,mbc0,de0,mbcr,der)>1);
  de->setVal(m_de); mbc->setVal(m_mbc);
  const double sig_val  = NSigPredicted[flv_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));//sigint;
  const double cmb_val  = NCmbPredicted[flv_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));//cmbint;
  const double part_val = NPrtPredicted[flv_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));//partint;
  const double f_bkg = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
  return f_bkg;
}

bool PurityFit::IsGoodNumbers(const ICPVEv& ev){
  if(abs(ev.bin) > 8){
    cout << "Bad Numbers: bin " << ev.bin << endl;
    return false;
  }
  return true;
}

void PurityFit::MixCPVTree(RooDataSet* ds){
  cout << "*** Preparing mixed CPV tree (SigMC + GenMC) ***" << endl;
  vector<ICPVEv> ev_sig_svd1_arr[2][16];
  vector<ICPVEv> ev_sig_svd2_arr[2][16];
//  const int SampleNum = 0;
  int sig_index_svd1[2][16];
  int sig_index_svd2[2][16];
//  int deltanum;
//  if(m_mode == 1)      deltanum = 50;
//  else if(m_mode == 2) deltanum = 20;
//  else if(m_mode == 3) deltanum = 50;
//  else if(m_mode == 5) deltanum = 20;
  for(int k=0; k<2; k++){
    for(int i=0; i<16; i++){
      sig_index_svd1[k][i] = m_cstr;//SampleNum*deltanum;
      sig_index_svd2[k][i] = m_cstr;//SampleNum*deltanum;
    }
  }
//  vector<ICPVEv> ev_bkg_svd1_arr[2][16];
//  vector<ICPVEv> ev_bkg_svd2_arr[2][16];

  FillFlvBinsArrays(ev_sig_svd1_arr,sig_ds_svd1);
  FillFlvBinsArrays(ev_sig_svd2_arr,sig_ds_svd2);
//  FillFlvBinsArrays(ev_bkg_svd1_arr,raw_back_ds_svd1);
//  FillFlvBinsArrays(ev_bkg_svd2_arr,raw_back_ds_svd2);

  ICPVEv ev, sev;
  stringstream out;
  out.str("");
  out << "data/";
  out << "mixcpvtree_m" << m_mode << "_mh0" << m_h0mode;
  out << "_ns" << m_nstr << "_cs" << m_cstr;
  out << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  cout << "Tree is ready: " << tree->GetEntries() << endl;

  const int NTot = ds->sumEntries();
  const double mbc0 = mbc_center->getVal();
  const double de0  = de_center->getVal();
  const double mbcr = mbc_radius->getVal();
  const double der  = de_radius->getVal();

  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();
  const double f_cont_in_comb = 1. - pdf_gen_cmb->Get_fbb()->getVal();

  int zerobin = 0;
  double sig_val,cmb_val,part_val,f_cmb;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    FillEvent(aset,ev,2);
    ev.sigarea = cuts->EllipsR2(ev.de,ev.mbc,mbc0,de0,mbcr,der)>1 ? 0 : 1;
    if(!ev.bin){
      zerobin++;
      continue;
    }
    de->setVal(ev.de); mbc->setVal(ev.mbc);
    const int bin_ind = cuts->bin_ind(ev.bin);
    const int flv_ind = cuts->flv_ind(cuts->q(ev.tag_LH));

    ev.fbb      = pdf_gen_cmb->Get_fbb()->getVal();
    ev.fbb_err  = pdf_gen_cmb->Get_fbb()->getError();
    ev.Nsig     = Nsig->getVal();//NSigPredicted[flv_ind][bin_ind];
    ev.Nsig_err = Nsig->getError();//NSigPredicted_err[flv_ind][bin_ind];
    ev.psig     = pdf_sig->getVal(RooArgSet(*de,*mbc));
    ev.pcmb     = pdf_comb->getVal(RooArgSet(*de,*mbc));
    ev.pprt     = pdf_part->getVal(RooArgSet(*de,*mbc));
    ev.pcnt     = pdf_cont->getVal(RooArgSet(*de,*mbc));
    ev.fprt     = f_p_f_bbc->getVal();
    ev.fprt_err = f_p_f_bbc->getError();
    ev.wrtag    = WrTagMap[bin_ind];
//    cout << ev.wrtag << " ";
    ev.NTot  = EventsMap[flv_ind][bin_ind];

    sig_val  = NSigPredicted[flv_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted[flv_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted[flv_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont         = ev.f_cont_in_comb*f_cmb;

    const int bin_mc_ind = cuts->bin_ind(ev.bin_mc);
    const int flv_mc_ind = cuts->flv_ind(ev.flv_mc);
    sig_val  = NSigPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg_bin_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb_bin_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_bin_mc[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont_bin_mc         = ev.f_cont_in_comb_bin_mc*f_cmb;

    sig_val  = NSigPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg_flv_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb_flv_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_flv_mc[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont_flv_mc         = ev.f_cont_in_comb_flv_mc*f_cmb;

    sig_val  = NSigPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg_mc = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_mc[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont_mc         = ev.f_cont_in_comb_mc*f_cmb;

    if(ev.b0f == 1 || ev.b0f == 5 || ev.b0f == 10){
      if(ev.exp > 30){
        CopyEvent(ev_sig_svd2_arr[flv_mc_ind][bin_mc_ind][sig_index_svd1[flv_mc_ind][bin_mc_ind]],sev);
        sig_index_svd1[flv_mc_ind][bin_mc_ind] += m_nstr;
      } else{
        CopyEvent(ev_sig_svd1_arr[flv_mc_ind][bin_mc_ind][sig_index_svd2[flv_mc_ind][bin_mc_ind]],sev);
        sig_index_svd2[flv_mc_ind][bin_mc_ind] += m_nstr;
      }
      CopyTimeEvent(sev,ev);
    }
    tree->Fill();
  }

//  for(int flv_ind=0; flv_ind<2; flv_ind++){
//    for(int bin_ind=0; bin_ind<16; bin_ind++){
//      int sig_index = 0;
//      int bkg_index = 0;
//      do{
//        if(!(Ncur%10000)) cout << "Event " << Ncur << endl;
//        const double f_bkg = GetFBkg(Ncur,m_de,m_mbc,flv_ind,bin_ind);
//        const bool sig_event = f_bkg < rndm.Rndm() ? true : false;
//        if(sig_event){
//          sig_index++;
//          CopyEvent(ev_sig_svd1_arr[flv_ind][bin_ind][sig_index],ev);
//          ev.de = m_de; ev.mbc = m_mbc; ev.f_bkg = f_bkg;
//          if(IsGoodNumbers(ev)) tree->Fill();
//        } else{
//          bkg_index++;
//          CopyEvent(ev_bkg_svd1_arr[flv_ind][bin_ind][bkg_index],ev);
//          ev.de = m_de; ev.mbc = m_mbc; ev.f_bkg = f_bkg;
//          if(IsGoodNumbers(ev)) tree->Fill();
//        }
//      } while(sig_index < ev_sig_svd1_arr[flv_ind][bin_ind].size() && bkg_index < ev_bkg_svd1_arr[flv_ind][bin_ind].size() && Ncur<NTot);
//      cout << cuts->flv(flv_ind) << " " << cuts->bin(bin_ind) << " " << sig_index << " (" << ev_sig_svd1_arr[flv_ind][bin_ind].size() << "), " << bkg_index << " (" << ev_bkg_svd1_arr[flv_ind][bin_ind].size() << "), " << Ncur << " (" << NTot << ")" << endl;
//    }
//  }

//  for(int flv_ind=0; flv_ind<2; flv_ind++){
//    for(int bin_ind=0; bin_ind<16; bin_ind++){
//      int sig_index = 0;
//      int bkg_index = 0;
//      do{
//        if(!(Ncur%10000)) cout << "Event " << Ncur << endl;
//        const double f_bkg = GetFBkg(Ncur,m_de,m_mbc,flv_ind,bin_ind);
//        const bool sig_event = f_bkg < rndm.Rndm() ? true : false;
//        if(sig_event){
//          sig_index++;
//          CopyEvent(ev_sig_svd2_arr[flv_ind][bin_ind][sig_index],ev);
//          ev.de = m_de; ev.mbc = m_mbc; ev.f_bkg = f_bkg;
//          if(IsGoodNumbers(ev)) tree->Fill();
//        } else{
//          bkg_index++;
//          CopyEvent(ev_bkg_svd2_arr[flv_ind][bin_ind][bkg_index],ev);
//          ev.de = m_de; ev.mbc = m_mbc; ev.f_bkg = f_bkg;
//          if(IsGoodNumbers(ev)) tree->Fill();
//        }
//      } while(sig_index < ev_sig_svd2_arr[flv_ind][bin_ind].size() && bkg_index < ev_bkg_svd2_arr[flv_ind][bin_ind].size() && Ncur<NTot);
//      cout << cuts->flv(flv_ind) << " " << cuts->bin(bin_ind) << " " << sig_index << " (" << ev_sig_svd2_arr[flv_ind][bin_ind].size() << "), " << bkg_index << " (" << ev_bkg_svd2_arr[flv_ind][bin_ind].size() << "), " << Ncur << " (" << NTot << ")" << endl;
//    }
//  }

//  cout << Ncur << " generated events are used" << endl;
  cout << "Writing a tree..." << endl;

  tree->Write();
  file->Close();
  return;
}

void PurityFit::MixCPVTree2(RooDataSet* ds){
  cout << "*** Preparing mixed CPV tree (SigMC + GenMC) ***" << endl;
  vector<ICPVEv> ev_sig_svd1_arr;
  vector<ICPVEv> ev_sig_svd2_arr;
  vector<int> index_svd1;
  vector<int> index_svd2;

  FillEvtVec(ev_sig_svd1_arr,sig_ds_svd1);
  FillEvtVec(ev_sig_svd2_arr,sig_ds_svd2);
  GetShuffledVector(ev_sig_svd1_arr.size(),index_svd1);
  GetShuffledVector(ev_sig_svd2_arr.size(),index_svd2);

  ICPVEv ev;
  stringstream out;
  out.str("");
  out << "data/";
  out << "mixcpvtree2_m" << m_mode << "_mh0" << m_h0mode;
  out << "_ns" << m_nstr << "_cs" << m_cstr;
  if(m_fixshape) out << "_fixshape";
  out << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  cout << "Initiating tree..." << endl;
  TTree* tree = GetCPVTree(ev);
  cout << "Tree is ready: " << tree->GetEntries() << endl;

  const int NTot    = ds->sumEntries();
  const double mbc0 = mbc_center->getVal();
  const double de0  = de_center->getVal();
  const double mbcr = mbc_radius->getVal();
  const double der  = de_radius->getVal();

  RooAbsPdf*  pdf_sig  = pdf_gen_sig->GetPdf();
  RooAddPdf*  pdf_comb = pdf_gen_cmb->GetPdf();
  RooProdPdf* pdf_part = pdf_gen_part->GetPdf();
  RooProdPdf* pdf_cont = pdf_gen_cmb->GetPdfQQ();

  int zerobin = 0;
  int cur_ind_svd1 = m_cstr;
  int cur_ind_svd2 = m_cstr;
  double sig_val,cmb_val,part_val,f_cmb;
  for(int i=0; i<NTot; i++){
    const RooArgSet* aset = ds->get(i);
    const int m_b0f = aset->getCatIndex("b0f");
    if(m_b0f == 1 || m_b0f == 5 || m_b0f == 10){
      // Signal event //
      if(aset->getCatIndex("exp")>30){
        CopyEvent(ev_sig_svd2_arr[index_svd2[cur_ind_svd2]],ev);
        cur_ind_svd2 += m_nstr;
        if(cur_ind_svd2>index_svd2.size()){
          cout << "MixCPVTree2: svd2 signal vector size is "  << index_svd2.size() << " while " << cur_ind_svd2 << " is requested." << endl;
          break;
        }
      } else{
        CopyEvent(ev_sig_svd1_arr[index_svd1[cur_ind_svd1]],ev);
        cur_ind_svd1 += m_nstr;
        if(cur_ind_svd1>index_svd1.size()){
          cout << "MixCPVTree2: svd1 signal vector size is "  << index_svd1.size() << " while " << cur_ind_svd1 << " is requested." << endl;
          break;
        }
      }
    } else{
      // Background event //
      FillEvent(aset,ev,2);
    }
    ev.sigarea = cuts->EllipsR2(ev.de,ev.mbc,mbc0,de0,mbcr,der)>1 ? 0 : 1;
    if(!ev.bin){
      zerobin++;
      continue;
    }
    de->setVal(ev.de); mbc->setVal(ev.mbc);
    const int bin_ind = cuts->bin_ind(ev.bin);
    const int flv_ind = cuts->flv_ind(ev.flv);
    const double f_cont_in_comb = 1. - fbbI[flv_ind][bin_ind]->getVal();

    ev.fbb      = fbbI[flv_ind][bin_ind]->getVal();
    ev.fbb_err  = fbbI[flv_ind][bin_ind]->getError();
    ev.Nsig     = Nsig->getVal();//NSigPredicted[flv_ind][bin_ind];
    ev.Nsig_err = Nsig->getError();//NSigPredicted_err[flv_ind][bin_ind];
    ev.psig     = pdf_sig->getVal(RooArgSet(*de,*mbc));
    ev.pcmb     = pdf_comb->getVal(RooArgSet(*de,*mbc));
    ev.pprt     = pdf_part->getVal(RooArgSet(*de,*mbc));
    ev.pcnt     = pdf_cont->getVal(RooArgSet(*de,*mbc));
    ev.fprt     = f_p_f_bbc->getVal();
    ev.fprt_err = f_p_f_bbc->getError();
    ev.wrtag    = WrTagMap[bin_ind];
    ev.NTot     = EventsMap[flv_ind][bin_ind];

    sig_val  = NSigPredicted[flv_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
    cmb_val  = NCmbPredicted[flv_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
    part_val = NPrtPredicted[flv_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

    ev.f_bkg          = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
    f_cmb             = cmb_val/(cmb_val + part_val + sig_val);
    ev.f_cont_in_comb = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted[flv_ind][bin_ind]/(cmb_val+part_val);
    ev.f_cont         = ev.f_cont_in_comb*f_cmb;

//    const int bin_mc_ind = cuts->bin_ind(ev.bin_mc);
//    const int flv_mc_ind = cuts->flv_ind(ev.flv_mc);
//    sig_val  = NSigPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
//    cmb_val  = NCmbPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
//    part_val = NPrtPredicted_bin_mc[flv_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

//    ev.f_bkg_bin_mc          = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
//    f_cmb                    = cmb_val/(cmb_val + part_val + sig_val);
//    ev.f_cont_in_comb_bin_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_bin_mc[flv_ind][bin_ind]/(cmb_val+part_val);
//    ev.f_cont_bin_mc         = ev.f_cont_in_comb_bin_mc*f_cmb;

//    sig_val  = NSigPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
//    cmb_val  = NCmbPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
//    part_val = NPrtPredicted_flv_mc[flv_mc_ind][bin_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

//    ev.f_bkg_flv_mc          = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
//    f_cmb                    = cmb_val/(cmb_val + part_val + sig_val);
//    ev.f_cont_in_comb_flv_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_flv_mc[flv_ind][bin_ind]/(cmb_val+part_val);
//    ev.f_cont_flv_mc         = ev.f_cont_in_comb_flv_mc*f_cmb;

//    sig_val  = NSigPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_sig->getVal(RooArgSet(*de,*mbc));
//    cmb_val  = NCmbPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_comb->getVal(RooArgSet(*de,*mbc));
//    part_val = NPrtPredicted_mc[flv_mc_ind][bin_mc_ind]*pdf_part->getVal(RooArgSet(*de,*mbc));

//    ev.f_bkg_mc          = (cmb_val + part_val)/(cmb_val + part_val + sig_val);
//    f_cmb                = cmb_val/(cmb_val + part_val + sig_val);
//    ev.f_cont_in_comb_mc = f_cont_in_comb*pdf_cont->getVal(RooArgSet(*de,*mbc))*NCmbPredicted_mc[flv_ind][bin_ind]/(cmb_val+part_val);
//    ev.f_cont_mc         = ev.f_cont_in_comb_mc*f_cmb;

    tree->Fill();
  }
  cout << "Writing a tree..." << endl;

  tree->Write();
  file->Close();
  return;
}

TTree* PurityFit::GetCPVTree(ICPVEv& ev){
  TTree* tree = new TTree("TEvent","TEvent");
  tree->Branch("exp",&ev.exp,"exp/I");
//  tree->Branch("run",&ev.run,"run/I");
//  tree->Branch("evtn",&ev.evtn,"evtn/I");
  tree->Branch("flv",&ev.flv,"flv/I");
  tree->Branch("bin",&ev.bin,"bin/I");

  tree->Branch("flv_mc",&ev.flv_mc,"flv_mc/I");
  tree->Branch("bin_mc",&ev.bin_mc,"bin_mc/I");

  tree->Branch("mp",&ev.mp,"mp/D");
  tree->Branch("mm",&ev.mm,"mm/D");
  tree->Branch("mp_mc",&ev.mp_mc,"mp_mc/D");
  tree->Branch("mm_mc",&ev.mm_mc,"mm_mc/D");

  tree->Branch("mode",&ev.mode,"mode/I");
  tree->Branch("h0mode",&ev.h0mode,"h0mode/I");
  if(m_mcflag) tree->Branch("b0f",&ev.b0f,"b0f/I");

  tree->Branch("costhBcms",&ev.costhBcms,"costhBcms/D");

  tree->Branch("z_sig",&ev.z_sig,"z_sig/D");
  tree->Branch("z_asc",&ev.z_asc,"z_asc/D");
  tree->Branch("sz_sig",&ev.sz_sig,"sz_sig/D");
  tree->Branch("sz_asc",&ev.sz_asc,"sz_asc/D");
  tree->Branch("chisq_sig",&ev.chisq_sig,"chisq_sig/D");
  tree->Branch("chisq_asc",&ev.chisq_asc,"chisq_asc/D");

  tree->Branch("z_sig_mc",&ev.z_sig_mc,"z_sig_mc/D");
  tree->Branch("z_asc_mc",&ev.z_asc_mc,"z_asc_mc/D");

  tree->Branch("ntrk_sig",&ev.ntrk_sig,"ntrk_sig/I");
  tree->Branch("ntrk_asc",&ev.ntrk_asc,"ntrk_asc/I");
  tree->Branch("ndf_sig",&ev.ndf_sig,"ndf_sig/I");
  tree->Branch("ndf_asc",&ev.ndf_asc,"ndf_asc/I");

  tree->Branch("f_bkg",&ev.f_bkg,"f_bkg/D");
  tree->Branch("f_cont",&ev.f_cont,"f_cont/D");
  tree->Branch("f_cont_in_comb",&ev.f_cont_in_comb,"f_cont_in_comb/D");

  tree->Branch("f_bkg_mc",&ev.f_bkg_mc,"f_bkg_mc/D");
  tree->Branch("f_cont_mc",&ev.f_cont_mc,"f_cont_mc/D");
  tree->Branch("f_cont_in_comb_mc",&ev.f_cont_in_comb_mc,"f_cont_in_comb_mc/D");

  tree->Branch("f_bkg_bin_mc",&ev.f_bkg_bin_mc,"f_bkg_bin_mc/D");
  tree->Branch("f_cont_bin_mc",&ev.f_cont_bin_mc,"f_cont_bin_mc/D");
  tree->Branch("f_cont_in_comb_bin_mc",&ev.f_cont_in_comb_bin_mc,"f_cont_in_comb_bin_mc/D");

  tree->Branch("f_bkg_flv_mc",&ev.f_bkg_flv_mc,"f_bkg_flv_mc/D");
  tree->Branch("f_cont_flv_mc",&ev.f_cont_flv_mc,"f_cont_flv_mc/D");
  tree->Branch("f_cont_in_comb_flv_mc",&ev.f_cont_in_comb_flv_mc,"f_cont_in_comb_flv_mc/D");

  tree->Branch("sigarea",&ev.sigarea,"sigarea/I");
  tree->Branch("de",&ev.de,"de/D");
  tree->Branch("mbc",&ev.mbc,"mbc/D");
  tree->Branch("tag_LH",&ev.tag_LH,"tag_LH/D");

  tree->Branch("psig",&ev.psig,"psig/D");
  tree->Branch("pcnt",&ev.pcnt,"pcnt/D");
  tree->Branch("pprt",&ev.pprt,"pprt/D");
  tree->Branch("pcmb",&ev.pcmb,"pcmb/D");
  tree->Branch("fbb",&ev.fbb,"fbb/D");
  tree->Branch("fbb_err",&ev.fbb_err,"fbb_err/D");
  tree->Branch("fprt",&ev.fprt,"fprt/D");
  tree->Branch("fprt_err",&ev.fprt_err,"fprt_err/D");
  tree->Branch("Nsig",&ev.Nsig,"Nsig/D");
  tree->Branch("Nsig_err",&ev.Nsig_err,"Nsig_err/D");

  tree->Branch("Ntot",&ev.NTot,"Ntot/I");
  tree->Branch("wrtag",&ev.wrtag,"wrtag/D");
  tree->Branch("Ntot_bin_mc",&ev.NTot_bin_mc,"Ntot_bin_mc/I");
  tree->Branch("wrtag_bin_mc",&ev.wrtag_bin_mc,"wrtag_bin_mc/D");
  tree->Branch("Ntot_flv_mc",&ev.NTot_flv_mc,"Ntot_flv_mc/I");
  tree->Branch("wrtag_flv_mc",&ev.wrtag_flv_mc,"wrtag_flv_mc/D");
  tree->Branch("Ntot_mc",&ev.NTot_mc,"Ntot_mc/I");
  tree->Branch("wrtag_mc",&ev.wrtag_mc,"wrtag_mc/D");

  return tree;
}

void PurityFit::CopyEvent(const ICPVEv& ev_from,ICPVEv& ev_to){
  ev_to.de        = ev_from.de;
  ev_to.mbc       = ev_from.mbc;
  ev_to.tag_LH    = ev_from.tag_LH;
  ev_to.costhBcms = ev_from.costhBcms;
  ev_to.flv       = ev_from.flv;
  ev_to.z_asc     = ev_from.z_asc;
  ev_to.z_sig     = ev_from.z_sig;
  ev_to.sz_asc    = ev_from.sz_asc;
  ev_to.sz_sig    = ev_from.sz_sig;
  ev_to.chisq_asc = ev_from.chisq_asc;
  ev_to.chisq_sig = ev_from.chisq_sig;
  ev_to.z_sig_mc  = ev_from.z_sig_mc;
  ev_to.z_asc_mc  = ev_from.z_asc_mc;
  ev_to.exp       = ev_from.exp;
  ev_to.bin       = ev_from.bin;
  ev_to.ntrk_asc  = ev_from.ntrk_asc;
  ev_to.ntrk_sig  = ev_from.ntrk_sig;
  ev_to.ndf_asc   = ev_from.ndf_asc;
  ev_to.ndf_sig   = ev_from.ndf_sig;
  ev_to.mode      = ev_from.mode;
  ev_to.h0mode    = ev_from.h0mode;
  ev_to.b0f       = ev_from.b0f;
  ev_to.mp        = ev_from.mp;
  ev_to.mm        = ev_from.mm;
  ev_to.mp_mc     = ev_from.mp_mc;
  ev_to.mm_mc     = ev_from.mm_mc;

  ev_to.psig      = ev_from.psig;
  ev_to.pcnt      = ev_from.pcnt;
  ev_to.pprt      = ev_from.pprt;
  ev_to.pcmb      = ev_from.pcmb;

  ev_to.fbb       = ev_from.fbb;
  ev_to.fbb_err   = ev_from.fbb_err;
  ev_to.fprt      = ev_from.fprt;
  ev_to.fprt_err  = ev_from.fprt_err;
  ev_to.Nsig      = ev_from.Nsig;
  ev_to.Nsig_err  = ev_from.Nsig_err;
  ev_to.NTot      = ev_from.NTot;
  ev_to.wrtag     = ev_from.wrtag;

  return;
}

void PurityFit::CopyTimeEvent(const ICPVEv& ev_from,ICPVEv& ev_to){
  ev_to.tag_LH    = ev_from.tag_LH;
  ev_to.costhBcms = ev_from.costhBcms;
  ev_to.z_asc     = ev_from.z_asc;
  ev_to.z_sig     = ev_from.z_sig;
  ev_to.sz_asc    = ev_from.sz_asc;
  ev_to.sz_sig    = ev_from.sz_sig;
  ev_to.chisq_asc = ev_from.chisq_asc;
  ev_to.chisq_sig = ev_from.chisq_sig;
  ev_to.z_sig_mc  = ev_from.z_sig_mc;
  ev_to.z_asc_mc  = ev_from.z_asc_mc;
  ev_to.ntrk_asc  = ev_from.ntrk_asc;
  ev_to.ntrk_sig  = ev_from.ntrk_sig;
  ev_to.ndf_asc   = ev_from.ndf_asc;
  ev_to.ndf_sig   = ev_from.ndf_sig;
  ev_to.b0f       = ev_from.b0f;
  return;
}

void PurityFit::FillEvent(const RooArgSet* aset, ICPVEv& ev, const int mctype){
  ev.de        = aset->getRealValue("de");
  ev.mbc       = aset->getRealValue("mbc");
  ev.tag_LH    = aset->getRealValue("tag_LH");
  ev.costhBcms = aset->getRealValue("costhBcms");
  ev.flv       = ev.tag_LH > 0 ? -1 : 1;
  ev.z_asc     = aset->getRealValue("z_asc");
  ev.z_sig     = aset->getRealValue("z_sig");
  ev.sz_asc    = aset->getRealValue("sz_asc");
  ev.sz_sig    = aset->getRealValue("sz_sig");
  ev.chisq_asc = aset->getRealValue("chisq_z_asc");
  ev.chisq_sig = aset->getRealValue("chisq_z_sig");

  ev.z_sig_mc  = mctype == 1 ? aset->getRealValue("z_sig_mc") : -99;
  ev.z_asc_mc  = mctype == 1 ? aset->getRealValue("z_asc_mc") : -99;
  ev.mp_mc     = mctype == 1 ? aset->getRealValue("mp_mc") : -99;
  ev.mm_mc     = mctype == 1 ? aset->getRealValue("mm_mc") : -99;
  ev.bin_mc    = mctype == 1 ? aset->getCatIndex("bin_mc") : 0;
  ev.flv_mc    = mctype == 1 ? aset->getCatIndex("flv_mc") : 0;

  ev.exp       = aset->getCatIndex("exp");
  ev.bin       = aset->getCatIndex("bin");
  ev.ntrk_asc  = aset->getCatIndex("ntrk_asc");
  ev.ntrk_sig  = aset->getCatIndex("ntrk_sig");
  ev.ndf_asc   = aset->getCatIndex("ndf_z_asc");
  ev.ndf_sig   = aset->getCatIndex("ndf_z_sig");
  ev.mode      = aset->getCatIndex("mode");
  ev.h0mode    = aset->getCatIndex("h0mode");
  ev.b0f       = m_mcflag ? aset->getCatIndex("b0f") : 0;
  ev.mp        = aset->getRealValue("mp");
  ev.mm        = aset->getRealValue("mm");
  return;
}

RooArgSet* PurityFit::GetArgSet(const int mctype, const bool include_bdt){
  cout << "Setting ArgSet... " << include_bdt << " " << m_svd << endl;
  RooArgSet* aset = new RooArgSet();
  aset = new RooArgSet();
  switch (m_svd) {
  case 1:
    aset->add(*exp_svd1);
    break;
  case 2:
    aset->add(*exp_svd2);
    break;
  default:
    aset->add(*exp);
    break;
  }
  if(m_mcflag){
    aset->add(*bin_mc);
    aset->add(*flv_mc);
  }
  if(mctype == 1){
    aset->add(*z_rec_mc);
    aset->add(*z_asc_mc);
    aset->add(*mp_mc);
    aset->add(*mm_mc);
  }
  aset->add(*mp);
  aset->add(*mm);
  aset->add(*b0f);
//  aset->add(*d0f);
//  aset->add(*h0f);
  aset->add(*rndm_pi0);
  aset->add(*mode); aset->add(*h0mode);
  aset->add(*flv);
  aset->add(*good_icpv);
  aset->add(*bin);
  aset->add(*de); aset->add(*mbc);
  aset->add(*md); aset->add(*mk); aset->add(*mh0);
  aset->add(*chi2_vtx_d0);
  aset->add(*costhBcms);
  if(!((m_mode == 2 || m_mode == 20 || m_mode == 5) && m_h0mode == 10)) aset->add(*mpi0);
  if(m_mode == 5 || m_mode == 50) aset->add(*dmetap);
  if(m_mode>9)                    aset->add(*dmdts0);

  aset->add(*pt_pip);
  aset->add(*pt_pim);
  aset->add(*z_pip);
  aset->add(*z_pim);
  aset->add(*r_pip);
  aset->add(*r_pim);

  if(m_h0mode == 20){
    aset->add(*pt_pi1);
    aset->add(*pt_pi2);
    aset->add(*p_pi0_h0);
    aset->add(*z_pi1);
    aset->add(*z_pi2);
    aset->add(*r_pi1);
    aset->add(*r_pi2);
  }
  if(m_mode == 3){
    aset->add(*cos_hel);
  }
  if(m_mode == 2 && m_h0mode == 10){
    aset->add(*e_g1);
  }

//  aset->add(*e_g2);
//  aset->add(*th_g1);
//  aset->add(*th_g2);
//  if(m_mode>9 || m_mode == 5){
//    aset->add(*e_g3);
//    aset->add(*e_g4);
//    aset->add(*th_g3);
//    aset->add(*th_g4);
//  }

  if(include_bdt){
    aset->add(*lh0);
    aset->add(*bdt);
  }
  aset->add(*tag_LH);
  aset->add(*ndf_asc); aset->add(*ntrk_asc);
  aset->add(*ndf_rec); aset->add(*ntrk_rec);
  aset->add(*z_rec); aset->add(*sz_rec); aset->add(*chisq_rec);
  aset->add(*z_asc); aset->add(*sz_asc); aset->add(*chisq_asc);
  cout << "done." << endl;
  return aset;
}

double PurityFit::de_line_size(void){
  const int numstr = 6./m_nstr;
  switch(m_mode){
  case 1:
    return !m_mcflag ? 60 : numstr*60;
  case 10:
    return !m_mcflag ? 17 : numstr*17;
  case 20:
    return !m_mcflag ? 3 : numstr*3;
  case 2:
    if(m_h0mode == 10) return !m_mcflag ? 15 : numstr*15;
    else               return !m_mcflag ? 8  : numstr*8;
  case 3:
    return !m_mcflag ? 40 : numstr*40;
  case 5:
    return !m_mcflag ? 4 : numstr*4;
  default:
    break;
  }
}

double PurityFit::mbc_line_size(void){
  const int numstr = 6./m_nstr;
  switch(m_mode){
  case 1:
    return !m_mcflag ? 85 : numstr*85;
  case 10:
    return !m_mcflag ? 25 : numstr*25;
  case 20:
    return !m_mcflag ? 5 : numstr*5;
  case 2:
    if(m_h0mode == 10) return !m_mcflag ? 25 : numstr*25;
    else               return !m_mcflag ? 13 : numstr*13;
  case 3:
    return !m_mcflag ? 45 : numstr*55;
  case 5:
    return !m_mcflag ? 6  : numstr*6;
  default:
    break;
  }
}
