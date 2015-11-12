#include "../rooksfw/ksfw-new-sl1.h"
#include "../rooksfw/rooksfw.h"
#include "../rooksfw/rooksfw.cc"

using namespace std;
void tmva_reader(const int type = 1){
// type = 1 -> Charged
// type = 2 -> Charm
// type = 0 -> Data
/////////////////
// TMVA Reader //
/////////////////
  TMVA::Reader* reader  = new TMVA::Reader("!Color:!Silent:V");

  TString fname;
  TString* tmststr;
  TList* treenames = new TList();
  TString fname;
  string ofname;
  switch(type){
    case 0:// Data
      fname = TString("fil_b2dpi_data_v2.root");
      ofname = string("Fil_b2dpi_data_v2.root");
      break;
    case 1:// Charged
      fname = TString("fil_b2dpi_charged_v2_0_10.root");
      ofname = string("Fil_b2dpi_charged_v2_0_10.root");
      break;
    case 2:// Charm
      fname = TString("fil_b2dpi_charm_0_10_v2.root");
      ofname = string("Fil_b2dpi_charm_v2_0_10.root");
      break;
    case 3:// Signal
      fname = TString("fil_bp2d0pip_sigmc.root");
      ofname = string("Fil_bp2d0pip_sigmc.root");
      break;
    default:
      cout << "Wrong type " << type << endl;
      return;
  }

  TFile ofile(ofname.c_str(),"RECREATE");
  TTree *TEvent = new TTree("TEvent","TEvent");

  Double_t chisq_vtx_d0,chisq_mass_d0,cos_thr,thr_sig,thr_oth;
  Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;
  Double_t k1mm2,k1et,k1hso00,k1hso02,k1hso04,k1hso10,k1hso12,k1hso14,k1hso20,k1hso22,k1hso24,k1hoo0,k1hoo1,k1hoo2,k1hoo3,k1hoo4;

  Double_t mbc,de,mp,mm,dz,dz_d0,atckpi_max,atckpi_pi,mk,md;
  Double_t bdtg;
  Int_t phsp,bin,exp,run,evtn;
  Double_t p_pi;
  Int_t good_icpv_mlt, good_icpv_sgl;

  Float_t m_costhBcms,m_chisq_vtx_d0,m_chisq_mass_d0,m_cos_thr,m_thr_sig,m_thr_oth, m_p_pi;
  Float_t m_k1mm2,m_k1et,m_k1hso00,m_k1hso02,m_k1hso04,m_k1hso10,m_k1hso12,m_k1hso14,m_k1hso20,m_k1hso22,m_k1hso24,m_k1hoo0,m_k1hoo1,m_k1hoo2,m_k1hoo3,m_k1hoo4;

  Int_t bpf,d0f,pif;
  Int_t flv;

  Double_t z_sig,z_sig_d0,z_asc;
  Double_t sz_sig,sz_sig_d0,sz_asc;
  Int_t ntrk_asc,ndf_z_sig,ndf_z_asc;
  Double_t chisq_z_sig,chisq_z_asc,cl_z_sig,cl_z_asc;
  Double_t costhBcms,Ecms;
  reader->AddVariable("abs(costhBcms)",&m_costhBcms);
  reader->AddVariable("log(chisq_vtx_d0)",&m_chisq_vtx_d0);
  reader->AddVariable("log(chisq_mass_d0)",&m_chisq_mass_d0);
  reader->AddVariable("abs(cos_thr)",&m_cos_thr);
  reader->AddVariable("thr_sig",&m_thr_sig);
  reader->AddVariable("thr_oth",&m_thr_oth);
  reader->AddVariable("p_pi",&m_p_pi);

  reader->AddVariable("k1mm2",&m_k1mm2);
  reader->AddVariable("k1et",&m_k1et);
  reader->AddVariable("k1hso00",&m_k1hso00);
  reader->AddVariable("k1hso02",&m_k1hso02);
  reader->AddVariable("k1hso04",&m_k1hso04);
  reader->AddVariable("k1hso10",&m_k1hso10);
  reader->AddVariable("k1hso12",&m_k1hso12);
  reader->AddVariable("k1hso14",&m_k1hso14);
  reader->AddVariable("k1hso20",&m_k1hso20);
  reader->AddVariable("k1hso22",&m_k1hso22);
  reader->AddVariable("k1hso24",&m_k1hso24);
  reader->AddVariable("k1hoo0",&m_k1hoo0);
  reader->AddVariable("k1hoo1",&m_k1hoo1);
  reader->AddVariable("k1hoo2",&m_k1hoo2);
  reader->AddVariable("k1hoo3",&m_k1hoo3);
  reader->AddVariable("k1hoo4",&m_k1hoo4);

  reader->BookMVA("Bp2D0pi","../TMVA/weights/MVA_Bp2D0pi_BDTG.weights.xml");

  TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree*)ifile->Get("TEvent");

  Int_t nptag,bin_mc,flv_asc,d0_flv_mc;
  Double_t mp_mc,mm_mc,d0_t_mc,z_mc_sig,z_mc_asc;

  if(type){
    tree->SetBranchAddress("bpf",&bpf);
    tree->SetBranchAddress("d0f",&d0f);
    tree->SetBranchAddress("pif",&pif);
    tree->SetBranchAddress("nptag",&nptag);
  }

  if(type == 3){
    tree->SetBranchAddress("mp_mc",&mp_mc);
    tree->SetBranchAddress("mm_mc",&mm_mc);
    tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
    tree->SetBranchAddress("z_sig_mc",&z_mc_sig);
    tree->SetBranchAddress("z_asc_mc",&z_mc_asc);
    tree->SetBranchAddress("bin_mc",&bin_mc);
    tree->SetBranchAddress("flv_mc",&flv_asc);
    tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
  }

  tree->SetBranchAddress("good_icpv_mlt",&good_icpv_mlt);
  tree->SetBranchAddress("good_icpv_sgl",&good_icpv_sgl);

  tree->SetBranchAddress("exp",&exp);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evtn",&evtn);
  tree->SetBranchAddress("mbc",&mbc);
  tree->SetBranchAddress("de",&de);
  tree->SetBranchAddress("mp",&mp);
  tree->SetBranchAddress("mm",&mm);
  tree->SetBranchAddress("bin",&bin);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("dz_d0",&dz_d0);
  tree->SetBranchAddress("atckpi_max",&atckpi_max);
  tree->SetBranchAddress("atckpi_pi",&atckpi_pi);
  tree->SetBranchAddress("phsp",&phsp);
  tree->SetBranchAddress("flv",&flv);

  tree->SetBranchAddress("mk",&mk);
  tree->SetBranchAddress("md",&md);

  tree->SetBranchAddress("z_sig",&z_sig);
  tree->SetBranchAddress("z_sig_d0",&z_sig_d0);
  tree->SetBranchAddress("z_asc",&z_asc);

  tree->SetBranchAddress("sz_sig",&sz_sig);
  tree->SetBranchAddress("sz_sig_d0",&sz_sig_d0);
  tree->SetBranchAddress("sz_asc",&sz_asc);

  tree->SetBranchAddress("ntrk_asc",&ntrk_asc);
  tree->SetBranchAddress("ndf_z_sig",&ndf_z_sig);
  tree->SetBranchAddress("ndf_z_asc",&ndf_z_asc);

  tree->SetBranchAddress("chisq_z_sig",&chisq_z_sig);
  tree->SetBranchAddress("chisq_z_asc",&chisq_z_asc);
  tree->SetBranchAddress("cl_z_sig",&cl_z_sig);
  tree->SetBranchAddress("cl_z_asc",&cl_z_asc);

  tree->SetBranchAddress("costhBcms",&costhBcms);
  tree->SetBranchAddress("Ecms",&Ecms);

  tree->SetBranchAddress("p_pi",&p_pi);
  tree->SetBranchAddress("chisq_vtx_d0",&chisq_vtx_d0);
  tree->SetBranchAddress("chisq_mass_d0",&chisq_mass_d0);

  tree->SetBranchAddress("cos_thr",&cos_thr);
  tree->SetBranchAddress("thr_sig",&thr_sig);
  tree->SetBranchAddress("thr_oth",&thr_oth);

  tree->SetBranchAddress("k0mm2",&k0mm2);
  tree->SetBranchAddress("k0et",&k0et);
  tree->SetBranchAddress("k0hso00",&k0hso00);
  tree->SetBranchAddress("k0hso02",&k0hso02);
  tree->SetBranchAddress("k0hso04",&k0hso04);
  tree->SetBranchAddress("k0hso10",&k0hso10);
  tree->SetBranchAddress("k0hso12",&k0hso12);
  tree->SetBranchAddress("k0hso14",&k0hso14);
  tree->SetBranchAddress("k0hso20",&k0hso20);
  tree->SetBranchAddress("k0hso22",&k0hso22);
  tree->SetBranchAddress("k0hso24",&k0hso24);
  tree->SetBranchAddress("k0hoo0",&k0hoo0);
  tree->SetBranchAddress("k0hoo1",&k0hoo1);
  tree->SetBranchAddress("k0hoo2",&k0hoo2);
  tree->SetBranchAddress("k0hoo3",&k0hoo3);
  tree->SetBranchAddress("k0hoo4",&k0hoo4);

  tree->SetBranchAddress("k1mm2",&k1mm2);
  tree->SetBranchAddress("k1et",&k1et);
  tree->SetBranchAddress("k1hso00",&k1hso00);
  tree->SetBranchAddress("k1hso02",&k1hso02);
  tree->SetBranchAddress("k1hso04",&k1hso04);
  tree->SetBranchAddress("k1hso10",&k1hso10);
  tree->SetBranchAddress("k1hso12",&k1hso12);
  tree->SetBranchAddress("k1hso14",&k1hso14);
  tree->SetBranchAddress("k1hso20",&k1hso20);
  tree->SetBranchAddress("k1hso22",&k1hso22);
  tree->SetBranchAddress("k1hso24",&k1hso24);
  tree->SetBranchAddress("k1hoo0",&k1hoo0);
  tree->SetBranchAddress("k1hoo1",&k1hoo1);
  tree->SetBranchAddress("k1hoo2",&k1hoo2);
  tree->SetBranchAddress("k1hoo3",&k1hoo3);
  tree->SetBranchAddress("k1hoo4",&k1hoo4);

  TEvent->Branch("good_icpv_mlt",&good_icpv_mlt,"good_icpv_mlt/I");
  TEvent->Branch("good_icpv_sgl",&good_icpv_sgl,"good_icpv_sgl/I");

  TEvent->Branch("exp",&exp,"exp/I");
  TEvent->Branch("run",&run,"run/I");
  TEvent->Branch("evtn",&evtn,"evtn/I");

  if(type){
    TEvent->Branch("bpf",&bpf,"bpf/I");
    TEvent->Branch("d0f",&d0f,"d0f/I");
    TEvent->Branch("pif",&pif,"pif/I");
    TEvent->Branch("nptag",&nptag,"nptag/I");
  }

  if(type == 3){
    TEvent->Branch("mp_mc",&mp_mc,"mp_mc/D");
    TEvent->Branch("mm_mc",&mm_mc,"mp_mc/D");
    TEvent->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
    TEvent->Branch("z_sig_mc",&z_mc_sig,"z_mc_sig/D");
    TEvent->Branch("z_asc_mc",&z_mc_asc,"z_mc_asc/D");
    TEvent->Branch("bin_mc",&bin_mc,"bin_mv/I");
    TEvent->Branch("flv_mc",&flv_asc,"flv_mc/I");
    TEvent->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
  }

  TEvent->Branch("bdtg",&bdtg,"bdtg/D");

  TEvent->Branch("mp",&mp,"mp/D");
  TEvent->Branch("mm",&mm,"mm/D");
  TEvent->Branch("bin",&bin,"bin/I");
  TEvent->Branch("flv",&flv,"flv/I");

  TEvent->Branch("mbc",&mbc,"mbc/D");
  TEvent->Branch("de",&de,"de/D");

  TEvent->Branch("mk",&mk,"mk/D");
  TEvent->Branch("md",&md,"md/D");

  TEvent->Branch("dz",&dz,"dz/D");
  TEvent->Branch("dz_d0",&dz_d0,"dz_d0/D");

  TEvent->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
  TEvent->Branch("atckpi_pi",&atckpi_pi,"atckpi_pi/D");

  TEvent->Branch("z_sig",&z_sig,"z_sig/D");
  TEvent->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
  TEvent->Branch("z_asc",&z_asc,"z_asc/D");

  TEvent->Branch("sz_sig",&sz_sig,"sz_sig/D");
  TEvent->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
  TEvent->Branch("sz_asc",&sz_asc,"sz_asc/D");

  TEvent->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
  TEvent->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
  TEvent->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");

  TEvent->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
  TEvent->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
  TEvent->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
  TEvent->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");

  TEvent->Branch("costhBcms",&costhBcms,"costhBcms/D");
  TEvent->Branch("Ecms",&Ecms,"Ecm/D");

  const int NTot = tree->GetEntries();
  for(int i=0; i<NTot; i++){
    if(!(i%10000)) cout << i << " events" << endl;
    tree->GetEvent(i);

    m_costhBcms   = (float)TMath::Abs(costhBcms);
    m_chisq_vtx_d0= (float)TMath::Log(chisq_vtx_d0);
    m_chisq_mass_d0= (float)TMath::Log(chisq_mass_d0);
    m_cos_thr     = (float)TMath::Abs(cos_thr);
    m_thr_sig     = (float)thr_sig;
    m_thr_oth     = (float)thr_oth;
    m_p_pi        = (float)p_pi;
    m_k1mm2       = (float)k1mm2;
    m_k1et        = (float)k1et;
    m_k1hso00     = (float)k1hso00;
    m_k1hso02     = (float)k1hso02;
    m_k1hso04     = (float)k1hso04;
    m_k1hso10     = (float)k1hso10;
    m_k1hso12     = (float)k1hso12;
    m_k1hso14     = (float)k1hso14;
    m_k1hso20     = (float)k1hso20;
    m_k1hso22     = (float)k1hso22;
    m_k1hso24     = (float)k1hso24;
    m_k1hoo0      = (float)k1hoo0;
    m_k1hoo1      = (float)k1hoo1;
    m_k1hoo2      = (float)k1hoo2;
    m_k1hoo3      = (float)k1hoo3;
    m_k1hoo4      = (float)k1hoo4;

    bdtg = reader->EvaluateMVA("Bp2D0pi");

    TEvent->Fill();
  }

  TEvent->Write();
  ofile.Write();
  ofile.Close();
}
