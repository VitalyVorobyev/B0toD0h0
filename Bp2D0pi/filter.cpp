#include "../BtoDh/cuts.h"
#include "../DalitzModelStudy/phasespace.cpp"

int min_decision(const vector<double>& v){
    double vmin = v[0];
    double imin = 0;
    for(int i=1; i<v.size(); i++){
        if(v[i]<vmin){
            vmin = v[i];
            imin = i;
        }
    }
    return imin;
}

int max_decision(const vector<double>& v){
    double vmax = v[0];
    double imax = 0;
    for(int i=1; i<v.size(); i++){
        if(v[i]>vmax){
            vmax = v[i];
            imax = i;
        }
    }
    return imax;
}

void filter(const int type){
  TChain* tree = new TChain("TEvent");
  string fname;
  switch(type){
  case 0://Data
    tree->Add("/home/vitaly/B0toDh0/Tuples/Bp2D0pi/b2d0pip_data_v2.root");
    fname = string("fil_b2dpi_data_v2.root");
    break;
  case 1://Charged
    tree->Add("/home/vitaly/B0toDh0/Tuples/Bp2D0pi/b2dpi_charged_v2_0_10.root");
    fname = string("fil_b2dpi_charged_v2_0_10.root");
    break;
  case 2://Charm
    tree->Add("/home/vitaly/B0toDh0/Tuples/Bp2D0pi/b2dpi_v2_charm_0_10.root");
    fname = string("fil_b2dpi_charm_0_10_v2.root");
    break;
  case 3://Signal
    tree->Add("/home/vitaly/B0toDh0/Tuples/Bp2D0pi/bp2d0pip_sigmc.root");
    fname = string("fil_bp2d0pip_sigmc.root");
    break;
  default:
    return;
  }
    Double_t p_pi,p_ks,cos_thr,thr_sig,thr_oth;
    Int_t phsp,bin,exp,run,evtn,flv;
    Int_t bpf,d0f,pif,nptag,d0_flv_mc;
    Double_t k0mm2,k0et,k0hso00,k0hso02,k0hso04,k0hso10,k0hso12,k0hso14,k0hso20,k0hso22,k0hso24,k0hoo0,k0hoo1,k0hoo2,k0hoo3,k0hoo4;
    Double_t k1mm2,k1et,k1hso00,k1hso01,k1hso02,k1hso03,k1hso04,k1hso10,k1hso12,k1hso14,k1hso20,k1hso22,k1hso24,k1hoo0,k1hoo1,k1hoo2,k1hoo3,k1hoo4;
    Double_t chisq_vtx_d0,chisq_mass_d0;

    Double_t mbc,de,mp,mm,atckpi_max,atckpi_pi,mk,md;

    Double_t z_sig,z_asc;
    Double_t sz_sig,sz_asc;
    Double_t z_sig_d0, sz_sig_d0;
    Int_t ntrk_asc,ndf_z_sig,ndf_z_asc;
    Double_t chisq_z_sig,chisq_z_asc;
    Double_t cl_z_sig,cl_z_asc;
    Double_t costhBcms,Ecms;

    Double_t dz,dz_d0;

    Double_t mp_mc,mm_mc,d0_t_mc,z_mc_sig,z_mc_asc;
    Int_t bin_mc,flv_mc;

    tree->Print();

    if(type == 3){
      tree->SetBranchAddress("mp_mc",&mp_mc);
      tree->SetBranchAddress("mm_mc",&mm_mc);
      tree->SetBranchAddress("d0_t_mc",&d0_t_mc);
      tree->SetBranchAddress("z_sig_mc",&z_mc_sig);
      tree->SetBranchAddress("z_asc_mc",&z_mc_asc);
      tree->SetBranchAddress("bin_mc",&bin_mc);
      tree->SetBranchAddress("flv_mc",&flv_mc);
      tree->SetBranchAddress("d0_flv_mc",&d0_flv_mc);
    }

    tree->SetBranchAddress("exp",&exp);
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("evtn",&evtn);

    tree->SetBranchAddress("chisq_vtx_d0",&chisq_vtx_d0);
    tree->SetBranchAddress("chisq_mass_d0",&chisq_mass_d0);

    tree->SetBranchAddress("p_pi",&p_pi);
    tree->SetBranchAddress("p_ks",&p_ks);
    tree->SetBranchAddress("cos_thr",&cos_thr);
    tree->SetBranchAddress("thr_sig",&thr_sig);
    tree->SetBranchAddress("thr_oth",&thr_oth);
    tree->SetBranchAddress("phsp",&phsp);
    tree->SetBranchAddress("bin",&bin);
    tree->SetBranchAddress("flv",&flv);

    if(type){
      tree->SetBranchAddress("bpf",&bpf);
      tree->SetBranchAddress("d0f",&d0f);
      tree->SetBranchAddress("pif",&pif);
      tree->SetBranchAddress("nptag",&nptag);
    }
    tree->SetBranchAddress("mbc",&mbc);
    tree->SetBranchAddress("de",&de);
    tree->SetBranchAddress("mp",&mp);
    tree->SetBranchAddress("mm",&mm);
    tree->SetBranchAddress("atckpi_max",&atckpi_max);
    tree->SetBranchAddress("atckpi_pi",&atckpi_pi);
    tree->SetBranchAddress("mks_raw",&mk);
    tree->SetBranchAddress("md0_raw",&md);

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
    tree->SetBranchAddress("k1hso01",&k1hso01);
    tree->SetBranchAddress("k1hso02",&k1hso02);
    tree->SetBranchAddress("k1hso03",&k1hso03);
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

    TFile *ofile = new TFile(fname.c_str(),"RECREATE");
    fname = string("TEvent");
    TTree* TEventTr = new TTree(fname.c_str(),fname.c_str());

    TEventTr->Branch("exp",&exp,"exp/I");
    TEventTr->Branch("run",&run,"run/I");
    TEventTr->Branch("evtn",&evtn,"evtn/I");

    if(type == 3){
      TEventTr->Branch("mp_mc",&mp_mc,"mp_mc/D");
      TEventTr->Branch("mm_mc",&mm_mc,"mp_mc/D");
      TEventTr->Branch("d0_t_mc",&d0_t_mc,"d0_t_mc/D");
      TEventTr->Branch("z_sig_mc",&z_mc_sig,"z_mc_sig/D");
      TEventTr->Branch("z_asc_mc",&z_mc_asc,"z_mc_asc/D");
      TEventTr->Branch("bin_mc",&bin_mc,"bin_mv/I");
      TEventTr->Branch("flv_mc",&flv_mc,"flv_mc/I");
      TEventTr->Branch("d0_flv_mc",&d0_flv_mc,"d0_flv_mc/I");
    }

    TEventTr->Branch("chisq_vtx_d0",&chisq_vtx_d0,"chisq_vtx_d0/D");
    TEventTr->Branch("chisq_mass_d0",&chisq_mass_d0,"chisq_mass_d0/D");

    TEventTr->Branch("p_pi",&p_pi,"p_pi/D");
    TEventTr->Branch("p_ks",&p_ks,"p_ks/D");
    TEventTr->Branch("cos_thr",&cos_thr,"cos_thr/D");
    TEventTr->Branch("thr_sig",&thr_sig,"thr_sig/D");
    TEventTr->Branch("thr_oth",&thr_oth,"thr_oth/D");
    TEventTr->Branch("thr_oth",&thr_oth,"thr_oth/D");
    TEventTr->Branch("phsp",&phsp,"phsp/I");
    TEventTr->Branch("bin",&bin,"bin/I");
    TEventTr->Branch("flv",&flv,"flv/I");

    if(type){
      TEventTr->Branch("bpf",&bpf,"bpf/I");
      TEventTr->Branch("d0f",&d0f,"d0f/I");
      TEventTr->Branch("pif",&pif,"pif/I");
      TEventTr->Branch("nptag",&nptag,"nptag/I");
    }
    TEventTr->Branch("mbc",&mbc,"mbc/D");
    TEventTr->Branch("de",&de,"de/D");
    TEventTr->Branch("mp",&mp,"mp/D");
    TEventTr->Branch("mm",&mm,"mm/D");
    TEventTr->Branch("dz",&dz,"dz/D");
    TEventTr->Branch("dz_d0",&dz_d0,"dz_d0/D");
    TEventTr->Branch("atckpi_max",&atckpi_max,"atckpi_max/D");
    TEventTr->Branch("atckpi_pi",&atckpi_pi,"atckpi_pi/D");
    TEventTr->Branch("mk",&mk,"mk/D");
    TEventTr->Branch("md",&md,"md/D");

    TEventTr->Branch("z_sig",&z_sig,"z_sig/D");
    TEventTr->Branch("z_sig_d0",&z_sig_d0,"z_sig_d0/D");
    TEventTr->Branch("z_asc",&z_asc,"z_asc/D");

    TEventTr->Branch("sz_sig",&sz_sig,"sz_sig/D");
    TEventTr->Branch("sz_sig_d0",&sz_sig_d0,"sz_sig_d0/D");
    TEventTr->Branch("sz_asc",&sz_asc,"sz_asc/D");

    TEventTr->Branch("ntrk_asc",&ntrk_asc,"ntrk_asc/I");
    TEventTr->Branch("ndf_z_sig",&ndf_z_sig,"ndf_z_sig/I");
    TEventTr->Branch("ndf_z_asc",&ndf_z_asc,"ndf_z_asc/I");

    TEventTr->Branch("chisq_z_sig",&chisq_z_sig,"chisq_z_sig/D");
    TEventTr->Branch("chisq_z_asc",&chisq_z_asc,"chisq_z_asc/D");
    TEventTr->Branch("cl_z_sig",&cl_z_sig,"cl_z_sig/D");
    TEventTr->Branch("cl_z_asc",&cl_z_asc,"cl_z_asc/D");

    TEventTr->Branch("costhBcms",&costhBcms,"costhBcms/D");
    TEventTr->Branch("Ecms",&Ecms,"Ecm/D");

    TEventTr->Branch("k0mm2",&k0mm2,"k0mm2/D");
    TEventTr->Branch("k0et",&k0et,"k0et/D");
    TEventTr->Branch("k0hso00",&k0hso00,"k0hso00/D");
    TEventTr->Branch("k0hso02",&k0hso02,"k0hso02/D");
    TEventTr->Branch("k0hso04",&k0hso04,"k0hso04/D");
    TEventTr->Branch("k0hso10",&k0hso10,"k0hso10/D");
    TEventTr->Branch("k0hso12",&k0hso12,"k0hso12/D");
    TEventTr->Branch("k0hso14",&k0hso14,"k0hso14/D");
    TEventTr->Branch("k0hso20",&k0hso20,"k0hso20/D");
    TEventTr->Branch("k0hso22",&k0hso22,"k0hso22/D");
    TEventTr->Branch("k0hso24",&k0hso24,"k0hso24/D");
    TEventTr->Branch("k0hoo0",&k0hoo0,"k0hoo0/D");
    TEventTr->Branch("k0hoo1",&k0hoo1,"k0hoo1/D");
    TEventTr->Branch("k0hoo2",&k0hoo2,"k0hoo2/D");
    TEventTr->Branch("k0hoo3",&k0hoo3,"k0hoo3/D");
    TEventTr->Branch("k0hoo4",&k0hoo4,"k0hoo4/D");

    TEventTr->Branch("k1mm2",&k1mm2,"k1mm2/D");
    TEventTr->Branch("k1et",&k1et,"k1et/D");
    TEventTr->Branch("k1hso00",&k1hso00,"k1hso00/D");
    TEventTr->Branch("k1hso01",&k1hso01,"k1hso01/D");
    TEventTr->Branch("k1hso02",&k1hso02,"k1hso02/D");
    TEventTr->Branch("k1hso03",&k1hso03,"k1hso03/D");
    TEventTr->Branch("k1hso04",&k1hso04,"k1hso04/D");
    TEventTr->Branch("k1hso10",&k1hso10,"k1hso10/D");
    TEventTr->Branch("k1hso12",&k1hso12,"k1hso12/D");
    TEventTr->Branch("k1hso14",&k1hso14,"k1hso14/D");
    TEventTr->Branch("k1hso20",&k1hso20,"k1hso20/D");
    TEventTr->Branch("k1hso22",&k1hso22,"k1hso22/D");
    TEventTr->Branch("k1hso24",&k1hso24,"k1hso24/D");
    TEventTr->Branch("k1hoo0",&k1hoo0,"k1hoo0/D");
    TEventTr->Branch("k1hoo1",&k1hoo1,"k1hoo1/D");
    TEventTr->Branch("k1hoo2",&k1hoo2,"k1hoo2/D");
    TEventTr->Branch("k1hoo3",&k1hoo3,"k1hoo3/D");
    TEventTr->Branch("k1hoo4",&k1hoo4,"k1hoo4/D");

    int good_icpv_mlt, good_icpv_sgl;
    TEventTr->Branch("good_icpv_mlt",&good_icpv_mlt,"good_icpv_mlt/I");
    TEventTr->Branch("good_icpv_sgl",&good_icpv_sgl,"good_icpv_sgl/I");

    PhaseSpace PHSP;
    const int NTot = tree->GetEntries();
    for(int i=0; i<NTot; i++){
      if(!(i%10000)){ cout << i << " events" << endl;}
      tree->GetEvent(i);
      if(type && (bpf<-1 || !bpf)) continue;
      if(md<1.85155 || md>1.87833) continue;
//      if(md<(1.865-0.013) || md>(1.865+0.013)) continue;
      if(mk<(0.4975-0.009) || mk>(0.4975+0.009)) continue;
//      if(chisq_z_sig>1000 || chisq_z_sig<0) continue;
//      if(chisq_z_asc>1000 || chisq_z_asc<0) continue;
      if(chisq_vtx_d0>500 || chisq_vtx_d0<0) continue;
//      if(chisq_mass_d0>1000 || chisq_mass_d0<0) continue;

      sz_sig    = 10.*TMath::Sqrt(sz_sig);
      sz_sig_d0 = 10.*TMath::Sqrt(sz_sig_d0);
      sz_asc    = 10.*TMath::Sqrt(sz_asc);
      z_sig    *= 10;
      z_asc    *= 10;
      z_sig_d0 *= 10;
      dz    = z_sig - z_asc;
      dz_d0 = z_sig_d0 - z_asc;

      // * Standatd ICPV cuts * //
      if(ndf_z_asc == 0 && sz_asc>0.5 || ndf_z_asc > 0 && (sz_asc>0.2 || chisq_z_asc/ndf_z_asc>50)){
        good_icpv_mlt = 0;
        good_icpv_sgl = 0;
      } else if(ndf_z_sig == 0 && sz_sig>0.5 || ndf_z_sig > 0 && (sz_sig>0.2 || chisq_z_sig/ndf_z_sig>50)){
        good_icpv_mlt = 0;
        good_icpv_sgl = 0;
      } else{
        good_icpv_mlt = 1;
        good_icpv_sgl = 1;
      }
      if(sz_sig_d0>0.5) good_icpv_sgl = 0;
      // * ////////////////// * //

      const int new_bin = PHSP.GetBin(mp,mm);
      bin = new_bin*bin>0 ? new_bin : -new_bin;
      TEventTr->Fill();
    }
    TEventTr->Print();
    TEventTr->Write();
    ofile->Close();
    return;
}
