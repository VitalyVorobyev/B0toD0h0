using namespace RooFit;
using namespace std;

void draw_smth(RooDataSet* ds, RooRealVar& var,string& label, string& cut, string& title,const double& ledge, const double& redge, const bool log, const int top){
  RooPlot* Frame = var.frame(Title(title.c_str()));
  ds->plotOn(Frame,DataError(RooAbsData::SumW2),MarkerSize(1),Cut(cut.c_str()));

  TCanvas* cm = new TCanvas(label.c_str(),"cm",600,400);
  cm->cd();
  ds->statOn(Frame,Layout(0.6,0.98,0.9),Cut(cut.c_str()));
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.00,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->SetGrid();
  pad1->SetLeftMargin(0.15);
  pad1->SetFillColor(0);
  if(log) pad1->SetLogy();
  Frame->GetXaxis()->SetTitleSize(0.05);
  Frame->GetXaxis()->SetTitleOffset(0.85);
  Frame->GetXaxis()->SetLabelSize(0.05);
  Frame->GetYaxis()->SetTitleOffset(1.6);
  Frame->Draw();
  TLine *lineLEFT = new TLine(ledge,0,ledge,top);
  lineLEFT->SetLineColor(kRed);
  lineLEFT->SetLineWidth(2);
  lineLEFT->SetLineStyle(1);
  lineLEFT->Draw();
  TLine *lineRIGHT = new TLine(redge,0,redge,top);
  lineRIGHT->SetLineColor(kRed);
  lineRIGHT->SetLineWidth(2);
  lineRIGHT->SetLineStyle(1);
  lineRIGHT->Draw();
  cm->Draw();

//  TPaveText *pt = new TPaveText(0.7,0.6,0.98,0.7,"brNDC");
//  pt->SetFillColor(0);
//  pt->SetTextAlign(12);
//  out.str("");
//  out << "Eff = " << std::setprecision(2) << int_m_etagg_sig/int_m_etagg_plot;
//  pt->AddText(out.str().c_str());
//  pt->Draw();

  stringstream out;
  out.str("");
  out << "pics/" << label << ".eps";
  cm->Print(out.str().c_str());
  out.str("");
  out << "pics/" << label << ".root";
  cm->Print(out.str().c_str());
  return;
}

void draw_rzhel(void){
  TChain* tree_sig_pi0   = new TChain("TEvent");
  TChain* tree_sig_omega = new TChain("TEvent");
  TChain* tree_bkg       = new TChain("TEvent");

  vector<int> streams;
  streams.push_back(0);
  streams.push_back(1);
  streams.push_back(2);
  streams.push_back(3);
  streams.push_back(4);
  streams.push_back(5);

  vector< string > types;
  types.push_back(string("uds"));
  types.push_back(string("charm"));
  types.push_back(string("charged"));
  types.push_back(string("mixed"));

  const string prefix("/home/vitaly/B0toDh0/Tuples/");
  stringstream out;
  out.str("");
  out << prefix << "Fil_b2dh_sigmcPi0_s8.root";
  tree_sig_pi0->Add(out.str().c_str());

  out.str("");
  out << prefix << "Fil_b2dh_sigmcOMEGA_s6.root";
  tree_sig_omega->Add(out.str().c_str());

  for(int i=0; i<streams.size(); i++){
    for(int j=0; j<types.size(); j++){
      out.str("");
      out << prefix <<"/Fil_b2dh_" << types[j] << "_" << streams[i] << "_" << streams[i]+10 << ".root";
      tree_bkg->Add(out.str().c_str());
    }
  }

  RooArgSet pi0_set;
  RooArgSet omega_set;

  RooRealVar de("de","de",-0.15,0.3,"GeV");
  de.setRange("pi0_signal",-0.1,0.0921343);
  de.setRange("omega_signal",-0.0565915,0.0467748);

  RooRealVar mbc("mbc","mbc",5.2,5.3,"GeV");
  de.setRange("pi0_signal",5.27093,5.28789);
  de.setRange("omega_signal",5.27198,5.28743);

  RooRealVar r_pip("r_pip","r_pip",-5,5,"mm");
  RooRealVar r_pim("r_pim","r_pim",-5,5,"mm");
  RooRealVar r_pi1("r_pi1","r_pi1",-5,5,"mm");
  RooRealVar r_pi2("r_pi2","r_pi2",-5,5,"mm");

  RooRealVar z_pip("z_pip","z_pip",-10,10,"mm");
  RooRealVar z_pim("z_pim","z_pim",-10,10,"mm");
  RooRealVar z_pi1("z_pi1","z_pi1",-10,10,"mm");
  RooRealVar z_pi2("z_pi2","z_pi2",-10,10,"mm");

  RooRealVar pt_pip("pt_pip","pt_pip",0,2.3,"GeV");
  RooRealVar pt_pim("pt_pim","pt_pim",0,2.3,"GeV");
  RooRealVar pt_pi1("pt_pi1","pt_pi1",0,2.3,"GeV");
  RooRealVar pt_pi2("pt_pi2","pt_pi2",0,2.3,"GeV");

  RooRealVar hel_pi0("hel_pi0","hel_pi0",-1.,1.);
  RooRealVar hel_h0("hel_h0","hel_h0",-1.,1.);

  RooRealVar th_g1("th_g1","Cos(#theta_{#gamma_{1}}})",-1.,1.,"GeV");
  th_g1.setRange("forward_g1",0.6,1.);
  th_g1.setRange("barrel_g1",-1,0.6);
  RooRealVar th_g2("th_g2","Cos(#theta_{#gamma_{2}}})",-1.,1.,"GeV");
  th_g2.setRange("forward_g2",0.6,1.);
  th_g2.setRange("barrel_g2",-1,0.6);

  RooRealVar e_g1("e_g1","E(#gamma_{1})",0,4.5,"GeV");
  RooRealVar e_g2("e_g2","E(#gamma_{2})",0,4.5,"GeV");

  RooCategory mode("mode","mode");
  mode.defineType("pi0",1);
  mode.defineType("eta",2);
  mode.defineType("omega",3);

  RooCategory h0mode("h0mode","h0mode");
  h0mode.defineType("gg",10);
  h0mode.defineType("ppp",20);

  RooCategory b0f("b0f","b0f");
  b0f.defineType("sig1",1);
  b0f.defineType("sig5",5);
  b0f.defineType("sig10",10);
  b0f.defineType("cmb",-1);
  b0f.defineType("part3",3);
  b0f.defineType("part4",4);

  pi0_set.add(de);
  pi0_set.add(mbc);
  pi0_set.add(r_pip);
  pi0_set.add(r_pim);
  pi0_set.add(z_pip);
  pi0_set.add(z_pim);
  pi0_set.add(pt_pip);
  pi0_set.add(pt_pim);
  pi0_set.add(hel_h0);
  pi0_set.add(e_g1);
  pi0_set.add(e_g2);
  pi0_set.add(th_g1);
  pi0_set.add(th_g2);
  pi0_set.add(mode);
  pi0_set.add(h0mode);
  pi0_set.add(b0f);

  omega_set.add(de);
  omega_set.add(mbc);
  omega_set.add(r_pip);
  omega_set.add(r_pim);
  omega_set.add(r_pi1);
  omega_set.add(r_pi2);
  omega_set.add(z_pip);
  omega_set.add(z_pim);
  omega_set.add(z_pi1);
  omega_set.add(z_pi2);
  omega_set.add(pt_pip);
  omega_set.add(pt_pim);
  omega_set.add(pt_pi1);
  omega_set.add(pt_pi2);
  omega_set.add(hel_pi0);
  omega_set.add(e_g1);
  omega_set.add(e_g2);
  omega_set.add(th_g1);
  omega_set.add(th_g2);
  omega_set.add(mode);
  omega_set.add(h0mode);
  omega_set.add(b0f);

  RooDataSet* pi0_sig_ds   = new RooDataSet("pi0_sig_ds","pi0_sig_ds",tree_sig_pi0,pi0_set,"mode == 1 && (b0f == 1 || b0f == 5 || b0f == 10) && de>-0.1 && de<0.0921343 && mbc>5.27093 && mbc<5.28789");
  RooDataSet* omega_sig_ds = new RooDataSet("omega_sig_ds","omega_sig_ds",tree_sig_omega,omega_set,"mode == 3 && (b0f == 1 || b0f == 5 || b0f == 10) && de>-0.1 && de<0.0921343 && mbc>5.27093 && mbc<5.28789");

  RooDataSet* pi0_bkg_ds   = new RooDataSet("pi0_bkg_ds","pi0_bkg_ds",tree_bkg,pi0_set,"mode == 1 && !(b0f == 1 || b0f == 5 || b0f == 10) && de>-0.0565915 && de<0.0467748 && mbc>5.27198 && mbc<5.28743");
  RooDataSet* omega_bkg_ds = new RooDataSet("omega_bkg_ds","omega_bkg_ds",tree_bkg,omega_set,"mode == 3 && !(b0f == 1 || b0f == 5 || b0f == 10) && de>-0.0565915 && de<0.0467748 && mbc>5.27198 && mbc<5.28743");
  RooDataSet* etagg_bkg_ds = new RooDataSet("etagg_bkg_ds","etagg_bkg_ds",tree_bkg,pi0_set,"mode == 2 && h0mode == 10 && !(b0f == 1 || b0f == 5 || b0f == 10) && abs(de)<0.1 && mbc>5.27");

  const string g1_frw_cut("th_g1>0.6");
  const string g1_bar_cut("th_g1<0.6");
  const string g2_frw_cut("th_g2>0.6");
  const string g2_bar_cut("th_g2<0.6");
  // eta -> gg //
  draw_smth(etagg_bkg_ds,hel_h0,string("hel_h0_bkg_etagg"),string(""),string("(E_{1}-E_{2})/(E_{1}+E_{2}), bkg #pi^{0}"),-1,1,false,220);
  draw_smth(etagg_bkg_ds,e_g1,string("e_g1_bkg_etagg_frw"),g1_frw_cut,string("Photon 1 energy (forward), bkg #eta#to#gamma#gamma"),0.08,10,false,220);
  draw_smth(etagg_bkg_ds,e_g1,string("e_g1_bkg_etagg_bar"),g1_bar_cut,string("Photon 1 energy (barrel), bkg #eta#to#gamma#gamma"), 0.08,10,false,500);
  draw_smth(etagg_bkg_ds,e_g2,string("e_g2_bkg_etagg_frw"),g2_frw_cut,string("Photon 2 energy (forward), bkg #eta#to#gamma#gamma"),0.08,10,false,220);
  draw_smth(etagg_bkg_ds,e_g2,string("e_g2_bkg_etagg_bar"),g2_bar_cut,string("Photon 2 energy (barrel), bkg #eta#to#gamma#gamma"), 0.08,10,false,500);

  // pi0 //
//  draw_smth(omega_sig_ds,pt_pip,string("pt_pip_sig_pi0"),string(""),string("p_{t}(#pi^{+}), signal #pi^{0}"),0.,5,false,2200);
//  draw_smth(omega_bkg_ds,pt_pip,string("pt_pip_bkg_pi0"),string(""),string("p_{t}(#pi^{+}), bkg #pi^{0}"),   0.,5,false,2200);

//  draw_smth(omega_sig_ds,pt_pim,string("pt_pim_sig_omega"),string(""),string("p_{t}(#pi^{-}), signal #pi^{0}"),0.,5,false,2200);
//  draw_smth(omega_bkg_ds,pt_pim,string("pt_pim_bkg_omega"),string(""),string("p_{t}(#pi^{-}), bkg #pi^{0}"),   0.,5,false,2200);

//  draw_smth(pi0_sig_ds,z_pip,string("z_pip_sig_pi0"),string(""),string("z(#pi^{+}), signal #pi^{0}"),-5,5,true,2200);
//  draw_smth(pi0_bkg_ds,z_pip,string("z_pip_bkg_pi0"),string(""),string("z(#pi^{+}), bkg #pi^{0}"),   -5,5,true,2200);

//  draw_smth(pi0_sig_ds,r_pip,string("r_pip_sig_pi0"),string(""),string("r(#pi^{+}), signal #pi^{0}"),-2,2,true,2200);
//  draw_smth(pi0_bkg_ds,r_pip,string("r_pip_bkg_pi0"),string(""),string("r(#pi^{+}), bkg #pi^{0}"),   -2,2,true,2200);

//  draw_smth(pi0_sig_ds,z_pim,string("z_pim_sig_pi0"),string(""),string("z(#pi^{-}), signal #pi^{0}"),-5,5,true,2200);
//  draw_smth(pi0_bkg_ds,z_pim,string("z_pim_bkg_pi0"),string(""),string("z(#pi^{-}), bkg #pi^{0}"),   -5,5,true,2200);

//  draw_smth(pi0_sig_ds,r_pim,string("r_pim_sig_pi0"),string(""),string("r(#pi^{-}), signal #pi^{0}"),-2,2,true,2200);
//  draw_smth(pi0_bkg_ds,r_pim,string("r_pim_bkg_pi0"),string(""),string("r(#pi^{-}), bkg #pi^{0}"),   -2,2,true,2200);

//  draw_smth(pi0_sig_ds,hel_h0,string("hel_h0_sig_pi0"),string(""),string("(E_{1}-E_{2})/(E_{1}+E_{2}), signal #pi^{0}"),-1,1,false,2200);
//  draw_smth(pi0_bkg_ds,hel_h0,string("hel_h0_bkg_pi0"),string(""),string("(E_{1}-E_{2})/(E_{1}+E_{2}), bkg #pi^{0}"),   -1,1,false,220);

//  draw_smth(pi0_sig_ds,e_g1,string("e_g1_sig_pi0_frw"),g1_frw_cut,string("Photon 1 energy (forward), signal #pi^{0}"),0.04,10,false,2200);
//  draw_smth(pi0_bkg_ds,e_g1,string("e_g1_bkg_pi0_frw"),g1_frw_cut,string("Photon 1 energy (forward), bkgd #pi^{0}"),  0.04,10,false,220);

//  draw_smth(pi0_sig_ds,e_g1,string("e_g1_sig_pi0_bar"),g1_bar_cut,string("Photon 1 energy (barrel), signal #pi^{0}"),0.04,10,false,5000);
//  draw_smth(pi0_bkg_ds,e_g1,string("e_g1_bkg_pi0_bar"),g1_bar_cut,string("Photon 1 energy (barrel), bkg #pi^{0}"),   0.04,10,false,500);

//  draw_smth(pi0_sig_ds,e_g2,string("e_g2_sig_pi0_frw"),g2_frw_cut,string("Photon 2 energy (forward), signal #pi^{0}"),0.04,10,false,2200);
//  draw_smth(pi0_bkg_ds,e_g2,string("e_g2_bkg_pi0_frw"),g2_frw_cut,string("Photon 2 energy (forward), bkgd #pi^{0}"),  0.04,10,false,220);

//  draw_smth(pi0_sig_ds,e_g2,string("e_g2_sig_pi0_bar"),g2_bar_cut,string("Photon 2 energy (barrel), signal #pi^{0}"),0.04,10,false,5000);
//  draw_smth(pi0_bkg_ds,e_g2,string("e_g2_bkg_pi0_bar"),g2_bar_cut,string("Photon 2 energy (barrel), bkg #pi^{0}"),   0.04,10,false,500);

  // omega //
//  draw_smth(omega_sig_ds,pt_pi1,string("pt_pi1_sig_omega"),string(""),string("p_{t}(#pi^{1}), signal #omega"),0.1,5,false,2200);
//  draw_smth(omega_bkg_ds,pt_pi1,string("pt_pi1_bkg_omega"),string(""),string("p_{t}(#pi^{1}), bkg #omega"),   0.1,5,false,2200);

//  draw_smth(omega_sig_ds,pt_pi2,string("pt_pi2_sig_omega"),string(""),string("p_{t}(#pi^{2}), signal #omega"),0.1,5,false,2200);
//  draw_smth(omega_bkg_ds,pt_pi2,string("pt_pi2_bkg_omega"),string(""),string("p_{t}(#pi^{2}), bkg #omega"),   0.1,5,false,2200);

//  draw_smth(omega_sig_ds,z_pi1,string("z_pi1_sig_omega"),string(""),string("z(#pi^{1}), signal #omega"),-1,1,true,2200);
//  draw_smth(omega_bkg_ds,z_pi1,string("z_pi1_bkg_omega"),string(""),string("z(#pi^{1}), bkg #omega"),   -1,1,true,2200);

//  draw_smth(omega_sig_ds,r_pi1,string("r_pi1_sig_omega"),string(""),string("r(#pi^{1}), signal #omega"),-0.1,0.1,true,2200);
//  draw_smth(omega_bkg_ds,r_pi1,string("r_pi1_bkg_omega"),string(""),string("r(#pi^{1}), bkg #omega"),   -0.1,0.1,true,2200);

//  draw_smth(omega_sig_ds,z_pi2,string("z_pi2_sig_omega"),string(""),string("z(#pi^{2}), signal #omega"),-1,1,true,2200);
//  draw_smth(omega_bkg_ds,z_pi2,string("z_pi2_bkg_omega"),string(""),string("z(#pi^{2}), bkg #omega"),   -1,1,true,2200);

//  draw_smth(omega_sig_ds,r_pi2,string("r_pi2_sig_omega"),string(""),string("r(#pi^{2}), signal #omega"),-0.1,0.1,true,2200);
//  draw_smth(omega_bkg_ds,r_pi2,string("r_pi2_bkg_omega"),string(""),string("r(#pi^{2}), bkg #omega"),   -0.1,0.1,true,2200);

//  draw_smth(omega_sig_ds,hel_pi0,string("hel_pi0_sig_omega"),string(""),string("(E_{1}-E_{2})/(E_{1}+E_{2}), signal #omega"),0.06,10,false,2200);
//  draw_smth(omega_bkg_ds,hel_pi0,string("hel_pi0_bkg_omega"),string(""),string("(E_{1}-E_{2})/(E_{1}+E_{2}), bkg #omega"),0.06,10,false,220);

//  draw_smth(omega_sig_ds,e_g1,string("e_g1_sig_omega_frw"),g1_frw_cut,string("Photon 1 energy (forward), signal #omega"),0.06,10,false,2200);
//  draw_smth(omega_bkg_ds,e_g1,string("e_g1_bkg_omega_frw"),g1_frw_cut,string("Photon 1 energy (forward), bkgd #omega"),0.06,10,false,220);

//  draw_smth(omega_sig_ds,e_g1,string("e_g1_sig_omega_bar"),g1_bar_cut,string("Photon 1 energy (barrel), signal #omega"),0.04,10,false,5000);
//  draw_smth(omega_bkg_ds,e_g1,string("e_g1_bkg_omega_bar"),g1_bar_cut,string("Photon 1 energy (barrel), bkg #omega"),0.04,10,false,500);

//  draw_smth(omega_sig_ds,e_g2,string("e_g2_sig_omega_frw"),g2_frw_cut,string("Photon 2 energy (forward), signal #omega"),0.06,10,false,2200);
//  draw_smth(omega_bkg_ds,e_g2,string("e_g2_bkg_omega_frw"),g2_frw_cut,string("Photon 2 energy (forward), bkgd #omega"),0.06,10,false,220);

//  draw_smth(omega_sig_ds,e_g2,string("e_g2_sig_omega_bar"),g2_bar_cut,string("Photon 2 energy (barrel), signal #omega"),0.04,10,false,5000);
//  draw_smth(omega_bkg_ds,e_g2,string("e_g2_bkg_omega_bar"),g2_bar_cut,string("Photon 2 energy (barrel), bkg #omega"),0.04,10,false,500);
  return;
}
