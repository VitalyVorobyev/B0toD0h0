int RooCPFit(void){
  const double cm2ps = 78.4857;
  ifile = TFile::Open("/home/vitaly/B0toDh0/TMVA/FIL_b2dh_sigOmega_s2_full.root");
  TTree* tree = (TTree*)ifile->Get("TEvent");

  RooRealVar dz_mc("dz_mc","#Deltaz",-10./cm2ps,10./cm2ps,"mm");
//  RooCategory tag("flv_mc","flv_mc");
//  tag->defineType("B0",1);
//  tag->defineType("anti-B0",-1);
//  RooCategory bin("bin_mc","bin_mc");
//  bin->defineType("1",1); bin->defineType("-1",-1);
//  bin->defineType("2",2); bin->defineType("-2",-2);
//  bin->defineType("3",3); bin->defineType("-3",-3);
//  bin->defineType("4",4); bin->defineType("-4",-4);
//  bin->defineType("5",5); bin->defineType("-5",-5);
//  bin->defineType("6",6); bin->defineType("-6",-6);
//  bin->defineType("7",7); bin->defineType("-7",-7);
//  bin->defineType("8",8); bin->defineType("-8",-8);
  RooDataSet d("data","data",tree,RooArgSet(dz_mc),"(dz_mc>0 || dz_mc<0)");
  d.Print();

  RooRealVar mean("mean","mean",0.0,-0.1,0.1,"cm");
  RooRealVar sigma("sigma","sigma",0.02,0.,1.5,"cm");
  RooGaussian pdf1("pdf1","pdf1",dz_mc,mean,sigma);

  pdf1.fitTo(d);

  ifile->Close();
  return 0;
}
