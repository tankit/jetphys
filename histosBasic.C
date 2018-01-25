// Purpose:  Create basic histograms for inclusive jets analysis
//
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 8, 2015
#include "histosBasic.h"
#include "settings.h"

#include "TMath.h"

#include <iostream>

using namespace std;

histosBasic::histosBasic(TDirectory *dir, string trigname, string cotrig,
			 double etamin, double etamax,
			 double pttrg, double ptmin, double ptmax,
			 bool ismcdir)
  : lumsum(0), lumsum2(0) {

  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // phase space
  this->trigname = trigname;
  this->cotrig = cotrig;
  this->etamin = etamin;
  this->etamax = etamax;
  this->pttrg = pttrg;
  this->ptmin = ptmin;
  this->ptmax = ptmax;

  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);

  // Binning agreed within JTF: pT>100 GeV from CaloJet resolutions,
  // pT<100 GeV to optimize bin widths for PFJets and b-tagging
  // (little higher than resolution, but fairly flat relative width)
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/QCDAnalysis/HighPtJetAnalysis/interface/DefaultPtBins.h?revision=1.2&view=markup
  const double x0[] =
    {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
     4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int nx0 = sizeof(x0)/sizeof(x0[0])-1;
  //
  const double *bx0 = &x0[0];
  const int nbx0 = nx0;

  // Wider version of the binning for less statistical scatter for b-jets
  const double xW[] =
    {1, 5, 15, 24, 37, 56, 84, 114, 153, 196, 245, 330, 430, 548, 686, 846,
     1032, 1248, 1497, 1784, 2116, 2500, 2941, 3450, 3637,
     4252, 4961, 5777, 6717, 7000};
  const int nxW = sizeof(xW)/sizeof(xW[0])-1;
  //
  const double *bxW = &xW[0];
  const int nbxW = nxW;


// Optimized binning created by optimizeBins.C ("MC"; lumi 1000/pb, eff 1e+10%)
// Using NLOxNP theory fit as input when available
  const int neta = 8;
  const int nbins = 65;
double vx[neta][nbins] =
  {{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492, 0}, // Eta_0.5-1.0
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832, 4037, 0, 0, 0, 0, 0}, // Eta_1.0-1.5
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_1.5-2.0
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_2.0-2.5
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_2.5-3.0
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_3.0-3.5
   {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; // Eta_3.5-4.0

  const double etarange[] =
  {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  const unsigned int netas = sizeof(etarange)/sizeof(etarange[0])-1;

  int ieta = int(0.5*(etamin+etamax)/0.5); assert(ieta<neta);
  vector<double> x;
  for (int i = 0; i != nbins && vx[ieta][i]!=0; ++i) {
    x.push_back(vx[ieta][i]);
  } // for i
  const int nx = x.size()-1;

  const double ay[] =
    {0, 0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411,
     2.5, 2.853, 2.964, 5.191};
  const int nay = sizeof(ay)/sizeof(ay[0]);

  vector<double> yW(nay);
  for (unsigned int i = 0; i != yW.size(); ++i) {
    yW[i] = ay[i];
  }
  const int nyW = yW.size()-1;

  vector<double> y(51);
  for (unsigned int i = 0; i != y.size(); ++i) {
    y[i] = -5. + 0.2*i;
  }
  const int ny = y.size()-1;

  vector<double> pv(26);
  for (unsigned int i = 0; i != pv.size(); ++i) {
    pv[i] = -0.5 + i;
  }
  const int npv = pv.size()-1;

  // new histograms for quark/gluon study (Ozlem)
  hgpt = new TH1D("hgpt","",nx,&x[0]);
  hgpt0 = new TH1D("hgpt0","",6500,0.,6500.);

  hqgl = new TH1D("hqgl","",101,0,1.01);
  hqgl_q = new TH1D("hqgl_q","",101,0,1.01);
  hqgl_g = new TH1D("hqgl_g","",101,0,1.01);
  hqgl_dq = new TH1D("hqgl_dq","",101,0,1.01);
  hqgl_dg = new TH1D("hqgl_dg","",101,0,1.01);

  // add TH2Ds for pT bins here
  hqgl2 =  new TH2D("hqgl2","",nx,&x[0],101,0,1.01);
  hqgl2_g =  new TH2D("hqgl2_g","",nx,&x[0],101,0,1.01);
  hqgl2_q =  new TH2D("hqgl2_q","",nx,&x[0],101,0,1.01);
  hqgl2_dg =  new TH2D("hqgl2_dg","",nx,&x[0],101,0,1.01);
  hqgl2_dq =  new TH2D("hqgl2_dq","",nx,&x[0],101,0,1.01);

  // raw spectrum
  hpt = new TH1D("hpt","",nx,&x[0]);
  hpt_pre = new TH1D("hpt_pre","",nx,&x[0]); // prescale weighed
  hpt_gtw = new TH1D("hpt_gtw","",nx,&x[0]); // _mc per trigger
  hpt_g0tw = new TH1D("hpt_g0tw","",nx,&x[0]); // _mc per trigger
  hgpt_g0tw = new TH1D("hgpt_g0tw","",nx,&x[0]);//Ozlem
  hqpt_g0tw = new TH1D("hqpt_g0tw","",nx,&x[0]);//Ozlem
  hgpt_g0 = new TH1D("hgpt_g0","",nx,&x[0]);//Ozlem
  hqpt_g0 = new TH1D("hqpt_g0","",nx,&x[0]);//Ozlem
  hcopt = new TH1D("hcopt","",nx,&x[0]);


  hpt_noid = new TH1D("hpt_noid","",nx,&x[0]);
  hpt_nojetid = new TH1D("hpt_nojetid","",nx,&x[0]);
  hpt_noevtid = new TH1D("hpt_noevtid","",nx,&x[0]);
  hptevt = new TH1D("hptevt","",nx,&x[0]);
  hpttmp = new TH1D("hpttmp","",nx,&x[0]);
  //hpt->Sumw2();
  //hcopt->Sumw2();
  //hpt_noid->Sumw2();
  //hpt_nojetid->Sumw2();
  //hpt_noevtid->Sumw2();
  //hptevt->Sumw2();

  // delete-m jackknife
  hpt_jk.resize(10);
  for (unsigned int i = 0; i != hpt_jk.size(); ++i) {
    hpt_jk[i] = new TH1D(Form("hpt_jk%d",i+1),"",nx,&x[0]);
  }
  h2jk = new TH2D("h2jk","Check of reshuffling",10,-0.5,9.5,10,-0.5,9.5);

  hpt_tmp = new TH1D("hpt_tmp","",nx,&x[0]);
  hpt_evtcount = new TH1D("hpt_evtcount","",nx,&x[0]);
  hpt_evt = new TH1D("hpt_evt","",nx,&x[0]);
  hpt_jet = new TH1D("hpt_jet","",nx,&x[0]);
  //
  const double nj = 3;
  vector<double> vnj(nj+1); vnj[0]=0; vnj[1]=1; vnj[2]=2; vnj[3]=3;

  // 1 GeV bins for localizing leading jets
  //hpt0 = new TH1D("hpt0","",int(_jp_emax),0.,_jp_emax);
  hpt0 = new TH1D("hpt0","",6500,0.,6500.);
  //hpt0->Sumw2(); 
  
  // leading and non-leading jets
  hpt1 = new TH1D("hpt1","",nx,&x[0]);
  hpt2 = new TH1D("hpt2","",nx,&x[0]);
  hpt3 = new TH1D("hpt3","",nx,&x[0]);

  // dijet mass
  hdjmass = new TH1D("hdjmass","",nx,&x[0]);
  hdjmass0 = new TH1D("hdjmass0","",int(_jp_sqrts),0.,_jp_sqrts);
  //hdjmass0 = new TH1D("hdjmass0","",13000,0.,13000.);
  pdjmass_ptratio = new TProfile("pdjmass_ptratio","",nx,&x[0]);
  pdjmass0_ptratio = new TProfile("pdjmass0_ptratio","",
				  int(_jp_sqrts),0.,_jp_sqrts);
  //pdjmass0_ptratio = new TProfile("pdjmass0_ptratio","",13000,0.,13000.);
  //hdjmass->Sumw2();
  //hdjmass0->Sumw2();
  //pdjmass_ptratio->Sumw2();
  //pdjmass0_ptratio->Sumw2();

  // basic properties
  ppt = new TProfile("ppt","",nx,&x[0]);
  pmass = new TProfile("pmass","",nx0,&x0[0]);
  hmass = new TH1D("hmass","",100,0.,0.5);
  pjec = new TProfile("pjec","",nx,&x[0]);
  pjec2 = new TProfile("pjec2","",nx,&x[0]);
  punc = new TProfile("punc","",nx,&x[0]);
  hnpvgood = new TH1D("hnpvgood","",100,-0.5,99.5);
  hrho = new TH1D("hrho","",200,0,100);
  //pmass->Sumw2();
  //hmass->Sumw2();
  //hnpvgood->Sumw2();
  //hrho->Sumw2();

  // JEC monitoring
  pjec_l1 = new TProfile("pjec_l1","",nx,&x[0]);
  pjec_l2l3 = new TProfile("pjec_l2l3","",nx,&x[0]);
  pjec_res = new TProfile("pjec_res","",nx,&x[0]);

  // pile-up information
  pa = new TProfile("pa","",nx,&x[0]);
  //
  ptrpu = new TProfile("ptrpu","",nx,&x[0]);
  pnpv = new TProfile("pnpv","",nx,&x[0]);
  pnpvall = new TProfile("pnpvall","",nx,&x[0]);
  prho = new TProfile("prho","",nx,&x[0]);
  prhovsnpv = new TProfile("prhovsnpv","",50,-0.5,49.5);
  pnpvvsrho = new TProfile("pnpvvsrho","",50,-0.5,49.5);
  prhovsnpvall = new TProfile("prhovsnpvall","",50,-0.5,49.5);
  h2rhovsnpv = new TH2D("h2rhovsnpv","",50,0.5,50.5,200,0,40);
  //
  prhovstrpu = new TProfile("prhovstrpu","",50,-0.5,49.5);
  pnpvvstrpu = new TProfile("pnpvvstrpu","",50,-0.5,49.5);
  pnpvallvstrpu = new TProfile("pnpvallvstrpu","",50,-0.5,49.5);
  pitpuvstrpu = new TProfile("itpuvstrpu","",50,-0.5,49.5);
  htrpu2 = new TH1D("htrpu2","",50,-0.5,49.5);
  hjet_vstrpu = new TH1D("hjet_vstrpu","",50,-0.5,49.5);
  hlumi_vstrpu = new TH1D("hlumi_vstrpu","",50,-0.5,49.5);
  //
  // luminosity
  hlumi = new TH1D("hlumi","",nx,&x[0]);
  hlumi2 = new TH1D("hlumi2","",nx,&x[0]);

  // inclusive efficiencies
  peff = new TProfile("peff","",nx,&x[0]);
  pideff = new TProfile("pideff","",nx,&x[0]);
  pvtxeff = new TProfile("pvtxeff","",nx,&x[0]);
  pdqmeff = new TProfile("pdqmeff","",nx,&x[0]);

  // control plots of components (JEC)
  pncand = new TProfile("pncand","",nx0,&x0[0]);
  pnch = new TProfile("pnch","",nx0,&x0[0]);
  pnne = new TProfile("pnne","",nx0,&x0[0]);
  pnnh = new TProfile("pnnh","",nx0,&x0[0]);
  pnce = new TProfile("pnce","",nx0,&x0[0]);
  pnmu = new TProfile("pnmu","",nx0,&x0[0]);
  pchf = new TProfile("pchf","",nx0,&x0[0]);
  pnef = new TProfile("pnef","",nx0,&x0[0]);
  pnhf = new TProfile("pnhf","",nx0,&x0[0]);
  pcef = new TProfile("pcef","",nx0,&x0[0]);
  pmuf = new TProfile("pmuf","",nx0,&x0[0]);
  pbeta = new TProfile("pbeta","",nx0,&x0[0]);
  pbetastar = new TProfile("pbetastar","",nx0,&x0[0]);
  hncand = new TH1D("hncand","",300,-0.5,299.5);
  hnch = new TH1D("hnch","",300,-0.5,299.5);
  hnne = new TH1D("hnne","",300,-0.5,299.5);
  hnnh = new TH1D("hnnh","",300,-0.5,299.5);
  hnce = new TH1D("hnce","",300,-0.5,299.5);
  hnmu = new TH1D("hnmu","",300,-0.5,299.5);
  hchf = new TH1D("hchf","",110,0.,1.1);
  hnef = new TH1D("hnef","",110,0.,1.1);
  hnhf = new TH1D("hnhf","",110,0.,1.1);
  hcef = new TH1D("hcef","",110,0.,1.1);
  hmuf = new TH1D("hmuf","",110,0.,1.1);
  hbeta = new TH1D("hbeta","",110,0.,1.1);
  hbetastar = new TH1D("hbetastar","",110,0.,1.1);
  // control plots of components (JEC tag-and-probe)
  pncandtp = new TProfile("pncandtp","",nx0,&x0[0]);
  pnchtp = new TProfile("pnchtp","",nx0,&x0[0]);
  pnnetp = new TProfile("pnnetp","",nx0,&x0[0]);
  pnnhtp = new TProfile("pnnhtp","",nx0,&x0[0]);
  pncetp = new TProfile("pncetp","",nx0,&x0[0]);
  pnmutp = new TProfile("pnmutp","",nx0,&x0[0]);
  pchftp = new TProfile("pchftp","",nx0,&x0[0]);
  pneftp = new TProfile("pneftp","",nx0,&x0[0]);
  pnhftp = new TProfile("pnhftp","",nx0,&x0[0]);
  pceftp = new TProfile("pceftp","",nx0,&x0[0]);
  pmuftp = new TProfile("pmuftp","",nx0,&x0[0]);
  pbetatp = new TProfile("pbetatp","",nx0,&x0[0]);
  pbetastartp = new TProfile("pbetastartp","",nx0,&x0[0]);
  //
  pchftp2 = new TProfile("pchftp2","",nx0,&x0[0]);
  pneftp2 = new TProfile("pneftp2","",nx0,&x0[0]);
  pnhftp2 = new TProfile("pnhftp2","",nx0,&x0[0]);
  pceftp2 = new TProfile("pceftp2","",nx0,&x0[0]);
  pmuftp2 = new TProfile("pmuftp2","",nx0,&x0[0]);
  //
  hncandtp = new TH1D("hncandtp","",300,-0.5,299.5);
  hnchtp = new TH1D("hnchtp","",300,-0.5,299.5);
  hnnetp = new TH1D("hnnetp","",300,-0.5,299.5);
  hnnhtp = new TH1D("hnnhtp","",300,-0.5,299.5);
  hncetp = new TH1D("hncetp","",300,-0.5,299.5);
  hnmutp = new TH1D("hnmutp","",300,-0.5,299.5);
  hchftp = new TH1D("hchftp","",110,0.,1.1);
  hneftp = new TH1D("hneftp","",110,0.,1.1);
  hnhftp = new TH1D("hnhftp","",110,0.,1.1);
  hceftp = new TH1D("hceftp","",110,0.,1.1);
  hmuftp = new TH1D("hmuftp","",110,0.,1.1);
  hbetatp = new TH1D("hbetatp","",110,0.,1.1);
  hbetastartp = new TH1D("hbetastartp","",110,0.,1.1);
  // control plots vs NPV
  pncandtp_vsnpv = new TProfile("pncandtp_vsnpv","",50,-0.5,49.5);
  pnchtp_vsnpv = new TProfile("pnchtp_vsnpv","",50,-0.5,49.5);
  pnnetp_vsnpv = new TProfile("pnnetp_vsnpv","",50,-0.5,49.5);
  pnnhtp_vsnpv = new TProfile("pnnhtp_vsnpv","",50,-0.5,49.5);
  pncetp_vsnpv = new TProfile("pncetp_vsnpv","",50,-0.5,49.5);
  pnmutp_vsnpv = new TProfile("pnmutp_vsnpv","",50,-0.5,49.5);
  pchftp_vsnpv = new TProfile("pchftp_vsnpv","",50,-0.5,49.5);
  pneftp_vsnpv = new TProfile("pneftp_vsnpv","",50,-0.5,49.5);
  pnhftp_vsnpv = new TProfile("pnhftp_vsnpv","",50,-0.5,49.5);
  pceftp_vsnpv = new TProfile("pceftp_vsnpv","",50,-0.5,49.5);
  pmuftp_vsnpv = new TProfile("pmuftp_vsnpv","",50,-0.5,49.5);
  pbetatp_vsnpv = new TProfile("pbetatp_vsnpv","",50,-0.5,49.5);
  pbetastartp_vsnpv = new TProfile("pbetastartp_vsnpv","",50,-0.5,49.5);
  //
  pchftp_vstrpu = new TProfile("pchftp_vstrpu","",50,-0.5,49.5);
  pneftp_vstrpu = new TProfile("pneftp_vstrpu","",50,-0.5,49.5);
  pnhftp_vstrpu = new TProfile("pnhftp_vstrpu","",50,-0.5,49.5);
  pceftp_vstrpu = new TProfile("pceftp_vstrpu","",50,-0.5,49.5);
  pmuftp_vstrpu = new TProfile("pmuftp_vstrpu","",50,-0.5,49.5);
  pbetatp_vstrpu = new TProfile("pbetatp_vstrpu","",50,-0.5,49.5);
  pbetastartp_vstrpu = new TProfile("pbetastartp_vstrpu","",50,-0.5,49.5);
  this->ismcdir = ismcdir;

  //pncand->Sumw2();
  //pnch->Sumw2();
  //pnne->Sumw2();
  //pnnh->Sumw2();
  //pchf->Sumw2();
  //pnef->Sumw2();
  //pnhf->Sumw2();
  //hncand->Sumw2();
  //hnch->Sumw2();
  //hnne->Sumw2();
  //hnnh->Sumw2();
  //hchf->Sumw2();
  //hnef->Sumw2();
  //hnhf->Sumw2();
  ////
  //pncandtp->Sumw2();
  //pnchtp->Sumw2();
  //pnnetp->Sumw2();
  //pnnhtp->Sumw2();
  //pchftp->Sumw2();
  //pneftp->Sumw2();
  //pnhftp->Sumw2();
  //hncandtp->Sumw2();
  //hnchtp->Sumw2();
  //hnnetp->Sumw2();
  //hnnhtp->Sumw2();
  //hchftp->Sumw2();
  //hneftp->Sumw2();
  //hnhftp->Sumw2();

  // control plots for topology (JEC)
  hselpt = new TH1D("hselpt","",nx,&x[0]);
  hy = new TH1D("hy","",100,-5.,5.); // May 11
  hy2 = new TH1D("hy2","", 100,-4.799,4.799);
  heta = new TH1D("heta","",100,-5.,5.);
  heta2 = new TH1D("heta2","",100,-4.799,4.799);
  hphi = new TH1D("hphi","",144,-TMath::TwoPi(),TMath::TwoPi());
  hdphi = new TH1D("hdphi","",36,0.,TMath::Pi());
  hdpt = new TH1D("hdpt","",100,0.,1.);
  pdpt = new TProfile("pdpt","",nx,&x[0]);
  hjet = new TH1D("hjet","",100,0.,1.);
  hmet = new TH1D("hmet","",100,0.,1.);
  hmetphi = new TH1D("hmetphi","",36,0.,TMath::Pi());
  // control plots for vertex
  hpvndof = new TH1D("hpvndof","",400,0.,400.);
  hpvx = new TH1D("hpvx","",400,-0.2,0.2);//-2.,2.);
  hpvy = new TH1D("hpvy","",400,-0.2,0.2);//-2.,2.);
  hpvz = new TH1D("hpvz","",400,-50.,50.);
  hpvr = new TH1D("hpvr","",200,0.,0.2);//2.);
  hpvrho = new TH1D("hpvrho","",200,0.,0.2);//2.);
  // closure plots for JEC
  hmpf = new TH1D("hmpf","",200,0.,2.);
  hmpf1 = new TH1D("hmpf1","",200,0.,2.);
  hmpf2 = new TH1D("hmpf2","",200,0.,2.);
  hmpfx = new TH1D("hmpfx","",200,0.,2.);
  hmpfy = new TH1D("hmpfy","",200,0.,2.);
  hmpfz = new TH1D("hmpfz","",200,0.,2.);
  pmpf = new TProfile("pmpf","",nx,&x[0]);
  pmpf1 = new TProfile("pmpf1","",nx,&x[0]);
  pmpf2 = new TProfile("pmpf2","",nx,&x[0]);
  pmpfx = new TProfile("pmpfx","",nx,&x[0]);
  pmpfy = new TProfile("pmpfy","",nx,&x[0]);
  pmpfz = new TProfile("pmpfz","",nx,&x[0]);

  const int n3 = 40;
  vector<double> v3(n3+1);
  for (unsigned int i = 0; i != n3+1; ++i) v3[i] = 0. + 1.*i/n3;
  const int na = 200;
  vector<double> va(na+1);
  for (unsigned int i = 0; i != na+1; ++i) va[i] = -1. + 2.*i/na;
  hdjasymm = new TH3D("hdjasymm",";p_{T,ave};p_{T,3rd}/p_{T,ave};Asymmetry",
		      nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjmpf = new TH3D("hdjmpf",";p_{T,ave};p_{T,3rd}/p_{T,ave};MPF",
		    nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjasymmtp = new TH3D("hdjasymmtp",";p_{T,tag};p_{T,3rd}/p_{T,tag};Asymmetry",
		      nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjmpftp = new TH3D("hdjmpftp",";p_{T,tag};p_{T,3rd}/p_{T,tag};MPF",
		      nx,&x[0],n3,&v3[0],na,&va[0]);

  hr21 = new TH1D("hr21",";pt2/pt1",50,0.,1.);
  hr31 = new TH1D("hr31",";pt3/pt1",50,0.,1.);
  hr32 = new TH1D("hr32",";pt3/pt2",50,0.,1.);
  pr21 = new TProfile("pr21",";npvgood;pt2/pt1",50,-0.5,49.5);
  pr31 = new TProfile("pr31",";npvgood;pt3/pt1",50,-0.5,49.5);
  pr32 = new TProfile("pr32",";npvgood;pt3/pt2",50,-0.5,49.5);
  px21 = new TProfile("px21",";npvgood;has pt2/pt1",50,-0.5,49.5);
  px31 = new TProfile("px31",";npvgood;has pt3/pt1",50,-0.5,49.5);
  px32 = new TProfile("px32",";npvgood;has pt3/pt2",50,-0.5,49.5);

  hyeta = new TH1D("hyeta","",100,-0.436/4.,0.436/4.);
  hyeta2 = new TH1D("hyeta2","",100,-0.436/4.,0.436/4.);
  hbetabetastar = new TH2D("hbetabetastar","",110,0.,1.1,110,0.,1.1);
  hetaphi = new TH2D("hetaphi","",100,-4.799,4.799,
		     144,-TMath::TwoPi(),TMath::TwoPi());
  //hyeta->Sumw2();
  //hyeta2->Sumw2();

  pdpt->Sumw2();
  pmpf->Sumw2();
  pmpf1->Sumw2();
  pmpf2->Sumw2();
  pmpfx->Sumw2();

  hdjasymm->Sumw2();
  hdjmpf->Sumw2();
  hdjasymmtp->Sumw2();
  hdjmpftp->Sumw2(); 


  // MC checks
  //htrpu = new TH1D("htrpu","",100,-0.5,99.5);
  htrpu = new TH1D("htrpu","",120,0.,60.); // for PU reweighing
  if (this->ismcdir) {
    //hpt_jt30 = new TH1D("hpt_jt30","",nx,&x[0]);
    //hpt_jt60 = new TH1D("hpt_jt60","",nx,&x[0]);
    //hpt_jt110 = new TH1D("hpt_jt110","",nx,&x[0]);
    //hpt_jt190 = new TH1D("hpt_jt190","",nx,&x[0]);
    //hpt_jt240 = new TH1D("hpt_jt240","",nx,&x[0]);
    //hpt_jt300 = new TH1D("hpt_jt300","",nx,&x[0]);
    //hpt_jt370 = new TH1D("hpt_jt370","",nx,&x[0]);
    ////
    //hpt0_jt30 = new TH1D("hpt0_jt30","",3450,0.,3450.);
    //hpt0_jt60 = new TH1D("hpt0_jt60","",3450,0.,3450.);
    //hpt0_jt110 = new TH1D("hpt0_jt110","",3450,0.,3450.);
    //hpt0_jt190 = new TH1D("hpt0_jt190","",3450,0.,3450.);
    //hpt0_jt240 = new TH1D("hpt0_jt240","",3450,0.,3450.);
    //hpt0_jt300 = new TH1D("hpt0_jt300","",3450,0.,3450.);
    //hpt0_jt370 = new TH1D("hpt0_jt370","",3450,0.,3450.);


    hpthat = new TH1D("hpthat","",nx,&x[0]);
    hpthatnlo = new TH1D("hpthatnlo","",nx,&x[0]);
    
    //unfolding studies (Mikael)
    //mT: (pTgen,ygen); (pTreco,yreco)
    mT = new TH2D("mT","mT(yjet);p_{T,gen};p_{T,reco}",nx,&x[0],nx,&x[0]);
    mTuw = new TH2D("mTuw","mTuw(yjet);p_{T,gen};p_{T,reco}",nx,&x[0],nx,&x[0]);
    mTf = new TH2D("mTf","mT(yjet);p_{T,gen};p_{T,reco}",
		   3485,15,3500, 3485,15,3500);
    mTfuw = new TH2D("mTfuw","mT(yjet);p_{T,gen};p_{T,reco}",
		   3485,15,3500, 3485,15,3500);
    mx = new TH1D("mx","mx(ygen);p_{T,gen}",nx,&x[0]); // pTgen, ygen
    mxuw = new TH1D("mxuw","mx(ygen);p_{T,gen}",nx,&x[0]); // pTgen, ygen
    mxf = new TH1D("mxf","mx(ygen);p_{T,gen}",3485,15,3500); // pTgen, ygen
    mxfuw = new TH1D("mxfuw","mx(ygen);p_{T,gen}",3485,15,3500); // pTgen, ygen
    my = new TH1D("my","my(yreco);p_{T,reco}",nx,&x[0]); // pTreco, yreco
    myuw = new TH1D("myuw","my(yreco);p_{T,reco}",nx,&x[0]); // pTreco, yreco
    myf = new TH1D("myf","my(yreco);p_{T,reco}",3485,15,3500); // pTreco, yreco
    myfuw = new TH1D("myfuw","my(yreco);p_{T,reco}",3485,15,3500); // pTreco, yreco

    //htrpu = new TH1D("htrpu","",100,-0.5,99.5);
    hitpu = new TH1D("hitpu","",100,-0.5,99.5);
    hootpuearly = new TH1D("hootpuearly","",100,-0.5,99.5);
    hootpulate = new TH1D("hootpulate","",100,-0.5,99.5);
    h2itvsoot = new TH2D("h2itvsoot","",25,-0.5,24.5,50,-0.5,49.5);
    //hitpu->Sumw2();
    //hootpuearly->Sumw2();
    //hootpulate->Sumw2();
    //h2itvsoot->Sumw2();

    hpt_noid_g = new TH1D("hpt_noid_g","",nx,&x[0]);
    hpt_nojetid_g = new TH1D("hpt_nojetid_g","",nx,&x[0]);
    hpt_noevtid_g = new TH1D("hpt_noevtid_g","",nx,&x[0]);

    hpt_r = new TH1D("hpt_r","",nx,&x[0]);
    hpt_g = new TH1D("hpt_g","",nx,&x[0]);
    hpt_gg = new TH1D("hpt_gg","",nx,&x[0]);
    hpt_gg0 = new TH1D("hpt_gg0","",nx,&x[0]);
    hpt_g0 = new TH1D("hpt_g0","",nx,&x[0]);
    hgpt_g = new TH1D("hgpt_g","",nx,&x[0]);//Ozlem
    hqpt_g = new TH1D("hqpt_g","",nx,&x[0]);//Ozlem
    hpt_g0_tmp = new TH1D("hpt_g0_tmp","",nx,&x[0]);
    ppt_r = new TProfile("ppt_r","",nx,&x[0]);
    ppt_g = new TProfile("ppt_g","",nx,&x[0]);
    const double nj = 3;
    vector<double> vnj(nj+1); vnj[0]=0; vnj[1]=1; vnj[2]=2; vnj[3]=3;

    // Response closure
    p3rvsnpv = new TProfile3D("p3rvsnpv","",nx,&x[0],ny,&y[0],npv,&pv[0]);
    p3rvsnpvW = new TProfile3D("p3rvsnpvW","",nxW,&xW[0],nyW,&yW[0],npv,&pv[0]);
    p2rvsnpv = new TProfile2D("p2rvsnpv","",nx,&x[0],50,-0.5,49.5);
    h2r_r = new TH2D("h2r_r","",nx,&x[0],600,0,3);
    h2r_g = new TH2D("h2r_g","",nx,&x[0],600,0,3);
    p2r_r = new TProfile("p2r_r","",nx,&x[0]);
    p2r_g = new TProfile("p2r_g","",nx,&x[0]);
    p2r_ruw = new TProfile("p2r_ruw","",nx,&x[0]);
    p2r_guw = new TProfile("p2r_guw","",nx,&x[0]);

    // Rapidity closure
    h2dy_r = new TH2D("h2dy_r","",nx,&x[0],200,-0.5,0.5);
    h2dy_g = new TH2D("h2dy_g","",nx,&x[0],200,-0.5,0.5);
    p2dy_r = new TProfile("p2dy_r","",nx,&x[0]);
    p2dy_g = new TProfile("p2dy_g","",nx,&x[0]);
    p2dy_ruw = new TProfile("p2dy_ruw","",nx,&x[0]);
    p2dy_guw = new TProfile("p2dy_guw","",nx,&x[0]);
    pdy_r = new TProfile2D("pdy_r","",nx,&x[0],144,0.,TMath::Pi());
    pdy_g = new TProfile2D("pdy_g","",nx,&x[0],144,0.,TMath::Pi());
    //hpthat->Sumw2();
    
    //hpt_noid_g->Sumw2();
    //hpt_nojetid_g->Sumw2();
    //hpt_noevtid_g->Sumw2();
    
    //hpt_r->Sumw2();
    //hpt_g->Sumw2();
    //hpt_gg->Sumw2();
    //hpt_g0->Sumw2();
    ////
    //p3rvsnpv->Sumw2();
    //p2rvsnpv->Sumw2();
    //h2r_r->Sumw2();
    //h2r_g->Sumw2();
    //p2r_r->Sumw2();
    //p2r_g->Sumw2();
    //// was missing from dec3 first attempt
    //h2dy_r->Sumw2();
    //h2dy_g->Sumw2();
    //p2dy_r->Sumw2();
    //p2dy_g->Sumw2();
    //pdy_r->Sumw2();
    //pdy_g->Sumw2();

  } // ismcdir

  curdir->cd();
}

void histosBasic::Write() {
  dir->cd();
  dir->Write();
}

histosBasic::~histosBasic() {
  Write();
};
