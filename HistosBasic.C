// Purpose:  Create basic histograms for inclusive jets analysis
//
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 8, 2015
#include "HistosBasic.h"
#include "settings.h"

#include "TMath.h"

#include <iostream>

using namespace std;

HistosBasic::HistosBasic(TDirectory *dir, string trigname, double etamin, double etamax,
                         double pttrg, double ptmin, double ptmax, bool ismcdir) {
  lumsum = 0;
  lumsum2 = 0;
  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // phase space
  this->trigname = trigname;
  this->etamin = etamin;
  this->etamax = etamax;
  this->pttrg = pttrg;
  this->ptmin = ptmin;
  this->ptmax = ptmax;
  this->ismcdir = ismcdir;

  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);

  unsigned int ieta = int(0.5*(etamin+etamax)/0.5); assert(ieta<_jp_noptranges);
  vector<double> x;
  // Chose ptrange according to eta
  for (unsigned int i = 0; i < _jp_nopts_eta and _jp_ptrangevseta[ieta][i]!=0; ++i) {
    if (_jp_ptrangevseta[ieta][i]<0.001) break; // There are zeros in the end of the array when we're out of barrel
    x.push_back(_jp_ptrangevseta[ieta][i]);
  }
  const int nx = x.size()-1;

  vector<double> y(51);
  for (unsigned int i = 0; i != y.size(); ++i) y[i] = -5. + 0.2*i;
  const int ny = y.size()-1;

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

  // delete-m jackknife
  hpt_jk.resize(10);
  for (unsigned int i = 0; i != hpt_jk.size(); ++i)
    hpt_jk[i] = new TH1D(Form("hpt_jk%d",i+1),"",nx,&x[0]);
  h2jk = new TH2D("h2jk","Check of reshuffling",10,-0.5,9.5,10,-0.5,9.5);

  hpt_tmp = new TH1D("hpt_tmp","",nx,&x[0]);
  hpt_evtcount = new TH1D("hpt_evtcount","",nx,&x[0]);
  hpt_evt = new TH1D("hpt_evt","",nx,&x[0]);
  hpt_jet = new TH1D("hpt_jet","",nx,&x[0]);

  // 1 GeV bins for localizing leading jets
  hpt0 = new TH1D("hpt0","",6500,0.,6500.);

  // leading and non-leading jets
  hpt1 = new TH1D("hpt1","",nx,&x[0]);
  hpt2 = new TH1D("hpt2","",nx,&x[0]);
  hpt3 = new TH1D("hpt3","",nx,&x[0]);

  // dijet mass
  hdjmass = new TH1D("hdjmass","",nx,&x[0]);
  hdjmass0 = new TH1D("hdjmass0","",static_cast<int>(_jp_sqrts),0.,_jp_sqrts);
  pdjmass_ptratio = new TProfile("pdjmass_ptratio","",nx,&x[0]);
  pdjmass0_ptratio = new TProfile("pdjmass0_ptratio","", static_cast<int>(_jp_sqrts),0.,_jp_sqrts);

  // basic properties
  ppt = new TProfile("ppt","",nx,&x[0]);
  pmass = new TProfile("pmass","",_jp_npts,&_jp_ptrange[0]);
  hmass = new TH1D("hmass","",100,0.,0.5);
  pjec = new TProfile("pjec","",nx,&x[0]);
  pjec2 = new TProfile("pjec2","",nx,&x[0]);
  punc = new TProfile("punc","",nx,&x[0]);
  hnpvgood = new TH1D("hnpvgood","",_jp_npvs,_jp_pvrange);
  hrho = new TH1D("hrho","",200,0,100);

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
  prhovsnpv = new TProfile("prhovsnpv","",_jp_npvs,_jp_pvrange);
  pnpvvsrho = new TProfile("pnpvvsrho","",_jp_npvs,_jp_pvrange);
  prhovsnpvall = new TProfile("prhovsnpvall","",_jp_npvs,_jp_pvrange);
  h2rhovsnpv = new TH2D("h2rhovsnpv","",_jp_npvs,_jp_pvrange,200,0,40);
  //
  prhovstrpu = new TProfile("prhovstrpu","",_jp_npvs,_jp_pvrange);
  pnpvvstrpu = new TProfile("pnpvvstrpu","",_jp_npvs,_jp_pvrange);
  pnpvallvstrpu = new TProfile("pnpvallvstrpu","",_jp_npvs,_jp_pvrange);
  pitpuvstrpu = new TProfile("itpuvstrpu","",_jp_npvs,_jp_pvrange);
  htrpu2 = new TH1D("htrpu2","",_jp_npvs,_jp_pvrange);
  hjet_vstrpu = new TH1D("hjet_vstrpu","",_jp_npvs,_jp_pvrange);
  hlumi_vstrpu = new TH1D("hlumi_vstrpu","",_jp_npvs,_jp_pvrange);

  // luminosity
  hlumi = new TH1D("hlumi","",nx,&x[0]);
  hlumi2 = new TH1D("hlumi2","",nx,&x[0]);

  // inclusive efficiencies
  peff = new TProfile("peff","",nx,&x[0]);
  pideff = new TProfile("pideff","",nx,&x[0]);
  pvtxeff = new TProfile("pvtxeff","",nx,&x[0]);
  pdqmeff = new TProfile("pdqmeff","",nx,&x[0]);

  // control plots of components (JEC)
  pncand = new TProfile("pncand","",_jp_npts,&_jp_ptrange[0]);
  pnch = new TProfile("pnch","",_jp_npts,&_jp_ptrange[0]);
  pnne = new TProfile("pnne","",_jp_npts,&_jp_ptrange[0]);
  pnnh = new TProfile("pnnh","",_jp_npts,&_jp_ptrange[0]);
  pnce = new TProfile("pnce","",_jp_npts,&_jp_ptrange[0]);
  pnmu = new TProfile("pnmu","",_jp_npts,&_jp_ptrange[0]);
  pchf = new TProfile("pchf","",_jp_npts,&_jp_ptrange[0]);
  pnef = new TProfile("pnef","",_jp_npts,&_jp_ptrange[0]);
  pnhf = new TProfile("pnhf","",_jp_npts,&_jp_ptrange[0]);
  pcef = new TProfile("pcef","",_jp_npts,&_jp_ptrange[0]);
  pmuf = new TProfile("pmuf","",_jp_npts,&_jp_ptrange[0]);
  pbeta = new TProfile("pbeta","",_jp_npts,&_jp_ptrange[0]);
  pbetastar = new TProfile("pbetastar","",_jp_npts,&_jp_ptrange[0]);
  pbetaprime = new TProfile("pbetaprime","",_jp_npts,&_jp_ptrange[0]);
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
  hbetaprime = new TH1D("hbetaprime","",110,0.,1.1);
  // control plots of components (JEC tag-and-probe)
  pncandtp = new TProfile("pncandtp","",_jp_npts,&_jp_ptrange[0]);
  pnchtp = new TProfile("pnchtp","",_jp_npts,&_jp_ptrange[0]);
  pnnetp = new TProfile("pnnetp","",_jp_npts,&_jp_ptrange[0]);
  pnnhtp = new TProfile("pnnhtp","",_jp_npts,&_jp_ptrange[0]);
  pncetp = new TProfile("pncetp","",_jp_npts,&_jp_ptrange[0]);
  pnmutp = new TProfile("pnmutp","",_jp_npts,&_jp_ptrange[0]);
  pchftp = new TProfile("pchftp","",_jp_npts,&_jp_ptrange[0]);
  pneftp = new TProfile("pneftp","",_jp_npts,&_jp_ptrange[0]);
  pnhftp = new TProfile("pnhftp","",_jp_npts,&_jp_ptrange[0]);
  pceftp = new TProfile("pceftp","",_jp_npts,&_jp_ptrange[0]);
  pmuftp = new TProfile("pmuftp","",_jp_npts,&_jp_ptrange[0]);
  pbetatp = new TProfile("pbetatp","",_jp_npts,&_jp_ptrange[0]);
  pbetastartp = new TProfile("pbetastartp","",_jp_npts,&_jp_ptrange[0]);
  pbetaprimetp = new TProfile("pbetaprimetp","",_jp_npts,&_jp_ptrange[0]);
  //
  pchftp2 = new TProfile("pchftp2","",_jp_npts,&_jp_ptrange[0]);
  pneftp2 = new TProfile("pneftp2","",_jp_npts,&_jp_ptrange[0]);
  pnhftp2 = new TProfile("pnhftp2","",_jp_npts,&_jp_ptrange[0]);
  pceftp2 = new TProfile("pceftp2","",_jp_npts,&_jp_ptrange[0]);
  pmuftp2 = new TProfile("pmuftp2","",_jp_npts,&_jp_ptrange[0]);
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
  hbetaprimetp = new TH1D("hbetaprimetp","",110,0.,1.1);
  // control plots vs NPV
  pncandtp_vsnpv = new TProfile("pncandtp_vsnpv","",_jp_npvs,_jp_pvrange);
  pnchtp_vsnpv = new TProfile("pnchtp_vsnpv","",_jp_npvs,_jp_pvrange);
  pnnetp_vsnpv = new TProfile("pnnetp_vsnpv","",_jp_npvs,_jp_pvrange);
  pnnhtp_vsnpv = new TProfile("pnnhtp_vsnpv","",_jp_npvs,_jp_pvrange);
  pncetp_vsnpv = new TProfile("pncetp_vsnpv","",_jp_npvs,_jp_pvrange);
  pnmutp_vsnpv = new TProfile("pnmutp_vsnpv","",_jp_npvs,_jp_pvrange);
  pchftp_vsnpv = new TProfile("pchftp_vsnpv","",_jp_npvs,_jp_pvrange);
  pneftp_vsnpv = new TProfile("pneftp_vsnpv","",_jp_npvs,_jp_pvrange);
  pnhftp_vsnpv = new TProfile("pnhftp_vsnpv","",_jp_npvs,_jp_pvrange);
  pceftp_vsnpv = new TProfile("pceftp_vsnpv","",_jp_npvs,_jp_pvrange);
  pmuftp_vsnpv = new TProfile("pmuftp_vsnpv","",_jp_npvs,_jp_pvrange);
  pbetatp_vsnpv = new TProfile("pbetatp_vsnpv","",_jp_npvs,_jp_pvrange);
  pbetastartp_vsnpv = new TProfile("pbetastartp_vsnpv","",_jp_npvs,_jp_pvrange);
  pbetaprimetp_vsnpv = new TProfile("pbetaprimetp_vsnpv","",_jp_npvs,_jp_pvrange);
  //
  pchftp_vstrpu = new TProfile("pchftp_vstrpu","",_jp_npvs,_jp_pvrange);
  pneftp_vstrpu = new TProfile("pneftp_vstrpu","",_jp_npvs,_jp_pvrange);
  pnhftp_vstrpu = new TProfile("pnhftp_vstrpu","",_jp_npvs,_jp_pvrange);
  pceftp_vstrpu = new TProfile("pceftp_vstrpu","",_jp_npvs,_jp_pvrange);
  pmuftp_vstrpu = new TProfile("pmuftp_vstrpu","",_jp_npvs,_jp_pvrange);
  pbetatp_vstrpu = new TProfile("pbetatp_vstrpu","",_jp_npvs,_jp_pvrange);
  pbetastartp_vstrpu = new TProfile("pbetastartp_vstrpu","",_jp_npvs,_jp_pvrange);
  pbetaprimetp_vstrpu = new TProfile("pbetaprimetp_vstrpu","",_jp_npvs,_jp_pvrange);
  //
  if (_jp_doPhiHistos) {
    pchfpostp_vsphi = new TProfile("pchfpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pnefpostp_vsphi = new TProfile("pnefpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pnhfpostp_vsphi = new TProfile("pnhfpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pcefpostp_vsphi = new TProfile("pcefpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pmufpostp_vsphi = new TProfile("pmufpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetapostp_vsphi = new TProfile("pbetapostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetastarpostp_vsphi = new TProfile("pbetastarpostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetaprimepostp_vsphi = new TProfile("pbetaprimepostp_vsphi","",_jp_nwphis,_jp_wphirange);
    pchfnegtp_vsphi = new TProfile("pchfnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pnefnegtp_vsphi = new TProfile("pnefnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pnhfnegtp_vsphi = new TProfile("pnhfnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pcefnegtp_vsphi = new TProfile("pcefnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pmufnegtp_vsphi = new TProfile("pmufnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetanegtp_vsphi = new TProfile("pbetanegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetastarnegtp_vsphi = new TProfile("pbetastarnegtp_vsphi","",_jp_nwphis,_jp_wphirange);
    pbetaprimenegtp_vsphi = new TProfile("pbetaprimenegtp_vsphi","",_jp_nwphis,_jp_wphirange);
  }

  // control plots for topology (JEC)
  hselpt = new TH1D("hselpt","",nx,&x[0]);
  hy = new TH1D("hy","",_jp_netas,_jp_etarange); // May 11
  hy2 = new TH1D("hy2","",_jp_netas,_jp_etarange);
  heta = new TH1D("heta","",_jp_netas,_jp_etarange);
  heta2 = new TH1D("heta2","",_jp_netas,_jp_etarange);
  hphi = new TH1D("hphi","",_jp_nphis,_jp_phirange);
  hdphi = new TH1D("hdphi","",_jp_nwposphis,_jp_wposphirange);
  hdpt = new TH1D("hdpt","",100,0.,1.);
  pdpt = new TProfile("pdpt","",nx,&x[0]);
  hjet = new TH1D("hjet","",100,0.,1.);
  hmet = new TH1D("hmet","",100,0.,1.);
  hmetphi = new TH1D("hmetphi","",_jp_nwposphis,_jp_wposphirange);
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
  pr21 = new TProfile("pr21",";npvgood;pt2/pt1",_jp_npvs,_jp_pvrange);
  pr31 = new TProfile("pr31",";npvgood;pt3/pt1",_jp_npvs,_jp_pvrange);
  pr32 = new TProfile("pr32",";npvgood;pt3/pt2",_jp_npvs,_jp_pvrange);
  px21 = new TProfile("px21",";npvgood;has pt2/pt1",_jp_npvs,_jp_pvrange);
  px31 = new TProfile("px31",";npvgood;has pt3/pt1",_jp_npvs,_jp_pvrange);
  px32 = new TProfile("px32",";npvgood;has pt3/pt2",_jp_npvs,_jp_pvrange);

  hyeta = new TH1D("hyeta","",100,-0.436/4.,0.436/4.);
  hyeta2 = new TH1D("hyeta2","",100,-0.436/4.,0.436/4.);
  hbetabetastar = new TH2D("hbetabetastar","",110,0.,1.1,110,0.,1.1);
  hbetabetaprime = new TH2D("hbetabetaprime","",110,0.,1.1,110,0.,1.1);
  hetaphi = new TH2D("hetaphi","",100,-4.799,4.799,_jp_nphis,_jp_phirange);

  pdpt->Sumw2();
  pmpf->Sumw2();
  pmpf1->Sumw2();
  pmpf2->Sumw2();
  pmpfx->Sumw2();

  hdjasymm->Sumw2();
  hdjmpf->Sumw2();
  hdjasymmtp->Sumw2();
  hdjmpftp->Sumw2();

  ppt_probepertag = new TProfile("ppt_probepertag","",nx,&x[0]);

  // MC checks
  htrpu = new TH1D("htrpu","",120,0.,60.); // for PU reweighing
  if (this->ismcdir) {
    hpthat = new TH1D("hpthat","",nx,&x[0]);
    hpthatnlo = new TH1D("hpthatnlo","",nx,&x[0]);

    //unfolding studies (Mikael)
    //mT: (pTgen,ygen); (pTreco,yreco)
    mT = new TH2D("mT","mT(yjet);p_{T,gen};p_{T,reco}",nx,&x[0],nx,&x[0]);
    mTuw = new TH2D("mTuw","mTuw(yjet);p_{T,gen};p_{T,reco}",nx,&x[0],nx,&x[0]);
    mTf = new TH2D("mTf","mT(yjet);p_{T,gen};p_{T,reco}",3485,15,3500, 3485,15,3500);
    mTfuw = new TH2D("mTfuw","mT(yjet);p_{T,gen};p_{T,reco}",3485,15,3500, 3485,15,3500);
    mx = new TH1D("mx","mx(ygen);p_{T,gen}",nx,&x[0]); // pTgen, ygen
    mxuw = new TH1D("mxuw","mx(ygen);p_{T,gen}",nx,&x[0]); // pTgen, ygen
    mxf = new TH1D("mxf","mx(ygen);p_{T,gen}",3485,15,3500); // pTgen, ygen
    mxfuw = new TH1D("mxfuw","mx(ygen);p_{T,gen}",3485,15,3500); // pTgen, ygen
    my = new TH1D("my","my(yreco);p_{T,reco}",nx,&x[0]); // pTreco, yreco
    myuw = new TH1D("myuw","my(yreco);p_{T,reco}",nx,&x[0]); // pTreco, yreco
    myf = new TH1D("myf","my(yreco);p_{T,reco}",3485,15,3500); // pTreco, yreco
    myfuw = new TH1D("myfuw","my(yreco);p_{T,reco}",3485,15,3500); // pTreco, yreco

    hitpu = new TH1D("hitpu","",100,-0.5,99.5);
    hootpuearly = new TH1D("hootpuearly","",100,-0.5,99.5);
    hootpulate = new TH1D("hootpulate","",100,-0.5,99.5);
    h2itvsoot = new TH2D("h2itvsoot","",25,-0.5,24.5,_jp_npvs,_jp_pvrange);

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

    // Response closure
    p2rvsnpv = new TProfile2D("p2rvsnpv","",nx,&x[0],_jp_npvs,_jp_pvrange);
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
    pdy_r = new TProfile2D("pdy_r","",nx,&x[0],_jp_nposphis,_jp_posphirange);
    pdy_g = new TProfile2D("pdy_g","",nx,&x[0],_jp_nposphis,_jp_posphirange);
  } // ismcdir

  curdir->cd();
}

void HistosBasic::Write() {
  dir->cd();
  dir->Write();
}

HistosBasic::~HistosBasic() {
  Write();
};
