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

  unsigned int ieta = int(0.5*(etamin+etamax)/0.5); assert(ieta<jp::noptranges);
  vector<double> x;
  // Chose ptrange according to eta
  for (unsigned int i = 0; i < jp::nopts_eta and jp::ptrangevseta[ieta][i]!=0; ++i) {
    if (jp::ptrangevseta[ieta][i]<0.001) break; // There are zeros in the end of the array when we're out of barrel
    x.push_back(jp::ptrangevseta[ieta][i]);
  }
  const int nx = x.size()-1;

  vector<double> y(51);
  for (unsigned int i = 0; i != y.size(); ++i) y[i] = -5. + 0.2*i;
  const int ny = y.size()-1;

  // raw spectrum
  hpt = new TH1D("hpt","",nx,&x[0]);
  hpt_pre = new TH1D("hpt_pre","",nx,&x[0]); // prescale weighed
  hpt_gtw = new TH1D("hpt_gtw","",nx,&x[0]); // _mc per trigger
  hpt_g0tw = new TH1D("hpt_g0tw","",nx,&x[0]); // _mc per trigger
  hgpt_g0tw = new TH1D("hgpt_g0tw","",nx,&x[0]);//Ozlem
  hqpt_g0tw = new TH1D("hqpt_g0tw","",nx,&x[0]);//Ozlem
  hupt_g0tw = new TH1D("hupt_g0tw","",nx,&x[0]);//// unmatched
  hgpt_g0 = new TH1D("hgpt_g0","",nx,&x[0]);//Ozlem
  hqpt_g0 = new TH1D("hqpt_g0","",nx,&x[0]);//Ozlem
  hupt_g0 = new TH1D("hupt_g0","",nx,&x[0]);// unmatched 
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
  hdjmass0 = new TH1D("hdjmass0","",static_cast<int>(jp::sqrts),0.,jp::sqrts);
  hdjmass_a01 = new TH1D("hdjmass_a01","",nx,&x[0]);
  hdjmass_a02 = new TH1D("hdjmass_a02","",nx,&x[0]);
  hdjmass_a03 = new TH1D("hdjmass_a03","",nx,&x[0]);
  pdjmass_ptratio = new TProfile("pdjmass_ptratio","",nx,&x[0]);
  pdjmass0_ptratio = new TProfile("pdjmass0_ptratio","", static_cast<int>(jp::sqrts),0.,jp::sqrts);
  // leading jet mass
  hjmass = new TH1D("hjmass","",nx,&x[0]);
  hjmass0 = new TH1D("hjmass0","",static_cast<int>(jp::sqrts),0.,jp::sqrts);
  hjmass_a01 = new TH1D("hjmass_a01","",nx,&x[0]);
  hjmass_a02 = new TH1D("hjmass_a02","",nx,&x[0]);
  hjmass_a03 = new TH1D("hjmass_a03","",nx,&x[0]);

  // basic properties
  ppt = new TProfile("ppt","",nx,&x[0]);
  pmass = new TProfile("pmass","",jp::npts,&jp::ptrange[0]);
  hmass = new TH1D("hmass","",100,0.,0.5);
  pjec = new TProfile("pjec","",nx,&x[0]);
  pjec2 = new TProfile("pjec2","",nx,&x[0]);
  punc = new TProfile("punc","",nx,&x[0]);
  hnpvgood = new TH1D("hnpvgood","",jp::npvs,jp::pvrange);
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
  prhovsnpv = new TProfile("prhovsnpv","",jp::npvs,jp::pvrange);
  pnpvvsrho = new TProfile("pnpvvsrho","",jp::npvs,jp::pvrange);
  prhovsnpvall = new TProfile("prhovsnpvall","",jp::npvs,jp::pvrange);
  h2rhovsnpv = new TH2D("h2rhovsnpv","",jp::npvs,jp::pvrange,200,0,40);
  //
  prhovstrpu = new TProfile("prhovstrpu","",jp::npvs,jp::pvrange);
  pnpvvstrpu = new TProfile("pnpvvstrpu","",jp::npvs,jp::pvrange);
  pnpvallvstrpu = new TProfile("pnpvallvstrpu","",jp::npvs,jp::pvrange);
  pitpuvstrpu = new TProfile("itpuvstrpu","",jp::npvs,jp::pvrange);
  htrpu2 = new TH1D("htrpu2","",jp::npvs,jp::pvrange);
  hjet_vstrpu = new TH1D("hjet_vstrpu","",jp::npvs,jp::pvrange);
  hlumi_vstrpu = new TH1D("hlumi_vstrpu","",jp::npvs,jp::pvrange);

  // luminosity
  hlumi = new TH1D("hlumi","",nx,&x[0]);
  hlumi2 = new TH1D("hlumi2","",nx,&x[0]);

  // inclusive efficiencies
  peff = new TProfile("peff","",nx,&x[0]);
  pideff = new TProfile("pideff","",nx,&x[0]);
  pvtxeff = new TProfile("pvtxeff","",nx,&x[0]);
  pdqmeff = new TProfile("pdqmeff","",nx,&x[0]);

  // control plots of components (JEC)
  pncand = new TProfile("pncand","",jp::npts,&jp::ptrange[0]);
  pnch = new TProfile("pnch","",jp::npts,&jp::ptrange[0]);
  pnne = new TProfile("pnne","",jp::npts,&jp::ptrange[0]);
  pnnh = new TProfile("pnnh","",jp::npts,&jp::ptrange[0]);
  pnce = new TProfile("pnce","",jp::npts,&jp::ptrange[0]);
  pnmu = new TProfile("pnmu","",jp::npts,&jp::ptrange[0]);
  pnhh = new TProfile("pnhh","",jp::npts,&jp::ptrange[0]);
  pnhe = new TProfile("pnhe","",jp::npts,&jp::ptrange[0]);
  pchf = new TProfile("pchf","",jp::npts,&jp::ptrange[0]);
  pnef = new TProfile("pnef","",jp::npts,&jp::ptrange[0]);
  pnhf = new TProfile("pnhf","",jp::npts,&jp::ptrange[0]);
  pcef = new TProfile("pcef","",jp::npts,&jp::ptrange[0]);
  pmuf = new TProfile("pmuf","",jp::npts,&jp::ptrange[0]);
  phhf = new TProfile("phhf","",jp::npts,&jp::ptrange[0]);
  phef = new TProfile("phef","",jp::npts,&jp::ptrange[0]);
  ppuf = new TProfile("ppuf","",jp::npts,&jp::ptrange[0]);

  ppttagptprobe = new TProfile("ppttagptprobe","",jp::npts,&jp::ptrange[0]);
  ppttageff = new TProfile("ppttageff","",jp::npts,&jp::ptrange[0]);
  ppttagmu = new TProfile("ppttagmu","",jp::npts,&jp::ptrange[0]);
  ppttagptprobe_noa = new TProfile("ppttagptprobe_noa","",jp::npts,&jp::ptrange[0]);
  ppttagmu_noa = new TProfile("ppttagmu_noa","",jp::npts,&jp::ptrange[0]);
  h2pttagptprobe = new TH2D("h2pttagptprobe","",jp::npts,&jp::ptrange[0],jp::npts,&jp::ptrange[0]);
  h2pttagmu = new TH2D("h2pttagmu","",jp::npts,&jp::ptrange[0],jp::npvs,jp::pvrange);
  h2pttagptprobe_noa = new TH2D("h2pttagptprobe_noa","",jp::npts,&jp::ptrange[0],jp::npts,&jp::ptrange[0]);
  h2pttagmu_noa = new TH2D("h2pttagmu_noa","",jp::npts,&jp::ptrange[0],jp::npvs,jp::pvrange);

  hncand = new TH1D("hncand","",300,-0.5,299.5);
  hnch = new TH1D("hnch","",300,-0.5,299.5);
  hnne = new TH1D("hnne","",300,-0.5,299.5);
  hnnh = new TH1D("hnnh","",300,-0.5,299.5);
  hnce = new TH1D("hnce","",300,-0.5,299.5);
  hnmu = new TH1D("hnmu","",300,-0.5,299.5);
  hnhh = new TH1D("hnhh","",300,-0.5,299.5);
  hnhe = new TH1D("hnhe","",300,-0.5,299.5);
  hchf = new TH1D("hchf","",110,0.,1.1);
  hnef = new TH1D("hnef","",110,0.,1.1);
  hnhf = new TH1D("hnhf","",110,0.,1.1);
  hcef = new TH1D("hcef","",110,0.,1.1);
  hmuf = new TH1D("hmuf","",110,0.,1.1);
  hhhf = new TH1D("hhhf","",110,0.,1.1);
  hhef = new TH1D("hhef","",110,0.,1.1);
  hpuf = new TH1D("hpuf","",110,0.,1.1);
  // control plots of components (JEC tag-and-probe)
  pncandtp = new TProfile("pncandtp","",jp::npts,&jp::ptrange[0]);
  pnchtp = new TProfile("pnchtp","",jp::npts,&jp::ptrange[0]);
  pnnetp = new TProfile("pnnetp","",jp::npts,&jp::ptrange[0]);
  pnnhtp = new TProfile("pnnhtp","",jp::npts,&jp::ptrange[0]);
  pncetp = new TProfile("pncetp","",jp::npts,&jp::ptrange[0]);
  pnmutp = new TProfile("pnmutp","",jp::npts,&jp::ptrange[0]);

  pchftp = new TProfile("pchftp","",jp::npts,&jp::ptrange[0]);
  pneftp = new TProfile("pneftp","",jp::npts,&jp::ptrange[0]);
  pnhftp = new TProfile("pnhftp","",jp::npts,&jp::ptrange[0]);
  pceftp = new TProfile("pceftp","",jp::npts,&jp::ptrange[0]);
  pmuftp = new TProfile("pmuftp","",jp::npts,&jp::ptrange[0]);
  phhftp = new TProfile("phhftp","",jp::npts,&jp::ptrange[0]);
  pheftp = new TProfile("pheftp","",jp::npts,&jp::ptrange[0]);
  ppuftp = new TProfile("ppuftp","",jp::npts,&jp::ptrange[0]);
  //
  hncandtp = new TH1D("hncandtp","",300,-0.5,299.5);
  hnchtp = new TH1D("hnchtp","",300,-0.5,299.5);
  hnnetp = new TH1D("hnnetp","",300,-0.5,299.5);
  hnnhtp = new TH1D("hnnhtp","",300,-0.5,299.5);
  hncetp = new TH1D("hncetp","",300,-0.5,299.5);
  hnmutp = new TH1D("hnmutp","",300,-0.5,299.5);
  hnhhtp = new TH1D("hnhhtp","",300,-0.5,299.5);
  hnhetp = new TH1D("hnhetp","",300,-0.5,299.5);

  hchftp = new TH1D("hchftp","",110,0.,1.1);
  hneftp = new TH1D("hneftp","",110,0.,1.1);
  hnhftp = new TH1D("hnhftp","",110,0.,1.1);
  hceftp = new TH1D("hceftp","",110,0.,1.1);
  hmuftp = new TH1D("hmuftp","",110,0.,1.1);
  hhhftp = new TH1D("hhhftp","",110,0.,1.1);
  hheftp = new TH1D("hheftp","",110,0.,1.1);
  hpuftp = new TH1D("hpuftp","",110,0.,1.1);
  // control plots vs NPV
  pncandtp_vsnpv = new TProfile("pncandtp_vsnpv","",jp::npvs,jp::pvrange);
  pnchtp_vsnpv = new TProfile("pnchtp_vsnpv","",jp::npvs,jp::pvrange);
  pnnetp_vsnpv = new TProfile("pnnetp_vsnpv","",jp::npvs,jp::pvrange);
  pnnhtp_vsnpv = new TProfile("pnnhtp_vsnpv","",jp::npvs,jp::pvrange);
  pncetp_vsnpv = new TProfile("pncetp_vsnpv","",jp::npvs,jp::pvrange);
  pnmutp_vsnpv = new TProfile("pnmutp_vsnpv","",jp::npvs,jp::pvrange);
  pnhhtp_vsnpv = new TProfile("pnhhtp_vsnpv","",jp::npvs,jp::pvrange);
  pnhetp_vsnpv = new TProfile("pnhetp_vsnpv","",jp::npvs,jp::pvrange);
  pchftp_vsnpv = new TProfile("pchftp_vsnpv","",jp::npvs,jp::pvrange);
  pneftp_vsnpv = new TProfile("pneftp_vsnpv","",jp::npvs,jp::pvrange);
  pnhftp_vsnpv = new TProfile("pnhftp_vsnpv","",jp::npvs,jp::pvrange);
  pceftp_vsnpv = new TProfile("pceftp_vsnpv","",jp::npvs,jp::pvrange);
  pmuftp_vsnpv = new TProfile("pmuftp_vsnpv","",jp::npvs,jp::pvrange);
  phhftp_vsnpv = new TProfile("phhftp_vsnpv","",jp::npvs,jp::pvrange);
  pheftp_vsnpv = new TProfile("pheftp_vsnpv","",jp::npvs,jp::pvrange);
  ppuftp_vsnpv = new TProfile("ppuftp_vsnpv","",jp::npvs,jp::pvrange);
  //
  pchftp_vstrpu = new TProfile("pchftp_vstrpu","",jp::npvs,jp::pvrange);
  pneftp_vstrpu = new TProfile("pneftp_vstrpu","",jp::npvs,jp::pvrange);
  pnhftp_vstrpu = new TProfile("pnhftp_vstrpu","",jp::npvs,jp::pvrange);
  pceftp_vstrpu = new TProfile("pceftp_vstrpu","",jp::npvs,jp::pvrange);
  pmuftp_vstrpu = new TProfile("pmuftp_vstrpu","",jp::npvs,jp::pvrange);
  phhftp_vstrpu = new TProfile("phhftp_vstrpu","",jp::npvs,jp::pvrange);
  pheftp_vstrpu = new TProfile("pheftp_vstrpu","",jp::npvs,jp::pvrange);
  ppuftp_vstrpu = new TProfile("ppuftp_vstrpu","",jp::npvs,jp::pvrange);
  //
  if (jp::doPhiHistos) {
    pchfpostp_vsphi = new TProfile("pchfpostp_vsphi","",jp::nwphis,jp::wphirange);
    pnefpostp_vsphi = new TProfile("pnefpostp_vsphi","",jp::nwphis,jp::wphirange);
    pnhfpostp_vsphi = new TProfile("pnhfpostp_vsphi","",jp::nwphis,jp::wphirange);
    pcefpostp_vsphi = new TProfile("pcefpostp_vsphi","",jp::nwphis,jp::wphirange);
    pmufpostp_vsphi = new TProfile("pmufpostp_vsphi","",jp::nwphis,jp::wphirange);
    phhfpostp_vsphi = new TProfile("phhfpostp_vsphi","",jp::nwphis,jp::wphirange);
    phefpostp_vsphi = new TProfile("phefpostp_vsphi","",jp::nwphis,jp::wphirange);
    ppufpostp_vsphi = new TProfile("ppufpostp_vsphi","",jp::nwphis,jp::wphirange);
    pchfnegtp_vsphi = new TProfile("pchfnegtp_vsphi","",jp::nwphis,jp::wphirange);
    pnefnegtp_vsphi = new TProfile("pnefnegtp_vsphi","",jp::nwphis,jp::wphirange);
    pnhfnegtp_vsphi = new TProfile("pnhfnegtp_vsphi","",jp::nwphis,jp::wphirange);
    pcefnegtp_vsphi = new TProfile("pcefnegtp_vsphi","",jp::nwphis,jp::wphirange);
    pmufnegtp_vsphi = new TProfile("pmufnegtp_vsphi","",jp::nwphis,jp::wphirange);
    phhfnegtp_vsphi = new TProfile("phhfnegtp_vsphi","",jp::nwphis,jp::wphirange);
    phefnegtp_vsphi = new TProfile("phefnegtp_vsphi","",jp::nwphis,jp::wphirange);
    ppufnegtp_vsphi = new TProfile("ppufnegtp_vsphi","",jp::nwphis,jp::wphirange);
  }

  //{ Histograms for QGL study (Ozlem)
  hgpt = new TH1D("hgpt","",nx,&x[0]);
  hgpt0 = new TH1D("hgpt0","",6500,0.,6500.);

  hqgl = new TH1D("hqgl","",101,0,1.01);
  hqgl_q = new TH1D("hqgl_q","",101,0,1.01);
  hqgl_g = new TH1D("hqgl_g","",101,0,1.01);
  hqgl_u = new TH1D("hqgl_u","",101,0,1.01);
  hqgl_dq = new TH1D("hqgl_dq","",101,0,1.01);
  hqgl_dg = new TH1D("hqgl_dg","",101,0,1.01);
  hqgl_du = new TH1D("hqgl_du","",101,0,1.01);

  // add TH2Ds for pT bins here
  hqgl2 =  new TH2D("hqgl2","",nx,&x[0],101,0,1.01);
  hqgl2_g =  new TH2D("hqgl2_g","",nx,&x[0],101,0,1.01);
  hqgl2_q =  new TH2D("hqgl2_q","",nx,&x[0],101,0,1.01);
  hqgl2_u =  new TH2D("hqgl2_u","",nx,&x[0],101,0,1.01);
  hqgl2_dg =  new TH2D("hqgl2_dg","",nx,&x[0],101,0,1.01);
  hqgl2_dq =  new TH2D("hqgl2_dq","",nx,&x[0],101,0,1.01);
  hqgl2_du =  new TH2D("hqgl2_du","",nx,&x[0],101,0,1.01);
  //} QGL study

  // control plots for topology (JEC)
  hselpt = new TH1D("hselpt","",nx,&x[0]);
  hy = new TH1D("hy","",jp::netas,jp::etarange); // May 11
  hy2 = new TH1D("hy2","",jp::netas,jp::etarange);
  heta = new TH1D("heta","",jp::netas,jp::etarange);
  heta2 = new TH1D("heta2","",jp::netas,jp::etarange);
  hphi = new TH1D("hphi","",jp::nphis,jp::phirange);
  hdphi = new TH1D("hdphi","",jp::nwposphis,jp::wposphirange);
  hdpt = new TH1D("hdpt","",100,0.,1.);
  pdpt = new TProfile("pdpt","",nx,&x[0]);
  hjet = new TH1D("hjet","",100,0.,1.);
  hmet = new TH1D("hmet","",100,0.,1.);
  hmetphi = new TH1D("hmetphi","",jp::nwposphis,jp::wposphirange);

  hnjet_pre       = new TH1D("hnjet_pre",";Njets;",15,0,15);
  hj0pt_pre       = new TH1D("hj0pt_pre",";p_{T,leading} (GeV);",nx,&x[0]);
  hj0phi_pre      = new TH1D("hj0phi_pre",";#phi_{leading};",jp::nphis,jp::phirange);
  hj1pt_pre       = new TH1D("hj1pt_pre",";p_{T,j2} (GeV);",200,0,200);
  hj2pt_pre       = new TH1D("hj2pt_pre",";p_{T,j3} (GeV);",200,0,200);
  hj0genpt_pre       = new TH1D("hj0genpt_pre",";p_{T,leading,gen} (GeV);",nx,&x[0]);
  hj0genr_pre        = new TH1D("hj0genr_pre",";;",100,0.0,0.5);

  hresp_jetsel       = new TH1D("hresp_jetsel",";Resp(jetsel,gen);",100,0,2.5);
  hdphi_recoiljetsel = new TH1D("hdphi_recoiljetsel",";#Delta#phi(leading,recoiljets);",50,1,TMath::Pi());
  //hweight_recoiljetsel  = new TH1D("hweight_recoiljetsel",";p_{T,recoil} #times cos(#Delta#phi(leading,recoiljet);",100,0,10);
  //h2etapt_recoiljetsel  = new TH2D("h2etapt_recoiljetsel",";#eta_{recoiljet};p_{T,recoiljet};",100,-5,5,nx,&x[0]);
  //h2etapt_recoiljetsel_weight  = new TH2D("h2etapt_recoiljetsel_weight",";#eta_{recoiljet};p_{T,recoiljet};",100,-5,5,nx,&x[0]);
  hdphi_recoilmin    = new TH1D("hdphi_recoilmin",";Min of #Delta#phi(leading,recoiljets);",jp::nwposphis,jp::wposphirange);
  hdphi_recoilmax    = new TH1D("hdphi_recoilmax",";Max of #Delta#phi(leading,recoiljets);",jp::nwposphis,jp::wposphirange);
  hdphi_recoilcol    = new TH2D("hdphi_recoilcol",";#Delta#phi_{min};#Delta#phi_{max}",jp::nwposphis,jp::wposphirange,jp::nwposphis,jp::wposphirange);
  pjtgenptsf_recoil  = new TProfile("pjtgenptsf_recoil",";p_{T,gen} of jets in the recoil;JER SF",nx,&x[0]);
  pjtptsf_recoil     = new TProfile("pjtptsf_recoil",";p_{T,reco} of jets in the recoil;JER SF",nx,&x[0]);

  //double minR = -3., maxR = 3.;
  for (unsigned int i = 0; i < 3; ++i) {

    hdphi_recoils[i]     = new TH1D(Form("hdphi_recoils%d",i),";#Delta#phi(recoils,recoiljets);",jp::nwposphis,jp::wposphirange); 
    hdphi_recoils_test[i]     = new TH1D(Form("hdphi_recoils_test%d",i),";#Delta#phi(recoils,recoiljets);",jp::nwposphis,jp::wposphirange);
    h2pt_dphi_recoils_test[i] = new TH2D(Form("h2pt_dphi_recoils_test%d",i),";p_{T,recoil};#Delta#phi(recoils,recoiljets)",nx,&x[0],jp::nwposphis,jp::wposphirange);
    h2pt_dphi_recoils[i] = new TH2D(Form("h2pt_dphi_recoils%d",i),";p_{T,recoil};#Delta#phi(recoils,recoiljets)",nx,&x[0],jp::nwposphis,jp::wposphirange);
    h2nrecoiljets[i]     = new TH2D(Form("h2nrecoiljets%d",i),";Nrecoiljets;Ngoodrecoils",15,0,15,15,0,15);
    hrecoilrat_good[i]   = new TH1D(Form("hrecoilrat_good%d",i),";p_{T,j2}/p_{T,recoil};",100,0.,2.);
    hrecoilrat_fail[i]   = new TH1D(Form("hrecoilrat_fail%d",i),";p_{T,j2}/p_{T,recoil};",100,0.,2.);
    if(jp::ismc) {
      hdphi_recoils_weight[i]     = new TH1D(Form("hdphi_recoils_weight%d",i),";#Delta#phi(recoils,recoiljets);",jp::nwposphis,jp::wposphirange);
      h2pt_dphi_recoils_weight[i] = new TH2D(Form("h2pt_dphi_recoils_weight%d",i),";p_{T,recoil};#Delta#phi(recoils,recoiljets)",nx,&x[0],jp::nwposphis,jp::wposphirange);
    }

    hnjet[i]          = new TH1D(Form("hnjet%d",i),";Njets;",15,0,15);
    hj0pt[i]          = new TH1D(Form("hj0pt%d",i),";p_{T,leading} (GeV);",nx,&x[0]);
    hj0phi[i]         = new TH1D(Form("hj0phi%d",i),";#phi_{leading};",jp::nphis,jp::phirange);
    hj1pt[i]          = new TH1D(Form("hj1pt%d",i),";p_{T,j2} (GeV);",nx,&x[0]);
    hj2pt[i]          = new TH1D(Form("hj2pt%d",i),";p_{T,j3} (GeV);",nx,&x[0]);
    hrecoilpt[i]      = new TH1D(Form("hrecoilpt%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hrecoilphi[i]     = new TH1D(Form("hrecoilphi%d",i),";#phi_{recoil};",jp::nphis,jp::phirange);
    hrecoilrat[i]     = new TH1D(Form("hrecoilrat%d",i),";p_{T,j2}/p_{T,recoil};",100,0.,2.);
    hrecoilrat2[i]    = new TH1D(Form("hrecoilrat2%d",i),";p_{T,j3}/p_{T,recoil};",100,0.,2.);
    hptleadrat[i]            = new TH1D(Form("hptleadrat%d",i),";p_{T,j2}/p_{T,lead};",100,0.,2.);
    hptaverat[i]             = new TH1D(Form("hptaverat%d",i),";p_{T,j2}/p_{T,ave};",100,0.,2.);
    hptleadrat_pre[i]        = new TH1D(Form("hptleadrat_pre%d",i),";p_{T,j2}/p_{T,lead};",100,0.,2.);
    hptaverat_pre[i]         = new TH1D(Form("hptaverat_pre%d",i),";p_{T,j2}/p_{T,ave};",100,0.,2.);
    h2ptrecoilrat[i]         = new TH2D(Form("h2ptrecoilrat%d",i),";p_{T,j2}/p_{T,recoil};p_{T,j3}/p_{T,recoil}",100,0.,2.,100,0.,2.);
    h2ptrecoilrat_pre[i]     = new TH2D(Form("h2ptrecoilrat_pre%d",i),";p_{T,j2}/p_{T,recoil};p_{T,j3}/p_{T,recoil}",100,0.,2.,100,0.,2.);
    hrecoildphi[i]    = new TH1D(Form("hrecoildphi%d",i),";#Delta#phi(leading,recoil);",jp::nwposphis,jp::wposphirange);
    hrecoiljets[i]    = new TH1D(Form("hrecoiljets%d",i),";N_{recoiljets};",15,0,15);
    hrecoilmjb[i]     = new TH1D(Form("hrecoilmjb%d",i),";p_{T,leading}/p_{T,recoil};",100,0.,4.0);
    hrecoilmjbinv[i]  = new TH1D(Form("hrecoilmjbinv%d",i),";p_{T,recoil}/p_{T,leading};",100,0.,4.0);
    hrecoilmpftpinv[i]  = new TH1D(Form("hrecoilmpftpinv%d",i),";1/MPF;",100,0.,4.);
    hrecoilmpftpinv2[i] = new TH1D(Form("hrecoilmpftpinv2%d",i),";1/MPF;",100,0.,4.);
    hrecoilmpf[i]     = new TH1D(Form("hrecoilmpf%d",i),"",100,-100,100);
    hrecoilmpftp[i]   = new TH1D(Form("hrecoilmpftp%d",i),";MPF(MET1);",100,-1.,3.);
    hrecoilmpf2[i]    = new TH1D(Form("hrecoilmpf2%d",i),"",100,-100,100);
    hrecoilmpftp2[i]  = new TH1D(Form("hrecoilmpftp2%d",i),";MPF(MET2);",100,-1.,3.);
    hleadmpf[i]       = new TH1D(Form("hleadmpf%d",i),"",100,-100,100);
    hleadmpftp[i]     = new TH1D(Form("hleadmpftp%d",i),";MPF(MET1,leading);",100,-1.,3.);
    hleadmpf2[i]      = new TH1D(Form("hleadmpf2%d",i),"",100,-100,100);
    hleadmpftp2[i]    = new TH1D(Form("hleadmpftp2%d",i),";MPF(MET2,leading);",100,-1.,3.);
    pmjb_recoil[i]    = new TProfile(Form("pmjb_recoil%d",i),";p_{T,recoil} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
    pmjb_leading[i]   = new TProfile(Form("pmjb_leading%d",i),";p_{T,leading} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
    pmjbinv_leading[i]  = new TProfile(Form("pmjbinv_leading%d",i),";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
    pmpfinv_leading[i]  = new TProfile(Form("pmpfinv_leading%d",i),";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
    pmpfinv_leading2[i] = new TProfile(Form("pmpfinv_leading2%d",i),";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
    pmpf_recoil[i]    = new TProfile(Form("pmpf_recoil%d",i),";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
    pmpf_recoil2[i]   = new TProfile(Form("pmpf_recoil2%d",i),";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
    pmpf_leading[i]   = new TProfile(Form("pmpf_leading%d",i),";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
    pmpf_leading2[i]  = new TProfile(Form("pmpf_leading2%d",i),";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
    pmpflead_leading[i]   = new TProfile(Form("pmpflead_leading%d",i),";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
    pmpflead_leading2[i]  = new TProfile(Form("pmpflead_leading2%d",i),";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
    pmpfleadinv_leading[i]   = new TProfile(Form("pmpfleadinv_leading%d",i),";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
    pmpfleadinv_leading2[i]  = new TProfile(Form("pmpfleadinv_leading2%d",i),";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
    pmjbinv_recoil[i]   = new TProfile(Form("pmjbinv_recoil%d",i),";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
    pmpflead_recoil[i]  = new TProfile(Form("pmpflead_recoil%d",i),";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
    pmpflead_recoil2[i] = new TProfile(Form("pmpflead_recoil2%d",i),";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
    h2mjb_recoil[i]       = new TH2D(Form("h2mjb_recoil%d",i),";p_{T,recoil} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
    h2mjb_leading[i]      = new TH2D(Form("h2mjb_leading%d",i),";p_{T,leading} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
    h2mjbinv_leading[i]   = new TH2D(Form("h2mjbinv_leading%d",i),";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
    h2mpfinv_leading[i]   = new TH2D(Form("h2mpfinv_leading%d",i),";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
    h2mpfinv_leading2[i]  = new TH2D(Form("h2mpfinv_leading2%d",i),";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
    h2mpf_recoil[i]       = new TH2D(Form("h2mpf_recoil%d",i),";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpf_recoil2[i]      = new TH2D(Form("h2mpf_recoil2%d",i),";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpf_leading[i]      = new TH2D(Form("h2mpf_leading%d",i),";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpf_leading2[i]     = new TH2D(Form("h2mpf_leading2%d",i),";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpflead_leading[i]  = new TH2D(Form("h2mpflead_leading%d",i),";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
    h2mpflead_leading2[i] = new TH2D(Form("h2mpflead_leading2%d",i),";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
    h2mpfleadinv_leading[i]  = new TH2D(Form("h2mpfleadinv_leading%d",i),";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
    h2mpfleadinv_leading2[i] = new TH2D(Form("h2mpfleadinv_leading2%d",i),";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
    h2mjbinv_recoil[i]   = new TH2D(Form("h2mjbinv_recoil%d",i),";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
    h2mpflead_recoil[i]  = new TH2D(Form("h2mpflead_recoil%d",i),";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
    h2mpflead_recoil2[i] = new TH2D(Form("h2mpflead_recoil2%d",i),";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
 
    hj0genpt[i]          = new TH1D(Form("hj0genpt%d",i),";p_{T,leading,gen} (GeV);",nx,&x[0]);
    hj0genr[i]           = new TH1D(Form("hj0genr%d",i),";;",100,0.0,0.5);
    hrecoilgenpt[i]      = new TH1D(Form("hrecoilgenpt%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hrecoilgenphi[i]     = new TH1D(Form("hrecoilgenphi%d",i),";#phi_{recoil};",jp::nphis,jp::phirange);
    hresp_leading[i]     = new TH1D(Form("hresp_leading%d",i),";Resp(leading,gen);",100,0,2.5);
    hresp_recoil[i]      = new TH1D(Form("hresp_recoil%d",i),";Resp(recoil,gen);",100,0,5);
    hresp_mjb[i]         = new TH1D(Form("hresp_mjb%d",i),";Resp(mjb,gen);",100,0,2.5);
    presp_leading[i]     = new TProfile(Form("presp_leading%d",i),"",nx,&x[0]);
    presp_recoil[i]      = new TProfile(Form("presp_recoil%d",i),"",nx,&x[0]);
    hdphi_jet12[i]       = new TH1D(Form("hdphi_jet12%d",i),";#Delta#phi(j1,j2);",jp::nwposphis,jp::wposphirange);
    hdphi_jet23[i]       = new TH1D(Form("hdphi_jet23%d",i),";#Delta#phi(j2,j3);",jp::nwposphis,jp::wposphirange);
    hCRecoil[i]          = new TH1D(Form("hCRecoil%d",i),";CRecoil;",100,0,1);
    hCRecoil_ptcl[i]     = new TH1D(Form("hCRecoil_ptcl%d",i),";CRecoil_ptcl;",100,0,1);
    pCRecoil[i]          = new TProfile(Form("pCRecoil%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    pCRecoil_leading[i]  = new TProfile(Form("pCRecoil_leading%d",i),";p_{T,leading} (GeV);",nx,&x[0]);
    pCRecoil_ptave[i]    = new TProfile(Form("pCRecoil_ptave%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    pCRecoil_ptcl[i]     = new TProfile(Form("pCRecoil_ptcl%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    pCRecoil_ptcl2[i]    = new TProfile(Form("pCRecoil_ptcl2%d",i),";p_{T,recoil_gen} (GeV);",nx,&x[0]);
    pCRecoil_ptcl2_leading[i] = new TProfile(Form("pCRecoil_ptcl2_leading%d",i),";p_{T,leading_gen} (GeV);",nx,&x[0]);
    pCRecoil_ptcl2_ptave[i]   = new TProfile(Form("pCRecoil_ptcl2_ptave%d",i),";p_{T,ave_gen} (GeV);",nx,&x[0]);
    pjtgenptsf[i]         = new TProfile(Form("pjtgenptsf%d",i),";p_{T,leading,gen};JER SF",nx,&x[0]);
    pjtptsf[i]            = new TProfile(Form("pjtptsf%d",i),";p_{T,leading,reco};JER SF",nx,&x[0]);
 
    pmjb_ptave[i]       = new TProfile(Form("pmjb_ptave%d",i),";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
    pmjbinv_ptave[i]    = new TProfile(Form("pmjbinv_ptave%d",i),";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
    pmpf_ptave[i]       = new TProfile(Form("pmpf_ptave%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
    pmpf_ptave2[i]      = new TProfile(Form("pmpf_ptave2%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
    pmpfinv_ptave[i]    = new TProfile(Form("pmpfinv_ptave%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
    pmpfinv_ptave2[i]   = new TProfile(Form("pmpfinv_ptave2%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
 
    pmpf_ptave_unc[i]    = new TProfile(Form("pmpf_ptave_unc%d",i),";p_{T,ave} (GeV);R_{MPFu}",nx,&x[0]);
    pmpf_ptave_one[i]   = new TProfile(Form("pmpf_ptave_one%d",i),";p_{T,ave} (GeV);R_{MPF1}",nx,&x[0]);
    pmpf_ptave_n[i]     = new TProfile(Form("pmpf_ptave_n%d",i),";p_{T,ave} (GeV);R_{MPFn}",nx,&x[0]);
    pmpfinv_ptave_unc[i]  = new TProfile(Form("pmpfinv_ptave_unc%d",i),";p_{T,ave} (GeV);R_{MPFu}",nx,&x[0]);
    pmpfinv_ptave_one[i] = new TProfile(Form("pmpfinv_ptave_one%d",i),";p_{T,ave} (GeV);R_{MPF1}",nx,&x[0]);
    pmpfinv_ptave_n[i]   = new TProfile(Form("pmpfinv_ptave_n%d",i),";p_{T,ave} (GeV);R_{MPFn}",nx,&x[0]);
 
    h2mjb_ptave[i]       = new TH2D(Form("h2mjb_ptave%d",i),";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
    h2mjbinv_ptave[i]    = new TH2D(Form("h2mjbinv_ptave%d",i),";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
    h2mpf_ptave[i]       = new TH2D(Form("h2mpf_ptave%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpf_ptave2[i]      = new TH2D(Form("h2mpf_ptave2%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpfinv_ptave[i]    = new TH2D(Form("h2mpfinv_ptave%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
    h2mpfinv_ptave2[i]   = new TH2D(Form("h2mpfinv_ptave2%d",i),";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
 
    h2mpf_ptave_unc[i]    = new TH2D(Form("h2mpf_ptave_unc%d",i),";p_{T,ave} (GeV);R_{MPFu}",nx,&x[0],100,-4,6);
    h2mpf_ptave_one[i]   = new TH2D(Form("h2mpf_ptave_one%d",i),";p_{T,ave} (GeV);R_{MPF1}",nx,&x[0],100,-4,6);
    h2mpf_ptave_n[i]     = new TH2D(Form("h2mpf_ptave_n%d",i),";p_{T,ave} (GeV);R_{MPFn}",nx,&x[0],100,-4,6);
 
    h2mpfinv_ptave_unc[i]  = new TH2D(Form("h2mpfinv_ptave_unc%d",i),";p_{T,ave} (GeV);R_{MPFu}",nx,&x[0],100,-4,6);
    h2mpfinv_ptave_one[i] = new TH2D(Form("h2mpfinv_ptave_one%d",i),";p_{T,ave} (GeV);R_{MPF1}",nx,&x[0],100,-4,6);
    h2mpfinv_ptave_n[i]   = new TH2D(Form("h2mpfinv_ptave_n%d",i),";p_{T,ave} (GeV);R_{MPFn}",nx,&x[0],100,-4,6);
 
    //hleadmpf_unc[i]          = new TH1D(Form("hleadmpf_unc%d",i),"",100,-100,100);
    //hleadmpf_one[i]         = new TH1D(Form("hleadmpf_one%d",i),"",100,-100,100);
    //hleadmpf_n[i]           = new TH1D(Form("hleadmpf_n%d",i),"",100,-100,100);
    //hleadmpftp_unc[i]       = new TH1D(Form("hleadmpftp_unc%d",i),";p_{T,leading} (GeV);MPF_{u}",100,-3.,3.);
    //hleadmpftp_one[i]       = new TH1D(Form("hleadmpftp_one%d",i),";p_{T,leading} (GeV);MPF_{1};",100,-3.,3.);
    //hleadmpftp_n[i]         = new TH1D(Form("hleadmpftp_n%d",i),";p_{T,leading} (GeV);MPF_{n};",100,-3.,3.);
 
    pmpflead_leading_unc[i]  = new TProfile(Form("pmpflead_leading_unc%d",i),";p_{T,leading} (GeV);MPF_{u}",nx,&x[0]);
    pmpflead_leading_one[i] = new TProfile(Form("pmpflead_leading_one%d",i),";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0]);
    pmpflead_leading_n[i]   = new TProfile(Form("pmpflead_leading_n%d",i),";p_{T,leading} (GeV);MPF_{n}",nx,&x[0]);
 
    h2mpflead_leading_unc[i]  = new TH2D(Form("h2mpflead_leading_unc%d",i),";p_{T,leading} (GeV);MPF_{u}",nx,&x[0],100,-5,5);
    h2mpflead_leading_one[i] = new TH2D(Form("h2mpflead_leading_one%d",i),";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0],100,-5,5);
    h2mpflead_leading_n[i]   = new TH2D(Form("h2mpflead_leading_n%d",i),";p_{T,leading} (GeV);MPF_{n}",nx,&x[0],100,-5,5);
 
    pmpf_recoil_unc[i]  = new TProfile(Form("pmpf_recoil_unc%d",i),";p_{T,recoil} (GeV);MPF_{u}",nx,&x[0]);
    pmpf_recoil_one[i] = new TProfile(Form("pmpf_recoil_one%d",i),";p_{T,recoil} (GeV);1 - MPF_{1}",nx,&x[0]);
    pmpf_recoil_n[i]   = new TProfile(Form("pmpf_recoil_n%d",i),";p_{T,recoil} (GeV);MPF_{n}",nx,&x[0]);
 
    h2mpf_recoil_unc[i]  = new TH2D(Form("h2mpf_recoil_unc%d",i),";p_{T,recoil} (GeV);MPF_{u}",nx,&x[0],100,-5,5);
    h2mpf_recoil_one[i] = new TH2D(Form("h2mpf_recoil_one%d",i),";p_{T,recoil} (GeV);1 - MPF_{1}",nx,&x[0],100,-5,5);
    h2mpf_recoil_n[i]   = new TH2D(Form("h2mpf_recoil_n%d",i),";p_{T,recoil} (GeV);MPF_{n}",nx,&x[0],100,-5,5);

    prho_recoil[i]       = new TProfile(Form("prho_recoil%d",i),";p_{T,recoil} (GeV);#rho",nx,&x[0]);
    prho_leading[i]      = new TProfile(Form("prho_leading%d",i),";p_{T,leading} (GeV);#rho",nx,&x[0]);
    prho_ptave[i]        = new TProfile(Form("prho_ptave%d",i),";p_{T,ave} (GeV);#rho",nx,&x[0]);
    h2rho_recoil[i]      = new TH2D(Form("h2rho_recoil%d",i),";p_{T,recoil} (GeV);#rho",nx,&x[0],200,0,80);
    h2rho_leading[i]     = new TH2D(Form("h2rho_leading%d",i),";p_{T,leading} (GeV);#rho",nx,&x[0],200,0,80);
    h2rho_ptave[i]       = new TH2D(Form("h2rho_ptave%d",i),";p_{T,ave} (GeV);#rho",nx,&x[0],200,0,80);

    //pmjb_ptave[i]->BuildOptions(minR,maxR,"i");
    //pmjbinv_ptave[i]->BuildOptions(minR,maxR,"i");
    //pmpf_ptave[i]->BuildOptions(minR,maxR,"i");       
    //pmpfinv_ptave[i]->BuildOptions(minR,maxR,"i");
    //pmpf_ptave2[i]->BuildOptions(minR,maxR,"i");
    //pmpfinv_ptave2[i]->BuildOptions(minR,maxR,"i");

  }

  // multijet flavour
  hptave_all.resize(6); hptlead_all.resize(6); hptrecoil_all.resize(6);
  hptave_weight_all.resize(6); hptlead_weight_all.resize(6); hptrecoil_weight_all.resize(6);

  hptave_b.resize(9); hptave_c.resize(9); hptave_g.resize(9); hptave_ud.resize(9); hptave_s.resize(9); //hptave_no.resize(9);
  hptave_weight_b.resize(9); hptave_weight_c.resize(9); hptave_weight_g.resize(9); hptave_weight_ud.resize(9); hptave_weight_s.resize(9); //hptave_weight_no.resize(9);
  hptlead_b.resize(9); hptlead_c.resize(9); hptlead_g.resize(9); hptlead_ud.resize(9); hptlead_s.resize(9); //hptlead_no.resize(9);
  hptlead_weight_b.resize(9); hptlead_weight_c.resize(9); hptlead_weight_g.resize(9); hptlead_weight_s.resize(9); hptlead_weight_ud.resize(9); //hptlead_weight_no.resize(9);
  hptrecoil_b.resize(9); hptrecoil_c.resize(9); hptrecoil_g.resize(9); hptrecoil_ud.resize(9); hptrecoil_s.resize(9); //hptrecoil_no.resize(9);
  hptrecoil_weight_b.resize(9); hptrecoil_weight_c.resize(9); hptrecoil_weight_g.resize(9); hptrecoil_weight_ud.resize(9); hptrecoil_weight_s.resize(9); //hptrecoil_weight_no.resize(9);

  hptave_phys_b.resize(9); hptave_phys_c.resize(9); hptave_phys_g.resize(9); hptave_phys_ud.resize(9); hptave_phys_s.resize(9); hptave_phys_no.resize(9);
  hptlead_phys_b.resize(9); hptlead_phys_c.resize(9); hptlead_phys_g.resize(9); hptlead_phys_ud.resize(9); hptlead_phys_s.resize(9); hptlead_phys_no.resize(9);
  hptrecoil_phys_b.resize(9); hptrecoil_phys_c.resize(9); hptrecoil_phys_g.resize(9); hptrecoil_phys_ud.resize(9); hptrecoil_phys_s.resize(9); hptrecoil_phys_no.resize(9);

  for (unsigned int i = 0; i < hptave_b.size(); ++i) {
    if(i<6) {
      hptave_all[i]            = new TH1D(Form("hptave_all%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
      hptlead_all[i]           = new TH1D(Form("hptlead_all%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
      hptrecoil_all[i]         = new TH1D(Form("hptrecoil_all%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
      hptave_weight_all[i]     = new TH1D(Form("hptave_weight_all%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
      hptlead_weight_all[i]    = new TH1D(Form("hptlead_weight_all%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
      hptrecoil_weight_all[i]  = new TH1D(Form("hptrecoil_weight_all%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);

      //h2ptave_all[i]         = new TH2D(Form("h2ptave_all%d",i),";p_{T,ave} (GeV);eta",nx,&x[0],36,jp::wetarange);
      //h2ptlead_all[i]        = new TH2D(Form("h2ptlead_all%d",i),";p_{T,lead} (GeV);eta",nx,&x[0],36,jp::wetarange);
      //h2ptrecoil_all[i]      = new TH2D(Form("h2ptrecoil_all%d",i),";p_{T,recoil} (GeV);eta",nx,&x[0],36,jp::etarange);
    }
    hptave_b[i]              = new TH1D(Form("hptave_b%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_c[i]              = new TH1D(Form("hptave_c%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_g[i]              = new TH1D(Form("hptave_g%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_s[i]              = new TH1D(Form("hptave_s%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_ud[i]            = new TH1D(Form("hptave_ud%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    //hptave_no[i]             = new TH1D(Form("hptave_no%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_weight_b[i]       = new TH1D(Form("hptave_weight_b%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_weight_c[i]       = new TH1D(Form("hptave_weight_c%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_weight_g[i]       = new TH1D(Form("hptave_weight_g%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_weight_s[i]       = new TH1D(Form("hptave_weight_s%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_weight_ud[i]     = new TH1D(Form("hptave_weight_ud%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    //hptave_weight_no[i]      = new TH1D(Form("hptave_weight_no%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptlead_b[i]             = new TH1D(Form("hptlead_b%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_c[i]             = new TH1D(Form("hptlead_c%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_g[i]             = new TH1D(Form("hptlead_g%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_s[i]             = new TH1D(Form("hptlead_s%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_ud[i]           = new TH1D(Form("hptlead_ud%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    //hptlead_no[i]            = new TH1D(Form("hptlead_no%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_weight_b[i]      = new TH1D(Form("hptlead_weight_b%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_weight_c[i]      = new TH1D(Form("hptlead_weight_c%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_weight_g[i]      = new TH1D(Form("hptlead_weight_g%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_weight_s[i]      = new TH1D(Form("hptlead_weight_s%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_weight_ud[i]    = new TH1D(Form("hptlead_weight_ud%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    //hptlead_weight_no[i]     = new TH1D(Form("hptlead_weight_no%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptrecoil_b[i]           = new TH1D(Form("hptrecoil_b%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_c[i]           = new TH1D(Form("hptrecoil_c%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_g[i]           = new TH1D(Form("hptrecoil_g%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_s[i]           = new TH1D(Form("hptrecoil_s%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_ud[i]         = new TH1D(Form("hptrecoil_ud%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    //hptrecoil_no[i]          = new TH1D(Form("hptrecoil_no%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_weight_b[i]    = new TH1D(Form("hptrecoil_weight_b%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_weight_c[i]    = new TH1D(Form("hptrecoil_weight_c%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_weight_g[i]    = new TH1D(Form("hptrecoil_weight_g%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_weight_s[i]    = new TH1D(Form("hptrecoil_weight_s%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_weight_ud[i]  = new TH1D(Form("hptrecoil_weight_ud%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    //hptrecoil_weight_no[i]   = new TH1D(Form("hptrecoil_weight_no%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);

    hptave_phys_b[i]       = new TH1D(Form("hptave_phys_b%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_phys_c[i]       = new TH1D(Form("hptave_phys_c%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_phys_g[i]       = new TH1D(Form("hptave_phys_g%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_phys_s[i]       = new TH1D(Form("hptave_phys_s%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_phys_ud[i]      = new TH1D(Form("hptave_phys_ud%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptave_phys_no[i]      = new TH1D(Form("hptave_phys_no%d",i),";p_{T,ave} (GeV);",nx,&x[0]);
    hptlead_phys_b[i]      = new TH1D(Form("hptlead_phys_b%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_phys_c[i]      = new TH1D(Form("hptlead_phys_c%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_phys_g[i]      = new TH1D(Form("hptlead_phys_g%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_phys_s[i]      = new TH1D(Form("hptlead_phys_s%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_phys_ud[i]     = new TH1D(Form("hptlead_phys_ud%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptlead_phys_no[i]     = new TH1D(Form("hptlead_phys_no%d",i),";p_{T,lead} (GeV);",nx,&x[0]);
    hptrecoil_phys_b[i]    = new TH1D(Form("hptrecoil_phys_b%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_phys_c[i]    = new TH1D(Form("hptrecoil_phys_c%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_phys_g[i]    = new TH1D(Form("hptrecoil_phys_g%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_phys_s[i]    = new TH1D(Form("hptrecoil_phys_s%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_phys_ud[i]   = new TH1D(Form("hptrecoil_phys_ud%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);
    hptrecoil_phys_no[i]   = new TH1D(Form("hptrecoil_phys_no%d",i),";p_{T,recoil} (GeV);",nx,&x[0]);

  }

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
  hmpfinv = new TH1D("hmpfinv","",200,0.,2.); //mpf inverse added on Mar 31, 2020
  hmpfinv1 = new TH1D("hmpfinv1","",200,0.,2.);
  hmpfinv2 = new TH1D("hmpfinv2","",200,0.,2.);
  hmpfx = new TH1D("hmpfx","",200,0.,2.);
  hmpfy = new TH1D("hmpfy","",200,0.,2.);
  hmpfz = new TH1D("hmpfz","",200,0.,2.);
  pmpf = new TProfile("pmpf","",nx,&x[0]);
  pmpf1 = new TProfile("pmpf1","",nx,&x[0]);
  pmpf2 = new TProfile("pmpf2","",nx,&x[0]);
  pmpfinv = new TProfile("pmpfinv","",nx,&x[0]);
  pmpfinv1 = new TProfile("pmpfinv1","",nx,&x[0]);
  pmpfinv2 = new TProfile("pmpfinv2","",nx,&x[0]);
  pmpfx = new TProfile("pmpfx","",nx,&x[0]);
  pmpfy = new TProfile("pmpfy","",nx,&x[0]);
  pmpfz = new TProfile("pmpfz","",nx,&x[0]);

  const int n3 = 40;
  vector<double> v3(n3+1);
  for (unsigned int i = 0; i != n3+1; ++i) v3[i] = 0. + 1.*i/n3;
  const int na = 200;
  vector<double> va(na+1);
  for (unsigned int i = 0; i != na+1; ++i) va[i] = -1. + 2.*i/na;
if (jp::do3dHistos) {
  hdjasymm = new TH3D("hdjasymm",";p_{T,ave};p_{T,3rd}/p_{T,ave};Asymmetry",
                      nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjmpf = new TH3D("hdjmpf",";p_{T,ave};p_{T,3rd}/p_{T,ave};MPF",
                    nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjasymmtp = new TH3D("hdjasymmtp",";p_{T,tag};p_{T,3rd}/p_{T,tag};Asymmetry",
                        nx,&x[0],n3,&v3[0],na,&va[0]);
  hdjmpftp = new TH3D("hdjmpftp",";p_{T,tag};p_{T,3rd}/p_{T,tag};MPF",
                      nx,&x[0],n3,&v3[0],na,&va[0]);
}
  hr21 = new TH1D("hr21",";pt2/pt1",50,0.,1.);
  hr31 = new TH1D("hr31",";pt3/pt1",50,0.,1.);
  hr32 = new TH1D("hr32",";pt3/pt2",50,0.,1.);
  pr21 = new TProfile("pr21",";npvgood;pt2/pt1",jp::npvs,jp::pvrange);
  pr31 = new TProfile("pr31",";npvgood;pt3/pt1",jp::npvs,jp::pvrange);
  pr32 = new TProfile("pr32",";npvgood;pt3/pt2",jp::npvs,jp::pvrange);
  px21 = new TProfile("px21",";npvgood;has pt2/pt1",jp::npvs,jp::pvrange);
  px31 = new TProfile("px31",";npvgood;has pt3/pt1",jp::npvs,jp::pvrange);
  px32 = new TProfile("px32",";npvgood;has pt3/pt2",jp::npvs,jp::pvrange);

  hyeta = new TH1D("hyeta","",100,-0.436/4.,0.436/4.);
  hyeta2 = new TH1D("hyeta2","",100,-0.436/4.,0.436/4.);
  hetaphi = new TH2D("hetaphi","",100,-4.799,4.799,jp::nphis,jp::phirange);

  pdpt->Sumw2();
  pmpf->Sumw2();
  pmpf1->Sumw2();
  pmpf2->Sumw2();
  pmpfx->Sumw2();
if (jp::do3dHistos) {
  hdjasymm->Sumw2();
  hdjmpf->Sumw2();
  hdjasymmtp->Sumw2();
  hdjmpftp->Sumw2();
}
  ppt_probepertag = new TProfile("ppt_probepertag","",nx,&x[0]);

  // MC checks
  htrpu = new TH1D("htrpu","",120,0.,60.); // for PU reweighing
  if (this->ismcdir) {
    hpthat = new TH1D("hpthat","",nx,&x[0]);
    hpthatnlo = new TH1D("hpthatnlo","",nx,&x[0]);
    hjtptsf = new TH1D("hjtptsf",";SF;",100,0,3.);
    hjtptsmearNew = new TH1D("hjtptsmearNew",";C_{JER} (smearing);",100,0.,2.);

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
    h2itvsoot = new TH2D("h2itvsoot","",25,-0.5,24.5,jp::npvs,jp::pvrange);

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
    hupt_g = new TH1D("hupt_g","",nx,&x[0]); // unmatched partons vs gen pT
    hpt_g0_tmp = new TH1D("hpt_g0_tmp","",nx,&x[0]);
    ppt_r = new TProfile("ppt_r","",nx,&x[0]);
    ppt_g = new TProfile("ppt_g","",nx,&x[0]);

    // Response closure
    p2rvsnpv = new TProfile2D("p2rvsnpv","",nx,&x[0],jp::npvs,jp::pvrange);
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
    pdy_r = new TProfile2D("pdy_r","",nx,&x[0],jp::nposphis,jp::posphirange);
    pdy_g = new TProfile2D("pdy_g","",nx,&x[0],jp::nposphis,jp::posphirange);
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
