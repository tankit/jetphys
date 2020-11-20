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
  hnjet0          = new TH1D("hnjet0",";Njets;",15,0,15);
  hj0pt0          = new TH1D("hj0pt0",";p_{T,leading} (GeV);",nx,&x[0]);
  hj0phi0         = new TH1D("hj0phi0",";#phi_{leading};",jp::nphis,jp::phirange);
  hj1pt0          = new TH1D("hj1pt0",";p_{T,j2} (GeV);",nx,&x[0]);
  hj2pt0          = new TH1D("hj2pt0",";p_{T,j3} (GeV);",nx,&x[0]);
  hrecoilpt0      = new TH1D("hrecoilpt0",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilphi0     = new TH1D("hrecoilphi0",";#phi_{recoil};",jp::nphis,jp::phirange);
  hrecoilrat0     = new TH1D("hrecoilrat0",";p_{T,j2}/p_{T,recoil};",100,0.,2.);
  hrecoilrat20    = new TH1D("hrecoilrat20",";p_{T,j3}/p_{T,recoil};",100,0.,2.);
  hrecoildphi0    = new TH1D("hrecoildphi0",";#Delta#phi(leading,recoil);",jp::nwposphis,jp::wposphirange);
  hrecoiljets0    = new TH1D("hrecoiljets0",";N_{recoiljets};",15,0,15);
  hrecoilmjb0     = new TH1D("hrecoilmjb0",";p_{T,leading}/p_{T,recoil};",100,0.,4.0);
  hrecoilmjbinv0  = new TH1D("hrecoilmjbinv0",";p_{T,recoil}/p_{T,leading};",100,0.,4.0);
  hrecoilmpftpinv0  = new TH1D("hrecoilmpftpinv0",";1/MPF;",100,0.,4.);
  hrecoilmpftpinv20 = new TH1D("hrecoilmpftpinv20",";1/MPF;",100,0.,4.);
  hrecoilmpf0     = new TH1D("hrecoilmpf0","",100,-100,100);
  hrecoilmpftp0   = new TH1D("hrecoilmpftp0",";MPF(MET1);",100,-1.,3.);
  hrecoilmpf20    = new TH1D("hrecoilmpf20","",100,-100,100);
  hrecoilmpftp20  = new TH1D("hrecoilmpftp20",";MPF(MET2);",100,-1.,3.);
  hleadmpf0       = new TH1D("hleadmpf0","",100,-100,100);
  hleadmpftp0     = new TH1D("hleadmpftp0",";MPF(MET1,leading);",100,-1.,3.);
  hleadmpf20      = new TH1D("hleadmpf20","",100,-100,100);
  hleadmpftp20    = new TH1D("hleadmpftp20",";MPF(MET2,leading);",100,-1.,3.);
  pmjb_recoil0    = new TProfile("pmjb_recoil0",";p_{T,recoil} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjb_leading0   = new TProfile("pmjb_leading0",";p_{T,leading} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_leading0  = new TProfile("pmjbinv_leading0",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpfinv_leading0  = new TProfile("pmpfinv_leading0",";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
  pmpfinv_leading20 = new TProfile("pmpfinv_leading20",";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
  pmpf_recoil0    = new TProfile("pmpf_recoil0",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_recoil20   = new TProfile("pmpf_recoil20",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading0   = new TProfile("pmpf_leading0",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading20  = new TProfile("pmpf_leading20",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpflead_leading0   = new TProfile("pmpflead_leading0",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpflead_leading20  = new TProfile("pmpflead_leading20",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading0   = new TProfile("pmpfleadinv_leading0",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading20  = new TProfile("pmpfleadinv_leading20",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmjbinv_recoil0   = new TProfile("pmjbinv_recoil0",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpflead_recoil0  = new TProfile("pmpflead_recoil0",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  pmpflead_recoil20 = new TProfile("pmpflead_recoil20",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  h2mjb_recoil0       = new TH2D("h2mjb_recoil0",";p_{T,recoil} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjb_leading0      = new TH2D("h2mjb_leading0",";p_{T,leading} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_leading0   = new TH2D("h2mjbinv_leading0",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpfinv_leading0   = new TH2D("h2mpfinv_leading0",";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
  h2mpfinv_leading20  = new TH2D("h2mpfinv_leading20",";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
  h2mpf_recoil0       = new TH2D("h2mpf_recoil0",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_recoil20      = new TH2D("h2mpf_recoil20",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading0      = new TH2D("h2mpf_leading0",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading20     = new TH2D("h2mpf_leading20",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpflead_leading0  = new TH2D("h2mpflead_leading0",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpflead_leading20 = new TH2D("h2mpflead_leading20",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading0  = new TH2D("h2mpfleadinv_leading0",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading20 = new TH2D("h2mpfleadinv_leading20",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mjbinv_recoil0   = new TH2D("h2mjbinv_recoil0",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil0  = new TH2D("h2mpflead_recoil0",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil20 = new TH2D("h2mpflead_recoil20",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  hnjet1          = new TH1D("hnjet1",";Njets;",15,0,15);
  hj0pt1          = new TH1D("hj0pt1",";p_{T,leading} (GeV);",nx,&x[0]);
  hj0phi1         = new TH1D("hj0phi1",";#phi_{leading};",jp::nphis,jp::phirange);
  hj1pt1          = new TH1D("hj1pt1",";p_{T,j2} (GeV);",nx,&x[0]);
  hj2pt1          = new TH1D("hj2pt1",";p_{T,j3} (GeV);",nx,&x[0]);
  hrecoilpt1      = new TH1D("hrecoilpt1",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilphi1     = new TH1D("hrecoilphi1",";#phi_{recoil};",jp::nphis,jp::phirange);
  hrecoilrat1     = new TH1D("hrecoilrat1",";p_{T,j2}/p_{T,recoil};",100,0.,2.);
  hrecoilrat21    = new TH1D("hrecoilrat21",";p_{T,j3}/p_{T,recoil};",100,0.,2.);
  hrecoildphi1    = new TH1D("hrecoildphi1",";#Delta#phi(leading,recoil);",jp::nwposphis,jp::wposphirange);
  hrecoiljets1    = new TH1D("hrecoiljets1",";N_{recoiljets};",15,0,15);
  hrecoilmjb1     = new TH1D("hrecoilmjb1",";p_{T,leading}/p_{T,recoil};",100,0.,4.0);
  hrecoilmjbinv1  = new TH1D("hrecoilmjbinv1",";p_{T,recoil}/p_{T,leading};",100,0.,4.0);
  hrecoilmpftpinv1  = new TH1D("hrecoilmpftpinv1",";1/MPF;",100,0.,4.);
  hrecoilmpftpinv21 = new TH1D("hrecoilmpftpinv21",";1/MPF;",100,0.,4.);
  hrecoilmpf1     = new TH1D("hrecoilmpf1","",100,-100,100);
  hrecoilmpftp1   = new TH1D("hrecoilmpftp1",";MPF(MET1);",100,-1.,3.);
  hrecoilmpf21    = new TH1D("hrecoilmpf21","",100,-100,100);
  hrecoilmpftp21  = new TH1D("hrecoilmpftp21",";MPF(MET2);",100,-1.,3.);
  hleadmpf1       = new TH1D("hleadmpf1","",100,-100,100);
  hleadmpftp1     = new TH1D("hleadmpftp1",";MPF(MET1,leading);",100,-1.,3.);
  hleadmpf21      = new TH1D("hleadmpf21","",100,-100,100);
  hleadmpftp21    = new TH1D("hleadmpftp21",";MPF(MET2,leading);",100,-1.,3.);
  pmjb_recoil1    = new TProfile("pmjb_recoil1",";p_{T,recoil} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjb_leading1   = new TProfile("pmjb_leading1",";p_{T,leading} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_leading1  = new TProfile("pmjbinv_leading1",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpfinv_leading1  = new TProfile("pmpfinv_leading1",";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
  pmpfinv_leading21 = new TProfile("pmpfinv_leading21",";p_{T,leading} (GeV);1 - MPF",nx,&x[0]);
  pmpf_recoil1    = new TProfile("pmpf_recoil1",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_recoil21   = new TProfile("pmpf_recoil21",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading1   = new TProfile("pmpf_leading1",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading21  = new TProfile("pmpf_leading21",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpflead_leading1   = new TProfile("pmpflead_leading1",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpflead_leading21  = new TProfile("pmpflead_leading21",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading1   = new TProfile("pmpfleadinv_leading1",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading21  = new TProfile("pmpfleadinv_leading21",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmjbinv_recoil1   = new TProfile("pmjbinv_recoil1",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpflead_recoil1  = new TProfile("pmpflead_recoil1",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  pmpflead_recoil21 = new TProfile("pmpflead_recoil21",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  h2mjb_recoil1       = new TH2D("h2mjb_recoil1",";p_{T,recoil} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjb_leading1      = new TH2D("h2mjb_leading1",";p_{T,leading} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_leading1   = new TH2D("h2mjbinv_leading1",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpfinv_leading1   = new TH2D("h2mpfinv_leading1",";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
  h2mpfinv_leading21  = new TH2D("h2mpfinv_leading21",";p_{T,leading} (GeV);1 - MPF",nx,&x[0],100,-4,6);
  h2mpf_recoil1       = new TH2D("h2mpf_recoil1",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_recoil21      = new TH2D("h2mpf_recoil21",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading1      = new TH2D("h2mpf_leading1",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading21     = new TH2D("h2mpf_leading21",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpflead_leading1  = new TH2D("h2mpflead_leading1",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpflead_leading21 = new TH2D("h2mpflead_leading21",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading1  = new TH2D("h2mpfleadinv_leading1",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading21 = new TH2D("h2mpfleadinv_leading21",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mjbinv_recoil1   = new TH2D("h2mjbinv_recoil1",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil1  = new TH2D("h2mpflead_recoil1",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil21 = new TH2D("h2mpflead_recoil21",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  hnjet           = new TH1D("hnjet",";Njets;",15,0,15);
  hj0pt           = new TH1D("hj0pt",";p_{T,leading} (GeV);",nx,&x[0]);
  hj0phi          = new TH1D("hj0phi",";#phi_{leading};",jp::nphis,jp::phirange);
  hj1pt           = new TH1D("hj1pt",";p_{T,j2} (GeV);",nx,&x[0]);
  hj2pt           = new TH1D("hj2pt",";p_{T,j3} (GeV);",nx,&x[0]);
  hrecoilpt       = new TH1D("hrecoilpt",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilphi      = new TH1D("hrecoilphi",";#phi_{T,recoil};",jp::nphis,jp::phirange);
  hrecoilrat      = new TH1D("hrecoilrat",";p_{T,j2}/p_{T,recoil};",100,0.,0.6);
  hrecoilrat2     = new TH1D("hrecoilrat2",";p_{T,j3}/p_{T,recoil};",100,0.,0.6);
  hrecoildphi     = new TH1D("hrecoildphi",";#Delta#phi(leading,recoil);",jp::nwposphis,jp::wposphirange);
  hrecoiljets     = new TH1D("hrecoiljets",";N_{recoiljets};",15,0,15);
  hrecoilmjb      = new TH1D("hrecoilmjb",";p_{T,leading}/p_{T,recoil};",100,0.,2.0);
  hrecoilmjb2     = new TH1D("hrecoilmjb2",";p_{T,leading}/p_{T,recoil};",100,0.,4.0);
  hrecoilmjbinv   = new TH1D("hrecoilmjbinv",";p_{T,recoil}/p_{T,leading};",100,0.,4.0);
  hrecoilmpftpinv  = new TH1D("hrecoilmpftpinv",";1/MPF;",100,0.,4.);
  hrecoilmpftpinv2 = new TH1D("hrecoilmpftpinv2",";1/MPF;",100,0.,4.);
  hrecoilmpf      = new TH1D("hrecoilmpf","",100,-100,100);
  hrecoilmpftp    = new TH1D("hrecoilmpftp",";MPF(MET1);",100,-1.,3.);
  hrecoilmpf2     = new TH1D("hrecoilmpf2","",100,-100,100);
  hrecoilmpftp2   = new TH1D("hrecoilmpftp2",";MPF(MET2);",100,-1.,3.);
  hleadmpf        = new TH1D("hleadmpf","",100,-100,100);
  hleadmpftp      = new TH1D("hleadmpftp",";MPF(MET1,leading);",100,-1.,3.);
  hleadmpf2       = new TH1D("hleadmpf2","",100,-100,100);
  hleadmpftp2     = new TH1D("hleadmpftp2",";MPF(MET2,leading);",100,-1.,3.);
  pmjb_recoil     = new TProfile("pmjb_recoil",";p_{T,recoil} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjb_leading    = new TProfile("pmjb_leading",";p_{T,leading} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_leading = new TProfile("pmjbinv_leading",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpfinv_leading  = new TProfile("pmpfinv_leading",";p_{T,leading} (GeV);1/MPF",nx,&x[0]);
  pmpfinv_leading2 = new TProfile("pmpfinv_leading2",";p_{T,leading} (GeV);1/MPF",nx,&x[0]);
  pmpf_recoil     = new TProfile("pmpf_recoil",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_recoil2    = new TProfile("pmpf_recoil2",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading    = new TProfile("pmpf_leading",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_leading2   = new TProfile("pmpf_leading2",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0]);
  pmpflead_leading   = new TProfile("pmpflead_leading",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpflead_leading2  = new TProfile("pmpflead_leading2",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading   = new TProfile("pmpfleadinv_leading",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmpfleadinv_leading2  = new TProfile("pmpfleadinv_leading2",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0]);
  pmjbinv_recoil   = new TProfile("pmjbinv_recoil",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0]);
  pmpflead_recoil  = new TProfile("pmpflead_recoil",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  pmpflead_recoil2 = new TProfile("pmpflead_recoil2",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0]);
  h2mjb_recoil       = new TH2D("h2mjb_recoil",";p_{T,recoil} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjb_leading      = new TH2D("h2mjb_leading",";p_{T,leading} (GeV);p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_leading   = new TH2D("h2mjbinv_leading",";p_{T,leading} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpfinv_leading   = new TH2D("h2mpfinv_leading",";p_{T,leading} (GeV);1/R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_leading2  = new TH2D("h2mpfinv_leading2",";p_{T,leading} (GeV);1/R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_recoil       = new TH2D("h2mpf_recoil",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_recoil2      = new TH2D("h2mpf_recoil2",";p_{T,recoil} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading      = new TH2D("h2mpf_leading",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_leading2     = new TH2D("h2mpf_leading2",";p_{T,leading} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpflead_leading  = new TH2D("h2mpflead_leading",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpflead_leading2 = new TH2D("h2mpflead_leading2",";p_{T,leading} (GeV);1 - MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading  = new TH2D("h2mpfleadinv_leading",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mpfleadinv_leading2 = new TH2D("h2mpfleadinv_leading2",";p_{T,leading} (GeV);1 + MPF_{leading}",nx,&x[0],100,-4,6);
  h2mjbinv_recoil   = new TH2D("h2mjbinv_recoil",";p_{T,recoil} (GeV);1/R_{MJB} = p_{T,recoil}/p_{T,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil  = new TH2D("h2mpflead_recoil",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  h2mpflead_recoil2 = new TH2D("h2mpflead_recoil2",";p_{T,recoil} (GeV);R_{MPF,leading}",nx,&x[0],100,-4,6);
  hdphi_recoiljet    = new TH1D("hdphi_recoiljet",";#Delta#phi(leading,recoiljet);",jp::nwposphis,jp::wposphirange);
  hdphi_recoiljetsel = new TH1D("hdphi_recoiljetsel",";#Delta#phi(leading,recoiljets);",50,1,TMath::Pi());
  hweight_recoiljetsel  = new TH1D("hweight_recoiljetsel",";p_{T,recoil} #times cos(#Delta#phi(leading,recoiljet);",100,0,10);
  h2etapt_recoiljetsel  = new TH2D("h2etapt_recoiljetsel",";#eta_{recoiljet};p_{T,recoiljet};",100,-5,5,nx,&x[0]);
  h2etapt_recoiljetsel_weight  = new TH2D("h2etapt_recoiljetsel_weight",";#eta_{recoiljet};p_{T,recoiljet};",100,-5,5,nx,&x[0]);
  hj0genpt_pre       = new TH1D("hj0genpt_pre",";p_{T,leading,gen} (GeV);",nx,&x[0]);
  hj0genr_pre        = new TH1D("hj0genr_pre",";;",100,0.0,0.5);
  hj0genpt0          = new TH1D("hj0genpt0",";p_{T,leading,gen} (GeV);",nx,&x[0]);
  hj0genr0           = new TH1D("hj0genr0",";;",100,0.0,0.5);
  hj0genpt1          = new TH1D("hj0genpt1",";p_{T,leading,gen} (GeV);",nx,&x[0]);
  hj0genr1           = new TH1D("hj0genr1",";;",100,0.0,0.5);
  hj0genpt           = new TH1D("hj0genpt",";p_{T,leading,gen} (GeV);",nx,&x[0]);
  hj0genr            = new TH1D("hj0genr",";;",100,0.0,0.5);
  hrecoilgenpt0      = new TH1D("hrecoilgenpt0",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilgenphi0     = new TH1D("hrecoilgenphi0",";#phi_{recoil};",jp::nphis,jp::phirange);
  hrecoilgenpt1      = new TH1D("hrecoilgenpt1",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilgenphi1     = new TH1D("hrecoilgenphi1",";#phi_{recoil};",jp::nphis,jp::phirange);
  hrecoilgenpt       = new TH1D("hrecoilgenpt",";p_{T,recoil} (GeV);",nx,&x[0]);
  hrecoilgenphi      = new TH1D("hrecoilgenphi",";#phi_{recoil};",jp::nphis,jp::phirange);
  hresp_jet          = new TH1D("hresp_jet",";Resp(jet,gen);",100,0,2.5);
  hresp_jetsel       = new TH1D("hresp_jetsel",";Resp(jetsel,gen);",100,0,2.5);
  hresp_tmp          = new TH1D("hresp_tmp",";Resp(jet,gen);",100,0,2.5);
  hresp_leading0     = new TH1D("hresp_leading0",";Resp(leading,gen);",100,0,2.5);
  hresp_recoil0      = new TH1D("hresp_recoil0",";Resp(recoil,gen);",100,0,5);
  hresp_mjb0         = new TH1D("hresp_mjb0",";Resp(mjb,gen);",100,0,2.5);
  hresp_leading1     = new TH1D("hresp_leading1",";Resp(leading,gen);",100,0,2.5);
  hresp_recoil1      = new TH1D("hresp_recoil1",";Resp(recoil,gen);",100,0,5);
  hresp_mjb1         = new TH1D("hresp_mjb1",";Resp(mjb,gen);",100,0,2.5);
  hresp_leading      = new TH1D("hresp_leading",";Resp(leading,gen);",100,0,2.5);
  hresp_recoil       = new TH1D("hresp_recoil",";Resp(recoil,gen);",100,0,5);
  hresp_mjb          = new TH1D("hresp_mjb",";Resp(mjb,gen);",100,0,2.5);
  presp_leading0     = new TProfile("presp_leading0","",nx,&x[0]);
  presp_recoil0      = new TProfile("presp_recoil0","",nx,&x[0]);
  presp_leading1     = new TProfile("presp_leading1","",nx,&x[0]);
  presp_recoil1      = new TProfile("presp_recoil1","",nx,&x[0]);
  presp_leading      = new TProfile("presp_leading","",nx,&x[0]);
  presp_recoil       = new TProfile("presp_recoil","",nx,&x[0]);
  hdphi_jet120       = new TH1D("hdphi_jet120",";#Delta#phi(j1,j2);",jp::nwposphis,jp::wposphirange);
  hdphi_jet230       = new TH1D("hdphi_jet230",";#Delta#phi(j2,j3);",jp::nwposphis,jp::wposphirange);
  hdphi_jet121       = new TH1D("hdphi_jet121",";#Delta#phi(j1,j2);",jp::nwposphis,jp::wposphirange);
  hdphi_jet231       = new TH1D("hdphi_jet231",";#Delta#phi(j2,j3);",jp::nwposphis,jp::wposphirange);
  hdphi_jet12        = new TH1D("hdphi_jet12",";#Delta#phi(j1,j2);",jp::nwposphis,jp::wposphirange);
  hdphi_jet23        = new TH1D("hdphi_jet23",";#Delta#phi(j2,j3);",jp::nwposphis,jp::wposphirange);
  hdphi_recoilmin    = new TH1D("hdphi_recoilmin",";Min of #Delta#phi(leading,recoiljets);",jp::nwposphis,jp::wposphirange);
  hdphi_recoilmax    = new TH1D("hdphi_recoilmax",";Max of #Delta#phi(leading,recoiljets);",jp::nwposphis,jp::wposphirange);
  hdphi_recoilcol    = new TH2D("hdphi_recoilcol",";#Delta#phi_{min};#Delta#phi_{max}",jp::nwposphis,jp::wposphirange,jp::nwposphis,jp::wposphirange);
  hCRecoil           = new TH1D("hCRecoil",";CRecoil;",100,0,1);
  hCRecoil_ptcl      = new TH1D("hCRecoil_ptcl",";CRecoil_ptcl;",100,0,1);
  pCRecoil           = new TProfile("pCRecoil",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_leading   = new TProfile("pCRecoil_leading",";p_{T,leading} (GeV);",nx,&x[0]);
  pCRecoil_ptave     = new TProfile("pCRecoil_ptave",";p_{T,ave} (GeV);",nx,&x[0]);
  pCRecoil_ptcl      = new TProfile("pCRecoil_ptcl",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2     = new TProfile("pCRecoil_ptcl2",";p_{T,recoil_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_leading = new TProfile("pCRecoil_ptcl2_leading",";p_{T,leading_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_ptave   = new TProfile("pCRecoil_ptcl2_ptave",";p_{T,ave_gen} (GeV);",nx,&x[0]);
  pjtgenptsf         = new TProfile("pjtgenptsf",";p_{T,leading,gen};JER SF",nx,&x[0]);
  pjtptsf            = new TProfile("pjtptsf",";p_{T,leading,reco};JER SF",nx,&x[0]);
  pjtgenptsf_recoil  = new TProfile("pjtgenptsf_recoil",";p_{T,gen} of jets in the recoil;JER SF",nx,&x[0]);
  pjtptsf_recoil     = new TProfile("pjtptsf_recoil",";p_{T,reco} of jets in the recoil;JER SF",nx,&x[0]);
  hCRecoil1          = new TH1D("hCRecoil1",";CRecoil;",100,0,1);
  hCRecoil_ptcl1     = new TH1D("hCRecoil_ptcl1",";CRecoil_ptcl;",100,0,1);
  pCRecoil1          = new TProfile("pCRecoil1",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_leading1  = new TProfile("pCRecoil_leading1",";p_{T,leading} (GeV);",nx,&x[0]);
  pCRecoil_ptave1    = new TProfile("pCRecoil_ptave1",";p_{T,ave} (GeV);",nx,&x[0]);
  pCRecoil_ptcl1     = new TProfile("pCRecoil_ptcl1",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_ptcl21    = new TProfile("pCRecoil_ptcl21",";p_{T,recoil_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_leading1 = new TProfile("pCRecoil_ptcl2_leading1",";p_{T,leading_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_ptave1   = new TProfile("pCRecoil_ptcl2_ptave1",";p_{T,ave_gen} (GeV);",nx,&x[0]);
  pjtgenptsf1         = new TProfile("pjtgenptsf1",";p_{T,leading,gen};JER SF",nx,&x[0]);
  pjtptsf1            = new TProfile("pjtptsf1",";p_{T,leading,reco};JER SF",nx,&x[0]);
  hCRecoil0          = new TH1D("hCRecoil0",";CRecoil;",100,0,1);
  hCRecoil_ptcl0     = new TH1D("hCRecoil_ptcl0",";CRecoil_ptcl;",100,0,1);
  pCRecoil0          = new TProfile("pCRecoil0",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_leading0  = new TProfile("pCRecoil_leading0",";p_{T,leading} (GeV);",nx,&x[0]);
  pCRecoil_ptave0    = new TProfile("pCRecoil_ptave0",";p_{T,ave} (GeV);",nx,&x[0]);
  pCRecoil_ptcl0     = new TProfile("pCRecoil_ptcl0",";p_{T,recoil} (GeV);",nx,&x[0]);
  pCRecoil_ptcl20    = new TProfile("pCRecoil_ptcl20",";p_{T,recoil_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_leading0 = new TProfile("pCRecoil_ptcl2_leading0",";p_{T,leading_gen} (GeV);",nx,&x[0]);
  pCRecoil_ptcl2_ptave0   = new TProfile("pCRecoil_ptcl2_ptave0",";p_{T,ave_gen} (GeV);",nx,&x[0]);
  pjtgenptsf0         = new TProfile("pjtgenptsf0",";p_{T,leading,gen};JER SF",nx,&x[0]);
  pjtptsf0            = new TProfile("pjtptsf0",";p_{T,leading,reco};JER SF",nx,&x[0]);

  pmjb_ptave0       = new TProfile("pmjb_ptave0",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_ptave0    = new TProfile("pmjbinv_ptave0",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmpf_ptave0       = new TProfile("pmpf_ptave0",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_ptave20      = new TProfile("pmpf_ptave20",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave0    = new TProfile("pmpfinv_ptave0",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave20   = new TProfile("pmpfinv_ptave20",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmjb_ptave1       = new TProfile("pmjb_ptave1",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_ptave1    = new TProfile("pmjbinv_ptave1",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmpf_ptave1       = new TProfile("pmpf_ptave1",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_ptave21      = new TProfile("pmpf_ptave21",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave1    = new TProfile("pmpfinv_ptave1",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave21   = new TProfile("pmpfinv_ptave21",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmjb_ptave        = new TProfile("pmjb_ptave",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmjbinv_ptave     = new TProfile("pmjbinv_ptave",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0]);
  pmpf_ptave        = new TProfile("pmpf_ptave",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpf_ptave2       = new TProfile("pmpf_ptave2",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave     = new TProfile("pmpfinv_ptave",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);
  pmpfinv_ptave2    = new TProfile("pmpfinv_ptave2",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0]);

  pmpf_ptave_ue0    = new TProfile("pmpf_ptave_ue0",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpf_ptave_one0   = new TProfile("pmpf_ptave_one0",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpf_ptave_n0     = new TProfile("pmpf_ptave_n0",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);
  pmpf_ptave_ue1    = new TProfile("pmpf_ptave_ue1",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpf_ptave_one1   = new TProfile("pmpf_ptave_one1",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpf_ptave_n1     = new TProfile("pmpf_ptave_n1",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);
  pmpf_ptave_ue     = new TProfile("pmpf_ptave_ue",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpf_ptave_one    = new TProfile("pmpf_ptave_one",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpf_ptave_n      = new TProfile("pmpf_ptave_n",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);

  pmpfinv_ptave_ue0  = new TProfile("pmpfinv_ptave_ue0",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpfinv_ptave_one0 = new TProfile("pmpfinv_ptave_one0",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpfinv_ptave_n0   = new TProfile("pmpfinv_ptave_n0",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);
  pmpfinv_ptave_ue1  = new TProfile("pmpfinv_ptave_ue1",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpfinv_ptave_one1 = new TProfile("pmpfinv_ptave_one1",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpfinv_ptave_n1   = new TProfile("pmpfinv_ptave_n1",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);
  pmpfinv_ptave_ue   = new TProfile("pmpfinv_ptave_ue",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0]);
  pmpfinv_ptave_one  = new TProfile("pmpfinv_ptave_one",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0]);
  pmpfinv_ptave_n    = new TProfile("pmpfinv_ptave_n",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0]);

  h2mjb_ptave0       = new TH2D("h2mjb_ptave0",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_ptave0    = new TH2D("h2mjbinv_ptave0",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mpf_ptave0       = new TH2D("h2mpf_ptave0",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_ptave20      = new TH2D("h2mpf_ptave20",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave0    = new TH2D("h2mpfinv_ptave0",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave20   = new TH2D("h2mpfinv_ptave20",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mjb_ptave1       = new TH2D("h2mjb_ptave1",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_ptave1    = new TH2D("h2mjbinv_ptave1",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mpf_ptave1       = new TH2D("h2mpf_ptave1",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_ptave21      = new TH2D("h2mpf_ptave21",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave1    = new TH2D("h2mpfinv_ptave1",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave21   = new TH2D("h2mpfinv_ptave21",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mjb_ptave        = new TH2D("h2mjb_ptave",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mjbinv_ptave     = new TH2D("h2mjbinv_ptave",";p_{T,ave} (GeV);R_{MJB} = p_{T,leading}/p_{T,recoil}",nx,&x[0],100,-4,6);
  h2mpf_ptave        = new TH2D("h2mpf_ptave",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpf_ptave2       = new TH2D("h2mpf_ptave2",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave     = new TH2D("h2mpfinv_ptave",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave2    = new TH2D("h2mpfinv_ptave2",";p_{T,ave} (GeV);R_{MPF}",nx,&x[0],100,-4,6);

  h2mpf_ptave_ue0    = new TH2D("h2mpf_ptave_ue0",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpf_ptave_one0   = new TH2D("h2mpf_ptave_one0",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpf_ptave_n0     = new TH2D("h2mpf_ptave_n0",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);
  h2mpf_ptave_ue1    = new TH2D("h2mpf_ptave_ue1",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpf_ptave_one1   = new TH2D("h2mpf_ptave_one1",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpf_ptave_n1     = new TH2D("h2mpf_ptave_n1",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);
  h2mpf_ptave_ue     = new TH2D("h2mpf_ptave_ue",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpf_ptave_one    = new TH2D("h2mpf_ptave_one",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpf_ptave_n      = new TH2D("h2mpf_ptave_n",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);

  h2mpfinv_ptave_ue0  = new TH2D("h2mpfinv_ptave_ue0",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_one0 = new TH2D("h2mpfinv_ptave_one0",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_n0   = new TH2D("h2mpfinv_ptave_n0",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_ue1  = new TH2D("h2mpfinv_ptave_ue1",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_one1 = new TH2D("h2mpfinv_ptave_one1",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_n1   = new TH2D("h2mpfinv_ptave_n1",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_ue   = new TH2D("h2mpfinv_ptave_ue",";p_{T,ave} (GeV);R_{MPF_u}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_one  = new TH2D("h2mpfinv_ptave_one",";p_{T,ave} (GeV);R_{MPF_1}",nx,&x[0],100,-4,6);
  h2mpfinv_ptave_n    = new TH2D("h2mpfinv_ptave_n",";p_{T,ave} (GeV);R_{MPF_n}",nx,&x[0],100,-4,6);

  hleadmpf_ue0          = new TH1D("hleadmpf_ue0","",100,-100,100);
  hleadmpf_ue1          = new TH1D("hleadmpf_ue1","",100,-100,100);
  hleadmpf_ue           = new TH1D("hleadmpf_ue","",100,-100,100);
  hleadmpf_one0         = new TH1D("hleadmpf_one0","",100,-100,100);
  hleadmpf_one1         = new TH1D("hleadmpf_one1","",100,-100,100);
  hleadmpf_one          = new TH1D("hleadmpf_one","",100,-100,100);
  hleadmpf_n0           = new TH1D("hleadmpf_n0","",100,-100,100);
  hleadmpf_n1           = new TH1D("hleadmpf_n1","",100,-100,100);
  hleadmpf_n            = new TH1D("hleadmpf_n","",100,-100,100);

  //hleadmpftp_ue0        = new TH1D("hleadmpftp_ue0",";MPF(ue,leading);",100,-3.,3.);
  //hleadmpftp_ue1        = new TH1D("hleadmpftp_ue1",";MPF(ue,leading);",100,-3.,3.);
  //hleadmpftp_ue         = new TH1D("hleadmpftp_ue",";MPF(ue,leading);",100,-3.,3.);
  //hleadmpftp_one0       = new TH1D("hleadmpftp_one0",";MPF(leading,leading);",100,-3.,3.);
  //hleadmpftp_one1       = new TH1D("hleadmpftp_one1",";MPF(leading,leading);",100,-3.,3.);
  //hleadmpftp_one        = new TH1D("hleadmpftp_one",";MPF(leading,leading);",100,-3.,3.);
  //hleadmpftp_n0         = new TH1D("hleadmpftp_n0",";MPF(recoil,leading);",100,-3.,3.);
  //hleadmpftp_n1         = new TH1D("hleadmpftp_n1",";MPF(recoil,leading);",100,-3.,3.);
  //hleadmpftp_n          = new TH1D("hleadmpftp_n",";MPF(recoil,leading);",100,-3.,3.);

  pmpflead_leading_ue0  = new TProfile("pmpflead_leading_ue0",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0]);
  pmpflead_leading_one0 = new TProfile("pmpflead_leading_one0",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0]);
  pmpflead_leading_n0   = new TProfile("pmpflead_leading_n0",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0]);
  pmpflead_leading_ue1  = new TProfile("pmpflead_leading_ue1",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0]);
  pmpflead_leading_one1 = new TProfile("pmpflead_leading_one1",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0]);
  pmpflead_leading_n1   = new TProfile("pmpflead_leading_n1",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0]);
  pmpflead_leading_ue   = new TProfile("pmpflead_leading_ue",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0]);
  pmpflead_leading_one  = new TProfile("pmpflead_leading_one",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0]);
  pmpflead_leading_n    = new TProfile("pmpflead_leading_n",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0]);

  h2mpflead_leading_ue0  = new TH2D("h2mpflead_leading_ue0",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0],100,-5,5);
  h2mpflead_leading_one0 = new TH2D("h2mpflead_leading_one0",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0],100,-5,5);
  h2mpflead_leading_n0   = new TH2D("h2mpflead_leading_n0",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0],100,-5,5);
  h2mpflead_leading_ue1  = new TH2D("h2mpflead_leading_ue1",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0],100,-5,5);
  h2mpflead_leading_one1 = new TH2D("h2mpflead_leading_one1",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0],100,-5,5);
  h2mpflead_leading_n1   = new TH2D("h2mpflead_leading_n1",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0],100,-5,5);
  h2mpflead_leading_ue   = new TH2D("h2mpflead_leading_ue",";p_{T,leading} (GeV);MPF_{ue}",nx,&x[0],100,-5,5);
  h2mpflead_leading_one  = new TH2D("h2mpflead_leading_one",";p_{T,leading} (GeV);1 + MPF_{1}",nx,&x[0],100,-5,5);
  h2mpflead_leading_n    = new TH2D("h2mpflead_leading_n",";p_{T,leading} (GeV);MPF_{n}",nx,&x[0],100,-5,5);

  prho_recoil        = new TProfile("prho_recoil",";p_{T,recoil} (GeV);#rho",nx,&x[0]);
  prho_leading       = new TProfile("prho_leading",";p_{T,leading} (GeV);#rho",nx,&x[0]);
  prho_ptave         = new TProfile("prho_ptave",";p_{T,ave} (GeV);#rho",nx,&x[0]);
  prho_recoil0       = new TProfile("prho_recoil0",";p_{T,recoil} (GeV);#rho",nx,&x[0]);
  prho_leading0      = new TProfile("prho_leading0",";p_{T,leading} (GeV);#rho",nx,&x[0]);
  prho_ptave0        = new TProfile("prho_ptave0",";p_{T,ave} (GeV);#rho",nx,&x[0]);
  prho_recoil1       = new TProfile("prho_recoil1",";p_{T,recoil} (GeV);#rho",nx,&x[0]);
  prho_leading1      = new TProfile("prho_leading1",";p_{T,leading} (GeV);#rho",nx,&x[0]);
  prho_ptave1        = new TProfile("prho_ptave1",";p_{T,ave} (GeV);#rho",nx,&x[0]);
  prho_recoil2       = new TProfile("prho_recoil2",";p_{T,recoil} (GeV);#rho",nx,&x[0]);
  prho_leading2      = new TProfile("prho_leading2",";p_{T,leading} (GeV);#rho",nx,&x[0]);
  prho_ptave2        = new TProfile("prho_ptave2",";p_{T,ave} (GeV);#rho",nx,&x[0]);

  h2rho_recoil       = new TH2D("h2rho_recoil",";p_{T,recoil} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_leading      = new TH2D("h2rho_leading",";p_{T,leading} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_ptave        = new TH2D("h2rho_ptave",";p_{T,ave} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_recoil0      = new TH2D("h2rho_recoil0",";p_{T,recoil} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_leading0     = new TH2D("h2rho_leading0",";p_{T,leading} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_ptave0       = new TH2D("h2rho_ptave0",";p_{T,ave} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_recoil1      = new TH2D("h2rho_recoil1",";p_{T,recoil} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_leading1     = new TH2D("h2rho_leading1",";p_{T,leading} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_ptave1       = new TH2D("h2rho_ptave1",";p_{T,ave} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_recoil2      = new TH2D("h2rho_recoil2",";p_{T,recoil} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_leading2     = new TH2D("h2rho_leading2",";p_{T,leading} (GeV);#rho",nx,&x[0],200,0,80);
  h2rho_ptave2       = new TH2D("h2rho_ptave2",";p_{T,ave} (GeV);#rho",nx,&x[0],200,0,80);

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
