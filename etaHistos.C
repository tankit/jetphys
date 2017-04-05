// Purpose:  Create eta-dependent histograms for inclusive jets analysis
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#include "etaHistos.h"

etaHistos::etaHistos(TDirectory *dir, string trigname) {

  TDirectory *curdir = gDirectory;
  assert(dir->cd());
  this->dir = dir;
  
  // phase space
  this->trigname = trigname;
  
  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);
  
  // Binning agreed within JTF: pT>100 GeV from CaloJet resolutions,
  // pT<100 GeV to optimize bin widths for PFJets and b-tagging
  // (little higher than resolution, but fairly flat relative width)
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/QCDAnalysis/HighPtJetAnalysis/interface/DefaultPtBins.h?revision=1.2&view=markup
  const double ptrange[] =
    {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
     4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int npts = sizeof(ptrange)/sizeof(ptrange[0])-1;
  
  const double etarange[] =
  {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  const unsigned int netas = sizeof(etarange)/sizeof(etarange[0])-1;
  
  const int na = 200;
  vector<double> va(na+1);
  for (unsigned int i = 0; i != na+1; ++i)
    va[i] = -1. + 2.*i/na;

  hdjasymm_a005 = new TH3D("hdjasymm_a005",";p_{T,ave};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymm_a01 =  new TH3D("hdjasymm_a01",";p_{T,ave};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymm_a015 = new TH3D("hdjasymm_a015",";p_{T,ave};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymm_a02 =  new TH3D("hdjasymm_a02",";p_{T,ave};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymm_a025 = new TH3D("hdjasymm_a025",";p_{T,ave};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymm_a03 =  new TH3D("hdjasymm_a03",";p_{T,ave};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);

  hdjasymmtp_a005 = new TH3D("hdjasymmtp_a005",";p_{T,tag};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymmtp_a01 =  new TH3D("hdjasymmtp_a01",";p_{T,tag};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymmtp_a015 = new TH3D("hdjasymmtp_a015",";p_{T,tag};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymmtp_a02 =  new TH3D("hdjasymmtp_a02",";p_{T,tag};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymmtp_a025 = new TH3D("hdjasymmtp_a025",";p_{T,tag};#eta;Asymmetry",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjasymmtp_a03 =  new TH3D("hdjasymmtp_a03",";p_{T,tag};#eta;Asymmetry", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);

  hdjmpf_a005 = new TH3D("hdjmpf_a005",";p_{T,ave};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpf_a01 =  new TH3D("hdjmpf_a01",";p_{T,ave};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpf_a015 = new TH3D("hdjmpf_a015",";p_{T,ave};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpf_a02 =  new TH3D("hdjmpf_a02",";p_{T,ave};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpf_a025 = new TH3D("hdjmpf_a025",";p_{T,ave};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpf_a03 =  new TH3D("hdjmpf_a03",";p_{T,ave};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);

  hdjmpftp_a005 = new TH3D("hdjmpftp_a005",";p_{T,tag};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpftp_a01 =  new TH3D("hdjmpftp_a01",";p_{T,tag};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpftp_a015 = new TH3D("hdjmpftp_a015",";p_{T,tag};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpftp_a02 =  new TH3D("hdjmpftp_a02",";p_{T,tag};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpftp_a025 = new TH3D("hdjmpftp_a025",";p_{T,tag};#eta;MPF",npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  hdjmpftp_a03 =  new TH3D("hdjmpftp_a03",";p_{T,tag};#eta;MPF", npts,&ptrange[0],netas,&etarange[0],na,&va[0]);
  
  // Weights:
  hdjasymm_a005->Sumw2();
  hdjasymm_a01->Sumw2();
  hdjasymm_a015->Sumw2();
  hdjasymm_a02->Sumw2();
  hdjasymm_a025->Sumw2();
  hdjasymm_a03->Sumw2();
  
  hdjasymmtp_a005->Sumw2();
  hdjasymmtp_a01->Sumw2();
  hdjasymmtp_a015->Sumw2();
  hdjasymmtp_a02->Sumw2();
  hdjasymmtp_a025->Sumw2();
  hdjasymmtp_a03->Sumw2();
  
  hdjmpf_a005->Sumw2();
  hdjmpf_a01->Sumw2();
  hdjmpf_a015->Sumw2();
  hdjmpf_a02->Sumw2();
  hdjmpf_a025->Sumw2();
  hdjmpf_a03->Sumw2();
  
  hdjmpftp_a005->Sumw2();
  hdjmpftp_a01->Sumw2();
  hdjmpftp_a015->Sumw2();
  hdjmpftp_a02->Sumw2();
  hdjmpftp_a025->Sumw2();
  hdjmpftp_a03->Sumw2();
  
  curdir->cd();
}

etaHistos::~etaHistos() {
  dir->cd();
  dir->Write();
  //delete dir;
};
