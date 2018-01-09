// Purpose:  Create eta-dependent histograms for inclusive jets analysis
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#include "mcHistos.h"

mcHistos::mcHistos(TDirectory *dir, string trigname) {

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
    {1,    5,    6,    8,    10,   12,   15,   18,   21,   24,   28,   32,
     37,   43,   49,   56,   64,   74,   84,   97,   114,  133,  153,  174,  
     196,  220,  245,  272,  300,  330,  362,  395,  430,  468,  507,  548,  
     592,  638,  686,  737,  790,  846,  905,  967,  1032, 1101, 1172, 1248, 
     1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 
     2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 
     4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int npts = sizeof(ptrange)/sizeof(ptrange[0])-1;
  
  const double etarange[] =
  {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  const unsigned int netas = sizeof(etarange)/sizeof(etarange[0])-1;
  
  const int na = 200;
  vector<double> va(na+1);
  for (unsigned int i = 0; i != na+1; ++i)
    va[i] = -1. + 2.*i/na;
  
  const int nb = 200;
  vector<double> vb(nb+1);
  for (unsigned int i = 0; i != nb+1; ++i)
    vb[i] = 0+2.*i/nb;

  // Loop over alpha entries of interest
  for (auto alpha : alpharange) {
    int major_alpha = 100*alpha;
    // Start by coming up with a nice number identifier
    int padding = TMath::Log10(major_alpha);
    padding = 2-padding;
    while (major_alpha%10==0)
      major_alpha /= 10; 
    string number = std::to_string(major_alpha);
    for (int i = 0; i < padding; ++i)
      number = string("0")+number;

    // Fill all histo types with a corresponding histogram
    hdjasymm.push_back(  new TH3D((string("hdjasymm_a")+number).c_str(),";p_{T,ave};#eta;Asymmetry",
                                  npts,&ptrange[0],netas,&etarange[0],na,&va[0]) );
    hdjasymmtp.push_back(new TH3D((string("hdjasymmtp_a")+number).c_str(),";p_{T,tag};#eta;Asymmetry",
                                  npts,&ptrange[0],netas,&etarange[0],na,&va[0]) );
    hdjresp_tag.push_back(new TH3D((string("hdjresp_tag_a")+number).c_str(),";p_{T,ave};#eta;R_{tag}",
                                   npts,&ptrange[0],netas,&etarange[0],nb,&vb[0]) );
    hdjresptp_tag.push_back(new TH3D((string("hdjresptp_tag_a")+number).c_str(),";p_{T,ave};#eta;R_{tag}",
                                      npts,&ptrange[0],netas,&etarange[0],nb,&vb[0]) );
    hdjresp_probe.push_back(new TH3D((string("hdjresp_probe_a")+number).c_str(),";p_{T,ave};#eta;R_{probe}",
                                     npts,&ptrange[0],netas,&etarange[0],nb,&vb[0]) );
    hdjresptp_probe.push_back(new TH3D((string("hdjresptp_probe_a")+number).c_str(),";p_{T,ave};#eta;R_{probe}",
                                       npts,&ptrange[0],netas,&etarange[0],nb,&vb[0]) );
  }

  // Weights:
  for (unsigned i = 0; i < alpharange.size(); ++i) {
    hdjasymm[i]->Sumw2();
    hdjasymmtp[i]->Sumw2();
    hdjresp_tag[i]->Sumw2();
    hdjresptp_tag[i]->Sumw2();
    hdjresp_probe[i]->Sumw2();
    hdjresptp_probe[i]->Sumw2();
  }
  
  curdir->cd();
}

mcHistos::~mcHistos() {
  dir->cd();
  dir->Write();
  for (unsigned i = 0; i < alpharange.size(); ++i) {
    delete hdjasymm[i];
    delete hdjasymmtp[i];
    delete hdjresp_tag[i];
    delete hdjresptp_tag[i];
    delete hdjresp_probe[i];
    delete hdjresptp_probe[i];
  }
  delete dir;
};
