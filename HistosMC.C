// Purpose:  Create eta-dependent histograms for inclusive jets analysis
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#include "HistosMC.h"

HistosMC::HistosMC(TDirectory *dir, string trigname) {

  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // phase space
  this->trigname = trigname;

  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);

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
                                  jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],na,&va[0]) );
    hdjasymmtp.push_back(new TH3D((string("hdjasymmtp_a")+number).c_str(),";p_{T,tag};#eta;Asymmetry",
                                  jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],na,&va[0]) );
    hdjresp_tag.push_back(new TH3D((string("hdjresp_tag_a")+number).c_str(),";p_{T,ave};#eta;R_{tag}",
                                   jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],nb,&vb[0]) );
    hdjresptp_tag.push_back(new TH3D((string("hdjresptp_tag_a")+number).c_str(),";p_{T,ave};#eta;R_{tag}",
                                      jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],nb,&vb[0]) );
    hdjresp_probe.push_back(new TH3D((string("hdjresp_probe_a")+number).c_str(),";p_{T,ave};#eta;R_{probe}",
                                     jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],nb,&vb[0]) );
    hdjresptp_probe.push_back(new TH3D((string("hdjresptp_probe_a")+number).c_str(),";p_{T,ave};#eta;R_{probe}",
                                       jp::npts,&jp::ptrange[0],jp::nwetas,&jp::wetarange[0],nb,&vb[0]) );
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

void HistosMC::Write() {
  dir->cd();
  dir->Write();
}


HistosMC::~HistosMC() {
  Write();
};
