// Purpose:  Create eta-dependent histograms for inclusive jets analysis
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#include "histosEta.h"

histosEta::histosEta(TDirectory *dir, string trigname, double pttrg, double ptmin, double ptmax) {

  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // phase space
  this->trigname = trigname;

  this->pttrg = pttrg;
  this->ptmin = ptmin;
  this->ptmax = ptmax;

  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);

  const int na = 200;
  vector<double> va(na+1);
  for (unsigned int i = 0; i != na+1; ++i)
    va[i] = -1. + 2.*i/na;

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
    hdjasymm.push_back(  new TH3D((string("hdjasymm_a")+number).c_str(),
                              ";p_{T,ave};#eta;Asymmetry",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
    hdjasymmtp.push_back(new TH3D((string("hdjasymmtp_a")+number).c_str(),
                              ";p_{T,tag};#eta;Asymmetry",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
    hdjasymmpt.push_back(new TH3D((string("hdjasymmpt_a")+number).c_str(),
                              ";p_{T,probe};#eta;Asymmetry",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
    hdjmpf.push_back(    new TH3D((string("hdjmpf_a")+number).c_str(),
                              ";p_{T,ave};#eta;MPF",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
    hdjmpftp.push_back(  new TH3D((string("hdjmpftp_a")+number).c_str(),
                              ";p_{T,tag};#eta;MPF",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
    hdjmpfpt.push_back(  new TH3D((string("hdjmpfpt_a")+number).c_str(),
                              ";p_{T,probe};#eta;MPF",
                              _jp_npts,&_jp_ptrange[0],_jp_nwetas,&_jp_wetarange[0],na,&va[0]) );
  }

  // Weights:
  for (unsigned i = 0; i < alpharange.size(); ++i) {
    hdjasymm[i]->Sumw2();
    hdjasymmtp[i]->Sumw2();
    hdjasymmpt[i]->Sumw2();
    hdjmpf[i]->Sumw2();
    hdjmpftp[i]->Sumw2();
    hdjmpfpt[i]->Sumw2();
  }

  // components vs eta (Ozlem)
  pncandtp_vseta = new TProfile("pcandtp_vseta","",_jp_netas,_jp_etarange);
  pnchtp_vseta = new TProfile("pnchtp_vseta","",_jp_netas,_jp_etarange);
  pnnetp_vseta = new TProfile("pnnetp_vseta","",_jp_netas,_jp_etarange);
  pnnhtp_vseta = new TProfile("pnnhtp_vseta","",_jp_netas,_jp_etarange);
  pncetp_vseta = new TProfile("pncetp_vseta","",_jp_netas,_jp_etarange);
  pnmutp_vseta = new TProfile("pnmutp_vseta","",_jp_netas,_jp_etarange);
  pnhhtp_vseta = new TProfile("pnhhtp_vseta","",_jp_netas,_jp_etarange);
  pnhetp_vseta = new TProfile("pnhetp_vseta","",_jp_netas,_jp_etarange);
  phhftp_vseta = new TProfile("phhftp_vseta","",_jp_netas,_jp_etarange);
  pheftp_vseta = new TProfile("pheftp_vseta","",_jp_netas,_jp_etarange);
  pchftp_vseta = new TProfile("pchftp_vseta","",_jp_netas,_jp_etarange);
  pneftp_vseta = new TProfile("pneftp_vseta","",_jp_netas,_jp_etarange);
  pnhftp_vseta = new TProfile("pnhftp_vseta","",_jp_netas,_jp_etarange);
  pceftp_vseta = new TProfile("pceftp_vseta","",_jp_netas,_jp_etarange);
  pmuftp_vseta = new TProfile("pmuftp_vseta","",_jp_netas,_jp_etarange);
  pbetatp_vseta = new TProfile("pbetatp_vseta","",_jp_netas,_jp_etarange);
  pbetastartp_vseta = new TProfile("pbetastartp_vseta","",_jp_netas,_jp_etarange);
  pbetaprimetp_vseta = new TProfile("pbetaprimetp_vseta","",_jp_netas,_jp_etarange);

  curdir->cd();
}

void histosEta::Write() {
  dir->cd();
  dir->Write();
}

histosEta::~histosEta() {
  Write();
};
