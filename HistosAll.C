// Purpose:  Create histos for all trigs
// Author:   hannu.siikonen@cern.ch
// Created:  November 19, 2018

#include "HistosAll.h"

HistosAll::HistosAll(TDirectory *dir) {

  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // Once and for all (even if few too many with Sumw2)
  TH1::SetDefaultSumw2(kTRUE);

  pmetave        = new TProfile("pmetave",       ";p_T^{tag};MET^{proj}_{ave};",jp::nwwpts,jp::wwptrange);
  pmetave_nol2l3 = new TProfile("pmetave_nol2l3",";p_T^{tag};MET^{proj}_{ave};",jp::nwwpts,jp::wwptrange);
  ppttagave        = new TProfile("ppttagave",       ";p_T^{tag};p^{proj}_{T,ave};",jp::nwwpts,jp::wwptrange);
  ppttagave_nol2l3 = new TProfile("ppttagave_nol2l3",";p_T^{tag};p^{proj}_{T,ave};",jp::nwwpts,jp::wwptrange);
  for (unsigned idx = 0; idx < jp::nwwpts; ++idx) {
    int num = jp::wwptrange[idx];
    string number = std::to_string(num);
    p2ptjet1.push_back           ( new TProfile2D((string("p2ptj1_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet2.push_back           ( new TProfile2D((string("p2ptj2_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet3.push_back           ( new TProfile2D((string("p2ptj3_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet4plus.push_back       ( new TProfile2D((string("p2ptj4p_")       +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet1_nol2l3.push_back    ( new TProfile2D((string("p2ptj1_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet2_nol2l3.push_back    ( new TProfile2D((string("p2ptj2_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet3_nol2l3.push_back    ( new TProfile2D((string("p2ptj3_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    p2ptjet4plus_nol2l3.push_back( new TProfile2D((string("p2ptj4p_nol2l3_")+number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet1.push_back           ( new TH2D((string("h2nj1_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet2.push_back           ( new TH2D((string("h2nj2_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet3.push_back           ( new TH2D((string("h2nj3_")        +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet4plus.push_back       ( new TH2D((string("h2nj4p_")       +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet1_nol2l3.push_back    ( new TH2D((string("h2nj1_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet2_nol2l3.push_back    ( new TH2D((string("h2nj2_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet3_nol2l3.push_back    ( new TH2D((string("h2nj3_nol2l3_") +number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
    h2njet4plus_nol2l3.push_back( new TH2D((string("h2nj4p_nol2l3_")+number+string("GeV")).c_str(),";p_T;#eta;",jp::npts_semi,jp::ptrange_semi,jp::netas,jp::etarange) );
  }
  h2ptjet1            = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet2            = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet3            = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet4plus        = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet1_nol2l3     = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet2_nol2l3     = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet3_nol2l3     = vector<TH2D*>(jp::nwwpts,0);
  h2ptjet4plus_nol2l3 = vector<TH2D*>(jp::nwwpts,0);

  curdir->cd();
}

void HistosAll::Write() {
  dir->cd();
  dir->Write();
}

HistosAll::~HistosAll() {
  Write();
};
