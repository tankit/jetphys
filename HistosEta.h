// Purpose:  Histograms that include the whole eta profiles
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#ifndef __HistosEta_h__
#define __HistosEta_h__

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectory.h"
#include "TMath.h"

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "settings.h"

using std::string;
using std::vector;

class HistosEta {
  public:
  // phase space
  string trigname;

  double pttrg;
  double ptmin;
  double ptmax;

  bool ismcdir;

  // control plots vs eta
  TProfile *pncandtp_vseta;
  TProfile *pnchtp_vseta;
  TProfile *pnnetp_vseta;
  TProfile *pnnhtp_vseta;
  TProfile *pncetp_vseta;
  TProfile *pnmutp_vseta;
  TProfile *pnhhtp_vseta;
  TProfile *pnhetp_vseta;
  TProfile *pchftp_vseta;
  TProfile *pneftp_vseta;
  TProfile *pnhftp_vseta;
  TProfile *pceftp_vseta;
  TProfile *pmuftp_vseta;
  TProfile *phhftp_vseta;
  TProfile *pheftp_vseta;
  TProfile *pbetatp_vseta;
  TProfile *pbetastartp_vseta;
  TProfile *ppuftp_vseta;

  // Control plots of resolutions in the pt-eta plane
  // tp: tag-probe (tag), pt: probe-tag (probe)
  vector<TH3D *> hdjasymm;
  vector<TH3D *> hdjasymmtp;
  //vector<TH3D *> hdjasymmpt;
  vector<TH3D *> hdjmpf;
  vector<TH3D *> hdjmpftp;
  //vector<TH3D *> hdjmpfpt;

  TProfile3D *p3rvsnpv;
  TProfile3D *p3rvsnpvW;

  const vector<float> alpharange = {0.05,0.10,0.15,0.20,0.25,0.30};

  HistosEta() {}
  HistosEta(TDirectory *dir, string trigname, double pttrg, double ptmin, double ptmax, bool ismcdir = false);
  ~HistosEta();

  void Write();

 private:
  TDirectory *dir;
};

#endif // __HistosEta_h__
