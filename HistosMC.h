// Purpose:  Histograms that include the whole eta profiles
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#ifndef __HistosMC_h__
#define __HistosMC_h__

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

class HistosMC {
 public:
  // phase space
  string trigname;
  
  // Control plots of resolutions in the pt-eta plane
  vector<TH3D *> hdjasymm;
  vector<TH3D *> hdjasymmtp;
  vector<TH3D *> hdjresp_tag; 
  vector<TH3D *> hdjresptp_tag; 
  vector<TH3D *> hdjresp_probe;
  vector<TH3D *> hdjresptp_probe;

  const vector<float> alpharange = {0.05,0.10,0.15,0.20,0.25,0.30};

  HistosMC() {}
  HistosMC(TDirectory *dir, string trigname);
  ~HistosMC();

  void Write();

 private:
  TDirectory *dir;
};

#endif // __HistosMC_h__
