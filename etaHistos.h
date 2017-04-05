// Purpose:  Histograms that include the whole eta profiles
// Author:   hannu.siikonen@cern.ch
// Created:  April 3, 2017

#ifndef __etaHistos_h__
#define __etaHistos_h__

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectory.h"

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "settings.h"

using std::string;

class etaHistos {
 public:
  // phase space
  string trigname;
  
  // Control plots of resolutions in the pt-eta plane
  TH3D *hdjasymm_a005;
  TH3D *hdjmpf_a005;
  TH3D *hdjasymmtp_a005;
  TH3D *hdjmpftp_a005;
  TH3D *hdjasymm_a01;
  TH3D *hdjmpf_a01;
  TH3D *hdjasymmtp_a01;
  TH3D *hdjmpftp_a01;
  TH3D *hdjasymm_a015;
  TH3D *hdjmpf_a015;
  TH3D *hdjasymmtp_a015;
  TH3D *hdjmpftp_a015;
  TH3D *hdjasymm_a02;
  TH3D *hdjmpf_a02;
  TH3D *hdjasymmtp_a02;
  TH3D *hdjmpftp_a02;
  TH3D *hdjasymm_a025;
  TH3D *hdjmpf_a025;
  TH3D *hdjasymmtp_a025;
  TH3D *hdjmpftp_a025;
  TH3D *hdjasymm_a03;
  TH3D *hdjmpf_a03;
  TH3D *hdjasymmtp_a03;
  TH3D *hdjmpftp_a03;

  etaHistos() {}
  etaHistos(TDirectory *dir, string trigname);
  ~etaHistos();

// private:
  TDirectory *dir;
};

#endif // __etaHistos_h__
