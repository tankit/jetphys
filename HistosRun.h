// Purpose:  Keep track of trigger rates and prescales per run
//
// Author:   mikko.voutilainen@cern.ch
// Created:  June 4, 2010
// Updated:  June 4, 2010

#ifndef __HistosRun_h__
#define __HistosRun_h__

#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TMath.h"

#include <string>
#include <map>

#include "settings.h"

class HistosRun {

 public:

  double etamin;
  double etamax;

  std::vector<std::string> trg;
  std::map<std::string, double> pt;

  double lumsum;
  double lumsum2;
  std::map<int, std::map<int, float> > lums;
  std::map<int, float> runlums;
  std::map<int, float> runlums2;
  std::map<std::string, std::map<int, float> > runlums_trg;

  // trigger rate counters
  std::map<std::string, std::map<int, int> > t_trg;
  std::map<std::string, std::map<int, int> > tw_trg; // prescale weight
  std::map<std::string, std::map<int, int> > t_trgtp;
  std::map<std::string, std::map<int, int> > c_trg; // calo

  // Average number of vertices
  std::map<std::string, std::map<int, int> > npv_trg;
  std::map<std::string, std::map<int, int> > npvgood_trg;

  // Average component fractions
  std::map<std::string, std::map<int, double> > c_chf;
  std::map<std::string, std::map<int, double> > c_nef;
  std::map<std::string, std::map<int, double> > c_nhf;
  std::map<std::string, std::map<int, double> > c_betastar;
  std::map<std::string, std::map<int, double> > c_betaprime;
  // Same with tag-and-probe method
  std::map<std::string, std::map<int, double> > c_chftp;
  std::map<std::string, std::map<int, double> > c_neftp;
  std::map<std::string, std::map<int, double> > c_nhftp;
  std::map<std::string, std::map<int, double> > c_betastartp;
  std::map<std::string, std::map<int, double> > c_betaprimetp;

  // prescale counters
  std::map<std::string, std::map<int, int> > p_trg;
  std::map<std::string, std::map<int, int> > p_trgpair;

  // Other non-lumi checks

  HistosRun(TDirectory *dir, double etamin = 0., double etamax = 3.);
  ~HistosRun();

 private:
  TDirectory *dir;
};

#endif // __HistosRun_h__
