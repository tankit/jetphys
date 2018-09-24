// Purpose: Estimation of inclusive jet cross section systematics
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: March 22, 2010

//  #include "settings.h"
 #include "tools.h"

#include "CondFormats/JetMETObjects/src/Utilities.cc"

// #include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// For JEC uncertainty
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/Utilities.cc+);
//  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
  // For JEC uncertainty
  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
  R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
  //  R__LOAD_LIBRARY(systematics.C)
  R__LOAD_LIBRARY(tools.C+g)



void mk_systematics() {

  // gROOT->ProcessLine(".exception"); // mikä

  // compile code
  //gROOT->ProcessLine(".L Config.cpp+");
  // gROOT->ProcessLine(".L tools.C+");
  //gROOT->ProcessLine(".L unfold.C+"); // for jer_systematics => include .C


  //gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  //cafe::Config *cfg = new cafe::Config("pfjet");
  //string type = cfg->get("type","DATA");
  //delete cfg;


  systematics(jp::type);
  //systematics(type, true); // 38X JEC; always run after 36X
  //jec_systematics();
  //jer_systematics();
  //jid_systematics();
  //lum_systematics();

  // sources(jp::type);
}
