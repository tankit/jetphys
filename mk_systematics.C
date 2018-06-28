// Purpose: Estimation of inclusive jet cross section systematics
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: March 22, 2010
{

  gROOT->ProcessLine(".exception");

  // compile code
  //gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L tools.C+");
  //gROOT->ProcessLine(".L unfold.C+"); // for jer_systematics => include .C

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L systematics.C+");

  //gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  //cafe::Config *cfg = new cafe::Config("pfjet");
  //string type = cfg->get("type","DATA");
  //delete cfg;

  #include "settings.h"

  systematics(jp::type);
  //systematics(type, true); // 38X JEC; always run after 36X
  //jec_systematics();
  //jer_systematics();
  //jid_systematics();
  //lum_systematics();

  sources(jp::type);
}
