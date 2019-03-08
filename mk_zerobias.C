#include "zerobias.h"

// JEC libraries
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)

//gROOT->ProcessLine(".L zerobias.C+");
R__LOAD_LIBRARY(zerobias.C+g);
R__LOAD_LIBRARY(zerobiasDraw.C+g);
//R__LOAD_LIBRARY(zerobiasDrawMu.C+g);


void mk_zerobias() {

  TChain *c = new TChain("T");
  
  // Desktop files
  //c->AddFile("/Volumes/LaCie_10TB/data/offset/Legacy_jets/Legacy16_RunG_unskimmed_jettuples.root");
  //c->AddFile("/Volumes/LaCie_10TB/data/offset/Legacy_jets/Legacy16_RunH_unskimmed_jettuples.root");
  //c->AddFile("/Volumes/LaCie_10TB/data/offset/Legacy_jets/SingleNeutrino_MC.root");
  
  //zerobias zb(c);
  //zb.Loop();

  //zerobiasDraw();
  //zerobiasDrawComp();
  zerobiasDrawMu();
}
