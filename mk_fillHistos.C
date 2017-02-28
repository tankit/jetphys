// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: June 1, 2015

{
  //gROOT->ProcessLine(".exception");

  // Jet energy corrections and uncertainties
  // Localized version available in github.com/miquork/jecsys/
  // Copied from there
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // compile code
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L runHistos.C+");
  gROOT->ProcessLine(".L basicHistos.C+");

  gROOT->ProcessLine(".L fillHistos.C+g"); // +g for assert to work

  // This is already included in the .C files above
  // Including it again breaks CLING in ROOT 6.04.10
  //#include "settings.h"

  string algo = "ak4";
  if (_jp_algo=="AK8") algo = "ak8";
  if (_jp_algo=="AK4") algo = "ak4";

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",algo.c_str()));
  TChain *d = (_jp_ak4ak8 ? new TChain("ak4/ProcessedTree") : 0);
  const char* p = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/Multijet/13TeV/";

  if (_jp_type=="DATA") {
    cout << "Running over DT" << endl;
    cout << "Load trees..." << endl;

    // Uncomment the one you want to use:

  //// Run2 2016 ReReco
    // 2016B
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart3.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart4.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart5.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart6.root",p));
    // 2016C
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunC-ReReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunC-ReReco-80Xpart2.root",p));
    // 2016D
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart3.root",p));
    // 2016E
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart3.root",p));
    // 2016Fearly
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart2.root",p));
    // 2016Flate
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFlate-ReReco-80X.root",p));
    // 2016G
    c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart1.root",p));
    c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart2.root",p));
    c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart3.root",p));
    c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart4.root",p));
    c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart5.root",p));

  //// Run2 2016 PromptReco
    // 2016B
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunB-part1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunB-part2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunB-part3.root",p));
    // 2016C
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunC-v2-part1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunC-v2-part2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunC-v2-part3.root",p));
    // 2016D
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunD-part1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunD-part2.root",p));
    // 2016E
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunE-part1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunE-part2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016-RunE-part3.root",p));
    // 2016F
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunF-PromptReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunF-PromptReco-80Xpart2.root",p));
    // 2016G
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart3.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart4.root",p));
    // 2016H
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart1.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart2.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart3.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart4.root",p));
    //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart5.root",p));
  } else if (_jp_type=="MC") {
    if (_jp_pthatbins) {
      cout << "Running over pthat bins" << endl;

      cout << "Load trees..." << endl;
      cout << "Files for a pthat bin run not specified" << endl;
    } else {
      cout << "Running over pythia flat sample" << endl;

      cout << "Load trees..." << endl;
      c->AddFile(Form("%sMC/Ntuples-MC-Pythia8-Flat15to7000-25ns-CUETM1-13TeV.root",p));
    }
  } else if (_jp_type=="HW") {
    if (_jp_pthatbins) {
      cout << "Running over pthat bins" << endl;

      cout << "Load trees..." << endl;
      cout << "Files for a pthat bin run not specified" << endl;
    } else {
      cout << "Running over Herwig flat sample" << endl;
      cout << "Load trees..." << endl;
      //c->AddFile(Form("%sMC/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_80X.root",p));
      c->AddFile(Form("%sMC/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_80X_v2.root",p));
    }
  } else {
    cout << "Enter a proper _jp_type" << endl;
  }

  int centries = c->GetEntries();
  cout << "Got " << centries << " entries" << endl;
  if (_jp_ak4ak8) {
    cout << "No ak8 trees currently specified" << endl;
    //cout << "Got " << d->GetEntries() << " entries" << endl;
  }

  if (centries > 0) {
    // ROOT 6 does not load the libraries before it thinks it has to.
    // Thus everything needs to be done in the initializer
    fillHistos(c, _jp_type, d);
  } else {
    cout << "Please specify files to be looped over!" << endl << "Exiting ..." << endl;
  }
  //gROOT->ProcessLine(".q");
}
