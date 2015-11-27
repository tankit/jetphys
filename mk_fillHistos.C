// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: June 1, 2015

//void mk_fillHistos() {
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

  gROOT->ProcessLine(".L fillHistos.C+");

  // This is already included in the .C files above
  // Including it again breaks CLING in ROOT 6.04.10
  // #include "settings.h"

  std::string algo = "ak4";
  if (_jp_algo=="AK8") algo = "ak8";
  if (_jp_algo=="AK4") algo = "ak4";
  const char *a = algo.c_str();

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",a));
  TChain *d = (_jp_ak4ak8 ? new TChain("ak4/ProcessedTree") : 0);

  if (_jp_type=="DATA") {
    cout << "Load trees..." << endl;

    // Sep 19 files
    c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
    //c->AddFile("moredata/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root");

    // Wrong files...
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-Cert246908-251883_13TeV_PromptReco_Collisions15-LowPileUp.root");
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFIle-254833-Run2015C.root");
    cout << "Got " << c->GetEntries() << " entries" << endl;

    if (_jp_ak4ak8) {
      cout << "Load friend trees..." << endl;

      c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
      c->AddFile("moredata/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root");

      cout << "Got " << d->GetEntries() << " entries" << endl;
    }
  }

  if (_jp_type=="MC") {
    if (_jp_pthatbins) {
      cout << "Running over pthat bins" << endl;
      cout << "Load trees..." << endl;

      if (_jp_ak4ak8) {
	cout << "Load friend trees..." << endl;
      }
    }
    else {
      cout << "Running over flat sample" << endl;
      cout << "Load trees..." << endl;

      c->AddFile("moremc/Pythia8-CUETP8M1-Ntuples-PFJets.root");

      cout << "Got " << c->GetEntries() << " entries" << endl;

      if (_jp_ak4ak8) {
	cout << "Load friend trees..." << endl;
	cout << "Got " << d->GetEntries() << " entries" << endl;
      }
    }
  }

  if (_jp_type=="HW") {
    cout << "Running over Herwig flat sample" << endl;
    cout << "Load trees..." << endl;
    cout << "Got " << c->GetEntries() << " entries" << endl;

    if (_jp_ak4ak8) {
      cout << "Load friend trees..." << endl;
      cout << "Got " << d->GetEntries() << " entries" << endl;
    }
  }

  // Awkward patch for ROOT6:
  // Call 'Loop()' and 'delete this;' inside constructor;
  fillHistos(c, _jp_type, d);

}
