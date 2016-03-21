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

  gROOT->ProcessLine(".L fillHistos.C+g"); // +g for assert to work

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
  const char *p = "/Volumes/Macintosh HD 2Tb/data/QCDJETS/";

  if (_jp_type=="DATA") {
    cout << "Load trees..." << endl;

    // 2015D
    //c->AddFile("moredata/ProcessedTree_JetHT_2015_PromptReco_goldJSON_upto-260627_SMPJ_ak4_withfilters.root");
    // 76X with botched prescales (and other problems?)
    //c->AddFile("%sNtuples-Data-RunC-ReReco16Dec2015.root"); // different trigs
//     c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-part1.root",p));
//     c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-part2.root",p));
//     c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-part3.root",p));
//     c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-part4.root",p));
//     c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-part5.root",p));
    // New 76X with fixed prescales
    c->AddFile(Form("%sNtuples-Data-RunC-ReReco16Dec2015-76X.root",p));
    c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-76X-part1.root",p));
    c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-76X-part2.root",p));
    c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-76X-part3.root",p));
    c->AddFile(Form("%sNtuples-Data-RunD-ReReco16Dec2015-76X-part4.root",p));

    // Sep 19 files
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
    //c->AddFile("moredata/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root");

    // Wrong files...
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-Cert246908-251883_13TeV_PromptReco_Collisions15-LowPileUp.root");
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
    //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFIle-254833-Run2015C.root");
    cout << "Got " << c->GetEntries() << " entries" << endl;

    if (_jp_ak4ak8) {
      cout << "Load friend trees..." << endl;

      //c->AddFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
      //c->AddFile("moredata/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root");

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

      //c->AddFile("moremc/ProcessedTree_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v3.root"); //  25 ns
      //c->AddFile("moremc/Pythia8-CUETP8M1-Ntuples-PFJets.root"); // 50 ns
      c->AddFile(Form("%sNtuples-MC-Pythia8-Flat15to7000-25ns-CUETM1-13TeV.root",p));//76

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
