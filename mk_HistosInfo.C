// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: June 15, 2017

#include "settings.h"

#include "CondFormats/JetMETObjects/src/Utilities.cc"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// For JEC uncertainty
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "HistosRun.h"
#include "HistosBasic.h"
#include "HistosEta.h"
#include "HistosMC.h"

#include "histosInfo.h"

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
#ifdef USEASSERT
R__LOAD_LIBRARY(tools.C+g)
R__LOAD_LIBRARY(histosInfo.C+g)
#else
R__LOAD_LIBRARY(tools.C+)
R__LOAD_LIBRARY(histosInfo.C)
#endif // USEASSERT or not?

void mk_histosInfo() {

  string algo = "ak4";
  if (strings_equal(jp::algo,"AK8")) algo = "ak8";
  if (strings_equal(jp::algo,"AK4")) algo = "ak4";

  // This is already included in the .C files above
  // Including it again breaks CLING in ROOT 6.04.10

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",algo.c_str()));
  const char* p = "/work/jet_tuples/";
  //const char* p = "/Volumes/LaCie/data_2017/SMPJ/data/";
  //const char* p = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/Multijet/13TeV/";

  if (jp::isdt) {
    cout << "Running over DT" << endl;
    cout << "Load trees..." << endl;

    vector<string> files;
    if (strings_equal(jp::run,"RunB17")) {
      const char* ps = "Data/2017/";
      files.push_back(Form("%s%sRunB-2017-part1.root",p,ps));
      files.push_back(Form("%s%sRunB-2017-part2.root",p,ps));
      files.push_back(Form("%s%sRunB-2017-part3.root",p,ps));
    } else if (strings_equal(jp::run,"RunB16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunB16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunC16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunC16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunD16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunD16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunE16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunE16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunFearly16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunFearly16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunFlate16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunFlate16_Sep23_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunG16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunG16_Sep23_0.root",p,ps));
      files.push_back(Form("%s%sRunG16_Sep23_1.root",p,ps));
      files.push_back(Form("%s%sRunG16_Sep23_2.root",p,ps));
      files.push_back(Form("%s%sRunG16_Sep23_3.root",p,ps));
      files.push_back(Form("%s%sRunG16_Sep23_4.root",p,ps));
    } else if (strings_equal(jp::run,"RunH16_Han")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunH16_07Aug17_0.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_1.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_2.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_3.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_4.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_5.root",p,ps));
    } else if (strings_equal(jp::run,"RunFearly16_Gold")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunFearly16_07Aug17_Gold_0.root",p,ps));
      files.push_back(Form("%s%sRunFearly16_07Aug17_Gold_1.root",p,ps));
      files.push_back(Form("%s%sRunFearly16_07Aug17_Gold_2.root",p,ps));
      files.push_back(Form("%s%sRunFearly16_07Aug17_Gold_3.root",p,ps));
    } else if (strings_equal(jp::run,"RunFlate16_HanNew")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunFlate16_07Aug17_New_0.root",p,ps));
      files.push_back(Form("%s%sRunFlate16_07Aug17_New_1.root",p,ps));
    } else if (strings_equal(jp::run,"RunG16_HanNew")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunH16_07Aug17_New_0.root",p,ps));
    } else if (strings_equal(jp::run,"RunH16_HanNew")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sRunH16_07Aug17_New_0.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_1.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_2.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_3.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_4.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_5.root",p,ps));
      files.push_back(Form("%s%sRunH16_07Aug17_New_6.root",p,ps));
    } else if (strings_equal(jp::run,"RunB16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart5.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-LegacyReReco-80Xpart6.root",p,ps));
    } else if (strings_equal(jp::run,"RunC16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunC-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunC-LegacyReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(jp::run,"RunD16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunD-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunD-LegacyReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunD-LegacyReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(jp::run,"RunE16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunE-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunE-LegacyReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunE-LegacyReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(jp::run,"RunFearly16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunFearly-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunFearly-LegacyReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(jp::run,"RunFlate16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunFlate-LegacyReReco-80Xpart1.root",p,ps));
    } else if (strings_equal(jp::run,"RunG16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunG-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-LegacyReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-LegacyReReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-LegacyReReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-LegacyReReco-80Xpart5.root",p,ps));
    } else if (strings_equal(jp::run,"RunH16_Eng")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunH-LegacyReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-LegacyReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-LegacyReReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-LegacyReReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-LegacyReReco-80Xpart5.root",p,ps));
    } else if (strings_equal(jp::run,"RunB16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart5.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart6.root",p,ps));
    } else if (strings_equal(jp::run,"RunC16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunC-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunC-ReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(jp::run,"RunD16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunD-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunD-ReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunD-ReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(jp::run,"RunE16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunE-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunE-ReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunE-ReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(jp::run,"RunFearly16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunFearly-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunFearly-ReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(jp::run,"RunFlate16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunFlate-ReReco-80X.root",p,ps));
    } else if (strings_equal(jp::run,"RunG16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart5.root",p,ps));
    } else if (strings_equal(jp::run,"RunH16")) {
      const char* ps = "Data/2016/";
      files.push_back(Form("%s%sNtuples-Data-2016RunH-PromptReco-80Xpart1.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-PromptReco-80Xpart2.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-PromptReco-80Xpart3.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-PromptReco-80Xpart4.root",p,ps));
      files.push_back(Form("%s%sNtuples-Data-2016RunH-PromptReco-80Xpart5.root",p,ps));
    } else {
      cout << "Enter a proper value for jp::run!" << endl;
    }
    cout << "Loading the files:" << endl;
    for (auto &fname : files) {
      cout << fname << endl;
      c->AddFile(fname.c_str());
    }
  } else if (jp::ispy) {
    if (jp::pthatbins) {
      cout << "Running over pthat bins in pythia8" << endl;
      cout << "Load trees..." << endl;

      if (jp::pthatfiles.size()!=jp::npthatbins or jp::pthatsigmas.size()!=jp::npthatbins or
        jp::pthatnevts.size()!=jp::npthatbins or jp::pthatranges.size()!=jp::npthatbins+1) {
        cout << "The pthat file information dimensions don't match! Aborting..." << endl;
        cout << "Bins: " << jp::npthatbins << endl;
        cout << "Files: " << jp::pthatfiles.size() << endl;
        cout << "Sigmas: " << jp::pthatsigmas.size() << endl;
        cout << "Evts: " << jp::pthatnevts.size() << endl;
        cout << "Ranges-1: " << jp::pthatranges.size()-1 << endl;
        return;
      }

      const char* ps = "MC/P825ns80X_Moriond17/";
      for (auto & mcfile : jp::pthatfiles) {
        c->AddFile(Form("%s%s%s",p,ps,mcfile.c_str()));
      }
    } else {
      cout << "Running over pythia flat sample" << endl;

      cout << "Load trees..." << endl;
      c->AddFile(Form("%sP825ns80X/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8.root",p));
    }
  } else if (jp::ishw) {
    if (jp::pthatbins) {
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
    cout << "Enter a proper jp::type" << endl;
  }

  int centries = c->GetEntries();
  cout << "Got " << centries << " entries" << endl;

  if (centries > 0) {
    histosInfo filler(c);
	filler.Loop();
  } else {
    cout << "Please specify files to be looped over!" << endl << "Exiting ..." << endl;
  }
  //gROOT->ProcessLine(".q");
}
