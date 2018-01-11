// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: June 15, 2017

#include "settings.h"

#include "CondFormats/JetMETObjects/src/Utilities.cc"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
// For JEC uncertainty
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)

#include "tools.h"
R__LOAD_LIBRARY(tools.C+)
#include "runHistos.h"
R__LOAD_LIBRARY(runHistos.C+)
#include "basicHistos.h"
R__LOAD_LIBRARY(basicHistos.C+)
#include "etaHistos.h"
R__LOAD_LIBRARY(etaHistos.C+)
#include "mcHistos.h"
R__LOAD_LIBRARY(mcHistos.C+)

#include "fillHistos.h"
R__LOAD_LIBRARY(fillHistos.C+g) // g for assert to work

void mk_fillHistos() {

  string algo = "ak4";
  if (strings_equal(_jp_algo,"AK8")) algo = "ak8";
  if (strings_equal(_jp_algo,"AK4")) algo = "ak4";

  // This is already included in the .C files above
  // Including it again breaks CLING in ROOT 6.04.10

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",algo.c_str()));
  const char* p = "/work/jet_tuples/";
  //const char* p = "/Volumes/LaCie/data_2017/SMPJ/data/";
  //const char* p = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/Multijet/13TeV/";

  if (_jp_isdt) {
    cout << "Running over DT" << endl;
    cout << "Load trees..." << endl;

    if (strings_equal(_jp_run,"RunB17")) {
      const char* ps = "Data/2017/";
      c->AddFile(Form("%s%sRunB-2017-part1.root",p,ps));
      c->AddFile(Form("%s%sRunB-2017-part2.root",p,ps));
      c->AddFile(Form("%s%sRunB-2017-part3.root",p,ps));
    } else if (strings_equal(_jp_run,"RunB16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart2.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart3.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart4.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart5.root",p,ps));
      //c->AddFile(Form("%s%sNtuples-Data-2016RunB-ReReco-80Xpart6.root",p,ps));
    } else if (strings_equal(_jp_run,"RunC16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sNtuples-Data-2016RunC-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunC-ReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(_jp_run,"RunD16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart2.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(_jp_run,"RunE16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart2.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart3.root",p,ps));
    } else if (strings_equal(_jp_run,"RunFearly16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart2.root",p,ps));
    } else if (strings_equal(_jp_run,"RunFlate16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunFlate-ReReco-80X.root",p,ps));
    } else if (strings_equal(_jp_run,"RunG16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart2.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart3.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart4.root",p,ps));
      c->AddFile(Form("%s%sNtuples-Data-2016RunG-ReReco-80Xpart5.root",p,ps));
    } else if (strings_equal(_jp_run,"RunH16")) {
      const char* ps = "Data/2016/";
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart1.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart2.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart3.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart4.root",p,ps));
      c->AddFile(Form("%s%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart5.root",p,ps));
    } else {
      cout << "Enter a proper value for _jp_run!" << endl;
    }
  } else if (_jp_ispy) {
    if (_jp_pthatbins) {
      cout << "Running over pthat bins in pythia8" << endl;
      cout << "Load trees..." << endl;

      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_30to50_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_50to80_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_80to120_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_120to170_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_170to300_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_300to470_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_470to600_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_600to800_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_800to1000_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_1000to1400_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_1400to1800_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_1800to2400_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_2400to3200_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sP825ns80X_Moriond17/QCD_Pt_3200toInf_TuneCUETP8M_13TeV_pythia8.root",p));

    } else {
      cout << "Running over pythia flat sample" << endl;

      cout << "Load trees..." << endl;
      c->AddFile(Form("%sP825ns80X/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8.root",p));
    }
  } else if (_jp_ishw) {
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

  if (centries > 0) {
    fillHistos filler(c);
	filler.Loop();
  } else {
    cout << "Please specify files to be looped over!" << endl << "Exiting ..." << endl;
  }
  //gROOT->ProcessLine(".q");
}
