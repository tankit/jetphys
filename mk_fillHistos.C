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
  if (_jp_algo=="AK8") algo = "ak8";
  if (_jp_algo=="AK4") algo = "ak4";

  // This is already included in the .C files above
  // Including it again breaks CLING in ROOT 6.04.10

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",algo.c_str()));
  const char* p = "/work/jet_tuples/";
  //const char* p = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/Multijet/13TeV/";

  if (_jp_type=="DATA") {
    cout << "Running over DT" << endl;
    cout << "Load trees..." << endl;

    // Run2 2016 Newest (updated March 7th, 2017)
    //switch (_jp_run) {
    if (_jp_run=="RunB") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart2.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart3.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart4.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart5.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunB-ReReco-80Xpart6.root",p));
    } else if (_jp_run=="RunC") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunC-ReReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunC-ReReco-80Xpart2.root",p));
    } else if (_jp_run=="RunD") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart2.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunD-ReReco-80Xpart3.root",p));
    } else if(_jp_run=="RunE") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart2.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunE-ReReco-80Xpart3.root",p));
    } else if(_jp_run=="RunFearly") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFearly-ReReco-80Xpart2.root",p));
    } else if(_jp_run=="RunFlate") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunFlate-ReReco-80X.root",p));
    } else if(_jp_run=="RunG") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart1.root",p));
      //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart2.root",p));
      //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart3.root",p));
      //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart4.root",p));
      //c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunG-ReReco-80Xpart5.root",p));
    } else if(_jp_run=="RunH") {
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart1.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart2.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart3.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart4.root",p));
      c->AddFile(Form("%sData/2016/Ntuples-Data-2016RunH-PromptReco-80Xpart5.root",p));
    } else {
      cout << "Enter a proper value for _jp_run!" << endl;
    }

  } else if (_jp_type=="MC") {
    if (_jp_pthatbins) {
      cout << "Running over pthat bins in pythia8" << endl;
      cout << "Load trees..." << endl;

      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_30to50_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_50to80_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_80to120_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_120to170_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_170to300_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_300to470_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_470to600_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_600to800_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_800to1000_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_1000to1400_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_1400to1800_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_1800to2400_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_2400to3200_TuneCUETP8M_13TeV_pythia8.root",p));
      c->AddFile(Form("%sMC/P825ns80X_Moriond17/QCD_Pt_3200toInf_TuneCUETP8M_13TeV_pythia8.root",p));

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

  if (centries > 0) {
    fillHistos filler(c);
	filler.Loop();
  } else {
    cout << "Please specify files to be looped over!" << endl << "Exiting ..." << endl;
  }
  //gROOT->ProcessLine(".q");
}
