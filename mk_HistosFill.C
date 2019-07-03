// Purpose: Fill jet physics analysis histograms
// Author:    mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: March 20, 2010

#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// For JEC uncertainty
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "settings.h"
#include "tools.h"
#include "HistosRun.h"
#include "HistosBasic.h"
#include "HistosEta.h"
#include "HistosMC.h"
#include "HistosAll.h"
#include "HistosFill.h"

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
#if JETRESO == 1
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetResolutionObject.cc+)
R__LOAD_LIBRARY(JetMETCorrections/Modules/src/JetResolution.cc+)
#endif
#ifdef USEASSERT
R__LOAD_LIBRARY(tools.C+g)
R__LOAD_LIBRARY(HistosRun.C+g)
R__LOAD_LIBRARY(HistosBasic.C+g)
R__LOAD_LIBRARY(HistosEta.C+g)
R__LOAD_LIBRARY(HistosMC.C+g)
R__LOAD_LIBRARY(HistosAll.C+g)
R__LOAD_LIBRARY(HistosFill.C+g)
#else
R__LOAD_LIBRARY(tools.C+)
R__LOAD_LIBRARY(HistosRun.C+)
R__LOAD_LIBRARY(HistosBasic.C+)
R__LOAD_LIBRARY(HistosEta.C+)
R__LOAD_LIBRARY(HistosMC.C+)
R__LOAD_LIBRARY(HistosAll.C+)
R__LOAD_LIBRARY(HistosFill.C+)
#endif // USEASSERT or not?

void mk_HistosFill() {

  bool fail = false;
  // Check the lengths of trigger infos & such
  if (jp::isdt) {
    if (jp::triggers.size()!=jp::notrigs) {
      cout << "jp::trigger size mismatch" << endl;
      fail = true;
    } else if (jp::trigthr.size()!=jp::notrigs) {
      cout << "jp::trigthr size mismatch" << endl;
      fail = true;
    } else if (jp::trigranges.size()!=jp::notrigs) {
      cout << "jp::trigranges size mismatch" << endl;
      fail = true;
    } else if (jp::triglumi.size()!=jp::notrigs) {
      cout << "jp::triglumi size mismatch" << endl;
      fail = true;
    } else if (jp::triglumi.size()!=jp::notrigs) {
      cout << "jp::triglumi size mismatch" << endl;
      fail = true;
    } else if (jp::triglumiera.size()!=jp::eras.size()) {
      cout << "Era count mismatch in jp::triglumiera and jp::eras" << endl;
      fail = true;
    } else if (jp::dtfiles.size()<jp::eras.size()) {
      cout << "Data files not provided for all eras!" << endl;
      fail = true;
    } else if (jp::IOVnames.size()!=jp::IOVranges.size()) {
      cout << "IOV count mismatch in jp::IOVnames and jp::IOVranges" << endl;
    } else {
      for (auto &lumis : jp::triglumiera) {
        if (lumis.size()!=jp::notrigs) {
          cout << "jp::triglumiera size mismatch" << endl;
          fail = true;
          break;
        }
      }
    }
  }
  if (jp::ispy) {
    if (jp::pthatbins) {
      // Check that pthat file dimensions are ok
      if (jp::pthatfiles.size()<jp::npthatbins) {
        cout << "The pthat info dimensions don't match! jp::pthatfiles vs jp::npthatbins" << endl;
        fail = true;
      } else if (jp::pthatsigmas.size()!=jp::npthatbins) {
        cout << "The pthat info dimensions don't match! jp::pthatsigmas vs jp::npthatbins" << endl;
        fail = true;
      } else if (jp::pthatranges.size()!=jp::npthatbins+1) {
        cout << "The pthat info dimensions don't match! jp::pthatranges vs jp::npthatbins" << endl;
        fail = true;
      }
    } else if (jp::htbins) {
      // Check that ht file dimensions are ok
      if (jp::htfiles.size()<jp::nhtbins) {
        cout << "The ht info dimensions don't match! jp::htfiles vs jp::nhtbins" << endl;
        fail = true;
      } else if (jp::htsigmas.size()!=jp::nhtbins) {
        cout << "The ht info dimensions don't match! jp::htsigmas vs jp::nhtbins" << endl;
        fail = true;
      } else if (jp::htranges.size()!=jp::nhtbins+1) {
        cout << "The ht info dimensions don't match! jp::htranges vs jp::nhtbins" << endl;
        fail = true;
      }
    }
  }
  if (fail) return;

  string algo = "ak4";
  if (jp::strings_equal(jp::algo,"AK8")) algo = "ak8";

  // connect trees
  TChain *c = new TChain(Form("%s/ProcessedTree",algo.c_str()));
  const char* p = jp::filepath;

  vector<const char*> files;
  if (jp::isdt) {
    cout << "Running over DT" << endl;
    cout << "Load trees..." << endl;

    assert(jp::dtfiles.size()>=jp::eras.size());
    const char* ps = jp::dtpath;

    for (auto &fname : jp::dtfiles.at(jp::run))
      files.push_back(Form("%s%s%s",p,ps,fname));
  } else if (jp::ispy) {
    if (jp::pthatbins) {
      cout << "Running over pthat binned files in pythia8" << endl;
      cout << "Load trees..." << endl;

      const char *ps = jp::pthatpath;

      for (auto &fname : jp::pthatfiles) {
        const char *name1 = Form("%s%s%s.root",p,ps,fname);
        const char *name2 = Form("%s%s%s_ext.root",p,ps,fname);
        std::ifstream stream1(name1);
        std::ifstream stream2(name2);
        if (stream1.good()) files.push_back(name1);
        if (stream2.good()) files.push_back(name2);
      }
    } else if (jp::htbins) {
      cout << "Running over ht binned files in pythia8" << endl;
      cout << "Load trees..." << endl;

      const char *ps = jp::htpath;

      for (auto &fname : jp::htfiles) {
        const char *name1 = Form("%s%s%s.root",p,ps,fname);
        const char *name2 = Form("%s%s%s_ext.root",p,ps,fname);
        std::ifstream stream1(name1);
        std::ifstream stream2(name2);
        if (stream1.good()) files.push_back(name1);
        if (stream2.good()) files.push_back(name2);
      }
    } else {
      cout << "Running over pythia flat sample" << endl;

      for (auto &fname : jp::mcfiles.at(jp::mcfile))
        files.push_back(Form("%s%s",p,fname));
    }
  } else if (jp::ishw) {
    if (jp::pthatbins) cout << "No pthat binned files exist for Herwig!" << endl;
    else if (jp::htbins) cout << "No ht binned files exist for Herwig!" << endl;

    cout << "Running over Herwig flat sample" << endl;

    for (auto &fname : jp::mcfiles.at(jp::hwfile)) files.push_back(Form("%s%s",p,fname));
  } else if (jp::isnu) {
    if (jp::pthatbins) cout << "No pthat binned files exist for Neutrino Gun!" << endl;
    else if (jp::htbins) cout << "No ht binned files exist for Neutrino Gun!" << endl;

    cout << "Running over HSingle Neutrino sample" << endl;

    for (auto &fname : jp::mcfiles.at(jp::nufile)) files.push_back(Form("%s%s",p,fname));
  } else {
    cout << "Enter a proper type!!" << endl;
  }

  if (files.size()==0) {
    if (jp::isdt) {
      cout << "Enter a proper value for run!" << endl;
      cout << "Entered: " << jp::run << endl << "Options:";
      for (auto &fname : jp::dtfiles) cout << " " << fname.first;
      cout << endl;
    } else if (jp::ispy or jp::ishw or jp::isnu) {
      if (jp::pthatbins) {
        cout << "Problems with pthat file logic!" << endl;
      } else if (jp::htbins) {
        cout << "Problems with ht file logic!" << endl;
      } else {
        cout << "Enter a proper value for MC tag!" << endl;
        cout << "Entered: " << jp::mcfile << endl << "Options:";
        for (auto &fname : jp::mcfiles) cout << " " << fname.first;
        cout << endl;
      }
    }
  } else {
    cout << "Loading the files:" << endl;
    for (auto &fname : files) {
      cout << fname << endl;
      c->AddFile(fname);
    }
  }

  int centries = c->GetEntries();
  cout << "Got " << centries << " entries" << endl;

  if (centries > 0) {
    try {
      HistosFill filler(c);
      filler.Loop();
    } catch (const exception& e) {
      cout << e.what() << endl;
    }
  } else {
    cout << "Please specify files to be looped over!" << endl << "Exiting ..." << endl;
  }
  gROOT->ProcessLine(".q");
}
