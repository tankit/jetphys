// JEC intervals of validity, used for DT
#ifndef JEC_IOV
#define JEC_IOV
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <map>
#include <string>
#include <iostream>
#include <utility>

#include "settings.h"

using namespace std;

namespace jec {
  
  struct IOVdata {
    vector<string> names;
    unsigned int low;
    unsigned int up;
    FactorizedJetCorrector* corr;
    FactorizedJetCorrector* l1rc;
    JetCorrectionUncertainty* unc;
  };

  // header part
  class IOV {

  public:
    IOV();
    ~IOV();
    void add(string id, string jecgt, string jecvers, unsigned int runmin, unsigned int runmax);
    bool setCorr(unsigned int,FactorizedJetCorrector**,FactorizedJetCorrector**,JetCorrectionUncertainty**);

  private:
    vector<IOVdata> _jecs;
    int _current;
  }; // class IOV

  IOV::IOV() : _current(-1) {}
  IOV::~IOV() {
    for (auto &jec: _jecs) {
      delete jec.corr;
      delete jec.l1rc;
      delete jec.unc;
    }
  }

  // body part (separate later to another file)
  void IOV::add(string id, string jecgt, string jecvers, unsigned int runmin, unsigned int runmax) {
    IOVdata dat;
    dat.low = runmin;
    dat.up = runmax;

    // sanity checks to avoid IOV overlaps
    assert(runmax>=runmin);
    for (auto it = _jecs.begin(); it != _jecs.end(); ++it) {
      assert(runmax<it->low or runmin>it->up);
    }

    const char *tmps;
    const char *path = "CondFormats/JetMETObjects/data/";
    string jectmp = jecgt + id + jecvers + (_jp_ismc ? "_MC_" : "_DATA_"); 
    const char *jecname = jecgt.c_str();

    vector<JetCorrectorParameters> vpar;

    // L1FastJet for AK*PF, L1Offset for others
    tmps = Form("%s%sL1FastJet_%s.txt",path,jecname,_jp_algo);
    if (access(tmps, F_OK)==-1) {
      cout << "The IOV files of " << tmps << " etc. were not found" << endl;
      return;
    }
    dat.names.push_back(string(tmps));
    vpar.push_back(JetCorrectorParameters(tmps));

    tmps = Form("%s%sL2Relative_%s.txt",path,jecname,_jp_algo);
    dat.names.push_back(string(tmps));
    vpar.push_back(JetCorrectorParameters(tmps));

    tmps = Form("%s%sL3Absolute_%s.txt",path,jecname,_jp_algo);
    dat.names.push_back(string(tmps));
    vpar.push_back(JetCorrectorParameters(tmps));

    if (_jp_isdt) {
      if (!_jp_skipl2l3res) {
        tmps = Form("%s%sL2L3Residual_%s.txt",path,jecname,_jp_algo);
        dat.names.push_back(string(tmps));
        vpar.push_back(JetCorrectorParameters(tmps));
      }
      tmps = Form("%s%sUncertainty_%s.txt",path,jecname,_jp_algo);
      dat.unc = new JetCorrectionUncertainty(tmps);
    }
    dat.corr = new FactorizedJetCorrector(vpar);

    // For type-I and type-II MET
    vector<JetCorrectorParameters> vrc;
    tmps = Form("%s%sL1RC_%s.txt",path,jecname,_jp_algo);
    vrc.push_back(JetCorrectorParameters(tmps));
    dat.l1rc = new FactorizedJetCorrector(vrc);

    _jecs.push_back(dat);
  } // add

  bool IOV::setCorr(unsigned int run,FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc) {
    assert(_jecs.size()!=0);
    auto cit = &_jecs[_current];
    if (cit->low <= run and cit->up >= run) // If no change is needed
      return true;

    // If the current IOV was not fine, search through all the IOVs
    for (int i = 0, N = int(_jecs.size()); i<N; ++i) {
      auto it = &_jecs[i];
      if (it->low <= run and it->up >= run) {
        if (i!=_current) {
          _current = i;
          *corr = it->corr;
          *l1rc = it->l1rc;
          *unc = it->unc;
          cout << endl << "IOV handling in use." << endl;
          for (auto &name: it->names) {
            cout << "Loading ... " << name << endl;
          }
        }
        return true;
      }
    }
    cout << "IOV for run " << run << " not found!!" << endl << flush;
    return false;
  } // getCorr


} // namescape jec
#endif
