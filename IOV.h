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
      assert( runmax<it->low || runmin>it->up );
    }
   
    jecgt = jecgt + id + jecvers + "_DATA_"; 
    const char *s;
    const char *p = "CondFormats/JetMETObjects/data/";
    const char *t = jecgt.c_str();;
    const char *a = _jp_algo.c_str();

    vector<JetCorrectorParameters> vpar;

    // L1FastJet for AK*PF, L1Offset for others
    s = Form("%s%sL1FastJet_%s.txt",p,t,a);
    dat.names.push_back(string(s));
    vpar.push_back(JetCorrectorParameters(s));

    s = Form("%s%sL2Relative_%s.txt",p,t,a);
    dat.names.push_back(string(s));
    vpar.push_back(JetCorrectorParameters(s));

    s = Form("%s%sL3Absolute_%s.txt",p,t,a);
    dat.names.push_back(string(s));
    vpar.push_back(JetCorrectorParameters(s));

    if (!_jp_skipl2l3res) {
      s = Form("%s%sL2L3Residual_%s.txt",p,t,a);
      dat.names.push_back(string(s));
      vpar.push_back(JetCorrectorParameters(s));
    }

    dat.corr = new FactorizedJetCorrector(vpar);

    // For type-I and type-II MET
    vector<JetCorrectorParameters> vrc;
    s = Form("%s%sL1RC_%s.txt",p,t,a);
    vrc.push_back(JetCorrectorParameters(s));
    dat.l1rc = new FactorizedJetCorrector(vrc);

    s = Form("%s%sUncertainty_%s.txt",p,t,a);
    dat.unc = new JetCorrectionUncertainty(s);

    _jecs.push_back(dat);
  } // add
  
  bool IOV::setCorr(unsigned int run,FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc) {
    assert(_jecs.size()!=0);
    for (int i = 0, N = int(_jecs.size()); i<N; ++i) {
      auto it = &_jecs[i];
      if (it->low <= run && it->up >= run) {
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
