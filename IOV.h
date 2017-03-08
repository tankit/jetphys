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

using namespace std;

namespace jec {
  
  struct IOVdata {
    unsigned int low;
    unsigned int up;
    FactorizedJetCorrector* corr;
    FactorizedJetCorrector* l1rc;
    JetCorrectionUncertainty* unc;
  };

  // header part
  class IOV {

  public:
    IOV(string algo = "AK4PF") : _algo(algo) {
    };
    ~IOV() {}; // leaks memory, fix later
    void add(string id, string jecgt, string jecvers, unsigned int runmin, unsigned int runmax);
    bool setCorr(unsigned int run,FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc);

  private:
    string _algo;
    vector<IOVdata> jecs;
  }; // class IOV

  // body part (separate later to another file)
  void IOV::add(string id, string jecgt, string jecvers, unsigned int runmin, unsigned int runmax) {
    IOVdata dat;
    dat.low = runmin;
    dat.up = runmax;

    // sanity checks to avoid IOV overlaps
    assert(runmax>=runmin);
    for (auto it = jecs.begin(); it != jecs.end(); ++it) {
      assert( runmax<it->low || runmin>it->up );
    }
   
    jecgt = jecgt + id + jecvers + "_DATA_"; 
    const char *s;
    const char *p = "CondFormats/JetMETObjects/data/";
    const char *t = jecgt.c_str();;
    const char *a = _algo.c_str();

    // L1FastJet for AK*PF, L1Offset for others
    s = Form("%s%sL1FastJet_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1 = new JetCorrectorParameters(s);
    s = Form("%s%sL2Relative_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2 = new JetCorrectorParameters(s);
    s = Form("%s%sL3Absolute_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l3 = new JetCorrectorParameters(s);
    s = Form("%s%sL2L3Residual_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> vpar;
    vpar.push_back(*par_l1);
    vpar.push_back(*par_l2);
    vpar.push_back(*par_l3);
    vpar.push_back(*par_l2l3res);
    dat.corr = new FactorizedJetCorrector(vpar);

    // For type-I and type-II MET
    s = Form("%s%sL1FastJet_%s.txt",p,t,a);
    JetCorrectorParameters *par_l1rc = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> vrc;
    vrc.push_back(*par_l1rc);
    dat.l1rc = new FactorizedJetCorrector(vrc);

    s = Form("%s%sUncertainty_%s.txt",p,t,a);
    dat.unc = new JetCorrectionUncertainty(s);

    jecs.push_back(dat);
  } // add
  
  bool IOV::setCorr(unsigned int run,FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc) {
    assert(jecs.size()!=0);
    for (auto it = jecs.begin(); it != jecs.end(); ++it) {
      if (it->low <= run && it->up >= run) {
        *corr = it->corr;
        *l1rc = it->l1rc;
        *unc = it->unc;
        return true;
      }
    }
    cout << "IOV for run " << run << " not found!!" << endl << flush;
    return false;
  } // getCorr


} // namescape jec
#endif
