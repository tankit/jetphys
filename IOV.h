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
#if JETRESO == 1
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#endif

struct JECIOV {
  vector<string> names;
  int low;
  int up;
  FactorizedJetCorrector* corr;
  FactorizedJetCorrector* l1rc;
  JetCorrectionUncertainty* unc;
};

// header part
class IOV {

public:
  IOV();
  ~IOV();
  void addJEC(string id, int runmin = -1, int runmax = -1);
  void addJRC(string id, int runmin = -1, int runmax = -1);
  bool setJEC(FactorizedJetCorrector**,FactorizedJetCorrector**,JetCorrectionUncertainty** unc=0,int run = -1);

private:
  void addJECLvl(JECIOV &jec, vector<JetCorrectorParameters> &vpar, string tag, bool unc);
  vector<JECIOV> _jecs;
  int _current;
  string _jecname;
  string _jecpath;
}; // class IOV

IOV::IOV() : _current(-1) {}
IOV::~IOV() {}

void IOV::addJECLvl(JECIOV &jec, vector<JetCorrectorParameters> &vpar, string tag, bool unc = false) {
  string path = _jecpath+_jecname+"_"+tag+"_"+jp::algo+".txt";
  if (access(path.c_str(), F_OK)==-1) {
    cout << "The IOV file" << path << " was not found" << endl;
  } else {
    jec.names.push_back(path);
    if (unc) jec.unc = new JetCorrectionUncertainty(path.c_str());
    else vpar.push_back(JetCorrectorParameters(path.c_str()));
  }
}

void IOV::addJEC(string id, int runmin, int runmax) {
  vector<JetCorrectorParameters> vpar, vrc;
  JECIOV jec;
  jec.low = runmin;
  jec.up = runmax;

  if (runmin != -1) {
    // sanity checks to avoid IOV overlaps
    assert(runmax>=runmin);
    for (auto it = _jecs.begin(); it != _jecs.end(); ++it) assert(runmax<it->low or runmin>it->up);
  }

  _jecname = jp::jecgt + id + (jp::ismc ? jp::jecversmc : jp::jecversdt) + (jp::ismc ? "_MC" : "_DATA");
  _jecpath = jp::useJECDB ? (jp::JECDBLoc+"/textFiles/"+_jecname+"/") : "CondFormats/JetMETObjects/data/";

  addJECLvl(jec,vpar,"L1FastJet"); // L1Offset for others than AK*PF
  addJECLvl(jec,vpar,"L2Relative");
  addJECLvl(jec,vpar,"L3Absolute");

  if (jp::isdt) {
    if (!jp::skipl2l3res) addJECLvl(jec,vpar,"L2L3Residual");
    addJECLvl(jec,vpar,"Uncertainty",true);
  }
  jec.corr = new FactorizedJetCorrector(vpar);

  // For type-I and type-II MET
  vrc.push_back(JetCorrectorParameters(string(_jecpath+_jecname+"_L1RC_"+jp::algo+".txt").c_str()));
  jec.l1rc = new FactorizedJetCorrector(vrc);

  _jecs.push_back(jec);
} // addJEC

void IOV::addJRC(string id, int runmin, int runmax) {
  vector<JetCorrectorParameters> vpar, vrc;
  JECIOV jec;
  jec.low = runmin;
  jec.up = runmax;

  if (runmin != -1) {
    // sanity checks to avoid IOV overlaps
    assert(runmax>=runmin);
    for (auto it = _jecs.begin(); it != _jecs.end(); ++it) assert(runmax<it->low or runmin>it->up);
  }

  _jecname = jp::jecgt + id + (jp::ismc ? jp::jecversmc : jp::jecversdt) + (jp::ismc ? "_MC" : "_DATA");
  _jecpath = jp::useJECDB ? (jp::JECDBLoc+"/textFiles/"+_jecname+"/") : "CondFormats/JetMETObjects/data/";

  addJECLvl(jec,vpar,"L1FastJet"); // L1Offset for others than AK*PF
  addJECLvl(jec,vpar,"L2Relative");
  addJECLvl(jec,vpar,"L3Absolute");

  if (jp::isdt) {
    if (!jp::skipl2l3res) addJECLvl(jec,vpar,"L2L3Residual");
    addJECLvl(jec,vpar,"Uncertainty",true);
  }
  jec.corr = new FactorizedJetCorrector(vpar);

  // For type-I and type-II MET
  vrc.push_back(JetCorrectorParameters(string(_jecpath+_jecname+"_L1RC_"+jp::algo+".txt").c_str()));
  jec.l1rc = new FactorizedJetCorrector(vrc);

  _jecs.push_back(jec);
} // addJRC

bool IOV::setJEC(FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc, int run) {
  assert(_jecs.size()!=0);

  if (_current != -1) { // Check if we are already in the correct IOV
    auto cit = &_jecs[_current];
    if (cit->low <= run and cit->up >= run)
      return true;
  }

  if (jp::ismc or run==-1) {
    cout << endl << "Loading a single JEC" << endl;
    *corr = _jecs[0].corr;
    *l1rc = _jecs[0].l1rc;
    if (jp::isdt) *unc = _jecs[0].unc;
    for (auto &name: _jecs[0].names)
      cout << "Loading ... " << name << endl;
    return true;
  } else {
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
  } // setJEC
  cout << "IOV for run " << run << " not found!!" << endl << flush;
  return false;
}

#endif
