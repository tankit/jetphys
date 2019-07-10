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
#include <random>

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

#if JETRESO == 1
struct JERIOV {
  vector<string> names;
  JME::JetResolution *Res;
  JME::JetResolutionScaleFactor *SF;
};
#endif

// header part
class IOV {
public:
  IOV();
  ~IOV();
  void addJEC(string id, int runmin = -1, int runmax = -1);
  bool setJEC(FactorizedJetCorrector**,FactorizedJetCorrector**,JetCorrectionUncertainty** unc=0,int run = -1);
#if JETRESO == 1
  double getJERCF(double jPt, double jEta, double jE, double rho, double jPtGen = -4);
#endif

private:
  void addJECLvl(JECIOV &jec, vector<JetCorrectorParameters> &vpar, string tag, bool unc);
  vector<JECIOV> _jecs;
  string _jecname;
  string _jecpath;
  int _currentJEC;
#if JETRESO == 1
  void addJERLvl(bool isSF=true);
  double getJER(double jPt, double jEta, double rho);
  JERIOV _jer;
  string _jername;
  string _jerpath;
  std::uint32_t _seed;
  std::mt19937 _mersennetwister;
  static constexpr const double MIN_JET_ENERGY = 0.01;
#endif
}; // class IOV

IOV::IOV() {
  _currentJEC = -1;
#if JETRESO == 1
  // At the value of _seed: the old question - should the seed of a rng be random itself?
  // Here we prefer stability, but the user can vary the seed if necessary. Moreover, https://xkcd.com/221/
  _seed = 4;
  _mersennetwister = std::mt19937(_seed);

  unsigned pos = jp::JERpos;
  if (pos>=jp::JERIOVnames.size()) pos = 0;
  _jername = jp::jergt + jp::JERIOVnames[pos] + jp::jervers + (jp::ismc ? "_MC" : "_DATA");
  _jerpath = jp::JERDBLoc+"/textFiles/"+_jername+"/";
  addJERLvl(false);
  addJERLvl(true);

#endif
}
IOV::~IOV() {}

void IOV::addJECLvl(JECIOV &jec, vector<JetCorrectorParameters> &vpar, string tag, bool unc = false) {
  string path = _jecpath+_jecname+"_"+tag+"_"+jp::algo+".txt";
  if (access(path.c_str(), F_OK)==-1) {
    cout << "The JEC IOV file " << path << " was not found" << endl;
  } else {
    jec.names.push_back(path);
    if (unc) jec.unc = new JetCorrectionUncertainty(path);
    else vpar.push_back(JetCorrectorParameters(path));
  }
}

#if JETRESO == 1
void IOV::addJERLvl(bool isSF) {
  string path = _jerpath+_jername+"_"+(isSF ? "SF" : "PtResolution")+"_"+jp::algo+".txt";
  if (access(path.c_str(), F_OK)==-1) {
    cout << "The JER IOV file " << path << " was not found" << endl;
  } else {
    _jer.names.push_back(path);
    cout << "Loading ... " << path << endl;
    if (isSF) _jer.SF  = new JME::JetResolutionScaleFactor(path);
    else _jer.Res = new JME::JetResolution(path);
  }
}
#endif

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

bool IOV::setJEC(FactorizedJetCorrector** corr,FactorizedJetCorrector** l1rc,JetCorrectionUncertainty** unc, int run) {
  assert(_jecs.size()!=0);

  if (_currentJEC != -1) { // Check if we are already in the correct IOV
    auto cit = &_jecs[_currentJEC];
    if (cit->low <= run and cit->up >= run)  return true;
  }

  if (jp::ismc or run==-1) {
    cout << endl << "Loading a single JEC" << endl;
    *corr = _jecs[0].corr;
    *l1rc = _jecs[0].l1rc;
    if (jp::isdt) *unc = _jecs[0].unc;
    for (auto &name: _jecs[0].names) cout << "Loading ... " << name << endl;
    return true;
  } else {
    // If the current IOV was not fine, search through all the IOVs
    for (int i = 0, N = int(_jecs.size()); i<N; ++i) {
      auto it = &_jecs[i];
      if (it->low <= run and it->up >= run) {
        if (i!=_currentJEC) {
          _currentJEC = i;
          *corr = it->corr;
          *l1rc = it->l1rc;
          *unc = it->unc;
          cout << endl << "IOV handling in use." << endl;
          for (auto &name: it->names)
            cout << "Loading ... " << name << endl;
        }
        return true;
      }
    }
  } // setJEC
  cout << "IOV for run " << run << " not found!!" << endl << flush;
  return false;
}

#if JETRESO == 1
double IOV::getJER(double jPt, double jEta, double rho) {
  return _jer.Res->getResolution({{JME::Binning::JetPt, jPt}, {JME::Binning::JetEta, jEta}, {JME::Binning::Rho, rho}});
}

// The method presented here can be found in https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
// and the corresponding code in https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
double IOV::getJERCF(double jPt, double jEta, double jE, double rho, double jPtGen) {
  double Reso = getJER(jPt,jEta,rho);
  double SF = _jer.SF->getScaleFactor({{JME::Binning::JetEta, jEta}, {JME::Binning::Rho, rho}}, Variation::NOMINAL);

  // Case 0: by default the JER correction factor is equal to 1
  double CF = 1.;
  // We see if the gen jet meets our requirements
  bool condPt = jPtGen>MIN_JET_ENERGY;
  double relDPt = condPt ? (jPt - jPtGen)/jPt : 0.0;
  bool condPtReso = fabs(relDPt) < 3*Reso;
  if (condPt and condPtReso) {
    // Case 1: we have a "good" gen jet matched to the reco jet (indicated by positive gen jet pt)
    CF += (SF - 1.)*relDPt;
  } else if (SF > 1) {
    // Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation
    double sigma = Reso*std::sqrt(SF*SF - 1);
    std::normal_distribution<> d(0, sigma);
    CF += d(_mersennetwister);
  }

  // Negative or too small smearFactor. Safety precautions.
  double CFLimit = MIN_JET_ENERGY / jE;
  if (CF < CFLimit) CF = CFLimit;

  return CF;
}
#endif

#endif
