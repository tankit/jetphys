// JEC intervals of validity
#ifndef JEC_IOV
#define JEC_IOV
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <map>
#include <string>
#include <iostream>
#include <utility>

using namespace std;

namespace jec {
  
  // header part
  class IOV {

  public:
     IOV(string algo = "AK5PF") : _algo(algo) {};
    ~IOV() {}; // leaks memory, fix later
    void add(string name, int runmin, int runmax, bool isdata = true);
    FactorizedJetCorrector *get(int run);

  private:
    typedef map<int, pair<int, FactorizedJetCorrector*> > IOVmap;
    string _algo;
    IOVmap jecs;
    
  }; // class IOV

  // body part (separate later to another file)
  void IOV::add(string id, int runmin, int runmax, bool isdata) {

    // sanity checks to avoid IOV overlaps
    assert(runmax>=runmin);
    assert(jecs.find(runmin)==jecs.end());
    for (IOVmap::const_iterator it = jecs.begin(); it != jecs.end(); ++it) {
      if (it->first < runmin) assert(it->second.first < runmin);
      if (it->first > runmin) assert(it->first > runmax);
    }
    
    const char *s;
    const char *p = "CondFormats/JetMETObjects/data/";
    const char *t = "GR_R_42_V23_";
    const char *a = _algo.c_str();

    // L1FastJet for AK*PF, L1Offset for others
    s = Form("%s%sL1FastJet_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1 = new JetCorrectorParameters(s);
      //(_algo=="AK5PF" || _algo=="AK7PF" ?
      //new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1FastJet_%s.txt",a)) :
      //new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L1Offset_%s.txt",a)));
      //
    s = Form("%s%sL2Relative_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2 = new JetCorrectorParameters(s);
    s = Form("%s%sL3Absolute_%s.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l3 = new JetCorrectorParameters(s);
    //JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_L2L3Residual_%s.txt",a));
    s = Form("%s%sL2L3Residual_%s.txt",p,t,a); cout<<s<<endl<<flush;
    // Switched off IOV handling for now
      //Form("CondFormats/JetMETObjects/data/Jec_V14_%s_L2L3Residual_%s.txt",id.c_str(),a));
    JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> vpar;
    vpar.push_back(*par_l1);
    vpar.push_back(*par_l2);
    vpar.push_back(*par_l3);
    if (isdata) vpar.push_back(*par_l2l3res);
    FactorizedJetCorrector *jec = new FactorizedJetCorrector(vpar);
    
    jecs[runmin] = pair<int, FactorizedJetCorrector*>(runmax, jec);
  } // add
  
  FactorizedJetCorrector *IOV::get(int run) {
    assert(jecs.size()!=0);
    for (IOVmap::const_iterator it = jecs.begin(); it != jecs.end(); ++it) {
      if (it->first <= run && run <= it->second.first)
	return it->second.second;
    }
    cout << "IOV for run " << run << " not found!!" << endl << flush;
    assert(false);
    return 0;
  } // get


} // namescape jec
#endif
