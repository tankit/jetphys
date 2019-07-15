// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: April 19, 2010

////////////////////////////////////////////////////////////////////////
// Notes:   Automatically created using TChain::MakeClass("HistosFill")
//          Keep variable declarations in the automatic order,
//          update array maximum sizes!!
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  4 11:57:12 2016 by ROOT version 6.06/00
// from TTree ProcessedTree/ProcessedTree
//////////////////////////////////////////////////////////

#ifndef HistosFill_h
#define HistosFill_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "settings.h"
#include "HistosBasic.h"
#include "HistosEta.h"
#include "HistosMC.h"
#include "HistosRun.h"
#include "HistosAll.h"
#include "tools.h"
#include "IOV.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#if JETRESO == 1
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#endif

class HistosFill {
public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  static const Int_t kMaxGenJets_   = 100;
  static const Int_t kMaxPFJetsCHS_ = 100;

  // Declaration of leaf types
  Bool_t          EvtHdr__mIsPVgood;
  Int_t           EvtHdr__mRun;
  Long64_t        EvtHdr__mEvent;
  Int_t           EvtHdr__mLumi;
  Int_t           EvtHdr__mBunch;
  Int_t           EvtHdr__mNVtx;
  Int_t           EvtHdr__mNVtxGood;
  Int_t           EvtHdr__mOOTPUEarly;
  Int_t           EvtHdr__mOOTPULate;
  Int_t           EvtHdr__mINTPU;
  Int_t           EvtHdr__mNBX;
  Float_t         EvtHdr__mPVndof;
  Float_t         EvtHdr__mTrPu;
  Float_t         EvtHdr__mPVx;
  Float_t         EvtHdr__mPVy;
  Float_t         EvtHdr__mPVz;
  Float_t         EvtHdr__mBSx;
  Float_t         EvtHdr__mBSy;
  Float_t         EvtHdr__mBSz;
  Float_t         EvtHdr__mPthat;
  Float_t         EvtHdr__mWeight;
  Float_t         EvtHdr__mCaloRho;
  Float_t         EvtHdr__mPFRho;
  Float_t         PFMet__et_;
  Float_t         PFMet__sumEt_;
  Float_t         PFMet__phi_;
#ifdef NEWMODE
  Float_t         PFMetT0__et_;
  Float_t         PFMetT0__sumEt_;
  Float_t         PFMetT0__phi_;
  Float_t         PFMetT0T1__et_;
  Float_t         PFMetT0T1__sumEt_;
  Float_t         PFMetT0T1__phi_;
#endif
  vector<int>     FilterDecision_;
  vector<int>     TriggerDecision_;
  vector<int>     L1Prescale_;
  vector<int>     HLTPrescale_;
  Int_t           GenJets__;
  Double_t        GenJets__fCoordinates_fX[kMaxGenJets_];   //[GenJets__]
  Double_t        GenJets__fCoordinates_fY[kMaxGenJets_];   //[GenJets__]
  Double_t        GenJets__fCoordinates_fZ[kMaxGenJets_];   //[GenJets__]
  Double_t        GenJets__fCoordinates_fT[kMaxGenJets_];   //[GenJets__]
  Int_t           PFJetsCHS__;
  Double_t        PFJetsCHS__P4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
#ifdef NEWMODE
  Int_t           PFJetsCHS__genIdx_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
#else
  Double_t        PFJetsCHS__genP4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
#endif
  Float_t         PFJetsCHS__genR_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__looseID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
#ifdef NEWMODE
  Float_t         PFJetsCHS__QGL_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__QGAx2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__QGMul_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__QGPtD_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__partonFlavourPhysicsDef_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  vector<int>     genFlavour_;
  vector<int>     genFlavourPartonPhysicsDef_;
  vector<int>     genFlavourHadron_;
#else
  Float_t         PFJetsCHS__QGtagger_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  vector<float>   genFlavour_;
  vector<float>   genFlavourHadron_;
#endif
  Float_t         PFJetsCHS__chf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__nhf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__nemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__cemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__muf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__hf_hf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__hf_phf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__hf_hm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__hf_phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__chm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__nhm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__elm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mum_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__ncand_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__betaPrime_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mpuTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mlvTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mjtTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

  /////////////////////////////////////////////////////////////////////////////
  // Following lines added by hand and must come *after* auto-generated header
  /////////////////////////////////////////////////////////////////////////////

  Int_t GetEntry(Long64_t entry) {
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
  }
  virtual void     Show(Long64_t entry = -1) {
    if (!fChain) return;
    fChain->Show(entry);
  }
  virtual Long64_t LoadTree(Long64_t entry);
  virtual Int_t    FindMCSliceIdx(string filename);
  virtual void     PrintInfo(string info, bool printcout = false);
  virtual void     PrintMemInfo(bool printcout = false);

  // We don't delete that much stuff here, since ROOT takes care of garbage collection (and gets very easily angry!!!)
  ~HistosFill() {}
  HistosFill(TChain *tree);
  virtual bool     Init(TChain *tree); // custom

  virtual void     Loop();
  virtual bool     PreRun();
  virtual bool     AcceptEvent();
  virtual void     Report();

  void             FillJetID(vector<bool> &id);
  bool             GetFilters();
  bool             GetTriggers();

  bool             LoadJSON();
  bool             LoadLumi();
  bool             LoadPuProfiles();

  void             InitBasic(string name);
  void             FillBasic(string name);
  void             FillSingleBasic(HistosBasic *h);
  void             WriteBasic();

  void             InitEta(string name);
  void             FillEta(string name, Float_t *pt, Float_t *eta, Float_t *phi);
  void             FillSingleEta(HistosEta *h, Float_t *pt, Float_t *eta, Float_t *phi);
  void             WriteEta();

  void             InitMC(string name);
  void             FillMC(string name, Float_t *recopt, Float_t *genpt, Float_t *pt, Float_t *eta, Float_t *phi);
  void             FillSingleMC(HistosMC *h, Float_t *recopt, Float_t *genpt, Float_t *pt, Float_t *eta, Float_t *phi);
  void             WriteMC();

  void             InitRun(string name, double etamin, double etamax);
  void             FillRun(string name);
  void             WriteRun();

  void             InitAll(string name);
  void             FillAll(string name);
  void             WriteAll();

  inline double DPhi(double phi1, double phi2) { // Return value between 0 and phi.
    double dphi = fabs(phi1 - phi2);
    return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
  }

private:
  bool   _pass;
  bool   _pass_qcdmet;
  bool   _initsuccess;

  int _eraIdx;
  int _nbadevts_dup;
  int _nbadevts_run;
  int _nbadevts_ls;
  int _nbadevts_lum;
  int _nbadevts_json;
  int _nbadevts_veto;
  int _bscounter_bad;
  int _bscounter_good;
  int _rhocounter_bad;
  int _rhocounter_good;
  int _trgcounter;
  int _evtcounter;
  int _totcounter;
  Long64_t _jentry;
  Long64_t _nentries;
  Long64_t _ntot ;
  Long64_t _nskip;;
  double _xsecMinBias;
  double _w, _w0;
  double _binnedmcweight;
  int    _binnedmcrepeats;

  vector<string> _availFlts;
  vector<string> _availTrigs;
  vector<unsigned int> _goodTrigs;
  vector<double> _goodWgts;
  vector<bool> _jetids;
  vector<int> _jkmore;

  set<string> _trigs;
  set<int> _badruns;
  set<pair<int, int> > _badlums;
  set<pair<int, int> > _nolums;
  set<pair<int, int> > _badjson;
  set<pair<int, int> > _passedlumis;

  map<int, map<int, set<int> > > _duplicates;
  map<string, map<int, int> > _prescales;
  map<int, map<int, int> > _json;
  map<int, map<int, float> > _avgpu;
  map<int, map<int, float> > _lums;
  map<int, map<int, float> > _lums2;
  map<string, double> _wt; // Trigger pileup and trigger weights
  map<string, vector<HistosBasic*> > _histos;
  map<string, vector<HistosEta*> > _etahistos;
  map<string, vector<HistosMC*> > _mchistos;
  map<string, HistosRun*> _runhistos;
  map<string, HistosAll*> _allhistos;
  map<string, TH1D*> _pudist;
  map<string, int> _cnt; // efficiency counters

  TH1D *_pumc;
  TH1D *hmcweight;
  TH2D *h2HotExcl;
  TH2F *_h2mu = 0;
  TH3D *_h3probg;

  TFile *_outfile;
  TFile *fHotExcl;

  TLorentzVector p4, gp4, genp4, _j1, _j2;
  IOV _iov;
  FactorizedJetCorrector *_JEC, *_L1RC;
  JetCorrectionUncertainty *_jecUnc;

  ofstream *ferr;
  MemInfo_t _info;

  /////////////////////////////////////////////////////////////////////////////
  // Auxiliary variables (declare after TTree variables only)
  /////////////////////////////////////////////////////////////////////////////
  Float_t         trpu;
  Float_t         pvrho;

  Float_t         &pthat;
  Float_t         &weight;
  Int_t           &run;
  Long64_t        &evt;
  Int_t           &lbn;

  Int_t           &itpu;
  Int_t           &ootpulate;
  Int_t           &ootpuearly;

  Int_t           &npv;
  Int_t           &npvgood;
  Float_t         &pvx;
  Float_t         &pvy;
  Float_t         &pvz;
  Float_t         &pvndof;
  Float_t         &bsx;
  Float_t         &bsy;
  //
  Int_t           &njt;
  Int_t           jt3leads[3]; // The three leading jets
  static const Int_t _maxnjt = 250;//100;//90;//kMaxPFJets_;
  Double_t        *jtp4x;//[_maxnjt];   //[njt]
  Double_t        *jtp4y;//[_maxnjt];   //[njt]
  Double_t        *jtp4z;//[_maxnjt];   //[njt]
  Double_t        *jtp4t;//[_maxnjt];   //[njt]
  Float_t         jte[_maxnjt];   //[njt]
  Float_t         jtpt[_maxnjt];   //[njt]
  Float_t         jtptu[_maxnjt];   //EXTRA
  Float_t         jteu[_maxnjt];   //EXTRA
  Float_t         jteta[_maxnjt];   //[njt]
  Float_t         jtphi[_maxnjt];   //[njt]
  Float_t         jty[_maxnjt];   //[njt]
  Float_t         jtpt_nol2l3[_maxnjt];   //[njt]
  Float_t         *jta;//[_maxnjt];   //[njt]
  Float_t         *jtjes;//[_maxnjt];   //[njt]
  Float_t         *jtbetaprime;//[_maxnjt];   //[njt]
  Float_t         jtjesnew[_maxnjt];   //[njt]
  Float_t         jtjes_l1[_maxnjt];   //[njt]
  Float_t         jtjes_l2l3[_maxnjt];   //[njt]
  Float_t         jtjes_res[_maxnjt];   //[njt]
  Bool_t          *jtidloose;//[_maxnjt];   //[njt]
  Bool_t          *jtidtight;//[_maxnjt];   //[njt]

  Float_t         *jtgenr;//[_maxnjt];   //[njt]
  Double_t        *jtgenp4x;//[_maxnjt];   //[njt]
  Double_t        *jtgenp4y;//[_maxnjt];   //[njt]
  Double_t        *jtgenp4z;//[_maxnjt];   //[njt]
  Double_t        *jtgenp4t;//[_maxnjt];   //[njt]
  Int_t           *jtgenidx;//[_maxnjt];   //[njt]
  Float_t         jtgenpt[_maxnjt];   //[njt]
  Float_t         jtgeny[_maxnjt];   //EXTRA
  Float_t         jtgeneta[_maxnjt];   //[njt]
  Float_t         jtgenphi[_maxnjt];   //[njt]

  Float_t         *qgl;//[_maxnjt];   //[njt]
  Int_t           *partonflavor;//[_maxnjt];   //[njt]
#ifdef NEWMODE
  Int_t           *partonflavorphys;
#endif

  Int_t           *jtn;//[_maxnjt];   //[njt]
  Int_t           *jtnch;//[_maxnjt];   //[njt]
  Int_t           *jtnnh;//[_maxnjt];   //[njt]
  Int_t           *jtnne;//[_maxnjt];   //[njt]
  Int_t           *jtnce;//[_maxnjt];   //[njt]
  Int_t           *jtnmu;//[_maxnjt];   //[njt]
  Int_t           *jtnhh;//[_maxnjt];   //[njt]
  Int_t           *jtnhe;//[_maxnjt];   //[njt]

  Float_t         *jtchf;//[_maxnjt];   //[njt]
  Float_t         *jtnhf;//[_maxnjt];   //[njt]
  Float_t         *jtnef;//[_maxnjt];   //[njt]
  Float_t         *jtcef;//[_maxnjt];   //[njt]
  Float_t         *jtmuf;//[_maxnjt];   //[njt]
  Float_t         *jthhf;//[_maxnjt];   //[njt]
  Float_t         *jthef;//[_maxnjt];   //[njt]

  Int_t           &gen_njt;
  Double_t        *gen_jtp4x;//[_maxnjt];   //[njt]
  Double_t        *gen_jtp4y;//[_maxnjt];   //[njt]
  Double_t        *gen_jtp4z;//[_maxnjt];   //[njt]
  Double_t        *gen_jtp4t;//[_maxnjt];   //[njt]
  Float_t         gen_jtpt[_maxnjt];   //[njt]
  Float_t         gen_jteta[_maxnjt];   //[njt]
  Float_t         gen_jtphi[_maxnjt];   //[njt]
  Float_t         gen_jty[_maxnjt];   //[njt]
  Int_t           gen_partonflavor[_maxnjt];

  Float_t         &rho;
  Float_t         &met;
  Float_t         &metphi;
  Float_t         &metsumet;
#ifdef NEWMODE
  Float_t         &met01;
  Float_t         &metsumet01;
#endif

  Float_t         met1;
  Float_t         metphi1;
  Float_t         met1_nol2l3;
  Float_t         metphi1_nol2l3;
  Float_t         met2;
  Float_t         metphi2;

  string         _runinfo;
};

#endif
