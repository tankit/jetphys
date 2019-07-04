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

#ifndef histosInfo_h
#define histosInfo_h

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
#include "tools.h"
#include "IOV.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

class histosInfo {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  static const Int_t kMaxGenJets_   = 100;
  static const Int_t kMaxPFJetsCHS_ = 100;

  // Declaration of leaf types
  //QCDEvent        *events;
  vector<vector<int> > filterIdList_;
  Int_t           EvtHdr__mRun;
  UInt_t          EvtHdr__mEvent; // Int_t -> UInt_t
  Int_t           EvtHdr__mLumi;
  Int_t           EvtHdr__mBunch;
  Float_t         EvtHdr__mTrPu;
  Float_t         EvtHdr__mWeight;
  Float_t         EvtHdr__mPFRho;
  Float_t         PFMet__et_;
  Float_t         PFMet__sumEt_;
  Float_t         PFMet__phi_;
  vector<int>     TriggerDecision_;
  vector<string>  triggerList_;
  vector<int>     L1Prescale_;
  vector<int>     HLTPrescale_;
  Int_t           PFJetsCHS__;
  Double_t        PFJetsCHS__P4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__P4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__unc_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__looseID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

  // List of branches
  TBranch        *b_events_EvtHdr__mRun;   //!
  TBranch        *b_events_EvtHdr__mEvent;   //!
  TBranch        *b_events_EvtHdr__mLumi;   //!
  TBranch        *b_events_EvtHdr__mBunch;   //!
  TBranch        *b_events_EvtHdr__mTrPu;   //!
  TBranch        *b_events_EvtHdr__mWeight;   //!
  TBranch        *b_events_EvtHdr__mPFRho;   //!
  TBranch        *b_events_PFMet__et_;   //!
  TBranch        *b_events_PFMet__sumEt_;   //!
  TBranch        *b_events_PFMet__phi_;   //!
  TBranch        *b_events_TriggerDecision_;   //!
  TBranch        *b_events_triggerList_;   //!
  TBranch        *b_events_L1Prescale_;   //!
  TBranch        *b_events_HLTPrescale_;   //!
  TBranch        *b_events_PFJetsCHS__;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fX;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fY;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fT;   //!
  TBranch        *b_PFJetsCHS__cor_;   //!
  TBranch        *b_PFJetsCHS__unc_;   //!
  TBranch        *b_PFJetsCHS__area_;   //!
  TBranch        *b_PFJetsCHS__looseID_;   //!
  TBranch        *b_PFJetsCHS__tightID_;   //!

  /////////////////////////////////////////////////////////////////////////////
  // Following lines added by hand and must come *after* auto-generated header
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // start manually added code
  /////////////////////////////////////////////////////////////////////////////

  // Auxiliary variables (declare after TTree variables only)
  Float_t         pthat;
  Float_t         weight;
  Int_t           run;
  UInt_t          evt;
  Int_t           lbn;

  Float_t         trpu;
  Int_t           itpu;
  Int_t           ootpulate;
  Int_t           ootpuearly;

  Int_t           npv;
  Int_t           npvgood;
  Float_t         pvx;
  Float_t         pvy;
  Float_t         pvz;
  Float_t         pvrho;
  //Float_t         pvchisq;
  Float_t         pvndof;
  //UInt_t          pvntrk;
  Float_t         bsx;
  Float_t         bsy;
  //
  Int_t           njt;
  Int_t           jt3leads[3]; // The three leading jets
  static const Int_t _maxnjt = 250;//100;//90;//kMaxPFJets_;
  Double_t        *jtp4x;//[_maxnjt];   //[njt]
  Double_t        *jtp4y;//[_maxnjt];   //[njt]
  Double_t        *jtp4z;//[_maxnjt];   //[njt]
  Double_t        *jtp4t;//[_maxnjt];   //[njt]
  Float_t         jte[_maxnjt];   //[njt]
  Float_t         jtpt[_maxnjt];   //[njt]
  //Float_t         jtptuchs[_maxnjt];   //EXTRA
  Float_t         jtptu[_maxnjt];   //EXTRA
  //Float_t         jteuchs[_maxnjt];   //EXTRA
  Float_t         jteu[_maxnjt];   //EXTRA
  Float_t         jteta[_maxnjt];   //[njt]
  Float_t         jtphi[_maxnjt];   //[njt]
  Float_t         jty[_maxnjt];   //[njt]
  Float_t         *jta;//[_maxnjt];   //[njt]
  Float_t         *jtjes;//[_maxnjt];   //[njt]
  Float_t         *jtbeta;
  Float_t         *jtbetastar;
  Float_t         *jtbetaprime;
  Float_t         jtjesnew[_maxnjt];   //[njt]
  Float_t         jtjes_l1[_maxnjt];   //[njt]
  Float_t         jtjes_l2l3[_maxnjt];   //[njt]
  Float_t         jtjes_res[_maxnjt];   //[njt]
  Bool_t          *jtidloose;//[_maxnjt];   //[njt]
  Bool_t          *jtidtight;//[_maxnjt];   //[njt]

  Float_t         rho;
  Float_t         met;
  Float_t         metphi;
  Float_t         metsumet;

  Float_t         met1;
  Float_t         metphi1;
  Float_t         met2;
  Float_t         metphi2;

  /////////////////////////////////////////////////////////////////////////////
  // end manually added code
  /////////////////////////////////////////////////////////////////////////////

  histosInfo();
  histosInfo(TChain *tree); // custom
  virtual ~histosInfo();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree); // custom
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1); // custom

  /////////////////////////////////////////////////////////////////////////////
  // start more manually added code
  /////////////////////////////////////////////////////////////////////////////

  // Additional code for processing
  ofstream *ferr;
  TFile *_outfile;
  map<string, TH1D*> _pudist;
  TH1D *_pumc;
  //TH1F *pudt;
  set<string> _trigs;
  TH1D *hmcweight;

  vector<string> _availTrigs;
  vector<unsigned int> _goodTrigs;
  vector<double> _goodWgts;
  vector<int> _goodNos;
  vector<int> _goodVNos;

  bool LoadJSON(const char* filename);
  bool LoadLumi(const char* filename);
  bool LoadPuProfiles(const char* datafile, const char* mcfile);

  bool GetTriggers();

  inline double DPhi(double phi1, double phi2) { // Return value between 0 and phi.
    double dphi = fabs(phi1 - phi2);
    return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
  }

private:
  Long64_t _entries;
  double _xsecMinBias;
  double _w, _w0;
  map<string, double> _wt; // Mainly for trigger pileup weigths

  vector<int> _jkmore;

  TLorentzVector p4, gp4, genp4;
  IOV _iov;
  FactorizedJetCorrector *_JEC, *_L1RC;
  JetCorrectionUncertainty *_jecUnc;
  //L3Corr *_jecUnc2;
  double _pthatweight;

  bool   _pass;
  bool   _pass_qcdmet;

  map<int, map<int, int> > _json;
  map<int, map<int, float> > _avgpu;
  map<int, map<int, float> > _lums;
  map<int, map<int, float> > _lums2;
  map<string, map<int, int> > _prescales;
  TH2F *_ecalveto;

  map<int, map<int, set<int> > > _duplicates;
  set<int> _badruns;
  set<pair<int, int> > _badlums;
  set<pair<int, int> > _nolums;
  set<pair<int, int> > _badjson;
  set<pair<int, int> > _passedlumis;
};

#endif

#ifdef histosInfo_cxx
histosInfo::histosInfo()
{
  cout << "This is not supposed to happen! Provide a tree" << endl;
}

histosInfo::histosInfo(TChain *tree) : fChain(0)
{
  assert(tree);
  Init(tree); // custom

  // Reset some pointers
  _outfile = 0;
  _outfile = new TFile("pileups.root","RECREATE");
}

histosInfo::~histosInfo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t histosInfo::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// Check that the correct tree is open in the chain
Long64_t histosInfo::LoadTree(Long64_t entry)
{
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;

  // A new tree is opened
  if (fChain->GetTreeNumber() != fCurrent) {

    fCurrent = fChain->GetTreeNumber();
    *ferr << "Opening tree number " << fChain->GetTreeNumber() << endl;

    if (jp::isdt) {
      // Reload the triggers and print them
      if (!GetTriggers()) {
        *ferr << "Failed to load DT triggers. Check that the SMPJ tuple has the required histograms. Aborting..." << endl;
        return -4;
      }
      *ferr << "Tree " << fCurrent << " triggers:" << endl;
      for (auto trigi = 0u; trigi < _availTrigs.size(); ++trigi) {
        auto str = _availTrigs[trigi];
        *ferr << str;
        auto trgplace = std::find(_goodTrigs.begin(),_goodTrigs.end(),trigi);
        if (trgplace!=_goodTrigs.end()) *ferr << "_y (" << _goodWgts[trgplace-_goodTrigs.begin()] << ") ";
        else *ferr << "_n ";
        if (trigi%(jp::notrigs+1)==jp::notrigs) *ferr << endl;
      }
      *ferr << endl << flush;
    } else if (jp::pthatbins) {
      TString filename = fChain->GetCurrentFile()->GetName();
      // Check the position of the current file in the list of file names
      unsigned currFile = std::find_if(jp::pthatfiles.begin(),jp::pthatfiles.end(),[&filename] (string s) { return filename.Contains(s); })-jp::pthatfiles.begin();
      if (jp::pthatnevts[currFile]<=0.0 or jp::pthatsigmas[currFile]<=0) {
        *ferr << "Suspicious pthat slice information for file " << jp::pthatfiles[currFile] << ". Aborting..." << endl;
        return -3;
      }
      _pthatweight = jp::pthatsigmas[currFile]/jp::pthatnevts[currFile];
      _pthatweight /= jp::pthatsigmas[jp::npthatbins-1]/jp::pthatnevts[jp::npthatbins-1]; // Normalize
      *ferr << "Pthat bin changing." << endl << "File " << currFile << " " << fChain->GetCurrentFile()->GetName();
      *ferr << " should correspond to the range [" << jp::pthatranges[currFile] << "," << jp::pthatranges[currFile+1] << "]";
      *ferr << endl << "Weight: " << _pthatweight << endl;
    }
    // slices with pthat bins
  }
  return centry;
}

void histosInfo::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
  fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
  fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
  fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch, &b_events_EvtHdr__mBunch);
  fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu, &b_events_EvtHdr__mTrPu);
  fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);
  fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
  fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);
  fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_, &b_events_PFMet__phi_);
  fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
  fChain->SetBranchAddress("triggerList_", &triggerList_, &b_events_triggerList_);
  fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
  fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
}

Bool_t histosInfo::Notify()
{
  return kTRUE;
}

void histosInfo::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t histosInfo::Cut(Long64_t entry)
{
  return 1;
}
#endif // #ifdef HistosFill_cxx
