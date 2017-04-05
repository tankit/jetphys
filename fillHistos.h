// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: April 19, 2010

////////////////////////////////////////////////////////////////////////
// Notes:   Automatically created using TChain::MakeClass("fillHistos")
//          Keep variable declarations in the automatic order,
//          update array maximum sizes!!
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  4 11:57:12 2016 by ROOT version 6.06/00
// from TTree ProcessedTree/ProcessedTree
//////////////////////////////////////////////////////////

#ifndef fillHistos_h
#define fillHistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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

#include <map>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>

using std::map;
using std::vector;
using std::string;
using std::cout;
using std::endl;

#include "settings.h"
#include "basicHistos.h"
#include "etaHistos.h"
#include "runHistos.h"
#include "tools.h"

#include "IOV.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Header file for the classes stored in the TTree if any.
//#include "Math/GenVector/PxPyPzE4D.h"

class fillHistos {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  static const Int_t kMaxGenJets_   = 100;
  static const Int_t kMaxPFJetsCHS_ = 100;

  // Declaration of leaf types
  //QCDEvent        *events;
  vector<vector<int> > filterIdList_;
  Bool_t          EvtHdr__mIsPVgood;
  Bool_t          EvtHdr__mHCALNoise;
  Bool_t          EvtHdr__mHCALNoiseNoMinZ;
  Int_t           EvtHdr__mRun;
  UInt_t          EvtHdr__mEvent; // Int_t -> UInt_t
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
  Float_t         CaloMet__et_;
  Float_t         CaloMet__CaloMetPt_;
  Float_t         CaloMet__sumEt_;
  Float_t         CaloMet__phi_;
  Float_t         PFMet__et_;
  Float_t         PFMet__CaloMetPt_;
  Float_t         PFMet__sumEt_;
  Float_t         PFMet__phi_;
  Float_t         MvaMet__et_;
  Float_t         MvaMet__CaloMetPt_;
  Float_t         MvaMet__sumEt_;
  Float_t         MvaMet__phi_;
  vector<int>     TriggerDecision_;
  vector<string>  triggerList_;
  vector<int>     L1Prescale_;
  vector<int>     HLTPrescale_;
  //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > HLTObj_;
  //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > L1Obj_;
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
  Double_t        PFJetsCHS__genP4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Double_t        PFJetsCHS__genP4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__genR_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  vector<double>  PFJetsCHS__jecLabels_[kMaxPFJetsCHS_];
  Float_t         PFJetsCHS__unc_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  vector<float>   PFJetsCHS__uncSrc_[kMaxPFJetsCHS_];
  Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__looseID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__CSVpfPositive_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__CSVpfNegative_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  //Float_t         PFJetsCHS__boosted_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__QGtagger_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__recommend1_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__recommend2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__recommend3_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  //Float_t         PFJetsCHS__pfCombinedCvsL_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  //Float_t         PFJetsCHS__pfCombinedCvsB_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
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
  //Int_t           PFJetsCHS__cm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__beta_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__betaStar_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mpuTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mlvTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Int_t           PFJetsCHS__mjtTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__hof_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__pujid_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__calojetpt_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  Float_t         PFJetsCHS__calojetef_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
  vector<float>   genFlavour_;
  vector<float>   genFlavourHadron_;

  // List of branches
  TBranch        *b_events_filterIdList_;   //!
  TBranch        *b_events_EvtHdr__mIsPVgood;   //!
  TBranch        *b_events_EvtHdr__mHCALNoise;   //!
  TBranch        *b_events_EvtHdr__mHCALNoiseNoMinZ;   //!
  TBranch        *b_events_EvtHdr__mRun;   //!
  TBranch        *b_events_EvtHdr__mEvent;   //!
  TBranch        *b_events_EvtHdr__mLumi;   //!
  TBranch        *b_events_EvtHdr__mBunch;   //!
  TBranch        *b_events_EvtHdr__mNVtx;   //!
  TBranch        *b_events_EvtHdr__mNVtxGood;   //!
  TBranch        *b_events_EvtHdr__mOOTPUEarly;   //!
  TBranch        *b_events_EvtHdr__mOOTPULate;   //!
  TBranch        *b_events_EvtHdr__mINTPU;   //!
  TBranch        *b_events_EvtHdr__mNBX;   //!
  TBranch        *b_events_EvtHdr__mPVndof;   //!
  TBranch        *b_events_EvtHdr__mTrPu;   //!
  TBranch        *b_events_EvtHdr__mPVx;   //!
  TBranch        *b_events_EvtHdr__mPVy;   //!
  TBranch        *b_events_EvtHdr__mPVz;   //!
  TBranch        *b_events_EvtHdr__mBSx;   //!
  TBranch        *b_events_EvtHdr__mBSy;   //!
  TBranch        *b_events_EvtHdr__mBSz;   //!
  TBranch        *b_events_EvtHdr__mPthat;   //!
  TBranch        *b_events_EvtHdr__mWeight;   //!
  TBranch        *b_events_EvtHdr__mCaloRho;   //!
  TBranch        *b_events_EvtHdr__mPFRho;   //!
  TBranch        *b_events_CaloMet__et_;   //!
  TBranch        *b_events_CaloMet__CaloMetPt_;   //!
  TBranch        *b_events_CaloMet__sumEt_;   //!
  TBranch        *b_events_CaloMet__phi_;   //!
  TBranch        *b_events_PFMet__et_;   //!
  TBranch        *b_events_PFMet__CaloMetPt_;   //!
  TBranch        *b_events_PFMet__sumEt_;   //!
  TBranch        *b_events_PFMet__phi_;   //!
  TBranch        *b_events_MvaMet__et_;   //!
  TBranch        *b_events_MvaMet__CaloMetPt_;   //!
  TBranch        *b_events_MvaMet__sumEt_;   //!
  TBranch        *b_events_MvaMet__phi_;   //!
  TBranch        *b_events_TriggerDecision_;   //!
  TBranch        *b_events_triggerList_;   //!
  TBranch        *b_events_L1Prescale_;   //!
  TBranch        *b_events_HLTPrescale_;   //!
  TBranch        *b_events_GenJets__;   //!
  TBranch        *b_GenJets__fCoordinates_fX;   //!
  TBranch        *b_GenJets__fCoordinates_fY;   //!
  TBranch        *b_GenJets__fCoordinates_fZ;   //!
  TBranch        *b_GenJets__fCoordinates_fT;   //!
  TBranch        *b_events_PFJetsCHS__;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fX;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fY;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ;   //!
  TBranch        *b_PFJetsCHS__P4__fCoordinates_fT;   //!
  TBranch        *b_PFJetsCHS__genP4__fCoordinates_fX;   //!
  TBranch        *b_PFJetsCHS__genP4__fCoordinates_fY;   //!
  TBranch        *b_PFJetsCHS__genP4__fCoordinates_fZ;   //!
  TBranch        *b_PFJetsCHS__genP4__fCoordinates_fT;   //!
  TBranch        *b_PFJetsCHS__genR_;   //!
  TBranch        *b_PFJetsCHS__cor_;   //!
  TBranch        *b_PFJetsCHS__jecLabels_;   //!
  TBranch        *b_PFJetsCHS__unc_;   //!
  TBranch        *b_PFJetsCHS__uncSrc_;   //!
  TBranch        *b_PFJetsCHS__area_;   //!
  TBranch        *b_PFJetsCHS__looseID_;   //!
  TBranch        *b_PFJetsCHS__tightID_;   //!
  TBranch        *b_PFJetsCHS__CSVpfPositive_;   //!
  TBranch        *b_PFJetsCHS__CSVpfNegative_;   //!
  //TBranch        *b_PFJetsCHS__boosted_;   //!
  TBranch        *b_PFJetsCHS__QGtagger_;   //!
  TBranch        *b_PFJetsCHS__partonFlavour_;   //!
  TBranch        *b_PFJetsCHS__hadronFlavour_;   //!
  TBranch        *b_PFJetsCHS__recommend1_;   //!
  TBranch        *b_PFJetsCHS__recommend2_;   //!
  TBranch        *b_PFJetsCHS__recommend3_;   //!
  //TBranch        *b_PFJetsCHS__pfCombinedCvsL_;   //!
  //TBranch        *b_PFJetsCHS__pfCombinedCvsB_;   //!
  TBranch        *b_PFJetsCHS__chf_;   //!
  TBranch        *b_PFJetsCHS__nhf_;   //!
  TBranch        *b_PFJetsCHS__nemf_;   //!
  TBranch        *b_PFJetsCHS__cemf_;   //!
  TBranch        *b_PFJetsCHS__muf_;   //!
  TBranch        *b_PFJetsCHS__hf_hf_;   //!
  TBranch        *b_PFJetsCHS__hf_phf_;   //!
  TBranch        *b_PFJetsCHS__hf_hm_;   //!
  TBranch        *b_PFJetsCHS__hf_phm_;   //!
  TBranch        *b_PFJetsCHS__chm_;   //!
  TBranch        *b_PFJetsCHS__nhm_;   //!
  TBranch        *b_PFJetsCHS__phm_;   //!
  TBranch        *b_PFJetsCHS__elm_;   //!
  TBranch        *b_PFJetsCHS__mum_;   //!
  TBranch        *b_PFJetsCHS__ncand_;   //!
  //TBranch        *b_PFJetsCHS__cm_;   //!
  TBranch        *b_PFJetsCHS__beta_;   //!
  TBranch        *b_PFJetsCHS__betaStar_;   //!
  TBranch        *b_PFJetsCHS__mpuTrk_;   //!
  TBranch        *b_PFJetsCHS__mlvTrk_;   //!
  TBranch        *b_PFJetsCHS__mjtTrk_;   //!
  TBranch        *b_PFJetsCHS__hof_;   //!
  TBranch        *b_PFJetsCHS__pujid_;   //!
  TBranch        *b_PFJetsCHS__calojetpt_;   //!
  TBranch        *b_PFJetsCHS__calojetef_;   //!
  TBranch        *b_events_genFlavour_;   //!
  TBranch        *b_events_genFlavourHadron_;   //!

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

  Bool_t          hlt_jet30u;
  Bool_t          hlt_jet60u;
  Bool_t          hlt_jet80u;
  Bool_t          hlt_jet110u;
  Bool_t          hlt_jet150u;
  Bool_t          hlt_jet190u;
  Bool_t          hlt_jet240u;
  Bool_t          hlt_jet300u;
  Bool_t          hlt_jet370u;

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
  static const Int_t _njt = 100;//90;//kMaxPFJets_;
  Double_t        *jtp4x;//[_njt];   //[njt]
  Double_t        *jtp4y;//[_njt];   //[njt]
  Double_t        *jtp4z;//[_njt];   //[njt]
  Double_t        *jtp4t;//[_njt];   //[njt]
  Float_t         jte[_njt];   //[njt]
  Float_t         jtpt[_njt];   //[njt]
  //Float_t         jtptuchs[_njt];   //EXTRA
  Float_t         jtptu[_njt];   //EXTRA
  //Float_t         jteuchs[_njt];   //EXTRA
  Float_t         jteu[_njt];   //EXTRA
  Float_t         jteta[_njt];   //[njt]
  Float_t         jtphi[_njt];   //[njt]
  Float_t         jty[_njt];   //[njt]
  Float_t         *jta;//[_njt];   //[njt]
  Float_t         *jtjes;//[_njt];   //[njt]
  Float_t         *jtbeta;
  Float_t         *jtbetastar;
  Float_t         jtjesnew[_njt];   //[njt]
  Float_t         jtjes_l1[_njt];   //[njt]
  Float_t         jtjes_l2l3[_njt];   //[njt]
  Float_t         jtjes_res[_njt];   //[njt]
  Bool_t          *jtidloose;//[_njt];   //[njt]
  Bool_t          *jtidtight;//[_njt];   //[njt]

  //Short_t         jtgenid[_njt];   //[njt]
  //Short_t         jtgenflv[_njt];   //[njt]
  Float_t         *jtgenr;//[_njt];   //[njt]
  Double_t        *jtgenp4x;//[_njt];   //[njt]
  Double_t        *jtgenp4y;//[_njt];   //[njt]
  Double_t        *jtgenp4z;//[_njt];   //[njt]
  Double_t        *jtgenp4t;//[_njt];   //[njt]
  Float_t         jtgenpt[_njt];   //[njt]
  Float_t         jtgeny[_njt];   //EXTRA

  Int_t           *jtn;//[_njt];   //[njt]
  Int_t           *jtnch;//[_njt];   //[njt]
  Int_t           *jtnnh;//[_njt];   //[njt]
  Int_t           *jtnne;//[_njt];   //[njt]
  Int_t           *jtnce;//[_njt];   //[njt]
  Int_t           *jtnmu;//[_njt];   //[njt]
  Float_t         *jtchf;//[_njt];   //[njt]
  Float_t         *jtnhf;//[_njt];   //[njt]
  Float_t         *jtnef;//[_njt];   //[njt]
  Float_t         *jtcef;//[_njt];   //[njt]
  Float_t         *jtmuf;//[_njt];   //[njt]

  Int_t           gen_njt;
  Double_t        *gen_jtp4x;//[_njt];   //[njt]
  Double_t        *gen_jtp4y;//[_njt];   //[njt]
  Double_t        *gen_jtp4z;//[_njt];   //[njt]
  Double_t        *gen_jtp4t;//[_njt];   //[njt]
  //Float_t         gen_jte[_njt];   //[njt]
  Float_t         gen_jtpt[_njt];   //[njt]
  Float_t         gen_jteta[_njt];   //[njt]
  Float_t         gen_jtphi[_njt];   //[njt]
  Float_t         gen_jty[_njt];   //[njt]

  Float_t         rho;
  Float_t         met;
  Float_t         metphi;
  Float_t         metsumet;

  Float_t         c_rho;
  Float_t         c_met;
  Float_t         c_metsumet;

  Float_t         met1;
  Float_t         metphi1;
  Float_t         met2;
  Float_t         metphi2;

  /////////////////////////////////////////////////////////////////////////////
  // end manually added code
  /////////////////////////////////////////////////////////////////////////////

  fillHistos(TTree *tree=0); // custom
  virtual ~fillHistos();
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
  string _type;
  TFile *_outfile;
  map<string, TH1D*> pudist;
  vector<string> _triggers;
  TH1F *pumc;
  //TH1F *pudt;
  set<string> _trigs;
  vector<bool> _jetids;
  map<string, vector<basicHistos*> > _histos;
  map<string, vector<etaHistos*> > _etahistos;
  map<string, runHistos*> _runhistos;
  TH1D *hmcweight;

  vector<string> _availTrigs;

  void loadJSON(const char* filename);
  void loadLumi(const char* filename);
  void loadPUProfiles(const char* datafile, const char* mcfile);
  void loadPrescales(const char* prescalefilename);

  void loadECALveto(const char* filename);

  void initBasics(string name);
  void fillBasics(string name);
  void fillBasic(basicHistos *h);
  void writeBasics();

  void initEtas(string name);
  void fillEtas(string name);
  void fillEta(etaHistos *h);
  void writeEtas();

  void initRunHistos(string name, double ymin, double ymax);
  void fillRunHistos(string name);
  void writeRunHistos();

  void fillJetID(vector<bool> &id);

  void getTriggers();

  // Move to simpleMath.h later
  inline double delta_phi(double phi1, double phi2) {
    
    double dphi = fabs(phi1 - phi2);
    return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
  }

private:
  bool _mc;
  bool _dt;
  Long64_t _entries;
  Long64_t _entry;
  double _xsecMinBias;
  double _w, _w0;
  map<string, double> _wt;

  vector<int> _jkmore;
  map<Int_t, Int_t> _outlist;

  TLorentzVector p4, gp4, genp4, cp4, pp4;
  jec::IOV *iov;
  FactorizedJetCorrector *_JEC, *_L1RC;
  JetCorrectionUncertainty *_jecUnc;
  //L3Corr *_jecUnc2;

  bool   _pass;
  bool   _evtid;

  set<int> _runveto;
  map<int, map<int, int> > _json;
  map<int, map<int, float> > _avgpu;
  map<int, map<int, float> > _lums;
  map<int, map<int, float> > _lums2;
  //map<string, map<int, int> > _premap;
  map<string, map<int, map<int, int> > > _premap;
  map<string, map<int, int> > _prescales;
  TH2F *ecalveto;

  map<int, map<int, set<int> > > _duplicates;
  set<int> _badruns;
  set<pair<int, int> > _badlums;
  set<pair<int, int> > _nolums;
  set<pair<int, int> > _badjson;
  set<pair<int, int> > _jt15lums;
  int _nbadevts_dup;
  int _nbadevts_run;
  int _nbadevts_ls;
  int _nbadevts_lum;
  int _nbadevts_json;
  int _nbadevts_veto;
  int _nbadevts_stream;
  int _bscounter_bad;
  int _bscounter_good;
  int _halocounter_bad;
  int _halocounter_good;
  int _ecalcounter_bad;
  int _ecalcounter_good;
  int _rhocounter_bad;
  int _rhocounter_good;
  int _trgcounter;
  int _evtcounter;
  int _totcounter;

  TLorentzVector _j1, _j2;

  /////////////////////////////////////////////////////////////////////////////
  // end more manually added code
  /////////////////////////////////////////////////////////////////////////////
};

#endif

#ifdef fillHistos_cxx
fillHistos::fillHistos(TTree *tree) : fChain(0)
{
  _type = _jp_type;

  assert(tree);
  Init(tree); // custom

  // Reset some pointers
  _outfile = 0;

  // ROOT 6 requires this
  Loop();
}

fillHistos::~fillHistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fillHistos::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// Check that the correct tree is open in the chain
Long64_t fillHistos::LoadTree(Long64_t entry)
{
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;

  // If a new tree is opened, reload the triggers and print them
  if (_dt && fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    getTriggers();
    *ferr << "Tree " << fCurrent << " triggers:" << endl;
    for (auto str : _availTrigs) {
      if (str.length()==0)
        *ferr << "x ";
      else
        *ferr << str << " ";
    }
    *ferr << endl << flush;
  }
  return centry;
}

void fillHistos::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  _dt = (_type=="DATA");
  _mc = !_dt;

  fChain->SetBranchAddress("filterIdList_", &filterIdList_, &b_events_filterIdList_);
  fChain->SetBranchAddress("EvtHdr_.mIsPVgood", &EvtHdr__mIsPVgood, &b_events_EvtHdr__mIsPVgood);
  fChain->SetBranchAddress("EvtHdr_.mHCALNoise", &EvtHdr__mHCALNoise, &b_events_EvtHdr__mHCALNoise);
  fChain->SetBranchAddress("EvtHdr_.mHCALNoiseNoMinZ", &EvtHdr__mHCALNoiseNoMinZ, &b_events_EvtHdr__mHCALNoiseNoMinZ);
  fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
  fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
  fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
  fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch, &b_events_EvtHdr__mBunch);
  fChain->SetBranchAddress("EvtHdr_.mNVtx", &EvtHdr__mNVtx, &b_events_EvtHdr__mNVtx);
  fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood, &b_events_EvtHdr__mNVtxGood);
  fChain->SetBranchAddress("EvtHdr_.mOOTPUEarly", &EvtHdr__mOOTPUEarly, &b_events_EvtHdr__mOOTPUEarly);
  fChain->SetBranchAddress("EvtHdr_.mOOTPULate", &EvtHdr__mOOTPULate, &b_events_EvtHdr__mOOTPULate);
  fChain->SetBranchAddress("EvtHdr_.mINTPU", &EvtHdr__mINTPU, &b_events_EvtHdr__mINTPU);
  fChain->SetBranchAddress("EvtHdr_.mNBX", &EvtHdr__mNBX, &b_events_EvtHdr__mNBX);
  fChain->SetBranchAddress("EvtHdr_.mPVndof", &EvtHdr__mPVndof, &b_events_EvtHdr__mPVndof);
  fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu, &b_events_EvtHdr__mTrPu);
  fChain->SetBranchAddress("EvtHdr_.mPVx", &EvtHdr__mPVx, &b_events_EvtHdr__mPVx);
  fChain->SetBranchAddress("EvtHdr_.mPVy", &EvtHdr__mPVy, &b_events_EvtHdr__mPVy);
  fChain->SetBranchAddress("EvtHdr_.mPVz", &EvtHdr__mPVz, &b_events_EvtHdr__mPVz);
  fChain->SetBranchAddress("EvtHdr_.mBSx", &EvtHdr__mBSx, &b_events_EvtHdr__mBSx);
  fChain->SetBranchAddress("EvtHdr_.mBSy", &EvtHdr__mBSy, &b_events_EvtHdr__mBSy);
  fChain->SetBranchAddress("EvtHdr_.mBSz", &EvtHdr__mBSz, &b_events_EvtHdr__mBSz);
  fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat, &b_events_EvtHdr__mPthat);
  fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);
  fChain->SetBranchAddress("EvtHdr_.mCaloRho", &EvtHdr__mCaloRho, &b_events_EvtHdr__mCaloRho);
  fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho, &b_events_EvtHdr__mPFRho);
  fChain->SetBranchAddress("CaloMet_.et_", &CaloMet__et_, &b_events_CaloMet__et_);
  fChain->SetBranchAddress("CaloMet_.CaloMetPt_", &CaloMet__CaloMetPt_, &b_events_CaloMet__CaloMetPt_);
  fChain->SetBranchAddress("CaloMet_.sumEt_", &CaloMet__sumEt_, &b_events_CaloMet__sumEt_);
  fChain->SetBranchAddress("CaloMet_.phi_", &CaloMet__phi_, &b_events_CaloMet__phi_);
  fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
  fChain->SetBranchAddress("PFMet_.CaloMetPt_", &PFMet__CaloMetPt_, &b_events_PFMet__CaloMetPt_);
  fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);
  fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_, &b_events_PFMet__phi_);
  fChain->SetBranchAddress("MvaMet_.et_", &MvaMet__et_, &b_events_MvaMet__et_);
  fChain->SetBranchAddress("MvaMet_.CaloMetPt_", &MvaMet__CaloMetPt_, &b_events_MvaMet__CaloMetPt_);
  fChain->SetBranchAddress("MvaMet_.sumEt_", &MvaMet__sumEt_, &b_events_MvaMet__sumEt_);
  fChain->SetBranchAddress("MvaMet_.phi_", &MvaMet__phi_, &b_events_MvaMet__phi_);
  fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
  fChain->SetBranchAddress("triggerList_", &triggerList_, &b_events_triggerList_);
  fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
  fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
  fChain->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fX", GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fY", GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fT", GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
  fChain->SetBranchAddress("PFJetsCHS_", &PFJetsCHS__, &b_events_PFJetsCHS__);
  fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fX", PFJetsCHS__P4__fCoordinates_fX, &b_PFJetsCHS__P4__fCoordinates_fX);
  fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fY", PFJetsCHS__P4__fCoordinates_fY, &b_PFJetsCHS__P4__fCoordinates_fY);
  fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fZ", PFJetsCHS__P4__fCoordinates_fZ, &b_PFJetsCHS__P4__fCoordinates_fZ);
  fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fT", PFJetsCHS__P4__fCoordinates_fT, &b_PFJetsCHS__P4__fCoordinates_fT);
  fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fX", PFJetsCHS__genP4__fCoordinates_fX, &b_PFJetsCHS__genP4__fCoordinates_fX);
  fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fY", PFJetsCHS__genP4__fCoordinates_fY, &b_PFJetsCHS__genP4__fCoordinates_fY);
  fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fZ", PFJetsCHS__genP4__fCoordinates_fZ, &b_PFJetsCHS__genP4__fCoordinates_fZ);
  fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fT", PFJetsCHS__genP4__fCoordinates_fT, &b_PFJetsCHS__genP4__fCoordinates_fT);
  fChain->SetBranchAddress("PFJetsCHS_.genR_", PFJetsCHS__genR_, &b_PFJetsCHS__genR_);
  fChain->SetBranchAddress("PFJetsCHS_.cor_", PFJetsCHS__cor_, &b_PFJetsCHS__cor_);
  fChain->SetBranchAddress("PFJetsCHS_.jecLabels_", PFJetsCHS__jecLabels_, &b_PFJetsCHS__jecLabels_);
  fChain->SetBranchAddress("PFJetsCHS_.unc_", PFJetsCHS__unc_, &b_PFJetsCHS__unc_);
  fChain->SetBranchAddress("PFJetsCHS_.uncSrc_", PFJetsCHS__uncSrc_, &b_PFJetsCHS__uncSrc_);
  fChain->SetBranchAddress("PFJetsCHS_.area_", PFJetsCHS__area_, &b_PFJetsCHS__area_);
  fChain->SetBranchAddress("PFJetsCHS_.looseID_", PFJetsCHS__looseID_, &b_PFJetsCHS__looseID_);
  fChain->SetBranchAddress("PFJetsCHS_.tightID_", PFJetsCHS__tightID_, &b_PFJetsCHS__tightID_);
  fChain->SetBranchAddress("PFJetsCHS_.CSVpfPositive_", PFJetsCHS__CSVpfPositive_, &b_PFJetsCHS__CSVpfPositive_);
  fChain->SetBranchAddress("PFJetsCHS_.CSVpfNegative_", PFJetsCHS__CSVpfNegative_, &b_PFJetsCHS__CSVpfNegative_);
  //fChain->SetBranchAddress("PFJetsCHS_.boosted_", PFJetsCHS__boosted_, &b_PFJetsCHS__boosted_);
  fChain->SetBranchAddress("PFJetsCHS_.QGtagger_", PFJetsCHS__QGtagger_, &b_PFJetsCHS__QGtagger_);
  fChain->SetBranchAddress("PFJetsCHS_.partonFlavour_", PFJetsCHS__partonFlavour_, &b_PFJetsCHS__partonFlavour_);
  fChain->SetBranchAddress("PFJetsCHS_.hadronFlavour_", PFJetsCHS__hadronFlavour_, &b_PFJetsCHS__hadronFlavour_);
  fChain->SetBranchAddress("PFJetsCHS_.recommend1_", PFJetsCHS__recommend1_, &b_PFJetsCHS__recommend1_);
  fChain->SetBranchAddress("PFJetsCHS_.recommend2_", PFJetsCHS__recommend2_, &b_PFJetsCHS__recommend2_);
  fChain->SetBranchAddress("PFJetsCHS_.recommend3_", PFJetsCHS__recommend3_, &b_PFJetsCHS__recommend3_);
  //fChain->SetBranchAddress("PFJetsCHS_.pfCombinedCvsL_", PFJetsCHS__pfCombinedCvsL_, &b_PFJetsCHS__pfCombinedCvsL_);
  //fChain->SetBranchAddress("PFJetsCHS_.pfCombinedCvsB_", PFJetsCHS__pfCombinedCvsB_, &b_PFJetsCHS__pfCombinedCvsB_);
  fChain->SetBranchAddress("PFJetsCHS_.chf_", PFJetsCHS__chf_, &b_PFJetsCHS__chf_);
  fChain->SetBranchAddress("PFJetsCHS_.nhf_", PFJetsCHS__nhf_, &b_PFJetsCHS__nhf_);
  fChain->SetBranchAddress("PFJetsCHS_.nemf_", PFJetsCHS__nemf_, &b_PFJetsCHS__nemf_);
  fChain->SetBranchAddress("PFJetsCHS_.cemf_", PFJetsCHS__cemf_, &b_PFJetsCHS__cemf_);
  fChain->SetBranchAddress("PFJetsCHS_.muf_", PFJetsCHS__muf_, &b_PFJetsCHS__muf_);
  fChain->SetBranchAddress("PFJetsCHS_.hf_hf_", PFJetsCHS__hf_hf_, &b_PFJetsCHS__hf_hf_);
  fChain->SetBranchAddress("PFJetsCHS_.hf_phf_", PFJetsCHS__hf_phf_, &b_PFJetsCHS__hf_phf_);
  fChain->SetBranchAddress("PFJetsCHS_.hf_hm_", PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
  fChain->SetBranchAddress("PFJetsCHS_.hf_phm_", PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
  fChain->SetBranchAddress("PFJetsCHS_.chm_", PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
  fChain->SetBranchAddress("PFJetsCHS_.nhm_", PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
  fChain->SetBranchAddress("PFJetsCHS_.phm_", PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
  fChain->SetBranchAddress("PFJetsCHS_.elm_", PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
  fChain->SetBranchAddress("PFJetsCHS_.mum_", PFJetsCHS__mum_, &b_PFJetsCHS__mum_);
  fChain->SetBranchAddress("PFJetsCHS_.ncand_", PFJetsCHS__ncand_, &b_PFJetsCHS__ncand_);
  //fChain->SetBranchAddress("PFJetsCHS_.cm_", PFJetsCHS__cm_, &b_PFJetsCHS__cm_);
  fChain->SetBranchAddress("PFJetsCHS_.beta_", PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
  fChain->SetBranchAddress("PFJetsCHS_.betaStar_", PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
  fChain->SetBranchAddress("PFJetsCHS_.mpuTrk_", PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
  fChain->SetBranchAddress("PFJetsCHS_.mlvTrk_", PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
  fChain->SetBranchAddress("PFJetsCHS_.mjtTrk_", PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
  fChain->SetBranchAddress("PFJetsCHS_.hof_", PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
  fChain->SetBranchAddress("PFJetsCHS_.pujid_", PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
  fChain->SetBranchAddress("PFJetsCHS_.calojetpt_", PFJetsCHS__calojetpt_, &b_PFJetsCHS__calojetpt_);
  fChain->SetBranchAddress("PFJetsCHS_.calojetef_", PFJetsCHS__calojetef_, &b_PFJetsCHS__calojetef_);
  fChain->SetBranchAddress("genFlavour_", &genFlavour_, &b_events_genFlavour_);
  fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_, &b_events_genFlavourHadron_);
}

Bool_t fillHistos::Notify()
{
  return kTRUE;
}

void fillHistos::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t fillHistos::Cut(Long64_t entry)
{
  return 1;
}
#endif // #ifdef fillHistos_cxx
