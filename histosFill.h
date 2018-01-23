// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: April 19, 2010

////////////////////////////////////////////////////////////////////////
// Notes:   Automatically created using TChain::MakeClass("histosFill")
//          Keep variable declarations in the automatic order,
//          update array maximum sizes!!
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  4 11:57:12 2016 by ROOT version 6.06/00
// from TTree ProcessedTree/ProcessedTree
//////////////////////////////////////////////////////////

#ifndef histosFill_h
#define histosFill_h

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
#include "histosBasic.h"
#include "histosEta.h"
#include "histosMC.h"
#include "histosRun.h"
#include "tools.h"
#include "IOV.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

class histosFill {
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
  Int_t           jt3leads[3]; // The three leading jets
  static const Int_t _njt = 250;//100;//90;//kMaxPFJets_;
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
  Float_t         jtgeneta[_njt];   //[njt]
  Float_t         jtgenphi[_njt];   //[njt]

  Float_t         *qgl;
  Float_t         *partonflavor;

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
  Int_t           gen_partonflavor[_njt];

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

  histosFill();
  histosFill(TChain *tree); // custom
  virtual ~histosFill();
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
  map<string, TH1D*> pudist;
  TH1F *pumc;
  //TH1F *pudt;
  set<string> _trigs;
  vector<bool> _jetids;
  map<string, vector<histosBasic*> > _histos;
  map<string, vector<histosEta*> > _etahistos;
  map<string, vector<histosMC*> > _mchistos;
  map<string, histosRun*> _runhistos;
  TH1D *hmcweight;
  TH2D *h2etaphiexcl;
  TFile *fetaphiexcl;
  TH3D *h3probg;

  vector<string> _availTrigs;

  void loadJSON(const char* filename);
  void loadLumi(const char* filename);
  void loadPUProfiles(const char* datafile, const char* mcfile);
  void loadPrescales(const char* prescalefilename);

  void loadECALveto(const char* filename);

  void initBasics(string name);
  void fillBasics(string name);
  void fillBasic(histosBasic *h);
  void writeBasics();

  void initEtas(string name);
  void fillEtas(string name, Float_t *pt, Float_t *eta, Float_t *phi);
  void fillEta(histosEta *h, Float_t *pt, Float_t *eta, Float_t *phi);
  void writeEtas();

  void initMcHistos(string name);
  void fillMcHistos(string name, Float_t *recopt, Float_t *genpt, Float_t *pt, Float_t *eta, Float_t *phi);
  void fillMcHisto(histosMC *h, Float_t *recopt, Float_t *genpt, Float_t *pt, Float_t *eta, Float_t *phi);
  void writeMcHistos();

  void initRunHistos(string name, double etamin, double etamax);
  void fillRunHistos(string name);
  void writeRunHistos();

  void fillJetID(vector<bool> &id);

  bool getTriggers();

  inline double delta_phi(double phi1, double phi2) { // Return value between 0 and phi.
    double dphi = fabs(phi1 - phi2);
    return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
  }

private:
  Long64_t _entries;
  double _xsecMinBias;
  double _w, _w0;
  map<string, double> _wt; // Mainly for trigger pileup weigths

  vector<int> _jkmore;
  map<Int_t, Int_t> _outlist;

  TLorentzVector p4, gp4, genp4, cp4, pp4;
  jec::IOV _iov;
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
  //map<string, map<int, int> > _premap;
  map<string, map<int, map<int, int> > > _premap;
  map<string, map<int, int> > _prescales;
  TH2F *_ecalveto;

  map<int, map<int, set<int> > > _duplicates;
  set<int> _badruns;
  set<pair<int, int> > _badlums;
  set<pair<int, int> > _nolums;
  set<pair<int, int> > _badjson;
  set<pair<int, int> > _passedlumis;
  int _nbadevts_dup;
  int _nbadevts_run;
  int _nbadevts_ls;
  int _nbadevts_lum;
  int _nbadevts_json;
  int _nbadevts_veto;
  int _bscounter_bad;
  int _bscounter_good;
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

#ifdef histosFill_cxx
histosFill::histosFill()
{
  cout << "This is not supposed to happen! Provide a tree" << endl;
}

histosFill::histosFill(TChain *tree) : fChain(0)
{
  assert(tree);
  Init(tree); // custom

  // Reset some pointers
  _outfile = 0;
}

histosFill::~histosFill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t histosFill::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// Check that the correct tree is open in the chain
Long64_t histosFill::LoadTree(Long64_t entry)
{
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;

  // A new tree is opened
  if (fChain->GetTreeNumber() != fCurrent) {
    // Reload the triggers and print them
    fCurrent = fChain->GetTreeNumber();
    if (_jp_isdt) {
      *ferr << "Opening tree number " << fChain->GetTreeNumber() << endl;
      if (!getTriggers()) {
        *ferr << "Failed to load DT triggers. Check that the SMPJ tuple has the required histograms. Aborting..." << endl;
        return -4;
      }
      *ferr << "Tree " << fCurrent << " triggers:" << endl;
      for (auto str : _availTrigs) {
        if (str.length()==0)
          *ferr << "x ";
        else
          *ferr << str << " ";
      }
      *ferr << endl << flush;
    } else {
      TString filename = fChain->GetCurrentFile()->GetName();
      // Check the position of the current file in the list of file names
      unsigned currFile = std::find_if(_jp_pthatfiles.begin(),_jp_pthatfiles.end(),[&filename] (string s) { return filename.Contains(s); })-_jp_pthatfiles.begin();
      if (_jp_pthatnevts[currFile]<=0.0 or _jp_pthatsigmas[currFile]<=0) {
        *ferr << "Suspicious pthat slice information for file " << _jp_pthatfiles[currFile] << ". Aborting..." << endl;
        return -3;
      }
      _pthatweight = _jp_pthatsigmas[currFile]/_jp_pthatnevts[currFile];
      _pthatweight /= _jp_pthatsigmas[_jp_npthatbins-1]/_jp_pthatnevts[_jp_npthatbins-1]; // Normalize
      *ferr << "Pthat bin changing." << endl << "File " << currFile << " " << fChain->GetCurrentFile()->GetName();
      *ferr << " should correspond to the range [" << _jp_pthatranges[currFile] << "," << _jp_pthatranges[currFile+1] << "]";
      *ferr << endl << "Weight: " << _pthatweight << endl;
    }
    // slices with pthat bins
  }
  return centry;
}

void histosFill::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

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
  fChain->SetBranchAddress(Form("PFJets%s_",_jp_chs), &PFJetsCHS__, &b_events_PFJetsCHS__);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fX",_jp_chs), PFJetsCHS__P4__fCoordinates_fX, &b_PFJetsCHS__P4__fCoordinates_fX);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fY",_jp_chs), PFJetsCHS__P4__fCoordinates_fY, &b_PFJetsCHS__P4__fCoordinates_fY);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fZ",_jp_chs), PFJetsCHS__P4__fCoordinates_fZ, &b_PFJetsCHS__P4__fCoordinates_fZ);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fT",_jp_chs), PFJetsCHS__P4__fCoordinates_fT, &b_PFJetsCHS__P4__fCoordinates_fT);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fX",_jp_chs), PFJetsCHS__genP4__fCoordinates_fX, &b_PFJetsCHS__genP4__fCoordinates_fX);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fY",_jp_chs), PFJetsCHS__genP4__fCoordinates_fY, &b_PFJetsCHS__genP4__fCoordinates_fY);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fZ",_jp_chs), PFJetsCHS__genP4__fCoordinates_fZ, &b_PFJetsCHS__genP4__fCoordinates_fZ);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fT",_jp_chs), PFJetsCHS__genP4__fCoordinates_fT, &b_PFJetsCHS__genP4__fCoordinates_fT);
  fChain->SetBranchAddress(Form("PFJets%s_.genR_",_jp_chs), PFJetsCHS__genR_, &b_PFJetsCHS__genR_);
  fChain->SetBranchAddress(Form("PFJets%s_.cor_",_jp_chs), PFJetsCHS__cor_, &b_PFJetsCHS__cor_);
  fChain->SetBranchAddress(Form("PFJets%s_.jecLabels_",_jp_chs), PFJetsCHS__jecLabels_, &b_PFJetsCHS__jecLabels_);
  fChain->SetBranchAddress(Form("PFJets%s_.unc_",_jp_chs), PFJetsCHS__unc_, &b_PFJetsCHS__unc_);
  fChain->SetBranchAddress(Form("PFJets%s_.uncSrc_",_jp_chs), PFJetsCHS__uncSrc_, &b_PFJetsCHS__uncSrc_);
  fChain->SetBranchAddress(Form("PFJets%s_.area_",_jp_chs), PFJetsCHS__area_, &b_PFJetsCHS__area_);
  fChain->SetBranchAddress(Form("PFJets%s_.looseID_",_jp_chs), PFJetsCHS__looseID_, &b_PFJetsCHS__looseID_);
  fChain->SetBranchAddress(Form("PFJets%s_.tightID_",_jp_chs), PFJetsCHS__tightID_, &b_PFJetsCHS__tightID_);
  fChain->SetBranchAddress(Form("PFJets%s_.CSVpfPositive_",_jp_chs), PFJetsCHS__CSVpfPositive_, &b_PFJetsCHS__CSVpfPositive_);
  fChain->SetBranchAddress(Form("PFJets%s_.CSVpfNegative_",_jp_chs), PFJetsCHS__CSVpfNegative_, &b_PFJetsCHS__CSVpfNegative_);
  //fChain->SetBranchAddress(Form("PFJets%s_.boosted_",_jp_chs), PFJetsCHS__boosted_, &b_PFJetsCHS__boosted_);
  fChain->SetBranchAddress(Form("PFJets%s_.QGtagger_",_jp_chs), PFJetsCHS__QGtagger_, &b_PFJetsCHS__QGtagger_);
  fChain->SetBranchAddress(Form("PFJets%s_.partonFlavour_",_jp_chs), PFJetsCHS__partonFlavour_, &b_PFJetsCHS__partonFlavour_);
  fChain->SetBranchAddress(Form("PFJets%s_.hadronFlavour_",_jp_chs), PFJetsCHS__hadronFlavour_, &b_PFJetsCHS__hadronFlavour_);
  fChain->SetBranchAddress(Form("PFJets%s_.recommend1_",_jp_chs), PFJetsCHS__recommend1_, &b_PFJetsCHS__recommend1_);
  fChain->SetBranchAddress(Form("PFJets%s_.recommend2_",_jp_chs), PFJetsCHS__recommend2_, &b_PFJetsCHS__recommend2_);
  fChain->SetBranchAddress(Form("PFJets%s_.recommend3_",_jp_chs), PFJetsCHS__recommend3_, &b_PFJetsCHS__recommend3_);
  //fChain->SetBranchAddress(Form("PFJets%s_.pfCombinedCvsL_",_jp_chs), PFJetsCHS__pfCombinedCvsL_, &b_PFJetsCHS__pfCombinedCvsL_);
  //fChain->SetBranchAddress(Form("PFJets%s_.pfCombinedCvsB_",_jp_chs), PFJetsCHS__pfCombinedCvsB_, &b_PFJetsCHS__pfCombinedCvsB_);
  fChain->SetBranchAddress(Form("PFJets%s_.chf_",_jp_chs), PFJetsCHS__chf_, &b_PFJetsCHS__chf_);
  fChain->SetBranchAddress(Form("PFJets%s_.nhf_",_jp_chs), PFJetsCHS__nhf_, &b_PFJetsCHS__nhf_);
  fChain->SetBranchAddress(Form("PFJets%s_.nemf_",_jp_chs), PFJetsCHS__nemf_, &b_PFJetsCHS__nemf_);
  fChain->SetBranchAddress(Form("PFJets%s_.cemf_",_jp_chs), PFJetsCHS__cemf_, &b_PFJetsCHS__cemf_);
  fChain->SetBranchAddress(Form("PFJets%s_.muf_",_jp_chs), PFJetsCHS__muf_, &b_PFJetsCHS__muf_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_hf_",_jp_chs), PFJetsCHS__hf_hf_, &b_PFJetsCHS__hf_hf_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_phf_",_jp_chs), PFJetsCHS__hf_phf_, &b_PFJetsCHS__hf_phf_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_hm_",_jp_chs), PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_phm_",_jp_chs), PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.chm_",_jp_chs), PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
  fChain->SetBranchAddress(Form("PFJets%s_.nhm_",_jp_chs), PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
  fChain->SetBranchAddress(Form("PFJets%s_.phm_",_jp_chs), PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.elm_",_jp_chs), PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
  fChain->SetBranchAddress(Form("PFJets%s_.mum_",_jp_chs), PFJetsCHS__mum_, &b_PFJetsCHS__mum_);
  fChain->SetBranchAddress(Form("PFJets%s_.ncand_",_jp_chs), PFJetsCHS__ncand_, &b_PFJetsCHS__ncand_);
  //fChain->SetBranchAddress("PFJets%s_.cm_",_jp_chs), PFJetsCHS__cm_, &b_PFJetsCHS__cm_);
  fChain->SetBranchAddress(Form("PFJets%s_.beta_",_jp_chs), PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
  fChain->SetBranchAddress(Form("PFJets%s_.betaStar_",_jp_chs), PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
  fChain->SetBranchAddress(Form("PFJets%s_.mpuTrk_",_jp_chs), PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.mlvTrk_",_jp_chs), PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.mjtTrk_",_jp_chs), PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.hof_",_jp_chs), PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
  fChain->SetBranchAddress(Form("PFJets%s_.pujid_",_jp_chs), PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
  fChain->SetBranchAddress(Form("PFJets%s_.calojetpt_",_jp_chs), PFJetsCHS__calojetpt_, &b_PFJetsCHS__calojetpt_);
  fChain->SetBranchAddress(Form("PFJets%s_.calojetef_",_jp_chs), PFJetsCHS__calojetef_, &b_PFJetsCHS__calojetef_);
  fChain->SetBranchAddress("genFlavour_", &genFlavour_, &b_events_genFlavour_);
  fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_, &b_events_genFlavourHadron_);
}

Bool_t histosFill::Notify()
{
  return kTRUE;
}

void histosFill::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t histosFill::Cut(Long64_t entry)
{
  return 1;
}
#endif // #ifdef histosFill_cxx
