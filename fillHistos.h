// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: Nov 28, 2015
////////////////////////////////////////////////////////////////////////
// Notes:   Automatically created using TChain::MakeClass("fillHistos")
//          Keep variable declarations in the automatic order,
//          update array maximum sizes!!
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 26 11:41:50 2015 by ROOT version 6.02/10
// from TTree ProcessedTree/ProcessedTree
// found on file: data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root
//////////////////////////////////////////////////////////

#ifndef fillHistos_h
#define fillHistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <map>
#include <vector>
#include <string>
#include <set>

#include "settings.h"
#include "basicHistos.h"
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
   TTree          *fChain2;  // _ak4ak8
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxfilterIdList = 1;
   const Int_t kMaxCaloMet__et = 1;
   const Int_t kMaxCaloMet__CaloMetPt = 1;
   const Int_t kMaxCaloMet__sumEt = 1;
   const Int_t kMaxCaloMet__phi = 1;
   const Int_t kMaxPFMet__et = 1;
   const Int_t kMaxPFMet__CaloMetPt = 1;
   const Int_t kMaxPFMet__sumEt = 1;
   const Int_t kMaxPFMet__phi = 1;
   const Int_t kMaxTriggerDecision = 1;
   const Int_t kMaxtriggerList = 1;
   const Int_t kMaxL1Prescale = 1;
   const Int_t kMaxHLTPrescale = 1;
   const Int_t kMaxHLTObj = 1;
   const Int_t kMaxL1Obj = 1;
   static const Int_t kMaxGenJets_ = 100;
   static const Int_t kMaxCaloJets_ = 100;
/*    const Int_t kMaxCaloJets__genR = 1; */
/*    const Int_t kMaxCaloJets__cor = 1; */
/*    const Int_t kMaxCaloJets__jecLabels = 1; */
/*    const Int_t kMaxCaloJets__unc = 1; */
/*    const Int_t kMaxCaloJets__uncSrc = 1; */
/*    const Int_t kMaxCaloJets__area = 1; */
/*    const Int_t kMaxCaloJets__looseID = 1; */
/*    const Int_t kMaxCaloJets__tightID = 1; */
/*    const Int_t kMaxCaloJets__emf = 1; */
/*    const Int_t kMaxCaloJets__fHPD = 1; */
/*    const Int_t kMaxCaloJets__fRBX = 1; */
/*    const Int_t kMaxCaloJets__n90hits = 1; */
/*    const Int_t kMaxCaloJets__nTrkCalo = 1; */
/*    const Int_t kMaxCaloJets__nTrkVtx = 1; */
   static const Int_t kMaxPFJets_ = 100;
/*    const Int_t kMaxPFJets__genR = 100; */
/*    const Int_t kMaxPFJets__cor = 100; */
/*    const Int_t kMaxPFJets__jecLabels = 100; */
/*    const Int_t kMaxPFJets__unc = 100; */
/*    const Int_t kMaxPFJets__uncSrc = 100; */
/*    const Int_t kMaxPFJets__area = 100; */
/*    const Int_t kMaxPFJets__looseID = 100; */
/*    const Int_t kMaxPFJets__tightID = 100; */
/*    const Int_t kMaxPFJets__chf = 100; */
/*    const Int_t kMaxPFJets__nhf = 100; */
/*    const Int_t kMaxPFJets__nemf = 100; */
/*    const Int_t kMaxPFJets__cemf = 100; */
/*    const Int_t kMaxPFJets__muf = 100; */
/*    const Int_t kMaxPFJets__hf_hf = 100; */
/*    const Int_t kMaxPFJets__hf_phf = 100; */
/*    const Int_t kMaxPFJets__hf_hm = 100; */
/*    const Int_t kMaxPFJets__hf_phm = 100; */
/*    const Int_t kMaxPFJets__chm = 100; */
/*    const Int_t kMaxPFJets__nhm = 100; */
/*    const Int_t kMaxPFJets__phm = 100; */
/*    const Int_t kMaxPFJets__elm = 100; */
/*    const Int_t kMaxPFJets__mum = 100; */
/*    const Int_t kMaxPFJets__ncand = 100; */
/*    const Int_t kMaxPFJets__beta = 100; */
/*    const Int_t kMaxPFJets__betaStar = 100; */
/*    const Int_t kMaxPFJets__mpuTrk = 100; */
/*    const Int_t kMaxPFJets__mlvTrk = 100; */
/*    const Int_t kMaxPFJets__mjtTrk = 100; */
/*    const Int_t kMaxPFJets__hof = 100; */
/*    const Int_t kMaxPFJets__pujid = 100; */
/*    const Int_t kMaxPFJets__calojetpt = 100; */
/*    const Int_t kMaxPFJets__calojetef = 100; */
   static const Int_t kMaxPFJetsCHS_ = 100;
/*    const Int_t kMaxPFJetsCHS__genR = 100; */
/*    const Int_t kMaxPFJetsCHS__cor = 100; */
/*    const Int_t kMaxPFJetsCHS__jecLabels = 100; */
/*    const Int_t kMaxPFJetsCHS__unc = 100; */
/*    const Int_t kMaxPFJetsCHS__uncSrc = 100; */
/*    const Int_t kMaxPFJetsCHS__area = 100; */
/*    const Int_t kMaxPFJetsCHS__looseID = 100; */
/*    const Int_t kMaxPFJetsCHS__tightID = 100; */
/*    const Int_t kMaxPFJetsCHS__chf = 100; */
/*    const Int_t kMaxPFJetsCHS__nhf = 100; */
/*    const Int_t kMaxPFJetsCHS__nemf = 100; */
/*    const Int_t kMaxPFJetsCHS__cemf = 100; */
/*    const Int_t kMaxPFJetsCHS__muf = 100; */
/*    const Int_t kMaxPFJetsCHS__hf_hf = 100; */
/*    const Int_t kMaxPFJetsCHS__hf_phf = 100; */
/*    const Int_t kMaxPFJetsCHS__hf_hm = 100; */
/*    const Int_t kMaxPFJetsCHS__hf_phm = 100; */
/*    const Int_t kMaxPFJetsCHS__chm = 100; */
/*    const Int_t kMaxPFJetsCHS__nhm = 100; */
/*    const Int_t kMaxPFJetsCHS__phm = 100; */
/*    const Int_t kMaxPFJetsCHS__elm = 100; */
/*    const Int_t kMaxPFJetsCHS__mum = 100; */
/*    const Int_t kMaxPFJetsCHS__ncand = 100; */
/*    const Int_t kMaxPFJetsCHS__beta = 100; */
/*    const Int_t kMaxPFJetsCHS__betaStar = 100; */
/*    const Int_t kMaxPFJetsCHS__mpuTrk = 100; */
/*    const Int_t kMaxPFJetsCHS__mlvTrk = 100; */
/*    const Int_t kMaxPFJetsCHS__mjtTrk = 100; */
/*    const Int_t kMaxPFJetsCHS__hof = 100; */
/*    const Int_t kMaxPFJetsCHS__pujid = 100; */
/*    const Int_t kMaxPFJetsCHS__calojetpt = 100; */
/*    const Int_t kMaxPFJetsCHS__calojetef = 100; */

   // Declaration of leaf types
 //QCDEvent        *events;
   vector<vector<int> > filterIdList_;
   Bool_t          EvtHdr__mIsPVgood;
   Bool_t          EvtHdr__mHCALNoise;
   Bool_t          EvtHdr__mHCALNoiseNoMinZ;
   Int_t           EvtHdr__mRun;
   Int_t           EvtHdr__mEvent;
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
   Int_t           CaloJets__;
   Double_t        CaloJets__P4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__genR_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__cor_[kMaxCaloJets_];   //[CaloJets__]
   vector<double>  CaloJets__jecLabels_[kMaxCaloJets_];
   Float_t         CaloJets__unc_[kMaxCaloJets_];   //[CaloJets__]
   vector<float>   CaloJets__uncSrc_[kMaxCaloJets_];
   Float_t         CaloJets__area_[kMaxCaloJets_];   //[CaloJets__]
   Bool_t          CaloJets__looseID_[kMaxCaloJets_];   //[CaloJets__]
   Bool_t          CaloJets__tightID_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__emf_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__fHPD_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__fRBX_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__n90hits_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__nTrkCalo_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__nTrkVtx_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           PFJets__;
   Double_t        PFJets__P4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__genR_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__cor_[kMaxPFJets_];   //[PFJets__]
   vector<double>  PFJets__jecLabels_[kMaxPFJets_];
   Float_t         PFJets__unc_[kMaxPFJets_];   //[PFJets__]
   vector<float>   PFJets__uncSrc_[kMaxPFJets_];
   Float_t         PFJets__area_[kMaxPFJets_];   //[PFJets__]
   Bool_t          PFJets__looseID_[kMaxPFJets_];   //[PFJets__]
   Bool_t          PFJets__tightID_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__chf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__nhf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__nemf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__cemf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__muf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hf_hf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hf_phf_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__hf_hm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__hf_phm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__chm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__nhm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__phm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__elm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mum_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__ncand_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__beta_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__betaStar_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mpuTrk_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mlvTrk_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mjtTrk_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hof_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__pujid_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__calojetpt_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__calojetef_[kMaxPFJets_];   //[PFJets__]
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
   Float_t         PFJetsCHS__beta_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__betaStar_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mpuTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mlvTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mjtTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hof_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pujid_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__calojetpt_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__calojetef_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

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
   TBranch        *b_events_TriggerDecision_;   //!
   TBranch        *b_events_triggerList_;   //!
   TBranch        *b_events_L1Prescale_;   //!
   TBranch        *b_events_HLTPrescale_;   //!
   TBranch        *b_events_GenJets__;   //!
   TBranch        *b_GenJets__fCoordinates_fX;   //!
   TBranch        *b_GenJets__fCoordinates_fY;   //!
   TBranch        *b_GenJets__fCoordinates_fZ;   //!
   TBranch        *b_GenJets__fCoordinates_fT;   //!
   TBranch        *b_events_CaloJets__;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fX;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fY;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fZ;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fT;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fX;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fY;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fZ;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fT;   //!
   TBranch        *b_CaloJets__genR_;   //!
   TBranch        *b_CaloJets__cor_;   //!
   TBranch        *b_CaloJets__jecLabels_;   //!
   TBranch        *b_CaloJets__unc_;   //!
   TBranch        *b_CaloJets__uncSrc_;   //!
   TBranch        *b_CaloJets__area_;   //!
   TBranch        *b_CaloJets__looseID_;   //!
   TBranch        *b_CaloJets__tightID_;   //!
   TBranch        *b_CaloJets__emf_;   //!
   TBranch        *b_CaloJets__fHPD_;   //!
   TBranch        *b_CaloJets__fRBX_;   //!
   TBranch        *b_CaloJets__n90hits_;   //!
   TBranch        *b_CaloJets__nTrkCalo_;   //!
   TBranch        *b_CaloJets__nTrkVtx_;   //!
   TBranch        *b_events_PFJets__;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fX;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fY;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fZ;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fT;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fX;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fY;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fZ;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fT;   //!
   TBranch        *b_PFJets__genR_;   //!
   TBranch        *b_PFJets__cor_;   //!
   TBranch        *b_PFJets__jecLabels_;   //!
   TBranch        *b_PFJets__unc_;   //!
   TBranch        *b_PFJets__uncSrc_;   //!
   TBranch        *b_PFJets__area_;   //!
   TBranch        *b_PFJets__looseID_;   //!
   TBranch        *b_PFJets__tightID_;   //!
   TBranch        *b_PFJets__chf_;   //!
   TBranch        *b_PFJets__nhf_;   //!
   TBranch        *b_PFJets__nemf_;   //!
   TBranch        *b_PFJets__cemf_;   //!
   TBranch        *b_PFJets__muf_;   //!
   TBranch        *b_PFJets__hf_hf_;   //!
   TBranch        *b_PFJets__hf_phf_;   //!
   TBranch        *b_PFJets__hf_hm_;   //!
   TBranch        *b_PFJets__hf_phm_;   //!
   TBranch        *b_PFJets__chm_;   //!
   TBranch        *b_PFJets__nhm_;   //!
   TBranch        *b_PFJets__phm_;   //!
   TBranch        *b_PFJets__elm_;   //!
   TBranch        *b_PFJets__mum_;   //!
   TBranch        *b_PFJets__ncand_;   //!
   TBranch        *b_PFJets__beta_;   //!
   TBranch        *b_PFJets__betaStar_;   //!
   TBranch        *b_PFJets__mpuTrk_;   //!
   TBranch        *b_PFJets__mlvTrk_;   //!
   TBranch        *b_PFJets__mjtTrk_;   //!
   TBranch        *b_PFJets__hof_;   //!
   TBranch        *b_PFJets__pujid_;   //!
   TBranch        *b_PFJets__calojetpt_;   //!
   TBranch        *b_PFJets__calojetef_;   //!
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
   TBranch        *b_PFJetsCHS__beta_;   //!
   TBranch        *b_PFJetsCHS__betaStar_;   //!
   TBranch        *b_PFJetsCHS__mpuTrk_;   //!
   TBranch        *b_PFJetsCHS__mlvTrk_;   //!
   TBranch        *b_PFJetsCHS__mjtTrk_;   //!
   TBranch        *b_PFJetsCHS__hof_;   //!
   TBranch        *b_PFJetsCHS__pujid_;   //!
   TBranch        *b_PFJetsCHS__calojetpt_;   //!
   TBranch        *b_PFJetsCHS__calojetef_;   //!

  // _ak4ak8 tree
   Int_t           t4_EvtHdr__mRun;
   UInt_t          t4_EvtHdr__mEvent;
   Int_t           t4_EvtHdr__mLumi;
   Float_t         t4_EvtHdr__mPthat;
   Int_t           t4_EvtHdr__mOOTPUEarly;
   Int_t           t4_EvtHdr__mOOTPULate;
   Int_t           t4_EvtHdr__mINTPU;
   //
   Int_t           t4_PFJets__;
   Double_t        t4_PFJets__P4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
   Double_t        t4_PFJets__P4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
   Double_t        t4_PFJets__P4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
   Double_t        t4_PFJets__P4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
   Float_t         t4_PFJets__cor_[kMaxPFJets_];   //[PFJets__]
   Float_t         t4_PFJets__area_[kMaxPFJets_];   //[PFJets__]
   Bool_t          t4_PFJets__looseID_[kMaxPFJets_];   //[PFJets__]
   Bool_t          t4_PFJets__tightID_[kMaxPFJets_];   //[PFJets__]
   //
   Int_t           t4_GenJets__;
   Double_t        t4_GenJets__fCoordinates_fX[kMaxGenJets_];
   Double_t        t4_GenJets__fCoordinates_fY[kMaxGenJets_];
   Double_t        t4_GenJets__fCoordinates_fZ[kMaxGenJets_];
   Double_t        t4_GenJets__fCoordinates_fT[kMaxGenJets_];

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
   //
   static const bool dofriends = false;//true;

   // _ak4ak8
   Int_t           t4_njt;
   Double_t        *t4_jtp4x;//[_njt];   //[njt]
   Double_t        *t4_jtp4y;//[_njt];   //[njt]
   Double_t        *t4_jtp4z;//[_njt];   //[njt]
   Double_t        *t4_jtp4t;//[_njt];   //[njt]
   Float_t         t4_jte[_njt];   //[njt]
   Float_t         t4_jtpt[_njt];   //[njt]
   Float_t         t4_jtptu[_njt];   //EXTRA
   Float_t         t4_jteu[_njt];   //EXTRA
   Float_t         t4_jteta[_njt];   //[njt]
   Float_t         t4_jtphi[_njt];   //[njt]
   Float_t         t4_jty[_njt];   //[njt]
   Float_t         *t4_jta;//[_njt];   //[njt]
   Float_t         *t4_jtjes;//[_njt];   //[njt]
   Float_t         t4_jtjesnew[_njt];   //[njt]
   Bool_t          *t4_jtidloose;//[_njt];   //[njt]
   Bool_t          *t4_jtidtight;//[_njt];   //[njt]
   //
   Int_t           t4gen_njt;
   Double_t        *t4gen_jtp4x;
   Double_t        *t4gen_jtp4y;
   Double_t        *t4gen_jtp4z;
   Double_t        *t4gen_jtp4t;
   Float_t         t4gen_jtpt[_njt];   //[njt]
   Float_t         t4gen_jteta[_njt];   //[njt]
   Float_t         t4gen_jtphi[_njt];   //[njt]
   Float_t         t4gen_jty[_njt];   //[njt]
   // For trigger matching
   Int_t           t4c_njt;
   Double_t        *t4c_jtp4x;//[_njt];   //[njt]
   Double_t        *t4c_jtp4y;//[_njt];   //[njt]
   Double_t        *t4c_jtp4z;//[_njt];   //[njt]
   Double_t        *t4c_jtp4t;//[_njt];   //[njt]
   Float_t         t4c_jte[_njt];   //[njt]
   Float_t         t4c_jtpt[_njt];   //[njt]
   Float_t         t4c_jtptu[_njt];   //EXTRA
   Float_t         t4c_jteu[_njt];   //EXTRA
   Float_t         t4c_jteta[_njt];   //[njt]
   Float_t         t4c_jtphi[_njt];   //[njt]
   Float_t         t4c_jty[_njt];   //[njt]
   //Float_t         *t4c_jta;//[_njt];   //[njt]
   Float_t         t4c_jta;//[_njt];   //[njt]
   Float_t         *t4c_jtjes;//[_njt];   //[njt]
   Float_t         t4c_jtjesnew[_njt];   //[njt]

   //
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

   // Composition corrected for PU
   Int_t           jtn2[_njt];   //[njt]
   Int_t           jtnch2[_njt];   //[njt]
   Int_t           jtnnh2[_njt];   //[njt]
   Int_t           jtnne2[_njt];   //[njt]
   Int_t           jtnce2[_njt];   //[njt]
   Int_t           jtnmu2[_njt];   //[njt]
   Float_t         jtchf2[_njt];   //[njt]
   Float_t         jtnhf2[_njt];   //[njt]
   Float_t         jtnef2[_njt];   //[njt]
   Float_t         jtcef2[_njt];   //[njt]
   Float_t         jtmuf2[_njt];   //[njt]

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

   fillHistos(TTree *tree=0, std::string type="DATA", TTree *tree2=0);
   virtual ~fillHistos();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tree2=0); // custom
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1, Long64_t entry4 = -1); // custom

   // Additional code for processing
   ofstream *ferr;
   std::string _type;
   TFile *_outfile;
   //TH1D *pudist;
   std::map<string, TH1D*> pudist;
   std::vector<std::string> _triggers;
   TH1D *pumc;
   std::set<std::string> _trigs;
   std::vector<bool> _jetids;
   std::vector<bool> _jetids5;
   std::map<std::string, std::vector<basicHistos*> > _histos;
   std::map<std::string, runHistos*> _runhistos;
   TH1D *hmcweight;

   void loadJSON(const char* filename);
   void loadLumi(const char* filename);
   //void loadPrescale(const char* filename, std::string trigger);
   void loadPUProfiles(const char* datafile, const char* mcfile);
   void loadECALveto(const char* filename);

   void initBasics(std::string name);
   void fillBasics(std::string name);
   void fillBasic(basicHistos *h);
   void writeBasics();

   void initRunHistos(std::string name, double ymin, double ymax);
   void fillRunHistos(std::string name);
   void writeRunHistos();

   void fillJetID(std::vector<bool> &id);
   void fillJetID5(std::vector<bool> &id);

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
   std::map<string, double> _wt;

   std::vector<int> _jkmore;
   std::map<Int_t, Int_t> _outlist;

   TLorentzVector p4, gp4, genp4, cp4, pp4;
   jec::IOV *iov;
   FactorizedJetCorrector *_JEC, *_L1RC;
   FactorizedJetCorrector *_JEC_ak4pf;
   JetCorrectionUncertainty *_jecUnc;
   //L3Corr *_jecUnc2;

   bool   _pass;

   std::set<int> _runveto;
   std::map<int, std::map<int, int> > _json;
   std::map<int, std::map<int, float> > _lums;
   std::map<int, std::map<int, float> > _lums2;
   std::map<std::string, std::map<int, int> > _prescales;
   TH2F *ecalveto;

   std::map<int, std::map<int, std::set<int> > > _duplicates;
   std::set<int> _badruns;
   std::set<std::pair<int, int> > _badlums;
   std::set<std::pair<int, int> > _nolums;
   std::set<std::pair<int, int> > _badjson;
   std::set<std::pair<int, int> > _jt15lums;
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

   TLorentzVector j1, j2;
};

#endif

#ifdef fillHistos_cxx
fillHistos::fillHistos(TTree *tree, std::string type, TTree *tree2)
{
  _type = type;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root:/ak4");
      dir->GetObject("ProcessedTree",tree);

   }
*/
  assert(tree);
  Init(tree, tree2); // custom

  // Reset some pointers
  _outfile = 0;

  // Awkward patch for ROOT6 to call Loop
  Loop();
  //delete this; // Needed?
}

fillHistos::~fillHistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fillHistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fillHistos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fillHistos::Init(TTree *tree, TTree *tree2)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain2 = tree2;
   assert(!_jp_ak4ak8 || fChain2);
   if (fChain2) fChain2->SetMakeClass(1);

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
   fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
   fChain->SetBranchAddress("triggerList_", &triggerList_, &b_events_triggerList_);
   fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
   fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
   fChain->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fX", &GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fY", &GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", &GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fT", &GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_", &CaloJets__, &b_events_CaloJets__);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fX", &CaloJets__P4__fCoordinates_fX, &b_CaloJets__P4__fCoordinates_fX);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fY", &CaloJets__P4__fCoordinates_fY, &b_CaloJets__P4__fCoordinates_fY);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fZ", &CaloJets__P4__fCoordinates_fZ, &b_CaloJets__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fT", &CaloJets__P4__fCoordinates_fT, &b_CaloJets__P4__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fX", &CaloJets__genP4__fCoordinates_fX, &b_CaloJets__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fY", &CaloJets__genP4__fCoordinates_fY, &b_CaloJets__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fZ", &CaloJets__genP4__fCoordinates_fZ, &b_CaloJets__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fT", &CaloJets__genP4__fCoordinates_fT, &b_CaloJets__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_.genR_", &CaloJets__genR_, &b_CaloJets__genR_);
   fChain->SetBranchAddress("CaloJets_.cor_", &CaloJets__cor_, &b_CaloJets__cor_);
   fChain->SetBranchAddress("CaloJets_.jecLabels_", &CaloJets__jecLabels_, &b_CaloJets__jecLabels_);
   fChain->SetBranchAddress("CaloJets_.unc_", &CaloJets__unc_, &b_CaloJets__unc_);
   fChain->SetBranchAddress("CaloJets_.uncSrc_", &CaloJets__uncSrc_, &b_CaloJets__uncSrc_);
   fChain->SetBranchAddress("CaloJets_.area_", &CaloJets__area_, &b_CaloJets__area_);
   fChain->SetBranchAddress("CaloJets_.looseID_", &CaloJets__looseID_, &b_CaloJets__looseID_);
   fChain->SetBranchAddress("CaloJets_.tightID_", &CaloJets__tightID_, &b_CaloJets__tightID_);
   fChain->SetBranchAddress("CaloJets_.emf_", &CaloJets__emf_, &b_CaloJets__emf_);
   fChain->SetBranchAddress("CaloJets_.fHPD_", &CaloJets__fHPD_, &b_CaloJets__fHPD_);
   fChain->SetBranchAddress("CaloJets_.fRBX_", &CaloJets__fRBX_, &b_CaloJets__fRBX_);
   fChain->SetBranchAddress("CaloJets_.n90hits_", &CaloJets__n90hits_, &b_CaloJets__n90hits_);
   fChain->SetBranchAddress("CaloJets_.nTrkCalo_", &CaloJets__nTrkCalo_, &b_CaloJets__nTrkCalo_);
   fChain->SetBranchAddress("CaloJets_.nTrkVtx_", &CaloJets__nTrkVtx_, &b_CaloJets__nTrkVtx_);
   fChain->SetBranchAddress("PFJets_", &PFJets__, &b_events_PFJets__);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fX", PFJets__P4__fCoordinates_fX, &b_PFJets__P4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fY", PFJets__P4__fCoordinates_fY, &b_PFJets__P4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fZ", PFJets__P4__fCoordinates_fZ, &b_PFJets__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fT", PFJets__P4__fCoordinates_fT, &b_PFJets__P4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fX", PFJets__genP4__fCoordinates_fX, &b_PFJets__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fY", PFJets__genP4__fCoordinates_fY, &b_PFJets__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fZ", PFJets__genP4__fCoordinates_fZ, &b_PFJets__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fT", PFJets__genP4__fCoordinates_fT, &b_PFJets__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJets_.genR_", PFJets__genR_, &b_PFJets__genR_);
   fChain->SetBranchAddress("PFJets_.cor_", PFJets__cor_, &b_PFJets__cor_);
   fChain->SetBranchAddress("PFJets_.jecLabels_", PFJets__jecLabels_, &b_PFJets__jecLabels_);
   fChain->SetBranchAddress("PFJets_.unc_", PFJets__unc_, &b_PFJets__unc_);
   fChain->SetBranchAddress("PFJets_.uncSrc_", PFJets__uncSrc_, &b_PFJets__uncSrc_);
   fChain->SetBranchAddress("PFJets_.area_", PFJets__area_, &b_PFJets__area_);
   fChain->SetBranchAddress("PFJets_.looseID_", PFJets__looseID_, &b_PFJets__looseID_);
   fChain->SetBranchAddress("PFJets_.tightID_", PFJets__tightID_, &b_PFJets__tightID_);
   fChain->SetBranchAddress("PFJets_.chf_", PFJets__chf_, &b_PFJets__chf_);
   fChain->SetBranchAddress("PFJets_.nhf_", PFJets__nhf_, &b_PFJets__nhf_);
   fChain->SetBranchAddress("PFJets_.nemf_", PFJets__nemf_, &b_PFJets__nemf_);
   fChain->SetBranchAddress("PFJets_.cemf_", PFJets__cemf_, &b_PFJets__cemf_);
   fChain->SetBranchAddress("PFJets_.muf_", PFJets__muf_, &b_PFJets__muf_);
   fChain->SetBranchAddress("PFJets_.hf_hf_", PFJets__hf_hf_, &b_PFJets__hf_hf_);
   fChain->SetBranchAddress("PFJets_.hf_phf_", PFJets__hf_phf_, &b_PFJets__hf_phf_);
   fChain->SetBranchAddress("PFJets_.hf_hm_", PFJets__hf_hm_, &b_PFJets__hf_hm_);
   fChain->SetBranchAddress("PFJets_.hf_phm_", PFJets__hf_phm_, &b_PFJets__hf_phm_);
   fChain->SetBranchAddress("PFJets_.chm_", PFJets__chm_, &b_PFJets__chm_);
   fChain->SetBranchAddress("PFJets_.nhm_", PFJets__nhm_, &b_PFJets__nhm_);
   fChain->SetBranchAddress("PFJets_.phm_", PFJets__phm_, &b_PFJets__phm_);
   fChain->SetBranchAddress("PFJets_.elm_", PFJets__elm_, &b_PFJets__elm_);
   fChain->SetBranchAddress("PFJets_.mum_", PFJets__mum_, &b_PFJets__mum_);
   fChain->SetBranchAddress("PFJets_.ncand_", PFJets__ncand_, &b_PFJets__ncand_);
   fChain->SetBranchAddress("PFJets_.beta_", PFJets__beta_, &b_PFJets__beta_);
   fChain->SetBranchAddress("PFJets_.betaStar_", PFJets__betaStar_, &b_PFJets__betaStar_);
   fChain->SetBranchAddress("PFJets_.mpuTrk_", PFJets__mpuTrk_, &b_PFJets__mpuTrk_);
   fChain->SetBranchAddress("PFJets_.mlvTrk_", PFJets__mlvTrk_, &b_PFJets__mlvTrk_);
   fChain->SetBranchAddress("PFJets_.mjtTrk_", PFJets__mjtTrk_, &b_PFJets__mjtTrk_);
   fChain->SetBranchAddress("PFJets_.hof_", PFJets__hof_, &b_PFJets__hof_);
   fChain->SetBranchAddress("PFJets_.pujid_", PFJets__pujid_, &b_PFJets__pujid_);
   fChain->SetBranchAddress("PFJets_.calojetpt_", PFJets__calojetpt_, &b_PFJets__calojetpt_);
   fChain->SetBranchAddress("PFJets_.calojetef_", PFJets__calojetef_, &b_PFJets__calojetef_);
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
   fChain->SetBranchAddress("PFJetsCHS_.beta_", PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
   fChain->SetBranchAddress("PFJetsCHS_.betaStar_", PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
   fChain->SetBranchAddress("PFJetsCHS_.mpuTrk_", PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mlvTrk_", PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mjtTrk_", PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.hof_", PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
   fChain->SetBranchAddress("PFJetsCHS_.pujid_", PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
   fChain->SetBranchAddress("PFJetsCHS_.calojetpt_", PFJetsCHS__calojetpt_, &b_PFJetsCHS__calojetpt_);
   fChain->SetBranchAddress("PFJetsCHS_.calojetef_", PFJetsCHS__calojetef_, &b_PFJetsCHS__calojetef_);
   Notify();
}

Bool_t fillHistos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fillHistos::Show(Long64_t entry, Long64_t entry4)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
   // careful here, showing entry will also load it
   if (fChain2) fChain2->Show(entry4);
}
Int_t fillHistos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fillHistos_cxx
