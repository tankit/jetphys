// Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: June 2, 2015
// Updated: Aug 31, 2016
// Further updates: see git log

#define HistosFill_cxx
#include "HistosFill.h"


// Set the shortcuts for variables
HistosFill::HistosFill(TChain *tree) :
  pthat(EvtHdr__mPthat),
  weight(EvtHdr__mWeight),
  run(EvtHdr__mRun),
  evt(EvtHdr__mEvent),
  lbn(EvtHdr__mLumi),
  itpu(EvtHdr__mINTPU),
  ootpulate(EvtHdr__mOOTPULate),
  ootpuearly(EvtHdr__mOOTPUEarly),
  npv(EvtHdr__mNVtx),
  npvgood(EvtHdr__mNVtxGood),
  pvx(EvtHdr__mPVx),
  pvy(EvtHdr__mPVy),
  pvz(EvtHdr__mPVz),
  pvndof(EvtHdr__mPVndof),
  bsx(EvtHdr__mBSx),
  bsy(EvtHdr__mBSy),
  njt(PFJetsCHS__),
  gen_njt(GenJets__),
  rho(EvtHdr__mPFRho),
  met(PFMet__et_),
  metphi(PFMet__phi_),
  metsumet(PFMet__sumEt_)
{
  assert(tree);
  _initsuccess = Init(tree);
}


// Print the same info on a file and conditionally to the output
void HistosFill::PrintInfo(std::__cxx11::string info, bool printcout)
{
  *ferr << info << endl << flush;
  if (printcout) cout << info << endl << flush;
}


// Report memory info
void HistosFill::PrintMemInfo(bool printcout)
{
  gSystem->GetMemInfo(&_info);
  PrintInfo(Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                 _info.fMemTotal,  _info.fMemUsed,  _info.fMemFree,
                 _info.fSwapTotal, _info.fSwapUsed, _info.fSwapFree),printcout);
}


// Mostly setting up the root tree and its branches
bool HistosFill::Init(TTree *tree)
{
  ferr = new ofstream(Form("reports/HistosFill-%s.log",_jp_type),ios::out);

  if (!tree) {
    PrintInfo("No tree given, processing makes no sense.",true);
    return false; // With no tree, Loop will be interrupted
  }
  _outfile = new TFile(Form("output-%s-1.root",_jp_type), "NEW");
  if (!_outfile or _outfile->IsZombie()) {
    PrintInfo(Form("Opening the output file output-%s-1.root failed. Check if the file already exists.",_jp_type),false);
    return false;
  }

  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  worryHCALHotExcl = false;
  for (auto i = 0u; i<4; ++i) rangesHCALHotExcl[i] = 0;

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
  fChain->SetBranchAddress(Form("PFJets%s_.chm_",_jp_chs), PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
  fChain->SetBranchAddress(Form("PFJets%s_.nhm_",_jp_chs), PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
  fChain->SetBranchAddress(Form("PFJets%s_.phm_",_jp_chs), PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.elm_",_jp_chs), PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
  fChain->SetBranchAddress(Form("PFJets%s_.mum_",_jp_chs), PFJetsCHS__mum_, &b_PFJetsCHS__mum_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_hm_",_jp_chs), PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_phm_",_jp_chs), PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.ncand_",_jp_chs), PFJetsCHS__ncand_, &b_PFJetsCHS__ncand_);
  //fChain->SetBranchAddress("PFJets%s_.cm_",_jp_chs), PFJetsCHS__cm_, &b_PFJetsCHS__cm_);
  fChain->SetBranchAddress(Form("PFJets%s_.beta_",_jp_chs), PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
  fChain->SetBranchAddress(Form("PFJets%s_.betaStar_",_jp_chs), PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
  fChain->SetBranchAddress(Form("PFJets%s_.betaPrime_",_jp_chs), PFJetsCHS__betaPrime_, &b_PFJetsCHS__betaPrime_);
  fChain->SetBranchAddress(Form("PFJets%s_.mpuTrk_",_jp_chs), PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.mlvTrk_",_jp_chs), PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.mjtTrk_",_jp_chs), PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
  fChain->SetBranchAddress(Form("PFJets%s_.hof_",_jp_chs), PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
  fChain->SetBranchAddress(Form("PFJets%s_.pujid_",_jp_chs), PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
  fChain->SetBranchAddress(Form("PFJets%s_.calojetpt_",_jp_chs), PFJetsCHS__calojetpt_, &b_PFJetsCHS__calojetpt_);
  fChain->SetBranchAddress(Form("PFJets%s_.calojetef_",_jp_chs), PFJetsCHS__calojetef_, &b_PFJetsCHS__calojetef_);
  fChain->SetBranchAddress("genFlavour_", &genFlavour_, &b_events_genFlavour_);
  fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_, &b_events_genFlavourHadron_);

  if (_jp_quick) { // Activate only some branches
    fChain->SetBranchStatus("*",0);
    // Luminosity calculation
    if (_jp_ismc) fChain->SetBranchStatus("EvtHdr_.mPthat",1); // pthat
    if (_jp_ismc) fChain->SetBranchStatus("EvtHdr_.mWeight",1); // weight
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mRun",1); // run
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mEvent",1); // evt
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mLumi",1); // lbn

    // Event properties
    fChain->SetBranchStatus("EvtHdr_.mNVtx",1); // npv
    fChain->SetBranchStatus("EvtHdr_.mNVtxGood",1); // npvgood
    fChain->SetBranchStatus("EvtHdr_.mPFRho",1); // rho

    // Jet properties (jtpt, jte, jteta, jty, jtphi etc.)
    fChain->SetBranchStatus(Form("PFJets%s_",_jp_chs),1); // njt
    fChain->SetBranchStatus(Form("PFJets%s_.P4_*",_jp_chs),1); // jtp4*
    fChain->SetBranchStatus(Form("PFJets%s_.cor_",_jp_chs),1); // jtjes
    fChain->SetBranchStatus(Form("PFJets%s_.area_",_jp_chs),1); // jta

    if (_jp_ismc) {
      fChain->SetBranchStatus(Form("PFJets%s_.genP4_*",_jp_chs),1); // jtgenp4*
      fChain->SetBranchStatus(Form("PFJets%s_.genR_",_jp_chs),1); // jtgenr
    }

    // for quark/gluon study (Ozlem)
    fChain->SetBranchStatus(Form("PFJets%s_.QGtagger_",_jp_chs),1); // qgl
    if (_jp_ismc) fChain->SetBranchStatus(Form("PFJets%s_.partonFlavour_",_jp_chs),1);

    // Component fractions
    fChain->SetBranchStatus(Form("PFJets%s_.chf_",_jp_chs),1); // jtchf
    fChain->SetBranchStatus(Form("PFJets%s_.nemf_",_jp_chs),1); // jtnef
    fChain->SetBranchStatus(Form("PFJets%s_.nhf_",_jp_chs),1); // jtnhf
    fChain->SetBranchStatus(Form("PFJets%s_.cemf_",_jp_chs),1); // jtcef !!
    fChain->SetBranchStatus(Form("PFJets%s_.muf_",_jp_chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_hf_",_jp_chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_phf_",_jp_chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.ncand_",_jp_chs),1); // jtn
    fChain->SetBranchStatus(Form("PFJets%s_.beta_",_jp_chs),1); // jtbeta
    fChain->SetBranchStatus(Form("PFJets%s_.betaStar_",_jp_chs),1); // jtbetastar
    fChain->SetBranchStatus(Form("PFJets%s_.betaPrime_",_jp_chs),1); // jtbetaprime
    fChain->SetBranchStatus(Form("PFJets%s_.chm_",_jp_chs),1); // jtnch
    fChain->SetBranchStatus(Form("PFJets%s_.phm_",_jp_chs),1); // jtnne
    fChain->SetBranchStatus(Form("PFJets%s_.nhm_",_jp_chs),1); // jtnnh
    fChain->SetBranchStatus(Form("PFJets%s_.elm_",_jp_chs),1); // jtnce !!
    fChain->SetBranchStatus(Form("PFJets%s_.mum_",_jp_chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_hm_",_jp_chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_phm_",_jp_chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.tightID_",_jp_chs),1); // jtidtight
    fChain->SetBranchStatus(Form("PFJets%s_.looseID_",_jp_chs),1); // jtidloose

    //fChain->SetBranchStatus("rho",1);
    fChain->SetBranchStatus("PFMet_.et_",1); // met
    fChain->SetBranchStatus("PFMet_.phi_",1); // metphi
    fChain->SetBranchStatus("PFMet_.sumEt_",1); // metsumet

    fChain->SetBranchStatus("TriggerDecision_",1);
    fChain->SetBranchStatus("L1Prescale_",1);
    fChain->SetBranchStatus("HLTPrescale_",1);

    // Event cleaning
    //fChain->SetBranchStatus("pvrho",1);
    fChain->SetBranchStatus("EvtHdr_.mPVx",1); // pvx
    fChain->SetBranchStatus("EvtHdr_.mPVy",1); // pvy
    fChain->SetBranchStatus("EvtHdr_.mPVz",1); // pvz
    fChain->SetBranchStatus("EvtHdr_.mPVndof",1); // pvndof
    fChain->SetBranchStatus("EvtHdr_.mBSx",1); // bsx
    fChain->SetBranchStatus("EvtHdr_.mBSy",1); // bsy
    //

    if (_jp_ismc) {
      fChain->SetBranchStatus("EvtHdr_.mTrPu",1); // trpu
      fChain->SetBranchStatus("EvtHdr_.mINTPU",1); // itpu
      fChain->SetBranchStatus("EvtHdr_.mOOTPULate",1); // ootpulate
      fChain->SetBranchStatus("EvtHdr_.mOOTPUEarly",1); // ootpuearly
      fChain->SetBranchStatus("GenJets_",1); // gen_njt
      fChain->SetBranchStatus("GenJets_.fCoordinates.fX",1); // gen_jtp4x
      fChain->SetBranchStatus("GenJets_.fCoordinates.fY",1); // gen_jtp4y
      fChain->SetBranchStatus("GenJets_.fCoordinates.fZ",1); // gen_jtp4z
      fChain->SetBranchStatus("GenJets_.fCoordinates.fT",1); // gen_jtp4t
    }
  } else {
    fChain->SetBranchStatus("*",1);
  } // quick/slow

  // Set pointers to branches. NOTE: These are here so that the PFJetsCHS stuff etc. can be easily updated
  // if the file format changes.
  jtp4x = &PFJetsCHS__P4__fCoordinates_fX[0];
  jtp4y = &PFJetsCHS__P4__fCoordinates_fY[0];
  jtp4z = &PFJetsCHS__P4__fCoordinates_fZ[0];
  jtp4t = &PFJetsCHS__P4__fCoordinates_fT[0];
  jta = &PFJetsCHS__area_[0];
  jtjes = &PFJetsCHS__cor_[0];
  jtbeta = &PFJetsCHS__beta_[0];
  jtbetastar = &PFJetsCHS__betaStar_[0];
  jtbetaprime = &PFJetsCHS__betaPrime_[0];
  jtidloose = &PFJetsCHS__looseID_[0];
  jtidtight = &PFJetsCHS__tightID_[0];
  //
  jtgenr = &PFJetsCHS__genR_[0];
  jtgenp4x = &PFJetsCHS__genP4__fCoordinates_fX[0];
  jtgenp4y = &PFJetsCHS__genP4__fCoordinates_fY[0];
  jtgenp4z = &PFJetsCHS__genP4__fCoordinates_fZ[0];
  jtgenp4t = &PFJetsCHS__genP4__fCoordinates_fT[0];
  //

  // for quark/gluon study (Ozlem)
  qgl = &PFJetsCHS__QGtagger_[0];
  partonflavor = &PFJetsCHS__partonFlavour_[0];

  jtn = &PFJetsCHS__ncand_[0];
  jtnch = &PFJetsCHS__chm_[0];
  jtnnh = &PFJetsCHS__nhm_[0];
  jtnne = &PFJetsCHS__phm_[0];
  jtnce = &PFJetsCHS__elm_[0];
  jtnmu = &PFJetsCHS__mum_[0];
  jtnhh = &PFJetsCHS__hf_hm_[0];
  jtnhe = &PFJetsCHS__hf_phm_[0];

  jtchf = &PFJetsCHS__chf_[0];
  jtnhf = &PFJetsCHS__nhf_[0];
  jtnef = &PFJetsCHS__nemf_[0];
  jtcef = &PFJetsCHS__cemf_[0];
  jtmuf = &PFJetsCHS__muf_[0];
  jthhf = &PFJetsCHS__hf_hf_[0];
  jthef = &PFJetsCHS__hf_phf_[0];

  gen_jtp4x = &GenJets__fCoordinates_fX[0];
  gen_jtp4y = &GenJets__fCoordinates_fY[0];
  gen_jtp4z = &GenJets__fCoordinates_fZ[0];
  gen_jtp4t = &GenJets__fCoordinates_fT[0];

  return true;
}


// Loop over events
void HistosFill::Loop()
{
  if (!_initsuccess) return;

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("Starting Loop() initialization:");
  PrintMemInfo(true);

  if (!PreRun()) return;

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("Beginning Loop() proper:");
  PrintMemInfo(true);

  TStopwatch stop;
  stop.Start();
  TDatime bgn;

  ///////////////
  // Event loop
  ///////////////
  Long64_t nbytes = 0, nb = 0;
  Long64_t hopval = 1000000;
  Long64_t repval = 100001;
  for (Long64_t djentry=0; djentry<_nentries;djentry+=1+_jp_skim) { // Event loop
    _jentry = djentry+_nskip; // Add a shift from beginning
    if (LoadTree(_jentry) < 0) break;
    nb = fChain->GetEntry(_jentry);   nbytes += nb;

    int testval = djentry;
    if (_jp_skim>0) testval /= 1+_jp_skim;

    if (testval%hopval==repval) { // 1M report (first report timed to be early)
      // Report memory usage to avoid malloc problems when writing file
      PrintInfo( Form("Doing Loop(), %dM events:",int(_jentry/1e6 + 0.5)) );
      PrintMemInfo();
      cout << endl << Form("Processed %lld events (%1.1f%%) in %1.0f sec.",
                           _jentry-1, 100.*djentry/_nentries, stop.RealTime()) << endl;
      cout << "BGN: ";
      bgn.Print();
      TDatime now;
      cout << "NOW: ";
      now.Print();
      TDatime nxt;
      nxt.Set(nxt.Convert()+static_cast<UInt_t>(stop.RealTime()*(static_cast<Double_t>(hopval)/static_cast<Double_t>(testval))));
      cout << "NXT: ";
      nxt.Print();
      now.Set(now.Convert()+static_cast<UInt_t>(stop.RealTime()*(static_cast<Double_t>(_nentries)/static_cast<Double_t>(djentry)-1.0)));
      cout << "ETA: ";
      now.Print();

      if (_jp_save) {
        for (auto &manyhists : _histos)
          for (auto &onehist : manyhists.second)
            onehist->Write();
          if (_jp_doEtaHistos) {
            for (auto &manyhists : _etahistos)
              for (auto & onehist : manyhists.second)
                onehist->Write();
            if (_jp_doEtaHistosMcResponse) {
              for (auto &manyhists : _mchistos)
                for (auto & onehist : manyhists.second)
                  onehist->Write();
            }
          }
      }
      stop.Continue();
    }

    if (!AcceptEvent()) continue;

    if (_jp_debug) cout << "Histos are being filled!" << endl;
    // Here can categorize events into different triggers, epochs, topologies etc.
    // Eta and pT binning are handled in the FillSingleBasic class
    if (_jp_doBasicHistos) {
      FillBasic("Standard");
    }

    if (_jp_doEtaHistos and _pass) {
      FillEta("FullEta_Reco", jtpt, jteta, jtphi);
      if (_jp_ismc and _jp_doEtaHistosMcResponse) {
        FillEta("FullEta_Gen", jtgenpt, jtgeneta, jtgenphi);
        FillMC("FullEta_RecoPerGen_vReco", jtpt, jtgenpt, jtpt,    jteta,    jtphi);
        FillMC("FullEta_RecoPerGen_vGen",  jtpt, jtgenpt, jtgenpt, jtgeneta, jtgenphi);
      }
    }

    // Run quality checks
    if (_jp_isdt and _jp_doRunHistos) {
      FillRun("Runs");
      FillRun("RunsBarrel");
      FillRun("RunsTransition");
      FillRun("RunsEndcap");
    }
  } // for jentry
  cout << endl;

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("Finished processing %lld entries:",_nentries),true);
  PrintMemInfo(true);

  if (_jp_doRunHistos)   WriteRun();
  if (_jp_doEtaHistos)   WriteEta();
  if (_jp_ismc and _jp_doEtaHistos and _jp_doEtaHistosMcResponse) WriteMC();
  if (_jp_doBasicHistos) WriteBasic(); // this needs to be last, output file closed

  Report();

  stop.Stop();
  TDatime now;
  cout << "Stopping at: ";
  now.Print();
  PrintInfo(Form("Processing used %f s CPU time (%f h)",stop.CpuTime(),stop.CpuTime()/3600.),true);
  PrintInfo(Form("Processing used %f s real time (%f h)",stop.RealTime(),stop.RealTime()/3600.),true);
}


// Setup before event loop
bool HistosFill::PreRun()
{
  _nentries = fChain->GetEntriesFast();
  _ntot = fChain->GetEntries();
  _nskip = _jp_nskip;//0;
  _nentries = (_jp_nentries==-1 ? _ntot-_nskip : min(_ntot-_nskip, _jp_nentries));
  assert(_nskip+_nentries);

  _nbadevts_json = _nbadevts_dup = _nbadevts_run = _nbadevts_ls = _nbadevts_lum = 0;
  _bscounter_bad = _bscounter_good = 0;
  _ecalcounter_good = _ecalcounter_bad = 0;
  _rhocounter_good = _rhocounter_bad = 0;
  _trgcounter = _evtcounter = _totcounter = 0;
  _ecalveto = 0;

  // Initialize _pthatweight. It will be loaded for each tree separately.
  if (_jp_pthatbins)
    _pthatweight = 0;

  PrintInfo("\nCONFIGURATION DUMP:");
  PrintInfo("-------------------");
  PrintInfo(Form("Running over %sPF",_jp_algo));
  PrintInfo(Form("Running over %s",_jp_ismc?"MC":"data"));
  PrintInfo(Form("%s time-dependent JEC (IOV)",_jp_useIOV ? "Applying" : "Not applying"));
  PrintInfo(Form("%s jets in bad ECAL towers",_jp_doVetoECAL ? "Vetoing" : "Not vetoing"));
  PrintInfo(Form("%s Hot ECAL tower exclusion in eta-phi plane",_jp_doVetoECALHot ? "Doing" : "Not doing"));
  PrintInfo(Form("%s Hot HCAL tower exclusion in eta-phi plane",_jp_doVetoHCALHot ? "Doing" : "Not doing"));
  *ferr << endl;

  if (_jp_isdt) {
    PrintInfo(Form("%s additional JSON selection",_jp_dojson ? "Applying" : "Not applying"));
    PrintInfo(Form("%s additional prescalees on run/lumi.",_jp_doprescale ? "Applying" : "Not applying"));
    PrintInfo(Form("%s luminosities.",_jp_doprescale ? "Recalculating" : "Not recalculating"));
    PrintInfo(Form("%s additional run-level histograms",_jp_doRunHistos ? "Storing" : "Not storing"));
    PrintInfo(Form("%s basic set of histograms",_jp_doBasicHistos ? "Storing" : "Not storing"));
    PrintInfo(Form("%s histograms with a full eta-range",_jp_doEtaHistos ? "Storing" : "Not storing"));
  } else if (_jp_ismc) {
    PrintInfo(Form("%s pileup profile in MC to data",_jp_reweighPU ? "Reweighing" : "Not reweighing"));
    PrintInfo(Form("Processing %s samples", _jp_pthatbins ? "pThat binned" : "\"flat\""));
  }
  *ferr << endl;

  _JEC = 0;
  _L1RC = 0;
  _jecUnc = 0;

  // Time dependent JEC (only for dt)
  if (_jp_isdt and _jp_useIOV) {
    // If multiple IOV's are used, we set _JEC etc for each event separately, checking that the IOV is correct
    for (unsigned iovidx=0; iovidx<_jp_nIOV; ++iovidx)
      _iov.add(_jp_IOVnames[iovidx],_jp_IOVranges[iovidx][0],_jp_IOVranges[iovidx][1]);
  } else {
    // If only one great IOV is used, we can set _JEC etc. directly here.
    _iov.add("");
    bool setcorrection = _iov.setCorr(&_JEC,&_L1RC,&_jecUnc);
    if (!setcorrection or !_JEC or !_L1RC or (_jp_isdt and !_jecUnc)) {
      cout << "Issues while loading JEC; aborting..." << endl;
      return false;
    }
  } // JEC redone

  // Load latest JSON selection
  if (_jp_isdt and _jp_dojson and !LoadJSON(_jp_json)) {
    cout << "Issues loading the JSON file; aborting..." << endl;
    return false;
  }

  // Load PU profiles for MC reweighing
  if (_jp_ismc and _jp_reweighPU and !LoadPuProfiles(_jp_pudata, _jp_pumc)) {
    cout << "Issues loading the PU histograms for reweighting; aborting..." << endl;
    return false;
  }

  // Load prescale information to patch 76X
  if (_jp_isdt and _jp_doprescale and !LoadPrescales(_jp_prescalefile)) {
    cout << "Issues loading the prescale information; aborting..." << endl;
    return false;
  }

  // load ECAL veto file for cleaning data
  if (_jp_doVetoECAL and !LoadVetoECAL(_jp_fECALVeto)) {
    cout << "Issues loading the ECAL veto; aborting..." << endl;
    return false;
  }

  // load luminosity tables (prescales now stored in event)
  if (_jp_isdt and _jp_dolumi and !LoadLumi(_jp_lumifile)) {
    cout << "Issues loading the Lumi file; aborting..." << endl;
    return false;
  }

  if (_jp_ismc) cout << Form("Running on MC produced with %1.3g nb-1 (%lld evts)",
    1000. * _ntot / _jp_xsecMinBias, _ntot) << endl;
  if (_jp_isdt) cout << Form("Running on %lld events of data",_ntot) << endl;

  // Initialize histograms for different epochs and DQM selections
  if (_jp_doBasicHistos) {
    InitBasic("Standard");
  }

  if (_jp_doEtaHistos) {
    InitEta("FullEta_Reco");
    if (_jp_ismc) {
      InitEta("FullEta_Gen");
      if (_jp_doEtaHistosMcResponse) {
        InitMC("FullEta_RecoPerGen_vReco");
        InitMC("FullEta_RecoPerGen_vGen");
      }
    }
  }

  // For debugging purposes, save the weights used for pileup profiles.
  if (_jp_ismc and _jp_reweighPU) {
    _outfile->cd();
    _outfile->mkdir("puwgt");
    _outfile->cd("puwgt");
    for (auto &puprof : _pudist)
      puprof.second->Write();
  }

  if (_jp_isdt and _jp_doRunHistos) {
    InitRun("Runs",0.,3.);
    InitRun("RunsBarrel",0.,1.);
    InitRun("RunsTransition",1.,2.);
    InitRun("RunsEndcap",2.,3.);
  }

  if (_jp_doVetoECALHot) {
    string ECALHotTag = "";
    if (std::regex_search(_jp_run,std::regex("Run[BCD]16"))) ECALHotTag = "BCD";
    if (std::regex_search(_jp_run,std::regex("RunE16"))) ECALHotTag = "EF";
    if (std::regex_search(_jp_run,std::regex("RunFearly16"))) ECALHotTag = "EF";
    if (std::regex_search(_jp_run,std::regex("RunFlate16"))) ECALHotTag = "GH";
    if (std::regex_search(_jp_run,std::regex("Run[GH]16"))) ECALHotTag = "GH";
    assert(ECALHotTag!="");
    fECALHotExcl = new TFile(Form("rootfiles/hotjets-run%s.root",ECALHotTag.c_str()),"READ");
    assert(fECALHotExcl and !fECALHotExcl->IsZombie() && Form("file rootfiles/hotjets-run%s.root missing",ECALHotTag.c_str()));
    h2ECALHotExcl = (TH2D*)fECALHotExcl->Get(Form("h2hot%s",_jp_ECALHotType));
    assert(h2ECALHotExcl and "erroneous eta-phi exclusion type");
    PrintInfo(Form("Loading ECAL corrections rootfiles/hotjets-run%s.root with h2hot %s",ECALHotTag.c_str(),_jp_ECALHotType));
  }

  if (_jp_dotrpufile) {
    // Map of mu per (run,lumi)
    TFile *fmu = new TFile(_jp_trpufile,"READ");
    assert(fmu && !fmu->IsZombie());
    _h2mu = (TH2F*)fmu->Get("hLSvsRUNxMU");
    assert(_h2mu and !_h2mu->IsZombie());
  }

  // Qgl: load quark/gluon probability histos (Ozlem)
  // 1. open the previous output-MC-1_iteration1.root in the beginning
  // 2. get the hqgl_q and hqgl_g for each eta bin, store in array
  // 3. find correct hqgl_q and hqgl_g from array (normalized)
  // 4. calculate probg = g / (q+g)
  if (_jp_doqglfile) {
    TFile *finmc = new TFile(_jp_qglfile,"READ");
    assert(finmc && !finmc->IsZombie());

    double veta[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7};
    const int neta = sizeof(veta)/sizeof(veta[0])-1;
    // same or vpt, vtrigpt, npt
    const int npt = _jp_notrigs;
    const double *vtrigpt = &_jp_trigthr[0];
    double vpt[npt+1]; vpt[0] = 0;
    for (int trigidx = 0; trigidx != npt; ++trigidx) vpt[trigidx+1] = _jp_trigranges[trigidx][1];
    const int nqgl = 101; double vqgl[nqgl+1];
    for (int qglidx = 0; qglidx != nqgl+1; ++qglidx) vqgl[qglidx] = 0. + 0.01*qglidx;

    _h3probg = new TH3D("_h3probg","Gluon prob.;#eta_{jet};p_{T,jet};QGL",
                        //6,0,3.0, 10,0,1000, 101,0,1.01);
                        neta,veta, npt,vpt, nqgl,vqgl); //101,0,1.01);

    // Loop over pt and eta bins in the gql histos
    for (int ieta = 1; ieta != _h3probg->GetNbinsX()+1; ++ieta) {
      for (int ipt = 1; ipt != _h3probg->GetNbinsY()+1; ++ipt) {
        // Gluons in the given eta/pt slice
        string sg = Form("Standard/Eta_%1.1f-%1.1f/jt%1.0f/hqgl_g",
                          veta[ieta-1], veta[ieta], vtrigpt[ipt-1]);
        cout << sg << endl << flush;
        TH1D *hqgl_g = dynamic_cast<TH1D*>(finmc->Get(sg.c_str())); assert(hqgl_g);
        hqgl_g->Scale(1./hqgl_g->Integral());

        // Quarks in the given eta/pt slice
        string sq = Form("Standard/Eta_%1.1f-%1.1f/jt%1.0f/hqgl_q", veta[ieta-1], veta[ieta], vtrigpt[ipt-1]);
        TH1D *hqgl_q = dynamic_cast<TH1D*>(finmc->Get(sq.c_str())); assert(hqgl_q);
        hqgl_q->Scale(1./hqgl_q->Integral());

        // Loop over qgl indices
        for (int iqgl = 1; iqgl != _h3probg->GetNbinsZ()+1; ++iqgl) {
          // Fetch gluon probability
          double probg = hqgl_g->GetBinContent(iqgl) /
          (hqgl_g->GetBinContent(iqgl) + hqgl_q->GetBinContent(iqgl));
          _h3probg->SetBinContent(ieta, ipt, iqgl, probg);
        } // for iqgl
      } // for ipt
    } // for ieta
  }
  PrintInfo("Finished pre-run processing!",true);
  return true;
}


// Routines and selections before histograms are filled
bool HistosFill::AcceptEvent()
{
  if (_jp_isdt) { // For DT fetch true pileup from the json or histogram info
    trpu = _avgpu[run][lbn];
    if (trpu==0 and _jp_dotrpufile and _h2mu) { // This is rarely used, as no new histogram is provided
      int irun = _h2mu->GetXaxis()->FindBin(run);
      int ilbn = _h2mu->GetYaxis()->FindBin(lbn);
      trpu = _h2mu->GetBinContent(irun, ilbn);
    }
  } else { // For MC the info exists in the SMPJ tree
    trpu = EvtHdr__mTrPu;
  }

  if (njt > _maxnjt) {
    PrintInfo(Form("Array overflow: njt = %d > njtmax= %d",njt,_maxnjt),true);
    assert(njt<_maxnjt);
  }

  if (_jp_ismc and gen_njt > _maxnjt) {
    PrintInfo(Form("Array overflow: gen_njt = %d > njtmax= %d",gen_njt,_maxnjt),true);
    assert(gen_njt<_maxnjt);
  }

  if (_jp_debug) {
    cout << endl << flush;
    Show(_jentry);
    cout << endl << endl << flush;

    cout << "***Checking basic event variables are read out:" << endl;
    cout << "isdata = " << _jp_isdt << " / ismc = " << _jp_ismc << endl;
    cout << "trpu = " << trpu << endl;
    cout << "pthat = " << pthat << endl;
    cout << "weight = " << weight << endl;
    cout << "njt = " << njt << endl;
    cout << "idloose[0] = " << (njt>0 ? jtidloose[0] : -1) << endl;
    cout << "idtight[0] = " << (njt>0 ? jtidtight[0] : -1) << endl;
    cout << "***end basic event variables" << endl;
    cout << endl << flush;
  }

  // Check if duplicate
  if (_jp_isdt and _jp_checkduplicates) {
    set<int>& events = _duplicates[run][lbn];
    if (events.find(evt)!=events.end()) {
      ++_nbadevts_dup;
      return false;
    }
    events.insert(evt);
  }
  ++_cnt["01all"];

  // Check if good run/LS, including JSON selection
  if (_jp_isdt) {
    if (_jp_dojson) {
      // Does the run/LS pass the latest JSON selection?
      if (_json[run][lbn]==0) {
        _badjson.insert(pair<int, int>(run, lbn));
        ++_nbadevts_json;
        return false;
      }
    } // _jp_dojson

    if (_jp_dolumi) {
      // Do we have the run listed in the .csv file?
      auto irun = _lums.find(run);
      if (irun==_lums.end()) {
        _badruns.insert(run);
        ++_nbadevts_run;
        return false;
      }

      // Do we have the LS listed in the .csv file?
      auto ils = irun->second.find(lbn);
      if (ils==irun->second.end()) {
        _badlums.insert(pair<int, int>(run,lbn));
        ++_nbadevts_ls;
        return false;
      }

      // Does the .csv file list a non-zero luminosity?
      if (ils->second==0) {
        _nolums.insert(pair<int, int>(run, lbn));
        ++_nbadevts_lum;
        //continue; // Could be Poisson fluctuation to zero
      }
    } // _jp_dolumi
  } // _jp_isdt

  // Keep track of LBNs
  _passedlumis.insert(pair<int, int>(run, lbn));

  ++_cnt["02ls"];

  // Reset event ID
  _pass = true;

  // Reject events with no vertex
  pvrho = tools::oplus(pvx, pvy);
  _pass = _pass and npvgood>0 and pvrho<2.;
  if (_pass) ++_cnt["03vtx"];

  // Event cuts against beam backgrounds
  if (_pass) {
    if (tools::oplus(pvx-bsx, pvy-bsy)>0.15 or pvndof<=4 or fabs(pvz) >= 24.) {
      ++_bscounter_bad;
      _pass = false;
    } else {
      ++_bscounter_good;
      ++_cnt["04bsc"];
    }
  }

  // load correct IOV for JEC
  if (_jp_isdt and _jp_useIOV) {
    bool setcorrection = _iov.setCorr(&_JEC,&_L1RC,&_jecUnc,run);
    if (!setcorrection or !_JEC or !_L1RC or !_jecUnc) {
      cout << "Issues while loading JEC; aborting..." << endl;
      assert(false);
    }
  }

  if (_jp_debug) cout << "JEC and MET calculation and leading jets info!" << endl;
  // Calculate jec and propagate jec to MET 1 and MET 2
  double mex = met * cos(metphi);
  double mey = met * sin(metphi);
  double ucx = mex;
  double ucy = mey;
  // Find leading jets (residual JEC may change ordering)
  // CAUTION: for less than 3 jets, we input -1 on the place of the index
  // NOTE: Above this we only use the initial corrections. This serves as a good approximation,
  // as the values before reapplying JEC are more like the ones we knew when collecting data.
  for (int i = 0; i<3; ++i)
    jt3leads[i] = -1;
  ///////////////
  // Jet loop
  ///////////////
  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    p4.SetPxPyPzE(jtp4x[jetidx],jtp4y[jetidx],jtp4z[jetidx],jtp4t[jetidx]);
    // Divide by the original JES
    if (_jp_debug) cout << "Entering jet loop!" << endl;
    if (_jp_undojes) p4 *= 1/jtjes[jetidx];

    jtptu[jetidx] = p4.Pt();
    jteu[jetidx] = p4.E();

    if (_jp_debug) cout << "Recalculating JEC!" << endl;
    // Recalculate JEC
    _JEC->setRho(rho);
    _JEC->setNPV(npvgood);
    _JEC->setJetA(jta[jetidx]);
    _JEC->setJetPt(jtptu[jetidx]);
    _JEC->setJetE(jteu[jetidx]);
    _JEC->setJetEta(p4.Eta());
    jtjesnew[jetidx] = _JEC->getCorrection();

    if (_jp_debug) cout << "Recalculating JEC (again)!" << endl;
    // Recalculate JEC (again to get subcorrections)
    _JEC->setRho(rho);
    _JEC->setNPV(npvgood);
    _JEC->setJetA(jta[jetidx]);
    _JEC->setJetPt(jtptu[jetidx]);
    _JEC->setJetE(jteu[jetidx]);
    _JEC->setJetEta(p4.Eta());
    //
    vector<float> v = _JEC->getSubCorrections();
    double jec_res = 1;
    if (_jp_ismc or _jp_skipl2l3res) {
      assert(v.size()==3);
    } else {
      assert(v.size()==4);
      jec_res = v[3]/v[2];
    }
    double jec_l1 = v[0];
    double jec_l2l3 = v[2]/v[0];
    jtjes_l1[jetidx] = jec_l1;
    jtjes_l2l3[jetidx] = jec_l2l3;
    jtjes_res[jetidx] = jec_res;
    assert(jtjesnew[jetidx] == v[v.size()-1]);

    if (_jp_debug) cout << "Reapplying JEC!" << endl;
    if (_jp_redojes) p4 *= jtjesnew[jetidx];

    jte[jetidx] = p4.E();
    jtpt[jetidx] = p4.Pt();
    jteta[jetidx] = p4.Eta();
    jtphi[jetidx] = p4.Phi();
    jty[jetidx] = p4.Rapidity();

    if (_jp_debug) cout << "Gen info!" << endl;
    // Calculate gen level info
    if (_jp_ismc) {
      gp4.SetPxPyPzE(jtgenp4x[jetidx],jtgenp4y[jetidx],jtgenp4z[jetidx],jtgenp4t[jetidx]);
      jtgenpt[jetidx] = gp4.Pt();
      jtgeny[jetidx] = gp4.Rapidity();
      jtgeneta[jetidx] = gp4.Eta();
      jtgenphi[jetidx] = gp4.Phi();
    }

    if (_jp_debug) cout << "Jet " << jetidx << " corrected!" << endl;

    // Only use jets with corr. pT>recopt GeV to equalize data and MC thresholds
    if (jtpt[jetidx] > _jp_recopt and fabs(jteta[jetidx])<4.7) {
      // MET 1: the one where JEC is applied. MET1 needs to be recalculated as JEC changes.
      // Subtract uncorrected jet pT from met, put back corrected & add L1RC offset to keep PU isotropic.
      _L1RC->setRho(rho);
      _L1RC->setJetA(jta[jetidx]);
      _L1RC->setJetPt(jtptu[jetidx]);
      _L1RC->setJetE(jteu[jetidx]);
      _L1RC->setJetEta(jteta[jetidx]);
      double l1corr = _L1RC->getCorrection();
      double dpt = - jtpt[jetidx] + l1corr*jtptu[jetidx];
      //double dpt = - jtpt[jetidx] + (l1chs - l1pf + l1corr)*jtptu[jetidx];
      mex += dpt * cos(jtphi[jetidx]);
      mey += dpt * sin(jtphi[jetidx]);

      // MET 2: record unclustered energy.
      // Keep track of remaining pT in unclustered energy, add to MET l1corr jets, from which ue is substracted.
      // Effectively this means substracting jets (without their PU and UE) from MET (=> homogeneous background).
      double ue = 1.068 * jta[jetidx]; // CAUTION: One should check that the magical coefficient here is good.
      double dptu = -ue + l1corr*jtptu[jetidx];
      ucx += dptu * cos(jtphi[jetidx]);
      ucy += dptu * sin(jtphi[jetidx]);
    }

    if (jt3leads[0]==-1 or jtpt[jt3leads[0]]<jtpt[jetidx]) {
      jt3leads[2] = jt3leads[1];
      jt3leads[1] = jt3leads[0];
      jt3leads[0] = jetidx;
    } else if (jt3leads[1]==-1 or jtpt[jt3leads[1]]<jtpt[jetidx]) {
      jt3leads[2] = jt3leads[1];
      jt3leads[1] = jetidx;
    } else if (jt3leads[2]==-1 or jtpt[jt3leads[2]]<jtpt[jetidx]) {
      jt3leads[2] = jetidx;
    }
  } // for jetidx

  // Type I MET (this is the best one we've got; works optimally if we keep T0Txy MET as the raw MET and apply here the newest JEC)
  // met1 = -jtpt[jetidx]+l1corr*jtptu[jetidx]+metraw
  met1 = tools::oplus(mex, mey);
  metphi1 = atan2(mey, mex);

  // Correct unclustered energy; jec for 10 GeV jets varies between 1.1-1.22 at |y|<2.5,
  // 2.5-3.0 even goes up to 1.35 => assume 1.15  => try 1.5 => to 1.25 (high pT threshold on jets)
  mex += 0.25*ucx;
  mey += 0.25*ucy;
  // Type II MET witch C = 1.25 (This is not recommended for pfJets.
  // met2 = met1 - C*uncl.ptsum = met1 + 0.25*(l1corr*jtptu[jetidx]-ue[jetidx]) = 1.25*metraw + 1.25*l1corr*jtptu[jetidx] - jtpt[jetidx] - 0.25ue[jetidx]
  met2 = tools::oplus(mex, mey);
  metphi2 = atan2(mey, mex);

  // The leading indices
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];

  // Zero jets not sensible
  _pass = _pass and i0>=0;
  if (_pass) ++_cnt["04njt"];
  else _pass = false;

  if (_jp_debug) {
    cout << "Indices for the three leading jets: " << jt3leads[0] << " " << jt3leads[1] << " " << jt3leads[2] << endl;
    cout << "Gen flav calculation!" << endl;
  }

  // Event cuts against beam backgrounds
  if (_pass and _ecalveto) {
    if ( (njt>=1 and _ecalveto->GetBinContent(_ecalveto->FindBin(jteta[i0],jtphi[i0]))!=0)
      or (njt>=2 and _ecalveto->GetBinContent(_ecalveto->FindBin(jteta[i1],jtphi[i1]))!=0) ) {
      ++_ecalcounter_bad;
      _pass = false;
    } else {
      ++_ecalcounter_good;
      ++_cnt["05ecal"];
    }
  } // ecal veto

  // Check rho
  if (_pass) {
    if (rho>40.) {
      ++_rhocounter_bad;
      _pass = false;
      if (_jp_debug)
        cout << Form("\nrun:ev:ls %d:%d:%d : rho=%1.1f njt=%d npv=%d jtpt0=%1.1f sumet=%1.1f met=%1.1f\n",
                    run, lbn, evt, rho, njt, npv, (njt>0 ? jtpt[i0] :0.), metsumet, met) << flush;
    } else {
      ++_rhocounter_good;
      ++_cnt["06rho"];
    }
  }

  // Reset prescales (dynamic can change within run)
  for (auto &scaleit : _prescales)
    scaleit.second[run] = 0;

  // Fill trigger information
  _trigs.clear();

  // Simulate other triggers for MC, if so wished
  // (this is slow, though)
  if (_jp_ismc) {
    // Always insert the generic mc trigger
    if (_jp_debug) cout << "Entering PU weight calculation!" << endl;
    if (_jp_domctrigsim and njt>0) {
      // Only add the greatest trigger present
      // Calculate trigger PU weight
      bool found = false;
      double wtrue = 1.0;
      for (int itrg = _jp_notrigs-1; itrg >= 0; --itrg) {
        if (jtpt[i0]>_jp_trigranges[itrg][0]) {
          string trg_name = _jp_triggers[itrg];
          _trigs.insert(trg_name);
          _wt[trg_name] = 1.;

          // Reweight in-time pile-up
          if (_jp_reweighPU) {
            int k = _pudist[trg_name]->FindBin(trpu);
            wtrue = _pudist[trg_name]->GetBinContent(k);
            _wt[trg_name] *= wtrue;
          }
          found = true;
          break; // Don't add lesser triggers
        }
      }
      _pass = _pass and found;
      if (_pass) ++_cnt["07mctrg"];
      else return false; // Leading jet is weak

      // check for non-zero PU weight
      _pass = _pass and wtrue!=0;
      if (_pass) ++_cnt["07puw"];
      else return false; // Bad pu areas with zero weight excluded
    } // _jp_domctrigsim
    _trigs.insert("mc");
    _wt["mc"] = 1.0;
    if (_jp_reweighPU) {
      int k = _pudist[_jp_reftrig]->FindBin(trpu);
      _wt["mc"] *= _pudist[_jp_reftrig]->GetBinContent(k);
    }
  } else if (_jp_isdt) {
    // For data, check trigger bits
    if (_jp_debug) {
      cout << "TriggerDecision_.size()=="<<TriggerDecision_.size()<<endl<<flush;
      cout << "_availTrigs.size()=="<<_availTrigs.size()<<endl<<flush;
      cout << "_goodTrigs.size()=="<<_goodTrigs.size()<<endl<<flush;
    }
    assert(TriggerDecision_.size() == _availTrigs.size());

    for (auto itrg = 0u; itrg != _jp_notrigs; ++itrg)
      _wt[_jp_triggers[itrg]] = 1.0;

    for (auto goodIdx = 0u; goodIdx < _goodTrigs.size(); ++goodIdx) {
      auto &itrg = _goodTrigs[goodIdx];
      if (TriggerDecision_[itrg]!=1) continue; // -1, 0, 1

      string trigname = _availTrigs[itrg];
      if (_jp_debug and TriggerDecision_[itrg]>0)
        cout << trigname << " " << TriggerDecision_[itrg]
             << " " << L1Prescale_[itrg] << " " << HLTPrescale_[itrg] << endl;

      // Set prescale from event for now
      //if (L1Prescale_[itrg]>0 and HLTPrescale_[itrg]>0) { There's trouble in 2017 L1, so we let it pass
      if (HLTPrescale_[itrg]>0) {
        double l1 = L1Prescale_[itrg];
        if (l1==0) l1 = 1;
        _prescales[trigname][run] = l1 * HLTPrescale_[itrg];
      } else {
        cout << "Error for trigger " << trigname << " prescales: "
              << "L1  =" << L1Prescale_[itrg]
              << "HLT =" << HLTPrescale_[itrg] << endl;
        _prescales[trigname][run] = 0;
        if (_jp_debug) { // check prescale
          double prescale = _prescales[trigname][run];
          if (L1Prescale_[itrg]*HLTPrescale_[itrg]!=prescale) {
            cout << "Trigger " << trigname << ", "
            << "Prescale(txt file) = " << prescale << endl;
            cout << "L1 = " << L1Prescale_[itrg] << ", "
            << "HLT = " << HLTPrescale_[itrg] << endl;
            assert(false);
          }
        } // debug
      }

      if (_prescales[trigname][run]>0) {
        // Set trigger only if prescale information is known
        _trigs.insert(trigname);
        _wt[trigname] = _goodWgts[goodIdx];
      } else {
        // Make sure all info is good! This is crucial if there is something odd with the tuples
        PrintInfo(Form("Missing prescale for %s in run %d",trigname.c_str(),run));
      }
    } // for itrg (_goodTrigs)

  }

  ++_totcounter;
  if (_pass) ++_evtcounter;
  _pass = _pass and _trigs.size()>0;
  if (_pass) {
    ++_trgcounter;
    if (_jp_isdt) ++_cnt["07trg"];
  }

  // Retrieve event weight. _w0 is static, _w is chanching with the trigger
  _w0 = 1.0;
  if (_jp_ismc) {
    _w0 *= weight;
    if (_jp_pthatbins)
      _w0 *= _pthatweight;
    assert(_w0>0);
  }
  _w = _w0;

  // TODO: implement reweighing for k-factor (NLO*NP/LOMC)

  if (_jp_ismc) {
    ///////////////
    // Gen Jet loop
    ///////////////
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      genp4.SetPxPyPzE(gen_jtp4x[gjetidx],gen_jtp4y[gjetidx],gen_jtp4z[gjetidx],gen_jtp4t[gjetidx]);
      gen_jtpt[gjetidx] = genp4.Pt();
      gen_jteta[gjetidx] = genp4.Eta(); // for matching
      gen_jtphi[gjetidx] = genp4.Phi(); // for matching
      gen_jty[gjetidx] = genp4.Rapidity();

      // Ozlem: loop for finding partonflavor by matching genjets and jets
      int ireco = -1;
      for (int j = 0; j != njt && ireco == -1; ++j) {
        p4.SetPxPyPzE(jtp4x[j],jtp4y[j],jtp4z[j],jtp4t[j]);
        if (genp4.DeltaR(p4) < 0.4)
          ireco = j;
      } // for j
      if (ireco!=-1)
        gen_partonflavor[gjetidx] = Int_t(partonflavor[ireco]+0.25);
      else
        gen_partonflavor[gjetidx] = -1;
    } // for gjetidx
  } // MC

  _jetids.resize(njt);
  for (unsigned int jetid = 0; jetid != _jetids.size(); ++jetid)
    _jetids[jetid] = true;
  FillJetID(_jetids);

  if (_jp_isdt and _jp_doVetoHCALHot and worryHCALHotExcl) {
    bool good0 = jteta[i0] < rangesHCALHotExcl[0] or jteta[i0] > rangesHCALHotExcl[1] or
                 jtphi[i0] < rangesHCALHotExcl[2] or jtphi[i0] > rangesHCALHotExcl[3];
    bool good1 = jteta[i1] < rangesHCALHotExcl[0] or jteta[i1] > rangesHCALHotExcl[1] or
                 jtphi[i1] < rangesHCALHotExcl[2] or jtphi[i1] > rangesHCALHotExcl[3];
    _pass = _pass and good0 and good1;
    if (_pass) ++_cnt["08etaphiexclHCAL"];
  }

  if (_jp_isdt and _jp_doVetoECALHot) {
    // Abort if one of the leading jets is in a difficult zone
    bool good0 = h2ECALHotExcl->GetBinContent(h2ECALHotExcl->FindBin(jteta[i0],jtphi[i0])) < 0;
    bool good1 = h2ECALHotExcl->GetBinContent(h2ECALHotExcl->FindBin(jteta[i1],jtphi[i1])) < 0;
    _pass = _pass and good0 and good1;
    if (_pass) ++_cnt["08etaphiexclECAL"];
  }

  if (_jetids[i0] and _pass) {
    // Check if overweight PU event
    if (_jp_ismc and _pass) {
      if (jtpt[i0] < 1.5*jtgenpt[i0]) ++_cnt["09ptgenlim"];
      else _pass = false;

      if (_pass) {
        if (_jp_pthatbins) {
          if (jtpt[i0] < _pthatuplim) ++_cnt["10pthatlim"];
          else _pass = false;
        } else { // Flat case
          if (jtpt[i0] < 1.5*pthat) ++_cnt["10pthatlim"];
          else _pass = false;
        }
      }
    }
    if (_pass) ++_cnt["11jtid"]; // Non-restrictive
  }

  // Equipped in FillBasic and FillRun
  _pass_qcdmet = (met < 0.4 * metsumet || met < 45.); // QCD-11-004

  return true;
}


// Report event stuff
void HistosFill::Report()
{
  // List bad runs
  PrintInfo(Form("Processed %d events in total",_totcounter));
  PrintInfo(Form("Processed %d events passing basic data quality and trigger cuts",_trgcounter));
  PrintInfo(Form("(out of %d passing data quality cuts)",_evtcounter));
  if (_jp_dojson)
    PrintInfo(Form("Found %d bad events according to new JSON (events cut)",_nbadevts_json));
  if (_jp_dolumi) {
    PrintInfo(Form("Found %d bad events according to lumi file.",_nbadevts_run));
    if (_badruns.size()>0) {
      *ferr << "The found " << _badruns.size() << " bad runs:";
      for (auto &runit : _badruns)
        *ferr << " " << runit;
      *ferr << endl;
    }
    PrintInfo(Form("Found %lu bad LS and %lu non-normalizable LS in good runs",_badlums.size(),_nolums.size()));
    PrintInfo(Form("These contained %d discarded events in bad LS and %d in non-normalizable LS",_nbadevts_ls,_nbadevts_lum));
    *ferr << endl;
  }
  if (_jp_checkduplicates)
    PrintInfo(Form("Found %d duplicate events, which were properly discarded",_nbadevts_dup));

  // Report beam spot cut efficiency
  PrintInfo(Form("Beam spot counter discarded %d events out of %d (%f %%)\nBeam spot expectation is less than 0.5%%",
                 _bscounter_bad,_bscounter_good,double(_bscounter_bad)/double(_bscounter_good)*100.));

  // Report ECAL hole veto efficiency
  if (_jp_doVetoECAL)
    PrintInfo(Form("ECAL hole veto counter discarded %d events out of %d (%f %%)\n"
                   "ECAL hole expectation is less than 2.6%% [=2*57/(60*72)]",
                   _ecalcounter_bad,_ecalcounter_good,double(_ecalcounter_bad)/double(_ecalcounter_good)*100.));

  // Report rho veto efficiency
  PrintInfo(Form("Rho<40 veto counter discarded %d events out of %d (%f %%)\nRho veto expectation is less than 1 ppm",
                _rhocounter_bad,_rhocounter_good,double(_rhocounter_bad)/double(_rhocounter_good)*100.));
  for (auto &cit : _cnt)
    PrintInfo(Form("%s: %d (%1.1f%%)",cit.first.c_str(),cit.second,100.*cit.second/max(1,_cnt["01all"])),true);

  if (_jp_isdt) PrintInfo("Note that for DT it is likely that we lose a large percentage of events for the trigger."
                          "Events triggered by JetHT and AK8PFJet are included in addition to AK4PFJet.");
  *ferr << endl;

  PrintInfo("Reporting lumis not discraded in reports/passedlumis.json",true);
  ofstream fout(Form("reports/passedlumis-%s.json",_jp_type), ios::out);
  for (auto &lumit : _passedlumis) fout << lumit.first << " " << lumit.second << endl;

  if (_jp_dojson and _badjson.size()>0) {
    PrintInfo("Reporting lumis discarded by json selection (_jp_dojson) in reports/badlumis_json",true);
    ofstream fout2(Form("reports/badlumis_json-%s.json",_jp_type), ios::out);
    for (auto &jsit : _badjson) fout2 << jsit.first << " " << jsit.second << endl;
  } // _jp_dojson

  if (_jp_dolumi) {
    if (_badlums.size()>0) {
      PrintInfo("Reporting lumis discarded by lumifile selection (_jp_dolumi) in reports/badlumis_lumi.json",true);
      ofstream fout2(Form("reports/badlumis_lumi-%s.json",_jp_type), ios::out);
      for (auto &jsit : _nolums) fout2 << jsit.first << " " << jsit.second << endl;
    } // _badlums
    if (_nolums.size()>0) {
      PrintInfo("Reporting lumis discarded fby lumifile selection (_jp_dolumi) @ zero luminosity in reports/badlumis_zerolumi.json",true);
      ofstream fout2(Form("reports/badlumis_zerolumi-%s.json",_jp_type), ios::out);
      for (auto &jsit : _nolums) fout2 << jsit.first << " " << jsit.second << endl;
    } // _nolums
  } // _jp_dolumi
}


// Initialize basic histograms for trigger and eta bins
void HistosFill::InitBasic(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitBasic(%s):",name.c_str()));
  PrintMemInfo();

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = _outfile;
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  bool enteroutdir = f->cd(name.c_str());
  assert(enteroutdir);
  //TDirectory *topdir = gDirectory;
  TDirectory *topdir = f->GetDirectory(name.c_str()); assert(topdir);
  topdir->cd();

  // Pseudorapidity bins + HF + barrel
  double etas[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 0., 1.3};
  const int netas = sizeof(etas)/sizeof(etas[0])-1;

  // define triggers
  vector<string> triggers;
  // define efficient pT ranges for triggers for control plots
  map<string, pair<double, double> > pt;
  // define pT values for triggers
  map<string, double> pttrg;
  if (_jp_ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(_jp_recopt, _jp_emax);
    pttrg["mc"] = _jp_recopt;
  }
  if (_jp_isdt or _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_notrigs; ++itrg) {
      string trg = _jp_triggers[itrg];
      triggers.push_back(trg);
      double pt1 = _jp_trigranges[itrg][0];
      double pt2 = _jp_trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = _jp_trigthr[itrg];
      pttrg[trg] = pt0;
    }
  }

  // Loop over pseudorapidity, trigger bins
  for (int etaidx = 0; etaidx != netas; ++etaidx) {
    if (etas[etaidx+1] > etas[etaidx]) { // create real bins only
      // subdirectory for rapidity bin
      const char *etaname = Form("Eta_%1.1f-%1.1f", etas[etaidx], etas[etaidx+1]);
      assert(topdir);
      //assert(topdir->mkdir(yname));
      topdir->mkdir(etaname);
      bool enteretadir = topdir->cd(etaname);
      assert(enteretadir);
      //TDirectory *ydir = gDirectory;
      TDirectory *etadir = topdir->GetDirectory(etaname); assert(etadir);
      etadir->cd();

      for (unsigned int j = 0; j != triggers.size(); ++j) {
        // subdirectory for trigger
        const char *trg = triggers[j].c_str();
        assert(etadir);
        //assert(ydir->mkdir(trg));
        etadir->mkdir(trg);
        bool entertrgdir = etadir->cd(trg);
        assert(entertrgdir);
        //TDirectory *dir = gDirectory;
        TDirectory *dir = etadir->GetDirectory(trg); assert(dir);
        dir->cd();

        // Initialize and store
        assert(dir);
        HistosBasic *h = new HistosBasic(dir, trg, etas[etaidx], etas[etaidx+1], pttrg[trg],
                                         pt[trg].first, pt[trg].second, triggers[j]=="mc");
        _histos[name].push_back(h);
      } // for j
    } // real bin
  } // for etaidx

  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitBasic(%s) finished:",name.c_str()));
  PrintMemInfo();
} // initBasic


// Loop over basic histogram containers to fill all
void HistosFill::FillBasic(string name)
{
  for (unsigned int histidx = 0; histidx != _histos[name].size(); ++histidx)
    FillSingleBasic(_histos[name][histidx]);
}


// Fill basic histograms after applying pt, y cuts
void HistosFill::FillSingleBasic(HistosBasic *h)
{
  assert(h);
  h->hpttmp->Reset();
  h->hpt_tmp->Reset();
  if (h->ismcdir) {
    assert(h->hpt_g0_tmp);
    h->hpt_g0_tmp->Reset();
  }

  _w = _w0 * _wt[h->trigname];
  if (_w <= 0) return;

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  // Luminosity information
  if (_jp_isdt and h->lums[run][lbn]==0) {
    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    double prescale(0);
    map<int, int>::const_iterator ip = _prescales[h->trigname].find(run);
    if (ip==_prescales[h->trigname].end()) {
      if (fired) {
        PrintInfo(Form("No prescale info for trigger %s in run %d!",
                       h->trigname.c_str(),run));
        assert(false);
      }
    } else prescale = ip->second;

    if (prescale==0 and fired) {
      PrintInfo(Form("Prescale zero for trigger %s in run %d!",
                     h->trigname.c_str(),run));
      prescale = 1.;
      assert(false);
    }

    h->lumsum += (prescale ? lum / prescale : 0.);
    h->lumsum2 += (prescale ? lum2 / prescale : 0.);
    h->lums[run][lbn] = 1;

    h->hlumi_vstrpu->Fill(trpu, prescale ? lum / prescale : 0.);
  }

  // For MC vs truePU
  if (_jp_ismc)
    h->hlumi_vstrpu->Fill(trpu, _w);

  if (_jp_debug) {
    if (h == _histos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit)
        cout << *trgit << ", ";
      cout << "Current: " << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired and passed
  if (!fired or !_pass) return;

  if (_jp_debug) cout << Form("Subdirectory Eta_%1.1f-%1.1f/%s",h->etamin,h->etamax,h->trigname.c_str()) << endl;
  if (_jp_debug) cout << "Calculate and fill dijet mass" << endl << flush;

  if (h->ismcdir) h->hpthat->Fill(pthat, _w);
  if (h->ismcdir) h->hpthatnlo->Fill(pthat);

  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>_jp_recopt) { // First leading jet
    if (i1>=0 and _jetids[i1] and jtpt[i1]>_jp_recopt) { // Second leading jet
      //{ Calculate and fill dijet mass.
      _j1.SetPtEtaPhiE(jtpt[i0],jteta[i0],jtphi[i0],jte[i0]);
      _j2.SetPtEtaPhiE(jtpt[i1],jteta[i1],jtphi[i1],jte[i1]);
      double djmass = (_j1+_j2).M();
      double etamaxdj = max(fabs(jteta[i0]),fabs(jteta[i1]));
      bool goodjets = (jtpt[i0]>30. and jtpt[i1]>30.);
      // The eta sectors are filled according to max eta
      // NOTE: alpha and dphi cuts are currently completely missing!
      if (goodjets and etamaxdj >= h->etamin and etamaxdj < h->etamax) {
        assert(h->hdjmass); h->hdjmass->Fill(djmass, _w);
        assert(h->hdjmass0); h->hdjmass0->Fill(djmass, _w);
        assert(h->pdjmass_ptratio); h->pdjmass_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
        assert(h->pdjmass0_ptratio); h->pdjmass0_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
      }
      //} Dijet mass


      //{ Tag & probe hoods: Tag in barrel and fires trigger, probe in eta bin unbiased
      if (_jp_debug) cout << "Calculate and fill dijet balance" << endl << flush;

      double dphi = DPhi(jtphi[i0], jtphi[i1]);

      if (dphi > 2.7) { // Back-to-back condition
        double pt3 = ((i2>=0 and jtpt[i2]>_jp_recopt) ? jtpt[i2] : 0.);
        double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
        double alpha = pt3/ptave;

        for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
          int itag = jt3leads[itag_lead];
          int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
          double etatag = jteta[itag];
          double etaprobe = jteta[iprobe];

          if (fabs(etatag) < 1.3 and fabs(etaprobe) >= h->etamin and fabs(etaprobe) < h->etamax) { // Eta sel.
            double pttag = jtpt[itag];
            double ptprobe = jtpt[iprobe];
            double phiprobe = jtphi[iprobe];

            double alphatp = pt3/pttag;
            //{ Dijet balance histograms
            if (_jp_do3dHistos) {
              double asymm = (ptprobe - pttag)/(2*ptave);
              double asymmtp = (ptprobe - pttag)/(2*pttag);
              double mpf = met1*cos(DPhi(metphi1,jtphi[itag]))/2;
              double mpftp = mpf/pttag;
              mpf /= ptave;

              // Asymmetry and mpf
              assert(h->hdjasymm);   h->hdjasymm->Fill(ptave, alpha, asymm, _w);
              assert(h->hdjasymmtp); h->hdjasymmtp->Fill(pttag, alphatp, asymmtp, _w);
              assert(h->hdjmpf);     h->hdjmpf->Fill(ptave, alpha, mpf, _w);
              assert(h->hdjmpftp);   h->hdjmpftp->Fill(pttag, alphatp, mpftp, _w);
            }
            //} // Dijet balance

            if (alphatp < 0.3) {
              //{ Composition vs pt tag pt
              // Fractions vs pt: we do pt selection later in HistosCombine
              assert(h->pncandtp);    h->pncandtp->Fill(pttag, jtn[iprobe], _w);
              assert(h->pnchtp);      h->pnchtp->Fill(pttag, jtnch[iprobe], _w);
              assert(h->pnnetp);      h->pnnetp->Fill(pttag, jtnne[iprobe]-jtnhe[iprobe], _w);
              assert(h->pnnhtp);      h->pnnhtp->Fill(pttag, jtnnh[iprobe]-jtnhh[iprobe], _w);
              assert(h->pncetp);      h->pncetp->Fill(pttag, jtnce[iprobe], _w);
              assert(h->pnmutp);      h->pnmutp->Fill(pttag, jtnmu[iprobe], _w);

              assert(h->pchftp);      h->pchftp->Fill(pttag, jtchf[iprobe], _w);
              assert(h->pneftp);      h->pneftp->Fill(pttag, (jtnef[iprobe]-jthef[iprobe]), _w);
              assert(h->pnhftp);      h->pnhftp->Fill(pttag, (jtnhf[iprobe]-jthhf[iprobe]), _w);
              assert(h->pceftp);      h->pceftp->Fill(pttag, jtcef[iprobe], _w);
              assert(h->pmuftp);      h->pmuftp->Fill(pttag, jtmuf[iprobe], _w);
              assert(h->phhftp);      h->phhftp->Fill(pttag, jthhf[iprobe], _w);
              assert(h->pheftp);      h->pheftp->Fill(pttag, jthef[iprobe], _w);
              assert(h->pbetatp);     h->pbetatp->Fill(pttag, jtbeta[iprobe], _w);
              assert(h->pbetastartp); h->pbetastartp->Fill(pttag, jtbetastar[iprobe], _w);
              assert(h->pbetaprimetp); h->pbetaprimetp->Fill(pttag, jtbetaprime[iprobe], _w);

              assert(h->ppt_probepertag); h->ppt_probepertag->Fill(pttag,ptprobe/pttag,_w);

              double metstuff = met1 * cos(DPhi(metphi1, phiprobe));
              assert(h->pmpfz); h->pmpfz->Fill(ptave, 1 + metstuff / ptave, _w);
              if (ptave >= h->ptmin and ptave < h->ptmax) { // Ave fires trigger
                assert(h->hmpfz); h->hmpfz->Fill(1 + metstuff / ptave, _w);
              }
              assert(h->pmpfy); h->pmpfy->Fill(ptprobe, 1 + metstuff / ptprobe, _w);
              if (ptprobe >= h->ptmin and ptprobe < h->ptmax) { // Probe fires trigger
                assert(h->hmpfy); h->hmpfy->Fill(1 + metstuff / ptprobe, _w);
              }
              assert(h->pmpfx); h->pmpfx->Fill(pttag, 1 + metstuff / pttag, _w);
              if (pttag >= h->ptmin and pttag < h->ptmax) { // Tag fires trigger
                assert(h->hmpfx); h->hmpfx->Fill(1 + metstuff / pttag, _w);

                // The distributions of event counts per fraction
                assert(h->hncandtp);    h->hncandtp->Fill(jtn[iprobe], _w);
                assert(h->hnchtp);      h->hnchtp->Fill(jtnch[iprobe], _w);
                assert(h->hnnetp);      h->hnnetp->Fill(jtnne[iprobe]-jtnhe[iprobe], _w);
                assert(h->hnnhtp);      h->hnnhtp->Fill(jtnnh[iprobe]-jtnhh[iprobe], _w);
                assert(h->hncetp);      h->hncetp->Fill(jtnce[iprobe], _w);
                assert(h->hnmutp);      h->hnmutp->Fill(jtnmu[iprobe], _w);
                assert(h->hnhhtp);      h->hnhhtp->Fill(jtnhh[iprobe], _w);
                assert(h->hnhetp);      h->hnhetp->Fill(jtnhe[iprobe], _w);
                assert(h->hchftp);      h->hchftp->Fill(jtchf[iprobe], _w);
                assert(h->hneftp);      h->hneftp->Fill((jtnef[iprobe]-jthef[iprobe]), _w);
                assert(h->hnhftp);      h->hnhftp->Fill((jtnhf[iprobe]-jthhf[iprobe]), _w);
                assert(h->hceftp);      h->hceftp->Fill(jtcef[iprobe], _w);
                assert(h->hmuftp);      h->hmuftp->Fill(jtmuf[iprobe], _w);
                assert(h->hhhftp);      h->hhhftp->Fill(jthhf[iprobe], _w);
                assert(h->hheftp);      h->hheftp->Fill(jthef[iprobe], _w);
                assert(h->hbetatp);     h->hbetatp->Fill(jtbeta[iprobe], _w);
                assert(h->hbetastartp); h->hbetastartp->Fill(jtbetastar[iprobe], _w);
                assert(h->hbetaprimetp); h->hbetaprimetp->Fill(jtbetaprime[iprobe], _w);

                // Fractions vs number of primary vertices
                assert(h->pncandtp_vsnpv);    h->pncandtp_vsnpv->Fill(npvgood, jtn[iprobe], _w);
                assert(h->pnchtp_vsnpv);      h->pnchtp_vsnpv->Fill(npvgood, jtnch[iprobe], _w);
                assert(h->pnnetp_vsnpv);      h->pnnetp_vsnpv->Fill(npvgood, jtnne[iprobe]-jtnhe[iprobe], _w);
                assert(h->pnnhtp_vsnpv);      h->pnnhtp_vsnpv->Fill(npvgood, jtnnh[iprobe]-jtnhh[iprobe], _w);
                assert(h->pncetp_vsnpv);      h->pncetp_vsnpv->Fill(npvgood, jtnce[iprobe], _w);
                assert(h->pnmutp_vsnpv);      h->pnmutp_vsnpv->Fill(npvgood, jtnmu[iprobe], _w);
                assert(h->pnhhtp_vsnpv);      h->pnhhtp_vsnpv->Fill(npvgood, jtnhh[iprobe], _w);
                assert(h->pnhetp_vsnpv);      h->pnhetp_vsnpv->Fill(npvgood, jtnhe[iprobe], _w);
                assert(h->pchftp_vsnpv);      h->pchftp_vsnpv->Fill(npvgood, jtchf[iprobe], _w);
                assert(h->pneftp_vsnpv);      h->pneftp_vsnpv->Fill(npvgood, (jtnef[iprobe]-jthef[iprobe]), _w);
                assert(h->pnhftp_vsnpv);      h->pnhftp_vsnpv->Fill(npvgood, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                assert(h->pceftp_vsnpv);      h->pceftp_vsnpv->Fill(npvgood, jtcef[iprobe], _w);
                assert(h->pmuftp_vsnpv);      h->pmuftp_vsnpv->Fill(npvgood, jtmuf[iprobe], _w);
                assert(h->phhftp_vsnpv);      h->phhftp_vsnpv->Fill(npvgood, jthhf[iprobe], _w);
                assert(h->pheftp_vsnpv);      h->pheftp_vsnpv->Fill(npvgood, jthef[iprobe], _w);
                assert(h->pbetatp_vsnpv);     h->pbetatp_vsnpv->Fill(npvgood, jtbeta[iprobe], _w);
                assert(h->pbetastartp_vsnpv); h->pbetastartp_vsnpv->Fill(npvgood, jtbetastar[iprobe], _w);
                assert(h->pbetaprimetp_vsnpv); h->pbetaprimetp_vsnpv->Fill(npvgood, jtbetaprime[iprobe], _w);

                // Fractions vs true pileup
                assert(h->pchftp_vstrpu);      h->pchftp_vstrpu->Fill(trpu, jtchf[iprobe], _w);
                assert(h->pneftp_vstrpu);      h->pneftp_vstrpu->Fill(trpu, (jtnef[iprobe]-jthef[iprobe]), _w);
                assert(h->pnhftp_vstrpu);      h->pnhftp_vstrpu->Fill(trpu, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                assert(h->pceftp_vstrpu);      h->pceftp_vstrpu->Fill(trpu, jtcef[iprobe], _w);
                assert(h->pmuftp_vstrpu);      h->pmuftp_vstrpu->Fill(trpu, jtmuf[iprobe], _w);
                assert(h->phhftp_vstrpu);      h->phhftp_vstrpu->Fill(trpu, jthhf[iprobe], _w);
                assert(h->pheftp_vstrpu);      h->pheftp_vstrpu->Fill(trpu, jthef[iprobe], _w);
                assert(h->pbetatp_vstrpu);     h->pbetatp_vstrpu->Fill(trpu, jtbeta[iprobe], _w);
                assert(h->pbetastartp_vstrpu); h->pbetastartp_vstrpu->Fill(trpu, jtbetastar[iprobe], _w);
                assert(h->pbetaprimetp_vstrpu); h->pbetaprimetp_vstrpu->Fill(trpu, jtbetaprime[iprobe], _w);

                if (_jp_doPhiHistos) {
                  if (etaprobe>0) {
                    assert(h->pchfpostp_vsphi);      h->pchfpostp_vsphi->Fill(phiprobe, jtchf[iprobe], _w);
                    assert(h->pnefpostp_vsphi);      h->pnefpostp_vsphi->Fill(phiprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
                    assert(h->pnhfpostp_vsphi);      h->pnhfpostp_vsphi->Fill(phiprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                    assert(h->pcefpostp_vsphi);      h->pcefpostp_vsphi->Fill(phiprobe, jtcef[iprobe], _w);
                    assert(h->pmufpostp_vsphi);      h->pmufpostp_vsphi->Fill(phiprobe, jtmuf[iprobe], _w);
                    assert(h->phhfpostp_vsphi);      h->phhfpostp_vsphi->Fill(phiprobe, jthhf[iprobe], _w);
                    assert(h->phefpostp_vsphi);      h->phefpostp_vsphi->Fill(phiprobe, jthef[iprobe], _w);
                    assert(h->pbetapostp_vsphi);     h->pbetapostp_vsphi->Fill(phiprobe, jtbeta[iprobe], _w);
                    assert(h->pbetastarpostp_vsphi); h->pbetastarpostp_vsphi->Fill(phiprobe, jtbetastar[iprobe], _w);
                    assert(h->pbetaprimepostp_vsphi); h->pbetaprimepostp_vsphi->Fill(phiprobe, jtbetaprime[iprobe], _w);
                  } else {
                    assert(h->pchfnegtp_vsphi);      h->pchfnegtp_vsphi->Fill(phiprobe, jtchf[iprobe], _w);
                    assert(h->pnefnegtp_vsphi);      h->pnefnegtp_vsphi->Fill(phiprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
                    assert(h->pnhfnegtp_vsphi);      h->pnhfnegtp_vsphi->Fill(phiprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                    assert(h->pcefnegtp_vsphi);      h->pcefnegtp_vsphi->Fill(phiprobe, jtcef[iprobe], _w);
                    assert(h->pmufnegtp_vsphi);      h->pmufnegtp_vsphi->Fill(phiprobe, jtmuf[iprobe], _w);
                    assert(h->phhfnegtp_vsphi);      h->phhfnegtp_vsphi->Fill(phiprobe, jthhf[iprobe], _w);
                    assert(h->phefnegtp_vsphi);      h->phefnegtp_vsphi->Fill(phiprobe, jthef[iprobe], _w);
                    assert(h->pbetanegtp_vsphi);     h->pbetanegtp_vsphi->Fill(phiprobe, jtbeta[iprobe], _w);
                    assert(h->pbetastarnegtp_vsphi); h->pbetastarnegtp_vsphi->Fill(phiprobe, jtbetastar[iprobe], _w);
                    assert(h->pbetaprimenegtp_vsphi); h->pbetaprimenegtp_vsphi->Fill(phiprobe, jtbetaprime[iprobe], _w);
                  }
                }
              } // Tag fires trigger
            //} Composition vs pt tag pt
            } // dijet system
          } // etatag < 1.3
        } // for itag (two leading jets)
      } // Back-to-back
    //} Tag and probe
    } // Second leading jet

    if (jtpt[i0]>=h->ptmin and jtpt[i0]<h->ptmax and fabs(jteta[i0]) < 1.3) { // Jet quality stats
      bool has2 = (i1>=0 and jtpt[i1] > _jp_recopt and fabs(jteta[i1])>=h->etamin and fabs(jteta[i1])<h->etamax);
      bool has3 = (i2>=0 and jtpt[i2] > _jp_recopt and fabs(jteta[i2])>=h->etamin and fabs(jteta[i2])<h->etamax and jtpt[i1] > 0.70 * jtpt[i0]);
      bool has32 = (has3 and fabs(jteta[i1]) < 1.3);

      h->hr21->Fill(has2  ? jtpt[i1] / jtpt[i0] : 0.);
      h->hr31->Fill(has3  ? jtpt[i2] / jtpt[i0] : 0.);
      h->hr32->Fill(has32 ? jtpt[i2] / jtpt[i1] : 0.);
      if (has2)  h->pr21->Fill(npvgood, has2 ? jtpt[i1] / jtpt[i0] : 0.);
      if (has3)  h->pr31->Fill(npvgood, has3 ? jtpt[i2] / jtpt[i0] : 0.);
      if (has32) h->pr32->Fill(npvgood, has3 ? jtpt[i2] / jtpt[i1] : 0.);
      h->px21->Fill(npvgood, has2  ? 1 : 0);
      h->px31->Fill(npvgood, has3  ? 1 : 0);
      h->px32->Fill(npvgood, has32 ? 1 : 0);
    } // Jet quality stats
  } // First leading jet

  // retrieve event-wide variables
  double dphi = (i1>=0 ? DPhi(jtphi[i0], jtphi[i1]) : 0.);
  double dpt = (i1>=0 ? fabs(jtpt[i0]-jtpt[i1])/(jtpt[i0]+jtpt[i1]) : 0.999);

  if (_jp_debug) cout << "Entering jet loop" << endl << flush;
  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    if (_jp_debug) cout << "Loop over jet " << jetidx << "/" << njt << endl << flush;

    // adapt variable names from different trees
    double pt = jtpt[jetidx];
    double eta = jteta[jetidx];
    double energy = jte[jetidx];
    double mass = sqrt(fabs(pow(energy,2) - pow(pt*cosh(eta),2)));
    double y = jty[jetidx];
    double phi = jtphi[jetidx];
    double jec = jtjesnew[jetidx];
    bool id = _jetids[jetidx];

    double jec2 = jtjesnew[jetidx]/jtjes[jetidx];
    bool etarange = fabs(eta) >= h->etamin and fabs(eta) < h->etamax;
    double ptgen = jtgenpt[jetidx];

    // Check effect of ID cuts
    if (etarange) { // Jet in eta range
      if (_jp_debug) {
        cout << "..." << h->trigname << " | " << " index " << jetidx << "/" << njt
          << " jet pt: " << pt << " eta : " << eta << " id: " << id << " jec: " << jec << endl;
        cout << "...evt id: " << _pass_qcdmet << " weight: " << _w << " met: " << met << " metsumet: " << metsumet << endl;
      }

      assert(h->hpt_noid); h->hpt_noid->Fill(pt, _w);
      assert(h->hpt_nojetid); if (_pass_qcdmet) h->hpt_nojetid->Fill(pt, _w);
      assert(h->hpt_noevtid); if (id) h->hpt_noevtid->Fill(pt, _w);
      // Same versus generator pT as MC extra
      // to decouple efficiency from JEC and JER
      if (h->ismcdir) {
        h->hpt_noid_g->Fill(ptgen, _w);
        if (_pass_qcdmet) h->hpt_nojetid_g->Fill(ptgen, _w);
        if (id)    h->hpt_noevtid_g->Fill(ptgen, _w);
      }
    } // Jet in eta range

    if (pt>_jp_recopt) { // pt visible
      // Flags for studying gen eta vs reco eta effects
      bool mcgendr = jtgenr[jetidx] < 0.25;
      bool mcgenetarange = fabs(jtgeneta[jetidx]) >= h->etamin && fabs(jtgeneta[jetidx]) < h->etamax;

      if (h->ismcdir and mcgendr and mcgenetarange) // GenJets matched to any reco jets in any events
        h->hpt_gg0->Fill(ptgen, _w);

      if (_pass_qcdmet and id) { // id OK
        if (h->ismcdir and mcgendr and mcgenetarange) // GenJets matched to good reco jets in good events
          h->hpt_gg->Fill(ptgen, _w);

        if (etarange) { // Correct jet eta range
          if (_jp_debug) cout << "..unfolding" << endl << flush;
          if (h->ismcdir) {
            // unfolding studies (Mikael)
            h->my->Fill(pt, _w);
            h->myf->Fill(pt, _w);
            h->myuw->Fill(pt);
            h->myfuw->Fill(pt);

            if (mcgendr) {
              if (ptgen>_jp_recopt) {
                double etagen = jtgeneta[jetidx];
                if (fabs(etagen) >= h->etamin and fabs(etagen) < h->etamax) {
                  h->mT->Fill(ptgen, jtpt[jetidx], _w);
                  h->mTf->Fill(ptgen, jtpt[jetidx], _w);
                  h->mTuw->Fill(ptgen, jtpt[jetidx]);
                  h->mTfuw->Fill(ptgen, jtpt[jetidx]);
                }
              }
            } // Delta r
          }

          if (_jp_debug) cout << "..jec uncertainty" << endl << flush;
          // Get JEC uncertainty
          double unc = 0.01; // default for MC
          if (_jp_isdt and _jecUnc) {
            _jecUnc->setJetEta(eta);
            _jecUnc->setJetPt(pt);
            unc = _jecUnc->getUncertainty(true);
            //_jecUnc2->Rjet(pt, unc); // use Fall10 absolute scale uncertainty
          }

          // calculate and/or retrieve efficiencies
          double ideff = 1.;
          double vtxeff = 1.;
          double dqmeff = 1.;
          double trigeff = 1.;
          double eff = ideff * vtxeff * dqmeff * trigeff;

          if (_jp_debug) cout << "..raw spectrum" << endl << flush;

          // REMOVED: "For trigger efficiency"

          // new histograms for quark/gluon study (Ozlem)
          double probg = 1.-qgl[jetidx]; // First approximation
          if (_jp_doqglfile) { // If we loaded a previous file to _h3probg, use this for better probg
            assert(_h3probg);
            probg = _h3probg->GetBinContent(_h3probg->FindBin(eta,pt,qgl[jetidx]));
          }
          if (probg>=0 and probg<=1) {
            assert(h->hgpt);  h->hgpt->Fill(pt,_w*probg);
            assert(h->hgpt0); h->hgpt0->Fill(pt, _w*probg);

            assert(h->hqgl);  h->hqgl->Fill(qgl[jetidx], _w);
            assert(h->hqgl2); h->hqgl2->Fill(pt, qgl[jetidx], _w);
            if (_jp_ismc) {
              assert(h->hqgl_g);
              assert(h->hqgl_q);
              bool isgluon = (fabs(partonflavor[jetidx]-21)<0.5);
              bool isquark = (fabs(partonflavor[jetidx])<7);
              assert(isgluon || isquark);

              // For data templates from scaling Pythia (wq & wg), see instructions at
              // https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#Systematics
              double x = qgl[jetidx];
              if (isgluon) {
                h->hqgl_g->Fill(x, _w);
                h->hqgl2_g->Fill(pt, x, _w);
                double wg = 1;//-55.7067*pow(x,7) + 113.218*pow(x,6) -21.1421*pow(x,5) -99.927*pow(x,4) + 92.8668*pow(x,3) -34.3663*x*x + 6.27*x + 0.612992;
                assert(wg>0);
                h->hqgl_dg->Fill(x, _w*wg);
                h->hqgl2_dg->Fill(pt, x, _w*wg);
              } else if (isquark) {
                h->hqgl_q->Fill(x, _w);
                h->hqgl2_q->Fill(pt, x, _w);
                double wq = 1;// -0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581;
                assert(wq>0);
                h->hqgl_dq->Fill(x, _w*wq);
                h->hqgl2_dq->Fill(pt, x, _w*wq);
              } else {
                PrintInfo("Quark/Gluon status missing from partonflavor");
              }
            }
          } // probg quark/gluon

          // raw spectrum
          assert(h->hpt); h->hpt->Fill(pt,_w);
          assert(h->hpt_tmp); h->hpt_tmp->Fill(pt); // Event statistics
          assert(h->hpt_pre);
          if (_jp_isdt) h->hpt_pre->Fill(pt, _w*_prescales[h->trigname][run]);
          if (_jp_ismc) h->hpt_pre->Fill(pt, _w0*_wt["mc"]);
          assert(h->hpt0); h->hpt0->Fill(pt, _w);
          // REMOVED: "h->hpt_plus_38x->Fill(pt, _w);" etc.
          // Do proper event statistics
          assert(h->hpttmp); assert(h->hptevt);
          if (h->hpttmp->GetBinContent(h->hpttmp->FindBin(pt))==0)
            h->hptevt->Fill(pt, _w);
          h->hpttmp->Fill(pt);

          // leading and non-leading jets
          if (jetidx==i0) {
            assert(h->hpt1); h->hpt1->Fill(pt, _w);
          } else if (jetidx==i1) {
            assert(h->hpt2); h->hpt2->Fill(pt,_w);
          } else if (jetidx==i2) {
            assert(h->hpt3); h->hpt3->Fill(pt,_w);
          }

          if (_jp_debug) cout << "..basic properties" << endl << flush;

          // basic properties
          assert(h->ppt); h->ppt->Fill(pt, pt, _w);
          assert(h->pmass); h->pmass->Fill(pt, mass/energy, _w);
          assert(h->pjec); h->pjec->Fill(pt, jec, _w);
          assert(h->pjec2); h->pjec2->Fill(pt, jec2, _w);
          assert(h->punc); h->punc->Fill(pt, unc, _w);
          // JEC monitoring
          assert(h->pjec_l1); h->pjec_l1->Fill(pt, jtjes_l1[jetidx], _w);
          assert(h->pjec_l2l3); h->pjec_l2l3->Fill(pt, jtjes_l2l3[jetidx], _w);
          assert(h->pjec_res); h->pjec_res->Fill(pt, jtjes_res[jetidx], _w);

          // efficiencies
          assert(h->peff); h->peff->Fill(pt, eff, _w);
          assert(h->pideff); h->pideff->Fill(pt, ideff, _w);
          assert(h->pvtxeff); h->pvtxeff->Fill(pt, vtxeff, _w);
          assert(h->pdqmeff); h->pdqmeff->Fill(pt, dqmeff, _w);

          if (_jp_debug) cout << "..control plots of components" << endl << flush;

          // Composition stuff without T&P (according to triggers)
          assert(h->pncand); h->pncand->Fill(pt, jtn[jetidx], _w);
          assert(h->pnch); h->pnch->Fill(pt, jtnch[jetidx], _w);
          assert(h->pnne); h->pnne->Fill(pt, jtnne[jetidx]-jtnhe[jetidx], _w);
          assert(h->pnnh); h->pnnh->Fill(pt, jtnnh[jetidx]-jtnhh[jetidx], _w);
          assert(h->pnce); h->pnce->Fill(pt, jtnce[jetidx], _w);
          assert(h->pnmu); h->pnmu->Fill(pt, jtnmu[jetidx], _w);
          assert(h->pnhh); h->pnhh->Fill(pt, jtnhh[jetidx], _w);
          assert(h->pnhe); h->pnhe->Fill(pt, jtnhe[jetidx], _w);
          //
          assert(h->pchf); h->pchf->Fill(pt, jtchf[jetidx], _w);
          assert(h->pnef); h->pnef->Fill(pt, (jtnef[jetidx]-jthef[jetidx]), _w);
          assert(h->pnhf); h->pnhf->Fill(pt, (jtnhf[jetidx]-jthhf[jetidx]), _w);
          assert(h->pcef); h->pcef->Fill(pt, jtcef[jetidx], _w);
          assert(h->pmuf); h->pmuf->Fill(pt, jtmuf[jetidx], _w);
          assert(h->phhf); h->phhf->Fill(pt, jthhf[jetidx], _w);
          assert(h->phhf); h->phef->Fill(pt, jthef[jetidx], _w);
          assert(h->pbeta); h->pbeta->Fill(pt, jtbeta[jetidx], _w);
          assert(h->pbetastar); h->pbetastar->Fill(pt, jtbetastar[jetidx], _w);
          assert(h->pbetaprime); h->pbetaprime->Fill(pt, jtbetaprime[jetidx], _w);

          // control plots for topology (JEC)
          h->pa->Fill(pt, jta[jetidx], _w);
          h->ptrpu->Fill(pt, trpu, _w);
          h->prho->Fill(pt, rho, _w);
          h->pnpv->Fill(pt, npvgood, _w);
          h->pnpvall->Fill(pt, npv, _w);
          if (pt >= h->ptmin && pt < h->ptmax) { // Trigger pt range
            h->htrpu2->Fill(trpu, _w);
            //
            h->pnpvvsrho->Fill(rho, npvgood, _w);
            h->prhovsnpv->Fill(npvgood, rho, _w);
            h->prhovsnpvall->Fill(npv, rho, _w);
            h->h2rhovsnpv->Fill(npvgood, rho, _w);
            //
            h->prhovstrpu->Fill(trpu, rho, _w);
            h->pnpvvstrpu->Fill(trpu, npvgood, _w);
            h->pnpvallvstrpu->Fill(trpu, npv, _w);
            h->pitpuvstrpu->Fill(trpu, itpu, _w);
            h->hjet_vstrpu->Fill(trpu, _w);

            if (_jp_debug) cout << "..control plots for topology" << endl << flush;

            h->htrpu->Fill(trpu, _w);
            if (h->ismcdir) {
              h->hitpu->Fill(itpu, _w);
              h->hootpuearly->Fill(ootpuearly, _w);
              h->hootpulate->Fill(ootpulate, _w);
              h->h2itvsoot->Fill(itpu, ootpulate, _w);
            }

            h->hnpvgood->Fill(npvgood, _w);
            h->hrho->Fill(rho, _w);
            h->hselpt->Fill(pt, _w);
            h->hmass->Fill(mass/energy, _w);
            h->hy->Fill(y, _w);
            h->hy2->Fill(y, _w);
            h->heta->Fill(eta, _w);
            h->heta2->Fill(eta, _w);
            h->hphi->Fill(phi, _w);
            h->hdphi->Fill(dphi, _w);
            h->hdpt->Fill(dpt, _w);
            h->hjet->Fill(pt / metsumet, _w);
            h->hmet->Fill(met / metsumet, _w);
            h->hmetphi->Fill(DPhi(metphi, phi), _w);
            // control plots for vertex
            h->hpvndof->Fill(pvndof);
            h->hpvx->Fill(pvx-bsx);
            h->hpvy->Fill(pvy-bsy);
            h->hpvz->Fill(pvz-0.);
            h->hpvr->Fill(tools::oplus(pvx-bsx, pvy-bsy));
            h->hpvrho->Fill(pvrho-tools::oplus(bsx, bsy));
            // closure plots for JEC
            h->hmpf->Fill(1 + met * cos(DPhi(metphi, phi)) / pt, _w);
            h->hmpf1->Fill(1 + met1 * cos(DPhi(metphi1, phi)) / pt, _w);
            h->hmpf2->Fill(1 + met2 * cos(DPhi(metphi2, phi)) / pt, _w);
            // Component fractions
            h->hncand->Fill(jtn[jetidx], _w);
            h->hnch->Fill(jtnch[jetidx], _w);
            h->hnne->Fill(jtnne[jetidx]-jtnhe[jetidx], _w);
            h->hnnh->Fill(jtnnh[jetidx]-jtnhh[jetidx], _w);
            h->hnce->Fill(jtnce[jetidx], _w);
            h->hnmu->Fill(jtnmu[jetidx], _w);
            //
            h->hchf->Fill(jtchf[jetidx], _w);
            h->hnef->Fill((jtnef[jetidx]-jthef[jetidx]), _w);
            h->hnhf->Fill((jtnhf[jetidx]-jthhf[jetidx]), _w);
            h->hcef->Fill(jtcef[jetidx], _w);
            h->hmuf->Fill(jtmuf[jetidx], _w);
            h->hhhf->Fill(jthhf[jetidx], _w);
            h->hhef->Fill(jthef[jetidx], _w);
            h->hbeta->Fill(jtbeta[jetidx], _w);
            h->hbetastar->Fill(jtbetastar[jetidx], _w);
            h->hbetaprime->Fill(jtbetaprime[jetidx], _w);

            h->hyeta->Fill(TMath::Sign(y-eta,y), _w);
            h->hyeta2->Fill(y-eta, _w);
            h->hbetabetastar->Fill(jtbeta[jetidx], jtbetastar[jetidx], _w);
            h->hbetabetaprime->Fill(jtbeta[jetidx], jtbetaprime[jetidx], _w);
            h->hetaphi->Fill(eta, phi, _w);
          } // within trigger pT range

          // closure plots for JEC
          h->pdpt->Fill(pt, dpt, _w);
          h->pmpf->Fill(pt, 1 + met * cos(DPhi(metphi, phi)) / pt, _w);
          h->pmpf1->Fill(pt, 1 + met1 * cos(DPhi(metphi1, phi)) / pt, _w);
          h->pmpf2->Fill(pt, 1 + met2 * cos(DPhi(metphi2, phi)) / pt, _w);

          if (h->ismcdir and mcgendr) { // MC extras
            if (_jp_debug)
              cout << "genmatch " << jetidx << " ptg="<<ptgen << " yg="<<jtgeny[jetidx] << " yr="<< y << endl;

            double r = (ptgen ? pt/ptgen : 0);
            //double resp = (jtjesnew[jetidx] ? r / jtjesnew[jetidx] : 0);
            double dy = (r ? TMath::Sign(jty[jetidx]-jtgeny[jetidx], jtgeny[jetidx]) : 0.);
            h->hpt_r->Fill(pt, _w);
            h->hpt_g->Fill(ptgen, _w);

            h->hpt_gtw->Fill(ptgen, _w);
            // Ozlem: (gluon vs quark)
            if (partonflavor[jetidx]==21) h->hgpt_g->Fill(ptgen, _w);
            else h->hqpt_g->Fill(ptgen, _w);

            h->ppt_r->Fill(pt, pt, _w);
            h->ppt_g->Fill(ptgen, ptgen, _w);

            if (r) {
              // Response closure vs NPV
              h->p2rvsnpv->Fill(ptgen, npvgood, r, _w);

              // Response closure
              h->h2r_r->Fill(pt, r, _w);
              h->h2r_g->Fill(ptgen, r, _w);
              h->p2r_r->Fill(pt, r, _w);
              h->p2r_g->Fill(ptgen, r, _w);
              h->p2r_ruw->Fill(pt, r); // unweighted!
              h->p2r_guw->Fill(ptgen, r); // unweighted!

              // Rapidity closure
              h->h2dy_r->Fill(pt, dy, _w);
              h->h2dy_g->Fill(ptgen, dy, _w);
              h->p2dy_r->Fill(pt, dy, _w);
              h->p2dy_g->Fill(ptgen, dy, _w);
              h->p2dy_ruw->Fill(pt, dy); // unweighted
              h->p2dy_guw->Fill(ptgen, dy); // unweighted
              h->pdy_r->Fill(pt, fabs(y), dy, _w);
              h->pdy_g->Fill(ptgen, fabs(y), dy, _w);
            }
          } // is MC
        } // etabin
      } // if id
    } // pt visible
  } // for jetidx

  // Event statistics
  for (int xidx = 1; xidx != h->hpt_tmp->GetNbinsX()+1; ++xidx) {
    if (h->hpt_tmp->GetBinContent(xidx)!=0) {
      double pt = h->hpt_tmp->GetBinCenter(xidx);
      int njet = h->hpt_tmp->GetBinContent(xidx);
      h->hpt_evtcount->Fill(pt);
      h->hpt_evt->Fill(pt, _w);
      h->hpt_jet->Fill(pt, _w*njet);
    }
  } // for xidx

  if (_jp_ismc) {
    if (_jp_debug) cout << "Truth loop:" << endl;
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) { // Unbiased gen spectrum (for each trigger)
      double etagen = gen_jteta[gjetidx];
      double ptgen = gen_jtpt[gjetidx];

      if (fabs(etagen) >= h->etamin and fabs(etagen) < h->etamax) {
        if (_jp_debug)
          cout << "genjet " << gjetidx << "/" << gen_njt << " ptg="<<ptgen << " etag="<<etagen << endl;

        h->hpt_g0tw->Fill(ptgen, _w);
        // Ozlem: (gluon vs quark)
        if (gen_partonflavor[gjetidx]==21) h->hgpt_g0tw->Fill(ptgen, _w);
        else h->hqpt_g0tw->Fill(ptgen, _w);

        // unfolding studies (Mikael)
        if (h->ismcdir) { // Only the 'mc' trigger
          if (gen_jtpt[gjetidx]>_jp_recopt) {
            h->mx->Fill(ptgen, _w);
            h->mxf->Fill(ptgen, _w);
            h->mxuw->Fill(ptgen);
            h->mxfuw->Fill(ptgen);
          } // recopt limit

          h->hpt_g0->Fill(gen_jtpt[gjetidx], _w);
          // Ozlem: (gluon vs quark)
          if (gen_partonflavor[gjetidx]==21) h->hgpt_g0->Fill(gen_jtpt[gjetidx], _w);
          else h->hqpt_g0->Fill(gen_jtpt[gjetidx], _w);

          assert(h->hpt_g0_tmp); h->hpt_g0_tmp->Fill(gen_jtpt[gjetidx]);
        } // mcdir (a subset of _jp_ismc)
      } // gen jet eta
    } // genjet loop
  } // MC
} // FillSingleBasic


// Write and delete histograms
void HistosFill::WriteBasic()
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteBasic():");
  PrintMemInfo();

  for (auto histit = _histos.begin(); histit != _histos.end(); ++histit) {
    for (auto histidx = 0u; histidx != histit->second.size(); ++histidx) {
      // Luminosity information
      HistosBasic *h = histit->second[histidx];
      for (int j = 0; j != h->hlumi->GetNbinsX()+1; ++j) {
        h->hlumi->SetBinContent(j, _jp_isdt ? h->lumsum : 1. );
        h->hlumi2->SetBinContent(j, _jp_isdt ? h->lumsum2 : 1. );
      }
      delete h; // histit->second[histidx];
    } // for histidx
  } // for histit

  cout << "\nOutput stored in " << _outfile->GetName() << endl;
  _outfile->Close();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("writeBasic() finished:");
  PrintMemInfo();
} // writeBasic


// Initialize eta histograms for trigger bins
void HistosFill::InitEta(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitEta(%s):",name.c_str()));
  PrintMemInfo();

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = _outfile;
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  bool enteroutdir = f->cd(name.c_str());
  assert(enteroutdir);
  TDirectory *topdir = f->GetDirectory(name.c_str()); assert(topdir);
  topdir->cd();

  // define triggers
  vector<string> triggers;
  // define efficient pT ranges for triggers for control plots
  map<string, pair<double, double> > pt;
  // define pT values for triggers
  map<string, double> pttrg;
  if (_jp_ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(_jp_recopt, _jp_emax);
    pttrg["mc"] = _jp_recopt;
  }
  if (_jp_isdt or _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_notrigs; ++itrg) {
      string trg = _jp_triggers[itrg];
      triggers.push_back(trg);
      double pt1 = _jp_trigranges[itrg][0];
      double pt2 = _jp_trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = _jp_trigthr[itrg];
      pttrg[trg] = pt0;
    }
  }

  assert(topdir);

  for (unsigned int j = 0; j != triggers.size(); ++j) {
    // subdirectory for trigger
    const char *trg = triggers[j].c_str();
    topdir->mkdir(trg);
    assert(topdir->cd(trg));
    TDirectory *dir = topdir->GetDirectory(trg);
    assert(dir);
    dir->cd();

    // Initialize and store
    assert(dir);
    HistosEta *h = new HistosEta(dir, trg,pttrg[trg],pt[trg].first, pt[trg].second, triggers[j]=="mc");
    _etahistos[name].push_back(h);
  } // for j

  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitEta(%s) finished:",name.c_str()));
  PrintMemInfo();
} // InitEta


// Loop over basic histogram containers to fill all
void HistosFill::FillEta(string name, Float_t* _pt, Float_t* _eta, Float_t* _phi)
{
  for (auto histidx = 0u; histidx != _etahistos[name].size(); ++histidx)
    FillSingleEta(_etahistos[name][histidx], _pt, _eta, _phi);
}


// Fill basic histograms after applying pt, y cuts
void HistosFill::FillSingleEta(HistosEta *h, Float_t* _pt, Float_t* _eta, Float_t* _phi)
{
  assert(h);

  _w = _w0 * _wt[h->trigname];
  if (_w <= 0) return;

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  if (_jp_debug) {
    if (h == _etahistos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit)
        cout << *trgit << ", ";
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired
  if (!fired or !_pass) return;

  // Tag & probe hoods
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>_jp_recopt and i1>=0 and _jetids[i1] and jtpt[i1]>_jp_recopt) { // Leading jets
    double dphi = DPhi(jtphi[i0], jtphi[i1]);
    if (dphi > 2.7) { // Back-to-back condition
      double pt3 = ((i2>=0 and jtpt[i2]>_jp_recopt) ? jtpt[i2] : 0.);
      double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
      double alpha = pt3/ptave;

      for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
        int itag = jt3leads[itag_lead];
        int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
        double etatag = jteta[itag];
        double etaprobe = jteta[iprobe];
        if (fabs(etatag) < 1.3) {
          double pttag = jtpt[itag];
          double ptprobe = _pt[iprobe];
          double alphatp = pt3/pttag;
          if (_jp_do3dHistos) {
            double asymm = (ptprobe - pttag)/(2*ptave);
//             double asymmtp = (ptprobe - pttag)/(2*pttag);
            double mpf = met1*cos(DPhi(metphi1,_phi[itag]))/(2*ptave);
//             double mpftp = met2*cos(DPhi(metphi2,_phi[itag]))/(2*pttag);
            for (auto alphaidx = 0u; alphaidx < h->alpharange.size(); ++alphaidx) {
              float alphasel = h->alpharange[alphaidx];
              if (alpha<alphasel) {
                // Val 10 = excluded, -10 = ok
                h->hdjasymm[alphaidx]  ->Fill(ptave,   etaprobe, asymm, _w);
                h->hdjmpf[alphaidx]    ->Fill(ptave,   etaprobe, mpf  , _w);
                h->hdjasymmtp[alphaidx]->Fill(pttag,   etaprobe, asymm, _w);
                h->hdjmpftp[alphaidx]  ->Fill(pttag,   etaprobe, mpf  , _w);
                h->hdjasymmpt[alphaidx]->Fill(ptprobe, etaprobe, asymm, _w);
                h->hdjmpfpt[alphaidx]  ->Fill(ptprobe, etaprobe, mpf  , _w);
              }
//               if (alphatp<alphasel) {
//                 h->hdjasymmtp[alphaidx]->Fill(pttag,   etaprobe, asymmtp, _w);
//                 h->hdjmpftp[alphaidx]  ->Fill(pttag,   etaprobe, mpftp  , _w);
//               }
            }
          }

          // for composition vs eta (Ozlem) || jetidx is the probe
          if (alphatp < 0.3 and pttag > h->ptmin and pttag < h->ptmax) { // Alpha and trigger
            assert(h->pchftp_vseta); h->pchftp_vseta->Fill(etaprobe, jtchf[iprobe], _w);
            assert(h->pneftp_vseta); h->pneftp_vseta->Fill(etaprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
            assert(h->pnhftp_vseta); h->pnhftp_vseta->Fill(etaprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
            assert(h->pceftp_vseta); h->pceftp_vseta->Fill(etaprobe, jtcef[iprobe], _w);
            assert(h->pmuftp_vseta); h->pmuftp_vseta->Fill(etaprobe, jtmuf[iprobe], _w);
            assert(h->phhftp_vseta); h->phhftp_vseta->Fill(etaprobe, jthhf[iprobe], _w);
            assert(h->pheftp_vseta); h->pheftp_vseta->Fill(etaprobe, jthef[iprobe], _w);
            assert(h->pbetatp_vseta); h->pbetatp_vseta->Fill(etaprobe, jtbeta[iprobe], _w);
            assert(h->pbetastartp_vseta); h->pbetastartp_vseta->Fill(etaprobe, jtbetastar[iprobe], _w);
            assert(h->pbetaprimetp_vseta); h->pbetaprimetp_vseta->Fill(etaprobe, jtbetaprime[iprobe], _w);
          } // select pt bin for profiles vseta
        } // etatag < 1.3
      } // tag & probe
    } // dphi > 2.7
  } // ids

  if (h->ismcdir and _pass_qcdmet) {
    for (int jetidx = 0; jetidx != njt; ++jetidx) {
      double pt = jtpt[jetidx];
      if (pt>_jp_recopt) { // pt visible
        bool mcgendr = jtgenr[jetidx] < 0.25;
        bool id = _jetids[jetidx];
        if (mcgendr and id) { // id OK
          double ptgen = jtgenpt[jetidx];
          if (ptgen>0 and jtjesnew[jetidx]!=0) {
            // CAUTION: Do we lay jtjesnew over the old correction?
            double resp = pt / (ptgen*jtjesnew[jetidx]);
            // Response closure vs NPV
            if (resp) {
              double eta = jteta[jetidx];
              h->p3rvsnpv->Fill(ptgen, eta, npvgood, resp, _w);
              h->p3rvsnpvW->Fill(ptgen, fabs(eta), npvgood, resp, _w);
            }
          }
        } // if id
      } // pt visible
    } // for jetidx
  } // if MC
} // FillSingleEta


// Write and delete histograms
void HistosFill::WriteEta()
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteEta():");
  PrintMemInfo();

  for (auto histit : _etahistos) {
    for (unsigned int histidx = 0; histidx != histit.second.size(); ++histidx) {
      HistosEta *h = histit.second[histidx];
      delete h;
    } // for histidx
  } // for histit

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteEta() finished:");
  PrintMemInfo();
} // WriteEta


// Initialize eta histograms for trigger bins
void HistosFill::InitMC(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("InitMC(%s):",name.c_str());
  PrintMemInfo();

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = _outfile;
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  bool entertopdir = f->cd(name.c_str());
  assert(entertopdir);
  TDirectory *topdir = f->GetDirectory(name.c_str()); assert(topdir);
  topdir->cd();

  // define triggers
  vector<string> triggers;
  // define efficient pT ranges for triggers for control plots
  map<string, pair<double, double> > pt;
  // define pT values for triggers
  map<string, double> pttrg;
  if (_jp_ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(_jp_recopt, _jp_emax);
    pttrg["mc"] = _jp_recopt;
  }
  if (_jp_isdt or _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_notrigs; ++itrg) {
      string trg = _jp_triggers[itrg];
      triggers.push_back(trg);
      double pt1 = _jp_trigranges[itrg][0];
      double pt2 = _jp_trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = _jp_trigthr[itrg];
      pttrg[trg] = pt0;
    }
  }

  assert(topdir);

  for (unsigned int j = 0; j != triggers.size(); ++j) {
    // subdirectory for trigger
    const char *trg = triggers[j].c_str();
    topdir->mkdir(trg);
    bool entertrgdir = topdir->cd(trg);
    assert(entertrgdir);
    TDirectory *dir = topdir->GetDirectory(trg);
    assert(dir);
    dir->cd();

    // Initialize and store
    assert(dir);
    HistosMC *h = new HistosMC(dir, trg);
    _mchistos[name].push_back(h);
  } // for j

  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitMC(%s) finished:",name.c_str()));
  PrintMemInfo();
} // InitMC


// Loop over basic histogram containers to fill all
void HistosFill::FillMC(string name,  Float_t* _recopt, Float_t* _genpt,
                              Float_t* _pt, Float_t* _eta,    Float_t* _phi)
{
  for (auto histidx = 0u; histidx != _mchistos[name].size(); ++histidx)
    FillSingleMC(_mchistos[name][histidx], _recopt, _genpt, _pt, _eta, _phi);
}


// Fill basic histograms after applying pt, y cuts
void HistosFill::FillSingleMC(HistosMC *h,  Float_t* _recopt,  Float_t* _genpt,
                             Float_t* _pt, Float_t* _eta,     Float_t* _phi)
{
  assert(h);

  _w = _w0 * _wt[h->trigname];
  if (_w <= 0) return;

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  if (_jp_debug) {
    if (h == _mchistos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit)
        cout << *trgit << ", ";
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired
  if (!fired or !_pass) return;

  // Tag & probe hoods
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>_jp_recopt and i1>=0 and _jetids[i1] and jtpt[i1]>_jp_recopt) { // Leading jets
    double dphi = DPhi(jtphi[i0], jtphi[i1]);
    if (dphi > 2.7) { // Back-to-back condition
      double pt3 = ((i2>=0 and jtpt[i2]>_jp_recopt) ? jtpt[i2] : 0.);
      double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
      double alpha = pt3/ptave;

      for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
        int itag = jt3leads[itag_lead];
        int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
        double etatag = jteta[itag];
        double etaprobe = jteta[iprobe];
        if (etatag < 1.3) {
          if (_jp_do3dHistos) {
            double pttag = _pt[itag];
            double ptprobe = _pt[iprobe];
            double alphatp = pt3/pttag;
            double asymm = (ptprobe - pttag)/(2*ptave);
            double asymmtp = (ptprobe - pttag)/(2*pttag);
            double ptresp_tag = _recopt[itag]/_genpt[itag];
            double ptresp_probe = _recopt[iprobe]/_genpt[iprobe];
            for (unsigned alphaidx = 0; alphaidx < h->alpharange.size(); ++alphaidx) {
              float alphasel = h->alpharange[alphaidx];
              if (alphatp<alphasel) {
                h->hdjasymmtp     [alphaidx]->Fill(pttag, etaprobe, asymmtp,      _w);
                h->hdjresptp_tag  [alphaidx]->Fill(pttag, etaprobe, ptresp_tag,   _w);
                h->hdjresptp_probe[alphaidx]->Fill(pttag, etaprobe, ptresp_probe, _w);
              } // alphatp sel.
              if (alpha<alphasel) {
                h->hdjasymm     [alphaidx]->Fill(ptave, etaprobe, asymm,        _w);
                h->hdjresp_tag  [alphaidx]->Fill(ptave, etaprobe, ptresp_tag,   _w);
                h->hdjresp_probe[alphaidx]->Fill(ptave, etaprobe, ptresp_probe, _w);
              } // alpha sel.
            } // alpha loop
          } // 3D histos?
        } // etatag < 1.3
      } // itag
    } // dphi > 2.7
  } // two or more jets, phase space
} // FillSingleEta


// Write and delete histograms
void HistosFill::WriteMC()
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteMC():");
  PrintMemInfo();

  for (auto histit : _mchistos) {
    for (auto histidx = 0u; histidx != histit.second.size(); ++histidx) {
      HistosMC *h = histit.second[histidx];
      delete h;
    } // for histidx
  } // for histit

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteMC() finished:");
  PrintMemInfo();
} // WriteMC


// Initialize basic histograms for trigger and eta bins
void HistosFill::InitRun(string name, double etamin, double etamax) {

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitRun(%s):",name.c_str()));
  PrintMemInfo();

  TDirectory *curdir = gDirectory;
  TFile *f = _outfile;
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  bool enteroutdir = f->cd(name.c_str());
  assert(enteroutdir);
  TDirectory *dir = f->GetDirectory(name.c_str()); assert(dir);
  dir->cd();

  HistosRun *h = new HistosRun(dir, etamin, etamax);
  _runhistos[name] = h;

  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("InitRun() finished:");
  PrintMemInfo();
} // InitRun


// Fill run histograms
void HistosFill::FillRun(string name)
{
  HistosRun *h = _runhistos[name];
  assert(h);

  // Luminosity information
  if (_jp_isdt and h->lums[run][lbn]==0) {

    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    // Let lum==0 pass, it can be a Poisson fluctuation for a valid LS

    h->lumsum += lum;
    h->lumsum2 += lum2;
    h->runlums[run] += lum;
    h->runlums2[run] += lum2;
    h->lums[run][lbn] = 1;

    for (auto trgidx = 0u; trgidx != h->trg.size(); ++trgidx) {
      string const& t = h->trg[trgidx];
      double prescale(0);
      if (_prescales[t].find(run)==_prescales[t].end()) {
        if (_trigs.find(t)!=_trigs.end()) {
          PrintInfo(Form("Prescale not found for trigger %s run %d",
                         t.c_str(),run));
          assert(false);
        }
      } else prescale = _prescales[t][run];
      h->runlums_trg[t][run] += (prescale ? lum / prescale : 0.);
    } // for trgidx

    // Initialize counters for a new run
    if (h->lums[run].size()==1) {
      for (auto trgidx = 0u; trgidx != h->trg.size(); ++trgidx) {
        string const& t = h->trg[trgidx];

        h->p_trg[t][run] = 0;
        h->t_trg[t][run] = 0;
        h->npv_trg[t][run] = 0;
        h->c_chf[t][run] = 0;
        h->c_nef[t][run] = 0;
        h->c_nhf[t][run] = 0;
        h->c_betastar[t][run] = 0;
        h->c_betaprime[t][run] = 0;
        h->t_trgtp[t][run] = 0;
        h->c_chftp[t][run] = 0;
        h->c_neftp[t][run] = 0;
        h->c_nhftp[t][run] = 0;
        h->c_betastartp[t][run] = 0;
        h->c_betaprimetp[t][run] = 0;

        for (unsigned int j = 0; j != h->trg.size(); ++j) {
          string const& t2 = h->trg[j];
          h->p_trgpair[t+t2][run] = 0;
        } // for j
      } // for trgidx
    } // new run
  } // data, luminosity zero

  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  double dphi = (njt>=2 ? DPhi(jtphi[i0], jtphi[i1]) : 0.);
  double pt3 = (njt>=3 ? jtpt[i2] : 0.);

  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    double pt = jtpt[jetidx];
    double eta = jteta[jetidx];

    if (h->etamin <= fabs(eta) and fabs(eta) < h->etamax and _pass and _jetids[jetidx] and _pass_qcdmet) {
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit) {
        string const& t = *trgit;

        if (pt > 18.) ++h->p_trg[t][run];
        if (pt > h->pt[t]) {
          ++h->t_trg[t][run]; // unweighted events
          h->tw_trg[t][run] += _prescales[t][run]; // prescale weighted events
          h->npv_trg[t][run] += npv;
          h->npvgood_trg[t][run] += npvgood;
          h->c_chf[t][run] += jtchf[jetidx];
          h->c_nef[t][run] += (jtnef[jetidx]-jthef[jetidx]);
          h->c_nhf[t][run] += (jtnhf[jetidx]-jthhf[jetidx]);
          h->c_betastar[t][run] += jtbetastar[jetidx];
          h->c_betaprime[t][run] += jtbetaprime[jetidx];
        }

        int itag = (jetidx==0 ? 1 : 0);
        if (jetidx<2 and dphi > 2.7 and pt3 < jtpt[itag] and fabs(jteta[itag]) < 1.3 and
            jtpt[itag] > h->pt[t] and _jetids[itag]) {
          ++h->t_trgtp[t][run];
          h->c_chftp[t][run] += jtchf[jetidx];
          h->c_neftp[t][run] += (jtnef[jetidx]-jthef[jetidx]);
          h->c_nhftp[t][run] += (jtnhf[jetidx]-jthhf[jetidx]);
          h->c_betastartp[t][run] += jtbetastar[jetidx];
          h->c_betaprimetp[t][run] += jtbetaprime[jetidx];
        }

        for (auto trgit2 = _trigs.begin(); trgit2 != _trigs.end(); ++trgit2) {
          string const& t2 = *trgit2;
          if (t!=t2) ++h->p_trgpair[t+t2][run];
        } // for trgit2
      } // for trgit
    } // conditions
  } // for jetidx
} // FillRun


// Write and delete histograms
void HistosFill::WriteRun()
{
  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteRun():");
  PrintMemInfo();

  for (auto histit = _runhistos.begin(); histit != _runhistos.end(); ++histit) {
    HistosRun *h = histit->second;
    delete h;
  } // for histit

  PrintInfo(Form("\nOutput (HistosRun) stored in %s",_outfile->GetName()),true);

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteRun() finished:");
  PrintMemInfo();
} // WriteRun


void HistosFill::FillJetID(vector<bool> &id)
{
  assert(int(id.size())==njt);

  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    id[jetidx] = ((fabs(jteta[jetidx])<2.5 ? jtidtight[jetidx] : jtidloose[jetidx]));

    if (_jp_doVetoECAL) {
      assert(_ecalveto);
      int ibin = _ecalveto->FindBin(jteta[jetidx],jtphi[jetidx]);
      id[jetidx] = (id[jetidx] and _ecalveto->GetBinContent(ibin)==0);
    }
  }
} // FillJetID


// Load good run and LS information
bool HistosFill::LoadJSON(const char* filename)
{
  cout << "Processing LoadJSON(\"" << filename << "\"..." << endl;
  ifstream file(filename, ios::in);
  if (!file.is_open()) return false;
  char c;
  string s, s2;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s && sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (_jp_debug)
      cout << "\"" << rn << "\": " << flush;

    while (file.get(c) && c==' ') {};
    if (_jp_debug) { cout << c << flush; assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun && file >> s >> s2 && sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (_jp_debug)
        cout << "["<<ls1<<","<<ls2<<"]"<<s1 << flush;

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (!endrun && s2!=",") {
        if (_jp_debug) { cout<<"s1: "<<s2<<endl<<flush; assert(s2==","); }
      }
    } // while ls
    if (_jp_debug)
      cout << endl;

    if (s2=="]}") continue;
    else if (s2!="],") {
      if (_jp_debug) { cout<<"s2: "<<s2<<endl<<flush; assert(s2=="],"); }
    }
  } // while run
  if (s2!="]}") { cout<<"s3: "<<s2<<endl<<flush; return false; }

  cout << "Called LoadJSON(\"" << filename << "\"):" << endl;
  cout << "Loaded " << nrun << " good runs and " << nls
       << " good lumi sections" << endl;
  return true;
} // LoadJSON


// Load luminosity information
bool HistosFill::LoadLumi(const char* filename)
{
  cout << "Processing LoadLumi(\"" << filename << "\")..." << endl;

  // Check lumi against the list of good runs
  const int a_goodruns[] = {};
  const int ngoodruns = sizeof(a_goodruns)/sizeof(a_goodruns[0]);
  set<int> goodruns;
  if (ngoodruns>0) { // This is an old remnant
    for (int runidx = 0; runidx != ngoodruns; ++runidx)
      goodruns.insert(a_goodruns[runidx]);

    for (auto runit = goodruns.begin(); runit != goodruns.end(); ++runit)
      cout << *runit << ", ";
    cout << endl;
  }
  set<pair<int, int> > nolums;

  ifstream f(filename, ios::in);
  if (!f.is_open()) return false;
  float secLS = 2.3310e+01;
  string s;
  int rn, fill, ls, ifoo;
  float del, rec, avgpu, energy;
  char sfoo[512];
  bool getsuccess1 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess1) return false;
  cout << endl << "string: " << s << " !" << endl << flush;

  // HOX: the lumi file format has been changing. Change the conditions when needed.
  if (s!="#Data tag : v2 , Norm tag: None") return false;

  bool getsuccess2 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess2) return false;
  cout << endl << "string: " << s << " !" << endl << flush;
  if (s!="#run:fill,ls,time,beamstatus,E(GeV),delivered(/ub),recorded(/ub),avgpu,source") return false;

  int nls(0);
  double lumsum(0);
  double lumsum_good(0);
  double lumsum_json(0);
  bool skip(false);
  while (getline(f, s, '\n')) {
    // Skip if not STABLE BEAMS or wrong number of arguments
    // STABLE BEAMS alts: ADJUST, BEAM DUMP, FLAT TOP, INJECTION PHYSICS BEAM, N/A, RAMP DOWN, SETUP, SQUEEZE
    if (sscanf(s.c_str(),"%d:%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,%f,%f,%f,%f,%s",
        &rn,&fill,&ls,&ifoo, &ifoo,&ifoo,&ifoo,&ifoo,&ifoo,&ifoo, &energy,&del,&rec,&avgpu,sfoo)!=15)
      skip=true;

    if (_jp_debug)
      cout << "Run " << rn << " ls " << ls << " lumi " << rec*1e-6 << "/pb" << endl;

    if (skip) { // The user should know if this happens, since we can choose to use only STABLE BEAMS
      if (skip) cout << "Skipping line (effects the recorded lumi):\n" << s << endl;
      skip = false;
      continue;
    }

    if (_lums[rn][ls]!=0) return false;
    if (_avgpu[rn][ls]!=0) return false;
    // lumiCalc.py returns lumi in units of mub-1 (=>nb-1=>pb-1)
    double lum = rec*1e-6;
    double lum2 = del*1e-6;
    if (lum==0 and goodruns.find(rn)!=goodruns.end() and (!_jp_dojson or _json[rn][ls]==1))
      nolums.insert(pair<int, int>(rn,ls));

    _avgpu[rn][ls] = avgpu; // * 69000. / 78400.; // brilcalc --minBiasXsec patch
    _lums[rn][ls] = lum;
    _lums2[rn][ls] = lum2;
    lumsum += lum;
    if (goodruns.find(rn)!=goodruns.end()) // Apr 17
      lumsum_good += lum;
    if ((!_jp_dojson || _json[rn][ls]))
      lumsum_json += lum;
    ++nls;
    if (nls>100000000) return false;
  }

  cout << "Called LoadLumi(\"" << filename << "\"):" << endl;
  cout << "Loaded " << _lums.size() << " runs with "
       << nls << " lumi sections containing "
       << lumsum << " pb-1 of data,\n of which "
       << lumsum_good << " pb-1 is in good runs ("
       << 100.*lumsum_good/lumsum << "%)"<< endl;
  cout << "This corresponds to " << nls*secLS/3600
       << " hours of data-taking" << endl;
  cout << "The JSON file contains "
       << lumsum_json << " pb-1 ("
       << 100.*lumsum_json/lumsum << "%)"<< endl;

  // Report any empty lumi section
  if (nolums.size()!=0) {
    cout << "Warning, found " << nolums.size() << " non-normalizable LS:";
    for (auto lumit = nolums.begin();
         lumit != nolums.end(); ++lumit) {

      cout << " ["<<lumit->first<<","<<lumit->second;
      auto lumit2 = lumit;
      ++lumit2;
      if (lumit2->first!=lumit->first or lumit2->second!=lumit->second+1)
        cout << "]";
      else {
        for (int lumadd = 0; lumit2!=nolums.end() and
                             lumit2->first==lumit->first and
                             lumit2->second==lumit->second+lumadd+1; ++lumadd, ++lumit2) {};
        lumit = --lumit2;
        cout << "-" << lumit->second << "]";
      }
    } // for lumit
    cout << endl;
  } // nolums
  return true;
} // LoadLumi


bool HistosFill::LoadPuProfiles(const char *datafile, const char *mcfile)
{
  cout << "Processing LoadPuProfiles(\"" << datafile << "\",\"" << mcfile << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  // Load pile-up files and hists from them
  TFile *f_pudist = new TFile(datafile, "READ");
  if (!f_pudist or f_pudist->IsZombie()) return false;
  TFile *fpumc = new TFile(mcfile,"READ");
  if (!fpumc or fpumc->IsZombie()) return false;

  _pumc = dynamic_cast<TH1D*>(fpumc->Get("pileupmc"));
  if (!_pumc) return false;
  double maxmcpu = _pumc->GetMaximum();
  int lomclim = _pumc->FindFirstBinAbove(maxmcpu/1000.0);
  int upmclim = _pumc->FindLastBinAbove(maxmcpu/1000.0);
  int maxmcbin = _pumc->FindFirstBinAbove(0.999*maxmcpu);
  PrintInfo(Form("Maximum bin: %d for MC",maxmcbin),true);
  PrintInfo(Form("Hazardous pu below & above: %f, %f",_pumc->GetBinLowEdge(lomclim),_pumc->GetBinLowEdge(upmclim+1)),true);
  int nbinsmc = _pumc->GetNbinsX();
  int kmc = _pumc->FindBin(33); // Check that pu=33 occurs at the same place as for data

  // For data, load each trigger separately
  for (auto itrg = 0u ; itrg != _jp_notrigs; ++itrg) {
    string t = _jp_triggers[itrg];
    _pudist[t] = dynamic_cast<TH1D*>(f_pudist->Get(t.c_str()));
    if (!_pudist[t]) return false;
    int nbinsdt = _pudist[t]->GetNbinsX();
    int kdt = _pudist[t]->FindBin(33);
    if (kdt!=kmc or nbinsdt!=nbinsmc) {
      cout << "The pileup histogram dt vs mc binning or range do not match (dt left mc right):" << endl;
      cout << " Bins: " << nbinsdt << " " << nbinsmc << endl;
      cout << " Pu=33 bin: " << kdt << " " << kmc << endl;
      return false;
    }
    double maxdtpu = _pudist[t]->GetMaximum();
    int lodtlim = _pudist[t]->FindFirstBinAbove(maxdtpu/1000.0);
    int updtlim = _pudist[t]->FindLastBinAbove(maxdtpu/1000.0);
    int maxdtbin = _pudist[t]->FindFirstBinAbove(0.999*maxdtpu);

    for (int bin = 0; bin < lomclim; ++bin) // Set fore-tail to zero
      _pudist[t]->SetBinContent(bin,0.0);
    PrintInfo(Form("Maximum bin: %d for DT trg %s",maxdtbin,t.c_str()),true);
    PrintInfo(Form("Hazardous pu below & above: %f, %f",_pudist[t]->GetBinLowEdge(lodtlim),_pudist[t]->GetBinLowEdge(updtlim+1)),true);
    _pudist[t]->Divide(_pumc);
  }
  PrintInfo("Finished processing pileup histos!",true);

  curdir->cd();
  return true;
} // LoadPuProfiles


bool HistosFill::LoadPrescales(const char *prescalefile)
{
  cout << "Processing LoadPrescales(\"" << prescalefile << "\")..." << endl;
  fstream fin(prescalefile);

  const int ns = 1024;
  char s[ns];

  fin.getline(s, ns);
  stringstream ss;
  ss << s;

  string srun, sls, strg;
  vector<string> trgs;
  ss >> srun;
  if (srun!="RUN") return false;

  while (ss >> strg) trgs.push_back(strg);

  int run, ls, pre;
  while(fin.getline(s, ns)) {
    stringstream ss;
    ss << s;
    ss >> run >> ls;
    if (_jp_debug) cout << " run " << run << " ls " << ls << ": ";

    for (unsigned int itrg = 0; ss >> pre; ++itrg) {
      if (itrg==trgs.size()) return false;
      if (_jp_debug) cout << pre << "/" << trgs[itrg] << " ";
    }
    if (_jp_debug) cout << endl;
  }
  return true;
} // LoadPrescales


bool HistosFill::LoadVetoECAL(const char *file)
{
  cout << "Processing LoadVetoECAL(\"" << file << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  TFile *fe = new TFile(file, "READ");
  if (!fe or fe->IsZombie()) return false;

  _ecalveto = (TH2F*)fe->Get("ecalveto");
  if (!_ecalveto) return false;

  curdir->cd();
  return true;
} // LoadVetoECAL


// Check that the correct tree is open in the chain
Long64_t HistosFill::LoadTree(Long64_t entry)
{
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;

  // A new tree is opened
  if (fChain->GetTreeNumber() != fCurrent) {

    fCurrent = fChain->GetTreeNumber();
    PrintInfo(Form("Opening tree number %d", fChain->GetTreeNumber()));

    if (_jp_isdt) {
      // Reload the triggers and print them
      if (!GetTriggers()) {
        PrintInfo("Failed to load DT triggers. Check that the SMPJ tuple has the required histograms. Aborting...");
        return -4;
      }
      PrintInfo(Form("Tree %d triggers:",fCurrent));
      for (auto trigi = 0u; trigi < _availTrigs.size(); ++trigi) {
        auto str = _availTrigs[trigi];
        *ferr << str;
        auto trgplace = std::find(_goodTrigs.begin(),_goodTrigs.end(),trigi);
        if (trgplace!=_goodTrigs.end()) *ferr << "_y (" << _goodWgts[trgplace-_goodTrigs.begin()] << ") ";
        else *ferr << "_n ";
        if (trigi%(_jp_notrigs+1)==_jp_notrigs) *ferr << endl;
      }
      *ferr << endl << flush;
      if (_jp_doVetoHCALHot) {
        // Reload hot HCAL zones
        bool worryHCALHotExcl = false;
        PrintInfo("Looking for HCAL exclusion zones...");
        for (auto era = 0u; era < _jp_HCALHotEras; ++era) {
          std::regex HCALEraName(_jp_HCALHotRuns[era]);
          if (std::regex_search(_jp_run,HCALEraName)) {
            PrintInfo(Form("Doing HCAL exclusion for era %s at : (etamin, etamax, phimin, phimax)",_jp_HCALHotRuns[era]));
            for (auto idxetaphi = 0u; idxetaphi < 4; ++idxetaphi) {
              *ferr << "  " << _jp_HCALHotRanges[era][idxetaphi];
              rangesHCALHotExcl[idxetaphi] =  _jp_HCALHotRanges[era][idxetaphi];
            }
            *ferr << endl;
            worryHCALHotExcl = true;
            break;
          }
        }
      }
    } else if (_jp_pthatbins) {
      TString filename = fChain->GetCurrentFile()->GetName();
      // Check the position of the current file in the list of file names
      unsigned currFile = std::find_if(_jp_pthatfiles.begin(),_jp_pthatfiles.end(),[&filename] (string s) { return filename.Contains(s); })-_jp_pthatfiles.begin();
      if (_jp_pthatnevts[currFile]<=0.0 or _jp_pthatsigmas[currFile]<=0) {
        PrintInfo(Form("Suspicious pthat slice information for file %s. Aborting...",_jp_pthatfiles[currFile].c_str()));
        return -3;
      }
      _pthatweight = _jp_pthatsigmas[currFile]/_jp_pthatnevts[currFile];
      _pthatweight /= (_jp_pthatsigmas[_jp_npthatbins-1]/_jp_pthatnevts[_jp_npthatbins-1]); // Normalize
      _pthatuplim = _jp_pthatuplims[currFile];
      PrintInfo(Form("Pthat bin changing.\nFile %d %s should correspond to the range [%f,%f]\nWeight: %f & uplim: %f",
                     currFile,fChain->GetCurrentFile()->GetName(),_jp_pthatranges[currFile],_jp_pthatranges[currFile+1],_pthatweight,_pthatuplim));;
    }
    // slices with pthat bins
  }
  return centry;
}


// Update the available trigger types for each new tree
bool HistosFill::GetTriggers()
{
  TH1F *triggers = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerNames")); assert(triggers);
  TAxis *xax = triggers->GetXaxis();

  regex pfjet("HLT_PFJet([0-9]*)_v([0-9]*)");
  regex ak8("HLT_AK8PFJet([0-9]*)_v[0-9]*");
  regex jetht("HLT_PFHT([0-9]*)_v[0-9]*");

  _availTrigs.clear();
  _goodTrigs.clear();
  _goodWgts.clear();
  for (int trgidx = xax->GetFirst(); trgidx <= xax->GetLast(); ++trgidx) {
//     cout << trgidx << " " << xax->GetLast() << endl;
    string trgName = xax->GetBinLabel(trgidx);
    if (trgName.compare("")==0) continue; // Ignore empty places on x-axis
    string trigger = "x"; // HLT_PFJet are given non-empty trigger names
    if (std::regex_match(trgName,pfjet)) {
      trigger=std::regex_replace(trgName, pfjet, "jt$1", std::regex_constants::format_no_copy);
      _goodTrigs.push_back(_availTrigs.size());
      double trigthr = std::stod(std::regex_replace(trgName, pfjet, "$1", std::regex_constants::format_no_copy));
      unsigned int thrplace = std::find(_jp_trigthr,_jp_trigthr+_jp_notrigs,trigthr)-_jp_trigthr;
      if (_jp_useversionlumi) {
        string trgdummy = std::regex_replace(trgName, pfjet, "$1_$2", std::regex_constants::format_no_copy);
        bool found = false;
        for (auto IDidx = 0u; IDidx < _jp_notrigIDs; ++IDidx) {
          for (auto &currTag : _jp_trigtags[IDidx]) {
            if (std::regex_match(trgdummy,currTag)) {
              if (found) {
                cerr << "Double match for " << trgName << ", check trigger naming in _jp_trigtags. Aborting..." << endl;
                return false;
              }
              found = true;
              if (thrplace < _jp_notrigs+1)
                _goodWgts.push_back(_jp_trigwgts[IDidx][thrplace]/_jp_trigwgts[_jp_notrigIDs][thrplace]);
              else {
                cerr << "Error searching the trigger weight! Aborting..." << endl;
                return false;
              }
            }
          }
        }
        if (!found) {
          cerr << "No info for " << trgName << " in _jp_trigtag. Could be a dummy trigger in the tuple." << endl;
          _goodWgts.push_back(0.0);
        }
      } else {
        _goodWgts.push_back(1.0);//_jp_triglumi[_jp_notrigs-1]/_jp_triglumi[thrplace]);
      }
    } else if (std::regex_match(trgName,ak8)) {
      trigger=std::regex_replace(trgName, ak8, "ak8jt$1", std::regex_constants::format_no_copy);
    } else if (std::regex_match(trgName,jetht)) {
      trigger=std::regex_replace(trgName, jetht, "jetht$1", std::regex_constants::format_no_copy);
    } else {
      trigger="x";
      cout << "Unknown trigger type " << trgName << endl;
    }
    _availTrigs.push_back(trigger);
  }

  return _availTrigs.size()>0;
} // GetTriggers
