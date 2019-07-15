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
#ifdef NEWMODE
  met(PFMetT0__et_),
  metphi(PFMetT0__phi_),
  metsumet(PFMetT0__sumEt_),
  met01(PFMetT0T1__et_),
  metsumet01(PFMetT0T1__sumEt_)
#else
  met(PFMet__et_),
  metphi(PFMet__phi_),
  metsumet(PFMet__sumEt_)
#endif
{
  assert(tree);
  _initsuccess = Init(tree);
}


// Print the same info on a file and conditionally to the output
void HistosFill::PrintInfo(string info, bool printcout)
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
bool HistosFill::Init(TChain *tree)
{
  ferr = 0;
  ferr = new ofstream(Form("reports/HistosFill-%s.log",jp::type),ios::out);

  if (!tree) {
    PrintInfo("No tree given, processing makes no sense.",true);
    return false; // With no tree, Loop will be interrupted
  }
  string type = jp::type;
  if (jp::pthatbins) type = "Pthat";
  else if (jp::htbins) type = "HT";
  _outfile = new TFile(Form("output-%s-1.root",type.c_str()), "NEW");
  if (!_outfile or _outfile->IsZombie()) {
    PrintInfo(Form("Opening the output file output-%s-1.root failed. Check if the file already exists.",jp::type),false);
    return false;
  }

  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("EvtHdr_.mIsPVgood", &EvtHdr__mIsPVgood);
  fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun);
  fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent);
  fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi);
  fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch);
  fChain->SetBranchAddress("EvtHdr_.mNVtx", &EvtHdr__mNVtx);
  fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood);
  fChain->SetBranchAddress("EvtHdr_.mOOTPUEarly", &EvtHdr__mOOTPUEarly);
  fChain->SetBranchAddress("EvtHdr_.mOOTPULate", &EvtHdr__mOOTPULate);
  fChain->SetBranchAddress("EvtHdr_.mINTPU", &EvtHdr__mINTPU);
  fChain->SetBranchAddress("EvtHdr_.mNBX", &EvtHdr__mNBX);
  fChain->SetBranchAddress("EvtHdr_.mPVndof", &EvtHdr__mPVndof);
  fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu);
  fChain->SetBranchAddress("EvtHdr_.mPVx", &EvtHdr__mPVx);
  fChain->SetBranchAddress("EvtHdr_.mPVy", &EvtHdr__mPVy);
  fChain->SetBranchAddress("EvtHdr_.mPVz", &EvtHdr__mPVz);
  fChain->SetBranchAddress("EvtHdr_.mBSx", &EvtHdr__mBSx);
  fChain->SetBranchAddress("EvtHdr_.mBSy", &EvtHdr__mBSy);
  fChain->SetBranchAddress("EvtHdr_.mBSz", &EvtHdr__mBSz);
  fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat);
  fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight);
  fChain->SetBranchAddress("EvtHdr_.mCaloRho", &EvtHdr__mCaloRho);
  fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho);
#ifdef NEWMODE
  fChain->SetBranchAddress("PFMetT0_.et_", &PFMetT0__et_);
  fChain->SetBranchAddress("PFMetT0_.sumEt_", &PFMetT0__sumEt_);
  fChain->SetBranchAddress("PFMetT0_.phi_", &PFMetT0__phi_);
  fChain->SetBranchAddress("PFMetT0T1_.et_", &PFMetT0T1__et_);
  fChain->SetBranchAddress("PFMetT0T1_.sumEt_", &PFMetT0T1__sumEt_);
  fChain->SetBranchAddress("PFMetT0T1_.phi_", &PFMetT0T1__phi_);
#endif
  fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_);
  fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_);
  fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_);
  if (jp::fetchMETFilters and jp::doMETFiltering) fChain->SetBranchAddress("FilterDecision_", &FilterDecision_);
  fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_);
  fChain->SetBranchAddress("L1Prescale_", &L1Prescale_);
  fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_);
  fChain->SetBranchAddress("GenJets_", &GenJets__);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fX", GenJets__fCoordinates_fX);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fY", GenJets__fCoordinates_fY);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", GenJets__fCoordinates_fZ);
  fChain->SetBranchAddress("GenJets_.fCoordinates.fT", GenJets__fCoordinates_fT);
#ifdef NEWMODE
  fChain->SetBranchAddress("PFJetsCHS_.genIdx_", PFJetsCHS__genIdx_);
#else
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fX",jp::chs), PFJetsCHS__genP4__fCoordinates_fX);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fY",jp::chs), PFJetsCHS__genP4__fCoordinates_fY);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fZ",jp::chs), PFJetsCHS__genP4__fCoordinates_fZ);
  fChain->SetBranchAddress(Form("PFJets%s_.genP4_.fCoordinates.fT",jp::chs), PFJetsCHS__genP4__fCoordinates_fT);
#endif
#ifdef NEWMODE
  fChain->SetBranchAddress(Form("PFJets%s_.QGL_",jp::chs), PFJetsCHS__QGL_); // qgl
#else
  fChain->SetBranchAddress(Form("PFJets%s_.QGtagger_",jp::chs), PFJetsCHS__QGtagger_); // qgl
#endif
  fChain->SetBranchAddress(Form("PFJets%s_",jp::chs), &PFJetsCHS__);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fX",jp::chs), PFJetsCHS__P4__fCoordinates_fX);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fY",jp::chs), PFJetsCHS__P4__fCoordinates_fY);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fZ",jp::chs), PFJetsCHS__P4__fCoordinates_fZ);
  fChain->SetBranchAddress(Form("PFJets%s_.P4_.fCoordinates.fT",jp::chs), PFJetsCHS__P4__fCoordinates_fT);
  fChain->SetBranchAddress(Form("PFJets%s_.genR_",jp::chs), PFJetsCHS__genR_);
  fChain->SetBranchAddress(Form("PFJets%s_.cor_",jp::chs), PFJetsCHS__cor_);
  fChain->SetBranchAddress(Form("PFJets%s_.area_",jp::chs), PFJetsCHS__area_);
  fChain->SetBranchAddress(Form("PFJets%s_.looseID_",jp::chs), PFJetsCHS__looseID_);
  fChain->SetBranchAddress(Form("PFJets%s_.tightID_",jp::chs), PFJetsCHS__tightID_);
  fChain->SetBranchAddress(Form("PFJets%s_.partonFlavour_",jp::chs), PFJetsCHS__partonFlavour_);
  fChain->SetBranchAddress(Form("PFJets%s_.hadronFlavour_",jp::chs), PFJetsCHS__hadronFlavour_);
#ifdef NEWMODE
  fChain->SetBranchAddress(Form("PFJets%s_.partonFlavourPhysicsDef_",jp::chs), PFJetsCHS__partonFlavourPhysicsDef_);
#endif
  fChain->SetBranchAddress(Form("PFJets%s_.chf_",jp::chs), PFJetsCHS__chf_);
  fChain->SetBranchAddress(Form("PFJets%s_.nhf_",jp::chs), PFJetsCHS__nhf_);
  fChain->SetBranchAddress(Form("PFJets%s_.nemf_",jp::chs), PFJetsCHS__nemf_);
  fChain->SetBranchAddress(Form("PFJets%s_.cemf_",jp::chs), PFJetsCHS__cemf_);
  fChain->SetBranchAddress(Form("PFJets%s_.muf_",jp::chs), PFJetsCHS__muf_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_hf_",jp::chs), PFJetsCHS__hf_hf_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_phf_",jp::chs), PFJetsCHS__hf_phf_);
  fChain->SetBranchAddress(Form("PFJets%s_.chm_",jp::chs), PFJetsCHS__chm_);
  fChain->SetBranchAddress(Form("PFJets%s_.nhm_",jp::chs), PFJetsCHS__nhm_);
  fChain->SetBranchAddress(Form("PFJets%s_.phm_",jp::chs), PFJetsCHS__phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.elm_",jp::chs), PFJetsCHS__elm_);
  fChain->SetBranchAddress(Form("PFJets%s_.mum_",jp::chs), PFJetsCHS__mum_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_hm_",jp::chs), PFJetsCHS__hf_hm_);
  fChain->SetBranchAddress(Form("PFJets%s_.hf_phm_",jp::chs), PFJetsCHS__hf_phm_);
  fChain->SetBranchAddress(Form("PFJets%s_.ncand_",jp::chs), PFJetsCHS__ncand_);
  fChain->SetBranchAddress(Form("PFJets%s_.betaPrime_",jp::chs), PFJetsCHS__betaPrime_);
  fChain->SetBranchAddress("genFlavour_", &genFlavour_);
  fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_);
#ifdef NEWMODE
  fChain->SetBranchAddress("genFlavourPartonPhysicsDef_", &genFlavourPartonPhysicsDef_);
#endif

  if (jp::quick) { // Activate only some branches
    fChain->SetBranchStatus("*",0);
    // Luminosity calculation
    if (jp::ismc) fChain->SetBranchStatus("EvtHdr_.mPthat",1); // pthat
    if (jp::ismc) fChain->SetBranchStatus("EvtHdr_.mWeight",1); // weight
    if (jp::isdt) fChain->SetBranchStatus("EvtHdr_.mRun",1); // run
    if (jp::isdt) fChain->SetBranchStatus("EvtHdr_.mEvent",1); // evt
    if (jp::isdt) fChain->SetBranchStatus("EvtHdr_.mLumi",1); // lbn

    // Event properties
    fChain->SetBranchStatus("EvtHdr_.mNVtx",1); // npv
    fChain->SetBranchStatus("EvtHdr_.mNVtxGood",1); // npvgood
    fChain->SetBranchStatus("EvtHdr_.mPFRho",1); // rho

    // Jet properties (jtpt, jte, jteta, jty, jtphi etc.)
    fChain->SetBranchStatus(Form("PFJets%s_",jp::chs),1); // njt
    fChain->SetBranchStatus(Form("PFJets%s_.P4_*",jp::chs),1); // jtp4*
    fChain->SetBranchStatus(Form("PFJets%s_.cor_",jp::chs),1); // jtjes
    fChain->SetBranchStatus(Form("PFJets%s_.area_",jp::chs),1); // jta

    if (jp::ismc) {

      fChain->SetBranchStatus(Form("PFJets%s_.genR_",jp::chs),1); // jtgenr
#ifdef NEWMODE
      fChain->SetBranchStatus(Form("PFJets%s_.genIdx_",jp::chs),1); // jtgenidx
#else
      fChain->SetBranchStatus(Form("PFJets%s_.genP4_*",jp::chs),1); // jtgenp4*
#endif
    }

    // for quark/gluon study (Ozlem)
#ifdef NEWMODE
    fChain->SetBranchStatus(Form("PFJets%s_.QGL_",jp::chs),1); // qgl
#else
    fChain->SetBranchStatus(Form("PFJets%s_.QGtagger_",jp::chs),1); // qgl
#endif
    if (jp::ismc) fChain->SetBranchStatus(Form("PFJets%s_.partonFlavour_",jp::chs),1);
    if (jp::ismc) fChain->SetBranchStatus(Form("PFJets%s_.partonFlavourPhysicsDef_",jp::chs),1);

    // Component fractions
    fChain->SetBranchStatus(Form("PFJets%s_.chf_",jp::chs),1); // jtchf
    fChain->SetBranchStatus(Form("PFJets%s_.nemf_",jp::chs),1); // jtnef
    fChain->SetBranchStatus(Form("PFJets%s_.nhf_",jp::chs),1); // jtnhf
    fChain->SetBranchStatus(Form("PFJets%s_.cemf_",jp::chs),1); // jtcef !!
    fChain->SetBranchStatus(Form("PFJets%s_.muf_",jp::chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_hf_",jp::chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_phf_",jp::chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.ncand_",jp::chs),1); // jtn
    fChain->SetBranchStatus(Form("PFJets%s_.betaPrime_",jp::chs),1); // jtbetaprime
    fChain->SetBranchStatus(Form("PFJets%s_.chm_",jp::chs),1); // jtnch
    fChain->SetBranchStatus(Form("PFJets%s_.phm_",jp::chs),1); // jtnne
    fChain->SetBranchStatus(Form("PFJets%s_.nhm_",jp::chs),1); // jtnnh
    fChain->SetBranchStatus(Form("PFJets%s_.elm_",jp::chs),1); // jtnce !!
    fChain->SetBranchStatus(Form("PFJets%s_.mum_",jp::chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_hm_",jp::chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.hf_phm_",jp::chs),1); // jtnmu !!
    fChain->SetBranchStatus(Form("PFJets%s_.tightID_",jp::chs),1); // jtidtight
    fChain->SetBranchStatus(Form("PFJets%s_.looseID_",jp::chs),1); // jtidloose

    //fChain->SetBranchStatus("rho",1);
#ifdef NEWMODE
    fChain->SetBranchStatus("PFMetT0_.et_",1); // met
    fChain->SetBranchStatus("PFMetT0_.phi_",1); // metphi
    fChain->SetBranchStatus("PFMetT0_.sumEt_",1); // metsumet
    fChain->SetBranchStatus("PFMetT0T1_.et_",1); // met
    fChain->SetBranchStatus("PFMetT0T1_.phi_",1); // metphi
    fChain->SetBranchStatus("PFMetT0T1_.sumEt_",1); // metsumet
#else
    fChain->SetBranchStatus("PFMet_.et_",1); // met
    fChain->SetBranchStatus("PFMet_.phi_",1); // metphi
    fChain->SetBranchStatus("PFMet_.sumEt_",1); // metsumet
#endif

    if (jp::fetchMETFilters and jp::doMETFiltering) fChain->SetBranchStatus("FilterDecision_",1);
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

    if (jp::ismc) {
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
  jtbetaprime = &PFJetsCHS__betaPrime_[0];
  jtidloose = &PFJetsCHS__looseID_[0];
  jtidtight = &PFJetsCHS__tightID_[0];
  //
  jtgenr = &PFJetsCHS__genR_[0];
#ifdef NEWMODE
  jtgenidx = &PFJetsCHS__genIdx_[0];
#else
  jtgenp4x = &PFJetsCHS__P4__fCoordinates_fX[0];
  jtgenp4y = &PFJetsCHS__P4__fCoordinates_fY[0];
  jtgenp4z = &PFJetsCHS__P4__fCoordinates_fZ[0];
  jtgenp4t = &PFJetsCHS__P4__fCoordinates_fT[0];
#endif
  //

  partonflavor = &PFJetsCHS__partonFlavour_[0];
#ifdef NEWMODE
  partonflavorphys = &PFJetsCHS__partonFlavourPhysicsDef_[0];
  // for quark/gluon study (Ozlem)
  qgl = &PFJetsCHS__QGL_[0];
#else
  // for quark/gluon study (Ozlem)
  qgl = &PFJetsCHS__QGtagger_[0];
#endif

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

  // Put some logging info into the file
  _runinfo = "\n";
  vector<string> infofiles = {"settings.h","HistosFill.h","HistosFill.C","HistosBasic.h","HistosBasic.C","HistosEta.h","HistosEta.C","HistosAll.h","HistosAll.C","HistosNormalize.C","HistosCombine.C"};
  if (jp::gitinfo) {
    gROOT->ProcessLine(".!git log > gitlog.txt");
    infofiles.push_back("gitlog.txt");
  }
  for (auto &infofile : infofiles) {
    std::ifstream tmpstrm(infofile.c_str());
    string tmpstr;

    tmpstrm.seekg(0, std::ios::end);
    tmpstr.reserve(tmpstrm.tellg());
    tmpstrm.seekg(0, std::ios::beg);
    tmpstr.assign((std::istreambuf_iterator<char>(tmpstrm)),std::istreambuf_iterator<char>());
    _runinfo += Form("\n###I\n###I %s Block\n###I\n\n",infofile.c_str());
    _runinfo += tmpstr;
  }
  if (jp::gitinfo) gROOT->ProcessLine(".!rm gitlog.txt");

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
  for (Long64_t djentry=0; djentry<_nentries;djentry+=1+jp::skim) { // Event loop
    _jentry = djentry+_nskip; // Add a shift from beginning
    if (LoadTree(_jentry) < 0) break;
    nb = fChain->GetEntry(_jentry);   nbytes += nb;

    int testval = djentry;
    if (jp::skim>0) testval /= 1+jp::skim;

    if (testval%hopval==repval) { // 1M report (first report timed to be early)
      // Report memory usage to avoid malloc problems when writing file
      PrintInfo(Form("Doing Loop(), %dM events:",int(_jentry/1e6 + 0.5)),true);
      PrintMemInfo();
      PrintInfo(Form("\nProcessed %lld events (%1.1f%%) in %1.0f sec.",_jentry-1, 100.*djentry/_nentries, stop.RealTime()),true);
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

      if (jp::save) {
        for (auto &manyhists : _histos)
          for (auto &onehist : manyhists.second)
            onehist->Write();
          if (jp::doEtaHistos) {
            for (auto &manyhists : _etahistos)
              for (auto & onehist : manyhists.second)
                onehist->Write();
            if (jp::doEtaHistosMcResponse) {
              for (auto &manyhists : _mchistos)
                for (auto & onehist : manyhists.second)
                  onehist->Write();
            }
          }
      }
      stop.Continue();
    }

    if (!AcceptEvent()) continue;

    if (jp::debug) PrintInfo("Histos are being filled!",true);
    // Here can categorize events into different triggers, epochs, topologies etc.
    // Eta and pT binning are handled in the FillSingleBasic class
    if (jp::doBasicHistos) {
      FillBasic("Standard");
    }

    if (jp::doEtaHistos and _pass) {
      FillEta("FullEta_Reco", jtpt, jteta, jtphi);
      if (jp::ismc and jp::doEtaHistosMcResponse) {
        FillEta("FullEta_Gen", jtgenpt, jtgeneta, jtgenphi);
        FillMC("FullEta_RecoPerGen_vReco", jtpt, jtgenpt, jtpt,    jteta,    jtphi);
        FillMC("FullEta_RecoPerGen_vGen",  jtpt, jtgenpt, jtgenpt, jtgeneta, jtgenphi);
      }
    }

    // Run quality checks
    if (jp::isdt and jp::doRunHistos) {
      FillRun("Runs");
      FillRun("RunsBarrel");
      FillRun("RunsTransition");
      FillRun("RunsEndcap");
    }

    if (jp::isdt and jp::doMpfHistos and _pass) {
      FillAll("AllTrigs");
    }
  } // for jentry
  PrintInfo("",true);

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("Finished processing %lld entries:",_nentries),true);
  PrintMemInfo(true);

  if (jp::isdt and jp::doRunHistos)   WriteRun();
  if (jp::doEtaHistos)   WriteEta();
  if (jp::ismc and jp::doEtaHistos and jp::doEtaHistosMcResponse) WriteMC();
  if (jp::isdt and jp::doMpfHistos)   WriteAll();
  if (jp::doBasicHistos) WriteBasic(); // this needs to be last, output file closed

  Report();

  PrintInfo("Stopping processing at: ",true);
  TDatime now1;
  now1.Print();
  PrintInfo(Form("Processing used %f s CPU time (%f h)",stop.CpuTime(),stop.CpuTime()/3600.),true);
  PrintInfo(Form("Processing used %f s real time (%f h)",stop.RealTime(),stop.RealTime()/3600.),true);

  // Gather info from the report file
  std::ifstream tmpstrm(Form("reports/HistosFill-%s.log",jp::type));
  string tmpstr;

  tmpstrm.seekg(0, std::ios::end);
  tmpstr.reserve(tmpstrm.tellg());
  tmpstrm.seekg(0, std::ios::beg);
  tmpstr.assign((std::istreambuf_iterator<char>(tmpstrm)),std::istreambuf_iterator<char>());
  _runinfo += "\n###I\n###I Run Info Block\n###I\n\n";
  _runinfo += tmpstr;

  // Write RunInfo and close the file
  _outfile->mkdir("infodir");
  _outfile->cd("infodir");
  TNamed runinfo("info",_runinfo);
  runinfo.Write();
  _outfile->Close();

  PrintInfo("File closed, exiting at: ",true);
  TDatime now2;
  now2.Print();
  PrintInfo(Form("Processing and logistics used %f s CPU time (%f h)",stop.CpuTime(),stop.CpuTime()/3600.),true);
  PrintInfo(Form("Processing and logistics used %f s real time (%f h)",stop.RealTime(),stop.RealTime()/3600.),true);
  stop.Continue();
}


// Setup before event loop
bool HistosFill::PreRun()
{
  _nentries = fChain->GetEntriesFast();
  _ntot = fChain->GetEntries();
  _nskip = jp::nskip;//0;
  _nentries = (jp::nentries==-1 ? _ntot-_nskip : min(_ntot-_nskip, jp::nentries));
  assert(_nskip+_nentries);

  _nbadevts_json = _nbadevts_dup = _nbadevts_run = _nbadevts_ls = _nbadevts_lum = 0;
  _bscounter_bad = _bscounter_good = 0;
  _rhocounter_good = _rhocounter_bad = 0;
  _trgcounter = _evtcounter = _totcounter = 0;

  // Initialize _binnedmcweight. It will be loaded for each tree separately.
  if (jp::pthatbins or jp::htbins) {
    _binnedmcweight = 0;
    _binnedmcrepeats = 0;
  }

  PrintInfo("\nCONFIGURATION DUMP:");
  PrintInfo("-------------------");
  PrintInfo(Form("Running over %sPF",jp::algo));
  PrintInfo(Form("Running over %s",jp::ismc?"MC":"data"));
  PrintInfo(Form("%s time-dependent JEC (IOV)",jp::useIOV ? "Applying" : "Not applying"));
  PrintInfo(Form("%s Hot zone exclusion in eta-phi plane",jp::doVetoHot ? "Doing" : "Not doing"));
  *ferr << endl;

  _eraIdx = -1;
  if (jp::isdt) {
    PrintInfo(Form("%s additional JSON selection",jp::dojson ? "Applying" : "Not applying"));
    PrintInfo(Form("%s luminosities.",jp::dolumi ? "Recalculating" : "Not recalculating"));
    PrintInfo(Form("%s additional run-level histograms",jp::doRunHistos ? "Storing" : "Not storing"));
    PrintInfo(Form("%s basic set of histograms",jp::doBasicHistos ? "Storing" : "Not storing"));
    PrintInfo(Form("%s histograms with a full eta-range",jp::doEtaHistos ? "Storing" : "Not storing"));
    // Find out era index
    int eraNo = 0;
    for (auto &eraMatch : jp::eras) {
      if (std::regex_search(jp::run,eraMatch)) {
        _eraIdx = eraNo;
        break;
      }
      ++eraNo;
    }
    if (_eraIdx==-1) {
      PrintInfo(Form("Era not found for %s. Aborting!",jp::run),true);
      return false;
    } else {
      PrintInfo(Form("Matched %s with era index %d!",jp::run,_eraIdx),true);
    }
  } else if (jp::ismc) {
    PrintInfo(Form("%s pileup profile in MC to data",jp::reweighPU ? "Reweighting" : "Not reweighting"));
    PrintInfo(Form("Processing %s samples", jp::pthatbins ? "pThat binned" : "\"flat\""));
  }
  *ferr << endl;

  _JEC = 0;
  _L1RC = 0;
  _jecUnc = 0;

  // Time dependent JEC (only for dt)
  if (jp::redojes) {
    if (jp::isdt and jp::useIOV) {
      // If multiple IOV's are used, we set _JEC etc for each event separately, checking that the IOV is correct
      for (unsigned iovidx=0; iovidx<jp::IOVnames.size(); ++iovidx)
        _iov.addJEC(jp::IOVnames[iovidx],jp::IOVranges.at(iovidx).at(0),jp::IOVranges.at(iovidx).at(1));
    } else {
      // If only one great IOV is used, we can set _JEC etc. directly here.
      _iov.addJEC("");
      bool setcorrection = _iov.setJEC(&_JEC,&_L1RC,&_jecUnc);
      if (!setcorrection or !_JEC or !_L1RC or (jp::isdt and !_jecUnc)) {
        PrintInfo("Issues while loading JEC; aborting...",true);
        return false;
      }
    } // JEC redone
  }

  // Load latest JSON selection
  if (jp::isdt and jp::dojson) {
    if (!LoadJSON()) {
      PrintInfo("Issues loading the JSON file; aborting...",true);
      return false;
    }
  }

  // Load PU profiles for MC reweighing
  if (jp::ismc and jp::reweighPU) {
    if (!LoadPuProfiles()) {
      PrintInfo("Issues loading the PU histograms for reweighting; aborting...",true);
      return false;
    }
  }

  // load luminosity tables (prescales now stored in event)
  if (jp::isdt and jp::dolumi) {
    if (!LoadLumi()) {
      PrintInfo("Issues loading the Lumi file; aborting...",true);
      return false;
    }
  }

  if (jp::ismc) PrintInfo(Form("Running on MC produced with %1.3g nb-1 (%lld evts)",
                          1000. * _ntot / jp::xsecMinBias, _ntot),true);
  if (jp::isdt) PrintInfo(Form("Running on %lld events of data",_ntot),true);

  // Initialize histograms for different epochs and DQM selections
  if (jp::doBasicHistos) {
    InitBasic("Standard");
  }

  if (jp::doEtaHistos) {
    InitEta("FullEta_Reco");
    if (jp::ismc) {
      InitEta("FullEta_Gen");
      if (jp::doEtaHistosMcResponse) {
        InitMC("FullEta_RecoPerGen_vReco");
        InitMC("FullEta_RecoPerGen_vGen");
      }
    }
  }

  if (jp::isdt and jp::doMpfHistos) {
    InitAll("AllTrigs");
  }

  if (jp::isdt and jp::doRunHistos) {
    InitRun("Runs",0.,3.);
    InitRun("RunsBarrel",0.,1.);
    InitRun("RunsTransition",1.,2.);
    InitRun("RunsEndcap",2.,3.);
  }

  // For debugging purposes, save the weights used for pileup profiles.
  if (jp::ismc and jp::reweighPU) {
    _outfile->cd();
    _outfile->mkdir("puwgt");
    _outfile->cd("puwgt");
    for (auto &puprof : _pudist)
      puprof.second->Write();
  }

  if (jp::doVetoHot) {
    string HotTag = "";
    string HotYr = "";
    if (jp::yid==0) {
      HotYr = "16";
      if      (std::regex_search(jp::run,regex("^Run[BCD]"))) HotTag = "BCD";
      else if (std::regex_search(jp::run,regex("^RunE")))     HotTag = "EF";
      else if (std::regex_search(jp::run,regex("^RunFe")))    HotTag = "EF";
      else if (std::regex_search(jp::run,regex("^RunFl")))    HotTag = "GH";
      else if (std::regex_search(jp::run,regex("^Run[GH]")))  HotTag = "GH";
    } else if (jp::yid==1) {
      HotYr = "17";
      if      (std::regex_search(jp::run,regex("^RunB"))) HotTag = "B";
      else if (std::regex_search(jp::run,regex("^RunC"))) HotTag = "C";
      else if (std::regex_search(jp::run,regex("^RunD"))) HotTag = "D";
      else if (std::regex_search(jp::run,regex("^RunE"))) HotTag = "E";
      else if (std::regex_search(jp::run,regex("^RunF"))) HotTag = "F";
    }
    assert(HotTag!="");
    fHotExcl = new TFile(Form("rootfiles/hotjets-%srun%s.root",HotYr.c_str(),HotTag.c_str()),"READ");
    assert(fHotExcl and !fHotExcl->IsZombie() and Form("file rootfiles/hotjets-run%s.root missing",HotTag.c_str()));
    h2HotExcl = (TH2D*)fHotExcl->Get(Form("h2hot%s",jp::HotType));
    assert(h2HotExcl and "erroneous eta-phi exclusion type");
    PrintInfo(Form("Loading hot zone corrections rootfiles/hotjets-run%s.root with h2hot %s",HotTag.c_str(),jp::HotType));
  }
  // Qgl: load quark/gluon probability histos (Ozlem)
  // 1. open the previous output-MC-1_iteration1.root in the beginning
  // 2. get the hqgl_q and hqgl_g for each eta bin, store in array
  // 3. find correct hqgl_q and hqgl_g from array (normalized)
  // 4. calculate probg = g / (q+g)
  if (jp::doqglfile) {
    TFile *finmc = new TFile(jp::qglfile,"READ");
    assert(finmc && !finmc->IsZombie());

    double veta[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7};
    const int neta = sizeof(veta)/sizeof(veta[0])-1;
    // same or vpt, vtrigpt, npt
    const int npt = jp::notrigs;
    const double *vtrigpt = &jp::trigthr[0];
    double vpt[npt+1]; vpt[0] = 0;
    for (int trigidx = 0; trigidx != npt; ++trigidx) vpt[trigidx+1] = jp::trigranges[trigidx][1];
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
        PrintInfo(sg,true);
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
  if (jp::isdt) // For DT fetch true pileup from the json or histogram info
    trpu = _avgpu[run][lbn];
  else // For MC the info exists in the SMPJ tree
    trpu = EvtHdr__mTrPu;

  if (njt > _maxnjt) {
    PrintInfo(Form("Array overflow: njt = %d > njtmax= %d",njt,_maxnjt),true);
    assert(njt<_maxnjt);
  }

  if (jp::ismc and gen_njt > _maxnjt) {
    PrintInfo(Form("Array overflow: gen_njt = %d > njtmax= %d",gen_njt,_maxnjt),true);
    assert(gen_njt<_maxnjt);
  }

  if (jp::debug) {
    PrintInfo("",true);
    Show(_jentry);
    PrintInfo("",true);
    PrintInfo("",true);

    PrintInfo("***Checking basic event variables are read out:",true);
    PrintInfo(Form("isdata = %d / ismc = %d",jp::isdt,jp::ismc),true);
    PrintInfo(Form("trpu = %f",trpu),true);
    PrintInfo(Form("pthat = %f",pthat),true);
    PrintInfo(Form("weight = %f",weight),true);
    PrintInfo(Form("njt = %d",njt),true);
    PrintInfo(Form("idloose[0] = %d",njt>0 ? jtidloose[0] : -1),true);
    PrintInfo(Form("idtight[0] = %d",njt>0 ? jtidtight[0] : -1),true);
    PrintInfo("***end basic event variables");
    PrintInfo("",true);
  }

  // Check if duplicate
  if (jp::isdt and jp::checkduplicates) {
    set<int>& events = _duplicates[run][lbn];
    if (events.find(evt)!=events.end()) {
      ++_nbadevts_dup;
      return false;
    }
    events.insert(evt);
  }
  ++_cnt["01all"];

  // Check if good run/LS, including JSON selection
  if (jp::isdt) {
    if (jp::dojson) {
      // Does the run/LS pass the latest JSON selection?
      if (_json[run][lbn]==0) {
        _badjson.insert(pair<int, int>(run, lbn));
        ++_nbadevts_json;
        return false;
      }
    } // jp::dojson

    if (jp::dolumi) {
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
    } // jp::dolumi
  } // jp::isdt

  // Keep track of LBNs
  _passedlumis.insert(pair<int, int>(run, lbn));

  ++_cnt["02ls"];

  // Reset event ID
  _pass = true;


  if (jp::fetchMETFilters) {
    if (jp::isnu or FilterDecision_.size()==0) ++_cnt["03METFlt"];
    else {
      // If we perform MET filtering, any filter firing will cause the event to be discarded.
      if (jp::doMETFiltering) return false;
      _pass = false;

      // If we don't perform MET filtering, we can do something more elaborate with the filter info
      for (unsigned FDec : FilterDecision_) {
        assert(FDec<_availFlts.size());
        auto &FName = _availFlts[FDec];

        // TODO: Do something interesting?

        if (jp::debug) PrintInfo(Form("%s %u",FName.c_str(),FDec),true);
      } // for FDec (_FilterDecision)
    }
  } else {
    // This is equal to the SMPJ GoodVertexFilter and is redundant with MET filters on
    pvrho = tools::oplus(pvx, pvy);
    _pass = _pass and pvndof > 4 and fabs(pvz) <= 24. and pvrho <= 2. and npvgood>0.;
    if (_pass) ++_cnt["03vtx"];
  }

  // Event cuts against beam backgrounds
  if (_pass) {
    if (tools::oplus(pvx-bsx, pvy-bsy)>0.15) {
      ++_bscounter_bad;
      _pass = false;
    } else {
      ++_bscounter_good;
      ++_cnt["04bsc"];
    }
  }

  // load correct IOV for JEC
  if (jp::redojes and jp::isdt and jp::useIOV) {
    bool setcorrection = _iov.setJEC(&_JEC,&_L1RC,&_jecUnc,run);
    if (!setcorrection or !_JEC or !_L1RC or !_jecUnc) {
      PrintInfo("Issues while loading JEC; aborting...",true);
      assert(false);
    }
  }

  if (jp::debug) PrintInfo("JEC and MET calculation and leading jets info!",true);
  // Calculate jec and propagate jec to MET 1 and MET 2
  double mex = met * cos(metphi);
  double mey = met * sin(metphi);
  double mex_nol2l3 = mex;
  double mey_nol2l3 = mey;
  double ucx = mex;
  double ucy = mey;
  // Find leading jets (residual JEC may change ordering)
  // CAUTION: for less than 3 jets, we input -1 on the place of the index
  // NOTE: Above this we only use the initial corrections. This serves as a good approximation,
  // as the values before reapplying JEC are more like the ones we knew when collecting data.
  for (int i = 0; i<3; ++i) jt3leads[i] = -1;
  ///////////////
  // Jet loop
  ///////////////
  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    p4.SetPxPyPzE(jtp4x[jetidx],jtp4y[jetidx],jtp4z[jetidx],jtp4t[jetidx]);
    // Divide by the original JES
    if (jp::debug) PrintInfo("Entering jet loop!",true);
    if (jp::undojes) p4 *= 1/jtjes[jetidx];

    jtptu[jetidx] = p4.Pt();
    jteu[jetidx] = p4.E();

    double jec_res = 1;
    if (jp::redojes) {
      if (jp::debug) PrintInfo("Recalculating JEC!",true);
      // Recalculate JEC
      _JEC->setRho(rho);
      _JEC->setNPV(npvgood);
      _JEC->setJetA(jta[jetidx]);
      _JEC->setJetPt(jtptu[jetidx]);
      _JEC->setJetE(jteu[jetidx]);
      _JEC->setJetEta(p4.Eta());
      vector<float> v = _JEC->getSubCorrections();
      jtjesnew[jetidx] = v.back();

      if (jp::ismc or jp::skipl2l3res) {
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
    } else {
      jtjesnew[jetidx] = 1.;
      jtjes_l1[jetidx] = 1.;
      jtjes_l2l3[jetidx] = 1.;
      jtjes_res[jetidx] = 1.;
    }

    // Calculate gen level info
    if (jp::ismc) {
      if (jp::debug) PrintInfo("Gen info!",true);
#ifdef NEWMODE
      auto &gjetidx = jtgenidx[jetidx];
      if (gjetidx>=0 and gjetidx<GenJets__)
        gp4.SetPxPyPzE(gen_jtp4x[gjetidx],gen_jtp4y[gjetidx],gen_jtp4z[gjetidx],gen_jtp4t[gjetidx]);
      else
        gp4.SetPxPyPzE(0.,0.,0.,0.);
#else
      gp4.SetPxPyPzE(jtgenp4x[jetidx],jtgenp4y[jetidx],jtgenp4z[jetidx],jtgenp4t[jetidx]);
#endif
      jtgenpt[jetidx] = gp4.Pt();
      jtgeny[jetidx] = gp4.Rapidity();
      jtgeneta[jetidx] = gp4.Eta();
      jtgenphi[jetidx] = gp4.Phi();
    }

    if (jp::debug) PrintInfo("Reapplying JEC!",true);
    if (jp::redojes) p4 *= jtjesnew[jetidx];

#if JETRESO == 1
    // For DATA, we have redundant scale factor files, which return a "1" for SF.
    // We save the processing time by skipping this.
    if (jp::ismc) {
      double genpt = gp4.Pt();
      // We enforce the gen match, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
      if (genpt>0) {
        double genDR = p4.DeltaR(gp4);
        // We require the match to be closer than half cone, otherwise gen Pt is worthless
        if (genDR>0.2) genpt = -4;
      }
      p4 *= _iov.getJERCF(p4.Pt(),p4.Eta(),p4.E(),rho,genpt);
    }
#endif

    jte[jetidx] = p4.E();
    jtpt[jetidx] = p4.Pt();
    jteta[jetidx] = p4.Eta();
    jtphi[jetidx] = p4.Phi();
    jty[jetidx] = p4.Rapidity();
    if (jp::doMpfHistos) jtpt_nol2l3[jetidx] = p4.Pt()/jec_res;

    if (jp::debug) PrintInfo(Form("Jet %d corrected!",jetidx),true);

    // Only use jets with corr. pT>recopt GeV to equalize data and MC thresholds
    if (jtpt[jetidx] > jp::recopt and fabs(jteta[jetidx])<4.7) {
      // MET 1: the one where JEC is applied. MET1 needs to be recalculated as JEC changes.
      // Subtract uncorrected jet pT from met, put back corrected & add L1RC offset to keep PU isotropic.
      double l1corr = 1.;
      if (jp::redojes) {
        _L1RC->setRho(rho);
        _L1RC->setJetA(jta[jetidx]);
        _L1RC->setJetPt(jtptu[jetidx]);
        _L1RC->setJetE(jteu[jetidx]);
        _L1RC->setJetEta(jteta[jetidx]);
        l1corr = _L1RC->getCorrection();
      } else {
        l1corr = jtjes[jetidx];
      }
      double dpt = - jtpt[jetidx] + l1corr*jtptu[jetidx]; // old: + (l1chs - l1pf + l1corr)*jtptu[jetidx];

      mex += dpt * cos(jtphi[jetidx]);
      mey += dpt * sin(jtphi[jetidx]);
      if (jp::doMpfHistos) {
        double dpt_nol2l3 = - jtpt_nol2l3[jetidx] + l1corr*jtptu[jetidx];
        mex_nol2l3 += dpt_nol2l3 * cos(jtphi[jetidx]);
        mey_nol2l3 += dpt_nol2l3 * sin(jtphi[jetidx]);
      }

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
  met1 = tools::oplus(mex, mey); // met1 = -jtpt[jetidx]+l1corr*jtptu[jetidx]+metraw
  metphi1 = atan2(mey, mex);
  if (jp::doMpfHistos) {
    met1_nol2l3 = tools::oplus(mex_nol2l3,mey_nol2l3);
    metphi1_nol2l3 = atan2(mey_nol2l3, mex_nol2l3);
  }

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
  if (_pass) ++_cnt["05njt"];

  if (jp::debug) PrintInfo(Form("Indices for the three leading jets: %d %d %d\nGen flav calculation!",i0,i1,i2),true);

  // Check rho
  if (_pass) {
    if (rho>500.) {
    //if (rho>40.) { // This was for run 1 settings
      ++_rhocounter_bad;
      _pass = false;
      if (jp::debug) PrintInfo(Form("\nrun:ev:ls %d:%d:%lld : rho=%1.1f njt=%d npv=%d jtpt0=%1.1f sumet=%1.1f met=%1.1f\n",run, lbn, evt, rho, njt, npv, (njt>0 ? jtpt[i0] :0.),metsumet,met),true);
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
  if (jp::ismc) {
    // Always insert the generic mc trigger
    if (jp::debug) PrintInfo("Entering PU weight calculation!",true);
#ifdef NEWMODE
    if (_pass and (jtgenidx[i0]!=-1 or jp::isnu)) ++_cnt["07mcgenjet"];
    else return false;
#endif
    if (jp::domctrigsim and njt>0) {
      // Only add the greatest trigger present
      // Calculate trigger PU weight
      bool found = false;
      bool wcond = false;
      for (unsigned itrg = 0; itrg < jp::triggers.size(); ++itrg) {
        double wtrue = 1.0;
        // We fire all the triggers up to the unfeasible turn-on point
        if (jtpt[i0]>jp::trigranges.at(itrg).at(0)) {
          string trg_name = jp::triggers.at(itrg);
          _trigs.insert(trg_name);
          _wt[trg_name] = 1.;

          // Reweight in-time pile-up
          if (jp::reweighPU) {
            int k = _pudist[trg_name]->FindBin(trpu);
            wtrue = _pudist[trg_name]->GetBinContent(k);
            _wt[trg_name] *= wtrue;
            wcond |= wtrue!=0;
          }
          found = true;
        } else {
          break;
        }
      }
      _pass = _pass and found;
      if (_pass) ++_cnt["07mctrg"];
      else return false; // Leading jet is weak

      // check for non-zero PU weight
      _pass = _pass and wcond;
      if (_pass) ++_cnt["07puw"];
      else return false; // Bad pu areas with zero weight excluded
    } // jp::domctrigsim
    _trigs.insert("mc");
    _wt["mc"] = 1.0;
    if (jp::reweighPU) {
      int k = _pudist[jp::reftrig]->FindBin(trpu);
      _wt["mc"] *= _pudist[jp::reftrig]->GetBinContent(k);
    }
  } else if (jp::isdt) {
    // For data, check trigger bits
    if (jp::debug) {
      PrintInfo(Form("TriggerDecision_.size()==%zu",TriggerDecision_.size()));
      PrintInfo(Form("_availTrigs.size()==%zu",_availTrigs.size()));
      PrintInfo(Form("_goodTrigs.size()==%zu",_goodTrigs.size()));
    }
    #ifndef NEWMODE
    assert(TriggerDecision_.size()==_availTrigs.size());
    #endif

    // New and old mode: TriggerDecision and L1/HLT Prescales have the same indexing.
    for (auto itrg = 0u; itrg != jp::notrigs; ++itrg) _wt[jp::triggers[itrg]] = 1.0;

    #ifdef NEWMODE
    for (unsigned itrg = 0; itrg<TriggerDecision_.size(); ++itrg) {
      unsigned TDec = TriggerDecision_[itrg]; // Location of the current place
      assert(TDec<_availTrigs.size());
      auto &TName = _availTrigs[TDec];

      auto trgPlace = std::find(_goodTrigs.begin(),_goodTrigs.end(),TDec);
      if (trgPlace==_goodTrigs.end()) continue;
      unsigned goodIdx = static_cast<unsigned int>(trgPlace-_goodTrigs.begin());
    #else
    for (auto goodIdx = 0u; goodIdx < _goodTrigs.size(); ++goodIdx) {
      auto &itrg = _goodTrigs[goodIdx];
      auto &TDec = TriggerDecision_[itrg]; // Trigger fired or not: -1, 0, 1
      if (TDec!=1) continue;
      assert(itrg<_availTrigs.size());
      auto &TName = _availTrigs[itrg];
    #endif // NEWMODE

      #ifdef NEWMODE
      if (jp::debug)
      #else
      if (jp::debug and TDec>0)
      #endif
      PrintInfo(Form("%s %d %d %d %d",TName.c_str(),itrg,TDec,L1Prescale_[itrg],HLTPrescale_[itrg]),true);

      // Set prescale from event for now
      //if (L1Prescale_[itrg]>0 and HLTPrescale_[itrg]>0) { There's trouble in 2017 L1, so we let it pass
      if (HLTPrescale_[itrg]>0 or L1Prescale_[itrg]>0) {
        double l1 = L1Prescale_[itrg];
        double hlt = HLTPrescale_[itrg];
        if (l1==0) l1 = 1;
        if (hlt==0) hlt = 1;
        _prescales[TName][run] = l1 * hlt;
      } else {
        PrintInfo(Form("Error for trigger %s prescales: L1 = %d HLT = %d",TName.c_str(),L1Prescale_[itrg],HLTPrescale_[itrg]));
        _prescales[TName][run] = 0;
        if (jp::debug) { // check prescale
          double prescale = _prescales[TName][run];
          if (L1Prescale_[itrg]*HLTPrescale_[itrg]!=prescale) {
            PrintInfo(Form("Trigger %s, Prescale (txt file) = %f",TName.c_str(),prescale),true);
            PrintInfo(Form("L1 = %d, HLT = %d",L1Prescale_[itrg],HLTPrescale_[itrg]),true);
            assert(false);
          }
        } // debug
      }

      if (_prescales[TName][run]>0) {
        // Set trigger only if prescale information is known
        _trigs.insert(TName);
        _wt[TName] = _goodWgts[goodIdx];
      } else {
        // Make sure all info is good! This is crucial if there is something odd with the tuples
        PrintInfo(Form("Missing prescale for %s in run %d",TName.c_str(),run),true);
      }
    } // for itrg (FilterDecision or _goodTrigs)
  } // if isdt

  ++_totcounter;
  if (_pass) ++_evtcounter;
  _pass = _pass and _trigs.size()>0;
  if (_pass) {
    ++_trgcounter;
    if (jp::isdt) ++_cnt["08trg"];
  }

  // Retrieve event weight. _w0 is static, _w is chanching with the trigger
  _w0 = 1.0;
  if (jp::ismc) {
    _w0 *= weight;
    if (jp::pthatbins or jp::htbins) _w0 *= _binnedmcweight;
    assert(_w0>0);
  }
  _w = _w0;

  // TODO: implement reweighing for k-factor (NLO*NP/LOMC)

  if (jp::ismc) {
    ///////////////
    // Gen Jet loop
    ///////////////
    bool doht = jp::isnu or (jp::htbins and !jp::pthatbins);
    double htsum = 0.0;
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      genp4.SetPxPyPzE(gen_jtp4x[gjetidx],gen_jtp4y[gjetidx],gen_jtp4z[gjetidx],gen_jtp4t[gjetidx]);
      gen_jtpt[gjetidx] = genp4.Pt();
      gen_jteta[gjetidx] = genp4.Eta(); // for matching
      gen_jtphi[gjetidx] = genp4.Phi(); // for matching
      gen_jty[gjetidx] = genp4.Rapidity();
      if (doht) htsum += genp4.Pt();

      // Ozlem: loop for finding partonflavor by matching genjets and jets
      int ireco = -1;
      for (int j = 0; j != njt && ireco == -1; ++j) {
        p4.SetPxPyPzE(jtp4x[j],jtp4y[j],jtp4z[j],jtp4t[j]);
        if (genp4.DeltaR(p4) < 0.4)
          ireco = j;
      } // for j
      if (ireco!=-1)
#ifdef NEWMODE
        gen_partonflavor[gjetidx] = Int_t(partonflavorphys[ireco]+0.25); //qglokan
#else
        gen_partonflavor[gjetidx] = Int_t(partonflavor[ireco]+0.25); //qglokan
#endif
      else
        gen_partonflavor[gjetidx] = -1;
    } // for gjetidx
    // In NuGun samples, we don't have gen jets
    for (int jetidx = 0; jetidx != njt; ++jetidx) htsum += jtpt[jetidx];
    if (doht) pthat = htsum/2.0;

    // Check if overweight PU event
    if (_pass) {
      if (!jp::isnu) {
        if (jtpt[i0] < 1.5*jtgenpt[i0] or jp::isnu) ++_cnt["09ptgenlim"];
        else _pass = false;
      }

      if (_pass) {
        if (doht) {
          if (jtpt[i0] < 2.0*pthat) ++_cnt["10htlim"];
          else _pass = false;
        } else {
          double lim = (pthat < 100) ? 2.0 : 1.5;
          if (jtpt[i0] < lim*pthat) ++_cnt["10pthatlim"];
          else _pass = false;
        }
      }
    }
  } // MC

  _jetids.resize(njt);
  for (unsigned int jetid = 0; jetid != _jetids.size(); ++jetid) _jetids[jetid] = true;
  FillJetID(_jetids);

  if (_pass and _jetids[i0]) ++_cnt["11jtid"]; // Non-restrictive
#ifdef NEWMODE
  // Equipped in FillBasic and FillRun
  _pass_qcdmet = met01 < 45. or met01 < 0.3 * metsumet01; // updated 4/2018
#else
  _pass_qcdmet = met < 45. or met < 0.3 * metsumet;
#endif

  return true;
}


// Report event stuff
void HistosFill::Report()
{
  // List bad runs
  PrintInfo(Form("Processed %d events in total",_totcounter));
  PrintInfo(Form("Processed %d events passing basic data quality and trigger cuts",_trgcounter));
  PrintInfo(Form("(out of %d passing data quality cuts)",_evtcounter));
  if (jp::dojson)
    PrintInfo(Form("Found %d bad events according to new JSON (events cut)",_nbadevts_json));
  if (jp::dolumi) {
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
  if (jp::checkduplicates)
    PrintInfo(Form("Found %d duplicate events, which were properly discarded",_nbadevts_dup));

  // Report beam spot cut efficiency
  PrintInfo(Form("Beam spot counter discarded %d events out of %d (%f %%)\nBeam spot expectation is less than 0.5%%",
                 _bscounter_bad,_bscounter_good,double(_bscounter_bad)/double(_bscounter_good)*100.));

  // Report rho veto efficiency
  PrintInfo(Form("Rho<40 veto counter discarded %d events out of %d (%f %%)\nRho veto expectation is less than 1 ppm",
                _rhocounter_bad,_rhocounter_good,double(_rhocounter_bad)/double(_rhocounter_good)*100.));
  for (auto &cit : _cnt)
    PrintInfo(Form("%s: %d (%1.1f%%)",cit.first.c_str(),cit.second,100.*cit.second/max(1,_cnt["01all"])),true);

  if (jp::isdt) PrintInfo("Note that for DT it is likely that we lose a large percentage of events for the trigger.\n"
                          "Events triggered by JetHT and AK8PFJet are included in addition to AK4PFJet.");
  *ferr << endl;

  PrintInfo("Reporting lumis not discarded in reports/passedlumis.json",true);
  ofstream fout(Form("reports/passedlumis-%s.json",jp::type), ios::out);
  for (auto &lumit : _passedlumis) fout << lumit.first << " " << lumit.second << endl;

  if (jp::dojson and _badjson.size()>0) {
    PrintInfo("Reporting lumis discarded by json selection (jp::dojson) in reports/badlumis_json",true);
    ofstream fout2(Form("reports/badlumis_json-%s.json",jp::type), ios::out);
    for (auto &jsit : _badjson) fout2 << jsit.first << " " << jsit.second << endl;
  } // jp::dojson

  if (jp::dolumi) {
    if (_badlums.size()>0) {
      PrintInfo("Reporting lumis discarded by lumifile selection (jp::dolumi) in reports/badlumis_lumi.json",true);
      ofstream fout2(Form("reports/badlumis_lumi-%s.json",jp::type), ios::out);
      for (auto &jsit : _nolums) fout2 << jsit.first << " " << jsit.second << endl;
    } // _badlums
    if (_nolums.size()>0) {
      PrintInfo("Reporting lumis discarded fby lumifile selection (jp::dolumi) @ zero luminosity in reports/badlumis_zerolumi.json",true);
      ofstream fout2(Form("reports/badlumis_zerolumi-%s.json",jp::type), ios::out);
      for (auto &jsit : _nolums) fout2 << jsit.first << " " << jsit.second << endl;
    } // _nolums
  } // jp::dolumi
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
  if (jp::ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(jp::recopt, jp::emax);
    pttrg["mc"] = jp::recopt;
  }
  if (jp::isdt or jp::domctrigsim) {
    // This is done both for data and MC, because why not?
    for (unsigned itrg = 0; itrg != jp::notrigs; ++itrg) {
      string trg = jp::triggers[itrg];
      triggers.push_back(trg);
      double pt1 = jp::trigranges[itrg][0];
      double pt2 = jp::trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = jp::trigthr[itrg];
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
  if (!fired) return;

  // Luminosity information
  if (jp::isdt and h->lums[run][lbn]==0) {
    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    double prescale(0);
    auto ip = _prescales[h->trigname].find(run);
    if (ip==_prescales[h->trigname].end()) {
      PrintInfo(Form("No prescale info for trigger %s in run %d!",h->trigname.c_str(),run));
      assert(false);
    } else {
      prescale = ip->second;
      if (prescale==0) {
        PrintInfo(Form("Prescale zero for trigger %s in run %d!",h->trigname.c_str(),run));
        prescale = 1.;
        assert(false);
      }
    }

    h->lumsum += lum / prescale;
    h->lumsum2 += lum2 / prescale;
    h->lums[run][lbn] = 1;

    h->hlumi_vstrpu->Fill(trpu, lum / prescale);
  }

  // For MC vs truePU
  if (jp::ismc)
    h->hlumi_vstrpu->Fill(trpu, _w);

  if (jp::debug) {
    if (h == _histos.begin()->second[0]) {
      PrintInfo(Form("Triggers size: %lu",_trigs.size()),true);
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit) cout << *trgit << ", ";
      PrintInfo(Form("Current: (%s)",h->trigname.c_str()),true);
    }
  }

  // Are all of our conditions met?
  if (!_pass) return;

  if (jp::debug) PrintInfo(Form("Subdirectory Eta_%1.1f-%1.1f/%s",h->etamin,h->etamax,h->trigname.c_str()),true);
  if (jp::debug) PrintInfo("Calculate and fill dijet mass",true);

  if (h->ismcdir) {
    h->hpthat->Fill(pthat, _w);
    h->hpthatnlo->Fill(pthat);
  }

  //{ Pre-calculate some nice garden tools
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (i0 < 0.) return; // This should not happen, but check just in case

  double ptave = (i1>=0 ? 0.5 * (jtpt[i0] + jtpt[i1]) : jtpt[i0]);
  double dphi = (i1>=0 ? DPhi(jtphi[i0], jtphi[i1]) : 0.);
  double dpt = (i1>=0 ? fabs(jtpt[i0]-jtpt[i1])/(2*ptave) : 0.999);
  // If the jetID is bad for the third jet (and the third jet is visible), we set pt3 to ptave (alpha = 1)
  double pt3 = ((i1>=0 and i2>=0 and jtpt[i2]>jp::recopt) ? (_jetids[i2] ? jtpt[i2] : ptave) : 0.);
  double alpha = pt3/ptave;
  //} Garden tools

  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>jp::recopt) { // First leading jet
    if (i1>=0 and _jetids[i1] and jtpt[i1]>jp::recopt) { // Second leading jet
      //{ Calculate and fill dijet mass.
      _j1.SetPtEtaPhiE(jtpt[i0],jteta[i0],jtphi[i0],jte[i0]);
      _j2.SetPtEtaPhiE(jtpt[i1],jteta[i1],jtphi[i1],jte[i1]);
      double djmass = (_j1+_j2).M();
      double etamaxdj = max(fabs(jteta[i0]),fabs(jteta[i1]));
      bool goodjets = (jtpt[i0]>30. and jtpt[i1]>30.);
      // The eta sectors are filled according to max eta
      if (goodjets and etamaxdj >= h->etamin and etamaxdj < h->etamax) {
        assert(h->hdjmass); h->hdjmass->Fill(djmass, _w);
        assert(h->hdjmass0); h->hdjmass0->Fill(djmass, _w);
        assert(h->pdjmass_ptratio); h->pdjmass_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
        assert(h->pdjmass0_ptratio); h->pdjmass0_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
        if (dphi > 2.7) { // Back-to-back condition
          if (alpha<0.1) { assert(h->hdjmass_a01); h->hdjmass_a01->Fill(djmass, _w); }
          if (alpha<0.2) { assert(h->hdjmass_a02); h->hdjmass_a02->Fill(djmass, _w); }
          if (alpha<0.3) { assert(h->hdjmass_a03); h->hdjmass_a03->Fill(djmass, _w); }
        }
      }
      //} Dijet mass
      //{ Calculate and fill jet mass.
      assert(h->hjmass);  h->hjmass->Fill(_j1.M(),weight);  h->hjmass->Fill(_j2.M(),weight);
      assert(h->hjmass0); h->hjmass0->Fill(_j1.M(),weight); h->hjmass0->Fill(_j2.M(),weight);
      if (dphi > 2.7) { // Back-to-back condition
        if (alpha<0.1) { assert(h->hjmass_a01); h->hjmass_a01->Fill(_j1.M(),weight); h->hjmass_a01->Fill(_j2.M(),weight); }
        if (alpha<0.2) { assert(h->hjmass_a02); h->hjmass_a02->Fill(_j1.M(),weight); h->hjmass_a02->Fill(_j2.M(),weight); }
        if (alpha<0.3) { assert(h->hjmass_a03); h->hjmass_a03->Fill(_j1.M(),weight); h->hjmass_a03->Fill(_j2.M(),weight); }
      }
      //}

      //{ Tag & probe hoods: Tag in barrel and fires trigger, probe in eta bin unbiased
      if (jp::debug) PrintInfo("Calculate and fill dijet balance",true);

      if (dphi > 2.7) { // Back-to-back condition
        for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
          int itag = jt3leads[itag_lead];
          int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
          double etatag = jteta[itag];
          double etaprobe = jteta[iprobe];

          // Eta sel: tag in barrel, probe in current slice
          if (fabs(etatag) < 1.3 and fabs(etaprobe) >= h->etamin and fabs(etaprobe) < h->etamax) { // Eta sel.
            double pttag = jtpt[itag];
            double ptprobe = jtpt[iprobe];
            double phiprobe = jtphi[iprobe];

            //{ Dijet balance histograms
            if (jp::do3dHistos) {
              double asymm = (ptprobe - pttag)/(2*ptave);
              double asymmtp = (ptprobe - pttag)/(2*pttag);
              double mpf = met1*cos(DPhi(metphi1,jtphi[itag]))/2;
              double mpftp = mpf/pttag;
              mpf /= ptave;

              // Asymmetry and mpf
              assert(h->hdjasymm);   h->hdjasymm->Fill(ptave, alpha, asymm, _w);
              assert(h->hdjasymmtp); h->hdjasymmtp->Fill(pttag, alpha, asymmtp, _w);
              assert(h->hdjmpf);     h->hdjmpf->Fill(ptave, alpha, mpf, _w);
              assert(h->hdjmpftp);   h->hdjmpftp->Fill(pttag, alpha, mpftp, _w);
            }
            assert(h->ppttagptprobe_noa); h->ppttagptprobe_noa->Fill(pttag,ptprobe, _w);
            assert(h->h2pttagptprobe_noa); h->h2pttagptprobe_noa->Fill(pttag,ptprobe, _w);
            assert(h->ppttagmu_noa); h->ppttagmu_noa->Fill(pttag,trpu, _w);
            assert(h->h2pttagmu_noa); h->h2pttagmu_noa->Fill(pttag,trpu, _w);

            //} // Dijet balance
            if (alpha >= 0.3) {
              assert(h->ppttageff); h->ppttageff->Fill(pttag, 0, _w);
            } else {
              assert(h->ppttagptprobe); h->ppttagptprobe->Fill(pttag,ptprobe, _w);
              assert(h->h2pttagptprobe); h->h2pttagptprobe->Fill(pttag,ptprobe, _w);
              assert(h->ppttageff); h->ppttageff->Fill(pttag, 1, _w);
              assert(h->ppttagmu); h->ppttagmu->Fill(pttag,trpu, _w);
              assert(h->h2pttagmu); h->h2pttagmu->Fill(pttag,trpu, _w);
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
              assert(h->ppuftp); h->ppuftp->Fill(pttag, jtbetaprime[iprobe]*jtchf[iprobe], _w);

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
                assert(h->hpuftp);      h->hpuftp->Fill(jtbetaprime[iprobe]*jtchf[iprobe], _w);

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
                assert(h->ppuftp_vsnpv);      h->ppuftp_vsnpv->Fill(npvgood, jtbetaprime[iprobe]*jtchf[iprobe], _w);

                // Fractions vs true pileup
                assert(h->pchftp_vstrpu);      h->pchftp_vstrpu->Fill(trpu, jtchf[iprobe], _w);
                assert(h->pneftp_vstrpu);      h->pneftp_vstrpu->Fill(trpu, (jtnef[iprobe]-jthef[iprobe]), _w);
                assert(h->pnhftp_vstrpu);      h->pnhftp_vstrpu->Fill(trpu, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                assert(h->pceftp_vstrpu);      h->pceftp_vstrpu->Fill(trpu, jtcef[iprobe], _w);
                assert(h->pmuftp_vstrpu);      h->pmuftp_vstrpu->Fill(trpu, jtmuf[iprobe], _w);
                assert(h->phhftp_vstrpu);      h->phhftp_vstrpu->Fill(trpu, jthhf[iprobe], _w);
                assert(h->pheftp_vstrpu);      h->pheftp_vstrpu->Fill(trpu, jthef[iprobe], _w);
                assert(h->ppuftp_vstrpu);      h->ppuftp_vstrpu->Fill(trpu, jtbetaprime[iprobe]*jtchf[iprobe], _w);

                if (jp::doPhiHistos) {
                  if (etaprobe>0) {
                    assert(h->pchfpostp_vsphi);      h->pchfpostp_vsphi->Fill(phiprobe, jtchf[iprobe], _w);
                    assert(h->pnefpostp_vsphi);      h->pnefpostp_vsphi->Fill(phiprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
                    assert(h->pnhfpostp_vsphi);      h->pnhfpostp_vsphi->Fill(phiprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                    assert(h->pcefpostp_vsphi);      h->pcefpostp_vsphi->Fill(phiprobe, jtcef[iprobe], _w);
                    assert(h->pmufpostp_vsphi);      h->pmufpostp_vsphi->Fill(phiprobe, jtmuf[iprobe], _w);
                    assert(h->phhfpostp_vsphi);      h->phhfpostp_vsphi->Fill(phiprobe, jthhf[iprobe], _w);
                    assert(h->phefpostp_vsphi);      h->phefpostp_vsphi->Fill(phiprobe, jthef[iprobe], _w);
                    assert(h->ppufpostp_vsphi);      h->ppufpostp_vsphi->Fill(phiprobe, jtbetaprime[iprobe]*jtchf[iprobe], _w);
                  } else {
                    assert(h->pchfnegtp_vsphi);      h->pchfnegtp_vsphi->Fill(phiprobe, jtchf[iprobe], _w);
                    assert(h->pnefnegtp_vsphi);      h->pnefnegtp_vsphi->Fill(phiprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
                    assert(h->pnhfnegtp_vsphi);      h->pnhfnegtp_vsphi->Fill(phiprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
                    assert(h->pcefnegtp_vsphi);      h->pcefnegtp_vsphi->Fill(phiprobe, jtcef[iprobe], _w);
                    assert(h->pmufnegtp_vsphi);      h->pmufnegtp_vsphi->Fill(phiprobe, jtmuf[iprobe], _w);
                    assert(h->phhfnegtp_vsphi);      h->phhfnegtp_vsphi->Fill(phiprobe, jthhf[iprobe], _w);
                    assert(h->phefnegtp_vsphi);      h->phefnegtp_vsphi->Fill(phiprobe, jthef[iprobe], _w);
                    assert(h->ppufnegtp_vsphi);      h->ppufnegtp_vsphi->Fill(phiprobe, jtbetaprime[iprobe]*jtchf[iprobe], _w);
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
      bool has2 = (i1>=0 and jtpt[i1] > jp::recopt and fabs(jteta[i1])>=h->etamin and fabs(jteta[i1])<h->etamax);
      bool has3 = (i2>=0 and jtpt[i2] > jp::recopt and fabs(jteta[i2])>=h->etamin and fabs(jteta[i2])<h->etamax and jtpt[i1] > 0.70 * jtpt[i0]);
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

  if (jp::debug) PrintInfo("Entering jet loop",true);
  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    if (jp::debug) PrintInfo(Form("Loop over jet %d/%d",jetidx,njt));

    // adapt variable names from different trees
    double pt = jtpt[jetidx];
    double eta = jteta[jetidx];
    double energy = jte[jetidx];
    double mass = sqrt(fabs(pow(energy,2) - pow(pt*cosh(eta),2)));
    double y = jty[jetidx];
    double phi = jtphi[jetidx];
    double jec = jtjesnew[jetidx];
    bool pass_id = _jetids[jetidx];

    double jec2 = jtjesnew[jetidx]/jtjes[jetidx];
    bool etarange = fabs(eta) >= h->etamin and fabs(eta) < h->etamax;
    double ptgen = jtgenpt[jetidx];

    // Check effect of ID cuts
    if (etarange) { // Jet in eta range
      if (jp::debug) {
        PrintInfo(Form("... %s |  index %d/%d jet pt: %f eta : %f id: %d jec: %f",h->trigname.c_str(),jetidx,njt,pt,eta,pass_id,jec),true);
        PrintInfo(Form("...evt id: %d weight: %f met: %f metsumet: %f",_pass_qcdmet,_w,met,metsumet),true);
      }

      assert(h->hpt_noid); h->hpt_noid->Fill(pt, _w);
      assert(h->hpt_nojetid); if (_pass_qcdmet) h->hpt_nojetid->Fill(pt, _w);
      assert(h->hpt_noevtid); if (pass_id) h->hpt_noevtid->Fill(pt, _w);
      // Same versus generator pT as MC extra
      // to decouple efficiency from JEC and JER
      if (h->ismcdir) {
        h->hpt_noid_g->Fill(ptgen, _w);
        if (_pass_qcdmet) h->hpt_nojetid_g->Fill(ptgen, _w);
        if (pass_id)    h->hpt_noevtid_g->Fill(ptgen, _w);
      }
    } // Jet in eta range

    if (pt>jp::recopt) { // pt visible
      // Flags for studying gen eta vs reco eta effects
      bool mcgendr = jtgenr[jetidx] < 0.25;
      bool mcgenetarange = fabs(jtgeneta[jetidx]) >= h->etamin && fabs(jtgeneta[jetidx]) < h->etamax;

      if (h->ismcdir and mcgendr and mcgenetarange) // GenJets matched to any reco jets in any events
        h->hpt_gg0->Fill(ptgen, _w);

      if (_pass_qcdmet and pass_id) { // id OK
        if (h->ismcdir and mcgendr and mcgenetarange) // GenJets matched to good reco jets in good events
          h->hpt_gg->Fill(ptgen, _w);

        if (etarange) { // Correct jet eta range
          if (jp::debug) PrintInfo("..unfolding",true);
          if (h->ismcdir) {
            // unfolding studies (Mikael)
            h->my->Fill(pt, _w);
            h->myf->Fill(pt, _w);
            h->myuw->Fill(pt);
            h->myfuw->Fill(pt);

            if (mcgendr) {
              if (ptgen>jp::recopt) {
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

          if (jp::debug) PrintInfo("..jec uncertainty",true);
          // Get JEC uncertainty
          double unc = 0.01; // default for MC
          if (jp::isdt and _jecUnc) {
            _jecUnc->setJetEta(eta);
            _jecUnc->setJetPt(pt);
            unc = _jecUnc->getUncertainty(true);
            //_jecUnc2->Rjet(pt, unc); // use Fall10 absolute scale uncertainty
          }

          if (jp::dotrigeffsimple) {
            // calculate and/or retrieve efficiencies
            // REMOVED: "For trigger efficiency" (only dummies left)
            double ideff = 1.;
            double vtxeff = 1.;
            double dqmeff = 1.;
            double trigeff = 1.;
            double eff = ideff * vtxeff * dqmeff * trigeff;

            if (jp::debug) PrintInfo("..raw spectrum",true);

            // efficiencies
            assert(h->peff); h->peff->Fill(pt, eff, _w);
            assert(h->pideff); h->pideff->Fill(pt, ideff, _w);
            assert(h->pvtxeff); h->pvtxeff->Fill(pt, vtxeff, _w);
            assert(h->pdqmeff); h->pdqmeff->Fill(pt, dqmeff, _w);
          }

          // raw spectrum
          assert(h->hpt); h->hpt->Fill(pt,_w);
          assert(h->hpt0); h->hpt0->Fill(pt, _w); // Finer binning
          // leading jets (1,2,3)
          if (jetidx==i0) {
            assert(h->hpt1); h->hpt1->Fill(pt, _w);
          } else if (jetidx==i1) {
            assert(h->hpt2); h->hpt2->Fill(pt,_w);
          } else if (jetidx==i2) {
            assert(h->hpt3); h->hpt3->Fill(pt,_w);
          }
          assert(h->hpt_tmp); h->hpt_tmp->Fill(pt); // Event statistics
          assert(h->hpt_pre);
          if (jp::isdt) h->hpt_pre->Fill(pt, _w*_prescales[h->trigname][run] / _wt[h->trigname]);
          if (jp::ismc) h->hpt_pre->Fill(pt, _w0*_wt["mc"]);
          // Do proper event statistics
          assert(h->hpttmp); assert(h->hptevt);
          if (h->hpttmp->GetBinContent(h->hpttmp->FindBin(pt))==0)
            h->hptevt->Fill(pt, _w);
          h->hpttmp->Fill(pt);

          if (jp::debug) PrintInfo("..basic properties",true);

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

          if (jp::debug) PrintInfo("..control plots of components",true);

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
          assert(h->ppuf); h->ppuf->Fill(pt, jtbetaprime[jetidx]*jtchf[jetidx], _w);

          // control plots for topology (JEC)
          h->pa->Fill(pt, jta[jetidx], _w);
          h->ptrpu->Fill(pt, trpu, _w);
          h->prho->Fill(pt, rho, _w);
          h->pnpv->Fill(pt, npvgood, _w);
          h->pnpvall->Fill(pt, npv, _w);
          if (pt >= h->ptmin and pt < h->ptmax) { // Trigger pt range
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

            if (jp::debug) PrintInfo("..control plots for topology",true);

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
            h->hpuf->Fill(jtbetaprime[jetidx]*jtchf[jetidx], _w);

            h->hyeta->Fill(TMath::Sign(y-eta,y), _w);
            h->hyeta2->Fill(y-eta, _w);
            h->hetaphi->Fill(eta, phi, _w);
          } // within trigger pT range

          // closure plots for JEC
          h->pdpt->Fill(pt, dpt, _w);
          h->pmpf->Fill(pt, 1 + met * cos(DPhi(metphi, phi)) / pt, _w);
          h->pmpf1->Fill(pt, 1 + met1 * cos(DPhi(metphi1, phi)) / pt, _w);
          h->pmpf2->Fill(pt, 1 + met2 * cos(DPhi(metphi2, phi)) / pt, _w);

          // Histograms for quark/gluon study (Ozlem)
          double probg = 1. - qgl[jetidx]; // First approximation
          if (jp::doqglfile) { // If we loaded a previous file to _h3probg, use this for a better probg value
            assert(_h3probg);
            probg = _h3probg->GetBinContent(_h3probg->FindBin(eta,pt,qgl[jetidx]));
          }
          if (probg>=0 and probg<=1) {
            assert(h->hgpt);  h->hgpt->Fill(pt,_w*probg);
            assert(h->hgpt0); h->hgpt0->Fill(pt, _w*probg);
            assert(h->hqgl);  h->hqgl->Fill(qgl[jetidx], _w);
            assert(h->hqgl2); h->hqgl2->Fill(pt, qgl[jetidx], _w);
            if (jp::ismc) {
              assert(h->hqgl_g);
              assert(h->hqgl_q);
#ifdef NEWMODE
              bool isgluon = (fabs(partonflavorphys[jetidx]-21)<0.5); //qglokan
              bool isquark = (fabs(partonflavorphys[jetidx])<6 and partonflavorphys[jetidx]!=0); //qglokan why we are ingnoring t(6) b'(7) t'(8)
              bool isunmatch=((fabs(partonflavorphys[jetidx])>=6 and partonflavorphys[jetidx]!=21) || partonflavorphys[jetidx]==0);
#else
              bool isgluon = (fabs(partonflavor[jetidx]-21)<0.5); //qglokan
              bool isquark = (fabs(partonflavor[jetidx])<6 and partonflavor[jetidx]!=0); //qglokan why we are ingnoring t(6) b'(7) t'(8)
              bool isunmatch=(fabs(partonflavor[jetidx])<7); //qglokan
#endif
              assert(isgluon or isquark or isunmatch);

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
              } else if(isunmatch){
                h->hqgl_u->Fill(x, _w);
                h->hqgl2_u->Fill(pt, x, _w);
                double wg = 1;//-55.7067*pow(x,7) + 113.218*pow(x,6) -21.1421*pow(x,5) -99.927*pow(x,4) + 92.8668*pow(x,3) -34.3663*x*x + 6.27*x + 0.612992;
                assert(wg>0);
                h->hqgl_du->Fill(x, _w*wg);
                h->hqgl2_du->Fill(pt, x, _w*wg);
              }else {
                PrintInfo("Quark/Gluon status missing from partonflavor");
              }
            }
          } // probg quark/gluon

          if (h->ismcdir and mcgendr) { // MC extras
            if (jp::debug) PrintInfo(Form("genmatch %d ptg=%f yg=%f yr=%f",jetidx,ptgen,jtgeny[jetidx],y),true);

            double r = (ptgen ? pt/ptgen : 0);
            //double resp = (jtjesnew[jetidx] ? r / jtjesnew[jetidx] : 0);
            double dy = (r ? TMath::Sign(jty[jetidx]-jtgeny[jetidx], jtgeny[jetidx]) : 0.);
            h->hpt_r->Fill(pt, _w);
            h->hpt_g->Fill(ptgen, _w);

            h->hpt_gtw->Fill(ptgen, _w);
            // Ozlem: (gluon vs quark)
#ifdef NEWMODE
            if (partonflavorphys[jetidx]==21) h->hgpt_g->Fill(ptgen, _w); //qglokan we need unmatched ones
            else if (abs(partonflavorphys[jetidx])<6 and partonflavorphys[jetidx]!=0) h->hqpt_g->Fill(ptgen, _w);
            else h->hupt_g->Fill(ptgen, _w);
#else
            if (partonflavor[jetidx]==21) h->hgpt_g->Fill(ptgen, _w); //qglokan
            else if (abs(partonflavor[jetidx])<6 and partonflavor[jetidx]!=0) h->hqpt_g->Fill(ptgen, _w);
            else h->hupt_g->Fill(ptgen, _w);
#endif

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

  if (jp::ismc) {
    if (jp::debug) PrintInfo("Truth loop:",true);
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) { // Unbiased gen spectrum (for each trigger)
      double etagen = gen_jteta[gjetidx];
      double ptgen = gen_jtpt[gjetidx];

      if (fabs(etagen) >= h->etamin and fabs(etagen) < h->etamax) {
        if (jp::debug) PrintInfo(Form("genjet %d/%d ptg=%f etag=%f",gjetidx,gen_njt,ptgen,etagen),true);

        h->hpt_g0tw->Fill(ptgen, _w);
        // Ozlem: (gluon vs quark)
        if (gen_partonflavor[gjetidx]==21) h->hgpt_g0tw->Fill(ptgen, _w); //qglokan
        else if(abs(gen_partonflavor[gjetidx])<6 and gen_partonflavor[gjetidx]!=0) h->hqpt_g0tw->Fill(ptgen, _w);
        else h->hupt_g0tw->Fill(ptgen, _w);

        // unfolding studies (Mikael)
        if (h->ismcdir) { // Only the 'mc' trigger
          if (gen_jtpt[gjetidx]>jp::recopt) {
            h->mx->Fill(ptgen, _w);
            h->mxf->Fill(ptgen, _w);
            h->mxuw->Fill(ptgen);
            h->mxfuw->Fill(ptgen);
          } // recopt limit

          h->hpt_g0->Fill(gen_jtpt[gjetidx], _w);
          // Ozlem: (gluon vs quark)
          if (gen_partonflavor[gjetidx]==21) h->hgpt_g0->Fill(gen_jtpt[gjetidx], _w); //qglokan
          else if(abs(gen_partonflavor[gjetidx])<6 and gen_partonflavor[gjetidx]!=0) h->hqpt_g0->Fill(gen_jtpt[gjetidx], _w);
          else h->hupt_g0->Fill(gen_jtpt[gjetidx], _w);

          assert(h->hpt_g0_tmp); h->hpt_g0_tmp->Fill(gen_jtpt[gjetidx]);
        } // mcdir (a subset of jp::ismc)
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

  for (auto &histrange : _histos) {
    for (auto &h : histrange.second) {
      // Luminosity information
      for (int j = 0; j != h->hlumi->GetNbinsX()+1; ++j) {
        h->hlumi->SetBinContent(j, jp::isdt ? h->lumsum : 1. );
        h->hlumi2->SetBinContent(j, jp::isdt ? h->lumsum2 : 1. );
      }
      h->Write();
      delete h;
    } // for histidx
  } // for histit

  PrintInfo(Form("\nOutput (HistosBasic) stored in %s",_outfile->GetName()),true);

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
  if (jp::ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(jp::recopt, jp::emax);
    pttrg["mc"] = jp::recopt;
  }
  if (jp::isdt or jp::domctrigsim) {
    // This is done both for data and MC, because why not?
    for (unsigned itrg = 0; itrg != jp::notrigs; ++itrg) {
      string trg = jp::triggers[itrg];
      triggers.push_back(trg);
      double pt1 = jp::trigranges[itrg][0];
      double pt2 = jp::trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = jp::trigthr[itrg];
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

  if (jp::debug) {
    if (h == _etahistos.begin()->second[0]) {
      PrintInfo(Form("Triggers size: %lu",_trigs.size()),true);
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit) cout << *trgit << ", ";
      PrintInfo(Form("(%s)",h->trigname.c_str()),true);
    }
  }

  // check if required trigger fired
  if (!fired or !_pass) return;

  // Tag & probe hoods
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>jp::recopt and i1>=0 and _jetids[i1] and jtpt[i1]>jp::recopt) { // Quality conditios for leading jets
    double dphi = DPhi(jtphi[i0], jtphi[i1]);
    if (dphi > 2.7) { // Back-to-back condition
      double pt3 = ((i2>=0 and jtpt[i2]>jp::recopt) ? jtpt[i2] : 0.);
      double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
      double alpha = pt3/ptave;

      for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
        int itag = jt3leads[itag_lead];
        int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
        double etatag = jteta[itag];
        double etaprobe = jteta[iprobe];
        if (fabs(etatag) < 1.3) { // Tag required to be in the barrel region
          double pttag = jtpt[itag];
          double ptprobe = _pt[iprobe];
          if (jp::do3dHistos) {
            double asymm = (ptprobe - pttag)/(2*ptave);
            double mpf = met1*cos(DPhi(metphi1,_phi[itag]))/(2*ptave);
            for (auto alphaidx = 0u; alphaidx < h->alpharange.size(); ++alphaidx) {
              float alphasel = h->alpharange[alphaidx];
              if (alpha<alphasel) {
                // Val 10 = excluded, -10 = ok
                h->hdjasymm[alphaidx]  ->Fill(ptave,   etaprobe, asymm, _w);
                h->hdjmpf[alphaidx]    ->Fill(ptave,   etaprobe, mpf  , _w);
                h->hdjasymmtp[alphaidx]->Fill(pttag,   etaprobe, asymm, _w);
                h->hdjmpftp[alphaidx]  ->Fill(pttag,   etaprobe, mpf  , _w);
                //h->hdjasymmpt[alphaidx]->Fill(ptprobe, etaprobe, asymm, _w);
                //h->hdjmpfpt[alphaidx]  ->Fill(ptprobe, etaprobe, mpf  , _w);
              }
            }
          }

          // for composition vs eta
          if (alpha < 0.3 and pttag >= h->ptmin and pttag < h->ptmax) { // Alpha and trigger
            assert(h->pchftp_vseta); h->pchftp_vseta->Fill(etaprobe, jtchf[iprobe], _w);
            assert(h->pneftp_vseta); h->pneftp_vseta->Fill(etaprobe, (jtnef[iprobe]-jthef[iprobe]), _w);
            assert(h->pnhftp_vseta); h->pnhftp_vseta->Fill(etaprobe, (jtnhf[iprobe]-jthhf[iprobe]), _w);
            assert(h->pceftp_vseta); h->pceftp_vseta->Fill(etaprobe, jtcef[iprobe], _w);
            assert(h->pmuftp_vseta); h->pmuftp_vseta->Fill(etaprobe, jtmuf[iprobe], _w);
            assert(h->phhftp_vseta); h->phhftp_vseta->Fill(etaprobe, jthhf[iprobe], _w);
            assert(h->pheftp_vseta); h->pheftp_vseta->Fill(etaprobe, jthef[iprobe], _w);
            assert(h->ppuftp_vseta); h->ppuftp_vseta->Fill(etaprobe, jtbetaprime[iprobe]*jtchf[iprobe], _w);
          } // select pt bin for profiles vseta
        } // etatag < 1.3
      } // tag & probe
    } // dphi > 2.7
  } // ids

  if (h->ismcdir and _pass_qcdmet) {
    for (int jetidx = 0; jetidx != njt; ++jetidx) {
      double pt = jtpt[jetidx];
      if (pt>jp::recopt) { // pt visible
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

  PrintInfo(Form("\nOutput (HistosEta) stored in %s",_outfile->GetName()),true);

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
  if (jp::ismc) {
    triggers.push_back("mc");
    pt["mc"] = pair<double, double>(jp::recopt, jp::emax);
    pttrg["mc"] = jp::recopt;
  }
  if (jp::isdt or jp::domctrigsim) {
    // This is done both for data and MC, because why not?
    for (unsigned itrg = 0; itrg != jp::notrigs; ++itrg) {
      string trg = jp::triggers[itrg];
      triggers.push_back(trg);
      double pt1 = jp::trigranges[itrg][0];
      double pt2 = jp::trigranges[itrg][1];
      pt[trg] = pair<double, double>(pt1, pt2);
      double pt0 = jp::trigthr[itrg];
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

  if (jp::debug) {
    if (h == _mchistos.begin()->second[0]) {
      PrintInfo(Form("Triggers size: %lu",_trigs.size()));
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit) cout << *trgit << ", ";
      PrintInfo(Form("(%s)",h->trigname.c_str()));
    }
  }

  // check if required trigger fired
  if (!fired or !_pass) return;

  // Tag & probe hoods
  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_pass_qcdmet and i0>=0 and _jetids[i0] and jtpt[i0]>jp::recopt and i1>=0 and _jetids[i1] and jtpt[i1]>jp::recopt) { // Leading jets
    double dphi = DPhi(jtphi[i0], jtphi[i1]);
    if (dphi > 2.7) { // Back-to-back condition
      double pt3 = ((i2>=0 and jtpt[i2]>jp::recopt) ? jtpt[i2] : 0.);
      double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
      double alpha = pt3/ptave;

      for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
        int itag = jt3leads[itag_lead];
        int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
        double etatag = jteta[itag];
        double etaprobe = jteta[iprobe];
        if (etatag < 1.3) {
          if (jp::do3dHistos) {
            double pttag = _pt[itag];
            double ptprobe = _pt[iprobe];
            double asymm = (ptprobe - pttag)/(2*ptave);
            double asymmtp = (ptprobe - pttag)/(2*pttag);
            double ptresp_tag = _recopt[itag]/_genpt[itag];
            double ptresp_probe = _recopt[iprobe]/_genpt[iprobe];
            for (unsigned alphaidx = 0; alphaidx < h->alpharange.size(); ++alphaidx) {
              float alphasel = h->alpharange[alphaidx];
              if (alpha<alphasel) {
                h->hdjasymmtp     [alphaidx]->Fill(pttag, etaprobe, asymmtp,      _w);
                h->hdjresptp_tag  [alphaidx]->Fill(pttag, etaprobe, ptresp_tag,   _w);
                h->hdjresptp_probe[alphaidx]->Fill(pttag, etaprobe, ptresp_probe, _w);
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

  PrintInfo(Form("\nOutput (HistosMC) stored in %s",_outfile->GetName()),true);

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
  if (jp::isdt and h->lums[run][lbn]==0) {

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
        h->c_betaprime[t][run] = 0;
        h->t_trgtp[t][run] = 0;
        h->c_chftp[t][run] = 0;
        h->c_neftp[t][run] = 0;
        h->c_nhftp[t][run] = 0;
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
          h->c_betaprime[t][run] += jtbetaprime[jetidx];
        }

        int itag = (jetidx==0 ? 1 : 0);
        if (jetidx<2 and dphi > 2.7 and pt3 < jtpt[itag] and fabs(jteta[itag]) < 1.3 and
            jtpt[itag] > h->pt[t] and _jetids[itag]) {
          ++h->t_trgtp[t][run];
          h->c_chftp[t][run] += jtchf[jetidx];
          h->c_neftp[t][run] += (jtnef[jetidx]-jthef[jetidx]);
          h->c_nhftp[t][run] += (jtnhf[jetidx]-jthhf[jetidx]);
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


// Initialize basic histograms for trigger and eta bins
void HistosFill::InitAll(string name) {

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo(Form("InitAll(%s):",name.c_str()));
  PrintMemInfo();

  TDirectory *curdir = gDirectory;
  TFile *f = _outfile;
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  bool enteroutdir = f->cd(name.c_str());
  assert(enteroutdir);
  TDirectory *dir = f->GetDirectory(name.c_str()); assert(dir);
  dir->cd();

  HistosAll *h = new HistosAll(dir);
  _allhistos[name] = h;

  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("InitAll() finished:");
  PrintMemInfo();
} // InitAll


// Fill all-trigs histograms
void HistosFill::FillAll(string name)
{
  HistosAll *h = _allhistos[name];
  assert(h);

  if (_pass_qcdmet and njt>1 and _jetids[0]) {

    // Tag & probes method
    // We do implicit binning or histogramming by tag pt, so no need to look at firing triggers.
    double pttag = jtpt[0];
    double etatag = jteta[0];
    if (fabs(etatag)<1.3 and pttag>jp::wwptrange[0]) {
      // We cannot use the indices i0, i1 and i2 here due to general mapping issues #SAD
      // That is, even if the ordering of jets might be different in the original scheme (0,1,2,...),
      // the final scheme (i1,i2,i3) and the nol2l3 scheme, the most straightforward method (and a
      // good approximation) is to assume that the original ordering is conserved.
      int ptbin = -1;
      for (unsigned ptloc = 0; ptloc < jp::nwwpts; ++ptloc) {
        if (pttag < jp::wwptrange[ptloc+1]) {
          ptbin = ptloc;
          break;
        }
      }
      // We want to be inside a supported bin and we check the quality of the four leading jets
      if (ptbin>0 and (njt<2 or _jetids[1]) and (njt<3 or _jetids[2]) and (njt<4 or _jetids[3]) and (njt<5 or _jetids[4])) {
        // This code sector is separated from standard weighting strategies: for data, no weight is applied
        double wt = 1.0;
        if (!jp::isdt) wt = _w;

        // Let's create a universal tag unit vector
        double ju_px = cos(jtphi[0]);
        double ju_py = sin(jtphi[0]);

        // Let's fill MET histos and pttag histos
        double pttag_nol2l3 = jtpt_nol2l3[0];
        double projmet = ju_px * cos(metphi1) + ju_py * sin(metphi1);
        double projmet_nol2l3 = ju_px * cos(metphi1_nol2l3) + ju_py * sin(metphi1_nol2l3);
        h->pmetave         ->Fill(pttag,       projmet       *met1,       wt);
        h->pmetave_nol2l3  ->Fill(pttag_nol2l3,projmet_nol2l3*met1_nol2l3,wt);
        h->ppttagave       ->Fill(pttag,       pttag,                     wt);
        h->ppttagave_nol2l3->Fill(pttag_nol2l3,pttag_nol2l3,              wt);

        // First jet
        int pidx = 1;
        if (pidx<njt and jtpt[pidx]>jp::recopt) {
          double proj1 = ju_px*cos(jtphi[pidx]) + ju_py*sin(jtphi[pidx]);
          double etaj1       = jteta      [pidx];
          double ptj1        = jtpt       [pidx];
          double ptj1_nol2l3 = jtpt_nol2l3[pidx];
          h->h2njet1        [ptbin]->Fill(ptj1,       etaj1,                  wt);
          h->h2njet1_nol2l3 [ptbin]->Fill(ptj1_nol2l3,etaj1,                  wt);
          h->p2ptjet1       [ptbin]->Fill(ptj1,       etaj1,proj1*ptj1,       wt);
          h->p2ptjet1_nol2l3[ptbin]->Fill(ptj1_nol2l3,etaj1,proj1*ptj1_nol2l3,wt);
        }
        // Second jet
        pidx = 2;
        if (pidx<njt and jtpt[pidx]>jp::recopt) {
          double proj2 = ju_px*cos(jtphi[pidx]) + ju_py*sin(jtphi[pidx]);
          double etaj2       = jteta      [pidx];
          double ptj2        = jtpt       [pidx];
          double ptj2_nol2l3 = jtpt_nol2l3[pidx];
          h->h2njet2        [ptbin]->Fill(ptj2,       etaj2,                  wt);
          h->h2njet2_nol2l3 [ptbin]->Fill(ptj2_nol2l3,etaj2,                  wt);
          h->p2ptjet2       [ptbin]->Fill(ptj2,       etaj2,proj2*ptj2,       wt);
          h->p2ptjet2_nol2l3[ptbin]->Fill(ptj2_nol2l3,etaj2,proj2*ptj2_nol2l3,wt);
        }
        // Third jet
        pidx = 3;
        if (pidx<njt and jtpt[pidx]>jp::recopt) {
          double proj3 = ju_px*cos(jtphi[pidx]) + ju_py*sin(jtphi[pidx]);
          double etaj3       = jteta      [pidx];
          double ptj3        = jtpt       [pidx];
          double ptj3_nol2l3 = jtpt_nol2l3[pidx];
          h->h2njet3        [ptbin]->Fill(ptj3,       etaj3,                  wt);
          h->h2njet3_nol2l3 [ptbin]->Fill(ptj3_nol2l3,etaj3,                  wt);
          h->p2ptjet3       [ptbin]->Fill(ptj3,       etaj3,proj3*ptj3,       wt);
          h->p2ptjet3_nol2l3[ptbin]->Fill(ptj3_nol2l3,etaj3,proj3*ptj3_nol2l3,wt);
        }
        // Fourth jet & friends
        pidx = 4;
        if (pidx<njt and jtpt[pidx]>jp::recopt) {
          // First, we seek for the jet<->bin connections
          map<int,int> jet2bin;
          map<int,int> jet2bin_nol2l3;
          for (int jidx = 4; jidx < njt; ++jidx) {
            jet2bin       [jidx] = h->h2njet4plus[ptbin]->FindBin(jtpt       [jidx],jteta[jidx]);
            jet2bin_nol2l3[jidx] = h->h2njet4plus[ptbin]->FindBin(jtpt_nol2l3[jidx],jteta[jidx]);
          }
          // Loop over the remaining jets (we do not consider jetids here).
          // For each jet, we fill the bin that corresponds to it with a sum over the jets within that bin.
          for (pidx = 4; pidx < njt; ++pidx) {
            int assoc        = jet2bin       [pidx];
            int assoc_nol2l3 = jet2bin_nol2l3[pidx];
            // This is an indicator that this jet has already been considered
            if (assoc==-1 and assoc_nol2l3==-1) continue;

            double cumul        = 0;
            double cumul_nol2l3 = 0;
            // We loop over jets, including the current jet, and take the cumulative sum corresponding to the bin of this jet
            // A jet already considered is marked with the tag '-1'
            for (int jidx = pidx; jidx < njt; ++jidx) {
              bool match        = (assoc == jet2bin[jidx]);
              bool match_nol2l3 = (assoc_nol2l3 == jet2bin_nol2l3[jidx]);
              if (match or match_nol2l3) {
                double proj4p = ju_px*cos(jtphi[jidx]) + ju_py*sin(jtphi[jidx]);
                if (match) {
                  cumul        += proj4p*jtpt       [jidx];
                  jet2bin       [jidx] = -1;
                }
                if (match_nol2l3) {
                  cumul_nol2l3 += proj4p*jtpt_nol2l3[jidx];
                  jet2bin_nol2l3[jidx] = -1;
                }
              }
            }
            double etaj4 = jteta[pidx];
            if (assoc!=-1) {
              double ptj4 = jtpt[pidx];
              h->h2njet4plus [ptbin]->Fill(ptj4,etaj4,      wt);
              h->p2ptjet4plus[ptbin]->Fill(ptj4,etaj4,cumul,wt);
            }
            if (assoc_nol2l3!=-1) {
              double ptj4_nol2l3 = jtpt_nol2l3[pidx];
              h->h2njet4plus_nol2l3 [ptbin]->Fill(ptj4_nol2l3,etaj4,             wt);
              h->p2ptjet4plus_nol2l3[ptbin]->Fill(ptj4_nol2l3,etaj4,cumul_nol2l3,wt);
            }
          } // 4th jet & friends jet loop
        } // 4th jet & friends
      }
    }
  }

} // FillAll


// Write and delete histograms
void HistosFill::WriteAll()
{
  for (auto &hh : _allhistos) {
    TDirectory *curdir = gDirectory;
    string name = hh.first;
    TFile *f = _outfile;
    assert(f && !f->IsZombie());
    TDirectory *dir = f->GetDirectory(name.c_str()); assert(dir);
    dir->cd();
    HistosAll *h = hh.second;
    // We create correctly scaled end user histograms to reduce future error rate
    for (unsigned idx = 0; idx < jp::nwwpts; ++idx) {
      // Total event counts are best fetched from the MET histos
      double sf        = h->pmetave       ->GetBinEntries(idx+1);
      double sf_nol2l3 = h->pmetave_nol2l3->GetBinEntries(idx+1);
      sf        = (sf>0        ? 1.0/sf        : 1.0);
      sf_nol2l3 = (sf_nol2l3>0 ? 1.0/sf_nol2l3 : 1.0);

      int num = jp::wwptrange[idx];
      string number = std::to_string(num);
      // Get the histos as projections
      h->h2ptjet1           [idx] = h->p2ptjet1           [idx]->ProjectionXY((string("h2ptj1_")        +number+string("GeV")).c_str(),"e");
      h->h2ptjet2           [idx] = h->p2ptjet2           [idx]->ProjectionXY((string("h2ptj2_")        +number+string("GeV")).c_str(),"e");
      h->h2ptjet3           [idx] = h->p2ptjet3           [idx]->ProjectionXY((string("h2ptj3_")        +number+string("GeV")).c_str(),"e");
      h->h2ptjet4plus       [idx] = h->p2ptjet4plus       [idx]->ProjectionXY((string("h2ptj4p_")       +number+string("GeV")).c_str(),"e");
      h->h2ptjet1_nol2l3    [idx] = h->p2ptjet1_nol2l3    [idx]->ProjectionXY((string("h2ptj1_nol2l3_") +number+string("GeV")).c_str(),"e");
      h->h2ptjet2_nol2l3    [idx] = h->p2ptjet2_nol2l3    [idx]->ProjectionXY((string("h2ptj2_nol2l3_") +number+string("GeV")).c_str(),"e");
      h->h2ptjet3_nol2l3    [idx] = h->p2ptjet3_nol2l3    [idx]->ProjectionXY((string("h2ptj3_nol2l3_") +number+string("GeV")).c_str(),"e");
      h->h2ptjet4plus_nol2l3[idx] = h->p2ptjet4plus_nol2l3[idx]->ProjectionXY((string("h2ptj4p_nol2l3_")+number+string("GeV")).c_str(),"e");
      // Multiply by bin event counts
      h->h2ptjet1           [idx]->Multiply(h->h2njet1           [idx]);
      h->h2ptjet2           [idx]->Multiply(h->h2njet2           [idx]);
      h->h2ptjet3           [idx]->Multiply(h->h2njet3           [idx]);
      h->h2ptjet4plus       [idx]->Multiply(h->h2njet4plus       [idx]);
      h->h2ptjet1_nol2l3    [idx]->Multiply(h->h2njet1_nol2l3    [idx]);
      h->h2ptjet2_nol2l3    [idx]->Multiply(h->h2njet2_nol2l3    [idx]);
      h->h2ptjet3_nol2l3    [idx]->Multiply(h->h2njet3_nol2l3    [idx]);
      h->h2ptjet4plus_nol2l3[idx]->Multiply(h->h2njet4plus_nol2l3[idx]);
      // Divide by total event counts
      h->h2ptjet1           [idx]->Scale(sf);
      h->h2ptjet2           [idx]->Scale(sf);
      h->h2ptjet3           [idx]->Scale(sf);
      h->h2ptjet4plus       [idx]->Scale(sf);
      h->h2ptjet1_nol2l3    [idx]->Scale(sf_nol2l3);
      h->h2ptjet2_nol2l3    [idx]->Scale(sf_nol2l3);
      h->h2ptjet3_nol2l3    [idx]->Scale(sf_nol2l3);
      h->h2ptjet4plus_nol2l3[idx]->Scale(sf_nol2l3);
    }
    curdir->cd();
  }

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteAll():");
  PrintMemInfo();

  for (auto histit = _allhistos.begin(); histit != _allhistos.end(); ++histit) {
    HistosAll *h = histit->second;
    delete h;
  } // for histit

  PrintInfo(Form("\nOutput (HistosAll) stored in %s",_outfile->GetName()),true);

  // Report memory usage to avoid malloc problems when writing file
  PrintInfo("WriteAll() finished:");
  PrintMemInfo();
} // WriteAll


void HistosFill::FillJetID(vector<bool> &id)
{
  assert(int(id.size())==njt);

  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    id[jetidx] = jtidtight[jetidx];
    float etabs = fabs(jteta[jetidx]);
    if (jp::yid==0) { // 2016 specialities
      // Loose ID is still in use
      if (etabs>2.5) id[jetidx] = jtidloose[jetidx];
      // TightLepVeto
      if (etabs<=2.7 and (jtcef[jetidx]>=0.90 or jtmuf[jetidx]>=0.80)) id[jetidx] = false;
    } else if (jp::yid==1) { // 2017 specialities
      if (etabs<=2.7 and (jtcef[jetidx]>=0.80 or jtmuf[jetidx]>=0.80)) id[jetidx] = false;
    } else if (jp::yid==2) {
      if (etabs<=2.7 and (jtcef[jetidx]>=0.80 or jtmuf[jetidx]>=0.80)) id[jetidx] = false;
    }

    if (jp::doVetoHot) {
      // Abort if one of the leading jets is in a difficult zone
      assert(h2HotExcl);
      bool good = h2HotExcl->GetBinContent(h2HotExcl->FindBin(jteta[jetidx],jtphi[jetidx])) <= 0;
      id[jetidx] = (id[jetidx] and good);
    }
  }
} // FillJetID


// Load good run and LS information
bool HistosFill::LoadJSON()
{
  PrintInfo(string("Processing LoadJSON() with ") + jp::json + " ...",true);
  ifstream file(jp::json, ios::in);
  if (!file.is_open()) return false;
  char c;
  string s, s2;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (jp::debug) PrintInfo(Form("\"%d\": ",rn),true);

    while (file.get(c) and c==' ') {};
    if (jp::debug) { PrintInfo(Form("%c",c),true); assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (jp::debug) PrintInfo(Form("[%d,%d,%s]",ls1,ls2,s1),true);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (jp::debug and !endrun and s2!=",") { PrintInfo(string("s1: ")+s2,true); assert(s2==","); }
    } // while ls
    if (jp::debug) PrintInfo("",true);

    if (s2=="]}") continue;
    else if (jp::debug and s2!="],") PrintInfo(string("s2: ")+s2,true); assert(s2=="],");
  } // while run
  if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }

  PrintInfo(string("Called LoadJSON() with ") + jp::json + ":",true);
  PrintInfo(Form("Loaded %d good runs and %d good lumi sections",nrun,nls));
  return true;
} // LoadJSON


// Load luminosity information
bool HistosFill::LoadLumi()
{
  PrintInfo(string("Processing LoadLumi() with ") + jp::lumifile + "...",true);

  // Check lumi against the list of good runs
  const int a_goodruns[] = {};
  const int ngoodruns = sizeof(a_goodruns)/sizeof(a_goodruns[0]);
  set<int> goodruns;
  if (ngoodruns>0) { // This is an old remnant
    for (int runidx = 0; runidx != ngoodruns; ++runidx)
      goodruns.insert(a_goodruns[runidx]);

    for (auto runit = goodruns.begin(); runit != goodruns.end(); ++runit) cout << *runit << ", ";
    cout << endl;
  }
  set<pair<int, int> > nolums;

  ifstream f(jp::lumifile, ios::in);
  if (!f.is_open()) return false;
  float secLS = 2.3310e+01;
  string s;
  int rn, fill, ls, ifoo;
  float del, rec, avgpu, energy;
  char sfoo[512];
  bool getsuccess1 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess1) return false;
  PrintInfo(string("\nstring: ") + s + " !",true);

  // HOX: the lumi file format has been changing. Change the conditions when needed.
  if (s!="#Data tag : 19v2 , Norm tag: None") return false;

  bool getsuccess2 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess2) return false;
  PrintInfo(string("\nstring: ") + s + " !",true);
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

    if (jp::debug) PrintInfo(Form("Run %d ls %d lumi %f/pb",rn,ls,rec*1e-6),true);

    if (skip) { // The user should know if this happens, since we can choose to use only STABLE BEAMS
      if (skip) PrintInfo(string("Skipping line (effects the recorded lumi):\n")+s,true);
      skip = false;
      continue;
    }

    if (_lums[rn][ls]!=0) return false;
    if (_avgpu[rn][ls]!=0) return false;
    // lumiCalc.py returns lumi in units of mub-1 (=>nb-1=>pb-1)
    double lum = rec*1e-6;
    double lum2 = del*1e-6;
    if (lum==0 and goodruns.find(rn)!=goodruns.end() and (!jp::dojson or _json[rn][ls]==1))
      nolums.insert(pair<int, int>(rn,ls));

    _avgpu[rn][ls] = avgpu; // * 69000. / 78400.; // brilcalc --minBiasXsec patch
    _lums[rn][ls] = lum;
    _lums2[rn][ls] = lum2;
    lumsum += lum;
    if (goodruns.find(rn)!=goodruns.end()) // Apr 17
      lumsum_good += lum;
    if ((!jp::dojson || _json[rn][ls]))
      lumsum_json += lum;
    ++nls;
    if (nls>100000000) return false;
  }

  PrintInfo(Form("Called LoadLumi() with %s:\nLoaded %lu runs with %d lumi sections containing %f"
                 " pb-1 of data,\n of which %f pb-1 is in good runs (%f%%)\nThis corresponds to %f"
                 " hours of data-taking\nThe JSON file contains %f pb-1 (%f%%)",
                 jp::lumifile,_lums.size(),nls,lumsum,lumsum_good,
                 100.*lumsum_good/lumsum,nls*secLS/3600,lumsum_json,100.*lumsum_json/lumsum),true);

  // Report any empty lumi section
  if (nolums.size()!=0) {
    PrintInfo(Form("Warning, found %lu non-normalizable LS:",nolums.size()),true);
    for (auto lumit = nolums.begin(); lumit != nolums.end(); ++lumit) {
      cout << " ["<<lumit->first<<","<<lumit->second;
      auto lumit2 = lumit;
      ++lumit2;
      if (lumit2->first!=lumit->first or lumit2->second!=lumit->second+1) cout << "]";
      else {
        for (int lumadd = 0; lumit2!=nolums.end() and lumit2->first==lumit->first and
                             lumit2->second==lumit->second+lumadd+1; ++lumadd, ++lumit2) {};
        lumit = --lumit2;
        cout << "-" << lumit->second << "]";
      }
    } // for lumit
    cout << endl;
  } // nolums
  return true;
} // LoadLumi


bool HistosFill::LoadPuProfiles()
{
  if (jp::PUIOVidx>=jp::PUIOVs.size()) {
    PrintInfo(Form("Selected PU IOV too great (%u vs. %zu), check the value of jp::PUIOVidx!",jp::PUIOVidx,jp::PUIOVs.size()));
    return false;
  }
  string datafile = jp::pudtpath + jp::PUIOVs[jp::PUIOVidx] + "pileup_DT.root";
  string mcfile   = jp::pumcpath;
  if (jp::isnu)      mcfile += "pileup_NU.root";
  else if (jp::ishw) mcfile += jp::puhwfile;
  else if (jp::ispy) {
    if (jp::pthatbins)   mcfile += "pileup_Pthat.root";
    else if (jp::htbins) mcfile += "pileup_P8MG.root";
    else                 mcfile += jp::pup8file;
  } else {
    PrintInfo("Problems with PU file types!",true);
    return false;
  }

  PrintInfo(Form("Processing LoadPuProfiles() using %s and %s ...",datafile.c_str(),mcfile.c_str()),true);

  TDirectory *curdir = gDirectory;
  // Load pile-up files and hists from them
  TFile *f_pudist = new TFile(datafile.c_str(), "READ");
  if (!f_pudist or f_pudist->IsZombie()) { PrintInfo(string("Unable to read PU data file")+datafile,true); return false; }
  TFile *fpumc = new TFile(mcfile.c_str(),"READ");
  if (!fpumc or fpumc->IsZombie()) { PrintInfo(string("Unable to read PU MC file")+mcfile,true); return false; }

  _pumc = dynamic_cast<TH1D*>(fpumc->Get("pileupmc")->Clone("pumchelp"));
  if (!_pumc) return false;
  if (jp::isnu) { // In the neutrino gun samples we look at the hardest PU event, so need to shift PU by -1
    _pumc->SetBinContent(0,_pumc->GetBinContent(0)+_pumc->GetBinContent(1));
    for (int idx = 1; idx < _pumc->GetNbinsX(); ++idx)
      _pumc->SetBinContent(idx,_pumc->GetBinContent(idx+1));
  }
  double maxmcpu = _pumc->GetMaximum();
  _pumc->Scale(1.0/maxmcpu);
  int lomclim = _pumc->FindFirstBinAbove(0.01);
  int upmclim = _pumc->FindLastBinAbove(0.01);
  int maxmcbin = _pumc->FindFirstBinAbove(0.999);
  PrintInfo(Form("Maximum bin: %d for MC",maxmcbin),true);
  PrintInfo(Form("Hazardous pu below & above: %f, %f",_pumc->GetBinLowEdge(lomclim),_pumc->GetBinLowEdge(upmclim+1)),true);
  int nbinsmc = _pumc->GetNbinsX();
  int kmc = _pumc->FindBin(33); // Check that pu=33 occurs at the same place as for data

  // For data, load each trigger separately
  for (auto &t : jp::triggers) {
    _pudist[t] = dynamic_cast<TH1D*>(f_pudist->Get(t));
    if (!_pudist[t]) {
      PrintInfo(Form("The trigger %s was not found in the DT pileup file!",t),true);
      return false;
    }
    int nbinsdt = _pudist[t]->GetNbinsX();
    int kdt = _pudist[t]->FindBin(33);
    if (kdt!=kmc or nbinsdt!=nbinsmc) {
      PrintInfo("The pileup histogram dt vs mc binning or range do not match (dt left mc right):",true);
      PrintInfo(Form(" Bins: dt:%d mc:%d",nbinsdt,nbinsmc),true);
      PrintInfo(Form(" Pu=33 bin: dt:%d mc:%d",kdt,kmc),true);
      return false;
    }
    double maxdtpu = _pudist[t]->GetMaximum();
    int lodtlim = _pudist[t]->FindFirstBinAbove(maxdtpu/100.0);
    int updtlim = _pudist[t]->FindLastBinAbove(maxdtpu/100.0);
    int maxdtbin = _pudist[t]->FindFirstBinAbove(0.999*maxdtpu);

    for (int bin = 0; bin < lomclim; ++bin) // Set fore-tail to zero
      _pudist[t]->SetBinContent(bin,0.0);
    PrintInfo(Form("Maximum bin: %d for DT trg %s",maxdtbin,t),true);
    PrintInfo(Form("Hazardous pu below & above: %f, %f",_pudist[t]->GetBinLowEdge(lodtlim),_pudist[t]->GetBinLowEdge(updtlim+1)),true);
    _pudist[t]->Divide(_pumc);
  }
  PrintInfo("Finished processing pileup histos!",true);

  curdir->cd();
  return true;
} // LoadPuProfiles


Int_t HistosFill::FindMCSliceIdx(string filename)
{
  int sliceIdx = 0;
  bool sliceFound = false;
  if (jp::pthatbins) {
    for (auto &fname : jp::pthatfiles) {
      regex rfile(fname);
      std::cmatch mfile;
      if (std::regex_search(filename.c_str(),mfile,rfile)) {
        sliceFound = true;
        break;
      }
      ++sliceIdx;
    }
  } else if (jp::htbins) {
    for (auto &fname : jp::htfiles) {
      regex rfile(fname);
      std::cmatch mfile;
      if (std::regex_search(filename.c_str(),mfile,rfile)) {
        sliceFound = true;
        break;
      }
      ++sliceIdx;
    }
  }
  if (sliceFound) return sliceIdx;
  return -1;
}


// Check that the correct tree is open in the chain
Long64_t HistosFill::LoadTree(Long64_t entry)
{
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;

  // A new tree is opened
  if (fChain->GetTreeNumber() != fCurrent) {

    fCurrent = fChain->GetTreeNumber();
    PrintInfo(Form("Opening tree number %d", fChain->GetTreeNumber()));

    if (jp::fetchMETFilters) {
      // Reload the MET filters and print them
      if (!GetFilters()) {
        PrintInfo("Failed to load DT filters. Check that the SMPJ tuple has the required histograms. Aborting...");
        return -4;
      }
      PrintInfo(Form("Tree %d MET filters:",fCurrent),true);
      for (auto &flt : _availFlts) PrintInfo(flt,true);
    }

    if (jp::isdt) {
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
        if (trigi%(jp::notrigs+1)==jp::notrigs) *ferr << endl;
      }
      *ferr << endl << flush;
    } else if (jp::pthatbins or jp::htbins) {
      // If there are two pthat files with the same pthat range, we convey this information through "prevweight"
      bool htmode = !jp::pthatbins;
      const char* bintag = htmode ? "HT" : "Pthat";
      Long64_t noevts = fChain->GetTree()->GetEntries();
      if (_binnedmcrepeats<=0) {
        const char* fname = fChain->GetCurrentFile()->GetName();
        // Check the position of the current file in the list of file names
        int sliceIdx = FindMCSliceIdx(fname);
        if (sliceIdx<0) {
          PrintInfo(Form("%s slice file name contradictory %s. Aborting...",bintag,fname));
          return -3;
        }
        if ((!htmode and jp::pthatsigmas[sliceIdx]<=0) or (htmode and jp::htsigmas[sliceIdx]<=0)) {
          PrintInfo(Form("Suspicious %s information for file %s. Aborting...",bintag,fname));
          return -3;
        }
        PrintInfo(Form("%s bin changing.\nFile %d %s, position %lld with %lld events.",bintag,sliceIdx,fname,centry,noevts),true);

        _binnedmcrepeats = 0;
        // We have a look if the next tree in the chain has the same range
        Long64_t next = entry-centry;
        while (fChain->GetTreeNumber()<fChain->GetNtrees()-1) {
          next = next+fChain->GetTree()->GetEntries();
          if (next>=_ntot) break;
          Long64_t nentry = fChain->LoadTree(next);
          if (nentry < 0) break;
          Long64_t nextevts = fChain->GetTree()->GetEntries();
          const char* nextname = fChain->GetCurrentFile()->GetName();
          int nextIdx = FindMCSliceIdx(nextname);
          if (nextIdx<0) {
            PrintInfo(Form("%s slice file name contradictory %s. Aborting...",bintag,nextname));
            return -3;
          }
          if (sliceIdx!=nextIdx) break;
          PrintInfo(Form("File extension %d %s, position %lld with %lld events",nextIdx,nextname,nentry,nextevts),true);
          noevts += nextevts;
          ++_binnedmcrepeats;
        }
        // Normalization with the amount of entries within the current tree
        double sigmacurr = htmode ? jp::htsigmas[sliceIdx  ] : jp::pthatsigmas[sliceIdx  ];
        double slicecurr = htmode ? jp::htranges[sliceIdx  ] : jp::pthatranges[sliceIdx  ];
        double slicenext = htmode ? jp::htranges[sliceIdx+1] : jp::pthatranges[sliceIdx+1];
        _binnedmcweight  = sigmacurr/noevts;
        // This is a normalization procedure by the luminosity of the furthest PtHat/HT bin.
        // In practice, it does not hurt if the normalevts number is arbitrary.
        _binnedmcweight /= (htmode ? (jp::htsigmas.back()/jp::htnormalevts) : (jp::pthatsigmas.back()/jp::pthatnormalevts));
        PrintInfo(Form("The given slice has the %s range [%f,%f]\nWeight: %f, with a total x-sec %f pb and %lld events.",
                  bintag,slicecurr,slicenext,_binnedmcweight,sigmacurr,noevts),true);
      } else {
        --_binnedmcrepeats;
        PrintInfo(Form("%s bin remains the same while file is changing.\nFile %s\nWeight: %f",
                  bintag,fChain->GetCurrentFile()->GetName(),_binnedmcweight),true);
      }
    }
    // slices with PtHat/HT bins
  }
  return centry;
}


// Update the available MET filter types for each new tree
bool HistosFill::GetFilters()
{
  TH1F *filters = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/FilterActive")); assert(filters);
  TAxis *xax = filters->GetXaxis();
  regex filter("Flag_([a-zA-Z0-9]+)");

  _availFlts.clear();
  bool empty  = false;
  bool passall = false;
  for (int fltidx = xax->GetFirst(); fltidx <= xax->GetLast(); ++fltidx) {
    string fltName = xax->GetBinLabel(fltidx);

    if (fltName=="") {
      empty = true;
    } else if (fltName=="PassAll") {
      if (empty) {
        PrintInfo("PassAll MET filter in the wrong position!",true);
        return false;
      }
      passall = true;
    } else if (std::regex_match(fltName,filter)) {
      if (empty or passall) {
        PrintInfo(Form("MET filter %s in the wrong position!",fltName.c_str()),true);
        return false;
      }
      string stripName = std::regex_replace(fltName, filter, "$1", std::regex_constants::format_no_copy);
      _availFlts.push_back(stripName);
    } else {
      PrintInfo(Form("Unknown filter type %s",fltName.c_str()),true);
      return false;
    }
  }

  return _availFlts.size()>0;
}

// Update the available trigger types for each new tree
bool HistosFill::GetTriggers()
{
  TH1F *triggers = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerNames")); assert(triggers);
  TAxis *xax = triggers->GetXaxis();
  TH1F *usedtrigs = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerPass")); assert(usedtrigs);
  TAxis *uxax = usedtrigs->GetXaxis();

  regex zbs("HLT_ZeroBias_v([0-9]+)");
  regex pfjet("HLT_PFJet([0-9]+)_v([0-9]+)");
  regex ak8("HLT_AK8PFJet([0-9]+)_v[0-9]+");
  regex jetht("HLT_PFHT([0-9]+)_v[0-9]+");
  regex hipfjet("HLT_HI[AK4]*PFJet([0-9]+)_v[0-9]+");
  regex hipfjetfwd("HLT_HI[AK4]*PFJet[FWwDd]*([0-9]+)[FWwDd]*_v[0-9]+");
  regex hiak8("HLT_HIAK8PFJet([0-9]+)_v[0-9]+");
  regex hiak8fwd("HLT_HIAK8PFJet[FWwDd]*([0-9]+)[FWwDd]*_v[0-9]+");

  // List triggers with actual contents
  map<string,unsigned int> utrigs;
  bool zbcase = false;
  for (int trgidx = uxax->GetFirst(); trgidx <= uxax->GetLast(); ++trgidx) {
    string trgName = uxax->GetBinLabel(trgidx);
    if (trgName=="") continue;
    utrigs[trgName] = usedtrigs->GetBinContent(trgidx);
    if (std::regex_match(trgName,zbs)) zbcase = true;
  }

  auto &eraLumis  = jp::triglumiera[_eraIdx];
  assert(eraLumis.size()==jp::notrigs);
  _availTrigs.clear();
  _goodTrigs.clear();
  _goodWgts.clear();
  for (int trgidx = xax->GetFirst(); trgidx <= xax->GetLast(); ++trgidx) {
    string trgName = xax->GetBinLabel(trgidx);
    if (trgName.compare("")==0) continue; // Ignore empty places on x-axis

    string trigger = "x"; // HLT_PFJet are given non-empty trigger names
    if (std::regex_match(trgName,zbs)) {
      trigger = "jt0"; // Zero for ZeroBias
      if (zbcase and utrigs.find(trgName) != utrigs.end()) {
        unsigned thrplace = 0;
        if (thrplace < jp::notrigs) {
          _goodTrigs.push_back(_availTrigs.size());
          if (jp::useeraweights) {
            // Get a weight for the current trig version normalized with the average of all triggers
            double eraLumiWgt  =     eraLumis.back()/    eraLumis[thrplace];
            double yearLumiWgt = jp::triglumi.back()/jp::triglumi[thrplace];
            // Ideally this is a constant common for all triggers.
            // Sadly, this is not true and thus this takes into account the time dependence of
            _goodWgts.push_back(eraLumiWgt/yearLumiWgt);
            PrintInfo(Form("Trigger %s responding loud and clear with %u events and relative weight %f!",trgName.c_str(),utrigs[trgName],_goodWgts.back()),true);
          } else {
            _goodWgts.push_back(1.0);
            PrintInfo(Form("Trigger %s responding loud and clear with %u events!",trgName.c_str(),utrigs[trgName]),true);
          }
        } else { // No trig era weighting: no relative weights
          PrintInfo(Form("The trigger %s is available, but not supported",trgName.c_str()),true);
        }
      }
    } else if (std::regex_match(trgName,hiak8)) {
      trigger=std::regex_replace(trgName, hiak8, "ak8jt$1", std::regex_constants::format_no_copy);
    } else if (std::regex_match(trgName,hiak8fwd)) {
      trigger=std::regex_replace(trgName, hiak8fwd, "ak8jt$1fwd", std::regex_constants::format_no_copy);
    } else if (std::regex_match(trgName,hipfjet)) {
      trigger=std::regex_replace(trgName, hipfjet, "jt$1", std::regex_constants::format_no_copy);
      if (_eraIdx==1 and jp::yid==2 and !zbcase and utrigs.find(trgName) != utrigs.end()) {
        double trigthr = std::stod(std::regex_replace(trgName, hipfjet, "$1", std::regex_constants::format_no_copy));
        unsigned thrplace = static_cast<unsigned>(std::find(jp::trigthr.begin()+1,jp::trigthr.end(),trigthr)-jp::trigthr.begin());
        if (thrplace < jp::notrigs) {
          _goodTrigs.push_back(_availTrigs.size());
          PrintInfo(Form("Trigger %s responding loud and clear with %u events!",trgName.c_str(),utrigs[trgName]),true);
          _goodWgts.push_back(1.0);
        } else {
          PrintInfo(Form("The trigger %s is available, but not supported",trgName.c_str()),true);
        }
      }
    } else if (std::regex_match(trgName,hipfjetfwd)) {
      trigger=std::regex_replace(trgName, hipfjetfwd, "jt$1fwd", std::regex_constants::format_no_copy);
      if (_eraIdx==0 and jp::yid==2 and !zbcase and utrigs.find(trgName) != utrigs.end()) {
        double trigthr = std::stod(std::regex_replace(trgName, hipfjetfwd, "$1", std::regex_constants::format_no_copy));
        unsigned thrplace = static_cast<unsigned>(std::find(jp::trigthr.begin()+1,jp::trigthr.end(),trigthr)-jp::trigthr.begin());
        if (thrplace < jp::notrigs) {
          _goodTrigs.push_back(_availTrigs.size());
          PrintInfo(Form("Trigger %s responding loud and clear with %u events!",trgName.c_str(),utrigs[trgName]),true);
          _goodWgts.push_back(1.0);
        } else {
          PrintInfo(Form("The trigger %s is available, but not supported",trgName.c_str()),true);
        }
      }
    } else if (std::regex_match(trgName,pfjet)) {
      trigger=std::regex_replace(trgName, pfjet, "jt$1", std::regex_constants::format_no_copy);
      if (!zbcase and utrigs.find(trgName) != utrigs.end()) {
        double trigthr = std::stod(std::regex_replace(trgName, pfjet, "$1", std::regex_constants::format_no_copy));
        unsigned thrplace = static_cast<unsigned>(std::find(jp::trigthr.begin()+1,jp::trigthr.end(),trigthr)-jp::trigthr.begin());
        if (thrplace < jp::notrigs) {
          _goodTrigs.push_back(_availTrigs.size());
          if (jp::useeraweights) {
            // Get a weight for the current trig version normalized with the average of all triggers
            double eraLumiWgt  =     eraLumis.back()/    eraLumis[thrplace];
            double yearLumiWgt = jp::triglumi.back()/jp::triglumi[thrplace];
            // Ideally this is a constant common for all triggers.
            // Sadly, this is not true and thus this takes into account the time dependence of
            _goodWgts.push_back(eraLumiWgt/yearLumiWgt);
            PrintInfo(Form("Trigger %s responding loud and clear with %u events and relative weight %f!",trgName.c_str(),utrigs[trgName],_goodWgts.back()),true);
          } else {
            PrintInfo(Form("Trigger %s responding loud and clear with %u events!",trgName.c_str(),utrigs[trgName]),true);
            _goodWgts.push_back(1.0);
          }
        } else { // No trig era weighting: no relative weights
          PrintInfo(Form("The trigger %s is available, but not supported",trgName.c_str()),true);
        }
      }
    } else if (std::regex_match(trgName,ak8)) {
      trigger=std::regex_replace(trgName, ak8, "ak8jt$1", std::regex_constants::format_no_copy);
    } else if (std::regex_match(trgName,jetht)) {
      trigger=std::regex_replace(trgName, jetht, "jetht$1", std::regex_constants::format_no_copy);
    } else {
      PrintInfo(Form("Unknown trigger type %s",trgName.c_str()),true);
    }

    _availTrigs.push_back(trigger);
  }

  return _availTrigs.size()>0;
} // GetTriggers
