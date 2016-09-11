// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: June 2, 2015
// Updated: Aug 31, 2016

#define fillHistos_cxx
#include "fillHistos.h"

void fillHistos::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t ntot = fChain->GetEntries();
  Long64_t nskip = _jp_nskip;//0;
  nentries = (_jp_nentries==-1 ? ntot-nskip : min(ntot-nskip, _jp_nentries));
  assert(nskip+nentries);

  //nentries = 10; // debug
  //nentries = 1000;//very short test runs
  //nentries = 100000;//short test runs
  //nentries = 1000000;//medium test runs
  //nentries = 5000000; // lunch-break run for MC (with trigsim off)

  map<string, int> cnt; // efficiency counters

  ferr = new ofstream("fillHistos.log",ios::out);

  // Report memory usage to avoid malloc problems when writing file
  *ferr << endl << "Starting Loop() initialization:" << endl << flush;
  cout << endl << "Starting Loop() initialization:" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
  cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  _entries = ntot;//(_mc ? fChain->GetEntries() : ntot);
  _nbadevts_dup = _nbadevts_run = _nbadevts_ls = _nbadevts_lum = 0;
  _nbadevts_veto = _nbadevts_stream = 0;
  _bscounter_bad = _bscounter_good = _halocounter_bad = _halocounter_good = 0;
  _ecalcounter_good = _ecalcounter_bad = 0;
  _rhocounter_good = _rhocounter_bad = 0;
  _trgcounter = _evtcounter = _totcounter = 0;
  ecalveto = 0;

  // Set cross section weights for pThat bins
  hmcweight = 0;
  if (_jp_pthatbins) {
    double bins[]={5,15,30,50,80,120,170,300,470,600,800,1000,1400,1800,3500};
    const int nbins = sizeof(bins)/sizeof(bins[0])-1;
    double lums[] = {4.49e-5,1.32e-2,1.22e-1,1.04,8.40,5.32e+1,2.57e+2,5.51e+3,
                     5.68e+4,2.73e+5,2.20e+6,6.30e+6,2.02e+8,8.20e+8};
    const int nlums = sizeof(lums)/sizeof(lums[0]);
    assert(nbins==nlums);
    hmcweight = new TH1D("hmcweight",";#hat{p}_{T} (GeV)",nbins,bins);
    for (int i = 0; i != nbins; ++i) {
      hmcweight->SetBinContent(i+1, 1./lums[i]);
    } // for i
  }

  // Map of mu per (run,lumi)
  TFile *fmu = new TFile("pileup/MUperLSvsRUN_MB.root","READ");
  assert(fmu && !fmu->IsZombie());
  TH2F *h2mu = (TH2F*)fmu->Get("hLSvsRUNxMU"); assert(h2mu);
  //TH2F *h2mu = (TH2F*)fmu->Get("hLSvsRuNxMU_cleaned"); assert(h2mu);

  // map AK4 events to AK8 events
  map<Int_t, map<Int_t, map<UInt_t, Long64_t> > > ak4entry;
  map<Double_t, map<Int_t, Long64_t> > ak4entrymc;
  if (_jp_ak4ak8) {
    if (_dt) {
      assert(fChain2);
      fChain2->SetBranchStatus("*",0);
      fChain2->SetBranchStatus("EvtHdr_.mRun",1);
      fChain2->SetBranchStatus("EvtHdr_.mLumi",1);
      fChain2->SetBranchStatus("EvtHdr_.mEvent",1);

      cout << "Mapping AK4 to AK8 events for data" << endl;
      Long64_t nentries5 = fChain2->GetEntriesFast();
      for (Long64_t jentry5=0; jentry5<nentries5;jentry5++) {
        fChain2->GetEntry(jentry5);
        // Check for duplicate entries
        assert(ak4entry[t4_EvtHdr__mRun][t4_EvtHdr__mLumi][t4_EvtHdr__mEvent]==0);
        ak4entry[t4_EvtHdr__mRun][t4_EvtHdr__mLumi][t4_EvtHdr__mEvent] = jentry5;
      } // for jentry5
      cout << "Found mapping for " << ak4entry.size() << " runs" <<endl<<flush;
    }

    if (_mc) {
      assert(fChain2);
      fChain2->SetBranchStatus("*",0);
      fChain2->SetBranchStatus("EvtHdr_.mPthat",1);
      fChain2->SetBranchStatus("EvtHdr_.mINTPU",1);
      fChain2->SetBranchStatus("EvtHdr_.mOOTPUEarly",1);
      fChain2->SetBranchStatus("EvtHdr_.mOOTPULate",1);

      cout << "Mapping AK4 to AK8 events for MC" << endl;
      Long64_t nentries5 = fChain2->GetEntriesFast();
      int ndup(0);
      for (Long64_t jentry5=0; jentry5<nentries5;jentry5++) {

        fChain2->GetEntry(jentry5);
        Int_t npu = t4_EvtHdr__mINTPU + 500*t4_EvtHdr__mOOTPUEarly
          + 250000*t4_EvtHdr__mOOTPULate;
        // Check for duplicate entries
        if (!(ak4entrymc[t4_EvtHdr__mPthat][npu]==0)) {
          cout << "Error: Duplicates in entries "
               << ak4entrymc[t4_EvtHdr__mPthat][npu]
               << " and " << jentry5 << ": "
               << " pThat="<<t4_EvtHdr__mPthat
               << " INTPU="<<t4_EvtHdr__mINTPU
               << " Early="<<t4_EvtHdr__mOOTPUEarly
               << " Late="<<t4_EvtHdr__mOOTPULate
               << endl << flush;
          //assert(ak4entrymc[t4_EvtHdr__mPthat][npu]==0);
          ++ndup;
        } else {
          ak4entrymc[t4_EvtHdr__mPthat][npu] = jentry5;
        }
      } // for jentry5
      cout << "Found mapping for "<<ak4entrymc.size()<<" pThats"<<endl<<flush;
      if (ndup!=0) {
        cout << "Found "<<ndup<<" duplicates, first ones kept"<<endl<<flush;
      }
    } // _mc

    cout << "/nCreating a shuffled list of jackknife removals based on AK4/n";

    // Modification pre-CWR, Sep 27, 2013
    // Randomly remove event from one of the histograms that was not
    // yet left out in the latest cycle (these are tracked in jkmore)
    // ...but how to make sure same event is removed from AK4 and AK8?
    // => Create list here at the beginning together with ak4entry mappin
    TRandom3 *rnd = new TRandom3();
    Long64_t nentries5 = (fChain2 ? fChain2->GetEntriesFast()
                          : fChain->GetEntriesFast());
    for (Long64_t jentry5=0; jentry5<nentries5; jentry5++) {

      int n = 10;//h->hpt_jk.size();
      if (_jkmore.size()==0) { // if jkmore is empty, reset it first
        for (int ijk = 0; ijk != n; ++ijk) _jkmore.push_back(ijk);
      }
      int m = _jkmore.size(); assert(m!=0);
      int mout = rnd->Integer(m);
      int iout = _jkmore[mout];
      _outlist[jentry5] = iout;
      _jkmore.erase(_jkmore.begin()+mout);
      if (jentry5<300) {
        cout << " " << iout;
        if (jentry5%10==9) cout << ",  ";
        if (jentry5%20==19) cout << endl;
      }
    }
    cout << endl;
  } // _jp_ak4ak8

  if (_jp_quick) {

    fChain->SetBranchStatus("*",0);
    if (fChain2) fChain2->SetBranchStatus("*",0);

    // Luminosity calculation
    if (_mc) fChain->SetBranchStatus("EvtHdr_.mPthat",1); // pthat
    if (_mc) fChain->SetBranchStatus("EvtHdr_.mWeight",1); // weight
    if (_dt) fChain->SetBranchStatus("EvtHdr_.mRun",1); // run
    if (_dt) fChain->SetBranchStatus("EvtHdr_.mEvent",1); // evt
    if (_dt) fChain->SetBranchStatus("EvtHdr_.mLumi",1); // lbn

    // Event properties
    fChain->SetBranchStatus("EvtHdr_.mNVtx",1); // npv
    fChain->SetBranchStatus("EvtHdr_.mNVtxGood",1); // npvgood
    fChain->SetBranchStatus("EvtHdr_.mPFRho",1); // rho

    // Jet properties (jtpt, jte, jteta, jty, jtphi etc.)
    fChain->SetBranchStatus("PFJets_",1); // njt
    fChain->SetBranchStatus("PFJets_.P4_*",1); // jtp4*
    fChain->SetBranchStatus("PFJets_.cor_",1); // jtjes
    fChain->SetBranchStatus("PFJets_.area_",1); // jta

    if (_mc) {
      fChain->SetBranchStatus("PFJets_.genP4_*",1); // jtgenp4*
      fChain->SetBranchStatus("PFJets_.genR_",1); // jtgenr
    }

    // Component fractions
    fChain->SetBranchStatus("PFJets_.chf_",1); // jtchf
    //fChain->SetBranchStatus("PFJets_.phf_",1); // jtnef
    fChain->SetBranchStatus("PFJets_.nemf_",1); // jtnef
    fChain->SetBranchStatus("PFJets_.nhf_",1); // jtnhf
    //fChain->SetBranchStatus("PFJets_.elf_",1); // jtcef !!
    fChain->SetBranchStatus("PFJets_.cemf_",1); // jtcef !!
    fChain->SetBranchStatus("PFJets_.muf_",1); // jtmuf !!
    fChain->SetBranchStatus("PFJets_.ncand_",1); // jtn
    fChain->SetBranchStatus("PFJets_.beta_",1); // jtbeta
    fChain->SetBranchStatus("PFJets_.betaStar_",1); // jtbetastar
    fChain->SetBranchStatus("PFJets_.chm_",1); // jtnch
    fChain->SetBranchStatus("PFJets_.phm_",1); // jtnne
    fChain->SetBranchStatus("PFJets_.nhm_",1); // jtnnh
    fChain->SetBranchStatus("PFJets_.elm_",1); // jtnce !!
    fChain->SetBranchStatus("PFJets_.mum_",1); // jtnmu !!
    fChain->SetBranchStatus("PFJets_.tightID_",1); // jtidtight
    fChain->SetBranchStatus("PFJets_.looseID_",1); // jtidloose

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

    if (_mc) {
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

    if (_jp_ak4ak8) {
      cout << "Activating friend branches" << endl;
      if (_dt) {
        fChain2->SetBranchStatus("EvtHdr_.mRun",1); // run
        fChain2->SetBranchStatus("EvtHdr_.mEvent",1); // evt
        fChain2->SetBranchStatus("EvtHdr_.mLumi",1); // lbn
      }
      if (_mc) {
        fChain2->SetBranchStatus("EvtHdr_.mPthat",1);
        fChain2->SetBranchStatus("EvtHdr_.mOOTPUEarly",1);
        fChain2->SetBranchStatus("EvtHdr_.mOOTPULate",1);
        fChain2->SetBranchStatus("EvtHdr_.mINTPU",1);
        fChain2->SetBranchStatus("GenJets_",1); // njt
        fChain2->SetBranchStatus("GenJets_.fCoordinates.fX",1);
        fChain2->SetBranchStatus("GenJets_.fCoordinates.fY",1);
        fChain2->SetBranchStatus("GenJets_.fCoordinates.fZ",1);
        fChain2->SetBranchStatus("GenJets_.fCoordinates.fT",1);
      }
      fChain2->SetBranchStatus("PFJets_",1); // njt
      fChain2->SetBranchStatus("PFJets_.P4_*",1); // jtp4x
      fChain2->SetBranchStatus("PFJets_.cor_",1); // jtjes
      fChain2->SetBranchStatus("PFJets_.area_",1); // jta
      fChain2->SetBranchStatus("PFJets_.tightID_",1); // jtidtight
      fChain2->SetBranchStatus("PFJets_.looseID_",1); // jtidloose
    }
  } else {
    fChain->SetBranchStatus("*",1);
    if (fChain2) fChain2->SetBranchStatus("*",1);
  } // quick/slow

  // Set pointers to branches
  jtp4x = &PFJets__P4__fCoordinates_fX[0];
  jtp4y = &PFJets__P4__fCoordinates_fY[0];
  jtp4z = &PFJets__P4__fCoordinates_fZ[0];
  jtp4t = &PFJets__P4__fCoordinates_fT[0];
  jta = &PFJets__area_[0];
  jtjes = &PFJets__cor_[0];
  jtbeta = &PFJets__beta_[0];
  jtbetastar = &PFJets__betaStar_[0];
  jtidloose = &PFJets__looseID_[0];
  jtidtight = &PFJets__tightID_[0];
  //
  jtgenr = &PFJets__genR_[0];
  jtgenp4x = &PFJets__genP4__fCoordinates_fX[0];
  jtgenp4y = &PFJets__genP4__fCoordinates_fY[0];
  jtgenp4z = &PFJets__genP4__fCoordinates_fZ[0];
  jtgenp4t = &PFJets__genP4__fCoordinates_fT[0];
  //
  jtn = &PFJets__ncand_[0];
  jtnch = &PFJets__chm_[0];
  jtnnh = &PFJets__nhm_[0];
  jtnne = &PFJets__phm_[0];
  jtnce = &PFJets__elm_[0];
  jtnmu = &PFJets__mum_[0];
  jtchf = &PFJets__chf_[0];
  jtnhf = &PFJets__nhf_[0];
  //jtnef = &PFJets__phf_[0];
  jtnef = &PFJets__nemf_[0];
  //jtcef = &PFJets__elf_[0];
  jtcef = &PFJets__cemf_[0];
  jtmuf = &PFJets__muf_[0];
  //
  gen_jtp4x = &GenJets__fCoordinates_fX[0];
  gen_jtp4y = &GenJets__fCoordinates_fY[0];
  gen_jtp4z = &GenJets__fCoordinates_fZ[0];
  gen_jtp4t = &GenJets__fCoordinates_fT[0];

  if (_jp_ak4ak8) {
    t4_jtp4x = &t4_PFJets__P4__fCoordinates_fX[0];
    t4_jtp4y = &t4_PFJets__P4__fCoordinates_fY[0];
    t4_jtp4z = &t4_PFJets__P4__fCoordinates_fZ[0];
    t4_jtp4t = &t4_PFJets__P4__fCoordinates_fT[0];
    t4_jta = &t4_PFJets__area_[0];
    t4_jtjes = &t4_PFJets__cor_[0];
    t4_jtidloose = &t4_PFJets__looseID_[0];
    t4_jtidtight = &t4_PFJets__tightID_[0];
    if (_mc) {
      t4gen_jtp4x = &t4_GenJets__fCoordinates_fX[0];
      t4gen_jtp4y = &t4_GenJets__fCoordinates_fY[0];
      t4gen_jtp4z = &t4_GenJets__fCoordinates_fZ[0];
      t4gen_jtp4t = &t4_GenJets__fCoordinates_fT[0];
    }
  }

  //assert(_jp_algo=="AK4" || _jp_algo=="AK8");
  const char *a = _jp_algo.c_str();
  cout << "\nCONFIGURATION DUMP:" << endl;
  cout << "-------------------" << endl;
  cout << Form("Running over %sPF",a) << endl;
  if(_mc) cout << "Running over MC" << endl;
  if(_dt) cout << "Running over data" << endl;
  cout << (_jp_redoJEC ? "Re-calculating JEC on the fly" :
           "Using stored JEC in the root tuple") << endl;
  cout << (_jp_useIOV ? "Applying" : "Not applying")
       << " time-dependent JEC (IOV)" << endl;
  cout << (_jp_doEras ? "Storing all " : "Not storing")
       << " eras separately" << endl;
  cout << (_jp_doECALveto ? "Vetoing" : "Not vetoing")
       << " jets in bad ECAL towers" << endl;
  cout << (_jp_doCHS ? "Applying" : "Not applying")
       << " CHS through betaStar" << endl;
  cout << endl;

  if (_dt) {
    cout << (_jp_dojson ? "Applying" : "Not applying")
         << " additional JSON selection" << endl;
    cout << (_jp_doRunHistos ? "Storing" : "Not storing")
         << " additional run-level histograms" << endl;
    cout << (_jp_doBasicHistos ? "Storing" : "Not storing")
         << " basic set of histograms" << endl;
  }
  if (_mc) {
    cout << (_jp_reweighPU ? "Reweighing" : "Not reweighing")
         << " pileup profile in MC to data" << endl;
    cout << (_jp_pthatbins ? "Processing pThat binned samples"
             : "Processing \"flat\" samples") << endl;
  }
  if (_jp_ak4ak8) {
    cout << "Adding extra information for AK4/AK8 ratio" << endl;
  }
  cout << endl;

  // Time dependent JEC
  if (_jp_useIOV) {
    bool isdata = _dt;
    iov = new jec::IOV(Form("%sPF",a));
    iov->add("RunA",160431,163869,isdata);
    iov->add("RunB",165088,167913,isdata);
  }

  // Full redoing of JEC
  const char *s;
  const char *p = "CondFormats/JetMETObjects/data/";
  //const char *t = "GR_R_42_V23_";
  string jecgt = _jp_jecgt + "_" + (_dt ? "DATA" : "MC") + "_";
  const char *t = jecgt.c_str();

  _JEC = 0;
  _L1RC = 0;
  {
    cout << "Loading "<<a<<"PF JEC" << endl;
    s = Form("%s%sL1FastJet_%sPF.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1 = new JetCorrectorParameters(s);
    s = Form("%s%sL2Relative_%sPF.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2 = new JetCorrectorParameters(s);
    s = Form("%s%sL3Absolute_%sPF.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l3 = new JetCorrectorParameters(s);
    s = Form("%s%sL2L3Residual_%sPF.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> vpar;
    vpar.push_back(*par_l1);
    vpar.push_back(*par_l2);
    vpar.push_back(*par_l3);
    if (_dt) vpar.push_back(*par_l2l3res);
    _JEC = new FactorizedJetCorrector(vpar);

    // For type-I MET
    s = Form("%s%sL1FastJet_%sPF.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1rc = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> vrc;
    vrc.push_back(*par_l1rc);
    _L1RC = new FactorizedJetCorrector(vrc);

  } // JEC redone
  assert(_JEC);
  assert(_L1RC);

  // Full redoing of JEC
  _JEC_ak4pf = 0;
  if (_jp_ak4ak8) {
    cout << "Loading "<<a<<"PF JEC" << endl;
    s = Form("%s%sL1FastJet_%sPF.txt",p,t,"AK4"); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1 = new JetCorrectorParameters(s);
    s = Form("%s%sL2Relative_%sPF.txt",p,t,"AK4"); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2 = new JetCorrectorParameters(s);
    s = Form("%s%sL3Absolute_%sPF.txt",p,t,"AK4"); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l3 = new JetCorrectorParameters(s);
    s = Form("%s%sL2L3Residual_%sPF.txt",p,t,"AK4"); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> vpar;
    vpar.push_back(*par_l1);
    vpar.push_back(*par_l2);
    vpar.push_back(*par_l3);
    if (_dt) vpar.push_back(*par_l2l3res);
    _JEC_ak4pf = new FactorizedJetCorrector(vpar);

    assert(_JEC_ak4pf);
  } // JEC redone

  _jecUnc = 0;
  if (_dt) {
    s = Form("%s%sUncertainty_%sPFchs.txt",p,t,a);
    cout<<"**"<<s<<endl<<flush;
    _jecUnc = new JetCorrectionUncertainty(s);
  }

  // Set list of triggers
  for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
    _triggers.push_back(_jp_triggers[itrg]);
  }

  // Load latest JSON selection
  if (_dt && _jp_dojson) loadJSON(_jp_json.c_str());

  // Load PU profiles for MC reweighing
  if (_mc && _jp_reweighPU)
    loadPUProfiles(_jp_pudata.c_str(), _jp_pumc.c_str());

  // Load prescale information to patch 76X
  if (_jp_prescalefile!="")
    loadPrescales(_jp_prescalefile.c_str());

   // load ECAL veto file for cleaning data
   if (_jp_doECALveto) loadECALveto(_jp_ecalveto.c_str());

   // Add these runs to the manual veto list
   if (_dt) {
     // Veto list for 38X 36/pb
     //_runveto.insert(142418);
     // _runveto.insert(162765); // rate 30% in all triggers for 1.1/fb
   }

   // load luminosity tables (prescales now stored in event)
   if (_dt && _jp_dolumi) loadLumi(_jp_lumifile.c_str());

   if (_mc) cout << Form("Running on MC produced with %1.3g nb-1 (%ld evts)",
                         1000. * _entries / _jp_xsecMinBias,
                         (long int)_entries) << endl;
   if (_dt) cout << Form("Running on %ld events of data",
                         (long int)_entries) << endl;

   // Initialize histograms for different epochs and DQM selections
   if (_jp_doBasicHistos) {
     initBasics("Standard");
     if (_dt && _jp_doEras) {
       initBasics("RunA");
       initBasics("RunB");
     }
   }

   if (_dt && _jp_doRunHistos) {
     initRunHistos("Runs",0.,3.);
     initRunHistos("RunsBarrel",0.,1.);
     initRunHistos("RunsTransition",1.,2.);
     initRunHistos("RunsEndcap",2.,3.);
   }

   // Report memory usage to avoid malloc problems when writing file
   *ferr << "Beginning Loop() proper:" << endl << flush;
   cout  << "Beginning Loop() proper:" << endl << flush;
   gSystem->GetMemInfo(&info);
   *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
   cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

   // Event loop
   TStopwatch stop;
   stop.Start();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=nskip; jentry<(nentries+nskip);jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      Long64_t jentry5 = jentry;
      if (_jp_ak4ak8 && _dt) {
        jentry5 = ak4entry[EvtHdr__mRun][EvtHdr__mLumi][EvtHdr__mEvent];
        if (jentry5==0 && jentry!=0) continue;
      }
      if (_jp_ak4ak8 && _mc) {
        Int_t npu = EvtHdr__mINTPU + 500*EvtHdr__mOOTPUEarly
          + 250000*EvtHdr__mOOTPULate;
        jentry5 = ak4entrymc[EvtHdr__mPthat][npu];
        if (jentry5==0 && jentry!=0) continue;
      }
      if (fChain2) fChain2->GetEntry(jentry5);

      // For jackknife in MC (and now data as well)
      _entry = jentry5;

      if (jentry%50000==0) cout << "." << flush;

      if (jentry==10000 || jentry==100000 || jentry==1000000 || jentry==5000000){
        cout << endl
             << Form("Processed %ld events (%1.1f%%) in %1.0f sec. ETA:",
                     (long int)jentry, 100.*jentry/ntot,
                     stop.RealTime()) << endl;
        TDatime now; now.Set(now.Convert()+stop.RealTime()*ntot/jentry);
        now.Print();
        stop.Continue();
      }

      // Set auxiliary event variables (jets, triggers later)
      assert(!_jp_pthatbins || hmcweight);
      pthat = EvtHdr__mPthat;
      weight = (_jp_pthatbins ?
                hmcweight->GetBinContent(hmcweight->FindBin(pthat)) :
                EvtHdr__mWeight);
      // START TEMP PATCH ///
      //if (_mc && !_jp_pthatbins && (fabs(weight-1)<1e-4 || weight==0)) {
      //weight = (pthat ? pow(pthat/15.0,-4.5) : 1);
      //}
      // END TEMP PATCH //
      run = EvtHdr__mRun;
      evt = EvtHdr__mEvent;
      lbn = EvtHdr__mLumi;

      trpu = EvtHdr__mTrPu;
      itpu = EvtHdr__mINTPU;
      //itpu = EvtHdr__mINTPU / 6.; // TEMP PATCH
      ootpulate = EvtHdr__mOOTPULate;
      ootpuearly = EvtHdr__mOOTPUEarly;

      if (_dt) {
        trpu = _avgpu[run][lbn];
        if (trpu==0) {
          //trpu = getTruePU(run, lbn);
          int irun = h2mu->GetXaxis()->FindBin(run);
          int ilbn = h2mu->GetYaxis()->FindBin(lbn);
          trpu = h2mu->GetBinContent(irun, ilbn);
        }
      }

      npv = EvtHdr__mNVtx;
      npvgood = EvtHdr__mNVtxGood;
      pvx = EvtHdr__mPVx;
      pvy = EvtHdr__mPVy;
      pvz = EvtHdr__mPVz;
      pvndof = EvtHdr__mPVndof;
      bsx = EvtHdr__mBSx;
      bsy = EvtHdr__mBSy;

      rho = EvtHdr__mPFRho;
      met = PFMet__et_;
      metphi = PFMet__phi_;
      metsumet = PFMet__sumEt_;

      njt = PFJets__;       //assert(njt < kMaxPFJets_);
      gen_njt = GenJets__;
      t4_njt = t4_PFJets__;
      t4gen_njt = t4_GenJets__;

      assert(njt<_njt);
      if (!(njt < _njt)) {
        *ferr << "Array overflow: njt = " << njt
              << " > njtmax=" << _njt << endl;
        cout << "Array overflow: njt = "<< njt
              << " > njtmax=" << _njt << endl;
        cout << flush;
        assert(njt<_njt);
      }

      if (_debug) {
        cout << endl << flush;
        Show(jentry, jentry5);
        cout << endl << endl << flush;

        cout << "***Checking basic event variables are read out:" << endl;
        cout << "isdata = " << _dt << " / ismc = " << _mc << endl;
        cout << "trpu = " << trpu << endl;
        cout << "pthat = " << pthat << endl;
        cout << "weight = " << weight << endl;
        cout << "njt = " << njt << endl;
        cout << "idloose[0] = " << (njt>0 ? jtidloose[0] : -1) << endl;
        cout << "idtight[0] = " << (njt>0 ? jtidtight[0] : -1) << endl;
        cout << "***end basic event variables" << endl;
        cout << endl << flush;

        if (_jp_ak4ak8) {
          cout << "AK4AK8PF: t4_njt="<<t4_njt<<"("<<t4_PFJets__<<")"<<endl;
          cout << "EVENT: " << jentry << ", " << jentry5 << endl;
          cout << "mRun: "<<EvtHdr__mRun<<" vs "<<t4_EvtHdr__mRun << endl;
          cout << "mLumi: "<<EvtHdr__mLumi<<" vs "<<t4_EvtHdr__mLumi << endl;
          cout << "mEvent: "<<EvtHdr__mEvent<<" vs "<<t4_EvtHdr__mEvent << endl;
          cout << endl << flush;
          assert(EvtHdr__mRun == t4_EvtHdr__mRun);
          assert(EvtHdr__mEvent == t4_EvtHdr__mEvent);
        }
      }

      // Check if duplicate
      if (_dt && _jp_checkduplicates) {
        set<int>& events = _duplicates[run][lbn];
        if (events.find(evt)!=events.end()) {
          ++_nbadevts_dup;
          continue;
        }
        events.insert(evt);
      }

      ++cnt["01all"];

      // Check if good run/LS, including JSON selection
      if (_dt && _jp_dojson) {

        // Does the run/LS pass the latest JSON selection?
        if (_json[run][lbn]==0) {
          _badjson.insert(pair<int, int>(run, lbn));
          ++_nbadevts_json;
          continue;
        }
      } // _dt && _jp_dojson

      if (_dt && _jp_dolumi) {
        // Do we have the run listed in the .csv file?
        map<int, map<int, float> >::const_iterator irun = _lums.find(run);
        if (irun==_lums.end()) {
          _badruns.insert(run);
          ++_nbadevts_run;
          continue;
        }
        // Do we have the LS listed in the .csv file?
        map<int, float>::const_iterator ils = irun->second.find(lbn);
        if (ils==irun->second.end()) {
          _badlums.insert(pair<int, int>(run,lbn));
          ++_nbadevts_ls;
          continue;
        }
        // Does the .csv file list a non-zero luminosity?
        if (ils->second==0) {
          _nolums.insert(pair<int, int>(run, lbn));
          ++_nbadevts_lum;
          //continue; // Could be Poisson fluctuation to zero
        }
      } // _dt && _jp_dolumi

      // Do we exercise run veto based on cross section stability?
      if (_runveto.find(run)!=_runveto.end()) {
        ++_nbadevts_veto;
        continue;
      }

      // Keep track of LBNs
      _jt15lums.insert(pair<int, int>(run, lbn));

      // Reset event ID
      _pass = true;

      if (_pass) ++cnt["02ls"];

      // Reject events with no vertex
      pvrho = tools::oplus(pvx, pvy);
      _pass = (_pass && npvgood>0 && pvrho<2.);

      if (_pass) ++cnt["03vtx"];

      // Event cuts against beam backgrounds
      if (_pass && (tools::oplus(pvx-bsx, pvy-bsy)>0.15 ||
                    pvndof<=4 || fabs(pvz) >= 24.)) {
        ++_bscounter_bad;
        _pass = false;
      }
      if (_pass) ++_bscounter_good;
      if (_pass) ++cnt["04bsc"];

      // Event cuts against beam backgrounds
      /* // Commented out on Aug 4 for testing
      if (_pass && ecalveto &&
          ((njt>=1 &&
            ecalveto->GetBinContent(ecalveto->FindBin(jteta[0],jtphi[0]))!=0) ||
           (njt>=2 &&
            ecalveto->GetBinContent(ecalveto->FindBin(jteta[1],jtphi[1]))!=0))){
        ++_ecalcounter_bad;
        _pass = false;
      } // ecal veto
      */
      if (_pass) ++_ecalcounter_good;
      if (_pass) ++cnt["05ecal"];

      // Check rho
      if (_pass && rho>40.) {
        ++_rhocounter_bad;
        _pass = false;
        if (_debug)
          cout << Form("\nrun:ev:ls %d:%d:%d : rho=%1.1f njt=%d npv=%d"
                       " jtpt0=%1.1f sumet=%1.1f met=%1.1f\n",
                       run, lbn, evt, rho, njt, npv,
                       (njt>0 ? jtpt[0] :0.), metsumet, met) << flush;
      }
      if (_pass) ++_rhocounter_good;
      if (_pass) ++cnt["06rho"];


      // Reset prescales (dynamic can change within run)
      for (map<std::string, std::map<int, int> >::iterator it
             = _prescales.begin(); it != _prescales.end(); ++it) {
        it->second[run] = 0;
      }

      // Fill trigger information
      _trigs.clear();

      // Add special "mc" trigger
      if (_mc) _trigs.insert("mc");
      // Simulate other triggers for MC, if so wished
      // (this is slow, though)
      if (_mc && _jp_domctrigsim) {
        for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
          if (njt>0 && jtpt[0]>_jp_trigthr[itrg]) {
            _trigs.insert(_jp_triggers[itrg]);
          }
        }
      } // _jp_domctrigsim

      // For data, check trigger bits
      if (_debug) cout << "TriggerDecision_.size()=="<<TriggerDecision_.size()<<endl<<flush;
      if (_debug) cout << "_availTrigs.size()=="<<_availTrigs.size()<<endl<<flush;
      assert(TriggerDecision_.size() == _availTrigs.size());

      for (unsigned int itrg = 0; itrg != TriggerDecision_.size(); ++itrg) {

        string strg = _availTrigs[itrg];
        bool pass = (TriggerDecision_[itrg]==1) && strg != ""; // -1, 0, 1

        if (pass) {

          // Set prescale from event for now
          if (L1Prescale_[itrg]>0 && HLTPrescale_[itrg]>0) {
            _prescales[strg][run] = L1Prescale_[itrg] * HLTPrescale_[itrg];
          } else {
            cout << "Error for trigger " << strg << " prescales: "
                 << "L1  =" << L1Prescale_[itrg]
                 << "HLT =" << HLTPrescale_[itrg] << endl;
            _prescales[strg][run] = 0;
          }

          // check prescale
          if (_debug) {
            double prescale = _prescales[strg][run];
            if (L1Prescale_[itrg]*HLTPrescale_[itrg]!=prescale) {
              cout << "Trigger " << strg << ", "
                   << "Prescale(txt file) = " << prescale << endl;
              cout << "L1 = " << L1Prescale_[itrg] << ", "
                   << "HLT = " << HLTPrescale_[itrg] << endl;
              assert(false);
            }
          } // debug

          if (_prescales[strg][run]!=0) {
            // Set trigger only if prescale information is known
            _trigs.insert(strg);
          } else {
            // Make sure all info is good!
            cout << "Missing prescale for " << strg
                 << " in run " << run << endl << flush;
          }
        }
      } // for itrg

      ++_totcounter;
      if (_pass) ++_evtcounter;
      if (_trigs.size()!=0 && _pass) ++_trgcounter;
      if (_trigs.size()!=0 && _pass && _dt) ++cnt["07trg"];

      // Retrieve event weight
      _w0 = (_mc ? weight : 1); assert(_w0);
      _w = _w0;

      // Calculate trigger PU weight
      for (unsigned int itrg = 0; itrg != _triggers.size(); ++itrg) {

        const char *t = _triggers[itrg].c_str();
        _wt[t] = 1.;

        // Reweigh in-time pile-up
        if (_mc && _jp_reweighPU) {

          //int k = pudist[t]->FindBin(trpu);
          //double w1 = pudist[t]->GetBinContent(k);
          int k = pudt->FindBin(trpu);
          double w1 = pudt->GetBinContent(k);
          double w2 = pumc->GetBinContent(k);
          Double_t wtrue = (w1==0 || w2==0 ? 1. : w1 / w2);
          _wt[t] *= wtrue;
        }
      } // for itrg
      _wt["mc"] = _wt[_jp_mctrig];

      // check for non-zero PU weight
      if (_pass && _mc && _jp_reweighPU) {
        _pass = (pudt->GetBinContent(pudt->FindBin(trpu))!=0);
      }
      if (_trigs.size()!=0 && _pass && _mc) ++cnt["07puw"];

      // To-do: implement reweighing for k-factor (NLO*NP/LOMC)

      // load correct IOV for JEC
      if (_dt && _jp_useIOV) {
        _JEC = iov->get(run); assert(_JEC);
      }

      // Calculate pT, eta, phi, y, E and uncorrected pT
      // oversmear jets and MET in MC
      double mex = met * cos(metphi);
      double mey = met * sin(metphi);
      for (int i = 0; i != njt; ++i) {

        // Kostas stores UNCORRECTED four-vector
        p4.SetPxPyPzE(jtp4x[i],jtp4y[i],jtp4z[i],jtp4t[i]);
        jtptu[i] = p4.Pt();
        jteu[i] = p4.E();

        // Recalculate JEC
        _JEC->setRho(rho);
        _JEC->setNPV(npvgood);
        _JEC->setJetA(jta[i]);
        _JEC->setJetPt(jtptu[i]);
        _JEC->setJetE(jteu[i]);
        _JEC->setJetEta(p4.Eta());
        jtjesnew[i] = _JEC->getCorrection();

        // Recalculate JEC (again to get subcorrections)
        _JEC->setRho(rho);
        _JEC->setNPV(npvgood);
        _JEC->setJetA(jta[i]);
        _JEC->setJetPt(jtptu[i]);
        _JEC->setJetE(jteu[i]);
        _JEC->setJetEta(p4.Eta());
        //
        vector<float> v = _JEC->getSubCorrections();
        assert((_dt && v.size()==4) || (_mc && v.size()==3));
        double jec_l1 = v[0];
        double jec_l2l3 = v[2]/v[0];
        double jec_res = (_dt ? v[3]/v[2] : 1.);
        jtjes_l1[i] = jec_l1;
        jtjes_l2l3[i] = jec_l2l3;
        jtjes_res[i] = jec_res;
        assert(jtjesnew[i] = v[v.size()-1]);

        // Correct jets
        p4 *= jtjesnew[i];
        jte[i] = p4.E();
        jtpt[i] = p4.Pt();
        jteta[i] = p4.Eta();
        jtphi[i] = p4.Phi();
        jty[i] = p4.Rapidity();

        // Calculate gen level info
        if (_mc) {
          gp4.SetPxPyPzE(jtgenp4x[i],jtgenp4y[i],jtgenp4z[i],jtgenp4t[i]);
          jtgenpt[i] = gp4.Pt();
          jtgeny[i] = gp4.Rapidity();
        }

        // Oversmear MC to match data
        // Results from Matthias Schroeder, JA July 21, 2011:
        // https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=148123
        /* // Turned off on Oct 4 for testing
        if (_mc && jtgenr[i] < 0.25 && jtgenpt[i] > 0. &&
            jtpt[i] > 0.5*jtgenpt[i] && jtpt[i] < 2.0*jtgenpt[i]) {
          double kover = 1.00;
          double x = fabs(jteta[i]);
          // updated to noResJEC on Aug29 (for |eta|>2.3 in particular)
          if (x < 0.5)              kover = 1.062;//1.052;
          if (x >= 0.5 && x < 1.1)  kover = 1.057;//1.057;
          if (x >= 1.1 && x < 1.7)  kover = 1.089;//1.096;
          if (x >= 1.7 && x < 2.3)  kover = 1.127;//1.134;
          if (x >= 2.3 && x < 5.0)  kover = 1.158;//1.288;
          double dpt = (kover - 1) * (jtpt[i] - jtgenpt[i]);
          double kf = 1 + dpt / jtpt[i];
          jte[i] *= kf;
          jtpt[i] *= kf;
          jtptu[i] *= kf;
          // To-do: propagate this to MET also
          mex -= dpt * cos(jtphi[i]) / jtjesnew[i];
          //mey -= dpt * cos(jtphi[i]) / jtjesnew[i]; // BUG!!!
          mey -= dpt * sin(jtphi[i]) / jtjesnew[i]; // BUG!!!
        }
        */

        met = tools::oplus(mex, mey);
        metphi = atan2(mey, mex);
      } // for i

      if (_mc) {
        for (int i = 0; i != gen_njt; ++i) {

          genp4.SetPxPyPzE(gen_jtp4x[i],gen_jtp4y[i],gen_jtp4z[i],gen_jtp4t[i]);
          gen_jtpt[i] = genp4.Pt();
          gen_jteta[i] = genp4.Eta(); // for matching
          gen_jtphi[i] = genp4.Phi(); // for matching
          gen_jty[i] = genp4.Rapidity();
          // for matching
        } // for i

        if (_jp_ak4ak8) {
          for (int i = 0; i != t4gen_njt; ++i) {

            genp4.SetPxPyPzE(t4gen_jtp4x[i],t4gen_jtp4y[i],
                             t4gen_jtp4z[i],t4gen_jtp4t[i]);
            t4gen_jtpt[i] = genp4.Pt();
            t4gen_jteta[i] = genp4.Eta(); // for matching
            t4gen_jtphi[i] = genp4.Phi(); // for matching
            t4gen_jty[i] = genp4.Rapidity();
            // for matching
          } // for i
        } // _jp_ak4ak8
      } // _mc

      // Propagate jec to MET
      double ucx = -mex;
      double ucy = -mey;
      for (int i = 0; i != njt; ++i) {

        // Only use jets with corr. pT>25 GeV to equalize data and MC thresholds
        if (jtpt[i] > _jp_recopt && fabs(jteta[i])<4.7) {

          // Subtract uncorrected jet pT from met, put back corrected
          // Also add RC offset to keep PU isotropic
          // Remember that MET is negative vector sum
          _L1RC->setRho(rho);
          _L1RC->setJetA(jta[i]);
          _L1RC->setJetPt(jtptu[i]);
          _L1RC->setJetE(jteu[i]);
          _L1RC->setJetEta(jteta[i]);
          double l1corr = _L1RC->getCorrection();
          double dpt = jtpt[i] - l1corr*jtptu[i];
          mex -= dpt * cos(jtphi[i]);
          mey -= dpt * sin(jtphi[i]);

          // Keep track of remaining pT in unclustered energy, i.e.
          // subtract jets from -MET to have the non-jet component
          // treat UE and PU underneath jets as unclustered in order
          // to keep the homogeneous
          double ue = 1.068 * jta[i];
          ucx -= (l1corr * jtptu[i] - ue) * cos(jtphi[i]);
          ucy -= (l1corr * jtptu[i] - ue) * sin(jtphi[i]);
        }
      } // for i
      // Type I MET

      met1 = tools::oplus(mex, mey);
      metphi1 = atan2(mey, mex);
      // Correct unclustered energy; jec for 10 GeV jets varies between
      // 1.1-1.22 at |y|<2.5, 2.5-3.0 even goes up to 1.35
      // => assume average correction of about 1.15 needed
      // => did not seem even nearly enough; try 1.5
      // => reduce down to 1.25 (high pT threshold on jets)
      mex -= 0.25*ucx;//0.5*ucx;
      mey -= 0.25*ucy;//0.5*ucy;
      // Type II MET
      met2 = tools::oplus(mex, mey);
      metphi2 = atan2(mey, mex);

      // Repeat for AK8 jets
      if (_jp_ak4ak8) {
        for (int i = 0; i != t4_njt; ++i) {

          // Kostas stores UNCORRECTED four-vector
          pp4.SetPxPyPzE(t4_jtp4x[i],t4_jtp4y[i],t4_jtp4z[i],t4_jtp4t[i]);
          t4_jtptu[i] = pp4.Pt();
          t4_jteu[i] = pp4.E();

          // Recalculate JEC
          _JEC_ak4pf->setRho(rho);
          _JEC_ak4pf->setNPV(npvgood);
          _JEC_ak4pf->setJetA(t4_jta[i]);
          _JEC_ak4pf->setJetPt(t4_jtptu[i]);
          _JEC_ak4pf->setJetE(t4_jteu[i]);
          _JEC_ak4pf->setJetEta(pp4.Eta());
          t4_jtjesnew[i] = _JEC_ak4pf->getCorrection();

          // Correct jets
          pp4 *= t4_jtjesnew[i];
          t4_jte[i] = pp4.E();
          t4_jtpt[i] = pp4.Pt();
          t4_jteta[i] = pp4.Eta();
          t4_jtphi[i] = pp4.Phi();
          t4_jty[i] = pp4.Rapidity();
        } // for i
      } // _jp_ak4ak8

      if (njt!=0 && _pass) ++cnt["08njt"];

      _jetids.resize(njt);
      for (unsigned int i = 0; i != _jetids.size(); ++i) _jetids[i] = true;
      fillJetID(_jetids);

      if (_jp_ak4ak8) {
        _jetids5.resize(t4_njt);
        for (unsigned int i = 0; i != _jetids5.size(); ++i) _jetids5[i] = true;
        fillJetID5(_jetids5);
      }

      if (njt!=0 && _jetids[0] && _pass) ++cnt["09jtid"];

      // Check if overweight PU event
      if (_mc && njt!=0 && _jetids[0] && _pass) {
        _pass = (jtpt[0] < 1.5*jtgenpt[0] || jtpt[0] < 1.5*pthat);
        //_pass = (jtpt[0] < 2.0*jtgenpt[0] || jtpt[0] < 2.0*pthat);
      }
      if (njt!=0 && _jetids[0] && _pass && _mc) ++cnt["10puw"];

      // Here can categorize events into different triggers, epochs,
      // topologies etc.
      // Eta and pT binning are handled in the fillBasic class
      if (_jp_doBasicHistos) {
        fillBasics("Standard");
        if (_dt && _jp_doEras) {
          //if (run>=160431 && run<=163869) fillBasics("RunA");
          //if (run>=165088 && run<=167913) fillBasics("RunB");
        }
      }

      // Run quality checks
      if (_dt && _jp_doRunHistos) {
        fillRunHistos("Runs");
        fillRunHistos("RunsBarrel");
        fillRunHistos("RunsTransition");
        fillRunHistos("RunsEndcap");
      }

      // Report memory usage to avoid malloc problems when writing file
      if (jentry%1000000==0) {
        *ferr << Form("Doing Loop(), %dM events:",
                      int(jentry/1e6 + 0.5)) << endl << flush;
        gSystem->GetMemInfo(&info);
        *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d,"
                      " Stot:%d, SUsed:%d, SFree:%d",
                      info.fMemTotal, info.fMemUsed, info.fMemFree,
                      info.fSwapTotal, info.fSwapUsed, info.fSwapFree)
              << endl << flush;
      } // 1M report

   } // for jentry
   cout << endl;

   // Report memory usage to avoid malloc problems when writing file
   *ferr << "Finished processing " << nentries << " entries:" << endl << flush;
   cout  << "Finished processing " << nentries << " entries:" << endl << flush;
   gSystem->GetMemInfo(&info);
   *ferr <<Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
   cout << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

   if (_jp_doRunHistos)   writeRunHistos(); // before writeBasics!
   if (_jp_doBasicHistos) writeBasics();

   // List bad runs
   cout << "Processed " << _totcounter << " events in total" << endl;
   cout << "Processed " << _trgcounter << " events passing "
        << " basic data quality and trigger cuts" << endl;
   cout << "(out of " << _evtcounter << " passing data quality cuts)" << endl;
   if (_badruns.size()!=0 || _badlums.size()!=0 || _nolums.size()!=0 ||
       _nbadevts_dup!=0 || _nbadevts_json!=0) {
     cout << "Found " << _badruns.size() << " bad runs:";
     for (set<int>::const_iterator it = _badruns.begin();
          it != _badruns.end(); ++it) {
       cout << " " << *it;
     } // for it
     cout << endl;
     cout << "These contained " << _nbadevts_run << " bad events" << endl;
     cout << "Found " << _nbadevts_json << " bad events according to new JSON"
          << (_jp_dojson ? " (events cut)" : "(events not cut)") << endl;
     cout << "Found " << _badlums.size() << " bad LS and "
          << _nolums.size() << " non-normalizable LS in good runs" << endl;
     cout << "These contained " << _nbadevts_ls << " discarded events"
          << " in bad LS and " << _nbadevts_lum << " in non-normalizable LS"
          << endl;
     cout << endl;
     cout << "Found " << _nbadevts_dup << " duplicate events, which were"
          << " properly discarded" << endl;
     cout << "The vetoed runs contained " << _nbadevts_veto
          << " events" << endl;
   } // has badruns
   cout << "Runs not in JetMETTau stream contained " << _nbadevts_stream
        << " events" << endl;

   // Report beam spot cut efficiency
   cout << "Beam spot counter discarded " << _bscounter_bad
        << " events out of " << _bscounter_good
        << " (" << double(_bscounter_bad)/double(_bscounter_good)*100.
        << "%)" << endl;
   cout << "Beam spot expectation is less than 0.5%" << endl;

   // Report ECAL hole veto efficiency
   cout << "ECAL hole veto counter discarded " << _ecalcounter_bad
        << " events out of " << _ecalcounter_good
        << " (" << double(_ecalcounter_bad)/double(_ecalcounter_good)*100.
        << "%)" << endl;
   cout << "ECAL hole expectation is less than 2.6% [=2*57/(60*72)]" << endl;

   // Report rho veto efficiency
   cout << "Rho<40 veto counter discarded " << _rhocounter_bad
        << " events out of " << _rhocounter_good
        << " (" << double(_rhocounter_bad)/double(_rhocounter_good)*100.
        << "%)" << endl;
   cout << "Rho veto expectation is less than 1 ppm" << endl;

   // Report beam halo efficiency
   cout << "Beam halo counter flagged (not discarded)" << _halocounter_bad
        << " events out of " << _halocounter_good
        << " (" << double(_halocounter_bad)/double(_halocounter_good)*100.
        << "%) " << endl;
   cout << "This is after the beam spot constraint" << endl;

   cout << endl;
   for (map<string,int>::const_iterator it = cnt.begin(); it != cnt.end(); ++it)
     cout << Form("%s: %d (%1.1f%%)", it->first.c_str(), it->second,
                  100. * it->second / max(1, cnt["01all"])) << endl;
   cout << endl;

   // Report LS actually used for Jet15U in the analysis
   // (not necessarily containing any Jet15U triggers, though)
   cout << "Reporting JetMETTau LS in fillhistos.json" << endl;
   ofstream fout("fillhistos.json", ios::out);
   for (set<pair<int, int> >::const_iterator it = _jt15lums.begin();
        it != _jt15lums.end(); ++it) {
     fout << it->first << " " << it->second << endl;
   }
   if (_jp_dojson) {
     cout << "Reporting LS marked newly bad in fillhistos.json.bad" << endl;
     ofstream fout2("fillhistos.json.bad", ios::out);
     for (set<pair<int, int> >::const_iterator it = _badjson.begin();
          it != _badjson.end(); ++it) {
       fout2 << it->first << " " << it->second << endl;
     }
   } // _jp_dojson

   stop.Stop();
   cout << "Processing used " << stop.CpuTime() << "s CPU time ("
        << stop.CpuTime()/3600. << "h)" << endl;
   cout << "Processing used " << stop.RealTime() << "s real time ("
        << stop.RealTime()/3600. << "h)" << endl;
   cout << endl << endl;

   delete ferr;
}


// Initialize basic histograms for trigger and eta bins
void fillHistos::initBasics(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initBasics("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile :
              new TFile(Form("output-%s-1.root",_type.c_str()), "RECREATE"));
  assert(f && !f->IsZombie());
  //assert(f->mkdir(name.c_str()));
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
  //TDirectory *topdir = gDirectory;
  TDirectory *topdir = f->GetDirectory(name.c_str()); assert(topdir);
  topdir->cd();

  // Rapidity bins + HF + barrel
  double y[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 0., 1.3};
  const int ny = sizeof(y)/sizeof(y[0])-1;

  // define triggers
  vector<string> triggers;
  if (_mc) triggers.push_back("mc");
  //triggers.push_back("jt30");
  //triggers.push_back("jt60");
  //triggers.push_back("jt110");
  //triggers.push_back("jt190");
  //triggers.push_back("jt240");
  //triggers.push_back("jt370");
   for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
     triggers.push_back(_jp_triggers[itrg]);
   }

  // define efficient pT ranges for triggers for control plots
  map<string, pair<double, double> > pt;
  if (_mc) pt["mc"] = pair<double, double>(_jp_recopt, _jp_emax);
//   pt["jt30"] = pair<double, double>(56, 97);
//   pt["jt60"] = pair<double, double>(97, 174);
//   pt["jt110"] = pair<double, double>(174,300);
//   pt["jt190"] = pair<double, double>(300,362);
//   pt["jt240"] = pair<double, double>(362,507);
//   pt["jt370"] = pair<double, double>(507, 2500);
  for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
    string trg = _jp_triggers[itrg];
    double pt1 = _jp_trigranges[itrg][0];
    double pt2 = _jp_trigranges[itrg][1];
    pt[trg] = pair<double, double>(pt1, pt2);
   }

  map<string, double> pttrg;
  if (_mc) pttrg["mc"] = _jp_recopt;
//   pttrg["jt30"] = 30.;
//   pttrg["jt60"] = 60.;
//   pttrg["jt110"] = 110.;
//   pttrg["jt190"] = 190.;
//   pttrg["jt240"] = 240.;
//   pttrg["jt300"] = 300.;
//   pttrg["jt370"] = 370.
  for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
    string trg = _jp_triggers[itrg];
    double pt0 = _jp_trigthr[itrg];
    pttrg[trg] = pt0;
   }

  // Loop over rapidity, trigger bins
  for (int i = 0; i != ny; ++i) {

    if (y[i+1] > y[i]) { // create real bins only

      // subdirectory for rapidity bin
      const char *yname = Form("Eta_%1.1f-%1.1f", y[i], y[i+1]);
      assert(topdir);
      //assert(topdir->mkdir(yname));
      topdir->mkdir(yname);
      assert(topdir->cd(yname));
      //TDirectory *ydir = gDirectory;
      TDirectory *ydir = topdir->GetDirectory(yname); assert(ydir);
      ydir->cd();

      for (unsigned int j = 0; j != triggers.size(); ++j) {

        // subdirectory for trigger
        const char *trg = triggers[j].c_str();
        assert(ydir);
        //assert(ydir->mkdir(trg));
        ydir->mkdir(trg);
        assert(ydir->cd(trg));
        //TDirectory *dir = gDirectory;
        TDirectory *dir = ydir->GetDirectory(trg); assert(dir);
        dir->cd();

        // Initialize and store
        assert(dir);
        basicHistos *h = new basicHistos(dir, trg, "", y[i], y[i+1], pttrg[trg],
                                         pt[trg].first, pt[trg].second,
                                         triggers[j]=="mc", dofriends, _jp_ak4ak8);
        _histos[name].push_back(h);
      } // for j
    } // real bin
  } // for i

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initBasics("<<name<<") finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initBasic


// Loop over basic histogram containers to fill all
void fillHistos::fillBasics(string name)
{
  for (unsigned int i = 0; i != _histos[name].size(); ++i)
    fillBasic(_histos[name][i]);
}


// Fill basic histograms after applying pt, y cuts
void fillHistos::fillBasic(basicHistos *h)
{
  assert(h);
  h->hpttmp->Reset();
  h->hpt_tmp->Reset();
  if (h->ismc) {
    assert(h->hpt_g0_tmp);
    h->hpt_g0_tmp->Reset();
  }
  if (h->ak4ak8) {
    h->hpt_tmp_ak4pf->Reset();
  }
  if (h->ak4ak8 && h->ismc) {
    assert(h->hpt_g0_tmp_ak4pf);
    h->hpt_g0_tmp_ak4pf->Reset();
  }

  _w = _w0 * _wt[h->trigname]; assert(_w);

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  // Luminosity information
  if (_dt && h->lums[run][lbn]==0) {
    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    if (lum==0) {
      // cerr << "Run " << run << " LS " << lbn << " empty!" << endl;
      // apparently these can be Poisson fluctuations in the .csv file
    }

    double prescale(0);
    map<int, int>::const_iterator ip = _prescales[h->trigname].find(run);
    if (ip==_prescales[h->trigname].end()) {
      if (fired) {
        *ferr << "No prescale info for trigger " << h->trigname
              << " in run " << run << "!" << endl << flush;
        assert(false);
      }
    }
    else
      prescale = ip->second;

    if (prescale==0 && fired) {
      *ferr << "Prescale zero for trigger " << h->trigname
            << " in run " << run << "!" << endl << flush;
      prescale = 1.;
      assert(false);
    }

    h->lumsum += (prescale ? lum / prescale : 0.);
    h->lumsum2 += (prescale ? lum2 / prescale : 0.);
    h->lums[run][lbn] = 1;

    h->hlumi_vstrpu->Fill(trpu, prescale ? lum / prescale : 0.);
  }
  // For MC vs truePU
  if (_mc)
    h->hlumi_vstrpu->Fill(trpu, _w);

  if (_debug) {
    if (h == _histos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (set<string>::iterator it = _trigs.begin();
          it != _trigs.end(); ++it) {
        cout << *it << ", ";
      }
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired
  if (!fired) return;

  // Check if event looks good
  bool evtid = (met <= 0.4 * metsumet || met < 45.); // QCD-11-004

  if (_debug) cout << Form("Subdirectory Eta_%1.1f-%1.1f/%s",
      h->ymin,h->ymax,h->trigname.c_str()) << endl;
  if (_debug) cout << "Calculate and fill dijet mass" << endl << flush;

  if (h->ismc) h->hpthat->Fill(pthat, _w);
  if (h->ismc) h->hpthatnlo->Fill(pthat);

  if (njt>=2) { // Calculate and fill dijet mass

    // Find leading jets (residual JEC may change ordering)
    map<double, int> ptorder;
    for (int i = 0; i != njt; ++i) {
      double idx = -jtpt[i]; // note minus
      while (ptorder.find(idx) != ptorder.end()) idx += 1e-5*jteta[i];
      assert(ptorder.find(idx)==ptorder.end());
      ptorder[idx] = i;
    }
    int i0 = (ptorder.begin())->second;
    int i1 = (++ptorder.begin())->second;

    //TLorentzVector j1, j2; // use class variables, save instantiation
    j1.SetPtEtaPhiE(jtpt[i0],jteta[i0],jtphi[i0],jte[i0]);
    j2.SetPtEtaPhiE(jtpt[i1],jteta[i1],jtphi[i1],jte[i1]);
    double djmass = (j1+j2).M();
    double ymaxdj = max(fabs(jty[i0]),fabs(jty[i1]));
    bool goodmass = (jtpt[i0]>30. && jtpt[i1]>30.);
    if (evtid && goodmass && _jetids[i0] && _jetids[i1] &&
        ymaxdj >= h->ymin && ymaxdj < h->ymax) {
      assert(h->hdjmass);
      h->hdjmass->Fill(djmass, _w);
      assert(h->hdjmass0);
      h->hdjmass0->Fill(djmass, _w);
      assert(h->pdjmass_ptratio);
      h->pdjmass_ptratio->Fill(djmass, j1.Pt()/j2.Pt(), _w);
      assert(h->pdjmass0_ptratio);
      h->pdjmass0_ptratio->Fill(djmass, j1.Pt()/j2.Pt(), _w);
    }

  } // dijet mass

  if (_debug) cout << "Calculate and fill dijet balance" << endl << flush;

  // Calculate and fill dijet balance histograms
  if (njt>=2 && jtpt[0]>0 && jtpt[1]>0) {
    // Could reorder in pT, but nah (for now)
    int i0 = 0;
    int i1 = 1;
    int i2 = (njt>=3 ? 2 : -1);
    int iref = (fabs(jty[i0]) < fabs(jty[i1]) ? i0 : i1);
    int iprobe = (iref==i0 ? i1 : i0);
    double yref = fabs(jty[iref]);
    double yprobe = fabs(jty[iprobe]);
    double pt1 = jtpt[iref];
    double pt2 = jtpt[iprobe];
    double pt3 = (i2>-1 ? jtpt[i2] : 0.);
    double ptave = 0.5 * (pt2 + pt1); assert(ptave);
    pt3 = (pt3>_jp_recopt ? pt3 / ptave : 0.);
    double asymm = (pt2 - pt1) / (pt2 + pt1);
    double mpf = 1 + met2*cos(delta_phi(metphi2,jtphi[iref]))/ptave;
    //
    double pt3tp = (pt3>_jp_recopt ? pt3 / pt1 : 0.);
    double asymmtp = (pt2 - pt1) / pt1;
    double mpftp = 1 + met2*cos(delta_phi(metphi2,jtphi[iref])) / pt1;

    if (evtid && _jetids[iref] && _jetids[iprobe] &&
        pt1>_jp_recopt && pt2>_jp_recopt
        && delta_phi(jtphi[iref],jtphi[iprobe]) > 2.8 &&
        yprobe >= h->ymin && yprobe < h->ymax && yref < 1.3) {
      assert(h->hdjasymm);
      h->hdjasymm->Fill(ptave, pt3, asymm, _w);
      assert(h->hdjmpf);
      h->hdjmpf->Fill(ptave, pt3, mpf-1, _w);
      //
      assert(h->hdjasymmtp);
      h->hdjasymmtp->Fill(pt1, pt3tp, asymmtp, _w);
      assert(h->hdjmpftp);
      h->hdjmpftp->Fill(pt1, pt3tp, mpftp-1, _w);
    } // first combo
    // If it so happens that both jets at |y|<1.3, fill also the other combo
    if (evtid && _jetids[iref] && _jetids[iprobe] &&
        pt1>_jp_recopt && pt2>_jp_recopt
        && delta_phi(jtphi[iref],jtphi[iprobe]) > 2.8 &&
        yref >= h->ymin && yref < h->ymax && yprobe < 1.3) {
      double mpf2 = 1 + met2*cos(delta_phi(metphi2,jtphi[iprobe]))/ptave;
      double asymm2tp = (pt1 - pt2) / pt2;
      double mpf2tp = 1 + met2*cos(delta_phi(metphi2,jtphi[iprobe])) / pt2;
      assert(h->hdjasymm);
      h->hdjasymm->Fill(ptave, pt3, -asymm, _w);
      assert(h->hdjmpf);
      h->hdjmpf->Fill(ptave, pt3, mpf2-1, _w);
      //
      assert(h->hdjasymmtp);
      h->hdjasymmtp->Fill(pt2, pt3tp, asymm2tp, _w);
      assert(h->hdjmpftp);
      h->hdjmpftp->Fill(pt2, pt3tp, mpf2tp-1, _w);
    } // second combo
  }

  // Fill jet pT ratios vs nvtxgood (pile-up)
  bool has2 = (njt>=2 && jtpt[1] > _jp_recopt &&
      fabs(jty[1])>=h->ymin && fabs(jty[1])<h->ymax);
  bool has3 = (njt>=3 && jtpt[2] > _jp_recopt &&
      fabs(jty[2])>=h->ymin && fabs(jty[2])<h->ymax
      && jtpt[1] > 0.70 * jtpt[0]);
  bool has32 = (has3 && fabs(jty[1]) < 1.3);
  if (_pass && evtid && _jetids[0] && jtpt[0]>=h->ptmin && jtpt[0]<h->ptmax &&
      fabs(jty[0]) < 1.3) {

    h->hr21->Fill(has2 ? jtpt[1] / jtpt[0] : 0.);
    h->hr31->Fill(has3 ? jtpt[2] / jtpt[0] : 0.);
    h->hr32->Fill(has32 ? jtpt[2] / jtpt[1] : 0.);
    if (has2) h->pr21->Fill(npvgood, has2 ? jtpt[1] / jtpt[0] : 0.);
    if (has3) h->pr31->Fill(npvgood, has3 ? jtpt[2] / jtpt[0] : 0.);
    if (has32) h->pr32->Fill(npvgood, has3 ? jtpt[2] / jtpt[1] : 0.);
    h->px21->Fill(npvgood, has2 ? 1 : 0);
    h->px31->Fill(npvgood, has3 ? 1 : 0);
    h->px32->Fill(npvgood, has32 ? 1 : 0);
  }


  if (_debug) cout << "Entering jet loop" << endl << flush;


  for (int i = 0; i != njt && _pass; ++i) {

    if (_debug) cout << "Loop over jet " << i << "/" << njt << endl << flush;

    // adapt variable names from different trees
    double pt = jtpt[i];
    double eta = jteta[i];
    double energy = jte[i];
    double mass = sqrt(fabs(pow(energy,2) - pow(pt*cosh(eta),2)));
    double y = jty[i];
    double phi = jtphi[i];
    double jec = jtjesnew[i];
    bool id = _jetids[i];

    double jec2 = jtjesnew[i]/jtjes[i];

    // Tag-and-probe for composition:
    // tag in barrel and fires trigger, probe in eta bin unbiased
    // only two leading jets back-to-back, third has less than 0.3*tag pT
    if (i<2 && njt>=2 && pt>_jp_recopt &&
        fabs(y) >= h->ymin && fabs(y) < h->ymax){

      int iref = (i==0 ? 1 : 0);
      double yref = jty[iref];
      double ptref = jtpt[iref];
      double dphi = delta_phi(phi, jtphi[iref]);
      double pt3 = (njt>=3 ? jtpt[2] : 0.);

      int ipf4(-1);
      double dr4min(999.);
      double ptprobepf4(0.);
      if (h->ak4ak8) {
        for (int j = 0; j != t4_njt; ++j) {
          double dr = tools::oplus(delta_phi(t4_jtphi[j], jtphi[i]),
              fabs(t4_jteta[j] - jteta[i]));
          if (dr < dr4min) {
            ipf4 = j;
            dr4min = dr;
          }
        }
        ptprobepf4 = (dr4min<0.8 ? t4_jtpt[ipf4] : 0.);
      }

      if (evtid && id && _jetids[iref] &&
          fabs(yref) < 1.3 && dphi > 2.7 && pt3 < 0.3*ptref) {

        assert(h->pncandtp);
        h->pncandtp->Fill(ptref, jtn[i], _w);
        assert(h->pnchtp);
        h->pnchtp->Fill(ptref, jtnch[i], _w);
        assert(h->pnnetp);
        h->pnnetp->Fill(ptref, jtnne[i], _w);
        assert(h->pnnhtp);
        h->pnnhtp->Fill(ptref, jtnnh[i], _w);
        assert(h->pncetp);
        h->pncetp->Fill(ptref, jtnce[i], _w);
        assert(h->pnmutp);
        h->pnmutp->Fill(ptref, jtnmu[i], _w);
        //
        assert(h->pchftp);
        h->pchftp->Fill(ptref, jtchf[i], _w);
        assert(h->pneftp);
        h->pneftp->Fill(ptref, jtnef[i], _w);
        assert(h->pnhftp);
        h->pnhftp->Fill(ptref, jtnhf[i], _w);
        assert(h->pceftp);
        h->pceftp->Fill(ptref, jtcef[i], _w);
        assert(h->pmuftp);
        h->pmuftp->Fill(ptref, jtmuf[i], _w);
        assert(h->pbetatp);
        h->pbetatp->Fill(ptref, jtbeta[i], _w);
        assert(h->pbetastartp);
        h->pbetastartp->Fill(ptref, jtbetastar[i], _w);
        //
        assert(h->pchftp2);
        h->pchftp2->Fill(ptref, jtchf2[i], _w);
        assert(h->pneftp2);
        h->pneftp2->Fill(ptref, jtnef2[i], _w);
        assert(h->pnhftp2);
        h->pnhftp2->Fill(ptref, jtnhf2[i], _w);
        assert(h->pceftp2);
        h->pceftp2->Fill(ptref, jtcef2[i], _w);
        assert(h->pmuftp2);
        h->pmuftp2->Fill(ptref, jtmuf2[i], _w);
        //
        if (dr4min<0.4) h->pak4ak8tp_50->Fill(ptref, ptprobepf4/jtpt[i]);
        if (dr4min<0.2) h->pak4ak8tp_25->Fill(ptref, ptprobepf4/jtpt[i]);

        if (ptref > h->ptmin && ptref < h->ptmax) {

          h->hncandtp->Fill(jtn[i], _w);
          h->hnchtp->Fill(jtnch[i], _w);
          h->hnnetp->Fill(jtnne[i], _w);
          h->hnnhtp->Fill(jtnnh[i], _w);
          h->hncetp->Fill(jtnce[i], _w);
          h->hnmutp->Fill(jtnmu[i], _w);
          h->hchftp->Fill(jtchf[i], _w);
          h->hneftp->Fill(jtnef[i], _w);
          h->hnhftp->Fill(jtnhf[i], _w);
          h->hceftp->Fill(jtcef[i], _w);
          h->hmuftp->Fill(jtmuf[i], _w);
          h->hbetatp->Fill(jtbeta[i], _w);
          h->hbetastartp->Fill(jtbetastar[i], _w);
          //
          if (dr4min<0.4) h->hak4ak8tp->Fill(ptprobepf4/jtpt[i]);
          //
          assert(h->pncandtp_vsnpv);
          h->pncandtp_vsnpv->Fill(npvgood, jtn[i], _w);
          assert(h->pnchtp_vsnpv);
          h->pnchtp_vsnpv->Fill(npvgood, jtnch[i], _w);
          assert(h->pnnetp_vsnpv);
          h->pnnetp_vsnpv->Fill(npvgood, jtnne[i], _w);
          assert(h->pnnhtp_vsnpv);
          h->pnnhtp_vsnpv->Fill(npvgood, jtnnh[i], _w);
          assert(h->pncetp_vsnpv);
          h->pncetp_vsnpv->Fill(npvgood, jtnce[i], _w);
          assert(h->pnmutp_vsnpv);
          h->pnmutp_vsnpv->Fill(npvgood, jtnmu[i], _w);
          //
          assert(h->pchftp_vsnpv);
          h->pchftp_vsnpv->Fill(npvgood, jtchf[i], _w);
          assert(h->pneftp_vsnpv);
          h->pneftp_vsnpv->Fill(npvgood, jtnef[i], _w);
          assert(h->pnhftp_vsnpv);
          h->pnhftp_vsnpv->Fill(npvgood, jtnhf[i], _w);
          assert(h->pceftp_vsnpv);
          h->pceftp_vsnpv->Fill(npvgood, jtcef[i], _w);
          assert(h->pmuftp_vsnpv);
          h->pmuftp_vsnpv->Fill(npvgood, jtmuf[i], _w);
          assert(h->pbetatp_vsnpv);
          h->pbetatp_vsnpv->Fill(npvgood, jtbeta[i], _w);
          assert(h->pbetastartp_vsnpv);
          h->pbetastartp_vsnpv->Fill(npvgood, jtbetastar[i], _w);
          //
          assert(h->pchftp_vstrpu);
          h->pchftp_vstrpu->Fill(trpu, jtchf2[i], _w);
          assert(h->pneftp_vstrpu);
          h->pneftp_vstrpu->Fill(trpu, jtnef2[i], _w);
          assert(h->pnhftp_vstrpu);
          h->pnhftp_vstrpu->Fill(trpu, jtnhf2[i], _w);
          assert(h->pceftp_vstrpu);
          h->pceftp_vstrpu->Fill(trpu, jtcef2[i], _w);
          assert(h->pmuftp_vstrpu);
          h->pmuftp_vstrpu->Fill(trpu, jtmuf2[i], _w);
          assert(h->pbetatp_vstrpu);
          h->pbetatp_vstrpu->Fill(trpu, jtbeta[i], _w);
          assert(h->pbetastartp_vstrpu);
          h->pbetastartp_vstrpu->Fill(trpu, jtbetastar[i], _w);
        }
      } // dijet system
    } // tag-and-probe

    // Check effect of ID cuts
    //if (pt>_jp_recopt && fabs(y) >= h->ymin && fabs(y) < h->ymax) {
    if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {

      if (_debug) {
        cout << "..." << h->trigname << " | " << " index " << i << "/" << njt
          << " jet pt: " << pt << " y : " << y
          << " id " << id << " jec: " << jec << endl;
        cout << "...evt id: " << evtid << " weight: " << _w
          << " met: " << met << " metsumet: " << metsumet << endl;
      }

      assert(h->hpt_noid);
      h->hpt_noid->Fill(pt, _w);
      assert(h->hpt_nojetid);
      if (evtid) h->hpt_nojetid->Fill(pt, _w);
      assert(h->hpt_noevtid);
      if (id)    h->hpt_noevtid->Fill(pt, _w);
      // Same versus generator pT as MC extra
      // to decouple efficiency from JEC and JER
      if (h->ismc) {
        h->hpt_noid_g->Fill(jtgenpt[i], _w);
        if (evtid) h->hpt_nojetid_g->Fill(jtgenpt[i], _w);
        if (id)    h->hpt_noevtid_g->Fill(jtgenpt[i], _w);
      }
    } // ID cuts

    // Check effect of reco y vs gen y binning
    if (h->ismc) {
      double ygen = jtgeny[i]; // use jtgeny, if available
      // GenJets matched to good reco jets in good events
      if (evtid && id && pt>_jp_recopt && jtgenr[i] < 0.25 &&
          fabs(ygen) >= h->ymin && fabs(ygen) < h->ymax) {
        h->hpt_gg->Fill(jtgenpt[i], _w);
      }
      // GenJets matched to any reco jets in any events
      if (pt>_jp_recopt && jtgenr[i] < 0.25 &&
          fabs(ygen) >= h->ymin && fabs(ygen) < h->ymax) {
        h->hpt_gg0->Fill(jtgenpt[i], _w);
      }
    }

    // Debugging JEC
    /*
       if (h->ismc && h->ymin==0 && h->ymax==0.5) {

       int igen = -1;
    // Find matching genjet
    double drmin = 999;
    for (int j = 0; j != gen_njt; ++j) {
    double dr = tools::oplus(delta_phi(gen_jtphi[j], jtphi[i]),
    fabs(gen_jteta[j] - jteta[i]));
    if (dr < drmin) {
    igen = j;
    drmin = dr;
    }
    }

    if (npvgood==1 && fabs(jteta[i])<1.3 &&
    jtgenpt[i]>=80 && jtgenpt[i]<120.) {
    cout << "jet " << i << " eta = " << jteta[i]
    << " y " << jty[i] << " ygen = " << jtgeny[i]
    << " pt = " << jtpt[i]
    << " ptuncorr = " << jtptu[i]
    << " (ratio = " << jtpt[i]/jtptu[i]
    << " , jec = " << jtjesnew[i]
    << " , old = " << jtjes[i] << ")"
    << " npvgood = " << npvgood
    << " jtgenpt = " << jtgenpt[i]
    << " dr = " << jtgenr[i]
    << " genpt = " << (igen!=-1 ? gen_jtpt[igen] : -1)
    << " drmin = " << drmin
    << " gen/reco = " << jtgenpt[i]/jtpt[i] << endl;
    }
    } // debugging JEC
     */

    // calculate efficiencies and fill histograms
    if (evtid && id && pt>_jp_recopt &&
        fabs(y) >= h->ymin && fabs(y) < h->ymax) {

      if (_debug) cout << "..jec uncertainty" << endl << flush;

      // Get JEC uncertainty
      double unc = 0.01; // default for MC
      if (_jecUnc) {
        _jecUnc->setJetEta(eta);
        _jecUnc->setJetPt(pt);
        unc = _jecUnc->getUncertainty(true);
        //_jecUnc2->Rjet(pt, unc); // use Fall10 absolute scale uncertainty
      }

      // retrieve event-wide variables
      double dphi = (njt>1 ? delta_phi(jtphi[0], jtphi[1]) : 0.);
      double dpt = (njt>1 ? fabs(jtpt[0]-jtpt[1])/(jtpt[0]+jtpt[1]) : 0.999);
      //double met = this->met;
      //double metphi = this->metphi;
      double sumet = this->metsumet;

      // calculate and/or retrieve efficiencies
      double ideff = 1.;
      double vtxeff = 1.;
      double dqmeff = 1.;
      double trigeff = 1.;
      double eff = ideff * vtxeff * dqmeff * trigeff;

      if (_debug) cout << "..raw spectrum" << endl << flush;

      // For trigger efficiency
      if (h->ismc) {
        h->hpt_jt30->Fill(pt, _w0 * _wt["jt30"]);
        h->hpt_jt60->Fill(pt, _w0 * _wt["jt60"]);
        //h->hpt_jt80->Fill(pt, _w0 * _wt["jt80"]);
        h->hpt_jt110->Fill(pt, _w0 * _wt["jt110"]);
        //h->hpt_jt150->Fill(pt, _w0 * _wt["jt150"]);
        h->hpt_jt190->Fill(pt, _w0 * _wt["jt190"]);
        h->hpt_jt240->Fill(pt, _w0 * _wt["jt240"]);
        h->hpt_jt300->Fill(pt, _w0 * _wt["jt300"]);
        h->hpt_jt370->Fill(pt, _w0 * _wt["jt370"]);
        //
        h->hpt0_jt30->Fill(pt, _w0 * _wt["jt30"]);
        h->hpt0_jt60->Fill(pt, _w0 * _wt["jt60"]);
        //h->hpt0_jt80->Fill(pt, _w0 * _wt["jt80"]);
        h->hpt0_jt110->Fill(pt, _w0 * _wt["jt110"]);
        //h->hpt0_jt150->Fill(pt, _w0 * _wt["jt150"]);
        h->hpt0_jt190->Fill(pt, _w0 * _wt["jt190"]);
        h->hpt0_jt240->Fill(pt, _w0 * _wt["jt240"]);
        h->hpt0_jt300->Fill(pt, _w0 * _wt["jt300"]);
        h->hpt0_jt370->Fill(pt, _w0 * _wt["jt370"]);
      }

      // raw spectrum
      assert(h->hpt);
      h->hpt->Fill(pt, _w);
      h->hpt_tmp->Fill(pt); // Event statistics
      assert(h->hpt_pre);
      if (_dt) h->hpt_pre->Fill(pt, _w*_prescales[h->trigname][run]);
      if (_mc) h->hpt_pre->Fill(pt, _w0*_wt["mc"]);
      assert(h->hpt0);
      h->hpt0->Fill(pt, _w);
      //if (y>=0) {
      //h->hpt_plus->Fill(pt, _w);
      //h->hpt0_plus->Fill(pt, _w);
      //h->hpt_plus_38x->Fill(pt, _w);
      //}
      //if (y<0) {
      //h->hpt_minus->Fill(pt, _w);
      //h->hpt0_minus->Fill(pt, _w);
      //h->hpt_minus_38x->Fill(pt, _w);
      //}
      // Do proper event statistics
      if (h->hpttmp->GetBinContent(h->hpttmp->FindBin(pt))==0)
        h->hptevt->Fill(pt, _w);
      h->hpttmp->Fill(pt);

      // delete-m jackknife for AK4/AK8 statistics
      /*
         if (_dt) {
         int n = h->hpt_jk.size();
         for (int iout = 0; iout != n; ++iout) {
      //if (evt%n!=iout) h->hpt_jk[iout]->Fill(pt, _w);
      if (_entry%n!=iout) h->hpt_jk[iout]->Fill(pt, _w); // March 28, 2013
      }
      }
      if (_mc) {
      int n = h->hpt_jk.size();
      for (int iout = 0; iout != n; ++iout) {
      if (_entry%n!=iout) h->hpt_jk[iout]->Fill(pt, _w);
      }
      }
       */
      //if (_dt || _mc) {
      if (_jp_ak4ak8) {
        // Modification pre-CWR, Sep 27, 2013
        // Randomly remove event from one of the histograms that was not
        // yet left out in the latest cycle (these are tracked in jkmore)
        // ...but how to make sure same event is removed from AK4 and AK8?
        // => move into beginning to pre-create removal list
        int n = h->hpt_jk.size(); assert(n==10);
        //if (_jkmore.size()==0) { // if jkmore is empty, reset it first
        //for (int ijk = 0; ijk != n; ++ijk) _jkmore.push_back(ijk);
        //}
        ///int m = _jkmore.size(); assert(m!=0);
        //int mout = gRandom3->Integer(m);
        //int iout = _jkmore[mout];
        int iout = _outlist[_entry];
        for (int ijk = 0; ijk != n; ++ijk) {
          if (ijk!=iout) h->hpt_jk[ijk]->Fill(pt, _w);
        }
        h->h2jk->Fill(_entry%n, iout, _w);
        //_jkmore.erase(_jkmore.begin()+mout);
      }

      // leading and non-leading jets
      assert(h->hpt1);
      if (i==0) { h->hpt1->Fill(pt, _w); }
      assert(h->hpt2);
      if (i==1) { h->hpt2->Fill(pt, _w); }
      assert(h->hpt3);
      if (i==2) { h->hpt3->Fill(pt, _w); }

      if (_debug) cout << "..basic properties" << endl << flush;

      // basic properties
      h->ppt->Fill(pt, pt, _w); assert(h->ppt);
      h->pmass->Fill(pt, mass/energy, _w); assert(h->pmass);
      h->pjec->Fill(pt, jec, _w); assert(h->pjec);
      h->pjec2->Fill(pt, jec2, _w); assert(h->pjec2);
      h->punc->Fill(pt, unc, _w); assert(h->punc);

      // JEC monitoring
      h->pjec_l1->Fill(pt, jtjes_l1[i], _w); assert(h->pjec_l1);
      h->pjec_l2l3->Fill(pt, jtjes_l2l3[i], _w); assert(h->pjec_l2l3);
      h->pjec_res->Fill(pt, jtjes_res[i], _w); assert(h->pjec_res);

      // Pile-up information
      h->pa->Fill(pt, jta[i], _w);
      h->ptrpu->Fill(pt, trpu, _w);
      h->prho->Fill(pt, rho, _w);
      h->pnpv->Fill(pt, npvgood, _w);
      h->pnpvall->Fill(pt, npv, _w);
      if (pt >= h->ptmin && pt < h->ptmax) {
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
      }

      // efficiencies
      assert(h->peff);
      h->peff->Fill(pt, eff, _w);
      assert(h->pideff);
      h->pideff->Fill(pt, ideff, _w);
      assert(h->pvtxeff);
      h->pvtxeff->Fill(pt, vtxeff, _w);
      assert(h->pdqmeff);
      h->pdqmeff->Fill(pt, dqmeff, _w);

      if (_debug) cout << "..control plots of components" << endl << flush;

      // control plots of jet components (JEC)
      assert(h->pncand);
      h->pncand->Fill(pt, jtn[i], _w);
      assert(h->pnch);
      h->pnch->Fill(pt, jtnch[i], _w);
      assert(h->pnne);
      h->pnne->Fill(pt, jtnne[i], _w);
      assert(h->pnnh);
      h->pnnh->Fill(pt, jtnnh[i], _w);
      assert(h->pnce);
      h->pnce->Fill(pt, jtnce[i], _w);
      assert(h->pnmu);
      h->pnmu->Fill(pt, jtnmu[i], _w);
      //
      assert(h->pchf);
      h->pchf->Fill(pt, jtchf[i], _w);
      assert(h->pnef);
      h->pnef->Fill(pt, jtnef[i], _w);
      assert(h->pnhf);
      h->pnhf->Fill(pt, jtnhf[i], _w);
      assert(h->pcef);
      h->pcef->Fill(pt, jtcef[i], _w);
      assert(h->pmuf);
      h->pmuf->Fill(pt, jtmuf[i], _w);
      assert(h->pbeta);
      h->pbeta->Fill(pt, jtbeta[i], _w);
      assert(h->pbetastar);
      h->pbetastar->Fill(pt, jtbetastar[i], _w);

      // Find AK4 match for AK4/AK8 ratio studies
      double ptpf4(0.), dr4min(999.);
      if (h->ak4ak8) {

        int ipf4 = -1;
        dr4min = 999.;
        for (int j = 0; j != t4_njt; ++j) {
          double dr = tools::oplus(delta_phi(t4_jtphi[j], phi),
                                   fabs(t4_jteta[j] - eta));
          if (dr < dr4min) {
            ipf4 = j;
            dr4min = dr;
          }
        }
        ptpf4 = (dr4min<0.4 ? t4_jtpt[ipf4] : 0.);
        if (dr4min<0.4) h->pak4ak8_50->Fill(pt, ptpf4/pt);
        if (dr4min<0.2) h->pak4ak8_25->Fill(pt, ptpf4/pt);

      } // _jp_ak4ak8

      // control plots for topology (JEC)
      if (pt >= h->ptmin && pt < h->ptmax) {
        if (_debug) cout << "..control plots for topology" << endl << flush;

        h->htrpu->Fill(trpu, _w);
        if (h->ismc) {
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
        h->hjet->Fill(pt / sumet, _w);
        h->hmet->Fill(met / sumet, _w);
        h->hmetphi->Fill(delta_phi(metphi, phi), _w);
        // control plots for vertex
        h->hpvndof->Fill(pvndof);
        h->hpvx->Fill(pvx-bsx);
        h->hpvy->Fill(pvy-bsy);
        h->hpvz->Fill(pvz-0.);
        h->hpvr->Fill(tools::oplus(pvx-bsx, pvy-bsy));
        h->hpvrho->Fill(pvrho-tools::oplus(bsx, bsy));
        // closure plots for JEC
        h->hmpf->Fill(1 + met * cos(delta_phi(metphi, phi)) / pt, _w);
        h->hmpf1->Fill(1 + met1 * cos(delta_phi(metphi1, phi)) / pt, _w);
        h->hmpf2->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);
        //
        if ((njt<3 || jtpt[2]<0.15*(jtpt[0]+jtpt[1]) || jtpt[2] < _jp_recopt) &&
            (njt>=2 && jtpt[1]>_jp_recopt && delta_phi(jtphi[0],jtphi[1]) > 2.7))
            h->hmpfy->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);

        // Component fractions
        h->hncand->Fill(jtn[i], _w);
        h->hnch->Fill(jtnch[i], _w);
        h->hnne->Fill(jtnne[i], _w);
        h->hnnh->Fill(jtnnh[i], _w);
        h->hnce->Fill(jtnce[i], _w);
        h->hnmu->Fill(jtnmu[i], _w);
        //
        h->hchf->Fill(jtchf[i], _w);
        h->hnef->Fill(jtnef[i], _w);
        h->hnhf->Fill(jtnhf[i], _w);
        h->hcef->Fill(jtcef[i], _w);
        h->hmuf->Fill(jtmuf[i], _w);
        h->hbeta->Fill(jtbeta[i], _w);
        h->hbetastar->Fill(jtbetastar[i], _w);
        //
        if (dr4min<0.4) h->hak4ak8->Fill(ptpf4/pt);

        h->hyeta->Fill(TMath::Sign(y-eta,y), _w);
        h->hyeta2->Fill(y-eta, _w);
        h->hbetabetastar->Fill(jtbeta[i], jtbetastar[i], _w);
        h->hetaphi->Fill(eta, phi, _w);
      } // within trigger pT range

      int iprobe = i;
      int itag = (iprobe==0 ? 1 : 0);
      double pttag = (njt>=2 ? jtpt[itag] : 0);
      if (iprobe<2 && pttag >= h->ptmin && pttag < h->ptmax &&
          fabs(jty[itag]) < 1.3) {

        if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2]<_jp_recopt) &&
            (njt>=2 && jtpt[1]>_jp_recopt && delta_phi(jtphi[0],jtphi[1])>2.7))
          h->hmpfx->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / pttag, _w);
      }
      double ptave = (njt>=2 ? 0.5*(pt+pttag) : 0);
      if (iprobe<2 && ptave >= h->ptmin && ptave < h->ptmax &&
          fabs(jty[itag]) < 1.3) {

        if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2]<_jp_recopt) &&
            (njt>=2 && jtpt[1]>_jp_recopt && delta_phi(jtphi[0],jtphi[1])>2.7))
          h->hmpfz->Fill(1 + met2 * cos(delta_phi(metphi2, phi)) / ptave, _w);
      }

      // closure plots for JEC
      h->pdpt->Fill(pt, dpt, _w);
      h->pmpf->Fill(pt, 1 + met * cos(delta_phi(metphi, phi)) / pt, _w);
      h->pmpf1->Fill(pt, 1 + met1 * cos(delta_phi(metphi1, phi)) / pt, _w);
      h->pmpf2->Fill(pt, 1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);
      //
      if ((njt<3 || jtpt[2] < 0.15*(jtpt[0]+jtpt[1]) || jtpt[2]<_jp_recopt) &&
          (njt>=2 && jtpt[1]>_jp_recopt && delta_phi(jtphi[0],jtphi[1])>2.7)) {

        h->pmpfx->Fill(pttag, 1 + met2 * cos(delta_phi(metphi2, phi))/pttag, _w);
        h->pmpfy->Fill(pt, 1 + met2 * cos(delta_phi(metphi2, phi))/pt, _w);
        h->pmpfz->Fill(ptave, 1 + met2 * cos(delta_phi(metphi2, phi))/ptave, _w);
      }

      // MC extras
      if (_mc && jtgenr[i]<0.25) {
        h->hpt_gtw->Fill(jtgenpt[i], _w);
      }
      //
      if (h->ismc) {

        if (jtgenr[i]<0.25) {

          //int flv = abs(jtgenflv[i]);
          double ptgen = jtgenpt[i];
          //double r = (jtgenpt[i] ? pt/jtgenpt[i] : 0);
          double r = (ptgen ? pt/ptgen : 0);
          //double resp = (jtjesnew[i] ? r / jtjesnew[i] : 0);
          double dy = (r ? TMath::Sign(jty[i]-jtgeny[i], jtgeny[i]) : 0.);
          h->hpt_r->Fill(pt, _w);
          h->hpt_g->Fill(ptgen, _w);
          h->ppt_r->Fill(pt, pt, _w);
          h->ppt_g->Fill(ptgen, ptgen, _w);

          // Response closure vs NPV
          //if (r) h->p3rvsnpv->Fill(ptgen, jteta[i], npvgood, resp, _w);
          //if (r) h->p3rvsnpvW->Fill(ptgen, fabs(jteta[i]), npvgood, resp, _w);
          if (r) h->p2rvsnpv->Fill(ptgen, npvgood, r, _w);

          // Response closure
          if (r) h->h2r_r->Fill(pt, r, _w);
          if (r) h->h2r_g->Fill(ptgen, r, _w);
          if (r) h->p2r_r->Fill(pt, r, _w);
          if (r) h->p2r_g->Fill(ptgen, r, _w);
          if (r) h->p2r_ruw->Fill(pt, r); // unweighted!
          if (r) h->p2r_guw->Fill(ptgen, r); // unweighted!

          // Rapidity closure
          if (r) h->h2dy_r->Fill(pt, dy, _w);
          if (r) h->h2dy_g->Fill(ptgen, dy, _w);
          if (r) h->p2dy_r->Fill(pt, dy, _w);
          if (r) h->p2dy_g->Fill(ptgen, dy, _w);
          if (r) h->p2dy_ruw->Fill(pt, dy); // unweighted
          if (r) h->p2dy_guw->Fill(ptgen, dy); // unweighted
          if (r) h->pdy_r->Fill(pt, fabs(y), dy, _w);
          if (r) h->pdy_g->Fill(ptgen, fabs(y), dy, _w);
        }
      } // is MC

    } // if id && etabin

    // Fill outside of eta bin also
    if (h->ismc && evtid && id && pt > _jp_recopt &&
        jtgenr[i]<0.25 && jtgenpt[i]!=0 && jtjesnew[i]!=0) {

      double ptgen = jtgenpt[i];
      double r = pt / ptgen;
      double resp = r / jtjesnew[i];

      // Response closure vs NPV
      h->p3rvsnpv->Fill(ptgen, jteta[i], npvgood, resp, _w);
      h->p3rvsnpvW->Fill(ptgen, fabs(jteta[i]), npvgood, resp, _w);
    } // if id && MC
  } // for i

  if (h->ak4ak8) { // _jp_ak4ak8

    bool evtid5 = (met < 0.4 * metsumet || met < 45.);

    for (int i = 0; i != t4_njt; ++i) {

      double y = t4_jty[i];
      if (evtid5 && t4_jtidtight[i] && _pass &&
          fabs(y) >= h->ymin && fabs(y) < h->ymax) {

        double pt = t4_jtpt[i];

        if (pt > _jp_recopt) {

          h->hpt_ak4pf->Fill(pt, _w);
          h->hpt_tmp_ak4pf->Fill(pt); // Event statistics
        } // reco pt
      } // y bin
    } // for i
  } // _jp_ak4ak8


  // Event statistics
  for (int i = 1; i != h->hpt_tmp->GetNbinsX()+1; ++i) {
    if (h->hpt_tmp->GetBinContent(i)!=0) {

      double pt = h->hpt_tmp->GetBinCenter(i);
      int njet = h->hpt_tmp->GetBinContent(i);
      h->hpt_evtcount->Fill(pt);
      h->hpt_evt->Fill(pt, _w);
      h->hpt_jet->Fill(pt, _w*njet);
    }
  } // for i

  //if (_jp_ak4ak8) {
  if (h->ak4ak8) {

    for (int i = 1; i != h->hpt_tmp_ak4pf->GetNbinsX()+1; ++i) {

      if (h->hpt_tmp_ak4pf->GetBinContent(i)!=0) {

        double pt = h->hpt_tmp_ak4pf->GetBinCenter(i);
        int njet = h->hpt_tmp_ak4pf->GetBinContent(i);
        h->hpt_evtcount_ak4pf->Fill(pt);
        h->hpt_evt_ak4pf->Fill(pt, _w);
        h->hpt_jet_ak4pf->Fill(pt, _w*njet);
      }

      // Jet number correlations
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          h->hpt_tmp_ak4pf->GetBinContent(i)!=0) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i);
        int n7 = h->hpt_tmp->GetBinContent(i);
        h->hpt_ak4ak8->Fill(pt, min(n7,2), min(n5,2), _w);
      }
      if (h->hpt_tmp->GetBinContent(i)!=0) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n7 = h->hpt_tmp->GetBinContent(i);
        h->hpt_ak8ak8->Fill(pt, min(n7,2), min(n7,2), _w);
      }
      if (h->hpt_tmp_ak4pf->GetBinContent(i)!=0) {
        double pt = h->hpt_tmp_ak4pf->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i);
        h->hpt_ak4ak4->Fill(pt, min(n5,2), min(n5,2), _w);
      }
      // Jet number correlations for nearby off-diagonals
      // Minus one
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          (i>1 && h->hpt_tmp_ak4pf->GetBinContent(i-1)!=0)) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i-1);
        int n7 = h->hpt_tmp->GetBinContent(i);
        h->hpt_ak4ak8m1->Fill(pt, min(n7,2), min(n5,2), _w);
      }
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          (i>1 && h->hpt_tmp->GetBinContent(i-1)!=0)) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n7 = h->hpt_tmp->GetBinContent(i);
        int n7m1 = h->hpt_tmp->GetBinContent(i-1);
        h->hpt_ak8ak8m1->Fill(pt, min(n7,2), min(n7m1,2), _w);
      }
      if (h->hpt_tmp_ak4pf->GetBinContent(i)!=0 ||
          (i>1 && h->hpt_tmp_ak4pf->GetBinContent(i-1)!=0)) {
        double pt = h->hpt_tmp_ak4pf->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i);
        int n5m1 = h->hpt_tmp_ak4pf->GetBinContent(i-1);
        h->hpt_ak4ak4m1->Fill(pt, min(n5,2), min(n5m1,2), _w);
      }
      // Minus two
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          (i>2 && h->hpt_tmp_ak4pf->GetBinContent(i-2)!=0)) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i-2);
        int n7 = h->hpt_tmp->GetBinContent(i);
        h->hpt_ak4ak8m2->Fill(pt, min(n7,2), min(n5,2), _w);
      }
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          (i>2 && h->hpt_tmp->GetBinContent(i-2)!=0)) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n7 = h->hpt_tmp->GetBinContent(i);
        int n7m2 = h->hpt_tmp->GetBinContent(i-2);
        h->hpt_ak8ak8m2->Fill(pt, min(n7,2), min(n7m2,2), _w);
      }
      if (h->hpt_tmp_ak4pf->GetBinContent(i)!=0 ||
          (i>2 && h->hpt_tmp_ak4pf->GetBinContent(i-2)!=0)) {
        double pt = h->hpt_tmp_ak4pf->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i);
        int n5m2 = h->hpt_tmp_ak4pf->GetBinContent(i-2);
        h->hpt_ak4ak4m2->Fill(pt, min(n5,2), min(n5m2,2), _w);
      }
      // Plus one
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          h->hpt_tmp_ak4pf->GetBinContent(i+1)!=0) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i+1);
        int n7 = h->hpt_tmp->GetBinContent(i);
        h->hpt_ak4ak8p1->Fill(pt, min(n7,2), min(n5,2), _w);
      }
      if (h->hpt_tmp->GetBinContent(i)!=0 ||
          h->hpt_tmp->GetBinContent(i+1)!=0) {
        double pt = h->hpt_tmp->GetBinCenter(i);
        int n7 = h->hpt_tmp->GetBinContent(i);
        int n7p1 = h->hpt_tmp->GetBinContent(i+1);
        h->hpt_ak8ak8p1->Fill(pt, min(n7,2), min(n7p1,2), _w);
      }
      if (h->hpt_tmp_ak4pf->GetBinContent(i)!=0 ||
          h->hpt_tmp_ak4pf->GetBinContent(i+1)!=0) {
        double pt = h->hpt_tmp_ak4pf->GetBinCenter(i);
        int n5 = h->hpt_tmp_ak4pf->GetBinContent(i);
        int n5p1 = h->hpt_tmp_ak4pf->GetBinContent(i+1);
        h->hpt_ak4ak4p1->Fill(pt, min(n5,2), min(n5p1,2), _w);
      }

    } // for i
  } // _jp_ak4ak8

  // Unbiased generator spectrum (for each trigger)
  if (_mc) {
    for (int i = 0; i != gen_njt; ++i) {
      double y = gen_jty[i];
      if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {
        h->hpt_g0tw->Fill(gen_jtpt[i], _w);
      }

      if (h->ak4ak8) {
        double pt = gen_jtpt[i];
        double phi = gen_jtphi[i];
        if (pt >= h->ptmin && pt < h->ptmax) {
          for (int j = 0; j != t4gen_njt; ++j) {
            double y5 = t4gen_jty[j];
            double phi5 = t4gen_jty[j];
            double dr = tools::oplus(y-y5, phi-phi5);
            if (dr<0.25) {
              assert(h->hak4ak8g0);
              h->hak4ak8g0->Fill(t4gen_jtpt[j]/gen_jtpt[i], _w);
            }
          } // for j
        }
      } // _jp_ak4ak8
    }
  }
  //
  if (h->ismc) {

    // unfolding studies (Mikael)
    for (int i = 0; i != gen_njt; ++i) {

      double ygen = fabs(gen_jty[i]);

      for (int j = 0; j != njt; ++j) {

        //double yreco = fabs(jty[j]);
        bool id = (_jetids[j] && evtid && _pass);

        if ((ygen >= h->ymin && ygen < h->ymax && gen_jtpt[i]>_jp_recopt) &&
            //(yreco >= h->ymin && yreco < h->ymax)
            (jtpt[j]>_jp_recopt && id)) {

          double dr = tools::oplus(delta_phi(gen_jtphi[i], jtphi[j]),
                                   fabs(gen_jteta[i] - jteta[j]));
          if (dr < 0.25) {
            h->mT->Fill(gen_jtpt[i], jtpt[j], _w);
            h->mTf->Fill(gen_jtpt[i], jtpt[j], _w);
            h->mTuw->Fill(gen_jtpt[i], jtpt[j]);
            h->mTfuw->Fill(gen_jtpt[i], jtpt[j]);
          }
        } // rapidity bin
      } // for j
    } // for i
    //
    for (int i = 0; i != gen_njt; ++i) {
      double ygen = fabs(gen_jty[i]);
      if (ygen >= h->ymin && ygen < h->ymax && gen_jtpt[i]>_jp_recopt) {
        h->mx->Fill(gen_jtpt[i], _w);
        h->mxf->Fill(gen_jtpt[i], _w);
        h->mxuw->Fill(gen_jtpt[i]);
        h->mxfuw->Fill(gen_jtpt[i]);
      }
    } // for i
    //
    for (int j = 0; j != njt; ++j) {
      double yreco = fabs(jty[j]);
      bool id  = (_jetids[j] && evtid && _pass && jtpt[j]>_jp_recopt);
      if (yreco >= h->ymin && yreco < h->ymax && id) {
        h->my->Fill(jtpt[j], _w);
        h->myf->Fill(jtpt[j], _w);
        h->myuw->Fill(jtpt[j]);
        h->myfuw->Fill(jtpt[j]);
      }
    } // for j

    for (int i = 0; i != gen_njt; ++i) {
      double y = gen_jty[i];
      if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {
        h->hpt_g0->Fill(gen_jtpt[i], _w);
        assert(h->hpt_g0_tmp);
        h->hpt_g0_tmp->Fill(gen_jtpt[i]);
      }
    } // for i

    if (h->ak4ak8) {
      for (int i = 0; i != t4gen_njt; ++i) {
        double y = t4gen_jty[i];
        if (fabs(y) >= h->ymin && fabs(y) < h->ymax) {
          assert(h->hpt_g0_ak4pf);
          h->hpt_g0_ak4pf->Fill(t4gen_jtpt[i], _w);
          assert(h->hpt_g0_tmp_ak4pf);
          h->hpt_g0_tmp_ak4pf->Fill(t4gen_jtpt[i]);
        }
      } // for i
    } // _ak4pf


    if (h->ak4ak8) {

      for (int i = 1; i != h->hpt_g0_tmp_ak4pf->GetNbinsX()+1; ++i) {

        // Jet number correlations
        if (h->hpt_g0_tmp->GetBinContent(i)!=0 ||
            h->hpt_g0_tmp_ak4pf->GetBinContent(i)!=0) {
          double pt = h->hpt_g0_tmp->GetBinCenter(i);
          int n5 = h->hpt_g0_tmp_ak4pf->GetBinContent(i);
          int n7 = h->hpt_g0_tmp->GetBinContent(i);
          h->hpt_g0_ak4ak8->Fill(pt, min(n7,2), min(n5,2), _w);
        }
        // Jet number correlations for nearby off-diagonals
        if (h->hpt_g0_tmp->GetBinContent(i)!=0 ||
          (i>1 && h->hpt_g0_tmp_ak4pf->GetBinContent(i-1)!=0)) {
          double pt = h->hpt_g0_tmp->GetBinCenter(i);
          int n5 = h->hpt_g0_tmp_ak4pf->GetBinContent(i-1);
          int n7 = h->hpt_g0_tmp->GetBinContent(i);
          h->hpt_g0_ak4ak8m1->Fill(pt, min(n7,2), min(n5,2), _w);
        }
        //
        if (h->hpt_g0_tmp->GetBinContent(i)!=0 ||
            (i>2 && h->hpt_g0_tmp_ak4pf->GetBinContent(i-2)!=0)) {
          double pt = h->hpt_g0_tmp->GetBinCenter(i);
          int n5 = h->hpt_g0_tmp_ak4pf->GetBinContent(i-2);
          int n7 = h->hpt_g0_tmp->GetBinContent(i);
          h->hpt_g0_ak4ak8m2->Fill(pt, min(n7,2), min(n5,2), _w);
        }
        //
        if (h->hpt_g0_tmp->GetBinContent(i)!=0 ||
            h->hpt_g0_tmp_ak4pf->GetBinContent(i+1)!=0) {
          double pt = h->hpt_g0_tmp->GetBinCenter(i);
          int n5 = h->hpt_g0_tmp_ak4pf->GetBinContent(i+1);
          int n7 = h->hpt_g0_tmp->GetBinContent(i);
          h->hpt_g0_ak4ak8p1->Fill(pt, min(n7,2), min(n5,2), _w);
        }
      } // for i
    } // _jp_ak4ak8
  } // gen spectrum

} // fillBasic


// Write and delete histograms
void fillHistos::writeBasics()
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeBasics():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (map<string, vector<basicHistos*> >::iterator it = _histos.begin();
       it != _histos.end(); ++it) {
    for (unsigned int i = 0; i != it->second.size(); ++i) {

      // Luminosity information
      basicHistos *h = it->second[i];
      for (int j = 0; j != h->hlumi->GetNbinsX()+1; ++j) {
        h->hlumi->SetBinContent(j, _dt ? h->lumsum : 1. );
        h->hlumi2->SetBinContent(j, _dt ? h->lumsum2 : 1. );
      }

      delete h;//it->second[i];
    } // for i
  } // for it

  cout << "\nOutput stored in " << _outfile->GetName() << endl;
  _outfile->Close();
  _outfile->Delete();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeBasic() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeBasic


// Initialize basic histograms for trigger and eta bins
 void fillHistos::initRunHistos(string name, double ymin, double ymax) {

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initRunHistos("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile :
              new TFile(Form("output-%s-1.root",_type.c_str()), "RECREATE"));
  assert(f && !f->IsZombie());
  //assert(f->mkdir(name.c_str()));
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
  //TDirectory *dir = gDirectory;
  TDirectory *dir = f->GetDirectory(name.c_str()); assert(dir);
  dir->cd();

  runHistos *h = new runHistos(dir, ymin, ymax);
  _runhistos[name] = h;

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initRunHistos() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initRunHistos


// Fill run histograms
void fillHistos::fillRunHistos(string name)
{
  runHistos *h = _runhistos[name];
  assert(h);

  // Luminosity information
  if (_dt && h->lums[run][lbn]==0) {

    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
    //if (lum==0) return; // no lumi no pass
    // let pass, lum==0 can be Poisson fluctuation for a valid LS

    h->lumsum += lum;
    h->lumsum2 += lum2;
    h->runlums[run] += lum;
    h->runlums2[run] += lum2;
    h->lums[run][lbn] = 1;

    for (unsigned int i = 0; i != h->trg.size(); ++i) {

      string const& t = h->trg[i];
      double prescale(0);
      if (_prescales[t].find(run)==_prescales[t].end()) {
        if (_trigs.find(t)!=_trigs.end()) {
          *ferr << "Prescale not found for trigger " << t
                << " run " << run << endl << flush;
          assert(false);
        }
      }
      else
        prescale = _prescales[t][run]; //assert(prescale);
      h->runlums_trg[t][run] += (prescale ? lum / prescale : 0.);
    } // for i

    // Initialize counters for a new run
    if (h->lums[run].size()==1) {
      for (unsigned int i = 0; i != h->trg.size(); ++i) {

        string const& t = h->trg[i];
        h->p_trg[t][run] = 0;
        h->t_trg[t][run] = 0;
        h->npv_trg[t][run] = 0;
        h->c_chf[t][run] = 0;
        h->c_nef[t][run] = 0;
        h->c_nhf[t][run] = 0;
        h->c_betastar[t][run] = 0;
        h->t_trgtp[t][run] = 0;
        h->c_chftp[t][run] = 0;
        h->c_neftp[t][run] = 0;
        h->c_nhftp[t][run] = 0;
        h->c_betastartp[t][run] = 0;

        for (unsigned int j = 0; j != h->trg.size(); ++j) {
          string const& t2 = h->trg[j];
          h->p_trgpair[t+t2][run] = 0;
        } // for j
      } // for i
    } // new run
  }

  //bool evtid = (met < 0.5 * metsumet || met < 100.);
  bool evtid = (met < 0.4 * metsumet || met < 45.); // QCD-11-004

  double dphi = (njt>=2 ? delta_phi(jtphi[0], jtphi[1]) : 0.);
  double pt3 = (njt>=3 ? jtpt[2] : 0.);

  for (int i = 0; i != njt; ++i) {

    double pt = jtpt[i];
    double y = jty[i];

    if (h->ymin <= fabs(y) && fabs(y) < h->ymax && _pass && _jetids[i]
        && evtid) {

      for (set<string>::const_iterator it = _trigs.begin(); it != _trigs.end(); ++it) {
        string const& t = *it;

        if (pt > 18.) ++h->p_trg[t][run];
        if (pt > h->pt[t]) {
          ++h->t_trg[t][run]; // unweighted events
          h->tw_trg[t][run] += _prescales[t][run]; // prescale weighted events
          h->npv_trg[t][run] += npv;
          h->npvgood_trg[t][run] += npvgood;
          h->c_chf[t][run] += jtchf[i];
          h->c_nef[t][run] += jtnef[i];
          h->c_nhf[t][run] += jtnhf[i];
          h->c_betastar[t][run] += jtbetastar[i];
        }

        int iref = (i==0 ? 1 : 0);
        if (i<2 && dphi > 2.7 && pt3 < jtpt[iref] && fabs(jty[iref]) < 1.3 &&
            jtpt[iref] > h->pt[t] && _jetids[iref]) {
          ++h->t_trgtp[t][run];
          h->c_chftp[t][run] += jtchf[i];
          h->c_neftp[t][run] += jtnef[i];
          h->c_nhftp[t][run] += jtnhf[i];
          h->c_betastartp[t][run] += jtbetastar[i];
        }

        for (set<string>::const_iterator jt = _trigs.begin();
             jt != _trigs.end(); ++jt) {

          string const& t2 = *jt;
          if (t!=t2) ++h->p_trgpair[t+t2][run];
        } // for jt
      } // for it
    }
  } // for i

} // fillRunHistos


// Write and delete histograms
void fillHistos::writeRunHistos()
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeRunHistos():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (map<string, runHistos*>::iterator it = _runhistos.begin();
       it != _runhistos.end(); ++it) {

    runHistos *h = it->second;
    delete h;
  } // for it

  cout << "\nOutput (runHistos) stored in " << _outfile->GetName() << endl;

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeRunHistos() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeRunHistos


void fillHistos::fillJetID(vector<bool> &id)
{
  assert(int(id.size())==njt);

  for (int i = 0; i != njt; ++i) {

    id[i] = ((fabs(jty[i])<2.5 ? jtidtight[i] : jtidloose[i]));
             //&& (_mc || jtbeta[i]!=0));
             //&& (_mc || (1-jtbetastar[i])>0.5));

    if (_jp_doECALveto) {
      assert(ecalveto);
      int ibin = ecalveto->FindBin(jteta[i],jtphi[i]);
      id[i] = (id[i] && ecalveto->GetBinContent(ibin)==0);
    }
  }

} // fillJetID


void fillHistos::fillJetID5(vector<bool> &id)
{
  assert(int(id.size())==t4_njt);

  for (int i = 0; i != t4_njt; ++i) {

    id[i] = ((fabs(t4_jty[i])<2.5 ? t4_jtidtight[i] : t4_jtidloose[i]));

    if (_jp_doECALveto) {
      assert(ecalveto);
      int ibin = ecalveto->FindBin(t4_jteta[i],t4_jtphi[i]);
      id[i] = (id[i] && ecalveto->GetBinContent(ibin)==0);
    }
  }

} // fillJetID5


// Load good run and LS information
void fillHistos::loadJSON(const char* filename)
{
  cout << "Processing loadJSON(\"" << filename << "\"..." << endl;
  ifstream file(filename, ios::in);
  assert(file.is_open());
  char c;
  string s, s2;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  assert(c=='{');
  while (file >> s && sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (_debug)
      cout << "\"" << rn << "\": " << flush;

    while (file.get(c) && c==' ') {};
    if (_debug)
      cout << c << flush; assert(c=='[');
    ++nrun;

    bool endrun = false;
    while (!endrun && file >> s >> s2 &&
           sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (_debug)
        cout << "["<<ls1<<","<<ls2<<"]"<<s1 << flush;

      for (int ls = ls1; ls != ls2+1; ++ls) {
        //assert(_json[rn].find(ls)==_json[rn].end()); // ok if 2 JSON files
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (!endrun && s2!=",") {
        if (_debug)
          cout<<"s1: "<<s2<<endl<<flush; assert(s2==",");
      }
    } // while ls
    if (_debug)
      cout << endl;

    if (s2=="]}") continue;
    else if (s2!="],") {
      if (_debug)
        cout<<"s2: "<<s2<<endl<<flush; assert(s2=="],");
    }
  } // while run
  if (s2!="]}") { cout<<"s3: "<<s2<<endl<<flush; assert(s2=="]}"); }

  cout << "Called loadJSON(\"" << filename << "\"):" << endl;
  cout << "Loaded " << nrun << " good runs and " << nls
       << " good lumi sections" << endl;

} // loadJSON


// Load luminosity information
void fillHistos::loadLumi(const char* filename)
{
  cout << "Processing loadLumi(\"" << filename << "\")..." << endl;

  // Check lumi against the list of good runs
  const int a_goodruns[] = {};
  const int ngoodruns = sizeof(a_goodruns)/sizeof(a_goodruns[0]);
  set<int> goodruns;
  for (int i = 0; i != ngoodruns; ++i) {
    goodruns.insert(a_goodruns[i]);
  }

  set<pair<int, int> > nolums;

  if (true) {
    for (set<int>::const_iterator it = goodruns.begin(); it != goodruns.end(); ++it) {
      cout << *it << ", ";
    }
    cout << endl;
  }

  ifstream f(filename, ios::in);
  assert(f.is_open());
  float secLS = 2.3310e+01;
  string s;
  int rn, fill, ls, ifoo;
  float del, rec, avgpu, ffoo;
  char sfoo[512];
  //assert(getline(f, s, '\r'));
  assert(getline(f, s, '\n'));
  //assert(f >> s);
  cout << endl << "string: " << s << " !" << endl << flush;
  //bool v1 = (s==string("run,ls,delivered,recorded"));
  //bool v2 = (s==string("Run,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)"));
  //bool v3 = (s==string("Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)"));
  //assert(v1 || v2 || v3);
  assert (s=="#Data tag : v1 , Norm tag: None");
  //assert (s=="#Data tag : online , Norm tag: None");
  //
  //assert(getline(f, s, '\r'));
  assert(getline(f, s, '\n'));
  //assert(f>>s);
  cout << endl << "string: " << s << " !" << endl << flush;
  assert (s=="#run:fill,ls,time,beamstatus,E(GeV),delivered(/ub),recorded(/ub),avgpu,source");

  int nls(0);
  double lumsum(0);
  double lumsum_good(0);
  double lumsum_json(0);
  bool skip(false);
  /*
  while ((v1 && f >> s &&
    sscanf(s.c_str(),"%d,%d,%f,%f",&rn,&ls,&del,&rec)==4) ||
   (v2 && getline(f, s, '\r') &&
    (sscanf(s.c_str(),"%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,"
      "%f,%f,%f", &rn,&ls,&ifoo, &ifoo,&ifoo,&ifoo,
      &ifoo,&ifoo,&ifoo, &ffoo,&del,&rec)==12 ||
     (skip=true))) ||
   (v3 && getline(f, s, '\r') &&
    (sscanf(s.c_str(),"%d:%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,"
      "%f,%f,%f", &rn,&fill,&ls,&ifoo, &ifoo,&ifoo,&ifoo,
      &ifoo,&ifoo,&ifoo, &ffoo,&del,&rec)==13 ||
     (skip=true)))) {
  */
  while (//getline(f, s, '\r') &&
         getline(f, s, '\n') &&
         //f >> s &&
         (sscanf(s.c_str(),"%d:%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,"
          "%f,%f,%f,%f,%s", &rn,&fill,&ls,&ifoo, &ifoo,&ifoo,&ifoo,
           &ifoo,&ifoo,&ifoo, &ffoo,&del,&rec,&avgpu,sfoo)==15 ||
         (skip=true))) {

    if (_debug) {
      if (skip) cout << "Skipping line:\n" << s << endl;
      cout << "Run " << run << " ls " << ls
           << " lumi " << rec*1e-6 << "/pb" << endl;
    }

    // LS is not STABLE BEAMS but something else:
    // ADJUST, BEAM DUMP, FLAT TOP, INJECTION PHYSICS BEAM, N/A, RAMP DOWN,
    // SETUP, SQUEEZE
    if (skip) {
      skip = false;
      continue;
    }

    assert(_lums[rn][ls]==0);
    assert(_avgpu[rn][ls]==0);
    // Try to get this in units of pb-1
    // apparently it is given in s^-1 cm^-2
    //double lum = lvtx * secLS * 1e-36 ;
    //if (lum==0) lum = lhf * secLS * 1e-36 ;
    // lumiCalc.py returns lumi in units of mub-1 (=>nb-1=>pb-1)
    double lum = rec*1e-6;
    //double lum2 = lhf * secLS * 1e-36 ;
    //if (lum2==0) lum2 = lvtx * secLS * 1e-36 ;
    double lum2 = del*1e-6;
    //assert(lum!=0);
    if (lum==0 && goodruns.find(rn)!=goodruns.end() &&
        (!_jp_dojson || _json[rn][ls]==1)) {
      //cerr << "Warning: run " << rn << " LS " << ls << " is zero!" << endl;
      nolums.insert(pair<int, int>(rn,ls));
    }

    _avgpu[rn][ls] = avgpu * 69000. / 78400.; // brilcalc --minBiasXsec patch
    _lums[rn][ls] = lum;
    _lums2[rn][ls] = lum2;
    lumsum += lum;
    if (goodruns.find(rn)!=goodruns.end()) // Apr 17
      lumsum_good += lum;
    if ((!_jp_dojson || _json[rn][ls]))
      lumsum_json += lum;
    ++nls;
    assert(nls<100000000);
  }

  cout << "Called loadLumi(\"" << filename << "\"):" << endl;
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
    for (set<pair<int, int> >::const_iterator it = nolums.begin();
         it != nolums.end(); ++it) {

      cout << " ["<<it->first<<","<<it->second;
      set<pair<int, int> >::const_iterator jt = it; ++jt;
      if (jt->first!=it->first || jt->second!=it->second+1)
        cout << "]";
      else {
        for (int i = 0; jt!=nolums.end() && jt->first==it->first
               && jt->second==it->second+i+1; ++i, ++jt) {};
        it = --jt;
        cout << "-" << it->second << "]";
      }
    } // for it
    cout << endl;
  } // nolums

} // loadLumi


void fillHistos::loadPUProfiles(const char *datafile, const char *mcfile)
{
  cout << "Processing loadPUProfiles(\"" << datafile << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  // Load pile-up files and hists from them
  TFile *fpudist = new TFile(datafile, "READ");
  assert(fpudist && !fpudist->IsZombie());
  TFile *fpumc = new TFile(mcfile,"READ");
  assert(fpumc && !fpumc->IsZombie());

  //pumc = (TH1D*)fpumc->Get("hpu370"); assert(pumc);
  pumc = (TH1F*)fpumc->Get("pileupmc"); assert(pumc);
  pudt = (TH1F*)fpudist->Get("pileupdt"); assert(pudt);

  // Normalize
  pumc->Scale(1./pumc->Integral());
  pudt->Scale(1./pudt->Integral());

  // For data, load each trigger separately
  /*
  for (unsigned int itrg = 0 ; itrg != _triggers.size(); ++itrg) {

    const char *t = _triggers[itrg].c_str();
    pudist[t] = (TH1D*)fpudist->Get(Form("pileup_%s",t)); assert(pudist[t]);
    pudist[t]->Scale(1./pudist[t]->Integral());
  }
  */

  curdir->cd();
} // loadPUProfiles


void fillHistos::loadPrescales(const char *prescalefile)
{
  cout << "Processing loadPrescales(\"" << prescalefile << "\")..." << endl;
  fstream fin(prescalefile);

  const int ns = 1024;
  char s[ns];

  fin.getline(s, ns);
  stringstream ss;
  ss << s;

  string srun, sls, strg;
  vector<string> trgs;
  ss >> srun; assert(srun=="RUN");

  while (ss >> strg) trgs.push_back(strg);

  int run, ls, pre;
  while(fin.getline(s, ns)) {

    stringstream ss;
    ss << s;
    ss >> run >> ls;
    if (_debug) cout << " run " << run << " ls " << ls << ": ";
    //ss >> run;
    //cout << " run " << run << ": ";

    int itrg(0);
    for (unsigned int itrg = 0; ss >> pre; ++itrg) {
      //cout << "trg" << trgs[itrg] << " run " << run << " ls " << ls;
      assert(itrg!=trgs.size());
      //_premap[trgs[itrg]][run] = pre;
      _premap[trgs[itrg]][run][ls] = pre;
      if (_debug) cout << pre << "/" << trgs[itrg] << " ";
    }
    if (_debug) cout << endl;
  }
} // loadPrescales


void fillHistos::loadECALveto(const char *file)
{
  cout << "Processing loadECALveto(\"" << file << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  TFile *fe = new TFile(file, "READ");
  assert(fe && !fe->IsZombie());

  ecalveto = (TH2F*)fe->Get("ecalveto"); assert(ecalveto);

  curdir->cd();
} // loadECALveto


// Update the available trigger types for each new tree
void fillHistos::getTriggers()
{
  TH1F *triggers = (TH1F*)fChain->GetCurrentFile()->Get("ak4/TriggerNames");
  TAxis *xax = triggers->GetXaxis();

  std::regex pfjet("HLT_PFJet([0-9]*)_v[0-9]*");

  _availTrigs.clear();
  for (int i = xax->GetFirst(); i <= xax->GetLast(); ++i) {
    string trgName = xax->GetBinLabel(i);
    if (trgName.compare("")==0) // Ignore empty places on x-axis
      continue;
    string trigger = ""; // HLT_PFJet are given non-empty trigger names
    if (std::regex_match(trgName,pfjet))
      trigger=std::regex_replace(trgName, pfjet, "jt$1", std::regex_constants::format_no_copy);
    _availTrigs.push_back(trigger);
  }
} // getTriggers
