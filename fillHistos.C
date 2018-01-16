// Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: June 2, 2015
// Updated: Aug 31, 2016
// Further updates: see git log

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

  _entries = ntot;//(_jp_ismc ? fChain->GetEntries() : ntot);
  _nbadevts_dup = _nbadevts_run = _nbadevts_ls = _nbadevts_lum = 0;
  _nbadevts_veto = _nbadevts_stream = 0;
  _bscounter_bad = _bscounter_good = _halocounter_bad = _halocounter_good = 0;
  _ecalcounter_good = _ecalcounter_bad = 0;
  _rhocounter_good = _rhocounter_bad = 0;
  _trgcounter = _evtcounter = _totcounter = 0;
  _ecalveto = 0;

  // Set cross section weights for pThat bins
  if (_jp_pthatbins)
    _pthatweight = 0;

  // Map of mu per (run,lumi)
  TFile *fmu = new TFile("pileup/MUperLSvsRUN_MB.root","READ");
  assert(fmu && !fmu->IsZombie());
  TH2F *h2mu = (TH2F*)fmu->Get("hLSvsRUNxMU");
  assert(h2mu && !h2mu->IsZombie());
  //TH2F *h2mu = (TH2F*)fmu->Get("hLSvsRuNxMU_cleaned"); assert(h2mu);

  if (_jp_doqgl) { // Qgl: load quark/gluon probability histos (Ozlem)
    TFile *finmc = new TFile(_jp_qglfile,"READ");
    assert(finmc && !finmc->IsZombie());

    double veta[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7};
    const int neta = sizeof(veta)/sizeof(veta[0])-1;
    // same or vpt, vtrigpt, npt
    const int npt = _jp_ntrigs;
    const double *vtrigpt = &_jp_trigthr[0];
    double vpt[npt+1]; vpt[0] = 0;
    for (int trigidx = 0; trigidx != npt; ++trigidx) vpt[trigidx+1] = _jp_trigranges[trigidx][1];
    const int nqgl = 101; double vqgl[nqgl+1];
    for (int qglidx = 0; qglidx != nqgl+1; ++qglidx) vqgl[qglidx] = 0. + 0.01*qglidx;

    h3probg = new TH3D("h3probg","Gluon prob.;#eta_{jet};p_{T,jet};QGL",
      //6,0,3.0, 10,0,1000, 101,0,1.01);
      neta,veta, npt,vpt, nqgl,vqgl); //101,0,1.01);

    // Loop over pt and eta bins in the gql histos
    for (int ieta = 1; ieta != h3probg->GetNbinsX()+1; ++ieta) {
      for (int ipt = 1; ipt != h3probg->GetNbinsY()+1; ++ipt) {
        // Gluons in the given eta/pt slice
        string sg = Form("Standard/Eta_%1.1f-%1.1f/jt%1.0f/hqgl_g",
                         veta[ieta-1], veta[ieta], vtrigpt[ipt-1]);
        cout << sg << endl << flush;
        TH1D *hqgl_g = (TH1D*)finmc->Get(sg.c_str());
        assert(hqgl_g);
        hqgl_g->Scale(1./hqgl_g->Integral());

        // Quarks in the given eta/pt slice
        string sq = Form("Standard/Eta_%1.1f-%1.1f/jt%1.0f/hqgl_q",
                         veta[ieta-1], veta[ieta], vtrigpt[ipt-1]);
        TH1D *hqgl_q = (TH1D*)finmc->Get(sq.c_str());
        assert(hqgl_q);
        hqgl_q->Scale(1./hqgl_q->Integral());

        // Loop over qgl indices
        for (int iqgl = 1; iqgl != h3probg->GetNbinsZ()+1; ++iqgl) {
          // Fetch gluon probability
          double probg = hqgl_g->GetBinContent(iqgl) /
                        (hqgl_g->GetBinContent(iqgl) + hqgl_q->GetBinContent(iqgl));
          h3probg->SetBinContent(ieta, ipt, iqgl, probg);
        } // for iqgl
      } // for ipt
    } // for ieta
  }

  if (_jp_quick) {

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
    //fChain->SetBranchStatus(Form("PFJets%s_.phf_",_jp_chs),1); // jtnef
    fChain->SetBranchStatus(Form("PFJets%s_.nemf_",_jp_chs),1); // jtnef
    fChain->SetBranchStatus(Form("PFJets%s_.nhf_",_jp_chs),1); // jtnhf
    //fChain->SetBranchStatus(Form("PFJets%s_.elf_",_jp_chs),1); // jtcef !!
    fChain->SetBranchStatus(Form("PFJets%s_.cemf_",_jp_chs),1); // jtcef !!
    fChain->SetBranchStatus(Form("PFJets%s_.muf_",_jp_chs),1); // jtmuf !!
    fChain->SetBranchStatus(Form("PFJets%s_.ncand_",_jp_chs),1); // jtn
    fChain->SetBranchStatus(Form("PFJets%s_.beta_",_jp_chs),1); // jtbeta
    fChain->SetBranchStatus(Form("PFJets%s_.betaStar_",_jp_chs),1); // jtbetastar
    fChain->SetBranchStatus(Form("PFJets%s_.chm_",_jp_chs),1); // jtnch
    fChain->SetBranchStatus(Form("PFJets%s_.phm_",_jp_chs),1); // jtnne
    fChain->SetBranchStatus(Form("PFJets%s_.nhm_",_jp_chs),1); // jtnnh
    fChain->SetBranchStatus(Form("PFJets%s_.elm_",_jp_chs),1); // jtnce !!
    fChain->SetBranchStatus(Form("PFJets%s_.mum_",_jp_chs),1); // jtnmu !!
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

  // Set pointers to branches
  jtp4x = &PFJetsCHS__P4__fCoordinates_fX[0];
  jtp4y = &PFJetsCHS__P4__fCoordinates_fY[0];
  jtp4z = &PFJetsCHS__P4__fCoordinates_fZ[0];
  jtp4t = &PFJetsCHS__P4__fCoordinates_fT[0];
  jta = &PFJetsCHS__area_[0];
  jtjes = &PFJetsCHS__cor_[0];
  jtbeta = &PFJetsCHS__beta_[0];
  jtbetastar = &PFJetsCHS__betaStar_[0];
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

  //
  jtn = &PFJetsCHS__ncand_[0];
  jtnch = &PFJetsCHS__chm_[0];
  jtnnh = &PFJetsCHS__nhm_[0];
  jtnne = &PFJetsCHS__phm_[0];
  jtnce = &PFJetsCHS__elm_[0];
  jtnmu = &PFJetsCHS__mum_[0];
  jtchf = &PFJetsCHS__chf_[0];
  jtnhf = &PFJetsCHS__nhf_[0];
  //jtnef = &PFJetsCHS__phf_[0];
  jtnef = &PFJetsCHS__nemf_[0];
  //jtcef = &PFJetsCHS__elf_[0];
  jtcef = &PFJetsCHS__cemf_[0];
  jtmuf = &PFJetsCHS__muf_[0];
  //
  gen_jtp4x = &GenJets__fCoordinates_fX[0];
  gen_jtp4y = &GenJets__fCoordinates_fY[0];
  gen_jtp4z = &GenJets__fCoordinates_fZ[0];
  gen_jtp4t = &GenJets__fCoordinates_fT[0];

  const char *a = _jp_algo;
  cout << "\nCONFIGURATION DUMP:" << endl;
  cout << "-------------------" << endl;
  cout << Form("Running over %sPF",a) << endl;
  if(_jp_ismc) cout << "Running over MC" << endl;
  if(_jp_isdt) cout << "Running over data" << endl;
  cout << (_jp_useIOV ? "Applying" : "Not applying")
       << " time-dependent JEC (IOV)" << endl;
  cout << (_jp_doECALveto ? "Vetoing" : "Not vetoing")
       << " jets in bad ECAL towers" << endl;
  cout << endl;

  if (_jp_isdt) {
    cout << (_jp_dojson ? "Applying" : "Not applying")
         << " additional JSON selection" << endl;
    cout << (_jp_doRunHistos ? "Storing" : "Not storing")
         << " additional run-level histograms" << endl;
    cout << (_jp_doBasicHistos ? "Storing" : "Not storing")
         << " basic set of histograms" << endl;
    cout << (_jp_doEtaHistos ? "Storing" : "Not storing")
         << " histograms with a full eta-range" << endl;
  }
  if (_jp_ismc) {
    cout << (_jp_reweighPU ? "Reweighing" : "Not reweighing")
         << " pileup profile in MC to data" << endl;
    cout << (_jp_pthatbins ? "Processing pThat binned samples"
             : "Processing \"flat\" samples") << endl;
  }
  cout << endl;

  _JEC = 0;
  _L1RC = 0;
  _jecUnc = 0;

  // Time dependent JEC (only for dt)
  if (_jp_isdt && _jp_useIOV) {
    for (unsigned iovidx=0; iovidx<_jp_nIOV; ++iovidx)
      _iov.add(_jp_IOVnames[iovidx],_jp_jecgt,_jp_jecvers,_jp_IOVranges[iovidx][0],_jp_IOVranges[iovidx][1]);
  } else {
    // At least a singular recalculation of JEC is always performed
    const char *s;
    const char *p = "CondFormats/JetMETObjects/data/";
    string jecgt = _jp_jecgt + _jp_jecvers + "_" + (_jp_isdt ? "DATA" : "MC") + "_";
    const char *t = jecgt.c_str();

    cout << "Loading "<<a<<"PF JEC" << endl;
    vector<JetCorrectorParameters> vpar;

    s = Form("%s%sL1FastJet_%s.txt",p,t,a);
    cout<<s<<endl<<flush;
    vpar.push_back(JetCorrectorParameters(s));

    s = Form("%s%sL2Relative_%s.txt",p,t,a);
    cout<<s<<endl<<flush;
    vpar.push_back(JetCorrectorParameters(s));

    s = Form("%s%sL3Absolute_%s.txt",p,t,a);
    cout<<s<<endl<<flush;
    vpar.push_back(JetCorrectorParameters(s));

    if (_jp_isdt) {
      if (!_jp_skipl2l3res) {
        s = Form("%s%sL2L3Residual_%s.txt",p,t,a);
        cout<<s<<endl<<flush;
        vpar.push_back(JetCorrectorParameters(s));
      }

      s = Form("%s%sUncertainty_%s.txt",p,t,a);
      cout<<"**"<<s<<endl<<flush;
      _jecUnc = new JetCorrectionUncertainty(s);
    }
    _JEC = new FactorizedJetCorrector(vpar);

    // For type-I and type-II MET
    vector<JetCorrectorParameters> vrc;
    s = Form("%s%sL1RC_%s.txt",p,t,a);
    cout<<s<<endl<<flush;
    vrc.push_back(JetCorrectorParameters(s));
    _L1RC = new FactorizedJetCorrector(vrc);

    assert(_JEC);
    assert(_L1RC);
  } // JEC redone


  // Set list of triggers
  for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {
    _triggers.push_back(_jp_triggers[itrg]);
  }

  // Load latest JSON selection
  if (_jp_isdt && _jp_dojson) loadJSON(_jp_json);

  // Load PU profiles for MC reweighing
  if (_jp_ismc && _jp_reweighPU)
    loadPUProfiles(_jp_pudata, _jp_pumc);

  // Load prescale information to patch 76X
  if (_jp_prescalefile!="")
    loadPrescales(_jp_prescalefile);

  // load ECAL veto file for cleaning data
  if (_jp_doECALveto) loadECALveto(_jp_ecalveto);

  // REMOVED: "manual veto list"

  // load luminosity tables (prescales now stored in event)
  if (_jp_isdt && _jp_dolumi) loadLumi(_jp_lumifile);

  if (_jp_ismc) cout << Form("Running on MC produced with %1.3g nb-1 (%ld evts)",
                        1000. * _entries / _jp_xsecMinBias,
                        (long int)_entries) << endl;
  if (_jp_isdt) cout << Form("Running on %ld events of data",
                        (long int)_entries) << endl;

  // Initialize histograms for different epochs and DQM selections
  if (_jp_doBasicHistos) {
    initBasics("Standard");
  }

  if (_jp_doEtaHistos) {
    initEtas("FullEta_Reco");
    if (_jp_ismc) {
      initEtas("FullEta_Gen");
      if (_jp_doEtaHistosMcResponse) {
        initMcHistos("FullEta_RecoPerGen_vReco");
        initMcHistos("FullEta_RecoPerGen_vGen");
      }
    }
  }

  if (_jp_isdt && _jp_doRunHistos) {
    initRunHistos("Runs",0.,3.);
    initRunHistos("RunsBarrel",0.,1.);
    initRunHistos("RunsTransition",1.,2.);
    initRunHistos("RunsEndcap",2.,3.);
  }

  if (_jp_etaphiexcl) {
    fetaphiexcl = new TFile("hotjets.root","READ");
    assert(fetaphiexcl && !fetaphiexcl->IsZombie() && "file hotjets.root missing");
    h2etaphiexcl = (TH2D*)fetaphiexcl->Get(_jp_etaphitype);
    assert(h2etaphiexcl && "erroneous eta-phi exclusion type");
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

     if (_js_skim != 0 && jentry%_js_skim!=0) continue;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    _entry = jentry;

    if (jentry%50000==0) cout << "." << flush;

    if (jentry==10000 || jentry==100000 || jentry==1000000 || jentry==5000000 || jentry==20000000 || jentry==40000000 || jentry==80000000){
      cout << endl << Form("Processed %ld events (%1.1f%%) in %1.0f sec. ETA:",
                      (long int)jentry, 100.*jentry/ntot, stop.RealTime()) << endl;
      TDatime now; now.Set(now.Convert()+stop.RealTime()*ntot/jentry);
      now.Print();
      stop.Continue();
    }

    // Set auxiliary event variables (jets, triggers later)
    assert(_jp_isdt || (_jp_pthatbins && _pthatweight));
    pthat = EvtHdr__mPthat;
    weight = EvtHdr__mWeight;
    if (_jp_ismc && _jp_pthatbins)
      weight *= _pthatweight;
    // REMOVED: "TEMP PATCH"
    run = EvtHdr__mRun;
    evt = EvtHdr__mEvent;
    lbn = EvtHdr__mLumi;

    trpu = EvtHdr__mTrPu;
    itpu = EvtHdr__mINTPU;
    ootpulate = EvtHdr__mOOTPULate;
    ootpuearly = EvtHdr__mOOTPUEarly;

    if (_jp_isdt) {
      trpu = _avgpu[run][lbn];
      if (trpu==0) {
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

    njt = PFJetsCHS__;       //assert(njt < kMaxPFJetsCHS_);
    gen_njt = GenJets__;
    gen_njt = GenJets__;

    //assert(njt<_njt);
    if (!(njt < _njt)) {
      *ferr << "Array overflow: njt = " << njt
           << " > njtmax=" << _njt << endl;
      cout << "Array overflow: njt = "<< njt
           << " > njtmax=" << _njt << endl;
      cout << flush;
      assert(njt<_njt);
    }

     //assert(_jp_isdt || gen_njt<_njt);
    if (_jp_ismc && !(gen_njt < _njt)) {
      *ferr << "Array overflow: gen_njt = " << gen_njt
           << " > njtmax=" << _njt << endl;
      cout << "Array overflow: gen_njt = "<< njt
           << " > njtmax=" << _njt << endl;
      cout << flush;
      assert(gen_njt<_njt);
    }

    if (_jt_debug) {
      cout << endl << flush;
      Show(jentry);
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
    if (_jp_isdt && _jp_checkduplicates) {
      set<int>& events = _duplicates[run][lbn];
      if (events.find(evt)!=events.end()) {
        ++_nbadevts_dup;
        continue;
      }
      events.insert(evt);
    }

    ++cnt["01all"];

    // Check if good run/LS, including JSON selection
    if (_jp_isdt && _jp_dojson) {

      // Does the run/LS pass the latest JSON selection?
      if (_json[run][lbn]==0) {
        _badjson.insert(pair<int, int>(run, lbn));
        ++_nbadevts_json;
        continue;
      }
    } // _jp_isdt && _jp_dojson

    if (_jp_isdt && _jp_dolumi) {
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
    } // _jp_isdt && _jp_dolumi

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
    if (_pass && (tools::oplus(pvx-bsx, pvy-bsy)>0.15 || pvndof<=4 || fabs(pvz) >= 24.))
    {
      ++_bscounter_bad;
      _pass = false;
    }
    if (_pass) ++_bscounter_good;
    if (_pass) ++cnt["04bsc"];

    // Event cuts against beam backgrounds
    if (_pass && _ecalveto &&
         ( (njt>=1 && _ecalveto->GetBinContent(_ecalveto->FindBin(jteta[0],jtphi[0]))!=0)
        || (njt>=2 && _ecalveto->GetBinContent(_ecalveto->FindBin(jteta[1],jtphi[1]))!=0) ))
    {
      ++_ecalcounter_bad;
      _pass = false;
    } // ecal veto
    if (_pass) ++_ecalcounter_good;
    if (_pass) ++cnt["05ecal"];

    // Check rho
    if (_pass && rho>40.) {
      ++_rhocounter_bad;
      _pass = false;
      if (_jt_debug)
        cout << Form("\nrun:ev:ls %d:%d:%d : rho=%1.1f njt=%d npv=%d"
                     " jtpt0=%1.1f sumet=%1.1f met=%1.1f\n",
                     run, lbn, evt, rho, njt, npv,
                     (njt>0 ? jtpt[0] :0.), metsumet, met) << flush;
    }
    if (_pass) ++_rhocounter_good;
    if (_pass) ++cnt["06rho"];


    // Reset prescales (dynamic can change within run)
    for (auto scaleit = _prescales.begin(); scaleit != _prescales.end(); ++scaleit)
      scaleit->second[run] = 0;

    // Fill trigger information
    _trigs.clear();

    // Simulate other triggers for MC, if so wished
    // (this is slow, though)
    if (_jp_ismc) {
      // Always insert the generic mc trigger
      _trigs.insert("mc");
      if (_jp_domctrigsim && njt>0) {
        // Only add the greatest trigger present
        for (int itrg = _jp_ntrigs; itrg > 0; --itrg) {
          if (jtpt[0]>_jp_trigranges[itrg-1][0]) {
            _trigs.insert(_jp_triggers[itrg-1]);
            break; // Don't add lesser triggers
          }
        }
      }
    } // _jp_domctrigsim

    if (_jp_isdt) {
      // For data, check trigger bits
      if (_jt_debug) cout << "TriggerDecision_.size()=="<<TriggerDecision_.size()<<endl<<flush;
      if (_jt_debug) cout << "_availTrigs.size()=="<<_availTrigs.size()<<endl<<flush;
      assert(TriggerDecision_.size() == _availTrigs.size());

      for (unsigned int itrg = 0; itrg != TriggerDecision_.size(); ++itrg) {

        string strg = _availTrigs[itrg];
        bool pass = TriggerDecision_[itrg]==1 && strg.length()!=0; // -1, 0, 1

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
          if (_jt_debug) {
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
            // Make sure all info is good! This is crucial if there is something odd with the tuples
            *ferr << "Missing prescale for " << strg
                  << " in run " << run << endl << flush;
          }
        }
      } // for itrg
    }

    ++_totcounter;
    if (_pass) ++_evtcounter;
    if (_trigs.size()!=0 && _pass) ++_trgcounter;
    if (_trigs.size()!=0 && _pass && _jp_isdt) ++cnt["07trg"];

    // Retrieve event weight
    _w0 = (_jp_ismc ? weight : 1);
    assert(_w0);
    _w = _w0;

    // Calculate trigger PU weight
    for (unsigned int itrg = 0; itrg != _triggers.size(); ++itrg) {
      const char *trg_name = _triggers[itrg].c_str();
      _wt[trg_name] = 1.;

      // Reweigh in-time pile-up
      if (_jp_ismc && _jp_reweighPU) {
        int k = pudist[trg_name]->FindBin(trpu);
        double w1 = pudist[trg_name]->GetBinContent(k);
        double w2 = pumc->GetBinContent(k);
        Double_t wtrue = (w1==0 || w2==0 ? 1. : w1 / w2);
        _wt[trg_name] *= wtrue;

        // check for non-zero PU weight
        if (_pass)
          _pass = (pudist[trg_name]->GetBinContent(pudist[trg_name]->FindBin(trpu))!=0);
      }
    } // for itrg
    _wt["mc"] = _wt[_jp_mctrig];
    if (_trigs.size()!=0 && _pass && _jp_ismc) ++cnt["07puw"];

    // TODO: implement reweighing for k-factor (NLO*NP/LOMC)

    // load correct IOV for JEC
    if (_jp_isdt && _jp_useIOV) {
      assert(_iov.setCorr(run,&_JEC,&_L1RC,&_jecUnc));
      assert(_JEC);
      assert(_L1RC);
      assert(_jecUnc);
    }

    // Calculate pT, eta, phi, y, E and uncorrected pT
    // oversmear jets and MET in MC
    double mex = met * cos(metphi);
    double mey = met * sin(metphi);
    for (int jetidx = 0; jetidx != njt; ++jetidx) {
      p4.SetPxPyPzE(jtp4x[jetidx],jtp4y[jetidx],jtp4z[jetidx],jtp4t[jetidx]);
      // Divide by the original JES
      if (_jp_undojes)
        p4 *= 1/jtjes[jetidx];

      jtptu[jetidx] = p4.Pt();
      jteu[jetidx] = p4.E();

      // Recalculate JEC
      _JEC->setRho(rho);
      _JEC->setNPV(npvgood);
      _JEC->setJetA(jta[jetidx]);
      _JEC->setJetPt(jtptu[jetidx]);
      _JEC->setJetE(jteu[jetidx]);
      _JEC->setJetEta(p4.Eta());
      jtjesnew[jetidx] = _JEC->getCorrection();

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
      if (_jp_ismc || _jp_skipl2l3res) {
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

      // Correct jets
      if (_jp_redojes)
        p4 *= jtjesnew[jetidx];
      jte[jetidx] = p4.E();
      jtpt[jetidx] = p4.Pt();
      jteta[jetidx] = p4.Eta();
      jtphi[jetidx] = p4.Phi();
      jty[jetidx] = p4.Rapidity();

      // Calculate gen level info
      if (_jp_ismc) {
        gp4.SetPxPyPzE(jtgenp4x[jetidx],jtgenp4y[jetidx],jtgenp4z[jetidx],jtgenp4t[jetidx]);
        jtgenpt[jetidx] = gp4.Pt();
        jtgeny[jetidx] = gp4.Rapidity();
        jtgeneta[jetidx] = gp4.Eta();
        jtgenphi[jetidx] = gp4.Phi();
      }

      // REMOVED: "Oversmear MC to match data"
      // TODO: Should we check, what MET is used nowadays?
      met = tools::oplus(mex, mey);
      metphi = atan2(mey, mex);
    } // for jetidx

    // Find leading jets (residual JEC may change ordering)
    // NOTE: for less than 3 jets, we input -1 on the place of the index
    for (int i = 0; i<3; ++i)
      jt3leads[i] = -1;
    for (int jetidx = 0; jetidx < njt; ++jetidx) {
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
    }

    if (_jp_ismc) {
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

    double ucx = -mex;
    double ucy = -mey;
    // Propagate jec to MET 1 and MET 2
    for (int jetidx = 0; jetidx != njt; ++jetidx) {

      // Only use jets with corr. pT>25 GeV to equalize data and MC thresholds
      if (jtpt[jetidx] > _jp_recopt && fabs(jteta[jetidx])<4.7) {

        // MET 1 stuff
        // Subtract uncorrected jet pT from met, put back corrected
        // Also add RC offset to keep PU isotropic
        // Remember that MET is negative vector sum
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

        // MET 2 stuff
        // Keep track of remaining pT in unclustered energy, i.e.
        // subtract jets from -MET to have the non-jet component
        // treat UE and PU underneath jets as unclustered in order
        // to keep the homogeneous
        double ue = 1.068 * jta[jetidx];
        ucx -= (l1corr * jtptu[jetidx] - ue) * cos(jtphi[jetidx]);
        ucy -= (l1corr * jtptu[jetidx] - ue) * sin(jtphi[jetidx]);
      }
    } // for jetidx
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

    if (njt!=0 && _pass) ++cnt["08njt"];

    _jetids.resize(njt);
    for (unsigned int jetid = 0; jetid != _jetids.size(); ++jetid)
      _jetids[jetid] = true;
    fillJetID(_jetids);

    if (njt!=0 && _jetids[0] && _pass) ++cnt["09jtid"];

    // Check if overweight PU event
    if (_jp_ismc && njt!=0 && _jetids[0] && _pass) {
      _pass = (jtpt[0] < 1.5*jtgenpt[0] || jtpt[0] < 1.5*pthat);
    }
    if (njt!=0 && _jetids[0] && _pass && _jp_ismc) ++cnt["10pthat"];

    // Equipped in fillBasics and fillRunHistos
    _evtid = (met < 0.4 * metsumet || met < 45.); // QCD-11-004

    if (_jp_etaphiexcl) {
      // Abort if one of the leading jets is in a difficult zone
      bool good0 = h2etaphiexcl->GetBinContent(h2etaphiexcl->FindBin(jteta[0],jtphi[0])) < 0;
      bool good1 = h2etaphiexcl->GetBinContent(h2etaphiexcl->FindBin(jteta[1],jtphi[1])) < 0;
      _pass = good0 and good1;
    }

    // Here can categorize events into different triggers, epochs,
    // topologies etc.
    // Eta and pT binning are handled in the fillBasic class
    if (_jp_doBasicHistos) {
      fillBasics("Standard");
    }

    if (_jp_doEtaHistos && _pass) {
      fillEtas("FullEta_Reco", jtpt, jteta, jtphi);
      if (_jp_ismc && _jp_doEtaHistosMcResponse) {
        fillEtas("FullEta_Gen", jtgenpt, jtgeneta, jtgenphi);
        fillMcHistos("FullEta_RecoPerGen_vReco", jtpt, jtgenpt, jtpt,    jteta,    jtphi);
        fillMcHistos("FullEta_RecoPerGen_vGen",  jtpt, jtgenpt, jtgenpt, jtgeneta, jtgenphi);
      }
    }

    // Run quality checks
    if (_jp_isdt && _jp_doRunHistos) {
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

  if (_jp_doRunHistos)   writeRunHistos();
  if (_jp_doEtaHistos)   writeEtas();
  if (_jp_ismc && _jp_doEtaHistos && _jp_doEtaHistosMcResponse) writeMcHistos();
  if (_jp_doBasicHistos) writeBasics(); // this needs to be last, output file closed

  // List bad runs
  cout << "Processed " << _totcounter << " events in total" << endl;
  cout << "Processed " << _trgcounter << " events passing "
       << " basic data quality and trigger cuts" << endl;
  cout << "(out of " << _evtcounter << " passing data quality cuts)" << endl;
  if (_badruns.size()!=0 || _badlums.size()!=0 || _nolums.size()!=0 ||
      _nbadevts_dup!=0 || _nbadevts_json!=0) {
    cout << "Found " << _badruns.size() << " bad runs:";
    for (auto runit = _badruns.begin(); runit != _badruns.end(); ++runit)
      cout << " " << *runit;
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
  for (auto cit = cnt.begin(); cit != cnt.end(); ++cit)
    cout << Form("%s: %d (%1.1f%%)", cit->first.c_str(), cit->second, 100. * cit->second / max(1, cnt["01all"])) << endl;
  cout << endl;

  // Report LS actually used for Jet15U in the analysis
  // (not necessarily containing any Jet15U triggers, though)
  cout << "Reporting JetMETTau LS in fillhistos.json" << endl;
  ofstream fout("fillhistos.json", ios::out);
  for (auto lumit = _jt15lums.begin(); lumit != _jt15lums.end(); ++lumit)
    fout << lumit->first << " " << lumit->second << endl;
  if (_jp_dojson) {
    cout << "Reporting LS marked newly bad in fillhistos.json.bad" << endl;
    ofstream fout2("fillhistos.json.bad", ios::out);
    for (auto jsit = _badjson.begin(); jsit != _badjson.end(); ++jsit)
      fout2 << jsit->first << " " << jsit->second << endl;
  } // _jp_dojson

  stop.Stop();
  cout << "Processing used " << stop.CpuTime() << "s CPU time ("
       << stop.CpuTime()/3600. << "h)" << endl;
  cout << "Processing used " << stop.RealTime() << "s real time ("
       << stop.RealTime()/3600. << "h)" << endl;
  cout << endl << endl;

  delete ferr;
  if (_jp_ismc || !_jp_useIOV) {
    delete _JEC;
    delete _L1RC;
    delete _jecUnc;
  }
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
  TFile *f = (_outfile ? _outfile : new TFile(Form("output-%s-1.root",_jp_type.c_str()), "RECREATE"));
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
  //TDirectory *topdir = gDirectory;
  TDirectory *topdir = f->GetDirectory(name.c_str()); assert(topdir);
  topdir->cd();

  // Pseudorapidity bins + HF + barrel
  double etas[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 0., 1.3};
  const int ny = sizeof(etas)/sizeof(etas[0])-1;

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
  if (_jp_isdt || _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {
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
  for (int etaidx = 0; etaidx != ny; ++etaidx) {

    if (etas[etaidx+1] > etas[etaidx]) { // create real bins only

      // subdirectory for rapidity bin
      const char *yname = Form("Eta_%1.1f-%1.1f", etas[etaidx], etas[etaidx+1]);
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
        basicHistos *h = new basicHistos(dir, trg, "", etas[etaidx], etas[etaidx+1], pttrg[trg],
                                         pt[trg].first, pt[trg].second, triggers[j]=="mc");
        _histos[name].push_back(h);
      } // for j
    } // real bin
  } // for etaidx

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
  for (unsigned int histidx = 0; histidx != _histos[name].size(); ++histidx)
    fillBasic(_histos[name][histidx]);
}


// Fill basic histograms after applying pt, y cuts
void fillHistos::fillBasic(basicHistos *h)
{
  assert(h);
  h->hpttmp->Reset();
  h->hpt_tmp->Reset();
  if (h->ismcdir) {
    assert(h->hpt_g0_tmp);
    h->hpt_g0_tmp->Reset();
  }

  _w = _w0 * _wt[h->trigname];
  assert(_w);

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  // Luminosity information
  if (_jp_isdt && h->lums[run][lbn]==0) {
    double lum = _lums[run][lbn];
    double lum2 = _lums2[run][lbn];
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
  if (_jp_ismc)
    h->hlumi_vstrpu->Fill(trpu, _w);

  if (_jt_debug) {
    if (h == _histos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit)
        cout << *trgit << ", ";
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired and passed
  if (!fired || !_pass) return;

  if (_jt_debug) cout << Form("Subdirectory Eta_%1.1f-%1.1f/%s",h->etamin,h->etamax,h->trigname.c_str()) << endl;
  if (_jt_debug) cout << "Calculate and fill dijet mass" << endl << flush;

  if (h->ismcdir) h->hpthat->Fill(pthat, _w);
  if (h->ismcdir) h->hpthatnlo->Fill(pthat);

  int i0 = jt3leads[0];
  int i1 = jt3leads[1];
  int i2 = jt3leads[2];
  if (_evtid and i0>=0 and _jetids[i0] and jtpt[i0]>_jp_recopt) { // First leading jet
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
        assert(h->hdjmass);
        h->hdjmass->Fill(djmass, _w);
        assert(h->hdjmass0);
        h->hdjmass0->Fill(djmass, _w);
        assert(h->pdjmass_ptratio);
        h->pdjmass_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
        assert(h->pdjmass0_ptratio);
        h->pdjmass0_ptratio->Fill(djmass, _j1.Pt()/_j2.Pt(), _w);
      }
      //} Dijet mass


      //{ Tag & probe hoods: Tag in barrel and fires trigger, probe in eta bin unbiased
      if (_jt_debug) cout << "Calculate and fill dijet balance" << endl << flush;

      double dphi = delta_phi(jtphi[i0], jtphi[i1]);

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
              double mpf = met*cos(delta_phi(metphi2,jtphi[itag]))/2;
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
              double metstuff = met2 * cos(delta_phi(metphi2, phiprobe));
              if (pttag >= h->ptmin and pttag < h->ptmax)
                h->hmpfx->Fill(1 + metstuff / pttag, _w);
              h->pmpfx->Fill(pttag, 1 + metstuff / pttag, _w);
              if (ptprobe >= h->ptmin and ptprobe < h->ptmax)
                h->hmpfy->Fill(1 + metstuff / ptprobe, _w);
              h->pmpfy->Fill(ptprobe, 1 + metstuff / ptprobe, _w);
              if (ptave >= h->ptmin and ptave < h->ptmax)
                h->hmpfz->Fill(1 + metstuff / ptave, _w);
              h->pmpfz->Fill(ptave, 1 + metstuff / ptave, _w);

              //{ Composition vs pt tag pt
              // Fractions vs pt: we do pt selection later in combineHistos
              assert(h->pncandtp);    h->pncandtp->Fill(pttag, jtn[iprobe], _w);
              assert(h->pnchtp);      h->pnchtp->Fill(pttag, jtnch[iprobe], _w);
              assert(h->pnnetp);      h->pnnetp->Fill(pttag, jtnne[iprobe], _w);
              assert(h->pnnhtp);      h->pnnhtp->Fill(pttag, jtnnh[iprobe], _w);
              assert(h->pncetp);      h->pncetp->Fill(pttag, jtnce[iprobe], _w);
              assert(h->pnmutp);      h->pnmutp->Fill(pttag, jtnmu[iprobe], _w);
              assert(h->pchftp);      h->pchftp->Fill(pttag, jtchf[iprobe], _w);
              assert(h->pneftp);      h->pneftp->Fill(pttag, jtnef[iprobe], _w);
              assert(h->pnhftp);      h->pnhftp->Fill(pttag, jtnhf[iprobe], _w);
              assert(h->pceftp);      h->pceftp->Fill(pttag, jtcef[iprobe], _w);
              assert(h->pmuftp);      h->pmuftp->Fill(pttag, jtmuf[iprobe], _w);
              assert(h->pbetatp);     h->pbetatp->Fill(pttag, jtbeta[iprobe], _w);
              assert(h->pbetastartp); h->pbetastartp->Fill(pttag, jtbetastar[iprobe], _w);

              if (pttag > h->ptmin and pttag < h->ptmax) { // Did tag fire the trigger?
                // The distributions of event counts per fraction
                assert(h->hncandtp);    h->hncandtp->Fill(jtn[iprobe], _w);
                assert(h->hnchtp);      h->hnchtp->Fill(jtnch[iprobe], _w);
                assert(h->hnnetp);      h->hnnetp->Fill(jtnne[iprobe], _w);
                assert(h->hnnhtp);      h->hnnhtp->Fill(jtnnh[iprobe], _w);
                assert(h->hncetp);      h->hncetp->Fill(jtnce[iprobe], _w);
                assert(h->hnmutp);      h->hnmutp->Fill(jtnmu[iprobe], _w);
                assert(h->hchftp);      h->hchftp->Fill(jtchf[iprobe], _w);
                assert(h->hneftp);      h->hneftp->Fill(jtnef[iprobe], _w);
                assert(h->hnhftp);      h->hnhftp->Fill(jtnhf[iprobe], _w);
                assert(h->hceftp);      h->hceftp->Fill(jtcef[iprobe], _w);
                assert(h->hmuftp);      h->hmuftp->Fill(jtmuf[iprobe], _w);
                assert(h->hbetatp);     h->hbetatp->Fill(jtbeta[iprobe], _w);
                assert(h->hbetastartp); h->hbetastartp->Fill(jtbetastar[iprobe], _w);

                // Fractions vs number of primary vertices
                assert(h->pncandtp_vsnpv);    h->pncandtp_vsnpv->Fill(npvgood, jtn[iprobe], _w);
                assert(h->pnchtp_vsnpv);      h->pnchtp_vsnpv->Fill(npvgood, jtnch[iprobe], _w);
                assert(h->pnnetp_vsnpv);      h->pnnetp_vsnpv->Fill(npvgood, jtnne[iprobe], _w);
                assert(h->pnnhtp_vsnpv);      h->pnnhtp_vsnpv->Fill(npvgood, jtnnh[iprobe], _w);
                assert(h->pncetp_vsnpv);      h->pncetp_vsnpv->Fill(npvgood, jtnce[iprobe], _w);
                assert(h->pnmutp_vsnpv);      h->pnmutp_vsnpv->Fill(npvgood, jtnmu[iprobe], _w);
                assert(h->pchftp_vsnpv);      h->pchftp_vsnpv->Fill(npvgood, jtchf[iprobe], _w);
                assert(h->pneftp_vsnpv);      h->pneftp_vsnpv->Fill(npvgood, jtnef[iprobe], _w);
                assert(h->pnhftp_vsnpv);      h->pnhftp_vsnpv->Fill(npvgood, jtnhf[iprobe], _w);
                assert(h->pceftp_vsnpv);      h->pceftp_vsnpv->Fill(npvgood, jtcef[iprobe], _w);
                assert(h->pmuftp_vsnpv);      h->pmuftp_vsnpv->Fill(npvgood, jtmuf[iprobe], _w);
                assert(h->pbetatp_vsnpv);     h->pbetatp_vsnpv->Fill(npvgood, jtbeta[iprobe], _w);
                assert(h->pbetastartp_vsnpv); h->pbetastartp_vsnpv->Fill(npvgood, jtbetastar[iprobe], _w);

                // Fractions vs true pileup
                assert(h->pchftp_vstrpu);      h->pchftp_vstrpu->Fill(trpu, jtchf[iprobe], _w);
                assert(h->pneftp_vstrpu);      h->pneftp_vstrpu->Fill(trpu, jtnef[iprobe], _w);
                assert(h->pnhftp_vstrpu);      h->pnhftp_vstrpu->Fill(trpu, jtnhf[iprobe], _w);
                assert(h->pceftp_vstrpu);      h->pceftp_vstrpu->Fill(trpu, jtcef[iprobe], _w);
                assert(h->pmuftp_vstrpu);      h->pmuftp_vstrpu->Fill(trpu, jtmuf[iprobe], _w);
                assert(h->pbetatp_vstrpu);     h->pbetatp_vstrpu->Fill(trpu, jtbeta[iprobe], _w);
                assert(h->pbetastartp_vstrpu); h->pbetastartp_vstrpu->Fill(trpu, jtbetastar[iprobe], _w);
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

      h->hr21->Fill(has2  ? jtpt[1] / jtpt[0] : 0.);
      h->hr31->Fill(has3  ? jtpt[2] / jtpt[0] : 0.);
      h->hr32->Fill(has32 ? jtpt[2] / jtpt[1] : 0.);
      if (has2)  h->pr21->Fill(npvgood, has2 ? jtpt[i1] / jtpt[i0] : 0.);
      if (has3)  h->pr31->Fill(npvgood, has3 ? jtpt[i2] / jtpt[i0] : 0.);
      if (has32) h->pr32->Fill(npvgood, has3 ? jtpt[i2] / jtpt[i1] : 0.);
      h->px21->Fill(npvgood, has2  ? 1 : 0);
      h->px31->Fill(npvgood, has3  ? 1 : 0);
      h->px32->Fill(npvgood, has32 ? 1 : 0);
    } // Jet quality stats
  } // First leading jet

  // retrieve event-wide variables
  double dphi = (i1>=0 ? delta_phi(jtphi[i0], jtphi[i1]) : 0.);
  double dpt = (i1>=0 ? fabs(jtpt[i0]-jtpt[i1])/(jtpt[i0]+jtpt[i1]) : 0.999);
  //double met = this->met;
  //double metphi = this->metphi;
  double sumet = this->metsumet;

  if (_jt_debug) cout << "Entering jet loop" << endl << flush;
  for (int jetidx = 0; jetidx != njt; ++jetidx) {
    if (_jt_debug) cout << "Loop over jet " << jetidx << "/" << njt << endl << flush;

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
      if (_jt_debug) {
        cout << "..." << h->trigname << " | " << " index " << jetidx << "/" << njt
          << " jet pt: " << pt << " eta : " << eta << " id: " << id << " jec: " << jec << endl;
        cout << "...evt id: " << _evtid << " weight: " << _w << " met: " << met << " metsumet: " << metsumet << endl;
      }

      assert(h->hpt_noid); h->hpt_noid->Fill(pt, _w);
      assert(h->hpt_nojetid); if (_evtid) h->hpt_nojetid->Fill(pt, _w);
      assert(h->hpt_noevtid); if (id) h->hpt_noevtid->Fill(pt, _w);
      // Same versus generator pT as MC extra
      // to decouple efficiency from JEC and JER
      if (h->ismcdir) {
        h->hpt_noid_g->Fill(ptgen, _w);
        if (_evtid) h->hpt_nojetid_g->Fill(ptgen, _w);
        if (id)    h->hpt_noevtid_g->Fill(ptgen, _w);
      }
    } // Jet in eta range

    // REMOVED: "Debugging JEC"

    if (pt>_jp_recopt) { // pt visible
      // Flags for studying gen eta vs reco eta effects
      bool mcstudy = h->ismcdir and jtgenr[jetidx] < 0.25;
      bool mcgenetarange = fabs(jtgeneta[jetidx]) >= h->etamin && fabs(jtgeneta[jetidx]) < h->etamax;

      if (mcstudy and mcgenetarange) // GenJets matched to any reco jets in any events
        h->hpt_gg0->Fill(ptgen, _w);

      if (_evtid and id) { // id OK
        if (mcstudy and mcgenetarange) // GenJets matched to good reco jets in good events
          h->hpt_gg->Fill(ptgen, _w);

        if (etarange) { // Correct jet eta range
          if (_jt_debug) cout << "..jec uncertainty" << endl << flush;

          // Get JEC uncertainty
          double unc = 0.01; // default for MC
          if (_jecUnc) {
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

          if (_jt_debug) cout << "..raw spectrum" << endl << flush;

          // REMOVED: "For trigger efficiency"

          if (_jp_doqgl) {
            // new histograms for quark/gluon study (Ozlem)
            // 1. open the previous output-MC-1_iteration1.root in the beginning
            // 2. get the hqgl_q and hqgl_g for each eta bin, store in array
            // 3. find correct hqgl_q and hqgl_g from array (normalized)
            // 4. calculate probg = g / (q+g)
            //double probg = 1.-qgl[jetidx]; // Later, map QGL to probability using hqgl_x
            assert(h3probg);
            double probg = h3probg->GetBinContent(h3probg->FindBin(jty[jetidx],jtpt[jetidx],qgl[jetidx]));
            if (probg>=0 and probg<=1) {
              assert(h->hgpt);
              h->hgpt->Fill(pt,_w*probg);
              assert(h->hgpt0);
              h->hgpt0->Fill(pt, _w*probg);

              assert(h->hqgl);
              h->hqgl->Fill(qgl[jetidx], _w);
              if (_jp_ismc) {
                assert(h->hqgl_g);
                assert(h->hqgl_q);
                bool isgluon = (fabs(partonflavor[jetidx]-21)<0.5);
                bool isquark = (fabs(partonflavor[jetidx])<7);
                assert(isgluon || isquark);
                if (isgluon) h->hqgl_g->Fill(qgl[jetidx], _w);
                if (isquark) h->hqgl_q->Fill(qgl[jetidx], _w);

                // For data templates from scaling Pythia, see instructions at
                // https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#Systematics
                double x = qgl[jetidx];
                double wq =  -0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581;
                double wg = -55.7067*pow(x,7) + 113.218*pow(x,6) -21.1421*pow(x,5) -99.927*pow(x,4) + 
                            92.8668*pow(x,3) -34.3663*x*x + 6.27*x + 0.612992;
                if (isgluon) h->hqgl_dg->Fill(qgl[jetidx], _w*wg);
                if (isquark) h->hqgl_dq->Fill(qgl[jetidx], _w*wq);
              }
            } // probg quark/gluon
          }

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

          if (_jt_debug) cout << "..basic properties" << endl << flush;

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

          if (_jt_debug) cout << "..control plots of components" << endl << flush;

          // Composition stuff without T&P (according to triggers)
          assert(h->pncand); h->pncand->Fill(pt, jtn[jetidx], _w);
          assert(h->pnch); h->pnch->Fill(pt, jtnch[jetidx], _w);
          assert(h->pnne); h->pnne->Fill(pt, jtnne[jetidx], _w);
          assert(h->pnnh); h->pnnh->Fill(pt, jtnnh[jetidx], _w);
          assert(h->pnce); h->pnce->Fill(pt, jtnce[jetidx], _w);
          assert(h->pnmu); h->pnmu->Fill(pt, jtnmu[jetidx], _w);
          //
          assert(h->pchf); h->pchf->Fill(pt, jtchf[jetidx], _w);
          assert(h->pnef); h->pnef->Fill(pt, jtnef[jetidx], _w);
          assert(h->pnhf); h->pnhf->Fill(pt, jtnhf[jetidx], _w);
          assert(h->pcef); h->pcef->Fill(pt, jtcef[jetidx], _w);
          assert(h->pmuf); h->pmuf->Fill(pt, jtmuf[jetidx], _w);
          assert(h->pbeta); h->pbeta->Fill(pt, jtbeta[jetidx], _w);
          assert(h->pbetastar); h->pbetastar->Fill(pt, jtbetastar[jetidx], _w);

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

            if (_jt_debug) cout << "..control plots for topology" << endl << flush;

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
            // Component fractions
            h->hncand->Fill(jtn[jetidx], _w);
            h->hnch->Fill(jtnch[jetidx], _w);
            h->hnne->Fill(jtnne[jetidx], _w);
            h->hnnh->Fill(jtnnh[jetidx], _w);
            h->hnce->Fill(jtnce[jetidx], _w);
            h->hnmu->Fill(jtnmu[jetidx], _w);
            //
            h->hchf->Fill(jtchf[jetidx], _w);
            h->hnef->Fill(jtnef[jetidx], _w);
            h->hnhf->Fill(jtnhf[jetidx], _w);
            h->hcef->Fill(jtcef[jetidx], _w);
            h->hmuf->Fill(jtmuf[jetidx], _w);
            h->hbeta->Fill(jtbeta[jetidx], _w);
            h->hbetastar->Fill(jtbetastar[jetidx], _w);

            h->hyeta->Fill(TMath::Sign(y-eta,y), _w);
            h->hyeta2->Fill(y-eta, _w);
            h->hbetabetastar->Fill(jtbeta[jetidx], jtbetastar[jetidx], _w);
            h->hetaphi->Fill(eta, phi, _w);
          } // within trigger pT range

          // closure plots for JEC
          h->pdpt->Fill(pt, dpt, _w);
          h->pmpf->Fill(pt, 1 + met * cos(delta_phi(metphi, phi)) / pt, _w);
          h->pmpf1->Fill(pt, 1 + met1 * cos(delta_phi(metphi1, phi)) / pt, _w);
          h->pmpf2->Fill(pt, 1 + met2 * cos(delta_phi(metphi2, phi)) / pt, _w);

          if (mcstudy) { // MC extras
            if (_jt_debug)
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

        // MC: Filling outside of eta bin
        if (mcstudy and ptgen!=0 and jtjesnew[jetidx]!=0) {
          double r = (ptgen ? pt / ptgen : 0);
          double resp = r / jtjesnew[jetidx];

          // Response closure vs NPV
          if (r) {
            h->p3rvsnpv->Fill(ptgen, jteta[jetidx], npvgood, resp, _w);
            h->p3rvsnpvW->Fill(ptgen, fabs(jteta[jetidx]), npvgood, resp, _w);
          }
        } // if MC 
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

  if (_jp_ismc) { // Unbiased generator spectrum (for each trigger)
    if (_jt_debug) cout << "Truth loop:" << endl;
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      double geneta = gen_jteta[gjetidx];
      double genpt = gen_jtpt[gjetidx];
      if (fabs(geneta) >= h->etamin && fabs(geneta) < h->etamax) {
        if (_jt_debug)
          cout << "genjet " << gjetidx << "/" << gen_njt << " ptg="<<genpt << " etag="<<geneta << endl;

        h->hpt_g0tw->Fill(genpt, _w);
        // Ozlem: (gluon vs quark)
        if (gen_partonflavor[gjetidx]==21) h->hgpt_g0tw->Fill(genpt, _w);
        else h->hqpt_g0tw->Fill(genpt, _w);
      }
    }
  }
  //
  if (h->ismcdir) {
    // unfolding studies (Mikael)
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      double ygen = fabs(gen_jty[gjetidx]);

      for (int j = 0; j != njt; ++j) {
        //double yreco = fabs(jty[j]);
        bool id = (_jetids[j] && _evtid);
        if ((ygen >= h->etamin && ygen < h->etamax && gen_jtpt[gjetidx]>_jp_recopt) &&
            //(yreco >= h->etamin && yreco < h->etamax)
            (jtpt[j]>_jp_recopt && id)) {

          double dr = tools::oplus(delta_phi(gen_jtphi[gjetidx], jtphi[j]),
                                   fabs(gen_jteta[gjetidx] - jteta[j]));
          if (dr < 0.25) {
            h->mT->Fill(gen_jtpt[gjetidx], jtpt[j], _w);
            h->mTf->Fill(gen_jtpt[gjetidx], jtpt[j], _w);
            h->mTuw->Fill(gen_jtpt[gjetidx], jtpt[j]);
            h->mTfuw->Fill(gen_jtpt[gjetidx], jtpt[j]);
          }
        } // rapidity bin
      } // for j
    } // for gjetidx
    //
    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      double ygen = fabs(gen_jty[gjetidx]);
      if (ygen >= h->etamin && ygen < h->etamax && gen_jtpt[gjetidx]>_jp_recopt) {
        h->mx->Fill(gen_jtpt[gjetidx], _w);
        h->mxf->Fill(gen_jtpt[gjetidx], _w);
        h->mxuw->Fill(gen_jtpt[gjetidx]);
        h->mxfuw->Fill(gen_jtpt[gjetidx]);
      }
    } // for gjetidx
    //
    for (int j = 0; j != njt; ++j) {
      double yreco = fabs(jty[j]);
      bool id  = (_jetids[j] && _evtid && jtpt[j]>_jp_recopt);
      if (yreco >= h->etamin && yreco < h->etamax && id) {
        h->my->Fill(jtpt[j], _w);
        h->myf->Fill(jtpt[j], _w);
        h->myuw->Fill(jtpt[j]);
        h->myfuw->Fill(jtpt[j]);
      }
    } // for j

    for (int gjetidx = 0; gjetidx != gen_njt; ++gjetidx) {
      double etagen = gen_jty[gjetidx];
      if (fabs(etagen) >= h->etamin && fabs(etagen) < h->etamax) {
        h->hpt_g0->Fill(gen_jtpt[gjetidx], _w);
        // Ozlem: (gluon vs quark)
        if (gen_partonflavor[gjetidx]==21) h->hgpt_g0->Fill(gen_jtpt[gjetidx], _w);
        else h->hqpt_g0->Fill(gen_jtpt[gjetidx], _w);

        assert(h->hpt_g0_tmp);
        h->hpt_g0_tmp->Fill(gen_jtpt[gjetidx]);
      }
    } // for gjetidx
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

  for (auto histit = _histos.begin(); histit != _histos.end(); ++histit) {
    for (auto histidx = 0u; histidx != histit->second.size(); ++histidx) {
      // Luminosity information
      basicHistos *h = histit->second[histidx];
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
  *ferr << "writeBasic() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeBasic


// Initialize eta histograms for trigger bins
void fillHistos::initEtas(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initEtas("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile : new TFile(Form("output-%s-1.root",_jp_type), "RECREATE"));
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
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
  if (_jp_isdt || _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {
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
    etaHistos *h = new etaHistos(dir, trg,pttrg[trg],pt[trg].first, pt[trg].second);
    _etahistos[name].push_back(h);
  } // for j

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initEtas("<<name<<") finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initEtas


// Loop over basic histogram containers to fill all
void fillHistos::fillEtas(string name, Float_t* _pt, Float_t* _eta, Float_t* _phi)
{
  for (auto histidx = 0u; histidx != _etahistos[name].size(); ++histidx)
    fillEta(_etahistos[name][histidx], _pt, _eta, _phi);
}


// Fill basic histograms after applying pt, y cuts
void fillHistos::fillEta(etaHistos *h, Float_t* _pt, Float_t* _eta, Float_t* _phi)
{
  assert(h);

  _w = _w0 * _wt[h->trigname];
  assert(_w);

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  if (_jt_debug) {
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
  if (_evtid and i0>=0 and _jetids[i0] and jtpt[i0]>_jp_recopt and i1>=0 and _jetids[i1] and jtpt[i1]>_jp_recopt) { // Leading jets
    double dphi = delta_phi(jtphi[i0], jtphi[i1]);
    if (dphi > 2.7) { // Back-to-back condition
      double pt3 = ((i2>=0 and jtpt[i2]>_jp_recopt) ? jtpt[i2] : 0.);
      double ptave = 0.5 * (jtpt[i0] + jtpt[i1]);
      double alpha = pt3/ptave;

      for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
        int itag = jt3leads[itag_lead];
        int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
        double etatag = jteta[itag];
        double etaprobe = jteta[iprobe];
        double pttag = jtpt[itag];
        double ptprobe = _pt[iprobe];
        double alphatp = pt3/pttag;
        if (fabs(etatag) < 1.3) {
          if (_jp_do3dHistos) {
            double asymm = (ptprobe - pttag)/(2*ptave);
            //double asymmtp = (ptprobe - pttag)/(2*pttag);
            double mpf = met*cos(delta_phi(metphi2,_phi[itag]))/(2*ptave);
            //double mpftp = met2*cos(delta_phi(metphi2,_phi[itag]))/(2*pttag);
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
            }
          }

          // for composition vs eta (Ozlem) || jetidx is the probe
          if (alphatp < 0.3 and pttag > h->ptmin and pttag < h->ptmax) { // Alpha and trigger
            assert(h->pchftp_vseta); h->pchftp_vseta->Fill(etaprobe, jtchf[iprobe], _w);
            assert(h->pneftp_vseta); h->pneftp_vseta->Fill(etaprobe, jtnef[iprobe], _w);
            assert(h->pnhftp_vseta); h->pnhftp_vseta->Fill(etaprobe, jtnhf[iprobe], _w);
            assert(h->pceftp_vseta); h->pceftp_vseta->Fill(etaprobe, jtcef[iprobe], _w);
            assert(h->pmuftp_vseta); h->pmuftp_vseta->Fill(etaprobe, jtmuf[iprobe], _w);
            assert(h->pbetatp_vseta); h->pbetatp_vseta->Fill(etaprobe, jtbeta[iprobe], _w);
            assert(h->pbetastartp_vseta); h->pbetastartp_vseta->Fill(etaprobe, jtbetastar[iprobe], _w);
          } // select pt bin for profiles vseta
        } // etatag < 1.3
      } // tag & probe
    } // dphi > 2.7
  } // ids
} // fillEta


// Write and delete histograms
void fillHistos::writeEtas()
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeEtas():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (auto histit : _etahistos) {
    for (unsigned int histidx = 0; histidx != histit.second.size(); ++histidx) {
      etaHistos *h = histit.second[histidx];
      delete h;
    } // for histidx
  } // for histit

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeEtas() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeEtas


// Initialize eta histograms for trigger bins
void fillHistos::initMcHistos(string name)
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initMcHistos("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = (_outfile ? _outfile : new TFile(Form("output-%s-1.root",_jp_type), "RECREATE"));
  assert(f && !f->IsZombie());
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
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
  if (_jp_isdt || _jp_domctrigsim) {
    // This is done both for data and MC, because why not?
    for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {
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
    mcHistos *h = new mcHistos(dir, trg);
    _mchistos[name].push_back(h);
  } // for j

  _outfile = f;
  curdir->cd();

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initMcHistos("<<name<<") finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // initMcHistos


// Loop over basic histogram containers to fill all
void fillHistos::fillMcHistos(string name,  Float_t* _recopt, Float_t* _genpt,
                              Float_t* _pt, Float_t* _eta,    Float_t* _phi)
{
  for (auto histidx = 0u; histidx != _mchistos[name].size(); ++histidx)
    fillMcHisto(_mchistos[name][histidx], _recopt, _genpt, _pt, _eta, _phi);
}


// Fill basic histograms after applying pt, y cuts
void fillHistos::fillMcHisto(mcHistos *h,  Float_t* _recopt,  Float_t* _genpt,
                             Float_t* _pt, Float_t* _eta,     Float_t* _phi)
{
  assert(h);

  _w = _w0 * _wt[h->trigname];
  assert(_w);

  bool fired = (_trigs.find(h->trigname)!=_trigs.end());

  if (_jt_debug) {
    if (h == _mchistos.begin()->second[0]) {
      cout << "Triggers size: " << _trigs.size() << endl;
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit)
        cout << *trgit << ", ";
      cout << "(" << h->trigname << ")" << endl;
    }
  }

  // check if required trigger fired
  if (!fired) return;

  // Calculate and fill dijet balance histograms
  if (njt>=2 && _evtid && delta_phi(_phi[0],_phi[1])>2.8
      && _jetids[0] && _jetids[1] && _pt[0]>_jp_recopt && _pt[1]>_jp_recopt) {
    // Two leading jets
    for (int itag = 0; itag<2; ++itag) {
      int iprobe = (itag==0 ? 1 : 0);
      double etatag = fabs(_eta[itag]); // This is usable without sign
      double etaprobe = _eta[iprobe];   // Here we need sign

      if (etatag < 1.3) {
        double pttag = _pt[itag];
        double ptprobe = _pt[iprobe];
        double pt3 = (njt>2 ? _pt[2] : 0.);
        double ptave = 0.5 * (pttag + ptprobe);
        assert(ptave);
        if (_jp_do3dHistos) {
          double alpha = pt3/ptave;
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
            }
            if (alpha<alphasel) {
              h->hdjasymm     [alphaidx]->Fill(ptave, etaprobe, asymm,        _w);
              h->hdjresp_tag  [alphaidx]->Fill(ptave, etaprobe, ptresp_tag,   _w);
              h->hdjresp_probe[alphaidx]->Fill(ptave, etaprobe, ptresp_probe, _w);
            }
          }
        }
      } // etatag < 1.3
    } // for itag (two leading jets)
  } // two or more jets, phase space
} // fillEta


// Write and delete histograms
void fillHistos::writeMcHistos()
{
  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeMcHistos():" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  for (auto histit : _mchistos) {
    for (auto histidx = 0u; histidx != histit.second.size(); ++histidx) {
      mcHistos *h = histit.second[histidx];
      delete h;
    } // for histidx
  } // for histit

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "writeMcHistos() finished:" << endl << flush;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
               info.fMemTotal, info.fMemUsed, info.fMemFree,
               info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;
} // writeMcHistos


// Initialize basic histograms for trigger and eta bins
void fillHistos::initRunHistos(string name, double etamin, double etamax) {

  // Report memory usage to avoid malloc problems when writing file
  *ferr << "initRunHistos("<<name<<"):" << endl << flush;
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  *ferr << Form("MemInfo(Tot:%d, Used:%d, Free:%d, Stot:%d, SUsed:%d, SFree:%d",
                info.fMemTotal, info.fMemUsed, info.fMemFree,
                info.fSwapTotal, info.fSwapUsed, info.fSwapFree) << endl<<flush;

  TDirectory *curdir = gDirectory;

  // open file for output
  TFile *f = _outfile ? _outfile : new TFile(Form("output-%s-1.root",_jp_type), "RECREATE");
  assert(f && !f->IsZombie());
  //assert(f->mkdir(name.c_str()));
  f->mkdir(name.c_str());
  assert(f->cd(name.c_str()));
  //TDirectory *dir = gDirectory;
  TDirectory *dir = f->GetDirectory(name.c_str()); assert(dir);
  dir->cd();

  runHistos *h = new runHistos(dir, etamin, etamax);
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
  if (_jp_isdt && h->lums[run][lbn]==0) {

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
          *ferr << "Prescale not found for trigger " << t
                << " run " << run << endl << flush;
          assert(false);
        }
      }
      else
        prescale = _prescales[t][run]; //assert(prescale);
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
        h->t_trgtp[t][run] = 0;
        h->c_chftp[t][run] = 0;
        h->c_neftp[t][run] = 0;
        h->c_nhftp[t][run] = 0;
        h->c_betastartp[t][run] = 0;

        for (unsigned int j = 0; j != h->trg.size(); ++j) {
          string const& t2 = h->trg[j];
          h->p_trgpair[t+t2][run] = 0;
        } // for j
      } // for trgidx
    } // new run
  }

  double dphi = (njt>=2 ? delta_phi(jtphi[0], jtphi[1]) : 0.);
  double pt3 = (njt>=3 ? jtpt[2] : 0.);

  for (int jetidx = 0; jetidx != njt; ++jetidx) {

    double pt = jtpt[jetidx];
    double y = jty[jetidx];
    double eta = jteta[jetidx];

    if (h->etamin <= fabs(eta) and fabs(eta) < h->etamax and _pass and _jetids[jetidx] and _evtid) {
      for (auto trgit = _trigs.begin(); trgit != _trigs.end(); ++trgit) {
        string const& t = *trgit;

        if (pt > 18.) ++h->p_trg[t][run];
        if (pt > h->pt[t]) {
          ++h->t_trg[t][run]; // unweighted events
          h->tw_trg[t][run] += _prescales[t][run]; // prescale weighted events
          h->npv_trg[t][run] += npv;
          h->npvgood_trg[t][run] += npvgood;
          h->c_chf[t][run] += jtchf[jetidx];
          h->c_nef[t][run] += jtnef[jetidx];
          h->c_nhf[t][run] += jtnhf[jetidx];
          h->c_betastar[t][run] += jtbetastar[jetidx];
        }

        int itag = (jetidx==0 ? 1 : 0);
        if (jetidx<2 and dphi > 2.7 and pt3 < jtpt[itag] and fabs(jty[itag]) < 1.3 and
            jtpt[itag] > h->pt[t] and _jetids[itag]) {
          ++h->t_trgtp[t][run];
          h->c_chftp[t][run] += jtchf[jetidx];
          h->c_neftp[t][run] += jtnef[jetidx];
          h->c_nhftp[t][run] += jtnhf[jetidx];
          h->c_betastartp[t][run] += jtbetastar[jetidx];
        }

        for (auto trgit2 = _trigs.begin(); trgit2 != _trigs.end(); ++trgit2) {
          string const& t2 = *trgit2;
          if (t!=t2) ++h->p_trgpair[t+t2][run];
        } // for trgit2
      } // for trgit
    } // conditions
  } // for jetidx
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

  for (auto histit = _runhistos.begin(); histit != _runhistos.end(); ++histit) {
    runHistos *h = histit->second;
    delete h;
  } // for histit

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

  for (int jetidx = 0; jetidx != njt; ++jetidx) {

    id[jetidx] = ((fabs(jty[jetidx])<2.5 ? jtidtight[jetidx] : jtidloose[jetidx]));
             //&& (_jp_ismc || jtbeta[jetidx]!=0));
             //&& (_jp_ismc || (1-jtbetastar[jetidx])>0.5));

    if (_jp_doECALveto) {
      assert(_ecalveto);
      int ibin = _ecalveto->FindBin(jteta[jetidx],jtphi[jetidx]);
      id[jetidx] = (id[jetidx] && _ecalveto->GetBinContent(ibin)==0);
    }
  }

} // fillJetID


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
    if (_jt_debug)
      cout << "\"" << rn << "\": " << flush;

    while (file.get(c) && c==' ') {};
    if (_jt_debug)
      cout << c << flush; assert(c=='[');
    ++nrun;

    bool endrun = false;
    while (!endrun && file >> s >> s2 &&
           sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (_jt_debug)
        cout << "["<<ls1<<","<<ls2<<"]"<<s1 << flush;

      for (int ls = ls1; ls != ls2+1; ++ls) {
        //assert(_json[rn].find(ls)==_json[rn].end()); // ok if 2 JSON files
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (!endrun && s2!=",") {
        if (_jt_debug)
          cout<<"s1: "<<s2<<endl<<flush; assert(s2==",");
      }
    } // while ls
    if (_jt_debug)
      cout << endl;

    if (s2=="]}") continue;
    else if (s2!="],") {
      if (_jt_debug)
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
  for (int runidx = 0; runidx != ngoodruns; ++runidx)
    goodruns.insert(a_goodruns[runidx]);

  set<pair<int, int> > nolums;

  for (auto runit = goodruns.begin(); runit != goodruns.end(); ++runit)
    cout << *runit << ", ";
  cout << endl;

  ifstream f(filename, ios::in);
  assert(f.is_open());
  float secLS = 2.3310e+01;
  string s;
  int rn, fill, ls, ifoo;
  float del, rec, avgpu, ffoo;
  char sfoo[512];
  assert(getline(f, s, '\n'));
  cout << endl << "string: " << s << " !" << endl << flush;

  // HOX: the lumi file format has been changing:
  //bool v1 = (s==string("run,ls,delivered,recorded"));
  //bool v2 = (s==string("Run,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)"));
  //bool v3 = (s==string("Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)"));
  //assert(v1 || v2 || v3);
  //assert (s=="#Data tag : online , Norm tag: None");
  assert (s=="#Data tag : v1 , Norm tag: None");

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
  // REMOVED: "while ((v1 && f >> s &&" etc.
  while (//getline(f, s, '\r') &&
        getline(f, s, '\n') &&
        //f >> s &&
        (sscanf(s.c_str(),"%d:%d,%d:%d,%d/%d/%d %d:%d:%d,STABLE BEAMS,"
        "%f,%f,%f,%f,%s", &rn,&fill,&ls,&ifoo, &ifoo,&ifoo,&ifoo,
        &ifoo,&ifoo,&ifoo, &ffoo,&del,&rec,&avgpu,sfoo)==15 ||
        (skip=true))) {

    if (_jt_debug) {
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

    _avgpu[rn][ls] = avgpu; // * 69000. / 78400.; // brilcalc --minBiasXsec patch
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

  pumc = (TH1F*)fpumc->Get("pileupmc"); assert(pumc);

  // Normalize
  pumc->Scale(1./pumc->Integral());

  // For data, load each trigger separately
  for (unsigned int itrg = 0 ; itrg != _triggers.size(); ++itrg) {
    const char *t = _triggers[itrg].c_str();
    pudist[t] = (TH1D*)fpudist->Get(Form("%s",t)); assert(pudist[t]);
    pudist[t]->Scale(1./pudist[t]->Integral());
  }
  // REMOVED: "data with only one histo:"

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
    if (_jt_debug) cout << " run " << run << " ls " << ls << ": ";
    //ss >> run;
    //cout << " run " << run << ": ";

    int itrg(0);
    for (unsigned int itrg = 0; ss >> pre; ++itrg) {
      //cout << "trg" << trgs[itrg] << " run " << run << " ls " << ls;
      assert(itrg!=trgs.size());
      _premap[trgs[itrg]][run][ls] = pre;
      if (_jt_debug) cout << pre << "/" << trgs[itrg] << " ";
    }
    if (_jt_debug) cout << endl;
  }
} // loadPrescales


void fillHistos::loadECALveto(const char *file)
{
  cout << "Processing loadECALveto(\"" << file << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  TFile *fe = new TFile(file, "READ");
  assert(fe && !fe->IsZombie());

  _ecalveto = (TH2F*)fe->Get("ecalveto"); assert(_ecalveto);

  curdir->cd();
} // loadECALveto


// Update the available trigger types for each new tree
void fillHistos::getTriggers()
{
  TH1F *triggers = (TH1F*)fChain->GetCurrentFile()->Get("ak4/TriggerNames");
  TAxis *xax = triggers->GetXaxis();

  std::regex pfjet("HLT_PFJet([0-9]*)_v[0-9]*");

  _availTrigs.clear();
  for (int trgidx = xax->GetFirst(); trgidx <= xax->GetLast(); ++trgidx) {
    string trgName = xax->GetBinLabel(trgidx);
    if (trgName.compare("")==0) // Ignore empty places on x-axis
      continue;
    string trigger = ""; // HLT_PFJet are given non-empty trigger names
    if (std::regex_match(trgName,pfjet))
      trigger=std::regex_replace(trgName, pfjet, "jt$1", std::regex_constants::format_no_copy);
    _availTrigs.push_back(trigger);
  }
} // getTriggers
