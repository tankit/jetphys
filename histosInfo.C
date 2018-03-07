// Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: April 19, 2010
// Updated: June 2, 2015
// Updated: Aug 31, 2016
// Further updates: see git log

#define histosInfo_cxx
#include "histosInfo.h"

void histosInfo::Loop()
{
  if (fChain == 0) return;
  Long64_t _nentries = fChain->GetEntriesFast();
  Long64_t _ntot = fChain->GetEntries();
  Long64_t nskip = _jp_nskip;//0;
  _nentries = (_jp_nentries==-1 ? _ntot-nskip : min(_ntot-nskip, _jp_nentries));
  assert(nskip+_nentries);

  map<string, int> _cnt; // efficiency counters

  ferr = new ofstream(Form("reports/histosInfo-%s.log",_jp_type),ios::out);

  MemInfo_t info;
  gSystem->GetMemInfo(&info);

  // Initialize _pthatweight. It will be loaded for each tree separately.
  if (_jp_pthatbins)
    _pthatweight = 0;

  if (_jp_quick) { // Activate only some branches
    fChain->SetBranchStatus("*",0);
    // Luminosity calculation
    if (_jp_ismc) fChain->SetBranchStatus("EvtHdr_.mWeight",1); // weight
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mRun",1); // run
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mEvent",1); // evt
    if (_jp_isdt) fChain->SetBranchStatus("EvtHdr_.mLumi",1); // lbn

    fChain->SetBranchStatus("PFMet_.et_",1); // met
    fChain->SetBranchStatus("PFMet_.phi_",1); // metphi
    fChain->SetBranchStatus("PFMet_.sumEt_",1); // metsumet

    fChain->SetBranchStatus("TriggerDecision_",1);
    fChain->SetBranchStatus("L1Prescale_",1);
    fChain->SetBranchStatus("HLTPrescale_",1);

    if (_jp_ismc) fChain->SetBranchStatus("EvtHdr_.mTrPu",1); // trpu
  } else {
    fChain->SetBranchStatus("*",1);
  } // quick/slow

  // Load latest JSON selection
  if (_jp_isdt and _jp_dojson and !loadJSON(_jp_json)) {
    cout << "Issues loading the JSON file; aborting..." << endl;
    return;
  }

  // Load PU profiles for MC reweighing
  if (_jp_ismc and _jp_reweighPU and !loadPUProfiles(_jp_pudata, _jp_pumc)) {
    cout << "Issues loading the PU histograms for reweighting; aborting..." << endl;
    return;
  }

  // load luminosity tables (prescales now stored in event)
  if (_jp_isdt and _jp_dolumi and !loadLumi(_jp_lumifile)) {
    cout << "Issues loading the Lumi file; aborting..." << endl;
    return;
  }

  // Event loop
  TStopwatch stop;
  stop.Start();
  TDatime bgn;

  vector<Long64_t> infojentrys = {5000,10000,50000,100000,500000,1000000,5000000,
                                  10000000,20000000,30000000,40000000,50000000,60000000,70000000,80000000,90000000,
                                  100000000,200000000,300000000,400000000,500000000,600000000,700000000,800000000,900000000,1000000000};
  Long64_t nbytes = 0, nb = 0;

  map<string,TH1D*> hltTrigs;
  map<string,TH1D*> aftTrigs;
  for (int ID = -1; ID < int(_jp_notrigIDs); ++ID) {
    for (auto i = 0u; i < _jp_notrigs; ++i) {
      string iname = Form("%s_%d",_jp_triggers[i],ID);
      string uname = Form("aft%s_%d",_jp_triggers[i],ID);
      hltTrigs[iname] = new TH1D(iname.c_str(),iname.c_str(),_jp_maxpu,0,_jp_maxpu);
      aftTrigs[uname] = new TH1D(uname.c_str(),uname.c_str(),_jp_maxpu,0,_jp_maxpu);
    }
    string aname = Form("jt500_%d",ID);
    string ename = Form("aftjt500_%d",ID);
    hltTrigs[aname] = new TH1D(aname.c_str(),aname.c_str(),_jp_maxpu,0,_jp_maxpu);
    aftTrigs[ename] = new TH1D(ename.c_str(),ename.c_str(),_jp_maxpu,0,_jp_maxpu);
  }
  map<string,TH1D*> versionTrigs;
  for (auto i = 0u; i < _jp_notrigs; ++i) {
    string iname = Form("vers%s",_jp_triggers[i]);
    versionTrigs[iname] = new TH1D(iname.c_str(),iname.c_str(),25,0,25);
  }
  versionTrigs["versjt500"] = new TH1D("versjt500","versjt500",25,0,25);

  for (Long64_t djentry=0; djentry<nentries;djentry+=1+_jp_skim) { // Event loop
    Long64_t jentry = djentry+nskip;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%50000==0) cout << "." << flush;

    auto infoplace = std::find(infojentrys.begin(),infojentrys.end(),djentry);
    if (infoplace!=infojentrys.end()) {
      // How long is it going to take?
      cout << endl << Form("Processed %ld events (%1.1f%%) in %1.0f sec.",
                      (long int)jentry, 100.*djentry/nentries, stop.RealTime()) << endl;
      cout << "BGN: ";
      bgn.Print();
      TDatime now;
      cout << "NOW: ";
      now.Print();
      if (++infoplace!=infojentrys.end()) {
        TDatime nxt;
        nxt.Set(nxt.Convert()+static_cast<UInt_t>(stop.RealTime()*(static_cast<Double_t>((*infoplace))/static_cast<Double_t>(djentry)-1.0)));
        cout << "NXT: ";
        nxt.Print();
      }
      now.Set(now.Convert()+static_cast<UInt_t>(stop.RealTime()*(static_cast<Double_t>(nentries)/static_cast<Double_t>(djentry)-1.0)));
      cout << "ETA: ";
      now.Print();

      stop.Continue();
    }

    // Set auxiliary event variables NOTE: These are kept here, so that the EvtHdr stuff etc.
    assert(_jp_isdt or !_jp_pthatbins or _pthatweight>0);
    weight = EvtHdr__mWeight;
    // Weight is most likely a redundant constant if pthatbins are used, but keep it here just in case.
    if (_jp_ismc and _jp_pthatbins) weight *= _pthatweight;
    // REMOVED: "TEMP PATCH"
    run = EvtHdr__mRun;
    evt = EvtHdr__mEvent;
    lbn = EvtHdr__mLumi;
    trpu = EvtHdr__mTrPu;

    if (_jp_isdt) trpu = _avgpu[run][lbn];

    met = PFMet__et_;
    metphi = PFMet__phi_;
    metsumet = PFMet__sumEt_;

    ++_cnt["01all"];

    // Check if good run/LS, including JSON selection
    if (_jp_isdt) {
      if (_jp_dojson) {
        // Does the run/LS pass the latest JSON selection?
        if (_json[run][lbn]==0)
          continue;
      } // _jp_dojson

      if (_jp_dolumi) {
        // Do we have the run listed in the .csv file?
        auto irun = _lums.find(run);
        if (irun==_lums.end())
          continue;
        // Do we have the LS listed in the .csv file?
        auto ils = irun->second.find(lbn);
        if (ils==irun->second.end())
          continue;
      } // _jp_dolumi
    } // _jp_isdt

    ++_cnt["02ls"];

    // Reset event ID
    _pass = true;

    // Reset prescales (dynamic can change within run)
    for (auto &scaleit : _prescales)
      scaleit.second[run] = 0;

    // Fill trigger information
    _trigs.clear();

    // Simulate other triggers for MC, if so wished
    // (this is slow, though)
    if (_jp_ismc) {
      // Always insert the generic mc trigger
      if (_jp_domctrigsim and njt>0) {
        // Only add the greatest trigger present
        for (int itrg = _jp_notrigs; itrg > 0; --itrg) {
          if (jtpt[0]>_jp_trigranges[itrg-1][0]) {
            _trigs.insert(_jp_triggers[itrg-1]);
            break; // Don't add lesser triggers
          }
        }
      } // _jp_domctrigsim
    } // _jp_ismc

    if (_jp_isdt) {
      // For data, check trigger bits
      if (_jp_debug) {
        cout << "TriggerDecision_.size()=="<<TriggerDecision_.size()<<endl<<flush;
        cout << "_availTrigs.size()=="<<_availTrigs.size()<<endl<<flush;
        cout << "_goodTrigs.size()=="<<_goodTrigs.size()<<endl<<flush;
      }
      assert(TriggerDecision_.size() == _availTrigs.size());

      for (auto itrg = 0u; itrg != _jp_notrigs; ++itrg)
        _wt[string(_jp_triggers[itrg])] = 1.0;

      for (auto goodIdx = 0u; goodIdx < _goodTrigs.size(); ++goodIdx) {
        auto &itrg = _goodTrigs[goodIdx];
        if (TriggerDecision_[itrg]!=1) continue; // -1, 0, 1
        double iwgt = 0.0;
        if (_jp_useversionlumi) iwgt = _goodWgts[goodIdx];

        string trigname = _availTrigs[itrg];
        if (_jp_debug and TriggerDecision_[itrg]>0)
          cout << trigname << " " << TriggerDecision_[itrg]
               << " " << L1Prescale_[itrg] << " " << HLTPrescale_[itrg] << endl;

        // Set prescale from event for now
        if (L1Prescale_[itrg]>0 and HLTPrescale_[itrg]>0) {
          _prescales[trigname][run] = L1Prescale_[itrg] * HLTPrescale_[itrg];
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

        if (_prescales[trigname][run]!=0) {
          // Set trigger only if prescale information is known
          string newname = Form("%s_%d",trigname.c_str(),_goodNos[goodIdx]);
          string versname = Form("vers%s",trigname.c_str());
          _trigs.insert(newname.c_str());
          if (_jp_useversionlumi) _wt[trigname] = iwgt;
          versionTrigs[versname]->Fill(_goodVNos[goodIdx],_w0 * _wt[trigname]);
          hltTrigs[newname]->Fill(trpu,_w0 * _wt[trigname]);
        } else {
          // Make sure all info is good! This is crucial if there is something odd with the tuples
          *ferr << "Missing prescale for " << trigname
                << " in run " << run << endl << flush;
        }
      } // for itrg (_goodTrigs)
    }

    _pass = _pass and _trigs.size()>0;

    // Retrieve event weight
    _w0 = (_jp_ismc ? weight : 1);
    assert(_w0>0);
    _w = _w0;

    // Calculate trigger PU weight
    if (_jp_ismc) {
      for (auto itrg = 0u; itrg != _jp_notrigs; ++itrg) {
        const char *trg_name = string(_jp_triggers[itrg]).c_str();
        _wt[trg_name] = 1.;

        // Reweight in-time pile-up
        if (_jp_reweighPU) {
          int k = _pudist[trg_name]->FindBin(trpu);
          double wtrue = _pudist[trg_name]->GetBinContent(k);
          _wt[trg_name] *= wtrue;

          // check for non-zero PU weight
          _pass = _pass and wtrue!=0;
        }
      } // for itrg
      if (_pass) ++_cnt["07puw"];
      else continue; // Fatal failure
    }

    // Equipped in fillBasics and fillRunHistos
    _pass_qcdmet = (met < 0.4 * metsumet || met < 45.); // QCD-11-004

    if (_pass and _pass_qcdmet) {
      regex dummyrgx("(jt[0-9]*)_([0-9]*)");
      for (auto &trg : _trigs) {
        string trgtrue = std::regex_replace(trg, dummyrgx, "$1", std::regex_constants::format_no_copy);
        string aftname = Form("aft%s",trg.c_str());
        aftTrigs[aftname]->Fill(trpu,_w0 * _wt[trgtrue]);
      }
    }
  } // for jentry
  cout << endl;

  stop.Stop();
  TDatime now;
  cout << "Stopping at: ";
  now.Print();

  for (auto &hist : hltTrigs) hist.second->Write();
  for (auto &hist : versionTrigs) hist.second->Write();
  for (auto &hist : aftTrigs) hist.second->Write();

  delete ferr;
}

// Load good run and LS information
bool histosInfo::loadJSON(const char* filename)
{
  cout << "Processing loadJSON(\"" << filename << "\"..." << endl;
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

  cout << "Called loadJSON(\"" << filename << "\"):" << endl;
  cout << "Loaded " << nrun << " good runs and " << nls
       << " good lumi sections" << endl;
  return true;
} // loadJSON


// Load luminosity information
bool histosInfo::loadLumi(const char* filename)
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
  if (s!="#Data tag : v1 , Norm tag: None") return false;

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
  return true;
} // loadLumi


bool histosInfo::loadPUProfiles(const char *datafile, const char *mcfile)
{
  cout << "Processing loadPUProfiles(\"" << datafile << "\",\"" << mcfile << "\")..." << endl;

  TDirectory *curdir = gDirectory;

  // Load pile-up files and hists from them
  TFile *f_pudist = new TFile(datafile, "READ");
  if (!f_pudist or f_pudist->IsZombie()) return false;
  TFile *fpumc = new TFile(mcfile,"READ");
  if (!fpumc or fpumc->IsZombie()) return false;

  _pumc = dynamic_cast<TH1D*>(fpumc->Get("pileupmc"));
  if (!_pumc) return false;
  _pumc->Scale(1./_pumc->Integral());
  double maxmcpu = _pumc->GetMaximum();
  int lomclim = _pumc->FindFirstBinAbove(maxmcpu/100.0);
  int upmclim = _pumc->FindLastBinAbove(maxmcpu/100.0);
  for (int bin = 0; bin < lomclim; ++bin)
    _pumc->SetBinContent(bin,0.0);
  for (int bin = upmclim+1; bin <= _pumc->GetNbinsX(); ++bin)
    _pumc->SetBinContent(bin,0.0);
  *ferr << "Discarding mc pu below & above: " << _pumc->GetBinLowEdge(lomclim) << " " << _pumc->GetBinLowEdge(upmclim+1) << endl;
  // Normalize
  int nbinsmc = _pumc->GetNbinsX();
  int kmc = _pumc->FindBin(33);

  // For data, load each trigger separately
  for (auto itrg = 0u ; itrg != _jp_notrigs; ++itrg) {
    string t = string(_jp_triggers[itrg]);
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
    _pudist[t]->Scale(1./_pudist[t]->Integral());
    double maxdtpu = _pudist[t]->GetMaximum();
    int lodtlim = _pudist[t]->FindFirstBinAbove(maxdtpu/100.0);
    int updtlim = _pudist[t]->FindLastBinAbove(maxdtpu/100.0);
    for (int bin = 0; bin < lodtlim; ++bin)
      _pudist[t]->SetBinContent(bin,0.0);
    for (int bin = updtlim+1; bin <= _pudist[t]->GetNbinsX(); ++bin)
      _pudist[t]->SetBinContent(bin,0.0);
    *ferr << "Discarding dt pu below & above: " << _pudist[t]->GetBinLowEdge(lodtlim) << " " << _pudist[t]->GetBinLowEdge(updtlim+1) << " " << t << endl;
    _pudist[t]->Divide(_pumc);
  }
  // REMOVED: "data with only one histo:"

  curdir->cd();
  return true;
} // loadPUProfiles



// Update the available trigger types for each new tree
bool histosInfo::getTriggers()
{
  TH1F *triggers = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerNames")); assert(triggers);
  TAxis *xax = triggers->GetXaxis();

  regex pfjet("HLT_PFJet([0-9]*)_v([0-9]*)");
  regex ak8("HLT_AK8PFJet([0-9]*)_v[0-9]*");
  regex jetht("HLT_PFHT([0-9]*)_v[0-9]*");

  _availTrigs.clear();
  _goodTrigs.clear();
  _goodWgts.clear();
  _goodNos.clear();
  _goodVNos.clear();
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
        int versdummy = stoi(std::regex_replace(trgName, pfjet, "$2", std::regex_constants::format_no_copy));
        bool found = false;
        for (auto IDidx = 0u; IDidx < _jp_notrigIDs; ++IDidx) {
          for (auto &currTag : _jp_trigtags[IDidx]) {
            if (std::regex_match(trgdummy,currTag)) {
              if (found) {
                cerr << "Double match for " << trgName << ", check trigger naming in _jp_trigtags. Aborting..." << endl;
                return false;
              }
              found = true;
              if (thrplace < _jp_notrigs+1) {
                _goodWgts.push_back(_jp_trigwgts[IDidx][thrplace]);
                _goodNos.push_back(IDidx);
                _goodVNos.push_back(versdummy);
              } else {
                cerr << "Error searching the trigger weight! Aborting..." << endl;
                return false;
              }
            }
          }
        }
        if (!found) {
          cout << "No info for " << trgName << " in _jp_trigtag. Could be a dummy trigger in the tuple." << endl;
          _goodWgts.push_back(0.0);
          _goodNos.push_back(-1.0);
          _goodVNos.push_back(-1.0);
        }
      } else {
        _goodWgts.push_back(_jp_triglumi[_jp_notrigs-1]/_jp_triglumi[thrplace]);
        _goodNos.push_back(-1.0);
        _goodVNos.push_back(-1.0);
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
} // getTriggers
