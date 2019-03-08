// Purpose: Combine different triggers into a single spectrum
// Authors:  mikko.voutilainen@cern.ch
//           hannu.siikonen@cern.ch
// Created: March 22, 2010
// Updated: November 21, 2018

// We throw a name to recurseFile and the function loops over all directories.
// As we reach a certain level, a new histogram is declared.
// Onto this histogram we add all the weighted products of the jtX subdirectories.
// The output directory is used as a temporary saving place for the histos, since pointers in ROOT are problematic.

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TRegexp.h"

#include <iostream>
#include <map>
#include <string>
#include <set>

#include "settings.h"

using namespace std;

void recurseFile(TDirectory *indir, TDirectory *outdir, string hname = "hpt", bool ptSelect = true,
                 bool othProf = false, int lvl = -1, double etamid = 0);
map<string, pair<double, double> > _ptranges;
vector<string> _ignoretrgs;


// global variables (not pretty, but works)
TDirectory *_top = 0;

void HistosCombine() {
  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-2a.root",jp::type),"READ");
  assert(fin and !fin->IsZombie());
  _top = gDirectory;

  TFile *fout = new TFile(Form("output-%s-2b.root",jp::type),"RECREATE");
  assert(fout and !fout->IsZombie());

  cout << "Calling HistosCombine("<<jp::type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursions. These may take a few seconds" << endl << flush;

  // Store pT ranges to a nice map
  if (jp::ismc) {
    if (jp::usemctrig) {
      _ptranges["mc"] = pair<double,double>(0., jp::sqrts/2.);
    } else {
      _ptranges["mc"] = pair<double,double>(0., 0.);
      _ignoretrgs.push_back("mc");
    }
  }
  for (auto itrg = 0u; itrg != jp::notrigs; ++itrg) {
    _ptranges[jp::triggers[itrg]] = pair<double, double>(jp::trigranges[itrg][0], jp::trigranges[itrg][1]);

    if (jp::ismc and jp::usemctrig) { // When in mc, we need to know how to use trigger weighting. Either mc or triggers are ignored
      _ptranges[jp::triggers[itrg]] = pair<double,double>(0.,0.);
      _ignoretrgs.push_back(jp::triggers[itrg]);
    }
  }

  // Loop over all the directories recursively
  // List here the histograms that need merging
  recurseFile(fin, fout, "hpt");
  recurseFile(fin, fout, "hpt_pre");
  recurseFile(fin, fout, "hpt_notrigeff");
  recurseFile(fin, fout, "htrigeff");
  recurseFile(fin, fout, "htrigeffmc");
  recurseFile(fin, fout, "htrigeffsf");
  recurseFile(fin, fout, "hpt_notimedep");
  recurseFile(fin, fout, "hpt_withtimedep");
  recurseFile(fin, fout, "htimedep");
  recurseFile(fin, fout, "htimefit");
  if (jp::isdt) recurseFile(fin, fout, "hlumi");
  if (jp::isdt) recurseFile(fin, fout, "hlumi_orig");

  recurseFile(fin, fout, "hpt_evt");
  recurseFile(fin, fout, "hpt_evtcount");
  recurseFile(fin, fout, "hpt_jet");

  if (jp::ismc) recurseFile(fin, fout, "hpt_g0tw");

  recurseFile(fin, fout, "hpt0");
  recurseFile(fin, fout, "hpt1");
  recurseFile(fin, fout, "hpt2");
  recurseFile(fin, fout, "hpt3");
  recurseFile(fin, fout, "hpt_noid");
  recurseFile(fin, fout, "hpt_noevtid");
  recurseFile(fin, fout, "hpt_nojetid");

  recurseFile(fin, fout, "pa");
  recurseFile(fin, fout, "pnpv");
  recurseFile(fin, fout, "pnpvall");
  recurseFile(fin, fout, "prho");
  //
  recurseFile(fin, fout, "pchs");
  //
  recurseFile(fin, fout, "pmpf");
  recurseFile(fin, fout, "pmpf1");
  recurseFile(fin, fout, "pmpf2");
  recurseFile(fin, fout, "pmpfx");
  recurseFile(fin, fout, "pmpfy");
  recurseFile(fin, fout, "pmpfz");
  recurseFile(fin, fout, "punc");
  //
  recurseFile(fin, fout, "pncand");
  recurseFile(fin, fout, "pnch");
  recurseFile(fin, fout, "pnne");
  recurseFile(fin, fout, "pnnh");
  recurseFile(fin, fout, "pnce");
  recurseFile(fin, fout, "pnmu");

  recurseFile(fin, fout, "pchf");
  recurseFile(fin, fout, "pnef");
  recurseFile(fin, fout, "pnhf");
  recurseFile(fin, fout, "pcef");
  recurseFile(fin, fout, "pmuf");
  recurseFile(fin, fout, "phhf");
  recurseFile(fin, fout, "phef");
  recurseFile(fin, fout, "pbeta");
  recurseFile(fin, fout, "pbetastar");
  recurseFile(fin, fout, "ppuf");
  recurseFile(fin, fout, "pjec");
  recurseFile(fin, fout, "pjec2");
  recurseFile(fin, fout, "pjec_res");
  //
  recurseFile(fin, fout, "pncandtp");
  recurseFile(fin, fout, "pnchtp");
  recurseFile(fin, fout, "pnnetp");
  recurseFile(fin, fout, "pnnhtp");
  recurseFile(fin, fout, "pncetp");
  recurseFile(fin, fout, "pnmutp");

  recurseFile(fin, fout, "pchftp");
  recurseFile(fin, fout, "pneftp");
  recurseFile(fin, fout, "pnhftp");
  recurseFile(fin, fout, "pceftp");
  recurseFile(fin, fout, "pmuftp");
  recurseFile(fin, fout, "phhftp");
  recurseFile(fin, fout, "pheftp");
  recurseFile(fin, fout, "pbetatp");
  recurseFile(fin, fout, "pbetastartp");
  recurseFile(fin, fout, "ppuftp");
  recurseFile(fin, fout, "pjectp");
  //
  recurseFile(fin, fout, "pchftp_vsnpv", false, true);
  recurseFile(fin, fout, "pneftp_vsnpv", false, true);
  recurseFile(fin, fout, "pnhftp_vsnpv", false, true);
  recurseFile(fin, fout, "pceftp_vsnpv", false, true);
  recurseFile(fin, fout, "pmuftp_vsnpv", false, true);
  recurseFile(fin, fout, "phhftp_vsnpv", false, true);
  recurseFile(fin, fout, "pheftp_vsnpv", false, true);
  recurseFile(fin, fout, "pbetatp_vsnpv", false, true);
  recurseFile(fin, fout, "pbetastartp_vsnpv", false, true);
  recurseFile(fin, fout, "ppuftp_vsnpv", false, true);
  //
  recurseFile(fin, fout, "pchftp_vstrpu", false, true);
  recurseFile(fin, fout, "pneftp_vstrpu", false, true);
  recurseFile(fin, fout, "pnhftp_vstrpu", false, true);
  recurseFile(fin, fout, "pceftp_vstrpu", false, true);
  recurseFile(fin, fout, "pmuftp_vstrpu", false, true);
  recurseFile(fin, fout, "phhftp_vstrpu", false, true);
  recurseFile(fin, fout, "pheftp_vstrpu", false, true);
  recurseFile(fin, fout, "pbetatp_vstrpu", false, true);
  recurseFile(fin, fout, "pbetastartp_vstrpu", false, true);
  recurseFile(fin, fout, "ppuftp_vstrpu", false, true);
  //
  recurseFile(fin, fout, "pchftp_vseta", false, true);
  recurseFile(fin, fout, "pneftp_vseta", false, true);
  recurseFile(fin, fout, "pnhftp_vseta", false, true);
  recurseFile(fin, fout, "pceftp_vseta", false, true);
  recurseFile(fin, fout, "pmuftp_vseta", false, true);
  recurseFile(fin, fout, "phhftp_vseta", false, true);
  recurseFile(fin, fout, "pheftp_vseta", false, true);
  recurseFile(fin, fout, "pbetatp_vseta", false, true);
  recurseFile(fin, fout, "pbetastartp_vseta", false, true);
  recurseFile(fin, fout, "ppuftp_vseta", false, true);
  //
  recurseFile(fin, fout, "pchfnegtp_vsphi", false, true);
  recurseFile(fin, fout, "pnefnegtp_vsphi", false, true);
  recurseFile(fin, fout, "pnhfnegtp_vsphi", false, true);
  recurseFile(fin, fout, "pcefnegtp_vsphi", false, true);
  recurseFile(fin, fout, "pmufnegtp_vsphi", false, true);
  recurseFile(fin, fout, "phhfnegtp_vsphi", false, true);
  recurseFile(fin, fout, "phefnegtp_vsphi", false, true);
  recurseFile(fin, fout, "pbetanegtp_vsphi", false, true);
  recurseFile(fin, fout, "pbetastarnegtp_vsphi", false, true);
  recurseFile(fin, fout, "ppufnegtp_vsphi", false, true);
  //
  recurseFile(fin, fout, "pchfpostp_vsphi", false, true);
  recurseFile(fin, fout, "pnefpostp_vsphi", false, true);
  recurseFile(fin, fout, "pnhfpostp_vsphi", false, true);
  recurseFile(fin, fout, "pcefpostp_vsphi", false, true);
  recurseFile(fin, fout, "pmufpostp_vsphi", false, true);
  recurseFile(fin, fout, "phhfpostp_vsphi", false, true);
  recurseFile(fin, fout, "phefpostp_vsphi", false, true);
  recurseFile(fin, fout, "pbetapostp_vsphi", false, true);
  recurseFile(fin, fout, "pbetastarpostp_vsphi", false, true);
  recurseFile(fin, fout, "ppufpostp_vsphi", false, true);

  recurseFile(fin, fout, "ppt_probepertag");

  recurseFile(fin, fout, "hdjasymm");
  recurseFile(fin, fout, "hdjasymmtp");
  recurseFile(fin, fout, "hdjasymmpt");
  recurseFile(fin, fout, "hdjmpf");
  recurseFile(fin, fout, "hdjmpftp");
  recurseFile(fin, fout, "hdjmpfpt");

  recurseFile(fin, fout, "hdjasymm_a005");
  recurseFile(fin, fout, "hdjasymmtp_a005");
  recurseFile(fin, fout, "hdjasymmpt_a005");
  recurseFile(fin, fout, "hdjmpf_a005");
  recurseFile(fin, fout, "hdjmpftp_a005");
  recurseFile(fin, fout, "hdjmpfpt_a005");
  recurseFile(fin, fout, "hdjasymm_a01");
  recurseFile(fin, fout, "hdjasymmtp_a01");
  recurseFile(fin, fout, "hdjasymmpt_a01");
  recurseFile(fin, fout, "hdjmpf_a01");
  recurseFile(fin, fout, "hdjmpftp_a01");
  recurseFile(fin, fout, "hdjmpfpt_a01");
  recurseFile(fin, fout, "hdjasymm_a015");
  recurseFile(fin, fout, "hdjasymmtp_a015");
  recurseFile(fin, fout, "hdjasymmpt_a015");
  recurseFile(fin, fout, "hdjmpf_a015");
  recurseFile(fin, fout, "hdjmpftp_a015");
  recurseFile(fin, fout, "hdjmpfpt_a015");
  recurseFile(fin, fout, "hdjasymm_a02");
  recurseFile(fin, fout, "hdjasymmtp_a02");
  recurseFile(fin, fout, "hdjasymmpt_a02");
  recurseFile(fin, fout, "hdjmpf_a02");
  recurseFile(fin, fout, "hdjmpftp_a02");
  recurseFile(fin, fout, "hdjmpfpt_a02");
  recurseFile(fin, fout, "hdjasymm_a025");
  recurseFile(fin, fout, "hdjasymmtp_a025");
  recurseFile(fin, fout, "hdjasymmpt_a025");
  recurseFile(fin, fout, "hdjmpf_a025");
  recurseFile(fin, fout, "hdjmpftp_a025");
  recurseFile(fin, fout, "hdjmpfpt_a025");
  recurseFile(fin, fout, "hdjasymm_a03");
  recurseFile(fin, fout, "hdjasymmtp_a03");
  recurseFile(fin, fout, "hdjasymmpt_a03");
  recurseFile(fin, fout, "hdjmpf_a03");
  recurseFile(fin, fout, "hdjmpftp_a03");
  recurseFile(fin, fout, "hdjmpfpt_a03");

  recurseFile(fin, fout, "hdjresp_tag_a005");
  recurseFile(fin, fout, "hdjresptp_tag_a005");
  recurseFile(fin, fout, "hdjresp_tag_a01");
  recurseFile(fin, fout, "hdjresptp_tag_a01");
  recurseFile(fin, fout, "hdjresp_tag_a015");
  recurseFile(fin, fout, "hdjresptp_tag_a015");
  recurseFile(fin, fout, "hdjresp_tag_a02");
  recurseFile(fin, fout, "hdjresptp_tag_a02");
  recurseFile(fin, fout, "hdjresp_tag_a025");
  recurseFile(fin, fout, "hdjresptp_tag_a025");
  recurseFile(fin, fout, "hdjresp_tag_a03");
  recurseFile(fin, fout, "hdjresptp_tag_a03");

  recurseFile(fin, fout, "hdjresp_probe_a005");
  recurseFile(fin, fout, "hdjresptp_probe_a005");
  recurseFile(fin, fout, "hdjresp_probe_a01");
  recurseFile(fin, fout, "hdjresptp_probe_a01");
  recurseFile(fin, fout, "hdjresp_probe_a015");
  recurseFile(fin, fout, "hdjresptp_probe_a015");
  recurseFile(fin, fout, "hdjresp_probe_a02");
  recurseFile(fin, fout, "hdjresptp_probe_a02");
  recurseFile(fin, fout, "hdjresp_probe_a025");
  recurseFile(fin, fout, "hdjresptp_probe_a025");
  recurseFile(fin, fout, "hdjresp_probe_a03");
  recurseFile(fin, fout, "hdjresptp_tag_a03");

  recurseFile(fin, fout, "hdjmass", false);
  recurseFile(fin, fout, "hdjmass0", false);

  if (!jp::ismc) recurseFile(fin, fout, "peff_new");

  curdir->cd();

  cout << endl;
  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();
  cout << "Output file closed" << endl;

  fin->Close();
  fin->Delete();
  cout << "Input file closed" << endl;
} // HistosNormalize

inline void binCp(TH1 *_hpt, TH1 *hpt, int bin) { // Copying a TH1 bin is simple
  _hpt->SetBinContent(bin, hpt->GetBinContent(bin));
  _hpt->SetBinError(bin, hpt->GetBinError(bin));
}

inline void binCp(TProfile *_hpt, TProfile *hpt, int bin) { // Copying a TProfile bin is complicated
  (*_hpt)[bin] += (*hpt)[bin]; // GetW
  (*_hpt->GetSumw2())[bin] += (*hpt->GetSumw2())[bin]; // GetW2
  _hpt->SetBinEntries(bin, hpt->GetBinEntries(bin)+_hpt->GetBinEntries(bin)); // GetB
  (*_hpt->GetBinSumw2())[bin] = (*hpt->GetBinSumw2())[bin]; // GetB2
}

template<typename T>
inline void fillHisto(T *hpt, TDirectory *outdir, TDirectory *indir, bool ptSelect, bool othProf) {
  const char* hname = hpt->GetName();
  T *_hpt = dynamic_cast<T*>(dynamic_cast<TObject*>(outdir->Get(hpt->GetName())));
  if (!_hpt) {
    cout << "Savings histogram fetch failed: " << hpt->GetName() << " in " << indir->GetName() << endl;
    if (jp::debug) cout << "Cloned _" << hname << endl;
    return;
  }

  if (ptSelect) {
    if (_ptranges.find(indir->GetName())==_ptranges.end())
      cout << "pT range not found for directory " << indir->GetName() << endl;
    assert(_ptranges.find(indir->GetName())!=_ptranges.end());
    double ptmin = _ptranges[indir->GetName()].first;
    double ptmax = _ptranges[indir->GetName()].second;

    bool isth3 = hpt->InheritsFrom("TH3");
    bool isth23 = isth3 or hpt->InheritsFrom("TH2");
    for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
      double pt = hpt->GetXaxis()->GetBinCenter(i);
      if (pt > ptmin and pt < ptmax) { // TODO: We could do better than a linear search!
        if (isth23) { // TH2 and TH3
          for (int j = 1; j != _hpt->GetNbinsY()+1; ++j) {
            if (isth3) {
              for (int k = 1; k != _hpt->GetNbinsZ()+1; ++k) {
                _hpt->SetBinContent(i,j,k, hpt->GetBinContent(i,j,k));
                _hpt->SetBinError(i,j,k, hpt->GetBinError(i,j,k));
              } // for k
            } else {
              _hpt->SetBinContent(i,j, hpt->GetBinContent(i,j));
              _hpt->SetBinError(i,j, hpt->GetBinError(i,j));
            }
          } // for j
        } else { // TH1 and TProfile
          binCp(_hpt,hpt,i);
        }
      } // in ptrange
    } // for i
  } else {
    _hpt->Add(hpt);
    outdir->cd();
  }
}


void recurseFile(TDirectory *indir, TDirectory *outdir, string hname, bool ptSelect, bool othProf, int lvl, double etamid)
{
  TDirectory *curdir = gDirectory;
  // Automatically go through the list of keys (directories)
  TObject *obj = 0;
  if (lvl==-1) {
    lvl++;
    cout << endl << hname << " ";
  }

  if (lvl<2) {
    TList *keys = indir->GetListOfKeys();
    TKey *key = 0;
    TListIter itkey(keys);
    while ((key = dynamic_cast<TKey*>(itkey.Next()))) {
      // We want the key to match to a directory or the selected histogram
      string classname = key->GetClassName();
      string kname = key->GetName();
      if (classname=="TDirectoryFile") { // All keys are looped over only in search for subdirs
        obj = key->ReadObj(); assert(obj);
      } else {
        continue;
      }
      if (obj->InheritsFrom("TDirectory")) {
        int loclvl = lvl;
        // Found a subdirectory: copy it to output and go deeper

        TDirectory *outdir2 = outdir;
        if (loclvl == 0) {
          if (outdir->FindKey(obj->GetName())==0) outdir->mkdir(obj->GetName());
          bool enterkeydir = outdir->cd(kname.c_str());
          assert(enterkeydir);
          outdir2 = outdir->GetDirectory(kname.c_str()); assert(outdir2);
          outdir2->cd();
          if (jp::debug) cout << kname.c_str() << endl;
        } else {
          if (jp::debug) cout << kname.c_str() << " (at bottom)" << endl;
          loclvl++; // Increase the level by one when we enter e.g. jt40
        }

        bool enterobjdir = indir->cd(obj->GetName());
        assert(enterobjdir);
        TDirectory *indir2 = indir->GetDirectory(obj->GetName()); assert(indir2);
        indir2->cd();

        // Check if directory name contains information on eta bin width
        // If yes, the next level is the bottom level with triggers
        // set flag and reset the combined histogram pointer
        float etamin, etamax;
        if (loclvl==0 and (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2) and (etamax>etamin) ) {
          outdir2->cd();
          TObject *inobj = indir2->Get(Form("%s/%s",jp::isdt ? jp::reftrig : "mc",hname.c_str()));
          if (inobj) {
            cout << "[" << etamin << "," << etamax << "]";
            TH1 *hpt = dynamic_cast<TH1*>(inobj->Clone(hname.c_str()));
            hpt->Reset();
            if (hpt) {
              recurseFile(indir2, outdir2, hname, ptSelect, othProf, 1, 0.5*(etamin+etamax));
              hpt->Write();
              hpt->Delete();
              indir2->cd();
            }
          }
        } else if (loclvl==0 and TString(indir2->GetName()).Contains("FullEta")) {
          outdir2->cd();
          TObject *inobj = indir2->Get(Form("%s/%s",jp::isdt ? jp::reftrig : "mc",hname.c_str()));
          if (inobj) {
            if (string(indir2->GetName())!="FullEta_Gen") cout << "FullEta";
            TH1 *hpt = dynamic_cast<TH1*>(inobj->Clone(hname.c_str()));
            hpt->Reset();
            if (hpt) {
              recurseFile(indir2, outdir2, hname, ptSelect, othProf, 1, etamid);
              hpt->Write();
              hpt->Delete();
              indir2->cd();
            }
          }
        } else {
          if (loclvl==0 or std::find(_ignoretrgs.begin(),_ignoretrgs.end(),indir2->GetName())==_ignoretrgs.end()) {
            string sid = indir->GetName();
            if (loclvl>0 and (sid=="Eta_0.0-1.3" or (TString(sid).Contains("FullEta") and sid!="FullEta_Gen"))) {
              if (string(indir2->GetName())=="jt40" or string(indir2->GetName())=="mc") cout << " ";
              cout << indir2->GetName();
              if (string(indir2->GetName())!=jp::reftrig) cout << ",";
            }
            recurseFile(indir2, outdir2, hname, ptSelect, othProf, loclvl, etamid);
          }
        }
      }
      obj->Delete();
    } // while key
  } else {
    obj = indir->Get(hname.c_str());
    if (obj) {
      // Copy over TProfile, TH3, TH2, TH1 histograms, in this precise order.
      if (obj->InheritsFrom("TProfile")) {
        TProfile *hpt = dynamic_cast<TProfile*>(obj);
        if (hpt) fillHisto(hpt, outdir, indir, ptSelect, othProf);
      } else {
        TH1 *hpt = 0;
        if (obj->InheritsFrom("TH3")) hpt = dynamic_cast<TH3D*>(obj);
        else if (obj->InheritsFrom("TH2")) hpt = dynamic_cast<TH2D*>(obj);
        else if (obj->InheritsFrom("TH1")) hpt = dynamic_cast<TH1D*>(obj);
        if (hpt) fillHisto(hpt, outdir, indir, ptSelect, othProf);
      }
      obj->Delete();
    }
  }
  curdir->cd();
} // recurseFile
