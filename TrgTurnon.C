// Purpose: Combine different triggers into a single spectrum
// Author:  hannu.siikonen@cern.ch
// Created: May 25, 2018

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
                 int lvl = -1, double etamid = 0);

map<string,string> _friends = {{"jt40","jt0"},{"jt60","jt40"},{"jt80","jt60"},{"jt140","jt80"},{"jt200","jt140"},{"jt260","jt200"},{"jt320","jt260"},{"jt400","jt320"},{"jt450","jt400"},{"jt500","jt450"}};

// global variables (not pretty, but works)
TDirectory *_top = 0;

void TrgTurnon() {
  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-2a.root",jp::type),"READ");
  assert(fin and !fin->IsZombie());
  _top = gDirectory;

  TFile *fout = new TFile(Form("trgeff-%s.root",jp::type),"RECREATE");
  assert(fout and !fout->IsZombie());

  cout << "Calling TrgEff ("<<jp::type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursions. These may take a few seconds" << endl << flush;

  // Loop over all the directories recursively
  // List here the histograms that need merging
  recurseFile(fin, fout, "hpt");
  recurseFile(fin, fout, "hpt0");
  recurseFile(fin, fout, "hpt1");
  recurseFile(fin, fout, "hpt2");
  recurseFile(fin, fout, "hpt3");

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

template<typename T>
inline void fillHisto(T *hpt, TDirectory *outdir, TDirectory *indir, bool ptSelect) {
  const char* hname = hpt->GetName();
}


void recurseFile(TDirectory *indir, TDirectory *outdir, string hname, bool ptSelect, int lvl, double etamid)
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

        bool enterobjdir = indir->cd(obj->GetName());
        assert(enterobjdir);
        TDirectory *indir2 = indir->GetDirectory(obj->GetName()); assert(indir2);
        indir2->cd();
        if (loclvl==0 and TString(indir2->GetName()).Contains("FullEta")) return;
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
        outdir2->cd();

        float etamin, etamax;
        if (loclvl==0 and (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2) and (etamax>etamin) ) {

          TObject *inobj = indir2->Get(Form("%s/%s",jp::isdt ? jp::reftrig : "mc",hname.c_str()));
          if (inobj) {
            cout << "[" << etamin << "," << etamax << "]";
            TH1 *hpt = dynamic_cast<TH1*>(inobj->Clone(hname.c_str()));
            hpt->Reset();
            if (hpt) {
              recurseFile(indir2, outdir2, hname, ptSelect, 1, 0.5*(etamin+etamax));
              hpt->Write();
              hpt->Delete();
              indir2->cd();
            }
          }
        } else {
          string sid = indir->GetName();
          if (loclvl>0 and (sid=="Eta_0.0-1.3" or (TString(sid).Contains("FullEta") and sid!="FullEta_Gen"))) {
            if (string(indir2->GetName())=="jt40" or string(indir2->GetName())=="mc") cout << " ";
            cout << indir2->GetName();
            if (string(indir2->GetName())!=jp::reftrig) cout << ",";
          }
          recurseFile(indir2, outdir2, hname, ptSelect, loclvl, etamid);
        }
      }
      obj->Delete();
    } // while key
  } else {
    obj = indir->Get(hname.c_str());
    if (obj and _friends.find(indir->GetName()) != _friends.end()) {
      // Copy over TProfile, TH3, TH2, TH1 histograms, in this precise order.
      TH1 *hpt = 0;
      hpt = dynamic_cast<TH1D*>(obj);
      hpt->SetName(Form("%s%s",hpt->GetName(),indir->GetName()));
      string friendname = Form("../%s/%s",_friends[indir->GetName()].c_str(),hname.c_str());
      TH1 *hptd = 0;
      hptd = dynamic_cast<TH1D*>(indir->Get(friendname.c_str()));
      hpt->Divide(hptd);
      hpt->Write();

      obj->Delete();
    }
  }
  curdir->cd();
} // recurseFile
