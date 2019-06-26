// A handy script for fetching the number of entries in a pthat sliced MC sample
// The user should replace fileformat and dirname according to present needs
// Aouthor Hannu Siikonen hannu.siikonen@cern.ch

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

void GetSliceEvts() {
  //TString dirname="/afs/cern.ch/user/h/hsiikone/eosw/cms/store/group/phys_smp/Multijet/13TeV/MC/P825ns80X_Moriond17";
  //std::regex fileformat("QCD_Pt_([0-9]*)to([0-9]*|Inf)_TuneCUETP8M_13TeV_pythia8.root");
  TString dirname="/work/jet_tuples/MC/2016/P8M1_Moriond17";
  std::regex fileformat("Pthat_([0-9]*)to([0-9]*|Inf).root");
  std::cmatch match;

  TSystemDirectory dir(dirname.Data(), dirname.Data());
  TList *files = dir.GetListOfFiles();
  set<pair<unsigned int,pair<unsigned int,TString> > > vals;
  if (files) {
    cout << "Found some files!" << endl;
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      cout << "Processing file " << fname << endl;
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        cout << "File is a .root file!" << endl;
        if (std::regex_match(fname.Data(), match, fileformat)) {
          cout << "File was matched with regex!" << endl;
          TFile *f = new TFile((dirname+"/"+fname).Data());
          TTree *t = (TTree*) f->Get("ak4/ProcessedTree");
          vals.insert(std::make_pair(stoi(match[1]),std::make_pair(t->GetEntries(),fname)));
        }
      }
    }
  }

  for (auto it = vals.begin(); it != vals.end(); ++it) {
    cout << "Lower pthat: " << it->first << " Entries: " << it->second.first << endl;
  }
  cout << "That is: " << endl;
  cout << "{";
  for (auto it = vals.begin(); it != vals.end(); ++it) {
    cout << it->second.first << ",";
  }
  cout << "}" << endl;
  for (auto it = vals.begin(); it != vals.end(); ++it) {
    cout << it->second.second << endl;
  }
}
