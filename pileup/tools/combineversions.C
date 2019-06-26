#include "../../../../settings.h"
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

void combineversions() {
//  vector<string> tags = {"500", "450", "400", "320", "260", "200", "140", "80", "60", "40"};
  vector<string> tags = {"450", "400", "320", "260", "200", "140", "80", "60", "40"};

  TString dirname="./";
  std::cmatch match;
  TFile f("pileup_DTu.root","RECREATE");

  TSystemDirectory dir(dirname.Data(), dirname.Data());
  TList *files = dir.GetListOfFiles();

  map<string,TH1D*> collection;
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() and fname.EndsWith(".root")) {
        std::regex fileformat("HLT_PFJet([0-9]*)_v([0-9]*).root");
        if (std::regex_match(fname.Data(), match, fileformat)) {
          string tnumber = match[1];
          string vnumber = match[2];
          TString name;
          name.Form("jt%s", tnumber.c_str());
          unsigned int thrplace = std::find(jp::trigthr,jp::trigthr+jp::notrigs,stoi(tnumber))-jp::trigthr;

          double wgt = 1.0;
          string trgdummy = Form("%s_%s",tnumber.c_str(),vnumber.c_str());
          bool found = false;
          for (auto IDidx = 0u; IDidx < jp::notrigIDs; ++IDidx) {
            for (auto &currTag : jp::trigtags[IDidx]) {
              if (std::regex_match(trgdummy,currTag)) {
                if (found) {
                  cout << "Double match for " << trgdummy << ", check trigger naming in jp::trigtags. Aborting..." << endl;
                  continue;
                }
                found = true;
                if (thrplace < jp::notrigs+1)
                  wgt = jp::trigwgts[IDidx][thrplace];
                else {
                  cout << "Error searching the trigger weight! Aborting..." << endl;
                  continue;
                }
              }
            }
          }
          if (!found) {
            cout << name.Data() << " v" << vnumber << " not found!" << endl;
            continue;
          }
          TFile *roothandle = new TFile(fname.Data(),"READ");
          roothandle->cd();
          TH1D* dhandle = dynamic_cast<TH1D*>(gROOT->FindObject("pileup"));
          f.cd();
          TH1D* handle = dynamic_cast<TH1D*>(dhandle->Clone(name.Data()));
          if (handle) {
            handle->Scale(wgt);
            handle->SetName(name.Data());
            handle->SetTitle(name.Data());
            auto place = collection.find(name.Data());
            if (place==collection.end())
              collection[name.Data()] = handle;
            else
              place->second->Add(handle,1.0);
            roothandle->Close();
          }
        }
      }
    }
  }

  f.cd();
  for (unsigned i = 0; i < tags.size(); ++i) {
    TString name;
    name.Form("jt%s", tags[i].c_str());
    collection[name.Data()]->Write();
  }
}
