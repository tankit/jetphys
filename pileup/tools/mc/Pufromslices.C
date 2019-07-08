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

#include "../../../settings.h"

void Pufromslices() {
  //TString dirname="/eos/cms/store/group/phys_smp/Multijet/13TeV/MC/P825ns80X_Moriond17";
  //TString dirname="/work/jet_tuples/MC/P825ns80X_Moriond17";
  //TString dirname="/work/jet_tuples/MC/2016/P8M1";
  TString dirname="/work/jet_tuples/SIM/2016/MGP8M1/";
  //std::regex fileformat("QCD_Pt_([0-9]*)to([0-9]*|Inf)_TuneCUETP8M_13TeV_pythia8.root");
  std::regex fileformat("HT_([0-9]*)to([0-9]*|Inf).root");
  vector<int> minval = {50,100,200,300,500,700,1000,1500,2000};
  //std::regex fileformat("Pthat_([0-9]*)to([0-9]*|Inf).root");
  //vector<int> minval = {30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200};
  bool debughistos = true;

  std::cmatch match;
  TFile *output = new TFile("pileup_MC.root","RECREATE");
  TH1D *summary = new TH1D("pileupmc","",jp::maxpu,0,jp::maxpu);
  vector<TH1D*> individuals;

  TSystemDirectory dir(dirname.Data(), dirname.Data());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        if (std::regex_match(fname.Data(), match, fileformat)) {
          TFile *f = new TFile((dirname+"/"+fname).Data());
          TTree *t = (TTree*) f->Get("ak4/ProcessedTree");
          int number = stoi(match[1]);
          const char *histname = Form("pileupmc%d",number);
          TH1D *hist = new TH1D(histname,"",jp::maxpu,0,jp::maxpu);
          t->Draw(Form("EvtHdr_.mTrPu>>%s",histname));
          int pos = 0;
          while (minval[pos]!=number) ++pos;
          hist->Scale(jp::pthatsigmas[pos]/t->GetEntries());
          summary->Add(hist);
          individuals.push_back(hist);
        }
      }
    }
  }
  output->cd();
  summary->Scale(1.0/summary->GetMaximum());
  summary->Write();
  if (debughistos) {
    for (auto &hist : individuals) {
      hist->Scale(1.0/hist->GetMaximum()); // The scales are for comparison, so these cannot be used for weighting
      hist->Write();
    }
  }
}
