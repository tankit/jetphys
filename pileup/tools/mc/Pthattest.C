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

void Pthattest() {
  int bins = 100;
  int max = 6000;
  TString filename="/work/jet_tuples/MC/2017/P8CP5/Flat-long.root";
  TFile *output = new TFile("pthat.root","RECREATE");
  TFile *f = new TFile((filename).Data());
  TTree *t = (TTree*) f->Get("ak4/ProcessedTree");
  TH1D *hist1 = new TH1D("pthat","",bins,0,max);
  TH1D *hist2 = new TH1D("pthat_wgt","",bins,0,max);
  TH2D *hist2 = new TH1D("pthat_wgt","",bins,0,max);
  t->Draw("EvtHdr_.mPthat>>pthat");
  t->Draw("EvtHdr_.mPthat>>pthat_wgt","EvtHdr_.mWeight");
  t->Draw("EvtHdr_.mWeight:EvtHdr_.mPthat>>wgts");
  output->cd();
  hist1->Write();
  hist2->Write();
  wgts->Write();
  output->Close();
  output->Delete();
}
