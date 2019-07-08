#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"

void normalize() {
  vector<string> tags = {"500", "450", "400", "320", "260", "200", "140", "80", "60", "40", "0"};
//  vector<string> tags = {"450", "400", "320", "260", "200", "140", "80", "60", "40", "0"};
  TFile f("pileup_DT.root","RECREATE");
  TFile g("pileup_DTu.root","READ");

  double scale = 1.0;
  for (unsigned i = 0; i < tags.size(); ++i) {
    TString name;
    name.Form("jt%s", tags[i].c_str());
    g.cd();
    TH1D* dhandle = dynamic_cast<TH1D*>(gROOT->FindObject(name.Data()));
    f.cd();
    TH1D* handle = dynamic_cast<TH1D*>(dhandle->Clone(name.Data()));
    if (handle) {
      if (tags[i]=="500") scale = 1.0/handle->GetMaximum(); // The reference tag needs to be first so that the weights line up
//      if (tags[i]=="450") scale = 1.0/handle->GetMaximum(); // The reference tag needs to be first so that the weights line up
      handle->Scale(scale);
      handle->SetTitle(name.Data());
      handle->Write();
    }
  }
  g.Close();
}
