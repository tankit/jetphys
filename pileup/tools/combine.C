#include "../../../../settings.h"

void combine() {
//  vector<string> tags = {"500", "450", "400", "320", "260", "200", "140", "80", "60", "40"};
  vector<string> tags = {"450", "400", "320", "260", "200", "140", "80", "60", "40"};
  TFile f("pileup_DTu.root","RECREATE");

  for (unsigned i = 0; i < tags.size(); ++i) {
    TString file, name;
    file.Form("HLT_PFJet%s.root", tags[i].c_str());
    name.Form("jt%s", tags[i].c_str());

    TFile g(file.Data(),"READ");
    g.cd();
    TH1D* handle = dynamic_cast<TH1D*>(gROOT->FindObject("pileup"));
    if (handle) {
      handle->Scale(_jp_triglumi[tags.size()-1]/_jp_triglumi[tags.size()-1-i]);
      f.cd();
      handle->SetName(name.Data());
      handle->SetTitle(name.Data());
      handle->Write();
      handle->Delete();
      g.Close();
    }
  }
}
