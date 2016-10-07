
{
    vector<const char*> tags = {"40", "60", "80", "140", "200", "260", "320", "400", "450"};
    TFile f("pileup_DT.root","RECREATE");

    for (unsigned i = 0; i < tags.size(); ++i) {
        TString file, name;
        file.Form("HLT_PFJet%s.root", tags[i]);
        name.Form("jt%s", tags[i]);

        TFile g(file.Data(),"READ");
        g.cd();
        TH1D* handle = (TH1D*)gROOT->FindObject("pileup");
        f.cd();
        handle->SetName(name.Data());
        handle->SetTitle(name.Data());
        handle->Write();
        handle->Delete();
        g.Close();
    }
}
