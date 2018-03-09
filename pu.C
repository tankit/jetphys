// Purpose: plot PF composition distributions for multiple variables
// - (CHF, NHF, NEF) x (data, MC) x (direct, TP) for each pT, eta bin
// - same for multiplicities
#include "TFile.h"
#include "TH1D.h"

#include "tdrstyle_mod15.C"

void pu() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  const char *crun = "";
  TFile *fm = new TFile("./output-MC-1.root");
  assert(fm && !fm->IsZombie());

  TFile *fd = new TFile("./output-DATA-1.root");
  assert(fd && !fd->IsZombie());

  curdir->cd();
  gSystem->MakeDirectory("pdf");

  //////////////////////
  // Plotter settings //
  //////////////////////

  // Choose fractions to be plotted
  const int nfrac = 3;//4;
  string fracs[nfrac] = {"chf", "nhf", "nef"};// "cef"};

  // Define colors and markers for each fraction
  map<string,int> colors;
  colors["chf"] = kRed;
  colors["nhf"] = kGreen+2;
  colors["nef"] = kBlue;
  colors["cef"] = kCyan+2;
  map<string, int> markers;
  markers["chf"] = kFullCircle;
  markers["nhf"] = kFullDiamond;
  markers["nef"] = kFullSquare;
  markers["cef"] = kFullStar;
  map<string, int> dpmarkers;
  dpmarkers["chf"] = kOpenCircle;
  dpmarkers["nhf"] = kOpenDiamond;
  dpmarkers["nef"] = kOpenSquare;
  dpmarkers["cef"] = kOpenStar;
  map<string, const char*> namesbb;
  namesbb["chf"] = "h^{#pm}";
  namesbb["nhf"] = "h^{0}";//"K^{0}_{had}";
  namesbb["nef"] = "#gamma";//"#pi^{0}+K^{0}_{EM}";
  namesbb["cef"] = "e^{#pm} #mu^{#pm}";
  map<string, const char*> namesec;
  namesec["chf"] = "h^{#pm}_{trk}";
  namesec["nhf"] = "h_{had}";
  namesec["nef"] = "#pi^{0}+h_{EM}";
  namesec["cef"] = "e^{#pm}";

  //////////////////////////////////////////////


  // Default bin, add loop later
  const char *cd = "Standard";
  //const char *ceta = "Eta_0.0-0.5";
  //const char *ctrg = "jt450";

  const int neta = 9;
  const char* etas[neta] = {"Eta_0.0-0.5","Eta_0.5-1.0","Eta_1.0-1.5","Eta_1.5-2.0","Eta_2.0-2.5","Eta_2.5-3.0","Eta_3.0-3.2","Eta_3.2-4.7","Eta_0.0-1.3"};
  const int ntrg = 9;
  const char* trgs[ntrg] = {"jt40","jt60","jt80","jt140","jt200",
			    "jt260","jt320","jt400","jt450"};

  vector<string> filenames;
  for (int ieta = 0; ieta != neta; ++ieta) {
    const char *ceta = etas[ieta];
    float y1(0),y2(0);
    assert(sscanf(ceta,"Eta_%f-%f",&y1,&y2)==2);
    double y = 0.5*(y1+y2);
    string etafile = Form("pu_y%02d-%02d",int(10.*y1),int(10*y2));
    for (int itrg = 0; itrg != ntrg; ++itrg) {
      const char *ctrg = trgs[itrg];

      // Find out the jet pT range for this trigger
      TH1D *hselpt = (TH1D*)fm->Get(Form("%s/%s/%s/hselpt",cd,ceta,ctrg));
      assert(hselpt);
      double minpt(0), maxpt(0);
      for (int i = 2; i != hselpt->GetNbinsX()+1; ++i) {
        if (hselpt->GetBinContent(i-1)==0 && hselpt->GetBinContent(i)!=0)
          if (minpt==0) minpt = hselpt->GetBinLowEdge(i);
        if (hselpt->GetBinContent(i)!=0 && hselpt->GetBinContent(i+1)==0)
          maxpt = hselpt->GetBinLowEdge(i+1);
      }

      // For loop, delete old histograms
      TObject *tob;
      tob = gROOT->FindObject("h"); if (tob) delete tob;
      tob = gROOT->FindObject("h2"); if (tob) delete tob;
      tob = gROOT->FindObject("c1"); if (tob) delete tob;

      TLegend *leg = tdrLeg(0.70,0.89-5*0.05,0.95,0.89);

      const char *cf = "trpu";
      // Least biased tag-and-probe (TP) fractions sensitive to JEC
      TH1D *hd = (TH1D*)fd->Get(Form("%s/%s/%s/h%s",cd,ceta,ctrg,cf));
      assert(hd);
      TH1D *hm = (TH1D*)fm->Get(Form("%s/%s/%s/h%s",cd,ceta,ctrg,cf));
      assert(hm);
      TCanvas *c1 = tdrCanvas("c1",hd,4,11);
      c1->cd(1);

      // Normalize data and MC to unit area
      hd->Scale(100./hd->Integral());
      hm->Scale(100./hm->Integral());
      double maxscale = hd->GetMaximum();
      hd->Scale(1.0/maxscale);
      hm->Scale(1.0/maxscale);

      // Regular TP fractions
      int color = kRed;
      int marker = kSquare;
      tdrDraw(hm,"HIST",kNone,0,kSolid,color,1001,color-9);
      hm->SetFillColorAlpha(color-9, 0.35); // 35% transparent
      tdrDraw(hd,"P",marker,color,kSolid,color,1001,color-9);
      if (marker!=kFullDiamond && marker!=kOpenDiamond) hd->SetMarkerSize(0.5);

      // Add basic labelling
      TLatex *tex = new TLatex();
      tex->SetTextSize(0.040);
      tex->SetNDC();

      double pl = gPad->GetLeftMargin();
      double pt = gPad->GetTopMargin();
      double pr = gPad->GetRightMargin();
      double pb = gPad->GetBottomMargin();
      double xpos = pl +0.7*(1-pl-pr);
      double ypos = 1-pt - 0.08*(1-pt-pb);
      double ypos2 = 1-pt - 0.13*(1-pt-pb);
      if (y1==0) tex->DrawLatex(xpos,ypos,Form("|#eta|<%1.1f (%s)",y2,ctrg));
      else tex->DrawLatex(xpos,ypos,Form("%1.1f<|#eta|<%1.1f (%s)",y1,y2,ctrg));
      tex->DrawLatex(xpos,ypos2,Form("%1.0f<p_{T}<%1.0f GeV",minpt,maxpt));

      const char *name = "trpu";

      TLine *l = new TLine();
      l->SetLineStyle(kDotted);
      l->DrawLine(0,1,1,1);

      // Add labelling for TP and DP methods
      TH1D *hl1 = (TH1D*)hd->Clone("hl1");
      hl1->SetMarkerStyle(kFullCircle);
      leg->AddEntry(hd,"tag&probe","P");

      gPad->RedrawAxis();

      string sfile = Form("pu_y%02d-%02d_%s.pdf",int(10.*y1),int(10*y2),ctrg);

      c1->cd(0);
      c1->SaveAs(Form("pdf/%s",sfile.c_str()));
      c1->Print(Form("pdf/%s.gif+100",etafile.c_str()));

      filenames.push_back(Form("\\includegraphics[width=0.30\\textwidth]{%s}",sfile.c_str()));
    } // itrg
  } // ieta

  for (unsigned int i = 0; i != filenames.size(); ++i) {
    cout << filenames[i] << endl;
  }
}
