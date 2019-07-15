#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"

#include "tdrstyle_mod18.C"
#include "ptresolution.h"
#include "settings.h"
#include "tools.h"

#include <iostream>

void drawDagostini(string type) {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile(Form("output-%s-3.root",type.c_str()),"READ");
  assert(f && !f->IsZombie());

  assert(f->cd("Standard"));
  f->cd("Standard");
  TDirectory *din = f->GetDirectory("Standard"); assert(din);
  curdir->cd();

  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);   // 800
  //c1->SetTopMargin(0.20);
  c1->Divide(3,3);//,0.1,0.00); 3,2

  TCanvas *c1b = new TCanvas("c1b","c1b",600,600);

  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->SetTopMargin(0.10);
  c2->Divide(3,3,-1,-1);

  TCanvas *c3 = new TCanvas("c3","c3",1200,1200);
  c3->SetTopMargin(0.10);
  c3->Divide(3,3,-1,-1);

  TH1D *h = new TH1D("h",";p_{T} (GeV);Unfolding correction",
		     int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.45);
  h->SetMaximum(1.15);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.25,0.74,0.45,0.97,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);

  const int ny = 7;
  double y1 = 0; double y2 = 0;
  for (int iy = 0; iy != ny; ++iy) {
    if (iy==6) {y1 = 3.2; y2 = 4.7;}
    else {y1 = 0.5*iy;  y2 = 0.5*(iy+1);}

    assert(din->cd(Form("Eta_%1.1f-%1.1f",y1,y2)));
    din->cd(Form("Eta_%1.1f-%1.1f",y1,y2));
    TDirectory *d = din->GetDirectory(Form("Eta_%1.1f-%1.1f",y1,y2)); assert(d);

    TH1D *hreco = (TH1D*)d->Get("hreco"); assert(hreco);
    TH2D *h2resp = (TH2D*)d->Get("mtu"); assert(h2resp);
    TGraphErrors *gfwd = (TGraphErrors*)d->Get("gfold_fwd"); assert(gfwd);
    TGraphErrors *gbin = (TGraphErrors*)d->Get("gfold_bin"); assert(gbin);
    TGraphErrors *gsvd = (TGraphErrors*)d->Get("gfold_svd"); assert(gsvd);
    TGraphErrors *gbayes = (TGraphErrors*)d->Get("gfold"); assert(gbayes);
    curdir->cd();

    // Find the low end of reconstructed data used in unfolding
    double ptmin(0);
    for (int i = 1; i != hreco->GetNbinsX()+1 && !ptmin; ++i) {
      if (hreco->GetBinContent(i)>0) ptmin = hreco->GetBinLowEdge(i);
    }

    c1->cd(iy+1);
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetPalette(1);
    gPad->SetRightMargin(0.10);
    h2resp->GetXaxis()->SetMoreLogLabels();
    h2resp->GetXaxis()->SetNoExponent();
    h2resp->GetXaxis()->SetTitle("Measured p_{T,reco} (GeV)");
    h2resp->GetXaxis()->SetTitleOffset(1.5);
    h2resp->GetYaxis()->SetMoreLogLabels();
    h2resp->GetYaxis()->SetNoExponent();
    h2resp->GetYaxis()->SetTitle("True p_{T,gen} (GeV)");
    h2resp->GetYaxis()->SetTitleOffset(2.0);
    h2resp->GetZaxis()->SetRangeUser(1e-3,0.9999);//1.05);
    h2resp->DrawClone("COLZ");
    //h2resp->Draw("SAMEBOX");

    tex->SetTextSize(0.045);
    tex->SetTextAlign(31); // align right
    tex->SetNDC(kTRUE);
    tex->DrawLatex(0.8, 0.2, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (y1==0) {
      c1b->cd();
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetRightMargin(0.13);//0.11);
      gPad->SetBottomMargin(0.14);//0.10);
      gPad->SetLeftMargin(0.18);//0.10);
      //h2resp->GetXaxis()->SetTitleOffset(1.4);
      //TH2D *h2 = (TH2D*)h2resp->DrawClone("COLZ");

      // copy TH2D over to a fresh one to get default graphics style
      TH2D *h2 = new TH2D("h2",";Measured Jet p_{T} (GeV);"
			  "True Jet p_{T} (GeV)",
			  h2resp->GetNbinsX(),
			  h2resp->GetXaxis()->GetXbins()->GetArray(),
			  h2resp->GetNbinsY(),
			  h2resp->GetYaxis()->GetXbins()->GetArray());
      for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
	for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	  h2->SetBinContent(i, j, h2resp->GetBinContent(i, j));
	  //h2->SetBinError(i, j, h2resp->GetBinError(i, j));
	} // for j
      } // for i
      h2->Draw("COLZ");
      h2->GetXaxis()->SetMoreLogLabels();
      h2->GetXaxis()->SetNoExponent();
      //h2->GetYaxis()->SetMoreLogLabels();
      h2->GetYaxis()->SetNoExponent();
      h2->GetYaxis()->SetTitleOffset(1.5);
      h2->GetZaxis()->SetRangeUser(1e-3,0.9999);

      //    h2->GetXaxis()->SetRangeUser(jp::fitptmin,jp::xmax);//1327.);
      // h2->GetYaxis()->SetRangeUser(jp::recopt,jp::xmax);
      
      h2->GetXaxis()->SetRangeUser(jp::unfptminreco,jp::xmax);//1327.);
      h2->GetYaxis()->SetRangeUser(jp::unfptmingen,jp::xmax);
      //h2resp->Draw("SAME COLZ");

      //tex->DrawLatex(0.8, 0.15, _algo=="AK7" ?
      tex->DrawLatex(0.8, 0.20, string(jp::algo)=="AK7" ?
		     "Anti-k_{T} R=0.8, |y|<0.5" :
		     "Anti-k_{T} R=0.4, |y|<0.5");
      //cmsPrel(_lumi);
      //cmsPrel(0); // simulation
    }

    c2->cd(iy+1);
    gPad->SetLogx();

    h->SetMinimum(0.1);
    h->SetMaximum(1.15);
    h->SetXTitle(iy==ny-1 ? "p_{T} (GeV)" : "");
    h->SetYTitle(iy==0 ? "Unfolding correction" : "");
    h->DrawClone("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDotted);
    l->DrawLine(jp::unfptminnlo,0.15,jp::unfptminnlo,1.15); // minimum pt in unfolding

    //  l->SetLineStyle(kDashed);
    //  l->DrawLine(jp::xminpas,0.15,jp::xminpas,1.15); // Line showing data points included in results

    gfwd->SetName("gfwd");
    gfwd->SetLineWidth(2);
    gfwd->SetLineColor(kRed);
    gfwd->Draw("SAMEL");

    gbin->SetName("gbin");
    gbin->SetMarkerStyle(kOpenSquare);
    gbin->SetMarkerColor(kGreen+1);
    gbin->SetLineColor(kGreen+1);
     gbin->Draw("SAMEPz");

    gbayes->SetName("gbayes");
    gbayes->SetMarkerStyle(kFullCircle);
    gbayes->Draw("SAMEP");

    gsvd->SetName("gsvd");
    gsvd->SetMarkerStyle(kOpenDiamond);
    gsvd->SetMarkerColor(kCyan+1);
    gsvd->SetLineColor(kCyan+1);
     gsvd->Draw("SAMEP");

    tex->SetTextSize(iy<3 ? 0.053 : 0.045);
    tex->SetTextAlign(21); // align middle
    tex->SetNDC(kFALSE);
    tex->DrawLatex(150, 0.5, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (iy==2) {
      leg->AddEntry(gbayes,"RooUnfoldBayes","P");
      leg->AddEntry(gbin,"RooUnfoldBinByBin","P");
       leg->AddEntry(gsvd,"RooUnfoldSvd","P");
      leg->AddEntry(gfwd,"Forward smearing","L");
      leg->Draw();
    }

    c3->cd(iy+1);
    gPad->SetLogx();

    h->SetMinimum(0.90+0.0001);
    h->SetMaximum(1.10-0.0001);
    h->SetXTitle(iy==ny-1 ? "p_{T} (GeV)" : "");
    h->SetYTitle(iy==0 ? "Correction / Forward smearing" : "");
    h->DrawClone("AXIS");

    l->SetLineStyle(kDotted);
    l->DrawLine(ptmin,0.9,ptmin,1.1);
    l->SetLineStyle(kDashed);
    //l->DrawLine(56,0.9,56,1.1); // v2
    l->DrawLine(jp::xminpas,0.9,jp::xminpas,1.1);

    TGraphErrors *grfwd = tools::ratioGraphs(gfwd, gfwd);
    grfwd->Draw("SAMEL");
    TGraphErrors *grbin = tools::ratioGraphs(gbin, gfwd);
    grbin->Draw("SAMEPz");
    TGraphErrors *grbayes = tools::ratioGraphs(gbayes, gfwd);
    grbayes->Draw("SAMEPz");
    TGraphErrors *grsvd = tools::ratioGraphs(gsvd, gfwd);
    grsvd->Draw("SAMEPz");

    tex->DrawLatex(150,0.91, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (iy==2) leg->Draw();

  }

  const char *a = string(jp::algo).c_str();
  const char *t = type.c_str();

  c1->cd(3);
  tex->SetTextSize(0.045);
  tex->SetNDC(kTRUE);
  tex->DrawLatex(0.35, 0.85, Form("%s %s",t,a));
  c1->cd(0);
  //cmsPrel(type=="DATA" ? _lumi : 0, true);
  c1->SaveAs(Form("pdf/roounfold_matrix_%s_%s.pdf",a,t));
  c1b->SaveAs(Form("pdf/roounfold_matrix0_%s_%s.pdf",a,t));

  c2->cd(2);
  tex->SetTextSize(0.053);
  tex->DrawLatex(0.50, 0.85, Form("%s %s",t,a));
  c2->cd(0);
  //cmsPrel(type=="DATA" ? jp::lumi : 0, true);
  c2->SaveAs(Form("pdf/roounfold_comparison_%s_%s.pdf",a,t));

  c3->cd(2);
  tex->DrawLatex(0.50, 0.85, Form("%s %s",t,a));
  c3->cd(0);
  //cmsPrel(type=="DATA" ? jp::lumi : 0, true);
  c3->SaveAs(Form("pdf/roounfold_ratiotofwd_%s_%s.pdf",a,t));


} // drawDagostini
