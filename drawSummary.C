// Purpose: Draw publication quality summary plots of analysis
// Author:  mikko.voutilainen@cern.ch
// Created: April 21, 2010
// Updated: April 21, 2010
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TROOT.h"
#include "TMath.h"
#include "TProfile.h"
#include "TBox.h"

#include <iostream>
#include <fstream>
#include <string>

#include "tools.h"
#include "tdrstyle_mod18.C"

using namespace std;

// ranges for drawing
#include "settings.h"
const bool _final = false; // cmsFinal

const int kRightHatch = 3654;
const int kLeftHatch = 3645;
const int kVertHatch = 3699;
const int kHorizHatch = 3600;
const int kCrossHatch = 3644;

// possible commands
void drawEtaSpectra(string type);
void drawDataTheoryRatioLog();
void drawDataTheoryRatio();
void drawDataTheoryRatio3x2();//bool ak7 = false, bool v38x = false);
void drawTheoryDataRatio(bool ak7 = false, bool v38x = false);
void drawUnfoldSummary(string type, string algo="PF");
void drawJetIDSummary(string type, string idtype = "JetID");
void drawNPCorr(string type);

// helper functions
//void GetPoint(TGraphErrors *g, int n, double &x, double &y,
//      double &ex, double &ey);
//void SetPoint(TGraphErrors *g, int n, double x, double y,
//      double ex, double ey);
void scaleGraph(TGraphErrors *g, double scale);
void offsetGraph(TGraphErrors *g, double offset);
void offsetGraph(TH1D *g, double offset);
void ratioGraph(TGraphErrors *g, TF1 *f);
void ratioGraph(TGraphErrors *g, TH1D *h);
void pruneGraph(TGraphErrors *g, double errmax, double ymin=0.,
		double xmin=-1e10, double xmax=1e10);
void rangeGraph(TGraph *g, double xmin, double xmax);
TGraphErrors *makeBand(TGraph *f, const TH1D *epl, const TH1D *emn,
		       double xmin, double xmax, bool hist = false);
TGraphErrors *makeBand(TF1 *f, const TH1D *epl, const TH1D *emn,
		       double xmin, double xmax);


void drawEtaSpectra(string type) {

  TDirectory *curdir = gDirectory;

  const char *t = type.c_str();

  // Unfolded data
  TFile *fin = new TFile(Form("output-%s-3.root",t),"READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  // Uncertainties
  TFile *fin2 = new TFile(Form("output-%s-4.root",t),"READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  // Theory
  TFile *fnlo = new TFile(Form("output-%s-2c.root",t),"READ");
  assert(fnlo && !fnlo->IsZombie());
  assert(fnlo->cd("Standard"));
  TDirectory *dnlo = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");

  double c = 10;
  double scale[] = {pow(c,5),pow(c,4),pow(c,3),pow(c,2),pow(c,1),pow(c,0.)};
  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();
  c1->SetLogy();

  double xmaxeta = jp::xmax;
  TH1D *h = new TH1D("h","",int(xmaxeta-jp::xmin),jp::xmin,xmaxeta);
  h->SetMinimum(1.0001e-4);
  h->SetMaximum(5e11);
  h->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
  h->GetYaxis()->SetTitleOffset(1.15);// for ^{2}
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->Draw("AXIS");

  TLegend *leg = new TLegend(0.58,0.69,0.90,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();

  TLegend *leg2 = new TLegend(0.19,0.22,0.53,0.30,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);
  leg2->Draw();

  TLatex *tjet = new TLatex(0.19,0.18,Form("Anti-k_{T} R=%1.1f",
					   string(jp::algo)=="AK7" ? 0.7 : 0.5));
  tjet->SetTextSize(0.04);
  tjet->SetNDC();
  tjet->Draw();

  //cmsPrel(jp::isdt ? jp::lumi : 0);

  vector<TH1D*> vhnlo(etas.size());
  vector<TGraphErrors*> vnlo(etas.size());
  vector<TGraphErrors*> vnloup(etas.size());
  vector<TGraphErrors*> vnlodw(etas.size());
  vector<TGraphErrors*> vgerr(etas.size());
  vector<TH1D*> vheup(etas.size());
  vector<TH1D*> vhedw(etas.size());
  vector<TGraphErrors*> vdata(etas.size());
  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;

    assert(dnlo->cd(etas[i].c_str()));
    TDirectory *dt = gDirectory;

    assert(din2->cd(etas[i].c_str()));
    TDirectory *d2 = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load and scale central values
    TGraphErrors *g = (TGraphErrors*)d->Get("gcorrpt"); assert(g);
    scaleGraph(g, scale[i]);
    TF1 *f = (TF1*)d->Get("fus"); assert(f);
    f->SetParameter(0, f->GetParameter(0)*scale[i]);
    vdata[i] = g;

    // Load theory prediction
    TH1D *htheory = (TH1D*)dt->Get("hnlo"); assert(htheory);
    TH1D *hthup = (TH1D*)dt->Get("hsysup"); assert(hthup);
    TH1D *hthdw = (TH1D*)dt->Get("hsysdw"); assert(hthdw);
    TGraphErrors *gnlo = new TGraphErrors(0);
    TGraphErrors *gnloup = new TGraphErrors(0);
    TGraphErrors *gnlodw = new TGraphErrors(0);
    for (int j = 0; j != g->GetN(); ++j) {
      double x = g->GetX()[j];
      if (x >= jp::xmin) {//jp::xminpas) {
	double y = htheory->GetBinContent(htheory->FindBin(x));
	int n = gnlo->GetN();
	gnlo->SetPoint(n, x, y * scale[i]);
	double yup = hthup->GetBinContent(hthup->FindBin(x));
	gnloup->SetPoint(n, x, y * (1 + yup) * scale[i]);
	double ydw = hthdw->GetBinContent(hthdw->FindBin(x));
	gnlodw->SetPoint(n, x, y * (1 + ydw) * scale[i]);
      }
    }

    htheory->GetXaxis()->SetRangeUser(g->GetX()[0],g->GetX()[g->GetN()-1]);
    htheory->SetLineWidth(2);
    htheory->SetLineColor(kRed);
    htheory->Scale(scale[i]);
    vhnlo[i] = htheory;

    gnlo->SetLineWidth(2);
    gnlo->SetLineColor(kRed);
    gnloup->SetLineColor(kRed);
    gnlodw->SetLineColor(kRed);
    vnlo[i] = gnlo;
    vnloup[i] = gnloup;
    vnlodw[i] = gnlodw;

    // Find last entry
    double ptmax, yptmax;
    g->GetPoint(g->GetN()-1, ptmax, yptmax);

    // Produce uncertainty band around ansatz fit
    TH1D *hpl = (TH1D*)d2->Get("htot_pl"); assert(hpl);
    TH1D *hmn = (TH1D*)d2->Get("htot_mn"); assert(hmn);
    ptmax = hpl->GetBinLowEdge(hpl->FindBin(ptmax)+1);

    TGraphErrors *gerr = (jp::centerOnAnsatz ? makeBand(f, hpl, hmn, jp::xmin, jp::xmax) :
			  makeBand(g, hpl, hmn, jp::xmin, jp::xmax));
    vgerr[i] = gerr;
    vheup[i] = hpl;
    vhedw[i] = hmn;

    // Find reasonable upper range for fit
    // at 1 TeV range, we require 1 event in a bin of O(100 GeV) width
    // => update to ptmax
    //f->SetRange(1., 3500./cosh(etamin));
    //double x = f->GetX(scale[i]*1./100.);
    f->SetRange(jp::xmin, ptmax);//x);

    g->SetMarkerStyle(marker[i]);
    gerr->SetFillStyle(1001);
    gerr->SetFillColor(kYellow+1);
    f->SetLineColor(kRed);

    gerr->Draw(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    if (gnlo) gnlo->Draw("SAME L");
    else      f->Draw("SAME");
    if (gnloup) gnloup->Draw("SAME L");
    if (gnlodw) gnlodw->Draw("SAME L");
    g->Draw("SAMEP");

    const char *seta = (etamin==0 ? Form("|y|<%1.2g    ",etamax) :
			Form("%1.2g#leq|y|<%1.2g",etamin,etamax));
    //const char *smult = (scale[i]==1 ? "" : Form(" (#times%1.0f)",scale[i]));
    const char *smult = Form(" (#times %1.0f^{%1.0f})",c,log(scale[i])/log(c));
    leg->AddEntry(g,Form("%s%s",seta,smult),"P");

    if (i==0) {
      assert(gnlo);
      leg2->AddEntry(gnlo, "NLO#otimesNP theory", "L"); // cmsFinal
      leg2->AddEntry(gerr, "Exp. uncertainty", "F");
    }
  } // for i


  cout << endl;
  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;
    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Print ansatz fits
    TF1 *fs = (TF1*)d->Get("fs"); assert(fs);
    TF1 *f = (TF1*)d->Get("fus"); assert(f);
    if (fs) {

      cout << "fs:" << endl;
      cout << Form("$%1.1f<|y|<%1.1f$ & %1.1f$\\pm$%1.1f & %1.2f$\\pm$%1.2f &"
		   " %1.1f$\\pm$%1.1f & %1.2f$\\pm$%1.2f & "
		   " %1.1f / %d \\\\", etamin, etamax,
		   fs->GetParameter(1)/1e12, fs->GetParError(1)/1e12,
		   fs->GetParameter(3), fs->GetParError(3),
		   fs->GetParameter(4), fs->GetParError(4),
		   fs->GetParameter(2), fs->GetParError(2),
		   fs->GetChisquare(), fs->GetNDF())
	   << endl;
      cout << "fus:" << endl;
      cout << Form("$%1.1f<|y|<%1.1f$ & %1.1f$\\pm$%1.1f & %1.2f$\\pm$%1.2f &"
		   " %1.1f$\\pm$%1.1f & %1.2f$\\pm$%1.2f & "
		   " %1.1f / %d \\\\", etamin, etamax,
		   f->GetParameter(1)/1e12, f->GetParError(1)/1e12,
		   f->GetParameter(3), f->GetParError(3),
		   f->GetParameter(4), f->GetParError(4),
		   f->GetParameter(2), f->GetParError(2),
		   f->GetChisquare(), f->GetNDF())
	   << endl;
    }
  } // for i
  cout << endl;
  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;
    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Print ansatz fits
    TF1 *fs = (TF1*)d->Get("fs"); assert(fs);
    cout << Form("   {%1.3g, %1.3g, %1.3g, %1.3g}, // %1.2g<|y|<%1.2g:"
		 " %1.1f/%d", fs->GetParameter(1),
		 fs->GetParameter(3), fs->GetParameter(4),
		 fs->GetParameter(2), etamin, etamax,
		 fs->GetChisquare(), fs->GetNDF()) << endl;
  } // for i

  const char *a = string(jp::algo).c_str();
  //if(_eps) c1->SaveAs(Form("eps/summaryEtaSpectra_%s_%s.eps",a,t));
  if(jp::pdf) c1->SaveAs(Form("pdf/summaryEtaSpectra_%s_%s.pdf",a,t));
  //if(_gif) c1->SaveAs(Form("gif/summaryEtaSpectra_%s_%s.gif",a,t));
  //if(_png) c1->SaveAs(Form("png/summaryEtaSpectra_%s_%s.png",a,t));

  ///////////////////////////////////////////////////////////////
  // for PAS, cut the NLO line at pT>18 GeV, color band yellow
  for (unsigned int i = 0; i != vnlo.size(); ++i) {
    pruneGraph(vnlo[i],1e10,-1e10,jp::xminpas,1e10);
    pruneGraph(vnloup[i],1e10,-1e10,jp::xminpas,1e10);
    pruneGraph(vnlodw[i],1e10,-1e10,jp::xminpas,1e10);
    pruneGraph(vgerr[i],1e10,-1e10,jp::xminpas,1e10);
    vgerr[i]->SetFillColor(kYellow+1);
  } // for i
  //h->SetMaximum(5e11);
  //h->SetMaximum(0.99999e11);
  h->SetMaximum(1e11);
  h->SetMinimum(1e-7);
  h->GetXaxis()->SetRangeUser(jp::xminpas,xmaxeta);

  // Redraw fully
  h->Draw("AXIS");
  for (unsigned int i = 0; i != vgerr.size(); ++i) vgerr[i]->Draw("SAME E3");
  for (unsigned int i = 0; i != vhnlo.size(); ++i) vhnlo[i]->Draw("SAME HIST][");
  //for (unsigned int i = 0; i != vnlo.size(); ++i) vnlo[i]->Draw("SAME L");
  for (unsigned int i = 0; i != vdata.size(); ++i) vdata[i]->DrawClone("SAME Pz");

  leg->Draw();
  leg2->Draw();
  tjet->Draw();

  TLatex *tex = new TLatex();
  //tex->SetTextSize(0.045);
  tex->SetTextSize(0.04);
  tex->SetNDC();
  //tex->DrawLatex(0.19,0.33+0.03,"#mu_{R} = #mu_{F} = p_{T}");
  tex->DrawLatex(0.19,0.30+0.03,"#mu_{R} = #mu_{F} = p_{T}");

  //if (jp::isdt and _final) cmsFinal(jp::lumi);
  //else cmsPrel(jp::isdt ? jp::lumi : 0);

  gPad->RedrawAxis();

  // Save the histograms for redrawing PAS plots (here, so pT<30 GeV clipped)
//   if (_pas) {
//     TDirectory *curdirpas = gDirectory;
//     TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
//     for (unsigned int i = 0; i != etas.size(); ++i) {

//       TGraphErrors *gnlo = (TGraphErrors*)vnlo[i]->Clone();
//       scaleGraph(gnlo, 1./scale[i]);
//       if (gnlo) gnlo->Write(Form("fig6_%d_gnlo",i),TObject::kOverwrite);
//       delete gnlo;

//       TGraphErrors *g = (TGraphErrors*)vdata[i]->Clone();
//       scaleGraph(g, 1./scale[i]);
//       if (g) g->Write(Form("fig6_%d_gpt",i),TObject::kOverwrite);
//       delete g;
//     }
//     fpas->Close();
//     curdirpas->cd();
//   } // _pas

  //if(_eps) c1->SaveAs(Form("eps/summaryEtaSpectraPAS_%s_%s.eps",a,t));
  if(jp::pdf) c1->SaveAs(Form("pdf/summaryEtaSpectraPAS_%s_%s.pdf",a,t));
  //if(_gif) c1->SaveAs(Form("gif/summaryEtaSpectraPAS_%s_%s.gif",a,t));
  //if(_png) c1->SaveAs(Form("png/summaryEtaSpectraPAS_%s_%s.png",a,t));

  // Print out the Durham tables
  gROOT->ProcessLine(".! mkdir txt");
  for (unsigned int i = 0; i != vdata.size(); ++i) {
    ofstream fout(Form("txt/CMS_incjet_y%1.1f_%1.1f_%s.txt",0.5*i,0.5*(i+1),
		       jp::algo),ios::out);
    fout << Form("Inclusive jet double-differential cross sections in the |rapidity| range %1.1f to %1.1f, using a jet resolution R value of %1.1f.",0.5*i,0.5*(i+1),string(jp::algo)=="AK7" ? 0.7 : 0.5) << endl;
    fout << "RE : P P --> JET X" << endl;
    fout << Form("ABS(YRAP) : %1.1f TO %1.1f",0.5*i,0.5*(i+1)) << endl;
    fout << Form("R : %1.1f",string(jp::algo)=="AK7" ? 0.7 : 0.5) << endl;
    fout << "SQRT(S) : 7 TeV" << endl;
    fout << "x-axis header: PT IN GEV" << endl;
    fout << "y-axis header: D2(SIG)/DPT/DYRAP IN PB/GEV" << endl;
    fout << endl;
    fout << Form("  %-6s %-5s %-5s    %-10s  %-21s %-21s",
		 "x","xlow","xhigh","y","dy(stat)","dy(sys)") << endl;
    fout << endl;
    TGraphErrors *g = vdata[i];
    TH1D *hup = vheup[i];
    TH1D *hdw = vhedw[i];
    for (int j = 0; j != g->GetN(); ++j) {

      double x = g->GetX()[j];
      if (x>18) {
	int ibin = hup->FindBin(x);
	double xminb = hup->GetBinLowEdge(ibin);
	double xmaxb = hup->GetBinLowEdge(ibin+1);
	double y = g->GetY()[j];
	double ey = g->GetEY()[j] / y * 100.;
	double sysup = hup->GetBinContent(ibin)*100.;
	double sysdw = hdw->GetBinContent(ibin)*100.;
	fout << Form("%6.1f %5.0f %5.0f  %10.3g  %+5.1f %+5.1f (percent)"
		     "  %+5.1f %+5.1f (percent)",
		     x, xminb, xmaxb, y/scale[i], ey, -ey, sysup, sysdw)
	     << endl;
      }
    } // for j
    fout.close();
  } // for i

  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawEtaSpectra


void drawDataTheoryRatioLog() {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile("output-DATA-4c.root","READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  TFile *fin2 = new TFile("output-DATA-4.root","READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  TFile *fnlo = new TFile("output-DATA-5.root","READ");
  assert(fnlo && !fnlo->IsZombie());
  assert(fnlo->cd("Standard"));
  TDirectory *dnlo = gDirectory;

  TFile *finmc = new TFile("output-MC-4c.root","READ");
  assert(finmc && !finmc->IsZombie());
  assert(finmc->cd("Standard"));
  TDirectory *dinmc = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");

  double c = 10.;
  double scale[] = {pow(c,5),pow(c,4),pow(c,3),pow(c,2),pow(c,1),pow(c,0.)};
  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();
  c1->SetLogy();

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(pow(c,-1.5));
  h->SetMaximum(pow(c,6));
  h->GetYaxis()->SetTitle("Data / NLO theory");
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->Draw("AXIS");

  TLegend *leg = new TLegend(11.,0.5,50.,8e5,"","br");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.033);

  TLegend *leg2 = new TLegend(100.,0.04,400.,0.7,"","br");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.033);

  double ymin = 0.5;
  double ymax = 5e5;
  TLine *lptmin = new TLine(50.,ymin,50.,ymax);
  lptmin->SetLineStyle(kDashed);
  lptmin->SetLineColor(kBlue);
  lptmin->Draw();

  TLatex *tptmin = new TLatex(1.1*50., 0.9*ymax,"p_{T} > 50 GeV");
  tptmin->SetTextSize(0.033);
  tptmin->SetTextColor(kBlue);
  tptmin->Draw("SAME");

  TLine *lxmin = new TLine(50./980.*4000.,ymin,50./980.*4000.,ymax);
  lxmin->SetLineStyle(kDotted);
  lxmin->SetLineColor(kBlue);
  lxmin->Draw();

  TLatex *txmin = new TLatex(1.1*50./980.*4000, 0.9*ymax,"x > 0.0051");
  txmin->SetTextSize(0.033);
  txmin->SetTextColor(kBlue);
  txmin->Draw();

  //cmsPrel(jp::lumi);

  TLatex *tjet = new TLatex(0.17,0.22,"Anti-k_{T} R=0.5");// PF");
  tjet->SetTextSize(0.045);
  tjet->SetNDC();
  tjet->Draw();

  vector<TF1*> vnlo(etas.size());
  //vector<TH1D*> verr(etas.size());
  vector<TGraphErrors*> vgerr(etas.size());
  vector<TGraphErrors*> vmc(etas.size());
  vector<TGraphErrors*> vdt(etas.size());
  vector<TH1D*> hpas(etas.size());

  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;

    assert(din2->cd(etas[i].c_str()));
    TDirectory *d2 = gDirectory;

    assert(dnlo->cd(etas[i].c_str()));
    TDirectory *dt = gDirectory;

    assert(dinmc->cd(etas[i].c_str()));
    TDirectory *dmc = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TGraphErrors *g = (TGraphErrors*)d->Get("gcorrpt"); assert(g);
    TGraphErrors *gmc = (TGraphErrors*)dmc->Get("gcorrpt"); assert(gmc);

    // Load theory
    TF1 *fus = (TF1*)d->Get("fus"); assert(fus);
    TH1D *htheory = (TH1D*)dt->Get("hnlo"); assert(htheory);

    // For PAS: evaluate ratio of ansatz and theory
    TH1D *hdtoth = (TH1D*)htheory->Clone(Form("fig7_%d_dataovertheory",i));
    for (int j = 1; j != hdtoth->GetNbinsX()+1; ++j) {
      double xminb = hdtoth->GetBinLowEdge(j);
      double xmaxb = min(4000./cosh(0.5*i), hdtoth->GetBinLowEdge(j+1));
      double y = htheory->GetBinContent(j);
      if (xminb*cosh(0.5*i) < 4000. && y) {
	double fy = fus->Integral(xminb, xmaxb) / (xmaxb - xminb);
	hdtoth->SetBinContent(j, fy / y);
	hdtoth->SetBinError(j, 0.);
      }
      else {
	hdtoth->SetBinContent(j ,0.);
	hdtoth->SetBinError(j, 0.);
      }
    }
    hpas[i] = hdtoth;

    if (htheory) ratioGraph(g, htheory);
    else         ratioGraph(g, fus);
    scaleGraph(g, scale[i]);
    if (htheory) ratioGraph(gmc, htheory);
    else         ratioGraph(gmc, fus);
    scaleGraph(gmc, scale[i]);
    pruneGraph(gmc, 0.50, 0., 10., jp::xmax);
    gmc->SetLineWidth(2);
    gmc->SetLineStyle(kDashed);

    // Produce uncertainty band around ansatz fit
    TH1D *hpl = (TH1D*)d2->Get("htot_pl"); assert(hpl);
    TH1D *hmn = (TH1D*)d2->Get("htot_mn"); assert(hmn);
    TF1 *f = new TF1(Form("f%d",i),"[0]", jp::xmin, jp::xmax);
    f->SetParameter(0, scale[i]);
    TGraphErrors *gerr = (jp::centerOnAnsatz ?
			  makeBand(f, hpl, hmn, jp::xmin, f->GetXmax()) :
			  makeBand(g, hpl, hmn, jp::xmin, f->GetXmax()));

    g->SetMarkerStyle(marker[i]);
    gmc->SetMarkerStyle(marker[i]);
    gmc->SetMarkerColor(kBlue);
    gmc->SetLineColor(kBlue);
    gerr->SetFillStyle(1001);
    gerr->SetFillColor(kYellow+1);
    f->SetLineColor(kRed);

    gerr->Draw(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    f->Draw("SAME");
    gmc->Draw("SAMEL");
    g->Draw("SAMEP");

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    const char *smult = (scale[i]==1 ? "" : Form(" (#times%1.0f^{%1.0f})", c,
						 log(scale[i])/log(c)));
    leg->AddEntry(g,Form("%s%s",seta,smult),"P");
    leg->AddEntry(gmc," ","");

    if (i==0) {
      if (htheory) leg2->AddEntry(f, "NLO theory", "L");
      else         leg2->AddEntry(f, "ansatz fit to data", "L");
      leg2->AddEntry(gmc,"Pythia D6T QCD", "L");
      leg2->AddEntry(gerr, "Exp. uncertainty", "F");
    }
    assert(gmc->GetN());

    vnlo[i] = f;
    vgerr[i] = gerr;
    vmc[i] = gmc;
    vdt[i] = g;
  } // for i

  // Get the tick marks on top of the bands
  h->Draw("AXIS SAME");
  leg->Draw();
  leg2->Draw();

  //if(_eps) c1->SaveAs("eps/summaryDataTheoryRatioLog.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryDataTheoryRatioLog.pdf");
  //if(_gif) c1->SaveAs("gif/summaryDataTheoryRatioLog.gif");
  //if(_png) c1->SaveAs("png/summaryDataTheoryRatioLog.png");

  /////////////////////////////////////////////////////
  // for PAS, cut at pT>30 GeV, color error band yellow
  h->GetXaxis()->SetRangeUser(30,jp::xmax);
  h->SetMinimum(2e-2);
  h->Draw("AXIS");
  for (unsigned int i = 0; i != etas.size(); ++i) {
    vgerr[i]->SetFillColor(kYellow+1);
    vgerr[i]->Draw(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    vnlo[i]->Draw("SAME");
    vmc[i]->Draw("SAME L");
    vdt[i]->Draw("SAME P");
  }
  h->Draw("AXIS SAME");
  tjet->Draw();
  leg->Draw();
  leg2->Draw();

  //if(_eps) c1->SaveAs("eps/summaryDataTheoryRatioLogPAS.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryDataTheoryRatioLogPAS.pdf");
  //if(_gif) c1->SaveAs("gif/summaryDataTheoryRatioLogPAS.gif");
  //if(_png) c1->SaveAs("png/summaryDataTheoryRatioLogPAS.png");

  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawDataTheoryRatioLog


void drawDataTheoryRatio() {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile("output-DATA-4c.root","READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  TFile *fin2 = new TFile("output-DATA-4.root","READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  TFile *fnlo = new TFile("output-DATA-5.root","READ");
  assert(fnlo && !fnlo->IsZombie());
  assert(fnlo->cd("Standard"));
  TDirectory *dnlo = gDirectory;

  //TFile *finmc = new TFile("output-MC-3a.root","READ");
  TFile *finmc = new TFile("output-MC-4c.root","READ");
  assert(finmc && !finmc->IsZombie());
  assert(finmc->cd("Standard"));
  TDirectory *dinmc = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");

  double off = 1.;
  double offset[] = {5*off, 4*off, 3*off, 2*off, off, 0.};
  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.0);
  h->SetMaximum(9.0*off);
  h->GetYaxis()->SetTitle("Data / NLO theory");
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetNdivisions(0,int(9*off)*2,0);
  h->Draw("AXIS");

  TLegend *leg = new TLegend(0.58,0.70,0.93,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);//0.045);

  TLegend *leg2 = new TLegend(0.18,0.78,0.48,0.93,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);

  TLatex *ta = new TLatex();
  ta->SetTextFont(42);
  ta->SetTextSize(0.04);

  vector<TF1*> vnlo(etas.size());
  vector<TGraphErrors*> vgerr(etas.size());
  vector<TGraphErrors*> vmc2(etas.size());
  vector<TGraphErrors*> vgnlo(etas.size());
  vector<TGraphErrors*> vdt(etas.size());

  TLatex *t = new TLatex();
  t->SetTextSize(0.05);
  t->SetTextFont(42);
  t->SetTextAlign(31); // align right

  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;

    assert(din2->cd(etas[i].c_str()));
    TDirectory *d2 = gDirectory;

    assert(dnlo->cd(etas[i].c_str()));
    TDirectory *dt = gDirectory;

    assert(dinmc->cd(etas[i].c_str()));
    TDirectory *dmc = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TGraphErrors *g = (TGraphErrors*)d->Get("gcorrpt"); assert(g);
    TGraphErrors *gmc = (TGraphErrors*)dmc->Get("gcorrpt"); assert(gmc);

    // Load theory
    TF1 *fus = (TF1*)d->Get("fus"); assert(fus);
    TH1D *htheory = (TH1D*)dt->Get("hnlo"); assert(htheory);
    TGraphErrors *gtheory = (TGraphErrors*)dt->Get("gnlofit"); assert(gtheory);
    TF1 *ftheory = (TF1*)dt->Get("fnlo"); assert(ftheory);

    // Ansatz ratio to theory
    TF1 *fdt = new TF1(Form("fdt%d",i),"[0]*exp([1]/x)*pow(x,[2])"
		       "*pow(1-x*cosh([4])/4000.,[3]) / "
		       "([5]*exp([6]/x)*pow(x,[7]+[9]*log(0.01*x))"
		       "*pow(1-x*cosh([4])/4000.,[8]))+[10]", jp::xmin, jp::xmax);
    fdt->SetParameters(fus->GetParameter(0), fus->GetParameter(1),
		       fus->GetParameter(2), fus->GetParameter(3),
		       fus->GetParameter(4),
		       ftheory->GetParameter(0), ftheory->GetParameter(1),
		       ftheory->GetParameter(2), ftheory->GetParameter(3),
		       ftheory->GetParameter(5));
    fdt->SetParameter(10, 0.);

    offsetGraph(gtheory, offset[i]);
    if (ftheory) ratioGraph(g, ftheory);
    else ratioGraph(g, fus);
    offsetGraph(g, offset[i]);
    if (ftheory) ratioGraph(gmc, ftheory);
    else         ratioGraph(gmc, fus);
    offsetGraph(gmc, offset[i]);
    gmc->SetLineWidth(2);
    gmc->SetLineStyle(kDashed);

    // Produce uncertainty band around ansatz fit
    TH1D *hpl = (TH1D*)d2->Get("htot_pl"); assert(hpl);
    TH1D *hmn = (TH1D*)d2->Get("htot_mn"); assert(hmn);
    TF1 *f = new TF1(Form("f%d",i),"[0]", jp::xmin, jp::xmax);
    f->SetParameter(0, 1+offset[i]);

    TGraphErrors *gerr = (jp::centerOnAnsatz ?
			  makeBand(fdt, hpl,hmn, jp::xmin,g->GetX()[g->GetN()-1]) :
			  makeBand(g, hpl,hmn, jp::xmin,g->GetX()[g->GetN()-1]));
    fdt->SetParameter(10, offset[i]);
    offsetGraph(gerr, offset[i]);

    g->SetMarkerStyle(marker[i]);
    gmc->SetMarkerStyle(marker[i]);
    gmc->SetMarkerColor(kBlue);
    gmc->SetLineColor(kBlue+2);
    gtheory->SetLineColor(kRed+2);
    gtheory->SetLineWidth(2);
    gerr->SetFillStyle(1001);
    gerr->SetFillColor(kYellow+1);
    f->SetLineWidth(2);
    f->SetLineColor(kRed);

    //TLine *line = new TLine(xmin, offset[i]+off, bxmax, offset[i]+off);
    TLine *line = new TLine(jp::xmin, offset[i]+off, jp::xmax, offset[i]+off);
    if (i!=0) line->Draw();

    gerr->Draw(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    f->Draw("SAME");
    gmc->Draw("SAMEL");
    g->Draw("SAMEP");

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    leg->AddEntry(g,seta,"P");
    leg->AddEntry(gmc," ","");

    if (i==0) {
      if (htheory) leg2->AddEntry(f, "NLO", "L");
      else         leg2->AddEntry(f, "ansatz fit to data", "L");
      leg2->AddEntry(gmc,"Pythia D6T", "L");
      leg2->AddEntry(gerr, "Exp. uncertainty", "F");
    }
    assert(gmc->GetN());

    vnlo[i] = f;
    vgnlo[i] = gtheory;
    vgerr[i] = gerr;
    vmc2[i] = gmc;
    vdt[i] = g;
  } // for i

  // Get the tick marks on top of the bands
  h->Draw("AXIS SAME");
  leg->Draw();
  leg2->Draw();
  // Theory scale uncertainty band
  TLine *lup = new TLine();
  lup->SetLineColor(kRed); lup->SetLineStyle(kDotted);
  lup->DrawLineNDC(0.19,0.92,0.245,0.92);
  lup->DrawLineNDC(0.19,0.89,0.245,0.89);

  //if(_eps) c1->SaveAs("eps/summaryDataTheoryRatio.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryDataTheoryRatio.pdf");
  //if(_gif) c1->SaveAs("gif/summaryDataTheoryRatio.gif");
  //if(_png) c1->SaveAs("png/summaryDataTheoryRatio.png");

  /////////////////////////////////////////////////////
  // for PAS, cut at pT>20 GeV, color error band yellow
  TH1D *hpas = (TH1D*)h->DrawClone("AXIS");
  hpas->GetXaxis()->SetRangeUser(jp::xmin, jp::xmax);
  const int nylab = 12;
  const double ylabels[nylab] = {0, 1, 2, 1, 2, 1, 0, 1, 2, 3, 4, 5};
  for (int i = 0; i != nylab; ++i)
    t->DrawTextNDC(0.14, 0.135+0.0665*i, Form("%1.0f",ylabels[i]));
  for (unsigned int i = 0; i != etas.size(); ++i) {
    vgerr[i]->SetFillColor(kYellow+1);
    vgerr[i]->Draw(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
  }
  for (unsigned int i = 0; i != etas.size(); ++i) {

    //TLine *line = new TLine(jp::xminpas, offset[i]+off, bxmax, offset[i]+off);
    TLine *line = new TLine(jp::xminpas, offset[i]+off, jp::xmax, offset[i]+off);
    if (i!=0) line->Draw();

    vnlo[i]->Draw("SAME");
    TGraphErrors *gmc = vmc2[i];
    for (int j = 0; j != gmc->GetN(); ++j)
      gmc->SetPointError(j, 0., 0.);
    gmc->SetLineColor(kBlue);
    gmc->Draw("SAME L");
    vdt[i]->Draw("SAME P");
  }
  hpas->Draw("AXIS SAME");
  //cmsPrel(jp::lumi);
  leg->Draw();
  leg2->Draw();
  if (jp::centerOnAnsatz)
    ta->DrawTextNDC(leg2->GetX1NDC()+0.07, leg2->GetY1NDC()-0.04,//0.02,
		    "(centered on ansatz)");
  lup->DrawLineNDC(0.19,0.92,0.245,0.92);
  lup->DrawLineNDC(0.19,0.89,0.245,0.89);

  //if(_eps) c1->SaveAs("eps/summaryDataTheoryRatioPAS.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryDataTheoryRatioPAS.pdf");
  //if(_gif) c1->SaveAs("gif/summaryDataTheoryRatioPAS.gif");
  //if(_png) c1->SaveAs("png/summaryDataTheoryRatioPAS.png");

  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawDataTheoryRatio

void drawDataTheoryRatio3x2() {//bool ak7, bool v38x) {

  //bool _usejpt = true;
  bool _ak7 = (string(jp::algo)=="AK7");

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile("output-DATA-4c.root","READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  TFile *fin2 = new TFile("output-DATA-4.root","READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  // EXTRA
  //TFile *fin2010 = new TFile("../pfjet/output-DATA-4c.root","READ");
  /*
  TFile *fin2010 = new TFile("fastnlo/data2010.root");
  assert(fin2010 && !fin2010->IsZombie());
  assert(fin2010->cd("Standard"));
  TDirectory *din2010 = gDirectory;
  */
  //TFile *fin2010t = new TFile("../pfjet/output-DATA-5.root","READ");
  //assert(fin2010t && !fin2010t->IsZombie());
  //assert(fin2010t->cd("Standard"));
  //TDirectory *din2010t = gDirectory;
  //EXTRA

  // EXTRA2
  //TFile *fin2011 = new TFile("backup/oct26/output-DATA-4c.root","READ");
  //TFile *fin2011 = new TFile("fastnlo/data2011v2.root","READ");
  //TFile *fin2011 = new TFile("GRV23_AK5/output-DATA-4c.root","READ");
  /*
  TFile *fin2011 = new TFile("GRV23_AK5_42x_v4/output-DATA-4c.root","READ");
  assert(fin2011 && !fin2011->IsZombie());
  assert(fin2011->cd("Standard"));
  TDirectory *din2011 = gDirectory;
  */
  //TFile *fin2011t = new TFile("backup/oct26/output-DATA-5.root","READ");
  //assert(fin2011t && !fin2011t->IsZombie());
  //assert(fin2011t->cd("Standard"));
  //TDirectory *din2011t = gDirectory;
  //EXTRA

  TFile *fmc = new TFile("output-MC-2a.root","READ");
  assert(fmc && !fmc->IsZombie());
  assert(fmc->cd("Standard"));
  TDirectory *dinmc = gDirectory;

  TFile *fnlo = new TFile("output-DATA-5.root","READ");
  assert(fnlo && !fnlo->IsZombie());
  assert(fnlo->cd("Standard"));
  TDirectory *dnlo = gDirectory;

  /*
  TFile *fnlo11 = new TFile("GRV23_AK5_42x_v4/output-DATA-5.root","READ");
  assert(fnlo11 && !fnlo11->IsZombie());
  assert(fnlo11->cd("Standard"));
  TDirectory *dnlo11 = gDirectory;
  */

  TFile *fout = new TFile("comparisons/graphsDataNLORatio.root","RECREATE");
  assert(fout && !fout->IsZombie());

  curdir->cd();

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");
  //string ssak7 = (ak7 ? "_AK7" : "");
  //string ssv38x = (v38x ? "_38X" : "");
  //string ssid = ssak7 + ssv38x;
  const char *sid = "";//ssid.c_str();
  //bool _ak7 = (string(jp::algo)=="AK7");

  TCanvas *c1a = new TCanvas("c1a","c1a",1800,1200);
  c1a->SetTopMargin(0.10);
  c1a->Divide(3,2,-1,-1);

  TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
  c1->SetTopMargin(0.10);
  c1->Divide(3,2,-1,-1);

  TCanvas *c2a = new TCanvas("c2a","c2a",600,600);

  TCanvas *c2b = new TCanvas("c2b","c2b",1800,1200);
  c2b->SetTopMargin(0.10);
  c2b->Divide(3,2,-1,-1);

  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.0001);
  h->SetMaximum(1.7999);
  //h->GetYaxis()->SetTitle("Data / NLO#otimesNP (PDF4LHC)"); // cmsFinal
  h->GetYaxis()->SetTitle("Data / NLO#otimesNP (CT10)"); // cmsFinal
  if (jp::herapdf) h->GetYaxis()->SetTitle("Data / NLO#otimesNP (HERAPDF1.5)");
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  TH1D *h2 = (TH1D*)h->Clone("h2");

  // Set text scale fixed in pixels
  int npix = int(min(c1->GetWw(),c1->GetWh())*0.045*0.7);

  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(npix);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(npix);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(npix);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(npix);
  h->GetYaxis()->SetTitleOffset(2.5);
  h->GetXaxis()->SetTitleOffset(2.5);

  TLegend *leg = new TLegend(0.58,0.70,0.93,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  for (unsigned int iy = 0; iy != etas.size(); ++iy) {

    assert(din->cd(etas[iy].c_str()));
    TDirectory *d = gDirectory;

    assert(din2->cd(etas[iy].c_str()));
    TDirectory *d2 = gDirectory;

    /*
    assert(din2010->cd(etas[iy].c_str()));
    TDirectory *d2010 = gDirectory;
    //assert(din2010t->cd(etas[iy].c_str()));
    //TDirectory *d2010t = gDirectory;

    assert(din2011->cd(etas[iy].c_str()));
    TDirectory *d2011 = gDirectory;
    */
    //assert(din2011t->cd(etas[iy].c_str()));
    //TDirectory *d2011t = gDirectory;

    assert(dnlo->cd(etas[iy].c_str()));
    TDirectory *dt = gDirectory;

    /*
    assert(dnlo11->cd(etas[iy].c_str()));
    TDirectory *dt11 = gDirectory;
    */

    assert(dinmc->cd(etas[iy].c_str()));
    TDirectory *dmc = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[iy].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TH1D *hg = (TH1D*)d->Get("hcorrpt"); assert(hg); // for global fit
    TGraphErrors *g = (TGraphErrors*)d->Get("gcorrpt"); assert(g);
    /*
    TGraphErrors *g2010 = (TGraphErrors*)d2010->Get("gcorrpt_38x"); assert(g2010);
    TGraphErrors *g2011 = (TGraphErrors*)d2011->Get("gcorrpt"); assert(g2011);
    */
    //TGraphErrors *g2010 = (TGraphErrors*)d2010->Get("gcorrpt"); assert(g2010);
    /*
    TGraphErrors *gc = (TGraphErrors*)d->Get("gcorrpt_ak5calo");
    TGraphErrors *gj = (TGraphErrors*)d->Get("gcorrpt_ak5jpt");
    if (v38x && !ak7) {
      g = (TGraphErrors*)d->Get("gcorrpt_38x"); assert(g);
      gc = (TGraphErrors*)d->Get("gcorrpt_ak5calo_38x");
      gj = (TGraphErrors*)d->Get("gcorrpt_ak5jpt_38x");
    }
    if (ak7 && !v38x) {
      g = (TGraphErrors*)d->Get("gcorrpt_ak7pf"); assert(g);
      gc = (TGraphErrors*)d->Get("gcorrpt_ak7calo");
      gj = (TGraphErrors*)d->Get("gcorrpt_ak7jpt");
    }
    if (ak7 && v38x) {
      g = (TGraphErrors*)d->Get("gcorrpt_ak7pf_38x"); assert(g);
      gc = (TGraphErrors*)d->Get("gcorrpt_ak7calo_38x");
      gj = (TGraphErrors*)d->Get("gcorrpt_ak7jpt_38x");
    }
    // turn-off calo
    gc = 0;
    */
    // Load theory
    TF1 *fus = (TF1*)d->Get("fus"); assert(fus);
    TH1D *htheory = (TH1D*)dt->Get("hnlo"); assert(htheory);
    /*
    TH1D *hdataref = (TH1D*)dt->Get(_ak7 ? "hdata2011" : "hdata2010");
    TH1D *hdataref_up = (TH1D*)dt->Get(_ak7 ? "hdata2011_up" : "hdata2010_up");
    TH1D *hdataref_dw = (TH1D*)dt->Get(_ak7 ? "hdata2011_dw" : "hdata2010_dw");
    */
    TH1D *hdataref = (TH1D*)dt->Get("hdata2012");
    TH1D *hdataref_up = 0;//(TH1D*)dt->Get("hdata2012");
    TH1D *hdataref_dw = 0;//(TH1D*)dt->Get("hdata2012");

    //assert(hdataref);
    //TH1D *htheory2010 = (TH1D*)dt->Get("hnlo2010"); assert(htheory2010);
    /*
    TH1D *htheory2010 = (TH1D*)dt11->Get("hnlo2010"); assert(htheory2010);
    //TH1D *htheory2011 = (TH1D*)dt->Get("hnlo2011"); assert(htheory2011);
    TH1D *htheory2011 = (TH1D*)dt11->Get("hnlo2011"); assert(htheory2011);
    */
    if (jp::herapdf) {
      htheory = (TH1D*)dt->Get("hnlo_hera"); assert(htheory); // HERAPDF
    }
    TGraphErrors *gtheory = (TGraphErrors*)dt->Get("gnlofit"); assert(gtheory);
    TH1D *hthsysup = (TH1D*)dt->Get("hsysup"); assert(hthsysup);
    TH1D *hthsysdw = (TH1D*)dt->Get("hsysdw"); assert(hthsysdw);
    /*
    if (ak7) {
      fus = (TF1*)d->Get("fus"); assert(fus);
      htheory = (TH1D*)dt->Get("hnlo_ak7"); assert(htheory);
      gtheory = (TGraphErrors*)dt->Get("gnlofit_ak7"); assert(gtheory);
      hthsysup = (TH1D*)dt->Get("hsysup"); assert(hthsysup);
      hthsysdw = (TH1D*)dt->Get("hsysdw"); assert(hthsysdw);
    }
    */

    // Generator level MC spectrum
    TH1D *hmc = (TH1D*)dmc->Get("mc/hpt_g0"); //assert(hmc);

    // Ratio of data and theory
    hg->Divide(htheory);
    ratioGraph(g, htheory);
    /*
    ratioGraph(g2010, htheory2010);
    ratioGraph(g2011, htheory2011);
    */
    if (hdataref) hdataref->Divide(htheory);
    if (hdataref_up) hdataref_up->Divide(htheory);
    if (hdataref_dw) hdataref_dw->Divide(htheory);
    //if (gc) ratioGraph(gc, htheory);
    //if (gj) ratioGraph(gj, htheory);

    // Remove low pT CaloJet and JPTJet points
    //if (gc)
    //for (int i = gc->GetN()-1; i != -1; --i)
    //if (gc->GetX()[i]<30.) gc->RemovePoint(i);
    //if (gj)
    //for (int i = gj->GetN()-1; i != -1; --i)
    //if (gj->GetX()[i]<20.) gj->RemovePoint(i);

    // Calculate ansatz at bin center
    TGraphErrors *gus = new TGraphErrors(0);
    TF1 *fusnew = new TF1(Form("fusnew_%d",iy),
			  "[0]*exp([1]/x)*pow(x,[2])"
			  "*pow(1-x*cosh([4])/4000.,[3])", 1., jp::xmax);
    for (int i = 0; i != fus->GetNpar(); ++i) {
      fusnew->SetParameter(i, fus->GetParameter(i));
    }
    for (int i = 0; i != g->GetN(); ++i) {
      double x = g->GetX()[i];
      double y = fusnew->Eval(x);
      if (x >= jp::xmin) {
	int n = gus->GetN();
	gus->SetPoint(n, x, y);
      }
    }
    ratioGraph(gus, htheory);
    gus->SetLineColor(kGray+2);

    // Produce experimental uncertainty band around ansatz fit, theory or data
    TGraphErrors *gth = new TGraphErrors(htheory);
    ratioGraph(gth, htheory);
    TH1D *hpl = (TH1D*)d2->Get("htot_pl"); assert(hpl);
    TH1D *hmn = (TH1D*)d2->Get("htot_mn"); assert(hmn);
    //if (v38x) {
    //hpl = (TH1D*)d2->Get("htot_38x_pl"); assert(hpl);
    //hmn = (TH1D*)d2->Get("htot_38x_mn"); assert(hmn);
    //}
    TGraphErrors *gerr = (jp::centerOnAnsatz ?
			  makeBand(gus, hpl, hmn, jp::xmin, jp::xmax) :
			  jp::centerOnTheory ?
			  makeBand(gth, hpl, hmn, jp::xmin, jp::xmax, true) :
			  makeBand(g, hpl, hmn, jp::xmin, jp::xmax));
    // Fix for new missing low pT theory
    for (int i = gerr->GetN()-1; i != -1; --i) {
      if (fabs(gerr->GetY()[i]-1)>0.5) gerr->RemovePoint(i);
    }

    // Produce theory uncertainty band
    TGraph *gthsysup = new TGraph(gus->GetN());
    TGraph *gthsysdw = new TGraph(gus->GetN());
    for (int i = 0; i != gus->GetN(); ++i) {
      double x = gus->GetX()[i];
      int jup = hthsysup->FindBin(x);
      int jdw = hthsysdw->FindBin(x);
      gthsysup->SetPoint(i, x, 1+hthsysup->GetBinContent(jup));
      gthsysdw->SetPoint(i, x, 1+hthsysdw->GetBinContent(jdw));
      hthsysup->SetBinContent(jup, 1+hthsysup->GetBinContent(jup));
      hthsysdw->SetBinContent(jdw, 1+hthsysdw->GetBinContent(jdw));
    }
    gthsysup->SetLineColor(kRed);
    gthsysdw->SetLineColor(kRed);
    hthsysup->SetLineColor(kRed);
    hthsysdw->SetLineColor(kRed);

    // Produce MC prediction
    TGraph *gmc = new TGraph(gus->GetN());
    for (int i = 0; i != gus->GetN() && hmc; ++i) {
      double x = gus->GetX()[i];
      gmc->SetPoint(i, x, hmc->GetBinContent(hmc->FindBin(x))
		    / htheory->GetBinContent(htheory->FindBin(x)));
    }
    gmc->SetLineWidth(2);
    gmc->SetLineColor(kBlue);
    gmc->SetLineStyle(kDashed);

    // 3x2 format for AN
    c1a->cd(iy+1);
    gPad->SetLogx();
    h->GetXaxis()->SetRangeUser(jp::xmin, jp::xmax);
    h->DrawClone("AXIS");

    //TObject *gerr_1a = gerr->DrawClone(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    TObject *gerr_1a = gerr->DrawClone("SAME E3");
    if (jp::centerOnAnsatz) gus->DrawClone("SAMEL");
    TLine *l = new TLine(jp::xmin, 1., jp::xmax, 1.);
    l->DrawClone();
    gthsysup->DrawClone("SAME L");
    gthsysdw->DrawClone("SAME L");
    if (jp::plotPythia) gmc->DrawClone("SAMEL");

    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(0.7);
    g->DrawClone("SAME PZ");

    //TLegend *leg = new TLegend(30, _plotPythia ? 1.2 : 1.3, 150,1.75,"","br");
    TLegend *leg = new TLegend(60, jp::plotPythia ? 1.1 : 1.2, 300, 1.8,"","br");
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->SetTextFont(43);
    leg->SetTextSize(npix);
    //leg->AddEntry(g, "Data / NLO#otimesNP theory", "P");
    leg->AddEntry(g, "Data", "P");
    //leg->AddEntry(g2010, "2010 Data / theory", "P"); // EXTRA
    if (hdataref) leg->AddEntry(hdataref, "Sanmay, 8 TeV", "P");
				// (_ak7 ? "QCD-11-004" : "QCD-10-010"), "P");
    leg->AddEntry(gthsysup, "Theory uncertainty","L");
    leg->AddEntry(gerr_1a, "Exp. uncertainty","F");
    if (jp::plotPythia) leg->AddEntry(gmc, "Pythia D6T","L");
    leg->Draw();

    TLegend *leg2010 = new TLegend(30, 1.3, 150, 1.75, "", "br");
    leg2010->SetBorderSize(0);
    leg2010->SetFillStyle(kNone);
    leg2010->SetTextFont(43);
    leg2010->SetTextSize(npix);
    //leg2010->AddEntry(g2011, "2011 Data(v1) / theory", "P"); // EXTRA
    /*
    leg2010->AddEntry(g2011, "2011(42X) Data / theory", "P"); // EXTRA
    leg2010->AddEntry(g2010, "2010 Data / theory", "P"); // EXTRA
    */
    //leg2010->AddEntry(hdataref, (_ak7 ? "QCD-11-004" : "QCD-10-010"), "P");
    leg2010->Draw();

    //TLatex *t = new TLatex(40, 0.15, etamin==0 ?
    //			   Form("|y| < %1.1f",etamax) :
    //			   Form("%1.1f #leq |y| < %1.1f",etamin,etamax));
    // Vivian
    TLatex *t = new TLatex(100, 0.65, etamin==0 ?
    			   Form("|y| < %1.1f",etamax) :
    			   Form("%1.1f #leq |y| < %1.1f",etamin,etamax));
    t->SetTextFont(43);
    t->SetTextSize(npix);
    t->DrawClone();

    h->DrawClone("AXIS SAME");


    // Single bins with different jet types
    c2a->cd();
    gPad->SetLogx();
    TH1D *h2a = (TH1D*)h2->DrawClone("AXIS");
    h2a->GetXaxis()->SetRangeUser(37, 1999.);//xmax);
    h2a->SetMinimum(0.5);
    h2a->SetMaximum(1.5);

    //gerr->DrawClone(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    gerr->DrawClone("SAME E3");
    if (jp::centerOnAnsatz) gus->DrawClone("SAMEL");
    l->DrawLine(37,1,jp::xmax,1);
    gthsysup->DrawClone("SAME L");
    gthsysdw->DrawClone("SAME L");

    /*
    if (g2010) {
      g2010->SetMarkerStyle(kOpenCircle);
      g2010->SetMarkerColor(kBlue);
      g2010->SetLineColor(kBlue);
      g2010->DrawClone("SAMEPz");
    }
    if (g2011) {
      g2011->SetMarkerStyle(kOpenSquare);
      g2011->SetMarkerColor(kRed);
      g2011->SetLineColor(kRed);
      //g2011->DrawClone("SAMEPz");
    }
    */
    //if (gc) {

      // Clean out bugged low pT end temporarily
      //for (int i = gc->GetN(); i != -1; --i) {
      //if (gc->GetX()[i] < 150.) gc->RemovePoint(i);
      //}

      //gc->SetMarkerStyle(kOpenSquare);//kOpenCircle);
      //gc->SetMarkerColor(kRed);//kBlue);
      //gc->SetLineColor(kRed);//kBlue);
      //gc->DrawClone("SAMEPz");
    //}
    //if (gj && _usejpt) {
    //gj->SetMarkerStyle(kOpenSquare);
    //gj->SetMarkerColor(kRed);
    //gj->SetLineColor(kRed);
    //gj->DrawClone("SAMEPZ");
    //}

    g->SetFillStyle(1001);
    g->SetFillColor(kYellow+1);
    g->SetMarkerSize(1.);
    g->DrawClone("SAME PZ");

    t->SetTextFont(42);
    t->SetTextSize(0.045);
    TLatex *t2a = (TLatex*)t->DrawClone();
    t2a->SetX(700.);
    t2a->SetY(1.40);

    //TLegend *leg2a = new TLegend(0.25,_usejpt?0.75:0.81,0.35,0.93,"","brNDC");
    //TLegend *leg2a = new TLegend(0.20,_usejpt?0.75:0.81,0.60,0.93,"","brNDC");
    TLegend *leg2a = new TLegend(0.20,0.81,0.60,0.93,"","brNDC");
    leg2a->SetBorderSize(0);
    leg2a->SetFillStyle(kNone);
    leg2a->SetTextSize(0.045);
    //leg2a->AddEntry(g,"PF","P");
    leg2a->AddEntry(g,"2011 PF (4.7/fb)","PF");
    //if (gc) leg2a->AddEntry(gc,"2011 Calo (4.7/fb)","P");
    //if (g2011) leg2a->AddEntry(g2011,"2011 (v2)","P");
    /*if (g2010) leg2a->AddEntry(g2010,"2010 PF (36/pb)","P");*/
    //if (gj && _usejpt) leg2a->AddEntry(gj,"JPT","P");
    leg2a->Draw();

    TLegend *leg2b = new TLegend(0.40,0.75,0.70,0.93,"","brNDC");
    leg2b->SetBorderSize(0);
    leg2b->SetFillStyle(kNone);
    leg2b->SetTextSize(0.045);
    leg2b->AddEntry(g, "Data / NLO#otimesNP theory","PF"); // cmsFinal
    leg2b->AddEntry(gthsysup, "Theory uncertainty","L");
    leg2b->AddEntry(gerr, "Exp. uncertainty","F");
    //leg2b->Draw();

    //cmsPrel(jp::lumi);
    h2->DrawClone("AXIS SAME");

    //if(_eps) c2a->SaveAs(Form("eps/compareReco_Rap%d%s.eps",iy,sid));
    if(jp::pdf) c2a->SaveAs(Form("pdf/compareReco_Rap%d%s.pdf",iy,sid));

    c2b->cd(iy+1);
    gPad->SetLogx();
    //h2->GetXaxis()->SetRangeUser(jp::xmin,xmax);
    TH2D *h2b = (TH2D*)h2->Clone(Form("%s_2b",h2->GetName()));
    double x7min = 97; double x7max = 2116;
    h2b->GetXaxis()->SetRangeUser(x7min, x7max);
    h2b->GetYaxis()->SetRangeUser(0.5+0.0001,1.5-0.0001);
    h2b->GetXaxis()->SetTitle(iy==5 ? "p_{T} (GeV)" : "");
    h2b->GetYaxis()->SetTitle(iy==0 ? "Data / NLO#otimesNP (CT10)" : "");
    h2b->DrawClone("AXIS");

    gerr->DrawClone(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    if (jp::centerOnAnsatz) gus->DrawClone("SAMEL");
    //l->Draw();
    l->DrawLine(x7min, 1, x7max, 1);
    gthsysup->DrawClone("SAME L");
    gthsysdw->DrawClone("SAME L");

    /*
    if (g2010) g2010->SetMarkerSize(1.3);//0.7);
    if (g2010) g2010->DrawClone("SAMEPz");
    */
    //if (gc) gc->SetMarkerSize(1.3);//0.7);
    //if (gc) gc->DrawClone("SAMEPz");
    //if (gj) gj->SetMarkerSize(0.7);
    //if (gj && _usejpt) gj->DrawClone("SAMEPz");
    g->SetMarkerSize(1.3);//0.7);
    g->DrawClone("SAME Pz");

    TLatex *t2b = new TLatex(110, 1.40, etamin==0 ?
			     Form(" |y| < %1.1f",etamax) :
			     Form("%1.1f < |y| < %1.1f",etamin,etamax));
    t2b->SetTextFont(43);
    t2b->SetTextSize(npix);
    t2b->Draw();

    TLegend *legc2b = new TLegend(400,1.30,2000,1.47,"","br");
    legc2b->SetBorderSize(0);
    legc2b->SetFillStyle(kNone);
    //legc2b->SetTextSize(0.045);
    legc2b->SetTextFont(43);
    legc2b->SetTextSize(npix);
    legc2b->AddEntry(g,"2011 (4.7/fb)","PF");
    /* if (g2010) legc2b->AddEntry(g2010,"2010 (36/pb)","P"); */
    if (iy==2) legc2b->Draw();

    //t->DrawClone();
    //if (iy==2) leg2a->Draw();
    //leg2b->Draw();
    //h2->DrawClone("AXIS SAME");

    gPad->RedrawAxis();

    fout->cd();
    g->Write(Form("DataOverNLO_Rap%1.2g-%1.2g%s",0.5*iy,0.5*(iy+1),sid));
    hg->Write(Form("HistDataOverNLO_Rap%1.2g-%1.2g%s",0.5*iy,0.5*(iy+1),sid));
    curdir->cd();


    // 3x2 format for PRL
    c1->cd(iy+1);
    gPad->SetLogx();
    //h->GetXaxis()->SetRangeUser(jp::xminpas, xmax);
    h->GetXaxis()->SetRangeUser(74., jp::xmax); // Vivian
    h->SetMinimum(0.5+0.001);//0.6001); // Vivian
    h->SetMaximum(1.5-0.001);//1.3999); // Vivian
    h->DrawClone("AXIS");
    h->SetMinimum(0.0001);
    h->SetMaximum(1.7999);

    rangeGraph(gerr,jp::xminpas,jp::xmax);
    rangeGraph(gus,jp::xminpas,jp::xmax);
    rangeGraph(gthsysup,jp::xminpas,jp::xmax);
    rangeGraph(gthsysdw,jp::xminpas,jp::xmax);
    hthsysup->GetXaxis()->SetRangeUser(jp::xminpas,g->GetX()[g->GetN()-1]);
    hthsysdw->GetXaxis()->SetRangeUser(jp::xminpas,g->GetX()[g->GetN()-1]);
    rangeGraph(gmc,jp::xminpas,jp::xmax);

    gerr->SetFillColor(kYellow+1);
    gus->SetLineColor(kYellow+2);
    //gerr->DrawClone(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    //gerr->DrawClone("SAME E3");
    gerr->DrawClone(jp::centerOnTheory ? "SAME E2" : "SAME E3");
    if (jp::centerOnAnsatz) gus->DrawClone("SAMEL");
    //l->DrawLine(jp::xminpas, 1., xmax, 1.);
    l->DrawLine(74, 1., jp::xmax, 1.);
    //gthsysup->DrawClone("SAME L");
    //gthsysdw->DrawClone("SAME L");
    hthsysup->DrawClone("SAME HIST][");
    hthsysdw->DrawClone("SAME HIST][");
    if (jp::plotPythia) gmc->DrawClone("SAME L");

    // Vivian
    l->SetLineStyle(kDotted);
    l->SetLineColor(kGray+1);
    /*
    l->DrawLine(114, 0.6, 114, 1.4); // 80
    l->DrawLine(220, 0.6, 220, 1.4); // 140
    l->DrawLine(300, 0.6, 300, 1.4); // 200
    l->DrawLine(362, 0.6, 362, 1.4); // 260
    l->DrawLine(430, 0.6, 430, 1.4); // 320
    l->DrawLine(507, 0.6, 507, 1.4); // 400
    */
    l->DrawLine(133, 0.6, 133, 1.4); // 80
    l->DrawLine(220, 0.6, 220, 1.4); // 140
    l->DrawLine(300, 0.6, 300, 1.4); // 200
    l->DrawLine(395, 0.6, 395, 1.4); // 260
    l->DrawLine(507, 0.6, 507, 1.4); // 320, was 468
    //l->DrawLine(592, 0.6, 592, 1.4); // 400
    l->SetLineColor(kBlack);
    l->SetLineStyle(kSolid);

    if (hdataref) {
      hdataref->SetMarkerStyle(kOpenCircle);
      hdataref->SetMarkerSize(1.0);//0.7);
      hdataref->SetMarkerColor(kRed);
      hdataref->SetLineColor(kRed);
      hdataref->DrawClone("SAME Pz");
    }
    if (hdataref_up) {
      TGraph *g = new TGraph(hdataref_up);
      for (int i = g->GetN()-1; i != -1; --i)
	if (g->GetY()[i]<0.1) g->RemovePoint(i);
      //hdataref_up->SetLineColor(kRed);
      //hdataref_up->SetLineStyle(kDotted);
      //hdataref_up->DrawClone("SAME HISTX");
      g->SetLineColor(kRed);
      g->SetLineStyle(kDotted);
      g->Draw("SAMEL");
    }
    if (hdataref_dw) {
      TGraph *g = new TGraph(hdataref_dw);
      for (int i = g->GetN()-1; i != -1; --i)
	if (g->GetY()[i]<0.1) g->RemovePoint(i);
      //hdataref_dw->SetLineColor(kRed);
      //hdataref_dw->SetLineStyle(kDotted);
      //hdataref_dw->DrawClone("SAME HISTX");
      g->SetLineColor(kRed);
      g->SetLineStyle(kDotted);
      g->Draw("SAMEL");
    }

    /*
    g2010->SetMarkerStyle(kOpenCircle);
    g2010->SetMarkerSize(1.0);//0.7);
    g2010->SetMarkerColor(kBlue);
    g2010->SetLineColor(kBlue);
    g2010->DrawClone("SAME Pz"); // EXTRA

    g2011->SetMarkerStyle(kOpenDiamond);
    g2011->SetMarkerSize(1.3);
    g2011->SetMarkerColor(kRed);
    g2011->SetLineColor(kRed);
    g2011->DrawClone("SAME Pz"); // EXTRA
    */

    g->SetMarkerSize(1.0);//0.7);
    g->DrawClone("PZ");

    //TBox *box = new TBox(xmin, 0.4, 56., 1.3);
    //TBox *box = new TBox(jp::xminpas, 0.4, _ak7 ? 64 : 56., 1.5);
    TBox *box = new TBox(jp::xminpas, 0.4, _ak7 ? 74 : 64., 1.8);
    box->SetFillStyle(kCrossHatch);
    box->SetLineColor(kGray);
    box->SetFillColor(kGray);
    //box->Draw("SAME");
    TLatex *tbox = new TLatex();
    tbox->SetTextFont(43);
    tbox->SetTextSize(npix);
    tbox->SetTextColor(kGray+1);
    //tbox->DrawLatex(jp::xminpas*1.05, 1.2, "trigger");

    t->SetTextFont(43);
    t->SetTextSize(npix);
    t->DrawClone();

    if (iy==1) leg2010->Draw();
    if (iy==2) leg->Draw();
    h->DrawClone("AXIS SAME");


    // Single bins with different jet types
    c2->cd();
    gPad->SetLogx();
    h2->GetXaxis()->SetRangeUser(jp::xminpas, jp::xmax);
    h2->DrawClone("AXIS");

    //gerr->DrawClone(jp::centerOnAnsatz ? "SAME E3" : "SAME E3");
    gerr->DrawClone("SAME E3");
    if (jp::centerOnAnsatz) gus->DrawClone("SAMEL");
    l->DrawLine(jp::xminpas, 1., jp::xmax, 1.);
    gthsysup->DrawClone("SAME L");
    gthsysdw->DrawClone("SAME L");

    //if (gc) gc->DrawClone("SAMEPZ");
    //if (gj && _usejpt) gj->DrawClone("SAMEPZ");
    /*
    g2010->Draw("SAMEPZ");
    g2011->Draw("SAMEPZ");
    */
    g->SetMarkerSize(1.0);
    g->DrawClone("SAME PZ");


    t->SetTextFont(42);
    t->SetTextSize(0.045);
    t->DrawClone();

    //leg2a->Draw();
    //leg2b->Draw();
    leg->Draw();
    TLegend *leg2010b = (TLegend*)leg2010->Clone();
    leg2010b->SetY1(0.25);
    leg2010b->SetY2(0.55);
    leg2010b->Draw();

    h2->DrawClone("AXIS SAME");
    //if (!ak7 && !v38x && iy==0 && _final) cmsFinal(jp::lumi);
    //else
    //cmsPrel(jp::lumi);


    //if(_eps) c2->SaveAs(Form("eps/compareRecoPAS_Rap%d%s.eps",iy,sid));
    if(jp::pdf) c2->SaveAs(Form("pdf/compareRecoPAS_Rap%d%s.pdf",iy,sid));
  } // for iy

  c1a->cd(0);

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);

  latex->SetTextAlign(11); // align left
  latex->DrawLatex(0.06,0.96,Form("CMS preliminary, %.2g pb^{-1}",jp::lumi));
  latex->SetTextAlign(21); // align middle
  //latex->DrawLatex(0.55,0.96,ak7 ? "Anti-k_{T} R=0.7" : "Anti-k_{T} R=0.5");
  latex->DrawLatex(0.55,0.96,_ak7 ? "Anti-k_{T} R=0.7" : "Anti-k_{T} R=0.5");
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.98,0.96,"#sqrt{s} = 7 TeV");

  //if(_eps) c1a->SaveAs(Form("eps/summaryDataTheoryRatio3x2%s.eps",sid));
  if(jp::pdf) c1a->SaveAs(Form("pdf/summaryDataTheoryRatio3x2%s.pdf",sid));

  c1->cd();

  latex->SetTextAlign(21); // align middle
  latex->DrawLatex(0.55,0.96,_ak7 ? "Anti-k_{T} R=0.7" : "Anti-k_{T} R=0.5");
  //if (!ak7 && !v38x && _final) cmsFinal(jp::lumi, true);
  //else
  //cmsPrel(jp::lumi, true);

  const char *t = string(jp::algo).c_str();
  //if(_eps) c1->SaveAs(Form("eps/summaryDataTheoryRatioPRL3x2%s%s.eps",sid,t));
  if(jp::pdf) c1->SaveAs(Form("pdf/summaryDataTheoryRatioPRL3x2%s%s.pdf",sid,t));

  c2b->cd(0);

  latex->SetTextAlign(11); // align left
  //latex->DrawLatex(0.06, 0.96, Form("CMS preliminary, %.2g fb^{-1}",jp::lumi/1e3));
  latex->DrawLatex(0.06, 0.96, Form("CMS preliminary"));//,jp::lumi/1e3));
  latex->SetTextAlign(21); // align middle
  //latex->DrawLatex(0.55, 0.96, ak7 ? "Anti-k_{T} R=0.7" : "Anti-k_{T} R=0.5");
  //latex->DrawLatex(0.55, 0.96, string(jp::algo)=="AK7" ?
  latex->DrawLatex(0.52, 0.96, _ak7 ?
		   "Anti-k_{T} R = 0.7" : "Anti-k_{T} R = 0.5");
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.98, 0.96, "#sqrt{s} = 7 TeV");

  //if(_eps) c2b->SaveAs(Form("eps/compareReco_AllRap%s.eps",sid));
  if(jp::pdf) c2b->SaveAs(Form("pdf/compareReco_AllRap%s.pdf",sid));

  fout->Write();
  fout->Close();
} // drawDataTheoryRatio3x2


void drawTheoryDataRatio(bool ak7, bool v38x) {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile("output-DATA-4c.root","READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  TFile *fin2 = new TFile("output-DATA-4.root","READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  TFile *fnlo = new TFile("output-DATA-5.root","READ");
  assert(fnlo && !fnlo->IsZombie());
  assert(fnlo->cd("Standard"));
  TDirectory *dnlo = gDirectory;

  curdir->cd();

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");
  string ssak7 = (ak7 ? "_AK7" : "");
  string ssv38x = (v38x ? "_38X" : "");
  string ssid = ssak7 + ssv38x;
  const char *sid = ssid.c_str();

  TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
  c1->SetTopMargin(0.10);
  c1->Divide(3,2,-1,-1);

  const double tscale = 1.0;
  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.8001*tscale);
  h->SetMaximum(1.6999*tscale);
  h->GetYaxis()->SetTitle("NLO#otimesNP / Data");
  if (tscale!=1)
    h->GetYaxis()->SetTitle(Form("NLO#otimesNP #times %1.2f / Data",tscale));
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  // Set text scale fixed in pixels
  int npix = int(min(c1->GetWw(),c1->GetWh())*0.045*0.7);

  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(npix);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(npix);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(npix);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(npix);
  h->GetYaxis()->SetTitleOffset(2.5);
  h->GetXaxis()->SetTitleOffset(2.5);

  TLegend *leg = new TLegend(0.38,0.57,0.73,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  TLegend *leg2 = new TLegend(0.20,0.75,0.55,0.93,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);

  for (unsigned int iy = 0; iy != etas.size(); ++iy) {

    assert(din->cd(etas[iy].c_str()));
    TDirectory *d = gDirectory;

    assert(din2->cd(etas[iy].c_str()));
    TDirectory *d2 = gDirectory;

    assert(dnlo->cd(etas[iy].c_str()));
    TDirectory *dt = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[iy].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TH1D *hdt = (TH1D*)d->Get("hcorrpt"); assert(hdt);
    if (v38x && !ak7) {
      hdt = (TH1D*)d->Get("hcorrpt_38x"); assert(hdt);
    }
    if (ak7 && !v38x) {
      hdt = (TH1D*)d->Get("hcorrpt_ak7pf"); assert(hdt);
    }
    if (ak7 && v38x) {
      hdt = (TH1D*)d->Get("hcorrpt_ak7pf_38x"); assert(hdt);
    }

    // Load theory
    TH1D *hnlo = (TH1D*)dt->Get("hnlo"); assert(hnlo);
    TH1D *hnlo_ct66 = (TH1D*)dt->Get("hnlo_ct66"); assert(hnlo_ct66);
    TH1D *hnlo_ct10 = (TH1D*)dt->Get("hnlo_ct10"); assert(hnlo_ct10);
    TH1D *hnlo_mstw = (TH1D*)dt->Get("hnlo_mstw"); assert(hnlo_mstw);
    TH1D *hnlo_nnpdf = (TH1D*)dt->Get("hnlo_nnpdf"); assert(hnlo_nnpdf);
    TH1D *hnlo_hera = (TH1D*)dt->Get("hnlo_hera"); assert(hnlo_hera);
    if (ak7) {
      hnlo_ct10 = (TH1D*)dt->Get("hnlo_ak7"); assert(hnlo_ct10);
    }

    // Ratio of data and theory
    TH1D *hdt0 = (TH1D*)hdt->Clone("hdt0");
    for (int i = 0; i != hdt0->GetNbinsX()+1; ++i) {
      hdt0->SetBinContent(i, hdt->GetBinContent(i)/tscale);
      hdt0->SetBinError(i, 0.);
    }
    hnlo->Divide(hdt0);
    hnlo_ct66->Divide(hdt0);
    hnlo_ct10->Divide(hdt0);
    hnlo_mstw->Divide(hdt0);
    hnlo_nnpdf->Divide(hdt0);
    hnlo_hera->Divide(hdt0);
    delete hdt0;

    // Produce experimental uncertainty band around data
    TH1D *hpl = (TH1D*)d2->Get("htot_pl"); assert(hpl);
    TH1D *hmn = (TH1D*)d2->Get("htot_mn"); assert(hmn);
    if (v38x) {
      hpl = (TH1D*)d2->Get("htot_38x_pl"); assert(hpl);
      hmn = (TH1D*)d2->Get("htot_38x_mn"); assert(hmn);
    }
    TH1D *herr = (TH1D*)hpl->Clone("herr");
    for (int i = 1; i != herr->GetNbinsX()+1; ++i) {
      herr->SetBinContent(i, 1 + 0.5*(hpl->GetBinContent(i)
				      + hmn->GetBinContent(i)));
      herr->SetBinError(i, 0.5*(hpl->GetBinContent(i)-hmn->GetBinContent(i)));
    }
    // Produce theoretical uncertainty band around data
    TH1D *hthsysup = (TH1D*)dt->Get("hsysup"); assert(hthsysup);
    TH1D *hthsysdw = (TH1D*)dt->Get("hsysdw"); assert(hthsysdw);
    TH1D *hthsys = (TH1D*)hthsysup->Clone("hthsys");
    TH1D *hthup = (TH1D*)hthsysup->Clone("hthup");
    TH1D *hthdw = (TH1D*)hthsysdw->Clone("hthdw");
    for (int i = 1; i != hthsys->GetNbinsX()+1; ++i) {
      hthsys->SetBinContent(i, 1 - 0.5*(hthsysup->GetBinContent(i)
					+ hthsysdw->GetBinContent(i)));
      hthsys->SetBinError(i, 0.5*(hthsysup->GetBinContent(i)
				  - hthsysdw->GetBinContent(i)));
      hthup->SetBinContent(i, 1-hthsysup->GetBinContent(i));
      hthdw->SetBinContent(i, 1-hthsysdw->GetBinContent(i));
    }
    // Produce combined systematics band
    TH1D *hsys = (TH1D*)hpl->Clone("hsys");
    for (int i = 1; i != hsys->GetNbinsX()+1; ++i) {
      int j = hthsys->FindBin(hsys->GetBinCenter(i));
      double eup = tools::oplus(hpl->GetBinContent(i),
				hthsysdw->GetBinContent(j));
      double edw = -tools::oplus(hmn->GetBinContent(i),
				 hthsysdw->GetBinContent(j));
      hsys->SetBinContent(i, 1 + 0.5*(eup + edw));
      hsys->SetBinError(i, 0.5*(eup - edw));
    }

    // 3x2 format (support for PRL)
    c1->cd(iy+1);
    gPad->SetLogx();
    h->GetXaxis()->SetRangeUser(jp::xminpas, jp::xmax);
    h->DrawClone("AXIS");

    hsys->SetFillStyle(1001);
    hsys->SetFillColor(kGray+1);
    hsys->Draw("SAME E3");
    herr->SetFillStyle(1001);
    herr->SetFillColor(kYellow);
    herr->Draw("SAME E3");
    hthsys->SetLineWidth(1);
    hthsys->SetLineColor(kRed+2);
    hthsys->SetFillColor(kRed+2);
    hthsys->SetFillStyle(3003);
    hthup->SetLineColor(kRed);
    hthup->Draw("SAME");
    hthdw->SetLineColor(kRed);
    hthdw->Draw("SAME");

    if (iy==0) leg2->AddEntry(herr,"Exp. uncertainty","F");
    //if (iy==0) leg2->AddEntry(hthup,"PDF4LHC+NP uncertainty","L");
    if (iy==0) leg2->AddEntry(hthup,"CT10+NP uncertainty","L");
    if (iy==0) leg2->AddEntry(hsys,"Exp. #oplus Th. uncertainty","F");

    TLine *l = new TLine();
    l->DrawLine(jp::xminpas, 1., jp::xmax, 1.);

    hnlo_hera->SetMarkerSize(0.5);
    hnlo_hera->SetMarkerStyle(kFullCircle);
    hnlo_hera->SetMarkerColor(kRed);
    hnlo_hera->Draw("SAME P");
    if (iy==0) leg->AddEntry(hnlo_hera, "HERAPDF1.5", "P");
    //
    hnlo_ct66->SetMarkerSize(0.5);
    hnlo_ct66->SetMarkerStyle(kFullSquare);
    hnlo_ct66->SetMarkerColor(kBlue);
    hnlo_ct66->Draw("SAME P");
    if (iy==0) leg->AddEntry(hnlo_ct66, "CTEQ 6.6", "P");
    //
    hnlo->SetMarkerSize(0.5);
    hnlo->SetMarkerStyle(kFullCircle);
    hnlo->SetMarkerColor(kGreen+2);
    hnlo->Draw("SAME P");
    //if (iy==0) leg->AddEntry(hnlo, "PDF4LHC", "P");
    if (iy==0) leg->AddEntry(hnlo, "CT10", "P");
    //
    hnlo_ct10->SetMarkerSize(0.5);
    hnlo_ct10->SetMarkerStyle(kOpenSquare);
    hnlo_ct10->SetMarkerColor(kGreen+2);
    hnlo_ct10->SetLineColor(kGreen+2);
    hnlo_ct10->SetLineStyle(kDashed);
    hnlo_ct10->Draw("SAME L");
    if (iy==0) leg->AddEntry(hnlo_ct10, "CTEQ 10", "L");
    //
    hnlo_mstw->SetMarkerSize(0.5);
    hnlo_mstw->SetMarkerStyle(kOpenCircle);
    hnlo_mstw->SetMarkerColor(kGreen+2);
    hnlo_mstw->SetLineColor(kGreen+2);
    hnlo_mstw->SetLineStyle(kDashDotted);
    hnlo_mstw->Draw("SAME L");
    if (iy==0) leg->AddEntry(hnlo_mstw, "MSTW", "L");
    //
    hnlo_nnpdf->SetMarkerSize(0.5);
    hnlo_nnpdf->SetMarkerStyle(kOpenDiamond);
    hnlo_nnpdf->SetMarkerColor(kGreen+2);
    hnlo_nnpdf->SetLineColor(kGreen+2);
    hnlo_nnpdf->SetLineStyle(kDotted);
    hnlo_nnpdf->Draw("SAME L");
    if (iy==0) leg->AddEntry(hnlo_nnpdf, "NNPDF1.0", "L");

    TLatex *t = new TLatex(40, 0.15, etamin==0 ?
			   Form("|y| < %1.1f",etamax) :
			   Form("%1.1f #leq |y| < %1.1f",etamin,etamax));
    t->SetTextFont(43);
    t->SetTextSize(npix);
    t->Draw();

    if (iy==0) leg->Draw();
    if (iy==1) leg2->Draw();
    gPad->RedrawAxis();
  } // for iy

  c1->cd();

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextAlign(21); // align middle
  latex->DrawLatex(0.55,0.96,ak7 ? "Anti-k_{T} R=0.7" : "Anti-k_{T} R=0.5");
  //cmsPrel(jp::lumi, true);

  //if(_eps) c1->SaveAs(Form("eps/summaryTheoryDataRatio%s.eps",sid));
  if(jp::pdf) c1->SaveAs(Form("pdf/summaryTheoryDataRatio%s.pdf",sid));
} // drawTheoryDataRatio


void drawUnfoldSummary(string type, string algo) {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-4c.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");

  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};
  int color[] = {kGray+2, kBlack, kRed, kRed+2, kBlue, kBlue+2};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.40);
  h->SetMaximum(1.00);
  h->GetYaxis()->SetTitle("Correction");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->DrawClone("AXIS");

  TLegend *leg = new TLegend(0.43,0.23,0.73,0.53,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->Draw();

  TLatex *tjet = new TLatex(0.43,0.18,"Anti-k_{T} R=0.5");
  tjet->SetTextSize(0.045);
  tjet->SetNDC();
  tjet->Draw();

  if (type!="PF") {
    TLatex *talgo = new TLatex(0.43,0.54,(algo+" jets").c_str());
    talgo->SetTextSize(0.045);
    talgo->SetNDC();
    talgo->Draw();
  }

  //cmsPrel(jp::isdt ? jp::lumi : 0);

  vector<TGraphErrors*> gfolds(etas.size());
  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load and scale centra values
    TGraphErrors *gfold = (TGraphErrors*)d->Get("gfold"); assert(gfold);
    if (algo=="CALO") {
      gfold = (TGraphErrors*)d->Get("gfold_ak5calo"); assert(gfold);
    }

    gfold->SetMarkerStyle(marker[i]);
    gfold->SetMarkerColor(color[i]);
    gfold->SetLineColor(color[i]);
    gfold->SetLineWidth(2);
    gfold->Draw("SAMEPL");
    gfolds[i] = gfold;

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    leg->AddEntry(gfold,seta,"PL");

  } // for i

  // Save the histograms for redrawing PAS plots
//   if (_pas && algo=="PF") {
//     TDirectory *curdir = gDirectory;
//     TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
//     for (unsigned int i = 0; i != etas.size(); ++i) {

//       TGraphErrors *g = gfolds[i];
//       g->Write(Form("fig3_%d_gpt",i),TObject::kOverwrite);
//     }
//     fpas->Close();
//     curdir->cd();
//   } // _pas

  const char *t = type.c_str();
  const char *a = algo.c_str();
  if (algo!="PF") {
    //if(_eps) c1->SaveAs(Form("eps/summaryUnfolding_%s_%s.eps",a,t));
    if(jp::pdf) c1->SaveAs(Form("pdf/summaryUnfolding_%s_%s.pdf",a,t));
    //if(_gif) c1->SaveAs(Form("gif/summaryUnfolding_%s_%s.gif",a,t));
    //if(_png) c1->SaveAs(Form("png/summaryUnfolding_%s_%s.png",a,t));
  }
  else {
    //if(_eps) c1->SaveAs(Form("eps/summaryUnfolding_%s.eps",t));
    if(jp::pdf) c1->SaveAs(Form("pdf/summaryUnfolding_%s.pdf",t));
    //if(_gif) c1->SaveAs(Form("gif/summaryUnfolding_%s.gif",t));
    //if(_png) c1->SaveAs(Form("png/summaryUnfolding_%s.png",t));
  }

  { // PAS plots
    h->SetMinimum(0.50); // new c-term :/
    h->GetXaxis()->SetRangeUser(jp::xminpas,jp::xmax);
    h->DrawClone("AXIS");
    for (unsigned int i = 0; i != etas.size(); ++i) {
      gfolds[i]->Draw("SAMEPL");
    }
    leg->SetX1NDC(0.33);
    leg->SetX2NDC(0.63);
    leg->Draw();
    tjet->DrawLatex(0.33,0.18,"Anti-k_{T} R=0.5");
    //if (jp::isdt and algo=="PF" && _final) cmsFinal(jp::lumi);
    //else cmsPrel(jp::isdt ? jp::lumi : 0);

    if (algo!="PF") {
      //if(_eps) c1->SaveAs(Form("eps/summaryUnfoldingPAS_%s_%s.eps",a,t));
      if(jp::pdf) c1->SaveAs(Form("pdf/summaryUnfoldingPAS_%s_%s.pdf",a,t));
      //if(_gif) c1->SaveAs(Form("gif/summaryUnfoldingPAS_%s_%s.gif",a,t));
      //if(_png) c1->SaveAs(Form("png/summaryUnfoldingPAS_%s_%s.png",a,t));
    }
    else {
      //if(_eps) c1->SaveAs(Form("eps/summaryUnfoldingPAS_%s.eps",t));
      if(jp::pdf) c1->SaveAs(Form("pdf/summaryUnfoldingPAS_%s.pdf",t));
      //if(_gif) c1->SaveAs(Form("gif/summaryUnfoldingPAS_%s.gif",t));
      //if(_png) c1->SaveAs(Form("png/summaryUnfoldingPAS_%s.png",t));
    }
  }

  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawUnfolding


void drawJetIDSummary(string type, string idtype) {

  TDirectory *curdir = gDirectory;

  assert(idtype=="JetID" || idtype=="MetCut" || idtype=="JetIDPlusMetCut");

  const char *t = type.c_str();
  TFile *fin = new TFile(Form("output-%s-2b.root",t),"READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din1 = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0");

  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};
  int color[] = {kGray+2, kBlack, kRed, kRed+2, kBlue, kBlue+2};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TH1D *h = new TH1D("h","",100,jp::xmin,jp::xmax);
  if (idtype=="JetID") {
    h->GetYaxis()->SetTitle("Fraction passing Jet ID [%]");
    h->SetMaximum(101.0001);
    h->SetMinimum(90.0000);
  }
  if (idtype=="MetCut") {
    h->GetYaxis()->SetTitle("Fraction passing MET cut [%]");
    h->SetMaximum(105-0.0001);
    h->SetMinimum(60.0000);
  }
  if (idtype=="JetIDPlusMetCut") {
    h->GetYaxis()->SetTitle("Fraction passing Jet ID + MET cut [%]");
    h->SetMaximum(110-0.0001);
    h->SetMinimum(0.0000);
  }
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->Draw("AXIS");

  TLegend *leg = new TLegend(0.40,0.23,0.70,0.58,t,"brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();

  TLatex *tjet = new TLatex(0.18,0.18,"Anti-k_{T} R=0.5");// PFJets");
  tjet->SetTextSize(0.04);
  tjet->SetNDC();
  tjet->Draw();

  //cmsPrel(type=="DATA" ? jp::lumi : 0);

  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din1->cd(etas[i].c_str()));
    TDirectory *d1 = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load JetID
    TH1D *h1 = (TH1D*)d1->Get("hpt");
    TH1D *h2 = 0;
    if (idtype=="JetID")
      h2 = (TH1D*)d1->Get("hpt_nojetid"); // no JetID, but MET cut
    if (idtype=="MetCut")
      h2 = (TH1D*)d1->Get("hpt_noevtid"); // no MET cut, but JetID
    if (idtype=="JetIDPlusMetCut")
      h2 = (TH1D*)d1->Get("hpt_noid"); // no JetID or MET cut
    assert(h1);
    assert(h2);

    // Add missing 3/pb patch for hpt_nojetid
    /*
    if (type=="DATA") {
      TH1D *peff = (TH1D*)d1->Get("peff_new"); assert(peff);
      assert(h2->GetNbinsX()==peff->GetNbinsX());
      for (int j = 1; j != peff->GetNbinsX()+1; ++j) {
	double eff = peff->GetBinContent(j);
	h2->SetBinContent(j, eff ? h2->GetBinContent(j)/eff : 0.);
	h2->SetBinError(j, eff ? h2->GetBinError(j)/eff : 0.);
      } // for j
    } // 3/pb patch
    */

    TH1D *hid = (TH1D*)h1->Clone(Form("hid%d",i));
    hid->Divide(h1, h2, 1, 1, "B");

    hid->Scale(100.);
    hid->SetMarkerStyle(marker[i]);
    hid->SetMarkerColor(color[i]);
    hid->SetLineColor(color[i]);
    hid->SetLineWidth(2);
    hid->Draw("SAMEPL");

    TF1 *f = new TF1(Form("fid%d",i),"[0]-[1]*exp([2]*x)",jp::xmin,jp::xmax);
    f->SetParameters(100,2,-0.1);
    f->SetParLimits(1,0.,10.);
    f->SetLineColor(color[i]);
    if (idtype=="JetIDPlusMetCut" || (idtype=="MetCut" && etamin>=2.5))
      f->SetRange(jp::xmin, 200.);
    hid->Fit(f,i==5 ? "QRNW" : "QRN");
    f->SetRange(jp::xmin,jp::xmax);
    f->Draw("SAME");

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    leg->AddEntry(hid,seta,"PL");

  } // for i

  const char *id = idtype.c_str();
  //if (_eps) c1->SaveAs(Form("eps/summary%s_%s.eps",id,t));
  if (jp::pdf) c1->SaveAs(Form("pdf/summary%s_%s.pdf",id,t));
  //if (_gif) c1->SaveAs(Form("gif/summary%s_%s.gif",id,t));
  //if (_png) c1->SaveAs(Form("png/summary%s_%s.png",id,t));

  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawJetID


void drawNPCorr(string type) {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-5.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  assert(fin->cd("Standard"));
  TDirectory *din = gDirectory;

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  etas.push_back("Eta_2.0-2.5");
  etas.push_back("Eta_2.5-3.0"); // missing NLO theory

  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};
  int color[] = {kGray+2, kBlack, kRed, kRed+2, kBlue, kBlue+2};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.8);
  h->SetMaximum(1.5);
  h->GetYaxis()->SetTitle("Non-perturbative corrections");
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->Draw("AXIS");

  TLegend *leg = new TLegend(0.58,0.63,0.93,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  TLatex *tjet = new TLatex(0.18,0.18,"Anti-k_{T} R=0.5 Jets");
  tjet->SetTextSize(0.05);
  tjet->SetNDC();

  vector<TGraph*> vnp(etas.size());
  vector<TH1D*> vhnp(etas.size());
  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din->cd(etas[i].c_str()));
    TDirectory *d = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TH1D *hnp = (TH1D*)d->Get("hnpcorr"); assert(hnp);
    TGraph *gnp = new TGraph(hnp);
    gnp->SetMarkerStyle(marker[i]);
    gnp->SetMarkerColor(color[i]);
    gnp->SetLineColor(color[i]);
    gnp->SetLineWidth(2);
    gnp->Draw("SAMEPL");
    vnp[i] = gnp;
    vhnp[i] = hnp;

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    leg->AddEntry(gnp,seta,"P");

  } // for i

  // Hide the biased pThat<15 GeV region (pt<20 GeV or so)
  TBox *b = new TBox(10,0.8,21.,1.5);
  b->SetFillStyle(kCrossHatch);
  b->SetFillColor(kGray);
  gStyle->SetHatchesLineWidth(2);
  b->Draw();

  // Get the tick marks on top of the bands
  h->Draw("AXIS SAME");

  tjet->Draw();
  //cmsPrel(0);
  leg->Draw();

  //if(_eps) c1->SaveAs("eps/summaryNPCorr.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryNPCorr.pdf");
  //if(_gif) c1->SaveAs("gif/summaryNPCorr.gif");
  //if(_png) c1->SaveAs("png/summaryNPCorr.png");

  /////////////////////////////////////////////////////////////////////////
  // for PAS, cut out the low pT part
  h->GetXaxis()->SetRangeUser(18.,jp::xmax);
  h->SetMaximum(1.45);
  h->SetMinimum(0.85);
  h->Draw("AXIS");
  double vptmax[6] = {846, 790, 686, 594, 468, 330};
  for (unsigned int i = 0; i != vnp.size(); ++i) {

    TGraph *g = vnp[i];
    for (int j = g->GetN()-1; j!=-1; --j) {
      if (g->GetX()[j]>vptmax[i]) g->RemovePoint(j);
    }

    vnp[i]->Draw("SAMEPL");
  }
  h->Draw("AXIS SAME");
  tjet->Draw();
  leg->Draw();

  c1->Update();

  //if(_eps) c1->SaveAs("eps/summaryNPCorrPAS.eps");
  if(jp::pdf) c1->SaveAs("pdf/summaryNPCorrPAS.pdf");
  //if(_gif) c1->SaveAs("gif/summaryNPCorrPAS.gif");
  //if(_png) c1->SaveAs("png/summaryNPCorrPAS.png");


  // Durham-style tables for Klaus
  // Print out the Durham tables
  gROOT->ProcessLine(".! mkdir txt");
  for (unsigned int i = 0; i != vnp.size(); ++i) {
    ofstream fout(Form("txt/CMS_incjet_npcorr_y%1.1f_%1.1f.txt",
		       0.5*i,0.5*(i+1)), ios::out);
    fout << Form("Non-perturbative corrections to inclusive jet double-differential cross sections in the |rapidity| range %1.1f to %1.1f, using a jet resolution R value of %1.1f.",0.5*i,0.5*(i+1),0.5) << endl;
    fout << "RE : P P --> JET X" << endl;
    fout << Form("ABS(YRAP) : %1.1f TO %1.1f",0.5*i,0.5*(i+1)) << endl;
    fout << Form("R : %1.1f",0.5) << endl;
    fout << "SQRT(S) : 7 TeV" << endl;
    fout << "x-axis header: PT IN GEV" << endl;
    fout << "y-axis header: CORRECTION TO THEORY" << endl;
    fout << endl;
    fout << Form("  %-6s %-5s %-5s    %-10s  %-21s %-21s",
		 "x","xlow","xhigh","y","dy(stat)","dy(sys)") << endl;
    fout << endl;
    TGraphErrors *g = new TGraphErrors(vhnp[i]);
    TH1D *h = vhnp[i];
    for (int j = 0; j != g->GetN(); ++j) {

      double x = g->GetX()[j];
      if (x>18) {
	int ibin = h->FindBin(x);
	double xmin = h->GetBinLowEdge(ibin);
	double xmax = h->GetBinLowEdge(ibin+1);
	double y = g->GetY()[j];
	double ey = g->GetEY()[j] / y * 100.;
	double sysup = h->GetBinError(ibin);//*100.;
	double sysdw = h->GetBinError(ibin);//*100.;
	fout << Form("%6.1f %5.0f %5.0f  %10.3g  %+5.1f %+5.1f (percent)"
		     //"  %+5.1f %+5.1f (percent)",
		     "  %+5.2f %+5.2f",
		     x, xmin, xmax, y, ey, -ey, sysup, sysdw) << endl;
      }
    } // for j
    fout.close();
  } // for i


  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawNPCorr


void drawClosure(string type, bool isb, string algo="PF", bool is38x=false) {

  TDirectory *curdir = gDirectory;

  assert(type=="MC" || type=="HW");
  assert(algo=="PF" || algo=="PFAK7" || algo=="JPT" || algo=="CALO");
  assert(!(isb && algo!="PF"));

  const char *t = type.c_str();

  TFile *fin1mc = new TFile(Form("output-%s-2a.root",type.c_str()),"READ");
  assert(fin1mc && !fin1mc->IsZombie());
  assert(fin1mc->cd("Standard"));
  TDirectory *din1mc = gDirectory;
  //
  TFile *fin1tw = new TFile(Form("output-%s-2b.root",t),"READ");
  assert(fin1tw && !fin1tw->IsZombie());
  assert(fin1tw->cd("Standard"));
  TDirectory *din1tw = gDirectory;

  //TFile *fin2 = new TFile(Form("output-%s-3a.root",type.c_str()),"READ");
  ////
  //TFile *fin2 = new TFile(Form("output-%s-3b.root",t),"READ");
  //TFile *fin2 = new TFile(Form("output-%s-3c.root",t),"READ");
  TFile *fin2 = new TFile(Form("output-%s-3.root",t),"READ");
  assert(fin2 && !fin2->IsZombie());
  assert(fin2->cd("Standard"));
  TDirectory *din2 = gDirectory;

  //TFile *fin2b = new TFile(Form("output-%s-4c.root",t),"READ");
  TFile *fin2b = new TFile(Form("output-%s-3.root",t),"READ");
  assert(fin2b && !fin2b->IsZombie());
  assert(fin2b->cd("Standard"));
  TDirectory *din2b = gDirectory;

  //TFile *fin3 = new TFile(Form("output-%s-4.root",t),"READ");
  //assert(fin3 && !fin3->IsZombie());
  //assert(fin3->cd("Standard"));
  //TDirectory *din3 = gDirectory;

  //TFile *frap = new TFile("rapiditybias.root","READ");
  //assert(frap && !frap->IsZombie());

  TFile *fclosure = (isb ? 0 :  new TFile(Form("closure-%s.root",t),
					  "RECREATE"));

  vector<string> etas;
  etas.push_back("Eta_0.0-0.5");
  etas.push_back("Eta_0.5-1.0");
  etas.push_back("Eta_1.0-1.5");
  etas.push_back("Eta_1.5-2.0");
  if (!isb) etas.push_back("Eta_2.0-2.5");
  if (!isb) etas.push_back("Eta_2.5-3.0");
  if (isb) etas.push_back("Eta_2.0-2.2");

  // maximum pT in data
  //double maxpt[6] = {1497, 1410, 1172, 905, 638, 430};
  double maxpt[6] = {2500, 2238, 2000, 1410, 967, 638};

  int marker[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		  kFullTriangleUp, kOpenTriangleUp};
  int color[] = {kGray+2, kBlack, kRed, kRed+2, kBlue, kBlue+2};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->SetLogx();

  TH1D *h = new TH1D("h","",int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h->SetMinimum(0.5);//isb ? 0.2 : 0.6);
  h->SetMaximum(1.4);//isb ? 1.6 : 1.3);
  //h->GetYaxis()->SetTitle(Form("Reconstructed (%s) / Generated",algo.c_str()));
  bool match = false;//true;
  h->GetYaxis()->SetTitle(Form("Reco / Gen (%s)",match ? "matched" : "all"));
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.28,0.18,0.63,0.48,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  TLegend *leg2 = new TLegend(0.18,0.87,0.53,0.92,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);

  c1->cd();
  h->Draw("AXIS");
  leg->Draw();
  leg2->Draw();

  c2->cd();
  h->Draw("AXIS");
  leg->Draw();
  leg2->Draw();

  // draw widest uncertainty band at the back
  // TH1D *hup = 0;
//   TH1D *hdw = 0;
//   for (unsigned int i = 0; i != etas.size(); ++i) {

//     assert(din3->cd(etas[i].c_str()));
//     TDirectory *d3 = gDirectory;

//     TH1D *up = (TH1D*)d3->Get("hjec_tot_pl"); assert(up);
//     TH1D *dw = (TH1D*)d3->Get("hjec_tot_mn"); assert(dw);

//     for (int j = 1; j != up->GetNbinsX()+1; ++j) {
//       if (up->GetBinCenter(j) < maxpt[i]) {
// 	up->SetBinContent(j, 1+up->GetBinContent(j));
// 	dw->SetBinContent(j, 1+dw->GetBinContent(j));
//       }
//       else {
// 	up->SetBinContent(j ,1.);
// 	dw->SetBinContent(j ,1.);
//       }
//     }
//     if (!hup) hup = (TH1D*)up->Clone("hup");
//     if (!hdw) hdw = (TH1D*)dw->Clone("hdw");

//     //if (i==0)
//     for (int j = 1; j != up->GetNbinsX()+1; ++j) {

//       if (up->GetBinContent(j) > hup->GetBinContent(j))
// 	hup->SetBinContent(j, up->GetBinContent(j));
//       if (dw->GetBinContent(j) < hdw->GetBinContent(j))
// 	hdw->SetBinContent(j, dw->GetBinContent(j));
//     }
//   }
//   TF1 *fcenter = new TF1("fcenter","1", jp::xmin, jp::xmax);
//   TH1D *hpl = (TH1D*)hup->Clone("hpl"); offsetGraph(hpl, -1);
//   TH1D *hmn = (TH1D*)hdw->Clone("hmn"); offsetGraph(hmn, -1);
//   TGraphErrors *gerr = makeBand(fcenter, hpl, hmn, jp::xmin, jp::xmax);

//   gerr->SetFillStyle(1001);
//   gerr->SetFillColor(kYellow);
//   hup->SetLineWidth(2);
//   hdw->SetLineWidth(2);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  tex->SetNDC();

  c1->cd();
  //gerr->Draw("SAME E3");
  //hup->Draw("SAMEL");
  //hdw->Draw("SAMEL");
  l->DrawLine(jp::xmin,1,jp::xmax,1);
  tex->DrawLatex(0.80,0.20,string(jp::algo).c_str());
  c2->cd();
  //gerr->Draw("SAME E3");
  //hup->Draw("SAMEL");
  //hdw->Draw("SAMEL");
  l->DrawLine(jp::xmin,1,jp::xmax,1);
  tex->DrawLatex(0.80,0.20,string(jp::algo).c_str());

  //leg2->AddEntry(gerr, "MC truth JEC", "F");

  for (unsigned int i = 0; i != etas.size(); ++i) {

    assert(din1mc->cd((etas[i]+"/mc").c_str()));
    TDirectory *d1mc = gDirectory;
    //
    assert(din1tw->cd(etas[i].c_str()));
    TDirectory *d1tw = gDirectory;

    assert(din2->cd(etas[i].c_str()));
    TDirectory *d2 = gDirectory;

    assert(din2b->cd(etas[i].c_str()));
    TDirectory *d2b = gDirectory;

    float etamin, etamax;
    assert(sscanf(etas[i].c_str(),"Eta_%f-%f",&etamin,&etamax)==2);

    // Load central values
    TH1D *hgen = 0;
    if (match)  hgen = (TH1D*)d1mc->Get("hpt_g"); // match+ID+yreco
    if (!match) hgen = (TH1D*)d1mc->Get("hpt_g0"); // all ygen
    assert(hgen);
    TH1D *hreco = 0;
    hreco = (TH1D*)d2->Get("hcorrpt");
    assert(hreco);
    hreco->Divide(hgen);
    //
    TH1D *hreco2 = 0;
    hreco2 = (TH1D*)d2b->Get("hcorrpt");
    assert(hreco2);
    hreco2->Divide(hgen);

    hreco->SetLineColor(color[i]);
    hreco->SetMarkerColor(color[i]);
    hreco->SetMarkerStyle(marker[i]);
    hreco2->SetLineColor(color[i]);
    hreco2->SetMarkerColor(color[i]);
    hreco2->SetMarkerStyle(marker[i]);

    // Load correction for central value, based on JEC residual
//     if (false) { // fix residual JEC

//       if (i==0)
// 	cout << "\n***ATTENTION: Correcting for JEC residual in closure!***\n"
// 	     << endl;

//       assert(din3->cd(etas[i].c_str()));
//       TDirectory *d3 = gDirectory;

//       TH1D *hc = (TH1D*)d3->Get("hptshift_pf"); assert(hc);
//       for (int j = 1; j != hreco2->GetNbinsX()+1; ++j) {
//       double c = 1./(1+hc->GetBinContent(j));
//       hreco2->SetBinContent(j, c*hreco2->GetBinContent(j));
//       hreco2->SetBinError(j, c*hreco2->GetBinError(j));
//       } // for j
//     } // fix residual JEC

    const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
    leg->AddEntry(hreco,seta,"P");

    // Do some cleaning on the closure to remove huge fluctuations
    // but for inclusive jets only (stats poor for b-jets anyway)
    if (!isb) {

      for (int j = 1+1; j != hreco->GetNbinsX()+1-1; ++j) {
	if (hreco->GetBinError(j)>0.05) {
	  double y0 = hreco->GetBinContent(j-1);
	  double ey0 = hreco->GetBinError(j-1);
	  double y1 = hreco->GetBinContent(j+1);
	  double ey1 = hreco->GetBinError(j+1);
	  hreco->SetBinContent(j, 0.5*(y0+y1));
	  hreco->SetBinError(j, sqrt(0.5*ey0*ey0+0.5*ey1*ey1+(y1-y0)*(y1-y0)));
	}
      } // for j
    } // !isb
    for (int j = 1; j != hreco->GetNbinsX()+1; ++j) {
      if (hreco->GetBinCenter(j) > maxpt[i]) {
	hreco->SetBinContent(j, 0.);
	hreco->SetBinError(j, 0.);
      }
    } // for j
    // Same for corrected numbers
    if (!isb) {

      for (int j = 1+1; j != hreco2->GetNbinsX()+1-1; ++j) {
	if (hreco2->GetBinError(j)>0.05) {
	  double y0 = hreco2->GetBinContent(j-1);
	  double ey0 = hreco2->GetBinError(j-1);
	  double y1 = hreco2->GetBinContent(j+1);
	  double ey1 = hreco2->GetBinError(j+1);
	  hreco2->SetBinContent(j, 0.5*(y0+y1));
	  hreco2->SetBinError(j, sqrt(0.5*ey0*ey0+0.5*ey1*ey1+(y1-y0)*(y1-y0)));
	}
      } // for j
    } // !isb
    for (int j = 1; j != hreco2->GetNbinsX()+1; ++j) {
      if (hreco2->GetBinCenter(j) > maxpt[i]) {
	hreco2->SetBinContent(j, 0.);
	hreco2->SetBinError(j, 0.);
      }
    } // for j

    c1->cd();
    hreco->GetXaxis()->SetRangeUser(h->GetXaxis()->GetXmin(),
				    h->GetXaxis()->GetXmax());
    hreco->DrawClone("SAME P");
    c2->cd();
    hreco2->GetXaxis()->SetRangeUser(h->GetXaxis()->GetXmin(),
				     h->GetXaxis()->GetXmax());
    hreco2->DrawClone("SAME P");

    if (!isb && algo=="PF" && is38x) {
      fclosure->cd();
      hreco->Write(Form("closure_Eta_%1.1f-%1.1f",etamin,etamax));
    }
  } // for i

  // Get the tick marks on top of the bands
  c1->cd();
  h->Draw("AXIS SAME");
  //cmsPrel(0);
  leg->Draw();

  c2->cd();
  h->Draw("AXIS SAME");
  //cmsPrel(0);
  leg->Draw();

  const char *a = algo.c_str();
  const char *v = (is38x ? "_38X" : "");
  const char *tc = string(jp::algo).c_str();
  if (isb) {
    //if(_eps) c1->SaveAs("eps/summaryBClosure.eps");
    if(jp::pdf) c1->SaveAs("pdf/summaryBClosure.pdf");
    //if(_gif) c1->SaveAs("gif/summaryBClosure.gif");
    //if(_png) c1->SaveAs("png/summaryBClosure.png");

    //if(_eps) c1->SaveAs("eps/summaryBClosureCorr.eps");
    if(jp::pdf) c1->SaveAs("pdf/summaryBClosureCorr.pdf");
    //if(_gif) c1->SaveAs("gif/summaryBClosureCorr.gif");
    //if(_png) c1->SaveAs("png/summaryBClosureCorr.png");
  }
  else {
    //if(_eps) c1->SaveAs(Form("eps/summaryClosure%s_%s%s%s.eps",t,a,v,tc));
    if(jp::pdf) c1->SaveAs(Form("pdf/summaryClosure%s_%s%s%s.pdf",t,a,v,tc));
    //if(_gif) c1->SaveAs(Form("gif/summaryClosure%s_%s%s%s.gif",t,a,v,tc));
    //if(_png) c1->SaveAs(Form("png/summaryClosure%s_%s%s%s.png",t,a,v,tc));

    //if(_eps) c2->SaveAs(Form("eps/summaryClosureCorr%s_%s%s%s.eps",t,a,v,tc));
    if(jp::pdf) c2->SaveAs(Form("pdf/summaryClosureCorr%s_%s%s%s.pdf",t,a,v,tc));
    //if(_gif) c2->SaveAs(Form("gif/summaryClosureCorr%s_%s%s%s.gif",t,a,v,tc));
    //if(_png) c2->SaveAs(Form("png/summaryClosureCorr%s_%s%s%s.png",t,a,v,tc));
  }

  if (gROOT->IsBatch()) delete c1;
  fclosure->Close();
  fclosure->Delete();

  curdir->cd();

} // drawClosure

// void GetPoint(TGraphErrors *g, int n, double &x, double &y,
// 	      double &ex, double &ey) {
//   g->GetPoint(n, x, y);
//   ex = g->GetErrorX(n);
//   ey = g->GetErrorY(n);
// }
// void SetPoint(TGraphErrors *g, int n, double x, double y,
// 	      double ex, double ey) {
//   g->SetPoint(n, x, y);
//   g->SetPointError(n, ex, ey);
// }
void scaleGraph(TGraphErrors *g, double scale) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(g, i, x, y, ex, ey);
    tools::SetPoint(g, i, x, scale*y, ex, scale*ey);
  }
}
void offsetGraph(TGraphErrors *g, double offset) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(g, i, x, y, ex, ey);
    tools::SetPoint(g, i, x, y+offset, ex, ey);
  }
}
void offsetGraph(TH1D *g, double offset) {

  for (int i = 1; i != g->GetNbinsX()+1; ++i) {
    g->SetBinContent(i, g->GetBinContent(i)+offset);
  }
}
void ratioGraph(TGraphErrors *g, TF1 *f) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(g, i, x, y, ex, ey);
    tools::SetPoint(g, i, x, y/f->Eval(x), ex, ey/f->Eval(x));
  }
}
void ratioGraph(TGraphErrors *g, TH1D *h) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(g, i, x, y, ex, ey);
    int j = h->FindBin(x);
    double y0 = h->GetBinContent(j);
    tools::SetPoint(g, i, x, y/y0, ex, ey/y0);
  }
}
void pruneGraph(TGraphErrors *g, double errmax, double ymin,
		double xmin, double xmax) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(g, i, x, y, ex, ey);
    if ((y!=0 && fabs(ey / y) > errmax) || fabs(y)<ymin ||
	x < xmin || x > xmax) {
      g->RemovePoint(i); --i;
    }
  }
}
void rangeGraph(TGraph *g, double xmin, double xmax) {

  for (int i = 0; i != g->GetN(); ++i) {
    double x = g->GetX()[i];
    if (x < xmin || x > xmax) {
      g->RemovePoint(i); --i;
    }
  } // for i
} // rangeGraph


TGraphErrors *makeBand(TGraph *g, const TH1D *epl, const TH1D *emn,
		       double xmin, double xmax, bool hist) {

  TGraphErrors *gerr = new TGraphErrors(0);
  gerr->SetName(Form("eband_%s",epl->GetName()));
  for (int i = 0; i != g->GetN(); ++i) {

    double x, y;
    g->GetPoint(i, x, y);
    if (x>=xmin && x<=xmax) {

      int jpl = ((TH1D*)epl)->FindBin(x);
      int jmn = ((TH1D*)emn)->FindBin(x);

      double ypl = (1. + fabs(epl->GetBinContent(jpl))) * y;
      double ymn = (1. - fabs(emn->GetBinContent(jmn))) * y;
      double ymid = 0.5 * (ypl + ymn);
      double ey = 0.5 * (ypl - ymn);
      double ex = 0;

      if (hist) {
	x = epl->GetBinCenter(jpl);
	ex = 0.5*epl->GetBinWidth(jpl);
      }


      tools::SetPoint(gerr, gerr->GetN(), x, ymid, ex, ey);
    }
  } // for i
  gerr->SetFillStyle(1001);
  gerr->SetFillColor(kYellow+1);

  return gerr;
} // makeBand(TGraph)

TGraphErrors *makeBand(TF1 *f, const TH1D *epl, const TH1D *emn,
		       double xmin, double xmax) {

  TH1D *h = (TH1D*)epl->Clone(Form("tmp_%s",epl->GetName()));
  assert(epl->GetNbinsX()==emn->GetNbinsX());

  int imin = h->FindBin(xmin);
  vector<double> vx(1,h->GetBinLowEdge(imin));
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    // evaluate uncertainty at bin center
    double x = h->GetBinCenter(i);
    double y = f->Eval(x);
    double ypl = (1. + fabs(epl->GetBinContent(i))) * y;
    double ymn = (1. - fabs(emn->GetBinContent(i))) * y;
    double ymid = 0.5 * (ypl + ymn);
    double ey = 0.5 * (ypl - ymn);

    h->SetBinContent(i, ymid);
    h->SetBinError(i, ey);
    if (h->GetBinLowEdge(i) >= xmin && h->GetBinLowEdge(i) < xmax) {
      vx.push_back(h->GetBinLowEdge(i+1));
    }
  } // for i

  TH1D *hout = new TH1D(Form("eband_%s",epl->GetName()),"",vx.size()-1,&vx[0]);
  for (int i = 1; i != hout->GetNbinsX()+1; ++i) {

    int ibin = h->FindBin(hout->GetBinCenter(i));
    hout->SetBinContent(i, h->GetBinContent(ibin));
    hout->SetBinError(i, h->GetBinError(ibin));
  }
  hout->SetMarkerStyle(kNone);
  TGraphErrors *gerr = new TGraphErrors(hout);
  gerr->SetFillStyle(1001);
  gerr->SetFillColor(kYellow+1);

  delete h;
  delete hout;
  return gerr;
} // makeBand(TF1)
