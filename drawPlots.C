// Purpose: Draw plots of analysis stages
// Author:  mikko.voutilainen@cern.ch
// Created: March 21, 2010
// Updated: March 21, 2010
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TLine.h"
#include "TROOT.h"
#include "TArrow.h"
#include "THStack.h"
#include "TMultiGraph.h"

#include "tools.h"
#include "tdrstyle_mod12.C"

#include <iostream>
#include <string>
#include <map>

using namespace std;

// ranges for drawing
#include "settings.h"


map<string, map<int, pair<double, double> > > _trgratios;
map<string, string> _trgrationames;

// possible commands
void drawJetID(string type);
void drawTriggerRatio(string type, string algo="PF");
void drawUnfolding(string type, string algo="PF");
void drawResolution(string type);
void drawMultijets(string type);

// b-tagging specific commands
//void drawBFraction(string type);
//void drawBTagPurity(string type);
//void drawBTagEfficiency(string type);
//void drawBTagNorm(string type);

// Homogenizing plot styles for PAS
//void drawMCPurity();
//void drawMCEfficiency();
//void drawTemplatePurity();
//void drawTemplate();

void drawJetIDBin(TDirectory *did, TDirectory *dnd, int ibin, string type);

void drawJetID(string type) {

  TFile *fin = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");

  // Select top categories for JEC uncertainty
  assert(fin->cd("Standard"));
  TDirectory *dtop0 = gDirectory;
  //assert(fin->cd("NoJetID"));
  //if(!fin->cd("NoJetID"))
  //assert(fin->cd("Standard"));
  TDirectory *dtop1 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = dtop0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")) {

      assert(dtop0->cd(obj->GetName()));
      TDirectory *din0 = gDirectory;
      assert(dtop1->cd(obj->GetName()));
      TDirectory *din1 = gDirectory;

      // Process subdirectory
      drawJetIDBin(din0, din1, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawPlots

// Draw comparison of raw pT spectra before and after JetID
void drawJetIDBin(TDirectory *did, TDirectory *dnd, int ibin, string type) {

  gStyle->SetOptStat(0);

  float etamin, etamax;
  assert(sscanf(did->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char* seta = (etamin==0 ? Form("|y| #leq %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.60,0.85,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  TH1D *hid = (TH1D*)did->Get("hpt"); assert(hid);
  //TH1D *hnd = (TH1D*)dnd->Get("hpt"); assert(hnd);
  TH1D *hnd = (TH1D*)did->Get("hpt_nojetid"); assert(hnd);
  //TH1D *hnd = (TH1D*)did->Get("hpt_noid"); assert(hnd);

  // Comparison of spectra
  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name,c1name,600,600);
  c1->SetLogy();
  c1->SetLogx();

  //hnd->GetYaxis()->SetTitle("number of jets");
  TH1D *h = new TH1D(Form("hid_%d",ibin),"",int(xmax-xmin),xmin,xmax);
  h->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
  h->GetYaxis()->SetTitleOffset(1.15);//2);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMinimum(1e-3);
  h->SetMaximum(1e11);
  h->DrawClone("AXIS");

  hnd->SetLineColor(kRed);
  hnd->SetLineWidth(2);
  hnd->Draw("SAME");

  hid->SetLineColor(kBlack);
  hid->SetLineWidth(2);
  hid->Draw("SAME");

  teta->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if(_eps) c1->SaveAs(Form("eps/jetID_ptSpectra_%s_Eta%d.eps",t,k));
  if(_pdf) c1->SaveAs(Form("pdf/jetID_ptSpectra_%s_Eta%d.pdf",t,k));
  if(_gif) c1->SaveAs(Form("gif/jetID_ptSpectra_%s_Eta%d.gif",t,k));
  if(_png) c1->SaveAs(Form("png/jetID_ptSpectra_%s_Eta%d.png",t,k));
  if (gROOT->IsBatch()) delete c1;


  // Ratio of spectra
  const char *c2name = Form("c2_%d",ibin);
  TCanvas *c2 = new TCanvas(c2name,c2name,600,600);
  c2->SetLogx();

  //TH1D *hr = (TH1D*)hnd->Clone("hr");
  TH1D *hr = (TH1D*)hid->Clone(Form("hr_%d",ibin));
  hr->Divide(hr, hnd, 1, 1, "B");

  h->GetYaxis()->SetTitle("Fraction of jets passing JetID");
  h->SetMaximum(1.02);
  h->SetMinimum(0.95);
  h->DrawClone("AXIS");

  hr->Draw("SAME");

  teta->Draw();

  if(_eps) c2->SaveAs(Form("eps/jetID_fracPass_%s_Eta%d.eps",t,k));
  if(_pdf) c2->SaveAs(Form("pdf/jetID_fracPass_%s_Eta%d.pdf",t,k));
  if(_gif) c2->SaveAs(Form("gif/jetID_fracPass_%s_Eta%d.gif",t,k));
  if(_png) c2->SaveAs(Form("png/jetID_fracPass_%s_Eta%d.png",t,k));
  if (gROOT->IsBatch()) delete c2;


  // Ratio of spectra
  const char *c3name = Form("c3_%d",ibin);
  TCanvas *c3 = new TCanvas(c3name,c3name,600,600);
  c3->SetLogx();

  //TH1D *hr2 = (TH1D*)hnd->Clone("hr2");
  TH1D *hr2 = (TH1D*)hnd->Clone(Form("hr2_%d",ibin));
  hr2->Divide(hnd, hid, 1, 1, "B");

  h->GetYaxis()->SetTitle("Ratio of jets before JetID");
  h->DrawClone("AXIS");

  hr2->Draw("SAME");

  TH1D *hr2s = (TH1D*)hr2->Clone(Form("hr2s_%d",ibin));
  hr2s->Smooth(10);
  hr2s->SetLineColor(kBlack);
  hr2s->SetMarkerStyle(kFullCircle);
  hr2s->Draw("SAMEP");

  teta->Draw();

  if(_eps) c3->SaveAs(Form("eps/jetID_excessBefore_%s_Eta%d.eps",t,k));
  if(_pdf) c3->SaveAs(Form("pdf/jetID_excessBefore_%s_Eta%d.pdf",t,k));
  if(_gif) c3->SaveAs(Form("gif/jetID_excessBefore_%s_Eta%d.gif",t,k));
  if(_png) c3->SaveAs(Form("png/jetID_excessBefore_%s_Eta%d.png",t,k));
  if (gROOT->IsBatch()) delete c3;


  // Calculate purity after JetID, assuming JetID is 80% efficient
  // in removing background, and 98% efficient for signal
  // before N0 = S + B = S'/s + B'/b = B'/s*(S'/B' + s/b)
  // after  N1 = s*S + b*B = S' + B' = B'*(S'/B' + 1)
  // r = N0 / N1 = 1/s * (p + s/b) / (p+1)
  // => r*p + r = 1/s*p + 1/b
  // => (r-1/s)*p = 1/b-r => p = (1/b-r)/(r-1/s)
  // X = (S+B)/S = 1 + 1/p
  const char *c4name = Form("c4_%d",ibin);
  TCanvas *c4 = new TCanvas(c4name,c4name,600,600);
  c4->SetLogx();

  TH1D *hp = (TH1D*)hr2s->Clone(Form("hp_%d",ibin));
  for (int i = 1; i != hp->GetNbinsX()+1; ++i) {

    double r = hr2s->GetBinContent(i);
    const double seff = 0.98;
    const double beff = 0.20;
    double p = (1./beff - r) / (r - 1./seff);
    double X = 1 + 1/p;
    hp->SetBinContent(i, X);
  }

  h->GetYaxis()->SetTitle("Increase due to impurity");
  h->DrawClone("AXIS");

  hp->Draw("SAME");
  teta->Draw();

  if(_eps) c4->SaveAs(Form("eps/jetID_incrEst_%s_Eta%d.eps",t,k));
  if(_pdf) c4->SaveAs(Form("pdf/jetID_incrEst_%s_Eta%d.pdf",t,k));
  if(_gif) c4->SaveAs(Form("gif/jetID_incrEst_%s_Eta%d.gif",t,k));
  if(_png) c4->SaveAs(Form("png/jetID_incrEst_%s_Eta%d.png",t,k));
  if (gROOT->IsBatch()) delete c4;

} // drawJetID


void drawTriggerRatioBin(TDirectory *din0, TDirectory *din1, int ibin,
			 string type, string algo,
			 TDirectory *dmc0 = 0, TDirectory *dmc1 = 0);

void drawTriggerRatio(string type, string algo) {

  TFile *fin0 = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin0 && !fin0->IsZombie());
  TFile *fin1 = new TFile(Form("output-%s-2a.root",type.c_str()),"READ");
  assert(fin1 && !fin1->IsZombie());

  TFile *fmc0 = new TFile("output-MC-2b.root","READ");
  assert(fmc0 && !fmc0->IsZombie());
  TFile *fmc1 = new TFile("output-MC-2a.root","READ");
  assert(fmc1 && !fmc1->IsZombie());

  // Select top categories for JEC uncertainty
  assert(fin0->cd("Standard"));
  TDirectory *dtop0 = gDirectory;
  assert(fin1->cd("Standard"));
  TDirectory *dtop1 = gDirectory;

  assert(fmc0->cd("Standard"));
  TDirectory *dmctop0 = gDirectory;
  assert(fmc1->cd("Standard"));
  TDirectory *dmctop1 = gDirectory;

  _trgratios.clear();

  // Automatically go through the list of keys (directories)
  TList *keys = dtop0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin1 = 0, ibin2 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(dtop0->cd(obj->GetName()));
      TDirectory *din0 = gDirectory;
      assert(dtop1->cd(obj->GetName()));
      TDirectory *din1 = gDirectory;

      assert(dmctop0->cd(obj->GetName()));
      TDirectory *dmc0 = gDirectory;
      assert(dmctop1->cd(obj->GetName()));
      TDirectory *dmc1 = gDirectory;

      // Process subdirectory
      //if (jp::isdt and algo=="PF")
      //drawTriggerRatioBin(din0, din1, ibin1++, type, algo, dmc0, dmc1);
      drawTriggerRatioBin(din0, din1, ibin2++, type, algo);
    } // inherits TDirectory
  } // while

} // drawTriggerRatio

void drawTriggerRatioBin(TDirectory *din0, TDirectory *din1, int ibin,
			 string type, string algo,
			 TDirectory *dmc0, TDirectory *dmc1) {

  const bool _mc = jp::ismc;
  const bool _x = (dmc0 && dmc1);
  assert(algo=="PF" || algo=="CALO");
  assert(!(_x && algo=="CALO"));

  // Combined reference spectrum
  TH1D *hpt = 0;
  //if (algo=="PF") hpt = (TH1D*)din0->Get("hpt");
  if (algo=="PF") hpt = (TH1D*)din0->Get("hpt_notrigeff");
  if (algo=="CALO") hpt = (TH1D*)din0->Get("hpt_ak5calo");
  assert(hpt);
  // No CALO for MC yet
  //TH1D *mpt = (_x ? (TH1D*)dmc0->Get("hpt") : 0);
  TH1D *mpt = (_x ? (TH1D*)dmc0->Get("hpt_ak5calo") : 0);

  // Trigger-wise spectra
  vector<string> trigs;
  double c = 1.00;
  trigs.push_back("jt40");
  trigs.push_back("jt80");
  trigs.push_back("jt140");
  trigs.push_back("jt200");
  trigs.push_back("jt260");
  trigs.push_back("jt320");
  //trigs.push_back("jt400");

  // Thresholds for triggers
  map<string, double> trgptmax;
  trgptmax["jt40"] = 133;//114;//144;
  trgptmax["jt80"] = 220;//196;
  trgptmax["jt140"] = 300;//245;
  trgptmax["jt200"] = 395;//362;//330;
  trgptmax["jt260"] = 507;//468;//430 ;
  trgptmax["jt320"] = 967;//592;//548;//507;
  //trgptmax["jt400"] = 967;

  // For trigger ratio table
  _trgrationames["jt40"] = "Jet40/Jet80";
  _trgrationames["jt80"] = "Jet80/Jet140";
  _trgrationames["jt140"] = "Jet140/Jet200";
  _trgrationames["jt200"] = "Jet200/Jet260";
  _trgrationames["jt260"] = "Jet260/Jet320";
  _trgrationames["jt320"] = "%Obsolete";//Jet320/Jet400";
  //_trgrationames["jt400"] = "%Obsolete";

  // Temporary rescaling to match triggers
  map<string, double> scale;
  scale["jt40"] = _mc ? 1. : c;
  scale["jt80"] = _mc ? 1 : c;
  scale["jt140"] = _mc ? 1 : c;
  scale["jt200"] = _mc ? 1. : c;
  scale["jt260"] = _mc ? 1. : c;
  scale["jt320"] = _mc ? 1. : c;
  scale["jt400"] = _mc ? 1. : c;

  // Labels
  map<string, const char*> labs;
  labs["jt40"] = "Jet40";
  labs["jt80"] = "Jet80";
  labs["jt140"] = "Jet140";
  labs["jt200"] = "Jet200";
  labs["jt260"] = "Jet260";
  labs["jt320"] = "Jet320";
  labs["jt400"] = "Jet400";

  // Line styles
  map<string, int> lstyle;
  lstyle["jt40"] = kSolid;
  lstyle["jt80"] = kSolid;
  lstyle["jt140"] = kSolid;
  lstyle["jt200"] = kSolid;
  lstyle["jt260"] = kSolid;
  lstyle["jt320"] = kSolid;
  lstyle["jt400"] = kSolid;

  // Line colors
  map<string, int> lcolor;
  lcolor["jt40"] = kBlack;
  lcolor["jt80"] = kGreen+2;
  lcolor["jt140"] = kRed;
  lcolor["jt200"] = kCyan+2;
  lcolor["jt260"] = kOrange+2;
  lcolor["jt320"] = kMagenta+2;
  lcolor["jt400"] = kBlue+2;

  // Marker styles
  map<string, int> mstyle;
  mstyle["jt40"] = kFullCircle;
  mstyle["jt80"] = kFullSquare;
  mstyle["jt140"] = kOpenSquare;
  mstyle["jt200"] = kOpenTriangleUp;
  mstyle["jt260"] = kFullTriangleDown;
  mstyle["jt320"] = kDiamond;
  mstyle["jt400"] = kFullStar;

  // Marker colors
  map<string, int> mcolor;
  mcolor["jt40"] = kBlack;
  mcolor["jt80"] = kGreen+2;
  mcolor["jt140"] = kRed;
  mcolor["jt200"] = kCyan+2;
  mcolor["jt260"] = kOrange+2;
  mcolor["jt320"] = kMagenta+2;
  mcolor["jt400"] = kBlue+2;

  map<string, TH1D*> hpts;
  map<string, TH1D*> mpts;
  for (unsigned int i = 0; i != trigs.size(); ++i) {

    const string& tt = trigs[i];
    const char* t = tt.c_str();
    TH1D *h = 0;
    //if (algo=="PF") h = (TH1D*)din1->Get(Form("%s/hpt",t));
    if (algo=="PF") h = (TH1D*)din1->Get(Form("%s/hpt_notrigeff",t));
    if (algo=="CALO") h = (TH1D*)din1->Get(Form("%s/hpt_ak5calo",t));
    assert(h);
    if (scale[tt]!=1.) {
      h = (TH1D*)h->Clone(Form("hpt%s",t));
      h->Scale(scale[tt]);
    }
    hpts[t] = h;
    // No CALO available for MC yet
    //TH1D *m = (_x ? (TH1D*)dmc1->Get(Form("%s/hpt",t)) : 0);
    TH1D *m = (_x ? (TH1D*)dmc1->Get(Form("%s/hpt_ak5calo",t)) : 0);
    if (!(m||!_x)) {
      cout << dmc1->GetPath() << endl << flush;
      cout << t << endl << flush;
      assert(m||!_x||string(t)=="jt140u");
    }
    mpts[t] = m;
  } // for i

  gStyle->SetOptStat(0);

  const char *c1name = Form("c1_%d", ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogy();
  c1->SetLogx();

  TH1D *h = new TH1D(Form("trig%d",ibin),"",int(xmax-xmin),xmin,xmax);
  h->GetYaxis()->SetTitle("Uncorrected d^{2}#sigma/dp_{T}dy (pb/GeV)");
  h->GetYaxis()->SetTitleOffset(1.15); // ^{2}
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(1e9);
  h->SetMinimum(1e-3);
  h->DrawClone("AXIS");

  float etamin, etamax;
  assert(sscanf(din0->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.25, 0.88, seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);
  teta->Draw();

  int nt = trigs.size();
  TLegend *leg = new TLegend(0.65,0.93-0.05*nt,0.90,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  TLatex *tjet = new TLatex(0.19, 0.18, "Anti-k_{T} R=0.5");
  tjet->SetTextSize(0.05);
  tjet->SetNDC();
  tjet->Draw();

  cmsPrel(_mc ? 0 : _lumi);

  // Underlay MC for comparison, if available
  if (_x) {

    mpt->SetLineWidth(2);
    mpt->SetFillStyle(1001);
    TLegend *leg1b = new TLegend(0.65, 0.57, 0.90, 0.63,"","brNDC");
    leg1b->SetFillStyle(kNone);
    leg1b->SetBorderSize(0);
    leg1b->SetTextSize(0.05);
    leg1b->AddEntry(mpt,"MC (color)","F");
    leg1b->Draw();

    for (unsigned int i = 0; i != trigs.size(); ++i) {

      const char* t = trigs[i].c_str();
      TH1D *m = mpts[t];
      if (m) {
	m->SetLineStyle(lstyle[t]);
	m->SetLineColor(lcolor[t]);
	m->SetLineWidth(2);
	m->Draw("SAME H");
      }
    } // for i
  } // if _x

  for (unsigned int i = 0; i != trigs.size(); ++i) {

    const char* t = trigs[i].c_str();
    TH1D *h = hpts[t];
    h->SetLineStyle(lstyle[t]);
    h->SetLineColor(lcolor[t]);
    h->SetMarkerStyle(mstyle[t]);
    h->SetMarkerColor(mcolor[t]);
    h->SetLineWidth(2);
    h->Draw("SAME");

    leg->AddEntry(h, labs[t], "P");
  } // for i

  leg->Draw();

  // Save the histograms for redrawing PAS plots
  if (_pas) {
    TDirectory *curdir = gDirectory;
    TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
    for (map<string, TH1D*>::const_iterator it = hpts.begin();
	 it != hpts.end(); ++it) {
      TH1D *h = it->second;
      if (h) h->Write(Form("fig1_%d_%s_hpt",ibin,it->first.c_str()),
		      TObject::kOverwrite);
    }
    fpas->Close();
    curdir->cd();
  } // _pas

  const char *t = type.c_str();
  const char *a = algo.c_str();
  const int k = ibin;
  if (_x) {
    if(_eps) c1->SaveAs(Form("eps/triggerSpectra_DATAMC_Eta%d.eps",k));
    if(_pdf) c1->SaveAs(Form("pdf/triggerSpectra_DATAMC_Eta%d.pdf",k));
    if(_gif) c1->SaveAs(Form("gif/triggerSpectra_DATAMC_Eta%d.gif",k));
    if(_png) c1->SaveAs(Form("png/triggerSpectra_DATAMC_Eta%d.png",k));
  }
  else if (algo!="PF") {
    if(_eps) c1->SaveAs(Form("eps/triggerSpectra_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c1->SaveAs(Form("pdf/triggerSpectra_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c1->SaveAs(Form("gif/triggerSpectra_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c1->SaveAs(Form("png/triggerSpectra_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/triggerSpectra_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/triggerSpectra_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/triggerSpectra_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/triggerSpectra_%s_Eta%d.png",t,k));
  }
  if (gROOT->IsBatch()) delete c1;


  // Ratios of triggers to the combined spectrum
  const char *c2name = Form("c2_%d", ibin);
  TCanvas *c2 = new TCanvas(c2name, c2name, 600, 600);
  c2->SetLogx();

  TH1D *hr = (TH1D*)hpt->Clone(Form("hr%d",ibin));
  hr->Divide(hr,hpt,1,1,"");//B");

  h->SetTitleOffset(1.25);
  h->SetMaximum(1.15);
  h->SetMinimum(0.00);
  if (algo=="PF") h->GetYaxis()->SetTitle("PFJet trigger efficiency");
  if (algo=="CALO") h->GetYaxis()->SetTitle("CaloJet trigger efficiency");
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  TH1D *hc = (TH1D*)h->DrawClone("AXIS");

  TLatex *teta2 = (TLatex*)teta->Clone();
  teta2->SetTextSize(0.05);
  teta2->SetX(0.18);
  teta2->SetY(0.88);
  teta2->Draw();

  cmsPrel(jp::isdt ? _lumi : 0);

  //TLegend *leg2 = new TLegend(0.65, 0.23, 0.90, 0.59,"","brNDC");
  TLegend *leg2 = new TLegend(0.70, 0.23, 0.95, 0.59,"","brNDC");
  leg2->SetFillColor(kWhite);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);

  TF1 *ftrig = new TF1(Form("ftrig%d",ibin),
		       "abs([0])+0.5*(1-abs([0]))*(1+TMath::Erf((x-[1])*[2]))",
		       xmin,xmax);
  ftrig->SetParameters(0.02,25.,0.1);
  double ftrig_pars_arr[6][4] =
    {{0.02, 40., 1./3., 0.1},
     {0.00, 80., 1./5., 0.1},
     {0.00, 140., 1./10., 0.1},
     {0.00, 200., 1./12., 0.1},
     {0.00, 260., 1./17., 0.1},
     {0.00, 320., 1./17., 0.1}};
  //{0.00, 400., 1./17., 0.1}};
  map<string, double*> ftrig_pars;
  ftrig_pars["jt40"] = &ftrig_pars_arr[0][0];
  ftrig_pars["jt80"] = &ftrig_pars_arr[1][0];
  ftrig_pars["jt140"] = &ftrig_pars_arr[2][0];
  ftrig_pars["jt200"] = &ftrig_pars_arr[3][0];
  ftrig_pars["jt260"] = &ftrig_pars_arr[4][0];
  ftrig_pars["jt320"] = &ftrig_pars_arr[5][0];
  //ftrig_pars["jt400"] = &ftrig_pars_arr[6][0];

  //{1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
  // 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  // 507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
  // 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

  /*
  // Move to beginning
  map<string, double> trgptmax;
  trgptmax["jt40"] = 144;
  trgptmax["jt80"] = 196;
  trgptmax["jt140"] = 245;
  trgptmax["jt200"] = 330;
  trgptmax["jt260"] = 430;
  trgptmax["jt320"] = 507;
  trgptmax["jt400"] = 967;
  */

  // Underlay MC for comparison, if available
  if (_x) {

    TLegend *leg2b = new TLegend(0.65, 0.17, 0.90, 0.23,"","brNDC");
    leg2b->SetFillColor(kWhite);
    leg2b->SetBorderSize(0);
    leg2b->SetTextSize(0.05);
    leg2b->AddEntry(mpt,"MC (color)","F");
    leg2b->Draw();

    for (unsigned int i = 0; i != trigs.size(); ++i) {

      if (trigs[i]!="mb") {
	string const& tt = trigs[i];
	bool binom = false;
	const char* t = tt.c_str();
	if (mpts[t]) {
	  TH1D *m = (TH1D*)mpts[t]->Clone(Form("mr%d_%s",ibin,t)); assert(m);
	  m->Divide(m,mpt,1,1, binom ? "B" : "");
	  m->SetLineStyle(lstyle[t]);
	  m->SetLineColor(lcolor[t]);
	  m->Draw("SAME H");
	}
      }
    } // for i
  } // if _x

  vector<TH1D*> vtrig(trigs.size());
  vector<TF1*> vfit(trigs.size());
  for (unsigned int i = 0; i != trigs.size(); ++i) {

    string const& t = trigs[i];
    const char* tr = t.c_str();
    TH1D *h = (TH1D*)hpts[t]->Clone(Form("hr%d_%s",ibin,tr));
    TH1D *h0 = (TH1D*)hpt->Clone(Form("h0%d_%s",ibin,tr));
    // Rebin above turn-on into a single bin
    {
      int i0 = min(h->FindBin(trgptmax[t]), h->GetNbinsX());
      if (_mc) i0 = h->GetNbinsX();
      double sum(0), esum(0), sum0(0), esum0(0);
      for (int i = i0; i != h->GetNbinsX()+1; ++i) {
	sum += h->GetBinContent(i);
	double err = h->GetBinError(i);
	esum += err*err;
	h->SetBinContent(i, 0.);
	h->SetBinError(i, 0.);
	sum0 += h0->GetBinContent(i);
	double err0 = h0->GetBinError(i);
	esum0 += err0*err0;
	h0->SetBinContent(i, 0.);
	h0->SetBinError(i, 0.);
      }
      h->SetBinContent(i0, sum);
      h->SetBinError(i0, sqrt(esum));
      h0->SetBinContent(i0, sum0);
      h0->SetBinError(i0, sqrt(esum));
      bool binom = false;
      h->Divide(h,h0,1,1,binom ? "B" : "");

      // Store trigger plateau ratios
      double eff = h->GetBinContent(i0);
      double err = h->GetBinError(i0);
      _trgratios[t][ibin] = make_pair<double, double>(eff, err);
    }

    h->SetLineStyle(lstyle[t]);
    h->SetLineColor(lcolor[t]);
    h->SetMarkerStyle(mstyle[t]);
    h->SetMarkerColor(mcolor[t]);
    h->Draw("SAME");
    vtrig[i] = h;

    if (trigs[i]!="mb") {
      assert(ftrig_pars[t]);
      ftrig->SetParameters(ftrig_pars[t][0],ftrig_pars[t][1],
			   ftrig_pars[t][2]);
      ftrig->SetRange(xmin, xmax);
      h->Fit(ftrig,"QRNW");
      ftrig->SetLineColor(lcolor[t]);
      ftrig->DrawClone("SAME");
      leg2->AddEntry(h, labs[t], "LP");
      vfit[i] = (TF1*)ftrig->Clone(Form("%s_fit",h->GetName()));
    }
    else {
      TF1 *ftrigmb = new TF1(Form("ftrigmb%d",ibin),
			     Form("[0]+log(x*%1.4g/40.)*[1]", cosh(etamin)),
			     xmin,xmax);
      ftrigmb->SetParameters(0.83,0.01);
      h->Fit(ftrigmb,"QRN");
      ftrigmb->SetLineColor(lcolor[t]);
      if (_mc) ftrigmb->DrawClone("SAME");
      leg2->AddEntry(h, labs[t], "LP");
      cout << Form("MinBias eff = %1.3g + %1.2g [%1.1f-%1.1f]",
		   ftrigmb->GetParameter(0), ftrigmb->GetParameter(1),
		   etamin, etamax) << endl;
      vfit[i] = ftrigmb;
    }
  } // for i

  TArrow *a60 = new TArrow(37,0.96,37,1.04,0.02,"<>");
  a60->SetX1(trgptmax["jt30"]); a60->SetX2(trgptmax["jt30"]);
  a60->SetLineColor(lcolor["jt60"]);//+2);
  a60->Draw();

  TArrow *a80 = (TArrow*)a60->Clone();
  a80->SetX1(trgptmax["jt60"]); a80->SetX2(trgptmax["jt60"]);
  a80->SetLineColor(lcolor["jt80"]);//+2);
  //a80->Draw();

  TArrow *a110 = (TArrow*)a60->Clone();
  //a110->SetX1(trgptmax["jt80"]); a110->SetX2(trgptmax["jt80"]);
  a110->SetX1(trgptmax["jt60"]); a110->SetX2(trgptmax["jt60"]);
  a110->SetLineColor(lcolor["jt110"]);//+2);
  a110->Draw();

  TArrow *a150 = (TArrow*)a60->Clone();
  a150->SetX1(trgptmax["jt110"]); a150->SetX2(trgptmax["jt110"]);
  a150->SetLineColor(lcolor["jt150"]);//+2);
  //a150->Draw();

  TArrow *a190 = (TArrow*)a60->Clone();
  //a190->SetX1(trgptmax["jt150"]); a190->SetX2(trgptmax["jt150"]);
  a190->SetX1(trgptmax["jt110"]); a190->SetX2(trgptmax["jt110"]);
  a190->SetLineColor(lcolor["jt190"]);//+2);
  a190->Draw();

  TArrow *a240 = (TArrow*)a60->Clone();
  a240->SetX1(trgptmax["jt190"]); a240->SetX2(trgptmax["jt190"]);
  a240->SetLineColor(lcolor["jt240"]);//+2);
  a240->Draw();

  TArrow *a300 = (TArrow*)a60->Clone();
  a300->SetX1(trgptmax["jt240"]); a300->SetX2(trgptmax["jt240"]);
  a300->SetLineColor(lcolor["jt300"]);//+2);
  a300->Draw();

  TArrow *a370 = (TArrow*)a60->Clone();
  a370->SetX1(trgptmax["jt300"]); a370->SetX2(trgptmax["jt300"]);
  a370->SetLineColor(lcolor["jt370"]);//+2);
  a370->Draw();

  //leg2->AddEntry(a60,"Turn-on","AL");

  leg2->Draw();

  // Save the histograms for redrawing PAS plots
  if (_pas) {
    TDirectory *curdir = gDirectory;
    TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
    for (unsigned int i = 0; i != trigs.size(); ++i) {

      const char *t = trigs[i].c_str();
      TH1D *h = vtrig[i]; assert(h);
      TGraphAsymmErrors *g = new TGraphAsymmErrors(h);
      // Little bit of cleaning to remove the points above the turn-on
      // these are often zero, sometimes >1 due to uncorrelation from
      // prescales and mess up the graph
      for (int j = g->GetN()-1; j != -1; --j) {
	if (g->GetY()[j]>1 || g->GetY()[j]==0 ||
	    //g->GetX()[j]>trgptthr[trigs[i]])
	    g->GetX()[j]>trgptmax[trigs[i]])
	  g->RemovePoint(j);
      }
      // Also reduce error bars that go beyond 1.00
      for (int j = 0; j != g->GetN(); ++j) {
	if (g->GetY()[j]+g->GetEYhigh()[j]>1.00) {
	  g->SetPointError(j, g->GetEXlow()[j], g->GetEXhigh()[j],
			   g->GetEYlow()[j],1.-g->GetY()[j]);
	}
      }
      g->Write(Form("fig2_%d_%s_rpt",ibin,t),TObject::kOverwrite);
      delete g;
      TF1 *fit = vfit[i]; assert(fit);
      if (trigs[i]=="mb")
	fit->SetParameters(1,0,0,1);
      fit->Write(Form("fig2_%d_%s_fit",ibin,t),TObject::kOverwrite);
    }
    fpas->Close();
    curdir->cd();
  } // _pas


  if (_x) {
    if(_eps) c2->SaveAs(Form("eps/triggerRatios_DATAMC_Eta%d.eps",k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatios_DATAMC_Eta%d.pdf",k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatios_DATAMC_Eta%d.gif",k));
    if(_png) c2->SaveAs(Form("png/triggerRatios_DATAMC_Eta%d.png",k));
  }
  else if (algo!="PF") {
    if(_eps) c2->SaveAs(Form("eps/triggerRatios_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatios_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatios_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c2->SaveAs(Form("png/triggerRatios_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c2->SaveAs(Form("eps/triggerRatios_%s_Eta%d.eps",t,k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatios_%s_Eta%d.pdf",t,k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatios_%s_Eta%d.gif",t,k));
    if(_png) c2->SaveAs(Form("png/triggerRatios_%s_Eta%d.png",t,k));
  }

  // Zoom up to the turn-on points
  hc->SetMinimum(0.74999);
  hc->SetMaximum(1.10001);//1.15001);//1.10001);
  hc->GetYaxis()->SetTitleOffset(1.2);//1.3);
  c2->Update();

  if (_x) {
    if(_eps) c2->SaveAs(Form("eps/triggerRatioZoom_DATAMC_Eta%d.eps",k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatioZoom_DATAMC_Eta%d.pdf",k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatioZoom_DATAMC_Eta%d.gif",k));
    if(_png) c2->SaveAs(Form("png/triggerRatioZoom_DATAMC_Eta%d.png",k));
  }
  else if (algo!="PF") {
    if(_eps) c2->SaveAs(Form("eps/triggerRatioZoom_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatioZoom_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatioZoom_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c2->SaveAs(Form("png/triggerRatioZoom_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c2->SaveAs(Form("eps/triggerRatioZoom_%s_Eta%d.eps",t,k));
    if(_pdf) c2->SaveAs(Form("pdf/triggerRatioZoom_%s_Eta%d.pdf",t,k));
    if(_gif) c2->SaveAs(Form("gif/triggerRatioZoom_%s_Eta%d.gif",t,k));
    if(_png) c2->SaveAs(Form("png/triggerRatioZoom_%s_Eta%d.png",t,k));
  }

  if (gROOT->IsBatch()) delete c2;

  // Print out all trigger ratios
  if (_trgratios.begin()->second.size()==6 && !_x) {
    cout << "**********************************" << endl;
    cout << "Results for " << algo << " jets" << endl;
    for (map<string, map<int, pair<double, double> > >::const_iterator
           it = _trgratios.begin(); it != _trgratios.end(); ++it) {
      //cout << it->first << endl << "    ";
      cout << _trgrationames[it->first] << endl << "    ";
      for (map<int, pair<double, double> >::const_iterator
             jt = it->second.begin(); jt != it->second.end(); ++jt) {
        double eff = jt->second.first;
        double err = jt->second.second;
        if (eff-1>2.*err)
          cout << Form(" & {\\color{blue} \\bf %1.3f $\\pm$ %1.3f}",
                       jt->second.first, jt->second.second);
        else if (eff-1>1.*err)
          cout << Form(" & {\\color{blue} %1.3f $\\pm$ %1.3f}",
                       jt->second.first, jt->second.second);
        else if (1-eff>2.*err)
          cout << Form(" & {\\color{red} \\bf %1.3f $\\pm$ %1.3f}",
                       jt->second.first, jt->second.second);
        else if (1-eff>1.*err)
          cout << Form(" & {\\color{red} %1.3f $\\pm$ %1.3f}",
                       jt->second.first, jt->second.second);
        else if (fabs(eff-1)>err)
          cout << Form(" & {\\bf %1.3f $\\pm$ %1.3f}",
                       jt->second.first, jt->second.second);
        else
          cout << Form(" & %1.3f $\\pm$ %1.3f",
                       jt->second.first, jt->second.second);
      } // for jt
      cout << " \\\\" << endl;
    } // for it

    cout << "**********************************" << endl;
  } // _trgratios

} // drawTriggerRatioBin


void drawUnfoldingBin(TDirectory *din, int ibin, string type, string algo);

void drawUnfolding(string type, string algo) {

  assert(algo=="PF" || algo=="CALO");

  //TFile *fin = new TFile(Form("output-%s-3a.root",type.c_str()),"READ");
  TFile *fin = new TFile(Form("output-%s-3b.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")) {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;

      // Process subdirectory
      drawUnfoldingBin(din, ibin++, type, algo);
    } // inherits TDirectory
  } // while

} // drawUnfolding

void drawUnfoldingBin(TDirectory *din, int ibin, string type, string algo) {

  gStyle->SetOptStat(0);

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.35,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.045);//5);

  if (etamin>2.7 || fabs(etamax-1.3)<0.2) return; // patch missing HF unfolding

  TH1D *hpt = (TH1D*)din->Get("hpt"); //assert(hpt);
  if (!hpt) { cout << din->GetName() << endl << flush; assert(false); }
  TGraphErrors *gpt = (TGraphErrors*)din->Get("gpt"); assert(gpt);
  TGraphErrors *gfold = (TGraphErrors*)din->Get("gfold"); assert(gfold);
  TGraphErrors *gratio = (TGraphErrors*)din->Get("gratio"); assert(gratio);
  TF1 *fs = (TF1*)din->Get("fs"); assert(fs);

  if (algo=="CALO") {
    hpt = (TH1D*)din->Get("hpt_ak5calo"); assert(hpt);
    gpt = (TGraphErrors*)din->Get("gpt_ak5calo"); assert(gpt);
    gfold = (TGraphErrors*)din->Get("gfold_ak5calo"); assert(gfold);
    gratio = (TGraphErrors*)din->Get("gratio_ak5calo"); assert(gratio);
    fs = (TF1*)din->Get("fs_ak5calo"); assert(fs);
  }

  // b-tagging specific results
  TGraphErrors *gpt_b = (TGraphErrors*)din->Get("gpt_b");
  TGraphErrors *gfold_b = (TGraphErrors*)din->Get("gfold_b");
  TGraphErrors *gratio_b = (TGraphErrors*)din->Get("gratio_b");
  TF1 *fs_b = (TF1*)din->Get("fs_b");

  const char *schi2 = Form("#chi^{2} / NDF = %1.1f / %d",
			   fs->GetChisquare(), fs->GetNDF());
  TLatex *tchi2 = new TLatex(0.35, 0.82, schi2);
  tchi2->SetNDC();
  tchi2->SetTextSize(0.045);//5);

  const char *schi2b = (fs_b ? Form("#chi^{2} / NDF = %1.1f / %d (b-jets)",
				    fs_b->GetChisquare(), fs_b->GetNDF()) : "");
  TLatex *tchi2b = new TLatex(0.35, 0.76, schi2b);
  tchi2b->SetNDC();
  tchi2b->SetTextSize(0.045);//5);
  tchi2b->SetTextColor(kBlue);

  TLatex *talgo = new TLatex(0.35, 0.76, (algo+" jets").c_str());
  talgo->SetNDC();
  talgo->SetTextSize(0.045);
  talgo->SetTextColor(kBlack);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogy();
  c1->SetLogx();

  TH1D *h = new TH1D(Form("hansatzfit%d",ibin),"",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
  h->GetYaxis()->SetTitleOffset(1.15); // for ^{2}
  h->SetMinimum(1e-4);//1e-2);
  h->SetMaximum(1e10);
  h->Draw("AXIS");

  fs->SetRange(xmin, fs->GetX(h->GetMinimum()));
  if (fs_b) fs_b->SetRange(xmin, fs_b->GetX(h->GetMinimum()));

  gpt->SetMarkerStyle(kFullCircle);
  gpt->Draw("SAMEP");
  fs->Draw("SAME");

  teta->Draw();
  tchi2->Draw();
  cmsPrel(_lumi);

  const char *a = algo.c_str();
  const char *t = type.c_str();
  const int k = ibin;
  if (algo!="PF") {
    talgo->Draw("SAME");
    if(_eps) c1->SaveAs(Form("eps/unfold_ansatzFits_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c1->SaveAs(Form("pdf/unfold_ansatzFits_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c1->SaveAs(Form("gif/unfold_ansatzFits_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c1->SaveAs(Form("png/unfold_ansatzFits_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/unfold_ansatzFits_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/unfold_ansatzFits_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/unfold_ansatzFits_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/unfold_ansatzFits_%s_Eta%d.png",t,k));
  }

  // Add b-jet spectrum on top
  if (gpt_b && fs_b) {

    if (fs_b->GetNDF()!=0) {
      gpt_b->SetMarkerStyle(kOpenCircle);
      gpt_b->SetMarkerColor(kBlue);
      gpt_b->SetLineColor(kBlue);
      fs_b->SetLineColor(kBlue);
      gpt_b->Draw("SAMEP");
      fs_b->Draw("SAME");

      tchi2b->Draw();
    }

    if (algo=="PF") {
      if(_eps) c1->SaveAs(Form("eps/unfoldPlusB_ansatzFits_%s_Eta%d.eps",t,k));
      if(_eps) c1->SaveAs(Form("pdf/unfoldPlusB_ansatzFits_%s_Eta%d.pdf",t,k));
      if(_eps) c1->SaveAs(Form("gif/unfoldPlusB_ansatzFits_%s_Eta%d.gif",t,k));
      if(_eps) c1->SaveAs(Form("png/unfoldPlusB_ansatzFits_%s_Eta%d.png",t,k));
    }
  }

  // Clean up afterwards
  if (gROOT->IsBatch()) delete c1;


  const char *c2name = Form("c2_%d",ibin);
  TCanvas *c2 = new TCanvas(c2name, c2name, 600, 600);
  c2->SetLogx();

  TH1D *h2 = new TH1D("h2","",int(xmax-xmin),xmin,xmax);
  h2->GetXaxis()->SetTitle("p_{T} (GeV)");
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->GetYaxis()->SetTitle("Unfolding correction");
  h2->GetYaxis()->SetRangeUser(0.60,1.15);
  h2->Draw("AXIS");

  gfold->SetMarkerStyle(kFullCircle);
  gfold->Draw("SAMEP");

  teta->Draw();
  cmsPrel(jp::isdt ? _lumi : 0);

  if (algo!="PF") {
    talgo->Draw();
    if(_eps) c2->SaveAs(Form("eps/unfold_correction_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c2->SaveAs(Form("pdf/unfold_correction_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c2->SaveAs(Form("gif/unfold_correction_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c2->SaveAs(Form("png/unfold_correction_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c2->SaveAs(Form("eps/unfold_correction_%s_Eta%d.eps",t,k));
    if(_pdf) c2->SaveAs(Form("pdf/unfold_correction_%s_Eta%d.pdf",t,k));
    if(_gif) c2->SaveAs(Form("gif/unfold_correction_%s_Eta%d.gif",t,k));
    if(_png) c2->SaveAs(Form("png/unfold_correction_%s_Eta%d.png",t,k));
  }

  // Add b-jet unfolding on top
  if (gfold_b) {

    if (gfold_b->GetN()!=0) {
      gfold_b->SetMarkerStyle(kOpenCircle);
      gfold_b->SetMarkerColor(kBlue);
      gfold_b->SetLineColor(kBlue);

      gfold_b->Draw("SAMEP");

      TLegend *leg = new TLegend(0.35,0.74,0.65,0.86,"","brNDC");
      leg->SetFillStyle(kNone);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(gfold,"Inclusive jets","P");
      leg->AddEntry(gfold_b,"b-jets","P");
      leg->Draw();
    }

    if (algo=="PF") {
      if(_eps) c2->SaveAs(Form("eps/unfoldPlusB_correction_%s_Eta%d.eps",t,k));
      if(_pdf) c2->SaveAs(Form("pdf/unfoldPlusB_correction_%s_Eta%d.pdf",t,k));
      if(_gif) c2->SaveAs(Form("gif/unfoldPlusB_correction_%s_Eta%d.gif",t,k));
      if(_png) c2->SaveAs(Form("png/unfoldPlusB_correction_%s_Eta%d.png",t,k));
    }
  }

  // Clean up afterwards
  if (gROOT->IsBatch()) delete c2;


  const char *c3name = Form("c3_%d",ibin);
  TCanvas *c3 = new TCanvas(c3name, c3name, 600, 600);
  c3->SetLogx();

  TH1D *h3 = new TH1D("h3","",int(xmax-xmin),xmin,xmax);
  h3->GetXaxis()->SetTitle("p_{T} (GeV)");
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  h3->GetXaxis()->SetRangeUser(xmin, xmax);
  h3->GetYaxis()->SetTitle("Ratio to ansatz");
  h3->GetYaxis()->SetRangeUser(0.5,2.);//0.,4.0);
  h3->Draw("AXIS");

  gratio->SetMarkerStyle(kFullCircle);
  gratio->Draw("SAMEP");

  TLine *l1 = new TLine(xmin,1.,xmax,1.);
  l1->Draw();

  teta->Draw();
  tchi2->Draw();
  cmsPrel(jp::isdt ? _lumi : 0);

  if (algo!="PF") {
    talgo->Draw();
    if(_eps) c3->SaveAs(Form("eps/unfold_ansatzRatio_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c3->SaveAs(Form("pdf/unfold_ansatzRatio_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c3->SaveAs(Form("gif/unfold_ansatzRatio_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c3->SaveAs(Form("png/unfold_ansatzRatio_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c3->SaveAs(Form("eps/unfold_ansatzRatio_%s_Eta%d.eps",t,k));
    if(_pdf) c3->SaveAs(Form("pdf/unfold_ansatzRatio_%s_Eta%d.pdf",t,k));
    if(_gif) c3->SaveAs(Form("gif/unfold_ansatzRatio_%s_Eta%d.gif",t,k));
    if(_png) c3->SaveAs(Form("png/unfold_ansatzRatio_%s_Eta%d.png",t,k));
  }

  // Zoom up a bit on the low pT range
  h3->GetYaxis()->SetRangeUser(xmin,100.);
  h3->SetMinimum(0.80);
  h3->SetMaximum(1.4);
  if (algo!="PF") {
    talgo->Draw();
    if(_eps) c3->SaveAs(Form("eps/unfold_ansatzZoom_%s_%s_Eta%d.eps",a,t,k));
    if(_pdf) c3->SaveAs(Form("pdf/unfold_ansatzZoom_%s_%s_Eta%d.pdf",a,t,k));
    if(_gif) c3->SaveAs(Form("gif/unfold_ansatzZoom_%s_%s_Eta%d.gif",a,t,k));
    if(_png) c3->SaveAs(Form("png/unfold_ansatzZoom_%s_%s_Eta%d.png",a,t,k));
  }
  else {
    if(_eps) c3->SaveAs(Form("eps/unfold_ansatzZoom_%s_Eta%d.eps",t,k));
    if(_pdf) c3->SaveAs(Form("pdf/unfold_ansatzZoom_%s_Eta%d.pdf",t,k));
    if(_gif) c3->SaveAs(Form("gif/unfold_ansatzZoom_%s_Eta%d.gif",t,k));
    if(_png) c3->SaveAs(Form("png/unfold_ansatzZoom_%s_Eta%d.png",t,k));
  }
  // Back to previous settings
  h3->GetYaxis()->SetRangeUser(xmin,xmax);
  h3->SetMinimum(0.5);//0.);
  h3->SetMaximum(2.0);//4.0);

  // Do the same plots for b-jets
  if (gratio_b) {

    if (gratio_b->GetN()!=0) {
      gratio_b->SetMarkerStyle(kOpenCircle);
      gratio_b->SetMarkerColor(kBlue);
      gratio_b->SetLineColor(kBlue);

      gratio_b->Draw("SAMEP");
      l1->Draw();

      tchi2b->Draw();
    }

    if (algo=="PF") {
      if(_eps) c3->SaveAs(Form("eps/unfoldPlusB_ansatzRatio_%s_Eta%d.eps",t,k));
      if(_pdf) c3->SaveAs(Form("pdf/unfoldPlusB_ansatzRatio_%s_Eta%d.pdf",t,k));
      if(_gif) c3->SaveAs(Form("gif/unfoldPlusB_ansatzRatio_%s_Eta%d.gif",t,k));
      if(_png) c3->SaveAs(Form("png/unfoldPlusB_ansatzRatio_%s_Eta%d.png",t,k));
    }
  }

  // Clean up afterwards
  if (gROOT->IsBatch()) delete c3;

} // drawUnfoldingBin


// Only one resolution bin for now, need to update
void drawResolution(string type) {

  TDirectory *curdir = gDirectory;

  TFile *fin1 = new TFile(Form("output-%s-3a.root",type.c_str()),"READ");
  assert(fin1 && !fin1->IsZombie());

  TFile *fin2 = new TFile("AbsJEC-MC-V8P-out.root","READ");
  assert(fin2 && !fin2->IsZombie());

  TF1 *fres = (TF1*)fin1->Get("Standard/Eta_0.0-0.5/fres"); assert(fres);
  TGraphErrors *gasc = (TGraphErrors*)fin2->Get("gasc"); assert(gasc);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TH1D *h = new TH1D("h","",45,5.,50.);
  h->SetMaximum(0.40);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("Jet p_{T} resolution #sigma_{p_{T}}/p_{T}");
  h->GetYaxis()->SetTitleOffset(1.2);

  h->Draw("AXIS");

  TLatex *tmc = new TLatex(0.40,0.85,"MinBias MC 900 GeV");
  tmc->SetNDC();
  tmc->SetTextSize(0.05);
  tmc->Draw();

  TLatex *teta = new TLatex(0.40,0.78,"|#eta| < 2.5");
  teta->SetNDC();
  teta->SetTextSize(0.05);
  teta->Draw();

  fres->Draw("SAME");

  gasc->SetMarkerStyle(kFullCircle);
  gasc->Draw("SAMEP");

  const char* t = type.c_str();
  if(_eps) c1->SaveAs(Form("eps/resolution_%s_Eta0.eps",t));
  if(_pdf) c1->SaveAs(Form("pdf/resolution_%s_Eta0.pdf",t));
  if(_gif) c1->SaveAs(Form("gif/resolution_%s_Eta0.gif",t));
  if(_png) c1->SaveAs(Form("png/resolution_%s_Eta0.png",t));
  if (gROOT->IsBatch()) delete c1;

  curdir->cd();
}


void drawMultijetsBin(TDirectory *din, int ibin, string type);

void drawMultijets(string type) {

  //TFile *fin = new TFile(Form("output-%s-2a.root",type.c_str()),"READ");
  TFile *fin = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")) {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;

      // Process subdirectory
      drawMultijetsBin(din, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawMultiJets


void drawMultijetsBin(TDirectory *din, int ibin, string type) {

  gStyle->SetOptStat(0);

  /*
  TH1D *hpt = (TH1D*)din->Get("mb/hpt"); assert(hpt);
  TH1D *hpt1 = (TH1D*)din->Get("mb/hpt1"); assert(hpt1);
  TH1D *hpt2 = (TH1D*)din->Get("mb/hpt2"); assert(hpt2);
  TH1D *hpt3 = (TH1D*)din->Get("mb/hpt3"); assert(hpt3);
  */
  TH1D *hpt = (TH1D*)din->Get("hpt"); assert(hpt);
  TH1D *hpt1 = (TH1D*)din->Get("hpt1"); assert(hpt1);
  TH1D *hpt2 = (TH1D*)din->Get("hpt2"); assert(hpt2);
  TH1D *hpt3 = (TH1D*)din->Get("hpt3"); assert(hpt3);
  TH1D *hpt0 = (TH1D*)hpt->Clone("hpt0");

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.18,0.87,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  // 3: 1st-3rd jets
  // 2: 1sr-2nd jets
  // 1: leading jet only
  hpt3->Add(hpt2);
  hpt3->Add(hpt1);
  hpt2->Add(hpt1);

  hpt0->Divide(hpt0, hpt, 1, 1, "B");
  hpt3->Divide(hpt3, hpt, 1, 1, "B");
  hpt2->Divide(hpt2, hpt, 1, 1, "B");
  hpt1->Divide(hpt1, hpt, 1, 1, "B");

  hpt0->SetFillColor(kBlack);
  hpt0->SetFillStyle(1001);
  hpt3->SetFillColor(kRed);
  hpt3->SetFillStyle(1001);
  hpt2->SetFillColor(kYellow);
  hpt2->SetFillStyle(1001);
  hpt1->SetFillColor(kGreen+2);
  hpt1->SetFillStyle(1001);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.30,0.20,0.60,0.50,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hpt0,"Inclusive jets","F");
  leg->AddEntry(hpt3,"Three leading jets","F");
  leg->AddEntry(hpt2,"Two leading jets","F");
  leg->AddEntry(hpt1,"Leading jet only","F");

  TH1D *h = new TH1D(Form("hmulti%d",ibin),"",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("Ratio of multijet spectra");
  h->SetMaximum(1.15);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  hpt0->Draw("SAME");
  hpt3->Draw("SAME HF");
  hpt2->Draw("SAME HF");
  hpt1->Draw("SAME HF");

  cmsPrel(jp::isdt ? _lumi : 0);
  teta->Draw();
  leg->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if(_eps) c1->SaveAs(Form("eps/multijetRatio_%s_Eta%d.eps",t,k));
  if(_pdf) c1->SaveAs(Form("pdf/multijetRatio_%s_Eta%d.pdf",t,k));
  if(_gif) c1->SaveAs(Form("gif/multijetRatio_%s_Eta%d.gif",t,k));
  if(_png) c1->SaveAs(Form("png/multijetRatio_%s_Eta%d.png",t,k));
  if (gROOT->IsBatch()) delete c1;

}


/*
void drawBFractionBin(TDirectory *din, int ibin, string type,
		      TDirectory *dmc = 0);

void drawBFraction(string type) {

  TFile *fin = new TFile(Form("output-%s-4.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *fmc = new TFile(Form("output-%s-4.root","MC"),"READ");
  assert(fmc && !fmc->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;
  assert(fmc->cd("Standard"));
  TDirectory *dmc0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0, ibin1 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;
      assert(dmc0->cd(obj->GetName()));
      TDirectory *dmc = gDirectory;

      // Process subdirectory
      if (jp::isdt)
	drawBFractionBin(din, ibin1++, type, dmc);
      drawBFractionBin(din, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawBFraction


void drawBFractionBin(TDirectory *din, int ibin, string type,
		      TDirectory *dmc) {

  gStyle->SetOptStat(0);

  bool _x = (din && dmc);

  TH1D *hbf = (TH1D*)din->Get("hbfrac"); assert(hbf);
  TH1D *hbfmc = (dmc ? (TH1D*)dmc->Get("hbfrac") : 0); assert(hbfmc||!_x);

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| #leq %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.35,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.70,0.70,0.90,0.90,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hbf,type.c_str(),"P");

  TH1D *h = new TH1D("h0_bfrac","",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("Fraction of b-jets");
  h->SetMaximum(0.07);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  if (_x) {
    hbfmc->SetMarkerStyle(kOpenCircle);
    hbfmc->Draw("PSAME");
    leg->AddEntry(hbfmc,"MC","P");
  }

  hbf->SetMarkerStyle(kFullCircle);
  hbf->Draw("PSAME");

  teta->Draw();
  cmsPrel(jp::isdt ? _lumi : 0);
  leg->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if (_x) {
    if(_eps) c1->SaveAs(Form("eps/bfraction_DATAMC_Eta%d.eps",k));
    if(_pdf) c1->SaveAs(Form("pdf/bfraction_DATAMC_Eta%d.pdf",k));
    if(_gif) c1->SaveAs(Form("gif/bfraction_DATAMC_Eta%d.gif",k));
    if(_png) c1->SaveAs(Form("png/bfraction_DATAMC_Eta%d.png",k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/bfraction_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/bfraction_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/bfraction_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/bfraction_%s_Eta%d.png",t,k));
  }
  if (gROOT->IsBatch()) delete c1;

}
*/
/*
void drawBTagPurityBin(TDirectory *din, int ibin, string type,
		       TDirectory *dmc = 0);

void drawBTagPurity(string type) {

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *fmc = new TFile(Form("output-%s-1.root","MC"),"READ");
  assert(fmc && !fmc->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;
  assert(fmc->cd("Standard"));
  TDirectory *dmc0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0, ibin1 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;
      assert(dmc0->cd(obj->GetName()));
      TDirectory *dmc = gDirectory;

      // Process subdirectory
      if (jp::isdt)
	drawBTagPurityBin(din, ibin1++, type, dmc);
      drawBTagPurityBin(din, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawBTagPurity


void drawBTagPurityBin(TDirectory *din, int ibin, string type,
		       TDirectory *dmc) {

  gStyle->SetOptStat(0);

  bool _x = (din && dmc);

  TProfile *pfb = static_cast<TProfile*>(din->Get("mb/pfb"); assert(pfb);
  TProfile *pfbmc = (_x ? static_cast<TProfile*>(dmc->Get("mb/pfb") : 0);
  assert(pfbmc||!_x);

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| #leq %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.35,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.70,0.70,0.90,0.90,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(pfb,type.c_str(),"P");

  TH1D *h = new TH1D("h0_fb","",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tagged sample purity");
  h->SetMaximum(1.0);
  h->SetMinimum(0.4);
  h->Draw("AXIS");

  if (_x) {
    pfbmc->SetMarkerStyle(kOpenCircle);
    pfbmc->Draw("SAMEP");
    leg->AddEntry(pfbmc, "MC", "P");
  }

  pfb->SetMarkerStyle(kFullCircle);
  pfb->Draw("PSAME");

  teta->Draw();
  cmsPrel(jp::isdt ? _lumi : 0);
  leg->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if (_x) {
    if(_eps) c1->SaveAs(Form("eps/btagpurity_DATAMC_Eta%d.eps",k));
    if(_pdf) c1->SaveAs(Form("pdf/btagpurity_DATAMC_Eta%d.pdf",k));
    if(_gif) c1->SaveAs(Form("gif/btagpurity_DATAMC_Eta%d.gif",k));
    if(_png) c1->SaveAs(Form("png/btagpurity_DATAMC_Eta%d.png",k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/btagpurity_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/btagpurity_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/btagpurity_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/btagpurity_%s_Eta%d.png",t,k));
  }
  if (gROOT->IsBatch()) delete c1;

}
*/
/*
void drawBTagEfficiencyBin(TDirectory *din, int ibin, string type,
			   TDirectory *dmc = 0);

void drawBTagEfficiency(string type) {

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *fmc = new TFile(Form("output-%s-1.root","MC"),"READ");
  assert(fmc && !fmc->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;
  assert(fmc->cd("Standard"));
  TDirectory *dmc0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0, ibin1 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;
      assert(dmc0->cd(obj->GetName()));
      TDirectory *dmc = gDirectory;

      // Process subdirectory
      if (jp::isdt)
	drawBTagEfficiencyBin(din, ibin1++, type, dmc);
      drawBTagEfficiencyBin(din, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawBTagEfficiency


void drawBTagEfficiencyBin(TDirectory *din, int ibin, string type,
			   TDirectory *dmc) {

  gStyle->SetOptStat(0);

  bool _x = (din && dmc);

  TH1D *peffb = (TH1D*)din->Get("mb/peffb"); assert(peffb);
  TH1D *peffbmc = (_x ? (TH1D*)dmc->Get("mb/peffb") : 0); assert(peffbmc||!_x);

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| #leq %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.35,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.70,0.70,0.90,0.90,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(peffb,type.c_str(),"P");

  TH1D *h = new TH1D("h0_bfrac","",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tagging efficiency");
  h->SetMaximum(1.0);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  if (_x) {
    peffbmc->SetMarkerStyle(kOpenCircle);
    peffbmc->Draw("SAMEP");
    leg->AddEntry(peffbmc,"MC","P");
  }

  peffb->SetMarkerStyle(kFullCircle);
  peffb->Draw("PSAME");

  teta->Draw();
  cmsPrel(jp::isdt ? _lumi : 0);
  leg->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if (_x) {
    if(_eps) c1->SaveAs(Form("eps/btagefficiency_DATAMC_Eta%d.eps",k));
    if(_pdf) c1->SaveAs(Form("pdf/btagefficiency_DATAMC_Eta%d.pdf",k));
    if(_gif) c1->SaveAs(Form("gif/btagefficiency_DATAMC_Eta%d.gif",k));
    if(_png) c1->SaveAs(Form("png/btagefficiency_DATAMC_Eta%d.png",k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/btagefficiency_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/btagefficiency_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/btagefficiency_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/btagefficiency_%s_Eta%d.png",t,k));
  }
  if (gROOT->IsBatch()) delete c1;

}
*/
/*
void drawBTagNormBin(TDirectory *din, int ibin, string type,
		     TDirectory *dmc = 0);

void drawBTagNorm(string type) {

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *fmc = new TFile(Form("output-%s-1.root","MC"),"READ");
  assert(fmc && !fmc->IsZombie());

  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;
  assert(fmc->cd("Standard"));
  TDirectory *dmc0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0, ibin1 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;
      assert(dmc0->cd(obj->GetName()));
      TDirectory *dmc = gDirectory;

      // Process subdirectory
      if (jp::isdt)
	drawBTagNormBin(din, ibin1++, type, dmc);
      drawBTagNormBin(din, ibin++, type);
    } // inherits TDirectory
  } // while

} // drawBTagNorm

void drawBTagNormBin(TDirectory *din, int ibin, string type,
		     TDirectory *dmc) {

  gStyle->SetOptStat(0);

  bool _x = (din && dmc);

  TH1D *hpt = (TH1D*)din->Get("mb/hbpt"); assert(hpt);
  TProfile *peffb = static_cast<TProfile*>(din->Get("mb/peffb"); assert(peffb);
  TProfile *pfb = static_cast<TProfile*>(din->Get("mb/pfb"); assert(pfb);

  TH1D *hptmc = (_x ? (TH1D*)dmc->Get("mb/hbpt") : 0); assert(hptmc||!_x);
  TProfile *peffbmc = (_x ? static_cast<TProfile*>(dmc->Get("mb/peffb") : 0);
  assert(peffbmc||!_x);
  TProfile *pfbmc = (_x ? static_cast<TProfile*>(dmc->Get("mb/pfb") : 0);
  assert(pfbmc||!_x);

  TH1D *hnorm = (TH1D*)hpt->Clone(Form("hnorm_%d",ibin));
  hnorm->Reset();
  assert(peffb->GetNbinsX()==pfb->GetNbinsX());
  assert(hpt->GetNbinsX()==pfb->GetNbinsX());
  for (int i = 1; i != hnorm->GetNbinsX()+1; ++i) {
    double effb = peffb->GetBinContent(i);
    double fb = pfb->GetBinContent(i);
    double y = (effb ? fb / effb : 0);
    double ey = (effb*fb ? tools::oplus(pfb->GetBinError(i)/fb,
					peffb->GetBinError(i)/effb) * y : 0);
    hnorm->SetBinContent(i, y);
    hnorm->SetBinError(i, ey);
  } // for i

  TH1D *hnormmc = (_x ? (TH1D*)hptmc->Clone(Form("hnormmc_%d",ibin)) : 0);
  if (_x ) {
    hnormmc->Reset();
    assert(peffbmc->GetNbinsX()==pfbmc->GetNbinsX());
    assert(hptmc->GetNbinsX()==pfbmc->GetNbinsX());
    for (int i = 1; i != hnormmc->GetNbinsX()+1; ++i) {
      double effb = peffbmc->GetBinContent(i);
      double fb = pfbmc->GetBinContent(i);
      double y = (effb ? fb / effb : 0);
      double ey = (effb*fb ? tools::oplus(pfbmc->GetBinError(i)/fb,
					  peffbmc->GetBinError(i)/effb)*y : 0);
      hnormmc->SetBinContent(i, y);
      hnormmc->SetBinError(i, ey);
    } // for i
  } // _x


  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| #leq %1.2g",etamax) :
		      Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
  TLatex *teta = new TLatex(0.35,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.05);

  const char *c1name = Form("c1_%d",ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.70,0.70,0.90,0.90,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hnorm,type.c_str(),"P");

  TH1D *h = new TH1D("h0_bnorm","",int(xmax-xmin),xmin,xmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tag rate normalization");
  h->SetMaximum(7.0);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  if (_x) {
    hnormmc->SetMarkerStyle(kOpenCircle);
    hnormmc->Draw("SAMEP");
    leg->AddEntry(hnormmc, "MC", "P");
  }

  hnorm->SetMarkerStyle(kFullCircle);
  hnorm->Draw("PSAME");

  teta->Draw();
  cmsPrel(jp::isdt ? 0 : _lumi);
  leg->Draw();

  const char* t = type.c_str();
  const int k = ibin;
  if (_x) {
    if(_eps) c1->SaveAs(Form("eps/btagnorm_DATAMC_Eta%d.eps",k));
    if(_pdf) c1->SaveAs(Form("pdf/btagnorm_DATAMC_Eta%d.pdf",k));
    if(_gif) c1->SaveAs(Form("gif/btagnorm_DATAMC_Eta%d.gif",k));
    if(_png) c1->SaveAs(Form("png/btagnorm_DATAMC_Eta%d.png",k));
  }
  else {
    if(_eps) c1->SaveAs(Form("eps/btagnorm_%s_Eta%d.eps",t,k));
    if(_pdf) c1->SaveAs(Form("pdf/btagnorm_%s_Eta%d.pdf",t,k));
    if(_gif) c1->SaveAs(Form("gif/btagnorm_%s_Eta%d.gif",t,k));
    if(_png) c1->SaveAs(Form("png/btagnorm_%s_Eta%d.png",t,k));
  }
  if (gROOT->IsBatch()) delete c1;

} // draw BTagNormBin
*/
/*
// Redraw the MC truth purity fit from Hauke
void drawMCPurity() {

  TFile *f = new TFile("fits/h1fbVsRefPt_JetY0to0.5_alleta_PAS.root");
  //TFile *f = new TFile("fits/fb_merged_g_f.root");
  assert(f && !f->IsZombie());

  TCanvas *c0 = (TCanvas*)f->Get("h1fbVsRefPt_JetY0to0.5_alleta_PAS;1");
  assert(c0);
  TMultiGraph *m = (TMultiGraph*)c0->FindObject("h1fbVsRefPt_JetY0to0.5_alleta_PAS");
  assert(m);
  TList *l = m->GetListOfGraphs();
  assert(l);


  const int ny = 4;//5;
  vector<TGraphErrors*> vmc(ny);
  for (int i = 0; i != ny; ++i) {
    const char *hname = Form("h1fbVsRefPt_JetY%1.2gto%1.2g", 0.5*i, 0.5*(i+1));
    //const char *hname = Form("btag/h1fbVsRefPt_JetY%1.2gto%1.2g",0.5*i,0.5*(i+1));
    TGraphErrors *g = (TGraphErrors*)l->FindObject(hname);
    //TGraphErrors *g = (TGraphErrors*)f->Get(hname);
    if (!g) {
      cout << hname << endl << flush;
      assert(g);
    }
    vmc[i] = (TGraphErrors*)g->Clone();
  }
  delete c0;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.30,0.19,0.70,0.19+5*0.06,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  double bxmax2 = bxmax;//272;
  TH1D *h = new TH1D("h0_mcpurity","",int(bxmax2-xminpas),xminpas,bxmax2);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tagged sample purity");
  h->SetMaximum(1.0);
  h->SetMinimum(0.3);//0.0);
  h->Draw("AXIS");

  cmsPrel(0);

  for (int i = 0; i != ny; ++i) {

    TGraphErrors *g = vmc[i];
    g->Draw("SAMEPZ");

    TF1 *ff = (TF1*)g->GetListOfFunctions()->First();
    assert(ff);
    ff->SetLineWidth(2);
    //g->GetListOfFunctions();

    const char *seta = (i==0 ? "|y| < 0.5" :
			Form("%1.1f #leq |y| < %1.1f", 0.5*i, 0.5*(i+1)));
    leg->AddEntry(g, seta, "LP");
  }

  leg->Draw();
  h->Draw("AXIS SAME");

  if(_eps) c1->SaveAs("eps/btagpurity_AllEta.eps");
  if(_pdf) c1->SaveAs("pdf/btagpurity_AllEta.pdf");
  if(_gif) c1->SaveAs("gif/btagpurity_AllEta.gif");
  if(_png) c1->SaveAs("png/btagpurity_AllEta.png");

} // drawMCPurity

// Redraw the MC efficiency fit from Hauke
void drawMCEfficiency() {


  TFile *f = new TFile("fits/BTagEff_merged_g_f.root");
  assert(f && !f->IsZombie());

  const int ny = 4;
  vector<TGraphErrors*> vmc(ny);
  for (int i = 0; i != ny; ++i) {
    const char *hname = Form("btag/h1BTagEffVsRefPt_JetY%1.2gto%1.2g",
			     0.5*i,0.5*(i+1));
    TGraphErrors *g = (TGraphErrors*)f->Get(hname);
    if (!g) {
      cout << hname << endl << flush;
      assert(g);
    }
    vmc[i] = (TGraphErrors*)g->Clone();
  }

  const int color[ny] = {kBlue+1, kRed+1, kGreen+2, kMagenta+2};
  const int marker[ny] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare};

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.50,0.19,0.90,0.19+5*0.06,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  //double bxmax2 = 272;//245.;
  TH1D *h = new TH1D("h0_mceff","",int(bxmax-xminpas),xminpas,bxmax);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tagging efficiency");
  h->SetMaximum(0.7);
  h->SetMinimum(0.);
  h->Draw("AXIS");

  cmsPrel(0);

  for (int i = 0; i != ny; ++i) {

    TGraphErrors *g = vmc[i];
    g->SetLineColor(color[i]);
    g->SetMarkerColor(color[i]);
    g->SetMarkerStyle(marker[i]);
    g->Draw("SAMEPZ");

    TF1 *ff = (TF1*)g->GetListOfFunctions()->First();
    assert(ff);
    ff->SetLineColor(color[i]);
    ff->SetLineWidth(2);

    const char *seta = (i==0 ? "|y| < 0.5" :
			Form("%1.1f #leq |y| < %1.1f", 0.5*i, 0.5*(i+1)));
    leg->AddEntry(g, seta, "LP");
  }

  leg->Draw();
  h->Draw("AXIS SAME");

  if(_eps) c1->SaveAs("eps/btagefficiency_AllEta.eps");
  if(_pdf) c1->SaveAs("pdf/btagefficiency_AllEta.pdf");
  if(_gif) c1->SaveAs("gif/btagefficiency_AllEta.gif");
  if(_png) c1->SaveAs("png/btagefficiency_AllEta.png");

} // drawMCEfficiency
*/
/*
// Redraw the template fit from Daniel in template/fb_plotsfix.root
void drawTemplatePurity() {

  //TFile *f = new TFile("template/fb_plotsfix.root","READ");
  //TFile *f = new TFile("template/fitresult_fix.root","READ");
  TFile *f = new TFile("template/fitresult-1_fix.root","READ");
  assert(f && !f->IsZombie());

  //TGraphAsymmErrors *gdata = (TGraphAsymmErrors*)f->Get("hist0");
  //TGraphAsymmErrors *gmc = (TGraphAsymmErrors*)f->Get("truthHist0");
  //assert(gdata);
  //assert(gmc);

  TCanvas *c0 = (TCanvas*)f->Get("c1;1");
  assert(c0);
  c0->Draw();
  TGraphAsymmErrors *gdata = (TGraphAsymmErrors*)c0->FindObject("hist0");
  TGraphAsymmErrors *gmc = (TGraphAsymmErrors*)c0->FindObject("truthHist0");
  //TGraphAsymmErrors *gsys = (TGraphAsymmErrors*)c0->FindObject("syshist0");
  assert(gdata);
  assert(gmc);
  //assert(gsys);
  //gsys = (TGraphAsymmErrors*)gsys->Clone();
  delete c0;

  cout << gdata->GetN() << endl;
  cout << gmc->GetN() << endl;

  // Calculate data / MC scale factor
  TGraphErrors *gr = new TGraphErrors(gdata->GetN());
  assert(gdata->GetN()==gmc->GetN());
  for (int i = 0; i != gr->GetN(); ++i) {
    double x1, y1, x2, y2;
    gmc->GetPoint(i, x1, y1);
    gdata->GetPoint(i, x2, y2);
    assert(fabs(x2-x1) < 0.01*fabs(x2+x1));
    gr->SetPoint(i, 0.5*(x2+x1), y1 ? y2 / y1 : 0.);
    double dy = (y1*y2 ? sqrt(pow(gmc->GetErrorY(i)/y1,2) +
			      pow(gdata->GetErrorY(i)/y2,2)) * y2 / y1 : 0.);
    gr->SetPointError(i, 0.5*fabs(x2-x1), dy);
  } // for i
  TF1 *fconst = new TF1("fconst","[0]",xminpas,100);
  gr->Fit(fconst, "QRN");
  cout << "data / MC scale factor = "  << fconst->GetParameter(0)
       << " +/- " << fconst->GetParError(0)
       << " ("<<fconst->GetXmin()<<","<<fconst->GetXmax()<<")"<< endl;
  cout << "chi2 / NDF = "  << fconst->GetChisquare() << " / "
       << fconst->GetNDF() << endl;
  fconst->SetRange(xminpas, bxmax);
  gr->Fit(fconst, "QRN");
  cout << "data / MC scale factor = "  << fconst->GetParameter(0)
       << " +/- " << fconst->GetParError(0)
       << " ("<<fconst->GetXmin()<<","<<fconst->GetXmax()<<")"<< endl;
  cout << "chi2 / NDF = "  << fconst->GetChisquare() << " / "
       << fconst->GetNDF() << endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogx();

  TLegend *leg = new TLegend(0.18,0.80,0.48,0.92,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  TLatex *teta = new TLatex(0.92,0.88,"|y| < 2.0");
  teta->SetTextAlign(31); // align right
  teta->SetNDC();
  teta->SetTextSize(0.05);

  TLatex *tfb = new TLatex(0.92,0.25, Form("Data / MC = %1.3f #pm %1.3f",
					   fconst->GetParameter(0),
					   fconst->GetParError(0)));
  tfb->SetTextAlign(31); // align right
  tfb->SetNDC();
  tfb->SetTextSize(0.05);
  TLatex *tchi2 = new TLatex(0.92,0.19, Form("#chi^{2} / NDF = %1.1f / %d",
					     fconst->GetChisquare(),
					     fconst->GetNDF()));
  tchi2->SetTextAlign(31); // align right
  tchi2->SetNDC();
  tchi2->SetTextSize(0.05);

  double bxmax2 = bxmax;//272.;
  TH1D *h = new TH1D("h0_template","",int(bxmax2-xminpas),xminpas,bxmax2);
  h->GetXaxis()->SetTitle("p_{T} (GeV)");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetYaxis()->SetTitle("b-tagged sample purity");
  h->SetMaximum(1.0);
  h->SetMinimum(0.3);//0.0);//0.5);
  h->Draw("AXIS");

  cmsPrel(60);//7.7);

  gmc->SetMarkerStyle(kOpenCircle);
  gmc->SetLineWidth(2);
  gmc->SetFillStyle(1001);
  gmc->SetFillColor(kYellow+2);

  gdata->SetMarkerStyle(kFullCircle);
  gdata->SetLineColor(kRed);
  gdata->SetLineWidth(2);
  //gsys->SetFillStyle(1001);
  //gsys->SetFillColor(kYellow);

  gmc->Draw("SAME E2");
  gdata->Draw("SAME P");

  leg->AddEntry(gdata, "Data", "LP");
  leg->AddEntry(gmc, "MC", "F");
  leg->Draw();

  teta->Draw();
  tfb->Draw();
  tchi2->Draw();
  h->Draw("AXIS SAME");

  if(_eps) c1->SaveAs("eps/templatePurity_Eta0_2.eps");
  if(_pdf) c1->SaveAs("pdf/templatePurity_Eta0_2.pdf");
  if(_gif) c1->SaveAs("gif/templatePurity_Eta0_2.gif");
  if(_png) c1->SaveAs("png/templatePurity_Eta0_2.png");

} // drawTemplatePurity

// Redraw the template fit from Daniel in template/fit_ptX_Y_etaA-B.root
void drawTemplate() {

  vector<string> bins;
  bins.push_back("pt37-56_eta0.0-2.0");
  //bins.push_back("pt56-97_eta0.0-2.0");
  bins.push_back("pt84-97_eta0.0-2.0");

  vector<string> trigs;
  trigs.push_back("Jet6U");
  //trigs.push_back("Jet15U");
  trigs.push_back("Jet30U");

  vector<double> nmax;
  nmax.push_back(160);//90);
  nmax.push_back(40);//140);

  vector<string> pts;
  pts.push_back("37 #leq p_{T} < 56 GeV");
  //pts.push_back("56 #leq p_{T} < 97 GeV");
  pts.push_back("84 #leq p_{T} < 97 GeV");

  for (unsigned int ibin = 0; ibin != bins.size(); ++ibin) {

  const char *cs = bins[ibin].c_str();

  TFile *f = new TFile(Form("template/fit_%s.root",cs),"READ");
  assert(f && !f->IsZombie());

  TCanvas *c0 = (TCanvas*)f->Get("c1"); assert(c0);
  TH1D *hdata = (TH1D*)c0->FindObject(Form("data_%s",cs));
  assert(hdata);  hdata->Sumw2(); hdata = (TH1D*)hdata->Clone("data");
  THStack *hs = (THStack*)c0->FindObject("hs");
  assert(hs);
  TList *lh = hs->GetHists();
  TH1D *hb = (TH1D*)lh->FindObject(Form("bTemplate_%s",cs));
  assert(hb); hb = (TH1D*)hb->Clone("hb");
  TH1D *hc = (TH1D*)lh->FindObject(Form("cTemplate_%s",cs));
  assert(hc); hc = (TH1D*)hc->Clone("hc");
  TH1D *hl = (TH1D*)lh->FindObject(Form("liGluTemplate_%s",cs));
  assert(hl); hl = (TH1D*)hl->Clone("hl");
  delete c0;

  // Calculate the agreement between data and theory
  int ndf = -2;
  double chi2 = 0;
  for (int i = 1; i != hdata->GetNbinsX()+1; ++i) {

    double sum = hb->GetBinContent(i) + hc->GetBinContent(i)
      + hl->GetBinContent(i);
    double data = hdata->GetBinContent(i);
    double err = hdata->GetBinError(i);
    chi2 += (data!=0 ? pow(data - sum, 2) / (err*err) : 0.);
    ndf += (data!=0 ? 1 : 0);
  }
  cout << "Chi2 / NDF = " << chi2 << " / " << ndf << endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  TLegend *leg = new TLegend(0.55,0.72,0.92,0.92,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  TLatex *teta = new TLatex(0.18,0.88,"|y| < 2.0");
  teta->SetTextAlign(11); // align left
  teta->SetNDC();
  teta->SetTextSize(0.045);

  TLatex *tpt = new TLatex(0.55,0.65,pts[ibin].c_str());
  tpt->SetTextAlign(11); // align right
  tpt->SetNDC();
  tpt->SetTextSize(0.045);

  TLatex *tchi2 = new TLatex(0.55,0.59, Form("#chi^{2} / NDF = %1.1f / %d",
					     chi2, ndf));
  tchi2->SetNDC();
  tchi2->SetTextSize(0.045);

  int nbins = hdata->GetNbinsX();
  TH1D *h = new TH1D("h0_template","",nbins, hdata->GetBinLowEdge(1),
		     hdata->GetBinLowEdge(nbins+1));
  h->GetXaxis()->SetTitle("Secondary vertex mass (GeV)");
  h->GetYaxis()->SetTitle("Number of jets / 0.25 GeV");
  h->SetMaximum(nmax[ibin]);
  h->SetMinimum(0);
  h->Draw("AXIS");

  cmsPrel(60);//7.7);

  hs->Draw("SAME");

  hdata->SetMarkerSize(1.5);
  hdata->SetMarkerStyle(kFullCircle);
  hdata->SetLineWidth(2);
  hdata->Draw("SAME LP");

  leg->AddEntry(hdata, "Data", "LP");
  leg->AddEntry(hb, "b template", "F");
  leg->AddEntry(hc, "c template", "F");
  leg->AddEntry(hl, "light template", "F");
  leg->Draw();

  teta->Draw();
  tpt->Draw();
  tchi2->Draw();
  h->Draw("AXIS SAME");

  const char *trg = trigs[ibin].c_str();
  if(_eps) c1->SaveAs(Form("eps/templateFit_%s_Eta0_2.eps",trg));
  if(_pdf) c1->SaveAs(Form("pdf/templateFit_%s_Eta0_2.pdf",trg));
  if(_gif) c1->SaveAs(Form("gif/templateFit_%s_Eta0_2.gif",trg));
  if(_png) c1->SaveAs(Form("png/templateFit_%s_Eta0_2.png",trg));
  }

} // drawTemplate
*/
