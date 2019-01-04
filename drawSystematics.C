// Purpose: Draw inclusive jet cross section systematics
// Author:  mikko.voutilainen@cern.ch
// Created: March 23, 2010
// Updated: March 23, 2010
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "tools.h"
#include "tdrstyle_mod.C"

using namespace std;

// pT range
#include "settings.h"
//const double xmin = 15.;//5.;
//const double xmax = 200.;

const int kRightHatch = 3654;
const int kLeftHatch = 3645;
const int kVertHatch = 3699;
const int kHorizHatch = 3600;
const int kCrossHatch = 3644;

void drawSystematicsBin(TDirectory *din, int ibin, string type, bool is38x);
void drawTheoryBin(TDirectory *din, int ibin, string type);
void drawPDFBin(TDirectory *din, TDirectory *dt5, TDirectory *dt7,
		int ibin, string type);


void drawSystematics(string type, bool is38x = false) {

  TFile *fin = new TFile(Form("output-%s-4.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *ft = new TFile(Form("output-%s-5.root",type.c_str()),"READ");
  assert(ft && !ft->IsZombie());
  TFile *fr = new TFile("fastnlo/lhc57.root","READ");
  assert(fr && !fr->IsZombie());

  // For AK5/AK7 PDF comparisons
  TFile *ft5 = new TFile("GRV23_AK5_42x_v4/output-DATA-5.root","READ");
  assert(ft5 && !ft5->IsZombie());
  TFile *ft5x = new TFile("GRV23_AK5_42x_v4/output-DATA-5-forPDF.root","READ");
  assert(ft5x && !ft5x->IsZombie());
  TFile *ft7 = new TFile("GRV23_AK7_42x_v4/output-DATA-5.root","READ");
  assert(ft7 && !ft7->IsZombie());
  TFile *ft7x = new TFile("GRV23_AK7_42x_v4/output-DATA-5-forPDF.root","READ");
  assert(ft7x && !ft7x->IsZombie());

  // Select top categories for JEC uncertainty
  assert(fin->cd("Standard"));
  TDirectory *din0 = gDirectory;
  assert(ft->cd("Standard"));
  TDirectory *dt0 = gDirectory;
  assert(fr->cd("Standard"));
  TDirectory *dr0 = gDirectory;

  assert(ft5->cd("Standard"));
  TDirectory *dt50 = gDirectory;
  assert(ft5x->cd("Standard"));
  TDirectory *dt50x = gDirectory;
  assert(ft7->cd("Standard"));
  TDirectory *dt70 = gDirectory;
  assert(ft7x->cd("Standard"));
  TDirectory *dt70x = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  int ibin = 0, ibin2 = 0, ibin3 = 0, ibin4 = 0, ibin5 = 0, ibin57 = 0;
  int ibinpdf57 = 0;
  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);
    
    // Found a subdirectory
    //if (obj->InheritsFrom("TDirectory")) {
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_0.0-1.3"
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      TDirectory *din = gDirectory;
      assert(dt0->cd(obj->GetName()));
      TDirectory *dt = gDirectory;
      assert(dr0->cd(obj->GetName()));
      TDirectory *dr = gDirectory;

      assert(dt50->cd(obj->GetName()));
      TDirectory *dt5 = gDirectory;
      assert(dt50x->cd(obj->GetName()));
      TDirectory *dt5x = gDirectory;
      assert(dt70->cd(obj->GetName()));
      TDirectory *dt7 = gDirectory;
      assert(dt70x->cd(obj->GetName()));
      TDirectory *dt7x = gDirectory;

      // Process subdirectory
      drawSystematicsBin(din, ibin++, type, is38x);
      drawTheoryBin(dt, ibin5++, type);
      drawTheoryBin(dr, ibin57++, "AK5AK7"); // AK5/AK7

      // need the output-DATA5-forPDF.root for this
      drawPDFBin(dr, dt5x, dt7x, ibinpdf57++, "AK5AK7");
    } // inherits TDirectory
  } // while

  fin->Close();
  fin->Delete();

} // drawSystematics


void drawSystematicsBin(TDirectory *din, int ibin, string type, bool is38x) {

  setTDRStyle();
  gStyle->SetOptStat(0);

  // List systematics to be drawn
  vector<string> sys;
  sys.push_back(is38x ? "tot_38x" : "tot");
  sys.push_back(is38x ? "jec_tot_38x" : "jec_tot");
  //sys.push_back("jec_abs");
  //sys.push_back("jec_rel");
  //sys.push_back("jec_pt");
  sys.push_back("jer");
  sys.push_back("lum");

  // Legend label
  map<string, const char*> lab;
  lab["tot"] = "Total uncertainty";
  lab["tot_38x"] = "Total uncertainty (38X)";
  lab["jec_tot"] = "Jet energy scale";// (%1.2g%%)",tools::oplus(5,0.5*ibin));
  lab["jec_tot_38x"] = "Jet energy scale (38X)";
  lab["jec_abs"] = Form("Absolute p_{T} (%1.1g%%)", 5.);
  lab["jec_rel"] = Form("Relative p_{T} (%1.2g%%)", 1.0*(ibin+1));
  lab["jec_pt"] = Form("Jet p_{T} extrapolation");//(%1.2g%%)", 1.0*(ibin+1));
  lab["jer"] = "Jet p_{T} resolution";// (10%)";
  lab["lum"] = "Luminosity";// (10%)";

  // Line styles for the systematics
  map<string, unsigned int> lstyle;
  lstyle["tot"] = kSolid;
  lstyle["tot_38x"] = kSolid;
  lstyle["jec_tot"] = kSolid;
  lstyle["jec_tot_38x"] = kSolid;
  lstyle["jec_abs"] = kSolid;
  lstyle["jec_rel"] = kDotted;
  lstyle["jec_pt"] = kDotted;
  lstyle["jer"] = kDashed;
  lstyle["lum"] = kSolid;

  // Line colors for the systematics
  map<string, unsigned int> lcolor;
  lcolor["tot"] = kBlack;
  lcolor["tot_38x"] = kBlack;
  lcolor["jec_tot"] = kYellow+2;
  lcolor["jec_tot_38x"] = kYellow+2;
  lcolor["jec_abs"] = kYellow+2;
  lcolor["jec_rel"] = kOrange+2;
  lcolor["jec_pt"] = kOrange+2;
  lcolor["jer"] = kGreen+2;//kBlack;
  lcolor["lum"] = kBlue;

  // Fill styles for the systematics
  map<string, unsigned int> fstyle;
  fstyle["tot"] = 1001;
  fstyle["tot_38x"] = 1001;
  fstyle["jec_tot"] = 1001;
  fstyle["jec_tot_38x"] = 1001;
  fstyle["jec_abs"] = 1001;
  fstyle["jec_rel"] = (ibin<4 ? 1001 : kCrossHatch);//3001;
  fstyle["jec_pt"] = 1001;//(ibin<4 ? 1001 : kCrossHatch);//3001;
  fstyle["jer"] = kNone;//1001;//kNone;
  fstyle["lum"] = kNone;
  gStyle->SetHatchesLineWidth(2);

  // Fill colors for the systematics
  map<string, unsigned int> fcolor;
  fcolor["tot"] = kGray;
  fcolor["tot_38x"] = kGray;
  fcolor["jec_tot"] = kYellow;
  fcolor["jec_tot_38x"] = kYellow;
  fcolor["jec_abs"] = kYellow;
  fcolor["jec_rel"] = kOrange;
  fcolor["jec_pt"] = kOrange;
  fcolor["jer"] = kGreen-3;//kBlack;
  fcolor["lum"] = kBlue;

  const char *c1name = Form("c%d", ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  const char *h1name = Form("h%d", ibin);
  TH1D *h1 = new TH1D(h1name,"",int(xmax-xmin),xmin,xmax);
  h1->SetMaximum(ibin<4 ? 40 : 60);
  h1->GetYaxis()->SetTitle("Uncertainty on d^{2}#sigma/dp_{T}dy (%)");
  h1->GetYaxis()->SetTitleOffset(1.15); // for ^{2}
  h1->GetXaxis()->SetTitle("p_{T} (GeV)");
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.17,0.68,0.60,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  h1->Draw("AXIS");

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g", etamax)
		      : Form("%1.2g #leq |y| < %1.2g", etamin, etamax));
  TLatex *teta = new TLatex(ibin==0 ? 0.73 : 0.68,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.045);
  teta->Draw();

  for (unsigned int i = 0; i != sys.size(); ++i) {
   
    const char *cs = sys[i].c_str();
    TH1D *hsys = (TH1D*)din->Get(Form("h%s_av",cs));
    if (!hsys) cout << "Systematic for " << cs << " not found!" << endl;
    assert(hsys);
    
    hsys->Scale(100.); // percentages
    hsys->SetLineStyle(lstyle[cs]);
    hsys->SetLineColor(lcolor[cs]);
    hsys->SetFillStyle(fstyle[cs]);
    hsys->SetFillColor(fcolor[cs]);
    hsys->SetLineWidth(2);
    hsys->Draw("SAME");

    bool fill = (fstyle[cs]!=kNone);

    leg->AddEntry(hsys, lab[cs], fill ? "F" : "L");
  } // for i

  leg->Draw();
  cmsPrel(_lumi);

  h1->Draw("AXIS SAME"); // other histos were on top of border

  string s = type + (is38x ? "_38X" : "");
  const char *t = s.c_str();
  const char *a = _algo.c_str();
  if(_eps) c1->SaveAs(Form("eps/systTot_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c1->SaveAs(Form("pdf/systTot_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c1->SaveAs(Form("gif/systTot_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c1->SaveAs(Form("png/systTot_%s_%s_Eta%d.png",a,t,ibin));
  if (gROOT->IsBatch()) delete c1;

  ////////////////////////////////////////////////////////////////////

  // Draw two-sided band
  const char *c2name = Form("c%d", ibin);
  TCanvas *c2 = new TCanvas(c2name, c2name, 600, 600);
  c2->SetLogx();
  
  h1->SetMinimum(-60);//ibin<4 ? -40 : -70);
  h1->SetMaximum(100);//ibin<4 ? 80 : 120);
  h1->Draw("AXIS");

  vector<TH1D*> herr(sys.size());
  vector<TH1D*> herr2(sys.size());
  for (unsigned int i = 0; i != sys.size(); ++i) {
   
    const char *cs = sys[i].c_str();
    TH1D *hsup = (TH1D*)din->Get(Form("h%s_pl",cs));
    if (!hsup) cout << "Systematic for " << cs << " plus not found!" << endl;
    assert(hsup);
    TH1D *hsdw = (TH1D*)din->Get(Form("h%s_mn",cs));
    if (!hsdw) cout << "Systematic for " << cs << " minus not found!" << endl;
    assert(hsdw);

    hsup->Scale(100.); // percentages
    hsup->SetLineStyle(lstyle[cs]);
    hsup->SetLineColor(lcolor[cs]);
    hsup->SetFillStyle(fstyle[cs]);
    hsup->SetFillColor(fcolor[cs]);
    hsup->SetLineWidth(2);
    hsup->Draw("SAME");

    hsdw->Scale(100.); // percentages
    hsdw->SetLineStyle(lstyle[cs]);
    hsdw->SetLineColor(lcolor[cs]);
    hsdw->SetFillStyle(fstyle[cs]);
    hsdw->SetFillColor(fcolor[cs]);
    hsdw->SetLineWidth(2);
    hsdw->Draw("SAME");

    //bool fill = (fstyle[cs]!=kNone);
    //leg->AddEntry(hsup, lab[cs], fill ? "F" : "L");
    // Save these for PAS
    const char *cs2 = cs;
    if (sys[i]=="jec_rel") cs2 = "reljec";
    if (sys[i]=="jec_abs") cs2 = "absjec";
    TH1D *h = (TH1D*)hsup->Clone(Form("fig4_%d_h%s",ibin,cs2));
    TH1D *h2 = (TH1D*)hsup->Clone(Form("fig7_%d_h%s",ibin,cs2));
    for (int j = 1; j != h->GetNbinsX()+1; ++j) {
      h->SetBinContent(j, 0.5*(hsup->GetBinContent(j)+hsdw->GetBinContent(j)));
      h->SetBinError(j, 0.5*(hsup->GetBinContent(j)-hsdw->GetBinContent(j)));
      h2->SetBinContent(j, 1+0.01*h->GetBinContent(j));
      h2->SetBinError(j, 0.01*h->GetBinError(j));
    }
    herr[i] = h;
    herr2[i] = h2;
  } // for i

  leg->Draw();
  teta->Draw();
  
  cmsPrel(_lumi);

  h1->Draw("AXIS SAME"); // other histos were on top of border


  // Save the histograms for redrawing PAS plots
  if (_pas) {
    TDirectory *curdir = gDirectory;
    TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
    for (unsigned int i = 0; i != sys.size(); ++i) {

      TH1D *h = herr[i]; assert(h);
      h->Write(h->GetName(),TObject::kOverwrite);
      delete h;
      TH1D *h2 = herr2[i]; assert(h);
      if (sys[i]=="tot" || sys[i]=="jec_rel")
	h2->Write(h2->GetName(),TObject::kOverwrite);
      delete h2;
    }
    fpas->Close();
    curdir->cd();
  } // _pas

  if(_eps) c2->SaveAs(Form("eps/systPAS_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c2->SaveAs(Form("pdf/systPAS_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c2->SaveAs(Form("gif/systPAS_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c2->SaveAs(Form("png/systPAS_%s_%s_Eta%d.png",a,t,ibin));
  if (gROOT->IsBatch()) delete c2;


} // drawSystematicsBin


void drawTheoryBin(TDirectory *din, int ibin, string type) {

  if (ibin>=6) return; // 3.0-3.2 and 3.2-4.7 not implemented yet

  setTDRStyle();
  gStyle->SetOptStat(0);

  // List systematics to be drawn
  vector<string> sys;
  sys.push_back("tot");
  sys.push_back("pdf");
  sys.push_back("npcorr");
  sys.push_back("scale");
  if (type!="AK5AK7") sys.push_back("alphas");
  if (type=="AK5AK7") sys.push_back("mc");

  // List systematics to be drawn
  map<string, const char*> name;
  name["tot"] = "hsys";
  name["pdf"] = "hpdf";
  name["mc"] = "hmc";
  name["scale"] = "hsc";
  name["npcorr"] = "hnp";
  name["alphas"] = "has";

  // Legend label
  map<string, const char*> lab;
  lab["tot"] = "Total uncertainty";
  lab["pdf"] = "PDF (CT10 68%)";//(CTEQ6.6)";
  lab["mc"] = "Statistics";//(CTEQ6.6)";
  lab["scale"] = "Scale (6-point scheme)";
  lab["npcorr"] = "NP (Pythia-Herwig++)";
  lab["alphas"] = "#alpha_{S}";

  // Line styles for the systematics
  map<string, unsigned int> lstyle;
  lstyle["tot"] = kSolid;
  lstyle["pdf"] = kDashed;
  lstyle["mc"] = kDashDotted;
  lstyle["scale"] = kDotted;
  lstyle["npcorr"] = kSolid;
  lstyle["alphas"] = kDashDotted;

  // Line colors for the systematics
  map<string, unsigned int> lcolor;
  lcolor["tot"] = kBlack;
  lcolor["pdf"] = kRed;
  lcolor["mc"] = kBlack;
  lcolor["scale"] = kGreen+2;
  lcolor["npcorr"] = kBlue;
  lcolor["alphas"] = kMagenta+2;

  // Fill styles for the systematics
  map<string, unsigned int> fstyle;
  fstyle["tot"] = 1001;
  fstyle["pdf"] = kNone;
  fstyle["mc"] = kNone;
  fstyle["scale"] = kNone;
  fstyle["npcorr"] = kNone;
  fstyle["alphas"] = kNone;

  // Fill colors for the systematics
  map<string, unsigned int> fcolor;
  fcolor["tot"] = kYellow;
  fcolor["pdf"] = kRed;
  fcolor["mc"] = kBlack;
  fcolor["scale"] = kGreen+2;
  fcolor["npcorr"] = kBlue;
  fcolor["alphas"] = kMagenta;

  const char *c1name = Form("cs%d", ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  const char *h1name = Form("ht%d", ibin);
  double ptmax = 1327;//1000;
  //if (ibin==4) ptmax = 800;
  //if (ibin==5) ptmax = 500;
  TH1D *h1 = new TH1D(h1name,"",int(ptmax-xmin),xmin,ptmax);
  h1->SetMaximum(ibin<4 ? +60 : 100);
  h1->SetMinimum(ibin<4 ? -40 : -80);
  if (type=="AK5AK7") {
    h1->SetMaximum(ibin<4 ? +6 : 10);
    h1->SetMinimum(ibin<4 ? -4 : -8);
  }
  h1->GetYaxis()->SetTitle("Theory uncertainty (%)");
  h1->GetXaxis()->SetTitle("p_{T} (GeV)");
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.38,0.63,0.78,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  h1->Draw("AXIS");

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g", etamax)
		      : Form("%1.2g#leq|y|<%1.2g", etamin, etamax));
  TLatex *teta = new TLatex(0.19,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.045);

  TLatex *ts = new TLatex(0.48,0.18,
			  type=="AK5AK7" ?
			  "Anti-k_{T} R=0.5 / R=0.7" :
			  _algo=="AK7" ?
			  "Anti-k_{T} R=0.7 Jets" :
			  "Anti-k_{T} R=0.5 Jets");
  ts->SetTextSize(0.05);
  ts->SetNDC();

  vector<TH1D*> vup(sys.size());
  vector<TH1D*> vdw(sys.size());
  vector<TH1D*> vsys(sys.size());
  for (unsigned int i = 0; i != sys.size(); ++i) {

    string c = sys[i];
    const char *cs = name[c];
    TH1D *hsysup = (TH1D*)din->Get(Form("%sup",cs));
    TH1D *hsysdw = (TH1D*)din->Get(Form("%sdw",cs));
    if (!hsysup) cout << "Systematic for " << cs << "up not found!" << endl;
    assert(hsysup);
    vup[i] = hsysup;
    if (!hsysdw) cout << "Systematic for " << cs << "dw not found!" << endl;
    assert(hsysdw);
    vdw[i] = hsysdw;
   
    hsysup->Scale(100.); // percentages
    hsysup->SetLineStyle(lstyle[c]);
    hsysup->SetLineColor(lcolor[c]);
    hsysup->SetFillStyle(fstyle[c]);
    hsysup->SetFillColor(fcolor[c]);
    hsysup->SetLineWidth(2);
    hsysup->Draw("SAME");

    hsysdw->Scale(100.); // percentages
    hsysdw->SetLineStyle(lstyle[c]);
    hsysdw->SetLineColor(lcolor[c]);
    hsysdw->SetFillStyle(fstyle[c]);
    hsysdw->SetFillColor(fcolor[c]);
    hsysdw->SetLineWidth(2);
    hsysdw->Draw("SAME");

    bool fill = (fstyle[c]!=kNone);

    leg->AddEntry(hsysup, lab[c], fill ? "F" : "L");

    // Save these for PAS
    TH1D *h = (TH1D*)hsysup->Clone(Form("fig7_%d_h%s",ibin,c.c_str()));
    for (int j = 1; j != h->GetNbinsX()+1; ++j) {
      h->SetBinContent(j, 1 + 0.01*0.5*(hsysup->GetBinContent(j)
					+hsysdw->GetBinContent(j)));
      h->SetBinError(j, 0.01*0.5*(hsysup->GetBinContent(j)
				  -hsysdw->GetBinContent(j)));
    }
    vsys[i] = h;
  } // for i

  // Hide the biased pThat<15 GeV region (pt<20 GeV or so)
  TBox *b = new TBox(10,h1->GetMinimum(),21.,h1->GetMaximum());
  b->SetFillStyle(kCrossHatch);
  b->SetFillColor(kGray);
  gStyle->SetHatchesLineWidth(2);
  b->Draw();

  teta->Draw();
  leg->Draw();

  h1->Draw("AXIS SAME"); // other histos were on top of border

  // Save the histograms for redrawing PAS plots
  if (_pas) {
    TDirectory *curdir = gDirectory;
    TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
    for (unsigned int i = 0; i != sys.size(); ++i) {

      TH1D *h = vsys[i]; assert(h);
      if (sys[i]=="tot")
	h->Write(TString(h->GetName()).ReplaceAll("htot","hth"),
		 TObject::kOverwrite);
      delete h;
    }
    fpas->Close();
    curdir->cd();
  } // _pas

  const char *t = (type=="AK5AK7" ? "DATA" : type.c_str());
  const char *a = (type=="AK5AK7" ? type.c_str() : _algo.c_str());
  if(_eps) c1->SaveAs(Form("eps/theory_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c1->SaveAs(Form("pdf/theory_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c1->SaveAs(Form("gif/theory_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c1->SaveAs(Form("png/theory_%s_%s_Eta%d.png",a,t,ibin));

  //////////////////////////////////////////////////////////////
  // for PAS, cut out pT<30 GeV
  h1->Draw("AXIS");
  //h1->GetXaxis()->SetRangeUser(28,xmax);
  h1->GetXaxis()->SetRangeUser(xminpas,xmax);
  //h1->SetMinimum(-30);
  //h1->SetMaximum(35);
  h1->SetMinimum(type=="AK5AK7" ? -4 : -39.999);//-55);
  h1->SetMaximum(type=="AK5AK7" ? 6 : 59.999);
  for (unsigned int i = 0; i != vup.size(); ++i) {

    double vptmax[6] = {1327, 1327, 1101, 846, 592, 362};//395};
    double maxpt = vptmax[ibin];
    //if (type=="AK5AK7") 
    vup[i]->GetXaxis()->SetRangeUser(xmin, maxpt);
    //if (type=="AK5AK7") 
    vdw[i]->GetXaxis()->SetRangeUser(xmin, maxpt);
    vup[i]->Draw("SAME");
    vdw[i]->Draw("SAME");
  } // for i
  teta->Draw();
  leg->Draw();
  h1->Draw("AXIS SAME");

  //cmsPrel(0); // not for PAS?
  ts->Draw();

  if(_eps) c1->SaveAs(Form("eps/theoryPAS_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c1->SaveAs(Form("pdf/theoryPAS_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c1->SaveAs(Form("gif/theoryPAS_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c1->SaveAs(Form("png/theoryPAS_%s_%s_Eta%d.png",a,t,ibin));

  if (gROOT->IsBatch()) delete c1;


  if (type=="AK5AK7") {
    
    TDirectory *curdir = gDirectory;
    TFile *fout = new TFile("ak5ak7ratio.root","UPDATE");

    const double vptmax[6] = {1327, 1327, 1101, 846, 592, 362};
    TH1D *hpt = vup[0];
    vector<double> x;
    for (int i = 0; i != hpt->GetNbinsX()+1; ++i) {
      double pt = hpt->GetBinLowEdge(i);
      if (pt>=xmin57 && pt<=vptmax[ibin])
	x.push_back(pt);
    }
    TH1D *href = new TH1D(Form("href_%d",ibin),"",x.size()-1,&x[0]);

    for (unsigned int isrc = 1; isrc != vup.size(); ++isrc) {
      TH1D *hup = tools::Rebin(vup[isrc],href);
      hup->SetTitle(sys[isrc].c_str());
      hup->SetName(Form("nloxnp_sysup%d_eta%d",isrc-1,ibin));
      hup->Write(Form("nloxnp_sysup%d_eta%d",isrc-1,ibin));
      TH1D *hdw = tools::Rebin(vdw[isrc],href);
      hdw->SetTitle(sys[isrc].c_str());
      hdw->SetName(Form("nloxnp_sysdw%d_eta%d",isrc-1,ibin));
      hdw->Write(Form("nloxnp_sysdw%d_eta%d",isrc-1,ibin));
    }

    fout->Close();
    curdir->cd();
  }
  
} // drawTheoryBin


void drawPDFBin(TDirectory *din, TDirectory *dt5, TDirectory *dt7,
		int ibin, string type) {

  if (ibin>=6) return; // 3.0-3.2 and 3.2-4.7 not implemented yet

  setTDRStyle();
  gStyle->SetOptStat(0);

  // List systematics to be drawn
  vector<string> sys;
  sys.push_back("stat");
  sys.push_back("pdf");
  sys.push_back("ct10");
  sys.push_back("mstw");
  sys.push_back("nnpdf");
  sys.push_back("hera");
  assert(type=="AK5AK7");

  // List systematics to be drawn
  map<string, const char*> name;
  name["stat"] = "hmc";
  name["pdf"] = "hpdf";
  name["ct10"] = "hnlo_ct10"; // by Gregory Soyez (ct10k by Klaus Rabbertz
  name["mstw"] = "hnlo_mstw";
  name["nnpdf"] = "hnlo_nnpdf";
  name["hera"] = "hnlo_hera";

  // Legend label
  map<string, const char*> lab;
  lab["stat"] = "Statistics";
  lab["pdf"] = "PDF (CT10 68%)";
  lab["ct10"] = "CT10 (G. Soyez)"; // by Gregory Soyez (central value)
  lab["mstw"] = "MSTW2008";
  lab["nnpdf"] = "NNPDF2.1";
  lab["hera"] = "HERA1.0"; // only 1.0 (not 1.5) available for both AK5 and AK7

  // Line styles for the systematics
  map<string, unsigned int> lstyle;
  lstyle["mc"] = kSolid;//kDashed;
  lstyle["pdf"] = kDashed;
  lstyle["mc"] = kDashDotted;
  lstyle["scale"] = kDotted;
  lstyle["npcorr"] = kSolid;
  lstyle["alphas"] = kDashDotted;

  // Line colors for the systematics
  map<string, unsigned int> lcolor;
  lcolor["stat"] = kGray+1;
  lcolor["pdf"] = kYellow+1;
  lcolor["ct10"] = kBlack;
  lcolor["mstw"] = kGreen+2;
  lcolor["nnpdf"] = kBlue;
  lcolor["hera"] = kMagenta+2;

  // Fill styles for the systematics
  map<string, unsigned int> fstyle;
  fstyle["stat"] = 1001;
  fstyle["pdf"] = 1001;
  fstyle["ct10"] = kNone;
  fstyle["mstw"] = kNone;
  fstyle["nnpdf"] = kNone;
  fstyle["hera"] = kNone;

  // Fill colors for the systematics
  map<string, unsigned int> fcolor;
  fcolor["stat"] = kGray;
  fcolor["pdf"] = kYellow;
  fcolor["ct10"] = kBlack;
  fcolor["mstw"] = kGreen+2;
  fcolor["nnpdf"] = kBlue;
  fcolor["hera"] = kMagenta;

  const char *c1name = Form("cs%d", ibin);
  TCanvas *c1 = new TCanvas(c1name, c1name, 600, 600);
  c1->SetLogx();

  const char *h1name = Form("hp%d", ibin);
  double ptmax = 1327;
  TH1D *h1 = new TH1D(h1name,"",int(ptmax-xmin),xmin,ptmax);
  h1->SetMaximum(+2.00*3);
  h1->SetMinimum(-1.00*3);
  h1->GetYaxis()->SetTitle("Theory / CT10 - 1 (%)");
  h1->GetXaxis()->SetTitle("p_{T} (GeV)");
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.38,0.63,0.78,0.93,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  h1->Draw("AXIS");

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  const char *seta = (etamin==0 ? Form("|y| < %1.2g", etamax)
		      : Form("%1.2g#leq|y|<%1.2g", etamin, etamax));
  TLatex *teta = new TLatex(0.19,0.88,seta);
  teta->SetNDC();
  teta->SetTextSize(0.045);

  TLatex *ts = new TLatex(0.48,0.18,
			  type=="AK5AK7" ?
			  "Anti-k_{T} R=0.5 / R=0.7" :
			  _algo=="AK7" ?
			  "Anti-k_{T} R=0.7 Jets" :
			  "Anti-k_{T} R=0.5 Jets");
  ts->SetTextSize(0.05);
  ts->SetNDC();

  //vector<TH1D*> vup(sys.size());
  //vector<TH1D*> vdw(sys.size());
  TH1D *hsys(0);
  vector<TH1D*> vsys(sys.size());
  for (unsigned int i = 0; i != sys.size(); ++i) {

    string c = sys[i];
    const char *cs = name[c];

    if (c=="pdf" || c=="stat") {
      
      TH1D *hsysup = (TH1D*)din->Get(Form("%sup",cs)); assert(hsysup);
      TH1D *hsysdw = (TH1D*)din->Get(Form("%sdw",cs)); assert(hsysdw);
      hsys = (TH1D*)hsysup->Clone("hsys"); assert(hsys);
      for (int j = 1; j != hsys->GetNbinsX()+1; ++j) {
	hsys->SetBinContent(j, 0.5*(hsysup->GetBinContent(j) +
				    hsysdw->GetBinContent(j)));
	hsys->SetBinError(j, 0.5*fabs(hsysup->GetBinContent(j) -
				      hsysdw->GetBinContent(j)));
      }
    }
    else {//if (c != "pdf") {

      TH1D *hsys5 = (TH1D*)dt5->Get(cs);
      assert(hsys5);
      TH1D *hsys7 = (TH1D*)dt7->Get(cs);
      assert(hsys7);
      TH1D *hsys57 = (TH1D*)hsys7->Clone(Form("%s_ak5ak7",cs));
      //hsys57->Divide(hsys5, hsys7);
      // Divide histograms with different number of bins
      for (int j = 1; j != hsys7->GetNbinsX()+1; ++j) {
	int k = hsys5->FindBin(hsys7->GetBinCenter(j));
	if (hsys7->GetBinContent(j) &&
	    hsys5->GetBinLowEdge(k)==hsys7->GetBinLowEdge(j) &&
	    hsys5->GetBinLowEdge(k+1)==hsys7->GetBinLowEdge(j+1))
	  hsys57->SetBinContent(j, hsys5->GetBinContent(k)
				/ hsys7->GetBinContent(j));
	else
	  hsys57->SetBinContent(j, 0);
	hsys57->SetBinError(j, 0);
      }

      TH1D *hsys50 = (TH1D*)dt5->Get("hnlo_ct10k");
      assert(hsys50);
      TH1D *hsys70 = (TH1D*)dt7->Get("hnlo_ct10k");
      assert(hsys70);
      TH1D *hsys570 = (TH1D*)hsys70->Clone(Form("%s_ak5ak7_0",cs));
      //hsys570->Divide(hsys50, hsys70);
      // Divide histograms with different number of bins
      for (int j = 1; j != hsys70->GetNbinsX()+1; ++j) {
	int k = hsys50->FindBin(hsys70->GetBinCenter(j));
	if (hsys70->GetBinContent(j) &&
	    hsys50->GetBinLowEdge(k)==hsys70->GetBinLowEdge(j) &&
	    hsys50->GetBinLowEdge(k+1)==hsys70->GetBinLowEdge(j+1))
	  hsys570->SetBinContent(j, hsys50->GetBinContent(k)
				 /  hsys70->GetBinContent(j));
	else
	  hsys570->SetBinContent(j, 0);
	hsys570->SetBinError(j, 0);
      }

      hsys = (TH1D*)hsys570->Clone(Form("%s_pdfsys",cs)); assert(hsys);
      //hsys->Divide(hsys57, hsys570);
      //hsys570->Divide(hsys570, hsys570);
      //hsys->Add(hsys, hsys570, 1, -1);
      // Divide histograms with different number of bins
      for (int j = 1; j != hsys570->GetNbinsX()+1; ++j) {
	
	int k = hsys57->FindBin(hsys570->GetBinCenter(j));
	if (hsys570->GetBinContent(j))
	  hsys->SetBinContent(j, 100.*(hsys57->GetBinContent(k)
				       /  hsys570->GetBinContent(j) - 1));
	hsys->SetBinError(j, 0);
	
	// AK5/AK7 ratio as a check
	if (false) {
	h1->GetYaxis()->SetTitle("R(0.5,0.7) - 1 (%)");
	h1->SetMinimum(-40);
	h1->SetMaximum(+40);
	int k1 = hsys5->FindBin(hsys570->GetBinCenter(j));
	int k2 = hsys7->FindBin(hsys570->GetBinCenter(j));
	if (hsys7->GetBinContent(k2) != 0)
	  hsys->SetBinContent(j, 100.*(hsys5->GetBinContent(k1)
				       / hsys7->GetBinContent(k2) - 1));
	hsys->SetBinError(j, 0);
	}
      }
    } // c!="pdf"

    vsys[i] = hsys;

    //hsys->Scale(100.); // percentages
    hsys->SetLineStyle(lstyle[c]);
    hsys->SetLineColor(lcolor[c]);
    hsys->SetFillStyle(fstyle[c]);
    hsys->SetFillColor(fcolor[c]);
    hsys->SetLineWidth(2);

    bool fill = (fstyle[c]!=kNone);
    if (fill) hsys->SetMarkerSize(0);
    hsys->DrawClone(fill ? "SAME E3" : "SAME");


    leg->AddEntry(hsys, lab[c], fill ? "F" : "L");

    /*
    TH1D *hsysup = (TH1D*)din->Get(Form("%sup",cs));
    TH1D *hsysdw = (TH1D*)din->Get(Form("%sdw",cs));
    if (!hsysup) cout << "Systematic for " << cs << "up not found!" << endl;
    assert(hsysup);
    vup[i] = hsysup;
    if (!hsysdw) cout << "Systematic for " << cs << "dw not found!" << endl;
    assert(hsysdw);
    vdw[i] = hsysdw;
   
    hsysup->Scale(100.); // percentages
    hsysup->SetLineStyle(lstyle[c]);
    hsysup->SetLineColor(lcolor[c]);
    hsysup->SetFillStyle(fstyle[c]);
    hsysup->SetFillColor(fcolor[c]);
    hsysup->SetLineWidth(2);
    hsysup->Draw("SAME");

    hsysdw->Scale(100.); // percentages
    hsysdw->SetLineStyle(lstyle[c]);
    hsysdw->SetLineColor(lcolor[c]);
    hsysdw->SetFillStyle(fstyle[c]);
    hsysdw->SetFillColor(fcolor[c]);
    hsysdw->SetLineWidth(2);
    hsysdw->Draw("SAME");

    bool fill = (fstyle[c]!=kNone);

    leg->AddEntry(hsysup, lab[c], fill ? "F" : "L");

    // Save these for PAS
    TH1D *h = (TH1D*)hsysup->Clone(Form("fig7_%d_h%s",ibin,c.c_str()));
    for (int j = 1; j != h->GetNbinsX()+1; ++j) {
      h->SetBinContent(j, 1 + 0.01*0.5*(hsysup->GetBinContent(j)
					+hsysdw->GetBinContent(j)));
      h->SetBinError(j, 0.01*0.5*(hsysup->GetBinContent(j)
				  -hsysdw->GetBinContent(j)));
    }
    vsys[i] = h;
    */
  } // for i

  /*
  // Hide the biased pThat<15 GeV region (pt<20 GeV or so)
  TBox *b = new TBox(10,h1->GetMinimum(),21.,h1->GetMaximum());
  b->SetFillStyle(kCrossHatch);
  b->SetFillColor(kGray);
  gStyle->SetHatchesLineWidth(2);
  b->Draw();
  */

  teta->Draw();
  leg->Draw();

  h1->Draw("AXIS SAME"); // other histos were on top of border

  // Save the histograms for redrawing PAS plots
  /*
  if (_pas) {
    TDirectory *curdir = gDirectory;
    TFile *fpas = new TFile("comparisons/pasPF.root", "UPDATE");
    for (unsigned int i = 0; i != sys.size(); ++i) {

      TH1D *h = vsys[i]; assert(h);
      if (sys[i]=="tot")
	h->Write(TString(h->GetName()).ReplaceAll("htot","hth"),
		 TObject::kOverwrite);
      delete h;
    }
    fpas->Close();
    curdir->cd();
  } // _pas
  */

  const char *t = (type=="AK5AK7" ? "DATA" : type.c_str());
  const char *a = (type=="AK5AK7" ? type.c_str() : _algo.c_str());
  if(_eps) c1->SaveAs(Form("eps/pdf_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c1->SaveAs(Form("pdf/pdf_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c1->SaveAs(Form("gif/pdf_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c1->SaveAs(Form("png/pdf_%s_%s_Eta%d.png",a,t,ibin));

  //////////////////////////////////////////////////////////////
  // for PAS, cut out pT<30 GeV
  h1->Draw("AXIS");
  h1->GetXaxis()->SetRangeUser(xminpas,xmax);
  //h1->SetMinimum(-0.3);
  //h1->SetMaximum(+0.3);
  //for (unsigned int i = 0; i != vup.size(); ++i) {
  for (unsigned int i = 0; i != vsys.size(); ++i) {

    double vptmax[6] = {1327, 13e27, 1101, 846, 592, 362};
    double maxpt = vptmax[ibin];
    /*
    vup[i]->GetXaxis()->SetRangeUser(xmin, maxpt);
    vdw[i]->GetXaxis()->SetRangeUser(xmin, maxpt);
    vup[i]->Draw("SAME");
    vdw[i]->Draw("SAME");
    */
    vsys[i]->GetXaxis()->SetRangeUser(xmin, maxpt);
    bool fill = (fstyle[sys[i]]!=kNone);
    vsys[i]->Draw(fill ? "SAME E3" : "SAME");
  } // for i
  teta->Draw();
  leg->Draw();
  h1->Draw("AXIS SAME");

  ts->Draw();

  if(_eps) c1->SaveAs(Form("eps/pdfPAS_%s_%s_Eta%d.eps",a,t,ibin));
  if(_pdf) c1->SaveAs(Form("pdf/pdfPAS_%s_%s_Eta%d.pdf",a,t,ibin));
  if(_gif) c1->SaveAs(Form("gif/pdfPAS_%s_%s_Eta%d.gif",a,t,ibin));
  if(_png) c1->SaveAs(Form("png/pdfPAS_%s_%s_Eta%d.png",a,t,ibin));

  if (gROOT->IsBatch()) delete c1;
  
} // drawTheoryBin
