// Purpose: Draw reports of trigger rates and prescales by run
// Author:  mikko.voutilainen@cern.ch
// Created: June 5, 2010
// Updated: June 5, 2010
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
#include "TMultiGraph.h"
#include "TLine.h"
#include "TROOT.h"
#include "TEllipse.h"
#include "TArrow.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include "tools.h"
#include "tdrstyle_mod18.C"
#include "settings.h"

using namespace std;

// ranges for drawing
//#include "settings.h"

// for correcting jet rates for PU smearing, set _npv
#include "ptresolution.h"

// Ansatz Kernel
int cnt_a = 0;
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) cout << "." << flush;

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[1]; // rapidity
  double res = res = ptresolution(pt, eta+1e-3) * pt; // PF
  const double s = TMath::Gaus(ptmeas, pt, res, kTRUE);
  const double f = p[2] * exp(p[3]/pt) * pow(pt, p[4])
    * pow(1 - pt*cosh(eta) / 4000., p[5]);

  return (f * s);
}

// Smeared Ansatz
double _epsilon = 1e-12;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  if (!_kernel) _kernel = new TF1("_kernel",smearedAnsatzKernel,1.,1000.,6);

  const double pt = x[0];
  const double eta = p[0];
  double res = ptresolution(pt, eta+1e-3) * pt;
  const double sigma = min(res, 0.30);
  double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  ptmin = max(1.,ptmin); // safety check
  double ptmax = pt / (1. - 2.*sigma); // xmax*(1-2*sigma)=x
  ptmax = min(4000./cosh(eta),ptmax); // safety check
  const double par[6] = {pt, eta, p[1], p[2], p[3], p[4]};
  _kernel->SetParameters(&par[0]);

  //if (p[5]>0 && p[5]<4000./cosh(eta)) ptmin = p[5]; // for smearing matrix
  //if (p[6]>0 && p[6]<4000./cosh(eta)) ptmax = p[6]; // for smearing matrix

  return ( _kernel->Integral(ptmin, ptmax, _epsilon) );
}

double unfold(double pt, double eta, double npv) {

  // ansatz parameters
  const int npar = 5;
  double p[npar] = {0.00, 1.657e14, -17.86, -5.055, 9.480};

  //_npv = npv;
  double fs = smearedAnsatz(&pt, &p[0]);
  double f = p[1] * exp(p[2]/pt) * pow(pt, p[3])
    * pow(1 - pt*cosh(p[0]) / 4000., p[4]);

  if (!(fs>0 && f>0)) {
    cout << "fs: " << fs << ", f: " << f << endl;
    assert(false);
  }

  return (f/fs);
}


void drawRunHistos(string type) {

  TDirectory *curdir = gDirectory;

  setTDRStyle();
  gStyle->SetOptStat(0);

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  ofstream fout("drawRunHistos.txt",ios::out);

  vector<string> etas;
  //etas.push_back("");
  etas.push_back("Barrel");
  //etas.push_back("Transition");
  //etas.push_back("Endcap");
  //etas.push_back("Endcap1");
  //etas.push_back("Endcap2");

  vector<string> types;
  types.push_back("PF");
  //types.push_back("Calo");

  vector<string> trigs;
  trigs.push_back("jt40");
  trigs.push_back("jt60");
  trigs.push_back("jt80");
  trigs.push_back("jt140");
  trigs.push_back("jt200");
  trigs.push_back("jt260");
  trigs.push_back("jt320");
  trigs.push_back("jt400");
  trigs.push_back("jt450");

  map<string, double> thr;
  thr["jt40"] = 49.;
  thr["jt60"] = 74.;
  thr["jt80"] = 97.;
  thr["jt140"] = 174.;
  thr["jt200"] = 245.;
  thr["jt260"] = 330.;
  thr["jt320"] = 430;
  thr["jt400"] = 507.;
  thr["jt450"] = 548.;

  map<string, double> meanpt;
  // from Eta_0.0-1.3 hselpt (upper limit also)
  meanpt["jt40"] = 60;
  meanpt["jt60"] = 80;
  meanpt["jt80"] = 120;
  meanpt["jt140"] = 210;
  meanpt["jt200"] = 290;
  meanpt["jt260"] = 370;
  meanpt["jt320"] = 460;
  meanpt["jt400"] = 530;
  meanpt["jt450"] = 1100;

  map<string, double> scale0;
  scale0["jt40"] =  8.269e+06;
  scale0["jt60"] =  8.269e+06;
  scale0["jt80"] =  3.391e+05;
  scale0["jt140"] =  1.535e+04;
  scale0["jt200"] =  3229;
  scale0["jt260"] =  1193;
  scale0["jt320"] =  0.1*462.4;
  scale0["jt400"] =  0.1*462.4;
  scale0["jt450"] =  0.1*462.4;

  ///////////////////////////////////////////////////

  map<string, map<string, map<string, double> > > scale;

  if (string(jp::algo)=="AK4") {

    /*9fb*/
    scale["jt40"]["Barrel"]["PF"] =  5e7;
    scale["jt60"]["Barrel"]["PF"] =  5e7;
    scale["jt80"]["Barrel"]["PF"] =  2e7;
    scale["jt140"]["Barrel"]["PF"] =  1e6;
    scale["jt200"]["Barrel"]["PF"] =  1e6;
    scale["jt260"]["Barrel"]["PF"] =  1e5;
    scale["jt320"]["Barrel"]["PF"] =  1e5;
    scale["jt400"]["Barrel"]["PF"] =  5e4;
    scale["jt450"]["Barrel"]["PF"] =  1e4;

    scale["jt40"]["Transition"]["PF"] =  3.918e+06;
    scale["jt60"]["Transition"]["PF"] =  3.918e+06;
    scale["jt80"]["Transition"]["PF"] =  1.748e+05;
    scale["jt140"]["Transition"]["PF"] =  7379;
    scale["jt200"]["Transition"]["PF"] =  1415;
    scale["jt260"]["Transition"]["PF"] =  481.1;
    scale["jt320"]["Transition"]["PF"] =  166;
    scale["jt400"]["Transition"]["PF"] =  55.3;
    scale["jt450"]["Transition"]["PF"] =  55.3;

    scale["jt40"]["Endcap"]["PF"] =  1.948e+06;
    scale["jt60"]["Endcap"]["PF"] =  1.948e+06;
    scale["jt80"]["Endcap"]["PF"] =  6.93e+04;
    scale["jt140"]["Endcap"]["PF"] =  1579;
    scale["jt200"]["Endcap"]["PF"] =  162.9;
    scale["jt260"]["Endcap"]["PF"] =  33.11;
    scale["jt320"]["Endcap"]["PF"] =  6.235;
    scale["jt400"]["Endcap"]["PF"] =  0.8653;
    scale["jt450"]["Endcap"]["PF"] =  0.8653;

  }
  if (string(jp::algo)=="AK7") {

    scale["jt40"]["Barrel"]["PF"] =  2.424e+06;
    scale["jt60"]["Barrel"]["PF"] =  2.424e+06;
    scale["jt80"]["Barrel"]["PF"] =  1.681e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.535e+04;
    scale["jt200"]["Barrel"]["PF"] =  3229;
    scale["jt260"]["Barrel"]["PF"] =  744.9;
    scale["jt320"]["Barrel"]["PF"] =  181.4;
    scale["jt400"]["Barrel"]["PF"] =  2.492;
    scale["jt450"]["Barrel"]["PF"] =  2.492;

    scale["jt40"]["Transition"]["PF"] =  6.453e+06;
    scale["jt60"]["Transition"]["PF"] =  6.453e+06;
    scale["jt80"]["Transition"]["PF"] =  2.323e+05;
    scale["jt140"]["Transition"]["PF"] =  9107;
    scale["jt200"]["Transition"]["PF"] =  1706;
    scale["jt260"]["Transition"]["PF"] =  569.7;
    scale["jt320"]["Transition"]["PF"] =  194.8;
    scale["jt400"]["Transition"]["PF"] =  64.02;
    scale["jt450"]["Transition"]["PF"] =  64.02;

    scale["jt40"]["Endcap"]["PF"] =  2.884e+06;
    scale["jt60"]["Endcap"]["PF"] =  2.884e+06;
    scale["jt80"]["Endcap"]["PF"] =  8.321e+04;
    scale["jt140"]["Endcap"]["PF"] =  1790;
    scale["jt200"]["Endcap"]["PF"] =  183;
    scale["jt260"]["Endcap"]["PF"] =  37.78;
    scale["jt320"]["Endcap"]["PF"] =  6.933;
    scale["jt400"]["Endcap"]["PF"] =  0.7181;
    scale["jt450"]["Endcap"]["PF"] =  0.7181;

    scale["jt40"]["Endcap1"]["PF"] =  1.753e+06;
    scale["jt60"]["Endcap1"]["PF"] =  1.753e+06;
    scale["jt80"]["Endcap1"]["PF"] =  5.691e+04;
    scale["jt140"]["Endcap1"]["PF"] =  1487;
    scale["jt200"]["Endcap1"]["PF"] =  168.3;
    scale["jt260"]["Endcap1"]["PF"] =  36.51;
    scale["jt320"]["Endcap1"]["PF"] =  6.858;
    scale["jt400"]["Endcap1"]["PF"] =  0.7181;
    scale["jt450"]["Endcap1"]["PF"] =  0.7181;

    scale["jt40"]["Endcap2"]["PF"] =  1.117e+06;
    scale["jt60"]["Endcap2"]["PF"] =  1.117e+06;
    scale["jt80"]["Endcap2"]["PF"] =  2.57e+04;
    scale["jt140"]["Endcap2"]["PF"] =  263.6;
    scale["jt200"]["Endcap2"]["PF"] =  12.96;
    scale["jt260"]["Endcap2"]["PF"] =  0.888;
    scale["jt320"]["Endcap2"]["PF"] =  0.06549;
    scale["jt400"]["Endcap2"]["PF"] =  0.01417;
    scale["jt450"]["Endcap2"]["PF"] =  0.01417;
  }

  /////////////////////////////////////////////////

  map<string, string> lab;
  lab["jt40"] = "HLT_Jet40";
  lab["jt60"] = "HLT_Jet60";
  lab["jt80"] = "HLT_Jet80";
  lab["jt140"] = "HLT_Jet140";
  lab["jt200"] = "HLT_Jet200";
  lab["jt260"] = "HLT_Jet260";
  lab["jt320"] = "HLT_Jet320";
  lab["jt400"] = "HLT_Jet400";
  lab["jt450"] = "HLT_Jet450";


  TCanvas *c1c = new TCanvas("c1c","c1c",600,600);
  TCanvas *c1d = new TCanvas("c1d","c1d",600,600);

  TCanvas *c1 = new TCanvas("c1","c1",1200,300);
  c1->SetLogy();
  // Adapt to very wide format
  c1->SetLeftMargin(0.05);
  c1->SetRightMargin(0.005);

  TH1D *runs = (TH1D*)fin->Get(Form("Runs%s/runs",etas[0].c_str()));
  assert(runs);
  int nruns = runs->GetNbinsX();
  float lumi;
  assert(sscanf(runs->GetTitle(),"runs %f pb-1",&lumi)==1);

  TH1D *h = new TH1D("h","",nruns,-0.5,nruns-0.5);
  h->SetMinimum(1e-2);
  h->SetMaximum(1e2);
  h->GetYaxis()->SetTitle("Trigger rate per run");
  h->GetYaxis()->SetTitleOffset(0.4);
  vector<int> vrun(nruns);
  for (int i = 1; i != nruns+1; ++i) {
    //if (i%2==1) // Reduce overlap on run numbers
    if (i%3==1) // Reduce overlap on run numbers
      h->GetXaxis()->SetBinLabel(i, Form("%d",int(runs->GetBinContent(i)+0.5)));
    vrun[i-1] = runs->GetBinContent(i);
    if (i>1) assert(vrun[i-2]<vrun[i-1] || vrun[i-1]==0);
  } // for i
  for (int i = nruns+1; i < h->GetNbinsX()+1; ++i) {
    h->GetXaxis()->SetBinLabel(i, "      ");
  }

  int iCa = TMath::BinarySearch(nruns, &vrun[0], 246908); // 2015C first
  int iDb = TMath::BinarySearch(nruns, &vrun[0], 260627); // 2015D last

  int ifirst = iCa;
  int i12b = 0.75*iCa+0.25*iDb;
  int i12c = 0.50*iCa+0.50*iDb;
  int i12d = 0.25*iCa+0.75*iDb;
  int ilast = iDb;

  cout << "ifirst = " << ifirst << endl;
  cout << "i12b = " << i12b << endl;
  cout << "i12c = " << i12c << endl;
  cout << "i12d = " << i12d << endl;
  cout << "ilast = " << ilast << endl;

  map<string, pair<int,int> > ranges;
  ranges["jt40"] = pair<int,int>(ifirst,ilast);
  ranges["jt60"] = pair<int,int>(ifirst,ilast);
  ranges["jt80"] = pair<int,int>(ifirst,ilast);
  ranges["jt140"] = pair<int,int>(ifirst,ilast);
  ranges["jt200"] = pair<int,int>(ifirst,ilast);
  ranges["jt260"] = pair<int,int>(ifirst,ilast);
  ranges["jt320"] = pair<int,int>(ifirst,ilast);
  ranges["jt400"] = pair<int,int>(ifirst,ilast);
  ranges["jt450"] = pair<int,int>(ifirst,ilast);

  const int nc = 3;
  int tricolor[nc] = {kCyan+1, kYellow+1, kMagenta+1};

  const int runbins[] =
    {iCa, i12b, i12c, i12d, ilast};
    // 20/fb divisions, 572 runs, epoch firsts: 1, 63, 209, 389, (569=0)
    //{/*A*/ 0, 31, /*B*/ 63, 111, 160, /*C*/ 209, 269, 329,
    // /*D*/ 389, 449, 509, 569};//572};

  const int nrunbins = sizeof(runbins)/sizeof(runbins[0])-1;

  map<string, int> mcolor;
  mcolor["jt40"] = kBlue+1;
  mcolor["jt80"] = kGreen+1;
  mcolor["jt140"] = kGreen+2;
  mcolor["jt200"] = kYellow+2;
  mcolor["jt260"] = kRed+1;
  mcolor["jt320"] = kRed+2;
  mcolor["jt400"] = kBlack;

  //cmsPrel(lumi);

  TLatex *tjet = new TLatex(0.19,0.19,"Anti-k_{T} R=0.5 PF, |y| < 3.0");
  tjet->SetTextSize(0.045);
  tjet->SetNDC();
  tjet->Draw();

  TArrow *arr = new TArrow(10,1.45,10,1.1,0.01);
  arr->SetLineColor(kRed);
  arr->SetLineWidth(2);
  TArrow *arr2 = new TArrow(10,1.25,10,1.1,0.01);
  arr2->SetLineColor(kMagenta+2);//kBlue);
  arr2->SetLineWidth(2);

  map<int, set<string> > mbadruns;

  // maps for type, eta, trig, values
  map<string, map<string, map<string, pair<double, double> > > > _rate;
  map<string, map<string, map<string, pair<double, double> > > > _ratebefore;
  map<string, map<string, map<string, pair<double, double> > > > _rateafter;

  for (unsigned int itype = 0; itype != types.size(); ++itype) {
  for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {


    map<string, pair<double, double> > _norm;

    string ttc = types[itype];
    const char* tc = ttc.c_str();
    string tte = etas[ieta];
    const char* te = tte.c_str();
    const char *sc = (ttc=="Calo" ? "c" : "");

    TMultiGraph *mrun = new TMultiGraph();
    vector<TGraphErrors*> vmrun(trigs.size());

    c1c->cd();
    TH1D *hc = new TH1D("hc",Form(";Run;Normalized rate in %s",te),
			nruns,-0.5,nruns+0.5);
    hc->SetMinimum(0.8);//tte=="Barrel" ? 0.9 : 0.8);
    hc->SetMaximum(1.4);//tte=="Barrel" ? 1.2 : 1.4);
    hc->DrawClone("AXIS");
    c1d->cd();
    hc->SetMaximum(1.2);
    hc->DrawClone("AXIS");

    c1->cd();

    for (unsigned int j = 0; j != trigs.size(); ++j) {

      // Load central values
      assert(fin->cd(Form("Runs%s",etas[ieta].c_str())));
      TDirectory *d = gDirectory;

      string tt = trigs[j];
      const char* t = tt.c_str();
      TH1D *hn = (TH1D*)d->Get(Form("npvgood_%s",t)); assert(hn);
      TH1D *r0 = (TH1D*)d->Get(Form("r0_%s%s","",t)); assert(r0); // no Calo
      TH1D *r = (TH1D*)d->Get(Form("r_%s%s",sc,t)); assert(r);

      // Correct rates for JER vs NPV
      for (int i = 1; i != r->GetNbinsX()+1; ++i) {
	double eta = 0.;
	if (etas[ieta]=="Transition") eta = 1.0;
	if (etas[ieta]=="Endcap") eta = 2.0;
	if (etas[ieta]=="Endcap1") eta = 2.0;
	if (etas[ieta]=="Endcap2") eta = 2.5;
	double pt = meanpt[t];
	double npv = hn->GetBinContent(i);
	double y = r->GetBinContent(i);
	double y0 = r0->GetBinContent(i);
	_ak7 = (string(jp::algo)=="AK7"); // change PU smearing
	double c = 1;//unfold(pt, eta, npv) / unfold(pt, eta, _npv0);
	r->SetBinContent(i, c * y);
	r0->SetBinContent(i, c * y0);
      }

      // Clone for marking deviating runs
      TH1D *r1 = (TH1D*)r->Clone(Form("r1_%s%s",sc,t));
      TH1D *r2 = (TH1D*)r->Clone(Form("r2_%s%s",sc,t));

      r0->Scale(1./scale0[t]);//[te][tc]);
      r->Scale(1./scale[t][te][tc]);
      r1->Scale(1./scale[t][te][tc]);
      r2->Scale(1./scale[t][te][tc]);


      // Coarse grain vs run
      TF1 *fit = new TF1("fit","[0]",0,350);
      TGraphErrors *grun = new TGraphErrors(0);//nrunbins);
      grun->SetName(Form("grun%s%s_%s",tc,te,t));
      for (int irun = 0; irun != nrunbins; ++irun) {

	double midrun = 0.5*(runbins[irun]+runbins[irun+1]-1);
	double drun = 0.5*(runbins[irun+1]-runbins[irun]);
	fit->SetRange(runbins[irun]-0.5, runbins[irun+1]-0.5);

	r->Fit(fit,"QRN");
	if (fit->GetParError(0)<0.10) {
	  int n = grun->GetN();
	  grun->SetPoint(n, midrun, fit->GetParameter(0));
	  grun->SetPointError(n, drun, fit->GetParError(0));
	}
      } // for irun
      grun->SetMarkerStyle(kOpenCircle);
      mrun->Add(grun);
      vmrun[j] = grun;
      delete fit;

      c1c->cd();
      if (!(ttc=="PF"
	    &&(tt=="jt320"||tt=="jt400"))) grun->SetMarkerSize(0);
      else grun->SetMarkerSize(1.5);
      grun->SetLineColor(mcolor[tt]);
      grun->SetMarkerColor(mcolor[tt]);
      grun->DrawClone("SAMEPz");
      c1d->cd();
      if (tt=="jt320") {
	grun->DrawClone("SAMEPz");

	//TF1 *f1d = new TF1("f1d","[0]",ifirst,its);
	TF1 *f1d = new TF1("f1d","[0]",ifirst,i12b);
	f1d->SetLineStyle(kDashed);
	f1d->SetLineColor(grun->GetLineColor());
	grun->Fit(f1d,"QRN");
	f1d->DrawClone("SAME");
	//f1d->SetRange(i12b,ijump);
	f1d->SetRange(i12b,i12c);
	grun->Fit(f1d,"QRN");
	f1d->DrawClone("SAME");
	//f1d->SetRange(ijump,i12c);
	f1d->SetRange(i12c,i12d);
	grun->Fit(f1d,"QRN");
	f1d->DrawClone("SAME");
	//f1d->SetRange(i12c,ilast);
	f1d->SetRange(i12d,ilast);
	grun->Fit(f1d,"QRN");
	f1d->DrawClone("SAME");
      }
      c1->cd();

      // Find bins that deviate from average
      assert(r1->GetNbinsX()==nruns);
      double ymax = 1.;
      double ymin = 1.;
      for (int k = 1; k != nruns+1; ++k) {

	if (r1->GetBinContent(k) > ymax) ymax = r1->GetBinContent(k);
	if (r1->GetBinContent(k) < ymin &&
	    r1->GetBinContent(k) > 0) ymin = r1->GetBinContent(k);

	if (fabs(r1->GetBinContent(k)-1) < 3.*r1->GetBinError(k)
	    || fabs(r1->GetBinContent(k)-1) < 0.05
	    ) {
	  if (r->GetBinContent(k)!=0)
	    r->SetBinContent(k, max(0.51,min(1.49,r->GetBinContent(k))));
	  r1->SetBinContent(k, 0.);
	  r1->SetBinError(k, 0.);
	}
	else {
	  r->SetBinContent(k, 0.);
	  r->SetBinError(k, 0.);
	  if (r1->GetBinContent(k)!=0) {
	    int run = vrun[k-1];
	    cout << Form( "Run %d trigger %s rate is %1.3g+/-%1.2g"
			  " (%1.3g+/-%1.2g) times normal",
			  run, t,
			  r1->GetBinContent(k), r1->GetBinError(k),
			  r0->GetBinContent(k), r0->GetBinError(k)) << endl;
	    //if (!(string(t)=="mb" && run>137028) &&
	    //!(string(t)=="jt6u") &&
	    //!(string(t)=="jt15u" && run<136066) &&
	    //!(string(t)=="jt30u" && run<136066) &&
	    //!(string(t)=="jt50u" && run<136066)
	    //)
	    //mbadruns[run].insert(t);
	  }
	}
	// Get the really low stats out from cluttering the picture
	double maxerr = 0.10;
	if (r0->GetBinError(k)>maxerr) {
	  r0->SetBinContent(k, 0.);
	  r0->SetBinError(k, 0.);
	}
	if (r->GetBinError(k)>maxerr) {
	  r->SetBinContent(k, 0.);
	  r->SetBinError(k, 0.);
	}
	if (r1->GetBinError(k)>maxerr) {
	  r1->SetBinContent(k, 0.);
	  r1->SetBinError(k, 0.);
	}
	if (r2->GetBinError(k)>maxerr) {
	  r2->SetBinContent(k, 0.);
	  r2->SetBinError(k, 0.);
	}
      } // for k
      cout << "ymax = " << ymax << " ymin = " << ymin << endl;
      ymax = pow(10.,1+int(log10(ymax)));
      ymin = pow(10.,-1+int(log10(ymin)));
      cout << "ymax_round = " << ymax << " ymin_round = " << ymin << endl;
      h->SetMaximum(ymax);
      h->SetMinimum(ymin);
      h->GetYaxis()->SetTitle(Form("%s normalized rate",lab[t].c_str()));

      TF1 *f = new TF1("f","[0]",-1,nruns);
      r->Fit(f,"QRN");
      cout << Form("Normalization for %s: %1.4g",t,f->GetParameter(0)) << endl;
      delete f;

      TLegend *leg = new TLegend(0.73, 0.77, 0.93, 0.92, "", "brNDC");
      leg->SetFillStyle(kNone);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.045);
      if (jp::debug) leg->AddEntry(r0,"Jet p_{t}>20 GeV","L");
      leg->AddEntry(r,Form("Jet p_{T}>%1.0f GeV (%s)",thr[t],t),"L");
      leg->AddEntry(r1,"Flagged runs","PL");

      TH1D *haxis = (TH1D*)h->DrawClone("AXIS");

      r0->SetLineWidth(2);
      if (jp::debug) r0->Draw("SAME");
      r->SetLineColor(kBlue);
      r->Draw("SAME");
      r1->SetLineColor(kRed);
      r1->SetMarkerColor(kRed);
      r1->DrawClone("SAME");
      r1->SetMarkerStyle(kCircle);
      r1->Draw("SAME");

      leg->Draw();
      //tjet->Draw();
      tjet->DrawLatex(0.19,0.19,Form("Anti-k_{T} R=%1.1f PF, %s",
				     string(jp::algo)=="AK4" ? 0.4 : 0.7, te));
      //cmsPrel(lumi);

      if(jp::pdf) c1->SaveAs(Form("pdf/RunHistos_%s_%s_%s.pdf",tc,te,t));

      // Same in linear scale
      c1->SetLogy(0);
      haxis->SetMaximum(1.5);
      haxis->SetMinimum(0.5);

      TLine *lup = new TLine(0,1.1,h->GetNbinsX()-1,1.1);
      lup->SetLineColor(kGreen+2);
      lup->SetLineStyle(kDashed);
      lup->Draw();
      TLine *lmid = new TLine(0,1,h->GetNbinsX()-1,1);
      lmid->SetLineColor(kGreen+2);
      lmid->Draw();
      TLine *ldw = new TLine(0,0.9,h->GetNbinsX()-1,0.9);
      ldw->SetLineColor(kGreen+2);
      ldw->SetLineStyle(kDashed);
      ldw->Draw();


      // Rate before and after tech stop
      TF1 *fpol0 = new TF1("fpol0","[0]",0,0.);
      fpol0->SetParameter(0, 1.);
      //
      fpol0->SetRange(ranges[t].first-0.5, ranges[t].second+0.5);
      r2->Fit(fpol0,"QRN");
      double chi = sqrt(max(1., fpol0->GetChisquare()/fpol0->GetNDF()));
      _rate[tc][te][t] =
	pair<double, double>(fpol0->GetParameter(0),
				  chi*fpol0->GetParError(0));
      _norm[t] =
	pair<double, double>(fpol0->GetParameter(0),
				  chi*fpol0->GetParError(0));
      //
      fpol0->SetRange(ranges[t].first-0.5,
		      //min(ranges[t].second+0.5,its));
		      ranges[t].second+0.5);
      r2->Fit(fpol0,"QRN");
      chi = sqrt(max(1., fpol0->GetChisquare()/fpol0->GetNDF()));
      _ratebefore[tc][te][t] =
	pair<double, double>(fpol0->GetParameter(0),
				  chi*fpol0->GetParError(0));
      fpol0->SetLineColor(kYellow+2);
      //if (ranges[t].first<its) fpol0->DrawClone("SAME");
      fpol0->DrawClone("SAME");
      //
      fpol0->SetRange(//max(its,ranges[t].first-0.5),
		      //min(ranges[t].second+0.5,nruns+0.5));
		      ranges[t].first-0.5,
		      ranges[t].second+0.5);
      r2->Fit(fpol0,"QRN");
      chi = sqrt(max(1., fpol0->GetChisquare()/fpol0->GetNDF()));
      _rateafter[tc][te][t] =
	pair<double, double>(fpol0->GetParameter(0),
				  chi*fpol0->GetParError(0));
      fpol0->SetLineColor(kRed+1);
      //if (ranges[t].second>its) fpol0->DrawClone("SAME");
      fpol0->DrawClone("SAME");

      int ic = -1;
      //
      fpol0->SetLineColor(kRed+1);
      //fpol0->SetRange(ifirst, its);
      fpol0->SetRange(ifirst, i12c);
      r2->Fit(fpol0,"QRN");
      fpol0->DrawClone("SAME");
      //
      //fpol0->SetRange(its, ilast);
      fpol0->SetRange(i12c, ilast);
      r2->Fit(fpol0,"QRN");
      fpol0->DrawClone("SAME");
      //

      // Redraw markers on top
      r->Draw("SAME");
      r1->Draw("SAME");

      //arr2->SetX1(its);
      //arr2->SetX2(its);
      //arr2->DrawClone();
      arr2->SetX1(i12b);
      arr2->SetX2(i12b);
      arr2->DrawClone();
      //arr2->SetX1(ijump);
      //arr2->SetX2(ijump);
      arr2->SetX1(i12c);
      arr2->SetX2(i12c);
      arr2->DrawClone();
      arr2->SetX1(i12d);
      arr2->SetX2(i12d);
      arr2->DrawClone();

      arr->SetX1(ranges[t].first);
      arr->SetX2(ranges[t].first);
      arr->DrawClone();
      arr->SetX1(ranges[t].second);
      arr->SetX2(ranges[t].second);
      arr->DrawClone();

      c1->Update();


      if (jp::pdf) c1->SaveAs(Form("pdf/RunHistos_%s_%s_%s_lin.pdf",tc,te,t));

      c1->SetLogy();

    } // for j


    // Coarse grain vs trigger (and run)
    TF1 *fit = new TF1("mfit","[0]",0,350);
    TGraphErrors *grun = new TGraphErrors(nrunbins);
    for (int irun = 0; irun != nrunbins; ++irun) {

      double midrun = 0.5*(runbins[irun]+runbins[irun+1]-1);
      double drun = 0.5*(runbins[irun+1]-runbins[irun]);
      //fit->SetRange(runbins[irun]+0.5, runbins[irun+1]+0.5);
      fit->SetRange(runbins[irun]-0.5, runbins[irun+1]-0.5);

      mrun->Fit(fit,"QRN");
      double k = sqrt(fit->GetChisquare() / max(1,fit->GetNDF()));
      grun->SetPoint(irun, midrun, fit->GetParameter(0));
      grun->SetPointError(irun, drun, k*fit->GetParError(0));
    } // for irun

    c1c->cd();

    TLine *lc = new TLine();
    double y2 = hc->GetMaximum();
    double y1 = hc->GetMinimum();
    lc->SetLineStyle(kDashed);
    //lc->DrawLine(its,y1,its,y2);
    lc->DrawLine(i12b,y1,i12b,y2);
    //lc->DrawLine(ijump,y1,ijump,y2);
    lc->DrawLine(i12c,y1,i12c,y2);
    lc->SetLineStyle(kDotted);
    //lc->DrawLine(its,y1,its,y2);
    //lc->DrawLine(ijump,y1,ijump,y2);
    lc->DrawLine(i12d,y1,i12d,y2);
    lc->SetLineStyle(kDotted);
    lc->DrawLine(-0.5,1.00,nruns-0.5,1.00);
    c1d->cd();
    lc->SetLineStyle(kDashed);
    //lc->DrawLine(its,y1,its,y2);
    lc->DrawLine(i12b,y1,i12b,y2);
    //lc->DrawLine(ijump,y1,ijump,y2);
    lc->DrawLine(i12c,y1,i12c,y2);
    lc->SetLineStyle(kDotted);
    //lc->DrawLine(its,y1,its,y2);
    //lc->DrawLine(ijump,y1,ijump,y2);
    lc->DrawLine(i12d,y1,i12d,y2);
    lc->SetLineStyle(kDotted);
    lc->DrawLine(-0.5,1.00,nruns-0.5,1.00);
    c1c->cd();
    lc->SetLineColor(kRed);
    lc->DrawLine(-0.5,1.05,nruns-0.5,1.05);
    lc->DrawLine(-0.5,0.95,nruns-0.5,0.95);


    grun->SetLineWidth(2);
    grun->SetMarkerSize(1.5);
    grun->SetLineColor(kBlack);
    grun->SetMarkerColor(kBlack);
    grun->SetMarkerStyle(kFullCircle);
    grun->Draw("SAMEP");

    bool hb = (tte=="Barrel");
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.040);
    //tex->DrawLatex(  5, 0.82, "2012A");
    //tex->DrawLatex(100, 0.82, "2012B");
    //tex->DrawLatex(220, 0.82, "2012C");
    tex->DrawLatex(ifirst+5, 0.82, "A");
    tex->DrawLatex(i12b+5, 0.82, "B");
    tex->DrawLatex(i12c+5, 0.82, "C");
    tex->DrawLatex(i12d+5, 0.82, "D");
    c1d->cd();
    grun->Draw("SAMEP");
    //tex->DrawLatex(  5, 0.82, "2012A");
    //tex->DrawLatex(100, 0.82, "2012B");
    //tex->DrawLatex(220, 0.82, "2012C");
    tex->DrawLatex(ifirst+5, 0.82, "A");
    tex->DrawLatex(i12b+5, 0.82, "B");
    tex->DrawLatex(i12c+5, 0.82, "C");
    tex->DrawLatex(i12d+5, 0.82, "D");

    //TF1 *f1d = new TF1("f1d2","[0]",ifirst,its);
    TF1 *f1d = new TF1("f1d2","[0]",ifirst,i12b);
    f1d->SetLineStyle(kDashed);
    //f1d->SetLineColor(grun->GetLineColor());
    f1d->SetLineColor(kMagenta+1);
    grun->Fit(f1d,"QRN");
    f1d->DrawClone("SAME");
    //f1d->SetRange(its,i12b);
    f1d->SetRange(i12b,i12c);
    f1d->SetLineColor(kYellow+2);
    grun->Fit(f1d,"QRN");
    f1d->DrawClone("SAME");
    //f1d->SetRange(i12b,ijump);
    f1d->SetRange(i12c,i12d);
    f1d->SetLineColor(kOrange+1);
    grun->Fit(f1d,"QRN");
    f1d->DrawClone("SAME");
    //f1d->SetRange(ijump,i12c);
    f1d->SetRange(i12d,ilast);
    f1d->SetLineColor(kCyan+1);
    grun->Fit(f1d,"QRN");
    f1d->DrawClone("SAME");
    //f1d->SetRange(i12c,ilast);
    //grun->Fit(f1d,"QRN");
    //f1d->DrawClone("SAME");
    c1c->cd();

    TF1 *fera = new TF1(Form("fera_%s",tte.c_str()),"1+[0]/100.",0,500);
    fera->SetLineColor(kRed);
    fera->SetLineWidth(2);

    if (ttc=="PF") {// && tte=="Barrel") {
      //cmsPrel(_lumi);
      if (tte=="Barrel") tex->DrawLatex(350, 1.35,//hb ? 1.175 : 1.35,
					"|#eta_{jet}| < 1.0");
      if (tte=="Transition") tex->DrawLatex(350, 1.35,//hb ? 1.175 : 1.35,
					    "1 < |#eta_{jet}| < 2");
      if (tte=="Endcap") tex->DrawLatex(350, 1.35,//1.175 : 1.35,
					"2 < |#eta_{jet}| < 3");
      if (tte=="Endcap1") tex->DrawLatex(350, 1.35,//1.175 : 1.35,
					"2 < |#eta_{jet}| < 2.5");
      if (tte=="Endcap2") tex->DrawLatex(350, 1.35,//1.175 : 1.35,
					"2.5 < |#eta_{jet}| < 3");

      c1d->cd();
      //cmsPrel(_lumi);
      tex->SetTextSize(0.045);
      if (tte=="Barrel") tex->DrawLatex(350, 1.16,"|#eta_{jet}| < 1.0");
      if (tte=="Transition") tex->DrawLatex(350, 1.16,"1 < |#eta_{jet}| < 2");
      if (tte=="Endcap") tex->DrawLatex(350, 1.16,"2 < |#eta_{jet}| < 3");
      if (tte=="Endcap1") tex->DrawLatex(350, 1.16,"2 < |#eta_{jet}| < 2.5");
      if (tte=="Endcap2") tex->DrawLatex(350, 1.16,"2.5 < |#eta_{jet}| < 3");
      TLegend *legd = new TLegend(0.18,0.80,0.38,0.92,"","brNDC");
      legd->SetFillStyle(kNone);
      legd->SetBorderSize(0);
      legd->SetTextSize(0.040);
      legd->AddEntry(grun, "Trigger average", "LP");
      legd->AddEntry(vmrun[vmrun.size()-2], "HLT_PFJet320", "LP");
      legd->Draw();
      c1c->cd();

      TLegend *leg = new TLegend(0.18,0.60,0.38,0.92,"","brNDC");
      leg->SetFillStyle(kNone);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.040);
      //leg->SetHeader("|#eta_{jet}| < 1.0");
      leg->AddEntry(grun, "Trigger average", "LP");
      for (unsigned int itrig = 0; itrig != vmrun.size(); ++itrig) {
	const char *t = trigs[itrig].c_str();
	leg->AddEntry(vmrun[itrig], Form("%s", t), "LP");

	if (string(t)=="jt370") {

	  tex->SetTextSize(0.030);

	  //fera->SetRange(ifirst, its);
	  fera->SetRange(ifirst, i12c);
	  vmrun[itrig]->Fit(fera, "QRN");
	  fera->DrawClone("SAME");
	  tex->DrawLatex( 40, hb ? 0.92 : 0.84,
			  Form("%+1.1f%%", fera->GetParameter(0)));

	  //fera->SetRange(its, ilast);
	  fera->SetRange(i12c, ilast);
	  vmrun[itrig]->Fit(fera, "QRN");
	  fera->DrawClone("SAME");
	  tex->DrawLatex(125, hb ? 0.92 : 0.84,
			 Form("%+1.1f%%", fera->GetParameter(0)));

	  fera->SetLineWidth(1);
	  fera->SetRange(ifirst, ilast); // full 2012
	  vmrun[itrig]->Fit(fera, "QRN");
	  fera->DrawClone("SAME");
 	  tex->DrawLatex(320, hb ? 1.15 : 1.30,
			 Form("dL=%+1.1f%%", fera->GetParameter(0)));

	  tex->SetTextSize(0.040);
	}
      }
      leg->Draw();
      hc->GetYaxis()->SetTitle("Normalized inclusive jet rate");
    }

    if (jp::pdf) c1c->SaveAs(Form("pdf/RunHistos_%s_%s_lin.pdf",tc,te));
    if (jp::pdf) c1d->SaveAs(Form("pdf/RunHistos_%s_%s_lin2.pdf",tc,te));

    c1->cd();

    cout << " Residual trigger rates normalization" << endl;
    for (unsigned int j = 0; j != trigs.size(); ++j) {
      const string& t = trigs[j];
      cout << Form("%s: %1.3f +/- %1.3f\n", t.c_str(),
		   _norm[t].first, _norm[t].second);
    } // for j
    cout << " New trigger rate normalization" << endl;
    for (unsigned int j = 0; j != trigs.size(); ++j) {
      const string& t = trigs[j];
      cout << Form("scale[\"%s\"][\"%s\"][\"%s\"] =  %1.4g;\n",
		   t.c_str(),te,tc, scale[t][te][tc]*_norm[t].first);
      fout << Form("    scale[\"%s\"][\"%s\"][\"%s\"] =  %1.4g;\n",
		   t.c_str(),te,tc, scale[t][te][tc]*_norm[t].first);
    } // for j
    fout << endl;

  } // for i
  } // for itype


  vector<vector<TGraphErrors*> > graphs(types.size());
  for (unsigned int i = 0; i != graphs.size(); ++i) {
    graphs[i].resize(etas.size());
    for (unsigned int j = 0; j != graphs[i].size(); ++j) {
      graphs[i][j] = new TGraphErrors(0);
    }
  }

  cout << "*******************************" << endl;
  cout << "Trigger rate before and after tech stop" << endl;
  cout << "Trigger" << endl;

  for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
    cout << " & " << etas[ieta];
  }
  cout << " \\\\" << endl;

  for (unsigned int itrig = 0; itrig != trigs.size(); ++itrig) {

    const string& t = trigs[itrig];
    cout << t;

    for (unsigned int itype = 0; itype != types.size(); ++itype) {

      const string& tc = types[itype];
      cout << (tc=="CALO" ? " (Calo)" : " (PF)");

      for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {

	const string& te = etas[ieta];
	pair<double, double> const& af = _rateafter[tc][te][t];
	pair<double, double> const& be = _ratebefore[tc][te][t];
	double delta = af.first - be.first;
	double err = tools::oplus(af.second, be.second);
	cout << Form(" & %1.3f $\\pm$ %1.3f", delta, err);

	TGraphErrors *g = graphs[itype][ieta]; assert(g);
	int n = g->GetN();
	g->SetPoint(n, thr[t]*(1+0.01*ieta-0.05*itype), 100.*delta);
	g->SetPointError(n, 0., 100.*err);
      } // for ieta
      cout << " \\\\" << endl;
    } // for itype
  } // for itrig

  {
    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    c2->SetLogx();
    double ptmin = 40;
    double ptmax = 600;
    TH1D *h = new TH1D("h2",";p_{T} threshold (GeV);"
		       "Rate change after tech stop (%)",
		       int(ptmax-ptmin), ptmin, ptmax);
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();
    h->SetMinimum(-30+0.0001);
    h->SetMaximum(+30-0.0001);
    h->Draw("AXIS");

    TLine *l = new TLine(ptmin, 0., ptmax, 0.);
    l->SetLineStyle(kDashed);
    l->Draw();

    TLegend *leg = new TLegend(0.2,0.92-0.05*etas.size(),0.4,0.92,"","brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->SetTextSize(0.045);
    leg->Draw("SAME");

    TLegend *leg2 = new TLegend(0.5,0.92-0.05*2,0.7,0.92,"","brNDC");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(kNone);
    leg2->SetTextSize(0.045);
    if (types.size()!=1) leg2->Draw("SAME");

    int colors[2] = {kBlack, kBlue};
    int markers[4] = {kFullCircle, kOpenCircle, kOpenSquare, kOpenDiamond};

    // Both overlaid
    for (unsigned int i0 = 0; i0 != graphs.size(); ++i0) {
      for (unsigned int j = 0; j != graphs[i0].size(); ++j) {

	int i = graphs.size()-i0-1; // reverse PF and Calo
	TGraphErrors *g = graphs[i][j]; assert(g);
	g->SetMarkerStyle(markers[j%4]);
	g->SetMarkerSize(2);
	g->SetMarkerColor(colors[i%2]);
	g->SetLineColor(colors[i%2]);
	g->Draw("SAME P");
	if (graphs.size()==2&& graphs[0].size()==1) {
	  if (j==0) leg2->AddEntry(g,types[i]=="" ? "PF":types[i].c_str(),"P");
	}
	else {
	  if (i==0) leg->AddEntry(g,etas[j].c_str(),"P");
	  if (j==0) leg2->AddEntry(g,types[i].c_str(),"P");
	}
      } // for j
    } // for i

    if(jp::pdf) c2->SaveAs("pdf/RunHistos_CaloPF_XDet_RateChange.pdf");

    // Calo and PF separately
    for (unsigned int i = 0; i != graphs.size(); ++i) {

      h->Draw("AXIS");
      l->Draw();
      leg->Clear(); leg->Draw();
      string t = types[i]; if (t=="") t = "PF";
      const char *ct = t.c_str();

      for (unsigned int j = 0; j != graphs[i].size(); ++j) {

	TGraphErrors *g = graphs[i][j]; assert(g);
	g->Draw("SAME P");
	leg->AddEntry(g, (etas[j]+" ("+t+")").c_str(), "P");

	//TF1 *f = new TF1(Form("fg%d%d",i,j),"pol1",70.,260);
	TF1 *f = new TF1(Form("fg%d%d",i,j),"[0]+[1]/x",45.,550.);
	f->SetLineColor(g->GetLineColor());
	if (etas[j]!="") f->SetLineStyle(kDashed);
	g->Fit(f, "QRNF");
	f->Draw("SAME");
      } // for j

      if(jp::pdf) c2->SaveAs(Form("pdf/RunHistos_%s_XDet_RateChange.pdf",ct));
    } // for i

    // Calo - PF difference
    if (types.size()==2) {

      h->GetYaxis()->SetTitle("Calo-PF rate change (%)");
      h->SetMinimum(-4);
      h->SetMaximum(6);
      h->Draw("AXIS");
      l->Draw();
      leg->Clear(); leg->Draw();

      for (unsigned int j = 0; j != graphs[0].size(); ++j) {

	TGraphErrors *g0 = graphs[0][j]; assert(g0);
	TGraphErrors *g1 = graphs[1][j]; assert(g1);
	assert(g0->GetN()==g1->GetN());
	TGraphErrors *g = (TGraphErrors*)g0->Clone();
	for (int i = 0; i != g->GetN(); ++i) {

	  double x = 0.5*(g1->GetX()[i]+g0->GetX()[i]);
	  double dx = 0.5*(g1->GetX()[i]-g0->GetX()[i]);
	  double y = 0.5*(g1->GetY()[i]-g0->GetY()[i]);
	  double dy = tools::oplus(g1->GetEY()[i], g0->GetEY()[i])/sqrt(2.);
	  g->SetPoint(i, x, y);
	  g->SetPointError(i, dx, dy);
	}
	g->Draw("SAME P");
	leg->AddEntry(g, etas[j].c_str(), "P");
      } // for j

      if(jp::pdf) c2->SaveAs("pdf/RunHistos_CaloMinusPF_XDet_RateChange.pdf");
    } // types==2
  }

  cout << "=============================================================" << endl
       << " Bad runs:" << endl
       << "=============================================================" << endl;
  typedef map<int, set<string> >::const_iterator IT;
  for (IT it = mbadruns.begin(); it != mbadruns.end(); ++it) {
    cout << "  _runveto.insert("<<it->first<<"); // ";

    for (set<string>::const_iterator jt = it->second.begin();
	 jt != it->second.end(); ++jt) {
      if (jt != it->second.begin()) cout << ",";
      cout << *jt;
    }
    cout << endl;
  }
  cout << "=============================================================" << endl;


  fout.close();
  if (gROOT->IsBatch()) delete c1;

  curdir->cd();

} // drawRunHistos

void drawRunLumi(string type) {

  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  vector<string> trigs;
  trigs.push_back("mb");
  trigs.push_back("jt6u");
  trigs.push_back("jt15u");

  for (unsigned int i = 0; i != trigs.size(); ++i) {

  } // for i
}

void drawRateVsNvtx(string type="DATA", string algo="PF") {

  assert(algo=="PF" || algo=="Calo");

  TDirectory *curdir = gDirectory;

  setTDRStyle();

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  // Loop not working properly yet, fix
  vector<string> etas;
  //etas.push_back("");
  etas.push_back("Barrel");
  //etas.push_back("Transition");
  //etas.push_back("Endcap");

  vector<string> trigs;
  trigs.push_back("jt40");
  trigs.push_back("jt80");
  trigs.push_back("jt140");
  trigs.push_back("jt200");
  trigs.push_back("jt260");
  trigs.push_back("jt320");
  //trigs.push_back("jt400");
  map<string, string> trigname;
  trigname["jt40"] = "Jet40";
  trigname["jt80"] = "Jet80";
  trigname["jt140"] = "Jet140";
  trigname["jt200"] = "Jet200";
  trigname["jt260"] = "Jet260";
  trigname["jt320"] = "Jet320";
  trigname["jt400"] = "Jet400";

  map<string, double> meanpt;
  // from Eta_0.0-1.3 hselpt (upper limit also)
  meanpt["jt40"] = 89.33;//69.51;//73.8;
  meanpt["jt80"] = 157.7;//136.7;//146.6;
  meanpt["jt140"] = 248.4;//215.4;//245.5;
  meanpt["jt200"] = 334.8;//275.6;//303.2;
  meanpt["jt260"] = 436.7;//366.9;//366.9;
  meanpt["jt320"] = 605.8;//461.8;//474.7;
  meanpt["jt400"] = 1104;//602.6;//556.1;

  map<string, map<string, map<string, double> > > scale;
  if (string(jp::algo)=="AK4") {

    // 9/fb
    // V4
    scale["jt40"]["Barrel"]["PF"] =  5.207e+06;
    scale["jt80"]["Barrel"]["PF"] =  2.519e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.219e+04;
    scale["jt200"]["Barrel"]["PF"] =  2633;
    scale["jt260"]["Barrel"]["PF"] =  986.1;
    scale["jt320"]["Barrel"]["PF"] =  386.9;
    scale["jt400"]["Barrel"]["PF"] =  151.8;
    /* //V1
    scale["jt40"]["Barrel"]["PF"] =  4.956e+06;
    scale["jt80"]["Barrel"]["PF"] =  2.451e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.198e+04;
    scale["jt200"]["Barrel"]["PF"] =  2594;
    scale["jt260"]["Barrel"]["PF"] =  972.5;
    scale["jt320"]["Barrel"]["PF"] =  381.8;
    scale["jt400"]["Barrel"]["PF"] =  149.9;
    */

    /*3fb*/
    /*
    scale["jt40"]["Barrel"]["PF"] =  5.268e+06;
    scale["jt80"]["Barrel"]["PF"] =  2.51e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.207e+04;
    scale["jt200"]["Barrel"]["PF"] =  2608;
    scale["jt260"]["Barrel"]["PF"] =  974.6;
    scale["jt320"]["Barrel"]["PF"] =  382.5;
    scale["jt400"]["Barrel"]["PF"] =  149.5;
    */
    scale["jt40"]["Transition"]["PF"] =  4.043e+06;
    scale["jt80"]["Transition"]["PF"] =  1.761e+05;
    scale["jt140"]["Transition"]["PF"] =  7341;
    scale["jt200"]["Transition"]["PF"] =  1411;
    scale["jt260"]["Transition"]["PF"] =  478.3;
    scale["jt320"]["Transition"]["PF"] =  164.6;
    scale["jt400"]["Transition"]["PF"] =  54.52;

    scale["jt40"]["Endcap"]["PF"] =  2.005e+06;
    scale["jt80"]["Endcap"]["PF"] =  6.902e+04;
    scale["jt140"]["Endcap"]["PF"] =  1618;
    scale["jt200"]["Endcap"]["PF"] =  175.8;
    scale["jt260"]["Endcap"]["PF"] =  34.86;
    scale["jt320"]["Endcap"]["PF"] =  6.757;
    scale["jt400"]["Endcap"]["PF"] =  0.7181;


    /*1.6fb*/
    /*
    scale["jt40"]["Barrel"]["PF"] =  5.067e+06;
    scale["jt80"]["Barrel"]["PF"] =  2.426e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.18e+04;
    scale["jt200"]["Barrel"]["PF"] =  2527;
    scale["jt260"]["Barrel"]["PF"] =  953.8;
    scale["jt320"]["Barrel"]["PF"] =  376.3;
    scale["jt400"]["Barrel"]["PF"] =  147.1;

    scale["jt40"]["Transition"]["PF"] =  3.911e+06;
    scale["jt80"]["Transition"]["PF"] =  1.709e+05;
    scale["jt140"]["Transition"]["PF"] =  7182;
    scale["jt200"]["Transition"]["PF"] =  1368;
    scale["jt260"]["Transition"]["PF"] =  467.6;
    scale["jt320"]["Transition"]["PF"] =  161;
    scale["jt400"]["Transition"]["PF"] =  53.37;

    scale["jt40"]["Endcap"]["PF"] =  1.985e+06;
    scale["jt80"]["Endcap"]["PF"] =  6.938e+04;
    scale["jt140"]["Endcap"]["PF"] =  1627;
    scale["jt200"]["Endcap"]["PF"] =  175.8;
    scale["jt260"]["Endcap"]["PF"] =  35.06;
    scale["jt320"]["Endcap"]["PF"] =  6.874;
    scale["jt400"]["Endcap"]["PF"] =  0.7181;
    */
  }
  if (string(jp::algo)=="AK7") {
    // 20/fb, new thresholds
    scale["jt40"]["Barrel"]["PF"] =  2.424e+06;
    scale["jt80"]["Barrel"]["PF"] =  1.681e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.535e+04;
    scale["jt200"]["Barrel"]["PF"] =  3229;
    scale["jt260"]["Barrel"]["PF"] =  744.9;
    scale["jt320"]["Barrel"]["PF"] =  181.4;
    scale["jt400"]["Barrel"]["PF"] =  2.492;
    // 20/fb,old thresholds
    /*
    scale["jt40"]["Barrel"]["PF"] =  8.269e+06;
    scale["jt80"]["Barrel"]["PF"] =  3.391e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.535e+04;
    scale["jt200"]["Barrel"]["PF"] =  3229;
    scale["jt260"]["Barrel"]["PF"] =  1193;
    scale["jt320"]["Barrel"]["PF"] =  462.4;
    scale["jt400"]["Barrel"]["PF"] =  179.3;
    */
    /*
    scale["jt40"]["Barrel"]["PF"] =  8.698e+06;
    scale["jt80"]["Barrel"]["PF"] =  3.361e+05;
    scale["jt140"]["Barrel"]["PF"] =  1.504e+04;
    scale["jt200"]["Barrel"]["PF"] =  3154;
    scale["jt260"]["Barrel"]["PF"] =  1163;
    scale["jt320"]["Barrel"]["PF"] =  450.9;
    scale["jt400"]["Barrel"]["PF"] =  174.3;
    */

    scale["jt40"]["Transition"]["PF"] =  6.453e+06;
    scale["jt80"]["Transition"]["PF"] =  2.323e+05;
    scale["jt140"]["Transition"]["PF"] =  9107;
    scale["jt200"]["Transition"]["PF"] =  1706;
    scale["jt260"]["Transition"]["PF"] =  569.7;
    scale["jt320"]["Transition"]["PF"] =  194.8;
    scale["jt400"]["Transition"]["PF"] =  64.02;

    scale["jt40"]["Endcap"]["PF"] =  2.884e+06;
    scale["jt80"]["Endcap"]["PF"] =  8.321e+04;
    scale["jt140"]["Endcap"]["PF"] =  1790;
    scale["jt200"]["Endcap"]["PF"] =  183;
    scale["jt260"]["Endcap"]["PF"] =  37.78;
    scale["jt320"]["Endcap"]["PF"] =  6.933;
    scale["jt400"]["Endcap"]["PF"] =  0.7181;

    scale["jt40"]["Endcap1"]["PF"] =  1.753e+06;
    scale["jt80"]["Endcap1"]["PF"] =  5.691e+04;
    scale["jt140"]["Endcap1"]["PF"] =  1487;
    scale["jt200"]["Endcap1"]["PF"] =  168.3;
    scale["jt260"]["Endcap1"]["PF"] =  36.46;
    scale["jt320"]["Endcap1"]["PF"] =  6.858;
    scale["jt400"]["Endcap1"]["PF"] =  0.7181;

    scale["jt40"]["Endcap2"]["PF"] =  1.117e+06;
    scale["jt80"]["Endcap2"]["PF"] =  2.569e+04;
    scale["jt140"]["Endcap2"]["PF"] =  259.3;
    scale["jt200"]["Endcap2"]["PF"] =  12.96;
    scale["jt260"]["Endcap2"]["PF"] =  0.888;
    scale["jt320"]["Endcap2"]["PF"] =  0.06549;
    scale["jt400"]["Endcap2"]["PF"] =  0.01417;

  }

  set<int> _runveto;
  //_runveto.insert(132959); //mb

  map<string, set<int> > _trigveto;
  //_trigveto["jt6u"].insert(141958); // jt6u

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TCanvas *c1ts = new TCanvas("c1ts","c1ts",600,600);

  const char *a = algo.c_str();
  TH1D *h = new TH1D("h",Form(";#LTN_{PV,good}#GT;Normalized rate per run (%s%s)",
			      a, jp::algo),
		     //100,0.8,15.);
		     100,0.,30.);
  h->SetMinimum(0.7);
  h->SetMaximum(1.3);

  double dldnpv(0), dodnpv(0);
  /*
  TGraphErrors *gro = new TGraphErrors(trigs.size());
  TGraphErrors *gro2 = new TGraphErrors(trigs.size());
  TGraphErrors *grl = new TGraphErrors(trigs.size());
  TGraphErrors *grl2 = new TGraphErrors(trigs.size());
  */
  for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {

    TGraphErrors *grl = new TGraphErrors(0);
    TGraphErrors *grl2 = new TGraphErrors(0);
    TGraphErrors *gro = new TGraphErrors(0);
    TGraphErrors *gro2 = new TGraphErrors(0);
    TGraphErrors *grs = new TGraphErrors(0);
    //TGraphErrors *grs2 = new TGraphErrors(0);

  const char *te = etas[ieta].c_str();
  double eta = 0;
  if (etas[ieta]=="Transition") eta = 1.0;
  if (etas[ieta]=="Endcap") eta = 2.0;

  for (unsigned int i = 0; i != trigs.size(); ++i) {

    const char *t = trigs[i].c_str();

    TH1D *hr = (TH1D*)fin->Get(Form("Runs%s/r_%s",te,t)); assert(hr);
    TH1D *hn = (TH1D*)fin->Get(Form("Runs%s/npvgood_%s",te,t)); assert(hn);
    TH1D *runs = (TH1D*)fin->Get(Form("Runs%s/runs",te)); assert(hn);

    if (algo=="Calo") {
      hr = (TH1D*)fin->Get(Form("Runs%s/r_c%s",te,t)); assert(hr);
    }

    TGraphErrors *g = new TGraphErrors(0);
    TGraphErrors *gu = new TGraphErrors(0);
    TGraphErrors *gx = new TGraphErrors(0);
    TGraphErrors *gy = new TGraphErrors(0);
    TGraphErrors *gb = new TGraphErrors(0); // before tech stop
    TGraphErrors *ga = new TGraphErrors(0); // after tech stop
    TGraphErrors *gby = new TGraphErrors(0); // before tech stop
    TGraphErrors *gay = new TGraphErrors(0); // after tech stop
    for (int j = 1; j != hr->GetNbinsX()+1; ++j) {

      double x = hn->GetBinContent(j);
      double ex = hn->GetBinError(j);
      double y = hr->GetBinContent(j);
      double ey = hr->GetBinError(j);
      int run = runs->GetBinContent(j);
      if (_runveto.find(run)==_runveto.end() &&
	  _trigveto[trigs[i]].find(run)==_trigveto[trigs[i]].end()
	  ) {

	double maxe = 0.05;
	double maxdy = 0.3;
	double sc = scale[t][te][a]; //assert(sc!=0);
	if (sc==0) {
	  cout << "t="<<t<<", te="<<te<<", a="<<a
	       << ",_algo="<<jp::algo<<endl<<flush;
	  assert(sc!=0);
	}
	if (x!=0 && ex/x<maxe && ex!=0 && y!=0 && ey/y<maxe && ey!=0 &&
	    fabs(y/sc-1) < maxdy) {

	  int nu = gu->GetN();
	  gu->SetPoint(nu, x, y/sc);
	  gu->SetPointError(nu, ex, ey/sc);

	  int n = g->GetN();
	  double pt = meanpt[t];
	  _ak7 = (string(jp::algo)=="AK7"); // change PU smearing
	  //double c = unfold(pt, eta, x) / unfold(pt, eta, 12.);
	  double c = unfold(pt, eta, x) / unfold(pt, eta, 13.5);
	  g->SetPoint(n, x, y/sc*c);
	  g->SetPointError(n, ex, ey/sc*c);

	  int n0 = gx->GetN();
	  gx->SetPoint(n0, x, y/sc*c);
	  gx->SetPointError(n0, ex, ey/sc*c);
	  int m0 = gy->GetN();
	  gy->SetPoint(m0, y/sc*c, x);
	  gy->SetPointError(m0, ey/sc*c, ex);

	  if (run>192000) {
	    int na = ga->GetN();
	    ga->SetPoint(na, x, y/sc*c);
	    ga->SetPointError(na, ex, ey/sc*c);
	    gay->SetPoint(na, y/sc*c, x);
	    gay->SetPointError(na, ey/sc*c, ex);
	  }
	  else {
	    int nb = gb->GetN();
	    gb->SetPoint(nb, x, y/sc*c);
	    gb->SetPointError(nb, ex, ey/sc*c);
	    gby->SetPoint(nb, y/sc*c, x);
	    gby->SetPointError(nb, ey/sc*c, ex);
	  }
	}
      }
    } // for j

    c1->cd();
    h->Draw("AXIS");

    TGraphErrors *gsys = new TGraphErrors(100);
    for (int j = 0; j != 100; ++j) {
      double npv = 1.0 + 9.*j/(100.-1);
      gsys->SetPoint(j, npv, 1.);
    }
    gsys->SetFillColor(kYellow+1);
    gsys->SetFillStyle(1001);
    gsys->Draw("SAME E3");

    g->Draw("SAMEP");

    gu->SetMarkerStyle(kOpenCircle);
    gu->Draw("SAMEP");

    TF1 *fit = new TF1(Form("fit_%s",t),"[0]+[1]*x",0,30);
    fit->SetParameters(1., 0.);
    fit->SetLineColor(kBlue);
    fit->SetLineWidth(2);
    gx->Fit(fit, "QRNF"); // F needed to get emat below
    fit->Draw("SAME");

    TMatrixD em(2,2);
    gMinuit->mnemat(&em[0][0],2);
    TF1 *efit0 = new TF1(Form("efit_%s",t),
			 "[0]+[1]*x+[2]*sqrt([3]+2*[4]*x+[5]*x*x)",0,30);
    efit0->SetParameters(fit->GetParameter(0),fit->GetParameter(1),
			 max(1.,sqrt(fit->GetChisquare()/fit->GetNDF())),
			 em[0][0], em[0][1], em[1][1]);
    efit0->SetLineColor(kBlue+2);
    efit0->SetLineStyle(kDotted);
    efit0->DrawClone("SAME");
    efit0->SetParameter(2, -efit0->GetParameter(2));
    efit0->DrawClone("SAME");

    TF1 *fit2 = new TF1(Form("fit2_%s",t),"pol0",0,30);
    gy->Fit(fit2, "QRN");
    double x = fit2->GetParameter(0);
    double y = fit->Eval(x);

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(h->GetXaxis()->GetXmin(), 1, h->GetXaxis()->GetXmax(), 1);
    TEllipse *circ = new TEllipse(x,y,0.15,0.01);
    circ->SetFillStyle(1001);
    circ->SetFillColor(kRed);
    circ->Draw("SAME");

    //cmsPrel(_lumi);
    TLatex *tex = new TLatex(0.75,0.85,trigname[trigs[i]].c_str());
    tex->SetTextSize(0.045);
    tex->SetNDC();
    tex->Draw();
    if (etas[ieta]=="Barrel")     tex->DrawLatex(0.75, 0.80, "|#eta|<1.0");
    if (etas[ieta]=="Transition") tex->DrawLatex(0.75, 0.80, "1<|#eta|<2");
    if (etas[ieta]=="Endcap")     tex->DrawLatex(0.75, 0.80, "2<|#eta|<3");

    double mpt = meanpt[trigs[i]]; assert(mpt!=0);
    double dsigmadjec = 5.; // change in xsec per change in JEC
    tex->DrawLatex(0.20, 0.25, Form("#DeltaO/#LTN_{PV}#GT #approx %1.2f#pm%1.2f GeV"
				    ", or",
				    fit->GetParameter(1)*mpt/dsigmadjec,
				    fit->GetParError(1)*mpt/dsigmadjec));
    tex->DrawLatex(0.20, 0.20, Form("#DeltaL/#LTN_{PV}#GT #approx %1.2f#pm%1.2f%%",
				    fit->GetParameter(1)*100.,
				    fit->GetParError(1)*100.));
    if (mpt*cosh(eta) < 2000.) {

      int n = grl->GetN();
      grl->SetPoint(n, mpt*0.99, 100.*fit->GetParameter(1));
      grl->SetPointError(n, 0., 100.*fit->GetParError(1));
      dodnpv = 0.17;
      grl2->SetPoint(n, mpt*1.01,100.*(fit->GetParameter(1)
				       - dodnpv/mpt*dsigmadjec));
      grl2->SetPointError(n, 0., 100.*fit->GetParError(1));
      //
      gro->SetPoint(n, mpt*0.99, fit->GetParameter(1)*mpt/dsigmadjec);
      gro->SetPointError(n, 0., fit->GetParError(1)*mpt/dsigmadjec);
      dldnpv = -0.12*0.01;
      dldnpv += pow(5.*1./mpt, 2.);
      gro2->SetPoint(n, mpt*1.01, (fit->GetParameter(1)-dldnpv)*mpt/dsigmadjec);
      gro2->SetPointError(n, 0., fit->GetParError(1)*mpt/dsigmadjec);
      //
      grs->SetPoint(n, mpt*0.99, fit->GetParameter(1)>0 ? sqrt(fit->GetParameter(1)*pow(mpt/dsigmadjec,2)) : 0);
      grs->SetPointError(n, 0., 0.5*sqrt(1./fabs(fit->GetParameter(1))*pow(mpt/dsigmadjec,2))*fit->GetParError(1));

    }

    gPad->RedrawAxis();

    if(jp::pdf) c1->SaveAs(Form("pdf/RunHistos_%s_%s_RateVsNvtx_%s.pdf",a,te,t));

    // Draw again with before/after separation
    c1ts->cd();
    h->Draw("AXIS");

    //TLine *l = new TLine(h->GetXaxis()->GetXmin(), 1.,
    //		 h->GetXaxis()->GetXmax(), 1.);
    l->DrawLine(h->GetXaxis()->GetXmin(), 1, h->GetXaxis()->GetXmax(), 1);
    l->SetLineStyle(kDashed);
    l->Draw();

    gb->SetMarkerColor(kBlue);
    gb->SetLineColor(kBlue);
    gb->Draw("SAMEP");

    ga->SetMarkerColor(kRed);
    ga->SetLineColor(kRed);
    ga->Draw("SAMEP");

    TF1 *fitb = new TF1(Form("fitb_%s",t),"pol1",0,30);
    fitb->SetParameters(1.,0.); // Needed with "F" below
    fitb->SetLineColor(kBlue+2);
    fitb->SetLineWidth(2);
    gb->Fit(fitb, "QRNF"); // "F" needed to get emat below
    fitb->Draw("SAME");

    TMatrixD emat(2,2);
    gMinuit->mnemat(&emat[0][0],2);

    TF1 *efit = new TF1(Form("efit_%s",t),
			"[0]+[1]*x+[2]*sqrt([3]+2*[4]*x+[5]*x*x)",0,30);
    efit->SetParameters(fitb->GetParameter(0),fitb->GetParameter(1),
			max(1.,sqrt(fitb->GetChisquare()/fitb->GetNDF())),
			emat[0][0], emat[0][1], emat[1][1]);
    efit->SetLineColor(kBlue+2);
    efit->SetLineStyle(kDotted);//kDashed);
    efit->DrawClone("SAME");
    efit->SetParameter(2, -efit->GetParameter(2));
    efit->DrawClone("SAME");

    TF1 *fitb2 = new TF1(Form("fitb2_%s",t),"pol0",0,30);
    gby->Fit(fitb2, "QRN");
    double bx = fitb2->GetParameter(0);
    double by = fitb->Eval(bx);
    l->DrawLine(h->GetXaxis()->GetXmin(), 1, h->GetXaxis()->GetXmax(), 1);
    TEllipse *circb = new TEllipse(bx,by,0.15,0.01);
    circb->SetFillStyle(1001);
    circb->SetFillColor(kBlue+2);
    circb->Draw("SAME");

    TF1 *fita = new TF1(Form("fita_%s",t),"pol1",0,30);
    fita->SetLineColor(kRed+2);
    fita->SetLineWidth(2);
    ga->Fit(fita, "QRN");
    fita->Draw("SAME");

    gMinuit->mnemat(&emat[0][0],2);

    efit->SetParameters(fita->GetParameter(0),fita->GetParameter(1),
			max(1.,sqrt(fita->GetChisquare()/fita->GetNDF())),
			emat[0][0], emat[0][1], emat[1][1]);
    efit->SetLineColor(kRed+2);
    efit->DrawClone("SAME");
    efit->SetParameter(2, -efit->GetParameter(2));
    efit->DrawClone("SAME");

    TF1 *fita2 = new TF1(Form("fita2_%s",t),"pol0",0,30);
    gay->Fit(fita2, "QRN");
    double ax = fita2->GetParameter(0);
    double ay = fita->Eval(ax);
    TEllipse *circa = new TEllipse(ax,ay,0.15,0.01);
    circa->SetFillStyle(1001);
    circa->SetFillColor(kRed+2);
    circa->Draw("SAME");

    fit->SetLineColor(kGreen+2);
    fit->Draw("SAME");
    circ->SetFillColor(kGreen+2);
    circ->Draw("SAME");

    tex->Draw();

    //if(jp::pdf) c1ts->SaveAs(Form("pdf/RunHistos_%s_RateVsNvtxTS_%s.pdf",a,t));
  } // for i


  TCanvas *c2a = new TCanvas("c2a","c2a",600,600);
  gPad->SetLogx();

  TH1D *h2a = new TH1D("h2a",";Mean trigger p_{T} (GeV);#DeltaL / #LTN_{PV}#GT (%)",
		       int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h2a->GetXaxis()->SetNoExponent();
  h2a->GetXaxis()->SetMoreLogLabels();
  h2a->SetMinimum(-2);//-1.5);
  h2a->SetMaximum(2);//1.5);

  h2a->Draw("HIST");//AXIS");

  grl->SetMarkerStyle(kFullCircle);
  grl->Draw("SAMEP");

  TF1 *f1a = new TF1("f1a","[0]",jp::xmin,min(2000./cosh(eta), jp::xmax));
  grl->Fit(f1a,"QN");
  f1a->SetLineColor(kBlack);
  f1a->DrawClone("SAME");

  grl2->SetLineColor(kBlue);
  grl2->SetMarkerColor(kBlue);
  grl2->SetMarkerStyle(kOpenCircle);
  //grl2->Draw("SAMEP");

  //grl2->Fit(f1a,"QN");
  f1a->SetLineColor(kBlue);
  //f1a->DrawClone("SAME");

  TLegend *leg2a = new TLegend(0.20,0.77,0.40,0.92,"","brNDC");
  leg2a->SetFillStyle(kNone);
  leg2a->SetBorderSize(0);
  leg2a->SetTextSize(0.045);
  if (etas[ieta]=="Barrel")     leg2a->SetHeader("Inclusive jets, |#eta|<1.0");
  if (etas[ieta]=="Transition") leg2a->SetHeader("Inclusive jets, 1<|#eta|<2");
  if (etas[ieta]=="Endcap")     leg2a->SetHeader("Inclusive jets, 2<|#eta|<3");
  leg2a->AddEntry(grl,"Data","PL");
  //leg2a->AddEntry(grl2,Form("shifted by dO/dN_{^{PV}} = %1.2f GeV",
  //		   dodnpv),"PL");
  leg2a->Draw();

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->SetLineColor(kBlue);
  tex->DrawLatex(0.20,0.20,Form("dL/d#LTN_{PV}#GT = %1.2f#pm%1.2f %% (#chi^{2}/NDF=%1.1f)",
				f1a->GetParameter(0),
				f1a->GetParError(0),
				f1a->GetChisquare()/f1a->GetNDF()));

  //cmsPrel(_lumi);

  if (jp::pdf) c2a->SaveAs(Form("pdf/RunHistos_%s_%s_DeltaLvsPt.pdf",a,te));


  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  gPad->SetLogx();

  TH1D *h2 = new TH1D("h2",";Mean trigger p_{T} (GeV);#DeltaO / #LTN_{PV}#GT (GeV)",
		      int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h2->GetXaxis()->SetNoExponent();
  h2->GetXaxis()->SetMoreLogLabels();
  h2->SetMinimum(-1);//-0.5);//-1.5);
  h2->SetMaximum(1);//0.5);//0.8);

  h2->Draw("HIST");//AXIS");

  gro->SetMarkerStyle(kFullCircle);
  gro->Draw("SAMEP");

  TF1 *f1 = new TF1("f1","[0]+[1]*x",jp::xmin,min(2000./cosh(eta), jp::xmax));
  f1->FixParameter(1,0);
  gro->Fit(f1,"QN");
  f1->SetLineColor(kBlack);
  f1->DrawClone("SAME");

  cout << "dL/dNpv = ";
  cout << f1->GetParameter(0) << "+" << f1->GetParameter(1) << "x" << endl;

  TF1 *f2 = new TF1("f2","[0]+[1]*log(x)",jp::xmin,jp::xmax);
  gro->Fit(f2,"QN");
  f2->SetLineStyle(kDashed);
  f2->SetLineColor(kBlack);
  //f2->DrawClone("SAME");

  gro2->SetLineColor(kBlue);
  gro2->SetMarkerColor(kBlue);
  gro2->SetMarkerStyle(kOpenCircle);
  //gro2->Draw("SAMEP");

  //gro2->Fit(f1,"QN");
  f1->SetLineColor(kBlue);
  //f1->DrawClone("SAME");

  gro2->Fit(f2,"QN");
  f2->SetLineColor(kBlue);
  //f2->DrawClone("SAME");

  TLegend *leg2 = new TLegend(0.20,0.77,0.40,0.92,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  if (etas[ieta]=="Barrel")     leg2->SetHeader("Inclusive jets, |#eta|<1.0");
  if (etas[ieta]=="Transition") leg2->SetHeader("Inclusive jets, 1<|#eta|<2");
  if (etas[ieta]=="Encap")      leg2->SetHeader("Inclusive jets, 2<|#eta|<3");
  leg2->AddEntry(gro,"Data","PL");
  //leg2->AddEntry(gro2,Form("shifted by dL/dN_{^{PV}} = %1.2f%%",
  //		   dldnpv*100.),"PL");
  leg2->Draw();

  //tex->DrawLatex(0.20,0.20,Form("dO/dN_{PV} = %1.2f#pm%1.2f GeV",
  //			f1->GetParameter(0),
  //			f1->GetParError(0)));
  tex->DrawLatex(0.20,0.20,Form("dO/d#LTN_{PV}#GT = %1.2f#pm%1.2f GeV (#chi^{2}/NDF=%1.1f)",
				f1->GetParameter(0),
				f1->GetParError(0),
				f1->GetChisquare()/f1->GetNDF()));

  //cmsPrel(_lumi);

  if (jp::pdf) c2->SaveAs(Form("pdf/RunHistos_%s_%s_DeltaOvsPt.pdf",a,te));


  TCanvas *c2c = new TCanvas("c2c","c2c",600,600);
  gPad->SetLogx();

  TH1D *h2c = new TH1D("h2c",";Mean trigger p_{T} (GeV);#DeltaJER / #LTN_{PV}#GT (GeV)",
		       int(jp::xmax-jp::xmin),jp::xmin,jp::xmax);
  h2c->GetXaxis()->SetNoExponent();
  h2c->GetXaxis()->SetMoreLogLabels();
  h2c->SetMinimum(-10);
  h2c->SetMaximum(10);

  h2c->Draw("HIST");

  grs->SetMarkerStyle(kFullCircle);
  grs->Draw("SAMEP");

  TF1 *f1c = new TF1("f1c","[0]",jp::xmin,min(2000./cosh(eta), jp::xmax));
  grs->Fit(f1c,"QN");
  f1c->SetLineColor(kBlack);
  f1c->DrawClone("SAME");

  TLegend *leg2c = new TLegend(0.20,0.77,0.40,0.92,"","brNDC");
  leg2c->SetFillStyle(kNone);
  leg2c->SetBorderSize(0);
  leg2c->SetTextSize(0.045);
  if (etas[ieta]=="Barrel")     leg2c->SetHeader("Inclusive jets, |#eta|<1.0");
  if (etas[ieta]=="Transition") leg2c->SetHeader("Inclusive jets, 1<|#eta|<2");
  if (etas[ieta]=="Endcap")     leg2c->SetHeader("Inclusive jets, 2<|#eta|<3");
  leg2c->AddEntry(grs,"Data","PL");
  leg2c->Draw();

  tex->DrawLatex(0.20,0.20,Form("dJER/d#LTN_{PV}#GT = %1.2f#pm%1.1f GeV (#chi^{2}/NDF=%1.1f)",
				f1c->GetParameter(0),
				f1c->GetParError(0),
				f1c->GetChisquare()/f1c->GetNDF()));

  //cmsPrel(_lumi);

  if (jp::pdf) c2c->SaveAs(Form("pdf/RunHistos_%s_%s_DeltaJERvsPt.pdf",a,te));
  } // for ieta
} // drawRateVsNvtx



// PFjet composition plots vs run
void drawRunComposition(string type = "DATA") {

  TDirectory *curdir = gDirectory;

  setTDRStyle();

  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  TFile *fmc = new TFile("output-MC-1.root","READ");
  //TFile *fmc = new TFile("backup/sep26/output-MC-1.root","READ");
  //cout << "*** Using backup MC!!! **** " << endl;
  assert(fmc && !fmc->IsZombie());

  ofstream fout("drawRunHistos12b.txt",ios::out);

  const bool _tp = true;//false;
  const bool _dm = !_tp;
  const char *sm = (_tp ? "tp" : "");

  vector<string> etas;
  //etas.push_back("");
  etas.push_back("Barrel");
  etas.push_back("Transition");
  etas.push_back("Endcap1");
  etas.push_back("Endcap2");
  etas.push_back("Endcap"); // must have last

  vector<string> trgs;
  //trgs.push_back("jt40");
  //trgs.push_back("jt80");
  //trgs.push_back("jt140");
  //trgs.push_back("jt200");
  //trgs.push_back("jt260");
  trgs.push_back("jt320");
  //trgs.push_back("jt400");

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  TCanvas *c4 = new TCanvas("c4","c4",600,600);
  TCanvas *c5 = new TCanvas("c5","c5",600,600);

  double refchf(0), refnef(0), refnhf(0);

  for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {

    const char *se = etas[ieta].c_str();

    fout << "// " << se << " jt140" << endl;

    // Mean for trigger jt140
    if (etas[ieta]=="Barrel") {
      if (_tp) {
        refchf = 0.6317;
        refnef = 0.2919;
        refnhf = 0.0678;
      }
      if (_dm) {
	refchf = 0.6274;
	refnef = 0.2867;
	refnhf = 0.0766;
      }
    }
    if (etas[ieta]=="Transition") {
      if (_tp) {
        refchf = 0.6408;
        refnef = 0.2659;
        refnhf = 0.0794;
      }
      if (_dm) {
	refchf = 0.6351;
	refnef = 0.2650;
	refnhf = 0.0851;
      }
    }
    if (etas[ieta]=="Endcap1") {
      if (_tp) {
        refchf = 0.5383;
        refnef = 0.3601;
        refnhf = 0.0981;
      }
      if (_dm) {
	refchf = 0.4361; // tp
	refnef = 0.3630; // tp
	refnhf = 0.1967; // tp
      }
    }
    if (etas[ieta]=="Endcap2") {
      if (_tp) {
        refchf = 0.1216;
        refnef = 0.4224;
        refnhf = 0.4427;
      }
      if (_dm) {
	refchf = 0.4361; // tp
	refnef = 0.3630; // tp
	refnhf = 0.1967; // tp
      }
    }
    if (etas[ieta]=="Endcap") {
      if (_tp) {
        refchf = 0.4305;
        refnef = 0.3762;
        refnhf = 0.1872;
      }
      if (_dm) {
	refchf = 0.4361; // tp
	refnef = 0.3630; // tp
	refnhf = 0.1967; // tp
      }
    }

    TH1D *runs = (TH1D*)fin->Get(Form("Runs%s/runs",se));
    assert(runs);
    int nruns = runs->GetNbinsX();
    float lumi;
    assert(sscanf(runs->GetTitle(),"runs %f pb-1",&lumi)==1);

    TH1D *h = new TH1D("h","",nruns,-0.5,nruns-0.5);
    h->SetMinimum(0.0);
    h->SetMaximum(0.9);
    h->GetYaxis()->SetTitle("PF composition per run");
    h->GetXaxis()->SetTitle("run");
    vector<int> vrun(nruns);
    for (int i = 1; i != nruns+1; ++i) {
      int run = runs->GetBinContent(i);
      if (run>0) {
	vrun[i-1] = run;
	if (i>1) assert(vrun[i-2]<vrun[i-1] || vrun[i-1]==0);
      }
      else {
	vrun[i-1] = vrun[max(0,i-2)];
      }
    } // for i


    //int i190456 = TMath::BinarySearch(nruns, &vrun[0], 190456);
    int i190645 = TMath::BinarySearch(nruns, &vrun[0], 190645);
    int i191859 = TMath::BinarySearch(nruns, &vrun[0], 191859); // before TS
    int i193093 = TMath::BinarySearch(nruns, &vrun[0], 193093); // after TS
    int i193621 = TMath::BinarySearch(nruns, &vrun[0], 193621); // last JetMon
    int i193834 = TMath::BinarySearch(nruns, &vrun[0], 193834); // first JetHT
    int i195530 = TMath::BinarySearch(nruns, &vrun[0], 195530);
    int i195540 = TMath::BinarySearch(nruns, &vrun[0], 195540);
    int i196531 = TMath::BinarySearch(nruns, &vrun[0], 196531);

    double ifirst = i190645;
    double its = 0.5*(i191859+i193093);
    double iht = 0.5*(i193621+i193834);
    double ijump = 0.5*(i195530+i195540);
    //double ilast = i196531;

    int i190456 = TMath::BinarySearch(nruns, &vrun[0], 190456); // 2012A first
    //int i193621 = TMath::BinarySearch(nruns, &vrun[0], 193621); // 2012A last
    //int i193834 = TMath::BinarySearch(nruns, &vrun[0], 193834); // 2012B first
    //int i196531 = TMath::BinarySearch(nruns, &vrun[0], 196531); // 2012B last
    int i198049 = TMath::BinarySearch(nruns, &vrun[0], 198049); // 2012C first
    int i203742 = TMath::BinarySearch(nruns, &vrun[0], 203742); // 2012C last
    int i203777 = TMath::BinarySearch(nruns, &vrun[0], 203777); // 2012D first
    int i208686 = TMath::BinarySearch(nruns, &vrun[0], 208686); // 2012D last

    //double ifirst = i190456;
    double i12b = 0.5*(i193621 + i193834);
    double i12c = 0.5*(i196531 + i198049);
    double i12d = 0.5*(i203742 + i203777);
    double ilast = i208686;

    cout << "ifirst = " << ifirst << endl;
    cout << "its = " << its << endl;
    cout << "ilast = " << ilast << endl;

    const int runbins[] =
    // 0.57/fb division
      //{0, 5, 10, 15, 20, 25, /*TS*/ 30, 35, 41};
    // 1.6/fb divisions
    //{0, 6, 12, 18, 23, 28, /*TS*/ 33, 37, 42, /*HT*/ 46,
    // 51, 56, 61, 66, 71, 76, 81,
    // /*3fb*/ 86, 91, 96, 101, 106, 111, 116, 122,
    // /*5fb*/127, 132, 137, 142, 147, 152, 157, 162, 167, 172, 178, 184};
    //
    // 20/fb divisions, 572 runs, epoch firsts: 1, 63, 209, 389, (569=0)
    {/*A*/ 0, 31, /*B*/ 63, 111, 160, /*C*/ 209, 269, 329,
     /*D*/ 389, 449, 509, 569};//572};

    const int nrunbins = sizeof(runbins)/sizeof(runbins[0])-1;

    TLegend *leg = new TLegend(0.20,0.78,0.40,0.92,"","brNDC");
    leg->SetFillStyle(kNone); leg->SetBorderSize(0); leg->SetTextSize(0.045);

    TF1 *f0 = new TF1("f0","[0]",ifirst,ilast);
    TF1 *f1 = new TF1("f1","[0]+[1]*x",ifirst,ilast);

    TMultiGraph *mchf = new TMultiGraph();
    TMultiGraph *mnef = new TMultiGraph();
    TMultiGraph *mnhf = new TMultiGraph();
    double mcchf(0);
    double mcnef(0);
    double mcnhf(0);

    double y1(0), y2(0), y3(0);
    if (etas[ieta]=="Barrel") {
      y1 = 0.0; y2 = 0.5; y3 = 1.0;
    }
    if (etas[ieta]=="Transition") {
      y1 = 1.0; y2 = 1.5; y3 = 2.0;
    }
    if (etas[ieta]=="Endcap") {
      y1 = 2.0; y2 = 2.5; y3 = 3.0;
    }
    if (etas[ieta]=="Endcap1") {
      y1 = 2.0; y2 = 2.5; y3 = 2.5;
    }
    if (etas[ieta]=="Endcap2") {
      y1 = 2.5; y2 = 3.0; y3 = 3.0;
    }

    for (unsigned int itrg = 0; itrg != trgs.size(); ++itrg) {

      string sst = trgs[itrg];
      const char *st = trgs[itrg].c_str();
      cout << endl << "Trigger " << st << endl;

      TH1D *hn = (TH1D*)fin->Get(Form("Runs%s/npvgood_%s",se,st)); assert(hn);

      TH1D *hchf = (TH1D*)fin->Get(Form("Runs%s/c_chf%s_%s",se,sm,st));
      assert(hchf);
      TH1D *hnef = (TH1D*)fin->Get(Form("Runs%s/c_nef%s_%s",se,sm,st));
      assert(hnef);
      TH1D *hnhf = (TH1D*)fin->Get(Form("Runs%s/c_nhf%s_%s",se,sm,st));
      assert(hnhf);

      TH1D *hmcchf0 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hchf%s",
					   y1,y2,st,sm));
      assert(hmcchf0);
      if (y3!=y2) {
      TH1D *hmcchf1 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hchf%s",
					   y2,y3,st,sm));
      assert(hmcchf1);
      hmcchf0->Add(hmcchf1);
      }
      mcchf = hmcchf0->GetMean();
      cout << "mcchf : " << mcchf << endl;
      if (sst=="jt140")
	fout << Form("        refchf = %1.4f;", mcchf) << endl;

      TH1D *hmcnef0 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hnef%s",
					   y1,y2,st,sm));
      assert(hmcnef0);
      if (y3!=y2) {
      TH1D *hmcnef1 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hnef%s",
					   y2,y3,st,sm));
      assert(hmcnef1);
      hmcnef0->Add(hmcnef1);
      }
      mcnef = hmcnef0->GetMean();
      cout << "mcnef : " << mcnef << endl;
      if (sst=="jt140")
	fout << Form("        refnef = %1.4f;", mcnef) << endl;

      TH1D *hmcnhf0 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hnhf%s",
					   y1,y2,st,sm));
      assert(hmcnhf0);
      if (y3!=y2) {
      TH1D *hmcnhf1 = (TH1D*)fmc->Get(Form("Standard/Eta_%1.1f-%1.1f/%s/hnhf%s",
					   y2,y3,st,sm));
      assert(hmcnhf1);
      hmcnhf0->Add(hmcnhf1);
      }
      mcnhf = hmcnhf0->GetMean();
      cout << "mcnhf : " << mcnhf << endl;
      if (sst=="jt140")
	fout << Form("        refnhf = %1.4f;", mcnhf) << endl;

      // Clean data
      for (int i = 1; i != hchf->GetNbinsX()+1; ++i) {
	if (fabs(hchf->GetBinContent(i)/refchf-1) > 0.5 ||
	    fabs(hnef->GetBinContent(i)/refnef-1) > 0.5 ||
	    fabs(hnhf->GetBinContent(i)/refnhf-1) > 0.5) {
	  hchf->SetBinContent(i,0); hchf->SetBinError(i,0);
	  hnef->SetBinContent(i,0); hnef->SetBinError(i,0);
	  hnhf->SetBinContent(i,0); hnhf->SetBinError(i,0);
	}
      } // for i

      TGraphErrors *gchf0 = new TGraphErrors(hchf);
      TGraphErrors *gnef0 = new TGraphErrors(hnef);
      TGraphErrors *gnhf0 = new TGraphErrors(hnhf);

      TGraphErrors *gchfvsnpv = new TGraphErrors(0);
      TGraphErrors *gnefvsnpv = new TGraphErrors(0);
      TGraphErrors *gnhfvsnpv = new TGraphErrors(0);
      for (int i = 1; i != hnhf->GetNbinsX()+1; ++i) {

	if (hchf->GetBinContent(i)!=0 &&
	    hchf->GetBinError(i)<0.1*hchf->GetBinContent(i)) {
	  int n = gchfvsnpv->GetN();
	  gchfvsnpv->SetPoint(n, hn->GetBinContent(i), hchf->GetBinContent(i));
	  gchfvsnpv->SetPointError(n, 0., hchf->GetBinError(i));
	}

	if (hnef->GetBinContent(i)!=0 &&
	    hnef->GetBinError(i)<0.1*hnef->GetBinContent(i)) {
	  int n = gnefvsnpv->GetN();
	  gnefvsnpv->SetPoint(n, hn->GetBinContent(i), hnef->GetBinContent(i));
	  gnefvsnpv->SetPointError(n, 0., hnef->GetBinError(i));
	}

	if (hnhf->GetBinContent(i)!=0 &&
	    hnhf->GetBinError(i)<0.1*hnhf->GetBinContent(i)) {
	  int n = gnhfvsnpv->GetN();
	  gnhfvsnpv->SetPoint(n, hn->GetBinContent(i), hnhf->GetBinContent(i));
	  gnhfvsnpv->SetPointError(n, 0., hnhf->GetBinError(i));
	}

      } // for i

      c1->cd();
      h->DrawClone("AXIS");

      f0->SetLineWidth(2);

      hchf->Fit(f0,"QRN");
      f0->SetLineColor(kRed);
      f0->DrawClone("SAME");
      hnef->Fit(f0,"QRN");
      f0->SetLineColor(kBlue);
      f0->DrawClone("SAME");
      hnhf->Fit(f0,"QRN");
      f0->SetLineColor(kGreen+2);
      f0->DrawClone("SAME");

      hchf->Fit(f1,"QRN");
      f1->SetLineColor(kRed);
      f1->DrawClone("SAME");
      hnef->Fit(f1,"QRN");
      f1->SetLineColor(kBlue);
      f1->DrawClone("SAME");
      hnhf->Fit(f1,"QRN");
      f1->SetLineColor(kGreen+2);
      f1->DrawClone("SAME");

      f0->SetLineStyle(kDashed);

      f0->SetLineColor(kRed+1);
      f0->SetParameter(0, mcchf);
      f0->DrawClone("SAME");
      f0->SetLineColor(kBlue+1);
      f0->SetParameter(0, mcnef);
      f0->DrawClone("SAME");
      f0->SetLineColor(kGreen+3);
      f0->SetParameter(0, mcnhf);
      f0->DrawClone("SAME");

      f0->SetLineStyle(kDotted);
      f0->SetRange(ifirst, ilast);

      hchf->Fit(f0,"QRN");
      f0->SetLineColor(kRed);
      f0->DrawClone("SAME");
      hnef->Fit(f0,"QRN");
      f0->SetLineColor(kBlue);
      f0->DrawClone("SAME");
      hnhf->Fit(f0,"QRN");
      f0->SetLineColor(kGreen+2);
      f0->DrawClone("SAME");

      hchf->SetLineColor(kRed);
      hchf->SetMarkerColor(kRed);
      hchf->SetMarkerSize(0.7);
      hchf->SetMarkerStyle(kFullCircle);
      hchf->Draw("SAME");

      hnef->SetLineColor(kBlue);
      hnef->SetMarkerColor(kBlue);
      hnef->SetMarkerSize(0.7);
      hnef->SetMarkerStyle(kOpenCircle);
      hnef->Draw("SAME");

      hnhf->SetLineColor(kGreen+2);
      hnhf->SetMarkerColor(kGreen+2);
      hnhf->SetMarkerSize(0.7);
      hnhf->SetMarkerStyle(kFullSquare);
      hnhf->Draw("SAME");

      //leg->SetHeader("2#leq|#eta|<3, 153#leqp_{T}<196 GeV");
      leg->SetHeader(Form("%1.1f #leq |#eta| < %1.1f, %s",y1,y3,st));
      if (itrg==0) {
	leg->AddEntry(hchf,"chf","PL");
	leg->AddEntry(hnef,"nef","PL");
	leg->AddEntry(hnhf,"nhf","PL");
      }
      leg->Draw();

      if (jp::pdf) c1->SaveAs(Form("pdf/RunComposition_PF_%s_%s.pdf",se,st));

      // Coarse grain vs run
      TF1 *fit = new TF1("fit","[0]",0,350);
      TGraphErrors *gchf = new TGraphErrors(nrunbins);
      TGraphErrors *gnef = new TGraphErrors(nrunbins);
      TGraphErrors *gnhf = new TGraphErrors(nrunbins);
      gchf->SetName(Form("gchf%s_%s",se,st));
      gnef->SetName(Form("gnef%s_%s",se,st));
      gnhf->SetName(Form("gnhf%s_%s",se,st));
      for (int irun = 0; irun != nrunbins; ++irun) {

	double midrun = 0.5*(runbins[irun]+runbins[irun+1]);
	double drun = 0.5*(runbins[irun+1]-runbins[irun]);
	fit->SetRange(runbins[irun]+0.5, runbins[irun+1]+0.5);

	gchf0->Fit(fit,"QRN");
	double kchf = refchf / mcchf;
	double k = 1.;//mcchf/fit->GetParameter(0); // XTRA
	//double k0 = 1. / (1. + (mcchf - fit->GetParameter(0))); // XTRA
	gchf->SetPoint(irun, midrun, k*kchf*fit->GetParameter(0));
	gchf->SetPointError(irun, drun, k*kchf*fit->GetParError(0));

	gnef0->Fit(fit,"QRN");
	double knef = refnef / mcnef;
	gnef->SetPoint(irun, midrun, k*knef*fit->GetParameter(0));
	gnef->SetPointError(irun, drun, k*knef*fit->GetParError(0));

	gnhf0->Fit(fit,"QRN");
	double knhf = refnhf / mcnhf;
	gnhf->SetPoint(irun, midrun, k*knhf*fit->GetParameter(0));
	gnhf->SetPointError(irun, drun, k*knhf*fit->GetParError(0));

      } // for irun

      mchf->Add(gchf);
      mnef->Add(gnef);
      mnhf->Add(gnhf);

      c2->cd();
      if (itrg==0) {

	h->DrawClone("AXIS");

	fit->SetLineStyle(kDashed);
	fit->SetLineWidth(2);
	fit->SetRange(ifirst, ilast);
	fit->SetLineColor(kRed+1);
	fit->SetParameter(0, refchf);
	fit->DrawClone("SAME");
	fit->SetLineColor(kBlue+1);
	fit->SetParameter(0, refnef);
	fit->DrawClone("SAME");
	fit->SetLineColor(kGreen+3);
	fit->SetParameter(0, refnhf);
	fit->DrawClone("SAME");
      }

      gchf->SetMarkerStyle(kFullCircle);
      gchf->SetMarkerColor(kRed);
      gchf->SetLineColor(kRed);
      gchf->Draw("SAME P");

      gnef->SetMarkerStyle(kOpenCircle);
      gnef->SetMarkerColor(kBlue);
      gnef->SetLineColor(kBlue);
      gnef->Draw("SAME P");

      gnhf->SetMarkerStyle(kFullSquare);
      gnhf->SetMarkerColor(kGreen+2);
      gnhf->SetLineColor(kGreen+2);
      gnhf->Draw("SAME P");

      TLine *l = new TLine();
      l->SetLineStyle(kDashed);
      //l->DrawLine(its,0.,its,0.9);
      l->DrawLine(i12b,0.,i12b,0.9);
      l->DrawLine(i12c,0.,i12c,0.9);
      l->DrawLine(i12d,0.,i12d,0.9);

      TLatex *tex = new TLatex();
      tex->SetTextSize(0.040);
      //tex->DrawLatex(ifirst+2,0.55,"Run2012A");
      //tex->DrawLatex(its+10,0.55,"Run2012B");
      tex->DrawLatex(ifirst+2,0.55,"A");
      tex->DrawLatex(i12b+2,0.55,"B");
      tex->DrawLatex(i12c+2,0.55,"C");
      tex->DrawLatex(i12d+2,0.55,"D");

      leg->SetHeader(Form("%1.1f #leq |#eta| < %1.1f",y1,y3));
      if (itrg==0) leg->Draw();

      if (ieta==0 && trgs[itrg]=="jt30") {
	TCanvas *c4 = new TCanvas("c4","c4",600,600);
	gnhfvsnpv->Draw("AP");
	//gnefvsnpv->Draw("AP");
	//gchfvsnpv->Draw("AP");
      }

    } // for itrg

    // Coarse grain vs trigger (and run)
    TF1 *fit = new TF1("mfit","[0]",0,350);
    TGraphErrors *gchf = new TGraphErrors(nrunbins);
    TGraphErrors *gnef = new TGraphErrors(nrunbins);
    TGraphErrors *gnhf = new TGraphErrors(nrunbins);
    for (int irun = 0; irun != nrunbins; ++irun) {

      double midrun = 0.5*(runbins[irun]+runbins[irun+1]);
      double drun = 0.5*(runbins[irun+1]-runbins[irun]);
      fit->SetRange(runbins[irun]+0.5, runbins[irun+1]+0.5);

      mchf->Fit(fit,"QRN");
      gchf->SetPoint(irun, midrun, fit->GetParameter(0));
      gchf->SetPointError(irun, drun, fit->GetParError(0));

      mnef->Fit(fit,"QRN");
      gnef->SetPoint(irun, midrun, fit->GetParameter(0));
      gnef->SetPointError(irun, drun, fit->GetParError(0));

      mnhf->Fit(fit,"QRN");
      gnhf->SetPoint(irun, midrun, fit->GetParameter(0));
      gnhf->SetPointError(irun, drun, fit->GetParError(0));

    } // for irun

    c3->cd();
    h->DrawClone("AXIS");

    // MC reference
    fit->SetLineStyle(kDashed);
    fit->SetLineWidth(2);
    fit->SetRange(ifirst, ilast);
    fit->SetLineColor(kRed+1);
    fit->SetParameter(0, refchf);
    fit->DrawClone("SAME");
    fit->SetLineColor(kBlue+1);
    fit->SetParameter(0, refnef);
    fit->DrawClone("SAME");
    fit->SetLineColor(kGreen+3);
    fit->SetParameter(0, refnhf);
    fit->DrawClone("SAME");

    // Data full range (without v4)
    fit->SetLineWidth(1);
    fit->SetLineStyle(kSolid);
    fit->SetRange(ifirst, ilast);

    gchf->Fit(fit,"QRN");
    fit->SetLineColor(kRed);
    fit->DrawClone("SAME");
    gnef->Fit(fit,"QRN");
    fit->SetLineColor(kBlue);
    fit->DrawClone("SAME");
    gnhf->Fit(fit,"QRN");
    fit->SetLineColor(kGreen+2);
    fit->DrawClone("SAME");

    gchf->Fit(fit,"QRN");
    fit->SetLineColor(kRed);
    fit->DrawClone("SAME");
    gnef->Fit(fit,"QRN");
    fit->SetLineColor(kBlue);
    fit->DrawClone("SAME");
    gnhf->Fit(fit,"QRN");
    fit->SetLineColor(kGreen+2);
    fit->DrawClone("SAME");

    gchf->SetMarkerStyle(kFullCircle);
    gchf->SetMarkerColor(kRed);
    gchf->SetLineColor(kRed);
    gchf->Draw("SAME P");

    gnef->SetMarkerStyle(kOpenCircle);
    gnef->SetMarkerColor(kBlue);
    gnef->SetLineColor(kBlue);
    gnef->Draw("SAME P");

    gnhf->SetMarkerStyle(kFullSquare);
    gnhf->SetMarkerColor(kGreen+2);
    gnhf->SetLineColor(kGreen+2);
    gnhf->Draw("SAME P");


    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    //l->DrawLine(its,0.,its,0.9);
    l->DrawLine(i12b,0.,i12b,0.9);
    l->DrawLine(i12c,0.,i12c,0.9);
    l->DrawLine(i12d,0.,i12d,0.9);

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.040);
    //tex->DrawLatex(ifirst+2,0.55,"Run2012A");
    //tex->DrawLatex(its+2,0.55,"Run2012B");
    tex->DrawLatex(ifirst+2,0.55,"A");
    tex->DrawLatex(i12b+2,0.55,"B");
    tex->DrawLatex(i12c+2,0.55,"C");
    tex->DrawLatex(i12d+2,0.55,"D");

    leg->SetHeader(Form("%1.1f #leq |#eta| < %1.1f",y1,y3));
    leg->Draw();

    //cmsPrel(_lumi);


    // Differences zoomed up
    TGraphErrors *dchf = (TGraphErrors*)gchf->Clone();
    TGraphErrors *dnef = (TGraphErrors*)gnef->Clone();
    TGraphErrors *dnhf = (TGraphErrors*)gnhf->Clone();
    for (int i = 0; i != gchf->GetN(); ++i) {
      double x, y, ex, ey;
      tools::GetPoint(gchf, i, x, y, ex, ey);
      tools::SetPoint(dchf, i, x, 100.*(y-refchf), ex, 100.*ey);
      //
      tools::GetPoint(gnef, i, x, y, ex, ey);
      tools::SetPoint(dnef, i, x, 100.*(y-refnef), ex, 100.*ey);
      //
      tools::GetPoint(gnhf, i, x, y, ex, ey);
      tools::SetPoint(dnhf, i, x, 100.*(y-refnhf), ex, 100.*ey);
    }
    // clean out odd points
    for (int i = dchf->GetN()-1; i != -1; --i) {
      //if (fabs(dchf->GetY()[i])>3.) dchf->RemovePoint(i);
    }
    for (int i = dnef->GetN()-1; i != -1; --i) {
      //if (fabs(dnef->GetY()[i])>3.) dnef->RemovePoint(i);
    }
    for (int i = dnhf->GetN()-1; i != -1; --i) {
      //if (fabs(dnhf->GetY()[i])>3.) dnhf->RemovePoint(i);
    }

    c5->cd();

    double ymax0 = h->GetMaximum();
    double ymin0 = h->GetMinimum();
    double ymin = (y3>2 ? -9 : -3);//(ieta==2 ? -7 : -1.0);
    double ymax = (y3>2 ? 9: 3);//(ieta==2 ? +7 : +1.0);
    h->SetMaximum(ymax);
    h->SetMinimum(ymin);
    h->GetYaxis()->SetTitle("Data-Ref for E fraction per run [%]");
    h->DrawClone("AXIS");
    h->SetMaximum(ymax0);
    h->SetMinimum(ymin0);

    l->SetLineStyle(kDotted);
    l->DrawLine(ifirst-0.5, 0, ilast+0.5, 0);
    l->SetLineStyle(kDashed);
    //l->DrawLine(its,ymin,its,ymax);
    //l->DrawLine(iht,ymin,iht,ymax);
    //l->DrawLine(ijump,ymin,ijump,ymax);
    l->DrawLine(i12b,ymin,i12b,ymax);
    l->DrawLine(i12c,ymin,i12c,ymax);
    l->DrawLine(i12d,ymin,i12d,ymax);

    TF1 *f5 = new TF1("f5", "[0]", ifirst-0.5, ilast+0.5);
    //TF1 *f5a = new TF1("f5a", "[0]", ifirst, its);
    TF1 *f5a = new TF1("f5a", "[0]", ifirst, i12b);
    f5a->SetLineStyle(kDashed);
    //TF1 *f5b = new TF1("f5b", "[0]", iht, ijump);
    TF1 *f5b = new TF1("f5b", "[0]", i12b, i12c);
    f5b->SetLineStyle(kDashed);
    //TF1 *f5c = new TF1("f5c", "[0]", ijump, ilast);
    TF1 *f5c = new TF1("f5c", "[0]", i12c, i12d);
    f5c->SetLineStyle(kDashed);
    TF1 *f5d = new TF1("f5d", "[0]", i12d, ilast);
    f5d->SetLineStyle(kDashed);

    dchf->SetMarkerStyle(kFullCircle);
    dchf->SetMarkerColor(kRed);
    dchf->SetLineColor(kRed);
    dchf->Draw("SAME P");

    dchf->Fit(f5,"QRN");
    f5->SetLineColor(dchf->GetLineColor());
    f5->DrawClone("SAME");
    dchf->Fit(f5a,"QRN");
    f5a->SetLineColor(dchf->GetLineColor());
    f5a->DrawClone("SAME");
    dchf->Fit(f5b,"QRN");
    f5b->SetLineColor(dchf->GetLineColor());
    f5b->DrawClone("SAME");
    dchf->Fit(f5c,"QRN");
    f5c->SetLineColor(dchf->GetLineColor());
    f5c->DrawClone("SAME");
    dchf->Fit(f5d,"QRN");
    f5d->SetLineColor(dchf->GetLineColor());
    f5d->DrawClone("SAME");

    dnef->SetMarkerStyle(kOpenCircle);
    dnef->SetMarkerColor(kBlue);
    dnef->SetLineColor(kBlue);
    dnef->Draw("SAME P");

    dnef->Fit(f5,"QRN");
    f5->SetLineColor(dnef->GetLineColor());
    f5->DrawClone("SAME");
    dnef->Fit(f5a,"QRN");
    f5a->SetLineColor(dnef->GetLineColor());
    f5a->DrawClone("SAME");
    dnef->Fit(f5b,"QRN");
    f5b->SetLineColor(dnef->GetLineColor());
    f5b->DrawClone("SAME");
    dnef->Fit(f5c,"QRN");
    f5c->SetLineColor(dnef->GetLineColor());
    f5c->DrawClone("SAME");
    dnef->Fit(f5d,"QRN");
    f5d->SetLineColor(dnef->GetLineColor());
    f5d->DrawClone("SAME");

    dnhf->SetMarkerStyle(kFullSquare);
    dnhf->SetMarkerColor(kGreen+2);
    dnhf->SetLineColor(kGreen+2);
    dnhf->Draw("SAME P");

    dnhf->Fit(f5,"QRN");
    f5->SetLineColor(dnhf->GetLineColor());
    f5->DrawClone("SAME");
    dnhf->Fit(f5a,"QRN");
    f5a->SetLineColor(dnhf->GetLineColor());
    f5a->DrawClone("SAME");
    dnhf->Fit(f5b,"QRN");
    f5b->SetLineColor(dnhf->GetLineColor());
    f5b->DrawClone("SAME");
    dnhf->Fit(f5c,"QRN");
    f5c->SetLineColor(dnhf->GetLineColor());
    f5c->DrawClone("SAME");
    dnhf->Fit(f5d,"QRN");
    f5d->SetLineColor(dnhf->GetLineColor());
    f5d->DrawClone("SAME");

    tex->SetTextSize(0.040);
    double y = ymin + 0.05*(ymax-ymin);
    //tex->DrawLatex(ifirst+10,y,"2012A");
    //tex->DrawLatex(its+1,y,"2012B");
    //tex->DrawLatex(iht+8,y,"2012B/JetHT");
    //tex->DrawLatex(ijump+8,y,"Post 3/fb");
    tex->DrawLatex(ifirst+2,y,"A");
    tex->DrawLatex(i12b+2,y,"B");
    tex->DrawLatex(i12c+2,y,"C");
    tex->DrawLatex(i12d+2,y,"D");

    leg->SetHeader(Form("%1.1f #leq |#eta| < %1.1f",y1,y3));
    leg->Draw();

    //cmsPrel(_lumi);


    if (jp::pdf) c2->SaveAs(Form("pdf/RunCompositionA_PF_%s.pdf",se));
    if (jp::pdf) c3->SaveAs(Form("pdf/RunCompositionB_PF_%s.pdf",se));
    if (jp::pdf) c5->SaveAs(Form("pdf/RunCompositionC_PF_%s.pdf",se));

  } // for ieta

  fout.close();
} // drawRunComposition

