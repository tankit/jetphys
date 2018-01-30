// Purpose: Plot PF energy fractions vs pT/PU/NPV (Run II cleaned code)
// Author: mikko.voutilainen(at)nospam.cern.ch
// Created: Dec 11, 2015
#include "TFile.h"
#include "THStack.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"

#include "tools.h"
#include "tdrstyle_mod15.C"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

using namespace std;
using namespace tools;

// Everyday settings:

bool _shiftJES = false;//true;

TF1 *_fjes(0);

// Vectors holding changing values: these are indexed as 'mode' (see above the drawFracs -function).
const constexpr char* _modenames[4] = {"","_vstrpu","_vsnpv","_vseta"};
const constexpr double _rangemin[4] = {37, 0.5, 0.5,-5};
const constexpr double _rangemax[4] = {3450, 38.5, 40.5};
const constexpr double _h2min[4] = {-4+1e-5 -6, -6+1e-5, -6+1e-5, -6+1e-5};
const constexpr double _h2max[4] = {+4-1e-5 +6, +6+10-1e-5, +6+10-1e-5, +6+10-1e-5};
const constexpr unsigned _nmodes = 4;

const constexpr char* _fracs[5] = {"betastar","chf", "nef", "nhf", "cef"};
const constexpr bool _dobeta=false;
const constexpr unsigned _nfrac = 5;
//const constexpr char* _fracs[] = {"betastar","chf", "beta", "nef", "nhf", "cef"};
//const constexpr bool _dobeta=true;
//const constexpr int _nfrac = 6;

const char ctp[] = "tp"; // "tp" for tag and probe, "" for nothing

// Wide binning for cleaner graphs and less jitter
const constexpr double xw[] =
//{1, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
//97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
//507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
//1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};
//2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
//4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
{1, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468, 548,
  686, 846, 1032, 1248, 1588, 2000, 2500, 3103, 3450, 3637, 3832,
  4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int nbins = sizeof(xw)/sizeof(xw[0]);


// CAUTION: We are using non-const global variables. Since this is just a macro, this is ok.
vector<pair<double,double> > _etas;
map<string, pair<int, int> > _style;
map<string, string> _name;


/* mode:
  *  0 - Pt
  *  1 - trpu
  *  2 - npv
  *
  * mc_type:
  *  "MC" (Pythia, default)
  *  "HW" (Herwig)
  *  "DT" (For dt vs dt mode)
  *
  * dt_type:
  *  "DT" (Data sample, default)
  *  "MC" (Pythia, for mc vs mc)
  *  "HW" (Herwig, for mc vs mc)
  */
void drawFracs(unsigned mode, string mc_path="./", string dt_path="./", string plot_title="RunH", string savedir=".", string mc_type = "MC", string dt_type="DT");

void makeProfile(unsigned mode, TDirectory *dmc, TDirectory *ddt, string mc_type, string dt_type, string plot_title, string savedir, unsigned tagnro = 0, double etamin = -5.0, double etamax = 5.0);

double jesShift(double x) {
  // shift by L3Residual
  if (!_fjes) {
    _fjes = new TF1("fjes","([0]+[1]*100./3.*(max(0.,1.03091-0.051154*pow(max(150.,min(340.,x)),-0.154227))-max(0.,1.03091-0.051154*pow(208.,-0.154227))) )",10,4000);
    _fjes->SetParameters( 0.9881, 0.2440 ); // Fall15_25nsV2
  }
  double jec = 1./_fjes->Eval(x);

  return jec;
} // jesShift

TF1 *fhb(0);
Double_t jesFit(Double_t *x, Double_t *p) {
  double pt = *x;

  // Initialize SinglePionHCAL and PileUpPt shapes
  if (!fhb) { 
    fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
    fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  }

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  double ptx = max(150.,min(340.,pt));

  return (p[0] + p[1]/3.*100*(fhb->Eval(ptx)-fhb->Eval(ptref)));
} // jesFit

void drawFracs(unsigned mode, string mc_path, string dt_path, string plot_title, string savedir, string mc_type, string dt_type) {
  setTDRStyle();

  bool all_DT = false, all_MC = false;
  if (mc_type!="MC" and mc_type!="HW")
    all_DT = true;
  if (dt_type=="MC" and dt_type=="HW")
    all_MC = true;
  assert(!(all_MC and all_DT));
  assert(mode>=0 and mode<_nmodes);

  setTDRStyle();

  // Opening the requested files {
  string dt_file = all_MC ? dt_type : "DATA";
  dt_file = "output-" + dt_file + "-2b.root";
  cout << Form("%s%s",dt_path.c_str(),dt_file.c_str()) << endl;
  TFile *fdt = new TFile(Form("%s%s",dt_path.c_str(),dt_file.c_str()),"READ");

  assert(fdt and !fdt->IsZombie());
  assert(fdt->cd("Standard"));
  TDirectory *ddt = gDirectory;

  cout << "Opened the data file: ";
  cout << Form("%s%s",dt_path.c_str(),dt_file.c_str()) << endl;

  string mc_file = all_DT ? "DATA" : mc_type;
  mc_file = "output-" + mc_file + "-2b.root";
  TFile *fmc = new TFile(Form("%s%s",mc_path.c_str(),mc_file.c_str()),"READ");

  assert(fmc and !fdt->IsZombie());
  assert(fmc->cd("Standard"));
  TDirectory *dmc = gDirectory;

  cout << "Opened the mc file: ";
  cout << Form("%s%s",mc_path.c_str(),mc_file.c_str()) << endl;
  // }

  _etas.push_back(make_pair<double, double>(0., 1.3));
  //_etas.push_back(make_pair<double, double>(0.0, 0.5));
  //_etas.push_back(make_pair<double, double>(0.5, 1.0));
  _etas.push_back(make_pair<double, double>(1.0, 1.5));
  _etas.push_back(make_pair<double, double>(1.5, 2.0));
  _etas.push_back(make_pair<double, double>(2.0, 2.5));
  _etas.push_back(make_pair<double, double>(2.5, 3.0));
  //_etas.push_back(make_pair<double, double>(3.0, 3.2));
  //_etas.push_back(make_pair<double, double>(3.2, 4.7));

  _style["betastar"] = make_pair<int, int>(kRed+2, kOpenCircle);
  _style["chf"] = make_pair<int, int>(_dobeta ? kRed+1 : kRed, _dobeta ? kOpenDiamond : kFullCircle);
  _style["beta"] = make_pair<int, int>(kRed, kFullCircle);
  _style["nef"] = make_pair<int, int>(kBlue, kFullSquare);
  _style["nhf"] = make_pair<int, int>(kGreen+1, kFullDiamond);
  _style["cef"] = make_pair<int, int>(kCyan+1, kOpenTriangleUp);
  _style["muf"] = make_pair<int, int>(kMagenta+1, kOpenTriangleDown);

  _name["betastar"] = "Charged pile-up";
  _name["chf"] = (_dobeta ? "Charged unassoc." : "Charged hadrons");
  _name["beta"] = "Charged hadrons";
  _name["nef"] = "Photons";
  _name["nhf"] = "Neutral hadrons";
  //_name["cef"] = "Electrons";
  _name["muf"] = "Muons";
  _name["cef"] = "Electrons+muons";

  for (unsigned int ieta = 0; ieta != _etas.size(); ++ieta) { // Loop through eta slices
    makeProfile(mode, dmc, ddt, mc_type, dt_type, plot_title, savedir, ieta, _etas[ieta].first, _etas[ieta].second);
  } // For ieta
} // drawFracs


void makeProfile(unsigned mode, TDirectory *dmc, TDirectory *ddt, string mc_type, string dt_type, string plot_title, string savedir, unsigned tagnro, double eta1, double eta2) {
  // List of differences
  map<string, TH1D*> mdf;
  // Bool flags for mode
  bool _vspt = false, _vspu = false, _vsnpv = false, _vseta = false;
  if (mode==0)
    _vspt=true;
  else if (mode==1)
    _vspu=true;
  else if (mode==2)
    _vsnpv=true;
  else if (mode==3)
    _vseta=true;
  else {
    cout << "Problematic mode " << mode << endl;
    return;
  }

  // Build appropriate wide binning
  vector<double> x;
  if (_vspt) {
    for (int i = 0; i != nbins && xw[i]!=0; ++i) x.push_back(xw[i]);
  } else if (_vspu) {
    for (int i = 0; i <= 25; ++i) x.push_back(-0.5+2*i);
  } else if (_vsnpv) {
    for (int i = 0; i <= 25; ++i) x.push_back(-0.5+2*i);
  }
  const int nx = x.size()-1;

  THStack *hsdt = new THStack(Form("hsdt%d",tagnro),"stacked histograms");
  THStack *hsmc = new THStack(Form("hsmc%d",tagnro),"stacked histograms");
  THStack *hsdf = new THStack(Form("hsdf%d",tagnro),"differences");

  TH1D *h = new TH1D("h",";p_{T} (GeV);PF energy fractions",nx,&x[0]);
  TH1D *h2 = new TH1D("h2",Form(";p_{T} (GeV);%s-%s (%%)",dt_type.c_str(),mc_type.c_str()),nx,&x[0]);

  if (_vspt) {
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
  } else if (_vspu) {
    h->SetXTitle("TruePU");
    h2->SetXTitle("TruePU");
  } else {
    h->SetXTitle("N_{PV,good}");
    h2->SetXTitle("N_{PV,good}");
  }
  h->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h2->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h->SetMaximum(1-1e-5);
  h->SetMinimum(0+1e-5);
  h2->SetMinimum(_h2min[mode]);
  h2->SetMaximum(_h2max[mode]);

  lumi_13TeV = plot_title;
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,0);

  c1->cd(1);
  TLegend *leg = tdrLeg(0.20,0.23-0.05,0.50,0.53-0.05);

  const char *dirname = Form("Eta_%1.1f-%1.1f",eta1,eta2);

  TH1D *href = new TH1D("href","", nx, &x[0]);
  map<string,TH1D*> mcHistos;
  map<string,TH1D*> dtHistos;
  for (int jfrac = 0; jfrac < _nfrac; ++jfrac) {
    const char *hname0 = Form("p%s%s%s",_fracs[jfrac],ctp,_modenames[mode]);
    bool enterdtdir = ddt->cd(dirname); assert(enterdtdir);
    TProfile *pdt = dynamic_cast<TProfile*>(gDirectory->Get(hname0)); assert(pdt);
    if (!pdt) {
      cout << hname0 << " not found in (dt) " << gDirectory->GetName() << endl << flush;
      gDirectory->ls();
      assert(false);
    }
    pdt->SetName(Form("%s_dt",pdt->GetName()));
    dtHistos[_fracs[jfrac]] = (_vspt ? tools::Rebin(pdt, href) : pdt->ProjectionX(Form("%s_px",pdt->GetName())));
    pdt->Delete();

    bool entermcdir = dmc->cd(dirname);
    assert(entermcdir);
    TProfile *pmc = dynamic_cast<TProfile*>(gDirectory->Get(hname0)); assert(pmc);
    if (!pmc) {
      cout << hname0 << " not found in (mc) " << gDirectory->GetName() << endl << flush;
      gDirectory->ls();
      assert(false);
    }
    pmc->SetName(Form("%s_mc",pmc->GetName()));
    mcHistos[_fracs[jfrac]] = (_vspt ? tools::Rebin(pmc, href) : pmc->ProjectionX(Form("%s_px",pmc->GetName())));
    pmc->Delete();
  }

  for (int jfrac = 0; jfrac != _nfrac; ++jfrac) {
    string sf = _fracs[jfrac];
    const char *cf = _fracs[jfrac];
    const char *cpu = _modenames[mode];

    bool enterdtdir = ddt->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)); assert(enterdtdir);
    const char *hname1 = Form("p%s%s%s",cf,ctp,cpu);
    TProfile *pdt = (TProfile*)(gDirectory->Get(hname1));

    pdt->SetName(Form("%s_dt",pdt->GetName()));

    bool entermcdir = dmc->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2));
    assert(entermcdir);
    const char *hname2 = Form("p%s%s%s",cf,ctp,cpu);
    TProfile *pmc = (TProfile*)(gDirectory->Get(hname2)); assert(pmc);
    if (!pmc) cout << hname2 << " not found in " << gDirectory->GetName() << endl << flush;
    pmc->SetName(Form("%s_mc",pmc->GetName()));

    // Rebin histogram
    TH1D *hmc = (!_vspt ? pmc->ProjectionX() : tools::Rebin(pmc, href));
    TH1D *hdt = (!_vspt ? pdt->ProjectionX() : tools::Rebin(pdt, href));

    // Scale data by response to MC-equivalent energy fractions
    if (_shiftJES) {
      for (int i = 1; i != hdt->GetNbinsX()+1; ++i) {
        double jec = jesShift(hdt->GetBinCenter(i));
        hdt->SetBinContent(i, hdt->GetBinContent(i)/jec);
      } // for i
    } // _shiftJES

    if (sf=="cef") { // For cef, add muf
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)+hdt2->GetBinContent(i));
      delete hdt2;

      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)+hmc2->GetBinContent(i));
      delete hmc2;
    } // cef -> cef + muf
    else if (sf=="betastar") { // For betastar, multiply by chf
      assert(ddt->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pdt2 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pchf%s%s",ctp,cpu))); assert(pdt2);
      TH1D *hdt2 = (!_vspt ? pdt2->ProjectionX() : tools::Rebin(pdt2, href));
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)*hdt2->GetBinContent(i));
      delete hdt2;

      assert(dmc->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pmc2 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pchf%s%s",ctp,cpu))); assert(pmc2);
      TH1D *hmc2 = (!_vspt ? pmc2->ProjectionX() : tools::Rebin(pmc2, href));
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)*hmc2->GetBinContent(i));
      delete hmc2;
    } // betastar -> chf * betastar
    else if (sf=="beta") { // For beta, multiply by chf
      assert(ddt->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pdt2 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pchf%s%s",ctp,cpu))); assert(pdt2);
      TH1D *hdt2 = (_vspu ? pdt2->ProjectionX() : tools::Rebin(pdt2, href));
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)*hdt2->GetBinContent(i));
      delete hdt2;

      assert(dmc->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pmc2 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pchf%s%s",ctp,cpu))); assert(pmc2);
      TH1D *hmc2 = (_vspu ? pmc2->ProjectionX() : tools::Rebin(pmc2, href));
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)*hmc2->GetBinContent(i));
      delete hmc2;
    } // beta -> chf * beta
    else if (sf=="chf") { // For chf, multiply by (1-beta-betastar)
      assert(ddt->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pdt2 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pbeta%s%s",ctp,cpu))); assert(pdt2);
      TH1D *hdt2 = (_vspu ? pdt2->ProjectionX() : tools::Rebin(pdt2, href));
      TProfile *pdt3 = dynamic_cast<TProfile*>(gDirectory->Get(Form("pbetastar%s%s",ctp,cpu))); assert(pdt3);
      TH1D *hdt3 = (_vspu ? pdt3->ProjectionX() : tools::Rebin(pdt3, href));
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i) {
        hdt->SetBinContent(i, hdt->GetBinContent(i) * (1 - (_dobeta ? hdt2->GetBinContent(i) : 0) - hdt3->GetBinContent(i)));
      }
      delete hdt2;
      delete hdt3;

      assert(dmc->cd(Form("Eta_%1.1f-%1.1f",eta1,eta2)));
      TProfile *pmc2 = (TProfile*)(gDirectory->Get(Form("pbeta%s%s",ctp,cpu))); assert(pmc2);
      TH1D *hmc2 = (_vspu ? pmc2->ProjectionX() : tools::Rebin(pmc2, href));
      TProfile *pmc3 = (TProfile*)(gDirectory->Get(Form("pbetastar%s%s",ctp,cpu))); assert(pmc3);
      TH1D *hmc3 = (_vspu ? pmc3->ProjectionX() : tools::Rebin(pmc3, href));
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i) {
        hmc->SetBinContent(i, hmc->GetBinContent(i) * (1 - (_dobeta ? hmc2->GetBinContent(i) : 0) - hmc3->GetBinContent(i)));
      }
      delete hmc2;
      delete hmc3;
    } // chf -> chf*(1-beta-betastar)
    delete href;

    hmc->SetMarkerStyle(kNone);
    hmc->SetFillStyle(1001);
    hmc->SetFillColor(_style[cf].first - 7);
    hmc->SetLineColor(_style[cf].first + 1);
    hmc->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    hsmc->Add(hmc, "SAME H");

    hdt->SetFillStyle(1001); // for legend
    hdt->SetFillColor(_style[cf].first - 7); // for legend
    hdt->SetLineColor(_style[cf].first + 1);
    hdt->SetMarkerStyle(_style[cf].second);
    hdt->SetMarkerSize(sf=="nhf" or (sf=="chf" and _dobeta) ? 1.3 : 1.0);
    if (_vspt)
      hdt->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    hsdt->Add(hdt, "SAME P");

    // Then, do the difference
    TH1D *hdf = dynamic_cast<TH1D*>(hdt->Clone(Form("hdf%d",tagnro))); assert(hdt);
    hdf->Add(hdt, hmc, 100, -100);
    hdf->SetLineColor(_style[cf].first + 1);
    hdf->SetMarkerColor(_style[cf].first + 1);
    hdf->SetMarkerStyle(_style[cf].second);
    hdf->SetLineWidth(tagnro==0 ? 2 : 0);
    hdf->SetMarkerSize(tagnro==0 ? (sf=="nhf" or (sf=="chf" and _dobeta) ? 1.8 : 1.3) : 0.7);
    hsdf->Add(hdf, "SAME P");

    c1->cd(2);

    if (jfrac==0) {
      TLine *l = new TLine();
      l->DrawLine(_rangemin[mode], 0,_rangemax[mode], 0);
      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextSize(h2->GetYaxis()->GetLabelSize());
      tex->DrawLatex(0.17,0.80,Form("Anti-k_{T} R=0.4%s",_shiftJES ? ", shifted by JES" : ""));
    }
    hdf->Draw("SAME");
    leg->AddEntry(hdt,_name[cf].c_str(),"PF");
    mdf[sf] = hdf;
  } // for jfrac

  c1->cd(1);
  if (_vspt) gPad->SetLogx();
  hsmc->Draw("SAME");
  hsdt->Draw("SAME");
  leg->Draw("SAME"); // redraw
  gPad->RedrawAxis();

  c1->cd(2);
  if (_vspt) gPad->SetLogx();
  //hsdf->Draw("SAME");
  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/drawFracs_%1.1f-%1.1f%s%s%s%s.pdf", savedir.c_str(), eta1, eta2, _shiftJES ? "_shiftJES" : "",
                  _vsnpv ? "_vsNPV" : "",_vspu ?  "_vsTRPU" : "", _vseta ? "_vseta" : ""));

  // Estimate jet response slope by analyzing composition
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(h2->GetYaxis()->GetLabelSize()*0.7);
  if (_dobeta) {
    TF1 *fchf = new TF1("fchf",jesFit,40,3000,2);
    fchf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
    TH1D *hchfa = dynamic_cast<TH1D*>(mdf["beta"]->Clone("hcfha")); assert(hchfa);
    hchfa->Add(mdf["chf"]);
    hchfa->Add(mdf["betastar"]);
    hchfa->Fit(fchf,"QRN");
    fchf->SetLineColor(kRed+3);
    fchf->Draw("SAME");
    hchfa->SetMarkerStyle(kFullStar);
    hchfa->SetMarkerColor(kRed+3);
    hchfa->SetLineColor(kRed+3);
    hchfa->Draw("SAMEP");
  } else {
    TF1 *fchf = new TF1("fchf",jesFit,40,3000,2);
    fchf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
    mdf["chf"]->Fit(fchf,"QRN");
    fchf->SetLineColor(kRed+2);
    fchf->Draw("SAME");

    tex->SetTextColor(kRed);
    tex->DrawLatex(0.17,0.40,Form("%1.2f#pm%1.2f%%, #chi^2/NDF=%1.1f/%d", fchf->GetParameter(0), 
                                  fchf->GetParError(0), fchf->GetChisquare(), fchf->GetNDF()));
  }

  TF1 *fnhf = new TF1("fnhf",jesFit,40,3000,2);
  fnhf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
  mdf["nhf"]->Fit(fnhf,"QRN");
  fnhf->SetLineColor(kGreen+2);
  fnhf->Draw("SAME");

  TF1 *fnef = new TF1("fnef",jesFit,40,3000,2);
  fnef->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
  mdf["nef"]->Fit(fnef,"QRN");
  fnef->SetLineColor(kBlue+1);
  fnef->Draw("SAME");

  TH1D *hall = dynamic_cast<TH1D*>(mdf["chf"]->Clone("hall")); assert(hall);
  if (_dobeta) hall->Add(mdf["beta"]);
  hall->Add(mdf["betastar"]);
  hall->Add(mdf["nef"]);
  hall->Add(mdf["nhf"]);
  hall->Add(mdf["cef"]);
  hall->SetMarkerStyle(kFullCross);
  hall->SetMarkerColor(kBlack);

  h2->SetMaximum(+5);//+3.0);
  h2->SetMinimum(-5);//-1.5);
  if (tagnro==0) c1->SaveAs(Form("%s/drawFracs_WithFit%s%s%s%s.pdf", savedir.c_str(), _shiftJES ? "_shiftJES" : "",
                                 _vsnpv ? "_vsNPV" : "", _vspu ? "_vsTRPU" : "", _vseta ? "_vseta" : ""));
}
