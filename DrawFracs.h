#ifndef DRAWFRACS_H
#define DRAWFRACS_H

#include "TFile.h"
#include "TSystem.h"
#include "THStack.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TError.h"

#include "tdrstyle_mod18.C"
#include "settings.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

class Fracs {

public:
  /* mode:
   *  0 - Pt
   *  1 - trpu
   *  2 - npv
   *  3 - eta */
  void DrawFracs(unsigned mode);

private:
  TDirectory *_dmc;
  TDirectory *_ddt;

  void makeProfile(unsigned mode, TDirectory *dmc, TDirectory *ddt, string trg, double etamin = -5.0, double etamax = 5.0);

  TF1 *_fhb;
  TF1 *_fjes;

  // Vectors holding mode info variables.
  vector<string> _modes;
  vector<double> _rangemin;
  vector<double> _rangemax;
  vector<double> _h2min;
  vector<double> _h2max;

  // Quite static settings {
  string _tp;

  vector <string> _fracs;

  static const constexpr bool _dofit = false;
  static const constexpr bool _shiftJES = false;
  // } Static settings

  // CAUTION: We are using non-const global variables. Since this is just a macro, this is ok.
  vector<pair<double,double> > _etas;
  map<string, pair<int, int> > _style;
  map<string, string> _name;

  double jesShift(double x) { // shift by L3Residual
    return 1./_fjes->Eval(x);
  } // jesShift

  // Bool flags for mode
  bool _vspt;
  bool _vspu;
  bool _vsnpv;
  bool _vseta;
  bool _vsphi;
  bool _phipos;

  bool _pertrg;
  bool _usetriglines;

  string _sphi;

  string _plot_title;
  string _lumi_title;
  string savedir;
  string _mc_type;
  string _dt_type;
  string _dt_typeverb;
  string _mc_typeverb;

  vector<double> _x;

  string _savedir;
  int _order;

  Fracs() {}
public:
 /* mc_type:
  *  "MC" (Pythia, default)
  *  "HW" (Herwig)
  *  "DT" (For dt vs dt mode)
  *
  * dt_type:
  *  "DT" (Data sample, default)
  *  "MC" (Pythia, for mc vs mc)
  *  "HW" (Herwig, for mc vs mc)
  * order:
  *  0 standard
  *  1 reshuffled */
  Fracs(const char *mc_path, const char *dt_path, string plot_title, string lumi_title, string savedir, bool pertrg = false, string mc_type = "MC", string dt_type = "RunG", int order = 0, string tp = "tp") {
    _plot_title = plot_title;
    _lumi_title = lumi_title;
    _savedir = savedir;
    assert(order==0 or order==1 or order==2 or order==3);
    _order = order;
    // This decides the ordering from bottom to top. Beta and muf are ignored.
    if (order==0)
      _fracs = {"puf","chf","nef","nhf","cef","muf","hhf","hef"};
    else if (order==1)
      _fracs = {"cef","muf","nef","puf","chf","nhf","hhf","hef"};
    else if (order==2)
      _fracs = {"nhf","puf","chf","nef","cef","muf","hhf","hef"};
    else if (order==3)
      _fracs = {"puf","chf","nhf","nef","cef","muf","hhf","hef"};

    gSystem->MakeDirectory(_savedir.c_str());

    _mc_type = mc_type;
    _dt_type = dt_type;
    if (_mc_type=="MCNU") _mc_type="MC";
    else if (_mc_type=="HWNU") _mc_type="HW";
    else if (_mc_type=="PthatNU") _mc_type="Pthat";
    else if (_mc_type=="HTNU") _mc_type="HT";
    else if (_mc_type=="NU") _mc_type="MC";
    if (_dt_type=="MCNU") _dt_type="MC";
    else if (_dt_type=="HWNU") _dt_type="HW";
    else if (_dt_type=="PthatNU") _dt_type="Pthat";
    else if (_dt_type=="HTNU") _dt_type="HT";
    else if (_dt_type=="NU") _dt_type="MC";
    _mc_typeverb = _mc_type;
    _dt_typeverb = _dt_type;
    if (_mc_typeverb=="DT") _mc_typeverb = "Data";
    else if (_mc_typeverb=="MC") _mc_typeverb = "Pythia8";
    else if (_mc_typeverb=="HW") _mc_typeverb = "Herwig++";
    else if (_mc_typeverb=="Pthat") _mc_typeverb = "P8 (binned)";
    else if (_mc_typeverb=="HT") _mc_typeverb = "Madgraph+P8";
    if (_dt_typeverb=="DT") _dt_typeverb = "Data";
    else if (_dt_typeverb=="MC") _dt_typeverb = "Pythia8";
    else if (_dt_typeverb=="HW") _dt_typeverb = "Herwig++";
    else if (_dt_typeverb=="Pthat") _dt_typeverb = "P8 (binned)";
    else if (_dt_typeverb=="HT") _dt_typeverb = "Madgraph+P8";

    _pertrg = pertrg;

    bool all_DT = false;
    bool all_MC = false;
    if (_mc_type!="MC" and _mc_type!="HW" and _mc_type!="Pthat" and _mc_type!="HT")
      all_DT = true;
    if (_dt_type=="MC" or _dt_type=="HW" or _dt_type=="Pthat" or _dt_type=="HT")
      all_MC = true;
    assert(!(all_MC and all_DT));

    // Opening the requested files {
    const char *cdt_type = all_MC ? _dt_type.c_str() : "DATA";
    TFile *fdt = new TFile(Form("%soutput-%s-%s.root",dt_path,cdt_type,pertrg?"1":"2b"),"READ");
    assert(fdt and !fdt->IsZombie());
    bool enterdtstandard = fdt->cd("Standard"); assert(enterdtstandard);
    _ddt = gDirectory;
    cout << "Opened the data file: ";
    cout << fdt->GetName() << endl;

    const char *cmc_type = all_DT ? "DATA" : mc_type.c_str();
    TFile *fmc = new TFile(Form("%soutput-%s-%s.root",mc_path,cmc_type,pertrg?"1":"2b"),"READ");
    assert(fmc and !fdt->IsZombie());
    bool entermcstandard = fmc->cd("Standard"); assert(entermcstandard);
    _dmc = gDirectory;
    cout << "Opened the mc file: ";
    cout << fmc->GetName() << endl;
    // }

    _etas.push_back(make_pair<double, double>(0., 1.3));
    _etas.push_back(make_pair<double, double>(0.0, 0.5));
    _etas.push_back(make_pair<double, double>(0.5, 1.0));
    _etas.push_back(make_pair<double, double>(1.0, 1.5));
    _etas.push_back(make_pair<double, double>(1.5, 2.0));
    _etas.push_back(make_pair<double, double>(2.0, 2.5));
    _etas.push_back(make_pair<double, double>(2.5, 3.0));
    _etas.push_back(make_pair<double, double>(3.0, 3.2));
    _etas.push_back(make_pair<double, double>(3.2, 4.7));

    _style["puf"] = make_pair<int, int>(kRed+2, kOpenCircle);
    _style["chf"] = make_pair<int, int>(kRed, kFullCircle);
    _style["nef"] = make_pair<int, int>(kBlue, kFullSquare);
    _style["nhf"] = make_pair<int, int>(kGreen+1, kFullDiamond);
    _style["cef"] = make_pair<int, int>(kCyan+1, kOpenTriangleUp);
    _style["muf"] = make_pair<int, int>(kMagenta+1, kOpenTriangleDown);
    _style["hhf"] = make_pair<int, int>(kViolet+2, kOpenDiamond);
    _style["hef"] = make_pair<int, int>(kOrange+2, kOpenSquare);

    _name["puf"] = "Charged pileup (overlaid)";
    _name["chf"] = "Charged hadrons";
    _name["nef"] = "Photons";
    _name["nhf"] = "Neutral hadrons";
    _name["cef"] = "Electrons, muons";
    _name["hhf"] = "Forward hadrons";
    _name["hef"] = "Forward photons";


    _fjes = new TF1("fjes","([0]+[1]*100./3.*(max(0.,1.03091-0.051154*pow(max(150.,min(340.,x)),-0.154227))-max(0.,1.03091-0.051154*pow(208.,-0.154227))) )",10,4000);
    _fjes->SetParameters( 0.9881, 0.2440 ); // Fall15_25nsV2
    // Initialize SinglePionHCAL and PileUpPt shapes
    _fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
    _fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH

    _modes = {"","_vstrpu","_vsnpv","_vseta","_vsphi","_vsphi"};
    _rangemin = {15, 0.5, 0.5,-5,-3.15,-3.15};
    _rangemax = {3450, 38.5, 40.5, 5,3.15,3.15};
    _h2min = {-4+1e-5, -4.5+1e-5, -6+1e-5, -5+1e-5, -4.5, -4.5};
    _h2max = { 4-1e-5, 4.5-1e-5, 16-1e-5, 5-1e-5, 4.5, 4.5};

    _tp = tp; // "tp" for tag and probe, "" for nothing (the latter is seldom used)
    _usetriglines = true;
  }
};

#endif // DRAWFRACS_H

#ifndef JESFIT
#define JESFIT

Double_t jesFit(Double_t *x, Double_t *p, TF1 *fhb) {
  double pt = *x;
  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  double ptx = max(150.,min(340.,pt));

  return (p[0] + p[1]/3.*100*(fhb->Eval(ptx)-fhb->Eval(ptref)));
} // jesFit

#endif // JESFIT
