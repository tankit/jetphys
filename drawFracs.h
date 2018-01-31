#ifndef DRAWFRACS_H
#define DRAWFRACS_H

#include "TFile.h"
#include "TSystem.h"
#include "THStack.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"

#include "tools.h"
#include "tdrstyle_mod18.C"
#include "settings.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

using namespace tools;

class Fracs {

public:
  /* mode:
   *  0 - Pt
   *  1 - trpu
   *  2 - npv
   *  3 - eta */
  void drawFracs(unsigned mode);

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

  static const constexpr bool _dobeta = false;
  static const constexpr bool _dofit = false;
  static const constexpr bool _shiftJES = false;
  // } Static settings

  // CAUTION: We are using non-const global variables. Since this is just a macro, this is ok.
  vector<pair<double,double> > _etas;
  map<string, pair<int, int> > _style;
  map<string, string> _name;

  // Wide binning for cleaner graphs and less jitter
  vector <double> xw;

  double jesShift(double x) { // shift by L3Residual
    return 1./_fjes->Eval(x);
  } // jesShift

  // Bool flags for mode
  bool _vspt;
  bool _vspu;
  bool _vsnpv;
  bool _vseta;

  bool _pertrg;

  string _plot_title;
  string savedir;
  string _mc_type;
  string _dt_type;

  vector<double> _x;

  string _savedir;

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
  *  "HW" (Herwig, for mc vs mc) */
  Fracs(const char *mc_path, const char *dt_path, string plot_title, string savedir, bool pertrg = false, string mc_type = "MC", string dt_type = "RunG") {
    _plot_title = plot_title;
    _savedir = savedir;
    gSystem->MakeDirectory(_savedir.c_str());

    _mc_type = mc_type;
    _dt_type = dt_type;

    _pertrg = pertrg;
    if (pertrg) {
      for (auto &trg : _jp_triggers) gSystem->MakeDirectory(Form("%s/%s",_savedir.c_str(),trg));
    }

    bool all_DT = false;
    bool all_MC = false;
    if (_mc_type!="MC" and _mc_type!="HW")
      all_DT = true;
    if (_dt_type=="MC" and _dt_type=="HW")
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

    _fjes = new TF1("fjes","([0]+[1]*100./3.*(max(0.,1.03091-0.051154*pow(max(150.,min(340.,x)),-0.154227))-max(0.,1.03091-0.051154*pow(208.,-0.154227))) )",10,4000);
    _fjes->SetParameters( 0.9881, 0.2440 ); // Fall15_25nsV2
    // Initialize SinglePionHCAL and PileUpPt shapes
    _fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
    _fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH

    _modes = {"","_vstrpu","_vsnpv","_vseta"};
    _rangemin = {37, 0.5, 0.5,-5};
    _rangemax = {3450, 38.5, 40.5, 5};
    _h2min = {-10+1e-5, -6+1e-5, -6+1e-5, -6+1e-5};
    _h2max = { 10-1e-5, 16-1e-5, 16-1e-5, 16-1e-5};

    _fracs = {"betastar","beta","chf","nef","nhf","cef","muf"};

    xw =
    //{1, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    //97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    //507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    //1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};
    //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
    //4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
    {1, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468, 548,
      686, 846, 1032, 1248, 1588, 2000, 2500, 3103, 3450, 3637, 3832,
      4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

    _tp = "tp"; // "tp" for tag and probe, "" for nothing (the latter is seldom used)
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
