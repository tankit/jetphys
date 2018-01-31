#ifndef DRAWTRIGFRACS_H
#define DRAWTRIGFRACS_H

#include "TFile.h"
#include "THStack.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"

#include "tools.h"
#include "tdrstyle_mod15.C"
#include "settings.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

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

// Quite static settings {
const char ctp[] = "tp"; // "tp" for tag and probe, "" for nothing (the latter is seldom used)

const constexpr char* _fracs[7] = {"betastar","beta","chf","nef","nhf","cef","muf"};
const constexpr bool _dobeta=false;
const constexpr bool _dofit=false;
bool _shiftJES = false;
// } Static settings

// CAUTION: We are using non-const global variables. Since this is just a macro, this is ok.
vector<pair<double,double> > _etas;
map<string, pair<int, int> > _style;
map<string, string> _name;

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

TF1 *_fjes(0);
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

#endif // DRAWTRIGFRACS_H
