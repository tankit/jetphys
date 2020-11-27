// Purpose: Combine different triggers into a single spectrum
// Authors:  mikko.voutilainen@cern.ch
//           hannu.siikonen@cern.ch
// Created: March 22, 2010
// Updated: November 21, 2018

// We throw a name to recurseCombFile and the function loops over all directories.
// As we reach a certain level, a new histogram is declared.
// Onto this histogram we add all the weighted products of the jtX subdirectories.
// The output directory is used as a temporary saving place for the histos, since pointers in ROOT are problematic.

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TRegexp.h"

#include <iostream>
#include <map>
#include <string>
#include <set>

#include "settings.h"

using namespace std;

void recurseCombFile(TDirectory *indir, TDirectory *outdir, bool isdt, string hname = "hpt", bool ptSelect = true,
                 bool othProf = false, int lvl = -1, double etamid = 0);
map<string, pair<double, double> > _ptranges;
vector<string> _ignoretrgs;


// global variables (not pretty, but works)
TDirectory *_top = 0;

void HistosCombine(string type = "") {
  TDirectory *curdir = gDirectory;
  bool isdt = jp::isdt;
  if (type=="") type = jp::type;
  else          isdt = (type=="DATA");

  TFile *fin = new TFile(Form("output-%s-2a.root",type.c_str()),"READ");
  assert(fin and !fin->IsZombie());
  _top = gDirectory;

  TFile *fout = new TFile(Form("output-%s-2b.root",type.c_str()),"RECREATE");
  assert(fout and !fout->IsZombie());

  cout << "Calling HistosCombine("<<type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursions. These may take a few seconds" << endl << flush;

  // Store pT ranges to a nice map
  if (!isdt) {
    if (jp::usemctrig) {
      _ptranges["mc"] = pair<double,double>(0., jp::sqrts/2.);
    } else {
      _ptranges["mc"] = pair<double,double>(0., 0.);
      _ignoretrgs.push_back("mc");
    }
  }
  for (auto itrg = 0u; itrg != jp::notrigs; ++itrg) {
    _ptranges[jp::triggers[itrg]] = pair<double, double>(jp::trigranges[itrg][0], jp::trigranges[itrg][1]);

    if (!isdt and jp::usemctrig) { // When in mc, we need to know how to use trigger weighting. Either mc or triggers are ignored
      _ptranges[jp::triggers[itrg]] = pair<double,double>(0.,0.);
      _ignoretrgs.push_back(jp::triggers[itrg]);
    }
  }

  // Loop over all the directories recursively
  // List here the histograms that need merging
  recurseCombFile(fin, fout, isdt, "hpt");
  recurseCombFile(fin, fout, isdt, "hpt_pre");
  recurseCombFile(fin, fout, isdt, "hpt_notrigeff");
  recurseCombFile(fin, fout, isdt, "htrigeff");
  recurseCombFile(fin, fout, isdt, "htrigeffmc");
  recurseCombFile(fin, fout, isdt, "htrigeffsf");
  recurseCombFile(fin, fout, isdt, "hpt_notimedep");
  recurseCombFile(fin, fout, isdt, "hpt_withtimedep");
  recurseCombFile(fin, fout, isdt, "htimedep");
  recurseCombFile(fin, fout, isdt, "htimefit");
  if (isdt) recurseCombFile(fin, fout, isdt, "hlumi");
  if (isdt) recurseCombFile(fin, fout, isdt, "hlumi_orig");

  recurseCombFile(fin, fout, isdt, "hpt_evt");
  recurseCombFile(fin, fout, isdt, "hpt_evtcount");
  recurseCombFile(fin, fout, isdt, "hpt_jet");

  if (!isdt) recurseCombFile(fin, fout, isdt, "hpt_g0tw");

  recurseCombFile(fin, fout, isdt, "hpt0");
  recurseCombFile(fin, fout, isdt, "hpt1");
  recurseCombFile(fin, fout, isdt, "hpt2");
  recurseCombFile(fin, fout, isdt, "hpt3");
  recurseCombFile(fin, fout, isdt, "hpt_noid");
  recurseCombFile(fin, fout, isdt, "hpt_noevtid");
  recurseCombFile(fin, fout, isdt, "hpt_nojetid");

  recurseCombFile(fin, fout, isdt, "pa");
  recurseCombFile(fin, fout, isdt, "pnpv");
  recurseCombFile(fin, fout, isdt, "pnpvall");
  recurseCombFile(fin, fout, isdt, "prho");
  recurseCombFile(fin, fout, isdt, "prhovsnpv");
  recurseCombFile(fin, fout, isdt, "prhovstrpu");
  recurseCombFile(fin, fout, isdt, "pnpvvstrpu");
  recurseCombFile(fin, fout, isdt, "htrpu");
  recurseCombFile(fin, fout, isdt, "htrpu2");
  recurseCombFile(fin, fout, isdt, "pmjb_recoil");
  recurseCombFile(fin, fout, isdt, "pmjb_leading");
  recurseCombFile(fin, fout, isdt, "pmjbinv_leading");
  recurseCombFile(fin, fout, isdt, "pmpfinv_leading");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil2");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_unc");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_one");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_n");
  recurseCombFile(fin, fout, isdt, "pmpf_leading");
  recurseCombFile(fin, fout, isdt, "pmpf_leading2");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading2");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_unc");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_one");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_n");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading2");
  recurseCombFile(fin, fout, isdt, "pmjbinv_recoil");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil2");
  recurseCombFile(fin, fout, isdt, "pmjb_recoil1");
  recurseCombFile(fin, fout, isdt, "pmjb_leading1");
  recurseCombFile(fin, fout, isdt, "pmjbinv_leading1");
  recurseCombFile(fin, fout, isdt, "pmpfinv_leading1");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil1");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil21");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_unc1");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_one1");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_n1");
  recurseCombFile(fin, fout, isdt, "pmpf_leading1");
  recurseCombFile(fin, fout, isdt, "pmpf_leading21");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading1");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading21");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_unc1");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_one1");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_n1");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading1");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading21");
  recurseCombFile(fin, fout, isdt, "pmjbinv_recoil1");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil1");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil21");
  recurseCombFile(fin, fout, isdt, "pmjb_recoil0");
  recurseCombFile(fin, fout, isdt, "pmjb_leading0");
  recurseCombFile(fin, fout, isdt, "pmjbinv_leading0");
  recurseCombFile(fin, fout, isdt, "pmpfinv_leading0");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil0");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil20");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_unc0");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_one0");
  recurseCombFile(fin, fout, isdt, "pmpf_recoil_n0");
  recurseCombFile(fin, fout, isdt, "pmpf_leading0");
  recurseCombFile(fin, fout, isdt, "pmpf_leading20");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading0");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading20");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_unc0");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_one0");
  recurseCombFile(fin, fout, isdt, "pmpflead_leading_n0");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading0");
  recurseCombFile(fin, fout, isdt, "pmpfleadinv_leading20");
  recurseCombFile(fin, fout, isdt, "pmjbinv_recoil0");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil0");
  recurseCombFile(fin, fout, isdt, "pmpflead_recoil20");
  recurseCombFile(fin, fout, isdt, "h2mjb_recoil");
  recurseCombFile(fin, fout, isdt, "h2mjb_leading");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_leading");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading2");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil2");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_unc");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_one");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_n");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading2");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading2");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_unc");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_one");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_n");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading2");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_recoil");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil2");
  recurseCombFile(fin, fout, isdt, "h2mjb_recoil1");
  recurseCombFile(fin, fout, isdt, "h2mjb_leading1");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_leading1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading21");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil1");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil21");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_unc1");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_one1");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_n1");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading1");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading21");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading1");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading21");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_unc1");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_one1");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_n1");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading1");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading21");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_recoil1");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil1");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil21");
  recurseCombFile(fin, fout, isdt, "h2mjb_recoil0");
  recurseCombFile(fin, fout, isdt, "h2mjb_leading0");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_leading0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_leading20");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil0");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil20");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_unc0");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_one0");
  recurseCombFile(fin, fout, isdt, "h2mpf_recoil_n0");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading0");
  recurseCombFile(fin, fout, isdt, "h2mpf_leading20");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading0");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading20");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_unc0");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_one0");
  recurseCombFile(fin, fout, isdt, "h2mpflead_leading_n0");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading0");
  recurseCombFile(fin, fout, isdt, "h2mpfleadinv_leading20");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_recoil0");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil0");
  recurseCombFile(fin, fout, isdt, "h2mpflead_recoil20");
  recurseCombFile(fin, fout, isdt, "pmjb_ptave");
  recurseCombFile(fin, fout, isdt, "pmjbinv_ptave");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave2");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave2");
  recurseCombFile(fin, fout, isdt, "h2mjb_ptave");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_ptave");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave2");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave2");
  recurseCombFile(fin, fout, isdt, "pmjb_ptave1");
  recurseCombFile(fin, fout, isdt, "pmjbinv_ptave1");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave1");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave21");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave1");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave21");
  recurseCombFile(fin, fout, isdt, "h2mjb_ptave1");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_ptave1");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave1");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave21");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave21");
  recurseCombFile(fin, fout, isdt, "pmjb_ptave0");
  recurseCombFile(fin, fout, isdt, "pmjbinv_ptave0");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave0");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave20");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave0");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave20");
  recurseCombFile(fin, fout, isdt, "h2mjb_ptave0");
  recurseCombFile(fin, fout, isdt, "h2mjbinv_ptave0");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave0");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave20");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave20");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_unc");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_one");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_n");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_unc1");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_one1");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_n1");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_unc0");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_one0");
  recurseCombFile(fin, fout, isdt, "pmpf_ptave_n0");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_unc");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_one");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_n");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_unc1");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_one1");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_n1");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_unc0");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_one0");
  recurseCombFile(fin, fout, isdt, "h2mpf_ptave_n0");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_unc");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_one");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_n");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_unc1");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_one1");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_n1");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_unc0");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_one0");
  recurseCombFile(fin, fout, isdt, "pmpfinv_ptave_n0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_unc");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_one");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_n");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_unc1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_one1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_n1");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_unc0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_one0");
  recurseCombFile(fin, fout, isdt, "h2mpfinv_ptave_n0");
  recurseCombFile(fin, fout, isdt, "hCRecoil");
  recurseCombFile(fin, fout, isdt, "pCRecoil");
  recurseCombFile(fin, fout, isdt, "pCRecoil_leading");
  recurseCombFile(fin, fout, isdt, "pCRecoil_ptave");
  recurseCombFile(fin, fout, isdt, "hCRecoil1");
  recurseCombFile(fin, fout, isdt, "pCRecoil1");
  recurseCombFile(fin, fout, isdt, "pCRecoil_leading1");
  recurseCombFile(fin, fout, isdt, "pCRecoil_ptave1");
  recurseCombFile(fin, fout, isdt, "hCRecoil0");
  recurseCombFile(fin, fout, isdt, "pCRecoil0");
  recurseCombFile(fin, fout, isdt, "pCRecoil_leading0");
  recurseCombFile(fin, fout, isdt, "pCRecoil_ptave0");
  recurseCombFile(fin, fout, isdt, "prho_recoil");
  recurseCombFile(fin, fout, isdt, "prho_leading");
  recurseCombFile(fin, fout, isdt, "prho_ptave");
  recurseCombFile(fin, fout, isdt, "prho_recoil0");
  recurseCombFile(fin, fout, isdt, "prho_leading0");
  recurseCombFile(fin, fout, isdt, "prho_ptave0");
  recurseCombFile(fin, fout, isdt, "prho_recoil1");
  recurseCombFile(fin, fout, isdt, "prho_leading1");
  recurseCombFile(fin, fout, isdt, "prho_ptave1");
  recurseCombFile(fin, fout, isdt, "prho_recoil2");
  recurseCombFile(fin, fout, isdt, "prho_leading2");
  recurseCombFile(fin, fout, isdt, "prho_ptave2");
  recurseCombFile(fin, fout, isdt, "h2rho_recoil");
  recurseCombFile(fin, fout, isdt, "h2rho_leading");
  recurseCombFile(fin, fout, isdt, "h2rho_ptave");
  recurseCombFile(fin, fout, isdt, "h2rho_recoil0");
  recurseCombFile(fin, fout, isdt, "h2rho_leading0");
  recurseCombFile(fin, fout, isdt, "h2rho_ptave0");
  recurseCombFile(fin, fout, isdt, "h2rho_recoil1");
  recurseCombFile(fin, fout, isdt, "h2rho_leading1");
  recurseCombFile(fin, fout, isdt, "h2rho_ptave1");
  recurseCombFile(fin, fout, isdt, "h2rho_recoil2");
  recurseCombFile(fin, fout, isdt, "h2rho_leading2");
  recurseCombFile(fin, fout, isdt, "h2rho_ptave2");
  if(jp::ismc) {
    recurseCombFile(fin, fout, isdt, "hCRecoil_ptcl");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_leading");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_ptave");
    recurseCombFile(fin, fout, isdt, "pjtgenptsf");
    recurseCombFile(fin, fout, isdt, "pjtptsf");
    recurseCombFile(fin, fout, isdt, "pjtgenptsf_recoil");
    recurseCombFile(fin, fout, isdt, "pjtptsf_recoil");
    recurseCombFile(fin, fout, isdt, "hCRecoil_ptcl1");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl1");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl21");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_leading1");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_ptave1");
    recurseCombFile(fin, fout, isdt, "pjtgenptsf1");
    recurseCombFile(fin, fout, isdt, "pjtptsf1");
    recurseCombFile(fin, fout, isdt, "hCRecoil_ptcl0");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl0");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl20");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_leading0");
    recurseCombFile(fin, fout, isdt, "pCRecoil_ptcl2_ptave0");
    recurseCombFile(fin, fout, isdt, "pjtgenptsf0");
    recurseCombFile(fin, fout, isdt, "pjtptsf0");
  }
  //
  recurseCombFile(fin, fout, isdt, "pchs");
  //
  recurseCombFile(fin, fout, isdt, "pmpf");
  recurseCombFile(fin, fout, isdt, "pmpf1");
  recurseCombFile(fin, fout, isdt, "pmpf2");
  recurseCombFile(fin, fout, isdt, "pmpfinv");
  recurseCombFile(fin, fout, isdt, "pmpfinv1");
  recurseCombFile(fin, fout, isdt, "pmpfinv2");
  recurseCombFile(fin, fout, isdt, "pmpfx");
  recurseCombFile(fin, fout, isdt, "pmpfy");
  recurseCombFile(fin, fout, isdt, "pmpfz");
  recurseCombFile(fin, fout, isdt, "punc");
  //
  recurseCombFile(fin, fout, isdt, "pncand");
  recurseCombFile(fin, fout, isdt, "pnch");
  recurseCombFile(fin, fout, isdt, "pnne");
  recurseCombFile(fin, fout, isdt, "pnnh");
  recurseCombFile(fin, fout, isdt, "pnce");
  recurseCombFile(fin, fout, isdt, "pnmu");

  recurseCombFile(fin, fout, isdt, "pchf");
  recurseCombFile(fin, fout, isdt, "pnef");
  recurseCombFile(fin, fout, isdt, "pnhf");
  recurseCombFile(fin, fout, isdt, "pcef");
  recurseCombFile(fin, fout, isdt, "pmuf");
  recurseCombFile(fin, fout, isdt, "phhf");
  recurseCombFile(fin, fout, isdt, "phef");
  recurseCombFile(fin, fout, isdt, "pbeta");
  recurseCombFile(fin, fout, isdt, "pbetastar");
  recurseCombFile(fin, fout, isdt, "ppuf");
  recurseCombFile(fin, fout, isdt, "pjec");
  recurseCombFile(fin, fout, isdt, "pjec2");
  recurseCombFile(fin, fout, isdt, "pjec_res");
  //
  recurseCombFile(fin, fout, isdt, "pncandtp");
  recurseCombFile(fin, fout, isdt, "pnchtp");
  recurseCombFile(fin, fout, isdt, "pnnetp");
  recurseCombFile(fin, fout, isdt, "pnnhtp");
  recurseCombFile(fin, fout, isdt, "pncetp");
  recurseCombFile(fin, fout, isdt, "pnmutp");

  recurseCombFile(fin, fout, isdt, "pchftp");
  recurseCombFile(fin, fout, isdt, "pneftp");
  recurseCombFile(fin, fout, isdt, "pnhftp");
  recurseCombFile(fin, fout, isdt, "pceftp");
  recurseCombFile(fin, fout, isdt, "pmuftp");
  recurseCombFile(fin, fout, isdt, "phhftp");
  recurseCombFile(fin, fout, isdt, "pheftp");
  recurseCombFile(fin, fout, isdt, "pbetatp");
  recurseCombFile(fin, fout, isdt, "pbetastartp");
  recurseCombFile(fin, fout, isdt, "ppuftp");
  recurseCombFile(fin, fout, isdt, "pjectp");
  //
  recurseCombFile(fin, fout, isdt, "pchftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pneftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pnhftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pceftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pmuftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "phhftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pheftp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pbetatp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "pbetastartp_vsnpv", false, true);
  recurseCombFile(fin, fout, isdt, "ppuftp_vsnpv", false, true);
  //
  recurseCombFile(fin, fout, isdt, "pchftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pneftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pnhftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pceftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pmuftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "phhftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pheftp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pbetatp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "pbetastartp_vstrpu", false, true);
  recurseCombFile(fin, fout, isdt, "ppuftp_vstrpu", false, true);
  //
  recurseCombFile(fin, fout, isdt, "pchftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pneftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pnhftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pceftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pmuftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "phhftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pheftp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pbetatp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "pbetastartp_vseta", false, true);
  recurseCombFile(fin, fout, isdt, "ppuftp_vseta", false, true);
  //
  recurseCombFile(fin, fout, isdt, "pchfnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pnefnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pnhfnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pcefnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pmufnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "phhfnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "phefnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pbetanegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pbetastarnegtp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "ppufnegtp_vsphi", false, true);
  //
  recurseCombFile(fin, fout, isdt, "pchfpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pnefpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pnhfpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pcefpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pmufpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "phhfpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "phefpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pbetapostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "pbetastarpostp_vsphi", false, true);
  recurseCombFile(fin, fout, isdt, "ppufpostp_vsphi", false, true);

  recurseCombFile(fin, fout, isdt, "ppt_probepertag");

  recurseCombFile(fin, fout, isdt, "hdjasymm");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt");
  recurseCombFile(fin, fout, isdt, "hdjmpf");
  recurseCombFile(fin, fout, isdt, "hdjmpftp");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt");

  recurseCombFile(fin, fout, isdt, "hdjasymm_a005");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a005");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a005");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a005");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a005");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a005");
  recurseCombFile(fin, fout, isdt, "hdjasymm_a01");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a01");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a01");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a01");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a01");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a01");
  recurseCombFile(fin, fout, isdt, "hdjasymm_a015");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a015");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a015");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a015");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a015");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a015");
  recurseCombFile(fin, fout, isdt, "hdjasymm_a02");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a02");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a02");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a02");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a02");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a02");
  recurseCombFile(fin, fout, isdt, "hdjasymm_a025");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a025");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a025");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a025");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a025");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a025");
  recurseCombFile(fin, fout, isdt, "hdjasymm_a03");
  recurseCombFile(fin, fout, isdt, "hdjasymmtp_a03");
  recurseCombFile(fin, fout, isdt, "hdjasymmpt_a03");
  recurseCombFile(fin, fout, isdt, "hdjmpf_a03");
  recurseCombFile(fin, fout, isdt, "hdjmpftp_a03");
  recurseCombFile(fin, fout, isdt, "hdjmpfpt_a03");

  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a005");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a005");
  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a01");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a01");
  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a015");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a015");
  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a02");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a02");
  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a025");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a025");
  recurseCombFile(fin, fout, isdt, "hdjresp_tag_a03");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a03");

  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a005");
  recurseCombFile(fin, fout, isdt, "hdjresptp_probe_a005");
  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a01");
  recurseCombFile(fin, fout, isdt, "hdjresptp_probe_a01");
  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a015");
  recurseCombFile(fin, fout, isdt, "hdjresptp_probe_a015");
  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a02");
  recurseCombFile(fin, fout, isdt, "hdjresptp_probe_a02");
  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a025");
  recurseCombFile(fin, fout, isdt, "hdjresptp_probe_a025");
  recurseCombFile(fin, fout, isdt, "hdjresp_probe_a03");
  recurseCombFile(fin, fout, isdt, "hdjresptp_tag_a03");

  recurseCombFile(fin, fout, isdt, "hdjmass", false);
  recurseCombFile(fin, fout, isdt, "hdjmass0", false);

  if (isdt) recurseCombFile(fin, fout, isdt, "peff_new");

  curdir->cd();

  cout << endl;
  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();
  cout << "Output file closed" << endl;

  fin->Close();
  fin->Delete();
  cout << "Input file closed" << endl;
} // HistosNormalize

inline void binCp(TH1 *_hpt, TH1 *hpt, int bin) { // Copying a TH1 bin is simple
  _hpt->SetBinContent(bin, hpt->GetBinContent(bin));
  _hpt->SetBinError(bin, hpt->GetBinError(bin));
}

inline void binCp(TProfile *_hpt, TProfile *hpt, int bin) { // Copying a TProfile bin is complicated
  (*_hpt)[bin] += (*hpt)[bin]; // GetW
  (*_hpt->GetSumw2())[bin] += (*hpt->GetSumw2())[bin]; // GetW2
  _hpt->SetBinEntries(bin, hpt->GetBinEntries(bin)+_hpt->GetBinEntries(bin)); // GetB
  (*_hpt->GetBinSumw2())[bin] = (*hpt->GetBinSumw2())[bin]; // GetB2
}

template<typename T>
inline void fillHisto(T *hpt, TDirectory *outdir, TDirectory *indir, bool ptSelect, bool othProf) {
  const char* hname = hpt->GetName();
  T *_hpt = dynamic_cast<T*>(dynamic_cast<TObject*>(outdir->Get(hpt->GetName())));
  if (!_hpt) {
    cout << "Savings histogram fetch failed: " << hpt->GetName() << " in " << indir->GetName() << endl;
    if (jp::debug) cout << "Cloned _" << hname << endl;
    return;
  }

  if (ptSelect) {
    if (_ptranges.find(indir->GetName())==_ptranges.end())
      cout << "pT range not found for directory " << indir->GetName() << endl;
    assert(_ptranges.find(indir->GetName())!=_ptranges.end());
    double ptmin = _ptranges[indir->GetName()].first;
    double ptmax = _ptranges[indir->GetName()].second;

    bool isth3 = hpt->InheritsFrom("TH3");
    bool isth23 = isth3 or hpt->InheritsFrom("TH2");
    for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
      double pt = hpt->GetXaxis()->GetBinCenter(i);
      if (pt > ptmin and pt < ptmax) { // TODO: We could do better than a linear search!
        if (isth23) { // TH2 and TH3
          for (int j = 1; j != _hpt->GetNbinsY()+1; ++j) {
            if (isth3) {
              for (int k = 1; k != _hpt->GetNbinsZ()+1; ++k) {
                _hpt->SetBinContent(i,j,k, hpt->GetBinContent(i,j,k));
                _hpt->SetBinError(i,j,k, hpt->GetBinError(i,j,k));
              } // for k
            } else {
              _hpt->SetBinContent(i,j, hpt->GetBinContent(i,j));
              _hpt->SetBinError(i,j, hpt->GetBinError(i,j));
            }
          } // for j
        } else { // TH1 and TProfile
          binCp(_hpt,hpt,i);
        }
      } // in ptrange
    } // for i
  } else {
    _hpt->Add(hpt);
    outdir->cd();
  }
}


void recurseCombFile(TDirectory *indir, TDirectory *outdir, bool isdt, string hname, bool ptSelect, bool othProf, int lvl, double etamid)
{
  TDirectory *curdir = gDirectory;
  // Automatically go through the list of keys (directories)
  TObject *obj = 0;
  if (lvl==-1) {
    lvl++;
    cout << endl << hname << " ";
  }

  if (lvl<2) {
    TList *keys = indir->GetListOfKeys();
    TKey *key = 0;
    TListIter itkey(keys);
    while ((key = dynamic_cast<TKey*>(itkey.Next()))) {
      // We want the key to match to a directory or the selected histogram
      string classname = key->GetClassName();
      string kname = key->GetName();
      if (classname=="TDirectoryFile") { // All keys are looped over only in search for subdirs
        obj = key->ReadObj(); assert(obj);
      } else {
        continue;
      }
      if (obj->InheritsFrom("TDirectory")) {
        int loclvl = lvl;
        // Found a subdirectory: copy it to output and go deeper

        TDirectory *outdir2 = outdir;
        if (loclvl == 0) {
          if (outdir->FindKey(obj->GetName())==0) outdir->mkdir(obj->GetName());
          bool enterkeydir = outdir->cd(kname.c_str());
          assert(enterkeydir);
          outdir2 = outdir->GetDirectory(kname.c_str()); assert(outdir2);
          outdir2->cd();
          if (jp::debug) cout << kname.c_str() << endl;
        } else {
          if (jp::debug) cout << kname.c_str() << " (at bottom)" << endl;
          loclvl++; // Increase the level by one when we enter e.g. jt40
        }

        bool enterobjdir = indir->cd(obj->GetName());
        assert(enterobjdir);
        TDirectory *indir2 = indir->GetDirectory(obj->GetName()); assert(indir2);
        indir2->cd();

        // Check if directory name contains information on eta bin width
        // If yes, the next level is the bottom level with triggers
        // set flag and reset the combined histogram pointer
        float etamin, etamax;
        if (loclvl==0 and (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2) and (etamax>etamin) ) {
          outdir2->cd();
          TObject *inobj = indir2->Get(Form("%s/%s",isdt ? jp::reftrig : "mc",hname.c_str()));
          if (inobj) {
            cout << "[" << etamin << "," << etamax << "]";
            TH1 *hpt = dynamic_cast<TH1*>(inobj->Clone(hname.c_str()));
            hpt->Reset();
            if (hpt) {
              recurseCombFile(indir2, outdir2, isdt, hname, ptSelect, othProf, 1, 0.5*(etamin+etamax));
              hpt->Write();
              hpt->Delete();
              indir2->cd();
            }
          }
        } else if (loclvl==0 and TString(indir2->GetName()).Contains("FullEta")) {
          outdir2->cd();
          TObject *inobj = indir2->Get(Form("%s/%s",isdt ? jp::reftrig : "mc",hname.c_str()));
          if (inobj) {
            if (string(indir2->GetName())!="FullEta_Gen") cout << "FullEta";
            TH1 *hpt = dynamic_cast<TH1*>(inobj->Clone(hname.c_str()));
            hpt->Reset();
            if (hpt) {
              recurseCombFile(indir2, outdir2, isdt, hname, ptSelect, othProf, 1, etamid);
              hpt->Write();
              hpt->Delete();
              indir2->cd();
            }
          }
        } else {
          if (loclvl==0 or std::find(_ignoretrgs.begin(),_ignoretrgs.end(),indir2->GetName())==_ignoretrgs.end()) {
            string sid = indir->GetName();
            if (loclvl>0 and (sid=="Eta_0.0-1.3" or (TString(sid).Contains("FullEta") and sid!="FullEta_Gen"))) {
              if (string(indir2->GetName())=="jt40" or string(indir2->GetName())=="mc") cout << " ";
              cout << indir2->GetName();
              if (string(indir2->GetName())!=jp::reftrig) cout << ",";
            }
            recurseCombFile(indir2, outdir2, isdt, hname, ptSelect, othProf, loclvl, etamid);
          }
        }
      }
      obj->Delete();
    } // while key
  } else {
    obj = indir->Get(hname.c_str());
    if (obj) {
      // Copy over TProfile, TH3, TH2, TH1 histograms, in this precise order.
      if (obj->InheritsFrom("TProfile")) {
        TProfile *hpt = dynamic_cast<TProfile*>(obj);
        if (hpt) fillHisto(hpt, outdir, indir, ptSelect, othProf);
      } else {
        TH1 *hpt = 0;
        if (obj->InheritsFrom("TH3")) hpt = dynamic_cast<TH3D*>(obj);
        else if (obj->InheritsFrom("TH2")) hpt = dynamic_cast<TH2D*>(obj);
        else if (obj->InheritsFrom("TH1")) hpt = dynamic_cast<TH1D*>(obj);
        if (hpt) fillHisto(hpt, outdir, indir, ptSelect, othProf);
      }
      obj->Delete();
    }
  }
  curdir->cd();
} // recurseCombFile
