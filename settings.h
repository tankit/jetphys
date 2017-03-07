// Purpose: Analysis settings for jet physics package
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: June 1, 2015
#ifndef __settings_h__
#define __settings_h__
#include <string>
using std::string;

// All the settings variables are in global name space
// To avoid conflicts, and to make them easier to find
// global variables all start with _jp_ (for JetPhysics)

// Print out debugging information (for this only, no _jp_)
const bool _debug = false;

////////////////////////////////////////
// Necessary analysis settings/inputs //
////////////////////////////////////////

// Algorithm to use ("AK4" or "AK8")
std::string _jp_algo = "AK4";
// Data type ("DATA", "MC", or "HW")
std::string _jp_type = "DATA";

// Number of events to process (-1 for all)
Long64_t _jp_nentries = -1; // all
//Long64_t _jp_nentries = 10; // debug
//Long64_t _jp_nentries = 100000; // short test run
//Long64_t _jp_nentries = 1000000; // for MC
// Number of events to skip from the beginning (for debugging)
Long64_t _jp_nskip = 0;

// PU profiles for data and MC
bool _jp_reweighPU = true;
std::string _jp_pudata = "pileup/pileup_DT.root";
std::string _jp_pumc   = "pileup/pileup_PY_80X.root";
//std::string _jp_pumc   = "pileup/pileup_PY_74X.root";
//std::string _jp_prescalefile = "pileup/prescales74x_run.txt";
std::string _jp_prescalefile = "";//pileup/prescales74x.txt";

//// JEC text files in CondFormats/JetMETObjects/data/
bool _jp_redoJEC = true;
//// Due to the HIP problem and other intricacies in Run2-16, the datasets are handled separately.
//// Note: different corrections for PromptReco and ReReco {

/// PromptReco
std::string _jp_jecgt_dt = "Spring16_25nsV10BCD"; // IOV [1,276811] B,C,D
//std::string _jp_jecgt_dt = "Spring16_25nsV10E";   // IOV [276831,277420] E
//std::string _jp_jecgt_dt = "Spring16_25nsV10F";   // IOV [277772,278801] Fearly
//std::string _jp_jecgt_dt = "Spring16_25nsV10p2";  // IOV [278802,Infinity] Flate,G,H
std::string _jp_jecgt_mc = "Spring16_25nsV10";
/// ReReco
//std::string _jp_jecgt_dt = "Spring16_23Sep2016BCDV2"; // IOV [1,276811] B,C,D
//std::string _jp_jecgt_dt = "Spring16_23Sep2016EFV2";  // IOV [276831,278801] E,Fearly
//std::string _jp_jecgt_dt = "Spring16_23Sep2016GV2";   // IOV [278802,280385] Flate,G
//std::string _jp_jecgt_dt = "Spring16_23Sep2016HV2";   // IOV [280919,Infinity] H
//std::string _jp_jecgt_mc = "Spring16_23Sep2016V2";

/// Old format: no distinct dt filenames
//std::string _jp_jecgt = "Spring16_25nsV6";
//std::string _jp_jecgt = "Fall15_25nsV2"; // 76X
//std::string _jp_jecgt = "Summer15_25nsV7"; // 74X
//std::string _jp_jecgt = "Summer15_25nsV6";
//std::string _jp_jecgt = "Summer15_50nsV5";

// }

// Veto jets near ECAL boundaries in JetID
const bool _jp_doECALveto = false;
string _jp_ecalveto = "lumicalc/ecalveto.root";

// Reapply json selection based on the latest one (check lumicalc if false!)
const bool _jp_dojson = true;
// Here: there are slight differences between PromptReco and ReReco in the 2016 run
//string _jp_json = "lumicalc/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
string _jp_json = "lumicalc/Cert_271036-282092_13TeV_PromptReco_Collisions16_JSON.txt";

// Calculate luminosity on the fly based on .csv file
const bool _jp_dolumi = true;
string _jp_lumifile = "lumicalc/brilcalc_lumibyls.csv"; // Run II

// This is the 13 TeV 25 ns list (Run2016BCDEFG)
// Check the recommended settings from https://twiki.cern.ch/twiki/bin/view/CMS/InclusiveJetsLegacy
const int _jp_ntrigger = 9; // jt450 unprescaled, so drop 500, but add Zero Bias
std::string _jp_triggers[_jp_ntrigger] =
  {"jt40",    "jt60",    "jt80",   "jt140",   "jt200",   "jt260",   "jt320",   "jt400",  "jt450"};
double _jp_trigthr[_jp_ntrigger] =
  {40,       60,        80,       140,       200,       260,       320,       400,      450};
double _jp_trigranges[_jp_ntrigger][2] =
  //{ {0,74}, {74,114}, {114,174}, {174,245}, {245,300}, {300,395}, {395,468}, {468,507}, {507,6500} }; // V[3,4], AK4
  { {0,84}, {84,114}, {114,174}, {174,245}, {245,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[3,4], AK4
// brilcalc lumi --hltpath "HLT_PFJet40_v*" -i [JSON]
bool _jp_usetriglumi = true; // use luminosity numbers below, in /ub
double _jp_triglumi[_jp_ntrigger] = // in /ub
//// 2016: ////
  // brilcalc lumi -i Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json --hltpath "HLT_PFJet*" -o brilcalc_lemibyls.csv --byls --minBiasXsec 69000
  {271792.436,738011.349,2806215.978,24609118.02,105577219.6,603409186.2,1800913935,5277989761,36810440678};

// Decide whether or not to simulate triggers from MC (this is slow)
bool _jp_domctrigsim = false;
// Use "mc" trigger for whole pT range instead of stiching triggers together
bool _jp_usemctrig = true;
std::string _jp_mctrig = "jt450"; // reference trigger (for PU profile)
// Unprescaled luminosity for plots
const double _jp_lumi = 36.810440678; // /fb from brilcalc for jt450
const double _jp_sqrts = 13000.; // GeV
const double _jp_emax = _jp_sqrts/2.; // GeV

//////////////////////////////////
// Additional analysis switches //
//////////////////////////////////

// Use Intervals-Of-Validity for JEC
const bool _jp_useIOV = false;
// Produce run-level histograms
const bool _jp_doRunHistos = false; // Set to false to save time
// Produce basic set of histograms
const bool _jp_doBasicHistos = true;

// Process pThatbins instead of flat sample
const bool _jp_pthatbins = false;
// Number of pthat bins
const unsigned int _jp_npthatbins = 14;
// The corresponding ranges
double _jp_pthatranges[_jp_npthatbins+1] = // The last number is ~inf
{30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200,20000}}
// The corresponding lumi
double _jp_pthatinvlums[_jp_npthatbins] = // Arbitrary scale
{140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445}
//{140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445}

// Correct for trigger efficiency based on MC
const bool _jp_dotrigeff = false;//true;
 // Correct pT<114 GeV only, if above _jp_dotrigeff=true
const bool _jp_dotrigefflowptonly = false;//true;
// Correct for time-dependence (prescales) in data
const bool _jp_dotimedep = false;//true;
// For creating smearing matrix
const bool _jp_doMatrix = false;

// Check that run / lumi section was listed in the .csv file
const bool _jp_dolumcheck = false;
//// Veto bad run list
//const bool _jp_dorunveto = false;
// Check for duplicates (warning: takes a lot of memory!)
const bool _jp_checkduplicates = false;//true;
// Only load selected branches (large speedup, but be careful!)
const bool _jp_quick = true;
// Center uncertainties around ansatz (true) or data (false)
const bool _jp_centerOnAnsatz = false;
const bool _jp_centerOnTheory = false;

// Plot Pythia for final PRL results
const bool _plotPythia = false;
// Minimum and maximum pT range to be plotted and fitted
const double _jp_recopt = 15;
const double _jp_fitptmin = 43;
// Changed on 2013-05-020: analysis from 49 GeV to 56 GeV
const double _jp_xmin57 = 56;
const double _jp_xminpas = 56;
const double _jp_xmin = 24.;//20.;
const double _jp_xmax = 1999.;

const double _jp_xsecMinBias = 7.126E+10;

//// Write histograms to PAS file
//const bool _jp_pas = true;
// Draw b-jets against MC@NLO instead of reco MC
//const bool _mcnlo = true;
//// Draw againts HERAPDF1.7 instead of PDF4LHC
const bool _jp_herapdf = false;

// Produce plots
const bool _jp_pdf = true;

#endif // __settings_h__
