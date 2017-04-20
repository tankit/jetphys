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

// Algorithm to use ("AK4PF" or "AK8PF")
string _jp_algo = "AK4PFchs";
// Data type ("DATA", "MC", or "HW")
string _jp_type = "DATA";
// In case of DATA, choose run ("RunB/C/D/E/Fearly/Flate/G/H")
string _jp_run = "RunG";
// Kostas stored UNCORRECTED four-vector, current status: CORRECTED
// HOX: this is a source of constant anxiety, should be rechecked from time to time
bool _jp_undojes = true;
// We can choose also not to apply the new jes onto a four-vector
bool _jp_redojes = true;
// For debugging
bool _jp_skipl2l3res = false;

// Number of events to process (-1 for all)
Long64_t _jp_nentries = 
-1; // all
//10; // debug
//100000; // short test run
//1000000; // for MC
// Number of events to skip from the beginning (for debugging)
Long64_t _jp_nskip = 0;

//// MC: PU profiles for data and MC
bool _jp_reweighPU = true;
string _jp_pudata = "pileup/pileup_DT.root";
string _jp_pumc   = "pileup/pu.root";
string _jp_prescalefile = "";//pileup/prescales74x.txt";

//// MC: Process pThatbins instead of flat sample
const bool _jp_pthatbins = true;
// Number of pthat bins
const unsigned int _jp_npthatbins = 14;
// The corresponding ranges
vector<double> _jp_pthatranges = // The last number is ~inf, the first has -1 to allow underflow
  {30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200,20000};
// The corresponding lumi
vector<double> _jp_pthatsigmas = // Arbitrary scale
  {140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
  //{140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
vector<unsigned int> _jp_pthatnevts = 
  {9699558,9948791,7742665,5748730,7838066,11701816,3959986,9628335,11915305,6992746,2477018,1584378,596904,391735};
// The filenames need to be given here and in mk_fillHistos, since ROOT is exceedingly stupid
vector<string> _jp_pthatfiles = {
  "MC/P825ns80X_Moriond17/QCD_Pt_30to50_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_50to80_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_80to120_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_120to170_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_170to300_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_300to470_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_470to600_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_600to800_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_800to1000_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_1000to1400_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_1400to1800_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_1800to2400_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_2400to3200_TuneCUETP8M_13TeV_pythia8.root",
  "MC/P825ns80X_Moriond17/QCD_Pt_3200toInf_TuneCUETP8M_13TeV_pythia8.root",
};

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
string _jp_triggers[_jp_ntrigger] =
  {"jt40",    "jt60",    "jt80",   "jt140",   "jt200",   "jt260",   "jt320",   "jt400",  "jt450"};
double _jp_trigthr[_jp_ntrigger] =
  {40,        60,        80,       140,       200,       260,       320,       400,      450};
double _jp_trigranges[_jp_ntrigger][2] =
  { {0,84}, {84,114}, {114,196}, {196,272}, {272,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[5,6], AK4
  //{ {0,84}, {84,114}, {114,174}, {174,245}, {245,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[3,4], AK4
// brilcalc lumi --hltpath "HLT_PFJet40_v*" -i [JSON]
bool _jp_usetriglumi = true; // use luminosity numbers below, in /ub
double _jp_triglumi[_jp_ntrigger] = // in /ub
//// 2016: ////
  // brilcalc lumi -i Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json --hltpath "HLT_PFJet*" -o brilcalc_lemibyls.csv --byls --minBiasXsec 69000
  {264854.065,719060.448,2734876.213,23973802.117,102890630.006,587903231.81,1754681848.318,5142568134.826,35863307558.002};

// Decide whether or not to simulate triggers from MC (this is slow)
bool _jp_domctrigsim = false;
// Use "mc" trigger for whole pT range instead of stiching triggers together
bool _jp_usemctrig = true;
string _jp_mctrig = "jt450"; // reference trigger (for PU profile)
// Unprescaled luminosity for plots
const double _jp_lumi = 36.810440678; // /fb from brilcalc for jt450
const double _jp_sqrts = 13000.; // GeV
const double _jp_emax = _jp_sqrts/2.; // GeV

//////////////////////////////////
// Additional analysis switches //
//////////////////////////////////

string _jp_jecgt = "Summer16_23Sep2016";
string _jp_jecvers = "V6"; // Summer16 // "V2" ; // Spring16

// Use Intervals-Of-Validity for JEC
const bool _jp_useIOV = true;
const unsigned int _jp_nIOV = 4;
string _jp_IOVnames[_jp_nIOV]
  {"BCD",    "EF",    "G",   "H"};
// Trigger IOVs: the 1 for -inf and 400000 for inf (currently)
double _jp_IOVranges[_jp_nIOV][2] =
  { {1,276811}, {276831,278801}, {278802,280385}, {280919,400000} }; // Spring/Summer16_23Sep2016

// Produce run-level histograms
const bool _jp_doRunHistos = false; // Set to false to save time
// Produce basic set of histograms
const bool _jp_doBasicHistos = true;
// Produce full-eta TH3 histograms
const bool _jp_doEtaHistos = true;
// Special reco/gen histos in mc
const bool _jp_doEtaHistosMcResponse = true;

// Correct for trigger efficiency based on MC
const bool _jp_dotrigeff = false;
 // Correct pT<114 GeV only, if above _jp_dotrigeff=true
const bool _jp_dotrigefflowptonly = false;
// Correct for time-dependence (prescales) in data
const bool _jp_dotimedep = false;
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

//// Draw againts HERAPDF1.7 instead of PDF4LHC
const bool _jp_herapdf = false;

// Produce plots
const bool _jp_pdf = true;

// Simple helper
bool _jp_isdt = _jp_type=="DATA";
bool _jp_ismc = !_jp_isdt;

#endif // __settings_h__
