// Purpose: Analysis settings for jet physics package
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: June 1, 2015

#ifndef __settings_h__
#define __settings_h__
#include <map>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>
#include <cassert>

using std::map;
using std::vector;
using std::string;
using std::cout;
using std::endl;

// All the settings variables are in global name space
// To avoid conflicts, and to make them easier to find
// global variables all start with _jp_ (for JetPhysics)

// Print out debugging information (for this only, no _jp_)
constexpr bool _jt_debug = false;

////////////////////////////////////////
// Necessary analysis settings/inputs //
////////////////////////////////////////

constexpr bool strings_equal(char const * a, char const * b) {
  return *a == *b && (*a == '\0' || strings_equal(a + 1, b + 1));
}

//{ BEGIN fundamental file settings
// Do we use CHS jets? ("CHS" for yes, "", for legacy/no)
const constexpr char _jp_chs[] = "CHS";
// Algorithm to use ("AK4PF" or "AK8PF" + "chs" for chs jets)
const constexpr char _jp_algo[] = "AK4PFchs";
// Data type ("DATA", "MC", or "HW")
const constexpr char _jp_type[] = "DATA";
// In case of DATA, choose run ("RunB/C/D/E/Fearly/Flate/G/H")
const constexpr char _jp_run[] = "RunG16";
// Simple helper
constexpr bool _jp_isdt = strings_equal(_jp_type,"DATA");
constexpr bool _jp_ismc = !_jp_isdt;
constexpr bool _jp_ispy = strings_equal(_jp_type,"MC");
constexpr bool _jp_ishw = strings_equal(_jp_type,"HW");
//} END fundamental settings

//{ BEGIN run settings
// Number of events to process (-1 for all)
constexpr Long64_t _jp_nentries =
//-1; // all
//10; // debug
100000; // short test run
// Number of events to skip from the beginning (for debugging)
constexpr Long64_t _jp_nskip = 0;
constexpr Long64_t _js_skim = 0; // "prescale", off if zero
//} END run settings

//{ BEGIN Switches to decide what goes into fillhistos
// If the output files grow too much, processing becomes more difficult. 
// Produce run-level histograms
constexpr bool _jp_doRunHistos = false; // Not needed very often
// Produce basic set of histograms
constexpr bool _jp_doBasicHistos = true;
// Produce full-eta TH3 histograms
constexpr bool _jp_doEtaHistos = true;
// Special reco/gen histos in mc
constexpr bool _jp_doEtaHistosMcResponse = false; // Special studies
// Toggle 3D asymmetry and mpf histos off, if too much space is consumed
constexpr bool _jp_do3dHistos = true;
//} END 

//{ BEGIN JES switches and JEC settings
// At some point of time Kostas stored UNCORRECTED four-vector.
// Current status: CORRECTED
// NOTE: this is a source of constant anxiety, should be rechecked from time to time
constexpr bool _jp_undojes = true;
// We can choose also not to apply the new jes onto a four-vector
constexpr bool _jp_redojes = true;
// For debugging
constexpr bool _jp_skipl2l3res = false;

// https://github.com/cms-jet/JECDatabase/tree/master/tarballs
// Summer16_03Feb2017G_V7_DATA.tar.gz [BCD, EF, G, H]
const constexpr char _jp_jecgt[] = "Summer16_03Feb2017";//BCD_";//"Summer15_50ns";// "Summer16_23Sep2016";
const constexpr char _jp_jecvers[] = "_V7";//"V4"; // Summer16_03Feb // "V6"; // Summer16_23Sep // "V2" ; // Spring16

//const constexpr char* const days[] = 
//{ "mon", "tue", "wed", "thur","fri", "sat", "sun" };
// Use Intervals-Of-Validity for JEC
constexpr const bool _jp_useIOV = true ;//false
constexpr const unsigned int _jp_nIOV = 4;
const constexpr char* _jp_IOVnames[_jp_nIOV] =
  {"BCD",    "EF",    "G",   "H"};
// Trigger IOVs: the 1 for -inf and 400000 for inf (currently)
const constexpr double _jp_IOVranges[_jp_nIOV][2] =
  { {1,276811}, {276831,278801}, {278802,280385}, {280919,400000} }; // Spring/Summer16_23Sep2016
//} END JES and JEC

// QGL studies by Ozlem
constexpr bool _jp_doqgl = false;
const constexpr char _jp_qglfile[] = "output-DATA_RunG_part1_16oct-1.root";

//{ BEGIN MC: PU profiles for data and MC
constexpr bool _jp_reweighPU = false;//true
const constexpr char _jp_pudata[] = "pileup/pileup_DT.root";
const constexpr char _jp_pumc[]   = "pileup/pu.root";
const constexpr char _jp_prescalefile[] = "";//pileup/prescales74x.txt";
//} END PU profiles

//{ BEGIN MC: Process pThatbins instead of flat sample
// All the following values need to be updated, if other source tuples are used.
// These can be obtained using getSliceEvts.C
constexpr bool _jp_pthatbins = true;//false;
// Number of pthat bins
constexpr unsigned int _jp_npthatbins = 14;
// The corresponding ranges, the last number is ~inf
const constexpr double _jp_pthatranges[_jp_npthatbins+1] = 
  {30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200,20000};
// The corresponding lumis in an arbitrary scale
const constexpr double _jp_pthatsigmas[_jp_npthatbins+1] =
  {140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
  //{140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
const constexpr double _jp_pthatnevts[_jp_npthatbins+1] =
  {9699558,9948791,7742665,5748730,7838066,11701816,3959986,9628335,11915305,6992746,2477018,1584378,596904,391735};
// The filenames need to be given here and in mk_fillHistos, since ROOT is exceedingly stupid
const vector<string> _jp_pthatfiles[_jp_npthatbins] = {
    "P825ns80X_Moriond17/QCD_Pt_30to50_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_50to80_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_80to120_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_120to170_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_170to300_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_300to470_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_470to600_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_600to800_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_800to1000_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_1000to1400_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_1400to1800_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_1800to2400_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_2400to3200_TuneCUETP8M_13TeV_pythia8.root",
    "P825ns80X_Moriond17/QCD_Pt_3200toInf_TuneCUETP8M_13TeV_pythia8.root",
};
//} END pthatbins stuff

// Reapply json selection based on the latest one (check lumicalc if false!)
constexpr bool _jp_dojson = false; //true;
// Here: there are slight differences between PromptReco and ReReco in the 2016 run
//string _jp_json = "PromptReco/Cert_299614-299617_13TeV_PromptReco_Collisions17_50ns_JSON.txt";
const constexpr char _jp_json[] = "lumicalc/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
//string _jp_json = "lumicalc/Cert_271036-282092_13TeV_PromptReco_Collisions16_JSON.txt";

// Calculate luminosity on the fly based on .csv file
constexpr bool _jp_dolumi = false; //true;
const constexpr char _jp_lumifile[] = "lumicalc/brilcalc_lumibyls.csv"; // Run II

//{ BEGIN Trigger stuff
// Decide whether or not to simulate triggers from MC (this is slow)
constexpr bool _jp_domctrigsim = true;
// Use "mc" trigger for whole pT range instead of stiching triggers together (requires trigsim)
constexpr bool _jp_usemctrig = true;
// reference trigger (for PU profile) in the mc folder
const constexpr char _jp_mctrig[] = "jt450";

// This is the 13 TeV 25 ns list (Run2016BCDEFG)
// Check the recommended settings from https://twiki.cern.ch/twiki/bin/view/CMS/InclusiveJetsLegacy
constexpr int _jp_ntrigs = 9; // jt450 unprescaled, so drop 500, but add Zero Bias
const constexpr char* _jp_triggers[_jp_ntrigs] =
  {"jt40",    "jt60",    "jt80",   "jt140",   "jt200",   "jt260",   "jt320",   "jt400",  "jt450"};
const constexpr double _jp_trigthr[_jp_ntrigs] =
  {40,        60,        80,       140,       200,       260,       320,       400,      450};
const constexpr double _jp_trigranges[_jp_ntrigs][2] =
  { {0,84}, {84,114}, {114,196}, {196,272}, {272,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[5,6], AK4
  //{ {0,84}, {84,114}, {114,174}, {174,245}, {245,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[3,4], AK4
// brilcalc lumi --hltpath "HLT_PFJet40_v*" -i [JSON]
constexpr bool _jp_usetriglumi = true; // use luminosity numbers below, in /ub
const constexpr double _jp_triglumi[_jp_ntrigs] = // in /ub
//// 2016: ///
// brilcalc lumi -i Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json --hltpath "HLT_PFJet*" -o brilcalc_lumibyls.csv --byls --minBiasXsec 69000
// {264854.065,719060.448,2734876.213,23973802.117,102890630.006,587903231.81,1754681848.318,5142568134.826,35863307558.002}; // Ozlem, add new numbers here for RunG-part1_JSON.txt
  {264314.035000,718756.543000,2727540.789000,23949603.660000,102671554.653000,587533626.778000,1753875820.267000,5139583558.868000,35918235379.677002};//Ozlem added new numbers here for RunG-part1_JSON.txt
// Unprescaled luminosity for plots
constexpr double _jp_lumi = 35.918235; // 36.810440678; (outdated)
//} END Trigger stuff

//{ BEGIN RUN2 params
constexpr double _jp_sqrts = 13000.; // GeV
constexpr double _jp_emax = _jp_sqrts/2.; // Max possible jet pt
constexpr double _jp_recopt = 15; // Min observed jet pt
//} END Run2

//{ BEGIN Outdated, probably shouldn't be touched
// Eta phi exclusion due to ecal problems
constexpr bool _jp_etaphiexcl = false;
const constexpr char _jp_etaphitype[] = "h2hotr"; // h2hotr (Robert) or h2hotm (Mikko)

// Veto jets near ECAL boundaries in JetID
constexpr bool _jp_doECALveto = false;//true;
const constexpr char _jp_ecalveto[] = "lumicalc/ecalveto.root";
//} END Outdated

/////////////////////////////////////////////////////
// Additional analysis switches, might be outdated //
/////////////////////////////////////////////////////

// Correct for trigger efficiency based on MC
constexpr bool _jp_dotrigeff = false;
 // Correct pT<114 GeV only, if above _jp_dotrigeff=true
constexpr bool _jp_dotrigefflowptonly = false;
// Correct for time-dependence (prescales) in data
constexpr bool _jp_dotimedep = false;
// For creating smearing matrix
constexpr bool _jp_doMatrix = false;

// Check that run / lumi section was listed in the .csv file
constexpr bool _jp_dolumcheck = false;
// Veto bad run list
constexpr bool _jp_dorunveto = false;
// Check for duplicates (warning: takes a lot of memory!)
constexpr bool _jp_checkduplicates = false;
// Only load selected branches (large speedup, but be careful!)
constexpr bool _jp_quick = true;
// Center uncertainties around ansatz (true) or data (false)
constexpr bool _jp_centerOnAnsatz = false;
constexpr bool _jp_centerOnTheory = false;

// Plot Pythia for final PRL results
constexpr bool _plotPythia = false;

// TODO: Do we still need these?
// Minimum and maximum pT range to be plotted and fitted
constexpr double _jp_fitptmin = 114;//43;
// Changed on 2013-05-020: analysis from 49 GeV to 56 GeV
constexpr double _jp_xmin57 = 114;//56;
constexpr double _jp_xminpas = 114;//56;
constexpr double _jp_xmin = 114;//24.;//20.;
constexpr double _jp_xmax = 1497;//TEMP PATCH for partial data //1999.;

constexpr double _jp_xsecMinBias = 7.126E+10;

// Draw againts HERAPDF1.7 instead of PDF4LHC (drawSummary)
constexpr bool _jp_herapdf = false;

// Produce plots (drawSummary and drawRunHistos)
constexpr bool _jp_pdf = true;

#endif // __settings_h__
