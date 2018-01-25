// Purpose: Analysis settings for jet physics package
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: June 1, 2015

// All the settings variables are in global name space
// To avoid conflicts, and to make them easier to find
// global variables all start with _jp_ (for JetPhysics)

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


constexpr bool strings_equal(char const * a, char const * b) {
  return *a == *b and (*a == '\0' or strings_equal(a + 1, b + 1));
}


// Debugging info
constexpr bool _jp_debug = false;
// Will assert be used in the programs?
//#define USEASSERT // Expert option: comment this if certain that no problems occur. More than 1/3 off the run time.

//{ BEGIN fundamental file settings
// Do we use CHS jets? ("CHS" for yes, "", for legacy/no)
const constexpr char _jp_chs[] = "CHS";
// Algorithm to use ("AK4PF" or "AK8PF" + "chs" for chs jets)
const constexpr char _jp_algo[] = "AK4PFchs";
// Data type ("DATA", "MC", or "HW")
const constexpr char _jp_type[] = "MC";
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
-1; // all
//10; // debug
//10000; // shorter test run
//100000; // short test run
//1000000; // shortish test run
// Number of events to skip from the beginning (for debugging)
constexpr Long64_t _jp_nskip = 0;
// Step between events
constexpr Long64_t _jp_skim = 0; // "prescale", off if zero
// Only load selected branches (large speedup, but be careful!)
constexpr bool _jp_quick = true;
// Save the infos at times. Setting this to true might help with long runs. (Not recommended
constexpr bool _jp_save = false;
//} END run settings


//{ BEGIN A group of constexpr switches and their settings to turn some parts of the code on/off
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

// QGL studies by Ozlem (use an external file for qgl studies)
constexpr bool _jp_doqglfile = false;
const constexpr char _jp_qglfile[] = "output-DATA_RunGfullpart-1.root";

// At some point of time Kostas stored UNCORRECTED four-vector. Current status: CORRECTED
// NOTE: this is a source of constant anxiety, should be rechecked from time to time
constexpr bool _jp_undojes = true;
// We can choose also not to apply the new jes onto a four-vector
constexpr bool _jp_redojes = true;
// For debugging
constexpr bool _jp_skipl2l3res = false;

// Reapply json selection based on the latest one (check lumicalc if false!).
constexpr bool _jp_dojson = true;
const constexpr char _jp_json[] = "lumicalc/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
//const constexpr char _jp_json[] = "lumicalc/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt";

// Calculate luminosity on the fly based on .csv file and take only events with non-zero luminosity.
constexpr bool _jp_dolumi = true;
const constexpr char _jp_lumifile[] = "lumicalc/brilcalc_lumibyls16.csv";
//const constexpr char _jp_lumifile[] = "lumicalc/brilcalc_lumibyls17.csv";

// Add prescale information from a file
constexpr bool _jp_doprescale = false; // CAUTION, this shouldn't be nowadays necessary
const constexpr char _jp_prescalefile[] = "pileup/prescales74x.txt";

constexpr bool _jp_dotrpufile = false; // CAUTION, the histo file in the repo is 2015
const constexpr char _jp_trpufile[] = "pileup/MUperLSvsRUN_MB.root";

// Decide whether or not to simulate triggers from MC (this is slow)
constexpr bool _jp_domctrigsim = true;
// Use "mc" trigger for whole pT range instead of stiching triggers together in histosCombine (requires trigsim)
constexpr bool _jp_usemctrig = true;
// This is the 13 TeV 25 ns list of triggers (Run2016BCDEFG)
// Check the recommended settings from https://twiki.cern.ch/twiki/bin/view/CMS/InclusiveJetsLegacy
constexpr int _jp_ntrigs = 9; // CAUTION: In 2016 jt450 is unprescaled but in 2017 it is not (jt500!)
const constexpr char* _jp_triggers[_jp_ntrigs] =
{"jt40",    "jt60",    "jt80",   "jt140",   "jt200",   "jt260",   "jt320",   "jt400",  "jt450"}; // "jt500"};
// reference trigger for PU profile in the mc folder and for trigger lumi weighting
const constexpr char _jp_reftrig[] = "jt450"; // "jt500";
// Thresholds for the corresponding triggers (same as in trigger name)
const constexpr double _jp_trigthr[_jp_ntrigs] =
{40,        60,        80,       140,       200,       260,       320,       400,      450}; // 500};
// Trigger ranges (differ from thresholds)
const constexpr double _jp_trigranges[_jp_ntrigs][2] =
{ {0,84}, {84,114}, {114,196}, {196,272}, {272,330}, {330,395}, {395,468}, {468,548}, {548,6500} }; // V[5,6], AK4

// Trigger lumis for weighting
constexpr bool _jp_usetriglumi = true; // use luminosity numbers below, in /ub
const constexpr double _jp_triglumi[_jp_ntrigs] = // in /ub
// 2016 final: brilcalc lumi -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --byls --normtag=$PWD/../../../Normtags/normtag_DATACERT.json --minBiasXsec 80000 --hltpath="HLT_PFJet*" -o brilcalc_lumibyls16.csv
{264821.835,718829.826,2733420.74,23966019.286,102854094.409,587728815.19,1753996573.885,5141160409.826,35860066023.522};
// 2017 final: brilcalc lumi -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt --byls --normtag=/eos/user/h/hsiikone/Normtags/normtag_BRIL.json --minBiasXsec 80000 --hltpath="HLT_PFJet*" -o brilcalc_lumibyls17.csv
// {480907.381,1135408.863,4300901.657,39784735.913,218494113.161,553531633.649,1402704877.698,4201460774.775,10434425834.42,41366526620.215};
// Unprescaled luminosity for plots
constexpr double _jp_lumi = _jp_triglumi[_jp_ntrigs-1]/1000000000.0; // in /fb
//} END Switches and their settings

//{ JEC and IOV settings. In the modern world we have a group of IOV's for which different corrections are applied.
// https://github.com/cms-jet/JECDatabase/tree/master/tarballs
// Summer16_03Feb2017G_V7_DATA.tar.gz [BCD, EF, G, H]
const constexpr char _jp_jecgt[] = "Summer16_07Aug2017";//"Summer16_03Feb2017";//BCD_";//"Summer15_50ns";// "Summer16_23Sep2016";
const constexpr char _jp_jecvers[] = "_V4";//"_V7";//"V4"; // Summer16_03Feb // "V6"; // Summer16_23Sep // "V2" ; // Spring16

// Use Intervals-Of-Validity for JEC
constexpr const bool _jp_useIOV = true ;//false
constexpr const unsigned int _jp_nIOV = 3; //4;
const constexpr char* _jp_IOVnames[_jp_nIOV] =
  {"BCD",    "EF",    "GH"};//{"BCD",    "EF",    "G",   "H"};
// Trigger IOVs: the 1 for -inf and 400000 for inf (currently)
const constexpr int _jp_IOVranges[_jp_nIOV][2] =
  { {1,276811}, {276831,278801}, {278802,400000} }; // Spring/Summer16_23Sep2016
//  { {1,276811}, {276831,278801}, {278802,280385}, {280919,400000} }; // Spring/Summer16_23Sep2016
//} END JES and JEC


//{ BEGIN For MC reweighting purposes, PU profiles for data and MC
constexpr bool _jp_reweighPU = true;
const constexpr char _jp_pudata[] = "pileup/pileup_DT.root";
const constexpr char _jp_pumc[]   = "pileup/pu.root";
//} END PU profiles


//{ BEGIN MC: Process pThatbins instead of flat sample
// All the following values need to be updated, if other source tuples are used.
// These can be obtained using getSliceEvts.C
constexpr bool _jp_pthatbins = true;//false;
// Number of pthat bins
constexpr unsigned int _jp_npthatbins = 14;
// The corresponding ranges, the last number is ~inf
const vector<double> _jp_pthatranges =
  {30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200,20000};
// The corresponding lumis in an arbitrary scale
const vector<double> _jp_pthatsigmas =
  {140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
  //{140932000,19204300,2762530,471100,117276,7823,648.2,186.9,32.293,9.4183,0.84265,0.114943,0.00682981,0.000165445};
const vector<double> _jp_pthatnevts =
  {9699558,9948791,7742665,5748730,7838066,11701816,3959986,9628335,11915305,6992746,2477018,1584378,596904,391735};
// The filenames need to be given here and in mk_histosFill, since ROOT is exceedingly stupid
const vector<string> _jp_pthatfiles = {
    "QCD_Pt_30to50_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_50to80_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_80to120_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_120to170_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_170to300_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_300to470_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_470to600_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_600to800_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_800to1000_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_1000to1400_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_1400to1800_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_1800to2400_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_2400to3200_TuneCUETP8M_13TeV_pythia8.root",
    "QCD_Pt_3200toInf_TuneCUETP8M_13TeV_pythia8.root",
};
//} END pthatbins stuff


//{ BEGIN RUN2 params
constexpr double _jp_sqrts = 13000.; // GeV
constexpr double _jp_emax = _jp_sqrts/2.; // Max possible jet pt
constexpr double _jp_recopt = 15; // Min observed jet pt
constexpr double _jp_xsecMinBias = 7.126E+10;
//} END Run2

//{ BEGIN Special histosFill analyses that are not typically needed
// Eta phi exclusion due to ecal problems
constexpr bool _jp_doetaphiexcl = false;
const constexpr char _jp_etaphitype[] = "h2hotr"; // h2hotr (Robert) or h2hotm (Mikko)

// Veto jets near ECAL boundaries in JetID
constexpr bool _jp_doECALveto = false;//true;
const constexpr char _jp_ecalveto[] = "lumicalc/ecalveto.root";

// Check for duplicates (warning: takes a lot of memory!)
constexpr bool _jp_checkduplicates = false;
//} END Special analyses that are not typically needed

//{ BEGIN histosNormalize (leave these off if not interested on details)
// Correct for trigger efficiency based on MC
constexpr bool _jp_dotrigeff = false; // CAUTION: Requires output-MC-1.root file
// Correct pT<114 GeV only, if above _jp_dotrigeff=true
constexpr bool _jp_dotrigefflowptonly = false; // CAUTION: Not needed
// Correct for time-dependence (prescales) in data
constexpr bool _jp_dotimedep = false; // CAUTION: Unknown consequences
//} END histosNormalize

//{ BEGIN drawSummary
// Center uncertainties around ansatz (true) or data (false)
constexpr bool _jp_centerOnAnsatz = false;
constexpr bool _jp_centerOnTheory = false;
// Plot Pythia for final PRL results
constexpr bool _jp_plotPythia = false;
// Draw againts HERAPDF1.7 instead of PDF4LHC (drawSummary)
constexpr bool _jp_herapdf = false;
// Produce plots (also drawRunHistos)
constexpr bool _jp_pdf = true;
//} END drawSummary

//{ BEGIN limits for histosNormalize, dagostini, drawSummary, drawRunHistos
// Minimum and maximum pT range to be plotted and fitted
constexpr double _jp_fitptmin = 114;//43;
// Changed on 2013-05-020: analysis from 49 GeV to 56 GeV
constexpr double _jp_xmin57 = 114;//56;
constexpr double _jp_xminpas = 114;//56;
constexpr double _jp_xmin = 114;//24.;//20.;
constexpr double _jp_xmax = 1497;//TEMP PATCH for partial data //1999.;
//} END limits

#endif // __settings_h__
