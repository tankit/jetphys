// Purpose: Jet physics analysis chain based on ROOT6 (6.02/10)
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010
// Updated: June 1, 2015
{

  // ***** BEFORE RUNNING THE ANALYSIS ****** //
  // NB: more detailed instructions are at the end
  // * Provide output JSON file or get the golden JSON from lxplus:
  //      /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/
  //   => settings.h : jp::json
  //
  // * Produce brilcalc_lumibyls.csv on lxplus:
  //   [http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html]
  //       cd public
  //       setenv PATH $HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
  //       brilcalc lumi -b "STABLE BEAMS" --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i [JSON] -u /fb
  //       brilcalc lumi -i [JSON] -b "STABLE BEAMS" -o brilcalc_lumibyls.csv --byls --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json  --minBiasXsec 69000
  //   [could add --hltpath to calculate luminosity per path]
  //   => settings.h : jp::lumifile
  //
  // * Produce pileup reweighing files for data and MC
  //   MC: ProcessedTree->Draw("EvtHdr_.TrPu>>pileupmc(500,0,50)");
  //   Data: use brilcalc_lumibyls.csv [need code for this]
  //   => settings.h : jp::pudata , jp::pumc
  //
  // * Run duplicate checker with new data files, but turn of later (slow)
  //   => settings.h : jp::checkduplicates
  //
  // * Point mk_HistosFill.C to correct data/MC files,
  //   mk_theory.C to correct theory, and HistosFill.C to correct JSON etc.
  //   [tbd: make these all configurable]
  // * Optimize bins for given luminosity (optimizeBins.C, HistosBasic.C)
  //
  // (Check that new runs are listed in JetMETTau stream)
  // (Update trigger thresholds in mk_HistosCombine.C)
  // (Update lumi and max pT in settings.h)
  // Reoptimize bins for given luminosity (optimizeBins.C, HistosBasic.C)
  // Instructions of getting the external data are at the end

  // ***** AFTER RUNNING THE ANALYSIS ****** //
  // (=> update his part)
  // Check the new results in the Analysis Note:
  // cd pas; cp -pi ../comparisons/pasPF.root .; root -l -b -q mk_pasPlots.C
  // cp -p fig7_Comparison_*_LUM.pdf ../latex/figs/
  // cp -p pflow_*_LUM.pdf ../latex/figs/
  // cd ../latex; sh copy_anplots.sh; pdflatex inclbpfjets.tex


  // Setup directories
  gROOT->ProcessLine(".! mkdir pdf");

  // Load settings
  #include "settings.h"


// Step 1:  - make histograms of measured jet properties (pt, dijet mass)
//          - make histograms of correction factors (JEC, JER)
//         (- for MC include generator truth information for analysis closure)
  cout << "\nStep 1: Histogram measured jet pt and correction factors"
       << "\n========================================================\n";
  gROOT->ProcessLine(".x mk_HistosFill.C");
  //cout << "NB: this step is now performed in JetAnalyzer" << endl;
  //cout << "use the output of 'analyzeJets Jets_FWL_cfg.py'" << endl;
  cout << "NB: this step is very slow (~4h) and is switched off" << endl;
  cout << "NB: to run, execute 'root -l -b -q mk_HistosFill.C'" << endl;

// Step 2a: - apply corrections, normalize luminosity and eta width
  cout << "\nStep 2a: Apply corrections and normalization factors"
       << "\n====================================================\n";
  gROOT->ProcessLine(".x mk_HistosNormalize.C+g");

// Step 2b: - stitch different triggers together
  cout << "\nStep 2b: Stitch different triggers together"
       << "\n===========================================\n";
  gROOT->ProcessLine(".x mk_HistosCombine.C+g");

// Step 3a: - unfold spectrum iteratively using Ansatz method
// Step 3b: - unfold spectrum using reweighted Pythia MC
  // need to comment out because causes crash when opening 2b.root again
  //cout << "\nStep 3a: Unfold spectrum using Ansatz method"
  //   << "\n============================================\n";
  //gROOT->ProcessLine(".x mk_unfold.C");
  gROOT->ProcessLine(".x mk_dagostini.C");

// Step 2c:  - reformat theory curves
  cout << "\nStep 2c: Reformat theory predictions"
       << "\n======================================\n";
  gROOT->ProcessLine(".x mk_theory.C+g");

// Step 3a: - unfold spectrum using NLO forward smearing
  //cout << "\nStep 3a: Unfold spectrum using NLO forward smearing"
  //   << "\n===================================================\n";
  //gROOT->ProcessLine(".x mk_forwardsmear.C");

// Step 3: - unfold spectrum using d'Agostini method
  cout << "\nStep 3: Unfold spectrum using d'Agostini method"
       << "\n===================================================\n";
  gROOT->ProcessLine(".x mk_dagostini.C+g");

// JEC/JER residuals: move by hand to systematics.C
//gROOT->ProcessLine(".x mk_resolution.C");

// Step 4:  - produce systematics:
//            (i) vary JEC
//            (ii) vary JER
//            (iii) vary JID
//            (iv) vary lumi
//            for all of these, use parameterized sources
  cout << "\nStep 4: Calculate systematics"
       << "\n============================================\n";
  gROOT->ProcessLine(".x mk_systematics.C+g");

// Step 4b: - calculate and plot rapidity bias (+reco efficiency)
  cout << "\nStep 4b: Calculate (and plot) rapidity bias"
       << "\n============================================\n";
  gROOT->ProcessLine(".x mk_rapiditybias.C+g");

// Sted 4c: - apply ad-hoc corrections for various biases
//            (i) rapidity bias
//            (ii) JEC residual (optional)
//            (iii) JER residual (not yet implemented)
  cout << "\nStep 4c: Apply ad-hoc bias correctios (rapidity, JEC)"
       << "\n=====================================================\n";
  gROOT->ProcessLine(".x mk_correctHistos.C+g");

// Step 6:  - produce pretty plots of analysis steps
  cout << "\nStep 6: Draw plots (raw spectra, triggers, unfolding)"
       << "\n=====================================================\n";
  gROOT->ProcessLine(".x mk_drawPlots.C");

// Step 7: - plot systematics
  cout << "\nStep 7: Draw plots (systematics)"
       << "\n================================\n";
  gROOT->ProcessLine(".x mk_drawSystematics.C+g");

// Step 8: - plot summary
  cout << "\nStep 8: Draw plots (summary)"
       << "\n================================\n";
  gROOT->ProcessLine(".x mk_drawSummary.C+g");

// Step 9: - plot summary
  cout << "\nStep 9a: Draw plots (comparison)"
       << "\n================================\n";
  gROOT->ProcessLine(".x mk_drawComparison.C");
  cout << "...skip obsolete code, is now run in drawSummary.C" << endl;

  if (jp::isdt) {
    cout << "\nStep 9b: Draw run quality checks"
	 << "\n================================\n";
    gROOT->ProcessLine(".x mk_drawRunHistos.C");
    cout << "\nStep 9c: Draw JEC checks"
	 << "\n================================\n";
    gROOT->ProcessLine(".x mk_jecChecks.C+g");
    gROOT->ProcessLine(".x mk_drawJEC.C");
  }


// Step 10: - produce a note draft
  //cout << "\nStep 10: Produce a draft note on latex"
  //     << "\n====================================\n";
  //gROOT->ProcessLine(".x mk_latex.C");


  // Getting the list of runs in the root file
  // ------------------------------------------
  // - Open the file and the ProcessedTree tree interactively
  // - command:
  //   ProcessedTree->Draw("EvtHdr_.mRun>>h(20000,160000,180000)")
  //   for (int i = 1; i != h->GetNbinsX()+1; ++i) { if(h->GetBinContent(i)!=0)
  //   cout << h->GetBinLowEdge(i) << ", "; } cout << endl;

  // Latest JSON file should be here:
  // --------------------------------
  // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions11/7TeV/
  // NB: but ideally will get the output JSON from Kostas!

  // Combine JSON files, if many:
  // ----------------------------
  // (- cut out the overlapping regions by hand (better methods?)
  // (compareJSON.py --or Cert_PromptReco_JSON.txt Cert_ReReco_JSON.txt Cert_PromptReco_OR_ReReco_JSON.txt)
  // - mergeJSON.py JSON1.json JSON2.json JSON3.json ... --output=JSONS.json

  // Get pixellumi_by_LS.csv file with pixelLumiCalc.py (new, Feb22)
  // --------------------------------------------------
  // https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/122.html
  // cmsrel CMSSW_5_0_1
  // cmsenv
  // cvs co -r HEAD RecoLuminosity/LumiDB
  // cd RecoLuminosity/LumiDB/; scram b
  // setenv FRONTIER_FORCELOAD long
  // cd lumicalc
  // - ../RecoLuminosity/LumiDB/scripts/pixelLumiCalc.py -i lumiSummary_Jan16th.json overview [total recorded for double-checking below: 4.943/fb]
  // - ../RecoLuminosity/LumiDB/scripts/pixelLumiCalc.py -i lumiSummary_Jan16th.json lumibyls -o pixellumi_by_LS.csv
  // - $HOME> rsync -rut lxplus.cern.ch:~/scratch0/CMSSW_4_2_8/src/lumicalc .

  // Create pile-up profiles for data, per trigger
  // ---------------------------------------------
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
  // - CMSSW_5_0_1 HEAD
  // Repeat for all triggers:
  // - pixelLumiCalc.py lumibyls -i lumiSummary_Jan16th.json --hltpath "HLT_Jet60_v*" -o pixellumi_by_LS_jt60.csv
  // - rfcp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt
  // (should use newer one?: pileup_2011_JSON_pixelLumi.txt)
  // use own script to merge pixel-lumibyls, pixel-lumibyls-bytrigger and
  // the HF luminosity in the official pileupJSON into pileupJSON_jt*.txt:
  // root -l -b -q mk_calcPileup.C (calcPileup.C::calcPileupJSON, ::mergeJSON)
  // (new tools to replace above, not tested:
  //   pileupReCalc_HLTpaths.py -i YourOutput.csv --inputLumiJSON pileup_2011_JSON_pixelLumi.txt -o My_HLT_corrected_PileupJSON.txt
  //   pileupCalc.py -i lumiSummary_Jan16th.json --inputLumiJSON pileupJSON_jt60.txt --calcMode true --minBiasXsec 73500 --maxPileupBin 50 --numPileupBins 50  pileup_jt30.root)
  // - MC profiles should be created by hand and are found at
  //   qcdjet/pileup/pileup_PY.root


  // ********* all sorts of obsolete commands ************

  // Create pile-up profiles for data, per trigger (obsolete, Feb22)
  // ---------------------------------------------
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
  // (Need HF lumi reweighed by pixel lumi to get good results per crossing,
  //  but this is not yet available on Feb22; Use antioniv's script instead)
  // Calculate eff lumi per LS for each trigger, then use single Xing file
  // - trying this in CMSSW_5_0_1 HEAD (not sure if right)
  // - [lumiCalc2.py lumibyls -b stable -i lumiSummary_Jan16th.json --norm pp7TeV --hltpath "HLT_Jet60_v*" -o lumicalc2v1_by_LS_jt60.csv [recorded 1.099/pb]]
  // Repeat for all triggers:
  // - [pixelLumiCalc.py lumibyls -i lumiSummary_Jan16th.json --hltpath "HLT_Jet60_v*" -o pixellumi_by_LS_jt60.csv]
  // - lumiCalc2.py lumibyls -i lumiSummary_Jan16th.json --hltpath "HLT_Jet60_v*" -o lumicalc2v1_by_LS_jt60.csv
  // - [lumiCalc2.py lumibylsXing --xingMinLum 0.1 -b stable -i lumiSummary_Jan16th.json --norm pp7TeV -o lumicalc2v1_by_LSXing.csv]
  // - rfcp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt
  // - [cp ~antoniov/public/analysis/ExclusiveDijets/computeInstLumiByXing.py .]

  // Per-trigger pile-up profile (obsolete)
  // --------------------------------------
  // Need profiles per run, was provided by Andre David previously
  // - ../RecoLuminosity/LumiDB/scripts/estimatePileup2.py -i Nov13th/lumiSummary.json --xingMinLum=0.0100 --minBiasXsec=71300 --maxPileupBin=60 --saveRuns pileup47fb.root
  // - mk_calcPileup.C (modify for all runs in a single file)

  // Estimate pile-up profile
  // ------------------------
  // (https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities)
  // - cd CMSSW_4_2_8/src/lumicalc
  // - ../RecoLuminosity/LumiDB/scripts/lumiCalc2.py -i Cert_*_JSON.txt -o Cert_*_JSON.csv -xingMinLum 5.0e-03 lumibylsXing
  // - ../RecoLuminosity/LumiDB/scripts/estimatePileup2.py --csvInput=Cert_*_JSON.csv pudist.root
  // - $HOME> rsync -rut lxplus.cern.ch:~/scratch0/CMSSW_4_2_8/src/lumicalc .

  // Get lumi_by_LS.csv file with lumiCalc2.py (replaced by pixelLumiCalc)
  // ----------------------------------------
  // (old: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc)
  // (old: https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc2011)
  // (new: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc)
  // - scramv1 p CMSSW CMSSW_4_2_8
  // - cd CMSSW_4_2_8/src
  // - addkpg RecoLuminosity/LumiDB; scramv1 b
  // - [cvs update -r V03-00-06 (Aug10->Sep2, what's in HEAD?); this is v1]
  // - [cvs update -r V03-03-00 (Oct 5, bug fix for dL/dt correction)]
  // - [cvs update -r V03-07-00 (Nov 15, bug fix for -hltpath)]
  // - cvs update -r HEAD (for --correctionv2 switch)
  // - scramv1 b
  // - cmsenv
  // - [$HOME> rsync -rut lumicalc lxplus.cern.ch:~/scratch0/CMSSW_4_2_8/src/]
  // - cd lumicalc
  // - [../RecoLuminosity/LumiDB/scripts/lumiCalc2.py -i Cert_*_JSON.txt lumibyls -o lumicalc2_by_LS.csv]
  // - [../RecoLuminosity/LumiDB/scripts/lumiCalc2.py -norm pp7TeV -i Cert_*_JSON.txt lumibyls -o lumicalc2_by_LS.txt]
  // - ../RecoLuminosity/LumiDB/scripts/lumiCalc2.py --correctionv2 -norm pp7TeV -i lumiSummary_Jan16th.json lumibyls -o lumicalc2v2_by_LS.txt
  // - $HOME> rsync -rut lxplus.cern.ch:~/scratch0/CMSSW_4_2_8/src/lumicalc .

  // obsolete:
  // - lumiCalc.py -i Cert_*_JSON.txt lumibyls -o lumicalc_by_LS.csv
  // - $HOME> rsync -rut lxplus.cern.ch:~/scratch0/CMSSW_4_2_2/src/lumicalc .

  // obsolete (prescales are now stored in the event record):
  // Get prescales with lumiCalc.py
  // ------------------------------
  // - cd lumicalc (as above)
  // - check that lumicalc_by_LS.csv is available for reading in the run list
  // - root -l -b -q mk_prescales.C
  // - root -l -b -q merge_prescales.C++
  // - $HOME> rsync -rut lxplus.cern.ch:~/scratch0/CMSSW_4_2_2/src/lumicalc .

  // Get prescales with lumiCalc.py (July 9, OBSOLETE!! New method above)
  // (- scramv1 project -n CMSSW_3_8_0_pre7 CMSSW CMSSW_3_8_0_pre7
  // (- cd CMSSW_3_8_0_pre7; cmsenv;
  // - setup CMSSW as above
  // - rsync -rut maccms316.cern.ch:/Users/voutila/pfjet/prescales .
  // - cd prescales
  // - sh db_prescales.sh (edit to add list of runs, and to set L1 prescales)
  //
  // Above script is based on these manual commmands:
  // edmConfigFromDB --runNumber 138924 --format summary.ascii --paths HLT_L1_BscMinBiasOR_BptxPlusORMinus
  // Except that the dynamic prescales just make a mess of this...
  // But, for now we are assuming the first number on the list is correct (4E29)
  // edmConfigFromDB --runNumber 138924 --format summary.ascii --paths HLT_Jet15U | grep HLT_Jet15U | tail -n1 | awk '{print $4}'


  // Getting the prescale files (OBSOLETE!! New method was above):
  // - copy prescales directory to lxplus
  // - scramv1 project -n CMSSW_3_3_4 CMSSW CMSSW_3_3_4
  // - cd CMSSW_3_3_4; cmsenv; cd ../prescales
  // - get run range for [132440,xxxxxx] from (or use dates):
  //   https://twiki.cern.ch/twiki/bin/view/CMS/LumiWiki2010Data
  // - HLT keys from https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/RunSummary
  //   or from (much easier, but not official)
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/Data2010HLTMenus
  // - run "python prescales.py <HLT key>" by adding entries to prescales.sh
  // - sh prescales.sh
  // - add hltfiles*.txt to hltfiles.txt
  // - add triggers to mk_prescales.C if needed
  // - root -l -b -q mk_prescales.C



  // Step eXtras: only to be run every now and then, and best by hand
  // gROOT->ProcessLine(".x mk_drawTagTables.C");
  // gROOT->ProcessLine(".x mk_optimizeBins.C");

  // Notes:
  //
  // Triggers
  //
  // The triggers may become a problem for PFJets. I checked the
  // efficiency of hlt_l1jet6u, hlt_jet15u etc. relative to MinBias
  // after tight ID at |eta|<2.0, and the efficiency barely rose up to
  // 50% in the last bin at pT~40 GeV. We would either need a looser
  // L1+HLT trigger, or significantly more MinBias statistics to cover
  // the pT gap between 30 and about 50-70 GeV.
  //
  // JEC
  //
  // The jets have been corrected when pT(uncorr)>4 GeV. For PFJets this
  // means a threshold of around 5 GeV. Below this the spectrum is not
  // reliable. I already observed this when trying to fit ansatz to the
  // smeared spectrum: fit from 5-50 GeV is good with chi2/NDF~1, while fit
  // from 4-50 GeV does not converge properly anymore. By eye the spectrum
  // at pT<5 GeV doesn't look that bad, but the log scale hides a lot.
  //
  // The measurement down to 5 GeV is probably needed for reliable unfolding
  // of the results at pT>10 GeV. Should make a plot of Gen pT in a bin of
  // pT reco at [10,11] GeV to see how low we need to measure data.
  //
  //
  // JetID
  //
  // Need to prove efficiency (~100%?) and rejection power (70-100%?) of
  // tight PFJetID. Make plots of pT spectrum and MET/sumET etc. before
  // and after ID to estimate possible contamination:
  // before = S+B, after = eff*S+rej*B. If B>>S, then purity after cuts
  // may not be 100%. The B is not necessarily "pure" background, it can
  // also be very badly measured jets, which are not well modeled by MC
  // and get boosted up in an expontetially falling pT spectrum. We still
  // want to reject these so that we can trust MC-based JEC and JER for
  // the rest of the jets (proven by comparing RD and MC distributions).
}
