{

  // For JEC (for uncertainty)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L drawJetPt.C+g");

  drawJetPt("0.0-0.5");
  drawJetPt("0.5-1.0");
  drawJetPt("1.0-1.5");
  drawJetPt("1.5-2.0");
  drawJetPt("2.0-2.5");
  drawJetPt("2.5-3.0");
  drawJetPt("3.2-4.7");
}
