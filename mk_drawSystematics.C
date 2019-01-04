// Purpose: Draw inclusive jet cross section systematics
// Author:  mikko.voutilainen@cern.ch
// Created: March 23, 2010
// Updated: March 23, 2010
{

  // compile code
  //gROOT->ProcessLine(".x cms_jes_style.C");
  gROOT->ProcessLine(".L tdrstyle_mod.C");
  setTDRStyle();
  gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L drawSystematics.C+");

  gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  cafe::Config *cfg = new cafe::Config("pfjet");
  string type = cfg->get("type","DATA");
  delete cfg;

  drawSystematics(type); // 36X JEC
  //drawSystematics(type, true); // 38X JEC
}
