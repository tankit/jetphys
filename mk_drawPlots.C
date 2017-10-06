// Purpose: Draw plots of analysis stages
// Author:  mikko.voutilainen@cern.ch
// Created: March 21, 2010
// Updated: March 21, 2010
{

  // compile code
  //gROOT->ProcessLine(".x cms_jes_style.C");
  gROOT->ProcessLine(".L tdrstyle_mod12.C");
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L drawPlots.C+");

  gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  cafe::Config *cfg = new cafe::Config("pfjet");
  string type = cfg->get("type","DATA");
  delete cfg;
  
  
  /*      
  drawJetID(type);
  drawTriggerRatio(type, "PF");
  if (type=="DATA") drawTriggerRatio(type, "CALO");
  drawUnfolding(type, "PF");
  if (type=="DATA") drawUnfolding(type, "CALO");
  drawResolution(type);
  drawMultijets(type);
  
  drawBFraction(type);
  drawBTagPurity(type);
  drawBTagEfficiency(type);
  drawBTagNorm(type);

  drawMCPurity();
  drawMCEfficiency();
  drawTemplatePurity();
  drawTemplate();
  */

  //drawUnfolding(type, "PF");
  //drawUnfolding(type, "CALO");
  //drawJetID(type);

  // For PAS (or check triggers)
  drawTriggerRatio(type, "PF");
  //drawTriggerRatio(type, "CALO");
  
  //drawMultijets(type);
  //drawMCPurity();
  //drawMCEfficiency();
  //drawTemplatePurity();
  //drawTemplate();
}
