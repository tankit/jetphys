{

  //gROOT->ProcessLine(".L tdrstyle_mod12.C");
  //setTDRStyle();

  //gStyle->SetOptStat(0);
  //gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L drawSummary.C+");

  //gSystem->Setenv("CAFE_CONFIG", "pfjet.config");
  //cafe::Config *cfg = new cafe::Config("pfjet");
  //string type = cfg->get("type","DATA");
  //delete cfg;

  //drawDataTheoryRatio3x2(false, false); // AK5 36X

  #include "settings.h"

    drawEtaSpectra(jp::type);
  //drawDataTheoryRatio();
  //drawDataTheoryRatioLog(); // graphs for pasPlots
  //drawDataTheoryRatio3x2(false, false); // AK5 36X
  //drawDataTheoryRatio3x2(true, false); // AK7 38X
  //drawDataTheoryRatio3x2(false, true); // AK5 38X
  //drawDataTheoryRatio3x2(true, true); // AK7 38X
    //drawDataTheoryRatio3x2();//_algo=="AK7" , false);
  //drawUnfoldSummary(type, "PF");
  //if (type=="DATA") drawUnfoldSummary(type, "CALO");

  //drawJetIDSummary(type, "JetID");
  //drawJetIDSummary(type, "MetCut");
  //drawJetIDSummary(type, "JetIDPlusMetCut");
  //drawNPCorr(type);

    /*
  if (jp::ismc)
    drawClosure(jp::type,false); // Incl jet (PF in 36x)
    */
  //drawClosure("HW",false); // Incl jet (PF in 36x)
  //drawClosure("MC",false,"PFAK7"); // Incl jet
  //drawClosure("MC",false,"JPT"); // Incl jet
  //drawClosure("MC",false,"CALO"); // Incl jet
  //drawClosure("MC",true); // b-jet
  //drawClosure("MC",false,"PF",true); // Incljet (PF in 38X)
  //drawClosure("MC",false,"PFAK7",true); // Incl jet

  //drawTheoryDataRatio(); // AK5, 36X JEC
  //drawTheoryDataRatio(false, true); // AK5, 38X JEC
  //drawTheoryDataRatio(true, true); // AK7, 38X JEC

  //drawDataTheoryRatio3x2(true); // AK7
}
