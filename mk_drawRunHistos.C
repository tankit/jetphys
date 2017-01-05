{
  //gROOT->ProcessLine(".L Config.cpp+");
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L drawRunHistos.C+g");

  drawRunHistos("DATA");
  //drawRunComposition("DATA");
}
