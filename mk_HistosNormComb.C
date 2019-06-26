// Purpose: Normalize and combine jet physics analysis histograms
// Author:  hannu.siikonen@cern.ch
// Created: December 3, 2018
{
  // Normkalization: compile code
  gROOT->ProcessLine(".L HistosNormalize.C+g");

  HistosNormalize();

  // Combination: compile code
  gROOT->ProcessLine(".L HistosCombine.C+g");

  HistosCombine();
}
