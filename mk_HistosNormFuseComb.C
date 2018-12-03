// Purpose: Normalize and combine jet physics analysis histograms
// Author:  hannu.siikonen@cern.ch
// Created: December 3, 2018
{
  // Normalization: compile code
  gROOT->ProcessLine(".L HistosNormalize.C+g");

  HistosNormalize();

  // Fuse: compile code
  gROOT->ProcessLine(".L HistosFuse.C+g");

  HistosFuse();

  // Combination: compile code
  gROOT->ProcessLine(".L HistosCombine.C+g");

  HistosCombine();
}
