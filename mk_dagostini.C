// Purpose: d'Agostini ("Baysian" or Richardson-Lucy) unfolding, including
//          response matrix generation from NLO theory and parameterized JER
// Author:  mikko.voutilainen@cern.ch
// Created: September 2, 2012
// Updated: June 4, 2015
{

  // compile code
  gROOT->ProcessLine(".L tools.C+");

  // Retrieve RooUnfold package from
  //   http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
  //   svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
  //   cd RooUnfold; make; cd ..
  gSystem->Load("RooUnfold/libRooUnfold"); // .so

  gROOT->ProcessLine(".L dagostini.C+");

  #include "settings.h"

  dagostiniUnfold(_jp_type);
  drawDagostini(_jp_type);
}
