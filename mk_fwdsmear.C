// Purpose: d'Agostini ("Baysian" or Richardson-Lucy) unfolding, including
//          response matrix generation from NLO theory and parameterized JER
// Author:  mikko.voutilainen@cern.ch
// Created: September 2, 2012
// Updated: June 4, 2015


  #include "tools.h"
   R__LOAD_LIBRARY(tools.C+)

  // Retrieve RooUnfold package from
  //   http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
  //   svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
  //   cd RooUnfold; make; cd ..

  #include "RooUnfold/src/RooUnfold.h"
  #include "RooUnfold/src/RooUnfoldResponse.h"

  R__LOAD_LIBRARY(RooUnfold/libRooUnfold.so)
  R__LOAD_LIBRARY(fwdsmear.C+)



  #include "settings.h"

 void mk_fwdsmear() {


  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+g");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+g");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+g");

     
 fwdsmearUnfold(jp::type);

}
