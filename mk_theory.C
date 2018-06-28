// Purpose: Reformat theory curves
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2010
// Updated: June 1, 2010
{

  // compile code
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L theory.C+g");

  #include "settings.h"

  theory(jp::type); // new generic, use _algo inside
}
