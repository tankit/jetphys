// Purpose: Reformat theory curves
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2010
// Updated: June 1, 2010
{

  // compile code
  gROOT->ProcessLine(".L tools.C+");
  gROOT->ProcessLine(".L theory.C+");

  #include "settings.h"

  theory(_jp_type); // new generic, use _algo inside
}
