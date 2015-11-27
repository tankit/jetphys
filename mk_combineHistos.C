// Purpose: Combine different triggers into a single spectrum
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: June 2, 2015
{

  // compile code
  gROOT->ProcessLine(".L combineHistos.C+");

  #include "settings.h"

  combineHistos(_jp_type);
}
