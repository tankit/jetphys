// Purpose: Reformat theory curves
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2010
// Updated: June 1, 2010

  #include "settings.h"
  #include "tools.h"
 
  R__LOAD_LIBRARY(tools.C+g)
  R__LOAD_LIBRARY(theory.C+g)
  
  void mk_theory() {

  theory(jp::type); // new generic, use _algo inside
}
