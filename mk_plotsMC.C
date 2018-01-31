// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "tools.h"
#include "drawFracs.h"
#include "drawTrigFracs.h"
R__LOAD_LIBRARY(tools.C+g)
R__LOAD_LIBRARY(drawFracs.C+g)
R__LOAD_LIBRARY(drawTrigFracs.C+g)

void mk_plotsMC() {
    const char *path="./";
    string title="RunG";

    Fracs fracs(path,path,title,"pypdf");
    for (unsigned i = 0; i<4; ++i) fracs.drawFracs(i);
//     for (unsigned i = 0; i<4; ++i) drawFracs(i,path,path,title,"hwpdf","HW");
}
