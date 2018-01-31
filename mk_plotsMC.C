// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "tools.h"
#include "drawFracs.h"
R__LOAD_LIBRARY(tools.C+g)
R__LOAD_LIBRARY(drawFracs.C+g)

void mk_plotsMC() {
    string path="./";
    string title="RunG";

    for (unsigned i = 0; i<4; ++i) drawFracs(i,path,path,title,"pypdf");
//     for (unsigned i = 0; i<4; ++i) drawFracs(i,path,path,title,"hwpdf","HW");
}
