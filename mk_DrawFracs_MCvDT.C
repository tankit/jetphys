// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

// fracs uses output-*-2b.root and tfracs output-*-1.root
void mk_DrawFracs_MCvDT() {
    const char *path="./";
    const char *savedir="pypdf";
    //const char *savedir="hwpdf";
    string title="RunG";

    Fracs fracs(path,path,title,savedir,false,"MC","RunG",3);
    Fracs tfracs(path,path,title,savedir,true,"MC","RunG",3);
    for (unsigned i = 0; i<6; ++i) {
      fracs.DrawFracs(i);
      tfracs.DrawFracs(i);
    }
}
