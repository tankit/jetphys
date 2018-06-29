// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

// fracs uses output-*-2b.root and tfracs output-*-1.root
void mk_DrawFracs() {
    const char *mcpath="./";
    const char *dtpath="./";
    const char *savedir="pypdf";
    const char *mctype="MC";
    //const char *savedir="hwpdf";
    //const char *mctype="HW";
    //string title="RunFlGH";
    string title="Run2016GH";
    string lumit="16.5 fb^{-1}";

    Fracs fracs(mcpath,dtpath,title,lumit,savedir,false,mctype,"DT",3);
    //Fracs tfracs(mcpath,dtpath,title,lumit,savedir,true,mctype,"DT",3);
    vector<unsigned int> order = {0,3};//,1,2,4,5};
    for (auto &i : order) {
      fracs.DrawFracs(i);
      //tfracs.DrawFracs(i);
    }
}
