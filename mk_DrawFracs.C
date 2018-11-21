// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

// fracs uses output-*-2b.root and tfracs output-*-1.root
void mk_DrawFracs() {
    const char *mcpath="./";
    const char *dtpath="./";
    vector<const char*> savedirs = {"pypdf","hwpdf"}; // pypdf, hwpdf
    vector<const char*> mctypes = {"MC","HW"}; // MC, HW

    // 2016 IOV's: tot. 35.86 fb^{-1}
    //string title="Run2016BCD";
    //string lumit="12.6 fb^{-1}";
    //string title="Run2016EF";
    //string lumit="6.7 fb^{-1}";
    string title="Run2016GH";
    string lumit="16.5 fb^{-1}";

    for (unsigned idx = 0; idx < mctypes.size(); ++idx) {
      const char *mctype = mctypes[idx];
      const char *savedir = savedirs[idx];
      Fracs fracs(mcpath,dtpath,title,lumit,savedir,false,mctype,"DT",3);
      Fracs tfracs(mcpath,dtpath,title,lumit,savedir,true,mctype,"DT",3);
      Fracs dfracs(mcpath,dtpath,title,lumit,savedir,false,mctype,"DT",3,"");
      vector<unsigned int> order = {0,3};//,1,2,4,5};
      for (auto &i : order) {
        fracs.DrawFracs(i);
        tfracs.DrawFracs(i);
      }
      dfracs.DrawFracs(0);
    }
}
