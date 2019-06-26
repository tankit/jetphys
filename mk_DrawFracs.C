// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

// fracs uses output-*-2b.root and tfracs output-*-1.root
void mk_DrawFracs() {
    const char *mcpath="./";
    const char *dtpath="./";
    vector<const char*> savedirs = {"pypdf","nupdf","pynupdf"};
    vector<const char*> mctypes = {"MC","NU","MCNU"};
    vector<bool> dotfracs = {true,true,false};
    //vector<const char*> savedirs = {"pypdf","hwpdf","nupdf","pynupdf","hwnupdf"};
    //vector<const char*> mctypes = {"MC","HW","NU","MCNU","HWNU"};
    //vector<bool> dotfracs = {true,true,true,false,false};

    // 2016 IOV's: tot. 35.86 fb^{-1}
    //string title="Run2016BCD";
    //string lumit="12.6 fb^{-1}";
    //string title="Run2016EF";
    //string lumit="6.7 fb^{-1}";
    //string title="Run2016GH";
    //string lumit="16.5 fb^{-1}";

    // 2017 IOV's: tot. 41.52 fb^{-1}
    //string title="Run2017B";
    //string lumit="4.8 fb^{-1}";
    //string title="Run2017C";
    //string lumit="9.6 fb^{-1}";
    //string title="Run2017DE";
    //string lumit="13.6 fb^{-1}";
    //string title="Run2017F";
    //string lumit="13.5 fb^{-1}";

    // 2018 IOV's: tot. 53.27 fb^{-1}
    //string title="Run2018A";
    //string lumit="14.1 fb^{-1}";
    //string title="Run2018B";
    //string lumit="7.1 fb^{-1}";
    //string title="Run2018C";
    //string lumit="6.9 fb^{-1}";
    string title="Run2018D";
    string lumit="29.2 fb^{-1}";

    for (unsigned idx = 0; idx < mctypes.size(); ++idx) {
      const char *mctype = mctypes[idx];
      const char *savedir = savedirs[idx];
      bool do_t = dotfracs[idx];
      cout << mctype << endl;
      Fracs fracs(mcpath,dtpath,title,lumit,savedir,false,mctype,"DT",3);
      Fracs dfracs(mcpath,dtpath,title,lumit,savedir,false,mctype,"DT",3,"");
      vector<unsigned int> order = {0,3,1,2};//,4,5};
      for (auto &i : order)
        fracs.DrawFracs(i);
      dfracs.DrawFracs(0);
      if (do_t) {
        Fracs tfracs(mcpath,dtpath,title,lumit,savedir,true,mctype,"DT",3);
        for (auto &i : order)
          tfracs.DrawFracs(i);
      }
    }
}
