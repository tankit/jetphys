// Purpose: Run the creation of output-*2*.root files and plotting
// Author:  hannu.siikonen@cern.ch
// Created: July 4, 2019

// Normalization: compile code
R__LOAD_LIBRARY(HistosNormalize.C+g)
// Fuse: compile code
R__LOAD_LIBRARY(HistosFuse.C+g)
// Combination: compile code
R__LOAD_LIBRARY(HistosCombine.C+g)
// DrawFracs: plot the results
#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

void mk_HistosFullCircle() {
  // Quick configuration settings:
  bool dodt    = false;
  bool dopy    = false;
  bool dohw    = false;
  bool donu    = true;
  bool doht    = true;
  bool doall   = true;
  bool doplots = true;
  // Choose one: 0:16BCD,1:16EF,2:16GH,3:17B,4:17C,5:17DE,6:17F,7:18A,8:18B,9:18C,10:19D
  int lumidx = 2;

  // Output-*-2*.root files:
  //if (dodt) HistosNormalize("DATA");
  //if (dopy) HistosNormalize("MC");
  //if (dohw) HistosNormalize("HW");
  //if (donu) HistosNormalize("NU");
  //if (doht) HistosNormalize("HT");
  //if (donu and dopy) HistosFuse("MC");
  //if (donu and dohw) HistosFuse("HW");
  //if (donu and doht) HistosFuse("HT");
  //if (dodt)                      HistosCombine("DATA");
  //if ((!donu or doall) and dopy) HistosCombine("MC");
  //if ((!donu or doall) and dohw) HistosCombine("HW");
  //if ((!donu or doall) and dohw) HistosCombine("HW");
  //if (donu and doall)            HistosCombine("NU");
  //if (donu and dopy)             HistosCombine("MCNU");
  //if (donu and dohw)             HistosCombine("HWNU");
  //if (donu and doht)             HistosCombine("HTNU");

  // Plotting:
  if (doplots) {
    vector<const char*> savedirs;
    vector<const char*> mctypes;
    vector<bool> dotfracs;
    if (donu) {
      if (doall) {
        savedirs.push_back("nupdf");
        mctypes.push_back("NU");
        dotfracs.push_back(true);
      }
      if (dopy) {
        savedirs.push_back("pynupdf");
        mctypes.push_back("MCNU");
        dotfracs.push_back(false);
      }
      if (dohw) {
        savedirs.push_back("hwnupdf");
        mctypes.push_back("HWNU");
        dotfracs.push_back(false);
      }
      if (doht) {
        savedirs.push_back("htnupdf");
        mctypes.push_back("HTNU");
        dotfracs.push_back(false);
      }
    }
    if (!donu or doall) {
      if (dopy) {
        savedirs.push_back("pypdf");
        mctypes.push_back("MC");
        dotfracs.push_back(true);
      }
      if (dohw) {
        savedirs.push_back("hwpdf");
        mctypes.push_back("HW");
        dotfracs.push_back(true);
      }
      if (dohw) {
        savedirs.push_back("htpdf");
        mctypes.push_back("HT");
        dotfracs.push_back(true);
      }
    }

    // Picking luminosity info for the plots.
    string title, lumit;
    if (lumidx==0) {
      // 2016 IOV's: tot. 35.86 fb^{-1}
      title="Run2016BCD";
      lumit="12.6 fb^{-1}";
    } else if (lumidx==1) {
      title="Run2016EF";
      lumit="6.7 fb^{-1}";
    } else if (lumidx==2) {
      title="Run2016GH";
      lumit="16.5 fb^{-1}";
    } else if (lumidx==3) {
      // 2017 IOV's: tot. 41.52 fb^{-1}
      title="Run2017B";
      lumit="4.8 fb^{-1}";
    } else if (lumidx==4) {
      title="Run2017C";
      lumit="9.6 fb^{-1}";
    } else if (lumidx==5) {
      title="Run2017DE";
      lumit="13.6 fb^{-1}";
    } else if (lumidx==6) {
      title="Run2017F";
      lumit="13.5 fb^{-1}";
    } else if (lumidx==7) {
      // 2018 IOV's: tot. 53.27 fb^{-1}
      title="Run2018A";
      lumit="14.1 fb^{-1}";
    } else if (lumidx==8) {
      title="Run2018B";
      lumit="7.1 fb^{-1}";
    } else if (lumidx==9) {
      title="Run2018C";
      lumit="6.9 fb^{-1}";
    } else if (lumidx==10) {
      title="Run2018D";
      lumit="29.2 fb^{-1}";
    } else {
      cout << "Lumi index erroneous, check your value: " << lumidx << endl;
      return;
    }

    for (unsigned idx = 0; idx < mctypes.size(); ++idx) {
      const char *mctype = mctypes[idx];
      const char *savedir = savedirs[idx];
      bool do_t = dotfracs[idx];
      cout << mctype << endl;
      Fracs fracs("./","./",title,lumit,savedir,false,mctype,"DT",3);
      Fracs dfracs("./","./",title,lumit,savedir,false,mctype,"DT",3,"");
      vector<unsigned int> order = {0,3,1,2};//,4,5};
      for (auto &i : order)
        fracs.DrawFracs(i);
      dfracs.DrawFracs(0);
      if (do_t) {
        Fracs tfracs("./","./",title,lumit,savedir,true,mctype,"DT",3);
        for (auto &i : order)
          tfracs.DrawFracs(i);
      }
    }
  }
}
