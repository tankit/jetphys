#include "drawFracs.h"
R__LOAD_LIBRARY(drawFracs.C+g)

void mk_drawFracs_MCvMC() {
  const char *path1="./";
  const char *path2="./";
  const char *savedir="pdf";
  string name1="G";
  string name2="G";
  // Title format: Second vs First
  string title1="MC";
  string title2="HW";

  Fracs fracs(path1,path2,title1,savedir,false,name1,name2);
  Fracs tfracs(path1,path2,title2,savedir,true,name1,name2);

  for (unsigned i = 0; i<4; ++i) {
    fracs.drawFracs(i);
    tfracs.drawFracs(i);
  }
}
