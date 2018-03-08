#include "DrawFracs.h"
R__LOAD_LIBRARY(DrawFracs.C+g)

void mk_DrawFracs_DTvDT() {
  const char *path1="./";
  const char *path2="./";
  const char *savedir="pdf";
  string name1="G";
  string name2="G";
  // Title format: Second vs First
  string title="RunG";

  Fracs fracs(path1,path2,title,savedir,false,name1,name2);
  Fracs tfracs(path1,path2,title,savedir,true,name1,name2);

  for (unsigned i = 0; i<4; ++i) {
    fracs.DrawFracs(i);
    tfracs.DrawFracs(i);
  }
}
