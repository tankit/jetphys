// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_mcPlots for further options.

{
    gROOT->ProcessLine(".L tools.C");
    gROOT->ProcessLine(".L drawFracs.C");

    string path1="./";
    string path2="./";
    string name1="G";
    string name2="G";
    // Title format: Second vs First
    string title="RunG";

    drawFracs(0,path1,path2,title,"pdf",name1,name2);
    drawFracs(1,path1,path2,title,"pdf",name1,name2);
    drawFracs(2,path1,path2,title,"pdf",name2,name2);
}
