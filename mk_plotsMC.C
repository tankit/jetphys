// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_generalPlots for further options.

{
    gROOT->ProcessLine(".L tools.C");
    gROOT->ProcessLine(".L drawFracs.C");

    string path="./";
    string title="RunFgG";

    //drawFracs(0,path,path,title,"hwpdf","HW");
    //drawFracs(1,path,path,title,"hwpdf","HW");
    //drawFracs(2,path,path,title,"hwpdf","HW");

    drawFracs(0,path,path,title,"pypdf");
    drawFracs(1,path,path,title,"pypdf");
    drawFracs(2,path,path,title,"pypdf");
}
