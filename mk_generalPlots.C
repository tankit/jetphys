// A script for producing MCvDT (pythia&herwig) plots quickly.
// Ref. mk_mcPlots for further options.

{
    gROOT->ProcessLine(".L tools.C");
    gROOT->ProcessLine(".L drawFracs.C");

    string path1="RR_BCD/";
    string path2="RR_EFearly/";
    string name1="BCD";
    string name2="EFe";
    // Title format: Second vs First
    string title="RunEFevsBCD";

    drawFracs(0,path1,path2,title,"pdf",name1,name2);
    drawFracs(1,path1,path2,title,"pdf",name1,name2);
    drawFracs(2,path1,path2,title,"pdf",name2,name2);
}
