#include "tdrstyle_mod20.C"

void drawMultijetUL18(TString dtfile="2018ABCD", int jer = 0, bool usel2res = false, int rebin = 2, TString pt="30", TString jt="")
{

  //gROOT->ForceStyle(kTRUE);
  //gStyle->SetImageScaling(3.);

  bool UL18IOVs = false; //for both 20201111 and 20201119 (before and after fixing MPFue_leading_MG and MPF0_leading_MG, respectively) 
  //bool UL18IOVs = (dtfile=="2018AB" or dtfile=="2018CD" or dtfile=="2018ABCD") ? false : true; //reweight each IOV separately
  if(jer!=0 and dtfile!="2018ABCD") { cout << " No MC inputs with JER up/down for " << dtfile << ". Use 2018ABCD! " << endl; return 0; }

  TString ver = "_JECV5JRV2";
  TString rebinVer = "_rebinV37";
  TString doVeto = "_doVetoHot"; if(dtfile=="2018C" || dtfile=="2018D" || dtfile=="2018CD" || dtfile=="2018ABCD") doVeto = "_doVetoHEM";

  bool met2 = false; //met2 = true; //metType-II for the plot (automatically store with metType-I)
  bool ptave_rebin = false; //ptave_rebin = true; //20200423 test version only
  bool ptlead_mcDir = true; //20200424 version

  int jtpt = 15; if(pt=="20") jtpt = 20; else if(pt=="30") jtpt = 30;

  TString mcfile = "MGP8CP5";
  TString sample = "MGP8CP5"+ver, suffix = "_MG";

  bool drawMC = false; //since Mar 27, 2020
  bool noRebin = false; //true; //for test purpose (Mar 27, 2020)

//  bool usemcDir = false; //NOT USE ANY MORE from Sep 13 !! //usemcDir = true;
  bool usemcDir = true; //This is for recoil, so recoil comes from mcDir (Mar 30, 2020) -- be careful!
  TString useSampleNo = "2b"; if(usemcDir) useSampleNo = "1";
  cout << "***********************************************" << endl;
  if(usemcDir) cout << "*** Using mcDir to get MC Recoil   (default)     ***" << endl;
  else          cout << "*** Using trigDir to get MC Recoil (old)  ***" << endl;
  cout << "***********************************************" << endl;
  cout << "*** Use mcDir to get MC ptave from Apr 20, 2020" << endl;

  bool leading_rebin = false; //leading_rebin = true; //for test only (no more)

  bool writeOutput = false, writeMC = true;
  //writeOutput = true;
  TString writeDate = "Rebin2_20201119"; //"Rebin2_20201111";

  if(writeOutput) {
    cout << " Making input files for global fits"; if(writeMC) cout << " and storing the histos from MC";
    cout << endl;
  } else cout << " Drawing plots only (If you want to make input files, then writeOutput = true) " << endl;

  TString fileName = "";

  TString jetCorrVer = "_jecV5";

  printf(" Run Smmmary\n");
  printf(" -----------\n");
  printf(" Jet ID: TightLepVeto with pT > %s GeV and |eta| < 2.5\n",pt.Data());
  TString jerVer = "Summer19UL18_JRV2";
  if(jer==0) printf(" 1) JER %s with SF = %s\n",jerVer.Data(),"Normal"); else printf(" JER Version: %s with SF = %s\n",jerVer.Data(),(jer>0) ? "Up" : "Down");
  printf(" 2) For Data only: %s\n",(usel2res) ? "only_L2Res" : "L2L3Res (default)");
  printf(" 3) MET type-%s\n",(met2) ? "II" : "I (default)");
  printf(" => jetCorrVer = %s\n",jetCorrVer.Data());
  printf("\n");

  //TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0.0001); 

  bool mjbinv = true, mjbinvDraw = false;
  bool mpfinv = false, mpfinvDraw = false;
  bool mpflead = true, mpfleadDraw = false;
  TString metType = "met1"; if(met2) metType = "met2";
  bool both = true; //both = false;
  bool recoil = true, leading = false;
  //leading = true, recoil = false;
  if (both) recoil = false, leading = false;
  bool usehist2 = false, usehist2Draw = false, mpfprofDraw = true; //true; //new from Mar 27, 2020
  usehist2 = true; //usehist2Draw = true; //false from Mar 28, 2020 !!!

  cout << " ################################################## " << endl;
  cout << " use h2mjbinv_leading   " << endl;
  cout << " use h2mpflead_leading  " << endl;
  cout << "     (pmpflead_leading) " << endl;
  if(mpfprofDraw) cout << endl << " *** CAUTION! usehist2Draw for MPF = false (new from Mar 27, 2020) " << endl << endl;
  cout << " ################################################## " << endl;

  TString ul18 = "_Summer19UL18";
  TString l2res = ""; if(usel2res) l2res = "_usel2res";
  TString jerSF = ""; if(jer>0) jerSF = "_jerUp"; else if(jer<0) jerSF = "_jerDown";
  if(UL18IOVs and dtfile!="2018ABCD") { jerSF = jerSF+"_run"+dtfile; cout << " MC Pilup reweighted by each IOV separately " << endl; } //different IOV (Apr 20, 2020);

  dtfile = "DATA-2b_Run"+dtfile+ver;
 
  TString ptCut = "152030", ipt = "";
  if(pt=="20") { ptCut = "15"+pt+"30"; ipt = "1"; } else if(pt=="15") { ptCut = pt+"2030"; ipt = "0"; }

  double xmin = 90., xmax = 3600.; //xmin = 80;
  xmin = 272; xmin = 245; xmin = 220; xmin = 196; xmin = 100;
  double hmin = 0.7, hmax = 1.3;
  double hminRat = 0.5, hmaxRat = 1.5; hminRat = 0.75, hmaxRat = 1.25; hminRat = 0.9; hmaxRat = 1.1; //hminRat = 0.95, hmaxRat = 1.05; //if change hminRat, then change LabelOffset, too!
  double leg_ymin = 0.16, leg_ymax = 0.31, leg_xmin = 0.52, leg_xmax = 0.75;
  TString str_mc = "MJB, MGP8", str_mc2 = "MPF, MGP8";
  TString str_dt = "2018A";
  if(dtfile.Contains("2018B")) str_dt = "2018B"; else if(dtfile.Contains("2018C")) str_dt = "2018C"; else if(dtfile.Contains("2018D")) str_dt = "2018D"; 
  if(dtfile.Contains("2018AB")) str_dt = "2018AB"; else if(dtfile.Contains("2018CD")) str_dt = "2018CD";
  if(dtfile.Contains("2018ABCD")) str_dt = "2018ABCD";

//  xmin = 300, leg_xmin = 0.4, leg_xmax = 0.63; //Oct 22, 2019
  xmin = 15; //32; //49; //64;
  if(xmin<100) {
    hmin = 0.68, hminRat = 0.8;
    //if(xmin<50) hmin = 0.66, hminRat = 0.65, hmaxRat = 1.25;
    if(xmin<50) hmax = 1.25, hmin = 0.45, hminRat = 0.68, hmaxRat = 1.32;
    hmin = 0.5;
    //if(xmin==32) hmin = 0.56, hminRat = 0.9;
    //leg_xmin = 0.54, leg_xmax = 0.77;
  }
  if(writeOutput) hmin = 0, hmax = 2;
  hminRat = 0.95, hmaxRat = 1.05;
  if(pt!="30") hminRat = 0.7, hmaxRat = 1.3;

  setTDRStyle(); leg_ymin = 0.06, leg_ymax = 0.23;
  if(!drawMC) leg_ymin = 0.14; 

  TString mcDir = "Standard/Eta_0.0-1.3/mc";
  TString trigDir = "Standard/Eta_0.0-1.3"+jt;
  cout << " gDirectory->cd(" << trigDir << ")" << endl;

//gDirectory->cd("Standard/Eta_0.0-1.3/jt500");
//gDirectory->cd(trigDir);
TProfile *pmjb_recoil = 0, *pmjb_leading = 0, *pmpf_recoil = 0, *pmpf_leading = 0, *pmpf_recoil2 = 0, *pmpf_leading2 = 0, *pmjbinv_leading = 0, *pmpfinv_leading = 0, *pmpflead_leading = 0, *pmpflead_leading2 = 0;
TProfile *pmpfleadinv_leading = 0, *pmpfleadinv_leading2 = 0; TH2D *h2mpfleadinv_leading = 0, *h2mpfleadinv_leading2 = 0;
//TProfile *h2mpfleadinv_leading_pfx = 0, *h2mpfleadinv_leading2_pfx = 0; TH1D *hmpfleadinv_leading = 0, *hmpfleadinv_leading2 = 0;
TProfile *pmjb_ptave = 0, *pmpf_ptave = 0, *pmpf_ptave2 = 0; TH2D *h2mjb_ptave = 0, *h2mpf_ptave = 0, *h2mpf_ptave2 = 0; TH1D *hdt_mjb_ptave = 0, *hdt_mpf_ptave = 0, *hdt_mpf_ptave2 = 0;
TProfile *pmjbinv_ptave = 0, *pmpfinv_ptave = 0, *pmpfinv_ptave2 = 0; TH2D *h2mjbinv_ptave = 0, *h2mpfinv_ptave = 0, *h2mpfinv_ptave2 = 0; //exist in HistosCombine from UL17new or UL17ave
TH1D *hdt_mjbinv_ptave = 0, *hdt_mpfinv_ptave = 0, *hdt_mpfinv_ptave2 = 0;
TProfile *pmpflead_leading_ue = 0; TH2D *h2mpflead_leading_ue = 0; TH1D *hdt_mpflead_leading_ue = 0;
TProfile *pmpfleadinv_leading_ue = 0; TH2D *h2mpfleadinv_leading_ue = 0; TH1D *hdt_mpfleadinv_leading_ue = 0;
////TProfile *pmpflead_ptave_ue = 0; TH2D *h2mpflead_ptave_ue = 0; TH1D *hdt_mpflead_ptave_ue = 0;
TH2D *h2mjb_recoil = 0, *h2mpf_recoil = 0, *h2mpf_recoil2 = 0, *h2mjb_leading = 0, *h2mjbinv_leading = 0, *h2mpf_leading = 0, *h2mpf_leading2 = 0;
TFile *input = 0;

TString inputName = "./jetphys_output/output-"+dtfile+"_dPhiPt"+ptCut+doVeto+"_MetJetPt15"+ul18+l2res+".root";

cout << "Data        : " << inputName << endl;
if (! gSystem->AccessPathName(inputName)) input = new TFile(inputName.Data());
else { cout << "Opening the input file failed. Check it the Data file exists." << endl; return 0; } 

input->GetObject(trigDir+"/pmjb_recoil"+ipt,pmjb_recoil); pmjb_recoil->Rebin(rebin);
input->GetObject(trigDir+"/pmjb_leading"+ipt,pmjb_leading); pmjb_leading->Rebin(rebin);
input->GetObject(trigDir+"/pmpf_recoil"+ipt,pmpf_recoil); pmpf_recoil->Rebin(rebin);
input->GetObject(trigDir+"/pmpf_recoil2"+ipt,pmpf_recoil2); pmpf_recoil2->Rebin(rebin);
input->GetObject(trigDir+"/pmjbinv_leading"+ipt,pmjbinv_leading); pmjbinv_leading->Rebin(rebin);
input->GetObject(trigDir+"/pmpfinv_leading"+ipt,pmpfinv_leading); pmpfinv_leading->Rebin(rebin); //no more using (confirmed on Sep 16, 2020)
input->GetObject(trigDir+"/pmpflead_leading"+ipt,pmpf_leading); pmpf_leading->Rebin(rebin); //testing -> using (confirmed on Sep 16, 2020)
input->GetObject(trigDir+"/pmpflead_leading2"+ipt,pmpf_leading2); pmpf_leading2->Rebin(rebin); //testing -> using (confirmed on Sep 16, 2020)

  input->GetObject(trigDir+"/pmpfleadinv_leading"+ipt,pmpfleadinv_leading); pmpfleadinv_leading->Rebin(rebin);
  input->GetObject(trigDir+"/pmpfleadinv_leading2"+ipt,pmpfleadinv_leading2); pmpfleadinv_leading2->Rebin(rebin);
  input->GetObject(trigDir+"/h2mpfleadinv_leading"+ipt,h2mpfleadinv_leading);
  input->GetObject(trigDir+"/h2mpfleadinv_leading2"+ipt,h2mpfleadinv_leading2);
  //h2mpfleadinv_leading_pfx = h2mpfleadinv_leading->ProfileX(); h2mpfleadinv_leading_pfx->Rebin(rebin); hmpfleadinv_leading = h2mpfleadinv_leading_pfx->ProjectionX();
  //h2mpfleadinv_leading2_pfx = h2mpfleadinv_leading2->ProfileX(); h2mpfleadinv_leading2_pfx->Rebin(rebin); hmpfleadinv_leading2 = h2mpfleadinv_leading2_pfx->ProjectionX();
  if(!ptave_rebin) {
    input->GetObject(trigDir+"/pmjb_ptave"+ipt,pmjb_ptave); pmjb_ptave->Rebin(rebin); hdt_mjb_ptave = pmjb_ptave->ProjectionX(); hdt_mjb_ptave->SetStats(0); hdt_mjb_ptave->SetMarkerStyle(1);
    input->GetObject(trigDir+"/pmpf_ptave"+ipt,pmpf_ptave); pmpf_ptave->Rebin(rebin); hdt_mpf_ptave = pmpf_ptave->ProjectionX(); hdt_mpf_ptave->SetStats(0); hdt_mpf_ptave->SetMarkerStyle(1);
    input->GetObject(trigDir+"/pmpf_ptave2"+ipt,pmpf_ptave2); pmpf_ptave2->Rebin(rebin); hdt_mpf_ptave2 = pmpf_ptave2->ProjectionX(); hdt_mpf_ptave2->SetStats(0); hdt_mpf_ptave2->SetMarkerStyle(1); 
    input->GetObject(trigDir+"/h2mjb_ptave"+ipt,h2mjb_ptave); //h2mjb_ptave->Rebin(rebin);
    input->GetObject(trigDir+"/h2mpf_ptave"+ipt,h2mpf_ptave); //h2mpf_ptave->Rebin(rebin);
    input->GetObject(trigDir+"/h2mpf_ptave2"+ipt,h2mpf_ptave2); //h2mpf_ptave2->Rebin(rebin);
    input->GetObject(trigDir+"/pmjbinv_ptave"+ipt,pmjbinv_ptave); pmjbinv_ptave->Rebin(rebin); hdt_mjbinv_ptave = pmjbinv_ptave->ProjectionX(); hdt_mjbinv_ptave->SetStats(0); hdt_mjbinv_ptave->SetMarkerStyle(1);
    input->GetObject(trigDir+"/pmpfinv_ptave"+ipt,pmpfinv_ptave); pmpfinv_ptave->Rebin(rebin); hdt_mpfinv_ptave = pmpfinv_ptave->ProjectionX(); hdt_mpfinv_ptave->SetStats(0); hdt_mpfinv_ptave->SetMarkerStyle(1);
    input->GetObject(trigDir+"/pmpfinv_ptave2"+ipt,pmpfinv_ptave2); pmpfinv_ptave2->Rebin(rebin); hdt_mpfinv_ptave2 = pmpfinv_ptave2->ProjectionX(); hdt_mpfinv_ptave2->SetStats(0); hdt_mpfinv_ptave2->SetMarkerStyle(1);
    input->GetObject(trigDir+"/h2mjbinv_ptave"+ipt,h2mjbinv_ptave); //h2mjbinv_ptave->Rebin(rebin);
    input->GetObject(trigDir+"/h2mpfinv_ptave"+ipt,h2mpfinv_ptave); //h2mpfinv_ptave->Rebin(rebin);
    input->GetObject(trigDir+"/h2mpfinv_ptave2"+ipt,h2mpfinv_ptave2); //h2mpfinv_ptave2->Rebin(rebin);
  }
  input->GetObject(trigDir+"/pmpflead_leading_ue"+ipt,pmpflead_leading_ue); pmpflead_leading_ue->Rebin(rebin); hdt_mpflead_leading_ue = pmpflead_leading_ue->ProjectionX();
  input->GetObject(trigDir+"/pmpfleadinv_leading_ue"+ipt,pmpfleadinv_leading_ue); pmpfleadinv_leading_ue->Rebin(rebin); hdt_mpfleadinv_leading_ue = pmpfleadinv_leading_ue->ProjectionX();
  input->GetObject(trigDir+"/h2mpflead_leading_ue"+ipt,h2mpflead_leading_ue);
  input->GetObject(trigDir+"/h2mpfleadinv_leading_ue"+ipt,h2mpfleadinv_leading_ue);
  hdt_mpflead_leading_ue->SetStats(0); hdt_mpflead_leading_ue->SetMarkerStyle(1); hdt_mpfleadinv_leading_ue->SetStats(0); hdt_mpfleadinv_leading_ue->SetMarkerStyle(1);
  ////input->GetObject(trigDir+"/pmpflead_ptave_ue"+ipt,pmpflead_ptave_ue); pmpflead_ptave_ue->Rebin(rebin); hdt_mpflead_ptave_ue = pmpflead_ptave_ue->ProjectionX();
  ////input->GetObject(trigDir+"/h2mpflead_ptave_ue"+ipt,h2mpflead_ptave_ue);
  ////hdt_mpflead_ptave_ue->SetStats(0); hdt_mpflead_ptave_ue->SetMarkerStyle(1); 

if(usehist2) {
  input->GetObject(trigDir+"/h2mjb_leading"+ipt,h2mjb_leading);
  input->GetObject(trigDir+"/h2mjbinv_leading"+ipt,h2mjbinv_leading); //testing
  input->GetObject(trigDir+"/h2mpflead_leading"+ipt,h2mpf_leading); //testing
  input->GetObject(trigDir+"/h2mpflead_leading2"+ipt,h2mpf_leading2); //testing
}

TProfile *pmjb_recoil_rebin = 0, *pmjb_leading_rebin = 0, *pmpf_recoil_rebin = 0, *pmpf_leading_rebin = 0, *pmpf_recoil2_rebin = 0, *pmpf_leading2_rebin = 0;
TProfile *pmjbinv_leading_rebin = 0;
TH2D *h2mjb_recoil_rebin = 0, *h2mpf_recoil_rebin = 0, *h2mpf_recoil2_rebin = 0;
TFile *inputRebin = 0;

TString inputRebinName = "./jetphys_output/output-"+dtfile+"_dPhiPt"+ptCut+doVeto+"_MetJetPt15"+ul18+l2res+rebinVer+".root";

cout << "Data Rebin  : " << inputRebinName << endl;
if (! gSystem->AccessPathName(inputRebinName)) inputRebin = new TFile(inputRebinName.Data());
else { cout << "Opening the input file failed. Check it the Data Rebin file exists." << endl; return 0; }

inputRebin->GetObject(trigDir+"/pmjb_recoil"+ipt,pmjb_recoil_rebin); pmjb_recoil_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmjb_leading"+ipt,pmjb_leading_rebin); pmjb_leading_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmpf_recoil"+ipt,pmpf_recoil_rebin); pmpf_recoil_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmpf_leading"+ipt,pmpf_leading_rebin); pmpf_leading_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmpf_recoil2"+ipt,pmpf_recoil2_rebin); pmpf_recoil2_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmpf_leading2"+ipt,pmpf_leading2_rebin); pmpf_leading2_rebin->Rebin(rebin);
inputRebin->GetObject(trigDir+"/pmjbinv_leading"+ipt,pmjbinv_leading_rebin); pmjbinv_leading_rebin->Rebin(rebin);
if(ptave_rebin) {
  inputRebin->GetObject(trigDir+"/pmjb_ptave"+ipt,pmjb_ptave); pmjb_ptave->Rebin(rebin); hdt_mjb_ptave = pmjb_ptave->ProjectionX(); hdt_mjb_ptave->SetStats(0); hdt_mjb_ptave->SetMarkerStyle(1);
  inputRebin->GetObject(trigDir+"/pmpf_ptave"+ipt,pmpf_ptave); pmpf_ptave->Rebin(rebin); hdt_mpf_ptave = pmpf_ptave->ProjectionX(); hdt_mpf_ptave->SetStats(0); hdt_mpf_ptave->SetMarkerStyle(1);
  inputRebin->GetObject(trigDir+"/pmpf_ptave2"+ipt,pmpf_ptave2); pmpf_ptave2->Rebin(rebin); hdt_mpf_ptave2 = pmpf_ptave2->ProjectionX(); hdt_mpf_ptave2->SetStats(0); hdt_mpf_ptave2->SetMarkerStyle(1);
  inputRebin->GetObject(trigDir+"/h2mjb_ptave"+ipt,h2mjb_ptave); //h2mjb_ptave->Rebin(rebin);
  inputRebin->GetObject(trigDir+"/h2mpf_ptave"+ipt,h2mpf_ptave); //h2mpf_ptave->Rebin(rebin);
  inputRebin->GetObject(trigDir+"/h2mpf_ptave2"+ipt,h2mpf_ptave2); //h2mpf_ptave2->Rebin(rebin);
  inputRebin->GetObject(trigDir+"/pmjbinv_ptave"+ipt,pmjbinv_ptave); pmjbinv_ptave->Rebin(rebin); hdt_mjbinv_ptave = pmjbinv_ptave->ProjectionX(); hdt_mjbinv_ptave->SetStats(0); hdt_mjbinv_ptave->SetMarkerStyle(1); 
  inputRebin->GetObject(trigDir+"/pmpfinv_ptave"+ipt,pmpfinv_ptave); pmpfinv_ptave->Rebin(rebin); hdt_mpfinv_ptave = pmpfinv_ptave->ProjectionX(); hdt_mpfinv_ptave->SetStats(0); hdt_mpfinv_ptave->SetMarkerStyle(1); 
  inputRebin->GetObject(trigDir+"/pmpfinv_ptave2"+ipt,pmpfinv_ptave2); pmpfinv_ptave2->Rebin(rebin); hdt_mpfinv_ptave2 = pmpfinv_ptave2->ProjectionX(); hdt_mpfinv_ptave2->SetStats(0); hdt_mpfinv_ptave2->SetMarkerStyle(1); 
  inputRebin->GetObject(trigDir+"/h2mjbinv_ptave"+ipt,h2mjbinv_ptave); //h2mjbinv_ptave->Rebin(rebin);
  inputRebin->GetObject(trigDir+"/h2mpfinv_ptave"+ipt,h2mpfinv_ptave); //h2mpfinv_ptave->Rebin(rebin);
  inputRebin->GetObject(trigDir+"/h2mpfinv_ptave2"+ipt,h2mpfinv_ptave2); //h2mpfinv_ptave2->Rebin(rebin);
}
if(usehist2) {
  inputRebin->GetObject(trigDir+"/h2mjb_recoil"+ipt,h2mjb_recoil_rebin);
  inputRebin->GetObject(trigDir+"/h2mpf_recoil"+ipt,h2mpf_recoil_rebin);
  inputRebin->GetObject(trigDir+"/h2mpf_recoil2"+ipt,h2mpf_recoil2_rebin);
}

//-- output-MC in the directory of jetphys_output since Apr 19, 2020
TString inputMCName     = "./jetphys_output/output-MC-1_"+sample+"_dPhiPt"+ptCut+"_smeared_doVetoHot_MetJetPt15"+ul18+jerSF+".root";
TString inputMCtrigName = "./jetphys_output/output-MC-2b_"+sample+"_dPhiPt"+ptCut+"_smeared_doVetoHot_MetJetPt15"+ul18+jerSF+".root";
TString inputCRecoilMCName = "./jetphys_output/output-MC-"+useSampleNo+"_"+sample+"_dPhiPt"+ptCut+"_smeared_doVetoHot_MetJetPt15"+ul18+jerSF+".root";

cout << "MC (CRecoil): " << inputCRecoilMCName << endl;
//--Do I need to update as follows? (Mar 30, 2020)
//cout << "MC (CRecoil_recoil ): " << inputMCName << endl;
//cout << "MC (CRecoil_leading): " << inputMCtrigName << endl;
cout << "MC          : " << inputMCName << endl;
cout << "MC Trigger  : " << inputMCtrigName << endl;

TProfile *pCRecoil_dt, *pCRecoil_mc, *pCRecoil_ptcl, *pCRecoil_ptcl2;
TProfile *pCRecoil_leading_dt, *pCRecoil_ptave_dt;
TProfile *pCRecoil_leading_mc, *pCRecoil_ptave_mc; //not yet
TFile *inputCRecoil = 0, *inputCRecoilMC = 0, *inputCRecoilMCtrig;

inputCRecoil = new TFile(inputName.Data()); //this is what we used in 2018 (test for UL17: /Users/minsuk/hip/test_20200330)
cout << " *** pCRecil uses pT rebin for data from Apr 20, 2020 and output-MC-1 (inputCRecoilMC/mcDir) for MC *** " << endl;
//-- pCRecoil_leading uses pT bin for data and output-MC-2b (inputMCtrig/trigDir) for MC
//-- pCRecoil_ptave uses pT bin for data and output-MC-1 (inputCRecoilMC/mcDir) for MC (hybrid?)
//inputCRecoil = new TFile(inputRebinName.Data()); //updated because we are using rebin for recoil (Mar 12, 2020) until 20200328 (trigDir) and test2_20200330 (usemcDir)
//inputCRecoil->GetObject(trigDir+"/pCRecoil"+ipt,pCRecoil_dt); pCRecoil_dt->Rebin(rebin);

  inputCRecoil->GetObject(trigDir+"/pCRecoil_leading"+ipt,pCRecoil_leading_dt); pCRecoil_leading_dt->Rebin(rebin);
  if(!ptave_rebin) { inputCRecoil->GetObject(trigDir+"/pCRecoil_ptave"+ipt,pCRecoil_ptave_dt), pCRecoil_ptave_dt->Rebin(rebin); }


//TProfile *pCRecoil_check = 0;
//TFile *inputCRecoilCheck = new TFile(inputName.Data());
TFile *inputRebinCRecoil = new TFile(inputRebinName.Data());
inputRebinCRecoil->GetObject(trigDir+"/pCRecoil"+ipt,pCRecoil_dt); pCRecoil_dt->Rebin(rebin);
if(ptave_rebin) { inputRebinCRecoil->GetObject(trigDir+"/pCRecoil_ptave"+ipt,pCRecoil_ptave_dt), pCRecoil_ptave_dt->Rebin(rebin); } //tested on Apr 23, 2020

TFile *inputMCtrig = 0;
if (! gSystem->AccessPathName(inputMCtrigName)) inputMCtrig = new TFile(inputMCtrigName.Data());
else { cout << "Opening the input file failed. Check it the MC Trigger file exists." << endl; return 0; }

if (! gSystem->AccessPathName(inputCRecoilMCName)) {
  inputCRecoilMC = new TFile(inputCRecoilMCName.Data());
  if(usemcDir) inputCRecoilMC->GetObject(mcDir+"/pCRecoil"+ipt,pCRecoil_mc); else inputCRecoilMC->GetObject(trigDir+"/pCRecoil"+ipt,pCRecoil_mc);
  if(usemcDir) inputCRecoilMC->GetObject(mcDir+"/pCRecoil_ptcl"+ipt,pCRecoil_ptcl); else inputCRecoilMC->GetObject(trigDir+"/pCRecoil_ptcl"+ipt,pCRecoil_ptcl);
  if(usemcDir) inputCRecoilMC->GetObject(mcDir+"/pCRecoil_ptcl2"+ipt,pCRecoil_ptcl2); else inputCRecoilMC->GetObject(trigDir+"/pCRecoil_ptcl2"+ipt,pCRecoil_ptcl2);
  if(usemcDir) inputCRecoilMC->GetObject(mcDir+"/pCRecoil_leading"+ipt,pCRecoil_leading_mc); else inputCRecoilMC->GetObject(trigDir+"/pCRecoil_leading"+ipt,pCRecoil_leading_mc);
  if(!ptlead_mcDir) inputMCtrig->GetObject(trigDir+"/pCRecoil_leading"+ipt,pCRecoil_leading_mc); //updated on Apr 20, 2020 -> is the previous way better? (Apr 24, 2020)
  if(usemcDir) inputCRecoilMC->GetObject(mcDir+"/pCRecoil_ptave"+ipt,pCRecoil_ptave_mc); else inputCRecoilMC->GetObject(trigDir+"/pCRecoil_ptave"+ipt,pCRecoil_ptave_mc);
  pCRecoil_mc->Rebin(rebin); pCRecoil_ptcl->Rebin(rebin); pCRecoil_ptcl2->Rebin(rebin); pCRecoil_leading_mc->Rebin(rebin); pCRecoil_ptave_mc->Rebin(rebin);
} else { cout << "Opening the input file failed. Check it the MC CRecoil file exists." << endl; return 0; }

TProfile *pmjb_recoil_mc = 0, *pmjb_leading_mc = 0, *pmpf_recoil_mc = 0, *pmpf_leading_mc = 0, *pmpf_recoil2_mc = 0, *pmpf_leading2_mc = 0, *pmjbinv_leading_mc = 0, *pmpfinv_leading_mc = 0;
TProfile *pmpfleadinv_leading_mc = 0, *pmpfleadinv_leading2_mc = 0; TH2D *h2mpfleadinv_leading_mc = 0, *h2mpfleadinv_leading2_mc = 0;
TProfile *pmjb_ptave_mc = 0, *pmpf_ptave_mc = 0, *pmpf_ptave2_mc = 0; TH2D *h2mjb_ptave_mc = 0, *h2mpf_ptave_mc = 0, *h2mpf_ptave2_mc = 0; TH1D *hmc_mjb_ptave = 0, *hmc_mpf_ptave = 0, *hmc_mpf_ptave2 = 0;
TProfile *pmjbinv_ptave_mc = 0, *pmpfinv_ptave_mc = 0, *pmpfinv_ptave2_mc = 0; TH2D *h2mjbinv_ptave_mc = 0, *h2mpfinv_ptave_mc = 0, *h2mpfinv_ptave2_mc = 0; //not in HistosCombine
TH1D *hmc_mjbinv_ptave = 0, *hmc_mpfinv_ptave = 0, *hmc_mpfinv_ptave2 = 0;
TProfile *pmpflead_leading_ue_mc = 0; TH2D *h2mpflead_leading_ue_mc = 0; TH1D *hmc_mpflead_leading_ue = 0;
TProfile *pmpfleadinv_leading_ue_mc = 0; TH2D *h2mpfleadinv_leading_ue_mc = 0; TH1D *hmc_mpfleadinv_leading_ue = 0;
////TProfile *pmpflead_ptave_ue_mc = 0; TH2D *h2mpflead_ptave_ue_mc = 0; TH1D *hmc_mpflead_ptave_ue = 0;
TH2D *h2mjb_recoil_mc = 0, *h2mpf_recoil_mc = 0, *h2mpf_recoil2_mc = 0, *h2mjb_leading_mc = 0, *h2mjbinv_leading_mc = 0, *h2mpf_leading_mc = 0, *h2mpf_leading2_mc = 0;
TFile *inputMC = 0;
if (! gSystem->AccessPathName(inputMCName)) {
  inputMC = new TFile(inputMCName.Data());
} else { cout << "Opening the input file failed. Check it the MC file exists." << endl; return 0; }
inputMC->GetObject(mcDir+"/pmjb_recoil"+ipt,pmjb_recoil_mc); pmjb_recoil_mc->Rebin(rebin);
inputMC->GetObject(mcDir+"/pmpf_recoil"+ipt,pmpf_recoil_mc); pmpf_recoil_mc->Rebin(rebin);
inputMC->GetObject(mcDir+"/pmpf_recoil2"+ipt,pmpf_recoil2_mc); pmpf_recoil2_mc->Rebin(rebin);
if(usehist2) {
  //--need for error calculation in global fit (wider bin) Sep 20
  inputMC->GetObject(mcDir+"/h2mjb_recoil"+ipt,h2mjb_recoil_mc);
  inputMC->GetObject(mcDir+"/h2mpf_recoil"+ipt,h2mpf_recoil_mc);
  inputMC->GetObject(mcDir+"/h2mpf_recoil2"+ipt,h2mpf_recoil2_mc);
}

//--updated on Apr 24, 2020
if(ptlead_mcDir) {
  inputMC->GetObject(mcDir+"/pmjb_leading"+ipt,pmjb_leading_mc); pmjb_leading_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/pmjbinv_leading"+ipt,pmjbinv_leading_mc); pmjbinv_leading_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/pmpfinv_leading"+ipt,pmpfinv_leading_mc); pmpfinv_leading_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/pmpflead_leading"+ipt,pmpf_leading_mc); pmpf_leading_mc->Rebin(rebin); //testing -> using (confirmed on Sep 16, 2020)
  inputMC->GetObject(mcDir+"/pmpflead_leading2"+ipt,pmpf_leading2_mc); pmpf_leading2_mc->Rebin(rebin); //testing -> using (confirmed on Sep 16, 2020)

  inputMC->GetObject(mcDir+"/pmpfleadinv_leading"+ipt,pmpfleadinv_leading_mc); pmpfleadinv_leading_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/pmpfleadinv_leading2"+ipt,pmpfleadinv_leading2_mc); pmpfleadinv_leading2_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading"+ipt,h2mpfleadinv_leading_mc);
  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading2"+ipt,h2mpfleadinv_leading2_mc);
  inputMC->GetObject(mcDir+"/pmpflead_leading_ue"+ipt,pmpflead_leading_ue_mc); pmpflead_leading_ue_mc->Rebin(rebin); hmc_mpflead_leading_ue = pmpflead_leading_ue_mc->ProjectionX();
  inputMC->GetObject(mcDir+"/h2mpflead_leading_ue"+ipt,h2mpflead_leading_ue_mc);
  inputMC->GetObject(mcDir+"/pmpfleadinv_leading_ue"+ipt,pmpfleadinv_leading_ue_mc); pmpfleadinv_leading_ue_mc->Rebin(rebin); hmc_mpfleadinv_leading_ue = pmpfleadinv_leading_ue_mc->ProjectionX();
  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading_ue"+ipt,h2mpfleadinv_leading_ue_mc);
  hmc_mpflead_leading_ue->SetStats(0); hmc_mpflead_leading_ue->SetMarkerStyle(1); hmc_mpfleadinv_leading_ue->SetStats(0); hmc_mpfleadinv_leading_ue->SetMarkerStyle(1);
} else {
  inputMCtrig->GetObject(trigDir+"/pmjb_leading"+ipt,pmjb_leading_mc); pmjb_leading_mc->Rebin(rebin);
  inputMCtrig->GetObject(trigDir+"/pmjbinv_leading"+ipt,pmjbinv_leading_mc); pmjbinv_leading_mc->Rebin(rebin);
  inputMCtrig->GetObject(trigDir+"/pmpfinv_leading"+ipt,pmpfinv_leading_mc); pmpfinv_leading_mc->Rebin(rebin);
  inputMCtrig->GetObject(trigDir+"/pmpflead_leading"+ipt,pmpf_leading_mc); pmpf_leading_mc->Rebin(rebin); //testing
  inputMCtrig->GetObject(trigDir+"/pmpflead_leading2"+ipt,pmpf_leading2_mc); pmpf_leading2_mc->Rebin(rebin); //testing

  inputMCtrig->GetObject(trigDir+"/pmpfleadinv_leading"+ipt,pmpfleadinv_leading_mc); pmpfleadinv_leading_mc->Rebin(rebin);
  inputMCtrig->GetObject(trigDir+"/pmpfleadinv_leading2"+ipt,pmpfleadinv_leading2_mc); pmpfleadinv_leading2_mc->Rebin(rebin);
  inputMCtrig->GetObject(trigDir+"/h2mpfleadinv_leading"+ipt,h2mpfleadinv_leading_mc);
  inputMCtrig->GetObject(trigDir+"/h2mpfleadinv_leading2"+ipt,h2mpfleadinv_leading2_mc);
  inputMCtrig->GetObject(trigDir+"/pmpflead_leading_ue"+ipt,pmpflead_leading_ue_mc); pmpflead_leading_ue_mc->Rebin(rebin); hmc_mpflead_leading_ue = pmpflead_leading_ue_mc->ProjectionX();
  inputMCtrig->GetObject(trigDir+"/h2mpflead_leading_ue"+ipt,h2mpflead_leading_ue_mc);
  inputMCtrig->GetObject(trigDir+"/pmpfleadinv_leading_ue"+ipt,pmpfleadinv_leading_ue_mc); pmpfleadinv_leading_ue_mc->Rebin(rebin); hmc_mpfleadinv_leading_ue = pmpfleadinv_leading_ue_mc->ProjectionX();
  inputMCtrig->GetObject(trigDir+"/h2mpfleadinv_leading_ue"+ipt,h2mpfleadinv_leading_ue_mc);
  hmc_mpflead_leading_ue->SetStats(0); hmc_mpflead_leading_ue->SetMarkerStyle(1); hmc_mpfleadinv_leading_ue->SetStats(0); hmc_mpfleadinv_leading_ue->SetMarkerStyle(1);
}

  inputMC->GetObject(mcDir+"/pmjb_ptave"+ipt,pmjb_ptave_mc); pmjb_ptave_mc->Rebin(rebin); hmc_mjb_ptave = pmjb_ptave_mc->ProjectionX(); hmc_mjb_ptave->SetStats(0); hmc_mjb_ptave->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/pmpf_ptave"+ipt,pmpf_ptave_mc); pmpf_ptave_mc->Rebin(rebin); hmc_mpf_ptave = pmpf_ptave_mc->ProjectionX(); hmc_mpf_ptave->SetStats(0); hmc_mpf_ptave->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/pmpf_ptave2"+ipt,pmpf_ptave2_mc); pmpf_ptave2_mc->Rebin(rebin); hmc_mpf_ptave2 = pmpf_ptave2_mc->ProjectionX(); hmc_mpf_ptave2->SetStats(0); hmc_mpf_ptave2->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/h2mjb_ptave"+ipt,h2mjb_ptave_mc); //h2mjb_ptave_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/h2mpf_ptave"+ipt,h2mpf_ptave_mc); //h2mpf_ptave_mc->Rebin(rebin);
  inputMC->GetObject(mcDir+"/h2mpf_ptave2"+ipt,h2mpf_ptave2_mc); //h2mpf_ptave2_mc->Rebin(rebin);
//--Perhaps from inputMCtrig/trigDir (Apr 24, 2020)
//  inputMC->GetObject(mcDir+"/pmpflead_leading_ue",pmpflead_leading_ue_mc); pmpflead_leading_ue_mc->Rebin(rebin); hmc_mpflead_leading_ue = pmpflead_leading_ue_mc->ProjectionX();
//  inputMC->GetObject(mcDir+"/h2mpflead_leading_ue",h2mpflead_leading_ue_mc);
//  hmc_mpflead_leading_ue->SetStats(0); hmc_mpflead_leading_ue->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/pmjbinv_ptave"+ipt,pmjbinv_ptave_mc); pmjbinv_ptave_mc->Rebin(rebin); hmc_mjbinv_ptave = pmjbinv_ptave_mc->ProjectionX(); hmc_mjbinv_ptave->SetStats(0); hmc_mjbinv_ptave->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/pmpfinv_ptave"+ipt,pmpfinv_ptave_mc); pmpfinv_ptave_mc->Rebin(rebin); hmc_mpfinv_ptave = pmpfinv_ptave_mc->ProjectionX(); hmc_mpfinv_ptave->SetStats(0); hmc_mpfinv_ptave->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/pmpfinv_ptave2"+ipt,pmpfinv_ptave2_mc); pmpfinv_ptave2_mc->Rebin(rebin); hmc_mpfinv_ptave2 = pmpfinv_ptave2_mc->ProjectionX(); hmc_mpfinv_ptave2->SetStats(0); hmc_mpfinv_ptave2->SetMarkerStyle(1);
  inputMC->GetObject(mcDir+"/h2mjbinv_ptave"+ipt,h2mjbinv_ptave_mc); //h2mjbinv_ptave->Rebin(rebin);
  inputMC->GetObject(mcDir+"/h2mpfinv_ptave"+ipt,h2mpfinv_ptave_mc); //h2mpfinv_ptave->Rebin(rebin);
  inputMC->GetObject(mcDir+"/h2mpfinv_ptave2"+ipt,h2mpfinv_ptave2_mc); //h2mpfinv_ptave2->Rebin(rebin);
//--Perhaps from inputMCtrig/trigDir (Apr 24, 2020)
//  inputMC->GetObject(mcDir+"/pmpfleadinv_leading",pmpfleadinv_leading_mc); pmpfleadinv_leading_mc->Rebin(rebin);
//  inputMC->GetObject(mcDir+"/pmpfleadinv_leading2",pmpfleadinv_leading2_mc); pmpfleadinv_leading2_mc->Rebin(rebin);
//  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading",h2mpfleadinv_leading_mc);
//  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading2",h2mpfleadinv_leading2_mc);
//  inputMC->GetObject(mcDir+"/pmpfleadinv_leading_ue",pmpfleadinv_leading_ue_mc); pmpfleadinv_leading_ue_mc->Rebin(rebin); hmc_mpfleadinv_leading_ue = pmpfleadinv_leading_ue_mc->ProjectionX();
//  inputMC->GetObject(mcDir+"/h2mpfleadinv_leading_ue",h2mpfleadinv_leading_ue_mc);
//  hmc_mpfleadinv_leading_ue->SetStats(0); hmc_mpfleadinv_leading_ue->SetMarkerStyle(1);
    ////inputMC->GetObject(mcDir+"/pmpflead_ptave_ue",pmpflead_ptave_ue_mc); pmpflead_ptave_ue_mc->Rebin(rebin); hmc_mpflead_ptave_ue = pmpflead_ptave_ue_mc->ProjectionX();
    ////inputMC->GetObject(mcDir+"/h2mpflead_ptave_ue",h2mpflead_ptave_ue_mc);
    ////hmc_mpflead_ptave_ue->SetStats(0); hmc_mpflead_ptave_ue->SetMarkerStyle(1);

//Revisit on Sep 11 2019, but it needs rebinned outputs...
//if(useTrigDir) {
//  inputMCtrig->GetObject(trigDir+"/pmjb_recoil",pmjb_recoil_mc);
//  inputMCtrig->GetObject(trigDir+"/pmpf_recoil",pmpf_recoil_mc);
//  inputMCtrig->GetObject(trigDir+"/pmpf_recoil2",pmpf_recoil2_mc);
//}
if(usehist2) {
  //--updated on Apr 24
  if(ptlead_mcDir) {
    inputMC->GetObject(mcDir+"/h2mjb_leading"+ipt,h2mjb_leading_mc);
    inputMC->GetObject(mcDir+"/h2mjbinv_leading"+ipt,h2mjbinv_leading_mc); //testing
    inputMC->GetObject(mcDir+"/h2mpflead_leading"+ipt,h2mpf_leading_mc); //testing
    inputMC->GetObject(mcDir+"/h2mpflead_leading2"+ipt,h2mpf_leading2_mc); //testing
  } else {
    inputMCtrig->GetObject(trigDir+"/h2mjb_leading"+ipt,h2mjb_leading_mc);
    inputMCtrig->GetObject(trigDir+"/h2mjbinv_leading"+ipt,h2mjbinv_leading_mc); //testing
    inputMCtrig->GetObject(trigDir+"/h2mpflead_leading"+ipt,h2mpf_leading_mc); //testing
    inputMCtrig->GetObject(trigDir+"/h2mpflead_leading2"+ipt,h2mpf_leading2_mc); //testing
  }
}
if(leading_rebin) {
    inputMC->GetObject(mcDir+"/pmjb_leading"+ipt,pmjb_leading_mc); pmjb_leading_mc->Rebin(rebin);
    inputMC->GetObject(mcDir+"/h2mjb_leading"+ipt,h2mjb_leading_mc);

    inputMC->GetObject(mcDir+"/pmjbinv_leading"+ipt,pmjbinv_leading_mc); pmjbinv_leading_mc->Rebin(rebin);
    inputMC->GetObject(mcDir+"/h2mjbinv_leading"+ipt,h2mjbinv_leading_mc);
    inputMC->GetObject(mcDir+"/pmpflead_leading"+ipt,pmpf_leading_mc); pmpf_leading_mc->Rebin(rebin);
    inputMC->GetObject(mcDir+"/pmpflead_leading2"+ipt,pmpf_leading2_mc); pmpf_leading2_mc->Rebin(rebin);
    inputMC->GetObject(mcDir+"/h2mpflead_leading"+ipt,h2mpf_leading_mc);
    inputMC->GetObject(mcDir+"/h2mpflead_leading2"+ipt,h2mpf_leading2_mc);
}

TH1::SetDefaultSumw2();
gStyle->SetErrorX(0.0001);

TH1D *hmjbinv = 0, *hmjbinv_mc = 0;
//TH1D *hdt_mjb0_leading = leading_rebin ? pmjbinv_leading_rebin->ProjectionX("hdt_mjb0_leading") : pmjbinv_leading->ProjectionX("hdt_mjb0_leading");
TH1D *hdt_mjb0_leading = pmjbinv_leading->ProjectionX("hdt_mjb0_leading"), *hmc_mjb0_leading = pmjbinv_leading_mc->ProjectionX("hmc_mjb0_leading");
hdt_mjb0_leading->SetStats(0); hdt_mjb0_leading->SetMarkerStyle(1); hmc_mjb0_leading->SetStats(0); hmc_mjb0_leading->SetMarkerStyle(1);
if(mjbinv) {
  if(leading_rebin) hmjbinv = pmjbinv_leading_rebin->ProjectionX("hmjbinv"); //Oct 22, 2019
  else hmjbinv = pmjbinv_leading->ProjectionX("hmjbinv");  
  hmjbinv_mc = pmjbinv_leading_mc->ProjectionX("hmjbinv_mc");

  TH1D *hmjbx_leading = pmjb_leading->ProjectionX("hmjbx_leading");
  for(int i=1; i!=hmjbinv->GetNbinsX()+1; ++i) {
    int ibin = hdt_mjb_ptave->GetBinLowEdge(i);
    double xdt = hmjbinv->GetBinContent(i), err_xdt = hmjbinv->GetBinError(i);
    double xmc = hmjbinv_mc->GetBinContent(i), err_xmc = hmjbinv_mc->GetBinError(i);
    if(hmjbinv->GetBinContent(i)!=0) { hmjbinv->SetBinContent(i,1/xdt); hmjbinv->SetBinError(i,err_xdt/(xdt*xdt)); } 
    if(hmjbinv_mc->GetBinContent(i)!=0) { hmjbinv_mc->SetBinContent(i,1/xmc); hmjbinv_mc->SetBinError(i,err_xmc/(xmc*xmc)); }
    //--important update for UL17mc and set the correct error value!
    //if(UL17mc) { if(i==1) cout << " *** mjb_ptave and mpf_ptave(2) are inversed and their errors are estimated from \"mjbinv_ptave\" and \"mpfinv_ptave(2)\", respectively (Apr 22, 2020)" << endl; }
    //else { cout << " *** mjbmc = 0 for MC/Pt20, so the current code doesn't work" << endl; return 0; }
    double mjbdt = hdt_mjb_ptave->GetBinContent(i), mjbmc = hmc_mjb_ptave->GetBinContent(i);
    double err0_mjbdt = hdt_mjb_ptave->GetBinError(i), err0_mjbmc = hmc_mjb_ptave->GetBinError(i);
    double err_mjbdt = hdt_mjbinv_ptave->GetBinError(i), err_mjbmc = hmc_mjbinv_ptave->GetBinError(i);
//    if(mjbdt!=0) hdt_mjb_ptave->SetBinContent(i,(1+mjbdt)/(1-mjbdt)); if(mjbmc!=0) hmc_mjb_ptave->SetBinContent(i,(1+mjbmc)/(1-mjbmc));
    if(ibin>=jtpt) hdt_mjb_ptave->SetBinContent(i,(1+mjbdt)/(1-mjbdt)), hmc_mjb_ptave->SetBinContent(i,(1+mjbmc)/(1-mjbmc));
    hdt_mjb_ptave->SetBinError(i,err_mjbdt), hmc_mjb_ptave->SetBinError(i,err_mjbmc);
    ////hdt_mjb_ptave->SetBinContent(i,(1+mjbdt)/(1-mjbdt)); hmc_mjb_ptave->SetBinContent(i,(1+mjbmc)/(1-mjbmc)); // 1 if mjb = 0
    //--test version for Run2017C with ptave (Apr 9, 2020)
//    double mpfdt = hdt_mpfinv_ptave->GetBinContent(i), mpfdt2 = hdt_mpfinv_ptave2->GetBinContent(i);
//    if(mpfdt!=0) hdt_mpfinv_ptave->SetBinContent(i,(1+mpfdt)/(1-mpfdt)); if(mpfdt!=0) hdt_mpfinv_ptave2->SetBinContent(i,(1+mpfdt2)/(1-mpfdt2));
    //--swapped: ptave2 version for Run2017D/E, so hdt_mpf_ptave is to be inversed (i.e., hdt_mpfinv_ptave has been already inversed in the code)
    double mpfdt = hdt_mpf_ptave->GetBinContent(i), mpfdt2 = hdt_mpf_ptave2->GetBinContent(i);
    double err0_mpfdt = hdt_mpf_ptave->GetBinError(i), err0_mpfdt2 = hdt_mpf_ptave2->GetBinError(i);
    double err_mpfdt = hdt_mpfinv_ptave->GetBinError(i), err_mpfdt2 = hdt_mpfinv_ptave2->GetBinError(i);
    if(ibin>=jtpt) hdt_mpf_ptave->SetBinContent(i,(1+mpfdt)/(1-mpfdt)), hdt_mpf_ptave2->SetBinContent(i,(1+mpfdt2)/(1-mpfdt2));
    //if(mpfdt!=0) hdt_mpf_ptave->SetBinContent(i,1-mpfdt); if(mpfdt2!=0) hdt_mpf_ptave2->SetBinContent(i,1-mpfdt2); //test (Apr 23, 2020)
    hdt_mpf_ptave->SetBinError(i,err_mpfdt), hdt_mpf_ptave2->SetBinError(i,err_mpfdt2);

      double mpfmc = hmc_mpf_ptave->GetBinContent(i), mpfmc2 = hmc_mpf_ptave2->GetBinContent(i);
      double err0_mpfmc = hmc_mpf_ptave->GetBinError(i), err0_mpfmc2 = hmc_mpf_ptave2->GetBinError(i);
      double err_mpfmc = hmc_mpfinv_ptave->GetBinError(i), err_mpfmc2 = hmc_mpfinv_ptave2->GetBinError(i);
      if(ibin>=jtpt) hmc_mpf_ptave->SetBinContent(i,(1+mpfmc)/(1-mpfmc)), hmc_mpf_ptave2->SetBinContent(i,(1+mpfmc2)/(1-mpfmc2));
      //if(mpfmc!=0) hmc_mpf_ptave->SetBinContent(i,1-mpfmc); if(mpfmc2!=0) hmc_mpf_ptave2->SetBinContent(i,1-mpfmc2); //test (Apr 23, 2020)
      hmc_mpf_ptave->SetBinError(i,err_mpfmc); hmc_mpf_ptave2->SetBinError(i,err_mpfmc2);

    double xinv = hdt_mjb_ptave->GetBinContent(i), err_xinv = hdt_mjb_ptave->GetBinError(i);
//    printf("    i = %2i (%4i- GeV):  x = %8.3f and   1/x  = %8.3f (%8.3f)\n",i,ibin,mjbdt,xinv,hdt_mjbinv_ptave->GetBinContent(i));
//    printf(" -> i = %2i (%4i- GeV): ex = %8.3f and e(1/x) = %8.3f (%8.3f)\n",i,ibin,err0_mjbdt,err_xinv,hdt_mjbinv_ptave->GetBinError(i));
    //double xinv = hdt_mpf_ptave->GetBinContent(i), err_xinv = hdt_mpf_ptave->GetBinError(i);
    //printf("    i = %2i (%4i- GeV):  x = %8.3f and   1/x  = %8.3f (%8.3f)\n",i,ibin,mpfdt,xinv,hdt_mpfinv_ptave->GetBinContent(i));
    //printf(" ->                   : ex = %8.3f and e(1/x) = %8.3f (%8.3f)\n",err0_mjbdt,err_xinv,hdt_mpfinv_ptave->GetBinError(i));
  }
}
TH1D *hmpfinv = 0, *hmpfinv_mc = 0;
if(mpfinv) {
  hmpfinv = pmpfinv_leading->ProjectionX("hmpfinv");
  hmpfinv_mc = pmpfinv_leading_mc->ProjectionX("hmpfinv_mc");
  for(int i=1; i!=hmpfinv->GetNbinsX()+1; ++i) {
    if(hmpfinv->GetBinContent(i)!=0) hmpfinv->SetBinContent(i,1/(double)hmpfinv->GetBinContent(i));
    if(hmpfinv_mc->GetBinContent(i)!=0) hmpfinv_mc->SetBinContent(i,1/(double)hmpfinv_mc->GetBinContent(i));
  }
}

TProfile *hmjb_leading_pfx = 0, *hmjb_leading_pfx_mc = 0;
TProfile *hmjbinv_leading_pfx = 0, *hmjbinv_leading_pfx_mc = 0, *hmpf_leading_pfx = 0, *hmpf_leading2_pfx = 0, *hmpf_leading_pfx_mc = 0, *hmpf_leading2_pfx_mc = 0;
TH1D *hmjb_leading = 0, *hmjb_leading_mc = 0, *hmjbinv_leading = 0, *hmjbinv_leading_mc = 0, *hmpf_leading = 0, *hmpf_leading2 = 0, *hmpf_leading_mc = 0, *hmpf_leading2_mc = 0;
//--added after discussion with Mikko about error calculation in global fit (Sep 21)
TProfile *hmjb_recoil_pfx = 0, *hmjb_recoil_pfx_mc = 0; TH1D *hmjb_recoil = 0, *hmjb_recoil_mc = 0;
TProfile *hmpf_recoil_pfx = 0, *hmpf_recoil_pfx_mc = 0; TH1D *hmpf_recoil = 0, *hmpf_recoil_mc = 0;
TProfile *hmpf_recoil2_pfx = 0, *hmpf_recoil2_pfx_mc = 0; TH1D *hmpf_recoil2 = 0, *hmpf_recoil2_mc = 0;

if(usehist2) {
  TCanvas *c0 = new TCanvas("c0","c0",300,300,600,600);
  c0->Divide(2,2);
  c0->cd(1); h2mjbinv_leading->Draw("zcol");
  hmjb_leading_pfx    = h2mjb_leading->ProfileX("hmjb_leading_pfx");                hmjb_leading_pfx->Rebin(rebin);       hmjb_leading       = hmjb_leading_pfx->ProjectionX("hmjb_leading");
  hmjb_leading_pfx_mc = h2mjb_leading_mc->ProfileX("hmjb_leading_pfx_mc");          hmjb_leading_pfx_mc->Rebin(rebin);    hmjb_leading_mc    = hmjb_leading_pfx_mc->ProjectionX("hmjb_leading_mc");
  hmjbinv_leading_pfx    = h2mjbinv_leading->ProfileX("hmjbinv_leading_pfx");       hmjbinv_leading_pfx->Rebin(rebin);    hmjbinv_leading    = hmjbinv_leading_pfx->ProjectionX("hmjbinv_leading");
  hmjbinv_leading_pfx_mc = h2mjbinv_leading_mc->ProfileX("hmjbinv_leading_pfx_mc"); hmjbinv_leading_pfx_mc->Rebin(rebin); hmjbinv_leading_mc = hmjbinv_leading_pfx_mc->ProjectionX("hmjbinv_leading_mc");
  hmpf_leading_pfx       = h2mpf_leading->ProfileX("hmpf_leading_pfx");             hmpf_leading_pfx->Rebin(rebin);       hmpf_leading       = hmpf_leading_pfx->ProjectionX("hmpf_leading");
  hmpf_leading2_pfx      = h2mpf_leading2->ProfileX("hmpf_leading2_pfx");           hmpf_leading2_pfx->Rebin(rebin);      hmpf_leading2      = hmpf_leading2_pfx->ProjectionX("hmpf_leading2");
  hmpf_leading_pfx_mc    = h2mpf_leading_mc->ProfileX("hmpf_leading_pfx_mc");       hmpf_leading_pfx_mc->Rebin(rebin);    hmpf_leading_mc    = hmpf_leading_pfx_mc->ProjectionX("hmpf_leading_mc");
  hmpf_leading2_pfx_mc   = h2mpf_leading2_mc->ProfileX("hmpf_leading2_pfx_mc");     hmpf_leading2_pfx_mc->Rebin(rebin);   hmpf_leading2_mc   = hmpf_leading2_pfx_mc->ProjectionX("hmpf_leading2_mc");
  c0->cd(2); hmjbinv_leading_pfx->Draw("hist"); pmjbinv_leading->Draw("histsame"); hmjbinv_leading_pfx->SetLineColor(2);

  for(int i=1; i!=hmjbinv_leading->GetNbinsX()+1; ++i) {
    if(hmjbinv_leading->GetBinContent(i)!=0) hmjbinv_leading->SetBinContent(i,1/(double)hmjbinv_leading->GetBinContent(i));
    //double x = hmjbinv_leading_pfx->GetBinContent(i), xinv = (x!=0) ? 1/x : 0; hmjbinv_leading->SetBinContent(i,xinv);
    //printf("i = %2i: x = %8.3f and 1/x = %8.3f (%8.3f)\n",i,x,xinv,hmjbinv_leading->GetBinContent(i));
    if(hmjbinv_leading_mc->GetBinContent(i)!=0) hmjbinv_leading_mc->SetBinContent(i,1/(double)hmjbinv_leading_mc->GetBinContent(i));
  }

  c0->cd(3); hmjbinv_leading->Draw("hist"); hmjbinv_leading->SetLineColor(2); if(mjbinv) hmjbinv->Draw("histsame");

  c0->cd(4);
  hmjb_recoil_pfx    = h2mjb_recoil_rebin->ProfileX("hmjb_recoil_pfx");    hmjb_recoil_pfx->Rebin(rebin);     hmjb_recoil       = hmjb_recoil_pfx->ProjectionX("hmjb_recoil");
  hmjb_recoil_pfx_mc = h2mjb_recoil_mc->ProfileX("hmjb_recoil_pfx_mc");    hmjb_recoil_pfx_mc->Rebin(rebin);  hmjb_recoil_mc    = hmjb_recoil_pfx_mc->ProjectionX("hmjb_recoil_mc");
  hmjb_recoil->Draw("hist"); pmjb_recoil->Draw("sames"); hmjb_recoil->SetLineColor(2); hmjb_recoil_mc->SetLineColor(4);
  //hmjb_recoil_mc->Draw("hist"); pmjb_recoil_mc->Darw("sames");
  //hmjb_recoil->Draw("hist"); hmjb_recoil_mc->Draw("histsame");
  hmpf_recoil_pfx     = h2mpf_recoil_rebin->ProfileX("hmpf_recoil_pfx");   hmpf_recoil_pfx->Rebin(rebin);     hmpf_recoil       = hmpf_recoil_pfx->ProjectionX("hmpf_recoil");
  hmpf_recoil_pfx_mc  = h2mpf_recoil_mc->ProfileX("hmpf_recoil_pfx_mc");   hmpf_recoil_pfx_mc->Rebin(rebin);  hmpf_recoil_mc    = hmpf_recoil_pfx_mc->ProjectionX("hmpf_recoil_mc");
  hmpf_recoil2_pfx    = h2mpf_recoil2_rebin->ProfileX("hmpf_recoil2_pfx"); hmpf_recoil2_pfx->Rebin(rebin);    hmpf_recoil2      = hmpf_recoil2_pfx->ProjectionX("hmpf_recoil2");
  hmpf_recoil2_pfx_mc = h2mpf_recoil2_mc->ProfileX("hmpf_recoil2_pfx_mc"); hmpf_recoil2_pfx_mc->Rebin(rebin); hmpf_recoil2_mc   = hmpf_recoil2_pfx_mc->ProjectionX("hmpf_recoil2_mc");
  //hmpf_recoil->Draw("hist"); pmjb_recoil->Darw("sames"); hmpf_recoil->SetLineColor(2); hmpf_recoil_mc->SetLineColor(4);

  float ymin = 0.6; ymin = 0.95;
  hmjbinv_leading->SetMinimum(ymin); hmjbinv_leading_mc->SetMinimum(ymin); hmpf_leading_mc->SetMinimum(ymin); if(usehist2) hmpf_leading2_mc->SetMinimum(ymin);

  pmjb_leading->SetLineColor(4);
  TH1 *hcomp_mjb_leading = pmjb_leading->ProjectionX(); hcomp_mjb_leading->Divide(hmjb_leading);
  hcomp_mjb_leading->SetLineColor(4); hcomp_mjb_leading->SetMarkerColor(4); hcomp_mjb_leading->SetMarkerStyle(2);

  TCanvas *c1 = new TCanvas("c1","c1",350,350,600,600);
  c1->Divide(2,2);
  c1->cd(1); hmjb_leading->Draw("hist"); pmjb_leading->DrawCopy("histsame"); hcomp_mjb_leading->Draw("esame");
  c1->cd(2); hmjb_leading_mc->Draw("hist"); pmjb_leading_mc->Draw("histsame");
  c1->cd(3); hmjbinv_leading->Draw("hist"); if(mjbinv) hmjbinv->Draw("histsame");
  c1->cd(4); hmjbinv_leading_mc->Draw("hist"); hmjbinv_leading->SetLineColor(2); if(mjbinv) hmjbinv_mc->Draw("histsame");

  TCanvas *c10 = new TCanvas("c10","c10",400,400,600,600);
  c10->Divide(2,2);
  c10->cd(1); hmpf_leading->Draw("hist"); pmpf_leading->DrawCopy("histsame"); hmpf_leading->SetLineColor(2);
  c10->cd(2); hmpf_leading2->Draw("hist"); pmpf_leading2->DrawCopy("histsame"); hmpf_leading2->SetLineColor(2);
  c10->cd(3); hmpf_leading_mc->Draw("hist"); pmpf_leading_mc->DrawCopy("histsame"); hmpf_leading_mc->SetLineColor(2);
  c10->cd(4); hmpf_leading2_mc->Draw("hist"); pmpf_leading2_mc->DrawCopy("histsame"); hmpf_leading2_mc->SetLineColor(2);

}

TCanvas* cvv = new TCanvas("cvv","cvv",0,0,600,600); gPad->SetLogx();
//pmjb_recoil->Draw();
pmjb_recoil->GetYaxis()->SetRangeUser(hmin,hmax); pmjb_recoil->GetXaxis()->SetRangeUser(xmin,xmax);
pmjb_recoil->GetYaxis()->SetTitleOffset(1.1);
pmjb_recoil->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,recoil}");
//pmjb_recoil->GetYaxis()->SetLabelSize(0.03); pmjb_recoil->GetXaxis()->SetLabelSize(0.03);
pmjb_recoil->SetStats(0); pmjb_leading->SetStats(0); pmpf_recoil->SetStats(0); pmpf_leading->SetStats(0);

pmjb_recoil_rebin->GetYaxis()->SetRangeUser(hmin,hmax); pmjb_recoil_rebin->GetXaxis()->SetRangeUser(xmin,xmax);
pmjb_recoil_rebin->GetYaxis()->SetTitleOffset(1.1);
pmjb_recoil_rebin->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,recoil}");
pmjb_recoil_rebin->GetYaxis()->SetLabelSize(0.03); pmjb_recoil_rebin->GetXaxis()->SetLabelSize(0.03);
pmjb_recoil_rebin->SetStats(0); pmjb_leading_rebin->SetStats(0); pmpf_recoil_rebin->SetStats(0); pmpf_leading_rebin->SetStats(0);

pmjb_recoil_mc->GetYaxis()->SetRangeUser(hmin,hmax); pmjb_recoil_mc->GetXaxis()->SetRangeUser(xmin,xmax);
pmjb_recoil_mc->GetYaxis()->SetTitleOffset(1.1);
pmjb_recoil_mc->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,recoil}");
pmjb_recoil_mc->GetYaxis()->SetLabelSize(0.03); pmjb_recoil_mc->GetXaxis()->SetLabelSize(0.03);
pmjb_recoil_mc->SetStats(0); pmjb_leading_mc->SetStats(0); pmpf_recoil_mc->SetStats(0); pmpf_leading_mc->SetStats(0);

pmjb_leading_mc->GetYaxis()->SetRangeUser(hmin,hmax); pmjb_leading_mc->GetXaxis()->SetRangeUser(xmin,xmax);
pmjb_leading_mc->GetYaxis()->SetTitleOffset(1.1);
pmjb_leading_mc->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,leading}");
pmjb_leading_mc->GetYaxis()->SetLabelSize(0.03); pmjb_leading_mc->GetXaxis()->SetLabelSize(0.03);
pmjb_leading_mc->SetStats(0); pmjb_leading_mc->SetStats(0); pmpf_leading_mc->SetStats(0); pmpf_leading_mc->SetStats(0);

hmjb_recoil->GetYaxis()->SetRangeUser(hmin,hmax); hmjb_recoil->GetXaxis()->SetRangeUser(xmin,xmax);
hmjb_recoil->GetYaxis()->SetTitleOffset(1.1);
hmjb_recoil->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,recoil}");
hmjb_recoil->GetYaxis()->SetLabelSize(0.03); hmjb_recoil->GetXaxis()->SetLabelSize(0.03);
hmjb_recoil->SetStats(0); pmjb_leading->SetStats(0); pmpf_recoil->SetStats(0); pmpf_leading->SetStats(0);

hmjb_recoil_mc->GetYaxis()->SetRangeUser(hmin,hmax); hmjb_recoil_mc->GetXaxis()->SetRangeUser(xmin,xmax);
hmjb_recoil_mc->GetYaxis()->SetTitleOffset(1.1);
hmjb_recoil_mc->GetYaxis()->SetTitle("R = p_{T,leading}/p_{T,recoil}");
hmjb_recoil_mc->GetYaxis()->SetLabelSize(0.03); hmjb_recoil_mc->GetXaxis()->SetLabelSize(0.03);
hmjb_recoil_mc->SetStats(0); pmjb_leading_mc->SetStats(0); pmpf_recoil_mc->SetStats(0); pmpf_leading_mc->SetStats(0);

gPad->SetTicks();

  pmjb_recoil->SetMarkerStyle(25); pmjb_recoil->SetLineColor(1); pmjb_recoil->SetMarkerColor(1); pmjb_recoil->Draw("e1");
  pmjb_leading->SetMarkerStyle(24); pmjb_leading->SetLineColor(1); pmjb_leading->SetMarkerColor(1); pmjb_leading->Draw("e1same");
  pmpf_leading->SetMarkerStyle(20); pmpf_leading->SetLineColor(kRed+1); pmpf_leading->SetMarkerColor(kRed+1); pmpf_leading->Draw("e1same");
  pmpf_recoil->SetMarkerStyle(21); pmpf_recoil->SetLineColor(kOrange+8); pmpf_recoil->SetMarkerColor(kOrange+8); pmpf_recoil->Draw("e1same");

  TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
  leg->SetBorderSize(0); leg->SetTextSize(0.03); leg->SetLineColor(1);
  leg->SetLineStyle(1); leg->SetLineWidth(1); leg->SetFillColor(0); leg->SetFillStyle(1001);
  if (both or recoil) {
    //leg->AddEntry(pmjb_recoil,"MJB, p_{T, recoil}","P");
    //leg->AddEntry(pmjb_leading,"MJB, p_{T, leading}","P");
    //leg->AddEntry(pmpf_recoil,"MPF, p_{T, recoil}","P");
    //leg->AddEntry(pmpf_leading,"MPF, p_{T, leading}","P");
    leg->AddEntry(pmjb_recoil,"MJB, p_{T, recoil}","P");
    if(drawMC) leg->AddEntry(pmjb_recoil_mc,str_mc,"P");
    leg->AddEntry(pmpf_recoil,"MPF, p_{T, recoil}","P");
    if(drawMC) leg->AddEntry(pmpf_recoil_mc,str_mc2,"P");
  }
  //leg->Draw();

  //double leg_xmin2 = 0.72, leg_xmax2 = 0.85;
  TLegend *leg2 = new TLegend(leg_xmin+0.2,leg_ymin,leg_xmax+0.15,leg_ymax);
  leg2->SetBorderSize(0); leg2->SetTextSize(0.03); leg2->SetLineColor(1);
  leg2->SetLineStyle(1); leg2->SetLineWidth(1); leg2->SetFillColor(0); leg2->SetFillStyle(1001);

  TCanvas* cvvOld = new TCanvas("cvvOld","cvvOld",100,100,600,600); gPad->SetLogx();

  if( ! jt.Contains("jt") and ! jt.Contains("mc")) {
    pmjb_recoil_rebin->SetMarkerStyle(25); pmjb_recoil_rebin->SetLineColor(1); pmjb_recoil_rebin->SetMarkerColor(1); pmjb_recoil_rebin->Draw("e1");
    pmjb_leading->SetMarkerStyle(24); pmjb_leading->SetLineColor(1); pmjb_leading->SetMarkerColor(1); pmjb_leading->Draw("e1same");
    pmpf_leading->SetMarkerStyle(20); pmpf_leading->SetLineColor(kRed+1); pmpf_leading->SetMarkerColor(kRed+1); pmpf_leading->Draw("e1same");
    pmpf_recoil_rebin->SetMarkerStyle(21); pmpf_recoil_rebin->SetLineColor(kOrange+8); pmpf_recoil_rebin->SetMarkerColor(kOrange+8); pmpf_recoil_rebin->Draw("e1same");
    leg->Draw();

    pmjb_recoil_mc->SetMarkerStyle(26); pmjb_recoil_mc->SetLineColor(kCyan+1); pmjb_recoil_mc->SetMarkerColor(kCyan+1); pmjb_recoil_mc->Draw("e1same");
    pmpf_recoil_mc->SetMarkerStyle(22); pmpf_recoil_mc->SetLineColor(kOrange+1); pmpf_recoil_mc->SetMarkerColor(kOrange+1); pmpf_recoil_mc->Draw("e1same");
    pmjb_leading_mc->SetMarkerStyle(27); pmjb_leading_mc->SetLineColor(kCyan+1); pmjb_leading_mc->SetMarkerColor(kCyan+1); pmjb_leading_mc->Draw("e1same");
    pmpf_leading_mc->SetMarkerStyle(33); pmpf_leading_mc->SetLineColor(kOrange+1); pmpf_leading_mc->SetMarkerColor(kOrange+1); pmpf_leading_mc->Draw("e1same");

    if(both or leading) {
      //leg2->AddEntry(pmjb_recoil_mc,str_mc,"P");
      //leg2->AddEntry(pmjb_leading_mc,str_mc,"P");
      //leg2->AddEntry(pmpf_recoil_mc,str_mc,"P");
      //leg2->AddEntry(pmpf_leading_mc,str_mc,"P");
      leg2->AddEntry(pmjb_leading,"MJB, p_{T, leading}","P");
      if(drawMC) leg2->AddEntry(pmjb_leading_mc,str_mc,"P");
      leg2->AddEntry(pmpf_leading,"MPF, p_{T, leading}","P");
      if(drawMC) leg2->AddEntry(pmpf_leading_mc,str_mc2,"P");
      leg2->Draw();
    }
  }

  TCanvas* cvvNew = new TCanvas("cvvNew","cvvNew",150,150,600,600); gPad->SetLogx();
  pmjb_recoil_rebin->SetMarkerStyle(25); pmjb_recoil_rebin->SetLineColor(1); pmjb_recoil_rebin->SetMarkerColor(1); pmjb_recoil_rebin->Draw("e1");
  if(mjbinv) { hmjbinv->SetMarkerStyle(24); hmjbinv->SetLineColor(1); hmjbinv->SetMarkerColor(1); hmjbinv->Draw("e1same"); }
  pmpf_leading->SetMarkerStyle(20); pmpf_leading->SetLineColor(kRed+1); pmpf_leading->SetMarkerColor(kRed+1); pmpf_leading->Draw("e1same");
  pmpf_recoil_rebin->SetMarkerStyle(21); pmpf_recoil_rebin->SetLineColor(kOrange+8); pmpf_recoil_rebin->SetMarkerColor(kOrange+8); pmpf_recoil_rebin->Draw("e1same");

  pmjb_recoil_mc->SetMarkerStyle(26); pmjb_recoil_mc->SetLineColor(kCyan+1); pmjb_recoil_mc->SetMarkerColor(kCyan+1); pmjb_recoil_mc->Draw("e1same");
  pmpf_recoil_mc->SetMarkerStyle(22); pmpf_recoil_mc->SetLineColor(kOrange+1); pmpf_recoil_mc->SetMarkerColor(kOrange+1); pmpf_recoil_mc->Draw("e1same");
  if(mjbinv) { hmjbinv_mc->SetMarkerStyle(27); hmjbinv_mc->SetLineColor(kCyan+1); hmjbinv_mc->SetMarkerColor(kCyan+1); hmjbinv_mc->Draw("e1same"); }
  pmpf_leading_mc->SetMarkerStyle(33); pmpf_leading_mc->SetLineColor(kOrange+1); pmpf_leading_mc->SetMarkerColor(kOrange+1); pmpf_leading_mc->Draw("e1same");

//return;
  TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.030); tex->SetTextColor(kBlack);
  tex->DrawLatex(0.20,0.80,"|#eta_{leading}| < 1.3");
  tex->SetTextSize(0.04); tex->DrawLatex(0.11,0.91,"CMS");
  tex->SetTextFont(42); tex->SetTextSize(0.038);
  tex->DrawLatex(0.6,0.91,str_dt+" (13 TeV)");
  //if(str_dt=="2018ABCD") tex->DrawLatex(0.5,0.91,str_dt+" (13 TeV)"); else tex->DrawLatex(0.6,0.91,str_dt+" (13 TeV)");


  TCanvas* cvv2 = new TCanvas("cvv2","cvv2",100,100,600,600); gPad->SetTicks();
  pmjb_recoil_mc->SetMarkerStyle(26); pmjb_recoil_mc->SetLineColor(kCyan+1); pmjb_recoil_mc->SetMarkerColor(kCyan+1); pmjb_recoil_mc->Draw("e1");
  //pmpf_recoil_mc->SetMarkerStyle(2); pmpf_recoil_mc->SetLineColor(kOrange+8); pmpf_recoil_mc->SetMarkerColor(kOrange+8); pmpf_recoil_mc->Draw("e1same");
  if (dtfile.Contains("MC")) {
    pmjb_recoil->SetMarkerStyle(25); pmjb_recoil->SetLineColor(kCyan+1); pmjb_recoil->SetMarkerColor(kCyan+1); pmjb_recoil->Draw("histsame");
    pmpf_recoil->SetMarkerStyle(21); pmpf_recoil->SetLineColor(kOrange+8); pmpf_recoil->SetMarkerColor(kOrange+8); pmpf_recoil->Draw("histsame");
  }
  pmjb_recoil_mc->GetXaxis()->SetRangeUser(100,3600); pmjb_recoil_mc->GetYaxis()->SetRangeUser(0.5,1.3); gPad->SetLogx();
  //leg2->Draw();

  TCanvas* cvv3 = new TCanvas("cvv3","cvv3",150,150,600,600); gPad->SetTicks();
  pmjb_recoil_rebin->SetMarkerStyle(25); pmjb_recoil_rebin->SetLineColor(kCyan+1); pmjb_recoil_rebin->SetMarkerColor(kCyan+1); pmjb_recoil_rebin->Draw("e1");
  pmpf_recoil_rebin->SetMarkerStyle(21); pmpf_recoil_rebin->SetLineColor(kOrange+8); pmpf_recoil_rebin->SetMarkerColor(kOrange+8); pmpf_recoil_rebin->Draw("e1same");
  //leg2->Draw();

  TCanvas* cvv4 = new TCanvas("cvv4","cvv4",200,200,600,600); gPad->SetTicks();
  //pmpf_recoil_mc->SetMarkerStyle(2); pmpf_recoil_mc->SetLineColor(kOrange+8); pmpf_recoil_mc->SetMarkerColor(kOrange+8); pmpf_recoil_mc->Draw("e1");

  pmjbinv_leading->Draw(""); pmjbinv_leading->SetMarkerStyle(24); pmjbinv_leading->SetLineColor(kCyan+3); pmjbinv_leading->SetMarkerColor(kCyan+3);

  TCanvas* cvv5 = new TCanvas("cvv5","cvv5",250,250,600,600); gPad->SetTicks();
  if(mjbinv) {
    hmjbinv->Draw(""); hmjbinv->SetMarkerStyle(24); hmjbinv->SetLineColor(kCyan+3); hmjbinv->SetMarkerColor(kCyan+3);
  }

  pmpf_recoil_mc->SetMinimum(0.5);

  TCanvas* cvv6 = new TCanvas("cvv6","cvv6",300,300,600,600); gPad->SetTicks();
  pmjbinv_leading->Draw(); pmjbinv_leading_rebin->Draw("histsames");

  //TH1D *hdt_mjb_recoil = pmjb_recoil->ProjectionX();
  TH1D *hdt_mjb_recoil = pmjb_recoil_rebin->ProjectionX("hdt_mjb_recoil"); //use rebin!
  TH1D *hmc_mjb_recoil = pmjb_recoil_mc->ProjectionX("hmc_mjb_recoil");
  TH1D *hdt_mjb_leading = pmjb_leading->ProjectionX("hdt_mjb_leading"), *hdt_mjbx_leading = pmjb_leading->ProjectionX("hdt_mjbx_leading");
  TH1D *hmc_mjb_leading = pmjb_leading_mc->ProjectionX("hmc_mjb_leading"), *hmc_mjbx_leading = pmjb_leading_mc->ProjectionX("hmc_mjbx_leading");
  hdt_mjbx_leading->SetStats(0); hdt_mjbx_leading->SetMarkerStyle(1); hmc_mjbx_leading->SetStats(0); hmc_mjbx_leading->SetMarkerStyle(1);
  if(mjbinv) {
    hdt_mjb_leading->Reset(); hmc_mjb_leading->Reset();
    hdt_mjb_leading = (TH1D*)hmjbinv->Clone(); hmc_mjb_leading = (TH1D*)hmjbinv_mc->Clone();
    mjbinvDraw = true;
  }
  if(usehist2Draw) {
    hdt_mjb_leading->Reset(); hmc_mjb_leading->Reset();
    hdt_mjb_leading = (TH1D*)hmjbinv_leading->Clone(); hmc_mjb_leading = (TH1D*)hmjbinv_leading_mc->Clone();
    if( !mjbinv ) { hdt_mjb_leading = (TH1D*)hmjb_leading->Clone(); hmc_mjb_leading = (TH1D*)hmjb_leading_mc->Clone(); }
    //MJB recoil updated on Sep 21, 2019, resulting in slightly different distributions from the previous ones, multijet_Rebin2_20190920
    hdt_mjb_recoil->Reset(); hmc_mjb_recoil->Reset();
    hdt_mjb_recoil = (TH1D*)hmjb_recoil->Clone(); hmc_mjb_recoil = (TH1D*)hmjb_recoil_mc->Clone();
  }

  TH1D *hdt_mpf_recoil = pmpf_recoil_rebin->ProjectionX("hdt_mpf_recoil"); //use rebin!
  TH1D *hmc_mpf_recoil = pmpf_recoil_mc->ProjectionX("hmc_mpf_recoil");
  TH1D *hdt_mpf_leading = pmpf_leading->ProjectionX("hdt_mpf_leading");
  TH1D *hmc_mpf_leading = pmpf_leading_mc->ProjectionX("hmc_mpf_leading");
  TH1D *hdt_mpfleadinv_leading = 0, *hmc_mpfleadinv_leading = 0;

  hdt_mpfleadinv_leading = pmpfleadinv_leading->ProjectionX("hdt_mpfleadinv_leading"); hdt_mpfleadinv_leading->SetStats(0); hdt_mpfleadinv_leading->SetMarkerStyle(1); 
  hmc_mpfleadinv_leading = pmpfleadinv_leading_mc->ProjectionX("hmc_mpfleadinv_leading"); hmc_mpfleadinv_leading->SetStats(0); hmc_mpfleadinv_leading->SetMarkerStyle(1);

  if (met2) {
    hdt_mpf_recoil->Reset();  hdt_mpf_recoil = pmpf_recoil2_rebin->ProjectionX("hdt_mpf_recoil"); //use rebin!
    hmc_mpf_recoil->Reset();  hmc_mpf_recoil = pmpf_recoil2_mc->ProjectionX("hmc_mpf_recoil");
    hdt_mpf_leading->Reset(); hdt_mpf_leading = pmpf_leading2->ProjectionX("hdt_mpf_leading");
    hmc_mpf_leading->Reset(); hmc_mpf_leading = pmpf_leading2_mc->ProjectionX("hmc_mpf_leading");

    hdt_mpfleadinv_leading->Reset(); hdt_mpfleadinv_leading = pmpfleadinv_leading2->ProjectionX("hdt_mpfleadinv_leading");
    hmc_mpfleadinv_leading->Reset(); hmc_mpfleadinv_leading = pmpfleadinv_leading2_mc->ProjectionX("hmc_mpfleadinv_leading");

  }
//  if(usehist2Draw) {
  if(usehist2Draw and ! mpfprofDraw) {
    hdt_mpf_leading->Reset(); hmc_mpf_leading->Reset();
    hdt_mpf_recoil->Reset(); hmc_mpf_recoil->Reset();
    if(met2) {
      hdt_mpf_leading = (TH1D*)hmpf_leading2->Clone(); hmc_mpf_leading = (TH1D*)hmpf_leading2_mc->Clone();
      hdt_mpf_recoil = (TH1D*)hmpf_recoil2->Clone(); hmc_mpf_recoil = (TH1D*)hmpf_recoil2_mc->Clone();
    } else {
      hdt_mpf_leading = (TH1D*)hmpf_leading->Clone(); hmc_mpf_leading = (TH1D*)hmpf_leading_mc->Clone();
      hdt_mpf_recoil = (TH1D*)hmpf_recoil->Clone(); hmc_mpf_recoil = (TH1D*)hmpf_recoil_mc->Clone();
    }
  }
  //--no more used
  if(mpfinv) {
    hdt_mpf_leading->Reset(); hmc_mpf_leading->Reset();
    hdt_mpf_leading = (TH1D*)hmpfinv->Clone(); hmc_mpf_leading = (TH1D*)hmpfinv_mc->Clone();
    mpfinvDraw = true;
  }

  hdt_mjb_recoil->GetYaxis()->SetTitle("#LTp_{T,leading}/p_{T,recoil}#GT");
  if(mjbinvDraw) hdt_mjb_recoil->GetYaxis()->SetTitle("#LTp_{T,leading}#GT / #LTp_{T,recoil}#GT");
  hdt_mjb_recoil->GetYaxis()->SetRangeUser(hmin,hmax); hdt_mjb_recoil->GetXaxis()->SetRangeUser(xmin,xmax);
  hmc_mjb_recoil->GetYaxis()->SetTitle("Data / MC"); hmc_mjb_recoil->GetXaxis()->SetRangeUser(xmin,xmax); hmc_mjb_recoil->GetYaxis()->SetRangeUser(hminRat,hmaxRat);
  hdt_mjb_recoil->SetMarkerStyle(25); hdt_mjb_recoil->SetLineColor(1); hdt_mjb_recoil->SetMarkerColor(1);
  hmc_mjb_recoil->SetMarkerStyle(26); hmc_mjb_recoil->SetLineColor(kCyan+1); hmc_mjb_recoil->SetMarkerColor(kCyan+1);

  hdt_mjb_leading->GetYaxis()->SetTitle("#LTp_{T,leading}/p_{T,recoil}#GT"); hdt_mjb_leading->GetYaxis()->SetRangeUser(hmin,hmax); hdt_mjb_leading->GetXaxis()->SetRangeUser(xmin,xmax);
  hmc_mjb_leading->GetYaxis()->SetTitle("Data / MC"); hmc_mjb_leading->GetXaxis()->SetRangeUser(xmin,xmax); //hmc_mjb_leading->GetYaxis()->SetRangeUser(hminRat,hmaxRat);
  hdt_mjb_leading->SetMarkerStyle(24); hdt_mjb_leading->SetLineColor(1); hdt_mjb_leading->SetMarkerColor(1);
  hmc_mjb_leading->SetMarkerStyle(27); hmc_mjb_leading->SetLineColor(kCyan+1); hmc_mjb_leading->SetMarkerColor(kCyan+1);

  hdt_mpf_recoil->GetYaxis()->SetTitle("#LTp_{T,leading}/p_{T,recoil}#GT"); hdt_mpf_recoil->GetYaxis()->SetRangeUser(hmin,hmax); hdt_mpf_recoil->GetXaxis()->SetRangeUser(xmin,xmax);
  hmc_mpf_recoil->GetYaxis()->SetTitle("Data / MC"); hmc_mpf_recoil->GetXaxis()->SetRangeUser(xmin,xmax); hmc_mpf_recoil->GetYaxis()->SetRangeUser(hminRat,hmaxRat);
  hdt_mpf_recoil->SetMarkerStyle(21); hdt_mpf_recoil->SetLineColor(kOrange+8); hdt_mpf_recoil->SetMarkerColor(kOrange+8);
  hmc_mpf_recoil->SetMarkerStyle(22); hmc_mpf_recoil->SetLineColor(kOrange+1); hmc_mpf_recoil->SetMarkerColor(kOrange+1);

  hdt_mpf_leading->GetYaxis()->SetTitle("#LTp_{T,leading}/p_{T,recoil}#GT"); hdt_mpf_leading->GetYaxis()->SetRangeUser(hmin,hmax); hdt_mpf_leading->GetXaxis()->SetRangeUser(xmin,xmax);
  hmc_mpf_leading->GetYaxis()->SetTitle("Data / MC"); hmc_mpf_leading->GetXaxis()->SetRangeUser(xmin,xmax); //hmc_mpf_leading->GetYaxis()->SetRangeUser(hminRat,hmaxRat);
  hdt_mpf_leading->SetMarkerStyle(20); hdt_mpf_leading->SetLineColor(kRed+1); hdt_mpf_leading->SetMarkerColor(kRed+1);
  hmc_mpf_leading->SetMarkerStyle(33); hmc_mpf_leading->SetLineColor(kOrange+1); hmc_mpf_leading->SetMarkerColor(kOrange+1);

  hdt_mjb_leading->SetMarkerSize(0.8); hdt_mpf_leading->SetMarkerSize(0.8);
  hdt_mjb_recoil->SetMarkerSize(0.8); hdt_mpf_recoil->SetMarkerSize(0.8);

  hmc_mjb_recoil->GetXaxis()->SetTitle("p_{T, recoil} or p_{T, leading} (GeV)");
  if (recoil) hmc_mjb_recoil->GetXaxis()->SetTitle("p_{T, recoil} (GeV)"); else if (leading) hmc_mjb_recoil->GetXaxis()->SetTitle("p_{T, leading} (GeV)");
  hmc_mjb_recoil->GetYaxis()->SetNdivisions(505);
  hmc_mjb_recoil->GetYaxis()->SetTitleOffset(0.7); hmc_mjb_recoil->GetXaxis()->SetTitleOffset(1.1);
  if(usehist2Draw) { hmc_mjb_recoil->GetYaxis()->SetTitleSize(0.09); hdt_mjb_recoil->GetYaxis()->SetTitleSize(0.045); }

  TH1D *hpad_ratio = (TH1D*)hmc_mjb_recoil->Clone();
  hpad_ratio->GetXaxis()->SetMoreLogLabels(); hpad_ratio->GetXaxis()->SetNoExponent();
  TCanvas *c2 = tdrDiCanvas("c2",hdt_mjb_recoil,hpad_ratio,4,11);

  c2->cd(1); gPad->SetLogx();
  
  int npix = int(min(c2->GetWw(),c2->GetWh())*0.045*0.7);
  TLatex latex2; latex2.SetNDC(); latex2.SetTextFont(43); latex2.SetTextSize(npix); latex2.SetTextColor(kBlack);
  //latex2.DrawLatex(0.19,0.7,"|#eta_{leading}^{jet}| < 1.3");
//  latex2.DrawLatex(0.19,0.7,"p_{T} > "+pt+" GeV, |#eta_{leading}| < 1.3");
  double leg_ymin2 = 0.8; //0.72;
  latex2.DrawLatex(0.19,leg_ymin2,"p_{T} > "+pt+" GeV, |#eta_{leading}| < 1.3");
  tex->SetTextSize(0.047);
  //tex->DrawLatex(0.6275,0.9356322,"Run2018C (13 TeV)");
  //tex->DrawLatex(0.62,0.9356322,"Run2018C");
  TString str_ul = "UL"+str_dt;
  if(dtfile.Contains("ABCD")) tex->DrawLatex(0.57,0.9356322,str_ul);
  else if(dtfile.Contains("AB")) tex->DrawLatex(0.61,0.9356322,str_ul);
  else if(dtfile.Contains("CD")) tex->DrawLatex(0.61,0.9356322,str_ul);
  else tex->DrawLatex(0.64,0.9356322,str_ul);

  //hdt_mjb_recoil->SetMarkerStyle(25); hdt_mjb_recoil->SetLineColor(kCyan+1); hdt_mjb_recoil->SetMarkerColor(kCyan+1);
  //hmc_mjb_recoil->SetMarkerStyle(24); hmc_mjb_recoil->SetLineColor(kCyan+3); hmc_mjb_recoil->SetMarkerColor(kCyan+3);
  //hdt_mpf_recoil->SetMarkerStyle(21); hdt_mpf_recoil->SetLineColor(kOrange+8); hdt_mpf_recoil->SetMarkerColor(kOrange+8);
  //hmc_mpf_recoil->SetMarkerStyle(22); hmc_mpf_recoil->SetLineColor(kOrange+1); hmc_mpf_recoil->SetMarkerColor(kOrange+1);

  if (both) {
    //hdt_mjb_recoil->DrawCopy("same"); if(drawMC) hmc_mjb_recoil->DrawCopy("same");
    //hdt_mpf_recoil->DrawCopy("same"); if(drawMC) hmc_mpf_recoil->DrawCopy("same");
    //hdt_mjb_leading->DrawCopy("same"); if(drawMC) hmc_mjb_leading->DrawCopy("same");
    //hdt_mpf_leading->DrawCopy("same"); if(drawMC) hmc_mpf_leading->DrawCopy("same");
    hdt_mpf_recoil->DrawCopy("same"); if(drawMC) hmc_mpf_recoil->DrawCopy("same");
    hdt_mjb_recoil->DrawCopy("same"); if(drawMC) hmc_mjb_recoil->DrawCopy("same");
    hdt_mpf_leading->DrawCopy("same"); if(drawMC) hmc_mpf_leading->DrawCopy("same");
    hdt_mjb_leading->DrawCopy("same"); if(drawMC) hmc_mjb_leading->DrawCopy("same");
  } else if (recoil) {
    hdt_mjb_recoil->Draw("same"); if(drawMC) hmc_mjb_recoil->Draw("same");
    hdt_mpf_recoil->Draw("same"); if(drawMC) hmc_mpf_recoil->Draw("same");
  } else if (leading) {
    hdt_mjb_leading->Draw("same"); if(drawMC) hmc_mjb_leading->Draw("same");
    hdt_mpf_leading->Draw("same"); if(drawMC) hmc_mpf_leading->Draw("same");
  }

  leg->Draw(); leg2->Draw();

  Float_t uxmin = gPad->GetUxmin();
  Float_t uxmax = gPad->GetUxmax();
  Float_t uymin = gPad->GetUymin();
  cout << "c2->cd(1) uxmin = " << uxmin << ", uxmax = " << uxmax << ", uymin = " << uymin << endl;

  Double_t xlow = 50, xup = 3000;
  if(xmin>100) xlow = 300;
  Int_t nbinsx = Int_t(xup - xlow);
/*
  if (gPad->GetLogx()) {
    cout << "gPad->GetLogx()" << endl;
    //hdt_mjb_recoil->GetXaxis()->SetNdivisions(0);

    TLine tick;
    tick.SetLineWidth(1); tick.SetLineColor(1);
    Float_t length = (gPad->GetLogy()) ? 0.02 : 0.013;
    for (Int_t i=xlow; i<=xup; i+=50) {
      Float_t xx = i;
      if (gPad->GetLogy()) tick.DrawLine(xx, pow(10,uymin), xx, pow(10,uymin) + (i%100 == 0 ? 2*length : length));
      else tick.DrawLine(xx, uymin, xx, uymin + (i%100 == 0 ? 2*length : length));
    }
  }
*/
  c2->cd(2); gPad->SetLogx();

  TH1D *hdiv_mjb_recoil = (TH1D*)hdt_mjb_recoil->Clone("hdiv_mjb_recoil");
  TH1D *hdiv_mjb_leading = (TH1D*)hdt_mjb_leading->Clone("hdiv_mjb_leading");
  TH1D *hdiv_mpf_recoil = (TH1D*)hdt_mpf_recoil->Clone("hdiv_mpf_recoil");
  TH1D *hdiv_mpf_leading = (TH1D*)hdt_mpf_leading->Clone("hdiv_mpf_leading");

  hdiv_mjb_recoil->Divide(hmc_mjb_recoil); hdiv_mjb_leading->Divide(hmc_mjb_leading);
  hdiv_mpf_recoil->Divide(hmc_mpf_recoil); hdiv_mpf_leading->Divide(hmc_mpf_leading);

  //hdiv_mpf_recoil->SetMarkerStyle(2); hdiv_mpf_leading->SetMarkerStyle(27);

  if (both) {
    hdiv_mpf_recoil->Draw("same"); hdiv_mpf_leading->Draw("same");
    hdiv_mjb_recoil->Draw("same"); hdiv_mjb_leading->Draw("same");
    //hdiv_mpf_recoil->Draw("same"); hdiv_mpf_leading->Draw("esame");
  } else if (recoil) {
    hdiv_mpf_recoil->Draw("same"); hdiv_mjb_recoil->Draw("same");
  } else if (leading) {
    hdiv_mpf_leading->Draw("same"); hdiv_mjb_leading->Draw("same");
  }

  TLegend *legMJB = new TLegend(0.5,0.4,0.7,0.52);
  legMJB->SetBorderSize(0); legMJB->SetTextSize(0.05); legMJB->SetLineColor(1);
  legMJB->SetLineStyle(1); legMJB->SetLineWidth(1); legMJB->SetFillColor(0); legMJB->SetFillStyle(1001);
  legMJB->AddEntry(hdiv_mjb_recoil,"MJB, p_{T, recoil}","P");
  legMJB->AddEntry(hdiv_mjb_leading,"MJB, p_{T, leading}","P");
  //legMJB->Draw();

  TLegend *legMPF = new TLegend(0.72,0.4,0.85,0.52);
  legMPF->SetBorderSize(0); legMPF->SetTextSize(0.05); legMPF->SetLineColor(1);
  legMPF->SetLineStyle(1); legMPF->SetLineWidth(1); legMPF->SetFillColor(0); legMPF->SetFillStyle(1001);
  legMPF->AddEntry(hdiv_mpf_recoil,"MPF, p_{T, recoil}","P");
  legMPF->AddEntry(hdiv_mpf_leading,"MPF, p_{T, leading}","P");
  //legMPF->Draw(); 

  Float_t ymax = 1;
  TLine *line = new TLine(xmin,ymax,5000,ymax);
  line->SetLineColor(kBlack); line->SetLineStyle(2);
  line->Draw();

  if(xmin==15 or xmin==32) {
    hmc_mjb_recoil->GetXaxis()->SetMoreLogLabels(); hmc_mjb_recoil->GetXaxis()->SetNoExponent();
  } else {

  uxmin = gPad->GetUxmin();
  uxmax = gPad->GetUxmax();
  uymin = gPad->GetUymin();
  cout << "c2->cd(2) uxmin = " << uxmin << ", uxmax = " << uxmax << ", uymin = " << uymin << endl;
  if(uymin>0.75 and uymin<0.85) { uymin *= 1.015; cout << "--------------------------> updated uymin = " << uymin << endl; }
  if(uymin>0.85) { uymin *= 1.012; cout << "--------------------------> updated uymin = " << uymin << endl; }

  //Double_t xlow = 100, xup = 1000;
  //Int_t nbinsx = Int_t(xup - xlow);

  if (gPad->GetLogx()) {
    cout << "gPad->GetLogx()" << endl;
    //hmc_mjb_recoil->GetXaxis()->SetNdivisions(0);
    hmc_mjb_recoil->GetXaxis()->SetLabelSize(0);

    Float_t length = (gPad->GetLogy()) ? 0.02 : 0.013;
//    TLine tick;
//    tick.SetLineWidth(1); tick.SetLineColor(1);
//    for (Int_t i=xlow; i<=xup; i+=50) {
//      Float_t xx = i;
//      if (gPad->GetLogy()) tick.DrawLine(xx, pow(10,uymin), xx, pow(10,uymin) + (i%100 == 0 ? 2*length : length));
//      else tick.DrawLine(xx, uymin, xx, uymin + (i%100 == 0 ? 2*length : length));
//    }
    Float_t ylatex = (gPad->GetLogy()) ? pow(10,uymin) - 1.3*length : uymin - 2*length;
    for (Int_t i=xlow; i<=xup; i+=25) {
      if(i>100) i+=75;
      if(i>400) i+=100; if(i>600) i+=200; if(i>1000) i+=600;
      Float_t xx = i;
    
      TLatex* latex = new TLatex(xx, ylatex, Form("%.0f", xx));
      
      latex->SetTextAlign(  22);
      latex->SetTextFont (  42);
      //latex->SetTextFont (  62);
      latex->SetTextSize (0.07);
    
      latex->Draw("same");
    }

  }
}
  //TString cvname = "Run"+str_dt+"_"+mcfile+"_"+metType+l2res; //cout << cvname << endl;
  TString CVname = "UL"+str_dt+"_"+sample+"_"+metType+l2res;
  cout << CVname << " " << jerSF << " " << jetCorrVer << endl;
  //c2->Print(CVname+".pdf");
  TString _drawMC = drawMC ? "_drawMC" : "";
  if(! usehist2Draw) c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_prof"+_drawMC+".pdf"), c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_prof"+_drawMC+".png");
  else if(mpfprofDraw) c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_MPFprof"+_drawMC+".pdf"), c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_MPFprof"+_drawMC+".png");
  else {
    if(xmin==300) c2->Print(CVname+jetCorrVer+jerSF+"_Pt"+pt+_drawMC+".pdf"), c2->Print(CVname+jetCorrVer+jerSF+"_Pt"+pt+_drawMC+".png");
    else {
      if(noRebin) c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_noRebin"+_drawMC+".pdf"), c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+"_noRebin"+_drawMC+".png");
      else c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+_drawMC+".pdf"), c2->Print(CVname+jetCorrVer+jerSF+"_xmin"+xmin+"_Pt"+pt+_drawMC+".png");
    }
  }

  TCanvas *ccv = new TCanvas("ccv","CRecoil",700,100,600,600);
  pCRecoil_dt->Draw(); pCRecoil_mc->Draw("same"); //pCRecoil_check->Draw("histsame");

  double xmin2 = 50;
  pCRecoil_dt->GetXaxis()->SetRangeUser(xmin2,xmax); pCRecoil_mc->GetXaxis()->SetRangeUser(xmin2,xmax); //pCRecoil_check->GetXaxis()->SetRangeUser(xmin2,xmax);
  gPad->SetLogx(); pCRecoil_dt->GetXaxis()->SetMoreLogLabels(); pCRecoil_dt->GetXaxis()->SetNoExponent(); pCRecoil_dt->GetYaxis()->SetRangeUser(0.35,0.55);

  pCRecoil_dt->GetYaxis()->SetTitle("C_{recoil}"); pCRecoil_dt->SetLineColor(kBlue-4); pCRecoil_dt->SetMarkerColor(kBlue-4); pCRecoil_dt->SetMarkerStyle(20);
  pCRecoil_mc->GetYaxis()->SetTitle("C_{recoil}"); pCRecoil_mc->SetLineColor(kBlue-4); pCRecoil_mc->SetMarkerColor(kBlue-4); pCRecoil_mc->SetMarkerStyle(24);
  pCRecoil_ptcl->GetYaxis()->SetTitle("C_{recoil ptcl}"); pCRecoil_ptcl2->GetYaxis()->SetTitle("C_{recoil ptcl}");

  //TLegend *legCRecoil = new TLegend(0.6,0.33,0.75,0.52,NULL,"brNDC");
  TLegend *legCRecoil = new TLegend(0.4,0.68,0.55,0.87,NULL,"brNDC");
  legCRecoil->SetBorderSize(0); legCRecoil->SetTextSize(0.04); legCRecoil->SetLineColor(1);
  legCRecoil->SetLineStyle(1); legCRecoil->SetLineWidth(1); legCRecoil->SetFillColor(0); legCRecoil->SetFillStyle(1001);
  if(! mcfile.Contains("MG")) legCRecoil->AddEntry(pCRecoil_mc,"MC ("+sample+")","ep");
  else legCRecoil->AddEntry(pCRecoil_mc,"MC (MGP8)","ep");
  legCRecoil->AddEntry(pCRecoil_dt,"data ("+str_dt+")","ep");
  //legCRecoil->AddEntry(pCRecoil_check,"data check","l");
  legCRecoil->Draw();

  TH1 *hdt_CRecoil = pCRecoil_dt->ProjectionX("hdt_CRecoil");
  hdt_CRecoil->SetLineColor(1); hdt_CRecoil->SetMarkerColor(1); hdt_CRecoil->SetMarkerStyle(1); hdt_CRecoil->SetStats(0);
  TH1 *hmc_CRecoil = pCRecoil_mc->ProjectionX("hmc_CRecoil"); 
  TH1 *hmc_CRecoil_ptcl = pCRecoil_ptcl->ProjectionX("hmc_CRecoil_ptcl");
  TH1 *hmc_CRecoil_ptcl2 = pCRecoil_ptcl2->ProjectionX("hmc_CRecoil_ptcl2");
  hmc_CRecoil->SetLineColor(1); hmc_CRecoil->SetMarkerColor(1); hmc_CRecoil->SetMarkerStyle(1); hmc_CRecoil->SetStats(0);
  hmc_CRecoil_ptcl->SetLineColor(1); hmc_CRecoil_ptcl->SetMarkerColor(1); hmc_CRecoil_ptcl->SetMarkerStyle(1); hmc_CRecoil_ptcl->SetStats(0);
  hmc_CRecoil_ptcl2->SetLineColor(1); hmc_CRecoil_ptcl2->SetMarkerColor(1); hmc_CRecoil_ptcl2->SetMarkerStyle(1); hmc_CRecoil_ptcl2->SetStats(0);

  TH1 *hdt_CRecoil_leading = 0, *hdt_CRecoil_ptave = 0, *hmc_CRecoil_leading = 0, *hmc_CRecoil_ptave = 0;
  hdt_CRecoil_leading = pCRecoil_leading_dt->ProjectionX("hdt_CRecoil_leading"), hdt_CRecoil_ptave = pCRecoil_ptave_dt->ProjectionX("hdt_CRecoil_ptave");
  hdt_CRecoil_leading->SetStats(0); hdt_CRecoil_ptave->SetStats(0); hdt_CRecoil_leading->SetMarkerStyle(1); hdt_CRecoil_ptave->SetMarkerStyle(1);
  hmc_CRecoil_leading = pCRecoil_leading_mc->ProjectionX("hmc_CRecoil_leading"), hmc_CRecoil_ptave = pCRecoil_ptave_mc->ProjectionX("hmc_CRecoil_ptave");
  hmc_CRecoil_leading->SetStats(0); hmc_CRecoil_ptave->SetStats(0); hmc_CRecoil_leading->SetMarkerStyle(1); hmc_CRecoil_ptave->SetMarkerStyle(1);

  TH1D *hdt_mpf_recoil2 = pmpf_recoil2_rebin->ProjectionX("hdt_mpf_recoil2"); //use rebin!
  TH1D *hmc_mpf_recoil2 = pmpf_recoil2_mc->ProjectionX("hmc_mpf_recoil2");
  //--Added on Sep 15 after finding a but of x-axis
  TH1D *hdt_mpf_leading2 = pmpf_leading2->ProjectionX("hdt_mpf_leading2");
  TH1D *hmc_mpf_leading2 = pmpf_leading2_mc->ProjectionX("hmc_mpf_leading2");
  if(usehist2) {
    hdt_mpf_leading2 = (TH1D*)hmpf_leading2->Clone("hdt_mpf_leading2");
    hmc_mpf_leading2 = (TH1D*)hmpf_leading2_mc->Clone("hmc_mpf_leading2");
  }
  TH1D *hdt_mpfleadinv_leading2 = 0, *hmc_mpfleadinv_leading2 = 0;
  hdt_mpfleadinv_leading2 = pmpfleadinv_leading2->ProjectionX("hdt_mpfleadinv_leading2"); hdt_mpfleadinv_leading2->SetStats(0); hdt_mpfleadinv_leading2->SetMarkerStyle(1);
  hmc_mpfleadinv_leading2 = pmpfleadinv_leading2_mc->ProjectionX("hmc_mpfleadinv_leading2"); hmc_mpfleadinv_leading2->SetStats(0); hmc_mpfleadinv_leading2->SetMarkerStyle(1);

  hdt_mpf_leading2->GetYaxis()->SetLabelSize(0.035); hmc_mpf_leading2->GetYaxis()->SetLabelSize(0.035);
  hdt_mpf_leading2->GetYaxis()->SetTitleSize(0.035); hmc_mpf_leading2->GetYaxis()->SetTitleSize(0.035);

  hdt_mpf_recoil2->GetYaxis()->SetTitleOffset(1); hmc_mpf_recoil2->GetYaxis()->SetTitleOffset(1);
  hdt_mpf_leading2->GetYaxis()->SetTitleOffset(1); hmc_mpf_leading2->GetYaxis()->SetTitleOffset(1);

  hdt_mpf_leading2->GetYaxis()->SetLabelOffset(0.005); hmc_mpf_leading2->GetYaxis()->SetLabelOffset(0.005);

  hdt_mpf_recoil2->SetLineColor(1); hmc_mpf_recoil2->SetLineColor(1);
  hdt_mpf_leading2->SetLineColor(1); hmc_mpf_leading2->SetLineColor(1);

  hdt_mpf_recoil2->SetMarkerColor(1); hmc_mpf_recoil2->SetMarkerColor(1);
  hdt_mpf_leading2->SetMarkerColor(1); hmc_mpf_leading2->SetMarkerColor(1);

  hdt_mpf_recoil2->SetMarkerStyle(1); hmc_mpf_recoil2->SetMarkerStyle(1);
  hdt_mpf_leading2->SetMarkerStyle(1); hmc_mpf_leading2->SetMarkerStyle(1);

  hdt_mpf_recoil2->SetMinimum(0); hmc_mpf_recoil2->SetMinimum(0);
  hdt_mpf_leading2->SetMinimum(0); hmc_mpf_leading2->SetMinimum(0);

  hdt_mpf_recoil2->SetStats(0); hmc_mpf_recoil2->SetStats(0);
  hdt_mpf_leading2->SetStats(0); hmc_mpf_leading2->SetStats(0);

if(writeOutput) {
  hmc_mjb_recoil->GetYaxis()->SetNdivisions(510);
  hmc_mjb_recoil->GetYaxis()->SetLabelSize(0.03); hmc_mjb_recoil->GetXaxis()->SetTickLength(0.03);
  hmc_mjb_recoil->GetYaxis()->SetLabelSize(0.035); hmc_mjb_recoil->GetXaxis()->SetLabelSize(0.035);
  hmc_mjb_recoil->GetYaxis()->SetTitleSize(0.035); hmc_mjb_recoil->GetXaxis()->SetTitleSize(0.035);
  hmc_mjb_recoil->GetYaxis()->SetTitleOffset(1); hmc_mjb_recoil->GetXaxis()->SetTitleOffset(1);

  hdt_mjb_recoil->GetYaxis()->SetTitleSize(0.035); hdt_mjb_recoil->GetYaxis()->SetLabelSize(0.035);
  hdt_mpf_recoil->GetYaxis()->SetTitleSize(0.035); hdt_mpf_recoil->GetYaxis()->SetLabelSize(0.035);
  hmc_mpf_recoil->GetYaxis()->SetTitleSize(0.035); hmc_mpf_recoil->GetYaxis()->SetLabelSize(0.035);

  hdt_mjb_leading->GetYaxis()->SetLabelSize(0.035); hdt_mjb_leading->GetYaxis()->SetTitleSize(0.035);
  hdt_mpf_leading->GetYaxis()->SetLabelSize(0.035); hdt_mpf_leading->GetYaxis()->SetTitleSize(0.035);
  hmc_mjb_leading->GetYaxis()->SetLabelSize(0.035); hmc_mjb_leading->GetYaxis()->SetTitleSize(0.035);
  hmc_mpf_leading->GetYaxis()->SetLabelSize(0.035); hmc_mpf_leading->GetYaxis()->SetTitleSize(0.035);

  hdt_mjb_recoil->GetYaxis()->SetTitleOffset(1); hdt_mjb_leading->GetYaxis()->SetTitleOffset(1);
  hdt_mpf_recoil->GetYaxis()->SetTitleOffset(1); hdt_mpf_leading->GetYaxis()->SetTitleOffset(1);
  hmc_mjb_recoil->GetYaxis()->SetTitleOffset(1); hmc_mjb_leading->GetYaxis()->SetTitleOffset(1);
  hmc_mpf_recoil->GetYaxis()->SetTitleOffset(1); hmc_mpf_leading->GetYaxis()->SetTitleOffset(1);

  hdt_mjb_leading->GetYaxis()->SetLabelOffset(0.005);
  hdt_mpf_leading->GetYaxis()->SetLabelOffset(0.005);
  hmc_mjb_leading->GetYaxis()->SetLabelOffset(0.005);
  hmc_mpf_leading->GetYaxis()->SetLabelOffset(0.005);

  hdt_mjb_recoil->SetName("hdt_mjb_recoil"); hdt_mjb_leading->SetName("hdt_mjb_leading");
  hdt_mpf_recoil->SetName("hdt_mpf_recoil"); hdt_mpf_leading->SetName("hdt_mpf_leading");
  hmc_mjb_recoil->SetName("hmc_mjb_recoil"); hmc_mjb_leading->SetName("hmc_mjb_leading");
  hmc_mpf_recoil->SetName("hmc_mpf_recoil"); hmc_mpf_leading->SetName("hmc_mpf_leading");

  hdt_mjb_recoil->SetLineColor(1); hdt_mjb_leading->SetLineColor(1);
  hdt_mpf_recoil->SetLineColor(1); hdt_mpf_leading->SetLineColor(1);
  hmc_mjb_recoil->SetLineColor(1); hmc_mjb_leading->SetLineColor(1);
  hmc_mpf_recoil->SetLineColor(1); hmc_mpf_leading->SetLineColor(1);

  hdt_mjb_recoil->SetMarkerColor(1); hdt_mjb_leading->SetMarkerColor(1);
  hdt_mpf_recoil->SetMarkerColor(1); hdt_mpf_leading->SetMarkerColor(1);
  hmc_mjb_recoil->SetMarkerColor(1); hmc_mjb_leading->SetMarkerColor(1);
  hmc_mpf_recoil->SetMarkerColor(1); hmc_mpf_leading->SetMarkerColor(1);

  hdt_mjb_recoil->SetMarkerStyle(1); hdt_mjb_leading->SetMarkerStyle(1);
  hdt_mpf_recoil->SetMarkerStyle(1); hdt_mpf_leading->SetMarkerStyle(1);
  hmc_mjb_recoil->SetMarkerStyle(1); hmc_mjb_leading->SetMarkerStyle(1);
  hmc_mpf_recoil->SetMarkerStyle(1); hmc_mpf_leading->SetMarkerStyle(1);

  hdt_mjb_recoil->SetMinimum(0); hdt_mjb_leading->SetMinimum(0);
  hdt_mpf_recoil->SetMinimum(0); hdt_mpf_leading->SetMinimum(0);
  hmc_mjb_recoil->SetMinimum(0); hmc_mjb_leading->SetMinimum(0);
  hmc_mpf_recoil->SetMinimum(0); hmc_mpf_leading->SetMinimum(0);

  hdt_mjb_recoil->SetStats(0); hdt_mjb_leading->SetStats(0);
  hdt_mpf_recoil->SetStats(0); hdt_mpf_leading->SetStats(0);
  hmc_mjb_recoil->SetStats(0); hmc_mjb_leading->SetStats(0);
  hmc_mpf_recoil->SetStats(0); hmc_mpf_leading->SetStats(0);

  TString metType2 = "_metType2";
  TString dataCorr = (usel2res) ? "_L1L2Res" : "_L1L2L3Res";
  TString jerScale = ""; if(jer>0) jerScale = "_JER-up"; else if(jer<0) jerScale = "_JER-down";
  TString option = (usel2res || jer!=0) ? "UPDATE" : "RECREATE";
  if(pt!="30") option = "UPDATE";

  if(writeOutput) {
  //hdt_mjb_recoil->GetYaxis()->SetTitle("MJB"+dataCorr); hdt_mjb_leading->GetYaxis()->SetTitle("MJB"+dataCorr);
  //hdt_mpf_recoil->GetYaxis()->SetTitle("MPF"+dataCorr); hdt_mpf_leading->GetYaxis()->SetTitle("MPF"+dataCorr);
  hdt_mjb_recoil->GetYaxis()->SetTitle("MJB"); hdt_mjb_leading->GetYaxis()->SetTitle("MJB");
  hdt_mpf_recoil->GetYaxis()->SetTitle("MPF"); hdt_mpf_leading->GetYaxis()->SetTitle("MPF");
  hmc_mjb_recoil->GetYaxis()->SetTitle("MJB"); hmc_mjb_leading->GetYaxis()->SetTitle("MJB");
  hmc_mpf_recoil->GetYaxis()->SetTitle("MPF"); hmc_mpf_leading->GetYaxis()->SetTitle("MPF");

  //hdt_mpf_recoil2->GetYaxis()->SetTitle("MPF"+metType2); hmc_mpf_recoil2->GetYaxis()->SetTitle("MPF"+metType2);
  //hdt_mpf_leading2->GetYaxis()->SetTitle("MPF"+metType2); hmc_mpf_leading2->GetYaxis()->SetTitle("MPF"+metType2);
  hdt_mpf_recoil2->GetYaxis()->SetTitle("MPF"); hmc_mpf_recoil2->GetYaxis()->SetTitle("MPF");
  hdt_mpf_leading2->GetYaxis()->SetTitle("MPF"); hmc_mpf_leading2->GetYaxis()->SetTitle("MPF");
  }

  TString foutput = "multijet_"+writeDate+"_Run"+str_dt+"_"+sample+jetCorrVer+".root";
  jetCorrVer = "_jecV5_jerV2";
  foutput = "multijet_"+writeDate+"_UL"+str_dt+jetCorrVer+".root";

  TFile *outfile = new TFile(foutput,option);

  TString dirName = "Data";
  if(option == "RECREATE") {
    cout << "RECREATE " << foutput << " with directories, e.g, Data/Pt" << pt << endl;
    outfile->mkdir(dirName); outfile->mkdir(dirName+"/Pt"+pt); 
    if(writeMC) {
      outfile->mkdir("MC"); outfile->mkdir("MC_JER-up"); outfile->mkdir("MC_JER-down");
      outfile->mkdir("MC/Pt"+pt); outfile->mkdir("MC_JER-up/Pt"+pt); outfile->mkdir("MC_JER-down/Pt"+pt);
    }
    if(pt=="30") {
      cout << "-> Adding Pt20 and Pt15 directories for both Data and MC" << endl;
      outfile->mkdir(dirName+"/Pt20"); outfile->mkdir(dirName+"/Pt15");
      if(writeMC) {
        outfile->mkdir("MC/Pt20"); outfile->mkdir("MC/Pt15");
        outfile->mkdir("MC_JER-up/Pt20"); outfile->mkdir("MC_JER-down/Pt20"); outfile->mkdir("MC_JER-up/Pt15"); outfile->mkdir("MC_JER-down/Pt15");
      }
    } else if(pt=="15"){
      outfile->mkdir(dirName+"/Pt20"); outfile->mkdir(dirName+"/Pt30");
      if(writeMC) {
        outfile->mkdir("MC/Pt20"); outfile->mkdir("MC/Pt30");
        outfile->mkdir("MC_JER-up/Pt20"); outfile->mkdir("MC_JER-down/Pt20"); outfile->mkdir("MC_JER-up/Pt30"); outfile->mkdir("MC_JER-down/Pt30");
      }
    }
  } else if(jer!=0) cout << "Updating " << foutput << " with " << jerScale << endl;
  else if(usel2res) cout << "Updating " << foutput << " with " << dataCorr << endl;
  else cout << "Updating " << foutput << " with " << sample << endl;
  if(option == "UPDATE") {
//    cout << "Updating " << foutput << " with addeing directories of Data/Pt" << pt << " and MC/Pt" << pt << endl;
//    outfile->mkdir("Data/Pt"+pt); outfile->mkdir("MC/Pt"+pt);
    cout << " Are there the directories for Data/Pt" << pt << " and MC/Pt" << pt << " ?? " << endl;
  }

  if(jer == 0) {
    cout << " -> " << dirName+"/Pt"+pt << endl;
    gDirectory->cd(dirName+"/Pt"+pt);

    h2mjb_recoil_rebin->Write("h2_MJB_recoil"+dataCorr); h2mpf_recoil_rebin->Write("h2_MPF_recoil"+dataCorr); h2mpf_recoil2_rebin->Write("h2_MPF_recoil"+dataCorr+metType2);
    h2mjbinv_leading->GetYaxis()->SetTitle("p_{T,recoil}/p_{T,leading}"); h2mjbinv_leading->Write("h2_MJB0_leading"+dataCorr);
    h2mjb_leading->Write("h2_MJBx_leading"+dataCorr);
    h2mpf_leading->Write("h2_MPF_leading"+dataCorr); h2mpf_leading2->Write("h2_MPF_leading"+dataCorr+metType2);
    h2mpfleadinv_leading->Write("h2_MPF0_leading"+dataCorr), h2mpfleadinv_leading2->Write("h2_MPF0_leading"+dataCorr+metType2);
    h2mjb_ptave->Write("h2_MJB0_ptave"+dataCorr), h2mpf_ptave->Write("h2_MPF0_ptave"+dataCorr), h2mpf_ptave2->Write("h2_MPF0_ptave"+dataCorr+metType2);
    h2mjbinv_ptave->Write("h2_MJBx_ptave"+dataCorr), h2mpfinv_ptave->Write("h2_MPFx_ptave"+dataCorr), h2mpfinv_ptave2->Write("h2_MPFx_ptave"+dataCorr+metType2);
    h2mpflead_leading_ue->Write("h2_MPFue_leading"+dataCorr); h2mpfleadinv_leading_ue->Write("h2_MPFue0_leading"+dataCorr); ////h2mpflead_ptave_ue->Write("h2_MPFue_ptave"+dataCorr);

    hdt_mjb_recoil->Write("MJB_recoil"+dataCorr); hdt_mjb_leading->Write("MJB_leading"+dataCorr); hdt_mjb0_leading->Write("MJB0_leading"+dataCorr); hdt_mjbx_leading->Write("MJBx_leading"+dataCorr);
    hdt_mpf_recoil->Write("MPF_recoil"+dataCorr); hdt_mpf_leading->Write("MPF_leading"+dataCorr); hdt_mpfleadinv_leading->Write("MPF0_leading"+dataCorr);
    hdt_mpf_recoil2->Write("MPF_recoil"+dataCorr+metType2); hdt_mpf_leading2->Write("MPF_leading"+dataCorr+metType2); hdt_mpfleadinv_leading2->Write("MPF0_leading"+dataCorr+metType2);
    hdt_mjb_ptave->Write("MJB_ptave"+dataCorr), hdt_mpf_ptave->Write("MPF_ptave"+dataCorr), hdt_mpf_ptave2->Write("MPF_ptave"+dataCorr+metType2);
    hdt_mjbinv_ptave->Write("MJBx_ptave"+dataCorr), hdt_mpfinv_ptave->Write("MPFx_ptave"+dataCorr), hdt_mpfinv_ptave2->Write("MPFx_ptave"+dataCorr+metType2);
    hdt_mpflead_leading_ue->Write("MPFue_leading"+dataCorr); hdt_mpfleadinv_leading_ue->Write("MPFue0_leading"+dataCorr); ////hdt_mpflead_ptave_ue->Write("MPFue_ptave"+dataCorr);
    hdt_CRecoil->Write("CRecoil"+dataCorr); hdt_CRecoil_leading->Write("CRecoil_leading"+dataCorr); hdt_CRecoil_ptave->Write("CRecoil_ptave"+dataCorr);
  }

  if(writeMC && !usel2res) {

    gDirectory->cd("/MC"+jerScale+"/Pt"+pt);
    //jerScale = ""; //updated on Sep 11, 2019
    jerScale = suffix; //updated for Mikko (Sep 12) 

    h2mjb_recoil_mc->Write("h2_MJB_recoil"+jerScale); h2mpf_recoil_mc->Write("h2_MPF_recoil"+jerScale); h2mpf_recoil2_mc->Write("h2_MPF_recoil"+jerScale+metType2);
    h2mjbinv_leading_mc->GetYaxis()->SetTitle("p_{T,recoil}/p_{T,leading}"); h2mjbinv_leading_mc->Write("h2_MJB0_leading"+jerScale);
    h2mjb_leading_mc->Write("h2_MJBx_leading"+jerScale);
    h2mpf_leading_mc->Write("h2_MPF_leading"+jerScale); h2mpf_leading2_mc->Write("h2_MPF_leading"+jerScale+metType2);
    h2mpfleadinv_leading_mc->Write("h2_MPF0_leading"+jerScale), h2mpfleadinv_leading2_mc->Write("h2_MPF0_leading"+jerScale+metType2);
    h2mjb_ptave_mc->Write("h2_MJB0_ptave"+jerScale), h2mpf_ptave_mc->Write("h2_MPF0_ptave"+jerScale), h2mpf_ptave2_mc->Write("h2_MPF0_ptave"+jerScale+metType2);
    h2mjbinv_ptave_mc->Write("h2_MJBx_ptave"+jerScale), h2mpfinv_ptave_mc->Write("h2_MPFx_ptave"+jerScale), h2mpfinv_ptave2_mc->Write("h2_MPFx_ptave"+jerScale+metType2);
    h2mpflead_leading_ue_mc->Write("h2_MPFue_leading"+jerScale); h2mpfleadinv_leading_ue_mc->Write("h2_MPFue0_leading"+jerScale); ////h2mpflead_ptave_ue_mc->Write("h2_MPFue_ptave"+jerScale);

    hmc_mjb_recoil->Write("MJB_recoil"+jerScale); hmc_mjb_leading->Write("MJB_leading"+jerScale); hmc_mjb0_leading->Write("MJB0_leading"+jerScale); hdt_mjbx_leading->Write("MJBx_leading"+jerScale);
    hmc_mpf_recoil->Write("MPF_recoil"+jerScale); hmc_mpf_leading->Write("MPF_leading"+jerScale); hmc_mpfleadinv_leading->Write("MPF0_leading"+jerScale);
    hmc_mpf_recoil2->Write("MPF_recoil"+jerScale+metType2); hmc_mpf_leading2->Write("MPF_leading"+jerScale+metType2); hmc_mpfleadinv_leading2->Write("MPF0_leading"+jerScale+metType2);
    hmc_mjb_ptave->Write("MJB_ptave"+jerScale), hmc_mpf_ptave->Write("MPF_ptave"+jerScale), hmc_mpf_ptave2->Write("MPF_ptave"+jerScale+metType2);
    hmc_mjbinv_ptave->Write("MJBx_ptave"+jerScale), hmc_mpfinv_ptave->Write("MPFx_ptave"+jerScale), hmc_mpfinv_ptave2->Write("MPFx_ptave"+jerScale+metType2);
    hmc_mpflead_leading_ue->Write("MPFue_leading"+jerScale); hmc_mpfleadinv_leading_ue->Write("MPFue0_leading"+jerScale); ////hmc_mpflead_ptave_ue->Write("MPFue_ptave"+jerScale);
 
    hmc_CRecoil->Write("CRecoil"+jerScale); hmc_CRecoil_ptcl->Write("CRecoil_ptcl"+jerScale); hmc_CRecoil_ptcl2->Write("CRecoil_ptcl2"+jerScale);
    hmc_CRecoil_leading->Write("CRecoil_leading"+jerScale); hmc_CRecoil_ptave->Write("CRecoil_ptave"+jerScale);
  }

  outfile->Write(); outfile->Close();

}

}

