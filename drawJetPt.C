#include "tdrstyle_mod18.C"
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF2.h"
#include "TGraphErrors.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

bool doIOVeff = true; // correct ECALprefire per IOV

// "Rebin(2)"
const int nptb = 51+10+6;
const float ptbins[nptb+1] =
  {15, 18, 21, 24, 28, 32, 37, 43, 49, 56,
   64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362,
   395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032,
   1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116,
   2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
   4037, 4252, 4477, 4713, 4961, 5220};

TH1F *rebinXsec(TH1F *h) {
  //assert(h->GetBinLowEdge(1)<=pt[0]);
  //assert(h->GetBinLowEdge(h->GetNbinsX()+1)>=pt[nptb]);
  TH1F *h2 = new TH1F(Form("%s_rb2",h->GetName()),"",nptb,ptbins);
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    float pt1 = h2->GetBinLowEdge(i);
    float pt2 = h2->GetBinLowEdge(i+1);
    int j1 = h->FindBin(pt1);
    int j2 = h->FindBin(pt2)-1;
    if (!(h2->GetBinLowEdge(i)==h->GetBinLowEdge(j1))) {
      cout << h2->GetBinLowEdge(i) << ", " << h->GetBinLowEdge(j1)
	   << endl << flush;
      assert(false);
    }
    if (!(h2->GetBinLowEdge(i+1)==h->GetBinLowEdge(j2+1))) {
      cout << h2->GetBinLowEdge(i+1) << ", " << h->GetBinLowEdge(j2+1)
	   << endl << flush;
      assert(false);
    }

    double sumx(0);
    double errx2(0);
    for (int j = j1; j != j2+1; ++j) {
      sumx += h->GetBinContent(j) * h->GetBinWidth(j);
      if (h->GetBinContent(j)!=0)
	errx2 += pow(h->GetBinError(j) * h->GetBinWidth(j),2);
    }
    h2->SetBinContent(i, sumx / h2->GetBinWidth(i));
    h2->GetBinError(i, sqrt(errx2) / h2->GetBinWidth(i));
  }

  return h2;
}


const double ptmin = 114;//64;
void drawJetPt(string sy = "0.0-0.5") {

  // string to char* for easier use with TForm later
  const char *cy = sy.c_str();

  // pdflatex cannot handle periods inside file name so remove them
  string sy2 = string(TString(cy).ReplaceAll(".",""));
  const char *cy2 = sy2.c_str();

  setTDRStyle();
  
  TDirectory *curdir = gDirectory;

  // Common format ROOT tuple from Engin Eren
  // https://gitlab.cern.ch/CMS-SMP-J/InclusiveJetsLegacy/blob/master/
  // common2016.root
  //TFile *fin = new TFile("rootfiles/common2016-20170204.root","READ");
  //TFile *fin = new TFile("rootfiles/common2016-03FebV3.root","READ");
  TFile *fin = new TFile("rootfiles/common2016-03FebV3-2.root","READ");
  assert(fin && !fin->IsZombie());
  
  //TFile *fin2 = new TFile("rootfiles/common2016_Legacy.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2016_Legacy_v2.root","READ");
  TFile *fin2 = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  assert(fin2 && !fin2->IsZombie());

  TFile *fin3 = new TFile("rootfiles/common2017_V11.root","READ");
  assert(fin3 && !fin3->IsZombie());

  TFile *fu = new TFile("rootfiles/unfold.root","READ");
  assert(fu && !fu->IsZombie());

  // ECAL prefire efficiency file
  TFile *fe = new TFile("rootfiles/ECALprefireEff.root");
  assert(fe && !fe->IsZombie());
  TH2D *heff = (TH2D*)fe->Get("h_JetHT_AllEras_inclusive"); assert(heff);
  //
  TFile *fe2 = new TFile("rootfiles/ECALprefireEffRebinned.root");
  assert(fe2 && !fe2->IsZombie());
  TH2D *heff2bcd = (TH2D*)fe2->Get("h_JetHT_Run2016BCD_spike_finor_pteta");
  assert(heff2bcd);
  TH2D *heff2ef = (TH2D*)fe2->Get("h_JetHT_Run2016EF_spike_finor_pteta");
  assert(heff2ef);
  TH2D *heff2gh = (TH2D*)fe2->Get("h_JetHT_Run2016GH_spike_finor_pteta");
  assert(heff2gh);
  TH2D *heff2 = (TH2D*)heff2bcd->Clone("heff2");
  heff2->Scale(12.9/36.5);
  heff2->Add(heff2ef,6.8/36.5);
  heff2->Add(heff2gh,16.8/36.5);
  
  const char *p = "CondFormats/JetMETObjects/data/";
  // 2016 uncertainty (25 ns)
  //const char *t = "Spring16_23Sep2016GHV1_DATA";
  //const char *t = "Summer16_23Sep2016GV3_DATA";
  //const char *t = "Summer16_03Feb2017G_V2_DATA";
  //const char *t = "Summer16_03Feb2017_V3_DATA";
  const char *t = "Summer16_07Aug2017GH_V10_DATA";
  const char *a = "AK4";
  const char *s = Form("%s%s_Uncertainty_%sPFchs.txt",p,t,a);
  cout<<"**"<<s<<endl<<flush;
  JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(s);

  const char *ss1 = Form("%s%s_UncertaintySources_%sPFchs.txt",p,t,a);
  const char *sss1 = "RelativeBal";//Sample";
  cout<<"**"<<ss1<<":"<<sss1<<endl<<flush;
  JetCorrectorParameters *ps1 = new JetCorrectorParameters(ss1,sss1);
  JetCorrectionUncertainty *uncs1 = new JetCorrectionUncertainty(*ps1);

  const char *ss2 = Form("%s%s_UncertaintySources_%sPFchs.txt",p,
			"Summer16_03Feb2017G_V2_DATA",a);
  const char *sss2a = "RelativePtEC1";
  cout<<"**"<<ss2<<":"<<sss2a<<endl<<flush;
  JetCorrectorParameters *ps2a = new JetCorrectorParameters(ss2,sss2a);
  JetCorrectionUncertainty *uncs2a = new JetCorrectionUncertainty(*ps2a);
  const char *sss2b = "RelativePtEC2";
  cout<<"**"<<ss2<<":"<<sss2b<<endl<<flush;
  JetCorrectorParameters *ps2b = new JetCorrectorParameters(ss2,sss2b);
  JetCorrectionUncertainty *uncs2b = new JetCorrectionUncertainty(*ps2b);

  // 2015 uncertainty (50 ns)
  const char *t0 = "Summer15_50nsV5_DATA";
  const char *a0 = "AK4";
  const char *s0 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t0,a0);
  cout<<"**"<<s0<<endl<<flush;
  JetCorrectionUncertainty *unc0 = new JetCorrectionUncertainty(s0);

  // 2012 uncertainty (8 TeV)
  const char *t1 = "Winter14_V8_DATA";
  const char *a1 = "AK5";
  const char *s1 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t1,a1);
  cout<<"**"<<s1<<endl<<flush;
  JetCorrectionUncertainty *unc1 = new JetCorrectionUncertainty(s1);

  // 5 eras: BCDEFGH, BCD, EF, G and H
  const int nera = 5;//5;//4;//5;
  //TFile *fins[nera] = {fin, fin, fin, fin, fin};
  //const char *eras[nera] = {"BCDEFGH","BCD","EFearly","FlateG","H"};
  //TFile *fins[nera] = {fin, fin, fin, fin, fin2};
  TFile *fins[nera] = {fin2, fin2, fin2, fin2, fin3};
  //const char *
  //string eras[nera] = {"CD","EF","G","H", "BCDEFGH"};
  string eras[nera] = {"BCD","EF","GH", "BCDEFGH", "B"};
  // luminosity re-normalization
  //double lumi[nera] = {1, 1, 1, 1, 1};
  //double lumi[nera] = {1, 1, 1, 1, 1};
  double lumi[nera] = {1, 1, 1, 1, 4800};

  // Mithat's settings
  /*
  int color[nera] = {kBlack, kRed+2, kOrange+2, kGreen+2, kBlue+1};
  int marker[nera] = {kFullSquare, kFullCircle, kFullDiamond,
		      kFullStar, kOpenStar};
  const char* label[nera] = {"BCDEFGH",
			     "BCD",
			     "EF^{1}",
			     "F^{2}G",
			     "H"};
  */
  //int color[nera] = {kRed+2, kOrange+2, kMagenta+2, kBlue+1, kBlack};
  int color[nera] = {kRed+2, kOrange+2, kBlue+1, kBlack, kBlack};//kGreen+2};
  //int marker[nera] = {kFullCircle, kFullDiamond,
  //		      kFullStar, kOpenStar,kFullSquare};
  int marker[nera] = {kFullCircle, kFullDiamond,
  		      kFullStar, kFullSquare, kOpenCircle};
  const char* label[nera] = {"BCD",//"CD",
			     "EF^{1}",
			     "F^{2}GH", "All","2017B"};
			     //"F^{2}G",
			     //"H", "BCDEFGH (AK4)"};


  // Reference data from the 50 ns analysis in 2015
  //TFile *fin50ns = new TFile("rootfiles/CommonPlots-74X-50ns.root","READ");
  //assert(fin50ns && !fin50ns->IsZombie());
  //curdir->cd();

  // Unfolded reference data from SMP-15-007 stored in HEPDATA
  TFile *fhepd = new TFile("rootfiles/HEPData-ins1459051-v1-root.root","READ");
  assert(fhepd && !fhepd->IsZombie());
  curdir->cd();

  // Unfolded reference data from SMP-15-007 directly from Paolo
  // InclusiveJetMeasurement_AK4_CrystalBall_JECSummer15_V5_Run2015BC_LowScale
  // InclusiveJetMeasurement_AK4_GaussianCore_JECSummer15_V5_Run2015BC_LowScale
  // InclusiveJetMeasurement_AK7_CrystalBall_JECSummer15_V5_Run2015BC_LowScale
  // InclusiveJetMeasurement_AK7_GaussianCore_JECSummer15_V5_Run2015BC_LowScale
  //TFile *f15007 = new TFile("rootfiles/InclusiveJetMeasurement_AK7_GaussianCore"//
  //			    "_JECSummer15_V5_Run2015BC_LowScale.root","READ");
  //assert(f15007 && !f15007->IsZombie());
  //curdir->cd();

  //TFile *f15007a=new TFile("rootfiles/InclusiveJetMeasurement_AK4_GaussianCore"
  //			   "_JECSummer15_V5_Run2015BC_LowScale.root","READ");
  //assert(f15007a && !f15007b->IsZombie());
  //curdir->cd();

  TH1F *h50ns4(0), *h50ns7(0);
  /*
    TH1F *h50ns;
  if (sy=="0.0-0.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_1bin");
  if (sy=="0.5-1.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_2bin");
  if (sy=="1.0-1.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_3bin");
  if (sy=="1.5-2.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_4bin");
  if (sy=="2.0-2.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_5bin");
  if (sy=="2.5-3.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_6bin");
  if (sy=="3.2-4.7") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_7bin");
  //
  if (sy=="0.0-0.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_1bin");
  if (sy=="0.5-1.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_2bin");
  if (sy=="1.0-1.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_3bin");
  if (sy=="1.5-2.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_4bin");
  if (sy=="2.0-2.5") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_5bin");
  if (sy=="2.5-3.0") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_6bin");
  if (sy=="3.2-4.7") h50ns = (TH1F*)fin50ns->Get("pt_DETInclJet_7bin");
*/

  // HEPData (1-7 AK7, 8-14 AK4; both unfolded)
  map<string, string> mhep7, mhep4;
  mhep7["0.0-0.5"] = "Table 1";
  mhep7["0.5-1.0"] = "Table 2";
  mhep7["1.0-1.5"] = "Table 3";
  mhep7["1.5-2.0"] = "Table 4";
  mhep7["2.0-2.5"] = "Table 5";
  mhep7["2.5-3.0"] = "Table 6";
  mhep7["3.2-4.7"] = "Table 7";
  //
  mhep4["0.0-0.5"] = "Table 8";
  mhep4["0.5-1.0"] = "Table 9";
  mhep4["1.0-1.5"] = "Table 10";
  mhep4["1.5-2.0"] = "Table 11";
  mhep4["2.0-2.5"] = "Table 12";
  mhep4["2.5-3.0"] = "Table 13";
  mhep4["3.2-4.7"] = "Table 14";
  //
  const char *ct7 = mhep7[sy].c_str();
  TH1F *hhepdy7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct7)); assert(hhepdy7);
  TH1F *hhepde7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct7)); assert(hhepde7);
  for (int i = 1; i != hhepdy7->GetNbinsX()+1; ++i) {
    hhepdy7->SetBinError(i, hhepde7->GetBinContent(i));
  }
  h50ns7 = hhepdy7;
  //
  const char *ct4 = mhep4[sy].c_str();
  TH1F *hhepdy4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct4)); assert(hhepdy4);
  TH1F *hhepde4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct4)); assert(hhepde4);
  for (int i = 1; i != hhepdy4->GetNbinsX()+1; ++i) {
    hhepdy4->SetBinError(i, hhepde4->GetBinContent(i));
  }
  h50ns4 = hhepdy4;

  /*
  // SMP-15-007 unfolded from Paolo
  // These contain both 'Detector' and 'Unfolded' spectra
  // but are not normalized by bin width
  // For some reason this des not seem to work properly
  map<string, string> mpg;
  mpg["0.0-0.5"] = "y0";
  mpg["0.5-1.0"] = "y1";
  mpg["1.0-1.5"] = "y2";
  mpg["1.5-2.0"] = "y3";
  mpg["2.0-2.5"] = "y4";
  mpg["2.5-3.0"] = "y5";
  mpg["3.2-4.7"] = "y6";
  const char *cr = mpg[sy].c_str();
  //TH1F *h15007 = (TH1F*)f15007->Get(Form("Detector_Spectra_%s",cr));
  TH1F *h15007 = (TH1F*)f15007->Get(Form("Unfolded_Spectra_%s",cr));
  assert(h15007);
  for (int i = 1; i != h15007->GetNbinsX()+1; ++i) {
    double Ldptdy = 100.*h15007->GetBinWidth(i);
    h15007->SetBinContent(i, h15007->GetBinContent(i) / Ldptdy);
    h15007->SetBinError(i, h15007->GetBinError(i) / Ldptdy);
  }
  // patch different pT range
  if (h15007->GetNbinsX() != h50ns->GetNbinsX()) {
    TH1F *htmp = h15007;
    h15007 = (TH1F*)h50ns->Clone(Form("hpt0_%s",cr));
    for (int i = 1; i != h15007->GetNbinsX()+1; ++i) {
      int j = htmp->FindBin(h15007->GetBinCenter(i));
      h15007->SetBinContent(i, htmp->GetBinContent(j));
      h15007->SetBinError(i, htmp->GetBinError(j));
    }
    delete htmp;
  }
  //h50ns = h15007;
  */
  
  // Settings for the spectrum fit and JEC-equivalent plot y-axis range
  double eta(0.), ymin(-4), ymax(+6); double emax(6000.); double ptmax(2600.);
  if (sy=="0.0-0.5") { eta = 0; ymin = -4; ymax = +6; }
  if (sy=="0.5-1.0") { eta = 0.5; ymin = -4; ymax = +6; }
  if (sy=="1.0-1.5") { eta = 1.0; ymin = -4; ymax = +6; }
  if (sy=="1.5-2.0") { eta = 1.5; ymin = -6; ymax = +9; }
  if (sy=="2.0-2.5") { eta = 2.0; ymin = -6; ymax = +9; }
  //if (sy=="2.5-3.0") { eta = 2.5; ymin = -12; ymax = +18; }
  if (sy=="2.5-3.0") { eta = 2.5; ymin = -15; ymax = +30; }
  if (sy=="3.2-4.7") { eta = 3.2; ymin = -12; ymax = +18; }
  assert(h50ns4);
  assert(h50ns7);

  // patch trigger biased 300 GeV point
  //int ipatch = h50ns->FindBin(320.);
  //h50ns->SetBinContent(ipatch, 0.);
  //h50ns->SetBinError(ipatch, 0.);

  // Replace 50 ns DET with HEPDATA unfolded
  /*
  if (false) {
    double p9161_d8x1y1_xval[] = { 123.5, 143.0, 163.5, 185.0, 208.0, 232.5, 258.5, 286.0, 315.0, 346.0, 378.5, 412.5, 449.0, 487.5, 527.5, 570.0, 615.0, 662.0, 711.5, 763.5, 818.0, 875.5, 936.0, 999.5, 1066.5, 1136.5, 1210.0, 1287.5, 1368.5, 1453.5, 1542.5, 1636.0, 1734.0, 1837.0, 1945.0 };
    double p9161_d8x1y1_yval[] = { 5349.33, 3064.67, 1777.82, 851.947, 426.18, 253.548, 148.523, 81.1185, 45.7049, 28.643, 17.0765, 9.79779, 6.29211, 4.13279, 2.59375, 1.59837, 1.04159, 0.657478, 0.412662, 0.264576, 0.169119, 0.109256, 0.0699437, 0.0439286, 0.0280389, 0.0196009, 0.0135767, 0.00681612, 0.00341421, 
    0.00225516, 0.00207889, 0.00176915, 0.00111131, 4.20276E-4, 1.88882E-4 };
    double p9161_d8x1y1_ystatminus[] = { 667.011, 211.904, 159.866, 97.3924, 55.0472, 15.8243, 10.3151, 5.76438, 1.50443, 0.996153, 0.631592, 0.403098, 0.179287, 0.130679, 0.0743668, 0.0486397, 0.0358208, 0.0173034, 0.0124171, 0.00809688, 0.0062119, 0.00484187, 0.00379823, 0.00291624, 0.00229541, 0.00183255, 0.0015417, 0.0010303, 7.29814E-4, 5.59543E-4, 5.53858E-4, 4.72921E-4, 3.88821E-4, 2.14123E-4, 1.49492E-4 };
    const int npoints = 35;
    
    h50ns->Clear();
    for (int i = 0; i != npoints; ++i) {
      int j = h50ns->FindBin(p9161_d8x1y1_xval[i]);
      h50ns->SetBinContent(j, p9161_d8x1y1_yval[i]);
      h50ns->SetBinError(j, p9161_d8x1y1_ystatminus[i]);
    }
  } // if true
  */
  
  //////////////////////////////////////////////////////////////////////////
  
  // Plot ECAL prefire inefficiency
  TH1D *h0 = new TH1D("h0",";Jet p_{T} (GeV);"
		      "ECAL prefire fraction",3445,55,3500);
  h0->SetMinimum(0.);
  h0->SetMaximum(1.);
  h0->GetXaxis()->SetNoExponent();
  h0->GetXaxis()->SetMoreLogLabels();

  TH1D *h02 = new TH1D("h02",";Jet p_{T} (GeV);"
		       "Ratio to fit",3445,55,3500);
  h02->SetMinimum(0.);//0.5);//0.8);
  h02->SetMaximum(2.0);//1.5);//1.2);
  h02->GetXaxis()->SetNoExponent();
  h02->GetXaxis()->SetMoreLogLabels();

  lumi_13TeV = "Run2016, 35 fb^{-1}";
  //TCanvas *c0 = tdrDiCanvas("c10",h0,h02,4,11);
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);

  //c0->cd(1);
  gPad->SetLogx();

  int ieta1 = heff->GetYaxis()->FindBin(eta+0.1);
  int ieta2 = heff->GetYaxis()->FindBin(eta+0.4);
  TH1D *heffeta1 = heff->ProjectionX("heffeta1",ieta1,ieta1,"e");
  TH1D *heffeta2 = heff->ProjectionX("heffeta2",ieta2,ieta2,"e");
  TH1D *heffeta12 = heff->ProjectionX("heffeta12",ieta1,ieta2,"e");
  heffeta12->Scale(1./(ieta2+1-ieta1));
  //
  int ietaRB = heff2->GetYaxis()->FindBin(eta+0.1);
  TH1D *heffetaRB = heff2->ProjectionX("heffetaRB",ietaRB,ietaRB,"e");
  TH1D *heffetaBCD = heff2bcd->ProjectionX("heffetaBCD",ietaRB,ietaRB,"e");
  TH1D *heffetaEF = heff2ef->ProjectionX("heffetaEF",ietaRB,ietaRB,"e");
  TH1D *heffetaGH = heff2gh->ProjectionX("heffetaGH",ietaRB,ietaRB,"e");
  
  // Attempt to properly weight the two eta bins using d2sigma/dpt/dy
  TF2 *f2 = new TF2("f2","[0]*pow(x,[1])"
                    "*pow(1-min(2*x*cosh(y),12998.)/13000.,[2]-2*y)"
                    "*(2*x*cosh(y)<12998.)",
                    60,4000.,0,5.2);
  f2->SetParameters(1.15e14,-5,10);
  TH1D *heffeta = (TH1D*)heffeta12->Clone("heffeta"); heffeta->Clear();
  TGraphErrors *geffeta = new TGraphErrors(0);
  TGraphErrors *geffetaRB = new TGraphErrors(0);
  TGraphErrors *geffetaBCD = new TGraphErrors(0);
  TGraphErrors *geffetaEF = new TGraphErrors(0);
  TGraphErrors *geffetaGH = new TGraphErrors(0);
  for (int i = 1; i != heffeta->GetNbinsX()+1; ++i) {
    double pt1 = heffeta->GetBinLowEdge(i);
    double pt2 = heffeta->GetBinLowEdge(i+1);
    double lumi = 35000*0.002; // /pb x unpre-fireable efficiency (0.2%)
    double n1 = f2->Integral(pt1, pt2, eta, eta+0.25)*lumi;
    double n2 = f2->Integral(pt1, pt2, eta+0.25, eta+0.5)*lumi;
    double eff1 = heffeta1->GetBinContent(i);
    double eff2 = heffeta2->GetBinContent(i);

    double eff = (n1+n2!=0 ? (eff1*n1 + eff2*n2)/(n1+n2) : 0);
    double eff_e = (n1+n2!=0 ? eff*(1-eff)/sqrt(n1+n2) : 0);
    heffeta->SetBinContent(i, eff);
    heffeta->SetBinError(i, eff_e);

    if (eff!=0 && eff_e!=0) {
      int n = geffeta->GetN();
      double mine = (eta<2.5 ? 0.003 : 0.01);
      geffeta->SetPoint(n, 0.5*(pt1+pt2), eff);
      geffeta->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
			     sqrt(mine*mine+eff_e*eff_e));
      //
      geffetaRB->SetPoint(n, 0.5*(pt1+pt2), heffetaRB->GetBinContent(i));
      geffetaRB->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
			       sqrt(mine*mine+eff_e*eff_e));
      //
      geffetaBCD->SetPoint(n, 0.5*(pt1+pt2), heffetaBCD->GetBinContent(i));
      geffetaBCD->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
				sqrt(mine*mine+eff_e*eff_e));
      geffetaEF->SetPoint(n, 0.5*(pt1+pt2), heffetaEF->GetBinContent(i));
      geffetaEF->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
			       sqrt(mine*mine+eff_e*eff_e));
      geffetaGH->SetPoint(n, 0.5*(pt1+pt2), heffetaGH->GetBinContent(i));
      geffetaGH->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
			       sqrt(mine*mine+eff_e*eff_e));
    }
  }

  TF1 *feff = new TF1("feff","[0]*0.5*(1+erf((x-[1])/([2]*sqrt(x))))",
		      55,4000./cosh(eta));
  TF1 *feffBCD = new TF1("feffBCD","[0]*0.5*(1+erf((x-[1])/([2]*sqrt(x))))",
                         55,4000./cosh(eta));
  TF1 *feffEF = new TF1("feffEF","[0]*0.5*(1+erf((x-[1])/([2]*sqrt(x))))",
		      55,4000./cosh(eta));
  TF1 *feffGH = new TF1("feffGH","[0]*0.5*(1+erf((x-[1])/([2]*sqrt(x))))",
		      55,4000./cosh(eta));
  feff->SetParameters(eta<2.5 ? 0.1 : 0.6, 180, 6);
  //heffeta->Fit(feff,"QRNW");
  //geffeta->Fit(feff,"QRNW");
  //geffeta->Fit(feff,"QRN"); feff->SetLineColor(kRed);
  geffetaRB->Fit(feff,"QRN"); feff->SetLineColor(kBlue);
  geffetaBCD->Fit(feffBCD,"QRN"); feffBCD->SetLineColor(kRed+2);
  geffetaEF->Fit(feffEF,"QRN"); feffEF->SetLineColor(kOrange+2);
  geffetaGH->Fit(feffGH,"QRN"); feffGH->SetLineColor(kBlue+1);
  
  tdrDraw(heffeta,"z",kFullSquare,kBlack,kSolid,kBlack);
  tdrDraw(geffeta,"Pz",kFullSquare,kBlack,kSolid,kBlack);
  tdrDraw(heffetaRB,"HIST Pz",kFullDiamond,kBlue,kSolid,kBlue);
  tdrDraw(geffetaRB,"Pz",kFullDiamond,kBlue,kSolid,kBlue);
  tdrDraw(heffeta12,"HIST Pz",kOpenSquare,kBlack,kSolid,kBlack);
  tdrDraw(heffeta1,"HIST PLz",kOpenTriangleUp,kBlack,kSolid,kBlack,kNone);
  tdrDraw(heffeta2,"HIST PLz",kOpenTriangleDown,kBlack,kSolid,kBlack,kNone);
  feff->Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kBlue);//kRed);
  double d = (eta<2.5 ? 0.18 : 0);
  tex->DrawLatex(0.50,0.40+d,"#frac{C}{2}#times#left(1 + "
		 "erf#left(#frac{p_{T}-#mu}{S#sqrt{p_{T}}}#right)#right)");
  tex->DrawLatex(0.50,0.30+d,Form("#chi^{2} / NDF = %1.1f / %d",
				  feff->GetChisquare(), feff->GetNDF()));
  tex->DrawLatex(0.50,0.25+d,Form("C = %1.4f #pm %1.4f",
				  feff->GetParameter(0), feff->GetParError(0)));
  tex->DrawLatex(0.50,0.20+d,Form("#mu = %1.1f #pm %1.1f GeV",
				  feff->GetParameter(1), feff->GetParError(1)));
  tex->DrawLatex(0.50,0.15+d,Form("S = %1.2f #pm %1.2f",
				  feff->GetParameter(2), feff->GetParError(2)));


  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.65,0.60,Form("%1.1f#leq|y|<%1.1f",eta,eta+0.5));
  tex->DrawLatex(0.19,0.75,Form("%1.1f#leq|y|<%1.1f",eta,eta+0.5));
  //
  
  TLegend *leg0 = tdrLeg(0.62,0.65,0.82,0.90);
  leg0->AddEntry(heffeta2,Form("%1.2f#leq|y|<%1.2f",eta,eta+0.25),"PL");
  leg0->AddEntry(heffeta1,Form("%1.2f#leq|y|<%1.2f",eta+0.25,eta+0.5),"PL");
  leg0->AddEntry(heffeta12,"direct avg.","P");
  leg0->AddEntry(heffeta,"weighted avg.","P");
  leg0->AddEntry(heffetaRB,"rebinned","P");

  //C0(2);
  //gPad->SetLogx();
  
  c0->SaveAs(Form("pdf/drawJetPt_ECAPprefire_%s.pdf",cy2));
  
  //////////////////////////////////////////////

  
  // Load 2016 results
  TH1F *h1s[nera];
  for (int iera = 0; iera != nera; ++iera) {

    const char *cera = eras[iera].c_str(); 
    //string hname = Form("ak4/y_%s/hpt_Run%s_rereco",cy,cera);
    //string hname = Form("ak7/y_%s/hptRun%s_Febr2017V3",cy,cera);
    //string hname = Form("ak7/y_%s/hptRun%s_Febr2017V3_unfolded",cy,cera);
    //string hname = Form("ak4/y_%s/hptData_%s2016_particle_%dbin",cy,cera,
    string hname = Form("ak4/y_%s/hptData_%s2016_detector_%dbin",cy,cera,
			eta==3.2 ? 7 : int(eta/0.5)+1);
    // PATCH start
    //if (string(cera)=="FlateG") hname = Form("ak7/y_%s/hptRunG_Febr2017V3",cy);
    if (eras[iera]=="BCDEFGH")
      //hname = Form("ak4/y_%s/hptData_%dbin",cy,
      //	   eta==3.2 ? 7 : int(eta/0.5)+1);
      //hname = Form("ak4/y_%s/hptData_full2016_particle_%dbin",cy,
      hname = Form("ak4/y_%s/hptData_full2016_detector_%dbin",cy,
		   eta==3.2 ? 7 : int(eta/0.5)+1);
    if (eras[iera]=="B")
      hname = Form("ak4/Eta_%s/hpt_data_2017_%s_det",cy,eras[iera].c_str());
    // PATCH end
    TH1F *hera = (TH1F*)fins[iera]->Get(hname.c_str());
    cout << hname << endl << flush;
    assert(hera);
    //hera->Scale(1.033); // for luminosity patch for 2016

    //if (eras[iera]=="B") {
    if (true) {
      TF1 *f1 = (TF1*)fu->Get(Form("fr_%s",cy)); assert(f1);
      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	if (hera->GetBinContent(i)!=0) {
	  double pt = hera->GetBinCenter(i);
	  double rdu = f1->Eval(pt);
	  hera->SetBinContent(i, hera->GetBinContent(i)/rdu);
	  hera->SetBinError(i, hera->GetBinError(i)/rdu);
	}
      }
    }

    // patch different pT binning
    /*
    if (hera->GetNbinsX() != h50ns4->GetNbinsX()) {
      TH1F *htmp = hera;
      hera = (TH1F*)h50ns4->Clone(Form("hpt0_%d",iera));
      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	int j = htmp->FindBin(hera->GetBinCenter(i));
	hera->SetBinContent(i, htmp->GetBinContent(j));
	hera->SetBinError(i, htmp->GetBinError(j));
      }
    }
    */
    hera->Scale(1./lumi[iera]);
    h1s[iera] = hera;

    // patch ECAL trigger inefficiency
    for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
      double pt = hera->GetBinCenter(i);
      int ix = heff->FindBin(pt, eta+0.1);
      //double ineff = heff->GetBinContent(ix);
      //double ineff = feff->Eval(pt);
      // Calculate weighted inefficiency for wide pT bin
      double xsec(0), xsecECAL(0);
      double pt1 = hera->GetBinLowEdge(i);
      double pt2 = hera->GetBinLowEdge(i+1);
      int nx = 9; // odd number
      double dx = (pt2-pt1)/nx;
      for (int j = 0; j != nx; ++j) {
	double x = pt1 + dx * (j+0.5);
	double xs = pow(x,-5)*pow(1-2*x*cosh(eta)/13000.,10) * dx;
	xsec += xs;
	double ineff = feff->Eval(x);
	if (doIOVeff) {
	  if (eras[iera]=="BCDEFGH") ineff = feff->Eval(x);
	  if (eras[iera]=="BCD") ineff = feffBCD->Eval(x);
	  if (eras[iera]=="EF") ineff = feffEF->Eval(x);
	  if (eras[iera]=="GH") ineff = feffGH->Eval(x);
	  if (eras[iera]=="B") ineff = feffGH->Eval(x);
	}
	xsecECAL += (1-ineff) * xs;
      }
      double ineff = 1 - xsecECAL / xsec;
      if (eta<2 || eta>3) ineff = 0;
      
      if (ineff != 1) {
	hera->SetBinContent(i, hera->GetBinContent(i)/(1-ineff));
	hera->SetBinError(i, hera->GetBinError(i)/(1-ineff));
      }
    }

    //if (eras[iera]=="BCDEFGH")
    //hera = rebinXsec(hera);
  } // for iera

  TF1 *f1s[nera];
  for (int iera = 0; iera != nera; ++iera) {
    TF1 *f1 = new TF1(Form("f1_%d",iera),
		      Form("[0]*pow(x,[1]+[2]*log(x))"
			   "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
		      ptmin,emax/cosh(eta));
    f1->SetParameters(1e12,-5,0.001,10);
    f1->SetNpx(2640.*2); // this limits precision of DeltaJEC
    h1s[iera]->Fit(f1,"RNI");
    f1->SetLineColor(color[iera]);
    f1s[iera] = f1;
  } // for iera

  TF1 *f50ns7 = new TF1("f50ns7",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns7->SetParameters(1e12,-5,0.001,10);
  f50ns7->SetNpx(2640.*2); // this limits precision of DeltaJEC
  h50ns7->Fit(f50ns7,"RNI");
  f50ns7->SetLineColor(kGreen+2);

  TF1 *f50ns4 = new TF1("f50ns4",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns4->SetParameters(1e12,-5,0.001,10);
  f50ns4->SetNpx(2640.*2); // this limits precision of DeltaJEC
  h50ns4->Fit(f50ns4,"RNI");
  f50ns4->SetLineColor(kGreen+2);

  TF1 *fref = f50ns7; // compare to 50 ns data

  // Calculate ratio of fits
  TF1 *f1rs[nera];
  for (int iera = 0; iera != nera; ++iera) {

    //fref = (eras[iera]=="BCDEFGH" ? f50ns4 : f50ns7);
    fref = f50ns4;

    TF1 *f1 = f1s[iera];
    TF1 *f1r = new TF1(Form("f1r_%d",iera),
		       Form("[0]*pow(x,[1]+[2]*log(x))"
			    "*pow(1-2.*x*cosh(%1.2f)/13000.,[3]) / "
			    "([4]*pow(x,[5]+[6]*log(x))"
			    "*pow(1-2.*x*cosh(%1.2f)/13000.,[7]))",
			    eta,eta),
		       ptmin,emax/cosh(eta));
    f1r->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		       f1->GetParameter(2),f1->GetParameter(3),
		       fref->GetParameter(0),fref->GetParameter(1),
		       fref->GetParameter(2),fref->GetParameter(3));
    f1r->SetLineColor(color[iera]);
    f1rs[iera] = f1r;
  } // for iera

  // Calculate ratio of fit (nominally should be 1 when 50 ns is fref)
  TF1 *f50nsr = new TF1("f50nsr",
			Form("[0]*pow(x,[1]+[2]*log(x))"
		       "*pow(1-2.*x*cosh(%1.2f)/13000.,[3]) / "
			     "([4]*pow(x,[5]+[6]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[7]))",
			     eta,eta),
			ptmin,emax/cosh(eta));
  f50nsr->SetParameters(f50ns4->GetParameter(0),f50ns4->GetParameter(1),
			f50ns4->GetParameter(2),f50ns4->GetParameter(3),
			fref->GetParameter(0),fref->GetParameter(1),
			fref->GetParameter(2),fref->GetParameter(3));
  f50nsr->SetLineColor(kGreen+2);

  // Divide 2016 data by reference fit
  TH1F *hrs[nera];
  for (int iera = 0; iera != nera; ++iera) {
    //fref = (eras[iera]=="BCDEFGH" ? f50ns4 : f50ns7);
    fref = f50ns4;
    TH1F *hr = (TH1F*)h1s[iera]->Clone(Form("hr_%d",iera));
    hr->Divide(fref);
    hrs[iera] = hr;
  } // for iera

  // Divide 2016 data by BCDEFGH for direct time stability
  TH1F *h1rs[nera];
  for (int iera = 1; iera != nera; ++iera) {
    TH1F *h1r = (TH1F*)h1s[iera]->Clone(Form("h1r_%d",iera));
    h1r->Divide(h1s[0]);
    h1rs[iera] = h1r;
  } // for iera

  // Ratio of 50 ns data to reference (50 ns) fit
  TH1F *hr50ns4 = (TH1F*)h50ns4->Clone("hr50ns4");
  hr50ns4->Divide(f50ns4);//fref);

  TH1F *hr50ns7 = (TH1F*)h50ns7->Clone("hr50ns7");
  hr50ns7->Divide(f50ns7);//fref);

  
  TH1D *h = new TH1D("h",";Jet p_{T} (GeV);"
		     "d#sigma^{2} / dy dp_{T} (pb / GeV)",3445,55,3500);
  h->SetMinimum(1e-6 *1.0001);
  h->SetMaximum(1e+8 *0.9999);

  TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);"
		     "Ratio to 50 ns",3445,55,3500);
  h2->SetMinimum(0.);//0.5);//0.8);
  h2->SetMaximum(2.0);//1.5);//1.2);
  h2->GetXaxis()->SetNoExponent();
  h2->GetXaxis()->SetMoreLogLabels();

  lumi_13TeV = "Run2016, 35 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,11);

  c1->cd(1);
  tdrDraw(h50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  tdrDraw(h50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    tdrDraw(h1s[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
  }
  gPad->SetLogx();
  gPad->SetLogy();
  
  TLegend *leg = tdrLeg(0.50,0.90-0.06*nera,0.80,0.90);
  for (int iera = 0; iera != nera; ++iera) {
    leg->AddEntry(h1s[iera],label[iera],"PL");
  }
  leg->AddEntry(h50ns4,"74X 50 ns AK4","PL");
  leg->AddEntry(h50ns7,"74X 50 ns AK7","PL");

  c1->cd(2);
  tdrDraw(hr50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  tdrDraw(hr50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    tdrDraw(hrs[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
  }
  gPad->SetLogx();

  c1->cd(1);
  //f50ns->SetRange(60,3500);
  f50ns4->SetRange(ptmin,3500);
  f50ns4->Draw("SAME");
  f50ns7->SetRange(ptmin,3500);
  f50ns7->Draw("SAME");
  // Draw 2016 fits as well (or don't)
  for (int iera = 0; iera != nera; ++iera) {
    //f1s[iera]->Draw("SAME");
  }

  //TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.20,0.10,Form("fit(50ns) AK7: #chi^{2}/NDF = %1.1f/%d",
				f50ns7->GetChisquare(),f50ns7->GetNDF()));
  tex->DrawLatex(0.20,0.05,Form("fit(50ns) AK4: #chi^{2}/NDF = %1.1f/%d",
				f50ns4->GetChisquare(),f50ns4->GetNDF()));
						
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.20,0.20,Form("y #in %s",cy)); 

  c1->cd(2);
  f50nsr->Draw("SAME");
  // Draw ratio of fits as well (or don't)
  for (int iera = 0; iera != nera; ++iera) {
    //f1rs[iera]->Draw("SAME");
  }

  c1->SaveAs(Form("pdf/drawJetPt_%s.pdf",cy2));


  TH1D *h3 = new TH1D("h3",";Jet p_{T} (GeV);#DeltaJEC-equivalent (%)",
		      3445,55,3500);
  h3->SetMinimum(ymin);
  h3->SetMaximum(ymax);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  
  // Calculate reference uncertainty
  TH1F *hunc = (TH1F*)h50ns4->Clone("hunc");
  TH1F *huncs1 = (TH1F*)h50ns4->Clone("huncs1");
  TH1F *huncs2 = (TH1F*)h50ns4->Clone("huncs2");
  TH1F *hunc0 = (TH1F*)h50ns4->Clone("hunc0");
  TH1F *hunc1 = (TH1F*)h50ns4->Clone("hunc1");
  for (int i = 1; i != hunc->GetNbinsX()+1; ++i) {
    
    double pt = hunc->GetBinCenter(i);

    unc->setJetEta(eta+0.25);
    unc->setJetPt(pt);
    hunc->SetBinContent(i, 0);//-1.5);
    hunc->SetBinError(i, 100.*unc->getUncertainty(true));

    uncs1->setJetEta(eta+0.25);
    uncs1->setJetPt(pt);
    huncs1->SetBinContent(i, 100.*uncs1->getUncertainty(true));
    huncs1->SetBinError(i, 0.);

    uncs2a->setJetEta(eta+0.25);
    uncs2a->setJetPt(pt);
    uncs2b->setJetEta(eta+0.25);
    uncs2b->setJetPt(pt);
    huncs2->SetBinContent(i, 100.*2.*uncs2a->getUncertainty(true) +
			  100.*2.*uncs2b->getUncertainty(true));
    huncs2->SetBinError(i, 0.);
    
    unc0->setJetEta(eta+0.25);
    unc0->setJetPt(pt);
    hunc0->SetBinContent(i, 0.0);
    hunc0->SetBinError(i, 100.*unc0->getUncertainty(true));

    unc1->setJetEta(eta+0.25);
    unc1->setJetPt(pt);
    hunc1->SetBinContent(i, 0.0);
    hunc1->SetBinError(i, 100.*unc1->getUncertainty(true));
  }

  TH1F *hds[nera];
  TH1F *huncl = (TH1F*)h50ns4->Clone("huncl");
  for (int iera = 0; iera != nera; ++iera) {
    
    //fref = (eras[iera]=="BCDEFGH" ? f50ns4 : f50ns7);
    fref = f50ns4;
    TH1F *hr = hrs[iera];
    TH1F *hd = (TH1F*)hr->Clone(Form("hd_%d",iera));
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {

      if (hr->GetBinContent(i)==0) continue;
      double x = hr->GetBinCenter(i);
      double y = fref->Eval(x);
      double y2 = y * hr->GetBinContent(i);
      double x2 = fref->GetX(y2,0.85*x,1.15*x,1e-3);
      hd->SetBinContent(i, 100.*(x/x2-1));

      // Luminosity uncertainty from BCDEFGH
      if (iera==0) {
	double y3 = y / 1.027; // lumi uncertainty for 2016, after 3.3% patch
	double x3 = fref->GetX(y3,0.85*x,1.15*x,1e-3);
	int j = huncl->FindBin(hd->GetBinCenter(i));
	huncl->SetBinContent(j, 0);
	huncl->SetBinError(j, 100.*fabs(x/x3-1));
      }

      // 2016 uncertainty around first entry (BCDEFGH)
      //if (iera==0) hunc->SetBinContent(hunc->FindBin(hd->GetBinCenter(i)),
      //			       hd->GetBinContent(i));

      double y2_up = y * (hr->GetBinContent(i) + hr->GetBinError(i));
      double x2_dw = fref->GetX(y2_up,0.85*x,1.15*x);
      double y2_dw = y * (hr->GetBinContent(i) - hr->GetBinError(i));
      double x2_up = fref->GetX(y2_dw,0.85*x,1.15*x);
      hd->SetBinError(i, 100.*sqrt(pow((x2_up-x2_dw)/x2,2) + pow(0.001,2)));
    } // for i
    hds[iera] = hd;

  } // for iera


  TCanvas *c2 = tdrCanvas("c2",h3,4,11,kSquare);

  tdrDraw(hunc0,"E3", kNone, kBlack,kSolid,-1,1001,kOrange-9);
  tdrDraw(hunc,"E3", kSolid, kBlack,kSolid,-1,1001,kYellow+1);
  tdrDraw(hunc1,"E3", kNone, kBlack,kSolid,-1,1001,kRed-9);
  //tdrDraw(huncl,"E3", kNone, kBlack,kSolid,-1,1001,kBlue-9);
  //tdrDraw(hunc,"E3", kSolid, kBlack,kSolid,-1,3001,kYellow+1);
  tdrDraw(huncs1,"HIST][", kSolid, kBlack, kSolid,-1,kNone,kOrange-9);
  tdrDraw(huncs2,"HIST][", kSolid, kBlue, kSolid,-1,kNone,kOrange-9);
  
  for  (int iera = 0; iera != nera; ++iera) {
    tdrDraw(hds[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
  }
  gPad->SetLogx();

  //TLegend *legd = tdrLeg(0.50,0.90-0.06*nera,0.80,0.90);
  TLegend *legd = tdrLeg(0.50,0.90-0.04*(nera+3),0.80,0.90);
  for  (int iera = 0; iera != nera; ++iera) {
    legd->AddEntry(hrs[iera],label[iera],"PL");
  }

  TLegend *legu = tdrLeg(0.50,0.15,0.80,0.23);
  //legd->AddEntry(hunc0,"2015 (50 ns) unc.","F");
  //legd->AddEntry(hunc,"2016 (25 ns) unc.","F");
  legd->AddEntry(hunc1,"2012 unc.","F");
  legd->AddEntry(hunc,"2016 unc.","F");
  legd->AddEntry(hunc0,"2015 unc.","F");
  //legd->AddEntry(huncl,"2.6% lum.","F");
  
  tex->SetTextColor(kBlack); tex->SetTextSize(0.040);
  //tex->DrawLatex(0.19,0.74,Form("Anti-k_{T} R=0.4"));
  tex->DrawLatex(0.19,0.74,Form("Anti-k_{T} R=0.7 (0.4)"));
  //tex->DrawLatex(0.19,0.69,Form("PF+CHS"));
  tex->DrawLatex(0.19,0.69,Form("PF+CHS unf."));
  tex->DrawLatex(0.19,0.64,Form("|y|#in %s",cy));
  tex->SetTextSize(0.025);
  tex->DrawLatex(0.19,0.59,Form("vs HEPData-ins1459051-v1"));
  tex->SetTextSize(0.040);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(60,0,3500,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(60,+1.5,3500,+1.5);
  l->DrawLine(60,-1.5,3500,-1.5);

  c2->RedrawAxis();
  c2->SaveAs(Form("pdf/drawJetPt_DeltaJEC_%s.pdf",cy2));

  TH1F *he = hds[0]; // BCDEFGH
  TF1 *fg = new TF1("fg","[0]+[1]*log(x)+"
		    "[2]*TMath::Gaus(x,[3],[4]*sqrt(x))",
		    60,2500);
  fg->SetParameters(-1.5,+0.2, -1.5, 400, 15.);
  he->Fit(fg, "RN");
  fg->DrawClone("SAME");
  fg->SetLineStyle(kDashed);
  double p0 = fg->GetParameter(0);
  double p1 = fg->GetParameter(1);
  double ecal = fg->GetParameter(2);
  fg->SetParameter(2,0);
  fg->DrawClone("SAME");
  fg->SetParameters(0,0,ecal,fg->GetParameter(3),fg->GetParameter(4));
  fg->SetLineStyle(kDotted);
  fg->DrawClone("SAME");
  
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d",
				fg->GetChisquare(),fg->GetNDF()));
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.20,0.260,Form("p_{0} = %1.2f#pm%1.2f%%", p0,
				 fg->GetParError(0)));
  tex->DrawLatex(0.20,0.225,Form("p_{1} = %1.2f#pm%1.2f%%", p1,
				 fg->GetParError(1)));

  tex->DrawLatex(0.60,0.225,Form("#Delta = %1.2f#pm%1.2f%%",
				fg->GetParameter(2),
				fg->GetParError(2)));
  tex->DrawLatex(0.60,0.190,Form("#mu = %1.0f#pm%1.0f GeV",
				fg->GetParameter(3),
				fg->GetParError(3)));
  tex->DrawLatex(0.60,0.155,Form("#sigma = %1.1f#pm%1.1f"
				 "#sqrt{GeV}",
				fg->GetParameter(4),
				fg->GetParError(4)));

  c2->SaveAs(Form("pdf/drawJetPt_DeltaJEC2_%s.pdf",cy2));


  TH1D *h4 = new TH1D("h4",";Jet p_{T} (GeV);RunX / RunBCDEFGH",3445,55,3500);
  bool isec2 = (fabs(eta)>2.2);
  h4->SetMinimum(isec2 ? 0.60 : 0.90);//0.80);
  h4->SetMaximum(isec2 ? 1.80 : 1.20);//1.40);
  h4->GetXaxis()->SetNoExponent();
  h4->GetXaxis()->SetMoreLogLabels();

  TCanvas *c3 = tdrCanvas("c3",h4,4,11,kSquare);
  gPad->SetLogx();
  for (int iera = 1; iera != nera; ++iera) {
    tdrDraw(h1rs[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
  }

  TLegend *leg3 = tdrLeg(0.50,0.90-0.06*nera,0.70,0.90);
  for  (int iera = 1; iera != nera; ++iera) {
    leg3->AddEntry(h1rs[iera],label[iera],"PL");
  }

  tex->SetTextColor(kBlack); tex->SetTextSize(0.040);
  tex->DrawLatex(0.19,0.74,Form("Anti-k_{T} R=0.4"));
  tex->DrawLatex(0.19,0.69,Form("PF+CHS"));
  tex->DrawLatex(0.19,0.64,Form("|y| #in %s",cy));

  
  c3->SaveAs(Form("pdf/drawJetPt_XsecoverBCDEFGH_%s.pdf",cy2));

}


// Tool to estimate unfolding corrections
void unfold(string sy = "0.0-0.5", int ieta = 1) {

  TFile *f = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *f = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  assert(f && !f->IsZombie());
  TFile *fout = new TFile("rootfiles/unfold.root",
			  ieta==1 ? "RECREATE" : "UPDATE");

  const char *cy = sy.c_str();
  TH1D *hd = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_detector_%dbin",
				cy,ieta));
  assert(hd);
  TH1D *hu = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_particle_%dbin",
				cy,ieta));
  assert(hu);

  fout->cd();
  TH1D *hr = (TH1D*)hu->Clone(Form("hr_%s",cy));
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    int i1 = hd->FindBin(hr->GetBinLowEdge(i));
    int i2 = hd->FindBin(hr->GetBinLowEdge(i+1)-0.5);
    double yd = (hd->GetBinContent(i1)*hd->GetBinWidth(i1) +
		 hd->GetBinContent(i2)*hd->GetBinWidth(i2)) /
      (hd->GetBinWidth(i1) + hd->GetBinWidth(i2));
    double yu = hu->GetBinContent(i);
    if (yu!=0 && yd!=0) {
      hr->SetBinContent(i, yd/yu);
      hr->SetBinError(i, hu->GetBinError(i)/hu->GetBinContent(i)
		      * hr->GetBinContent(i));
    }
    else {
      hr->SetBinContent(i, 0);
      hr->SetBinError(i, 0);
    }
  }

  //TF1 *f1 = new TF1(Form("fr_%d",ieta),"[0]+[1]*log(x)+[2]*log(x)*log(x)",
  TF1 *f1 = new TF1(Form("fr_%s",cy),"[0]+[1]*log(x)+[2]*log(x)*log(x)",
		    114,3000);
  f1->SetParameters(1,0.01,-0.001);
  hr->Fit(f1,"QRN");
  f1->Write();

  fout->Write();
  fout->Close();
} // unfold
