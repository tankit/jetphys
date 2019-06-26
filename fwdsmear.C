
// Purpose: d'Agostini ("Bayesian" or Richardson-Lucy) unfolding, including
//          response matrix generation from NLO theory and parameterized JER
// Author:  mikko.voutilainen@cern.ch
// Created: September 2, 2012
// Updated: June 5, 2015
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"

#include "RooUnfold/src/RooUnfold.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldBinByBin.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
//#include "RooUnfold.h"

#include "tdrstyle_mod18.C"
#include "ptresolution.h"
#include "settings.h"
#include "tools.h"

#include <iostream>

using namespace std;

// Resolution function
int _jk = 0; // global variable
bool _jet = false; // global variable

Double_t fPtRes(Double_t *x, Double_t *p) {

   return ptresolution(x[0], p[0]);
   // return effectiveJER(x[0], p[0]);
}


// Ansatz Kernel 
int cnt_a = 0;
//const int nk = 4; // number of kernel parameters (excluding pt, eta)
const int nk = 3; // number of kernel parameters (excluding pt, eta)
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) {
    cout << "+" << flush;
  }

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[1]; // rapidity

  //cout << Form("true pt: %10.5f p[0]: %10.5f  %10.5f  %10.5f  %10.5f  %10.5f ", pt, p[0], p[1], p[2], p[3], p[4]) << endl;
  //  double resold = ptresolution(pt, eta+1e-3) * pt;
  //  double res = effectiveJER(pt, eta) * pt;    // Switching this on breaks the integration if the kernel. Resolutions look okay
  double res = ptresolution(pt, eta) * pt;

  const double s = TMath::Gaus(ptmeas, pt, res, kTRUE);
 
  //   const double f = p[2] * exp(p[3]/pt) * pow(pt, p[4])
  //    * pow(1 - pt*cosh(eta) / jp::emax, p[5]);
  const double f = p[2] * pow(pt, p[3])
    * pow(1 - pt*cosh(eta) / jp::emax, p[4]);
  
  return (f * s);
}

// Smeared Ansatzz
// double _epsilon = 1e-12; // Tolerance too small when using effective JER...
 double _epsilon = 1e-6;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  const double pt = x[0];
  const double eta = p[0];
  //const double eta = 0.0;

  if (!_kernel) _kernel = new TF1("_kernel", smearedAnsatzKernel,
				  1., jp::emax/cosh(eta), nk+2);

  double res = ptresolution(pt, eta+1e-3) * pt;
  //  double res = effectiveJER(pt, eta) * pt;
  
  const double sigma = max(0.10, min(res/pt, 0.30));
  double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  ptmin = max(1.,ptmin); // safety check
  double ptmax = pt / (1. - 3.*sigma); // xmax*(1-3*sigma)=x
  //  cout << Form("1pt %10.5f sigma %10.5f ptmin %10.5f ptmax %10.5f eta %10.5f",pt, sigma, ptmin, ptmax, eta) << endl << flush;
  ptmax = min(jp::emax/cosh(eta), ptmax); // safety check
  // cout << Form("2pt %10.5f sigma %10.5f ptmin %10.5f ptmax %10.5f eta %10.5f",pt, sigma, ptmin, ptmax, eta) << endl << flush;
  if(ptmin>ptmax)return 1.0;//ptmax=ptmin;

  const double par[nk+2] = {pt, eta, p[1], p[2], p[3]};
  _kernel->SetParameters(&par[0]);

  // Set pT bin limits needed in smearing matrix generation
  if (p[4]>0 && p[4]<jp::emax/cosh(eta)) ptmin = p[4];
  if (p[5]>0 && p[5]<jp::emax/cosh(eta)) ptmax = p[5];
				            
   return ( _kernel->Integral(ptmin, ptmax, _epsilon) );
   //  return ( 1.0); // integral fails due to nan ptmin ptmax
}
void recurseFile(TDirectory *indir, TDirectory *indir2, TDirectory *outdir,
                 bool ismc);
void Unfold_histo(TH1D *hpt, TH1D *hpt2, TDirectory *outdir,
			   bool ismc, bool kscale = false, string id = "");


void fwdsmearUnfold(string type) {

  TFile *fin = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  // TFile *fin2 = new TFile(Form("output-%s-2c.root",type.c_str()),"READ");
  // TFile *fin2 = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  // TFile *fin2 = new TFile(Form("output-%s-2c.root","MC"),"READ");
  TFile *fin2 = new TFile(jp::dagfile1 ? "output-MC-1.root" : "output-MC-2b.root","READ");

  assert(fin2 && !fin2->IsZombie());

  TFile *fout = new TFile(Form("output-%s-3.root",type.c_str()),"RECREATE");
  assert(fout && !fout->IsZombie());

    bool ismc = jp::ismc;

  recurseFile(fin, fin2, fout, ismc);

  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();

  fin->Close();
  fin->Delete();
}


void recurseFile(TDirectory *indir, TDirectory *indir2, TDirectory *outdir,
                 bool ismc) {

  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory: copy it to output and go deeper
    if (obj->InheritsFrom("TDirectory")) {

      if (jp::debug) cout << key->GetName() << endl;

      assert(outdir->mkdir(obj->GetName()));
      outdir->mkdir(obj->GetName());
      assert(outdir->cd(obj->GetName()));
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName()); assert(outdir2);
      outdir2->cd();

      assert(indir->cd(obj->GetName()));
      TDirectory *indir2a = indir->GetDirectory(obj->GetName()); assert(indir2a);
      indir2a->cd();

      if (indir2->cd(obj->GetName())) {
        TDirectory *indir2b = indir2->GetDirectory(obj->GetName()); assert(indir2b);

        recurseFile(indir2a, indir2b, outdir2, ismc);
      }
    } // inherits from TDirectory

    // Found hpt plot: call unfolding routine
    if (obj->InheritsFrom("TH1") &&
        (string(obj->GetName())=="hpt" /* ||
	 string(obj->GetName())=="hpt_jet" ||
	 string(obj->GetName())=="hpt_jk1" ||
	 string(obj->GetName())=="hpt_jk2" ||
	 string(obj->GetName())=="hpt_jk3" ||
	 string(obj->GetName())=="hpt_jk4" ||
	 string(obj->GetName())=="hpt_jk5" ||
	 string(obj->GetName())=="hpt_jk6" ||
	 string(obj->GetName())=="hpt_jk7" ||
	 string(obj->GetName())=="hpt_jk8" ||
	 string(obj->GetName())=="hpt_jk9" ||
	 string(obj->GetName())=="hpt_jk10"*/
	 )) {

      cout << "+" << flush;

      _jk = 0;
      if (TString(obj->GetName()).Contains("hpt_jk")) {
	assert( sscanf(obj->GetName(), "hpt_jk%d", &_jk) == 1);
      }
      _jet = TString(obj->GetName()).Contains("hpt_jet");

      TH1D *hpt = (TH1D*)obj;
      //TH1D *hpt2 = (TH1D*)indir2->Get("hnlo"); assert(hpt2);
      //  TH1D *hpt2 = (TH1D*)indir2->Get(jp::dagfile1 ? "mc/hpt_g" : "hgpt"); assert(hpt2);
          TH1D *hpt2 = (TH1D*)indir2->Get(jp::dagfile1 ? "mc/hpt_g" : "hpt"); assert(hpt2);

       if (hpt2)
        Unfold_histo(hpt, hpt2, outdir, ismc);

    } // hpt

    // Try to process friends similarly
    /*
    if (obj->InheritsFrom("TH1") &&
        (string(obj->GetName())=="hpt_ak5calo")) {

      cout << "-" << flush;

      _jk = 0; _jet = false;
      TH1D *hpt = (TH1D*)obj;
      TH1D *hpt2 = (TH1D*)indir2->Get("hnlo"); assert(hpt2);
      if (hpt2)
        dagostiniUnfold_histo(hpt, hpt2, outdir, ismc, false, "_ak5calo");
    } // hpt
    */
  } // while key

  curdir->cd();
} // recurseFile


void Unfold_histo(TH1D *hpt, TH1D *hnlo, TDirectory *outdir,
			   bool ismc, bool kscale, string id) {

  float y1, y2;
  assert(sscanf(outdir->GetName(),"Eta_%f-%f",&y1,&y2)==2);
  sscanf(outdir->GetName(),"Eta_%f-%f",&y1,&y2);
  cout << outdir->GetName() << " y1:" << y1 << " y2: " << y2 << endl;
  const char *c = id.c_str();
  if (_jk) c = Form("_jk%d",_jk);
  if (_jet) c = "_jet";


  _ismcjer = ismc;

  // initial fit of the NLO curve to a histogram
  TF1 *fnlo = new TF1(Form("fus%s",c),
		      //   "[0]*exp([1]/x)*pow(x,[2])"
		           "[0]*pow(x,[1])"
                      "*pow(1-x*cosh([3])/[4],[2])", //10., 1000.);
		      jp::unfptminnlo, min(jp::xmax, jp::emax/cosh(y1)));

  
  fnlo->SetParameters(5e10,-5.2,8.9,y1,jp::emax);
  //  fnlo->SetParameters(2e14*2e-10,-18,-5,10,y1,jp::emax);
  fnlo->FixParameter(3,y1);
  fnlo->FixParameter(4,jp::emax);

  //hnlo->Fit(fnlo,"QRN");
  //hnlo->Scale(2e-10); // TEMP PATCH
  fnlo->SetRange(max(60.,jp::unfptminnlo), min(jp::xmax, jp::emax/cosh(y1)));
  cout << "fit hnlo" << endl;
  hnlo->Fit(fnlo,"RN");  // There seems to be abnormal terminations
  fnlo->SetRange(jp::unfptminnlo, min(jp::xmax, jp::emax/cosh(y1)));

  // Graph of theory points with centered bins
  const double minerr = 0.02;
  TGraphErrors *gnlo = new TGraphErrors(0);
  TGraphErrors *gnlo2 = new TGraphErrors(0); // above + minerr
  gnlo->SetName("gnlo");
  gnlo2->SetName("gnlo2");
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {
  
    double y = hnlo->GetBinContent(i);
    double dy = hnlo->GetBinError(i);
    double ptmin = hnlo->GetBinLowEdge(i);
    double ptmax = hnlo->GetBinLowEdge(i+1);

    double y0 = fnlo->Integral(ptmin, ptmax) / (ptmax - ptmin);
    double x = fnlo->GetX(y0, ptmin, ptmax);

    int n = gnlo->GetN();
    tools::SetPoint(gnlo, n, x, y, 0, dy);
    tools::SetPoint(gnlo2, n, x, y, 0, tools::oplus(dy, minerr*y));
  }

  // Second fit to properly centered graph
  //gnlo2->Fit(fnlo,"QRN");
  fnlo->SetRange(max(60.,jp::unfptminnlo), min(jp::xmax, jp::emax/cosh(y1)));
  cout << "fit to gnlo2" << endl;
  gnlo2->Fit(fnlo,"RN");
  fnlo->SetRange(jp::unfptminnlo, min(jp::xmax, jp::emax/cosh(y1)));

  // Bin-centered data points
  TGraphErrors *gpt = new TGraphErrors(0);
  gpt->SetName(Form("gpt%s",c));
    for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

    double ptmin = hpt->GetBinLowEdge(i);
    double ptmax = hpt->GetBinLowEdge(i+1);
    double y = fnlo->Integral(ptmin, ptmax) / (ptmax - ptmin);
    double x = fnlo->GetX(y, ptmin, ptmax);
    double ym = hpt->GetBinContent(i);
    double ym_err = hpt->GetBinError(i);
    if (ym>0) {
      tools::SetPoint(gpt, gpt->GetN(), x, ym, 0., ym_err);
    }
  } // for i

  // Create smeared theory curve
  double maxpt = jp::emax/cosh(y1);
  cout << "y1 "<< y1 << " c "<< c <<endl<<flush;
  TF1 *fnlos = new TF1(Form("fs%s",c),smearedAnsatz,jp::unfptminnlo,maxpt,nk+3); 
  fnlos->SetParameters(y1, fnlo->GetParameter(0), fnlo->GetParameter(1),
                       fnlo->GetParameter(2), 0, 0);
  cout << "par0 "<< fnlos->GetParameter(0) << " y1 "<< y1<<endl<<flush; 

  // if (jp::debug)
    cout << "Calculate forward smearing and unfold hpt" << endl << flush;

  TGraphErrors *gfold_fwd = new TGraphErrors(0);
  gfold_fwd->SetName(Form("gfold_fwd%s",c));


  TGraphErrors *gcorrpt_fwd = new TGraphErrors(0);
  gcorrpt_fwd->SetName(Form("gcorrpt_fwd%s",c));
  TH1D *hcorrpt_fwd = (TH1D*)hpt->Clone(Form("hcorrpt_fwd%s",c));

 
  
  for (int i = 0; i != gpt->GetN(); ++i) {
    //  for (int i = 1; i != gpt->GetN()+1; ++i) {
    double x, y, ex, ey;

    tools::GetPoint(gpt, i, x, y, ex, ey);
    double k = fnlo->Eval(x)/fnlos->Eval(x); // fnlos causes a seg fault

    if (!TMath::IsNaN(k)) {

      tools::SetPoint(gfold_fwd, gfold_fwd->GetN(), x, k, ex, 0.);
      tools::SetPoint(gcorrpt_fwd, gcorrpt_fwd->GetN(), x, k*y, ex, k*ey);

      int j = hpt->FindBin(x);
      hcorrpt_fwd->SetBinContent(j, k*hpt->GetBinContent(j));
      hcorrpt_fwd->SetBinError(j, k*hpt->GetBinError(j));
    }
  }

  cout << "Forward smearing done!" << endl;

    
  outdir->cd();

  // Save resolution function
    TF1 *fres = new TF1(Form("fres%s",c), fPtRes, jp::xmin, jp::xmax, 1);
  fres->SetParameter(0, y1);

  // Store NLO ratio to (unsmeared) fit
  TGraphErrors *grationlo = new TGraphErrors(0);
  grationlo->SetName(Form("grationlo%s",c));
  for (int i = 0; i != gnlo2->GetN(); ++i) {

    double x, y, ex, ey;
    tools::GetPoint(gnlo2, i, x, y, ex, ey);
    double ys = fnlo->Eval(x);
    if (!TMath::IsNaN(ys))
      tools::SetPoint(grationlo, grationlo->GetN(), x, y / ys, ex, ey / ys);
  }

   
  // Store data ratio to (smeared) fit (this is also time consuming)
  TGraphErrors *gratio = new TGraphErrors(0);
  gratio->SetName(Form("gratio%s",c));
  for (int i = 0; i != gpt->GetN(); ++i) {

    double x, y, ex, ey;
    tools::GetPoint(gpt, i, x, y, ex, ey);
    double ys = fnlos->Eval(x);
    if (!TMath::IsNaN(ys))
      tools::SetPoint(gratio, gratio->GetN(), x, y / ys, ex, ey / ys);
      }

    outdir->cd();


 cout << "Smeared ratio stored" << endl;   
  if (!_jk && !_jet) {

 cout << "Saving plots" << endl;
    // Inputs and central method results
    hpt->Write("hpt");
    hnlo->Write("hnlo");
    gpt->Write("gpt");
    gnlo2->Write("gnlo");

     // Fit functions

    //cout << "Saving fit functions" << endl;
    //uResp->Write();
    fnlo->Write();
    // Calculating points for fs is taking significant time,
    // and even got stuck at some point when too far out of the range
    fnlos->SetRange(jp::unfptminnlo, min(jp::xmax, jp::emax/cosh(y1)));
     //     fnlos->SetNpx(1000); // otherwise ugly on log x-axis after write
     fnlos->Write();
     // Also calculating fnlos seems to take a lot of time when using effective JER... 
    
    fres->Write();
    grationlo->Write();
    gratio->Write();

    gfold_fwd->Write();
    gcorrpt_fwd->Write();
    hcorrpt_fwd->Write();

    cout << "Plots saved" << endl;
  }


} //Unfold_histo
