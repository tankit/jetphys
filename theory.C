// Purpose: Reformat theory curves
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2010
// Updated: June 2, 2015
//#include "ptresolution.h"
#include "tools.h"
#include "settings.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include <string>
#include <fstream>
#include <sstream>

using namespace std;

// Return largest of three numbers
inline double max3(double a, double b, double c) {
  return max(max(a,b),c);
}
// Return largest of three numbers
inline double min3(double a, double b, double c) {
  return min(min(a,b),c);
}

// Ansatz Kernel (for jer_systematics)
/*
int cnt_a = 0;
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) cout << "." << flush;

  const double pt = x[0]; // true pT, p[0] is measured pT
  const double eta = p[5];
  const double res = ptresolution(pt, eta+1e-3) * pt * (1. + p[6]);
  const double s = TMath::Gaus(p[0], pt, res, kTRUE);
  const double f = p[1] * exp(p[2]/pt) * pow(pt, p[3])
    * pow(1 - pt*cosh(eta)/3500., p[4]);

  return (f * s);
}

// Smeared Ansatz
double _ieps = 1e-12;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  if (!_kernel) _kernel = new TF1("_kernel", smearedAnsatzKernel, 1.,1000.,7);

  const double pt = x[0];
  const double eta = p[4];
  const double res = ptresolution(pt, eta+1e-3) * pt * (1 + p[5]);
  const double sigma = min(res, 0.30);
  const double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  const double ptmax = pt / (1. - 2.*sigma); // xmax*(1-2*sigma)=x
  const double par[7] = {pt, p[0], p[1], p[2], p[3], p[4], p[5]};
  _kernel->SetParameters(&par[0]);

  return ( _kernel->Integral(ptmin, ptmax, &par[0], _ieps) );
}
*/

//void theoryBin(TDirectory *din, TDirectory *dnlo, TDirectory *dnp,
//       TDirectory *dnlo1, TDirectory *dnlo2, TDirectory *dnlo3,
//       TDirectory *dnlo4,
//       TDirectory *das, TDirectory *dout);

// Load various theory predictions (FastNLO, LO MC)
void theoryBin(TDirectory *din, TDirectory *dth, TDirectory *dout);

// Also load old experimental data
void dataBin(TDirectory *din, TDirectory *dout);

void theory(string type) {


  //TFile *fin = new TFile(Form("output-%s-3a.root",type.c_str()),"READ");
  TFile *fin = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  // LO Pythia MC prediction
  TFile *fmc = new TFile("output-MC-2b.root","READ");
  assert(fmc && !fmc->IsZombie());

  // NB: central PDF is chosen later
  /*
  TFile *fnlo_cteq = (_algo=="AK7" ?
		      new TFile("fastnlo/fnl2332d_ct10-nlo_aspdf.root",
				"READ") :
		      new TFile("fastnlo/fnl2342b_ct10_aspdf.root","READ"));
  assert(fnlo_cteq && !fnlo_cteq->IsZombie());

  // CT10 from Gregory Soyez
  TFile *fnlo_ct10 = (_algo=="AK7" ?
		      new TFile("fastnlo/lhc7.root","READ") :
		      new TFile("fastnlo/lhc5.root","READ"));
  assert(fnlo_ct10 && !fnlo_ct10->IsZombie());

  TFile *fnlo_mstw = (_algo=="AK7" ?
		      new TFile("fastnlo/fnl2332d_mstw2008-nlo_aspdf.root",
				"READ") :
		      new TFile("fastnlo/fnl2342b_mstw2008nlo_aspdf.root",
				"READ"));
  assert(fnlo_mstw && !fnlo_mstw->IsZombie());

  TFile *fnlo_nnpdf = (_algo=="AK7" ?
		       new TFile("fastnlo/fnl2332d_nnpdf21-nlo_aspdf.root",
				 "READ") :
		       new TFile("fastnlo/fnl2342b_nnpdf21100_aspdf.root",
				 "READ"));
  assert(fnlo_nnpdf && !fnlo_nnpdf->IsZombie());

  // HERA10,15 empty for AK7?
  TFile *fnlo_hera = (_algo=="AK7" ?
		      new TFile("fastnlo/fnl2332d_hera10all_aspdf.root",
				"READ") :
		      new TFile("fastnlo/fnl2342b_hera10all_aspdf.root",
				"READ"));
  assert(fnlo_hera && !fnlo_hera->IsZombie());

  // Alpha-s variation
  TFile *fnlo_as = new TFile(_algo=="AK7" ?
			     "pdf4lhc/fnl2342a_ct10as_as-2.root" :
			     "pdf4lhc/fnl2342a_ct10as_as-2.root", "READ");
  assert(fnlo_as && !fnlo_as->IsZombie());

  TFile *fnp = new TFile(Form("fastnlo/NPcor%s.root",_algo.c_str()));
  assert(fnp && !fnp->IsZombie());
  */

  TFile *fout = new TFile(Form("output-%s-2c.root",type.c_str()),"RECREATE");
  assert(fout && !fout->IsZombie());

  // Select top category
  assert(fin->cd("Standard"));
  fin->cd("Standard");
  TDirectory *din0 = gDirectory;

  assert(fmc->cd("Standard"));
  fmc->cd("Standard");
  TDirectory *dmc0 = gDirectory;

  fout->mkdir("Standard");
  assert(fout->cd("Standard"));
  fout->cd("Standard");
  TDirectory *dout0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = din0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")
	&& string(obj->GetName())!="Eta_0.0-1.3"
	&& string(obj->GetName())!="Eta_3.0-3.2"
	&& string(obj->GetName())!="Eta_3.2-4.7") {

      assert(din0->cd(obj->GetName()));
      din0->cd(obj->GetName());
      TDirectory *din = gDirectory;

      assert(dmc0->cd(obj->GetName()));
      dmc0->cd(obj->GetName());
      TDirectory *dmc = gDirectory;
      /*
      assert(fnlo_cteq->cd());
      TDirectory *dnlo_cteq = gDirectory;
      assert(fnlo_ct10->cd());
      TDirectory *dnlo_ct10 = gDirectory;
      assert(fnlo_mstw->cd());
      TDirectory *dnlo_mstw = gDirectory;
      assert(fnlo_nnpdf->cd());
      TDirectory *dnlo_nnpdf = gDirectory;
      assert(fnlo_hera->cd());
      TDirectory *dnlo_hera = gDirectory;
      //
      assert(fnlo_as->cd());
      TDirectory *das = gDirectory;

      assert(fnp->cd());
      TDirectory *dnp = gDirectory;
      */

      dout0->mkdir(obj->GetName());
      assert(dout0->cd(obj->GetName()));
      dout0->cd(obj->GetName());
      TDirectory *dout = gDirectory;

      // Process subdirectory
      //theoryBin(din, dnlo_cteq, dnp, dnlo_ct10, dnlo_mstw, dnlo_nnpdf,
      //	dnlo_hera,
      //	das, dout);
      theoryBin(din, dmc, dout);
      //dataBin(din, dout);
    } // inherits TDirectory
  } // while

  cout << "Output stored in " << fout->GetName() << endl;
  fout->Write();
  fout->Close();
  fout->Delete();

  fin->Close();
  fin->Delete();

  fmc->Close();
  /*
  fnlo_cteq->Close();
  fnlo_cteq->Delete();
  fnlo_ct10->Close();
  fnlo_ct10->Delete();
  fnlo_mstw->Close();
  fnlo_mstw->Delete();
  fnlo_nnpdf->Close();
  fnlo_nnpdf->Delete();
  fnlo_as->Close();
  fnlo_as->Delete();
  */

  // Empty work space so that next steps won't get messed up
  //if (gROOT->IsBatch()) {
  //gROOT->Clear();
  //}
} // theory


//void theoryBin(TDirectory *din, TDirectory *dnlo0, TDirectory *dnp,
//       TDirectory *dnlo1, TDirectory *dnlo2, TDirectory *dnlo3,
//       TDirectory *dnlo4,
//       TDirectory *das, TDirectory *dout){//, bool ak7) {

void theoryBin(TDirectory *din, TDirectory *dth, TDirectory *dout) {

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);
  sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax);
  int ieta = int((0.5*(etamin+etamax))/0.5);
  int jeta = (string(jp::algo)=="AK7" ? min(4,ieta) : ieta);

  // inclusive jets
  TH1D *hpt = (TH1D*)din->Get("hpt"); assert(hpt);

  // theory curves
  // http://www-ekp.physik.uni-karlsruhe.de/~rabbertz/fastNLO_LHC/InclusiveJets/
  // -> fnl2342_cteq66_aspdf_full.root
  // Numbering scheme explained in
  // https://twiki.cern.ch/twiki/bin/view/CMS/CMSfastNLO
  // 2-point to 6-point theory uncertainty:
  // h200X00->h300X09, h400X00->h300X08
  TH1D *hnlo(0);//, *hpdfup(0), *hpdfdw(0), *hscup(0), *hscdw(0);//, *herr(0);

  cout << din->GetName() << " ieta="<<ieta<< " jeta="<<jeta << endl;

  TH1D *hmc = (TH1D*)dth->Get("hpt_g0tw");
  assert(hmc);

  /*
  int S = (string(jp::algo)=="AK7" ? 1 : 3);
  TH1D *hnlo_cteq = (TH1D*)dnlo0->Get(Form("h%d00%d00",S,jeta+1));
  assert(hnlo_cteq);

  hnlo_cteq->Scale(1000.);

  // CTEQ10
  TH1D *hnlo_ct10 = (TH1D*)dnlo1->Get(Form("h300%d00",jeta+1));
  TH1D *hpdfup_ct10 = (TH1D*)dnlo1->Get(Form("h300%d02",jeta+1));
  TH1D *hpdfdw_ct10 = (TH1D*)dnlo1->Get(Form("h300%d01",jeta+1));
  TH1D *hscup_ct10 = (TH1D*)dnlo1->Get(Form("h300%d09",jeta+1));
  TH1D *hscdw_ct10 = (TH1D*)dnlo1->Get(Form("h300%d08",jeta+1));
  assert(hnlo_ct10); assert(hpdfup_ct10); assert(hpdfdw_ct10);
  assert(hscup_ct10); assert(hscdw_ct10);

  hnlo_ct10->Scale(1000.);
  hpdfup_ct10->Scale(1./1.65);
  hpdfdw_ct10->Scale(1./1.65);

  // MSTW2008
  TH1D *hnlo_mstw = (TH1D*)dnlo2->Get(Form("h%d00%d00",S,jeta+1));
  assert(hnlo_mstw);

  hnlo_mstw->Scale(1000.);

  // NNPDF2010
  TH1D *hnlo_nnpdf = (TH1D*)dnlo3->Get(Form("h%d00%d00",S,jeta+1));
  assert(hnlo_nnpdf);

  hnlo_nnpdf->Scale(1000.);

  // HERA10
  TH1D *hnlo_hera = (TH1D*)dnlo4->Get(Form("h%d00%d00",S,jeta+1));
  TH1D *hpdfup_hera = (TH1D*)dnlo4->Get(Form("h%d00%d02",S,jeta+1));
  TH1D *hpdfdw_hera = (TH1D*)dnlo4->Get(Form("h%d00%d01",S,jeta+1));
  assert(hnlo_hera); assert(hpdfup_hera); assert(hpdfdw_hera);

  hnlo_hera->Scale(1000.);

  // alpha-S
  TH1D *hasup = (TH1D*)das->Get(Form("h300%d02",ieta+1));
  TH1D *hasdw = (TH1D*)das->Get(Form("h300%d01",ieta+1));
  assert(hasup); assert(hasdw);

  TH1D *hnp = (TH1D*)dnp->Get(Form("corr%d",min(int(etamin/0.5+0.5),5)));

  if (!hnp) cout << "eta: " << etamin << " " << etamax << endl;
  assert(hnp);
  */


  // make sure new histograms get created in the output file
  dout->cd();

  /*
  TH1D *hnpup(0), *hnpdw(0), *hsysup(0), *hsysdw(0);
  {
    // Move CTEQ6.6
    hnlo_cteq = (TH1D*)hnlo_cteq->Clone("hnlo0_ct10k");
    hnlo_ct10 = (TH1D*)hnlo_ct10->Clone("hnlo0_ct10");
    hnlo_mstw = (TH1D*)hnlo_mstw->Clone("hnlo0_mstw");
    hnlo_nnpdf = (TH1D*)hnlo_nnpdf->Clone("hnlo0_nnpdf");
    hnlo_hera = (TH1D*)hnlo_hera->Clone("hnlo0_hera");

    // Determine PDF4LHC from CT10, MSTW08, NNPDF10
    hnlo = (TH1D*)hnlo_ct10->Clone("hnlo0"); // central PDF again later
    hpdfup = (TH1D*)hpdfup_ct10->Clone("hpdfup"); // central PDF
    hpdfdw = (TH1D*)hpdfdw_ct10->Clone("hpdfdw"); // central PDF
    for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

      // Sanity checks
      assert(hpdfup_ct10->GetBinContent(i)>=0);
      assert(hpdfdw_ct10->GetBinContent(i)<=0);
    }
    hscup = (TH1D*)hscup_ct10->Clone("hscup"); // central PDF
    hscdw = (TH1D*)hscdw_ct10->Clone("hscdw"); // central PDF
    hasup = (TH1D*)hasup->Clone("hasup");
    hasdw = (TH1D*)hasdw->Clone("hasdw");

    hnp = (TH1D*)hnp->Clone("hnpcorr");

    hnpup = (TH1D*)hnlo->Clone("hnpup");
    hnpdw = (TH1D*)hnlo->Clone("hnpdw");

    hsysup = (TH1D*)hnlo->Clone("hsysup");
    hsysdw = (TH1D*)hnlo->Clone("hsysdw");
  }
  */

  /*
  // Patch up the hnp to adapt the wrong binning at low pT
  if (true) { // patch hnp

    TH1D *hnp_tmp = (TH1D*)hnlo->Clone("hnp_tmp");
    for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

      double x = hnlo->GetBinCenter(i);
      int j1 = hnp->FindBin(x);
      double x1 = hnp->GetBinCenter(j1);
      double y1 = hnp->GetBinContent(j1);
      double ey1 = hnp->GetBinError(j1);
      int j2 = (x>x1? j1+1 : j1-1);
      double x2 = hnp->GetBinCenter(j2);
      double y2 = hnp->GetBinContent(j2);
      double ey2 = hnp->GetBinError(j2);
      double y = y1 + (y2-y1) / (x2-x1) * (x-x1);
      double ey = ey1 + (ey2-ey1) / (x2-x1) * (x-x1);
      hnp_tmp->SetBinContent(i, y);
      hnp_tmp->SetBinError(i, ey);
    } // for i
    delete hnp;
    hnp = hnp_tmp;
    hnp->SetName("hnpcorr");
  } // patch hnp
  */

  /*
  // Turn NP uncertainty into uncertainty histogram
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

    int j = hnp->FindBin(hnlo->GetBinCenter(i));
    hnpup->SetBinContent(i, +hnp->GetBinError(j)/hnp->GetBinContent(j));
    hnpup->SetBinError(i, 0.);
    hnpdw->SetBinContent(i, -hnp->GetBinError(j)/hnp->GetBinContent(j));
    hnpdw->SetBinError(i, 0.);
    if (!(hnp->GetBinLowEdge(j)>=hnlo->GetBinLowEdge(i)) ||
	!(hnp->GetBinLowEdge(j+1)<=hnlo->GetBinLowEdge(i+1))) {
      cout << Form("hnlo = [%1.0f,%1.0f]; hnp = [%1.0f,%1.0f]",
		   hnlo->GetBinLowEdge(i), hnlo->GetBinLowEdge(i+1),
		   hnp->GetBinLowEdge(j), hnp->GetBinLowEdge(j+1)) << endl;
      if (hnlo->GetBinLowEdge(i)>100 && hnlo->GetBinLowEdge(i)<2500.) {
	assert(hnp->GetBinLowEdge(j)>=hnlo->GetBinLowEdge(i));
	assert(hnp->GetBinLowEdge(j+1)<=hnlo->GetBinLowEdge(i+1));
      }
    }
  } // for i
  */

  /*
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

    double x = hnlo->GetBinCenter(i);
    int i1 = hnpup->FindBin(x);
    double y1a = hnpup->GetBinContent(i1);
    double y1b = hnpdw->GetBinContent(i1);
    int i2 = hscup->FindBin(x);
    double y2a = hscup->GetBinContent(i2);
    double y2b = hscdw->GetBinContent(i2);
    int i3 = hpdfup->FindBin(x);
    double y3a = hpdfup->GetBinContent(i3);
    double y3b = hpdfdw->GetBinContent(i3);
    int i4 = hasup->FindBin(x);
    double y4a = hasup->GetBinContent(i4);
    double y4b = hasdw->GetBinContent(i4);
    double errup = sqrt(pow(max3(y1a,y1b,0),2) + pow(max3(y2a,y2b,0),2)
			+ pow(max3(y3a,y3b,0),2) + pow(max3(y4a,y4b,0),2));
    double errdw = sqrt(pow(min3(y1a,y1b,0),2) + pow(min3(y2a,y2b,0),2)
			+ pow(min3(y3a,y3b,0),2) + pow(min3(y4a,y4b,0),2));
    hsysup->SetBinContent(i, errup);
    hsysup->SetBinError(i, 0.);
    hpdfup->SetBinError(i, 0.);
    hscup->SetBinError(i, 0.);
    hasup->SetBinError(i, 0);
    hsysdw->SetBinContent(i, -errdw);
    hsysdw->SetBinError(i, 0.);
    hpdfdw->SetBinError(i, 0.);
    hscdw->SetBinError(i, 0.);
    hasdw->SetBinError(i, 0);
  } // for i
  */

  /*
  for (int i = 1; i != hnlo_ct10->GetNbinsX()+1; ++i) {

    double x = hnlo_ct10->GetBinCenter(i);
    int jnp = hnp->FindBin(x);
    double np = hnp->GetBinContent(jnp);
    int j = hnlo_ct10->FindBin(x);
    hnlo->SetBinContent(j, 1e-3*hnlo->GetBinContent(j)*np);
    int jc = hnlo_ct10->FindBin(x);
    hnlo_ct10->SetBinContent(jc, 1e-3*hnlo_ct10->GetBinContent(jc)*np);
    int jk = hnlo_cteq->FindBin(x);
    hnlo_cteq->SetBinContent(jk, 1e-3*hnlo_cteq->GetBinContent(jk)*np);
    int jm = hnlo_mstw->FindBin(x);
    hnlo_mstw->SetBinContent(jm, 1e-3*hnlo_mstw->GetBinContent(jm)*np);
    int jn = hnlo_nnpdf->FindBin(x);
    hnlo_nnpdf->SetBinContent(jn, 1e-3*hnlo_nnpdf->GetBinContent(jn)*np);
    int jh = hnlo_hera->FindBin(x);
    hnlo_hera->SetBinContent(jh, 1e-3*hnlo_hera->GetBinContent(jh)*np);
    //
    hnlo->SetBinError(j, 0.);
    hnlo_ct10->SetBinError(jc, 0.);
    hnlo_cteq->SetBinError(jk, 0.);
    hnlo_mstw->SetBinError(jm, 0.);
    hnlo_nnpdf->SetBinError(jn, 0.);
    hnlo_hera->SetBinError(jh, 0.);
  }
  */

  // initial fit of the NLO curve
  TF1 *fnlo = new TF1("fnlo",
		      "[0]*exp([1]/x)*pow(x,[2])"
		      "*pow(1-x*cosh([4])/3500.,[3])", 10., 1000.);
  fnlo->SetParameters(2e14,-18,-5.2,8.9,0.5*ieta);
  fnlo->FixParameter(4,0.5*ieta);

  hnlo = hmc;
  hnlo->Fit(fnlo,"QRN");

  // Graph of theory points with centered bins
  const double minerr = 0.02;
  TGraphErrors *gnlo = new TGraphErrors(0);
  TGraphErrors *gnlo2 = new TGraphErrors(0); // above + minerr
  TGraphErrors *gnlocut = new TGraphErrors(0); // only up to expected pT
  gnlo->SetName("gnlo");
  gnlo2->SetName("gnlo2");
  gnlocut->SetName("gnlocut");
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
    if (y*(ptmax-ptmin) > 0.0001) { // for 10/fb
      int m = gnlocut->GetN();
      tools::SetPoint(gnlocut, m, x, y, 0, 0);
    }
  }

  gnlo2->Fit(fnlo,"QRN");

  // Divide graph with fit to check stability
  TGraphErrors *gnlofit = new TGraphErrors(0);
  gnlofit->SetName("gnlofit");
  for (int i = 0; i != gnlo->GetN(); ++i) {

    double x, y, ex, ey;
    tools::GetPoint(gnlo, i, x, y, ex, ey);
    double f = fnlo->Eval(x);
    tools::SetPoint(gnlofit, i, x, y/f, ex, ey/f);
  }

  // Rebin theory to match data bins
  TH1D *hnlo0 = hnlo;
  hnlo0->SetName("hnlo0");
  hnlo = tools::Rebin(hnlo0, hpt);
  hnlo->SetName("hnlo");

  // Same for earlier 2010 and 2011 studies to get a matching theory set
  // (moved from PDF4LHC to CT10 since 2010)
  if (false) {

    // Match this to drawSummary.C inputs
    TDirectory *curdir = gDirectory;
    TFile *fin2010 = new TFile("../pfjet/output-DATA-4c.root","READ");
    assert(fin2010 && !fin2010->IsZombie());
    assert(fin2010->cd("Standard"));
    fin2010->cd("Standard");
    assert(gDirectory->cd(din->GetName()));
    gDirectory->cd(din->GetName());
    TDirectory *din2010 = gDirectory;
    TH1D *hpt2010 = (TH1D*)din2010->Get("hpt");
    curdir->cd();
    TH1D *hnlo2010 = tools::Rebin(hnlo0, hpt2010);
    hnlo2010->SetName("hnlo2010");
    hnlo2010->Write();
  }
  if (false) {
    // Match this to drawSummary.C inputs
    TDirectory *curdir = gDirectory;
    TFile *fin2011 = new TFile("backup/oct26/output-DATA-4c.root","READ");
    assert(fin2011 && !fin2011->IsZombie());
    assert(fin2011->cd("Standard"));
    fin2011->cd("Standard");
    assert(gDirectory->cd(din->GetName()));
    gDirectory->cd(din->GetName());
    TDirectory *din2011 = gDirectory;
    TH1D *hpt2011 = (TH1D*)din2011->Get("hpt");
    curdir->cd();
    TH1D *hnlo2011 = tools::Rebin(hnlo0, hpt2011);
    hnlo2011->SetName("hnlo2011");
    hnlo2011->Write();
  }

  /*
  {
    TH1D *hnlo0_cteq = hnlo_cteq;
    hnlo0_cteq->SetName("hnlo0_ct10k");
    hnlo_cteq = tools::Rebin(hnlo0_cteq, hpt);
    hnlo_cteq->SetName("hnlo_ct10k");
    //
    TH1D *hnlo0_ct10 = hnlo_ct10;
    hnlo0_ct10->SetName("hnlo0_ct10");
    hnlo_ct10 = tools::Rebin(hnlo0_ct10, hpt);
    hnlo_ct10->SetName("hnlo_ct10");
    //
    TH1D *hnlo0_mstw = hnlo_mstw;
    hnlo0_mstw->SetName("hnlo0_mstw");
    hnlo_mstw = tools::Rebin(hnlo0_mstw, hpt);
    hnlo_mstw->SetName("hnlo_mstw");
    //
    TH1D *hnlo0_nnpdf = hnlo_nnpdf;
    hnlo0_nnpdf->SetName("hnlo0_nnpdf");
    hnlo_nnpdf = tools::Rebin(hnlo0_nnpdf, hpt);
    hnlo_nnpdf->SetName("hnlo_nnpdf");
    //
    TH1D *hnlo0_hera = hnlo_hera;
    hnlo0_hera->SetName("hnlo0_hera");
    hnlo_hera = tools::Rebin(hnlo0_hera, hpt);
    hnlo_hera->SetName("hnlo_hera");
  }
  */

  dout->cd();
  gnlo->Write();
  gnlocut->Write();
  gnlofit->Write();
  fnlo->Write();

  din->cd();
} // theoryBin


void dataBin(TDirectory *din, TDirectory *dout) {

  float y1, y2;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&y1,&y2)==2);
  sscanf(din->GetName(),"Eta_%f-%f",&y1,&y2);
  int iy = int((0.5*(y1+y2))/0.5);

  TH1D *hpt = (TH1D*)din->Get("hpt"); assert(hpt);

  // Read in text file
  if (string(jp::algo)=="AK7") {
    const int nchr = 2048;
    char chr[nchr];
    ifstream fs("fastnlo/InclusiveJets_Table_postCWR.txt", ios::in);
    assert(fs.is_open());
    fs.getline(chr, nchr);

    dout->cd();
    TH1D *h = (TH1D*)hpt->Clone("hdata2011");
    h->Reset();
    TH1D *hdw = (TH1D*)hpt->Clone("hdata2011_dw");
    hdw->Reset();
    TH1D *hup = (TH1D*)hpt->Clone("hdata2011_up");
    hup->Reset();
    while(fs.getline(chr, nchr)) {

      //cout << chr << endl;

      float ymin, ymax, xsec, err, err2, uncorr, np;
      int ptmin, ptmax;
      assert(sscanf(chr, "%f %f %d %d %f %f %f %f %f",
		    &ymin, &ymax, &ptmin, &ptmax,
		    &xsec, &err, &uncorr, &err2, &np)==9);
      sscanf(chr, "%f %f %d %d %f %f %f %f %f",
	     &ymin, &ymax, &ptmin, &ptmax,
	     &xsec, &err, &uncorr, &err2, &np);

      // sum up the uncertainty
      float foo;
      stringstream str(chr);
      for (int i = 0; i != 9; ++i) str >> foo; assert(foo==np);
      //assert(str >> foo);
      str >> foo;
      //assert(str >> foo); // lumierr
      str >> foo; // lumierr
      float eup, edw, sumedw(foo*foo), sumeup(foo*foo);
      while (str >> edw >> eup) {
	sumeup += pow(max(eup, -edw),2);
	sumedw += pow(min(eup, -edw),2);
      }
      eup = sqrt(sumeup);
      edw = sqrt(sumedw);

      if (fabs(ymin-y1)<0.1 && fabs(ymax-y2)<0.1) {
	int i = h->FindBin(0.5*(ptmin+ptmax));
	if (h->GetBinLowEdge(i)==ptmin &&
	    h->GetBinLowEdge(i+1)==ptmax) {
	  h->SetBinContent(i, xsec);
	  h->SetBinError(i, xsec*err);
	  hup->SetBinContent(i, xsec*(1+eup));
	  hup->SetBinError(i, xsec*(1+eup)*err);
	  hdw->SetBinContent(i, xsec*(1-edw));
	  hdw->SetBinError(i, xsec*(1-edw)*err);
	}
	else {
	  cout << "ymin = " << y1 << " ymax = " << y2
	       << "ptmin = " << ptmin << " ptmax = " << ptmax
	       << " pt1 = " << h->GetBinLowEdge(i)
	       << " pt2 = " << h->GetBinLowEdge(i+1) << endl;
	}
      }

    } // while
  } // ak7

  if (string(jp::algo)=="AK5") {

    const int nchr = 2048;
    char chr[nchr];
    ifstream fs(Form("fastnlo/InclusiveJets2010_Table%d.txt",iy), ios::in);
    assert(fs.is_open());
    // remove header lines
    for (int i = 0; i != 10; ++i) fs.getline(chr, nchr);

    dout->cd();
    TH1D *h = (TH1D*)hpt->Clone("hdata2010");
    h->Reset();
    TH1D *hup = (TH1D*)hpt->Clone("hdata2010_up");
    hup->Reset();
    TH1D *hdw = (TH1D*)hpt->Clone("hdata2010_dw");
    hdw->Reset();
    while(fs.getline(chr, nchr)) {

      float xsec, err, err2, sys, sys2;
      float pt, ptmin, ptmax;
      assert(string(chr)=="" ||
	     sscanf(chr, "%f %f %f %f %f %f %f %f",
		    &pt, &ptmin, &ptmax,
		    &xsec, &err, &err2, &sys, &sys2)==8);
      if(string(chr)!="")
	sscanf(chr, "%f %f %f %f %f %f %f %f",
	       &pt, &ptmin, &ptmax,
	       &xsec, &err, &err2, &sys, &sys2);

      int i = h->FindBin(0.5*(ptmin+ptmax));
      if (h->GetBinLowEdge(i)==ptmin &&
	  h->GetBinLowEdge(i+1)==ptmax) {
	h->SetBinContent(i, xsec);
	h->SetBinError(i, err);
	hup->SetBinContent(i, xsec+sys);
	hup->SetBinError(i, err);
	hdw->SetBinContent(i, xsec+sys2);
	hdw->SetBinError(i, err);
      }
      else {
	cout << "ymin = " << y1 << " ymax = " << y2
	     << "ptmin = " << ptmin << " ptmax = " << ptmax
	     << " pt1 = " << h->GetBinLowEdge(i)
	     << " pt2 = " << h->GetBinLowEdge(i+1) << endl;
      }
    } // while
  } // ak5


  // Read in text file (SMP-13-002)
  if (true) {
    const int nchr = 2048;
    char chr[nchr];
    ifstream fs(Form("fastnlo/InclusiveJets_Table_%s.txt",jp::algo));
    assert(fs.is_open());
    fs.getline(chr, nchr);

    dout->cd();
    TH1D *h = (TH1D*)hpt->Clone("hdata2013");
    h->Reset();
    TH1D *hdw = (TH1D*)hpt->Clone("hdata2013_dw");
    hdw->Reset();
    TH1D *hup = (TH1D*)hpt->Clone("hdata2013_up");
    hup->Reset();
    while(fs.getline(chr, nchr)) {

      float ymin, ymax, xsec, err, err2, uncorr, np;
      int ptmin, ptmax;
      assert(sscanf(chr, "%f %f %d %d %f %f %f %f %f",
		    &ymin, &ymax, &ptmin, &ptmax,
		    &xsec, &err, &uncorr, &err2, &np)==9);
      sscanf(chr, "%f %f %d %d %f %f %f %f %f",
	     &ymin, &ymax, &ptmin, &ptmax,
	     &xsec, &err, &uncorr, &err2, &np);

      // sum up the uncertainty
      float foo;
      stringstream str(chr);
      for (int i = 0; i != 9; ++i) str >> foo; assert(foo==np);
      //assert(str >> foo);
      str >> foo;
      //assert(str >> foo); // lumierr
      str >> foo; // lumierr
      float eup, edw, sumedw(foo*foo), sumeup(foo*foo);
      while (str >> edw >> eup) {
	sumeup += pow(max(eup, -edw),2);
	sumedw += pow(min(eup, -edw),2);
      }
      eup = sqrt(sumeup);
      edw = sqrt(sumedw);

      if (fabs(ymin-y1)<0.1 && fabs(ymax-y2)<0.1) {
	int i = h->FindBin(0.5*(ptmin+ptmax));
	if (h->GetBinLowEdge(i)==ptmin &&
	    h->GetBinLowEdge(i+1)==ptmax) {
	  h->SetBinContent(i, xsec);
	  h->SetBinError(i, xsec*err);
	  hup->SetBinContent(i, xsec*(1+eup));
	  hup->SetBinError(i, xsec*(1+eup)*err);
	  hdw->SetBinContent(i, xsec*(1-edw));
	  hdw->SetBinError(i, xsec*(1-edw)*err);
	}
	else {
	  cout << "ymin = " << y1 << " ymax = " << y2
	       << "ptmin = " << ptmin << " ptmax = " << ptmax
	       << " pt1 = " << h->GetBinLowEdge(i)
	       << " pt2 = " << h->GetBinLowEdge(i+1) << endl;
	}
      }

    } // while
  } // ak7

}
