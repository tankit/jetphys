// Purpose: Estimation of inclusive jet cross section systematics
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: June 11, 2015
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include "ptresolution.h"
#include "tools.h"
#include "settings.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Ansatz Kernel (for jer_systematics
int cnt_a = 0;
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) cout << "." << flush;

  const double pt = x[0]; // true pT, p[0] is measured pT
  const double eta = p[5];
  const double res = ptresolution(pt, eta+1e-3) * pt * (1. + p[6]);
  const double s = TMath::Gaus(p[0], pt, res, kTRUE);
  const double f = p[1] * exp(p[2]/pt) * pow(pt, p[3])
    * pow(1 - pt*cosh(eta)/jp::emax., p[4]);

  return (f * s);
}

// Smeared Ansatz
double _epsilon = 1e-12;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  if (!_kernel) _kernel = new TF1("_kernel", smearedAnsatzKernel, 1.,1000.,7);

  const double pt = x[0];
  const double eta = p[4];
  const double res = ptresolution(pt, eta+1e-3) * pt * (1 + p[5]);
  const double sigma = min(res, 0.30);
  double xmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  xmin = max(1.,xmin); // safety check
  double xmax = pt / (1. - 2.*sigma); // xmax*(1-2*sigma)=x
  xmax = min(jp::emax/cosh(eta),xmax); // safety check
  const double par[7] = {pt, p[0], p[1], p[2], p[3], p[4], p[5]};
  _kernel->SetParameters(&par[0]);

  return ( _kernel->Integral(xmin, xmax, &par[0], _epsilon) );
}

// JEC systematics are evaluated in two ways:
// - vary JEC by uncertainty and rescale original data => ratio of histos
// - shift Ansatzt by JEC uncertainty => ratio of functions
// The former usually suffers from low statistics fluctuations at high pT,
// but still works as a good cross check for the latter method. The Ansatz
// approach also decouples unfolding from JEC uncertainty, which could be
// good or bad (not quite sure yet how the correlation should be treated).
// The Ansatz method is very fast and because it does not suffer from low
// statistics, can be used to break JEC uncertainty into uncorrelated sources
// (D0 measurement of inclusive jets used 91 sources for JEC uncertainty...)

// systematics container for calculating total
struct sysc {

  TH1D *plus;
  TH1D *minus;
  TH1D *av;

  sysc(TH1D *pl = 0, TH1D *mn = 0, TH1D *avg = 0) :
    plus(pl), minus(mn), av(avg) { };
};

sysc *jec_systematics(TDirectory *dzr, TDirectory *dunc,
		      TDirectory *dpl, TDirectory *dmn,
		      TDirectory *dout, string type,
		      string jectype = "abs");

void jec_shifts(TDirectory *dzr, TDirectory *dout, string type,
		string algo = "jpt");

sysc *jer_systematics(TDirectory *dzr, TDirectory *dout,
		      string type, string jectype="inc");

sysc *lum_systematics(TDirectory *dzr, TDirectory *dout);

void tot_systematics(TDirectory *dzr, TDirectory *dout,
		     sysc *cjec1, sysc *cjec2, sysc *cjec3,
		     sysc *cjer1, sysc *cjer2,
		     sysc *clum);

void statistics(TDirectory *dzr, TDirectory *dout);
void sourceBin(TDirectory *dth, TDirectory *dout);

bool _ismc = false; // implement as non-global later
void systematics(string type) {

  _ismc = (type=="MC" || type=="HW");

  // Could do ratio of JEC histos before or after unfolding, i.e.
  // output2.root or output3.root. The answer could be slightly different
  // due to correlations between unfolding and JEC

  // Unfolded data
  TFile *fin = new TFile(Form("output-%s-3.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  // Raw data
  TFile *fin2 = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin2 && !fin2->IsZombie());

  TFile *fout = new TFile(Form("output-%s-4.root",type.c_str()),
			  "RECREATE");
  assert(fout && !fout->IsZombie());

  // Select top categories for JEC uncertainty
  if (!fin->cd("JECPlus"))
    assert(fin->cd("Standard"));
  TDirectory *dpl0 = gDirectory;
  if (!fin->cd("JECMinus"))
    assert(fin->cd("Standard"));
  TDirectory *dmn0 = gDirectory;
  if (!fin->cd("ResJEC"))
    assert(fin->cd("Standard"));
  TDirectory *dzr0 = gDirectory;
  if (!fin2->cd("ResJEC"))
    assert(fin2->cd("Standard"));
  TDirectory *dunc0 = gDirectory;


  fout->mkdir("Standard");
  assert(fout->cd("Standard"));
  TDirectory *dout0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = dzr0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")) {

      assert(dzr0->cd(obj->GetName()));
      TDirectory *dzr = gDirectory;
      assert(dpl0->cd(obj->GetName()));
      TDirectory *dpl = gDirectory;
      assert(dmn0->cd(obj->GetName()));
      TDirectory *dmn = gDirectory;
      assert(dunc0->cd(obj->GetName()));
      TDirectory *dunc = gDirectory;

      dout0->mkdir(obj->GetName());
      assert(dout0->cd(obj->GetName()));
      TDirectory *dout = gDirectory;

      // Process subdirectory
      sysc *cjec  = jec_systematics(dzr,dunc,dpl,dmn, dout, type, "tot");
      if(cjec);
      sysc *cjec1 = jec_systematics(dzr,dunc,dpl,dmn, dout, type, "abs");
      sysc *cjec2 = jec_systematics(dzr,dunc,dpl,dmn, dout, type, "rel");
      sysc *cjec3 = jec_systematics(dzr,dunc,dpl,dmn, dout, type, "bjt");
      sysc *cjec4 = jec_systematics(dzr,dunc,dpl,dmn, dout, type, "pt");
      sysc *cjer1 = jer_systematics(dzr, dout, type, "inc");
      sysc *cjer2 = jer_systematics(dzr, dout, type, "bjt");
      sysc *clum = lum_systematics(dzr, dout);
      tot_systematics(dzr, dout, cjec1, cjec2, cjec3, cjer1, cjer2, clum);
      statistics(dzr, dout);

      jec_shifts(dzr, dout, type, "pf");
    } // inherits TDirectory
  } // while

  cout << "Output stored in " << fout->GetName() << endl;
  fout->Write();
  fout->Close();
  fout->Delete();

  fin->Close();
  fin->Delete();
} // systematics


// JEC systematics
sysc *jec_systematics(TDirectory *dzr, TDirectory *dunc,
		      TDirectory *dpl, TDirectory *dmn,
		      TDirectory *dout, string type,
		      string jectype) {

  float etamin, etamax;
  assert(sscanf(dzr->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);

  // Load the uncertainty
  TH1D *hunc = (TH1D*)dunc->Get("punc"); assert(hunc);
  //JetCorrectionUncertainty *func = new JetCorrectionUncertainty(Form("CondFormats/JetMETObjects/data/GR_R_42_V23_Uncertainty_%sPF.txt",jp::algo));
  string s = Form("CondFormats/JetMETObjects/data/%s_%s_Uncertainty_%sPF.txt", jp::jecgt, jp::type, jp::algo);
  cout << s << endl << flush;
  JetCorrectionUncertainty *func = new JetCorrectionUncertainty(s.c_str());
  const char *jt = jectype.c_str();

  // inclusive jets
  TH1D *hzr = (TH1D*)dzr->Get("hpt"); assert(hzr);
  TH1D *hpl = (TH1D*)dpl->Get("hpt"); assert(hpl);
  TH1D *hmn = (TH1D*)dmn->Get("hpt"); assert(hmn);

  // make sure new histograms get created in the output file
  dout->cd();

  // inclusive jets
  TF1 *fpt0 = (TF1*)dzr->Get("fus"); assert(fpt0); fpt0->SetName("fpt0");
  TF1 *fpt = new TF1("fpt1","[0]*exp([1]/(x*(1-[6])))*pow(x*(1-[6]),[2])"
		     "*pow(1-x*(1-[6])*cosh([4])/[5], [3])",
		     jp::recopt, jp::emax/cosh(etamin));
  fpt->SetParameters(fpt0->GetParameter(0), fpt0->GetParameter(1),
		     fpt0->GetParameter(2), fpt0->GetParameter(3),
		     fpt0->GetParameter(4), fpt0->GetParameter(5),
		     0.);

  // first estimate is just ratio of histograms with shifted JEC
  // inclusive jets
  TH1D *hjpl0 = (TH1D*)hpl->Clone(Form("hjec_%s_pl0", jt));
  hjpl0->Divide(hzr);
  TH1D *hjmn0 = (TH1D*)hmn->Clone(Form("hjec_%s_mn0", jt));
  hjmn0->Divide(hzr);
  TH1D *hjav0 = (TH1D*)hzr->Clone(Form("hjec_%s_av0", jt));
  hjav0->Reset();
  // Center around zero, calculate average uncertainty
  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {

    if (hzr->GetBinContent(i)!=0) {

      double pl = hjpl0->GetBinContent(i)-1;
      double epl = hjpl0->GetBinError(i);
      hjpl0->SetBinContent(i, pl);
      double mn = hjmn0->GetBinContent(i)-1;
      double emn = hjmn0->GetBinError(i);
      hjmn0->SetBinContent(i, mn);

      // Average error should be weighted by statistics
      double av = (epl&&emn ? (pl+mn*epl*epl/(emn*emn)) / (1+(epl*epl)/(emn*emn))
		   : 0.);
      double eav = (epl&&emn ? epl / sqrt(1.+epl*epl/(emn*emn)) : 0.);
      hjav0->SetBinContent(i, av);
      hjav0->SetBinError(i, eav);
    }
  }

  // for second estimate calculate ratio of Ansatzes with shifted pT
  TH1D *hjpl = (TH1D*)hzr->Clone(Form("hjec_%s_pl", jt));
  hjpl->Reset();
  TH1D *hjmn = (TH1D*)hzr->Clone(Form("hjec_%s_mn", jt));
  hjmn->Reset();
  TH1D *hjav = (TH1D*)hzr->Clone(Form("hjec_%s_av", jt));
  hjav->Reset();

  string s = Form("CondFormats/JetMETObjects/data/%s_%s_Uncertainty_%sPF.txt",
                  jp::jecgt, jp::type, jp::algo);
  cout << s << endl << flush;
  JetCorrectionUncertainty *rjet = new JetCorrectionUncertainty(s.c_str());

  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {

    if (hzr->GetBinContent(i)!=0) {

	double xmin = hzr->GetBinLowEdge(i);
	double xmax = hzr->GetBinLowEdge(i+1);

	double unc = hunc->GetBinContent(i);
	{
	  double eta = 0.5*(etamin+etamax);
	  func->setJetPt(hzr->GetBinCenter(i));
	  func->setJetEta(eta);
	  double u1 = func->getUncertainty(true);
	  func->setJetPt(hzr->GetBinCenter(i));
	  func->setJetEta(-eta);
	  double u2 = func->getUncertainty(true);
	  unc = 0.5*(u1+u2);
	}

	// Estimate the JEC uncertainty
	assert(jectype=="tot" || jectype=="abs" || jectype=="rel" ||
	       jectype=="bjt" || jectype=="pt");
	double djec = 0.;
	if (jectype=="abs") djec = unc;//(_ismc ? 0.01 : unc);
	if (jectype=="rel") djec = unc;//(_ismc ? 0.00 : 0.00);
	if (jectype=="tot") djec = unc;//(_ismc ? 0.01 : unc);
	if (jectype=="pt") {
	  rjet->setJetPt(hzr->GetBinCenter(i));
	  rjet->setJetEta(0.);
	  djec = rjet->getUncertainty(true);
	}
	// Add low pt uncertainty for MC
	if (_ismc && (jectype=="abs" || jectype=="tot")) {
	  djec = tools::oplus(djec, 0.04*pow(hzr->GetBinCenter(i)/10.,-2));
	}

	double djec0 = djec;
	double djec1 = djec;
	if (type==""); // suppress warning

	assert(xmax!=xmin);
	double ypl = fpt->Integral(xmin/(1+djec0), xmax/(1+djec1)) / (xmax-xmin);
	double ymn = fpt->Integral(xmin*(1+djec0), xmax*(1+djec1)) / (xmax-xmin);
	double yzr = fpt->Integral(xmin, xmax) / (xmax-xmin);

	// sanity checks
	//cout << Form("xmin %1.3g xmax %1.3g", xmin, xmax) << endl;
	//cout << Form("ymin %1.3g ymax %1.3g yzr %1.3g",
	//	 fpt->Eval(xmin), fpt->Eval(xmax), yzr) << endl;
	if (yzr > fpt->Eval(xmin) || yzr < fpt->Eval(xmax)) {

	  double ymin = fpt->Eval(xmin);
	  double ymax = fpt->Eval(xmax);
	  if (ymin < ymax)
	    cerr << Form("Warning: range [%3.3g, %3.3g] is not falling,",
			 xmin, xmax);
	  else
	    cerr << Form("Warning: range [%3.3g, %3.3g] may not be monotonous,",
			 xmin, xmax);
	  cerr << Form(" yzr=%1.3g, ymin=%1.3g, ymax=%1.3g",yzr,ymin,ymax)<<endl;
	}

	if (TMath::IsNaN(yzr) || yzr<0) yzr = 0;
	hjpl->SetBinContent(i, yzr ? ypl / yzr - 1 : 0.);
	hjmn->SetBinContent(i, yzr ? ymn / yzr - 1 : 0.);
	hjav->SetBinContent(i, yzr ? 0.5*(fabs(ypl/yzr-1)+fabs(ymn/yzr-1)) : 0.);
    }
  } // for i


  dzr->cd();

  return ( new sysc(hjpl, hjmn, hjav) );
} // jec_systematics


// Estimate impact of know JEC shifts (from MC truth, jet matching)
// in an ad-hoc way
void jec_shifts(TDirectory *dzr, TDirectory *dout, string type, string algo) {
  //,
  //	bool is38x) {

  float etamin, etamax;
  assert(sscanf(dzr->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);

  const int neta = 8;
  const int npar = 5;
// Produced by jetmatching.C on 20110224 213705;
// Note that parameters are multiplied by 100
const double p_ak5calo[8][5] =
  {{ 99.4, -2.27,  0.61,  0.261, -0.049}, // 1.1
   { 99.5, -1.30,  1.29, -0.256, -0.027}, // 0.8
   { 99.4, -2.24,  1.23,  0.219, -0.072}, // 1.2
   { 98.9, -4.04,  1.85,  0.866, -0.286}, // 1.9
   { 99.3, -2.33,  2.26,  0.391, -0.227}, // 2.1
   { 99.1,  0.31,  0.19, -0.706,  0.006}, // 0.3
   {103.5,  2.51, -4.28,  0.983,  2.349}, // 0.5
   { 98.0, -0.58, -2.01,  0.712,  0.932}}; // 1.8
const double p_ak5jpt[8][5] =
  {{ 97.2, -0.46,  0.17,  0.581, -0.157}, // 1.9
   { 96.5, -0.01,  0.85,  0.491, -0.341}, // 0.9
   { 96.1, -0.00,  0.52,  0.309, -0.092}, // 1.4
   { 96.3, -2.34,  1.20,  0.767, -0.317}, // 1.3
   { 96.2, -0.09,  1.49,  0.243, -0.378}, // 1.4
   { 96.4,  2.08, -0.54, -0.325,  0.041}, // 0.3
   {104.3,  4.71, -3.19,  1.051,  1.842}, // 0.4
   { 98.6,  2.94, -1.79, -0.304,  0.370}}; // 2.4

  // Select right table for the job
  int ieta = int(0.5*(etamin+etamax));
  const double (*p)[npar](0);
  if (type=="DATA") {
    if (algo=="calo") p = &p_ak5calo[0];
    if (algo=="jpt") p = &p_ak5jpt[0];
  }

  TF1 *fref = new TF1("fref","[0]+log(0.01*x)*([1]+log(0.01*x)*([2]+log(0.01*x)"
		      "*([3]+log(0.01*x)*[4])))", 0., 2000.);
  fref->SetParameters(1., 0., 0., 0., 0.);
  if (p) {
    assert(fref->GetNpar()==npar);
    for (int i = 0; i != npar; ++i) {
      fref->SetParameter(i, 0.01*p[min(ieta,neta-1)][i]);
    }
  }

  // inclusive jets
  TH1D *hzr = (TH1D*)dzr->Get("hpt"); assert(hzr);

  // make sure new histograms get created in the output file
  dout->cd();

  // inclusive jets
  TF1 *fpt0 = (TF1*)dzr->Get("fus"); assert(fpt0); fpt0->SetName("fpt0");
  TF1 *fpt = new TF1("fpt1","[0]*exp([1]/(x*(1-[6])))*pow(x*(1-[6]),[2])"
		     "*pow(1-x*(1-[6])*cosh([4])/[5], [3])",
		     jp::recopt, jp::emax);
  fpt->SetParameters(fpt0->GetParameter(0), fpt0->GetParameter(1),
		     fpt0->GetParameter(2), fpt0->GetParameter(3),
		     fpt0->GetParameter(4), fpt0->GetParameter(5),
		     0.);

  // for second estimate calculate ratio of Ansatzes with shifted pT
  TH1D *hs = (TH1D*)hzr->Clone(Form("hptshift_%s%s", algo.c_str(),""));
  hs->Reset();

  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {

    double xmin = hzr->GetBinLowEdge(i);
    double xmax = hzr->GetBinLowEdge(i+1);

    // Use function
    double djec0 = fref->Eval(xmin) - 1.;
    double djec1 = fref->Eval(xmax) - 1.;

    assert(xmax!=xmin);
    double ys = fpt->Integral(xmin/(1+djec0), xmax/(1+djec1)) / (xmax - xmin);
    double yzr = fpt->Integral(xmin, xmax) / (xmax - xmin);

    if (TMath::IsNaN(yzr) || yzr<0) yzr = 0;
    hs->SetBinContent(i, yzr ? ys / yzr - 1 : 0.);
  } // for i

  dzr->cd();

  return;
} // jec_shifts


sysc *jer_systematics(TDirectory *din, TDirectory *dout,
		      string type, string jertype) {

  assert(type=="MC" || type=="HW" || type=="DATA" || type=="TH");
  assert(jertype=="inc" || jertype=="bjt");
  bool _b = (jertype=="bjt");
  bool _th = (type=="TH");

  float etamin, etamax;
  assert(sscanf(din->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);

  TH1D *hzr = (TH1D*)din->Get(_th ? "hnlo" : "hpt"); assert(hzr);

  TF1 *fpt0 = (TF1*)din->Get(_th ? "fnlo" : "fus"); assert(fpt0);
  fpt0->SetName("fpt0");
  TF1 *fpt = new TF1("fpt2","[0]*exp([1]/x)*pow(x,[2])"
		     "*pow(1-x*cosh([4])/[5].,[3])",xmin,xmax);
  fpt->SetParameters(fpt0->GetParameter(0), fpt0->GetParameter(1),
		     fpt0->GetParameter(2), fpt0->GetParameter(3),
		     fpt0->GetParameter(4), fpt0->GetParameter(5));

  TF1 *fs = new TF1("fs2",smearedAnsatz,xmin,xmax,6);
  fs->SetParameters(fpt->GetParameter(0), fpt->GetParameter(1),
		    fpt->GetParameter(2), fpt->GetParameter(3),
		    fpt->GetParameter(4), 0.);

  // make sure new histograms get created in the output file
  dout->cd();

  // shift kernel resolutions up and down 20%
  TH1D *hspl = (TH1D*)hzr->Clone(_b ? "hjer_bjt_pl" : "hjer_pl"); hspl->Reset();
  TH1D *hsmn = (TH1D*)hzr->Clone(_b ? "hjer_bjt_mn" : "hjer_mn"); hsmn->Reset();
  TH1D *hsav = (TH1D*)hzr->Clone(_b ? "hjer_bjt_av" : "hjer_av"); hsav->Reset();

  _epsilon = 1e-4; // default 1e-12
  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {

    double ptmin = hzr->GetBinLowEdge(i);
    double ptmax = hzr->GetBinLowEdge(i+1);

    if (ptmin>=xmin && ptmax<=xmax && hzr->GetBinContent(i)!=0) {

      double y = fpt->Integral(ptmin, ptmax) / (ptmax - ptmin);
      double x = fpt->GetX(y, ptmin, ptmax);

      // sanity checks
      //cout << Form("xmin %1.0f xmax %1.0f x %1.0f", xmin, xmax, x) << endl;
      //if (x < xmin || x > xmax) {

      //cerr << Form("bin center x=%1.3f not in bin [%1.3g,%1.3g]"
      //	   ", reverting to bin midpoint", x, xmin, xmax) << endl;
      //x = 0.5*(xmin + xmax);
      assert(x >= ptmin);
      assert(x <= ptmax);
      //}

      // Estimate the JER uncertainty
      double djer = 0.;//0.10;//0.20;
      if (jertype=="inc") djer = (_ismc ? 0.05 : 0.10);
      // above overridden later

      // Over-ride with Jet Algorithm group's recommendations (interpolated)
      if (jertype=="inc" && !_ismc) {
	int iy = min(int((etamin+0.001)/0.5),5);
	djer = kpar[iy][1];
      }

      fs->SetParameter(5, +djer);
      double ypl = fs->Eval(x);
      fs->SetParameter(5, -djer);
      double ymn = fs->Eval(x);
      fs->SetParameter(5, 0.);
      double yzr = fs->Eval(x);

      hspl->SetBinContent(i, yzr ? max(ypl / yzr - 1, 0.) : 0.);
      hsmn->SetBinContent(i, yzr ? min(ymn / yzr - 1, 0.) : 0.);
      hsav->SetBinContent(i, yzr && ypl && ymn ?
			  0.5*(fabs(ypl/yzr-1) + fabs(ymn/yzr-1)) : 0);
    }
  } // for i


  din->cd();

  return ( new sysc(hspl, hsmn, hsav) );
} // jer_systematics


sysc *lum_systematics(TDirectory *din, TDirectory *dout) {

  TH1D *hzr = (TH1D*)din->Get("hpt"); assert(hzr);

  // make sure new histograms get created in the output file
  dout->cd();
  TH1D *hlpl = (TH1D*)hzr->Clone("hlum_pl"); hlpl->Reset();
  TH1D *hlmn = (TH1D*)hzr->Clone("hlum_mn"); hlmn->Reset();
  TH1D *hlav = (TH1D*)hzr->Clone("hlum_av"); hlav->Reset();

  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {

    double lumsys = (_ismc ? 0. : 0.022);
    if (hzr->GetBinContent(i)!=0) {
      hlpl->SetBinContent(i, +lumsys);
      hlmn->SetBinContent(i, -lumsys);
      hlav->SetBinContent(i, lumsys);
    }
  }

  din->cd();

  return ( new sysc(hlpl, hlmn, hlav) );
} // lum_systematics

void tot_systematics(TDirectory *din, TDirectory *dout,
		     sysc *cjec1, sysc *cjec2, sysc *cjec3,
		     sysc *cjer1, sysc *cjer2,
		     sysc *clum) {

  TDirectory *curdir = gDirectory;

  TH1D *hpt = (TH1D*)din->Get("hpt"); assert(hpt);
  TGraphErrors *gpt = (TGraphErrors*)din->Get("gpt");
  if (!gpt) {
    cout << "gpt not found in " << din->GetName() << endl;
    assert(gpt);
  }

  string ss = "";
  const char *s = ss.c_str();

  // Create new graphs in the output
  dout->cd();
  TGraphErrors *gpt_stat = (TGraphErrors*)gpt->Clone();
  gpt_stat->SetName(Form("gpt%s_stat",s));
  TGraphErrors *gpt_syst = (TGraphErrors*)gpt->Clone();
  gpt_syst->SetName(Form("gpt%s_syst",s));

  TH1D *htpl = (TH1D*)cjec1->plus->Clone(Form("htot%s_pl",s)); htpl->Reset();
  TH1D *htmn = (TH1D*)cjec1->minus->Clone(Form("htot%s_mn",s)); htmn->Reset();
  TH1D *htav = (TH1D*)cjec1->av->Clone(Form("htot%s_av",s)); htav->Reset();

  for (int i = 1; i != cjec1->av->GetNbinsX()+1; ++i) {

    double jec1_pl = cjec1->plus->GetBinContent(i);
    double jec1_mn = cjec1->minus->GetBinContent(i);
    double jec2_pl = cjec2->plus->GetBinContent(i);
    double jec2_mn = cjec2->minus->GetBinContent(i);
    double jer1_pl = cjer1->plus->GetBinContent(i);
    double jer1_mn = cjer1->minus->GetBinContent(i);
    double lum_pl = clum->plus->GetBinContent(i);

    double tot_pl = sqrt(pow(max(jec1_pl,0.),2) + pow(max(jec1_mn,0.),2) +
			 pow(max(jec2_pl,0.),2) + pow(max(jec2_mn,0.),2) +
			 pow(max(jer1_pl,0.),2) + pow(max(jer1_mn,0.),2) +
			 lum_pl*lum_pl);
    double tot_mn = sqrt(pow(min(jec1_pl,0.),2) + pow(min(jec1_mn,0.),2) +
			 pow(min(jec2_pl,0.),2) + pow(min(jec2_mn,0.),2) +
			 pow(min(jer1_pl,0.),2) + pow(min(jer1_mn,0.),2) +
			 lum_pl*lum_pl);
    double tot_av = 0.5*(tot_pl + tot_mn);

    if (hpt->GetBinContent(i)!=0) {
      htpl->SetBinContent(i, tot_pl);
      htmn->SetBinContent(i, -tot_mn);
      htav->SetBinContent(i, tot_av);
    }

    // Not the most efficient way to find a point, but should work...
    double xmin = cjec1->plus->GetBinLowEdge(i);
    double xmax = cjec1->plus->GetBinLowEdge(i+1);
    for (int j = 0; j != gpt_syst->GetN(); ++j) {
      double x, y;
      gpt_syst->GetPoint(j, x, y);
      if (x>=xmin && x<xmax && hpt->GetBinContent(hpt->FindBin(x))!=0) {
	double ex = gpt_syst->GetErrorX(j);
	gpt_syst->SetPointError(j, ex, tot_av);
      }
    } // for j
  } // for i

  dout->Add(gpt_stat);
  dout->Add(gpt_syst);

  curdir->cd();
} // tot_systematics


void statistics(TDirectory *din, TDirectory *dout) {

  TH1D *hzr = (TH1D*)din->Get("hpt"); assert(hzr);

  // make sure new histograms get created in the output file
  dout->cd();

  TF1 *fpt0 = (TF1*)din->Get("fus"); assert(fpt0); fpt0->SetName("fpt0");
  TF1 *fpt = new TF1("fpt3","[0]*exp([1]/x)*pow(x,[2])"
		     "*pow(1-x*cosh([4])/[5],[3])",
		     jp::recopt, jp::emax);
  fpt->SetParameters(fpt0->GetParameter(0), fpt0->GetParameter(1),
		     fpt0->GetParameter(2), fpt0->GetParameter(3),
		     fpt0->GetParameter(4), fpt0->GetParameter(5));

  TF1 *fs = new TF1("fs3",smearedAnsatz, jp::recopt, jp::emax, 6);
  fs->SetParameters(fpt->GetParameter(0), fpt->GetParameter(1),
		    fpt->GetParameter(2), fpt->GetParameter(3),
		    fpt->GetParameter(4), 0.);

  // Estimate statistical uncertainty with 10 mub-1 and 1pb-1
  TH1D *hmub = (TH1D*)hzr->Clone("hsta_mub"); hmub->Reset();
  TH1D *hnb = (TH1D*)hzr->Clone("hsta_nb"); hnb->Reset();
  TH1D *hpb = (TH1D*)hzr->Clone("hsta_pb"); hpb->Reset();

  _epsilon = 1e-4; // default 1e-12
  for (int i = 1; i != hzr->GetNbinsX()+1 && false; ++i) {

    double xmin = hzr->GetBinLowEdge(i);
    double xmax = hzr->GetBinLowEdge(i+1);
    double y = fpt->Integral(xmin, xmax) / (xmax - xmin);
    double x = fpt->GetX(y, xmin, xmax);

    // sanity checks
    assert(x >= xmin);
    assert(x <= xmax);

    double nd = fs->Eval(x)*(xmax-xmin); // 1/bin/pb-1
    double stat_mub = 1./sqrt(nd*0.0001);
    double stat_nb  = 1./sqrt(nd*0.01);
    double stat_pb  = 1./sqrt(nd*1.);

    if (stat_mub<1.) hmub->SetBinContent(i, stat_mub);
    if (stat_nb<1.) hnb->SetBinContent(i, stat_nb);
    if (stat_pb<1.) hpb->SetBinContent(i, stat_pb);
  } // for i

  din->cd();

} // statistics


void sources(string type = "DATA") {

  assert(type=="DATA" || type=="MC" || type=="HW");
  const char *t = type.c_str();

  cout << "*** Determining uncertainty correlations" << endl;
  //cout << "!!! NB: Running over GRV23_AK7 directory" << endl;

  // NLO theory prediction with NP corrections as baseline
  TFile *fth = new TFile(Form("output-%s-5.root",t),"READ");
  assert(fth && !fth->IsZombie());

  assert(fth->cd("Standard"));
  TDirectory *dth0 = gDirectory;

  TFile *fout = new TFile(Form("output-%s-4.root",t), "UPDATE");
  assert(fout && !fout->IsZombie());

  fout->mkdir("Standard");
  assert(fout->cd("Standard"));
  TDirectory *dout0 = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = dth0->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory
    if (obj->InheritsFrom("TDirectory")) {

      assert(dth0->cd(obj->GetName()));
      TDirectory *dth = gDirectory;

      dout0->mkdir(obj->GetName());
      assert(dout0->cd(obj->GetName()));
      TDirectory *dout = gDirectory;

      sourceBin(dth, dout);
    } // inherits TDirectory
  } // while

  cout << "Output stored in " << fout->GetName() << endl;
  fout->Write();
  fout->Close();
  fout->Delete();

  fth->Close();
  fth->Delete();
} // systematics

// Calculate uncertainty sources based on JEC
void sourceBin(TDirectory *dth, TDirectory *dout) {

  float etamin, etamax;
  assert(sscanf(dth->GetName(),"Eta_%f-%f",&etamin,&etamax)==2);

  // inclusive jets
  TH1D *hnlo = (TH1D*)dth->Get("hnlo"); assert(hnlo);
  TF1 *fnlo0 = (TF1*)dth->Get("fnlo"); assert(fnlo0); fnlo0->SetName("fnlo0");
  TF1 *fnlo = new TF1("fnlo","[0]*exp([1]/x)*pow(x,[2])"
		     "*pow(1-x*cosh([4])/[5], [3])",
		      jp::recopt, jp::emax);
  fnlo->SetParameters(fnlo0->GetParameter(0), fnlo0->GetParameter(1),
		      fnlo0->GetParameter(2), fnlo0->GetParameter(3),
		      fnlo0->GetParameter(4), fnlo0->GetParamater(5));

  // make sure new histograms get created in the output file
  dout->cd();

  string srcnames[] =
    {"Absolute",  "HighPtExtra", "SinglePion", "Flavor", "Time",
     "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
     "RelativeStatEC2", "RelativeStatHF", "RelativeFSR",
     "PileUpDataMC", "PileUpOOT", "PileUpPt", "PileUpBias", "PileUpJetRate"};
  const int nsrc = sizeof(srcnames)/sizeof(srcnames[0]);

  TH1D *hup0 = (TH1D*)hnlo->Clone("src_up"); hup0->Reset();
  TH1D* hdw0 = (TH1D*)hnlo->Clone("src_dw"); hdw0->Reset();

  TH1D *horig0 = (TH1D*)hnlo->Clone("orig"); horig0->Reset();

  for (int isrc = 0; isrc != nsrc; ++isrc) {

    string s = Form("CondFormats/JetMETObjects/data/%s_%s_UncertaintySources_%sPF.txt",
                    jp::jecgt, jp::type, jp::algo);
    cout << s << endl << flush;
    JetCorrectorParameters *p = new JetCorrectorParameters(s.c_str(),
							   srcnames[isrc]);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);

    TH1D *hup = (TH1D*)hnlo->Clone(Form("src_%d_up",isrc)); hup->Reset();
    TH1D* hdw = (TH1D*)hnlo->Clone(Form("src_%d_dw",isrc)); hdw->Reset();
    TH1D *horig = (TH1D*)hnlo->Clone(Form("orig_%d",isrc)); horig->Reset();

    hup->SetTitle(srcnames[isrc].c_str());
    hdw->SetTitle(srcnames[isrc].c_str());
    horig->SetTitle(srcnames[isrc].c_str());

    for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

      if (hnlo->GetBinContent(i)!=0) {

	double xmin = hnlo->GetBinLowEdge(i);
	double xmax = hnlo->GetBinLowEdge(i+1);

	double pt = 0.5*(xmin+xmax);
	double eta = 0.5*(etamin+etamax);
	unc->setJetPt(pt);
	unc->setJetEta(eta);
	double u1 = unc->getUncertainty(true);
	unc->setJetPt(pt);
	unc->setJetEta(-eta);
	double u2 = unc->getUncertainty(true);
	double djec = 0.5*(u1+u2);

	double djec0 = djec;
	double djec1 = djec;

	assert(xmax!=xmin);
	double ypl = fnlo->Integral(xmin/(1+djec0), xmax/(1+djec1))/(xmax-xmin);
	double ymn = fnlo->Integral(xmin*(1+djec0), xmax*(1+djec1))/(xmax-xmin);
	double yzr = fnlo->Integral(xmin, xmax) / (xmax-xmin);

	// sanity checks
	double ymin = fnlo->Eval(xmin);
	double ymax = fnlo->Eval(xmax);
	if (TMath::IsNaN(yzr) || yzr<0) {
	  cerr << Form("Warning: yzer=%1.3g for range [%3.3g, %3.3g]",
		       yzr, xmin, xmax);
	  yzr = 0;
	}
	else if (yzr > ymin || yzr < ymax) {

	  if (ymin < ymax)
	    cerr << Form("Warning: range [%3.3g, %3.3g] is not falling,",
			 xmin, xmax);
	  else
	    cerr << Form("Warning: range [%3.3g, %3.3g] may not be monotonous,",
			 xmin, xmax);
	  cerr << Form(" yzr=%1.3g, ymin=%1.3g, ymax=%1.3g",
		       yzr, ymin, ymax) << endl;
	}

	hup->SetBinContent(i, yzr ? ypl / yzr - 1 : 0.);
	hdw->SetBinContent(i, yzr ? ymn / yzr - 1 : 0.);
	horig->SetBinContent(i, djec);

	// Add sources in quadrature
	double eup = (yzr ? max(ypl/yzr-1, ymn/yzr-1) : 0);
	double edw = (yzr ? min(ypl/yzr-1, ymn/yzr-1) : 0);
	hup0->SetBinContent(i, +tools::oplus(hup0->GetBinContent(i), eup));
	hdw0->SetBinContent(i, -tools::oplus(hdw0->GetBinContent(i), edw));
	horig0->SetBinContent(i, tools::oplus(horig0->GetBinContent(i), djec));
      }
    } // for i
  } // for isrc

  // Lumi source
  TH1D *hlup = (TH1D*)hnlo->Clone(Form("src_%d_up",nsrc)); hlup->Reset();
  TH1D *hldw = (TH1D*)hnlo->Clone(Form("src_%d_dw",nsrc)); hldw->Reset();
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {
    hlup->SetBinContent(i, 1./(1-0.022)-1);
    hldw->SetBinContent(i, 1./(1+0.022)-1);
  }
  hlup->SetTitle("Luminosity");
  hldw->SetTitle("Luminosity");

  // Resolution (unfolding) source
  sysc *sys = jer_systematics(dth, dout, "TH", "inc");
  dout->cd();

  TH1D *hjerup = (TH1D*)hnlo->Clone(Form("src_%d_up",nsrc+1)); hjerup->Reset();
  TH1D *hjerdw = (TH1D*)hnlo->Clone(Form("src_%d_dw",nsrc+1)); hjerdw->Reset();
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {
    int j = sys->av->FindBin(hnlo->GetBinCenter(i));
    hjerup->SetBinContent(i, sys->plus->GetBinContent(j));
    hjerdw->SetBinContent(i, sys->minus->GetBinContent(j));
  }
  hjerup->SetTitle("Unfolding");
  hjerdw->SetTitle("Unfolding");

  dth->cd();
} // sources
