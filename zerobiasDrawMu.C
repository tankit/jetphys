#include "tdrstyle_mod15.C"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLine.h"

#include "TLatex.h"

#include <vector>
#include <string>

using namespace std;

void zerobiasDrawMu(string mode="raw", string type="MC") {

  assert(mode=="raw" || mode=="corr" || mode=="rc" || mode=="mc");
  assert(type=="MC" || type=="DATA");

  setTDRStyle();

  //TFile *fin = new TFile("rootfiles/zerobias-data-v17c.root","READ");
  //TFile *find = new TFile("rootfiles/zerobias-data-v18.root","READ");
  TFile *find = new TFile("rootfiles/zerobias-data-v18v14fix.root","READ");
  assert(find && !find->IsZombie());
  TFile *finm = new TFile("rootfiles/zerobias-mc-v15.root","READ");
  assert(finm && !finm->IsZombie());
  TFile *fin = (type=="DATA" ? find : finm);

  int vpt[] = {10, 15, 20, 25, 30, 40, 60, 80, 120};
  const int npt = sizeof(vpt)/sizeof(vpt[0]);
  TH1D* vh[npt];
  TH1D* vh2[npt];

  int marker[npt] = {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle,
		     kFullDiamond, kOpenDiamond, kDot, kDot, kDot};
  int color[npt] = {kRed, kRed, kGreen+2, kGreen+2,
		    kBlue, kBlue, kBlack, kBlack, kBlack};

  TH1D *hmu0 = (TH1D*)fin->Get("hmu0"); assert(hmu0);
  hmu0->Rebin(2);
  hmu0->Scale(1./hmu0->Integral());

  for (int ipt = 0; ipt != npt; ++ipt) {

    int pt = vpt[ipt];

    TH1D *h = (TH1D*)fin->Get(Form("ppt%sovptchsvsmu_pt%d",mode.c_str(),pt));
    assert(h);
    vh[ipt] = h;

    string mode2 = mode;
    if (mode=="corr") mode2 = "res";
    TH1D *h2 = (TH1D*)fin->Get(Form("h%svsmu_pt%d",mode2.c_str(),pt));
    h2->Rebin(2);
    h2->Scale(1./h2->Integral());
    h2->Divide(hmu0);
    assert(h2);
    vh2[ipt] = h2;
  }

  TH1D *hup = new TH1D("hup",Form(";#mu;0.6#timesp_{T,%s}/p_{T,ch}",mode.c_str()),45,0,45);
  hup->SetMaximum(1.65);
  hup->SetMinimum(0.65);
  TH1D *hdw = new TH1D("hdw",";#mu;Ratio",45,0,45);
  hdw->SetMaximum(1.5);
  hdw->SetMinimum(0.5);

  lumi_13TeV = (mode=="MC" ? "2016 Simulation V15" : "Run2016GH V18 16.8 fb^{-1}");
  //if (mode=="MC") extraText = "Simulation Preliminary";

  //TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);//33);//11);
  TCanvas *c1 = tdrCanvas("c1",hup,4,11,kSquare);
  
  c1->cd(1);

  //TLegend *leg = tdrLeg(0.17,0.45,0.37,0.90);
  //TLegend *leg = tdrLeg(0.17,0.45,0.37,0.75); // DiCanvas
  TLegend *leg = tdrLeg(0.17,0.56,0.37,0.80); // Canvas
  leg->SetTextSize(0.040);

  for (int ipt = 0; ipt != npt; ++ipt) {

    int ptmin = vpt[ipt];
    int ptmax = (ipt!=npt-1 ? vpt[ipt+1] : 6500.);

    //if (ptmin>60) continue;
    if (ptmin>40) continue;

    TH1D *h = vh[ipt];

    TF1 *f1 = new TF1(Form("f1_%d",ptmin),"[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,45);
    if (ptmin>20) f1->FixParameter(3,0);
    f1->SetLineColor(color[ipt]);
    f1->SetLineStyle(ipt%2==0 ? kSolid : kDotted);
    //f1->SetLineStyle(ipt%2==0 ? kDotted : kSolid);
    f1->SetLineWidth(2);

    h->Fit(f1,"QRN");
    f1->Draw("SAME");

    TF1 *f1b = new TF1(Form("f1b_%d",ptmin),"[0]+[1]*x",0,45);
    f1b->SetLineColor(color[ipt]);
    f1b->SetLineStyle(ipt%2==0 ? kDotted : kSolid);

    h->Fit(f1b,"QRN");
    //f1b->Draw("SAME");

    double mu = 25;
    double ptmid = 0.5*(ptmin+ptmax);
    double off = (f1->Eval(mu)/f1->Eval(0) - 1)/(0.7*mu) * ptmid;
    tdrDraw(h,"P",marker[ipt],color[ipt]);
    if (ipt!=npt-1) leg->AddEntry(h,Form("%d<p_{T,ch}/0.6<%d (%1.2f)",ptmin,ptmax,off),"PL");
    else            leg->AddEntry(h,Form("%d<p_{T,ch}/0.6 (%1.2f)",ptmin,off),"PL");
  }

  c1->cd(2);


  //c1->SaveAs(Form("pdf/zerobiasDrawMu_bal_%s.pdf",mode.c_str()));
  c1->SaveAs(Form("pdf/zerobiasDrawMu_bal_%s_%s-v14fix.pdf",
		  type.c_str(),mode.c_str()));


  TH1D *h2 = new TH1D("h2",";#mu;Normalized N_{jet} / #mu (a.u.)",45,0,45);
  h2->SetMaximum(1.65);
  h2->SetMinimum(0.65);

  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);

  const int nskip2 = (mode=="raw" ? (type=="MC" ? 4 : 4) :
		      (type=="MC" ? 3 : 3));

  TLegend *leg2 = tdrLeg(0.17,0.56,0.37,0.80); // Canvas
  leg2->SetTextSize(0.040);


  for (int ipt = 0; ipt != npt; ++ipt) {

    int ptmin = vpt[ipt];
    int ptmax = (ipt!=npt-1 ? vpt[ipt+1] : 6500.);

    //if (ptmin>40) continue;
    //if(ptmin<20) continue;
    if (ipt<nskip2) continue;
    if (mode=="raw" ? (type=="MC" ? ptmin>60 : ptmin>80) :
	(type=="MC" ? ptmin>60 : ptmin>60 )) continue;

    TH1D *h = vh2[ipt];
    
    // Clean out biggest error bars
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      if (h->GetBinError(i)>0.1 || (mode!="raw" && h->GetBinContent(i)<0.9)) {
	h->SetBinContent(i,0);
	h->SetBinError(i,0);
      }
    }

    TF1 *f1 = new TF1(Form("f2_%d",ptmin),"[0]+[1]*x+abs([2])*x*x-[3]*x*x*x",0,45);
    //if (ptmin>40) f1->FixParameter(3,0);
    //if (ptmin>30) f1->FixParameter(2,0);
    if (mode=="raw" && type=="MC" && ptmin>40) f1->FixParameter(3,0);
    f1->SetLineColor(color[ipt-nskip2]);
    f1->SetLineStyle((ipt-nskip2)%2==0 ? kSolid : kDotted);
    //f1->SetLineStyle(ipt%2==0 ? kDotted : kSolid);
    f1->SetLineWidth(1);//2);

    h->Fit(f1,"QRN");
    f1->Draw("SAME");

    TF1 *f1b = new TF1(Form("f2b_%d",ptmin),"[0]+[1]*x",0,45);
    f1b->SetLineColor(color[ipt-nskip2]);
    f1b->SetLineStyle(ipt%2==0 ? kDotted : kSolid);

    h->Fit(f1b,"QRN");
    //f1b->Draw("SAME");

    double mu = 25;
    double ptmid = 0.5*(ptmin+ptmax);
    double off = (f1->Eval(mu)/f1->Eval(0) - 1)/(0.7*mu) * ptmid;
    tdrDraw(h,"P",marker[ipt-nskip2],color[ipt-nskip2]);
    if (ipt!=npt-1) leg2->AddEntry(h,Form("%d<p_{T,%s}<%d (%1.2f)",ptmin,mode.c_str(),ptmax,off),"PL");
    else            leg2->AddEntry(h,Form("%d<p_{T,%s} (%1.2f)",ptmin,mode.c_str(),off),"PL");


    // Estimate impact from JER
    // f(pT)=N0*pT^-a, a~5; approx. f'(pT)=N0*exp(-alpha*pT)
    // df/dpT=df'/dpT => alpha = a/pT
    // sigma=const => F'(pT)=f'(pT)*exp(a^2*sigma^2/2)
    // F(pT) ~ f(pT)*exp(a^2/2 * (sigma/pT)^2)
    // (sigma/pT)^2 = (sigma0/pT)^2 + mu*sigmaPU^2/pT^2
    // => F(pT) = f(pT) * exp(a^2/2*(sigma/pT)^2) * exp((a*sigmaPU/pT)^2/2 * mu)
    // => F(pT) = const1 * exp( const2 * mu),
    // where const2 ~ 0.5(5*sigmaPU/pT)^2 
    TF1 *fjer = new TF1(Form("fjer_%d",ptmin),"exp([0]*x)/exp([0]*25.)",0,50);
    // Chad Harrington, MC truth JER:
    // https://indico.cern.ch/event/706518/#3-jet-energy-resolution-with-f
    // Feb 28, 2018, slide 5: JER=0.4 at pT=15 GeV for mu=60-70, and 0.18 at mu=0-10
    // => sigmaPU = sqrt( (0.40^2 - 0.18^2)/(65-5) ) * 15 = 0.69
    // or, JER=0.21 at mu=10-20 and 0.34 at mu=50-60
    // => sigmaPU = sqrt( (0.34^2 - 0.21^2)/(55-15) ) * 15 = 0.63
    // (mu:60-10 => 0.56);
    // or, JER=0.13 at pT=30 GeV and mu=0-10, 0.23 at mu=60-70
    // => sigmaPU = sqrt( (0.23^2 - 0.12^2)/(65-5) ) * 30 = 0.76
    //double sigmapu = 0.76;//0.56;//0.63;//0.69; 
    double sigmapu = (type=="DATA" ? 1.15 : 1.0)*0.66; //0.76 for data
    double ptavg = 0.5*(ptmin+ptmax);
    fjer->SetParameter(0, 0.5 * pow(5*sigmapu/ptavg,2));
    fjer->SetLineColor(color[ipt-nskip2]);//+2);
    fjer->SetLineStyle((ipt-nskip2)%2==0 ? kSolid : kDotted);
    fjer->SetLineWidth(3);
    fjer->Draw("SAME");

  } // for ipt

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.75,0.85,type.c_str());

  //c2->SaveAs(Form("pdf/zerobiasDrawMu_rate_%s.pdf",mode.c_str()));
  c2->SaveAs(Form("pdf/zerobiasDrawMu_rate_%s_%s-v14fix.pdf",
		  type.c_str(),mode.c_str()));

}
