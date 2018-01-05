#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include <map>

#include "tdrstyle_mod15.C"

using namespace std;

TH1D *_hq(0), *_hg(0);
Double_t fitTemplate(Double_t *x, Double_t *p) {

  assert(_hq);
  assert(_hg);

  double qgl = x[0];
  double fg = p[0];
  double fq = 1-fg;

  int iqgl = _hg->FindBin(qgl);
  double ng = _hg->GetBinContent(iqgl);
  double nq = _hq->GetBinContent(iqgl);

  return (fg*ng + nq*fq);
}

struct fitvals {
  double pt;
  double fggen;
  double fqgen;
  double fgreco;
  double fqreco;
  double fgdata;
  double fqdata;
} ;

// Run with 'root -l drawGluon.C+g'
fitvals drawQGLbin(double etamin = 0, double etamax = 1.3, int trg = 260) {
  //double ptmin = 330, double ptmax = 395) {

  //double etamin = 0.;
  //double etamax = 1.3;
  //double ptmin = 330;
  //double ptmax = 395;
  //int trg = 260;

  // Define pT ranges for each trigger (or take from settings.h)
  map<int, pair<int,int> > trigpts;
  trigpts[40]  = make_pair<int,int>(0,   84);
  trigpts[60]  = make_pair<int,int>(84,  114);
  trigpts[80]  = make_pair<int,int>(114, 196);
  trigpts[140] = make_pair<int,int>(196, 272);
  trigpts[200] = make_pair<int,int>(272, 330);
  trigpts[260] = make_pair<int,int>(330, 395);
  trigpts[320] = make_pair<int,int>(395, 468);
  trigpts[400] = make_pair<int,int>(468, 548);
  trigpts[450] = make_pair<int,int>(548, 6500);

  double ptmin = trigpts[trg].first;
  double ptmax = trigpts[trg].second;

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fd = new TFile("output-DATA_31Octlongprobp_RunGpart1-1.root","READ");
  assert(fd && !fd->IsZombie());
  TFile *fm = new TFile("output-MC_30Octlongprobg-1.root","READ");
  assert(fm && !fm->IsZombie());

  //string dir = "Standard/Eta_0.0-0.5/jt80";//mc";
  //string dir = "Standard/Eta_0.0-0.5/jt320";//mc";
  //string dir = "Standard/Eta_3.2-4.7/jt450";//450";//mc";
  //string dir = "Standard/Eta_0.0-1.3/jt260";
  string dir = Form("Standard/Eta_%1.1f-%1.1f/jt%d",etamin,etamax,trg);
  TH1D *hg = (TH1D*)fm->Get(Form("%s/hqgl_g",dir.c_str()));
  assert(hg);
  TH1D *hq = (TH1D*)fm->Get(Form("%s/hqgl_q",dir.c_str()));
  assert(hq);
  //TH1D *hg = (TH1D*)fm->Get(Form("%s/hgpt_g0tw",dir.c_str()));
  //assert(hg);
  //TH1D *hq = (TH1D*)fm->Get(Form("%s/hqpt_g0tw",dir.c_str()));                                                                        
  //assert(hq);
  //TH1D *ha = (TH1D*)fm->Get(Form("%s/hpt_g0tw",dir.c_str()));
  //assert(ha);
  //TH1D *hd = (TH1D*)fd->Get(Form("%s/hpt_g0tw",dir.c_str())); 
  //assert(hd);  
  TH1D *ha = (TH1D*)fm->Get(Form("%s/hqgl",dir.c_str()));
  assert(ha);
  TH1D *hd = (TH1D*)fd->Get(Form("%s/hqgl",dir.c_str()));
  assert(hd);

  curdir->cd();

  hg->Scale(1./ha->Integral());
  hq->Scale(1./ha->Integral());
  ha->Scale(1./ha->Integral());
  hd->Scale(1./hd->Integral());

  // Set MC shapes for template fit
  TH1D *hqm = (TH1D*)hq->Clone("hqm");
  hqm->Scale(1./hqm->Integral());
  TH1D *hgm = (TH1D*)hg->Clone("hgm");
  hgm->Scale(1./hgm->Integral());
  
  // Apply shaping to create data templates
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#Systematics
  TH1D *hqd = (TH1D*)hq->Clone("hqd");
  TF1 *fwq = new TF1("fwq","-0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581",
		 0,1);
  hqd->Multiply(fwq);
  hqd->Scale(1./hqd->Integral());
  TH1D *hgd = (TH1D*)hg->Clone("hgd");
  TF1 *fwg = new TF1("fwg","-55.7067*x^7 + 113.218*x^6 -21.1421*x^5 -99.927*x^4 + 92.8668*x^3 -34.3663*x^2 + 6.27*x + 0.612992",0,1);
  hgd->Multiply(fwg);
  hgd->Scale(1./hgd->Integral());

  // Add MC uncertainty to data to avoid biasing chi2
   for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
            hd->SetBinError(i, sqrt(pow(hd->GetBinError(i),2)+
  			    pow(ha->GetBinError(i),2)));
   }

  // Do a template fit
  TF1 *f1m = new TF1("f1m",fitTemplate,0,1,1);
  f1m->SetLineWidth(2);
  f1m->SetNpx(1000.);
  _hq = hqm;
  _hg = hgm;
  ha->Fit(f1m,"N");

  TF1 *f1d = new TF1("f1d",fitTemplate,0,1,1);
  f1d->SetLineWidth(2);
  f1d->SetNpx(1000.);
  _hq = hqd;
  _hg = hgd;
  hd->Fit(f1d,"N");
  gStyle->SetOptFit(kFALSE);


  // Calculate MC truth quark and gluon fractions
  double efg(0), efq(0);
  double fg = hg->IntegralAndError(-1,hg->GetNbinsX()+1,efg);
  double fq = hq->IntegralAndError(-1,hq->GetNbinsX()+1,efq);

  TH1D *hsum = (TH1D*)hg->Clone("hsum");
  hsum->Add(hq);

  TH1D *hup = new TH1D("hup",";Quark-Gluon likelihood;1/#sigma d#sigma/dQGL",
		       100,0,1.00);
  hup->SetMinimum(0.);
  hup->SetMaximum(0.085);//0.06);

  TH1D *hdw = new TH1D("hdw",";Quark-Gluon likelihood;Data/MC",
		       100,0,1.00);
  hdw->SetMinimum(0.75);
  hdw->SetMaximum(1.25);


  lumi_13TeV = "Ozlem's analysis, 36 fb^{-1}";
  //second plot
  //second plot

  // also tdrDiCanvas
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  c1->cd(1);
  tdrDraw(hd,"HP",kFullCircle,kBlack);
  tdrDraw(ha,"P",kOpenDiamond,kBlack);
  tdrDraw(hgm,"P",kOpenCircle,kBlue);
  tdrDraw(hqm,"P",kOpenSquare,kRed);
  tdrDraw(hgd,"P",kOpenCircle,kBlue+1); hgd->SetMarkerSize(0.7);
  tdrDraw(hqd,"P",kOpenSquare,kRed+2);  hqd->SetMarkerSize(0.7);
  
  // Draw fits on top
  _hq = hqm;
  _hg = hgm;
  f1m->SetLineColor(kBlue);
  f1m->DrawClone("SAME");

  _hq = hqd;
  _hg = hgd;
  f1d->SetLineColor(kBlue);
  f1d->DrawClone("SAME");


  TLegend *leg = tdrLeg(0.60, 0.66, 0.85, 0.90);
  leg->AddEntry(hd,"Data","FP");
  leg->AddEntry(ha,"MC","P");
  leg->AddEntry(hgm,"Gluons (MC)","LP");
  leg->AddEntry(hqm,"Quarks (MC)","LP");
  //c2 
  //c2
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();

  //tex->DrawLatex(0.20,0.74,"3.2 < |y_{jet}| < 4.7");
  tex->DrawLatex(0.20,0.74,Form("%1.1f #leq |y_{jet}| < %1.1f",etamin,etamax));
  //tex->DrawLatex(0.20,0.68,Form("%1.0f < p_{T,jet} < %1.0f GeV",548.,6500.));
  tex->DrawLatex(0.20,0.68,Form("%1.0f < p_{T,jet} < %1.0f GeV",ptmin,ptmax));
  tex->DrawLatex(0.55,0.42,"anti-k_{T} R=0.4 PF+CHS");
  tex->DrawLatex(0.20,0.60,Form("f_{g}^{reco} = %1.3f #pm %1.3f",f1m->GetParameter(0),f1m->GetParError(0)));
  tex->DrawLatex(0.20,0.54,Form("f_{q}^{reco} = %1.3f #pm %1.3f",1-f1m->GetParameter(0), f1m->GetParError(0)));
  tex->DrawLatex(0.20,0.48,Form("f_{g}^{gen} = %1.3f #pm %1.3f",fg,efg));
  tex->DrawLatex(0.20,0.42,Form("f_{q}^{gen} = %1.3f #pm %1.3f",fq,efq));

  tex->DrawLatex(0.55,0.60,Form("f_{g}^{data} = %1.3f #pm %1.3f",f1d->GetParameter(0),f1d->GetParError(0)));
  tex->DrawLatex(0.55,0.54,Form("f_{q}^{data} = %1.3f #pm %1.3f",1-f1d->GetParameter(0), f1d->GetParError(0)));
  tex->DrawLatex(0.55,0.48,Form("#chi^{2} / NDF = %1.1f / %d",f1d->GetChisquare(),f1d->GetNDF()));

  gPad->RedrawAxis();

  c1->cd(2);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,1,1);

  TH1D *hrd = (TH1D*)hd->Clone("hrd");
  _hq = hqd;
  _hg = hgd;
  hrd->Divide(f1d);
  tdrDraw(hrd,"P",kFullCircle,kRed);

  TH1D *hrm = (TH1D*)ha->Clone("hrm");
  _hq = hqm;
  _hg = hgm;
  hrm->Divide(f1m);
  tdrDraw(hrm,"P",kOpenSquare,kBlue);

  //c1->SaveAs("pdf/drawQGL_hpt_eta3.2_4.7_450.pdf");
  c1->SaveAs(Form("pdf/drawQGL_hpt_eta%1.1f_%1.1f_%d.pdf",etamin,etamax,trg));

  if (gROOT->IsBatch()) {
    delete hup;
    delete hdw;
  }

  fitvals fv;
  fv.pt = 0.5*(ptmax+ptmin);
  fv.fggen = fg;
  fv.fqgen = fq;
  fv.fgreco = f1m->GetParameter(0);
  fv.fqreco = 1-f1m->GetParameter(0);
  fv.fgdata = f1d->GetParameter(0);
  fv.fqdata = 1-f1d->GetParameter(0);

  return fv;
}


#include "TGraphErrors.h"
#include <vector>

void drawQGL() {

  setTDRStyle();
  
  TGraphErrors *ggdata = new TGraphErrors(0);
  TGraphErrors *gqdata = new TGraphErrors(0);
  TGraphErrors *ggreco = new TGraphErrors(0);
  TGraphErrors *gqreco = new TGraphErrors(0);
  TGraphErrors *gggen = new TGraphErrors(0);
  TGraphErrors *gqgen = new TGraphErrors(0);

  std::vector<int> trgs;
  trgs.push_back(40);
  trgs.push_back(60);
  trgs.push_back(80);
  trgs.push_back(140);
  trgs.push_back(200);
  trgs.push_back(260);
  trgs.push_back(320);
  trgs.push_back(400);
  trgs.push_back(450);

  for (int itrg = 0; itrg != trgs.size(); ++itrg) {

    fitvals fv = drawQGLbin(0, 1.3, trgs[itrg]);

    ggdata->SetPoint(itrg, fv.pt, fv.fgdata);
    gqdata->SetPoint(itrg, fv.pt, fv.fqdata);
    ggreco->SetPoint(itrg, fv.pt, fv.fgreco);
    gqreco->SetPoint(itrg, fv.pt, fv.fqreco);
    gggen->SetPoint(itrg, fv.pt, fv.fggen);
    gqgen->SetPoint(itrg, fv.pt, fv.fqgen);
  }

  TH1D *h = new TH1D("h",";p_{T,jet};d#sigma_{gluon/quark}/d#sigma_{inc}",
		     3970,30,4000);
  h->SetMinimum(0);
  h->SetMaximum(0);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  //TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLogx();

  ggdata->SetMarkerStyle(kFullCircle);
  ggdata->SetMarkerColor(kBlue);
  ggdata->Draw("SAMEP");
  gqdata->SetMarkerStyle(kFullSquare);
  gqdata->SetMarkerColor(kRed);
  gqdata->Draw("SAMEP");

  ggreco->SetMarkerStyle(kOpenCircle);
  ggreco->Draw("SAMEP");
  gqreco->SetMarkerStyle(kOpenSquare);
  gqreco->Draw("SAMEP");

  gggen->SetMarkerStyle(kOpenCircle);
  gggen->SetMarkerColor(kBlue);
  gggen->Draw("SAMEP");
  gqgen->SetMarkerStyle(kOpenSquare);
  gqgen->SetMarkerColor(kRed);
  gqgen->Draw("SAMEP");

}
