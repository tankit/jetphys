#include "TFile.h"
#include "TH1D.h"

#include "tdrstyle_mod18.C"

// Code for template fit
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

// Loop for template fit over all the pT bins
/*
void tempHist(TH1D *hg, TH1D *hq, double eta, bool isData = true) {

  assert(hg);
  assert(hq);

  for (int ipt = 1; ipt != hg->GetNbinsX()+1; ++ipt) {

    double pt = hg->GetBinCenter(ipt);
    // select correct trigger based on pt
    // retrieve hqgl (hm), hqgl_q (hq), hqgl_g (hg) for MC and hqgl (hd) for data
    // retrieve the overall cross section also (hpt)
    // do template fit f1d or f1m as in drawGluon
    if (isData) hd->Fit(ft,"QRN");
    else        hm->Fit(ft,"QRN");
    

    // set template fit values to histograms, multiply by pT spectrum
    double sigma = hpt->GetBinContent(ipt);
    double esigma = hpt->GetBinError(ipt);
    hg->SetBinContent(ipt, sigma*ft->GetParameter(0));
    hg->SetBinError(ipt, pow(sigma*ft->GetParError(0), esigma*ft->GetParameter(0),2));

    hq->SetBinContent(ipt, sigma*(1-ft->GetParameter(0)));
    hq->SetBinError(ipt, pow(sigma*ft->GetParError(0), esigma*(1-ft->GetParameter(0)),2));
  }

  return;
}
*/


// Run with 'root -l drawGluon.C+g'
void drawGluon() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fd = new TFile("output-DATA_RunGfull_part-1.root","READ");
  assert(fd && !fd->IsZombie());
  TFile *fm = new TFile("output-MC_30Octlongprobg-1.root","READ");
  assert(fm && !fm->IsZombie());

  string dir = "Standard/Eta_0.0-1.3/mc";
  //TH1D *hg = (TH1D*)f->Get(Form("%s/hgpt",dir.c_str()));
  TH1D *hgpt = (TH1D*)fm->Get(Form("%s/hgpt_g0tw",dir.c_str()));
  assert(hgpt);
  //TH1D *ha = (TH1D*)f->Get(Form("%s/hpt",dir.c_str()));
  TH1D *hapt = (TH1D*)fm->Get(Form("%s/hpt_g0tw",dir.c_str()));
  assert(hapt);

  curdir->cd();

  TH1D *hg = (TH1D*)hgpt->Clone("hg");
  hg->Divide(hapt);

  TH1D *hq = (TH1D*)hgpt->Clone("hq");
  hq->Divide(hgpt);
  hq->Add(hg,-1);

  TH1D *h = new TH1D("h",";p_{T,jet} (GeV);"
		     "d#sigma_{gluon/quark} / d#sigma_{inclusive}",
		     100,30,4000);
  h->SetMinimum(0.);
  h->SetMaximum(1.1);
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetMoreLogLabels();


  lumi_13TeV = "Ozlem's analysis, 36 fb^{-1}";
  //second plot
  //second plot

  // also tdrDiCanvas
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLogx();
  hg->GetXaxis()->SetRangeUser(30,4000.); // for ROOT plotting bug
  hq->GetXaxis()->SetRangeUser(30,4000.); // for ROOT plotting bug
  tdrDraw(hg,"P",kFullCircle,kBlue);
  tdrDraw(hq,"P",kOpenSquare,kRed);

  TLegend *leg = tdrLeg(0.60, 0.78, 0.80, 0.90);
  leg->AddEntry(hg,"Gluons","LP");
  leg->AddEntry(hq,"Quarks","LP");
  //c2 
  //c2
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();

  tex->DrawLatex(0.20,0.74,"0.0 < |y_{jet}| < 1.3");
  tex->DrawLatex(0.20,0.68,"anti-k_{T} R=0.4 PF+CHS");


  c1->RedrawAxis();

  c1->SaveAs("pdf/drawGluon_eta0.0_1.3.pdf");
}


