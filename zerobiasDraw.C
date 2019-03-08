#include "tdrstyle_mod15.C"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TProfile.h"

#include "TLatex.h"

void zerobiasDraw() {

  setTDRStyle();

  //TFile *fin = new TFile("rootfiles/zerobias.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobiasL1Data.root","READ");
  ///TFile *fin = new TFile("rootfiles/zerobiasL1MC.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobiasL1RC.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobiasL1RCMC.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobiasL1DataGH.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobiasL1MCMC.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobias-data15.root","READ");
  //
  //TFile *fin = new TFile("rootfiles/zerobias-mc-v15.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobias-data-v16m.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobias-data-v17.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobias-data-v17c.root","READ");
  //TFile *fin = new TFile("rootfiles/zerobias-data-v18.root","READ");
  TFile *fin = new TFile("rootfiles/zerobias-data-v18v14fix.root","READ");
  assert(fin && !fin->IsZombie());

  TH1D *hraw = (TH1D*)fin->Get("hraw"); assert(hraw);
  TH1D *hl1 = (TH1D*)fin->Get("hl1"); assert(hl1);
  TH1D *hjec = (TH1D*)fin->Get("hjec"); assert(hjec);
  TH1D *hres = (TH1D*)fin->Get("hres"); assert(hres);

  TProfile *prho = (TProfile*)fin->Get("prho"); assert(prho);

  TH1D *hmu00up = (TH1D*)fin->Get("hmu00up"); assert(hmu00up);
  TH1D *hmu15dw = (TH1D*)fin->Get("hmu15dw"); assert(hmu15dw);
  TH1D *hmu1520 = (TH1D*)fin->Get("hmu1520"); assert(hmu1520);
  TH1D *hmu2025 = (TH1D*)fin->Get("hmu2025"); assert(hmu2025);
  TH1D *hmu2530 = (TH1D*)fin->Get("hmu2530"); assert(hmu2530);
  TH1D *hmu30up = (TH1D*)fin->Get("hmu30up"); assert(hmu30up);


  //TH1D *h = new TH1D("h",";#eta;N_{jet} / N_{events} #times 1000",104,-5.2,5.2);
  TH1D *h = new TH1D("h",";#eta;N_{jet} / N_{events}",104,-5.2,5.2);
  //lumi_13TeV = "ZeroBias 2016G Legacy";
  lumi_13TeV = "ZeroBias 2016GH Legacy";
  //lumi_13TeV = "Neutrino Gun";

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  //h->SetMaximum(2.1e-3*1000.);
  //h->SetMaximum(2.1e-3); // 30?
  h->SetMaximum(3e-2); // 
  //h->SetMaximum(1e-13); // PATCH
  
  //hraw->Scale(1000.);
  hraw->SetLineColor(kBlack);
  hraw->SetLineStyle(kDotted);
  hraw->Draw("SAME");

  //hl1->Scale(1000.);
  hl1->SetLineColor(kRed);
  hl1->SetLineStyle(kDashed);
  hl1->Draw("SAME");

  //hjec->Scale(1000.);
  hjec->SetLineColor(kBlue);
  hjec->SetLineStyle(kSolid);//kDashDotted);
  hjec->Draw("SAME");

  //hres->Scale(1000.);
  hres->SetLineColor(kGreen+2);
  hres->SetLineStyle(kSolid);
  hres->SetLineWidth(2);
  hres->Draw("SAME");

  TLegend *leg = tdrLeg(0.40,0.72,0.60,0.90);
  leg->AddEntry(hres,"L1L2L3+Res","L");
  leg->AddEntry(hjec,"L1L2L3","L");
  leg->AddEntry(hl1,"L1","L");
  leg->AddEntry(hraw,"No JEC","L");

  //double mu = 20.93;
  double mu = prho->GetMean();
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  //tex->DrawLatex(0.33,0.22,"AK4PFchs p_{T}>30 GeV");
  tex->DrawLatex(0.33,0.22,"AK4PFchs p_{T}>15 GeV");
  tex->DrawLatex(0.33,0.17,Form("#sigma_{jet}/#sigma_{MB} (#times#mu)"
				" = %1.2g (%1.2g)",
				hres->Integral() / mu,// * 0.001,
				hres->Integral()));// * 0.001));
  // https://arxiv.org/pdf/1106.0208.pdf Fig.1
  // d^2sigma/dpT/dy = 2*10^7 pb/GeV at sqrt(s)=7 TeV, pT=18 GeV, |y|<0.5
  // sigma_{MB} ~ 69 mb = 6.9*10^10 pb
  // => sigma_{jet} / sigma_{MB} = 3*10^-4
  // corresponds to pT~33 GeV at sqrt(s)=13 TeV

  //TF1 *f1 = new TF1("f1","[0]*pow(cosh(x),[1])",-5.2,5.2);
  // for pT^-5 spectrum integrated from x0 to infinity, <pT>=4/3*x0
  // => use 4/3*x0 in place of x0=30 GeV for kinematic cutoff term
  TF1 *f1 = new TF1("f1","[0]*pow(cosh(x),[1])"
		    "*pow(1-2*4/3*30*cosh(x)/13000.,[2])",-5.2,5.2);
  f1->SetParameters(hjec->GetBinContent(hjec->FindBin(0)), 0., 14);
  f1->SetLineColor(kBlue);
  f1->FixParameter(1,0);
  hjec->Fit(f1,"RNW ");
  f1->DrawClone("SAME");

  hl1->Fit(f1,"RWN");
  f1->SetLineStyle(kDashed);
  f1->SetLineColor(kRed);
  //f1->DrawClone("SAME");

  hres->Fit(f1,"RWN");
  f1->SetLineStyle(kSolid);
  f1->SetLineColor(kGreen+2);
  f1->DrawClone("SAME");
  

  
  //c1->SaveAs("pdf/zerobiasDraw.pdf");
  //
  //c1->SaveAs("pdf/zerobiasDraw-mc-v15.pdf");
  //c1->SaveAs("pdf/zerobiasDraw-data-v16m.pdf");
  //c1->SaveAs("pdf/zerobiasDraw-data-v17.pdf");
  //c1->SaveAs("pdf/zerobiasDraw-data-v17c.pdf");
  c1->SaveAs("pdf/zerobiasDraw-data-v18.pdf");


  TH1D *h2up = new TH1D("h2up",";#eta;N_{jet} / (N_{events}#times#mu)"
			//" #times 10000",104,-5.2,5.2);
			" #times 2000",104,-5.2,5.2);
  h2up->SetMaximum(3.5);//5);//1.5);
  TH1D *h2dw = new TH1D("h2dw",";#eta;N_{jet,#mu} / N_{jet}",
			104,-5.2,5.2);
  h2dw->SetMaximum(2.0);//1.3);
  h2dw->SetMinimum(0.5);//0.7);

  TCanvas *c2 = tdrDiCanvas("c2",h2up,h2dw,4,11);//,kSquare);
  
  c2->cd(1);

  hmu00up->Scale(2000);//00.);
  hmu00up->SetLineColor(kGreen+2);
  hmu00up->SetLineStyle(kSolid);
  hmu00up->SetLineWidth(2);
  hmu00up->Draw("SAME");  

  hmu15dw->Scale(2000);//00.);
  hmu15dw->SetLineColor(kBlue);
  hmu15dw->SetLineStyle(kSolid);
  ///hmu15dw->SetLineWidth(2);
  hmu15dw->Draw("SAME");  

  hmu1520->Scale(2000);//00.);
  hmu1520->SetLineColor(kMagenta+2);
  hmu1520->SetLineStyle(kDotted);
  ///hmu1520->SetLineWidth(2);
  hmu1520->Draw("SAME");  

  hmu2025->Scale(2000);//00.);
  hmu2025->SetLineColor(kYellow+2);
  hmu2025->SetLineStyle(kDashDotted);
  ///hmu2025->SetLineWidth(2);
  hmu2025->Draw("SAME");  

  hmu2530->Scale(2000);//00.);
  hmu2530->SetLineColor(kOrange+1);
  hmu2530->SetLineStyle(kDotted);
  ///hmu2530->SetLineWidth(2);
  hmu2530->Draw("SAME");  

  hmu30up->Scale(2000);//00.);
  hmu30up->SetLineColor(kRed);
  hmu30up->SetLineStyle(kSolid);
  ///hmu30up->SetLineWidth(2);
  hmu30up->Draw("SAME");  

  TLegend *leg2 = tdrLeg(0.40,0.60,0.60,0.90);
  leg2->AddEntry(hmu15dw,"#mu < 15","L");
  leg2->AddEntry(hmu1520,"15 #geq #mu < 20","L");
  leg2->AddEntry(hmu2025,"20 #geq #mu < 25","L");
  //leg2->AddEntry(hmu00up,"All #mu (#approx 21)","L");
  leg2->AddEntry(hmu00up,Form("All #mu (#approx %1.1f)",mu),"L");
  leg2->AddEntry(hmu2530,"25 #geq #mu < 30","L");
  leg2->AddEntry(hmu30up,"30 #geq #mu","L");

  //tex->DrawLatex(0.37,0.10,"AK4PFchs p_{T} > 30 GeV");
  tex->DrawLatex(0.37,0.10,"AK4PFchs p_{T} > 15 GeV");
  tex->DrawLatex(0.37,0.05,Form("#sigma_{jet}/#sigma_{MB} (#times#mu)"
				" = %1.2g (%1.2g)",
				hmu00up->Integral() * mu / 2000.,
				hmu00up->Integral() / 2000.));

  c2->cd(2);

  TH1D *hmu00up2 = (TH1D*)hmu00up->Clone("hmu00up2");
  hmu00up2->Divide(hmu00up);
  hmu00up2->Draw("SAME");

  TH1D *hmu15dw2 = (TH1D*)hmu15dw->Clone("hmu15dw2");
  hmu15dw2->Divide(hmu00up);
  hmu15dw2->Draw("SAME");

  TH1D *hmu30up2 = (TH1D*)hmu30up->Clone("hmu30up2");
  hmu30up2->Divide(hmu00up);
  hmu30up2->Draw("SAME");

  //c2->SaveAs("pdf/zerobiasDraw.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1RC.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1Data.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1MC.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1RCMC.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1DataGH.pdf");
  //c2->SaveAs("pdf/zerobiasDraw_L1MCMC.pdf");
  //
  //c2->SaveAs("pdf/zerobiasDrawTwoPad-mc-v15.pdf");
  //c2->SaveAs("pdf/zerobiasDrawTwoPad-data-v16m.pdf");
  //c2->SaveAs("pdf/zerobiasDrawTwoPad-data-v17.pdf");
  //c2->SaveAs("pdf/zerobiasDrawTwoPad-data-v17c.pdf");
  //c2->SaveAs("pdf/zerobiasDrawTwoPad-data-v18.pdf");
  c2->SaveAs("pdf/zerobiasDrawTwoPad-data-v18v14fix.pdf");
}
