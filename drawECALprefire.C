// Original code copied here from drawJetPt.C for easier manipulation and extensio
#include "tdrstyle_mod18.C"
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF2.h"
#include "TGraphErrors.h"

#include <string>
using namespace std;

void drawECALprefire(double eta = 2.0) {

  setTDRStyle();

  string sy = Form("%02d-%02d",int(10*eta+0.5),int(10*(eta+0.5)+0.5));
  const char *cy2 = sy.c_str();

  // Original ECAL prefire efficiency file:
  // fine eta bins, only some triggers and only AllEras
  TFile *fe = new TFile("rootfiles/ECALprefireEff.root");
  assert(fe && !fe->IsZombie());
  TH2D *heff = (TH2D*)fe->Get("h_JetHT_AllEras_inclusive"); assert(heff);

  // New ECAL prefire efficiency file:
  // wide eta bins, all triggers, and BCD, EF, GH separately
  TFile *fe2 = new TFile("rootfiles/ECALprefireEffRebinned.root");
  assert(fe2 && !fe2->IsZombie());
  TH2D *heff2bcd = (TH2D*)fe2->Get("h_JetHT_Run2016BCD_spike_finor_pteta");
  assert(heff2bcd);
  TH2D *heff2ef = (TH2D*)fe2->Get("h_JetHT_Run2016EF_spike_finor_pteta");
  assert(heff2ef);
  TH2D *heff2gh = (TH2D*)fe2->Get("h_JetHT_Run2016GH_spike_finor_pteta");
  assert(heff2gh);
  // Calculate combined efficiency using _unprescaled luminosities_
  // => fix later using per-trigger weights for each pT bin
  TH2D *heff2 = (TH2D*)heff2bcd->Clone("heff2");
  heff2->Scale(12.9/36.5);
  heff2->Add(heff2ef,6.8/36.5);
  heff2->Add(heff2gh,16.8/36.5);


  // Plot ECAL prefire inefficiency
  TH1D *h0 = new TH1D("h0",";Jet p_{T} (GeV);"
		      "ECAL prefire fraction",3445,55,3500);
  h0->SetMinimum(0.);
  h0->SetMaximum(1.);
  h0->GetXaxis()->SetNoExponent();
  h0->GetXaxis()->SetMoreLogLabels();

  TH1D *h02 = new TH1D("h02",";Jet p_{T} (GeV);"
		       "Ratio to fit",3445,55,3500);
  h02->SetMinimum(0.);
  h02->SetMaximum(2.0);
  h02->GetXaxis()->SetNoExponent();
  h02->GetXaxis()->SetMoreLogLabels();

  lumi_13TeV = "Run2016, 35 fb^{-1}";
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);

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

    //heff2->Scale(12.9/36.5);
    //heff2->Add(heff2ef,6.8/36.5);
    //heff2->Add(heff2gh,16.8/36.5);
    double effrb = heffetaRB->GetBinContent(i);
    double nrb = f2->Integral(pt1, pt2, eta, eta+0.5)*36500*0.002;
    double effrb_e = (nrb != 0 ? effrb*(1-effrb)/sqrt(nrb) : 0);
    //
    double effbcd = heffetaBCD->GetBinContent(i);
    double nbcd = f2->Integral(pt1, pt2, eta, eta+0.5)*12900*0.002;
    double effbcd_e = (nbcd != 0 ? effbcd*(1-effbcd)/sqrt(nbcd) : 0);
    //
    double effef = heffetaEF->GetBinContent(i);
    double nef = f2->Integral(pt1, pt2, eta, eta+0.5)*6800*0.002;
    double effef_e = (nef != 0 ? effef*(1-effef)/sqrt(nef) : 0);
    //
    double effgh = heffetaGH->GetBinContent(i);
    double ngh = f2->Integral(pt1, pt2, eta, eta+0.5)*16800*0.002;
    double effgh_e = (ngh != 0 ? effgh*(1-effgh)/sqrt(ngh) : 0);

    if (eff!=0 && eff_e!=0) {
      int n = geffeta->GetN();
      double mine = (eta<2.5 ? 0.003 : 0.01);
      geffeta->SetPoint(n, 0.5*(pt1+pt2), eff);
      geffeta->SetPointError(n, 0.5*(pt2-pt1)/sqrt(n1+n2),
			     sqrt(mine*mine+eff_e*eff_e));
      //
      geffetaRB->SetPoint(n, 0.5*(pt1+pt2), heffetaRB->GetBinContent(i));
      geffetaRB->SetPointError(n, 0.5*(pt2-pt1)/sqrt(nrb),
			       sqrt(mine*mine+effrb_e*effrb_e));
      //
      geffetaBCD->SetPoint(n, 0.5*(pt1+pt2), heffetaBCD->GetBinContent(i));
      geffetaBCD->SetPointError(n, 0.5*(pt2-pt1)/sqrt(nbcd),
				sqrt(mine*mine+effbcd_e*effbcd_e));
      //
      geffetaEF->SetPoint(n, 0.5*(pt1+pt2), heffetaEF->GetBinContent(i));
      geffetaEF->SetPointError(n, 0.5*(pt2-pt1)/sqrt(nef),
			       sqrt(mine*mine+effef_e*effef_e));
      //
      geffetaGH->SetPoint(n, 0.5*(pt1+pt2), heffetaGH->GetBinContent(i));
      geffetaGH->SetPointError(n, 0.5*(pt2-pt1)/sqrt(ngh),
			       sqrt(mine*mine+effgh_e*effgh_e));
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

  // Use combined fit as starting point for others, for better convergence
  feffBCD->SetParameters(feff->GetParameter(0),feff->GetParameter(1),feff->GetParameter(2));
  feffEF->SetParameters(feff->GetParameter(0),feff->GetParameter(1),feff->GetParameter(2));
  feffGH->SetParameters(feff->GetParameter(0),feff->GetParameter(1),feff->GetParameter(2));
  // Reduce low stat EF freedom a bit to better match the shape of others
  // And while at it, can could the others as well
  //feffBCD->FixParameter(2,feff->GetParameter(2));
  //if (eta<2.5) 
  feffEF->FixParameter(2,feff->GetParameter(2));
  //feffGH->FixParameter(2,feff->GetParameter(2));

  geffetaBCD->Fit(feffBCD,"QRN"); feffBCD->SetLineColor(kGreen+2);
  geffetaEF->Fit(feffEF,"QRN"); feffEF->SetLineColor(kOrange+2);
  geffetaGH->Fit(feffGH,"QRN"); feffGH->SetLineColor(kRed);
  
  //tdrDraw(heffeta,"z",kFullSquare,kBlack,kSolid,kBlack);
  //tdrDraw(geffeta,"Pz",kFullSquare,kBlack,kSolid,kBlack);
  tdrDraw(heffetaRB,"HIST Pz",kFullDiamond,kBlue,kSolid,kBlue);
  tdrDraw(geffetaRB,"Pz",kFullDiamond,kBlue,kSolid,kBlue);
  //tdrDraw(heffeta12,"HIST Pz",kOpenSquare,kBlack,kSolid,kBlack);
  //tdrDraw(heffeta1,"HIST PLz",kOpenTriangleUp,kBlack,kSolid,kBlack,kNone);
  //tdrDraw(heffeta2,"HIST PLz",kOpenTriangleDown,kBlack,kSolid,kBlack,kNone);
  feff->Draw("SAME");
  feffBCD->Draw("SAME");
  feffEF->Draw("SAME");
  feffGH->Draw("SAME");

  tdrDraw(heffetaBCD,"HIST Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  tdrDraw(geffetaBCD,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  tdrDraw(heffetaEF,"HIST Pz",kOpenDiamond,kOrange+2,kSolid,kOrange+2);
  tdrDraw(geffetaEF,"Pz",kOpenDiamond,kOrange+2,kSolid,kOrange+2);
  tdrDraw(heffetaGH,"HIST Pz",kOpenDiamond,kRed,kSolid,kRed);
  tdrDraw(geffetaGH,"Pz",kOpenDiamond,kRed,kSolid,kRed);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kBlue);
  //double d = (eta<2.5 ? 0.18 : 0);
  double dy = (eta<2.5 ? 0.54 : 0.26);
  double dx = (eta<2.5 ? 0.20 : 0.45);
  tex->DrawLatex(dx,0.10+dy,"#frac{C}{2}#times#left(1 + "
		 "erf#left(#frac{p_{T}-#mu}{S#sqrt{p_{T}}}#right)#right)");
  tex->SetTextSize(0.030); double dd = 0.035;//0.05;
  for (int i = 0; i != 4; ++i) {
    TF1 *f(0); int color(kBlack); double dx2(dx); double dy2(dy);
    if (i==0) { f = feff; color = kBlue; dx2 += (eta<2.5 ? 0 : 0.22); if(eta>=2.5) dy2=0.62;}
    if (i==1) { f = feffGH; color = kRed; dx2 += (eta<2.5 ? 0.22 : 0);}
    if (i==2) { f = feffEF; color = kOrange+2; dx2 += (eta<2.5 ? 0.40 : 0.18); }
    if (i==3) { f = feffBCD; color = kGreen+2; dx2 += (eta<2.5 ? 0.58 : 0.36); }
    tex->SetTextColor(color);
    tex->DrawLatex(dx2,dy2-0*dd,Form("%s%1.1f/%d",
				     i==0 ? "#chi^{2}/NDF=" : "",
				     f->GetChisquare(), f->GetNDF()));
    tex->DrawLatex(dx2,dy2-1*dd,Form("%s%1.3f#pm%1.3f", i==0 ? "C=" : "",
				     f->GetParameter(0), f->GetParError(0)));
    tex->DrawLatex(dx2,dy2-2*dd,Form("%s%1.1f#pm%1.1f", i==0 ? "#mu=" : "",
				     f->GetParameter(1), f->GetParError(1)));
    tex->DrawLatex(dx2,dy2-3*dd,Form("%s%1.1f#pm%1.1f", i==0 ? "S=" : "",
				     f->GetParameter(2), f->GetParError(2)));
  }
  

  tex->SetTextSize(0.045);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.19,0.75,Form("%1.1f#leq|y|<%1.1f",eta,eta+0.5));
  
  //TLegend *leg0 = tdrLeg(0.62,0.65,0.82,0.90);
  //leg0->AddEntry(heffeta2,Form("%1.2f#leq|y|<%1.2f",eta,eta+0.25),"PL");
  //leg0->AddEntry(heffeta1,Form("%1.2f#leq|y|<%1.2f",eta+0.25,eta+0.5),"PL");
  //leg0->AddEntry(heffeta12,"direct avg.","P");
  //leg0->AddEntry(heffeta,"weighted avg.","P");
  //leg0->AddEntry(heffetaRB,"rebinned","P");

  TLegend *leg0 = tdrLeg(0.62,0.70,0.82,0.90);
  leg0->AddEntry(heffetaRB,"Weighted avg.","PL");
  leg0->AddEntry(heffetaGH,"GH","PL");
  leg0->AddEntry(heffetaEF,"EF","PL");
  leg0->AddEntry(heffetaBCD,"BCD","PL");

  c0->SaveAs(Form("pdf/drawECALprefire_%s.pdf",cy2));

}
