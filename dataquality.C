#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "tdrstyle_mod18.C"

// Draw 2D plot of jet rates in (eta,phi) to spot issues
void dataquality() {
  const char *eras = "GH";
  const double maxgoodsig = 2.;
  const double minsig = 3.;
  const double maxsig = 8.5;
  const double minsumsig = 9.;

  TDirectory *curdir = gDirectory;
  setTDRStyle();
  gStyle->SetOptTitle();
  gStyle->SetPalette(1);

  TFile *fd = new TFile("./output-DATA-1.root","READ");
  assert(fd && !fd->IsZombie());

  TFile *fm = new TFile("./output-MC-1.root","READ");
  assert(fm && !fm->IsZombie());

  string triggers[] = {"jt40","jt60","jt80","jt140","jt200","jt260","jt320","jt400","jt450"};
  const int ntrg = sizeof(triggers)/sizeof(triggers[0]);
  cout << "Using " << ntrg << " triggers!" << endl;

  double etabins[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;
  cout << "Using " << neta << " etabins!" << endl;

  vector<double> etalims = {4.2,4.0,3.8,3.4,3.0,2.6};
  vector<int> etatrgs = {1,2,3,4,5,9};
  vector<double> etaws = {9.0,4.5,3.0,2.25,1.8,1.0}; // ngoodtrig/nalltrigs

  TH2D *h2hotr(0);
  TH2D *h2hotm(0);
  TH2D *h2hotHCAL(0);

  // Create map of known hot ECAL regions from Robert Schoefbeck:
  // https://github.com/schoef/JetMET/blob/master/JEC/python/L2res/jet_cleaning.py#L2-L8
  const int nhot = 7;
  double thr[nhot][4] = {
    { -2.650, -2.500,  -1.35,    -1.05 },
    { -2.964, -2.650,  -1.10,    -0.80 },
    { -2.964, -2.650,  -0.25,     0.10 },
    { -2.964, -2.650,  -3.14159, -2.80 },
    { -2.964, -2.650,   2.9,      3.14159 },
    {  2.650,  2.964,  -2.00,    -1.60 },
    {  2.650,  3.139,   0.00,     0.25 }
  };
  const int nhot2 = 10;
  double dx = TMath::TwoPi()/72.;
  double thr2[nhot2][4] = {
  // Mikko's extra regions first
    { 2.60,  2.60+dx,   -2.70-dx, -2.70+dx},
    { 1.40,  1.40+2*dx, -0.80-dx, -0.80+dx},
    {-3.20, -3.20+dx,    0.20,     0.20+2*dx},
    { 2.60,  2.60+dx,    0.35,     0.35+3*dx},
    { 2.40,  2.40+2*dx, -0.10,    -0.10+2*dx },
  // Robert's sub-regions next
    { 2.80,     2.80+2*dx, 0.00,  0+2*dx},
    { 2.70,     2.70+dx,  -1.80, -1.8+2*dx},
    {-2.80-dx, -2.80+dx,   3.00,  3.0+2*dx},
    {-2.70,    -2.70+dx,  -1.20, -1.2+2*dx},
    {-2.80,    -2.80+dx,  -3.15, -3.15+2*dx}
  };
  const int nhotHCAL = 3;
  double thrHCAL[nhotHCAL][4] = {
    {-2.250, -1.930, 2.200, 2.500},
    {-3.489, -3.139, 2.237, 2.475},
    {-3.600, -3.139, 2.237, 2.475}
  };

  TH2D *h2s[ntrg], *h2as[ntrg], *h2bs[ntrg];
  TH2D *h2hots[ntrg], *h2colds[ntrg];
  TH1D *hrmss[ntrg], *hrms2s[ntrg];

  for (int dtmc = 0; dtmc <= 1; ++dtmc) {
    TFile *f = (dtmc==0) ? fd : fm;
    assert(f->cd("Standard"));
    TDirectory *din = gDirectory;
    string nametag = (dtmc==0) ? "data" : "mc";
    string roottag = (dtmc==0) ? "" : "mc";
    for (int itrg = 0; itrg < ntrg; ++itrg) {
      string strg = triggers[itrg];
      const char *ctrg = triggers[itrg].c_str();

      TH2D *h2 = 0;
      for (int ieta = 0; ieta < neta; ++ieta) {
        double etamin = etabins[ieta];
        double etamax = etabins[ieta+1];

        assert(din->cd(Form("Eta_%1.1f-%1.1f",etamin,etamax)));
        assert(gDirectory->cd(ctrg));
        TDirectory *d = gDirectory;

        TH2D *h = (TH2D*)d->Get("hetaphi");
        assert(h);
        if (!h2) h2 = (TH2D*)h->Clone("h2");
        else h2->Add(h);
      } // for ieta

      TH2D *h2hot = (TH2D*)h2->Clone("h2hot");
      TH2D *h2cold = (TH2D*)h2->Clone("h2cold");
      if (itrg==0) {
        h2hotr = (TH2D*)h2->Clone("h2hotr");
        h2hotm = (TH2D*)h2->Clone("h2hotm");
        h2hotHCAL = (TH2D*)h2->Clone("h2hotHCAL");
        for (int idxeta = 1; idxeta != h2hot->GetNbinsX()+1; ++idxeta) {
          for (int idxphi = 1; idxphi != h2hot->GetNbinsY()+1; ++idxphi) {
            double eta = h2hot->GetXaxis()->GetBinCenter(idxeta);
            double phi = h2hot->GetYaxis()->GetBinCenter(idxphi);

            // Reset hot region map (fill later from h2a)
            h2hot->SetBinContent(idxeta, idxphi, 0);
            h2hot->SetBinError(idxeta, idxphi, 0);
            h2cold->SetBinContent(idxeta, idxphi, 0);
            h2cold->SetBinError(idxeta, idxphi, 0);

            // Produce map for Robert's hot spots
            h2hotr->SetBinContent(idxeta, idxphi, -10);
            h2hotr->SetBinError(idxeta, idxphi, 0);
            for (int k = 0; k != nhot; ++k) {
              if (eta >= thr[k][0] and eta <= thr[k][1] and phi >= thr[k][2] and phi <= thr[k][3])
                h2hotr->SetBinContent(idxeta, idxphi, 10);
            } // for k

            // Produce map for Mikko's manual hot spots
            h2hotm->SetBinContent(idxeta, idxphi, -10);
            h2hotm->SetBinError(idxeta, idxphi, 0);
            for (int k = 0; k != nhot2; ++k) {
              if (eta >= thr2[k][0] and eta <= thr2[k][1] and phi >= thr2[k][2] and phi <= thr2[k][3])
                h2hotm->SetBinContent(idxeta, idxphi, 10);
            } // for k

            // Produce map for Mikko's manual hot spots
            h2hotHCAL->SetBinContent(idxeta, idxphi, -10);
            h2hotHCAL->SetBinError(idxeta, idxphi, 0);
            for (int k = 0; k != nhotHCAL; ++k) {
              if (eta >= thrHCAL[k][0] and eta <= thrHCAL[k][1] and phi >= thrHCAL[k][2] and phi <= thrHCAL[k][3])
                h2hotHCAL->SetBinContent(idxeta, idxphi, 10);
            } // for k
          } // for j
        } // for i
      }
      vector<double> xbins;
      for (unsigned idxeta = 1u; idxeta < h2->GetNbinsX()+1; ++idxeta)
        xbins.push_back(h2->GetXaxis()->GetBinLowEdge(idxeta));
      xbins.push_back(h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX()));

      TH2D *h2a = (TH2D*)h2->Clone("h2a"); // statistical fluctuation
      TH2D *h2b = (TH2D*)h2->Clone("h2b"); // relative fluctuation
      TH2D *h2dummy = (TH2D*)h2->Clone("h2dummy"); // Helper
      TH1D *hrms = new TH1D("hrms","Relative RMS of #eta_{jet} strip;#eta_{jet}Relative RMS",xbins.size()-1,&xbins[0]);
      TH1D *hrms2 = new TH1D("hrms2","Relative RMS of #eta_{jet} strip;#eta_{jet};Expected relative RMS",xbins.size()-1,&xbins[0]);
      // For each eta value we average a "pedestal" over phi and calculate the relative fluctuation.
      for (int idxeta = 1; idxeta < h2->GetNbinsX()+1; ++idxeta) {
        // 1. Calculate pedestal [ped] (mean value of populated bins)
        int nbins(0), nbinsgood(0);
        double rms2good(0), rms2(0), sumgood(0), sum(0);
        for (int idxphi = 1; idxphi < h2->GetNbinsY()+1; ++idxphi) {
          if (h2->GetBinContent(idxeta,idxphi)>0) {
            sum += h2->GetBinContent(idxeta,idxphi);
            ++nbins;
          }
        }
        double ped = sum/nbins;

        // 2. Subtract pedestal, calculate raw rms [rms] (prerequisite: ped)
        for (int idxphi = 1; idxphi < h2->GetNbinsY()+1; ++idxphi) {
          double diff = 0.0;
          if (h2->GetBinContent(idxeta, idxphi)>0) {
            diff = h2->GetBinContent(idxeta,idxphi)-ped;
            rms2 += pow(diff, 2);
          }
          h2dummy->SetBinContent(idxeta, idxphi, fabs(diff)/maxgoodsig);
        } // for j
        double rms = sqrt(rms2/nbins);

        // 3. Calculate cleaned pedestal with only [-maxsig,+maxsig] range [pedgood] (prerequisite: h2dummy, rms)
        for (int idxphi = 1; idxphi < h2->GetNbinsY()+1; ++idxphi) {
          if (h2->GetBinContent(idxeta,idxphi)>0 and h2dummy->GetBinContent(idxeta,idxphi) < rms) {
            sumgood += h2->GetBinContent(idxeta,idxphi);
            ++nbinsgood;
          }
        }
        double pedgood = sumgood/nbinsgood;

        // 4. Substract good pedastal, calculate good rms [rmsgood] (prerequisite: pedgood, h2dummy, rms)
        for (int idxphi = 1; idxphi < h2->GetNbinsY()+1; ++idxphi) {
          if (h2->GetBinContent(idxeta,idxphi)>0 and h2dummy->GetBinContent(idxeta,idxphi) < rms)
            rms2good += pow(h2->GetBinContent(idxeta,idxphi)-pedgood, 2);
        }
        double rmsgood = sqrt(rms2good/nbinsgood);

        // One could repeat the steps 1-4 until required convergence
        hrms->SetBinContent(idxeta, rmsgood/pedgood);
        hrms2->SetBinContent(idxeta, 1./sqrt(pedgood));
        // For MC, replace sqrt(ped) with rms, because not using event counts
        for (int idxphi = 1; idxphi < h2->GetNbinsY()+1; ++idxphi) {
          if (h2->GetBinContent(idxeta,idxphi)>0) {
            h2a->SetBinContent(idxeta, idxphi, (h2->GetBinContent(idxeta,idxphi)-pedgood)/rmsgood);
            h2b->SetBinContent(idxeta, idxphi, (h2->GetBinContent(idxeta,idxphi)-pedgood)/pedgood);
          }
        } // for idxphi
      } // for idxeta

      // Keep track of found hot and cold regions
      for (int idxeta = 1; idxeta < h2hot->GetNbinsX()+1; ++idxeta) {
        for (int idxphi = 1; idxphi < h2hot->GetNbinsY()+1; ++idxphi) {
          if (h2a->GetBinContent(idxeta,idxphi)>=minsig)
            h2hot->SetBinContent(idxeta,idxphi, min(maxsig, h2a->GetBinContent(idxeta,idxphi)));
          else
            h2hot->SetBinContent(idxeta,idxphi, 0);
          if (h2a->GetBinContent(idxeta,idxphi)<=-minsig)
            h2cold->SetBinContent(idxeta,idxphi, min(maxsig, -h2a->GetBinContent(idxeta,idxphi)));
          else
            h2cold->SetBinContent(idxeta,idxphi, 0);
        }
      }

      h2s[itrg] = h2;
      h2as[itrg] = h2a;
      h2bs[itrg] = h2b;
      hrmss[itrg] = hrms;
      hrms2s[itrg] = hrms2;
      h2hots[itrg] = h2hot;
      h2colds[itrg] = h2cold;
    } // itrg

    // Sum up hot region maps
    TH2D *h2cumul = (TH2D*)h2as[0]->Clone("h2cumul");
    TH2D *h2hot = (TH2D*)h2hots[0]->Clone("h2hot");
    TH2D *h2hot2 = (TH2D*)h2hots[0]->Clone("h2hot2");
    TH2D *h2cold = (TH2D*)h2colds[0]->Clone("h2cold");
    TH2D *h2cold2 = (TH2D*)h2colds[0]->Clone("h2cold2");
    for (int idxeta = 1; idxeta < h2hot->GetNbinsX()+1; ++idxeta) {
      double abseta = h2hot->GetXaxis()->GetBinCenter(idxeta);
      int breaker = -1;
      for (auto ie = 0u; ie < etalims.size(); ++ie) {
        if (abseta>etalims[ie]) {
          breaker = ie;
          break;
        }
      }
      int maxtrg = 0;
      double scale = 1.0;
      if (breaker<0) maxtrg = ntrg;
      else if (breaker>0) {
        maxtrg = etatrgs[breaker-1];
        scale = etaws[breaker-1];
      }

      for (int idxphi = 1; idxphi < h2hot->GetNbinsY()+1; ++idxphi) {
        // Hot regions
        h2hot->SetBinContent(idxeta, idxphi, 0);
        h2hot2->SetBinContent(idxeta, idxphi, 0);
        // Cold regions
        h2cold->SetBinContent(idxeta, idxphi, 0);
        h2cold2->SetBinContent(idxeta, idxphi, 0);
        for (int itrg = 0; itrg < maxtrg; ++itrg) {
          h2hot->SetBinContent(idxeta,idxphi, h2hot->GetBinContent(idxeta,idxphi) + scale*h2hots[itrg]->GetBinContent(idxeta,idxphi));
          h2cold->SetBinContent(idxeta,idxphi, h2cold->GetBinContent(idxeta,idxphi) + scale*h2colds[itrg]->GetBinContent(idxeta,idxphi));
          h2cumul->SetBinContent(idxeta,idxphi, h2cumul->GetBinContent(idxeta,idxphi) + scale*h2as[itrg]->GetBinContent(idxeta,idxphi));
        } // for itrg

        // Select regions with excesses in at least two triggers
        if (h2hot->GetBinContent(idxeta,idxphi)>=minsumsig)
          h2hot2->SetBinContent(idxeta,idxphi,10);
        else
          h2hot2->SetBinContent(idxeta,idxphi,-10);

        // Select regions with deficits in at least two triggers
        if (h2cold->GetBinContent(idxeta,idxphi)>=minsumsig)
          h2cold2->SetBinContent(idxeta,idxphi,10);
        else
          h2cold2->SetBinContent(idxeta,idxphi,-10);
      } // for idxphi
    } // for idxeta

    h2hotr->SetFillStyle(0);
    h2hotr->SetLineColor(kBlack);
    h2hotr->SetFillColor(kBlack);
    h2hotr->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2hotr->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
    h2hotr->GetZaxis()->SetRangeUser(0,10);
    h2hotm->SetFillStyle(0);
    h2hotm->SetLineColor(kRed);
    h2hotm->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2hotm->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
    h2hotm->GetZaxis()->SetRangeUser(0,10);
    h2hotHCAL->SetFillStyle(0);
    h2hotHCAL->SetLineColor(kMagenta+1);
    h2hotHCAL->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2hotHCAL->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
    h2hotHCAL->GetZaxis()->SetRangeUser(0,10);

    h2hotr->Scale(h2s[ntrg-1]->Integral());
    h2hotm->Scale(h2s[ntrg-1]->Integral());
    h2hotHCAL->Scale(h2s[ntrg-1]->Integral());
    for (int itrg = 0; itrg != ntrg; ++itrg) {
      const char *ctrg = triggers[itrg].c_str();

      TCanvas *c1 = new TCanvas("c1","c1",600,600);
      gPad->SetLeftMargin(0.10);
      gPad->SetRightMargin(0.15);
      gPad->SetTopMargin(0.10);
      gPad->SetBottomMargin(0.10);

      h2s[itrg]->SetTitle("Number of jets;#eta_{jet};#phi_{jet}");
      h2s[itrg]->GetXaxis()->SetRangeUser(-4.79,4.79);
      h2s[itrg]->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

      h2hotHCAL->SetFillStyle(0);
      h2s[itrg]->DrawClone("COLZ");
      h2hotr->DrawClone("SAMEBOX");
      h2hotHCAL->DrawClone("SAMEBOX");
      h2hotm->DrawClone("SAMEBOX");
      gPad->Update();
      c1->SaveAs(Form("pdf/%squality_njet_%s.pdf",nametag.c_str(),ctrg));

      TCanvas *c2a = new TCanvas("c2a","c2a",600,600);
      gPad->SetLeftMargin(0.10);
      gPad->SetRightMargin(0.15);
      gPad->SetTopMargin(0.10);
      gPad->SetBottomMargin(0.10);

      // For drawing, mark deficit <-8 as -8 (so blue box instead of white)
      h2as[itrg]->SetTitle(Form("Significance of excess/deficit (%s);#eta_{jet};#phi_{jet}",ctrg));

      h2as[itrg]->SetMinimum(-8);
      h2as[itrg]->SetMaximum(+8);
      h2as[itrg]->GetXaxis()->SetRangeUser(-4.79,4.79);
      h2as[itrg]->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

      h2as[itrg]->DrawClone("COLZ");
      h2hotr->DrawClone("SAMEBOX");
      h2hotHCAL->DrawClone("SAMEBOX");
      h2hotm->DrawClone("SAMEBOX");

      TCanvas *c2b = new TCanvas("c2b","c2b",600,600);
      gPad->SetLeftMargin(0.10);
      gPad->SetRightMargin(0.15);
      gPad->SetTopMargin(0.10);
      gPad->SetBottomMargin(0.10);

      h2bs[itrg]->SetTitle("Relative fluctuation;#eta_{jet};#phi_{jet}");
      h2bs[itrg]->SetMinimum(-1.0);
      h2bs[itrg]->SetMaximum(+1.0);
      h2bs[itrg]->GetXaxis()->SetRangeUser(-4.79,4.79);
      h2bs[itrg]->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

      h2bs[itrg]->Draw("COLZ");
      h2hotr->DrawClone("BOX SAME");
      h2hotHCAL->DrawClone("BOX SAME");
      h2hotm->DrawClone("BOX SAME");

      //TCanvas *c3 = new TCanvas("c3","c3",600,600);
      //gPad->SetTopMargin(0.10);
      ////hrmss[itrg]->SetMinimum(0.);
      ////hrmss[itrg]->SetMaximum(0.5);
      //hrmss[itrg]->SetLineColor(kRed);
      //hrmss[itrg]->DrawClone();
      //hrms2s[itrg]->SetLineColor(kBlue);
      //hrms2s[itrg]->DrawClone("SAME");

      c2a->SaveAs(Form("pdf/%squality_significance_%s.pdf",nametag.c_str(),ctrg));
      c2b->SaveAs(Form("pdf/%squality_relfluctuation_%s.pdf",nametag.c_str(),ctrg));
      //c3->SaveAs(Form("pdf/%squality_relRMS_%s.pdf",nametag.c_str(),ctrg));
    }
    TCanvas *c0 = new TCanvas("c0","c0",600,600);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);

    h2cumul->SetTitle("Cumulative h2as;#eta_{jet};#phi_{jet}");
    h2cumul->SetMinimum(-45);
    h2cumul->SetMaximum(+45);
    h2cumul->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2cumul->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

    h2cumul->Draw("COLZ");
    h2hotr->Draw("SAMEBOX");
    h2hotHCAL->Draw("SAMEBOX");
    h2hotm->Draw("SAMEBOX");

    TCanvas *c0hot = new TCanvas("c0hot","c0hot",600,600);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);

    h2hot->SetTitle("Cumulative h2hot;#eta_{jet};#phi_{jet}");
    h2hot->SetMinimum(-20);
    h2hot->SetMaximum(+20);
    h2hot->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2hot->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

    h2hot->DrawClone("COLZ");
    h2hotr->DrawClone("SAMEBOX");
    h2hotHCAL->DrawClone("SAMEBOX");
    h2hotm->DrawClone("SAMEBOX");

    TCanvas *c0cold = new TCanvas("c0cold","c0cold",600,600);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);

    h2cold->SetTitle("Cumulative h2cold;#eta_{jet};#phi_{jet}");
    h2cold->SetMinimum(-30);
    h2cold->SetMaximum(+30);
    h2cold->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2cold->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());

    h2cold->DrawClone("COLZ");
    h2hotr->DrawClone("SAMEBOX");
    h2hotHCAL->DrawClone("SAMEBOX");
    h2hotm->DrawClone("SAMEBOX");

    c0->SaveAs(Form("pdf/%squality_cumulation.pdf",nametag.c_str()));
    c0hot->SaveAs(Form("pdf/%squality_hots.pdf",nametag.c_str()));
    c0cold->SaveAs(Form("pdf/%squality_colds.pdf",nametag.c_str()));

    h2hot2->GetXaxis()->SetRangeUser(-4.79,4.79);
    h2hot2->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
    TFile *fouthot = new TFile(Form("rootfiles/hotjets%s-run%s.root",roottag.c_str(),eras),"RECREATE");
    h2hot->Write("h2hot");
    h2hot2->Write("h2hotfilter");
    //h2hotr->Write("h2hotrobert");
    //h2hotm->Write("h2hotmikko");
    //h2hotHCAL->Write("h2hotHCAL");
    fouthot->Close();
    TFile *foutcold = new TFile(Form("rootfiles/coldjets%s-run%s.root",roottag.c_str(),eras),"RECREATE");
    h2cold->Write("h2cold");
    h2cold2->Write("h2hole");
    foutcold->Close();
  }
}
