// Purpose: Normalize inclusive jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Co-author: hannu.siikonen@cern.ch
// Created: March 21, 2010
// Updated: April 4, 2017
// Archaic blocks denoted by wide commenting


#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TRegexp.h"

#include <iostream>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

#include "settings.h"

// A list of pt dependent stuff (affects how we process)
vector<string> hptlike = {"hpt",
                          "hpt_evt",
                          "hpt_jet",
                          "hpt_pre",
                          "hpt0",
                          "hpt1",
                          "hpt2",
                          "hpt3",
                          "hpt_jk1",
                          "hpt_jk2",
                          "hpt_jk3",
                          "hpt_jk4",
                          "hpt_jk5",
                          "hpt_jk6",
                          "hpt_jk7",
                          "hpt_jk8",
                          "hpt_jk9",
                          "hpt_jk10",
                          "hpt_noid",
                          "hpt_noevtid",
                          "hpt_nojetid",
                          "hselpt",
                          "hpt_r",
                          "hpt_g",
                          "hpt_gg",
                          "hpt_g0",
                          "hpt_g0tw",
                          "hdjmass",
                          "hdjmass0"};

void recurseFile(TDirectory *indir, TDirectory *outdir,
                 double etawid = 1., double etamid = 0., int lvl = 0);

// Use this to fix luminosity
double _lumiscale = 1.00;
std::map<std::string, double> triglumi;

void HistosNormalize() {
  TFile *fin = new TFile(Form("output-%s-1.root",_jp_type),"READ");
  assert(fin && !fin->IsZombie());

  TFile *fout = new TFile(Form("output-%s-2a.root",_jp_type),"RECREATE");
  assert(fout and !fout->IsZombie());

  if (_lumiscale!=1 and _jp_isdt)
    cout << "Attention! : Scaling luminosity to the new estimate"
         << " by multiplying with " << _lumiscale << endl;

  if (_jp_usetriglumi) { // Setting up lumis
    cout << "Reading trigger luminosity from settings.h" << endl;
    for (unsigned int i = 0; i < _jp_notrigs; ++i) {
      double lumi = _jp_triglumi[i]/1e6; // /ub to /pb
      cout << Form(" *%s: %1.3f /pb", _jp_triggers[i],lumi) << endl;
      triglumi[_jp_triggers[i]] = lumi;
    }
    if (_jp_ismc) triglumi["mc"] = _jp_triglumi[_jp_notrigs-1]/1e6;
  }

  cout << "Calling HistosNormalize("<<_jp_type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Loop over all the directories recursively
  recurseFile(fin, fout);

  cout << endl;
  cout << "Recursive loop done." << endl;
  cout << "Output written in " << fout->GetName() << endl;
  fout->Close();
  cout << "Output file closed" << endl;
  fout->Delete();
  cout << "Output file pointer deleted" << endl << flush;

  fin->Close();
  fin->Delete();
} // HistosNormalize


void recurseFile(TDirectory *indir, TDirectory *outdir, double etawid, double etamid, int lvl) {
  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *obj;
  TKey *key;

  while ( (key = dynamic_cast<TKey*>(itkey.Next())) ) {
    if (_jp_debug) cout << key->GetName() << endl << flush;
    obj = key->ReadObj(); assert(obj);

    // Found a subdirectory: copy it to output and go deeper
    if (obj->InheritsFrom("TDirectory")) {
      int loclvl = lvl;
      outdir->mkdir(obj->GetName());
      bool enteroutdir = outdir->cd(obj->GetName());
      assert(enteroutdir);
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName()); assert(outdir2);
      outdir2->cd();

      bool enterindir = indir->cd(obj->GetName());
      assert(enterindir);
      TDirectory *indir2 = indir->GetDirectory(obj->GetName());
      indir2->cd();

      // Check if directory name contains information on eta bin width
      float etamin, etamax;
      if ( (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2) and (etamax>etamin) ) {
        etawid = 2.*(etamax-etamin);
        etamid = 0.5*(etamax+etamin);
        loclvl++;
      } else if (TString(indir2->GetName()).Contains("FullEta")) {
        loclvl++;
      } else if (loclvl>0) {
        loclvl++;
        cout << endl << " trg: " << indir2->GetName();
      }
      if (loclvl==1) cout << endl << "Entering: " << indir2->GetName();

      recurseFile(indir2, outdir2, etawid, etamid, loclvl);
    } // inherits from TDirectory

    if (obj->InheritsFrom("TH1")) { // Normalize (or just pass on) histogrammish objects
      outdir->cd();
      string name = obj->GetName();
      string trgname = indir->GetName();

      TObject *obj2 = obj->Clone(obj->GetName()); // Copy the input object to output

      if (name=="hlumi" and _jp_usetriglumi) { // hlumi handling
        TH1D *hlumi = dynamic_cast<TH1D*>(indir->Get("hlumi")); assert(hlumi);
        TH1D *hlumi0 = dynamic_cast<TH1D*>(indir->Get(Form("../%s/hlumi",_jp_reftrig))); assert(hlumi0);

        TH1D *hlumi_orig = dynamic_cast<TH1D*>(outdir->FindObject("hlumi_orig"));
        if (!hlumi_orig) hlumi_orig = dynamic_cast<TH1D*>(hlumi->Clone("hlumi_orig"));

        TH1D *hlumi_new = dynamic_cast<TH1D*>(outdir->FindObject("hlumi"));
        if (hlumi_new) hlumi = hlumi_new;

        double lumi = triglumi[trgname];
        for (int i = 1; i != hlumi->GetNbinsX()+1; ++i) hlumi->SetBinContent(i, lumi);

        double lumi0 = triglumi[_jp_reftrig];
        for (int i = 1; i != hlumi0->GetNbinsX()+1; ++i) hlumi0->SetBinContent(i, lumi0);
        continue;
      } // hlumi

      double lumi = 1;
      double lumiref = 1;

      if (!_jp_useversionlumi) {
        if (_jp_usetriglumi) {
          lumi = triglumi[trgname];
          lumiref = triglumi[_jp_reftrig];
        } else if (_jp_isdt or !_jp_usemctrig) {
          TH1D* lumihisto = dynamic_cast<TH1D*>(indir->Get("hlumi"));
          TH1D* lumihistoref = dynamic_cast<TH1D*>(indir->Get(Form("../%s/hlumi",_jp_reftrig)));
          assert(lumihisto);
          assert(lumihistoref);
          lumi = lumihisto->GetBinContent(1);
          lumiref = lumihistoref->GetBinContent(1);
        }
      }

      if (std::find(hptlike.begin(), hptlike.end(), name) != hptlike.end()) { // A block for hpt objects
        TH1D *hpt = dynamic_cast<TH1D*>(obj2);
        bool isgen = TString(obj2->GetName()).Contains("pt_g");
        bool isoth = (TString(obj2->GetName()).Contains("pt_no") ||
                      TString(obj2->GetName()).Contains("djmass") ||
                      TString(obj2->GetName()).Contains("hpt0"));
        bool ispre = (TString(obj2->GetName()).Contains("_pre"));
        bool isjk = (TString(obj2->GetName()).Contains("hpt_jk"));
        bool isjet = (TString(obj2->GetName()).Contains("hpt_jet"));

        TProfile *peff = dynamic_cast<TProfile*>(indir->Get("peff")); assert(peff);

        //// Test MC-based normalization for trigger efficiency
        bool dotrigeff = ((string(obj2->GetName())=="hpt") or isjk or isjet);
        TH1D *htrigeff = dynamic_cast<TH1D*>(outdir->FindObject("htrigeff"));
        TH1D *hpt_notrigeff = 0;
        if (!htrigeff and _jp_dotrigeff and dotrigeff) { // This is done only once
          TH1D *htrigeffmc = dynamic_cast<TH1D*>(outdir->FindObject("htrigeffmc"));
          TH1D *htrigeffsf = dynamic_cast<TH1D*>(outdir->FindObject("htrigeffsf"));
          TFile *fmc = new TFile("output-MC-1.root","READ");
          assert(fmc and !fmc->IsZombie());
          assert(fmc->cd("Standard"));
          fmc->cd("Standard");
          TDirectory *dmc0 = fmc->GetDirectory("Standard");
          TDirectory *dmc = dmc0->GetDirectory(Form("Eta_%1.1f-%1.1f", etamid-0.25*etawid,etamid+0.25*etawid));
          assert(dmc);
          dmc->cd();

          // Add MC truth based trigger efficiency
          if(!htrigeffmc and dmc->cd(indir->GetName())) {
            TDirectory *dir1 = dmc->GetDirectory(indir->GetName()); assert(dir1);
            TH1D *hpty = dynamic_cast<TH1D*>(dir1->Get("hpt")); assert(hpty);
            bool entermcdir = dmc->cd("mc");
            assert(entermcdir);
            TDirectory *dir2 = dmc->GetDirectory("mc"); assert(dir2);
            TH1D *hptx = dynamic_cast<TH1D*>(dir2->Get(Form("hpt_%s",indir->GetName())));

            outdir->cd();
            if (hpty && hptx) htrigeffmc = dynamic_cast<TH1D*>(hpty->Clone("htrigeffmc"));
            if (hpty && hptx) htrigeffmc->Divide(hpty,hptx,1,1,"B");
          }

          // Add data/MC scale factor for trigger efficiency
          if (!(_jp_ismc) && !htrigeffsf) {
            bool enterdir = dmc->cd(indir->GetName());
            assert(enterdir);
            TDirectory *dirmc = dmc->GetDirectory(indir->GetName()); assert(dirmc);
            TProfile *pm = dynamic_cast<TProfile*>(dirmc->Get("ptrigefftp"));
            TProfile *pd = dynamic_cast<TProfile*>(indir->Get("ptrigefftp"));

            outdir->cd();
            if (pm and pd) {
              htrigeffsf = pm->ProjectionX("htrigeffsf");
              htrigeffsf->Divide(pd,pm,1);
            }
          }

          // Combine MC trigger efficiency and scalefactor
          if (htrigeffmc) { // not available for 'mc' directory
            outdir->cd();
            htrigeff = dynamic_cast<TH1D*>(htrigeffmc->Clone("htrigeff"));
            assert(_jp_ismc || htrigeffsf);
            if (_jp_isdt) htrigeff->Multiply(htrigeffsf);

            TH1D *h = dynamic_cast<TH1D*>(indir->Get("hpt"));
            assert(outdir->FindObject("hpt_notrigeff")==0);
            outdir->cd();
            hpt_notrigeff = dynamic_cast<TH1D*>(h->Clone("hpt_notrigeff"));
          }

          fmc->Close();
        } //// dotrigeff

        //// Scale data to account for time dependence
        bool dotimedep = ((string(obj2->GetName())=="hpt") or isjk or isjet);
        TH1D *htimedep = dynamic_cast<TH1D*>(outdir->FindObject("htimedep"));
        TH1D *hpt_notimedep = 0, *hpt_withtimedep = 0;
        double ktime = 1.;
        if (_jp_dotimedep and dotimedep) { // This is done only once
          if (!htimedep) {
            TH1D *htimefit = dynamic_cast<TH1D*>(outdir->FindObject("htimefit"));

            TH1D *h = dynamic_cast<TH1D*>(indir->Get("hpt"));
            TH1D *hsel = dynamic_cast<TH1D*>(indir->Get("hselpt"));
            TH1D *hpre = dynamic_cast<TH1D*>(indir->Get("hpt_pre"));

            outdir->cd();
            if (h) hpt_notimedep = dynamic_cast<TH1D*>(h->Clone("hpt_notimedep"));
            if (hpre and h) htimedep = dynamic_cast<TH1D*>(hpre->Clone("htimedep"));
            if (hpre and h) htimedep->Divide(hpre,h);//,1,1,"B");

            if (htimedep and lumi>0 and lumiref>0) {
              htimedep->Scale(lumi/lumiref); // This sounds incorrect, idea from the old code
            }

            // Find proper pT range and fit
            double minpt = 0.;
            double maxpt = 6500.;
            if (hsel) {
              for (int i = 1; i != hsel->GetNbinsX()+1; ++i) {
                if (hsel->GetBinContent(i)!=0 &&
                    hsel->GetBinLowEdge(i)>=_jp_xmin57) {
                  if (minpt<20) minpt = hsel->GetBinLowEdge(i);
                  maxpt = hsel->GetBinLowEdge(i+1);
                }
              }
            }
            TF1 *ftmp = new TF1("ftmp","[0]",minpt,maxpt);
            ftmp->SetParameter(0,1);
            if (htimedep and htimedep->Integral()>0) htimedep->Fit(ftmp,"QRN");

            if (htimedep and ftmp->GetParameter(0)>0)
              ktime = 1./ftmp->GetParameter(0);

            if (htimedep) {
              outdir->cd();
              htimefit = dynamic_cast<TH1D*>(hsel->Clone("htimefit"));
              hpt_withtimedep = dynamic_cast<TH1D*>(h->Clone("hpt_withtimedep"));

              for (int i = 1; i != htimefit->GetNbinsX()+1; ++i) {

                if (hsel->GetBinContent(i)!=0) {
                  htimefit->SetBinContent(i, ftmp->GetParameter(0));
                  htimefit->SetBinError(i, ftmp->GetParError(0));
                }

                // Calculate with time dependence here to add ktime fit error
                hpt_withtimedep->SetBinContent(i, hpt_notimedep->GetBinContent(i)
                                              * htimefit->GetBinContent(i));
                double err1 = hpt_notimedep->GetBinError(i)
                  / hpt_notimedep->GetBinContent(i);
                double err2 = htimefit->GetBinError(i)
                  / htimefit->GetBinContent(i);
                hpt_withtimedep->SetBinError(i, hpt_notimedep->GetBinContent(i)
                                            * sqrt(pow(err1,2) + pow(err2,2)));
              }
            }
          } else { // We need ktime even if all this has already been calculated
            TH1D *hsel = dynamic_cast<TH1D*>(indir->Get("hselpt"));

            // Find proper pT range and fit
            double minpt = 0.;
            double maxpt = 6500.;
            if (hsel) {
              for (int i = 1; i != hsel->GetNbinsX()+1; ++i) {
                if (hsel->GetBinContent(i)!=0 &&
                  hsel->GetBinLowEdge(i)>=_jp_xmin57) {
                  if (minpt<20) minpt = hsel->GetBinLowEdge(i);
                  maxpt = hsel->GetBinLowEdge(i+1);
                  }
              }
            }

            TF1 *ftmp = new TF1("ftmp","[0]",minpt,maxpt);
            ftmp->SetParameter(0,1);
            if (htimedep and htimedep->Integral()>0) htimedep->Fit(ftmp,"QRN");

            if (htimedep and ftmp->GetParameter(0)>0)
              ktime = 1./ftmp->GetParameter(0);
          }
        } // dotimedep


        if (!(hpt->GetNbinsX()==peff->GetNbinsX() || isoth || isgen) || isoth || isgen) {
          assert(hpt->GetNbinsX()==peff->GetNbinsX() || isoth);
        }

        // Lumi weighting checked in each bin separately
        for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
          // Normalization for bin width in y, pT
          double norm = hpt->GetBinWidth(i) * etawid;
          double trigeff = 1.;
          double pt = hpt->GetBinCenter(i);
          // Normalization for all the common efficiencies
          if (peff->GetBinContent(i)!=0 && !isgen)
            norm *= peff->GetBinContent(i);
          // Test MC-based normalization for trigger efficiency
          if (dotrigeff and htrigeff and _jp_dotrigeff) {
            if (htrigeff->GetBinContent(i)!=0) {
              trigeff = min(1.,max(0.,htrigeff->GetBinContent(i)));
              if (_jp_dotrigefflowptonly and pt>=114) trigeff = 1;
              norm *= trigeff;
            }
          }

          // Normalization for luminosity
          if (lumiref>0 and !isgen) {
            if (ispre) {
              norm *= lumiref;
            } else if (lumi>0)
              norm *= lumi/lumiref;
          }

          // Fix luminosity from .csv VTX to lumiCalc vdM
          if (!_jp_ismc) norm *= _lumiscale;
          // Scale normalization for jackknife
          if (isjk) norm *= 0.9;

          // Let's divide by a magical number. This is totally OK.
          if (_jp_ismc and !_jp_pthatbins) norm /= 2500.; //(xsecw / (sumw * adhocw) ); // equals 2551.;

          // Correct data for time-dependence
          double norm_notime = norm;
          if (dotimedep and htimedep and _jp_dotimedep) norm *= ktime;

          if (!(peff->GetBinContent(i)!=0||hpt->GetBinContent(i)==0 || isgen || isoth || hpt->GetBinCenter(i)<_jp_recopt || hpt->GetBinCenter(i)*cosh(etamid)>3500.)) {
            cerr << "Hist " << hpt->GetName() << " " << indir->GetName() << " pt=" << hpt->GetBinCenter(i) << " etamid = " << etamid << endl << flush;
            assert(peff->GetBinContent(i)!=0||hpt->GetBinContent(i)==0||isgen||hpt->GetBinCenter(i)<_jp_recopt);
          }

          assert(norm!=0);
          hpt->SetBinContent(i, hpt->GetBinContent(i) / norm);
          hpt->SetBinError(i, hpt->GetBinError(i) / norm);
          if (hpt_notrigeff) {
            hpt_notrigeff->SetBinContent(i, hpt_notrigeff->GetBinContent(i)/ norm * trigeff);
            hpt_notrigeff->SetBinError(i, hpt_notrigeff->GetBinError(i)/ norm * trigeff);
          }
          if (hpt_notimedep) {
            hpt_notimedep->SetBinContent(i, hpt_notimedep->GetBinContent(i)/ norm_notime);
            hpt_notimedep->SetBinError(i, hpt_notimedep->GetBinError(i)/ norm_notime);
          }
          if (hpt_withtimedep) { // ktime already applied => use norm_notime
            hpt_withtimedep->SetBinContent(i, hpt_withtimedep->GetBinContent(i)/ norm_notime);
            hpt_withtimedep->SetBinError(i, hpt_withtimedep->GetBinError(i)/ norm_notime);
          }
        } // for i
      } else if (_jp_isdt or (_jp_reweighPU and trgname!="mc")) { // TH3, TH2 and TH1 that is not a hpt object
        // We need to weight also the MC trigsim profiles if PU reweighting is used, as the PU histos contain weight info.
        TString namehandle(name);
        TRegexp re_profile("^p");
        bool isprofile = namehandle.Contains(re_profile);
        if (isprofile) {
          if (obj->InheritsFrom("TProfile")) {
            // For TProfile scale acts in a problematic way but add has scaling properties
            TProfile *handle = dynamic_cast<TProfile*>(obj2);
            handle->Reset();
            handle->Add(dynamic_cast<TProfile*>(obj),lumiref/lumi);
          }
        } else if (!namehandle.Contains("hlumi")) {
          if (obj->InheritsFrom("TH3")) {
            TH3D *handle = dynamic_cast<TH3D*>(obj2);
            handle->Scale(lumiref/lumi);
          } else if (obj->InheritsFrom("TH2")) {
            TH2D *handle = dynamic_cast<TH2D*>(obj2);
            handle->Scale(lumiref/lumi);
          } else if (obj->InheritsFrom("TH1")) {
            TH1D *handle = dynamic_cast<TH1D*>(obj2);
            handle->Scale(lumiref/lumi);
          }
        }
      }
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    } // inherits from TH1
    obj->Delete();
  } // while key

  curdir->cd();
} // recurseFile
