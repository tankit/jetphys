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

void recurseNormFile(TDirectory *indir, TDirectory *outdir, bool isdt,
                 double etawid = 1., double etamid = 0., int lvl = 0);

// Use this to fix luminosity
std::map<std::string, double> triglumi;

void HistosNormalize(string type = "")
{
  bool isdt = jp::isdt;
  if (type=="") type = jp::type;
  else          isdt = (type=="DATA");
  TFile *fin = new TFile(Form("output-%s-1.root",type.c_str()),"READ");
  assert(fin and !fin->IsZombie());

  TFile *fout = new TFile(Form("output-%s-2a.root",type.c_str()),"RECREATE");
  assert(fout and !fout->IsZombie());

  if (isdt and jp::usetriglumi) { // Setting up lumis
    cout << "Reading trigger luminosity from settings.h" << endl;
    int eraIdx = -1;
    if (jp::usetriglumiera) {
      int eraNo = 0;
      for (auto &eraMatch : jp::eras) {
        if (std::regex_search(jp::run,eraMatch)) {
          eraIdx = eraNo;
          break;
        }
        ++eraNo;
      }
      if (eraIdx!=-1) cout << "Using weights according to the run era!" << endl;
      else cout << "Could not locate the given era! :(" << endl;
    }
    for (unsigned int i = 0; i < jp::notrigs; ++i) {
      double lumi = (jp::usetriglumiera ? jp::triglumiera[eraIdx][i]/1e6 : jp::triglumi[i]/1e6); // /ub to /pb
      cout << Form(" *%s: %1.3f /pb", jp::triggers[i],lumi) << endl;
      triglumi[jp::triggers[i]] = lumi;
    }
  }

  cout << "Calling HistosNormalize("<<type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Loop over all the directories recursively
  recurseNormFile(fin, fout, isdt);

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


void recurseNormFile(TDirectory *indir, TDirectory *outdir, bool isdt, double etawid, double etamid, int lvl) {
  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *obj;
  TKey *key;

  while ( (key = dynamic_cast<TKey*>(itkey.Next())) ) {
    if (jp::debug) cout << key->GetName() << endl << flush;
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

      recurseNormFile(indir2, outdir2, isdt, etawid, etamid, loclvl);
      // inherits from TDirectory
    } else if (obj->InheritsFrom("TH1")) { // Normalize (or just pass on) histogrammish objects
      outdir->cd();
      string name = obj->GetName();
      string trgname = indir->GetName();

      TObject *obj2 = obj->Clone(obj->GetName()); // Copy the input object to output

      double lumi = 1;
      double lumiref = 1;

      if (isdt) {
        if (jp::usetriglumi) { // Use lumi info from settings.h
          if (name=="hlumi") {
            // Overwrite hlumi
            TH1D *hlumi = dynamic_cast<TH1D*>(indir->Get("hlumi")); assert(hlumi);
            TH1D *hlumi0 = dynamic_cast<TH1D*>(indir->Get(Form("../%s/hlumi",jp::reftrig))); assert(hlumi0);

            TH1D *hlumi_orig = dynamic_cast<TH1D*>(outdir->FindObject("hlumi_orig"));
            if (!hlumi_orig) hlumi_orig = dynamic_cast<TH1D*>(hlumi->Clone("hlumi_orig"));

            TH1D *hlumi_new = dynamic_cast<TH1D*>(outdir->FindObject("hlumi"));
            if (hlumi_new) hlumi = hlumi_new;

            double lumi = triglumi[trgname];
            for (int i = 1; i != hlumi->GetNbinsX()+1; ++i) hlumi->SetBinContent(i, lumi);

            double lumi0 = triglumi[jp::reftrig];
            for (int i = 1; i != hlumi0->GetNbinsX()+1; ++i) hlumi0->SetBinContent(i, lumi0);
            continue;
          } // hlumi

          // Set lumi weights
          lumi = triglumi[trgname];
          lumiref = triglumi[jp::reftrig];
        } else { // Use lumi info from hlumi
          TH1D* lumihisto = dynamic_cast<TH1D*>(indir->Get("hlumi"));
          TH1D* lumihistoref = dynamic_cast<TH1D*>(indir->Get(Form("../%s/hlumi",jp::reftrig)));
          assert(lumihisto);
          assert(lumihistoref);

          // Set lumi weights
          lumi = lumihisto->GetBinContent(1);
          lumiref = lumihistoref->GetBinContent(1);
        }
      }

      if (std::find(hptlike.begin(), hptlike.end(), name) != hptlike.end()) { // A block for hpt objects
        TH1D *hpt = dynamic_cast<TH1D*>(obj2);
        double norm0 = etawid;

        // Let's divide by a magical number. This is totally OK (disabled).
        // if (!isdt and !jp::pthatbins) norm0 /= 2500.; //(xsecw / (sumw * adhocw) ); // equals 2551.;

        bool isgen = TString(obj2->GetName()).Contains("pt_g");
        bool isoth = (TString(obj2->GetName()).Contains("pt_no") ||
                      TString(obj2->GetName()).Contains("djmass") ||
                      TString(obj2->GetName()).Contains("hpt0"));
        bool isjk = (TString(obj2->GetName()).Contains("hpt_jk"));
        bool isjet = (TString(obj2->GetName()).Contains("hpt_jet"));
        bool hptstuff = (string(obj2->GetName())=="hpt") or isjk or isjet;

        // Normalization for luminosity
        if (isdt and lumiref>0) {
          if (TString(obj2->GetName()).Contains("_pre")) norm0 *= lumiref;
          else if (lumi>0)                               norm0 *= lumi/lumiref;
        }

        // Scale normalization for jackknife (but why?)
        if (isjk) norm0 *= 0.9;

        TProfile *peff;
        if (jp::dotrigeffsimple) {
          peff = dynamic_cast<TProfile*>(indir->Get("peff"));
          assert(peff);
        }

        //// Test MC-based normalization for trigger efficiency
        TH1D *htrigeff = dynamic_cast<TH1D*>(outdir->FindObject("htrigeff"));
        TH1D *hpt_notrigeff = 0;
        if (jp::dotrigeff) { // Requires output-MC-1.root, seldom in use
          if (hptstuff) {
            if (!htrigeff) { // This is done only once
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
            if (!htrigeffmc and dmc->cd(indir->GetName())) {
              TDirectory *dir1 = dmc->GetDirectory(indir->GetName()); // Use the given trigger
              assert(dir1);
              TH1D *hpty = dynamic_cast<TH1D*>(dir1->Get("hpt")); assert(hpty);

              bool entermcdir = dmc->cd("mc"); // Use the 'mc' trigger
              assert(entermcdir);
              TDirectory *dir2 = dmc->GetDirectory("mc");
              assert(dir2);
              TH1D *hptx = dynamic_cast<TH1D*>(dir2->Get(Form("hpt_%s",indir->GetName())));

              outdir->cd();
              if (hpty and hptx) {
                htrigeffmc = dynamic_cast<TH1D*>(hpty->Clone("htrigeffmc"));
                htrigeffmc->Divide(hpty,hptx,1,1,"B");
              }
            }

            // Add data/MC scale factor for trigger efficiency
            if (isdt and !htrigeffsf) {
              bool enterdir = dmc->cd(indir->GetName());
              assert(enterdir);
              TDirectory *dirmc = dmc->GetDirectory(indir->GetName());
              assert(dirmc);
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
              assert(!isdt or htrigeffsf);
              if (isdt) htrigeff->Multiply(htrigeffsf);

              TH1D *h = dynamic_cast<TH1D*>(indir->Get("hpt"));
              assert(outdir->FindObject("hpt_notrigeff")==0);
              outdir->cd();
              hpt_notrigeff = dynamic_cast<TH1D*>(h->Clone("hpt_notrigeff"));
            }

            fmc->Close();
            } // set trigeff (once)
          } // trigeff only for some histos
        } // jp::dotrigeff

        //// Scale data to account for time dependence
        TH1D *htimedep = dynamic_cast<TH1D*>(outdir->FindObject("htimedep"));
        TH1D *hpt_notimedep = 0, *hpt_withtimedep = 0;
        double ktime = 0.0;
        if (jp::dotimedep) { // Seldom used
          if (hptstuff) {
            if (!htimedep) { // This is done only once
              TH1D *htimefit = dynamic_cast<TH1D*>(outdir->FindObject("htimefit"));

              TH1D *h = dynamic_cast<TH1D*>(indir->Get("hpt"));
              TH1D *hsel = dynamic_cast<TH1D*>(indir->Get("hselpt"));
              TH1D *hpre = dynamic_cast<TH1D*>(indir->Get("hpt_pre"));

              outdir->cd();
              if (h) hpt_notimedep = dynamic_cast<TH1D*>(h->Clone("hpt_notimedep"));
              if (hpre and h) htimedep = dynamic_cast<TH1D*>(hpre->Clone("htimedep"));
              if (hpre and h) htimedep->Divide(hpre,h);//,1,1,"B");

              if (htimedep and lumi>0 and lumiref>0)
                htimedep->Scale(lumi/lumiref); // This sounds incorrect, idea from the old code

              // Find proper pT range and fit
              double minpt = 0.;
              double maxpt = 6500.;
              if (hsel) {
                for (int i = 1; i != hsel->GetNbinsX()+1; ++i) {
                  if (hsel->GetBinContent(i)!=0 and hsel->GetBinLowEdge(i)>=jp::xmin57) {
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
                    hsel->GetBinLowEdge(i)>=jp::xmin57) {
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
            } // fetch ktime
          } // timedep only for some histos
        } // jp::dotimedep

        if (jp::dotrigeffsimple) assert(hpt->GetNbinsX()==peff->GetNbinsX() or isoth);

        // Lumi weighting checked in each bin separately
        for (int binidx = 1; binidx != hpt->GetNbinsX()+1; ++binidx) {
          // Normalization for bin width in y, pT
          double norm = hpt->GetBinWidth(binidx) * norm0;
          double trigeff = 1.;
          double pt = hpt->GetBinCenter(binidx);
          // Normalization for all the common efficiencies (by default, off)
          if (jp::dotrigeffsimple and !isgen and peff->GetBinContent(binidx)!=0) norm *= peff->GetBinContent(binidx);

          // Test MC-based normalization for trigger efficiency
          if (htrigeff) {
            if (htrigeff->GetBinContent(binidx)!=0) {
              trigeff = min(1.,max(0.,htrigeff->GetBinContent(binidx))); // Efficiency between 0 and 1
              if (jp::dotrigefflowptonly and pt>=114) trigeff = 1;
              norm *= trigeff;
            }
          }

          // Correct data for time-dependence
          double norm_notime = norm;
          if (ktime>0.0) norm *= ktime;

          assert(norm!=0);
          hpt->SetBinContent(binidx, hpt->GetBinContent(binidx) / norm);
          hpt->SetBinError(binidx, hpt->GetBinError(binidx) / norm);

          if (hpt_notrigeff) {
            hpt_notrigeff->SetBinContent(binidx, (hpt_notrigeff->GetBinContent(binidx)/ norm) * trigeff);
            hpt_notrigeff->SetBinError(binidx, (hpt_notrigeff->GetBinError(binidx)/ norm) * trigeff);
          }
          if (hpt_notimedep) {
            hpt_notimedep->SetBinContent(binidx, hpt_notimedep->GetBinContent(binidx)/ norm_notime);
            hpt_notimedep->SetBinError(binidx, hpt_notimedep->GetBinError(binidx)/ norm_notime);
          }
          if (hpt_withtimedep) { // ktime already applied => use norm_notime
            hpt_withtimedep->SetBinContent(binidx, hpt_withtimedep->GetBinContent(binidx)/ norm_notime);
            hpt_withtimedep->SetBinError(binidx, hpt_withtimedep->GetBinError(binidx)/ norm_notime);
          }

          // By default off, done for some data histos, require Zero peff but non-zero hpt and a good pt & eta region
          if (jp::dotrigeffsimple and !isgen and !isoth and peff->GetBinContent(binidx)==0 and hpt->GetBinContent(binidx)!=0 and hpt->GetBinCenter(binidx)>jp::recopt and hpt->GetBinCenter(binidx)*cosh(etamid)<3500.)
            cerr << "Hist " << hpt->GetName() << " " << indir->GetName() << " pt=" << hpt->GetBinCenter(binidx) << " etamid = " << etamid << endl << flush;
        } // for binidx
      } else if (isdt) { // TH3, TH2 and TH1 that is not a hpt object
        // This we do only for data - there are no lumi weights to be applied for MC
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
      // Save the stuff into an identical directory
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
      // inherits from TH1
    } else { // Others
      // Save the stuff into an identical directory
      TObject *obj2 = obj->Clone(obj->GetName()); // Copy the input object to output
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    }
    obj->Delete();
  } // while key

  curdir->cd();
} // recurseNormFile
