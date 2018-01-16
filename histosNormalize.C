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
                          "hpt_l1off",
                          "hpt_l1fast",
                          "hpt_plus",
                          "hpt_minus",
                          "hpt0_plus",
                          "hpt0_minus",
                          "hpt_noid",
                          "hpt_noevtid",
                          "hpt_nojetid",
                          "hpt_ak5calo",
                          "hpt_ak5pf",
                          "hpt_evt_ak5pf",
                          "hpt_jet_ak5pf",
                          "hselpt",
                          "hpt_r",
                          "hpt_g",
                          "hpt_gg",
                          "hpt_g0",
                          "hpt_g0tw",
                          "hdjmass",
                          "hdjmass0",
                          "hdjmass0_hgg"};

void recurseFile(TDirectory *indir, TDirectory *outdir,
                 double etawid = 1., double etamid = 0.);

// Use this to fix luminosity
double _lumiscale = 1.00;
std::map<std::string, double> triglumi;

void histosNormalize() {

  TFile *fin = new TFile(Form("output-%s-1.root",_jp_type),"READ");
  assert(fin && !fin->IsZombie());

  TFile *fout = new TFile(Form("output-%s-2a.root",_jp_type),"RECREATE");
  assert(fout && !fout->IsZombie());

  if (_lumiscale!=1 && !_jp_ismc)
    cout << "Attention! : Scaling luminosity to the new estimate"
         << " by multiplying with " << _lumiscale << endl;

  if (!(_jp_ismc) && _jp_usetriglumi) {
    cout << "Reading trigger luminosity from settings.h" << endl;
    for (int i = 0; i != _jp_ntrigs; ++i) {
      double lumi = _jp_triglumi[i]/1e6; // /ub to /pb
      cout << Form(" *%s: %1.3f /pb", _jp_triggers[i],lumi) << endl;
      triglumi[_jp_triggers[i]] = lumi;
    }
  }

  cout << "Calling histosNormalize("<<_jp_type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Loop over all the directories recursively
  recurseFile(fin, fout);

  cout << endl;
  cout << "Recursive loop done." << endl;
  cout << "Writing output to " << fout->GetName() << endl;
  cout << "This may again take a minute" << endl << flush;
  fout->Write();
  cout << "Output written in " << fout->GetName() << endl;
  fout->Close();
  cout << "Output file closed" << endl;
  fout->Delete();
  cout << "Output file pointer deleted" << endl << flush;


  fin->Close();
  fin->Delete();

} // histosNormalize


void recurseFile(TDirectory *indir, TDirectory *outdir,
                 double etawid, double etamid) {

  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;
  TDirectory *dir;

  while ( (key = itkey.Next()) ) {

    if (_jp_debug) cout << key->GetName() << endl << flush;
    obj = ((TKey*)key)->ReadObj(); assert(obj);
    dir = indir;

    // Found a subdirectory: copy it to output and go deeper
    if (obj->InheritsFrom("TDirectory")) {

      //assert(outdir->mkdir(obj->GetName()));
      outdir->mkdir(obj->GetName());
      assert(outdir->cd(obj->GetName()));
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName()); assert(outdir2);
      outdir2->cd();

      assert(indir->cd(obj->GetName()));
      TDirectory *indir2 = indir->GetDirectory(obj->GetName());
      indir2->cd();

      // Check if directory name contains information on eta bin width
      float etamin, etamax;
      if ( (sscanf(indir->GetName(),"Eta_%f-%f",&etamin,&etamax)==2)
           && (etamax>etamin) ) {
        etawid = 2.*(etamax-etamin);
        etamid = 0.5*(etamax+etamin);
        //cout << "Eta bin width: " << etawid << flush << endl;
      }

      recurseFile(indir2, outdir2, etawid, etamid);
      //outdir2->Write(); // does this speedup or slow down?
    } // inherits from TDirectory

    // Normalize all histograms
    if (obj->InheritsFrom("TH1")) {

      outdir->cd();
      string name = obj->GetName();
      string trgname = dir->GetName();

//////// Special handling for hlumi
      if (name=="hlumi" && !_jp_ismc && _jp_usetriglumi) {
        TH1D *hlumi = (TH1D*)dir->Get("hlumi");
        assert(hlumi);

        TH1D *hlumi0 = (TH1D*)dir->Get("../jt450/hlumi");
        assert(hlumi0);

        TH1D *hlumi_orig = (TH1D*)outdir->FindObject("hlumi_orig");
        if (!hlumi_orig) hlumi_orig = (TH1D*)hlumi->Clone("hlumi_orig");

        TH1D *hlumi_new = (TH1D*)outdir->FindObject("hlumi");
        if (hlumi_new) hlumi = hlumi_new;

        double lumi = triglumi[trgname];
        for (int i = 1; i != hlumi->GetNbinsX()+1; ++i) {
          hlumi->SetBinContent(i, lumi);
        }

        double lumi0 = triglumi["jt450"];
        for (int i = 1; i != hlumi0->GetNbinsX()+1; ++i) {
          hlumi0->SetBinContent(i, lumi0);
        }

        continue;
      }
//////// hlumi

      TObject *obj2 = obj->Clone(obj->GetName());
      
      double lumi = 1;
      double lumiref = 1;
      if (!_jp_ismc) {
        if (_jp_usetriglumi) {
          lumi = triglumi[trgname];
          lumiref = triglumi["jt450"];
        } else {
          TH1D* lumihisto = (TH1D*)dir->Get("hlumi");
          TH1D* lumihistoref = (TH1D*)dir->Get("../jt450/hlumi");
          assert(lumihisto);
          assert(lumihistoref);
          lumi = lumihisto->GetBinContent(1);
          lumiref = lumihistoref->GetBinContent(1);
        }
      }

      cout << "." << flush;

//////// A block for hpt objects
      if (std::find(hptlike.begin(), hptlike.end(), name) != hptlike.end()) {

        TH1D *hpt = (TH1D*)obj2;
        bool isgen = TString(obj2->GetName()).Contains("pt_g");
        bool isoth = (TString(obj2->GetName()).Contains("pt_no") ||
                      TString(obj2->GetName()).Contains("djmass") ||
                      TString(obj2->GetName()).Contains("hpt0") ||
                      TString(obj2->GetName()).Contains("l1off") ||
                      TString(obj2->GetName()).Contains("l1fast"));
        bool iscalo = (TString(obj2->GetName()).Contains("_ak5calo"));
        bool ispf5 = (TString(obj2->GetName()).Contains("_ak5pf"));
        bool ispre = (TString(obj2->GetName()).Contains("_pre"));
        bool isjk = (TString(obj2->GetName()).Contains("hpt_jk"));
        bool isjet = (TString(obj2->GetName()).Contains("hpt_jet"));

        TProfile *peff = (TProfile*)dir->Get("peff");
        assert(peff);

        
        //// Test MC-based normalization for trigger efficiency
        bool dotrigeff = ((string(obj2->GetName())=="hpt") || isjk || isjet);
        TH1D *htrigeff = (TH1D*)outdir->FindObject("htrigeff");
        TH1D *htrigeffmc = (TH1D*)outdir->FindObject("htrigeffmc");
        TH1D *htrigeffsf = (TH1D*)outdir->FindObject("htrigeffsf");
        TH1D *hpt_notrigeff = 0;

        if (!htrigeff && _jp_dotrigeff && dotrigeff) {

          TFile *fmc = new TFile("output-MC-1.root","READ");
          assert(fmc && !fmc->IsZombie());
          assert(fmc->cd("Standard"));
          fmc->cd("Standard");
          TDirectory *dmc0 = fmc->GetDirectory("Standard");
          TDirectory *dmc = dmc0->GetDirectory(Form("Eta_%1.1f-%1.1f",
                                                    etamid-0.25*etawid,etamid+0.25*etawid));
          assert(dmc);
          dmc->cd();

          // Add MC truth based trigger efficiency
          if(!htrigeffmc && dmc->cd(dir->GetName())) {

            TDirectory *dir1 = dmc->GetDirectory(dir->GetName()); assert(dir1);
            TH1D *hpty = (TH1D*)dir1->Get("hpt"); assert(hpty);
            assert(dmc->cd("mc"));
            dmc->cd("mc");
            TDirectory *dir2 = dmc->GetDirectory("mc"); assert(dir2);
            TH1D *hptx = (TH1D*)dir2->Get(Form("hpt_%s",dir->GetName()));

            outdir->cd();
            if (hpty && hptx) htrigeffmc = (TH1D*)hpty->Clone("htrigeffmc");
            if (hpty && hptx) htrigeffmc->Divide(hpty,hptx,1,1,"B");
          }

          // Add data/MC scale factor for trigger efficiency
          if (!(_jp_ismc) && !htrigeffsf) {

            assert(dmc->cd(dir->GetName()));
            dmc->cd(dir->GetName());
            TDirectory *dirmc = dmc->GetDirectory(dir->GetName()); assert(dirmc);
            TProfile *pm = (TProfile*)dirmc->Get("ptrigefftp");
            TProfile *pd = (TProfile*)dir->Get("ptrigefftp");

            outdir->cd();
            if (pm && pd) htrigeffsf = pm->ProjectionX("htrigeffsf");
            if (pm && pd) htrigeffsf->Divide(pd,pm,1);
          }

          // Combine MC trigger efficiency and scalefactor
          if (htrigeffmc) { // not available for 'mc' directory
            outdir->cd();
            htrigeff = (TH1D*)htrigeffmc->Clone("htrigeff");
            assert(!!(_jp_ismc) || htrigeffsf);
            if (!(_jp_ismc)) htrigeff->Multiply(htrigeffsf);

            TH1D *h = (TH1D*)dir->Get("hpt");
            assert(outdir->FindObject("hpt_notrigeff")==0);
            outdir->cd();
            hpt_notrigeff = (TH1D*)h->Clone("hpt_notrigeff");
          }

          fmc->Close();
        } //// dotrigeff


        //// Scale data to account for time dependence
        bool dotimedep = ((string(obj2->GetName())=="hpt") || isjk || isjet);
        TH1D *htimedep = (TH1D*)outdir->FindObject("htimedep");
        TH1D *htimefit = (TH1D*)outdir->FindObject("htimefit");
        TH1D *hpt_notimedep = 0, *hpt_withtimedep = 0;
        double ktime = 1.;

        if (!htimedep && _jp_dotimedep && dotimedep) {

          TH1D *h = (TH1D*)dir->Get("hpt");
          TH1D *hsel = (TH1D*)dir->Get("hselpt");
          TH1D *hpre = (TH1D*)dir->Get("hpt_pre");

          outdir->cd();
          if (h) hpt_notimedep = (TH1D*)h->Clone("hpt_notimedep");
          if (hpre && h) htimedep = (TH1D*)hpre->Clone("htimedep");
          if (hpre && h) htimedep->Divide(hpre,h);//,1,1,"B");

          if (htimedep && lumi>0 && lumiref>0) {
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
          if (htimedep && htimedep->Integral()>0) htimedep->Fit(ftmp,"QRN");

          if (htimedep && ftmp->GetParameter(0)>0)
            ktime = 1./ftmp->GetParameter(0);

          if (htimedep) {
            outdir->cd();
            htimefit = (TH1D*)hsel->Clone("htimefit");
            hpt_withtimedep = (TH1D*)h->Clone("hpt_withtimedep");

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
        } // dotimedep


        if (!(hpt->GetNbinsX()==peff->GetNbinsX() || isoth || isgen) || isoth || isgen) {
          assert(hpt->GetNbinsX()==peff->GetNbinsX() || isoth);
        }

        // Lumi weighting and unnecessary stuff
        for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

          // Normalization for bin width in y, pT
          double norm = hpt->GetBinWidth(i) * etawid;
          double trigeff = 1.;
          double pt = hpt->GetBinCenter(i);
          // Normalization for all the common efficiencies
          if (peff->GetBinContent(i)!=0 && !isgen)
            norm *= peff->GetBinContent(i);
          // Test MC-based normalization for trigger efficiency
          if (dotrigeff && htrigeff && _jp_dotrigeff) {
            if (htrigeff->GetBinContent(i)!=0) {
              trigeff = min(1.,max(0.,htrigeff->GetBinContent(i)));
              if (_jp_dotrigefflowptonly && pt>=114) trigeff = 1;
              norm *= trigeff;
            }
          }

          // Normalization for luminosity
          // Normalization for luminosity
          if (lumi>0 && lumiref>0 && !isoth && !isgen && !ispre)
            norm *= lumi/lumiref;
          if (lumi>0 && isoth && !isgen && !ispre)
            norm *= lumi;
          if (lumiref && !isoth && !isgen && ispre)
            norm *= lumiref;

          // Fix luminosity from .csv VTX to lumiCalc vdM
          if (!_jp_ismc) norm *= _lumiscale;
          // Scale normalization for jackknife
          if (isjk) norm *= 0.9;

          if (_jp_ismc && _jp_pthatbins) norm *= 1.;
          // Let's divide by a magical number. This is totally OK.
          if (_jp_ismc && !_jp_pthatbins) {
            norm /= 2500.; //(xsecw / (sumw * adhocw) ); // equals 2551.;
          }

          // Correct data for time-dependence
          double norm_notime = norm;
          if (dotimedep && htimedep && _jp_dotimedep) {
            norm *= ktime;
          }

          if (!(peff->GetBinContent(i)!=0||hpt->GetBinContent(i)==0 || isgen || iscalo || ispf5 || isoth || hpt->GetBinCenter(i)<_jp_recopt || hpt->GetBinCenter(i)*cosh(etamid)>3500.)) {
            cerr << "Hist " << hpt->GetName() << " " << dir->GetName() << " pt=" << hpt->GetBinCenter(i) << " etamid = " << etamid << endl << flush;
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
//////// hpt objects
      } else {
        // This is enough. For MC, lumi = lumiref = 1.
        TString namehandle(name);
        TRegexp re_profile("^p");
        if (!namehandle.Contains(re_profile) && !namehandle.Contains("hlumi")) {
          if (obj->InheritsFrom("TH3")) {
            TH3D *handle = (TH3D*)obj2;
            handle->Scale(lumiref/lumi);
          } else if (obj->InheritsFrom("TH2")) {
            TH2D *handle = (TH2D*)obj2;
            handle->Scale(lumiref/lumi);
          } else if (obj->InheritsFrom("TH1")) {
            TH1D *handle = (TH1D*)obj2;
            handle->Scale(lumiref/lumi);
          }
        }
      }

      dir->cd();
    } // inherits from TH1

  } // while key

  curdir->cd();
} // recurseFile
