// Purpose: Combine different triggers into a single spectrum
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: June 2, 2015

// We throw a name to recurseFile and the function loops over all directories.
// As we reach a certain level, a new histogram is declared.
// Onto this histogram we add all the weighted products of the jtX subdirectories.
// The output directory is used as a temporary saving place for the histos, since pointers in ROOT are problematic.

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TRegexp.h"

#include <iostream>
#include <map>
#include <string>
#include <set>

#include "settings.h"

using namespace std;

void recurseFile(TDirectory *indir, TDirectory *outdir, string hname = "hpt", bool ptSelect = true, 
                 bool othProf = false, int lvl = 0, double etamid = 0);
map<string, pair<double, double> > _ptranges;
map<string, double> _trigwgts;

// global variables (not pretty, but works)
TDirectory *_top = 0;

void histosCombine() {
  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-2a.root",_jp_type),"READ");
  assert(fin and !fin->IsZombie());
  _top = gDirectory;

  TFile *fout = new TFile(Form("output-%s-2b.root",_jp_type),"RECREATE");
  assert(fout and !fout->IsZombie());

  cout << "Calling histosCombine("<<_jp_type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursions. These may take a few seconds" << endl << flush;

  // Store pT ranges to a nice map
  for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {
    _ptranges[_jp_triggers[itrg]] = pair<double, double>(_jp_trigranges[itrg][0], _jp_trigranges[itrg][1]);

    if (_jp_ismc) { // When in mc, we need to know how to use trigger weighting. Either mc or triggers are ignored
      _ptranges["mc"] = pair<double,double>(0., _jp_sqrts/2.);
      if (_jp_usemctrig) {
        _ptranges[_jp_triggers[itrg]] = pair<double,double>(0.,0.);
      } else if (_jp_domctrigsim) {
        _ptranges["mc"] = pair<double,double>(0., 0.);
      }
    } else if (_jp_usetriglumi) {
      _trigwgts[_jp_triggers[itrg]] = _jp_triglumi[_jp_ntrigs-1]/_jp_triglumi[itrg];
    }
  }

  // Loop over all the directories recursively
  // List here the histograms that need merging
  recurseFile(fin, fout, "hpt");
  recurseFile(fin, fout, "hpt_pre");
  recurseFile(fin, fout, "hpt_notrigeff");
  recurseFile(fin, fout, "htrigeff");
  recurseFile(fin, fout, "htrigeffmc");
  recurseFile(fin, fout, "htrigeffsf");
  recurseFile(fin, fout, "hpt_notimedep");
  recurseFile(fin, fout, "hpt_withtimedep");
  recurseFile(fin, fout, "htimedep");
  recurseFile(fin, fout, "htimefit");
  if (_jp_isdt) recurseFile(fin, fout, "hlumi");
  if (_jp_isdt) recurseFile(fin, fout, "hlumi_orig");

  recurseFile(fin, fout, "hpt_evt");
  recurseFile(fin, fout, "hpt_evtcount");
  recurseFile(fin, fout, "hpt_jet");

  if (_jp_ismc) recurseFile(fin, fout, "hpt_g0tw");

  recurseFile(fin, fout, "hpt0");
  recurseFile(fin, fout, "hpt1");
  recurseFile(fin, fout, "hpt2");
  recurseFile(fin, fout, "hpt3");
  recurseFile(fin, fout, "hpt_noid");
  recurseFile(fin, fout, "hpt_noevtid");
  recurseFile(fin, fout, "hpt_nojetid");

  recurseFile(fin, fout, "pa");
  recurseFile(fin, fout, "pnpv");
  recurseFile(fin, fout, "pnpvall");
  recurseFile(fin, fout, "prho");
  //
  recurseFile(fin, fout, "pchs");
  //
  recurseFile(fin, fout, "pmpf");
  recurseFile(fin, fout, "pmpf1");
  recurseFile(fin, fout, "pmpf2");
  recurseFile(fin, fout, "pmpfx");
  recurseFile(fin, fout, "pmpfy");
  recurseFile(fin, fout, "pmpfz");
  recurseFile(fin, fout, "punc");
  //
  recurseFile(fin, fout, "pncand");
  recurseFile(fin, fout, "pnch");
  recurseFile(fin, fout, "pnne");
  recurseFile(fin, fout, "pnnh");
  recurseFile(fin, fout, "pnce");
  recurseFile(fin, fout, "pnmu");
  recurseFile(fin, fout, "pchf");
  recurseFile(fin, fout, "pnef");
  recurseFile(fin, fout, "pnhf");
  recurseFile(fin, fout, "pcef");
  recurseFile(fin, fout, "pmuf");
  recurseFile(fin, fout, "pbeta");
  recurseFile(fin, fout, "pbetastar");
  recurseFile(fin, fout, "pjec");
  recurseFile(fin, fout, "pjec2");
  recurseFile(fin, fout, "pjec_res");
  //
  recurseFile(fin, fout, "pncandtp");
  recurseFile(fin, fout, "pnchtp");
  recurseFile(fin, fout, "pnnetp");
  recurseFile(fin, fout, "pnnhtp");
  recurseFile(fin, fout, "pncetp");
  recurseFile(fin, fout, "pnmutp");
  recurseFile(fin, fout, "pchftp");
  recurseFile(fin, fout, "pneftp");
  recurseFile(fin, fout, "pnhftp");
  recurseFile(fin, fout, "pceftp");
  recurseFile(fin, fout, "pmuftp");
  recurseFile(fin, fout, "pbetatp");
  recurseFile(fin, fout, "pbetastartp");
  recurseFile(fin, fout, "pjectp");
  //
  recurseFile(fin, fout, "pchftp_vsnpv", false, true);
  recurseFile(fin, fout, "pneftp_vsnpv", false, true);
  recurseFile(fin, fout, "pnhftp_vsnpv", false, true);
  recurseFile(fin, fout, "pceftp_vsnpv", false, true);
  recurseFile(fin, fout, "pmuftp_vsnpv", false, true);
  recurseFile(fin, fout, "pbetatp_vsnpv", false, true);
  recurseFile(fin, fout, "pbetastartp_vsnpv", false, true);
  //
  recurseFile(fin, fout, "pchftp_vstrpu", false, true);
  recurseFile(fin, fout, "pneftp_vstrpu", false, true);
  recurseFile(fin, fout, "pnhftp_vstrpu", false, true);
  recurseFile(fin, fout, "pceftp_vstrpu", false, true);
  recurseFile(fin, fout, "pmuftp_vstrpu", false, true);
  recurseFile(fin, fout, "pbetatp_vstrpu", false, true);
  recurseFile(fin, fout, "pbetastartp_vstrpu", false, true);
  //
  recurseFile(fin, fout, "pchftp_vseta", false, true);
  recurseFile(fin, fout, "pneftp_vseta", false, true);
  recurseFile(fin, fout, "pnhftp_vseta", false, true);
  recurseFile(fin, fout, "pceftp_vseta", false, true);
  recurseFile(fin, fout, "pmuftp_vseta", false, true);
  recurseFile(fin, fout, "pbetatp_vseta", false, true);
  recurseFile(fin, fout, "pbetastartp_vseta", false, true);

  recurseFile(fin, fout, "hdjasymm");
  recurseFile(fin, fout, "hdjasymmtp");
  recurseFile(fin, fout, "hdjasymmpt");
  recurseFile(fin, fout, "hdjmpf");
  recurseFile(fin, fout, "hdjmpftp");
  recurseFile(fin, fout, "hdjmpfpt");

  recurseFile(fin, fout, "hdjasymm_a005");
  recurseFile(fin, fout, "hdjasymmtp_a005");
  recurseFile(fin, fout, "hdjasymmpt_a005");
  recurseFile(fin, fout, "hdjmpf_a005");
  recurseFile(fin, fout, "hdjmpftp_a005");
  recurseFile(fin, fout, "hdjmpfpt_a005");
  recurseFile(fin, fout, "hdjasymm_a01");
  recurseFile(fin, fout, "hdjasymmtp_a01");
  recurseFile(fin, fout, "hdjasymmpt_a01");
  recurseFile(fin, fout, "hdjmpf_a01");
  recurseFile(fin, fout, "hdjmpftp_a01");
  recurseFile(fin, fout, "hdjmpfpt_a01");
  recurseFile(fin, fout, "hdjasymm_a015");
  recurseFile(fin, fout, "hdjasymmtp_a015");
  recurseFile(fin, fout, "hdjasymmpt_a015");
  recurseFile(fin, fout, "hdjmpf_a015");
  recurseFile(fin, fout, "hdjmpftp_a015");
  recurseFile(fin, fout, "hdjmpfpt_a015");
  recurseFile(fin, fout, "hdjasymm_a02");
  recurseFile(fin, fout, "hdjasymmtp_a02");
  recurseFile(fin, fout, "hdjasymmpt_a02");
  recurseFile(fin, fout, "hdjmpf_a02");
  recurseFile(fin, fout, "hdjmpftp_a02");
  recurseFile(fin, fout, "hdjmpfpt_a02");
  recurseFile(fin, fout, "hdjasymm_a025");
  recurseFile(fin, fout, "hdjasymmtp_a025");
  recurseFile(fin, fout, "hdjasymmpt_a025");
  recurseFile(fin, fout, "hdjmpf_a025");
  recurseFile(fin, fout, "hdjmpftp_a025");
  recurseFile(fin, fout, "hdjmpfpt_a025");
  recurseFile(fin, fout, "hdjasymm_a03");
  recurseFile(fin, fout, "hdjasymmtp_a03");
  recurseFile(fin, fout, "hdjasymmpt_a03");
  recurseFile(fin, fout, "hdjmpf_a03");
  recurseFile(fin, fout, "hdjmpftp_a03");
  recurseFile(fin, fout, "hdjmpfpt_a03");

  recurseFile(fin, fout, "hdjresp_tag_a005");
  recurseFile(fin, fout, "hdjresptp_tag_a005");
  recurseFile(fin, fout, "hdjresp_tag_a01");
  recurseFile(fin, fout, "hdjresptp_tag_a01");
  recurseFile(fin, fout, "hdjresp_tag_a015");
  recurseFile(fin, fout, "hdjresptp_tag_a015");
  recurseFile(fin, fout, "hdjresp_tag_a02");
  recurseFile(fin, fout, "hdjresptp_tag_a02");
  recurseFile(fin, fout, "hdjresp_tag_a025");
  recurseFile(fin, fout, "hdjresptp_tag_a025");
  recurseFile(fin, fout, "hdjresp_tag_a03");
  recurseFile(fin, fout, "hdjresptp_tag_a03");
  
  recurseFile(fin, fout, "hdjresp_probe_a005");
  recurseFile(fin, fout, "hdjresptp_probe_a005");
  recurseFile(fin, fout, "hdjresp_probe_a01");
  recurseFile(fin, fout, "hdjresptp_probe_a01");
  recurseFile(fin, fout, "hdjresp_probe_a015");
  recurseFile(fin, fout, "hdjresptp_probe_a015");
  recurseFile(fin, fout, "hdjresp_probe_a02");
  recurseFile(fin, fout, "hdjresptp_probe_a02");
  recurseFile(fin, fout, "hdjresp_probe_a025");
  recurseFile(fin, fout, "hdjresptp_probe_a025");
  recurseFile(fin, fout, "hdjresp_probe_a03");
  recurseFile(fin, fout, "hdjresptp_tag_a03");

  recurseFile(fin, fout, "hdjmass", false);
  recurseFile(fin, fout, "hdjmass0", false);

  if (!_jp_ismc) recurseFile(fin, fout, "peff_new");

  curdir->cd();

  cout << endl;
  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();
  cout << "Output file closed" << endl;

  fin->Close();
  fin->Delete();
  cout << "Input file closed" << endl;
} // histosNormalize


void recurseFile(TDirectory *indir, TDirectory *outdir, string hname, bool ptSelect,
                 bool othProf, int lvl, double etamid) {

  TDirectory *curdir = gDirectory;
  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TObject *obj = 0;
  TKey *key = 0;

  if (lvl<2) {
    TListIter itkey(keys);
    while ((key = dynamic_cast<TKey*>(itkey.Next()))) {
      // We want the key to match to a directory or the selected histogram
      string classname = key->GetClassName();
      string kname = key->GetName();
      if (classname=="TDirectoryFile") { // All keys are looped over only in search for subdirs
        obj = key->ReadObj(); assert(obj);
      } else {
        continue;
      }
      if (classname=="TDirectoryFile" and obj->InheritsFrom("TDirectory")) {
        // Found a subdirectory: copy it to output and go deeper

        TDirectory *outdir2 = outdir;
        if (lvl == 0) {
          if (outdir->FindKey(obj->GetName())==0) outdir->mkdir(obj->GetName());
          bool enterkeydir = outdir->cd(key->GetName());
          assert(enterkeydir);
          outdir2 = outdir->GetDirectory(key->GetName()); assert(outdir2);
          outdir2->cd();
          if (_jp_debug) cout << key->GetName() << endl;
        } else {
          if (_jp_debug) cout << key->GetName() << " (at bottom)" << endl;
          lvl += 1;
        }

        bool enterobjdir = indir->cd(obj->GetName());
        assert(enterobjdir);
        TDirectory *indir2 = indir->GetDirectory(obj->GetName()); assert(indir2);
        indir2->cd();

        // Check if directory name contains information on eta bin width
        // If yes, the next level is the bottom level with triggers
        // set flag and reset the combined histogram pointer
        float etamin, etamax;
        if (lvl==0 and (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2) and (etamax>etamin) )
          recurseFile(indir2, outdir2, hname, ptSelect, othProf, 1, 0.5*(etamin+etamax));
        else if (lvl==0 and TString(indir2->GetName()).Contains("FullEta"))
          recurseFile(indir2, outdir2, hname, ptSelect, othProf, 1, etamid);
        else
          recurseFile(indir2, outdir2, hname, ptSelect, othProf, lvl, etamid);
      }
      obj->Delete();
    } // while key
  } else {
    // Search for the histogram. This is more cost-effective and avoids using multiple cycles of the same object.
    key = dynamic_cast<TKey*>(keys->FindObject(hname.c_str()));
    if (key) {
      obj = key->ReadObj();
      assert(obj);
      if (obj) {
        int prevcycle = -1;

        if (othProf) {
          if (obj->InheritsFrom("TProfile")) {
            TProfile *ppt = dynamic_cast<TProfile*>(obj);
            TProfile *_ppt = dynamic_cast<TProfile*>(outdir->FindObject(hname.c_str()));

            if (_ppt) {
              TKey* juu = dynamic_cast<TKey*>(outdir->FindKey(hname.c_str()));
              prevcycle = juu->GetCycle();
            } else {
              outdir->cd();
              _ppt = dynamic_cast<TProfile*>(ppt->Clone(obj->GetName())); assert(_ppt);
              _ppt->Reset();
              indir->cd();
              if (_jp_debug) cout << "Cloned _" << obj->GetName() << endl;
            }
            double wgt = 1.0;
            if (_jp_isdt and _jp_usetriglumi) wgt = _trigwgts[indir->GetName()];
            _ppt->Add(ppt,wgt);
            outdir->cd();
            _ppt->Write();
          }
        } else {
          // Copy over TH3, TH2, TH1 histograms, in this precise order
          // Careful with if-then, because TH3 inherits from TH1+TH2
          // Clone and reset histogram the first time it is seen

          if (obj->InheritsFrom("TH3")) {
            TH3D *hpt3 = dynamic_cast<TH3D*>(obj);
            TH3D *_hpt3 = dynamic_cast<TH3D*>(outdir->FindObject(hname.c_str()));

            if (_hpt3) {
              TKey* juu = dynamic_cast<TKey*>(outdir->FindKey(hname.c_str()));
              prevcycle = juu->GetCycle();
            } else {
              outdir->cd();
              _hpt3 = dynamic_cast<TH3D*>(hpt3->Clone(obj->GetName())); assert(_hpt3);
              _hpt3->Reset();
              indir->cd();
              if (_jp_debug) cout << "Cloned _" << obj->GetName() << endl;
            }

            if (ptSelect) {
              if (_ptranges.find(indir->GetName())==_ptranges.end())
                cout << "pT range not found for directory " << indir->GetName() << endl;
              assert(_ptranges.find(indir->GetName())!=_ptranges.end());
              double ptmin = _ptranges[indir->GetName()].first;
              double ptmax = _ptranges[indir->GetName()].second;

              for (int i = 1; i != _hpt3->GetNbinsX()+1; ++i) {
                double pt = _hpt3->GetXaxis()->GetBinCenter(i); 
                if (pt > ptmin && pt < ptmax) { // TODO: We could do better than a linear search!
                  for (int j = 1; j != _hpt3->GetNbinsY()+1; ++j) {
                    for (int k = 1; k != _hpt3->GetNbinsZ()+1; ++k) {
                      _hpt3->SetBinContent(i,j,k, hpt3->GetBinContent(i,j,k));
                      _hpt3->SetBinError(i,j,k, hpt3->GetBinError(i,j,k));
                    } // for l
                  } // for j
                  break;
                } // in ptrange
              } // for i
            } else {
              _hpt3->Add(hpt3);
            }
            outdir->cd();
            _hpt3->Write();
          } else if (obj->InheritsFrom("TH2")) {
            TH2D *hpt2 = dynamic_cast<TH2D*>(obj);
            TH2D *_hpt2 = dynamic_cast<TH2D*>(outdir->FindObject(hname.c_str()));

            if (_hpt2) {
              TKey* juu = dynamic_cast<TKey*>(outdir->FindKey(hname.c_str()));
              prevcycle = juu->GetCycle();
            } else {
              outdir->cd();
              _hpt2 = dynamic_cast<TH2D*>(hpt2->Clone(obj->GetName())); assert(_hpt2);
              _hpt2->Reset();
              indir->cd();
              if (_jp_debug) cout << "Cloned _" << obj->GetName() << endl;
            }

            if (ptSelect) {
              if (_ptranges.find(indir->GetName())==_ptranges.end())
                cout << "pT range not found for directory " << indir->GetName() << endl;
              assert(_ptranges.find(indir->GetName())!=_ptranges.end());
              double ptmin = _ptranges[indir->GetName()].first;
              double ptmax = _ptranges[indir->GetName()].second;

              for (int i = 1; i != _hpt2->GetNbinsX()+1; ++i) {
                double pt = _hpt2->GetXaxis()->GetBinCenter(i);
                if (pt > ptmin && pt < ptmax) { // TODO: We could do better than a linear search!
                  for (int j = 1; j != _hpt2->GetNbinsY()+1; ++j) {
                    _hpt2->SetBinContent(i,j, hpt2->GetBinContent(i,j));
                    _hpt2->SetBinError(i,j, hpt2->GetBinError(i,j));
                  } // for j
                } // in ptrange
              } // for i
            } else {
              _hpt2->Add(hpt2);
            }
            outdir->cd();
            _hpt2->Write();
          } else if (obj->InheritsFrom("TH1")) {
            TH1D *hpt = 0;
            TH1D *_hpt1 = dynamic_cast<TH1D*>(outdir->FindObject(hname.c_str()));
            if (obj->InheritsFrom("TProfile")) {
              hpt = (dynamic_cast<TProfile*>(obj))->ProjectionX(Form("%s_%s",obj->GetName(),indir->GetName()));
            } else {
              hpt = dynamic_cast<TH1D*>(obj);
            }
            if (_hpt1) {
              TKey* juu = dynamic_cast<TKey*>(outdir->FindKey(hname.c_str()));
              prevcycle = juu->GetCycle();
            } else {
              outdir->cd();
              _hpt1 = dynamic_cast<TH1D*>(hpt->Clone(obj->GetName())); assert(_hpt1);
              _hpt1->Reset();
              indir->cd();
              if (_jp_debug) cout << "Cloned _" << obj->GetName() << endl;
            }

            if (ptSelect) {
              if (_ptranges.find(indir->GetName())==_ptranges.end())
                cout << "pT range not found for directory " << indir->GetName() << endl;
              assert(_ptranges.find(indir->GetName())!=_ptranges.end());
              double ptmin = _ptranges[indir->GetName()].first;
              double ptmax = _ptranges[indir->GetName()].second;

              for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
                double pt = hpt->GetBinCenter(i);
                if (pt > ptmin and pt < ptmax) { // TODO: We could do better than a linear search!
                  _hpt1->SetBinContent(i, hpt->GetBinContent(i));
                  _hpt1->SetBinError(i, hpt->GetBinError(i));
                } // in ptrange
              } // for i
            } else {
              _hpt1->Add(hpt);
            }
            if (obj->InheritsFrom("TProfile")) outdir->Delete(hpt->GetName());
            outdir->cd();
            _hpt1->Write();
          } // TH
        } // othProf
        if (prevcycle!=-1) outdir->Delete(Form("%s;%d",hname.c_str(),prevcycle));
        obj->Delete();
      }
    }
  }
  if (indir==_top)
    cout << "." << flush;
  curdir->cd();
} // recurseFile
