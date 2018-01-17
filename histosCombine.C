// Purpose: Combine different triggers into a single spectrum
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: June 2, 2015
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

#include <iostream>
#include <map>
#include <string>

#include "settings.h"

using namespace std;

TH1D *recurseFile(TDirectory *indir, TDirectory *outdir, string hname = "hpt",
                  bool ptSelect = true,bool atBottom = false, TH1D *_hpt = 0, double etamid = 0);
map<string, pair<double, double> > _ptranges;
map<string, map<int, pair<double, double> > > _massranges;

// global variables (not pretty, but works)
TDirectory *_top = 0;

void histosCombine() {
  
  TDirectory *curdir = gDirectory;

  TFile *fin = new TFile(Form("output-%s-2a.root",_jp_type),"READ");
  assert(fin && !fin->IsZombie());
  _top = gDirectory;

  TFile *fout = new TFile(Form("output-%s-2b.root",_jp_type),"RECREATE");
  assert(fout && !fout->IsZombie());

  cout << "Calling histosCombine("<<_jp_type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursions. These may take a few seconds" << endl << flush;

  // Store pT ranges to a nice map
  for (int itrg = 0; itrg != _jp_ntrigs; ++itrg) {

    _ptranges[_jp_triggers[itrg]] =
      pair<double, double>(_jp_trigranges[itrg][0], _jp_trigranges[itrg][1]);

    // 
    if (_jp_ismc && _jp_usemctrig) {
      _ptranges["mc"] = pair<double,double>(0., _jp_sqrts/2.);
      _ptranges[_jp_triggers[itrg]] = pair<double,double>(0.,0.);
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
  if (!_jp_ismc) recurseFile(fin, fout, "hlumi");
  if (!_jp_ismc) recurseFile(fin, fout, "hlumi_orig");

  recurseFile(fin, fout, "hpt_evt");
  recurseFile(fin, fout, "hpt_evtcount");
  recurseFile(fin, fout, "hpt_evt");
  recurseFile(fin, fout, "hpt_jet");

  if (_jp_ismc) recurseFile(fin, fout, "hpt_g0tw");

  recurseFile(fin, fout, "hpt0");
  recurseFile(fin, fout, "hpt1");
  recurseFile(fin, fout, "hpt2");
  recurseFile(fin, fout, "hpt3");
  recurseFile(fin, fout, "hpt_noid");
  recurseFile(fin, fout, "hpt_noevtid");
  recurseFile(fin, fout, "hpt_nojetid");
  recurseFile(fin, fout, "hpt_ak5calo");

  recurseFile(fin, fout, "pa");
  recurseFile(fin, fout, "pnpv");
  recurseFile(fin, fout, "pnpvall");
  recurseFile(fin, fout, "prho");
  recurseFile(fin, fout, "prho1");
  recurseFile(fin, fout, "prho2");
  recurseFile(fin, fout, "prho3");
  //
  recurseFile(fin, fout, "pchs");
  recurseFile(fin, fout, "pchsx");
  recurseFile(fin, fout, "pchs1");
  recurseFile(fin, fout, "pchs2");
  recurseFile(fin, fout, "pchs3");
  //
  recurseFile(fin, fout, "pa_ak5calo");
  recurseFile(fin, fout, "prho_ak5calo");
  recurseFile(fin, fout, "prho1_ak5calo");
  recurseFile(fin, fout, "prho2_ak5calo");
  recurseFile(fin, fout, "prho3_ak5calo");

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
  //recurseFile(fin, fout, "pjec_res");
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
  // TODO: We could add a combination method for the vsnpv stuff.
  // At the moment jt40 is used as an approximation, since its weight dominates

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

  recurseFile(fin, fout, "hdjmass0");
  recurseFile(fin, fout, "hdjmass0_hgg");

  // JEC stability checks
  recurseFile(fin, fout, "hpt_plus");
  recurseFile(fin, fout, "hpt0_plus");
  recurseFile(fin, fout, "hpt_minus");
  recurseFile(fin, fout, "hpt0_minus");
  recurseFile(fin, fout, "hpt_l1fast");
  recurseFile(fin, fout, "hpt_l1off");

  if (!_jp_ismc) recurseFile(fin, fout, "peff_new");

  recurseFile(fin, fout, "pemf");
  recurseFile(fin, fout, "pemftp");
  recurseFile(fin, fout, "pcalopf");
  recurseFile(fin, fout, "pcalopftp");
  recurseFile(fin, fout, "pcalotp");
  recurseFile(fin, fout, "ppftp");
  recurseFile(fin, fout, "pemftp25");
  recurseFile(fin, fout, "pcalopftp25");
  recurseFile(fin, fout, "pcalotp25");
  recurseFile(fin, fout, "ppftp25");

  recurseFile(fin, fout, "pak4ak8_25");
  recurseFile(fin, fout, "pak4ak8_50");
  recurseFile(fin, fout, "pak4ak8tp_25");
  recurseFile(fin, fout, "pak4ak8tp_50");
  recurseFile(fin, fout, "hak4ak8");
  recurseFile(fin, fout, "hak4ak8tp");

  curdir->cd();

  cout << endl;
  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();
  cout << "Output file closed" << endl;

  // For some reason closing input is very slow. Quit instead of waiting
  //return; 

  fin->Close();
  fin->Delete();
  cout << "Input file closed" << endl;

} // histosNormalize


TH1D* recurseFile(TDirectory *indir, TDirectory *outdir, string hname,
                  bool ptSelect, bool atBottom, TH1D *_hpt, double etamid) {

  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    // We want the key to match to a directory or the selected histogram
    string classname = ((TKey*)key)->GetClassName();
    string kname = key->GetName();
    if (classname=="TDirectoryFile" or kname==hname) {
      obj = ((TKey*)key)->ReadObj();
      assert(obj);
    } else {
      continue;
    }

    if (classname=="TDirectoryFile" and obj->InheritsFrom("TDirectory")) {
      // Found a subdirectory: copy it to output and go deeper

      TDirectory *outdir2 = outdir;
      if (!atBottom) {
        if (outdir->FindKey(obj->GetName())==0) outdir->mkdir(obj->GetName());
        assert(outdir->cd(key->GetName()));
        outdir2 = outdir->GetDirectory(key->GetName()); assert(outdir2);
        outdir2->cd();
        if (_jp_debug) cout << key->GetName() << endl;
      } else if (_jp_debug) {
        cout << key->GetName() << " (at bottom)" << endl;      
      }

      assert(indir->cd(obj->GetName()));
      TDirectory *indir2 = indir->GetDirectory(obj->GetName()); assert(indir2);
      indir2->cd();

      // Check if directory name contains information on eta bin width
      // If yes, the next level is the bottom level with triggers
      // set flag and reset the combined histogram pointer
      float etamin, etamax;
      if ( (sscanf(indir2->GetName(),"Eta_%f-%f",&etamin,&etamax)==2)
         && (etamax>etamin) ) {

        _hpt = recurseFile(indir2, outdir2, hname, ptSelect, true, _hpt,
               0.5*(etamin+etamax));
        if (_hpt) {
          outdir2->cd();
          // Write out to save ram
          _hpt->Write();
          _hpt = 0;
        }
      } else if (TString(indir2->GetName()).Contains("FullEta")) {
        _hpt = recurseFile(indir2, outdir2, hname, ptSelect, true, _hpt, etamid);
        if (_hpt) {
          outdir2->cd();
          // Write out to save ram
          _hpt->Write();
          _hpt = 0;
        }
      } else {
        _hpt = recurseFile(indir2, outdir2, hname, ptSelect, atBottom, _hpt, etamid);
      }
      // inherits from TDirectory
    } else if (kname==hname) {
      // Flatten TProfile to TH1D for later processing
      if (obj->InheritsFrom("TProfile")) {
        if (_hpt==0) {
          outdir->cd();
          _hpt = ((TProfile*)obj)->ProjectionX(obj->GetName());
          indir->cd();
        }
      } // TProfile
      // Flatter TProfile3D to TH3D for later processing
      if (obj->InheritsFrom("TProfile3D")) {
        if (_hpt==0) {
          outdir->cd();
          _hpt = (TH1D*)((TProfile3D*)obj)->ProjectionXYZ(obj->GetName());
          indir->cd();
        }
      } // TProfile3D

      // Copy over TH3, TH2, TH1 histograms, in this precise order
      // Careful with if-then, because TH3 inherits from TH1+TH2
      // Clone and reset histogram the first time it is seen
      if (obj->InheritsFrom("TH3")) {
        TH3D *hpt3 = (TH3D*)obj;
        
        if (_hpt==0) {
          outdir->cd();
          _hpt = (TH1D*)hpt3->Clone(hpt3->GetName());
          _hpt->Reset();
          indir->cd();
          assert(_hpt);
          if (_jp_debug) cout << "Cloned _" << hpt3->GetName() << endl;
        }
        
        if (_ptranges.find(indir->GetName())==_ptranges.end())
          cout << "pT range not found for directory " << indir->GetName() << endl;
        assert(_ptranges.find(indir->GetName())!=_ptranges.end());
        double ptmin = _ptranges[indir->GetName()].first;
        double ptmax = _ptranges[indir->GetName()].second;
        
        TH3D *_hpt3 = (TH3D*)_hpt;
        if (ptSelect) {
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
        // TH3
      } else if (obj->InheritsFrom("TH2")) {
        TH2D *hpt2 = (TH2D*)obj;
        
        if (_hpt==0) {
          outdir->cd();
          _hpt = (TH1D*)hpt2->Clone(hpt2->GetName());
          _hpt->Reset();
          indir->cd();
          assert(_hpt);
          if (_jp_debug) cout << "Cloned _" << hpt2->GetName() << endl;
        }
        
        if (_ptranges.find(indir->GetName())==_ptranges.end())
          cout << "pT range not found for directory " << indir->GetName() << endl;
        assert(_ptranges.find(indir->GetName())!=_ptranges.end());
        double ptmin = _ptranges[indir->GetName()].first;
        double ptmax = _ptranges[indir->GetName()].second;
        
        TH2D *_hpt2 = (TH2D*)_hpt;
        if (ptSelect) {
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
        // TH2
      } else if (obj->InheritsFrom("TH1")) {
        TH1D *hpt = (TH1D*)obj;
        
        if (_hpt==0) {
          outdir->cd();
          _hpt = (TH1D*)hpt->Clone(hpt->GetName());
          _hpt->Reset();
          indir->cd();
          assert(_hpt);
          if (_jp_debug) cout << "Cloned _" << hpt->GetName() << endl;
        }
        
        if (_ptranges.find(indir->GetName())==_ptranges.end())
          cout << "pT range not found for directory " << indir->GetName() << endl;
        assert(_ptranges.find(indir->GetName())!=_ptranges.end());
        double ptmin = _ptranges[indir->GetName()].first;
        double ptmax = _ptranges[indir->GetName()].second;
        
        // Replace ranges for mass histograms
        if (TString(hname.c_str()).Contains("hdjmass")) {
          assert(etamid!=0);
          int ieta = int(etamid/0.5); assert(ieta<=7);
          ptmin = _massranges[indir->GetName()][ieta].first;
          ptmax = _massranges[indir->GetName()][ieta].second;
          if (ptmin==0) ptmin = 10.;
          if (ptmax==0) ptmax = 3000.;
        } // mass histo
        
        if (ptSelect) {
          for (int i = 1; i != _hpt->GetNbinsX()+1; ++i) {
            double pt = _hpt->GetBinCenter(i);
            if (pt > ptmin && pt < ptmax) { // TODO: We could do better than a linear search!
              _hpt->SetBinContent(i, hpt->GetBinContent(i));
              _hpt->SetBinError(i, hpt->GetBinError(i));
            } // in ptrange
          } // for i
        } else {
          _hpt->Add(hpt);
        }
      } // TH1D
    }

    // Free memory, avoid malloc error
    obj->Delete();
  } // while key

  // Free memory
  if (indir==_top) {
  //outdir->Write();
    cout << "." << flush;
  }
  curdir->cd();

  return _hpt;
} // recurseFile
