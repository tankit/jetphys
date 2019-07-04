// Purpose: Fuse Single Neutrino and standard MC jet analysis histograms
// Author: hannu.siikonen@cern.ch
// Created: December 3, 2018


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

using namespace std;

#include "settings.h"


void FuseRecurseFile(TDirectory *indir, TDirectory *indir_nu, TDirectory *outdir);

void HistosFuse(string type = "")
{
  if (type=="") {
    type = jp::type;
    if (!jp::ispy and !jp::ishw) return;
  }

  TFile *fin = new TFile(Form("output-%s-2a.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());
  TFile *fin_nu = new TFile("output-NU-2a.root","READ");
  assert(fin_nu && !fin_nu->IsZombie());

  TFile *fout = new TFile(Form("output-%sNU-2a.root",type.c_str()),"RECREATE");
  assert(fout and !fout->IsZombie());

  cout << "Calling HistosFuse("<<type<<");" << endl;
  cout << "Input file " << fin->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Loop over all the directories recursively
  FuseRecurseFile(fin, fin_nu, fout);

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


void FuseRecurseFile(TDirectory *indir, TDirectory *indir_nu, TDirectory *outdir) {
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
      outdir->mkdir(obj->GetName());
      bool enteroutdir = outdir->cd(obj->GetName());
      assert(enteroutdir);
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName());
      assert(outdir2);
      outdir2->cd();

      bool enterindir = indir->cd(obj->GetName());
      assert(enterindir);
      bool enterindir_nu = indir_nu->cd(obj->GetName());
      assert(enterindir_nu);

      TDirectory *indir2;
      TDirectory *indir2_nu;

      cout << obj->GetName();
      // For the ZeroBias trigger we do a goode ole switcheroo'
      if (TString(obj->GetName()).Contains("jt0")) {
        indir2 = indir_nu->GetDirectory(obj->GetName());
        indir2_nu = indir->GetDirectory(obj->GetName());
        cout << " Using NU info!";
      } else {
        indir2 = indir->GetDirectory(obj->GetName());
        indir2_nu = indir_nu->GetDirectory(obj->GetName());
        if (TString(obj->GetName()).Contains("jt")) cout << " Using original info!";
      }
      cout << endl;
      indir2->cd();

      FuseRecurseFile(indir2, indir2_nu, outdir2);
    } else {
      outdir->cd();

      TObject *obj2 = obj->Clone(obj->GetName()); // Copy the input object to output

      // Save the stuff into an identical directory
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    } // inherits from TH1
    obj->Delete();
  } // while key

  curdir->cd();
} // FuseRecurseFile
