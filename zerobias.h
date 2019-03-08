//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 23 22:10:14 2018 by ROOT version 6.04/10
// from TTree T/Offset Tree
// found on file: moredata/offset/Legacy_jets/Legacy16_RunG_unskimmed_jettuples.root
//////////////////////////////////////////////////////////

#ifndef zerobias_h
#define zerobias_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class zerobias {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           bx;
   ULong64_t       event;
   Float_t         mu;
   Float_t         rho;
   Float_t         rhoC0;
   Float_t         rhoCC;
   Int_t           nPVall;
   Int_t           nPV;
   Float_t         pv_ndof[86];   //[nPVall]
   Float_t         pv_z[86];   //[nPVall]
   Float_t         pv_rho[86];   //[nPVall]
   Float_t         ht;
   Int_t           nJets;
   Float_t         jet_eta[4];   //[nJets]
   Float_t         jet_phi[4];   //[nJets]
   Float_t         jet_pt[4];   //[nJets]
   Float_t         jet_area[4];   //[nJets]
   Float_t         jet_ch[4];   //[nJets]
   Float_t         jet_nh[4];   //[nJets]
   Float_t         jet_ne[4];   //[nJets]
   Float_t         jet_hfh[4];   //[nJets]
   Float_t         jet_hfe[4];   //[nJets]
   Float_t         jet_lep[4];   //[nJets]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_event;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoC0;   //!
   TBranch        *b_rhoCC;   //!
   TBranch        *b_nPVall;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_rho;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_ch;   //!
   TBranch        *b_jet_nh;   //!
   TBranch        *b_jet_ne;   //!
   TBranch        *b_jet_hfh;   //!
   TBranch        *b_jet_hfe;   //!
   TBranch        *b_jet_lep;   //!

   zerobias(TTree *tree=0);
   virtual ~zerobias();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef zerobias_cxx
zerobias::zerobias(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("moredata/offset/Legacy_jets/Legacy16_RunG_unskimmed_jettuples.root");
     //if (!f || !f->IsOpen()) {
     //  f = new TFile("moredata/offset/Legacy_jets/Legacy16_RunG_unskimmed_jettuples.root");
     //}
     //f->GetObject("T",tree);
     exit(0);
   }
   Init(tree);
}

zerobias::~zerobias()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t zerobias::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t zerobias::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void zerobias::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoC0", &rhoC0, &b_rhoC0);
   fChain->SetBranchAddress("rhoCC", &rhoCC, &b_rhoCC);
   fChain->SetBranchAddress("nPVall", &nPVall, &b_nPVall);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("pv_ndof", pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_z", pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_rho", pv_rho, &b_pv_rho);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_ch", jet_ch, &b_jet_ch);
   fChain->SetBranchAddress("jet_nh", jet_nh, &b_jet_nh);
   fChain->SetBranchAddress("jet_ne", jet_ne, &b_jet_ne);
   fChain->SetBranchAddress("jet_hfh", jet_hfh, &b_jet_hfh);
   fChain->SetBranchAddress("jet_hfe", jet_hfe, &b_jet_hfe);
   fChain->SetBranchAddress("jet_lep", jet_lep, &b_jet_lep);
   Notify();

}

Bool_t zerobias::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void zerobias::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t zerobias::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef zerobias_cxx
