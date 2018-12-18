#define zerobias_cxx
#include "zerobias.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

#include "TFile.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <iostream>

using namespace std;

struct histos {
  TProfile *pmu;
  TProfile *prho;
  TProfile *pnpv;
  
  TH1D *hrawvseta;
  TH1D *hl1vseta;
  TH1D *hjecvseta;
  TH1D *hresvseta;
  TH1D *hrcvseta;
  TH1D *hmcvseta;

  TH1D *hrawvspt;
  TH1D *hl1vspt;
  TH1D *hjecvspt;
  TH1D *hresvspt;
  TH1D *hrcvspt;
  TH1D *hmcvspt;

  TH1D *hrawvsmu;
  TH1D *hl1vsmu;
  TH1D *hjecvsmu;
  TH1D *hresvsmu;
  TH1D *hrcvsmu;
  TH1D *hmcvsmu;
  
  TH1D *hmu00up;
  TH1D *hmu15dw;
  TH1D *hmu1520;
  TH1D *hmu2025;
  TH1D *hmu2530;
  TH1D *hmu30up;

  TProfile *pptcorrovptchvsmu;
  TProfile *pptrawovptchvsmu;
  TProfile *pptrcovptchvsmu;
  TProfile *pptmcovptchvsmu;

  TH2D *h2ptcorrovptchvsmu;
  TH2D *h2ptrawovptchvsmu;
  TH2D *h2ptrcovptchvsmu;
  TH2D *h2ptmcovptchvsmu;
};

void zerobias::Loop()
{
//   In a ROOT session, you can do:
//      root> .L zerobias.C
//      root> zerobias t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();

   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("nJets",1);  // activate branchname
   fChain->SetBranchStatus("jet_pt",1);
   fChain->SetBranchStatus("jet_eta",1);
   fChain->SetBranchStatus("jet_area",1);
   fChain->SetBranchStatus("jet_ch",1);
   fChain->SetBranchStatus("rho",1);
   fChain->SetBranchStatus("mu",1);
   fChain->SetBranchStatus("nPV",1);

   string s = "";
   FactorizedJetCorrector *_JEC;
   FactorizedJetCorrector *_JECRC;
   FactorizedJetCorrector *_JECMC;
   {
     const char *cd = "CondFormats/JetMETObjects/data/";
     //s = Form("%s/Summer16_07Aug2017%s_V5_DATA_L1FastJet_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V5_DATA_L1RC_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017_V5_MC_L1FastJet_AK4PFchs.txt",cd);
     //s = Form("%s/Summer16_07Aug2017_V5_MC_L1RC_AK4PFchs.txt",cd);
     //
     //s = Form("%s/Summer16_07Aug2017%s_V15_MC_L1FastJet_AK4PFchs.txt",cd,"");
     //s = Form("%s/Summer16_07Aug2017%s_V16_DATA_L1FastJet_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V17_DATA_L1FastJet_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V18_DATA_L1FastJet_AK4PFchs.txt",cd,"GH");
     s = Form("%s/Summer16_07Aug2017%s_V14Mikko_L1fix_Data_L1FastJet_AK4PFchs.txt",cd,"GH");
     cout << s << endl << flush;
     JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
     //s = Form("%s/Summer16_07Aug2017%s_V15_MC_L1RC_AK4PFchs.txt",cd,"");
     //s = Form("%s/Summer16_07Aug2017%s_V17_DATA_L1RC_AK4PFchs.txt",cd,"GH");
     s = Form("%s/Summer16_07Aug2017%s_V18_DATA_L1RC_AK4PFchs.txt",cd,"GH");
     cout << s << endl << flush;
     JetCorrectorParameters *l1rc = new JetCorrectorParameters(s);
     s = Form("%s/Summer16_07Aug2017%s_V15_MC_L1FastJet_AK4PFchs.txt",cd,"");
     cout << s << endl << flush;
     JetCorrectorParameters *l1mc = new JetCorrectorParameters(s);
     s = Form("%s/Summer16_07Aug2017%s_V15_MC_L2Relative_AK4PFchs.txt",cd,"");
     //s = Form("%s/Summer16_07Aug2017%s_V16_DATA_L2Relative_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V17_DATA_L2Relative_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V18_DATA_L2Relative_AK4PFchs.txt",cd,"GH");
     cout << s << endl << flush;
     JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
     //s = Form("%s/Summer16_07Aug2017%s_V15_MC_L2L3Residual_AK4PFchs.txt",cd,""); //==1
     //s = Form("%s/Summer16_07Aug2017%s_V16M_DATA_L2L3Residual_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V17_DATA_L2L3Residual_AK4PFchs.txt",cd,"GH");
     //s = Form("%s/Summer16_07Aug2017%s_V17C_DATA_L2L3Residual_AK4PFchs.txt",cd,"GH");
     s = Form("%s/Summer16_07Aug2017%s_V18_DATA_L2L3Residual_AK4PFchs.txt",cd,"GH");
     cout << s << endl << flush;
     JetCorrectorParameters *res = new JetCorrectorParameters(s);
     vector<JetCorrectorParameters> v;
     v.push_back(*l1);
     v.push_back(*l2);
     v.push_back(*res);
     _JEC = new FactorizedJetCorrector(v);
     vector<JetCorrectorParameters> vrc;
     vrc.push_back(*l1rc);
     vrc.push_back(*l2);
     vrc.push_back(*res);
     _JECRC = new FactorizedJetCorrector(vrc);
     vector<JetCorrectorParameters> vmc;
     vmc.push_back(*l1mc);
     vmc.push_back(*l2);
     vmc.push_back(*res);
     _JECMC = new FactorizedJetCorrector(vmc);
   }

   //TFile *fout = new TFile("rootfiles/zerobias.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-data15.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-mc-v15.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-data-v16m.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-data-v17.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-data-v17c.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-mc-v15.root","RECREATE");
   //TFile *fout = new TFile("rootfiles/zerobias-data-v18.root","RECREATE");
   TFile *fout = new TFile("rootfiles/zerobias-data-v18v14fix.root","RECREATE");

   const double x[] = 
     //{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
     {5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40,
      43, 49, 56, 64, 74, 84, 97, 114,
      133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
      548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
      2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}; // Eta_0.0-0.5
   const int nx = sizeof(x)/sizeof(x[0])-1;
   
   TH1D *hptlead = new TH1D("hptlead",";p_{T,lead};N_{jet}",nx,x);
   TH1D *hpttag = new TH1D("hpttag",";p_{T,tag};N_{jet}",nx,x);

   TH1D *hmu0  = new TH1D("hmu0",";#mu;1/N_{evt}",60,0,60);
   TH1D *hrho0 = new TH1D("hrho0",";#rho;1/N_{evt}",60,0,60);
   TH1D *hnpv0 = new TH1D("hnpv0",";N_{PV};1/N_{evt}",60,0,60);

   TProfile *pmu0 = new TProfile("pmu0",";#mu;#LT#mu#GT (GeV)",60,0,60);
   TProfile *prho0 = new TProfile("prho0",";#mu;#LT#rho#GT (GeV)",60,0,60);
   TProfile *pnpv0 = new TProfile("pnpv0",";#mu;#LTN_{PV,good}#GT",60,0,60);

   // Repeat these plots for each pT cut
   int vpt[] = {5, 10, 15, 20, 25, 30, 40, 60, 80, 120};
   const int npt = sizeof(vpt)/sizeof(vpt[0]);
   vector<histos*> vh(npt);
   for (int ipt = 0; ipt != npt; ++ipt) {

     // Create new histograms
     int pt = vpt[ipt];
     TProfile *pmu = new TProfile(Form("pmu_pt%d",pt),";#mu;#LT#mu#GT (GeV)",60,0,60);
     TProfile *prho = new TProfile(Form("prho_pt%d",pt),";#mu;#LT#rho#GT (GeV)",60,0,60);
     TProfile *pnpv = new TProfile(Form("pnpv_pt%d",pt),";#mu;#LTN_{PV,good}#GT",60,0,60);
     
     TH1D *hrawvseta = new TH1D(Form("hrawvseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hl1vseta = new TH1D(Form("hl1vseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hjecvseta = new TH1D(Form("hjecvseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hresvseta = new TH1D(Form("hresvseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hrcvseta = new TH1D(Form("hrcvseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmcvseta = new TH1D(Form("hmcvseta_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);

     TH1D *hrawvspt = new TH1D(Form("hrawvspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     TH1D *hl1vspt = new TH1D(Form("hl1vspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     TH1D *hjecvspt = new TH1D(Form("hjecvspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     TH1D *hresvspt = new TH1D(Form("hresvspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     TH1D *hrcvspt = new TH1D(Form("hrcvspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     TH1D *hmcvspt = new TH1D(Form("hmcvspt_pt%d",pt),";p_{T} (GeV);N_{jet}",nx,x);
     
     TH1D *hrawvsmu = new TH1D(Form("hrawvsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);
     TH1D *hl1vsmu = new TH1D(Form("hl1vsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);
     TH1D *hjecvsmu = new TH1D(Form("hjecvsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);
     TH1D *hresvsmu = new TH1D(Form("hresvsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);
     TH1D *hrcvsmu = new TH1D(Form("hrcvsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);
     TH1D *hmcvsmu = new TH1D(Form("hmcvsmu_pt%d",pt),";#mu;N_{jet}",60,0,60);

     TH1D *hmu00up = new TH1D(Form("hmu00up_%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmu15dw = new TH1D(Form("hmu15dw_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmu1520 = new TH1D(Form("hmu1520_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmu2025 = new TH1D(Form("hmu2025_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmu2530 = new TH1D(Form("hmu2530_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);
     TH1D *hmu30up = new TH1D(Form("hmu30up_pt%d",pt),";#eta;N_{jet}",104,-5.2,5.2);

     TProfile *pptcorrovptchvsmu = new TProfile(Form("pptcorrovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,corr}/p_{T,ch}",60,0,60);
     TProfile *pptrawovptchvsmu = new TProfile(Form("pptrawovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,raw}/p_{T,ch}",60,0,60);
     TProfile *pptrcovptchvsmu = new TProfile(Form("pptrcovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,RC}/p_{T,ch}",60,0,60);
     TProfile *pptmcovptchvsmu = new TProfile(Form("pptmcovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,MC}/p_{T,ch}",60,0,60);

     TH2D *h2ptcorrovptchvsmu = new TH2D(Form("h2ptcorrovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,corr}/p_{T,ch}",60,0,60,200,0,4);
     TH2D *h2ptrawovptchvsmu = new TH2D(Form("h2ptrawovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,raw}/p_{T,ch}",60,0,60,200,0,4);
     TH2D *h2ptrcovptchvsmu = new TH2D(Form("h2ptrcovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,RC}/p_{T,ch}",60,0,60,200,0,4);
     TH2D *h2ptmcovptchvsmu = new TH2D(Form("h2ptmcovptchsvsmu_pt%d",pt),";#mu;0.6#timesp_{T,MC}/p_{T,ch}",60,0,60,200,0,4);

     // Store in a vector array for easy handing;
     histos *h = new histos();
     vh[ipt] = h;

     h->pmu = pmu;
     h->prho = prho;
     h->pnpv = pnpv;
  
     h->hrawvseta = hrawvseta;
     h->hl1vseta = hl1vseta;
     h->hjecvseta = hjecvseta;
     h->hresvseta = hresvseta;
     h->hrcvseta = hrcvseta;
     h->hmcvseta = hmcvseta;

     h->hrawvspt = hrawvspt;
     h->hl1vspt = hl1vspt;
     h->hjecvspt = hjecvspt;
     h->hresvspt = hresvspt;
     h->hrcvspt = hrcvspt;
     h->hmcvspt = hmcvspt;
  
     h->hrawvsmu = hrawvsmu;
     h->hl1vsmu = hl1vsmu;
     h->hjecvsmu = hjecvsmu;
     h->hresvsmu = hresvsmu;
     h->hrcvsmu = hrcvsmu;
     h->hmcvsmu = hmcvsmu;

     h->hmu00up = hmu00up;
     h->hmu15dw = hmu15dw;
     h->hmu1520 = hmu1520;
     h->hmu2025 = hmu2025;
     h->hmu2530 = hmu2530;
     h->hmu30up = hmu30up;

     h->pptcorrovptchvsmu = pptcorrovptchvsmu;
     h->pptrawovptchvsmu = pptrawovptchvsmu;
     h->pptrcovptchvsmu = pptrcovptchvsmu;
     h->pptmcovptchvsmu = pptmcovptchvsmu;

     h->h2ptcorrovptchvsmu = h2ptcorrovptchvsmu;
     h->h2ptrawovptchvsmu = h2ptrawovptchvsmu;
     h->h2ptrcovptchvsmu = h2ptrcovptchvsmu;
     h->h2ptmcovptchvsmu = h2ptmcovptchvsmu;
   }
   double sum(0), sum15dw(0), sum1520(0), sum2025(0), sum2530(0), sum30up(0);
   double nmu(0), nmu15dw(0), nmu1520(0), nmu2025(0), nmu2530(0), nmu30up(0);

   //const double ptmin = 15.;//30

   Long64_t ntot = fChain->GetEntries();
   cout << "Read in " << ntot << " entries" << endl << flush;
   cout << "Process " << 1.*ntot/100000 << " blocks" << endl << flush;

   // Set number of entries
   //nentries = 1000000;
   //nentries = 100;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (ientry%100000==0) cout << "." << flush;

      ++sum; nmu += mu;
      if (mu<15)           { ++sum15dw; nmu15dw += mu; }
      if (mu>=15 && mu<20) { ++sum1520; nmu1520 += mu; }
      if (mu>=20 && mu<25) { ++sum2025; nmu2025 += mu; }
      if (mu>=25 && mu<30) { ++sum2530; nmu2530 += mu; }
      if (mu>=30)          { ++sum30up; nmu30up += mu; }
      
      hmu0->Fill(mu);
      hrho0->Fill(rho);
      hnpv0->Fill(nPV);

      pmu0->Fill(mu,mu);
      prho0->Fill(mu,rho);
      pnpv0->Fill(mu,nPV);

      for (int i = 0; i != nJets; ++i) {

	_JEC->setRho(rho);
	_JEC->setJetA(jet_area[i]);
	_JEC->setJetPt(jet_pt[i]);
	_JEC->setJetEta(jet_eta[i]);
	//double jec = _JEC->getCorrection();
	//double ptcorr = jec*jet_pt[i];
	vector<float> jecs = _JEC->getSubCorrections();

	_JECRC->setRho(rho);
	_JECRC->setJetA(jet_area[i]);
	_JECRC->setJetPt(jet_pt[i]);
	_JECRC->setJetEta(jet_eta[i]);
	double jecrc = _JECRC->getCorrection();

	_JECMC->setRho(rho);
	_JECMC->setJetA(jet_area[i]);
	_JECMC->setJetPt(jet_pt[i]);
	_JECMC->setJetEta(jet_eta[i]);
	double jecmc = _JECMC->getCorrection();

	double ptcorr = jecs[2]*jet_pt[i];
	double ptch = jet_ch[i]*jet_pt[i];

	for (int ipt = 0; ipt != npt; ++ipt) {
	  
	  double ptmin = vpt[ipt];
	  double ptmax = (ipt!=npt-1 ? vpt[ipt+1] : 6500.);
	  histos *h = vh[ipt];

	  // Rates vs eta (for given pT)
	  if (jet_pt[i]>ptmin && jet_pt[i]<ptmax)
	    h->hrawvseta->Fill(jet_eta[i]);
	  if (jecs[0]*jet_pt[i]>ptmin && jecs[0]*jet_pt[i]<ptmax)
	    h->hl1vseta->Fill(jet_eta[i]);
	  if (jecs[1]*jet_pt[i]>ptmin && jecs[1]*jet_pt[i]>ptmin)
	    h->hjecvseta->Fill(jet_eta[i]);
	  if (jecs[2]*jet_pt[i]>ptmin && jecs[2]*jet_pt[i]>ptmin)
	    h->hresvseta->Fill(jet_eta[i]);
	  if (jecrc*jet_pt[i]>ptmin && jecrc*jet_pt[i]>ptmin)
	    h->hrcvseta->Fill(jet_eta[i]);
	  if (jecmc*jet_pt[i]>ptmin && jecmc*jet_pt[i]>ptmin)
	    h->hmcvseta->Fill(jet_eta[i]);

	  if (fabs(jet_eta[i])<1.3) {

	    // Rates vs pT (for barrel eta)
	    h->hrawvspt->Fill(jet_pt[i]);
	    h->hl1vspt->Fill(jecs[0]*jet_pt[i]);
	    h->hjecvspt->Fill(jecs[1]*jet_pt[i]);
	    h->hresvspt->Fill(jecs[2]*jet_pt[i]);
	    h->hrcvspt->Fill(jecrc*jet_pt[i]);
	    h->hmcvspt->Fill(jecmc*jet_pt[i]);

	    // Rates vs mu (for given pT and barrel eta)
	    if (mu>0.1) {
	      if (jet_pt[i]>ptmin && jet_pt[i]<ptmax)
		h->hrawvsmu->Fill(mu, 1./mu);
	      if (jecs[0]*jet_pt[i]>ptmin && jecs[0]*jet_pt[i]<ptmax)
		h->hl1vsmu->Fill(mu, 1./mu);
	      if (jecs[1]*jet_pt[i]>ptmin && jecs[1]*jet_pt[i]>ptmin)
		h->hjecvsmu->Fill(mu, 1./mu);
	      if (jecs[2]*jet_pt[i]>ptmin && jecs[2]*jet_pt[i]>ptmin)
		h->hresvsmu->Fill(mu, 1./mu);
	      if (jecrc*jet_pt[i]>ptmin && jecrc*jet_pt[i]>ptmin)
		h->hrcvsmu->Fill(mu, 1./mu);
	      if (jecmc*jet_pt[i]>ptmin && jecmc*jet_pt[i]>ptmin)
		h->hmcvsmu->Fill(mu, 1./mu);
	    }
	  }

	  // Rates vs eta (for given pT and mu range)
	  if (ptcorr>ptmin && ptcorr<ptmax) {

	    h->pmu->Fill(mu,mu);
	    h->prho->Fill(mu,rho);
	    h->pnpv->Fill(mu,nPV);
	    
	    h->hmu00up->Fill(jet_eta[i]);
	    if (mu<15)           h->hmu15dw->Fill(jet_eta[i]);
	    if (mu>=15 && mu<20) h->hmu1520->Fill(jet_eta[i]);
	    if (mu>=20 && mu<25) h->hmu2025->Fill(jet_eta[i]);
	    if (mu>=25 && mu<30) h->hmu2530->Fill(jet_eta[i]);
	    if (mu>=30)          h->hmu30up->Fill(jet_eta[i]);
	  } // pt>ptmin

	  if (ptch>0.6*ptmin && ptch<0.6*ptmax) {

	    if (fabs(jet_eta[i])<1.3 ) {

	      //if (ptch > 0.2*jet_pt[i]) {
	      h->pptcorrovptchvsmu->Fill(mu, 0.6*ptcorr/ptch);
	      h->pptrawovptchvsmu->Fill(mu, 0.6*jet_pt[i]/ptch);
	      h->pptrcovptchvsmu->Fill(mu, 0.6*jecrc*jet_pt[i]/ptch);
	      h->pptmcovptchvsmu->Fill(mu, 0.6*jecmc*jet_pt[i]/ptch);
	      //}
	      h->h2ptcorrovptchvsmu->Fill(mu, 0.6*ptcorr/ptch);
	      h->h2ptrawovptchvsmu->Fill(mu, 0.6*jet_pt[i]/ptch);
	      h->h2ptrcovptchvsmu->Fill(mu, 0.6*jecrc*jet_pt[i]/ptch);
	      h->h2ptmcovptchvsmu->Fill(mu, 0.6*jecmc*jet_pt[i]/ptch);
	    } // |eta|<1.3
	  } // ptch>0.6*ptmin
	} // for ipt

	if (fabs(jet_eta[i])<1.3) {

	  if (i==0) hptlead->Fill(ptcorr);
	  if (i<2)  hpttag->Fill(ptcorr);
	} // for |jet_eta|<1.3

      } // for i
   }

   hmu0->Scale(1./nentries);
   hrho0->Scale(1./nentries);
   hnpv0->Scale(1./nentries);

   for (int ipt = 0; ipt != npt; ++ipt) {

     histos *h = vh[ipt];
     h->hrawvseta->Scale(1./nmu);
     h->hl1vseta->Scale(1./nmu);
     h->hjecvseta->Scale(1./nmu);
     h->hresvseta->Scale(1./nmu);
     h->hrcvseta->Scale(1./nmu);
     h->hmcvseta->Scale(1./nmu);

     h->hrawvspt->Scale(1./nmu);
     h->hl1vspt->Scale(1./nmu);
     h->hjecvspt->Scale(1./nmu);
     h->hresvspt->Scale(1./nmu);
     h->hrcvspt->Scale(1./nmu);
     h->hmcvspt->Scale(1./nmu);

     h->hrawvsmu->Scale(1./sum);
     h->hl1vsmu->Scale(1./sum);
     h->hjecvsmu->Scale(1./sum);
     h->hresvsmu->Scale(1./sum);
     h->hrcvsmu->Scale(1./sum);
     h->hmcvsmu->Scale(1./sum);

     h->hmu00up->Scale(1./nmu);
     h->hmu15dw->Scale(1./nmu15dw);
     h->hmu1520->Scale(1./nmu1520);
     h->hmu2025->Scale(1./nmu2025);
     h->hmu2530->Scale(1./nmu2530);
     h->hmu30up->Scale(1./nmu30up);
   }

   cout << "Processed nentries="<<nentries<<" (sum="<<sum
	<<", nmu="<<nmu<<", <mu>="<<nmu/sum<<")"<<endl;

   fout->Write();
   fout->Close();
}
