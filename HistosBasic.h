// Purpose:  Define basic histograms for jet physics studies
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 9, 2015
#ifndef __HistosBasic_h__
#define __HistosBasic_h__

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectory.h"

#include <string>
#include <map>
#include <vector>

class HistosBasic {

 public:

  // phase space
  std::string trigname;
  double etamin;
  double etamax;
  double pttrg;
  double ptmin;
  double ptmax;
  bool ismcdir;

  // raw spectrum
  TH1D *hpt;
  TH1D *hpt_pre;
  TH1D *hcopt;
  TH1D *hpt_noid;
  TH1D *hpt_nojetid;
  TH1D *hpt_noevtid;
  TH1D *hptevt;
  TH1D *hpttmp;

  // delete-m jackknife
  std::vector<TH1D*> hpt_jk;
  TH2D *h2jk;

  TH1D *hpt_tmp;
  TH1D *hpt_evtcount;
  TH1D *hpt_evt;
  TH1D *hpt_jet;

  // 1 GeV ins
  TH1D *hpt0;

  // leading and non-leading jets
  TH1D *hpt1;
  TH1D *hpt2;
  TH1D *hpt3;

  // dijet mass
  TH1D *hdjmass;
  TH1D *hdjmass0;
  TH1D *hdjmass_a01;
  TH1D *hdjmass_a02;
  TH1D *hdjmass_a03;
  TProfile *pdjmass_ptratio;
  TProfile *pdjmass0_ptratio;
  // jet mass (for the two leading jets)
  TH1D *hjmass;
  TH1D *hjmass0;
  TH1D *hjmass_a01;
  TH1D *hjmass_a02;
  TH1D *hjmass_a03;

  // Ozlem: new histograms for quark/gluon study
  TH1D *hgpt;
  TH1D *hgpt0;

  TH1D *hgpt_g0tw;
  TH1D *hqpt_g0tw;
  TH1D *hupt_g0tw; // unmatched
  TH1D *hgpt_g0;
  TH1D *hqpt_g0;
  TH1D *hupt_g0; // unmatched 
  TH1D *hgpt_g; // gluons vs get pT
  TH1D *hqpt_g;
  TH1D *hupt_g; // unmatched vs gen pT
  //
  TH1D *hqgl;
  TH1D *hqgl_q;
  TH1D *hqgl_g;
  TH1D *hqgl_u;
  TH1D *hqgl_dq;
  TH1D *hqgl_dg;
  TH1D *hqgl_du;
  //
  TH2D *hqgl2;
  TH2D *hqgl2_g;
  TH2D *hqgl2_q;
  TH2D *hqgl2_u;
  TH2D *hqgl2_dg;
  TH2D *hqgl2_dq;
  TH2D *hqgl2_du;

  // basic properties
  TProfile *ppt;
  TProfile *pmass;
  TH1D *hmass;
  TProfile *pjec;
  TProfile *pjec2;
  TProfile *punc;
  TH1D *hnpvgood;
  TH1D *hrho;

  TProfile *pjec_l1;
  TProfile *pjec_l2l3;
  TProfile *pjec_res;

  TH1D *htrpu;
  TH1D *hitpu;
  TH1D *hootpuearly;
  TH1D *hootpulate;
  TH2D *h2itvsoot;

  // pile-up information
  TProfile *pa;
  TProfile *ptrpu;
  TProfile *pnpv;
  TProfile *pnpvall;
  TProfile *prho;
  TProfile *pnpvvsrho;
  TProfile *prhovsnpv;
  TProfile *prhovsnpvall;
  TH2D *h2rhovsnpv;
  //
  TProfile *prhovstrpu;
  TProfile *pnpvvstrpu;
  TProfile *pnpvallvstrpu;
  TProfile *pitpuvstrpu;
  TH1D *htrpu2;
  TH1D *hjet_vstrpu;
  TH1D *hlumi_vstrpu;

  // luminosity
  TH1D *hlumi;
  TH1D *hlumi2;
  std::map<int, std::map<int, float> > lums;
  double lumsum;
  double lumsum2;

  // inclusive efficiencies
  TProfile *peff;
  TProfile *pideff;
  TProfile *pvtxeff;
  TProfile *pdqmeff;

  // control plots of components (JEC)
  TProfile *pncand;
  TProfile *pnch;
  TProfile *pnne;
  TProfile *pnnh;
  TProfile *pnce;
  TProfile *pnmu;
  TProfile *pnhh;
  TProfile *pnhe;
  TProfile *pchf;
  TProfile *pnef;
  TProfile *pnhf;
  TProfile *pcef;
  TProfile *pmuf;
  TProfile *phhf;
  TProfile *phef;
  TProfile *ppuf;
  TH1D *hncand;
  TH1D *hnch;
  TH1D *hnne;
  TH1D *hnnh;
  TH1D *hnce;
  TH1D *hnmu;
  TH1D *hnhh;
  TH1D *hnhe;
  TH1D *hchf;
  TH1D *hnef;
  TH1D *hnhf;
  TH1D *hcef;
  TH1D *hmuf;
  TH1D *hhhf;
  TH1D *hhef;
  TH1D *hpuf;
  // control plots of components (JEC tag-and-probe)
  TProfile *pncandtp;
  TProfile *pnchtp;
  TProfile *pnnetp;
  TProfile *pnnhtp;
  TProfile *pncetp;
  TProfile *pnmutp;
  TProfile *pnhhtp;
  TProfile *pnhetp;
  TProfile *pchftp;
  TProfile *pneftp;
  TProfile *pnhftp;
  TProfile *pceftp;
  TProfile *pmuftp;
  TProfile *phhftp;
  TProfile *pheftp;
  TProfile *ppuftp;
  //
  TProfile *ppttagptprobe;
  TProfile *ppttageff;
  TProfile *ppttagmu;
  TProfile *ppttagptprobe_noa;
  TProfile *ppttagmu_noa;
  TH2D *h2pttagptprobe;
  TH2D *h2pttagmu;
  TH2D *h2pttagptprobe_noa;
  TH2D *h2pttagmu_noa;
  //
  TH1D *hncandtp;
  TH1D *hnchtp;
  TH1D *hnnetp;
  TH1D *hnnhtp;
  TH1D *hncetp;
  TH1D *hnmutp;
  TH1D *hnhhtp;
  TH1D *hnhetp;
  TH1D *hchftp;
  TH1D *hneftp;
  TH1D *hnhftp;
  TH1D *hceftp;
  TH1D *hmuftp;
  TH1D *hhhftp;
  TH1D *hheftp;
  TH1D *hpuftp;
  // Npv
  TProfile *pncandtp_vsnpv;
  TProfile *pnchtp_vsnpv;
  TProfile *pnnetp_vsnpv;
  TProfile *pnnhtp_vsnpv;
  TProfile *pncetp_vsnpv;
  TProfile *pnmutp_vsnpv;
  TProfile *pnhhtp_vsnpv;
  TProfile *pnhetp_vsnpv;
  TProfile *pchftp_vsnpv;
  TProfile *pneftp_vsnpv;
  TProfile *pnhftp_vsnpv;
  TProfile *pceftp_vsnpv;
  TProfile *pmuftp_vsnpv;
  TProfile *phhftp_vsnpv;
  TProfile *pheftp_vsnpv;
  TProfile *ppuftp_vsnpv;
  // Pileup
  TProfile *pchftp_vstrpu;
  TProfile *pneftp_vstrpu;
  TProfile *pnhftp_vstrpu;
  TProfile *pceftp_vstrpu;
  TProfile *pmuftp_vstrpu;
  TProfile *phhftp_vstrpu;
  TProfile *pheftp_vstrpu;
  TProfile *ppuftp_vstrpu;
  // control plots vs phi ("pos" and "neg" for positive and negative eta values)
  TProfile *pchfpostp_vsphi;
  TProfile *pnefpostp_vsphi;
  TProfile *pnhfpostp_vsphi;
  TProfile *pcefpostp_vsphi;
  TProfile *pmufpostp_vsphi;
  TProfile *phhfpostp_vsphi;
  TProfile *phefpostp_vsphi;
  TProfile *ppufpostp_vsphi;
  TProfile *pchfnegtp_vsphi;
  TProfile *pnefnegtp_vsphi;
  TProfile *pnhfnegtp_vsphi;
  TProfile *pcefnegtp_vsphi;
  TProfile *pmufnegtp_vsphi;
  TProfile *phhfnegtp_vsphi;
  TProfile *phefnegtp_vsphi;
  TProfile *ppufnegtp_vsphi;

  // control plots for topology (JEC)
  TH1D *hselpt;
  TH1D *hy;
  TH1D *hy2;
  TH1D *heta;
  TH1D *heta2;
  TH1D *hphi;
  TH1D *hdphi;
  TH1D *hdpt;
  TProfile *pdpt;
  TH1D *hjet;
  TH1D *hmet;
  TH1D *hmetphi;
  // control plots for vertex
  TH1D *hpvndof;
  TH1D *hpvx;
  TH1D *hpvy;
  TH1D *hpvz;
  TH1D *hpvr;
  TH1D *hpvrho;
  // closure plots for JEC
  TH1D *hmpf;
  TH1D *hmpf1;
  TH1D *hmpf2;
  TH1D *hmpfx;
  TH1D *hmpfy;
  TH1D *hmpfz;
  TProfile *pmpf;
  TProfile *pmpf1;
  TProfile *pmpf2;
  TProfile *pmpfx;
  TProfile *pmpfy;
  TProfile *pmpfz;

  // Control plots of resolutions:
  // dijet asymmetry binned in dijet pT and 3rd jet pT
  TH3D *hdjasymm;
  TH3D *hdjmpf;
  TH3D *hdjasymmtp;
  TH3D *hdjmpftp;

  TH1D *hr21;
  TH1D *hr31;
  TH1D *hr32;
  TProfile *pr21;
  TProfile *pr31;
  TProfile *pr32;
  TProfile *px21;
  TProfile *px31;
  TProfile *px32;

  TH1D *hyeta;
  TH1D *hyeta2;
  TH2D *hetaphi;

  // MC checks
  TH1D *hpt_jt30;
  TH1D *hpt_jt60;
  //TH1D *hpt_jt80;
  TH1D *hpt_jt110;
  //TH1D *hpt_jt150;
  TH1D *hpt_jt190;
  TH1D *hpt_jt240;
  TH1D *hpt_jt300;
  TH1D *hpt_jt370;
  //
  TH1D *hpt0_jt30;
  TH1D *hpt0_jt60;
  //TH1D *hpt0_jt80;
  TH1D *hpt0_jt110;
  //TH1D *hpt0_jt150;
  TH1D *hpt0_jt190;
  TH1D *hpt0_jt240;
  TH1D *hpt0_jt300;
  TH1D *hpt0_jt370;

  // unfolding studies (Mikael)
  TH2D *mT;
  TH2D *mTf;
  TH2D *mTuw;
  TH2D *mTfuw;
  TH1D *mx;
  TH1D *mxf;
  TH1D *mxuw;
  TH1D *mxfuw;
  TH1D *my;
  TH1D *myuw;
  TH1D *myf;
  TH1D *myfuw;

  TH1D *hpthat;
  TH1D *hpthatnlo;
  TH1D *hpt_noid_g;
  TH1D *hpt_nojetid_g;
  TH1D *hpt_noevtid_g;
  TH1D *hpt_r;
  TH1D *hpt_g;
  TH1D *hpt_gtw;
  TH1D *hpt_gg;
  TH1D *hpt_gg0;
  TH1D *hpt_g0;
  TH1D *hpt_g0tw;
  TH1D *hpt_g0_tmp;
  //
  TProfile *ppt_r;
  TProfile *ppt_g;

  TProfile *ppt_probepertag;

  // Response closure
  TProfile2D *p2rvsnpv;
  TH2D *h2r_r;
  TH2D *h2r_g;
  TProfile *p2r_r;
  TProfile *p2r_g;
  TProfile *p2r_ruw;
  TProfile *p2r_guw;

  //Rapidity closure
  TH2D *h2dy_r;
  TH2D *h2dy_g;
  TProfile *p2dy_r;
  TProfile *p2dy_g;
  TProfile *p2dy_ruw;
  TProfile *p2dy_guw;
  TProfile2D *pdy_r;
  TProfile2D *pdy_g;

  HistosBasic(TDirectory *dir, string trigname, double etamin = 0., double etamax = 2.0,
              double pttrg = 10., double ptmin = 10., double ptmax = 50., bool ismcdir = false);
  ~HistosBasic();

  void Write();

// private:
  TDirectory *dir;
};

#endif // __HistosBasic_h__
