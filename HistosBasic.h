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
  TH1D *hnjet_pre;
  TH1D *hj0pt_pre;
  TH1D *hj0phi_pre;
  TH1D *hj1pt_pre;
  TH1D *hj2pt_pre;
  TH1D *hnjet0;
  TH1D *hj0pt0;
  TH1D *hj0phi0;
  TH1D *hj1pt0;
  TH1D *hj2pt0;
  TH1D *hrecoilpt0;
  TH1D *hrecoilphi0;
  TH1D *hrecoilrat0;
  TH1D *hrecoilrat20;
  TH1D *hrecoildphi0;
  TH1D *hrecoiljets0;
  TH1D *hrecoilmjb0;
  TH1D *hrecoilmjbinv0;
  TH1D *hrecoilmpftpinv0;
  TH1D *hrecoilmpftpinv20;
  TH1D *hrecoilmpf0;
  TH1D *hrecoilmpftp0;
  TH1D *hrecoilmpf20;
  TH1D *hrecoilmpftp20;
  TH1D *hleadmpf0;
  TH1D *hleadmpftp0;
  TH1D *hleadmpf20;
  TH1D *hleadmpftp20;
  TProfile *pmjb_recoil0;
  TProfile *pmjb_leading0;
  TProfile *pmjbinv_leading0;
  TProfile *pmpfinv_leading0;
  TProfile *pmpfinv_leading20;
  TProfile *pmpf_recoil0;
  TProfile *pmpf_recoil20;
  TProfile *pmpf_leading0;
  TProfile *pmpf_leading20;
  TProfile *pmpflead_leading0;
  TProfile *pmpflead_leading20;
  TProfile *pmpfleadinv_leading0;
  TProfile *pmpfleadinv_leading20;
  TProfile *pmjbinv_recoil0;
  TProfile *pmpflead_recoil0;
  TProfile *pmpflead_recoil20;
  TH2D *h2mjb_recoil0;
  TH2D *h2mjb_leading0;
  TH2D *h2mjbinv_leading0;
  TH2D *h2mpfinv_leading0;
  TH2D *h2mpfinv_leading20;
  TH2D *h2mpf_recoil0;
  TH2D *h2mpf_recoil20;
  TH2D *h2mpf_leading0;
  TH2D *h2mpf_leading20;
  TH2D *h2mpflead_leading0;
  TH2D *h2mpflead_leading20;
  TH2D *h2mpfleadinv_leading0;
  TH2D *h2mpfleadinv_leading20;
  TH2D *h2mjbinv_recoil0;
  TH2D *h2mpflead_recoil0;
  TH2D *h2mpflead_recoil20;

  TH1D *hnjet1;
  TH1D *hj0pt1;
  TH1D *hj0phi1;
  TH1D *hj1pt1;
  TH1D *hj2pt1;
  TH1D *hrecoilpt1;
  TH1D *hrecoilphi1;
  TH1D *hrecoilrat1;
  TH1D *hrecoilrat21;
  TH1D *hrecoildphi1;
  TH1D *hrecoiljets1;
  TH1D *hrecoilmjb1;
  TH1D *hrecoilmjbinv1;
  TH1D *hrecoilmpftpinv1;
  TH1D *hrecoilmpftpinv21;
  TH1D *hrecoilmpf1;
  TH1D *hrecoilmpftp1;
  TH1D *hrecoilmpf21;
  TH1D *hrecoilmpftp21;
  TH1D *hleadmpf1;
  TH1D *hleadmpftp1;
  TH1D *hleadmpf21;
  TH1D *hleadmpftp21;
  TProfile *pmjb_recoil1;
  TProfile *pmjb_leading1;
  TProfile *pmjbinv_leading1;
  TProfile *pmpfinv_leading1;
  TProfile *pmpfinv_leading21;
  TProfile *pmpf_recoil1;
  TProfile *pmpf_recoil21;
  TProfile *pmpf_leading1;
  TProfile *pmpf_leading21;
  TProfile *pmpflead_leading1;
  TProfile *pmpflead_leading21;
  TProfile *pmpfleadinv_leading1;
  TProfile *pmpfleadinv_leading21;
  TProfile *pmjbinv_recoil1;
  TProfile *pmpflead_recoil1;
  TProfile *pmpflead_recoil21;
  TH2D *h2mjb_recoil1;
  TH2D *h2mjb_leading1;
  TH2D *h2mjbinv_leading1;
  TH2D *h2mpfinv_leading1;
  TH2D *h2mpfinv_leading21;
  TH2D *h2mpf_recoil1;
  TH2D *h2mpf_recoil21;
  TH2D *h2mpf_leading1;
  TH2D *h2mpf_leading21;
  TH2D *h2mpflead_leading1;
  TH2D *h2mpflead_leading21;
  TH2D *h2mpfleadinv_leading1;
  TH2D *h2mpfleadinv_leading21;
  TH2D *h2mjbinv_recoil1;
  TH2D *h2mpflead_recoil1;
  TH2D *h2mpflead_recoil21;

  TH1D *hnjet;
  TH1D *hj0pt;
  TH1D *hj0phi;
  TH1D *hj1pt;
  TH1D *hj2pt;
  TH1D *hrecoilpt;
  TH1D *hrecoilphi;
  TH1D *hrecoilrat;
  TH1D *hrecoilrat2;
  TH1D *hrecoildphi;
  TH1D *hrecoiljets;
  TH1D *hrecoilmjb;
  TH1D *hrecoilmjb2;
  TH1D *hrecoilmjbinv;
  TH1D *hrecoilmpftpinv;
  TH1D *hrecoilmpftpinv2;
  TH1D *hrecoilmpf;
  TH1D *hrecoilmpftp;
  TH1D *hrecoilmpf2;
  TH1D *hrecoilmpftp2;
  TH1D *hleadmpf;
  TH1D *hleadmpftp;
  TH1D *hleadmpf2;
  TH1D *hleadmpftp2;
  TProfile *pmjb_recoil;
  TProfile *pmjb_leading;
  TProfile *pmjbinv_leading;
  TProfile *pmpfinv_leading;
  TProfile *pmpfinv_leading2;
  TProfile *pmpf_recoil;
  TProfile *pmpf_recoil2;
  TProfile *pmpf_leading;
  TProfile *pmpf_leading2;
  TProfile *pmpflead_leading;
  TProfile *pmpflead_leading2;
  TProfile *pmpfleadinv_leading;
  TProfile *pmpfleadinv_leading2;
  TProfile *pmjbinv_recoil;
  TProfile *pmpflead_recoil;
  TProfile *pmpflead_recoil2;
  TH2D *h2mjb_recoil;
  TH2D *h2mjb_leading;
  TH2D *h2mjbinv_leading;
  TH2D *h2mpfinv_leading;
  TH2D *h2mpfinv_leading2;
  TH2D *h2mpf_recoil;
  TH2D *h2mpf_recoil2;
  TH2D *h2mpf_leading;
  TH2D *h2mpf_leading2;
  TH2D *h2mpflead_leading;
  TH2D *h2mpflead_leading2;
  TH2D *h2mpfleadinv_leading;
  TH2D *h2mpfleadinv_leading2;
  TH2D *h2mjbinv_recoil;
  TH2D *h2mpflead_recoil;
  TH2D *h2mpflead_recoil2;
  TH1D *hdphi_recoiljet;
  TH1D *hdphi_recoiljetsel;
  TH1D *hweight_recoiljetsel;
  TH2D *h2etapt_recoiljetsel;
  TH2D *h2etapt_recoiljetsel_weight;
  TH1D *hj0genpt_pre;
  TH1D *hj0genr_pre;
  TH1D *hj0genpt0;
  TH1D *hj0genr0;
  TH1D *hj0genpt1;
  TH1D *hj0genr1;
  TH1D *hj0genpt;
  TH1D *hj0genr;
  TH1D *hrecoilgenpt0;
  TH1D *hrecoilgenphi0;
  TH1D *hrecoilgenpt1;
  TH1D *hrecoilgenphi1;
  TH1D *hrecoilgenpt;
  TH1D *hrecoilgenphi;
  TH1D *hresp_jet;
  TH1D *hresp_jetsel;
  TH1D *hresp_tmp;
  TH1D *hresp_leading0;
  TH1D *hresp_recoil0;
  TH1D *hresp_mjb0;
  TH1D *hresp_leading1;
  TH1D *hresp_recoil1;
  TH1D *hresp_mjb1;
  TH1D *hresp_leading;
  TH1D *hresp_recoil;
  TH1D *hresp_mjb;
  TProfile *presp_leading0;
  TProfile *presp_recoil0;
  TProfile *presp_leading1;
  TProfile *presp_recoil1;
  TProfile *presp_leading;
  TProfile *presp_recoil;
  TH1D *hjtptsf;
  TH1D *hjtptsmearNew;

  TProfile *pmjb_ptave0;
  TProfile *pmjbinv_ptave0;
  TProfile *pmpf_ptave0;
  TProfile *pmpf_ptave20;
  TProfile *pmpfinv_ptave0;
  TProfile *pmpfinv_ptave20;
  TProfile *pmjb_ptave1;
  TProfile *pmjbinv_ptave1;
  TProfile *pmpf_ptave1;
  TProfile *pmpf_ptave21;
  TProfile *pmpfinv_ptave1;
  TProfile *pmpfinv_ptave21;
  TProfile *pmjb_ptave;
  TProfile *pmjbinv_ptave;
  TProfile *pmpf_ptave;
  TProfile *pmpf_ptave2;
  TProfile *pmpfinv_ptave;
  TProfile *pmpfinv_ptave2;
  TH2D *h2mjb_ptave0;
  TH2D *h2mjbinv_ptave0;
  TH2D *h2mpf_ptave0;
  TH2D *h2mpf_ptave20;
  TH2D *h2mpfinv_ptave0;
  TH2D *h2mpfinv_ptave20;
  TH2D *h2mjb_ptave1;
  TH2D *h2mjbinv_ptave1;
  TH2D *h2mpf_ptave1;
  TH2D *h2mpf_ptave21;
  TH2D *h2mpfinv_ptave1;
  TH2D *h2mpfinv_ptave21;
  TH2D *h2mjb_ptave;
  TH2D *h2mjbinv_ptave;
  TH2D *h2mpf_ptave;
  TH2D *h2mpf_ptave2;
  TH2D *h2mpfinv_ptave;
  TH2D *h2mpfinv_ptave2;

  TH1D *hleadmpf_ue0;
  TH1D *hleadmpftp_ue0;
  TH2D *h2mpflead_leading_ue0;
  TH2D *h2mpfleadinv_leading_ue0;
  TProfile *pmpflead_leading_ue0;
  TProfile *pmpfleadinv_leading_ue0;
  TH1D *hleadmpf_ue1;
  TH1D *hleadmpftp_ue1;
  TH2D *h2mpflead_leading_ue1;
  TH2D *h2mpfleadinv_leading_ue1;
  TProfile *pmpflead_leading_ue1;
  TProfile *pmpfleadinv_leading_ue1;
  TH1D *hleadmpf_ue;
  TH1D *hleadmpftp_ue;
  TH2D *h2mpflead_leading_ue;
  TH2D *h2mpfleadinv_leading_ue;
  TProfile *pmpflead_leading_ue;
  TProfile *pmpfleadinv_leading_ue;

  TH1D *hdphi_jet120;
  TH1D *hdphi_jet230;
  TH1D *hdphi_jet121;
  TH1D *hdphi_jet231;
  TH1D *hdphi_jet12;
  TH1D *hdphi_jet23;
  TH1D *hdphi_recoilmin;
  TH1D *hdphi_recoilmax;
  TH2D *hdphi_recoilcol;
  TH1D *hCRecoil;
  TH1D *hCRecoil_ptcl;
  TProfile *pCRecoil;
  TProfile *pCRecoil_leading;
  TProfile *pCRecoil_ptave;
  TProfile *pCRecoil_ptcl;
  TProfile *pCRecoil_ptcl2;
  TProfile *pCRecoil_ptcl2_leading;
  TProfile *pCRecoil_ptcl2_ptave;
  TProfile *pjtgenptsf;
  TProfile *pjtptsf;
  TProfile *pjtgenptsf_recoil;
  TProfile *pjtptsf_recoil;
  TH1D *hCRecoil1;
  TH1D *hCRecoil_ptcl1;
  TProfile *pCRecoil1;
  TProfile *pCRecoil_leading1;
  TProfile *pCRecoil_ptave1;
  TProfile *pCRecoil_ptcl1;
  TProfile *pCRecoil_ptcl21;
  TProfile *pCRecoil_ptcl2_leading1;
  TProfile *pCRecoil_ptcl2_ptave1;
  TProfile *pjtgenptsf1;
  TProfile *pjtptsf1;
  TH1D *hCRecoil0;
  TH1D *hCRecoil_ptcl0;
  TProfile *pCRecoil0;
  TProfile *pCRecoil_leading0;
  TProfile *pCRecoil_ptave0;
  TProfile *pCRecoil_ptcl0;
  TProfile *pCRecoil_ptcl20;
  TProfile *pCRecoil_ptcl2_leading0;
  TProfile *pCRecoil_ptcl2_ptave0;
  TProfile *pjtgenptsf0;
  TProfile *pjtptsf0;

  TProfile *prho_recoil;
  TProfile *prho_leading;
  TProfile *prho_ptave;
  TProfile *prho_recoil0;
  TProfile *prho_leading0;
  TProfile *prho_ptave0;
  TProfile *prho_recoil1;
  TProfile *prho_leading1;
  TProfile *prho_ptave1;
  TProfile *prho_recoil2;
  TProfile *prho_leading2;
  TProfile *prho_ptave2;

  TH2D *h2rho_recoil;
  TH2D *h2rho_leading;
  TH2D *h2rho_ptave;
  TH2D *h2rho_recoil0;
  TH2D *h2rho_leading0;
  TH2D *h2rho_ptave0;
  TH2D *h2rho_recoil1;
  TH2D *h2rho_leading1;
  TH2D *h2rho_ptave1;
  TH2D *h2rho_recoil2;
  TH2D *h2rho_leading2;
  TH2D *h2rho_ptave2;

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
  TH1D *hmpfinv;
  TH1D *hmpfinv1;
  TH1D *hmpfinv2;
  TH1D *hmpfx;
  TH1D *hmpfy;
  TH1D *hmpfz;
  TProfile *pmpf;
  TProfile *pmpf1;
  TProfile *pmpf2;
  TProfile *pmpfinv;
  TProfile *pmpfinv1;
  TProfile *pmpfinv2;
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
