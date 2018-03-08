// Purpose:  Keep track of trigger rates and prescales per run
//
// Author:   mikko.voutilainen@cern.ch
// Created:  June 4, 2010
// Updated:  August 9, 2011
#include "HistosRun.h"

HistosRun::HistosRun(TDirectory *dir, double etamin, double etamax)
  : lumsum(0) {

  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  // phase space
  this->etamin = etamin;
  this->etamax = etamax;

  this->trg.push_back("jt40");
  this->trg.push_back("jt60");
  this->trg.push_back("jt80");
  this->trg.push_back("jt140");
  this->trg.push_back("jt200");
  this->trg.push_back("jt260");
  this->trg.push_back("jt320");
  this->trg.push_back("jt400");
  this->trg.push_back("jt450");
  //this->trg.push_back("jt500");

  this->pt["jt40"] = 49;//64;
  this->pt["jt60"] = 74;//97;
  this->pt["jt80"] = 97;//114;
  this->pt["jt140"] = 174;//196;
  this->pt["jt200"] = 300;//272;
  this->pt["jt260"] = 362;//330;
  this->pt["jt320"] = 430;//395;
  this->pt["jt400"] = 507;//548;
 this->pt["jt450"] = 548;//638;
 //this->pt["jt500"] = 737;

  curdir->cd();
}

HistosRun::~HistosRun() {

  dir->cd();

  int nruns = lums.size();

  // Runs
  TH1I *runs = new TH1I("runs",Form("runs %1.3g pb-1",lumsum),nruns,-0.5,nruns-0.5);
  TH1D *runlumi = new TH1D("runlumi",Form("runlumi %1.3g pb-1",lumsum),nruns,-0.5,nruns-0.5);
  TH1D *runlumi2 = new TH1D("runlumi2",Form("runlumi2 %1.3g pb-1",lumsum2),nruns,-0.5,nruns-0.5);

  { int i = 1;
    for (map<int, float>::const_iterator it = runlums.begin();
	 it != runlums.end(); ++it, ++i) {
      runlumi->SetBinContent(i, it->second);
      runlumi->GetXaxis()->SetBinLabel(i, Form("%d", it->first));
    } // for it
  }
  { int i = 1;
    for (map<int, float>::const_iterator it = runlums2.begin();
	 it != runlums2.end(); ++it, ++i) {
      runlumi2->SetBinContent(i, it->second);
      runlumi2->GetXaxis()->SetBinLabel(i, Form("%d", it->first));
    } // for it
  }

  map<string, TH1D*> npv_trg;
  map<string, TH1D*> npvgood_trg;
  map<string, TH1D*> c_chf;
  map<string, TH1D*> c_nef;
  map<string, TH1D*> c_nhf;
  map<string, TH1D*> c_betastar;
  map<string, TH1D*> c_betaprime;
  map<string, TH1D*> c_chftp;
  map<string, TH1D*> c_neftp;
  map<string, TH1D*> c_nhftp;
  map<string, TH1D*> c_betastartp;
  map<string, TH1D*> c_betaprimetp;
  map<string, TH1D*> h0_trg;
  map<string, TH1D*> h_trg;
  map<string, TH1D*> hw_trg; // prescale weights
  map<string, TH1D*> r0_trg;
  map<string, TH1D*> r_trg;
  map<string, TH1D*> rw_trg; // prescale weight
  map<string, TH1D*> r_ctrg; // calo
  map<string, TH1D*> pl_trg;
  map<string, TH1D*> pt_trg;

  for (unsigned int j = 0; j != trg.size(); ++j) {

    string const& t1 = trg[j];
    const char* t = t1.c_str();
    npv_trg[t] = new TH1D(Form("npv_%s",t),Form("npv_%s",t),nruns,-0.5,nruns-0.5);
    npvgood_trg[t] = new TH1D(Form("npvgood_%s",t),Form("npvgood_%s",t),nruns,-0.5,nruns-0.5);
    c_chf[t] = new TH1D(Form("c_chf_%s",t),Form("c_chf_%s",t),nruns,-0.5,nruns-0.5);
    c_nef[t] = new TH1D(Form("c_nef_%s",t),Form("c_nef_%s",t),nruns,-0.5,nruns-0.5);
    c_nhf[t] = new TH1D(Form("c_nhf_%s",t),Form("c_nhf_%s",t),nruns,-0.5,nruns-0.5);
    c_betastar[t] = new TH1D(Form("c_betastar_%s",t),Form("c_betastar_%s",t),nruns,-0.5,nruns-0.5);
    c_betaprime[t] = new TH1D(Form("c_betaprime_%s",t),Form("c_betaprime_%s",t),nruns,-0.5,nruns-0.5);
    c_chftp[t] = new TH1D(Form("c_chftp_%s",t),Form("c_chftp_%s",t),nruns,-0.5,nruns-0.5);
    c_neftp[t] = new TH1D(Form("c_neftp_%s",t),Form("c_neftp_%s",t),nruns,-0.5,nruns-0.5);
    c_nhftp[t] = new TH1D(Form("c_nhftp_%s",t),Form("c_nhftp_%s",t),nruns,-0.5,nruns-0.5);
    c_betastartp[t] = new TH1D(Form("c_betastartp_%s",t),Form("c_betastartp_%s",t),nruns,-0.5,nruns-0.5);
    c_betaprimetp[t] = new TH1D(Form("c_betaprimetp_%s",t),Form("c_betaprimetp_%s",t),nruns,-0.5,nruns-0.5);
    h0_trg[t] = new TH1D(Form("h0_%s",t),Form("h0_%s",t),nruns,-0.5,nruns-0.5);
    h_trg[t] = new TH1D(Form("h_%s",t),Form("h_%s",t),nruns,-0.5,nruns-0.5);
    hw_trg[t] = new TH1D(Form("hw_%s",t),Form("hw_%s",t),nruns,-0.5,nruns-0.5);
    r0_trg[t] = new TH1D(Form("r0_%s",t),Form("r0_%s",t),nruns,-0.5,nruns-0.5);
    r_trg[t] = new TH1D(Form("r_%s",t),Form("r_%s",t),nruns,-0.5,nruns-0.5);
    rw_trg[t] = new TH1D(Form("rw_%s",t),Form("rw_%s",t),nruns,-0.5,nruns-0.5);
    r_ctrg[t] = new TH1D(Form("r_c%s",t),Form("r_c%s",t),nruns,-0.5,nruns-0.5);
    pl_trg[t] = new TH1D(Form("pl_%s",t),Form("pl_%s",t),nruns,-0.5,nruns-0.5);
    pt_trg[t] = new TH1D(Form("pt_%s",t),Form("pt_%s",t),nruns,-0.5,nruns-0.5);

    int i = 1;
    for (map<int, float>::const_iterator it = runlums.begin(); it != runlums.end(); ++it, ++i) {

      int run = it->first;
      runs->SetBinContent(i, run);

      double ntrig = t_trg[t][run];
      double npv = this->npv_trg[t][run];
      npv_trg[t]->SetBinContent(i, ntrig ? npv / ntrig : 0.);
      npv_trg[t]->SetBinError(i, ntrig ? npv / pow(ntrig,1.5) : 0.);

      double npvgood = this->npvgood_trg[t][run];
      npvgood_trg[t]->SetBinContent(i, ntrig ? npvgood / ntrig : 0.);
      npvgood_trg[t]->SetBinError(i, ntrig ? npvgood / pow(ntrig,1.5) : 0.);

      double chf = this->c_chf[t][run];
      c_chf[t]->SetBinContent(i, ntrig ? chf / ntrig : 0.);
      c_chf[t]->SetBinError(i, ntrig ? chf / pow(ntrig,1.5) : 0.);
      //
      double nef = this->c_nef[t][run];
      c_nef[t]->SetBinContent(i, ntrig ? nef / ntrig : 0.);
      c_nef[t]->SetBinError(i, ntrig ? nef / pow(ntrig,1.5) : 0.);
      //
      double nhf = this->c_nhf[t][run];
      c_nhf[t]->SetBinContent(i, ntrig ? nhf / ntrig : 0.);
      c_nhf[t]->SetBinError(i, ntrig ? nhf / pow(ntrig,1.5) : 0.);
      //
      double betastar = this->c_betastar[t][run];
      c_betastar[t]->SetBinContent(i, ntrig ? betastar / ntrig : 0.);
      c_betastar[t]->SetBinError(i, ntrig ? betastar / pow(ntrig,1.5) : 0.);
      //
      double betaprime = this->c_betaprime[t][run];
      c_betaprime[t]->SetBinContent(i, ntrig ? betaprime / ntrig : 0.);
      c_betaprime[t]->SetBinError(i, ntrig ? betaprime / pow(ntrig,1.5) : 0.);

      h0_trg[t]->SetBinContent(i, p_trg[t][run]);
      h_trg[t]->SetBinContent(i, t_trg[t][run]);
      hw_trg[t]->SetBinContent(i, tw_trg[t][run]);
      //
      // Same for TP method
      double ntrigtp = t_trgtp[t][run];
      double chftp = this->c_chftp[t][run];
      c_chftp[t]->SetBinContent(i, ntrigtp ? chftp / ntrigtp : 0.);
      c_chftp[t]->SetBinError(i, ntrigtp ? chftp / pow(ntrigtp,1.5) : 0.);
      //
      double neftp = this->c_neftp[t][run];
      c_neftp[t]->SetBinContent(i, ntrigtp ? neftp / ntrigtp : 0.);
      c_neftp[t]->SetBinError(i, ntrigtp ? neftp / pow(ntrigtp,1.5) : 0.);
      //
      double nhftp = this->c_nhftp[t][run];
      c_nhftp[t]->SetBinContent(i, ntrigtp ? nhftp / ntrigtp : 0.);
      c_nhftp[t]->SetBinError(i, ntrigtp ? nhftp / pow(ntrigtp,1.5) : 0.);
      //
      double betastartp = this->c_betastartp[t][run];
      c_betastartp[t]->SetBinContent(i, ntrigtp ? betastartp / ntrigtp : 0.);
      c_betastartp[t]->SetBinError(i, ntrigtp ? betastartp / pow(ntrigtp,1.5) : 0.);
      //
      double betaprimetp = this->c_betaprimetp[t][run];
      c_betaprimetp[t]->SetBinContent(i, ntrigtp ? betaprimetp / ntrigtp : 0.);
      c_betaprimetp[t]->SetBinError(i, ntrigtp ? betaprimetp / pow(ntrigtp,1.5) : 0.);
      //
      double lum = runlums_trg[t][run];
      r0_trg[t]->SetBinContent(i, lum ? p_trg[t][run] / lum : 0.);
      r0_trg[t]->SetBinError(i, lum ? sqrt(p_trg[t][run]) / lum : 0.);
      r_trg[t]->SetBinContent(i, lum ? t_trg[t][run] / lum : 0.);
      r_trg[t]->SetBinError(i, lum ? sqrt(t_trg[t][run]) / lum : 0.);

      double lum0 = runlums[run];
      rw_trg[t]->SetBinContent(i, lum0 ? tw_trg[t][run] / lum0 : 0.);
      rw_trg[t]->SetBinError(i, lum0 ? sqrt(tw_trg[t][run]) / lum0 : 0.);

      r_ctrg[t]->SetBinContent(i, lum ? c_trg[t][run] / lum : 0.);
      r_ctrg[t]->SetBinError(i, lum ? sqrt(c_trg[t][run]) / lum : 0.);

      pl_trg[t]->SetBinContent(i, lum ? runlums[run] / lum : 0.);

      int k = min(i+1, int(trg.size())-1);
      string const& t2 = trg[k];
      double p = (p_trg[t2][run] ? double(p_trgpair[t1+t2][run])/p_trg[t2][run]
		  : 0.);
      double err_p = (p ? sqrt(p*(1-p)/p_trg[t2][run]) : 1.);
      double prescale = (p ? 1./p : 0.);
      double err_prescale = pow(prescale,2) * err_p;
      pt_trg[t]->SetBinContent(i, prescale);
      pt_trg[t]->SetBinError(i, err_prescale);

      // Set bin labels
      npv_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      npvgood_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      h0_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      h_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      hw_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      r0_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      r_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      rw_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      r_ctrg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run)); // calo
      pl_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
      pt_trg[t]->GetXaxis()->SetBinLabel(i,Form("%d",run));
    } // for ir
  } // for j

  dir->Write();
};
