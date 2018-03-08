// Purpose: Plot PF energy fractions vs pT/PU/NPV (Run II cleaned code)
// Author: mikko.voutilainen(at)nospam.cern.ch
// Created: Dec 11, 2015
#include "DrawFracs.h"

void Fracs::DrawFracs(unsigned mode) {
  setTDRStyle();

  // Bool flags for mode
  _vspt = false;
  _vspu = false;
  _vsnpv = false;
  _vseta = false;
  _vsphi = false;
  _sphi = "";
  if (mode==0)
    _vspt=true;
  else if (mode==1)
    _vspu=true;
  else if (mode==2)
    _vsnpv=true;
  else if (mode==3)
    _vseta=true;
  else if (mode==4) {
    _vsphi = true;
    _sphi = "pos";
  } else if (mode==5) {
    _vsphi = true;
    _sphi = "neg";
  } else {
    cout << "Problematic mode " << mode << endl;
    return;
  }

  // Build appropriate wide binning
  _x.clear();
  if (_vspt) {
    for (auto &ptwi : _jp_wptrange) _x.push_back(ptwi);
  } else if (_vspu) {
    for (int i = 0; i <= 25; ++i) _x.push_back(-0.5+2*i);
  } else if (_vsnpv) {
    for (int i = 0; i <= 25; ++i) _x.push_back(-0.5+2*i);
  } else if (_vseta) {
    for (auto &ewi : _jp_wetarange) _x.push_back(ewi);
  } else if (_vsphi) {
    for (auto &phidx : _jp_phirange) _x.push_back(phidx);
  }

  assert(mode<_modes.size()+1);

  setTDRStyle();

  _dmc->cd();
  TDirectory *dmc = gDirectory;
  _ddt->cd();
  TDirectory *ddt = gDirectory;
  if (mode==3) {
    dmc->cd("../FullEta_Reco"); dmc = gDirectory;
    ddt->cd("../FullEta_Reco"); ddt = gDirectory;
    if (_pertrg) {
      for (unsigned itrg = 0; itrg < _jp_notrigs; ++itrg) {
        auto ctrg = _jp_triggers[itrg];
        makeProfile(mode, dmc, ddt, ctrg);
      }
    } else {
      makeProfile(mode, dmc, ddt, "");
    }
  } else {
    if (_pertrg) {
      for (unsigned int ieta = 0; ieta != _etas.size(); ++ieta) {
        for (unsigned itrg = 0; itrg < _jp_notrigs; ++itrg) {
          auto ctrg = _jp_triggers[itrg];
          makeProfile(mode, dmc, ddt, ctrg, _etas[ieta].first, _etas[ieta].second);
        }
      }
    } else {
      for (unsigned int ieta = 0; ieta != _etas.size(); ++ieta)
        makeProfile(mode, dmc, ddt, "", _etas[ieta].first, _etas[ieta].second);
    }
  }
} // DrawFracs


void Fracs::makeProfile(unsigned mode, TDirectory *dmc, TDirectory *ddt, string trg, double eta1, double eta2) {
  // List of differences
  string taguniq = string(_vseta?"_vseta":"") + string(_vsnpv?"_vsNPV":"") + string(_vspu?"_vsTRPU":"") + string(_vsphi?"_vsphi":"") + _sphi
                 + string(_shiftJES ? "_shiftJES" : "") + (_pertrg?"_"+trg:"") + string(_vseta?"":Form("_%1.1f-%1.1f", eta1, eta2));
  map<string, TH1D*> mdf;

  THStack *hsdt = new THStack(Form("hsdt%s",taguniq.c_str()),"stacked histograms");
  THStack *hsmc = new THStack(Form("hsmc%s",taguniq.c_str()),"stacked histograms");
  THStack *hsdf = new THStack(Form("hsdf%s",taguniq.c_str()),"differences");

  TH1D *h = new TH1D(Form("h%s",taguniq.c_str()),Form(";p_{T%s} (GeV);PF energy fractions",_tp=="tp"?",tag":""),_x.size()-1,&_x[0]);
  TH1D *h2 = new TH1D(Form("h2%s",taguniq.c_str()),Form(";p_{T%s} (GeV);%s-%s (%%)",_tp=="tp"?",tag":"",_dt_type.c_str(),_mc_type.c_str()),_x.size()-1,&_x[0]);

  if (_vspt) {
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
  } else if (_vspu) {
    h->SetXTitle("#mu");
    h2->SetXTitle("#mu");
  } else if (_vsnpv) {
    h->SetXTitle("N_{PV,good}");
    h2->SetXTitle("N_{PV,good}");
  } else if (_vseta) {
    h->SetXTitle("Eta");
    h2->SetXTitle("Eta");
  }
  h->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h2->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h->SetMaximum(1-1e-5); // One should occasionally check that there is no overflow
  h->SetMinimum(0+1e-5);
  h2->SetMinimum(_h2min[mode]);
  h2->SetMaximum(_h2max[mode]);

  lumi_13TeV = _plot_title;
  TCanvas *c1 = tdrDiCanvas(Form("c1%s",taguniq.c_str()),h,h2,4,0);

  c1->cd(1);
  TLegend *leg = 0;
  if (_order==0) {
    if (_vseta)
      leg = tdrLeg(0.38,0.1,0.58,0.4);
    else if (_vspt)
      leg = tdrLeg(0.20,0.18,0.50,0.48);
    else
      leg = tdrLeg(0.20,0.18,0.50,0.48);
  } else if (_order==1) {
    if (_vseta)
      leg = tdrLeg(0.38,0.45,0.58,0.75);
    else if (_vspt)
      leg = tdrLeg(0.20,0.45,0.50,0.75);
    else
      leg = tdrLeg(0.20,0.5,0.50,0.8);
  } else if (_order==2) {
    if (_vseta)
      leg = tdrLeg(0.38,0.1,0.58,0.4);
    else if (_vspt)
      leg = tdrLeg(0.20,0.25,0.50,0.55);
    else
      leg = tdrLeg(0.20,0.18,0.50,0.48);
  } else if (_order==3) {
    if (_vseta)
      leg = tdrLeg(0.38,0.1,0.58,0.4);
    else if (_vspt)
      leg = tdrLeg(0.20,0.18,0.50,0.48);
    else
      leg = tdrLeg(0.20,0.18,0.50,0.48);
  }

  const char *dirname = Form("Eta_%1.1f-%1.1f%s",eta1,eta2,_pertrg?("/"+trg).c_str():"");

  TH1D *href = new TH1D("href","", _x.size()-1, &_x[0]);
  map<string,TH1D*> mcHistos;
  map<string,TH1D*> dtHistos;
  for (auto &frc : _fracs) {
    const char *hname0 = Form("p%s%s%s%s",frc.c_str(),_sphi.c_str(),_tp.c_str(),_modes[mode].c_str());
    if (!_vseta) { bool enterdtdir = ddt->cd(dirname); assert(enterdtdir); }
    else if (_pertrg) { bool enterdtdir = ddt->cd(trg.c_str()); assert(enterdtdir); }
    else { bool enterdtdir = ddt->cd(); assert(enterdtdir); }
    TProfile *pdt = dynamic_cast<TProfile*>(gDirectory->Get(hname0));
    if (!pdt) {
      cout << hname0 << " not found in (dt) " << gDirectory->GetName() << endl << flush;
      gDirectory->ls();
      assert(false);
    }
    pdt->SetName(Form("%s_dt",pdt->GetName()));
    if (_vspt or _vseta) {
      pdt->Sumw2(kFALSE);
      TProfile *pdtmp = dynamic_cast<TProfile*>(pdt->Rebin(_x.size()-1,Form("%s_rb",pdt->GetName()),&_x[0]));
      pdt->Delete(); pdt = 0;
      pdt = pdtmp;
    }
    dtHistos[frc] = pdt->ProjectionX(Form("%s_px",pdt->GetName()));
    pdt->Delete();

    if (!_vseta) { bool entermcdir = dmc->cd(dirname); assert(entermcdir); }
    else if (_pertrg) { bool entermcdir = dmc->cd(trg.c_str()); assert(entermcdir); }
    else { bool entermcdir = dmc->cd(); assert(entermcdir); }
    TProfile *pmc = dynamic_cast<TProfile*>(gDirectory->Get(hname0));
    if (!pmc) {
      cout << hname0 << " not found in (mc) " << gDirectory->GetName() << endl << flush;
      gDirectory->ls();
      assert(false);
    }
    pmc->SetName(Form("%s_mc",pmc->GetName()));
    if (_vspt or _vseta) {
      pmc->Sumw2(kFALSE);
      TProfile *pmcmp = dynamic_cast<TProfile*>(pmc->Rebin(_x.size()-1,Form("%s_rb",pmc->GetName()),&_x[0]));
      pmc->Delete(); pdt = 0;
      pmc = pmcmp;
    }
    mcHistos[frc] = pmc->ProjectionX(Form("%s_px",pmc->GetName()));
    pmc->Delete();
  } // for frc in _fracs
  href->Delete();

  vector<pair<TH1D*,string>> hdts;
  for (auto &frc : _fracs) {
    if (!_dobeta and frc=="beta") continue;
    else if (frc=="muf") continue;

    TH1D *hmc = dynamic_cast<TH1D*>(mcHistos[frc]->Clone(Form("%sc",mcHistos[frc]->GetName()))); assert(hmc);
    TH1D *hdt = dynamic_cast<TH1D*>(dtHistos[frc]->Clone(Form("%sc",dtHistos[frc]->GetName()))); assert(hmc);

    // Scale data by response to MC-equivalent energy fractions
    if (_shiftJES) {
      for (int i = 1; i != hdt->GetNbinsX()+1; ++i) {
        double jec = jesShift(hdt->GetBinCenter(i));
        hdt->SetBinContent(i, hdt->GetBinContent(i)/jec);
      } // for i
    } // _shiftJES

    if (frc=="cef") { // For cef, add muf
      TH1D *hmc2 = mcHistos["muf"]; assert(hmc2);
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)+hmc2->GetBinContent(i));  // cef -> cef + muf

      TH1D *hdt2 = dtHistos["muf"]; assert(hdt2);
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)+hdt2->GetBinContent(i));  // cef -> cef + muf
    } else if (frc=="betastar") { // For betastar, multiply by chf
      TH1D *hmc2 = mcHistos["chf"]; assert(hmc2);
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)*hmc2->GetBinContent(i));  // betastar -> chf * betastar

      TH1D *hdt2 = dtHistos["chf"]; assert(hdt2);
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)*hdt2->GetBinContent(i));  // betastar -> chf * betastar
    } else if (frc=="beta") { // For beta, multiply by chf
      TH1D *hmc2 = mcHistos["chf"]; assert(hmc2);
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)*hmc2->GetBinContent(i));  // beta -> chf * beta

      TH1D *hdt2 = dtHistos["chf"]; assert(hdt2);
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)*hdt2->GetBinContent(i));  // beta -> chf * beta
    } else if (frc=="chf") { // For chf, multiply by (1-beta-betastar)
      TH1D *hmc2 = mcHistos["beta"]; assert(hmc2);
      TH1D *hmc3 = mcHistos["betastar"]; assert(hmc3);
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i) * (1 - (_dobeta ? hmc2->GetBinContent(i) : 0) - hmc3->GetBinContent(i)));  // chf -> chf*(1-beta-betastar)

      TH1D *hdt2 = dtHistos["beta"]; assert(hdt2);
      TH1D *hdt3 = dtHistos["betastar"]; assert(hdt3);
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i) * (1 - (_dobeta ? hdt2->GetBinContent(i) : 0) - hdt3->GetBinContent(i)));  // chf -> chf*(1-beta-betastar)
    } // others we do not touch

    hmc->SetMarkerStyle(kNone);
    hmc->SetFillStyle(1001);
    hmc->SetFillColor(_style[frc].first - 7);
    hmc->SetLineColor(_style[frc].first + 1);
    hsmc->Add(hmc, "SAME H");

    hdt->SetFillStyle(1001); // for legend
    hdt->SetFillColor(_style[frc].first - 7); // for legend
    hdt->SetLineColor(_style[frc].first + 1);
    hdt->SetMarkerStyle(_style[frc].second);
    hdt->SetMarkerSize(frc=="nhf" or (frc=="chf" and _dobeta) ? 1.3 : 1.0);
    if (_vspt or _vseta) hdt->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    hsdt->Add(hdt, "SAME P");

    // Then, do the difference
    TH1D *hdf = dynamic_cast<TH1D*>(hdt->Clone(Form("hdf%s",taguniq.c_str()))); assert(hdt);
    hdf->Add(hdt, hmc, 100, -100);
    hdf->SetLineColor(_style[frc].first + 1);
    hdf->SetMarkerColor(_style[frc].first + 1);
    hdf->SetMarkerStyle(_style[frc].second);
    hdf->SetLineWidth(2);
    hdf->SetMarkerSize(frc=="nhf" or (frc=="chf" and _dobeta) ? 1.8 : 1.3);
    hsdf->Add(hdf, "SAME P");

    c1->cd(2);

    hdf->Draw("SAME");
    mdf[frc] = hdf;
    hdts.push_back(std::make_pair(hdt,_name[frc]));
  } // for frc
  TLine *l = new TLine();
  l->DrawLine(_rangemin[mode], 0,_rangemax[mode], 0);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(h2->GetYaxis()->GetLabelSize());
  tex->DrawLatex(0.17,0.80,Form("Anti-k_{T} R=0.4%s%s",_shiftJES ?", shifted by JES":"",_pertrg?Form(" Trg=%s",trg.c_str()):""));

  for (auto rhdti = hdts.rbegin(); rhdti != hdts.rend(); ++rhdti)
    leg->AddEntry(rhdti->first, rhdti->second.c_str(),"PF");

  c1->cd(1);
  if (_vspt) gPad->SetLogx();
  hsmc->Draw("SAME");
  hsdt->Draw("SAME");
  leg->Draw("SAME"); // redraw
  if (_pertrg) {
    auto trigidx = std::find(_jp_triggers,_jp_triggers+_jp_notrigs,trg)-_jp_triggers;
    tex->SetTextSize(h2->GetYaxis()->GetLabelSize()/3.0);
    if (_order==0) {
      if (_vspt)
        tex->DrawLatex(0.2,0.65,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",",trg",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
      else
        tex->DrawLatex(0.38,0.65,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",_tp=="tp"?",tag":"",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
    } else if (_order==1) {
      if (_vspt)
        tex->DrawLatex(0.2,0.1,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",",trg",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
      else
        tex->DrawLatex(0.38,0.1,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",_tp=="tp"?",tag":"",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
    } else if (_order==2) {
      if (_vspt)
        tex->DrawLatex(0.2,0.80,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",",trg",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
      else
        tex->DrawLatex(0.38,0.65,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",_tp=="tp"?",tag":"",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
    } else if (_order==3) {
      if (_vspt)
        tex->DrawLatex(0.2,0.75,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",",trg",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
      else
        tex->DrawLatex(0.38,0.75,Form("p_{T%s} #in [%1.1f,%1.1f] GeV",_tp=="tp"?",tag":"",_jp_trigranges[trigidx][0],_jp_trigranges[trigidx][1]));
    }
  }
  gPad->RedrawAxis();

  c1->cd(2);
  if (_vspt) gPad->SetLogx();
  gPad->RedrawAxis();

  c1->SaveAs(Form("%s/DrawFracs%s.pdf",_savedir.c_str(), taguniq.c_str()));

  if (_dofit) { // Estimate jet response slope by analyzing composition
    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(h2->GetYaxis()->GetLabelSize()*0.7);
    if (_dobeta) {
      TF1 *fchf = new TF1("fchf",[this](Double_t *x, Double_t *p){return jesFit(x,p,_fhb);},40,3000,2);
      fchf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
      TH1D *hchfa = dynamic_cast<TH1D*>(mdf["beta"]->Clone("hcfha")); assert(hchfa);
      hchfa->Add(mdf["chf"]);
      hchfa->Add(mdf["betastar"]);
      hchfa->Fit(fchf,"QRN");
      fchf->SetLineColor(kRed+3);
      fchf->Draw("SAME");
      hchfa->SetMarkerStyle(kFullStar);
      hchfa->SetMarkerColor(kRed+3);
      hchfa->SetLineColor(kRed+3);
      hchfa->Draw("SAMEP");
    } else {
      TF1 *fchf = new TF1("fchf",[this](Double_t *x, Double_t *p){return jesFit(x,p,_fhb);},40,3000,2);
      fchf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
      mdf["chf"]->Fit(fchf,"QRN");
      fchf->SetLineColor(kRed+2);
      fchf->Draw("SAME");

      tex->SetTextColor(kRed);
      tex->DrawLatex(0.17,0.40,Form("%1.2f#pm%1.2f%%, #chi^2/NDF=%1.1f/%d", fchf->GetParameter(0),
                                    fchf->GetParError(0), fchf->GetChisquare(), fchf->GetNDF()));
    }

    TF1 *fnhf = new TF1("fnhf",[this](Double_t *x, Double_t *p){return jesFit(x,p,_fhb);},40,3000,2);
    fnhf->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
    mdf["nhf"]->Fit(fnhf,"QRN");
    fnhf->SetLineColor(kGreen+2);
    fnhf->Draw("SAME");

    TF1 *fnef = new TF1("fnef",[this](Double_t *x, Double_t *p){return jesFit(x,p,_fhb);},40,3000,2);
    fnef->SetParameters( 0.9881-1, 0.2440 ); // Fall15_25nsV2
    mdf["nef"]->Fit(fnef,"QRN");
    fnef->SetLineColor(kBlue+1);
    fnef->Draw("SAME");

    TH1D *hall = dynamic_cast<TH1D*>(mdf["chf"]->Clone("hall")); assert(hall);
    if (_dobeta) hall->Add(mdf["beta"]);
    hall->Add(mdf["betastar"]);
    hall->Add(mdf["nef"]);
    hall->Add(mdf["nhf"]);
    hall->Add(mdf["cef"]);
    hall->SetMarkerStyle(kFullCross);
    hall->SetMarkerColor(kBlack);

    h2->SetMaximum(+5);//+3.0);
    h2->SetMinimum(-5);//-1.5);
    c1->SaveAs(Form("%s/DrawFracs_WithFit%s.pdf",savedir.c_str(),taguniq.c_str()));
  }
}
