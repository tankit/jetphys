// Purpose: Plot PF energy fractions vs pT/PU/NPV (Run II cleaned code)
// Author: mikko.voutilainen(at)nospam.cern.ch
// Created: Dec 11, 2015
#include "DrawFracs.h"

void Fracs::DrawFracs(unsigned mode) {
  setTDRStyle();
  TH1::SetDefaultSumw2(kTRUE);

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
    for (auto &ptwi : jp::wptrange) _x.push_back(ptwi);
  } else if (_vspu) {
    for (int i = 0; i <= 25; ++i) _x.push_back(-0.5+2*i);
  } else if (_vsnpv) {
    for (int i = 0; i <= 25; ++i) _x.push_back(-0.5+2*i);
  } else if (_vseta) {
    for (auto &ewi : jp::wetarange) _x.push_back(ewi);
  } else if (_vsphi) {
    for (auto &phidx : jp::phirange) _x.push_back(phidx);
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
      for (unsigned itrg = 0; itrg < jp::notrigs; ++itrg) {
        auto ctrg = jp::triggers[itrg];
        makeProfile(mode, dmc, ddt, ctrg);
      }
    } else {
      makeProfile(mode, dmc, ddt, "");
    }
  } else {
    if (_pertrg) {
      for (unsigned int ieta = 0; ieta != _etas.size(); ++ieta) {
        for (unsigned itrg = 0; itrg < jp::notrigs; ++itrg) {
          auto ctrg = jp::triggers[itrg];
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
                 + string(_shiftJES ? "_shiftJES" : "") + (_pertrg?"_"+trg:"") + string(_vseta?"":Form("_%.0f-%.0f", 10*eta1, 10*eta2));
  map<string, TH1D*> mdf;

  THStack *hsdt = new THStack(Form("hsdt%s",taguniq.c_str()),"stacked histograms");
  THStack *hsmc = new THStack(Form("hsmc%s",taguniq.c_str()),"stacked histograms");
  THStack *hsdf = new THStack(Form("hsdf%s",taguniq.c_str()),"differences");

  TH1D *h = new TH1D(Form("h%s",taguniq.c_str()),Form(";p_{T}^{%s} (GeV);PF energy fractions",_tp=="tp"?"tag":""),_x.size()-1,&_x[0]);
  TH1D *h2 = new TH1D(Form("h2%s",taguniq.c_str()),Form(";p_{T}^{%s} (GeV);%s-%s #scale[1.15]{(}10^{-2}#scale[1.15]{)}",_tp=="tp"?"tag":"",_dt_typeverb.c_str(),_mc_type.c_str()),_x.size()-1,&_x[0]);

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
    h->SetXTitle("#eta");
    h2->SetXTitle("#eta");
  }
  h->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h2->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
  h->SetMaximum(1-1e-5); // One should occasionally check that there is no overflow
  h->SetMinimum(0+1e-5);
  h2->SetMinimum(_h2min[mode]);
  h2->SetMaximum(_h2max[mode]);

  lumi_13TeV = _lumi_title;
  lumi_bonus = _plot_title;
  TCanvas *c1 = tdrDiCanvas(Form("c1%s%s",taguniq.c_str(),_tp.c_str()),h,h2,4,33);

  c1->cd(1);

  double xLo = 0.0, xHi = 0.0, yLo = 0.0, yHi = 0.0;
  if (_order==0) {
    if (_vseta) {
      xLo = 0.38; yLo = 0.1; xHi = 0.58; yHi = 0.4;
    } else if (_vspt) {
      xLo = 0.20; yLo = 0.18; xHi = 0.50; yHi = 0.48;
    } else {
      xLo = 0.20; yLo = 0.18; xHi = 0.50; yHi = 0.48;
    }
  } else if (_order==1) {
    if (_vseta) {
      xLo = 0.38; yLo = 0.45; xHi = 0.58; yHi = 0.75;
    } else if (_vspt) {
      xLo = 0.20; yLo = 0.45; xHi = 0.50; yHi = 0.75;
    } else {
      xLo = 0.20; yLo = 0.5; xHi = 0.50; yHi = 0.8;
    }
  } else if (_order==2) {
    if (_vseta) {
      xLo = 0.38; yLo = 0.1; xHi = 0.58; yHi = 0.4;
    } else if (_vspt) {
      xLo = 0.20; yLo = 0.25; xHi = 0.50; yHi = 0.55;
    } else {
      xLo = 0.20; yLo = 0.18; xHi = 0.50; yHi = 0.48;
    }
  } else if (_order==3) {
    if (_vseta) {
      xLo = 0.36; yLo = 0.15; xHi = 0.56; yHi = 0.45;
    } else if (_vspt) {
      if (eta2 != 3.0) {
        xLo = 0.18; yLo = 0.18; xHi = 0.48; yHi = 0.48;
      } else {
        xLo = 0.18; yLo = 0.14; xHi = 0.48; yHi = 0.44;
      }
    } else {
      xLo = 0.20; yLo = 0.18; xHi = 0.50; yHi = 0.48;
    }
  }
  TLegend *leg = tdrLeg(xLo,yLo,xHi,yHi);
  double scalefac = min(gPad->XtoPixel(gPad->GetX2()),gPad->YtoPixel(gPad->GetY1()));
  leg->SetTextFont(63);
  leg->SetTextSize(leg->GetTextSize()*scalefac);

  const char *dirname = Form("Eta_%1.1f-%1.1f%s",eta1,eta2,_pertrg?("/"+trg).c_str():"");

  map<string,TH1D*> mcHistos;
  map<string,TH1D*> dtHistos;
  for (auto &frc : _fracs) {
    const char *hname0 = Form("p%s%s%s%s",frc.c_str(),_sphi.c_str(),_tp.c_str(),_modes[mode].c_str());
    if (!_vseta) {
      if (!_vspt and (frc=="hef" or frc=="hhf")) continue;
//      if (_tp=="" and (frc=="hef" or frc=="hhf")) continue;
      bool enterdtdir = ddt->cd(dirname); assert(enterdtdir);
    } else if (_pertrg) { bool enterdtdir = ddt->cd(trg.c_str()); assert(enterdtdir); }
    else { bool enterdtdir = ddt->cd(); assert(enterdtdir); }
    TProfile *pdt = dynamic_cast<TProfile*>(gDirectory->Get(hname0));
    if (!pdt) {
      cout << hname0 << " not found in (dt) " << gDirectory->GetName() << endl << flush;
      gDirectory->ls();
      assert(false);
    }
    pdt->SetName(Form("%s_dt",pdt->GetName()));
    if (_vspt or _vseta) {
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
      TProfile *pmcmp = dynamic_cast<TProfile*>(pmc->Rebin(_x.size()-1,Form("%s_rb",pmc->GetName()),&_x[0]));
      pmc->Delete(); pdt = 0;
      pmc = pmcmp;
    }
    mcHistos[frc] = pmc->ProjectionX(Form("%s_px",pmc->GetName()));
    pmc->Delete();
  } // for frc in _fracs

  vector<pair<TH1D*,string>> hdts;
  for (auto &frc : _fracs) {
    if (frc=="muf") continue;
    else if (!_vseta and (!_vspt or eta2<3.0)) {
        if (frc=="hef" or frc=="hhf") continue;
    }

    TH1D *hmc = dynamic_cast<TH1D*>(mcHistos[frc]->Clone(Form("%sc",mcHistos[frc]->GetName()))); assert(hmc);
    TH1D *hdt = dynamic_cast<TH1D*>(dtHistos[frc]->Clone(Form("%sc",dtHistos[frc]->GetName()))); assert(hmc);
    hmc->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    hdt->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    TH1D *hdf = dynamic_cast<TH1D*>(hdt->Clone(Form("hdf%s",taguniq.c_str()))); assert(hdf);

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
    } else if (frc=="beta") { // For beta, multiply by chf
      TH1D *hmc2 = mcHistos["chf"]; assert(hmc2);
      for (int i = 1; i != hmc2->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i)*hmc2->GetBinContent(i));  // beta -> chf * beta

      TH1D *hdt2 = dtHistos["chf"]; assert(hdt2);
      for (int i = 1; i != hdt2->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i)*hdt2->GetBinContent(i));  // beta -> chf * beta
    } else if (frc=="chf") { // For chf, multiply by (1-beta-betastar)
      hdf->Add(hdt, hmc, 100, -100);
      TH1D *hmc3 = mcHistos["puf"]; assert(hmc3);
      for (int i = 1; i != hmc->GetNbinsX()+1; ++i)
        hmc->SetBinContent(i, hmc->GetBinContent(i) - hmc3->GetBinContent(i));  // chf -> chf*(1-beta-betastar)

      TH1D *hdt3 = dtHistos["puf"]; assert(hdt3);
      for (int i = 1; i != hdt->GetNbinsX()+1; ++i)
        hdt->SetBinContent(i, hdt->GetBinContent(i) - hdt3->GetBinContent(i));  // chf -> chf*(1-beta-betastar)
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
    hdt->SetMarkerSize(frc=="nhf" ? 1.3 : 1.0);
    if (_vspt or _vseta) hdt->GetXaxis()->SetRangeUser(_rangemin[mode],_rangemax[mode]);
    hsdt->Add(hdt, "SAME P");

    // Then, do the difference
    if (frc!="chf") hdf->Add(hdt, hmc, 100, -100);
    hdf->SetLineColor(_style[frc].first + 1);
    hdf->SetMarkerColor(_style[frc].first + 1);
    hdf->SetMarkerStyle(_style[frc].second);
    hdf->SetLineWidth(2);
    hdf->SetMarkerSize(frc=="nhf" ? 1.8 : 1.3);
    hsdf->Add(hdf,"SAME P");

    c1->cd(2);

    hdf->Draw("SAME P");
    mdf[frc] = hdf;
    hdts.push_back(std::make_pair(hdt,_name[frc]));
  } // for frc
  TLine *l = new TLine();
  l->DrawLine(_rangemin[mode], 0,_rangemax[mode], 0);
  c1->cd(2);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(63);
  tex->SetTextSize(leg->GetTextSize());
  //tex->DrawLatex(0.15,0.306,Form("Anti-k_{T} R=0.4 PF+CHS%s%s",_shiftJES ?", shifted by JES":"",_pertrg?Form(" Trg=PFJet%s",trg.substr(2).c_str()):""));
  int trigidx = _pertrg ? int(std::find(jp::triggers.begin(),jp::triggers.end(),trg)-jp::triggers.begin()) : 0;
  if (_pertrg) tex->DrawLatex(0.164,0.33,Form("%s",Form("Trg=PFJet%s, %1.1f#leq p_{T}[GeV]<%1.1f",trg.substr(2).c_str(),jp::trigranges[trigidx][0],jp::trigranges[trigidx][1])));
  tex->DrawLatex(0.164,0.83,Form("Anti-k_{T} R=0.4 PF+CHS"));

  for (auto rhdti = hdts.rbegin(); rhdti != hdts.rend(); ++rhdti)
    leg->AddEntry(rhdti->first, rhdti->second.c_str(),"PF");

  c1->cd(1);
  if (_vspt) gPad->SetLogx();
  hsmc->Draw("SAME");
  hsdt->Draw("SAME");
  leg->Draw("SAME"); // redraw

  double yD = 0.06;
  double yP = 0.85;
  if (_vspt) {
    if (eta2 <= 2.5) {
      if (eta1 == 0.0) {
        yP = 0.83;
        tex->DrawLatex(xLo,yP,Form("|#eta^{%s}|<%1.1f",_tp=="tp"?"probe":"",eta2));
      } else
        tex->DrawLatex(xLo,yP,Form("%1.1f#leq |#eta^{%s}|<%1.1f",eta1,_tp=="tp"?"probe":"",eta2));
      tex->DrawLatex(xLo,yP-yD,Form("%s: markers",_dt_typeverb.c_str()));
      tex->DrawLatex(xLo,yP-2*yD,Form("%s: histograms",_mc_typeverb.c_str()));
    } else {
      yP = 0.75;
      tex->DrawLatex(xLo,yP,Form("%1.1f#leq |#eta^{%s}|<%1.1f",eta1,_tp=="tp"?"probe":"",eta2));
      tex->DrawLatex(xLo,yP-yD,Form("%s: markers",_dt_typeverb.c_str()));
      tex->DrawLatex(xLo,yP-2*yD,Form("%s: histograms",_mc_typeverb.c_str()));
    }
    //tex->DrawLatex(0.75,yP,"19.7 fb^{-1}");
  } else if (_vseta) {
    yP = 0.81;
    tex->DrawLatex(xLo,yP,Form("%s: markers",_dt_typeverb.c_str()));
    tex->DrawLatex(xLo,yP-yD,Form("%s: histograms",_mc_typeverb.c_str()));
  } else {
    yP = 0.81;
    tex->DrawLatex(xLo,yP,Form("%s: markers",_dt_typeverb.c_str()));
    tex->DrawLatex(xLo,yP-yD,Form("%s: histograms",_mc_typeverb.c_str()));
  }
  gPad->RedrawAxis();
  if (_usetriglines and _vspt and !_pertrg) {
    l->SetLineColor(kGray+2);
    for (auto i = 0u; i < jp::notrigs-1; ++i) {
      l->DrawLine(jp::trigranges[i][1],0.0,jp::trigranges[i][1],1.0);
    }
  }
  c1->cd(1);
  CMS_lumi( (TCanvas*)gPad, 4, 33 );

  c1->cd(2);
  if (_vspt) gPad->SetLogx();
  gPad->RedrawAxis();

  gErrorIgnoreLevel = kWarning;
  c1->SaveAs(Form("%s/drawFracs%s%s.pdf",_savedir.c_str(),taguniq.c_str(),_tp=="tp"?"":"_dm"));
  cout << Form("\\includegraphics[width=0.30\\textwidth]{../%s/drawFracs%s%s.pdf}",_savedir.c_str(),taguniq.c_str(),_tp=="tp"?"":"_dm") << endl;
}
