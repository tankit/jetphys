// Draw quark/gluon likelihood (QGL) vs pTgen
// Fit 2D shape for quarks, gluons, unmatched

void normalizeTH2D(TH2D* h2) {

  for (int ix = 1; ix != h2->GetNbinsX()+1; ++ix) {
    TH1D *h1 = h2->ProjectionY("_py",ix,ix);
    double norm = h1->Integral();
    for (int iy = 1; iy != h2->GetNbinsY()+1; ++iy) {
      h2->SetBinContent(ix, iy, h2->GetBinContent(ix, iy) / norm);
    } // for iy
    delete h1;
  } // for ix

} // void normalizeTH2D


void drawQGL() {

  TFile *f = new TFile("rootfiles/Pythia_qginclusive.root","READ");
  assert(f && !f->IsZombie());

  TH2D *h2q = (TH2D*)f->Get("hpTGenvsQGL_quark_0"); assert(h2q);
  normalizeTH2D(h2q);
  TH2D *h2g = (TH2D*)f->Get("hpTGenvsQGL_gluon_0"); assert(h2g);
  normalizeTH2D(h2g);

  // Draw 2D histogram
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gPad->SetLogx();
  gStyle->SetOptStat(kFALSE);

  h2q->Draw("COLZ");
  h2q->GetXaxis()->SetRangeUser(64.,3103.);
  h2q->GetYaxis()->SetRangeUser(0.,1.);
  h2q->GetZaxis()->SetRangeUser(0.,0.3);

  // Draw (doubly) normalized 1D projection
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  TH1D *hpt = h2q->ProjectionX("hpt");

  TH1D *h1q = h2q->ProjectionY("h1q",hpt->FindBin(64.),
			       hpt->FindBin(2000.));
  h1q->Scale(1./h1q->Integral());
  TH1D *h1q_lo = h2q->ProjectionY("h1q_lo",hpt->FindBin(64.),
				   hpt->FindBin(200.));
  h1q_lo->Scale(1./h1q_lo->Integral());
  TH1D *h1q_hi = h2q->ProjectionY("h1q_hi",hpt->FindBin(1000.),
				   hpt->FindBin(2000.));
  h1q_hi->Scale(1./h1q_hi->Integral());

  TH1D *h1g = h2g->ProjectionY("h1g",hpt->FindBin(64.),
			       hpt->FindBin(2000.));
  h1g->Scale(1./h1g->Integral());
  TH1D *h1g_lo = h2g->ProjectionY("h1g_lo",hpt->FindBin(64.),
				   hpt->FindBin(200.));
  h1g_lo->Scale(1./h1g_lo->Integral());
  TH1D *h1g_hi = h2g->ProjectionY("h1g_hi",hpt->FindBin(1000.),
				   hpt->FindBin(2000.));
  h1g_hi->Scale(1./h1g_hi->Integral());


  h1q->SetFillStyle(3445);
  h1q->SetFillColor(kBlue);
  h1q->Draw("HIST");
  h1q->GetYaxis()->SetRangeUser(0,0.3);
  h1q_lo->SetLineColor(kBlue+2);
  h1q_lo->Draw("SAME HIST");
  h1q_hi->SetLineColor(kBlue-9);
  h1q_hi->Draw("SAME HIST");
  h1q->SetLineColor(kBlue);
  h1q->SetLineWidth(2);
  h1q->Draw("SAME");

  h1g->SetFillStyle(3454);
  h1g->SetFillColor(kRed);
  h1g->SetLineColor(kRed);
  h1g->Draw("SAME HIST");
  h1g_lo->SetLineColor(kRed+2);
  h1g_lo->Draw("SAME HIST");
  h1g_hi->SetLineColor(kRed-9);
  h1g_hi->Draw("SAME HIST");

  //TF1 *fg = new TF1("fg","[0]*exp([1]*x)+[2]",0,0.5);
  //fg->SetParameters(0.20,-20.,0.01);
  //TF1 *fg = new TF1("fg","[0]*exp([1]/(1-x))+[2]",0,0.5);
  //fg->SetParameters(0.50,-20.,0.01);
  //TF1 *fg = new TF1("fq","[0]*exp([1]*x)+[2]*exp([3]*x)+[4]",0.,0.5);
  //fg->SetParameters(0.16,-20.,0.04,-10.,0.01);
  //TF1 *fg = new TF1("fg","[0]*exp([1]*x)+[2]*exp([3]*x) + [4] +"
  //	    "[5]*exp([6]*(1-x))+[7]*exp([8]*(1-x))",0.,1.0);
  TF1 *fg = new TF1("fg","[0]*exp([1]*sqrt(x*x+[9]*[9]))"
		    "+[2]*exp([3]*sqrt(x*x+[9]*[9])) + [4] +"
		    "[5]*exp([6]*(1-x))+[7]*exp([8]*(1-x))",0.,1.0);
  fg->SetParameters(0.16,-20.,0.04,-10., 0.01, 0.08,-40,0.02,-20, 0.001);
  fg->SetLineColor(kOrange+2);
  h1g->Fit(fg,"QRN");
  fg->Draw("SAME");

  //TF1 *fq = new TF1("fq","[0]*exp([1]*(1-x))+[2]",0.5,1.0);
  //fq->SetParameters(0.20,-20.,0.01);
  //TF1 *fq = new TF1("fq","[0]*pow(x,[1])+[2]",0.5,1.0);
  //fq->SetParameters(0.20,2.,0.01);
  //TF1 *fq = new TF1("fq","[0]*exp([1]*(1-x))+[2]*exp([3]*(1-x))+[4]",0.5,1.0);
  //fq->SetParameters(0.40,-40.,0.10,-20.,0.01);
  //TF1 *fq = new TF1("fq","[0]*exp([1]*x)+[2]*exp([3]*x) + [4] +"
  //	    "[5]*exp([6]*(1-x))+[7]*exp([8]*(1-x))",0.,1.0);
  TF1 *fq = new TF1("fq","[0]*exp([1]*sqrt(x*x+[9]*[9]))"
		    "+[2]*exp([3]*sqrt(x*x+[9]*[9])) + [4] +"
		    "[5]*exp([6]*(1-x))+[7]*exp([8]*(1-x))",0.,1.0);
  fq->SetParameters(0.08,-20.,0.02,-10., 0.01, 0.40,-40,0.10,-20, 0.001);
  fq->SetLineColor(kMagenta+2);
  h1q->Fit(fq,"QRN");
  fq->Draw("SAME");

  c2->SaveAs("pdf/drawQGL_1D.pdf");
}

