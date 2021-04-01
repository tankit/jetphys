void drawPU(TString ul="16_UL_new", int hlt=450, bool useRMS = 0)
{

  TString pu_path = "pileup/DT/"+ul+"/";

  std::vector<string> pu_run;
  if(ul.Contains("16_UL")) pu_run = {"B", "C", "D", "E", "Fe", "Fl", "G", "H"};
  //if(ul.Contains("16_UL")) pu_run = {"B", "C", "D", "E", "Fe", "Fl", "G", "H", "BCD", "EF", "GH"};
  else if(ul.Contains("17_UL")) pu_run = {"B", "C", "D", "E", "F", "BCDEF", "DE"};
  else if(ul.Contains("18_UL")) pu_run = {"A", "B", "C", "D", "AB", "CD", "ABCD"};
  else { cout << "No pileup data! Exit! " << endl; return 0; }

  int hmin = 0, hmax = 60;
  if(ul.Contains("16_UL")) hmin = 16, hmax = 30;
  else if(ul.Contains("17_UL")) hmin = 16, hmax = 45;
  else if(ul.Contains("18_UL")) hmin = 16, hmax = 45;

  const int m = pu_run.size();

  double x[m], ex[m];
  double pu_mean[m], pu_error[m];

  TFile* f[m];
  TH1D* jetTrig[m];
  for(int i=0; i<m; ++i) {
    x[i] = i+0.5, ex[i] = 0;
    TString input = pu_path+pu_run[i]+"/pileup_DT.root"; 
    //cout << input << ": " << x[i] << endl;
    f[i] = new TFile(input);
    jetTrig[i] = (TH1D*)f[i]->Get(Form("jt%i",hlt)); jetTrig[i]->SetName(Form("jt%i_%i",hlt,i));
    pu_mean[i] = jetTrig[i]->GetMean(); pu_error[i] = useRMS ? pu_error[i] = jetTrig[i]->GetRMS() : jetTrig[i]->GetMeanError();
    printf(" Run %7s with %5.2f +/- %5.2f\n",pu_run[i].c_str(),pu_mean[i],pu_error[i]);
  }

  TH1F* hpileup = new TH1F("hpileup","",m,0,m);
  hpileup->SetMinimum(hmin); hpileup->SetMaximum(hmax);
  hpileup->SetStats(0);
  hpileup->GetXaxis()->SetTitle("Era");
  hpileup->GetYaxis()->SetTitle("Pileup");
  //hpileup->GetXaxis()->SetNdivisions(205);
  hpileup->GetXaxis()->SetLabelSize(0.05);
//  hpileup->GetYaxis()->SetTitleOffset(1.5);

  TGraphErrors *gr = new TGraphErrors(m,x,pu_mean,ex,pu_error);
  gr->SetHistogram(hpileup);
  gr->SetMarkerStyle(20); gr->SetMarkerSize(0.8);

  for(int i=0; i<m; ++i) {
    gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(x[i]), pu_run[i].c_str());
    hpileup->GetXaxis()->SetBinLabel(hpileup->GetXaxis()->FindBin(x[i]), pu_run[i].c_str());
  }
  gr->GetHistogram()->GetXaxis()->LabelsOption("h");

  TCanvas *cv = new TCanvas("cv","cv",0,0,600,600);
  gr->Draw("AP");

  gPad->SetTicks();

  TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.035); tex->SetTextColor(kBlack);
  if(hlt==0) tex->DrawLatex(0.14,0.91,"HLT_ZeroBias"); else tex->DrawLatex(0.14,0.91,Form("HLT_PFJet%i",hlt));

  if(ul.Contains("16_UL")) tex->DrawLatex(0.77,0.91,"UL2016");
  else if(ul.Contains("17_UL")) tex->DrawLatex(0.77,0.91,"UL2017");
  else if(ul.Contains("18_UL")) tex->DrawLatex(0.77,0.91,"UL2018");

  double leg_xmin = 0.12, leg_ymin = 0.75, leg_xmax = 0.30, leg_ymax = 0.85;  
  TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
  leg->SetBorderSize(0); leg->SetTextSize(0.03); leg->SetLineColor(1);
  leg->SetLineStyle(1); leg->SetLineWidth(1); leg->SetFillColor(0); leg->SetFillStyle(1001);
  if(hlt==0) leg->AddEntry(gr,"HLT_ZeroBias","P"); else leg->AddEntry(gr,Form("HLT_PFJet%i",hlt),"PL");
  //leg->Draw();

}
