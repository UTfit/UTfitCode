void calcSign(TH1D* histo, double val) {
  double yval = histo->GetBinContent(histo->FindBin(val));
  double area = 0.;
  
  for (int i=1; i<histo->GetNbinsX(); i++) {
    double yval_i = histo->GetBinContent(i);
    if(yval_i > yval) area += yval_i;
  }

  area *= 1./histo->Integral();

  // look for the closest n sigmas
  double prob[100];
  for(int i=0; i<100; i++) {
    double x = (i+1.)/100.*10.;
    prob[i] = TMath::Erf(x/sqrt(2.));
  }

  int igood = -99;
  for(int i=0; i<100; i++) {
    if(prob[i] <= area && prob[i+1] > area) igood = i;
  }

  double mysigma = (igood+1.)/100.*10.;

  cout << "Bayes ratio = " <<  histo->GetMaximum()/yval << endl;
  cout << "Probability of the CI " << area << endl;
  cout << "corresponding to " << mysigma  << "sigmas" << endl;
}

void calcSignificance(TString filename, TString histoname, double val) {

  TFile* file = new TFile(filename);
  TH1D*  histo = (TH1D*) file->Get(histoname);

  calcSign(histo, val);
}

