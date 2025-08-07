#include "calcSignificance.cxx"

void makeHistoDiff(TString filename, TString histoname,
		   double mean, double sigma) {

  TFile* file = new TFile(filename);
  TH1D*  histo = (TH1D*) file->Get(histoname);

  TH1D* diff = new TH1D("diff", "diff", 1000, -20.*sigma, 20.*sigma);
  //  TH1D* diff = new TH1D("diff", "diff", 1000, -1.5, 1.5);

  for(int i=0; i<10000000; i++) 
    diff->Fill(histo->GetRandom()-gRandom->Gaus(mean, sigma));

  calcSign(diff, 0.);
}

