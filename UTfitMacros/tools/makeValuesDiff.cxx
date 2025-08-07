#include "calcSignificance.cxx"

void makeValuesDiff(double mean1, double sigma1,
		    double mean2, double sigma2) {

  TH1D* diff = new TH1D("diff", "diff", 1000, -20.*sigma1, 20.*sigma1);
  for(int i=0; i<10000000; i++) 
    diff->Fill(gRandom->Gaus(mean1, sigma1)-gRandom->Gaus(mean2, sigma2));

  calcSign(diff, 0.);
}

