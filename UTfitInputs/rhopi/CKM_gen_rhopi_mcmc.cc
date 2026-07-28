#include <iostream>
#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TComplex.h>
#include <TDatime.h>
#include <TF1.h>
#include <string>
#include <vector>
#include <TString.h>
#include <TMatrixF.h>
#include <TMatrixD.h>
#include <TVector.h>
#include <TComplex.h>
#include "CKM_gen_rhopi_mcmc.hh"

using namespace std;

double Pi = acos(-1.);
bool lfirst = true;

// ---------------------------------------------------------------------
// Parallel-tempering sampling of the 5 free complex rho-pi amplitudes.
//
// The amplitude convention (Apm/Amp/Azz/Apm_bar/Amp_bar/Azz_bar, Apm fixed
// at (1,0) as the normalization/global-phase convention), the 26 U/I
// observable formulas, and the alpha-extraction formula are unchanged from
// the original generator and were verified against arXiv:hep-ex/0701015
// and arXiv:1304.3503. Only the sampling method changes: the original code
// drew i.i.d. flat-prior amplitude configurations and importance-weighted
// them by exp(-chi2/2) against the 26-parameter Belle+BaBar measurement.
// Because several of those 26 observables are tightly constrained, the
// high-likelihood region is a tiny fraction of the flat-prior volume:
// effective sample size (Kish ESS, via TH1::GetEffectiveEntries()) stayed
// ~O(1-10) even at tens of millions of draws, with the output histogram
// dominated by a handful of "lucky" draws and not reproducible across
// random seeds. Parallel tempering (replica exchange MCMC) fixes this.
// ---------------------------------------------------------------------

struct Walker {
  double magAmp, phAmp;
  double magAzz, phAzz;
  double magApmB, phApmB;
  double magAmpB, phAmpB;
  double magAzzB, phAzzB;
  double chisq;
};

// Fold x into [lo,hi] by reflection. Done with a single fmod rather than an
// iterative loop: a proposal far outside the range would otherwise cost one
// iteration per 2*(hi-lo) of excursion.
double ReflectBound(double x, double lo, double hi) {
  if (hi<=lo) return lo;
  if (!std::isfinite(x)) return 0.5*(lo+hi);
  const double range = hi-lo, period = 2.*range;
  double t = fmod(x-lo, period);
  if (t<0.) t += period;
  return (t<=range) ? lo+t : hi-(t-range);
}

// Wrap x into (-Pi,Pi], likewise in closed form.
double WrapPhase(double x) {
  if (!std::isfinite(x)) return 0.;
  double t = fmod(x+Pi, 2.*Pi);
  if (t<=0.) t += 2.*Pi;
  return t-Pi;
}

Walker RandomWalker(TRandom3 *rng) {
  Walker w;
  w.magAmp  = rng->Rndm()*Amp_max;  w.phAmp  = (2.*rng->Rndm()-1.)*Pi;
  w.magAzz  = rng->Rndm()*Azz_max;  w.phAzz  = (2.*rng->Rndm()-1.)*Pi;
  w.magApmB = rng->Rndm()*Apm_max;  w.phApmB = (2.*rng->Rndm()-1.)*Pi;
  w.magAmpB = rng->Rndm()*Amp_max;  w.phAmpB = (2.*rng->Rndm()-1.)*Pi;
  w.magAzzB = rng->Rndm()*Azz_max;  w.phAzzB = (2.*rng->Rndm()-1.)*Pi;
  w.chisq = 0.;
  return w;
}

// Builds the 6 rho-pi amplitudes from a walker, computes the 26 U/I
// observables (formulas unchanged from the original generator), and
// returns chi2 against the enabled Belle/BaBar measurements.
double ChiSqFromWalker(const Walker &w, const TMatrixDSym &Vneu_BaBar,
			const TMatrixDSym &Vneu_Belle, double *alphaOut) {

  TComplex Apm(1., 0., kTRUE);
  TComplex Amp(w.magAmp, w.phAmp, kTRUE);
  TComplex Azz(w.magAzz, w.phAzz, kTRUE);
  TComplex Apm_bar(w.magApmB, w.phApmB, kTRUE);
  TComplex Amp_bar(w.magAmpB, w.phAmpB, kTRUE);
  TComplex Azz_bar(w.magAzzB, w.phAzzB, kTRUE);

  TComplex e21alpha = (Apm_bar+Amp_bar+2.*Azz_bar)/(Apm+Amp+2.*Azz);
  double alpha = e21alpha.Theta()/2.;
  if (alpha<0.) alpha += Pi;
  if (alphaOut) *alphaOut = alpha;

  double Upm = pow(Apm.Rho(),2.) + pow(Apm_bar.Rho(),2.);
  double Ump = (pow(Amp.Rho(),2.) + pow(Amp_bar.Rho(),2.))/Upm;
  double Uzz = (pow(Azz.Rho(),2.) + pow(Azz_bar.Rho(),2.))/Upm;

  double Vpm = (pow(Apm.Rho(),2.) - pow(Apm_bar.Rho(),2.))/Upm;
  double Vmp = (pow(Amp.Rho(),2.) - pow(Amp_bar.Rho(),2.))/Upm;
  double Vzz = (pow(Azz.Rho(),2.) - pow(Azz_bar.Rho(),2.))/Upm;

  double Re_Upm = ((Apm.Re()*Amp.Re() + Apm.Im()*Amp.Im()) + (Apm_bar.Re()*Amp_bar.Re() + Apm_bar.Im()*Amp_bar.Im()))/Upm;
  double Re_Upz = ((Apm.Re()*Azz.Re() + Apm.Im()*Azz.Im()) + (Apm_bar.Re()*Azz_bar.Re() + Apm_bar.Im()*Azz_bar.Im()))/Upm;
  double Re_Umz = ((Amp.Re()*Azz.Re() + Amp.Im()*Azz.Im()) + (Amp_bar.Re()*Azz_bar.Re() + Amp_bar.Im()*Azz_bar.Im()))/Upm;

  double Im_Upm = ((Apm.Im()*Amp.Re() - Apm.Re()*Amp.Im()) + (Apm_bar.Im()*Amp_bar.Re() - Apm_bar.Re()*Amp_bar.Im()))/Upm;
  double Im_Upz = ((Apm.Im()*Azz.Re() - Apm.Re()*Azz.Im()) + (Apm_bar.Im()*Azz_bar.Re() - Apm_bar.Re()*Azz_bar.Im()))/Upm;
  double Im_Umz = ((Amp.Im()*Azz.Re() - Amp.Re()*Azz.Im()) + (Amp_bar.Im()*Azz_bar.Re() - Amp_bar.Re()*Azz_bar.Im()))/Upm;

  double Re_Vpm = ((Apm.Re()*Amp.Re() + Apm.Im()*Amp.Im()) - (Apm_bar.Re()*Amp_bar.Re() + Apm_bar.Im()*Amp_bar.Im()))/Upm;
  double Re_Vpz = ((Apm.Re()*Azz.Re() + Apm.Im()*Azz.Im()) - (Apm_bar.Re()*Azz_bar.Re() + Apm_bar.Im()*Azz_bar.Im()))/Upm;
  double Re_Vmz = ((Amp.Re()*Azz.Re() + Amp.Im()*Azz.Im()) - (Amp_bar.Re()*Azz_bar.Re() + Amp_bar.Im()*Azz_bar.Im()))/Upm;

  double Im_Vpm = ((Apm.Im()*Amp.Re() - Apm.Re()*Amp.Im()) - (Apm_bar.Im()*Amp_bar.Re() - Apm_bar.Re()*Amp_bar.Im()))/Upm;
  double Im_Vpz = ((Apm.Im()*Azz.Re() - Apm.Re()*Azz.Im()) - (Apm_bar.Im()*Azz_bar.Re() - Apm_bar.Re()*Azz_bar.Im()))/Upm;
  double Im_Vmz = ((Amp.Im()*Azz.Re() - Amp.Re()*Azz.Im()) - (Amp_bar.Im()*Azz_bar.Re() - Amp_bar.Re()*Azz_bar.Im()))/Upm;

  double Ipm = (Apm_bar.Im()*Apm.Re() - Apm_bar.Re()*Apm.Im())/Upm;
  double Imp = (Amp_bar.Im()*Amp.Re() - Amp_bar.Re()*Amp.Im())/Upm;
  double Izz = (Azz_bar.Im()*Azz.Re() - Azz_bar.Re()*Azz.Im())/Upm;

  double Re_Ipm = ((Apm_bar.Re()*Amp.Re() + Apm_bar.Im()*Amp.Im()) - (Amp_bar.Re()*Apm.Re() + Amp_bar.Im()*Apm.Im()))/Upm;
  double Re_Ipz = ((Apm_bar.Re()*Azz.Re() + Apm_bar.Im()*Azz.Im()) - (Azz_bar.Re()*Apm.Re() + Azz_bar.Im()*Apm.Im()))/Upm;
  double Re_Imz = ((Amp_bar.Re()*Azz.Re() + Amp_bar.Im()*Azz.Im()) - (Azz_bar.Re()*Amp.Re() + Azz_bar.Im()*Amp.Im()))/Upm;

  double Im_Ipm = ((Apm_bar.Im()*Amp.Re() - Apm_bar.Re()*Amp.Im()) + (Amp_bar.Im()*Apm.Re() - Amp_bar.Re()*Apm.Im()))/Upm;
  double Im_Ipz = ((Apm_bar.Im()*Azz.Re() - Apm_bar.Re()*Azz.Im()) + (Azz_bar.Im()*Apm.Re() - Azz_bar.Re()*Apm.Im()))/Upm;
  double Im_Imz = ((Amp_bar.Im()*Azz.Re() - Amp_bar.Re()*Azz.Im()) + (Azz_bar.Im()*Amp.Re() - Azz_bar.Re()*Amp.Im()))/Upm;

  double xneu[26];
  xneu[0]=Izz; xneu[1]=Imp; xneu[2]=Im_Imz; xneu[3]=Re_Imz; xneu[4]=Ipm;
  xneu[5]=Im_Ipz; xneu[6]=Re_Ipz; xneu[7]=Im_Ipm; xneu[8]=Re_Ipm; xneu[9]=Vzz;
  xneu[10]=Uzz; xneu[11]=Im_Vmz; xneu[12]=Re_Vmz; xneu[13]=Im_Umz; xneu[14]=Re_Umz;
  xneu[15]=Vmp; xneu[16]=Ump; xneu[17]=Im_Vpz; xneu[18]=Re_Vpz; xneu[19]=Im_Upz;
  xneu[20]=Re_Upz; xneu[21]=Im_Vpm; xneu[22]=Re_Vpm; xneu[23]=Im_Upm; xneu[24]=Re_Upm;
  xneu[25]=Vpm;

  double chisq = 0.;
  if (UseBabar==1) {
    for (int ix=0; ix<26; ix++)
      for (int iy=0; iy<26; iy++)
	chisq += (xneu[ix]-CenNeu_BaBar[ix])*Vneu_BaBar(ix,iy)*(xneu[iy]-CenNeu_BaBar[iy]);
  }
  if (UseBelle==1) {
    double chisqBelle = 0.;
    for (int ix=0; ix<26; ix++)
      for (int iy=0; iy<26; iy++)
	chisqBelle += (xneu[ix]-CenNeu_Belle[ix])*Vneu_Belle(ix,iy)*(xneu[iy]-CenNeu_Belle[iy]);
    chisq += chisqBelle;
  }
  return chisq;
}

Walker Propose(const Walker &w, double stepMag, double stepPhase, TRandom3 *rng) {
  Walker p = w;
  p.magAmp  = ReflectBound(w.magAmp  + rng->Gaus(0.,stepMag),   0., Amp_max);
  p.phAmp   = WrapPhase(w.phAmp + rng->Gaus(0.,stepPhase));
  p.magAzz  = ReflectBound(w.magAzz  + rng->Gaus(0.,stepMag),   0., Azz_max);
  p.phAzz   = WrapPhase(w.phAzz + rng->Gaus(0.,stepPhase));
  p.magApmB = ReflectBound(w.magApmB + rng->Gaus(0.,stepMag),   0., Apm_max);
  p.phApmB  = WrapPhase(w.phApmB + rng->Gaus(0.,stepPhase));
  p.magAmpB = ReflectBound(w.magAmpB + rng->Gaus(0.,stepMag),   0., Amp_max);
  p.phAmpB  = WrapPhase(w.phAmpB + rng->Gaus(0.,stepPhase));
  p.magAzzB = ReflectBound(w.magAzzB + rng->Gaus(0.,stepMag),   0., Azz_max);
  p.phAzzB  = WrapPhase(w.phAzzB + rng->Gaus(0.,stepPhase));
  return p;
}

// Gelman-Rubin R-hat from per-chain Welford mean/M2 accumulators.
double GelmanRubin(const vector<double>& mean, const vector<double>& M2, const vector<long>& n) {
  int J = (int) mean.size();
  double grand=0.; for (int j=0;j<J;j++) grand += mean[j]; grand/=J;
  double B=0.; for (int j=0;j<J;j++) B += pow(mean[j]-grand,2.);
  double L=0.; for (int j=0;j<J;j++) L += n[j]; L/=J;
  if (J<2 || L<2.) return -1.;
  B *= L/(J-1);
  double W=0.; for (int j=0;j<J;j++) W += (n[j]>1) ? M2[j]/(n[j]-1) : 0.; W/=J;
  if (W<=0.) return -1.;
  double varHat = ((L-1.)/L)*W + B/L;
  return sqrt(varHat/W);
}

int main(int argc, char *argv[]){
  if (argc>3) {
    cout << "CKM_gen: input error: please specify  " << endl;
    cout << "just one .conf file name (w/o extension) "   << endl;
    cout << "and(optional) out file name (w/o extension) "   << endl;
    cout << " (e.g. CKM_gen CKM_gen)" << endl;
    return 0;
  }

  TString conffile("CKM_gen_rhopi_mcmc");
  TString outfile("CKM_gen_rhopi_mcmc");
  if (argc==2) {
    conffile = TString((const char *)argv[1],strlen(argv[1]));
    outfile  = TString((const char *)argv[1],strlen(argv[1]));
  }

  if (argc==3) {
    conffile = TString((const char *)argv[1],strlen(argv[1]));
    outfile  = TString((const char *)argv[2],strlen(argv[2]));
  }

  TString conffile_dot_conf = conffile + ".conf";
  TString ouput_dot_root = outfile + ".root";

  const char* filename = conffile_dot_conf;
  Read_Parameters(filename);

  TFile *output = new TFile(ouput_dot_root,"RECREATE");
  output->cd();
  output->mkdir("Input","some input distributions");
  output->cd("Input");

  TH1D *histo1D[12];

  histo1D[0]  = new TH1D("sin2alpharhopi", " sin2alpha",  100,   -1., 1.);
  histo1D[1]  = new TH1D("alpharhopi",     " alpha",      180,    0., 180.);
  histo1D[2]  = new TH1D("AbsAmp",         "AbsAmp",      100,    0., Amp_max);
  histo1D[3]  = new TH1D("AbsAzz",         "AbszzA",      100,    0., Azz_max);
  histo1D[4]  = new TH1D("AbsApm_bar",     "AbsApm_bar",  100,    0., Apm_max);
  histo1D[5]  = new TH1D("AbsAmp_bar",     "AbsAmp_bar",  100,    0., Amp_max);
  histo1D[6]  = new TH1D("AbsAzz_bar",     "AbsAzz_bar",  100,    0., Azz_max);
  histo1D[7]  = new TH1D("ArgAmp",         "ArgAmp",      180, -180., 180.);
  histo1D[8]  = new TH1D("ArgAzz",         "ArgzzA",      180, -180., 180.);
  histo1D[9]  = new TH1D("ArgApm_bar",     "ArgApm_bar",  180, -180., 180.);
  histo1D[10] = new TH1D("ArgAmp_bar",     "ArgAmp_bar",  180, -180., 180.);
  histo1D[11] = new TH1D("ArgAzz_bar",     "ArgAzz_bar",  180, -180., 180.);

  cout << "=========> number of extractions " << NExtractions << endl;

  TMatrixDSym Vneu_BaBar(26);
  TMatrixDSym Vneu_Belle(26);
  Vneu_BaBar = (covStat_BaBar+covSys_BaBar).Invert();
  Vneu_Belle = cov_Belle.Invert();

  int jobpid = getpid();
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int seed = today+clock+jobpid;
  if (Seed != 0) seed = (int) Seed;
  TRandom3 *gRandom = new TRandom3(seed);
  cout << " Seed = " << seed << endl;

  // ---- parallel tempering ----

  long NStepsPerReplica = NExtractions/NReplicas;
  if (NStepsPerReplica<1) NStepsPerReplica=1;

  vector<double> beta(NTemps);
  for (int m=0; m<NTemps; m++)
    beta[m] = (NTemps>1) ? pow(BetaMin, ((double)m)/(NTemps-1)) : 1.;

  double maxRange = std::max(std::max(Apm_max,Amp_max),Azz_max);
  vector<double> stepMagLevel(NTemps), stepPhaseLevel(NTemps);

  const int AdaptWindow = 200;
  vector<long> accCount(NTemps,0), propCount(NTemps,0);           // sampling-phase only, for final report
  vector<long> accWindow(NTemps,0), propWindow(NTemps,0);         // burn-in only, reset every AdaptWindow steps
  vector<long> swapAcc(NTemps-1,0), swapProp(NTemps-1,0);         // sampling-phase only

  vector<double> chainMeanAlpha(NReplicas,0.), chainM2Alpha(NReplicas,0.);
  vector<double> chainMeanAzz(NReplicas,0.),   chainM2Azz(NReplicas,0.);
  vector<long>   chainCount(NReplicas,0);

  cout << "=========> " << NReplicas << " replicas x " << NStepsPerReplica
       << " cold-chain samples each (" << NBurnin << " burn-in steps), "
       << NTemps << " temperature levels down to beta=" << BetaMin << endl;

  for (int rep=0; rep<NReplicas; rep++) {

    // Each replica is an independent chain: it re-runs burn-in from a fresh
    // random start, so the adaptive step sizes and their acceptance counters
    // must start from the nominal values too. Carrying them over from the
    // previous replica lets the burn-in growth factor compound across
    // replicas without bound.
    for (int m=0; m<NTemps; m++) {
      double boost = 1./sqrt(beta[m]);
      stepMagLevel[m]   = std::min(StepMag*boost,   0.5*maxRange);
      stepPhaseLevel[m] = std::min(StepPhase*boost, Pi);
      accWindow[m] = propWindow[m] = 0;
    }

    vector<Walker> state(NTemps);
    for (int m=0; m<NTemps; m++) {
      state[m] = RandomWalker(gRandom);
      state[m].chisq = ChiSqFromWalker(state[m], Vneu_BaBar, Vneu_Belle, 0);
    }

    bool burning = true;
    long totalSteps = NBurnin + NStepsPerReplica;

    for (long step=1; step<=totalSteps; step++) {

      if (burning && step>NBurnin) burning = false;

      // per-level Metropolis update
      for (int m=0; m<NTemps; m++) {
	Walker prop = Propose(state[m], stepMagLevel[m], stepPhaseLevel[m], gRandom);
	prop.chisq = ChiSqFromWalker(prop, Vneu_BaBar, Vneu_Belle, 0);
	double logA = 0.5*beta[m]*(state[m].chisq - prop.chisq);
	bool accept = (logA>=0. || gRandom->Rndm() < exp(logA));
	if (accept) state[m] = prop;
	if (burning) {
	  propWindow[m]++;
	  if (accept) accWindow[m]++;
	} else {
	  propCount[m]++;
	  if (accept) accCount[m]++;
	}
      }

      // adaptive step-size tuning (burn-in only), target 20-40% acceptance
      if (burning && step%AdaptWindow==0) {
	for (int m=0; m<NTemps; m++) {
	  if (propWindow[m]==0) continue;
	  double rate = ((double)accWindow[m])/propWindow[m];
	  double factor = (rate<0.2) ? 0.8 : ((rate>0.4) ? 1.25 : 1.0);
	  stepMagLevel[m]   = std::min(stepMagLevel[m]*factor, 0.5*maxRange);
	  stepPhaseLevel[m] = std::min(stepPhaseLevel[m]*factor, Pi);
	  accWindow[m]=0; propWindow[m]=0;
	}
      }

      // replica-exchange (swap) moves between adjacent temperature levels
      if (step%SwapEvery==0) {
	for (int m=0; m<NTemps-1; m++) {
	  double logA = 0.5*(beta[m]-beta[m+1])*(state[m].chisq-state[m+1].chisq);
	  bool accept = (logA>=0. || gRandom->Rndm() < exp(logA));
	  if (accept) std::swap(state[m], state[m+1]);
	  if (!burning) {
	    swapProp[m]++;
	    if (accept) swapAcc[m]++;
	  }
	}
      }

      if (!burning) {
	// fill output histograms from the cold (beta=1) chain only, weight 1
	double alpha;
	ChiSqFromWalker(state[0], Vneu_BaBar, Vneu_Belle, &alpha);
	double alphaDeg = alpha*180./Pi;

	histo1D[0]->Fill(sin(2.*alpha));
	histo1D[1]->Fill(alphaDeg);
	histo1D[2]->Fill(state[0].magAmp);
	histo1D[3]->Fill(state[0].magAzz);
	histo1D[4]->Fill(state[0].magApmB);
	histo1D[5]->Fill(state[0].magAmpB);
	histo1D[6]->Fill(state[0].magAzzB);
	histo1D[7]->Fill(state[0].phAmp*180./Pi);
	histo1D[8]->Fill(state[0].phAzz*180./Pi);
	histo1D[9]->Fill(state[0].phApmB*180./Pi);
	histo1D[10]->Fill(state[0].phAmpB*180./Pi);
	histo1D[11]->Fill(state[0].phAzzB*180./Pi);

	// Welford running mean/variance, per replica, for the Gelman-Rubin check
	chainCount[rep]++;
	double dA = alphaDeg - chainMeanAlpha[rep];
	chainMeanAlpha[rep] += dA/chainCount[rep];
	chainM2Alpha[rep]   += dA*(alphaDeg-chainMeanAlpha[rep]);

	double dZ = state[0].magAzz - chainMeanAzz[rep];
	chainMeanAzz[rep] += dZ/chainCount[rep];
	chainM2Azz[rep]   += dZ*(state[0].magAzz-chainMeanAzz[rep]);
      }
    }

    cout << "replica " << rep+1 << "/" << NReplicas << " done" << endl;
  }

  cout << "=== per-level acceptance rate (sampling phase) ===" << endl;
  for (int m=0;m<NTemps;m++)
    cout << "  level " << m << " (beta=" << beta[m] << "): "
	 << (propCount[m]? ((double)accCount[m])/propCount[m] : 0.) << endl;

  cout << "=== adjacent-pair swap acceptance rate (sampling phase) ===" << endl;
  for (int m=0;m<NTemps-1;m++)
    cout << "  " << m << "<->" << m+1 << ": "
	 << (swapProp[m]? ((double)swapAcc[m])/swapProp[m] : 0.) << endl;

  cout << "=== Gelman-Rubin R-hat across " << NReplicas << " replicas (want ~1.0-1.1) ===" << endl;
  cout << "  alpha : " << GelmanRubin(chainMeanAlpha, chainM2Alpha, chainCount) << endl;
  cout << "  |Azz| : " << GelmanRubin(chainMeanAzz,   chainM2Azz,   chainCount) << endl;

  for(int h=0; h<12; h++) {
    //normalize
    if (histo1D[h]) {
      cout << "here: " << h << endl;
      double norm = 1./(histo1D[h]->Integral());
      cout << norm << endl;
      histo1D[h]->Scale(norm);
      histo1D[h]->Write();
    }
  }
  output->Close();
}


void Read_Parameters(const char* filename) {
  char buffer[256];
  ifstream reader (filename);
  map<string, double> data;
  string label, equal;
  double value;

  if ( ! reader.is_open())
    { cout << "Error opening file"; exit (1); }

  while ( reader.getline (buffer,256) )
    {
      istringstream i(buffer);
      i >> label
        >> equal
        >> value;
      data[label] = value;
    }

  map<string, double>::const_iterator iter;
  for (iter=data.begin(); iter != data.end(); iter++) {
    cout << iter->second << " " << iter->first << endl;
  }

  Assign_Parameters(data);

  return;
}

void Assign_Parameters(map<string, double> data) {

//     Parameters

  NExtractions = int (data["NExtractions"]);
  Opt_NP = int (data["Opt_NP"]);
  Seed = long (data["Seed"]);

  // parallel-tempering parameters, with sensible fallback defaults so
  // existing confs that predate this scheme still work.
  NReplicas = int (data["NReplicas"]); if (NReplicas<=0) NReplicas=10;
  NTemps    = int (data["NTemps"]);    if (NTemps<=1)    NTemps=10;
  BetaMin   = data["BetaMin"];         if (BetaMin<=0. || BetaMin>=1.) BetaMin=1.e-3;
  NBurnin   = int (data["NBurnin"]);   if (NBurnin<=0)   NBurnin=2000;
  SwapEvery = int (data["SwapEvery"]); if (SwapEvery<=0) SwapEvery=20;
  StepMag   = data["StepMag"];         if (StepMag<=0.)  StepMag=0.1;
  StepPhase = data["StepPhase"];       if (StepPhase<=0.) StepPhase=0.3;

  UseBabar = int (data["UseBabar"]);
  UseBelle = int (data["UseBelle"]);

  Apm_max = data["Apm_max"];
  Amp_max = data["Amp_max"];
  Azz_max = data["Azz_max"];

  CenNeu_BaBar[0] = data[("CenIzz_BaBar")];
  CenNeu_BaBar[1] = data[("CenImp_BaBar")];
  CenNeu_BaBar[2] = data[("CenIm_Imz_BaBar")];
  CenNeu_BaBar[3] = data[("CenRe_Imz_BaBar")];
  CenNeu_BaBar[4] = data[("CenIpm_BaBar")];
  CenNeu_BaBar[5] = data[("CenIm_Ipz_BaBar")];
  CenNeu_BaBar[6] = data[("CenRe_Ipz_BaBar")];
  CenNeu_BaBar[7] = data[("CenIm_Ipm_BaBar")];
  CenNeu_BaBar[8] = data[("CenRe_Ipm_BaBar")];
  CenNeu_BaBar[9] = data[("CenVzz_BaBar")];
  CenNeu_BaBar[10] = data[("CenUzz_BaBar")];
  CenNeu_BaBar[11] = data[("CenIm_Vmz_BaBar")];
  CenNeu_BaBar[12] = data[("CenRe_Vmz_BaBar")];
  CenNeu_BaBar[13] = data[("CenIm_Umz_BaBar")];
  CenNeu_BaBar[14] = data[("CenRe_Umz_BaBar")];
  CenNeu_BaBar[15] = data[("CenVmp_BaBar")];
  CenNeu_BaBar[16] = data[("CenUmp_BaBar")];
  CenNeu_BaBar[17] = data[("CenIm_Vpz_BaBar")];
  CenNeu_BaBar[18] = data[("CenRe_Vpz_BaBar")];
  CenNeu_BaBar[19] = data[("CenIm_Upz_BaBar")];
  CenNeu_BaBar[20] = data[("CenRe_Upz_BaBar")];
  CenNeu_BaBar[21] = data[("CenIm_Vpm_BaBar")];
  CenNeu_BaBar[22] = data[("CenRe_Vpm_BaBar")];
  CenNeu_BaBar[23] = data[("CenIm_Upm_BaBar")];
  CenNeu_BaBar[24] = data[("CenRe_Upm_BaBar")];
  CenNeu_BaBar[25] = data[("CenVpm_BaBar")];

  SigNeu_BaBar[0] = data[("SigIzz_BaBar")];
  SigNeu_BaBar[1] = data[("SigImp_BaBar")];
  SigNeu_BaBar[2] = data[("SigIm_Imz_BaBar")];
  SigNeu_BaBar[3] = data[("SigRe_Imz_BaBar")];
  SigNeu_BaBar[4] = data[("SigIpm_BaBar")];
  SigNeu_BaBar[5] = data[("SigIm_Ipz_BaBar")];
  SigNeu_BaBar[6] = data[("SigRe_Ipz_BaBar")];
  SigNeu_BaBar[7] = data[("SigIm_Ipm_BaBar")];
  SigNeu_BaBar[8] = data[("SigRe_Ipm_BaBar")];
  SigNeu_BaBar[9] = data[("SigVzz_BaBar")];
  SigNeu_BaBar[10] = data[("SigUzz_BaBar")];
  SigNeu_BaBar[11] = data[("SigIm_Vmz_BaBar")];
  SigNeu_BaBar[12] = data[("SigRe_Vmz_BaBar")];
  SigNeu_BaBar[13] = data[("SigIm_Umz_BaBar")];
  SigNeu_BaBar[14] = data[("SigRe_Umz_BaBar")];
  SigNeu_BaBar[15] = data[("SigVmp_BaBar")];
  SigNeu_BaBar[16] = data[("SigUmp_BaBar")];
  SigNeu_BaBar[17] = data[("SigIm_Vpz_BaBar")];
  SigNeu_BaBar[18] = data[("SigRe_Vpz_BaBar")];
  SigNeu_BaBar[19] = data[("SigIm_Upz_BaBar")];
  SigNeu_BaBar[20] = data[("SigRe_Upz_BaBar")];
  SigNeu_BaBar[21] = data[("SigIm_Vpm_BaBar")];
  SigNeu_BaBar[22] = data[("SigRe_Vpm_BaBar")];
  SigNeu_BaBar[23] = data[("SigIm_Upm_BaBar")];
  SigNeu_BaBar[24] = data[("SigRe_Upm_BaBar")];
  SigNeu_BaBar[25] = data[("SigVpm_BaBar")];

  SysNeu_BaBar[0] = data[("SysIzz_BaBar")];
  SysNeu_BaBar[1] = data[("SysImp_BaBar")];
  SysNeu_BaBar[2] = data[("SysIm_Imz_BaBar")];
  SysNeu_BaBar[3] = data[("SysRe_Imz_BaBar")];
  SysNeu_BaBar[4] = data[("SysIpm_BaBar")];
  SysNeu_BaBar[5] = data[("SysIm_Ipz_BaBar")];
  SysNeu_BaBar[6] = data[("SysRe_Ipz_BaBar")];
  SysNeu_BaBar[7] = data[("SysIm_Ipm_BaBar")];
  SysNeu_BaBar[8] = data[("SysRe_Ipm_BaBar")];
  SysNeu_BaBar[9] = data[("SysVzz_BaBar")];
  SysNeu_BaBar[10] = data[("SysUzz_BaBar")];
  SysNeu_BaBar[11] = data[("SysIm_Vmz_BaBar")];
  SysNeu_BaBar[12] = data[("SysRe_Vmz_BaBar")];
  SysNeu_BaBar[13] = data[("SysIm_Umz_BaBar")];
  SysNeu_BaBar[14] = data[("SysRe_Umz_BaBar")];
  SysNeu_BaBar[15] = data[("SysVmp_BaBar")];
  SysNeu_BaBar[16] = data[("SysUmp_BaBar")];
  SysNeu_BaBar[17] = data[("SysIm_Vpz_BaBar")];
  SysNeu_BaBar[18] = data[("SysRe_Vpz_BaBar")];
  SysNeu_BaBar[19] = data[("SysIm_Upz_BaBar")];
  SysNeu_BaBar[20] = data[("SysRe_Upz_BaBar")];
  SysNeu_BaBar[21] = data[("SysIm_Vpm_BaBar")];
  SysNeu_BaBar[22] = data[("SysRe_Vpm_BaBar")];
  SysNeu_BaBar[23] = data[("SysIm_Upm_BaBar")];
  SysNeu_BaBar[24] = data[("SysRe_Upm_BaBar")];
  SysNeu_BaBar[25] = data[("SysVpm_BaBar")];

  CenNeu_Belle[0] = data[("CenIzz_Belle")];
  CenNeu_Belle[1] = data[("CenImp_Belle")];
  CenNeu_Belle[2] = data[("CenIm_Imz_Belle")];
  CenNeu_Belle[3] = data[("CenRe_Imz_Belle")];
  CenNeu_Belle[4] = data[("CenIpm_Belle")];
  CenNeu_Belle[5] = data[("CenIm_Ipz_Belle")];
  CenNeu_Belle[6] = data[("CenRe_Ipz_Belle")];
  CenNeu_Belle[7] = data[("CenIm_Ipm_Belle")];
  CenNeu_Belle[8] = data[("CenRe_Ipm_Belle")];
  CenNeu_Belle[9] = data[("CenVzz_Belle")];
  CenNeu_Belle[10] = data[("CenUzz_Belle")];
  CenNeu_Belle[11] = data[("CenIm_Vmz_Belle")];
  CenNeu_Belle[12] = data[("CenRe_Vmz_Belle")];
  CenNeu_Belle[13] = data[("CenIm_Umz_Belle")];
  CenNeu_Belle[14] = data[("CenRe_Umz_Belle")];
  CenNeu_Belle[15] = data[("CenVmp_Belle")];
  CenNeu_Belle[16] = data[("CenUmp_Belle")];
  CenNeu_Belle[17] = data[("CenIm_Vpz_Belle")];
  CenNeu_Belle[18] = data[("CenRe_Vpz_Belle")];
  CenNeu_Belle[19] = data[("CenIm_Upz_Belle")];
  CenNeu_Belle[20] = data[("CenRe_Upz_Belle")];
  CenNeu_Belle[21] = data[("CenIm_Vpm_Belle")];
  CenNeu_Belle[22] = data[("CenRe_Vpm_Belle")];
  CenNeu_Belle[23] = data[("CenIm_Upm_Belle")];
  CenNeu_Belle[24] = data[("CenRe_Upm_Belle")];
  CenNeu_Belle[25] = data[("CenVpm_Belle")];

  SigTot_Belle[0] = sqrt(pow(data[("SigIzz_Belle")],2.)+pow(data[("SysIzz_Belle")],2.));
  SigTot_Belle[1] = sqrt(pow(data[("SigImp_Belle")],2.)+pow(data[("SysImp_Belle")],2.));
  SigTot_Belle[2] = sqrt(pow(data[("SigIm_Imz_Belle")],2.)+pow(data[("SysIm_Imz_Belle")],2.));
  SigTot_Belle[3] = sqrt(pow(data[("SigRe_Imz_Belle")],2.)+pow(data[("SysRe_Imz_Belle")],2.));
  SigTot_Belle[4] = sqrt(pow(data[("SigIpm_Belle")],2.)+pow(data[("SysIpm_Belle")],2.));
  SigTot_Belle[5] = sqrt(pow(data[("SigIm_Ipz_Belle")],2.)+pow(data[("SysIm_Ipz_Belle")],2.));
  SigTot_Belle[6] = sqrt(pow(data[("SigRe_Ipz_Belle")],2.)+pow(data[("SysRe_Ipz_Belle")],2.));
  SigTot_Belle[7] = sqrt(pow(data[("SigIm_Ipm_Belle")],2.)+pow(data[("SysIm_Ipm_Belle")],2.));
  SigTot_Belle[8] = sqrt(pow(data[("SigRe_Ipm_Belle")],2.)+pow(data[("SysRe_Ipm_Belle")],2.));
  SigTot_Belle[9] = sqrt(pow(data[("SigVzz_Belle")],2.)+pow(data[("SysVzz_Belle")],2.));
  SigTot_Belle[10] = sqrt(pow(data[("SigUzz_Belle")],2.)+pow(data[("SysUzz_Belle")],2.));
  SigTot_Belle[11] = sqrt(pow(data[("SigIm_Vmz_Belle")],2.)+pow(data[("SysIm_Vmz_Belle")],2.));
  SigTot_Belle[12] = sqrt(pow(data[("SigRe_Vmz_Belle")],2.)+pow(data[("SysRe_Vmz_Belle")],2.));
  SigTot_Belle[13] = sqrt(pow(data[("SigIm_Umz_Belle")],2.)+pow(data[("SysIm_Umz_Belle")],2.));
  SigTot_Belle[14] = sqrt(pow(data[("SigRe_Umz_Belle")],2.)+pow(data[("SysRe_Umz_Belle")],2.));
  SigTot_Belle[15] = sqrt(pow(data[("SigVmp_Belle")],2.)+pow(data[("SysVmp_Belle")],2.));
  SigTot_Belle[16] = sqrt(pow(data[("SigUmp_Belle")],2.)+pow(data[("SysUmp_Belle")],2.));
  SigTot_Belle[17] = sqrt(pow(data[("SigIm_Vpz_Belle")],2.)+pow(data[("SysIm_Vpz_Belle")],2.));
  SigTot_Belle[18] = sqrt(pow(data[("SigRe_Vpz_Belle")],2.)+pow(data[("SysRe_Vpz_Belle")],2.));
  SigTot_Belle[19] = sqrt(pow(data[("SigIm_Upz_Belle")],2.)+pow(data[("SysIm_Upz_Belle")],2.));
  SigTot_Belle[20] = sqrt(pow(data[("SigRe_Upz_Belle")],2.)+pow(data[("SysRe_Upz_Belle")],2.));
  SigTot_Belle[21] = sqrt(pow(data[("SigIm_Vpm_Belle")],2.)+pow(data[("SysIm_Vpm_Belle")],2.));
  SigTot_Belle[22] = sqrt(pow(data[("SigRe_Vpm_Belle")],2.)+pow(data[("SysRe_Vpm_Belle")],2.));
  SigTot_Belle[23] = sqrt(pow(data[("SigIm_Upm_Belle")],2.)+pow(data[("SysIm_Upm_Belle")],2.));
  SigTot_Belle[24] = sqrt(pow(data[("SigRe_Upm_Belle")],2.)+pow(data[("SysRe_Upm_Belle")],2.));
  SigTot_Belle[25] = sqrt(pow(data[("SigVpm_Belle")],2.)+pow(data[("SysVpm_Belle")],2.));

  // BaBar statistical correlation
  for(int ix=0;ix<26;ix++) {
    for(int iy=ix;iy<26;iy++) {
      char nameneu[256] ;
      sprintf(nameneu,"corStatBaBar_%d_%d",ix+1,iy+1);
      covStat_BaBar(ix,iy) = data[(nameneu)]*SigNeu_BaBar[ix]*SigNeu_BaBar[iy];
      cout << " read stat " << covStat_BaBar(ix,iy) << " index fixing for " << ix << " and " << iy << endl;
    }
  }

  // to make it symmetric without inputting same numbers twice
  for(int ix=0;ix<26;ix++) {
    for(int iy=0;iy<ix;iy++) {
      covStat_BaBar(iy,ix) = covStat_BaBar(iy,ix);
    }
  }

  // BaBar systematic correlation
  for(int ix=0;ix<26;ix++) {
    for(int iy=ix;iy<26;iy++) {
      char nameneu[256] ;
      sprintf(nameneu,"corSysBaBar_%d_%d",ix+1,iy+1);
      covSys_BaBar(ix,iy) = data[(nameneu)]*SysNeu_BaBar[ix]*SysNeu_BaBar[iy];
      cout << " read sys " << covSys_BaBar(ix,iy) << " index fixing for " << ix << " and " << iy << endl;
    }
  }

  // to make it symmetric without inputting same numbers twice
  for(int ix=0;ix<26;ix++) {
    for(int iy=0;iy<ix;iy++) {
      covSys_BaBar(ix,iy) = covSys_BaBar(iy,ix);
    }
  }


  for(int ix=0;ix<26;ix++) {
    for(int iy=ix;iy<26;iy++) {
      char nameneu[256] ;
      sprintf(nameneu,"corBelle_%d_%d",ix+1,iy+1);
      cov_Belle(ix,iy) = data[(nameneu)]*SigTot_Belle[ix]*SigTot_Belle[iy];
    }
  }

  for(int ix=0;ix<26;ix++) {
    for(int iy=0;iy<ix;iy++) {
      // to make it symmetric without inputting same numbers twice
      cov_Belle(ix,iy) = cov_Belle(iy,ix);
      }
    }

}

Double_t myGauss(Double_t *x, Double_t *par) {

  double xx =x[0];
  return exp(-0.5*(xx-par[0])*(xx-par[0])/(par[1]*par[1]));
}
