#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#include <iostream>
#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TDatime.h>
#include <TF1.h>
#include <string>
#include <vector>
#include <TString.h>
#include <TComplex.h>
#include <RooRealVar.h>
#include <RooMultiVarGaussian.h>
#include <RooDataSet.h>
#include "CKM_gen.hh"

using namespace std;
double Pi = acos(-1.);
bool lfirst = true;

int main(int argc, char *argv[]){
  if (argc>3) {
    cout << "CKM_gen: input error: please specify  " << endl;
    cout << "just one .conf file name (w/o extension) "   << endl;
    cout << "and(optional) out file name (w/o extension) "   << endl;
    cout << " (e.g. CKM_gen CKM_gen)" << endl;
    return 0;
  } 

  TString conffile("CKM_gen");
  TString outfile("CKM_gen");
  if (argc==2) {
    conffile = TString((const char *)argv[1],strlen(argv[1]));
    outfile  = TString((const char *)argv[1],strlen(argv[1]));
  }

  if (argc==3) {
    conffile = TString((const char *)argv[1],strlen(argv[1]));
    outfile  = TString((const char *)argv[2],strlen(argv[2]));
  }

  TString conffile_dot_conf = conffile + ".conf";
  TString output_dot_root = outfile + ".root";
  
  cout << output_dot_root << endl;

  const char* filename = conffile_dot_conf;
  Read_Parameters(filename);

  //  double alphaeff, stwoalpha;;

  //TH1D * distr_skepC =0;
  //TH1D * distr_skepS =0;

  // lets put input values by hand... waiting for nerd's suggestions
  int nrho = 192, neta = 192;
  double rho_min(-1.25), rho_max(1.25);
  double eta_min(-0.25), eta_max(1.25);

  std::vector<TH1D*> histo1D;
  std::vector<TH1D*> prior1D;
  std::vector<TH2D*> histo2D;

  /*
  TH1D *histo1D[15];
  TH1D *prior1D[15];
  TH2D *histo2D[27];
  for (int h=0; h<15; h++) {
    //histo1D[h]=NULL;
    prior1D[h]=NULL;
  }
  for (int h=0; h<27; h++) {
    histo2D[h]=NULL;
  }
  */
  
  /*
  if ( skepC != "0"){
      TString string_skepC = (const char*) skepC.c_str();
      string_skepC += ".root";
      TFile *input_skepC = new TFile(string_skepC);
      if ((TH1D*) input_skepC->Get("Input/input_C")) 
	distr_skepC = (TH1D*) input_skepC->Get("Input/input_C")->Clone();
      else 
	cout << "no C distribution histogram in the given rootfile" << endl;      
  }
  
  if ( skepS != "0"){
    TString string_skepS = (const char*) skepS.c_str();
    string_skepS += ".root";
    TFile *input_skepS = new TFile(string_skepS);
    if ((TH1D*) input_skepS->Get("Input/input_S") ) 
      distr_skepS = (TH1D*) input_skepS->Get("Input/input_S")->Clone();
    else 
      cout << "no S distribution histogram in the given rootfile" << endl;
  }
  */

  TH1D* priorPs = new TH1D("prior_Ps", "PriorPs", 100, 0., 10.);
  // create histogram for P (from Bs->KK)
  // to use in the weight
  if(SigPs != 0. && FlatPs != 0.) {
    cout << "Bs->KK prior for P applied" << endl;
    cout << "Change P range from [" << minP << ","<< maxP <<"] to [0.,CenPs+4*(SigPs+FlatPs)]" << endl;
    minP = 0.;
    maxP = CenPs+4*(SigPs+FlatPs);
    TF1  *Ps_dist = new TF1("Ps_dist",Hfun,minP,maxP,3);
    Ps_dist->SetParameters(CenPs,SigPs,FlatPs);
    priorPs->FillRandom("Ps_dist",200000);
    priorPs->Scale(1./priorPs->Integral());
  }

  histo1D.push_back(new TH1D((channel+"_input_sin2alpha").c_str(), (channel+"input sin2alpha").c_str(),  100, -1., 1.));
  histo1D.push_back(new TH1D((channel+"_input_alpha").c_str(),     (channel+" input alpha").c_str(),      180, 0., 180.));
  histo1D.push_back(new TH1D((channel+"_input_BRpm").c_str(),      (channel+" input BRpm").c_str(),    
			 100, MAX(0.,CenBRpm-6.*SigBRpm), CenBRpm+6.*SigBRpm));
  histo1D.push_back(new TH1D((channel+"_input_BRpz").c_str(),      (channel+" input BRpz").c_str(),    
			 100, MAX(0.,CenBRpz-6.*SigBRpz), CenBRpz+6.*SigBRpz)); 
  histo1D.push_back(new TH1D((channel+"_input_BRzz").c_str(),      (channel+" input BRzz").c_str(),    
			 100, MAX(0.,CenBRzz-6.*SigBRzz), CenBRzz+6.*SigBRzz));;
  histo1D.push_back(new TH1D((channel+"_input_C").c_str(),         (channel+" input C").c_str(),       100, -1., 1.));
  histo1D.push_back(new TH1D((channel+"_input_S").c_str(),         (channel+" input S").c_str(),       100, -1., 1.));
  histo1D.push_back(new TH1D((channel+"_input_Czz").c_str(),       (channel+" input Czz").c_str(),     100, -1., 1.));
  histo1D.push_back(new TH1D((channel+"_input_Szz").c_str(),       (channel+" input Szz").c_str(),     100, -1., 1.));
  histo1D.push_back(new TH1D((channel+"_input_P").c_str(),         (channel+" input P").c_str(),       100, minP,   maxP));
  histo1D.push_back(new TH1D((channel+"_input_T").c_str(),         (channel+" input T").c_str(),       100, minTPM, maxTPM));
  histo1D.push_back(new TH1D((channel+"_input_Tc").c_str(),        (channel+" input Tc").c_str(),      100, minT00, maxT00));
  histo1D.push_back(new TH1D((channel+"_input_PoverT").c_str(),    (channel+" input PoverT").c_str(),  100, 0., 5.));
  histo1D.push_back(new TH1D((channel+"_input_deltaP").c_str(),    (channel+" input deltaP").c_str(),  100, 0., 360.));
  histo1D.push_back(new TH1D((channel+"_input_deltaTc").c_str(),   (channel+" input deltaTc").c_str(), 100, 0., 360.));

  histo2D.push_back(new TH2D("alpha", "eta vs rho alpha", nrho, rho_min, rho_max, neta, eta_min, eta_max));

  histo2D.push_back(new TH2D("BRpm_vs_BRpz", "BRpm vs BRpz", 
			250, MAX(0.,CenBRpz-6.*SigBRpz), CenBRpz+6.*SigBRpz, 
			250, MAX(0.,CenBRpm-6.*SigBRpm), CenBRpm+6.*SigBRpm));
  histo2D.push_back(new TH2D("BRpm_vs_BRzz", "BRpm vs BRzz", 
			250, MAX(0.,CenBRzz-6.*SigBRzz), CenBRzz+6.*SigBRzz, 
			250, MAX(0.,CenBRpm-6.*SigBRpm), CenBRpm+6.*SigBRpm));
  histo2D.push_back(new TH2D("BRpz_vs_BRzz", "BRpz vs BRzz", 
			250, MAX(0.,CenBRzz-6.*SigBRzz), CenBRzz+6.*SigBRzz, 
			250, MAX(0.,CenBRpz-6.*SigBRpz), CenBRpz+6.*SigBRpz));
  histo2D.push_back(new TH2D("alpha_vs_C", "alpha vs C",                 100, -1.,   1., 100, 0., 180.)); //4
  histo2D.push_back(new TH2D("alpha_vs_ToverP", "alpha_vs_ToverP",       180,  0.,  18., 180, 0., 180.));
  histo2D.push_back(new TH2D("alpha_vs_T", "alpha vs T",                 180,  0., 180., 100, minTPM, maxTPM));
  histo2D.push_back(new TH2D("alpha_vs_P", "alpha vs P",                 180,  0., 180., 100, minP,   maxP));
  histo2D.push_back(new TH2D("alpha_vs_Tc", "alpha vs Tc",               180,  0., 180., 100, minT00, maxT00));
  histo2D.push_back(new TH2D("alpha_vs_deltaP", "alpha vs deltaP",       180,  0., 180., 180, 0., 180.)); //9
  //histo2D.push_back(new TH2D("alpha_vs_deltaTc", "alpha vs deltaTc",     180,  0., 180., 180, 0., 180.));
  //histo2D.push_back(new TH2D("deltaPb_vs_T", "deltaPb vs T",             180,  0., 180., 100, minTPM, maxTPM));
  //histo2D.push_back(new TH2D("deltaPb_vs_P", "deltaPb vs P",             180,  0., 180., 100, minP,   maxP));
  //histo2D.push_back(new TH2D("deltaPb_vs_Tc", "deltaPb vs Tc",           180,  0., 180., 100, minT00, maxT00));
  //histo2D.push_back(new TH2D("deltaPb_vs_deltaP", "deltaPb vs deltaP",   180,  0., 180., 180, 0., 180.));
  //histo2D.push_back(new TH2D("deltaPb_vs_deltaTc", "deltaPb vs deltaTc", 180,  0., 180., 180, 0., 180.));
  histo2D.push_back(new TH2D("T_vs_P", "T vs P",                         100,  minTPM, maxTPM, 100, minP,   maxP)); //10
  histo2D.push_back(new TH2D("T_vs_Tc", "T vs Tc",                       100,  minTPM, maxTPM, 100, minT00, maxT00));
  histo2D.push_back(new TH2D("T_vs_deltaP", "T vs deltaP",               100,  minTPM, maxTPM, 180, 0., 180.));
  histo2D.push_back(new TH2D("T_vs_deltaTc", "T vs deltaTc",             100,  minTPM, maxTPM, 180, 0., 180.));
  histo2D.push_back(new TH2D("P_vs_Tc", "P vs Tc",                       100,  minT00, maxT00, 100, minP, maxP));
  histo2D.push_back(new TH2D("P_vs_deltaP", "P vs deltaP",               100,  minP, maxP, 180, 0., 180.));
  histo2D.push_back(new TH2D("P_vs_deltaTc", "P vs deltaTc",             100,  minP, maxP, 180, 0., 180.));
  histo2D.push_back(new TH2D("Tc_vs_deltaP", "Tc vs deltaP",             100,  minT00, maxT00, 180, 0., 180.));
  histo2D.push_back(new TH2D("Tc_vs_deltaTc", "Tc vs deltaTc",           100,  minT00, maxT00, 180, 0., 180.));
  histo2D.push_back(new TH2D("deltaP_vs_deltaTc", "deltaP vs deltaTc",   180,  minT00, maxT00, 180, 0., 180.));//19
  
  prior1D.push_back(new TH1D("prior_sin2alpha","prior sin2alpha", 100, -1., 1.));
  prior1D.push_back(new TH1D("prior_alpha","prior alpha", 180, 0., 180.));
  prior1D.push_back(new TH1D("prior_BRpm", "prior BRpm", 100, 0., 50.));
  prior1D.push_back(new TH1D("prior_BRpz", "prior BRpz", 100, 0., 50.));
  prior1D.push_back(new TH1D("prior_BRzz", "prior BRzz", 100, 0., 5.));
  prior1D.push_back(new TH1D("prior_C", "prior C", 100, -1., 1.));
  prior1D.push_back(new TH1D("prior_S", "prior S", 100, -1., 1.));
  prior1D.push_back(new TH1D("prior_Czz", "prior Czz", 100, -1., 1.));
  prior1D.push_back(new TH1D("prior_Szz", "prior Szz", 100, -1., 1.));
  prior1D.push_back(new TH1D("prior_P", "prior P", 100, minP, maxP));
  prior1D.push_back(new TH1D("prior_T", "prior T", 100, minTPM, maxTPM));
  prior1D.push_back(new TH1D("prior_Tc", "prior Tc", 100, minT00, maxT00));
  prior1D.push_back(new TH1D("prior_deltaP", "prior deltaP", 100, 0., 360.));
  prior1D.push_back(new TH1D("prior_deltaTc", "prior deltaTc", 100, 0., 360.));
  //prior1D.push_back(new TH1D("prior_deltaPb", "prior deltaPb", 100, 0., 180.));

  
  // correlation for C and S
  RooRealVar  Scp("Scp","Scp",-0.9,0,"") ;
  RooRealVar  Ccp("Ccp","Ccp",-0.9,0,"") ;

  TVectorD CS_combo(2) ;
  CS_combo(0) = CenS2a ; CS_combo(1) = CenC2a ;

  TMatrixDSym CS_combo_V(2) ;
  double CS_combo_sigma1 = SigS2a;
  double CS_combo_sigma2 = SigC2a;
  double CS_combo_rho1 = 0.22;
  
  CS_combo_V(0,0) = CS_combo_sigma1*CS_combo_sigma1 ; CS_combo_V(0,1) = CS_combo_rho1*CS_combo_sigma1*CS_combo_sigma2 ;
  CS_combo_V(1,0) = CS_combo_rho1*CS_combo_sigma1*CS_combo_sigma2  ; CS_combo_V(1,1) = CS_combo_sigma2*CS_combo_sigma2;

  RooMultiVarGaussian mvg("mvg","mvg",RooArgList(Scp,Ccp),CS_combo,CS_combo_V) ;
  RooDataSet* d_CS_combo = mvg.generate(RooArgSet(Scp,Ccp),1000000) ;
  Int_t bins_a = 1000;
  Int_t bins_b = 1000;

  TH2F * myhist_CS_combo = (TH2F*) d_CS_combo->createHistogram(Scp,Ccp, bins_a, bins_b);

  
  // seed initialization 
  cout << "=========> number of extractions " << NExtractions << endl;
  int jobpid = getpid();
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int seed = today+clock+jobpid;
  TRandom3 *gRandom = new TRandom3(seed);
  cout << " Seed = " << seed << endl;  
  for(int k=1; k<=NExtractions; k++) {
    if(fmod(k,1000000.) == 0) cout << k << " events generated" << endl;

    // random generation of BR's and CP asy
    BRpm = gRandom->Gaus(CenBRpm,SigBRpm);
    BRpz = gRandom->Gaus(CenBRpz,SigBRpz);
    BRzz = gRandom->Gaus(CenBRzz,SigBRzz);
    //C = gRandom->Gaus(CenC2a,SigC2a);
    //S = gRandom->Gaus(CenS2a,SigS2a);
    myhist_CS_combo->GetRandom2(S,C);
    if(SigCzz2a != 0)
      Czz = gRandom->Gaus(CenCzz2a,SigCzz2a);
    else
      Czz = gRandom->Rndm()*2.-1.;
    // random generation of polarizations
    FLong_pm = gRandom->Gaus(CenFLpm,SigFLpm);
    FLong_pz = gRandom->Gaus(CenFLpz,SigFLpz);
    FLong_zz = gRandom->Gaus(CenFLzz,SigFLzz);
    
    // Sanity check of physical boundaries
    if (S*S+C*C >1. || fabs(Czz) > 1. ||
	fabs(FLong_pm) > 1.   ||
	fabs(FLong_pz) > 1.   ||
	fabs(FLong_zz) > 1.) continue; 
    
    // I seguenti valori servono a eliminare tB0 dal b->pi0 pi0,
    // mantenendo la struttura della rottura di SU(2)
    // infatti, seguendo hep-ph 0607246, ho scritto B00=tB0/2 (|A|^2+|Ab|^2)/2
    // mettendo tB0=2, equivale ad avere B00=(|A|^2+|Ab|^2)/2
    LifetimeRatio = gRandom->Gaus(CenLtRatio,SigLtRatio);
    tB0 = 2;
    tBP = LifetimeRatio*tB0;
    
    B0Ratio = gRandom->Gaus(CenB0Ratio,SigB0Ratio);
    BchRatio = 1- B0Ratio;
    if(B0Ratio < 0. || B0Ratio > 1.) continue; 
    
    // scale the BR to the polarized ones
    // including Y factors
    BRpm *= FLong_pm*0.5/B0Ratio;
    BRpz *= FLong_pz*0.5/BchRatio;
    BRzz *= FLong_zz*0.5/B0Ratio;
    
    if(sol() == 0) continue;

    for(int j=0;j<8;j++){
      
      if(prior == 1) 
	xprob = 1.;
      else 
	xprob = J(alpha[j],  TPM[j],  T00[j],  P[j],  dP[j],  dT[j]);
      if(prior == 3) xprob *= P[j]*T00[j];

      // for alpha_SM fit
      // xprob *= exp(-0.5*pow((alpha[j]*180./M_PI-92.9)/5.7,2.));

      if(TPM[j] < 0 ||
	 T00[j] < 0 ||
	 P[j] < 0) continue;
      
      if(cutrange == 1 && 
	 (TPM[j] < minTPM || TPM[j] > maxTPM ||
	  T00[j] < minT00 || T00[j] > maxT00 ||
	  P[j]   < minP   || P[j]   > maxP)) continue;
      
      //use prior from Bs->KK in the weight
      if(SigPs != 0. && FlatPs != 0.) 
	xprob *= priorPs->GetBinContent(priorPs->FindBin(P[j]));

      // this is assuming that babar fits for
      // f(B) = 1 - C cos() + S sin()
      // with S = 2Im(lambda)/(1+|lamnda|^2)
      TComplex A = TComplex(T00[j]*cos(dT[j]-alpha[j]),T00[j]*sin(dT[j]-alpha[j])) -
        TComplex(P[j]*cos(dP[j]),P[j]*sin(dP[j]));
      TComplex Abarst = TComplex(T00[j]*cos(dT[j]+alpha[j]),-T00[j]*sin(dT[j]+alpha[j])) -
        TComplex(P[j]*cos(dP[j]),-P[j]*sin(dP[j]));
      TComplex SzzComplex = -2.*Abarst*A/(pow(A.Rho(),2.)+pow(Abarst.Rho(),2.));
      Szz = SzzComplex.Im();

      if(SigSzz2a != 0.) {
        xprob *= exp(-0.5*pow((Szz-CenSzz2a)/SigSzz2a,2.));
      }

      // impose |P|<2*|T+-|
      //      if(cutPoverT != 0. && TPM[j] != 0.)
      //	if(P[j]/TPM[j] > cutPoverT) continue;
	
      // Output Parameters
      histo1D[0]->Fill(sin(2*alpha[j]),xprob);
      histo1D[1]->Fill(alpha[j]*180./Pi,xprob);

      histo1D[2]->Fill(BRpm,xprob);
      histo1D[3]->Fill(BRpz,xprob);
      histo1D[4]->Fill(BRzz,xprob);
      histo1D[5]->Fill(C,xprob);
      histo1D[6]->Fill(S,xprob);
      histo1D[7]->Fill(Czz,xprob);
      histo1D[8]->Fill(Szz,xprob);
      histo1D[9]->Fill(P[j],xprob);
      histo1D[10]->Fill(TPM[j],xprob);
      histo1D[11]->Fill(T00[j],xprob);
      if (TPM[j]>0) 
	histo1D[12]->Fill(P[j]/TPM[j],xprob);
      histo1D[13]->Fill(dP[j]*180./Pi,xprob);
      histo1D[14]->Fill(dT[j]*180./Pi,xprob);
      //  histo1D[16]->Fill(deltaPb*180./Pi,xprob); E LA NUOVA FISICA???
      //  histo1D[17]->Fill(Pb,xprob);
      
      
      // ce serve alpha eff????
      //	if(1-pow(C,2.) !=0) {
      //	  double sin2alphaeff = S/sqrt(1-pow(C,2.));
      //	  double alphaeff = 0.5*asin(sin2alphaeff)*180./Pi;
      //	  if(alphaeff<0) alphaeff += 180.;
      //	  histo1D[18]->Fill(alpha*180./Pi-alphaeff,xprob);
      //	  histo1D[19]->Fill(alpha*180./Pi-alphaeff,xprob);
      //	}
      
      histo2D[1]->Fill(BRpz,BRpm,xprob);
      histo2D[2]->Fill(BRzz,BRpm,xprob);
      histo2D[3]->Fill(BRzz,BRpz,xprob);
      
      histo2D[4]->Fill(C,alpha[j]*180./Pi,xprob);
      if (TPM[j]>0) histo2D[5]->Fill(P[j]/TPM[j],alpha[j]*180./Pi,xprob);
      
      histo2D[7]->Fill(alpha[j]*180./Pi,TPM[j],xprob);
      histo2D[8]->Fill(alpha[j]*180./Pi,P[j],xprob);
      histo2D[9]->Fill(alpha[j]*180./Pi,T00[j],xprob);
      histo2D[10]->Fill(TPM[j],P[j],xprob);
      histo2D[11]->Fill(TPM[j],T00[j],xprob);
      histo2D[12]->Fill(TPM[j],dP[j]*180./Pi,xprob);
      histo2D[13]->Fill(TPM[j],dT[j]*180./Pi,xprob);
      histo2D[14]->Fill(P[j],T00[j],xprob);
      histo2D[15]->Fill(P[j],dP[j]*180./Pi,xprob);
      histo2D[16]->Fill(P[j],dT[j]*180./Pi,xprob);
      histo2D[17]->Fill(T00[j],dP[j]*180./Pi,xprob);
      histo2D[18]->Fill(T00[j],dT[j]*180./Pi,xprob);
      histo2D[19]->Fill(dP[j]*180./Pi,dT[j]*180./Pi,xprob);
      
      prior1D[0]->Fill(sin(2*alpha[j]),1.);
      prior1D[1]->Fill(alpha[j]*180./Pi,1.);
      prior1D[2]->Fill(BRpm,1.);
      prior1D[3]->Fill(BRpz,1.);
      prior1D[4]->Fill(BRzz,1.);
      prior1D[5]->Fill(C,1.);
      prior1D[6]->Fill(S,1.);
      prior1D[7]->Fill(Czz,1.);
      prior1D[8]->Fill(Szz,1.);
      prior1D[9]->Fill(P[j],1.);
      prior1D[10]->Fill(TPM[j],1.);
      prior1D[11]->Fill(T00[j],1.);
      prior1D[12]->Fill(dP[j]*180./Pi,1.);
      prior1D[13]->Fill(dT[j]*180./Pi,1.);
    }
  }
  
  for (int ix=1; ix<nrho+1; ix++) {
    for (int iy=1; iy<neta+1; iy++) {
      double irho = rho_min + (ix+0.5)*(rho_max-rho_min)/nrho;
      double ieta = eta_min + (iy+0.5)*(eta_max-eta_min)/neta;
      double ialpha = (acos(-1.)- atan2(ieta,1-irho) - atan2(ieta,irho))*180./Pi;
      if(ialpha < 0.) ialpha += 180.;
      if(ialpha > 180.) ialpha -= 180.;
      histo2D[0]->SetBinContent(ix,iy,histo1D[1]->GetBinContent(histo1D[2]->FindBin(ialpha)));
    }
  }

  TFile *output = new TFile(output_dot_root,"RECREATE");
  output->cd();
  output->mkdir("Priors","prior distributions");
  output->cd("Priors");
  for(int i=0; i<14; i++) {
    //normalize
    if (prior1D[i] && prior1D[i]->Integral() != 0) {
      double norm = 1./(prior1D[i]->Integral());
      prior1D[i]->Scale(norm);
      prior1D[i]->Write();
    }
  }

  if(SigPs != 0. && FlatPs != 0.) 
    priorPs->Write();

  output->cd();
  output->mkdir("Input","some input distributions");
  output->cd("Input");  
  for(int i=0; i<15; i++) {
    //normalize
    if (histo1D[i]) {
      double norm = 1./(histo1D[i]->Integral());
      histo1D[i]->Scale(norm);
      histo1D[i]->Write();
    }
  }
    
  for(int ind=0;ind<20;ind++) {
    if (histo2D[ind]) {
      double Sum = histo2D[ind]->GetSum();
      histo2D[ind]->Scale(1./Sum);
      histo2D[ind]->Write();
    }
  }
  output->Close();
}


void Read_Parameters(const char* filename) {
  char buffer[256];  
  ifstream reader (filename);
  map<string, double> data; 
  map<string, string> files;
  string label, equal, svalue;
  char cvalue[40];
  double value;

  if ( ! reader.is_open())
    { cout << "Error opening file"; exit (1); }

  while ( reader.getline (buffer,256) )
    {
      istringstream i(buffer);
      i >> label
	>> equal
	>> cvalue;
      if (equal == "=") {
	value = atof(cvalue);
	data[label] = value;
	//cout << "value: " << value << endl;
      } else if (equal == "==") {
	svalue = (string) cvalue;
	files[label] = svalue;
	//cout << "svalue: " << svalue << endl;
      }
      equal = "";
    }

  map<string, double>::const_iterator iter;
  for (iter=data.begin(); iter != data.end(); iter++) {
    cout << iter->second << " " << iter->first << endl;
  }

  Assign_Parameters(data);
  Assign_Parameters(files);

  return;
}

void Assign_Parameters(map<string, string> files) {
  skepC = string (files["skepC"]);
  skepS = string (files["skepS"]);
  channel = string (files["channel"]);
}

void Assign_Parameters(map<string, double> data) {

//     Parameters         
  NExtractions = int (data["NExtractions"]);
  Opt_NP   = int (data["Opt_NP"]);
  prior    = int (data["prior"]);
  cutrange = int (data["cutrange"]);

  cutPoverT = data["cutrange"];

  CenLtRatio = double (data["CenLtRatio"]);
  SigLtRatio = double (data["SigLtRatio"]);  
  
  CenB0Ratio = double (data["CenB0Ratio"]);
  SigB0Ratio = double (data["SigB0Ratio"]);  
  CenBchRatio = double (data["CenBchRatio"]);
  SigBchRatio = double (data["SigBchRatio"]);  
  
  CenS2a = double (data["CenS2a"]);
  SigS2a = double (data["SigS2a"]);
  FlatS2a = double (data["FlatS2a"]);

  SCcorr = double (data["SCcorr"]);
  
  CenC2a = double (data["CenC2a"]);
  SigC2a = double (data["SigC2a"]);
  FlatC2a = double (data["FlatC2a"]);
  
  CenCzz2a = double (data["CenCzz2a"]);
  SigCzz2a = double (data["SigCzz2a"]);
  FlatCzz2a = double (data["FlatCzz2a"]);

  CenSzz2a = double (data["CenSzz2a"]);
  SigSzz2a = double (data["SigSzz2a"]);
  FlatSzz2a = double (data["FlatSzz2a"]);
  
  CenBRpm = double (data["CenBRpm"]);
  SigBRpm = double (data["SigBRpm"]);
  CenFLpm = double (data["CenFLpm"]);
  SigFLpm = double (data["SigFLpm"]);
  
  CenBRpz = double (data["CenBRpz"]);
  SigBRpz = double (data["SigBRpz"]);
  CenFLpz = double (data["CenFLpz"]);
  SigFLpz = double (data["SigFLpz"]);
  
  CenBRzz = double (data["CenBRzz"]);
  SigBRzz = double (data["SigBRzz"]);
  CenFLzz = double (data["CenFLzz"]);
  SigFLzz = double (data["SigFLzz"]);

  minTPM = double (data["minTPM"]);
  maxTPM = double (data["maxTPM"]);
  minT00 = double (data["minT00"]);
  maxT00 = double (data["maxT00"]);
  minP   = double (data["minP"]);
  maxP   = double (data["maxP"]);

  CenPs   = double (data["CenPs"]);
  SigPs   = double (data["SigPs"]);
  FlatPs  = double (data["FlatPs"]);
   
}

double J(double alp, double tPM, double t00, double p, double dp, double dt) {

  double res=(-32.*p*t00*tPM*(pow(tPM, 2.)*cos(2.*alp) + p*(p + 2.*tPM*cos(alp)*cos(dp)))*
	      (pow(t00, 2.)*pow(tPM, 2.) +
	       pow(p, 2.)*(pow(t00, 2.) + pow(tPM, 2.))*cos(2.*alp) +
	       pow(p, 2.)*t00*tPM*cos(2.*alp - dt) -
	       pow(t00, 2.)*pow(tPM, 2.)*cos(2.*dt) +
	       pow(p, 2.)*t00*tPM*cos(2.*alp + dt) -
	       2.*pow(p, 2.)*t00*tPM*cos(dt - 2.*dp) +
	       p*pow(t00, 2.)*tPM*cos(alp - dp) +
	       p*t00*pow(tPM, 2.)*cos(alp - dt - dp) -
	       pow(p, 2.)*pow(t00, 2.)*cos(2.*(dt - dp)) -
	       p*t00*pow(tPM, 2.)*cos(alp + dt - dp) -
	       p*pow(t00, 2.)*tPM*cos(alp + 2.*dt - dp) -
	       pow(p, 2.)*pow(tPM, 2.)*cos(2.*dp) +
	       p*pow(t00, 2.)*tPM*cos(alp + dp) -
	       p*pow(t00, 2.)*tPM*cos(alp - 2.*dt + dp) -
	       p*t00*pow(tPM, 2.)*cos(alp - dt + dp) +
	       p*t00*pow(tPM, 2.)*cos(alp + dt + dp))*
	      pow(sin(alp), 2.))/
    ((pow(p, 2.) + pow(t00, 2.) - 2.*p*t00*cos(alp)*cos(dt - dp))*
     pow(pow(p, 2.) + pow(tPM, 2.) + 2.*p*tPM*cos(alp)*cos(dp), 2.));
  return(1./fabs(res));
}

int sol(){
  
  double c = pow(LifetimeRatio,0.5)*
    (BRpz/LifetimeRatio+BRpm*(1.+C)/2.-BRzz*(1.+Czz))/(pow(2.*BRpm*BRpz*(1.+C),0.5));
  if (fabs(c)>1.) return(0);
  double cbar = pow(LifetimeRatio,0.5)*
    (BRpz/LifetimeRatio+BRpm*(1.-C)/2.-BRzz*(1.-Czz))/(pow(2.*BRpm*BRpz*(1.-C),0.5));
  if (fabs(cbar)>1) return(0);
  
  double sin2a = S/pow(1.-C*C,0.5);
  double cos2a = pow(1.-sin2a*sin2a,0.5);
  
  double r;
  double s = pow(1-c*c,0.5);
  double sbar = pow(1-cbar*cbar,0.5);

  r = cos2a*cbar-sin2a*sbar+c;
  if(r==0.) alpha[0] = M_PI/2.; else alpha[0]=atan((sin2a*cbar+cos2a*sbar+s)/r);
  if(r==0.) alpha[1] = M_PI/2.; else alpha[1]=atan((sin2a*cbar+cos2a*sbar-s)/r);

  r = cos2a*cbar+sin2a*sbar+c;
  if(r==0.) alpha[2] = M_PI/2.; else alpha[2]=atan((sin2a*cbar-cos2a*sbar+s)/r);
  if(r==0.) alpha[3] = M_PI/2.; else alpha[3]=atan((sin2a*cbar-cos2a*sbar-s)/r);
  
  r = -cos2a*cbar-sin2a*sbar+c;
  if(r==0.) alpha[4] = M_PI/2.; else alpha[4]=atan((sin2a*cbar-cos2a*sbar+s)/r);
  if(r==0.) alpha[5] = M_PI/2.; else alpha[5]=atan((sin2a*cbar-cos2a*sbar-s)/r);

  r = -cos2a*cbar+sin2a*sbar+c;
  if(r==0.) alpha[6] = M_PI/2.; else alpha[6]=atan((sin2a*cbar+cos2a*sbar+s)/r);
  if(r==0.) alpha[7] = M_PI/2.; else alpha[7]=atan((sin2a*cbar+cos2a*sbar-s)/r);
  
  for(int al=0;al<8;al++){
    TPM[al] = -1.;
    T00[al] = -1.;
    P[al]   = -1.;

    if(alpha[al]<0) alpha[al]=alpha[al]+M_PI;
    double xtpm = 0.;
    double xp   = 0.;
    double xt   = 0.;
    double xdp  = 0.;
    double xd0  = 0.;
    int    nsol = 0;

    for(int j1=-1; j1<2; j1=j1+2){
      if(nsol>0) break;
      xtpm = pow(BRpm/pow(sin(alpha[al]),2.)*(1.+((double) j1)*pow(1.-C*C-S*S,0.5))/tB0,0.5);
      xp   = xtpm*xtpm*(2.*pow(cos(alpha[al]),2.)-1.)+2.*BRpm/tB0*(1.-S/tan(alpha[al]));
      if(xp<0.) continue;
      xp=pow(xp,0.5);
      xdp=atan2(-BRpm*C/(sin(alpha[al])*xp*xtpm*tB0),
		-((xp*xp+xtpm*xtpm)/2.-BRpm/tB0)/(cos(alpha[al])*xp*xtpm));
      if(xdp<0.) xdp=xdp+2.*M_PI;
      for(int j2 = -1; j2<2;j2=j2+2){
	if(nsol>0) break;
	xt=pow(xp,4.)-16.*pow(BRzz*Czz/tB0/sin(2.*alpha[al]),2.)+pow(xp/cos(alpha[al]),2.)*(4.*BRzz/tB0-xp*xp);
	if(xt<0.) continue;
	xt=xp*xp*cos(2.*alpha[al])+4.*BRzz/tB0+2.*((double) j2)*pow(cos(alpha[al]),2.)*pow(xt,0.5);
	if(xt<0.) continue;
	xt=pow(xt,0.5);
	xd0=xdp+atan2(-2*BRzz*Czz/(sin(alpha[al])*xp*xt*tB0),
		      ((xp*xp+xt*xt)/2.-2.*BRzz/tB0)/(cos(alpha[al])*xp*xt));
	if(xd0<0.) xd0=xd0+2.*M_PI;
	if(xd0>2.*M_PI) xd0=xd0-2.*M_PI;
	double bp0 = tBP/4.*(xt*xt+xtpm*xtpm+2.*xt*xtpm*cos(xd0));
	if(fabs(bp0-BRpz)>.000001) continue;
	TPM[al] = xtpm;
	T00[al] = xt;
	P[al]   = xp;
	dP[al]  = xdp;
	dT[al]  = xd0;
	//	cout << al << " " << j1 << " " << j2 << endl;
	nsol++;
      }
    }
    
  }
  return(1);
}   

//////////////////////////////////////////////////////////////
//  convolution between a Gaussian and a flat distribution  //
//////////////////////////////////////////////////////////////

Double_t Hfun(Double_t *x, Double_t *par) {

  double xx =x[0];
  return xfun(xx,par[0],par[1],par[2]);
}

double xfun(double x, double a, double sig, double flat)
{
  double fun;
  fun = fun_conv(x,a,sig,flat);
  return fun;
}

double fun_conv(double x, double a, double sig, double flat)
{
  double   xup,xlow,range,cc1,cc2,rde,fun;

  if (flat == 0) {
    fun  =  exp(-0.5*(x-a)*(x-a)/(sig*sig));
  } else {
    xup  = a+flat;
    xlow = a-flat;
    range = xup-xlow;
    cc1   = (x-xlow)/(sqrt(2.)*sig);
    cc2   = (x-xup)/(sqrt(2.)*sig);
    rde   = 1/(2*range);
    fun  = rde*(erf(cc1)-erf(cc2)) ;
  }

  return fun;
}

