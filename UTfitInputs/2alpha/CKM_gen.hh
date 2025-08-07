#ifndef CKM_gen_HH
#define CKM_gen_HH

#include <string>
#include <vector>

using namespace std;

double xcen, xsig,xfla;
double xfun(double x, double a, double sig, double flat);
double J(double alp, double tPM, double t00, 
	 double p, double dp, double dt);
void Read_Parameters(const char* );
void Assign_Parameters(map<string, double>);
void Assign_Parameters(map<string, string>);
int sol();

Double_t Hfun(Double_t *x, Double_t *par);
double xfun(double x, double a, double sig, double flat);
double fun_conv(double x, double a, double sig, double flat);

vector< vector <int> > Opt;
vector<string> chan;
vector<string> fit;

double CenLtRatio, SigLtRatio;
double CenB0Ratio, SigB0Ratio;
double CenBchRatio, SigBchRatio;
int NExtractions;
int Opt_NP;

double CenBRpm, SigBRpm;
double CenFLpm, SigFLpm;
double CenBRpz, SigBRpz;
double CenFLpz, SigFLpz;
double CenBRzz, SigBRzz;
double CenFLzz, SigFLzz;

double CenS2a, SigS2a, FlatS2a;
double CenC2a, SigC2a, FlatC2a;
double SCcorr;
double CenCzz2a, SigCzz2a, FlatCzz2a;
double CenSzz2a, SigSzz2a, FlatSzz2a;

double minP, maxP;
double minTPM, maxTPM;
double minT00, maxT00;
double CenPs, SigPs, FlatPs;

int prior;
int cutrange;
string skepC, skepS;
string channel;

double B0Ratio, BchRatio;
double cutPoverT;

double alpha[8];
double TPM[8]; 
double T00[8]; 
double P[8]; 
double dP[8]; 
double dT[8]; 
double xprob;

double LifetimeRatio, tBP, tB0;
double FLong_pm, FLong_pz, FLong_zz; 
double BRpm, BRpz, BRzz;
double C, S, Czz, Szz;

#endif 

