#ifndef CKM_gen_rhopi_HH
#define CKM_gen_rhopi_HH

#include <iostream>
#include <string>
#include <vector>

using namespace std;

double xcen, xsig,xfla;
double xfun(double x, double a, double sig, double flat);
double fun_hock(double var, double cen, double gaus, double flat);
double fun_conv(double x, double a, double sig, double flat);
double Get_Val(double mycen, double mysig, double myfla);
double F_Sin2bpg(double);
Double_t myGauss(Double_t *x, Double_t *par);
void Fill_Strategies(double , double, double[10]);
void Read_Parameters(const char* );
void Assign_Parameters(map<string, double>);

double Ftt;

int NOpt, NChan;
double Rb,Rt;
vector< vector <int> > Opt;

vector<string> chan;
int NExtractions;
int Opt_NP;

int UseBabar;
int UseBelle;

double Apm_max;
double Amp_max;
double Azz_max;

// separated correlation matrices for BaBar
TMatrixF corStat_BaBar(26,26);
TMatrixF corSys_BaBar(26,26);
// single correlation matrix for Belle
TMatrixF cor_Belle(26,26);

// central values of parameters and total errors
// obtained by summing in quadrature the stat and syst
double CenNeu_BaBar[26], SigNeu_BaBar[26], SysNeu_BaBar[26]; //SigTot_BaBar[26];
double CenNeu_Belle[26], SigTot_Belle[26];
double xneu[26];


#endif 

