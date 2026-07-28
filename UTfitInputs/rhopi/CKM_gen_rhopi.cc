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
#include "CKM_gen_rhopi.hh"

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

  TString conffile("CKM_gen_rhopi");
  TString outfile("CKM_gen_rhopi");
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

  
  // 1D Histos
  double alpha;
  TComplex Apm, Amp, Azz;

  double Upm, Ump, Uzz;
  double Vpm, Vmp, Vzz;

  double Re_Upm, Re_Upz, Re_Umz;
  double Re_Vpm, Re_Vpz, Re_Vmz;

  double Im_Upm, Im_Upz, Im_Umz;
  double Im_Vpm, Im_Vpz, Im_Vmz;

  double Ipm, Imp, Izz;

  double Re_Ipm, Re_Ipz, Re_Imz;
  double Im_Ipm, Im_Ipz, Im_Imz;

  TH1D *histo1D[12];
  //for (int h=0; h<12; h++) {
  // histo1D[h]=0;
  //}

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

  // seed initialization 
  cout << "=========> number of extractions " << NExtractions << endl;
  
  TMatrixDSym Vneu_BaBar(26);
  TMatrixDSym Vneu_Belle(26);
  Vneu_BaBar = (covStat_BaBar+covSys_BaBar).Invert();
  //Vneu_BaBar = covStat_BaBar.Invert();
  Vneu_Belle = cov_Belle.Invert();

  for(int ix=0;ix<26;ix++) {
    for(int iy=0;iy<26;iy++) {
      cout << "BaBar " << Vneu_BaBar(ix,iy) << endl;
    }
  }

  for(int ix=0;ix<26;ix++) {
    for(int iy=0;iy<26;iy++) {
      cout << "Belle " << Vneu_Belle(ix,iy) << endl;
    }
  }
  
  int jobpid = getpid();
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int seed = today+clock+jobpid;
  TRandom3 *gRandom = new TRandom3(seed);
  cout << " Seed = " << seed << endl;  
  
  for(int k=1; k<=NExtractions; k++) {
    if(fmod(k,1000000.) == 0) cout << k << " events generated" << endl;

     // one arbitrary abs value and phase
    TComplex Apm(1., 0., kTRUE);
    TComplex Amp(gRandom->Rndm()*Amp_max, (2.*gRandom->Rndm()-1.)*Pi, kTRUE);
    TComplex Azz(gRandom->Rndm()*Azz_max, (2.*gRandom->Rndm()-1.)*Pi, kTRUE);

    TComplex Apm_bar(gRandom->Rndm()*Apm_max, (2.*gRandom->Rndm()-1.)*Pi, kTRUE);
    TComplex Amp_bar(gRandom->Rndm()*Amp_max, (2.*gRandom->Rndm()-1.)*Pi, kTRUE);
    TComplex Azz_bar(gRandom->Rndm()*Azz_max, (2.*gRandom->Rndm()-1.)*Pi, kTRUE);

    TComplex e21alpha = (Apm_bar+Amp_bar+2.*Azz_bar)/(Apm+Amp+2.*Azz);

    alpha = e21alpha.Theta()/2.;
    if(alpha<0.) alpha += Pi;

    
    Upm = pow(Apm.Rho(),2.) + pow(Apm_bar.Rho(),2.);
    Ump = (pow(Amp.Rho(),2.) + pow(Amp_bar.Rho(),2.))/Upm;
    Uzz = (pow(Azz.Rho(),2.) + pow(Azz_bar.Rho(),2.))/Upm;
    
    Vpm = (pow(Apm.Rho(),2.) - pow(Apm_bar.Rho(),2.))/Upm;
    Vmp = (pow(Amp.Rho(),2.) - pow(Amp_bar.Rho(),2.))/Upm;
    Vzz = (pow(Azz.Rho(),2.) - pow(Azz_bar.Rho(),2.))/Upm;
    
    Re_Upm = ((Apm.Re()*Amp.Re() + Apm.Im()*Amp.Im()) + (Apm_bar.Re()*Amp_bar.Re() + Apm_bar.Im()*Amp_bar.Im()))/Upm;
    Re_Upz = ((Apm.Re()*Azz.Re() + Apm.Im()*Azz.Im()) + (Apm_bar.Re()*Azz_bar.Re() + Apm_bar.Im()*Azz_bar.Im()))/Upm;
    Re_Umz = ((Amp.Re()*Azz.Re() + Amp.Im()*Azz.Im()) + (Amp_bar.Re()*Azz_bar.Re() + Amp_bar.Im()*Azz_bar.Im()))/Upm;
    
    Im_Upm = ((Apm.Im()*Amp.Re() - Apm.Re()*Amp.Im()) + (Apm_bar.Im()*Amp_bar.Re() - Apm_bar.Re()*Amp_bar.Im()))/Upm;
    Im_Upz = ((Apm.Im()*Azz.Re() - Apm.Re()*Azz.Im()) + (Apm_bar.Im()*Azz_bar.Re() - Apm_bar.Re()*Azz_bar.Im()))/Upm;
    Im_Umz = ((Amp.Im()*Azz.Re() - Amp.Re()*Azz.Im()) + (Amp_bar.Im()*Azz_bar.Re() - Amp_bar.Re()*Azz_bar.Im()))/Upm;
    
    Re_Vpm = ((Apm.Re()*Amp.Re() + Apm.Im()*Amp.Im()) - (Apm_bar.Re()*Amp_bar.Re() + Apm_bar.Im()*Amp_bar.Im()))/Upm;
    Re_Vpz = ((Apm.Re()*Azz.Re() + Apm.Im()*Azz.Im()) - (Apm_bar.Re()*Azz_bar.Re() + Apm_bar.Im()*Azz_bar.Im()))/Upm;
    Re_Vmz = ((Amp.Re()*Azz.Re() + Amp.Im()*Azz.Im()) - (Amp_bar.Re()*Azz_bar.Re() + Amp_bar.Im()*Azz_bar.Im()))/Upm;
    
    Im_Vpm = ((Apm.Im()*Amp.Re() - Apm.Re()*Amp.Im()) - (Apm_bar.Im()*Amp_bar.Re() - Apm_bar.Re()*Amp_bar.Im()))/Upm;
    Im_Vpz = ((Apm.Im()*Azz.Re() - Apm.Re()*Azz.Im()) - (Apm_bar.Im()*Azz_bar.Re() - Apm_bar.Re()*Azz_bar.Im()))/Upm;
    Im_Vmz = ((Amp.Im()*Azz.Re() - Amp.Re()*Azz.Im()) - (Amp_bar.Im()*Azz_bar.Re() - Amp_bar.Re()*Azz_bar.Im()))/Upm;
    
    Ipm = (Apm_bar.Im()*Apm.Re() - Apm_bar.Re()*Apm.Im())/Upm;
    Imp = (Amp_bar.Im()*Amp.Re() - Amp_bar.Re()*Amp.Im())/Upm;
    Izz = (Azz_bar.Im()*Azz.Re() - Azz_bar.Re()*Azz.Im())/Upm;
    
    Re_Ipm = ((Apm_bar.Re()*Amp.Re() + Apm_bar.Im()*Amp.Im()) - (Amp_bar.Re()*Apm.Re() + Amp_bar.Im()*Apm.Im()))/Upm;
    Re_Ipz = ((Apm_bar.Re()*Azz.Re() + Apm_bar.Im()*Azz.Im()) - (Azz_bar.Re()*Apm.Re() + Azz_bar.Im()*Apm.Im()))/Upm;
    Re_Imz = ((Amp_bar.Re()*Azz.Re() + Amp_bar.Im()*Azz.Im()) - (Azz_bar.Re()*Amp.Re() + Azz_bar.Im()*Amp.Im()))/Upm;
    
    Im_Ipm = ((Apm_bar.Im()*Amp.Re() - Apm_bar.Re()*Amp.Im()) + (Amp_bar.Im()*Apm.Re() - Amp_bar.Re()*Apm.Im()))/Upm;
    Im_Ipz = ((Apm_bar.Im()*Azz.Re() - Apm_bar.Re()*Azz.Im()) + (Azz_bar.Im()*Apm.Re() - Azz_bar.Re()*Apm.Im()))/Upm;
    Im_Imz = ((Amp_bar.Im()*Azz.Re() - Amp_bar.Re()*Azz.Im()) + (Azz_bar.Im()*Amp.Re() - Azz_bar.Re()*Amp.Im()))/Upm;

    xneu[0] = Izz;
    xneu[1] = Imp;
    xneu[2] = Im_Imz;
    xneu[3] = Re_Imz;
    xneu[4] = Ipm;
    xneu[5] = Im_Ipz;
    xneu[6] = Re_Ipz;
    xneu[7] = Im_Ipm;
    xneu[8] = Re_Ipm;
    xneu[9] = Vzz;
    xneu[10] = Uzz;
    xneu[11] = Im_Vmz;
    xneu[12] = Re_Vmz;
    xneu[13] = Im_Umz;
    xneu[14] = Re_Umz;
    xneu[15] = Vmp;
    xneu[16] = Ump;
    xneu[17] = Im_Vpz;
    xneu[18] = Re_Vpz;
    xneu[19] = Im_Upz;
    xneu[20] = Re_Upz;
    xneu[21] = Im_Vpm;
    xneu[22] = Re_Vpm;
    xneu[23] = Im_Upm;
    xneu[24] = Re_Upm;
    xneu[25] = Vpm;
    
    double xprob = 1.;
    
    double chisq = 0.;
    if(UseBabar == 1){
      for(int ix=0;ix<26;ix++) {
	for(int iy=0;iy<26;iy++) {
	  chisq += (xneu[ix]-CenNeu_BaBar[ix])*Vneu_BaBar(ix,iy)*(xneu[iy]-CenNeu_BaBar[iy]);
	}
      }
      xprob *= exp(-0.5*chisq);
    }
    
    if(UseBelle == 1){
      chisq = 0.;
      for(int ix=0;ix<26;ix++) {
	for(int iy=0;iy<26;iy++) {
	  chisq += (xneu[ix]-CenNeu_Belle[ix])*Vneu_Belle(ix,iy)*(xneu[iy]-CenNeu_Belle[iy]);
	}
      }
      xprob *= exp(-0.5*chisq);
    }
    // Output Parameters
    histo1D[0]->Fill(sin(2.*alpha),xprob);
    histo1D[1]->Fill(alpha*180./Pi,xprob);

    histo1D[2]->Fill(Amp.Rho(),xprob);
    histo1D[3]->Fill(Azz.Rho(),xprob);
    histo1D[4]->Fill(Apm_bar.Rho(),xprob);
    histo1D[5]->Fill(Amp_bar.Rho(),xprob);
    histo1D[6]->Fill(Azz_bar.Rho(),xprob);
    histo1D[7]->Fill(Amp.Theta()*180./Pi,xprob);
    histo1D[8]->Fill(Azz.Theta()*180./Pi,xprob);
    histo1D[9]->Fill(Apm_bar.Theta()*180./Pi,xprob);
    histo1D[10]->Fill(Amp_bar.Theta()*180./Pi,xprob);
    histo1D[11]->Fill(Azz_bar.Theta()*180./Pi,xprob);

  }

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

  /*
  SigTot_BaBar[0] = sqrt(pow(data[("SigIzz_BaBar")],2.)+pow(data[("SysIzz_BaBar")],2.));
  SigTot_BaBar[1] = sqrt(pow(data[("SigImp_BaBar")],2.)+pow(data[("SysImp_BaBar")],2.));
  SigTot_BaBar[2] = sqrt(pow(data[("SigIm_Imz_BaBar")],2.)+pow(data[("SysIm_Imz_BaBar")],2.));
  SigTot_BaBar[3] = sqrt(pow(data[("SigRe_Imz_BaBar")],2.)+pow(data[("SysRe_Imz_BaBar")],2.));
  SigTot_BaBar[4] = sqrt(pow(data[("SigIpm_BaBar")],2.)+pow(data[("SysIpm_BaBar")],2.));
  SigTot_BaBar[5] = sqrt(pow(data[("SigIm_Ipz_BaBar")],2.)+pow(data[("SysIm_Ipz_BaBar")],2.));
  SigTot_BaBar[6] = sqrt(pow(data[("SigRe_Ipz_BaBar")],2.)+pow(data[("SysRe_Ipz_BaBar")],2.));
  SigTot_BaBar[7] = sqrt(pow(data[("SigIm_Ipm_BaBar")],2.)+pow(data[("SysIm_Ipm_BaBar")],2.));
  SigTot_BaBar[8] = sqrt(pow(data[("SigRe_Ipm_BaBar")],2.)+pow(data[("SysRe_Ipm_BaBar")],2.));
  SigTot_BaBar[9] = sqrt(pow(data[("SigVzz_BaBar")],2.)+pow(data[("SysVzz_BaBar")],2.));
  SigTot_BaBar[10] = sqrt(pow(data[("SigUzz_BaBar")],2.)+pow(data[("SysUzz_BaBar")],2.));
  SigTot_BaBar[11] = sqrt(pow(data[("SigIm_Vmz_BaBar")],2.)+pow(data[("SysIm_Vmz_BaBar")],2.));
  SigTot_BaBar[12] = sqrt(pow(data[("SigRe_Vmz_BaBar")],2.)+pow(data[("SysRe_Vmz_BaBar")],2.));
  SigTot_BaBar[13] = sqrt(pow(data[("SigIm_Umz_BaBar")],2.)+pow(data[("SysIm_Umz_BaBar")],2.));
  SigTot_BaBar[14] = sqrt(pow(data[("SigRe_Umz_BaBar")],2.)+pow(data[("SysRe_Umz_BaBar")],2.));
  SigTot_BaBar[15] = sqrt(pow(data[("SigVmp_BaBar")],2.)+pow(data[("SysVmp_BaBar")],2.));
  SigTot_BaBar[16] = sqrt(pow(data[("SigUmp_BaBar")],2.)+pow(data[("SysUmp_BaBar")],2.));
  SigTot_BaBar[17] = sqrt(pow(data[("SigIm_Vpz_BaBar")],2.)+pow(data[("SysIm_Vpz_BaBar")],2.));
  SigTot_BaBar[18] = sqrt(pow(data[("SigRe_Vpz_BaBar")],2.)+pow(data[("SysRe_Vpz_BaBar")],2.));
  SigTot_BaBar[19] = sqrt(pow(data[("SigIm_Upz_BaBar")],2.)+pow(data[("SysIm_Upz_BaBar")],2.));
  SigTot_BaBar[20] = sqrt(pow(data[("SigRe_Upz_BaBar")],2.)+pow(data[("SysRe_Upz_BaBar")],2.));
  SigTot_BaBar[21] = sqrt(pow(data[("SigIm_Vpm_BaBar")],2.)+pow(data[("SysIm_Vpm_BaBar")],2.));
  SigTot_BaBar[22] = sqrt(pow(data[("SigRe_Vpm_BaBar")],2.)+pow(data[("SysRe_Vpm_BaBar")],2.));
  SigTot_BaBar[23] = sqrt(pow(data[("SigIm_Upm_BaBar")],2.)+pow(data[("SysIm_Upm_BaBar")],2.));
  SigTot_BaBar[24] = sqrt(pow(data[("SigRe_Upm_BaBar")],2.)+pow(data[("SysRe_Upm_BaBar")],2.));
  SigTot_BaBar[25] = sqrt(pow(data[("SigVpm_BaBar")],2.)+pow(data[("SysVpm_BaBar")],2.));
  */

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
