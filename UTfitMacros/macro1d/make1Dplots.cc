#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TImage.h>
#include <TColor.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TText.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include "FuncFromHisto.cc"
#define ncoo       50
#define nvalrebin   4

using namespace std;

TString fit1,fit2,fit3,dir1,dir2,dir3,season,addtxt="";
bool lsup1=false,lsup2=false,lsup3=false;
int colsup=4;

class hdata{
public:
  hdata()
         {m_index=0; m_cont=0;};
  hdata(int index,double cont)
         {m_index=index; m_cont=cont;};
  bool   operator<(const hdata &x)
         {return m_cont>x.m_cont;};
  int    Index()
         {return m_index;};
  double Cont()
         {return m_cont;};
private:
  int    m_index;
  double m_cont;
};

TString xlab,ylab;
bool lxlab=false,lylab=false;

TApplication   *theapp=NULL;
TFile          *datafile=NULL;
TCanvas        *tela=NULL;
TH1D           *PlotToBeDrawn,*Plot68,*Plot95;
TString        titletmp="none";
TString      lname="none";

bool   lgraphic=true,lmakeps=false,lprecise=false,lsigma=false,lmediana=false,ldebug=false,llogo=true;
int    nlevel=3;
int    rebin=0;
int    smooth=0;
int    col1=2;
int    col2=5;
int    val[nvalrebin]={1,2,5,10};
// 1,2,3,3.5,4 sigma
double valref[5]={68.27,95.45,99.73,99.95,99.99};
// 2.2,2.4,2.6,2.8,3 sigma
//double valref[5]={97.22,98.36,99.07,99.49,99.73};
// 2.1,2.2,2.3,2.4,2.5 sigma
//double valref[5]={96.43,97.22,97.86,98.36,98.76};
// 2.3,2.4,2.5,2.6,2.7 sigma
//double valref[5]={97.86,98.36,98.76,99.07,99.31};
// 2.6,2.7,2.8,2.9,3.0 sigma
//double valref[5]={99.07,99.31,99.49,99.63,99.73};
// 2.6,2.65,2.7,2.75,3 sigma
//double valref[5]={99.07,99.20,99.31,99.40,99.49};
// 1.5, 1.6, 1.7, 1.8, 1.9 sigma
//double valref[5]={86.64,89.04,91.09,92.81,94.26};
// 1.1, 1.2, 1.3, 1.4, 1.5 sigma
//double valref[5]={72.87,76.99,80.64,83.85,86.64};

void DrawImage(double xl, double yl, double xu, double yu){
  TImage *img = TImage::Open("../common/logo.png");
  TPad *pad = new TPad("pad","pad",xl,yl,xu,yu,0,0,0);
  img->SetConstRatio(kTRUE);
  pad->SetFillStyle(4000);
  pad->Draw();
  pad->cd();
  img->Draw();
}

void PutLogo(){

  double yl=0.77,yu=0.88, xl=0.21, xu=0.39;  

  
  double xc;
  int    nbin = PlotToBeDrawn->GetNbinsX();
  double hmax = PlotToBeDrawn->GetBinContent(1);
  double hmin = PlotToBeDrawn->GetBinContent(1);

  for (int i=0;i<nbin;i++){
    if (PlotToBeDrawn->GetBinContent(i+1)>hmax){
      hmax = PlotToBeDrawn->GetBinContent(i+1);
      xc   = PlotToBeDrawn->GetBinCenter(i+1);
    }
    if (PlotToBeDrawn->GetBinContent(i+1)<hmin) hmin = PlotToBeDrawn->GetBinContent(1);
  }

  double xmax = PlotToBeDrawn->GetXaxis()->GetXmax();
  double xmin = PlotToBeDrawn->GetXaxis()->GetXmin();
  
  bool left=true,up=false;   
  if (xc<(xmax+xmin)/2) left=false;
  if (fabs(hmax-hmin)/hmax<0.2) up = true;
  
  if (!up && !left){
    xl=.70; xu=.88; 
    //xl=.71; xu=.89; 
  }
  if (up){
    PlotToBeDrawn->SetMaximum(hmax/0.75);
  }

  DrawImage(xl,yl,xu,yu);

  TFile *output = new TFile("dump.root","RECREATE");
  output->cd();
  output->mkdir("dump");
  output->cd("dump");
  //PlotToBeDrawn->SetName("input_alpha_all");
  //PlotToBeDrawn->SetTitle("input_alpha_all");
  PlotToBeDrawn->Write();
  output->Close();

  gPad->Update();
  tela->cd();

  double newyl=yl-(yu-yl);
  double newxl=xl-(xu-xl)/2;
  double newyu=yu;
  double newxu=xu+(xu-xl)/2;

  TPad *padSeason = new TPad("padSeason","padSeason",newxl,newyl,newxu,newyu,0,0,0);

  padSeason->SetFillStyle(4000);
  padSeason->Draw();
  padSeason->cd();
  TLatex* tex = new TLatex;
  //TColor *col = new TColor(999,42/255.,44/255.,114/255.);
  tex->SetTextSize(0.16);
  tex->SetTextAlign(21);
  tex->SetTextColor(602);
  tex->DrawTextNDC(0.50,0.43, season);

  if (addtxt!=""){
    TLatex* addtex = new TLatex;
    addtex->SetTextSize(0.18);
    addtex->SetTextAlign(21);
    if (!up && !left){
      addtex->DrawTextNDC(0.38,0.1, addtxt);
    } else {
      addtex->DrawTextNDC(0.57,0.1, addtxt);
    }
  }

}

void SetXTitle(TH1D *id, TString tit){

  TString str;

  if (titletmp!="none") tit = titletmp;

  if      (strncmp(tit,"sin2bpg",7)==0)    str = "sin(2#beta+#gamma)";
  else if (strncmp(tit,"sin2b",5)==0)      str = "sin2#beta";
  else if (strncmp(tit,"sin2a",5)==0)      str = "sin2#alpha";
  else if (strncmp(tit,"gamma",5)==0)      str = "#gamma[^{o}]";
  else if (strncmp(tit,"beta",4)==0)       str = "#beta[^{o}]";
  else if (strncmp(tit,"alpha",5)==0)      str = "#alpha[^{o}]";
  else if (strncmp(tit,"rhof",4)==0)       str = "#bar{#rho}";
  else if (strncmp(tit,"etaf",4)==0)       str = "#bar{#eta}";
  else if (strncmp(tit,"Lambda",6)==0)     str = "#lambda";
  else if (strncmp(tit,"Imlamt",6)==0)     str = "Im#lambda_{t}[10^{-5}]";
  else if (strncmp(tit,"Relamt",6)==0)     str = "Re#lambda_{t}[10^{-3}]";
  else if (strncmp(tit,"Rb",2)==0)         str = "R_{b}";
  else if (strncmp(tit,"Rt",2)==0)         str = "R_{t}";
  else if (strncmp(tit,"Vub",3)==0)        str = "V_{ub}";
  else if (strncmp(tit,"Vcb",3)==0)        str = "V_{cb}";
  else if (strncmp(tit,"Vtd",3)==0)        str = "V_{td}#times10^{3}";
  else if (strncmp(tit,"DeltaMs",7)==0)    str = "#Delta m_{s}[ps^{-1}]";
  else if (strncmp(tit,"FbSqrtb",7)==0)    str = "F_{B_{s}}#sqrt{B_{s}}";
  else if (strncmp(tit,"Bk",2)==0)         str = "B_{k}";
  else if (strncmp(tit,"Xi",2)==0)         str = "#xi";
  else if (strncmp(tit,"Ftt",3)==0)        str = "F_{tt}";
  else                           str = tit;
  
  if (!lxlab){
    id->GetXaxis()->SetTitle(str);
  } else {
    id->GetXaxis()->SetTitle(xlab);
  }

}


int EvaluateSignChanges(TH1D *his){

  int nbin = his->GetNbinsX();
  double *vect = new double[nbin];
  for (int i=0;i<nbin;i++){
    vect[i] = his->GetBinContent(i+1);
  }
 
  double xmax=vect[0];
  double xmin=vect[0];
  for (int i=0;i<nbin;i++){
    if (vect[i]>=xmax) xmax=vect[i];
    if (vect[i]<=xmin) xmin=vect[i];
  }

  int n=0;
  if ((xmax-xmin)/xmax>=0.01){
    for (int i=0;i<=nbin-3;i++){
      if (vect[i]>xmax/10.0){
	if ((vect[i+2]-vect[i+1])*(vect[i+1]-vect[i])<0) n++;
      }
    }
  }

  return n;

}

int Do_Cont_Mediana(const TString plot, TH1D* Plot, double *GArray, int nbin, double *NArray, int nrebin){

  double *chk,*xcoo;
  int *icls;
  //TMath m;

  // delete all the underflow and overflow, otherwise fucked root 
  // take them into account in the normalization
  Plot->SetBinContent(0,0.);
  Plot->SetBinContent(1+Plot->GetNbinsX(),0.);

  chk  = new double[nlevel*2+1];
  xcoo = new double[nlevel*2+1];
  icls = new int[nlevel*2+1];

  for (int i=0;i<nlevel;i++){
    chk[i] = (1.-valref[nlevel-1-i]/100.0)/2;
    if (lsigma) chk[i] = TMath::Erfc((nlevel-i)/sqrt(2.))/2;
  }
  chk[nlevel] = 0.5;
  for (int i=0;i<nlevel;i++){
    chk[i+1+nlevel] = 1.-(1.-valref[i]/100.0)/2;
    if (lsigma) chk[i+1+nlevel] = 1.-TMath::Erfc((i+1)/sqrt(2.))/2;
  }  

  double low = Plot->GetXaxis()->GetXmin();
  double up  = Plot->GetXaxis()->GetXmax();

  double norm=pow(10.,floor(log10(up))-3);

  int newbin = (int) ceil((up-low)/norm);

  if (ldebug) cout << "up " << up << " low " << low << " norm " << norm << endl;
  if (ldebug) cout << "nrebin " << nrebin << " newbin " << newbin << endl;
  if (nrebin<newbin){
    nrebin = newbin;

    delete[] NArray;
    NArray = new double[nrebin];

    //    TF1 *fh;
    //fh = Plot->GetFunction("xx");
    double sum=0;
    for (int i=0;i<nrebin;i++){
      double xc = low+(up-low)/nrebin*(i+0.5);      
      NArray[i] = FuncFromHisto(Plot,xc,smooth);
      sum += NArray[i];
    }

    Plot->Reset();
    Plot->SetBins(nrebin,low,up);
    
    Plot68->Reset();
    Plot68->SetBins(nrebin,low,up);
    Plot95->Reset();
    Plot95->SetBins(nrebin,low,up);
 
    for (int i=0;i<nrebin;i++){
      NArray[i]/=sum;
      double xc = Plot->GetBinCenter(i+1);
      Plot->Fill(xc,NArray[i]);
      Plot68->Fill(xc,NArray[i]);
      Plot95->Fill(xc,NArray[i]);
    }
  }

  //compute coordinates
  double *coord = new double[nrebin];
  for (int i=0;i<nrebin;i++){
    coord[i] = low+(up-low)/nrebin*i;
  }

  double width=(low-up)/nrebin;
  int k=0;
  double sumold,sum=0;

  // compute CLs and the corresponding coordinates
  //
  for (int i=0;i<nlevel;i++)              xcoo[i]=low;
  for (int i=nlevel+1;i<nlevel*2+1;i++)   xcoo[i]=up;

  sum = NArray[0];
  for (int i=0;i<nlevel*2+1;i++){
    if (sum>chk[i]) {
      xcoo[i]=low;
      icls[i]=0;      
      k++;
    }
  }
  for (int i=1;i<nrebin;i++){
    sumold = sum;
    sum += NArray[i];
    if ( sum>=chk[k] && sumold<chk[k]){
      double tobeadded=(chk[k]-sumold)/NArray[i];
      xcoo[k]=coord[i]+tobeadded*width;
      k++;
    }
  }
  if (ldebug) cout << "sum: " << sum << endl;
  k=0;
  sum=0;
  for (int i=0;i<nlevel;i++){
    icls[i]=0;
    icls[nlevel*2-1-i]=nrebin;
  }

  int j=0;
  
  sum = NArray[0];
  for (int i=0;i<nlevel*2+1;i++){
    if (sum>chk[i]) {
      icls[i]=0;      
      k++;
    }
  }
  
  for (int i=0;i<nrebin;i++){
    sumold = sum;
    sum+= NArray[i];
    if ( sum>chk[k] && sumold<chk[k]){
      icls[k] = i;
      k++;
    }
  }

  // round interval coordinates
  norm=pow(10.,3-floor(log10(fabs(xcoo[nlevel]))));
  if(ldebug) cout << xcoo[nlevel] << " and " << norm << endl;

  double xcoo_up  = fabs(xcoo[nlevel+1]-xcoo[nlevel]);
  double xcoo_low = fabs(xcoo[nlevel-1]-xcoo[nlevel]);

  for (int i=0;i<nlevel*2+1;i++){
    xcoo[i]=round(xcoo[i]*norm)/norm;
  }
  xcoo_up  = round(xcoo_up *norm)/norm;
  xcoo_low = round(xcoo_low*norm)/norm;

  // print the results
  cout << " Results for " << plot << endl;
  cout <<"   "<< plot << " = " << xcoo[nlevel] << " +"<< xcoo_up << " -" << xcoo_low << endl;
  for (int j=0;j<nlevel;j++){
    if (lsigma){
      cout <<"   "<< plot << " in [" << xcoo[nlevel-1-j] << "," << xcoo[nlevel+1+j] <<"]" <<" at "<< (j+1) << " sigma";
      if (j==0)
	cout <<" corresponding to "<< (xcoo[nlevel-1-j]+xcoo[nlevel+1+j])/2 << " +/- " <<  (xcoo[nlevel+1+j]-xcoo[nlevel-1-j])/2 <<endl;
      else
	cout << endl;
    } else {
      cout <<"   "<< plot << " in [" << xcoo[nlevel-1-j] << "," << xcoo[nlevel+1+j] <<"]" <<" at "<< valref[j] << "%";
      if (j==0)
	cout <<" corresponding to "<< (xcoo[nlevel-1-j]+xcoo[nlevel+1+j])/2 << " +/- " << (xcoo[nlevel+1+j]-xcoo[nlevel-1-j])/2 <<endl;
      else
	cout << endl;
    }
  }

  //PlotToBeDrawn->Set(nrebin,NArray);
  for (int i=0;i<nrebin;i++){
    if (i<icls[nlevel-2] || i>icls[nlevel+2]){
      Plot95->SetBinContent(i+1,0.0);
    }
  }
  for (int i=0;i<nrebin;i++){
    if (i<icls[nlevel-1] || i>icls[nlevel+1]){
      Plot68->SetBinContent(i+1,0.0);
    }
  }
  
  return 0;
}

int Do_Cont_Moda(const TString plot, TH1D* &Plot, double *GArray, int nbin, double *NArray, int &nrebin){

  double low = Plot->GetXaxis()->GetXmin();
  double up  = Plot->GetXaxis()->GetXmax();
  // delete all the underflow and overflow, otherwise fucked root 
  // take them into account in the normalization
  Plot->SetBinContent(0,0.);
  Plot->SetBinContent(1+Plot->GetNbinsX(),0.);

  // double norm=1;
  // for (int i=6;i>-6;i--){
  //   if (fabs(up)>=pow(10.,i)  && fabs(up)<pow(10.,i+1)) norm = pow(10.,i-3);
  // }

  double norm=pow(10.,floor(log10(up))-3);


  int newbin = (int) ceil((up-low)/norm);

  if (ldebug) cout << "up " << up << " low " << low << " norm " << norm << endl;
  if (ldebug) cout << "nrebin " << nrebin << " newbin " << newbin << endl;
  if (nrebin<newbin){
    nrebin = newbin;

    delete[] NArray;
    NArray = new double[nrebin];

    //    TF1 *fh;
    //fh = Plot->GetFunction("xx");
    double sum=0;
    for (int i=0;i<nrebin;i++){
      double xc = low+(up-low)/nrebin*(i+0.5);

//       if(i<4) {
// 	NArray[i] = Plot->GetBinContent(i);
//       } else { 
	
	NArray[i] = FuncFromHisto(Plot,xc,smooth);
//       }
      //NArray[i] = fh->Eval(xc);
      sum += NArray[i];
    }

    Plot->Reset();
    Plot->SetBins(nrebin,low,up);
        
    Plot68->Reset();
    Plot68->SetBins(nrebin,low,up);
    Plot95->Reset();
    Plot95->SetBins(nrebin,low,up);
    

    for (int i=0;i<nrebin;i++){
      NArray[i]/=sum;
      double xc = Plot->GetBinCenter(i+1);
      Plot->Fill(xc,NArray[i]);
      Plot68->Fill(xc,NArray[i]);
      Plot95->Fill(xc,NArray[i]);
    }

  }


  list <hdata> myhis;
  for (int i=0;i<nrebin;i++){
    hdata xx(i,NArray[i]);
    myhis.push_back(xx);    
  }

  myhis.sort();

  list <hdata>::iterator it;
  double bsum=0;
  double *vbin[5];
  for (int ib=0;ib<5;ib++){
    vbin[ib] = new double[nrebin];
  }
  for (it=myhis.begin();it!=myhis.end();it++){
    int ind = it->Index();
    bsum += it->Cont();
    for (int il=0;il<nlevel;il++){
      vbin[il][ind]=0;
      if (bsum<=valref[il]/100) vbin[il][ind]=1;
    }
  }
  double vsum[5]={0.0,0.0,0.0,0.0,0.0};
  for (int i=0;i<nrebin;i++){
    for (int j=0;j<nlevel;j++){
      if (vbin[j][i]==1) vsum[j]+=NArray[i];
    }
  }
  cout << "%make1Dplots: computed CLs - ";
  for (int j=0;j<nlevel;j++){
    cout << vsum[j] << " " ;
  }
  cout  << endl;

  //Compute Cls..

  double center=Plot->GetBinCenter((myhis.begin())->Index()+1);

  double xlow[5][ncoo],xup[5][ncoo];

  double dx = Plot->GetBinCenter(2)-Plot->GetBinCenter(1);
  int nint1[5],nint2[5];
  for (int j=0;j<nlevel;j++){
    nint1[j]=0; nint2[j]=0;
    if (vbin[j][0]==1)  xlow[j][nint1[j]++]=low;
    for (int i=0;i<nrebin-1;i++){
      if (vbin[j][i]==0 && vbin[j][i+1]==1) xlow[j][nint1[j]++] = Plot->GetBinCenter(i+2)-dx/2;
      if (vbin[j][i]==1 && vbin[j][i+1]==0) xup[j][nint2[j]++]  = Plot->GetBinCenter(i+1)+dx/2;
    }
    if (vbin[j][nrebin-1]==1)  xup[j][nint2[j]++]=up;
  }

  
  cout << " Results for " << plot << endl;
  double sigp=0,sigm=0;
  if (nlevel!=0){
    for (int i=0;i<nint1[0];i++){
      if (center>xlow[0][i] && center<xup[0][i]){
	sigp = xup[0][i]-center;
	sigm = center-xlow[0][i];
      }
    }
  }
  cout <<"   "<< plot << " = " << center << " +"<< sigp << " -" << sigm << endl;
  for (int j=0;j<nlevel;j++){
    if (lsigma) {
      cout << "     at " << (j+1) <<" sigma ";
    } else {
      cout << "     at " << valref[j]<<"% ";
    }
    for (int i=0;i<nint1[j];i++){
      cout << "["<< xlow[j][i] << ",";
      cout << xup[j][i] <<"] ";
      if (j==0)
	cout <<" corresponding to "<< (xlow[j][i]+xup[j][i])/2 << " +/- " <<  (xup[j][i]-xlow[j][i])/2;
    }
    cout << endl;
  }

  //    PlotToBeDrawn->Set(nrebin,NArray);
  for (int i=0;i<nrebin;i++){
    if (vbin[1][i]==0){
      Plot95->SetBinContent(i+1,0.0);
    }
    if (vbin[0][i]==0){
      Plot68->SetBinContent(i+1,0.0);
    }
  }

  return 0;
}

void DefineColors() {
  TColor *color[120];
  color[0]  = new TColor(1300, 1.0, 1.0, 1.0, "");
  color[1]  = new TColor(1301, 0.0, 0.0, 0.0, "");
  color[2]  = new TColor(1302, 1.0, 0.0, 0.0, "");
  color[3]  = new TColor(1303, 0.0, 1.0, 0.0, "");
  color[4]  = new TColor(1304, 0.0, 0.0, 1.0, "");
  color[11] = new TColor(1311, 0.6, 0.4, 1.0, "");
  color[12] = new TColor(1312, 1.0, 0.5, 0.5, "");
  color[13] = new TColor(1313, 0.5, 1.0, 0.4, "");
  color[14] = new TColor(1314, 0.8, 0.6, 1.0, "");
  color[15] = new TColor(1315, 1.0, 0.90, 0.82, ""); 
  color[16] = new TColor(1316, 1.0, 0.75, 0.67, "");
  color[17] = new TColor(1317, 0.5, 1.0, 0.5, "");
  color[18] = new TColor(1318, 0.95, 1.0, 1.0, "");
  color[19] = new TColor(1319, 0.7, 1.0, 0.8, "");
  color[20] = new TColor(1320, 0.97, 0.87, 0.97, "");
  color[21] = new TColor(1321, 0.85, 0.44, 0.84, "");
  color[22] = new TColor(1322, 0.73, 0.33, 0.83, "");
  color[23] = new TColor(1323, 0.60, 0.20, 0.80, "");
  color[24] = new TColor(1324, 0.58, 0.00, 0.83, "");
  color[25] = new TColor(1325, 0.54, 0.17, 0.89, "");
  color[26] = new TColor(1326, 0.63, 0.13, 0.94, "");
  color[27] = new TColor(1327, 0.58, 0.44, 0.86, "");
  color[28] = new TColor(1328, 0.95, 0.85, 0.95, "");
  color[29] = new TColor(1329, 0.93, 0.51, 0.93, "");
  color[30] = new TColor(1330, 1.00, 0.98, 0.98, "");      
  color[31] = new TColor(1331, 0.97, 0.97, 1.00, "");
  color[32] = new TColor(1332, 1.00, 1.00, 1.00, ""); 
  color[33] = new TColor(1333, 0.52, 0.44, 1.00, ""); 
  color[34] = new TColor(1334, 0.12, 0.56, 1.00, "");
  color[35] = new TColor(1335, 0.00, 0.75, 1.00, "");
  color[36] = new TColor(1336, 0.88, 1.00, 1.00, "");
  color[37] = new TColor(1337, 0.37, 0.62, 0.63, "");
  color[38] = new TColor(1338, 0.40, 0.80, 0.67, "");
  color[39] = new TColor(1339, 0.50, 1.00, 0.83, "");
  color[40] = new TColor(1340, 0.00, 0.39, 0.00, "");
  color[41] = new TColor(1341, 0.33, 0.42, 0.18, "");
  color[42] = new TColor(1342, 0.56, 0.74, 0.56, "");
  color[43] = new TColor(1343, 0.18, 0.55, 0.34, "");
  color[44] = new TColor(1344, 0.24, 0.70, 0.44, "");
  color[45] = new TColor(1345, 0.13, 0.70, 0.67, "");
  color[46] = new TColor(1346, 0.60, 0.98, 0.60, "");
  color[47] = new TColor(1347, 0.00, 1.00, 0.50, "");
  color[48] = new TColor(1348, 0.49, 0.99, 0.00, "");
  color[49] = new TColor(1349, 0.00, 1.00, 0.00, "");
  color[50] = new TColor(1350, 0.50, 1.00, 0.00, "");
  color[51] = new TColor(1351, 0.00, 0.98, 0.60, "");
  color[52] = new TColor(1352, 0.68, 1.00, 0.18, "");
  color[53] = new TColor(1353, 0.20, 0.80, 0.20, "");
  color[54] = new TColor(1354, 0.60, 0.80, 0.20, "");
  color[55] = new TColor(1355, 0.13, 0.55, 0.13, "");
  color[56] = new TColor(1356, 0.42, 0.56, 0.14, "");
  color[57] = new TColor(1357, 0.74, 0.72, 0.42, "");
  color[58] = new TColor(1358, 0.94, 0.90, 0.55, "");
  color[59] = new TColor(1359, 0.93, 0.91, 0.67, "");
  color[60] = new TColor(1360, 0.99, 0.99, 0.75, "");
  color[61] = new TColor(1361, 1.00, 0.70, 0.00, "");
  color[62] = new TColor(1362, 1.00, 1.00, 0.00, "");
  color[63] = new TColor(1363, 1.00, 0.84, 0.00, "");
  color[64] = new TColor(1364, 1.00, 0.95, 0.60, "");
  color[65] = new TColor(1365, 0.85, 0.65, 0.13, "");
  color[66] = new TColor(1366, 0.72, 0.53, 0.04, "");
  color[67] = new TColor(1367, 0.74, 0.56, 0.56, "");
  color[68] = new TColor(1368, 0.80, 0.36, 0.36, "");
  color[69] = new TColor(1369, 0.55, 0.27, 0.07, "");
  color[70] = new TColor(1370, 0.63, 0.32, 0.18, "");
  color[71] = new TColor(1371, 0.80, 0.52, 0.25, "");
  color[72] = new TColor(1372, 0.99, 0.85, 0.70, "");
  color[73] = new TColor(1373, 0.96, 0.96, 0.86, "");
  color[74] = new TColor(1374, 0.99, 0.90, 0.75, "");
  color[75] = new TColor(1375, 0.96, 0.64, 0.38, "");
  color[76] = new TColor(1376, 0.82, 0.71, 0.55, "");
  color[77] = new TColor(1377, 0.82, 0.41, 0.12, "");
  color[78] = new TColor(1378, 0.70, 0.13, 0.13, "");
  color[79] = new TColor(1379, 0.65, 0.16, 0.16, "");
  color[80] = new TColor(1380, 0.91, 0.59, 0.48, "");
  color[81] = new TColor(1381, 0.98, 0.50, 0.45, "");
  color[82] = new TColor(1382, 1.00, 0.85, 0.80, "");
  color[83] = new TColor(1383, 1.00, 0.65, 0.00, "");
  color[84] = new TColor(1384, 1.00, 0.55, 0.00, "");
  color[85] = new TColor(1385, 1.00, 0.50, 0.31, "");
  color[86] = new TColor(1386, 0.94, 0.50, 0.50, "");
  color[87] = new TColor(1387, 1.00, 0.39, 0.28, "");
  color[88] = new TColor(1388, 1.00, 0.27, 0.00, "");
  color[89] = new TColor(1389, 1.00, 0.00, 0.00, "");

  color[90] = new TColor(1390, 0.90, 0.60, 0.60, ""); //red
  color[91] = new TColor(1391, 0.70, 0.25, 0.25, "");
  color[92] = new TColor(1392, 0.87, 0.87, 0.91, ""); //blue
  color[93] = new TColor(1393, 0.59, 0.58, 0.91, ""); 
  color[94] = new TColor(1394, 0.65, 0.55, 0.85, ""); //violet (gamma)
  color[95] = new TColor(1395, 0.49, 0.26, 0.64, ""); 
  color[96] = new TColor(1396, 0.95, 0.95, 0.45, ""); // yellow (alpha)
  color[97] = new TColor(1397, 0.95, 0.95, 0.05, "");
  color[98] = new TColor(1398, 0.75, 0.92, 0.68, ""); //green (2beta+gamma)
  color[99] = new TColor(1399, 0.36, 0.57, 0.30, "");
  color[100] = new TColor(1400, 0.97, 0.50, 0.09, ""); // orange
  color[101] = new TColor(1401, 0.76, 0.34, 0.09, "");
  color[102] = new TColor(1402, 0.97, 0.52, 0.75, ""); // pink 
  color[103] = new TColor(1403, 0.76, 0.32, 0.51, "");
  color[104] = new TColor(1404, 0.49, 0.60, 0.82, ""); // dark blue (kpnn)
  color[105] = new TColor(1405, 0.43, 0.48, 0.52, "");
  color[106] = new TColor(1406, 0.70, 0.70, 0.70, "");  // black
  color[107] = new TColor(1407, 0.40, 0.40, 0.40, "");
  color[108] = new TColor(1408, 0.10, 0.60, 0.10, "");  // dark green
  color[109] = new TColor(1409, 0.10, 0.40, 0.10, "");
}

int make1Dplots(TString sfile, TString plot, TString fit="none"){


  TH1D*  Plot=NULL;

  int nrebin;

  sfile+=".root";
  datafile = new TFile(sfile);
  if (!datafile){
    cout << "%make1Dplots: error opening file" << endl;
  }
  TDirectory *dir=NULL;
  dir = (TDirectory *)datafile->Get(fit);
  if (dir==NULL){
    if (fit!="none") {
      cout << "%make1Dplots: fit " << fit << " does not exist " << endl;
      return -1;
    } else {
      cout << "%make1Dplots: fit " << fit << " does not exist: working in the main directory... " << endl;
    }
  }

  TH1D* his;
  int ind=-1;
  if (fit != "none"){
    his = (TH1D*)datafile->Get(fit+"/"+plot);
  } else {
    ind = plot.Index("=");
    if (ind!=-1) {
      TString tmp(plot(ind+1,plot.Length()-ind+1));
      titletmp = tmp;
      TString tmpplot(plot(0,ind));
      his = (TH1D*)datafile->Get(tmpplot);
      plot = plot(ind+1,plot.Length());
    } else {
      his = (TH1D*)datafile->Get(plot);
    }
  }

  if (his==NULL){
    cout << "%make1Dplots: plot " << plot << " does not exist..." << endl;
    return -1;
  }

  Plot = (TH1D*)his->Clone();
  PlotToBeDrawn = new TH1D(*Plot);
  
  int nbin   = Plot->GetNbinsX();
  double low = Plot->GetXaxis()->GetXmin();
  double up  = Plot->GetXaxis()->GetXmax();
  
  double *GArray = new double[nbin];
  double *NArray = new double[nbin];

  for (int i=0;i<nbin;i++) {
    GArray[i]=0.0;
    if (Plot->GetBinContent(i+1)>0) GArray[i]=Plot->GetBinContent(i+1);
    NArray[i]=GArray[i];
  }

  // defines the histogram to be plotted (rebinned w.r.t.
  // original one)

  int nSignChange=100;
  int i=0;

  nrebin = nbin;
  if (rebin==0){
    while (nSignChange>5 && i<nvalrebin){
      PlotToBeDrawn = (TH1D*)Plot->Rebin(val[i],"hnew");
      nrebin        = PlotToBeDrawn->GetNbinsX();
      nSignChange   = EvaluateSignChanges(PlotToBeDrawn);
      for (int j=0;j<nrebin;j++) {
	NArray[j]=PlotToBeDrawn->GetBinContent(j+1);
      }
      i++;
    }
  } else {
    PlotToBeDrawn = (TH1D*)Plot->Rebin(rebin,"hnew");    
    nrebin        = PlotToBeDrawn->GetNbinsX();
    for (int j=0;j<nrebin;j++) {
      NArray[j]=PlotToBeDrawn->GetBinContent(j+1);
    }
  }
  
  if (ldebug) cout << "Number of Sign Changes " << nSignChange << endl;

  PlotToBeDrawn = new TH1D("plot","",nrebin,low,up);
  Plot68        = new TH1D("plot68","",nrebin,low,up);
  Plot95        = new TH1D("plot95","",nrebin,low,up);

  // compute overall sum
  double xsum=0;  
  for (int i=0;i<nbin;i++){
    xsum += GArray[i];
  }

  for (int i=0;i<nbin;i++){
    GArray[i]/=xsum;
  }
  for (int i=0;i<nrebin;i++){
    NArray[i]/=xsum;
  }
  
  for (int i=0;i<nrebin;i++){
    PlotToBeDrawn->SetBinContent(i+1,NArray[i]);
    Plot68->SetBinContent(i+1,NArray[i]);
    Plot95->SetBinContent(i+1,NArray[i]);
  }
  gStyle->SetOptStat(0);

  DefineColors();
  // SetDefaultPalette();
  // gStyle->SetTitleSize(0.06);
  tela = new TCanvas("his1D","1D histograms",3);
  tela->SetLeftMargin(0.20);
  tela->SetBottomMargin(0.15);
  tela->SetFillColor(0);
  tela->SetFrameFillColor(0);
  tela->SetHighLightColor(0);
  tela->SetBorderMode(0);

  // PlotToBeDrawn->SetTitle("");
  // PlotToBeDrawn->GetXaxis()->SetTitleSize(0.06);
  // PlotToBeDrawn->GetYaxis()->SetTitleSize(0.06);
  // PlotToBeDrawn->GetXaxis()->SetTitleFont(50);
  // PlotToBeDrawn->GetYaxis()->SetTitleFont(50);
  // PlotToBeDrawn->GetXaxis()->SetNdivisions(505);
  // PlotToBeDrawn->GetYaxis()->SetNdivisions(505);
  // PlotToBeDrawn->GetYaxis()->SetTitleOffset(1.5);

  //PlotToBeDrawn->SetMinimum(0.);

  if (!lylab){
    PlotToBeDrawn->GetYaxis()->SetTitle("Probability density");
  } else {
    PlotToBeDrawn->GetYaxis()->SetTitle(ylab);
  }
  SetXTitle(PlotToBeDrawn,plot);

  PlotToBeDrawn->Draw("HIST");

  if (lmediana){
    Do_Cont_Mediana(plot, PlotToBeDrawn, GArray, nbin, NArray, nrebin);
  } else {
    Do_Cont_Moda(plot, PlotToBeDrawn, GArray, nbin, NArray, nrebin);
  }

  PlotToBeDrawn->SetTitle("");
  PlotToBeDrawn->GetXaxis()->SetTitleSize(0.06);
  PlotToBeDrawn->GetYaxis()->SetTitleSize(0.06);
  //  PlotToBeDrawn->GetXaxis()->SetTitleFont(50);
  //  PlotToBeDrawn->GetYaxis()->SetTitleFont(50);
  PlotToBeDrawn->GetXaxis()->SetNdivisions(505);
  PlotToBeDrawn->GetYaxis()->SetNdivisions(505);
  PlotToBeDrawn->GetYaxis()->SetTitleOffset(1.5);

  //Plot95->SetFillColor(kYellow);
  // Plot95->SetFillColor(398);
  Plot95->SetFillColor(col2);
  Plot95->SetLineStyle(3);
  Plot95->Scale(1/Plot95->GetBinWidth(1));
  Plot95->Draw("HIST same");

  Plot68->SetLineStyle(3);
  //Plot68->SetFillColor(kRed);
  // Plot68->SetFillColor(399);
  Plot68->SetFillColor(col1);
  Plot68->Scale(1/Plot68->GetBinWidth(1));
  Plot68->Draw("HIST same");
  PlotToBeDrawn->SetLineStyle(1);
  PlotToBeDrawn->Scale(1/PlotToBeDrawn->GetBinWidth(1));
  PlotToBeDrawn->Draw("HIST same");
  PlotToBeDrawn->Draw("HIST sameaxis");
  
  TH1D *hsup1;
  TH1D *hsup2;
  TH1D *hsup3;

  if (lsup1){
    double norm = PlotToBeDrawn->GetSum()*PlotToBeDrawn->GetBinWidth(1);
    hsup1 = (TH1D*)datafile->Get(dir1+"/"+fit1);
    double norm1 = hsup1->GetSum()*hsup1->GetBinWidth(1);
    hsup1->Scale(norm/norm1);
    hsup1->SetLineWidth(3);
    hsup1->SetLineColor(colsup);
    hsup1->Draw("HIST SAME");
  }

  if (lsup2){
    double norm = PlotToBeDrawn->GetSum()*PlotToBeDrawn->GetBinWidth(1);
    hsup2 = (TH1D*)datafile->Get(dir2+"/"+fit2);
    double norm2 = hsup2->GetSum()*hsup2->GetBinWidth(1);
    hsup2->Scale(norm/norm2);
    hsup2->SetLineWidth(3);
    hsup2->SetLineColor(colsup);
    hsup2->SetLineStyle(2);
    hsup2->Draw("HIST SAME");
  }

  if (lsup3){
    double norm = PlotToBeDrawn->GetSum()*PlotToBeDrawn->GetBinWidth(1);
    hsup3 = (TH1D*)datafile->Get(dir3+"/"+fit3);
    double norm3 = hsup3->GetSum()*hsup3->GetBinWidth(1);
    hsup3->Scale(norm/norm3);
    hsup3->SetLineWidth(3);
    hsup3->SetLineColor(colsup);
    hsup3->SetLineStyle(3);
    hsup3->Draw("HIST SAME");
  }

  if(llogo)
    PutLogo();
  gPad->Update();
  tela->cd();
  tela->RedrawAxis("same");

  //  delete[] NArray;
  //  delete[] GArray;

  return 0;
}

int main(int argc, char *argv[]){

  if (argc<4 && argc!=3){
    cout << "####################################################"   << endl;
    cout << "%make1Dplots:                                       "   << endl;
    cout << " - needed parameters                                "   << endl;
    cout << "  'filename' (w/o extension - .root assumed)        "   << endl;
    cout << "  'plotname' (e.g. sin2b   )                        "   << endl;
    cout << "     for files converted from paw use h232=\"sin2b\""  << endl;              
    cout << "  'fitname-directory'  (=none no dir struct.)"   << endl;
    cout << endl;
    cout << " e.g. " << endl;
    cout << "  > make1Dplots output sin2b fit1           "   << endl;
    cout << "    read from file output.root and plot the histogram " << endl;
    cout << "    for sin2b resident in the fit directory fit1      " << endl;
    cout << endl;
    cout << " - to check the containts of the filename type " << endl;
    cout << "  > make1Dplots filename (dir) info"    << endl;
    cout << endl;
    cout << " - optional parameters                           " << endl;
    cout << "   --nograf  -> no graphical output " << endl;
    cout << "   --makeps  -> print .eps file     " << endl;
    cout << "   --sigma   -> express CLs in term of sigma (0.68,0.954,0.997)" << endl;
    cout << "                (default is 0.68,0.95,0.99,0.999, etc..)   " << endl;
    cout << "   --median  -> compute CLs w.r.t. median (default from most prob. value)" << endl;
    cout << "   --nologo      -> don't plot UTfit logo" << endl;
    cout << "   -season=winter10  -> attach the date information" << endl;
    cout << "   -addtext=test  -> attach additional information" << endl;
    cout << "   -nlevel=n -> specify the number of CLs (default n=3)   " << endl;
    cout << "   -xlab=namex  -> x label (default read from histogram)" << endl;
    cout << "   -ylab=namey  -> y label (default 'Probability density')" << endl;
    cout << "   -col1=index  -> color index for the 1 sigma interval" << endl;
    cout << "   -col2=index  -> color index for the 2 sigma interval" << endl;
    cout << "   -rebin=ndiv  -> rebin histogram to nbins/ndiv    " << endl;
    cout << "   -smooth=ntime  -> smooth with ROOT (ntime = number of iterative smoothing)" << endl;
    cout << "   -fit1=fitname  -> name of the first fit to be superimposed" << endl;
    cout << "   -fit2=fitname  -> name of the second fit to be superimposed" << endl;
    cout << "   -dir1=dirname  -> name of the dir containing fit1 " << endl;
    cout << "   -dir2=dirname  -> name of the dir containing fit2 " << endl;
    cout << "   -colsup=solsup -> color to be used for the superimposed plots" << endl;
    cout << "   -name=plotname -> output plot name, with extension" << endl;
    cout << endl;
    cout << " How to use it in a interactive root session:  " << endl;
    cout << "  root [0] .L make1Dplots.cc                    " << endl;
    cout << "  root [1] make1Dplots(filename,plotname,fitname) "<<endl;
    cout << "####################################################"   << endl;
    return 0;
  }

  char str[100];
  for (int i=1;i<argc;i++){

    if (strncmp(argv[i],"-fit1",5)==0){
      sscanf(argv[i],"-fit1=%s",str);
      TString fit(str,strlen(str));
      fit1  = fit;
      lsup1 = true;
    }
    if (strncmp(argv[i],"-dir1",5)==0){
      sscanf(argv[i],"-dir1=%s",str);
      TString dir(str,strlen(str));
      dir1 = dir;
    }

    if (strncmp(argv[i],"-fit2",5)==0){
      sscanf(argv[i],"-fit2=%s",str);
      TString fit(str,strlen(str));
      fit2  = fit;
      lsup2 = true;
    }
    if (strncmp(argv[i],"-dir2",5)==0){
      sscanf(argv[i],"-dir2=%s",str);
      TString dir(str,strlen(str));
      dir2 = dir;
    }

    if (strncmp(argv[i],"-fit3",5)==0){
      sscanf(argv[i],"-fit3=%s",str);
      TString fit(str,strlen(str));
      fit3  = fit;
      lsup3 = true;
    }
    if (strncmp(argv[i],"-dir3",5)==0){
      sscanf(argv[i],"-dir3=%s",str);
      TString dir(str,strlen(str));
      dir3 = dir;
    }

    if (strncmp(argv[i],"-season",7)==0){
      sscanf(argv[i],"-season=%s",str);
      TString label(str,strlen(str));
      season = label;
    }

    if (strncmp(argv[i],"-addtext",8)==0){
      sscanf(argv[i],"-addtext=%s",str);
      TString label(str,strlen(str));
      addtxt = label;
    }

    if (strncmp(argv[i],"--nologo",8)==0)  llogo = false;
    if (strncmp(argv[i],"--median",8)==0)  lmediana = true;
    if (strncmp(argv[i],"--nograf",8)==0)  lgraphic = false;
    if (strncmp(argv[i],"--makeps",8)==0)  lmakeps  = true;
    if (strncmp(argv[i],"--debug",7)==0)  ldebug   = true;
    if (strncmp(argv[i],"--precise",9)==0) lprecise = true;
    if (strncmp(argv[i],"--sigma",7)==0) lsigma = true;
    if (strncmp(argv[i],"-nlevel",7)==0) sscanf(argv[i],"-nlevel=%d",&nlevel);
    if (strncmp(argv[i],"-name",5)==0) {
      sscanf(argv[i],"-name=%s",str);
      lname= TString(str,strlen(str));
    }
    if (strncmp(argv[i],"-xlab",5)==0){
      sscanf(argv[i],"-xlab=%s",str);
      //TString sx(str,strlen(str));
      //xlab  = sx;

      char tmp[40];
      char sx1[40];
      sprintf(tmp,"%s",str);
      for (int c=1; c<=40; c++) {
	if (tmp[c]=='_')
	  if (tmp[c+1]!='{')
	    tmp[c]=' ';
      }
      sprintf(sx1,"%s",tmp);
      xlab  = sx1;

      lxlab = true;
    }
    if (strncmp(argv[i],"-ylab",5)==0){ 
      sscanf(argv[i],"-ylab=%s",str);
      //TString sy(str,strlen(str));
      //ylab  = sy;

      char tmp[40];
      char sy1[40];
      sprintf(tmp,"%s",str);
      for (int c=1; c<=40; c++) {
	if (tmp[c]=='_') 
	  if (tmp[c+1]!='{')
	    tmp[c]=' ';
      }
      sprintf(sy1,"%s",tmp);
      ylab  = sy1;
      
      lylab = true;
    }
    if (strncmp(argv[i],"-col1",5)==0)  sscanf(argv[i],"-col1=%d",&col1);
    if (strncmp(argv[i],"-col2",5)==0)  sscanf(argv[i],"-col2=%d",&col2);
    if (strncmp(argv[i],"-rebin",6)==0) sscanf(argv[i],"-rebin=%d",&rebin);
    if (strncmp(argv[i],"-smooth",7)==0) sscanf(argv[i],"-smooth=%d",&smooth);

    if (strncmp(argv[i],"-colsup",5)==0)  sscanf(argv[i],"-colsup=%d",&colsup);

  }
  //  if (lmakeps) lgraphic = true;

  theapp = new TApplication("app",0,NULL);


  TString sfile((const char *)argv[1],strlen(argv[1])); 
  if (argc==3 && strncmp(argv[2],"info",4)==0){
    sfile+=".root";
    datafile = new TFile(sfile);
    datafile->ls();
    return 0;
  }
  if (argc==4 && strncmp(argv[3],"info",4)==0){
    sfile+=".root";
    datafile = new TFile(sfile);
    datafile->cd(argv[2]);
    datafile->ls();
    return 0;
  }

  TString plot((const char *)argv[2],strlen(argv[2])); 
  TString fit((const char *)argv[3],strlen(argv[3])); 

  if (dir1=="") dir1 = fit;
  if (dir2=="") dir2 = fit;
  if (dir3=="") dir3 = fit;

  int ret=make1Dplots(sfile,plot,fit);

  if (lgraphic) {
    cout <<"%make1Dplots: select Quit ROOT from the menu file to continue " << endl;
    theapp->Run(true);
  }
  if (lmakeps){
    int ind=-1;
    for(int i=0;i<plot.Length();i++){
      ind = plot.Index("=");
    }
    if (ind!=-1) {
      plot = plot(ind+1,plot.Length());
    }
    if(lname!="none"){
      tela->SaveAs("./plots/"+lname);
    } else {
      TString psfile = "./plots/"+sfile+"_"+fit+"_"+plot+".eps";
      tela->SaveAs(psfile);
    }
  }
  if (datafile!=NULL) datafile->Close();
  tela->Close();

  return ret;

}

