#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <TH1.h>
#include <TF2.h>
#include <TH2.h>
#include <TROOT.h>
#include <TFile.h>
#include <TImage.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TColor.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TLine.h>

using namespace std;

TString season,addtxt="";
TString xlab,ylab;
bool lxlab=false,lylab=false;

TApplication   *theapp=NULL;
TFile        *datafile=NULL;
TCanvas          *tela=NULL;
TH1D             *Plot=NULL;
TH2D           *Plot2D=NULL;
TString     titletmp="none";
//TMath                     m;
TF1                     *f1;

bool lmakeps=false;
bool llogo=true;

int nx=0,ny=0;
int    xnbin;
double x_low=0,x_up=0,y_low=0,y_up=0,hlow,hup;

int logoposition = 3; 
double finverfc(double);
double fdelta(double x);
double quadra(double x,double x1,double x2,double x3,double y1,double y2,double y3);
double integrale(double f(double), double xmin, double xmax);
double fzero(double x, double area);
double f2(double x, double y);
double fconv(double x);
void SetStyle();
void PutLogo();

double xval=-999,xerr;
double xval2=-999,xerr2;
double xval3=-999,xerr3;
double mean,sigma,delta;

double quadra(double x,double x1,double x2,double x3,double y1,double y2,double y3){
  double a = ((y3-y1)/(x3-x1)*(x2-x1)-(y2-y1))/((x3*x3-x1*x1)/(x3-x1)*(x2-x1)-(x2*x2-x1*x1));
  double b = (y2-y1-a*(x2*x2-x1*x1))/(x2-x1);
  double c = y1-b*x1-a*x1*x1;
  return a*x*x+b*x+c;
}

double fhisto(double x){

  double val;

  int ibin = Plot->FindBin(x);
  if (ibin<0 || ibin>=xnbin) 
    val = 0.;
  else
    val = Plot->GetBinContent(ibin); 

  return val;

}

double fdelta(double x){
  delta = x;
  return integrale(fconv,hlow,hup);  
}

double integrale(double f(double),double xmin, double xmax){
  
  int n=150;
  double val=0;
  double x0,x1,x2,w=(xmax-xmin)/n;
  for (int i=0;i<n;i+=2){
    x0 = xmin+i*w;
    x1 = x0+w;
    x2 = x1+w;
    val+= w/3*(f(x0)+4*f(x1)+f(x2));
  }
  return val;
}

double fconv(double x){
  double f1,f2;
  f1 = fhisto(x);
  f2 = TMath::Gaus(x-delta,mean,sigma,true);
  if (f1<=0) f1=0;
  if (f2<=0) f2=0;
  return f1*f2;
}


double f2(double x, double y){

  double sump,val,tot;

  mean  = x;
  sigma = y;

  //tot = 1.; sump = 0.2;
  tot  = integrale(fdelta,-(x_up-x_low),0.0);
  sump = integrale(fdelta,0.0,x_up-x_low);
  sump/=(sump+tot);

  if (sump>0.5) sump = 1-sump;
  val = finverfc(sump);

  if (val>=5.9) val = 6.01;
  return val;

}


void SetStyle(){
  tela->SetLeftMargin(0.20);
  tela->SetRightMargin(0.20);
  tela->SetBottomMargin(0.15);
  tela->SetFillColor(0);
  tela->SetBorderMode(0);
  Plot2D->GetXaxis()->SetTitleSize(0.06);
  Plot2D->GetYaxis()->SetTitleSize(0.06);
  Plot2D->GetXaxis()->SetNdivisions(506);
  Plot2D->GetYaxis()->SetNdivisions(510);
  Plot2D->GetXaxis()->SetTitleOffset(0.85);
  Plot2D->GetYaxis()->SetTitleOffset(1.40);
  //Plot2D->GetYaxis()->SetTitleOffset(1.00);
}

/*
void PutLogo(){
  TImage *img = TImage::Open("../common/logo.png");
  // we should add an option to put the logo in the corner we want
  TPad *pad;
  if(logoposition==1) 
    //bottom-left
    pad = new TPad("pad","pad",0.22,0.16,0.37,0.27,0,0,0);
  else if(logoposition==2) 
    //top-left
    pad = new TPad("pad","pad",0.22,0.77,0.37,0.88,0,0,0);
  else if(logoposition==3) 
    //top-right
    //pad = new TPad("pad","pad",0.64,0.77,0.79,0.88,0,0,0);
    pad = new TPad("pad","pad",0.61,0.75,0.78,0.88,0,0,0);
  else if(logoposition==4) 
    // bottom-right
    pad = new TPad("pad","pad",0.64,0.16,0.79,0.27,0,0,0);

  //img->SetConstRatio(kFALSE);
  pad->Draw();
  pad->cd();
  img->Draw();

}
*/

void PutLogo(){
  TImage *img = TImage::Open("../common/logo.png");
  // we should add an option to put the logo in the corner we want
  TPad *pad;
  double img_x0=0,img_x1=0,img_y0=0,img_y1=0;
  if(logoposition==1){
    //bottom-left
    img_x0 = 0.22;
    img_x1 = 0.37;
    img_y0 = 0.16;
    img_y1 = 0.27;
  }else if(logoposition==2){
    //top-left
    img_x0 = 0.22;
    img_x1 = 0.37;
    img_y0 = 0.77;
    img_y1 = 0.88;
  }else if(logoposition==3){
    //top-right
    img_x0 = 0.64;
    img_x1 = 0.79;
    img_y0 = 0.77;
    img_y1 = 0.88;
  }else if(logoposition==4){
    // bottom-right
    img_x0 = 0.64;
    img_x1 = 0.79;
    img_y0 = 0.16;
    img_y1 = 0.27;
  }
  pad = new TPad("pad","pad",img_x0,img_y0,img_x1,img_y1,0,0,0);

  img->SetConstRatio(kFALSE);
  pad->Draw();
  pad->cd();
  img->Draw();
  gPad->Modified();
  gPad->Update();
  tela->cd();
  double x0 = img_x0;

  double x1 = img_x1;
  double y0 = img_y0-0.03;
  double y1 = y0+0.05;
  
  TPad *padSeason = new TPad("padSeason","padSeason",x0,y0,x1,y1,0,0,0);
  padSeason->SetFillStyle(1001);
  padSeason->SetFillColor(0);
  padSeason->SetBorderSize(0);
  padSeason->Draw();
  padSeason->cd();
  TLatex* tex = new TLatex;
  TColor *col = gROOT->GetColor(1);
  col->SetRGB(0.21569,0.22745,0.47843);

  tex->SetTextSize(0.6);
  tex->SetTextAlign(21);
  tex->SetTextColor(1);
  tex->DrawTextNDC(0.48,0.43, season);
  padSeason->Update();
  
  TLatex* addtex = new TLatex;
  addtex->SetTextSize(0.5);
  addtex->SetTextAlign(21);
  addtex->DrawTextNDC(0.48,0.0,addtxt);
  
}





void SetXTitle(TH2D *id, TString tit){

  TString str;

  if (titletmp!="none") tit = titletmp;

  if      (strncmp(tit,"sin2b",5)==0)      str = "sin2#beta";
  else if (strncmp(tit,"sin2a",5)==0)      str = "sin2#alpha";
  else if (strncmp(tit,"alpha",5)==0)      str = "#alpha[^{o}]";
  else if (strncmp(tit,"beta",4)==0)       str = "#beta[^{o}]";
  else if (strncmp(tit,"gamma",5)==0)      str = "#gamma[^{o}]";
  else if (strncmp(tit,"rho",3)==0)        str = "#bar{#rho}";
  else if (strncmp(tit,"eta",3)==0)        str = "#bar{#eta}";
  else if (strncmp(tit,"Lambda",6)==0)     str = "#lambda";
  else if (strncmp(tit,"Imlamt",6)==0)     str = "Im#lambda_{t}";
  else if (strncmp(tit,"Rb",2)==0)         str = "R_{b}";
  else if (strncmp(tit,"Rt",2)==0)         str = "R_{t}";
  else if (strncmp(tit,"Vub",3)==0)        str = "|V_{ub}|";
  else if (strncmp(tit,"Vcb",3)==0)        str = "|V_{cb}|";
  else if (strncmp(tit,"VtdOverVts",10)==0)   str = "V_{td}/V_{ts}";
  else if (strncmp(tit,"Vtd",3)==0)        str = "V_{td}#times10^{3}";
  else if (strncmp(tit,"DeltaMs",7)==0)    str = "#Delta m_{s}[ps^{-1}]";
  else if (strncmp(tit,"FbSqrtb",7)==0)    str = "F_{B}#sqrt{B}";
  else if (strncmp(tit,"Bk",2)==0)         str = "B_{k}";
  else if (strncmp(tit,"Xi",2)==0)         str = "#xi";
  else if (strncmp(tit,"Ftt",3)==0)        str = "F_{tt}";
  else if (strncmp(tit,"sin2bpg",7)==0)    str = "sin(2#beta+#gamma)";
  else if (strncmp(tit,"2bpg",4)==0)       str = "2#beta+#gamma[^{o}]";
  else if (strncmp(tit,"Epsk",4)==0)       str = "#epsilon_{K}";
  else if (strncmp(tit,"btaunu",6)==0)     str = "BR(B#rightarrow#tau#nu)[10^{-4}]";
  else if (strncmp(tit,"BRbarbsll",9)==0)  str = "#bar{BR}(B_{s}#rightarrowll)[10^{-9}]";
  else if (strncmp(tit,"BRbdll",6)==0)     str = "BR(B_{d}#rightarrowll)[10^{-9}]";
  else if (strncmp(tit,"BdOverBsll",10)==0)  str = "BR(B_{d}#rightarrowll)/BR(B_{s}#rightarrowll)";
  else if (strncmp(tit,"ASL_d",5)==0)      str = "A_{SL}^{d}";
  else if (strncmp(tit,"ASL_s",5)==0)      str = "A_{SL}^{s}";
  else if (strncmp(tit,"A_HC",4)==0)       str = "A_{#mu#mu}";
  else str = tit;
  
  TString laby;
  if (!lxlab){
    // b tau nu
    //TString labx2 = str+"[10^{-4}]";
    //id->GetXaxis()->SetTitle(labx2);
    //laby = "#sigma("+str+")[10^{-4}]";
    // bs to mumu
    //TString labx2 = str+"[10^{-9}]";
    //id->GetXaxis()->SetTitle(labx2);
    //laby = "#sigma("+str+")[10^{-9}]";
    // default
    id->GetXaxis()->SetTitle(str);
    laby = "#sigma("+str+")";
  } else {
    id->GetXaxis()->SetTitle(xlab);
    laby = "#sigma("+xlab+")";
  }
  id->GetYaxis()->SetTitle(laby);

}


int makepull(TString sfile, TString plot, TString fit){


  sfile+=".root";
  datafile = new TFile(sfile);

  TDirectory *dir=NULL;
  dir = (TDirectory *)datafile->Get(fit);
  if (dir==NULL){
    if (fit!="none") {
      cout << "%make1Dplots: fit " << fit << " does not exist " << endl;
      return -1;
    } else {
      cout << "%make1Dplots: fit " << fit << " does not exist " << endl;
      cout << " ... reading histograms from the main directory ... " << endl;
    }
  }

  TH1D* his;
  if (fit != "none"){
    his = (TH1D*)datafile->Get(fit+"/"+plot);  
  } else {
    int ind=-1;
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
    cout << "%makepull: plot " << plot << " does not exist " << endl;
    return -1;
  }

  gStyle->SetOptStat(0);
  tela = new TCanvas("his1D","1D histograms",10,10,700,700);

  Plot = new TH1D;
  Plot = (TH1D*)his->Clone();

  xnbin = Plot->GetNbinsX();
  double *v  = new double[xnbin];
  double  vtot = 0;
  for (int i=0;i<xnbin;i++){
    v[i] = Plot->GetBinContent(i+1);
    vtot += v[i];
  }
  for (int i=0;i<xnbin;i++){
    v[i]/= vtot;
  }
  Plot->Reset();
  for (int i=0;i<xnbin;i++){
    Plot->Fill(Plot->GetBinCenter(i+1),v[i]);
  }
  delete[] v;

  hup  = Plot->GetXaxis()->GetXmax();
  hlow = Plot->GetXaxis()->GetXmin();
  if (x_low==0 && x_up==0 && y_low==0 && y_up == 0){
    x_low  = hlow;
    x_up   = hup;
    y_low  = 0.0;
    double rms  = Plot->GetRMS(); 
    y_up  = rms*3;
  } 
  if (nx==0 && ny==0){
    nx=100;
    ny=20;
  }

  Plot2D = new TH2D("myfunc","",nx,x_low,x_up,ny,y_low,y_up);
  double stepx=(x_up-x_low)/nx,stepy=(y_up-y_low)/ny;
  double x,y;


  for (int i=0;i<nx;i++){
    x = (i+0.5)*stepx+x_low;
    for (int j=0;j<ny;j++){
      y = (j+0.5)*stepy+y_low;
      Plot2D->Fill(x,y,f2(x,y));
    }
  }

  SetStyle();  

  Plot2D->SetTitle("");
  SetXTitle(Plot2D,plot);

  //TColor col85(65,0.00,0.75,1.00);
  //TColor col84(64,0.12,0.56,1.00);

  /*Set by Maurizio
//Set by Guido
  TColor *col7 = gROOT->GetColor(7);
  col7->SetRGB(0.41569,0.95294,0.95686);
  TColor* col65 = gROOT->GetColor(65);
  col65->SetRGB(0.57647,0.87451,0.23922);
  TColor* col64 = gROOT->GetColor(64);
  col64->SetRGB(0.95294,0.86667,0.29020);
  TColor* col4 = gROOT->GetColor(4);
  col4->SetRGB(0.96078,0.64706,0.17255);
  TColor* col2 = gROOT->GetColor(2);
  col2->SetRGB(0.92941,0.21569,0.20000);
  */
  
  /* Set by Maurizio
//Option 1
  TColor *col7 = gROOT->GetColor(7);
  col7->SetRGB(224/255.,243/255.,248/255.);
  TColor* col65 = gROOT->GetColor(65);
  col65->SetRGB(145/255.,191/255.,219/255.);
  TColor* col64 = gROOT->GetColor(64);
  col64->SetRGB(69/255.,117/255.,180/255.);
  TColor* col4 = gROOT->GetColor(4);
  col4->SetRGB(254/255.,224/255.,144/255.);
  TColor* col2 = gROOT->GetColor(2);
  col2->SetRGB(252/255.,141/255.,89/255.);
  */
  

//Option 2
  TColor *col7 = gROOT->GetColor(7);
  col7->SetRGB(224/255.,243/255.,248/255.);
  TColor* col65 = gROOT->GetColor(65);
  col65->SetRGB(145/255.,191/255.,219/255.);
  TColor* col64 = gROOT->GetColor(64);
  col64->SetRGB(254/255.,224/255.,144/255.);
  TColor* col4 = gROOT->GetColor(4);
  col4->SetRGB(252/255.,141/255.,89/255.);
  TColor* col2 = gROOT->GetColor(2);
  col2->SetRGB(215/255.,48/255.,39/255.);

  
/*Option 4
  TColor *col7 = gROOT->GetColor(7);
  col7->SetRGB(127/255.,191/255.,123/255.);
  TColor* col65 = gROOT->GetColor(65);
  col65->SetRGB(145/255.,191/255.,219/255.);
  TColor* col64 = gROOT->GetColor(64);
  col64->SetRGB(254/255.,224/255.,144/255.);
  TColor* col4 = gROOT->GetColor(4);
  col4->SetRGB(252/255.,141/255.,89/255.);
  TColor* col2 = gROOT->GetColor(2);
  col2->SetRGB(215/255.,48/255.,39/255.);
*/


  Plot2D->SetMaximum(6.0);
  Plot2D->SetMinimum(0.0);
  int colors[6]={0,7,65,64,4,2};
  
  Plot2D->SetMaximum(6.0);
  Plot2D->SetContour(6);
  gStyle->SetPalette(6,colors);
  gPad->SetGrid();
  gROOT->ForceStyle();

  Plot2D->Draw("colz");
  Plot2D->SetContour(6);
  Plot2D->SetLineWidth(3);
  Plot2D->SetLineColor(1);
  Plot2D->Draw("cont3 same");

  //Plot2D->Draw("COLZ");

  TLatex t;
  t.SetTextSize(0.06);
  double xadd = (x_up-x_low)*0.15;
  double yadd = (y_up-y_low)*0.01;
  t.DrawLatex(x_up+xadd,
              y_up-yadd,"#sigma");

  if (xval!=-999){
    double lw = (y_up-y_low)/40;
    TLine *lx = new TLine();
    lx->SetLineWidth(4);

    cout << xval-0.1*(y_up-y_low) << endl;

    lx->DrawLine(xval-0.025*(x_up-x_low),xerr,xval+0.025*(x_up-x_low),xerr);
    lx->DrawLine(xval,xerr+0.025*(y_up-y_low),xval,xerr-0.025*(y_up-y_low));
        
  }

  if (xval2!=-999){
    double lw = (y_up-y_low)/40;
    TLine *ly = new TLine();
    ly->SetLineWidth(4);

    ly->DrawLine(xval2-0.025*(x_up-x_low)/sqrt(2.),xerr2-0.025*(y_up-y_low)/sqrt(2.),
		 xval2+0.025*(x_up-x_low)/sqrt(2.),xerr2+0.025*(y_up-y_low)/sqrt(2.));
    ly->DrawLine(xval2-0.025*(x_up-x_low)/sqrt(2.),xerr2+0.025*(y_up-y_low)/sqrt(2.),
		 xval2+0.025*(x_up-x_low)/sqrt(2.),xerr2-0.025*(y_up-y_low)/sqrt(2.));
    //ly->DrawLine(xval2-0.025*(x_up-x_low),xerr2,xval2+0.025*(x_up-x_low),xerr2);
    //ly->DrawLine(xval2,xerr2+0.025*(y_up-y_low),xval2,xerr2-0.025*(y_up-y_low));
    
  }

  if (xval3!=-999){
    double lw = (y_up-y_low)/40;
    TLine *ly = new TLine();
    ly->SetLineWidth(4);

    ly->DrawLine(xval3-0.025*(x_up-x_low)/sqrt(2.),xerr3-0.025*(y_up-y_low)/sqrt(2.),
		 xval3+0.025*(x_up-x_low)/sqrt(2.),xerr3+0.025*(y_up-y_low)/sqrt(2.));
    ly->DrawLine(xval3-0.025*(x_up-x_low)/sqrt(2.),xerr3+0.025*(y_up-y_low)/sqrt(2.),
		 xval3+0.025*(x_up-x_low)/sqrt(2.),xerr3-0.025*(y_up-y_low)/sqrt(2.));
    ly->DrawLine(xval3-0.025*(x_up-x_low),xerr3,xval3+0.025*(x_up-x_low),xerr3);
    ly->DrawLine(xval3,xerr3+0.025*(y_up-y_low),xval3,xerr3-0.025*(y_up-y_low));
    
  }

  if(llogo) PutLogo();

  gPad->Modified();
  gPad->Update();

  /*
  tela->cd();

  season="summer22";
  //double yl=0.77, yu=0.88, xl=0.21, xu=0.39;
  double yl=0.77, yu=0.88, xl=0.61, xu=0.79;
  double newyl=yl-(yu-yl);
  double newxl=xl-(xu-xl)/2;
  double newyu=yu;
  double newxu=xu+(xu-xl)/2;

  TPad *padSeason = new TPad("padSeason","padSeason",newxl,newyl,newxu,newyu,0,0,0);
  padSeason->SetFillStyle(4000);
  padSeason->Draw();
  padSeason->cd();
  TLatex* tex = new TLatex;
  TColor *col = new TColor(999,42/255.,44/255.,114/255.);
  tex->SetTextSize(0.16);
  tex->SetTextAlign(21);
  tex->SetTextColor(602);
  tex->DrawTextNDC(0.49,0.43, season);

  addtxt="SM fit";
  //addtxt="NP fit";
  if (addtxt!=""){
    TLatex* addtex = new TLatex;
    addtex->SetTextSize(0.18);
    addtex->SetTextAlign(21);
    //addtex->DrawTextNDC(0.38,0.1, addtxt);
    addtex->DrawTextNDC(0.57,0.1, addtxt);
  }
  */
  
  return 0;

}

int main(int argc, char *argv[]){



  
  theapp = new TApplication("app",0,NULL);

  if (argc<4 && argc!=3){
    cout << "####################################################"   << endl;
    cout << "%makepull:                                          "   << endl;
    cout << " - needed parameters                                "   << endl;
    cout << "  'filename' (w/o extension - .root assumed)        "   << endl;
    cout << "  'plotname' (e.g. sin2b   )                        "   << endl;
    cout << "     for files converted from paw use h232=\"sin2b\""  << endl;              
    cout << "  'fitname-directory'  (=none no dir struct.)"   << endl;
    cout << endl;
    cout << " e.g. " << endl;
    cout << "  > makepull output sin2b fit1           "   << endl;
    cout << "    read from file output.root and pull for the histogram " << endl;
    cout << "    for sin2b resident in the fit directory fit1      " << endl;
    cout << endl;
    cout << " - to check the containts of the filename type " << endl;
    cout << "  > make1Dplots filename (dir) info"    << endl;
    cout << endl;
    cout << " - optional parameters                           " << endl;
    cout << "   --makeps  -> print .eps file     " << endl;
    cout << "   --nologo  -> don't show UTfit logo     " << endl;
    cout << "   -val=xval -> xval is the measured point (meas.) to superimpose " << endl;
    cout << "   -err=xerr -> xerr is the error of xval                " << endl;
    cout << "   -range=``[xmin,xmax]x[ymin,ymax]'' -> define the graph range    "   <<endl;
    cout << "   -bins=``[xbins]x[ybins]'' -> define the graph binning    "   <<endl;
    cout << "   -logo=position  -> 1=bottom-left,          2=top-left" << endl;
    cout << "                      3=top-right(default),   4=bottom-right" << endl;
    cout << endl;
    cout << " How to use it in a interactive root session:  " << endl;
    cout << "  root [0] .L makepull.cc                    " << endl;
    cout << "  root [1] makepull(filename,plotname,fitname,range) "<<endl;
    cout << "####################################################"   << endl;

    cout << "   -season=winter10  -> attach the date information" << endl;
    cout << "   -addtext=test  -> attach additional information" << endl;

    return 0;
  }

  
  TString srange="null";
  for (int i=1;i<argc;i++){
    if (strncmp(argv[i],"-",1)==0){
      if (strncmp(argv[i],"--makeps",8)==0)  lmakeps  = true;
      if (strncmp(argv[i],"--nologo",8)==0)  llogo  = false;
      if (strncmp(argv[i],"-val",4)==0) sscanf(argv[i],"-val=%lf",&xval);
      if (strncmp(argv[i],"-err",4)==0) sscanf(argv[i],"-err=%lf",&xerr);
      if (strncmp(argv[i],"-skval",5)==0) sscanf(argv[i],"-skval=%lf",&xval2);
      if (strncmp(argv[i],"-skerr",5)==0) sscanf(argv[i],"-skerr=%lf",&xerr2);
      if (strncmp(argv[i],"-skval2",6)==0) sscanf(argv[i],"-skval2=%lf",&xval3);
      if (strncmp(argv[i],"-skerr2",6)==0) sscanf(argv[i],"-skerr2=%lf",&xerr3);
      if (strncmp(argv[i],"-logo",5)==0) sscanf(argv[i],"-logo=%d",&logoposition);
      if (strncmp(argv[i],"-range",6)==0){
        TString stmp(argv[i]+7);
        sscanf(stmp.Data(),"[%lf,%lf]x[%lf,%lf]",&x_low,&x_up,&y_low,&y_up);    
      }
      if (strncmp(argv[i],"-bins",5)==0){
        TString stmp(argv[i]+6);
        sscanf(stmp.Data(),"[%d]x[%d]",&nx,&ny);    
      }
      char str[100];
      if (strncmp(argv[i],"-xlab",5)==0){
	sscanf(argv[i],"-xlab=%s",str);
	TString sx(str,strlen(str));
	xlab  = sx;
	lxlab = true;
      }

      char inputstr[100];

      if (strncmp(argv[i],"-season",7)==0){
        sscanf(argv[i],"-season=%s",inputstr);
        TString label(inputstr,strlen(inputstr));
        season = label;
      }

      if (strncmp(argv[i],"-addtext",8)==0){
        string str(argv[i]);
        addtxt = str.substr(9,str.length()-1);
      }
      
    }
  }


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

  int ret=makepull(sfile,plot,fit);

  if (ret!=0) return -1;

  cout <<"%makepull: select Quit ROOT from the menu file to continue " << endl;
  theapp->Run(true);

  if (lmakeps){
    int ind=-1;
    for(int i=0;i<plot.Length();i++){
      ind = plot.Index("=");
    }
    if (ind!=-1) {
      plot = plot(ind+1,plot.Length());
    }
    TString psfile = "pull_"+plot+".eps";
    tela->Print(psfile);
  }
  if (datafile!=NULL) datafile->Close();
  tela->Close();

  return 0;

}


double fzero(double x, double area){
  //  TMath m;
  return 2*area-TMath::Erfc(x/sqrt(2.0));
}

double finverfc(double x){

  double w=3,eps=0.001;
  double diff,val,s=0.0;
  
  if (x<=1e-100) return 20.00;
  if (x>=0.5) return 0.00;

  do {
    if (fzero(s,x)==0)   val = s;
    if (fzero(s+w,x)==0) val = s+w;
    diff = fzero(s,x)*fzero(s+w,x);
    if (diff<0){
      w/=2;
      val = s+w; 
    } else if (diff>0) {
      s+=w;
    }
  } while (w>eps || diff == 0);

  return val;

}

