using namespace std;
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TImage.h>
#include <TASImage.h>
#include <TColor.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TLatex.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMarker.h>
#include <TApplication.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include "analysis.hh"

double Pi = acos(-1.);
double xval=-999,xerr;
double yval=-999,yerr;

void DrawImage(double xl, double yl, double xu, double yu){

  TASImage *img = new TASImage("../common/logo.png");
  if (nocol) 
    img = new TASImage("../common/logo_bandw.png");
  TPad *pad = new TPad("pad","pad",xl,yl,xu,yu,0,0,0);
  img->SetConstRatio(kFALSE);
  img->SetImageQuality(TAttImage::kImgBest);
  c1->cd();
  pad->Draw();
  pad->cd();
  img->Draw();
}

void PutLogo(int logoposition, int zoom, int doNP){

  double xl = 0.195;
  if (zoom!=0) xl = 0.180;
  else if (doNP == 1 ) xl = 0.260;
  double yl = 0.270;
  if (zoom!=0) yl = 0.255;

  if (logoposition==3 || logoposition==4) {
    xl = 0.735;
    if (zoom!=0) xl = 0.140;
    else if (doNP == 1) xl = 0.750;
  }
  if (logoposition==2 || logoposition==3) {
    yl = 0.770;
    if (zoom!=0) yl = 0.790;
    else if (doNP == 1) yl = 0.800;
    else if (lsquare == 1) yl = 0.790;
  }
  double dx = (doNP == 1 || zoom!=0 || lsquare ==1 ? 0.095 : 0.08);
  double dy = (doNP == 1 || zoom!=0 || lsquare ==1 ? 0.065 : 0.08);
  double scaling = 1.5;
  if (zoom!=0) scaling = 1.6;
  else if (doNP == 1 || lsquare ==1) scaling = 1.32;
  double xu = xl+dx*scaling;
  double yu = yl+dy*scaling;
 
 //  TPad *pad = new TPad("pad","pad",xl,yl,xu,yu);
  
//    double x0 = 0.180;
//     double y0 = 0.160;
//     if (logoposition==3 || logoposition==4)
//       x0 = 0.780; //right part of the plot
//     if (logoposition==2 || logoposition==3)
//       y0 = 0.820; //top part of the plot
//     double x1 = x0+(lsquare == 1? 0.095 :0.11);
//     double y1 = y0+(lsquare == 1? 0.075 :0.11);
//     TPad *pad = new TPad("pad","pad",x0,y0,x1,y1);
 
  
//   double xc;
//   int    nbin = PlotToBeDrawn->GetNbinsX();
//   double hmax = PlotToBeDrawn->GetBinContent(1);
//   double hmin = PlotToBeDrawn->GetBinContent(1);
  
//   for (int i=0;i<nbin;i++){
//     if (PlotToBeDrawn->GetBinContent(i+1)>hmax){
//       hmax = PlotToBeDrawn->GetBinContent(i+1);
//       xc   = PlotToBeDrawn->GetBinCenter(i+1);
//     }
//     if (PlotToBeDrawn->GetBinContent(i+1)<hmin) hmin = PlotToBeDrawn->GetBinContent(1);
//   }
  
//   double xmax = PlotToBeDrawn->GetXaxis()->GetXmax();
//   double xmin = PlotToBeDrawn->GetXaxis()->GetXmin();
    
//   bool left=true,up=false;   
//   if (xc<(xmax+xmin)/2) left=false;
//   if (fabs(hmax-hmin)/hmax<0.2) up = true;
  
//   if (!up && !left){
//     xl=.70; xu=.90; 
//   }
//   if (up){
//     PlotToBeDrawn->SetMaximum(hmax/0.75);
//   }
  
  DrawImage(xl,yl,xu,yu);
  
  gPad->Update();
  c1->cd();

  if(markSM) {
    //TMarker* mark = new TMarker(1., 1., 34);
    //TMarker* mark = new TMarker(1., 0., 34);
    TMarker* mark = new TMarker(0.0001, 0., 34);
    mark->SetMarkerSize(2.0);
    mark->SetMarkerColor(2);
    mark->Draw("SAME");
  }
  
  //season="NP fit";
  season="summer25";
  if (season!=""){  
    double newyl=yl-(yu-yl)*1.2;
    double newxl=xl-(xu-xl)*1.9;
    double newyu=yu;
    double newxu=xu+(xu-xl)/1.9;

    TPad *padSeason = new TPad("padSeason","padSeason",newxl,newyl,newxu,newyu,0,0,0);
        
    padSeason->SetFillStyle(4000);
    padSeason->Draw();
    padSeason->cd();

    TLatex *tex = new TLatex;
    //TColor *col = new TColor(999,42/255.,44/255.,114/255.);
    if (zoom != 0) tex->SetTextSize(0.03);
    else tex->SetTextSize(0.16);
    tex->SetTextAlign(21);
    tex->SetTextColor(602);

    if (zoom !=0) {
      if (logoposition==1 || logoposition==2)
	tex->DrawTextNDC(0.255,0.78, season);
      else 
	tex->DrawTextNDC(0.55,0.48, season);
    } else {
      if (logoposition==1 || logoposition==2){
	tex->DrawTextNDC(0.70,0.48, season);
      } else {
	tex->DrawTextNDC(0.69,0.48, season);
      }
    }
    
    addtxt = "SM fit";
    //addtxt = "CKM24";
    //addtxt = "sides&#epsilon_{K}";
    //addtxt = "CP conserving";
    addtxt2 = "";
    addtxt3 = "";
    if (addtxt!=""){

      TLatex* addtex = new TLatex;
      TLatex* addtex2 = new TLatex;
      TLatex* addtex3 = new TLatex;
      addtex->SetTextSize(0.18);
      addtex2->SetTextSize(0.18);
      addtex3->SetTextSize(0.18);
      addtex->SetTextColor(2);
      addtex3->SetTextColor(802);


      if (zoom !=0) {
	addtex->SetTextSize(0.03);
	addtex->SetTextColor(2);
	if (logoposition==1 || logoposition==2)
	  addtex->DrawLatex(0.30,0.73, addtxt);
	else 
	  addtex->DrawLatex(0.75,0.43, addtxt);
	//addtex->DrawTextNDC(0.75,0.43, addtxt);
      } else {
	if (logoposition==1 || logoposition==2){
	  addtex->SetTextAlign(11);
	  addtex->DrawLatex(0.8,0.20, addtxt);
	} else {
	  addtex->SetTextAlign(21);
	  addtex->DrawLatex(0.8,0.20, addtxt);
	  //addtex->DrawTextNDC(0.69,0.1, addtxt);
	  addtex2->DrawTextNDC(0.36,0.1, addtxt2);
	  addtex3->DrawTextNDC(0.01,0.1, addtxt3);
	}
      }

    }
    c1->cd();
  }
  
}

int main(int argc, char *argv[]){

  if ( argc<3) {
    cout << "####################################################"   << endl;
    cout << "%2Dplot:                                            "   << endl;
    cout << " - needed parameters                                "   << endl;
    cout << "  'filename' (w/o extension - .root assumed)        "   << endl;
    cout << "  'confname' (w/o extension - .conf assumed)        "   << endl;
    cout << "  'output file name'  (none -> same than conf file) "   << endl;
    cout << endl;
    cout << " e.g. " << endl;
    cout << "  > 2Dplot output analysis outplot               " << endl;
    cout << "    read from file output.root and make the plot " << endl;
    cout << "    as specified in analysis in plots/outplot.eps" << endl;
    cout << endl;
    cout << " you can also plot a single histo giving:        " << endl;
    cout << "  'filename' (w/o extension - .root assumed)     "   << endl;
    cout << "  'plotname' (e.g. sin2b   )                     "   << endl;
    cout << "  'fitname-directory'                            "   << endl;
    cout << "  'output file name'                             "   << endl;
    cout << " - optional parameters                           " << endl;
    cout << "   -xlab=namex  -> x label (default read from histogram)" << endl;
    cout << "   -ylab=namey  -> y label (default 'Probability density')" << endl;
    cout << "   -col1=index  -> color index for the 1 sigma interval" << endl;
    cout << "   -col2=index  -> color index for the 2 sigma interval" << endl;
    cout << "   -logo=position  -> 1=bottom-left(default), 2=top-left" << endl;
    cout << "                      3=top-right,            4=bottom-right" << endl;
    cout << "   --square     -> square shaped canvas" << endl;
    cout << "   --nologo     -> do not plot UTfit logo" << endl;
    cout << "   --drawlines  -> draw contour lines" << endl;
    cout << "   -season=winter10  -> attach the date information" << endl;
    cout << "   -addtext=test  -> attach additional information" << endl;
    cout << "   -plot2=nameplot2  -> plot to overimpose" << endl;
    cout << "   -dir2=namedir2  -> dir of the plot to overimpose" << endl;
    cout << "   -col3=index  -> color index for the 1 sigma interval" << endl;
    cout << "   -col4=index  -> color index for the 2 sigma interval" << endl;
    cout << "   -smooth=ntime -> smooth (ntime = number of iterative smoothing)" << endl;
    cout << endl;
    cout << "####################################################"   << endl;

    return 0;  
  }

  TString sfile((const char *)argv[1],strlen(argv[1]));
  TString input_dot_root = sfile + ".root";
  datafile= new TFile(input_dot_root);
  datafile2= new TFile("JpsiphiErrorStandard.root");

  TString conffile((const char *)argv[2],strlen(argv[2]));  
  TString outfile = conffile;

  if(argc == 4 ) {
    const int n = argc -1;
    outfile = TString((const char *)argv[n],strlen(argv[n]));
  }

  if(argc >= 5) {
    const int n = 4;
    outfile = TString((const char *)argv[n],strlen(argv[n]));
  }

  char inputstr[100];
  logoposition=1;
  int other=0;
  if(argc >= 5) {
    plotname = TString((const char *)argv[2],strlen(argv[2]));
    dirname  = TString((const char *)argv[3],strlen(argv[3]));
    for (int i=1;i<argc;i++){
      if (strncmp(argv[i],"-",1)==0){
	if (strncmp(argv[i],"-xval",4)==0) sscanf(argv[i],"-xval=%lf",&xval);
	if (strncmp(argv[i],"-xerr",4)==0) sscanf(argv[i],"-xerr=%lf",&xerr);
	if (strncmp(argv[i],"-yval",4)==0) sscanf(argv[i],"-yval=%lf",&yval);
	if (strncmp(argv[i],"-yerr",4)==0) sscanf(argv[i],"-yerr=%lf",&yerr);
	if (strncmp(argv[i],"-smooth",7)==0) sscanf(argv[i],"-smooth=%d",&smooth);
      }
    }
    graphic = 1;
    makeps = 1;
    makepdf = 0;
    other=1;
    col1=4;
    col2=5;
    col3=3;
    col4=4;
    col5=3;
    col6=4;
    col7=3;
    col8=4;
    for (int i=1;i<argc;i++){
      if (strncmp(argv[i],"-xlab",5)==0){
	sscanf(argv[i],"-xlab=%s",inputstr);
	TString sx(inputstr,strlen(inputstr));
	//	xlab  = sx;
	char tmp[40];
	char sx1[40];
	char inputstr[40];
	snprintf(tmp, sizeof(tmp), "%s", inputstr);
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
	sscanf(argv[i],"-ylab=%s",inputstr);
	TString sy(inputstr,strlen(inputstr));
	//	ylab  = sy;
	char tmp[40];
	char sx1[40];
	char inputstr[40];
	snprintf(tmp, sizeof(tmp), "%s", inputstr);
	for (int c=1; c<=40; c++) {
	  if (tmp[c]=='_')
	    if (tmp[c+1]!='{')
	      tmp[c]=' ';
	}
	sprintf(sx1,"%s",tmp);
	ylab  = sx1;
	
	lylab = true;
      }      
      if (strncmp(argv[i],"-plot2",6)==0){
	sscanf(argv[i],"-plot2=%s",inputstr);
	TString p2(inputstr,strlen(inputstr));
	plotname2  = p2;
	lplot2 = true;
      }
      if (strncmp(argv[i],"-plot3",6)==0){
	sscanf(argv[i],"-plot3=%s",inputstr);
	TString p3(inputstr,strlen(inputstr));
	plotname3  = p3;
	lplot3 = true;
      }
      if (strncmp(argv[i],"-plot4",6)==0){
	sscanf(argv[i],"-plot4=%s",inputstr);
	TString p4(inputstr,strlen(inputstr));
	plotname4  = p4;
	lplot4 = true;
      }

      if (strncmp(argv[i],"-dir2",5)==0){
	sscanf(argv[i],"-dir2=%s",inputstr);
	TString d2(inputstr,strlen(inputstr));
	dirname2  = d2;
      }
      if (strncmp(argv[i],"-dir3",5)==0){
	sscanf(argv[i],"-dir3=%s",inputstr);
	TString d3(inputstr,strlen(inputstr));
	dirname3  = d3;
      }
      if (strncmp(argv[i],"-dir4",5)==0){
	sscanf(argv[i],"-dir4=%s",inputstr);
	TString d4(inputstr,strlen(inputstr));
	dirname4  = d4;
      }
      
      if (strncmp(argv[i],"-season",7)==0){
	sscanf(argv[i],"-season=%s",inputstr);
	TString label(inputstr,strlen(inputstr));
	season = label;
      }
      
      if (strncmp(argv[i],"-addtext",8)==0){
	string str(argv[i]);
	addtxt = str.substr(9,str.length()-1);
      }

      if (strncmp(argv[i],"-col1",5)==0) sscanf(argv[i],"-col1=%d",&col1);
      if (strncmp(argv[i],"-col2",5)==0) sscanf(argv[i],"-col2=%d",&col2);
      if (strncmp(argv[i],"-col3",5)==0) sscanf(argv[i],"-col3=%d",&col3);
      if (strncmp(argv[i],"-col4",5)==0) sscanf(argv[i],"-col4=%d",&col4);
      if (strncmp(argv[i],"-col5",5)==0) sscanf(argv[i],"-col5=%d",&col5);
      if (strncmp(argv[i],"-col6",5)==0) sscanf(argv[i],"-col6=%d",&col6);
      if (strncmp(argv[i],"-col7",5)==0) sscanf(argv[i],"-col7=%d",&col7);
      if (strncmp(argv[i],"-col8",5)==0) sscanf(argv[i],"-col8=%d",&col8);
      if (strncmp(argv[i],"-logo",5)==0) sscanf(argv[i],"-logo=%d",&logoposition);
      if (strncmp(argv[i],"--square",5)==0) lsquare = true;
      if (strncmp(argv[i],"--nologo",8)==0) llogo = false;
      if (strncmp(argv[i],"--drawlines",11)==0) llines = true;
      if (strncmp(argv[i],"--nocol",8)==0) nocol = true;
      if (strncmp(argv[i],"--markSM",9)==0) markSM = true;
    }
  } else {  
    const int namelenght = strlen(argv[2])+5;
    char filename[namelenght];  
    sprintf(filename,conffile + ".conf",conffile + ".conf");
    Read_Parameters(filename);
    if(zoom!=0) {
      plotname = "etavsrhozoom";
      xmax=1.25;
      xmin=-0.25;
      ymax=1.25;
      ymin=0.0;
      etamin=0.0;
    }
  }

  gStyle = new TStyle;
  DefineColors();
  SetDefaultPalette(other);
  gStyle->SetOptStat(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetLineWidth(2);
  //  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);

  cout << "graphic " << graphic << endl;
  if (graphic){
    theapp = new TApplication("app",0,NULL);
  }

  if (zoom != 0) {
    c1 =  new TCanvas("his2D","2D histograms",600,500);
    c1->SetLeftMargin(0.16);
  } else  if( (doNP ==1 && plotname == "etavsrho") || lsquare) {
    c1 =  new TCanvas("his2D","2D histograms",600,600);
    c1->SetLeftMargin(0.20);
    //    c1->SetRightMargin(0.08);
  } else {
    c1 =  new TCanvas("his2D","2D histograms",600,352);
    c1->SetLeftMargin(0.12);
  }

  c1->SetFillColor(0);
  c1->SetBorderMode(0);

  char number[15];
  sprintf(number, "fit%d",fitnumber);
  string dir = number;

  if(plotname == "etavsrho" || (zoom !=0 && plotname == "etavsrhozoom")) {
    makeCKMplot(dir);
  } else {
    c1->SetBottomMargin(0.15);
    /////////////////////////////////////////////////////////////////
    // calculate all the TGraphs, that go as argument of makeOtherPlot
    // to move to a subroutine later?????????
    /////////////////////////////////////////////////////////////////
    // get the levels
    get_data(plotname.c_str(),dirname.c_str(),2);
    //TH2D* Constrain_1 = Constrain;
    TH2D* Constrain_1 = (TH2D*)Constrain->Clone();
    Constrain_1->SetName("Constrain_1"); 
    double level68_1 = level[4];
    double level95_1 = level[3];
    //TObjArray* histo_1 = GraphFromHisto(Constrain_1,level68_1,level95_1);
    TObjArray* histo_1 = (TObjArray*)GraphFromHisto(Constrain_1,level68_1,level95_1)->Clone();
    
    TH2D* Constrain_2 = NULL;
    TObjArray* histo_2 = NULL;
    double level68_2 = 0;
    double level95_2 = 0;
    if (lplot2) {
      // get the levels of the second plot
      get_data(plotname2.c_str(),dirname2.c_str(),2);
      //get_data(plotname2.c_str(),dirname2.c_str(),2,2);
      Constrain_2 = (TH2D*)Constrain->Clone();
      Constrain_2->SetName("Constrain_2"); 
      level68_2 = level[4];
      level95_2 = level[3];
      histo_2 = (TObjArray*)GraphFromHisto(Constrain_2,level68_2,level95_2)->Clone();
    }
    TH2D* Constrain_3 = NULL;
    TObjArray* histo_3 = NULL;
    double level68_3 = 0;
    double level95_3 = 0;
    if (lplot3) {
      // get the levels of the second plot
      get_data(plotname3.c_str(),dirname3.c_str(),2,2);
      Constrain_3 = (TH2D*)Constrain->Clone();
      Constrain_3->SetName("Constrain_3"); 
      level68_3 = level[4];
      level95_3 = level[3];
      histo_3 = (TObjArray*)GraphFromHisto(Constrain_3,level68_3,level95_3)->Clone();
    }

    TH2D* Constrain_4 = NULL;
    TObjArray* histo_4 = NULL;
    double level68_4 = 0;
    double level95_4 = 0;
    if (lplot4) {
      // get the levels of the second plot
      get_data(plotname4.c_str(),dirname4.c_str(),2,2);
      Constrain_4 = (TH2D*)Constrain->Clone();
      Constrain_4->SetName("Constrain_4"); 
      level68_4 = level[4];
      level95_4 = level[3];
      histo_4 = (TObjArray*)GraphFromHisto(Constrain_4,level68_4,level95_4)->Clone();
    }

    //plot
    if (lplot2&&lplot3&&lplot4) {
      makeOtherPlot(Constrain_4,histo_4,0,col7,col8);
      makeOtherPlot(Constrain_3,histo_3,1,col5,col6);
      makeOtherPlot(Constrain_2,histo_2,1,col3,col4);
      makeOtherPlot(Constrain_1,histo_1,1,col1,col2);
    } else if (lplot2&&lplot3) {
      makeOtherPlot(Constrain_3,histo_3,0,col5,col6);
      makeOtherPlot(Constrain_2,histo_2,1,col3,col4);
      makeOtherPlot(Constrain_1,histo_1,1,col1,col2);
    } else if (lplot2) {
      makeOtherPlot(Constrain_2,histo_2,0,col3,col4);
      makeOtherPlot(Constrain_1,histo_1,1,col1,col2);
    } else {
      makeOtherPlot(Constrain_1,histo_1,0,col1,col2);
    }
        
  }

  if (graphic && !makeps && !makepdf ){
    cout <<"Select Quit ROOT from the menu File to continue " << endl;
    theapp->Run(true);
  }

  TFile *output = new TFile("dump.root","RECREATE");
  if (makeps){
    TString output_dot_eps = "plots/" + outfile + ".eps";
    TString output_dot_root = "plots/" + outfile + ".root";
    c1->Print(output_dot_eps);
    c1->Print(output_dot_root);
    // output_dot_eps = "plots/" + outfile + ".C";
    // c1->SaveAs(output_dot_eps);
    //    TString output_dot_gif = "plots/" + outfile + ".gif";
    //    c1->Print(output_dot_gif);

    //TFile *output = new TFile("dump.root","RECREATE");
    output->cd();
    output->mkdir("dumped");
    output->cd("dumped");
    c1->Write();
  }
  output->Close();

  if (makepdf){
    TString output_dot_pdf = "plots/" + outfile + ".pdf";
    c1->Print(output_dot_pdf);
    // output_dot_eps = "plots/" + outfile + ".C";
    // c1->SaveAs(output_dot_eps);
    //    TString output_dot_gif = "plots/" + outfile + ".gif";
    //    c1->Print(output_dot_gif);

    TFile *output = new TFile("dump.root","RECREATE");
    output->cd();
    output->mkdir("dumped");
    output->cd("dumped");
    c1->Write();
    output->Close();

  }

  if (graphic) c1->Close();
  return 0;
}

void makeCKMplot(string dir) {  
  set_data(dir);
  prepare_epsk();
  prepare_Vub();
  prepare_b2taunu();
  prepare_dmd();
  prepare_dms();
  prepare_k2pnn();
  prepare_btovg();

  plot_rho_eta(dir);
  
  TLatex* tex = new TLatex;
  if(zoom != 0) {
    tex->SetTextSize(0.05);
    tex->SetTextAlign(22);
    if (set_gamma[0]&&set_gamma[5]) tex->DrawLatex(0.47,1.06, "#phi_{3}/#gamma");
    if (set_s2b[0]&&set_s2b[5]) tex->DrawLatex(-0.12, 0.55, "#phi_{1}/#beta");
    if (set_c2b[0]&&set_c2b[5]) tex->DrawLatex(1.12, 0.2, "cos2#beta");
    if (set_D0p0[0]&&set_D0p0[5]) tex->DrawLatex(-0.64, 1.1, "D^{0}#pi^{0}");
    if (set_s2bpg[0]&&set_s2bpg[5]) tex->DrawLatex( 0.68,-0.18, "sin(2#beta+#gamma)");
    if (set_s2a[0]&&set_s2a[5]) tex->DrawLatex(.75, 0.40, "#phi_{2}/#alpha");
    if (set_eps[0]&&set_eps[5]) tex->DrawLatex(0.91, 1.10, "#varepsilon_{K}");

    tex->SetTextSize(0.045);
    if (set_dms[0]&&set_dms[5]) tex->DrawLatex( 1.15, 0.95, "#frac{#Deltam_{d}}{#Deltam_{s}}");
    if (set_dmd[0]&&set_dmd[5]) tex->DrawLatex( 1.0, 0.85, "#Deltam_{d}");
    if (set_vub[0]&&set_vub[5]) tex->DrawLatex(0.34,0.16, "#left|#frac{V_{ub}}{V_{cb}}#right|");
    if (set_b2taunu[0]&&set_b2taunu[5]) tex->DrawLatex(-0.36,-0.5, "BR(B#rightarrow#tau#nu)");
    if (set_btovg[0]&&set_btovg[5]) {
      if(btovg_use_his == 2)
	tex->DrawLatex( 0.79, 0.72, "BR(B^{0}#rightarrow#rho^{0}#gamma)");
      else
	tex->DrawLatex( 0.79, 0.72, "BR(B#rightarrow#rho/#omega#gamma)");
      tex->DrawLatex( 0.79, 0.58, "BR(B#rightarrowK^{*}#gamma)");
      TLine* Axisratio = new TLine(0.58,0.65,1.00,0.65);
      Axisratio->SetLineWidth(1);
      Axisratio->Draw();
    }
  } else if (doNP ==1){
    tex->SetTextSize(0.035);
    tex->SetTextAlign(22);
    if (set_gamma[0]&&set_gamma[5]) tex->DrawLatex(0.3,1.1, "#phi_{3}/#gamma");
    if (set_s2b[0]&&set_s2b[5]) tex->DrawLatex(-0.80, 0.72, "#beta");
    if (set_c2b[0]&&set_c2b[5]) tex->DrawLatex(1.12, 0.2, "cos2#beta");
    if (set_D0p0[0]&&set_D0p0[5]) tex->DrawLatex(-0.64, 1.1, "D^{0}#pi^{0}");
    if (set_s2bpg[0]&&set_s2bpg[5]) tex->DrawLatex( 0.68,-0.18, "sin(2#beta+#gamma)");
    if (set_s2a[0]&&set_s2a[5]) tex->DrawLatex(.8, -0.45, "#phi_{2}/#alpha");
    if (set_eps[0]&&set_eps[5]) tex->DrawLatex(-0.90, 0.22, "#varepsilon_{K}");

    tex->SetTextSize(0.03);
    if (set_dms[0]&&set_dms[5]) tex->DrawLatex( 1.15, 0.95, "#frac{#Deltam_{d}}{#Deltam_{s}}");
    if (set_dmd[0]&&set_dmd[5]) tex->DrawLatex( 0.99, 0.85, "#Deltam_{d}");
    if (set_vub[0]&&set_vub[5]) tex->DrawLatex(-0.38,0.12, "#left|#frac{V_{ub}}{V_{cb}}#right|");
    if (set_b2taunu[0]&&set_b2taunu[5]) tex->DrawLatex(-0.36,-0.5, "BR(B#rightarrow#tau#nu)");
    if (set_btovg[0]&&set_btovg[5]) {
      if(btovg_use_his == 2)
	tex->DrawLatex( 0.79, 0.72, "BR(B^{0}#rightarrow#rho^{0}#gamma)");
      else
	tex->DrawLatex( 0.79, 0.72, "BR(B#rightarrow#rho/#omega#gamma)");
      tex->DrawLatex( 0.79, 0.58, "BR(B#rightarrowK^{*}#gamma)");
      TLine* Axisratio = new TLine(0.58,0.65,1.00,0.65);
      Axisratio->SetLineWidth(1);
      Axisratio->Draw();
    }
  } else {
    tex->SetTextSize(0.063);
    tex->SetTextAlign(22);
    if (set_gamma[0]&&set_gamma[5]) tex->DrawLatex(0.37,1.1, "#phi_{3}/#gamma");
    if (set_s2b[0]&&set_s2b[5]) tex->DrawLatex(-1.05, 0.83, "#beta");
    if (set_c2b[0]&&set_c2b[5]) tex->DrawLatex(1.12, 0.2, "cos2#beta");
    if (set_D0p0[0]&&set_D0p0[5]) tex->DrawLatex(-0.64, 1.1, "D^{0}#pi^{0}");
    if (set_s2bpg[0]&&set_s2bpg[5]) tex->DrawLatex( 0.65,-0.17, "2#beta+#gamma");
    if (set_s2a[0]&&set_s2a[5]) tex->DrawLatex(.8, 0.37, "#phi_{2}/#alpha");
    if (set_eps[0]&&set_eps[5]) tex->DrawLatex(-0.60, 0.25, "#varepsilon_{K}");

    tex->SetTextSize(0.053);
    if (set_dms[0]&&set_dms[5]) tex->DrawLatex( 1.16, 0.95, "#frac{#Deltam_{d}}{#Deltam_{s}}");
    if (set_dmd[0]&&set_dmd[5]) tex->DrawLatex( 0.98, 0.86, "#Deltam_{d}");
    if (set_vub[0]&&set_vub[5]) tex->DrawLatex( -0.36, 0.13, "#left|#frac{V_{ub}}{V_{cb}}#right|");
    if (set_b2taunu[0]&&set_b2taunu[5]) tex->DrawLatex(0.65,-0.2, "BR(B#rightarrow#tau#nu)");

    tex->SetTextSize(0.055);
    if (set_k2pnn[0]&&set_k2pnn[5]) tex->DrawLatex(0.5,1.1,"BR(K^{+}#rightarrow#pi^{+}#nu#bar{#nu})");
    if (set_kl2p0nn[0]&&set_kl2p0nn[5]) tex->DrawLatex(0.8,0.4,"BR(K_{L}#rightarrow#pi^{0}#nu#bar{#nu})");
    if (set_btovg[0]&&set_btovg[5]) {
      //      if(btovg_use_his == 2)
      //	tex->DrawLatex( 0.69, 0.72, "BR(B^{0}#rightarrow#rho^{0}#gamma)");
      //      else
	tex->DrawLatex( 0.69, 0.72, "BR(B#rightarrow#rho/#omega#gamma)");
      tex->DrawLatex( 0.69, 0.58, "BR(B#rightarrowK^{*}#gamma)");
      TLine* Axisratio = new TLine(0.48,0.65,0.90,0.65);
      Axisratio->SetLineWidth(1);
      Axisratio->Draw();
    }
  }


  if(logoposition!=0) {
    PutLogo(logoposition,zoom,doNP);
  }

//   if(logoposition!=0) {
//     TImage *img = TImage::Open("../common/logo.png");
//     img->SetConstRatio(kFALSE);
//     double x0 = (doNP == 1 || zoom!=0 ? 0.170 : 0.125);
//     double y0 = 0.110;
//     if (logoposition==3 || logoposition==4)
//       x0 = (doNP == 1 || zoom!=0 ? 0.800 : 0.785); //right part of the plot
//     if (logoposition==2 || logoposition==3)
//       y0 = (doNP == 1 || zoom!=0 ? 0.820 : 0.770); //top part of the plot
//     double x1 = x0+(doNP == 1 || zoom!=0 ? 0.095 : 0.11);
//     double y1 = y0+(doNP == 1 || zoom!=0 ? 0.075 : 0.11);
//     TPad *pad = new TPad("pad","pad",x0,y0,x1,y1);
//     pad->Draw();
//     pad->cd();
//     img->Draw();
//     gPad->Update();
//   }

  if (nocol) c1->GetCanvas()->SetGrayscale();
  c1->cd();
  c1->RedrawAxis("same");

  // drawing the triangle
  if(drawtriangle == 1) {
    string rhoplot="rho"+dir;
    string etaplot="eta"+dir;
    double drawrho = ((TH1D*)(((TDirectory*) 
 			       datafile->Get(dir.c_str()))->Get(rhoplot.c_str())))->GetMean();
    double draweta = ((TH1D*)(((TDirectory*) 
 			       datafile->Get(dir.c_str()))->Get(etaplot.c_str())))->GetMean();
    
    double xx[3],yy[3];
    xx[0] = 0.; yy[0] = 0.;
    xx[1] = drawrho; yy[1] = draweta;
    TGraph* Rb = new TGraph(2,xx,yy);
    Rb->SetLineWidth(2);
    Rb->Draw();
    xx[0] = 0.; yy[0] = 0.;
    xx[1] = 1.; yy[1] = 0.;
    TGraph* Rc = new TGraph(2,xx,yy);
    Rc->SetLineWidth(2);
    Rc->Draw();
    xx[0] = 1.; yy[0] = 0.;
    xx[1] = drawrho; yy[1] = draweta;
    TGraph* Rt = new TGraph(2,xx,yy);
    Rt->SetLineWidth(2);
    Rt->Draw();
  }
}

void  plot_rho_eta(string dir) {
  etamin = (doNP == 1 ? -1.25 : -0.25);
  if(zoom != 0) etamin = 0.;

  // calculate s2b contour levels
  if (set_s2b[0] == 1) {
    get_data("s2b","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      s2blevel[ind] = level[ind];
    }
    s2bConstrain = Constrain;
  }

  // calculate c2b contour levels
  if (set_c2b[0] == 1) {
    get_data("c2b","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      c2blevel[ind] = level[ind];
    }
    c2bConstrain = Constrain;
  }

  // calculate dmd contour levels
  if (set_dmd[0] == 1) {
    get_data("dmd","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      dmdlevel[ind] = level[ind];
    }
    dmdConstrain = Constrain;
  }

  // calculate dms contour levels
  if (set_dms[0] == 1) {
    get_data("dms","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      dmslevel[ind] = level[ind];
    }
    dmsConstrain = Constrain;
  }

  // calculate vub contour levels
  if (set_vub[0] == 1) {
    get_data("vub","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      vublevel[ind] = level[ind];
    }
    vubConstrain = Constrain;
  }

  // calculate eps contour levels
  if (set_eps[0] == 1) {
    get_data("epsk","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      epslevel[ind] = level[ind];
    }
    epsConstrain = Constrain;
  }

  // calculate s2bpg contour levels
  if (set_s2bpg[0] == 1) {
    TString s2bpg = "2bpg_all";
    if(s2bpg_use_his == 2) s2bpg = "2bpg_D";
    if(s2bpg_use_his == 3) s2bpg = "2bpg_Dstar";
    if(s2bpg_use_his == 4) s2bpg = "2bpg_Drho";
    get_data(s2bpg,"rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      s2bpglevel[ind] = level[ind];
    }
    s2bpgConstrain = Constrain;
  }

  // calculate s2a  contour levels
  if (set_s2a[0] == 1) {
    TString s2a = "s2a";
    if(s2a_use_his == 2) s2a = "pipi_sin2alpha";
    if(s2a_use_his == 3) s2a = "rhorho_sin2alpha";
    if(s2a_use_his == 4) s2a = "rhopi_sin2alpha";
    get_data(s2a,"rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      s2alevel[ind] = level[ind];
    }
    s2aConstrain = Constrain;
  }

  // calculate k2pnn  contour levels
  if (set_k2pnn[0] == 1) {
    get_data("k2pnn","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      k2pnnlevel[ind] = level[ind];
    }
    k2pnnConstrain = Constrain;
  }

  // calculate gamma contour levels
  if (set_gamma[0] == 1) {
    TString gamma = "gamma_all";
    if(gamma_use_his == 2) gamma = "gamma_GLW";
    if(gamma_use_his == 3) gamma = "gamma_ADS";
    if(gamma_use_his == 4) gamma = "gamma_DAL";
    get_data(gamma,"rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      gammalevel[ind] = level[ind];
    }
    gammaConstrain = Constrain;
  }

  // calculate b2taunu contour levels
  if (set_b2taunu[0] == 1) {
    get_data("b2taunu","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      b2taunulevel[ind] = level[ind];
    }
    b2taunuConstrain = Constrain;
  } 

  // calculate btovg contour levels
  if (set_btovg[0] == 1) {
    get_data("BtoVg","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      btovglevel[ind] = level[ind];
    }
    btovgConstrain = Constrain;
  }

  // calculate kl->p0nn contour levels
  if (set_kl2p0nn[0] == 1) {
    get_data("klp0nn","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      kl2p0nnlevel[ind] = level[ind];
    }
    kl2p0nnConstrain = Constrain;
  }

  // calculate BtoKpp contour levels
  if (set_BtoKpp[0] == 1) {
    get_data("BtoKpp","rhoetahisto",2);
    for(int ind=0; ind<5;ind++) {
      BtoKpplevel[ind] = level[ind];
    }
    BtoKppConstrain = Constrain;
  }

  int vubd     = set_vub[0];
  int dmdd     = set_dmd[0];
  int dmsd     = set_dms[0];
  int epsd     = set_eps[0];
  int s2bd     = set_s2b[0];
  int c2bd     = set_c2b[0];
  int s2bpgd   = set_s2bpg[0];
  int s2ad     = set_s2a[0];
  int k2pnnd   = set_k2pnn[0];
  int gammad   = set_gamma[0];
  int b2taunud = set_b2taunu[0];
  int btovgd   = set_btovg[0];
  int D0p0d    = set_D0p0[0];
  int kl2p0nnd = set_kl2p0nn[0];
  int BtoKppd  = set_BtoKpp[0];

  for (int i = 1; i<6; i++) {
    vubsel[i]     = 0;
    dmdsel[i]     = 0;
    dmssel[i]     = 0;
    epssel[i]     = 0;
    s2bsel[i]     = 0;
    c2bsel[i]     = 0;
    s2bpgsel[i]   = 0;
    s2asel[i]     = 0;
    k2pnnsel[i]   = 0;
    gammasel[i]   = 0;
    b2taunusel[i] = 0;
    btovgsel[i]   = 0;
    D0p0sel[i]    = 0;
    kl2p0nnsel[i] = 0;
    BtoKppsel[i] = 0;
  }

  prepareopt(set_vub, vubsel);
  prepareopt(set_dmd, dmdsel);
  prepareopt(set_dms, dmssel);
  prepareopt(set_s2b, s2bsel);
  prepareopt(set_c2b, c2bsel);
  prepareopt(set_s2bpg, s2bpgsel);
  prepareopt(set_s2a, s2asel);
  prepareopt(set_k2pnn, k2pnnsel);
  prepareopt(set_eps, epssel);
  prepareopt(set_gamma, gammasel);
  prepareopt(set_b2taunu, b2taunusel);
  prepareopt(set_btovg, btovgsel);
  prepareopt(set_D0p0, D0p0sel);
  prepareopt(set_kl2p0nn, kl2p0nnsel);
  prepareopt(set_BtoKpp, BtoKppsel);

  // preliminary plots for 2alpha and 2bpg, if needed


  if(set_s2a[0] == 2) {
    s2ahisto1 =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha1")).Clone();
    s2ahisto2 =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha2")).Clone();
    s2ahisto3 =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha3")).Clone();
    s2ahisto4 =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha4")).Clone();
    s2ahisto5 =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha5")).Clone();
    s2ahisto =  (TObjArray*) (*CalcGraph(s2a1D,s2as1,s2as2,"alpha")).Clone();
  }

  if (set_s2bpg[0] == 2) {
    s2bpghisto1 =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg1")).Clone();
    s2bpghisto2 =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg2")).Clone();
    s2bpghisto3 =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg3")).Clone();
    s2bpghisto4 =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg4")).Clone();
    s2bpghisto  =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg")).Clone();
  }

  //if (set_s2bpg[0] == 2) {
  //s2bpghisto =  (TObjArray*) (*CalcGraph(s2bpg1D,s2bpgs1,s2bpgs2,"2bpg")).Clone();
  //}

  TH2D* null(0);
  null = new TH2D("null","",128,xmin,xmax,128,etamin,ymax);
  null->SetXTitle("#bar{#rho}");
  null->SetYTitle("#bar{#eta}");
  null->SetTitleSize(0.06,"X");
  null->SetTitleSize(0.06,"Y");
  null->SetTitleOffset(0.7,"X");
  if(zoom!=0) 
    null->SetTitleOffset(0.8,"Y");
  else
    null->SetTitleOffset(0.6,"Y");
    
  null->SetStats(0);
  null->Draw();

  // Draw constraints as areas
  tria(vubd,dmdd,epsd,dmsd,s2bd,c2bd,s2bpgd,s2ad,k2pnnd,gammad,b2taunud,btovgd,D0p0d,kl2p0nnd,BtoKppd,"AREA");

  string plot=plotname+dir;
  if(drawtotarea == 1) {
    get_data(plot.c_str(),dir.c_str(),2);
    int sel[6];
    sel[0] = 0;
    sel[1] = 0;
    sel[2] = 0;
    sel[3] = 1;
    sel[4] = 1;
    sel[5] = 1;
    setdrawlevel(Constrain,level,sel, "AREA");
  }
  c1->Update();

  // Draw constraints as lines
  tria(vubd,dmdd,epsd,dmsd,s2bd,c2bd,s2bpgd,s2ad,k2pnnd,gammad,b2taunud,btovgd,D0p0d,kl2p0nnd,BtoKppd,"CONT");
  TLine* RefAxis = new TLine(xmin,0.,xmax,0.);
  RefAxis->Draw();
  RefAxis->DrawLine(0.,etamin,0.,ymax);

  int vubt     = set_vub[5];
  int dmdt     = set_dmd[5];
  int dmst     = set_dms[5];
  int epst     = set_eps[5];
  int s2bt     = set_s2b[5];
  int c2bt     = set_c2b[5];
  int s2bpgt   = set_s2b[5];
  int s2at     = set_s2a[5];
  int k2pnnt   = set_k2pnn[5];
  int gammat   = set_gamma[5];
  int b2taunut = set_b2taunu[5];
  int btovgt   = set_btovg[5];
  int D0p0t    = set_D0p0[5];
  int kl2p0nnt = set_kl2p0nn[5];
  int BtoKppt  = set_BtoKpp[5];

  if(drawtotarea == 0) {
    get_data(plot.c_str(),dir.c_str(),1);
  }
  tria_tit(vubt,dmdt,epst,dmst,s2bt,c2bt,s2bpgt,s2at,k2pnnt,gammat,b2taunut,btovgt,D0p0t,kl2p0nnt,BtoKppt);
  null->Draw("same");
  TLine* Bound = new TLine(xmin,ymax,xmax,ymax);
  Bound->Draw();
  Bound->DrawLine(xmax,etamin,xmax,ymax);
}

void tria_tit(int vubt, int dmdt, int epst, int dmst, int s2bt, 
	      int c2bt, int s2bpgt, int s2at, int k2pnnt, int gammat, 
	      int b2taunut, int btovgt, int D0p0t, int kl2p0nnt, int BtoKppt) {
}

void tria(int vubd, int dmdd, int epsd, int dmsd, int s2bd, int c2bd, 
	  int s2bpgd, int s2ad, int k2pnnd, int gammad, int b2taunud, 
	  int btovgd, int D0p0d, int kl2p0nnd, int BtoKppd, const string DrawOpts) {

  // Draw 2bpg
  if (s2bpgd == 1)
    setdrawlevel(s2bpgConstrain,s2bpglevel,s2bpgsel, DrawOpts);
  if (s2bpgd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (s2bpgsel[4] > 0)
        if (s2bpgsel[3] > 0) {
          drawFromGraph(s2bpghisto1, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto2, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto3, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto4, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto1, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto2, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto3, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto4, 1, "AREA", s2bpgsel[3]);
          //        drawFromGraph(s2bpghisto, 0, "AREA", s2bpgsel[4]);
          //        drawFromGraph(s2bpghisto, 1, "AREA", s2bpgsel[3]);
        } else {
          drawFromGraph(s2bpghisto1, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto2, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto3, 0, "AREA", s2bpgsel[4]);
          drawFromGraph(s2bpghisto4, 0, "AREA", s2bpgsel[4]);
          //        drawFromGraph(s2bpghisto, 0, "AREA", s2bpgsel[4]);
        } else if (s2bpgsel[3] > 0) {
          drawFromGraph(s2bpghisto1, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto2, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto3, 1, "AREA", s2bpgsel[3]);
          drawFromGraph(s2bpghisto4, 1, "AREA", s2bpgsel[3]);
          //        drawFromGraph(s2bpghisto, 1, "AREA", s2bpgsel[3]);
        }
    } else if (DrawOpts.compare("CONT")==0) {
      if(s2bpgsel[1] > 0)
        drawFromGraph(s2bpghisto, 0, "CONT", s2bpgsel[1]);
      if(s2bpgsel[0] > 0)
        drawFromGraph(s2bpghisto, 1, "CONT", s2bpgsel[0]);
    }
  }

  // Draw c2b
  if (c2bd == 1) 
    setdrawlevel(c2bConstrain,c2blevel,c2bsel,DrawOpts);
  if (c2bd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (c2bsel[4] > 0)
	if (c2bsel[3] > 0) {
	  draw_cos2b(b1s2_c,c2bsel[4]);
	  draw_cos2b(b1s1_c,c2bsel[3]);
	} else
	  draw_cos2b(b1s2_c,c2bsel[4]);
      else if (c2bsel[3] > 0)
 	draw_cos2b(b1s1_c,c2bsel[3]);
    } else if (DrawOpts.compare("CONT")==0) {
      if (c2bsel[0] > 0)
 	cos2b_lines(b1s1_c,c2bsel[0]);
      if (c2bsel[1] > 0)
 	cos2b_lines(b1s2_c,c2bsel[1]);
    }
  }

  // Draw D0p0
  // questa va sistemata....
  if (D0p0d == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (D0p0sel[4] > 0)
	if (D0p0sel[3] > 0) {
	  D0p0_areas(b1s2mn_d0p0,b1s2mx_d0p0,D0p0sel[4],"LEFT");
	  D0p0_areas(b1s2mn_d0p0,b1s2mx_d0p0,D0p0sel[4],"RIGHT");
	  D0p0_areas(b1s1mn_d0p0,b1s1mx_d0p0,D0p0sel[3],"LEFT");
	  D0p0_areas(b1s1mn_d0p0,b1s1mx_d0p0,D0p0sel[3],"RIGHT");
	} else {
	  D0p0_areas(b1s2mn_d0p0,b1s2mx_d0p0,D0p0sel[4],"LEFT");
	  D0p0_areas(b1s2mn_d0p0,b1s2mx_d0p0,D0p0sel[4],"RIGHT");

	}
      else if (D0p0sel[3] > 0) {
	  D0p0_areas(b1s1mn_d0p0,b1s1mx_d0p0,D0p0sel[3],"LEFT");
	  D0p0_areas(b1s1mn_d0p0,b1s1mx_d0p0,D0p0sel[3],"RIGHT");

      }
    }
    else if (DrawOpts.compare("CONT")==0) {
      if (D0p0sel[0] > 0) {
	D0p0_lines(b1s1mn_d0p0,b1s1mx_d0p0,D0p0sel[0]);
      }
      if (D0p0sel[1] > 0) {
	D0p0_lines(b1s2mn_d0p0,b1s2mx_d0p0,D0p0sel[1]);
      }
    }
  }

  // Draw Gamma
  if (gammad == 1)
    setdrawlevel(gammaConstrain,gammalevel,gammasel, DrawOpts);
  if (gammad == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if(gammasel[4] > 0)
	if(gammasel[3] > 0) {
	  drawgamma(g1s2mn, g1s2mx, gammasel[4]);
	  drawgamma(g1s1mn, g1s1mx, gammasel[3]);
	  if(gamma_use_his == 2) {
	    drawgamma(g1s1mn1, g1s1mx1, gammasel[3]);
	    drawgamma(g1s1mn2, g1s1mx2, gammasel[3]);
	  }
	} else 
	  drawgamma(g1s2mn, g1s2mx, gammasel[4]);
      else if(gammasel[3] > 0)
	drawgamma(g1s1mn, g1s1mx, gammasel[3]);
    } else if (DrawOpts.compare("CONT")==0) {
      if(gammasel[0] > 0) {
	gamma_lines(g1s1mn,g1s1mx,gammasel[0]);
	if(gamma_use_his == 2) {
	  gamma_lines(g1s1mn1,g1s1mx1,gammasel[0]);
	  gamma_lines(g1s1mn2,g1s1mx2,gammasel[0]);
	}
      }
      if(gammasel[1] > 0)
	gamma_lines(g1s2mn,g1s2mx,gammasel[1]);
    }
  }
  
  // Draw 2a
  if (s2ad == 1)
    setdrawlevel(s2aConstrain,s2alevel,s2asel, DrawOpts);
  if (s2ad == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (s2asel[4] > 0)
	if (s2asel[3] > 0) {
	  drawFromGraph(s2ahisto1, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto2, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto3, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto4, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto5, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto1, 1, "AREA", s2asel[3]);	  
	  drawFromGraph(s2ahisto2, 1, "AREA", s2asel[3]);	  
	  drawFromGraph(s2ahisto3, 1, "AREA", s2asel[3]);	  
	  drawFromGraph(s2ahisto4, 1, "AREA", s2asel[3]);	  
	  drawFromGraph(s2ahisto5, 1, "AREA", s2asel[3]);	  
	} else {
	  drawFromGraph(s2ahisto1, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto2, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto3, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto4, 0, "AREA", s2asel[4]);
	  drawFromGraph(s2ahisto5, 0, "AREA", s2asel[4]);
	}
      else if (s2asel[3] > 0) {
	drawFromGraph(s2ahisto1, 1,"AREA", s2asel[3]);
	drawFromGraph(s2ahisto2, 1,"AREA", s2asel[3]);
	drawFromGraph(s2ahisto3, 1,"AREA", s2asel[3]);
	drawFromGraph(s2ahisto4, 1,"AREA", s2asel[3]);
	drawFromGraph(s2ahisto5, 1,"AREA", s2asel[3]);
      }
    } else if (DrawOpts.compare("CONT")==0) {
      if(s2asel[1] > 0) {
	drawFromGraph(s2ahisto, 0, "CONT", s2asel[1]);
      }
      if(s2asel[0] > 0) {
	drawFromGraph(s2ahisto, 1, "CONT", s2asel[0]);
      }
    }
  }

  // Draw kpnn
  if (k2pnnd == 1)
    setdrawlevel(k2pnnConstrain,k2pnnlevel,k2pnnsel, DrawOpts);
  if (k2pnnd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (k2pnnsel[4] > 0)
	if (k2pnnsel[3] > 0) {
	  drawfunc(rho2_k2pnn, eta2_k2pnn, k2pnnsel[4], DrawOpts);
	  drawfunc(rho1_k2pnn, eta1_k2pnn, k2pnnsel[3], DrawOpts);
	} else 
	  drawfunc(rho2_k2pnn, eta2_k2pnn, k2pnnsel[4], DrawOpts);
      else if (k2pnnsel[3] > 0)
	drawfunc(rho1_k2pnn, eta1_k2pnn, k2pnnsel[3], DrawOpts);
    } else if (DrawOpts.compare("CONT")==0) {
      if (k2pnnsel[0] > 0)
	drawfunc(rho1_k2pnn, eta1_k2pnn, k2pnnsel[0], DrawOpts);
      if (k2pnnsel[1] > 0)
 	drawfunc(rho2_k2pnn, eta2_k2pnn, k2pnnsel[1], DrawOpts);
    }
  }

  // Draw klp0nn
  if (kl2p0nnd == 1)
    setdrawlevel(kl2p0nnConstrain,kl2p0nnlevel,kl2p0nnsel, DrawOpts);
  if (kl2p0nnd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (kl2p0nnsel[4] > 0)
	if (kl2p0nnsel[3] > 0) {
	  kl2p0nn_areas(kl2p0nns2mn,kl2p0nns2mx,kl2p0nnsel[4]);
	  kl2p0nn_areas(kl2p0nns1mn,kl2p0nns1mx,kl2p0nnsel[3]);
	  kl2p0nn_areas(-kl2p0nns2mn,-kl2p0nns2mx,kl2p0nnsel[4]);
	  kl2p0nn_areas(-kl2p0nns1mn,-kl2p0nns1mx,kl2p0nnsel[3]);
	} else {
	  kl2p0nn_areas(kl2p0nns2mn,kl2p0nns2mx,kl2p0nnsel[4]);
	  kl2p0nn_areas(-kl2p0nns2mn,-kl2p0nns2mx,kl2p0nnsel[4]);
	} else if (kl2p0nnsel[3] > 0) {
	  kl2p0nn_areas(kl2p0nns1mn,kl2p0nns1mx,kl2p0nnsel[3]);
	  kl2p0nn_areas(-kl2p0nns1mn,-kl2p0nns1mx,kl2p0nnsel[3]);
	}
    } else if (DrawOpts.compare("CONT")==0) {
      if (kl2p0nnsel[0] > 0) {
 	kl2p0nn_lines(kl2p0nns1mn,kl2p0nns1mx,kl2p0nnsel[0]);
 	kl2p0nn_lines(-kl2p0nns1mn,-kl2p0nns1mx,kl2p0nnsel[0]);
      }
      if (kl2p0nnsel[1] > 0) {
 	kl2p0nn_lines(kl2p0nns2mn,kl2p0nns2mx,kl2p0nnsel[1]);
 	kl2p0nn_lines(-kl2p0nns2mn,-kl2p0nns2mx,kl2p0nnsel[1]);
      }
    }
  }
  
  // Draw BtoKpp
  if (BtoKppd == 1) 
    setdrawlevel(BtoKppConstrain,BtoKpplevel,BtoKppsel,DrawOpts);
  if (BtoKppd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (BtoKppsel[4] > 0)
	if (BtoKppsel[3] > 0) {
	  drawBtoKpp(tan(BtoKpps2mn),tan(BtoKpps2mx),BtoKppsel[4]);
	  drawBtoKpp(tan(BtoKpps1mn),tan(BtoKpps1mx),BtoKppsel[3]);
	} else {
	  drawBtoKpp(tan(BtoKpps2mn),tan(BtoKpps2mx),BtoKppsel[4]);
	}
      else if (BtoKppsel[3] > 0) {
 	drawBtoKpp(tan(BtoKpps1mn),tan(BtoKpps1mx),BtoKppsel[3]);
      }
    } else if (DrawOpts.compare("CONT")==0) {
      if (BtoKppsel[0] > 0) {
 	BtoKpp_lines(tan(BtoKpps1mn),tan(BtoKpps1mx),BtoKppsel[0]);
      }
      if (BtoKppsel[1] > 0) {
 	BtoKpp_lines(tan(BtoKpps2mn),tan(BtoKpps2mx),BtoKppsel[1]);
      }
    }
  }

  // Draw btovg
  if (btovgd == 1) 
    setdrawlevel(btovgConstrain,btovglevel,btovgsel,DrawOpts);
  if (btovgd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (btovgsel[4] > 0)
	if (btovgsel[3] > 0) {
	  drawfunc(rho2_btovg, eta2_btovg, btovgsel[4], DrawOpts);
	  drawfunc(rho1_btovg, eta1_btovg, btovgsel[3], DrawOpts);
	} else
	  drawfunc(rho2_btovg, eta2_btovg, btovgsel[4], DrawOpts);
      else if (btovgsel[3] > 0)
	drawfunc(rho1_btovg, eta1_btovg, btovgsel[3], DrawOpts);
    } else if (DrawOpts.compare("CONT")==0) {
      if (btovgsel[0] > 0)
 	drawfunc(rho1_btovg, eta1_btovg, btovgsel[0], DrawOpts);
      if (btovgsel[1] > 0)
 	drawfunc(rho2_btovg, eta2_btovg, btovgsel[1], DrawOpts);
    }
  }

  // Draw dmd
  if (dmdd == 1) 
    setdrawlevel(dmdConstrain,dmdlevel,dmdsel,DrawOpts);
  if (dmdd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (dmdsel[4] > 0)
	if (dmdsel[3] > 0) {
	  drawfunc(rho2_dmd, eta2_dmd, dmdsel[4], DrawOpts);
	  drawfunc(rho1_dmd, eta1_dmd, dmdsel[3], DrawOpts);
	} else
	  drawfunc(rho2_dmd, eta2_dmd, dmdsel[4], DrawOpts);
      else if (dmdsel[3] > 0)
	drawfunc(rho1_dmd, eta1_dmd, dmdsel[3], DrawOpts);
    } else if (DrawOpts.compare("CONT")==0) {
      if (dmdsel[0] > 0)
 	drawfunc(rho1_dmd, eta1_dmd, dmdsel[0], DrawOpts);
      if (dmdsel[1] > 0)
 	drawfunc(rho2_dmd, eta2_dmd, dmdsel[1], DrawOpts);
    }
  }

  // Draw dms
  if (dmsd == 1) { 
    dmsConstrain->SetLineStyle(3);
    setdrawlevel(dmsConstrain,dmslevel,dmssel,DrawOpts);
  }
  if ( dmsd==2 ) {
    if(doMFV != 0) {
      if (DrawOpts.compare("AREA")==0) {
	if (dmssel[4] > 0)
	  if (dmssel[3] > 0) {
	    drawfunc(rho2_dms, eta2_dms, dmssel[4], DrawOpts);
	    drawfunc(rho1_dms, eta1_dms, dmssel[3], DrawOpts);
	  } else
	    drawfunc(rho2_dms, eta2_dms, dmssel[4], DrawOpts);
	else if (dmssel[3] > 0)
	  drawfunc(rho1_dms, eta1_dms, dmssel[3], DrawOpts);
      } else if (DrawOpts.compare("CONT")==0) {
	if (dmssel[0] > 0)
	  drawfunc(rho1_dms, eta1_dms, dmssel[0], DrawOpts);
	if (dmssel[1] > 0)
	  drawfunc(rho2_dms, eta2_dms, dmssel[1], DrawOpts);
      }
    } else {
      if (DrawOpts.compare("CONT")==0) {
	if (dmssel[1] > 0)
	  drawfunc(rho2_dms, eta2_dms, dmssel[1], "DASHED");
      } else if (DrawOpts.compare("AREA")==0) {
	if (dmssel[4] > 0)
	  drawfunc(rho2_dms, eta2_dms, dmssel[4], DrawOpts);
      }
    }
  }

  // Draw b -> tau nu
  if (b2taunud == 1) 
    setdrawlevel(b2taunuConstrain,b2taunulevel,b2taunusel,DrawOpts);
  if (b2taunud == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (b2taunusel[4] > 0)
	if (b2taunusel[3] > 0) {
	  drawfunc(rho2_b2taunu, eta2_b2taunu, b2taunusel[4], DrawOpts);
	  drawfunc(rho1_b2taunu, eta1_b2taunu, b2taunusel[3], DrawOpts);
	} else
	  drawfunc(rho2_b2taunu, eta2_b2taunu, b2taunusel[4], DrawOpts);
      else if (b2taunusel[3] > 0)
 	drawfunc(rho1_b2taunu, eta1_b2taunu, b2taunusel[3], DrawOpts);
    } else if (DrawOpts.compare("CONT")==0) {
      if (b2taunusel[0] > 0)
 	drawfunc(rho1_b2taunu, eta1_b2taunu, b2taunusel[0], DrawOpts);
      if (b2taunusel[1] > 0)
 	drawfunc(rho2_b2taunu, eta2_b2taunu, b2taunusel[1], DrawOpts);
    }
  }

  // Draw Vub
  if (vubd == 1) 
    setdrawlevel(vubConstrain,vublevel,vubsel,DrawOpts);
  if (vubd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (vubsel[4] > 0)
	if (vubsel[3] > 0) {
	  drawfunc(rho2_vub, eta2_vub, vubsel[4], DrawOpts);
	  drawfunc(rho1_vub, eta1_vub, vubsel[3], DrawOpts);
	} else
	  drawfunc(rho2_vub, eta2_vub, vubsel[4], DrawOpts);
      else if (vubsel[3] > 0)
 	drawfunc(rho1_vub, eta1_vub, vubsel[3], DrawOpts);
    } else if (DrawOpts.compare("CONT")==0) {
      if (vubsel[0] > 0)
 	drawfunc(rho1_vub, eta1_vub, vubsel[0], DrawOpts);
      if (vubsel[1] > 0)
 	drawfunc(rho2_vub, eta2_vub, vubsel[1], DrawOpts);
    }
  }

  // Draw epsk
  if (epsd == 1) 
    setdrawlevel(epsConstrain,epslevel,epssel,DrawOpts);
  if ( epsd==2 ) {
    if (DrawOpts.compare("AREA")==0) {
      if (epssel[4] > 0) {
	if (epssel[3] > 0) {
	  drawfunc200(rho2_eps, eta2_eps, epssel[4], DrawOpts);
	  drawfunc200(rho1_eps, eta1_eps, epssel[3], DrawOpts);
	  drawfunc200(rho2m_eps, eta2m_eps, epssel[4], DrawOpts);
	  drawfunc200(rho1m_eps, eta1m_eps, epssel[3], DrawOpts);
	} else {
	  drawfunc200(rho2_eps, eta2_eps, epssel[4], DrawOpts);
	  drawfunc200(rho2m_eps, eta2m_eps, epssel[4], DrawOpts);
	} 
      } else if (epssel[3] > 0) {
	drawfunc200(rho1_eps, eta1_eps, epssel[3], DrawOpts);
	drawfunc200(rho1m_eps, eta1m_eps, epssel[3], DrawOpts);
      }
    } else if (DrawOpts.compare("CONT")==0) {
      if (epssel[0] > 0) {
	drawfunc200(rho1_eps, eta1_eps, epssel[0], DrawOpts);
	drawfunc200(rho1m_eps, eta1m_eps, epssel[0], DrawOpts);
      }
      if (epssel[1] > 0) {
	drawfunc200(rho2_eps, eta2_eps, epssel[1], DrawOpts);
	drawfunc200(rho2m_eps, eta2m_eps, epssel[1], DrawOpts);
      }
    }
  }

  // Draw s2b
  if (s2bd == 1) 
    setdrawlevel(s2bConstrain,s2blevel,s2bsel,DrawOpts);
  if (s2bd == 2) {
    if (DrawOpts.compare("AREA")==0) {
      if (s2bsel[4] > 0)
	if (s2bsel[3] > 0) {
	  //	  sin2b_areas(b1s2mn,b1s2mx,s2bsel[4]);
	  sin2b_areas(b2s2mn,b2s2mx,s2bsel[4]);
	  //	  sin2b_areas(b1s1mn,b1s1mx,s2bsel[3]);
	  sin2b_areas(b2s1mn,b2s1mx,s2bsel[3]);
	} else {
	  //	  sin2b_areas(b1s2mn,b1s2mx,s2bsel[4]);
	  sin2b_areas(b2s2mn,b2s2mx,s2bsel[4]);
	}
      else if (s2bsel[3] > 0) {
	//	sin2b_areas(b1s1mn,b1s1mx,s2bsel[3]);
	sin2b_areas(b2s1mn,b2s1mx,s2bsel[3]);
      }
    } else if (DrawOpts.compare("CONT")==0) {
      if (s2bsel[0] > 0) {
	//	sin2b_lines(b1s1mn,b1s1mx,s2bsel[0]);
	sin2b_lines(b2s1mn,b2s1mx,s2bsel[0]);      
      }
      if (s2bsel[1] > 0) {
	//	sin2b_lines(b1s2mn,b1s2mx,s2bsel[1]);
	sin2b_lines(b2s2mn,b2s2mx,s2bsel[1]);
      }
    }
  }

}

void BtoKpp_lines(double bound1, double bound2,int color) {    
    
  double xx[3],yy[3];
  xx[0] = rho0CPS;   yy[0] = 0.;
  xx[1] = -2.0; yy[1] = bound1*(xx[1]-rho0CPS);
  TGraph* BtoKpp1 = new TGraph(2,xx,yy);
  BtoKpp1->SetLineColor(color);
  BtoKpp1->Draw();
  xx[1] = -2.0; yy[1] = bound2*(xx[1]-rho0CPS);
  BtoKpp1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = bound1*(xx[1]-rho0CPS);
  BtoKpp1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = bound2*(xx[1]-rho0CPS);
  BtoKpp1->DrawGraph(2,xx,yy);
}		       

void BtoKpp_areas(double bound1, double bound2,int color) {
     
  double xx[4],yy[4];
  xx[0] = rho0CPS;   yy[0] = 0.;
  yy[1] = -1.25; xx[1] = (yy[1]/bound1)+rho0CPS;
  yy[2] = -1.25; xx[2] = -1.25;
  yy[3] = -1.25; xx[3] = (yy[3]/bound2)+rho0CPS;

  TGraph* BtoKpp1 = new TGraph(4,xx,yy);
  BtoKpp1->SetFillColor(color);
  BtoKpp1->Draw("F");
     
  yy[1] =  1.25; xx[1] = (yy[1]/bound1)+rho0CPS;
  yy[2] =  1.25; xx[2] = 1.25;
  yy[3] =  1.25; xx[3] = (yy[3]/bound2)+rho0CPS;
  BtoKpp1->DrawGraph(4,xx,yy,"F");
 
}
void sin2b_lines(double bound1, double bound2,int color) {    
    
  double xx[3],yy[3];
  xx[0] = 1.;   yy[0] = 0.;
  xx[1] = -2.; yy[1] = 3.0/tan(bound1);
  TGraph* sen2b1 = new TGraph(2,xx,yy);
  sen2b1->SetLineColor(color);
  sen2b1->Draw();
  xx[1] = -2.0; yy[1] = 3.0/tan(bound2);
  sen2b1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = -3.0/tan(bound1);
  sen2b1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = -3.0/tan(bound2);
  sen2b1->DrawGraph(2,xx,yy);
}		       

void sin2b_areas(double bound1, double bound2,int color) {
     
  double xx[3],yy[3];
  xx[0] = 1.;   yy[0] = 0.;
  xx[1] = -2.0; yy[1] = 3.0/tan(bound1);
  xx[2] = -2.0; yy[2] = 3.0/tan(bound2);

  TGraph* sen2b1 = new TGraph(3,xx,yy);
  sen2b1->SetFillColor(color);
  sen2b1->Draw("F");
     
  xx[1] =  4.; yy[1] = -3.0/tan(bound1);
  xx[2] =  4.; yy[2] = -3.0/tan(bound2);
  sen2b1->DrawGraph(3,xx,yy,"F");
 
}

void gamma_lines(double bound1, double bound2,int color) {    
    
  double xx[3],yy[3];
  xx[0] = 0.;   yy[0] = 0.;
  xx[1] = -2.; yy[1] = -2.0*tan(bound1);
  TGraph* gamma1 = new TGraph(2,xx,yy);
  gamma1->SetLineColor(color);
  gamma1->Draw();
  xx[1] = -2.0; yy[1] = -2.0*tan(bound2);
  gamma1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = 4.0*tan(bound1);
  gamma1->DrawGraph(2,xx,yy);
  xx[1] =  4.; yy[1] = 4.0*tan(bound2);
  gamma1->DrawGraph(2,xx,yy);
}		       


void gamma_areas(double bound1, double bound2,int color) {
     
  double xx[3],yy[3];
  xx[0] = 0.;   yy[0] = 0.;
  xx[1] = -2.0/tan(bound1); yy[1] = -2.0;
  xx[2] = -2.0/tan(bound2); yy[2] = -2.0;

  TGraph* gamma1 = new TGraph(3,xx,yy);
  gamma1->SetFillColor(color);
  gamma1->Draw("F");
     
  yy[0] =  4.; xx[0] = 4.0/tan(bound2);
  yy[1] =  4.; xx[1] = 4.0/tan(bound1);
  yy[2] = 0.;  xx[2] = 0.;
    gamma1->DrawGraph(3,xx,yy,"F");
}

void cos2b_lines(double bound,int color) {    
    
  double xx[3],yy[3];
  double twob = acos(-1.)-acos(bound);
  xx[0] = -2.; yy[0] = 3.0/tan(twob/2);
  xx[1] =  1.; yy[1] = 0.;
  xx[2] =  2.; yy[2] = 1.0/tan(twob/2);
  TGraph* sen2b1 = new TGraph(3,xx,yy);
  sen2b1->SetLineColor(color);
  sen2b1->Draw();
  xx[0] = -2.; yy[0] = -3.0/tan(twob/2);
  xx[1] =  1.; yy[1] = 0.;
  xx[2] =  2.; yy[2] = -1.0/tan(twob/2);
  sen2b1->DrawGraph(3,xx,yy);  
}		       

void draw_cos2b(double bound, int color) {
  // left side
  double twob = acos(bound);
  if(tan(twob/2.) < 1.25/2.25) {
    if(tan(twob/2.) < 0.25/2.25) {      
      cos2b_areas(twob, twob, color, "LEFT");
    } else {
      cos2b_areas(twob, 2*atan2(0.25,2.25), color, "LEFT");
      cos2b_areas(2*atan2(0.25,2.25), twob, color, "LEFT");
    }
  } else {
    if(tan(twob/2.) < 0.25/2.25) {    
      cos2b_areas(2*atan2(1.25,2.25), twob, color, "LEFT");
      cos2b_areas(twob,2*atan2(1.25,2.25), color, "LEFT");
    } else {
      cos2b_areas(twob, 2*atan2(0.25,2.25), color, "LEFT");
      cos2b_areas(2*atan2(0.25,2.25), 2*atan2(1.25,2.25), color, "LEFT");
      cos2b_areas(2*atan2(1.25,2.25), twob,color, "LEFT");
    }
  }

  // Right side
  if(tan(twob/2.) < 1.25/0.25) {
    if(tan(twob/2.) < 0.25/0.25) {
      cout << "sto al 1" << endl;
      cos2b_areas(twob, twob, color, "RIGHT");
    } else {
      cout << "sto al 1" << endl;
      cos2b_areas(acos(-1.)/2., twob, color, "RIGHT");
      cos2b_areas(twob,acos(-1.)/2., color, "RIGHT");
    }
  } else {
    if(tan(twob/2.) < 0.25/0.25) {    
      cout << "sto al 3" << endl;
      cos2b_areas(2*atan2(1.25,0.25), twob, color, "RIGHT");
      cos2b_areas(twob,2*atan2(1.25,0.25), color, "RIGHT");
    } else {
      cout << "sto al 4" << endl;
      cos2b_areas(twob, acos(-1.)/2., color, "RIGHT");
      cos2b_areas(2*atan2(1.25,0.25), acos(-1.)/2., color, "RIGHT");
      cos2b_areas(twob,2*atan2(1.25,0.25), color, "RIGHT");
    }
  }
}

void cos2b_areas(double twob_min, double twob_max, int color, string sides) {    
  double xx[3],yy[3];
  if(sides.compare("RIGHT")==0) {
    xx[0] =  1.25; yy[0] = 0.25*tan(twob_max/2);
    xx[1] =  1.; yy[1] = 0;
    xx[2] =  1.25; yy[2] = -0.25*tan(twob_min/2);
  } else {
    xx[0] = -1.25;   yy[0] = 2.25*tan(twob_max/2);
    xx[1] =  1.;   yy[1] = 0;
    xx[2] = -1.25;   yy[2] = -2.25*tan(twob_min/2);
  }
  TGraph* sen2b1 = new TGraph(3,xx,yy);
  sen2b1->SetFillColor(color);
  sen2b1->Draw("F");

}		       

void D0p0_lines(double beta_min, double beta_max ,int color) {    
    
  double xx[3],yy[3];
//   double twob = acos(-1.)-acos(bound);
//   xx[0] = -2.; yy[0] = 3.0/tan(twob/2);
//   xx[1] =  1.; yy[1] = 0.;
//   xx[2] =  2.; yy[2] = 1.0/tan(twob/2);
 
  xx[0] =  1.25; yy[0] = 0.25*tan(beta_min);
  xx[1] =  1.; yy[1] = 0;
  xx[2] =  1.25; yy[2] = -0.25*tan(beta_max);

  TGraph* sen2b1 = new TGraph(3,xx,yy);
  sen2b1->SetLineColor(color);
  sen2b1->Draw();
//   xx[0] = -2.; yy[0] = -3.0/tan(twob/2);
//   xx[1] =  1.; yy[1] = 0.;
//   xx[2] =  2.; yy[2] = -1.0/tan(twob/2);
    xx[0] = -1.25;   yy[0] = 2.25*tan(beta_max);
    xx[1] =  1.;   yy[1] = 0;
    xx[2] = -1.25;   yy[2] = -2.25*tan(beta_min);
  sen2b1->DrawGraph(3,xx,yy);  
}		       

void D0p0_areas(double beta_min, double beta_max, int color, string sides) {    
    
  double xx[5],yy[5];
  if(sides.compare("RIGHT")==0) {
    xx[0] =  1.25; yy[0] = 0.25*tan(beta_min);
    xx[1] =  1.; yy[1] = 0;
    xx[2] =  1.25; yy[2] = -0.25*tan(beta_max);
    if(yy[2]<-0.25) {
      xx[3] =  1.25; yy[3] = -0.25;
    } else {
      xx[3] =  1.25; yy[3] = -0.25*tan(beta_max);
    }
    if(yy[0]>1.25) {
      xx[4] =  1.25; yy[4] = 1.25;
    } else {
      xx[4] =  1.25; yy[4] = 0.25*tan(beta_min);
    }
  } else {
    xx[0] = -1.25;   yy[0] = 2.25*tan(beta_max);
    xx[1] =  1.;   yy[1] = 0;
    xx[2] = -1.25;   yy[2] = -2.25*tan(beta_min);
    if(yy[2]<-0.25) {
      xx[3] =  -1.25; yy[3] = -0.25;
    } else {
      xx[3] =  -1.25; yy[3] = -0.25*tan(beta_min);
    }
    if(yy[0]>1.25) {
      xx[4] =  -1.25; yy[4] = 1.25;
    } else {
      xx[4] =  -1.25; yy[4] = 0.25*tan(beta_max);
    }
  }
  TGraph* D0p0b1 = new TGraph(3,xx,yy);
  D0p0b1->SetFillColor(color);
  D0p0b1->Draw("F");
}		       

void kl2p0nn_lines(double kl2p0nnsmn, double kl2p0nnsmx, int color) {
    
  double xx[4],yy[4];
  xx[0] = -1.25; yy[0] = kl2p0nnsmn;
  xx[1] =  1.25; yy[1] = kl2p0nnsmn;
  xx[2] =  1.25; yy[2] = kl2p0nnsmx;
  xx[3] = -1.25; yy[3] = kl2p0nnsmx;
  
  TGraph* kl2p0nn1 = new TGraph(4,xx,yy);
  kl2p0nn1->SetLineColor(color);
  kl2p0nn1->Draw();
}		       

void kl2p0nn_areas(double kl2p0nnsmn, double kl2p0nnsmx, int color) {
     
  double xx[4],yy[4];
  xx[0] = -1.25; yy[0] = kl2p0nnsmn;
  xx[1] =  1.25; yy[1] = kl2p0nnsmn;
  xx[2] =  1.25; yy[2] = kl2p0nnsmx;
  xx[3] = -1.25; yy[3] = kl2p0nnsmx;

  TGraph* kl2p0nn1 = new TGraph(4,xx,yy);
  kl2p0nn1->SetFillColor(color);
  kl2p0nn1->Draw("F");
     
}

void drawgamma(double min, double max, int col) {
  if (min*max > 0) {
    if(tan(min)>0)
      if(tan(min)<0.25/1.25) {
	gamma_areas(min,atan2(0.25,1.25),col);
	if(tan(min)<1.) {
	  gamma_areas(atan2(0.25,1.25),acos(-1.)/4.,col);
	  if(tan(max)<-1)
	    gamma_areas(acos(-1.)/4.,max,col);
	  else {
	    gamma_areas(acos(-1.)/4,3*acos(-1.)/4,col);
	    gamma_areas(3*acos(-1.)/4.,max,col);
	  }
	}
      } else if(tan(min)<1.) {
	gamma_areas(min,acos(-1.)/4.,col);
	gamma_areas(acos(-1.)/4.,max,col);
      } else 
	gamma_areas(min,max,col);
    else 
      gamma_areas(min,max,col);
  }
  else 
    if (tan(min)<0) {
      if(tan(min)<-1) {
	gamma_areas(min,3*acos(-1.)/4,col);
	gamma_areas(3*acos(-1.)/4,atan2(-0.25,1.25),col);
	gamma_areas(atan2(-0.25,1.25),atan2(0.,-1.),col);
      } else if(tan(min)<-0.25/1.25) {
	gamma_areas(min,atan2(-0.25,1.25),col);
	gamma_areas(atan2(-0.25,1.25),atan2(0.,-1.),col);
      } else 
	gamma_areas(min,0,col);
      if(tan(max)<0.25/1.25) {
	gamma_areas(atan2(0.,1.),max,col);
      } else if(tan(max)<1.) {
	gamma_areas(acos(0.9999999),atan2(0.25,1.25),col);
	gamma_areas(atan2(0.25,1.25),max,col);	
      } else {
	gamma_areas(acos(0.9999999),atan2(0.25,1.25),col);
	gamma_areas(atan2(0.25,1.25),acos(-1.)/4.,col);
	gamma_areas(acos(-1.)/4.,max,col);
      }
    }
}

void drawBtoKpp(double min, double max, int col) {
  if(max>0.25/(1.25-rho0CPS)) {
    if(min<-0.25/(1.25-rho0CPS)) {
      BtoKpp_areas(min,-0.25/(1.25-rho0CPS),col);
      min = -0.25/(1.25-rho0CPS);
    }
    if(max>0.25/(1.25+rho0CPS)) {
      if(min < 0.25/(1.25+rho0CPS)) {
	BtoKpp_areas(min,0.25/(1.25+rho0CPS),col);
	min = 0.25/(1.25+rho0CPS);
      }
      if(max>1.25/(1.25+rho0CPS)) {
	if(min < 1.25/(1.25+rho0CPS)) {
	  BtoKpp_areas(min,1.25/(1.25+rho0CPS),col);	  
	  min = 1.25/(1.25+rho0CPS);
	}
	if(max>1.25/(1.25-rho0CPS)) {
	  if(min < 1.25/(1.25-rho0CPS)) {
	    BtoKpp_areas(min,1.25/(1.25-rho0CPS),col);	  
	    min = 1.25/(1.25-rho0CPS);
	  }
	}
      }
    }
  }
  BtoKpp_areas(min,max,col);
}

void set_data(string dir) {  
//   Vub
  if(set_vub[0] == 2) {

    TH1D* tmp_vub = (TH1D*) datafile->Get("Input/input_Rb");
    double vubOvcb = tmp_vub->GetMean();
    double err_vubOvcb = tmp_vub->GetRMS();

    vub1mn=vubOvcb-err_vubOvcb;
    vub1mx=vubOvcb+err_vubOvcb;
    vub2mn=vubOvcb-2*err_vubOvcb;
    vub2mx=vubOvcb+2*err_vubOvcb;

  }

  //  BtoVg
  if(set_btovg[0] == 2) {
    TString btovg_tit = "Input/input_Rt_all";
    if(btovg_use_his == 2) btovg_tit = "Input/input_Rt_neutral";
    if(btovg_use_his == 3) btovg_tit = "Input/input_Rt_charged";
    TH1D* tmp_btovg = (TH1D*) datafile->Get(btovg_tit)->Clone();

    btovg1mn=CalcLevel_Mediana(tmp_btovg,0.68,"Low");
    btovg1mx=CalcLevel_Mediana(tmp_btovg,0.68,"Up");
    btovg2mn=CalcLevel_Mediana(tmp_btovg,0.95,"Low");
    btovg2mx=CalcLevel_Mediana(tmp_btovg,0.95,"Up");
    
    delete tmp_btovg;    
  }

//   B2taunu
  if(set_b2taunu[0] == 2) {
    TH1D* tmp_b2taunu = (TH1D*) datafile->Get("Input/input_Rb_btaunu")->Clone();

    b2taunu1mn=CalcLevel(tmp_b2taunu,0.68,"Low");
    b2taunu1mx=CalcLevel(tmp_b2taunu,0.68,"Up");
    b2taunu2mn=CalcLevel(tmp_b2taunu,0.95,"Low");
    b2taunu2mx=CalcLevel(tmp_b2taunu,0.95,"Up");

    delete tmp_b2taunu;
  }

//  delta_md  
  if(set_dmd[0] == 2) {
    TH1D* tmp_dmd = (TH1D*) datafile->Get("Input/input_dmd")->Clone();

    dmd=tmp_dmd->GetMean();
    err_dmd =tmp_dmd->GetRMS();
    delete tmp_dmd;
    
    dmd1mn=dmd-err_dmd;
    dmd1mx=dmd+err_dmd;
    dmd2mn=dmd-2*err_dmd;
    dmd2mx=dmd+2*err_dmd;
  }

//  epsilon  
  if(set_eps[0] == 2) {
    
    TH1D* tmp_y0 = (TH1D*) (datafile->Get("Input/input_epsk_y0"))->Clone();
    TH1D* tmp_x0 = (TH1D*) (datafile->Get("Input/input_epsk_x0"))->Clone();
/*      
    epsk2mn[0] = CalcLevel_Mediana(tmp_y0,0.95,"Low");
    epsk2mn[1] = CalcLevel_Mediana(tmp_x0,0.95,"Up");
    epsk2mn[2] = 2.8;
    
    epsk1mn[0] = CalcLevel_Mediana(tmp_y0,0.68,"Low");
    epsk1mn[1] = CalcLevel_Mediana(tmp_x0,0.68,"Up");
    epsk1mn[2] = 2.8;
    
    epsk2mx[0] = CalcLevel_Mediana(tmp_y0,0.95,"Up");
    epsk2mx[1] = CalcLevel_Mediana(tmp_x0,0.95,"Low");
    epsk2mx[2] = 2.8;
    
    epsk1mx[0] = CalcLevel_Mediana(tmp_y0,0.68,"Up");
    epsk1mx[1] = CalcLevel_Mediana(tmp_x0,0.68,"Low");
    epsk1mx[2] = 2.8;
 */
    epsk2mn[0] = tmp_y0->GetMean()-2.*tmp_y0->GetRMS();
    epsk2mn[1] = tmp_x0->GetMean()-2.*tmp_x0->GetRMS();
    epsk2mn[2] = 2.8;
    
    epsk1mn[0] = tmp_y0->GetMean()-tmp_y0->GetRMS();
    epsk1mn[1] = tmp_x0->GetMean()-tmp_x0->GetRMS();
    epsk1mn[2] = 2.8;
    
    epsk2mx[0] = tmp_y0->GetMean()+2.*tmp_y0->GetRMS();
    epsk2mx[1] = tmp_x0->GetMean()+2.*tmp_x0->GetRMS();
    epsk2mx[2] = 2.8;
    
    epsk1mx[0] = tmp_y0->GetMean()+tmp_y0->GetRMS();
    epsk1mx[1] = tmp_x0->GetMean()+tmp_x0->GetRMS();
    epsk1mx[2] = 2.8;
   
    delete tmp_y0;
    delete tmp_x0;
  }
  
  //  B to Kpp
  if(set_BtoKpp[0] == 2) {

    TH1D* tmp_BtoKpp_rho0  = (TH1D*) datafile->Get("Input/input_rho0_CPS")->Clone();
    TH1D* tmp_BtoKpp_coeff = (TH1D*) datafile->Get("Input/input_coeffang_CPS")->Clone();
    
    rho0CPS = tmp_BtoKpp_rho0->GetMean();

    BtoKpps2mn = CalcLevel(tmp_BtoKpp_coeff,0.95,"Low");
    BtoKpps1mn = CalcLevel(tmp_BtoKpp_coeff,0.68,"Low");
    BtoKpps2mx = CalcLevel(tmp_BtoKpp_coeff,0.95,"Up");
    BtoKpps1mx = CalcLevel(tmp_BtoKpp_coeff,0.68,"Up");

    delete tmp_BtoKpp_rho0;
    delete tmp_BtoKpp_coeff;

  }

  //  sen2beta
  if(set_s2b[0] == 2) {
    //TH1D* tmp_sin2b = (TH1D*) datafile->Get("Input/input_sin2b")->Clone();
    TH1D* tmp_sin2b = (TH1D*) datafile->Get("Input/sin2b_tot")->Clone();
    sin2b = tmp_sin2b->GetMean();
    err_sin2b = tmp_sin2b->GetRMS();
    delete tmp_sin2b;
    
    b1s1mn=asin(sin2b-err_sin2b)/2.;
    b1s1mx=asin(sin2b+err_sin2b)/2.;
    b1s2mn=asin(sin2b-2*err_sin2b)/2.;
    b1s2mx=asin(sin2b+2*err_sin2b)/2.;
    b2s1mn=(acos(-1.)-asin(sin2b-err_sin2b))/2.;
    b2s1mx=(acos(-1.)-asin(sin2b+err_sin2b))/2.;
    b2s2mn=(acos(-1.)-asin(sin2b-2*err_sin2b))/2.;
    b2s2mx=(acos(-1.)-asin(sin2b+2*err_sin2b))/2.;
  }

  // gamma
  if(set_gamma[0] == 2) {
    //TString gamma_tit = "Input/input_gamma_deg_all";
    TString gamma_tit = "Input/input_gamma";
    if(gamma_use_his == 2) gamma_tit = "Input/input_gamma_deg_glw";
    if(gamma_use_his == 3) gamma_tit = "Input/input_gamma_deg_ads";
    if(gamma_use_his == 4) gamma_tit = "Input/input_gamma_deg_dal";

    TFile* gamma_root(0);
    TH1D* tmp_gamma(0);
    if(gamma_use_root != "") {
      gamma_root = new TFile(gamma_use_root.c_str());
      tmp_gamma = (TH1D*) gamma_root->Get(gamma_tit)->Clone();
    } else {
      tmp_gamma = (TH1D*) datafile->Get(gamma_tit)->Clone();
    }
    
    g1s1mn = Pi/180*CalcLevel(tmp_gamma,0.68,"Low");
    g1s1mx = Pi/180*CalcLevel(tmp_gamma,0.68,"Up");
    
    if(gamma_use_his == 2) {
      g1s1mn1 = Pi/180*CalcLevel(tmp_gamma,0.68,"Low2")+Pi;
      g1s1mx1 = Pi/180*CalcLevel(tmp_gamma,0.68,"Up2");
      g1s1mn2 = Pi/180*CalcLevel(tmp_gamma,0.68,"Low3")+Pi;
      g1s1mx2 = Pi/180*CalcLevel(tmp_gamma,0.68,"Up3");
    }

    g1s2mn = Pi/180*CalcLevel(tmp_gamma,0.95,"Low");
    g1s2mx = Pi/180*CalcLevel(tmp_gamma,0.95,"Up");

    g2s1mn = g1s1mn + Pi;
    g2s1mx = g1s1mx + Pi;
    if(gamma_use_his == 2) {
      g2s1mn1 = g1s1mn1 + Pi;
      g2s1mx1 = g1s1mx1 + Pi;
      g2s1mn2 = g1s1mn2 + Pi;
      g2s1mx2 = g1s1mx2 + Pi;
    }
    g2s2mn = g1s2mn + Pi;
    g2s2mx = g1s2mx + Pi;

  }

  // 2bpg
  if(set_s2bpg[0] == 2) {
    if(s2bpg_use_his < 3)
      s2bpg1D = (TH1D*) datafile->Get("Input/input_2bpg")->Clone();

    if(s2bpg_use_his == 3) 
      s2bpg1D = (TH1D*) datafile->Get("Input/input_2bpg_star")->Clone();
    if(s2bpg_use_his == 4)
      s2bpg1D = (TH1D*) datafile->Get("Input/input_2bpg_rho")->Clone();
    
    if(s2bpg_use_his == 1) {
      s2bpg1D->Multiply((TH1D*) datafile->Get("Input/input_2bpg_star")->Clone());
      s2bpg1D->Multiply((TH1D*) datafile->Get("Input/input_2bpg_rho")->Clone());
    }

    s2bpgs1= CalcThreshold(s2bpg1D,0.68);
    s2bpgs2= CalcThreshold(s2bpg1D,0.95);

  }

  // 2alpha
  if(set_s2a[0] == 2) {
    // usare  TString s2a_tit = "Input/input_alpha";
    //TString s2a_tit = "Input/all_input_alpha"; 
    TString s2a_tit = "Input/input_alpha";
    if(s2a_use_his == 2) s2a_tit = "Input/pipi_input_alpha";
    if(s2a_use_his == 3) s2a_tit = "Input/rhorho_input_alpha";
    if(s2a_use_his == 4) s2a_tit = "Input/rhopi_input_alpha";
    s2a1D = (TH1D*) datafile->Get(s2a_tit)->Clone();

    s2as1 = CalcThreshold(s2a1D,0.68);
    s2as2 = CalcThreshold(s2a1D,0.95);
  }
  
//  cos2beta
  if(set_c2b[0] == 2) {
    TH1D* tmp_cos2b = (TH1D*) datafile->Get("Input/input_cos2b")->Clone();
    b1s1_c=CalcLevel(tmp_cos2b,0.68,"Low");
    b1s2_c=CalcLevel(tmp_cos2b,0.95,"Low");
    delete tmp_cos2b;
  }


  //D0p0  
  if(set_D0p0[0] == 2) {
  // sono stanco, je le dico a mano
    cout << "WARNING:the input is in cos2b now!!!!" << endl;
  }

//   dm_s 
  if(set_dms[0] == 2) {
    TH1D* tmp_dms(0);
    if(doMFV==2) {
      tmp_dms= (TH1D*)(datafile->Get("Input/input_dmsodmd"))->Clone();
   } else
     tmp_dms= (TH1D*)(datafile->Get("Input/input_dms"))->Clone();
    dms=tmp_dms->GetMean();
    err_dms =tmp_dms->GetRMS();
    dms2mx=dms+2*err_dms;
    dms2mn=dms-2*err_dms;
    dms1mx=dms+err_dms;
    dms1mn=dms-err_dms;
  }

  // k p nn
  if(set_k2pnn[0] == 2) {
    
    TH1D* tmp_x0 =(TH1D*) datafile->Get("Input/input_k2pnn_BRfunc")->Clone();
    TH1D* tmp_x1 =(TH1D*) datafile->Get("Input/input_k2pnn_sigmasq")->Clone();
    TH1D* tmp_x2 =(TH1D*) datafile->Get("Input/input_k2pnn_rho0bar")->Clone();
    
    k2pnn2mn[0] = CalcLevel_Mediana(tmp_x0,0.95,"Low");
    k2pnn2mn[1] = CalcLevel_Mediana(tmp_x1,0.95,"Low");
    k2pnn2mn[2] = tmp_x2->GetMean(); // check that this is reasonable

    cout << k2pnn2mn[0] << "  "  << k2pnn2mn[1] << "  " << k2pnn2mn[2] << endl;

    k2pnn2mx[0] = CalcLevel_Mediana(tmp_x0,0.95,"Up");
    k2pnn2mx[1] = CalcLevel_Mediana(tmp_x1,0.95,"Up");
    k2pnn2mx[2] = k2pnn2mn[2];

    cout << k2pnn2mx[0] << "  "  << k2pnn2mx[1] << "  " << k2pnn2mx[2] << endl;

    k2pnn1mn[0] = CalcLevel_Mediana(tmp_x0,0.68,"Low");
    k2pnn1mn[1] = CalcLevel_Mediana(tmp_x1,0.68,"Low");
    k2pnn1mn[2] = k2pnn2mn[2];

    cout << k2pnn1mn[0] << "  "  << k2pnn1mn[1] << "  " << k2pnn1mn[2] << endl;

    k2pnn1mx[0] = CalcLevel_Mediana(tmp_x0,0.68,"Up");
    k2pnn1mx[1] = CalcLevel_Mediana(tmp_x1,0.68,"Up");
    k2pnn1mx[2] = k2pnn2mn[2];

    cout << k2pnn1mx[0] << "  "  << k2pnn1mx[1] << "  " << k2pnn1mx[2] << endl;
  }

  // kl -> p0nn
  if(set_kl2p0nn[0] == 2) {
    TH1D* tmp_kl2p0nn = (TH1D*)(((TDirectory*) 
    				 datafile->Get("Input"))->Get("input_kl2p0nn"))->Clone();
    
    kl2p0nns1mn=CalcLevel_Mediana(tmp_kl2p0nn,0.68,"Low");
    kl2p0nns1mx=CalcLevel_Mediana(tmp_kl2p0nn,0.68,"Up");
    kl2p0nns2mn=CalcLevel_Mediana(tmp_kl2p0nn,0.95,"Low");
    kl2p0nns2mx=CalcLevel_Mediana(tmp_kl2p0nn,0.95,"Up");

    delete tmp_kl2p0nn;    
  }

}

// DrawOpt  = 1  ----> Overimpose
// DrawOpt  = 0  ----> Draw in a new Canvas
// DrawOpt !=0,1 ----> Do not Draw

void get_data(const char* HistoName, const char* DirName,int DrawOpt,int filenum) {

  if ( filenum==2 ) 
    if(strcmp(DirName,"none"))
      Constrain = (TH2D*)(((TDirectory*) datafile2->Get(DirName))->Get(HistoName))->Clone();
    else
      Constrain = (TH2D*)(datafile2->Get(HistoName))->Clone();
  else 
    if(strcmp(DirName,"none"))
      Constrain = (TH2D*)(((TDirectory*) datafile->Get(DirName))->Get(HistoName))->Clone();
    else
      Constrain = (TH2D*)(datafile->Get(HistoName))->Clone();
  
  int nx = Constrain->GetNbinsX();
  int ny = Constrain->GetNbinsY();
  // delete all the underflow and overflow, otherwise fucked root 
  // take them into account in the normalization
  for(int ix = 0; ix <= nx+1; ix++) {
    Constrain->SetBinContent(ix,0,0.);
    Constrain->SetBinContent(ix,ny+1,0.);
  }
  for(int iy = 0; iy <= ny+1; iy++) {
    Constrain->SetBinContent(0,iy,0.);
    Constrain->SetBinContent(nx+1,iy,0.);
  }    
  Double_t Sum = Constrain->GetSum();
  Constrain->Scale(1./Sum);

  for (int ix=1; ix<=nx; ix++) {
    for (int iy=1; iy<=ny; iy++) {
      if ( Constrain->GetBinContent(ix,iy) < 0.00000001 ) {
	Constrain->SetBinContent(ix,iy,0.0000000);
      }
    }
  }

  Sum = Constrain->GetSum();
  Constrain->Scale(1./Sum);

  vector<Double_t> conf; 
  conf.push_back(0.9999);
  conf.push_back(0.999);
  conf.push_back(0.99);
  conf.push_back(0.95);
  conf.push_back(0.68);

  for (int j=0; j<5; j++) {
    level[j] = 0.;
    area[j] = 0.;
  }

  vector< Double_t > OrderedGArray;
  for (int ix=1; ix<=nx; ix++) {
    for (int iy=1; iy<=ny; iy++) {
      OrderedGArray.push_back(Constrain->GetBinContent(ix,iy));
    }
  }

  sort( OrderedGArray.begin(), OrderedGArray.end(),greater<double>() );

  vector< Double_t> SumArray;
  SumArray.push_back(OrderedGArray[0]);
  for (int i=1; i<= (int) OrderedGArray.size(); i++) {
    SumArray.push_back(SumArray[i-1] + OrderedGArray[i]);
    for (int j=0; j<5; j++) {
      if (SumArray[i] <= conf[j]) {
	level[j] = OrderedGArray[i];
	area[j] = SumArray[i];
      }
    }
  }

//  if(DrawOpt == 1 && positiveeta == 1) {
//    for (int ix=1; ix<=nx; ix++) {
//      for (int iy=1; iy<=32; iy++) {
//	Constrain->SetBinContent(ix,iy,Constrain->GetBinContent(ix,65-iy));	
//      }
//    }
//  }

  if(DrawOpt == 1 && positiveeta == 1) {
    for (int ix=1; ix<=nx; ix++) {
      for (int iy=1; iy<=ny; iy++) {
	if(2*iy<ny) 
	  Constrain->SetBinContent(ix,iy,Constrain->GetBinContent(ix,ny-iy));	
      }
    }
  }

  for (int k=0;k<smooth;k++)
    Constrain->Smooth();

  for (int j=0; j<5; j++) {
    cout << "level: " << level[j] << endl;
    cout << "area: " << area[j] << endl;
  }

  if (abs(DrawOpt) <= 1 ) {
    Constrain->SetContour(2);
    for (int k=3; k<5 ; k++) {
      int l=k-3;
      Constrain->SetContourLevel(l,level[k]);
    }
    if (DrawOpt == 1) {
      Constrain->SetLineColor(1);
      Constrain->SetLineWidth(2);
      Constrain->Draw("cont3same");
    }
    if (DrawOpt == 0) {
      Constrain->SetLineColor(1);
      Constrain->SetLineWidth(1);
      Constrain->Draw("cont3");
    }
  }
}

double Epsk(double rho, int coni) {
  double x1(0.),x2(0.4),y1,y2,A,B;
  if (coni == 1) {
    y1   = epsk1mx[0];
    y2   = epsk1mx[1];
  }
  if (coni == 2) {
    y1   = epsk2mx[0];
    y2   = epsk2mx[1];
  }
  if (coni == -1) {
    y1   = epsk1mn[0];
    y2   = epsk1mn[1];
  }
  if (coni == -2) {
    y1   = epsk2mn[0];
    y2   = epsk2mn[1];
  }
  
//  B = -(1.-y1/y2)/((1.-x2)-y1/y2*(1.-x1));
//  A = (y1-y2)/(1./(1+B*(1-x1))-1./(1+B*(1-x2)));
  B = (y1-y2)/((1-x2)*y2+(x1-1)*y1);
  A = (y1*(x1*y2-x2*y2))/((1-x2)*y2+(x1-1)*y1);
  return A*(1./(1+B*(1-rho)));
}

void prepare_epsk() {
  double lowx(-3.0),upx(1.3);
  double rho;
  for(int i=0; i<100; i++) {
    rho = ((double) i+0.5)/100.*(upx+fabs(lowx))+lowx;
    eta1_eps[i] = fabs(Epsk(rho,1));
    rho1_eps[i] = rho;
    eta1m_eps[i] = -fabs(Epsk(rho,1));
    rho1m_eps[i] = -rho+2;
  }
  for(int i=100; i<200; i++) {
    rho = upx-((double) i-99.5)/100.*(upx+fabs(lowx));
    eta1_eps[i] = fabs(Epsk(rho,-1));
    rho1_eps[i] = rho;
    eta1m_eps[i] = -fabs(Epsk(rho,-1));
    rho1m_eps[i] = -rho+2;
  }
  for(int i=0; i<100; i++) {
    rho = ((double) i+0.5)/100.*(upx+fabs(lowx))+lowx;
    eta2_eps[i] = fabs(Epsk(rho,2));
    rho2_eps[i] = rho;
    eta2m_eps[i] = -fabs(Epsk(rho,2));
    rho2m_eps[i] = -rho+2;
  }
  for(int i=100; i<200; i++) {
    rho = upx-((double) i-99.5)/100.*(upx+fabs(lowx));
    eta2_eps[i] = fabs(Epsk(rho,-2));
    rho2_eps[i] = rho;
    eta2m_eps[i] = -fabs(Epsk(rho,-2));
    rho2m_eps[i] = -rho+2;
  }  
}

double K2pnn(double eta, int coni, int sign) {
  double brfunc,sigmasq,rho0bar;
  
  if (coni == 1) {
    brfunc   = k2pnn1mx[0];
    sigmasq  = k2pnn1mx[1];
    rho0bar  = k2pnn1mx[2];
  }
  if (coni == 2) {
    brfunc    = k2pnn2mx[0];
    sigmasq   = k2pnn2mx[1];
    rho0bar   = k2pnn2mx[2];
  }
  if (coni == -1) {
    brfunc    = k2pnn1mn[0];
    sigmasq   = k2pnn1mn[1];
    rho0bar   = k2pnn1mn[2];
  }
  if (coni == -2) {
    brfunc    = k2pnn2mn[0];
    sigmasq   = k2pnn2mn[1];
    rho0bar   = k2pnn2mn[2];
  }


  if(brfunc >= sigmasq*pow(eta,2.)) {
    if(sign > 0) 
      return rho0bar + sqrt(brfunc - sigmasq*pow(eta,2.));
    else
      return rho0bar - sqrt(brfunc - sigmasq*pow(eta,2.));
  } else {
    return 1.25;
  }
}
  
void prepare_k2pnn() {
  double lowy,upy(1.25);
  if(doNP==1) 
    lowy=-1.25;
  else
    lowy=-0.25;
    
  double eta;
  for(int i=0; i<100; i++) {
    eta = ((double) i)/99.*(upy-lowy)+lowy;
    rho1_k2pnn[i] = K2pnn(eta,-1,-1);
    eta1_k2pnn[i] = eta;
    rho1_k2pnn[199-i] = K2pnn(eta,-1,1);
    eta1_k2pnn[199-i] = eta;
  }
  for(int i=0; i<100; i++) {
    eta = ((double) i)/99.*(upy-lowy)+lowy;
    rho1_k2pnn[200+i] = K2pnn(eta,1,1);
    eta1_k2pnn[200+i] = eta;
    rho1_k2pnn[399-i] = K2pnn(eta,1,-1);
    eta1_k2pnn[399-i] = eta;
  }

  for(int i=0; i<100; i++) {
    eta = ((double) i)/99.*(upy-lowy)+lowy;
    rho2_k2pnn[i] = K2pnn(eta,-2,-1);
    eta2_k2pnn[i] = eta;
    rho2_k2pnn[199-i] = K2pnn(eta,-2,1);
    eta2_k2pnn[199-i] = eta;
  }

  for(int i=0; i<100; i++) {
    eta = ((double) i)/99.*(upy-lowy)+lowy;
    rho2_k2pnn[200+i] = K2pnn(eta,2,1);
    eta2_k2pnn[200+i] = eta;
    rho2_k2pnn[399-i] = K2pnn(eta,2,-1);
    eta2_k2pnn[399-i] = eta;
  }  
}

////////////////////////////////////////////////////
//   Calculates circles for dmd, corresponding to //
//   dmd mean +/- 1sigma and 2sigma               //
////////////////////////////////////////////////////

void prepare_dmd() {
  double rho;
  // 1 sigma circles
  for(int i=0; i<100; i++) {
    // top quarter of circle, dmd + 1sigma
    rho = (1-dmd1mx)+i*(xmax-(1-dmd1mx))/99.;
    eta1_dmd[i]     =  sqrt(fabs(dmd1mx*dmd1mx-(rho-1)*(rho-1)));
    rho1_dmd[i]     =  rho;
    eta1_dmd[i+200] = -sqrt(fabs(dmd1mx*dmd1mx-(rho-1)*(rho-1)));
    rho1_dmd[i+200] =  rho;
    // top quarter of circle, dmd - 1sigma
    rho = xmax+i*(1-dmd1mn-xmax)/99.;
    eta1_dmd[i+100] =  sqrt(fabs(dmd1mn*dmd1mn-(rho-1)*(rho-1)));
    rho1_dmd[i+100] =  rho;
    eta1_dmd[i+300] = -sqrt(fabs(dmd1mn*dmd1mn-(rho-1)*(rho-1)));
    rho1_dmd[i+300] =  rho;
  }
  // 2 sigmas circles
  for(int i=0; i<100; i++) {
    // top quarter of circle, dmd + 2sigma
    rho = (1-dmd2mx)+i*(xmax-(1-dmd2mx))/99.;
    eta2_dmd[i]     =  sqrt(fabs(dmd2mx*dmd2mx-(rho-1)*(rho-1)));
    rho2_dmd[i]     =  rho;
    eta2_dmd[i+200] = -sqrt(fabs(dmd2mx*dmd2mx-(rho-1)*(rho-1)));
    rho2_dmd[i+200] =  rho;    
    // top quarter of circle, dmd - 2sigma
    rho = xmax+i*(1-dmd2mn-xmax)/99.;
    eta2_dmd[i+100] =  sqrt(fabs(dmd2mn*dmd2mn-(rho-1)*(rho-1)));
    rho2_dmd[i+100] =  rho;
    eta2_dmd[i+300] = -sqrt(fabs(dmd2mn*dmd2mn-(rho-1)*(rho-1)));
    rho2_dmd[i+300] =  rho;
  }
}

////////////////////////////////////////////////
//   Calculates circle for dms, corresponding //
//   to dms mean + 2sigma                     //
////////////////////////////////////////////////

void prepare_dms() {
  double rho;
  TH1D* tmp_lam =(TH1D*) datafile->Get("Input/lambda_2")->Clone();
  double lam2 = pow(tmp_lam->GetMean(),2.);
  delete tmp_lam;
  // 1 sigma circles
  double rho0mx1=1.+lam2*dms1mx*dms1mx-sqrt(dms1mx*dms1mx*(1.+lam2)+lam2*lam2*pow(dms1mx,4.));
  double rho0mn1=1.+lam2*dms1mn*dms1mn-sqrt(dms1mn*dms1mn*(1.+lam2)+lam2*lam2*pow(dms1mn,4.));
  for(int i=0; i<100; i++) {
    // top quarter of circle, dms + 1sigma
    rho = rho0mx1+i*(xmax-rho0mx1)/99.;
    eta1_dms[i]     =  sqrt(fabs(dms1mx*dms1mx*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho1_dms[i]     =  rho;
    eta1_dms[i+200] = -sqrt(fabs(dms1mx*dms1mx*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho1_dms[i+200] =  rho;
    // top quarter of circle, dms - 1sigma
    rho = xmax+i*(rho0mn1-xmax)/99.;
    eta1_dms[i+100] =  sqrt(fabs(dms1mn*dms1mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho1_dms[i+100] =  rho;
    eta1_dms[i+300] = -sqrt(fabs(dms1mn*dms1mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho1_dms[i+300] =  rho;
  }
  // 2 sigmas circles
  double rho0mx2=1.+lam2*dms2mx*dms2mx-sqrt(dms2mx*dms2mx*(1.+lam2)+lam2*lam2*pow(dms2mx,4.));
  double rho0mn2=1.+lam2*dms2mn*dms2mn-sqrt(dms2mn*dms2mn*(1.+lam2)+lam2*lam2*pow(dms2mn,4.));
  double rho0mn2d=1.+lam2*dmd2mn*dmd2mn-sqrt(dmd2mn*dmd2mn*(1.+lam2)+lam2*lam2*pow(dmd2mn,4.));
  for(int i=0; i<100; i++) {
    // top quarter of circle, dms + 2sigma
    rho = rho0mx2+i*(xmax-rho0mx2)/99.;
    eta2_dms[i]     =  sqrt(fabs(dms2mx*dms2mx*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho2_dms[i]     =  rho;
    eta2_dms[i+200] = -sqrt(fabs(dms2mx*dms2mx*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
    rho2_dms[i+200] =  rho;
    if(doMFV!=0) {
    // top quarter of circle, dms - 2sigma
      rho = xmax+i*(rho0mn2-xmax)/99.;
      eta2_dms[i+100] =  sqrt(fabs(dms2mn*dms2mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
      rho2_dms[i+100] =  rho;
      eta2_dms[i+300] = -sqrt(fabs(dms2mn*dms2mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
      rho2_dms[i+300] =  rho;
    } else {
    // top quarter of circle, dmd - 2sigma
      rho = xmax+i*(rho0mn2d-xmax)/99.;
      eta2_dms[i+100] =  sqrt(fabs(dmd2mn*dmd2mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
      rho2_dms[i+100] =  rho;
      eta2_dms[i+300] = -sqrt(fabs(dmd2mn*dmd2mn*(1.-lam2*(1.-2.*rho))-(rho-1.)*(rho-1.)));
      rho2_dms[i+300] =  rho;
    }
  }
}
//////////////////////////////////////////////////////
//   Calculates circles for btovg, corresponding to //
//   btovg mean +/- 1sigma and 2sigma               //
//////////////////////////////////////////////////////

void prepare_btovg() {
  double rho;
  // 1 sigma circles
  for(int i=0; i<100; i++) {
    // top quarter of circle, btovg + 1sigma
    rho = (1-btovg1mx)+i*(xmax-(1-btovg1mx))/99.;
    eta1_btovg[i]     =  sqrt(fabs(btovg1mx*btovg1mx-(rho-1)*(rho-1)));
    eta1_btovg[i+200] = -sqrt(fabs(btovg1mx*btovg1mx-(rho-1)*(rho-1)));
    if(rho > xmax) rho = xmax;
    if(rho < xmin) rho = xmin;
    rho1_btovg[i]     =  rho;
    rho1_btovg[i+200] =  rho;
    // top quarter of circle, btovg - 1sigma
    rho = xmax+i*(1-btovg1mn-xmax)/99.;
    eta1_btovg[i+100] =  sqrt(fabs(btovg1mn*btovg1mn-(rho-1)*(rho-1)));
    eta1_btovg[i+300] = -sqrt(fabs(btovg1mn*btovg1mn-(rho-1)*(rho-1)));
    if(rho > xmax) rho = xmax;
    if(rho < xmin) rho = xmin;
    rho1_btovg[i+100] =  rho;
    rho1_btovg[i+300] =  rho;
  }
  // 2 sigmas circles
  for(int i=0; i<100; i++) {
    // top quarter of circle, btovg + 2sigma
    rho = (1-btovg2mx)+i*(xmax-(1-btovg2mx))/99.;
    if(btovg2mx*btovg2mx-(rho-1)*(rho-1) >= 0.) {
      eta2_btovg[i]     =  sqrt(btovg2mx*btovg2mx-(rho-1)*(rho-1));
      eta2_btovg[i+200] = -sqrt(btovg2mx*btovg2mx-(rho-1)*(rho-1));
    } else {
      eta2_btovg[i]     =  0.;
      eta2_btovg[i+200] =  0.;
    }      
    if(rho > xmax) rho = xmax;
    if(rho < xmin) rho = xmin;
    rho2_btovg[i]     =  rho;
    rho2_btovg[i+200] =  rho;    
    // top quarter of circle, btovg - 2sigma
    rho = xmax+i*(1-btovg2mn-xmax)/99.;
    if(btovg2mn*btovg2mn-(rho-1)*(rho-1) >= 0.) {
      eta2_btovg[i+100] =  sqrt(btovg2mn*btovg2mn-(rho-1)*(rho-1));
      eta2_btovg[i+300] = -sqrt(btovg2mn*btovg2mn-(rho-1)*(rho-1));
    } else {
      eta2_btovg[i+100] = 0.;
      eta2_btovg[i+300] = 0.;
    }
    if(rho > xmax) rho = xmax;   
    if(rho < xmin) rho = xmin;
    rho2_btovg[i+100] =  rho;
    rho2_btovg[i+300] =  rho;

  }
}

////////////////////////////////////////////////////
//   Calculates circles for B2taunu, corresponding to //
//   B2taunu mean +/- 1sigma and 2sigma               //
////////////////////////////////////////////////////

void prepare_b2taunu() {
  double rho;
  // 1 sigma circles
  for(int i=0; i<100; i++) {
    // top half circle, B2taunu + 1sigma
    rho = -b2taunu1mx+i*2*b2taunu1mx/99.;
    eta1_b2taunu[i] = sqrt(fabs(b2taunu1mx*b2taunu1mx-rho*rho));
    rho1_b2taunu[i] = rho;
    // bottom half circle, B2taunu + 1sigma
    rho = b2taunu1mx-i*2*b2taunu1mx/99.;
    eta1_b2taunu[i+100] = -sqrt(fabs(b2taunu1mx*b2taunu1mx-rho*rho));
    rho1_b2taunu[i+100] = rho;
    // top half circle, B2taunu - 1sigma
    rho = -b2taunu1mn+i*2*b2taunu1mn/99.;
    eta1_b2taunu[i+200] = sqrt(fabs(b2taunu1mn*b2taunu1mn-rho*rho));
    rho1_b2taunu[i+200] = rho;
    // bottom half circle, B2taunu - 1sigma
    rho = b2taunu1mn-i*2*b2taunu1mn/99.;
    eta1_b2taunu[i+300] = -sqrt(fabs(b2taunu1mn*b2taunu1mn-rho*rho));
    rho1_b2taunu[i+300] = rho;
  }

  // 2 sigma circles
  for(int i=0; i<100; i++) {
    // top half circle, B2taunu + 2sigma
    rho = -b2taunu2mx+i*2*b2taunu2mx/99.;
    eta2_b2taunu[i] = sqrt(fabs(b2taunu2mx*b2taunu2mx-rho*rho));
    rho2_b2taunu[i] = rho;
    // bottom half circle, B2taunu + 2sigma
    rho = b2taunu2mx-i*2*b2taunu2mx/99.;
    eta2_b2taunu[i+100] = -sqrt(fabs(b2taunu2mx*b2taunu2mx-rho*rho));
    rho2_b2taunu[i+100] = rho;
    // top half circle, B2taunu - 2sigma
    rho = -b2taunu2mn+i*2*b2taunu2mn/99.;
    eta2_b2taunu[i+200] = sqrt(fabs(b2taunu2mn*b2taunu2mn-rho*rho));
    rho2_b2taunu[i+200] = rho;
    rho = b2taunu2mn-i*2*b2taunu2mn/99.;
    // bottom half circle, B2taunu - 2sigma
    eta2_b2taunu[i+300] = -sqrt(fabs(b2taunu2mn*b2taunu2mn-rho*rho));
    rho2_b2taunu[i+300] = rho;
  }
}

////////////////////////////////////////////////////
//   Calculates circles for Vub, corresponding to //
//   Vub mean +/- 1sigma and 2sigma               //
////////////////////////////////////////////////////

void prepare_Vub() {
  double rho;
  // 1 sigma circles
  for(int i=0; i<100; i++) {
    // top half circle, Vub + 1sigma
    rho = -vub1mx+i*2*vub1mx/99.;
    eta1_vub[i] = sqrt(fabs(vub1mx*vub1mx-rho*rho));
    rho1_vub[i] = rho;
    // bottom half circle, Vub + 1sigma
    rho = vub1mx-i*2*vub1mx/99.;
    eta1_vub[i+100] = -sqrt(fabs(vub1mx*vub1mx-rho*rho));
    rho1_vub[i+100] = rho;
    // top half circle, Vub - 1sigma
    rho = -vub1mn+i*2*vub1mn/99.;
    eta1_vub[i+200] = sqrt(fabs(vub1mn*vub1mn-rho*rho));
    rho1_vub[i+200] = rho;
    // bottom half circle, Vub - 1sigma
    rho = vub1mn-i*2*vub1mn/99.;
    eta1_vub[i+300] = -sqrt(fabs(vub1mn*vub1mn-rho*rho));
    rho1_vub[i+300] = rho;
  }

  // 2 sigma circles
  for(int i=0; i<100; i++) {
    // top half circle, Vub + 2sigma
    rho = -vub2mx+i*2*vub2mx/99.;
    eta2_vub[i] = sqrt(fabs(vub2mx*vub2mx-rho*rho));
    rho2_vub[i] = rho;
    // bottom half circle, Vub + 2sigma
    rho = vub2mx-i*2*vub2mx/99.;
    eta2_vub[i+100] = -sqrt(fabs(vub2mx*vub2mx-rho*rho));
    rho2_vub[i+100] = rho;
    // top half circle, Vub - 2sigma
    rho = -vub2mn+i*2*vub2mn/99.;
    eta2_vub[i+200] = sqrt(fabs(vub2mn*vub2mn-rho*rho));
    rho2_vub[i+200] = rho;
    rho = vub2mn-i*2*vub2mn/99.;
    // bottom half circle, Vub - 2sigma
    eta2_vub[i+300] = -sqrt(fabs(vub2mn*vub2mn-rho*rho));
    rho2_vub[i+300] = rho;
  }
}

void prepareopt(int optin[6], int optout[6]) {
  if (optin[1] && optin[1]< 10) {
    optout[optin[1]-1] = optin[2];
  }
  if (optin[1] == 12) {
    optout[0] = optin[2];//+1;
    optout[1] = optin[2];
  }
  if (optin[1] == 23) {
    optout[1] = optin[2]+1;
    optout[2] = optin[2];
  }
  if (optin[1] == 123) {
    optout[0] = optin[2]+1;
    optout[1] = optin[2];
    optout[2] = optin[2]+2;
  }

  int color = optin[4];
  //if(optin[0] == 1) color = GetLevel(300+optin[4]);

  if (optin[3] && optin[3]< 10) {
    optout[optin[3]+2] = color;
  }
  if (optin[3] == 12) {
    optout[3] = color+1;
    optout[4] = color;
  }
  if (optin[3] == 23) {
    optout[4] = color+1;
    optout[5] = color;
  }
  if (optin[3] == 123) {
    optout[3] = color+1;
    optout[4] = color;
    optout[5] = color+2;
  }
}

void drawfunc(double rho[400], double eta[400], int nsind, string DrawOpts) {

  if (DrawOpts.compare("CONT")==0) {
    TGraph * line = new TGraph(400, rho, eta);
    line->SetLineColor(nsind);
    line->Draw();
  } else if (DrawOpts.compare("AREA")==0) {
    TGraph * area = new TGraph(400, rho, eta);
    area->SetFillColor(nsind);
    area->Draw("F");
  } else if (DrawOpts.compare("DASHED")==0) {
    TGraph * line = new TGraph(400, rho, eta);
    line->SetLineColor(nsind);
    line->SetLineStyle(2);
    line->Draw();
  }    
}

void drawfunc200(double rho[400], double eta[400], int nsind, string DrawOpts) {
  if (DrawOpts.compare("CONT")==0) {
    TGraph * line = new TGraph(200, rho, eta);
    line->SetLineColor(nsind);
    line->Draw();
  } else if (DrawOpts.compare("AREA")==0) {
    TGraph * area = new TGraph(200, rho, eta);
    area->SetFillColor(nsind);
    area->Draw("F");
  }
}

void setdrawlevel(TH2D* Constrain, double level[5], int sel[6], string DrawOpts) {
  if ( DrawOpts.compare("AREA")==0 ) {
    if (sel[5] > 0) {
      if (sel[3] > 0 && sel[4] > 0) {
	drawcont(Constrain,level,2,17,"AREA");
	drawcont(Constrain,level,1,18,"AREA");
	drawcont(Constrain,level,0,19,"AREA");
      } else if (sel[4] > 0) {
	drawcont(Constrain,level,2,sel[5],"AREA");
	drawcont(Constrain,level,1,sel[4],"AREA");
      }	else {
	drawcont(Constrain,level,2,sel[5],"AREA");
      }
    } else if (sel[4] > 0) {
      if (sel[3] > 0) {
	drawcont(Constrain,level,1,sel[4],"AREA");
	drawcont(Constrain,level,0,sel[3],"AREA");
      } else {
	drawcont(Constrain,level,1,sel[4],"AREA");
      }
    } else if (sel[3] > 0) {
      drawcont(Constrain,level,0,sel[3],"AREA");
    }
  } else if ( DrawOpts.compare("CONT")==0 ) {
    if (sel[0] > 0) {
      drawcont(Constrain,level,0,sel[0],"CONT");
    }
    if (sel[1] > 0) {
      drawcont(Constrain,level,1,sel[1],"CONT");
    }
    if (sel[2] > 0) {
      drawcont(Constrain,level,2,sel[2],"CONT");
    }
  }
}

//////////////////////////////////////////////////
//   given a rho vs eta histo and a level on z  //    
//     coordinate, returns the region having    //
//   z > level. It can be used to get 1sigma    //
//   or 2sigma region for not analitical bounds //
//////////////////////////////////////////////////

void drawcont(TH2D* ThisConstrain, double level[5], int nlevel, int color, string DrawOpts) {
  // set the level
  // read Plot infos
  double nrho   = ThisConstrain->GetNbinsX();
  double neta   = ThisConstrain->GetNbinsY(); 

  int pal = 0;
  if (DrawOpts.compare("AREA")==0) {
    if (color >=300) {
      for (int kol=0; kol<20; kol++) {
	if (color == Palette[kol]) pal = kol;
      }
    } else
      pal = color;
  } else
    pal = color;
  
  TH2D* ThisConstrain2 = smooth2D(ThisConstrain,smooth, level, nlevel);
  
  if(DrawOpts.compare("AREA")==0) {
    // set to 1 (0) z coordinate for cells (not) satisfying
    // the requirement z(rho,eta) > level
    for(int ix=1; ix<=nrho; ix++) {
      for(int iy=1; iy<=neta; iy++) {
	if(ThisConstrain2->GetBinContent(ix,iy) <= level[4-nlevel]) {
	  ThisConstrain2->SetBinContent(ix,iy,0.);
	} else {
	  ThisConstrain2->SetBinContent(ix,iy,pal);
	}
       }
    }
    
    // set the color for the bound and overimpose it to the plot
    ThisConstrain2->SetMaximum(20.);
    ThisConstrain2->Draw("COLSAME");
  } else if(DrawOpts.compare("CONT")==0) {
    ThisConstrain2->SetLineColor(color);
    ThisConstrain2->SetContour(1);
    ThisConstrain2->SetContourLevel(0,level[4-nlevel]);
    //ThisConstrain2->SetLineStyle(4);
    ThisConstrain2->SetLineWidth(2);
    ThisConstrain2->Draw("cont3same");
  } else {
    cout << "Do you want Area or Contours?" << endl;
  }
  return;
}

void SetDefaultPalette(int other) {

  for(int i=0;i<20;i++) {
    Palette[i] = 0;
  }
  if (other==0) {
    Palette[1] = 300+set_vub[4];
    Palette[2] = 300+set_dmd[4];
    Palette[3] = 300+set_eps[4];
    Palette[4] = 300+set_dms[4];
    Palette[5] = 300+set_s2b[4];
    Palette[6] = 300+set_s2bpg[4];
    Palette[7] = 300+set_s2a[4];
    Palette[8] = 300+set_k2pnn[4];
    Palette[9] = 300+set_gamma[4];
    Palette[10]= 300+set_c2b[4]; 

    Palette[11] = 300+set_vub[4]+1;
    Palette[12] = 300+set_dmd[4]+1;
    Palette[13] = 300+set_eps[4]+1;
    Palette[14] = 300+set_dms[4]+1;
    Palette[15] = 300+set_s2b[4]+1;
    Palette[16] = 300+set_s2bpg[4]+1;
    Palette[17] = 300+set_s2a[4]+1;
    Palette[18] = 300+set_k2pnn[4]+1;
    Palette[19] = 300+set_gamma[4]+1;
    
    // for 123 plots
    //Palette[17] = 8; // green  = 99%
    //Palette[18] = 5; // yellow = 95%
    //Palette[19] = 2; // red    = 68%
    //Palette[18] = 14; // grey   = 95%
    //Palette[19] = 4; // blue    = 68%
    //Palette[18] = 398; // customized   = 95%
    //Palette[19] = 399; // customized   = 68%
  } else {    
    Palette[1] = 1;
    Palette[2] = 2;
    Palette[3] = 3;
    Palette[4] = 4;
    Palette[5] = 5;
    Palette[6] = 6;
    Palette[7] = 7;
    Palette[8] = 8;
    Palette[9] = 9;
    Palette[10]= 390;
    Palette[11] = 391;
    Palette[12] = 392;
    Palette[13] = 393;
    Palette[14] = 394;
    Palette[15] = 395;
    Palette[16] = 396;
    Palette[17] = 397;
    Palette[18] = 398;
    Palette[19] = 399;
    
    // for 123 plots
    //Palette[17] = 8; // green  = 99%
    //Palette[18] = 5; // yellow = 95%
    //Palette[19] = 2; // red    = 68%
    //Palette[18] = 14; // grey   = 95%
    //Palette[19] = 4; // blue    = 68%
    //Palette[18] = 398; // customized   = 95%
    //Palette[19] = 399; // customized   = 68%
  }
  gStyle->SetPalette(20,Palette);
}

void SetOldPalette() {
   for(int i=0;i<20;i++) {
    Palette[i] = i;
  }
   gStyle->SetPalette(20,Palette);
}


int GetLevel(int indcol) {
  for(int i=0;i<20;i++) {
    if (Palette[i] == indcol)
      return i;
  }
  return 0;
}

void Read_Parameters(const char* filename) {
  char buffer[256];  
  ifstream reader (filename);
  map<string, int> data; 
  map<string, string> files; 
  string label, equal, svalue;
  char cvalue[40];
  int value;
  
  if ( ! reader.is_open())
    { cout << "Error opening file"; exit (1); }
  
  while ( reader.getline (buffer,256) ) {
    istringstream i(buffer);
    i >> label
      >> equal
      >> cvalue;
    if (equal == "=") {
      value = atoi(cvalue);
      data[label] = value;
    } else if (equal == "==") {
      svalue = (string) cvalue;
      files[label] = svalue;
    }
    equal = "";
  }
  map<string, int>::const_iterator iter;
  for (iter=data.begin(); iter != data.end(); iter++) {
    cout << iter->first << " " << iter->second << endl;
  }
  
  Assign_Parameters(data);
  Assign_Parameters(files);
  
  return;
}

void Assign_Parameters(map<string, string> files) {
  // not used yet... but it'here if we need it
  //plotname   = files["plotname"];
  //dirname   = files["dirname"];
  gamma_use_root = files["gamma_use_root"];
  season         = files["season"];
  addtxt         = files["addtxt"];
}

void Assign_Parameters(map<string, int> data) {
  
  graphic       = data["interactive_root_window"];
  makeps        = data["create_output_file"];
  makepdf       = data["create_output_file_pdf"];
  drawtotarea   = data["draw_total_area_color"];
  nocol         = data["nocol"];
  logoposition  = data["logo_position"];
  positiveeta   = data["positiveeta"];
  fitnumber     = data["fitnumber"];
  s2a_use_his   = data["s2a_use_his"];
  gamma_use_his = data["gamma_use_his"];
  s2bpg_use_his = data["s2bpg_use_his"];
  btovg_use_his = data["btovg_use_his"];
  drawtriangle  = data["drawtriangle"];
  doNP          = data["doNP"];
  doMFV         = data["doMFV"];
  zoom          = data["zoom"];
  smooth        = data["smooth"];

  set_dmd[0]  = data["dmd_draw_const"];
  set_dmd[1]  = data["dmd_sigma_line"];
  set_dmd[2]  = data["dmd_color_line"];
  set_dmd[3]  = data["dmd_sigma_area"];
  set_dmd[4]  = data["dmd_color_area"];
  set_dmd[5]  = data["dmd_draw_title"];
  set_eps[0]  = data["eps_draw_const"];
  set_eps[1]  = data["eps_sigma_line"];
  set_eps[2]  = data["eps_color_line"];
  set_eps[3]  = data["eps_sigma_area"];
  set_eps[4]  = data["eps_color_area"];
  set_eps[5]  = data["eps_draw_title"];
  set_dms[0]  = data["dms_draw_const"];
  set_dms[1]  = data["dms_sigma_line"];
  set_dms[2]  = data["dms_color_line"];
  set_dms[3]  = data["dms_sigma_area"];
  set_dms[4]  = data["dms_color_area"];
  set_dms[5]  = data["dms_draw_title"];
  set_vub[0]  = data["vub_draw_const"];
  set_vub[1]  = data["vub_sigma_line"];
  set_vub[2]  = data["vub_color_line"];
  set_vub[3]  = data["vub_sigma_area"];
  set_vub[4]  = data["vub_color_area"];
  set_vub[5]  = data["vub_draw_title"];
  set_s2b[0]  = data["s2b_draw_const"];
  set_s2b[1]  = data["s2b_sigma_line"];
  set_s2b[2]  = data["s2b_color_line"];
  set_s2b[3]  = data["s2b_sigma_area"];
  set_s2b[4]  = data["s2b_color_area"];
  set_s2b[5]  = data["s2b_draw_title"];
  set_c2b[0]  = data["c2b_draw_const"];
  set_c2b[1]  = data["c2b_sigma_line"];
  set_c2b[2]  = data["c2b_color_line"];
  set_c2b[3]  = data["c2b_sigma_area"];
  set_c2b[4]  = data["c2b_color_area"];
  set_c2b[5]  = data["c2b_draw_title"];
  set_s2bpg[0] = data["s2bpg_draw_const"];
  set_s2bpg[1] = data["s2bpg_sigma_line"];
  set_s2bpg[2] = data["s2bpg_color_line"];
  set_s2bpg[3] = data["s2bpg_sigma_area"];
  set_s2bpg[4] = data["s2bpg_color_area"];
  set_s2bpg[5] = data["s2bpg_draw_title"];
  set_s2a[0] = data["s2a_draw_const"];
  set_s2a[1] = data["s2a_sigma_line"];
  set_s2a[2] = data["s2a_color_line"];
  set_s2a[3] = data["s2a_sigma_area"];
  set_s2a[4] = data["s2a_color_area"];
  set_s2a[5] = data["s2a_draw_title"];
  set_k2pnn[0] = data["k2pnn_draw_const"];
  set_k2pnn[1] = data["k2pnn_sigma_line"];
  set_k2pnn[2] = data["k2pnn_color_line"];
  set_k2pnn[3] = data["k2pnn_sigma_area"];
  set_k2pnn[4] = data["k2pnn_color_area"];
  set_k2pnn[5] = data["k2pnn_draw_title"];
  set_gamma[0] = data["gamma_draw_const"];
  set_gamma[1] = data["gamma_sigma_line"];
  set_gamma[2] = data["gamma_color_line"];
  set_gamma[3] = data["gamma_sigma_area"];
  set_gamma[4] = data["gamma_color_area"];
  set_gamma[5] = data["gamma_draw_title"];
  set_b2taunu[0] = data["b2taunu_draw_const"];
  set_b2taunu[1] = data["b2taunu_sigma_line"];
  set_b2taunu[2] = data["b2taunu_color_line"];
  set_b2taunu[3] = data["b2taunu_sigma_area"];
  set_b2taunu[4] = data["b2taunu_color_area"];
  set_b2taunu[5] = data["b2taunu_draw_title"];
  set_btovg[0]  = data["btovg_draw_const"];
  set_btovg[1]  = data["btovg_sigma_line"];
  set_btovg[2]  = data["btovg_color_line"];
  set_btovg[3]  = data["btovg_sigma_area"];
  set_btovg[4]  = data["btovg_color_area"];
  set_btovg[5]  = data["btovg_draw_title"];
  set_D0p0[0]  = data["D0p0_draw_const"];
  set_D0p0[1]  = data["D0p0_sigma_line"];
  set_D0p0[2]  = data["D0p0_color_line"];
  set_D0p0[3]  = data["D0p0_sigma_area"];
  set_D0p0[4]  = data["D0p0_color_area"];
  set_D0p0[5]  = data["D0p0_draw_title"];
  set_kl2p0nn[0] = data["kl2p0nn_draw_const"];
  set_kl2p0nn[1] = data["kl2p0nn_sigma_line"];
  set_kl2p0nn[2] = data["kl2p0nn_color_line"];
  set_kl2p0nn[3] = data["kl2p0nn_sigma_area"];
  set_kl2p0nn[4] = data["kl2p0nn_color_area"];
  set_kl2p0nn[5] = data["kl2p0nn_draw_title"];
  set_BtoKpp[0] = data["BtoKpp_draw_const"];
  set_BtoKpp[1] = data["BtoKpp_sigma_line"];
  set_BtoKpp[2] = data["BtoKpp_color_line"];
  set_BtoKpp[3] = data["BtoKpp_sigma_area"];
  set_BtoKpp[4] = data["BtoKpp_color_area"];
  set_BtoKpp[5] = data["BtoKpp_draw_title"];
}

void DefineColors() {
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
  //  color[18] = new TColor(1318, 0.95, 1.0, 1.0, "");
  //  color[19] = new TColor(1319, 0.7, 1.0, 0.8, "");

  color[18] = new TColor(1318, 0.95, 1.0, 1.0, "");
  color[19] = new TColor(1319, 0.7, 1.0, 0.8, "");

  //  color[18] = new TColor(1318, 0.85, 1.0, 0.90, "");
  //  color[19] = new TColor(1319, 0.57, 1.0, 0.62, "");
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
  //color[90] = new TColor(1390, 1.00, 0.41, 0.71, "");
  //color[91] = new TColor(1391, 1.00, 0.08, 0.58, "");
  //color[92] = new TColor(1392, 1.00, 0.75, 0.80, "");
  //color[93] = new TColor(1393, 1.00, 0.71, 0.76, "");
  //color[94] = new TColor(1394, 0.86, 0.44, 0.58, "");
  //color[95] = new TColor(1395, 0.69, 0.19, 0.38, "");
  //color[96] = new TColor(1396, 0.78, 0.08, 0.52, "");
  //color[97] = new TColor(1397, 0.82, 0.13, 0.56, "");

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

}

double CalcThreshold(TH1D* InputPlot, double conf_level) {

  TH1D* PLot1D = (TH1D*) InputPlot->Clone();
  int nx = PLot1D->GetNbinsX();
  PLot1D->SetBinContent(0,0.);
  PLot1D->SetBinContent(1+nx,0.);
  conf_level *= PLot1D->GetSum();

  vector< Double_t > OrderedGArray1D;
  for (int ix=1; ix<=nx; ix++) {
    OrderedGArray1D.push_back(PLot1D->GetBinContent(ix));
  }
  
  double level1D = 0.;
  //double area1D = 0.;
  sort( OrderedGArray1D.begin(), OrderedGArray1D.end(),greater<double>() );
  vector< Double_t> SumArray1D;
  SumArray1D.push_back(OrderedGArray1D[0]);
  for (int i=1; i< (int) OrderedGArray1D.size(); i++) {
    SumArray1D.push_back(SumArray1D[i-1] + (double) OrderedGArray1D[i]);
    if (SumArray1D[i] <= conf_level) {
      level1D = OrderedGArray1D[i];
      //area1D = SumArray1D[i];
    }
  }
  delete PLot1D;
  return level1D;
}

double CalcLevel(TH1D* InputPlot, double conf_level, string sides) {

  TH1D* PLot1D = (TH1D*) InputPlot->Clone();
  int nx = PLot1D->GetNbinsX();
  PLot1D->SetBinContent(0,0.);
  PLot1D->SetBinContent(1+nx,0.);
  
  Double_t Sum = PLot1D->GetSum();
  
  PLot1D->Scale(1./Sum);
  double level1D = 0.;
  level1D = CalcThreshold(PLot1D, conf_level);
  
  int maxbin = PLot1D->GetMaximumBin();
  int levelbin = 0;
  
  if(sides.compare("Low")==0) {
    int found = 0;
    for(int ix=1; ix<maxbin; ix++) {
      int ibin = maxbin-ix;
      if(found == 0 &&
	 (double) PLot1D->GetBinContent(ibin) > level1D &&
	 (double) PLot1D->GetBinContent(ibin-1) <= level1D) {
	levelbin = ibin;
	found = 1;
      }
    }
  } else if(sides.compare("Up")==0) {
    int found = 0;
    for(int ix=maxbin; ix<nx; ix++) {
      if(found == 0 &&
	 (double) PLot1D->GetBinContent(ix-1) > level1D &&
	 (double) PLot1D->GetBinContent(ix) <= level1D) {
	levelbin = ix-1;
	found = 1;
      }
    }
  } else if(sides.compare("Low2")==0) {
    int found = 0;
    for(int ix=1; ix<maxbin; ix++) {
      int ibin = maxbin-ix;
      if(found < 3 &&
	 (double) PLot1D->GetBinContent(ibin) > level1D &&
	 (double) PLot1D->GetBinContent(ibin-1) <= level1D) {
	levelbin = ibin;
	found ++;
      }
    }
  } else if(sides.compare("Up2")==0) {
    int found = 0;
    for(int ix=maxbin; ix<nx; ix++) {
      if(found < 2 &&
	 (double) PLot1D->GetBinContent(ix-1) > level1D &&
	 (double) PLot1D->GetBinContent(ix) <= level1D) {
	levelbin = ix-1;
	found ++;
      }
    }
  } else if(sides.compare("Low3")==0) {
    int found = 0;
    for(int ix=1; ix<180; ix++) {
      int ibin = 180-ix;
      if(found < 5 &&
	 (double) PLot1D->GetBinContent(ibin) > level1D &&
	 (double) PLot1D->GetBinContent(ibin-1) <= level1D) {
	levelbin = ibin;
	found ++;
      }
    }
  } else if(sides.compare("Up3")==0) {
    int found = 0;
    for(int ix=maxbin; ix<nx; ix++) {
      if(found < 4 &&
	 (double) PLot1D->GetBinContent(ix-1) > level1D &&
	 (double) PLot1D->GetBinContent(ix) <= level1D) {
	levelbin = ix-1;
	found ++;
      }
    }
  } else {
    cout << "Error in CalcLevel; specify Upper or Lower level" << endl;
    return 0.;
  }
  
  TAxis* tmp = PLot1D->GetXaxis();
  double x_min = tmp->GetXmin();
  double bin_width = tmp->GetBinWidth(3);

  if (sides.compare("Low")==0)
    return x_min + (levelbin-1)*bin_width + bin_width*0.01;
  else
    return x_min + (levelbin-1)*bin_width + bin_width*0.99;
}


void makeOtherPlot(TH2D* histo2D, TObjArray* histo, int superimpose, int col68, int col95) {

  TAxis* tmp = (TAxis*) histo2D->GetXaxis();
  xmin = tmp->GetXmin();
  xmax = tmp->GetXmax();
  binx =  tmp->GetNbins();
  
  tmp = (TAxis*) histo2D->GetYaxis();
  ymin = tmp->GetXmin();
  // to avoid problems with common code
  etamin=ymin;
  ymax = tmp->GetXmax();
  biny =  tmp->GetNbins();
  
  TH2D* null2 = new TH2D("null","",binx,xmin,xmax,biny,ymin,ymax);
  //  TH2D* null2 = new TH2D("null","",binx,xmin,xmax,biny,-180.,180.);
  //  TH2D* null2 = new TH2D("null","",100,0.0,3.0,100,0.,0.6);
  if(strncmp(plotname.c_str(),"gammavsRb",9)==0) {
    null2->SetXTitle("R_{b}");
    null2->SetYTitle("#gamma");
  }
  if(strncmp(plotname.c_str(),"FbvsBk",6)==0) {
    null2->SetXTitle("B_{k}");
    null2->SetYTitle("f_{B_{d}}#sqrt{B_{B_{d}}}[GeV]");
  }
  if(strncmp(plotname.c_str(),"FbsvsBk",7)==0) {
    null2->SetXTitle("B_{k}");
    null2->SetYTitle("f_{B_{s}}#sqrt{B_{B_{s}}}[GeV]");
  }
  if(strncmp(plotname.c_str(),"BkvsXi",6)==0) {
    null2->SetXTitle("#xi");
    null2->SetYTitle("B_{k}");
  }
  if(strncmp(plotname.c_str(),"FbvsXi",6)==0) {
    null2->SetXTitle("#xi");
    null2->SetYTitle("f_{B_{d}}#sqrt{B_{B_{d}}}[GeV]");
  }
  if(strncmp(plotname.c_str(),"FbsvsXi",7)==0) {
    null2->SetXTitle("#xi");
    null2->SetYTitle("f_{B_{s}}#sqrt{B_{B_{s}}}[GeV]");
  }
  
  if (lxlab) null2->SetXTitle(xlab);
  if (lylab) null2->SetYTitle(ylab);
  
  null2->SetTitleSize(0.07,"X");
  null2->SetTitleSize(0.07,"Y");
  null2->SetNdivisions(505, "X");
  if (lsquare) {
    null2->SetTitleOffset(0.85,"X");
    null2->SetTitleOffset(1.20,"Y");
  } else {
    null2->SetTitleOffset(0.8,"X");
    null2->SetTitleOffset(0.8,"Y");
  }
  null2->SetStats(0);
  if (!superimpose){
    null2->Draw();
    gPad->SetFillColor(0);
    gPad->SetFrameFillColor(0);
  }

  if (drawtotarea == 1) {
    drawFromGraph(histo2D,histo, 0, "AREA", col95);
    drawFromGraph(histo2D,histo, 1, "AREA", col68);
    //drawFromGraph(histo, 0, "AREA", col95);
    //drawFromGraph(histo, 1, "AREA", col68);
  }

  if(llines) {
    //    drawFromGraph(histo, 0, "CONT", col95);
    drawFromGraph(histo2D,histo, 0, "CONT", col68);
    drawFromGraph(histo2D,histo, 1, "CONT", col68);
    //drawFromGraph(histo, 0, "CONT", col68);
    //drawFromGraph(histo, 1, "CONT", col68);
  }
  
  null2->Draw("same");
  //theapp->Run(true);

  // draw error bars
  if (xval!=-999 && yval!=-999 && !superimpose){
    TLine *lx = new TLine();
    lx->SetLineWidth(2);
    double zero=0,err;
    err = xerr;
    TGraphErrors *g1 = new TGraphErrors(1,&xval,&yval,&err,&zero);
    g1->SetLineWidth(2);
    g1->SetLineStyle(1);
    g1->SetMarkerStyle(20);
    g1->Draw("P");
    
    err = xerr;
    double min_x = max(xval-err,xmin);
    double max_x = min(xval+err,xmax);
    lx->DrawLine(min_x,yval,max_x,yval);

    TLine *ly = new TLine();
    ly->SetLineWidth(2);
    err = yerr;
    TGraphErrors* g2 = new TGraphErrors(1,&xval,&yval,&zero,&err);
    g2->SetLineWidth(2);
    g2->SetLineStyle(1);
    g2->SetMarkerStyle(20);
    g2->Draw("P");
    
    double min_y = max(yval-err,ymin);
    double max_y = min(yval+err,ymax);
    ly->DrawLine(xval,min_y,xval,max_y);
  }

  if(llogo && !superimpose) {
    PutLogo(logoposition,0,0);
  }

   TLine* Bound = new TLine(xmin,ymax,xmax,ymax);
   Bound->Draw();
   Bound->DrawLine(xmax,ymin,xmax,ymax);

   //   TLine* Bound = new TLine(0.0,0.6,3.0,0.6);
   //   Bound->Draw();
   //   Bound->DrawLine(3.0,0.0,3.0,0.6);

  //  Bound->DrawLine(xmin,0.,xmax,0.);
  //  Bound->DrawLine(0.,ymin,0.,ymax);

  //TLine* Bound = new TLine(xmin,180.,xmax,180.);
  //Bound->Draw();
  //Bound->DrawLine(xmax,-180.,xmax,180.);
  
  c1->RedrawAxis("same");
}


double CalcLevel_Mediana(TH1D* InputPlot, double conf_level, string sides) {
  
  TH1D* Plot = (TH1D*) InputPlot->Clone();
  int nbin = Plot->GetNbinsX();
  double low = Plot->GetXaxis()->GetXmin();
  double up  = Plot->GetXaxis()->GetXmax();
  double width=(up-low)/nbin;

  // delete all the underflow and overflow, otherwise fucked root 
  // take them into account in the normalization
  Plot->SetBinContent(0,0.);
  Plot->SetBinContent(1+nbin,0.);

  double *NArray = new double[nbin];
  for (int i=0;i<nbin;i++) {
    NArray[i]=Plot->GetBinContent(i+1);
  }

  //compute coordinates
  double *coord = new double[nbin];
  for (int i=0;i<nbin;i++){
    coord[i] = low+width*i;
  }

  cout << coord[10] << endl;

  double sumold,sum=0;
  double chk; 
  double xcoo = -99.;

  // compute CLs and the corresponding coordinates
  //
  xcoo = low;
  if(sides.compare("Low")) {
    chk = (1.-conf_level)/2.*Plot->GetSum();
  } else {
    chk = (1.-(1.-conf_level)/2.)*Plot->GetSum();
    xcoo = up;
  }

  sum = NArray[0];
  for (int i=1;i<nbin;i++){
    sumold = sum;
    sum += NArray[i];
    if ( sum>=chk  && sumold<chk ){
      //      double tobeadded=(chk -sumold)/NArray[i];
      //      xcoo =coord[i] +tobeadded*width;
      xcoo =coord[i] +0.5*width;
    }
  }
    
  return xcoo;

}

TH2D* smooth2D(TH2D* ThisConstrain, int smooth, double level[5], int nlevel) {
  double nrho   = ThisConstrain->GetNbinsX();
  double neta   = ThisConstrain->GetNbinsY(); 
  TH2D* ThisConstrain2 = (TH2D*) ThisConstrain->Clone();

  for(int jj=0; jj<smooth; jj++) {
    for (int ix=2; ix<nrho; ix++) {
      for (int iy=2; iy<neta; iy++) {
	double new_weight = (ThisConstrain2->GetBinContent(ix+1,iy+1)+
			     ThisConstrain2->GetBinContent(ix+1,iy)+
			     ThisConstrain2->GetBinContent(ix+1,iy-1)+
			     ThisConstrain2->GetBinContent(ix,iy+1)+
			     ThisConstrain2->GetBinContent(ix,iy)+
			     ThisConstrain2->GetBinContent(ix,iy-1)+
			     ThisConstrain2->GetBinContent(ix-1,iy+1)+
			     ThisConstrain2->GetBinContent(ix-1,iy)+
			     ThisConstrain2->GetBinContent(ix-1,iy-1))/9.;
	ThisConstrain2->SetBinContent(ix,iy,new_weight);
      }
    }
    
    for (int ix=2; ix<nrho; ix++) {
      for (int iy=2; iy<neta; iy++) {
	if(ThisConstrain2->GetBinContent(ix,iy) <= level[4-nlevel]) {
	  int icont = 0;
	  if(ThisConstrain2->GetBinContent(ix+1,iy+1) > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix+1,iy)   > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix+1,iy-1) > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix,iy+1)   > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix,iy-1)   > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy+1) > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy)   > level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy-1) > level[4-nlevel]) icont++;
	  if(icont >4)
	    ThisConstrain2->SetBinContent(ix,iy,level[4-nlevel]*1.01);  
	}
      }
    }
    
    for (int ix=2; ix<nrho; ix++) {
      for (int iy=2; iy<neta; iy++) {
	if(ThisConstrain2->GetBinContent(ix,iy)     > level[4-nlevel]) {
	  int icont = 0;
	  if(ThisConstrain2->GetBinContent(ix+1,iy+1) <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix+1,iy)   <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix+1,iy-1) <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix,iy+1)   <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix,iy-1)   <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy+1) <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy)   <= level[4-nlevel]) icont++;
	  if(ThisConstrain2->GetBinContent(ix-1,iy-1) <= level[4-nlevel]) icont++;
	  if(icont>6) 
	    ThisConstrain2->SetBinContent(ix,iy,level[4-nlevel]*0.99);
	}
      }
    }  
  }

  return ThisConstrain2;
}


// NEW SCHOOL!!!!!!!!!!!!!!!!!!!!

/////////////////////////////////////////////////////////////////////////
// given rho, eta, a 1D histo and a constraint
// returns the value of the 1D histogram for 
// those values of rho and eta.
// only 2alpha and 2bpg are implemented (for the moment)
/////////////////////////////////////////////////////////////////////////

double findLevel(double rho, double eta, TH1D* Histo, TString name) {
  double var=0;
  if (strncmp(name,"2bpg",4)==0)
    var = acos(cos(2*atan2(eta,1-rho)+atan2(eta,rho)))*180./Pi;
  else if (strncmp(name,"alpha",5)==0) {
    var = (acos(-1.)- atan2(eta,1-rho) - atan2(eta,rho))*180./Pi;
    if(var < 0.) var += 180.;
    if(var > 180.) var -= 180.;
  }
  else if (strncmp(name,"s2a",3)==0)
    var = sin(2*(acos(-1.)- atan2(eta,1-rho) - atan2(eta,rho)));
  else if (strncmp(name,"dms",3)==0)
    var = sqrt(pow(1.-rho,2.) + pow(eta,2.));
  else
    cout << "you should implement the relation with rho and eta" << endl;
  double thislevel = 0.;

  if ((strncmp(name,"2bpg",4)==0) && (var == 180.))
    thislevel = Histo->GetBinContent(Histo->FindBin(var)-1);
  else
    thislevel = Histo->GetBinContent(Histo->FindBin(var));
    
  return thislevel;
}

/////////////////////////////////////////////////////////////////////////
// interface function that produces a 2d plot for rho and eta
// from a 1D plot on a given constraint and calls the general function
// GraphFromHisto to return the TObjArray containing the contours
/////////////////////////////////////////////////////////////////////////

TObjArray* CalcGraph(TH1D* Histo, double level1, double level2, TString name) {
  int rhobin= 1000;
  int imin = 0;
  int imax = 1000;
  int etabin= 600;
  int jmin = 0;
  int jmax = 600;

  double min = xmin;
  double max = xmax;
  double bottom = etamin;
  double top = ymax;

  // this probably works only for doNP==1 and should
  // be changed if you want to create another sise of plot
  // MB has done this change only for alpha (see below)
  if(name == "2bpg1") {
    min = xmin; max = 1.0;  imin = 0;   imax = 900;
    bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300;
  }
  if(name == "2bpg2") {
    min = 1.0;  max = xmax; imin = 900; imax = 1000;
    bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300;
  }
  if(name == "2bpg3") {
    min = xmin; max = 1.0;  imin = 0;   imax = 900;
    bottom = 0; top = ymax;  jmin = 300;   jmax = 600;
  }
  if(name == "2bpg4") {
    min = 1.0;  max = xmax; imin = 900; imax = 1000;
    bottom = 0; top = ymax;  jmin = 300;   jmax = 600;
  }

  if (zoom != 0) {

    if(name == "alpha1") {
      min = xmin; max = 0.0;  imin = 0;   imax = 166;
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300;
    }
    if(name == "alpha2") {
      min = 0.0;  max =1.0;   imin = 166; imax = 834;
    }
    if(name == "alpha3") {
      min = 1.0;  max = xmax; imin = 834; imax = 1000;
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300;
    }
    if(name == "alpha4") {
      min = xmin; max = 0.0;  imin = 0;   imax = 166;
      bottom = 0; top = ymax;  jmin = 300;   jmax = 600;
    }
    if(name == "alpha5") {
      min = 1.0;  max = xmax; imin = 934; imax = 1000;
      bottom = 0; top = ymax;  jmin = 300;   jmax = 600;
    }
  }else if (doNP == 1) {

    if(name == "alpha1") {
      min = xmin; max = 0.0;  imin = 0;   imax = 500; 
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300; 
    }
    if(name == "alpha2") {
      min = 0.0;  max =1.0;   imin = 500; imax = 900;
    }
    if(name == "alpha3") {
      min = 1.0;  max = xmax; imin = 900; imax = 1000;
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 300; 
    }
    if(name == "alpha4") {
      min = xmin; max = 0.0;  imin = 0;   imax = 500;
      bottom = 0; top = ymax;  jmin = 300;   jmax = 600; 
    }
    if(name == "alpha5") {
      min = 1.0;  max = xmax; imin = 900; imax = 1000;
      bottom = 0; top = ymax;  jmin = 300;   jmax = 600; 
    }
  } else {

    if(name == "alpha1") {
      min = xmin; max = 0.0;  imin = 0;   imax = 500; 
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 100; 
    }
    if(name == "alpha2") {
      min = 0.0;  max =1.0;   imin = 500; imax = 900;
    }
    if(name == "alpha3") {
      min = 1.0;  max = xmax; imin = 900; imax = 1000;
      bottom = etamin; top = 0.0;  jmin = 0;   jmax = 100; 
    }
    if(name == "alpha4") {
      min = xmin; max = 0.0;  imin = 0;   imax = 500;
      bottom = 0; top = ymax;  jmin = 100;   jmax = 600; 
    }
    if(name == "alpha5") {
      min = 1.0;  max = xmax; imin = 900; imax = 1000;
      bottom = 0; top = ymax;  jmin = 100;   jmax = 600; 
    }
  }

  if(name == "alpha") {
    etabin= 1000;
    jmax=1000;
  }
  // we want the graph to be larger than the hist, to avoid 
  // problems at the boundaries
  TH2D* tmp2Dhisto = new TH2D(name, name, rhobin, xmin,xmax, etabin,etamin,ymax);

  for (int ix=imin+1; ix<=imax; ix++) {
    double irho = min + double(ix-imin-0.5)/(imax-imin)*(max-min);
    for (int iy=jmin+1; iy<=jmax; iy++) {
      double ieta = bottom + double(iy-jmin-0.5)/(jmax-jmin)*(top-bottom);
      tmp2Dhisto->SetBinContent(ix,iy,findLevel(irho,ieta,Histo,name));
    }
  }
  cout << "GETTING CONTOUR LIST FOR " << name << endl;
  TObjArray* tmpArray = (TObjArray*) GraphFromHisto(tmp2Dhisto,level1,level2);
  delete tmp2Dhisto;
  return tmpArray;
}

/////////////////////////////////////////////////////////////////////////
// Return the contours of a 2D histogram tmp2Dhisto, given the two levels
// l1 and l2 (usually corresponding to 68% and 95%
/////////////////////////////////////////////////////////////////////////

TObjArray* GraphFromHisto(TH2D* tmp2Dhisto, double level1, double level2) {

  double levels[3];  
  //  levels[0] = level2/2; // useless!!! But it works
  levels[0] = level2;
  levels[1] = level1;

  tmp2Dhisto->SetContour(3,levels);
  tmp2Dhisto->Draw("contlist");
  c1->Update();

  TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

  if (contours == NULL){
    printf("*** WARNING: No Contours Were Extracted!\n");
  }

  printf("TotalConts = %d\n", contours->GetSize());
  return contours;
}

void drawFromGraph(TObjArray* contours, int ind, string DrawOpts, int col) {
  TList* contour = (TList*)contours->At(ind);
  printf("TotalGraphs = %d\n", contour->GetSize());

  TGraph* curv = (TGraph*)contour->First();
  for (int i=0; i<contour->GetSize();i++) {

    TGraph* tmpTGraph = CloseTGraph(curv);
    if(DrawOpts.compare("AREA")==0) {
      tmpTGraph->SetLineWidth(2);
      tmpTGraph->SetLineColor(col);
      tmpTGraph->SetFillColor(col);
      tmpTGraph->Draw("F");
    }
    if( DrawOpts.compare("CONT")==0) {
      curv->SetLineWidth(2);
      curv->SetLineColor(col);
      curv->SetFillColor(col);
      curv->Draw();
    }
    curv = (TGraph*)contour->After(curv); // Get Next graph
  }

  //theapp->Run(true);
  c1->Update();
  return;
}

class Point{
public:
  Point(double x, double y){m_x = x; m_y = y;}
  void R(double r){m_r = r;}
  bool operator<(const Point& b) const { return m_r<b.m_r;}
  double m_x;
  double m_y;
  double m_r;
};

void drawFromGraph(TH2D* h2, TObjArray* contours, int ind, string DrawOpts, int col) {

  TList* contour = (TList*)contours->At(ind);
  printf("TotalGraphs = %d\n", contour->GetSize());

  double epsp=1;  
  double epsm=1-1e-6;  
  double x01 = h2->GetBinContent(h2->GetXaxis()->FindBin(xmin*epsp),h2->GetYaxis()->FindBin(ymax*epsm));
  double x11 = h2->GetBinContent(h2->GetXaxis()->FindBin(xmax*epsm),h2->GetYaxis()->FindBin(ymax*epsm));
  double x00 = h2->GetBinContent(h2->GetXaxis()->FindBin(xmin*epsp),h2->GetYaxis()->FindBin(ymin*epsp));
  double x10 = h2->GetBinContent(h2->GetXaxis()->FindBin(xmax*epsm),h2->GetYaxis()->FindBin(ymin*epsp));

  TGraph* curv = (TGraph*)contour->First();
  for (int i=0; i<contour->GetSize();i++) {

    bool l01=false;
    bool l11=false;
    bool l00=false;
    bool l10=false;
    double binwx = h2->GetXaxis()->GetBinWidth(1);
    double binwy = h2->GetYaxis()->GetBinWidth(1);

    double val;
    //double val00=-1,val01=-1,val10=-1,val11=-1,val;

    for(int j=0; j< curv->GetN(); j++) {
      double xx,yy;
      curv->GetPoint(j,xx,yy);
      if (xx>=xmax-binwx/2) curv->SetPoint(j,xmax,yy);
      if (xx<=xmin+binwx/2) curv->SetPoint(j,xmin,yy);
      if (yy>=ymax-binwy/2) curv->SetPoint(j,xx,ymax);
      if (yy<=ymin+binwy/2) curv->SetPoint(j,xx,ymin);
      val = h2->GetBinContent(h2->GetXaxis()->FindBin(xx),h2->GetYaxis()->FindBin(yy));
      if (val<x01){
	l01 = true;
        //val01 = val;
      }
      if (val<x00){
	l00 = true;
        //val00 = val;
      }
      if (val<x10){
	l10 = true;
        //val10 = val;
      }
      if (val<x11){
	l11 = true;
        //val11 = val;
      }
    }

    //  cout << "--------------" << endl;
    //  cout << "x00  " << x00 << endl;
    //  cout << "x01  " << x01 << endl;
    //  cout << "x10  " << x10 << endl;
    //  cout << "x11  " << x11 << endl;
    //  cout << "val  " << val << endl;

    vector<Point> vp;
    double xl,yl;
    curv->GetPoint(curv->GetN()-1,xl,yl);
    if (l01) {
      Point p(xmin,ymax);
      p.R(sqrt(pow(xmin-xl,2)+pow(ymax-yl,2)));
      vp.push_back(p);
    }
    if (l11) {
      Point p(xmax,ymax);
      p.R(sqrt(pow(xmax-xl,2)+pow(ymax-yl,2)));
      vp.push_back(p);
    }
    if (l10) {
      Point p(xmax,ymin);
      p.R(sqrt(pow(xmax-xl,2)+pow(ymin-yl,2)));
      vp.push_back(p);
    }
    if (l00) {
      Point p(xmin,ymin);
      p.R(sqrt(pow(xmin-xl,2)+pow(ymin-yl,2)));
      vp.push_back(p);
    }
    sort(vp.begin(),vp.end());

    if (DrawOpts=="AREA"){
      for (unsigned int i=0;i<vp.size();i++){
	curv->SetPoint(curv->GetN(),vp[i].m_x,vp[i].m_y);
      }
    }

    TGraph* tmpTGraph = CloseTGraph(curv);
    if(DrawOpts.compare("AREA")==0) {
      tmpTGraph->SetLineWidth(2);
      tmpTGraph->SetLineColor(col);
      tmpTGraph->SetFillColor(col);
      tmpTGraph->Draw("F");
    }
    if( DrawOpts.compare("CONT")==0) {
      curv->SetLineWidth(2);
      curv->SetLineColor(col);
      curv->SetFillColor(col);
      curv->Draw();
    }
    curv = (TGraph*)contour->After(curv); // Get Next graph
  }

  c1->Update();
  //  h2->Draw();
  //  theapp->Run(true);
  return;
}


TGraph* CloseTGraph(TGraph* inputgraph) {

  double x_i, x_j, y_i, y_j;
  inputgraph->GetPoint(0,x_i,y_i);
  inputgraph->GetPoint(inputgraph->GetN()-1,x_j,y_j);

  double xleft_i(0), xright_i(0);
  double ybottom_i(0), ytop_i(0);
  double xleft_j(0), xright_j(0);
  double ybottom_j(0), ytop_j(0);

  double deltax = (xmax-xmin)/binx;
  double deltay = (ymax-etamin)/biny;

  if(fabs(x_i-xmin)<deltax) xleft_i   = 1.;
  if(fabs(x_i-xmax)<deltax) xright_i  = 1.;
  if(fabs(y_i-etamin)<deltay) ybottom_i = 1.;
  if(fabs(y_i-ymax)<deltay) ytop_i    = 1.;

  if(fabs(x_j-xmin)<deltax) xleft_j   = 1.;
  if(fabs(x_j-xmax)<deltax) xright_j  = 1.;
  if(fabs(y_j-etamin)<deltay) ybottom_j = 1.;
  if(fabs(y_j-ymax)<deltay) ytop_j    = 1.;

  double xnew[inputgraph->GetN()+3];
  double ynew[inputgraph->GetN()+3];

  for(int i=0; i< inputgraph->GetN(); i++) {
    inputgraph->GetPoint(i,xnew[i],ynew[i]);
  }

  if(xleft_i ==1 && ybottom_j == 1) { 
    // we go from bottom to left
    // passing through the left-bottom corner
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
    xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmin;   ynew[inputgraph->GetN()+2] = y_i;
  } else if(xleft_j ==1 && ybottom_i == 1) { 
    // we go from left to bottom
    // passing through the left-bottom corner
    xnew[inputgraph->GetN()]   = xmin;   ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = x_i;    ynew[inputgraph->GetN()+2] = etamin;    
  } else if(xleft_j ==1 && ytop_i == 1) { 
    // we go from left to top
    // passing through the left-top corner
    xnew[inputgraph->GetN()]   = xmin;   ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = ymax;
    xnew[inputgraph->GetN()+2] = x_i;    ynew[inputgraph->GetN()+2] = ymax;    
  } else if(xleft_i ==1 && ytop_j == 1) { 
    // we go from left to top
    // passing through the left-top corner
    //xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = ymax;
    //xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = ymax;
    //xnew[inputgraph->GetN()+2] = xmin;   ynew[inputgraph->GetN()+2] = y_i;
  } else if(xright_i == 1 && ytop_j ==1) {
    // we go from right to top
    // passing through the tight-top corner
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = ymax;
    xnew[inputgraph->GetN()+1] = xmax;   ynew[inputgraph->GetN()+1] = ymax;
    xnew[inputgraph->GetN()+2] = xmax;   ynew[inputgraph->GetN()+2] = y_i;
  } else if(xright_j == 1 && ytop_i == 1) {
    // we go from top to right
    // passing through the tight-top corner
    xnew[inputgraph->GetN()]   = xmax;   ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmax;   ynew[inputgraph->GetN()+1] = ymax;
    xnew[inputgraph->GetN()+2] = x_i;    ynew[inputgraph->GetN()+2] = ymax;
  } else if(xright_i == 1 && ybottom_j ==1) {
    // we go from right to bottom
    // passing through the right-bottom corner
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
    xnew[inputgraph->GetN()+1] = xmax;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmax;   ynew[inputgraph->GetN()+2] = y_i;
  } else if(xright_j == 1 && ybottom_i ==1) {
    // we go from bottom to right
    // passing through the right-bottom corner
    xnew[inputgraph->GetN()]   = xmax;   ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmax;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = x_i;    ynew[inputgraph->GetN()+2] = etamin;
  } else if(xleft_i == 1 && ybottom_j ==1) {
    // we go from left to bottom
    // passing through the left-bottom corner
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
    xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmin;   ynew[inputgraph->GetN()+2] = y_i;
  } else if(xleft_j == 1 && ybottom_i ==1) {
    // we go from bottom to left
    // passing through the left-bottom corner
    xnew[inputgraph->GetN()]   = xmin;   ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = x_i;    ynew[inputgraph->GetN()+2] = etamin;
//    NP C_Bs vs phi_Bs
  } else if (ybottom_i ==1 && ytop_j ==1){ //from top to bottom
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = ymax;
    xnew[inputgraph->GetN()+1] = 1.0;    ynew[inputgraph->GetN()+1] = ymax;
    xnew[inputgraph->GetN()+2] = 1.0;    ynew[inputgraph->GetN()+2] = etamin;
  } else if (ybottom_j ==1 && ytop_i ==1){ //from bottom to top
    xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
    xnew[inputgraph->GetN()+1] = xmin;    ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmin;    ynew[inputgraph->GetN()+2] = ymax;
    //xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
    //xnew[inputgraph->GetN()+1] = 1.0;    ynew[inputgraph->GetN()+1] = etamin;
    //xnew[inputgraph->GetN()+2] = 1.0;    ynew[inputgraph->GetN()+2] = ymax;

// achilleplot Bd
  } else if (xleft_i ==1 && xright_j ==1){ //from left to right
    xnew[inputgraph->GetN()]   = xmax;    ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmax;    ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmin;    ynew[inputgraph->GetN()+2] = etamin;
  } else if (xleft_j ==1 && xright_i ==1){ //from right to left
    xnew[inputgraph->GetN()]   = xmin;    ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = xmin;    ynew[inputgraph->GetN()+1] = etamin;
    xnew[inputgraph->GetN()+2] = xmax;    ynew[inputgraph->GetN()+2] = etamin;	
  } else {
    // nominal
    xnew[inputgraph->GetN()]   = x_j; ynew[inputgraph->GetN()]   = y_j;
    xnew[inputgraph->GetN()+1] = x_j; ynew[inputgraph->GetN()+1] = y_j;
    xnew[inputgraph->GetN()+2] = x_j; ynew[inputgraph->GetN()+2] = y_j;
    
    // NP paper, achilleplot Bs
       /*
     if(x_j > 175.) {
       xnew[inputgraph->GetN()]   = xmax; ynew[inputgraph->GetN()]   = y_j;
       xnew[inputgraph->GetN()+1] = xmax; ynew[inputgraph->GetN()+1] = etamin;
       xnew[inputgraph->GetN()+2] = 0.; ynew[inputgraph->GetN()+2] = etamin;
     } else if(x_j < 10) {
       xnew[inputgraph->GetN()]   = xmin; ynew[inputgraph->GetN()]   = y_j;
       xnew[inputgraph->GetN()+1] = xmin; ynew[inputgraph->GetN()+1] = etamin;
       xnew[inputgraph->GetN()+2] = 0.; ynew[inputgraph->GetN()+2] = etamin;      
     } else {      
       xnew[inputgraph->GetN()]   = x_j; ynew[inputgraph->GetN()]   = y_j;
       xnew[inputgraph->GetN()+1] = x_j; ynew[inputgraph->GetN()+1] = y_j;
       xnew[inputgraph->GetN()+2] = x_j; ynew[inputgraph->GetN()+2] = y_j;
     }
       */
    // latticeQCD x SM paper
 
//     if(y_j > 0.) {
//       xnew[inputgraph->GetN()]   = x_j; ynew[inputgraph->GetN()]  = y_j;
//       xnew[inputgraph->GetN()+1] = 1.; ynew[inputgraph->GetN()+1] =  90.;
//       xnew[inputgraph->GetN()+2] = 1.; ynew[inputgraph->GetN()+2] = -90.;
//     } else {
//       xnew[inputgraph->GetN()]   = x_j; ynew[inputgraph->GetN()]  = y_j;
//       xnew[inputgraph->GetN()+1] = 1.; ynew[inputgraph->GetN()+1] = -90.;
//       xnew[inputgraph->GetN()+2] = 1.; ynew[inputgraph->GetN()+2] =  90.;
//     }
 
  }
  
  TGraph* newTGraph = new TGraph(inputgraph->GetN()+3,xnew,ynew);
  return newTGraph;
}
