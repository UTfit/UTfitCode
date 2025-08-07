{

//TCanvas *c1 = new TCanvas("c1","Stat. Sig",200,10,700,500);
c1 = new TCanvas("his1D","1D histograms",2);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.05);
c1->SetBottomMargin(0.15);
c1->SetFillColor(0);
c1->SetHighLightColor(0);
c1->SetBorderMode(0);
gPad->SetTickx(0);
gStyle->SetTickLength(0,"x");
gStyle->SetOptStat(0);
gStyle->SetLabelSize(0.1,"x");
gStyle->SetLabelOffset(.016,"x");
gStyle->SetTitleOffset(1.5,"X");
gStyle->SetTitleSize(0.05,"X");
gStyle->SetTitleSize(0.05,"Y");
gStyle->SetTitleOffset(1.1,"Y");
//gStyle->SetLabelSize(0.06,"y");

int color[5];
 color[0] =  416+3; // kGreen+3;
 color[1] =  416+1; // kGreen+1;
 //color[2] =  400; //kYellow
 color[2] =  800-3; // kOrange -3
 color[3] =  632; //kRed
 color[4] =  600; //kBlue

//--------- Insert here the vectors out of scales.py ------------
Double_t ReCK[5] = {0.4, 3.1, 1.0, 6.0, 1.3};
Double_t ImCK[5] = {7.5, 55.8, 17.4, 108.5, 23.7};
Double_t ReCD[5] = {0.6, 1.2, 0.4, 2.2, 0.7};
Double_t ImCD[5] = {3.6, 6.9, 2.0, 12.3, 4.2};
Double_t AbsCBd[5] = {8.7, 14.0, 7.2, 22.9, 14.0};
Double_t AbsCBs[5] = {10.3, 17.1, 8.8, 27.6, 17.2};
//---------------------------------------------------------------

double CiVal[5][5];
double ind=0.59;
for(int i=0;i<5;i++) {
  for(int j=0;j<5;j++) {
    //cout << i << " "<< j << " "<< index+j*0.90+i*0.125 << endl;
    CiVal[i][j] = ind+j*0.90+i*0.125;
  }
}

double Scale[5][5];
for(int j=0;j<5;j++) {
  Scale[0][j] = ReCK[j]; 
  Scale[1][j] = ImCK[j]; 
  //Scale[2][j] = ReCD[j]; 
  Scale[2][j] = ImCD[j]; 
  Scale[3][j] = AbsCBd[j]; 
  Scale[4][j] = AbsCBs[j]; 
}

TH1D* gr[5];
gr[0] = new TH1D("gr0","",36,0.45,4.95);
gr[1] = new TH1D("gr1","",36,0.45,4.95);
//gr[2] = new TH1D("gr2","",36,0.45,4.95);
gr[2] = new TH1D("gr3","",36,0.45,4.95);
gr[3] = new TH1D("gr4","",36,0.45,4.95);
gr[4] = new TH1D("gr5","",36,0.45,4.95);

for(int i=0;i<5;i++) {
  for(int j=0;j<5;j++) {
    gr[i]->SetBinContent(gr[i]->FindBin(CiVal[i][j]),Scale[i][j]);
    //cout << "sto mettendo: " << Scale[i][j] << " for " << CiVal[i][j] << " nel bin " << i << "," << j <<  endl;
  }
}

gr[1]->GetXaxis()->SetTitle("");
gr[1]->GetYaxis()->SetTitle("NP scale #Lambda (TeV)");

const char *wil[5]  = {"C_{1}","C_{2}","C_{3}","C_{4}","C_{5}"};
for (int j=1;j<=5;j++) gr[1]->GetXaxis()->SetBinLabel(gr[1]->FindBin(CiVal[1][j-1]),wil[j-1]);

gr[1]->SetMarkerColor(color[1]);
gr[1]->SetFillColor(color[1]);
gr[1]->SetMarkerStyle(21);
gr[1]->SetMarkerSize(2.);
gr[1]->SetMinimum(0.1);
gr[1]->SetMaximum(550);
gr[1]->Draw();

gr[0]->SetMarkerColor(color[0]);
gr[0]->SetFillColor(color[0]);
gr[0]->SetMarkerStyle(21);
gr[0]->SetMarkerSize(2.);
gr[0]->Draw("same");

//gr[2]->SetMarkerColor(color[2]);
//gr[2]->SetFillColor(color[2]);
//gr[2]->SetMarkerStyle(21);
//gr[2]->SetMarkerSize(2.);
//gr[2]->Draw("same");

gr[2]->SetMarkerColor(color[2]);
gr[2]->SetFillColor(color[2]);
gr[2]->SetMarkerStyle(21);
gr[2]->SetMarkerSize(2.);
gr[2]->Draw("same");

gr[3]->SetMarkerColor(color[3]);
gr[3]->SetFillColor(color[3]);
gr[3]->SetMarkerStyle(21);
gr[3]->SetMarkerSize(2.);
gr[3]->Draw("same");

gr[4]->SetMarkerColor(color[4]);
gr[4]->SetFillColor(color[4]);
gr[4]->SetMarkerStyle(21);
gr[4]->SetMarkerSize(2.);
gr[4]->Draw("same");

TLegend *legend=new TLegend(0.15,0.58,0.30,0.88);
legend->SetTextFont(72);
legend->SetTextSize(0.04);
// char *wilson[6]  = {"Re C_{K}","Im C_{K}","Re C_{D}","Im C_{D}","C_{Bd}","C_{Bs}"};
const char *wilson[5]  = {"Re C_{K}","Im C_{K}","Im C_{D}","C_{Bd}","C_{Bs}"};
for (int i=0; i<5; i++) {
  legend->AddEntry(gr[i],wilson[i],"p");
}
legend->Draw();

c1->SetLogy();

TASImage *img = new TASImage("../common/logo.png");
TPad *pad = new TPad("pad","pad",0.74,0.75,0.94,0.88,0,0,0);
img->SetConstRatio(kFALSE);
img->SetImageQuality(TAttImage::kImgBest);
c1->cd();
pad->Draw();
pad->cd();
img->Draw();

TPad *padSeason = new TPad("padSeason","padSeason",0.74,0.65,0.94,0.78,0,0,0);
padSeason->SetFillStyle(4000);
c1->cd();
padSeason->Draw();
padSeason->cd();

TLatex *tex = new TLatex;
//tex->SetTextSize(0.03);
tex->SetTextSize(0.28);
tex->SetTextAlign(21);
tex->SetTextColor(602);
tex->DrawTextNDC(0.50,0.63, "summer23");

/*
TLatex *tex = new TLatex;
//tex->SetTextSize(0.03);
tex->SetTextSize(0.035);
tex->SetTextAlign(21);
tex->SetTextColor(602);
tex->DrawTextNDC(0.85,0.72, "summer21");
*/
 
c1->SaveAs("df2-npscale-nmfv.eps");
}
