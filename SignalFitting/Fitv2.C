#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "../headerFiles/cutsAndBinUpsilonV2.h"

void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=18, float textFont=43){
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(textFont);
   //   if(bold)tex->SetTextFont(43);
   tex->SetTextSize(textSize);
   tex->SetTextColor(textColor);
   tex->SetLineWidth(1);
   tex->SetNDC();
   tex->Draw();
}

void Fitv2(

  int collId = kAADATA,
  float ptLow = 0, float ptHigh = 30,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 100,
  float muPtCut = 3.5,
  double* v2ValPtr=0, double* v2ErrPtr=0,
  int whichUpsilon=1 )

{

  float muEtaCut = 2.4;

  gStyle->SetOptStat(0);

  TString fileLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //Load the yields vs phi
  TFile* yieldsFile = TFile::Open(Form("ExtractedYields/yieldsVsPhi_%iS_%s.root",whichUpsilon,fileLabel.Data()),"READ");
  TH1D* yieldsVsPhi = (TH1D*)yieldsFile->Get("yieldsVsPhi;1");
  yieldsVsPhi->Sumw2();
  float binwidth = yieldsVsPhi->GetBinWidth(1);
  cout << "binwidth = " << binwidth << endl;

  //Plot it
  TCanvas* c1 = new TCanvas("c1","c1",50,50,550,550);
  yieldsVsPhi->Draw();
  yieldsVsPhi->GetXaxis()->SetTitle("|#Delta#phi|/#pi");
  double totalintegral = yieldsVsPhi->Integral(1,4);
  yieldsVsPhi->Scale(1.0/totalintegral);

  //Fit it
  /*TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265) + 2*[2]*cos(3*x*3.14159265) + 2*[3]*cos(4*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2","v3","v4");*/
  TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2");
  /*TF1* fitfunc = new TF1("fitfunc"," 1 + 2*[1]*cos(2*x*3.14159265)",0,0.5);
  fitfunc->SetParNames("v2");*/
  yieldsVsPhi->Fit("fitfunc");

  double v2Val = fitfunc->GetParameter(1)/binwidth;
  double v2Err = fitfunc->GetParError(1)/binwidth;
  cout << "v2 = " << v2Val << " +/- " << v2Err << endl;
  TLatex latex;
  latex.SetTextSize(0.05);
  //latex.SetTextAlign(13);
  latex.DrawLatex(0.25,50,Form("v_{2} = %.3f #pm %.3f",v2Val,v2Err));

  TString perc = "%";
  float pos_text_x = 0.15;
  float pos_text_y = 0.45;
  float pos_y_diff = 0.06;
  float text_size = 16;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < %.1f GeV/c", muEtaCut ), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %i-%i%s", cLow,cHigh, perc.Data() ), pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
  drawText(Form("v_{2} = %.3f #pm %.3f", v2Val,v2Err), pos_text_x,pos_text_y-pos_y_diff*5,2,text_size);

  c1->Update();

  c1->SaveAs(Form("Plots/v2_fit_pt%.1f-%.1f_y%.2f-%.2f_cent%i-%i.png",ptLow,ptHigh, yLow,yHigh, cLow, cHigh));

  *v2ValPtr = v2Val;
  *v2ErrPtr = v2Err;
  
  yieldsFile->Close();
  delete fitfunc;
  delete c1;
}
