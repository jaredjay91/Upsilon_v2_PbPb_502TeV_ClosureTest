#include <iostream>
#include "../rootFitHeaders.h"
//#include "../commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../PsetCollection.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"
#include "../Style_jaebeom.h"

const double pi = 3.14159265;

using namespace std;
using namespace RooFit;

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

void PlotMCmassAndEpAngle(
       int collId = kAADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,//Run 1 has p going in -z direction
       int cLow=20, int cHigh=120,
       float muPtCut=3.5
			) 
{

  float muEtaCut = 2.4;
  float eta_low = -muEtaCut;
  float eta_high = muEtaCut;
  
  gStyle->SetEndErrorSize(0);

  float dphiEp2Low = -1; 
  float dphiEp2High = 1;

  float ep2Low = -1; 
  float ep2High = 1;

  float ep2LowForPlot = ep2Low;    
  float ep2HighForPlot = ep2High;

  int   nEp2Bin  = (ep2High-ep2Low)*20;

  TFile* f1;
  bool flatten = kFALSE;
  TString flatString = "";
  if (flatten) {
    f1 = new TFile("../OniaTree_Skim_UpsTrig_MC_flattened_190712_v2_0point5.root");
    flatString = "_flattened";
  }
  else {
    f1 = new TFile("../OniaTree_Skim_UpsTrig_MC_190624.root"); 
  }
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i",ptLow, ptHigh, yLow, yHigh, eta_high,eta_low, eta_high,eta_low, cLow,cHigh);

  float massLow = 8.0;
  float massHigh = 10.0;
  int   nMassBin  = (massHigh-massLow)*20;

  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("ep2")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
  ws->var("ep2")->setRange(ep2Low*pi, ep2High*pi);
  ws->var("ep2")->Print();

  TCanvas* c1 = new TCanvas("c1","c1",0,0,1150,550);
  c1->Divide(2);
  c1->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  myPlot2->SetMinimum(0);
  myPlot2->SetMarkerStyle(kFullCircle);
  myPlot2->SetTitleSize(0.0);
  myPlot2->SetTitle("PbPb 1.6nb^{-1} (5.02 TeV)");
  myPlot2->SetTitleOffset(5.0);

  myPlot2->GetYaxis()->SetTitleOffset(1.5);
  myPlot2->GetYaxis()->SetTitle("Counts");
  myPlot2->GetYaxis()->SetTitleSize(0.05);
  myPlot2->GetYaxis()->SetLabelSize(0.05);
  myPlot2->GetYaxis()->CenterTitle();

  myPlot2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2->GetXaxis()->SetTitleOffset(1.5);
  myPlot2->GetXaxis()->SetLabelOffset(0.04);
  myPlot2->GetXaxis()->SetLabelSize(0.05); //23
  myPlot2->GetXaxis()->SetTitleSize(0.05);  //28
  myPlot2->GetXaxis()->CenterTitle();

  myPlot2->GetYaxis()->SetTickSize(0.04);
  myPlot2->GetYaxis()->SetNdivisions(404);
  myPlot2->GetXaxis()->SetTickSize(0.03);
  myPlot2->Draw();


  TString perc = "%";
  float pos_text_x = 0.2;
  float pos_text_y = 0.7;
  float pos_y_diff = 0.06;
  float text_size = 16;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < %.1f GeV/c", muEtaCut ), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %i-%i%s", cLow/2,cHigh/2, perc.Data() ), pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);

  float CMS_text_size = 30;
  float CMS_text_font = 63;
  float CMS_y = 0.83;
  float CMS_x = 0.78;
  drawText("CMS", CMS_x,CMS_y,text_color,CMS_text_size,CMS_text_font);
  float Pre_text_size = 20;
  float Pre_text_font = 53;
  drawText("Preliminary", CMS_x-CMS_text_size/500,CMS_y-Pre_text_size/500,text_color,Pre_text_size,Pre_text_font);

  c1->Update();

  //Plot something else
  c1->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  RooPlot* myPlotEp = ws->var("ep2")->frame(nEp2Bin); // bins
  ws->data("reducedDS")->plotOn(myPlotEp,Name("dataHist"));

  RooPlot* myPlotEp2 = (RooPlot*)myPlotEp->Clone();
  ws->data("reducedDS")->plotOn(myPlotEp2,Name("dataOS_FIT"),MarkerSize(.8));

  myPlotEp2->SetMinimum(0);
  myPlotEp2->SetMaximum(12000);
  myPlotEp2->GetXaxis()->SetTitle("#psi_{ep}");
  myPlotEp2->GetXaxis()->SetTitleOffset(1.5);
  myPlotEp2->GetXaxis()->SetLabelOffset(0.04);
  myPlotEp2->GetXaxis()->SetLabelSize(0.05); //23
  myPlotEp2->GetXaxis()->SetTitleSize(0.05);  //28
  //myPlotEp2->GetXaxis()->CenterTitle();
  myPlotEp2->GetYaxis()->SetTitleOffset(1.5);
  myPlotEp2->GetYaxis()->SetTitle("Counts");
  myPlotEp2->GetYaxis()->SetTitleSize(0.05);
  myPlotEp2->GetYaxis()->SetLabelSize(0.05);
  myPlotEp2->GetYaxis()->CenterTitle();
  myPlotEp2->Draw();

  drawText("CMS", CMS_x,CMS_y,text_color,CMS_text_size,CMS_text_font);
  drawText("Preliminary", CMS_x-CMS_text_size/500,CMS_y-Pre_text_size/500,text_color,Pre_text_size,Pre_text_font);

  c1->Update();

  TString plotTitle = Form("EpAng_pt%.1f-%.1f_y%.1f-%.1f_cent%i-%i%s",ptLow,ptHigh,yLow,yHigh,cLow/2,cHigh/2,flatString.Data());
  c1->SaveAs(Form("Plots/%s.pdf",plotTitle.Data()));
  c1->SaveAs(Form("Plots/%s.png",plotTitle.Data()));

}
