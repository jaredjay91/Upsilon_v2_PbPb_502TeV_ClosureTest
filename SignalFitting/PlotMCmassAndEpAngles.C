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
#include "fileList.h"

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


void PlotMCmassAndEpAngles(
       int collId = kAADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=10,
       float muPtCut=3.5
			) 
{

  float muEtaCut = 2.4;
  float eta_low = -muEtaCut;
  float eta_high = muEtaCut;
  
  TGaxis::SetMaxDigits(3);
  gStyle->SetEndErrorSize(0);

  float dphiEp2Low = -1; 
  float dphiEp2High = 1;

  float ep2Low = -1; 
  float ep2High = 1;

  float ep2LowForPlot = ep2Low;    
  float ep2HighForPlot = ep2High;

  int   nEp2Bin  = (ep2High-ep2Low)*20;

  //TFile* fflat = new TFile("../OniaTree_Skim_UpsTrig_MC_flattened_190624.root");
  //TFile* fflat = new TFile("../OniaTree_Skim_UpsTrig_MC_flattened_190626.root");
  //TFile* fraw = new TFile("../OniaTree_Skim_UpsTrig_MC_190624.root"); 

  TFile* fflat = new TFile(flatSkimFileName0point5);
  TFile* fraw = new TFile(rawSkimFileName); 

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i",ptLow, ptHigh, yLow, yHigh, eta_high,eta_low, eta_high,eta_low, cLow*2,cHigh*2);

  float massLow = 8.0;
  float massHigh = 10.0;
  int   nMassBin  = (massHigh-massLow)*20;

  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)fflat->Get("dataset");
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

//  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
//  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = ws->var("mass")->frame(nMassBin); // bins
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
  drawText(Form("Centrality %i-%i%s", cLow,cHigh, perc.Data() ), pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);

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

  //RooPlot* myPlotEp = ws->var("ep2")->frame(nEp2Bin); // bins
  //ws->data("reducedDS")->plotOn(myPlotEp,Name("dataHist"));

  //Plot raw event plane as well
  RooDataSet *datasetraw = (RooDataSet*)fraw->Get("dataset");
  RooWorkspace *wsraw = new RooWorkspace("workspaceraw");
  wsraw->import(*datasetraw);
  cout << "####################################" << endl;
  RooDataSet *reducedDSraw = (RooDataSet*)datasetraw->reduce(RooArgSet(*(wsraw->var("mass")), *(wsraw->var("ep2")), *(wsraw->var("pt")), *(wsraw->var("y"))), kineCut.Data() );
  reducedDSraw->SetName("reducedDSraw");
  wsraw->import(*reducedDSraw);
  wsraw->var("ep2")->setRange(ep2Low*pi, ep2High*pi);
  wsraw->var("ep2")->Print();

  //RooPlot* myPlotEpraw = wsraw->var("ep2")->frame(nEp2Bin); // bins
  //wsraw->data("reducedDSraw")->plotOn(myPlotEpraw,Name("dataHistraw"));

  RooPlot* myPlotEpraw2 = wsraw->var("ep2")->frame(nEp2Bin); // bins
  wsraw->data("reducedDSraw")->plotOn(myPlotEpraw2,Name("dataOS_FITraw"),MarkerSize(.8));

  myPlotEpraw2->SetMinimum(0);
  //myPlotEpraw2->SetMaximum(6000);
  myPlotEpraw2->GetXaxis()->SetTitle("#psi_{ep}");
  myPlotEpraw2->GetXaxis()->SetTitleOffset(1.5);
  myPlotEpraw2->GetXaxis()->SetLabelOffset(0.04);
  myPlotEpraw2->GetXaxis()->SetLabelSize(0.05); //23
  myPlotEpraw2->GetXaxis()->SetTitleSize(0.05);  //28
  //myPlotEpraw2->GetXaxis()->CenterTitle();
  myPlotEpraw2->GetYaxis()->SetTitleOffset(1.5);
  myPlotEpraw2->GetYaxis()->SetTitle("Counts");
  myPlotEpraw2->GetYaxis()->SetTitleSize(0.05);
  myPlotEpraw2->GetYaxis()->SetLabelSize(0.05);
  myPlotEpraw2->GetYaxis()->CenterTitle();
  myPlotEpraw2->Draw();

  //Plot flattened event plane
  RooPlot* myPlotEp2 = ws->var("ep2")->frame(nEp2Bin); // bins
  ws->data("reducedDS")->plotOn(myPlotEp2,Name("dataOS_FIT"),MarkerSize(.8),MarkerColor(kRed));

  /*myPlotEp2->SetMinimum(0);
  //myPlotEp2->SetMaximum(6000);
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
  myPlotEp2->GetYaxis()->CenterTitle();*/
  myPlotEp2->Draw("same");

  drawText("CMS", CMS_x,CMS_y,text_color,CMS_text_size,CMS_text_font);
  drawText("Preliminary", CMS_x-CMS_text_size/500,CMS_y-Pre_text_size/500,text_color,Pre_text_size,Pre_text_font);

  c1->Update();



  TLegend fitleg = TLegend(0.66,0.4,0.8,0.6); fitleg.SetTextSize(19);
  fitleg.SetTextFont(43);
  fitleg.SetBorderSize(0);
  fitleg.AddEntry(myPlotEp2->findObject("dataOS_FIT"),"flattened","pe");
  fitleg.AddEntry(myPlotEpraw2->findObject("dataOS_FITraw"),"raw","pe");
  fitleg.Draw("same");

  c1->Update();

  TString plotTitle = Form("EpAng_pt%.1f-%.1f_y%.1f-%.1f_cent%i-%i",ptLow,ptHigh,yLow,yHigh,cLow,cHigh);
  c1->SaveAs(Form("FlatteningPlots/%s.pdf",plotTitle.Data()));
  c1->SaveAs(Form("FlatteningPlots/%s.png",plotTitle.Data()));

}
