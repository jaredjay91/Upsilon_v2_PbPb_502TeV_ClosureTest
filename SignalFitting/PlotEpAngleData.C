//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "../rootFitHeaders.h"
#include "../commonUtility.h"
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
void PlotEpAngleData( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,//Run 1 has p going in -z direction
       int cLow=20, int cHigh=120,
       float muPtCut=4.0
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

  int   nEp2Bin  = (ep2High-ep2Low)*10;


  TFile* f1 = new TFile("OniaTree_Skim_UpsTrig_MC_190619.root");
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i",ptLow, ptHigh, yLow, yHigh, eta_high,eta_low, eta_high,eta_low, cLow,cHigh);


  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ep2")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("ep2")->setRange(ep2Low*pi, ep2High*pi);
  ws->var("ep2")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  
  RooPlot* myPlot = ws->var("ep2")->frame(nEp2Bin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(ep2LowForPlot, ep2HighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.04);
  myPlot2->GetXaxis()->SetLabelSize(0.04);
  myPlot2->GetXaxis()->SetTitleSize(0.05);
  myPlot2->GetXaxis()->SetTitle("#psi");
  myPlot2->GetXaxis()->SetRangeUser(ep2Low*pi, ep2High*pi);
  myPlot2->Draw();
  TString perc = "%";


  c1->Update();

  c1->SaveAs("Ep2Plot.png");

  cout << "if ( binMatched( "<<muPtCut<<",  " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ; 

} 
 
