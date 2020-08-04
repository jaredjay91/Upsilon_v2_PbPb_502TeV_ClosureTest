
#include <iostream>
#include "../headerFiles/rootFitHeaders.h"
#include "../headerFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../headerFiles/cutsAndBinUpsilonV2.h"
#include "../headerFiles/PsetCollection.h"
#include "../headerFiles/CMS_lumi_square.C"
#include "../headerFiles/tdrstyle.C"
#include "../headerFiles/Style_jaebeom.h"
#include "fileList.h"

const static float pi = 3.14159265;

using namespace std;
using namespace RooFit;

bool isAbout(float num1=0.0, float num2=0.0) {

  float numavg = (num1+num2)/2;
  float percdiff = abs(num1-num2)/numavg;
  if (percdiff>0.01) return kFALSE;
  else return kTRUE;

}

void PlotDataset( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,
       int cLow=0, int cHigh=200,
       float muPtCut=3.5,
       float whichv2=0.5,
       int whichVar=0
//RooDataSet::dataset[mass, pt, y, pt1, pt2, eta1, eta2, weight, cBin, ep2, recoQQsign, dphiEp2, qx, qy] = 4023601 entries
			) 
{

  TString vars[14] = {"mass", "pt", "y", "pt1", "pt2", "eta1", "eta2", "weight", "cBin", "ep2", "recoQQsign", "dphiEp2", "qx", "qy"};
  float varLow[14]  = { 8,  0, -3,  0,  0, -3, -3, 0,   0, -pi, -1, -pi, -50, -50};
  float varHigh[14] = {14, 30,  3, 30, 30,  3,  3, 2, 200,  pi,  1,  pi,  50,  50};

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

  int   nVarBin  = (varHigh[whichVar]-varLow[whichVar])*10;

  TString fileName;
  if (isAbout(whichv2,0.5)) fileName = flatSkimFileName0point5;
  else if (isAbout(whichv2,0.2)) fileName = flatSkimFileName0point2;
  else if (isAbout(whichv2,0.1)) fileName = flatSkimFileName0point1;
  else if (isAbout(whichv2,0.05)) fileName = flatSkimFileName0point05;

  TFile* f1 = new TFile(fileName);

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i",ptLow, ptHigh, yLow, yHigh, eta_high,eta_low, eta_high,eta_low, cLow,cHigh);

  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;

  //unweighted
  /*RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var(vars[whichVar])), *(ws->var("weight")), kineCut.Data() ));
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);*/

  //weighted
  RooDataSet *reducedDStmp = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var(vars[whichVar])), *(ws->var("weight"))), kineCut.Data() );
  reducedDStmp->SetName("reducedDStmp");
  ws->import(*reducedDStmp);
  RooDataSet *reducedDS = new RooDataSet("reducedDS", "A sample", *reducedDStmp->get(), Import(*reducedDStmp), WeightVar(*ws->var("weight")) );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);

  ws->var(vars[whichVar])->setRange(varLow[whichVar], varHigh[whichVar]);
  ws->var(vars[whichVar])->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  
  RooPlot* myPlot = ws->var(vars[whichVar])->frame(nVarBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(varLow[whichVar], varHigh[whichVar],"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.04);
  myPlot2->GetXaxis()->SetLabelSize(0.04);
  myPlot2->GetXaxis()->SetTitleSize(0.05);
  myPlot2->GetXaxis()->SetTitle(vars[whichVar]);
  myPlot2->GetXaxis()->SetRangeUser(varLow[whichVar], varHigh[whichVar]);
  myPlot2->Draw();
  TString perc = "%";

  c1->Update();

  c1->SaveAs(Form("DatasetPlots_weighted/%sPlot.png",vars[whichVar].Data()));

  cout << "if ( binMatched( "<<muPtCut<<",  " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ; 

} 
 
