#include <iostream>
#include "../headerFiles/rootFitHeaders.h"
#include <RooWorkspace.h>
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "../headerFiles/cutsAndBinUpsilonV2.h"
#include "TCanvas.h"
#include "fileList.h"

const double pi = 3.14159265;

using namespace std;
using namespace RooFit;

bool isAbout(float num1=0.0, float num2=0.0) {

  float numavg = (num1+num2)/2;
  float percdiff = abs(num1-num2)/numavg;
  if (percdiff>0.01) return kFALSE;
  else return kTRUE;

}


void GetYieldsvsPhi( int collId = kAADATA,
  float ptLow = 0, float ptHigh = 30,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 100,
  float muPtCut = 3.5,
       float whichv2=0.5
 ) {

  //const int numphibins = 4;
  //float phibins[5] = {0.0, 0.125, 0.25, 0.375, 0.5};

  const int numphibins = 4;
  const int phiarraysize = numphibins+1;
  float phibins[phiarraysize] = {0};
  float binsize = 0.5*pi/((float)phiarraysize);
  for (int u=0; u<=phiarraysize; u++) {
    phibins[u] = ((float)u)*binsize;
  }

  float muEtaCut = 2.4;
  float eta_low = -muEtaCut;
  float eta_high = muEtaCut;
  
  TGaxis::SetMaxDigits(3);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);

  TString fileName;
  if (isAbout(whichv2,0.5)) fileName = flatSkimFileName0point5;
  else if (isAbout(whichv2,0.2)) fileName = flatSkimFileName0point2;
  else if (isAbout(whichv2,0.1)) fileName = flatSkimFileName0point1;
  else if (isAbout(whichv2,0.05)) fileName = flatSkimFileName0point05;
  else if (isAbout(whichv2,0.0)) fileName = flatSkimFileName0;

  TFile* funweighted = new TFile(fileName);

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i",ptLow, ptHigh, yLow, yHigh, eta_high,eta_low, eta_high,eta_low, cLow*2,cHigh*2);

  float massLow = 8.0;
  float massHigh = 10.0;
  int   nMassBin  = (massHigh-massLow)*20;

  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  //Create histogram
  TH1D* yieldsVsPhi = new TH1D("yieldsVsPhi", "yieldsVsPhi",numphibins,phibins);

  //Extract yields from MC file
  float dphiEp2Low, dphiEp2High;
  for (int iphi=0; iphi<numphibins; iphi++) {

    dphiEp2Low = phibins[iphi];
    dphiEp2High = phibins[iphi+1];

    //TString cutSpec = kineCut + Form(" && dphiEp2>%f && dphiEp2<%f", dphiEp2Low, dphiEp2High);
    TString cutSpec = kineCut + Form(" && abs(dphiEp2)>=%f && abs(dphiEp2)<%f", dphiEp2Low, dphiEp2High);

    //import and merge datasets
    RooDataSet *dataset = (RooDataSet*)funweighted->Get("dataset");
    RooWorkspace *ws = new RooWorkspace("workspace");
    ws->import(*dataset);
    cout << "####################################" << endl;
    RooDataSet *reducedDStmp = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("dphiEp2")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight"))), cutSpec.Data() );
    reducedDStmp->SetName("reducedDStmp");
    ws->import(*reducedDStmp);
    RooDataSet *reducedDS = new RooDataSet("reducedDS", "A sample", *reducedDStmp->get(), Import(*reducedDStmp), WeightVar(*ws->var("weight")) );
    reducedDS->SetName("reducedDS");
    ws->import(*reducedDS);

    float yield = reducedDS->sumEntries();
    //float yielderr = sqrt(yield);
    float yielderr = 0.052*yield;
    //float dyield = yield/binsize;
    //float dyielderr = yielderr/binsize;
    delete reducedDStmp;
    delete reducedDS;
    delete ws;

    cout << yield << " +/- " << yielderr << endl;
    //cout << dyield << " +/- " << dyielderr << endl;
    yieldsVsPhi->SetBinContent(iphi+1, yield);
    yieldsVsPhi->SetBinError(iphi+1, yielderr);
  }

  TCanvas* c1 = new TCanvas("c1","c1",50,50,550,550);
  yieldsVsPhi->Draw();
  yieldsVsPhi->SetMinimum(0);

  c1->SaveAs(Form("ExtractedYields/yieldsVsPhi_1S_%s.png",kineLabel.Data()));

  TFile* outFile = new TFile(Form("ExtractedYields/yieldsVsPhi_1S_%s.root",kineLabel.Data()),"RECREATE");
  yieldsVsPhi->Write();
  outFile->Close();

  delete yieldsVsPhi;
  delete c1;
  delete outFile;
  funweighted->Close("R");
  delete funweighted;

  cout << endl << "Here's what's in memory right now" << endl;
  gDirectory->ls("-m");
  cout << "that's all." << endl << endl;

}
