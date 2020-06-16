#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "Fitv2.C"


void Get_v2_vs_var(int whichUpsilon=1) {

  gStyle->SetOptStat(0);

  float ptbins[4] = {0,3,6,30};
  float cbins[5] = {0,10,30,50,100};

  int numptbins = sizeof(ptbins)/sizeof(float)-1;
  int numcbins = sizeof(cbins)/sizeof(float)-1;

  int collId = kAADATA;
  float muPtCut = 3.5;

  double v2Val; double v2Err;
  double* v2ValPtr; double* v2ErrPtr;
  v2ValPtr = &v2Val; v2ErrPtr = &v2Err;

  float ptLow; float ptHigh;
  int cLow = 0; int cHigh = 100;
  float yLow = 0.0;
  float yHigh = 2.4;

  TH1D* hv2pt = new TH1D("hv2pt","hv2pt",numptbins,ptbins);
  TH1D* hv2c = new TH1D("hv2c","hv2c",numcbins,cbins);

  //PT BINS
  for (int ipt=0; ipt<numptbins; ipt++) {
    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    cout << "[" << ptLow << "," << ptHigh << "]" << endl;
    //cLow = 10; cHigh = 60;

    Fitv2(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, v2ValPtr, v2ErrPtr, whichUpsilon);
    hv2pt->SetBinContent(ipt+1,v2Val);
    hv2pt->SetBinError(ipt+1,v2Err);
  }
  
  TCanvas* cv2 = new TCanvas("cv2","cv2",800,400);
  cv2->Divide(2,1);
  cv2->cd(1);
  hv2pt->Draw();

  //Centrality BINS
  for (int ic=0; ic<numcbins; ic++) {
    cLow = cbins[ic];
    cHigh = cbins[ic+1];
    cout << "[" << cLow << "," << cHigh << "]" << endl;
    ptLow = 0; ptHigh = 30;

    Fitv2(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, v2ValPtr, v2ErrPtr, whichUpsilon);
    hv2c->SetBinContent(ic+1,v2Val);
    hv2c->SetBinError(ic+1,v2Err);
  }
  cv2->cd(2);
  hv2c->Draw();

  TString outFileName = Form("Plots/Ups_%i_v2.root",whichUpsilon);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hv2pt->Write();
  hv2c->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
