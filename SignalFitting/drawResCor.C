#include "fileList.h"

void drawResCor() {

  gStyle->SetOptStat(0);

  TFile* inFile = TFile::Open(resCorFileName,"READ");
  TH1D* hRpt = (TH1D*)inFile->Get("hRpt;1");
  TH1D* hRcent = (TH1D*)inFile->Get("hRcent;1");

  hRpt->GetXaxis()->SetTitle("p_T");
  hRpt->SetMinimum(0.4);
  hRpt->SetMaximum(1.0);
  //hRpt->GetXaxis()->SetTitleSize(0.06);
  hRpt->GetYaxis()->SetTitleSize(0.06);

  hRcent->GetXaxis()->SetTitle("centrality");
  hRcent->SetMinimum(0.4);
  hRcent->SetMaximum(1.0);
  //hRcent->GetXaxis()->SetTitleSize(0.06);
  hRcent->GetYaxis()->SetTitleSize(0.06);

  hRpt->Sumw2();
  hRcent->Sumw2();

  TCanvas* cpt = new TCanvas("cpt","cpt",0,0,500,500);
  hRpt->Draw();

  TCanvas* ccent = new TCanvas("ccent","ccent",0,0,500,500);
  hRcent->Draw();

  cpt->SaveAs("cpt.pdf");
  ccent->SaveAs("ccent.pdf");
}
