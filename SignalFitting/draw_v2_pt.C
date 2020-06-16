#include "../headerFiles/SONGKYO.h"
#include "../headerFiles/tdrstyle.C"
#include "../headerFiles/CMS_lumi_square.C"
#include "../headerFiles/cutsAndBinUpsilonV2.h"
#include "fileList.h"



void draw_v2_pt(int whichUpsilon=1, float whichv2=0.5) {

  float ptbins[4] = {0,3,6,30};
  float ptpredval[3] = {0.05,0.05,0.05};
  float cbins[5] = {0,10,30,50,100};
  float cpredval[4] = {0.05,0.05,0.05,0.05};

  int numptbins = sizeof(ptbins)/sizeof(float)-1;
  int numcbins = sizeof(cbins)/sizeof(float)-1;

  setTDRStyle();
  writeExtraText = true;

  //Get v2 histograms
  TFile* inFile = new TFile(Form("Plots/Ups_%i_v2.root",whichUpsilon),"READ");
  TH1D* hv2pt = (TH1D*)inFile->Get("hv2pt;1");

  //Make it perfect with shrunk error bars
  bool makePerfect = kFALSE;
  if (makePerfect) {
    for (int i=0; i<numptbins; i++) {
      float binerr = hv2pt->GetBinError(i+1);
      hv2pt->SetBinContent(i+1,ptpredval[i]);
      hv2pt->SetBinError(i+1,binerr/sqrt(2));
    }
  }

  //Apply event plane resolution correction
  bool EpResCorrection = kFALSE;
  TFile* EPRFile = new TFile(resCorFileName,"READ");
  TH1D* hRpt = (TH1D*)EPRFile->Get("hRpt;1");
  hv2pt->Sumw2(); hRpt->Sumw2();
  if (EpResCorrection) {
    cout << " Applying resoulution correction" << endl;
    hv2pt->Divide(hRpt);
  }

  //make histograms of systematics
  TH1D* hsyspt = new TH1D("hsyspt","hsyspt",numptbins,ptbins);
  for (int i=0; i<numptbins; i++) {
    hsyspt->SetBinContent(i+1,0);//0.036 is the largest systematic uncertainty tabulated in the J/Psi analysis.
  }

  //make TGraphs
  TGraphErrors* gv2pt = new TGraphErrors(hv2pt);
  TGraphErrors* gv2pt_sys = new TGraphErrors(hv2pt);
  for (int ipt=0; ipt<numptbins; ipt++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2pt->GetPoint(ipt, pxtmp, pytmp);
    extmp=gv2pt->GetErrorX(ipt);
    eytmp=gv2pt->GetErrorY(ipt);
    relsys=hsyspt->GetBinContent(ipt+1);
    //remove ex
    gv2pt->SetPointError(ipt, 0, eytmp);
    //set ey for systematic error
    gv2pt_sys->SetPointError(ipt, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2pt, 1, 1); 
  SetGraphStyleSys(gv2pt_sys, 1);

  float ptmin=0; float ptmax=30;
  gv2pt_sys->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  gv2pt_sys->GetXaxis()->CenterTitle();
  gv2pt_sys->GetXaxis()->SetTitleOffset(1.);
  gv2pt_sys->GetXaxis()->SetLimits(ptmin,ptmax);
  gv2pt_sys->GetXaxis()->SetTitleSize(0.06);
  gv2pt_sys->GetYaxis()->SetTitle("v_{2}");
  gv2pt_sys->GetYaxis()->CenterTitle();
  gv2pt_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2pt_sys->GetYaxis()->SetTitleSize(0.06);
  gv2pt_sys->SetMinimum(-0.05);
  gv2pt_sys->SetMaximum(2*whichv2*1.5);
  if (whichUpsilon==3) {
    gv2pt_sys->SetMinimum(-5);
    gv2pt_sys->SetMaximum(200*whichv2);
  }
  if (makePerfect) {
    gv2pt_sys->SetMinimum(0);
    gv2pt_sys->SetMaximum(0.5*whichv2);
  }

  TCanvas* c2 = new TCanvas("c2","c2",40,40,600,600);

  gv2pt_sys->Draw("A5");
  gv2pt->Draw("P");

  TLegend *leg= new TLegend(0.7, 0.2, 0.9, 0.3);
  SetLegendStyle(leg);
  leg->AddEntry(gv2pt,Form(" #Upsilon(%iS)",whichUpsilon),"lp");

  leg->Draw("same");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);

  int collId = kAADATA;
  int cLow = 0;
  if(collId == kPPDATA) CMS_lumi_square(c2, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi_square(c2, 2 ,33);
  else if(collId == kPADATA) CMS_lumi_square(c2, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi_square(c2, 21 ,33);

  c2->SaveAs(Form("Plots/v2_vs_pt_%is.png",whichUpsilon));
  c2->SaveAs(Form("Plots/v2_vs_pt_%is.pdf",whichUpsilon));
}
