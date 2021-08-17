#include "../headerFiles/SONGKYO.h"
#include "../headerFiles/tdrstyle.C"
#include "../headerFiles/CMS_lumi_square.C"
#include "../headerFiles/cutsAndBinUpsilonV2.h"
#include "fileList.h"



void draw_v2_centrality(int whichUpsilon=1, float whichv2=0.5) {

  float ptbins[4] = {0,3,6,30};
  float ptpredval[3] = {0.05,0.05,0.05};
  float cbins[5] = {0,10,30,50,100};
  float cpredval[4] = {0.05,0.05,0.05,0.05};

  int numptbins = sizeof(ptbins)/sizeof(float)-1;
  int numcbins = sizeof(cbins)/sizeof(float)-1;

  setTDRStyle();
  writeExtraText = true;
  extraText = "Simulation";

  //Get v2 histograms
  TFile* inFile = new TFile(Form("Plots/Ups_%i_v2.root",whichUpsilon),"READ");
  TH1D* hv2c = (TH1D*)inFile->Get("hv2c;1");

  //Make it perfect with shrunk error bars
  bool makePerfect = kFALSE;
  if (makePerfect) {
    for (int i=0; i<numcbins; i++) {
      float binerr = hv2c->GetBinError(i+1);
      hv2c->SetBinContent(i+1,cpredval[i]);
      hv2c->SetBinError(i+1,binerr/sqrt(2));
    }
  }

  //Apply event plane resolution correction
  bool EpResCorrection = kFALSE;
  TFile* EPRFile = new TFile(resCorFileName,"READ");
  TH1D* hRcent = (TH1D*)EPRFile->Get("hRcent;1");
  hv2c->Sumw2(); hRcent->Sumw2();
  if (EpResCorrection) {
    cout << " Applying resoulution correction" << endl;
    hv2c->Divide(hRcent);
  }

  //make histograms of systematics
  TH1D* hsysc = new TH1D("hsysc","hsysc",numcbins,cbins);
  for (int i=0; i<numcbins; i++) {
    hsysc->SetBinContent(i+1,0);//0.036 is the largest systematic uncertainty tabulated in the J/Psi analysis.
    //print out percent difference from true v2.
    float recov2 = hv2c->GetBinContent(i+1);
    cout << "Difference in v2 = " << recov2-whichv2 << endl;
    cout << "Percent difference in v2 = " << fabs(recov2-whichv2)/whichv2*100 << "%" << endl;
  }

  //make TGraphs
  TGraphErrors* gv2c = new TGraphErrors(hv2c);
  TGraphErrors* gv2c_sys = new TGraphErrors(hv2c);
  for (int ic=0; ic<numcbins; ic++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2c->GetPoint(ic, pxtmp, pytmp);
    extmp=gv2c->GetErrorX(ic);
    eytmp=gv2c->GetErrorY(ic);
    relsys=hsysc->GetBinContent(ic+1);
    //remove ex
    gv2c->SetPointError(ic, 0, eytmp);
    //set ey for systematic error
    gv2c_sys->SetPointError(ic, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2c, 1, 1); 
  SetGraphStyleSys(gv2c_sys, 1);

  float cmin=0; float cmax=100;
  gv2c_sys->GetXaxis()->SetTitle("Centrality (%)");
  gv2c_sys->GetXaxis()->CenterTitle();
  gv2c_sys->GetXaxis()->SetTitleOffset(1.);
  gv2c_sys->GetXaxis()->SetLimits(cmin,cmax);
  gv2c_sys->GetXaxis()->SetTitleSize(0.06);
  gv2c_sys->GetYaxis()->SetTitle("v_{2}");
  gv2c_sys->GetYaxis()->CenterTitle();
  gv2c_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2c_sys->GetYaxis()->SetTitleSize(0.06);
  gv2c_sys->SetMinimum(-0.05);
  gv2c_sys->SetMaximum(2*whichv2*1.5);
  if (whichUpsilon==3) {
    gv2c_sys->SetMinimum(-5);
    gv2c_sys->SetMaximum(200*whichv2);
  }
  if (makePerfect) {
    gv2c_sys->SetMinimum(0);
    gv2c_sys->SetMaximum(0.5*whichv2);
  }

  TCanvas* c2 = new TCanvas("c2","c2",40,40,600,600);

  gv2c_sys->Draw("A5");
  gv2c->Draw("P");

  TLine* l1 = new TLine(0,whichv2,100,whichv2);
  l1->SetLineStyle(kDashed);
  l1->Draw("same");

  TLegend *leg= new TLegend(0.7, 0.2, 0.9, 0.3);
  SetLegendStyle(leg);
  leg->AddEntry(gv2c,Form(" #Upsilon(%iS)",whichUpsilon),"lp");

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

  c2->SaveAs(Form("Plots/v2_vs_c_%is.png",whichUpsilon));
  c2->SaveAs(Form("Plots/v2_vs_c_%is.pdf",whichUpsilon));
}
