#include "GetYieldsvsPhi.C"

const int numybins = 11;
float ybins[14] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 0.0, 0.8, 1.6, 2.4, 0.0, 1.2, 2.4};

float ptbins[4] = {0,3,6,30};
const int numptbins = sizeof(ptbins)/sizeof(float)-1;

int cbins[5] = {0,10,30,50,100};
const int numcbins = sizeof(cbins)/sizeof(int)-1;

void GetAllYieldsvsPhi(float whichv2=0.5) {

  int collId = kAADATA;
  int whichModel = 0;
  float muPtCut=3.5;
  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;

  bool INTBIN = kFALSE;
  bool PTBINS = kFALSE;
  bool YBINS = kFALSE;
  bool CBINS = kTRUE;

  //integrated bin
  if (INTBIN) GetYieldsvsPhi(collId, 0, 30, 0, 2.4, 0, 100, muPtCut, whichv2);

  //Pt bins
  if (PTBINS) {
    cout << endl << "********** STARTING PT BINS **********" << endl;
    yLow = 0.0;
    yHigh = 2.4;
    cLow = 0;
    cHigh = 100;
    for (int ipt=1; ipt<numptbins; ipt++) {
      ptLow = ptbins[ipt];
      ptHigh = ptbins[ipt+1];
      cout << endl << "[" << ptLow << "," << ptHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichv2);
    }
  }

  //Rapidity bins
  if (YBINS) {
    ptLow = 0;
    ptHigh = 30;
    cLow = 0;
    cHigh = 100;
    for (int iy=0; iy<numybins; iy++) {
      yLow = ybins[iy];
      yHigh = ybins[iy+1];
      cout << endl << "[" << yLow << "," << yHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichv2);
    }
  }
 
  //Centrality bins
  if (CBINS) {
    cout << endl << "********** STARTING CENTRALITY BINS **********" << endl;
    ptLow = 0;
    ptHigh = 30;
    yLow = 0.0;
    yHigh = 2.4;
    for (int ic=0; ic<numcbins; ic++) {
      cLow = cbins[ic];
      cHigh = cbins[ic+1];
      cout << endl << "[" << cLow << "," << cHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichv2);
    }
  }
}



