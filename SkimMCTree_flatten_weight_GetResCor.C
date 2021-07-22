#include <ctime>

#include <TLorentzVector.h>
#include "headerFiles/commonUtility.h"
#include "headerFiles/HiEvtPlaneList.h"
#include "headerFiles/cutsAndBinUpsilonV2.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"

const Double_t pi = 3.141592653589;

Double_t getDPHI_Jared( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
 
  if ( dphi > pi ) dphi = dphi - 2*pi;
  if ( dphi <= -pi ) dphi = dphi + 2*pi;
  if ( fabs(dphi) > pi ) {
    //cout << " getDPHI error!!! dphi=" << phi1 << "-" << phi2 << " is bigger than 3.141592653589 " << endl;
    dphi = -10;
  }

  return dphi;
}

Double_t restrictToPiOver2(Double_t phi) {

  Double_t newphi = phi;
  if (phi>pi/2) newphi = phi-pi;
  else if (phi<-pi/2) newphi = phi+pi;
  return newphi;
}

bool isAbout(float num1=0.0, float num2=0.0) {

  float numavg = (num1+num2)/2;
  float percdiff = abs(num1-num2)/numavg;
  if (percdiff>0.01) return kFALSE;
  else return kTRUE;

}

static const long MAXTREESIZE = 1000000000000;

void SkimMCTree_flatten_weight_GetResCor(int nevt=-1,
      int dateStr=20210721,
      bool flattenBinByBin=kTRUE,
      Double_t v2 = 0.5)
{

  Double_t v2_present_pt[3] = {0.00311049, 0.00856134, 0.0298847};
  Double_t v2_present_c[4] = {0.00698024, 0.000933421, 0.00656082, 0.0387734};
  Double_t v2_present_average = 0.022685505;

  const int flatOrder = 21;

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  Double_t ptbins[4] = {0,3,6,30};
  const int numptbins = sizeof(ptbins)/sizeof(double)-1;

  Double_t cbins[5] = {0,10,30,50,100};
  const int numcbins = sizeof(cbins)/sizeof(double)-1;

  //TH1D* hCosAC = new TH1D("hCosAC","cos(2*(psiA-psiC))",50,-1.2,1.2);
  TH1D* hRpt = new TH1D("hRpt","EP Resolution factor vs pt",numptbins,ptbins);
  TH1D* hRcent = new TH1D("hRcent","EP Resolution factor vs cent",numcbins,cbins);

  //Event planes
  TH1D* hEpHF2old = new TH1D("hEpHF2old","hEpHF2old",50,-2,2);
  TH1D* hEpHF2new = new TH1D("hEpHF2new","hEpHF2new",50,-2,2);
  TH1D* hEpHFm2old = new TH1D("hEpHFm2old","hEpHFm2old",50,-2,2);
  TH1D* hEpHFm2new = new TH1D("hEpHFm2new","hEpHFm2new",50,-2,2);
  TH1D* hEpHFp2old = new TH1D("hEpHFp2old","hEpHFp2old",50,-2,2);
  TH1D* hEpHFp2new = new TH1D("hEpHFp2new","hEpHFp2new",50,-2,2);
  TH1D* hEptrackmid2old = new TH1D("hEptrackmid2old","hEptrackmid2old",50,-2,2);
  TH1D* hEptrackmid2new = new TH1D("hEptrackmid2new","hEptrackmid2new",50,-2,2);

  TH1D* hEpHF2 = new TH1D("hEpHF2","hEpHF2",50,-2,2);
  TH1D* hEpHFm2 = new TH1D("hEpHFm2","hEpHFm2",50,-2,2);
  TH1D* hEpHFp2 = new TH1D("hEpHFp2","hEpHFp2",50,-2,2);
  TH1D* hEptrackmid2 = new TH1D("hEptrackmid2","hEptrackmid2",50,-2,2);

  TH1D* hphi = new TH1D("hphi","hphi",50,-2,2);
  TH1D* hdphi = new TH1D("hdphi","hdphi",50,-2,2);
  TH1D* hdphiw = new TH1D("hdphiw","hdphiw",50,-2,2);

  //TString fnameData1 = "../Oniatree_Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_MC_190610.root";
  //TFile* MCfile = TFile::Open(fnameData1,"READ");
  //TTree *mytree = (TTree*)MCfile->Get("myTree");
  //TTree* tree = (TTree*)MCfile->Get("tree");
  //mytree->AddFriend(tree);

  TString inputMC1 = "../DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_v1.root";
  TString inputMC2 = "../DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_ext-v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC1.Data());
  mytree->Add(inputMC2.Data());

  TChain* tree = new TChain("tree"); 
  tree->Add(inputMC1.Data());
  tree->Add(inputMC2.Data());

  mytree->AddFriend(tree);

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  Int_t           Reco_mu_whichGen[maxBranchSize];
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_mu_whichGen;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);


  const int nEP = 29;  // number of event planes in the tree
  Double_t qx[nEP]; 
  Double_t qy[nEP]; 
  TBranch *b_qx;
  TBranch *b_qy;
  tree->SetBranchAddress("qx", qx, &b_qx);
  tree->SetBranchAddress("qy", qy, &b_qy);

  Double_t    	  epang[29];
  TBranch         *b_epang;
  tree->SetBranchAddress("epang", epang, &b_epang);
  
  TString binByBinStr = "";
  if (flattenBinByBin) binByBinStr = "BinByBin";

  TString v2Str;
  if (isAbout(v2,0.5)) v2Str = "0point5";
  else if (isAbout(v2,0.2)) v2Str = "0point2";
  else if (isAbout(v2,0.1)) v2Str = "0point1";
  else if (isAbout(v2,0.05)) v2Str = "0point05";
  else if (isAbout(v2,0.0)) v2Str = "0";

  TFile* newfile;
  newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_MC_flattened%s_order%i_n%i_%i_v2_%s_fixed.root",binByBinStr.Data(),flatOrder,nevt,dateStr,v2Str.Data()),"recreate");

  DiMuon dm;
  TTree* mmtree = new TTree("mmep","dimuonAndEventPlanes");
  mmtree->SetMaxTreeSize(MAXTREESIZE);
  mmtree->Branch("mm",&dm.run,branchString.Data());
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","pt weight", 0, 10000,"");
  RooRealVar* recoQQsign = new RooRealVar("recoQQsign","qq sign",-1,3,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");
  RooRealVar* qxVar   = new RooRealVar("qx","x-comp of 2nd order q vector", -100,100,"");
  RooRealVar* qyVar   = new RooRealVar("qy","y-comp of 2nd order q vector", -100,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar); argSet->add(*ep2Var); argSet->add(*recoQQsign); argSet->add(*dphiEp2Var); argSet->add(*qxVar); argSet->add(*qyVar);

  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  int kTrigSel = 13;
  int DIMUIDPASS = 0;
  int ALLPASS = 0;

  if(nevt == -1) nevt = mytree->GetEntries();

  //averages for resolution correction
  Double_t avgCosAB = 0;
  Double_t avgCosAC = 0;
  Double_t avgCosBC = 0;
  Double_t avgCosABpt[3] = {0};
  Double_t avgCosACpt[3] = {0};
  Double_t avgCosBCpt[3] = {0};
  Double_t avgCosABcent[4] = {0};
  Double_t avgCosACcent[4] = {0};
  Double_t avgCosBCcent[4] = {0};
  //errors on averages
  Double_t sumsqrsCosAB = 0;
  Double_t sumsqrsCosAC = 0;
  Double_t sumsqrsCosBC = 0;
  Double_t sumsqrsCosABpt[3] = {0};
  Double_t sumsqrsCosACpt[3] = {0};
  Double_t sumsqrsCosBCpt[3] = {0};
  Double_t sumsqrsCosABcent[4] = {0};
  Double_t sumsqrsCosACcent[4] = {0};
  Double_t sumsqrsCosBCcent[4] = {0};
  int ptPASS[3] = {0};
  int centPASS[4] = {0};

  //get values for flattening
  TFile* avgFile = TFile::Open(Form("averages/avgEpFile%i.root",nevt),"READ");
  TH1D* havgCosEp = (TH1D*)avgFile->Get("havgCosEp;1");
  TH1D* havgSinEp = (TH1D*)avgFile->Get("havgSinEp;1");
  TH1D* havgCosEpHFm2 = (TH1D*)avgFile->Get("havgCosEpHFm2;1");
  TH1D* havgSinEpHFm2 = (TH1D*)avgFile->Get("havgSinEpHFm2;1");
  TH1D* havgCosEpHFp2 = (TH1D*)avgFile->Get("havgCosEpHFp2;1");
  TH1D* havgSinEpHFp2 = (TH1D*)avgFile->Get("havgSinEpHFp2;1");
  TH1D* havgCosEptrackmid2 = (TH1D*)avgFile->Get("havgCosEptrackmid2;1");
  TH1D* havgSinEptrackmid2 = (TH1D*)avgFile->Get("havgSinEptrackmid2;1");
  Double_t avgCosEp[flatOrder] = {0};
  Double_t avgSinEp[flatOrder] = {0};
  Double_t avgCosEpHFm2[flatOrder] = {0};
  Double_t avgSinEpHFm2[flatOrder] = {0};
  Double_t avgCosEpHFp2[flatOrder] = {0};
  Double_t avgSinEpHFp2[flatOrder] = {0};
  Double_t avgCosEptrackmid2[flatOrder] = {0};
  Double_t avgSinEptrackmid2[flatOrder] = {0};
  for (int n=1; n<=flatOrder; n++) {
    avgCosEp[n-1] = havgCosEp->GetBinContent(n);
    avgSinEp[n-1] = havgSinEp->GetBinContent(n);
    avgCosEpHFm2[n-1] = havgCosEpHFm2->GetBinContent(n);
    avgSinEpHFm2[n-1] = havgSinEpHFm2->GetBinContent(n);
    avgCosEpHFp2[n-1] = havgCosEpHFp2->GetBinContent(n);
    avgSinEpHFp2[n-1] = havgSinEpHFp2->GetBinContent(n);
    avgCosEptrackmid2[n-1] = havgCosEptrackmid2->GetBinContent(n);
    avgSinEptrackmid2[n-1] = havgSinEptrackmid2->GetBinContent(n);
    cout << "avg(Cos(2*" << n << "*Psi)) = " << avgCosEp[n-1] << endl;
    cout << "avg(Sin(2*" << n << "*Psi)) = " << avgSinEp[n-1] << endl;
  }

  TH1D* havgCosEppt[3];
  TH1D* havgSinEppt[3];
  TH1D* havgCosEpHFm2pt[3];
  TH1D* havgSinEpHFm2pt[3];
  TH1D* havgCosEpHFp2pt[3];
  TH1D* havgSinEpHFp2pt[3];
  TH1D* havgCosEptrackmid2pt[3];
  TH1D* havgSinEptrackmid2pt[3];

  TH1D* havgCosEpcent[4];
  TH1D* havgSinEpcent[4];
  TH1D* havgCosEpHFm2cent[4];
  TH1D* havgSinEpHFm2cent[4];
  TH1D* havgCosEpHFp2cent[4];
  TH1D* havgSinEpHFp2cent[4];
  TH1D* havgCosEptrackmid2cent[4];
  TH1D* havgSinEptrackmid2cent[4];

  Double_t avgCosEppt[3][flatOrder] = {0};
  Double_t avgSinEppt[3][flatOrder] = {0};
  Double_t avgCosEpHFm2pt[3][flatOrder] = {0};
  Double_t avgSinEpHFm2pt[3][flatOrder] = {0};
  Double_t avgCosEpHFp2pt[3][flatOrder] = {0};
  Double_t avgSinEpHFp2pt[3][flatOrder] = {0};
  Double_t avgCosEptrackmid2pt[3][flatOrder] = {0};
  Double_t avgSinEptrackmid2pt[3][flatOrder] = {0};

  Double_t avgCosEpcent[4][flatOrder] = {0};
  Double_t avgSinEpcent[4][flatOrder] = {0};
  Double_t avgCosEpHFm2cent[4][flatOrder] = {0};
  Double_t avgSinEpHFm2cent[4][flatOrder] = {0};
  Double_t avgCosEpHFp2cent[4][flatOrder] = {0};
  Double_t avgSinEpHFp2cent[4][flatOrder] = {0};
  Double_t avgCosEptrackmid2cent[4][flatOrder] = {0};
  Double_t avgSinEptrackmid2cent[4][flatOrder] = {0};

  for (int i=0; i<numptbins; i++) {
    havgCosEppt[i] = (TH1D*)avgFile->Get(Form("havgCosEppt[%i];1",i));
    havgSinEppt[i] = (TH1D*)avgFile->Get(Form("havgSinEppt[%i];1",i));
    havgCosEpHFm2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFm2pt[%i];1",i));
    havgSinEpHFm2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFm2pt[%i];1",i));
    havgCosEpHFp2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFp2pt[%i];1",i));
    havgSinEpHFp2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFp2pt[%i];1",i));
    havgCosEptrackmid2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEptrackmid2pt[%i];1",i));
    havgSinEptrackmid2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEptrackmid2pt[%i];1",i));
    for (int n=1; n<=flatOrder; n++) {
      avgCosEppt[i][n-1] = havgCosEppt[i]->GetBinContent(n);
      avgSinEppt[i][n-1] = havgSinEppt[i]->GetBinContent(n);
      avgCosEpHFm2pt[i][n-1] = havgCosEpHFm2pt[i]->GetBinContent(n);
      avgSinEpHFm2pt[i][n-1] = havgSinEpHFm2pt[i]->GetBinContent(n);
      avgCosEpHFp2pt[i][n-1] = havgCosEpHFp2pt[i]->GetBinContent(n);
      avgSinEpHFp2pt[i][n-1] = havgSinEpHFp2pt[i]->GetBinContent(n);
      avgCosEptrackmid2pt[i][n-1] = havgCosEptrackmid2pt[i]->GetBinContent(n);
      avgSinEptrackmid2pt[i][n-1] = havgSinEptrackmid2pt[i]->GetBinContent(n);
    }
  }

  for (int i=0; i<numcbins; i++) {
    havgCosEpcent[i] = (TH1D*)avgFile->Get(Form("havgCosEpcent[%i];1",i));
    havgSinEpcent[i] = (TH1D*)avgFile->Get(Form("havgSinEpcent[%i];1",i));
    havgCosEpHFm2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFm2cent[%i];1",i));
    havgSinEpHFm2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFm2cent[%i];1",i));
    havgCosEpHFp2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFp2cent[%i];1",i));
    havgSinEpHFp2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFp2cent[%i];1",i));
    havgCosEptrackmid2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEptrackmid2cent[%i];1",i));
    havgSinEptrackmid2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEptrackmid2cent[%i];1",i));
    for (int n=1; n<=flatOrder; n++) {
      avgCosEpcent[i][n-1] = havgCosEpcent[i]->GetBinContent(n);
      avgSinEpcent[i][n-1] = havgSinEpcent[i]->GetBinContent(n);
      avgCosEpHFm2cent[i][n-1] = havgCosEpHFm2cent[i]->GetBinContent(n);
      avgSinEpHFm2cent[i][n-1] = havgSinEpHFm2cent[i]->GetBinContent(n);
      avgCosEpHFp2cent[i][n-1] = havgCosEpHFp2cent[i]->GetBinContent(n);
      avgSinEpHFp2cent[i][n-1] = havgSinEpHFp2cent[i]->GetBinContent(n);
      avgCosEptrackmid2cent[i][n-1] = havgCosEptrackmid2cent[i]->GetBinContent(n);
      avgSinEptrackmid2cent[i][n-1] = havgSinEptrackmid2cent[i]->GetBinContent(n);
    }
  }

  avgFile->Close();
  delete avgFile;

  //get values for re-centering
  cout << "getting values for re-centering..." << endl;
  //TFile* avgQFile = TFile::Open("avgQFile_fulldataset_2020_05_20.root","READ");
  TFile* avgQFile = TFile::Open(Form("averages/avgQFile%i.root",nevt),"READ");
  TH1D* havgqx = (TH1D*)avgQFile->Get("hqx;1");
  TH1D* havgqy = (TH1D*)avgQFile->Get("hqy;1");
  TH1D* havgqxHFm2 = (TH1D*)avgQFile->Get("hqxHFm2;1");
  TH1D* havgqyHFm2 = (TH1D*)avgQFile->Get("hqyHFm2;1");
  TH1D* havgqxHFp2 = (TH1D*)avgQFile->Get("hqxHFp2;1");
  TH1D* havgqyHFp2 = (TH1D*)avgQFile->Get("hqyHFp2;1");
  TH1D* havgqxtrackmid2 = (TH1D*)avgQFile->Get("hqxtrackmid2;1");
  TH1D* havgqytrackmid2 = (TH1D*)avgQFile->Get("hqytrackmid2;1");
  Double_t avgqx = havgqx->GetBinContent(1);
  Double_t avgqy = havgqy->GetBinContent(1);
  Double_t avgqxHFm2 = havgqxHFm2->GetBinContent(1);
  Double_t avgqyHFm2 = havgqyHFm2->GetBinContent(1);
  Double_t avgqxHFp2 = havgqxHFp2->GetBinContent(1);
  Double_t avgqyHFp2 = havgqyHFp2->GetBinContent(1);
  Double_t avgqxtrackmid2 = havgqxtrackmid2->GetBinContent(1);
  Double_t avgqytrackmid2 = havgqytrackmid2->GetBinContent(1);
  cout << "avg(qx) = " << avgqx << endl;
  cout << "avg(qy) = " << avgqy << endl;

  TH1D* havgqxpt = (TH1D*)avgQFile->Get("hqxpt;1");
  TH1D* havgqypt = (TH1D*)avgQFile->Get("hqypt;1");
  TH1D* havgqxHFm2pt = (TH1D*)avgQFile->Get("hqxHFm2pt;1");
  TH1D* havgqyHFm2pt = (TH1D*)avgQFile->Get("hqyHFm2pt;1");
  TH1D* havgqxHFp2pt = (TH1D*)avgQFile->Get("hqxHFp2pt;1");
  TH1D* havgqyHFp2pt = (TH1D*)avgQFile->Get("hqyHFp2pt;1");
  TH1D* havgqxtrackmid2pt = (TH1D*)avgQFile->Get("hqxtrackmid2pt;1");
  TH1D* havgqytrackmid2pt = (TH1D*)avgQFile->Get("hqytrackmid2pt;1");
  Double_t avgqxpt[3];
  Double_t avgqypt[3];
  Double_t avgqxHFm2pt[3];
  Double_t avgqyHFm2pt[3];
  Double_t avgqxHFp2pt[3];
  Double_t avgqyHFp2pt[3];
  Double_t avgqxtrackmid2pt[3];
  Double_t avgqytrackmid2pt[3];
  for (int i=0; i<3; i++) {
    avgqxpt[i] = havgqxpt->GetBinContent(i+1);
    avgqypt[i] = havgqypt->GetBinContent(i+1);
    avgqxHFm2pt[i] = havgqxHFm2pt->GetBinContent(i+1);
    avgqyHFm2pt[i] = havgqyHFm2pt->GetBinContent(i+1);
    avgqxHFp2pt[i] = havgqxHFp2pt->GetBinContent(i+1);
    avgqyHFp2pt[i] = havgqyHFp2pt->GetBinContent(i+1);
    avgqxtrackmid2pt[i] = havgqxtrackmid2pt->GetBinContent(i+1);
    avgqytrackmid2pt[i] = havgqytrackmid2pt->GetBinContent(i+1);
  }

  TH1D* havgqxcent = (TH1D*)avgQFile->Get("hqxcent;1");
  TH1D* havgqycent = (TH1D*)avgQFile->Get("hqycent;1");
  TH1D* havgqxHFm2cent = (TH1D*)avgQFile->Get("hqxHFm2cent;1");
  TH1D* havgqyHFm2cent = (TH1D*)avgQFile->Get("hqyHFm2cent;1");
  TH1D* havgqxHFp2cent = (TH1D*)avgQFile->Get("hqxHFp2cent;1");
  TH1D* havgqyHFp2cent = (TH1D*)avgQFile->Get("hqyHFp2cent;1");
  TH1D* havgqxtrackmid2cent = (TH1D*)avgQFile->Get("hqxtrackmid2cent;1");
  TH1D* havgqytrackmid2cent = (TH1D*)avgQFile->Get("hqytrackmid2cent;1");
  Double_t avgqxcent[4];
  Double_t avgqycent[4];
  Double_t avgqxHFm2cent[4];
  Double_t avgqyHFm2cent[4];
  Double_t avgqxHFp2cent[4];
  Double_t avgqyHFp2cent[4];
  Double_t avgqxtrackmid2cent[4];
  Double_t avgqytrackmid2cent[4];
  for (int i=0; i<4; i++) {
    avgqxcent[i] = havgqxcent->GetBinContent(i+1);
    avgqycent[i] = havgqycent->GetBinContent(i+1);
    avgqxHFm2cent[i] = havgqxHFm2cent->GetBinContent(i+1);
    avgqyHFm2cent[i] = havgqyHFm2cent->GetBinContent(i+1);
    avgqxHFp2cent[i] = havgqxHFp2cent->GetBinContent(i+1);
    avgqyHFp2cent[i] = havgqyHFp2cent->GetBinContent(i+1);
    avgqxtrackmid2cent[i] = havgqxtrackmid2cent->GetBinContent(i+1);
    avgqytrackmid2cent[i] = havgqytrackmid2cent->GetBinContent(i+1);
  }
  avgQFile->Close();
  cout << "Done." << endl;

  //check the recentered averages, they should be zero.
  Double_t avgqxrec = 0;
  Double_t avgqyrec = 0;
  Double_t avgqxHFm2rec = 0;
  Double_t avgqyHFm2rec = 0;
  Double_t avgqxHFp2rec = 0;
  Double_t avgqyHFp2rec = 0;
  Double_t avgqxtrackmid2rec = 0;
  Double_t avgqytrackmid2rec = 0;

  //Calculate average q-vector products for the scalar-product method.
  Double_t avgqAqB = 0;//averaged over all events
  Double_t avgqAqC = 0;
  Double_t avgqBqC = 0;
  Double_t avgqqAUps = 0;//averaged over upsilon 1s candidates

  newfile->cd();

  // event loop start

  cout << "Total events = " << nevt << ", : " << mytree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
  
    //Before any cuts, find averages of q-vector products.
    //Nominal: HF2 (8)
    //A: HFm2 (6), B: HFp2 (7), C: trackmid2 (9)
    Double_t qAx = qx[6]; Double_t qAy = qy[6];
    Double_t qBx = qx[7]; Double_t qBy = qy[7];
    Double_t qCx = qx[9]; Double_t qCy = qy[9];
    Double_t qAqB = qAx*qBx-qAy*qBy;
    Double_t qAqC = qAx*qCx-qAy*qCy;
    Double_t qBqC = qBx*qCx-qBy*qCy;
    avgqAqB += qAqB;
    avgqAqC += qAqC;
    avgqBqC += qBqC;

    //Start cutting down to the upsilons.
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      dm.clear();      // clear the output tree: 
      dm.run = runNb;
      dm.lumi = LS ;
      dm.event = eventNb ;
      dm.vz = zVtx;
      dm.cBin = Centrality ;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      
      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;

      bool muplSoft = ( passMuonTypePl && //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) 
          ) ; 

      bool mumiSoft = ( passMuonTypeMi && //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  
          ) ; 

      if ( !(muplSoft && mumiSoft) ) 
        continue;   
      
      DIMUIDPASS++;

      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;
   
      Double_t qxrec = qx[8]-avgqx;
      Double_t qyrec = qy[8]-avgqy;
      Double_t qxHFm2rec = qx[6]-avgqxHFm2;
      Double_t qyHFm2rec = qy[6]-avgqyHFm2;
      Double_t qxHFp2rec = qx[7]-avgqxHFp2;
      Double_t qyHFp2rec = qy[7]-avgqyHFp2;
      Double_t qxtrackmid2rec = qx[9]-avgqxtrackmid2;
      Double_t qytrackmid2rec = qy[9]-avgqytrackmid2;

      Double_t epHFm2rec = atan2(qyHFm2rec,qxHFm2rec)/2;
      Double_t epHFp2rec = atan2(qyHFp2rec,qxHFp2rec)/2;
      Double_t epHF2rec = atan2(qyrec,qxrec)/2;
      Double_t eptrackmid2rec = atan2(qytrackmid2rec,qxtrackmid2rec)/2;

      Double_t epHFm2 = epHFm2rec;
      Double_t epHFp2 = epHFp2rec;
      Double_t epHF2 = epHF2rec;
      Double_t eptrackmid2 = eptrackmid2rec;

      if (abs(epHF2)>5)
        continue;
      if (abs(epHFm2)>5)
        continue;
      if (abs(epHFp2)>5)
        continue;
      if (abs(eptrackmid2)>5)
        continue;

      ALLPASS++;

      //Find averages of q-vector products.
      Double_t q0x = qx[8]; Double_t q0y = qy[8];
      Double_t qqA = q0x*qAx-q0y*qAy;
      avgqqAUps += qqA;

      //Double_t epHFm2 = epang[6];
      //Double_t epHFp2 = epang[7];
      //Double_t eptrackmid2 = epang[9];

      //flatten event plane angles:
      //epang[0] = HFm1
      //epang[1] = HFp1
      //epang[2] = HF1
      //epang[3] = trackm1
      //epang[4] = trackp1
      //epang[5] = Castor1
      //epHFm2 = HFm2*  -5<eta<-3
      //epHFp2 = HFp2*   3<eta<5
      //epang[8] = HF2
      //eptrackmid2 = trackmid2*   -0.8<eta<0.8
      //epang[10] = trackm2
      //epang[11] = trackp2
      //epang[12] = Castor2
      Double_t DeltaEp2 = 0.0;
      Double_t DeltaEpHFm2 = 0.0;
      Double_t DeltaEpHFp2 = 0.0;
      Double_t DeltaEptrackmid2 = 0.0;
      for (int n=1; n<=flatOrder; n++) {
        DeltaEp2 += -cos(2*n*epHF2)*avgSinEp[n-1]/n;
        DeltaEp2 += sin(2*n*epHF2)*avgCosEp[n-1]/n;
        DeltaEpHFm2 += -cos(2*n*epHFm2)*avgSinEpHFm2[n-1]/n;
        DeltaEpHFm2 += sin(2*n*epHFm2)*avgCosEpHFm2[n-1]/n;
        DeltaEpHFp2 += -cos(2*n*epHFp2)*avgSinEpHFp2[n-1]/n;
        DeltaEpHFp2 += sin(2*n*epHFp2)*avgCosEpHFp2[n-1]/n;
        DeltaEptrackmid2 += -cos(2*n*eptrackmid2)*avgSinEptrackmid2[n-1]/n;
        DeltaEptrackmid2 += sin(2*n*eptrackmid2)*avgCosEptrackmid2[n-1]/n;
      }
      if (!flattenBinByBin) {
        epHF2 = epHF2 + DeltaEp2;
        epHFm2 = epHFm2 + DeltaEpHFm2;
        epHFp2 = epHFp2 + DeltaEpHFp2;
        eptrackmid2 = eptrackmid2 + DeltaEptrackmid2;
      }
      if (epHF2<-pi/2) epHF2 = epHF2+pi;
      if (epHF2>pi/2) epHF2 = epHF2-pi;
      if (epHFm2<-pi/2) epHFm2 = epHFm2+pi;
      if (epHFm2>pi/2) epHFm2 = epHFm2-pi;
      if (epHFp2<-pi/2) epHFp2 = epHFp2+pi;
      if (epHFp2>pi/2) epHFp2 = epHFp2-pi;
      if (eptrackmid2<-pi/2) eptrackmid2 = eptrackmid2+pi;
      if (eptrackmid2>pi/2) eptrackmid2 = eptrackmid2-pi;

      //flatten event plane angles bin by bin in centrality:
      int whichptbin = 0;
      int whichcentbin = 0;
      dm.pt     = JP_Reco->Pt();
      if (dm.pt<3) {
        whichptbin = 0;
      }
      else if (dm.pt>=3 && dm.pt<6) {
        whichptbin = 1;
      }
      else if (dm.pt>=6 && dm.pt<30) {
        whichptbin = 2;
      }
      if (dm.cBin<20) {
        whichcentbin = 0;
      }
      else if (dm.cBin>=20 && dm.cBin<60) {
        whichcentbin = 1;
      }
      else if (dm.cBin>=60 && dm.cBin<100) {
        whichcentbin = 2;
      }
      else if (dm.cBin>=100 && dm.cBin<200) {
        whichcentbin = 3;
      }
      Double_t epHF2cent = epHF2;
      Double_t epHFm2cent = epHFm2;
      Double_t epHFp2cent = epHFp2;
      Double_t eptrackmid2cent = eptrackmid2;
      Double_t DeltaEp2cent = 0.0;
      Double_t DeltaEpHFm2cent = 0.0;
      Double_t DeltaEpHFp2cent = 0.0;
      Double_t DeltaEptrackmid2cent = 0.0;
      for (int n=1; n<=flatOrder; n++) {
        DeltaEp2cent += -cos(2*n*epHF2cent)*avgSinEpcent[whichcentbin][n-1]/n;
        DeltaEp2cent += sin(2*n*epHF2cent)*avgCosEpcent[whichcentbin][n-1]/n;
        DeltaEpHFm2cent += -cos(2*n*epHFm2cent)*avgSinEpHFm2cent[whichcentbin][n-1]/n;
        DeltaEpHFm2cent += sin(2*n*epHFm2cent)*avgCosEpHFm2cent[whichcentbin][n-1]/n;
        DeltaEpHFp2cent += -cos(2*n*epHFp2cent)*avgSinEpHFp2cent[whichcentbin][n-1]/n;
        DeltaEpHFp2cent += sin(2*n*epHFp2cent)*avgCosEpHFp2cent[whichcentbin][n-1]/n;
        DeltaEptrackmid2cent += -cos(2*n*eptrackmid2cent)*avgSinEptrackmid2cent[whichcentbin][n-1]/n;
        DeltaEptrackmid2cent += sin(2*n*eptrackmid2cent)*avgCosEptrackmid2cent[whichcentbin][n-1]/n;
      }
      epHF2cent = epHF2cent + DeltaEp2cent;
      epHFm2cent = epHFm2cent + DeltaEpHFm2cent;
      epHFp2cent = epHFp2cent + DeltaEpHFp2cent;
      eptrackmid2cent = eptrackmid2cent + DeltaEptrackmid2cent;
      if (epHF2cent<-pi/2) epHF2cent = epHF2cent+pi;
      if (epHF2cent>pi/2) epHF2cent = epHF2cent-pi;
      if (epHFm2cent<-pi/2) epHFm2cent = epHFm2cent+pi;
      if (epHFm2cent>pi/2) epHFm2cent = epHFm2cent-pi;
      if (epHFp2cent<-pi/2) epHFp2cent = epHFp2cent+pi;
      if (epHFp2cent>pi/2) epHFp2cent = epHFp2cent-pi;
      if (eptrackmid2cent<-pi/2) eptrackmid2cent = eptrackmid2cent+pi;
      if (eptrackmid2cent>pi/2) eptrackmid2cent = eptrackmid2cent-pi;

      //Choose flattening method.
      if (flattenBinByBin) {
        epHF2 = epHF2cent;
        epHFm2 = epHFm2cent;
        epHFp2 = epHFp2cent;
        eptrackmid2 = eptrackmid2cent;
      }

      //Fill event plane histograms.
      hEpHF2->Fill(epHF2);
      hEpHFm2->Fill(epHFm2);
      hEpHFp2->Fill(epHFp2);
      hEptrackmid2->Fill(eptrackmid2);

      //set dimuon attributes
      //dm.phi    = JP_Reco->Phi();
      dm.phi    = restrictToPiOver2(JP_Reco->Phi());
      //dm.ep2 = epang[8];
      dm.ep2 = epHF2;

      //dm.dphiEp2 = getDPHI_Jared( dm.phi, dm.ep2);
      dm.dphiEp2 = restrictToPiOver2(dm.phi-dm.ep2);
      dm.mass   = JP_Reco->M();
      dm.y      = JP_Reco->Rapidity();
      dm.eta    = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();

      //Add weight
      Double_t v2_effective = v2;
      if (v2>0.01) {
        //Double_t v2_present = (v2_present_pt[whichptbin]+v2_present_c[whichcentbin])/2;
        Double_t v2_present = v2_present_pt[whichptbin]+v2_present_c[whichcentbin]-v2_present_average;
        //Double_t v2_present = v2_present_c[whichcentbin];
        //Double_t v2_present = v2_present_pt[whichptbin];
        v2_effective = v2-v2_present;
      }
      if (abs(dm.dphiEp2)<3.1416) dm.weight = 1 + 2*v2_effective*cos(2*dm.dphiEp2);
      else dm.weight = 1.0;

      hphi->Fill(dm.phi);
      hdphi->Fill(dm.dphiEp2);
      hdphiw->Fill(dm.dphiEp2,dm.weight);

      Double_t epHF2old = atan2(qy[8],qx[8])/2;
      Double_t epHFm2old = atan2(qy[6],qx[6])/2;
      Double_t epHFp2old = atan2(qy[7],qx[7])/2;
      Double_t eptrackmid2old = atan2(qy[9],qx[9])/2;
      hEpHF2old->Fill(epHF2old);
      hEpHF2new->Fill(epHF2rec);
      hEpHFm2old->Fill(epHFm2old);
      hEpHFm2new->Fill(epHFm2rec);
      hEpHFp2old->Fill(epHFp2old);
      hEpHFp2new->Fill(epHFp2rec);
      hEptrackmid2old->Fill(eptrackmid2old);
      hEptrackmid2new->Fill(eptrackmid2rec);

      avgqxrec += qxrec;
      avgqyrec += qyrec;
      avgqxHFm2rec += qxHFm2rec;
      avgqyHFm2rec += qyHFm2rec;
      avgqxHFp2rec += qxHFp2rec;
      avgqyHFp2rec += qyHFp2rec;
      avgqxtrackmid2rec += qxtrackmid2rec;
      avgqytrackmid2rec += qytrackmid2rec;

      recoQQsign->setVal((int)Reco_QQ_sign[irqq]);     
      massVar->setVal( (double)dm.mass ) ;
      ptVar->setVal(   (double)dm.pt   ) ;
      yVar->setVal(    (double)dm.y    ) ;
      pt1Var->setVal(  (double)dm.pt1  ) ;
      eta1Var->setVal( (double)dm.eta1 ) ;
      pt2Var->setVal(  (double)dm.pt2  ) ;
      eta2Var->setVal( (double)dm.eta2 ) ;
      ep2Var->setVal( (double)dm.ep2 ) ;
      dphiEp2Var->setVal(   (double)dm.dphiEp2  ) ;
      cBinVar->setVal( (double)dm.cBin ) ;
      evtWeight->setVal( (double)dm.weight ) ;
      dataSet->add( *argSet);

      //get values to use for event plane resolution.
      Double_t psiA = 0;
      Double_t psiB = 0;
      Double_t psiC = 0;
      if (dm.eta<0.0) {
        psiA = epHFp2;
        psiB = epHFm2;
        psiC = eptrackmid2;
      }
      else {
        psiA = epHFm2;
        psiB = epHFp2;
        psiC = eptrackmid2;
      }

      Double_t deltaAB = psiA-psiB;
      Double_t deltaAC = psiA-psiC;
      Double_t deltaBC = psiB-psiC;

      avgCosAB += cos(2*(deltaAB));
      avgCosAC += cos(2*(deltaAC));
      avgCosBC += cos(2*(deltaBC));

      //cout << "CosAB = " << cos(2*(deltaAB)) << endl;
      //hCosAC->Fill(cos(2*(deltaAC)));

      sumsqrsCosAB += pow(cos(2*(deltaAB)),2);
      sumsqrsCosAC += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBC += pow(cos(2*(deltaBC)),2);

      //binned averages:
      avgCosABpt[whichptbin] += cos(2*(deltaAB));
      avgCosACpt[whichptbin] += cos(2*(deltaAC));
      avgCosBCpt[whichptbin] += cos(2*(deltaBC));
      sumsqrsCosABpt[whichptbin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACpt[whichptbin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCpt[whichptbin] += pow(cos(2*(deltaBC)),2);
      ptPASS[whichptbin] += 1;

      avgCosABcent[whichcentbin] += cos(2*(deltaAB));
      avgCosACcent[whichcentbin] += cos(2*(deltaAC));
      avgCosBCcent[whichcentbin] += cos(2*(deltaBC));
      sumsqrsCosABcent[whichcentbin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACcent[whichcentbin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCcent[whichcentbin] += pow(cos(2*(deltaBC)),2);
      centPASS[whichcentbin] += 1;

      mmtree->Fill();

    } // end of dimuon loop
  } //end of event loop

  mmtree->Write();  // Don't need to call Write() for trees
  dataSet->Write(); 
  newfile->Close();

  //Averaged for scalar product method:
  TFile* avgSPFile = new TFile(Form("averages/avgSPFile_n%i_%i.root",nevt,dateStr),"RECREATE");
  TH1D* havgSP = new TH1D("havgSP","avgSP",4,0,4);
  cout << "sumqqAUps = " << avgqqAUps << endl;
  cout << "sumqAqB = " << avgqAqB << endl;
  cout << "sumqAqC = " << avgqAqC << endl;
  cout << "sumqBqC = " << avgqBqC << endl;
  avgqAqB = avgqAqB/nevt;
  avgqAqC = avgqAqC/nevt;
  avgqBqC = avgqBqC/nevt;
  avgqqAUps = avgqqAUps/ALLPASS;
  cout << "avgqqAUps = " << avgqqAUps << endl;
  cout << "avgqAqB = " << avgqAqB << endl;
  cout << "avgqAqC = " << avgqAqC << endl;
  cout << "avgqBqC = " << avgqBqC << endl;
  havgSP->SetBinContent(1,avgqqAUps);
  havgSP->SetBinContent(2,avgqAqB);
  havgSP->SetBinContent(3,avgqAqC);
  havgSP->SetBinContent(4,avgqBqC);
  havgSP->Write();
  TH1D* hv2SP = new TH1D("hv2SP","v2SP",1,0,1);
  Double_t v2SP = avgqqAUps/sqrt(avgqAqB*avgqAqC/avgqBqC);
  cout << "v2SP = " << v2SP << endl;
  hv2SP->SetBinContent(1,v2SP);
  hv2SP->Write();
  avgSPFile->Close();

  //include the averages used to calculate the event plane resolution correction.
  TFile* resCorFile = new TFile(Form("averages/resCorFile_n%i_%s.root",nevt,binByBinStr.Data()),"RECREATE");
  TH1D* hRint = new TH1D("hRint","EP Resolution factor",1,0,1);
  avgCosAB = avgCosAB/ALLPASS;
  avgCosAC = avgCosAC/ALLPASS;
  avgCosBC = avgCosBC/ALLPASS;
  Double_t rmsCosAB = sqrt(sumsqrsCosAB/ALLPASS - pow(avgCosAB,2));
  Double_t rmsCosAC = sqrt(sumsqrsCosAC/ALLPASS - pow(avgCosAC,2));
  Double_t rmsCosBC = sqrt(sumsqrsCosBC/ALLPASS - pow(avgCosBC,2));
  cout << "avg(Cos(2*(deltaAB))) = " << avgCosAB << " +/- " << rmsCosAB << endl;
  cout << "avg(Cos(2*(deltaAC))) = " << avgCosAC << " +/- " << rmsCosAC << endl;
  cout << "avg(Cos(2*(deltaBC))) = " << avgCosBC << " +/- " << rmsCosBC << endl;
  rmsCosAB = rmsCosAB/sqrt(ALLPASS);
  rmsCosAC = rmsCosAC/sqrt(ALLPASS);
  rmsCosBC = rmsCosBC/sqrt(ALLPASS);
  Double_t RA = sqrt(avgCosAB*avgCosAC/avgCosBC);
  Double_t RAerr = 0.5*RA*sqrt(pow(rmsCosAB/avgCosAB,2) + pow(rmsCosAC/avgCosAC,2) + pow(rmsCosBC/avgCosBC,2));
  cout << "Event plane resolution factor = " << RA << " +/- " << RAerr << endl;
  hRint->SetBinContent(1,RA);
  hRint->SetBinError(1,RAerr);
  hRint->Write();

  cout << "filling pt histogram" << endl;
  for (int i=0; i<numptbins; i++) {
    avgCosABpt[i] = avgCosABpt[i]/ptPASS[i];
    avgCosACpt[i] = avgCosACpt[i]/ptPASS[i];
    avgCosBCpt[i] = avgCosBCpt[i]/ptPASS[i];
    Double_t rmsCosABpt = sqrt(sumsqrsCosABpt[i]/ptPASS[i] - pow(avgCosABpt[i],2));
    Double_t rmsCosACpt = sqrt(sumsqrsCosACpt[i]/ptPASS[i] - pow(avgCosACpt[i],2));
    Double_t rmsCosBCpt = sqrt(sumsqrsCosBCpt[i]/ptPASS[i] - pow(avgCosBCpt[i],2));
    rmsCosABpt = rmsCosABpt/sqrt(ptPASS[i]);
    rmsCosACpt = rmsCosACpt/sqrt(ptPASS[i]);
    rmsCosBCpt = rmsCosBCpt/sqrt(ptPASS[i]);
    Double_t RApt = sqrt(avgCosABpt[i]*avgCosACpt[i]/avgCosBCpt[i]);
    Double_t RApterr = 0.5*RApt*sqrt(pow(rmsCosABpt/avgCosABpt[i],2) + pow(rmsCosACpt/avgCosACpt[i],2) + pow(rmsCosBCpt/avgCosBCpt[i],2));
    hRpt->SetBinContent(i+1,RApt);
    hRpt->SetBinError(i+1,RApterr);
  }
  hRpt->Write();

  cout << "filling cent histogram" << endl;
  for (int i=0; i<numcbins; i++) {
    avgCosABcent[i] = avgCosABcent[i]/centPASS[i];
    avgCosACcent[i] = avgCosACcent[i]/centPASS[i];
    avgCosBCcent[i] = avgCosBCcent[i]/centPASS[i];
    Double_t rmsCosABcent = sqrt(sumsqrsCosABcent[i]/centPASS[i] - pow(avgCosABcent[i],2));
    Double_t rmsCosACcent = sqrt(sumsqrsCosACcent[i]/centPASS[i] - pow(avgCosACcent[i],2));
    Double_t rmsCosBCcent = sqrt(sumsqrsCosBCcent[i]/centPASS[i] - pow(avgCosBCcent[i],2));
    rmsCosABcent = rmsCosABcent/sqrt(centPASS[i]);
    rmsCosACcent = rmsCosACcent/sqrt(centPASS[i]);
    rmsCosBCcent = rmsCosBCcent/sqrt(centPASS[i]);
    Double_t RAcent = sqrt(avgCosABcent[i]*avgCosACcent[i]/avgCosBCcent[i]);
    Double_t RAcenterr = 0.5*RAcent*sqrt(pow(rmsCosABcent/avgCosABcent[i],2) + pow(rmsCosACcent/avgCosACcent[i],2) + pow(rmsCosBCcent/avgCosBCcent[i],2));
    hRcent->SetBinContent(i+1,RAcent);
    hRcent->SetBinError(i+1,RAcenterr);
  }
  hRcent->Write();

  hRpt->Sumw2();
  hRcent->Sumw2();

  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,300);
  c1->Divide(2);
  c1->cd(1);
  hRpt->Draw();
  c1->cd(2);
  hRcent->Draw();

  hEpHF2->Write();
  hEpHFm2->Write();
  hEpHFp2->Write();
  hEptrackmid2->Write();

  TCanvas* c2 = new TCanvas("c2","c2",0,0,800,800);
  c2->Divide(2,2);
  c2->cd(1);
  hEpHF2->Draw();
  c2->cd(2);
  hEpHFm2->Draw();
  c2->cd(3);
  hEpHFp2->Draw();
  c2->cd(4);
  hEptrackmid2->Draw();

  c2->SaveAs(Form("plots/epHistos%s_n%i_%i.pdf",binByBinStr.Data(),nevt,dateStr));
  c2->SaveAs(Form("plots/epHistos%s_n%i_%i.png",binByBinStr.Data(),nevt,dateStr));

  TCanvas* c3 = new TCanvas("c3","c3",0,0,400,400);
  c3->cd();
  hEpHF2old->SetTitle("Event plane HF2");
  hEpHF2old->Draw();
  hEpHF2new->SetLineColor(2);
  hEpHF2new->Draw("same");
  hEpHF2->SetLineColor(3);
  hEpHF2->Draw("same");
  TLegend* legc3 = new TLegend(0.4,0.15,0.6,0.3); legc3->SetTextSize(12);
  legc3->SetTextFont(43);
  legc3->SetBorderSize(0);
  legc3->AddEntry(hEpHF2old,"ep raw","l");
  legc3->AddEntry(hEpHF2new,"ep recentered","l");
  legc3->AddEntry(hEpHF2,"ep rec+flattened","l");
  legc3->Draw("same");
  c3->SaveAs(Form("plots/EventPlaneHF2Change_n%i.pdf",nevt));
  c3->SaveAs(Form("plots/EventPlaneHF2Change_n%i.png",nevt));

  TCanvas* c4 = new TCanvas("c4","c4",0,0,400,400);
  c4->cd();
  hEpHFm2old->SetTitle("Event plane HFm2");
  hEpHFm2old->Draw();
  hEpHFm2new->SetLineColor(2);
  hEpHFm2new->Draw("same");
  hEpHFm2->SetLineColor(3);
  hEpHFm2->Draw("same");
  TLegend* legc4 = new TLegend(0.4,0.15,0.6,0.3); legc4->SetTextSize(12);
  legc4->SetTextFont(43);
  legc4->SetBorderSize(0);
  legc4->AddEntry(hEpHFm2old,"ep raw","l");
  legc4->AddEntry(hEpHFm2new,"ep recentered","l");
  legc4->AddEntry(hEpHFm2,"ep rec+flattened","l");
  legc4->Draw("same");
  c4->SaveAs(Form("plots/EventPlaneHFm2Change_n%i.pdf",nevt));
  c4->SaveAs(Form("plots/EventPlaneHFm2Change_n%i.png",nevt));

  TCanvas* c5 = new TCanvas("c5","c5",0,0,400,400);
  c5->cd();
  hEpHFp2old->SetTitle("Event plane HFp2");
  hEpHFp2old->Draw();
  hEpHFp2new->SetLineColor(2);
  hEpHFp2new->Draw("same");
  hEpHFp2->SetLineColor(3);
  hEpHFp2->Draw("same");
  TLegend* legc5 = new TLegend(0.4,0.15,0.6,0.3); legc5->SetTextSize(12);
  legc5->SetTextFont(43);
  legc5->SetBorderSize(0);
  legc5->AddEntry(hEpHFp2old,"ep raw","l");
  legc5->AddEntry(hEpHFp2new,"ep recentered","l");
  legc5->AddEntry(hEpHFp2,"ep rec+flattened","l");
  legc5->Draw("same");
  c5->SaveAs(Form("plots/EventPlaneHFp2Change_n%i.pdf",nevt));
  c5->SaveAs(Form("plots/EventPlaneHFp2Change_n%i.png",nevt));

  TCanvas* c6 = new TCanvas("c6","c6",0,0,400,400);
  c6->cd();
  hEptrackmid2old->SetTitle("Event plane trackmid2");
  hEptrackmid2old->Draw();
  hEptrackmid2new->SetLineColor(2);
  hEptrackmid2new->Draw("same");
  hEptrackmid2->SetLineColor(3);
  hEptrackmid2->Draw("same");
  TLegend* legc6 = new TLegend(0.4,0.15,0.6,0.3); legc6->SetTextSize(12);
  legc6->SetTextFont(43);
  legc6->SetBorderSize(0);
  legc6->AddEntry(hEptrackmid2old,"ep raw","l");
  legc6->AddEntry(hEptrackmid2new,"ep recentered","l");
  legc6->AddEntry(hEptrackmid2,"ep rec+flattened","l");
  legc6->Draw("same");
  c6->SaveAs(Form("plots/EventPlanetrackmid2Change_n%i.pdf",nevt));
  c6->SaveAs(Form("plots/EventPlanetrackmid2Change_n%i.png",nevt));

  resCorFile->Close();

  TCanvas* c7 = new TCanvas("c7","c7",0,0,400,400);
  c7->cd();
  hdphiw->Draw();
  hdphiw->SetLineColor(2);
  hdphi->Draw("same");

  TCanvas* c8 = new TCanvas("c8","c8",400,0,400,400);
  c8->cd();
  hphi->Draw();

  cout << "ptPASS = {" << ptPASS[0] << "," << ptPASS[1] << "," << ptPASS[2] << "}" << endl;
  cout << "centPASS = {" << centPASS[0] << "," << centPASS[1] << "," << centPASS[2] << "," << centPASS[3] << "}" << endl;

  cout << endl;
  cout << "avgqx = " << avgqx << endl;
  cout << "avgqy = " << avgqy << endl;
  cout << "avgqxHFm2 = " << avgqxHFm2 << endl;
  cout << "avgqyHFm2 = " << avgqyHFm2 << endl;
  cout << "avgqxHFp2 = " << avgqxHFp2 << endl;
  cout << "avgqyHFp2 = " << avgqyHFp2 << endl;
  cout << "avgqxtrackmid2 = " << avgqxtrackmid2 << endl;
  cout << "avgqytrackmid2 = " << avgqytrackmid2 << endl;

  cout << endl << "All these averages should be zero after recentering:" << endl;
  avgqxrec = avgqxrec/ALLPASS;
  avgqyrec = avgqyrec/ALLPASS;
  cout << "avgqx-recentered = " << avgqxrec << endl;
  cout << "avgqy-recentered = " << avgqyrec << endl;
  avgqxHFm2rec = avgqxHFm2rec/ALLPASS;
  avgqyHFm2rec = avgqyHFm2rec/ALLPASS;
  cout << "avgqxHFm2-recentered = " << avgqxHFm2rec << endl;
  cout << "avgqyHFm2-recentered = " << avgqyHFm2rec << endl;
  avgqxHFp2rec = avgqxHFp2rec/ALLPASS;
  avgqyHFp2rec = avgqyHFp2rec/ALLPASS;
  cout << "avgqxHFp2-recentered = " << avgqxHFp2rec << endl;
  cout << "avgqyHFp2-recentered = " << avgqyHFp2rec << endl;
  avgqxtrackmid2rec = avgqxtrackmid2rec/ALLPASS;
  avgqytrackmid2rec = avgqytrackmid2rec/ALLPASS;
  cout << "avgqxtrackmid2-recentered = " << avgqxtrackmid2rec << endl;
  cout << "avgqytrackmid2-recentered = " << avgqytrackmid2rec << endl;
} 

