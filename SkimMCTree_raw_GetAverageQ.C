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


static const long MAXTREESIZE = 1000000000000;

void SkimMCTree_raw_GetAverageQ(int nevt=-1, int dateStr=20210107) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  Double_t ptbins[4] = {0,3,6,30};
  const int numptbins = sizeof(ptbins)/sizeof(Double_t)-1;

  Double_t centbins[5] = {0,10,30,50,100};
  const int numcentbins = sizeof(centbins)/sizeof(Double_t)-1;

  //TH1D* hPseudoQx = new TH1D("hPseudoQx","cos(2*psi)",50,-1.2,1.2);
  //TH1D* hPseudoQy = new TH1D("hPseudoQy","sin(2*psi)",50,-1.2,1.2);

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
  
  TFile* newfile;
  newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_MC_RAW_%i.root",dateStr),"recreate");

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

  //averages for re-centering
  Double_t avgqx = 0;
  Double_t avgqy = 0;
  Double_t avgqxHFm2 = 0;
  Double_t avgqyHFm2 = 0;
  Double_t avgqxHFp2 = 0;
  Double_t avgqyHFp2 = 0;
  Double_t avgqxtrackmid2 = 0;
  Double_t avgqytrackmid2 = 0;

  Double_t avgqxpt[3] = {0};
  Double_t avgqypt[3] = {0};
  Double_t avgqxHFm2pt[3] = {0};
  Double_t avgqyHFm2pt[3] = {0};
  Double_t avgqxHFp2pt[3] = {0};
  Double_t avgqyHFp2pt[3] = {0};
  Double_t avgqxtrackmid2pt[3] = {0};
  Double_t avgqytrackmid2pt[3] = {0};

  Double_t avgqxcent[4] = {0};
  Double_t avgqycent[4] = {0};
  Double_t avgqxHFm2cent[4] = {0};
  Double_t avgqyHFm2cent[4] = {0};
  Double_t avgqxHFp2cent[4] = {0};
  Double_t avgqyHFp2cent[4] = {0};
  Double_t avgqxtrackmid2cent[4] = {0};
  Double_t avgqytrackmid2cent[4] = {0};


  //errors on averages
  Double_t sumsqrsqx = 0;
  Double_t sumsqrsqy = 0;

  int ptPASS[3] = {0};
  int centPASS[4] = {0};

  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
  
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    if( Centrality<20 || Centrality>180 ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {

      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;

      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      
      bool passMuonTypePl = true;
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

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

      dm.clear();      // clear the output tree: 
      dm.run = runNb;
      dm.lumi = LS ;
      dm.event = eventNb ;
      dm.vz = zVtx;
      dm.cBin = Centrality ;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      Double_t qxHF2 = qx[8];
      Double_t qyHF2 = qy[8];
      Double_t qxHFm2 = qx[6];
      Double_t qyHFm2 = qy[6];
      Double_t qxHFp2 = qx[7];
      Double_t qyHFp2 = qy[7];
      Double_t qxtrackmid2 = qx[9];
      Double_t qytrackmid2 = qy[9];

      Double_t epHFm2 = atan2(qyHFm2,qxHFm2)/2;
      Double_t epHFp2 = atan2(qyHFp2,qxHFp2)/2;
      Double_t epHF2 = atan2(qyHF2,qxHF2)/2;
      Double_t eptrackmid2 = atan2(qytrackmid2,qxtrackmid2)/2;

      if (abs(epHF2)>5)
        continue;
      if (abs(epHFm2)>5)
        continue;
      if (abs(epHFp2)>5)
        continue;
      if (abs(eptrackmid2)>5)
        continue;

      ALLPASS++;

      dm.phi    = JP_Reco->Phi();
      //dm.ep2 = epang[8];
      dm.ep2 = epHF2;
      dm.dphiEp2 = getDPHI_Jared( dm.phi, dm.ep2);

      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();

      dm.y      = JP_Reco->Rapidity();
      dm.eta      = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();
      dm.weight = 1.;

      //cout << "qx[8] = " << qx[8] << endl;
      //a little test:
      //cout << "epang[8] = " << epang[8] << endl;
      //cout << "atan(qy/qx)/2 = " << atan2(qy[8],qx[8])/2 << endl;

      recoQQsign->setVal((int)Reco_QQ_sign[irqq]);     
      massVar->setVal( (Double_t)dm.mass ) ;
      ptVar->setVal(   (Double_t)dm.pt   ) ;
      yVar->setVal(    (Double_t)dm.y    ) ;
      pt1Var->setVal(  (Double_t)dm.pt1  ) ;
      eta1Var->setVal( (Double_t)dm.eta1 ) ;
      pt2Var->setVal(  (Double_t)dm.pt2  ) ;
      eta2Var->setVal( (Double_t)dm.eta2 ) ;
      ep2Var->setVal( (Double_t)dm.ep2 ) ;
      dphiEp2Var->setVal(   (Double_t)dm.dphiEp2  ) ;
      qxVar->setVal(   (Double_t)qx[8]  ) ;
      qyVar->setVal(   (Double_t)qy[8]  ) ;
      cBinVar->setVal( (Double_t)dm.cBin ) ;
      evtWeight->setVal( (Double_t)dm.weight ) ;
      dataSet->add( *argSet);

      //get average values for recentering later on.
      //epang[0] = HFm1
      //epang[1] = HFp1
      //epang[2] = HF1
      //epang[3] = trackm1
      //epang[4] = trackp1
      //epang[5] = Castor1
      //epang[6] = HFm2*  -5<eta<-3
      //epang[7] = HFp2*   3<eta<5
      //epang[8] = HF2
      //epang[9] = trackmid2*   -0.8<eta<0.8
      //epang[10] = trackm2
      //epang[11] = trackp2
      //epang[12] = Castor2
      avgqx += qx[8];
      avgqy += qy[8];
      avgqxHFm2 += qx[6];
      avgqyHFm2 += qy[6];
      avgqxHFp2 += qx[7];
      avgqyHFp2 += qy[7];
      avgqxtrackmid2 += qx[9];
      avgqytrackmid2 += qy[9];

      sumsqrsqx += pow(qx[8],2);
      sumsqrsqy += pow(qy[8],2);

      //binned averages:
      int whichptbin = 0;
      int whichcentbin = 0;
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

      avgqxpt[whichptbin] += qx[8];
      avgqypt[whichptbin] += qy[8];
      avgqxHFm2pt[whichcentbin] += qx[6];
      avgqyHFm2pt[whichcentbin] += qy[6];
      avgqxHFp2pt[whichcentbin] += qx[7];
      avgqyHFp2pt[whichcentbin] += qy[7];
      avgqxtrackmid2pt[whichcentbin] += qx[9];
      avgqytrackmid2pt[whichcentbin] += qy[9];

      avgqxcent[whichcentbin] += qx[8];
      avgqycent[whichcentbin] += qy[8];
      avgqxHFm2cent[whichcentbin] += qx[6];
      avgqyHFm2cent[whichcentbin] += qy[6];
      avgqxHFp2cent[whichcentbin] += qx[7];
      avgqyHFp2cent[whichcentbin] += qy[7];
      avgqxtrackmid2cent[whichcentbin] += qx[9];
      avgqytrackmid2cent[whichcentbin] += qy[9];

      ptPASS[whichptbin] += 1;

      centPASS[whichcentbin] += 1;

      mmtree->Fill();

    } // end of dimuon loop
  } //end of event loop

  mmtree->Write();  // Don't need to call Write() for trees
  dataSet->Write(); 
  newfile->Close();

  //TCanvas* c2 = new TCanvas("c2","c2",0,0,400,400);
  //hCosAC->Draw();
  //hCosAC->SaveAs("hCosAC.pdf");

  //Save file with averages for flattening.
  TFile* avgFile = new TFile(Form("averages/avgQFile%i.root",nevt),"RECREATE");

  TH1D* hqx = new TH1D("hqx","Average of qx",1,0,30);
  TH1D* hqy = new TH1D("hqy","Average of qy",1,0,30);
  TH1D* hqxHFm2 = new TH1D("hqxHFm2","Average of qx",1,0,30);
  TH1D* hqyHFm2 = new TH1D("hqyHFm2","Average of qy",1,0,30);
  TH1D* hqxHFp2 = new TH1D("hqxHFp2","Average of qx",1,0,30);
  TH1D* hqyHFp2 = new TH1D("hqyHFp2","Average of qy",1,0,30);
  TH1D* hqxtrackmid2 = new TH1D("hqxtrackmid2","Average of qx",1,0,30);
  TH1D* hqytrackmid2 = new TH1D("hqytrackmid2","Average of qy",1,0,30);

  TH1D* hqxpt = new TH1D("hqxpt","Average of qx",numptbins,ptbins);
  TH1D* hqypt = new TH1D("hqypt","Average of qy",numptbins,ptbins);
  TH1D* hqxHFm2pt = new TH1D("hqxHFm2pt","Average of qx",numptbins,ptbins);
  TH1D* hqyHFm2pt = new TH1D("hqyHFm2pt","Average of qy",numptbins,ptbins);
  TH1D* hqxHFp2pt = new TH1D("hqxHFp2pt","Average of qx",numptbins,ptbins);
  TH1D* hqyHFp2pt = new TH1D("hqyHFp2pt","Average of qy",numptbins,ptbins);
  TH1D* hqxtrackmid2pt = new TH1D("hqxtrackmid2pt","Average of qx",numptbins,ptbins);
  TH1D* hqytrackmid2pt = new TH1D("hqytrackmid2pt","Average of qy",numptbins,ptbins);

  TH1D* hqxcent = new TH1D("hqxcent","Average of qx",numcentbins,centbins);
  TH1D* hqycent = new TH1D("hqycent","Average of qy",numcentbins,centbins);
  TH1D* hqxHFm2cent = new TH1D("hqxHFm2cent","Average of qx",numcentbins,centbins);
  TH1D* hqyHFm2cent = new TH1D("hqyHFm2cent","Average of qy",numcentbins,centbins);
  TH1D* hqxHFp2cent = new TH1D("hqxHFp2cent","Average of qx",numcentbins,centbins);
  TH1D* hqyHFp2cent = new TH1D("hqyHFp2cent","Average of qy",numcentbins,centbins);
  TH1D* hqxtrackmid2cent = new TH1D("hqxtrackmid2cent","Average of qx",numcentbins,centbins);
  TH1D* hqytrackmid2cent = new TH1D("hqytrackmid2cent","Average of qy",numcentbins,centbins);

    cout << "sum(qx) = " << avgqx << endl;
    cout << "sum(qx) = " << avgqy << endl;
    avgqx = avgqx/ALLPASS;
    avgqy = avgqy/ALLPASS;
    avgqxHFm2 = avgqxHFm2/ALLPASS;
    avgqyHFm2 = avgqyHFm2/ALLPASS;
    avgqxHFp2 = avgqxHFp2/ALLPASS;
    avgqyHFp2 = avgqyHFp2/ALLPASS;
    avgqxtrackmid2 = avgqxtrackmid2/ALLPASS;
    avgqytrackmid2 = avgqytrackmid2/ALLPASS;
    cout << "avg(qx) = " << avgqx << endl;
    cout << "avg(qy) = " << avgqy << endl;
    hqx->SetBinContent(1,avgqx);
    hqy->SetBinContent(1,avgqy);
    hqxHFm2->SetBinContent(1,avgqxHFm2);
    hqyHFm2->SetBinContent(1,avgqyHFm2);
    hqxHFp2->SetBinContent(1,avgqxHFp2);
    hqyHFp2->SetBinContent(1,avgqyHFp2);
    hqxtrackmid2->SetBinContent(1,avgqxtrackmid2);
    hqytrackmid2->SetBinContent(1,avgqytrackmid2);
  
  hqx->Write();
  hqy->Write();
  hqxHFm2->Write();
  hqyHFm2->Write();
  hqxHFp2->Write();
  hqyHFp2->Write();
  hqxtrackmid2->Write();
  hqytrackmid2->Write();

  for (int i=0; i<numptbins; i++) {
      avgqxpt[i] = avgqxpt[i]/ptPASS[i];
      avgqypt[i] = avgqypt[i]/ptPASS[i];
      avgqxHFm2pt[i] = avgqxHFm2pt[i]/ptPASS[i];
      avgqyHFm2pt[i] = avgqyHFm2pt[i]/ptPASS[i];
      avgqxHFp2pt[i] = avgqxHFp2pt[i]/ptPASS[i];
      avgqyHFp2pt[i] = avgqyHFp2pt[i]/ptPASS[i];
      avgqxtrackmid2pt[i] = avgqxtrackmid2pt[i]/ptPASS[i];
      avgqytrackmid2pt[i] = avgqytrackmid2pt[i]/ptPASS[i];
      hqxpt->SetBinContent(i+1,avgqxpt[i]);
      hqypt->SetBinContent(i+1,avgqypt[i]);
      hqxHFm2pt->SetBinContent(i+1,avgqxHFm2pt[i]);
      hqyHFm2pt->SetBinContent(i+1,avgqyHFm2pt[i]);
      hqxHFp2pt->SetBinContent(i+1,avgqxHFp2pt[i]);
      hqyHFp2pt->SetBinContent(i+1,avgqyHFp2pt[i]);
      hqxtrackmid2pt->SetBinContent(i+1,avgqxtrackmid2pt[i]);
      hqytrackmid2pt->SetBinContent(i+1,avgqytrackmid2pt[i]);
  }
  hqxpt->Write();
  hqypt->Write();
  hqxHFm2pt->Write();
  hqyHFm2pt->Write();
  hqxHFp2pt->Write();
  hqyHFp2pt->Write();
  hqxtrackmid2pt->Write();
  hqytrackmid2pt->Write();

  for (int i=0; i<numcentbins; i++) {
      avgqxcent[i] = avgqxcent[i]/centPASS[i];
      avgqycent[i] = avgqycent[i]/centPASS[i];
      avgqxHFm2cent[i] = avgqxHFm2cent[i]/centPASS[i];
      avgqyHFm2cent[i] = avgqyHFm2cent[i]/centPASS[i];
      avgqxHFp2cent[i] = avgqxHFp2cent[i]/centPASS[i];
      avgqyHFp2cent[i] = avgqyHFp2cent[i]/centPASS[i];
      avgqxtrackmid2cent[i] = avgqxtrackmid2cent[i]/centPASS[i];
      avgqytrackmid2cent[i] = avgqytrackmid2cent[i]/centPASS[i];
      hqxcent->SetBinContent(i+1,avgqxcent[i]);
      hqycent->SetBinContent(i+1,avgqycent[i]);
      hqxHFm2cent->SetBinContent(i+1,avgqxHFm2cent[i]);
      hqyHFm2cent->SetBinContent(i+1,avgqyHFm2cent[i]);
      hqxHFp2cent->SetBinContent(i+1,avgqxHFp2cent[i]);
      hqyHFp2cent->SetBinContent(i+1,avgqyHFp2cent[i]);
      hqxtrackmid2cent->SetBinContent(i+1,avgqxtrackmid2cent[i]);
      hqytrackmid2cent->SetBinContent(i+1,avgqytrackmid2cent[i]);
  }
  hqxcent->Write();
  hqycent->Write();
  hqxHFm2cent->Write();
  hqyHFm2cent->Write();
  hqxHFp2cent->Write();
  hqyHFp2cent->Write();
  hqxtrackmid2cent->Write();
  hqytrackmid2cent->Write();

  avgFile->Close();

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

}

