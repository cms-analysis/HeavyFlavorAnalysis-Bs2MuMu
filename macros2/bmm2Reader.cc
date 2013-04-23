#include "bmm2Reader.hh"
#include <cmath>
#include <string>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "MuScleFitCorrector.hh"

#include "candAna.hh"
#include "candAnaMuMu.hh"
#include "candAnaBu2JpsiK.hh"
#include "candAnaBs2JpsiPhi.hh"
#include "candAnaDstar.hh"
#include "candAnaHh.hh"
#include "candAnaBd2DstarPi.hh"

using namespace std;

// ----------------------------------------------------------------------
bmm2Reader::bmm2Reader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> bmm2Reader: constructor..." << endl;
  fVerbose = 0; 
  msc = 0; 
}


// ----------------------------------------------------------------------
bmm2Reader::~bmm2Reader() {
  cout << "==> bmm2Reader: destructor..." << endl;
}


// ----------------------------------------------------------------------
void bmm2Reader::startAnalysis() {
  cout << "==> bmm2Reader: fVerbose = " << fVerbose << endl;
  fpJSON = new JSON(JSONFILE.c_str(), (fVerbose!=0?1:0)); 
}


// ----------------------------------------------------------------------
void bmm2Reader::eventProcessing() {
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(0); 

  bool json = false;  
   
  if (fIsMC) {
    json = 1; 
    processType(); 
  } else {
    json = fpJSON->good(fRun, fLS); 
    if (fVerbose > 1 && !json) {
      cout << "JSON = 0 for run = " << fRun << " and LS = " << fLS << endl;
    }
    fProcessType = -98; 
  }
  
  // -- fill a few basic histograms
  TSimpleTrack *pT(0); 
  double x(0.);
  ((TH1D*)fpHistFile->Get("ntracks"))->Fill(fpEvt->nSimpleTracks()); 
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    x = pT->getP().Perp();
    ((TH1D*)fpHistFile->Get("pt0"))->Fill(x); 
    ((TH1D*)fpHistFile->Get("pt1"))->Fill(x); 
    x = pT->getP().Eta();
    ((TH1D*)fpHistFile->Get("eta"))->Fill(x); 
    x = pT->getP().Phi();
    ((TH1D*)fpHistFile->Get("phi"))->Fill(x); 
  }
  
  // -- call candidate analyses
  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->fIsMC        = fIsMC;
    lCandAnalysis[i]->fJSON        = json; 
    lCandAnalysis[i]->fRun         = fRun; 
    lCandAnalysis[i]->fEvt         = fEvt; 
    lCandAnalysis[i]->fLS          = fLS; 
    lCandAnalysis[i]->fEvent       = fEvent; 
    lCandAnalysis[i]->fProcessType = fProcessType; 
    lCandAnalysis[i]->fCandTau     = -1.; 
    lCandAnalysis[i]->fGenLifeTime = -1.; 

    //cout<<" call evtanalysis "<<i<<endl;
    lCandAnalysis[i]->evtAnalysis(fpEvt);
  }

}


// ----------------------------------------------------------------------
void bmm2Reader::bookHist() {
  fpHistFile->cd();
  TH1D *h(0);
  h = new TH1D("monEvents", "monEvents", 20, 0., 20.);
  h = new TH1D("ntracks", "ntracks", 100, 0., 1000.);
  h = new TH1D("pt0", "pt(tracks)", 100, 0., 10.);
  h = new TH1D("pt1", "pt(tracks)", 100, 0., 50.);
  h = new TH1D("eta", "eta(tracks)", 60, -3.0, 3.0);
  h = new TH1D("phi", "phi(tracks)", 50, -3.15, 3.15);

  (void)h; // make compiler warning go away

  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->bookHist();
  }

}


// ----------------------------------------------------------------------
void bmm2Reader::readCuts(TString filename, int dump) {
  if (dump) cout << "==> bmm2Reader: Reading " << filename << " for classes setup" << endl;
  
  ifstream is(filename.Data());
  char buffer[1000]; 
  char className[200], cutFile[200]; 
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "%s %s", className, cutFile);

    // -- set up candidate analyzer classes
    if (!strcmp(className, "candAnaMuMu")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu", cutFile); 
      a->BLIND = BLIND; 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaMuMu1301")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu1301", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaMuMu1302")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu1302", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaMuMu1313")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu1313", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaBu2JpsiK")) {
      candAna *a = new candAnaBu2JpsiK(this, "candAnaBu2JpsiK", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaBs2JpsiPhi")) {
      candAna *a = new candAnaBs2JpsiPhi(this, "candAnaBs2JpsiPhi", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaDstar")) {
      candAna *a = new candAnaDstar(this, "candAnaDstar", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaHh")) {
      candAna *a = new candAnaHh(this, "candAnaHh", cutFile); 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaBd2DstarPi")) {
      candAna *a = new candAnaBd2DstarPi(this, "candAnaBd2DstarPi", cutFile); 
      lCandAnalysis.push_back(a); 
    }


    // -- all the rest ...
    if (!strcmp(className, "JSON")) {
      char json[1000];
      sscanf(buffer, "%s %s", className, json);
      JSONFILE = string(json);
      if (dump) cout << "JSON FILE:           " << JSONFILE << endl;
    }

//     if (!strcmp(className, "ptSgMUID")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptSgMUID = new PidTable(name); 
//       if (dump) cout << "Seagulls MUID:           " << name << endl;
//     }

//     if (!strcmp(className, "ptCbMUID")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptCbMUID = new PidTable(name); 
//       if (dump) cout << "Cowboys MUID:            " << name << endl;
//     }

//     if (!strcmp(className, "ptSgMUT1")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptSgMUT1 = new PidTable(name); 
//       if (dump) cout << "Seagulls MUT1:           " << name << endl;
//     }

//     if (!strcmp(className, "ptCbMUT1")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptCbMUT1 = new PidTable(name); 
//       if (dump) cout << "Cowboys MUT1:            " << name << endl;
//     }

//     if (!strcmp(className, "ptSgMUT2")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptSgMUT2 = new PidTable(name); 
//       if (dump) cout << "Seagulls MUT2:           " << name << endl;
//     }

//     if (!strcmp(className, "ptCbMUT2")) {
//       char name[1000];
//       sscanf(buffer, "%s %s", className, name);
//       ptCbMUT2 = new PidTable(name); 
//       if (dump) cout << "Cowboys MUT2:            " << name << endl;
//     }

    if (!strcmp(className, "muScaleCorrector")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      TString tname(name);  
      msc = new MuScleFitCorrector(tname);
      if (dump) cout << "MuScaleCorrector:            " << name << endl;
    }


  }


  // -- in case the muon scale corrections were not initialized: 
  if (0 == msc) {
    TString tname;  
    if (2011 == fYear) {
      if (fIsMC) {
	tname = "MuScleFitCorrector_v3/MuScleFit_2011_MC_44X.txt";
      } else {
	tname = "MuScleFitCorrector_v3/MuScleFit_2011_DATA_44X.txt";
      }
    }
    if (2012 == fYear) { 
      if (fIsMC) {
	tname = "MuScleFitCorrector_v3/MuScleFit_2012_MC_52X.txt";
      } else {
	tname = "MuScleFitCorrector_v3/MuScleFit_2012_DATA_53X.txt";
      }
    }
    cout << "==> Initializing MuScleFitCorrector with " << tname << endl;
    msc = new  MuScleFitCorrector(tname);
  }  
}


// ----------------------------------------------------------------------
void bmm2Reader::processType() {

  TGenCand *pG;
  
  // documentation line partons (entries { d, u, s, c, b, t } )
  double docPartCnt[6];
  double docAntiCnt[6];

  // partons
  double parPartCnt[6];
  double parAntiCnt[6];    
    
  for (int i = 0; i < 6; i++) {
    docPartCnt[i] = 0; 
    docAntiCnt[i] = 0; 
    parPartCnt[i] = 0; 
    parAntiCnt[i] = 0; 
  }

  int aid(0);
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {

    pG = fpEvt->getGenCand(i);

    aid = TMath::Abs(pG->fID); 
    if ( aid == 1 || aid == 2 ||
         aid == 3 || aid == 4 || 
         aid == 5 || aid == 6 || 
         aid == 21) {
      if ( pG->fStatus == 3 ) {
        //      cout << "quark/gluon from documentation #" << i << "(ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 2 &&  TMath::Abs(pG->fID) != 21) {
        //      cout << "decayed quark/gluon #" << i << " (ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 1 ) {
        //      cout << "undecayed (?) quark/gluon #" << i << " (ID: " << pG->fID  << ")" << endl;
      }
    }

    for (int j = 0; j < 6; j++) {

      if ( pG->fStatus == 3 ) {
        if ( pG->fID == j+1 ) {  
          docPartCnt[j]++;
        }
        if ( pG->fID == -(j+1) ) {  
          docAntiCnt[j]++;
        }
      }

      if ( pG->fStatus == 2 ) {
        if ( pG->fID == j+1 ) {  
          parPartCnt[j]++;
        }
        if ( pG->fID == -(j+1) ) {  
          parAntiCnt[j]++;
        }
      }
    }
  }

  fProcessType = -99;
  // -- top 
  if (docPartCnt[5] >= 1 && docAntiCnt[5] >= 1) {
    fProcessType = 50; 
    //    printf("====> t: GGF (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[5] >= 1 && docAntiCnt[5] == 0) || (docPartCnt[5] == 0 && docAntiCnt[5] >= 1) ) {
    fProcessType = 51; 
    //    printf("====> t: FEX (%i)\n", fProcessType);
    return;
  }

  if (docPartCnt[5] == 0 && docAntiCnt[5] == 0 && (parPartCnt[5] >= 1 || parAntiCnt[5] >= 1)) {
    fProcessType = 52;
    //    printf("====> t: GSP (%i)\n", fProcessType); 
    return;
  }
  
  // -- beauty
  if (docPartCnt[4] >= 1 && docAntiCnt[4] >= 1) {
    fProcessType = 40; 
   //    printf("====> b: GGF (%i)\n", fProcessType);
    return;
  } 
  
  if ((docPartCnt[4] >= 1 && docAntiCnt[4] == 0) || (docPartCnt[4] == 0 && docAntiCnt[4] >= 1) ) {
    fProcessType = 41; 
    //    printf("====> b: FEX (%i)\n", fProcessType);
    return;
  }

  if (docPartCnt[4] == 0 && docAntiCnt[4] == 0 && (parPartCnt[4] >= 1 || parAntiCnt[4] >= 1)) {
    fProcessType = 42; 
    //    printf("====> b: GSP (%i)\n", fProcessType);

    return;
  }

  if (docPartCnt[3] >= 1 && docAntiCnt[3] >= 1) {
    fProcessType = 30; 
    //    printf("====> c: GGF (%i)\n", fProcessType);
    return;
  }
  
 
  if ((docPartCnt[3] >= 1 && docAntiCnt[3] == 0) || (docPartCnt[3] == 0 && docAntiCnt[3] >= 1) ) {
    fProcessType = 31; 
    //    printf("====> c: FEX (%i)\n", fProcessType);
    return;
  }
  
  if (docPartCnt[3] == 0 && docAntiCnt[3] == 0 && (parPartCnt[3] >= 1 || parAntiCnt[3] >= 1)) {
    fProcessType = 32; 
    //    printf("====> c: GSP (%i)\n", fProcessType);
    return;
  }
  
  // light flavors
  if ((docPartCnt[5] == 0 && docAntiCnt[5] == 0) && (parPartCnt[5] == 0 && parAntiCnt[5] == 0)
      && (docPartCnt[4] == 0 && docAntiCnt[4] == 0) && (parPartCnt[4] == 0 && parAntiCnt[4] == 0)
      && (docPartCnt[3] == 0 && docAntiCnt[3] == 0) && (parPartCnt[3] == 0 && parAntiCnt[3] == 0)
      ) {
    fProcessType = 1; 
    //    printf("====> UDS: light flavors (%i)\n", fProcessType);
    return;
  }

  // if no process type was determined
  //  printf("====> Could not determine process type !!!\n");

  fpEvt->dumpGenBlock();     

}
