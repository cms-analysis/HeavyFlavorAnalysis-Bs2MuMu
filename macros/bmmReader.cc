#include "bmmReader.hh"
#include "TRandom.h"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using std::string;
using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/bmmReader.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
bmmReader::bmmReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> bmmReader: constructor..." << endl;
  cout << "==> Defining analysis cuts" << endl;
  MASSMIN = 4.5;
  MASSMAX = 6.5; 
  fVerbose = 0; 

  cout << "reading events from  " << "evts" << endl;
  char  buffer[200];
  ifstream is("evts");
  int event(0); 
  while (is.getline(buffer, 200, '\n')) {
    event = atoi(buffer); 
    cout << event << endl;
    fEventVector.push_back(event);
  }
  is.close();
  
}


// ----------------------------------------------------------------------
bmmReader::~bmmReader() {
  cout << "==> bmmReader: destructor..." << endl;
}

// ----------------------------------------------------------------------
bool bmmReader::evtFoundInCN(int evt) {
  for (unsigned int i = 0; i < fEventVector.size(); ++i) {
    if (fEventVector[i] == evt) {
      return true;
    } 
  }
  return false; 
}

// ----------------------------------------------------------------------
void bmmReader::startAnalysis() {
  cout << "==>bmmReader: setup PidTables" << endl;

  // -- Note that the return value is -99 for ranges not covered by the PidTables!
  fpMuonID = new PidTable("pidtables/jpsi/data/muid0.both.dat"); 
  fpMuonTr = new PidTable("pidtables/jpsi/data/mutrig0.both.dat"); 

  fpJSON = new JSON(JSONFILE.c_str()); 
  cout << "==> bmmReader: Starting analysis..." << endl;

}


// ----------------------------------------------------------------------
void bmmReader::eventProcessing() {

  if (fVerbose > 1) cout << "event: " << fEvent << endl;
  
//   if (218505 == fEvt) {
//     fVerbose = 100; 
//   } else {
//     fVerbose = 0; 
//   }
  //  fVerbose = 100; 
  
  //   cout << "------------------------" 
  //        << fRun << " " << fEvt
  //        << endl;
  //  if (fRun == 1 && fEvt == 218505) fpEvt->dump();

  // -- initialize all variables
  initVariables(); 

  if (fIsMC) {
    processType();
    genMatch(); 
    recoMatch(); 
    candMatch(); 
  }

  if (fIsMC) MCKinematics();
  L1TSelection();
  HLTSelection();
  trackSelection(); 
  muonSelection();

  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(0); 
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(1, fpEvt->nCands()); 
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(2, fCands.size()); 

  if (0 == IGNORETRIGGER) {
    if (!fIsMC && !fGoodHLT) {
      if (fVerbose > 5) {
	cout << "no HLT trigger!" << endl;
      }
      return;
    }
  }

  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(3, fpEvt->nCands()); 
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(4, fCands.size()); 

  ((TH1D*)fpHistFile->Get("monAllCands"))->Fill(fpEvt->nCands()); 
  ((TH1D*)fpHistFile->Get("monTypeCands"))->Fill(fCands.size()); 

  if (SELMODE < 0) {
    // -- Fill ALL candidates
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      fpCand = fCands[iC];
      fCandIdx = iC; 
      if (fVerbose > 4) cout << "Filling candidate " << iC 
			     << " with sig tracks " << fpCand->fSig1 << ".." << fpCand->fSig2
			     << endl;
      fillCandidateVariables();
      fillCandidateHistograms();
    }
  } else {
    candidateSelection(SELMODE); 
    // -- Fill only 'best' candidate
    if (fVerbose > 4) cout << "Filling candidate " << fCandIdx << endl;
    fillCandidateVariables();
    fillCandidateHistograms();
  }

  if (fIsMC) efficiencyCalculation();

  //  studyL1T();
  
  fillHist();

}


// ----------------------------------------------------------------------
void bmmReader::initVariables() {

  // -- Fill candidates vector
  fCands.clear();
  fGoodMCKinematics = true; 
  fGoodL1T = fGoodHLT = true; 

  fvGoodMuonsID.clear();
  fvGoodMuonsPt.clear();
  fvGoodMuonsEta.clear();
  fvGoodTracks.clear();
  fvGoodTracksPt.clear();
  fvGoodTracksEta.clear();

  fvGoodCand.clear(); 
  fvGoodCandPt.clear();

  fGoodEvent = true;
  fPreselection = false; 
  
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) {
      if (fVerbose > 19) cout << "Skipping candidate at " << iC << " which is of type " << pCand->fType << endl;
      continue;
    }
    if (fVerbose > 2) cout << "Adding candidate at " << iC << " which is of type " << TYPE 
			   << " with sig tracks: " << pCand->fSig1 << " .. " << pCand->fSig2
			   << " and rec tracks: " 
			   << fpEvt->getSigTrack(pCand->fSig1)->fIndex
			   << " .. " 
			   << fpEvt->getSigTrack(pCand->fSig2)->fIndex
			   << endl;
    insertCand(pCand);
  }


  // -- clear reduced tree variables
  fCandPt = fCandM = -1.;
  fGenBpartial = 0; 
}


// ----------------------------------------------------------------------
void bmmReader::insertCand(TAnaCand* pCand) {
  if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) return;

  //  cout << "bmmReader::insertCand" << endl;
   fCands.push_back(pCand); 

  // -- this is just to initialize the vector variables
  fvGoodMuonsID.push_back(true); 
  fvGoodMuonsPt.push_back(true);
  fvGoodMuonsEta.push_back(true);
  fvGoodTracks.push_back(true);
  fvGoodTracksPt.push_back(true);
  fvGoodTracksEta.push_back(true);
  fvGoodCandPt.push_back(true);
  
  fvGoodCand.push_back(false); 
}



// ----------------------------------------------------------------------
void bmmReader::MCKinematics() {

}

// ----------------------------------------------------------------------
void bmmReader::L1TSelection() {
  fGoodL1T = false; 
  fL1TMu0 = fL1TMu3 = false; 
  TString a; 
  int ps(0); 
  bool result(false), error(false); 
  if (0) cout << " ----------------------------------------------------------------------" << endl;
  for (int i = 0; i < NL1T; ++i) {
    result = error = false;
    a = fpEvt->fL1TNames[i]; 
    ps = fpEvt->fL1TPrescale[i]; 
    result = fpEvt->fL1TResult[i]; 
    error  = fpEvt->fL1TError[i]; 
    if (0 && a.Contains("Mu")) { cout << a << endl; }
    for (unsigned int j = 0; j < L1TPath.size(); ++j) {
      if (a.Contains(L1TPath[j].c_str())) {
	if (result) {
	  if (!a.CompareTo("L1_DoubleMuOpen")) fL1TMu0 = true; 
	  if (!a.CompareTo("L1_DoubleMu3")) fL1TMu3 = true; 
	  fGoodL1T = true; 
	}
      }
    }

  }

  
}


// ----------------------------------------------------------------------
void bmmReader::studyL1T() {
  TTrgObj *p; 

  ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(1); 

  if (fL1TMu0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(10); 
  if (fL1TMu3) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(11); 

  if (fGoodTracks && fGoodTracksEta && fGoodMuonsID && fGoodMuonsEta) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(20); 
    if (fL1TMu0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(21); 
    if (fL1TMu3) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(22); 
  }

  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt &&  fGoodMuonsEta) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(30); 
    if (fL1TMu0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(31); 
    if (fL1TMu3) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(32); 
  }

  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt && fGoodMuonsEta &&fGoodPt ) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(40); 
    if (fL1TMu0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(41); 
    if (fL1TMu3) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(42); 
  }



  TTrgObj *pM1(0), *pM2(0), *pM(0); 

  for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
    p = fpEvt->getTrgObj(i); 
    //    cout << "= " << p->fLabel << endl;
    if (!p->fLabel.CompareTo("hltL1sL1DoubleMuOpen:HLT::")) {
      if (0 == pM1) {
	pM1 = p; 
      } else {
	pM2 = p; 
      }
      if (0) cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
    }

    if (0 && !p->fLabel.CompareTo("hltSingleMu3L3Filtered3:HLT::")) {
      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
    }

    if (0 && !p->fLabel.CompareTo("hltMu0L2Mu0L3Filtered0:HLT::")) {
      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
    }
  }

  if (0 == pM1 || 0 == pM2) return;
  if (pM1->fP.Perp() < pM2->fP.Perp()) {
    pM = pM1; 
    pM1 = pM2; 
    pM2 = pM; 
  }

  
  if (0) cout << "===> " << pM1->fLabel << " pT = " << pM1->fP.Perp() << " eta = " << pM1->fP.Eta() << " phi = " << pM1->fP.Phi() <<  endl;
  if (0) cout << "===> " << pM2->fLabel << " pT = " << pM2->fP.Perp() << " eta = " << pM2->fP.Eta() << " phi = " << pM2->fP.Phi() <<  endl;

  bool mmEta0(false), mmEta1(false);
  bool mmL1Eta0(false), mmL1Eta1(false);
  // -- L1T muons
  if ((TMath::Abs(pM1->fP.Eta()) < 2.4) && (TMath::Abs(pM2->fP.Eta()) < 2.4)) mmL1Eta0 = true; 
  if ((TMath::Abs(pM1->fP.Eta()) < 2.0) && (TMath::Abs(pM2->fP.Eta()) < 2.0)) mmL1Eta1 = true; 
  // -- candidate muons
  if ((TMath::Abs(fMu1Eta) < 2.4) && (TMath::Abs(fMu2Eta) < 2.4)) mmEta0 = true; 
  if ((TMath::Abs(fMu1Eta) < 2.0) && (TMath::Abs(fMu2Eta) < 2.0)) mmEta1 = true; 

  ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(50); 
  if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(51); 
  if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(52); 
  if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(53); 
  if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(54); 

  if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(55); 
  if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(56); 
  if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(57); 
  if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(58); 
  

  if (fGoodTracks && fGoodTracksEta && fGoodMuonsID && fGoodMuonsEta) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(60); 
    if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(61); 
    if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(62); 
    if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(63); 
    if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(64); 

    if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(65); 
    if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(66); 
    if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(67); 
    if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(68); 
    
  }

  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt &&  fGoodMuonsEta) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(70); 
    if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(71); 
    if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(72); 
    if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(73); 
    if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(74); 

    if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(75); 
    if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(76); 
    if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(77); 
    if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(78); 

  }

  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt && fGoodMuonsEta &&fGoodPt ) {
    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(80); 
    if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(81); 
    if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(82); 
    if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(83); 
    if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(84); 

    if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(85); 
    if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(86); 
    if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(87); 
    if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(88); 
  }

  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt && fGoodMuonsEta 
      && (fCandFLS3d > 3) && (fCandCosA > 0.95) 
      ) {

    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(90); 
    if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(91); 
    if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(92); 
    if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(93); 
    if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(94); 

    if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(95); 
    if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(96); 
    if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(97); 
    if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(98); 
    
  }


  if (fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID &&  fGoodMuonsPt && fGoodMuonsEta 
      && (fCandFLS3d > 3) && (fCandCosA > 0.95) 
      && fGoodPt ) {

    ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(100); 
    if (fL1TMu0 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(101); 
    if (fL1TMu0 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(102); 
    if (fL1TMu3 && mmL1Eta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(103); 
    if (fL1TMu3 && mmL1Eta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(104); 

    if (fHLTMu0 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(105); 
    if (fHLTMu0 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(106); 
    if (fHLTMu3 && mmEta0) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(107); 
    if (fHLTMu3 && mmEta1) ((TH1D*)fpHistFile->Get("l1tStudy"))->Fill(108); 
    
  }


}


// ----------------------------------------------------------------------
void bmmReader::HLTSelection() {
  fGoodHLT = false; 
  fHLTMu0 = fHLTMu3 = false; 
  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 
  //  cout << " ----------------------------------------------------------------------" << endl;
  if (HLTPath.end() != find(HLTPath.begin(), HLTPath.end(), "NOTRIGGER"))  {
    //    cout << "NOTRIGGER requested... " << endl;
    fGoodHLT = true; 
    return;
  }

  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 
    //     if (result && a.Contains("Mu")) //&& a.Contains("Bs")
    //       cout << a << ": wasRun = " << wasRun << " result: " << result << " prescale: " << ps << " error: " << error << endl;
    for (unsigned int j = 0; j < HLTPath.size(); ++j) {
      if (a.Contains(HLTPath[j].c_str())) {
	if (wasRun && result) {
	  if (!a.CompareTo(HLTPath[j].c_str())) {
	    fGoodHLT = true; 
	  }
	  if (!a.CompareTo("HLT_DoubleMu0")) {
	    fHLTMu0 = true; 
	  }
	  if (!a.CompareTo("HLT_DoubleMu0_Quarkonium_v1")) {
	    fHLTMu0 = true; 
	  }
	  if (!a.CompareTo("HLT_DoubleMu3")) {
	    fHLTMu3 = true; 
	  }
	  if (fVerbose > 0 && fGoodHLT) cout << "HLT fired " << a << " required: " << HLTPath[j].c_str() << endl;
	}
	//cout << a << ": wasRun = " << wasRun << " result: " << result << " error: " << error << endl;
      }
    }

  }

}

// ----------------------------------------------------------------------
void bmmReader::trackSelection() {
  
  TAnaTrack *pt, *ps; 
  
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    for (int it = fCands[iC]->fSig1; it <= fCands[iC]->fSig2; ++it) {
      ps = fpEvt->getSigTrack(it); 
      pt = fpEvt->getRecTrack(ps->fIndex); 

      if (TRACKQUALITY > 0 && (0 == (pt->fTrackQuality & TRACKQUALITY))) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed track quality: " << pt->fTrackQuality << endl;
	fvGoodTracks[iC] = false; 
      }

      if (TMath::Abs(pt->fTip) > TRACKTIP) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed tip: " << pt->fTip << endl;
	fvGoodTracks[iC] = false; 
      }
      
      if (TMath::Abs(pt->fLip) > TRACKLIP) { 
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed lip: " << pt->fLip << endl;          
	fvGoodTracks[iC] = false; 
      }

      if (pt->fPlab.Eta() < TRACKETALO) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed eta: " << pt->fPlab.Eta() << endl;          
	fvGoodTracksEta[iC] = false; 
      }

      if (pt->fPlab.Eta() > TRACKETAHI) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed eta: " << pt->fPlab.Eta() << endl;          
	fvGoodTracksEta[iC] = false; 
      }

      if (pt->fPlab.Perp() < TRACKPTLO) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed pt: " << pt->fPlab.Perp() << endl;          
	fvGoodTracksPt[iC] = false; 
      }
      
      if (pt->fPlab.Perp() > TRACKPTHI) {
	if (fVerbose > 5) cout << "track " << ps->fIndex << " failed pt: " << pt->fPlab.Perp() << endl;          
	fvGoodTracksPt[iC] = false; 
      }
    }
  }

}


// ----------------------------------------------------------------------
void bmmReader::muonSelection() {

  TAnaTrack *pt, *ps; 
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    for (int it = fCands[iC]->fSig1; it <= fCands[iC]->fSig2; ++it) {
      ps = fpEvt->getSigTrack(it); 
      if (TMath::Abs(ps->fMCID) != 13) continue;
      pt = fpEvt->getRecTrack(ps->fIndex); 
      if (false == muonID(pt)) fvGoodMuonsID[iC] = false; 
      if (fVerbose > 4)
	cout << "muon " << ps->fIndex
	     << hex << " with muid = " << pt->fMuID << " MUIDMASK = " << MUIDMASK << " MUIDRESULT = " << MUIDRESULT
	     << " ===>  muonID(track) = " << muonID(pt) << dec << endl;

      if (pt->fPlab.Perp() < MUPTLO) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUPTLO: " << pt->fPlab.Perp() << " (" << ps->fPlab.Perp() << ")" << endl;
	fvGoodMuonsPt[iC] = false; 
      }
      if (pt->fPlab.Perp() > MUPTHI) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUPTHI: " << pt->fPlab.Perp() << endl;          
	fvGoodMuonsPt[iC] = false; 
      }
      if (pt->fPlab.Eta() < MUETALO) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUETALO: " << pt->fPlab.Eta() << endl;          
	fvGoodMuonsEta[iC] = false; 
      }
      if (pt->fPlab.Eta() > MUETAHI) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUETAHI: " << pt->fPlab.Eta() << endl;          
	fvGoodMuonsEta[iC] = false; 
      }
    }
  }                                     
}


// ----------------------------------------------------------------------
void bmmReader::candidateSelection(int mode) {
  cout << "==>bmmReader::candidateSelection: nothing implemented" << endl;
}

// ----------------------------------------------------------------------
void bmmReader::fillCandidateVariables() {
  if (0 == fpCand) return;
  if (-1 == fCandIdx) return;

  int goodSV(0); 
  TAnaVertex *pVtx; 
  fPvN = 0; 
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    pVtx = fpEvt->getPV(i); 
    if (0 == pVtx->fStatus) {
      ++fPvN;
    } else {
      //      cout << "skipping  fake vertex" << endl;
    }
  }
  //  cout << "fPvN = " << fPvN << endl;

  if (fpCand->fPvIdx > -1 && fpCand->fPvIdx < fpEvt->nPV()) {
    goodSV = 1; 
    TAnaVertex *pv = fpEvt->getPV(fpCand->fPvIdx); 
    fPvX = pv->fPoint.X(); 
    fPvY = pv->fPoint.Y(); 
    fPvZ = pv->fPoint.Z(); 
  } else {
    fPvX = -99.;
    fPvY = -99.;
    fPvZ = -99.;
  }

  fJSON = 0; 
  if (fIsMC) {
    fJSON = 1; 
  } else {
    fJSON = fpJSON->good(fRun, fLS); 
    if (fVerbose > 1 && !fJSON) {
      cout << "JSON = 0 for run = " << fRun << " and LS = " << fLS << endl;
    }
  }
  
  if (fIsMC) {
    fCandTM    = tmCand(fpCand); 
  } else {
    fCandTM = 0; 
  }
  
  fCandType  = fpCand->fType;
  fCandPt    = fpCand->fPlab.Perp();
  fCandEta   = fpCand->fPlab.Eta();
  fCandPhi   = fpCand->fPlab.Phi();
  fCandM     = fpCand->fMass;
  fCandPvTip = fpCand->fPvTip;
  fCandPvTipE= fpCand->fPvTipE;
  fCandPvLip = fpCand->fPvLip;
  fCandPvLipE= fpCand->fPvLipE;

  TAnaTrack *p0; 
  TAnaTrack *p1(0), *ps1(0);
  TAnaTrack *p2(0), *ps2(0); 

  fCandQ    = 0;

  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    fCandQ += fpEvt->getRecTrack(p0->fIndex)->fQ;
    if (TMath::Abs(p0->fMCID) != 13) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }

  // -- for rare backgrounds there are no "true" muons
  if (fpCand->fType > 1000000) {
    p1 = fpEvt->getSigTrack(fpCand->fSig1); 
    p2 = fpEvt->getSigTrack(fpCand->fSig2); 
  }

  // -- switch to RecTracks!
  ps1= p1; 
  p1 = fpEvt->getRecTrack(ps1->fIndex);
  ps2= p2; 
  p2 = fpEvt->getRecTrack(ps2->fIndex);

  if (1301 == fpCand->fType ||1302 == fpCand->fType || 1313 == fpCand->fType) {
    // do nothing, these types are for efficiency/TNP studies
  } else {
    if (p1->fPlab.Perp() < p2->fPlab.Perp()) {
      p0 = p1; 
      p1 = p2; 
      p2 = p0; 
      
      p0 = ps1; 
      ps1 = ps2; 
      ps2 = p0;
    }
  }

  fMu1Id        = muonID(p1); 
  fMu1Pt        = p1->fPlab.Perp(); 
  fMu1Eta       = p1->fPlab.Eta(); 
  fMu1Phi       = p1->fPlab.Phi(); 
  fMu1PtNrf     = ps1->fPlab.Perp();
  fMu1EtaNrf    = ps1->fPlab.Eta();
  fMu1TkQuality = p1->fTrackQuality & TRACKQUALITY;
  fMu1W8Mu      = fpMuonID->effD(fMu1Pt, fMu1Eta, fMu1Phi);
  fMu1W8Tr      = fpMuonTr->effD(fMu1Pt, fMu1Eta, fMu1Phi);
  fMu1Q         = p1->fQ;
  fMu1Pix       = numberOfPixLayers(p1);
  fMu1BPix      = numberOfBPixLayers(p1);
  fMu1BPixL1    = numberOfBPixLayer1Hits(p1);

  if (fCandTM && fGenM1Tmi < 0) fpEvt->dump();

  //  cout << "  " << p1->fGenIndex << "  " << fCandTM << " fCandTmi: " << fCandTmi << " gen m " << fGenM1Tmi << " " << fGenM2Tmi << endl;

  if (fCandTM) {
    TGenCand *pg1 = fpEvt->getGenCand(p1->fGenIndex);
    fMu1PtGen     = pg1->fP.Perp();
    fMu1EtaGen    = pg1->fP.Eta();
    //cout << "bmmReader: m = " << fCandM << " from cand " << fpCand << endl;
  } else {
    fMu1PtGen     = -99.;
    fMu1EtaGen    = -99.;
  }
  
  fMu2Id        = muonID(p2); 
  fMu2Pt        = p2->fPlab.Perp(); 
  fMu2Eta       = p2->fPlab.Eta(); 
  fMu2Phi       = p2->fPlab.Phi(); 
  fMu2PtNrf     = ps2->fPlab.Perp();
  fMu2EtaNrf    = ps2->fPlab.Eta();
  fMu2TkQuality = p2->fTrackQuality & TRACKQUALITY;
  fMu2W8Mu      = fpMuonID->effD(fMu2Pt, fMu2Eta, fMu2Phi);
  fMu2W8Tr      = fpMuonTr->effD(fMu2Pt, fMu2Eta, fMu2Phi);
  fMu2Q         = p2->fQ;
  fMu2Pix       = numberOfPixLayers(p2);
  fMu2BPix      = numberOfBPixLayers(p2);
  fMu2BPixL1    = numberOfBPixLayer1Hits(p2);

  if (fCandTM) {
    TGenCand *pg2 = fpEvt->getGenCand(p2->fGenIndex);
    fMu2PtGen     = pg2->fP.Perp();
    fMu2EtaGen    = pg2->fP.Eta();
  } else {
    fMu2PtGen     = -99.;
    fMu2EtaGen    = -99.;
  }


  fCandW8Mu     = fMu1W8Mu*fMu2W8Mu;
  if (TMath::Abs(fCandW8Mu) > 1.) fCandW8Mu = 0.2; // FIXME correction for missing entries at low pT
  fCandW8Tr     = fMu1W8Tr*fMu2W8Tr;
  if (TMath::Abs(fCandW8Tr) > 1.) fCandW8Tr = 0.2; // FIXME correction for missing entries at low pT 

  // -- Add HLT muons
  TTrgObj *p(0), *h1(0), *h2(0); 
  if (fHLTMu0) {
    for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
      p = fpEvt->getTrgObj(i); 
      if (fHLTMu0 && (!p->fLabel.CompareTo("hltDiMuonL3PreFiltered0:HLT::") || !p->fLabel.CompareTo("hltDiMuonL3PreFiltered0:HLT::"))) {
	if (0 == h1) {
	  h1 = p; 
	} else {
	  h2 = p; 
	}
      }
    }

    if (0 != h1 && 0 != h2 && h2->fP.Perp() > h1->fP.Perp()) {
      p = h1; 
      h1 = h2; 
      h2 = p;
    }
  }

  if (fHLTMu0 && 0 != h1 && 0 != h2) {
    //     cout << "m1 " << p1->fPlab.Perp() << " " << p1->fPlab.Phi() << " " << p1->fPlab.Eta()  << endl;
    //     cout << "h1 " << h1->fP.Perp() << " " << h1->fP.Phi() << " " << h1->fP.Eta() << endl;
    //     cout << "m2 " << p2->fPlab.Perp() << " " << p2->fPlab.Phi() << " " << p2->fPlab.Eta()  << endl;
    //     cout << "h2 " << h2->fP.Perp() << " " << h2->fP.Phi() << " " << h2->fP.Eta() << endl;

    fHltMu1Pt  = h1->fP.Perp(); 
    fHltMu1Eta = h1->fP.Eta(); 
    fHltMu1Phi =  h1->fP.Phi(); 

    fHltMu2Pt  = h1->fP.Perp();  
    fHltMu2Eta = h1->fP.Eta(); 
    fHltMu2Phi = h1->fP.Phi(); 

    TLorentzVector lh1, lh2; 
    lh1.SetPtEtaPhiM(fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, MMUON);
    lh2.SetPtEtaPhiM(fHltMu2Pt, fHltMu2Eta, fHltMu2Phi, MMUON);
    double mReco(fpCand->fMass), mHLT((lh1+lh2).M());

    ((TH1D*)fpHistFile->Get("hltMass"))->Fill((mReco - mHLT)/mReco); 

    ((TH1D*)fpHistFile->Get("hltPt"))->Fill((p1->fPlab.Perp() - h1->fP.Perp())/p1->fPlab.Perp()); 
    ((TH1D*)fpHistFile->Get("hltPt"))->Fill((p2->fPlab.Perp() - h2->fP.Perp())/p2->fPlab.Perp()); 
    ((TH1D*)fpHistFile->Get("hltEta"))->Fill((p1->fPlab.Eta() - h1->fP.Eta())/TMath::Abs(p1->fPlab.Eta())); 
    ((TH1D*)fpHistFile->Get("hltEta"))->Fill((p2->fPlab.Eta() - h2->fP.Eta())/TMath::Abs(p2->fPlab.Eta())); 
  } else {
    fHltMu1Pt  = -99.;
    fHltMu1Eta = -99.;
    fHltMu1Phi = -99.;
    
    fHltMu2Pt  = -99.;
    fHltMu2Eta = -99.;
    fHltMu2Phi = -99.;
  }


  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0); 
  TVector3 svpv(fpCand->fVtx.fPoint - fpEvt->getPV(pvidx)->fPoint); 
  double alpha = svpv.Angle(fpCand->fPlab);
  fCandCosA  = TMath::Cos(alpha);
  fCandA     = alpha; 

  double iso = isoClassic(fpCand); 
  double iso1= isoClassicOnePv(fpCand); 
  double iso2= isoClassicWithDOCA(fpCand, 0.03); // 300um DOCA cut
  double iso3= isoClassicWithDOCA(fpCand, 0.04); // 400um DOCA cut
  double iso4= isoClassicWithDOCA(fpCand, 0.05); // 500um DOCA cut
  fCandIso   = iso; 
  fCandIso1  = iso1; 
  fCandIso2  = iso2; 
  fCandIso3  = iso3; 
  fCandIso4  = iso4; 
  fCandChi2  = fpCand->fVtx.fChi2;
  fCandDof   = fpCand->fVtx.fNdof;
  fCandProb  = fpCand->fVtx.fProb;
  fCandFL3d  =fpCand->fVtx.fD3d;
  fCandFL3dE =fpCand->fVtx.fD3dE;
  fCandFLS3d = fpCand->fVtx.fD3d/fpCand->fVtx.fD3dE; 
  fCandFLSxy = fpCand->fVtx.fDxy/fpCand->fVtx.fDxyE; 

  if (fpCand->fNstTracks.size() == 0) {
    //    cout << "HHHHEEEELLLLPPPP" << endl;
    fCandDocaTrk = 99.;
  } else {
    fCandDocaTrk = fpCand->fNstTracks[0].second.first;
  }
  // ??
  TAnaTrack *t1 = p1; 
  TAnaTrack *t2 = p2; 
  double bmu1   = TMath::Sin(fpCand->fPlab.Angle(t1->fPlab));
  double bmu2   = TMath::Sin(fpCand->fPlab.Angle(t2->fPlab));
  fMu1IP        = fpCand->fVtx.fD3d*bmu1/t1->fTip;
  fMu2IP        = fpCand->fVtx.fD3d*bmu2/t2->fTip;

  // -- fill cut variables
  fWideMass = ((fpCand->fMass > MASSMIN) && (fpCand->fMass < MASSMAX)); 

  fGoodTracks = fvGoodTracks[fCandIdx];
  fGoodTracksPt = fvGoodTracksPt[fCandIdx];
  fGoodTracksEta = fvGoodTracksEta[fCandIdx];
  fGoodMuonsID  = fvGoodMuonsID[fCandIdx];
  fGoodMuonsPt  = fvGoodMuonsPt[fCandIdx]; 
  fGoodMuonsEta  = fvGoodMuonsEta[fCandIdx]; 

  fGoodQ = (fMu1Q*fMu2Q < 0); 
  fGoodPt = (fCandPt > CANDPTLO);
  fGoodEta = ((fCandEta > CANDETALO) && (fCandEta < CANDETAHI)); 
  fGoodCosA = (fCandCosA > CANDCOSALPHA); 
  fGoodIso = (fCandIso1 > CANDISOLATION); 
  fGoodChi2 = (fCandChi2/fCandDof < CANDVTXCHI2);
  fGoodFLS =  ((fCandFLS3d > CANDFLS3D) && (fCandFLSxy > CANDFLSXY)); 
  if (TMath::IsNaN(fCandFLS3d)) fGoodFLS = false;

  fGoodDocaTrk = (fCandDocaTrk > CANDDOCATRK);
  fGoodIP = true; 

  fAnaCuts.update(); 

  fPreselection = fWideMass && fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID && fGoodMuonsPt && fGoodMuonsEta; 
  fPreselection = fPreselection && fGoodQ && (fCandPt > 4) && (fCandCosA > 0.9) && (fCandFLS3d > 2) && (fCandChi2/fCandDof < 10); 
}


// ----------------------------------------------------------------------
void bmmReader::processType() {

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

// ----------------------------------------------------------------------
int bmmReader::partialReco(TAnaCand *pCand) {
  cout << "wrong function" << endl;
  return -1;
}


// ----------------------------------------------------------------------
void bmmReader::genMatch() {
  cout << "wrong function" << endl;
}

// ----------------------------------------------------------------------
void bmmReader::recoMatch() {
  cout << "wrong function" << endl;
}

// ----------------------------------------------------------------------
void bmmReader::candMatch() {
  cout << "wrong function" << endl;
}


// ----------------------------------------------------------------------
int bmmReader::tmCand2(TAnaCand *pC) {
  cout << "wrong function" << endl;
  return 0; 
}


// ----------------------------------------------------------------------
int bmmReader::tmCand(TAnaCand *pC) {

  cout << "wrong function" << endl;
  return 0; 
}


// ----------------------------------------------------------------------
double bmmReader::isoClassic(TAnaCand *pC) {
  double iso(-1.), pt(0.), sumPt(0.), candPt(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    cIdx.push_back(pT->fIndex); 
    candPt += pT->fPlab.Perp(); 
  }
  
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  continue;
    pT = fpEvt->getRecTrack(i); 
    pt = pT->fPlab.Perp(); 
    if (pt < 0.9) continue;
    if (pT->fPlab.DeltaR(pC->fPlab) < 1.0) sumPt += pt; 
  }

  iso = candPt/(candPt + sumPt); 

  return iso; 
}


// ----------------------------------------------------------------------
double bmmReader::isoClassicOnePv(TAnaCand *pC) {
  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx; 
  int pvIdx = pC->fPvIdx;

  //   int verbose(0); 
  //   if (fRun == 148862 &&  fEvt == 1219569139) { verbose = 1; }
  //   if (fRun == 147284 &&  fEvt == 186315253) { verbose = 1; }

  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    cIdx.push_back(pT->fIndex); 
    candPtScalar += pT->fPlab.Perp(); 
    //     if (verbose) {
    //       int tIdx = fpEvt->getRecTrack(pT->fIndex)->fPvIdx;
    //       if (pvIdx != tIdx) {
    // 	cout << "Signal track pointing to PV " << tIdx << " instead of " << pvIdx << endl;
    //       }
    //     }
  }

  candPt = pC->fPlab.Perp(); 
  
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  continue;
    pT = fpEvt->getRecTrack(i); 
    //     if (verbose) {
    //       cout << "   track " << i 
    //  	   << " with pT = " << pT->fPlab.Perp()
    //  	   << " eta = " << pT->fPlab.Eta()
    //  	   << " pointing at PV " << pT->fPvIdx;
    //     }
    if (pvIdx != pT->fPvIdx) {
      //      if (verbose) cout << " skipped because of PV index mismatch" << endl;
      continue;
    }
    pt = pT->fPlab.Perp(); 
    if (pt < 0.9) {
      //      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (pT->fPlab.DeltaR(pC->fPlab) < 1.0) {
      sumPt += pt; 
      //      if (verbose) cout << endl;
    } 
    //    else {
      //      if (verbose) cout << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
    //    }
  }

  iso = candPt/(candPt + sumPt); 

  //   if (verbose) cout << "--> iso = " << candPt << " .. " << sumPt << " = " << iso << endl;
  //   if (verbose) cout << "--> iso = " << pC->fPlab.Perp() << " .. " << sumPt << " = " << pC->fPlab.Perp()/(pC->fPlab.Perp() + sumPt) << endl;
  return iso; 
}

// ----------------------------------------------------------------------
double bmmReader::isoClassicWithDOCA(TAnaCand *pC, float docaCut) {
  const double ptCut = 0.9, coneSize=1.0;
  const bool verbose=false;

  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx; 
  int pvIdx = pC->fPvIdx;

 
  //   if (fRun == 148862 &&  fEvt == 1219569139) { verbose = 1; }
  //   if (fRun == 147284 &&  fEvt == 186315253) { verbose = 1; }

  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    cIdx.push_back(pT->fIndex); 
    candPtScalar += pT->fPlab.Perp(); 
        if (verbose) {
          int tIdx = fpEvt->getRecTrack(pT->fIndex)->fPvIdx;
          if (pvIdx != tIdx) {
    	cout << "Signal track pointing to PV " << tIdx << " instead of " << pvIdx << endl;
          }
        }
  }

  candPt = pC->fPlab.Perp(); 
  
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  continue;
    pT = fpEvt->getRecTrack(i); 
        if (verbose) {
          cout << "   track " << i 
     	   << " with pT = " << pT->fPlab.Perp()
     	   << " eta = " << pT->fPlab.Eta()
     	   << " pointing at PV " << pT->fPvIdx;
        }
    if (pvIdx != pT->fPvIdx) {
      if (verbose) cout << " skipped because of PV index mismatch" << endl;
      continue;
    }
    pt = pT->fPlab.Perp(); 
    if (pt < ptCut) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (pT->fPlab.DeltaR(pC->fPlab) < coneSize) {
      sumPt += pt; 
      if (verbose) cout << endl;
    } 
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
    }
  }



  // Now consider the DOCA tracks, but only those which have PV=-1
  // DOCA of close tracks
  int nsize = pC->fNstTracks.size(); 
  if(nsize>0) {
    for(int i = 0; i<nsize;++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      double docaE = pC->fNstTracks[i].second.second;

      if(doca > docaCut) continue; // check the doca cut
      
      pT = fpEvt->getRecTrack(trkId);
      // Consider only tracks from an undefined PV
      if (pT->fPvIdx >= 0) { 
	if (verbose) cout << " doca track skipped because it has a defined  PV " << pT->fPvIdx <<endl;
	continue;
      }

      pt = pT->fPlab.Perp();  // cut on track pt
      if (pt < ptCut) {
	if (verbose) cout << " doca skipped because of pt = " << pt << endl;
	continue;
      }

      if (pT->fPlab.DeltaR(pC->fPlab) < coneSize) {
	sumPt += pt; 
	if (verbose) cout << " doaa track included "<<doca<<" "<<pt<<endl;
      } 
      else {
	if (verbose) cout << " doca track skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
      }

    } // for loop over tracks
  } // end if 

  iso = candPt/(candPt + sumPt); 

  //   if (verbose) cout << "--> iso = " << candPt << " .. " << sumPt << " = " << iso << endl;
  //   if (verbose) cout << "--> iso = " << pC->fPlab.Perp() << " .. " << sumPt << " = " << pC->fPlab.Perp()/(pC->fPlab.Perp() + sumPt) << endl;

  return iso; 
}


// ----------------------------------------------------------------------
void bmmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("b1"))->Fill(fpEvt->nRecTracks()); 
}


// ----------------------------------------------------------------------
void bmmReader::fillCandidateHistograms() {
  //  cout << "bmmReader::fillCandidateHistograms() " << endl;

  //   fAnaCuts.update(); 
  //   fAnaCuts.dumpAll(); 

  // -- only candidate histograms below
  if (0 == fpCand) return;
  
  // -- historic remnant?!
  //   for (int i = 0; i < fAnaCuts.ncuts(); ++i) {
  //     if (fAnaCuts.singleCutTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dSi", i)))->Fill(fCandM);
  //     if (fAnaCuts.cumulativeCutTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dCu", i)))->Fill(fCandM);
  //     if (fAnaCuts.allOtherCutsTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dAo", i)))->Fill(fCandM);
  //     if (fAnaCuts.nMinus1CutsTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dAo", i)))->Fill(fCandM);
  //   }

  // -- Fill distributions
  fpTracksPt->fill(fMu1Pt, fCandM);
  fpTracksPt->fill(fMu2Pt, fCandM);
  fpTracksEta->fill(fMu1Eta, fCandM);
  fpTracksEta->fill(fMu2Eta, fCandM);
  fpTracksQual->fill(fMu1TkQuality, fCandM);
  fpTracksQual->fill(fMu2TkQuality, fCandM);

  fpMuonsID->fill(1, fCandM);
  fpMuonsID->fill(1, fCandM);
  fpMuonsPt->fill(fMu1Pt, fCandM);
  fpMuonsPt->fill(fMu2Pt, fCandM);
  fpMuonsEta->fill(fMu1Eta, fCandM);
  fpMuonsEta->fill(fMu2Eta, fCandM);
  fpMuon1Eta->fill(fMu2Eta, fCandM);
  fpMuon2Eta->fill(fMu2Eta, fCandM);
  fpMuon1Pt->fill(fMu1Pt, fCandM);
  fpMuon2Pt->fill(fMu2Pt, fCandM);
  fpMuon1Eta->fill(fMu1Eta, fCandM);
  fpMuon2Eta->fill(fMu2Eta, fCandM);

  fpAllEvents->fill(1, fCandM); 
  fpHLT->fill(1, fCandM); 
  fpPvZ->fill(fPvZ, fCandM); 
  fpPvN->fill(fPvN, fCandM); 
  fpQ->fill(fCandQ, fCandM); 
  fpPt->fill(fCandPt, fCandM); 
  fpEta->fill(fCandEta, fCandM); 
  fpAlpha->fill(fCandA, fCandM);
  fpCosA->fill(fCandCosA, fCandM);
  fpCosA0->fill(fCandCosA, fCandM);
  fpIso->fill(fCandIso, fCandM);
  fpIso1->fill(fCandIso1, fCandM);
  fpIso2->fill(fCandIso2, fCandM);
  fpIso3->fill(fCandIso3, fCandM);
  fpIso4->fill(fCandIso4, fCandM);

  // -- N(PV) dependent
  if (fPvN > 0 && fPvN <= 2) {
    fpIsoPv1->fill(fCandIso, fCandM);
    fpIso1Pv1->fill(fCandIso1, fCandM);
    fpIso4Pv1->fill(fCandIso4, fCandM);
    fpFLS3dPv1->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv1->fill(fCandFLSxy, fCandM); 
    fpAlphaPv1->fill(fCandA, fCandM); 
  } else if (fPvN > 2 && fPvN <= 4) {
    fpIsoPv2->fill(fCandIso, fCandM);
    fpIso1Pv2->fill(fCandIso1, fCandM);
    fpIso4Pv2->fill(fCandIso4, fCandM);
    fpFLS3dPv2->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv2->fill(fCandFLSxy, fCandM); 
    fpAlphaPv2->fill(fCandA, fCandM); 
  } else if (fPvN > 4 && fPvN <= 6) {
    fpIsoPv3->fill(fCandIso, fCandM);
    fpIso1Pv3->fill(fCandIso1, fCandM);
    fpIso4Pv3->fill(fCandIso4, fCandM);
    fpFLS3dPv3->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv3->fill(fCandFLSxy, fCandM); 
    fpAlphaPv3->fill(fCandA, fCandM); 
  } else if  (fPvN > 6 && fPvN <= 8) {
    fpIsoPv4->fill(fCandIso, fCandM);
    fpIso1Pv4->fill(fCandIso1, fCandM);
    fpIso4Pv4->fill(fCandIso4, fCandM);
    fpFLS3dPv4->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv4->fill(fCandFLSxy, fCandM); 
    fpAlphaPv4->fill(fCandA, fCandM); 
  } else if  (fPvN > 8 && fPvN <= 10) {
    fpIsoPv5->fill(fCandIso, fCandM);
    fpIso1Pv5->fill(fCandIso1, fCandM);
    fpIso4Pv5->fill(fCandIso4, fCandM);
    fpFLS3dPv5->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv5->fill(fCandFLSxy, fCandM); 
    fpAlphaPv5->fill(fCandA, fCandM); 
  } else if  (fPvN > 10 ) {
    fpIsoPv6->fill(fCandIso, fCandM);
    fpIso1Pv6->fill(fCandIso1, fCandM);
    fpIso4Pv6->fill(fCandIso4, fCandM);
    fpFLS3dPv6->fill(fCandFLS3d, fCandM); 
    fpFLSxyPv6->fill(fCandFLSxy, fCandM); 
    fpAlphaPv6->fill(fCandA, fCandM); 
  }

  fpChi2->fill(fCandChi2, fCandM);
  fpChi2Dof->fill(fCandChi2/fCandDof, fCandM); 
  fpProb->fill(fCandProb, fCandM);   
  fpFLS3d->fill(fCandFLS3d, fCandM); 
  fpFL3d->fill(fCandFL3d, fCandM); 
  fpFL3dE->fill(fCandFL3dE, fCandM); 
  fpFLSxy->fill(fCandFLSxy, fCandM); 
  fpDocaTrk->fill(fCandDocaTrk, fCandM); 
  fpIP1->fill(fMu1IP, fCandM); 
  fpIP2->fill(fMu2IP, fCandM); 

  fTree->Fill(); 


}



// ---------------------------------------------------------------------- 
void bmmReader::bookHist() {
  cout << "==> bmmReader: bookHist " << endl;

  TH1D *h; 

  h = new TH1D("l1tStudy", "l1t study", 120, 0., 120.); 
  h = new TH1D("genStudy", "gen study", 100, 0., 100.); 
  h = new TH1D("acceptance", "acceptance", 100, 0., 100.); 
  h = new TH1D("presel", "presel", 100, 0., 100.); 
  h = new TH1D("efficiency", "efficiency", 100, 0., 100.); 
  h = new TH1D("effMass", "", 100, 4., 6.0); 
  h = new TH1D("effMass1", "", 100, 0., 20.0); 
  h = new TH1D("effMass2", "", 100, 0., 20.0); 
  h = new TH1D("effMass3", "", 100, 0., 20.0); 
  h = new TH1D("effMass4", "", 100, 0., 20.0); 

  h = new TH1D("hltMass", "hlt mass resolution", 100, -0.1, 0.1); 
  h = new TH1D("hltPt", "hlt pt resolution", 100, -0.1, 0.1); 
  h = new TH1D("hltEta", "hlt eta resolution", 100, -0.1, 0.1); 

 
  // -- mass histograms for efficiencies
  //   for (int i = 0; i < fAnaCuts.ncuts(); ++i) {
  //     h = new TH1D(Form("c%dSi", i), fAnaCuts.getName(i), 40, MASSMIN, MASSMAX); 
  //     h = new TH1D(Form("c%dAo", i), fAnaCuts.getName(i), 40, MASSMIN, MASSMAX); 
  //     h = new TH1D(Form("c%dCu", i), fAnaCuts.getName(i), 40, MASSMIN, MASSMAX); 
  //     h = new TH1D(Form("c%dNm", i), fAnaCuts.getName(i), 40, MASSMIN, MASSMAX); 
  //   }

  h = new TH1D("analysisDistributions", "analysisDistributions", 100, 0., 100.); 
  fpAllEvents= bookDistribution("allevents", "allevents", "fWideMass", 10, 0., 10.);           
  fpHLT      = bookDistribution("hlt", "hlt", "fGoodHLT", 10, 0., 10.);           
  fpPvZ      = bookDistribution("pvz", "z_{PV} [cm]", "fGoodHLT", 50, -25., 25.);           
  fpPvN      = bookDistribution("pvn", "N(PV) ", "fGoodHLT", 20, 0., 20.);           
  fpTracksQual= bookDistribution("tracksqual", "Qual(tracks)", "fGoodTracks", 20, -10., 10.);
  fpTracksPt = bookDistribution("trackspt", "p_{T} [GeV]", "fGoodTracksPt", 25, 0., 25.);
  fpTracksEta= bookDistribution("trackseta", "#eta_{T}", "fGoodTracksEta", 25, -2.5, 2.5);
  fpMuonsID  = bookDistribution("muonsid", "muon id", "fGoodMuonsID", 10, 0., 10.); 
  fpMuonsPt  = bookDistribution("muonspt", "p_{T, #mu} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
  fpMuon1Pt  = bookDistribution("muon1pt", "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
  fpMuon2Pt  = bookDistribution("muon2pt", "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
  fpMuonsEta = bookDistribution("muonseta", "#eta_{#mu}", "fGoodMuonsEta", 20, -2.5, 2.5); 
  fpMuon1Eta = bookDistribution("muon1eta", "#eta_{#mu1}", "fGoodMuonsEta", 20, -2.5, 2.5); 
  fpMuon2Eta = bookDistribution("muon2eta", "#eta_{#mu2}", "fGoodMuonsEta", 20, -2.5, 2.5); 
  fpQ        = bookDistribution("q", "q_{1} q_{2}", "fGoodQ", 3, -1., 2.); 
  fpPt       = bookDistribution("pt", "p_{T, B} [GeV]", "fGoodPt", 20, 0., 40.); 
  fpEta      = bookDistribution("eta", "#eta_{B}", "fGoodEta", 20, -2.5, 2.5); 
  fpCosA     = bookDistribution("cosa", "cos(#alpha)", "fGoodCosA", 30, 0.97, 1.); 
  fpCosA0    = bookDistribution("cosa0", "cos(#alpha)", "fGoodCosA", 44, -1.1, 1.1); 
  fpAlpha    = bookDistribution("alpha", "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpIso      = bookDistribution("iso",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIso1     = bookDistribution("iso1", "I", "fGoodIso", 22, 0., 1.1); 
  fpIso2     = bookDistribution("iso2", "I2", "fGoodIso", 22, 0., 1.1); 
  fpIso3     = bookDistribution("iso3", "I3", "fGoodIso", 22, 0., 1.1); 
  fpIso4     = bookDistribution("iso4", "I4", "fGoodIso", 22, 0., 1.1); 

  fpAlphaPv1 = bookDistribution("alpha1",  "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpAlphaPv2 = bookDistribution("alpha2",  "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpAlphaPv3 = bookDistribution("alpha3",  "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpAlphaPv4 = bookDistribution("alpha4",  "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpAlphaPv5 = bookDistribution("alpha5",  "#alpha", "fGoodCosA", 20, 0., 0.2); 
  fpAlphaPv6 = bookDistribution("alpha6",  "#alpha", "fGoodCosA", 20, 0., 0.2); 

  fpIsoPv1   = bookDistribution("isopv1",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIsoPv2   = bookDistribution("isopv2",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIsoPv3   = bookDistribution("isopv3",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIsoPv4   = bookDistribution("isopv4",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIsoPv5   = bookDistribution("isopv5",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIsoPv6   = bookDistribution("isopv6",  "I (old)", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv1   = bookDistribution("iso1pv1",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv2   = bookDistribution("iso1pv2",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv3   = bookDistribution("iso1pv3",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv4   = bookDistribution("iso1pv4",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv5   = bookDistribution("iso1pv5",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso1Pv6   = bookDistribution("iso1pv6",  "I", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv1   = bookDistribution("iso4pv1",  "I4", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv2   = bookDistribution("iso4pv2",  "I4", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv3   = bookDistribution("iso4pv3",  "I4", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv4   = bookDistribution("iso4pv4",  "I4", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv5   = bookDistribution("iso4pv5",  "I4", "fGoodIso", 22, 0., 1.1); 
  fpIso4Pv6   = bookDistribution("iso4pv6",  "I4", "fGoodIso", 22, 0., 1.1); 

  fpChi2     = bookDistribution("chi2",  "#chi^{2}", "fGoodChi2", 30, 0., 30.);              
  fpChi2Dof  = bookDistribution("chi2dof",  "#chi^{2}/dof", "fGoodChi2", 30, 0., 3.);       
  fpProb     = bookDistribution("pchi2dof",  "P(#chi^{2},dof)", "fGoodChi2", 25, 0., 1.);    
  fpFLS3d    = bookDistribution("fls3d", "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.);  
  fpFL3d     = bookDistribution("fl3d",  "l_{3d}", "fGoodFLS", 25, 0., 5.);  
  fpFL3dE    = bookDistribution("fl3dE", "#sigma(l_{3d})", "fGoodFLS", 25, 0., 0.5);  
  fpFLSxy    = bookDistribution("flsxy", "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.);  
  fpDocaTrk  = bookDistribution("docatrk", "d_{ca}(track)", "fGoodDocaTrk", 35, 0., 0.14);   
  fpIP1      = bookDistribution("ip1", "IP_{1}/lsin(#beta)", "fGoodIP", 40, -4., 4.);        
  fpIP2      = bookDistribution("ip2", "IP_{2}/lsin(#beta)", "fGoodIP", 40, -4., 4.);        

  fpFLS3dPv1  = bookDistribution("fls3dpv1",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 
  fpFLS3dPv2  = bookDistribution("fls3dpv2",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 
  fpFLS3dPv3  = bookDistribution("fls3dpv3",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 
  fpFLS3dPv4  = bookDistribution("fls3dpv4",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 
  fpFLS3dPv5  = bookDistribution("fls3dpv5",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 
  fpFLS3dPv6  = bookDistribution("fls3dpv6",  "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.); 

  fpFLSxyPv1  = bookDistribution("flsxypv1",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  fpFLSxyPv2  = bookDistribution("flsxypv2",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  fpFLSxyPv3  = bookDistribution("flsxypv3",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  fpFLSxyPv4  = bookDistribution("flsxypv4",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  fpFLSxyPv5  = bookDistribution("flsxypv5",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  fpFLSxyPv6  = bookDistribution("flsxypv6",  "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.); 
  
  h = new TH1D("b1", "Ntrk", 200, 0., 200.);
  h = new TH1D("bnc0", "NCand before selection", 20, 0., 20.);
  h = new TH1D("bnc1", "NCand after selection", 20, 0., 20.);
  h = new TH1D("monEvents", "monEvents", 10, 0., 10.);
  h = new TH1D("monAllCands", "monAllCands", 1000, 0., 1000.);
  h = new TH1D("monTypeCands", "monTypeCands", 1000, 0., 1000.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun,               "run/I");
  fTree->Branch("json",   &fJSON,              "json/O");
  fTree->Branch("evt",    &fEvt,               "evt/I");
  fTree->Branch("ls",     &fLS,                "ls/I");
  fTree->Branch("mck",    &fGoodMCKinematics,  "mck/O");
  fTree->Branch("tm",     &fCandTM,            "tm/I");
  fTree->Branch("pr",     &fGenBpartial,       "pr/I"); 
  fTree->Branch("procid", &fProcessType,       "procid/I");
  fTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fTree->Branch("l1t",    &fGoodL1T,           "l1t/O");
  fTree->Branch("pvn",    &fPvN,               "pvn/I");

  // -- global cuts and weights
  fTree->Branch("gmuid",  &fGoodMuonsID,       "gmuid/O");
  fTree->Branch("gmupt",  &fGoodMuonsPt,       "gmupt/O");
  fTree->Branch("gmueta", &fGoodMuonsEta,      "gmueta/O");
  fTree->Branch("gtqual", &fGoodTracks,        "gtqual/O");
  fTree->Branch("gtpt",   &fGoodTracksPt,      "gtpt/O");
  fTree->Branch("gteta",  &fGoodTracksEta,     "gteta/O");
  fTree->Branch("w8mu",   &fCandW8Mu,          "w8mu/D");
  fTree->Branch("w8tr",   &fCandW8Tr,          "w8tr/D");
  // -- cand
  fTree->Branch("q",      &fCandQ,             "q/I");
  fTree->Branch("type",   &fCandType,          "type/I");
  fTree->Branch("pt",     &fCandPt,            "pt/D");
  fTree->Branch("eta",    &fCandEta,           "eta/D");
  fTree->Branch("phi",    &fCandPhi,           "phi/D");
  fTree->Branch("m",      &fCandM,             "m/D");
  fTree->Branch("cosa",   &fCandCosA,          "cosa/D");
  fTree->Branch("alpha",  &fCandA,             "alpha/D");
  fTree->Branch("iso",    &fCandIso,           "iso/D");
  fTree->Branch("iso1",   &fCandIso1,          "iso1/D");
  fTree->Branch("iso2",   &fCandIso2,          "iso2/D");
  fTree->Branch("iso3",   &fCandIso3,          "iso3/D");
  fTree->Branch("iso4",   &fCandIso4,          "iso4/D");
  fTree->Branch("chi2",   &fCandChi2,          "chi2/D");
  fTree->Branch("dof",    &fCandDof,           "dof/D");
  fTree->Branch("fls3d",  &fCandFLS3d,         "fls3d/D");
  fTree->Branch("fl3d",   &fCandFL3d,          "fl3d/D");
  fTree->Branch("fl3dE",  &fCandFL3dE,         "fl3dE/D");
  fTree->Branch("flsxy",  &fCandFLSxy,         "flsxy/D");
  fTree->Branch("docatrk",&fCandDocaTrk,       "docatrk/D");
  fTree->Branch("lip",    &fCandPvLip,         "lip/D");
  fTree->Branch("lipE",   &fCandPvLipE,        "lipE/D");
  fTree->Branch("tip",    &fCandPvTip,         "tip/D");
  fTree->Branch("tipE",   &fCandPvTipE,        "tipE/D");
  // -- muons
  fTree->Branch("m1q",     &fMu1Q,              "m1q/I");
  fTree->Branch("m1id",    &fMu1Id,             "m1id/O");
  fTree->Branch("m1pt",    &fMu1Pt,             "m1pt/D");
  fTree->Branch("m1eta",   &fMu1Eta,            "m1eta/D");
  fTree->Branch("m1phi",   &fMu1Phi,            "m1phi/D");
  fTree->Branch("m1ip",    &fMu1IP,             "m1ip/D");
  fTree->Branch("m1gt",    &fMu1TkQuality,      "m1gt/I");
  fTree->Branch("m1pix",   &fMu1Pix,            "m1pix/I");
  fTree->Branch("m1bpix",  &fMu1BPix,           "m1bpix/I");
  fTree->Branch("m1bpixl1",&fMu1BPixL1,         "m1bpixl1/I");
  fTree->Branch("m2q",     &fMu2Q,              "m2q/I");
  fTree->Branch("m2id",    &fMu2Id,             "m2id/O");
  fTree->Branch("m2pt",    &fMu2Pt,             "m2pt/D");
  fTree->Branch("m2eta",   &fMu2Eta,            "m2eta/D");
  fTree->Branch("m2phi",   &fMu2Phi,            "m2phi/D");
  fTree->Branch("m2ip",    &fMu2IP,             "m2ip/D");
  fTree->Branch("m2gt",    &fMu2TkQuality,      "m2gt/I");
  fTree->Branch("m2pix",   &fMu2Pix,            "m2pix/I");
  fTree->Branch("m2bpix",  &fMu2BPix,           "m2bpix/I");
  fTree->Branch("m2bpixl1",&fMu2BPixL1,         "m2bpixl1/I");

  fTree->Branch("g1pt",    &fMu1PtGen,          "g1pt/D");
  fTree->Branch("g2pt",    &fMu2PtGen,          "g2pt/D");
  fTree->Branch("g1eta",   &fMu1EtaGen,         "g1eta/D");
  fTree->Branch("g2eta",   &fMu2EtaGen,         "g2eta/D");

  fTree->Branch("t1pt",    &fMu1PtNrf,          "t1pt/D");
  fTree->Branch("t1eta",   &fMu1EtaNrf,         "t1eta/D");
  fTree->Branch("t2pt",    &fMu2PtNrf,          "t2pt/D");
  fTree->Branch("t2eta",   &fMu2EtaNrf,         "t2eta/D");

  fTree->Branch("fHltMu1Pt",  &fHltMu1Pt,  "hm1pt/D");    
  fTree->Branch("fHltMu1Eta", &fHltMu1Eta, "hm1eta/D");  
  fTree->Branch("fHltMu1Phi", &fHltMu1Phi, "hm1phi/D");  
  fTree->Branch("fHltMu2Pt",  &fHltMu2Pt,  "hm2pt/D");    
  fTree->Branch("fHltMu2Eta", &fHltMu2Eta, "hm2eta/D");  
  fTree->Branch("fHltMu2Phi", &fHltMu2Phi, "hm2phi/D");  


  // -- Efficiency/Acceptance Tree
  fEffTree = new TTree("effTree", "effTree");
  fEffTree->Branch("run",    &fRun,               "run/I");
  fEffTree->Branch("evt",    &fEvt,               "evt/I");
  fEffTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fEffTree->Branch("procid", &fProcessType,       "procid/I");

  fEffTree->Branch("gpt",    &fETgpt,             "gpt/F");
  fEffTree->Branch("geta",   &fETgeta,            "geta/F");

  fEffTree->Branch("m1pt",   &fETm1pt,            "m1pt/F");
  fEffTree->Branch("g1pt",   &fETg1pt,            "g1pt/F");
  fEffTree->Branch("m1eta",  &fETm1eta,           "m1eta/F");
  fEffTree->Branch("g1eta",  &fETg1eta,           "g1eta/F");
  fEffTree->Branch("m1q",    &fETm1q,             "m1q/I");
  fEffTree->Branch("m1gt",   &fETm1gt,            "m1gt/O");
  fEffTree->Branch("m1id",   &fETm1id,            "m1id/O");

  fEffTree->Branch("m2pt",   &fETm2pt,            "m2pt/F");
  fEffTree->Branch("g2pt",   &fETg2pt,            "g2pt/F");
  fEffTree->Branch("m2eta",  &fETm2eta,           "m2eta/F");
  fEffTree->Branch("g2eta",  &fETg2eta,           "g2eta/F");
  fEffTree->Branch("m2q",    &fETm2q,             "m2q/I");
  fEffTree->Branch("m2gt",   &fETm2gt,            "m2gt/O");
  fEffTree->Branch("m2id",   &fETm2id,            "m2id/O");

  fEffTree->Branch("m",      &fETcandMass,        "m/F");
}


// ---------------------------------------------------------------------- 
AnalysisDistribution* bmmReader::bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn, ht, nbins, lo, hi); 
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX); 
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX); 
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX); 
  p->setAnalysisCuts(&fAnaCuts, hc); 
  p->setPreselCut(&fPreselection); 

  TH1 *h = (TH1D*)gFile->Get("analysisDistributions"); 
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (!strcmp(h->GetXaxis()->GetBinLabel(i), "")) {
      // cout << "adding at bin " << i << " the label " << hn << endl;
      h->GetXaxis()->SetBinLabel(i, hn);
      break;
    }
  }
  return p; 
}


// ----------------------------------------------------------------------
void bmmReader::efficiencyCalculation() {
  cout << "bmmReader::efficiencyCalculation() missing implementation" << endl;
}


// ----------------------------------------------------------------------
void bmmReader::basicCuts() {
  cout << "    bmmReader basic cuts" << endl;
  fAnaCuts.addCut("fWideMass", "m(B candidate) [GeV]", fWideMass); 
  //  fAnaCuts.addCut("fGoodL1T", "L1T", fGoodL1T); 
  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT); 
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID); 
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt); 
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta); 
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks); 
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt); 
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta); 
}


// ----------------------------------------------------------------------
void bmmReader::moreBasicCuts() {
  cout << "    bmmReader more basic cuts?" << endl;

}


// ----------------------------------------------------------------------
void bmmReader::candidateCuts() {
  cout << "    bmmReader candidate cuts" << endl;
  fAnaCuts.addCut("fGoodQ", "q_{1} 1_{2}", fGoodQ); 
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt); 
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta); 
  fAnaCuts.addCut("fGoodCosA", "cos(#alpha)", fGoodCosA); 
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS); 
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2); 
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso); 
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk); 
  fAnaCuts.addCut("fGoodIP", "sin#beta*l_{3d}/IP", fGoodIP); 
}


// ----------------------------------------------------------------------
void bmmReader::moreCandidateCuts() {
  cout << "    bmmReader more candidate cuts?" << endl;

}


// ----------------------------------------------------------------------
void bmmReader::readCuts(TString filename, int dump) {

  // -- set up cut sequence for analysis
  basicCuts(); 
  moreBasicCuts(); 
  candidateCuts(); 
  moreCandidateCuts(); 

  
  fCutFile = filename;
  if (dump) cout << "==> bmmReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(string(fCutFile.Data()), cutLines);

  char CutName[100];
  float CutValue;
  int ok(0);

  char  buffer[200];
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.Data());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
      if (1313 == TYPE) cstring = "#mu^{+}#mu^{-}";
      if (301313 == TYPE) cstring = "#mu^{+}#mu^{-}";
      if (200521 == TYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      if (300521 == TYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      ibin = 1;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "SELMODE")) {
      SELMODE = int(CutValue); ok = 1;
      if (dump) cout << "SELMODE:           " << SELMODE << endl;
      ibin = 2;
      hcuts->SetBinContent(ibin, SELMODE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Selection Mode :: %i", CutName, SELMODE));
    }

    if (!strcmp(CutName, "JSON")) {
      char json[1000]; ok = 1; 
      sscanf(buffer, "%s %s", CutName, json);
      JSONFILE = string(json); ok = 1; 
      if (dump) cout << "JSON FILE:           " << JSONFILE << endl;
      ibin = 3;
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: JSON File :: %s", CutName, JSONFILE.c_str()));
    }

    if (!strcmp(CutName, "TRIGGER")) {
      char triggerlist[1000]; ok = 1; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      string::size_type idx = tl.find_first_of(","); 
      while (idx != string::npos) {
 	string hlt = tl.substr(0, idx); 
	// cout << "idx: " << idx << " hlt: " << hlt << endl;
 	tl.erase(0, idx+1); 
 	HLTPath.push_back(hlt); 
 	idx = tl.find_first_of(","); 
      }
      HLTPath.push_back(tl); 
      if (dump) {
	cout << "TRIGGER:          "; 
	for (unsigned int i = 0; i < HLTPath.size(); ++i) cout << HLTPath[i] << " "; cout << endl;
      }
      ibin = 4; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }

    if (!strcmp(CutName, "L1TRIGGER")) {
      char triggerlist[1000]; ok = 1; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      string::size_type idx = tl.find_first_of(","); 
      while (idx != string::npos) {
 	string hlt = tl.substr(0, idx); 
	// cout << "idx: " << idx << " hlt: " << hlt << endl;
 	tl.erase(0, idx+1); 
 	L1TPath.push_back(hlt); 
 	idx = tl.find_first_of(","); 
      }
      L1TPath.push_back(tl); 
      if (dump) {
	cout << "L1TRIGGER:          "; 
	for (unsigned int i = 0; i < L1TPath.size(); ++i) cout << L1TPath[i] << " "; cout << endl;
      }
      ibin = 5; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }

    if (!strcmp(CutName, "TRUTHCAND")) {
      TRUTHCAND = int(CutValue); ok = 1;
      if (dump) cout << "TRUTHCAND:           " << TRUTHCAND << endl;
      ibin = 6;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "IGNORETRIGGER")) {
      IGNORETRIGGER = int(CutValue); ok = 1;
      if (dump) cout << "IGNORETRIGGER      " << IGNORETRIGGER << endl;
      ibin = 7;
      hcuts->SetBinContent(ibin, IGNORETRIGGER);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore trigger :: %i", CutName, IGNORETRIGGER));
    }

    if (!strcmp(CutName, "CANDPTLO")) {
      CANDPTLO = CutValue; ok = 1;
      if (dump) cout << "CANDPTLO:           " << CANDPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, CANDPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(%s) :: %3.1f", CutName, cstring.c_str(), CANDPTLO));
    }

    if (!strcmp(CutName, "CANDETALO")) {
      CANDETALO = CutValue; ok = 1;
      if (dump) cout << "CANDETALO:           " << CANDETALO << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, CANDETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(%s) :: %3.1f", CutName, cstring.c_str(), CANDETALO));
    }

    if (!strcmp(CutName, "CANDETAHI")) {
      CANDETAHI = CutValue; ok = 1;
      if (dump) cout << "CANDETAHI:           " << CANDETAHI << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, CANDETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(%s) :: %3.1f", CutName, cstring.c_str(), CANDETAHI));
    }

    if (!strcmp(CutName, "CANDCOSALPHA")) {
      CANDCOSALPHA = CutValue; ok = 1;
      if (dump) cout << "CANDCOSALPHA:           " << CANDCOSALPHA << endl;
      ibin = 14;
      hcuts->SetBinContent(ibin, CANDCOSALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: cos#alpha :: %5.4f", CutName, CANDCOSALPHA));
    }

    if (!strcmp(CutName, "CANDFLS3D")) {
      CANDFLS3D = CutValue; ok = 1;
      if (dump) cout << "CANDFLS3D:           " << CANDFLS3D << endl;
      ibin = 15;
      hcuts->SetBinContent(ibin, CANDFLS3D);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{3d}/#sigma(l_{3d}) :: %3.1f", CutName, CANDFLS3D));
    }

    if (!strcmp(CutName, "CANDFLSXY")) {
      CANDFLSXY = CutValue; ok = 1;
      if (dump) cout << "CANDFLSXY:           " << CANDFLSXY << endl;
      ibin = 16;
      hcuts->SetBinContent(ibin, CANDFLSXY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{xy}/#sigma(l_{xy}) :: %3.1f", CutName, CANDFLSXY));
    }

    if (!strcmp(CutName, "CANDVTXCHI2")) {
      CANDVTXCHI2 = CutValue; ok = 1;
      if (dump) cout << "CANDVTXCHI2:           " << CANDVTXCHI2 << endl;
      ibin = 17;
      hcuts->SetBinContent(ibin, CANDVTXCHI2);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #chi^{2} :: %3.1f", CutName, CANDVTXCHI2));
    }

    if (!strcmp(CutName, "CANDISOLATION")) {
      CANDISOLATION = CutValue; ok = 1;
      if (dump) cout << "CANDISOLATION:           " << CANDISOLATION << endl;
      ibin = 18;
      hcuts->SetBinContent(ibin, CANDISOLATION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: I_{trk} :: %4.2f", CutName, CANDISOLATION));
    }

    if (!strcmp(CutName, "CANDDOCATRK")) {
      CANDDOCATRK = CutValue; ok = 1;
      if (dump) cout << "CANDDOCATRK:           " << CANDDOCATRK << endl;
      ibin = 19;
      hcuts->SetBinContent(ibin, CANDDOCATRK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{trk} :: %4.3f", CutName, CANDDOCATRK));
    }



    if (!strcmp(CutName, "SIGBOXMIN")) {
      SIGBOXMIN = CutValue; ok = 1;
      if (dump) cout << "SIGBOXMIN:           " << SIGBOXMIN << endl;
      ibin = 90;
      hcuts->SetBinContent(ibin, SIGBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMIN :: %6.3f", CutName, SIGBOXMIN));
    }

    if (!strcmp(CutName, "SIGBOXMAX")) {
      SIGBOXMAX = CutValue; ok = 1;
      if (dump) cout << "SIGBOXMAX:           " << SIGBOXMAX << endl;
      ibin = 91;
      hcuts->SetBinContent(ibin, SIGBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMAX :: %6.3f", CutName, SIGBOXMAX));
    }

    if (!strcmp(CutName, "BGLBOXMIN")) {
      BGLBOXMIN = CutValue; ok = 1;
      if (dump) cout << "BGLBOXMIN:           " << BGLBOXMIN << endl;
      ibin = 92;
      hcuts->SetBinContent(ibin, BGLBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMIN :: %6.3f", CutName, BGLBOXMIN));
    }

    if (!strcmp(CutName, "BGLBOXMAX")) {
      BGLBOXMAX = CutValue; ok = 1;
      if (dump) cout << "BGLBOXMAX:           " << BGLBOXMAX << endl;
      ibin = 93;
      hcuts->SetBinContent(ibin, BGLBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMAX :: %6.3f", CutName, BGLBOXMAX));
    }

    if (!strcmp(CutName, "BGHBOXMIN")) {
      BGHBOXMIN = CutValue; ok = 1;
      if (dump) cout << "BGHBOXMIN:           " << BGHBOXMIN << endl;
      ibin = 94;
      hcuts->SetBinContent(ibin, BGHBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMIN :: %6.3f", CutName, BGHBOXMIN));
    }

    if (!strcmp(CutName, "BGHBOXMAX")) {
      BGHBOXMAX = CutValue; ok = 1;
      if (dump) cout << "BGHBOXMAX:           " << BGHBOXMAX << endl;
      ibin = 95;
      hcuts->SetBinContent(ibin, BGHBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMAX :: %6.3f", CutName, BGHBOXMAX));
    }

    // -- Tracks
    if (!strcmp(CutName, "TRACKQUALITY")) {
      TRACKQUALITY = CutValue; ok = 1;
      if (dump) cout << "TRACKQUALITY:           " << TRACKQUALITY << " " << endl;
      ibin = 100;
      hcuts->SetBinContent(ibin, TRACKQUALITY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: track quality :: %d", CutName, TRACKQUALITY));
    }

    if (!strcmp(CutName, "TRACKPTLO")) {
      TRACKPTLO = CutValue; ok = 1;
      if (dump) cout << "TRACKPTLO:           " << TRACKPTLO << " GeV" << endl;
      ibin = 101;
      hcuts->SetBinContent(ibin, TRACKPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(track) :: %3.1f", CutName, TRACKPTLO));
    }

    if (!strcmp(CutName, "TRACKPTHI")) {
      TRACKPTHI = CutValue; ok = 1;
      if (dump) cout << "TRACKPTHI:           " << TRACKPTHI << " GeV" << endl;
      ibin = 102;
      hcuts->SetBinContent(ibin, TRACKPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(track) :: %3.1f", CutName, TRACKPTHI));
    }

    if (!strcmp(CutName, "TRACKTIP")) {
      TRACKTIP = CutValue; ok = 1;
      if (dump) cout << "TRACKTIP:           " << TRACKTIP << " cm" << endl;
      ibin = 103;
      hcuts->SetBinContent(ibin, TRACKTIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{xy}(track) :: %3.1f", CutName, TRACKTIP));
    }

    if (!strcmp(CutName, "TRACKLIP")) {
      TRACKLIP = CutValue; ok = 1;
      if (dump) cout << "TRACKLIP:           " << TRACKLIP << " cm" << endl;
      ibin = 104;
      hcuts->SetBinContent(ibin, TRACKLIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{z}(track) :: %3.1f", CutName, TRACKLIP));
    }

    if (!strcmp(CutName, "TRACKETALO")) {
      TRACKETALO = CutValue; ok = 1;
      if (dump) cout << "TRACKETALO:           " << TRACKETALO << " " << endl;
      ibin = 105;
      hcuts->SetBinContent(ibin, TRACKETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{min}(track) :: %3.1f", CutName, TRACKETALO));
    }

    if (!strcmp(CutName, "TRACKETAHI")) {
      TRACKETAHI = CutValue; ok = 1;
      if (dump) cout << "TRACKETAHI:           " << TRACKETAHI << " " << endl;
      ibin = 106;
      hcuts->SetBinContent(ibin, TRACKETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{max}(track) :: %3.1f", CutName, TRACKETAHI));
    }

    // -- Muons
    if (!strcmp(CutName, "MUIDMASK")) {
      MUIDMASK = int(CutValue); ok = 1;
      if (dump) cout << "MUIDMASK:           " << MUIDMASK << endl;
      ibin = 200;
      hcuts->SetBinContent(ibin, MUIDMASK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDMask :: %d", CutName, MUIDMASK));
    }

    if (!strcmp(CutName, "MUIDRESULT")) {
      // MUIDRESULT == 0: compare result of & with ">=0"
      // MUIDRESULT != 0: compare result of & with "==MUIDRESULT"
      MUIDRESULT = int(CutValue); ok = 1;
      if (dump) cout << "MUIDRESULT:           " << MUIDRESULT << endl;
      ibin = 201;
      hcuts->SetBinContent(ibin, MUIDRESULT);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDResult :: %d", CutName, MUIDRESULT));
    }

    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; ok = 1;
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 202;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(#mu) :: %3.1f", CutName, MUPTLO));
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; ok = 1;
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 203;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(#mu) :: %3.1f", CutName, MUPTHI));
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; ok = 1;
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 204;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(#mu) :: %3.1f", CutName, MUETALO));
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; ok = 1;
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 205;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(#mu) :: %3.1f", CutName, MUETAHI));
    }

    if (!strcmp(CutName, "MUIP")) {
      MUIP = CutValue; ok = 1;
      if (dump) cout << "MUIP:           " << MUIP << endl;
      ibin = 206;
      hcuts->SetBinContent(ibin, MUIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: IP(#mu) :: %3.1f", CutName, MUIP));
    }

    if (!strcmp(CutName, "JPSITYPE")) {
      JPSITYPE = CutValue; ok = 1;
      if (dump) cout << "JPSITYPE:      " << JPSITYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, JPSITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: J/#psi ID :: %d", CutName, JPSITYPE));
    }

    if (!strcmp(CutName, "JPSIMASSLO")) {
      JPSIMASSLO = CutValue; ok = 1;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSLO << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, JPSIMASSLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSLO));
    }

    if (!strcmp(CutName, "JPSIMASSHI")) {
      JPSIMASSHI = CutValue; ok = 1;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSHI << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, JPSIMASSHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSHI));
    }


    //    if (!ok) cout << "==> bmmReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}

// ----------------------------------------------------------------------
void bmmReader::readFile(string filename, vector<string> &lines) {
  cout << "readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  char input[1000]; 
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] != '+') {
      lines.push_back(string(buffer));
    } else {
      sscanf(buffer, "+input %s", input);
      readFile(input, lines); 
    }
  }

}


// ----------------------------------------------------------------------
int bmmReader::checkCut(const char *s, TH1D *hcuts) {
  TString a; 
  int ok(0); 
  for (int i = 1; i <= hcuts->GetNbinsX(); ++i) {
    a = hcuts->GetXaxis()->GetBinLabel(i); 
    if (a.Contains(s)) {
      ok = 1; 
      break;
    }
  }  
  return ok; 
}


// ----------------------------------------------------------------------
bool bmmReader::muonID(TAnaTrack *pT) {
  int result = pT->fMuID & MUIDMASK;
  if (HLTPath.end() != find(HLTPath.begin(), HLTPath.end(), "NOTRIGGER"))  {
    return true; 
  }

  if (0 == MUIDRESULT ) {
    // -- result must be larger than zero for successful muon ID, i.e., the OR of the bits
    if (0 == result){
      if (fVerbose > 4) cout << "muon " << pT->fIndex << " failed MUID: " << hex << pT->fMuID << dec << endl;          
      return false; 
    } else {
      return true; 
    }
  } else {
    // -- result must be equal to the mask for successful muon ID, i.e., the AND of the bits
    if (MUIDRESULT != result){
      if (fVerbose > 4) cout << "muon " << pT->fIndex << " failed MUID: " << hex << pT->fMuID << dec << endl;          
      return false; 
    } else {
      return true;
    }
  }
}
