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
  fAnaCuts.addCut("fWideMass", "m(B candidate) [GeV]", fWideMass); 
  fAnaCuts.addCut("fGoodL1T", "L1T", fGoodL1T); 
  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT); 
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks); 
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt); 
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta); 
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID); 
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt); 
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta); 
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
bmmReader::~bmmReader() {
  cout << "==> bmmReader: destructor..." << endl;

}

// ----------------------------------------------------------------------
void bmmReader::startAnalysis() {
  cout << "==> bmmReader: Starting analysis..." << endl;

}


// ----------------------------------------------------------------------
void bmmReader::eventProcessing() {

  if (fVerbose > 1) cout << "event: " << fEvent << endl;

  // -- initialize all variables
  initVariables(); 

  if (fIsMC) MCKinematics();
  L1TSelection();
  HLTSelection();
  trackSelection(); 
  muonSelection();

  candidateSelection(SELMODE); 

  if (SELMODE < 0) {
    // -- Fill ALL candidates
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      if (fvGoodCand[iC]) {
	fpCand = fCands[iC];
	if (fVerbose > 4) cout << "Filling candidate " << iC 
			       << " with sig tracks " << fpCand->fSig1 << ".." << fpCand->fSig2
			       << endl;
	fillCandidateVariables();
	fillCandidateHistograms();
      }
    }
  } else {
    // -- Fill only 'best' candidate
    if (fVerbose > 4) cout << "Filling candidate " << fCandIdx << endl;
    fillCandidateVariables();
    fillCandidateHistograms();
  }

  studyL1T();
  
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
  
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) {
      if (fVerbose > 19) cout << "Skipping candidate at " << iC << " which is of type " << pCand->fType << endl;
      continue;
    }
    if (fVerbose > 2) cout << "Adding candidate at " << iC << " which is of type " << TYPE 
			   << " with sig tracks: " << pCand->fSig1 << " .. " << pCand->fSig2
			   << endl;
    insertCand(pCand);
  }


  // -- clear reduced tree variables
  fCandPt = fCandM = -1.;

}


// ----------------------------------------------------------------------
void bmmReader::insertCand(TAnaCand* pCand) {
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
  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 
    //    if (result) cout << a << ": wasRun = " << wasRun << " result: " << result << " error: " << error << endl;
    for (unsigned int j = 0; j < HLTPath.size(); ++j) {
      if (a.Contains(HLTPath[j].c_str())) {
	if (wasRun && result) {
	  if (!a.CompareTo("HLT_DoubleMu0")) fHLTMu0 = true; 
	  if (!a.CompareTo("HLT_DoubleMu0_Quarkonium_v1")) fHLTMu0 = true; 
	  if (!a.CompareTo("HLT_DoubleMu3")) fHLTMu3 = true; 
	  if (0) cout << "HLT fired " << a << endl;
	  fGoodHLT = true; 
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

      // FIXME: replace (pt->fTrackQuality < 0) with something meaningful!
      if (TRACKQUALITY > 0 && (pt->fTrackQuality < 0)) {
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
      if (0 == (pt->fMuID & MUID)) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUID: " << pt->fMuID << endl;          
	fvGoodMuonsID[iC] = false; 
      }
      if (pt->fPlab.Perp() < MUPTLO) {
	if (fVerbose > 4) cout << "muon " << ps->fIndex << " failed MUPTLO: " << pt->fPlab.Perp() << endl;          
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

  if (fpCand->fPvIdx > -1) {
    TAnaVertex *pv = fpEvt->getPV(fpCand->fPvIdx); 
    fPvX = pv->fPoint.X(); 
    fPvY = pv->fPoint.Y(); 
    fPvZ = pv->fPoint.Z(); 
  } else {
    fPvX = -99.;
    fPvY = -99.;
    fPvZ = -99.;
  }

  fCandTM    = tmCand(fpCand); 
  fCandPt    = fpCand->fPlab.Perp();
  fCandEta   = fpCand->fPlab.Eta();
  fCandPhi   = fpCand->fPlab.Phi();
  fCandM     = fpCand->fMass;

  TAnaTrack *p0; 
  TAnaTrack *p1 = fpEvt->getSigTrack(fpCand->fSig1); 
  TAnaTrack *p2 = fpEvt->getSigTrack(fpCand->fSig1+1); 
  if (p1->fPlab.Perp() < p2->fPlab.Perp()) {
    p0 = p1; 
    p1 = p2; 
    p2 = p0; 
  }

  fMu1Pt     = p1->fPlab.Perp(); 
  fMu1Eta    = p1->fPlab.Eta(); 
  fMu1Phi    = p1->fPlab.Phi(); 
  fMu2Pt     = p2->fPlab.Perp(); 
  fMu2Eta    = p2->fPlab.Eta(); 
  fMu2Phi    = p2->fPlab.Phi(); 
  
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0); 
  TVector3 svpv(fpCand->fVtx.fPoint - fpEvt->getPV(pvidx)->fPoint); 
  double alpha = svpv.Angle(fpCand->fPlab);
  fCandCosA  = TMath::Cos(alpha);

  double iso = isoClassic(fpCand); 
  fCandIso   = iso; 
  fCandChi2  = fpCand->fVtx.fChi2;
  fCandDof   = fpCand->fVtx.fNdof;
  fCandProb  = fpCand->fVtx.fProb;
  fCandFLS3d = fpCand->fVtx.fD3d/fpCand->fVtx.fD3dE; 
  fCandFLSxy = fpCand->fVtx.fDxy/fpCand->fVtx.fDxyE; 

  if (fpCand->fNstTracks.size() == 0) {
    cout << "HHHHEEEELLLLPPPP" << endl;
    fCandDocaTrk = 99.;
  } else {
    fCandDocaTrk = fpCand->fNstTracks[0].second.first;
  }
  // ??
  TAnaTrack *t1 = fpEvt->getRecTrack(p1->fIndex); 
  TAnaTrack *t2 = fpEvt->getRecTrack(p2->fIndex); 
  double bmu1   = TMath::Sin(fpCand->fPlab.Angle(t1->fPlab));
  double bmu2   = TMath::Sin(fpCand->fPlab.Angle(t2->fPlab));
  fMu1IP        = fpCand->fVtx.fD3d*bmu1/t1->fTip;
  fMu2IP        = fpCand->fVtx.fD3d*bmu2/t2->fTip;

  // -- fill cut variables
  fWideMass = ((fpCand->fMass > 4.8) && (fpCand->fMass < 6.0)); 
  
  fGoodTracks = fvGoodTracks[fCandIdx];
  fGoodTracksPt = fvGoodTracksPt[fCandIdx];
  fGoodTracksEta = fvGoodTracksEta[fCandIdx];
  fGoodMuonsID  = fvGoodMuonsID[fCandIdx];
  fGoodMuonsPt  = fvGoodMuonsPt[fCandIdx]; 
  fGoodMuonsEta  = fvGoodMuonsEta[fCandIdx]; 

  fGoodPt = (fCandPt > CANDPTLO);
  fGoodEta = ((fCandEta > CANDETALO) && (fCandEta < CANDETAHI)); 
  fGoodCosA = (fCandCosA > CANDCOSALPHA); 
  fGoodIso = (fCandIso > CANDISOLATION); 
  fGoodChi2 = (fCandChi2 < CANDVTXCHI2);
  fGoodFLS =  ((fCandFLS3d > CANDFLS3D) && (fCandFLSxy > CANDFLSXY)); 

  fGoodDocaTrk = (fCandDocaTrk > CANDDOCATRK);
  fGoodIP = true; 

  fAnaCuts.update(); 
  
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
void bmmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("b1"))->Fill(fpEvt->nRecTracks()); 
}


// ----------------------------------------------------------------------
void bmmReader::fillCandidateHistograms() {

  // -- only candidate histograms below
  if (0 == fpCand) return;
  
  for (int i = 0; i < fAnaCuts.ncuts(); ++i) {
    if (fAnaCuts.singleCutTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dSi", i)))->Fill(fCandM);
    if (fAnaCuts.cumulativeCutTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dCu", i)))->Fill(fCandM);
    if (fAnaCuts.allOtherCutsTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dAo", i)))->Fill(fCandM);
    if (fAnaCuts.nMinus1CutsTrue(i)) ((TH1D*)fpHistFile->Get(Form("c%dAo", i)))->Fill(fCandM);
  }

  // -- Fill distributions
  fpTracksPt->fill(fMu1Pt, fCandM);
  fpTracksPt->fill(fMu2Pt, fCandM);

  fpMuonsID->fill(1, fCandM);
  fpMuonsID->fill(1, fCandM);
  fpMuonsPt->fill(fMu1Pt, fCandM);
  fpMuonsPt->fill(fMu2Pt, fCandM);
  fpMuonsEta->fill(fMu1Eta, fCandM);
  fpMuonsEta->fill(fMu2Eta, fCandM);

  fpAllEvents->fill(1, fCandM); 
  fpHLT->fill(1, fCandM); 
  fpPvZ->fill(fPvZ, fCandM); 
  fpPt->fill(fCandPt, fCandM); 
  fpEta->fill(fCandEta, fCandM); 
  fpCosA->fill(fCandCosA, fCandM);
  fpCosA0->fill(fCandCosA, fCandM);
  fpIso->fill(fCandIso, fCandM);
  fpChi2->fill(fCandChi2, fCandM);
  fpChi2Dof->fill(fCandChi2/fCandDof, fCandM); 
  fpProb->fill(fCandProb, fCandM);   
  fpFLS3d->fill(fCandFLS3d, fCandM); 
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

 
  // -- mass histograms for efficiencies
  for (int i = 0; i < fAnaCuts.ncuts(); ++i) {
    h = new TH1D(Form("c%dSi", i), fAnaCuts.getName(i), 30, 4.8, 6.0); 
    h = new TH1D(Form("c%dAo", i), fAnaCuts.getName(i), 30, 4.8, 6.0); 
    h = new TH1D(Form("c%dCu", i), fAnaCuts.getName(i), 30, 4.8, 6.0); 
    h = new TH1D(Form("c%dNm", i), fAnaCuts.getName(i), 30, 4.8, 6.0); 
  }

  h = new TH1D("analysisDistributions", "analysisDistributions", 100, 0., 100.); 
  fpAllEvents= bookDistribution("allevents", "allevents", "fWideMass", 10, 0., 10.);           
  fpHLT      = bookDistribution("hlt", "hlt", "fGoodHLT", 10, 0., 10.);           
  fpPvZ      = bookDistribution("pvz", "z_{PV} [cm]", "fGoodHLT", 50, -10., 10.);           
  fpTracksPt = bookDistribution("trackspt", "p_{T} [GeV]", "fGoodTracksPt", 50, 0., 25.);           
  fpMuonsID  = bookDistribution("muonsid", "muon id", "fGoodMuonsID", 10, 0., 10.); 
  fpMuonsPt  = bookDistribution("muonspt", "p_{T, #mu} [GeV]", "fGoodMuonsPt", 50, 0., 25.); 
  fpMuonsEta = bookDistribution("muonseta", "#eta_{#mu}", "fGoodMuonsEta", 50, -2.5, 2.5); 
  fpPt       = bookDistribution("pt", "p_{T, B} [GeV]", "fGoodPt", 50, 0., 25.); 
  fpEta      = bookDistribution("eta", "#eta_{B}", "fGoodEta", 50, -2.5, 2.5); 
  fpCosA     = bookDistribution("cosa", "cos(#alpha)", "fGoodCosA", 60, 0.97, 1.); 
  fpCosA0    = bookDistribution("cosa0", "cos(#alpha)", "fGoodCosA", 101, -1.01, 1.01); 
  fpIso      = bookDistribution("iso",  "I", "fGoodIso", 101, 0., 1.01); 
  fpChi2     = bookDistribution("chi2",  "#chi^{2}", "fGoodChi2", 100, 0., 50.);              
  fpChi2Dof  = bookDistribution("chi2dof",  "#chi^{2}/dof", "fGoodChi2", 100, 0., 50.);       
  fpProb     = bookDistribution("pchi2dof",  "P(#chi^{2},dof)", "fGoodChi2", 100, 0., 1.);    
  fpFLS3d    = bookDistribution("fls3d", "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 100, 0., 50.);  
  fpFLSxy    = bookDistribution("flsxy", "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 100, 0., 50.);  
  fpDocaTrk  = bookDistribution("docatrk", "d_{ca}(track))", "fGoodDocaTrk", 100, 0., 0.1);   
  fpIP1      = bookDistribution("ip1", "IP_{1}/lsin(#beta)", "fGoodIP", 100, -1., 3.);        
  fpIP2      = bookDistribution("ip2", "IP_{1}/lsin(#beta)", "fGoodIP", 100, -1., 3.);        

  
  h = new TH1D("b1", "Ntrk", 200, 0., 200.);
  h = new TH1D("bnc0", "NCand before selection", 20, 0., 20.);
  h = new TH1D("bnc1", "NCand after selection", 20, 0., 20.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun,               "run/I");
  fTree->Branch("evt",    &fEvt,               "evt/I");
  fTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fTree->Branch("l1t",    &fGoodL1T,           "l1t/O");
  fTree->Branch("mck",    &fGoodMCKinematics,  "mck/O");
  fTree->Branch("tm",     &fCandTM,            "tm/I");
  fTree->Branch("pt",     &fCandPt,            "pt/D");
  fTree->Branch("eta",    &fCandEta,           "eta/D");
  fTree->Branch("phi",    &fCandPhi,           "phi/D");
  fTree->Branch("m",      &fCandM,             "m/D");
  fTree->Branch("cosa",   &fCandCosA,          "cosa/D");
  fTree->Branch("iso",    &fCandIso,           "iso/D");
  fTree->Branch("chi2",   &fCandChi2,          "chi2/D");
  fTree->Branch("fls3d",  &fCandFLS3d,         "fls3d/D");
  fTree->Branch("flsxy",  &fCandFLSxy,         "flsxy/D");
  fTree->Branch("docatrk",&fCandDocaTrk,       "docatrk/D");
  fTree->Branch("m1pt",   &fMu1Pt,             "m1pt/D");
  fTree->Branch("m1eta",  &fMu1Eta,            "m1eta/D");
  fTree->Branch("m1phi",  &fMu1Phi,            "m1phi/D");
  fTree->Branch("m2pt",   &fMu2Pt,             "m2pt/D");
  fTree->Branch("m2eta",  &fMu2Eta,            "m2eta/D");
  fTree->Branch("m2phi",  &fMu2Phi,            "m2phi/D");
  fTree->Branch("m1ip",   &fMu1IP,             "m1ip/D");
  fTree->Branch("m2ip",   &fMu2IP,            "m2ip/D");
}


// ---------------------------------------------------------------------- 
AnalysisDistribution* bmmReader::bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn, ht, nbins, lo, hi); 
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX); 
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX); 
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX); 
  p->setAnalysisCuts(&fAnaCuts, hc); 

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
void bmmReader::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> bmmReader: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin; 
  string cstring = "B cand"; 
  while (is.getline(buffer, 200, '\n')) {
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
      ibin = 3; 
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
      ibin = 3; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
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
      CANDISOLATION = CutValue; ok = 1;
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
    if (!strcmp(CutName, "MUID")) {
      MUID = int(CutValue); ok = 1;
      if (dump) cout << "MUID:           " << MUID << endl;
      ibin = 200;
      hcuts->SetBinContent(ibin, MUID);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuID :: %d", CutName, MUID));
    }

    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; ok = 1;
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 201;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(#mu) :: %3.1f", CutName, MUPTLO));
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; ok = 1;
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 202;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(#mu) :: %3.1f", CutName, MUPTHI));
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; ok = 1;
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 203;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(#mu) :: %3.1f", CutName, MUETALO));
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; ok = 1;
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 204;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(#mu) :: %3.1f", CutName, MUETAHI));
    }

    if (!strcmp(CutName, "MUIP")) {
      MUIP = CutValue; ok = 1;
      if (dump) cout << "MUIP:           " << MUIP << endl;
      ibin = 205;
      hcuts->SetBinContent(ibin, MUIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: IP(#mu) :: %3.1f", CutName, MUIP));
    }


    if (!ok) cout << "==> bmmReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

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
