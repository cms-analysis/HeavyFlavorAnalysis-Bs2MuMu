#include "candAnaBs2JpsiPhi.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

#include "../macros/AnalysisDistribution.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::candAnaBs2JpsiPhi(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fGenK1Tmi = fGenK2Tmi = fRecK1Tmi = fRecK2Tmi = -1; 
  BLIND = 0; 
  cout << "==> candAnaBs2JpsiPhi: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::~candAnaBs2JpsiPhi() {
  cout << "==> candAnaBs2JpsiPhi: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::candAnalysis() {

  if (0 == fpCand) return;

  // -- Check for J/psi mass
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  double chi2(0.), ndof(0.); 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    //    cout << "i = " << i << " pD = " << pD << endl;
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp();
      fJpsiEta  = pD->fPlab.Eta();
      fJpsiPhi  = pD->fPlab.Phi(); 

      chi2 = pD->fVtx.fChi2;
      ndof = pD->fVtx.fNdof;
    }
    if (pD->fType == PHITYPE) {
      if ((MKKLO < pD->fMass) && (pD->fMass < MKKHI)) fGoodMKK = true;
      fMKK     = pD->fMass;
      fPhiPt   = pD->fPlab.Perp();
      fPhiEta  = pD->fPlab.Eta();
      fPhiPhi  = pD->fPlab.Phi(); 
    }
  }


  // -- Get Kaons
  TAnaTrack *p0; 
  TAnaTrack *p1(0);
  TAnaTrack *p2(0); 
  
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (TMath::Abs(p0->fMCID) != 321) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }

  if (0 == p1) {
    cout << "candAnaBs2JpsiPhi::candAnalysis  no kaon 1 found " << endl;
    return; 
  }
  if (0 == p2) {
    cout << "candAnaBs2JpsiPhi::candAnalysis  no kaon 2 found " << endl;
    return; 
  }
 
  fKa1Pt        = p1->fRefPlab.Perp(); 
  fKa1Eta       = p1->fRefPlab.Eta(); 
  fKa1Phi       = p1->fRefPlab.Phi(); 
  fKa1TkQuality = highPurity(p1);
  fKa1PtNrf     = p1->fPlab.Perp();
  fKa1EtaNrf    = p1->fPlab.Eta();

  fKa2Pt        = p2->fRefPlab.Perp(); 
  fKa2Eta       = p2->fRefPlab.Eta(); 
  fKa2Phi       = p2->fRefPlab.Phi(); 
  fKa2TkQuality = highPurity(p2);
  fKa2PtNrf     = p2->fPlab.Perp();
  fKa2EtaNrf    = p2->fPlab.Eta();

  fKa1Missid = tightMuon(p1);  // true for tight  muons 
  fKa2Missid = tightMuon(p2);
  fKa1MuMatch = doTriggerMatching(p1); // see if it matches HLT muon 
  fKa2MuMatch = doTriggerMatching(p2); // see if it matches HLT muon 

  if(0) { // special tests d.k.
    double mva=0;

    fKa1Missid2 = mvaMuon(p1,mva);  // true for tight  muons 
    fKa1MuMatch = doTriggerMatching(p1,false); // see if it matches HLT muon 
    fKa1MuMatch2 = doTriggerMatching(p1,true); // see if it matches HLT muon 
    fKa1MuMatchR = doTriggerMatchingR(p1,false); // see if it matches HLT muon 
    fKa1MuMatchR2 = doTriggerMatchingR(p1,true); // see if it matches HLT muon 
    fKa1MuMatchR3 = matchToMuon(p1,true); // see if it matches HLT muon 
    
    // 2nd Kaon
    fKa2Missid2 = mvaMuon(p2,mva);  // true for tight  muons 
    fKa2MuMatch = doTriggerMatching(p2,false); // see if it matches HLT muon 
    fKa2MuMatch2 = doTriggerMatching(p2,true); // see if it matches HLT muon 
    fKa2MuMatchR = doTriggerMatchingR(p2,false); // see if it matches HLT muon 
    fKa2MuMatchR2 = doTriggerMatchingR(p2,true); // see if it matches HLT muon 
    fKa2MuMatchR3 = matchToMuon(p2,true); // see if it matches HLT muon 

    //if(fKa1Missid || fKa2Missid) cout<<"missid "<<fKa1Missid<<" "<<fKa2Missid<<" "<<fKa1MuMatch<<" "<<fKa2MuMatch<<endl;
  }

  if (fCandTmi > -1) {
    TGenCand *pg1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p1->fIndex)->getGenIndex());
    fKa1PtGen     = pg1->fP.Perp();
    fKa1EtaGen    = pg1->fP.Eta();
    TGenCand *pg2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p2->fIndex)->getGenIndex());
    fKa2PtGen     = pg2->fP.Perp();
    fKa2EtaGen    = pg2->fP.Eta();
  } else {
    fKa1PtGen     = -99.;
    fKa1EtaGen    = -99.;
    fKa2PtGen     = -99.;
    fKa2EtaGen    = -99.;
  }

  fDeltaR  = p1->fPlab.DeltaR(p2->fPlab); 

  TLorentzVector ka1, ka2; 
  ka1.SetPtEtaPhiM(fKa1Pt, fKa1Eta, fKa1Phi, MKAON); 
  ka2.SetPtEtaPhiM(fKa2Pt, fKa2Eta, fKa2Phi, MKAON); 

  TLorentzVector phiCand = ka1 + ka2; 
  fMKK     = phiCand.M();
  fPhiPt   = phiCand.Pt();
  fPhiEta  = phiCand.Eta();
  fPhiPhi  = phiCand.Phi();

  fGoodDeltaR = (fDeltaR < DELTAR);
  fGoodMKK    = ((MKKLO < fMKK ) && (fMKK < MKKHI)); 
  
  candAna::candAnalysis();

  fPreselection = fPreselection && fGoodJpsiMass && fGoodMKK && fGoodDeltaR; 
  fPreselection = fPreselection && fWideMass;

  // -- overwrite specific variables
  fCandTau     = fCandFL3d*MBPLUS/fCandP/TMath::Ccgs();
  fCandChi2    = chi2; 
  fCandDof     = ndof;
  fCandChi2Dof = chi2/ndof;
  
  
  if(1) { // misid test d.k.
    if( (p1->fIndex == fpMuon1->fIndex) || (p1->fIndex ==fpMuon2->fIndex) ) 
      cout<<" Kaon1 is a MUON "<<fEvt<<" "<<fpCand<<" "<<p1->fIndex<<" "<<fpMuon1->fIndex<<" "<<fpMuon2->fIndex<<" "<<fEvt<<endl;
    if( (p2->fIndex == fpMuon1->fIndex) || (p2->fIndex ==fpMuon2->fIndex) ) 
      cout<<" Kaon2 is a MUON "<<fEvt<<" "<<fpCand<<" "<<p2->fIndex<<" "<<fpMuon1->fIndex<<" "<<fpMuon2->fIndex<<" "<<fEvt<<endl;

    TVector3 muonMom1 = fpMuon1->fPlab;
    TVector3 muonMom2 = fpMuon2->fPlab;

    TVector3 trackMom = p1->fPlab;  // test track momentum
    double dR1 = muonMom1.DeltaR(trackMom);
    double dR2 = muonMom2.DeltaR(trackMom);

    if(dR1<dR2) { fKa1MuMatchR4 = dR1; fKa1MuMatchR5 = dR2;}
    else        { fKa1MuMatchR4 = dR2; fKa1MuMatchR5 = dR1;}

    trackMom = p2->fPlab;  // test track momentum
    dR1 = muonMom1.DeltaR(trackMom);
    dR2 = muonMom2.DeltaR(trackMom);

    if(dR1<dR2) { fKa2MuMatchR4 = dR1; fKa2MuMatchR5 = dR2;}
    else        { fKa2MuMatchR4 = dR2; fKa2MuMatchR5 = dR1;}
  } // end special test 

  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(4); 
}

// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::moreBasicCuts() {
  cout << "   candAnaBs2JpsiPhi: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
  fAnaCuts.addCut("fGoodDeltaR", "Delta R(KK)", fGoodDeltaR); 
  fAnaCuts.addCut("fGoodMKK", "m(KK) [GeV]", fGoodMKK); 
}

// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::genMatch() {
  
  fGenM1Tmi = fGenM2Tmi = fGenK1Tmi = -1; 
  fNGenPhotons = 0; 

  TGenCand *pC(0), *pB(0), *pPsi(0), *pPhi(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  int nb(0), ngamma(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1; 
      if (nb > 2) continue; // skip B decays where more than J/psi and phi came from B
      ngamma = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id); 
	if (22 == TMath::Abs(pC->fID)) ++ngamma;
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC; 
	  pM1 = pM2 = 0;
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenTWithIndex(idd); 
	      } else {
		pM2 = fpEvt->getGenTWithIndex(idd); 
	      }
	    }
	  }
	} else if (333 == TMath::Abs(pC->fID)) {
	  pPhi = fpEvt->getGenTWithIndex(id); 
	  pK1 = pK2 = 0;
	  for (int idd = pPhi->fDau1; idd <= pPhi->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (321 == TMath::Abs(pC->fID)) {
	      if (0 == pK1) {
		pK1 = fpEvt->getGenTWithIndex(idd); 
	      } else {
		pK2 = fpEvt->getGenTWithIndex(idd); 
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK1 && 0 != pK2 
	  && (pPsi->fMom1 == pPhi->fMom1)
	  ) {
	goodMatch = true; 
	fNGenPhotons = ngamma;
	break;
      }
    }
  }

  if (!goodMatch) {
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }

  fGenBTmi = -1; 
  fKa1GenID = -99999; 
  fKa2GenID = -99999; 
  if (goodMatch) {
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fKa1GenID = pK1->fID;
    fKa2GenID = pK2->fID;
    fGenBTmi = pB->fNumber; 
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1); 
    // the meson is the original except if it oscillated
    if (531 != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag(); 
    fGenLifeTime = x*m/p/TMath::Ccgs();
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber; 
      fGenM2Tmi = pM2->fNumber; 
    } else {
      fGenM1Tmi = pM2->fNumber; 
      fGenM2Tmi = pM1->fNumber; 
    }
    if (pK1->fP.Perp() > pK2->fP.Perp()) {
      fGenK1Tmi = pK1->fNumber; 
      fGenK2Tmi = pK2->fNumber; 
    } else {
      fGenK1Tmi = pK2->fNumber; 
      fGenK2Tmi = pK1->fNumber; 
    }
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
    fGenK1Tmi = -1;  
    fGenK2Tmi = -1;  
  }


  // -- check that only one reco track is matched to each gen particle
  //    else skip the *event*!
  TSimpleTrack *pT(0);
  int cntM1(0), cntM2(0), cntK1(0), cntK2(0); 
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (fGenM1Tmi > -1 && fGenM1Tmi == pT->getGenIndex()) ++cntM1;
    if (fGenM2Tmi > -1 && fGenM2Tmi == pT->getGenIndex()) ++cntM2;
    if (fGenK1Tmi > -1 && fGenK1Tmi == pT->getGenIndex()) ++cntK1;
    if (fGenK2Tmi > -1 && fGenK2Tmi == pT->getGenIndex()) ++cntK2;
  }

  static int cntBadEvents = 0; 
  if (cntM1 > 1 || cntM2 > 1 || cntK1 > 1 || cntK2 > 1) {
    cout << "BAD BAD event: multiple reco tracks matched to the same gen particle! " << ++cntBadEvents 
	 << ": " << cntM1 << " .. " << cntM2 << " .. " << cntK1 << " .. " << cntK2
	 << " (gen-reco matches) " << endl;
    fBadEvent = true; 
  }

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = fRecK2Tmi =-1; 
  TSimpleTrack *pT(0);
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (pT->getGenIndex() < 0) continue;

    // -- muon 1
    if (fGenM1Tmi > -1 && pT->getGenIndex() == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (fGenM2Tmi > -1 && pT->getGenIndex() == fGenM2Tmi) {
      fRecM2Tmi = i; 
    }

    // -- kaon 1
    if (fGenK1Tmi > -1 && pT->getGenIndex() == fGenK1Tmi) {
      fRecK1Tmi = i; 
    }

    // -- kaon 2
    if (fGenK2Tmi > -1 && pT->getGenIndex() == fGenK2Tmi) {
      fRecK2Tmi = i; 
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecK1Tmi > -1 && fRecK2Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fRecK1Tmi = " << fRecK1Tmi << " matched to fGenK1Tmi = " << fGenK1Tmi << endl;
    cout << "fRecK2Tmi = " << fRecK2Tmi << " matched to fGenK2Tmi = " << fGenK2Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::candMatch() {

  fCandTmi = -1;   
  int idx(-1), type(-1); 
  int d1Matched(0), d2Matched(0), d3Matched(0), d4Matched(0); 
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 
    if (TYPE != pCand->fType) continue;
    
    d1Matched = d2Matched = d3Matched = d4Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
      //       if (fGenM1Tmi > -1) cout << "  --> " << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
      if (fVerbose > 10) {
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
      }
      if (fRecM1Tmi > -1 && type == 13 && idx == fRecM1Tmi) {
	d1Matched = 1; 
      }
      if (fRecM2Tmi > -1 && type == 13 && idx == fRecM2Tmi) {
	d2Matched = 1; 
      }
      if (fRecK1Tmi > -1 && type == 321 && idx == fRecK1Tmi) {
	d3Matched = 1; 
      }
      if (fRecK2Tmi > -1 && type == 321 && idx == fRecK2Tmi) {
	d4Matched = 1; 
      }
    }
    
    if (d1Matched && d2Matched && d3Matched && d4Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi  << " " << fRecK2Tmi 
	 << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::bookHist() {
  cout << "==>candAnaBs2JpsiPhi: bookHist" << endl;
  candAna::bookHist();

  moreReducedTree(fTree);
  moreReducedTree(fAmsTree);

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");

  fEffTree->Branch("k2pt",   &fETk2pt,            "k2pt/F");
  fEffTree->Branch("g4pt",   &fETg4pt,            "g4pt/F");
  fEffTree->Branch("k2eta",  &fETk2eta,           "k2eta/F");
  fEffTree->Branch("g4eta",  &fETg4eta,           "g4eta/F");
  fEffTree->Branch("k2q",    &fETk2q,             "k2q/I");
  fEffTree->Branch("k2gt",   &fETk2gt,            "k2gt/O");

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::moreReducedTree(TTree *t) {

  // -- Additional reduced tree variables
  t->Branch("mpsi",  &fJpsiMass, "mpsi/D");
  t->Branch("psipt", &fJpsiPt,   "psipt/D");
  t->Branch("psieta",&fJpsiEta,  "psieta/D");
  t->Branch("psiphi",&fJpsiPhi,  "psiphi/D");
  t->Branch("mkk",   &fMKK,      "mkk/D");
  t->Branch("phipt", &fPhiPt,    "phipt/D");
  t->Branch("phieta",&fPhiEta,   "phieta/D");
  t->Branch("phiphi",&fPhiPhi,   "phiphi/D");
  t->Branch("dr",    &fDeltaR,   "dr/D");

  t->Branch("k1pt",  &fKa1Pt,    "k1pt/D");
  t->Branch("k1eta", &fKa1Eta,   "k1eta/D");
  t->Branch("k1phi", &fKa1Phi,   "k1phi/D");
  t->Branch("k1gt",  &fKa1TkQuality,"k1gt/I");
  t->Branch("k2pt",  &fKa2Pt,    "k2pt/D");
  t->Branch("k2eta", &fKa2Eta,   "k2eta/D");
  t->Branch("k2phi", &fKa2Phi,   "k2phi/D");
  t->Branch("k2gt",  &fKa2TkQuality,"k2gt/I");



  t->Branch("t3pt",  &fKa1PtNrf, "t3pt/D");
  t->Branch("t3eta", &fKa1EtaNrf,"t3eta/D");

  t->Branch("t4pt",  &fKa2PtNrf, "t4pt/D");
  t->Branch("t4eta", &fKa2EtaNrf,"t4eta/D");

  t->Branch("g3pt", &fKa1PtGen,  "g3pt/D");
  t->Branch("g3eta",&fKa1EtaGen, "g3eta/D");
  t->Branch("g3id", &fKa1GenID,  "g3id/I");

  t->Branch("g4pt", &fKa2PtGen,  "g4pt/D");
  t->Branch("g4eta",&fKa2EtaGen, "g4eta/D");
  t->Branch("g4id", &fKa2GenID,  "g4id/I");

  t->Branch("k1missid",  &fKa1Missid,    "k1missid/O");
  t->Branch("k2missid",  &fKa2Missid,    "k2missid/O");
  t->Branch("k1mumatch", &fKa1MuMatch,  "k1mumatch/O");
  t->Branch("k2mumatch", &fKa2MuMatch,  "k2mumatch/O");

  if(0) { // for testing d.k.

    t->Branch("k1missid2",  &fKa1Missid2,    "k1missid2/O");
    t->Branch("k1mumatch2", &fKa1MuMatch2,    "k1mumatch2/O");
    
    t->Branch("k1mumatchr", &fKa1MuMatchR,    "k1mumatchr/F");
    t->Branch("k1mumatchr2", &fKa1MuMatchR2,    "k1mumatchr2/F");
    t->Branch("k1mumatchr3", &fKa1MuMatchR3,    "k1mumatchr3/F");
    t->Branch("k1mumatchr4", &fKa1MuMatchR4,    "k1mumatchr4/F");
    t->Branch("k1mumatchr5", &fKa1MuMatchR5,    "k1mumatchr5/F");
    
    t->Branch("k2missid2",  &fKa2Missid2,    "k2missid2/O");
    t->Branch("k2mumatch2", &fKa2MuMatch2,   "k2mumatch2/O");
    
    t->Branch("k2mumatchr",  &fKa2MuMatchR,   "k2mumatchr/F");
    t->Branch("k2mumatchr2", &fKa2MuMatchR2,  "k2mumatchr2/F");
    t->Branch("k2mumatchr3", &fKa2MuMatchR3,  "k2mumatchr3/F");
    t->Branch("k2mumatchr4", &fKa2MuMatchR4,  "k2mumatchr4/F");
    t->Branch("k2mumatchr5", &fKa2MuMatchR5,  "k2mumatchr5/F");
  } // end testing 

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::fillCandidateHistograms(int offset) {
  candAna::fillCandidateHistograms(offset); 
}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::efficiencyCalculation() {

  fGoodEffCand = false;

  // -- gen level 
  TGenCand *pB(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenK1Tmi || -1 == fGenK2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenTWithIndex(fGenBTmi); 
  pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi); 
  pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi); 
  pK1 = fpEvt->getGenTWithIndex(fGenK1Tmi); 
  pK2 = fpEvt->getGenTWithIndex(fGenK2Tmi); 

  // -- reco level
  TSimpleTrack *prM1(0), *prM2(0), *prK1(0), *prK2(0); 
  double bla(0); 
  int m1Matched(0), m2Matched(0), k1Matched(0), k2Matched(0), m1ID(0), m1tmID(0), m1mvaID(0), m2ID(0), m2tmID(0), m2mvaID(0), 
    m1GT(0), m2GT(0), k1GT(0), k2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1; 
    prM1 = fpEvt->getSimpleTrack(fRecM1Tmi); 
    if (tightMuon(prM1)) m1tmID = 1; 
    if (mvaMuon(prM1, bla)) m1mvaID = 1; 
    if (prM1->getHighPurity()) {
      m1GT = 1; 
    } else {
      m1GT = 0;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1; 
    prM2 = fpEvt->getSimpleTrack(fRecM2Tmi); 
    if (tightMuon(prM2)) m2tmID = 1; 
    if (mvaMuon(prM2, bla)) m2mvaID = 1; 
    if (prM2->getHighPurity()) {
      m2GT = 1; 
    } else {
      m2GT = 0;
    }
  } 

  if (fRecK1Tmi > -1) {
    k1Matched = 1; 
    prK1 = fpEvt->getSimpleTrack(fRecK1Tmi); 
    if (prK1->getHighPurity()) {
      k1GT = 1; 
    } else {
      k1GT = 0;
    }
  } 

  if (fRecK2Tmi > -1) {
    k2Matched = 1; 
    prK2 = fpEvt->getSimpleTrack(fRecK2Tmi); 
    if (prK2->getHighPurity()) {
      k2GT = 1; 
    } else {
      k2GT = 0;
    }
  } 

  // -- cand level 
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fpEvt->getCand(fCandTmi);
  }

  m1ID = m1tmID; 
  m2ID = m2tmID; 

  // -- EffTree filling for all events with a signal decay
  fETgpt   = pB->fP.Perp(); 
  fETgeta  = pB->fP.Eta(); 
  fETg1pt  = pM1->fP.Perp(); 
  fETg1eta = pM1->fP.Eta(); 
  fETg2pt  = pM2->fP.Perp(); 
  fETg2eta = pM2->fP.Eta(); 
  fETg3pt  = pK1->fP.Perp(); 
  fETg3eta = pK1->fP.Eta(); 
  fETg4pt  = pK2->fP.Perp(); 
  fETg4eta = pK2->fP.Eta(); 
  if (m1Matched) {
    fETm1pt  = prM1->getP().Perp(); 
    fETm1eta = prM1->getP().Eta(); 
    fETm1q   = prM1->getCharge();
    fETm1gt  = (m1GT>0?true:false); 
    fETm1id  = (m1ID>0?true:false);
    fETm1tmid  = (m1tmID>0?true:false);
    fETm1mvaid = (m1mvaID>0?true:false);
  } else {
    fETm1pt  = -99.; 
    fETm1eta = -99.; 
    fETm1q   = -99;
    fETm1gt  = false; 
    fETm1id  = false;
    fETm1tmid  = false;
    fETm1mvaid = false;
  }
  if (m2Matched) {
    fETm2pt  = prM2->getP().Perp(); 
    fETm2eta = prM2->getP().Eta(); 
    fETm2q   = prM2->getCharge();
    fETm2gt  = (m2GT>0?true:false); 
    fETm2id  = (m2ID>0?true:false);
    fETm2tmid  = (m2tmID>0?true:false);
    fETm2mvaid = (m2mvaID>0?true:false);
  } else {
    fETm2pt  = -99.; 
    fETm2eta = -99.; 
    fETm2q   = -99;
    fETm2gt  = false; 
    fETm2id  = false;
    fETm2tmid  = false;
    fETm2mvaid = false;
  }
  if (k1Matched) {
    fETk1pt  = prK1->getP().Perp(); 
    fETk1eta = prK1->getP().Eta(); 
    fETk1q   = prK1->getCharge();
    fETk1gt  = (k1GT>0?true:false); 
  } else {
    fETk1pt  = -99.; 
    fETk1eta = -99.; 
    fETk1q   = -99;
    fETk1gt  = false; 
  }
  if (k2Matched) {
    fETk2pt  = prK2->getP().Perp(); 
    fETk2eta = prK2->getP().Eta(); 
    fETk2q   = prK2->getCharge();
    fETk2gt  = (k2GT>0?true:false); 
  } else {
    fETk2pt  = -99.; 
    fETk2eta = -99.; 
    fETk2q   = -99;
    fETk2gt  = false; 
  }
  if (pCand) {
    fETcandMass = pCand->fMass; 
  } else {
    fETcandMass = -99.;
  }

  fEffTree->Fill(); 

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaBs2JpsiPhi: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "JPSITYPE")) {
      JPSITYPE = static_cast<int>(CutValue);
      if (dump) cout << "JPSITYPE:      " << JPSITYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, JPSITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: J/#psi ID :: %d", CutName, JPSITYPE));
    }

    if (!strcmp(CutName, "JPSIMASSLO")) {
      JPSIMASSLO = CutValue; 
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSLO << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, JPSIMASSLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSLO));
    }

    if (!strcmp(CutName, "JPSIMASSHI")) {
      JPSIMASSHI = CutValue; 
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSHI << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, JPSIMASSHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSHI));
    }

    if (!strcmp(CutName, "PHITYPE")) {
      PHITYPE = static_cast<int>(CutValue); 
      if (dump) cout << "PHITYPE:      " << PHITYPE << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, PHITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #phi ID :: %d", CutName, PHITYPE));
    }

    if (!strcmp(CutName, "MKKLO")) {
      MKKLO = CutValue; 
      if (dump) cout << "MKKLO:           " << MKKLO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(KK) :: %3.1f", CutName, MKKLO));
    }

    if (!strcmp(CutName, "MKKHI")) {
      MKKHI = CutValue; 
      if (dump) cout << "MKKHI:           " << MKKHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(KK) :: %3.1f", CutName, MKKHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }


    
  }

}
