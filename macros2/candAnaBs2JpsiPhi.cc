#include "candAnaBs2JpsiPhi.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

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
      //       cout << "type = " << pD->fType 
      // 	   << " with mass = " << pD->fMass 
      // 	   << " fGoodJpsiMass = " << fGoodJpsiMass 
      // 	   << endl;
    }
    if (pD->fType == PHITYPE) {
      if ((MKKLO < pD->fMass) && (pD->fMass < MKKHI)) fGoodMKK = true;
      fMKK     = pD->fMass;
      fPhiPt   = pD->fPlab.Perp();
      fPhiEta  = pD->fPlab.Eta();
      fPhiPhi  = pD->fPlab.Phi(); 
      //       cout << "type = " << pD->fType 
      // 	   << " with mass = " << pD->fMass 
      // 	   << " fGoodMKK = " << fGoodMKK 
      // 	   << endl;
    }
  }

  //   // -- special case for truth candidates (which have no daughter cands)
  //   if (fpCand->fType > 999999) {
  //     TAnaTrack *p0; 
  //     TAnaTrack *p1(0);
  //     TAnaTrack *p2(0); 
  
  //     for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
  //       p0 = fpEvt->getSigTrack(it);     
  //       if (TMath::Abs(p0->fMCID) != 13) continue;
  //       if (0 == p1) {
  // 	p1 = p0; 
  //       } else {
  // 	p2 = p0; 
  //       }
  //     }
  
  //     if (0 == p1) {
  //       cout << "candAnaBs2JpsiPhi::candAnalysis  no muon 1 found " << endl;
  //       return; 
  //     }
  //     if (0 == p2) {
  //       cout << "candAnaBs2JpsiPhi::candAnalysis  no muon 2 found " << endl;
  //       return; 
  //     }
  
  //     TLorentzVector mu1, mu2; 
  //     mu1.SetPtEtaPhiM(p1->fPlab.Perp(), p1->fPlab.Eta(), p1->fPlab.Phi(), MMUON); 
  //     mu2.SetPtEtaPhiM(p2->fPlab.Perp(), p2->fPlab.Eta(), p2->fPlab.Phi(), MMUON); 
  
  //     TLorentzVector psi = mu1 + mu2; 
  //     if ((JPSIMASSLO < psi.M()) && (psi.M() < JPSIMASSHI)) fGoodJpsiMass = true;
  //     fJpsiMass = psi.M();
  //     fJpsiPt   = psi.Pt();
  //     fJpsiEta  = psi.Eta();
  //     fJpsiPhi  = psi.Phi(); 
  //   }    

  // -- Get Kaons
  TAnaTrack *p0; 
  TAnaTrack *p1(0), *ps1(0);
  TAnaTrack *p2(0), *ps2(0); 
  
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
 
  // -- switch to RecTracks!
  ps1 = p1; 
  p1 = fpEvt->getRecTrack(p1->fIndex);
  ps2 = p2; 
  p2 = fpEvt->getRecTrack(p2->fIndex);


  fKa1Pt        = p1->fPlab.Perp(); 
  fKa1Eta       = p1->fPlab.Eta(); 
  fKa1Phi       = p1->fPlab.Phi(); 
  fKa1TkQuality = p1->fTrackQuality & TRACKQUALITY;
  fKa1PtNrf     = ps1->fPlab.Perp();
  fKa1EtaNrf    = ps1->fPlab.Eta();

  fKa2Pt        = p2->fPlab.Perp(); 
  fKa2Eta       = p2->fPlab.Eta(); 
  fKa2Phi       = p2->fPlab.Phi(); 
  fKa2TkQuality = p2->fTrackQuality & TRACKQUALITY;
  fKa2PtNrf     = ps2->fPlab.Perp();
  fKa2EtaNrf    = ps2->fPlab.Eta();

  if (fCandTmi > -1) {
    TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(p1->fIndex)->fGenIndex);
    fKa1PtGen     = pg1->fP.Perp();
    fKa1EtaGen    = pg1->fP.Eta();
    TGenCand *pg2 = fpEvt->getGenCand(fpEvt->getRecTrack(p2->fIndex)->fGenIndex);
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
  
  fPreselection = fPreselection && fGoodJpsiMass && fGoodMKK && fGoodDeltaR; 


  candAna::candAnalysis();
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(3); 
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
  int npsi(0), nphi(0), nb(0), ngamma(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1; 
      if (nb > 2) continue; // skip B decays where more than J/psi and phi came from B
      ngamma = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (22 == TMath::Abs(pC->fID)) ++ngamma;
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC; 
	  npsi = pPsi->fDau2 - pPsi->fDau1 + 1; 
	  pM1 = pM2 = 0;
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenCand(idd); 
	      } else {
		pM2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	} else if (333 == TMath::Abs(pC->fID)) {
	  pPhi = fpEvt->getGenCand(id); 
	  nphi = pPhi->fDau2 - pPhi->fDau1 + 1; 
	  pK1 = pK2 = 0;
	  for (int idd = pPhi->fDau1; idd <= pPhi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (321 == TMath::Abs(pC->fID)) {
	      if (0 == pK1) {
		pK1 = fpEvt->getGenCand(idd); 
	      } else {
		pK2 = fpEvt->getGenCand(idd); 
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
  if (goodMatch) {
    fGenBTmi = pB->fNumber; 
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

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = fRecK2Tmi =-1; 
  TAnaTrack *pT(0);
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i); 
    if (pT->fGenIndex < 0) continue;

    // -- muon 1
    if (fGenM1Tmi > -1 && pT->fGenIndex == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (fGenM2Tmi > -1 && pT->fGenIndex == fGenM2Tmi) {
      fRecM2Tmi = i; 
    }

    // -- kaon 1
    if (fGenK1Tmi > -1 && pT->fGenIndex == fGenK1Tmi) {
      fRecK1Tmi = i; 
    }

    // -- kaon 2
    if (fGenK2Tmi > -1 && pT->fGenIndex == fGenK2Tmi) {
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

  // -- Additional reduced tree variables
  fTree->Branch("mpsi",  &fJpsiMass, "mpsi/D");
  fTree->Branch("psipt", &fJpsiPt,   "psipt/D");
  fTree->Branch("psieta",&fJpsiEta,  "psieta/D");
  fTree->Branch("psiphi",&fJpsiPhi,  "psiphi/D");
  fTree->Branch("mkk",   &fMKK,      "mkk/D");
  fTree->Branch("phipt", &fPhiPt,    "phipt/D");
  fTree->Branch("phieta",&fPhiEta,   "phieta/D");
  fTree->Branch("phiphi",&fPhiPhi,   "phiphi/D");
  fTree->Branch("dr",    &fDeltaR,   "dr/D");

  fTree->Branch("k1pt",  &fKa1Pt,    "k1pt/D");
  fTree->Branch("k1eta", &fKa1Eta,   "k1eta/D");
  fTree->Branch("k1phi", &fKa1Phi,   "k1phi/D");
  fTree->Branch("k1gt",  &fKa1TkQuality,"k1gt/I");
  fTree->Branch("k2pt",  &fKa2Pt,    "k2pt/D");
  fTree->Branch("k2eta", &fKa2Eta,   "k2eta/D");
  fTree->Branch("k2phi", &fKa2Phi,   "k2phi/D");
  fTree->Branch("k2gt",  &fKa2TkQuality,"k2gt/I");

  fTree->Branch("t3pt",  &fKa1PtNrf, "t3pt/D");
  fTree->Branch("t3eta", &fKa1EtaNrf,"t3eta/D");

  fTree->Branch("t4pt",  &fKa2PtNrf, "t4pt/D");
  fTree->Branch("t4eta", &fKa2EtaNrf,"t4eta/D");

  fTree->Branch("g3pt", &fKa1PtGen,  "g3pt/D");
  fTree->Branch("g3eta",&fKa1EtaGen, "g3eta/D");
  fTree->Branch("g4pt", &fKa2PtGen,  "g4pt/D");
  fTree->Branch("g4eta",&fKa2EtaGen, "g4eta/D");

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


  string name; 
  int i(0); 
  for (map<string, int>::iterator imap = fRegion.begin(); imap != fRegion.end(); ++imap) {  
    i    = imap->second;
    name = imap->first + "_";

    //    cout << "booking with offset = " << i << Form("   %smpsi", name.c_str()) << endl;

    fpMpsi[i]    = bookDistribution(Form("%smpsi", name.c_str()), "m(J/#psi) [GeV]", "fGoodJpsiMass", 40, 2.8, 3.4);           
    fpDeltaR[i]  = bookDistribution(Form("%sdeltar", name.c_str()), "#Delta R", "fGoodDeltaR", 50, 0., 1.);           
    fpMKK[i]     = bookDistribution(Form("%smkk", name.c_str()), "m(KK) [GeV]", "fGoodMKK", 50, 0.95, 1.15);           
    
    fpKaonsPt[i] = bookDistribution(Form("%skaonspt", name.c_str()), "p_{T, K} [GeV]", "fGoodTracksPt", 25, 0., 25.);           
    fpKaonsEta[i]= bookDistribution(Form("%skaonseta", name.c_str()), "#eta_{K}", "fGoodTracksEta", 25, -2.5, 2.5);
    fpPsiPt[i]   = bookDistribution(Form("%spsipt", name.c_str()), "p_{T, J/#psi} [GeV]", "fGoodTracksPt", 25, 0., 50.);           
    fpPsiEta[i]  = bookDistribution(Form("%spsieta", name.c_str()), "#eta_{J/#psi}", "fGoodTracksEta", 25, -2.5, 2.5);  
    fpPhiPt[i]   = bookDistribution(Form("%sphipt", name.c_str()), "p_{T, #phi} [GeV]", "fGoodTracksPt", 25, 0., 25.);
    fpPhiEta[i]  = bookDistribution(Form("%sphieta", name.c_str()), "#eta_{#phi}", "fGoodTracksEta", 25, -2.5, 2.5);  
  }

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::fillCandidateHistograms(int offset) {

  fpDeltaR[offset]->fill(fDeltaR, fCandM); 
  fpMKK[offset]->fill(fMKK, fCandM); 

  fpMpsi[offset]->fill(fJpsiMass, fCandM); 
  fpTracksPt[offset]->fill(fKa1Pt, fCandM);
  fpTracksPt[offset]->fill(fKa2Pt, fCandM);
  fpTracksEta[offset]->fill(fKa1Eta, fCandM);
  fpTracksEta[offset]->fill(fKa2Eta, fCandM);

  fpKaonsPt[offset]->fill(fKa1Pt, fCandM);
  fpKaonsPt[offset]->fill(fKa2Pt, fCandM);
  fpKaonsEta[offset]->fill(fKa1Eta, fCandM);
  fpKaonsEta[offset]->fill(fKa2Eta, fCandM);

  fpPsiPt[offset]->fill(fJpsiPt, fCandM); 
  fpPsiEta[offset]->fill(fJpsiEta, fCandM); 
  fpPhiPt[offset]->fill(fPhiPt, fCandM); 
  fpPhiEta[offset]->fill(fPhiEta, fCandM); 

  candAna::fillCandidateHistograms(offset); 

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::efficiencyCalculation() {

  fGoodEffCand = false;

//   ((TH1D*)fpHistFile->Get("efficiency"))->Fill(0); 
//   ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(1, "all events"); 


  // -- gen level 
  TGenCand *pB(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenK1Tmi || -1 == fGenK2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenCand(fGenBTmi); 
  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 
  pK1 = fpEvt->getGenCand(fGenK1Tmi); 
  pK2 = fpEvt->getGenCand(fGenK2Tmi); 

//   ((TH1D*)fpHistFile->Get("efficiency"))->Fill(1); 
//   ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(2, "gen signal decays"); 


  // -- reco level
  TAnaTrack *prM1(0), *prM2(0), *prK1(0), *prK2(0); 
  int m1Matched(0), m2Matched(0), k1Matched(0), k2Matched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0), k1GT(0), k2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1; 
    prM1 = fpEvt->getRecTrack(fRecM1Tmi); 
    if (goodMuon(prM1)) m1ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM1->fTrackQuality & TRACKQUALITY))) {
      m1GT = 0; 
    } else {
      m1GT = 1;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1; 
    prM2 = fpEvt->getRecTrack(fRecM2Tmi); 
    if (goodMuon(prM2)) m2ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM2->fTrackQuality & TRACKQUALITY))) {
      m2GT = 0; 
    } else {
      m2GT = 1;
    }
  } 

  if (fRecK1Tmi > -1) {
    k1Matched = 1; 
    prK1 = fpEvt->getRecTrack(fRecK1Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK1->fTrackQuality & TRACKQUALITY))) {
      k1GT = 0; 
    } else {
      k1GT = 1;
    }
  } 

  if (fRecK2Tmi > -1) {
    k2Matched = 1; 
    prK2 = fpEvt->getRecTrack(fRecK2Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK2->fTrackQuality & TRACKQUALITY))) {
      k2GT = 0; 
    } else {
      k2GT = 1;
    }
  } 

  // -- cand level 
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fpEvt->getCand(fCandTmi);
  }


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
    fETm1pt  = prM1->fPlab.Perp(); 
    fETm1eta = prM1->fPlab.Eta(); 
    fETm1q   = prM1->fQ;
    fETm1gt  = (m1GT>0?true:false); 
    fETm1id  = (m1ID>0?true:false);
  } else {
    fETm1pt  = -99.; 
    fETm1eta = -99.; 
    fETm1q   = -99;
    fETm1gt  = false; 
    fETm1id  = false;
  }
  if (m2Matched) {
    fETm2pt  = prM2->fPlab.Perp(); 
    fETm2eta = prM2->fPlab.Eta(); 
    fETm2q   = prM2->fQ;
    fETm2gt  = (m2GT>0?true:false); 
    fETm2id  = (m2ID>0?true:false);
  } else {
    fETm2pt  = -99.; 
    fETm2eta = -99.; 
    fETm2q   = -99;
    fETm2gt  = false; 
    fETm2id  = false;
  }
  if (k1Matched) {
    fETk1pt  = prK1->fPlab.Perp(); 
    fETk1eta = prK1->fPlab.Eta(); 
    fETk1q   = prK1->fQ;
    fETk1gt  = (k1GT>0?true:false); 
  } else {
    fETk1pt  = -99.; 
    fETk1eta = -99.; 
    fETk1q   = -99;
    fETk1gt  = false; 
  }
  if (k2Matched) {
    fETk2pt  = prK2->fPlab.Perp(); 
    fETk2eta = prK2->fPlab.Eta(); 
    fETk2q   = prK2->fQ;
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
  //   cout << Form("%4d", fEvent) << " bmBs2JpsiPhiReader: m = " << fETcandMass << " from cand " << pCand 
  //        << " mu gen: " << fGenM1Tmi << " " << fGenM2Tmi 
  //        << " mu gen: " << fETg1pt << " " << fETg2pt 
  //        << endl;

  fEffTree->Fill(); 


//   // -- results...
//   if ((TMath::Abs(pM1->fP.Eta()) < 2.5)
//       && (TMath::Abs(pM2->fP.Eta()) < 2.5) 
//       && (TMath::Abs(pK1->fP.Eta()) < 2.5)
//       && (TMath::Abs(pK2->fP.Eta()) < 2.5)
//       ) {
//     ((TH1D*)fpHistFile->Get("efficiency"))->Fill(2); 
//     ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(3, "+ eta cuts"); 
    
//     if ((pM1->fP.Perp() > 1.0)  
// 	&& (pM2->fP.Perp() > 1.0)
// 	&& (pK1->fP.Perp() > 0.4)
// 	&& (pK2->fP.Perp() > 0.4)
// 	) {
//       ((TH1D*)fpHistFile->Get("efficiency"))->Fill(3); 
//       ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(4, "+ pT cuts"); 

//       if (m1Matched && m2Matched && k1Matched && k2Matched
// 	  && (prM1->fPlab.Perp() > 3.0) && (TMath::Abs(prM1->fPlab.Eta()) < 2.4)
// 	  && (prM2->fPlab.Perp() > 3.0) && (TMath::Abs(prM2->fPlab.Eta()) < 2.4)
// 	  && (prK1->fPlab.Perp() > 0.5) && (TMath::Abs(prK1->fPlab.Eta()) < 2.4)
// 	  && (prK2->fPlab.Perp() > 0.5) && (TMath::Abs(prK2->fPlab.Eta()) < 2.4)
// 	  && (prM1->fQ*prM2->fQ < 0)
// 	  && m1GT && m2GT && k1GT && k2GT
// 	  ) {
// 	((TH1D*)fpHistFile->Get("efficiency"))->Fill(4); 
// 	((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(5, "+ reco tracks/cuts"); 
	
// 	if (m1ID && m2ID) {
// 	  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(5); 
// 	  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(6, "+ muon ID"); 
	  
// 	  if (fGoodHLT) {
// 	    ((TH1D*)fpHistFile->Get("efficiency"))->Fill(6); 
// 	    ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(7, "+ trigger"); 
	    
// 	    if (pCand) {
// 	      ((TH1D*)fpHistFile->Get("efficiency"))->Fill(7); 
// 	      ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(8, "+ candidate"); 
// 	      fGoodEffCand = true;
// 	      ((TH1D*)fpHistFile->Get("effMass"))->Fill(pCand->fMass);
// 	    } 
// 	  }
// 	} 
//       } 
//     } 
//   } 

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
  int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

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

    if (!strcmp(CutName, "PHITYPE")) {
      PHITYPE = CutValue; ok = 1;
      if (dump) cout << "PHITYPE:      " << PHITYPE << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, PHITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #phi ID :: %d", CutName, PHITYPE));
    }

    if (!strcmp(CutName, "MKKLO")) {
      MKKLO = CutValue; ok = 1;
      if (dump) cout << "MKKLO:           " << MKKLO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(KK) :: %3.1f", CutName, MKKLO));
    }

    if (!strcmp(CutName, "MKKHI")) {
      MKKHI = CutValue; ok = 1;
      if (dump) cout << "MKKHI:           " << MKKHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(KK) :: %3.1f", CutName, MKKHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue; ok = 1;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }


    
  }

}
