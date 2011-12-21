#include "candAnaBu2JpsiK.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBu2JpsiK::candAnaBu2JpsiK(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fGenK1Tmi = fRecK1Tmi = -1; 

  cout << "==> candAnaBu2JpsiK: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaBu2JpsiK::~candAnaBu2JpsiK() {
  cout << "==> candAnaBu2JpsiK: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::candAnalysis() {

  if (0 == fpCand) return;

  TAnaTrack *p0, *pk(0), *pks(0); 
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (321 == TMath::Abs(p0->fMCID)) {
      pks = p0; 
    }
  }

  if (0 == pks) {
    cout << "candAnaBu2JpsiK::candAnalysis:  no kaon found " << endl;
    return;
  }

  pk = fpEvt->getRecTrack(pks->fIndex);

  fKaonPt        = pk->fPlab.Perp(); 
  fKaonEta       = pk->fPlab.Eta();  
  fKaonPhi       = pk->fPlab.Phi(); 
  fKaonTkQuality = pk->fTrackQuality & TRACKQUALITY;
  fKaonPtNrf     = pks->fPlab.Perp();
  fKaonEtaNrf    = pks->fPlab.Eta();

//FIXME  
//   if (fIsMC) {
//     fGenBpartial = partialReco(fpCand); 
//     if (fpCand->fMass > 5.5 && fGenBpartial == 511) fpEvt->dumpGenBlock();
//   }

  if (fCandTmi > -1) {
    TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(pk->fIndex)->fGenIndex);
    fKPtGen     = pg1->fP.Perp();
    fKEtaGen    = pg1->fP.Eta();
  } else {
    fKPtGen     = -99.;
    fKEtaGen    = -99.;
  }
  

  // -- Check for J/psi mass
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp(); 
      fJpsiEta  = pD->fPlab.Eta(); 
      fJpsiPhi  = pD->fPlab.Phi(); 
      //       cout << "type = " << pD->fType 
      //        	   << " with mass = " << pD->fMass 
      //        	   << " fGoodJpsiMass = " << fGoodJpsiMass 
      //        	   << endl;
      break;
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
//       cout << "bmmNormalizationReader::fillCandidateVariables:  no muon 1 found " << endl;
//       return; 
//     }
//     if (0 == p2) {
//       cout << "bmmNormalizationReader::fillCandidateVariables:  no muon 2 found " << endl;
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

  
  candAna::candAnalysis();
  fPreselection = fPreselection && fGoodJpsiMass;

  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(3); 

}

// ----------------------------------------------------------------------
void candAnaBu2JpsiK::moreBasicCuts() {
  cout << "   candAnaBu2JpsiK: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::genMatch() {

  fGenM1Tmi = fGenM2Tmi = fGenK1Tmi = -1; 
  fNGenPhotons = 0; 

  TGenCand *pC(0), *pM1(0), *pM2(0), *pK(0), *pB(0), *pPsi(0); 
  int nb(0), npsi(0), nphotons(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (521 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (443 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pPsi = pC; 
	  npsi = pPsi->fDau2 - pPsi->fDau1 + 1; 
	  pM1 = pM2 = 0; 
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++nphotons;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenCand(idd); 
	      } else {
		pM2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	} else if (321 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pK = fpEvt->getGenCand(id); 
	} else 	if (22 == TMath::Abs(pC->fID)) {
	  ++nphotons;
	} else {
	  ++nb;
	}
      }
      if (nb > 2) {
	pM1 = pM2 = pK = pPsi = 0; 
	continue; // skip B decays where more than J/psi and kaon came from B
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK && (pPsi->fMom1 == pK->fMom1)) {
	goodMatch = true; 
	fNGenPhotons = nphotons;
	break;
      }
    }
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
    fGenK1Tmi = pK->fNumber; 
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
    fGenK1Tmi = -1;  
  }
  
  if (fVerbose > 10) {
    cout << "fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fGenK1Tmi = " << fGenK1Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = -1; 
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

    // -- kaon
    if (fGenK1Tmi > -1 && pT->fGenIndex == fGenK1Tmi) {
      fRecK1Tmi = i; 
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecK1Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fRecK1Tmi = " << fRecK1Tmi << " matched to fGenK1Tmi = " << fGenK1Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::candMatch() {

  fCandTmi = -1;   
  int idx(-1), type(-1); 
  int d1Matched(0), d2Matched(0), d3Matched(0); 
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 
    if (TYPE != pCand->fType) continue;
    
    d1Matched = d2Matched = d3Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
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
    }
    
    if (d1Matched && d2Matched && d3Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
  }

}




// ----------------------------------------------------------------------
void candAnaBu2JpsiK::bookHist() {
  cout << "==>candAnaBu2JpsiK: bookHist" << endl;
  candAna::bookHist();

  // -- Additional reduced tree variables
  fTree->Branch("mpsi", &fJpsiMass,  "mpsi/D");
  fTree->Branch("psipt",&fJpsiPt,    "psipt/D");
  fTree->Branch("kpt",  &fKaonPt,    "kpt/D");
  fTree->Branch("keta", &fKaonEta,   "keta/D");
  fTree->Branch("kphi", &fKaonPhi,   "kphi/D");
  fTree->Branch("kgt",  &fKaonTkQuality,"kgt/I");
  fTree->Branch("t3pt", &fKaonPtNrf, "t3pt/D");
  fTree->Branch("t3eta",&fKaonEtaNrf,"t3eta/D");
  fTree->Branch("g3pt", &fKPtGen,    "g3pt/D");
  fTree->Branch("g3eta",&fKEtaGen,   "g3eta/D");

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");

  string name; 
  int i(0); 
  for (map<string, int>::iterator imap = fRegion.begin(); imap != fRegion.end(); ++imap) {  
    i    = imap->second; 
    name = imap->first + "_";

    //    cout << "  booking analysis distributions for " << name << " at offset i = " << i << endl;

    fpMpsi[i]     = bookDistribution(Form("%smpsi", name.c_str()), "m(J/#psi) [GeV]", "fGoodJpsiMass", 40, 2.8, 3.4);           
    fpKaonPt[i]   = bookDistribution(Form("%skaonpt", name.c_str()), "p_{T, K} [GeV]", "fGoodTracksPt", 25, 0., 25.);           
    fpKaonEta[i]  = bookDistribution(Form("%skaoneta", name.c_str()), "#eta_{K}", "fGoodTracksEta", 25, -2.5, 2.5);
    fpPsiPt[i]    = bookDistribution(Form("%spsipt", name.c_str()), "p_{T, J/#psi} [GeV]", "fGoodTracksPt", 25, 0., 50.);           
    fpPsiEta[i]   = bookDistribution(Form("%spsieta", name.c_str()), "#eta_{J/#psi}", "fGoodTracksEta", 25, -2.5, 2.5);  
  }


}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::fillCandidateHistograms(int offset) {

  fpMpsi[offset]->fill(fJpsiMass, fCandM); 
  fpTracksPt[offset]->fill(fKaonPt, fCandM);
  fpTracksEta[offset]->fill(fKaonEta, fCandM);

  fpKaonPt[offset]->fill(fKaonPt, fCandM);
  fpKaonEta[offset]->fill(fKaonEta, fCandM);

  fpPsiPt[offset]->fill(fJpsiPt, fCandM); 
  fpPsiEta[offset]->fill(fJpsiEta, fCandM); 

  candAna::fillCandidateHistograms(offset); 

}



// ----------------------------------------------------------------------
void candAnaBu2JpsiK::efficiencyCalculation() {
  fGoodEffCand = false;

//   ((TH1D*)fpHistFile->Get("efficiency"))->Fill(0); 
//   ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(1, "all events"); 


  // -- gen level 
  TGenCand *pB(0), *pM1(0), *pM2(0), *pK(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenK1Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenCand(fGenBTmi); 
  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 
  pK  = fpEvt->getGenCand(fGenK1Tmi); 

//   ((TH1D*)fpHistFile->Get("efficiency"))->Fill(1); 
//   ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(2, "gen signal decays"); 


  // -- reco level
  TAnaTrack *prM1(0), *prM2(0), *prK(0); 
  int m1Matched(0), m2Matched(0), kMatched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0), kGT(0);
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
    kMatched = 1; 
    prK = fpEvt->getRecTrack(fRecK1Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK->fTrackQuality & TRACKQUALITY))) {
      kGT = 0; 
    } else {
      kGT = 1;
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
  fETg3pt  = pK->fP.Perp(); 
  fETg3eta = pK->fP.Eta(); 
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
  if (kMatched) {
    fETk1pt  = prK->fPlab.Perp(); 
    fETk1eta = prK->fPlab.Eta(); 
    fETk1q   = prK->fQ;
    fETk1gt  = (kGT>0?true:false); 
  } else {
    fETk1pt  = -99.; 
    fETk1eta = -99.; 
    fETk1q   = -99;
    fETk1gt  = false; 
  }
  if (pCand) {
    fETcandMass = pCand->fMass; 
  } else {
    fETcandMass = -99.;
  }

  //   cout << Form("%4d", fEvent) << " bmmNormalizationReader: m = " << fETcandMass << " from cand " << pCand 
  //        << " mu gen: " << fGenM1Tmi << " " << fGenM2Tmi 
  //        << " mu gen: " << fETg1pt << " " << fETg2pt 
  //        << endl;

  fEffTree->Fill(); 


//   // -- results...
//   if ((TMath::Abs(pM1->fP.Eta()) < 2.5) && (TMath::Abs(pM2->fP.Eta()) < 2.5) && (TMath::Abs(pK->fP.Eta()) < 2.5)) {
//     ((TH1D*)fpHistFile->Get("efficiency"))->Fill(2); 
//     ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(3, "+ eta cuts"); 
    
//     if ((pM1->fP.Perp() > 1.) && (pM2->fP.Perp() > 1.) && (pK->fP.Perp() > 0.4)) {
//       ((TH1D*)fpHistFile->Get("efficiency"))->Fill(3); 
//       ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(4, "+ pT cuts"); 
      
//       if (m1Matched && m2Matched && kMatched
// 	  && (prM1->fPlab.Perp() > 3.0) && (prM2->fPlab.Perp() > 3.0) && (prK->fPlab.Perp() > 0.5)
// 	  && (TMath::Abs(prM1->fPlab.Eta()) < 2.4) && (TMath::Abs(prM2->fPlab.Eta()) < 2.4) && (TMath::Abs(prK->fPlab.Eta()) < 2.4)
// 	  && (prM1->fQ*prM2->fQ < 0)
// 	  && m1GT && m2GT && kGT
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
void candAnaBu2JpsiK::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaBu2JpsiK: Reading " << fCutFile << " for cut settings" << endl;
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
    
  }

}
