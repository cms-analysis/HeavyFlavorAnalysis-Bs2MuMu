#include "candAnaBu2JpsiK.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBu2JpsiK::candAnaBu2JpsiK(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fGenK1Tmi = fRecK1Tmi = -1; 
  BLIND = 0; 
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
    //    cout << it << " " << p0->fMCID << endl;
    if (321 == TMath::Abs(p0->fMCID)) {
      pks = p0; 
    }
    if (211 == TMath::Abs(p0->fMCID)) {
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

  if (fCandTmi > -1) {
    TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(pk->fIndex)->fGenIndex);
    fKPtGen     = pg1->fP.Perp();
    fKEtaGen    = pg1->fP.Eta();
  } else {
    fKPtGen     = -99.;
    fKEtaGen    = -99.;
  }
  
  fKa1Missid = tightMuon(pk);  // true for tight  muons 
  fKa1MuMatch = doTriggerMatching(pk); // see if it matches HLT muon 
  //if(fKa1Missid) cout<<"missid "<<fKa1Missid<<" "<<fKa1MuMatch<<endl;

  // -- Check for J/psi mass
  //cout<<" check jpsi "<<endl;

  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  double chi2(0.);
  int ndof(0); 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    //cout<<i<<" "<<pD->fType<<" "<<JPSITYPE<<endl;

    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp(); 
      fJpsiEta  = pD->fPlab.Eta(); 
      fJpsiPhi  = pD->fPlab.Phi(); 

      chi2 = pD->fVtx.fChi2;
      ndof = pD->fVtx.fNdof;
      break;
    }
  }
 
  candAna::candAnalysis();
  // -- overwrite specific variables
  fCandTau   = fCandFL3d*MBPLUS/fCandP/TMath::Ccgs();
  fCandChi2  = chi2; 
  fCandDof   = ndof;

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
      //cout<<" found B0 "<<TYPE<<endl;
      pB = pC;
      nb = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (443 == TMath::Abs(pC->fID)) {
	  // cout<<" found JPsi "<<endl;
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
	      //cout<<" mu "<<pM1<<" "<<pM2<<endl;
	    }
	  }
	} else if ( (TYPE%1000)==66 && 211 == TMath::Abs(pC->fID)) { // get Bu2JpsiPi
	  ++nb;
	  pK = fpEvt->getGenCand(id); 
	  //cout<<" pi "<<pK<<" "<<nb<<endl;
	} else if (321 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pK = fpEvt->getGenCand(id); 
	  //cout<<" K "<<pK<<" "<<nb<<endl;
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
	//cout<<" goodMatch "<<goodMatch<<endl;
	fNGenPhotons = nphotons;
	break;
      }
    }
  }

  fGenBTmi = -1; 
  if (goodMatch) {
    fGenBTmi = pB->fNumber; 
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = pB; 
    double x = (pM1->fV - pM->fV).Mag(); 
    fGenLifeTime = x*m/p/TMath::Ccgs();
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


  // -- check that only one reco track is matched to each gen particle
  //    else skip the *event*!
  static int cntBadEvents = 0; 
  if (fGenM1Tmi > -1 && fGenM2Tmi > -1 && fGenK1Tmi > -1) {
    TAnaTrack *pT(0);
    int cntM1(0), cntM2(0), cntK1(0); 
    for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
      pT = fpEvt->getRecTrack(i); 
      if (fGenM1Tmi == pT->fGenIndex) ++cntM1;
      if (fGenM2Tmi == pT->fGenIndex) ++cntM2;
      if (fGenK1Tmi == pT->fGenIndex) ++cntK1;
    }
    
    if (cntM1 > 1 || cntM2 > 1 || cntK1 > 1) {
      if (fVerbose > 1) cout << "BAD BAD event: multiple reco tracks matched to the same gen particle! " << ++cntBadEvents 
			     << ": " << cntM1 << " .. " << cntM2 << " .. " << cntK1
			     << " (gen-reco matches) " << endl;
      fBadEvent = true; 
    }
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

  //cout<<TYPE<<endl;

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 

    //cout<<TYPE<<" "<<pCand->fType<<" "<<iC<<endl;

    if (TYPE != pCand->fType) continue;
    
    d1Matched = d2Matched = d3Matched = 0; 
    //cout<< TYPE;
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
      //cout<<" "<<idx<<" "<<type;
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
    
    //cout<<endl;

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
  fTree->Branch("k1missid",  &fKa1Missid,    "k1missid/O");

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");

//  string name; 
//  int i(0); 
//  for (map<string, int>::iterator imap = fRegion.begin(); imap != fRegion.end(); ++imap) {  
//    i    = imap->second; 
//    name = imap->first + "_";
//
//    //    cout << "  booking analysis distributions for " << name << " at offset i = " << i << endl;
//
//    fpMpsi[i]     = bookDistribution(Form("%smpsi", name.c_str()), "m(J/#psi) [GeV]", "fGoodJpsiMass", 40, 2.8, 3.4);           
//    fpKaonPt[i]   = bookDistribution(Form("%skaonpt", name.c_str()), "p_{T, K} [GeV]", "fGoodTracksPt", 25, 0., 25.);           
//    fpKaonEta[i]  = bookDistribution(Form("%skaoneta", name.c_str()), "#eta_{K}", "fGoodTracksEta", 25, -2.5, 2.5);
//    fpPsiPt[i]    = bookDistribution(Form("%spsipt", name.c_str()), "p_{T, J/#psi} [GeV]", "fGoodTracksPt", 25, 0., 50.);           
//    fpPsiEta[i]   = bookDistribution(Form("%spsieta", name.c_str()), "#eta_{J/#psi}", "fGoodTracksEta", 25, -2.5, 2.5);  
//  }


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

  fEffTree->Fill(); 

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
