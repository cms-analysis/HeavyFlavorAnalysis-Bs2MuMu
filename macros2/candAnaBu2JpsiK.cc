#include "candAnaBu2JpsiK.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

#include "../macros/AnalysisDistribution.hh"

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

  TAnaTrack *p0, *pk(0); 
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (321 == TMath::Abs(p0->fMCID)) {
      pk = p0; 
    }
    if (211 == TMath::Abs(p0->fMCID)) {
      pk = p0; 
    }
  }
  
  if (0 == pk) {
    cout << "candAnaBu2JpsiK::candAnalysis:  no kaon found " << endl;
    return;
  }

  fKaonPt        = pk->fRefPlab.Perp(); 
  fKaonEta       = pk->fRefPlab.Eta();  
  fKaonPhi       = pk->fRefPlab.Phi(); 
  fKaonTkQuality = highPurity(pk);
  fKaonPtNrf     = pk->fPlab.Perp();
  fKaonEtaNrf    = pk->fPlab.Eta();

  if (fCandTmi > -1) {
    TGenCand *pg1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(pk->fIndex)->getGenIndex());
    fKPtGen     = pg1->fP.Perp();
    fKEtaGen    = pg1->fP.Eta();
  } else {
    fKPtGen     = -99.;
    fKEtaGen    = -99.;
  }
  
  //cout<<" match kaon "<<fKaonPt<<endl;
  fKa1Missid = tightMuon(pk);  // true for tight  muons 
  fKa1MuMatch = doTriggerMatching(pk, true); // see if it matches HLT muon 

  if(0) { // for testing d.k.

    double mva=0;
    fKa1Missid2 = mvaMuon(pk,mva);  // true for tight  muons 
    
    fKa1MuMatch2 = doTriggerMatching(pk,false); // see if it matches HLT muon 

    fKa1MuMatchR = doTriggerMatchingR(pk,false); // matches to Bs/Jpsi-disp HLT muon 
    fKa1MuMatchR2 = doTriggerMatchingR(pk,true); // matches to fired HLT muon 

    fKa1MuMatchR5 = doTriggerMatchingR_OLD(pk,true); // matches to any "mu" HLT muon 
    //fKa1MuMatchR7 = doTriggerMatchingR_OLD(pk,false); // same as R
    
    fKa1MuMatchR3 = matchToMuon(pk,true); // matches muon, ignore self muon 
    //if(fKa1Missid) cout<<"missid "<<fKa1Missid<<" "<<fKa1MuMatch<<endl;
  } // end testing 




  // -- Check for J/psi mass
  //cout<<" check jpsi "<<endl;

  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  double chi2(0.);
  double ndof(0.), masse(0.); 
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
      masse = pD->fMassE;
      break;
    }
  }
 
  candAna::candAnalysis();
  // -- overwrite specific variables
  fCandTau   = fCandFL3d*MBPLUS/fCandP/TMath::Ccgs();
  fCandChi2  = chi2; 
  fCandDof   = ndof;
  fCandChi2Dof = chi2/ndof;
  fCandME      = masse;

  fPreselection = fPreselection && fGoodJpsiMass;
  fPreselection = fPreselection && fWideMass;

  if(0) { // special misid tests d.k.
    if( (pk->fIndex == fpMuon1->fIndex) || (pk->fIndex ==fpMuon2->fIndex) ) 
      cout<<" Kaon is a MUON "<<fEvt<<" "<<fpCand<<" "<<pk->fIndex<<" "<<fpMuon1->fIndex<<" "<<fpMuon2->fIndex<<" "<<fEvt<<endl;

    TVector3 trackMom = pk->fPlab;  // test track momentum
    TVector3 muonMom;
    muonMom = fpMuon1->fPlab;
    double dR1 = muonMom.DeltaR(trackMom);
    muonMom = fpMuon2->fPlab;
    double dR2 = muonMom.DeltaR(trackMom);

    if(dR1<dR2) { fKa1MuMatchR4 = dR1; fKa1MuMatchR6 = dR2;}
    else        { fKa1MuMatchR4 = dR2; fKa1MuMatchR6 = dR1;}
  } // end testing 


  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
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
  int nb(0), nphotons(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i); 
    if (521 == TMath::Abs(pC->fID)) {
      //cout<<" found B0 "<<TYPE<<endl;
      pB = pC;
      nb = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id); 
	if (443 == TMath::Abs(pC->fID)) {
	  // cout<<" found JPsi "<<endl;
	  ++nb;
	  pPsi = pC; 
	  pM1 = pM2 = 0; 
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++nphotons;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenTWithIndex(idd); 
	      } else {
		pM2 = fpEvt->getGenTWithIndex(idd); 
	      }
	      //cout<<" mu "<<pM1<<" "<<pM2<<endl;
	    }
	  }
	} else if ( (TYPE%1000)==66 && 211 == TMath::Abs(pC->fID)) { // get Bu2JpsiPi
	  ++nb;
	  pK = fpEvt->getGenTWithIndex(id); 
	  //cout<<" pi "<<pK<<" "<<nb<<endl;
	} else if (321 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pK = fpEvt->getGenTWithIndex(id); 
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
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fKGenID = pK->fID;
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


//   // -- check that only one reco track is matched to each gen particle
//   //    else skip the *event*!
//   static int cntBadEvents = 0; 
//   if (fGenM1Tmi > -1 && fGenM2Tmi > -1 && fGenK1Tmi > -1) {
//     TSimpleTrack *pT(0);
//     int cntM1(0), cntM2(0), cntK1(0); 
//     for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
//       pT = fpEvt->getSimpleTrack(i); 
//       if (fGenM1Tmi == pT->getGenIndex()) ++cntM1;
//       if (fGenM2Tmi == pT->getGenIndex()) ++cntM2;
//       if (fGenK1Tmi == pT->getGenIndex()) ++cntK1;
//     }
    
//     if (cntM1 > 1 || cntM2 > 1 || cntK1 > 1) {
//       cout << "BAD BAD event: multiple reco tracks matched to the same gen particle! " << ++cntBadEvents 
// 	   << ": " << cntM1 << " .. " << cntM2 << " .. " << cntK1
// 	   << " (gen-reco matches) " << endl;
//       fBadEvent = true; 
//     }
//   }
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = -1; 
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

    // -- kaon
    if (fGenK1Tmi > -1 && pT->getGenIndex() == fGenK1Tmi) {
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

  moreReducedTree(fTree);
  moreReducedTree(fAmsTree);

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");

}



// ----------------------------------------------------------------------
void candAnaBu2JpsiK::moreReducedTree(TTree *t) {
  // -- Additional reduced tree variables
  t->Branch("mpsi", &fJpsiMass,  "mpsi/D");
  t->Branch("psipt",&fJpsiPt,    "psipt/D");
  t->Branch("kpt",  &fKaonPt,    "kpt/D");
  t->Branch("keta", &fKaonEta,   "keta/D");
  t->Branch("kphi", &fKaonPhi,   "kphi/D");
  t->Branch("kgt",  &fKaonTkQuality,"kgt/I");
  t->Branch("t3pt", &fKaonPtNrf, "t3pt/D");
  t->Branch("t3eta",&fKaonEtaNrf,"t3eta/D");
  t->Branch("g3pt", &fKPtGen,    "g3pt/D");
  t->Branch("g3eta",&fKEtaGen,   "g3eta/D");
  t->Branch("g3id", &fKGenID,    "g3id/I"); 
  t->Branch("k1missid",  &fKa1Missid,    "k1missid/O");
  t->Branch("k1mumatch", &fKa1MuMatch,    "k1mumatch/O");

  if(0) { // testing d.k.
    t->Branch("k1missid2",  &fKa1Missid2,    "k1missid2/O");
    t->Branch("k1mumatch2", &fKa1MuMatch2,    "k1mumatch2/O");

    t->Branch("k1mumatchr", &fKa1MuMatchR,    "k1mumatchr/F");
    t->Branch("k1mumatchr2", &fKa1MuMatchR2,    "k1mumatchr2/F");
    t->Branch("k1mumatchr3", &fKa1MuMatchR3,    "k1mumatchr3/F");
    t->Branch("k1mumatchr4", &fKa1MuMatchR4,    "k1mumatchr4/F");
    t->Branch("k1mumatchr5", &fKa1MuMatchR5,    "k1mumatchr5/F");
    t->Branch("k1mumatchr6", &fKa1MuMatchR6,    "k1mumatchr6/F");
    //t->Branch("k1mumatchr7", &fKa1MuMatchR7,    "k1mumatchr7/F");
  }

}

// ----------------------------------------------------------------------
void candAnaBu2JpsiK::fillCandidateHistograms(int offset) {
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
  pB  = fpEvt->getGenTWithIndex(fGenBTmi); 
  pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi); 
  pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi); 
  pK  = fpEvt->getGenTWithIndex(fGenK1Tmi); 

  // -- reco level
  TSimpleTrack *prM1(0), *prM2(0), *prK(0); 
  double bla(0); 
  int m1Matched(0), m2Matched(0), kMatched(0), m1ID(0), m1tmID(0), m1mvaID(0), m2ID(0), m2tmID(0), m2mvaID(0), m1GT(0), m2GT(0), kGT(0);
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
    kMatched = 1; 
    prK = fpEvt->getSimpleTrack(fRecK1Tmi); 
    if (prK->getHighPurity()) {
      kGT = 1; 
    } else {
      kGT = 0;
    }
  } 

  m1ID = m1tmID; 
  m2ID = m2tmID; 

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
  if (kMatched) {
    fETk1pt  = prK->getP().Perp(); 
    fETk1eta = prK->getP().Eta(); 
    fETk1q   = prK->getCharge();
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
    
  }

}
