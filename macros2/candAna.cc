#include "candAna.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAna::candAna(bmm2Reader *pReader, string name, string cutsFile) {
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", reading cutsfile " << cutsFile << endl;
  fpReader = pReader; 
  fName = name; 

  MASSMIN = 4.5;
  MASSMAX = 6.5; 

  fHistDir = gFile->mkdir(fName.c_str());
  //   fHistDir->cd();
  //   cout << "pwd(): "; fHistDir->pwd();

  fRegion.insert(make_pair("A", 0)); // all
  fRegion.insert(make_pair("B", 1)); // barrel
  fRegion.insert(make_pair("E", 2)); // endcap

  fRegion.insert(make_pair("AR0", 3)); // run range 0: HLT_DoubleMu3_Jpsi_v1 
                                       //              HLT_DoubleMu3_Jpsi_v2
  fRegion.insert(make_pair("AR1", 4)); // run range 1: HLT_Dimuon6p5_Jpsi_Displaced_v1
  fRegion.insert(make_pair("AR2", 5)); // run range 2: HLT_Dimuon7_Jpsi_Displaced_v1 
                                       //              HLT_Dimuon7_Jpsi_Displaced_v3
  fRegion.insert(make_pair("AR3", 6)); // run range 3: HLT_DoubleMu3p5_Jpsi_Displaced_v2
  fRegion.insert(make_pair("AR4", 7)); // run range 4: HLT_DoubleMu5_Jpsi_Displaced_v1 
                                       //              HLT_DoubleMu5_Jpsi_Displaced_v2
  fRegion.insert(make_pair("AR5", 8)); // run range 5: HLT_DoubleMu5_Jpsi_Displaced_v4
  fRegion.insert(make_pair("AR6", 9)); // should be empty


}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}



// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {
  fpEvt = evt; 

  if (fIsMC) {
    processType();
    genMatch(); 
    recoMatch(); 
    candMatch(); 
  }

  triggerSelection();

  fRunRange = 6;
  if ((fRun >= 160329) && (fRun <= 163261)) {
    fRunRange = 0; 
  } 
  if ((fRun >= 163269) && (fRun <= 163869)) {
    fRunRange = 1; 
  } 
  if ((fRun >= 165088) && (fRun <= 167913)) {
    fRunRange = 2; 
  } 
  if ((fRun >= 170249) && (fRun <= 173198)) {
    fRunRange = 3; 
  } 
  if ((fRun >= 173236) && (fRun <= 178380)) {
    fRunRange = 4; 
  } 
  if ((fRun >= 178420) && (fRun <= 999999)) {
    fRunRange = 5; 
  } 


  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (fVerbose > 29) cout << "candidate at " << iC << " which is of type " << pCand->fType << endl;
    if (TYPE != pCand->fType) {
      if (fVerbose > 19) cout << "Skipping candidate at " << iC << " which is of type " << pCand->fType << endl;
      continue;
    }
    if (fVerbose > 2) cout << "Analyzing candidate at " << iC << " which is of type " << TYPE 
			   << " with sig tracks: " << pCand->fSig1 << " .. " << pCand->fSig2
			   << " and rec tracks: " 
			   << fpEvt->getSigTrack(pCand->fSig1)->fIndex
			   << " .. " 
			   << fpEvt->getSigTrack(pCand->fSig2)->fIndex
			   << endl;

    fpCand = pCand;
    fCandIdx = iC; 
    // -- call derived functions
    candAnalysis();
    if (fIsMC) {
      fTree->Fill(); 
    } else {
      if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX) {
	// do nothing
      } else {
	if (fPreselection) {
	  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(12); 
	  fTree->Fill(); 
	}         
      }
    }

    // -- fill histograms
    fillCandidateHistograms(fRegion[Form("AR%i", fRunRange)]);
    if (fRunRange > 0) {
      fillCandidateHistograms(fRegion["A"]);
      if (fBarrel) {
	fillCandidateHistograms(fRegion["B"]);
      } else {
	fillCandidateHistograms(fRegion["E"]);
      }
    }
  }

}

// ----------------------------------------------------------------------
void candAna::candAnalysis() {

  if (0 == fpCand) return;

  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(1); 

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
    fPvNtrk = pv->getNtracks();
  } else {
    fPvX = -99.;
    fPvY = -99.;
    fPvZ = -99.;
    fPvNtrk = -99;
  }

  if (fIsMC) {
    if (fCandTmi == fCandIdx) {
      if (fNGenPhotons) {
	fCandTM = 2; 
      } else {
	fCandTM = 1; 
      }
    } else {
      fCandTM = 0; 
    }
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
  
  // -- find mass-constrained candidate with the same rectracks as this one
  fCandM2 = -1.;
  vector<int> rectracks;
  int bla;
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    bla = fpEvt->getSigTrack(it)->fIndex; 
    rectracks.push_back(bla);
  } 

  TAnaCand *pC(0), *pDau(0), *pMcCand(0); 
  int mctype = TYPE + 100000; 
  unsigned int nmatch(0); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pC = fpEvt->getCand(iC);
    if (pC->fType != mctype) continue;
    nmatch = 0; 
    for (int it = pC->fSig1; it <= pC->fSig2; ++it) {
      bla = fpEvt->getSigTrack(it)->fIndex; 
      if (rectracks.end() != find(rectracks.begin(), rectracks.end(), bla)) {
	++nmatch; 
      }
    }
    if (pC->fDau1 > 0) {
      pDau = fpEvt->getCand(pC->fDau1); 
      for (int it = pDau->fSig1; it <= pDau->fSig2; ++it) {
	bla = fpEvt->getSigTrack(it)->fIndex; 
	if (rectracks.end() != find(rectracks.begin(), rectracks.end(), bla)) {
	  ++nmatch; 
	}
      }
    }
    if (nmatch == rectracks.size()) {
      pMcCand = pC; 
      break;
    }
  }
  if (pMcCand) {
    fCandM2 = pMcCand->fMass;
  } else {
    fCandM2 = -2.;
  }

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
  if (fpCand->fType > 1000000 && fpCand->fType < 2000000) {
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

  fMu1Id        = goodMuon(p1); 
  fMu1Pt        = p1->fPlab.Perp(); 
  fMu1Eta       = p1->fPlab.Eta(); 
  fMu1Phi       = p1->fPlab.Phi(); 
  fMu1PtNrf     = ps1->fPlab.Perp();
  fMu1EtaNrf    = ps1->fPlab.Eta();
  fMu1TkQuality = p1->fTrackQuality & TRACKQUALITY;
  fMu1Q         = p1->fQ;
  fMu1Pix       = fpReader->numberOfPixLayers(p1);
  fMu1BPix      = fpReader->numberOfBPixLayers(p1);
  fMu1BPixL1    = fpReader->numberOfBPixLayer1Hits(p1);

  if (p1->fMuIndex > -1) {
    fMu1Chi2      = fpEvt->getMuon(p1->fMuIndex)->fMuonChi2;
  } else {
    fMu1Chi2 = -98.;
  }

  if (fCandTM && fGenM1Tmi < 0) fpEvt->dump();

  //  cout << "  " << p1->fGenIndex << "  " << fCandTM << " fCandTmi: " << fCandTmi << " gen m " << fGenM1Tmi << " " << fGenM2Tmi << endl;

  if (fCandTM) {
    TGenCand *pg1 = fpEvt->getGenCand(p1->fGenIndex);
    fMu1PtGen     = pg1->fP.Perp();
    fMu1EtaGen    = pg1->fP.Eta();
  } else {
    fMu1PtGen     = -99.;
    fMu1EtaGen    = -99.;
  }
  
  fMu2Id        = goodMuon(p2); 
  fMu2Pt        = p2->fPlab.Perp(); 
  fMu2Eta       = p2->fPlab.Eta(); 
  fMu2Phi       = p2->fPlab.Phi(); 
  fMu2PtNrf     = ps2->fPlab.Perp();
  fMu2EtaNrf    = ps2->fPlab.Eta();
  fMu2TkQuality = p2->fTrackQuality & TRACKQUALITY;
  fMu2Q         = p2->fQ;
  fMu2Pix       = fpReader->numberOfPixLayers(p2);
  fMu2BPix      = fpReader->numberOfBPixLayers(p2);
  fMu2BPixL1    = fpReader->numberOfBPixLayer1Hits(p2);

  if  (fMu1Id && fMu2Id) {
    //  if (p1->fMuIndex > -1 && p2->fMuIndex > -1) {
    TVector3 rm1  = fpEvt->getMuon(p1->fMuIndex)->fPositionAtM2;
    TVector3 rm2  = fpEvt->getMuon(p2->fMuIndex)->fPositionAtM2;

    if (rm1.Mag() > 0.1 && rm2.Mag() > 0.1) {
      TVector3 rD   = rm2-rm1; 
      fMuDist   = rD.Mag(); 
      fMuDeltaR =  rm1.DeltaR(rm2); 
    } else {
      fMuDist   = -99.; 
      fMuDeltaR = -99.; 
    }

    if (fMuDeltaR > 100) {
      cout << "r1 = (" << rm1.X() << ", " << rm1.Y() << ", " << rm1.Z() 
	   << ") r2 = (" << rm2.X() << ", " << rm2.Y() << ", " << rm2.Z() << ")"
	   << endl;
    }
    if (fVerbose > 10) cout << "dist: " << fMuDist << " dr = " << fMuDeltaR << endl;
  } else {
    fMuDist   = -99.; 
    fMuDeltaR = -99.; 
  }

  if (p2->fMuIndex > -1) {
    fMu2Chi2      = fpEvt->getMuon(p2->fMuIndex)->fMuonChi2;
  } else {
    fMu2Chi2 = -98.;
  }

  if ((TMath::Abs(fMu1Eta) < 1.4) && (TMath::Abs(fMu2Eta) < 1.4)) {
    fBarrel = true; 
  } else {
    fBarrel = false; 
  }    

  double dphi = p1->fPlab.DeltaPhi(p2->fPlab);
  fCowboy = (p1->fQ*dphi > 0); 

  // -- Muon weights
  PidTable *pT, *pT1, *pT2; 
  if (fCowboy) {
    pT = fpReader->ptCbMUID; 
  } else {
    pT = fpReader->ptSgMUID; 
  }
  fMu1W8Mu      = pT->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi);
  fMu2W8Mu      = pT->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi);
  
  if (fCowboy) {
    pT1 = fpReader->ptCbMUT1;
    pT2 = fpReader->ptCbMUT2;
  } else {
    pT1 = fpReader->ptSgMUT1;
    pT2 = fpReader->ptSgMUT2;
  }
  fMu1W8Tr      = pT1->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi)*pT2->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi);
  fMu2W8Tr      = pT1->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi)*pT2->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi);


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

  // -- FIXME ????
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0); 
  // -- this is from the full candidate
  TAnaVertex sv = fpCand->fVtx;
  // -- this is from the dimuon vertex
  TAnaCand *pD; 
  //  cout << "looking at daughters " << fpCand->fDau1  << " .. " << fpCand->fDau2 << endl;
  for (int id = fpCand->fDau1; id <= fpCand->fDau2; ++id) {
    if (id < 0) break;
    pD = fpEvt->getCand(id); 
    //    cout << "looking at daughter " <<  id << " with type = " << pD->fType << endl;
    if (300443 == pD->fType) {
      sv = pD->fVtx;
      //      cout << "  Found J/psi vertex" << endl;
      break;
    }
  }

  // -- go back to original!
  sv = fpCand->fVtx;


  TVector3 svpv(sv.fPoint - fpEvt->getPV(pvidx)->fPoint); 
  double alpha = svpv.Angle(fpCand->fPlab);
  fCandCosA   = TMath::Cos(alpha);
  fCandA      = alpha; 

  double iso  = isoClassicWithDOCA(fpCand, 0.05); // 500um DOCA cut
  fCandIsoTrk = fCandItrk;
  fCandIso    = iso; 
  
  // -- the matrix
//   fIsoR05Pt03 = isoWithDOCA(fpCand, 0.05, 0.5, 0.3);
//   fIsoR05Pt05 = isoWithDOCA(fpCand, 0.05, 0.5, 0.5);
//   fIsoR05Pt07 = isoWithDOCA(fpCand, 0.05, 0.5, 0.7);
//   fIsoR05Pt09 = isoWithDOCA(fpCand, 0.05, 0.5, 0.9);
//   fIsoR05Pt11 = isoWithDOCA(fpCand, 0.05, 0.5, 1.1);

//   fIsoR07Pt03 = isoWithDOCA(fpCand, 0.05, 0.7, 0.3);
//   fIsoR07Pt05 = isoWithDOCA(fpCand, 0.05, 0.7, 0.5);
//   fIsoR07Pt07 = isoWithDOCA(fpCand, 0.05, 0.7, 0.7);
//   fIsoR07Pt09 = isoWithDOCA(fpCand, 0.05, 0.7, 0.9);
//   fIsoR07Pt11 = isoWithDOCA(fpCand, 0.05, 0.7, 1.1);

//   fIsoR10Pt03 = isoWithDOCA(fpCand, 0.05, 1.0, 0.3);
//   fIsoR10Pt05 = isoWithDOCA(fpCand, 0.05, 1.0, 0.5);
//   fIsoR10Pt07 = isoWithDOCA(fpCand, 0.05, 1.0, 0.7);
//   fIsoR10Pt09 = isoWithDOCA(fpCand, 0.05, 1.0, 0.9);
//   fIsoR10Pt11 = isoWithDOCA(fpCand, 0.05, 1.0, 1.1);

  fCandChi2  = sv.fChi2;
  fCandDof   = sv.fNdof;
  fCandProb  = sv.fProb;
  fCandFL3d  = sv.fD3d;
  fCandFL3dE = sv.fD3dE;
  fCandFLS3d = sv.fD3d/sv.fD3dE; 
  if (TMath::IsNaN(fCandFLS3d)) fCandFLS3d = -1.;
  fCandFLSxy = sv.fDxy/sv.fDxyE; 
  if (TMath::IsNaN(fCandFLSxy)) fCandFLSxy = -1.;

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
  fMu1IP        = sv.fD3d*bmu1/t1->fTip;
  fMu2IP        = sv.fD3d*bmu2/t2->fTip;

  // -- fill cut variables
  fWideMass = ((fpCand->fMass > MASSMIN) && (fpCand->fMass < MASSMAX)); 

  fGoodMuonsID   = (fMu1Id && fMu2Id);
  fGoodMuonsPt   = ((fMu1Pt > MUPTLO) && (fMu1Pt < MUPTHI) && (fMu2Pt > MUPTLO) && (fMu2Pt < MUPTHI));
  fGoodMuonsEta  = ((fMu1Eta > MUETALO) && (fMu1Eta < MUETAHI) && (fMu2Eta > MUETALO) && (fMu2Eta < MUETAHI));
  fGoodTracks    = (goodTrack(p1) && goodTrack(p2));
  fGoodTracksPt  = ((fMu1Pt > TRACKPTLO) && (fMu1Pt < TRACKPTHI) && (fMu2Pt > TRACKPTLO) && (fMu2Pt < TRACKPTHI));
  fGoodTracksEta = ((fMu1Eta > TRACKETALO) && (fMu1Eta < TRACKETAHI) && (fMu2Eta > TRACKETALO) && (fMu2Eta < TRACKETAHI));

  fGoodQ = (fMu1Q*fMu2Q < 0); 
  fGoodPt = (fCandPt > CANDPTLO);
  fGoodEta = ((fCandEta > CANDETALO) && (fCandEta < CANDETAHI)); 
  fGoodAlpha = (fCandA < CANDALPHA); 
  fGoodIso = (fCandIso > CANDISOLATION); 
  fGoodChi2 = (fCandChi2/fCandDof < CANDVTXCHI2);
  fGoodFLS =  ((fCandFLS3d > CANDFLS3D) && (fCandFLSxy > CANDFLSXY)); 
  if (TMath::IsNaN(fCandFLS3d)) fGoodFLS = false;

  fGoodDocaTrk = (fCandDocaTrk > CANDDOCATRK);
  fGoodLastCut = true; 

  fAnaCuts.update(); 

  fPreselection = fWideMass && fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsID && fGoodMuonsPt && fGoodMuonsEta; 
  fPreselection = fPreselection && fGoodQ && (fCandPt > 4) && (fCandA < 0.4) && (fCandFLS3d > 2) && (fCandChi2/fCandDof < 10); 

  //  fPreselection = true; 
}



// ----------------------------------------------------------------------
void candAna::fillCandidateHistograms(int offset) {
 
  // -- only candidate histograms below
  if (0 == fpCand) return;
  
  // -- Fill distributions
  fpTracksPt[offset]->fill(fMu1Pt, fCandM);
  fpTracksPt[offset]->fill(fMu2Pt, fCandM);
  fpTracksEta[offset]->fill(fMu1Eta, fCandM);
  fpTracksEta[offset]->fill(fMu2Eta, fCandM);

  fpMuonsPt[offset]->fill(fMu1Pt, fCandM);
  fpMuonsPt[offset]->fill(fMu2Pt, fCandM);
  fpMuonsEta[offset]->fill(fMu1Eta, fCandM);
  fpMuonsEta[offset]->fill(fMu2Eta, fCandM);
  fpMuon1Eta[offset]->fill(fMu1Eta, fCandM);
  fpMuon2Eta[offset]->fill(fMu2Eta, fCandM);
  fpMuon1Pt[offset]->fill(fMu1Pt, fCandM);
  fpMuon2Pt[offset]->fill(fMu2Pt, fCandM);

  fpPvZ[offset]->fill(fPvZ, fCandM); 
  fpPvN[offset]->fill(fPvN, fCandM); 
  fpPvNtrk[offset]->fill(fPvNtrk, fCandM); 
  fpPt[offset]->fill(fCandPt, fCandM); 
  fpEta[offset]->fill(fCandEta, fCandM); 
  fpAlpha[offset]->fill(fCandA, fCandM);
  fpCosA[offset]->fill(fCandCosA, fCandM);
  fpIso[offset]->fill(fCandIso, fCandM);

  fpChi2[offset]->fill(fCandChi2, fCandM);
  fpChi2Dof[offset]->fill(fCandChi2/fCandDof, fCandM); 
  fpProb[offset]->fill(fCandProb, fCandM);   
  fpFLS3d[offset]->fill(fCandFLS3d, fCandM); 
  fpFL3d[offset]->fill(fCandFL3d, fCandM); 
  fpFL3dE[offset]->fill(fCandFL3dE, fCandM); 
  fpFLSxy[offset]->fill(fCandFLSxy, fCandM); 
  fpDocaTrk[offset]->fill(fCandDocaTrk, fCandM); 

  int ipv = (fPvN<30?fPvN/2:NADPV-1); 
  if (fpNpvPvN[ipv][offset]) {
    fpNpvPvN[ipv][offset]->fill(fPvN, fCandM); 
    fpNpvChi2Dof[ipv][offset]->fill(fCandChi2/fCandDof, fCandM); 
    fpNpvProb[ipv][offset]->fill(fCandProb, fCandM);   
    fpNpvFLS3d[ipv][offset]->fill(fCandFLS3d, fCandM); 
    fpNpvFLSxy[ipv][offset]->fill(fCandFLSxy, fCandM); 
    fpNpvDocaTrk[ipv][offset]->fill(fCandDocaTrk, fCandM); 
    fpNpvIso[ipv][offset]->fill(fCandIso, fCandM);
    fpNpvIsoTrk[ipv][offset]->fill(fCandIsoTrk, fCandM);
  }
}


// ----------------------------------------------------------------------
void candAna::basicCuts() {
  cout << "    candAna basic cuts" << endl;
  fAnaCuts.addCut("fWideMass", "m(B candidate) [GeV]", fWideMass); 
  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT); 
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID); 
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt); 
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta); 
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks); 
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt); 
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta); 
}


// ----------------------------------------------------------------------
void candAna::moreBasicCuts() {
  cout << "    candAna more basic cuts?" << endl;

}


// ----------------------------------------------------------------------
void candAna::candidateCuts() {
  cout << "    candAna candidate cuts" << endl;
  fAnaCuts.addCut("fGoodQ", "q_{1} 1_{2}", fGoodQ); 
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt); 
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta); 
  fAnaCuts.addCut("fGoodAlpha", "#alpha", fGoodAlpha); 
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS); 
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2); 
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso); 
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk); 
  fAnaCuts.addCut("fGoodLastCut", "lastCut", fGoodLastCut); 
}


// ----------------------------------------------------------------------
void candAna::moreCandidateCuts() {
  cout << "    candAna more candidate cuts?" << endl;

}


// ---------------------------------------------------------------------- 
AnalysisDistribution* candAna::bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn, ht, nbins, lo, hi); 
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX); 
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX); 
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX); 
  p->setAnalysisCuts(&fAnaCuts, hc); 
  p->setPreselCut(&fPreselection); 

  TH1 *h = (TH1D*)fHistDir->Get("analysisDistributions"); 
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
bool candAna::goodTrack(TAnaTrack *pt) {

  if (TRACKQUALITY > 0 && (0 == (pt->fTrackQuality & TRACKQUALITY))) {
    if (fVerbose > 5) cout << "track " << pt->fIndex << " failed track quality: " << pt->fTrackQuality << endl;
    return false; 
  }
  
  if (TMath::Abs(pt->fTip) > TRACKTIP) {
    if (fVerbose > 5) cout << "track " << pt->fIndex << " failed tip: " << pt->fTip
			   << " pointing to PV = "  << pt->fPvIdx  << endl;
    return false; 
  }
  
  if (TMath::Abs(pt->fLip) > TRACKLIP) { 
    if (fVerbose > 5) cout << "track " << pt->fIndex << " failed lip: " << pt->fLip 
			   << " pointing to PV = "  << pt->fPvIdx  << endl;          
    return false; 
  }

  return true; 
}


// ----------------------------------------------------------------------
void candAna::processType() {

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
void candAna::genMatch() {
  cout << "candAna::genMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::recoMatch() {
  cout << "candAna::recoMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::candMatch() {
  cout << "candAna::candMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::triggerSelection() {

  fGoodHLT = false; 
  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 
  //  cout << " ----------------------------------------------------------------------" << endl;

  //   FIXME
  //   if (HLTPath.end() != find(HLTPath.begin(), HLTPath.end(), "NOTRIGGER"))  {
  //     //    cout << "NOTRIGGER requested... " << endl;
  //     fGoodHLT = true; 
  //     return;
  //   }
  
  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 
    
    //    if (-32 == fVerbose) cout << "path " << i << ": " << a << endl;
    if (wasRun && result) {
      if (-32 == fVerbose) {
	if ((a != "digitisation_step") 
	    && (a != "L1simulation_step") 
	    && (a != "digi2raw_step") 
	    && (a != "HLTriggerFinalPath") 
	    && (a != "raw2digi_step") 
	    && (a != "reconstruction_step") 
	    ) {
	  cout << "run and fired: " << a << endl;
	}
      }

      string spath; 
      int rmin, rmax; 
      for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin(); imap != HLTRANGE.end(); ++imap) {  
	spath = imap->first; 
	rmin = imap->second.first; 
	rmax = imap->second.second; 
	// 	if ((a != "digitisation_step") 
	// 	    && (a != "L1simulation_step") 
	// 	    && (a != "digi2raw_step") 
	// 	    && (a != "HLTriggerFinalPath") 
	// 	    && (a != "raw2digi_step") 
	// 	    && (a != "reconstruction_step") 
	// 	    ) {
	// 	  cout << "path: " << a << ", comparing to " << spath << " for run range " << rmin << " .. " << rmax << endl;
	// 	}
	if (!a.CompareTo(imap->first.c_str())) {
 	  fGoodHLT = true; 
	  if (fVerbose > 0) cout << "exact match: " << imap->first.c_str() << " HLT: " << a << " result: " << result << endl;
	}

 	if (a.Contains(spath.c_str()) && (rmin <= fRun) && (fRun <= rmax)) {
 	  fGoodHLT = true; 
	  if (fVerbose > 0) cout << "close match: " << imap->first.c_str() << " HLT: " << a 
				 << " result: " << result 
				 << " in run " << fRun 
				 << endl;
 	}
      }
    }      
  }


}



// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun,               "run/I");
  fTree->Branch("json",   &fJSON,              "json/O");
  fTree->Branch("evt",    &fEvt,               "evt/I");
  fTree->Branch("ls",     &fLS,                "ls/I");
  fTree->Branch("tm",     &fCandTM,            "tm/I");
  fTree->Branch("pr",     &fGenBpartial,       "pr/I"); 
  fTree->Branch("procid", &fProcessType,       "procid/I");
  fTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fTree->Branch("pvn",    &fPvN,               "pvn/I");
  fTree->Branch("cb",     &fCowboy,            "cb/O");

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
  fTree->Branch("cm",     &fCandM2,            "cm/D");
  fTree->Branch("cosa",   &fCandCosA,          "cosa/D");
  fTree->Branch("alpha",  &fCandA,             "alpha/D");
  fTree->Branch("iso",    &fCandIso,           "iso/D");
  fTree->Branch("isotrk", &fCandIsoTrk,        "isotrk/I");
  fTree->Branch("chi2",   &fCandChi2,          "chi2/D");
  fTree->Branch("dof",    &fCandDof,           "dof/D");
  fTree->Branch("prob",   &fCandProb,          "prob/D");
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
  fTree->Branch("m1chi2",  &fMu1Chi2,           "m1chi2/D");
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
  fTree->Branch("m2chi2",  &fMu2Chi2,           "m2chi2/D");

  fTree->Branch("mudist",  &fMuDist,            "mudist/D");
  fTree->Branch("mudeltar",&fMuDeltaR,          "mudeltar/D");

  fTree->Branch("g1pt",    &fMu1PtGen,          "g1pt/D");
  fTree->Branch("g2pt",    &fMu2PtGen,          "g2pt/D");
  fTree->Branch("g1eta",   &fMu1EtaGen,         "g1eta/D");
  fTree->Branch("g2eta",   &fMu2EtaGen,         "g2eta/D");

  fTree->Branch("t1pt",    &fMu1PtNrf,          "t1pt/D");
  fTree->Branch("t1eta",   &fMu1EtaNrf,         "t1eta/D");
  fTree->Branch("t2pt",    &fMu2PtNrf,          "t2pt/D");
  fTree->Branch("t2eta",   &fMu2EtaNrf,         "t2eta/D");

  fTree->Branch("hm1pt",  &fHltMu1Pt,  "hm1pt/D");    
  fTree->Branch("hm1eta", &fHltMu1Eta, "hm1eta/D");  
  fTree->Branch("hm1phi", &fHltMu1Phi, "hm1phi/D");  
  fTree->Branch("hm2pt",  &fHltMu2Pt,  "hm2pt/D");    
  fTree->Branch("hm2eta", &fHltMu2Eta, "hm2eta/D");  
  fTree->Branch("hm2phi", &fHltMu2Phi, "hm2phi/D");  


  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.); 
  h = 0; 

  TDirectory *pD;
  string name, dname; 
  int i(0); 
  for (map<string, int>::iterator imap = fRegion.begin(); imap != fRegion.end(); ++imap) {  
    i    = imap->second; 
    name = imap->first + "_";

    fpPvZ[i]       = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 50, -25., 25.);           
    fpPvN[i]       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 20, 0., 20.);           
    fpPvNtrk[i]    = bookDistribution(Form("%spvntrk", name.c_str()), "N_{trk}^{PV} ", "fGoodHLT", 20, 0., 200.);           
    fpTracksPt[i]  = bookDistribution(Form("%strackspt", name.c_str()), "p_{T} [GeV]", "fGoodTracksPt", 25, 0., 25.);
    fpTracksEta[i] = bookDistribution(Form("%strackseta", name.c_str()), "#eta_{T}", "fGoodTracksEta", 25, -2.5, 2.5);
    fpMuonsPt[i]   = bookDistribution(Form("%smuonspt", name.c_str()), "p_{T, #mu} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
    fpMuon1Pt[i]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
    fpMuon2Pt[i]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 25, 0., 25.); 
    fpMuonsEta[i]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", 20, -2.5, 2.5); 
    fpMuon1Eta[i]  = bookDistribution(Form("%smuon1eta", name.c_str()), "#eta_{#mu1}", "fGoodMuonsEta", 20, -2.5, 2.5); 
    fpMuon2Eta[i]  = bookDistribution(Form("%smuon2eta", name.c_str()), "#eta_{#mu2}", "fGoodMuonsEta", 20, -2.5, 2.5); 
    fpPt[i]        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 20, 0., 50.); 
    fpEta[i]       = bookDistribution(Form("%seta", name.c_str()), "#eta(B)", "fGoodEta", 20, -2.5, 2.5); 
    fpCosA[i]      = bookDistribution(Form("%scosa", name.c_str()), "cos(#alpha_{3D})", "fGoodAlpha", 30, 0.97, 1.); 
    fpAlpha[i]     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", "fGoodAlpha", 20, 0., 0.2); 
    fpIso[i]       = bookDistribution(Form("%siso", name.c_str()),  "isolation", "fGoodIso", 22, 0., 1.1); 
    fpIsoTrk[i]    = bookDistribution(Form("%sisotrk", name.c_str()),  "N_{trk}^{I}", "fGoodIso", 20, 0., 20.); 
    
    fpChi2[i]      = bookDistribution(Form("%schi2", name.c_str()),  "#chi^{2}", "fGoodChi2", 30, 0., 30.);              
    fpChi2Dof[i]   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", "fGoodChi2", 30, 0., 3.);       
    fpProb[i]      = bookDistribution(Form("%spchi2dof", name.c_str()),  "P(#chi^{2},dof)", "fGoodChi2", 25, 0., 1.);    
    fpFLS3d[i]     = bookDistribution(Form("%sfls3d", name.c_str()), "l_{3D}/#sigma(l_{3D})", "fGoodFLS", 25, 0., 100.);  
    fpFL3d[i]      = bookDistribution(Form("%sfl3d", name.c_str()),  "l_{3D} [cm]", "fGoodFLS", 25, 0., 5.);  
    fpFL3dE[i]     = bookDistribution(Form("%sfl3dE", name.c_str()), "#sigma(l_{3D}) [cm]", "fGoodFLS", 25, 0., 0.5);  
    fpFLSxy[i]     = bookDistribution(Form("%sflsxy", name.c_str()), "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.);  
    fpDocaTrk[i]   = bookDistribution(Form("%sdocatrk", name.c_str()), "d_{ca}^{min} [cm]", "fGoodDocaTrk", 35, 0., 0.14);   
    
    if (name == "A_") {
      for (int ipv = 0; ipv < NADPV; ++ipv) {
	pD = fHistDir->mkdir(Form("%sNpv%i", name.c_str(), ipv));
	pD->cd();
	dname = Form("%snpv%i_", name.c_str(), ipv);
	fpNpvPvN[ipv][i]     = bookDistribution(Form("%spvn", dname.c_str()), "N(PV) ", "fGoodHLT", 30, 0., 30.);           
	fpNpvChi2Dof[ipv][i] = bookDistribution(Form("%schi2dof", dname.c_str()), "#chi^{2}/dof", "fGoodChi2", 30, 0., 3.);       
	fpNpvProb[ipv][i]    = bookDistribution(Form("%spchi2dof", dname.c_str()), "P(#chi^{2},dof)", "fGoodChi2", 25, 0., 1.);    
	fpNpvFLS3d[ipv][i]   = bookDistribution(Form("%sfls3d", dname.c_str()), "l_{3d}/#sigma(l_{3d})", "fGoodFLS", 25, 0., 100.);  
	fpNpvFLSxy[ipv][i]   = bookDistribution(Form("%sflsxy", dname.c_str()), "l_{xy}/#sigma(l_{xy})", "fGoodFLS", 25, 0., 100.);  
	fpNpvDocaTrk[ipv][i] = bookDistribution(Form("%sdocatrk", dname.c_str()), "d_{ca}^{0} [cm]", "fGoodDocaTrk", 35, 0., 0.14);   
	fpNpvIso[ipv][i]     = bookDistribution(Form("%siso", dname.c_str()),  "I", "fGoodIso", 22, 0., 1.1); 
	fpNpvIsoTrk[ipv][i]  = bookDistribution(Form("%sisotrk", dname.c_str()),  "N_{trk}^{I}", "fGoodIso", 20, 0., 20.); 
	
      }
    } else {
      for (int ipv = 0; ipv < NADPV; ++ipv) {
	fpNpvPvN[ipv][i] = 0; 
      }
    }
  }
}

// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  // -- set up cut sequence for analysis
  basicCuts(); 
  moreBasicCuts(); 
  candidateCuts(); 
  moreCandidateCuts(); 

  fCutFile = fileName;

  if (dump) cout << "==> candAna: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;
  int ok(0);

  char  buffer[200];
  fHistDir->cd();
  if (dump) cout << "gDirectory: "; fHistDir->pwd();
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.c_str());
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

    if (!strcmp(CutName, "TRIGRANGE")) {
      char triggerlist[1000]; ok = 1; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      int r1(0), r2(0); 
      string hlt = splitTrigRange(tl, r1, r2); 
      HLTRANGE.insert(make_pair(hlt, make_pair(r1, r2))); 
      if (dump) {
	cout << "HLTRANGE:       " << hlt << " from " << r1 << " to " << r2 << endl; 
      }
      ibin = 5; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }


    if (!strcmp(CutName, "TRUTHCAND")) {
      TRUTHCAND = int(CutValue); ok = 1;
      if (dump) cout << "TRUTHCAND:           " << TRUTHCAND << endl;
      ibin = 7;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "IGNORETRIGGER")) {
      IGNORETRIGGER = int(CutValue); ok = 1;
      if (dump) cout << "IGNORETRIGGER      " << IGNORETRIGGER << endl;
      ibin = 8;
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

    if (!strcmp(CutName, "CANDALPHA")) {
      CANDALPHA = CutValue; ok = 1;
      if (dump) cout << "CANDALPHA:           " << CANDALPHA << endl;
      ibin = 15;
      hcuts->SetBinContent(ibin, CANDALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #alpha :: %5.4f", CutName, CANDALPHA));
    }

    if (!strcmp(CutName, "CANDFLS3D")) {
      CANDFLS3D = CutValue; ok = 1;
      if (dump) cout << "CANDFLS3D:           " << CANDFLS3D << endl;
      ibin = 16;
      hcuts->SetBinContent(ibin, CANDFLS3D);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{3d}/#sigma(l_{3d}) :: %3.1f", CutName, CANDFLS3D));
    }

    if (!strcmp(CutName, "CANDFLSXY")) {
      CANDFLSXY = CutValue; ok = 1;
      if (dump) cout << "CANDFLSXY:           " << CANDFLSXY << endl;
      ibin = 17;
      hcuts->SetBinContent(ibin, CANDFLSXY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{xy}/#sigma(l_{xy}) :: %3.1f", CutName, CANDFLSXY));
    }

    if (!strcmp(CutName, "CANDVTXCHI2")) {
      CANDVTXCHI2 = CutValue; ok = 1;
      if (dump) cout << "CANDVTXCHI2:           " << CANDVTXCHI2 << endl;
      ibin = 18;
      hcuts->SetBinContent(ibin, CANDVTXCHI2);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #chi^{2} :: %3.1f", CutName, CANDVTXCHI2));
    }

    if (!strcmp(CutName, "CANDISOLATION")) {
      CANDISOLATION = CutValue; ok = 1;
      if (dump) cout << "CANDISOLATION:           " << CANDISOLATION << endl;
      ibin = 19;
      hcuts->SetBinContent(ibin, CANDISOLATION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: I_{trk} :: %4.2f", CutName, CANDISOLATION));
    }

    if (!strcmp(CutName, "CANDDOCATRK")) {
      CANDDOCATRK = CutValue; ok = 1;
      if (dump) cout << "CANDDOCATRK:           " << CANDDOCATRK << endl;
      ibin = 20;
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

    //    if (!ok) cout << "==> candAna: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void candAna::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    exit(1);
  }
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
bool candAna::goodMuon(TAnaTrack *pT) {
  int result = pT->fMuID & MUIDMASK;

//   FIXME
//   if (HLTPath.end() != find(HLTPath.begin(), HLTPath.end(), "NOTRIGGER"))  {
//     return true; 
//   }

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


// ----------------------------------------------------------------------
string candAna::splitTrigRange(string tl, int &r1, int &r2) {

  string::size_type id1 = tl.find_first_of("("); 
  string::size_type id2 = tl.find_first_of(":"); 
  string::size_type id3 = tl.find_first_of(")"); 

  //cout << "tl: " << tl << endl;
  string hlt = tl.substr(0, id1);
  //cout << "hlt: " << hlt << endl;
  string a   = tl.substr(id1+1, id2-id1-1);
  r1 = atoi(a.c_str());
  //cout << "1st a: " << a << " -> r1 = " << r1 << endl;
  a  = tl.substr(id2+1, id3-id2-1); 
  r2 = atoi(a.c_str());
  //cout << "2nd a: " << a << " -> r2 = " << r2 << endl;

  return hlt; 

}


// ----------------------------------------------------------------------
double candAna::isoClassic(TAnaCand *pC) {
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
double candAna::isoClassicOnePv(TAnaCand *pC, double r, double ptmin) {
  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx; 
  int pvIdx = pC->fPvIdx;

  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    cIdx.push_back(pT->fIndex); 
    candPtScalar += pT->fPlab.Perp(); 
  }

  candPt = pC->fPlab.Perp(); 
  
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  continue;
    pT = fpEvt->getRecTrack(i); 
    if (pvIdx != pT->fPvIdx) {
      continue;
    }
    pt = pT->fPlab.Perp(); 
    if (pt < ptmin) {
      continue;
    }
    if (pT->fPlab.DeltaR(pC->fPlab) < r) {
      sumPt += pt; 
    } 
  }

  iso = candPt/(candPt + sumPt); 

  return iso; 
}

// ----------------------------------------------------------------------
double candAna::isoClassicWithDOCA(TAnaCand *pC, float docaCut, double r, double ptmin) {
  const double ptCut=ptmin, coneSize=r;
  const bool verbose=false;

  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx, pIdx; 
  int pvIdx = pC->fPvIdx;

  fCandItrk = 0; 
 
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
      pIdx.push_back(i); 
      ++fCandItrk;
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
      //      double docaE = pC->fNstTracks[i].second.second;

      if(doca > docaCut) continue; // check the doca cut
      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue; // skip tracks already included above
      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;
      
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
	++fCandItrk;
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
double candAna::isoWithDOCA(TAnaCand *pC, float docaCut, double r, double ptmin) {
  double ptCut=ptmin, coneSize=r;
  bool verbose=false;
  
  double iso(-1.), pt(0.), sumPt(0.), candPt(0.); 
  TAnaTrack *pT; 
  vector<int> cIdx, pIdx; 
  int pvIdx = pC->fPvIdx;
  
  fCandItrk = 0; 
  
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    cIdx.push_back(pT->fIndex); 
    if (verbose) {
      int tIdx = fpEvt->getRecTrack(pT->fIndex)->fPvIdx;
      cout << "Signal track " << pT->fIndex << endl;
      if (pvIdx != tIdx) {
	cout << "    pointing to PV " << tIdx << " instead of " << pvIdx << endl;
      }
    }
  }
  
  candPt = pC->fPlab.Perp(); 
  
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  continue;
    pT = fpEvt->getRecTrack(i); 
    if (pvIdx != pT->fPvIdx) {
      //      if (verbose) cout << " track " << i << " skipped because of PV index mismatch" << endl;
      continue;
    }
    pIdx.push_back(i); 
    pt = pT->fPlab.Perp(); 
    if (pt < ptCut) {
      //      if (verbose) cout << " track " << i << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (pT->fPlab.DeltaR(pC->fPlab) < coneSize) {
      ++fCandItrk;
      sumPt += pt; 
      if (verbose) 
	cout << "track " << i 
	     << " with pT = " << pT->fPlab.Perp()
	     << " eta = " << pT->fPlab.Eta()
	     << " dR = " << pT->fPlab.DeltaR(pC->fPlab)
	     << " pointing at PV " << pT->fPvIdx 
	     << endl;
    } 
    else {
      //      if (verbose) cout << " track " << i << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
    }
  }
  
  
  
  // Now consider ALL the DOCA tracks
  int nsize = pC->fNstTracks.size(); 
  if (nsize>0) {
    for(int i = 0; i<nsize;++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      //      double docaE = pC->fNstTracks[i].second.second;
      
      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue; // skip tracks already included above
      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;
      if(doca > docaCut) continue; // check doca cut
   
      pT = fpEvt->getRecTrack(trkId);

      pt = pT->fPlab.Perp();  // cut on track pt
      if (pt < ptCut) {
	//	if (verbose) cout << " close track " << i << " skipped because of pt = " << pt << endl;
	continue;
      }

      if (pT->fPlab.DeltaR(pC->fPlab) < 99.) {
	++fCandItrk;
	sumPt += pt; 
	if (verbose) cout << "close track " << i 
			  << " with pT = " << pT->fPlab.Perp()
			  << " eta = " << pT->fPlab.Eta()
			  << " pointing at PV " << pT->fPvIdx 
			  << endl;
      } 
      else {
	//	if (verbose) cout << " close track " << i << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
      }

    } // for loop over tracks
  } // end if 

  iso = candPt/(candPt + sumPt); 

  return iso; 
}

