#include "treeBmm.hh"

#include "TRandom.h"

#include "treeBmm.icc"


// Run with: ./bmm -c chains/bg-test -D root ; ./bmm -c chains/sg-test -D root ;

// ----------------------------------------------------------------------
void treeBmm::startAnalysis() {
  cout << "startAnalysis: ..." << endl;
}

// ----------------------------------------------------------------------
void treeBmm::eventProcessing() {

  if (fDebug & 1) cout << "==> Event: " << fEvent << endl;

  fER1->Fill(0.1); 

  fpHistFile->cd();
  ((TH1D*)fpHistFile->Get("runs"))->Fill(fpEvt->fRunNumber);

  // -- get candidate and cand./track properties
  initVariables();

  // -- Generator-level process discrimation
  processDiscrimination();

  // -- Muon efficiency / fake rate
  muonEfficiency();

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // -- Generator-level preselection (fGoodKinematics)
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  kinematicSelection(0);
  if (0 == fGoodKinematics) {
    if (fDebug & 1) cout << " --> no good kinematics ... " << endl;
  } else {
    if (fDebug & 1) cout << " --> good kinematics ... " << endl;
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // -- L1 trigger (fGoodL1)
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  L1Selection(100);
  if (0 == fGoodL1) {
    if (fDebug & 1) cout << " --> no good L1 trigger ... " << endl;
  } else {
    if (fDebug & 1) cout << " --> good L1 trigger ... " << endl;
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // -- High-level trigger (fGoodHLT)
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  HLTSelection(200);
  if (0 == fGoodHLT) {
    if (fDebug & 1) cout << " --> no good HLT trigger ... " << endl;
  } else {
    if (fDebug & 1) cout << " --> good HLT trigger ... " << endl;
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // -- Check indeces of cand. tracks and PV in event (fGoodEvent)
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  candidateSelection(300);
  if (0 == fGoodEvent) {
    if (fDebug & 1) cout << " --> no good event, RETURN ... " << "-> Cand "  << fGoodCand << ", PV " << fGoodPV << endl;
    return;                                                    
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RETURN !!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  } else {
    if (fDebug & 1) cout << " --> good event, PROCEED ... " << endl;
  }  

  // -- B-Candidate
  candidateProperties();

  // -- analysis efficiency
  fillAnalysisEff();

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  // -- Muon reco. tracks cuts: pT, eta & tip (fGoodTT)
  trackSelection(400);
  if (fGoodTT < 0) {
    if (fDebug & 1)  cout << " --> no good tracks ... " << fGoodTT << endl;
    fGoodTT = 0; 
  } else {
    if (fDebug & 1) cout << " --> good tracks ... " << endl;
  } 

  // -- Check muon ID of cand. track
  trackProperties();

  // -- offline HLT emulation
  fillHist();

  // -- residual plots
  ptResiduals();

  vertexResiduals();

  fTree->Fill();

}



// ----------------------------------------------------------------------
void treeBmm::initVariables() {

  if (fDebug & 2) { cout << "initVariables> Start" << endl; }

  fpHistFile->cd();

  // -- Find candidates for signal or norm. channel
  fBCI    = -1;
  fSTI[0] = fSTI[1] = fSTI[2] = -1;
  fPtL0   = fPtL1   = fPtK    = -99.;
  fDptL0  = fDptL1  = fDptK   = -99.;
  fEtaL0  = fEtaL1  = fEtaK   = -99.;
  fTipL0  = fTipL1  = fTipK   = -99.;
  fQL0    = fQL1    = fQK     = -99;
  fChi2   = fNdof   = fProb   = -99;
  fMass   = fMassJ  = -99.;
  fPt     = fP      = -99.;
  fMcL0   = fMcL1   = fMcK    = -99;
  fMuL0   = fMuL1   = fMuK    = -99.;
  fMomL0  = fMomL1  = fMomK   = -99;
  fGMoL0  = fGMoL1  = fGMoK   = -99;
  fTruthL0= fTruthL1= fTruthK = -99;
  fTruthB = fTruthJ =-99;
  fEta    = fTau    = fTxy    = -99.;
  fL3d    = fS3d    = fLxy    = fSxy = -99.;
  fRMM    = fDeta   = fDphi   = -99.; 
  fCosAngle = fCosAngle3 = -99.;

  fSig1   = fSig2   = fSig3   = -1;
  fgPtL0  = fgPtL1  = -99.;
  fgEtaL0 = fgEtaL1 = -99.;
  fgMuL0  = fgMuL1  = -99.;
  fgQL0   = fgQL1    = 99;
  fgRMM   = fgChi2   = fgS3d = -99.; 


  fnB     = fnJ     = 0;
  fgrB    = fgrJ     = 0;
  fgB     = fgJ     = 0;
  fgBmm   = fgJmm   = 0;
  fgMu    = frMu    = fnMu   = 0; 
  fnSigMu = frSigMu = 0;  
  fgK     = frK     = 0;
  fnTmMu  = fnTmK   = 0;
  frSigK  = 0;

  int bs_cand(-1), bplus_cand(-1);
  int sub_i(-1);

  for (int i = 0; i < NSUBSEL; i++ ) {

    sub_i = i + 1;

    if ( fNorm ) {
      
      fSubCand[i] = getNormCand(521, 443, fSel, sub_i);

      if ( sub_i == fSubSel ) {

	bplus_cand = fSubCand[i];
      }
      
    } else {
      
      fSubCand[i] = getSigCand(531, fSel, sub_i);

      if ( sub_i == fSubSel ) {

	bs_cand = fSubCand[i];
      }
    }
  }
  
  // -- Bs -> mu mu
  if ( bs_cand > -1 ) {
   
    fBCI = bs_cand;
    fpB   = fpEvt->getCand(bs_cand);
    
    // - B cand.
    fChi2 = fpB->fVtx.fChi2;
    fNdof = fpB->fVtx.fNdof;
    fProb = fpB->fVtx.fProb;
    
    fL3d  = fpB->fVtx.fD3d;
    fS3d  = fpB->fVtx.fD3dE;
    fLxy  = fpB->fVtx.fDxy;
    fSxy  = fpB->fVtx.fDxyE;
    
    fMass = fpB->fMass;

    fPt   = fpB->fPlab.Pt();
    fP    = fpB->fPlab.Mag();
    fEta  = fpB->fPlab.Eta();
    fTau  = fL3d*fMass/fP;
    fTxy  = fLxy*fMass/fPt;
    // fTxy  = fLxy*fCosAngle*fMass/fPt;

    // - muon cand.
    fSTI[0] = fpB->fSig1;
    fSTI[1] = fpB->fSig2;
    fpL1    = fpEvt->getRecTrack(fSTI[0]);
    fpL2    = fpEvt->getRecTrack(fSTI[1]);
    
    fPtL0 = fpL1->fPlab.Pt();
    fPtL1 = fpL2->fPlab.Pt();
    
    fDptL0 = fpL1->fPlab.Pt();
    fDptL1 = fpL2->fPlab.Pt();
    
    fEtaL0 = fpL1->fPlab.Eta();
    fEtaL1 = fpL2->fPlab.Eta();
    
    fTipL0 = fpL1->fTip;
    fTipL1 = fpL2->fTip;
    
    fQL0   = fpL1->fQ;
    fQL1   = fpL2->fQ;
    
  }

  // -- B+ -> mu mu K

  if ( bplus_cand > -1 ) {                                           
        
    fBCI = bplus_cand;

    // -- J/Psi -> mu mu                                            
    fpJpsi  = fpEvt->getCand(bplus_cand-1);

    // - J/Psi cand.   
    fPtJ   = fpJpsi->fPlab.Pt();
    fPJ    = fpJpsi->fPlab.Mag();
    fEtaJ  = fpJpsi->fPlab.Eta();

    fMassJ = fpJpsi->fMass;

    // - muon cand.
    fSTI[0] = fpJpsi->fSig1;
    fSTI[1] = fpJpsi->fSig2;
    fpL1    = fpEvt->getRecTrack(fSTI[0]);
    fpL2    = fpEvt->getRecTrack(fSTI[1]);
    
    fPtL0 = fpL1->fPlab.Pt();
    fPtL1 = fpL2->fPlab.Pt();
    
    fDptL0 = fpL1->fPlab.Pt();
    fDptL1 = fpL2->fPlab.Pt();
    
    fEtaL0 = fpL1->fPlab.Eta();
    fEtaL1 = fpL2->fPlab.Eta();
    
    fTipL0 = fpL1->fTip;
    fTipL1 = fpL2->fTip;
    
    fQL0   = fpL1->fQ;
    fQL1   = fpL2->fQ;      
    
    // -- B+ -> J/Psi K  
    fpB = fpEvt->getCand(bplus_cand);
    
    // -- B cand
    fChi2 = fpB->fVtx.fChi2;
    fNdof = fpB->fVtx.fNdof;
    fProb = fpB->fVtx.fProb;
    
    fL3d  = fpB->fVtx.fD3d;
    fS3d  = fpB->fVtx.fD3dE;
    fLxy  = fpB->fVtx.fDxy;
    fSxy  = fpB->fVtx.fDxyE;
    
    fMass = fpB->fMass;

    fPt   = fpB->fPlab.Pt();
    fP    = fpB->fPlab.Mag();
    fEta  = fpB->fPlab.Eta();
    fTau  = fL3d*fMass/fP;
    fTxy  = fLxy*fMass/fPt;
  
    if ( (fpJpsi->fSig1 != fpB->fSig1) || (fpB->fSig2 > fpEvt->nRecTracks()) ) {
 
      if (fpJpsi->fSig1 != fpB->fSig1) {
	cout << "==> Event: " << fEvent 
	     << " *** PROBLEM!!! J/Psi candidate points to different muon track than B+ candidate!!! ***" << endl;
	cout << " m1 = " << fpJpsi->fSig1 << " and m1 = " << fpB->fSig1 << endl;
	cout << " m2 = " << fpJpsi->fSig2 << " and m2 = " << fpB->fSig2 << endl;
      } 
      
      if (fpB->fSig2 > fpEvt->nRecTracks()) {
	cout << "==> Event: " << fEvent 
	     <<" *** PROBLEM!!! Track index of kaon is outside reco. tracks!!! ***" << endl;
	cout << " Kaon index = " << fpB->fSig2 << " and n_tracks = " << fpEvt->nRecTracks() << endl;
      }
      
      fSTI[0] = fSTI[1] = fSTI[2] = -1;
      fPtL0  = fPtL1  = fPtK  = -99.;
      fEtaL0 = fEtaL1 = fEtaK = -99.;
      fTipL0 = fTipL1 = fTipK = -99.;
      fQL0   = fQL1   = fQK   = -99.;

      return;
      
    } else {
      

      fSTI[2] = fpB->fSig2;
      fpK   = fpEvt->getRecTrack(fSTI[2]);
      fPtK  = fpK->fPlab.Pt();
      fDptK = fpK->fPlab.Pt();
      fEtaK = fpK->fPlab.Eta();
      fTipK = fpK->fTip;
      fQK   = fpK->fQ;
    }
  }
  
  // -- debug
  if (fDebug & 2) {

    if ( bs_cand > -1 ) {

      cout << "==> Event: " << fEvent << endl
	   << "  Choosing signal cand " << bs_cand << " which is of type " << fpB->fType
	   << " with pt,phi,eta = "
	   << fpB->fPlab.Pt()  << ", "
	   << fpB->fPlab.Phi() << ", "
	   << fpB->fPlab.Eta()
	   << " and has signal tracks at " << endl
	   << fSTI[0]
	   << " with pt,phi,eta = "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Pt()  << ", "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Phi() << ", "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Eta()
	   << " and " << endl
	   << fSTI[1]
	   << " with pt,phi,eta = "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Pt()  << ", "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Phi() << ", "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Eta()
	   << endl;
    }
    
    if ( bplus_cand > -1 ) {
	
      cout << "==> Event: " << fEvent << endl
	   << "  Choosing norm. cand " << bplus_cand << " which is of type " << fpJpsi->fType
	   << " with pt,phi,eta = "
	   << fpJpsi->fPlab.Pt()  << ", "
	   << fpJpsi->fPlab.Phi() << ", "
	   << fpJpsi->fPlab.Eta()
	   << " and has signal tracks at " << endl
	   << fSTI[0]
	   << " with pt,phi,eta = "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Pt()  << ", "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Phi() << ", "
	   << fpEvt->getRecTrack(fSTI[0])->fPlab.Eta()
	   << " and " << endl
	   << fSTI[1]
	   << " with pt,phi,eta = "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Pt()  << ", "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Phi() << ", "
	   << fpEvt->getRecTrack(fSTI[1])->fPlab.Eta()
	   << endl;
	
      cout << " and candidate " << bplus_cand+1 << " which is of type " << fpB->fType
	   << " with pt,phi,eta = "
	   << fpB->fPlab.Pt()  << ", "
	   << fpB->fPlab.Phi() << ", "
	   << fpB->fPlab.Eta() << endl
	   << "  and has signal track at "
	   << fSTI[2]
	   << " with pt,phi,eta = "
	   << fpEvt->getRecTrack(fSTI[2])->fPlab.Pt()  << ", "
	   << fpEvt->getRecTrack(fSTI[2])->fPlab.Phi() << ", "
	   << fpEvt->getRecTrack(fSTI[2])->fPlab.Eta() 
	   << endl;
    }
  }

  // -- Indeces to B+ candidates of all kaon candidates
  TAnaCand *pCand;
  fKTI.clear();
  fKCI.clear();

  for (int i = 0; i < fpEvt->nCands(); ++i) {
    
    if (fpEvt->getCand(i)->fType == 0 ) {  

      pCand = fpEvt->getCand(i);
      fKCI.push_back(i);
      fKTI.push_back(pCand->fSig2);
    }

    if (int(fpEvt->getCand(i)->fType/10) != 0 && int(fpEvt->getCand(i)->fType/10) != 443 && 
	int(fpEvt->getCand(i)->fType/10) != 531 && int(fpEvt->getCand(i)->fType/10) != 521  ) { 

      pCand = fpEvt->getCand(i); 
      cout << "Candidate " << i << " has unknown type: " << pCand->fType << endl;
    }
  }

  // ======================================================================
  // -- Systematics: Change some variables 
  // ======================================================================
  // -- muon ID
  if (0 && gRandom->Rndm() < 0.01) {
    fpEvt->fL1Decision = 1; 
  }

  // -- L1
  if (0 && fpEvt->fL1Decision == 0 && gRandom->Rndm() < 0.05) {
    fpEvt->fL1Decision = 1; 
  }

  // -- Tracking efficiency 2*1%
  if (0) {
    if (gRandom->Rndm() < 0.005) {
      fSTI[0] = -1; 
    }

    if (gRandom->Rndm() < 0.005) {
      fSTI[1] = -1; 
    }
  }

  // -- Tracking resolution
  if (0) {

    double k0 = 1./fPtL0;
    double k1 = 1./fPtL1;

    double smear, sigma; 
    sigma = 0.0005;  // momentum scale
    //    sigma = 0.0004;  // misalignment
    //    sigma = 0.0003;  // B field

    smear = gRandom->Gaus(0., sigma);
    k0 += smear;
    smear = gRandom->Gaus(0., sigma);
    k1 += smear;

    fPtL0 = 1./k0;
    fPtL1 = 1./k1;
  } 

  if (fDebug & 2) { cout << "initVariables> End" << endl; }   
}

// ---------------------------------------------------------------------- 
int treeBmm::getSigCand(int cand, int sel, int crit) {

  // -- criteria: 1 = leading muon + second highest pT muon
  //              2 = closest muon
  //              3 = leading muon + muon closest to it 
  //              4 = best chi2
 
  int sigCand(-1), sigMu1(-1), sigMu2(-1);

  int id(-1), bmmsel(-1);
  int i1(-1), i2(-1), itmp(-1);
  int nB(0), gB(0);

  // -- signal cand. criteria
  double pt1(0.), pt2(0.), pt_max(0.), pt_max2(0.);

  double dphi(0.), deta(0.);
  double rmm(99999.),  rmm_min(99999.);
  double chi2(99999.), chi2_min(99999.);

  // -- signal cand. quality 
  double s3d(100.);
  
  TAnaCand   *B;
  TAnaTrack  *L1, *L2;

  for (int i = 0; i < fpEvt->nCands(); ++i) {
    
    id     = int(fpEvt->getCand(i)->fType/10);
    bmmsel = fpEvt->getCand(i)->fType - 10*id;

    if ( id == cand && bmmsel == sel ) {
 
      nB++;

      B  = fpEvt->getCand(i);
      
      i1   = B->fSig1;
      i2   = B->fSig2;

      chi2 = B->fVtx.fChi2;
      s3d  = B->fVtx.fD3dE;

      if ( i1 >= 0 && i2 >= 0 ) { 

 	L1    = fpEvt->getRecTrack(i1);
	L2    = fpEvt->getRecTrack(i2);

	if ( (L1->fPlab.Pt() < fMinPt) ||
	     (L2->fPlab.Pt() < fMinPt) ||
	     (TMath::Abs(L1->fPlab.Eta()) > fMaxEta) || 
	     (TMath::Abs(L2->fPlab.Eta()) > fMaxEta) || 
	     (L1->fQ*L2->fQ > -0.5) ||
	     (s3d > 1.0)
	     
	     ) {
	  
	  continue;
	}

	// -- Re-arrange muons in descending pT-order --------------
	// ---------------------------------------------------------
	// -- :::: This should be obsolete, as the muons are already 
	// -- :::: filled in descending pT order in the ED Analyser!
	// ---------------------------------------------------------
	if ( L2->fPlab.Pt() > L1->fPlab.Pt() ) {
	  
	  itmp = i1;
	  i1 = i2;
	  i2 = itmp;

	  L1    = fpEvt->getRecTrack(i1);
	  L2    = fpEvt->getRecTrack(i2);
	}

	// ----------------------------------------------------------


	// -- two closest muons
	if (  crit == 1 ) {
	    
	  dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	  deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	  rmm  = TMath::Sqrt(dphi*dphi + deta*deta);

	  if ( rmm < rmm_min ) {
	    
	    sigMu1  = i1;
	    sigMu2  = i2;
	    sigCand = i;
	    
	    rmm_min = rmm;
	  }
	}

	// -- leading muon
	if (  crit == 2 ) {

	  pt1 = L1->fPlab.Pt();
	  pt2 = L2->fPlab.Pt();

	  if (pt1 > pt_max) {

	    pt_max  = pt1; 
	    pt_max2 = pt2; 

	    sigMu1  = i1;
	    sigMu2  = i2;
	    sigCand = i;

	    
	  } else if ( i1 == sigMu1 ) {

	    if ( pt2 > pt_max2 ) {

	      sigMu2  = i2;
	      sigCand = i;

	      pt_max2 = pt2;
	    }
	  }
	}

	// -- leading muon & closest muon
	if ( crit == 3 ) {

	  pt1 = L1->fPlab.Pt();

	  if ( pt1 > pt_max ) {
	    
	    dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	    deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	    rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	    
	    pt_max = pt1; 
	    rmm_min = rmm;

	    sigMu1  = i1;
	    sigMu2  = i2;
	    sigCand = i;
	    
	  } else if ( i1 == sigMu1 ) {
	    
	    dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	    deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	    rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	    
	    if ( rmm < rmm_min ) {

	      sigMu2  = i2;
	      sigCand = i;
	      
	      rmm_min = rmm;
	    }
	  }
	}

	// -- best chi2
	if (  crit == 4 ) {
	    
	  if ( chi2 < chi2_min ) {

	    sigMu1  = i1;
	    sigMu2  = i2;
	    sigCand = i;
	    
	    chi2_min = chi2;
	  }
	}
      }
    }

    // -- Get the truth-matched signal candidate (takes last one) 
    if ( id == cand && bmmsel == 1 ) {
      
      B  = fpEvt->getCand(i);
      
      fSig1 = B->fSig1;
      fSig2 = B->fSig2;

      fSigB = i;

      if (fSig1 >= 0 && fSig2 >= 0) {
      
	gB++;
      }
    }
  }

  fgrB = gB; 
  fnB = nB; 
  
  return sigCand;
}

// ---------------------------------------------------------------------- 
int treeBmm::genPartType(int index) {  

  if ( index < 0 ) {

    return -1;
  }

  int type(-1);
  int nD(0), subidx(-1);

  TGenCand *pG, *pD1, *pD2, *pGD1, *pGD2;

  pG = fpEvt->getGenCand(index);
  nD = getNrOfDaughters(index);

  if ( pG->fDau1 > 1 && pG->fDau2 > 1 ) { 
      
    pD1   = fpEvt->getGenCand(pG->fDau1);
    pD2   = fpEvt->getGenCand(pG->fDau2);
    
    if ( TMath::Abs(pD1->fID) == 13 && TMath::Abs(pD2->fID) == 13) {
	

      // -- generated J/Psi+ -> mu+ mu-
      if (TMath::Abs(pG->fID) == 443) {

	((TH1D*)gDirectory->Get("m702"))->Fill(nD);

	if (nD == 2) {
	  
	  type = 443;
	}
      }

      // -- generated Bs0 -> mu+ mu-
      if (TMath::Abs(pG->fID) == 531) {

	((TH1D*)gDirectory->Get("m700"))->Fill(nD);
	
	if (nD == 2) {

	  type = 531;
	}
      }
    }

    // -- generated B+ -> J/Psi K+
    // - D1 = J/Psi
    if ( TMath::Abs(pD1->fID) == 443 && TMath::Abs(pD2->fID) == 321 ) {	
      
      if ( pD1->fDau1 > 1 && pD1->fDau2 > 1 ) { 
	  
	subidx = pD1->fNumber;
	pGD1   = fpEvt->getGenCand(pD1->fDau1);
	pGD2   = fpEvt->getGenCand(pD1->fDau2);
	  
	if ( TMath::Abs(pGD1->fID) == 13 && TMath::Abs(pGD2->fID) == 13 &&
	     getNrOfDaughters(subidx) == 2 ) {
	    
	  ((TH1D*)gDirectory->Get("m701"))->Fill(nD);

	  if (nD == 2) {

	    type = 521;
	  }
	}
      }
    }

    // -- D2 = J/Psi
    if ( TMath::Abs(pD1->fID) == 321 && TMath::Abs(pD2->fID) == 443 ) {

      if ( pD2->fDau1 > 1 && pD2->fDau2 > 1 ) { 
	
	subidx = pD2->fNumber;
	pGD1   = fpEvt->getGenCand(pD2->fDau1);
	pGD2   = fpEvt->getGenCand(pD2->fDau2);
	
	if ( TMath::Abs(pGD1->fID) == 13 && TMath::Abs(pGD2->fID) == 13 &&
	     getNrOfDaughters(subidx) == 2 ) {
	    
	  ((TH1D*)gDirectory->Get("m701"))->Fill(nD);
	  
	  if (nD == 2) {

	    type = 521;
	  }
	}
      }
    }
  }
  
  return type;
}
      

// ----------------------------------------------------------------------
int treeBmm::getNrOfDaughters(int index) { 

  // -- #muons / per gen. B
  TGenCand *pG;

  int nDau(0);

  for (int j = 2; j < fpEvt->nGenCands(); ++j) {
    
    pG = fpEvt->getGenCand(j); 
    
    if ( TMath::Abs(pG->fMom1) == index ) {
      
      nDau++;
    }
  }

  return nDau;

}

// ----------------------------------------------------------------------
int treeBmm::getNormCand(int cand, int subcand, int sel, int crit) {

   // -- criteria: 1 = leading muon + second highest pT muon
  //              2 = closest muon
  //              3 = leading muon + muon closest to it 
  //              4 = best chi2
 
  int normCand(-1), normMu1(-1), normMu2(-1);

  int id(-1), ids(-1), bmmsel(-1);
  int i1(-1), i2(-1), itmp(-1);
  int nB(0), nJ(0), gB(0), gJ(0);

  // -- signal cand. criteria
  double pt1(0.), pt2(0.), pt_max(0.), pt_max2(0.);

  double dphi(0.), deta(0.);
  double rmm(99999.),  rmm_min(99999.);
  double chi2(99999.), chi2_min(99999.);

  // -- signal cand. quality 
  double s3d(100.);

  TAnaCand   *B, *J;
  TAnaTrack  *L1, *L2;

  for (int i = 0; i < fpEvt->nCands(); ++i) {
    
    id     = int(fpEvt->getCand(i)->fType/10);
    bmmsel = fpEvt->getCand(i)->fType - 10*id;
 
    if ( id == cand && bmmsel == sel ) {
 
      nB++;

      ids    = int(fpEvt->getCand(i - 1)->fType/10);
      bmmsel = fpEvt->getCand(i - 1)->fType - 10*ids;

      if ( ids == subcand && bmmsel == sel ) {

	fnJ++;
	
	J  = fpEvt->getCand(i - 1);
	
	i1 = J->fSig1;
	i2 = J->fSig2;
	
	chi2 = J->fVtx.fChi2;
	s3d  = J->fVtx.fD3dE;

	if ( i1 >= 0 && i2 >= 0 ) { 
	  
   	  L1    = fpEvt->getRecTrack(i1);
	  L2    = fpEvt->getRecTrack(i2);
	  
	  if ( (L1->fPlab.Pt() < fMinPt) ||
	       (L2->fPlab.Pt() < fMinPt) ||
	       (TMath::Abs(L1->fPlab.Eta()) > fMaxEta) || 
	       (TMath::Abs(L2->fPlab.Eta()) > fMaxEta) || 
	       (L1->fQ*L2->fQ > -0.5) ||
	       (s3d > 1.0)
	       
	       ) {
	    
	    continue;
	  }
	 	 
	  // -- Re-arrange muons in descending pT-order --------------
	  // ---------------------------------------------------------
	  // -- :::: This should be obsolete, as the muons are already 
	  // -- :::: filled in descending pT order in the ED Analyser!
	  // ---------------------------------------------------------
	  if ( L2->fPlab.Pt() > L1->fPlab.Pt() ) {
	    
	    itmp = i1;
	    i1 = i2;
	    i2 = itmp;
	    
	    L1    = fpEvt->getRecTrack(i1);
	    L2    = fpEvt->getRecTrack(i2);
	  }
	  
	  // ----------------------------------------------------------
	  
	  
	  // -- two closest muons
	  if (  crit == 1 ) {
	 
	    dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	    deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	    rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	    
	    if ( rmm < rmm_min ) {
	      
	      normMu1  = i1;
	      normMu2  = i2;
	      normCand = i;
	      
	      rmm_min = rmm;
	    }
	  }
	  
	  // -- leading muon
	  if (  crit == 2 ) {
	    
	    pt1 = L1->fPlab.Pt();
	    pt2 = L2->fPlab.Pt();
	    
	    if (pt1 > pt_max) {
	      
	      pt_max  = pt1; 
	      pt_max2 = pt2; 
	      
	      normMu1  = i1;
	      normMu2  = i2;
	      normCand = i;
	      
	      
	    } else if ( i1 == normMu1 ) {
	      
	      if ( pt2 > pt_max2 ) {
		
		normMu2  = i2;
		normCand = i;
		
		pt_max2 = pt2;
	      }
	    }
	  }
	  
	  // -- leading muon & closest muon
	  if ( crit == 3 ) {
	    
	    pt1 = L1->fPlab.Pt();
	    
	    if ( pt1 > pt_max ) {
	      
	      dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	      deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	      rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	      
	      pt_max = pt1; 
	      rmm_min = rmm;
	      
	      normMu1  = i1;
	      normMu2  = i2;
	      normCand = i;
	      
	    } else if ( i1 == normMu1 ) {
	      
	      dphi = L1->fPlab.DeltaPhi(L2->fPlab);
	      deta = L1->fPlab.Eta() - L2->fPlab.Eta();
	      rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	      
	      if ( rmm < rmm_min ) {
		
		normMu2  = i2;
		normCand = i;
		
		rmm_min = rmm;
	      }
	    }
	  }
	  
	  // -- best chi2
	  if (  crit == 4 ) {
	    
	    if ( chi2 < chi2_min ) {
	      
	      normMu1  = i1;
	      normMu2  = i2;
	      normCand = i;
	      
	      chi2_min = chi2;
	    }
	  }
	}
      }
    }

    // -- Get the truth-matched signal candidate (takes last one)
    if ( ids == subcand && bmmsel == 1 ) {
      
      gJ++;
      
      J  = fpEvt->getCand(i);
      
      fSig1 = J->fSig1;
      fSig2 = J->fSig2;

      fSigJ = i;
   
    }

    if ( id == cand && bmmsel == 1 ) {
      
      gB++;
      
      B  = fpEvt->getCand(i);
      
      fSig3 = B->fSig2;

      fSigB = i;
    } 
  }

  fgrB = gB; 
  fgrJ = gJ;
  fnB = nB; 
  fnJ = nJ;
  
  return normCand;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- This is on the GENERATOR level
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::kinematicSelection(int o) {
  // offset: 0

  if (fDebug & 2) { cout << "kinematicSelection> Start" << endl; }
  
  fGoodKinematics = 1;

  TGenCand *pG;

  double pT(0.);
  double eta(0.);

  int cnt(0); 

  int m0(-1), m1(-1); 
  int pid0(-1), pid1(-1), pidTmp(-1);
  double pT0(0.), pT1(0.);
 

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %% The first two particles in gen. block are the protons and %%
  // %% have no valid pT / eta entry. Therefore start at two! To  %%
  // %% extract gen. level info only use tracks with gIndex > 1.  %%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for (int i = 2; i < fpEvt->nGenCands(); ++i) {

    pG = fpEvt->getGenCand(i);

    if ( pG->fP.Pt() < fMinPt || 
	    TMath::Abs(pG->fP.Eta()) > fMaxEta ) { 
      
      continue; 
    }

    if ( TMath::Abs(pG->fID) == fD0 || TMath::Abs(pG->fID) == fD1 || !strcmp("qcd", fChannel) ) {
      
      ++cnt;
      
      pT  = pG->fP.Pt();
      eta = pG->fP.Eta();
      
      if (pT > pT0) {
	
	if (m0 > 1) {
	  
	  m1 = m0;
	  pid1 = pid0;
	  pT1 = pT0;
	}
	
	pT0 = pT;
	m0 = i;
	pid0 = TMath::Abs(pG->fID);

      }  else  if (pT > pT1) {
	
	pT1 = pT;
	m1 = i;
	pid1 = TMath::Abs(pG->fID);
      }
    }
    
    // -- Rare BG with different particles
    
    if ( fD0 != fD1  && pid0 == pid1 && strcmp("qcd", fChannel) ) {  
      
      cnt = 1;
      
      if ( pid1 == fD0 ) { 
	
	pidTmp = fD1; 
      }
      
      if ( pid1 == fD1 ) {
	
	pidTmp = fD0; 
      }
      
      m1 = -1;
      pid1 = -1;
      pT1 = 0.;
      
      for (int i = 2; i < fpEvt->nGenCands(); ++i) {
	
	pG = fpEvt->getGenCand(i);
	
	if ( pG->fP.Pt() < fMinPt || 
	     TMath::Abs(pG->fP.Eta()) > fMaxEta ) { 
	  
	  continue; 
	}	

	if ( TMath::Abs(pG->fID) == pidTmp ) {
	  
	  ++cnt;
	  
	  pT  = pG->fP.Pt();
	  eta = pG->fP.Eta();
	  
	  if (pT > pT1) {
	    
	    pT1 = pT;
	    m1 = i;
	    pid1 = TMath::Abs(pG->fID);
	    
	  }
	}
      }
    }
  }


  // -- Eff. histogram
  fER1->Fill(o+0.1); 
  if (cnt >= 2 ) {
    fER1->Fill(o+1.1); 
    
  } else {
    fGoodKinematics = 0;
    fER1->Fill(o+2.1);
  }

  
  // -- Gen. muon statistics ----------------------------------------

  int genType(-1);

  for (int i = 2; i < fpEvt->nGenCands(); ++i) {
    

    pG = fpEvt->getGenCand(i);
    genType = genPartType(i);

    if ( genType == 443 ) {
      fgJmm++;
    }

    if ( !fNorm && genType == 531 ) {
      fgBmm++;
    }

    if (  fNorm && genType == 521 ) {
      fgBmm++;
    }

    // -- #total gen. Bs0
    if ( !fNorm && TMath::Abs(pG->fID) == 531 ) {
      
      fgB++;
    }

    // -- #total gen. B+/-
    if ( fNorm && TMath::Abs(pG->fID) == 521 ) {
      
      fgB++;
    }

    // -- #total gen. J/Psi's
    if ( TMath::Abs(pG->fID) == 443 ) {
      
      fgJ++;
    }

    // -- #total gen. muons
    if ( TMath::Abs(pG->fID) == 13 ) {
      
      fgMu++;
    }

    // -- #total gen. kaons
    if ( TMath::Abs(pG->fID) ==321 ) {
      
      fgK++;
    }
  }

  if (fDebug & 2) { cout << "kinematicSelection> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- L1 trigger decision
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::L1Selection(int o) {
  // offset: 100

  if (fDebug & 2) { cout << "L1Selection> Start" << endl; }

  fGoodL1 = 1;
  
  int L1 = fpEvt->fL1Decision;

  // -- Rare backgrounds
  if ( SETL1 > 0) { 
    
    L1 = 1; 
  }

  // -- Eff. histogram
  fER1->Fill(o+0.1);
  if (!L1) {
    fGoodL1 = 0; 
    fER1->Fill(o+12.1);
  } else {
    fER1->Fill(o+11.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodL1) {
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }

  if (fDebug & 2) { cout << "L1Selection> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- HLT trigger decision
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::HLTSelection(int o) {
  // offset: 200

  if (fDebug & 2) { cout << "HLTSelection> Start" << endl; }

  fGoodHLT = 1;

  int HLT = fpEvt->fHLTDecision;

  // -- Rare backgrounds
  if ( SETHLT > 0) { 

    HLT = 1; 
  }

  // -- Eff. histograms
  fER1->Fill(o+0.1);
  if (!HLT) {
    fGoodHLT = 0;
    fER1->Fill(o+12.1);
  } else {
    fER1->Fill(o+11.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodHLT) { 
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1)) {
    if (0 == fGoodHLT) { 
      fER1->Fill(o+4.1);
    } else {
      fER1->Fill(o+3.1);
    }
  }

  if (fDebug & 2) { cout << "HLTSelection> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- Reconstructed candidate / Primary Vertex found?
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::candidateSelection(int o) {
  // offset: 300

  if (fDebug & 2) { cout << "candidateSelection> Start" << endl; }

  fGoodEvent = 1;
  fGoodCand = 1;
  fGoodPV = 1;

  if (fBCI < 0 ) {
    if (fDebug & 2) cout << "cand. selection> no candidate found " << endl;
    fGoodCand = -1;
    fER1->Fill(o+22.1);
  } else {
    fER1->Fill(o+21.1);

    // -- Signal tracks
    if (fSTI[0] < 0) {
      if (fDebug & 2) cout << "cand. selection> no index for signal track 0 " << endl;
      fGoodCand = -2;
      fER1->Fill(o+24.1);
    } else {
      fER1->Fill(o+23.1);
      
    }
    
    if (fSTI[1] < 0) {
      if (fDebug & 2)    cout << "cand. selection> no index for signal track 1 "  << endl;
      fGoodCand = -3;
      fER1->Fill(o+26.1);
    } else {
      fER1->Fill(o+25.1);
    }
    
    if (fD2 > 0 && fSTI[2] < 0) {
      if (fDebug & 2)    cout << "cand. selection> no index for signal track 2 "  << endl;
      fGoodCand = -4;
      fER1->Fill(o+28.1);
    } else {
      fER1->Fill(o+27.1);
    }
  }
  
  if ( !fNorm && fgBmm > 1 ) {
    if (fDebug & 2)    cout << "cand. selection> more than one Bs0 -> mu mu generated "  << endl;
    fGoodCand = fgBmm;
  }

  if (fChainFileName.Contains("cbg-00") && fgBmm > 0 ) {
    if (fDebug & 2)    cout << "cand. selection> this background events contain a signal !!!"  << endl;
    fGoodCand = -99;
  }

  if (fChainFileName.Contains("csg-004n") && fgBmm < 1 ) {
    if (fDebug & 2)    cout << "cand. selection> this signal event does not contain a signal !!!"  << endl;
    fGoodCand = -66;
  }

  // -- Primary vertex
  if (-99 == fpEvt->fPrimaryVertex2.fType) {
    if (fDebug & 2)    cout << "cand. selection> no primary vertex "  << endl;
    fGoodPV = 0;
    fER1->Fill(o+32.1);
  } else {
    fER1->Fill(o+31.1);
  }


  // -- Eff. histogram
  fER1->Fill(o+0.1);
  if (fGoodCand != 1 || 0 == fGoodPV) {
    fGoodEvent = 0;
    fER1->Fill(o+12.1);
  } else {
    fER1->Fill(o+11.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodEvent) { 
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1)) {
    if (0 == fGoodEvent) { 
      fER1->Fill(o+4.1);
    } else {
      fER1->Fill(o+3.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1) &&  (1 == fGoodHLT)) {
    if (0 == fGoodEvent) {
      fER1->Fill(o + 6.1);
    } else {
      fER1->Fill(o + 5.1);
    }
  }

  // -- This methods have to be called before RETURN statement ----
  fillEventStats();

  fillTriggerEff();

  int bs_i(-1), sub_i(-1);
  for (int i = 0; i < NSUBSEL; i++) { 

    sub_i = i + 1;
    bs_i  = fSubCand[i];

    if ( fOffCand > 0 ) { 

      offlineEff(bs_i, sub_i);
    }
  }

  if (fDebug & 2) { cout << "candidateSelection> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- Candidate properties
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::candidateProperties() {  

  if (fDebug & 2) { cout << "candidateProperties> Start" << endl; }

  TH1D *h;

  if (fChainFileName.Contains("sg-001")) {

    fMass *= 5.374/5.29;  //  <----- ?????
  }

  // -- RMM:        Opening cone of signal muons  
  fDeta = fpL1->fPlab.DeltaPhi(fpL2->fPlab);
  fDphi = fpL1->fPlab.Eta() - fpL2->fPlab.Eta();
  fRMM  = TMath::Sqrt(fDphi*fDphi + fDeta*fDeta);

  // -- COSALSPHA:   3D version: l is the vector from the PV to the SV
  TVector3 l  = fpB->fVtx.fPoint - fpEvt->fPrimaryVertex2.fPoint;
  fCosAngle3  = TMath::Cos(fpB->fPlab.Angle(l));

  // -- COSALSPHA:   2D version in transverse plane
  TVector2 tPV = fpEvt->fPrimaryVertex2.fPoint.XYvector();
  TVector2 tSV = fpB->fVtx.fPoint.XYvector();
  TVector2 tB  = fpB->fPlab.XYvector();
  TVector2 tl  = tSV - tPV;
  fCosAngle = TMath::Cos(tB.DeltaPhi(tl));
  
  
  // -- ISOLATION:   Compute isolation with tracks
  TVector3 ptv;
  int absQ(0);
  double df2, de2;
  double dr(0.), drmin(99.);
  double sum(0.), sum05(0.), sum06(0.), sum08(0.), sum10(0.), sum11(0.), sum12(0.), sum14(0.);

  double cone = 0.5*fRMM + 0.4;
  
  fIsoVeto = fIsoVeto05 = fIsoVeto06 = fIsoVeto08 = fIsoVeto10 = fIsoVeto11 = fIsoVeto12 = fIsoVeto14 = 0;

  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {

    ptv    = fpEvt->getRecTrack(it)->fPlab;
    absQ   = TMath::Abs(fpEvt->getRecTrack(it)->fQ);

    if ( fNorm ) {

      df2 = fpJpsi->fPlab.DeltaPhi(ptv);
      de2 = (ptv.Eta() - fpJpsi->fPlab.Eta());
      dr  = TMath::Sqrt(df2*df2 + de2*de2);

      if (it == fSTI[0]) { fRMJ1 = dr; continue; }
      if (it == fSTI[1]) { fRMJ2 = dr; continue; }
      if (it == fSTI[2]) { fRKJ  = dr; continue; } 

      ((TH1D*)gDirectory->Get("i300"))->Fill(df2);
      ((TH1D*)gDirectory->Get("i301"))->Fill(ptv.Pt());

      if ((ptv.Pt() > ISOPTLO) && (dr < cone) && (absQ > 0) ) {
	((TH1D*)gDirectory->Get("i302"))->Fill(ptv.Pt());
      }
    
    } else {

      if (it == fSTI[0]) { continue; }
      if (it == fSTI[1]) { continue; }
    }
    
    df2 = fpB->fPlab.DeltaPhi(ptv);
    de2 = (ptv.Eta() - fpB->fPlab.Eta());
    dr  = TMath::Sqrt(df2*df2 + de2*de2);

    ((TH1D*)gDirectory->Get("i200"))->Fill(df2);
    ((TH1D*)gDirectory->Get("i201"))->Fill(ptv.Pt());

    if ((dr < cone) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) {
      sum += ptv.Pt();
      ((TH1D*)gDirectory->Get("i202"))->Fill(ptv.Pt());
    }


    // -- Systematics
    double lostTrack = (gRandom->Rndm() < 0.01? 0. : 1.);
    lostTrack = 1.0; // this is equivalent to NO smearing
    
    if ((dr < 0.5) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum05 += lostTrack*ptv.Pt(); ++fIsoVeto05; }
    if ((dr < 0.6) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum06 += lostTrack*ptv.Pt(); ++fIsoVeto06; }
    if ((dr < 0.8) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum08 += lostTrack*ptv.Pt(); ++fIsoVeto08; }
    if ((dr < 1.0) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum10 += lostTrack*ptv.Pt(); ++fIsoVeto10; }
    if ((dr < 1.1) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum11 += lostTrack*ptv.Pt(); ++fIsoVeto11; }
    if ((dr < 1.2) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum12 += lostTrack*ptv.Pt(); ++fIsoVeto12; }
    if ((dr < 1.4) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { sum14 += lostTrack*ptv.Pt(); ++fIsoVeto14; }

    if ((ptv.Pt() > ISOPTLO) && (dr < cone)  && (absQ > 0) ) {
      ++fIsoVeto;
      //      cout << "  ... isolation veto triggered" << endl;
    }

    if ((ptv.Pt() > ISOPTLO) && (dr < drmin) && (absQ > 0) ) {
      drmin = dr;
    }
  }

  fIso05 = fPt/(fPt + sum05);
  fIso06 = fPt/(fPt + sum06);
  fIso08 = fPt/(fPt + sum08);
  fIso10 = fPt/(fPt + sum10);
  fIso11 = fPt/(fPt + sum11);
  fIso12 = fPt/(fPt + sum12);
  fIso14 = fPt/(fPt + sum14);

  fIso = fPt/(fPt + sum);

  // -- Choose cone size for isolation 
  if (TMath::Abs(ISOCONE - 0.5) < 0.001) {
    fIso = fIso05;
  } else if (TMath::Abs(ISOCONE - 0.6) < 0.001) {
    fIso = fIso06;
  } else if (TMath::Abs(ISOCONE - 0.8) < 0.001) {
    fIso = fIso08;
  } else if (TMath::Abs(ISOCONE - 1.0) < 0.001) {
    fIso = fIso10;
  } else if (TMath::Abs(ISOCONE - 1.2) < 0.001) {
    fIso = fIso11;
  } else if (TMath::Abs(ISOCONE - 1.1) < 0.001) {
    fIso = fIso12;
  } else if (TMath::Abs(ISOCONE - 1.4) < 0.001) {
    fIso = fIso14;
  } else {
    cout << "don't know about cone size " << ISOCONE << endl;
    exit(1);
  }

  // -- Systematics: Vertexing 
  if (0) {
    double smear, sigma;
    sigma = 0.3*fLxy;
    
    smear = gRandom->Gaus(0., sigma);
    fLxy += smear;
  }

  // -- Isolation plots
  ((TH2D*)gDirectory->Get("I100"))->Fill(fMass, fRMM);
  ((TH1D*)gDirectory->Get("i100"))->Fill(drmin);

  ((TH1D*)gDirectory->Get("i104"))->Fill(fDphi);
  ((TH1D*)gDirectory->Get("i105"))->Fill(fDeta);
  ((TH1D*)gDirectory->Get("i106"))->Fill(fRMM);

  ((TH1D*)gDirectory->Get("I104"))->Fill(fDeta, fDphi);
  ((TH1D*)gDirectory->Get("I105"))->Fill(fDphi, fRMM);
  ((TH1D*)gDirectory->Get("I106"))->Fill(fDeta, fRMM);


  if (fIsoVeto == 0) {

    ((TH2D*)gDirectory->Get("I101"))->Fill(fMass, fRMM);
    ((TH2D*)gDirectory->Get("I102"))->Fill(cone, drmin);

    ((TH1D*)gDirectory->Get("i101"))->Fill(drmin);

    if (fMass > 4.) {

      ((TH1D*)gDirectory->Get("i102"))->Fill(drmin);
    }
    if (fMass > 2.8 && fMass < 3.2) {

      ((TH1D*)gDirectory->Get("i103"))->Fill(drmin);
    }
  }

  // -- Fill histograms
  h = (TH1D*)gDirectory->Get("v100"); h->Fill(fChi2);
  h = (TH1D*)gDirectory->Get("v101"); h->Fill(fProb);

  h = (TH1D*)gDirectory->Get("v110"); h->Fill(fLxy);
  h = (TH1D*)gDirectory->Get("v111"); h->Fill(fSxy);
  h = (TH1D*)gDirectory->Get("v112"); h->Fill(fLxy/fSxy);

  h = (TH1D*)gDirectory->Get("v120"); h->Fill(fL3d);
  h = (TH1D*)gDirectory->Get("v121"); h->Fill(fS3d);
  h = (TH1D*)gDirectory->Get("v122"); h->Fill(fL3d/fS3d);

  h = (TH1D*)gDirectory->Get("b100"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b101"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b102"); h->Fill(fPt);
  h = (TH1D*)gDirectory->Get("b103"); h->Fill(fEta);
  h = (TH1D*)gDirectory->Get("b104"); h->Fill(fTau);

  if (fNorm) {

    h = (TH1D*)gDirectory->Get("j100"); h->Fill(fMassJ);
    h = (TH1D*)gDirectory->Get("j101"); h->Fill(fMassJ);
    h = (TH1D*)gDirectory->Get("j102"); h->Fill(fPtJ);
    h = (TH1D*)gDirectory->Get("j103"); h->Fill(fEtaJ);
  }

  // -- Angle between B (or J/Psi) momentum and sv direction
  h = (TH1D*)gDirectory->Get("b110"); h->Fill(fCosAngle);

  if (fPt > 5) {
    h = (TH1D*)gDirectory->Get("b120"); h->Fill(fCosAngle);
  }

  if (fL3d/fS3d > 4.) {
    h = (TH1D*)gDirectory->Get("b130"); h->Fill(fCosAngle);
  }

  if (fL3d/fS3d > 10.) {
    h = (TH1D*)gDirectory->Get("b131"); h->Fill(fCosAngle);
  }

  if (fLxy/fSxy > 4.) {
    h = (TH1D*)gDirectory->Get("b140"); h->Fill(fCosAngle);
  }

  if (fLxy/fSxy > 10.) {
    h = (TH1D*)gDirectory->Get("b141"); h->Fill(fCosAngle);
  }

  h = (TH1D*)gDirectory->Get("b150"); h->Fill(fRMM);
  h = (TH1D*)gDirectory->Get("b200"); h->Fill(fIso);
  h = (TH1D*)gDirectory->Get("b201"); h->Fill(fIsoVeto);

  if (fNorm) {

    h = (TH1D*)gDirectory->Get("j150"); h->Fill(fRMM);
    h = (TH1D*)gDirectory->Get("j151"); h->Fill(fRKJ);
    h = (TH1D*)gDirectory->Get("j152"); h->Fill(fRMJ1);
    h = (TH1D*)gDirectory->Get("j153"); h->Fill(fRMJ2);
  }

  if (fNorm) {
    kaonCandidates();
  }

  if (fDebug & 2) { cout << "candidateProperties> End" << endl; }
}

// ----------------------------------------------------------------------
void treeBmm::kaonCandidates() {

  if (fDebug & 2) { cout << "kaonCandidates> Start" << endl; }

  TAnaCand   *pCand;  
  TAnaTrack  *pTrack;

  TGenCand *pG;
  TGenCand *pM;

  int gIndex(-1), gPDG(-1), mPDG(-1);
  int tmKaons(0), tmkIndex(-1);

  for (unsigned int i = 0; i < fKCI.size(); i++) {

    pCand  = fpEvt->getCand(fKCI[i]);
    pTrack = fpEvt->getRecTrack(fKTI[i]);
    gIndex = pTrack->fGenIndex;

    if ( gIndex > 1 ) { 

      pG  = fpEvt->getGenCand(gIndex);

      if( genPartType(pG->fMom1) == 521 ) {

	tmkIndex = i;
	tmKaons++;
      } 
    }
  }

  ((TH1D*)gDirectory->Get("keff"))->Fill(0.1 + tmKaons);

  if ( fSTI[2] >= 0 ) {  

    ((TH1D*)gDirectory->Get("keff"))->Fill(10.1);

    if ( tmkIndex >= 0) {

      if ( tmkIndex == fSTI[2] ) {
	((TH1D*)gDirectory->Get("keff"))->Fill(11.1);
      } else {
	((TH1D*)gDirectory->Get("keff"))->Fill(12.1);
      }

      pCand  = fpEvt->getCand(fKCI[tmkIndex]);
      pTrack = fpEvt->getRecTrack(fKTI[tmkIndex]);

      ((TH1D*)gDirectory->Get("k100"))->Fill(pCand->fMass);
      
      ((TH1D*)gDirectory->Get("k101"))->Fill(pCand->fVtx.fChi2);
      ((TH1D*)gDirectory->Get("k102"))->Fill(pCand->fPlab.Pt());
      ((TH1D*)gDirectory->Get("k103"))->Fill(pCand->fPlab.Eta());
      ((TH1D*)gDirectory->Get("k104"))->Fill(pTrack->fPlab.Pt());
      ((TH1D*)gDirectory->Get("k105"))->Fill(pTrack->fPlab.Eta());
      
      gIndex = fpK->fGenIndex;
      
      if ( gIndex > 1 ) { 
	
	pG  = fpEvt->getGenCand(gIndex); 
	gPDG = pG->fID; 
	
	pM  = fpEvt->getGenCand(pG->fMom1);
	mPDG = pM->fID;
	
	if ( tmkIndex == fSTI[2] ) {

	  ((TH2D*)gDirectory->Get("K100"))->Fill(gPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K101"))->Fill(gPDG, pCand->fMass  - fpB->fMass);
	  ((TH2D*)gDirectory->Get("K102"))->Fill(mPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K103"))->Fill(mPDG, pCand->fMass  - fpB->fMass);

	} else {

	  ((TH2D*)gDirectory->Get("K200"))->Fill(gPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K201"))->Fill(gPDG, pCand->fMass  - fpB->fMass);
	  ((TH2D*)gDirectory->Get("K202"))->Fill(mPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K203"))->Fill(mPDG, pCand->fMass  - fpB->fMass);
	}
      }
    }
  }

  if (fDebug & 2) { cout << "kaonCandidates> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- Cuts on tracks (pt, eta, tip etc)
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::trackSelection(int o) {
  // offset: 400

  if (fDebug & 2) { cout << "trackSelection> Start" << endl; }

  fGoodTT = 1;

  fpHistFile->cd();
  TH1D *h;
  
  int ntrk = 2;
  if ( fD2 > 0 ) {
    
    ntrk = 3;
  }

  for (int it = 0; it < ntrk; ++it) {

    TAnaTrack *tt = fpEvt->getRecTrack(fSTI[it]);
    double x(0.);
    
    fER1->Fill(o + 20 + it*10. + 0.1);

    // -- pT
    x = tt->fPlab.Pt();
    h = (TH1D*)gDirectory->Get("t100"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t10%d", it+1)); h->Fill(x);
    if (x < PTLO) fGoodTT = -1;
    if (x > PTHI) fGoodTT = -2;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 1.1);
    if ((fGoodTT < 0) && (fDebug & 0x8)) {
      break;
    }

    // -- Eta
    x = tt->fPlab.Eta();
    h = (TH1D*)gDirectory->Get("t110"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t11%d", it+1)); h->Fill(x);
    if (x < ETALO) fGoodTT = -3;
    if (x > ETAHI) fGoodTT = -4;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 2.1);
    if ((fGoodTT < 0) && (fDebug & 0x8)) {
      break;
    }

    // -- Tip
    x = tt->fTip;
    h = (TH1D*)gDirectory->Get("t120"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t12%d", it+1)); h->Fill(x);
    if (x > TIPHI) fGoodTT = -5;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 3.1);
    if ((fGoodTT < 0) && (fDebug & 0x8)) {
      break;
    }
  }

  // -- Eff. histogram
  fER1->Fill(o + 0.1);
  if (fGoodTT < 0) {
    fER1->Fill(o + 12.1);
  } else {
    fER1->Fill(o + 11.1);
  }

  if (fDebug & 2) { cout << "trackSelection> End" << endl; }
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// -- Track properties
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void treeBmm::trackProperties() {

  if (fDebug & 2) { cout << "trackProperties> Start" << endl; }


  // -- fMuID = 1 if recontructed as muons, 0 else
  fMuL0 = fpL1->fMuID;
  fMuL1 = fpL2->fMuID;
  if (fNorm) fMuK = fpK->fMuID;

  // -- fMcL0 = particle PDG ID
  fMcL0 = fpL1->fMCID;
  fMcL1 = fpL2->fMCID;
  if (fNorm) fMcK = fpK->fMCID;

  ((TH1D*)gDirectory->Get("mid"))->Fill(fMuL0);
  ((TH1D*)gDirectory->Get("mid"))->Fill(fMuL1);
  if (fNorm) ((TH1D*)gDirectory->Get("kid"))->Fill(fMuK);

  // --Mother and grand-mother particles for muons
  if ( fpL1->fGenIndex > 1 && fpL2->fGenIndex > 1 ) {
    
    TGenCand *p1 = fpEvt->getGenCand(fpL1->fGenIndex);
    TGenCand *p2 = fpEvt->getGenCand(fpL2->fGenIndex);

    int type1(-1),   type2(-2);
    int gtype1(-1), gtype2(-2);
    int mom_i1(-1), mom_i2(-2);
    int gmo_i1(-1), gmo_i2(-2);
    
    if (p1->fMom1 >= 0) {

      TGenCand *m1    = fpEvt->getGenCand(p1->fMom1);
      fMomL0          = m1->fID;;
      mom_i1          = m1->fNumber;

      type1           = genPartType(mom_i1);

      if (m1->fMom1 >= 0) {

	TGenCand *gm1    = fpEvt->getGenCand(m1->fMom1);
	fGMoL0           = gm1->fID;
	gmo_i1           = gm1->fNumber;

	gtype1           = genPartType(gmo_i1);
      }
    }

    if (p2->fMom1 >= 0) {

      TGenCand *m2    = fpEvt->getGenCand(p2->fMom1);
      fMomL1          = m2->fID;
      mom_i2          = m2->fNumber;

      type2           = genPartType(mom_i2);

      if (m2->fMom1 >= 0) {

	TGenCand *gm2    = fpEvt->getGenCand(m2->fMom1);
	fGMoL1           = gm2->fID;
	gmo_i2           = gm2->fNumber;

	gtype2          = genPartType(gmo_i2);
      }
    }    
    
    // --Mother and grand-mother particles for muons
    int type3(-3);
    int gtype3(-3);
    int mom_i3(-3);
    int gmo_i3(-3);

    if (fNorm && fpK->fGenIndex > 1 ) {
      
      TGenCand *p3 = fpEvt->getGenCand(fpK->fGenIndex);
      
      
      if (p3->fMom1 >= 0) {
	
	TGenCand *m3    = fpEvt->getGenCand(p3->fMom1);
	fMomK           = m3->fID;;
	mom_i3          = m3->fNumber;
	
	type3           = genPartType(mom_i3);
	
	if (m3->fMom1 >= 0) {
	  
	  TGenCand *gm3    = fpEvt->getGenCand(m3->fMom1);
	  fGMoK            = gm3->fID;
	  gmo_i3           = gm3->fNumber;
	  
	  gtype3           = genPartType(gmo_i3);
	}
      }
    }
      
    // -- Check if truth-matched 
    if (!fNorm) {
      
      if (type1 == 531) {
	fTruthL0 = 1;
      } else {
	fTruthL0 = 0;
      }
      
      if (type2 == 531) {
	fTruthL1 = 1;
      } else {
	fTruthL1 = 0;
      }
      
      if (fTruthL0 && fTruthL1 && (mom_i1 == mom_i2) ) {
	fTruthB = 1;
      }  else {
	fTruthB = 0;
      }
    }

    if (fNorm) {

      if (type1 == 443 && gtype1 == 521) {
	fTruthL0 = 1;
      } else {
	fTruthL0 = 0;
      }
      
      if (type2 == 443 && gtype2 == 521) {
	fTruthL1 = 1;
      } else {
	fTruthL1 = 0;
      }

      if (type3 == 521) {
	fTruthK = 1;
      } else {
	fTruthK = 0;
      } 
      
      if (fTruthL0 && fTruthL1 && (gmo_i1 == gmo_i2) ) {
	fTruthJ = 1;
      }  else {
	fTruthJ = 0;
      }

      if (fTruthJ && fTruthK &&  (gmo_i1 == mom_i3) ) {
	fTruthB = 1;
      }  else {
	fTruthB = 0;
      }
    }
  }

  if (fDebug & 2) { cout << "trackProperties> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::fillHist() {

  if (fDebug & 2) { cout << "fillHist> Start" << endl; }
  
  TH1D *ho = (TH1D*)gDirectory->Get("oeff");

  // -- offline HLT emulation
  double HLTlPT(4.0);
  double HLTlETA(2.4);
  double HLTlTIP(10000.1);
  double HLTCHI2(20.);
  double HLTL3D(0.015);
  double HLTMASSLO(4.8);
  double HLTMASSHI(6.0);

  ho->Fill(0.1);
  if (1 == fGoodKinematics) {
    ho->Fill(1.1);
    if ((1 == fGoodKinematics) && (1 == fGoodL1)) {
      ho->Fill(2.1);
      if ((fPtL0 > HLTlPT) && (fPtL1 > HLTlPT)) {
	ho->Fill(3.1);
	if ((fEtaL0 > -HLTlETA) && (fEtaL0 < HLTlETA) && (fEtaL1 > -HLTlETA) && (fEtaL1 < HLTlETA)) {
	  ho->Fill(4.1);
	  if (fQL0*fQL1 < -0.5) {
	    ho->Fill(5.1);
	    if ((fTipL0 < HLTlTIP) && (fTipL1 < HLTlTIP)) {
	      ho->Fill(6.1);
	      if (fChi2 < HLTCHI2) {
		ho->Fill(7.1);
		if (fL3d > HLTL3D) {
		  ho->Fill(8.1);
		  if ((fMass > HLTMASSLO) && (fMass < HLTMASSHI)) {
		    ho->Fill(9.1);
		  }      
		}
	      }
	    }
	  }
	}
      }
    }    
  }

  if (fDebug & 2) { cout << "fillHist> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::fillTriggerEff() {

  if (fDebug & 2) { cout << "fillTriggerEff> Start" << endl; }

  if (1) {
    fAR1->Fill(0.1);
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(0.1), "no cuts (0.1)");
  }

  // -- generator kinematics
  if (1 == fGoodKinematics) {
    fAR1->Fill(1.1);  
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(1.1), "kinematics (1.1)");
  } else {
    fAR1->Fill(2.1);
  }

  // -- L1-trigger
  if (1 == fGoodKinematics && 1 == fGoodL1) {
    fAR1->Fill(11.1); 
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(11.1), "L1 (11.1)");
  } else {
    fAR1->Fill(12.1);
  }

  // -- High-level trigger
  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT) {
    fAR1->Fill(21.1);
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(21.1), "HLT (21.1)");
  } else {
    fAR1->Fill(22.1);
  }

  // -- Primary vertex
  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT  && 1 == fGoodPV) {
    fAR1->Fill(31.1);
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(31.1), "PV (31.1)");
  } else {
    fAR1->Fill(32.1);
  }

  // -- Candidate
  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT  && 1 == fGoodCand) {
    fAR1->Fill(41.1);
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(41.1), "Cand (41.1)");
  } else {
    fAR1->Fill(42.1);
  }

  // -- Event
  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT  && 1 == fGoodEvent) {
    fAR1->Fill(51.1);
    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(51.1), "Event (51.1)");
  } else {
    fAR1->Fill(52.1);
  }


  if (fDebug & 2) { cout << "fillTriggerEff> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::fillAnalysisEff() {

  if (fDebug & 2) { cout << "fillAnalysisEff> Start" << endl; }

  fGoodAna  = fGoodAnaF = fGoodVtxF = fGoodIsoF   = 0;
  fGoodPtMu = fGoodRmm  = fGoodMass = fGoodWindow = 0;
  fGoodPtB  = fGoodEtaB = fGoodCosa = fGoodLength = fGood3DLength = 0;
  fGoodIso  = fGoodVtx  = 0;

  //--------------------------------------------------------------------
  // -- selection histogram (0-3)
  //--------------------------------------------------------------------
  if (1) {
    histogram(0);
  }

  // -- generator kinematics (histo: 1)
  if (1 == fGoodKinematics) {
    histogram(1); 
  }

  // -- L1-trigger (histo: 2)
  if (1 == fGoodKinematics && 1 == fGoodL1) {
    histogram(2);
  }

  // -- High-level trigger (histo: 3)
  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT) {
    histogram(3);
  }

  //--------------------------------------------------------------------
  // -- single cut effiencies (offline)
  //--------------------------------------------------------------------

  if ((fPtL0 > PTLO) && (fPtL1 > PTLO))       { fGoodPtMu     = 1; }
  if ((RMMLO < fRMM) && (fRMM < RMMHI))       { fGoodRmm      = 1; }
  if ((MASSLO < fMass) && (fMass < MASSHI))   { fGoodMass     = 1; }
  if ( TMath::Abs(fMass - fMassB) < MASSWI)   { fGoodWindow   = 1; }
  if (fPt > PTBS)                             { fGoodPtB      = 1; }
  if (ETALO < fEta && fEta < ETAHI)           { fGoodEtaB     = 1; }
  if (fCosAngle > COSALPHA)                   { fGoodCosa     = 1; }
  if (fLxy/fSxy > LXYSXYLO)                   { fGoodLength   = 1; }
  if (fL3d/fS3d > L3DS3DLO)                   { fGood3DLength = 1; }
  if (fIso > ISOLATION)                       { fGoodIso      = 1; }
  if (fChi2 < VTXCHI)                         { fGoodVtx      = 1; }
  if ( fLxy/fSxy > 7 )                        { fGoodPresel   = 1; }


  //--------------------------------------------------------------------
  // -- Set Triggers to true for rare BG (for next part)
  //--------------------------------------------------------------------
  int ignoreTrigger = 0;
  if ( fD0 != 13 || fD1 != 13 ) {
    if (fDebug & 1) cout << "Ignore L1/HLT decision " << endl;
    ignoreTrigger = 1;
  }

  //--------------------------------------------------------------------
  // -- cumulative cut effiencies (offline)
  //--------------------------------------------------------------------
  if ( (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT) || ignoreTrigger ) {

    ((TH1D*)gDirectory->Get("PTLO"))->Fill(fPt);

    if ((fPtL0 > PTLO) && (fPtL1 > PTLO)) {
      fAR1->Fill(100.1); // histogram(10);
      fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(100.1), "p_{T}(#mu) (100.1)");


      ((TH1D*)gDirectory->Get("RMM"))->Fill(fRMM);

      if ((RMMLO < fRMM) && (fRMM < RMMHI)) {
	fAR1->Fill(101.1); // histogram(11);
	fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(101.1), "#Delta R(#mu #mu) (101.1)");

	((TH1D*)gDirectory->Get("MASSBAND"))->Fill(fMass);
	
	if ( (MASSLO < fMass) && (fMass < MASSHI) ) {
	  fAR1->Fill(110.1); 
	  fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(110.1), "mass (110.1)");
	  	    
	  // -- offline cuts (w/o factorization) 
	  // ----------------------------------

	  ((TH1D*)gDirectory->Get("PTBS"))->Fill(fPt);

	  if (fPt > PTBS) {
	    fAR1->Fill(120.1); // histogram(12);
	    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(120.1), "pT(B) (120.1)");
	    
	    ((TH1D*)gDirectory->Get("ETABS"))->Fill(fEta);
	      
	    if (ETALO < fEta && fEta < ETAHI) {
	      fAR1->Fill(121.1); // histogram(13);
	      fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(121.1), "eta(B) (121.1)");

	      ((TH1D*)gDirectory->Get("COSALPHA"))->Fill(fCosAngle);
		
	      if (fCosAngle > COSALPHA) {
		fAR1->Fill(122.1); // histogram(14);
		fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(122.1), "cos(#alpha_{xy}) (122.1)");

		((TH1D*)gDirectory->Get("LXYSXY"))->Fill(fLxy/fSxy);
		((TH1D*)gDirectory->Get("L3DS3D"))->Fill(fL3d/fS3d);
		  
		if (fL3d/fS3d > L3DS3DLO && fLxy/fSxy > LXYSXYLO) {
		  fAR1->Fill(123.1); //histogram(15); 
		  fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(123.1), "l_{xy}/#sigma_{xy} (123.1)");
		  fGoodAnaF = 1;  // <-------- good ana.     (before factorizing)
		  histogram(5);   // <-------- distributions (before factorizing)

		  ((TH1D*)gDirectory->Get("MASSWI_F"))->Fill(fMass);
		    
		  // -- signal window before 'factorizing' cuts
		  if ( TMath::Abs(fMass - fMassB) < MASSWI ) {
		    fAR1->Fill(210.1);
		    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(210.1), "m_{B} #pm 100 MeV (w/o fact) (210.1)");
		  }
		  
		  ((TH1D*)gDirectory->Get("ISOLATION"))->Fill(fIso);

		  ((TH1D*)gDirectory->Get("ISOLATION_F"))->Fill(fIso);
		  ((TH1D*)gDirectory->Get("VTXCHI_F"))->Fill(fChi2);


		  // -- 'factorizing' cuts
		  if (fIso > ISOLATION) {
		    fAR1->Fill(124.1); // histogram(16);
		    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(124.1), "Iso (124.1)");

		    ((TH1D*)gDirectory->Get("VTXCHI"))->Fill(fChi2);

		    if (fChi2 < VTXCHI) {
		      fAR1->Fill(125.1);  
		      fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(125.1), "Vertex (125.1)");
		      fGoodAna = 1;  // <-------- good ana.
		      histogram(4);  // <-------- distributions
		      
		      ((TH1D*)gDirectory->Get("MASSWI"))->Fill(fMass);

		      if ( TMath::Abs(fMass - fMassB) < MASSWI ) {
			fAR1->Fill(126.1);
			fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(126.1), "m_{B} #pm 100 MeV (126.1)");
		      }
		    }
		  }
		}
	      }
	    }
	  }

	  // -- offline cuts (w/ factorization)
	  // ----------------------------------
	    
	  // -- somewhat tightened preselection 
	  ((TH1D*)gDirectory->Get("LXYSXY_PRE"))->Fill(fLxy/fSxy);
	  ((TH1D*)gDirectory->Get("L3DS3D_PRE"))->Fill(fL3d/fS3d);

	  if ( (L3DS3DLO > 0 && fL3d/fS3d > PRESEL) || 
	       (LXYSXYLO > 0 && fLxy/fSxy > PRESEL) ) { // && fCosAngle > 0.9  
	      
	    fAR1->Fill(220.1); // histogram(20);
	    fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(220.1), "preselection (fact) (220.1)");
	      
	    ((TH1D*)gDirectory->Get("ISOLATION_PRE"))->Fill(fIso);
	    ((TH1D*)gDirectory->Get("VTXCHI_PRE"))->Fill(fChi2);

	    // -- 'factorizing' cuts
	    if (fIso > ISOLATION) {
	      fAR1->Fill(224.1); // histogram(21);
	      fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(224.1), "Iso (fact) (224.1)");
	      fGoodIsoF = 1;  // <-------- good iso. (factorizing)
	    }
	      
	    if (fChi2 < VTXCHI) {
	      fAR1->Fill(225.1); // histogram(22);
	      fAR1->GetXaxis()->SetBinLabel(fAR1->FindBin(225.1), "Vertex (fact) (225.1)");
	      fGoodVtxF = 1;  // <-------- good vtx. (factorizing)
	    }
	  }
	}
      }
    }
  }

  if (fDebug & 2) { cout << "fillAnalysisEff> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::offlineEff(int bs_index, int sel) {

  if (fDebug & 2) { cout << "offlineEff> Start" << endl; }

  if ( bs_index < 0 ) {

    return;
  }

  TAnaTrack  *L1, *L2;
  TAnaCand   *B;
  int        STI[3];         // Index of signal tracks

  int hist_i = sel -1;

  // -- Secondary vertex variables
  double L3d(-99.), S3d(-99.), Lxy(-99.), Sxy(-99.), Chi2(-99.);

  // -- Variables for muon tracks / kaon track
  double PtL0(-99.), PtL1(-99.);

  // -- Variables for the chosen B candidate
  double Mass(-99.), Pt(-99.), Eta(-99.), RMM(-99.), CosAngle(-99.), Iso(-99.);

  B   = fpEvt->getCand(bs_index);
    
  // - B cand.
  Chi2 = B->fVtx.fChi2;
  
  L3d  = B->fVtx.fD3d;
  S3d  = B->fVtx.fD3dE;
  Lxy  = B->fVtx.fDxy;
  Sxy  = B->fVtx.fDxyE;
  
  Mass = B->fMass;
  Pt   = B->fPlab.Pt();

  // - muon cand.
  STI[0] = B->fSig1;
  STI[1] = B->fSig2;
  L1     = fpEvt->getRecTrack(STI[0]);
  L2     = fpEvt->getRecTrack(STI[1]);
   
  PtL0 = L1->fPlab.Pt();
  PtL1 = L2->fPlab.Pt();

  double dphi = L1->fPlab.DeltaPhi(L2->fPlab);
  double deta = L1->fPlab.Eta() - L2->fPlab.Eta();
  
  RMM  = TMath::Sqrt(dphi*dphi + deta*deta);
  
  // - B cand.
  Pt   = B->fPlab.Pt();
  Eta  = B->fPlab.Eta();
  
  L3d  = B->fVtx.fD3d;
  S3d  = B->fVtx.fD3dE;
  Lxy  = B->fVtx.fDxy;
  Sxy  = B->fVtx.fDxyE;
  
  // -- 3D version: l is the vector from the PV to the SV
  if (-99 == fpEvt->fPrimaryVertex2.fType) {
    
    return;
  }

  // -- Pointing angle
  TVector2 tPV = fpEvt->fPrimaryVertex2.fPoint.XYvector();
  TVector2 tSV = B->fVtx.fPoint.XYvector();
  TVector2 tB  = B->fPlab.XYvector();
  TVector2 tl  = tSV - tPV;
  CosAngle = TMath::Cos(tB.DeltaPhi(tl));

  
  // -- Isolation
  TVector3 ptv;
  int absQ(0);
  double df2, de2;
  double dr(0.);
  double sum(0.);
  
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
    
    ptv   = fpEvt->getRecTrack(it)->fPlab;
    absQ   = TMath::Abs(fpEvt->getRecTrack(it)->fQ);
    
    if (it == STI[0]) { continue; }
    if (it == STI[1]) { continue; }
    if ( fNorm && it == STI[2]) { continue; }
    
    df2 = B->fPlab.DeltaPhi(ptv);
    de2 = (ptv.Eta() - B->fPlab.Eta());
    dr  = TMath::Sqrt(df2*df2 + de2*de2);
    
    if ( (dr < ISOCONE) && (ptv.Pt() > ISOPTLO) && (absQ > 0) ) { 

      sum += ptv.Pt();
    }
  }
  
  Iso = Pt/(Pt + sum);

  //--------------------------------------------------------------------
  // -- cumulative cut effiencies (offline)
  //--------------------------------------------------------------------

  if (1 == fGoodKinematics && 1 == fGoodL1  && 1 == fGoodHLT) {

    if ((PtL0 > PTLO) && (PtL1 > PTLO)) {
      fOR1[hist_i]->Fill(100.1);
      fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(100.1), "p_{T}(#mu) (100.1)");

      if ((RMMLO < RMM) && (RMM < RMMHI)) {
	fOR1[hist_i]->Fill(101.1); 
	fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(101.1), "#Delta R(#mu #mu) (101.1");
	
	if ( (MASSLO < Mass) && (Mass < MASSHI) ) {
	  fOR1[hist_i]->Fill(110.1); 
	  fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(110.1), "mass (110.1)");

	  // -- offline cuts (w/o factorization) 
	  // ----------------------------------

	  if (Pt > PTBS) {
	    fOR1[hist_i]->Fill(120.1);
	    fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(120.1), "pT(B) (120.1)");
	  
	    if (ETALO < Eta && Eta < ETAHI) {
	      fOR1[hist_i]->Fill(121.1);
	      fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(121.1), "eta(B) (121.1)");
	    
	      if (CosAngle > COSALPHA) {
		fOR1[hist_i]->Fill(122.1);
		fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(122.1), "cos(#alpha) (122.1)");
	      
		if (L3d/S3d > L3DS3DLO && Lxy/Sxy > LXYSXYLO) {
		  fOR1[hist_i]->Fill(123.1);
		  fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(123.1), "l_{xy}/#sigma_{xy} (123.1)");
		
		  // -- signal window before 'factorizing' cuts
		  if ( TMath::Abs(Mass - fMassB) < MASSWI ) {
		    fOR1[hist_i]->Fill(210.1);
		    fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(210.1),
							  "m_{B} #pm 100 MeV (w/o fact) (210.1)");
		  }
		
		  // -- 'factorizing' cuts
		  if (Iso > ISOLATION) {
		    fOR1[hist_i]->Fill(124.1);
		    fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(124.1), "Iso (124.1)");
		    if (Chi2 < VTXCHI) {
		      fOR1[hist_i]->Fill(125.1);
		      fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(125.1), "Vertex (125.1)");
		      if ( TMath::Abs(Mass - fMassB) < MASSWI ) {
			fOR1[hist_i]->Fill(126.1);
			fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(126.1), 
							      "m_{B} #pm 100 MeV (126.1)");
		      }
		    }
		  }
		}
	      }
	    }
	  }
	
	  // -- offline cuts (w/ factorization)
	  // ----------------------------------
	
	  // -- somewhat tightened preselection 
	  if ( Lxy/Sxy > 7 ) {    // && fCosAngle > 0.9

	    fOR1[hist_i]->Fill(220.1);
	    fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(220.1), "preselection (fact) (220.1)");

	    // -- 'factorizing' cuts
	    if (Iso > ISOLATION) {
	      fOR1[hist_i]->Fill(224.1);
	      fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(224.1), "Iso (fact) (224.1)");
	    }
	    	  
	    if (Chi2 < VTXCHI) {
	      fOR1[hist_i]->Fill(225.1);
	      fOR1[hist_i]->GetXaxis()->SetBinLabel(fOR1[hist_i]->FindBin(225.1), "Vertex (fact) (225.1)");
	    }
	  }
	}
      }
    }  
  }

  if (fDebug & 2) { cout << "offlineEff> End" << endl; }
}


// **************************************************************************************************
// --------------------------------------------------------------------------------------------------
void treeBmm::bookHist() {

  if (fDebug & 2) { cout << "bookHist> Start" << endl; }

  TH1 *h;
  TH2 *h2;
  TProfile *p;

  cout << "-->bookHist> " << endl;

  book(0);   // histogram(0) reco level, before any cuts
  book(1);   // histogram(1) reco level, after kinematic cuts on generator level have been applied
  book(2);   // histogram(2) reco level, after L1 cuts
  book(3);   // histogram(3) reco level, after HLT cuts
  book(4);   // histogram(4) reco level, after all offline cuts
  book(5); // histogram(15) reco level, after offline cuts without factorizing cuts (Vertex&Isolation)

  book2(0);   // histogram(0) reco level, before any cuts
  book2(1);   // histogram(1) reco level, after kinematic cuts on generator level have been applied
  book2(2);   // histogram(2) reco level, after L1 cuts
  book2(3);   // histogram(3) reco level, after HLT cuts
  book2(4);   // histogram(4) reco level, after all offline cuts
  book2(5);   // histogram(15) reco level, after offline cuts without factorizing cuts (Vertex&Isolation)

  // -- single cuts histograms
//   book(10);
//   book(11);
//   book(12);
//   book(13);
//   book(14);
//   book(15); book(5); // histogram(15) reco level, after offline cuts without factorizing cuts (Vertex&Isolation)
//   book(16);
//   book(20);
//   book(21);
//   book(22);

//   book2(10);
//   book2(11);
//   book2(12);
//   book2(13);
//   book2(14);
//   book2(15); book2(5);   // histogram(15) reco level, after offline cuts without factorizing cuts (Vertex&Isolation)
//   book2(16);
//   book2(20);
//   book2(21);
//   book2(22);

  //*******************************************************************************
  // Reduced Tree 
  //*******************************************************************************

  fTree = new TTree("events", "events");

  fTree->Branch("run",            &fRun ,            "run/I");

  fTree->Branch("process",        &fProcessType ,    "process/I");

  fTree->Branch("goodKinematics", &fGoodKinematics , "goodKinematics/I");
  fTree->Branch("goodL1",         &fGoodL1 ,         "goodL1/I");
  fTree->Branch("goodHLT",        &fGoodHLT ,        "goodHLT/I");
  fTree->Branch("goodPV",         &fGoodPV ,         "goodPV/I");
  fTree->Branch("goodCand",       &fGoodCand ,       "goodCand/I");
  fTree->Branch("goodEvent",      &fGoodEvent ,      "goodEvent/I");

  fTree->Branch("goodAna",        &fGoodAna  ,       "goodAna/I");
  fTree->Branch("goodAnaF",       &fGoodAnaF ,       "goodAnaF/I");
  fTree->Branch("goodVtxF",       &fGoodVtxF ,       "goodVtxF/I");
  fTree->Branch("goodIsoF",       &fGoodIsoF ,       "goodIsoF/I");

  fTree->Branch("goodTT",         &fGoodTT ,         "goodTT/I");

  fTree->Branch("goodPtMu",       &fGoodPtMu  ,       "goodPtMu/I");
  fTree->Branch("goodRmm",        &fGoodRmm  ,        "goodRmm/I");
  fTree->Branch("goodMass",       &fGoodMass  ,       "goodMass/I");
  fTree->Branch("goodPtB",        &fGoodPtB  ,        "goodPtB/I");
  fTree->Branch("goodCosa",       &fGoodCosa  ,       "goodCosa/I");
  fTree->Branch("goodLength",     &fGoodLength  ,     "goodLength/I");
  fTree->Branch("good3DLength",   &fGood3DLength  ,   "good3DLength/I");
  fTree->Branch("goodVtx",        &fGoodVtx  ,        "goodVtx/I");
  fTree->Branch("goodIso",        &fGoodIso  ,        "goodIso/I");
  fTree->Branch("goodWindow",     &fGoodWindow  ,     "goodWindow/I");

  fTree->Branch("dVtxPerp",       &fDvtxPerp,        "dVtxPerp/D");
  fTree->Branch("dVtxPar",        &fDvtxPar,         "dVtxPar/D");

  fTree->Branch("dPVtxX",         &fDpVtxX,          "dPVtxX/D");
  fTree->Branch("dPVtxY",         &fDpVtxY,          "dPVtxY/D");
  fTree->Branch("dPVtxZ",         &fDpVtxZ,          "dPVtxZ/D");

  fTree->Branch("dSVtxX",         &fDsVtxX,          "dSVtxX/D");
  fTree->Branch("dSVtxY",         &fDsVtxY,          "dSVtxY/D");
  fTree->Branch("dSVtxZ",         &fDsVtxZ,          "dSVtxZ/D");

  fTree->Branch("rflt",           &fFltR,            "rflt/D");
  fTree->Branch("gflt",           &fFltG,            "gflt/D");
  fTree->Branch("dFltRes",        &fFltRes,          "dFltRes/D");

  fTree->Branch("rflt3d",           &fFltR3D,            "rflt3d/D");
  fTree->Branch("gflt3d",           &fFltG3D,            "gflt3d/D");
  fTree->Branch("dFltRes3d",        &fFltRes3D,          "dFltRes3d/D");

  fTree->Branch("mass",           &fMass,            "mass/D");
  fTree->Branch("pt",             &fPt,              "pt/D");
  fTree->Branch("p",              &fP,               "p/D");
  fTree->Branch("eta",            &fEta,             "eta/D");

  fTree->Branch("ptl0",           &fPtL0,            "ptl0/D");
  fTree->Branch("etal0",          &fEtaL0,           "etal0/D");
  fTree->Branch("ql0",            &fQL0,             "ql0/D");
  fTree->Branch("tipl0",          &fTipL0,           "tipl0/D");

  fTree->Branch("ptl1",           &fPtL1,            "ptl1/D");
  fTree->Branch("etal1",          &fEtaL1,           "etal1/D");
  fTree->Branch("ql1",            &fQL1,             "ql1/D");
  fTree->Branch("tipl1",          &fTipL1,           "tipl1/D");

  fTree->Branch("ptk",            &fPtK,            "ptk/D");
  fTree->Branch("etak",           &fEtaK,           "etak/D");
  fTree->Branch("qk",             &fQK,             "qk/D");
  fTree->Branch("tipk",           &fTipK,           "tipk/D");

  fTree->Branch("glbl0",          &fMuL0,            "glbl0/D");
  fTree->Branch("glbl1",          &fMuL1,            "glbl1/D");
  fTree->Branch("glbk",           &fMuK,              "glbk/D");

  fTree->Branch("pdgl0",          &fMcL0,           "pdgl0/I");
  fTree->Branch("pdgl1",          &fMcL1,           "pdgl1/I");
  fTree->Branch("pdgk",           &fMcK,             "pdgk/I");

  fTree->Branch("moml0",          &fMomL0,           "moml0/I");
  fTree->Branch("moml1",          &fMomL1,           "moml1/I");
  fTree->Branch("momk",           &fMomK,            "momk/I");

  fTree->Branch("gmol0",          &fGMoL0,           "gmol0/I");
  fTree->Branch("gmol1",          &fGMoL1,           "gmol1/I");
  fTree->Branch("gmok",           &fGMoK,           "  gmok/I");

  fTree->Branch("truthl0",        &fTruthL0,       "truthl0/I");
  fTree->Branch("truthl1",        &fTruthL1,       "truthl1/I");
  fTree->Branch("truthk",         &fTruthK,         "truthk/I");

  fTree->Branch("truthJ",         &fTruthJ,         "truthJ/I");
  fTree->Branch("truthB",         &fTruthB,         "truthB/I");

  fTree->Branch("rmm",            &fRMM,            "rmm/D");

  fTree->Branch("rmj1",           &fRMJ1,           "rmj1/D");
  fTree->Branch("rmj2",           &fRMJ2,           "rmj2/D");
  fTree->Branch("rjk",            &fRKJ,            "rkj/D");

  fTree->Branch("chi2",           &fChi2,           "chi2/D");
  fTree->Branch("ndof",           &fNdof,           "ndof/D");
  fTree->Branch("prob",           &fProb,           "prob/D");
  fTree->Branch("cosa3",          &fCosAngle3,      "cosa3/D");
  fTree->Branch("cosa",           &fCosAngle,       "cosa/D");
  fTree->Branch("l3d",            &fL3d,            "l3d/D");
  fTree->Branch("s3d",            &fS3d,            "s3d/D");
  fTree->Branch("lxy",            &fLxy,            "lxy/D");
  fTree->Branch("sxy",            &fSxy,            "sxy/D");
  fTree->Branch("tau",            &fTau,            "tau/D");
  fTree->Branch("txy",            &fTxy,            "txy/D");
  fTree->Branch("dPhi",           &fDphi,           "dPhi/D");
  fTree->Branch("dEta",           &fDeta,           "dEta/D");
  fTree->Branch("i05",            &fIso05,          "i05/D");
  fTree->Branch("i06",            &fIso06,          "i06/D");
  fTree->Branch("i08",            &fIso08,          "i08/D");
  fTree->Branch("i10",            &fIso10,          "i10/D");
  fTree->Branch("i12",            &fIso12,          "i12/D");
  fTree->Branch("i14",            &fIso14,          "i14/D");
  fTree->Branch("iso",            &fIso,            "iso/D");
  fTree->Branch("isoveto",        &fIsoVeto,        "isoveto/I");
  fTree->Branch("iv05",           &fIsoVeto05,      "iv05/I");
  fTree->Branch("iv06",           &fIsoVeto06,      "iv06/I");
  fTree->Branch("iv08",           &fIsoVeto08,      "iv08/I");
  fTree->Branch("iv10",           &fIsoVeto10,      "iv10/I");
  fTree->Branch("iv12",           &fIsoVeto12,      "iv12/I");
  fTree->Branch("iv14",           &fIsoVeto14,      "iv14/I");

  // -- muon counters
  fTree->Branch("nGenMu",         &fgMu,            "nGenMu/I");
  fTree->Branch("nRecMu",         &frMu,            "nRecMu/I");
  fTree->Branch("nGlbMu",         &fnMu,            "nGlbMu/I");

  fTree->Branch("nRecSigMu",      &frSigMu,         "nRecSigMu/I");
  fTree->Branch("nGlbSigMu",      &fnSigMu,         "nGlbSigMu/I");

  fTree->Branch("nGenK",          &fgK,              "nGenK/I");
  fTree->Branch("nRecK",          &frK,              "nRecK/I");

  fTree->Branch("nRecSigK",       &frSigK,           "nRecSigK/I");

  fTree->Branch("nTmMu",          &fnTmMu,              "nTmMu/I");
  fTree->Branch("nTmK",           &fnTmK,                "nTmK/I");

  // -- candidate counters
  fTree->Branch("nJpsiCand",      &fnJ,             "nJpsiCand/I");
  fTree->Branch("nJpsiGenRec",    &fgrJ,          "nJpsiGenRec/I");
  fTree->Branch("nJpsiGen",       &fgJ,              "nJpsiGen/I");
  fTree->Branch("nJmmGen",        &fgJmm,            "nJmmGen/I");
  fTree->Branch("nBcand",         &fnB,                "nBcand/I");
  fTree->Branch("nBgenRec",       &fgrB,             "nBgenRec/I");
  fTree->Branch("nBgen",          &fgB,                 "nBgen/I");
  fTree->Branch("nBmmGen",        &fgBmm,             "nBmmGen/I");

  // -- Gen. signal
  fTree->Branch("gptl0",          &fgPtL0,           "gptl0/D");
  fTree->Branch("getal0",         &fgEtaL0,          "getal0/D");
  fTree->Branch("gql0",           &fgQL0,            "gql0/D");

  fTree->Branch("gptl1",          &fgPtL1,           "gptl1/D");
  fTree->Branch("getal1",         &fgEtaL1,          "getal1/D");
  fTree->Branch("gql1",           &fgQL1,            "gql1/D");

  fTree->Branch("grmm",           &fgRMM,            "grmm/D");
  fTree->Branch("gchi2",          &fgChi2,           "gchi2/D");
  fTree->Branch("gs3d",           &fgS3d,            "gs3d/D");


  //*******************************************************************************
  // Histograms
  //*******************************************************************************
  // -- event
  h = new TH1D("runs", "runs ",         100000, 0., 100000);
  h = new TH1D("mass", "mass ",             70, 5.,    5.7);

  // -- cuts
  h = new TH1D("PTLO", " ",              50,  0.,    25.); h->Sumw2(); setTitles(h,  "p_{T, #mu} [GeV]", "events/bin");
  h = new TH1D("RMM", " ",               50,  0.,     5.); h->Sumw2(); setTitles(h,  "#Delta R(#mu#mu)", "events/bin");
  h = new TH1D("MASSBAND", " ",         400,  2.,    10.); h->Sumw2(); setTitles(h,  "m_{#mu#mu} [GeV]", "events/bin");
  h = new TH1D("PTBS", " ",              60,  0.,    30.); h->Sumw2(); setTitles(h,  "p_{T, B} [GeV]", "events/bin");
  h = new TH1D("ETABS", " ",             50, -5.,     5.); h->Sumw2(); setTitles(h,  "#eta_{B} [GeV]", "events/bin");
  h = new TH1D("COSALPHA", " ",         200,  0.97,   1.); h->Sumw2(); setTitles(h,  "cos #alpha_{xy}", "events/bin");
  h = new TH1D("LXYSXY", " ",            50,  0.,    50.); h->Sumw2(); setTitles(h, "l_{xy}/#sigma_{xy}", "events/bin");
  h = new TH1D("L3DS3D", " ",            50,  0.,    50.); h->Sumw2(); setTitles(h, "l_{3D}/#sigma_{3D}", "events/bin");
  h = new TH1D("ISOLATION", " ",         55,  0.,    1.1); h->Sumw2(); setTitles(h,  "I", "events/bin");
  h = new TH1D("VTXCHI", " ",           200,  0.,    20.); h->Sumw2(); setTitles(h,  "#chi^{2}", "events/bin");
  h = new TH1D("MASSWI", " ",           120, 4.8,     6.); h->Sumw2(); setTitles(h,  "m_{#mu#mu} [GeV]", "events/bin");
 
  h = new TH1D("MASSWI_F", " ",          120, 4.8,     6.); h->Sumw2(); setTitles(h,  "m_{#mu#mu} [GeV]", "events/bin");
  h = new TH1D("ISOLATION_F", " ",        55,  0.,    1.1); h->Sumw2(); setTitles(h,  "I", "events/bin");
  h = new TH1D("VTXCHI_F", " ",          200,  0.,    20.); h->Sumw2(); setTitles(h,  "#chi^{2}", "events/bin");

  h = new TH1D("LXYSXY_PRE", " ",         50,  0.,    50.); h->Sumw2(); setTitles(h,"l_{xy}/#sigma_{xy}", "events/bin");
  h = new TH1D("L3DS3D_PRE", " ",         50,  0.,    50.); h->Sumw2(); setTitles(h,"l_{3D}/#sigma_{3D}", "events/bin");
  h = new TH1D("ISOLATION_PRE", " ",      55,  0.,    1.1); h->Sumw2(); setTitles(h,  "I", "events/bin");
  h = new TH1D("VTXCHI_PRE", " ",        200,  0.,    20.); h->Sumw2(); setTitles(h,  "#chi^{2}", "events/bin");

  // -- data sample check
  fNgen  = new TH1D("ngen", "N_{#mu}^{gen}",   10, 0., 10. );    fNgen->Sumw2();
  fNrec  = new TH1D("nrec", "N_{#mu}^{rec}",   10, 0., 10. );    fNrec->Sumw2();
  fNglb  = new TH1D("nglb", "N_{#mu}^{glb}",   10, 0., 10. );    fNglb->Sumw2();

  fNgenJ  = new TH1D("ngenJ", "N_{J/#psi}^{gen}",   10, 0., 10. );    fNgenJ->Sumw2();
  fNdecJ  = new TH1D("ndecJ", "N_{J/#psi #rightarrow #mu #mu}^{gen}",   10, 0., 10. );    fNdecJ->Sumw2();

  fNgenB  = new TH1D("ngenB", "N_{B}^{gen}",   10, 0., 10. );    fNgenB->Sumw2();
  fNdecB  = new TH1D("ndecB", "N_{B #rightarrow #mu #mu}^{gen}",   10, 0., 10. );    fNdecB->Sumw2();

  fErec  = new TH1D("erec", "#epsilon_{#mu}^{rec}",   10, 0., 1. );    fErec->Sumw2();
  fEglb = new TH1D("eglb", "#epsilon_{#mu}^{glb}",   10, 0., 1. );    fEglb->Sumw2();

  fNR0  = new TH2D("NR0", "Good Event",                     1000, 0., 1000., 1000, 0., 1000.);    fNR0->Sumw2();
  fNR1  = new TH2D("NR1", "N_{B_{s}^{0}}^{cand} / event ",  1000, 0., 1000., 1000, 0., 1000.);    fNR1->Sumw2();
  fNR2  = new TH2D("NR2", "N_{B_{s}^{0}}^{gen} / event ",   1000, 0., 1000., 1000, 0., 1000.);    fNR2->Sumw2();
  fNR3  = new TH2D("NR3", "N_{#mu}^{rec} / event",          1000, 0., 1000., 1000, 0., 1000.);    fNR3->Sumw2();
  fNR4  = new TH2D("NR4", "N_{#mu}^{glb} / event  ",        1000, 0., 1000., 1000, 0., 1000.);    fNR4->Sumw2();
  fNR5  = new TH2D("NR5", "L1 trigger / event  ",           1000, 0., 1000., 1000, 0., 1000.);    fNR5->Sumw2();
  fNR6  = new TH2D("NR6", "HLT trigger / event  ",          1000, 0., 1000., 1000, 0., 1000.);    fNR6->Sumw2();

  // -- event reduction
  fER1  = new TH1D("ER1", "Event Reduction ",                 1000, 0., 1000.);    fER1->Sumw2();
  fAR1  = new TH1D("AR1", "Analysis Reduction ",              1000, 0., 1000.);    fAR1->Sumw2();

  // -- offline hlt-emulation and efficiencies
  h = new TH1D("oeff", "Offline HLT efficiencies",            1000, 0., 1000.);    h->Sumw2();

  // -- offline event reduction for different candidate selection
  for (int i = 0; i < NSUBSEL; i++ ) {

    fOR1[i] = new TH1D(Form("OR1_%i",i+1), Form("Offline Reduction for subsel = %i",i+1),  1000, 0., 1000.);   fOR1[i]->Sumw2();

  }


  // -- generator Level (only signal)
  h = new TH1D("g100", "generator pT",                          40, 0., 80.);    h->Sumw2();   setTitles(h, "pT(#mu)", "events/bin");   
  h = new TH1D("g102", "generator eta",                       60, -3.0, 3.0);    h->Sumw2();   setTitles(h, "eta(#mu)", "events/bin");
  h = new TH1D("g104", "generator status",                     11, -1., 10.);    h->Sumw2();   setTitles(h, "generator status", "events/bin");
  h = new TH1D("g105", "generator pT, #mu_{1}",                 80, 0., 80.);    h->Sumw2();   setTitles(h, "pT(#mu_{1})", "events/bin");
  h = new TH1D("g106", "generator pT, #mu_{2}",                 80, 0., 80.);    h->Sumw2();   setTitles(h, "pT(#mu_{2})", "events/bin");
  h = new TH1D("g107", "generator eta, #mu_{1}",              60, -3.0, 3.0);    h->Sumw2();   setTitles(h, "eta(#mu_{1})", "events/bin");
  h = new TH1D("g108", "generator eta  #mu_{2}",              60, -3.0, 3.0);    h->Sumw2();   setTitles(h, "eta(#mu_{2})", "events/bin");
  h = new TH1D("g110", "generator pT, kaon ",                   80, 0., 40.);    h->Sumw2();   setTitles(h, "pT(K)", "events/bin");
  h = new TH1D("g112", "generator eta, kaon ",                 60, -3.0, 3.);    h->Sumw2();   setTitles(h, "eta(K)", "events/bin");
  h = new TH1D("g120", "deltaRmm",                              50, 0., 5.0);    h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");
  h = new TH1D("g121", "deltaRkj",                              50, 0., 5.0);    h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi K)", "events/bin");
  h = new TH1D("g122", "deltaRmj1",                             50, 0., 5.0);    h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{1})", "events/bin");
  h = new TH1D("g123", "deltaRmj2",                             50, 0., 5.0);    h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{2})", "events/bin");
  h = new TH1D("g130", "generator mass, m_{#mu#mu}",          1000, 0., 10.);    h->Sumw2();   setTitles(h, "m_{#mu#mu}", "events/bin");
  h = new TH1D("g140", "generator mass, m_{#mu#muK}",        1000, 0., 10.);    h->Sumw2();   setTitles(h, "m_{J/#Psi K}", "events/bin");

  // -- signal plots
  // - 1D
  h = new TH1D("s100", "#Delta #chi^{2}",                    100, -10., 10.);    h->Sumw2();    setTitles(h, "#Delta #chi^{2}", "events/bin");
  h = new TH1D("s101", "#Delta p_{T,leading} [GeV]",         100, -25., 25.);    h->Sumw2();    setTitles(h, "#Delta p_{T,leading}", "events/bin");
  h = new TH1D("s102", "#DeltaR(#mu_{1}#mu_{2}) - #DeltaR_{min}(#mu#mu)",  
                                                               100, -5., 5.);    h->Sumw2();    setTitles(h, "#DeltaR(#mu_{1}#mu_{2}) - #DeltaR_{min}(#mu#mu)", "events/bin");
  h = new TH1D("s103", "#chi^{2} rank",                         10, 0., 10.);    h->Sumw2();    setTitles(h, "#Delta #chi^{2}", "events/bin");

  h = new TH1D("s104", "Nr. of Bs",                             10, 0., 10.);    h->Sumw2();    setTitles(h, "nr. of B", "events/bin");

  // - 2D
  h2 = new TH2D("S106", "Nr. of cand. vs. #sigma_{3D}",              10, 0., 10., 10, 0., 10.);  h2->Sumw2();
  setTitles2(h2, "#sigma_{3D}", "Nr. of cand.");

  h2 = new TH2D("S107", "#mu_{2} vs. #mu_{1} selection",              4, -2., 2., 4, -2., 2.);  h2->Sumw2();
  setTitles2(h2, "#mu_{1} is selected", "#mu_{2} is selected");

  h2 = new TH2D("S108", "#mu_{2} vs. #mu_{1} pT ranking",           10, 0., 10., 10, 0., 10.);  h2->Sumw2();
  setTitles2(h2, "#mu_{1} rank", "#mu_{2} rank");

  h2 = new TH2D("S109", "#mu_{2} vs. #mu_{1} leading #mu",            4, -2., 2., 4, -2., 2.);  h2->Sumw2();
  setTitles2(h2, "#mu_{1} is leading #mu", "#mu_{2} is leading #mu");

  h2 = new TH2D("S110", "#mu_{2} vs. #mu_{1} #DeltaR_{min}",          4, -2., 2., 4, -2., 2.);  h2->Sumw2();
  setTitles2(h2, "#mu_{1} is #DeltaR_{min}(#mu_{1}#mu)", "#mu_{2} is #DeltaR_{min}(#mu#mu_{2})");

  h2 = new TH2D("S111", "#DeltaR_{min} vs. leading #mu",              6, -3., 3., 6, -3., 3.);  h2->Sumw2();
  setTitles2(h2, "#mu_{1 or 2} is #DeltaR_{min}", "#mu_{1 or 2} is leading");


  // -- Tracking plots
  h = new TH1D("teff", "Track Cut efficiencies",                20, 0., 20.);    h->Sumw2();

  h = new TH1D("t100", "pt leptons",                            50, 0., 25.);    h->Sumw2();    setTitles(h, "p_{T}^{#mu}", "events/bin");
  h = new TH1D("t101", "pt l1",                                 50, 0., 25.);    h->Sumw2();    setTitles(h, "p_{T}^{#mu 1}", "events/bin");
  h = new TH1D("t102", "pt l2",                                 50, 0., 25.);    h->Sumw2();    setTitles(h, "p_{T}^{#mu 2}", "events/bin");             
  h = new TH1D("t103", "pt kaon",                               50, 0., 25.);    h->Sumw2();    setTitles(h, "p_{T}^{kaon}", "events/bin");

  h = new TH1D("t110", "eta leptons",                           50, -5., 5.);    h->Sumw2();    setTitles(h, "#eta^{#mu}", "events/bin");
  h = new TH1D("t111", "eta l1",                                50, -5., 5.);    h->Sumw2();    setTitles(h, "#eta^{#mu 1}", "events/bin");
  h = new TH1D("t112", "eta l2",                                50, -5., 5.);    h->Sumw2();    setTitles(h, "#eta^{#mu 2}", "events/bin");
  h = new TH1D("t113", "eta kaon",                              50, -5., 5.);    h->Sumw2();    setTitles(h, "#eta^{kaon}", "events/bin");

  h = new TH1D("t120", "tip leptons",                           50, 0., 0.5);    h->Sumw2();    setTitles(h, "tip^{#mu}_{xy}", "events/bin");   
  h = new TH1D("t121", "tip l1",                                50, 0., 0.5);    h->Sumw2();    setTitles(h, "tip^{#mu 1}_{xy}", "events/bin");
  h = new TH1D("t122", "tip l2",                                50, 0., 0.5);    h->Sumw2();    setTitles(h, "tip^{#mu 2}_{xy}", "events/bin");
  h = new TH1D("t123", "tip kaon",                              50, 0., 0.5);    h->Sumw2();    setTitles(h, "tip^{kaon}_{xy}", "events/bin");

  // -- Muon IDs of muon & kaon tracks
  h = new TH1D("mid", "Muon ID of the 2 muon tracks",        200, -20., 20.);    h->Sumw2();
  h = new TH1D("kid", "Muon ID of the kaon track",           200, -20., 20.);    h->Sumw2();


  // ============= Candidate plots ==================================================
  // -- isolation plots
  h = new TH1D("i100", "isolation veto",                       11, -1., 10.);    h->Sumw2();    setTitles(h, "isolation veto", "events/bin");
  h = new TH1D("i101", "isolation veto, fIsoVeto == 0",       50, -0.1, 10.);    h->Sumw2();    setTitles(h, "isolation veto", "events/bin");
  h = new TH1D("i102", "isolation veto, m > 4",               50, -0.1, 10.);    h->Sumw2();    setTitles(h, "isolation veto", "events/bin");
  h = new TH1D("i103", "isolation veto, m ~ J/psi",           50, -0.1, 10.);    h->Sumw2();    setTitles(h, "isolation veto", "events/bin");

  h = new TH1D("i104", "#Delta #phi(#mu#mu)",                  100, -5., 5.);    h->Sumw2();    setTitles(h,"#Delta #phi(#mu#mu)", "events/bin");
  h = new TH1D("i105", "#Delta #eta(#mu#mu)",                  100, -5., 5.);    h->Sumw2();    setTitles(h,"#Delta #eta(#mu#mu)", "events/bin");
  h = new TH1D("i106", "#Delta R(#mu#mu)",                       50, 0., 5.);    h->Sumw2();    setTitles(h,"#Delta R(#mu#mu)", "events/bin");


  h = new TH1D("i200", "#Delta #phi",                          100, -5., 5.);    h->Sumw2();
  h = new TH1D("i201", "p_{T, Trk} [GeV]",                     100, 0., 10.);    h->Sumw2();
  h = new TH1D("i202", "p_{T, Trk in cone} [GeV]",             100, 0., 10.);    h->Sumw2();

  // -- Jpsi
  h = new TH1D("i300", "#Delta #phi",                          100, -5., 5.);    h->Sumw2();
  h = new TH1D("i301", "p_{T, Trk} [GeV]",                     100, 0., 10.);    h->Sumw2();
  h = new TH1D("i302", "p_{T, Trk in cone} [GeV]",             100, 0., 10.);    h->Sumw2();

  // -- 2D-histograms
  h2 = new TH2D("I100", "#Delta R(#mu#mu) vs. fMass",          60, 0., 6., 50, 0., 10.);    h2->Sumw2();
  h2 = new TH2D("I101", "#Delta R(#mu#mu) vs. fMass",          60, 0., 6., 50, 0., 10.);    h2->Sumw2();
  h2 = new TH2D("I102", "drmin vs. cone size",                50, 0., 10., 50, 0., 10.);    h2->Sumw2();

  h2 = new TH2D("I104", "#Delta #phi(#mu#mu) vs. #Delta #eta(#mu#mu)", 100, -5., 5., 100, -5., 5.); h2->Sumw2();
  h2 = new TH2D("I105", "#Delta R(#mu#mu) vs. #Delta #phi(#mu#mu)", 50, 0., 5., 100, -5., 5.); h2->Sumw2();
  h2 = new TH2D("I106", "#Delta R(#mu#mu) vs. #Delta #eta(#mu#mu)", 50, 0., 5., 100, -5., 5.); h2->Sumw2();


  // -- secondary vertex plots
  h = new TH1D("v100", "chi2",                                 100, 0., 20.);    h->Sumw2();    setTitles(h, "#chi^{2}", "events/bin");
  h = new TH1D("v101", "prob(chi2, ndof)",                       40, 0., 1.);    h->Sumw2();    setTitles(h, "P(#chi^{2},ndof)", "events/bin");

  h = new TH1D("v110", "lxy",                                   50, 0., 2.0);    h->Sumw2();    setTitles(h, "l_{xy} [cm]", "events/bin");                                 
  h = new TH1D("v111", "sxy",                                  50, 0., 0.05);    h->Sumw2();    setTitles(h, "#sigma_{xy} [cm]", "events/bin");
  h = new TH1D("v112", "lxy/sxy",                               50, 0., 50.);    h->Sumw2();    setTitles(h, "l_{xy}/#sigma_{xy}", "events/bin");

  h = new TH1D("v120", "l3d",                                   50, 0., 2.0);    h->Sumw2();    setTitles(h, "l_{3d} [cm]", "events/bin");
  h = new TH1D("v121", "s3d",                                  50, 0., 0.05);    h->Sumw2();    setTitles(h, "#sigma_{3d} [cm]", "events/bin");
  h = new TH1D("v122", "l3d/s3d",                               50, 0., 50.);    h->Sumw2();    setTitles(h, "l_{3d}/#sigma_{3d}", "events/bin");

  // -- B candidate plots
  h = new TH1D("b100", "mass",                                 500, 0., 50.);    h->Sumw2();    setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b101", "mass",                                120, 4.8, 6.0);    h->Sumw2();    setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b102", "pT",                                    60, 0., 30.);    h->Sumw2();    setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("b103", "eta",                                   50, -5., 5.);    h->Sumw2();    setTitles(h, "#eta", "events/bin");
  h = new TH1D("b104", "proper decay time",                    100, 0., 0.5);    h->Sumw2();    setTitles(h, "#tau_{B} [...]", "events/bin");
  h = new TH1D("b110", "cos(angle)",                         100, 0.99, 1.0);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b120", "cos(angle), pT(B) > 5",               100, 0.99, 1.);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b130", "cos(angle), l3d/s3d > 4",             100, 0.99, 1.);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b131", "cos(angle), l3d/s3d > 10",            100, 0.99, 1.);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b140", "cos(angle), lxy/sxy > 4",             100, 0.99, 1.);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b141", "cos(angle), lxy/sxy > 10",            100, 0.99, 1.);    h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");

  h = new TH1D("b150", "deltaRmm",                              50, 0., 5.0);    h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");
  h = new TH1D("b200", "isolation",                              40, 0., 1.);    h->Sumw2();   setTitles(h, "Isolation", "events/bin");
  h = new TH1D("b201", "isolation veto",                        10, 0., 10.);    h->Sumw2();   setTitles(h, "Isolation veto", "events/bin");


  // -- Jpsi candidate plots
  h = new TH1D("j100", "mass",                                 500, 0., 50.);    h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("j101", "mass", 60, 2.8, 3.4); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("j102", "pT", 60, 0., 30.); h->Sumw2();   setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("j103", "eta", 50, -5., 5.); h->Sumw2();  setTitles(h, "#eta", "events/bin");

  h = new TH1D("j150", "deltaRmm", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");
  h = new TH1D("j151", "deltaRkj", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi K)", "events/bin");
  h = new TH1D("j152", "deltaRmj1", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{1})", "events/bin");
  h = new TH1D("j153", "deltaRmj2", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{2})", "events/bin");

  // -- Kaon candidate plots for B+ -> J/Psi K+
  // -- 1D-histograms
  h = new TH1D("keff", "number of truth-matched kaons", 100, 0., 100.);

  h = new TH1D("k100", "mass", 200, 0., 10.);
  h = new TH1D("k101", "#chi^{2}", 500, 0., 50.);
  h = new TH1D("k102", "p_{T}^{B}", 60, 0., 30.);
  h = new TH1D("k103", "#eta^{B}",  50, -5., 5.);
  h = new TH1D("k104", "p_{T}^{K}", 50, 0., 25.);
  h = new TH1D("k105", "#eta^{K}",  50, -5., 5.);

  // -- 2D-histograms
  h2 = new TH2D("K100", "#Delta m vs. gen. PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K101", "#Delta #chi^{2} vs. gen. PDG", 5000, 0., 5000., 500, -25., 25.);
  h2 = new TH2D("K102", "#Delta m vs. gen. mother PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K103", "#Delta #chi^{2} vs. gen. mother PDG", 5000, 0., 5000., 500, -25., 25.);

  h2 = new TH2D("K200", "#Delta m vs. gen. PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K201", "#Delta #chi^{2} vs. gen. PDG", 5000, 0., 5000., 500, -25., 25.);
  h2 = new TH2D("K202", "#Delta m vs. gen. mother PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K203", "#Delta #chi^{2} vs. gen. mother PDG", 5000, 0., 5000., 500, -25., 25.);


  // =================== residual plots ===============================================

  // -- track origins (mother, gran-mother)
  h2 = new TH2D("D100", "Mother of 2. muon vs. Mother of 1. muon", 5000, 0., 5000., 5000, 0., 5000.);  h2->Sumw2();
  h2 = new TH2D("D101", "Grand-Mother of 2. muon vs. Grand-Mother of 1. muon", 5000, 0., 5000., 5000, 0., 5000.);  h2->Sumw2();

  h = new TH1D("d100", "Mother of 1. muon", 5000, -1., 4999.);  h->Sumw2();
  h = new TH1D("d101", "Mother of 2. muon", 5000, -1., 4999.);  h->Sumw2();
  h = new TH1D("d102", "Mother of kaon", 5000, -1., 4999.);  h->Sumw2();
  h = new TH1D("d103", "Grand-Mother of 1. muon", 5000, -1., 4999.);  h->Sumw2();
  h = new TH1D("d104", "Grand-Mother of 2. muon", 5000, -1., 4999.);  h->Sumw2();
  h = new TH1D("d105", "Grand-Mother of kaon", 5000, -1., 4999.);  h->Sumw2();  
 
  // -- pt resolution 
  // only global muons
  h2 = new TH2D("R100", "#sigma_{p_{T}, #mu} vs p_{T}", 25, 0., 25., 100, -1., 1.);
  h2 = new TH2D("R101", "#sigma_{p_{T}, #mu}/p_{T} vs p_{T}", 25, 0., 25., 100, -1., 1.);
  h2 = new TH2D("R102", "#sigma_{p_{T}, #mu} vs #eta", 25, 0., 2.5, 100, -1., 1.);
  h2 = new TH2D("R103", "#sigma_{p_{T}, #mu}/p_{T}  vs #eta", 25, 0., 2.5, 100, -1., 1.);

  // all tracks
  h2 = new TH2D("R200", "#sigma_{p_{T}, tracks} vs p_{T}", 25, 0., 25., 100, -1., 10.);
  h2 = new TH2D("R201", "#sigma_{p_{T}, tracks}/p_{T} vs p_{T}", 25, 0., 25., 100, -1., 10.);
  h2 = new TH2D("R202", "#sigma_{p_{T}, tracks} vs #eta", 25, 0., 2.5, 100, -1., 10.);
  h2 = new TH2D("R203", "#sigma_{p_{T}, tracks}/p_{T}  vs #eta", 25, 0., 2.5, 100, -1., 10.);

  // -- Profiles
  // only global muons
  p = new TProfile("r100","#sigma_{p_{T},#mu} vs p_{T}", 25, 0., 25., "s");
  p = new TProfile("r101","#sigma_{p_{T}, #mu}/p_{T} vs p_{T}", 25, 0., 25., "s");
  p = new TProfile("r102","#sigma_{p_{T}, #mu} vs #eta", 25, 0., 2.5, "s");
  p = new TProfile("r103","#sigma_{p_{T}, #mu}/p_{T}  vs #eta", 25, 0., 2.5, "s");

  h  = new TH1D("r100C", "Counter for profile r100, r101", 25, 0., 25.);  h->Sumw2();
  h  = new TH1D("r102C", "Counter for profile r102, r103", 25, 0., 2.5);  h->Sumw2();

  // all tracks
  p = new TProfile("r200","#sigma_{p_{T}, tracks} vs p_{T}", 25, 0., 25., "s");
  p = new TProfile("r201","#sigma_{p_{T}, tracks}/p_{T} vs p_{T}", 25, 0., 25., "s");
  p = new TProfile("r202","#sigma_{p_{T}, tracks} vs #eta", 25, 0., 2.5, "s");
  p = new TProfile("r203","#sigma_{p_{T}, tracks}/p_{T}  vs #eta", 25, 0., 2.5, "s");

  h  = new TH1D("r200C", "Counter for profile r200, r201", 25, 0., 25.);  h->Sumw2();
  h  = new TH1D("r202C", "Counter for profile r202, r203", 25, 0., 2.5);  h->Sumw2();

  // -- sec. vertex residuals
  h = new TH1D("v200", "Decay length resolution (2D)", 500, -0.05, 0.05); h->Sumw2(); setTitles(h, "l_{xy}^{rec} - l_{xy}^{sim} [cm]", "events/bin");
  h = new TH1D("v201", "Decay length resolution (3D)", 500, -0.1, 0.1); h->Sumw2(); setTitles(h, "l_{3D}^{rec} - l_{3D}^{sim} [cm]", "events/bin");
  // for Norm change to:
//   h = new TH1D("v201", "Decay length resolution (3D)", 500, -0.01, 0.01); h->Sumw2(); setTitles(h, "l_{3D}^{rec} - l_{3D}^{sim} [cm]", "events/bin");
  h = new TH1D("v202", "Proper decay time resolution (2D)", 500, -0.05, 0.05); h->Sumw2(); setTitles(h, "#tau_{xy}^{rec} - #tau_{xy}^{sim} [cm/c]", "events/bin");
  h = new TH1D("v203", "Proper decay time resolution (3D)", 500, -0.05, 0.05); h->Sumw2(); setTitles(h, "#tau_{3D}^{rec} - #tau_{3D}^{sim} [cm/c]", "events/bin");
  h = new TH1D("v204", "Proper decay time resolution (2D)", 500, -1000., 1000.); h->Sumw2(); setTitles(h, "#tau_{xy}^{rec} - #tau_{xy}^{sim} [fs]", "events/bin");
  h = new TH1D("v205", "Proper decay time resolution (3D)", 500, -1000., 1000.); h->Sumw2(); setTitles(h, "#tau_{3D}^{rec} - #tau_{3D}^{sim} [fs]", "events/bin");

  h = new TH1D("v300", "Primary Vtx resolution (2D)", 100, 0., 0.02); h->Sumw2(); setTitles(h, "PV_{xy}^{rec} - PV_{xy}^{sim} [cm]", "events/bin");
  h = new TH1D("v301", "Primary Vtx resolution (3D)", 100, 0., 0.02); h->Sumw2(); setTitles(h, "PV_{3D}^{rec} - PV_{3D}^{sim} [cm]", "events/bin");
  h = new TH1D("v400", "Secondary Vtx resolution (2D)", 100, 0., 0.05); h->Sumw2(); setTitles(h, "SV_{xy}^{rec} - SV_{xy}^{sim} [cm]", "events/bin");
  h = new TH1D("v401", "Secondary Vtx resolution (3D)", 100, 0., 0.05); h->Sumw2(); setTitles(h, "SV_{3D}^{rec} - SV_{3D}^{sim} [cm]", "events/bin");

  // =================== Muon efficieny / fake rate plots =================================

  fMisID = new TH1D("MisID", "Muon mis-ID Rate", 50, 0., 50.);      fMisID->Sumw2();

  // 1D-histograms
  //---------------

  // -- Counter histograms
  h = new TH1D("m700", "N_{#mu} in decay B_{s}^{0} #rightarrow #mu #mu)", 100, 0., 100.);  h->Sumw2();
  h = new TH1D("m701", "N_{#mu} in decay B #rightarrow J/#psi(#rightarrow mu mu) K", 100, 0., 100.);  h->Sumw2();
  h = new TH1D("m702", "N_{#mu} in decay J/#psi #rightarrow mu mu", 100, 0., 100.);  h->Sumw2();

  // -- pions
  h = new TH1D("m000", "PDG id all tracks", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("m010", "number of associated tracks per muons", 100, 0., 100.);  h->Sumw2();

  h = new TH1D("m100", "p_{T} (all #pi) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m101", "p_{T} (#pi = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m102", "p_{T} (#pi != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m110", "#eta (all #pi)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m111", "#eta (#pi = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m112", "#eta (#pi != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- kaons
  h = new TH1D("m200", "p_{T} (all K) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m201", "p_{T} (K = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m202", "p_{T} (K != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m210", "#eta (all K)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m211", "#eta (K = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m212", "#eta (K != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- protons
  h = new TH1D("m300", "p_{T} (all p) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m301", "p_{T} (p = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m302", "p_{T} (p != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m310", "#eta (all p)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m311", "#eta (p = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m312", "#eta (p != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- global muons
  h = new TH1D("m400", "p_{T} (all #mu, global) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m401", "p_{T} (#mu = #mu, global) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m402", "p_{T} (#mu != #mu, global) [GeV]", 50, 0., 25.);  h->Sumw2();setTitles(h, "p_{T} [GeV]", "events/bin");

  h = new TH1D("m410", "#eta (all #mu, global)", 50, -5., 5.);  h->Sumw2();  setTitles(h, "#eta", "events/bin");
  h = new TH1D("m411", "#eta (#mu = #mu, global)", 50, -5., 5.);  h->Sumw2();  setTitles(h, "#eta", "events/bin");
  h = new TH1D("m412", "#eta (#mu != #mu, global)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");

  // -- tracker muons
  h = new TH1D("m500", "p_{T} (all #mu, tracker) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m501", "p_{T} (#mu = #mu, tracker) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m502", "p_{T} (#mu != #mu, tracker) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");

  h = new TH1D("m510", "#eta (all #mu, tracker)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");
  h = new TH1D("m511", "#eta (#mu = #mu, tracker)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");
  h = new TH1D("m512", "#eta (#mu != #mu, tracker)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");

  // -- l1 muons
  h = new TH1D("m600", "p_{T} (all #mu, l1) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m601", "p_{T} (#mu = #mu, l1) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("m602", "p_{T} (#mu != #mu, l1) [GeV]", 50, 0., 25.);  h->Sumw2(); setTitles(h, "p_{T} [GeV]", "events/bin");

  h = new TH1D("m610", "#eta (all #mu, l1)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");
  h = new TH1D("m611", "#eta (#mu = #mu, l1)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");
  h = new TH1D("m612", "#eta (#mu != #mu, l1)", 50, -5., 5.);  h->Sumw2(); setTitles(h, "#eta", "events/bin");


  // 2D-histograms
  //--------------- 

  // -- pions
  h2 = new TH2D("M100", "p_{T} vs. #eta (all #pi) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M101", "p_{T}  vs. #eta (#pi = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M102", "p_{T} vs. #eta  (#pi != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- kaons
  h2 = new TH2D("M200", "p_{T} vs. #eta (all K) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M201", "p_{T}  vs. #eta (K = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M202", "p_{T} vs. #eta  (K != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- protons
  h2 = new TH2D("M300", "p_{T} vs. #eta (all protons) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M301", "p_{T}  vs. #eta (p = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M302", "p_{T} vs. #eta  (p != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- global muons
  h2 = new TH2D("M400", "p_{T} vs. #eta (all #mu, global) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M401", "p_{T}  vs. #eta (#mu = #mu, global) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M402", "p_{T} vs. #eta  (#mu != #mu, global) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- tracker muons
  h2 = new TH2D("M500", "p_{T} vs. #eta (all #mu, tracker) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M501", "p_{T}  vs. #eta (#mu = #mu, tracker) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M502", "p_{T} vs. #eta  (#mu != #mu, tracker) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- l1 muons
  h2 = new TH2D("M600", "p_{T} vs. #eta (all #mu, l1) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M601", "p_{T}  vs. #eta (#mu = #mu, l1) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M602", "p_{T} vs. #eta  (#mu != #mu, l1) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  if (fDebug & 2) { cout << "bookHist> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::book(int offset) {

  if (fDebug & 2) { cout << "book> Start" << endl; }

  TH1D *h;
  h = new TH1D(Form("c%d00", offset), "pT(B)",         60, 0., 30.);    h->Sumw2();     setTitles(h,  "p_{T, B} [GeV]", "events/bin");
  h = new TH1D(Form("c%d01", offset), "eta(B)",        50, -5., 5.);    h->Sumw2();     setTitles(h,  "#eta_{B}", "events/bin");

  h = new TH1D(Form("c%d10", offset), "leading pT leptons", 50, 0., 25.); h->Sumw2(); setTitles(h, "p_{T, #mu}^{max} [GeV]", "events/bin");
  h = new TH1D(Form("c%d11", offset), "eta leptons",    50, -5., 5.);   h->Sumw2();     setTitles(h, "#eta_{#mu}", "events/bin");
  h = new TH1D(Form("c%d12", offset), "pT leptons",     50, 0., 25.);   h->Sumw2();     setTitles(h, "p_{T, #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d13", offset), "TIP leptons",    50, 0., 0.2);   h->Sumw2();     setTitles(h, "TIP_{#mu} [cm]", "events/bin");

  h = new TH1D(Form("c%d14", offset), "s3d",            100, 0., 0.2);   h->Sumw2();     setTitles(h, "#sigma_{3D} [cm]", "events/bin");
  h = new TH1D(Form("c%d15", offset), "sxy",            100, 0., 0.1);  h->Sumw2();     setTitles(h, "#sigma_{xy} [cm]", "events/bin");
  h = new TH1D(Form("c%d16", offset), "l3d/s3d",         50, 0., 50.);  h->Sumw2();     setTitles(h, "l_{3D}/#sigma_{3D}", "events/bin");
  h = new TH1D(Form("c%d18", offset), "delta phi",      100, -5., 5.);  h->Sumw2();     setTitles(h, "#Delta #phi(#mu#mu)", "events/bin");
  h = new TH1D(Form("c%d19", offset), "delta eta",      100, -5., 5.);  h->Sumw2();     setTitles(h, "#Delta #eta(#mu#mu)", "events/bin");

  h = new TH1D(Form("c%d20", offset), "delta rmm",       50, 0., 5.0);  h->Sumw2();     setTitles(h, "#Delta R(#mu#mu)", "events/bin");
  h = new TH1D(Form("c%d21", offset), "cos(angle) 2D",    200,0.95,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{xy}", "events/bin");
  h = new TH1D(Form("c%d22", offset), "lxy/sxy",        50, 0., 50.);   h->Sumw2();     setTitles(h, "l_{xy}/#sigma_{xy}", "events/bin");
  h = new TH1D(Form("c%d23", offset), "lxy",           100, 0., 0.4);   h->Sumw2();     setTitles(h, "l_{xy} [cm]", "events/bin");
  h = new TH1D(Form("c%d24", offset), "isolation veto", 20, 0., 20.);   h->Sumw2();     setTitles(h, "I_{V}", "events/bin");
  h = new TH1D(Form("c%d25", offset), "isolation",      55, 0., 1.1);   h->Sumw2();     setTitles(h, "I", "events/bin");
  h = new TH1D(Form("c%d26", offset), "l3d",           100, 0., 0.4);   h->Sumw2();     setTitles(h, "l_{3D} [cm]", "events/bin");
  h = new TH1D(Form("c%d27", offset), "chi2",           50, 0., 5.0);   h->Sumw2();     setTitles(h, "#chi^{2}", "events/bin");
  h = new TH1D(Form("c%d28", offset), "cos(angle) 3D",    200,0.95,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{3D}", "events/bin");
  h = new TH1D(Form("c%d29", offset), "process",       52, 0., 52.);    h->Sumw2();     setTitles(h, "Process (GGF,FEX,GSP for t,b,c) ", "events/bin");

  h = new TH1D(Form("c%d30", offset), "mass",           120, 4.8, 6.0); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d31", offset), "mass",           500, 0., 50.); h->Sumw2();      setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d32", offset), "mass",           200, 4., 6.0); h->Sumw2();      setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d33", offset), "mass",           160, 2., 10.); h->Sumw2();      setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");

  h = new TH1D(Form("c%d35", offset), "NTRK",          100, 0., 200.);    h->Sumw2();   setTitles(h, "N_{Track}", "events/bin");

  h = new TH1D(Form("c%d36", offset), "prob(chi2, ndof)",    50, 0., 1.);   h->Sumw2();     setTitles(h, "P(#chi^{2}, ndof)", "events/bin");
  h = new TH1D(Form("c%d37", offset), "chi2/ndof",    50, 0., 5.0);   h->Sumw2();     setTitles(h, "#chi^{2}/ndof", "events/bin");

  h = new TH1D(Form("c%d40", offset), "cos(angle) 2D",     200,-1.0,1.0);    h->Sumw2();   setTitles(h, "cos #alpha_{xy}", "events/bin");
  h = new TH1D(Form("c%d41", offset), "cos(angle) 2D",     200,0.95,1.0);    h->Sumw2();   setTitles(h, "cos #alpha_{xy}", "events/bin"); 
  h = new TH1D(Form("c%d42", offset), "cos(angle) 2D",      50,0.97,1.0);    h->Sumw2();   setTitles(h, "cos #alpha_{xy}", "events/bin"); 
  h = new TH1D(Form("c%d43", offset), "cos(angle) 2D",     200,0.99,1.0);    h->Sumw2();   setTitles(h, "cos #alpha_{xy}", "events/bin");

  h = new TH1D(Form("c%d50", offset), "flight length (if sim. Vtx exists)", 100,  0.0, 0.1);   h->Sumw2();   setTitles(h, "t_{rec}", "events/bin");
  h = new TH1D(Form("c%d51", offset), "flight length (if sim. Vtx exists)", 100,  0.0, 0.1);   h->Sumw2();   setTitles(h, "t_{gen}", "events/bin");
  h = new TH1D(Form("c%d52", offset), "flight length (if sim. Vtx exists)", 100, -0.02, 0.02); h->Sumw2();   setTitles(h, "t_{rec} - t_{gen}", "events/bin");


  h = new TH1D(Form("c%d60", offset), "cos(angle) 3D",    200,-1.0,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{3D}", "events/bin");
  h = new TH1D(Form("c%d61", offset), "cos(angle) 3D",    200,0.95,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{3D}", "events/bin");
  h = new TH1D(Form("c%d62", offset), "cos(angle) 3D",    50,0.97,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{3D}", "events/bin");
  h = new TH1D(Form("c%d63", offset), "cos(angle) 3D",    200,0.99,1.0);   h->Sumw2();     setTitles(h, "cos #alpha_{3D}", "events/bin");

  h = new TH1D(Form("c%d70", offset), "proper decay length 3D", 100,  0.0, 0.5);   h->Sumw2();   setTitles(h, "ct_{3D}", "events/bin");
  h = new TH1D(Form("c%d71", offset), "proper decay length 2D", 100,  0.0, 0.5);   h->Sumw2();   setTitles(h, "ct_{xy}", "events/bin");
  h = new TH1D(Form("c%d72", offset), "proper decay length 2D (cos)", 100,  0.0, 0.5);   h->Sumw2();   setTitles(h, "ct_{xy}", "events/bin");

  h = new TH1D(Form("c%d73", offset), "proper decay length 3D", 200,  0.0, 0.2);   h->Sumw2();   setTitles(h, "ct_{3D}", "events/bin");
  h = new TH1D(Form("c%d74", offset), "proper decay length 2D", 200,  0.0, 0.2);   h->Sumw2();   setTitles(h, "ct_{xy}", "events/bin");
  h = new TH1D(Form("c%d75", offset), "proper decay length 2D (cos)", 200,  0.0, 0.2);   h->Sumw2();   setTitles(h, "ct_{xy}", "events/bin");


  if (fDebug & 2) { cout << "book> End" << endl; }
}

// ----------------------------------------------------------------------

void treeBmm::book2(int offset) {

  if (fDebug & 2) { cout << "book2> Start" << endl; }

  TH2D *h;

  h = new TH2D(Form("C%d00", offset), "dpT vs. sxy", 50, 0., 0.05, 50, 0., 0.05); h->Sumw2();  
  setTitles2(h, "#sigma_{pT}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d01", offset), "pT vs. eta",   50, 0., 25., 50, -5., 5); h->Sumw2(); 
  setTitles2(h, "p_{T, #mu} [GeV]", "#eta_{#mu} [cm]");

  h = new TH2D(Form("C%d02", offset), "pT vs. sxy", 50, 0., 25., 50, 0., 0.05); h->Sumw2(); 
  setTitles2(h, "p_{T,#mu}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d03", offset), "pT vs. lxy/sxy", 50, 0., 25., 50, 0., 50.); h->Sumw2();  
  setTitles2(h, "p_{T,#mu}", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d04", offset), "pT vs. TIP", 50, 0., 25., 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "TIP [mm]");

  h = new TH2D(Form("C%d05", offset), "pT vs. mass", 50, 0., 25., 500, 0., 50.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "m_{#mu #mu} [GeV]");

  h = new TH2D(Form("C%d06", offset), "pT vs. deltaR", 50, 0., 25., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d07", offset), "pT vs. cos(angle)", 50, 0., 25., 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d08", offset), "pT vs. isolation veto", 50, 0., 25., 10, 0., 10.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "Isolation veto");

  h = new TH2D(Form("C%d09", offset), "pT vs. chi2", 50, 0., 25., 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#chi^{2}");

  h = new TH2D(Form("C%d10", offset), "pT vs. pT(B)", 50, 0., 25., 60, 0., 30.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d11", offset), "pT vs. eta(B)", 50, 0., 25., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#eta_{B}");

  h = new TH2D(Form("C%d12", offset), "eta vs. sxy", 50, -5., 5, 50, 0., 0.05); h->Sumw2(); 
  setTitles2(h, "#eta_{#mu}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d13", offset), "eta vs. lxy/sxy", 50, -5., 5,  50, 0., 50.); h->Sumw2();  
  setTitles2(h, "#eta_{#mu}", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d14", offset), "eta vs. TIP", 50, -5., 5, 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "TIP [mm]");

  h = new TH2D(Form("C%d15", offset), "eta vs. mass", 50, -5., 5, 500, 0., 50.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "m_{#mu #mu} [GeV]");

  h = new TH2D(Form("C%d16", offset), "eta vs. deltaR", 50, -5., 5, 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d17", offset), "eta vs.cos(angle)", 50, -5., 5, 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d18", offset), "eta vs. isolation veto", 50, -5., 5, 10, 0., 10.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "Isolation veto");

  h = new TH2D(Form("C%d19", offset), "eta vs. chi2", 50, -5., 5, 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#chi^{2}");

  h = new TH2D(Form("C%d20", offset), "eta vs. pT(B)", 50, -5., 5, 60, 0., 30.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d21", offset), "eta vs. eta(B)", 50, -5., 5, 50, -5., 5.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#eta_{B}");

  h = new TH2D(Form("C%d22", offset), "mass vs. sxy", 500, 0., 50, 50, 0., 0.05); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d23", offset), "mass vs. lxy/sxy", 500, 0., 50, 50, 0., 50.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d24", offset), "mass vs. TIP", 500, 0., 50, 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "TIP [mm]");

  h = new TH2D(Form("C%d25", offset), "mass vs. deltaR", 500, 0., 50., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d26", offset), "mass vs. cos(angle)", 500, 0., 50, 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d27", offset), "mass vs. isolation veto", 500, 0., 50, 10, 0., 10.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "Isolation veto");

  h = new TH2D(Form("C%d28", offset), "mass vs. chi2", 500, 0., 50,  100, 0., 5.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#chi^{2}");

  h = new TH2D(Form("C%d29", offset), "mass vs. pT(B)", 500, 0., 50.,  60, 0., 30.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d30", offset), "mass vs. eta(B)", 500, 0., 50., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#eta_{B}");

  h = new TH2D(Form("C%d31", offset), "pT(B) vs. sxy",            60, 0., 30., 50, 0., 0.05); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#sigma_{xy} [cm]");
  
  h = new TH2D(Form("C%d32", offset), "pT(B) vs. lxy/sxy",        60, 0., 30., 50, 0., 50.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d33", offset), "pT(B) vs. TIP",            60, 0., 30., 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "TIP [mm]");

  h = new TH2D(Form("C%d34", offset), "pT(B) vs. deltaR",         60, 0., 30., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d35", offset), "pT(B) vs. cos(angle)",     60, 0., 30., 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "cos #alpha(p,v) ");

  h = new TH2D(Form("C%d36", offset), "pT(B) vs. isolation veto", 60, 0., 30., 10, 0., 10.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "Isolation veto");

  h = new TH2D(Form("C%d37", offset), "pT(B) vs. chi2",           60, 0., 30., 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#chi^{2}");

  h = new TH2D(Form("C%d38", offset), "pT(B) vs. eta(B)",         60, 0., 30., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#eta_{B}");

  h = new TH2D(Form("C%d40", offset), "delta phi vs. delta eta", 100, -5., 5., 100, -5., 5.); h->Sumw2();
  setTitles2(h, "#Delta #eta(#mu#mu)", "#Delta #phi(#mu#mu)");

  h = new TH2D(Form("C%d41", offset), "delta rmm vs. delta phi", 100, -5., 5., 50, 0., 5.); h->Sumw2();
  setTitles2(h, "#Delta #phi(#mu#mu)", "#Delta R(#mu#mu)");

  h = new TH2D(Form("C%d42", offset), "delta rmm vs. delta eta", 100, -5., 5., 50, 0., 5.); h->Sumw2();
  setTitles2(h, "#Delta #eta(#mu#mu)", "#Delta R(#mu#mu)");


  if (fDebug & 2) { cout << "book2> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::histogram(int offset) {

  if (fDebug & 2) { cout << "histogram> Start" << endl; }

  ((TH1D*)gDirectory->Get(Form("c%d00", offset)))->Fill(fPt);
  ((TH1D*)gDirectory->Get(Form("c%d01", offset)))->Fill(fEta);
  
  ((TH1D*)gDirectory->Get(Form("c%d10", offset)))->Fill(fPtL0);
  ((TH1D*)gDirectory->Get(Form("c%d11", offset)))->Fill(fEtaL0);
  ((TH1D*)gDirectory->Get(Form("c%d11", offset)))->Fill(fEtaL1);
  ((TH1D*)gDirectory->Get(Form("c%d12", offset)))->Fill(fPtL0);
  ((TH1D*)gDirectory->Get(Form("c%d12", offset)))->Fill(fPtL1);
  ((TH1D*)gDirectory->Get(Form("c%d13", offset)))->Fill(fTipL0);
  ((TH1D*)gDirectory->Get(Form("c%d13", offset)))->Fill(fTipL1);

  ((TH1D*)gDirectory->Get(Form("c%d14", offset)))->Fill(fS3d);
  ((TH1D*)gDirectory->Get(Form("c%d15", offset)))->Fill(fSxy);
  ((TH1D*)gDirectory->Get(Form("c%d16", offset)))->Fill(fL3d/fS3d);

  ((TH1D*)gDirectory->Get(Form("c%d18", offset)))->Fill(fDphi);
  ((TH1D*)gDirectory->Get(Form("c%d19", offset)))->Fill(fDeta);

  ((TH1D*)gDirectory->Get(Form("c%d20", offset)))->Fill(fRMM);
  ((TH1D*)gDirectory->Get(Form("c%d21", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d22", offset)))->Fill(fLxy/fSxy);
  ((TH1D*)gDirectory->Get(Form("c%d23", offset)))->Fill(fLxy);
  ((TH1D*)gDirectory->Get(Form("c%d24", offset)))->Fill(fIsoVeto);
  ((TH1D*)gDirectory->Get(Form("c%d25", offset)))->Fill(fIso);
  ((TH1D*)gDirectory->Get(Form("c%d26", offset)))->Fill(fL3d);
  ((TH1D*)gDirectory->Get(Form("c%d27", offset)))->Fill(fChi2);
  ((TH1D*)gDirectory->Get(Form("c%d28", offset)))->Fill(fCosAngle3);
  ((TH1D*)gDirectory->Get(Form("c%d29", offset)))->Fill(fProcessType);
  ((TH1D*)gDirectory->Get(Form("c%d30", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d31", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d32", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d33", offset)))->Fill(fMass);

  ((TH1D*)gDirectory->Get(Form("c%d35", offset)))->Fill(fpEvt->nRecTracks());

  ((TH1D*)gDirectory->Get(Form("c%d36", offset)))->Fill(fProb);
  ((TH1D*)gDirectory->Get(Form("c%d37", offset)))->Fill(fChi2/fNdof);

  ((TH1D*)gDirectory->Get(Form("c%d40", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d41", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d42", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d43", offset)))->Fill(fCosAngle);

  ((TH1D*)gDirectory->Get(Form("c%d50", offset)))->Fill(fFltR);
  ((TH1D*)gDirectory->Get(Form("c%d51", offset)))->Fill(fFltG);
  ((TH1D*)gDirectory->Get(Form("c%d52", offset)))->Fill(fFltRes);

  ((TH1D*)gDirectory->Get(Form("c%d60", offset)))->Fill(fCosAngle3);
  ((TH1D*)gDirectory->Get(Form("c%d61", offset)))->Fill(fCosAngle3);
  ((TH1D*)gDirectory->Get(Form("c%d62", offset)))->Fill(fCosAngle3);
  ((TH1D*)gDirectory->Get(Form("c%d63", offset)))->Fill(fCosAngle3);

  ((TH1D*)gDirectory->Get(Form("c%d70", offset)))->Fill(fTau);
  ((TH1D*)gDirectory->Get(Form("c%d71", offset)))->Fill(fTxy);
  ((TH1D*)gDirectory->Get(Form("c%d72", offset)))->Fill(fTxy*fCosAngle);

  ((TH1D*)gDirectory->Get(Form("c%d73", offset)))->Fill(fTau);
  ((TH1D*)gDirectory->Get(Form("c%d74", offset)))->Fill(fTxy);
  ((TH1D*)gDirectory->Get(Form("c%d75", offset)))->Fill(fTxy*fCosAngle);

  ((TH2D*)gDirectory->Get(Form("C%d00", offset)))->Fill(fDptL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d01", offset)))->Fill(fPtL0, fEtaL0);
  ((TH2D*)gDirectory->Get(Form("C%d01", offset)))->Fill(fPtL1, fEtaL1);
  ((TH2D*)gDirectory->Get(Form("C%d02", offset)))->Fill(fPtL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d03", offset)))->Fill(fPtL0, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d04", offset)))->Fill(fPtL0, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d04", offset)))->Fill(fPtL1, fTipL1);
  ((TH2D*)gDirectory->Get(Form("C%d05", offset)))->Fill(fPtL0, fMass);
  ((TH2D*)gDirectory->Get(Form("C%d06", offset)))->Fill(fPtL0, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d07", offset)))->Fill(fPtL0, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d08", offset)))->Fill(fPtL0, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d09", offset)))->Fill(fPtL0, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d10", offset)))->Fill(fPtL0, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d11", offset)))->Fill(fPtL0, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d12", offset)))->Fill(fEtaL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d13", offset)))->Fill(fEtaL0, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d14", offset)))->Fill(fEtaL0, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d15", offset)))->Fill(fEtaL0, fMass);
  ((TH2D*)gDirectory->Get(Form("C%d16", offset)))->Fill(fEtaL0, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d17", offset)))->Fill(fEtaL0, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d18", offset)))->Fill(fEtaL0, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d19", offset)))->Fill(fEtaL0, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d20", offset)))->Fill(fEtaL0, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d21", offset)))->Fill(fEtaL0, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d22", offset)))->Fill(fMass, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d23", offset)))->Fill(fMass, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d24", offset)))->Fill(fMass, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d25", offset)))->Fill(fMass, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d26", offset)))->Fill(fMass, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d27", offset)))->Fill(fMass, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d28", offset)))->Fill(fMass, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d29", offset)))->Fill(fMass, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d30", offset)))->Fill(fMass, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d31", offset)))->Fill(fPt, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d32", offset)))->Fill(fPt, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d33", offset)))->Fill(fPt, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d34", offset)))->Fill(fPt, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d35", offset)))->Fill(fPt, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d36", offset)))->Fill(fPt, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d37", offset)))->Fill(fPt, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d38", offset)))->Fill(fPt, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d40", offset)))->Fill(fDeta, fDphi);
  ((TH2D*)gDirectory->Get(Form("C%d41", offset)))->Fill(fDphi, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d42", offset)))->Fill(fDeta, fRMM);

  if (fDebug & 2) { cout << "histogram> End" << endl; }
}

// **************************************************************************************************
// --------------------------------------------------------------------------------------------------
void treeBmm::processDiscrimination() {

  if (fDebug & 2) { cout << "processDiscrimination> Start" << endl; }

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

    ((TH1D*)gDirectory->Get("g104"))->Fill(pG->fStatus);
    
    aid = TMath::Abs(pG->fID); 
    if ( aid == 1 || aid == 2 ||
 	 aid == 3 || aid == 4 || 
 	 aid == 5 || aid == 6 || 
 	 aid == 21) {
      if ( pG->fStatus == 3 ) {
	// 	cout << "quark/gluon from documentation #" << i << "(ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 2 &&  TMath::Abs(pG->fID) != 21) {
	// 	cout << "decayed quark/gluon #" << i << " (ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 1 ) {
	// 	cout << "undecayed (?) quark/gluon #" << i << " (ID: " << pG->fID  << ")" << endl;
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
  // GGF heavy flavour
  if (docPartCnt[5] == 1 && docAntiCnt[5] == 1) {
    fProcessType = 50; 
    //    printf("====> t: GGF (%i)\n", fProcessType);
    return;
  }
  
  if (docPartCnt[4] == 1 && docAntiCnt[4] == 1) {
    fProcessType = 40; 
    //    printf("====> b: GGF (%i)\n", fProcessType);
    return;
  } 
  
  if (docPartCnt[3] == 1 && docAntiCnt[3] == 1) {
    fProcessType = 30; 
    //    printf("====> c: GGF (%i)\n", fProcessType);
    return;
  }
  
  // FEX heavy flavour
  if ((docPartCnt[5] >= 1 && docAntiCnt[5] == 0) || (docPartCnt[5] == 0 && docAntiCnt[5] >= 1) ) {
    fProcessType = 51; 
    //    printf("====> t: FEX (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[4] >= 1 && docAntiCnt[4] == 0) || (docPartCnt[4] == 0 && docAntiCnt[4] >= 1) ) {
    fProcessType = 41; 
    //    printf("====> b: FEX (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[3] >= 1 && docAntiCnt[3] == 0) || (docPartCnt[3] == 0 && docAntiCnt[3] >= 1) ) {
    fProcessType = 31; 
    //    printf("====> c: FEX (%i)\n", fProcessType);
    return;
  }
  
  // GSP heavy flavour
  if (docPartCnt[5] == 0 && docAntiCnt[5] == 0 && (parPartCnt[5] >= 1 || parAntiCnt[5] >= 1)) {
    fProcessType = 52;
    //    printf("====> t: GSP (%i)\n", fProcessType); 
    return;
  }
  
  if (docPartCnt[4] == 0 && docAntiCnt[4] == 0 && (parPartCnt[4] >= 1 || parAntiCnt[4] >= 1)) {
    fProcessType = 42; 
    //    printf("====> b: GSP (%i)\n", fProcessType);
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



  if (fDebug & 2) { cout << "processDiscrimination> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::muonEfficiency() { 

  if (fDebug & 2) { cout << "muonEfficiency> Start" << endl; }

  double recMu(-1.), trkMu(-1.), l1Mu(-1.);
  int mcID(-1), genIdx(-1), momIdx(-1), gmoIdx(-1);
  double pT(0.), eta(0.);
  int mu_per_track(1);

  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {

    recMu  = fpEvt->getRecTrack(it)->fMuID;
    trkMu  = fpEvt->getRecTrack(it)->fKaID;
    l1Mu  = fpEvt->getRecTrack(it)->fElID;

    mcID   = fpEvt->getRecTrack(it)->fMCID;
    genIdx = fpEvt->getRecTrack(it)->fGenIndex;

    pT  = fpEvt->getRecTrack(it)->fPlab.Pt();
    eta = fpEvt->getRecTrack(it)->fPlab.Eta();
    
    ((TH1D*)gDirectory->Get("m000"))->Fill(mcID);

    if ( recMu > -0.5 ) {
      
      fnMu++;
    }

    if ( TMath::Abs(mcID) == 13 ) {
      
      frMu++;
    }

    if ( TMath::Abs(mcID) == 321 ) {
      
      frK++;
    }

    // -- number of tracks associated to same muon
    if ( TMath::Abs(mcID) == 13 ) {
      mu_per_track = 1;
      for (int jt = it + 1; jt < fpEvt->nRecTracks(); ++jt) {
	if  ( (genIdx > 1) 
	      && (genIdx == fpEvt->getRecTrack(jt)->fGenIndex) ) {
	  mu_per_track++;
	}
      }

      ((TH1D*)gDirectory->Get("m010"))->Fill(mu_per_track); 
    }
    
    // -- Muon statistics for sign. & norm. channel
    if ( TMath::Abs(mcID) == 13 ) {

      if ( genIdx > 1 ) {
	
	momIdx = fpEvt->getGenCand(genIdx)->fMom1;

	if ( momIdx >= 0 ) {

	  gmoIdx = fpEvt->getGenCand(momIdx)->fMom1;
	}
      }

      if (!fNorm && genPartType(momIdx) == 531) {
	
	frSigMu++;

	if (recMu > -0.5) {
	  
	  fnSigMu++;
	}

	if ( it == fSig1 || it == fSig2 ) {

	  fnTmMu++;
	}
      }
      
      if (fNorm && genPartType(gmoIdx) == 521) {
	
	frSigMu++;

	if (recMu > -0.5) {
	  
	  fnSigMu++;
	}

	if ( it == fSig1 || it == fSig2 ) {

	  fnTmMu++;
	}
      }
    }

    // -- Kaon statistics for norm. channel
    if ( TMath::Abs(mcID) == 321 ) {

      if ( genIdx > 1 ) {
	
	momIdx = fpEvt->getGenCand(genIdx)->fMom1;
      }
      
      if (fNorm && genPartType(momIdx) == 521) {
	
	frSigK++;
      }
      
      if ( it == fSig3 ) {
	
	fnTmK++;
      }
    }

    // -- Pions
    if ( TMath::Abs(mcID) == 211 ) { 

      if ( pT > 3. ) { 
	
	fMisID->Fill(0.1);
	
	((TH1D*)gDirectory->Get("m100"))->Fill(pT);   // all Pions
	((TH1D*)gDirectory->Get("m110"))->Fill(eta);
	((TH2D*)gDirectory->Get("M100"))->Fill(pT, eta);
	
	if ( recMu > -0.5 ) {
	  
	  fMisID->Fill(1.1);
	  
	  ((TH1D*)gDirectory->Get("m101"))->Fill(pT);   // mis-id. Pions
	  ((TH1D*)gDirectory->Get("m111"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M101"))->Fill(pT, eta);
	  
	}
	
	if ( recMu < -0.5 ) {
	  
	  fMisID->Fill(2.1);
	  
	  ((TH1D*)gDirectory->Get("m102"))->Fill(pT);   // correct-id. Pions
	  ((TH1D*)gDirectory->Get("m112"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M102"))->Fill(pT, eta);
	  
	}
      }
    }    

    // -- Kaons
    if ( TMath::Abs(mcID) == 321 ) { 

      if ( pT > 3. ) { 
	
	fMisID->Fill(10.1);

      ((TH1D*)gDirectory->Get("m200"))->Fill(pT);      // all Kaons
      ((TH1D*)gDirectory->Get("m210"))->Fill(eta);
      ((TH2D*)gDirectory->Get("M200"))->Fill(pT, eta);

	if ( recMu > -0.5 ) {

	  fMisID->Fill(11.1);
 
	  ((TH1D*)gDirectory->Get("m201"))->Fill(pT);   // mis-id. Kaons
	  ((TH1D*)gDirectory->Get("m211"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M201"))->Fill(pT, eta);

	}

	if ( recMu < -0.5 ) {

	  fMisID->Fill(12.1);

	  ((TH1D*)gDirectory->Get("m202"))->Fill(pT);   // correct-id. Kaons
	  ((TH1D*)gDirectory->Get("m212"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M202"))->Fill(pT, eta);

	}
      }
    }    

    // -- Protons
    if ( TMath::Abs(mcID) == 2212 ) { 

      if ( pT > 3. ) { 

	fMisID->Fill(20.1);

	((TH1D*)gDirectory->Get("m300"))->Fill(pT);      // all Protons
	((TH1D*)gDirectory->Get("m310"))->Fill(eta);
	((TH2D*)gDirectory->Get("M300"))->Fill(pT, eta);
	
	if ( recMu > -0.5 ) {

	  fMisID->Fill(21.1);
 
	  ((TH1D*)gDirectory->Get("m301"))->Fill(pT);   // mis-id. Protons
	  ((TH1D*)gDirectory->Get("m311"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M301"))->Fill(pT, eta);

	}

	if ( recMu < -0.5 ) {

	  fMisID->Fill(22.1);

	  ((TH1D*)gDirectory->Get("m302"))->Fill(pT);   // correct-id. Protons
	  ((TH1D*)gDirectory->Get("m312"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M302"))->Fill(pT, eta);

	}
      }
    }    

    // -- Global Muons
    if ( TMath::Abs(mcID) == 13 ) { 

      if ( pT > 3. ) { 

	fMisID->Fill(30.1);
	
	((TH1D*)gDirectory->Get("m400"))->Fill(pT);      // all Muons
	((TH1D*)gDirectory->Get("m410"))->Fill(eta);
	((TH2D*)gDirectory->Get("M400"))->Fill(pT, eta);
	
	if ( recMu > -0.5 ) {
	  
	  fMisID->Fill(31.1);
	  
	  ((TH1D*)gDirectory->Get("m401"))->Fill(pT);   // correct-id. Muons
	  ((TH1D*)gDirectory->Get("m411"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M401"))->Fill(pT, eta);
	  
	}

	if ( recMu < -0.5 ) {

	  fMisID->Fill(32.1);

	  ((TH1D*)gDirectory->Get("m402"))->Fill(pT);   // mis-id. Muons
	  ((TH1D*)gDirectory->Get("m412"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M402"))->Fill(pT, eta);

	}
      }
    }

    // -- Tracker Muons
    if ( TMath::Abs(mcID) == 13 ) { 

      if ( pT > 3. ) { 
	
	fMisID->Fill(40.1);
	
	((TH1D*)gDirectory->Get("m500"))->Fill(pT);      // all Muons
	((TH1D*)gDirectory->Get("m510"))->Fill(eta);
	((TH2D*)gDirectory->Get("M500"))->Fill(pT, eta);
	
	if ( trkMu > -0.5 ) {
	  
	  fMisID->Fill(41.1);
	  
	  ((TH1D*)gDirectory->Get("m501"))->Fill(pT);   // correct-id. Muons
	  ((TH1D*)gDirectory->Get("m511"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M501"))->Fill(pT, eta);
	  
	}

	if ( trkMu < -0.5 ) {

	  fMisID->Fill(42.1);

	  ((TH1D*)gDirectory->Get("m502"))->Fill(pT);   // mis-id. Muons
	  ((TH1D*)gDirectory->Get("m512"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M502"))->Fill(pT, eta);

	}
      }
    }

    // -- L1 Muons
    if ( TMath::Abs(mcID) == 13 ) { 

      if ( pT > 3. ) { 

	fMisID->Fill(50.1);
	
	((TH1D*)gDirectory->Get("m600"))->Fill(pT);      // all Muons
	((TH1D*)gDirectory->Get("m610"))->Fill(eta);
	((TH2D*)gDirectory->Get("M600"))->Fill(pT, eta);
	
	if ( l1Mu > -0.5 ) {
	  
	  fMisID->Fill(51.1);
	  
	  ((TH1D*)gDirectory->Get("m601"))->Fill(pT);   // correct-id. Muons
	  ((TH1D*)gDirectory->Get("m611"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M601"))->Fill(pT, eta);
	  
	}

	if ( l1Mu < -0.5 ) {

	  fMisID->Fill(52.1);

	  ((TH1D*)gDirectory->Get("m602"))->Fill(pT);   // mis-id. Muons
	  ((TH1D*)gDirectory->Get("m612"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M602"))->Fill(pT, eta);

	}
      }
    }
  }

  if (fDebug & 2) { cout << "muonEfficiency> End" << endl; }
}


// ---------------------------------------------------------------------
void treeBmm::fillEventStats() {

  if (fDebug & 2) { cout << "fillEventStats> Start" << endl; }

  // -- plots for the real signal
  signalPlots();

  // -- data sample control plots
  fNgen->Fill(fgMu);
  fNrec->Fill(frMu);
  fNglb->Fill(fnMu);


  fNgenJ->Fill(fgJ);
  fNdecJ->Fill(fgJmm);

  fNgenB->Fill(fgB);
  fNdecB->Fill(fgBmm);

  fErec->Fill(frMu/(1.*fgMu));
  fEglb->Fill(fnMu/(1.*frMu));

  if ( fEvent < 1000000 ) {

    fNR0->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fGoodEvent);
    fNR1->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fnB);
    fNR2->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fgB);
    fNR3->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), frMu);
    fNR4->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fnMu);
    fNR5->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fGoodL1);
    fNR6->SetBinContent(fEvent - 1000*int(fEvent/1000), int(fEvent/1000), fGoodHLT);
  }

  if (fDebug & 2) { cout << "fillEventStats> End" << endl; }
}


// ---------------------------------------------------------------------
void treeBmm::signalPlots() {

  if (fDebug & 2) { cout << "signalPlots> Start" << endl; }
  
  TAnaTrack  *L1, *L2;
  TAnaCand   *B;

  if ( fSig1 >= 0 && fSig2 >= 0 && fSigB >= 0) {

    L1 = fpEvt->getRecTrack(fSig1);
    L2 = fpEvt->getRecTrack(fSig2);

    fgMuL0 = L1->fMuID;
    fgMuL1 = L2->fMuID;

    fgQL0 = L1->fQ;
    fgQL1 = L2->fQ;

    fgPtL0  = L1->fPlab.Pt();
    fgPtL1  = L2->fPlab.Pt();
    
    fgEtaL0 = L1->fPlab.Eta();
    fgEtaL1 = L2->fPlab.Eta();

    double dphi0 = L1->fPlab.DeltaPhi(L2->fPlab);
    double deta0 = L1->fPlab.Eta() - L2->fPlab.Eta();

    fgRMM = TMath::Sqrt(dphi0*dphi0 + deta0*deta0);

    B  = fpEvt->getCand(fSigB);

    fgChi2 = B->fVtx.fChi2;
    fgS3d  = B->fVtx.fD3dE;
  
  } else {

    return;

  }

  // -- Leading Muon / best RMM 
  double pT_max(0.), pT(0.), delta_pT(-20);
  int leading(-1);
  
  double rmm_min(99.), rmm(99.), delta_rmm(-20);
  double dphi(0.), deta(0.);
  int rank1(0), rank2(0), rank_chi2(0);
  int best1(-1), best2(-1);

  int id(-1), bmmsel(-1), cand_s3d(0);
  double chi2(-1.), chi2_min(99.), delta_chi2(-20);
  double s3d(99.);

  TAnaTrack  *tt, *ts;
   
  if (fGoodEvent) {

    // -- Selection by cand. properties
    for (int i = 0; i < fpEvt->nCands(); ++i) {
    
      id     = int(fpEvt->getCand(i)->fType/10);
      bmmsel = fpEvt->getCand(i)->fType - 10*id;
      
      if ( ((!fNorm && id == 531) || (fNorm && id == 521))
	   && (bmmsel == fSel) ) {

	cand_s3d++;	

	B  = fpEvt->getCand(i);
	
	chi2 = B->fVtx.fChi2;
	s3d  = B->fVtx.fD3dE;

	if (chi2 < chi2_min) {
	  chi2_min = chi2;
	}

	if (chi2 < fgChi2) {
	  rank_chi2++;
	}

	if (s3d > 1.0) {

	  cand_s3d--;
	}
      }
    }
	
    // -- Selection by muon properties
    for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
      
      tt = fpEvt->getRecTrack(it);
      
      if ( (fSel == 3 && tt->fMuID > -0.5) ||
	   (fSel == 2 && TMath::Abs(tt->fMCID) == 13) ) {
	
	// -- Leading muon of event
	pT = tt->fPlab.Pt();
	
	if ( pT > pT_max ) { 
	  
	  pT_max = pT;
	  leading = it;
	}

	if ( pT > fgPtL0 ) { rank1++; }
	if ( pT > fgPtL1 ) { rank2++; }

        
	// -- Best RMM of event
	for (int is = it + 1; is < fpEvt->nRecTracks(); ++is) {
	  
	  ts  = fpEvt->getRecTrack(is);
	  
	  if ( (fSel == 3 && ts->fMuID > -0.5) ||
	       (fSel == 2 && TMath::Abs(ts->fMCID) == 13) ) {
	  
	    dphi = ts->fPlab.DeltaPhi(tt->fPlab);
	    deta = ts->fPlab.Eta() - tt->fPlab.Eta();
	    rmm  = TMath::Sqrt(dphi*dphi + deta*deta);
	    
	    if ( rmm <  rmm_min) { 

	      rmm_min  = rmm;
	      best1 = it;
	      best2 = is;
	    }
	  }
	}
      }
    }


    int l1(0), l2(0);
    if ( fSig1 == fSTI[0] || fSig1 == fSTI[1] ) {
      l1 = 1;
    }
    
    if ( fSig2 == fSTI[0] || fSig2 == fSTI[1] ) {
      l2 = 1;
    }
    
    int isLeading1(0), isLeading2(0);
    if ( fSig1 == leading ) {
      isLeading1 = 1;
    }
    
    if ( fSig2 == leading ) {
      isLeading2 = 1;
    }
    
    int isBest1(0), isBest2(0);
    if ( fSig1 == best1 || fSig1 == best2 ) {
      isBest1 = 1;
    }

    if ( fSig2 == best1 || fSig2 == best2 ) {
      isBest2 = 1;
    }


    if ( fgPtL0 > fgPtL1 ) {

      delta_pT = pT_max - fgPtL0;

    } else {

      delta_pT = pT_max - fgPtL1;
    }

    delta_rmm =  fgRMM - rmm_min;   
    delta_chi2 = fgChi2 - chi2_min;    
    
    // -- Signal Muon porperties
    ((TH2D*)gDirectory->Get("S106"))->Fill(fgS3d, cand_s3d);
    ((TH2D*)gDirectory->Get("S107"))->Fill(l1 - 0.5, l2 - 0.5);
    
    ((TH1D*)gDirectory->Get("s100"))->Fill(delta_chi2);
    ((TH1D*)gDirectory->Get("s101"))->Fill(delta_pT);
    ((TH1D*)gDirectory->Get("s102"))->Fill(delta_rmm);
    ((TH1D*)gDirectory->Get("s103"))->Fill(rank_chi2);
    ((TH1D*)gDirectory->Get("s104"))->Fill(fgB);
    
    ((TH2D*)gDirectory->Get("S108"))->Fill(rank1  + 0.1, rank2  + 0.1);
    ((TH2D*)gDirectory->Get("S109"))->Fill(isLeading1  - 0.5, isLeading2  - 0.5);
    ((TH2D*)gDirectory->Get("S110"))->Fill(isBest1  - 0.5, isBest2  - 0.5);
    ((TH2D*)gDirectory->Get("S111"))->Fill(isBest1 + isBest2  - 0.5, isLeading1 + isLeading2  - 0.5);
   
  }

    
  // ================================================================================
  // -- Generator level histograms
  // ================================================================================

  // -- Muons gen. index
  int m0(-99), m1(-99), m2(-99);

  if ( fSig1 >= 0 && fSig2 >= 0 ) {
    m0 = fpEvt->getRecTrack(fSig1)->fGenIndex;
    m1 = fpEvt->getRecTrack(fSig2)->fGenIndex;
  }

  if ( fSig3 >= 0 ) { 
    m2 = fpEvt->getRecTrack(fSig3)->fGenIndex;
  }

  if ( m0 > 1 && m1 > 1 ) {

    TGenCand *gm0 = fpEvt->getGenCand(m0);
    TGenCand *gm1 = fpEvt->getGenCand(m1);
    
    double dphi = gm0->fP.DeltaPhi(gm1->fP);
    double deta = gm0->fP.Eta() - gm1->fP.Eta();
    double gRMM = TMath::Sqrt(dphi*dphi + deta*deta);

    ((TH1D*)gDirectory->Get("g120"))->Fill(gRMM);
    
    double gpT  = gm0->fP.Pt();
    ((TH1D*)gDirectory->Get("g100"))->Fill(gpT);
    ((TH1D*)gDirectory->Get("g105"))->Fill(gpT);
    
    gpT  = gm1->fP.Pt();
    ((TH1D*)gDirectory->Get("g100"))->Fill(gpT);
    ((TH1D*)gDirectory->Get("g106"))->Fill(gpT);
        
    double geta  = gm0->fP.Eta();
    ((TH1D*)gDirectory->Get("g102"))->Fill(geta);
    ((TH1D*)gDirectory->Get("g107"))->Fill(geta);

    geta  = gm1->fP.Eta();
    ((TH1D*)gDirectory->Get("g102"))->Fill(geta);
    ((TH1D*)gDirectory->Get("g108"))->Fill(geta);
    
    TLorentzVector gmm  = gm0->fP + gm1->fP;
    double gmass = gmm.M();
    ((TH1D*)gDirectory->Get("g130"))->Fill(gmass);  

    // -- Kaon
    if ( m2 > 1 ) {

      TGenCand *gm0 = fpEvt->getGenCand(m0);
      TGenCand *gm1 = fpEvt->getGenCand(m1);
      TGenCand *gm2 = fpEvt->getGenCand(m2);
      
      TLorentzVector gjpsi  = gm0->fP + gm1->fP;
      
      double dphi = gjpsi.DeltaPhi(gm2->fP);
      double deta = gjpsi.Eta() - gm2->fP.Eta();
      double gRKJ = TMath::Sqrt(dphi*dphi + deta*deta);
      
      dphi = gjpsi.DeltaPhi(gm0->fP);
      deta = gjpsi.Eta() - gm0->fP.Eta();
      double gRMJ1 = TMath::Sqrt(dphi*dphi + deta*deta);
      
      dphi = gjpsi.DeltaPhi(gm1->fP);
      deta = gjpsi.Eta() - gm1->fP.Eta();
      double gRMJ2 = TMath::Sqrt(dphi*dphi + deta*deta);

      ((TH1D*)gDirectory->Get("g121"))->Fill(gRKJ);
      ((TH1D*)gDirectory->Get("g122"))->Fill(gRMJ1);
      ((TH1D*)gDirectory->Get("g123"))->Fill(gRMJ2);
      
      double gpT  = gm2->fP.Pt();
      ((TH1D*)gDirectory->Get("g110"))->Fill(gpT);
      
      double geta  = gm2->fP.Eta();
      ((TH1D*)gDirectory->Get("g112"))->Fill(geta);
      
      TLorentzVector gkj  = gm0->fP + gm1->fP + gm2->fP;
      double gmass = gkj.M();
      ((TH1D*)gDirectory->Get("g140"))->Fill(gmass);
    }
  }

  if (fDebug & 2) { cout << "signalPlots> End" << endl; }
}

// ----------------------------------------------------------------------
void treeBmm::ptResiduals() { 

  if (fDebug & 2) { cout << "ptResiduals> Start" << endl; }

  // -- pT resolution for misalignment studies

  TGenCand *pG, *pM, *pGM;
  TVector3 ptv;

  int gIndex(-1), mom1(-1), gmo1(-1);
  int mtmp1(-1), gtmp1(-1), mtmp2(-1), gtmp2(-1);
  double muonID(-1);
  double genPt(0.), genEta(0.),trkPt(0.), res0(0.), res1(0.);

  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {   

    ptv   = fpEvt->getRecTrack(it)->fPlab;
    trkPt = ptv.Pt();

    gIndex = fpEvt->getRecTrack(it)->fGenIndex;
    muonID = fpEvt->getRecTrack(it)->fMuID;

    if ( gIndex > 1 ) {
       
      pG      = fpEvt->getGenCand(gIndex);
      genPt   = pG->fP.Pt();
      genEta  = TMath::Abs(pG->fP.Eta());
      
      res0    = trkPt - genPt;
      if ( genPt > 0 ) res1    = (trkPt - genPt)/genPt;

       // -- pT Resolution of GlobalMuons
      if ( muonID > -0.5 ) {

	((TProfile*)gDirectory->Get("r100"))->Fill(genPt, res0, 1);
	((TProfile*)gDirectory->Get("r101"))->Fill(genPt, res1, 1);
	((TH1D*)gDirectory->Get("r100C"))->Fill(genPt);

	((TH2D*)gDirectory->Get("R100"))->Fill(genPt, res0);
	((TH2D*)gDirectory->Get("R101"))->Fill(genPt, res1);
	
	if ( genPt > 8. && genPt < 12. ) {
 
	  ((TProfile*)gDirectory->Get("r102"))->Fill(genEta, res0, 1);
	  ((TProfile*)gDirectory->Get("r103"))->Fill(genEta, res1, 1);

	  ((TH1D*)gDirectory->Get("r102C"))->Fill(genEta);
 
	  ((TH2D*)gDirectory->Get("R102"))->Fill(genEta, res0);
	  ((TH2D*)gDirectory->Get("R103"))->Fill(genEta, res1);
	}
      }

      // -- pT Resolution of all tracks
      ((TProfile*)gDirectory->Get("r200"))->Fill(genPt, res0, 1);
      ((TProfile*)gDirectory->Get("r201"))->Fill(genPt, res1, 1);
      
      ((TH1D*)gDirectory->Get("r200C"))->Fill(genPt);

      ((TH2D*)gDirectory->Get("R200"))->Fill(genPt, res0);
      ((TH2D*)gDirectory->Get("R201"))->Fill(genPt, res1);
      
      if ( genPt > 8. && genPt < 12. ) {
	
	((TProfile*)gDirectory->Get("r202"))->Fill(genEta, res0, 1); 
	((TProfile*)gDirectory->Get("r203"))->Fill(genEta, res1, 1); 
	
	((TH1D*)gDirectory->Get("r202C"))->Fill(genEta);
	
	((TH2D*)gDirectory->Get("R202"))->Fill(genEta, res0); 
	((TH2D*)gDirectory->Get("R203"))->Fill(genEta, res1); 
      }

      // -- origin of signal track (mother / grand-mother)
      if (pG->fMom1 >= 0) { pM  = fpEvt->getGenCand(pG->fMom1); mom1 = pM->fID;  }
      if (pM->fMom1 >= 0) { pGM = fpEvt->getGenCand(pM->fMom1); gmo1 = pGM->fID; }
   
      if (it == fSTI[0]) {
	
	((TH1D*)gDirectory->Get("d100"))->Fill(mom1);
	((TH1D*)gDirectory->Get("d103"))->Fill(gmo1);

	mtmp1 = mom1; gtmp1 = gmo1;
      }
      
      if (it == fSTI[1]) {
	
	((TH1D*)gDirectory->Get("d101"))->Fill(mom1);
	((TH1D*)gDirectory->Get("d104"))->Fill(gmo1);

	mtmp2 = mom1; gtmp2 = gmo1;
      }
      
      if (it == fSTI[2] && fD2 > 0) {
	
	((TH1D*)gDirectory->Get("d102"))->Fill(mom1);
	((TH1D*)gDirectory->Get("d105"))->Fill(gmo1);
      }
    }   
  }
  
  ((TH2D*)gDirectory->Get("D100"))->Fill(mtmp1, mtmp2);
  ((TH2D*)gDirectory->Get("D101"))->Fill(gtmp1, gtmp2);

  if (fDebug & 2) { cout << "ptResiduals> End" << endl; }
}

// ----------------------------------------------------------------------
void treeBmm::vertexResiduals() {


  if (fDebug & 2) { cout << "vertexResiduass> Start" << endl; }

  fDpVtxX = fDpVtxY = fDpVtxZ = -99999.;
  fDsVtxX = fDsVtxY = fDsVtxZ = -99999.;
  fDvtxPerp = fDvtxPar = -99999.;
  fFltR   = fFltG   = fFltRes   = -99999.;
  fFltR3D = fFltG3D = fFltRes3D = -99999.;

  TH1D *h;

  // -- Resolution of primary and secondary Vertex (not working in CMSSW_1_3_X)

  TVector3 resPV = fpEvt->fPrimaryVertex2.fSimPoint - fpEvt->fPrimaryVertex2.fPoint;

//   fDpVtxX      = resPV.X(); 
//   fDpVtxY      = resPV.Y(); 
//   fDpVtxZ      = resPV.Z(); 
 

  TVector3 resSV = fpB->fVtx.fSimPoint - fpB->fVtx.fPoint;

//   fDsVtxX      = resSV.X(); 
//   fDsVtxY      = resSV.Y(); 
//   fDsVtxZ      = resSV.Z(); 


  // -- Find generated prim. and secondary Vertex

  TGenCand *pJ; // --> J/Psi-meson
  TGenCand *pM; // --> B-meson
  TGenCand *pG; // --> Parent of B-meson

  TVector3 simSecVtx;
  TVector3 simPrimVtx;

  int simVtx(0), truth(0);
  double genPt(0.), genPlab(0.), mass(0.);

  if ( fpL1->fGenIndex > 1 && fpL2->fGenIndex > 1 ) {

    TGenCand *g1 = fpEvt->getGenCand(fpL1->fGenIndex);
    TGenCand *g2 = fpEvt->getGenCand(fpL2->fGenIndex);

    if ( !fNorm ) {

      // -- Signal
      if ( (g1->fMom1 >= 0) && genPartType(g1->fMom1) == 531 && (g2->fMom1 == g1->fMom1) ) {
	
	pM    = fpEvt->getGenCand(g1->fMom1);
	
	if ( pM->fMom1 >= 0 ) { 
	  
	  pG    = fpEvt->getGenCand(pM->fMom1);
	  truth = 1;
	}
      }

    } else {
      
      // -- Normalization
      if ( (g1->fMom1 >= 0) && genPartType(g1->fMom1) == 443 && (g2->fMom1 == g1->fMom1) ) {
      
	pJ    = fpEvt->getGenCand(g1->fMom1);
	
	if ( fpK->fGenIndex > 1 ) {
	  
	  TGenCand *g3 = fpEvt->getGenCand(fpK->fGenIndex);
	  
	  if ( (g3->fMom1 >= 0) && genPartType(g3->fMom1) == 521 && (g3->fMom1 == pJ->fMom1 ) ) {
	    
	    pM    = fpEvt->getGenCand(g3->fMom1);
	    
	    if ( pM->fMom1 >= 0 ) { 
	      
	      pG    = fpEvt->getGenCand(pM->fMom1);
	      truth = 1;
	    }
	  }
	}
      }
    }
	  
    if ( truth && fpL1->fPlab.Pt() > 3. && fpL2->fPlab.Pt() > 3. && fLxy/fSxy > 2. ) { 
	    
      simPrimVtx.SetXYZ( 0.1*pG->fV.x(),
			 0.1*pG->fV.y(),
			 0.1*pG->fV.z() );
	    
      simSecVtx.SetXYZ( 0.1*pM->fV.x(),
			0.1*pM->fV.y(),
			0.1*pM->fV.z() );
	    
      genPt    = pM->fP.Perp();
      genPlab  = pM->fP.P();
      mass     = pM->fMass;

      if ( mass < 0.5 ) mass = fMassB;
	    
      simVtx = 1;
	    
    }
  }

  // -------------------------------------------------------------------
 
  if (fDebug & 4) {

    cout << "PrimVtx (gen) " << " " << simPrimVtx.x() << " " << simPrimVtx.y() << " " << simPrimVtx.z() << endl;
    cout << "PrimVtx (rec) " << " " << fpEvt->fPrimaryVertex2.fPoint.x() 
	                     << " " << fpEvt->fPrimaryVertex2.fPoint.y() 
	                     << " " << fpEvt->fPrimaryVertex2.fPoint.z() << endl;

    cout << "SecVtx (gen) " << " " << simSecVtx.x() << " " << simSecVtx.y() << " " << simSecVtx.z() << endl;
    cout << "SecVtx (rec) " << " " << fpB->fVtx.fPoint.x() 
                            << " " << fpB->fVtx.fPoint.y() 
                            << " " << fpB->fVtx.fPoint.z() << endl;

  }


  // -------------------------------------------------------------------

  //**   TVector3 Res = fpB->fVtx.fSimPoint -  fpB->fVtx.fPoint; 
  TVector3 Res = simSecVtx -  fpB->fVtx.fPoint; 
  double theta = fpB->fPlab.Angle(Res);

  fDvtxPerp    = Res.Mag()*TMath::Sin(theta);  
  fDvtxPar     = Res.Mag()*TMath::Cos(theta);

  // -------------------------------------------------------------------

  if ( !simVtx ) {
    
    return;
  }


  // -- Vertex resolution of primary and secondary Vertex

  TVector3 primVtxRes = simPrimVtx - fpEvt->fPrimaryVertex2.fPoint;
  TVector3 secVtxRes  = simSecVtx  - fpB->fVtx.fPoint;

  fDpVtxX      = primVtxRes.X(); 
  fDpVtxY      = primVtxRes.Y(); 
  fDpVtxZ      = primVtxRes.Z(); 

  fDsVtxX      = secVtxRes.X(); 
  fDsVtxY      = secVtxRes.Y(); 
  fDsVtxZ      = secVtxRes.Z(); 


  // -- Proper time resolution

  //**   TVector3 fltGen = fpB->fVtx.fSimPoint - fpEvt->fPrimaryVertex2.fSimPoint;
  //**   TVector3 fltGen = simSecVtx            - fpEvt->fPrimaryVertex2.fSimPoint;
  TVector3 fltGen = simSecVtx - simPrimVtx;
  TVector3 fltRec = fpB->fVtx.fPoint    - fpEvt->fPrimaryVertex2.fPoint;

  if ( fltGen.Mag() == 0 || fltRec.Mag() == 0 ) { 

    return; 
  }

  if ( fPt > 0 && genPt > 0 ) {
              
    fFltR   = fltRec.Perp()*fMass/fPt;
    fFltG   = fltGen.Perp()*mass/genPt;
    fFltRes = fFltR - fFltG;

  } else {

    return;
  }

  if (fDebug & 4) {

    cout <<" FltRec: " << fFltR << endl;
    cout <<" FltGen: " << fFltG << endl;
    cout <<" FltRes: " << fFltRes << endl;
  }

  if ( fP > 0 && genPlab > 0 ) {

    fFltR3D   = fltRec.Mag()*fMass/fP;
    fFltG3D   = fltGen.Mag()*mass/genPlab;
    fFltRes3D = fFltR - fFltG;

  } else {

    return;
  }


  if (fDebug & 4) {

    cout <<" FltRec3D: " << fFltR3D << endl;
    cout <<" FltGen3D: " << fFltG3D << endl;
    cout <<" FltRes3D: " << fFltRes3D << endl;
  }

  double DPVxy = primVtxRes.Perp();
  double DPV3d = primVtxRes.Mag();
  double DSVxy = secVtxRes.Perp();
  double DSV3d = secVtxRes.Mag();

  double DLxy = fltRec.Perp() - fltGen.Perp();
  double DL3d = fltRec.Mag()  - fltGen.Mag();

  double DTxy = fltRec.Perp()*fMass/fPt - fltGen.Perp()*mass/genPt;
  double DT3d = fltRec.Mag()*fMass/fP - fltGen.Mag()*mass/genPlab;

  double DTxy_s = (0.01*1.E15/299792458.)*(fltRec.Perp()*fMass/fPt - fltGen.Perp()*mass/genPt);
  double DT3d_s = (0.01*1.E15/299792458.)*(fltRec.Mag()*fMass/fP - fltGen.Mag()*mass/genPlab);

  if (fDebug & 2) {

    cout << " DLxy = " << DLxy
	 << " DL3d = " << DL3d << endl;

    cout << " DTxy = " << DTxy
	 << " DT3d = " << DT3d << endl;

    cout << " DTxy_s = " << DTxy_s
	 << " DT3d_s = " << DT3d_s << endl;
  }
  
  h = (TH1D*)gDirectory->Get("v200"); h->Fill(DLxy);
  h = (TH1D*)gDirectory->Get("v201"); h->Fill(DL3d);
  h = (TH1D*)gDirectory->Get("v202"); h->Fill(DTxy);
  h = (TH1D*)gDirectory->Get("v203"); h->Fill(DT3d);
  h = (TH1D*)gDirectory->Get("v204"); h->Fill(DTxy_s);
  h = (TH1D*)gDirectory->Get("v205"); h->Fill(DT3d_s);
  
  h = (TH1D*)gDirectory->Get("v300"); h->Fill(DPVxy);
  h = (TH1D*)gDirectory->Get("v301"); h->Fill(DPV3d);
  h = (TH1D*)gDirectory->Get("v400"); h->Fill(DSVxy);
  h = (TH1D*)gDirectory->Get("v401"); h->Fill(DSV3d);


  if (fDebug & 2) { cout << "vertexResiduals> End" << endl; }
}


// **************************************************************************************************
// --------------------------------------------------------------------------------------------------
void treeBmm::readCuts(TString filename, int dump, double ptMin, double etaMax) {

  if (fDebug & 2) { cout << "readCuts> Start" << endl; }

  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];

  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());
  fn.ReplaceAll("bjk", "");
  fn.ReplaceAll("bmm", "");
  fn.ReplaceAll("cuts", "");
  fn.ReplaceAll(".", "");
  fn.ReplaceAll("/", "");

  fChannel = fn;

  if (dump) {
    cout << "Setting generator thresholds to pT = " << ptMin <<  " and eta = " << etaMax << endl;
    cout << "====================================" << endl;
    cout << "Cut file  " << fCutFile.Data() << endl;
    cout << "Cut label " << fn.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  fMinPt = ptMin;
  fMaxEta = etaMax;

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin; 

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!strcmp(CutName, "PTLO")) {
      PTLO = CutValue; ok = 1;
      if (dump) cout << "PTLO:           " << PTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, PTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(l) [GeV]");
      ((TH1D*)gDirectory->Get("PTLO"))->SetTitle(Form("p_{T}^{#mu} > %4.2f", PTLO));
    }

    if (!strcmp(CutName, "PTHI")) {
      PTHI = CutValue; ok = 1;
      if (dump) cout << "PTHI:           " << PTHI << " GeV" << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, PTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(l) [GeV]");
    }

    if (!strcmp(CutName, "ETALO")) {
      ETALO = CutValue; ok = 1;
      if (dump) cout << "ETALO:           " << ETALO << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, ETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta_{B}^{min}");
    }

    if (!strcmp(CutName, "ETAHI")) {
      ETAHI = CutValue; ok = 1;
      if (dump) cout << "ETAHI:           " << ETAHI << endl;
      ibin = 14;
      hcuts->SetBinContent(ibin, ETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta_{B}^{max}");
      ((TH1D*)gDirectory->Get("ETABS"))->SetTitle(Form("%4.2f < #eta_{B} < %4.2f", ETALO, ETAHI));
    }

    if (!strcmp(CutName, "RMMLO")) {
      RMMLO = CutValue; ok = 1;
      if (dump) cout << "RMMLO:           " << RMMLO << endl;
      ibin = 15;
      hcuts->SetBinContent(ibin, RMMLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{#mu#mu}^{min}");
    }

    if (!strcmp(CutName, "RMMHI")) {
      RMMHI = CutValue; ok = 1;
      if (dump) cout << "RMMHI:           " << RMMHI << endl;
      ibin = 16;
      hcuts->SetBinContent(ibin, RMMHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{#mu#mu}^{max}");
      ((TH1D*)gDirectory->Get("RMM"))->SetTitle(Form("%4.2f < R_{#mu#mu} < %4.2f", RMMLO, RMMHI));
    }

    if (!strcmp(CutName, "TIPHI")) {
      TIPHI = CutValue; ok = 1;
      if (dump) cout << "TIPHI:           " << TIPHI << endl;
      ibin = 17;
      hcuts->SetBinContent(ibin, TIPHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "TIP(l) [cm]");
    }

    if (!strcmp(CutName, "PTBS")) {
      PTBS = CutValue; ok = 1;
      if (dump) cout << "PTBS:           " << PTBS << " GeV" << endl;
      ibin = 100; 
      hcuts->SetBinContent(ibin, PTBS);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}(B_{s}) [GeV]");
      ((TH1D*)gDirectory->Get("PTBS"))->SetTitle(Form("p_{T}(B_{s}) > %4.2f", PTBS));
    }

    if (!strcmp(CutName, "VTXCHI")) {
      VTXCHI = CutValue; ok = 1;
      if (dump) cout << "VTXCHI:           " << VTXCHI << endl;
      ibin = 102;
      hcuts->SetBinContent(ibin, VTXCHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#chi^2");
      ((TH1D*)gDirectory->Get("VTXCHI"))->SetTitle(Form("#chi^2 < %4.2f", VTXCHI));
      ((TH1D*)gDirectory->Get("VTXCHI_F"))->SetTitle(Form("#chi^2 < %4.2f", VTXCHI));
      ((TH1D*)gDirectory->Get("VTXCHI_PRE"))->SetTitle(Form("#chi^2 < %4.2f", VTXCHI));
    }

    if (!strcmp(CutName, "L3DS3DLO")) {
      L3DS3DLO = CutValue; ok = 1;
      if (dump) cout << "L3DS3DLO:           " << L3DS3DLO << endl;
      ibin = 112;
      hcuts->SetBinContent(ibin, L3DS3DLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{3D}/#sigma_{3D}");
      ((TH1D*)gDirectory->Get("L3DS3D"))->SetTitle(Form("l_{3D}/#sigma_{3D} > %4.2f", L3DS3DLO));
      ((TH1D*)gDirectory->Get("L3DS3D_PRE"))->SetTitle("l_{3D}/#sigma_{3D} > 7");
    }
    
    if (!strcmp(CutName, "L3DLO")) {
      L3DLO = CutValue; ok = 1;
      if (dump) cout << "L3DLO:           " << L3DLO << endl;
      ibin = 103;
      hcuts->SetBinContent(ibin, L3DLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{3d} [cm]");
    }

    if (!strcmp(CutName, "LXYSXYLO")) {
      LXYSXYLO = CutValue; ok = 1;
      if (dump) cout << "LXYSXYLO:           " << LXYSXYLO << endl;
      ibin = 104;
      hcuts->SetBinContent(ibin, LXYSXYLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{xy}/#sigma_{xy}");
      ((TH1D*)gDirectory->Get("LXYSXY"))->SetTitle(Form("l_{xy}/#sigma_{xy} > %4.2f", LXYSXYLO));
      ((TH1D*)gDirectory->Get("LXYSXY_PRE"))->SetTitle("l_{xy}/#sigma_{xy} > 7");
    }

    if (!strcmp(CutName, "LXYLO")) {
      LXYLO = CutValue; ok = 1;
      if (dump) cout << "LXYLO:           " << LXYLO << endl;
      ibin = 105;
      hcuts->SetBinContent(ibin, LXYLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{xy} [cm]");
    }

    if (!strcmp(CutName, "COSALPHA")) {
      COSALPHA = CutValue; ok = 1;
      if (dump) cout << "COSALPHA:        " << COSALPHA << endl;
      ibin = 106;
      hcuts->SetBinContent(ibin, COSALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, "cos(#alpha)");
      ((TH1D*)gDirectory->Get("COSALPHA"))->SetTitle(Form("cos(#alpha) > %4.2f", COSALPHA));
    }

    if (!strcmp(CutName, "ISOVETO")) {
      ISOVETO = int(CutValue); ok = 1;
      if (dump) cout << "ISOVETO:           " << ISOVETO << endl;
      ibin = 107;
      hcuts->SetBinContent(ibin, ISOVETO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "I_{veto}");
    }

    if (!strcmp(CutName, "ISOLATION")) {
      ISOLATION = CutValue; ok = 1;
      if (dump) cout << "ISOLATION:           " << ISOLATION << endl;
      ibin = 108;
      hcuts->SetBinContent(ibin, ISOLATION);
      hcuts->GetXaxis()->SetBinLabel(ibin, "I");
      ((TH1D*)gDirectory->Get("ISOLATION"))->SetTitle(Form("Isolation < %4.2f", ISOLATION));
      ((TH1D*)gDirectory->Get("ISOLATION_F"))->SetTitle(Form("Isolation < %4.2f", ISOLATION));
      ((TH1D*)gDirectory->Get("ISOLATION_PRE"))->SetTitle(Form("Isolation < %4.2f", ISOLATION));
    }

    if (!strcmp(CutName, "ISOCONE")) {
      ISOCONE = CutValue; ok = 1;
      if (dump) cout << "ISOCONE:           " << Form("%3.2f", ISOCONE) << endl;
      ibin = 109;
      hcuts->SetBinContent(ibin, ISOCONE);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{I}");
    }

    if (!strcmp(CutName, "BMMSEL")) {
      BMMSEL = int(CutValue); ok = 1;
      if (dump) cout << "BMMSEL:           " << BMMSEL << endl;
      ibin = 110;
      hcuts->SetBinContent(ibin, BMMSEL);
      hcuts->GetXaxis()->SetBinLabel(ibin, "BMMSEL");
    }
    
    if (!strcmp(CutName, "SUBSEL")) {
      SUBSEL = int(CutValue); ok = 1;
      if (dump) cout << "SUBSEL:           " << SUBSEL << endl;
      ibin = 111;
      hcuts->SetBinContent(ibin, SUBSEL);
      hcuts->GetXaxis()->SetBinLabel(ibin, "SUBSEL");
    }
    
    if (!strcmp(CutName, "MASSLO")) {
      MASSLO = CutValue; ok = 1;
      if (dump) cout << "MASSLO:           " << MASSLO << endl;
      ibin = 113;
      hcuts->SetBinContent(ibin, MASSLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "m_{min}^{#mu #mu}");
    }
    
    if (!strcmp(CutName, "MASSHI")) {
      MASSHI = CutValue; ok = 1;
      if (dump) cout << "MASSHI:           " << MASSHI << endl;
      ibin = 114;
      hcuts->SetBinContent(ibin, MASSHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "m_{max}^{#mu #mu}");
      ((TH1D*)gDirectory->Get("MASSBAND"))->SetTitle(Form("%4.2f < m_{#mu #mu} < %4.2f", MASSLO, MASSHI));
    }
    
    if (!strcmp(CutName, "ISOPTLO")) {
      ISOPTLO = CutValue; ok = 1;
      if (dump) cout << "ISOPTLO:           " << ISOPTLO << endl;
      ibin = 115;
      hcuts->SetBinContent(ibin, ISOPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "iso. p_{T}^{min}");
    }
    
    if (!strcmp(CutName, "MASSWI")) {
      MASSWI = CutValue; ok = 1;
      if (dump) cout << "MASSWI:           " << MASSWI << endl;
      ibin = 118;
      hcuts->SetBinContent(ibin, 1000.*MASSWI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#Delta m{#mu #mu}");
      ((TH1D*)gDirectory->Get("MASSWI"))->SetTitle(Form("#Delta m{#mu #mu} < %4.2f", MASSWI));
      ((TH1D*)gDirectory->Get("MASSWI_F"))->SetTitle(Form("#Delta m{#mu #mu} < %4.2f", MASSWI));
    }

    if (!strcmp(CutName, "PRESEL")) {
      PRESEL = CutValue; ok = 1;
      if (dump) cout << "PRESEL(L/S): " << PRESEL << endl;
      ibin = 119;
      hcuts->SetBinContent(ibin, PRESEL);
      hcuts->GetXaxis()->SetBinLabel(ibin, "presel");
    }
    
    if (!ok) cout << "==> ERROR: Don't know about variable " << CutName << endl;
  }
  
  if (fDebug & 2) { cout << "readCuts> End" << endl; }
}


// ----------------------------------------------------------------------
void treeBmm::decayChannel(TString ch, int dump) {

  if (fDebug & 2) { cout << "decayChannel> Start" << endl; }

  fD0 = -1;  fD1 = -1;  fD2 = -1;
  fNorm = 0; fSel = 3; fSubSel = 1;
  fMassB = 5.369;

  if ( !strcmp("bd2pi", ch.Data()) ) {

    fD0 = 211;  fD1 = 211; fSel = 1;
    fChannel = TString("bd2pi");
    cout << "Selected decay mode: bd2pi." << endl;

  } else if ( !strcmp("bdpik", ch.Data()) ) {

    fD0 = 321;  fD1 = 211; fSel = 1;
    fChannel = TString("bdpik");
    cout << "Selected decay mode: bdpik." << endl;

  } else if ( !strcmp("bdpimunu", ch.Data()) ) {

    fD0 = 13;  fD1 = 211; fSel = 1;
    fChannel = TString("bdpimunu");
    cout << "Selected decay mode: bdpimunu." << endl;

  } else if ( !strcmp("bs2pi", ch.Data()) ) {

    fD0 = 211;  fD1 = 211; fSel = 1;
    fChannel = TString("bs2pi");
    cout << "Selected decay mode: bs2pi." << endl;
  }
  else if ( !strcmp("bskk", ch.Data()) ) {

    fD0 = 321;  fD1 = 321; fSel = 1;
    fChannel = TString("bskk");
    cout << "Selected decay mode: bskk." << endl;

  } else if ( !strcmp("bskpi", ch.Data()) ) {

    fD0 = 321;  fD1 = 211; fSel = 1;
    fChannel = TString("bskpi");
    cout << "Selected decay mode: bskpi." << endl;

  } else if ( !strcmp("bskmunu", ch.Data()) ) {

    fD0 = 13;  fD1 = 321; fSel = 1;
    fChannel = TString("bskmunu");
    cout << "Selected decay mode: bskmunu." << endl;

  } else if ( !strcmp("lbpk", ch.Data()) ) {

    fD0 = 2212;  fD1 = 321; fSel = 1;
    fChannel = TString("lbpk");
    cout << "Selected decay mode: lbpk." << endl;

  } else if ( !strcmp("lbppi", ch.Data()) ) {

    fD0 = 2212;  fD1 = 211; fSel = 1;
    fChannel = TString("lbppi");
    cout << "Selected decay mode: lbppi." << endl;

  } else if ( !strcmp("qcd", ch.Data()) ) {

    fD0 = 321;  fD1 = 211; fSel = 1;
    fChannel = TString("qcd");
    cout << "Selected decay mode: qcd." << endl;

  } else if ( !strcmp("bsmumug", ch.Data())  || 
	      !strcmp("bsmumup0", ch.Data()) ||
	      !strcmp("bdmumup0", ch.Data()) ||
	      !strcmp("bu3munu", ch.Data())  ||
	      !strcmp("bc3munu", ch.Data())  ||
	      !strcmp("bcjpmunu", ch.Data())    ) {
    
    fD0 = 13;  fD1 = 13;
    fChannel = TString("2mu_rare");
    cout << "Selected decay mode: 2mu (from rare decays)." << endl;

  } else if ( !strcmp("bjk", ch.Data()) ) {

    fNorm = 1; fMassB = 5.279;
    fD0 = 13;  fD1 = 13;  fD2 = 321;
    fChannel = TString("bdjpk");
    cout << "Selected decay mode: bdjpk." << endl;

  } else {
    
    fD0 = 13;  fD1 = 13;
    fChannel = TString("2mu");
    cout << "Selected decay mode: 2mu." << endl;
  }

  if ( BMMSEL > 0. ) {

    fSel = int(BMMSEL);

  }

  if ( SUBSEL > 0. ) {

    fSubSel = int(SUBSEL);

  }

  cout << "-----------------------------------------------------" << endl;
  if ( fNorm ) {
  cout << "  LOOKING AT THE NORMALIZATION CHANNEL " << endl;
  } else {
  cout << "        LOOKING AT THE SIGNAL CHANNEL" << endl;
  }

  cout << "-----------------------------------------------------" << endl;
  cout << "       Muon selection with SEL = " << fSel << "." << endl;
  cout << "  Candidate selection with sub-SEL = " << fSubSel << "." << endl;
  cout << "-----------------------------------------------------" << endl;
  if ( SETL1  > 0 ) cout << "             ***  L1 = 1 ***" << endl;
  if ( SETHLT > 0 ) cout << "             *** HLT = 1 ***" << endl;
  if ( SETL1  > 0 || SETHLT > 0 ) cout << "-----------------------------------------------------" << endl;

  if ( fOffCand > 0 ) {
    cout << " ... candidate selection 1 - 4 activated" << endl;
  } else {
    cout << " ... not doing candidate selection 1 - 4!" << endl;
  }
  cout << "-----------------------------------------------------" << endl;

  if (fDebug & 2) { cout << "decayChannel> End" << endl; }
}


// **************************************************************************************************
// --------------------------------------------------------------------------------------------------
void treeBmm::setTitles(TH1 *h, const char *sx, const char *sy, float size,
			float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy);
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}


// ----------------------------------------------------------------------
void treeBmm::setTitles2(TH2 *h, const char *sx, const char *sy, float size,
			float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy);
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}


// ----------------------------------------------------------------------
void treeBmm::boxes() {
  //  cout << "  boxes " << endl;
}
