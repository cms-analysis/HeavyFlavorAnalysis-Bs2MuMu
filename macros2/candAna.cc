#include "candAna.hh"

#include "../interface/HFMasses.hh"
#include "../macros/AnalysisDistribution.hh"
#include "MuScleFitCorrector.hh"


using namespace std;
//bool select_print; //dk


struct near_track_t {
  int ix;
  float doca;
  float p;
  float pt;
  float pt_rel;
  float deltaR;
};

// ----------------------------------------------------------------------
candAna::candAna(bmm2Reader *pReader, string name, string cutsFile) {
  fpReader = pReader; 
  fVerbose = fpReader->fVerbose;	 
  fYear    = fpReader->fYear; 
  fName    = name; 
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", reading cutsfile " << cutsFile << " setup for year " << fYear << endl;

  MASSMIN = 4.9;
  MASSMAX = 5.9; 
  BLIND = 0;

  fGenBTmi = fGenM1Tmi = fGenM2Tmi = fNGenPhotons = fRecM1Tmi = fRecM2Tmi = fCandTmi = -1; 

  fHistDir = gFile->mkdir(fName.c_str());

  cout << "======================================================================" << endl;	 

}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}



// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {
  
  fpEvt = evt; 
  fBadEvent = false;


  // TESTING d.k.
  //static int count = 0; //dk
  //select_print = false; //dk
  //if(fRun==208686 && fEvt==677184169) select_print = true; //dk
  //select_print = true; //dk
  //if(!select_print) return;

  //  play(); 
  //  return;

  if (fIsMC) {
    genMatch(); 
    recoMatch(); 
    candMatch(); 
    if (fBadEvent) {
      cout << "XXXXXXX BAD EVENT XXXXXX SKIPPING XXXXX" << endl;
      return;
    }
    efficiencyCalculation();
  } 

  triggerSelection();
  runRange(); 

  //cout<<" event "<<fEvt<<" run "<<fRun<<" cands "<<fpEvt->nCands()<<" "<<fVerbose<<" "<<fGoodHLT<<" "<<fIsMC<<" "<<count<<endl;
  //return;

  // Skip data events where there was no valid trigger
  // NO!  if(!fIsMC && !fGoodHLT) {return;}  

  TAnaCand *pCand(0);
  if (fVerbose == -66) { cout << "----------------------------------------------------------------------" << endl;}
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);


    if (fVerbose == -66) {
      cout << Form("%4d", iC) << " cand -> " << pCand->fType << endl;
      continue;
    }

    if (TYPE != pCand->fType) {
      if (fVerbose > 39) cout << "  skipping candidate at " << iC << " which is of type " << pCand->fType <<endl;
      continue;
    }

    if (fVerbose > 19) 
      cout<<"--------------- found candidate ------------------- " << iC << " type = " << pCand->fType << " evt = " << fEvt << endl;

    if (fVerbose > 99) {
      
      int gen1(-1), gen2(-1), gen0(-1);
 
      if (fIsMC) {
	gen1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex)->getGenIndex())->fID;
	gen2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex)->getGenIndex())->fID;
	gen0 = fpEvt->getGenTWithIndex(fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex)->getGenIndex())->fMom1)->fID;
      }

      cout << "Analyzing candidate at " << iC << " which is of type " << TYPE << " mass "<< pCand->fMass
	   << " with sig tracks: " << pCand->fSig1 << " .. " << pCand->fSig2
	   << " and rec tracks: " 
	   << fpEvt->getSigTrack(pCand->fSig1)->fIndex
	   << " .. " 
	   << fpEvt->getSigTrack(pCand->fSig2)->fIndex
	   << " gen IDs =  " << gen1 << " " << gen2 << " from " << gen0
	   << " tm = " << fGenM1Tmi << " " << fGenM2Tmi 
	   << " ctm " << fCandTmi 
	   << endl;

      if (TMath::Abs(gen1) != 13 || TMath::Abs(gen2) != 13) fpEvt->dumpGenBlock(); 

    }



    fpCand = pCand;
    fCandIdx = iC; 

    //cout<<" call analysis "<<fpCand<<" "<<iC<<endl;

    // -- call derived functions
    candAnalysis();

 
    if(fpMuon1 != NULL && fpMuon2 != NULL) { // do only when 2 muons exist
      fHLTmatch  = doTriggerMatching(false);      // match only with the specific objects 
      fHLTmatch2 = doTriggerMatching(true); // =true - check matching to all Mu L3 HLT FOR TESTING
      //fHLTmatch2 = doTriggerMatching(); // OLD
    }

    // TESTING d.k.
    //     if(fGoodHLT) {
    //      cout<<fHLTmatch<<" "<<fHLTmatch2<<endl;
    //     cout << " candidate at " << iC << " which is of type " << TYPE << " mass "<< pCand->fMass<<endl;
    //     }
    //     //if((fHLTmatch!=fHLTmatch2)) 
    //      if(fGoodHLT && (fHLTmatch!=fHLTmatch2)) 
    //        cout<<fHLTmatch<<" "<<fHLTmatch2
    //  	  << " candidate at " << iC << " which is of type " << TYPE << " mass "<< pCand->fMass<<" "<<fGoodHLT<<endl;

    if (0 && fCandM > 4.99 && fCandM < 5.02 && fCandFLS3d > 13 && fCandA < 0.05 && fMu1Pt > 4.5 && fMu2Pt > 4.5) {
      cout << "----------------------------------------------------------------------" << endl;
      int gen1(-1), gen2(-1), gen0(-1);
      gen1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex)->getGenIndex())->fID;
      gen2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex)->getGenIndex())->fID;
      gen0 = fpEvt->getGenTWithIndex(fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex)->getGenIndex())->fMom1)->fID;

      cout << "Analyzing candidate at " << iC << " which is of type " << TYPE 
	   << " with sig tracks: " << pCand->fSig1 << " .. " << pCand->fSig2
	   << " and rec tracks: " 
	   << fpEvt->getSigTrack(pCand->fSig1)->fIndex
	   << " .. " 
	   << fpEvt->getSigTrack(pCand->fSig2)->fIndex 
	   << " gen IDs =  " << gen1 << " " << gen2 << " from " << gen0
	   << " tm = " << fGenM1Tmi << " " << fGenM2Tmi 
	   << " ctm " << fCandTmi 
	   << endl;
      cout << fpCand->fMass << endl;

      TLorentzVector gm1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex)->getGenIndex())->fP;
      TLorentzVector gm2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex)->getGenIndex())->fP;
      cout << "gen level mass: " << (gm1+gm2).M() << endl;
      cout << "----------------------------------------------------------------------" << endl;
      fpEvt->dumpGenBlock();
    }

    if (fIsMC) {
      fTree->Fill(); 
      if (!fGoodMuonsID) fAmsTree->Fill();

      //if (fVerbose > 10) cout<<" write "<<fpCand->fType<<" "<<fEvt<<" "<<fGoodHLT<<" "<<fHLTmatch<<endl;

      ((TH1D*)fHistDir->Get("test3"))->Fill(6.); 
      if(fJSON)       ((TH1D*)fHistDir->Get("test3"))->Fill(7.); 
      if(fJSON&&fGoodHLT)       ((TH1D*)fHistDir->Get("test3"))->Fill(8.); 
      if(fJSON&&fGoodHLT&&fHLTmatch)       ((TH1D*)fHistDir->Get("test3"))->Fill(9.);  

    } else {  // DATA
      if (NOPRESELECTION) fPreselection = true; 
      if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX) {
	if (fPreselection && !fGoodMuonsID) fAmsTree->Fill();
	// do nothing
	//cout<<" blinded "<<BLIND<<" "<<fpCand->fMass<<" "<<fCandM<<" "<<SIGBOXMIN<<" "<<SIGBOXMAX<<" "<<fCandIso<<" "<<fPreselection<<endl;;

      } else {

	//if (fVerbose > 1)
	//cout<<" select "<<fRun<<" "<<fLS<<" "<<fEvt<<" "<<fJSON<<" "<<fGoodHLT<<" "<<fpCand->fType<<" "<<fPreselection<<" "<<fHLTmatch<<endl;

	if (fPreselection) { 
	  ((TH1D*)fHistDir->Get("test3"))->Fill(6.); 
	  if (fVerbose > 5) cout << " filling this cand into the tree" << endl;	  
	  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(12); 

 	  //if (fVerbose > 10)
	  //cout<<" write "<<fpCand->fType<<" "<<fRun<<" "<<fLS<<" "<<fEvt<<" "<<fJSON<<" "<<fGoodHLT<<" "
	  //<<fHLTmatch<<endl;

	  ((TH1D*)fHistDir->Get("run1"))->Fill(fRun); 
	  if(fJSON) ((TH1D*)fHistDir->Get("run2"))->Fill(fRun);
	  if(fJSON&&fGoodHLT) ((TH1D*)fHistDir->Get("run3"))->Fill(fRun);
	  if(fJSON&&fGoodHLT&&fHLTmatch) ((TH1D*)fHistDir->Get("run4"))->Fill(fRun); 

	  fTree->Fill(); 
	  if (!fGoodMuonsID) fAmsTree->Fill();
	  
	} else {	 

           if ( fVerbose > 5 ) cout << " failed preselection" << endl;        

	} // if preselection
      } // if blind
    } // if MC  
  }  // loop over cands

}

// ----------------------------------------------------------------------
void candAna::candAnalysis() {

  //cout<<" cand "<<fpCand<<endl;

  if (0 == fpCand) return;

  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(1); 

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

  //cout << "fPvN = " << fPvN << " "<<fpCand->fPvIdx <<" "<<fpCand->fPvIdx2 <<" "<<fpEvt->nPV()<<endl;

  if (fpCand->fPvIdx > -1 && fpCand->fPvIdx < fpEvt->nPV()) {
    TAnaVertex *pv = fpEvt->getPV(fpCand->fPvIdx); 
    fPvX = pv->fPoint.X(); 
    fPvY = pv->fPoint.Y(); 
    fPvZ = pv->fPoint.Z(); 
    fPvNtrk = pv->getNtracks();
    fPvNdof = pv->fNdof;
    fPvAveW8 = ((fPvNdof+2.)/2.)/fPvNtrk;
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

  fCandType     = fpCand->fType;
  fCandPt       = fpCand->fPlab.Perp();
  fCandP        = fpCand->fPlab.Mag();
  fCandEta      = fpCand->fPlab.Eta();
  fCandPhi      = fpCand->fPlab.Phi();
  fCandM        = fpCand->fMass;
  fCandME       = fpCand->fMassE;
  fCandDoca     = fpCand->fMaxDoca;

  // -- values of cand wrt PV
  fCandPvTip    = fpCand->fPvTip;
  fCandPvTipE   = fpCand->fPvTipE;
  fCandPvTipS   = fCandPvTip/fCandPvTipE;
  fCandPvLip    = fpCand->fPvLip;
  fCandPvLipE   = fpCand->fPvLipE;
  fCandPvLipS   = fCandPvLip/fCandPvLipE;
  fCandPvLip2   = fpCand->fPvLip2;
  fCandPvLipS2  = fpCand->fPvLip2/fpCand->fPvLipE2;
  if (fpCand->fPvLip2 > 999) fCandPvLipS2 = 999.; // in this case no 2nd best PV was found, reset to 'good' state
  fCandPvLip12  = fCandPvLip/fpCand->fPvLip2;
  fCandPvLipE12 = fCandPvLipE/fpCand->fPvLipE2;
  fCandPvLipS12 = fCandPvLipS/(fpCand->fPvLip2/fpCand->fPvLipE2);

  // -- 3d impact parameter wrt PV
  if (0) {
    fCandPvIp     = TMath::Sqrt(fCandPvLip*fCandPvLip + fCandPvTip*fCandPvTip); 
    fCandPvIpE    = (fCandPvLip*fCandPvLip*fCandPvLipE*fCandPvLipE + fCandPvTip*fCandPvTip*fCandPvTipE*fCandPvTipE)/(fCandPvIp*fCandPvIp); 
    fCandPvIpE    = TMath::Sqrt(fCandPvIpE); 
    fCandPvIpS    = fCandPvIp/fCandPvIpE;
    if (TMath::IsNaN(fCandPvIpS)) fCandPvIpS = -1.;
  }
  // -- new version directly from CMSSSW and no longer patched...
  fCandPvIp     = fpCand->fPvIP3d;
  fCandPvIpE    = fpCand->fPvIP3dE;
  fCandPvIpS    = fpCand->fPvIP3d/fpCand->fPvIP3dE;
  if (TMath::IsNaN(fCandPvIpS)) fCandPvIpS = -1.;

  fCandPvIp3D   = fpCand->fPvIP3d; 
  fCandPvIpE3D  = fpCand->fPvIP3dE; 
  fCandPvIpS3D  = fpCand->fPvIP3d/fpCand->fPvIP3dE; 
  if (TMath::IsNaN(fCandPvIpS3D)) fCandPvIpS3D = -1.;

  fCandM2 = constrainedMass();

  // -- new variables
  fCandPvDeltaChi2 = fpCand->fDeltaChi2;
  
  TAnaTrack *p0; 
  TAnaTrack *p1(0);
  TAnaTrack *p2(0); 
  
  fCandQ    = 0;

  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    fCandQ += p0->fQ;
    if (TMath::Abs(p0->fMCID) != 13) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }
  
  //cout<<p1<<" "<<p2<<endl;
  //cout<<p1->fMuIndex<<" "<<p2->fMuIndex<<endl;

  // -- for rare backgrounds there are no "true" muons
  if (fpCand->fType > 1000000 && fpCand->fType < 2000000) {
    p1 = fpEvt->getSigTrack(fpCand->fSig1); 
    p2 = fpEvt->getSigTrack(fpCand->fSig2); 
  }

  // Special case for Dstar (prompt and from Bd)
  if (fpCand->fType == 54 || fpCand->fType == 300054 || fpCand->fType == 300031) {
    return;
  }

  // -- Bd2DstarPi
  if (fpCand->fType == 600030 || fpCand->fType == 3000030) {
    if (fpCand->fDau1 < 0 || fpCand->fDau1 > fpEvt->nCands()) return;
    TAnaCand *pD = fpEvt->getCand(fpCand->fDau1); 
    pD = fpEvt->getCand(pD->fDau1); 
    // -- the pion and kaon of the D0 assume the "muon" roles (??)
    p1 = fpEvt->getSigTrack(pD->fSig1); 
    p2 = fpEvt->getSigTrack(pD->fSig2); 
  }


  if (1301 == fpCand->fType ||1302 == fpCand->fType || 1313 == fpCand->fType) {
    // do nothing, these types are for efficiency/TNP studies
  } else {
    if (p1->fPlab.Perp() < p2->fPlab.Perp()) {
      p0 = p1; 
      p1 = p2; 
      p2 = p0; 
    }
  }

  //cout<<p1->fMuIndex<<" "<<p2->fMuIndex<<endl;

  fpMuon1 = p1; 
  fpMuon2 = p2; 
  muScaleCorrectedMasses(); 

  fMu1TrkLayer  = fpReader->numberOfTrackerLayers(p1);
  fMu1TmId      = tightMuon(p1); 
  fMu1MvaId     = mvaMuon(p1, fMu1BDT);        
  fMu1rTmId     = tightMuon(p1, false); 
  fMu1rBDT      = -1.;
  fMu1rMvaId    = mvaMuon(p1, fMu1rBDT, false);
  fTrigMatchDeltaPt = 99.;
  fMu1TrigM     = doTriggerMatchingR(p1, true);
  if (fTrigMatchDeltaPt > 0.1) fMu1TrigM *= -1.;
  fMu1Id        = fMu1MvaId && (fMu1TrigM < 0.1) && (fMu1TrigM > 0); 
  fMu1Id        = fMu1MvaId;
  if (HLTRANGE.begin()->first == "NOTRIGGER") fMu1Id = true; 


  fMu1Pt        = p1->fRefPlab.Perp(); 
  fMu1Eta       = p1->fRefPlab.Eta(); 
  fMu1Phi       = p1->fRefPlab.Phi(); 
  fMu1PtNrf     = p1->fPlab.Perp();
  fMu1EtaNrf    = p1->fPlab.Eta();
  fMu1TkQuality = highPurity(p1); 
  fMu1Q         = p1->fQ;
  fMu1Pix       = fpReader->numberOfPixLayers(p1);
  fMu1BPix      = fpReader->numberOfBPixLayers(p1);
  fMu1BPixL1    = fpReader->numberOfBPixLayer1Hits(p1);
  fMu1PV        = p1->fPvIdx;
  fMu1IP        = p1->fBsTip;
  fMu1IPE       = p1->fBsTipE;

  if (p1->fMuIndex > -1) {
    TAnaMuon *pm = fpEvt->getMuon(p1->fMuIndex);
    fMu1Chi2     = pm->fMuonChi2;
    fMu1Iso      = isoMuon(fpCand, pm); 
    fMu1VtxProb  = pm->fVtxProb;
  } else {
    fMu1Chi2 = -98.;
    fMu1Iso  = -98.; 
    fMu1VtxProb = 99.;
  }

  if (fCandTM && fGenM1Tmi < 0) fpEvt->dump();
  
  TGenCand *pg1(0), *pg2(0);   
  if (fCandTmi > -1) {
    pg1           = fpEvt->getGenTWithIndex(p1->fGenIndex);
    fMu1PtGen     = pg1->fP.Perp();
    fMu1EtaGen    = pg1->fP.Eta();
    fMu1PhiGen    = pg1->fP.Phi();
  } else {
    fMu1PtGen     = -99.;
    fMu1EtaGen    = -99.;
    fMu1PhiGen    = -99.;
  }
  
  //cout<<" check muon 2"<<endl;
  //  fMu2Id        = goodMuon(p2); 
  fMu2TrkLayer  = fpReader->numberOfTrackerLayers(p2);
  fMu2TmId      = tightMuon(p2); 
  fMu2MvaId     = mvaMuon(p2, fMu2BDT);        
  fMu2rTmId     = tightMuon(p2, false); 
  fMu2rBDT      = -1.;
  fMu2rMvaId    = mvaMuon(p2, fMu2rBDT, false);
  fTrigMatchDeltaPt = 99.;
  fMu2TrigM     = doTriggerMatchingR(p2, true);
  if (fTrigMatchDeltaPt > 0.1) fMu2TrigM *= -1.;
  fMu2Id        = fMu2MvaId && (fMu2TrigM < 0.1) && (fMu2TrigM > 0); 
  fMu2Id        = fMu2MvaId;
  if (HLTRANGE.begin()->first == "NOTRIGGER") fMu2Id = true; 

  fMu2Pt        = p2->fRefPlab.Perp(); 
  fMu2Eta       = p2->fRefPlab.Eta(); 
  fMu2Phi       = p2->fRefPlab.Phi(); 
  fMu2PtNrf     = p2->fPlab.Perp();
  fMu2EtaNrf    = p2->fPlab.Eta();
  fMu2TkQuality = highPurity(p2);
  fMu2Q         = p2->fQ;
  fMu2Pix       = fpReader->numberOfPixLayers(p2);
  fMu2BPix      = fpReader->numberOfBPixLayers(p2);
  fMu2BPixL1    = fpReader->numberOfBPixLayer1Hits(p2);
  fMu2PV        = p2->fPvIdx;
  fMu2IP        = p2->fBsTip;
  fMu2IPE       = p2->fBsTipE;

  // -- fill tree for muon id MVA studies
  if (1 && (fMu1rTmId || fMu2rTmId)) {
    TAnaMuon *pt(0);
    for (int i = 0; i < 2; ++i) {
      pt = 0; 
      if (0 == i) {
	if (!fMu1rTmId) continue;
	int idx = p1->fMuIndex;
	if (idx > -1 && idx < fpEvt->nMuons()) {
	  pt = fpEvt->getMuon(idx);
	}
      }
      if (1 == i) {
	if (!fMu2rTmId) continue;
	int idx = p2->fMuIndex;
	if (idx > -1 && idx < fpEvt->nMuons()) {
	  pt = fpEvt->getMuon(idx);
	}
      }

      if (0 == pt) {
	cout << "no TAnaMuon found despite fMu1TmId != 0!!" << endl;
	continue;
      }

      fMuonData.pt            = pt->fPlab.Perp(); 
      fMuonData.eta           = pt->fPlab.Eta(); 

      fMuonData.validMuonHits    = 0; 
      fMuonData.glbNChi2         = pt->fGtrkNormChi2; 
      fMuonData.nMatchedStations = pt->fNmatchedStations; 
      fMuonData.validPixelHits   = fpReader->numberOfPixLayers(pt); // FIXME, kind of correct
      fMuonData.trkLayerWithHits = fpReader->numberOfTrackerLayers(pt);

      fMuonData.trkValidFract = pt->fItrkValidFraction; 
      fMuonData.segComp       = pt->fSegmentComp; 
      fMuonData.chi2LocMom    = pt->fChi2LocalMomentum;
      fMuonData.chi2LocPos    = pt->fChi2LocalPosition;
      fMuonData.glbTrackProb  = pt->fGtrkProb;
      fMuonData.NTrkVHits     = static_cast<float>(pt->fNumberOfValidTrkHits);
      fMuonData.NTrkEHitsOut  = static_cast<float>(pt->fNumberOfLostTrkHits);
      
      fMuonData.kink          = pt->fMuonChi2;

      fMuonData.dpt           = pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag();
      fMuonData.dptrel        = TMath::Abs(pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag())/pt->fInnerPlab.Mag();
      fMuonData.deta          = pt->fInnerPlab.Eta() - pt->fOuterPlab.Eta();
      fMuonData.dphi          = pt->fInnerPlab.DeltaPhi(pt->fOuterPlab);
      fMuonData.dr            = pt->fInnerPlab.DeltaR(pt->fOuterPlab);
      
      fMuonIdTree->Fill();
    }
  }
      


  if (fMu1Id) {
    ((TH1D*)fHistDir->Get("tm_pt"))->Fill(fMu1Pt); 
    ((TH1D*)fHistDir->Get("tm_eta"))->Fill(fMu1Eta); 
  }

  if (fMu2Id) {
    ((TH1D*)fHistDir->Get("tm_pt"))->Fill(fMu2Pt); 
    ((TH1D*)fHistDir->Get("tm_eta"))->Fill(fMu2Eta); 
  }

  if (fMu1MvaId) {
    ((TH1D*)fHistDir->Get("bm_pt"))->Fill(fMu1Pt); 
    ((TH1D*)fHistDir->Get("bm_eta"))->Fill(fMu1Eta); 
  }

  if (fMu2MvaId) {
    ((TH1D*)fHistDir->Get("bm_pt"))->Fill(fMu2Pt); 
    ((TH1D*)fHistDir->Get("bm_eta"))->Fill(fMu2Eta); 
  }

  // -- cut on fMuIndex so that fake muons (from rare backgrounds) can be treated above as real muons
  if (p1->fMuIndex > -1 && p2->fMuIndex > -1) {
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

    if (fVerbose > 10) cout << "dist: " << fMuDist << " dr = " << fMuDeltaR << endl;
  } else {
    fMuDist   = -99.; 
    fMuDeltaR = -99.; 
  }

  if (p2->fMuIndex > -1) {
    TAnaMuon *pm = fpEvt->getMuon(p2->fMuIndex);
    fMu2Chi2      = pm->fMuonChi2;
    fMu2Iso       = isoMuon(fpCand, pm); 
    fMu2VtxProb   = pm->fVtxProb;
  } else {
    fMu2Chi2 = -98.;
    fMu2Iso  = -98.; 
    fMu2VtxProb = 99.;
  }


  xpDistMuons();    

  if ((TMath::Abs(fMu1Eta) < 1.4) && (TMath::Abs(fMu2Eta) < 1.4)) {
    fBarrel = true; 
  } else {
    fBarrel = false; 
  }    

  fChan = detChan(fMu1Eta, fMu2Eta); 

  double dphi = p1->fPlab.DeltaPhi(p2->fPlab);
  fCowboy = (p1->fQ*dphi > 0); 

  // -- Muon weights
  //   PidTable *pT, *pT1, *pT2; 
  //   if (fCowboy) {
  //     pT = fpReader->ptCbMUID; 
  //   } else {
  //     pT = fpReader->ptSgMUID; 
  //   }
  
  //   fMu1W8Mu      = pT->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi);
  //   fMu2W8Mu      = pT->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi);
  
  //   if (fCowboy) {
  //     pT1 = fpReader->ptCbMUT1;
  //     pT2 = fpReader->ptCbMUT2;
  //   } else {
  //     pT1 = fpReader->ptSgMUT1;
  //     pT2 = fpReader->ptSgMUT2;
  //   }
  //   fMu1W8Tr      = pT1->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi)*pT2->effD(fMu1Pt, TMath::Abs(fMu1Eta), fMu1Phi);
  //   fMu2W8Tr      = pT1->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi)*pT2->effD(fMu2Pt, TMath::Abs(fMu2Eta), fMu2Phi);

  fMu1W8Mu = fMu2W8Mu = fMu1W8Tr = fMu2W8Tr = -2.;

  if (fCandTmi > -1) {
    pg2           = fpEvt->getGenTWithIndex(p2->fGenIndex);
    fMu2PtGen     = pg2->fP.Perp();
    fMu2EtaGen    = pg2->fP.Eta();
    fMu2PhiGen    = pg2->fP.Phi(); 
  } else {
    fMu2PtGen     = -99.;
    fMu2EtaGen    = -99.;
    fMu2PhiGen    = -99.;
  }

  fGenMass = -99.;
  if (0 != pg1 && 0 != pg2) {
    TLorentzVector gendimuon = pg1->fP + pg2->fP; 
    fGenMass = gendimuon.M(); 
  }
  //  cout << "m(mu,mu) = " << fGenMass << " n(photons) = " << fNGenPhotons << endl;

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
  TAnaVertex sv2m;
  bool good2m(false); 
  //  cout << "looking at daughters " << fpCand->fDau1  << " .. " << fpCand->fDau2 << endl;
  for (int id = fpCand->fDau1; id <= fpCand->fDau2; ++id) {
    if (id < 0) break;
    pD = fpEvt->getCand(id); 
    //    cout << "looking at daughter " <<  id << " with type = " << pD->fType << endl;
    if (300443 == pD->fType) {
      sv2m = pD->fVtx;
      good2m = true; 
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

  //  virtual double      isoClassicWithDOCA(TAnaCand*, double dca, double r = 1.0, double ptmin = 0.9); 
  double iso = isoClassicWithDOCA(fpCand, 0.05, 0.7, 0.9); // 500um DOCA cut
  fCandIso      = iso; 
  fCandPvTrk    = fCandI0trk;
  fCandIsoTrk   = fCandI2trk;
  pair<int, int> pclose; 
  pclose = nCloseTracks(fpCand, 0.03, 1, 0.5);
  fCandCloseTrk = pclose.first;
  fCandCloseTrkS1 = pclose.second;   
  pclose = nCloseTracks(fpCand, 0.03, 2, 0.5);
  fCandCloseTrkS2 = pclose.second;   
  pclose = nCloseTracks(fpCand, 0.03, 3, 0.5);
  fCandCloseTrkS3 = pclose.second;   

  fCandChi2    = sv.fChi2;
  fCandDof     = sv.fNdof;
  fCandProb    = sv.fProb;
  fCandChi2Dof = fCandChi2/fCandDof;

  fCandOtherVtx = TMath::Max(fMu1VtxProb, fMu2VtxProb) - sv.fProb;

  fCandFL3d  = sv.fD3d;
  fCandFL3dE = sv.fD3dE;
  fCandFLS3d = sv.fD3d/sv.fD3dE; 
  if (TMath::IsNaN(fCandFLS3d)) fCandFLS3d = -1.;
  fCandFLxy  = sv.fDxy;
  fCandFLSxy = sv.fDxy/sv.fDxyE; 
  if (TMath::IsNaN(fCandFLSxy)) fCandFLSxy = -1.;

  fCandTau   = fCandFL3d*MBS/fCandP/TMath::Ccgs();

  // -- variables for production mechanism studies
  //  fpOsCand      = osCand(fpCand);
  fOsIso        = osIsolation(fpCand, 1.0, 0.9); 
  fOsRelIso     = fOsIso/fCandPt; 
  int osm       = osMuon(fpCand, 0.);
  fOsMuonPt     = (osm > -1? fpEvt->getSimpleTrack(osm)->getP().Perp():-1.);
  fOsMuonPtRel  = (osm > -1? fpEvt->getSimpleTrack(osm)->getP().Perp(fpCand->fPlab):-1);
  fOsMuonDeltaR =  (osm > -1? fpCand->fPlab.DeltaR(fpEvt->getSimpleTrack(osm)->getP()):-1.);

//   if (301313 == fCandType) 
//     cout << fRun << " " << fEvt 
// 	 << " from PV = " << fpCand->fPvIdx << " with ntrk = " << fPvNtrk << " av w8 = " << fPvAveW8
// 	 << " cand tracks mindoca = " << fpCand->fMinDoca << " maxdoca = " << fpCand->fMaxDoca << " and chi2/dof = " << fCandChi2/fCandDof 
// 	 << endl;
  
  // -- dimuon vertex version
  if (good2m) { 
    f2MChi2  = sv2m.fChi2;
    f2MDof   = sv2m.fNdof;
    f2MProb  = sv2m.fProb;
    f2MFL3d  = sv2m.fD3d;
    f2MFL3dE = sv2m.fD3dE;
    f2MFLS3d = sv2m.fD3d/sv2m.fD3dE; 
    if (TMath::IsNaN(f2MFLS3d)) f2MFLS3d = -1.;
    f2MFLSxy = sv2m.fDxy/sv2m.fDxyE; 
    if (TMath::IsNaN(f2MFLSxy)) f2MFLSxy = -1.;
  } else {
    f2MChi2  = -1.;
    f2MDof   = -1;
    f2MProb  = -1.;
    f2MFL3d  = -1.;
    f2MFL3dE = -1.;
    f2MFLS3d = -1.; 
    f2MFLSxy = -1.; 
  }

  if (fpCand->fNstTracks.size() == 0) {
    //    cout << "HHHHEEEELLLLPPPP" << endl;
    fCandDocaTrk = 99.;
    fCandDocaTrkBdt = 99.;
  } else {
    fCandDocaTrk    = fpCand->fNstTracks[0].second.first;
    fCandDocaTrkBdt = fpCand->fNstTracks[0].second.first;
    
    int nsize(fpCand->fNstTracks.size());
    TSimpleTrack *ps; 
    for (int i = 0; i<nsize; ++i) {
      int trkId = fpCand->fNstTracks[i].first;
      ps = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((ps->getPvIndex() > -1) && (ps->getPvIndex() != pvidx)) continue;
      
      fCandDocaTrk = fpCand->fNstTracks[i].second.first;
      break;
    }
    
  }

  fillRedTreeData();

  if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX  && fCandIso < 0.7) {
    calcBDT();
  } else {
    calcBDT();
  }  

  fWideMass       = ((fpCand->fMass > MASSMIN) && (fpCand->fMass < MASSMAX)); 

  fGoodMuonsID    = (fMu1Id && fMu2Id);
  fGoodMuonsTmID  = (fMu1TmId && fMu2TmId);
  fGoodMuonsMvaID = (fMu1MvaId && fMu2MvaId); 
  fGoodMuonsPt    = ((fMu1Pt > MUPTLO) && (fMu1Pt < MUPTHI) && (fMu2Pt > MUPTLO) && (fMu2Pt < MUPTHI));
  fGoodMuonsEta   = ((fMu1Eta > MUETALO) && (fMu1Eta < MUETAHI) && (fMu2Eta > MUETALO) && (fMu2Eta < MUETAHI));
  fGoodTracks     = (highPurity(p1) && highPurity(p2));
  fGoodTracksPt   = ((fMu1Pt > TRACKPTLO) && (fMu1Pt < TRACKPTHI) && (fMu2Pt > TRACKPTLO) && (fMu2Pt < TRACKPTHI));
  fGoodTracksEta  = ((fMu1Eta > TRACKETALO) && (fMu1Eta < TRACKETAHI) && (fMu2Eta > TRACKETALO) && (fMu2Eta < TRACKETAHI));

  fGoodQ          = (fMu1Q*fMu2Q < 0); 
  fGoodPvAveW8    = (fPvAveW8 > PVAVEW8);
  fGoodPvLip      = (TMath::Abs(fCandPvLip) < CANDLIP); 
  fGoodPvLipS     = (TMath::Abs(fCandPvLipS) < CANDLIPS); 
  fGoodPvLip2     = (TMath::Abs(fCandPvLip2) > CANDLIP2); 
  fGoodPvLipS2    = (TMath::Abs(fCandPvLipS2) > CANDLIPS2); 
  fGoodMaxDoca    = (TMath::Abs(fCandDoca) < CANDDOCA); 
  fGoodIp         = (TMath::Abs(fCandPvIp) < CANDIP); 
  fGoodIpS        = (TMath::Abs(fCandPvIpS) < CANDIPS); 
    
  fGoodPt         = (fCandPt > CANDPTLO);
  fGoodEta        = ((fCandEta > CANDETALO) && (fCandEta < CANDETAHI)); 
  fGoodAlpha      = (fCandA < CANDALPHA); 
  fGoodChi2       = (fCandChi2/fCandDof < CANDVTXCHI2);
  fGoodFLS        =  ((fCandFLS3d > CANDFLS3D) && (fCandFLSxy > CANDFLSXY)); 
  if (TMath::IsNaN(fCandFLS3d)) fGoodFLS = false;

  fGoodCloseTrack = (fCandCloseTrk < CANDCLOSETRK); 
  fGoodIso        = (fCandIso > CANDISOLATION); 
  fGoodDocaTrk    = (fCandDocaTrk > CANDDOCATRK);
  fGoodLastCut    = true; 

  fAnaCuts.update(); 

  //   fPreselection = fWideMass && fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsPt && fGoodMuonsEta && fGoodQ;  
  //   fPreselection = fPreselection && (fCandPt > 5) && (fCandA < 0.3) && (fCandChi2/fCandDof < 10) && (fCandFL3d < 2); 

  // -- to be consistent with the BDT traning
  ((TH1D*)fHistDir->Get("test3"))->Fill(1.); 

  fPreselection = preselection(fRTD, fChan);
  //fPreselection = preselection(fRTD, fChan, ((TH1D*)fHistDir->Get("test3")) );
  if(fPreselection) ((TH1D*)fHistDir->Get("test3"))->Fill(2.); 

  fPreselection = fPreselection && fGoodHLT;
  if(fPreselection) ((TH1D*)fHistDir->Get("test3"))->Fill(3.); 

  //  fPreselection = true; 

}

// ----------------------------------------------------------------------
void candAna::fillCandidateHistograms(int offset) {
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
  fAnaCuts.addCut("fGoodPvAveW8", "<w8>", fGoodPvAveW8); 
  fAnaCuts.addCut("fGoodPvLip", "LIP(PV)", fGoodPvLip); 
  fAnaCuts.addCut("fGoodPvLipS", "LIPS(PV)", fGoodPvLipS); 
  fAnaCuts.addCut("fGoodPvLip2", "LIP2(PV)", fGoodPvLip2); 
  fAnaCuts.addCut("fGoodPvLipS2", "LIPS2(PV)", fGoodPvLipS2); 
  fAnaCuts.addCut("fGoodIp", "IP", fGoodIp); 
  fAnaCuts.addCut("fGoodIpS", "IPS", fGoodIpS); 
  fAnaCuts.addCut("fGoodMaxDoca", "MAXDOCA", fGoodMaxDoca); 
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt); 
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta); 
  fAnaCuts.addCut("fGoodAlpha", "#alpha", fGoodAlpha); 
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS); 
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2); 
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso); 
  fAnaCuts.addCut("fGoodCloseTrack", "close track veto", fGoodCloseTrack); 
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk); 
  fAnaCuts.addCut("fGoodLastCut", "lastCut", fGoodLastCut); 
}


// ----------------------------------------------------------------------
void candAna::moreCandidateCuts() {
  cout << "    candAna more candidate cuts?" << endl;

}



// ----------------------------------------------------------------------
bool candAna::highPurity(TAnaTrack *pt) {

  if (pt->fTrackQuality & 4) return true; 
  return false;

//   if (TRACKQUALITY > 0 && (0 == (pt->fTrackQuality & TRACKQUALITY))) {
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed track quality: " << pt->fTrackQuality << endl;
//     return false; 
//   }
  
//   if (TMath::Abs(pt->fTip) > TRACKTIP) {
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed tip: " << pt->fTip
// 			   << " pointing to PV = "  << pt->fPvIdx  << endl;
//     return false; 
//   }
  
//   if (TMath::Abs(pt->fLip) > TRACKLIP) { 
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed lip: " << pt->fLip 
// 			   << " pointing to PV = "  << pt->fPvIdx  << endl;          
//     return false; 
//   }

  return true; 
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
void candAna::efficiencyCalculation() {
  cout << "candAna::efficiencyCalculation()  wrong function" << endl;
}

// ----------------------------------------------------------------------
void candAna::triggerSelection() {

  fGoodHLT = false; 
  fhltType = 0; 

  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 

  if (0) {
    cout << " ----------------------------------------------------------------------" << endl;
    TAnaMuon *pM(0);
    for (int i = 0; i < fpEvt->nMuons(); ++i) {
      pM = fpEvt->getMuon(i);
      if (pM->fMuID & 1) {
	cout << "STA pt = " << pM->fPlab.Pt() << " eta = " << pM->fPlab.Eta() << " phi = " << pM->fPlab.Phi() << endl;
      }
    }

    TTrgObj *p; 
    
    for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
      p = fpEvt->getTrgObj(i); 
      p->dump();
    }
  }


  if (0) {
    TTrgObj *p; 
    
    for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
      p = fpEvt->getTrgObj(i); 
      //      cout << p->fLabel << endl;
      //      if (!p->fLabel.CompareTo("hltL1sL1DoubleMu33HighQ:HLT::")) cout << "= " << p->fLabel << endl;
      if (!p->fLabel.CompareTo("hltL1sL1DoubleMuOpen:HLT::")) {
        if (0) cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
      }
      
      if (1 && !p->fLabel.CompareTo("hltL1sL1DoubleMu33HighQ:HLT::")) {
        //cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_0"))->Fill(p->fP.Eta());
      }
      
      if (1 && !p->fLabel.CompareTo("hltL1sL1DoubleMu0or33HighQ:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_1"))->Fill(p->fP.Eta());
      }

      if (1 && !p->fLabel.CompareTo("hltDimuon33L1Filtered0:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_2"))->Fill(p->fP.Eta());
      }

      if (1 && !p->fLabel.CompareTo("hltDimuonL2PreFiltered0:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_3"))->Fill(p->fP.Eta());
      }

      if (1 && !p->fLabel.CompareTo("hltDimuonL2PreFiltered0:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_4"))->Fill(p->fP.Eta());
      }

      if (1 && !p->fLabel.CompareTo("hltDoubleDisplacedMu4L3PreFiltered:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_5"))->Fill(p->fP.Eta());
      }

      if (1 && !p->fLabel.CompareTo("hltDoubleMu3p5LowMassDisplacedL3Filtered:HLT::")) {
        //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
        ((TH1D*)fHistDir->Get("L1_6"))->Fill(p->fP.Eta());
      }
    }
    
  }


  if ( HLTRANGE.begin()->first == "NOTRIGGER" ) { // keep NOTRIGGER here and place ALLTRIGGER at the end.
    //    cout << "NOTRIGGER requested... " << endl;
    fGoodHLT = true; 
    return;
  }

  if (fVerbose == -32) {
    cout << "--------------------  L1" << endl;
    for (int i = 0; i < NL1T; ++i) {
      result = wasRun = error = false;
      a = fpEvt->fL1TNames[i]; 
      ps = fpEvt->fL1TPrescale[i]; 
      result = fpEvt->fL1TResult[i]; 
      error  = fpEvt->fL1TMask[i]; 
      //if (a.Contains("Mu")) {
      if (result ) {
	cout << a <<  " mask: " << error << " result: " << result << " ps: " << ps << endl;
      }
    }
  }
  
  ((TH1D*)fHistDir->Get("test9"))->Fill(0.);
  fhltType=0; // reset this variable 
  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 
    
    //    if (-32 == fVerbose) cout << "path " << i << ": " << a << endl;
    if (wasRun && result) {
      if ((a != "digitisation_step") 
	  && (a != "L1simulation_step") 
	  && (a != "digi2raw_step") 
	  && (a != "HLTriggerFinalPath") 
	  && (a != "raw2digi_step") 
	  && (a != "reconstruction_step") 
	  ) {

	if (-32 == fVerbose) cout << "run and fired: " << a << endl;
	
	// special hlt test, check and store which trigger has fired 
        if (a.Contains("Mu") || a.Contains("mu")) {
          fhltType = fhltType | 0x10;
          if (fVerbose > 2) cout << a << " " << fhltType<<" "<<wasRun << " "<< result <<  endl;
          ((TH1D*)fHistDir->Get("test9"))->Fill(9.);
	  if (a.Contains("L3") ) {
	    fhltType = fhltType | 0x20;
	    ((TH1D*)fHistDir->Get("test9"))->Fill(8.);
	  }
        }

	if (a.Contains("HLT_DoubleMu3_4_Dimuon5_Bs_Central")) {
          fhltType = fhltType | 0x1;
          //cout << a << " " << fhltType<<" "<<wasRun << " "<< result <<  endl;
          ((TH1D*)fHistDir->Get("test9"))->Fill(1.);
        }
        if (a.Contains("HLT_DoubleMu3p5_4_Dimuon5_Bs_Central")) {
          fhltType = fhltType | 0x2;
          //cout << a << " " << fhltType<<" "<< wasRun << " "<< result <<  endl;
          ((TH1D*)fHistDir->Get("test9"))->Fill(2.);
        } else if (a.Contains("HLT_DoubleMu4_Dimuon7_Bs_Forward")) {
          fhltType = fhltType | 0x4;
          //cout << a << " " << fhltType<<" "<< wasRun << " "<< result <<  endl;
          ((TH1D*)fHistDir->Get("test9"))->Fill(3.);
        } else if (a.Contains("HLT_DoubleMu4_Jpsi_Displaced")) {
          fhltType = fhltType | 0x8;
          //cout << a << " " << fhltType<<" "<< wasRun << " "<< result <<  endl;
          ((TH1D*)fHistDir->Get("test9"))->Fill(4.);
        } else if (a.Contains("LowMass")) {
          ((TH1D*)fHistDir->Get("test9"))->Fill(5.);
        } else if (a.Contains("SameSign")) {
          ((TH1D*)fHistDir->Get("test9"))->Fill(6.);
        } else if (a.Contains("Upsilon")) {
          ((TH1D*)fHistDir->Get("test9"))->Fill(7.);
        }
      } // if step

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
	  fHLTPath = a.Data();
	  if (fVerbose > 1) cout << "exact match: " << imap->first.c_str() << " HLT: " << a << " result: " << result << endl;
	}

 	if (a.Contains(spath.c_str()) && (rmin <= fRun) && (fRun <= rmax)) {
	  fHLTPath = a.Data();
 	  fGoodHLT = true; 
	  if (fVerbose > 1) cout << "close match: " << imap->first.c_str() << " HLT: " << a 
				 << " result: " << result 
				 << " in run " << fRun 
				 << endl;
 	}
      }
    }      
  }

  // I prefer to have it at the end. So the trigger informatin is processed correctly before it is overwritten
  if ( HLTRANGE.begin()->first == "ALLTRIGGER") { // extend to ALLTRIGGER 16/1/13 d.k.
    if(fVerbose >1) cout << "ALLTRIGGER requested... " << endl;
    // For data make always true, the event has been triggered anyway
    if(!fIsMC) {fGoodHLT = true;} 
    else { // for MC check if there was any HLT firing

      for (int i = 0; i < NHLT; ++i) {
	result = wasRun = false;
	//a = fpEvt->fHLTNames[i]; 
	//ps = fpEvt->fHLTPrescale[i]; 
	wasRun = fpEvt->fHLTWasRun[i]; 
	result = fpEvt->fHLTResult[i]; 
	//error  = fpEvt->fHLTError[i]; 
	if (wasRun && result) {fGoodHLT=true; break;}
      } // end HLT for loop 

    } // if-else

  } // ALLTRIGGER

  // Test trigger object for trigger matching TESTING ONLY
  if (true && fGoodHLT) {    
    ((TH1D*)fHistDir->Get("test10"))->Fill(0.);
    bool passed[20];
    for(int i=0;i<20;i++) passed[i]=false;
    TTrgObj *p;     
    for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
      p = fpEvt->getTrgObj(i); 

      //cout<<i<<"  "<< p->fLabel << " number " << p->fNumber <<" ID = " << p->fID << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() << " "<<p->fID << endl;

      if(false && p->fNumber > -1) {
	cout<<" Selected object: label "<< p->fLabel << " number " << p->fNumber <<" ID = " << p->fID << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
      }

      if ( p->fLabel.Contains("mu") || p->fLabel.Contains("Mu") ) {passed[0]=true; if( p->fLabel.Contains("L3")) passed[10]=true;}
      if ( p->fLabel.Contains("Displaced") ) {passed[1]=true; if( p->fLabel.Contains("L3")) passed[11]=true;}
      if ( p->fLabel.Contains("Jpsi") )      {passed[2]=true; if( p->fLabel.Contains("L3")) passed[12]=true;}
      if ( p->fLabel.Contains("Bs") )        {passed[3]=true; if( p->fLabel.Contains("L3")) passed[13]=true;}
      if ( p->fLabel.Contains("LowMass") )   {passed[4]=true; if( p->fLabel.Contains("L3")) passed[14]=true;}
      if ( p->fLabel.Contains("Vertex") )    {passed[5]=true; if( p->fLabel.Contains("L3")) passed[15]=true;}
      if ( p->fLabel.Contains("SameSign") )  {passed[6]=true; if( p->fLabel.Contains("L3")) passed[16]=true;}
      if ( p->fNumber > -1 )  {passed[7]=true; if( p->fLabel.Contains("L3")) passed[17]=true;}
      
    } // for look 
    bool passedFinal=false, passedL3=false;
    for(int i=0;i<20;i++) {
      if(passed[i]) {
	passedFinal=true;
	((TH1D*)fHistDir->Get("test10"))->Fill(float(i+1));
	if(i>=10) passedL3=true; // count only with L3
      }
    }
    if(passedFinal) ((TH1D*)fHistDir->Get("test10"))->Fill(9.);
    if(passedL3) ((TH1D*)fHistDir->Get("test10"))->Fill(10.);

  } // if true

  if (false == fGoodHLT) {
    if (fVerbose > 1) cout << "------->  event NOT triggered!" << endl;
  }

}

// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  TH1D *h11(0); 
  (void)h11; 
  h11 = new TH1D("tm_pt", "tight muon pT", 50, 0., 25.); 
  h11 = new TH1D("bm_pt", "BDT muon pT", 50, 0., 25.); 

  h11 = new TH1D("tm_eta", "tight muon eta", 50, -2.5, 2.5); 
  h11 = new TH1D("bm_eta", "BDT muon eta", 50, -2.5, 2.5); 


  h11 = new TH1D("L1_0", "hltL1sL1DoubleMu33HighQ", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_1", "hltL1sL1DoubleMu0or33HighQ", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_2", "hltDimuon33L1Filtered0", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_3", "hltDimuon33L2PreFiltered0", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_4", "hltDimuonL2PreFiltered0", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_5", "hltDoubleDisplacedMu4L3PreFiltered", 50, -2.5, 2.5); 
  h11 = new TH1D("L1_6", "hltDoubleMu3p5LowMassDisplacedL3Filtered", 50, -2.5, 2.5); 

  const double firstRun = 190400.5, lastRun = 210400.5; const int numRuns = 20000;
  h11 = new TH1D("run1", "runs all", numRuns, firstRun, lastRun); 
  h11 = new TH1D("run2", "runs json", numRuns, firstRun, lastRun); 
  h11 = new TH1D("run3", "runs json&hlt", numRuns, firstRun, lastRun); 
  h11 = new TH1D("run4", "runs json&hltmatched", numRuns, firstRun, lastRun); 
  //h11 = new TH1D("run5", "runs jpsik", numRuns, firstRun, lastRun); 
  //h11 = new TH1D("run6", "runs jpsiphi", numRuns, firstRun, lastRun); 
  //h11 = new TH1D("run7", "runs mumu", numRuns, firstRun, lastRun); 
  h11 = new TH1D("test1", "test1",50, 0., 50.); 
  h11 = new TH1D("test2", "test2",100, 0., 0.2); 
  h11 = new TH1D("test3", "test3",50, 0., 50.); 
  h11 = new TH1D("test4", "test4",100, 0., 0.2); 
  h11 = new TH1D("test5", "test5",100, 0., 0.2); 
  h11 = new TH1D("test6", "test6",100, 0., 0.2); 
  h11 = new TH1D("test7", "test7", 400, 0., 4.);
  h11 = new TH1D("test8", "test8",1000, 0., 2.); 
  h11 = new TH1D("test9", "test9",10, -0.5, 9.5); 
  h11 = new TH1D("test10", "test10",20, -0.5, 19.5); 


  h11 = new TH1D("gp1cms", "p1cms", 50, 0, 10.); 
  h11 = new TH1D("gp2cms", "p2cms", 50, 0, 10.); 
  h11 = new TH1D("gt1cms", "t1cms", 50, -1, 1.); 
  h11 = new TH1D("gt2cms", "t2cms", 50, -1, 1.); 

  h11 = new TH1D("gp1cmsg", "p1cms (with photons)", 50, 0, 10.); 
  h11 = new TH1D("gp2cmsg", "p2cms (with photons)", 50, 0, 10.); 
  h11 = new TH1D("gt1cmsg", "t1cms (with photons)", 50, 0, 1.); 
  h11 = new TH1D("gt2cmsg", "t2cms (with photons)", 50, 0, 1.); 

  h11 = new TH1D("rp1cms", "p1cms", 50, 0, 10.); 
  h11 = new TH1D("rp2cms", "p2cms", 50, 0, 10.); 
  h11 = new TH1D("rt1cms", "t1cms", 50, -1, 1.); 
  h11 = new TH1D("rt2cms", "t2cms", 50, -1, 1.); 
  h11 = new TH1D("rt3cms", "t3cms", 50, -1, 1.); 
  TH2D *h22(0); 
  (void)h22; 
  h22 = new TH2D("tvsm",   "tvsm", 50, 4.9, 5.9, 50, -1., 1.); 

  h11 = new TH1D("rp1cmsg", "p1cms (with photons)", 50, 0, 10.); 
  h11 = new TH1D("rp2cmsg", "p2cms (with photons)", 50, 0, 10.); 
  h11 = new TH1D("rt1cmsg", "t1cms (with photons)", 50, 0, 1.); 
  h11 = new TH1D("rt2cmsg", "t2cms (with photons)", 50, 0, 1.); 
  h11 = new TH1D("rt3cmsg", "t2cms (with photons)", 50, 0, 1.); 

  h11 = new TH1D("gt1", "gt1", 50, -1, 1.);   
  h11 = new TH1D("gt2", "gt2", 50, -1, 1.); 

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree(fTree); 
  
  // -- tree for AMS
  fAmsTree = new TTree("amsevents", "amsevents");
  setupReducedTree(fAmsTree); 

  // -- tree for AMS
  fMuonIdTree = new TTree("muonidtree", "muonidtree");
  setupMuonIdTree(fMuonIdTree); 

  // -- Efficiency/Acceptance Tree
  fEffTree = new TTree("effTree", "effTree");
  fEffTree->Branch("run",    &fRun,               "run/L");
  fEffTree->Branch("evt",    &fEvt,               "evt/L");
  fEffTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fEffTree->Branch("procid", &fProcessType,       "procid/I");
  fEffTree->Branch("bidx",   &fGenBTmi,           "bidx/I");

  fEffTree->Branch("gpt",    &fETgpt,             "gpt/F");
  fEffTree->Branch("geta",   &fETgeta,            "geta/F");

  fEffTree->Branch("m1pt",   &fETm1pt,            "m1pt/F");
  fEffTree->Branch("g1pt",   &fETg1pt,            "g1pt/F");
  fEffTree->Branch("m1eta",  &fETm1eta,           "m1eta/F");
  fEffTree->Branch("g1eta",  &fETg1eta,           "g1eta/F");
  fEffTree->Branch("m1q",    &fETm1q,             "m1q/I");
  fEffTree->Branch("m1gt",   &fETm1gt,            "m1gt/O");
  fEffTree->Branch("m1id",   &fETm1id,            "m1id/O");
  fEffTree->Branch("m1tmid", &fETm1tmid,          "m1tmid/O");
  fEffTree->Branch("m1mvaid",&fETm1mvaid,         "m1mvaid/O");

  fEffTree->Branch("m2pt",   &fETm2pt,            "m2pt/F");
  fEffTree->Branch("g2pt",   &fETg2pt,            "g2pt/F");
  fEffTree->Branch("m2eta",  &fETm2eta,           "m2eta/F");
  fEffTree->Branch("g2eta",  &fETg2eta,           "g2eta/F");
  fEffTree->Branch("m2q",    &fETm2q,             "m2q/I");
  fEffTree->Branch("m2gt",   &fETm2gt,            "m2gt/O");
  fEffTree->Branch("m2id",   &fETm2id,            "m2id/O");
  fEffTree->Branch("m2tmid", &fETm2tmid,          "m2tmid/O");
  fEffTree->Branch("m2mvaid",&fETm2mvaid,         "m2mvaid/O");

  fEffTree->Branch("m",      &fETcandMass,        "m/F");


  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.); 
  (void)h;

}


// ----------------------------------------------------------------------
void candAna::setupMuonIdTree(TTree *t) {

  t->Branch("pt",                   &fMuonData.pt,  "pt/F");
  t->Branch("eta",                  &fMuonData.eta, "eta/F");

  t->Branch("intvalidmuonhits",     &fMuonData.validMuonHits, "validmuonhits/I");
  t->Branch("intnmatchedstations",  &fMuonData.nMatchedStations, "nmatchedstations/I");
  t->Branch("intvalidpixelhits",    &fMuonData.validPixelHits, "validpixelhits/I");
  t->Branch("inttrklayerswithhits", &fMuonData.trkLayerWithHits, "trklayerswithhits/I");

  t->Branch("gchi2",                &fMuonData.glbNChi2, "gchi2/F");

  t->Branch("itrkvalidfraction",    &fMuonData.trkValidFract, "itrkvalidfraction/F");
  t->Branch("segcomp",              &fMuonData.segComp, "segcomp/F");
  t->Branch("chi2lmom",             &fMuonData.chi2LocMom, "chi2lmom/F");
  t->Branch("chi2lpos",             &fMuonData.chi2LocPos, "chi2lpos/F");
  t->Branch("gtrkprob",             &fMuonData.glbTrackProb, "gtrkprob/F");
  t->Branch("ntrkvhits",            &fMuonData.NTrkVHits, "ntrkvhits/F");
  t->Branch("ntrkehitsout",         &fMuonData.NTrkEHitsOut, "ntrkehitsout/F");

  t->Branch("kink",                 &fMuonData.kink, "kink/F");
  t->Branch("dpt",                  &fMuonData.dpt, "dpt/F");
  t->Branch("dptrel",               &fMuonData.dptrel, "dptrel/F");
  t->Branch("deta",                 &fMuonData.deta, "deta/F");
  t->Branch("dphi",                 &fMuonData.dphi, "dphi/F");
  t->Branch("dr",                   &fMuonData.dr, "dr/F");

}

// ----------------------------------------------------------------------
void candAna::setupReducedTree(TTree *t) {

  t->Branch("run",     &fRun,               "run/L");
  t->Branch("json",    &fJSON,              "json/O");
  t->Branch("evt",     &fEvt,               "evt/L");
  t->Branch("ls",      &fLS,                "ls/I");
  t->Branch("tm",      &fCandTM,            "tm/I");
  t->Branch("pr",      &fGenBpartial,       "pr/I"); 
  t->Branch("procid",  &fProcessType,       "procid/I");
  t->Branch("hlt",     &fGoodHLT,           "hlt/O");
  t->Branch("pvn",     &fPvN,               "pvn/I");
  t->Branch("cb",      &fCowboy,            "cb/O");
  t->Branch("rr",      &fRunRange,          "rr/I");
  t->Branch("bdt",     &fBDT,               "bdt/D");
  t->Branch("npv",     &fPvN,               "npv/I");
  t->Branch("pvw8",    &fPvAveW8,           "pvw8/D");

  // -- global cuts and weights
  t->Branch("gmuid",   &fGoodMuonsID,       "gmuid/O");
  t->Branch("gmutmid", &fGoodMuonsTmID,     "gmutmid/O");
  t->Branch("gmumvaid",&fGoodMuonsMvaID,    "gmumvaid/O");
  t->Branch("gmupt",   &fGoodMuonsPt,       "gmupt/O");
  t->Branch("gmueta",  &fGoodMuonsEta,      "gmueta/O");
  t->Branch("gtqual",  &fGoodTracks,        "gtqual/O");
  t->Branch("gtpt",    &fGoodTracksPt,      "gtpt/O");
  t->Branch("gteta",   &fGoodTracksEta,     "gteta/O");

  // -- PV
  t->Branch("pvlip",    &fCandPvLip,        "pvlip/D");
  t->Branch("pvlips",   &fCandPvLipS,       "pvlips/D");
  t->Branch("pvlip2",   &fCandPvLip2,       "pvlip2/D");
  t->Branch("pvlips2",  &fCandPvLipS2,      "pvlips2/D");
  t->Branch("pvip",     &fCandPvIp,         "pvip/D");
  t->Branch("pvips",    &fCandPvIpS,        "pvips/D");
  t->Branch("pvip3d",   &fCandPvIp3D,       "pvip3d/D");
  t->Branch("pvips3d",  &fCandPvIpS3D,      "pvips3d/D");

  // -- cand
  t->Branch("q",       &fCandQ,             "q/I");
  t->Branch("type",    &fCandType,          "type/I");
  t->Branch("pt",      &fCandPt,            "pt/D");
  t->Branch("eta",     &fCandEta,           "eta/D");
  t->Branch("phi",     &fCandPhi,           "phi/D");
  t->Branch("tau",     &fCandTau,           "tau/D");
  t->Branch("m",       &fCandM,             "m/D");
  t->Branch("me",      &fCandME,            "me/D");
  t->Branch("cm",      &fCandM2,            "cm/D");
  t->Branch("m3",      &fCandM3,            "m3/D");
  t->Branch("m4",      &fCandM4,            "m4/D");
  t->Branch("cosa",    &fCandCosA,          "cosa/D");
  t->Branch("alpha",   &fCandA,             "alpha/D");
  t->Branch("iso",     &fCandIso,           "iso/D");
  t->Branch("isotrk",  &fCandIsoTrk,        "isotrk/I");
  t->Branch("closetrk",&fCandCloseTrk,      "closetrk/I");
  t->Branch("chi2",    &fCandChi2,          "chi2/D");
  t->Branch("dof",     &fCandDof,           "dof/D");
  t->Branch("chi2dof", &fCandChi2Dof,       "chi2dof/D");
  t->Branch("prob",    &fCandProb,          "prob/D");
  t->Branch("fls3d",   &fCandFLS3d,         "fls3d/D");
  t->Branch("fl3d",    &fCandFL3d,          "fl3d/D");
  t->Branch("flxy",    &fCandFLxy,          "flxy/D");
  t->Branch("fl3dE",   &fCandFL3dE,         "fl3dE/D");
  t->Branch("flsxy",   &fCandFLSxy,         "flsxy/D");
  t->Branch("docatrk", &fCandDocaTrk,       "docatrk/D");  
  t->Branch("docatrkbdt", &fCandDocaTrkBdt, "docatrkbdt/D");  
  t->Branch("maxdoca", &fCandDoca,          "maxdoca/D");
  t->Branch("lip",     &fCandPvLip,         "lip/D");
  t->Branch("lipE",    &fCandPvLipE,        "lipE/D");
  t->Branch("tip",     &fCandPvTip,         "tip/D");
  t->Branch("tipE",    &fCandPvTipE,        "tipE/D");

  t->Branch("pvdchi2",   &fCandPvDeltaChi2, "pvdchi2/D");
  t->Branch("closetrks1", &fCandCloseTrkS1,   "closetrks1/I");
  t->Branch("closetrks2", &fCandCloseTrkS2,   "closetrks2/I");
  t->Branch("closetrks3", &fCandCloseTrkS3,   "closetrks3/I");
  t->Branch("othervtx", &fCandOtherVtx, "othervtx/D"); 

  t->Branch("osiso",   &fOsIso,             "osiso/D");
  t->Branch("osreliso",&fOsRelIso,          "osreliso/D");
  t->Branch("osmpt",   &fOsMuonPt,          "osmpt/D");
  t->Branch("osmptrel",&fOsMuonPtRel,       "osmptrel/D");
  t->Branch("osmdr",   &fOsMuonDeltaR,      "osmdr/D");

  // -- muons
  t->Branch("m1q",     &fMu1Q,              "m1q/I");
  t->Branch("m1id",    &fMu1Id,             "m1id/O");
  t->Branch("m1tmid",  &fMu1TmId,           "m1tmid/O");
  t->Branch("m1mvaid", &fMu1MvaId,          "m1mvaid/O");
  t->Branch("m1rtmid", &fMu1rTmId,          "m1rtmid/O");
  t->Branch("m1rmvaid",&fMu1rMvaId,          "m1rmvaid/O");
  t->Branch("m1mvabdt",&fMu1BDT,            "m1mvabdt/D");
  t->Branch("m1rmvabdt",&fMu1rBDT,          "m1rmvabdt/D");
  t->Branch("m1trigm", &fMu1TrigM,          "m1trigm/D");

  t->Branch("m1pt",    &fMu1Pt,             "m1pt/D");
  t->Branch("m1eta",   &fMu1Eta,            "m1eta/D");
  t->Branch("m1phi",   &fMu1Phi,            "m1phi/D");
  t->Branch("m1ip",    &fMu1IP,             "m1ip/D");
  t->Branch("m1gt",    &fMu1TkQuality,      "m1gt/I");
  t->Branch("m1pix",   &fMu1Pix,            "m1pix/I");
  t->Branch("m1bpix",  &fMu1BPix,           "m1bpix/I");
  t->Branch("m1bpixl1",&fMu1BPixL1,         "m1bpixl1/I");
  t->Branch("m1chi2",  &fMu1Chi2,           "m1chi2/D");
  t->Branch("m1pv",    &fMu1PV,             "m1pv/I");
  t->Branch("m1vtxprob",&fMu1VtxProb,       "m1vtxprob/D");
  t->Branch("m1xpdist",&fMu1XpDist,         "m1xpdist/D");
  t->Branch("m1iso",   &fMu1Iso,            "m1iso/D"); 
  
  t->Branch("m2q",     &fMu2Q,              "m2q/I");
  t->Branch("m2id",    &fMu2Id,             "m2id/O");
  t->Branch("m2tmid",  &fMu2TmId,           "m2tmid/O");
  t->Branch("m2mvaid", &fMu2MvaId,          "m2mvaid/O");
  t->Branch("m2rtmid", &fMu2rTmId,          "m2rtmid/O");
  t->Branch("m2rmvaid",&fMu2rMvaId,         "m2rmvaid/O");
  t->Branch("m2mvabdt",&fMu2BDT,            "m2mvabdt/D");
  t->Branch("m2rmvabdt",&fMu2rBDT,          "m2rmvabdt/D");
  t->Branch("m2trigm", &fMu2TrigM,          "m2trigm/D");
  
  t->Branch("m2pt",    &fMu2Pt,             "m2pt/D");
  t->Branch("m2eta",   &fMu2Eta,            "m2eta/D");
  t->Branch("m2phi",   &fMu2Phi,            "m2phi/D");
  t->Branch("m2ip",    &fMu2IP,             "m2ip/D");
  t->Branch("m2gt",    &fMu2TkQuality,      "m2gt/I");
  t->Branch("m2pix",   &fMu2Pix,            "m2pix/I");
  t->Branch("m2bpix",  &fMu2BPix,           "m2bpix/I");
  t->Branch("m2bpixl1",&fMu2BPixL1,         "m2bpixl1/I");
  t->Branch("m2chi2",  &fMu2Chi2,           "m2chi2/D");
  t->Branch("m2pv",    &fMu2PV,             "m2pv/I");
  t->Branch("m2vtxprob",&fMu2VtxProb,       "m2vtxprob/D");
  t->Branch("m2xpdist",&fMu2XpDist,         "m2xpdist/D");
  t->Branch("m2iso",   &fMu2Iso,            "m2iso/D"); 

  t->Branch("mudist",  &fMuDist,            "mudist/D");
  t->Branch("mudeltar",&fMuDeltaR,          "mudeltar/D");
  t->Branch("hltm",    &fHLTmatch,          "hltm/O");
  t->Branch("hltt",    &fhltType,           "hltt/I");

  t->Branch("g1pt",    &fMu1PtGen,          "g1pt/D");
  t->Branch("g2pt",    &fMu2PtGen,          "g2pt/D");
  t->Branch("g1eta",   &fMu1EtaGen,         "g1eta/D");
  t->Branch("g2eta",   &fMu2EtaGen,         "g2eta/D");
  t->Branch("g1phi",   &fMu1PhiGen,         "g1phi/D");
  t->Branch("g2phi",   &fMu2PhiGen,         "g2phi/D");
  t->Branch("gmass",   &fGenMass,           "gmass/D");
  t->Branch("gtau",    &fGenLifeTime,       "gtau/D");
  t->Branch("g1id",    &fMu1GenID,          "g1id/I");
  t->Branch("g2id",    &fMu2GenID,          "g2id/I");

  t->Branch("t1pt",    &fMu1PtNrf,          "t1pt/D");
  t->Branch("t1eta",   &fMu1EtaNrf,         "t1eta/D");
  t->Branch("t2pt",    &fMu2PtNrf,          "t2pt/D");
  t->Branch("t2eta",   &fMu2EtaNrf,         "t2eta/D");

  t->Branch("hm1pt",  &fHltMu1Pt,  "hm1pt/D");    
  t->Branch("hm1eta", &fHltMu1Eta, "hm1eta/D");  
  t->Branch("hm1phi", &fHltMu1Phi, "hm1phi/D");  
  t->Branch("hm2pt",  &fHltMu2Pt,  "hm2pt/D");    
  t->Branch("hm2eta", &fHltMu2Eta, "hm2eta/D");  
  t->Branch("hm2phi", &fHltMu2Phi, "hm2phi/D");  

  // for testing only 
  t->Branch("hltm2",    &fHLTmatch2,          "hltm2/O");

}


// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  // -- define default values for some cuts
  NOPRESELECTION = 0; 
  IGNORETRIGGER  = 0; 

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

  char  buffer[1000], XmlName[1000];
  fHistDir->cd();
  if (dump) cout << "gDirectory: "; fHistDir->pwd();
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); 
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
      SELMODE = int(CutValue); 
      if (dump) cout << "SELMODE:           " << SELMODE << endl;
      ibin = 2;
      hcuts->SetBinContent(ibin, SELMODE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Selection Mode :: %i", CutName, SELMODE));
    }

    if (!strcmp(CutName, "TRIGRANGE")) {
      char triggerlist[1000]; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      int r1(0), r2(0); 
      string hlt = splitTrigRange(tl, r1, r2); 
      HLTRANGE.insert(make_pair(hlt, make_pair(r1, r2))); 
      if (dump) {
	cout << "HLTRANGE:       " << hlt << " from " << r1 << " to " << r2 << endl; 
      }
      ibin = 3; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }


    if (!strcmp(CutName, "TRUTHCAND")) {
      TRUTHCAND = int(CutValue); 
      if (dump) cout << "TRUTHCAND:           " << TRUTHCAND << endl;
      ibin = 4;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "IGNORETRIGGER")) {
      IGNORETRIGGER = int(CutValue); 
      if (dump) cout << "IGNORETRIGGER      " << IGNORETRIGGER << endl;
      ibin = 5;
      hcuts->SetBinContent(ibin, IGNORETRIGGER);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore trigger :: %i", CutName, IGNORETRIGGER));
    }

    if (!strcmp(CutName, "NOPRESELECTION")) {
      NOPRESELECTION = int(CutValue); 
      if (dump) cout << "NOPRESELECTION     " << NOPRESELECTION << endl;
      ibin = 6;
      hcuts->SetBinContent(ibin, NOPRESELECTION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore preselection :: %i", CutName, NOPRESELECTION));
    }

    if (!strcmp(CutName, "CANDPTLO")) {
      CANDPTLO = CutValue; 
      if (dump) cout << "CANDPTLO:           " << CANDPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, CANDPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(%s) :: %3.1f", CutName, cstring.c_str(), CANDPTLO));
    }

    if (!strcmp(CutName, "CANDETALO")) {
      CANDETALO = CutValue; 
      if (dump) cout << "CANDETALO:           " << CANDETALO << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, CANDETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(%s) :: %3.1f", CutName, cstring.c_str(), CANDETALO));
    }

    if (!strcmp(CutName, "CANDETAHI")) {
      CANDETAHI = CutValue; 
      if (dump) cout << "CANDETAHI:           " << CANDETAHI << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, CANDETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(%s) :: %3.1f", CutName, cstring.c_str(), CANDETAHI));
    }

    if (!strcmp(CutName, "CANDCOSALPHA")) {
      CANDCOSALPHA = CutValue; 
      if (dump) cout << "CANDCOSALPHA:           " << CANDCOSALPHA << endl;
      ibin = 20;
      hcuts->SetBinContent(ibin, CANDCOSALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: cos#alpha :: %5.4f", CutName, CANDCOSALPHA));
    }

    if (!strcmp(CutName, "CANDALPHA")) {
      CANDALPHA = CutValue; 
      if (dump) cout << "CANDALPHA:           " << CANDALPHA << endl;
      ibin = 21;
      hcuts->SetBinContent(ibin, CANDALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #alpha :: %5.4f", CutName, CANDALPHA));
    }

    if (!strcmp(CutName, "CANDFLS3D")) {
      CANDFLS3D = CutValue; 
      if (dump) cout << "CANDFLS3D:           " << CANDFLS3D << endl;
      ibin = 22;
      hcuts->SetBinContent(ibin, CANDFLS3D);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{3d}/#sigma(l_{3d}) :: %3.1f", CutName, CANDFLS3D));
    }

    if (!strcmp(CutName, "CANDFLSXY")) {
      CANDFLSXY = CutValue; 
      if (dump) cout << "CANDFLSXY:           " << CANDFLSXY << endl;
      ibin = 23;
      hcuts->SetBinContent(ibin, CANDFLSXY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{xy}/#sigma(l_{xy}) :: %3.1f", CutName, CANDFLSXY));
    }

    if (!strcmp(CutName, "CANDVTXCHI2")) {
      CANDVTXCHI2 = CutValue; 
      if (dump) cout << "CANDVTXCHI2:           " << CANDVTXCHI2 << endl;
      ibin = 24;
      hcuts->SetBinContent(ibin, CANDVTXCHI2);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #chi^{2} :: %3.1f", CutName, CANDVTXCHI2));
    }

    if (!strcmp(CutName, "CANDISOLATION")) {
      CANDISOLATION = CutValue; 
      if (dump) cout << "CANDISOLATION:           " << CANDISOLATION << endl;
      ibin = 25;
      hcuts->SetBinContent(ibin, CANDISOLATION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: I_{trk} :: %4.2f", CutName, CANDISOLATION));
    }

    if (!strcmp(CutName, "CANDDOCATRK")) {
      CANDDOCATRK = CutValue; 
      if (dump) cout << "CANDDOCATRK:           " << CANDDOCATRK << endl;
      ibin = 26;
      hcuts->SetBinContent(ibin, CANDDOCATRK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{trk} :: %4.3f", CutName, CANDDOCATRK));
    }

    if (!strcmp(CutName, "CANDCLOSETRK")) {
      CANDCLOSETRK = CutValue; 
      if (dump) cout << "CANDCLOSETRK:           " << CANDCLOSETRK << endl;
      ibin = 27;
      hcuts->SetBinContent(ibin, CANDCLOSETRK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: N_{close tracks} :: %4.2f", CutName, CANDCLOSETRK));
    }

    if (!strcmp(CutName, "PVAVEW8")) {
      PVAVEW8 = CutValue; 
      if (dump) cout << "PVAVEW8:           " << PVAVEW8 << endl;
      ibin = 28;
      hcuts->SetBinContent(ibin, PVAVEW8);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: <w^{PV}_{trk}> :: %4.3f", CutName, PVAVEW8));
    }

    if (!strcmp(CutName, "CANDLIP")) {
      CANDLIP = CutValue; 
      if (dump) cout << "CANDLIP:           " << CANDLIP << endl;
      ibin = 40;
      hcuts->SetBinContent(ibin, CANDLIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{z} :: %4.3f", CutName, CANDLIP));
    }

    if (!strcmp(CutName, "CANDLIPS")) {
      CANDLIPS = CutValue; 
      if (dump) cout << "CANDLIPS:          " << CANDLIPS << endl;
      ibin = 41;
      hcuts->SetBinContent(ibin, CANDLIPS);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{z}/#sigma(l_{z}) :: %4.3f", CutName, CANDLIPS));
    }

    if (!strcmp(CutName, "CAND2LIP")) {
      CANDLIP2 = CutValue; 
      if (dump) cout << "CAND2LIP:           " << CANDLIP2 << endl;
      ibin = 42;
      hcuts->SetBinContent(ibin, CANDLIP2);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{z,2} :: %4.3f", CutName, CANDLIP2));
    }

    if (!strcmp(CutName, "CAND2LIPS")) {
      CANDLIPS2 = CutValue; 
      if (dump) cout << "CAND2LIPS:          " << CANDLIPS2 << endl;
      ibin = 43;
      hcuts->SetBinContent(ibin, CANDLIPS2);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{z,2}/#sigma(l_{z,2}) :: %4.3f", CutName, CANDLIPS2));
    }

    if (!strcmp(CutName, "CANDIP")) {
      CANDIP = CutValue; 
      if (dump) cout << "CANDIP:          " << CANDIP << endl;
      ibin = 44;
      hcuts->SetBinContent(ibin, CANDIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{3d} :: %4.3f", CutName, CANDIP));
    }

    if (!strcmp(CutName, "CANDIPS")) {
      CANDIPS = CutValue; 
      if (dump) cout << "CANDIPS:          " << CANDIPS << endl;
      ibin = 45;
      hcuts->SetBinContent(ibin, CANDIPS);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: l_{3d}/#sigma(l_{3d}) :: %4.3f", CutName, CANDIPS));
    }

    if (!strcmp(CutName, "MAXDOCA")) {
      CANDDOCA = CutValue; 
      if (dump) cout << "MAXDOCA:         " << CANDDOCA << endl;
      ibin = 46;
      hcuts->SetBinContent(ibin, CANDDOCA);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: d :: %4.3f", CutName, CANDDOCA));
    }


    if (!strcmp(CutName, "SIGBOXMIN")) {
      SIGBOXMIN = CutValue; 
      if (dump) cout << "SIGBOXMIN:           " << SIGBOXMIN << endl;
      ibin = 90;
      hcuts->SetBinContent(ibin, SIGBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMIN :: %6.3f", CutName, SIGBOXMIN));
    }

    if (!strcmp(CutName, "SIGBOXMAX")) {
      SIGBOXMAX = CutValue; 
      if (dump) cout << "SIGBOXMAX:           " << SIGBOXMAX << endl;
      ibin = 91;
      hcuts->SetBinContent(ibin, SIGBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMAX :: %6.3f", CutName, SIGBOXMAX));
    }

    if (!strcmp(CutName, "BGLBOXMIN")) {
      BGLBOXMIN = CutValue; 
      if (dump) cout << "BGLBOXMIN:           " << BGLBOXMIN << endl;
      ibin = 92;
      hcuts->SetBinContent(ibin, BGLBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMIN :: %6.3f", CutName, BGLBOXMIN));
    }

    if (!strcmp(CutName, "BGLBOXMAX")) {
      BGLBOXMAX = CutValue; 
      if (dump) cout << "BGLBOXMAX:           " << BGLBOXMAX << endl;
      ibin = 93;
      hcuts->SetBinContent(ibin, BGLBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMAX :: %6.3f", CutName, BGLBOXMAX));
    }

    if (!strcmp(CutName, "BGHBOXMIN")) {
      BGHBOXMIN = CutValue; 
      if (dump) cout << "BGHBOXMIN:           " << BGHBOXMIN << endl;
      ibin = 94;
      hcuts->SetBinContent(ibin, BGHBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMIN :: %6.3f", CutName, BGHBOXMIN));
    }

    if (!strcmp(CutName, "BGHBOXMAX")) {
      BGHBOXMAX = CutValue; 
      if (dump) cout << "BGHBOXMAX:           " << BGHBOXMAX << endl;
      ibin = 95;
      hcuts->SetBinContent(ibin, BGHBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMAX :: %6.3f", CutName, BGHBOXMAX));
    }

    // -- Tracks
    if (!strcmp(CutName, "TRACKQUALITY")) {
      TRACKQUALITY = static_cast<int>(CutValue); 
      if (dump) cout << "TRACKQUALITY:           " << TRACKQUALITY << " " << endl;
      ibin = 100;
      hcuts->SetBinContent(ibin, TRACKQUALITY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: track quality :: %d", CutName, TRACKQUALITY));
    }

    if (!strcmp(CutName, "TRACKPTLO")) {
      TRACKPTLO = CutValue; 
      if (dump) cout << "TRACKPTLO:           " << TRACKPTLO << " GeV" << endl;
      ibin = 101;
      hcuts->SetBinContent(ibin, TRACKPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(track) :: %3.1f", CutName, TRACKPTLO));
    }

    if (!strcmp(CutName, "TRACKPTHI")) {
      TRACKPTHI = CutValue; 
      if (dump) cout << "TRACKPTHI:           " << TRACKPTHI << " GeV" << endl;
      ibin = 102;
      hcuts->SetBinContent(ibin, TRACKPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(track) :: %3.1f", CutName, TRACKPTHI));
    }

    if (!strcmp(CutName, "TRACKTIP")) {
      TRACKTIP = CutValue; 
      if (dump) cout << "TRACKTIP:           " << TRACKTIP << " cm" << endl;
      ibin = 103;
      hcuts->SetBinContent(ibin, TRACKTIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{xy}(track) :: %3.1f", CutName, TRACKTIP));
    }

    if (!strcmp(CutName, "TRACKLIP")) {
      TRACKLIP = CutValue; 
      if (dump) cout << "TRACKLIP:           " << TRACKLIP << " cm" << endl;
      ibin = 104;
      hcuts->SetBinContent(ibin, TRACKLIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{z}(track) :: %3.1f", CutName, TRACKLIP));
    }

    if (!strcmp(CutName, "TRACKETALO")) {
      TRACKETALO = CutValue; 
      if (dump) cout << "TRACKETALO:           " << TRACKETALO << " " << endl;
      ibin = 105;
      hcuts->SetBinContent(ibin, TRACKETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{min}(track) :: %3.1f", CutName, TRACKETALO));
    }

    if (!strcmp(CutName, "TRACKETAHI")) {
      TRACKETAHI = CutValue; 
      if (dump) cout << "TRACKETAHI:           " << TRACKETAHI << " " << endl;
      ibin = 106;
      hcuts->SetBinContent(ibin, TRACKETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{max}(track) :: %3.1f", CutName, TRACKETAHI));
    }

    // -- Muons
    if (!strcmp(CutName, "MUIDMASK")) {
      MUIDMASK = int(CutValue); 
      if (dump) cout << "MUIDMASK:           " << MUIDMASK << endl;
      ibin = 200;
      hcuts->SetBinContent(ibin, MUIDMASK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDMask :: %d", CutName, MUIDMASK));
    }

    if (!strcmp(CutName, "MUIDRESULT")) {
      // MUIDRESULT == 0: compare result of & with ">=0"
      // MUIDRESULT != 0: compare result of & with "==MUIDRESULT"
      MUIDRESULT = int(CutValue); 
      if (dump) cout << "MUIDRESULT:           " << MUIDRESULT << endl;
      ibin = 201;
      hcuts->SetBinContent(ibin, MUIDRESULT);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDResult :: %d", CutName, MUIDRESULT));
    }

    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; 
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 202;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(#mu) :: %3.1f", CutName, MUPTLO));
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; 
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 203;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(#mu) :: %3.1f", CutName, MUPTHI));
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; 
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 204;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(#mu) :: %3.1f", CutName, MUETALO));
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; 
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 205;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(#mu) :: %3.1f", CutName, MUETAHI));
    }

    if (!strcmp(CutName, "MUIP")) {
      MUIP = CutValue; 
      if (dump) cout << "MUIP:           " << MUIP << endl;
      ibin = 206;
      hcuts->SetBinContent(ibin, MUIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: IP(#mu) :: %3.1f", CutName, MUIP));
    }

    if (!strcmp(CutName, "MUBDTXML")) {
      char xml[1000]; 
      sscanf(buffer, "%s %s", CutName, xml);
      string tl(xml); 
      if (dump) {
	cout << "MUBDTXML:       " << xml << endl; 
      }
      ibin = 207; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, xml));
      fMvaMuonID = setupMuonMvaReader(string(xml), mrd); 
    }
    

    if (!strcmp(CutName, "MUBDT")) {
      MUBDT = CutValue; 
      if (dump) cout << "MUBDT:           " << MUBDT << endl;
      ibin = 208;
      hcuts->SetBinContent(ibin, MUBDT);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BDT(#mu) :: %3.1f", CutName, MUBDT));
    }


    sscanf(buffer, "%s %s", CutName, XmlName);
    string ctmp = CutName; 
    string sXmlName;
    replaceAll(ctmp, " ", ""); 
    // -- barrel
    if (!strcmp(ctmp.c_str(), "xml0")) {
      sXmlName = "weights/" + string(XmlName) + "-Events0_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents0.push_back(setupReader(sXmlName, frd)); 
      sXmlName = "weights/" + string(XmlName) + "-Events1_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents1.push_back(setupReader(sXmlName, frd)); 
      sXmlName = "weights/" + string(XmlName) + "-Events2_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents2.push_back(setupReader(sXmlName, frd)); 
    }

    // -- endcap
    if (!strcmp(ctmp.c_str(), "xml1")) {
      sXmlName = "weights/" + string(XmlName) + "-Events0_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents0.push_back(setupReader(sXmlName, frd)); 
      sXmlName = "weights/" + string(XmlName) + "-Events1_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents1.push_back(setupReader(sXmlName, frd)); 
      sXmlName = "weights/" + string(XmlName) + "-Events2_BDT.weights.xml"; 
      if (dump) cout << "xml:                   " << sXmlName << endl;
      fReaderEvents2.push_back(setupReader(sXmlName, frd)); 
    }


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
bool candAna::tightMuon(TAnaTrack *pT, bool hadronsPass) {

  const int verbose(0); 

  if (verbose) cout << fYear << " --------- pT = " << pT->fPlab.Perp() << endl;

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    //cout << "NOTRIGGER requested... " << endl;
    return true;
  }

  //             654 3210
  //80 = 0x50 = 0101 0000
  bool muflag = ((pT->fMuID & 80) == 80);
  //  bool muflag = ((pT->fMuID & 16) == 16);
  if (verbose) cout << "muflag: " << hex << pT->fMuID << dec << " -> " << muflag << endl;

  bool mucuts(false); 
  if (verbose) cout << "mu index: " << pT->fMuIndex << " track index: " << pT->fIndex << endl;
  if (pT->fMuIndex > -1) {
    TAnaMuon *pM = fpEvt->getMuon(pT->fMuIndex);
    if (pM->fNmatchedStations > 1) mucuts = true; 
    if (verbose) cout << "matched muon stations: " << pM->fNmatchedStations << " -> " << mucuts << endl;
  }

  bool trackcuts(true); 

  if (TMath::Abs(pT->fBsTip) > 0.2) trackcuts = false;
  if (verbose)  cout << "fBsTip: " << pT->fBsTip << " -> " << trackcuts << endl;
  //if(pT->fBsTip != -99) cout << "fBsTip: " << pT->fBsTip << " -> " << trackcuts << endl;
  if (fpReader->numberOfPixLayers(pT) < 1) trackcuts = false;
  if (verbose)  cout << "pixel layers: " << fpReader->numberOfPixLayers(pT) << " -> " << trackcuts << endl;

  if (fYear == 2011) {
    if (pT->fValidHits < 11) trackcuts = false; 
    if (verbose)  cout << "valid hits: " << pT->fValidHits << " -> " << trackcuts << endl;
  } else if (fYear == 2012) {
    int trkHits = fpReader->numberOfTrackerLayers(pT);
    if (trkHits < 6) trackcuts = false; 
    if (verbose)  cout << "number of tracker layers: " << trkHits << " -> " << trackcuts << endl;
  } else {
    if (pT->fValidHits < 11) trackcuts = false; 
    if (verbose)  cout << "valid hits: " << pT->fValidHits << " -> " << trackcuts << endl;
  }

  if (muflag && mucuts && trackcuts) {
    if (verbose) cout << " +++ passed "<<endl;
    return true; 
  } else {
    //cout<<" failed "<<endl;
    return false;
  }



}


// ----------------------------------------------------------------------
bool candAna::tightMuon(TSimpleTrack *pT, bool hadronsPass) {
  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    return true;
  }
  if (0 == pT->getMuonID()) {
    return false; 
  }

  TAnaMuon *pM(0); 
  int idx = pT->getIndex();
  for (int i = 0; i < fpEvt->nMuons(); ++i) {
    pM = fpEvt->getMuon(i);
    if (idx == pM->fIndex) {
      return tightMuon(pM, hadronsPass); 
    }
  }
  return false; 
}




// ----------------------------------------------------------------------
bool candAna::mvaMuon(TAnaMuon *pt, double &result, bool hadronsPass) {
  
  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    return true;
  }

  if (!tightMuon(pt)) {
    result = -2.;
    return false; 
  }

  mrd.trkValidFract    = pt->fItrkValidFraction; 
  mrd.glbNChi2         = pt->fGtrkNormChi2; 
  mrd.pt               = pt->fPlab.Perp(); 
  mrd.eta              = pt->fPlab.Eta(); 
  mrd.segComp          = pt->fSegmentComp; 
  mrd.chi2LocMom       = pt->fChi2LocalMomentum;
  mrd.chi2LocPos       = pt->fChi2LocalPosition;
  mrd.glbTrackProb     = pt->fGtrkProb;
  mrd.NTrkVHits        = static_cast<float>(pt->fNumberOfValidTrkHits);
  mrd.NTrkEHitsOut     = static_cast<float>(pt->fNumberOfLostTrkHits);

  mrd.dpt                  = pt->fNmatchedStations; 
  mrd.intvalidpixelhits    = fpReader->numberOfPixLayers(pt); // FIXME, kind of correct
  mrd.inttrklayerswithhits = fpReader->numberOfTrackerLayers(pt);
  mrd.intnmatchedstations  = pt->fNmatchedStations;

  mrd.kink             = pt->fMuonChi2;
  
  mrd.dpt              = pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag();
  mrd.dptrel           = TMath::Abs(pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag())/pt->fInnerPlab.Mag();
  mrd.deta             = pt->fInnerPlab.Eta() - pt->fOuterPlab.Eta();
  mrd.dphi             = pt->fInnerPlab.DeltaPhi(pt->fOuterPlab);
  mrd.dr               = pt->fInnerPlab.DeltaR(pt->fOuterPlab);


  result = fMvaMuonID->EvaluateMVA("BDT"); 
  if (result > MUBDT) return true; 
  return false; 
}




// ----------------------------------------------------------------------
bool candAna::mvaMuon(TSimpleTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->getMuonID()) {
    result = -3.; 
    return false; 
  }

  TAnaMuon *pM(0); 
  int idx = pt->getIndex();
  for (int i = 0; i < fpEvt->nMuons(); ++i) {
    pM = fpEvt->getMuon(i);
    if (idx == pM->fIndex) {
      return mvaMuon(pM, result, hadronsPass); 
    }
  }
  return false; 

}


// ----------------------------------------------------------------------
bool candAna::mvaMuon(TAnaTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->fMuID) {
    result = -3.;
    return false; 
  }


  int idx = pt->fMuIndex;
  if (idx > -1 && idx < fpEvt->nMuons()) {
    TAnaMuon *pM = fpEvt->getMuon(idx);
    return mvaMuon(pM, result, hadronsPass); 
  } else {
    cout << "muon index out of range!!!!!" << endl;
  }
  return false; 

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
pair<int, int> candAna::nCloseTracks(TAnaCand *pC, double dcaCut, double dcaCutS, double ptCut) {
  int cnt(0), cnts(0); 
  int nsize = pC->fNstTracks.size(); 
  int pvIdx = pC->fPvIdx;
  int pvIdx2= nearestPV(pvIdx, 0.1);
  if (TMath::Abs(fCandPvLipS2) > 2) pvIdx2 = -1;
  
  TSimpleTrack *pT; 
  double pt(0.); 
  if (nsize > 0) {
    for (int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      
      if (doca > dcaCut) continue; // check the doca cut
      
      pT = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pvIdx)) continue;

      pt = pT->getP().Perp();  
      if (pt < ptCut) continue;
      
      ++cnt;
    }
  }

  if (nsize > 0) {
    for (int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double docas = pC->fNstTracks[i].second.first/pC->fNstTracks[i].second.second;
      
      if (docas > dcaCutS) continue; // check the doca cut
      
      pT = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pvIdx)) continue;
      
      pt = pT->getP().Perp();  
      if (pt < ptCut) continue;
      
      ++cnts;
    }
  }

  return make_pair(cnt, cnts);
}


// ----------------------------------------------------------------------
TAnaCand* candAna::osCand(TAnaCand *pC) {
  TAnaCand *a = new TAnaCand(); 
  a->fPvIdx = pC->fPvIdx; 
  a->fType = pC->fType; 
  a->fSig1 = a->fSig2 = -1; 
  a->fPlab = TVector3(-pC->fPlab.X(), -pC->fPlab.Y(), -pC->fPlab.Z()); //???
  return a; 
}


// ----------------------------------------------------------------------
double candAna::osIsolation(TAnaCand *pC, double r, double ptmin) {
  double iso(0.); 
  TSimpleTrack *pT(0); 
  int overlap(0), verbose(0); 

  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (verbose) {
      cout << "   track " << i 
     	   << " with pT = " << pT->getP().Perp()
     	   << " eta = " << pT->getP().Eta()
     	   << " pointing at PV " << pT->getPvIndex();
    }

    // -- check against overlap to primary candidate
    overlap = 0; 
    for (int j = pC->fSig1; j <= pC->fSig2; ++j) {
      if (i == fpEvt->getSigTrack(j)->fIndex) {
	overlap = 1;
	break;
      }
    }
    if (1 == overlap) continue;
    
    
// ???
//     if ((pT->fPvIdx > -1) && (pT->fPvIdx != pvIdx)) { 	 
//       if (verbose) cout << " doca track " << trkId << " skipped because it is from a different PV " << pT->fPvIdx <<endl; 	 
//       continue; 	 
//     }


    if (pT->getPvIndex() != pC->fPvIdx) { 	 
      if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
      continue;
    }
    
    double pt = pT->getP().Perp(); 
    if (pt < ptmin) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }

    if (pT->getP().DeltaR(pC->fPlab) > r) {
      iso += pt; 
      if (verbose) cout << endl;
    } 
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->getP().DeltaR(pC->fPlab) << endl;
    }
  }
  
  return iso;
 
}

// ----------------------------------------------------------------------
int candAna::osMuon(TAnaCand *pC, double r) {

  double mpt(-1.); 
  int idx (-1); 
  TSimpleTrack *pT(0); 
  int overlap(0), verbose(0); 
  
  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (!tightMuon(pT)) continue;

    // -- check against overlap to primary candidate
    overlap = 0; 
    for (int j = pC->fSig1; j <= pC->fSig2; ++j) {
      if (i == fpEvt->getSigTrack(j)->fIndex) {
	overlap = 1;
	break;
      }
    }
    if (1 == overlap) continue;

    if (verbose) {
      cout << "   track " << i 
     	   << " with pT = " << pT->getP().Perp()
     	   << " eta = " << pT->getP().Eta()
     	   << " pointing at PV " << pT->getPvIndex();
    }
    
    // ???
    //     if (pT->fPvIdx != pvIdx) { 	 
    //       if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
    //       continue;
    //     }
    
    if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pC->fPvIdx)) { 	 
      if (verbose) cout << " track " << i << " skipped because it is from a different PV " << pT->getPvIndex() <<endl; 	 
      continue; 	 
    }
    
    if (pT->getP().DeltaR(pC->fPlab) > r) {
      if (pT->getP().Perp() > mpt) {
	mpt = pT->getP().Perp(); 
	idx = i; 
      }
      if (verbose) cout << endl;
    } 
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->getP().DeltaR(pC->fPlab) << endl;
    }
  }
  
  return idx;


}




// ----------------------------------------------------------------------
double candAna::isoClassicWithDOCA(TAnaCand *pC, double docaCut, double r, double ptmin) {
  const double ptCut(ptmin), coneSize(r); 
  const bool verbose(false);

  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.); 
  TSimpleTrack *ps; 
  vector<int> cIdx, pIdx; 
  int pvIdx = pC->fPvIdx;
  
  fCandI0trk = 0; 
  fCandI1trk = 0; 
  fCandI2trk = 0; 

  getSigTracks(cIdx, pC); 
  for (unsigned int i = 0; i < cIdx.size(); ++i) {
    ps = fpEvt->getSimpleTrack(i); 
 
    if (verbose) cout << " track idx = " << ps->getIndex() << " with ID = " << fpEvt->getSimpleTrackMCID(ps->getIndex()) << endl;
    candPtScalar += ps->getP().Perp(); 
    if (verbose) {
      int tIdx = fpEvt->getSimpleTrack(ps->getIndex())->getPvIndex();
      if (pvIdx != tIdx) {
    	cout << "Signal track pointing to PV " << tIdx << " instead of " << pvIdx << endl;
      }
    }
  }
  
  candPt = pC->fPlab.Perp(); 
  
  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    ps = fpEvt->getSimpleTrack(i); 
    if (verbose) {
      cout << "   track " << i 
     	   << " with pT = " << ps->getP().Perp()
     	   << " eta = " << ps->getP().Eta()
     	   << " pointing at PV " << ps->getPvIndex();
    }
    

    if (ps->getPvIndex() != pvIdx) { 	 
      if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
      continue;
    }
    

    pt = ps->getP().Perp(); 
    if (pt < ptCut) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  {
      if (verbose) cout << " skipped because it is a sig track " << endl;
      continue;
    }
    if (ps->getP().DeltaR(pC->fPlab) < coneSize) {
      pIdx.push_back(i); 
      ++fCandI0trk;
      sumPt += pt; 
      if (verbose) cout << endl;
    } 
    else {
      if (verbose) cout << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
    }
  }

  // -- Now consider the DOCA tracks
  int nsize = pC->fNstTracks.size(); 
  if (nsize>0) {
    for(int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      // double docaE = pC->fNstTracks[i].second.second;

      if(doca > docaCut) continue; // check the doca cut

      ps = fpEvt->getSimpleTrack(trkId);


      if ((ps->getPvIndex() > -1) && (ps->getPvIndex() != pvIdx)) { 	 
	if (verbose) cout << " doca track " << trkId << " skipped because it is from a different PV " << ps->getPvIndex() <<endl; 	 
	continue; 	 
      }

      pt = ps->getP().Perp();  
      if (pt < ptCut) {
	if (verbose) cout << " doca track " << trkId << " skipped because of pt = " << pt << endl;
	continue;
      }

      if (ps->getP().DeltaR(pC->fPlab) > coneSize) {
	if (verbose) cout << " doca track " << trkId << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
	continue;
      }

      // -- Skip tracks already included above
      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue;
      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;

      ++fCandI1trk;
      sumPt += pt; 
      if (verbose) cout << " doca track " << trkId << " included "<<doca<<" "<<pt<<endl;

    } // for loop over tracks
  } // end if 

  fCandI2trk = fCandI0trk + fCandI1trk;

  iso = candPt/(candPt + sumPt); 

  //   if (verbose) cout << "--> iso = " << candPt << " .. " << sumPt << " = " << iso << endl;
  //   if (verbose) cout << "--> iso = " << pC->fPlab.Perp() << " .. " << sumPt << " = " << pC->fPlab.Perp()/(pC->fPlab.Perp() + sumPt) << endl;

  return iso; 
}


// ----------------------------------------------------------------------
void candAna::xpDistMuons() {

  fMu1XpDist = -99.;
  fMu2XpDist = -99.;

  if (fpMuon1->fMuIndex < 0) return;
  if (fpMuon2->fMuIndex < 0) return;

  TAnaMuon *pm[2]; 
  pm[0] = fpEvt->getMuon(fpMuon1->fMuIndex); 
  pm[1] = fpEvt->getMuon(fpMuon2->fMuIndex); 

  //int set0(-99), set1(-99); 
  for (int m = 0; m < 2; ++m) {
    for (int i = 0; i < 10; ++i) {
      
      int trackIdx = pm[m]->fXpTracks[i].idx; 
      if (trackIdx < 0 || trackIdx > fpEvt->nSimpleTracks()) continue;

      // -- ignore this track if it is the other muon of this decay
      if (trackIdx == pm[1-m]->fIndex) continue;

      // -- check that the track's PV is the current one of it is set
      int pvIdx = fpEvt->getSimpleTrack(trackIdx)->getPvIndex();
      if (pvIdx > -1) {
	if (pvIdx != fpCand->fPvIdx) continue;
      }

      // -- now set the mu i xp dist and break for this one
      if (0 == m) {
	fMu1XpDist = pm[m]->fXpTracks[i].dist; 
	//set0 = i; 
	break;
      } 
      if (1 == m) {
	fMu2XpDist = pm[m]->fXpTracks[i].dist; 
	//set1 = i; 
	break;
      } 
    }
  }

  //  cout << "used tracks " << set0 << " and " << set1 << endl;

}



// ----------------------------------------------------------------------
double candAna::isoMuon(TAnaCand *pCand, TAnaMuon *pMuon) {

  TSimpleTrack *sTrack;
  std::vector<near_track_t> nearTracks;
  std::map<int,float>::const_iterator it;
  std::map<int,int> cand_tracks;
  near_track_t nt;
  size_t k;
	
  // get the candidate structure
  findAllTrackIndices(pCand, &cand_tracks);
	
  TVector3 plabMu = pMuon->fPlab;
  for (it = pMuon->fNstTracks.begin(); it != pMuon->fNstTracks.end(); ++it) {
    
    // no tracks from candidate...
    if (cand_tracks.count(it->first) > 0)
      continue;
    
    sTrack = fpEvt->getSimpleTrack(it->first);
    
    // no tracks from foreign primary vertex
    if (sTrack->getPvIndex() >= 0 && sTrack->getPvIndex() != pCand->fPvIdx)
      continue;
    
    nt.ix = it->first;
    nt.doca = it->second;
    nt.p = sTrack->getP().Mag();
    nt.pt = sTrack->getP().Perp();
    TVector3 recTrack = sTrack->getP();
    nt.pt_rel = (recTrack - (recTrack * plabMu) * plabMu).Mag() / plabMu.Mag2();
    nt.deltaR = plabMu.DeltaR(recTrack);
    nearTracks.push_back(nt);
  }
	
  // compute isolation variable
  double result = 0;
  for (k = 0; k < nearTracks.size(); k++) {
    if (nearTracks[k].deltaR < 0.5 && nearTracks[k].pt > 0.5 && nearTracks[k].doca < 0.1) // 1 mm
      result += nearTracks[k].p;
  }
  result = plabMu.Mag()/(plabMu.Mag() + result);
  
  return result;
}


// ----------------------------------------------------------------------
void candAna::findAllTrackIndices(TAnaCand* pCand, map<int,int> *indices) {
  int j;
	
  // iterate through all own tracks. has to be done first, so the duplicate signal tracks
  // won't be added in the daughter anymore
  for (j = pCand->fSig1; j <= pCand->fSig2 && j>=0; j++)
    indices->insert(make_pair(fpEvt->getSigTrack(j)->fIndex,j));
	
  for (j = pCand->fDau1; j <= pCand->fDau2 && j>=0; j++)
    findAllTrackIndices(fpEvt->getCand(j),indices);
}


// ----------------------------------------------------------------------
void candAna::muScaleCorrectedMasses() {
  fCandM3 = fCandM4 = -99.;
  TLorentzVector myNegMuon, myPosMuon; 
  if (fpMuon1->fQ < 0) { 
    myNegMuon.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
    myPosMuon.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  } else {
    myPosMuon.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
    myNegMuon.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  }

  double mass0 = (myNegMuon + myPosMuon).Mag(); 
  fpReader->msc->applyPtCorrection(myNegMuon,-1);
  fpReader->msc->applyPtCorrection(myPosMuon,1);
  double mass1 = (myNegMuon + myPosMuon).Mag(); 
  //  cout << "candMass: " << fCandM << " TLV mass = " << mass0 << " msc mass = " << mass1 << endl;
  fCandM3 = mass1;
  fCandM4 = (mass1/mass0)*fCandM; 

}

 
// ----------------------------------------------------------------------
double candAna::constrainedMass() {
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
    return pMcCand->fMass;
  } else {
    return -2.;
  }
}


// ----------------------------------------------------------------------
void candAna::runRange() {
  fRunRange = 6;
  if ((fRun >= 160329) && (fRun <= 163261)) {
    fRunRange = 0; // the beginning
  } 
  if ((fRun >= 163269) && (fRun <= 163869)) {
    fRunRange = 1; // displaced J/psi
  } 
  if ((fRun >= 165088) && (fRun <= 167913)) {
    fRunRange = 2; // HLTDimuon7
  } 
  if ((fRun >= 170249) && (fRun <= 173198)) {
    fRunRange = 3; // 2e33
  } 
  if ((fRun >= 173236) && (fRun <= 178380)) {
    fRunRange = 4; // 3e33 WITH the ETA cut!
  } 
  if ((fRun >= 178420) && (fRun <= 999999)) {
    fRunRange = 5; // 5e33
  } 

  if (fIsMC) {
    if (string::npos != fHLTPath.find("HLT_Dimuon7_Jpsi_Displaced_v1")) {
      fRunRange = 2;
    }
    if (string::npos != fHLTPath.find("HLT_DoubleMu3p5_Jpsi_Displaced_v2_Bs")) {
      fRunRange = 3;
    }
    if (string::npos != fHLTPath.find("HLT_DoubleMu4_Jpsi_Displaced_v1")) {
      fRunRange = 4;
    }
  }

}




// ----------------------------------------------------------------------
int candAna::nearestPV(int pvIdx, double maxDist) {

  if (pvIdx==-1) return -1;  // add protection d.k. 9/6/12

  TAnaVertex *v0 = fpEvt->getPV(pvIdx); 

  if (0 == v0) return -1;  // add protection

  double zV0 = v0->fPoint.Z(); 
  
  int idx(-1); 
  double z(0.), zabs(0.), zmin(99.);
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    if (i == pvIdx) continue;
    z = fpEvt->getPV(i)->fPoint.Z();
    zabs = TMath::Abs(zV0 - z); 
    if (zabs < zmin) {
      idx = i; 
      zmin = zabs;
    }
  }
  
  if (zmin < maxDist) {
    //     cout << "pvIdx = " << pvIdx << " at = " << zV0 
    // 	 << ", nearest other PV with idx = " << idx << " at z = " << fpEvt->getPV(idx)->fPoint.Z() 
    // 	 << " and delta(z) = " << zmin << endl;
    return idx; 
  } else {
    return -1;
  }
}


// ----------------------------------------------------------------------
void candAna::getSigTracks(vector<int> &v, TAnaCand *pC) {
  TAnaCand *pD; 
  TAnaTrack *pT; 
  vector<int> bla; 

  // -- loop over daughters
  if (pC->fDau1 > -1) {
    for (int j = pC->fDau1; j <= pC->fDau2; ++j) {
      pD = fpEvt->getCand(j); 
      getSigTracks(bla, pD); 
    }

    for (unsigned j = 0; j < bla.size(); ++j) v.push_back(bla[j]);
  }

  // -- add direct sigtracks
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    if (v.end() == find(v.begin(), v.end(), pT->fIndex)) {
      v.push_back(pT->fIndex); 
    }
  }

}


// ----------------------------------------------------------------------
void candAna::calcBDT() {
  fBDT = -99.;
  //??  if (5 == mode && 5.2 < mass && mass < 5.45 && fb.iso < 0.7) continue; 
  if (fChan < 0) return;
  if (0 == fReaderEvents0.size()) {
    cout << "no BDT defined" << endl;
    return;
  }

  if (!preselection(fRTD, fChan)) return;

  //   if (fCandPt > 100) return;
  //   if (fCandPt < 6) return;
  //   if (fMu1Pt < 4) return;
  //   if (fMu2Pt < 4) return;
  //   if (fCandFL3d > 1.5) return;
  //   if (fCandFL3d < 0.) return;
  //   if (fCandM > 5.9) return;
  //   if (fCandM < 4.9) return;
  
  //   if (!fb.hlt) return;
  //   if (!fb.gmuid) return;
  
  frd.pt = fCandPt; 
  frd.eta = fCandEta; 
  frd.m1eta = fMu1Eta; 
  frd.m2eta = fMu2Eta; 
  frd.m1pt = fMu1Pt; 
  frd.m2pt = fMu2Pt;
  frd.fls3d = fCandFLS3d; 
  frd.alpha = fCandA; 
  frd.maxdoca = fCandDoca;
  frd.pvip = fCandPvIp; 
  frd.pvips = fCandPvIpS; 
  frd.iso = fCandIso; 
  frd.docatrk = fCandDocaTrk; 
  frd.chi2dof = fCandChi2/fCandDof; 
  frd.closetrk = fCandCloseTrk; 
  
  frd.m  = fCandM; 
  //  cout << "Evt = " << fEvt << " %3 = " << fEvt%3 << " chan = " << fChan << " " << " etas = " << fMu1Eta << " " << fMu2Eta;
  if (0 == fEvt%3) {
    fBDT   = fReaderEvents0[fChan]->EvaluateMVA("BDT"); 
  } else if (1 == fEvt%3) {
    fBDT   = fReaderEvents1[fChan]->EvaluateMVA("BDT"); 
  } else if (2 == fEvt%3) {
    fBDT   = fReaderEvents2[fChan]->EvaluateMVA("BDT"); 
  } else {
    cout << "all hell break loose" << endl;
  }
  //  cout << " bdt = " << fBDT << endl;
}


// ----------------------------------------------------------------------
int candAna::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < 1.4 && TMath::Abs(m2eta) < 1.4) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
}


// ----------------------------------------------------------------------
TMVA::Reader* candAna::setupMuonMvaReader(string xmlFile, mvaMuonIDData &d) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

  TString dir    = "weights/";
  TString methodNameprefix = "BDT";

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  string weightFile = "weights/TMVA-" + xmlFile + ".weights.xml"; 
  cout << "setupMuonMvaReader, open file " << weightFile << endl;
  ifstream is(weightFile.c_str()); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1); 
  string::size_type m1, m2;
  string stype; 
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	if (stype == "trkValidFract") {
	  cout << "  adding trkValidFract" << endl;
	  reader->AddVariable( "trkValidFract", &d.trkValidFract);
	  continue;
	}
	if (stype == "glbNChi2") {
	  reader->AddVariable("glbNChi2", &d.glbNChi2);
	  cout << "  adding glbNChi2" << endl;
	  continue;
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable("eta", &d.eta);
	  continue;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable("pt", &d.pt);
	  continue;
	}
	if (stype == "segComp") {
	  cout << "  adding segComp" << endl;
	  reader->AddVariable("segComp", &d.segComp);
	  continue;
	}
	if (stype == "chi2LocMom") {
	  cout << "  adding chi2LocMom" << endl;
	  reader->AddVariable("chi2LocMom", &d.chi2LocMom);
	  continue;
	}
	if (stype == "chi2LocPos") {
	  cout << "  adding chi2LocPos" << endl;
	  reader->AddVariable("chi2LocPos", &d.chi2LocPos);
	  continue;
	}
	if (stype == "glbTrackProb") {
	  cout << "  adding glbTrackProb" << endl;
	  reader->AddVariable("glbTrackProb", &d.glbTrackProb);
	  continue;
	}
	if (stype == "NTrkVHits") {
	  cout << "  adding NTrkVHits" << endl;
	  reader->AddVariable("NTrkVHits", &d.NTrkVHits);
	  continue;
	}
	if (stype == "NTrkEHitsOut") {
	  cout << "  adding NTrkEHitsOut" << endl;
	  reader->AddVariable("NTrkEHitsOut", &d.NTrkEHitsOut);
	  continue;
	}
	if (stype == "NTrkEHitsOut") {
	  cout << "  adding NTrkEHitsOut" << endl;
	  reader->AddVariable("NTrkEHitsOut", &d.NTrkEHitsOut);
	  continue;
	}

	// -- new convention for UL's BDT
	if (stype == "intnmatchedstations") {
	  cout << "  adding intnmatchedstations" << endl;
	  reader->AddVariable("intnmatchedstations", &d.intnmatchedstations);
	  continue;
	}
	if (stype == "intvalidpixelhits") {
	  cout << "  adding intvalidpixelhits" << endl;
	  reader->AddVariable("intvalidpixelhits", &d.intvalidpixelhits);
	  continue;
	}
	if (stype == "inttrklayerswithhits") {
	  cout << "  adding inttrklayerswithhits" << endl;
	  reader->AddVariable("inttrklayerswithhits", &d.inttrklayerswithhits);
	  continue;
	}
	if (stype == "gchi2") {
	  cout << "  adding gchi2" << endl;
	  reader->AddVariable("gchi2", &d.glbNChi2);
	  continue;
	}
	if (stype == "itrkvalidfraction") {
	  cout << "  adding itrkvalidfraction" << endl;
	  reader->AddVariable("itrkvalidfraction", &d.trkValidFract);
	  continue;
	}
	if (stype == "segcomp") {
	  cout << "  adding segcomp" << endl;
	  reader->AddVariable("segcomp", &d.segComp);
	  continue;
	}
	if (stype == "chi2lmom") {
	  cout << "  adding chi2lmom" << endl;
	  reader->AddVariable("chi2lmom", &d.chi2LocMom);
	  continue;
	}
	if (stype == "chi2lpos") {
	  cout << "  adding chi2lpos" << endl;
	  reader->AddVariable("chi2lpos", &d.chi2LocPos);
	  continue;
	}
	if (stype == "gtrkprob") {
	  cout << "  adding gtrkprob" << endl;
	  reader->AddVariable("gtrkprob", &d.glbTrackProb);
	  continue;
	}
	if (stype == "ntrkvhits") {
	  cout << "  adding ntrkvhits" << endl;
	  reader->AddVariable("ntrkvhits", &d.NTrkVHits);
	  continue;
	}
	if (stype == "inttrklayerswithhits") {
	  cout << "  adding inttrklayerswithhits" << endl;
	  reader->AddVariable("inttrklayerswithhits", &d.inttrklayerswithhits);
	  continue;
	}
	if (stype == "ntrkehitsout") {
	  cout << "  adding ntrkehitsout" << endl;
	  reader->AddVariable("ntrkehitsout", &d.NTrkEHitsOut);
	  continue;
	}
	if (stype == "kink") {
	  cout << "  adding kink" << endl;
	  reader->AddVariable("kink", &d.kink);
	  continue;
	}
	if (stype == "dpt") {
	  cout << "  adding dpt" << endl;
	  reader->AddVariable("dpt", &d.dpt);
	  continue;
	}
	if (stype == "deta") {
	  cout << "  adding deta" << endl;
	  reader->AddVariable("deta", &d.deta);
	  continue;
	}
	if (stype == "dphi") {
	  cout << "  adding dphi" << endl;
	  reader->AddVariable("dphi", &d.dphi);
	  continue;
	}
	if (stype == "dr") {
	  cout << "  adding dr" << endl;
	  reader->AddVariable("dr", &d.dr);
	  continue;
	}
	if (stype == "dptrel") {
	  cout << "  adding dptrel" << endl;
	  reader->AddVariable("dptrel", &d.dptrel);
	  continue;
	}


      }
      break;
    }
  }

  reader->BookMVA("BDT", TString(weightFile.c_str())); 
  return reader; 
}


// ----------------------------------------------------------------------
TMVA::Reader* candAna::setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  
  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  //  TString methodName = TString(fBdt) + TString(" method");
  //  TString weightfile = dir + fBdt + "_" + methodNameprefix + TString(".weights.xml");
  TString weightfile = xmlFile;

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  cout << "setupReader, open file " << weightfile << endl;
  ifstream is(weightfile); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1); 
  string::size_type m1, m2;
  string stype; 
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	//	cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
	if (stype == "m1pt") {
	  cout << "  adding m1pt" << endl;
	  reader->AddVariable( "m1pt", &rd.m1pt);
	}
	if (stype == "m2pt") {
	  cout << "  adding m2pt" << endl;
	  reader->AddVariable( "m2pt", &rd.m2pt);
	}
	if (stype == "m1eta") {
	  cout << "  adding m1eta" << endl;
	  reader->AddVariable( "m1eta", &rd.m1eta);
	}
	if (stype == "m2eta") {
	  reader->AddVariable( "m2eta", &rd.m2eta);
	  cout << "  adding m2eta" << endl;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable( "pt", &rd.pt);
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable( "eta", &rd.eta);
	}
	if (stype == "fls3d") {
	  cout << "  adding fls3d" << endl;
	  reader->AddVariable( "fls3d", &rd.fls3d);
	}
	if (stype == "alpha") {
	  cout << "  adding alpha" << endl;
	  reader->AddVariable( "alpha", &rd.alpha);
	}
	if (stype == "maxdoca") {
	  cout << "  adding maxdoca" << endl;
	  reader->AddVariable( "maxdoca", &rd.maxdoca);
	}
	if (stype == "pvip") {
	  cout << "  adding pvip" << endl;
	  reader->AddVariable( "pvip", &rd.pvip);
	}
	if (stype == "pvips") {
	  cout << "  adding pvips" << endl;
	  reader->AddVariable( "pvips", &rd.pvips);
	}
	if (stype == "iso") {
	  cout << "  adding iso" << endl;
	  reader->AddVariable( "iso", &rd.iso);
	}
	if (stype == "docatrk") {
	  cout << "  adding docatrk" << endl;
	  reader->AddVariable( "docatrk", &rd.docatrk);
	}
	if (stype == "closetrk") {
	  cout << "  adding closetrk" << endl;
	  reader->AddVariable( "closetrk", &rd.closetrk);
	}
	if (stype == "closetrks1") {
	  cout << "  adding closetrks1" << endl;
	  reader->AddVariable( "closetrks1", &rd.closetrks1);
	}
	if (stype == "closetrks2") {
	  cout << "  adding closetrks2" << endl;
	  reader->AddVariable( "closetrks2", &rd.closetrks2);
	}
	if (stype == "closetrks3") {
	  cout << "  adding closetrks3" << endl;
	  reader->AddVariable( "closetrks3", &rd.closetrks3);
	}
	if (stype == "chi2dof") {
	  cout << "  adding chi2dof" << endl;
	  reader->AddVariable( "chi2dof", &rd.chi2dof);
	}
	if (stype == "m1iso") {
	  cout << "  adding m1iso" << endl;
	  reader->AddVariable( "m1iso", &rd.m1iso);
	}
	if (stype == "m2iso") {
	  cout << "  adding m2iso" << endl;
	  reader->AddVariable( "m2iso", &rd.m2iso);
	}
	if (stype == "pvdchi2") {
	  cout << "  adding pvdchi2" << endl;
	  reader->AddVariable( "pvdchi2", &rd.pvdchi2);
	}
	if (stype == "othervtx") {
	  cout << "  adding othervtx" << endl;
	  reader->AddVariable( "othervtx", &rd.othervtx);
	}
	if (stype == "pvlip2") {
	  cout << "  adding pvlip2" << endl;
	  reader->AddVariable( "pvlip2", &rd.pvlip2);
	}
	if (stype == "pvlips2") {
	  cout << "  adding pvlips2" << endl;
	  reader->AddVariable( "pvlips2", &rd.pvlips2);
	}
      }
      break;
    }
  }
  
  nvars = -1; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
	if (stype == "m") {
	  cout << "  adding m as spectator" << endl;
	  reader->AddSpectator( "m", &rd.m);  
	}
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", weightfile); 
  return reader; 
}


// ----------------------------------------------------------------------
void candAna::replaceAll(std::string &s, std::string a, std::string b) {
  
  TString ts(s.c_str()); 
  ts.ReplaceAll(a.c_str(), b.c_str()); 
  s = ts.Data(); 

}



// ----------------------------------------------------------------------
void candAna::fillRedTreeData() {
  // -- this only fills the variables that are needed for the preselection() function
  fRTD.hlt       = fGoodHLT;
  fRTD.gmuid     = fGoodMuonsID;

  fRTD.pt        = fCandPt;
  fRTD.eta       = fCandEta;
  fRTD.m         = fCandM;

  fRTD.m1pt      = fMu1Pt;
  fRTD.m2pt      = fMu2Pt;

  fRTD.m1eta     = fMu1Eta;
  fRTD.m2eta     = fMu2Eta;

  fRTD.pvip      = fCandPvIp; 
  fRTD.pvips     = fCandPvIpS; 

  fRTD.pvlip     = fCandPvLip; 
  fRTD.pvlips    = fCandPvLipS; 

  fRTD.closetrk  = fCandCloseTrk;
  fRTD.iso       = fCandIso;

  fRTD.flsxy     = fCandFLSxy;
  fRTD.fl3d      = fCandFL3d;
  fRTD.fls3d     = fCandFLS3d;

  fRTD.chi2      = fCandChi2;
  fRTD.dof       = fCandDof;

  fRTD.alpha     = fCandA;
	
  fRTD.docatrk   = fCandDocaTrk;
  fRTD.maxdoca   = fCandDoca;

}

// ----------------------------------------------------------------------
// A simple trigger matcher based on deltaR (from Frank)
// Three versions exist, 
// (void) - checks the 2 muons, use  objects selected by the trigger path 
// (bool)  - check 2 muons, flag for HLT object selection 
// (TAnaTrack, bool) - check a single track matching any muon object, 
//
bool candAna::doTriggerMatching() { // call the normal version with (true)
  bool HLTmatch = doTriggerMatching(true);
  return HLTmatch;
}

//   bool HLTmatch = false;
//   const double deltaRthrsh0(0.2); // initial cone for testing 
//   const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020
//   int mu1match(-1), mu2match(-1);
//   double deltaRminMu1(100), deltaRminMu2(100);
//   string hlt1, hlt2;
//   TTrgObj *tto;
//   TLorentzVector tlvMu1, tlvMu2;
//   const bool localPrint = false;
//   //bool localPrint = false;
//   //localPrint = select_print;

//   if (fVerbose > 15 || localPrint) {
//     cout << "dump trigger objects ----------------------" << fRun << " " << fEvt <<" "<<fCandType<<endl;
//     //cout << "mu1: pt,eta,phi: " << fMu1Pt << " " << fMu1Eta << " " << fMu1Phi << " q: " << fMu1Q << endl;
//     //cout << "mu2: pt,eta,phi: " << fMu2Pt << " " << fMu2Eta << " " << fMu2Phi << " q: " << fMu2Q << endl;
//     cout << "mu1: pt,eta,phi: " << fpMuon1->fPlab.Perp() << " " << fpMuon1->fPlab.Eta() << " " <<fpMuon1->fPlab.Phi() 
// 	 << " q: " << fMu1Q << endl;
//     cout << "mu2: pt,eta,phi: " << fpMuon2->fPlab.Perp() << " " << fpMuon2->fPlab.Eta() << " " <<fpMuon2->fPlab.Phi() 
// 	 << " q: " << fMu2Q << endl;
//   }
  
//   tlvMu1.SetPtEtaPhiM(fpMuon1->fPlab.Perp(),fpMuon1->fPlab.Eta(),fpMuon1->fPlab.Phi(),MMUON);
//   tlvMu2.SetPtEtaPhiM(fpMuon2->fPlab.Perp(),fpMuon2->fPlab.Eta(),fpMuon2->fPlab.Phi(),MMUON);
//   ((TH1D*)fHistDir->Get("test1"))->Fill(1.); 

//   for(int i=0; i!=fpEvt->nTrgObj(); i++) {
//     tto = fpEvt->getTrgObj(i);

//     if (fVerbose > 97 || localPrint) {
//       cout << "i: " << i << " ";
//       ////cout << tto->fLabel << " ";
//       cout << tto->fP.DeltaR(tlvMu1) << ":" << tto->fP.DeltaR(tlvMu2)<<" ";
//       tto->dump();
//     }

//     // WARNING: this works only for 2012 data
//     bool selected = false;
//     if( fIsMC ) { //MC

//       if ( fYear==2012) {

// 	if( (fCandType==3000068 || fCandType==3000067) && tto->fLabel == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::") 
// 	  selected = true; // 2012 data, psik&psiphi
// 	else if ( (fCandType==1000080 || fCandType==1000082|| fCandType==1000091 ) &&  //BsMuMu, BsKK, Bdpipi 
// 		  ( tto->fLabel == "hltVertexmumuFilterBs345:HLT::"   // 2012 data, mumu, central 34
// 		 || tto->fLabel == "hltVertexmumuFilterBs3p545:HLT::" // 2012 data, mumu, central 3p54
// 		 || tto->fLabel == "hltVertexmumuFilterBs47:HLT::") ) // 2012 data, mumu, forward
// 	  selected = true;

//       } else if (fYear == 2011) {


//       } // year 2012


//     } else { // DATA

//       if ( fYear==2012) {

// 	if( (fCandType==300521 || fCandType==300531) && tto->fLabel == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::") // 2012 data, psik&psiphi
// 	  selected = true;
// 	else if ( (fCandType==301313 ||fCandType==1313 || fCandType==211211)&&  // mumu and HH 
// 		  ( tto->fLabel == "hltVertexmumuFilterBs345:HLT::"   // 2012 data, mumu, central 34
// 		 || tto->fLabel == "hltVertexmumuFilterBs3p545:HLT::" // 2012 data, mumu, central 3p54
// 	         || tto->fLabel == "hltVertexmumuFilterBs47:HLT::") ) // 2012 data, mumu, forward
// 	  selected = true;

//       } else if (fYear == 2011) {


//       } // year 2012

//     } // data 

    
//    //|| (fIsMC && (tto->fLabel.BeginsWith("hltDisplacedmumuFilter") || tto->fLabel.BeginsWith("hltVertexmumuFilterDimuon"))) // MC selection needs adjustment according to sample used
      
//     if(selected) {
      
//       double deltaR1 = tto->fP.DeltaR(tlvMu1);
//       double deltaR2 = tto->fP.DeltaR(tlvMu2);
//       if (fVerbose > 16 || localPrint) cout << i << " selected "<< tto->fLabel << deltaR1 << ":" << deltaR2 << endl;
//       if (deltaR1<deltaRthrsh0 && deltaR1<deltaRminMu1) {
// 	deltaRminMu1 = deltaR1;
// 	mu1match = i;
// 	hlt1 = tto->fLabel;
// 	if (fVerbose > 16 || localPrint) cout << " matched 1 " << i << endl;
//       }
//       if (deltaR2<deltaRthrsh0 && deltaR2<deltaRminMu2) {
// 	deltaRminMu2 = deltaR2;
// 	mu2match = i;
// 	hlt2 = tto->fLabel;
// 	if (fVerbose > 16 || localPrint) cout << " matched 2 " << i << endl;
//       }
//     }
      
    
//   } // end for loop 
  
//   HLTmatch = (mu1match>=0 && mu2match >=0 && deltaRminMu1<deltaRthrsh && deltaRminMu2<deltaRthrsh);

//   if (fVerbose > 15 || localPrint) cout << "trigger matched: " << HLTmatch
// 			 << " - result Mu1: " << mu1match << " dR: " << deltaRminMu1
// 			 << " Mu2: " << mu2match << " dR: " << deltaRminMu2 << endl;


//   ((TH1D*)fHistDir->Get("test8"))->Fill(deltaRminMu1); 
//   ((TH1D*)fHistDir->Get("test8"))->Fill(deltaRminMu2); 

//   if(deltaRminMu1<deltaRthrsh) {
//     ((TH1D*)fHistDir->Get("test1"))->Fill(3.); 
//     if (     hlt1 == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::" ) ((TH1D*)fHistDir->Get("test1"))->Fill(6.); 
//     else if (hlt1 == "hltVertexmumuFilterBs345:HLT::")             ((TH1D*)fHistDir->Get("test1"))->Fill(7.);
//     else if (hlt1 == "hltVertexmumuFilterBs3p545:HLT::")           ((TH1D*)fHistDir->Get("test1"))->Fill(8.);
//     else if (hlt1 == "hltVertexmumuFilterBs47:HLT::")              ((TH1D*)fHistDir->Get("test1"))->Fill(9.);
//   }

//   if(deltaRminMu2<deltaRthrsh) {
//     ((TH1D*)fHistDir->Get("test1"))->Fill(4.); 
//     if (     hlt2 == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::" )  ((TH1D*)fHistDir->Get("test1"))->Fill(10.); 
//     else if (hlt2 == "hltVertexmumuFilterBs345:HLT::")              ((TH1D*)fHistDir->Get("test1"))->Fill(11.);
//     else if (hlt2 == "hltVertexmumuFilterBs3p545:HLT::")            ((TH1D*)fHistDir->Get("test1"))->Fill(12.);
//     else if (hlt2 == "hltVertexmumuFilterBs47:HLT::")               ((TH1D*)fHistDir->Get("test1"))->Fill(13.);
//   }

//   if(HLTmatch) ((TH1D*)fHistDir->Get("test1"))->Fill(2.); 

//   //cout<<HLTmatch<<" "<<deltaRminMu1<<" "<<deltaRminMu2<<endl;

//   return HLTmatch;
// }

//------------------------------------
// Checks both muons from the triggered dimuon. NEW
// calls doTriggerMatchingR()
bool candAna::doTriggerMatching(bool anyTrig) {
  bool HLTmatch = false;
  const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020

  double dR1 = doTriggerMatchingR(fpMuon1, anyTrig);
  double dR2 = doTriggerMatchingR(fpMuon2, anyTrig);

  ((TH1D*)fHistDir->Get("test4"))->Fill(dR1); 
  ((TH1D*)fHistDir->Get("test4"))->Fill(dR2); 

  HLTmatch = (dR1<deltaRthrsh && dR2<deltaRthrsh);

  //cout<<HLTmatch<<" "<<dR1<<" "<<dR2<<endl;

  return HLTmatch;
}
// ---------------------------------------------------------------------------------
// To match a single track to a trigger object (selected or all)
// pt - track
// anyTrig - if true use all trigger objects
// calls doTriggerMatchingR()
bool candAna::doTriggerMatching(TAnaTrack *pt, bool anyTrig) {

  const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020
  double dR = doTriggerMatchingR(pt,anyTrig);

  ((TH1D*)fHistDir->Get("test6"))->Fill(dR); 

  bool HLTmatch = (dR<deltaRthrsh );
  return HLTmatch;

}

//-------------------------------------------------------------------------------
// To match a single track to a trigger object (selected or all)
// pt - track
// anyTrig - if true use all trigger objects
double candAna::doTriggerMatchingR(TAnaTrack *pt, bool anyTrig) {
  fTrigMatchDeltaPt = 99.;

  const bool localPrint = false;
  //bool localPrint = true;
  //localPrint = select_print;

  //bool HLTmatch = false;
  //const double deltaRthrsh0(0.2); // initial cone for testing 
  const double deltaRthrsh0(2.0); // initial cone for testing 
  const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020
  int mu1match(-1);
  string hlt1;
  double deltaRminMu1(100);
  TTrgObj *tto;
  TLorentzVector tlvMu1;
  
 
  if (fVerbose > 15 || localPrint) {
    cout << "dump trigger objects ----------------------" << fEvt<< endl;
    cout << "pt,eta,phi: " << pt->fPlab.Perp() << " " << pt->fPlab.Eta() << " " << pt->fPlab.Phi() << endl;
  }
  
  ((TH1D*)fHistDir->Get("test1"))->Fill(20.); 
  tlvMu1.SetPtEtaPhiM(pt->fPlab.Perp(),pt->fPlab.Eta(),pt->fPlab.Phi(),MMUON); // assume a muon
  
  for(int i=0; i!=fpEvt->nTrgObj(); i++) {
    tto = fpEvt->getTrgObj(i);
    
    if (fVerbose > 97  || localPrint ) {
      cout << "i: " << i << " "; 
      //cout << tto->fLabel << tto->fP.DeltaR(tlvMu1)<<" ";
      cout << tto->fP.DeltaR(tlvMu1)<<" ";
      tto->dump();
    }

    // WARNING: this works only for 2012 data    
    // the label changes according to datataking era, so we need to distinguish them
    bool selected = false;

#define OLD_MATCH
    if ( anyTrig ) {
#ifndef OLD_MATCH
      // Use the already selected objects
      if ( tto->fNumber > -1 ) {selected = true;}
#else

      if ( fYear==2012) {
	 //{if( tto->fLabel.Contains("L3") && (tto->fLabel.Contains("mu")  || tto->fLabel.Contains("Mu")) ) selected = true;}
	if( (tto->fLabel.Contains("mu")  || tto->fLabel.Contains("Mu")) ) {
	  selected = true;
	}
      } else { // 2011 cannot do Mu.AND.L3, some triger do not have it in the name
	if( (tto->fLabel.Contains("mu")||tto->fLabel.Contains("Mu")||tto->fLabel.Contains("Jpsi")||
	     tto->fLabel.Contains("Displaced")||tto->fLabel.Contains("Vertex")||tto->fLabel.Contains("LowMass"))
	    ) {
	  selected = true;
	}
      }
#endif      
      // Use the already selected objects 
      if ( tto->fNumber > -1 ) {selected = true;}
      
      // select objects with "mu" or "Mu" in the name.  Add also the "L3" selection. 
      //cout<<tto->fLabel.Contains("mu")<<" "<<tto->fLabel.Contains("Mu")<<endl;
      //if( (tto->fLabel.Contains("mu")  || tto->fLabel.Contains("Mu")) ) selected = true;

//       if ( fYear==2012) 
// 	{if( tto->fLabel.Contains("L3") && (tto->fLabel.Contains("mu")  || tto->fLabel.Contains("Mu")) ) selected = true;}

//       else  // 2011 cannot do Mu.AND.L3, some triger do not have it in the name
// 	{if( (tto->fLabel.Contains("mu")||tto->fLabel.Contains("Mu")||tto->fLabel.Contains("Jpsi")||
// 	      tto->fLabel.Contains("Displaced")||tto->fLabel.Contains("Vertex")||tto->fLabel.Contains("LowMass")) 
// 	     ) selected = true;}


      // select really all 
      //selected = true; 

      // Select all bsmm anaysis triggers
//       if( (tto->fLabel == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::") || //2012 data jpsi disp. 
//        (tto->fLabel == "hltVertexmumuFilterBs345:HLT::") ||   // 2012 data, mumu, central 34
//        (tto->fLabel == "hltVertexmumuFilterBs3p545:HLT::") || // 2012 data, mumu, central 3p54
//        (tto->fLabel == "hltVertexmumuFilterBs47:HLT::")  // 2012 data, mumu, forward
//        ) selected = true;

    } else {  // check triggers

      // data same as MC
      // if( fIsMC ) { //MC
	
      if ( fYear==2012) {
	
	if( (fCandType==300521 || fCandType==300531 || // data
	     fCandType==3000068 || fCandType==3000067) ) {  // MC
	 
	  // Two options, they are similat for trigered events but can differ for nontriggers,
	  // the first one is more restructive, use first 
	  if(tto->fLabel == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::") selected = true;
	  //if(tto->fLabel == "hltDoubleMu4JpsiDisplacedL3Filtered:HLT::") selected = true;
	  
	} else if ( fCandType==301313 ||fCandType==1313 || fCandType==211211 ||  // mumu and HH 
		    fCandType==1000080 || fCandType==1000082|| fCandType==1000091 ) {  //MC BsMuMu, BsKK, Bdpipi 
	  // Same comment as above 
	  if( tto->fLabel == "hltVertexmumuFilterBs345:HLT::"  // 2012 data, mumu, central 34
	      || tto->fLabel == "hltVertexmumuFilterBs3p545:HLT::" // 2012 data, mumu, central 3p54
	      || tto->fLabel == "hltVertexmumuFilterBs47:HLT::")   // 2012 data, mumu, forward
	    //if( tto->fLabel    == "hltDoubleMu34Dimuon5CentralBsL3Filtered:HLT::"   // 2012 data, mumu, central 34
	    //|| tto->fLabel == "hltDoubleMu3p54Dimuon5CentralBsL3Filtered:HLT::" // 2012 data, mumu, central 3p54
	    //|| tto->fLabel == "hltDoubleMu4Dimuon7ForwardBsL3Filtered:HLT::")   // 2012 data, mumu, forward
	    selected = true;
	} // if fCandType
	
      } else if (fYear == 2011) {
	
	if( (fCandType==300521 || fCandType==300531 ||  // psik & psiphi
	     fCandType==3000068 || fCandType==3000067) ) {  // MC
	  if(   tto->fLabel   == "hltDoubleMu3JpsiL3Filtered:HLT::"          //5E32
		||tto->fLabel == "hltDisplacedmumuFilterJpsi:HLT::"           // 5E32, 1E33, 1.4E33,2E33
		||tto->fLabel == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::") // 3E33, 5E33
	    selected = true;
	} else if ( fCandType==301313 || fCandType==1313 || fCandType==211211 ||                        // data mumu and HH 
		    fCandType==1000080 || fCandType==1000082|| fCandType==1000091 ) {  //MC BsMuMu, BsKK, Bdpipi 

	  if ( tto->fLabel   == "hltDoubleMu3BsL3Filtered:HLT::"         // 5E32
	       ||tto->fLabel == "hltDoubleMu2BsL3Filtered:HLT::"         // 1E33, 1.4E33
	       ||tto->fLabel == "hltDoubleMu2BarrelBsL3Filtered:HLT::"   // 2E33 central
	       ||tto->fLabel == "hltDoubleMu2Dimuon6BsL3Filtered:HLT::"  // 2E33 forward
	       ||tto->fLabel == "hltVertexmumuFilterBs4:HLT::"           // 3E33, 5E33, central
	       ||tto->fLabel == "hltVertexmumuFilterBs6:HLT::" )         // 3E33, 5E33, forward
	    selected = true;
	} // if cand type 
	
      } // year 
      
    } // allTrig
    
    if(selected) {
      double deltaR1 = tto->fP.DeltaR(tlvMu1);
      ((TH1D*)fHistDir->Get("test8"))->Fill(deltaR1); 
      if (fVerbose > 16 || localPrint) cout << i<<" "<<tto->fLabel << " "<<deltaR1 << endl;
      
      if (deltaR1<deltaRthrsh0 && deltaR1<deltaRminMu1) {
	if (fVerbose > 16 || localPrint) {
	  cout << " selected "<< deltaR1 <<" ";
	  tto->dump();
          //cout<<endl;
        }
        deltaRminMu1 = deltaR1;
        mu1match = i;
        hlt1 = tto->fLabel;
	if (pt->fPlab.Mag() > 0.) {
	  fTrigMatchDeltaPt = TMath::Abs(tto->fP.Rho() - pt->fPlab.Mag())/pt->fPlab.Mag(); 
	} else {
	  fTrigMatchDeltaPt = 98.;
	}
      } // if delta 
    } // selected 
    
  } // end for loop 
  
  
  if (fVerbose > 15 || localPrint) 
    cout << "best trigger matching: " << mu1match << " dR: " << deltaRminMu1 << " "<<hlt1<<endl;
  
  ((TH1D*)fHistDir->Get("test2"))->Fill(deltaRminMu1); 
  
  if ( mu1match>=0 && deltaRminMu1<deltaRthrsh ) {
    if(localPrint) cout<<" matched "<<hlt1<<endl;
    ((TH1D*)fHistDir->Get("test1"))->Fill(21.); 
    if (     hlt1 == "hltDisplacedmumuFilterDoubleMu4Jpsi:HLT::" ) ((TH1D*)fHistDir->Get("test1"))->Fill(22.); 
    else if (hlt1 == "hltVertexmumuFilterBs345:HLT::")   ((TH1D*)fHistDir->Get("test1"))->Fill(23.);
    else if (hlt1 == "hltVertexmumuFilterBs3p545:HLT::") ((TH1D*)fHistDir->Get("test1"))->Fill(24.);
    else if (hlt1 == "hltVertexmumuFilterBs47:HLT::")    ((TH1D*)fHistDir->Get("test1"))->Fill(25.);
  }

  return deltaRminMu1;
}

// ----------------------------------------------------------------------
void candAna::boostGames() {

  if(fpMuon1 == NULL || fpMuon2 == NULL) return; // protection for DSTAR d.k. 15/1/2013

  double gcosTheta, gcosTheta2, rcosTheta, rcosTheta2;
  TVector3 pvec = TVector3(0., 0., 1.);
  if ((fGenBTmi > -1) && (fGenM1Tmi > -1) && (fGenM2Tmi > -1 )) {
    TGenCand *pB(0), *pM1(0), *pM2(0); 

    pB  = fpEvt->getGenTWithIndex(fGenBTmi); 
    pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi); 
    pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi); 
    
    TVector3 boost = pB->fP.Vect();
    boost.SetMag(boost.Mag()/pB->fP.E());
    TLorentzVector pM1Cms = pM1->fP; 
    pM1Cms.Boost(-boost);
    TLorentzVector pM2Cms = pM2->fP; 
    pM2Cms.Boost(-boost);

    //    N = P_{beam} x P_b / | P_{beam} x P_b |
    TVector3 bvec = pB->fP.Vect();
    TVector3 nvec = pvec.Cross(bvec);     
    TLorentzVector nvec4; nvec4.SetXYZM(nvec.X(), nvec.Y(), nvec.Z(), 0); 
    nvec4.Boost(-boost); 

    if (pM1->fQ > 0) {
      gcosTheta = pM1Cms.CosTheta();
      gcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM1Cms.Vect()));
      
    } else {
      gcosTheta = pM2Cms.CosTheta();
      gcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM2Cms.Vect()));
    }

    if (0 == fNGenPhotons) {
      ((TH1D*)fHistDir->Get("gp1cms"))->Fill(pM1Cms.Rho()); 
      ((TH1D*)fHistDir->Get("gp2cms"))->Fill(pM2Cms.Rho()); 
      ((TH1D*)fHistDir->Get("gt1cms"))->Fill(gcosTheta); 
      ((TH1D*)fHistDir->Get("gt2cms"))->Fill(gcosTheta2); 
    } else {
      ((TH1D*)fHistDir->Get("gp1cmsg"))->Fill(pM1Cms.Rho()); 
      ((TH1D*)fHistDir->Get("gp2cmsg"))->Fill(pM2Cms.Rho()); 
      
      ((TH1D*)fHistDir->Get("gt1cmsg"))->Fill(gcosTheta); 
      ((TH1D*)fHistDir->Get("gt2cmsg"))->Fill(gcosTheta2); 
    }
  }

  


  // -- reco version
  TVector3 boost = fpCand->fPlab;
  double eboost  = TMath::Sqrt(fpCand->fPlab*fpCand->fPlab + fpCand->fMass*fpCand->fMass); 
  boost.SetMag(boost.Mag()/eboost); 
  TLorentzVector pM1Cms; pM1Cms.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
  pM1Cms.Boost(-boost);
  TLorentzVector pM2Cms; pM2Cms.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  pM2Cms.Boost(-boost);

  //    N = P_{beam} x P_b / | P_{beam} x P_b |
  TVector3 bvec = fpCand->fPlab;
  TVector3 nvec = pvec.Cross(bvec);     
  TLorentzVector nvec4; 
  //  nvec4.SetXYZM(nvec.X(), nvec.Y(), nvec.Z(), 0); 
  nvec4.SetXYZT(nvec.X(), nvec.Y(), nvec.Z(), 0); 
  nvec4.Boost(-boost); 

  
  if (fMu1Q > 0) {
    rcosTheta  = pM1Cms.CosTheta();
    rcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM1Cms.Vect()));
  } else {
    rcosTheta = pM2Cms.CosTheta();
    rcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM2Cms.Vect()));
  }


  
  ((TH2D*)fHistDir->Get("tvsm"))->Fill(fpCand->fMass, rcosTheta2); 

  ((TH1D*)fHistDir->Get("rp1cms"))->Fill(pM1Cms.Rho()); 
  ((TH1D*)fHistDir->Get("rp2cms"))->Fill(pM2Cms.Rho()); 
  
  ((TH1D*)fHistDir->Get("rt1cms"))->Fill(rcosTheta); 
  ((TH1D*)fHistDir->Get("rt2cms"))->Fill(rcosTheta2); 
  if (fGoodHLT) {
    ((TH1D*)fHistDir->Get("rt3cms"))->Fill(rcosTheta2); 
  }
  ((TH1D*)fHistDir->Get("gt1"))->Fill(rcosTheta-gcosTheta); 
  ((TH1D*)fHistDir->Get("gt2"))->Fill(rcosTheta2-gcosTheta2); 

  if (0)  cout << "muon 1:  p = " << fpMuon1->fPlab.Mag() << " " << pM1Cms.Rho()
	       << " muon 2: p = " << fpMuon2->fPlab.Mag() << " " << pM2Cms.Rho()  
	       << " theta g = " << gcosTheta  << " r = " << rcosTheta 
	       << " theta g2 = " << gcosTheta2  << " r = " << rcosTheta2 
	       << endl;
  

}
//
//-----------------------------------------------------------------------------------
// Loops over all muons, returns dR of the closests muon, excluding the 
// same track muon (if exists)
double candAna::matchToMuon(TAnaTrack *pt, bool skipSame) {
  //bool print = true;
  bool print = false;

  int numMuons = fpEvt->nMuons();
  if(print) cout << "Found " << numMuons << " rec muons in event" << endl;

  TVector3 trackMom = pt->fPlab;  // test track momentum
  int it0 = pt->fIndex;
  if(print) cout<<" check track "<<it0<<endl;

  TVector3 muonMom;
  TAnaMuon * muon = 0;
  //TSimpleTrack * pTrack = 0;
  double ptMuon=0., dRMin=9999.;
  int select = -1;
  for (int it = 0; it<numMuons; ++it) { // loop over muons
    muon = fpEvt->getMuon(it);
    //if(print) muon->dump();

    // check if this is a nice  muon, accept only global and tracker muons
    int muonId = muon->fMuID;
    //     if ( (muonId & 0x6) == 0 ) continue;  // skip muons which are not global/tracker
 
    int itrk = muon->fIndex;
    // Eliminate pure standalone muons and calo muons. Skip same track comparion (only of skipSame=true) 
    if(itrk<0 || (skipSame && itrk == it0)) {
      if(print) {
	if(itrk<0) cout<<"standalone only or calo muon? "<<hex<<muonId<<dec<<" "<<(muon->fPlab).Perp()<<endl;
	else       cout<<"skip, same track "<<endl;
      }
      continue;  // skip same track comparion and standalone/calo muons
    }

    //if(itrk>0) cout<<" tracker muon? "<<hex<<(muonId&0x6)<<dec<<" "<<(muon->fPlab).Perp()<<endl;
  
    // Use direct access, withour going through SimpleTracks
    muonMom = muon->fPlab;
    ptMuon  = muonMom.Perp();
    double dR = muonMom.DeltaR(trackMom);
    //double etaMuon = muonMom.Eta();
    //double phiMuon = muonMom.Phi();
    if(print) cout<<it<<" "<<ptMuon<<" "<<dR<<endl;
   
    // Go through reco track, Find the reco track  NOT NEEDED 
//     if(itrk>=0 && itrk< (fpEvt->nSimpleTracks()) ) {  // if the simple track exists 
//       //pTrack = fpEvt->getRecTrack(itrk);
//       pTrack = fpEvt->getSimpleTrack(itrk);
//       //cout<<it<<" "<<itrk<<" "<<pTrack<<endl;
//       if(pTrack != 0) {
// 	muonMom = pTrack->getP();
//  	ptMuon  = muonMom.Perp();
// 	dR = muonMom.DeltaR(trackMom);
// 	if(print) cout<<it<<" "<<ptMuon<<" "<<dR<<endl;
//       } //if ptrack
//     } // if track

    if(dR<dRMin) {dRMin=dR; select=it;} // select the best fit

  } // loop over muons

  if(select>-1) {
    int idx =  (fpEvt->getMuon(select))->fIndex;  // find the track index 
    if(print) cout<<"final muon match "<<dRMin<<" "<<select<<" "<<idx<<endl;
    
    ((TH1D*)fHistDir->Get("test7"))->Fill(dRMin);
  }

  return dRMin;
}



// ----------------------------------------------------------------------
void candAna::play() {

}
