#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFTruthCandidate.h"

#include <algorithm>

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TSimpleTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>


#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

using reco::Track;
using reco::Vertex;
using reco::TrackBaseRef;

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace edm;
using namespace std;
using namespace reco;

// ----------------------------------------------------------------------
HFTruthCandidate::HFTruthCandidate(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", string("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<InputTag>("BeamSpotLabel", InputTag("offlineBeamSpot"))),
  fMuonsLabel(iConfig.getUntrackedParameter<edm::InputTag>("muonsLabel", edm::InputTag("muons"))),
  fPartialDecayMatching(iConfig.getUntrackedParameter<bool>("partialDecayMatching", false)), 
  fMotherID(iConfig.getUntrackedParameter("motherID", 0)), 
  fType(iConfig.getUntrackedParameter("type", 67)),
  fGenType(iConfig.getUntrackedParameter("GenType", -67)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 99.)),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)) {

  vector<int> defaultIDs;
  defaultIDs.push_back(0);
  fDaughtersID = iConfig.getUntrackedParameter<vector<int> >("daughtersID", defaultIDs);
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTruthCandidate constructor" << endl;
  cout << "--- verbose:               " << fVerbose << endl;
  cout << "--- tracksLabel:           " << fTracksLabel << endl;
  cout << "--- PrimaryVertexLabel:    " << fPrimaryVertexLabel << endl;
  cout << "--- BeamSpotLabel:         " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel            " << fMuonsLabel << endl;
  cout << "--- motherID:              " << fMotherID << endl;
  cout << "--- type:                  " << fType << endl;
  cout << "--- GenType:               " << fGenType << endl;
  cout << "--- partialDecayMatching:  " << fPartialDecayMatching << endl;
  cout << "--- maxDoca:               " << fMaxDoca << endl;
  fDaughtersSet.clear(); 
  fStableDaughters = 0; 
  for (unsigned int i = 0; i < fDaughtersID.size(); ++i) {
    cout << "---   daughterID:              " << fDaughtersID[i] << endl;
    if (TMath::Abs(fDaughtersID[i]) == 11)   ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 13)   ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 211)  ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 321)  ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 2212) ++fStableDaughters; 
    fDaughtersSet.insert(TMath::Abs(fDaughtersID[i])); 
    fDaughtersGammaSet.insert(TMath::Abs(fDaughtersID[i])); 
    fDaughtersGamma2Set.insert(TMath::Abs(fDaughtersID[i])); 
  }    
  cout << "---    total stable particles: " << fStableDaughters << endl;
  fDaughtersGammaSet.insert(22); 
  fDaughtersGamma2Set.insert(22); 
  fDaughtersGamma2Set.insert(22); 
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
HFTruthCandidate::~HFTruthCandidate() {  
}


// ----------------------------------------------------------------------
void HFTruthCandidate::beginJob() {
}

// ----------------------------------------------------------------------
void HFTruthCandidate::endJob() {
}


// ----------------------------------------------------------------------
void HFTruthCandidate::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- load the muons
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(fMuonsLabel, muonHandle);
  const reco::MuonCollection *muonCollection(0);
  if (muonHandle.isValid()) {
    muonCollection = muonHandle.product();
  }

  // -- In generator block, find mother with declared decay channel
  multiset<int> genDaughters; 
  multiset<int> genIndices; 
  multiset<pair<int, int> > genMap; 
  TGenCand *pGen, *pDau, *pTmp;
  int matchedDecay(0);
  int iMom(-1), motherIndex(-1); 

  vector<int> bla(100); 
  vector<int>::iterator blaIt; 

  //   cout << "----------------------------------------------------------------------" << endl;
  //  cout << " ngenCands: " << gHFEvent->nGenCands() << endl;
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (TMath::Abs(pGen->fID) == fMotherID) {
      motherIndex = ig; 
      if (fVerbose > 1) {
	cout << "mother ";
	pGen->dump(); 
      }
      genDaughters.clear(); 
      genIndices.clear();
      genMap.clear();

      // -- version with descendants
      for (int id = ig+1; id < gHFEvent->nGenCands(); ++id) {
	pDau = gHFEvent->getGenCand(id);
	iMom = pDau->fMom1;
	while (iMom > ig) {
	  pTmp = gHFEvent->getGenCand(iMom);
	  iMom = pTmp->fMom1;
	}
	if (iMom == ig) {
	  if (fVerbose > 1) {
	    cout << "  daug: ";
	    pDau->dump(); 
	  }
	  if (fPartialDecayMatching) {
	    if (fDaughtersSet.find(TMath::Abs(pDau->fID)) != fDaughtersSet.end()) {
	      genDaughters.insert(TMath::Abs(pDau->fID)); 
	      genIndices.insert(id); 
	      genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
	    }
	  } else {
	    genDaughters.insert(TMath::Abs(pDau->fID)); 
	    genIndices.insert(id); 
	    genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
	  }
	}
      }

      // -- now check whether this is PARTIALLY the decay channel in question
      if (fPartialDecayMatching) {
	blaIt = set_intersection(genDaughters.begin(), genDaughters.end(), fDaughtersSet.begin(), fDaughtersSet.end(), bla.begin()); 
	if (static_cast<unsigned int>(blaIt - bla.begin()) == fDaughtersSet.size()) {
	  matchedDecay = 1; 
	  if (fVerbose > 0) {
	    cout << "matched partial decay: ";
	    for (vector<int>::iterator it = bla.begin(); it != blaIt; ++it) cout << *it << " "; 
	    cout << endl;
	  }
	  break;
	}
      }

      // -- now check whether this is the decay channel in question
      if (fDaughtersSet == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay" << endl;
	break;
      }
      if (fDaughtersGammaSet == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay with bremsstrahlung photon" << endl;
	break;
      }
      if (fDaughtersGamma2Set == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay with 2 bremsstrahlung photons" << endl;
	break;
      }
    }
  }


  // -- Dump generator candidate made from stable charged particles
  if (matchedDecay > 0) {
    int id(-1), idx(-1); 
    TLorentzVector comp; 
    for (multiset<pair<int, int> >::iterator i = genMap.begin(); i != genMap.end(); ++i) {
      idx = i->first; 
      id  = i->second; 
      if (id == 11 || id == 13 || id == 211 || id == 321 || id ==2212)  {
	comp += gHFEvent->getGenCand(idx)->fP ; 
      }
    }

    TAnaCand *pCand = gHFEvent->addCand();
    pCand->fPlab = comp.Vect();
    pCand->fMass = comp.M();
    pCand->fType = fGenType;
    pCand->fIndex= motherIndex;
    if (fVerbose > 1) {
      char line[200];
      sprintf(line, "p=%8.3f(%+9.3f,%+9.3f,%+9.3f), mass = %f", 
	      pCand->fPlab.Mag(), 
	      pCand->fPlab.X(), pCand->fPlab.Y(), pCand->fPlab.Z(), 
	      pCand->fMass);
      cout << line << endl;
    }

  }


  if (fVerbose > 2 && 0 == matchedDecay)  {
    cout << "Did not match decay" << endl;
    for (multiset<int>::iterator i = genDaughters.begin(); i != genDaughters.end(); ++i) {
      cout << " unmatched genDaughter: " << *i << endl;
    }
  }

  // -- Construct and dump reconstructed candidates matched to generator particles
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFTruthCandidate> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  Vertex dummy; 
  HFKalmanVertexFit  aKal(0, dummy, 1, fType); 
  map<int,int> trackIxMap; // map: genIx -> trkIx
  map<int,double> trackMassesMap; // map: genIx -> mass

  TSimpleTrack *pTrack; 
  if (matchedDecay > 0) {
    for (int it = 0; it < gHFEvent->nSimpleTracks(); ++it) {
      pTrack = gHFEvent->getSimpleTrack(it); 
      int gidx = pTrack->getGenIndex(); 
      if (genIndices.find(gidx) != genIndices.end()) {
	if (fVerbose > 2) {
	  cout << "Found simple track: " << it; 
	  pTrack->dump(); 
	}
	
	double mass = MMUON; 
	int mcid = TMath::Abs(gHFEvent->getGenCand(gidx)->fID); 
	if (321  == mcid) mass = MKAON;
	if (211  == mcid) mass = MPION;
	if (13   == mcid) mass = MMUON;
	if (2212 == mcid) mass = MPROTON;
				
	if (trackIxMap.count(gidx) > 0) {
	  // -- this code is not longer called since the introduction of the truth matching cleanup
	  ESHandle<MagneticField> magfield;
	  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
	  AnalyticalImpactPointExtrapolator ipExt(magfield.product());
	  GlobalPoint vtx(0,0,0);
	  FreeTrajectoryState fts;
	  TrajectoryStateOnSurface tsof;
	  TrackBaseRef trackView;
	  TVector3 ipGen,ipThis,ipOld;
					
	  // GENERATOR IMPACT POINT
	  pGen = gHFEvent->getGenCand(gidx);
	  fts = FreeTrajectoryState(GlobalPoint(pGen->fV.X(),pGen->fV.Y(),pGen->fV.Z()),
				    GlobalVector(pGen->fP.X(),pGen->fP.Y(),pGen->fP.Z()),
				    TrackCharge(pGen->fQ),
				    magfield.product());
	  tsof = ipExt.extrapolate(fts,vtx);
	  ipGen.SetXYZ(tsof.globalPosition().x(), tsof.globalPosition().y(), tsof.globalPosition().z());
					
	  // NEW TRACK IMPACT POINT
	  trackView = TrackBaseRef(hTracks,it);
	  fts = FreeTrajectoryState(GlobalPoint(trackView->vx(),trackView->vy(),trackView->vz()),
				    GlobalVector(trackView->px(),trackView->py(),trackView->pz()),
				    trackView->charge(),
				    magfield.product());
	  tsof = ipExt.extrapolate(fts,vtx);
	  ipThis.SetXYZ(tsof.globalPosition().x(), tsof.globalPosition().y(), tsof.globalPosition().z());
					
	  // OLD TRACK IMPACT POINT
	  trackView = TrackBaseRef(hTracks,trackIxMap[gidx]);
	  fts = FreeTrajectoryState(GlobalPoint(trackView->vx(),trackView->vy(),trackView->vz()),
				    GlobalVector(trackView->px(),trackView->py(),trackView->pz()),
				    trackView->charge(),
				    magfield.product());
	  tsof = ipExt.extrapolate(fts,vtx);
	  ipOld.SetXYZ(tsof.globalPosition().x(), tsof.globalPosition().y(), tsof.globalPosition().z());
					
	  // compare the better tracks
	  if ( (ipGen - ipThis).Mag() < (ipGen - ipOld).Mag() ) {
	    trackIxMap[gidx] = it;
	    trackMassesMap[gidx] = mass;
	    if (fVerbose > 0) cout << "-> overwriting to trackList: " << it << " with ID = " << mcid << endl;
	  } else {
	    if (fVerbose > 2) cout << "-> rejecting to trackList: " << it << " with ID = " << mcid << endl;
	  }
	} else {
	  trackIxMap.insert(make_pair(gidx, it));
	  trackMassesMap.insert(make_pair(gidx, mass));
	  if (fVerbose > 2) cout << "-> adding to trackList: " << it << " with ID = " << mcid << endl;
	}
      }
    }


    if (static_cast<int>(trackIxMap.size()) == fStableDaughters) {
			
      vector<Track> trackList;
      vector<int> trackIndices;
      vector<double> trackMasses;
      for (map<int,int>::const_iterator it = trackIxMap.begin(); it != trackIxMap.end(); ++it) {
	TrackBaseRef trackView(hTracks,it->second);
	Track track(*trackView);
	trackList.push_back(track);
	trackIndices.push_back(it->second);
	trackMasses.push_back(trackMassesMap[it->first]);
      }
      aKal.doNotFit(trackList, trackIndices, trackMasses, fType); 

      // -- Vertexing, with Kinematic Particles!
      ESHandle<MagneticField> magfield;
      iSetup.get<IdealMagneticFieldRecord>().get(magfield);
      const MagneticField *field = magfield.product();

      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
      if (!fTTB.isValid()) {
	cout << " -->HFTruthCandidate: Error: no TransientTrackBuilder found."<<endl;
	return;
      }

      Handle<VertexCollection> recoPrimaryVertexCollection;
      iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
      if(!recoPrimaryVertexCollection.isValid()) {
	cout << "==>HFTruthCandidate> No primary vertex collection found, skipping" << endl;
	return;
      }
      const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
      if (vertices.size() == 0) {
	cout << "==>HFTruthCandidate> No primary vertex found, skipping" << endl;
	return;
      }
	  
      // -- get the beamspot
      Handle<BeamSpot> bspotHandle;
      iEvent.getByLabel(fBeamSpotLabel, bspotHandle);
      if (!bspotHandle.isValid()) {
	cout << "==>HFbs2JpsiPhi> No beamspot found, skipping" << endl;
	return;
      }
      BeamSpot bspot = *bspotHandle;		

      HFSequentialVertexFit aSeq(hTracks, muonCollection, fTTB.product(), recoPrimaryVertexCollection, field, bspot, fVerbose);
      // -- setup with (relevant) muon hypothesis
      HFDecayTree theTree(1000000+fType, true, 0, false); 
      int ID(0), IDX(0); 
      for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	IDX = trackIndices[ii];
	ID  = 13;
	if (fVerbose > 2) cout << "-> adding track " << IDX << " with ID = " << ID << " to the tree" << endl; 
	theTree.addTrack(IDX, ID);
      }
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      aSeq.doFit(&theTree);


      // -- setup with (correct) truth hypothesis
      HFDecayTree theTree2(2000000+fType, true, 0, false);
      for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	IDX = trackIndices[ii];
	ID  = gHFEvent->getSimpleTrackMCID(IDX);
	if (fVerbose > 2) cout << "-> adding track " << IDX << " with ID = " << ID << " to the tree" << endl; 
	theTree2.addTrack(IDX, ID);
      }
      theTree2.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      aSeq.doFit(&theTree2);

      // -- special case for the normalization sample
      if (68 == fType || 66 == fType) {
	HFDecayTree theTree3(3000000 + fType, true, MBPLUS, false, -1.0, true);

	HFDecayTreeIterator iterator = theTree3.addDecayTree(300443, false, MJPSI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  =  gHFEvent->getSimpleTrackMCID(IDX);
	  if (13 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 13);
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (321 == TMath::Abs(ID)) {
	    theTree3.addTrack(IDX, 321);
	  }
	  if (211 == TMath::Abs(ID)) {
	    theTree3.addTrack(IDX, 321); // assign the kaon mass!
	  }
	}
	theTree3.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit for Bu2JpsiKp" << endl;
	aSeq.doFit(&theTree3);
      }
			
      // -- special case for the control sample
      if (67 == fType) {
	HFDecayTree theTree4(3000000 + fType, true, MBS, false, -1.0, true);

	HFDecayTreeIterator iterator = theTree4.addDecayTree(300443, false, MJPSI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (13 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 13);
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	iterator = theTree4.addDecayTree(300333, false, MPHI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (321 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 321);
	  }
	  if (211 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 321); // assign the kaon mass! do we ever get here?
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
	theTree4.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	if (fVerbose > 5) cout << "==>HFTruthCandidate> sequential fit for Bs2JpsiPhi" << endl;
	aSeq.doFit(&theTree4);
      }

      // -- special case for B0 -> J/psi Kstar reconstructed as B+ -> J/psi K (leaving out the pion)
      if (10 == fType) {
	HFDecayTree theTree5(3000000 + fType, true, MBPLUS, false, -1.0, true);

	HFDecayTreeIterator iterator = theTree5.addDecayTree(300443, false, MJPSI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (13 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 13);
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (321 == TMath::Abs(ID)) {
		  theTree5.addTrack(IDX, 321);
	  }
		// Comment to leave out the pion
//	  if (211 == TMath::Abs(ID)) {
//	    //do nothing and ignore the pion from the K*!
//		  theTree5.addTrack(IDX, 321); // assign the kaon mass!
//	  }
	}
	theTree5.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	if (fVerbose > 5) cout << "==>HFTruthCandidate> sequential fit for Bu2JpsiK of Bd2JpsiKstar" << endl;
	aSeq.doFit(&theTree5);
      }


      // -- special case for B0 -> J/psi Kstar reconstructed as Bs -> J/psi KK 
      if (11 == fType) {
	HFDecayTree theTree6(3000000 + fType, true, MBS, false, -1.0, true);

	HFDecayTreeIterator iterator = theTree6.addDecayTree(300443, false, MJPSI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (13 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 13);
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	iterator = theTree6.addDecayTree(300333, false, MPHI, false);
	for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
	  IDX = trackIndices[ii];
	  ID  = gHFEvent->getSimpleTrackMCID(IDX);
	  if (321 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 321);
	  }
	  if (211 == TMath::Abs(ID)) {
	    iterator->addTrack(IDX, 321);
	  }
	}
	iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
	theTree6.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

	if (fVerbose > 5) cout << "==>HFTruthCandidate> sequential fit for Bs2JpsiPhi of Bd2JpsiKstar" << endl;
	aSeq.doFit(&theTree6);
      }

      // -- special case for B0 -> D*-pi+
      int iPion1(-1), iPion2(-1), iPion(-1), iKaon(-1); 
      if (30 == fType) {
        for (unsigned int ii = 0; ii < trackIndices.size(); ++ii) {
          IDX = trackIndices[ii];
          ID  = gHFEvent->getSimpleTrackMCID(IDX);
          int mcidx   = gHFEvent->getSimpleTrack(IDX)->getGenIndex();
          int mcmoidx = gHFEvent->getGenCand(mcidx)->fMom1;
          int mcmoID  = gHFEvent->getGenCand(mcmoidx)->fID;
          //      cout << "IDX = " << IDX << " ID = " << ID << " mcmoID = " << mcmoID << endl;
          if (421 == TMath::Abs(mcmoID)) {
            if (321 == TMath::Abs(ID)) iKaon = IDX; 
            if (211 == TMath::Abs(ID)) iPion = IDX; 
          } else if (413 == TMath::Abs(mcmoID)) {
            iPion1 = IDX;
          } else {
            iPion2 = IDX; 
          }
        }

        if (iPion > -1 && iKaon > -1 && iPion1 > -1 && iPion2 > -1) {
          
          HFDecayTree theTree(3000000 + fType, true, MB_0, false, -1.0, true);
          theTree.addTrack(iPion2,211); // add the fast pion to the B0 candidate
          theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
          
          HFDecayTreeIterator iterator = theTree.addDecayTree(300031, false, MDSTARPLUS, false); // D*-
          iterator->addTrack(iPion1,211); // add the slow pion to the D*+ candidate
          iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
          
          iterator = iterator->addDecayTree(300032, true, MD0, false); // D0
          iterator->addTrack(iKaon,321);
          iterator->addTrack(iPion,211);
          iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
          
          aSeq.doFit(&theTree);
        }
      }

    }

    // -- now copy reduced part of gen block and clear the original
    gHFEvent->fillGenT(motherIndex); 
  }

  


}

//define this as a plug-in
DEFINE_FWK_MODULE(HFTruthCandidate);
