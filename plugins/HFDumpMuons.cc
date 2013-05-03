#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpMuons.h"

#include <iostream>
#include <cmath>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpUtilities.hh"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;

// -- sort the vector with xpTracks
static bool dist_less(const xpTrack &x, const xpTrack &y) {
  return (x.dist < y.dist); 
}


// ----------------------------------------------------------------------
HFDumpMuons::HFDumpMuons(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("ctfWithMaterialTracks"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fBeamSpotLabel(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpotLabel", edm::InputTag("offlineBeamSpot"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertexLabel", edm::InputTag("offlinePrimaryVertices"))),
  fCaloMuonsLabel(iConfig.getUntrackedParameter<InputTag>("calomuonsLabel")),
  fMaxTrackDistToStore(iConfig.getUntrackedParameter<double>("maxTrackDist",0.1)),
  fDocaVertex(iConfig.getUntrackedParameter<double>("docaVertex",0.05)),
  fKeepBest(iConfig.getUntrackedParameter<int>("keepBest",3)),
  fMaxCandTracks(iConfig.getUntrackedParameter<int>("maxCandTracks",3)),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fRunOnAOD(iConfig.getUntrackedParameter<bool>("runOnAOD",false)),
  fpropM1(iConfig.getParameter<edm::ParameterSet>("propM1")),
  fpropM2(iConfig.getParameter<edm::ParameterSet>("propM2"))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpMuons constructor" << endl;
  cout << "---  fTracksLabel             " << fTracksLabel.encode() << endl;
  cout << "---  fMuonsLabel:             " << fMuonsLabel.encode() << endl;
  cout << "---  BeamSpotLabel            " << fBeamSpotLabel << endl;
  cout << "---  PrimaryVertexLabel       " << fPrimaryVertexLabel << endl;
  cout << "---  fCaloMuonsLabel:         " << fCaloMuonsLabel.encode() << endl;
  cout << "---  fDoTruthMatching:        " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "---  fVerbose:                " << fVerbose << endl;
  cout << "---  fRunOnAOD:               " << fRunOnAOD << endl;
  cout << "---  fMaxTrackDistToStore:    " << fMaxTrackDistToStore << endl;
  cout << "---  docaVertex:              " << fDocaVertex << endl;
  cout << "---  keepBest:                " << fKeepBest << endl;
  cout << "---  maxCandTracks:           " << fMaxCandTracks << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpMuons::~HFDumpMuons() {

}

void HFDumpMuons::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  fpropM1.init(iSetup);
  fpropM2.init(iSetup);
}

// ----------------------------------------------------------------------
void HFDumpMuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // -- load the beam spot
  edm::Handle<reco::BeamSpot> bspotHandle;
  iEvent.getByLabel(fBeamSpotLabel, bspotHandle);
  fBeamSpot = 0; 
  if (bspotHandle.isValid()) fBeamSpot = bspotHandle.product();
  
  // -- load the primary vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(fPrimaryVertexLabel, vertexHandle);
  fVertexCollection = 0; 
  if (vertexHandle.isValid()) {
    fVertexCollection = vertexHandle.product();
  }


  // -- tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if (hTracks.isValid()) {
    fhTracks = &hTracks; // to be used in fillMuon()
  } else {
    cerr << "==> HFDumpMuons> ERROR loading the tracks" << endl;
    throw std::string("==> HFDumpMuons> ERROR loading the tracks");
    fhTracks = NULL;
  }
  
  // Load the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cerr << "==> HFDumpMuons> ERROR loading the transient track builder" << endl;
    throw std::string("==> HFDumpMuons> ERROR loading the transient track builder");
  }

  extrapolateTracks();
  
  // -- global muons
  Handle<MuonCollection> hMuons;
  if (fVerbose > 0) cout << "==> HFDumpMuons> " << fMuonsLabel << endl;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  fMuonCollection = hMuons.product();

  int im(0);
  for (reco::MuonCollection::const_iterator iMuon = hMuons->begin(); iMuon != hMuons->end();  iMuon++) {
    fillMuon(*iMuon, im); 
    ++im;
  }

  // -- calo muons
  if(!fRunOnAOD) {
    Handle<CaloMuonCollection> cMuons;
    if (fVerbose > 0) cout << "==> HFDumpMuons> " << fCaloMuonsLabel << endl;
    iEvent.getByLabel(fCaloMuonsLabel, cMuons);
    
    for (reco::CaloMuonCollection::const_iterator cMuon = cMuons->begin(); cMuon != cMuons->end();  cMuon++) {
      fillCaloMuon(*cMuon, im); 
      ++im;
    }
  }

  if (fVerbose > 0) {
    for (int im = 0; im < gHFEvent->nMuons(); ++im) {
      gHFEvent->getMuon(im)->dump();
    }
  }
  
  // make sure the pointer does not point to
  // invalid location
  fhTracks = NULL;
}


// ----------------------------------------------------------------------
void HFDumpMuons::fillMuon(const reco::Muon& rm, int im) {

  TrackRef gTrack = rm.globalTrack();
  TrackRef iTrack = rm.innerTrack();
  TrackRef oTrack = rm.outerTrack();

  TAnaMuon *pM = gHFEvent->addMuon();    

  
  if (rm.innerTrack().isNonnull()) {
    Track trk(*iTrack);
    fillAnaTrack(pM, trk, rm.innerTrack().index(), -2, fVertexCollection, fMuonCollection, fBeamSpot); 
  } else {
    pM->fIndex = -23;
  }
  pM->fMuIndex = im;
  pM->fMuID    = muonID(rm);
  pM->fQ       = rm.charge();

  pM->fMuonChi2         = rm.combinedQuality().trkKink;
  pM->fTimeInOut        = rm.time().timeAtIpInOut; 
  pM->fTimeInOutE       = rm.time().timeAtIpInOutErr; 
  pM->fTimeOutIn        = rm.time().timeAtIpOutIn; 
  pM->fTimeOutInE       = rm.time().timeAtIpOutInErr; 
  pM->fTimeNdof         = rm.time().nDof;
  pM->fNmatchedStations = rm.numberOfMatchedStations();

  bool isGlobalMuon = muon::isGoodMuon(rm, muon::AllGlobalMuons);

  // -- variables for MVA muon ID
  if (gTrack.isNonnull() && iTrack.isNonnull()) {
    const HitPattern track_hp  = iTrack->hitPattern();
    const HitPattern exp_track_out_hp = iTrack->trackerExpectedHitsOuter();
    reco::MuonQuality muQuality = rm.combinedQuality();
    
    pM->fItrkValidFraction       = iTrack->validFraction(); //1
    pM->fGtrkNormChi2            = gTrack->normalizedChi2(); //2
    pM->fChi2LocalPosition       = muQuality.chi2LocalPosition; //3
    pM->fNumberOfLostTrkHits     = exp_track_out_hp.numberOfLostTrackerHits(); //4
    pM->fSegmentComp             = muon::segmentCompatibility(rm); //5
    pM->fGtrkProb                = muQuality.glbTrackProbability; //6
    pM->fChi2LocalMomentum       = muQuality.chi2LocalMomentum; //6
    pM->fNumberOfValidTrkHits    = track_hp.numberOfValidTrackerHits(); //7
  }


  if (gTrack.isNonnull()) {
    Track trk(*gTrack);
    pM->fGlobalPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    pM->fNvalidMuonHits = gTrack->hitPattern().numberOfValidMuonHits();

    if (!fRunOnAOD) {
      vector<unsigned int> hits = muonStatHits(trk);
      pM->fNhitsDT  = hits.at(0); 
      pM->fNhitsCSC = hits.at(1); 
      pM->fNhitsRPC = hits.at(2); 
    } else {
      pM->fNhitsDT  = -1;
      pM->fNhitsCSC = -1;
      pM->fNhitsRPC = -1;
    }
  }

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    //    cout << " ipt = " << trk.pt() << " " <<  trk.eta() << " "  <<  trk.phi() << endl;
  }

  if (oTrack.isNonnull()) {
    Track trk(*oTrack);
    pM->fOuterPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  }

  // -- propagate muons to muon system to get their impact point
  TVector3 muPosM1; 
  bool validM1(false); 

  if (isGlobalMuon && doExtrapolate(rm.pt(), rm.eta())) {
    TrajectoryStateOnSurface prop_M1 = fpropM1.extrapolate(rm);
    TrajectoryStateOnSurface prop_M2 = fpropM2.extrapolate(rm);
    
    if (prop_M1.isValid()) {
      pM->fPositionAtM1.SetXYZ(prop_M1.globalPosition().x(), prop_M1.globalPosition().y(), prop_M1.globalPosition().z());
    }
    if (prop_M2.isValid()) {
      pM->fPositionAtM2.SetXYZ(prop_M2.globalPosition().x(), prop_M2.globalPosition().y(), prop_M2.globalPosition().z());
    }
    
    if (oTrack.isNonnull()) {
      TrajectoryStateOnSurface propOuter = fpropM1.extrapolate(*oTrack);
      if (propOuter.isValid()) {
	validM1 = true; 
	muPosM1.SetXYZ(propOuter.globalPosition().x(),propOuter.globalPosition().y(),propOuter.globalPosition().z());
	pM->fMuonTrackPosAtM1.SetXYZ(propOuter.globalPosition().x(),propOuter.globalPosition().y(),propOuter.globalPosition().z());
	pM->fMuonTrackPlabAtM1.SetXYZ(propOuter.globalMomentum().x(),propOuter.globalMomentum().y(),propOuter.globalMomentum().z());
      }
    }
  }

  // -- doca of close tracks to muon
  if (isGlobalMuon && iTrack.isNonnull() && fTTB.isValid()) {
    TwoTrackMinimumDistance md(TwoTrackMinimumDistance::SlowMode);
    Track trkMuon(*iTrack);
    TransientTrack transTrkMuon = fTTB->build(trkMuon);
	  
    for (size_t k = 0; k < (*fhTracks)->size(); k++) {
      if (k == iTrack.index()) continue; // own track
		  
      TrackBaseRef bRefTrk(*fhTracks,k);
      Track trk(*bRefTrk);
      TransientTrack transTrk = fTTB->build(trk);
		  
      md.calculate(transTrkMuon.initialFreeState(),transTrk.initialFreeState());
      if (md.distance() < fMaxTrackDistToStore) {
	pM->fNstTracks.insert(std::make_pair(k,md.distance()));
      }
    }
  }


  if (isGlobalMuon && validM1) {
    vector<xpTrack> xvec; 
    for (unsigned int i = 0; i < fXpTracks.size(); ++i) {
      xpTrack x = fXpTracks[i]; 
      if (x.idx == static_cast<int>(iTrack.index())) continue;
      x.dist = (muPosM1 - x.r).Mag();  
      xvec.push_back(x); 
    }
    
    // -- sort the vector & keep only the first TAnaMuon::NXPTRACKS
    sort(xvec.begin(), xvec.end(), dist_less);
    if (TAnaMuon::NXPTRACKS < xvec.size()) {
      for (unsigned int ii = 0; ii < TAnaMuon::NXPTRACKS; ++ii) pM->fXpTracks[ii] = xvec[ii]; 
    } else {
      for (unsigned int ii = 0; ii < xvec.size(); ++ii) pM->fXpTracks[ii] = xvec[ii]; 
    }
  }

  // do the muon vertex analysis
  if (isGlobalMuon && iTrack.isNonnull()) {
    std::set<unsigned> trks;
    double prob = -1.0;
    trks.insert(iTrack.index());
    findVertex(pM,&trks,&prob);
    pM->fVtxProb = prob;
    pM->fVtxTracks = trks;
  }
}

// ----------------------------------------------------------------------
void HFDumpMuons::fillCaloMuon(const reco::CaloMuon& rm, int im) {

  TAnaMuon *pM = gHFEvent->addMuon();    
  pM->fMuID    = 0;  // this assumes that fillCaloMuon is independent from the main muons d.k.

  if (rm.innerTrack().isNonnull()) {
    pM->fIndex = rm.innerTrack().index();
    pM->fQ        = rm.charge();
  } else {
    pM->fIndex = -23;
    pM->fQ     = 0;
  }
  pM->fMuIndex = im; 
  pM->fMuID   |= 0x1<<15;

  pM->fNhitsDT  = 0; 
  pM->fNhitsCSC = 0; 
  pM->fNhitsRPC = 0; 
  

  TrackRef iTrack = rm.innerTrack();

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  } 

}

void HFDumpMuons::findVertex(TAnaMuon *anaMu, std::set<unsigned> *trkIcs, double *prob)
{
	std::vector<TransientTrack> transTracks;
	std::vector<std::pair<double,unsigned> > bestTracks;
	std::map<int,float>::const_iterator mapIt;
	std::set<unsigned>::const_iterator it;
	std::set<unsigned> resultIcs;
	KalmanVertexFitter kvf;
	double best;
	unsigned ix;
	
	// build the transient tracks with 'trkIcs'
	for (it = trkIcs->begin(); it != trkIcs->end(); ++it) {
		TrackBaseRef baseRef(*fhTracks,*it);
		Track trk(*baseRef);
		TransientTrack ttrack = fTTB->build(trk);
		transTracks.push_back(ttrack);
	}
	
	for (mapIt = anaMu->fNstTracks.begin(); mapIt != anaMu->fNstTracks.end(); ++mapIt) {
		
		if (mapIt->second >= fDocaVertex)
			continue;
		
		if (trkIcs->count(mapIt->first) > 0)
			continue; // already included
		
		trkIcs->insert(mapIt->first);
		
		TrackBaseRef baseRef(*fhTracks,mapIt->first);
		Track trk(*baseRef);
		TransientTrack ttrack = fTTB->build(trk);
		transTracks.push_back(ttrack);
		
		TransientVertex vtx = kvf.vertex(transTracks);
		ChiSquared chi(vtx.totalChiSquared(), vtx.degreesOfFreedom());
		best = chi.probability();
		if (!TMath::IsNaN(best))
			bestTracks.push_back(make_pair(chi.probability(),mapIt->first));
		
		trkIcs->erase(mapIt->first);
		transTracks.pop_back();
	}
	
	// only iterate the most promosing 'keep'
	std::sort(bestTracks.begin(),bestTracks.end());
	if (bestTracks.size() > fKeepBest) bestTracks.erase(bestTracks.begin(),bestTracks.end()-fKeepBest);
	
	best = *prob;
	resultIcs = *trkIcs;
	for (ix = 0; ix < bestTracks.size(); ix++) {
		
		std::set<unsigned> curTracks = *trkIcs;
		double result = bestTracks[ix].first;
		
		curTracks.insert(bestTracks[ix].second);
		
		if (curTracks.size() < fMaxCandTracks)
			findVertex(anaMu,&curTracks,&result);
		
		if (best < result) {
			best = result;
			resultIcs = curTracks;
		}
	}
	
	// save
	*prob = best;
	*trkIcs = resultIcs;
} // findVertex()

// ----------------------------------------------------------------------
vector<unsigned int> HFDumpMuons::muonStatHits(const reco::Track& tr) {
  vector<unsigned int> theMuonHits;
  unsigned int nRecHitDT(0), nRecHitCSC(0), nRecHitRPC(0);

  for (trackingRecHit_iterator recHit = tr.recHitsBegin(); recHit != tr.recHitsEnd(); ++recHit){
    DetId detIdHit = (*recHit)->geographicalId();
    if (detIdHit.det() == DetId::Muon ){
      if (detIdHit.subdetId() == MuonSubdetId::DT ) nRecHitDT++;
      else if (detIdHit.subdetId() == MuonSubdetId::CSC ) nRecHitCSC++;
      else if (detIdHit.subdetId() == MuonSubdetId::RPC ) nRecHitRPC++;
    }
  }
  
  theMuonHits.push_back(nRecHitDT); 
  theMuonHits.push_back(nRecHitCSC);
  theMuonHits.push_back(nRecHitRPC);
  return theMuonHits; 
}

// ----------------------------------------------------------------------
void HFDumpMuons::extrapolateTracks() {

  fXpTracks.clear(); 

  for (unsigned int k = 0; k < (*fhTracks)->size(); ++k) {
    
    TrackBaseRef bRefTrk(*fhTracks, k);
    Track trk(*bRefTrk);

    double eta = TMath::Abs(trk.eta()); 
    double pt  = trk.pt();
    
    // -- skip regions where extrapolation is futile
    if (!doExtrapolate(pt, eta)) continue;

    TrajectoryStateOnSurface tsos = fpropM1.extrapolate(trk);
    if (tsos.isValid()) {
      xpTrack x; 
      x.idx = k; 
      x.r = TVector3(tsos.globalPosition().x(),tsos.globalPosition().y(),tsos.globalPosition().z()); 
      x.p = TVector3(tsos.globalMomentum().x(),tsos.globalMomentum().y(),tsos.globalMomentum().z());
      x.dist = 9999.;
      fXpTracks.push_back(x); 
    }
  }

}


// ----------------------------------------------------------------------
bool HFDumpMuons::doExtrapolate(double pt, double eta) { 
  if (pt < 0.8) return false;
  if (eta > 2.4) return false;
  if (eta < 0.8 && pt < 3.4) return false;
  if (eta > 0.8 && eta < 1.4 && (pt < -eta + 4.2)) return false;
  if (eta > 1.4 && eta < 2.4 && (pt < -eta + 2.8)) return false;
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpMuons::beginJob() {

  gHFFile->cd();
  // H1D *h1 = new TH1D("h2", "h2", 20, 0., 20.);

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpMuons::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpMuons);
