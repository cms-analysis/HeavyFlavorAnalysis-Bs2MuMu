#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpMuons.h"

#include <iostream>

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
  fCaloMuonsLabel(iConfig.getUntrackedParameter<InputTag>("calomuonsLabel")),
  fMaxTrackDistToStore(iConfig.getUntrackedParameter<double>("maxTrackDist",0.2)),
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
  cout << "---  fCaloMuonsLabel:         " << fCaloMuonsLabel.encode() << endl;
  cout << "---  fDoTruthMatching:        " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "---  fVerbose:                " << fVerbose << endl;
  cout << "---  fRunOnAOD:               " << fRunOnAOD << endl;
  cout << "---  fMaxTrackDistToStore:    " << fMaxTrackDistToStore << endl;
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

  TAnaMuon *pM = gHFEvent->addMuon();    
  pM->clear(); 
  if (rm.innerTrack().isNonnull()) {
    pM->fIndex = rm.innerTrack().index();
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
  //  pM->fTimeNdof         = rm.numberOfMatchedStations(); //FIXME!!!
  pM->fNmatchedStations = rm.numberOfMatchedStations();

  TrackRef gTrack = rm.globalTrack();
  TrackRef iTrack = rm.innerTrack();
  TrackRef oTrack = rm.outerTrack();


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
    //    cout << " gpt = " << trk.pt() << " " <<  trk.eta() << " "  <<  trk.phi() << endl;


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
  TrajectoryStateOnSurface prop_M1 = fpropM1.extrapolate(rm);
  TrajectoryStateOnSurface prop_M2 = fpropM2.extrapolate(rm);
  TVector3 muPosM1; 
  bool validM1(false); 
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
  
  // -- doca of close tracks to muon
  if (iTrack.isNonnull() && fTTB.isValid()) {
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


  if (validM1) {
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
