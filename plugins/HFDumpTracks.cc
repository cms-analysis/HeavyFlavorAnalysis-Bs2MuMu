#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTracks.h"
#include "HFDumpMuons.h"

#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TSimpleTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpUtilities.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpTracks::HFDumpTracks(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("ctfWithMaterialTracks"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("primaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<InputTag>("beamSpotLabel", InputTag("offlineBeamSpot"))),
  fGenEventLabel(iConfig.getUntrackedParameter<InputTag>("generatorEventLabel", InputTag("source"))),
  fSimTracksLabel(iConfig.getUntrackedParameter<InputTag>("simTracksLabel", InputTag("famosSimHits"))),
  fAssociatorLabel(iConfig.getUntrackedParameter<InputTag>("associatorLabel", InputTag("TrackAssociatorByChi2"))), 
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<InputTag>("trackingParticlesLabel", InputTag("trackingParticles"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fDumpSimpleTracks(iConfig.getUntrackedParameter<bool>("dumpSimpleTracks", true)),
  fDumpRecTracks(iConfig.getUntrackedParameter<bool>("dumpRecTracks", false)),
  fPropMuon(iConfig.getParameter<edm::ParameterSet>("propMuon"))
    {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTracks constructor  " << endl;
  cout << "---  tracksLabel:             " << fTracksLabel << endl;
  cout << "---  primaryVertexLabel:      " << fPrimaryVertexLabel << endl;
  cout << "---  beamSpotLabel:           " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel:              " << fMuonsLabel << endl;
  cout << "---  generatorEventLabel:     " << fGenEventLabel << endl;
  cout << "---  simTracksLabel:          " << fSimTracksLabel << endl;
  cout << "---  associatorLabel:         " << fAssociatorLabel << endl;
  cout << "---  trackingParticlesLabel:  " << fTrackingParticlesLabel << endl;
  cout << "---  dumpSimpleTracks:        " << fDumpSimpleTracks << endl;
  cout << "---  dumpRecTracks:           " << fDumpRecTracks << endl;
  cout << "---  doTruthMatching:         " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS, 3 = TAna01Event
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpTracks::~HFDumpTracks() {
  
}


// ----------------------------------------------------------------------
void HFDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (fVerbose > 0) cout << "==>HFDumpTracks> new event " << endl;
  // -- get the collection of RecoTracks 
  edm::Handle<edm::View<reco::Track> > tracksView;
  iEvent.getByLabel(fTracksLabel, tracksView);

  // -- get the primary vertex
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDumpTracks> No primary vertex collection found, skipping" << endl;
    return;
  }
  const reco::VertexCollection *vc = recoPrimaryVertexCollection.product();
  if (vc->size() == 0) {
    cout << "==>HFDumpTracks> No primary vertex found, skipping" << endl;
    return;
  }



  // -- get the beam spot
  const reco::BeamSpot *bs;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(fBeamSpotLabel, beamSpotHandle);

  if (beamSpotHandle.isValid()) {
    bs = beamSpotHandle.product();
  } else {
    bs = 0; 
    cout << "==>HFDumpTracks> No beam spot available from EventSetup" << endl;
  }


  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices, muonCollectionIndices; 
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  const reco::MuonCollection *mc = (MuonCollection*)hMuons.product();

  int muonIndex(0); 
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    TrackRef track = muon->innerTrack();
    muonIndices.push_back(track.index());
    muonCollectionIndices.push_back(muonIndex); 
    ++muonIndex; 
  }
 
  if (fVerbose > 0) cout << "==>HFDumpTracks> nMuons = " << hMuons->size() << endl;
  TH1D *h2 = (TH1D*)gHFFile->Get("h2");
  if (h2) h2->Fill(hMuons->size());
 
  if (fVerbose > 0) cout << "===> Tracks " << tracksView->size() << endl;

  TH1D *h1 = (TH1D*)gHFFile->Get("h1");
  if (h1) h1->Fill(tracksView->size());

  // -- fill association map between tracks and primary vertices
  //  tracksAndPv(iEvent); 
  //   if (fVerbose > 20) {
  //     for (unsigned int i = 0; i < tracksView->size(); ++i) {
  //       TrackBaseRef rTrackView(tracksView,i);
  //       Track trackView(*rTrackView);
  //       cout << Form("%4d: %d", i, fTrack2Pv[i]) << " pt = " << trackView.pt() <<endl; 
  //     }
  //   }
  
  int genIdx(-1); 
  for (unsigned int i = 0; i < tracksView->size(); ++i){    

    TrackBaseRef rTrackView(tracksView,i);
    const Track trackView(*rTrackView);

    // -- Muon?
    int mid = 0; 
    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (i == muonIndices[im]) {
	mid = 1; 
	break;
      }
    }
    
    // -- truth matching with deltaR comparision
    genIdx = -1; 
    if (3 == fDoTruthMatching) {    
      genIdx  = gHFEvent->getGenIndexWithDeltaR(trackView.pt(), trackView.eta(), trackView.phi(), trackView.charge()); 
    }
    
    // -- fill the tracks
    if (fDumpSimpleTracks) {
      TSimpleTrack *st = gHFEvent->addSimpleTrack();
      fillSimpleTrack(st, trackView, i, mid, genIdx, vc); 
    } 

    if (fDumpRecTracks) {
      TAnaTrack *at = gHFEvent->addRecTrack();
      fillAnaTrack(at, trackView, i, genIdx, vc, mc, bs); 
    }
  }


  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  cleanupTruthMatching(tracksView, magfield);


  if (fVerbose > 0) {
    if (fDumpSimpleTracks) {
      for (int it = 0; it < gHFEvent->nSimpleTracks(); ++it) {
	gHFEvent->getSimpleTrack(it)->dump();
      }
    }

    if (fDumpRecTracks) {
      for (int it = 0; it < gHFEvent->nRecTracks(); ++it) {
	gHFEvent->getRecTrack(it)->dump();
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTracks::beginJob() {
  gHFFile->cd();
  //  TH1D *h1 = new TH1D("h2", "h2", 20, 0., 20.);

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTracks::endJob() {
}

void HFDumpTracks::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	fPropMuon.init(iSetup);
} // beginRun()



// // ----------------------------------------------------------------------
// void HFDumpTracks::tracksAndPv(const edm::Event& iEvent) {

//   for (int i = 0; i < HFMAXTRACK; ++i) {
//     fTrack2Pv[i] = -1; 
//   }

//   edm::Handle<edm::View<Track> > tracksView;
//   iEvent.getByLabel(fTracksLabel, tracksView);
  
//   edm::Handle<VertexCollection> recoPrimaryVertexCollection;
//   iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  
//   int cnt(0); 
//   for (VertexCollection::const_iterator iv = recoPrimaryVertexCollection->begin(); iv != recoPrimaryVertexCollection->end(); ++iv) {
//     Vertex::trackRef_iterator v1TrackIter;
//     Vertex::trackRef_iterator v1TrackBegin = iv->tracks_begin();
//     Vertex::trackRef_iterator v1TrackEnd   = iv->tracks_end();
//     for (v1TrackIter = v1TrackBegin; v1TrackIter != v1TrackEnd; v1TrackIter++) {
//       if (fVerbose > 100) {
// 	cout << "vtx " << cnt << " with trk: " << " " << v1TrackIter->key()
// 	     << "  " << (**v1TrackIter).pt() 
// 	     << "  " << (**v1TrackIter).phi();

// 	TrackBaseRef rTrackView(tracksView, v1TrackIter->key());
// 	Track track(*rTrackView);
// 	cout << " -> track " << "  " << track.pt() << "  " << track.phi();
// 	cout << endl;
//       }
//       if (v1TrackIter->key() < HFMAXTRACK) fTrack2Pv[v1TrackIter->key()] = cnt;
//     }
//     ++cnt;
//   }

// }





// // ----------------------------------------------------------------------
// void HFDumpTracks::fillAnaTrack(TAnaTrack *pTrack, const reco::Track &trackView, int tidx, int mid, int gidx,
// 				const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs) {
  
//   math::XYZPoint refPt; 
//   double xE(0.), yE(0), zE(0.);
//   int pvIdx = getPv(tidx, vc);
//   if (pvIdx > -1) {
//     refPt = vc->at(pvIdx).position();
//     xE = vc->at(pvIdx).xError();
//     yE = vc->at(pvIdx).yError();
//     zE = vc->at(pvIdx).zError();
//   } else {
//     if (bs) {
//       refPt = bs->position(); 
//       xE = bs->BeamWidthX();
//       yE = bs->BeamWidthY();
//       zE = bs->sigmaZ();
//     } else {
//       xE = 99.;
//       yE = 99.;
//       zE = 99.;
//     } 
//   }
  
//   pTrack->fIndex = tidx;
//   pTrack->fPlab.SetPtEtaPhi(trackView.pt(),
// 			    trackView.eta(),
// 			    trackView.phi()
// 			    );
//   pTrack->fPtE   = trackView.ptError();
//   pTrack->fPhiE  = trackView.phiError();
//   pTrack->fEtaE  = trackView.etaError();
  
//   pTrack->fTip   = trackView.dxy(refPt);
//   pTrack->fTipE  = sqrt(trackView.d0Error()*trackView.d0Error() + 0.5*xE*xE + 0.5*yE*yE); 
  
//   pTrack->fLip   = trackView.dz(refPt); 
//   pTrack->fLipE  = sqrt(trackView.dzError()*trackView.dzError() + zE*zE); 
//   pTrack->fPvIdx = pvIdx;
  

//   math::XYZPoint bsPt(bs->x0(), bs->y0(), bs->z0());
//   pTrack->fBsTip  = trackView.dxy(bsPt);
//   pTrack->fBsTipE = sqrt(trackView.d0Error()*trackView.d0Error() 
//  			 + 0.5*bs->BeamWidthX()*bs->BeamWidthX() 
//  			 + 0.5*bs->BeamWidthY()*bs->BeamWidthY()); 
//   pTrack->fBsLip  = trackView.dz(bsPt);
//   pTrack->fBsLipE = sqrt(trackView.dzError()*trackView.dzError() + bs->sigmaZ()*bs->sigmaZ());
  
//   pTrack->fdxy  = trackView.dxy();
//   pTrack->fdxyE = trackView.dxyError();
//   pTrack->fd0   = trackView.d0();
//   pTrack->fd0E  = trackView.d0Error();
//   pTrack->fdz   = trackView.dz();
//   pTrack->fdzE  = trackView.dzError();
//   pTrack->fdsz  = trackView.dsz();
//   pTrack->fdszE = trackView.dszError();
  
//   pTrack->fQ = trackView.charge();
//   pTrack->fChi2 = trackView.chi2();
//   pTrack->fDof = int(trackView.ndof());
//   pTrack->fValidHits = trackView.numberOfValidHits();  
//   pTrack->fAlgorithm = trackView.algo(); 
  
//   // -- see https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=123067
//   pTrack->fValidHitFraction = trackView.validFraction();
  
//   //fixme move into HFDumpMuons and store with TAnaMuon
//   //   // store entrance to muon chamber
//   //   TrajectoryStateOnSurface tsos = fPropMuon.extrapolate(trackView);
//   //   if (tsos.isValid()) {
//   //     pTrack->fPosMuonChamber.SetXYZ(tsos.globalPosition().x(),tsos.globalPosition().y(),tsos.globalPosition().z());
//   //     pTrack->fPlabMuonChamber.SetXYZ(tsos.globalMomentum().x(),tsos.globalMomentum().y(),tsos.globalMomentum().z());
//   //   }
  
//   // -- from: RecoBTag/TrackProbability/src/TrackClassFilter.cc
//   reco::TrackBase::TrackQuality trackQualityUndef               =  reco::TrackBase::qualityByName("undefQuality");
//   reco::TrackBase::TrackQuality trackQualityLoose               =  reco::TrackBase::qualityByName("loose");
//   reco::TrackBase::TrackQuality trackQualityTight               =  reco::TrackBase::qualityByName("tight");
//   reco::TrackBase::TrackQuality trackQualityhighPur             =  reco::TrackBase::qualityByName("highPurity");
//   reco::TrackBase::TrackQuality trackQualityConfirmed           =  reco::TrackBase::qualityByName("confirmed");
//   reco::TrackBase::TrackQuality trackQualityGoodIterative       =  reco::TrackBase::qualityByName("goodIterative");
//   reco::TrackBase::TrackQuality trackQualityLooseSetWithPV      =  reco::TrackBase::qualityByName("looseSetWithPV");
//   reco::TrackBase::TrackQuality trackQualityHighPuritySetWithPV =  reco::TrackBase::qualityByName("highPuritySetWithPV");
  
//   int trakQuality  = 0;
//   if (trackView.quality(trackQualityUndef))               trakQuality |= 0x1<<10;
//   if (trackView.quality(trackQualityLoose))               trakQuality |= 0x1<<0;
//   if (trackView.quality(trackQualityTight))               trakQuality |= 0x1<<1;
//   if (trackView.quality(trackQualityhighPur))             trakQuality |= 0x1<<2;
//   if (trackView.quality(trackQualityConfirmed))           trakQuality |= 0x1<<3;
//   if (trackView.quality(trackQualityGoodIterative))       trakQuality |= 0x1<<4;
//   if (trackView.quality(trackQualityLooseSetWithPV))      trakQuality |= 0x1<<5;
//   if (trackView.quality(trackQualityHighPuritySetWithPV)) trakQuality |= 0x1<<6;
//   pTrack->fTrackQuality = trakQuality; 
  
//   // -- Muon ID
//   pTrack->fMuIndex = -1; 
//   pTrack->fMuID    = 0; 
//   for (unsigned int i = 0; i < mc->size(); ++i) {
//     const reco::Muon muon = mc->at(i); 
//     TrackRef track = muon.innerTrack();
//     if (static_cast<unsigned int>(tidx) == track.index()) {
//       pTrack->fMuIndex = i; 
//       pTrack->fMuID    = HFDumpMuons::muonID(muon);
//       break;
//     }
//   }

  
//   // -- hits of the track
//   const reco::HitPattern& p = trackView.hitPattern();
//   for (int i=0; i<p.numberOfHits(); i++) {
//     uint32_t hit = p.getHitPattern(i);
//     if (i < 20) pTrack->fHitPattern[i] = hit; 
//   }
  
//   //fixme add truth matching

// }




//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTracks);
