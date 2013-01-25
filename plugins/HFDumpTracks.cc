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

#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

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
  fCaloMuonsLabel(iConfig.getUntrackedParameter<InputTag>("calomuonsLabel")),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fLoadCalomuons(iConfig.getUntrackedParameter<bool>("loadCalomuons",true)),
  fPropMuon(iConfig.getParameter<edm::ParameterSet>("propMuon"))
    {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTracks constructor  " << endl;
  cout << "---  tracksLabel:             " << fTracksLabel << endl;
  cout << "---  primaryVertexLabel:      " << fPrimaryVertexLabel << endl;
  cout << "---  beamSpotLabel:           " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel:              " << fMuonsLabel << endl;
  cout << "---  calomuonsLabel:          " << fCaloMuonsLabel << endl;
  cout << "---  generatorEventLabel:     " << fGenEventLabel << endl;
  cout << "---  simTracksLabel:          " << fSimTracksLabel << endl;
  cout << "---  associatorLabel:         " << fAssociatorLabel << endl;
  cout << "---  trackingParticlesLabel:  " << fTrackingParticlesLabel << endl;
  cout << "---  doTruthMatching:         " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS, 3 = TAna01Event
  cout << "---  loadCalomuons:           " << fLoadCalomuons << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpTracks::~HFDumpTracks() {
  
}


// ----------------------------------------------------------------------
void HFDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  static int first(1); 
  if (1 == first) {
    first = 0; 
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    iSetup.get<TrackAssociatorRecord>().get(fAssociatorLabel.encode(), theAssociator);
    fAssociator = theAssociator.product();
    cout << "==>HFDumpTracks> fAssociator = " << fAssociator << endl;
  }


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
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFDumpTracks> No primary vertex found, skipping" << endl;
    return;
  }

  // -- get the beam spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(fBeamSpotLabel, beamSpotHandle);

  if (beamSpotHandle.isValid()) {
    beamSpot = *beamSpotHandle;
  } else {
    cout << "==>HFDumpTracks> No beam spot available from EventSetup" << endl;
  }




  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices, muonCollectionIndices;
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  int muonIndex(0); 
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    TrackRef track = muon->innerTrack();
    muonIndices.push_back(track.index());
    muonCollectionIndices.push_back(muonIndex); 
    ++muonIndex; 
  }
  
  // -- get the collection of calo muons and store their corresponding track indices
  vector<unsigned int> CalomuonIndices, CalomuonCollectionIndices;  
  if (fLoadCalomuons) {
	  Handle<CaloMuonCollection> cMuons;
	  iEvent.getByLabel(fCaloMuonsLabel, cMuons);
	  for (CaloMuonCollection::const_iterator cmuon = cMuons->begin(); cmuon != cMuons->end(); ++cmuon) {
		TrackRef track = cmuon->innerTrack();
		CalomuonIndices.push_back(track.index());
		CalomuonCollectionIndices.push_back(muonIndex); 
		++muonIndex; 
	  }
  }
  
  if (fVerbose > 0) cout << "==>HFDumpTracks> nMuons = " << hMuons->size() << endl;
  TH1D *h2 = (TH1D*)gHFFile->Get("h2");
  if (h2) h2->Fill(hMuons->size());
 
  // -- get the tracking particle collection needed for truth matching. Only on RECO data tier?
  RecoToSimCollection recSimColl;
  const RecoToSimCollection recSimColl2;

  if (1 == fDoTruthMatching) {
    try {
      edm::Handle<TrackingParticleCollection> trackingParticles;
      iEvent.getByLabel(fTrackingParticlesLabel, trackingParticles);
      recSimColl = fAssociator->associateRecoToSim(tracksView, trackingParticles, &iEvent);
    } catch (cms::Exception &ex) {
      cout << ex.explainSelf() << endl;
    }
  }


  // -- Get the stuff needed for FAMOS truth matching
  Handle<HepMCProduct> hepmc;
  const HepMC::GenEvent *genEvent = 0;
  edm::Handle<std::vector<SimTrack> > simTracks;
  if (2 == fDoTruthMatching) {
    iEvent.getByLabel(fGenEventLabel, hepmc);
    genEvent = hepmc->GetEvent();
    iEvent.getByLabel(fSimTracksLabel, simTracks); 
  }      


  // -- dE/dx (from Loic)
  Handle<ValueMap<DeDxData> > dEdxTrackHandle;
  try { iEvent.getByLabel("dedxHarmonic2", dEdxTrackHandle); } catch (...) {;}
  const ValueMap<DeDxData> dEdxTrack = *dEdxTrackHandle.product();
  

  if (fVerbose > 0) cout << "===> Tracks " << tracksView->size() << endl;
  TAnaTrack *pTrack; 
  TH1D *h1 = (TH1D*)gHFFile->Get("h1");
  if (h1) h1->Fill(tracksView->size());

  // -- fill association map between tracks and primary vertices
  tracksAndPv(iEvent); 
  if (fVerbose > 20) {
    for (unsigned int i = 0; i < tracksView->size(); ++i) {
      TrackBaseRef rTrackView(tracksView,i);
      Track trackView(*rTrackView);
      cout << Form("%4d: %d", i, fTrack2Pv[i]) << " pt = " << trackView.pt() <<endl; 
    }
  }
  
  for (unsigned int i = 0; i < tracksView->size(); ++i){    

    TrackBaseRef rTrackView(tracksView,i);
    Track trackView(*rTrackView);

    math::XYZPoint refPt; 
    double xE(0.), yE(0), zE(0.);
    if (fTrack2Pv[i] > -1) {
      //      fPV   = vertices[fTrack2Pv[i]]; 
      refPt = vertices[fTrack2Pv[i]].position();
      xE = vertices[fTrack2Pv[i]].xError();
      yE = vertices[fTrack2Pv[i]].yError();
      zE = vertices[fTrack2Pv[i]].zError();
    } else {
      // FIXME!!! Change to beamspot
      refPt = beamSpot.position(); 
      xE = beamSpot.BeamWidthX();
      yE = beamSpot.BeamWidthY();
      zE = beamSpot.sigmaZ();
    }

    pTrack = gHFEvent->addRecTrack();
    pTrack->fIndex = i;
    pTrack->fPlab.SetPtEtaPhi(trackView.pt(),
			      trackView.eta(),
			      trackView.phi()
			      );
    pTrack->fPtE   = trackView.ptError();
    pTrack->fPhiE  = trackView.phiError();
    pTrack->fEtaE  = trackView.etaError();

    pTrack->fTip   = trackView.dxy(refPt);
    pTrack->fTipE  = sqrt(trackView.d0Error()*trackView.d0Error() + 0.5*xE*xE + 0.5*yE*yE); 
    
    pTrack->fLip   = trackView.dz(refPt); 
    pTrack->fLipE  = sqrt(trackView.dzError()*trackView.dzError() + zE*zE); 
    pTrack->fPvIdx = fTrack2Pv[i];

    math::XYZPoint bsPt(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
    pTrack->fBsTip  = trackView.dxy(bsPt);
    pTrack->fBsTipE = sqrt(trackView.d0Error()*trackView.d0Error() 
			   + 0.5*beamSpot.BeamWidthX()*beamSpot.BeamWidthX() 
			   + 0.5*beamSpot.BeamWidthY()*beamSpot.BeamWidthY()); 
    pTrack->fBsLip  = trackView.dz(bsPt);
    pTrack->fBsLipE = sqrt(trackView.dzError()*trackView.dzError() + beamSpot.sigmaZ()*beamSpot.sigmaZ());

    pTrack->fdxy  = trackView.dxy();
    pTrack->fdxyE = trackView.dxyError();
    pTrack->fd0   = trackView.d0();
    pTrack->fd0E  = trackView.d0Error();
    pTrack->fdz   = trackView.dz();
    pTrack->fdzE  = trackView.dzError();
    pTrack->fdsz  = trackView.dsz();
    pTrack->fdszE = trackView.dszError();

    pTrack->fQ = trackView.charge();
    pTrack->fChi2 = trackView.chi2();
    pTrack->fDof = int(trackView.ndof());
    pTrack->fValidHits = trackView.numberOfValidHits();  
    pTrack->fAlgorithm = trackView.algo(); 

    // -- see https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=123067
    pTrack->fValidHitFraction = trackView.validFraction();
	
	// store entrance to muon chamber
	TrajectoryStateOnSurface tsos = fPropMuon.extrapolate(trackView);
	if (tsos.isValid()) {
		pTrack->fPosMuonChamber.SetXYZ(tsos.globalPosition().x(),tsos.globalPosition().y(),tsos.globalPosition().z());
		pTrack->fPlabMuonChamber.SetXYZ(tsos.globalMomentum().x(),tsos.globalMomentum().y(),tsos.globalMomentum().z());
	}
	
    // -- from: RecoBTag/TrackProbability/src/TrackClassFilter.cc
    reco::TrackBase::TrackQuality trackQualityUndef               =  reco::TrackBase::qualityByName("undefQuality");
    reco::TrackBase::TrackQuality trackQualityLoose               =  reco::TrackBase::qualityByName("loose");
    reco::TrackBase::TrackQuality trackQualityTight               =  reco::TrackBase::qualityByName("tight");
    reco::TrackBase::TrackQuality trackQualityhighPur             =  reco::TrackBase::qualityByName("highPurity");
    reco::TrackBase::TrackQuality trackQualityConfirmed           =  reco::TrackBase::qualityByName("confirmed");
    reco::TrackBase::TrackQuality trackQualityGoodIterative       =  reco::TrackBase::qualityByName("goodIterative");
    reco::TrackBase::TrackQuality trackQualityLooseSetWithPV      =  reco::TrackBase::qualityByName("looseSetWithPV");
    reco::TrackBase::TrackQuality trackQualityHighPuritySetWithPV =  reco::TrackBase::qualityByName("highPuritySetWithPV");

    int trakQuality  = 0;
    if (trackView.quality(trackQualityUndef))               trakQuality |= 0x1<<10;
    if (trackView.quality(trackQualityLoose))               trakQuality |= 0x1<<0;
    if (trackView.quality(trackQualityTight))               trakQuality |= 0x1<<1;
    if (trackView.quality(trackQualityhighPur))             trakQuality |= 0x1<<2;
    if (trackView.quality(trackQualityConfirmed))           trakQuality |= 0x1<<3;
    if (trackView.quality(trackQualityGoodIterative))       trakQuality |= 0x1<<4;
    if (trackView.quality(trackQualityLooseSetWithPV))      trakQuality |= 0x1<<5;
    if (trackView.quality(trackQualityHighPuritySetWithPV)) trakQuality |= 0x1<<6;
    pTrack->fTrackQuality = trakQuality; 

    // -- filling dE/dx information (from Loic)
    pTrack->fDeDx      = dEdxTrack[rTrackView].dEdx();
    pTrack->fDeDxNmeas = dEdxTrack[rTrackView].numberOfSaturatedMeasurements();
    pTrack->fDeDxNsat  = dEdxTrack[rTrackView].numberOfMeasurements();

    // -- Muon ID
    pTrack->fMuIndex = -1; 
    pTrack->fMuID    = 0;
    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (i == muonIndices[im]) {
		if (fVerbose > 3) cout << " ==>HFDumpTracks> Found a muon!!" << endl;
		const Muon &rMuon  = hMuons->at(muonCollectionIndices[im]);

		pTrack->fMuIndex = muonCollectionIndices[im]; 
		pTrack->fMuID    = HFDumpMuons::muonID(rMuon);
      }
    }
    
    for (unsigned int im = 0; im < CalomuonIndices.size(); ++im) {
      if (i == CalomuonIndices[im]) {
		if (fVerbose > 3) cout << " ==>HFDumpTracks> Found a calo muon!!" << endl;
		
		if (pTrack->fMuIndex < 0)
			pTrack->fMuIndex = CalomuonCollectionIndices[im];
		pTrack->fMuID   |= 0x1<<15;
      }
    }


    // -- hits of the track
    const reco::HitPattern& p = trackView.hitPattern();
    for (int i=0; i<p.numberOfHits(); i++) {
      uint32_t hit = p.getHitPattern(i);
      if (i < 20) pTrack->fHitPattern[i] = hit; 
    }

    int gen_pdg_id(-99999), gen_id(-99999), gen_cnt(0);

    // -- simple truth matching
    //    this procedure can be redone offline on the ntuple level!
    if (0 == fDoTruthMatching) {
      gen_id     = gHFEvent->getGenIndex(trackView.px(), trackView.py(), trackView.pz(), trackView.charge()); 
      if (gen_id > -1) gen_pdg_id = gHFEvent->getGenCand(gen_id)->fID; 
      if (13 == TMath::Abs(gen_pdg_id)) {
	if (fVerbose > 4) {cout << "Simple TM: "; pTrack->dump(); }
      }
    }

    // -- RECO truth matching with TrackingParticle
    if (1 == fDoTruthMatching) {
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[rTrackView];
 	TrackingParticleRef tpr = tp.begin()->first;  
 	const HepMC::GenParticle *genPar = 0; 
 	for (TrackingParticle::genp_iterator pit = tpr->genParticle_begin(); 
 	     pit != tpr->genParticle_end(); 
 	     ++pit){
 	  genPar     = pit->get();
 	  gen_pdg_id = (*genPar).pdg_id();
 	  gen_id     = (*genPar).barcode()-1;
 	  //gen_pt     = (*genPar).momentum().perp();
 	  //gen_phi    = (*genPar).momentum().phi();
 	  //gen_eta    = (*genPar).momentum().pseudoRapidity();
 	  gen_cnt++;
	}
      } catch (Exception event) {
	// 	cout << "%%>   Rec. Track #" << setw(2) << rTrack.index() << " pT: " 
	// 	     << setprecision(2) << setw(6) << trackView.pt() 
	// 	     <<  " matched to 0 MC Tracks" << endl;
      }
    }

    // -- FAMOS truth matching via SimHit
    if (2 == fDoTruthMatching) {
      for (trackingRecHit_iterator it = trackView.recHitsBegin();  it != trackView.recHitsEnd(); it++) {
	if ((*it)->isValid()) {
	  
	  int currentId(-1);
	  if (const SiTrackerGSRecHit2D *rechit = dynamic_cast<const SiTrackerGSRecHit2D *> (it->get())) {
	    currentId = rechit->simtrackId();          
	  }
	  
	  for (SimTrackContainer::const_iterator simTrack = simTracks->begin(); 
	       simTrack != simTracks->end(); 
	       simTrack++)   { 
	    
	    if (static_cast<int>(simTrack->trackId()) == currentId) {
	      int igen = simTrack->genpartIndex();
	      HepMC::GenParticle *genPar = genEvent->barcode_to_particle(igen);
	      if (genPar) {
		gen_pdg_id = (*genPar).pdg_id();
		gen_id     = (*genPar).barcode()-1;  // BC(HepMC) = BC(reco)+1
		if (fVerbose > 4) cout << Form("id = %i, bc = %i", gen_pdg_id, gen_id) << endl;
		goto done;
	      }
	    }
	  }
	}
      }
    done:;
    }            
    
    // -- truth matching with deltaR comparision
    if (3 == fDoTruthMatching) {    
      gen_id     = gHFEvent->getGenIndexWithDeltaR(trackView.pt(), trackView.eta(), trackView.phi(), trackView.charge()); 
      if (gen_id > -1) gen_pdg_id = gHFEvent->getGenCand(gen_id)->fID; 
      if (13 == TMath::Abs(gen_pdg_id)) {
	if (fVerbose > 4) {cout << "Simple TM: "; pTrack->dump(); }
      }
    }


    pTrack->fGenIndex = gen_id;
    pTrack->fMCID     = gen_pdg_id;
   }

  if (fVerbose > 0) {
    for (int it = 0; it < gHFEvent->nRecTracks(); ++it) {
      gHFEvent->getRecTrack(it)->dump();
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

// ----------------------------------------------------------------------
void HFDumpTracks::tracksAndPv(const edm::Event& iEvent) {

  for (int i = 0; i < HFMAXTRACK; ++i) {
    fTrack2Pv[i] = -1; 
  }

  edm::Handle<edm::View<Track> > tracksView;
  iEvent.getByLabel(fTracksLabel, tracksView);
  
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  
  int cnt(0); 
  for (VertexCollection::const_iterator iv = recoPrimaryVertexCollection->begin(); iv != recoPrimaryVertexCollection->end(); ++iv) {
    Vertex::trackRef_iterator v1TrackIter;
    Vertex::trackRef_iterator v1TrackBegin = iv->tracks_begin();
    Vertex::trackRef_iterator v1TrackEnd   = iv->tracks_end();
    for (v1TrackIter = v1TrackBegin; v1TrackIter != v1TrackEnd; v1TrackIter++) {
      if (fVerbose > 100) {
	cout << "vtx " << cnt << " with trk: " << " " << v1TrackIter->key()
	     << "  " << (**v1TrackIter).pt() 
	     << "  " << (**v1TrackIter).phi();

	TrackBaseRef rTrackView(tracksView, v1TrackIter->key());
	Track track(*rTrackView);
	cout << " -> track " << "  " << track.pt() << "  " << track.phi();
	cout << endl;
      }
      if (v1TrackIter->key() < HFMAXTRACK) fTrack2Pv[v1TrackIter->key()] = cnt;
    }
    ++cnt;
  }

}


//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTracks);
