#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2JpsiKp.h"

#include <iostream>
#include <utility>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatoricsNew.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTrackListBuilder.hh"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;


// ----------------------------------------------------------------------
HFBu2JpsiKp::HFBu2JpsiKp(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("generalTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)), 
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)), 
  fBuWindow(iConfig.getUntrackedParameter<double>("BuWindow", 0.8)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)), 
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)), 
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fType(iConfig.getUntrackedParameter<int>("type", 521))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2JpsiKp constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt :                  " << fMuonPt << endl;
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  BuWindow:                 " << fBuWindow << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  maxD0:                    " << fMaxD0 << endl;
  cout << "---  maxDz:                    " << fMaxDz << endl;
  cout << "---  vertexing:                " << fVertexing << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFBu2JpsiKp::~HFBu2JpsiKp() {
  
}


// ----------------------------------------------------------------------
void HFBu2JpsiKp::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
	// -- get the magnetic field
	ESHandle<MagneticField> magfield;
	iSetup.get<IdealMagneticFieldRecord>().get(magfield);
	const MagneticField *field = magfield.product();

	// -- get the primary vertex
	Handle<VertexCollection> recoPrimaryVertexCollection;
	iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
	if(!recoPrimaryVertexCollection.isValid()) {
		cout << "==>HFBu2JpsiKp> No primary vertex collection found, skipping" << endl;
		return;
	}
	const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
	if (vertices.size() == 0) {
		cout << "==>HFBu2JpsiKp> No primary vertex found, skipping" << endl;
		return;
	}
	fPV = vertices[gHFEvent->fBestPV]; 
	if (fVerbose > 0) {
		cout << "HFDimuons: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
	}

	// -- get the collection of muons
	Handle<MuonCollection> hMuons;
	iEvent.getByLabel(fMuonsLabel, hMuons);
	if (!hMuons.isValid()) {
		cout << "==>HFBu2JpsiKp> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
		return;
	}
	const MuonCollection *muonColl = hMuons.product();

	// -- get the collection of tracks
	Handle<View<Track> > hTracks;
	iEvent.getByLabel(fTracksLabel, hTracks);
	if(!hTracks.isValid()) {
		cout << "==>HFBu2JpsiKp> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
		return;
	}

	// -- Transient track builder for vertexing
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
	if (!fTTB.isValid()) {
		cout << " -->HFBu2JpsiKp: Error: no TransientTrackBuilder found."<<endl;
		return;
	}
	
	HFTrackListBuilder listBuilder(hTracks,muonColl,fTTB.product(),fVerbose);
	listBuilder.setMaxD0(fMaxD0);
	listBuilder.setMaxDz(fMaxDz);
	listBuilder.setMinPt(fMuonPt);
	vector<int> muonList = listBuilder.getMuonList();
	listBuilder.setMinPt(fTrackPt);
	listBuilder.setMaxDocaToTracks(fMaxDoca);
	listBuilder.setCloseTracks(&muonList);
	vector<int> trkList = listBuilder.getTrackList();

	if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
	
	HFTwoParticleCombinatoricsNew a(hTracks,fVerbose);
	HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList), MMUON,
													  (fPsiMuons < 2 ? trkList : muonList), MMUON,
													  MJPSI - fPsiWindow, MJPSI + fPsiWindow, 1 );

	if (fVerbose > 0) cout << "==>HFBu2JpsiKp> J/psi list size: " << psiList.size() << endl;

	HFSequentialVertexFit aSeq(hTracks, hMuons.product(), fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);

	// -- Build J/psi + track
	TLorentzVector psi, cpsi, m1, m2, ka, bu;
	for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
		int iMuon1 = psiIt->first; 
		int iMuon2 = psiIt->second; 

		TrackBaseRef mu1TrackView(hTracks, iMuon1);
		Track tMuon1(*mu1TrackView);
		if (tMuon1.pt() < fMuonPt)  continue;
		m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 

		TrackBaseRef mu2TrackView(hTracks, iMuon2);
		Track tMuon2(*mu2TrackView);
		if (tMuon2.pt() < fMuonPt)  continue;
		m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON); 

		psi = m1 + m2; 
		if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
		for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
			if (*trkIt == iMuon1 || *trkIt == iMuon2) continue; 
			TrackBaseRef rTrackView(hTracks, *trkIt);
			Track tKaon(*rTrackView);
			if (tKaon.d0() > fMaxD0) continue;
			if (tKaon.dz() > fMaxDz) continue;
			if (tKaon.pt() < fTrackPt) continue;
			ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON); 
			if (psi.DeltaR(ka) > fDeltaR) continue; 

			bu = ka + psi; 
			if (TMath::Abs(bu.M() - MBPLUS) > fBuWindow) continue;

			// -- sequential fit: J/Psi kaon
			if (fVerbose > 5) cout << "==>HFBu2JpsiKp> going to sequential fit" << endl;
			HFDecayTree theTree(300521, true, MBPLUS, false, -1.0, true);

			HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false);
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			theTree.addTrack(*trkIt,321);
			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit without mass constraint" << endl;
			aSeq.doFit(&theTree);

			// -- sequential fit: J/Psi kaon
			theTree.clear(400521, true, MBPLUS, false, -1.0, true);

			iterator = theTree.addDecayTree(400443, true, MJPSI, true);
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit with mass constraint" << endl;

			theTree.addTrack(*trkIt,321);
			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			aSeq.doFit(&theTree);
			if (fVerbose > 5) cout << "==>HFBu2JpsiKp> done with fitting for track " << *trkIt << endl;

			// -- global fit: J/Psi kaon
			theTree.clear(500521, true, MBPLUS, false, -1.0, true);
			iterator = theTree.addDecayTree(500443, false, MJPSI, false);
			iterator->addTrack(iMuon1,13,true);
			iterator->addTrack(iMuon2,13,true);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			theTree.addTrack(*trkIt,321,false);
			theTree.set_mass_tracks(MJPSI);
			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			aSeq.doFit(&theTree);
		}
	}
}

// ------------ method called once each job just before starting event loop  ------------
void  HFBu2JpsiKp::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBu2JpsiKp::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2JpsiKp);
