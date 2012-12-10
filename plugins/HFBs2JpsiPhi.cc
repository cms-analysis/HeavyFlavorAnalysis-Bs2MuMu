#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBs2JpsiPhi.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
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
HFBs2JpsiPhi::HFBs2JpsiPhi(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("generalTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)), 
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)), 
  fPhiWindow(iConfig.getUntrackedParameter<double>("phiWindow", 0.2)), 
  fBsWindow(iConfig.getUntrackedParameter<double>("BsWindow", 0.8)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fType(iConfig.getUntrackedParameter<int>("type", 531)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBs2JpsiPhi constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt :                  " << fMuonPt << endl;
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  phiWindow:                " << fPhiWindow << endl;
  cout << "---  BsWindow:                 " << fBsWindow << endl;
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
HFBs2JpsiPhi::~HFBs2JpsiPhi() {
  
}


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::analyze(const Event& iEvent, const EventSetup& iSetup)
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
		cout << "==>HFBs2JpsiPhi> No primary vertex collection found, skipping" << endl;
		return;
	}
	const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
	if (vertices.size() == 0) {
		cout << "==>HFBs2JpsiPhi> No primary vertex found, skipping" << endl;
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
		cout << "==>HFBs2JpsiPhi> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
		return;
	}
	const MuonCollection *muonColl = hMuons.product();

	// -- get the collection of tracks
	Handle<View<Track> > hTracks;
	iEvent.getByLabel(fTracksLabel, hTracks);
	if(!hTracks.isValid()) {
		cout << "==>HFBs2JpsiPhi> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
		return;
	}

	// -- Transient tracks for vertexing
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
	if (!fTTB.isValid()) {
		cout << " -->HFBs2JpsiPhi: Error: no TransientTrackBuilder found."<<endl;
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
	HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList),MMUON,
														  (fPsiMuons < 2 ? trkList : muonList),MMUON,
														 MJPSI-fPsiWindow,MJPSI+fPsiWindow, 1);
	HFTwoParticleCombinatoricsSet phiList = a.combine(trkList,MKAON,trkList,MKAON,MPHI-fPhiWindow,MPHI+fPhiWindow, 1);
	
	if (fVerbose > 0) cout << "==>HFBs2JpsiKp> J/psi list size: " << psiList.size() << endl;
	if (fVerbose > 0) cout << "==>HFBs2JpsiKp> phi list size: " << phiList.size() << endl;

	HFSequentialVertexFit aSeq(hTracks, fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);

	// -- Build J/psi + phi
	TLorentzVector psi, phi, m1, m2, ka1, ka2, bs;
	for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
		unsigned int iMuon1 = psiIt->first;
		unsigned int iMuon2 = psiIt->second;

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
		
		for (HFTwoParticleCombinatoricsNew::iterator phiIt = phiList.begin(); phiIt != phiList.end(); ++phiIt) {
			unsigned int iKaon1 = phiIt->first;
			unsigned int iKaon2 = phiIt->second;
			if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue; 
			if (iKaon2 == iMuon1 || iKaon2 == iMuon2) continue; 

			TrackBaseRef rTrackView1(hTracks, iKaon1);
			Track tKaon1(*rTrackView1);
			if (tKaon1.pt() < fTrackPt) continue;
			ka1.SetXYZM(tKaon1.px(), tKaon1.py(), tKaon1.pz(), MKAON); 
			if (psi.DeltaR(ka1) > fDeltaR) continue; 

			TrackBaseRef rTrackView2(hTracks, iKaon2);
			Track tKaon2(*rTrackView2);
			if (tKaon2.pt() < fTrackPt) continue;
			ka2.SetXYZM(tKaon2.px(), tKaon2.py(), tKaon2.pz(), MKAON); 
			if (psi.DeltaR(ka2) > fDeltaR) continue; 

			phi = ka1 + ka2; 
			if ((TMath::Abs(phi.M() - MPHI) > fPhiWindow)) continue;

			bs = psi + phi; 
			if (TMath::Abs(bs.M() - MBS) > fBsWindow) continue;
			
			// -- sequential fit: J/Psi kaons
			HFDecayTree theTree(300531, true, MBS, false, -1.0, true);

			HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false); // Don't use kinematic particle for the Psi
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			iterator = theTree.addDecayTree(300333, false, MPHI, false);
			iterator->addTrack(iKaon1,321);
			iterator->addTrack(iKaon2,321);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			aSeq.doFit(&theTree);

			// -- sequential fit: J/Psi (constraint) phi (unconstraint)
			theTree.clear(400531, true, MBS, false, -1.0, true);

			iterator = theTree.addDecayTree(400443, true, MJPSI, true);
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			iterator = theTree.addDecayTree(400333, false, MPHI, false);
			iterator->addTrack(iKaon1,321);
			iterator->addTrack(iKaon2,321);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			aSeq.doFit(&theTree);


			// -- global fit: J/Psi (constraint) phi (unconstraint)
			theTree.clear(500531, true, MBS, false, -1.0, true);
			iterator = theTree.addDecayTree(500443, false, MJPSI, false);
			iterator->addTrack(iMuon1,13,true);
			iterator->addTrack(iMuon2,13,true);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			iterator = theTree.addDecayTree(400333, false, MPHI, false);
			iterator->addTrack(iKaon1,321,false);
			iterator->addTrack(iKaon2,321,false);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			theTree.set_mass_tracks(MJPSI);
			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			aSeq.doFit(&theTree);
		}
	}
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBs2JpsiPhi::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBs2JpsiPhi::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBs2JpsiPhi);
