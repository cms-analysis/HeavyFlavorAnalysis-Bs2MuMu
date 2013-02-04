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
	HFVirtualDecay(iConfig),
	fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
	fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
	fBuWindow(iConfig.getUntrackedParameter<double>("BuWindow", 0.8))
{
	dumpConfiguration();
}

void HFBu2JpsiKp::dumpConfiguration()
{
	using namespace std;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "--- HFBu2JpsiKp configuration" << endl;
	HFVirtualDecay::dumpConfiguration();
	cout << "---  psiMuons:                 " << fPsiMuons << endl;
	cout << "---  psiWindow:                " << fPsiWindow << endl;
	cout << "---  BuWindow:                " << fBuWindow << endl;
	cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()

// ----------------------------------------------------------------------
void HFBu2JpsiKp::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
	try {
	  HFVirtualDecay::analyze(iEvent,iSetup);
	} catch(HFSetupException e) {
	  cout << "==>HFBu2JpsiKp> " << e.fMsg << endl;
	  return;
	}

	fListBuilder->setMinPt(fMuonPt); // work with muon pt
	vector<int> muonList = fListBuilder->getMuonList();
	if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
	
	fListBuilder->setMinPt(fTrackPt);
	fListBuilder->setMaxDocaToTracks(fMaxDoca);
	fListBuilder->setCloseTracks(&muonList);
	vector<int> trkList = fListBuilder->getTrackList();
	
	HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
	HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList), MMUON,
													  (fPsiMuons < 2 ? trkList : muonList), MMUON,
													  MJPSI - fPsiWindow, MJPSI + fPsiWindow, 1 );

	if (fVerbose > 0) cout << "==>HFBu2JpsiKp> J/psi list size: " << psiList.size() << endl;

	// -- Build J/psi + track
	TLorentzVector psi, cpsi, m1, m2, ka, bu;
	for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
		int iMuon1 = psiIt->first; 
		int iMuon2 = psiIt->second; 

		reco::TrackBaseRef mu1TrackView(fTracksHandle, iMuon1);
		reco::Track tMuon1(*mu1TrackView);
		if (tMuon1.pt() < fMuonPt)  continue;
		m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 

		TrackBaseRef mu2TrackView(fTracksHandle, iMuon2);
		Track tMuon2(*mu2TrackView);
		if (tMuon2.pt() < fMuonPt)  continue;
		m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON); 

		psi = m1 + m2; 
		if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
		for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
			if (*trkIt == iMuon1 || *trkIt == iMuon2) continue; 
			TrackBaseRef rTrackView(fTracksHandle, *trkIt);
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
			fSequentialFitter->doFit(&theTree);

			// -- sequential fit: J/Psi kaon
			theTree.clear(400521, true, MBPLUS, false, -1.0, true);

			iterator = theTree.addDecayTree(400443, true, MJPSI, true);
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit with mass constraint" << endl;

			theTree.addTrack(*trkIt,321);
			theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

			fSequentialFitter->doFit(&theTree);
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

			fSequentialFitter->doFit(&theTree);
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2JpsiKp);
