#include "HFBs2JpsiPhi.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatoricsNew.hh"

#include <iostream>

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFBs2JpsiPhi::HFBs2JpsiPhi(const edm::ParameterSet& iConfig) :
	HFVirtualDecay(iConfig),
	fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
	fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
	fPhiWindow(iConfig.getUntrackedParameter<double>("phiWindow", 0.2)),
	fBsWindow(iConfig.getUntrackedParameter<double>("BsWindow", 0.8))
{
	dumpConfiguration();
} // HFBs2JpsiPhi()

void HFBs2JpsiPhi::dumpConfiguration()
{
	using namespace std;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "--- HFBs2JpsiPhi configuration" << endl;
	HFVirtualDecay::dumpConfiguration();
	cout << "---  psiMuons:                 " << fPsiMuons << endl;
	cout << "---  psiWindow:                " << fPsiWindow << endl;
	cout << "---  phiWindow:                " << fPhiWindow << endl;
	cout << "---  BsWindow:                 " << fBsWindow << endl;
	cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace std;
	typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
	try {
		HFVirtualDecay::analyze(iEvent,iSetup);
	}
	catch (HFSetupException e) {
		cout << "==>HFBs2JpsiPhi> " << e.fMsg << endl;
		return;
	}
	
	fListBuilder->setMinPt(fMuonPt); // work with muon pt and not with track pt
	vector<int> muonList = fListBuilder->getMuonList();
	fListBuilder->setMinPt(fTrackPt);
	fListBuilder->setMaxDocaToTracks(fMaxDoca);
	fListBuilder->setCloseTracks(&muonList);
	vector<int> trkList = fListBuilder->getTrackList();
	
	if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
	HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
	HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList),MMUON,
													  (fPsiMuons < 2 ? trkList : muonList),MMUON,
													  MJPSI-fPsiWindow,MJPSI+fPsiWindow, 1);
	HFTwoParticleCombinatoricsSet phiList = a.combine(trkList,MKAON,trkList,MKAON,MPHI-fPhiWindow,MPHI+fPhiWindow, 1);
	
	if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> J/psi list size: " << psiList.size() << endl;
	if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> phi list size: " << phiList.size() << endl;
	
	// -- Build J/psi + phi
	TLorentzVector psi, phi, m1, m2, ka1, ka2, bs;
	for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
		unsigned int iMuon1 = psiIt->first;
		unsigned int iMuon2 = psiIt->second;
		
		reco::TrackBaseRef mu1TrackView(fTracksHandle, iMuon1);
		reco::Track tMuon1(*mu1TrackView);
		if (tMuon1.pt() < fMuonPt)  continue;
		m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);
		
		reco::TrackBaseRef mu2TrackView(fTracksHandle, iMuon2);
		reco::Track tMuon2(*mu2TrackView);
		if (tMuon2.pt() < fMuonPt)  continue;
		m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON);
		
		psi = m1 + m2;
		if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
		for (HFTwoParticleCombinatoricsNew::iterator phiIt = phiList.begin(); phiIt != phiList.end(); ++phiIt) {
			unsigned int iKaon1 = phiIt->first;
			unsigned int iKaon2 = phiIt->second;
			
			if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue;
			if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue;
			
			if (iKaon2 == iMuon1 || iKaon2 == iMuon2) continue;
			reco::TrackBaseRef rTrackView1(fTracksHandle, iKaon1);
			reco::Track tKaon1(*rTrackView1);
			if (tKaon1.pt() < fTrackPt) continue;
			ka1.SetXYZM(tKaon1.px(), tKaon1.py(), tKaon1.pz(), MKAON); 
			if (psi.DeltaR(ka1) > fDeltaR) continue;
			
			reco::TrackBaseRef rTrackView2(fTracksHandle, iKaon2);
			reco::Track tKaon2(*rTrackView2);
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
			
			fSequentialFitter->doFit(&theTree);
			
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
			
			fSequentialFitter->doFit(&theTree);
			
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
			
			fSequentialFitter->doFit(&theTree);
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBs2JpsiPhi);
