#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDiTracks.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatoricsNew.hh"

using std::cout;
using std::endl;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFDiTracks::HFDiTracks(const edm::ParameterSet& iConfig) :
	HFVirtualDecay(iConfig),
	fTrack1Mass(iConfig.getUntrackedParameter<double>("track1Mass", MMUON)),
	fTrack2Mass(iConfig.getUntrackedParameter<double>("track2Mass", MMUON)),
	fMassLow(iConfig.getUntrackedParameter<double>("massLow", 0.0)),
	fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 12.0)),
	fNbrMuons(iConfig.getUntrackedParameter<int>("nbrMuons",2)),
	fCloseToMuons(iConfig.getUntrackedParameter<bool>("closeToMuons",false))
{
	dumpConfiguration();
} // HFDiTracks()

void HFDiTracks::dumpConfiguration()
{
	cout << "----------------------------------------------------------------------" << endl;
	cout << "--- HFDiTracks constructor" << endl;
	HFVirtualDecay::dumpConfiguration();
	cout << "---  track1Mass:               " << fTrack1Mass << endl;
	cout << "---  track2Mass:               " << fTrack2Mass << endl;
	cout << "---  massLow:                  " << fMassLow << endl;
	cout << "---  massHigh:                 " << fMassHigh << endl;
	cout << "---  nbrMuons:                 " << fNbrMuons << endl;
	cout << "---  closeToMuons:             " << fCloseToMuons << endl;
	cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()

void HFDiTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
	try {
		HFVirtualDecay::analyze(iEvent,iSetup);
	}
	catch (HFSetupException e) {
		cout << "==>HFDiTracks> " << e.fMsg << endl;
		return;
	}
	
	std::vector<int> trkList;
	std::vector<int> muonList;
	
	if (fNbrMuons > 0 || fCloseToMuons) {
		fListBuilder->setMinPt(fMuonPt);
		muonList = fListBuilder->getMuonList();
	}
	if (fNbrMuons < 2) {
		fListBuilder->setMinPt(fTrackPt);
		if (fCloseToMuons) {
			fListBuilder->setMaxDocaToTracks(fMaxDoca);
			fListBuilder->setCloseTracks(&muonList);
		}
		trkList = fListBuilder->getTrackList();
	}
	
	HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
	HFTwoParticleCombinatoricsSet candSet = a.combine( (fNbrMuons < 1 ? trkList : muonList), fTrack1Mass,
													  (fNbrMuons < 2 ? trkList : muonList), fTrack2Mass,
													  fMassLow, fMassHigh,
													  TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
	
	if (fVerbose > 0) cout << "==>HFDiTracks> candidate list size: " << candSet.size() << endl;
	
	for (HFTwoParticleCombinatoricsNew::iterator trkIt = candSet.begin(); trkIt != candSet.end(); ++trkIt) {
		
		// -- Vertexing, with Kinematic Particles
		HFDecayTree theTree(fType, true, 0, false);
		theTree.addTrack(trkIt->first, idFromMass(fTrack1Mass));
		theTree.addTrack(trkIt->second, idFromMass(fTrack2Mass));
		theTree.setNodeCut(RefCountedHFNodeCut(new HFPvWeightCut(fMaxDoca,fPvWeight)));
		
		fSequentialFitter->doFit(&theTree);
	}
}

int HFDiTracks::idFromMass(double mass) {
  if (TMath::Abs(mass - MELECTRON) < 0.0001) return 11;
  if (TMath::Abs(mass - MMUON) < 0.0001) return 13;
  if (TMath::Abs(mass - MPION) < 0.0001) return 211;
  if (TMath::Abs(mass - MKAON) < 0.0001) return 321;
  if (TMath::Abs(mass - MPROTON) < 0.0001) return 2212;
  cout << "%%%> HFDiTracks: mass " << mass << " not associated to any stable particle" << endl;
  return 0; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDiTracks);
