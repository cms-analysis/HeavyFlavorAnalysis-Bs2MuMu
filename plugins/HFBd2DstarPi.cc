#include "HFBd2DstarPi.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatoricsNew.hh"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <vector>

// -- Yikes!
extern TAna01Event *gHFEvent;

using std::cout;
using std::endl;

// ----------------------------------------------------------------------
HFBd2DstarPi::HFBd2DstarPi(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fSlowPionPt(iConfig.getUntrackedParameter<double>("slowPionPt", 0.4)),
  fD0Window(iConfig.getUntrackedParameter<double>("D0Window", 0.1)),
  fDeltaM(iConfig.getUntrackedParameter<double>("deltaM", 0.03))
{
  dumpConfiguration();
}

void HFBd2DstarPi::dumpConfiguration()
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBd2DstarPi configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  D0Window:                 " << fD0Window << endl;
  cout << "---  DeltaM:                   " << fDeltaM << endl;
  cout << "---  slowPionPt:               " << fSlowPionPt << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()

// ----------------------------------------------------------------------
void HFBd2DstarPi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch (HFSetupException e) {
    cout << "==>HFBd2DstarPi> " << e.fMsg << endl;
    return; // problem with setup
  }
	
  // -- Build lists
  std::vector<int> trkHighList = fListBuilder->getTrackList(); // to be used for pion, kaon of D0 and fast pion of B0, default d0,dz,pt settings
  fListBuilder->setMaxD0(999.0);
  fListBuilder->setMaxD0(999.0);
  fListBuilder->setMaxDz(999.0);
  fListBuilder->setMinPt(fSlowPionPt); // reduce minimum pt requirement for slow pion
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&trkHighList);
  std::vector<int> trkLowList = fListBuilder->getTrackList();
	
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet d0Set = a.combine(trkHighList, MKAON, trkHighList, MPION, MD0-fD0Window, MD0+fD0Window, 0);
  if (fVerbose > 0) cout << "==>HFBd2DstarPi> d0 set size: " << d0Set.size() << endl;
	
  if (d0Set.size() == 0)
    return;
	
  HFTwoParticleCombinatoricsSet piSet = a.combine(trkLowList, MPION, trkHighList, MPION, 0.0, 6.0, 0);
  if (fVerbose > 0) cout << "==>HFBd2DstarPi> pion set size: " << piSet.size() << endl;
	
  TLorentzVector ka,pi,d0,pis,pif,dstar,b0;
  for (HFTwoParticleCombinatoricsNew::iterator d0It = d0Set.begin(); d0It != d0Set.end(); ++d0It) {
    unsigned int iKaon = d0It->first;
    unsigned int iPion = d0It->second;
		
    reco::TrackBaseRef kaonTrackView(fTracksHandle, iKaon);
    reco::Track tKaon(*kaonTrackView);
		
    ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON);
		
    reco::TrackBaseRef pionTrackView(fTracksHandle, iPion);
    reco::Track tPion(*pionTrackView);
    if (tPion.pt() < fTrackPt)  continue;
    pi.SetPtEtaPhiM(tPion.pt(), tPion.eta(), tPion.phi(), MPION);
		
    d0 = ka + pi;
    for (HFTwoParticleCombinatoricsNew::iterator piIt = piSet.begin(); piIt != piSet.end(); ++piIt) {
		
      unsigned int iPionSlow = piIt->first;
      unsigned int iPionFast = piIt->second;
			
      if (iPionSlow == iKaon || iPionFast == iKaon) continue;
      if (iPionSlow == iPion || iPionFast == iPion) continue;
			
      reco::TrackBaseRef rTrackView1(fTracksHandle, iPionSlow);
      reco::Track tPionSlow(*rTrackView1);
      pis.SetXYZM(tPionSlow.px(), tPionSlow.py(), tPionSlow.pz(), MPION);
			
      reco::TrackBaseRef rTrackView2(fTracksHandle, iPionFast);
      reco::Track tPionFast(*rTrackView2);
      pif.SetXYZM(tPionFast.px(), tPionFast.py(), tPionFast.pz(), MPION);
			
      // -- the slow D* pion has the same charge like the fast D0 pion, and other charge correlations
      if (tPionSlow.charge()*tPion.charge() < 0) continue;
      if (tPionSlow.charge()*tPionFast.charge() > 0) continue;
      if (tKaon.charge()*tPion.charge() > 0) continue;
			
      if (d0.DeltaR(pis) > fDeltaR) continue;
			
      dstar = d0 + pis; 
      double dm = dstar.M() - d0.M();
      if (TMath::Abs(0.145 - dm) > fDeltaM) continue;  // keep 0.135 .. 0.155 
      if (dstar.DeltaR(pif) > fDeltaR) continue;
			
      b0 = dstar + pif;
      if (TMath::Abs(b0.M() - MB_0) > 0.6) continue;
      if (b0.DeltaR(pif) > fDeltaR) continue;
      if (b0.DeltaR(d0) > fDeltaR) continue;
			
      // -- sequential fit: D0 pi_slow pi_fast
      if (fVerbose > 5) {
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "==>HFBd2DstarPi> going to sequential fit with track indices: " 
	     << iKaon << " " << iPion << " " << iPionSlow << " " << iPionFast
	     << endl;
      }
			
      HFDecayTree theTree(300000 + fType, true, MB_0, false, -1.0, true);
      theTree.addTrack(iPionFast,211); // add the fast pion to the B0 candidate
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      HFDecayTreeIterator iterator = theTree.addDecayTree(300001 + fType, false, MDSTARPLUS, false); // D*-
      iterator->addTrack(iPionSlow,211); // add the slow pion to the D*+ candidate
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      iterator = iterator->addDecayTree(300002 + fType, true, MD0, false); // D0
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      fSequentialFitter->doFit(&theTree);
    }
  }
} // analyze()

//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2DstarPi);
