#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpPhotons.h"

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

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"


///
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
///

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
HFDumpPhotons::HFDumpPhotons(const edm::ParameterSet& iConfig):
  fPFLabel(iConfig.getUntrackedParameter<InputTag>("pfLabel")),
  fPhotonsLabel(iConfig.getUntrackedParameter<InputTag>("photonsLabel")),  
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fRunOnAOD(iConfig.getUntrackedParameter<bool>("runOnAOD",false))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "---  HFDumpPhotons constructor" << endl;
  cout << "---  fPFLabel:                " << fPFLabel.encode() << endl;
  cout << "---  fPhotonsLabel:           " << fPhotonsLabel.encode() << endl;
  cout << "---  fDoTruthMatching:        " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS
  cout << "---  fVerbose:                " << fVerbose << endl;
  cout << "---  fRunOnAOD:               " << fRunOnAOD << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpPhotons::~HFDumpPhotons() {

}


// ----------------------------------------------------------------------
void HFDumpPhotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  Handle<PFCandidateCollection>  pfCandidates;
  iEvent.getByLabel(fPFLabel, pfCandidates);
  if(!pfCandidates.isValid()) {
    cout << "==>HFDumpPhotons> No PF collection found, skipping" << endl;
    return;
  }
  
  if ( fVerbose > 0 ){
    if ( pfCandidates->size() > 0  ) cout << "PF Collection Size " << pfCandidates->size() << endl;
  }
  
  int i(-1);
  TAnaTrack *pTrack;
  TLorentzVector PFGamma;
  TLorentzVector GenGamma;
  double dR(999);
  double dR_min(999);
  int genIndex(-1), mcid(-1);
  bool match = false;
  for( reco::PFCandidateCollection::const_iterator  iPF = pfCandidates->begin(); iPF != pfCandidates->end(); iPF++) {
    ++i;
    if ( iPF->particleId() == 4 ) {
    //if ( iPF->pdgId() == 22 ) {
      if ( fVerbose > 0 ) cout << " PF Photon Energy  " << iPF->energy() << endl;
    
      pTrack = gHFEvent->addRecTrack();
      pTrack->fIndex  = i;
      pTrack->fQ      = iPF->charge();
      pTrack->fPlab.SetPtEtaPhi(iPF->pt(),
				iPF->eta(),
				iPF->phi()
				);
      
      PFGamma.SetPtEtaPhiM(pTrack->fPlab.Pt(), pTrack->fPlab.Eta(), pTrack->fPlab.Phi(), 0.);
      for (int s = 0; s < gHFEvent->nGenCands() ; ++s) {
	if (gHFEvent->getGenCand(s)->fID == 22 ){
	  GenGamma.SetPxPyPzE(gHFEvent->getGenCand(s)->fP.Px(), gHFEvent->getGenCand(s)->fP.Py(), gHFEvent->getGenCand(s)->fP.Pz(), gHFEvent->getGenCand(s)->fP.E());
	  dR = PFGamma.DeltaR(GenGamma);
	  if ( fVerbose > 0 ) cout << " dR " << dR << endl;
	  if ( (dR < 0.2) && (dR < dR_min) ){
	    genIndex = s;
	    mcid = gHFEvent->getGenCand(s)->fID;
	    dR_min = dR;
	    match = true;
	    if ( fVerbose > 0 ) cout << "Match  dR " << dR << endl;
	  }
	}
      }
      
      if ( match  ){
	pTrack->fGenIndex = genIndex;
	pTrack->fMCID = mcid;
	if ( fVerbose > 0 ) cout << " mcid " << mcid << " genIndex " << genIndex << " dR_min " << dR_min << endl;
      }
      match = false;
      dR=999; dR_min=999; genIndex=-1; mcid=-1;
    }
  }
  
  /*
  Handle<PhotonCollection>  photonCandidates;
  iEvent.getByLabel(fPhotonsLabel, photonCandidates);
  if(!photonCandidates.isValid()) {
    cout << "==>HFDumpPhotons> No PhotonCollection found, skipping" << endl;
    return;
  }
  
  if ( photonCandidates->size() > 0  ) cout << "Photon Collection Size " << photonCandidates->size() << endl;
  
  vector<reco::Photon> localPhotons;
  
  for( reco::PhotonCollection::const_iterator  iPho = photonCandidates->begin(); iPho != photonCandidates->end(); iPho++) {
    
    if ( iPho->energy() < 5. ) cout << " Photon Energy  " << iPho->energy() << endl;
    
  }
  */
  
  
}



// ------------ method called once each job just before starting event loop  ------------
void  HFDumpPhotons::beginJob() {

  gHFFile->cd();
  // H1D *h1 = new TH1D("h2", "h2", 20, 0., 20.);

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpPhotons::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpPhotons);
