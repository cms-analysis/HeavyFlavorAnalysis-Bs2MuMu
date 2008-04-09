#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpStuff.h"

#include "DataFormats/Common/interface/Handle.h"
// #include "DataFormats/METReco/interface/CaloMET.h"
// #include "DataFormats/METReco/interface/CaloMETCollection.h"
// #include "DataFormats/METReco/interface/GenMETCollection.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include <TVector3.h>

// -- Yikes!
extern TAna00Event *gHFEvent;

using namespace std;
using namespace edm;


// ----------------------------------------------------------------------
HFDumpStuff::HFDumpStuff(const edm::ParameterSet& iConfig):
  fGenEventScaleLabel(iConfig.getUntrackedParameter<string>("GenEventScaleLabel", string("genEventScale"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<string>("PrimaryVertexLabel", string("offlinePrimaryVerticesFromCTFTracks")))
//   fCandidates1Label(iConfig.getUntrackedParameter<string>("Candidates1Label", string("JPsiToMuMu"))),
//   fCandidates2Label(iConfig.getUntrackedParameter<string>("Candidates2Label", string("JPsiToMuMu"))),
//   fCandidates3Label(iConfig.getUntrackedParameter<string>("Candidates3Label", string("JPsiToMuMu")))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpStuff constructor" << endl;
  cout << "--- GenEventScaleLabel: " << fGenEventScaleLabel.c_str() << endl;
  cout << "--- PrimaryVertexLabel: " << fPrimaryVertexLabel.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpStuff::~HFDumpStuff() {
  
}


// ----------------------------------------------------------------------
void HFDumpStuff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // -- Primary vertex
  try {
    edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel.c_str(), recoPrimaryVertexCollection);
    const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
    const reco::Vertex pV = vertices[0]; // ???? 
    ChiSquared chi2(pV.chi2(), pV.ndof());
    
    TAnaVertex *pVtx = new TAnaVertex();
    pVtx->setInfo(chi2.value(), int(chi2.degreesOfFreedom()), chi2.probability(), 0, 0);
    pVtx->fPoint.SetXYZ(pV.position().x(),
                        pV.position().y(),
                        pV.position().z());

    gHFEvent->fPrimaryVertex = *pVtx;

  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    cout << "==>HFDumpStuff> primaryVertex " << fPrimaryVertexLabel.c_str() << " not found " << endl;
  } 


  // -- Candidates list
//   try {
//     //    Handle<CandidateCollection> candidates1Handle;
//     Handle<reco::CandidateView> candidates1Handle;
//     iEvent.getByLabel(fCandidates1Label.c_str(), candidates1Handle);
//     for (int i = 0; i < candidates1Handle->size(); ++ i ) {
//       const Candidate &p = (*candidates1Handle)[i];
//       cout << "==>HFDumpStuff>  candidates1 " << i << " id = " << p.pdgId() 
// 	   << " mass = " << p.mass()
// 	   << " pt = " << p.pt() 
// 	   << " phi = " << p.phi() 
// 	   << " eta = " << p.eta() 
// 	   << endl;
//     }
//   } catch (cms::Exception &ex) {
//     //    cout << ex.explainSelf() << endl;
//     cout << "==>HFDumpStuff> Candidate list " << fCandidates1Label.c_str() << " not found " << endl;
//   }

//   try {
//     //    Handle<CandidateCollection> candidates2Handle;
//     Handle<reco::CandidateView> candidates2Handle;
//     iEvent.getByLabel(fCandidates2Label.c_str(), candidates2Handle);
//     for (int i = 0; i < candidates2Handle->size(); ++ i ) {
//       const Candidate &p = (*candidates2Handle)[i];
//       cout << "==> candidates2 " << i << " id = " << p.pdgId() 
// 	   << " mass = " << p.mass()
// 	   << " pt = " << p.pt() 
// 	   << " phi = " << p.phi() 

}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpStuff::beginJob(const edm::EventSetup& setup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpStuff::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFDumpStuff);
