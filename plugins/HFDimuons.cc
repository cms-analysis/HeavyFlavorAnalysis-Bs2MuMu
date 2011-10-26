#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDimuons.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

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
HFDimuons::HFDimuons(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fMassLow(iConfig.getUntrackedParameter<double>("massLow", 8.7)), 
  fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 11.2)), 
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.05)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fType(iConfig.getUntrackedParameter<int>("type", 1313)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDimuons constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  massLow:                  " << fMassLow << endl;
  cout << "---  massHigh:                 " << fMassHigh << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "---  vertexing:                " << fVertexing << endl;
  cout << "---  maxDoca                   " << fMaxDoca << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDimuons::~HFDimuons() {
  
}


// ----------------------------------------------------------------------
void HFDimuons::analyze(const Event& iEvent, const EventSetup& iSetup) {
  // -- get the magnetic field
  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  const MagneticField *field = magfield.product();
 
  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDimuons> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFDimuons> No primary vertex found, skipping" << endl;
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
  cout << "==>HFDimuons> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFDimuons> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }


  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im >= 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFDimuons> nMuons = " << hMuons->size() << endl;
    cout << "==>HFDimuons> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < 2) return;
  
  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFDimuons: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- set up vertex fitter 
  HFSequentialVertexFit aSeq(hTracks,fTTB.product(),recoPrimaryVertexCollection, field, fVerbose);
  //   HFKalmanVertexFit a(fTTB.product(), fPV, 1313, 0);
  //   a.setNoCuts();
  //   a.fMaxDoca = fMaxDoca;
  //   vector<Track> trackList; 
  //   vector<int> trackIndices;
  //   vector<double> trackMasses;

  TLorentzVector dimuon, m1, m2;
  int iMuon1(-1), iMuon2(-1); 
  for (unsigned int imuon1 = 0; imuon1 < muonIndices.size()-1; ++imuon1) {
    TrackBaseRef mu1TrackView(hTracks, muonIndices[imuon1]);
    Track tMuon1(*mu1TrackView);
    iMuon1 = muonIndices[imuon1]; 
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 
    if (tMuon1.pt() < fMuonPt)  continue;
    for (unsigned int imuon2 = imuon1 + 1; imuon2 < muonIndices.size(); ++imuon2) {
      TrackBaseRef mu2TrackView(hTracks, muonIndices[imuon2]);
      Track tMuon2(*mu2TrackView);
      iMuon2 = muonIndices[imuon2]; 
      if (tMuon2.pt() < fMuonPt)  continue;

      m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON); 
      dimuon = m1 + m2; 

      if (dimuon.M() < fMassLow || dimuon.M() > fMassHigh) {
	if (fVerbose > 0) {
	  cout << "==>HFDimuons> dimuon mass = " << dimuon.M() << ", skipping" << endl;
	}
	continue; 
      }
      
      //       // -- Vertexing, new style
      //       trackList.clear();
      //       trackIndices.clear(); 
      //       trackMasses.clear(); 
      
      //       trackList.push_back(tMuon1); 
      //       trackIndices.push_back(iMuon1); 
      //       trackMasses.push_back(MMUON);
      
      //       trackList.push_back(tMuon2); 
      //       trackIndices.push_back(iMuon2); 
      //       trackMasses.push_back(MMUON);
      //       if (fVertexing > 0) {
      // 	a.doFit(trackList, trackIndices, trackMasses, fType); 
      //       } else {
      // 	a.doNotFit(trackList, trackIndices, trackMasses, fType); 	
      //       }

      // -- Vertexing, with Kinematic Particles
      HFDecayTree theTree(301313, true, 0, false); // TODO adjust mass to meaningful value
      theTree.addTrack(iMuon1,13);
      theTree.addTrack(iMuon2,13);
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      aSeq.doFit(&theTree);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDimuons::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDimuons::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDimuons);
