#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBd2JpsiKstar.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

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

#define MMUON 0.10566
#define MKAON 0.49368
#define MJPSI 3.096916
#define MPION 0.139570
#define MPROTON 0.93827
#define MKSTAR 0.8960
#define MKSHORT 0.497614
#define MB_0 5.27953
#define MBPLUS 5.27927


// ----------------------------------------------------------------------
HFBd2JpsiKstar::HFBd2JpsiKstar(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fKaonPt(iConfig.getUntrackedParameter<double>("kaonPt", 0.4)), 
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.)),
  fType(iConfig.getUntrackedParameter<int>("type", 511))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBd2JpsiKstar constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  kaonPt:                   " << fKaonPt << endl;
  cout << "---  pionPt:                   " << fPionPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFBd2JpsiKstar::~HFBd2JpsiKstar() {
  
}


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- get the primary vertex
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFBd2JpsiKstar> No primary vertex collection found, skipping" << endl;
    return;
  }
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFBd2JpsiKstar> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fEventTag]; 
  if (fVerbose > 0) {
    cout << "HFDimuons: Taking vertex " << gHFEvent->fEventTag << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of muons
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  if (!hMuons.isValid()) {
  cout << "==>HFBd2JpsiKstar> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  edm::Handle<edm::View<reco::Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBd2JpsiKstar> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im > 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFBd2JpsiKstar> nMuons = " << hMuons->size() << endl;
    cout << "==>HFBd2JpsiKstar> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < 2) return;

  // -- reconstruct J/Psi, based on 'best' dimuon mass
  TLorentzVector psi, cpsi, m1, m2;
  double dmass(1000.0);
  Track tMuon1, tMuon2; 
  int iMuon1(-1), iMuon2(-1); 
  for (unsigned int imuon1 = 0; imuon1 < muonIndices.size()-1; ++imuon1) {
    TrackBaseRef mu1TrackView(hTracks, muonIndices[imuon1]);
    Track muon1(*mu1TrackView);
    if (muon1.pt() < fMuonPt)  continue;

    for (unsigned int imuon2 = imuon1 + 1; imuon2 < muonIndices.size(); ++imuon2) {
      TrackBaseRef mu2TrackView(hTracks, muonIndices[imuon2]);
      Track muon2(*mu2TrackView);
      if (muon2.pt() < fMuonPt)  continue;

      m1.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), MMUON); 
      m2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), MMUON); 
      cpsi = m1 + m2; 
      if ((TMath::Abs(psi.M()-MJPSI) < dmass) 
	  //&& (muon1.charge()*muon2.charge() < 0) 
	  ) {
	psi    = cpsi; 
	dmass  = TMath::Abs(psi.M()-MJPSI);
	tMuon1 = muon1;
	iMuon1 = muonIndices[imuon1]; 
	tMuon2 = muon2;
	iMuon2 = muonIndices[imuon2]; 
      }
    }
  }

  if (iMuon1 < 0 || iMuon2 < 0) {
    //if (fVerbose > 0) cout << " --> HFBd2JpsiKstar: No unlike charge dimuons found, dmass = " << dmass << endl;
    return;
  }
  if (dmass > 0.5) {
    if (fVerbose > 0) cout << " --> HFBd2JpsiKstar: No J/Psi candidate found: " << psi.M() << endl;
    return;
  }

  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFBd2JpsiKstar: Error: no TransientTrackBuilder found."<<endl;
    return;
  }
  std::vector<reco::Track> fitTracks;
  
  TAnaVertex BdVtx;
  TAnaTrack *pTrack; 

  TLorentzVector v1, v2, pi, ka, kstar1, kstar2, kstar, b0; 
  Track tKaon, tPion; 
  int   iKaon, iPion; 

  // -- Loop over all K* daughters: include both K* possibilities 
  // -- kaon!
  for (unsigned int iTrack1 = 0; iTrack1 < hTracks->size(); ++iTrack1){    
    if (static_cast<int>(iTrack1) == iMuon1 || static_cast<int>(iTrack1) == iMuon2) continue; 
    TrackBaseRef rTrackView(hTracks, iTrack1);
    Track track1(*rTrackView);
    
    if (track1.pt() < fKaonPt) continue;
    
    // -- pion!
    for (unsigned int iTrack2 = 0; iTrack2 < hTracks->size(); ++iTrack2){    
      if (iTrack2 == iTrack1) continue; 
      if (static_cast<int>(iTrack2) == iMuon1 || static_cast<int>(iTrack2) == iMuon2) continue; 
      TrackBaseRef sTrackView(hTracks, iTrack2);
      Track track2(*sTrackView);

      if (track2.pt() < fPionPt) continue;
      if (track1.charge()*track2.charge() > 0) continue;
	
      tPion = track2; iPion = iTrack2; 
      tKaon = track1; iKaon = iTrack1; 
      v1.SetXYZM(track1.px(), track1.py(), track1.pz(), MKAON); 
      v2.SetXYZM(track2.px(), track2.py(), track2.pz(), MPION); 
      kstar1 = v1 + v2;

      //       v1.SetXYZM(track1.px(), track1.py(), track1.pz(), MPION); 
      //       v2.SetXYZM(track2.px(), track2.py(), track2.pz(), MKAON); 
      //       kstar2 = v1 + v2;
      
      //       if (TMath::Abs(kstar1.M() - MKSTAR) > TMath::Abs(kstar2.M() - MKSTAR)) {
      // 	tPion = track1; iPion = iTrack1; 
      // 	tKaon = track2; iKaon = iTrack2; 
      // 	pi = v1; 
      // 	ka = v2; 
      // 	kstar = kstar2; 
      //       } else {
      // 	tPion = track2; iPion = iTrack2; 
      // 	tKaon = track1; iKaon = iTrack1; 
      // 	pi = v2; 
      // 	ka = v1; 
      // 	kstar = kstar1; 
      //       }

      kstar = kstar1; 
      if (TMath::Abs(kstar.M() - MKSTAR) >  0.7) continue; 
      b0 = kstar + psi; 
      if (TMath::Abs(b0.M() - MB_0) > 1.2) continue; 

      fitTracks.clear();
      fitTracks.push_back(tMuon1); 
      fitTracks.push_back(tMuon2); 
      fitTracks.push_back(tKaon); 
      fitTracks.push_back(tPion); 

      RefCountedKinematicTree bdTree;
      int OK = doVertexFit(fitTracks, bdTree); 

      if (OK == 0 || !bdTree->isValid()) {
	if (fVerbose > 0) {
	  cout << "----------------------------------------" << endl;
	  cout << "==> HFBd2JpsiKstar: Empty/Invalid kinematic fit tree !! " << endl;
	  cout << tMuon1.pt() << ", "  << tMuon1.d0() << ", "  << tMuon1.dz() << endl;
	  cout << tMuon2.pt() << ", "  << tMuon2.d0() << ", "  << tMuon2.dz() << endl;
	  cout << tKaon.pt() << ", "  << tKaon.d0() << ", "  << tKaon.dz() << endl;
	  cout << tPion.pt() << ", "  << tPion.d0() << ", "  << tPion.dz() << endl;
	  cout << "----------------------------------------" << endl;
	}
	continue; 
      }

      bdTree->movePointerToTheTop();
      RefCountedKinematicVertex bdVertex = bdTree->currentDecayVertex();
      RefCountedKinematicParticle bd = bdTree->currentParticle();

      TVector3 bdPlab = (bd->currentState().globalMomentum().x(), 
			 bd->currentState().globalMomentum().y(), 
			 bd->currentState().globalMomentum().z()
			 );
      AlgebraicVector7 bdParameters = bd->currentState().kinematicParameters().vector();
      double bdMass = bdParameters[6];

      if (bdVertex->vertexIsValid() && bdMass > 4.8 && bdMass < 6.0) {

	if (fVerbose > 0) {
	  cout << "----------------------------------------" << endl;
	  cout << "==> HFBd2JpsiKstar: Filling B0 candidate with mass = " << bdMass << endl;
	  cout << tMuon1.pt() << ", "  << tMuon1.d0() << ", "  << tMuon1.dz() << endl;
	  cout << tMuon2.pt() << ", "  << tMuon2.d0() << ", "  << tMuon2.dz() << endl;
	  cout << tKaon.pt() << ", "  << tKaon.d0() << ", "  << tKaon.dz() << endl;
	  cout << tPion.pt() << ", "  << tPion.d0() << ", "  << tPion.dz() << endl;
	  cout << "----------------------------------------" << endl;
	}

	TAnaVertex anaVt;

	ChiSquared chi(bd->chiSquared(), bd->degreesOfFreedom());
	anaVt.setInfo(bd->chiSquared(), bd->degreesOfFreedom(), chi.probability(), 0, 0);
	
	anaVt.fPoint.SetXYZ(bdVertex->position().x(), 
			    bdVertex->position().y(), 
			    bdVertex->position().z()
			    );
	
	// -- Distance to primary vertex
	VertexDistanceXY axy;
	double dXY      = axy.distance(fPV, bdVertex->vertexState()).value();
	double dXYE     = axy.distance(fPV, bdVertex->vertexState()).error();
	
	VertexDistance3D a3d;
	double d3d      = a3d.distance(fPV, bdVertex->vertexState()).value();
	double d3dE     = a3d.distance(fPV, bdVertex->vertexState()).error();
	
	anaVt.fDxy  = dXY; 
	anaVt.fDxyE = dXYE; 
	
	anaVt.fD3d  = d3d; 
	anaVt.fD3dE = d3dE; 
	anaVt.addTrack(iMuon1); 
	anaVt.addTrack(iMuon2); 
	anaVt.addTrack(iKaon); 
	anaVt.addTrack(iPion); 


	// -- fill candidate
	TAnaCand  *pCand = gHFEvent->addCand();
	
	pCand->fPlab = bdPlab;
	pCand->fMass = bdMass;
	pCand->fVtx  = anaVt;    
	pCand->fType = fType;
	pCand->fDau1 = -1;
	pCand->fDau2 = -1;
	pCand->fSig1 = gHFEvent->nSigTracks();
	pCand->fSig2 = pCand->fSig1 + 3;

	// -- fill (not yet refitted) sig tracks
	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tMuon1.charge()*13; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tMuon1.charge();
	pTrack->fPlab.SetPtEtaPhi(tMuon1.pt(),
				  tMuon1.eta(),
				  tMuon1.phi()
				  ); 
	pTrack->fIndex  = iMuon1;

	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tMuon2.charge()*13; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tMuon2.charge();
	pTrack->fPlab.SetPtEtaPhi(tMuon2.pt(),
				  tMuon2.eta(),
				  tMuon2.phi()
				  ); 
	pTrack->fIndex  = iMuon2;

	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tKaon.charge()*321; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tKaon.charge();
	pTrack->fPlab.SetPtEtaPhi(tKaon.pt(),
				  tKaon.eta(),
				  tKaon.phi()
				  ); 
	pTrack->fIndex  = iKaon;
	
	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tPion.charge()*211; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tPion.charge();
	pTrack->fPlab.SetPtEtaPhi(tPion.pt(),
				  tPion.eta(),
				  tPion.phi()
				  ); 
	pTrack->fIndex  = iPion;
	
	// -- Vertex the two muons only
	pCand = gHFEvent->addCand(); 
	fitTracks.clear();
	fitTracks.push_back(tMuon1); 
	fitTracks.push_back(tMuon2); 
	doJpsiVertexFit(fitTracks, iMuon1, iMuon2, pCand); 

      }
    }

  }
  
}


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::doJpsiVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2, TAnaCand *pCand){

  Track tMuon1 = Tracks[0]; 
  Track tMuon2 = Tracks[1]; 
  
  std::vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  RecoTransientTrack.push_back(fTTB->build(Tracks[0]));
  RecoTransientTrack.push_back(fTTB->build(Tracks[1]));
  
  // -- Do the vertexing
  KalmanVertexFitter theFitter(true);
  TransientVertex TransSecVtx = theFitter.vertex(RecoTransientTrack); 
  if (TransSecVtx.isValid() ) {
    if (isnan(TransSecVtx.position().x()) 
	|| isnan(TransSecVtx.position().y()) 
	|| isnan(TransSecVtx.position().z()) ) {
      cout << "==>HFDimuons> Something went wrong! SecVtx nan - continue ... " << endl;
      pCand->fType = -1;
      return; 
    }
  } else {
    cout << "==>HFDimuons> KVF failed! continue ..." << endl;
    pCand->fType = -1;
    return; 
  }
  
  // -- Get refitted tracks
  std::vector<reco::TransientTrack> refTT = TransSecVtx.refittedTracks();
  std::vector<reco::Track> refT; refT.clear(); 
  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); i != refTT.end(); i++) {
    const Track &ftt = i->track();
    refT.push_back(ftt);
  }
  
  // -- Build composite
  TLorentzVector comp, M1, M2;
  M1.SetXYZM(refT[0].px(), refT[0].py(), refT[0].pz(), MMUON); 
  M2.SetXYZM(refT[1].px(), refT[1].py(), refT[1].pz(), MMUON); 
  comp = M1 + M2;
  
  
  // -- Build vertex for ntuple
  TAnaVertex anaVtx;
  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());
  anaVtx.setInfo(chi.value(), chi.degreesOfFreedom(), chi.probability(), 0, 0);
  anaVtx.fPoint.SetXYZ(TransSecVtx.position().x(), 
		       TransSecVtx.position().y(), 
		       TransSecVtx.position().z());
  
  anaVtx.addTrack(iMuon1);
  anaVtx.addTrack(iMuon2);
  
  VertexDistanceXY axy;
  anaVtx.fDxy     = axy.distance(fPV, TransSecVtx).value();
  anaVtx.fDxyE    = axy.distance(fPV, TransSecVtx).error();
  VertexDistance3D a3d;
  anaVtx.fD3d     = a3d.distance(fPV, TransSecVtx).value();
  anaVtx.fD3dE    = a3d.distance(fPV, TransSecVtx).error();
  
  
  // -- fill candidate
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();
  pCand->fVtx  = anaVtx;    
  pCand->fType = fType*1000 + 443;
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + 1;
  
  // -- fill refitted sig tracks
  TAnaTrack *pTrack = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon1.charge()*13; 
  pTrack->fGenIndex = -1; 
  pTrack->fQ        = tMuon1.charge();
  pTrack->fPlab.SetXYZ(refT[0].px(),
		       refT[0].py(),
		       refT[0].pz()
		       ); 
  pTrack->fIndex  = iMuon1;
  
  pTrack            = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon2.charge()*13; 
  pTrack->fGenIndex = -1; 
  pTrack->fQ        = tMuon2.charge();
  pTrack->fPlab.SetXYZ(refT[1].px(),
		       refT[1].py(),
		       refT[1].pz()
		       ); 
  pTrack->fIndex  = iMuon2;
}


// ----------------------------------------------------------------------
int HFBd2JpsiKstar::doVertexFit(std::vector<reco::Track> &Tracks, RefCountedKinematicTree &bdTree){

  KinematicParticleFactoryFromTransientTrack pFactory;
  
  ParticleMass muon_mass = MMUON;
  ParticleMass kaon_mass = MKAON;
  ParticleMass pion_mass = MPION;
  float muon_sigma = 0.0000001;
  float kaon_sigma = 0.000016;
  float pion_sigma = 0.000016;

  float chi = 0.;
  float ndf = 0.;

  TransientTrack ttMuon1 = fTTB->build(Tracks[0]);
  TransientTrack ttMuon2 = fTTB->build(Tracks[1]);
  TransientTrack ttKaon  = fTTB->build(Tracks[2]);
  TransientTrack ttPion  = fTTB->build(Tracks[3]);
  
  vector<RefCountedKinematicParticle> muonParticles;
  vector<RefCountedKinematicParticle> kstarParticles;
  vector<RefCountedKinematicParticle> allParticles;
  muonParticles.push_back(pFactory.particle(ttMuon1, muon_mass, chi, ndf, muon_sigma));
  muonParticles.push_back(pFactory.particle(ttMuon2, muon_mass, chi, ndf, muon_sigma));
  allParticles.push_back(pFactory.particle(ttKaon,   muon_mass, chi, ndf, muon_sigma));
  allParticles.push_back(pFactory.particle(ttPion,   muon_mass, chi, ndf, muon_sigma));
  
  kstarParticles.push_back(pFactory.particle(ttKaon, kaon_mass, chi, ndf, kaon_sigma));
  kstarParticles.push_back(pFactory.particle(ttPion, pion_mass, chi, ndf, pion_sigma));
  allParticles.push_back(pFactory.particle(ttKaon,   kaon_mass, chi, ndf, kaon_sigma));
  allParticles.push_back(pFactory.particle(ttPion,   pion_mass, chi, ndf, pion_sigma));
  
  KinematicParticleVertexFitter kpvFitter;
  RefCountedKinematicTree jpTree; 
  try {
    jpTree = kpvFitter.fit(muonParticles);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBd2JpsiKstar: exception caught from kin fit 1" << endl;
  }
  if (jpTree->isEmpty()) {
    return 0;
  } 

  KinematicParticleFitter csFitter;
  ParticleMass jpsi = MJPSI;
  float jp_m_sigma = 0.00004;
  KinematicConstraint *jpsi_c2 = new MassKinematicConstraint(jpsi, jp_m_sigma);
  try {
    jpTree = csFitter.fit(jpsi_c2, jpTree);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBd2JpsiKstar: exception caught from kin fit 2" << endl;
    delete jpsi_c2; 
    return 0; 
  }

  if (jpTree->isEmpty()) return 0; 
  
  try {
    jpTree->movePointerToTheTop();
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBd2JpsiKstar: exception caught from kin fit " << endl;
    delete jpsi_c2; 
    return 0; 
  }

  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
  kstarParticles.push_back(jpsi_part);
  
  bdTree = kpvFitter.fit(kstarParticles);
  delete jpsi_c2; 
  return 1;
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBd2JpsiKstar::beginJob(const edm::EventSetup& setup) {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBd2JpsiKstar::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2JpsiKstar);
