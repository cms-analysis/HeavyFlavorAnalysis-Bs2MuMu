#ifndef _BSDUMPCANDIDATES_h_
#define _BSDUMPCANDIDATES_h_

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// user include files

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 

#include "RecoPixelVertexing/PixelVertexFinding/interface/PVPositionBuilder.h"
#include "RecoPixelVertexing/PixelVertexFinding/interface/PVClusterComparer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "RecoVertex/VertexPrimitives/interface/VertexReconstructor.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

#include <iostream>
#include <string>
#include <map>
#include <set>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class TFile;
class TTree;
class TH1D;
class TAna00Event;

class TrackAssociatorBase;
/* class VertexAssociatorBase; */

// ----------------------------------------------------------------------
class BSDumpCandidates : public edm::EDAnalyzer {

public:
  explicit BSDumpCandidates(const edm::ParameterSet&);
  ~BSDumpCandidates();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void decayChannel(const char *name);
  virtual void clearTracks();
  virtual void clearCandidateTracks();

  virtual void bmmTracks1(const edm::Event&);
  virtual void bmmTracks2(const edm::Event&);
  virtual void bmmTracks3(const edm::Event&);
  virtual void bmmTracks4(const edm::Event&);

  int primaryVertex(const edm::Event&);
  virtual void secondaryVertex(const edm::Event&, const edm::EventSetup&);
  double kalmanVertexFit(const edm::Event&, const edm::EventSetup&, int type, unsigned int ntracks);

  virtual void fillVertex(const edm::Event&, const edm::EventSetup&, TransientVertex *v, int type, unsigned int ntracks);
  virtual void truthCandTracks(const edm::Event&, const edm::EventSetup&);
  virtual void muonCandTracks(const edm::Event&, const edm::EventSetup&, double m_cand1, double m_cand2);
  virtual void kaonCandTracks(const edm::Event&, const edm::EventSetup&, double cone);

  int massMuonCand(const edm::Event&, const edm::EventSetup&, double m_cand1, double m_cand2);
  double rmmKaonCand(const edm::Event&, const edm::EventSetup&, double cone);
  int kaonCandidate(const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  std::string fGenEventLabel;
  std::string fTracksLabel;
  std::string fVtxAssociatorLabel;
  std::string fAssociatorLabel;
  std::string fTrackingParticlesLabel;
  std::string fTrackingVertexLabel;

  std::string fPrimaryVertexLabel;
  std::string fChannel;
  std::string fPrintChannel, fPrintChannel2;

  edm::InputTag fMuonsLabel1;
  edm::InputTag fMuonsLabel2;

  double fMass, fMass2;
  double fMassRange, fMassRange2, fMassKaon;

  int fBmmSel;
  int fVerbose;
  int fDoTruthMatching;

  int fNevt, fNgen, fNrec;

  int fTruthMC_I, fTruthMC_II, fTruthMC_mom, fTruthMC_gmo;
  int fTruthMC2, fTruthMC2_mom, fTruthMC2_gmo;

  TFile *fFile; 

  TrackAssociatorBase  *fAssociator;
/*   VertexAssociatorBase *fVtxAssociator; */
  

  // -- Stuff -----------------------------------------------------

  std::vector<double> SecVtxChi2;
  std::vector<double> InvMass;

  // ---------------
  // -- Collections
  // ---------------
 
  const reco::MuonCollection       *theGlobalMuonCollection;
  const reco::MuonCollection       *theTrackerMuonCollection;
  const reco::TrackCollection      *theTkCollection;
  std::vector<HepMC::GenParticle*> theGenCollection;

  reco::RecoToSimCollection        *recSimCollection;
  reco::VertexRecoToSimCollection  *recSimCollectionVertex;

  // ----------
  // -- Tracks
  // ----------

  // -- Muons (TM=2 or GM=3) 
  std::vector<const reco::Track*> MuonRecTracks; 
  std::vector<int> MuonRecTracksIndex;

  // -- Muons (full TM=1)
  std::vector<const reco::Track*> BmmRecTracks;
  std::vector<int> BmmRecTracksIndex;          
  std::vector<int> BmmRecTracksB;                 // Index of B-mother

  std::vector<const reco::Track*> JpsiRecTracks; 
  std::vector<int> JpsiRecTracksIndex;
  std::vector<int> JpsiRecTracksB;                // Index of B-grandmother

  // -- Kaons (full TM=1)
  std::vector<const reco::Track*> KaonRecTracks; 
  std::vector<int> KaonRecTracksIndex;
  std::vector<int> KaonRecTracksB;                // Index of B-mother


  // ---------------------------
  // -- Track candidates (pairs)
  // ----------------------------

  // -- Muons signal
  std::vector< std::pair<const reco::Track*, const reco::Track*> > BmmPairTracks;
  std::vector< std::pair<int, int> > BmmPairTracksIndex;

  // -- Muons norm
  std::vector< std::pair<const reco::Track*, const reco::Track*> > JpsiPairTracks;
  std::vector< std::pair<int, int> > JpsiPairTracksIndex;

  // -- Kaons norm
  std::vector<const reco::Track*> KaonTrack; 
  std::vector<int> KaonTrackIndex;

  // -----------------
  // -- Vertex tracks
  // -----------------
  std::vector<const reco::Track*> RecTracks; 
  std::vector<int> RecTracksIndex;
  std::vector<reco::Track> RefittedTracks;


  // ------------------
  // -- Primary Vertex
  // ------------------

  reco::Vertex primaryVertex2;

};

#endif
