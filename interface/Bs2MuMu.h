#ifndef Bs2MuMu_h
#define Bs2MuMu_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "SimTracker/TrackAssociation/test/testTrackAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
//#include "SimTracker/Records/interface/VertexAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
//#include "SimTracker/VertexAssociation/interface/VertexAssociatorByTracks.h"
#include "SimGeneral/HepPDT/interface/HepPDTable.h" 

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

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
#include "RecoVertex/VertexTools/interface/BeamSpot.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"

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

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

// ----------------------------------------------------------------------

using namespace edm;
using namespace reco;

class TFile;
class TTree;
class TH1D;

class TAna00Event;
class anaStuff;
class candStuff;

class Bs2MuMu : public edm::EDAnalyzer {

  //class TrackAssociatorBase;
  //class VertexAssociatorBase;


public:
  explicit Bs2MuMu(const edm::ParameterSet&);
  ~Bs2MuMu();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void printGenTracks(const edm::Event&);
  virtual void printSimTracks(const edm::Event&);
  virtual void printRecTracks(const edm::Event&);
  virtual void printMuonTracks(const edm::Event&);
  virtual void printReco2Sim(const edm::Event&, const char *option);

  virtual void decayChannel(const char *name);
  virtual void clearTracks();
  virtual void clearCandidateTracks();

  virtual void fillGeneratorBlock(const edm::Event&);
  virtual void fillRecTracks(const edm::Event&);

  virtual void bmmTracks1(const edm::Event&);
  virtual void bmmTracks2(const edm::Event&);
  virtual void bmmTracks3(const edm::Event&);

  virtual void trimBmmTracks(const edm::Event&);

  int primaryVertex(const edm::Event&);
  virtual void secondaryVertex(const edm::Event&, const edm::EventSetup&);
  double kalmanVertexFit(const edm::Event&, const edm::EventSetup&, int type, unsigned int ntracks);

  virtual void fillTrack(const edm::Event&, TAnaTrack *pTrack, reco::Track *it, int idx, int verb);
  virtual void fillVertex(const edm::Event&, const edm::EventSetup&, TransientVertex *v, int type, unsigned int ntracks);
  virtual void truthCandTracks(const edm::Event&, const edm::EventSetup&);
  virtual void muonCandTracks(const edm::Event&, const edm::EventSetup&, double m_cand1, double m_cand2);
  virtual void kaonCandTracks(const edm::Event&, const edm::EventSetup&, double cone);
  int massMuonCand(const edm::Event&, const edm::EventSetup&, double m_cand1, double m_cand2);
  double rmmKaonCand(const edm::Event&, const edm::EventSetup&, double cone);
  int kaonCandidate(const edm::Event&, const edm::EventSetup&);
  
  int idRecTrack(const reco::Track *track);

  // ----------member data ---------------------------
  string fLabel, fSourceLabel, fTracksLabel, fMuonLabel, fAssocLabel
    , fChannel, fPrintChannel, fPrintChannel2;

  double fMass, fMass2;

  int fVerbose, fGenVerbose, fSimVerbose, fRecVerbose, fGlbVerbose, fR2SVerbose;

  int fNevt, fNgen, fNrec;

  int fTruthMC_I, fTruthMC_II, fTruthMC_mom, fTruthMC_gmo;
  int fTruthMC2, fTruthMC2_mom, fTruthMC2_gmo;

  TFile *fFile; 
  TTree        *fTree;
  TAna00Event  *fEvent;

  anaStuff  *fStuff; // contains most of the class data members

  TrackAssociatorBase*  associatorByChi2;
  TrackAssociatorBase*  associatorByHits;
  //  VertexAssociatorBase* associatorByTracks;
  
  TH1D *fEff;  

  TH1D  *fM000[3], *fM100[3], *fM200[3], *fM300[3];
  TH2D  *fK100, *fK200;
  TH1D  *fPT300, *fPT310, *fPT320;


};

#endif
