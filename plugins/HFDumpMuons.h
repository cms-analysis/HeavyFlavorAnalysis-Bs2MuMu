#ifndef _HFDUMPMUONS_h_
#define _HFDUMPMUONS_h_

#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"


class TFile;
class TTree;
class TAna01Event;

// ----------------------------------------------------------------------
class HFDumpMuons : public edm::EDAnalyzer {
 public:
  explicit HFDumpMuons(const edm::ParameterSet&);
  ~HFDumpMuons();

  
  
  static int                muonID(const reco::Muon &);
  
 private:
  virtual void              beginJob();
  virtual void              beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void              analyze(const edm::Event&, const edm::EventSetup&);
  virtual void              endJob();
  void                      fillMuon(const reco::Muon& tr, int type);
  void                      fillCaloMuon(const reco::CaloMuon& tr, int type);
  std::vector<unsigned int> muonStatHits(const reco::Track& tr);
  edm::InputTag             fMuonsLabel;
  edm::InputTag             fCaloMuonsLabel;

  int                       fVerbose, fDoTruthMatching; 
  bool                      fRunOnAOD;
  
  PropagateToMuon           fpropM1, fpropM2;

};

#endif
