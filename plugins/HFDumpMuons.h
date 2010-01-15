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


class TFile;
class TTree;
class TAna01Event;

using namespace std;
using namespace reco;

// ----------------------------------------------------------------------
class HFDumpMuons : public edm::EDAnalyzer {
 public:
  explicit HFDumpMuons(const edm::ParameterSet&);
  ~HFDumpMuons();

  static int           muonID(const Muon &);
  
 private:
  virtual void         beginJob(const edm::EventSetup&) ;
  virtual void         analyze(const edm::Event&, const edm::EventSetup&);
  virtual void         endJob();
  void                 fillMuon(const Muon& tr, int type);
  vector<unsigned int> muonStatHits(const Track& tr);
  edm::InputTag        fMuonsLabel;

  int                  fVerbose, fDoTruthMatching; 

};

#endif
