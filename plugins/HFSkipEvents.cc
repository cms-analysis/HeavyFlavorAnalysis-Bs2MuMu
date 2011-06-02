#include <memory>

#include "DataFormats/TrackReco/interface/Track.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

class HFSkipEvents : public EDFilter {
public:
  explicit HFSkipEvents(const ParameterSet&);
  ~HFSkipEvents();

private:
  virtual void beginJob();
  virtual void endJob() ;
  virtual bool filter(Event&, const EventSetup&);

  int          fVerbose;

  int          filterOnPrimaryVertex; 
  InputTag     fPrimaryVertexCollectionLabel;

  int          filterOnTracks; 
  InputTag     fTrackCollectionLabel;

  int          filterOnMuons; 
  InputTag     fMuonCollectionLabel;

  int fNpv, fNtk, fNmu; 
  int fNfailed, fNpassed; 
  int fEvent; 

};

// ----------------------------------------------------------------------
HFSkipEvents::HFSkipEvents(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  filterOnPrimaryVertex(iConfig.getUntrackedParameter<int>("filterOnPrimaryVertex", 1)),
  fPrimaryVertexCollectionLabel(iConfig.getUntrackedParameter<InputTag>("primaryVertexCollectionLabel", edm::InputTag("offlinePrimaryVertices"))),
  filterOnTracks(iConfig.getUntrackedParameter<int>("filterOnTrackMaximum", 1)),
  fTrackCollectionLabel(iConfig.getUntrackedParameter<InputTag>("TrackCollectionLabel", edm::InputTag("generalTracks"))),
  filterOnMuons(iConfig.getUntrackedParameter<int>("filterOnMuonMinimum", 1)),
  fMuonCollectionLabel(iConfig.getUntrackedParameter<InputTag>("MuonCollectionLabel", edm::InputTag("muons")))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFSkipEvents constructor" << endl;
  cout << "---  verbose:                         " << fVerbose << endl;
  cout << "---  filterOnPrimaryVertex:           " << filterOnPrimaryVertex << endl;
  cout << "---  primaryVertexCollectionLabel:    " << fPrimaryVertexCollectionLabel << endl;
  cout << "---  filterOnTrackMaximum:            " << filterOnTracks << endl;
  cout << "---  filterOnMuonMinimum:             " << filterOnMuons << endl;
  cout << "---  TrackCollectionLabel:            " << fTrackCollectionLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fNpv = fNtk = fNmu = 0; 
  fEvent = fNfailed = fNpassed = 0; 
  
}


// ----------------------------------------------------------------------
HFSkipEvents::~HFSkipEvents() {

}



// ----------------------------------------------------------------------
bool HFSkipEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {    
  bool result(true);
  ++fEvent; 

  edm::Handle<reco::VertexCollection> hVertices;
  int goodVertices(-1);
  try{ 
    iEvent.getByLabel(fPrimaryVertexCollectionLabel, hVertices);
    const reco::VertexCollection vertices = *(hVertices.product());
    goodVertices = 0; 
    for (unsigned int i = 0; i < vertices.size(); ++i) {
      if (vertices[i].isFake()) {
      } else {
	++goodVertices;
      }
    }
  } catch (cms::Exception &ex) {
    if (fVerbose > 1) cout << "No Vertex collection with label " << fPrimaryVertexCollectionLabel << endl;
  }

  edm::Handle<MuonCollection> hMuons;
  int goodMuons(-1);
  try{ 
    iEvent.getByLabel(fMuonCollectionLabel, hMuons);
    goodMuons = (*(hMuons.product())).size();
  } catch (cms::Exception &ex)  {
    if (fVerbose > 1) cout << "No Muon collection with label " << fMuonCollectionLabel << endl;
  }
  bool goodMu = (goodMuons >= filterOnMuons);
  if (filterOnMuons > 0) {
    if (goodMu) {
      //      result = true;
      ++fNmu;
    } else {
      result = false;
    }
  }


  edm::Handle<reco::TrackCollection> hTracks;
  int goodTracks(-1);
  try{ 
    iEvent.getByLabel(fTrackCollectionLabel, hTracks);
    goodTracks = (*(hTracks.product())).size();
  } catch (cms::Exception &ex) {
    if (fVerbose > 1) cout << "No Track collection with label " << fTrackCollectionLabel << endl;
  }


  bool goodPv = (goodVertices >= filterOnPrimaryVertex);
  if (filterOnPrimaryVertex > 0)  {
    if (goodPv) {
      //      result = true;
      ++fNpv;
    } else {
      result = false; 
    }
  }

  bool goodTk = (goodTracks < filterOnTracks);
  if (filterOnTracks > 0) {
    if (goodTk) {
      //      result = true;
      ++fNtk;
    } else {
      result = false;
    }
  }


  if (fVerbose > 0) {
    char line[20]; 
    sprintf(line, "%7d", fEvent);
    cout << "HFSkipEvents: " << line
	 << " result: " << (result? "true ":"false")
	 << " PV: " << (goodPv?1:0) << "(" << goodVertices << ")"
	 << " Tk: " << (goodTk?1:0) << "(" << goodTracks << ")"
	 << " Mu: " << (goodMu?1:0) << "(" << goodMuons << ")"
	 << endl;
  }

  if (result) {
    ++fNpassed; 
  } else {
    ++fNfailed;
  }
  
  return result;
}


// --------------------------------------------------------------------------------
void  HFSkipEvents::beginJob() {
}


// --------------------------------------------------------------------------------
void  HFSkipEvents::endJob() {

    cout << "HFSkipEvents: " 
	 << " passed events: " << fNpassed
	 << " failed events: " << fNfailed
	 << " fNpv: " << fNpv
	 << " fNtk: " << fNtk
	 << endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFSkipEvents);
