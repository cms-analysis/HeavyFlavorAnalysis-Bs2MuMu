#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTree.h"


#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

// -- Yikes!
TAna00Event  *gHFEvent;
TFile        *gHFFile;

using namespace::std;

// ----------------------------------------------------------------------
HFTree::HFTree(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)) {

  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTree constructor" << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fFile = new TFile(iConfig.getParameter<string>("fileName").c_str(), "RECREATE");
  fTree = new TTree("T1","CMSSW HF tree");
  fEvent = new TAna00Event(0);

  fTree->Branch("TAna00Event", "TAna00Event", &fEvent, 256000/8, 1);

  gHFEvent = fEvent;
  gHFFile  = fFile;

  fNevt = 0;

}


// ----------------------------------------------------------------------
HFTree::~HFTree() {
  
  // -- Save output
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
  delete fFile;
}


// ----------------------------------------------------------------------
void HFTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  fNevt++;

  if (fVerbose > 0) {
    
    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
    cout << "==>HFTree> Start event: " << fNevt << endl;
  }

  gHFEvent->fRunNumber   = iEvent.id().run();
  gHFEvent->fEventNumber = iEvent.id().event();


  if (fVerbose > 0)  cout << "HFTree> filling tree for run: " << gHFEvent->fRunNumber
			  << ", event: "  << gHFEvent->fEventNumber 
			  << endl;
  
  fTree->Fill();

  gHFEvent->Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void  HFTree::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFTree::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFTree);
