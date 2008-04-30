#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSTree.h"


#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

// -- Yikes!
TAna00Event  *gBSEvent;
TFile        *gBSFile;

using namespace::std;

// ----------------------------------------------------------------------
BSTree::BSTree(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)) {

  using namespace std;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- BSTree constructor" << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fFile = new TFile(iConfig.getParameter<string>("fileName").c_str(), "RECREATE");
  fTree = new TTree("T1","CMSSW BS tree");
  fEvent = new TAna00Event(0);

  fTree->Branch("TAna00Event", "TAna00Event", &fEvent, 256000/8, 1);

  gBSEvent = fEvent;
  gBSFile  = fFile;

  TH1D *h1;
  TH2D *h2;

  // -- histograms filled in BSDumpCandidates 
  h1 = new TH1D("eff", "Efficiencies", 1000, 0., 1000. );

  for (int i = 0; i < 4; i++) {

    h1 = new TH1D(Form("m000_%i", i+1), Form("inv. Mass Cand0 (all kaon cand.), sel = %i", i+1), 1000, 0., 10. );
    h1 = new TH1D(Form("m531_%i", i+1), Form("inv. Mass Cand1 (bmm), sel = %i", i+1), 1000, 0., 10. );
    h1 = new TH1D(Form("m443_%i", i+1), Form("inv. Mass Cand2 (jpsi), sel = %i", i+1), 1000, 0., 10. );
    h1 = new TH1D(Form("m521_%i", i+1), Form("inv. Mass Cand3 (bjk), sel = %i", i+1), 1000, 0., 10. );
  }


  // -- histograms filled in L1/HLTriggerReport
  h1 = new TH1D("L1" , "L1-Trigger names",        200, 0., 200.);
  h1 = new TH1D("HLT", "HL-Trigger names",        200, 0., 200.);
 
  // -- histograms filled in BSDumpMuons
  h2  = new TH2D("nmuons", "N_{#mu} / event  ",        100, 0., 100., 100, 0., 100.);

  fNevt = 0;

}


// ----------------------------------------------------------------------
BSTree::~BSTree() {
  
  // -- Save output
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
  delete fFile;
}


// ----------------------------------------------------------------------
void BSTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  fNevt++;

  if (fVerbose > 0) {
    
    cout << endl;
    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
    cout << "==>BSTree> Start event: " << fNevt << endl;
  }

  gBSEvent->fRunNumber   = iEvent.id().run();
  gBSEvent->fEventNumber = iEvent.id().event();


  if (fVerbose > 0)  cout << "BSTree> filling tree for run: " << gBSEvent->fRunNumber
			  << ", event: "  << gBSEvent->fEventNumber 
			  << endl;
  
  fTree->Fill();

  gBSEvent->Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void  BSTree::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  BSTree::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(BSTree);
