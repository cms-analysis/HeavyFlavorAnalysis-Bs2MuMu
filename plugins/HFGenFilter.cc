#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFGenFilter.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>


#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace edm;
using namespace std;

// ----------------------------------------------------------------------
HFGenFilter::HFGenFilter(const edm::ParameterSet& iConfig) {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFGenFilter constructor" << endl;
  cout << "----------------------------------------------------------------------" << endl;
  
  fNgood= 0;
  fNtot = 0; 
}


// ----------------------------------------------------------------------
HFGenFilter::~HFGenFilter() {  
  cout << "==> HFGenFilter> Total number of accepted/total events = " << fNgood << "/" << fNtot << endl;
}


// ----------------------------------------------------------------------
bool HFGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bool goodEvent(false);

  TGenCand *pGen;

  ++fNtot;
  //cout << "==> HFGenFilter> new  event  " << fNtot << endl;
  //cout << "==> HFGenFilter> nGenCands = " << gHFEvent->nGenCands() << endl;

  // -- filter events with at least one muon and pT>2.5GeV
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if ((13 == TMath::Abs(pGen->fID)) && (pGen->fP.Perp() > 2.5)) {
      goodEvent = true;
    }
  }

  if (goodEvent) {
    ++fNgood;
    //cout << " accepted" << endl;
  } else {
    //cout << " not accepted" << endl;
  }
  return goodEvent;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFGenFilter);
