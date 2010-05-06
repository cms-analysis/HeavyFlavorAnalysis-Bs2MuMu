#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTrigger.h"

#include <iostream>
#include <bitset>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TTrgObj.hh"

#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"


// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


// ----------------------------------------------------------------------
HFDumpTrigger::HFDumpTrigger(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fHLTProcessName(iConfig.getUntrackedParameter<string>("HLTProcessName")),
  fL1GTReadoutRecordLabel(iConfig.getUntrackedParameter<InputTag>("L1GTReadoutRecordLabel", edm::InputTag("gtDigis"))),
  fL1GTmapLabel(iConfig.getUntrackedParameter<InputTag>("hltL1GtObjectMap")),
  fL1MuonsLabel(iConfig.getUntrackedParameter<InputTag>("L1MuonsLabel")),
  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel"))
{

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrigger constructor" << endl;
  cout << "--- Verbose                     : " << fVerbose << endl;
  cout << "--- HLT process name            : " << fHLTProcessName << endl;
  cout << "--- L1 GT Readout Record Label  : " << fL1GTReadoutRecordLabel << endl;
  cout << "--- L1 GT Object Map Label      : " << fL1GTmapLabel << endl;
  cout << "--- L1 Muons Label              : " << fL1MuonsLabel << endl;
  cout << "--- HLTResultsLabel             : " << fHLTResultsLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0; 
  
}


// ----------------------------------------------------------------------
HFDumpTrigger::~HFDumpTrigger() {
  
}


// ----------------------------------------------------------------------
void HFDumpTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  // ----------------------------------------------------------------------
  // -- L1 results: physics and technical triggers
  // ----------------------------------------------------------------------
  gHFEvent->fL1Decision = 0; 
  gHFEvent->fL1TWords[0]=0; gHFEvent->fL1TWords[1]=0; gHFEvent->fL1TWords[2]=0; gHFEvent->fL1TWords[3]=0; 
  gHFEvent->fL1TWasRun[0]=0;gHFEvent->fL1TWasRun[1]=0;gHFEvent->fL1TWasRun[2]=0;gHFEvent->fL1TWasRun[3]=0; 
  gHFEvent->fL1TTWords[0]=0;gHFEvent->fL1TTWords[1]=0;

  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(fL1GTReadoutRecordLabel,L1GTRR);
  Handle<L1GlobalTriggerObjectMapRecord> hL1GTmap; 
  iEvent.getByLabel("hltL1GtObjectMap", hL1GTmap);

  L1GtUtils l1GtUtils;
  l1GtUtils.retrieveL1EventSetup(iSetup);
  // cout << "L1 trigger menu: ";
  // cout << l1GtUtils.l1TriggerMenu() << endl;

  edm::ESHandle<L1GtTriggerMenu> hL1GtMenu;
  iSetup.get<L1GtTriggerMenuRcd>().get(hL1GtMenu);
  const L1GtTriggerMenu* l1GtMenu = hL1GtMenu.product();

  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  string algoname; 
  int    algobit(-1); 
  bool   result(false); 
  int    prescale(0); 
  int    mask(0); 
  int    itrig(-1), iword(-1), iErrorCode; 

  for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);

    iword = algobit/32;
    itrig = algobit%32;
    
    gHFEvent->fL1TNames[algobit] = TString(algoname);
    if (result)    gHFEvent->fL1TWords[iword]  |= (0x1 << itrig);
    if (0 == mask) gHFEvent->fL1TWasRun[iword] |= (0x1 << itrig);
    gHFEvent->fL1TPrescale[algobit] = prescale;
  }

  for (CItAlgo algo = menu->gtTechnicalTriggerMap().begin(); algo != menu->gtTechnicalTriggerMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);
    
    iword = algobit/32;
    itrig = algobit%32;
    
    gHFEvent->fL1TTNames[algobit] = TString(algoname);
    if (result)    gHFEvent->fL1TTWords[iword]  |= (0x1 << itrig);
    if (0 == mask) gHFEvent->fL1TTWasRun[iword] |= (0x1 << itrig);
    gHFEvent->fL1TTPrescale[algobit] = prescale;

  }

  // ----------------------------------------------------------------------
  // -- HLT results
  // ----------------------------------------------------------------------
  gHFEvent->fHLTDecision = 0; 
  gHFEvent->fHLTWords[0]=0; gHFEvent->fHLTWords[1]=0; gHFEvent->fHLTWords[2]=0; gHFEvent->fHLTWords[3]=0; 
  gHFEvent->fHLTWords[4]=0; gHFEvent->fHLTWords[5]=0; gHFEvent->fHLTWords[6]=0; gHFEvent->fHLTWords[7]=0; 

  gHFEvent->fHLTWasRun[0]=0;gHFEvent->fHLTWasRun[1]=0;gHFEvent->fHLTWasRun[2]=0;gHFEvent->fHLTWasRun[3]=0; 
  gHFEvent->fHLTWasRun[4]=0;gHFEvent->fHLTWasRun[5]=0;gHFEvent->fHLTWasRun[6]=0;gHFEvent->fHLTWasRun[7]=0; 


  // -- Read HLT configuration and names
  HLTConfigProvider hltConfig;
  bool hltConfigInitSuccess = hltConfig.init(fHLTProcessName);

  vector<string> validTriggerNames;
  if (hltConfigInitSuccess) validTriggerNames = hltConfig.triggerNames();

  if (validTriggerNames.size() < 1) {
    cout << "==>HFDumpTrigger: NO valid trigger names returned by HLT config provided!!??" << endl;
    return;
  }

  Handle<TriggerResults> hHLTresults;
  bool hltF = true;
  try {
    iEvent.getByLabel(fHLTResultsLabel, hHLTresults);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDumpTrigger> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    hltF = false;
  }
  
  if (hltF) {
    TriggerNames trigName;
    trigName.init(*hHLTresults);
    gHFEvent->fHLTDecision = hHLTresults->accept();
    if (fVerbose > 1) cout << "hHLTresults->size() = " << hHLTresults->size() << " and HLT accept = " << gHFEvent->fHLTDecision << endl;

    unsigned int index(999); 
    int wasrun(0), result(0), error(0);
    for (unsigned int it = 0; it < validTriggerNames.size(); ++it) {
      index = trigName.triggerIndex(validTriggerNames[it]); 
      result = (index < validTriggerNames.size() && hHLTresults->accept(index)) ? 1 : 0;
      wasrun = (index < validTriggerNames.size() && hHLTresults->wasrun(index)) ? 1 : 0;

      iword = it/32;
      itrig = it%32;

      gHFEvent->fHLTNames[it] = TString(validTriggerNames[it]); 
      if (result) gHFEvent->fHLTWords[iword] |= (0x1 << itrig);
      if (wasrun) gHFEvent->fHLTWasRun[iword]|= (0x1 << itrig);
    }
  }
  

  // ----------------------------------------------------------------------
  // -- Get trigger objects
  // ----------------------------------------------------------------------

  Handle<trigger::TriggerEvent> trgEvent;
  hltF = true;
  try {
    iEvent.getByLabel(fTriggerEventLabel, trgEvent);
  } catch (const cms::Exception& e) {
    hltF = false;
    cout<<"Error!! No TriggerEvent with label " << fTriggerEventLabel << endl;
  }

  if (hltF) {
    TriggerObjectCollection allObjects = trgEvent->getObjects();
    for (int i=0; i < trgEvent->sizeFilters(); i++){         
      Keys keys = trgEvent->filterKeys(i);
      if (keys.size() > 0) {

	if (fVerbose > 1) cout << trgEvent->filterTag(i) << endl; 

	TString label = TString(trgEvent->filterTag(i).label())
	  + TString(":") 
	  + TString(trgEvent->filterTag(i).process())
	  + TString(":") 
	  + TString(trgEvent->filterTag(i).instance())
	  + TString(":");

	for (unsigned int j=0; j<keys.size(); j++){
	  TTrgObj *pTO = gHFEvent->addTrgObj();
	  pTO->fP.SetPtEtaPhiE(allObjects[keys[j]].pt(), 
			  allObjects[keys[j]].eta(), 
			  allObjects[keys[j]].phi(), 
			  allObjects[keys[j]].energy()
			  ); 
	  pTO->fID     = allObjects[keys[j]].id(); 
	  pTO->fLabel  = label;
	  // pTO->fNumber = i;
	  if (fVerbose > 1) 
	    cout << " pt = " <<  allObjects[keys[j]].pt() 
		 << " eta = " << allObjects[keys[j]].eta() 
		 << " phi = " << allObjects[keys[j]].phi() 
		 << " e = " << allObjects[keys[j]].energy() 
		 << " id = " << allObjects[keys[j]].id() 
		 << " label: " << pTO->fLabel
		 << endl;
	}
      }
    }
  }  
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTrigger::beginRun(const Run &run, const EventSetup &iSetup) {

}


// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTrigger::beginJob() {

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTrigger::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTrigger);
