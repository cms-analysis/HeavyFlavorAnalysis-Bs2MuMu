#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTrigger.h"

#include <iostream>
#include <bitset>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

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
  for (int i = 0; i < NL1T; ++i) {
    gHFEvent->fL1TPrescale[i] = gHFEvent->fL1TResult[i] = gHFEvent->fL1TMask[i] = gHFEvent->fL1TError[i] = 0; 
  }

  for (int i = 0; i < NLTT; ++i) {
    gHFEvent->fLTTPrescale[i] = gHFEvent->fLTTResult[i] = gHFEvent->fLTTMask[i] = gHFEvent->fLTTError[i] = 0; 
  }

  for (int i = 0; i < NHLT; ++i) {
    gHFEvent->fHLTPrescale[i] = gHFEvent->fHLTResult[i] = gHFEvent->fHLTWasRun[i] = gHFEvent->fHLTError[i] = 0; 
  }
    
  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(fL1GTReadoutRecordLabel,L1GTRR);
  Handle<L1GlobalTriggerObjectMapRecord> hL1GTmap; 
  iEvent.getByLabel("hltL1GtObjectMap", hL1GTmap);

  L1GtUtils l1GtUtils;
  l1GtUtils.retrieveL1EventSetup(iSetup);
  // cout << "L1 trigger menu: ";
  // cout << l1GtUtils.l1TriggerMenu() << endl;

  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  string algoname; 
  int    algobit(-1); 
  bool   result(false); 
  int    prescale(0); 
  int    mask(0); 
  int    iErrorCode(0); 

  for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);

    gHFEvent->fL1TNames[algobit]    = TString(algoname);
    gHFEvent->fL1TResult[algobit]   = result;
    gHFEvent->fL1TMask[algobit]     = mask;
    gHFEvent->fL1TError[algobit]    = iErrorCode;
    gHFEvent->fL1TPrescale[algobit] = prescale;
  }

  for (CItAlgo algo = menu->gtTechnicalTriggerMap().begin(); algo != menu->gtTechnicalTriggerMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);
    
    gHFEvent->fLTTNames[algobit]    = TString(algoname);
    gHFEvent->fLTTResult[algobit]   = result;
    gHFEvent->fLTTMask[algobit]     = mask;
    gHFEvent->fLTTError[algobit]    = iErrorCode;
    gHFEvent->fLTTPrescale[algobit] = prescale;
  }

  // ----------------------------------------------------------------------
  // -- HLT results
  // ----------------------------------------------------------------------

  // -- Read HLT configuration and names
  HLTConfigProvider hltConfig;
  bool hltConfigInitSuccess = hltConfig.init(fHLTProcessName);

  vector<string> validTriggerNames;
  if (hltConfigInitSuccess) validTriggerNames = hltConfig.triggerNames();
  //can assert?!  hltConfig.dump("PrescaleTable");

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
    const TriggerNames &trigName = iEvent.triggerNames(*hHLTresults);
    gHFEvent->fHLTDecision = hHLTresults->accept();
    if (fVerbose > 1) cout << "hHLTresults->size() = " << hHLTresults->size() << " and HLT accept = " << gHFEvent->fHLTDecision << endl;

    unsigned int index(999); 
    bool wasrun(false), result(false), error(false);
    int prescale(1); 
    int psSet = -1; //hltConfig.prescaleSet(iEvent, iSetup);
    for (unsigned int it = 0; it < validTriggerNames.size(); ++it) {
      index    = trigName.triggerIndex(validTriggerNames[it]); 
      result   = (index < validTriggerNames.size() && hHLTresults->accept(index));
      wasrun   = (index < validTriggerNames.size() && hHLTresults->wasrun(index));
      error    = (index < validTriggerNames.size() && hHLTresults->error(index));
      if (psSet > -1) {
	prescale = hltConfig.prescaleValue(psSet, validTriggerNames[it]);
      } else {
	//	cout << "==>HFDumpTrigger> error in prescale set!?" << endl;
	prescale = 0;
      }

      gHFEvent->fHLTNames[index]    = TString(validTriggerNames[it]);
      gHFEvent->fHLTResult[index]   = result;
      gHFEvent->fHLTWasRun[index]   = wasrun;
      gHFEvent->fHLTError[index]    = error;
      gHFEvent->fHLTPrescale[index] = prescale;
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
