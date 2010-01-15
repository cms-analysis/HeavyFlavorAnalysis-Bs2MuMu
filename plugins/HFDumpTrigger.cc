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

///////My includes ////
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
  fL1GTReadoutRecordLabel(iConfig.getUntrackedParameter<InputTag>("L1GTReadoutRecordLabel", edm::InputTag("gtDigis"))),
  fL1GTmapLabel(iConfig.getUntrackedParameter<InputTag>("hltL1GtObjectMap")),
  fL1MuonsLabel(iConfig.getUntrackedParameter<InputTag>("L1MuonsLabel")),
  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel"))
{

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrigger constructor" << endl;
  cout << "--- Verbose                     : " << fVerbose << endl;
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

  // -- L1 trigger bits
  gHFEvent->fL1Decision = 0; 
  gHFEvent->fL1TWords[0]=0; gHFEvent->fL1TWords[1]=0; gHFEvent->fL1TWords[2]=0; gHFEvent->fL1TWords[3]=0; 
  gHFEvent->fL1TWasRun[0]=0;gHFEvent->fL1TWasRun[1]=0;gHFEvent->fL1TWasRun[2]=0;gHFEvent->fL1TWasRun[3]=0; 
  gHFEvent->fL1TTWords[0]=0;gHFEvent->fL1TTWords[1]=0;

  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(fL1GTReadoutRecordLabel,L1GTRR);
  Handle<L1GlobalTriggerObjectMapRecord> hL1GTmap; 
  iEvent.getByLabel("hltL1GtObjectMap", hL1GTmap);

  edm::ESHandle<L1GtTriggerMenu> hL1GtMenu;
  iSetup.get<L1GtTriggerMenuRcd>().get(hL1GtMenu);
  const L1GtTriggerMenu* l1GtMenu = hL1GtMenu.product();


  if (L1GTRR.isValid()) {
    gHFEvent->fL1Decision = (L1GTRR->decision()? 1: 0);

    const AlgorithmMap& algorithmMap = l1GtMenu->gtAlgorithmMap();
    for (CItAlgo itAlgo = algorithmMap.begin(); itAlgo != algorithmMap.end(); itAlgo++) {
      std::string aName = itAlgo->first;
      int algBitNumber = (itAlgo->second).algoBitNumber();
      if (fVerbose > 2) cout << "i = " << algBitNumber << " -> " << aName << endl;
      gHFEvent->fL1TNames[algBitNumber] = TString(aName);
    }

    const AlgorithmMap& algorithmTTMap = l1GtMenu->gtTechnicalTriggerMap();
    for (CItAlgo itAlgo = algorithmTTMap.begin(); itAlgo != algorithmTTMap.end(); itAlgo++) {
      std::string aName = itAlgo->first;
      int algBitNumber = (itAlgo->second).algoBitNumber();
      if (fVerbose > 2) cout << "i = " << algBitNumber << " -> " << aName << endl;
      gHFEvent->fL1TTNames[algBitNumber] = TString(aName);
    }

    int itrig(0);
    for (unsigned int iTrig = 0; iTrig < L1GTRR->decisionWord().size(); ++iTrig) {
      int l1flag = L1GTRR->decisionWord()[iTrig]; 
      int t1flag = L1GTRR->technicalTriggerWord()[iTrig]; 
      itrig = iTrig%32;
      if (iTrig < 32) {
	if (l1flag) gHFEvent->fL1TWords[0]  |= (0x1 << itrig);//note: changin iTrig to itrig = iTrig%32
	if (t1flag) gHFEvent->fL1TTWords[0] |= (0x1 << itrig);
      } else if (iTrig < 64) {
	if (l1flag) gHFEvent->fL1TWords[1]  |= (0x1 << itrig);//note: changin iTrig to itrig = iTrig%32
	if (t1flag) gHFEvent->fL1TTWords[1] |= (0x1 << itrig);
      } else if (iTrig < 96) {
	if (l1flag) gHFEvent->fL1TWords[2]  |= (0x1 << itrig);//note: changin iTrig to itrig = iTrig%32
      } else if (iTrig < 128 && l1flag) {
	if (l1flag) gHFEvent->fL1TWords[3]  |= (0x1 << itrig);//note: changin iTrig to itrig = iTrig%32
      }
    }

  } 

  if (fVerbose > 1) {
    cout << "L1 trigger: " << endl
	 << std::bitset<32>(gHFEvent->fL1TWords[0]) << endl
	 << std::bitset<32>(gHFEvent->fL1TWords[1]) << endl
	 << std::bitset<32>(gHFEvent->fL1TWords[2]) << endl
	 << std::bitset<32>(gHFEvent->fL1TWords[3]) << endl;
    cout << "L1 technical trigger: " << endl
	 << std::bitset<32>(gHFEvent->fL1TTWords[0]) << endl
	 << std::bitset<32>(gHFEvent->fL1TTWords[1]) << endl;
  }

  // -- L1 muons
  Handle<l1extra::L1MuonParticleCollection> hL1Muons;
  try {
    iEvent.getByLabel(fL1MuonsLabel, hL1Muons);
    l1extra::L1MuonParticleCollection::const_iterator l1muon;
    for (l1muon = hL1Muons->begin(); l1muon != hL1Muons->end(); ++l1muon ) {
      TLorentzVector a(0.0,0.0,0.0,0.0);
      a.SetPxPyPzE(l1muon->px(),l1muon->py(),l1muon->pz(),l1muon->energy());
      if (fVerbose > 3) cout << "L1 muon: pT = " << a.Pt() << " phi = " << a.Phi() << " eta = " << a.Eta() << " q = " << l1muon->charge() << endl;
    }
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFDumpTrigger>  l1extra::L1MuonParticleCollection " << fL1MuonsLabel.encode() << " not found " << endl;
  }


  // -- HLT bits
  gHFEvent->fHLTDecision = 0; 
  gHFEvent->fHLTWords[0]=0; gHFEvent->fHLTWords[1]=0; gHFEvent->fHLTWords[2]=0; gHFEvent->fHLTWords[3]=0; 
  gHFEvent->fHLTWasRun[0]=0;gHFEvent->fHLTWasRun[1]=0;gHFEvent->fHLTWasRun[2]=0;gHFEvent->fHLTWasRun[3]=0; 

  Handle<TriggerResults> hHLTresults;
  bool hltF = true;
  try {
    iEvent.getByLabel(fHLTResultsLabel, hHLTresults);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFDumpTrigger> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    hltF = false;
  }

  if (hltF) {
    TriggerNames trigName;
    trigName.init(*hHLTresults);
    //    HLTGlobal_wasrun=HLTR->wasrun();
    //    HLTGlobal_error=HLTR->error();
    gHFEvent->fHLTDecision = hHLTresults->accept();
    if (fVerbose > 1) cout << "hHLTresults->size() = " << hHLTresults->size() << " and HLT accept = " << gHFEvent->fHLTDecision << endl;

    for (unsigned int iTrig = 0; iTrig < hHLTresults->size(); ++iTrig) {
      int hltflag = hHLTresults->accept(iTrig); 
      int wasrun = hHLTresults->wasrun(iTrig); 

      // This is the HLT name for each path: 
      if (fVerbose > 2) cout << iTrig << " " << trigName.triggerName(iTrig) << endl;

      gHFEvent->fHLTNames[iTrig] = TString(trigName.triggerName(iTrig)); 

      if (iTrig < 32 && hltflag) {
	gHFEvent->fHLTWords[0] |= (0x1 << iTrig);
      } else if (iTrig < 64 && hltflag) {
	gHFEvent->fHLTWords[1] |= (0x1 << iTrig);
      } else if (iTrig < 96 && hltflag) {
	gHFEvent->fHLTWords[2] |= (0x1 << iTrig);
      } else if (iTrig < 128 && hltflag) {
	gHFEvent->fHLTWords[3] |= (0x1 << iTrig);
      }

      if (iTrig < 32 && wasrun) {
	gHFEvent->fHLTWasRun[0] |= (0x1 << iTrig);
      } else if (iTrig < 64 && hltflag) {
	gHFEvent->fHLTWasRun[1] |= (0x1 << iTrig);
      } else if (iTrig < 96 && hltflag) {
	gHFEvent->fHLTWasRun[2] |= (0x1 << iTrig);
      } else if (iTrig < 128 && hltflag) {
	gHFEvent->fHLTWasRun[3] |= (0x1 << iTrig);
      }
    }


    if (fVerbose > 1)  cout << "HLT trigger accept/run: " << endl
			    << std::bitset<32>(gHFEvent->fHLTWords[0]) << endl
			    << std::bitset<32>(gHFEvent->fHLTWasRun[0]) << endl << endl
			    << std::bitset<32>(gHFEvent->fHLTWords[1]) << endl
			    << std::bitset<32>(gHFEvent->fHLTWasRun[1]) << endl<< endl
			    << std::bitset<32>(gHFEvent->fHLTWords[2]) << endl
			    << std::bitset<32>(gHFEvent->fHLTWasRun[2]) << endl << endl
			    << std::bitset<32>(gHFEvent->fHLTWords[3]) << endl
			    << std::bitset<32>(gHFEvent->fHLTWasRun[3]) << endl;
  }   
  
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
void  HFDumpTrigger::beginJob(const edm::EventSetup& setup) {

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTrigger::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTrigger);
