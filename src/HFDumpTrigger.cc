#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpTrigger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

// -- Yikes!
extern TAna00Event *gHFEvent;

using namespace std;
using namespace edm;


// ----------------------------------------------------------------------
HFDumpTrigger::HFDumpTrigger(const edm::ParameterSet& iConfig):
  fTriggerLabel(iConfig.getParameter<edm::InputTag>("triggerLabel")),
  fL1MuLabel(iConfig.getUntrackedParameter< std::string > ("L1MuLabel")),
  fparticleMap(iConfig.getUntrackedParameter< std::string > ("particleMap")),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fL1TriggerName(iConfig.getUntrackedParameter<string>("L1TriggerName", string("L1_DoubleMu3"))),
  fHLTriggerName(iConfig.getUntrackedParameter<string>("triggerName",  string("HLTBJPsiMuMu"))) {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrigger constructor" << endl;
  cout << "--- Verbose: " << fVerbose << endl;
  cout << "--- " << fTriggerLabel.encode() << endl;
  cout << "--- " << fparticleMap.c_str() << endl;
  cout << "--- " << fL1MuLabel.c_str() << endl;
  cout << "--- " << fL1TriggerName.c_str() << endl;
  cout << "--- " << fHLTriggerName.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl; 
  
}


// ----------------------------------------------------------------------
HFDumpTrigger::~HFDumpTrigger() {
  
}


// ----------------------------------------------------------------------
void HFDumpTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  nevt++;

  // ===================================================================================
  // -- L1 trigger
  // ===================================================================================

  cout << endl << "============================ L1 Trigger ===================================" << endl; 

  gHFEvent->fL1Decision=0;
  gHFEvent->fL1w1=0; gHFEvent->fL1w2=0; gHFEvent->fL1w3=0; gHFEvent->fL1w4=0; 

  edm::Handle<l1extra::L1ParticleMapCollection> l1mapcoll; 
  try {iEvent.getByLabel(fparticleMap,l1mapcoll );} catch (Exception event) { cout << "  -- No L1 Map Collection" << endl;}
  const l1extra::L1ParticleMapCollection& L1MapColl = *l1mapcoll; 

  if (&L1MapColl) {
 
    // Fill L1 trigger bits
    for (int itrig = 0; itrig != l1extra::L1ParticleMap::kNumOfL1TriggerTypes; ++itrig){

      const l1extra::L1ParticleMap& map = ( L1MapColl )[ itrig ] ;
      TString trigName = map.triggerName();
      int l1flag = map.triggerDecision();

      cout << itrig << " : " << trigName << " :" << l1flag << endl;
   
      if (itrig < 32 && l1flag) {
	gHFEvent->fL1w1 |= (0x1 << itrig);
      } else if (itrig < 64 && l1flag) {
	gHFEvent->fL1w2 |= (0x1 << itrig);
      } else if (itrig < 96 && l1flag) {
	gHFEvent->fL1w3 |= (0x1 << itrig);
      } else if (itrig < 128 && l1flag) {
	gHFEvent->fL1w4 |= (0x1 << itrig);
      }
    }
   
   //      if ( !strcmp(trigName, fL1TriggerName1.c_str() ) ) {
   // 	cout << "!!! changing " << gHFEvent->fL1_mu3;
   //       	gHFEvent->fL1_mu3 = l1flag;
   // 	cout << " to " << gHFEvent->fL1_mu3 << endl;
   //       }
 
  } else {
    std::cout << "%HFDumpTrigger -- No L1 Map Collection" << std::endl;
  }

  edm::Handle<l1extra::L1MuonParticleCollection> l1extmu;
  try {iEvent.getByLabel(fL1MuLabel,l1extmu);} catch (Exception event) { cout << "  -- No L1Mu objects" << endl;}
  const l1extra::L1MuonParticleCollection& L1ExtMu = *l1extmu;
 
  TAnaTrack *pTrack;

  if (&L1ExtMu) {
    l1extra::L1MuonParticleCollection myl1mus;
    int il1exmu = 0;
    for (l1extra::L1MuonParticleCollection::const_iterator muItr = L1ExtMu.begin(); muItr != L1ExtMu.end(); ++muItr) {
      

      cout << "==> L1Muon " << il1exmu << endl;
      cout << "pt " << muItr->pt() << " e " << muItr->energy() << " eta " << muItr->eta() 
	   << " phi " << muItr->phi() << endl;
      
      cout << "iso " << muItr->isIsolated() << " mip " << muItr->isMip() 
	   << " forward " << muItr->isForward() << " RPC " << muItr->isRPC() << endl;
          
      L1MuGMTExtendedCand gmtCand = muItr->gmtMuonCand();

      cout << "quality " << gmtCand.quality() << endl;


      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 11;
      pTrack->fMCID     = muItr->charge()*-13; 
      pTrack->fMuID     = gmtCand.quality();
      pTrack->fIndex    = il1exmu;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = muItr->charge();
      pTrack->fPlab.SetPtEtaPhi(muItr->pt(),
				muItr->eta(),
				muItr->phi()
				);

//       pL1Muon->fIsIso     = muItr->isIsolated();
//       pL1Muon->fIsMIP     = muItr->isMip();
//       pL1Muon->fIsForward = muItr->isForward(); 
//       pL1Muon->fIsRPC     = muItr->isRPC();

      
      il1exmu++;
    }
  }
  else {
    std::cout << "%HFDumpTrigger -- No L1 MU object" << std::endl;
  }

  // ===================================================================================
  // -- HLT trigger
  // ===================================================================================
  
  gHFEvent->fHLTDecision=0;
  gHFEvent->fHLTw1=0; gHFEvent->fHLTw2=0; gHFEvent->fHLTw3=0; gHFEvent->fHLTw4=0; 

  cout << endl << "============================ HLT Trigger ===================================" << endl;
  edm::Handle<edm::TriggerResults> hltresults;
    
  try {
    iEvent.getByLabel(fTriggerLabel,hltresults);
    
  } catch (Exception event) {
    cout << "%HFDumpTrigger -- Couldn't get handle on HLT Trigger!" << endl;
  }
   
  if (!hltresults.isValid()) {
    cout << "%HFDumpTrigger -- No Trigger Result!" << endl;
  } 
  else {
    int ntrigs=hltresults->size();
    if (ntrigs==0){
	cout << "%HFDumpTrigger -- No trigger name given in TriggerResults of the input " << endl;
    } 

    // get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    triggerNames_.init(*hltresults); 
    
    for (int itrig=0; itrig< ntrigs; itrig++) {

      TString trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      int wasrun  = (*hltresults)[itrig].wasrun();

      cout << "%HLT-Report --  #" << itrig << "    "  << trigName
	   << " decision " << hltflag  << " was run " << wasrun << endl;
         

      if (itrig < 32 && hltflag) {
	gHFEvent->fHLTw1 |= (0x1 << itrig);
      } else if (itrig < 64 && hltflag) {
	gHFEvent->fHLTw2 |= (0x1 << itrig);
      } else if (itrig < 96 && hltflag) {
	gHFEvent->fHLTw3 |= (0x1 << itrig);
      } else if (itrig < 128 && hltflag) {
	gHFEvent->fHLTw4 |= (0x1 << itrig);
      }

//       if ( !strcmp(triggerNames_.triggerName(i).c_str(), fTriggerName.c_str() ) ) {
// 	cout << "!!! changing " << gHFEvent->fTriggerDecision;
//       	gHFEvent->fTriggerDecision = (*hltresults)[i].accept();
// 	cout << " to " << gHFEvent->fTriggerDecision << endl;
//       }
    } 
   
  }

  /////////////////////////////////////////////////
 
  // get hold of filtered candidates
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands0;
  try {
    iEvent.getByLabel ("SingleMuNoIsoLevel1Seed",filtercands0); 
    for (unsigned int i=0; i<filtercands0->size(); i++) {

      float pt= filtercands0->getParticleRef(i).get()->pt();
      float phi= filtercands0->getParticleRef(i).get()->phi();
      float eta= filtercands0->getParticleRef(i).get()->eta();
      cout << endl << "==>SingleMuNoIsoLevel1Seed:    Muon " << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      
      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 20;
      pTrack->fMCID     = filtercands0->getParticleRef(i).get()->charge()*-13; 
      pTrack->fIndex    = i;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands0->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 
     
    }
//     cout << gHFEvent->nHLTMuonsL1S() << " L1 SEED MUONS IN THE EVENT" << endl;
  } 
  catch (Exception event) {
    cout << "  -- No HLT SingleMuNoIsoLevel1Seed Collection" << endl;
  }

  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands1;
  try {
    iEvent.getByLabel ("SingleMuNoIsoL1Filtered",filtercands1);
    for (unsigned int i=0; i<filtercands1->size(); i++) {
      float pt= filtercands1->getParticleRef(i).get()->pt();
      float phi= filtercands1->getParticleRef(i).get()->phi();
      float eta= filtercands1->getParticleRef(i).get()->eta();
      cout << endl << "==>SingleMuNoIsoL1Filtered:   Muon " << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl;
    
      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 21;
      pTrack->fMCID     = filtercands0->getParticleRef(i).get()->charge()*-13; 
      pTrack->fIndex    = i;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands0->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

    }
//     cout << gHFEvent->nHLTMuonsL1F() << " L1 FILTERED MUONS IN THE EVENT" << endl;
  } 
  catch (Exception event) {
    cout << "  -- No HLT SingleMuNoIsoL1Filtered Collection" << endl;
  }

  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands2;
  try {
    iEvent.getByLabel ("SingleMuNoIsoL2PreFiltered",filtercands2);
    for (unsigned int i=0; i<filtercands2->size(); i++) {
      float pt= filtercands2->getParticleRef(i).get()->pt();
      float phi= filtercands2->getParticleRef(i).get()->phi();
      float eta= filtercands2->getParticleRef(i).get()->eta();
      cout << endl << "==>SingleMuNoIsoL2PreFiltered: Muon " << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 22;
      pTrack->fMCID     = filtercands0->getParticleRef(i).get()->charge()*-13; 
      pTrack->fIndex    = i;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands0->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

    }
//     cout << gHFEvent->nHLTMuonsL2F() << " L2 FILTERED MUONS IN THE EVENT" << endl;
  }
  catch (Exception event) {
    cout << "  -- No HLT SingleMuNoIsoL2PreFiltered Collection" << endl;
  }

  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands3;
  try {
    iEvent.getByLabel ("SingleMuNoIsoL3PreFiltered",filtercands3);
    for (unsigned int i=0; i<filtercands3->size(); i++) {
      float pt= filtercands3->getParticleRef(i).get()->pt();
      float phi= filtercands3->getParticleRef(i).get()->phi();
      float eta= filtercands3->getParticleRef(i).get()->eta();
      cout << endl << "==>SingleMuNoIsoL3PreFiltered: Muon " << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl;

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 23;
      pTrack->fMCID     = filtercands0->getParticleRef(i).get()->charge()*-13; 
      pTrack->fIndex    = i;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands0->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 
    }
//     cout << gHFEvent->nHLTMuonsL3F() << " L3 FILTERED MUONS IN THE EVENT" << endl;
  }
  catch (Exception event) {
    cout << "  -- No HLT SingleMuNoIsoL3PreFiltered Collection" << endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////

}

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTrigger::beginJob(const edm::EventSetup& setup) {

  nevt=0;
  
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTrigger::endJob() {

  cout << "HFDumpTrigger> Summary: Events processed: " << nevt << endl;
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFDumpTrigger);
