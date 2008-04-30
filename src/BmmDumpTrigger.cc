#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpTrigger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

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
using namespace l1extra ;


// ----------------------------------------------------------------------
BmmDumpTrigger::BmmDumpTrigger(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fparticleMap(iConfig.getUntrackedParameter< std::string > ("particleMap")),
  fL1MuLabel(iConfig.getUntrackedParameter< std::string > ("L1MuLabel")),
  fL1TriggerName(iConfig.getUntrackedParameter<string>("L1TriggerName", string("L1_DoubleMu3"))),
  fHLTriggerLabel(iConfig.getParameter<edm::InputTag>("hltLabel")),
  fHLTriggerName(iConfig.getUntrackedParameter<string>("HLTriggerName",  string("HLTBJPsiMuMu"))),
  fHLTFilterObject0(iConfig.getUntrackedParameter<string>("HLTfiltObj0",  string("DiMuonNoIsoLevel1Seed"))),
  fHLTFilterObject1(iConfig.getUntrackedParameter<string>("HLTfiltObj1",  string("DiMuonNoIsoLevel1Filtered"))),
  fHLTFilterObject2(iConfig.getUntrackedParameter<string>("HLTfiltObj2",  string("JpsitoMumuL1Seed"))),
  fHLTFilterObject3(iConfig.getUntrackedParameter<string>("HLTfiltObj3",  string("JpsitoMumuL2Filtered"))),
  fHLTFilterObject4(iConfig.getUntrackedParameter<string>("HLTfiltObj4",  string("displacedJpsitoMumuFilter"))) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- BmmDumpTrigger constructor" << endl;
  cout << "--- Verbose          : " << fVerbose << endl;
  cout << "--- particle Map     : " << fparticleMap.c_str() << endl;
  cout << "--- L1 Mu Label      : " << fL1MuLabel.c_str() << endl;
  cout << "--- L1 Trigger Name  : " << fL1TriggerName.c_str() << endl;
  cout << "--- HLT Label        : " << fHLTriggerLabel.encode() << endl;
  cout << "--- HLT Trigger Name : " << fHLTriggerName.c_str() << endl;
  cout << "--- HLT Filter Objects : " 
       << fHLTFilterObject0.c_str() << ", "
       << fHLTFilterObject1.c_str() << ", "
       << fHLTFilterObject2.c_str() << ", "
       << fHLTFilterObject3.c_str() << " and "
       << fHLTFilterObject4.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0; 
  
}


// ----------------------------------------------------------------------
BmmDumpTrigger::~BmmDumpTrigger() {
  
}


// ----------------------------------------------------------------------
void BmmDumpTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  TAnaTrack *pTrack;

  // ===================================================================================
  // -- L1 trigger
  // ===================================================================================

  // -- adding L1 muons to signal block

  edm::Handle<l1extra::L1MuonParticleCollection> l1extmu;

  try {
    
    iEvent.getByLabel(fL1MuLabel,l1extmu);
  
  } catch (Exception event) { 
  
    if ( fVerbose ) cout << "  -- No L1Mu objects" << endl;
  }
  
  const l1extra::L1MuonParticleCollection& L1ExtMu = *l1extmu;
 
  if ( fVerbose ) cout << endl << endl << " -- L1MuonParticleCollection -- " << endl;

  if (&L1ExtMu) {

    l1extra::L1MuonParticleCollection myl1mus;
    int il1exmu = 0;

    for (l1extra::L1MuonParticleCollection::const_iterator muItr = L1ExtMu.begin(); muItr != L1ExtMu.end(); ++muItr) {
      
      L1MuGMTExtendedCand gmtCand = muItr->gmtMuonCand();

      if ( fVerbose ) {

	cout << "==> L1Muon " << il1exmu << endl;
	cout << "pt " << muItr->pt() << " e " << muItr->energy() << " eta " << muItr->eta() 
	     << " phi " << muItr->phi() << endl;
      
	cout << "iso " << muItr->isIsolated() << " mip " << muItr->isMip() 
	     << " forward " << muItr->isForward() << " RPC " << muItr->isRPC() << endl;

 	cout << " .... quality " << gmtCand.quality() << endl;
      }


      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 11;
      pTrack->fChi2     = gmtCand.quality();
      pTrack->fMuID     = il1exmu;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fMCID     = muItr->charge()*-13; 
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
    if ( fVerbose ) cout << "%BmmDumpTrigger -- No L1 MU object" << endl;
  }


  // -- Getting L1 trigger bits

  if ( fVerbose ) cout << endl << "============================ L1 Trigger ===================================" << endl; 

  gHFEvent->fL1Decision=0;
  gHFEvent->fL1w1=0; gHFEvent->fL1w2=0; gHFEvent->fL1w3=0; gHFEvent->fL1w4=0; 

  edm::Handle<l1extra::L1ParticleMapCollection> l1mapcoll; 
  try {

    iEvent.getByLabel(fparticleMap,l1mapcoll );

  } catch (Exception event) { 

    if ( fVerbose ) cout << "  -- No L1 Map Collection" << endl;
  }

  const l1extra::L1ParticleMapCollection& L1MapColl = *l1mapcoll; 

  if (&L1MapColl) {
 
    // Fill L1 trigger bits
    for (int itrig = 0; itrig != l1extra::L1ParticleMap::kNumOfL1TriggerTypes; ++itrig){
      
      const l1extra::L1ParticleMap& map = ( L1MapColl )[ itrig ] ;
      TString trigName = map.triggerName();
      int l1flag = map.triggerDecision();
      if ( fVerbose ) cout  << endl << "%L1-Report --  #"  << setw(3)<< itrig << " : "  
			    << setw(30) << trigName << " = " << l1flag;
   
      if (itrig < 32 && l1flag) {
	gHFEvent->fL1w1 |= (0x1 << itrig);
      } else if (itrig < 64 && l1flag) {
	gHFEvent->fL1w2 |= (0x1 << itrig);
      } else if (itrig < 96 && l1flag) {
	gHFEvent->fL1w3 |= (0x1 << itrig);
      } else if (itrig < 128 && l1flag) {
	gHFEvent->fL1w4 |= (0x1 << itrig);
      }

      if ( !strcmp(trigName, fL1TriggerName.c_str()) && l1flag ) {
	
	if ( fVerbose ) cout << " =============>  L1 decision set to " << l1flag;
	gHFEvent->fL1Decision = l1flag;      
      }

      // -- adding different L1 muons to signal block
      int muon_type(-1);
      if ( !strcmp(trigName, "L1_DoubleMu3" ) ) muon_type = 12;
      if ( !strcmp(trigName, "L1_SingleMu3" ) ) muon_type = 13;
      if ( !strcmp(trigName, "L1_SingleMu5" ) ) muon_type = 14;
      if ( !strcmp(trigName, "L1_SingleMu7" ) ) muon_type = 15;
      if ( !strcmp(trigName, "L1_SingleMu10") ) muon_type = 16;
      if ( !strcmp(trigName, "L1_SingleMu14") ) muon_type = 17;
      if ( !strcmp(trigName, "L1_SingleMu20") ) muon_type = 18;
      if ( !strcmp(trigName, "L1_SingleMu25") ) muon_type = 19;
 
     
      if( l1flag && (muon_type > 0) ) {

	const L1MuonParticleVectorRef& muons = map.muonParticles();
	
	if (fVerbose) cout << endl << " ==>  Muon particles pT, eta:" << endl;
	
	int l1cnt(0);
	for( L1MuonParticleVectorRef::const_iterator muonItr = muons.begin() ;
	     muonItr != muons.end() ; ++muonItr ) {
	  

	  if (fVerbose) cout << "\t pT = " << (*muonItr)->pt() << ", eta = " << (*muonItr)->eta() << endl ;
	  
	  
	  pTrack   = gHFEvent->addSigTrack();
	  pTrack->fMuType   = muon_type;
	  pTrack->fMCID     = (*muonItr)->charge()*-13; 
	  pTrack->fMuID     = l1cnt;
	  pTrack->fIndex    = -1;
	  pTrack->fGenIndex = -1; 
	  pTrack->fQ        = (*muonItr)->charge();
	  pTrack->fPlab.SetPtEtaPhi((*muonItr)->pt(),
				    (*muonItr)->eta(),
				    (*muonItr)->phi()
				    );
	  
	  l1cnt++;
	}
      }
    }
    
  } else {
    
    if ( fVerbose ) cout << "%BmmDumpTrigger -- No L1 Map Collection" << endl;
  }



  // ===================================================================================
  // -- HLT trigger
  // ===================================================================================
  
  gHFEvent->fHLTDecision=0;
  gHFEvent->fHLTw1=0; gHFEvent->fHLTw2=0; gHFEvent->fHLTw3=0; gHFEvent->fHLTw4=0; 

  if ( fVerbose ) cout << endl << "============================ HLT Trigger ===================================" << endl;
  edm::Handle<edm::TriggerResults> hltresults;
    
  try {
    iEvent.getByLabel(fHLTriggerLabel,hltresults);
    
  } catch (Exception event) {
    if ( fVerbose ) cout << "%BmmDumpTrigger -- Couldn't get handle on HLT Trigger!" << endl;
  }
   
  if (!hltresults.isValid()) {
    if ( fVerbose ) cout << "%BmmDumpTrigger -- No Trigger Result!" << endl;
  } 
  else {
    int ntrigs=hltresults->size();
    if (ntrigs==0){
	if ( fVerbose ) cout << "%BmmDumpTrigger -- No trigger name given in TriggerResults of the input " << endl;
    } 

    // get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    triggerNames_.init(*hltresults); 
    
    for (int itrig=0; itrig< ntrigs; itrig++) {

      TString trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      int wasrun  = (*hltresults)[itrig].wasrun();

      if ( fVerbose ) cout << endl << "%HLT-Report --  #"  << setw(3) << itrig << setw(30)  <<  trigName
			   << " was run " << wasrun  << ", decision " << hltflag;
         
      if (itrig < 32 && hltflag) {
	gHFEvent->fHLTw1 |= (0x1 << itrig);
      } else if (itrig < 64 && hltflag) {
	gHFEvent->fHLTw2 |= (0x1 << itrig);
      } else if (itrig < 96 && hltflag) {
	gHFEvent->fHLTw3 |= (0x1 << itrig);
      } else if (itrig < 128 && hltflag) {
	gHFEvent->fHLTw4 |= (0x1 << itrig);
      }

      if ( !strcmp(trigName, fHLTriggerName.c_str() ) && hltflag ) {     
	if ( fVerbose ) cout << " =============>  HLT decision set to " << hltflag;
	gHFEvent->fHLTDecision = hltflag;
      }
    } 
  }


  /////////////////////////////////////////////////
 
  // get hold of filtered candidates

  if ( fVerbose ) cout << endl << endl << " -- HLT Filter Objects -- " << endl;

  // -- filter object 0
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands0;
  try {

    iEvent.getByLabel(fHLTFilterObject0.c_str(), filtercands0); 
    if ( fVerbose ) cout << "  -- Found HLT " << fHLTFilterObject0.c_str() <<  " Collection, size " 
			 << filtercands0->size() << endl;

    int hltcnt(0);
    for (unsigned int i=0; i<filtercands0->size(); i++) {

      float pt  = filtercands0->getParticleRef(i).get()->pt();
      float phi = filtercands0->getParticleRef(i).get()->phi();
      float eta = filtercands0->getParticleRef(i).get()->eta();

      if ( fVerbose ) cout << endl << "==> " << fHLTFilterObject0.c_str() << setw(20) << "    Muon " 
			   << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 20;
      pTrack->fMCID     = filtercands0->getParticleRef(i).get()->charge()*-13; 
      pTrack->fMuID     = hltcnt;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands0->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

      hltcnt++;
     
    }

  } catch (Exception event) {
   
    if ( fVerbose ) cout << "  -- No HLT " << fHLTFilterObject0.c_str() <<  " Collection" << endl;
  }

  // -- filter object 1
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands1;
  try {

    iEvent.getByLabel(fHLTFilterObject1.c_str(), filtercands1); 
    if ( fVerbose ) cout << "  -- Found HLT " << fHLTFilterObject1.c_str() <<  " Collection, size " 
			 << filtercands1->size() << endl;

    int hltcnt(0);
    for (unsigned int i=0; i<filtercands1->size(); i++) {

      float pt  = filtercands1->getParticleRef(i).get()->pt();
      float phi = filtercands1->getParticleRef(i).get()->phi();
      float eta = filtercands1->getParticleRef(i).get()->eta();

      if ( fVerbose ) cout << endl << "==> " << fHLTFilterObject1.c_str() << setw(20) << "    Muon " 
			   << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 21;
      pTrack->fMCID     = filtercands1->getParticleRef(i).get()->charge()*-13; 
      pTrack->fMuID     = hltcnt;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands1->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

      hltcnt++;
     
    }

  } catch (Exception event) {
   
    if ( fVerbose ) cout << "  -- No HLT " << fHLTFilterObject1.c_str() <<  " Collection" << endl;
  }

  // -- filter object 2
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands2;
  try {

    iEvent.getByLabel(fHLTFilterObject2.c_str(), filtercands2);
    if ( fVerbose ) cout << "  -- Found HLT " << fHLTFilterObject2.c_str() <<  " Collection, size " 
			 << filtercands2->size() << endl;

    int hltcnt(0);
    for (unsigned int i=0; i<filtercands2->size(); i++) {

      float pt  = filtercands2->getParticleRef(i).get()->pt();
      float phi = filtercands2->getParticleRef(i).get()->phi();
      float eta = filtercands2->getParticleRef(i).get()->eta();

      if ( fVerbose ) cout << endl << "==> " << fHLTFilterObject2.c_str() << setw(20) << "    Muon " 
			   << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 
    
      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 22;
      pTrack->fMCID     = filtercands2->getParticleRef(i).get()->charge()*-13; 
      pTrack->fMuID     = hltcnt;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands2->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

      hltcnt++;

    }
  } catch (Exception event) {
   
    if ( fVerbose ) cout << "  -- No HLT " << fHLTFilterObject2.c_str() <<  " Collection" << endl;
  }

  // -- filter object 3
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands3;
  try {

    iEvent.getByLabel(fHLTFilterObject3.c_str(), filtercands3);
    if ( fVerbose ) cout << "  -- Found HLT " << fHLTFilterObject3.c_str() <<  " Collection, size " 
			 << filtercands3->size() << endl;

    int hltcnt(0);
    for (unsigned int i=0; i<filtercands3->size(); i++) {

      float pt  = filtercands3->getParticleRef(i).get()->pt();
      float phi = filtercands3->getParticleRef(i).get()->phi();
      float eta = filtercands3->getParticleRef(i).get()->eta();

      if ( fVerbose ) cout << endl << "==> " << fHLTFilterObject3.c_str() << setw(20) << "    Muon " 
			   << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 23;
      pTrack->fMCID     = filtercands3->getParticleRef(i).get()->charge()*-13; 
      pTrack->fMuID     = hltcnt;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands3->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

      hltcnt++;

    }
  } catch (Exception event) {

    if ( fVerbose ) cout << "  -- No HLT " << fHLTFilterObject3.c_str() <<  " Collection" << endl;
  }

  // -- filter object 4
  edm::Handle<reco::HLTFilterObjectWithRefs> filtercands4;
  try {

    iEvent.getByLabel(fHLTFilterObject4.c_str(), filtercands4);
    if ( fVerbose ) cout << "  -- Found HLT " << fHLTFilterObject4.c_str() <<  " Collection, size " 
			 << filtercands4->size() << endl;

    int hltcnt(0);
    for (unsigned int i=0; i<filtercands4->size(); i++) {
      float pt= filtercands4->getParticleRef(i).get()->pt();
      float phi= filtercands4->getParticleRef(i).get()->phi();
      float eta= filtercands4->getParticleRef(i).get()->eta();

      if ( fVerbose ) cout << endl << "==> " << fHLTFilterObject4.c_str() << setw(20) << "    Muon " 
			   << i << ": pt " << pt << " phi " << phi << " eta " << eta << endl; 

      pTrack   = gHFEvent->addSigTrack();
      pTrack->fMuType   = 24;
      pTrack->fMCID     = filtercands4->getParticleRef(i).get()->charge()*-13; 
      pTrack->fMuID     = hltcnt;
      pTrack->fIndex    = -1;
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = filtercands4->getParticleRef(i).get()->charge();

      pTrack->fPlab.SetPtEtaPhi(pt,
				eta,
				phi
				); 

      hltcnt++;
    }
  } catch (Exception event) {
   
    if ( fVerbose ) cout << "  -- No HLT " << fHLTFilterObject4.c_str() <<  " Collection " << endl;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void  BmmDumpTrigger::beginJob(const edm::EventSetup& setup) {

}

// ------------ method called once each job just after ending the event loop  ------------
void  BmmDumpTrigger::endJob() {

}

//define this as a plug-in
//DEFINE_FWK_MODULE(BmmDumpTrigger);
