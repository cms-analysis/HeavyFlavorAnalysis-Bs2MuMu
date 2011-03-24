#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpGenerator.h"

#include <iostream>

#include <TRandom.h>

#include "HepMC/GenVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HFDumpGenerator::HFDumpGenerator(const ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fGenCandidatesLabel(iConfig.getUntrackedParameter<string>("generatorCandidates", string("MCCandidate"))), 
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEvent", string("EvtGenProducer")))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpGenerator constructor: " << fGenCandidatesLabel << "  " << fGenEventLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpGenerator::~HFDumpGenerator() {
  
}


// ----------------------------------------------------------------------
void HFDumpGenerator::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- CAREFUL!
  gHFEvent->Clear();

  // Das Clear ist dafuer da, wenn ein GenFilter (der mit gHFEvent arbeitet
  // anstatt mit dem HEPMC EDM-product direkt) da ist, der nicht alle
  // Events durchlaesst. Wenn das Clear() nicht da waere, wuerde bei einem
  // neuen Event der alte (weggefilterte) GenBlock auch noch drin sein.


  static int nevt(0); 
  ++nevt;

  if (fVerbose > 3) {
    cout << "=================HEPMC===================" << endl;
    Handle<HepMCProduct> evt;
    iEvent.getByLabel(fGenEventLabel.c_str(), evt);
    const HepMC::GenEvent *genEvent = evt->GetEvent();
    genEvent->print();
    
    cout << "=================HEPMC===================" << endl;
  }

  TGenCand  *pGen;
  // -- From PhysicsTools/HepMCCandAlgos/plugins/ParticleListDrawer.cc
  int iMo1(-1), iMo2(-1), iDa1(-1), iDa2(-1); 

  std::vector<const GenParticle *> cands;
  cands.clear();
  vector<const GenParticle *>::const_iterator found = cands.begin();

  edm::Handle<GenParticleCollection> genParticlesH;
  genParticlesH.clear();
  try {
    // iEvent.getByLabel ("genParticles", genParticlesH);
    iEvent.getByLabel (fGenCandidatesLabel.c_str(), genParticlesH);
  } catch(cms::Exception ce) {
    cout << "==> HFDumpGenerator caught std::exception " << ce.what() << endl;
  }

  for (GenParticleCollection::const_iterator p = genParticlesH->begin(); p != genParticlesH->end(); ++p) {
    cands.push_back( & * p );
  }

  if (fVerbose > 1) cout << Form("Number of genParticles = %i", (int)genParticlesH->size()) << endl;

  int i(-1); 
  for(GenParticleCollection::const_iterator p  = genParticlesH->begin(); p != genParticlesH->end();  p++) {
    ++i; 
    pGen = gHFEvent->addGenCand();
    pGen->fID     = p->pdgId();
    pGen->fStatus = p->status();  
    pGen->fNumber = i; 
    pGen->fQ      = p->charge(); 

    double vx = p->vx(), vy = p->vy(), vz = p->vz();
    pGen->fP.SetXYZT(p->px(), 
		     p->py(), 
		     p->pz(), 
		     p->energy());
    pGen->fV.SetXYZ(vx, vy, vz);

    // Particles Mothers and Daighters
    iMo1 = -1;
    iMo2 = -1;
    iDa1 = -1;
    iDa2 = -1;
    int nMo = p->numberOfMothers();
    int nDa = p->numberOfDaughters();

    pGen->fDau1 = 99999;
    pGen->fDau2 = -1;

    found = find(cands.begin(), cands.end(), p->mother(0));
    if (found != cands.end()) {
      iMo1 = found - cands.begin();
      pGen->fMom1 = iMo1;
    } else {
      pGen->fMom1 = -1;
    }

    found = find(cands.begin(), cands.end(), p->mother(nMo-1));
    if (found != cands.end()) {
      iMo2 = found - cands.begin();
      pGen->fMom2 = iMo2;
    } else {
      pGen->fMom2 = -1;
    }

    found = find(cands.begin(), cands.end(), p->daughter(0));
    if (found != cands.end()) {
      iDa1 = found - cands.begin();
      pGen->fDau1 = iDa1;
    }

    found = find(cands.begin(), cands.end(), p->daughter(nDa-1));
    if (found != cands.end()) {
      iDa2 = found - cands.begin();
      pGen->fDau2 = iDa2;
    }

    if (fVerbose > 2) pGen->dump();
  }
 
  genParticlesH.clear(); // WHY?

  if (fVerbose > 0) cout << "==> HFDumpGenerator> Event " << nevt << ", dumped  " << gHFEvent->nGenCands() << " generator cands" << endl;
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDumpGenerator::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpGenerator::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpGenerator);
