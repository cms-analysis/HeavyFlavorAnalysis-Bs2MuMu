/*
 *  HFGenDimuon.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 18.07.12.
 *
 */

#include "HFGenDimuon.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <iostream>
#include <utility>

using std::cout;
using std::endl;
using std::set;
using std::multiset;
using std::map;

extern TAna01Event *gHFEvent;

HFGenDimuon::HFGenDimuon(const edm::ParameterSet& iConfig) :
	fVerbose(iConfig.getUntrackedParameter<int>("Verbose",0)),
	fGenType(iConfig.getUntrackedParameter<int>("GenType",0))
{
	std::vector<int> defaultIDs(2,13);
	std::vector<int> dauVec = iConfig.getUntrackedParameter<std::vector<int> >("DaughtersID",defaultIDs);
	
	// convert to multiset
	for (std::vector<int>::const_iterator vit = dauVec.begin(); vit != dauVec.end(); ++vit)
		fDaughtersID.insert(*vit);
	
	// Dump configuration
	cout << "--- HFGenDimuon constructor" << endl;
	cout << "--- Verbose:		" << fVerbose << endl;
	cout << "--- GenType:		" << fGenType << endl;
	cout << "--- Daughters:	";
	for(std::multiset<int>::const_iterator sit = fDaughtersID.begin(); sit != fDaughtersID.end(); ++sit)
		cout << '\t' <<  *sit;
	cout << endl;
} // HFGenDimuon()	

HFGenDimuon::~HFGenDimuon()
{} // ~HFGenDimuon

void HFGenDimuon::beginJob()
{
} // beginJob()

void HFGenDimuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	// -- do all combinations with the daughters
	map<int,int> particles;
	set<int> indices;
	TGenCand *gen;
	TAnaVertex *vtx;
	int ix;
	
	// add the simulation Primary Vertex = Production of first daughter of proton
	for (ix = 0; ix < gHFEvent->nGenCands(); ix++) {
		gen = gHFEvent->getGenCand(ix);
		if (gen->fMom1 == 0) {
			vtx = gHFEvent->addPV();
			vtx->fPoint = gen->fV;
			break;
		}
	}
	
	for (ix = 0; ix < gHFEvent->nGenCands(); ix++) {
		gen = gHFEvent->getGenCand(ix);
		if (fDaughtersID.count(gen->fID) > 0)
			particles.insert(std::make_pair(ix,gen->fID));
	}
	
	if (fVerbose > 6) {
		cout << "==> HFGenDimuon: Found the following important particles in generator block" << endl;
		for (map<int,int>::const_iterator it = particles.begin(); it != particles.end(); ++it)
			cout << "\t(" << it->first << ", " << it->second << ")";
		cout << endl;
	}
	
	// try to find all daughters...
	combine(fDaughtersID.begin(),fDaughtersID.end(),&particles,&indices);
} // analyze()

void HFGenDimuon::endJob()
{
} // endJob()

void HFGenDimuon::combine(std::multiset<int>::iterator b, std::multiset<int>::iterator e, std::map<int,int> *part, std::set<int> *indices)
{
	TAnaCand *cand;
	TGenCand *gen;
	TAnaTrack *trk;
	set<int>::const_iterator it;
	multiset<int>::iterator old;
	map<int,int>::const_iterator mapIt;
	TLorentzVector p;
	set<int> bIx;
	
	if (b == e) {
		// found a valid combination
		if (fVerbose > 5) {
			cout << "==> HFGenDimuon: Found candidates with indices:" << endl;
			for (it = indices->begin(); it != indices->end(); ++it)
				cout << '\t' << *it;
			cout << endl;
		}
		
		cand = gHFEvent->addCand();
		cand->fType = fGenType;
		cand->fQ = 0;
		cand->fSig1 = gHFEvent->nSigTracks();
		cand->fSig2 = cand->fSig1 - 1;
		
		p.SetXYZT(0,0,0,0); // initialize to zero
		for (it = indices->begin(); it != indices->end(); ++it) {
			gen = gHFEvent->getGenCand(*it);
			cand->fQ = cand->fQ + gen->fQ;
			p = p + gen->fP;
			
			cand->fSig2 = gHFEvent->nSigTracks(); // enlarge the sig tracks array
			trk = gHFEvent->addSigTrack();
			trk->fIndex = -1; // no rec track as MC only
			trk->fGenIndex = *it; // generator index
			trk->fMCID = cand->fType;
			trk->fPlab = gen->fP.Vect();
			trk->fQ = gen->fQ;
		}
		
		cand->fMass = p.M();
		cand->fPlab = p.Vect();
		
		if (fVerbose > 6) {
			cout << "===> HFGenDimuon: Add candidate with:" << endl;
			cout << "	GenType:	" << cand->fType << endl;
			cout << "	Q:			" << cand->fQ << endl;
			cout << "	Mass:		" << cand->fMass << endl;
			cout << "	P			(" << cand->fPlab.X() << "," << cand->fPlab.Y() << "," << cand->fPlab.Z() << ")" << endl;
		}
		
		if (fVerbose > 7) {
			cout << "===> HFGenDimuon: Daughter production vertices:" << endl;
			for (it = indices->begin(); it != indices->end(); ++it) {
				gen = gHFEvent->getGenCand(*it);
				cout << '\t' << gen->fID << ":	(" << gen->fV.X() << "," << gen->fV.Y() << "," << gen->fV.Z() << ")" << endl;
			}
		}
	}
	else {
		// search for particle 'b' in part
		for (mapIt = part->begin(); mapIt != part->end(); ++mapIt) {
			if (mapIt->second == *b)
				bIx.insert(mapIt->first);
		}
		
		// process all of them
		old = b++;
		for (it = bIx.begin(); it != bIx.end(); ++it) {
			indices->insert(*it);
			part->erase(*it);
			combine(b,e,part,indices);
			part->insert(std::make_pair(*it,*old));
			indices->erase(*it);
		}
	}
} // combine()

//define this as a plug-in
DEFINE_FWK_MODULE(HFGenDimuon);
