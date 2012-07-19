/*
 *  PythiaVFilter.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.07.12.
 *
 */

#include "PythiaVFilter.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <iostream>

using std::set;
using std::vector;
using std::endl;
using std::cout;

template<class T>
void showVec(T v) {cout << '\t' << v;}

PythiaVFilter::PythiaVFilter(const edm::ParameterSet& iConfig) :
	fVerbose(iConfig.getUntrackedParameter("verbose",0)),
	mcLabel(iConfig.getUntrackedParameter("MCLabel",std::string("generator")))
{
	// default initializations
	double defMinPt = -9999, defMaxPt = 9999;
	double defMinEta = -9999, defMaxEta = 9999;
	vector<int> defParticleID;
	vector<double> defMinEtaV;
	vector<double> defMaxEtaV;
	vector<double> defMinPtV;
	vector<double> defMaxPtV;
	defParticleID.push_back(0);
	defMinEtaV.push_back(defMinEta);
	defMaxEtaV.push_back(defMaxEta);
	defMinPtV.push_back(defMinPt);
	defMaxPtV.push_back(defMaxPt);
	
	particleIDs = iConfig.getUntrackedParameter<vector<int> >("ParticleIDs",defParticleID);
	minEtas = iConfig.getUntrackedParameter<vector<double> >("MinEtas",defMinEtaV);
	maxEtas = iConfig.getUntrackedParameter<vector<double> >("MaxEtas",defMaxEtaV);
	minPts = iConfig.getUntrackedParameter<vector<double> >("MinPts",defMinPtV);
	maxPts = iConfig.getUntrackedParameter<vector<double> >("MaxPts",defMaxPtV);
	
	// sanity check
	while (minEtas.size() < particleIDs.size()) minEtas.push_back(defMinEta);
	while (maxEtas.size() < particleIDs.size()) maxEtas.push_back(defMaxEta);
	while (minPts.size() < particleIDs.size()) minPts.push_back(defMinPt);
	while (maxPts.size() < particleIDs.size()) maxPts.push_back(defMaxPt);
	
	// construct the particle lookup Table
	for (vector<int>::const_iterator it = particleIDs.begin(); it != particleIDs.end(); ++it)
		lookFor.insert(*it);
	
	// show the configuration
	cout << "-----------------------------" << endl;
	cout << "--- PythiaVFilter" << endl;
	cout << "verbose: " << fVerbose << endl;
	cout << "particleIDs:"; std::for_each(particleIDs.begin(),particleIDs.end(),showVec<int>); cout << endl;
	cout << "minEtas:"; std::for_each(minEtas.begin(),minEtas.end(),showVec<double>); cout << endl;
	cout << "maxEtas:"; std::for_each(maxEtas.begin(),maxEtas.end(),showVec<double>); cout << endl;
	cout << "minPts:"; std::for_each(minPts.begin(),minPts.end(),showVec<double>); cout << endl;
	cout << "maxPts:"; std::for_each(maxPts.begin(),maxPts.end(),showVec<double>); cout << endl;
} // PythiaVFilter()


PythiaVFilter::~PythiaVFilter()
{} // ~PythiaVFilter()

bool PythiaVFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	Handle<HepMCProduct> evt;
	iEvent.getByLabel(mcLabel,evt);
	HepMC::GenEvent genEvent(*(evt->GetEvent()));
	vector<double> currentEtas;
	vector<double> currentPts;
	vector<int> currentIDs;
	set<int> used;
	bool accepted = false;
	
	if (fVerbose > 0) {
		cout << "Looking for particles:";
		std::for_each(particleIDs.begin(),particleIDs.end(),showVec<int>);
		cout << endl;
	}
	
	for (HepMC::GenEvent::particle_iterator p = genEvent.particles_begin(); p != genEvent.particles_end(); ++p) {
		
		if(lookFor.count( (*p)->pdg_id() ) == 0)
			continue; // not interested in that guy
		
		if (fVerbose > 5)
			cout << "Found particle " << (*p)->pdg_id() << " with pt = " << (*p)->momentum().perp() << " and eta = " << (*p)->momentum().eta() << endl;
		
		// save this particles and its interessted properties
		currentIDs.push_back( (*p)->pdg_id() );
		currentPts.push_back( (*p)->momentum().perp() );
		currentEtas.push_back( (*p)->momentum().eta() );
	}
	
	if (fVerbose > 5) {
		cout << "current IDs: "; std::for_each(currentIDs.begin(), currentIDs.end(), showVec<int>); cout	<< endl;
		cout << "current Pts: "; std::for_each(currentPts.begin(), currentPts.end(), showVec<double>); cout << endl;
		cout << "current Etas: "; std::for_each(currentEtas.begin(), currentEtas.end(), showVec<double>); cout << endl;
	}
	
	// try to match
	accepted = match(0,&currentIDs,&currentEtas,&currentPts,&used);
	
	if (fVerbose > 5)
		cout << "accepted = " << accepted << endl;
	
	return accepted;
} // filter()

bool PythiaVFilter::match(unsigned ix, vector<int> *ids, vector<double> *etas, vector<double> *pts, set<int> *used)
{
	vector<int>::iterator it = ids->begin();
	bool matched = false;
	int v;
	
	if (ix >= particleIDs.size()) {
		matched = true;
		goto bail;
	}
	
	if (fVerbose > 10)
		cout << "Matching particle " << particleIDs[ix] << " with pt in [" << minPts[ix] << ", " << maxPts[ix] << "] and eta in [" << minEtas[ix] << ", " << maxEtas[ix] << "]" << endl;
	
	// try to find a candidate for ix
	while( !matched && (it = std::find(it, ids->end(), particleIDs[ix])) != ids->end() ) {
		
		v = it - ids->begin(); // index
		if (fVerbose > 10)
			cout << "Found particle " <<  (*ids)[v] << " pt = " << (*pts)[v] << ", eta = " << (*etas)[v] << endl;
		
		if ( minEtas[ix] < (*etas)[v] && (*etas)[v] < maxEtas[ix]) {
			// eta matches
			if ( minPts[ix] < (*pts)[v] && (*pts)[v] < maxPts[ix] ) {
				// pt matches
				if (used->count(v) == 0) {
					// not yet taken
					used->insert(v); // take it
					matched = match(ix+1, ids, etas, pts, used);
					used->erase(v); // and remove it
				}
			}
		}
		++it; // jump over current one
	}
	
bail:
	return matched;
} // match()

//define this as a plug-in
DEFINE_FWK_MODULE(PythiaVFilter);
