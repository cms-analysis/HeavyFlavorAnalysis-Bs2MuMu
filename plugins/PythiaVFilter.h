/*
 *  PythiaVFilter.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.07.12.
 *
 */

#ifndef PYTHIAVFILTER_H
#define PYTHIAVFILTER_H

#include <set>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Specify a list of particles required to be available in the generated event
class PythiaVFilter : public edm::EDFilter {
	
	public:
		explicit PythiaVFilter(const edm::ParameterSet&);
		~PythiaVFilter();
		
		virtual bool filter(edm::Event&, const edm::EventSetup&);
	
	private:
		bool match(unsigned ix, std::vector<int> *ids, std::vector<double> *etas, std::vector<double> *pts, std::set<int> *used);
	
	private:
		int fVerbose;
		std::string mcLabel;
		std::vector<int>	particleIDs;
		std::vector<double>	minEtas;
		std::vector<double> maxEtas;
		std::vector<double> minPts;
		std::vector<double> maxPts;
		
		std::set<int> lookFor;
};

#endif
