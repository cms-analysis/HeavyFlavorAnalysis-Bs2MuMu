/*
 *  HFGenDimuon.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 18.07.12.
 *
 */

#ifndef HFGENDIMUON_H
#define HFGENDIMUON_H

#include <set>
#include <map>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HFGenDimuon : public edm::EDAnalyzer {
	
	public:
		explicit HFGenDimuon(const edm::ParameterSet&);
		~HFGenDimuon();
	
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
	
		void combine(std::multiset<int>::iterator b, std::multiset<int>::iterator e, std::map<int,int> *part, std::set<int> *indices);
	
	private:
		int fVerbose;
		int fGenType;
		std::multiset<int> fDaughtersID;
};

#endif
