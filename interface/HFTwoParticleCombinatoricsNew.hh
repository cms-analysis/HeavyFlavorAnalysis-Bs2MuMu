/*
 *  HFTwoParticleCombinatoricsNew.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.11.12.
 *
 */

#ifndef HFTWOPARTICLECOMBINATORICSNEW_H
#define HFTWOPARTICLECOMBINATORICSNEW_H

#include <set>
#include <utility>

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/View.h"

class HFTwoParticleCombinatoricsCompare
{
	public:
		HFTwoParticleCombinatoricsCompare(int rmDuplicate = 0) : fRmDuplicate(rmDuplicate) {}
		
		bool operator()(const std::pair<int,int> &p1, const std::pair<int,int> &p2) {
			double v1,v2;
			bool result = false;
			if (fRmDuplicate) {
				v1 = std::min(p1.first,p1.second);
				v2 = std::min(p2.first,p2.second);
				if (v1 < v2) {
					result = true;
				} else if (v2 < v1) {
					result = false;
				} else {
					v1 = std::max(p1.first,p1.second);
					v2 = std::max(p2.first,p2.second);
					if (v1 < v2) {
						result = true;
					} else {
						result = false;
					}
				}
			} else {
				result = p1 < p2;
			}
			
			return result;
		}
	private:
		int fRmDuplicate;
};

class HFTwoParticleCombinatoricsNew {
	
	public:
		typedef std::set<std::pair<int,int>, HFTwoParticleCombinatoricsCompare> HFTwoParticleCombinatoricsSet;
		typedef HFTwoParticleCombinatoricsSet::iterator iterator;
	
	public:
		HFTwoParticleCombinatoricsNew(edm::Handle<edm::View<reco::Track> > &hTracks, int verbose) : fhTracks(hTracks), fVerbose(verbose) {}
		
		HFTwoParticleCombinatoricsSet combine(std::vector<int> &tlist1, double mass1, std::vector<int> &tlist2, double mass2, double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);
	
	private:
		edm::Handle<edm::View<reco::Track> > &fhTracks;
		int fVerbose;
};

#endif
