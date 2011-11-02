#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::make_pair;
using namespace reco;


// ----------------------------------------------------------------------
HFTwoParticleCombinatorics::HFTwoParticleCombinatorics(int verbose) : ttBuilder(NULL){  
  fVerbose  = verbose;
}

HFTwoParticleCombinatorics::HFTwoParticleCombinatorics(int verbose, const TransientTrackBuilder *fTTB ) :
  ttBuilder(fTTB) {  
  fVerbose  = verbose;
}

// ----------------------------------------------------------------------
HFTwoParticleCombinatorics::~HFTwoParticleCombinatorics() {
  if (fVerbose > 2) cout << "This is the end" << endl;
}

// ----------------------------------------------------------------------
void HFTwoParticleCombinatorics::combine(vector<pair<int, int> > &combList, 
					 vector<pair<int, TLorentzVector> > &tlist1, 
					 vector<pair<int, TLorentzVector> > &tlist2, 
					 double loMass, double hiMass, int rmDuplicate) {

  if (fVerbose > 3) {
    cout << "tlist1: " << tlist1.size() << endl;
    for (unsigned int i = 0; i < tlist1.size(); ++i) {
      cout << "i = " << i << " index = " << tlist1[i].first 
	   << " pT = " << tlist1[i].second.Pt() 
	   << " phi = " << tlist1[i].second.Phi() 
	   << " eta = " << tlist1[i].second.Eta() 
	   << endl;
    }
    cout << "tlist2: " << tlist2.size() << endl;
    for (unsigned int i = 0; i < tlist2.size(); ++i) {
      cout << "i = " << i << " index = " << tlist2[i].first 
	   << " pT = " << tlist2[i].second.Pt() 
	   << " phi = " << tlist2[i].second.Phi() 
	   << " eta = " << tlist2[i].second.Eta() 
	   << endl;
    }
  }

  TLorentzVector mom; 
  double mass(0.); 
  int duplicate(0); 
  for (unsigned int i = 0; i < tlist1.size(); ++i) {
    for (unsigned int j = 0; j < tlist2.size(); ++j) {
      if (tlist2[j].first == tlist1[i].first) continue;
      mom = tlist1[i].second +  tlist2[j].second;
      mass = mom.M(); 
      if (loMass < mass && mass < hiMass) {
	duplicate = 0; 
	for (unsigned int k = 0; k < combList.size(); ++k) {
	  if (combList[k].second == tlist1[i].first && combList[k].first == tlist2[j].first) { 
	    duplicate = 1; 
	    break;
	  }
	}
	if (1 == rmDuplicate) {
	  if (0 == duplicate) combList.push_back(make_pair(tlist1[i].first, tlist2[j].first));
	} else {
	   combList.push_back(make_pair(tlist1[i].first, tlist2[j].first));
	}
      }
    }
  }

  if (fVerbose > 2) {
    for (unsigned int k = 0; k < combList.size(); ++k) {
      cout << "combList. 1: " << combList[k].first << " 2: " << combList[k].second << endl;
    }
  }

}

// ----------------------------------------------------------------------
// This is a variant of the combinatorics using HFTwoParticleState as the return list
// Uses tracks (tf. needs masses). Suitable for use with V0 collection
void HFTwoParticleCombinatorics::combine(vector<HFTwoParticleState > &combList, 
					 vector<pair<int, reco::Track> > &tlist1, double mass1,
					 vector<pair<int, reco::Track> > &tlist2, double mass2,
					 double loMass, double hiMass, double maxDoca, int rmDuplicate) {
  if(ttBuilder==NULL) return;

  TwoTrackMinimumDistance md;
  TransientTrack tt1, tt2;
  FreeTrajectoryState freestates[2];
  TLorentzVector part1, part2, mom; 
  double mass(0.); 
  int duplicate(0); 

  for (unsigned int i = 0; i < tlist1.size(); ++i) {
    for (unsigned int j = 0; j < tlist2.size(); ++j) {
      if (tlist2[j].first == tlist1[i].first) continue; 
      if (tlist1[i].second.charge() == tlist2[j].second.charge()) continue;
      tt1=ttBuilder->build(tlist1[i].second);
      freestates[0]=tt1.impactPointTSCP().theState();
      tt2=ttBuilder->build(tlist2[j].second);
      freestates[1]=tt2.impactPointTSCP().theState();
      md.calculate(freestates[0], freestates[1]);

      if(!md.status()) continue;
      if(md.distance()<0.0) continue;
      if(md.distance()>maxDoca) continue;
      GlobalPoint vtx=md.crossingPoint();
      TrajectoryStateClosestToPoint p1TSCP = tt1.trajectoryStateClosestToPoint( vtx );
      TrajectoryStateClosestToPoint p2TSCP = tt2.trajectoryStateClosestToPoint( vtx );
      if( !p1TSCP.isValid() || !p2TSCP.isValid() ) continue;

      GlobalVector p1= p1TSCP.momentum();
      GlobalVector p2= p2TSCP.momentum();

      part1.SetXYZM(p1.x(), p1.y(), p1.z(), mass1);
      part2.SetXYZM(p2.x(), p2.y(), p2.z(), mass2);

      mom = part1 + part2;
      mass = mom.M();
      HFTwoParticleState p;
      if (loMass < mass && mass < hiMass) {
	duplicate = 0; 
	for (unsigned int k = 0; k < combList.size(); ++k) {
	  if (combList[k].id2 == tlist1[i].first && combList[k].id1 == tlist2[j].first) {
	    duplicate = 1; 
	    break;
	  }
	}
	if (1 == rmDuplicate) {  
	  if (0 == duplicate) {
	    p.id1=tlist1[i].first;
	    p.id2=tlist2[j].first;
	    p.mother = mom;
	    combList.push_back(p);
	  }
	} else {
	    p.id1=tlist1[i].first;
	    p.id2=tlist2[j].first;
	    p.mother = mom;
	    combList.push_back(p);
	}
      }
    }
  }

  if (fVerbose > 2) {
    for (unsigned int k = 0; k < combList.size(); ++k) {
      cout << "combList. 1: " << combList[k].id1 << " 2: " << combList[k].id2 << endl;
    }
  }

}
