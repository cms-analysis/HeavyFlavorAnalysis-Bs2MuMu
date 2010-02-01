#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"


// ----------------------------------------------------------------------
HFTwoParticleCombinatorics::HFTwoParticleCombinatorics(int verbose) {  
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

  if (fVerbose > 2) {
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

  if (fVerbose > 0) {
    for (unsigned int k = 0; k < combList.size(); ++k) {
      cout << "combList. 1: " << combList[k].first << " 2: " << combList[k].second << endl;
    }
  }

}


