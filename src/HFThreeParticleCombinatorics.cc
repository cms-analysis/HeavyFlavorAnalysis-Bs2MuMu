#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFThreeParticleCombinatorics.hh"


// ----------------------------------------------------------------------
HFThreeParticleCombinatorics::HFThreeParticleCombinatorics(int verbose) {  
  fVerbose  = verbose;
}

// ----------------------------------------------------------------------
HFThreeParticleCombinatorics::~HFThreeParticleCombinatorics() {
  if (fVerbose > 2) cout << "This is the end" << endl;
}

// ----------------------------------------------------------------------
void HFThreeParticleCombinatorics::combine(vector<triplet> &combList, 
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
  for (unsigned int i = 0; i < tlist1.size(); ++i) {                    // Kaon
    for (unsigned int j = 0; j < tlist2.size(); ++j) {                  // first pion
      if (tlist2[j].first == tlist1[i].first) continue;
      for(unsigned int k = j+1; k<tlist2.size(); k++) {                 // second pion
         if (tlist2[k].first == tlist1[i].first) continue;
	 triplet t=triplet(tlist1[i].first, tlist2[j].first, tlist2[k].first);        // triplet with three different tracks
         mom = tlist1[i].second +  tlist2[j].second + tlist2[k].second;
         mass = mom.M(); 
         if (loMass < mass && mass < hiMass) {
	    duplicate = 0; 
	    for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	       if(t==*it){
	          duplicate = 1; 
	          break;
	       }
	    }
      	    if (1 == rmDuplicate) {
	      if (0 == duplicate) combList.push_back(t);                   // add only if not already in the list
	    } else {
	      combList.push_back(t);                                       // add anyway 
	    }
         } 
      }
    }
  }

  if (fVerbose > 0) {
      for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
      cout << "combList. 1: " << it->ka << " 2: " << it->pi1 << " 3: " << it->pi2<< endl;
    }
  }

}


