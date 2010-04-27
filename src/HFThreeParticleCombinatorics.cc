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
// For D+ -> K pi pi
// First is Kaon, second and third are pions with ordered indices
// ----------------------------------------------------------------------
void HFThreeParticleCombinatorics::combine(vector<triplet> &combList, 
					 vector<pair<int, TLorentzVector> > &kaList, 
					 vector<pair<int, TLorentzVector> > &piList, 
					 double loMass, double hiMass, int rmDuplicate) {

  if (fVerbose > 2) {
    cout << "kaList: " << kaList.size() << endl;
    for (unsigned int i = 0; i < kaList.size(); ++i) {
      cout << "i = " << i << " index = " << kaList[i].first 
	   << " pT = " << kaList[i].second.Pt() 
	   << " phi = " << kaList[i].second.Phi() 
	   << " eta = " << kaList[i].second.Eta() 
	   << endl;
    }
    cout << "piList: " << piList.size() << endl;
    for (unsigned int i = 0; i < piList.size(); ++i) {
      cout << "i = " << i << " index = " << piList[i].first 
	   << " pT = " << piList[i].second.Pt() 
	   << " phi = " << piList[i].second.Phi() 
	   << " eta = " << piList[i].second.Eta() 
	   << endl;
    }
  }
  
  TLorentzVector mom; 
  double mass(0.); 
  int duplicate(0); 
  vector<pair<int, TLorentzVector> >::iterator ka, pi1, pi2;
  for (ka = kaList.begin(); ka != kaList.end(); ka++) {                 // Kaon
    for (pi1 = piList.begin(); pi1 != piList.end(); pi1++) {            // first pion
      if (pi1->first == ka->first) continue;                            // don't use the same track twice
      for(pi2 = pi1+1; pi2 != piList.end(); pi2++) {                    // second pion
         if (pi2->first == ka->first) continue;                         // don't use the same track twice
	 triplet t=triplet(ka->first, pi1->first, pi2->first);          // triplet with three different tracks
         mom = ka->second +  pi1->second + pi2->second;                 // mother is D+ candidate
         mass = mom.M(); 
         cout << "K pi pi mass: " << ka->first<<" / "<<pi1->first << " / " << pi2->first<< "  = "<<mass<<endl;
         if (loMass < mass && mass < hiMass) {                          // continue only when in mass window
	   duplicate = 0;                                               // check for duplication 
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

  if (fVerbose > 1) {
      for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	cout << "K pi pi List: " << it->ka() << " / " << it->pi1() << " / " << it->pi2()<< endl;
    }
  }

}




// ----------------------------------------------------------------------
// For Kshort pi --> 3 pi
// First two pions are from Kshort candidate, third comes from D+ vertex
// ----------------------------------------------------------------------
void HFThreeParticleCombinatorics::combine(vector<triplet> &combList, 
					 vector<pair<int, TLorentzVector> > &tlist, 
					 double loMass, double hiMass, int rmDuplicate) {

  if (fVerbose > 2) {
    cout << "tlist: " << tlist.size() << endl;
    for (unsigned int i = 0; i < tlist.size(); ++i) {
      cout << "i = " << i << " index = " << tlist[i].first 
	   << " pT = " << tlist[i].second.Pt() 
	   << " phi = " << tlist[i].second.Phi() 
	   << " eta = " << tlist[i].second.Eta() 
	   << endl;
    }
  }

  TLorentzVector mom; 
  double mass(0.); 
  int duplicate(0); 
  vector<pair<int, TLorentzVector> >::iterator i,j,k;
  for (i = tlist.begin(); i != tlist.end(); ++i) {                      // first pion from Kshort
    for (j = i+1; j != tlist.end(); ++j) {                              // second pion from Kshort
      mom = i->second +  j->second;                                     // Kshort candidate
      mass = mom.M(); 
      cout << "pi pi mass: " << i->first << " / " << j->first<< "  = "<<mass<<endl;
      if (loMass > mass || mass > hiMass) continue;                     // continue only if in mass window
      for(k = tlist.begin(); k!=tlist.end(); k++) {                     // third pion, from D+
	if (k==i || k==j) continue;                                     // don't use same track twice
	triplet t=triplet(k->first, i->first, j->first);                // triplet with three different tracks
	duplicate = 0;                                                  // check for duplication                                      
	for(vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	    if(t==*it){
	       duplicate = 1; 
	       break;
	    }
	}
      	if(1 == rmDuplicate) {
	   if (0 == duplicate) combList.push_back(t);                   // add only if not already in the list
	} else {
	   combList.push_back(t);                                       // add anyway 
	}
      }
    }
  }

  if (fVerbose > 0) {
      for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	cout << "Triple Pion List: " << it->pi1() << " /  " << it->pi2() << " / " << it->pi3()<< endl;
    }
  }

}
