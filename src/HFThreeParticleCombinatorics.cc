#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFThreeParticleCombinatorics.hh"

using std::cout;
using std::endl;
using std::vector;
using std::pair;


// ----------------------------------------------------------------------
HFThreeParticleCombinatorics::HFThreeParticleCombinatorics(int verbose) {  
  fVerbose  = verbose;
}

// ----------------------------------------------------------------------
HFThreeParticleCombinatorics::~HFThreeParticleCombinatorics() {
  if (fVerbose > 2) cout << "This is the end" << endl;
}


// ----------------------------------------------------------------------
// Combining particle from kaList with two particles from piList into a combination
// with mass between loMass and hiMass
//
// E.g.: D+ -> K pi pi
// First is Kaon, second and third are pions with ordered indices
// ----------------------------------------------------------------------
void HFThreeParticleCombinatorics::combine(vector<triplet> &combList, 
					 vector<pair<int, TLorentzVector> > &kaList, 
					 vector<pair<int, TLorentzVector> > &piList, 
					 double loMass, double hiMass) {

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
         if (loMass < mass && mass < hiMass) {                          // continue only when in mass window
	   duplicate = 0;                                               // check for duplication 
	    for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	       if(t==*it){
	          duplicate = 1; 
	          break;
	       }
	    }
	    if (0 == duplicate) combList.push_back(t);                   // add only if not already in the list
         } 
      }
    }
  }

  if (fVerbose > 0) {
      for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	cout << "K pi pi List: " << it->ka() << " / " << it->pi1() << " / " << it->pi2()<< endl;
    }
  }

}




// ----------------------------------------------------------------------
// Combine three particles from piList into a combination with mass between loMass and hiMass
// if hiResMass >0.0, the first two particles form a resonance
// with mass between loResMass and hiResMass
//
// E.g.: Kshort pi --> 3 pi
// First two pions are from Kshort candidate, third comes from D+ vertex
// ----------------------------------------------------------------------
void HFThreeParticleCombinatorics::combine(vector<triplet> &combList, 
					 vector<pair<int, TLorentzVector> > &piList, 
					   double loMass, double hiMass, double loResMass, double hiResMass) {

  if (fVerbose > 2) {
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
  vector<pair<int, TLorentzVector> >::iterator i,j,k;
  for (i = piList.begin(); i != piList.end(); ++i) {                         // first pion from Kshort
    for (j = i+1; j != piList.end(); ++j) {                                  // second pion from Kshort
      mom = i->second +  j->second;                                          // Resonance candidate
      mass = mom.M(); 
      if(hiResMass>0. && (loResMass > mass || mass > hiResMass)) continue;   // ignore if resonance outside mass window
      for(k = piList.begin(); k!=piList.end(); k++) {                        // third pion, from D+
	if (k==i || k==j) continue;                                          // don't use same track twice
	mom = i->second + j->second + k->second;                             // three particle candidate
        mass = mom.M();
        if(loMass > mass || mass > hiMass) continue;                        // ignore if outside mass window
	triplet t=triplet(k->first, i->first, j->first);                     // triplet with three different tracks
	duplicate = 0;                                                       // check for duplication                 
	for(vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	    if(t==*it){
	       duplicate = 1; 
	       break;
	    }
	}
	if (0 == duplicate) combList.push_back(t);                   // add only if not already in the list
      }
    }
  }

  if (fVerbose > 0) {
      for (vector<triplet>::iterator it=combList.begin(); it != combList.end(); it++) {
	cout << "Triple Pion List: " << it->pi1() << " /  " << it->pi2() << " / " << it->pi3()<< endl;
    }
  }

}
