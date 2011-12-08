#ifndef DANEKUTILS_H
#define DANEKUTILS_H

#include <iostream>
#include <utility>

#define DR      57.29577951

class danekUtils {

public:
  danekUtils();
  ~danekUtils();

  static inline double twoBodyDecayAngle(TVector3 v1, TVector3 v2) {
    float deta = v1.Eta() - v2.Eta();
    float dphi = v1.Phi() - v2.Phi();
    float dr = sqrt(deta*deta+dphi*dphi);
    return dr;    
    //return v1.Angle(v2);
  };

  static inline double twoBodyDecayMomPerp(TVector3 v1, TVector3 v2) {
    TVector3 pSum = v1 + v2;
    return pSum.Perp();
  };

  static double twoBodyDecayMass(TVector3 v1, TVector3 v2, double m1, double m2) {

    TVector3 pSum = v1 + v2;
    
    //double angle = v1.Angle(v2);      
    float p1 = v1.Mag();
    float p2 = v2.Mag();
    float e1 = sqrt(p1*p1+m1*m1);
    float e2 = sqrt(p2*p2+m2*m2);
    float eSum = e1 + e2;
    //float m = sqrt( eSum*eSum - pSum.Mag2() );
    
    TLorentzVector lv(pSum, eSum); 
    float m = lv.M();
    return m;
    
  };

};

#endif
