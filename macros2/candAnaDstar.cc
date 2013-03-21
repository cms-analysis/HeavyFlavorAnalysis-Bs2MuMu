#include "candAnaDstar.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"
#include "danekUtils.h"

using namespace std;

namespace {
  TVector3 DSVertex(0,0,0), DZVertex(0,0,0), PV(0,0,0);  
  TVector3 DSMom(0,0,0), DZMom(0,0,0), PiSlowMom(0,0,0), PiMom(0,0,0), KMom(0,0,0);  
}

// ----------------------------------------------------------------------
candAnaDstar::candAnaDstar(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaDstar: name = " << name << ", reading cutsfile " << cutsFile << endl;

  readCuts(cutsFile, 1); 

}


// ----------------------------------------------------------------------
candAnaDstar::~candAnaDstar() {
  cout << "==> candAnaDstar: destructor..." << endl;
  tree->Write();
}


// ----------------------------------------------------------------------
// To analyze the MC event
bool candAnaDstar::anaMC(TAna01Event *evt) {
  const bool print = false;

  fpEvt = evt; 
  bool foundDs = false, foundDz = false, foundPiSlow = false, foundPi=false, foundK=false, foundPV=false;
  //TVector3 DSVertex(0,0,0), DSMom(0,0,0), DZVertex(0,0,0), DZMom(0,0,0), PV(0,0,0);  
  int qds =0, qk=0, qpi=0, qpis=0;
  int pC0 = 0;

  
  int numGenCands = fpEvt->nGenCands();
  if(print) cout << "Found " << numGenCands << " gen cands in event" << endl;
  for (int it = 0; it < numGenCands; ++it) {  // loop over all gen candidates 
    foundDs = false, foundDz = false, foundPiSlow = false, foundPi=false, foundK=false;

    TGenCand * pCand = fpEvt->getGenCand(it);

    if( !foundPV && (abs(pCand->fID) == 5) ) {foundPV=true; PV=pCand->fV; if(print) cout<<" PV "<<PV.Z()<<endl;}   // get PV 

    if( ( abs(pCand->fID) != 413) ) continue;        // skip others 

    if(print) cout <<" DS "<<it<< " " << pCand->fNumber << " "<<pCand->fID<<" "<<pCand->fQ<<" "<<pCand->fStatus<<" "
		   <<pCand->fMom1<<" "<<pCand->fMom2<<" "<<pCand->fDau1<<" "
		   <<pCand->fDau2<<" "<<pCand->fP.Perp()<<" "<<pCand->fV.Z()<<endl;
    qds = pCand->fQ;
    DSVertex = (pCand->fV);
    DSMom = (pCand->fP.Vect());
    foundDs = true;
    pC0 = it;

    int i1 = (pCand->fDau2)-(pCand->fDau1)+1;
    if(i1!=2) {continue;} // fpEvt->dumpGenBlock();}
    //if(i1!=2) {cout<<" number of daughters1 "<<i1<<endl;continue;} // fpEvt->dumpGenBlock();}
    for(int id=(pCand->fDau1);id<=(pCand->fDau2);++id) {
      TGenCand * dau = fpEvt->getGenCand(id);  // check daughters

      if( abs(dau->fID) == 421 ) { //  D0

	foundDz=true;
	if(print) cout <<" D0 "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;
	//TVector3 v1 = dau->fP.Vect();
	DZMom = (dau->fP.Vect());

       
	int i2 = (dau->fDau2)-(dau->fDau1)+1;
	if(i2!=2) {continue;} // fpEvt->dumpGenBlock();}
	//if(i2!=2) {cout<<" number of daughters2 "<<i2<<endl;continue;} // fpEvt->dumpGenBlock();}
	for(int igd=(dau->fDau1);igd<=(dau->fDau2);++igd) {
	  TGenCand * gdau = fpEvt->getGenCand(igd);  // check grand-daughters
	  //TVector3 v2 = gdau->fP.Vect();
          
	  if( abs(gdau->fID) == 321) {  // kaon  
	    foundK = true;
	    KMom = (gdau->fP.Vect());
	    //float pt = v2.Perp();
	    //float eta = v2.Eta();
	    //float phi = v2.Phi();
	    if(print) cout <<" K "<<gdau->fNumber << " "<<gdau->fID<<" "<<gdau->fQ<<" "
		 <<gdau->fMom1<<" "<<gdau->fMom2<<" "
		 <<gdau->fP.Perp()<<" "<<gdau->fV.Z()<<endl;
	    
	    DZVertex = (gdau->fV);
	    qk = gdau->fQ;

	  } else if( abs(gdau->fID) == 211) {  // pion  
	    foundPi = true;
	    PiMom = (gdau->fP.Vect());
	    qpi = gdau->fQ;
	    //float pt = v2.Perp();
	    //float eta = v2.Eta();
	    //float phi = v2.Phi();
	    if(print) cout <<" Pi "<<gdau->fNumber << " "<<gdau->fID<<" "<<gdau->fQ<<" "
			   <<gdau->fMom1<<" "<<gdau->fMom2<<" "
			   <<gdau->fP.Perp()<<" "<<gdau->fV.Z()<<endl;

	  } // 

	  if( foundPi && foundK) break;

	} // end granddaughter loop 

      } else if( (abs(dau->fID)) == 211) { // slow pion

	foundPiSlow = true;
	if(print) cout <<" Slow pi "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;
	PiSlowMom = dau->fP.Vect();
	qpis = dau->fQ;

      }  // 

      if(foundPiSlow && foundDz) break; 

    } // daugther loop 

    // exit if we find the right candidate
    if(foundPV && foundDs && foundDz && foundPiSlow && foundPi && foundK) break;

  } // gen part loop 


  bool ok = foundPV && foundDs && foundDz && foundPiSlow && foundPi && foundK;
  if(ok) {
    
    if( (qds != qpis) || (qds != -qk) || ( qpis != -qk) || (qpis != qpi) ) {
      cout<<pC0<<" wrong charge Ds,K,pi,pi_slow ";
      cout<<qds<<" "<<qk<<" "<<qpi<<" "<<qpis<<endl;
      //fpEvt->dumpGenBlock();
    }
    int tmp = qpi + qpis;
    ((TH1D*)fHistDir->Get("h23"))->Fill(float(tmp));

    TVector3 t1(DSVertex-PV), t2(DZVertex-PV), t3(DZVertex-DSVertex);
    double a1 = t1.Angle(DSMom);  // D* pointing angle
    double a2 = t2.Angle(DZMom);  // D0 pointing angle
    double a3 = t3.Angle(DZMom);  // D0 pointing angle with respect PV
    double a4 = t1.Angle(t3);     // SV1 versus SV2

    ((TH1D*)fHistDir->Get("h11"))->Fill(t1.Mag());
    ((TH1D*)fHistDir->Get("h12"))->Fill(t2.Mag());
    ((TH1D*)fHistDir->Get("h13"))->Fill(t3.Mag());
    ((TH1D*)fHistDir->Get("h20"))->Fill(a1);
    ((TH1D*)fHistDir->Get("h17"))->Fill(a2);
    ((TH1D*)fHistDir->Get("h18"))->Fill(a3);
    ((TH1D*)fHistDir->Get("h19"))->Fill(a4);

    ((TH1D*)fHistDir->Get("h14"))->Fill(PiSlowMom.Perp());
    ((TH1D*)fHistDir->Get("h15"))->Fill(PiSlowMom.Angle(DZMom));
    ((TH1D*)fHistDir->Get("h16"))->Fill(PiSlowMom.Angle(DSMom));

    ((TH1D*)fHistDir->Get("h21"))->Fill(DZMom.Angle(DSMom));
    //((TH1D*)fHistDir->Get("h23"))->Fill(t3.Angle(DSMom));
    ((TH1D*)fHistDir->Get("h25"))->Fill(t2.Angle(DSMom));

    ((TH1D*)fHistDir->Get("h28"))->Fill(DSMom.Perp());
    ((TH1D*)fHistDir->Get("h29"))->Fill(DZMom.Perp());
    ((TH1D*)fHistDir->Get("h30"))->Fill(KMom.Perp());
    ((TH1D*)fHistDir->Get("h30"))->Fill(PiMom.Perp());

    double angle = danekUtils::twoBodyDecayAngle(KMom, PiMom);
    double pt    = danekUtils::twoBodyDecayMomPerp(KMom, PiMom);
    double m1    = danekUtils::twoBodyDecayMass(KMom, PiMom, MKAON, MPION);
    double m2    = danekUtils::twoBodyDecayMass(KMom, PiMom, MPION, MKAON);
    TVector3 t4(KMom+PiMom);

    ((TH1D*)fHistDir->Get("h41"))->Fill(pt); 
    ((TH1D*)fHistDir->Get("h43"))->Fill(t4.Perp()); 

    ((TH1D*)fHistDir->Get("h42"))->Fill(m2); 
    ((TH1D*)fHistDir->Get("h44"))->Fill(m1); 

    tmp = DSMom.Perp();
    ((TH2D*)fHistDir->Get("h51"))->Fill(tmp,m2); 
    ((TH2D*)fHistDir->Get("h52"))->Fill(angle,m2); 
    ((TH2D*)fHistDir->Get("h53"))->Fill(pt,m2); 
    tmp = KMom.Perp();
    ((TH2D*)fHistDir->Get("h54"))->Fill(tmp,m2); 
    tmp = PiMom.Perp();
    ((TH2D*)fHistDir->Get("h55"))->Fill(tmp,m2); 
    tmp = PiSlowMom.Perp();
    ((TH2D*)fHistDir->Get("h56"))->Fill(tmp,m2); 
    tmp = PiSlowMom.Angle(DZMom);
    ((TH2D*)fHistDir->Get("h57"))->Fill(tmp,m2); 
      

    double m21    = danekUtils::twoBodyDecayMass(PiSlowMom, DZMom, MPION, m1);
    double m22    = danekUtils::twoBodyDecayMass(PiSlowMom, DZMom, MPION, m2);
    double angle2 = danekUtils::twoBodyDecayAngle(DZMom, PiSlowMom);
    double pt2    = danekUtils::twoBodyDecayMomPerp(DZMom, PiSlowMom);

    ((TH1D*)fHistDir->Get("h45"))->Fill(m21); 
    ((TH1D*)fHistDir->Get("h46"))->Fill(m22); 
    ((TH1D*)fHistDir->Get("h47"))->Fill(angle2); 
    ((TH1D*)fHistDir->Get("h48"))->Fill(pt2); 
    ((TH1D*)fHistDir->Get("h49"))->Fill(m21-m1); 
    ((TH1D*)fHistDir->Get("h50"))->Fill(m22-m2); 


  }
  return ok;
}

//---------------------------------------------------------------------------------------------------------------
// Main candidate processing for Dstar
void candAnaDstar::candAnalysis() {
  static int count0=0, count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  const bool useMyTrees = false;

  fPreselection = false;  //  reset

  if (0 == fpCand) return;

  dumpHFTruthCand(fpCand); // testing


  //cout<<" Call canAna::candAnalysis "<<endl;
  candAna::candAnalysis();  // call the main analysis 
  //cout<<" Call canAnaDstar::candAnalysis "<<endl;

  // now do Dstar specific

  count0++;

  ((TH1D*)fHistDir->Get("Status"))->Fill(0.);
  //cout<<" Dstar candidate "<<fpCand->fType<<" in event "<<fEvt<<" "<<fVerbose<<endl;

  if (fVerbose>10) {
    cout << "DUMP HFDstarCandidate with mass = " << fpCand->fMass <<" in event "<<fEvt<<" cand num "<<count0<<endl;
    dumpHFDstarCand(fpCand); 
    //dumpHFTruthCand(pCand); 
  }


  TAnaCand *pC(0);
  TAnaTrack *pK, *pPi, *pPis; 
  double fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), dr(0.); 
  //int indexPi=-1, indexK=-1;
  // -- D0 
  if (fpCand->fDau1 < 0  || fpCand->fDau2 < 0) {
    if(fVerbose>1) {cout << "pCand->fDauX = -1!!! " << fpCand->fType << " skip event "<<endl; fpCand->dump();}    
    return;  // skip if no daughters 
  }
  ((TH1D*)fHistDir->Get("Status"))->Fill(1.);

  pC = fpEvt->getCand(fpCand->fDau1);  // D0 candidate 
  pK = 0;
  pPi =0; 

  //int pPisId = fpEvt->getSigTrack(fpCand->fSig1)->fIndex; // slow pion index
  //pPis = fpEvt->getRecTrack(pPisId); // slow pi from D*
  int pPisId = fpCand->fSig1; // slow pion index
  pPis = fpEvt->getSigTrack(pPisId); // slow pi from D*
  TVector3 piSlowMom = pPis->fPlab; // slow pi momentum vector 

  if (fVerbose>10 ) cout<<" slow pi "<<pPisId<<" "<< pPis->fPlab.Perp()  << " "<<pPis->fMCID<<endl;
  if (fVerbose>10 ) cout<<" D0 "<<pC->fType<<endl;

  // loop over D0 tracks 
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    //int index = fpEvt->getSigTrack(id)->fIndex;
    int index = id;
    // 	if (fVerbose>0 && tm) cout<<id<<" "<<index<<" ";

    if (211 == fpEvt->getSigTrack(id)->fMCID) {  // pion 
      //pPi = fpEvt->getRecTrack(index);
      pPi = fpEvt->getSigTrack(index);
      //indexPi = index;
      if (fVerbose>10) cout<< pPi->fMCID<<" "<<pPi->fPlab.Perp() << endl;
    } else {  // kaon
      //pK = fpEvt->getRecTrack(index);
      pK = fpEvt->getSigTrack(index);
      //indexK = index;
      if (fVerbose>10 ) cout<< pK->fMCID<<" "<<pK->fPlab.Perp() << endl;
    }
  }
	
  if(pPi == 0 || pK==0) {
    if(fVerbose>0) {cout << " pi or K not found " << fpCand->fType << endl; fpCand->dump();}    
    return;  // skip if no daughters 
  }
  ((TH1D*)fHistDir->Get("Status"))->Fill(2.);


//       tm = 0; 
//       // truthMatch return an interger, 0-no match, 1-correct match, -1-pi&K swapped
//       if(fIsMC) tm = truthMatch(pCand,fVerbose); // check truth matching
//       if ( tm == 1 ) { // do only for MC
// 	ok = anaMC(evt);
// 	if (fVerbose>0) {
// 	  cout << " -> " << pCand->fType;
// 	  cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
// 	  cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
// 	  cout << " mass: " << pCand->fMass << " " << ok<<endl;
	
// 	  cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
// 	  dumpHFDstarCand(pCand); 
// 	  fpEvt->dumpGenBlock(); 
// 	}
//       } // end if

  // Pt, eta
  double pt    = fpCand->fPlab.Perp(); // cand pt
  double ptdz  = pC->fPlab.Perp();
  double ptPis = pPis->fPlab.Perp();
  double ptPi  = pPi->fPlab.Perp();
  double ptK   = pK->fPlab.Perp();

  double eta   = fpCand->fPlab.Eta();
  //double etaPis= pPis->fPlab.Eta();
  double etaPi = pPi->fPlab.Eta();
  double etaK  = pK->fPlab.Eta();

  // charge 
  int qk   = pK->fQ;
  int qpi  = pPi->fQ;
  int qpis = pPis->fQ;

  // masses
  double mdstar = fpCand->fMass; 
  double mdz = pC->fMass;
  double dm = mdstar - mdz;

  // D* vertex
  TAnaVertex sv = fpCand->fVtx;  // D* vertex
  TVector3 sv1 = sv.fPoint;  // D* decay vertex position 
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0);  // D* PV index
  TVector3 pv =  fpEvt->getPV(pvidx)->fPoint;  // Dstar vertex 

  // D0 vertex
  TAnaVertex svD0 = pC->fVtx;  // D0 vertex
  TVector3 sv2 = svD0.fPoint; // D0 decay vertex position
  //int pvidx2 = (pC->fPvIdx > -1? pC->fPvIdx : 0);  // D0 PV index, DOES NOT HAVE A PV
  //TVector3 pv2 =  fpEvt->getPV(pvidx2)->fPoint;    // It is always -1 

  //cout<< fpCand->fPvIdx <<" "<< pC->fPvIdx <<endl;

  if(fpCand->fPvIdx==-1 ) ((TH1D*)fHistDir->Get("Status"))->Fill(3.);
 

  //((TH1D*)fHistDir->Get("Status"))->Fill(4.);

  ((TH1D*)fHistDir->Get("pvidx"))->Fill(float(pvidx));

  // Calculate angles
  TVector3 t1(sv1-pv), t2(sv2-pv), t3(sv2-sv1);
  //cout<<t1.Z()<<" "<<t2.Z()<<" "<<t3.Z()<<endl;

      
  //fls3d = sv.fD3d/sv.fD3dE; // D*
  //flsxy = sv.fDxy/sv.fDxyE; 
  fls3d = svD0.fD3d/svD0.fD3dE; //  use D0
  flsxy = svD0.fDxy/svD0.fDxyE; 
  prob  = svD0.fProb;
  chi2  = svD0.fChi2;
      
  //alpha = t1.Angle(pCand->fPlab);  // D* pointing angle
  alpha  = t2.Angle(pC->fPlab); // D0 angle
  falpha2 = t3.Angle(pC->fPlab); // D0 angle versus SV2-SV1
  dr = piSlowMom.Angle(fpCand->fPlab); // pislow openinig
  
  if(fVerbose>10 ) {
    cout<<" PVs "<<pvidx<<" "
	<<pv.X()<<" "<<pv.Y()<<" "<<pv.Z()<<" "
	<<sv1.X()<<" "<<sv1.Y()<<" "<<sv1.Z()<<" "
	<<sv2.X()<<" "<<sv2.Y()<<" "<<sv2.Z()<<" "
	<<t1.Mag()<<" "<<t2.Mag()<<" "<<t3.Mag()<<endl;
    cout<<sv.fD3d<<" "<<sv.fDxy<<" "<<svD0.fD3d<<" "<<svD0.fDxy<<endl;
    cout<<alpha<<" "<<t2.Angle(pC->fPlab)<<" "<<t3.Angle(pC->fPlab)<<" "<<t1.Angle(t3)<<endl;
  }

  const bool doHisto = true;
  if(doHisto) {  // Do for all data and for matched MC
    
    //((TH1D*)fHistDir->Get("h61"))->Fill( float(qpi+qpis));
    
    ((TH1D*)fHistDir->Get("h1"))->Fill(t1.Mag());
    ((TH1D*)fHistDir->Get("h2"))->Fill(t2.Mag());
    ((TH1D*)fHistDir->Get("h3"))->Fill(t3.Mag());
    
    ((TH1D*)fHistDir->Get("h4"))->Fill( (pC->fPlab).Angle(fpCand->fPlab) );
    ((TH1D*)fHistDir->Get("h5"))->Fill(sv.fD3d);
    ((TH1D*)fHistDir->Get("h6"))->Fill(svD0.fD3d);
    
    ((TH1D*)fHistDir->Get("h7"))->Fill(t2.Angle(pC->fPlab));
    ((TH1D*)fHistDir->Get("h8"))->Fill(t3.Angle(pC->fPlab));
    ((TH1D*)fHistDir->Get("h9"))->Fill(t1.Angle(t3));
    ((TH1D*)fHistDir->Get("h10"))->Fill(t1.Angle(fpCand->fPlab));
    
    
    ((TH1D*)fHistDir->Get("mds"))->Fill(mdstar);
    ((TH1D*)fHistDir->Get("mdz"))->Fill(mdz);
    ((TH2D*)fHistDir->Get("h2d"))->Fill(mdz, dm);
    ((TH1D*)fHistDir->Get("dm"))->Fill(dm);
    ((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
    ((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
    ((TH1D*)fHistDir->Get("prob"))->Fill(prob);
    ((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
    ((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
    ((TH1D*)fHistDir->Get("alpha2"))->Fill(falpha2);
    ((TH1D*)fHistDir->Get("pt"))->Fill(pt);
    ((TH1D*)fHistDir->Get("dr"))->Fill(dr);
    
    ((TH1D*)fHistDir->Get("ptdz"))->Fill(ptdz);
    ((TH1D*)fHistDir->Get("ptPis"))->Fill(ptPis);
    ((TH1D*)fHistDir->Get("ptPi"))->Fill(ptPi);
    ((TH1D*)fHistDir->Get("ptK"))->Fill(ptK);
    
  }  // if 
  
  
  ((TH1D*)fHistDir->Get("Status"))->Fill(10.);
  count1++;

  // Now the selection cuts cut 
  if(fVerbose>10 ) cout<<"Check pre-cuts "<<endl;
  if(ptPi<4. || ptK<4.) return;  //  does not do aything for data, cand reco is already with 4
  ((TH1D*)fHistDir->Get("Status"))->Fill(11.);
  
  if( (qpi+qpis)==0 ) return; // skip wrong sign decys 
  ((TH1D*)fHistDir->Get("Status"))->Fill(12.);
  
  //if (prob < 0.05) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(13.);
  
  if (ptPis < 0.4) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(14.);
  
  if (dr > 0.25) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(15.);
  
  if (chi2 > 2.0) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(16.);
  
  if (pt < 2) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(17.);
  
  if (alpha > 0.4) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(18.);
  
  if (fls3d < 2) return;
  ((TH1D*)fHistDir->Get("Status"))->Fill(19.);
  
  if( dm<0.130 || dm>0.160 ) return; 
  ((TH1D*)fHistDir->Get("Status"))->Fill(20.);
  
  if( (qk+qpi)!=0 ) return; 
  ((TH1D*)fHistDir->Get("Status"))->Fill(21.);
  
  if(fVerbose>10 ) cout<<"Passed pre-cuts "<<endl;
  count2++;
  
  // Look at muid
  //bool muid1 = goodMuon(pPi);  // true for good  muons 
  //bool muid2 = goodMuon(pK);
  //cout<<" call tight for pion "<<endl;
  bool muid1 = tightMuon(pPi);  // true for good/tight  muons 
  //cout<<" call tight for kaon "<<endl;
  bool muid2 = tightMuon(pK);
  
  //bool mumatch1old = doTriggerMatching(pPi,true); // see if it matches HLT muon
  //bool mumatch2old = doTriggerMatching(pK,true); // see if it matches HLT muon
  
  match1dr = doTriggerMatchingR(pPi,true); // see if it matches HLT muon
  match2dr = doTriggerMatchingR(pK,true);  // see if it matches HLT muon
  bool mumatch1 = (match1dr<0.02); // see if it matches HLT muon
  bool mumatch2 = (match2dr<0.02); // see if it matches HLT muon

  ((TH1D*)fHistDir->Get("dr7"))->Fill(match1dr);
  ((TH1D*)fHistDir->Get("dr7"))->Fill(match2dr);
  if(muid1) ((TH1D*)fHistDir->Get("dr8"))->Fill(match1dr);
  if(muid2) ((TH1D*)fHistDir->Get("dr8"))->Fill(match2dr);
  if(mumatch1) ((TH1D*)fHistDir->Get("dr1"))->Fill(match1dr);
  if(mumatch2) ((TH1D*)fHistDir->Get("dr1"))->Fill(match2dr);
  if(muid1&&mumatch1) ((TH1D*)fHistDir->Get("dr2"))->Fill(match1dr);
  if(muid2&&mumatch2) ((TH1D*)fHistDir->Get("dr2"))->Fill(match2dr);

  //if(mumatch1 != mumatch1old) cout<<" ERROR "<<muid1<<" "<<mumatch1<<" "<<mumatch1old<<" "<<match1dr<<endl;
  //if(mumatch2 != mumatch2old) cout<<" ERROR "<<muid2<<" "<<mumatch2<<" "<<mumatch2old<<" "<<match2dr<<endl;


  if(!useMyTrees) fPreselection = true;  // select this event for the standrad redtree
  //cout<<fPreselection<<" "<<fGoodHLT<<endl;
  
  // kink finder 
  double chiPi = -99.; // pPi->fChi2;
  double chiK = -99.; // pK->fChi2;
  //int mid1 = 0, mid2= 0;
  if (pPi->fMuIndex > -1) {
    chiPi= fpEvt->getMuon(pPi->fMuIndex)->fMuonChi2;
    //mid1 = fpEvt->getMuon(pPi->fMuIndex)->fMuID;
  }
  if (pK->fMuIndex > -1)  {
    chiK = fpEvt->getMuon(pK->fMuIndex)->fMuonChi2;
    //mid2 = fpEvt->getMuon(pK->fMuIndex)->fMuID;
  }
  
  //cout<<fRun<<" "<<fEvt<<" "<<fLS<<" "<<fJSON<<" "<<fGoodHLT<<" "<<fHLTmatch<<endl;
  double dr1 = matchToMuon(pPi,true); // skip same track muons
  double dr2 = matchToMuon(pK,true);

  if( (muid1 || muid2) ) {
    if(fVerbose>0) cout<<"fake muon "<<muid1<<" "<<muid2<<" "<<mumatch1<<" "<<mumatch2<<" "<<match1dr<<" "<<match2dr<<" "<<dr1<<" "<<dr2<<endl;
    
    if(muid1) { // pi misid 
      count3++; if(!mumatch1) count5++;
    }
    if(muid2) { // Kmisid
      count4++; if(!mumatch2) count6++;
    }
  } // if muid
  
  if(muid1) ((TH1D*)fHistDir->Get("dr5"))->Fill(dr1);
  else      ((TH1D*)fHistDir->Get("dr6"))->Fill(dr1);
  if(muid2) ((TH1D*)fHistDir->Get("dr5"))->Fill(dr2);
  else      ((TH1D*)fHistDir->Get("dr6"))->Fill(dr2);
  
  //cout<<" pi mis-id "<<muid1<<" "<<match1dr<<" "<<mumatch1<<" "<<dr1<<" "<<fhltType<<endl;
  //cout<<" k  mis-id "<<muid2<<" "<<match2dr<<" "<<mumatch2<<" "<<dr2<<" "<<fhltType<<endl;

  // Isolation 
  //                       dcaCut(cm) ptCut(GeV)         
  //int close1 = nCloseTracks(fpCand,0.03, 0.5); // around D*
  pair<int, int> pclose; 
  pclose = nCloseTracks(pC,    0.03, 2, 0.5); // around D0
  int close2 = pclose.first; 
  
  //                                         dca   R    Pt
  //double iso1 = isoClassicWithDOCA(fpCand, 0.05,0.7, 0.9); // D*
  double iso2 = isoClassicWithDOCA(pC,     0.05,0.7, 0.9); // D0
  
  
  ftm= 0; // tm;
  fmds=mdstar;
  fmdz=mdz;
  fchi2=chi2;
  falpha=alpha;
  // falpha2 = // not used 
  ffls3d=fls3d;
  fqpis=qpis;
  fdr=dr;
  
  feta=eta;
  fetapi=etaPi;
  fetak=etaK;
  
  fpt=pt;
  fptdz=ptdz;
  fptpis=ptPis;
  fptpi=ptPi;
  fptk=ptK;
  
  fmuid1 = muid1;
  fmuid2 = muid2;
  fmumat1 = mumatch1;
  fmumat2 = mumatch2;
  
  mudr1 = dr1;
  mudr2 = dr2;
  
  fchipi = chiPi;
  fchik = chiK;
  
  fiso = iso2;
  fnclose = close2;
  
  fpvd = pv.Z();
  
  if(fVerbose>0) if(count0%10 == 0) cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;
  else           if(count0%100 == 0) cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;
  
  if(useMyTrees) tree->Fill();
  
}

// ----------------------------------------------------------------------
void candAnaDstar::dumpHFTruthCand(TAnaCand *pC) {
  const bool PRINT = false;

  //TAnaTrack *pT(0);  // needs work 
  //if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  //for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
  //pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
  //pT->dump(); 
  //}

  // Check simple tracks 
  if(PRINT) cout<<"List all tracks "<<fpEvt->nRecTracks()<<"  "<<fpEvt->nSimpleTracks()<<endl;
  // loop over simple tracks  
  for(int itrk=0; itrk<(fpEvt->nSimpleTracks()); itrk++ ) {  // if the simple track exists 
    TSimpleTrack *pTrack = fpEvt->getSimpleTrack(itrk);
    //cout<<itrk<<" "<<pTrack<<endl;
    if(pTrack == 0) continue;

    //pTrack->dump();

    TVector3 mom = pTrack->getP(); // momentum 
    int index = pTrack->getIndex(); // same as itrk
    int q     = pTrack->getCharge();
    int qual  = pTrack->getHighPurity();
    int muonId= pTrack->getMuonID(); // muon id 
    int pvidx = pTrack->getPvIndex();
    int genidx=pTrack->getGenIndex();
    int inds  = pTrack->getIndices();
    int bits  = pTrack->getBits();
    double pt = mom.Perp();

    if(PRINT) cout<<" track "<<itrk<<" "<<index<<" "<<q<<" "<<qual<<" "<<muonId<<" "<<pvidx<<" "<<genidx<<" "
		  <<hex<<inds<<" "<<bits<<dec<<" "<<pt<<" "<<mom.Phi()<<" "<<mom.Eta()<<endl;
 
    ((TH1D*)fHistDir->Get("htrackpt1"))->Fill(pt);

    if(muonId>0) {
      ((TH1D*)fHistDir->Get("htrackpt2"))->Fill(pt);
      if(PRINT) cout<<" muon "<<muonId;
      bool muid  = tightMuon(pTrack);  // does this work on simpleTracks?
      if(muid) {
	if(PRINT) cout<<", tight muon"; 
	((TH1D*)fHistDir->Get("htrackpt3"))->Fill(pt);
      }
      if(PRINT) cout<<endl;
    } // muon 
    
  } // end track loop 


  // Check signal tracks
  if(PRINT) cout<<"List signal tracks "<<pC->fSig1<<" "<<pC->fSig2<<endl;
  for(int it = pC->fSig1; it<=pC->fSig2; it++) {  // loop over index of signal tracka
    TAnaTrack *pT = fpEvt->getSigTrack(it); // this gives TAnaTrack
    //pT->dump(); // so all normal TAnaTracks work 
    int itrk = pT->fIndex; // index of the SimpleTrack
    double pt=pT->fPlab.Perp();
    int muonId = pT->fMuID; // muon id
    if(PRINT)cout <<" signal "<<it<<" "<<itrk <<" "<<pT->fMCID<<" "<<pT->fGenIndex<<" "<<pt<<" "<<pT->fQ<<" "
		  <<pT->fChi2 <<" "<<pT->fDof<<" "<<pT->fTrackQuality<<" "<<pT->fAlgorithm<<" "<<muonId<<" "
		  <<pT->fMuIndex << " "<<pT->fPvIdx<<endl; 

    ((TH1D*)fHistDir->Get("hsigpt1"))->Fill(pt);

    // Check tight muon only if muonId>0
    if(muonId>0) {
      if(PRINT) cout<<" muon "<<muonId;
      ((TH1D*)fHistDir->Get("hsigpt2"))->Fill(pt);

      bool muid  = tightMuon(pT);  // 
      if(muid) {
	if(PRINT) cout<<", tight muon"; 
	((TH1D*)fHistDir->Get("hsigpt3"))->Fill(pt);
      }
      if(PRINT) cout<<endl;
    } // if muon 
  }  // end signal loop 


  // Check muons 
  if(PRINT) cout<<"List muon tracks "<< fpEvt->nMuons()<<endl;
  const float dRMin = 0.1;
  int muHlt = 0, muTightHlt = 0;
  for (int it = 0; it< fpEvt->nMuons(); ++it) { // loop over muons
    TAnaMuon *muon = fpEvt->getMuon(it); // get muon

    // some direct muon methods
    //muon->dump();

    // some direct muon methods
    int muonId = muon->fMuID; // muon ID bits
    int muonIdx = muon->fMuIndex;  // index in the muon list = it
 
    // get track index
    int itrk = muon->fIndex; // index of the SimpleTrack

    TVector3 muonMom = muon->fPlab;
    double ptMuon  = muonMom.Perp();

    ((TH1D*)fHistDir->Get("hmupt1"))->Fill(ptMuon);

    if(itrk<0) continue; // skip muons without inner tracker info
    ((TH1D*)fHistDir->Get("hmupt2"))->Fill(ptMuon);

    if(PRINT) cout<<" muon "<<it<<" index "<<itrk<<" "<<muon->fQ<<" "<<muon->fGlobalPlab.Perp()<<" "<<ptMuon<<" "<<muon->fValidHits
		  <<" "<<muon->fMuonChi2<<" m-id "<<hex<<muonId<<dec<<" "<<muonIdx<<endl;
 
    bool muid  = tightMuon(muon);  // does it work on muons? 
    if(muid) {
      if(PRINT) cout<<" A tight muon"<<endl;  
      ((TH1D*)fHistDir->Get("hmupt3"))->Fill(ptMuon);
    }

    // Not do the rigger matching
    double dR = 0.; // doTriggerMatchingR(muon,true); // see if it matches HLT muon
    ((TH1D*)fHistDir->Get("dr3"))->Fill(dR);
    if(muid) ((TH1D*)fHistDir->Get("dr4"))->Fill(dR);

    if(dR<dRMin) {
      muHlt++;
      if(muid) muTightHlt++;
      if(PRINT) cout<<" Matched to HLT "<<dR<<endl;
    }
    // Going through the SimpleTrack gives the same info
    //TSimpleTrack *pTrack = fpEvt->getSimpleTrack(itrk);
    //cout<<it<<" "<<itrk<<" "<<pTrack->getP().Perp()<<endl; // same as direct muon

  } // END MUON LOOP 

  if(PRINT) cout<<" HLT matched muons in this event "<<muHlt<<" tight "<<muTightHlt<<endl;
  ((TH1D*)fHistDir->Get("h11"))->Fill(float(muHlt));
  ((TH1D*)fHistDir->Get("h12"))->Fill(float(muTightHlt));

}

// ----------------------------------------------------------------------
void candAnaDstar::dumpHFDstarCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 


  // -- D0 daughters
  if (pC->fDau1 < 0) {
    cout << "XXXXXXXXX cannot get daughter cand of " << pC->fType << endl;
    return;
  }

  TAnaCand *pD = fpEvt->getCand(pC->fDau1); 
  cout << "HFDstarCand: idx = " << pC->fIndex << " type = " << pC->fType
       << " m* = " << pC->fMass << " m0 = " << pD->fMass << " dm = " << pC->fMass-pD->fMass << endl;

  // Check D0 tracks
  cout<<" loop over D0 tracks "<<endl; 
  for (int id = pD->fSig1; id <= pD->fSig2; ++id) {

    pT = fpEvt->getSigTrack(id); // this gives TAnaTrack
    //pT->dump(); // so all normal TAnaTracks work 
    //cout <<id<<" "<<pT->fIndex <<" "<<pT->fMCID<<" "<<pT->fGenIndex<<" "<<pT->fPlab.Perp()<<" "<<pT->fQ<<" "
    // <<pT->fChi2 <<" "<<pT->fDof<<" "<<pT->fTrackQuality<<" "<<pT->fAlgorithm<<" "<<pT->fMuID<<" "
    // <<pT->fMuIndex << " "<<pT->fPvIdx<<endl; 

    int index = pT->fIndex; // index of the SimpleTrack
    if (211 == pT->fMCID) {  // pion 
      cout<<" pion index "<<index<<" id "<< pT->fMCID<<" pt "<<pT->fPlab.Perp() << endl;
    } else {  // kaon
      cout<<" kaon index "<<index<<" id "<< pT->fMCID<<" pt "<<pT->fPlab.Perp() << endl;
    }
  }


  // -- slow pion
  //cout<<" simpler way "<<endl;
  int it = pC->fSig1;  // index of the signal track
  pT = fpEvt->getSigTrack(it); // this give TAnaTrack
  cout << " slow pion index "<<it<<" id "<<fpEvt->getSigTrack(pC->fSig1)->fMCID << " pt "<<pT->fPlab.Perp()<<endl; 
  pT->dump(); 

}

// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
int candAnaDstar::truthMatch(TAnaCand *pCand, int verbose) {

  // -- check slow pion
  TAnaTrack *pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
  if(fVerbose>0) cout<<fpEvt->getSigTrack(pCand->fSig1)->fIndex<<" "
		<<pT->fGenIndex<<" "
		<<(fpEvt->getRecTrack(pT->fIndex)->fGenIndex)<<" "<<verbose<<endl;  // same as above 

  if (pT->fGenIndex < 0) {
    if (verbose > 0) cout << "pT->fGenIndex < 0" << endl;
    return 0; 
  }
  TGenCand  *pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex); 
  if (0 == pG) {
    if (verbose > 0) cout << "0 == pG" << endl;
    return 0;
  }
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return 0;
  }

  if(fVerbose>0) cout<< pG->fID<<" "<<endl;  // gen id


  int moSlowPion = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 413 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 413 != pG->fID, pG->fID = " << pG->fID << " moSlowPion = " << moSlowPion << endl;
    return 0;
  }
  if (pG->fDau2 - pG->fDau1 > 1) {
    if (verbose > 0) cout << "pG->fDau2 - pG->fDau1 > 1" << endl;
    return 0; 
  }
  if (verbose > 0) cout << "slow pion " << pT->fIndex << " OK, fGenIndex = " << pT->fGenIndex << endl;

  if(fVerbose>0) cout<< moSlowPion<<" "<<pG->fID<<" "<<endl;


  // -- check D0 
  if (pCand->fDau1 < 0) {
    if (verbose > 0) cout << "no pCand->fDau1" << endl;
    return 0;
  }

  TAnaCand *pC = fpEvt->getCand(pCand->fDau1);
  
  // -- check D0 daughters
  int type(0), moIdx(-1); 
  int count=0, countPi=0, countK=0, missK=0, missPi=0;
  //int daus = (pC->fSig2) - (pC->fSig1);
  //cout<<" D0 with daugthers "<<(daus+1)<<endl;
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    pT = fpEvt->getSigTrack(id);
    type = pT->fMCID;
    pT = fpEvt->getRecTrack(pT->fIndex);
    if (pT->fGenIndex < 0) {
      if (verbose > 0) cout << "no pT->fGenIndex" << endl;
      return 0;
    }
    pG = fpEvt->getGenCand(pT->fGenIndex);
    if (moIdx < 0) {
      moIdx = pG->fMom1;
    } else {
      if (moIdx != pG->fMom1) {
	if (verbose > 0) cout << "moIdx != pG->fMom1" << endl;
	return 0;
      }
    }
    if (verbose > 0) 
      cout << "dau cand sigtrack " << id 
 	   << " with type = " << type 
 	   << " and gen ID = " << pG->fID 
 	   << " at gen idx = " << pT->fGenIndex 
 	   << endl;


    //if (TMath::Abs(type) !=  TMath::Abs(pG->fID)) {
    //if (verbose > 0) cout << "TMath::Abs(type) != TMath::Abs(pG->fID), type = " << type << " pG->fID = " << pG->fID 
    //			    << " track " << pT->fIndex << endl;
    //  return 0;
    // }

    count++;
    if( TMath::Abs(pG->fID) == 321 ) { // Kaon
      countK++;
      //cout <<count<<" "<<countK; 
      if( TMath::Abs(type) != 321 ) {
	//cout << " Kaon identified as pion, type = " << type << " pG->fID = " << pG->fID 
	//   << " track " << pT->fIndex;
	missK++;
	//return 0;
      }
      //cout<<endl;
    }
    if( TMath::Abs(pG->fID) == 211 ) { // Pion
      countPi++;
      //cout <<count<<" "<<countPi; 
      if( TMath::Abs(type) != 211 ) {
	//cout << " Pion identified as kaon, type = " << type << " pG->fID = " << pG->fID 
	//   << " track " << pT->fIndex;
	missPi++;
	//return 0;
      }
      //cout<<endl;
    }

  }


  // -- Get gen-level D0
  if (pG->fMom1 < 0) {
    if (verbose > 0) cout << "pG->fMom1 < 0" << endl;
    return 0; 
  }
  pG = fpEvt->getGenCand(pG->fMom1);
  if (pG->fMom1 != moSlowPion) {
    if (verbose > 0) cout << "pG->fMom1 != moSlowPion" << endl;
    return 0; 
  }
  if (pG->fDau2 - pG->fDau1 > 1) {
    if (verbose > 0) cout << "pG->fDau2 - pG->fDau1 > 1" << endl;
    return 0; 
  }

  if (verbose > 0)   cout << "===> truth matching OK" << endl;


  //cout<<" Found kaons:"<< countK<<" Found pions "<<countPi<<" Found "<<count<<"/"<<(daus+1)
  //  <<" Miss K/Pi "<<missK<<"/"<<missPi<<endl;
  //if(missK!=0 || missPi!=0) return false; // select righ combination 
  //if(missK!=1 || missPi!=1) return false; // select wrong combination 

  if     (missK==0 && missPi==0) return  1; // select righ combination 
  else if(missK==1 && missPi==1) return -1; // select wrong combination 

  return 0; 
}
  

// ----------------------------------------------------------------------
void candAnaDstar::moreBasicCuts() {
  cout << "   candAnaDstar: more basic cuts" << endl;
}


// ----------------------------------------------------------------------
void candAnaDstar::bookHist() {
  //  candAna::bookHist();
  cout << "==>candAnaDstar: bookHist" << endl;

  candAna::bookHist();

  //fHistDir->cd();

  TH1 *h = new TH1D("Status", "Status", 50, -0.5, 49.5);

  // Dstar histos
  h = new TH1D("mds", "m(dstar)", 50, 1.8, 2.5);
  h = new TH1D("mdz", "m(d0)", 70, 1.8, 2.5);
  h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("alpha2", "alpha2", 50, 0., 1.0);
  h = new TH1D("dr","dr",100, 0., 1);
  h = new TH1D("pt", "pT", 50, 0., 25);
  h = new TH1D("ptdz", "pT", 50, 0., 25);
  h = new TH1D("ptK",   "pT", 50, 0., 10);
  h = new TH1D("ptPi",  "pT", 50, 0., 10);
  h = new TH1D("ptPis", "pT", 50, 0., 5);
  TH2D *h2 = new TH2D("h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);


  // overview histos 
  h = new TH1D("htrackpt1","all track pT", 100, 0., 100);
  h = new TH1D("htrackpt2","muon track pT", 100, 0., 100);
  h = new TH1D("htrackpt3","tight muon track pT", 100, 0., 100);
  h = new TH1D("hsigpt1","all signal pT", 100, 0., 100);
  h = new TH1D("hsigpt2","muon signal pT", 100, 0., 100);
  h = new TH1D("hsigpt3","tight muon signal pT", 100, 0., 100);
  h = new TH1D("hmupt1","all muon pT", 100, 0., 100);
  h = new TH1D("hmupt2","tracker muon pT", 100, 0., 100);
  h = new TH1D("hmupt3","tight muon pT", 100, 0., 100);

  // test histos 
  h = new TH1D("dr1", "dr1", 400, 0., 2.);
  h = new TH1D("dr2", "dr2", 400, 0., 2.);
  h = new TH1D("dr3", "dr3", 400, 0., 2.);
  h = new TH1D("dr4", "dr4", 400, 0., 2.);
  h = new TH1D("dr5", "dr5", 400, 0., 2.);
  h = new TH1D("dr6", "dr6", 400, 0., 2.);
  h = new TH1D("dr7", "dr7", 400, 0., 2.);
  h = new TH1D("dr8", "dr8", 400, 0., 2.);
  h = new TH1D("dr9", "dr9", 400, 0., 2.);

  h = new TH1D("pvidx", "PVidx", 50, 0., 50.);
  h = new TH1D("h1", "h1", 100, 0., 1);
  h = new TH1D("h2", "h2", 100, 0., 1);
  h = new TH1D("h3", "h3", 100, 0., 1);
  h = new TH1D("h4", "h4", 350, 0., 3.5);
  h = new TH1D("h5", "h5", 100, 0., 1);
  h = new TH1D("h6", "h6", 100, 0., 1);
  h = new TH1D("h7", "h7", 350, 0., 3.5);
  h = new TH1D("h8", "h8", 350, 0., 3.5);
  h = new TH1D("h9", "h9", 350, 0., 3.5);
  h = new TH1D("h10","h10",100, 0., 1);
  h = new TH1D("h11", "h11", 10, 0., 10);
  h = new TH1D("h12", "h12", 10, 0., 10);

  // Ntuples 
  tree = new TTree("dstar","dstar");
  // special to dstar
  tree->Branch("ftm",&ftm,"ftm/I");
  tree->Branch("fmuid1",&fmuid1,"fmuid1/O");
  tree->Branch("fmuid2",&fmuid2,"fmuid2/O");
  tree->Branch("fmumat1",&fmumat1,"fmumat1/O");
  tree->Branch("fmumat2",&fmumat2,"fmumat2/O");
  tree->Branch("fmds",&fmds,"fmds/F");
  tree->Branch("fmdz",&fmdz,"fmdz/F");

  tree->Branch("ffls3d",&ffls3d,"ffls3d/F");
  tree->Branch("fchi2",&fchi2,"fchi2/F");
  tree->Branch("falpha",&falpha,"falpha/F");
  tree->Branch("falpha2",&falpha2,"falpha2/F");
  tree->Branch("fqpis",&fqpis,"fqpis/F");
  tree->Branch("fdr",&fdr,"fdr/F");

  tree->Branch("fpt",&fpt,"fpt/F");
  tree->Branch("fptdz",&fptdz,"fptdz/F");
  tree->Branch("fptpis",&fptpis,"fptpis/F");
  tree->Branch("fptpi",&fptpi,"fptpi/F");
  tree->Branch("fptk",&fptk,"fptk/F");

  tree->Branch("feta",&feta,"feta/F");
  tree->Branch("fetapi",&fetapi,"fetapi/F");
  tree->Branch("fetak",&fetak,"fetak/F");


  tree->Branch("fchipi",&fchipi,"fchipi/F");
  tree->Branch("fchik",&fchik,"fchik/F");
  tree->Branch("fiso",&fiso,"fiso/F");
  tree->Branch("fnclose",&fnclose,"fnclose/I");

  tree->Branch("mudr1",&mudr1,"fmudr1/F");
  tree->Branch("mudr2",&mudr2,"fmudr2/F");

  tree->Branch("hltdr1",    &match1dr,      "hltdr1/F");
  tree->Branch("hltdr2",    &match2dr,      "hltdr2/F");


  // general 
  tree->Branch("fpvd",&fpvd,"fpvd/F");
  tree->Branch("run",     &fRun,               "run/L");
  tree->Branch("json",    &fJSON,              "json/O");
  tree->Branch("evt",     &fEvt,               "evt/L");
  tree->Branch("ls",      &fLS,                "ls/I");
  //t->Branch("tm",      &fCandTM,            "tm/I");
  //t->Branch("pr",      &fGenBpartial,       "pr/I"); 
  //t->Branch("procid",  &fProcessType,       "procid/I");
  tree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  tree->Branch("pvn",    &fPvN,               "pvn/I");
  //t->Branch("cb",      &fCowboy,            "cb/O");
  //t->Branch("rr",      &fRunRange,          "rr/I");
  //t->Branch("bdt",     &fBDT,               "bdt/D");
  //t->Branch("npv",     &fPvN,               "npv/I");
  //t->Branch("pvw8",    &fPvAveW8,           "pvw8/D");
  tree->Branch("hltm",    &fHLTmatch,          "hltm/O");
  tree->Branch("hlttype",   &fhltType,        "hltype/I");

  // Add variables to the standard redtrees 
  fTree->Branch("ftm",&ftm,"ftm/I");
  fTree->Branch("fmuid1",&fmuid1,"fmuid1/O");
  fTree->Branch("fmuid2",&fmuid2,"fmuid2/O");
  fTree->Branch("fmumat1",&fmumat1,"fmumat1/O");
  fTree->Branch("fmumat2",&fmumat2,"fmumat2/O");
  fTree->Branch("fmds",&fmds,"fmds/F");
  fTree->Branch("fmdz",&fmdz,"fmdz/F");

  fTree->Branch("ffls3d",&ffls3d,"ffls3d/F");
  fTree->Branch("fchi2",&fchi2,"fchi2/F");
  fTree->Branch("falpha",&falpha,"falpha/F");
  //fTree->Branch("falpha2",&falpha2,"falpha2/F");
  fTree->Branch("fqpis",&fqpis,"fqpis/F");
  fTree->Branch("fdr",&fdr,"fdr/F");

  fTree->Branch("fpt",&fpt,"fpt/F");
  fTree->Branch("fptdz",&fptdz,"fptdz/F");
  fTree->Branch("fptpis",&fptpis,"fptpis/F");
  fTree->Branch("fptpi",&fptpi,"fptpi/F");
  fTree->Branch("fptk",&fptk,"fptk/F");

  fTree->Branch("feta",&feta,"feta/F");
  fTree->Branch("fetapi",&fetapi,"fetapi/F");
  fTree->Branch("fetak",&fetak,"fetak/F");

  fTree->Branch("fchipi",&fchipi,"fchipi/F");
  fTree->Branch("fchik",&fchik,"fchik/F");
  fTree->Branch("fiso",&fiso,"fiso/F");
  fTree->Branch("fnclose",&fnclose,"fnclose/I");

  fTree->Branch("mudr1",&mudr1,"fmudr1/F");
  fTree->Branch("mudr2",&mudr2,"fmudr2/F");

  fTree->Branch("hltdr1",    &match1dr,      "hltdr1/F");
  fTree->Branch("hltdr2",    &match2dr,      "hltdr2/F");

}


// ----------------------------------------------------------------------
void candAnaDstar::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaDstar: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;
  //int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    //ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

   
  }

}

