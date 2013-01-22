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
void candAnaDstar::candAnalysis() {
  static int count0=0, count1=0, count2=0, count3=0, count4=0, count5=0;

  if (0 == fpCand) return;

  candAna::candAnalysis();

  count0++;

  ((TH1D*)fHistDir->Get("Status"))->Fill(0.);
  //cout<<" Dstar candidate "<<fpCand->fType<<" in event "<<fEvt<<endl;

  if (fVerbose>10) {
    cout << "DUMP HFDstarCandidate with mass = " << fpCand->fMass << endl;
    //dumpHFTruthCand(pCand); 
    dumpHFDstarCand(fpCand); 
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

  int pPisId = fpEvt->getSigTrack(fpCand->fSig1)->fIndex; // slow pion index
  pPis = fpEvt->getRecTrack(pPisId); // slow pi from D*
  TVector3 piSlowMom = pPis->fPlab; // slow pi momentum vector 

//       if (fVerbose>0 && tm) cout<<" slow pi "<<pPisId<<" "<< pPis->fPlab.Perp()  << " "<<pPis->fMCID<<endl;
//       if (fVerbose>0 && tm) cout<<" D0 "<<pC->fType<<endl;

  // loop over D0 tracks 
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    int index = fpEvt->getSigTrack(id)->fIndex;
    // 	if (fVerbose>0 && tm) cout<<id<<" "<<index<<" ";

    if (211 == fpEvt->getSigTrack(id)->fMCID) {  // pion 
      pPi = fpEvt->getRecTrack(index);
      //indexPi = index;
      //if (fVerbose>0 && tm) cout<< pPi->fMCID<<" "<<pPi->fPlab.Perp() << endl;
    } else {  // kaon
      pK = fpEvt->getRecTrack(index);
      //indexK = index;
      // 	  if (fVerbose>0 && tm) cout<< pK->fMCID<<" "<<pK->fPlab.Perp() << endl;
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

  // direction vector from PV to D0 decay vertex 
  //TVector3 svpv1(sv.fPoint - fpEvt->getPV(pvidx)->fPoint); 
  //TVector3 svpv(sv1 - pv1); 
  //cout<<svpv.Z()<<" "<<svpv1.Z()<<endl;

//  ONLY FOR MC
//       if(ok) {  // Do here comparison between RECO and GEN quantities
// 	double c1 = (PV-pv1).Mag();
// 	double c2 = (DSVertex-sv1).Mag();
// 	double c3 = (DZVertex-sv2).Mag();

// 	//cout<<" RECO-MV vertex distance "<<c1<<" "<<c2<<" "<<c3<<endl;
// 	((TH1D*)fHistDir->Get("h31"))->Fill(c1);
// 	((TH1D*)fHistDir->Get("h32"))->Fill(c2);
// 	((TH1D*)fHistDir->Get("h33"))->Fill(c3);

// 	((TH1D*)fHistDir->Get("h60"))->Fill( float(qpi+qpis));

// 	double a10 = PiSlowMom.Angle(pPis->fPlab);
// 	double a11 = PiMom.Angle(pPi->fPlab);
// 	double a12 = KMom.Angle(pK->fPlab);
// 	double a13 = DSMom.Angle(pCand->fPlab);
// 	double a14 = DZMom.Angle(pC->fPlab);

// 	((TH1D*)fHistDir->Get("h34"))->Fill(a10);
// 	((TH1D*)fHistDir->Get("h35"))->Fill(a11);
// 	((TH1D*)fHistDir->Get("h36"))->Fill(a12);
// 	((TH1D*)fHistDir->Get("h37"))->Fill(a13);
// 	((TH1D*)fHistDir->Get("h38"))->Fill(a14);

// 	double t1len = t1.Mag();
// 	double t2len = t2.Mag();
// 	double t1lenMC = (DSVertex-PV).Mag();
// 	double t2lenMC = (DZVertex-PV).Mag();
// 	((TH1D*)fHistDir->Get("h26"))->Fill( abs(t1len-t1lenMC) );
// 	((TH1D*)fHistDir->Get("h27"))->Fill( abs(t2len-t2lenMC) );

//       } // if OK


      
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

      if(fVerbose>9 ) {
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

      // Now the selection cuts cut 
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

      count1++;

      // Look at muid
      //bool muid1 = goodMuon(pPi);  // true for good  muons 
      //bool muid2 = goodMuon(pK);
      bool muid1 = tightMuon(pPi);  // true for good/tight  muons 
      bool muid2 = tightMuon(pK);

      bool mumatch1 = doTriggerMatching(pPi,true); // see if it matches HLT muon
      bool mumatch2 = doTriggerMatching(pK,true); // see if it matches HLT muon

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
      if( muid1 || muid2) {
	//cout<<"fake muon "<<muid1<<" "<<muid2<<" "<<mumatch1<<" "<<mumatch2<<" "<<mid1<<" "<<mid2
	//  <<" "<<indexPi<<" "<<indexK<<endl;
	if(muid1) { // pi misid 
	  count2++; if(!mumatch1) count4++;
	}
	if(muid2) { // Kmisid
	  count3++; if(!mumatch2) count5++;
	}
      } // if muid

      double dr1 = matchToMuon(pPi);
      double dr2 = matchToMuon(pK);
      ((TH1D*)fHistDir->Get("mu_match_dr1"))->Fill(dr1);
      ((TH1D*)fHistDir->Get("mu_match_dr1"))->Fill(dr2);
      if(muid1) ((TH1D*)fHistDir->Get("mu_match_dr2"))->Fill(dr1);
      if(muid2) ((TH1D*)fHistDir->Get("mu_match_dr2"))->Fill(dr2);

      // Isolation 
      //                       dcaCut(cm) ptCut(GeV)         
      //int close1 = nCloseTracks(fpCand,0.03, 0.5); // around D*
      int close2 = nCloseTracks(pC,    0.03, 0.5); // around D0
      //                                      dca   R    Pt
      //double iso1 = isoClassicWithDOCA(fpCand, 0.05,0.7, 0.9); // D*
      double iso2 = isoClassicWithDOCA(pC,     0.05,0.7, 0.9); // D0


      ftm= 0; // tm;
      fmds=mdstar;
      fmdz=mdz;
      fchi2=chi2;
      falpha=alpha;
      // falpha2 = 
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

      if(count0%100 == 0) 
	cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<endl;

      tree->Fill();

}

// ----------------------------------------------------------------------
void candAnaDstar::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 
  if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    pT->dump(); 
  }
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
  //  if (fVerbose > -1) cout << "   mdz = " << pD->fMass << endl;

  // -- slow pion
  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pC->fSig1)->fIndex); 
  cout << fpEvt->getSigTrack(pC->fSig1)->fMCID << " " ; 
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
//   h = new TH1D("full_mds", "m(dstar)", 50, 1.8, 2.5);
//   h = new TH1D("full_mdz", "m(d0)", 50, 1.8, 2.5);
//   h = new TH1D("full_dm", "delta(m)", 60, 0.13, 0.16);
//   h = new TH1D("full_ncand", "ncand", 200, 0., 200);
//   h = new TH1D("full_fls3d", "fls3d", 60, 0., 20);
//   h = new TH1D("full_flsxy", "flsxy", 60, 0., 20);
//   h = new TH1D("full_prob", "prob(chi2/dof)", 100, 0., 1.);
//   h = new TH1D("full_chi2", "chi2", 100, 0., 10.);
//   h = new TH1D("full_alpha", "alpha", 350, 0., 3.5);
//   h = new TH1D("full_dr","dr",100, 0., 1);
//   h = new TH1D("full_pt", "pT", 50, 0., 25); 
//   h = new TH1D("full_ptdz", "pT", 50, 0., 25);
//   h = new TH1D("full_ptK",   "pT", 50, 0., 10);
//   h = new TH1D("full_ptPi",  "pT", 50, 0., 10);
//   h = new TH1D("full_ptPis", "pT", 50, 0., 5);
//   TH2D *h2 = new TH2D("full_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

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

  h = new TH1D("mu_match_dr", "mu_match_dr", 100, 0., 1.);
  h = new TH1D("mu_match_pt", "mu_match_pt", 100, 0., 20.);

  h = new TH1D("mu_match_dr1", "mu_match_dr1", 100, 0., 1.);
  h = new TH1D("mu_match_dr2", "mu_match_dr2", 100, 0., 1.);

//   h = new TH1D("all_mds", "m(dstar)", 70, 1.8, 2.5);
//   h = new TH1D("all_mdz", "m(d0)", 70, 1.8, 2.5);
//   h = new TH1D("all_dm", "delta(m)", 60, 0.13, 0.16);
//   h = new TH1D("all_ncand", "ncand", 200, 0., 200);
//   h = new TH1D("all_fls3d", "fls3d", 60, 0., 20);
//   h = new TH1D("all_flsxy", "flsxy", 60, 0., 20);
//   h = new TH1D("all_prob", "prob(chi2/dof)", 100, 0., 1.);
//   h = new TH1D("all_chi2", "chi2", 100, 0., 10.);
//   h = new TH1D("all_alpha", "alpha", 50, 0., 1.0);
//   h = new TH1D("all_dr","dr",100, 0., 1);
//   h = new TH1D("all_pt", "pT", 50, 0., 25);
//   h = new TH1D("all_ptdz", "pT", 50, 0., 25);
//   h = new TH1D("all_ptK",   "pT", 50, 0., 10);
//   h = new TH1D("all_ptPi",  "pT", 50, 0., 10);
//   h = new TH1D("all_ptPis", "pT", 50, 0., 5);

//   h2 = new TH2D("all_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);


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

//   h = new TH1D("h11", "h11", 100, 0., 1);
//   h = new TH1D("h12", "h12", 100, 0., 1);
//   h = new TH1D("h13", "h13", 100, 0., 1);
//   h = new TH1D("h14", "h14", 100, 0., 10);
//   h = new TH1D("h15", "h15", 100, 0., 1);
//   h = new TH1D("h16", "h16", 100, 0., 1);
//   h = new TH1D("h17", "h17", 350, 0., 3.5);
//   h = new TH1D("h18", "h18", 350, 0., 3.5);
//   h = new TH1D("h19", "h19", 350, 0., 3.5);
//   h = new TH1D("h20", "h20", 350, 0., 3.5);
//   h = new TH1D("h21", "h21", 350, 0., 3.5);
//   h = new TH1D("h22", "h22", 350, 0., 3.5);
//   h = new TH1D("h23", "h23",   5,-2.5, 2.5);
//   h = new TH1D("h24", "h24", 350, 0., 3.5);
//   h = new TH1D("h25", "h25", 350, 0., 3.5);
//   h = new TH1D("h26", "h26", 100, 0., 1);
//   h = new TH1D("h27", "h27", 100, 0., 1);
//   h = new TH1D("h28", "h28", 50, 0., 20);
//   h = new TH1D("h29", "h29", 50, 0., 20);
//   h = new TH1D("h30", "h30", 50, 0., 10);

//   h = new TH1D("h31", "h31", 100, 0., 1);
//   h = new TH1D("h32", "h32", 100, 0., 1);
//   h = new TH1D("h33", "h33", 100, 0., 1);
//   h = new TH1D("h34", "h34", 350, 0., 0.35);
//   h = new TH1D("h35", "h35", 350, 0., 0.35);
//   h = new TH1D("h36", "h36", 350, 0., 0.35);
//   h = new TH1D("h37", "h37", 350, 0., 0.35);
//   h = new TH1D("h38", "h38", 350, 0., 0.35);

//   h = new TH1D("h39", "h39", 50, 0., 10);

//   h = new TH1D("h41", "h41",50, 0.0, 20.);
//   h = new TH1D("h42", "h42",100, 1.0, 2.5);
//   h = new TH1D("h43", "h43",50, 0.0, 20.);
//   h = new TH1D("h44", "h44",100, 1.0, 2.5);

//   h = new TH1D("h45", "h45",100,  1.,2.5);
//   h = new TH1D("h46", "h46",100,  1.,2.5);
//   h = new TH1D("h47", "h47",100,  0., 1.0);
//   h = new TH1D("h48", "h48",50,  0., 20.);
//   h = new TH1D("h49", "h49",200,  0.13, 0.18);
//   h = new TH1D("h50", "h50",200,  0.13, 0.18);

//   h2 = new TH2D("h51", "h51", 40, 0., 20., 80, 1.5, 2.3);
//   h2 = new TH2D("h52", "h52", 35, 0., 3.5, 80, 1.5, 2.3);
//   h2 = new TH2D("h53", "h53", 40, 0., 20., 80, 1.5, 2.3);
//   h2 = new TH2D("h54", "h53", 40, 0., 10., 80, 1.5, 2.3);
//   h2 = new TH2D("h55", "h55", 40, 0., 10., 80, 1.5, 2.3);
//   h2 = new TH2D("h56", "h56", 40, 0.,  2., 80, 1.5, 2.3);
//   h2 = new TH2D("h57", "h57", 40, 0.,0.4,  80, 1.5, 2.3);

//   h = new TH1D("h60", "h60", 5,-2.5,2.5);
  //h = new TH1D("h61", "h61", 5,-2.5,2.5);
//   h = new TH1D("h62", "h62", 5,-2.5,2.5);
//   //h = new TH1D("h63", "h63", 5,-2.5,2.5);
//   //h = new TH1D("h64", "h64", 5,-2.5,2.5);
//   //h = new TH1D("h65", "h65", 5,-2.5,2.5);
//   //h = new TH1D("h66", "h66", 5,-2.5,2.5);
//   //h = new TH1D("h67", "h67", 5,-2.5,2.5);

//   h = new TH1D("h71", "dm",40,0.135,0.155);
//   h = new TH1D("h72", "dm",40,0.135,0.155);
//   h = new TH1D("h73", "dm",40,0.135,0.155);
//   h = new TH1D("h74", "dm",40,0.135,0.155);
//   h = new TH1D("h75", "dm",40,0.135,0.155);
//   h = new TH1D("h76", "dm",40,0.135,0.155);
//   h = new TH1D("h77", "dm",40,0.135,0.155);

  tree = new TTree("dstar","dstar");
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


  tree->Branch("fpvd",&fpvd,"fpvd/F");
  tree->Branch("fchipi",&fchipi,"fchipi/F");
  tree->Branch("fchik",&fchik,"fchik/F");
  tree->Branch("fiso",&fiso,"fiso/F");
  tree->Branch("fnclose",&fnclose,"fnclose/I");

  tree->Branch("mudr1",&mudr1,"fmudr1/F");
  tree->Branch("mudr2",&mudr2,"fmudr2/F");

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

