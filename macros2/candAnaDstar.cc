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
  int qds =0, qk =0, qpi=0, qpis=0;
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

void candAnaDstar::evtAnalysis(TAna01Event *evt) {
  fpEvt = evt; 
  fcands=0;
  int count = 0;

  TAnaCand *pCand(0), *pC(0), *pC1(0);
  TAnaTrack *pK, *pPi, *pPis; 
  //cout << "Evt: " << fEvt << " ----------------------------------------------------------------------" << endl;
  // -- loop over all seq vtx fit candidates for D*
  double mdstar(0.), mdz(0.), dm(0.), fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), pt(0.), dr(0.); 
  int qk =0, qpi=0, qpis=0;
  int tm=0; 
  int ncand(0); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {  // just count candidates
    pCand = fpEvt->getCand(iC);
    if (300054 == pCand->fType) ++ncand; 
    else if( pCand->fType == 54 || pCand->fType == -54) {
      if (fVerbose>0) {
	TVector3 s = pCand->fVtx.fPoint;
	cout << " -> " << iC <<" "<< pCand->fType<<" "<<pCand->fMom;
	cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	cout << " mass: " << pCand->fMass << " " << s.Z() <<" "<<pCand->fPlab.Mag()<<endl;
	cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
	dumpHFTruthCand(pCand); 
      }
    }  
  }


  if(fVerbose>0) cout<<" num of cands "<<ncand<<" "<<fVerbose<<" "<<fIsMC<<endl;

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    bool ok = false;
    //     if (54 == pCand->fType) {
    //       cout << "DUMP HFTruthCandidate with mass = " << pCand->fMass << endl;
    //       dumpHFTruthCand(pCand); 
    //     }

    if (300054 == pCand->fType) {  //select D*+

      if(fVerbose>0) dumpHFDstarCand(pCand); 

      // -- D0 
      if (pCand->fDau1 < 0) {
	if (pCand->fDau2 < 0) {
	  if(fVerbose>0) {cout << "pCand->fDauX = -1!!! " << pCand->fType << endl; pCand->dump();}
	}
	continue;  // skip if no daughters 
      }

      pC = fpEvt->getCand(pCand->fDau1);  // D0?
      pK = 0;
      pPi =0; 

      int pPisId = fpEvt->getSigTrack(pCand->fSig1)->fIndex;
      pPis = fpEvt->getRecTrack(pPisId); // slow pi from D*
      TVector3 piSlowMom = pPis->fPlab; // slow pi momentum vector 

      if (fVerbose>0 && tm) cout<<" slow pi "<<pPisId<<" "<< pPis->fPlab.Perp()  << " "<<pPis->fMCID<<endl;
      if (fVerbose>0 && tm) cout<<" D0 "<<pC->fType<<endl;

      for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
	int index = fpEvt->getSigTrack(id)->fIndex;
	if (fVerbose>0 && tm) cout<<id<<" "<<index<<" ";

	if (211 == fpEvt->getSigTrack(id)->fMCID) {
	  pPi = fpEvt->getRecTrack(index);
	  if (fVerbose>0 && tm) cout<< pPi->fMCID<<" "<<pPi->fPlab.Perp() << endl;
	} else {
	  pK = fpEvt->getRecTrack(index);
	  if (fVerbose>0 && tm) cout<< pK->fMCID<<" "<<pK->fPlab.Perp() << endl;
	}
      }
	

      tm = 0; 
      // truthMatch return an interger, 0-no match, 1-correct match, -1-pi&K swapped
      if(fIsMC) tm = truthMatch(pCand,fVerbose); // check truth matching
      if ( tm == 1 ) { // do only for MC
	ok = anaMC(evt);
	if (fVerbose>0) {
	  cout << " -> " << pCand->fType;
	  cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	  cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	  cout << " mass: " << pCand->fMass << " " << ok<<endl;
	
	  cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
	  dumpHFDstarCand(pCand); 
	  fpEvt->dumpGenBlock(); 
	}
      } // end if

      double ptdz  = pC->fPlab.Perp();
      double ptPis = pPis->fPlab.Perp();
      double ptPi  = pPi->fPlab.Perp();
      double ptK   = pK->fPlab.Perp();

      qk   = pK->fQ;
      qpi  = pPi->fQ;
      qpis = pPis->fQ;

      mdstar = pCand->fMass; 
      mdz = pC->fMass;
      dm = mdstar - mdz;

      // D* vertex
      TAnaVertex sv = pCand->fVtx;  // D* vertex
      int pvidx = (pCand->fPvIdx > -1? pCand->fPvIdx : 0);  // D* PV index

      // D0 vertex
      TAnaVertex svD0 = pC->fVtx;  // D0 vertex
      int pvidx2 = (pC->fPvIdx > -1? pC->fPvIdx : 0);  // D0 PV index

      ((TH1D*)fHistDir->Get("pvidx"))->Fill(float(pvidx));

      TVector3 pv1 =  fpEvt->getPV(pvidx)->fPoint;
      TVector3 pv2 =  fpEvt->getPV(pvidx2)->fPoint;

      TVector3 sv1 = sv.fPoint;
      TVector3 sv2 = svD0.fPoint;
      TVector3 t1(sv1-pv1), t2(sv2-pv1), t3(sv2-sv1), t4(pv1-pv2);

      if(ok) {  // Do here comparison between RECO and GEN quantities
	double c1 = (PV-pv1).Mag();
	double c2 = (DSVertex-sv1).Mag();
	double c3 = (DZVertex-sv2).Mag();

	//cout<<" RECO-MV vertex distance "<<c1<<" "<<c2<<" "<<c3<<endl;
	((TH1D*)fHistDir->Get("h31"))->Fill(c1);
	((TH1D*)fHistDir->Get("h32"))->Fill(c2);
	((TH1D*)fHistDir->Get("h33"))->Fill(c3);

	((TH1D*)fHistDir->Get("h60"))->Fill( float(qpi+qpis));

	double a10 = PiSlowMom.Angle(pPis->fPlab);
	double a11 = PiMom.Angle(pPi->fPlab);
	double a12 = KMom.Angle(pK->fPlab);
	double a13 = DSMom.Angle(pCand->fPlab);
	double a14 = DZMom.Angle(pC->fPlab);

	((TH1D*)fHistDir->Get("h34"))->Fill(a10);
	((TH1D*)fHistDir->Get("h35"))->Fill(a11);
	((TH1D*)fHistDir->Get("h36"))->Fill(a12);
	((TH1D*)fHistDir->Get("h37"))->Fill(a13);
	((TH1D*)fHistDir->Get("h38"))->Fill(a14);

	double t1len = t1.Mag();
	double t2len = t2.Mag();
	double t1lenMC = (DSVertex-PV).Mag();
	double t2lenMC = (DZVertex-PV).Mag();
	((TH1D*)fHistDir->Get("h26"))->Fill( abs(t1len-t1lenMC) );
	((TH1D*)fHistDir->Get("h27"))->Fill( abs(t2len-t2lenMC) );

      } // if OK


      TVector3 svpv(sv.fPoint - fpEvt->getPV(pvidx)->fPoint); 
      
      //fls3d = sv.fD3d/sv.fD3dE; 
      //flsxy = sv.fDxy/sv.fDxyE; 
      fls3d = svD0.fD3d/svD0.fD3dE; 
      flsxy = svD0.fDxy/svD0.fDxyE; 
      prob  = svD0.fProb;
      chi2  = svD0.fChi2;
      pt    = pCand->fPlab.Perp(); 

      //alpha = svpv.Angle(pCand->fPlab);  // D* pointing angle
      alpha = t2.Angle(pC->fPlab); // D0 angle
      dr = piSlowMom.Angle(pCand->fPlab); // pislow openinig

      if(fVerbose>9 && tm==1) {
	cout<<" PVs "<<pvidx<<" "<<pvidx2<<" "
	    <<pv1.X()<<" "<<pv1.Y()<<" "<<pv1.Z()<<" "
	    <<sv1.X()<<" "<<sv1.Y()<<" "<<sv1.Z()<<" "
	    <<sv2.X()<<" "<<sv2.Y()<<" "<<sv2.Z()<<" "
	    <<t1.Mag()<<" "<<t2.Mag()<<" "<<t3.Mag()<<" "<<t4.Mag()<<endl;
	cout<<sv.fD3d<<" "<<sv.fDxy<<" "<<svD0.fD3d<<" "<<svD0.fDxy<<endl;
	cout<<alpha<<" "<<t2.Angle(pC->fPlab)<<" "<<t3.Angle(pC->fPlab)<<" "<<t1.Angle(t3)<<endl;
      }

 

      if(!fIsMC || tm==1) {  // Do for all data and for matched MC

	((TH1D*)fHistDir->Get("h61"))->Fill( float(qpi+qpis));

	((TH1D*)fHistDir->Get("h1"))->Fill(t1.Mag());
	((TH1D*)fHistDir->Get("h2"))->Fill(t2.Mag());
	((TH1D*)fHistDir->Get("h3"))->Fill(t3.Mag());
	((TH1D*)fHistDir->Get("h4"))->Fill( (pC->fPlab).Angle(pCand->fPlab) );
	((TH1D*)fHistDir->Get("h5"))->Fill(sv.fD3d);
	((TH1D*)fHistDir->Get("h6"))->Fill(svD0.fD3d);
	((TH1D*)fHistDir->Get("h7"))->Fill(t2.Angle(pC->fPlab));
	((TH1D*)fHistDir->Get("h8"))->Fill(t3.Angle(pC->fPlab));
	((TH1D*)fHistDir->Get("h9"))->Fill(t1.Angle(t3));

	//((TH1D*)fHistDir->Get("h10"))->Fill( piSlowMom.Angle(pC->fPlab));
	//((TH1D*)fHistDir->Get("h22"))->Fill(t3.Angle(pCand->fPlab));
	//((TH1D*)fHistDir->Get("h24"))->Fill(t2.Angle(pCand->fPlab));

	((TH1D*)fHistDir->Get("full_mds"))->Fill(mdstar);
	((TH1D*)fHistDir->Get("full_mdz"))->Fill(mdz);
	((TH2D*)fHistDir->Get("full_h2d"))->Fill(mdz, dm);
	((TH1D*)fHistDir->Get("full_dm"))->Fill(dm);
	((TH1D*)fHistDir->Get("full_fls3d"))->Fill(fls3d);
	((TH1D*)fHistDir->Get("full_flsxy"))->Fill(flsxy);
	((TH1D*)fHistDir->Get("full_prob"))->Fill(prob);
	((TH1D*)fHistDir->Get("full_chi2"))->Fill(chi2);
	((TH1D*)fHistDir->Get("full_alpha"))->Fill(alpha);
	((TH1D*)fHistDir->Get("full_pt"))->Fill(pt);
	((TH1D*)fHistDir->Get("full_dr"))->Fill(dr);
		
	((TH1D*)fHistDir->Get("full_ptdz"))->Fill(ptdz);
	((TH1D*)fHistDir->Get("full_ptPis"))->Fill(ptPis);
	((TH1D*)fHistDir->Get("full_ptPi"))->Fill(ptPi);
	((TH1D*)fHistDir->Get("full_ptK"))->Fill(ptK);
	
      }  // if 


      // Now the selection cuts cut 
      if(ptPi<4. || ptK<4.) continue;  //  does not do aything for data, cand reco is already with 4
      if( (qpi+qpis)==0 ) continue; // skip wrong sign decys 
      //if (prob < 0.05) continue;
      if (ptPis < 0.4) continue;
      if (dr > 0.25) continue;
      if (chi2 > 2.0) continue;
      if (pt < 2) continue;
      if (alpha > 0.4) continue;
      if (fls3d < 2) continue;

      count++; // count selected candidates 

      // Look at muid
      bool muid1 = goodMuon(pPi);  // true for good  muons 
      bool muid2 = goodMuon(pK);
      bool mutid1 = tightMuon(pPi);  // true for good/tight  muons 
      bool mutid2 = tightMuon(pK);

      //if( muid1 || muid2) cout<<"fake muon "<<muid1<<" "<<muid2;

      // Save in a tree, save only masses between 130-160MeV
      if( (fcands<10) && (dm>0.130) && (dm<0.160) ) {
	ftm[fcands] = tm;
	fmds[fcands]=mdstar;
	fmdz[fcands]=mdz;

	fchi2[fcands]=chi2;
	falpha[fcands]=alpha;
	ffls3d[fcands]=fls3d;
	fqpis[fcands]=qpis;
	fdr[fcands]=dr;

	fpt[fcands]=pt;
	fptdz[fcands]=ptdz;
	fptpis[fcands]=ptPis;
	fptpi[fcands]=ptPi;
	fptk[fcands]=ptK;

	int mu = 0;
	if(muid1)  mu += 1;
	if(muid2)  mu += 2;
	if(mutid1) mu += 100;
	if(mutid2) mu += 200;

	fmu[fcands] = mu;

	fcands++;
      } // if fcands

      //if( muid1 || muid2) cout<<endl;
      // Histogram 
      ((TH1D*)fHistDir->Get("all_mds"))->Fill(mdstar);
      ((TH1D*)fHistDir->Get("all_mdz"))->Fill(mdz);
      ((TH2D*)fHistDir->Get("all_h2d"))->Fill(mdz, dm);
      ((TH1D*)fHistDir->Get("all_dm"))->Fill(dm);
      ((TH1D*)fHistDir->Get("all_fls3d"))->Fill(fls3d);
      ((TH1D*)fHistDir->Get("all_flsxy"))->Fill(flsxy);
      ((TH1D*)fHistDir->Get("all_prob"))->Fill(prob);
      ((TH1D*)fHistDir->Get("all_chi2"))->Fill(chi2);
      ((TH1D*)fHistDir->Get("all_alpha"))->Fill(alpha);
      ((TH1D*)fHistDir->Get("all_pt"))->Fill(pt);
      ((TH1D*)fHistDir->Get("all_dr"))->Fill(dr);
      ((TH1D*)fHistDir->Get("all_ptdz"))->Fill(ptdz);
      ((TH1D*)fHistDir->Get("all_ptPis"))->Fill(ptPis);
      ((TH1D*)fHistDir->Get("all_ptPi"))->Fill(ptPi);
      ((TH1D*)fHistDir->Get("all_ptK"))->Fill(ptK);

      if (tm ==1) {  // histogram truth matched candidates
	//if( muid1 || muid2) cout<<" tm "<<muid1<<" "<<muid2;
	((TH1D*)fHistDir->Get("h62"))->Fill( float(qpi+qpis));
	((TH1D*)fHistDir->Get("mds"))->Fill(mdstar);
	((TH1D*)fHistDir->Get("mdz"))->Fill(mdz);
	((TH1D*)fHistDir->Get("dm"))->Fill(dm);
	((TH2D*)fHistDir->Get("h2d"))->Fill(mdz, dm);
	((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
	((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
	((TH1D*)fHistDir->Get("prob"))->Fill(prob);
	((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
	((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
	((TH1D*)fHistDir->Get("pt"))->Fill(pt);
	((TH1D*)fHistDir->Get("dr"))->Fill(dr);

	((TH1D*)fHistDir->Get("ptdz"))->Fill(ptdz);
	((TH1D*)fHistDir->Get("ptPis"))->Fill(ptPis);
	((TH1D*)fHistDir->Get("ptPi"))->Fill(ptPi);
	((TH1D*)fHistDir->Get("ptK"))->Fill(ptK);
      } // it tm


      // Now do the final selection
      bool cut10 = (chi2<2.0) && (pt>6.) && (alpha<0.2) && (fls3d>2.0) && (dr<0.08) && (ptPis>0.5);
      bool cut3 = (dm>=0.135) && (dm<0.155);  // window
      //TCut cut32 = "(fmds-fmdz)>=0.1405&&(fmds-fmdz)<0.1430";  // left 1/2
      //TCut cut33 = "(fmds-fmdz)>=0.1480&&(fmds-fmdz)<0.1505";  // right 1/2
      //TCut cut34 = cut32||cut33;
      //TCut cut35 = cut34||cut31;
      
      // Mass cuts
      bool cut1 = (mdz>1.82) && (mdz<1.91);
      bool cut2 = (mdstar>1.97) && (mdstar<2.05);
      
      //if( muid1 || muid2) cout<<" tm "<<muid1<<" "<<muid2;
      //bool  cut5 = muid1 || muid2; // problem event
      //bool  cut6 = muid1;  // false pi->mu
      //bool  cut7 = muid2;  // false K->mu
      bool  cut8 = !(muid1 || muid2); // clean event
      
      if( cut10 && cut3 && cut1 && cut2 ) {  // final cuts passed 

	((TH1D*)fHistDir->Get("h71"))->Fill(dm);
	if(muid1) ((TH1D*)fHistDir->Get("h72"))->Fill(dm);
	if(muid2) ((TH1D*)fHistDir->Get("h73"))->Fill(dm);
	if(cut8) ((TH1D*)fHistDir->Get("h74"))->Fill(dm);
	if(mutid1) ((TH1D*)fHistDir->Get("h75"))->Fill(dm);
	if(mutid2) ((TH1D*)fHistDir->Get("h76"))->Fill(dm);
	if( !(mutid1||mutid2) ) ((TH1D*)fHistDir->Get("h77"))->Fill(dm);
      } // if final

//       if (0) {
// 	if (dm < 0.147 && dm > 0.145 && !tm) {
// 	  cout << "%%%%%%%%  ??" << endl;
// 	  tm = truthMatch(pCand, 1); 
// 	  cout << "  %%%" << endl;
// 	  dumpHFDstarCand(pCand); 
// 	  cout << "  %%%" << endl;

// 	  for (int i = 0; i < fpEvt->nCands(); ++i) {
// 	    pC1 = fpEvt->getCand(i);
// 	    if (300054 == pC1->fType) {
// 	      dumpHFDstarCand(pC1); 
// 	    }
// 	  }
// 	}
//       }


    } // if select candidate
  }  // candidate loop

  ((TH1D*)fHistDir->Get("all_ncand"))->Fill(ncand);
  if(count>0) ((TH1D*)fHistDir->Get("cands"))->Fill(float(count));
  //  cout << fpCand->fType << " -> mass = " << fpCand->fMass << endl;

  
  if(fcands>0) tree->Fill();
  
  
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
  int daus = (pC->fSig2) - (pC->fSig1);
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
void candAnaDstar::candAnalysis() {

  if (0 == fpCand) return;

  //  candAna::candAnalysis();

}

// ----------------------------------------------------------------------
void candAnaDstar::moreBasicCuts() {
  cout << "   candAnaDstar: more basic cuts" << endl;
}


// ----------------------------------------------------------------------
void candAnaDstar::bookHist() {
  //  candAna::bookHist();
  cout << "==>candAnaDstar: bookHist" << endl;

  fHistDir->cd();

  TH1 *h = new TH1D("mds", "m(dstar)", 70, 1.8, 2.5);
  h = new TH1D("mdz", "m(d0)", 70, 1.8, 2.5);
  h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("dr","dr",100, 0., 1);
  h = new TH1D("pt", "pT", 50, 0., 25);
  h = new TH1D("ptdz", "pT", 50, 0., 25);
  h = new TH1D("ptK",   "pT", 50, 0., 10);
  h = new TH1D("ptPi",  "pT", 50, 0., 10);
  h = new TH1D("ptPis", "pT", 50, 0., 5);
  TH2D *h2 = new TH2D("h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("all_mds", "m(dstar)", 70, 1.8, 2.5);
  h = new TH1D("all_mdz", "m(d0)", 70, 1.8, 2.5);
  h = new TH1D("all_dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("all_ncand", "ncand", 200, 0., 200);
  h = new TH1D("all_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("all_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("all_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("all_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("all_alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("all_dr","dr",100, 0., 1);
  h = new TH1D("all_pt", "pT", 50, 0., 25);
  h = new TH1D("all_ptdz", "pT", 50, 0., 25);
  h = new TH1D("all_ptK",   "pT", 50, 0., 10);
  h = new TH1D("all_ptPi",  "pT", 50, 0., 10);
  h = new TH1D("all_ptPis", "pT", 50, 0., 5);

  h2 = new TH2D("all_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("full_mds", "m(dstar)", 50, 1.8, 2.5);
  h = new TH1D("full_mdz", "m(d0)", 50, 1.8, 2.5);
  h = new TH1D("full_dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("full_ncand", "ncand", 200, 0., 200);
  h = new TH1D("full_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("full_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("full_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("full_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("full_alpha", "alpha", 350, 0., 3.5);
  h = new TH1D("full_dr","dr",100, 0., 1);
  h = new TH1D("full_pt", "pT", 50, 0., 25); 
  h = new TH1D("full_ptdz", "pT", 50, 0., 25);
  h = new TH1D("full_ptK",   "pT", 50, 0., 10);
  h = new TH1D("full_ptPi",  "pT", 50, 0., 10);
  h = new TH1D("full_ptPis", "pT", 50, 0., 5);
  h2 = new TH2D("full_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);


  h = new TH1D("cands", "pT", 100, 0., 100);
  h = new TH1D("pvidx", "PVidx", 100, 0., 100);

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

  h = new TH1D("h11", "h11", 100, 0., 1);
  h = new TH1D("h12", "h12", 100, 0., 1);
  h = new TH1D("h13", "h13", 100, 0., 1);
  h = new TH1D("h14", "h14", 100, 0., 10);
  h = new TH1D("h15", "h15", 100, 0., 1);
  h = new TH1D("h16", "h16", 100, 0., 1);
  h = new TH1D("h17", "h17", 350, 0., 3.5);
  h = new TH1D("h18", "h18", 350, 0., 3.5);
  h = new TH1D("h19", "h19", 350, 0., 3.5);
  h = new TH1D("h20", "h20", 350, 0., 3.5);
  h = new TH1D("h21", "h21", 350, 0., 3.5);
  h = new TH1D("h22", "h22", 350, 0., 3.5);
  h = new TH1D("h23", "h23",   5,-2.5, 2.5);
  h = new TH1D("h24", "h24", 350, 0., 3.5);
  h = new TH1D("h25", "h25", 350, 0., 3.5);
  h = new TH1D("h26", "h26", 100, 0., 1);
  h = new TH1D("h27", "h27", 100, 0., 1);
  h = new TH1D("h28", "h28", 50, 0., 20);
  h = new TH1D("h29", "h29", 50, 0., 20);
  h = new TH1D("h30", "h30", 50, 0., 10);

  h = new TH1D("h31", "h31", 100, 0., 1);
  h = new TH1D("h32", "h32", 100, 0., 1);
  h = new TH1D("h33", "h33", 100, 0., 1);
  h = new TH1D("h34", "h34", 350, 0., 0.35);
  h = new TH1D("h35", "h35", 350, 0., 0.35);
  h = new TH1D("h36", "h36", 350, 0., 0.35);
  h = new TH1D("h37", "h37", 350, 0., 0.35);
  h = new TH1D("h38", "h38", 350, 0., 0.35);

  h = new TH1D("h39", "h39", 50, 0., 10);

  h = new TH1D("h41", "h41",50, 0.0, 20.);
  h = new TH1D("h42", "h42",100, 1.0, 2.5);
  h = new TH1D("h43", "h43",50, 0.0, 20.);
  h = new TH1D("h44", "h44",100, 1.0, 2.5);

  h = new TH1D("h45", "h45",100,  1.,2.5);
  h = new TH1D("h46", "h46",100,  1.,2.5);
  h = new TH1D("h47", "h47",100,  0., 1.0);
  h = new TH1D("h48", "h48",50,  0., 20.);
  h = new TH1D("h49", "h49",200,  0.13, 0.18);
  h = new TH1D("h50", "h50",200,  0.13, 0.18);

  h2 = new TH2D("h51", "h51", 40, 0., 20., 80, 1.5, 2.3);
  h2 = new TH2D("h52", "h52", 35, 0., 3.5, 80, 1.5, 2.3);
  h2 = new TH2D("h53", "h53", 40, 0., 20., 80, 1.5, 2.3);
  h2 = new TH2D("h54", "h53", 40, 0., 10., 80, 1.5, 2.3);
  h2 = new TH2D("h55", "h55", 40, 0., 10., 80, 1.5, 2.3);
  h2 = new TH2D("h56", "h56", 40, 0.,  2., 80, 1.5, 2.3);
  h2 = new TH2D("h57", "h57", 40, 0.,0.4,  80, 1.5, 2.3);

  h = new TH1D("h60", "h60", 5,-2.5,2.5);
  h = new TH1D("h61", "h61", 5,-2.5,2.5);
  h = new TH1D("h62", "h62", 5,-2.5,2.5);
  //h = new TH1D("h63", "h63", 5,-2.5,2.5);
  //h = new TH1D("h64", "h64", 5,-2.5,2.5);
  //h = new TH1D("h65", "h65", 5,-2.5,2.5);
  //h = new TH1D("h66", "h66", 5,-2.5,2.5);
  //h = new TH1D("h67", "h67", 5,-2.5,2.5);

  h = new TH1D("h71", "dm",40,0.135,0.155);
  h = new TH1D("h72", "dm",40,0.135,0.155);
  h = new TH1D("h73", "dm",40,0.135,0.155);
  h = new TH1D("h74", "dm",40,0.135,0.155);
  h = new TH1D("h75", "dm",40,0.135,0.155);
  h = new TH1D("h76", "dm",40,0.135,0.155);
  h = new TH1D("h77", "dm",40,0.135,0.155);

  tree = new TTree("dstar","dstar");
  tree->Branch("fcands",&fcands,"fcands/I");
  tree->Branch("ftm",ftm,"ftm[fcands]/I");
  tree->Branch("fmu",fmu,"fmu[fcands]/I");
  tree->Branch("fmds",fmds,"fmds[fcands]/F");
  tree->Branch("fmdz",fmdz,"fmdz[fcands]/F");

  tree->Branch("ffls3d",ffls3d,"ffls3d[fcands]/F");
  tree->Branch("fchi2",fchi2,"fchi2[fcands]/F");
  tree->Branch("falpha",falpha,"falpha[fcands]/F");
  tree->Branch("fqpis",fqpis,"fqpis[fcands]/F");
  tree->Branch("fdr",fdr,"fdr[fcands]/F");

  tree->Branch("fpt",fpt,"fpt[fcands]/F");
  tree->Branch("fptdz",fptdz,"fptdz[fcands]/F");
  tree->Branch("fptpis",fptpis,"fptpis[fcands]/F");
  tree->Branch("fptpi",fptpi,"fptpi[fcands]/F");
  tree->Branch("fptk",fptk,"fptk[fcands]/F");


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
  int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

   
  }

}
