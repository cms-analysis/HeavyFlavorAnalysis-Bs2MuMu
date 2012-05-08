#include "candAnaHh.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"
#include "danekUtils.h"

using namespace std;

namespace {
  TVector3 BdVertexGen(0,0,0), PVGen(0,0,0);  
  TVector3 BdMomGen(0,0,0), Pi1MomGen(0,0,0), Pi2MomGen(0,0,0);  
}

// ----------------------------------------------------------------------
candAnaHh::candAnaHh(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaHh: name = " << name << ", reading cutsfile " << cutsFile << endl;

  readCuts(cutsFile, 1); 

}


// ----------------------------------------------------------------------
candAnaHh::~candAnaHh() {
  cout << "==> candAnaHh: destructor..." << endl;
  tree->Write();
}


// ----------------------------------------------------------------------
// To analyze the MC event
bool candAnaHh::anaMC(TAna01Event *evt) {
  const bool print = false;

  fpEvt = evt; 
  bool foundPV=false, foundBd = false, foundPi1 = false, foundPi2=false;
  int pC0 = 0;

  
  int numGenCands = fpEvt->nGenCands();
  if(print) cout << "Found " << numGenCands << " gen cands in event" << endl;
  for (int it = 0; it < numGenCands; ++it) {  // loop over all gen candidates 

    foundBd = false, foundPi1 = false, foundPi2=false;

    TGenCand * pCand = fpEvt->getGenCand(it);

    if( !foundPV && (abs(pCand->fID) == 5) ) {foundPV=true; PVGen=pCand->fV; if(print) cout<<" PV "<<PVGen.Z()<<endl;}   // get PV 

    if( ( abs(pCand->fID) != 511) ) continue;        // skip others 

    if(print) cout <<" Bd "<<it<< " " << pCand->fNumber << " "<<pCand->fID<<" "<<pCand->fQ<<" "<<pCand->fStatus<<" "
		   <<pCand->fMom1<<" "<<pCand->fMom2<<" "<<pCand->fDau1<<" "
		   <<pCand->fDau2<<" "<<pCand->fP.Perp()<<" "<<pCand->fV.Z()<<endl;

    BdMomGen = (pCand->fP.Vect());
    foundBd = true;
    pC0 = it;

    // Look at daugthers
    int i1 = (pCand->fDau2)-(pCand->fDau1)+1;
    int i2=0;
    //if(i1!=2) {continue;} // fpEvt->dumpGenBlock();}
    if(i1!=2) { if(print) cout<<" number of daughters wrong skip "<<i1<<endl; continue;} // fpEvt->dumpGenBlock();}

    for(int id=(pCand->fDau1);id<=(pCand->fDau2);++id) {
      TGenCand * dau = fpEvt->getGenCand(id);  // check daughters

      if( abs(dau->fID) == 211 ) { //  pions
	if(print) cout <<" Pi "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;

	BdVertexGen = (dau->fV);  // Bd decay vertex

	i2++;
	if(i2==1)      {foundPi1=true; Pi1MomGen = (dau->fP.Vect());}
	else if(i2==2) {foundPi2=true; Pi2MomGen = (dau->fP.Vect());}

      }      

    } // daugther loop 

      if(foundBd && foundPi1 && foundPi2) break; 

  } // gen part loop 


  bool ok = foundPV && foundBd && foundPi1 && foundPi2;
  if(ok) {

    TVector3 t1(BdVertexGen-PVGen);
    double a1 = t1.Angle(BdMomGen);  // Bd pointing angle

    ((TH1D*)fHistDir->Get("h11"))->Fill(t1.Mag());
    ((TH1D*)fHistDir->Get("h20"))->Fill(a1);

    ((TH1D*)fHistDir->Get("h28"))->Fill(BdMomGen.Perp());
    ((TH1D*)fHistDir->Get("h29"))->Fill(Pi1MomGen.Perp());
    ((TH1D*)fHistDir->Get("h30"))->Fill(Pi2MomGen.Perp());

    double dr    = danekUtils::twoBodyDecayAngle(Pi1MomGen, Pi2MomGen);
    //double pt    = danekUtils::twoBodyDecayMomPerp(Pi1MomGen, Pi2MomGen);
    //double m1    = danekUtils::twoBodyDecayMass(Pi1MomGen, Pi2MomGen, MPION, MPION);
    //double m2    = danekUtils::twoBodyDecayMass(Pi1MomGen, Pi2MomGen, MMUON, MMUON);

    ((TH1D*)fHistDir->Get("h47"))->Fill(dr); 
    //((TH1D*)fHistDir->Get("h41"))->Fill(pt); 
    //((TH1D*)fHistDir->Get("h42"))->Fill(m2); 
    //((TH1D*)fHistDir->Get("h44"))->Fill(m1); 
  } else {
    if(print) {cout<<" Bd->pipi not found in Gen"<<endl; fpEvt->dumpGenBlock(); }

  } // if OK

  return ok;
}
//------------------------------------------------------------------------------------------
void candAnaHh::evtAnalysis(TAna01Event *evt) {
  fpEvt = evt; 
  fcands=0;
  int count = 0;

  TAnaCand *pCand(0);
  TAnaTrack *pPi1, *pPi2; 
  double fl3d(0), fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), dr(0.); 
  int tm=0; 
  int ncand(0); 

  if(fVerbose>10) 
    cout << "Evt: " << fEvt << " ----------------------------------------------------------------------" << endl;


  // -- loop over all seq vtx fit candidates for D*
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {  // just count candidates
    pCand = fpEvt->getCand(iC);
    if (211211 == pCand->fType || pCand->fType == 91 || pCand->fType == -91) {

      if(pCand->fType == 211211) ++ncand;
      if (fIsMC>0 && fVerbose>0) {
	TVector3 s = pCand->fVtx.fPoint;
	cout << " -> " << iC <<" "<< pCand->fType<<" "<<pCand->fMom;
	cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	cout << " mass: " << pCand->fMass << " " << s.Z() <<" "<<pCand->fPlab.Mag()<<endl;
	cout << "DUMP HFDHh with mass = " << pCand->fMass << endl;
	dumpHFTruthCand(pCand); 
      }
    }  
  }
  

  if(fVerbose>0) cout<<" num of cands "<<ncand<<" "<<fVerbose<<" "<<fIsMC<<endl;
  ((TH1D*)fHistDir->Get("all_ncand"))->Fill(ncand);

  ((TH1D*)fHistDir->Get("status"))->Fill(0.);

  // Check MC Gen 
  bool ok = false;
  if(fIsMC) {
      ok = anaMC(evt);
      if(fVerbose>10) cout<<" Correct candidate = "<<ok<<endl;
      if(ok) ((TH1D*)fHistDir->Get("status"))->Fill(1.);
  }

  if(ncand>0) ((TH1D*)fHistDir->Get("status"))->Fill(2.);

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (211211 != pCand->fType) continue;  //select HH
    ((TH1D*)fHistDir->Get("status"))->Fill(10.);

    if(fVerbose>0) {cout<<" Dump candidate "<<endl; dumpHFHhCand(pCand); pCand->dump();}
    
    double candMass = pCand->fMass;
    TVector3 pBdMom  = pCand->fPlab;
    double candPt  = pBdMom.Perp();
    double candEta = pBdMom.Eta();
    double doca = pCand->fMaxDoca; // doca between the 2 pions

    // PV 
    int pvidx = (pCand->fPvIdx > -1? pCand->fPvIdx : 0);  // PV index
    TAnaVertex *pv =  fpEvt->getPV(pvidx);
    TVector3 pvPos =  pv->fPoint;
    double pvNtrk = pv->getNtracks();
    double pvNdof = pv->fNdof;
    double pvAveW = ((pvNdof+2.)/2.)/pvNtrk;
    if(fVerbose>10) cout<<" PV "<<pvPos.Z()<<" "<<pvNtrk<<" "<<pvNdof<<" "<<pvAveW<<" "<<pv->fChi2<<" "<<pv->fProb<<endl;

    // HH decay vertex
    TAnaVertex sv = pCand->fVtx;  //  HH vertex
    TVector3 svPos = sv.fPoint;

    TVector3 svpv(svPos-pvPos);

    fl3d = sv.fD3d; 
    fls3d = sv.fD3d/sv.fD3dE; 
    flsxy = sv.fDxy/sv.fDxyE; 
    prob  = sv.fProb;
    chi2  = sv.fChi2;

    alpha = svpv.Angle(pCand->fPlab);  // pointing angle
    //dr = piSlowMom.Angle(pCand->fPlab); // pislow openinig

    // Check that only 2 tracks come from teh candidate 
    int indx1 = (pCand->fSig1); // pion1
    int indx2 = (pCand->fSig2); // pion2

    if(fVerbose>10) cout<<" signal tracks "<<(indx2-indx1)<<" "<<indx1<<" "<<indx2<<endl;
    if(indx1<0 || indx2<0 || (indx2-indx1+1)!=2 ) {
      if(fVerbose>0) cout << " Wrong number of candidate isgnal tracks " << indx1<<" "<<indx2 << endl; 
      continue;
    }
    ((TH1D*)fHistDir->Get("status"))->Fill(11.);

    // Get pion tracks
    int pi1Id = fpEvt->getSigTrack(indx1)->fIndex;
    int pi2Id = fpEvt->getSigTrack(indx2)->fIndex;
    pPi1 = fpEvt->getRecTrack( pi1Id ); // track 1 
    pPi2 = fpEvt->getRecTrack( pi2Id ); // track 2
    //pPi1->dump(); 
    //pPi2->dump(); 
 
    TVector3 pi1Mom = pPi1->fPlab; // pi 1 momentum vector 
    TVector3 pi2Mom = pPi2->fPlab; // pi 2 momentum vector 
    double pt1 = pi1Mom.Perp();
    double pt2 = pi2Mom.Perp();

    if(fVerbose>10) {
      cout<<" pion1 "<<pi1Id<<" "<< pt1<<" "<<pPi1->fQ<<" "<<fpEvt->getSigTrack(indx1)->fMCID <<endl; 
      cout<<" pion2 "<<pi2Id<<" "<< pt2<<" "<<pPi2->fQ<<" "<<fpEvt->getSigTrack(indx2)->fMCID <<endl; 
    }

    // skip candidates with the same charge pions
    if(pPi1->fQ == pPi2->fQ) {continue;}
    ((TH1D*)fHistDir->Get("status"))->Fill(12.);

    // Get the di-pion opening angle
    dr = pi1Mom.Angle(pi2Mom);

    // 

    // Histogram 
    ((TH1D*)fHistDir->Get("all_fl3d"))->Fill(fl3d);
    ((TH1D*)fHistDir->Get("all_fls3d"))->Fill(fls3d);
    ((TH1D*)fHistDir->Get("all_flsxy"))->Fill(flsxy);
    ((TH1D*)fHistDir->Get("all_prob"))->Fill(prob);
    ((TH1D*)fHistDir->Get("all_chi2"))->Fill(chi2);
    ((TH1D*)fHistDir->Get("all_alpha"))->Fill(alpha);
    ((TH1D*)fHistDir->Get("all_pt"))->Fill(candPt);
    ((TH1D*)fHistDir->Get("all_m"))->Fill(candMass);
    ((TH1D*)fHistDir->Get("all_eta"))->Fill(candEta);
    ((TH1D*)fHistDir->Get("all_dr"))->Fill(dr);
    ((TH1D*)fHistDir->Get("all_ptPi"))->Fill(pt1);
    ((TH1D*)fHistDir->Get("all_ptPi"))->Fill(pt2);
    ((TH1D*)fHistDir->Get("all_doca"))->Fill(doca);
    ((TH1D*)fHistDir->Get("all_pvW"))->Fill(pvAveW);
    ((TH1D*)fHistDir->Get("all_pvid"))->Fill(float(pvidx));
    ((TH2D*)fHistDir->Get("h2d"))->Fill(candPt,candMass);


    // truthMatch return an interger, 0-no match, 1-correct match,
    tm = 0; 
    if(fIsMC) {

      tm = truthMatch(pCand,fVerbose); // check truth matching
      if(fVerbose>10) cout<<" Truth matching = "<<tm<<endl;

      if ( (tm == 1) && fVerbose>20 ) {
	  cout << " Truth matched cand -> " << pCand->fType;
	  cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	  cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	  cout << " mass: " << pCand->fMass << " " << tm<<endl;
	  cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
	  dumpHFHhCand(pCand); 
	  fpEvt->dumpGenBlock(); 
      } // end if
    } // if MC
       
    if(!fIsMC || tm==1) {  // Do for all data and for matched MC
      
      ((TH1D*)fHistDir->Get("status"))->Fill(21.);
      if(ok) ((TH1D*)fHistDir->Get("status"))->Fill(22.);

      // Histogram 
      ((TH1D*)fHistDir->Get("mc_fl3d"))->Fill(fl3d);
      ((TH1D*)fHistDir->Get("mc_fls3d"))->Fill(fls3d);
      ((TH1D*)fHistDir->Get("mc_flsxy"))->Fill(flsxy);
      ((TH1D*)fHistDir->Get("mc_prob"))->Fill(prob);
      ((TH1D*)fHistDir->Get("mc_chi2"))->Fill(chi2);
      ((TH1D*)fHistDir->Get("mc_alpha"))->Fill(alpha);
      ((TH1D*)fHistDir->Get("mc_pt"))->Fill(candPt);
      ((TH1D*)fHistDir->Get("mc_m"))->Fill(candMass);
      ((TH1D*)fHistDir->Get("mc_eta"))->Fill(candEta);
      ((TH1D*)fHistDir->Get("mc_dr"))->Fill(dr);
      ((TH1D*)fHistDir->Get("mc_ptPi"))->Fill(pt1);
      ((TH1D*)fHistDir->Get("mc_ptPi"))->Fill(pt2);
      ((TH1D*)fHistDir->Get("mc_doca"))->Fill(doca);
      ((TH1D*)fHistDir->Get("mc_pvW"))->Fill(pvAveW);
      ((TH1D*)fHistDir->Get("mc_pvid"))->Fill(float(pvidx));
    
    }  // if 

    if(ok && tm==1) {  // Do here comparison between RECO and GEN quantities
      double c1 = (PVGen-pvPos).Mag();
      double c2 = (BdVertexGen-svPos).Mag();

      //cout<<" RECO-MV vertex distance "<<c1<<" "<<c2<<" "<<c3<<endl;
      ((TH1D*)fHistDir->Get("h31"))->Fill(c1);
      ((TH1D*)fHistDir->Get("h32"))->Fill(c2);
      
      double a11 = Pi1MomGen.Angle(pi1Mom);
      double a12 = Pi2MomGen.Angle(pi2Mom);
      double a13 = BdMomGen.Angle(pBdMom);
      
      ((TH1D*)fHistDir->Get("h35"))->Fill(a11);
      ((TH1D*)fHistDir->Get("h36"))->Fill(a12);
      ((TH1D*)fHistDir->Get("h37"))->Fill(a13);
      
    } // if OK

    // Do the selection cuts
    if (doca > 0.025) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(13.); 
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(23.);

    if (pvAveW < 0.70) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(14.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(24.);

    if(pt1<4. || pt2<4.) continue;  //  does not do aything for data, cand reco is already with 4
    ((TH1D*)fHistDir->Get("status"))->Fill(15.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(25.);

    //if (dr > 0.25) continue;
    //if (chi2 > 2.0) continue;
    //if (pt < 2) continue;
    //if (alpha > 0.4) continue;
    //if (fls3d < 2) continue;
    
    count++; // count selected candidates 
    ((TH1D*)fHistDir->Get("status"))->Fill(20.);

    // Look at muid
    //bool muid1 = goodMuon(pPi);  // true for good  muons 
    //bool muid2 = goodMuon(pK);
    //bool mutid1 = tightMuon(pPi);  // true for good/tight  muons 
    //bool mutid2 = tightMuon(pK);
    
    //if( muid1 || muid2) cout<<"fake muon "<<muid1<<" "<<muid2;
    
    // Final histos after cuts
    if ( !fIsMC || tm ==1 ) {  // histogram truth matched candidates
      
      // Histogram 
      ((TH1D*)fHistDir->Get("fl3d"))->Fill(fl3d);
      ((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
      ((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
      ((TH1D*)fHistDir->Get("prob"))->Fill(prob);
      ((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
      ((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
      ((TH1D*)fHistDir->Get("pt"))->Fill(candPt);
      ((TH1D*)fHistDir->Get("m"))->Fill(candMass);
      ((TH1D*)fHistDir->Get("eta"))->Fill(candEta);
      ((TH1D*)fHistDir->Get("dr"))->Fill(dr);
      ((TH1D*)fHistDir->Get("ptPi"))->Fill(pt1);
      ((TH1D*)fHistDir->Get("ptPi"))->Fill(pt2);
      ((TH1D*)fHistDir->Get("doca"))->Fill(doca);
      ((TH1D*)fHistDir->Get("pvW"))->Fill(pvAveW);
      ((TH1D*)fHistDir->Get("pvid"))->Fill(float(pvidx));
      
    } // it tm
    
      // Save in a tree, save only masses between 130-160MeV
    if( (fcands<10) ) {
      ftm[fcands] = tm;
      
      fm[fcands]=candMass;
      fpt[fcands]=candPt;
      fchi2[fcands]=chi2;
      falpha[fcands]=alpha;
      ffls3d[fcands]=fls3d;
      fdr[fcands]=dr;
      fptpi1[fcands]=pt1;
      fptpi2[fcands]=pt2;
      fdoca[fcands]=doca;
      fweight[fcands]=pvAveW;
      
      fcands++;
    } // if fcands
    
  }  // candidate loop

  if(count>0) {
    ((TH1D*)fHistDir->Get("cands"))->Fill(float(count));
    ((TH1D*)fHistDir->Get("status"))->Fill(3.);
    if(fcands>0) tree->Fill();
  }

  
    
}

// ----------------------------------------------------------------------
void candAnaHh::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 
  if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    pT->dump(); 
  }
}


// ----------------------------------------------------------------------
void candAnaHh::dumpHFHhCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 

  cout << "HFHhCand: idx = " << pC->fIndex << " type = " << pC->fType<< " m = " << pC->fMass <<endl;
  cout<<"DOCA "<<pC->fMinDoca<<" "<<pC->fMaxDoca<<endl;


//   int nsize = pC->fNstTracks.size(); 
//   cout<<" tracks "<<nsize<<endl;
//   if (nsize>0) {
//     for(int i = 0; i<nsize; ++i) {
//       int trkId = pC->fNstTracks[i].first;
//       double doca = pC->fNstTracks[i].second.first;
//       cout<<i<<" "<<trkId<<" "<<doca<<endl;
//     }
//   }

  // -- D0 daughters
  if ( pC->fSig1<0  || pC->fSig2<0 ) {
    cout << "XXXXXXXXX cannot get signal cand of " << pC->fType << endl;
    return;
  }

  // -- pion 1
  int indx = (pC->fSig1); 
  pT = fpEvt->getRecTrack( fpEvt->getSigTrack(indx)->fIndex ); 
  cout << fpEvt->getSigTrack(indx)->fMCID << " " ; 
  pT->dump(); 

  // -- pion 2
  indx = (pC->fSig2); 
  pT = fpEvt->getRecTrack( fpEvt->getSigTrack(indx)->fIndex ); 
  cout << fpEvt->getSigTrack(indx)->fMCID << " " ; 
  pT->dump(); 



}


// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
int candAnaHh::truthMatch(TAnaCand *pCand, int verbose) {

  // -- check pion1
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
  if(fVerbose>0) cout<< pG->fID<<" "<<endl;  // gen id
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return 0;
  }

  // Check the mother
  int mom = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 511 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 511 != pG->fID, pG->fID = " << pG->fID  << endl;
    return 0;
  }


  // -- check pion2
  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex);
  if(fVerbose>0) cout<<fpEvt->getSigTrack(pCand->fSig2)->fIndex<<" "
		     <<pT->fGenIndex<<" "
		     <<(fpEvt->getRecTrack(pT->fIndex)->fGenIndex)<<" "<<verbose<<endl;  // same as above 
  
  if (pT->fGenIndex < 0) {
    if (verbose > 0) cout << "pT->fGenIndex < 0" << endl;
    return 0; 
  }

  pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex); 
  if (0 == pG) {
    if (verbose > 0) cout << "0 == pG" << endl;
    return 0;
  }
  if(fVerbose>0) cout<< pG->fID<<" "<<endl;  // gen id
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return 0;
  }

  // Check the mother
  mom = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 511 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 511 != pG->fID, pG->fID = " << pG->fID  << endl;
    return 0;
  }

  
  // -- check daughters
  int type(0), moIdx(-1); 
  int count=0, countPi=0, countK=0, missK=0, missPi=0;
  //int daus = (pC->fSig2) - (pC->fSig1);
  //cout<<" Bd daugthers "<<(daus+1)<<endl;
  for (int id = pCand->fSig1; id <= pCand->fSig2; ++id) {
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

    if( TMath::Abs(pG->fID) == 211 ) { // Pion
      countPi++;
      //cout <<count<<" "<<countPi; 
      // The check below does not work because we assing the muon mass to all tracks
      //if( TMath::Abs(type) != 211 ) {
      //cout << " Pion identified as ?, type = " << type << " pG->fID = " << pG->fID 
      //     << " track " << pT->fIndex;
      //missPi++;
      //return 0;
      //}
      //cout<<endl;
    }

  }

  if     (countPi==2) return  1; // select righ combination 

  return 0; 
}
  

// ----------------------------------------------------------------------
void candAnaHh::candAnalysis() {

  if (0 == fpCand) return;

  //  candAna::candAnalysis();

}

// ----------------------------------------------------------------------
void candAnaHh::moreBasicCuts() {
  cout << "   candAnaHh: more basic cuts" << endl;
}


// ----------------------------------------------------------------------
void candAnaHh::bookHist() {
  //  candAna::bookHist();
  cout << "==>candAnaHh: bookHist" << endl;

  fHistDir->cd();

  TH1 *h = new TH1D("status", "status", 100, -0.5, 99.5);
//   h = new TH1D("mdz", "m(d0)", 70, 1.8, 2.5);
//   h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);
//   h = new TH1D("ncand", "ncand", 200, 0., 200);
//   h = new TH1D("fls3d", "fls3d", 60, 0., 20);
//   h = new TH1D("flsxy", "flsxy", 60, 0., 20);
//   h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
//   h = new TH1D("chi2", "chi2", 100, 0., 10.);
//   h = new TH1D("alpha", "alpha", 50, 0., 1.0);
//   h = new TH1D("dr","dr",100, 0., 1);
//   h = new TH1D("pt", "pT", 50, 0., 25);
//   h = new TH1D("ptdz", "pT", 50, 0., 25);
//   h = new TH1D("ptK",   "pT", 50, 0., 10);
//   h = new TH1D("ptPi",  "pT", 50, 0., 10);
//   h = new TH1D("ptPis", "pT", 50, 0., 5);
  TH2D *h2 = new TH2D("h2d", "m vs pt", 50, 0., 25, 120, 0., 12.);

  h = new TH1D("all_m", "cand mass", 100, 0.0, 10.);
  h = new TH1D("all_ncand", "ncand", 200, 0., 200);
  h = new TH1D("all_fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("all_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("all_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("all_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("all_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("all_alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("all_dr","dr",100, 0., 3.);
  h = new TH1D("all_pt", "cand pT", 100, 0., 50.);
  h = new TH1D("all_ptPi","pions pt", 100, 0., 20);
  h = new TH1D("all_doca","pions doca", 100, 0., 0.2);
  h = new TH1D("all_pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("all_eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("all_pvid", "cand PVidx", 100, 0., 100);

  h = new TH1D("mc_m", "cand mass", 150, 0.0, 15.);
  h = new TH1D("mc_ncand", "ncand", 200, 0., 200);
  h = new TH1D("mc_fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("mc_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("mc_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("mc_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("mc_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("mc_alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("mc_dr","dr",100, 0., 3.);
  h = new TH1D("mc_pt", "cand pT", 100, 0., 50.);
  h = new TH1D("mc_ptPi",  "pions pt", 100, 0., 20);
  h = new TH1D("mc_doca","pions doca", 100, 0., 0.2);
  h = new TH1D("mc_pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("mc_eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("mc_pvid", "cand PVidx", 100, 0., 100);

  h = new TH1D("m", "cand mass", 150, 0.0, 15.);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("dr","dr",100, 0., 3.);
  h = new TH1D("pt", "cand pT", 100, 0., 50.);
  h = new TH1D("ptPi",  "pions pt", 100, 0., 20);
  h = new TH1D("doca","pions doca", 100, 0., 0.2);
  h = new TH1D("pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("pvid", "cand PVidx", 100, 0., 100);

  //h2 = new TH2D("all_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("cands", "pT", 100, 0., 100);

//   h = new TH1D("h1", "h1", 100, 0., 1);
//   h = new TH1D("h2", "h2", 100, 0., 1);
//   h = new TH1D("h3", "h3", 100, 0., 1);
//   h = new TH1D("h4", "h4", 350, 0., 3.5);
//   h = new TH1D("h5", "h5", 100, 0., 1);
//   h = new TH1D("h6", "h6", 100, 0., 1);
//   h = new TH1D("h7", "h7", 350, 0., 3.5);
//   h = new TH1D("h8", "h8", 350, 0., 3.5);
//   h = new TH1D("h9", "h9", 350, 0., 3.5);
//   h = new TH1D("h10","h10",100, 0., 1);

  h = new TH1D("h11", "h11", 100, 0., 1);
//   h = new TH1D("h12", "h12", 100, 0., 1);
//   h = new TH1D("h13", "h13", 100, 0., 1);
//   h = new TH1D("h14", "h14", 100, 0., 10);
//   h = new TH1D("h15", "h15", 100, 0., 1);
//   h = new TH1D("h16", "h16", 100, 0., 1);
//   h = new TH1D("h17", "h17", 350, 0., 3.5);
//   h = new TH1D("h18", "h18", 350, 0., 3.5);
//   h = new TH1D("h19", "h19", 350, 0., 3.5);
  h = new TH1D("h20", "h20", 100, 0.,0.2);
//   h = new TH1D("h21", "h21", 350, 0., 3.5);
//   h = new TH1D("h22", "h22", 350, 0., 3.5);
//   h = new TH1D("h23", "h23",   5,-2.5, 2.5);
//   h = new TH1D("h24", "h24", 350, 0., 3.5);
//   h = new TH1D("h25", "h25", 350, 0., 3.5);
//   h = new TH1D("h26", "h26", 100, 0., 1);
//   h = new TH1D("h27", "h27", 100, 0., 1);
  h = new TH1D("h28", "h28", 80, 0., 40);
  h = new TH1D("h29", "h29", 80, 0., 40);
  h = new TH1D("h30", "h30", 80, 0., 40);

  h = new TH1D("h31", "h31", 100, 0., 1);
  h = new TH1D("h32", "h32", 100, 0., 1);
//   h = new TH1D("h33", "h33", 100, 0., 1);
//   h = new TH1D("h34", "h34", 350, 0., 0.35);
  h = new TH1D("h35", "h35", 350, 0., 0.35);
  h = new TH1D("h36", "h36", 350, 0., 0.35);
  h = new TH1D("h37", "h37", 350, 0., 0.35);
//   h = new TH1D("h38", "h38", 350, 0., 0.35);

//   h = new TH1D("h39", "h39", 50, 0., 10);

  h = new TH1D("h41", "h41",80, 0.0, 40.);
  h = new TH1D("h42", "h42",100, 0.0, 20.);
//   h = new TH1D("h43", "h43",50, 0.0, 20.);
  h = new TH1D("h44", "h44",100, 0.0, 20.);

//   h = new TH1D("h45", "h45",100,  1.,2.5);
//   h = new TH1D("h46", "h46",100,  1.,2.5);
  h = new TH1D("h47", "h47",100,  0., 3.0);
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
//   h = new TH1D("h61", "h61", 5,-2.5,2.5);
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

  tree = new TTree("hh","hh");
  tree->Branch("fcands",&fcands,"fcands/I");
  tree->Branch("ftm",ftm,"ftm[fcands]/I");
  tree->Branch("fm",fm,"fm[fcands]/F");
  tree->Branch("ffls3d",ffls3d,"ffls3d[fcands]/F");
  tree->Branch("fchi2",fchi2,"fchi2[fcands]/F");
  tree->Branch("falpha",falpha,"falpha[fcands]/F");
  tree->Branch("fdr",fdr,"fdr[fcands]/F");
  tree->Branch("fpt",fpt,"fpt[fcands]/F");
  tree->Branch("fptpi1",fptpi1,"fptpi1[fcands]/F");
  tree->Branch("fptpi2",fptpi2,"fptpi2[fcands]/F");
  tree->Branch("fdoca",fdoca,"fdoca[fcands]/F");
  tree->Branch("fweight",fweight,"fweight[fcands]/F");


}


// ----------------------------------------------------------------------
void candAnaHh::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaHh: Reading " << fCutFile << " for cut settings" << endl;
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
