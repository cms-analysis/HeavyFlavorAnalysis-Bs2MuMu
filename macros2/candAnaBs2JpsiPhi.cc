#include "candAnaBs2JpsiPhi.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::candAnaBs2JpsiPhi(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaBs2JpsiPhi: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::~candAnaBs2JpsiPhi() {
  cout << "==> candAnaBs2JpsiPhi: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::candAnalysis() {

  if (0 == fpCand) return;

  // -- Check for J/psi mass (but see below for truth candidates)
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    //    cout << "i = " << i << " pD = " << pD << endl;
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp();
      fJpsiEta  = pD->fPlab.Eta();
      fJpsiPhi  = pD->fPlab.Phi(); 
      //       cout << "type = " << pD->fType 
      //        	   << " with mass = " << pD->fMass 
      // 	   << " fGoodJpsiMass = " << fGoodJpsiMass 
      // 	   << endl;
    }
  }

  // -- special case for truth candidates (which have no daughter cands)
  if (fpCand->fType > 999999) {
    TAnaTrack *p0; 
    TAnaTrack *p1(0);
    TAnaTrack *p2(0); 
    
    for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
      p0 = fpEvt->getSigTrack(it);     
      if (TMath::Abs(p0->fMCID) != 13) continue;
      if (0 == p1) {
	p1 = p0; 
      } else {
	p2 = p0; 
      }
    }
    
    if (0 == p1) {
      cout << "candAnaBs2JpsiPhi::candAnalysis  no muon 1 found " << endl;
      return; 
    }
    if (0 == p2) {
      cout << "candAnaBs2JpsiPhi::candAnalysis  no muon 2 found " << endl;
      return; 
    }

    TLorentzVector mu1, mu2; 
    mu1.SetPtEtaPhiM(p1->fPlab.Perp(), p1->fPlab.Eta(), p1->fPlab.Phi(), MMUON); 
    mu2.SetPtEtaPhiM(p2->fPlab.Perp(), p2->fPlab.Eta(), p2->fPlab.Phi(), MMUON); 
    
    TLorentzVector psi = mu1 + mu2; 
    if ((JPSIMASSLO < psi.M()) && (psi.M() < JPSIMASSHI)) fGoodJpsiMass = true;
    fJpsiMass = psi.M();
    fJpsiPt   = psi.Pt();
    fJpsiEta  = psi.Eta();
    fJpsiPhi  = psi.Phi(); 
  }    

  // -- Get Kaons
  TAnaTrack *p0; 
  TAnaTrack *p1(0), *ps1(0);
  TAnaTrack *p2(0), *ps2(0); 
  
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (TMath::Abs(p0->fMCID) != 321) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }

  if (0 == p1) {
    cout << "candAnaBs2JpsiPhi::candAnalysis  no kaon 1 found " << endl;
    return; 
  }
  if (0 == p2) {
    cout << "candAnaBs2JpsiPhi::candAnalysis  no kaon 2 found " << endl;
    return; 
  }
 
  // -- switch to RecTracks!
  ps1 = p1; 
  p1 = fpEvt->getRecTrack(p1->fIndex);
  ps2 = p2; 
  p2 = fpEvt->getRecTrack(p2->fIndex);


  fKa1Pt        = p1->fPlab.Perp(); 
  fKa1Eta       = p1->fPlab.Eta(); 
  fKa1Phi       = p1->fPlab.Phi(); 
  fKa1TkQuality = p1->fTrackQuality & TRACKQUALITY;
  fKa1PtNrf     = ps1->fPlab.Perp();
  fKa1EtaNrf    = ps1->fPlab.Eta();

  fKa2Pt        = p2->fPlab.Perp(); 
  fKa2Eta       = p2->fPlab.Eta(); 
  fKa2Phi       = p2->fPlab.Phi(); 
  fKa2TkQuality = p2->fTrackQuality & TRACKQUALITY;
  fKa2PtNrf     = ps2->fPlab.Perp();
  fKa2EtaNrf    = ps2->fPlab.Eta();

// FIXME
//   if (tmCand(fpCand)) {
//     TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(p1->fIndex)->fGenIndex);
//     fKa1PtGen     = pg1->fP.Perp();
//     fKa1EtaGen    = pg1->fP.Eta();
//     TGenCand *pg2 = fpEvt->getGenCand(fpEvt->getRecTrack(p2->fIndex)->fGenIndex);
//     fKa2PtGen     = pg2->fP.Perp();
//     fKa2EtaGen    = pg2->fP.Eta();
//   } else {
//     fKa1PtGen     = -99.;
//     fKa1EtaGen    = -99.;
//     fKa2PtGen     = -99.;
//     fKa2EtaGen    = -99.;
//   }

  fDeltaR  = p1->fPlab.DeltaR(p2->fPlab); 

  TLorentzVector ka1, ka2; 
  ka1.SetPtEtaPhiM(fKa1Pt, fKa1Eta, fKa1Phi, MKAON); 
  ka2.SetPtEtaPhiM(fKa2Pt, fKa2Eta, fKa2Phi, MKAON); 

  TLorentzVector phiCand = ka1 + ka2; 
  fMKK     = phiCand.M();
  fPhiPt   = phiCand.Pt();
  fPhiEta  = phiCand.Eta();
  fPhiPhi  = phiCand.Phi();

  fGoodDeltaR = (fDeltaR < DELTAR);
  fGoodMKK    = ((MKKLO < fMKK ) && (fMKK < MKKHI)); 
  
  fPreselection = fPreselection && fGoodJpsiMass && fGoodMKK && fGoodDeltaR; 


  candAna::candAnalysis();
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(3); 

  if (fIsMC) {
    fTree->Fill(); 
  } else {
    if (fPreselection) {
      ((TH1D*)fHistDir->Get("../monEvents"))->Fill(13); 
      fTree->Fill(); 
    }         
  }

}

// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::moreBasicCuts() {
  cout << "   candAnaBs2JpsiPhi: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
  fAnaCuts.addCut("fGoodDeltaR", "Delta R(KK)", fGoodDeltaR); 
  fAnaCuts.addCut("fGoodMKK", "m(KK) [GeV]", fGoodMKK); 
}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::processType() {

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::genMatch() {

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::candMatch() {

}




// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::bookHist() {
  candAna::bookHist();
  cout << "==>candAnaBs2JpsiPhi: bookHist" << endl;

  // -- Additional reduced tree variables
  fTree->Branch("mpsi",  &fJpsiMass, "mpsi/D");
  fTree->Branch("psipt", &fJpsiPt,   "psipt/D");
  fTree->Branch("psieta",&fJpsiEta,  "psieta/D");
  fTree->Branch("psiphi",&fJpsiPhi,  "psiphi/D");
  fTree->Branch("mkk",   &fMKK,      "mkk/D");
  fTree->Branch("phipt", &fPhiPt,    "phipt/D");
  fTree->Branch("phieta",&fPhiEta,   "phieta/D");
  fTree->Branch("phiphi",&fPhiPhi,   "phiphi/D");
  fTree->Branch("dr",    &fDeltaR,   "dr/D");

  fTree->Branch("k1pt",  &fKa1Pt,    "k1pt/D");
  fTree->Branch("k1eta", &fKa1Eta,   "k1eta/D");
  fTree->Branch("k1phi", &fKa1Phi,   "k1phi/D");
  fTree->Branch("k1gt",  &fKa1TkQuality,"k1gt/I");
  fTree->Branch("k2pt",  &fKa2Pt,    "k2pt/D");
  fTree->Branch("k2eta", &fKa2Eta,   "k2eta/D");
  fTree->Branch("k2phi", &fKa2Phi,   "k2phi/D");
  fTree->Branch("k2gt",  &fKa2TkQuality,"k2gt/I");

  fTree->Branch("t3pt",  &fKa1PtNrf, "t3pt/D");
  fTree->Branch("t3eta", &fKa1EtaNrf,"t3eta/D");

  fTree->Branch("t4pt",  &fKa2PtNrf, "t4pt/D");
  fTree->Branch("t4eta", &fKa2EtaNrf,"t4eta/D");

  fTree->Branch("g3pt", &fKa1PtGen,  "g3pt/D");
  fTree->Branch("g3eta",&fKa1EtaGen, "g3eta/D");
  fTree->Branch("g4pt", &fKa2PtGen,  "g4pt/D");
  fTree->Branch("g4eta",&fKa2EtaGen, "g4eta/D");

}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaBs2JpsiPhi: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;
  int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "JPSITYPE")) {
      JPSITYPE = CutValue; ok = 1;
      if (dump) cout << "JPSITYPE:      " << JPSITYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, JPSITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: J/#psi ID :: %d", CutName, JPSITYPE));
    }

    if (!strcmp(CutName, "JPSIMASSLO")) {
      JPSIMASSLO = CutValue; ok = 1;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSLO << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, JPSIMASSLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSLO));
    }

    if (!strcmp(CutName, "JPSIMASSHI")) {
      JPSIMASSHI = CutValue; ok = 1;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSHI << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, JPSIMASSHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSHI));
    }

    if (!strcmp(CutName, "MKKLO")) {
      MKKLO = CutValue; ok = 1;
      if (dump) cout << "MKKLO:           " << MKKLO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(KK) :: %3.1f", CutName, MKKLO));
    }

    if (!strcmp(CutName, "MKKHI")) {
      MKKHI = CutValue; ok = 1;
      if (dump) cout << "MKKHI:           " << MKKHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(KK) :: %3.1f", CutName, MKKHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue; ok = 1;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }


    
  }

}
