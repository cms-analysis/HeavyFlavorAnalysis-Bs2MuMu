#include "candAnaBu2JpsiK.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBu2JpsiK::candAnaBu2JpsiK(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaBu2JpsiK: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaBu2JpsiK::~candAnaBu2JpsiK() {
  cout << "==> candAnaBu2JpsiK: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::candAnalysis() {

  if (0 == fpCand) return;

  TAnaTrack *p0, *pk(0), *pks(0); 
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (321 == TMath::Abs(p0->fMCID)) {
      pks = p0; 
    }
  }

  if (0 == pks) {
    cout << "candAnaBu2JpsiK::candAnalysis:  no kaon found " << endl;
    return;
  }

  pk = fpEvt->getRecTrack(pks->fIndex);

  fKaonPt        = pk->fPlab.Perp(); 
  fKaonEta       = pk->fPlab.Eta();  
  fKaonPhi       = pk->fPlab.Phi(); 
  fKaonTkQuality = pk->fTrackQuality & TRACKQUALITY;
  fKaonPtNrf     = pks->fPlab.Perp();
  fKaonEtaNrf    = pks->fPlab.Eta();

//FIXME  
//   if (fIsMC) {
//     fGenBpartial = partialReco(fpCand); 
//     if (fpCand->fMass > 5.5 && fGenBpartial == 511) fpEvt->dumpGenBlock();
//   }

//FIXME
//   if (tmCand(fpCand)) {
//     TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(pk->fIndex)->fGenIndex);
//     fKPtGen     = pg1->fP.Perp();
//     fKEtaGen    = pg1->fP.Eta();
//   } else {
//     fKPtGen     = -99.;
//     fKEtaGen    = -99.;
//   }
  

  // -- Check for J/psi mass
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp(); 
      fJpsiEta  = pD->fPlab.Eta(); 
      fJpsiPhi  = pD->fPlab.Phi(); 
      break;
      //       cout << "type = " << pD->fType 
      // 	   << " with mass = " << pD->fMass 
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
      cout << "bmmNormalizationReader::fillCandidateVariables:  no muon 1 found " << endl;
      return; 
    }
    if (0 == p2) {
      cout << "bmmNormalizationReader::fillCandidateVariables:  no muon 2 found " << endl;
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

  
  candAna::candAnalysis();
  fPreselection = fPreselection && fGoodJpsiMass;

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
void candAnaBu2JpsiK::moreBasicCuts() {
  cout << "   candAnaBs2JpsiPhi: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::genMatch() {

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::candMatch() {

}




// ----------------------------------------------------------------------
void candAnaBu2JpsiK::bookHist() {
  candAna::bookHist();
  cout << "==>candAnaBu2JpsiK: bookHist" << endl;

  // -- Additional reduced tree variables
  fTree->Branch("mpsi", &fJpsiMass,  "mpsi/D");
  fTree->Branch("kpt",  &fKaonPt,    "kpt/D");
  fTree->Branch("keta", &fKaonEta,   "keta/D");
  fTree->Branch("kphi", &fKaonPhi,   "kphi/D");
  fTree->Branch("kgt",  &fKaonTkQuality,"kgt/I");
  fTree->Branch("t3pt", &fKaonPtNrf, "t3pt/D");
  fTree->Branch("t3eta",&fKaonEtaNrf,"t3eta/D");
  fTree->Branch("g3pt", &fKPtGen,    "g3pt/D");
  fTree->Branch("g3eta",&fKEtaGen,   "g3eta/D");

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiK::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaBu2JpsiK: Reading " << fCutFile << " for cut settings" << endl;
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
    
  }

}
