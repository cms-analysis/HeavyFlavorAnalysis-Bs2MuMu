#include "candAnaBd2DstarPi.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBd2DstarPi::candAnaBd2DstarPi(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  BLIND = 0; 
  cout << "==> candAnaBd2DstarPi: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaBd2DstarPi::~candAnaBd2DstarPi() {
  cout << "==> candAnaBd2DstarPi: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::candAnalysis() {

  if (0 == fpCand) return;


  // -- Check for D0  mass
  TAnaCand *pD = 0; 
  fGoodMD0 = fGoodDeltaM = fGoodMDs = false; 
  TLorentzVector ppi, ppipi, ppislow, pd0, pdstar; 
  if (fpCand->fDau1 < 0 || fpCand->fDau1 > fpEvt->nCands()) return;
  pD = fpEvt->getCand(fpCand->fDau1); 
  if (0 == pD) return;
  //   cout << "----------------------------------------------------------------------" << endl;
  //   cout << " fpCand = " << fpCand << " type = " << fpCand->fType << " mass = " << fpCand->fMass << endl;
  //   cout << " 1 pD   = " << pD << " type = " << pD->fType << " mass = " << pD->fMass << endl;
  // -- get DSTAR (sub)cand

  if (pD->fType == DSTARTYPE) {
    fMDs =  pD->fMass;  
  }
  
  // -- get D0 (sub)cand
  pD = fpEvt->getCand(pD->fDau1); 
  //   cout << " 2 pD   = " << pD << " type = " << pD->fType << " mass = " << pD->fMass << endl;
  if (pD->fType == D0TYPE) {
    fMD0 = pD->fMass; 
    fPTD0 = pD->fPlab.Perp();
    if ((MD0LO < fMD0) && (fMD0 < MD0HI)) {
      fGoodMD0 = true;
    }
  }

  fDeltaM = fMDs - fMD0; 
  //   cout << "  deltaM = " << fDeltaM << endl;
  if (TMath::Abs(fDeltaM - 0.145) < DELTAM) {
    fGoodDeltaM = true;
  }

  fCandM = fpCand->fMass;

  ((TH1D*)fHistDir->Get("hdm"))->Fill(fDeltaM); 
  ((TH1D*)fHistDir->Get("hmdz"))->Fill(fMD0); 
  ((TH1D*)fHistDir->Get("hmds"))->Fill(fMDs); 
  ((TH1D*)fHistDir->Get("hmb0"))->Fill(fCandM); 

  if (fGoodMD0)    ((TH1D*)fHistDir->Get("chdm"))->Fill(fDeltaM); 
  if (fGoodDeltaM) ((TH1D*)fHistDir->Get("chmdz"))->Fill(fMD0); 
  if (fGoodDeltaM && fGoodMD0) ((TH1D*)fHistDir->Get("chmds"))->Fill(fMDs); 
  if (fGoodDeltaM && fGoodMD0) ((TH1D*)fHistDir->Get("chmb0"))->Fill(fCandM); 

  //  fPTPiS, fPTPi2, fPTPi, fPTKa;

  candAna::candAnalysis();
  fCandTau   = fCandFL3d*MBPLUS/fCandP/TMath::Ccgs();

  fPreselection = fWideMass && fGoodTracks && fGoodTracksPt && fGoodTracksEta && fGoodMuonsPt && fGoodMuonsEta; 
  fPreselection = fPreselection && fGoodQ;
  fPreselection = fPreselection && (fCandA < 0.2) && (fCandFLS3d > 2) && (fCandChi2/fCandDof < 5); 
  fPreselection = fPreselection && fGoodMD0 && fGoodDeltaM;
  fPreselection = true;
  return;

}

// ----------------------------------------------------------------------
void candAnaBd2DstarPi::moreBasicCuts() {
  cout << "   candAnaBd2DstarPi: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodMD0", "m(D0)", fGoodMD0); 
  fAnaCuts.addCut("fGoodDeltaM", "Delta M", fGoodDeltaM); 
  fAnaCuts.addCut("fGoodDeltaR", "Delta R", fGoodDeltaR); 
}

// ----------------------------------------------------------------------
void candAnaBd2DstarPi::genMatch() {
  

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::candMatch() {

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::bookHist() {
  cout << "==>candAnaBd2DstarPi: bookHist" << endl;
  candAna::bookHist();

  TH1D *hdm = new TH1D("hdm", "delta(m)", 50, 0.135, 0.16); 
  hdm = new TH1D("hmdz", "mdz", 50, 1.75, 1.95); 
  hdm = new TH1D("hmds", "mds", 50, 1.9, 2.1); 
  hdm = new TH1D("hmb0", "mb0", 50, 4.9, 5.9); 

  hdm = new TH1D("chdm", "delta(m)", 50, 0.135, 0.16); 
  hdm = new TH1D("chmdz", "mdz", 50, 1.75, 1.95); 
  hdm = new TH1D("chmds", "mds", 50, 1.9, 2.1); 
  hdm = new TH1D("chmb0", "mb0", 50, 4.9, 5.9); 
  (void)hdm; // remove compiler warning

  moreReducedTree(fTree);
  moreReducedTree(fAmsTree);

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::moreReducedTree(TTree *t) {
  // -- Additional reduced tree variables
  t->Branch("ptd0", &fPTD0, "ptd0/D");
  t->Branch("ptpis", &fPTPiS, "ptpis/D");
  t->Branch("ptpi2", &fPTPi2, "ptpi2/D");
  t->Branch("ptpi",  &fPTPi,  "ptpi/D");
  t->Branch("ptka",  &fPTKa,  "ptka/D");
  t->Branch("md0", &fMD0, "md0/D");
  t->Branch("dm",  &fDeltaM, "dm/D");
}

// ----------------------------------------------------------------------
void candAnaBd2DstarPi::fillCandidateHistograms(int offset) {

  candAna::fillCandidateHistograms(offset); 

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::efficiencyCalculation() {

}


// ----------------------------------------------------------------------
void candAnaBd2DstarPi::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaBd2DstarPi: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "D0TYPE")) {
      D0TYPE = static_cast<int>(CutValue); 
      if (dump) cout << "D0TYPE:      " << D0TYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, D0TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: D0 ID :: %d", CutName, D0TYPE));
    }

    if (!strcmp(CutName, "DSTARTYPE")) {
      DSTARTYPE = static_cast<int>(CutValue); 
      if (dump) cout << "DSTARTYPE:      " << DSTARTYPE << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, DSTARTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: DSTAR ID :: %d", CutName, DSTARTYPE));
    }

    if (!strcmp(CutName, "MD0LO")) {
      MD0LO = CutValue; 
      if (dump) cout << "MD0LO:      " << MD0LO << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, MD0LO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MD0LO :: %3.1f", CutName, MD0LO));
    }

    if (!strcmp(CutName, "MD0HI")) {
      MD0HI = CutValue; 
      if (dump) cout << "MD0HI:      " << MD0HI << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, MD0HI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MD0HI :: %3.1f", CutName, MD0HI));
    }

    if (!strcmp(CutName, "DELTAM")) {
      DELTAM = CutValue; 
      if (dump) cout << "DELTAM:      " << DELTAM << endl;
      ibin = 214;
      hcuts->SetBinContent(ibin, DELTAM);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: deltaM :: %f", CutName, DELTAM));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue; 
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }


    
  }

}
