#include "candAnaMuMu.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaMuMu::candAnaMuMu(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candMuMuAna: constructor..." << endl;
}


// ----------------------------------------------------------------------
candAnaMuMu::~candAnaMuMu() {
  cout << "==> candMuMuAna: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaMuMu::candAnalysis() {
  candAna::candAnalysis();
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(2); 

  if (fIsMC) {
    fTree->Fill(); 
  } else {
    if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX) {
      // do nothing
    } else {
      if (fPreselection) {
	((TH1D*)fHistDir->Get("../monEvents"))->Fill(12); 
	fTree->Fill(); 
      }         
    }
  }

}

// ----------------------------------------------------------------------
void candAnaMuMu::processType() {

}


// ----------------------------------------------------------------------
void candAnaMuMu::genMatch() {

}


// ----------------------------------------------------------------------
void candAnaMuMu::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaMuMu::candMatch() {

}



// ----------------------------------------------------------------------
void candAnaMuMu::bookHist() {
  candAna::bookHist();
  cout << "==>candAnaMuMu: bookHist" << endl;
  cout << "pwd(): "; fHistDir->pwd();
  cout << "fTree: " << fTree << endl;
  fTree->Branch("mmspecial",      &fMMSpecial,  "mmspecial/D");

}
