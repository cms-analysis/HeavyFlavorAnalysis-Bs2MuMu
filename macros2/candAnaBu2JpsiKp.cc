#include "candAnaBu2JpsiKp.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBu2JpsiKp::candAnaBu2JpsiKp(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaBu2JpsiKp: name = " << name << ", reading cutsfile " << cutsFile << endl;
}


// ----------------------------------------------------------------------
candAnaBu2JpsiKp::~candAnaBu2JpsiKp() {
  cout << "==> candAnaBu2JpsiKp: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBu2JpsiKp::candAnalysis() {
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
void candAnaBu2JpsiKp::processType() {

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiKp::genMatch() {

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiKp::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaBu2JpsiKp::candMatch() {

}




// ----------------------------------------------------------------------
void candAnaBu2JpsiKp::bookHist() {
  candAna::bookHist();
  cout << "==>candAnaBu2JpsiKp: bookHist" << endl;
  fTree->Branch("mmkspecial",      &fMMKSpecial,  "mmkspecial/D");

}
