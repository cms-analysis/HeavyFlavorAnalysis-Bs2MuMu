#include "candAnaBs2JpsiPhi.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::candAnaBs2JpsiPhi(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaBs2JpsiPhi: name = " << name << ", reading cutsfile " << cutsFile << endl;
}


// ----------------------------------------------------------------------
candAnaBs2JpsiPhi::~candAnaBs2JpsiPhi() {
  cout << "==> candAnaBs2JpsiPhi: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBs2JpsiPhi::candAnalysis() {
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
  fTree->Branch("mmkspecial",      &fMMKSpecial,  "mmkspecial/D");

}
