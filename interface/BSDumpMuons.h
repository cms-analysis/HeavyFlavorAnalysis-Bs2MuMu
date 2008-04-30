#ifndef _BSDUMPMUONS_h_
#define _BSDUMPMUONS_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <string>
#include <map>
#include <set>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>
#include <TH2.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class TFile;
class TTree;
class TH2D;
class TAna00Event;

// ----------------------------------------------------------------------
class BSDumpMuons : public edm::EDAnalyzer {
 public:
  explicit BSDumpMuons(const edm::ParameterSet&);
  ~BSDumpMuons();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose;
  int           fNevt;

  edm::InputTag fMuonsLabel1;
  edm::InputTag fMuonsLabel2;

  edm::InputTag fTracksLabel1;
  edm::InputTag fTracksLabel2;
  edm::InputTag fTracksLabel3;
  edm::InputTag fTracksLabel4;
  edm::InputTag fTracksLabel5;
  edm::InputTag fTracksLabel6;
  edm::InputTag fTracksLabel7;

};

#endif
