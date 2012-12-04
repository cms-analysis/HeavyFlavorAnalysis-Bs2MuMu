#ifndef HFPYTHIADAUVFILTER_h
#define HFPYTHIADAUVFILTER_h
// -*- C++ -*-
//
// Package:    HFPythiaDauVFilter
// Class:      HFPythiaDauVFilter
// 
/**\class HFPythiaDauVFilter HFPythiaDauVFilter.cc 

 Description: Filter events using MotherId and ChildrenIds infos

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Daniele Pedrini
//         Created:  Apr 29 2008
// $Id: PythiaDauVFilter.h,v 1.4 2011/09/12 12:39:54 ursl Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class HFPythiaDauVFilter : public edm::EDFilter {
 public:
  explicit HFPythiaDauVFilter(const edm::ParameterSet&);
  ~HFPythiaDauVFilter();
  
  
  virtual bool filter(edm::Event&, const edm::EventSetup&);
 private:
  int fVerbose;  
  std::string label_;
  std::vector<int> dauIDs;
  int particleID;
  int motherID;
  bool chargeconju; 
  int ndaughters;
  std::vector<double> minptcut;
  double maxptcut;
  std::vector<double> minetacut;
  std::vector<double> maxetacut;
};
#define PYCOMP pycomp_
extern "C" {
  int PYCOMP(int& ip);
} 
#endif
