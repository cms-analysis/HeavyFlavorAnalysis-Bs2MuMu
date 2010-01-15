#ifndef _HFGENFILTER_h_
#define _HFGENFILTER_h_

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


class HFGenFilter : public edm::EDFilter {

public:

  explicit HFGenFilter(const edm::ParameterSet&);
  ~HFGenFilter();
  
  virtual bool filter(edm::Event&, const edm::EventSetup&);

private:
  
  int            fNgood, fNtot, fEvt;

};

#endif
