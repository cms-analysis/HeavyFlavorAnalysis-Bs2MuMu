/** \class L1TrigReport
 *
 * See header file for documentation
 *
 *  $Date: 2007/12/18 08:28:02 $
 *  $Revision: 1.11 $
 *
 *  \author Martin Grunewald
 *
 */

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/L1TrigReport.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TFile.h>
#include <TH1.h>

#include <iomanip>
#include <iostream>


// -- Yikes!
extern TAna00Event *gHFEvent;
extern TFile       *gHFFile;

//
// constructors and destructor
//
L1TrigReport::L1TrigReport(const edm::ParameterSet& iConfig) :
  l1GTReadoutRecTag_(iConfig.getParameter<edm::InputTag> ("L1GTReadoutRecord")),
  nEvents_(0),
  nErrors_(0),
  nAccepts_(0),
  l1Accepts_(0),
  l1Names_(0),
  init_(false),
  nSize_(0)
{
  LogDebug("") << "Level-1 Global Trigger Readout Record: " + l1GTReadoutRecTag_.encode();
  fL1 = new TH1D("l1", "L1 names", 128, 0., 128. );
}

L1TrigReport::~L1TrigReport()
{ 
  // -- Save output
  gHFFile->cd();

  fL1->Write();

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void L1TrigReport::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // accumulation of statistics event by event

  using namespace std;
  using namespace edm;

  // get hold of L1GlobalReadoutRecord
  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(l1GTReadoutRecTag_,L1GTRR); 

  if (L1GTRR.isValid()) {
    const unsigned int n(L1GTRR->decisionWord().size());
    // initialisation
    if ( (!init_) || (nSize_!=n) ) {
      if (!init_) {
	init_=true;
      } else {
	endJob();
	nEvents_=0;
	nErrors_=0;
	nAccepts_=0;
      }
      nSize_=n;
      l1Names_.resize(n);
      l1Accepts_.resize(n);
      for (unsigned int i=0; i!=n; ++i) {
	l1Accepts_[i]=0;
	l1Names_[i]="NameNotAvailable";
      }
    }
    const bool accept(L1GTRR->decision());
    LogDebug("") << "L1GlobalTriggerReadoutRecord decision: " << accept;
    nEvents_++;
    if (accept) ++nAccepts_;
    // decision for each L1 algorithm
    for (unsigned int i=0; i!=n; ++i) {
      if (L1GTRR->decisionWord()[i]) l1Accepts_[i]++;
    }
  } else {
    LogDebug("") << "L1GlobalTriggerReadoutRecord with label ["+l1GTReadoutRecTag_.encode()+"] not found!";
    nEvents_++;
    nErrors_++;
  }

  return;

}
void  L1TrigReport::beginJob() {

  gHFFile->cd();

}

void L1TrigReport::endJob()
{
  // final printout of accumulated statistics

  using namespace std;

    cout << dec << endl;
    cout << "L1T-Report " << "---------- Event  Summary ------------\n";
    cout << "L1T-Report"
	 << " Events total = " << nEvents_
	 << " passed = " << nAccepts_
	 << " failed = " << nEvents_-nErrors_-nAccepts_
	 << " errors = " << nErrors_
	 << "\n";

    cout << endl;
    cout << "L1T-Report " << "---------- L1Trig Summary ------------\n";
    cout << "L1T-Report "
	 << right << setw(10) << "L1T  Bit#" << " "
	 << right << setw(10) << "Passed" << " "
	 << right << setw(10) << "Failed" << " "
	 << right << setw(10) << "Errors" << " "
	 << "Name" << "\n";

  if (init_) {
    for (unsigned int i=0; i!=nSize_; ++i) {
      cout << "L1T-Report "
	   << right << setw(10) << i << " "
	   << right << setw(10) << l1Accepts_[i] << " "
	   << right << setw(10) << nEvents_-nErrors_-l1Accepts_[i] << " "
	   << right << setw(10) << nErrors_ << " "
	   << l1Names_[i] << "\n";

      fL1->SetBinContent(i, l1Accepts_[i]);
      fL1->GetXaxis()->SetBinLabel(i, Form("%s", l1Names_[i].c_str()));
    }
  } else {
    cout << "L1T-Report - No L1 GTRRs found!" << endl;
  }

    cout << endl;
    cout << "L1T-Report end!" << endl;
    cout << endl;

    return;
}
