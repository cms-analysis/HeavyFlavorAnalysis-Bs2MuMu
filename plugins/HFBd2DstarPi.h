#ifndef _HFBD2DSTARPI_h_
#define _HFBD2DSTARPI_h_

#include "HeavyFlavorAnalysis/Bs2MuMu/plugins/HFVirtualDecay.h"


// ----------------------------------------------------------------------

class HFBd2DstarPi : public HFVirtualDecay {
	public:
		explicit HFBd2DstarPi(const edm::ParameterSet&);
	
	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		
		virtual void dumpConfiguration();
	private:
		// additional variables needed for this decay
		double fSlowPionPt;
		double fD0Window,fDeltaM;
};

#endif
