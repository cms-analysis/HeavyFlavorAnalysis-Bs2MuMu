#ifndef _HFDSTAR_h_
#define _HFDSTAR_h_

#include "HeavyFlavorAnalysis/Bs2MuMu/plugins/HFVirtualDecay.h"

class HFDstar : public HFVirtualDecay {
	public:
		explicit HFDstar(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();

		double fSlowPionPt; 
		double fD0Window;
		double fDeltaM;
};

#endif
