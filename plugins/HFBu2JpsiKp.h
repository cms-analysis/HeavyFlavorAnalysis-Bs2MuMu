#ifndef _HFBU2JPSIKP_h_
#define _HFBU2JPSIKP_h_

#include "HeavyFlavorAnalysis/Bs2MuMu/plugins/HFVirtualDecay.h"

class HFBu2JpsiKp : public HFVirtualDecay {

	public:
		explicit HFBu2JpsiKp(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();

		int           fPsiMuons;
		double        fPsiWindow, fBuWindow;
};

#endif
