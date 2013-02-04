#ifndef _HFBS2JPSIPHI_h_
#define _HFBS2JPSIPHI_h_

#include "HeavyFlavorAnalysis/Bs2MuMu/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFBs2JpsiPhi : public HFVirtualDecay {
	public:
		explicit HFBs2JpsiPhi(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();
		
		int           fPsiMuons;
		double        fPsiWindow, fPhiWindow, fBsWindow;
};

#endif
