#ifndef _HFDITRACKS_h_
#define _HFDITRACKS_h_

#include <string>

#include "HeavyFlavorAnalysis/Bs2MuMu/plugins/HFVirtualDecay.h"

class HFDiTracks : public HFVirtualDecay {
	public:
		explicit HFDiTracks(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();
		
		virtual int  idFromMass(double mass);

		double fTrack1Mass;
		double fTrack2Mass;
		double fMassLow;
		double fMassHigh;
		
		int fNbrMuons;
		bool fCloseToMuons;
};

#endif
