#ifndef _HFBD2JPSIKS_H_
#define _HFBD2JPSIKS_H_

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class HFLambdas : public edm::EDAnalyzer {
	public:
		explicit HFLambdas(const edm::ParameterSet&);
		~HFLambdas();
	
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		
		int fVerbose;
		edm::InputTag	fTracksLabel,fPrimaryVertexLabel;
		edm::InputTag	fMuonsLabel;

		double fMuonPt,fPionPt;
		int fPsiMuons;
		double fPsiWindow,fL0Window,fLbWindow;
		double fDeltaR;
		int fVertexing;
		
		double fMaxDoca;
		double fPAngleL0;
		
		reco::Vertex	fPV;
		edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
