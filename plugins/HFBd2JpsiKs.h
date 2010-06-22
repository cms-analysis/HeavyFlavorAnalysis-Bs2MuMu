#ifndef _HFBD2JPSIKS_H_
#define _HFBD2JPSIKS_H_

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class HFBd2JpsiKs : public edm::EDAnalyzer {
	public:
		explicit HFBd2JpsiKs(const edm::ParameterSet&);
		~HFBd2JpsiKs();
	
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		
		int fVerbose;
		edm::InputTag	fTracksLabel,fPrimaryVertexLabel;
		edm::InputTag	fMuonsLabel;

		double fMuonPt,fPionPt;
		int fPsiMuons;
		double fPsiWindow,fKsWindow,fBdWindow;
		double fDeltaR;
		int fVertexing;
		
		double fMaxDoca;
		double fPAngleKs;
		
		reco::Vertex	fPV;
		edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
