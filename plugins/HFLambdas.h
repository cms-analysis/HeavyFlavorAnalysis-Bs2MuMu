#ifndef _HFBD2JPSIKS_H_
#define _HFBD2JPSIKS_H_

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class HFLambdas : public edm::EDAnalyzer
{
public:
    explicit HFLambdas(const edm::ParameterSet&);
    ~HFLambdas();

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    int fVerbose;
    edm::InputTag	fTracksLabel,fPrimaryVertexLabel;
    edm::InputTag	fMuonsLabel, fBeamSpotLabel;
    edm::InputTag       fMuonType;

    double fMuonPt,fPionPt,fProtonPt;
    muon::SelectionType fMuonSelType;
    double fTrackNormChi2;
    int fPsiMuons;
    double fPsiWindow,fksWindow,fL0Window,fLbWindow,fB0Window; // mass windows for event selection
    bool fUseAnalysisValuesForEff; // use same values for mass windows and node cuts as in analysis candidates
    double fPsiEffWindow,fL0EffWindow; // windows for efficiencies of J/Psi and Lambda0
    double fEffMaxChi2, fEffMin3d, fEffMaxDoca; // cuts for efficiencies of J/Psi and Lambda0
    double fDeltaR;
    double fMuonEtaMax;
    bool fDoB0, fDoLb;
    //int fVertexing;

    double fMaxDoca; // for LambdaCut
    double fMaxVtxChi2;
    double fPAngle;

    bool fUseV0producer;
    bool fDoVcands;
    bool fRemoveCandTracksFromVertex;

    reco::Vertex	fPV;
    edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
