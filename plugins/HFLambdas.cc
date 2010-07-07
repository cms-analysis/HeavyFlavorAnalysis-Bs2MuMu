#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFLambdas.h"

#include <iostream>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFThreeParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Yikes!
extern TAna01Event *gHFEvent;

//using namespace std;   // strange enough these namespaces are declared elsewhere....
//using namespace reco;
//using namespace edm;

// ----------------------------------------------------------------------
HFLambdas::HFLambdas(const ParameterSet& iConfig) :
    fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
    fMaxTracks(iConfig.getUntrackedParameter<int>("maxTracks", 1000)),
    fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))),
    fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
    fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
    fUseMuon(iConfig.getUntrackedParameter<int>("useMuon", 0)),
    fJPsiWindow(iConfig.getUntrackedParameter<double>("JPsiWindow", 0.4)),
    fL0Window(iConfig.getUntrackedParameter<double>("L0Window", 0.4)),
    fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
    fProtonPt(iConfig.getUntrackedParameter<double>("protonPt", 1.0)),
    fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 1.0)),
    fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)),
    fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)),
    fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 9999.0)),
    fMaxVtxChi2(iConfig.getUntrackedParameter<double>("maxVtxChi2", 9999.0)),
    fMinVtxSigXY(iConfig.getUntrackedParameter<double>("minVtxSigXY", -1.)),
    fMinVtxSig3d(iConfig.getUntrackedParameter<double>("minVtxSig3d", -1.)),
    fMinCosAngle(iConfig.getUntrackedParameter<double>("minCosAngle", -1.)),
    fMinPtCand(iConfig.getUntrackedParameter<double>("minPtCand", -99.)),
    fType(iConfig.getUntrackedParameter<int>("type", 1)) {
    std::cout << "----------------------------------------------------------------------" << endl;
    std::cout << "--- HFLambdas constructor" << std::endl;
    std::cout << "---  verbose:                  " << fVerbose << std::endl;
    std::cout << "---  tracksLabel:              " << fTracksLabel << std::endl;
    std::cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << std::endl;
    std::cout << "---  muonsLabel:               " << fMuonsLabel << std::endl;
    std::cout << "---  useMuon:                  " << fUseMuon << std::endl;
    std::cout << "---  phiWindow:                " << fPhiWindow << std::endl;
    std::cout << "---  JPsiWindow:               " << fJPsiWindow << std::endl;
    std::cout << "---  L0Window:                 " << fL0Window << std::endl;
    std::cout << "---  muonPt:                   " << fMuonPt << std::endl;
    std::cout << "---  protonPt:                 " << fProtonPt << std::endl;
    std::cout << "---  pionPt:                   " << fPionPt << std::endl;
    std::cout << "---  trackPt:                  " << fTrackPt << std::endl;
    std::cout << "---  deltaR:                   " << fDeltaR << std::endl;

    std::cout << "---  maxDoca:                  " << fMaxDoca << std::endl;
    std::cout << "---  maxVtxChi2:               " << fMaxVtxChi2 << std::endl;
    std::cout << "---  minVtxSigXY:              " << fMinVtxSigXY << std::endl;
    std::cout << "---  minVtxSig3d:              " << fMinVtxSig3d << std::endl;
    std::cout << "---  minCosAngle:              " << fMinCosAngle << std::endl;
    std::cout << "---  minPtCand:                " << fMinPtCand << std::endl;

    std::cout << "---  type:                     " << fType << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;

}


// ----------------------------------------------------------------------
HFLambdas::~HFLambdas() {

}


// ----------------------------------------------------------------------
void HFLambdas::analyze(const Event& iEvent, const EventSetup& iSetup) {

    if (fVerbose > 0) {
        std::cout << "-------------------------------------------------------------" << std::endl;
        std::cout << "==>HFLambdas: beginning of analyze():" << std::endl;
        std::string line("ps -F " + getpid());
        system(line.c_str());
    }

    // get the primary vertex
    Handle<VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
    if(!recoPrimaryVertexCollection.isValid()) {
        std::cout << "==>HFLambdas> No primary vertex collection found, skipping" << std::endl;
        return;
    }
    const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
    if (vertices.size() == 0) {
        std::cout << "==>HFLambdas> No primary vertex found, skipping" << std::endl;
        return;
    }
    fPV = vertices[gHFEvent->fBestPV];
    if (fVerbose > 0) {
        std::cout << "HFDimuons: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << std::endl;
    }

    // get the collection of muons
    Handle<MuonCollection> hMuons;
    iEvent.getByLabel(fMuonsLabel, hMuons);
    if (!hMuons.isValid()) {
        std::cout << "==>HFLambdas> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << std::endl;
        return;
    }

    // get the collection of tracks
    Handle<View<Track> > hTracks;
    iEvent.getByLabel(fTracksLabel, hTracks);
    if(!hTracks.isValid()) {
        std::cout << "==>HFLambdas> No valid TrackCollection with label " << fTracksLabel << " found, skipping" << std::endl;
        return;
    }

    if (hTracks->size() > static_cast<unsigned int>(fMaxTracks)) {
        std::cout << "==>HFLambdas> Too many tracks " << hTracks->size() << ", skipping" << std::endl;
        return;
    }
    if (fVerbose > 0) {
        std::cout << "==>HFLambdas> ntracks = " << hTracks->size() << std::endl;
    }

    // Transient tracks for vertexing
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
    if (!fTTB.isValid()) {
        std::cout << " -->HFLambdas: Error: no TransientTrackBuilder found."<<std::endl;
        return;
    }

    // get the collection of muons and store their corresponding track indices
    std::vector<unsigned int> muonIndices;
    for (MuonCollection::const_iterator muonIt = hMuons->begin(); muonIt != hMuons->end(); ++muonIt) {
        const int im = muonIt->track().index();
	if (fVerbose > 0) std::cout << "==>HFLambdas> Muon index: " << im << " Ptr: " << muonIt->track().get() << std::endl;
        if (im >= 0) muonIndices.push_back(im);
    }
    if (fVerbose > 0) {
        std::cout << "==>HFLambdas> nMuons = " << hMuons->size() << std::endl;
        std::cout << "==>HFLambdas> nMuonIndices = " << muonIndices.size() << std::endl;
    }

    // Build lists
    trackList_t prList, piList, trackMuonList;
    trackMuonList.reserve(200);
    piList.reserve(2000);
    prList.reserve(2000);

    for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack) {
        TrackBaseRef rTrackView(hTracks, itrack);
        Track tTrack(*rTrackView);
	//if (fVerbose > 0) std::cout << "==>HFLambdas> Track Ptr: " << rTrackView.get() << std::endl;

	// fill pion list
        if (tTrack.pt() > fPionPt)  {
	    TLorentzVector tlv;
            tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPION);
            piList.push_back(std::make_pair(itrack, tlv));
        }

	// fill proton list
        if (tTrack.pt() > fProtonPt) {
	    TLorentzVector tlv;
            tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPROTON);
            prList.push_back(std::make_pair(itrack, tlv));
        }

	// fill trackMuonList
        if (tTrack.pt() > fMuonPt) {
	    //if (fVerbose > 0) std::cout << "==>HFLambdas> added " << rTrackView.get() << " to trackMuonList" << std::std::endl;
	    TLorentzVector tlv;
            tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MMUON);
            trackMuonList.push_back(std::make_pair(itrack, tlv));
        }

    }

    // Now create the combinatorics for the J/psi
    // For each muon from the muon chambers combine with all tracks in trackMuonList
    // trackMuonList contain some of the muonsystem-muons as well
    std::vector<duplet_t> jpsiList;
    for(std::vector<unsigned int>::const_iterator itm=muonIndices.begin(); itm!=muonIndices.end(); itm++) {
	for(trackList_t::const_iterator ittrm=trackMuonList.begin(); ittrm!=trackMuonList.end(); ittrm++) {
	    if( (*itm) != ittrm->first ) { // then we have two distinct tracks
		jpsiList.push_back(std::make_pair( (*itm), ittrm->first));
	    }
	}
    }

    if (fVerbose > 0) std::cout << "==>HFLambdas> jpsiList size: " << jpsiList.size() << std::endl;

    HFKalmanVertexFit  hkvfitter(fTTB.product(), fPV, 0, fVerbose+10);
    std::vector<Track> trackList;
    std::vector<int> trackIndices;
    std::vector<double> trackMasses;

    hkvfitter.setNoCuts();
    hkvfitter.fMaxDoca     = fMaxDoca;
    hkvfitter.fVtxChi2     = fMaxVtxChi2;
    hkvfitter.fVtxSigXY    = fMinVtxSigXY;
    hkvfitter.fVtxSig3d    = fMinVtxSig3d;
    hkvfitter.fCosAngle    = fMinCosAngle;
    hkvfitter.fPtCand      = fMinPtCand;

    for (std::vector<duplet_t>::iterator it=jpsiList.begin(); it!=jpsiList.end(); ++it) {

        TrackBaseRef mu1TrackView(hTracks, it->first);
        Track tMu1(*mu1TrackView);
        TLorentzVector tlvMu1;
	tlvMu1.SetPtEtaPhiM(tMu1.pt(), tMu1.eta(), tMu1.phi(), MMUON);

        TrackBaseRef mu2TrackView(hTracks, it->second);
        Track tMu2(*mu2TrackView);
	TLorentzVector tlvMu2;
        tlvMu2.SetPtEtaPhiM(tMu2.pt(), tMu2.eta(), tMu2.phi(), MMUON);

        if (tMu2.charge() == tMu1.charge()) continue;          // muons must have opposite charge to be from a J/Psi

        trackList.clear();
        trackIndices.clear();
        trackMasses.clear();

        trackList.push_back(tMu1);
        trackIndices.push_back(it->first);
        trackMasses.push_back(MMUON);

        trackList.push_back(tMu2);
        trackIndices.push_back(it->second);
        trackMasses.push_back(MMUON);

        TLorentzVector tlvJPsi = tlvMu1 + tlvMu2;
        if ((TMath::Abs(tlvJPsi.M() - MJPSI) < fJPsiWindow )) {
	    if (fVerbose > 0) std::cout << "==>HFLambdas> added to tlvJPsi with mass " << tlvJPsi.M() << std::endl;
            hkvfitter.doFit(trackList, trackIndices, trackMasses, fType*10000+443, 2);
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFLambdas::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFLambdas::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFLambdas);
