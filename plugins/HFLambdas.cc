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

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

// Yikes!
extern TAna01Event *gHFEvent;

//using namespace std;   // strange enough these namespaces are declared elsewhere....
//using namespace reco;
//using namespace edm;

// ----------------------------------------------------------------------
HFLambdas::HFLambdas(const edm::ParameterSet& iConfig) :
    fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
    fMaxTracks(iConfig.getUntrackedParameter<int>("maxTracks", 300)),
    fMinTracks(iConfig.getUntrackedParameter<int>("minTracks", 4)),
    fTracksLabel(iConfig.getUntrackedParameter<edm::InputTag>("tracksLabel", edm::InputTag("goodTracks"))),
    fPrimaryVertexLabel(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertexLabel", edm::InputTag("offlinePrimaryVertices"))),
    fMuonsLabel(iConfig.getUntrackedParameter<edm::InputTag>("muonsLabel")),
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
    fMinPocaJpsi(iConfig.getUntrackedParameter<double>("minPocaJpsi", 0.)),
    fMinPocaL0(iConfig.getUntrackedParameter<double>("minPocaL0", 0.)),
    fType(iConfig.getUntrackedParameter<int>("type", 1)) {
    std::cout << "----------------------------------------------------------------------" << endl;
    std::cout << "--- HFLambdas constructor" << std::endl;
    std::cout << "---  verbose:                  " << fVerbose << std::endl;
    std::cout << "---  maxTracks:                " << fMaxTracks << std::endl;
    std::cout << "---  minTracks:                " << fMinTracks << std::endl;
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
void HFLambdas::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    if (fVerbose > 0) {
        std::cout << "-------------------------------------------------------------" << std::endl;
        std::cout << "==>HFLambdas: beginning of analyze():" << std::endl;
        std::string line("ps -F " + getpid());
        //system(line.c_str());
    }

    // get the primary vertex
    edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
    if(!recoPrimaryVertexCollection.isValid()) {
        std::cout << "==>HFLambdas> No primary vertex collection found, skipping" << std::endl;
        return;
    }
    const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
    if (vertices.size() == 0) {
        std::cout << "==>HFLambdas> No primary vertex found, skipping" << std::endl;
        return;
    }
    fPV = vertices[gHFEvent->fBestPV];
    if (fVerbose > 0) {
        std::cout << "HFDimuons: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << std::endl;
    }

    // get the collection of muons
    edm::Handle<reco::MuonCollection> hMuons;
    iEvent.getByLabel(fMuonsLabel, hMuons);
    if (!hMuons.isValid()) {
        std::cout << "==>HFLambdas> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << std::endl;
        return;
    }

    // get the collection of tracks
    edm::Handle<edm::View<reco::Track> > hTracks;
    iEvent.getByLabel(fTracksLabel, hTracks);
    if(!hTracks.isValid()) {
        std::cout << "==>HFLambdas> No valid TrackCollection with label " << fTracksLabel << " found, skipping" << std::endl;
        return;
    }

    if (hTracks->size() > static_cast<size_t>(fMaxTracks)) {
        std::cout << "==>HFLambdas> Too many tracks " << hTracks->size() << ", skipping" << std::endl;
        return;
    }
    if (hTracks->size() < static_cast<size_t>(fMinTracks)) {
	std::cout << "==>HFLambdas> Not enough tracks " << hTracks->size() << ", skipping" << std::endl;
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
    std::vector<index_t> muonIndices;
    for (reco::MuonCollection::const_iterator muonIt = hMuons->begin(); muonIt != hMuons->end(); ++muonIt) {
        const int im = muonIt->track().index(); // tried to change this to const index_t (=unsigned int) but leads to segviol
	if (fVerbose > 0) std::cout << "==>HFLambdas> Muon index: " << im << " Ptr: " << muonIt->track().get() << std::endl;
        if (im >= 0) muonIndices.push_back(im);
    }
    if (fVerbose > 0) {
        std::cout << "==>HFLambdas> nMuons = " << hMuons->size() << std::endl;
        std::cout << "==>HFLambdas> nMuonIndices = " << muonIndices.size() << std::endl;
    }

    // Build lists
    trackList_t prList, piList, trackMuonList;
    trackMuonList.reserve(fMaxTracks);
    piList.reserve(fMaxTracks);
    prList.reserve(fMaxTracks);

    for (index_t itrack = 0; itrack < hTracks->size(); ++itrack) {
	reco::TrackBaseRef rTrackView(hTracks, itrack);
	reco::Track tTrack(*rTrackView);
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

    // now do the combinatorics
    for(std::vector<index_t>::const_iterator itm1=muonIndices.begin(); itm1!=muonIndices.end(); itm1++) {
	if (fVerbose > 10) std::cout << "==>HFLambdas>This is muon " << (*itm1) << std::endl;
	for(trackList_t::const_iterator itm2=trackMuonList.begin(); itm2!=trackMuonList.end(); itm2++) {
	    if( (*itm1)!=itm2->first ) { // then we have two distinct muons
		reco::TrackBaseRef tbrMu1(hTracks, (*itm1));
		reco::TrackBaseRef tbrMu2(hTracks, itm2->first);
		reco::Track tMu1(*tbrMu1);
		reco::Track tMu2(*tbrMu2);
		if (tMu1.charge()==tMu2.charge()) continue; // charges must be opposite
		if (fVerbose > 10) std::cout << "==>HFLambdas>This is track muon " << itm2->first << std::endl;
		for(trackList_t::const_iterator itpr=prList.begin(); itpr!=prList.end(); itpr++) {
		    if( (*itm1)!=itpr->first && itm2->first!=itpr->first ) { // prevent from using the same track for two different particles
			if (fVerbose > 10) std::cout << "==>HFLambdas>This is proton " << itpr->first << std::endl;
			for(trackList_t::const_iterator itpi=piList.begin(); itpi!=piList.end(); itpi++) {
			    if( (*itm1)!=itpi->first && itm2->first!=itpi->first && itpr->first!=itpi->first ) { // and check again
				reco::TrackBaseRef tbrPi(hTracks, itpi->first);
				reco::TrackBaseRef tbrPr(hTracks, itpr->first);
				reco::Track tPi(*tbrPi);
				reco::Track tPr(*tbrPr);
				if (tPr.charge()==tPi.charge()) continue; // charges must be opposite
				// now we create a J/Psi

				// find points of closest approach
				TwoTrackMinimumDistance ttmdJpsi = calculatePoca(hTracks, (*itm1), itm2->first);
				TwoTrackMinimumDistance ttmdL0 = calculatePoca(hTracks, itpr->first, itpi->first);
				if (!(ttmdJpsi.status() && ttmdL0.status() )) { 
				    std::cout << "ttmd status invalid" << std::endl; continue;
				} // if something went wrong, discard this pair

				if (fVerbose > 10) {
				    std::cout << "status(): true"
					      << "  distance JPsi: " << ttmdJpsi.distance()
					      << "  distance L0: " << ttmdL0.distance()
					      << std::endl;
				}
				if (ttmdJpsi.distance() > fMaxDoca) continue;
				if (ttmdL0.distance() > fMaxDoca) continue;

				//std::cout << "mag " << ttmdJpsi.crossingPoint().mag() << " " << ttmdJpsi.crossingPoint().perp() << std::endl;
				const double pocaToVtxJPsi=calculateDistToPV(ttmdJpsi.crossingPoint(), fPV);
				const double pocaToVtxL0  =calculateDistToPV(ttmdL0.crossingPoint(), fPV);
				if (pocaToVtxL0<pocaToVtxJPsi) continue; // L0 vtx should be more far away than Jpsi vtx
				if (pocaToVtxJPsi<fMinPocaJpsi) continue;
				if (pocaToVtxL0<fMinPocaL0) continue;

				//reco::TrackBaseRef tbrMu1(hTracks, (*itm1));
				//reco::Track tMu1(*tbrMu1);
				TLorentzVector tlvMu1;
				tlvMu1.SetPtEtaPhiM(tMu1.pt(), tMu1.eta(), tMu1.phi(), MMUON);
				const TLorentzVector tlvJPsi = tlvMu1 + itm2->second;
				// and a Lambda_0
				const TLorentzVector tlvLambda = itpr->second + itpi->second;
				// and now we can make a Lambda_b
				const TLorentzVector tlvLambdaB = tlvJPsi + tlvLambda;
				// check mass windows
				if ( inMassWindow(tlvJPsi.M(), MJPSI, fJPsiWindow)
				     && inMassWindow(tlvLambda.M(), MLAMBDA_0, fL0Window)) {
				    // add candidates
				    const int indxCandJpsi=fillCand(tlvJPsi, ttmdJpsi, (*itm1), itm2->first, 10443);
				    const int indxCandL0=fillCand(tlvLambda, ttmdL0, itpr->first, itpi->first, 13122);
				    const int indxCandLb=fillCand(tlvLambdaB, indxCandJpsi, indxCandL0, 15122);
				    // some output
				    if (fVerbose > 20) {
					std::cout << "==>HFLambdas>" 
					          << "Masses: JPsi: " << tlvJPsi.M()
						  << " Lamdba0: " << tlvLambda.M()
					          << " LambdaB: " << tlvLambdaB.M()
					          << std::endl; 
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

}


// ------------ method called once each job just before starting event loop  ------------
void  HFLambdas::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFLambdas::endJob() {
}

bool HFLambdas::inMassWindow(const value_t& v, const value_t& m, const value_t& dm)
{
    const value_t t = fabs(v-m);
    return (t<dm) ? true : false;
}

int HFLambdas::fillCand(const TLorentzVector& tlvCand, const int& trk1, const int& trk2, const int& type)
{
    TAnaCand *cand = gHFEvent->addCand();
    cand->fPlab = tlvCand.Vect();
    cand->fMass = tlvCand.M();
    cand->fType = type;
    cand->fDau1 = -1;
    cand->fDau2 = -1;
    cand->fSig1 = trk1;
    cand->fSig2 = trk2;
    return gHFEvent->nCands()-1;
}

int HFLambdas::fillCand(const TLorentzVector& tlvCand, const TwoTrackMinimumDistance& ttmd, const int& trk1, const int& trk2, const int& type)
{
    TAnaCand *cand = gHFEvent->addCand();
    cand->fPlab = tlvCand.Vect();
    cand->fMass = tlvCand.M();
    cand->fType = type;
    cand->fDau1 = -1;
    cand->fDau2 = -1;
    cand->fSig1 = trk1;
    cand->fSig2 = trk2;
    cand->fMinDoca = ttmd.distance();
    cand->fPoca = TVector3(ttmd.crossingPoint().x(), ttmd.crossingPoint().y(), ttmd.crossingPoint().z());
    return gHFEvent->nCands()-1;
}

TwoTrackMinimumDistance HFLambdas::calculatePoca(const edm::Handle<edm::View<reco::Track> >& tracks, int track1, int track2)
{
    reco::TrackBaseRef tbr1(tracks, track1);
    reco::Track t1(*tbr1);
    reco::TransientTrack tt1(fTTB->build(t1));

    reco::TrackBaseRef tbr2(tracks, track2);
    reco::Track t2(*tbr2);
    reco::TransientTrack tt2(fTTB->build(t2));

    TwoTrackMinimumDistance ttmd;
    ttmd.calculate(tt1.initialFreeState(),tt2.initialFreeState());

    return ttmd;
    //return std::make_pair(ttmd.crossingPoint(), ttmd.distance());
}

double HFLambdas::calculateDistToPV(const GlobalPoint& pt, const reco::Vertex& vtx)
{
    const TVector3 tv3pt(pt.x(), pt.y(), pt.z());
    const TVector3 tv3vtx(vtx.position().x(), vtx.position().y(), vtx.position().z());
    const TVector3 tv3diff=tv3pt-tv3vtx;
    return tv3diff.Mag();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFLambdas);
