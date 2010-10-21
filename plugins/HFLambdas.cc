#include "HFLambdas.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"

#include<vector>

extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;

class HFLambdaCut : public HFMaxDocaCut {
	public:
		HFLambdaCut(double maxPAngle = 3.2, HFNodeCut *ksCut = NULL, double docaCut = 1E7) :
		  HFMaxDocaCut(docaCut), fKsCut(ksCut), fMaxPAngle(maxPAngle) {}
		virtual bool operator()();
		
		void setKshorCut(HFLambdaCut* newKsCut) {fKsCut = newKsCut;}
		void setMaxPAngle(double newAngle) {fMaxPAngle = newAngle;}
	protected:
		HFNodeCut *fKsCut; // cross reference for pointing angle calculation!
		double fMaxPAngle; // max pointing angle (in radians)
};

bool HFLambdaCut::operator()()
{
  bool result = HFMaxDocaCut::operator()();
  double pAngle;
  TVector3 v;
  
  if(result && fKsCut) {
    
    // check the pointing angle  
    pAngle = fKsCut->fPtCand.Angle(fKsCut->fVtxPos - this->fVtxPos);
    result = pAngle < fMaxPAngle;
  }
  
  return result;
}

HFLambdas::HFLambdas(const ParameterSet& iConfig) :
	fVerbose(iConfig.getUntrackedParameter<int>("verbose",0)),
	fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))),
	fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
	fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
	fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
	fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.4)),
	fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 1)),
	fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow",0.3)),
	fL0Window(iConfig.getUntrackedParameter<double>("ksWindow",0.3)),
	fLbWindow(iConfig.getUntrackedParameter<double>("bdWindow",0.8)),
	fDeltaR(iConfig.getUntrackedParameter<double>("deltaR",99.0)),
	fVertexing(iConfig.getUntrackedParameter<int>("vertexing",1)),
	fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
	fPAngleL0(iConfig.getUntrackedParameter<double>("pAngleKs",0.1))
{
	using namespace std;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "--- HFLambdas constructor" << endl;
	cout << "---  verbose:					" << fVerbose << endl;
	cout << "---  tracksLabel:                              " << fTracksLabel << endl;
	cout << "---  PrimaryVertexLabel:                       " << fPrimaryVertexLabel << endl;
	cout << "---  muonsLabel:				" << fMuonsLabel << endl;
	cout << "---  muonPt:                                   " << fMuonPt << endl;
	cout << "---  pionPt:                                   " << fPionPt << endl;
	cout << "---  psiMuons:					" << fPsiMuons << endl;
	cout << "---  psiWindow:                                " << fPsiWindow << endl;
	cout << "---  fL0Window:                                " << fL0Window << endl;
	cout << "---  bdWindow:                                 " << fLbWindow << endl;
	cout << "---  deltaR:                                   " << fDeltaR << endl;
	cout << "---  vertexing:                                " << fVertexing << endl;
	cout << "---  maxDoca:			       		" << fMaxDoca << endl;
	cout << "---  pAngleKs:                                 " << fPAngleL0 << endl;

} // HFLambdas()

HFLambdas::~HFLambdas()
{}

void HFLambdas::beginJob() {} // beginJob()
void HFLambdas::endJob() {} // endJob()

void HFLambdas::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	// -- get the primary vertex
        Handle<VertexCollection> recoPrimaryVertexCollection;
	iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
	if(!recoPrimaryVertexCollection.isValid()) {
		cout << "==>HFLambdas> No primary vertex collection found, skipping" << endl;
		return;
	}
	
	const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
	if (vertices.size() == 0) {
		cout << "==>HFLambdas> No primary vertex found, skipping" << endl;
		return;
	}
	
	fPV = vertices[gHFEvent->fBestPV];
	if(fVerbose > 0)
		cout << "==>HFLambdas: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
	
	// -- get the collection of muons
	Handle<MuonCollection> hMuons;
	iEvent.getByLabel(fMuonsLabel, hMuons);
	if(!hMuons.isValid()) {
		cout << "==>HFLambdas> No valid MuonCollection with label "<< fMuonsLabel << " found, skipping" << endl;
		return;
	}
	
	// -- get the collection of tracks
	Handle<View<Track> > hTracks;
	iEvent.getByLabel(fTracksLabel, hTracks);
	if(!hTracks.isValid()) {
		cout << "==>HFLambdas> No valid TrackCollection with label "<< fTracksLabel << " found, skipping" << endl;
		return;
	}
	
	// -- Transient tracks for vertexing
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
	if(!fTTB.isValid()) {
		cout << "==>HFLambdas: Error: no TransientTrackBuilder found."<<endl;
		return;
	}
	
	// -- get the collection of muons and store their corresponding track indices
	vector<unsigned int> muonIndices;
	for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
		int im = muon->track().index();
		if(im >= 0) muonIndices.push_back(im);
	}
	
	if (fVerbose > 0) {
		cout << "==>HFLambdas> nMuons = " << hMuons->size() << endl;
		cout << "==>HFLambdas> nMuonIndices = " << muonIndices.size() << endl;
	}
	
	if(muonIndices.size() < (unsigned)fPsiMuons) return;
	
	// -- Build muon list
	TLorentzVector tlv;
	vector<pair<int,TLorentzVector> > pilist,plist,tlist1,tlist2;
	int isMuon(0);
	if(2 == fPsiMuons) {
		tlist1.reserve(10);
		tlist2.reserve(10);
	} else {
		tlist1.reserve(100);
		tlist2.reserve(100);
	}
	pilist.reserve(100);
	plist.reserve(100);
	
	for(unsigned int itrack = 0; itrack < hTracks->size(); ++itrack) {
		TrackBaseRef rTrackView(hTracks,itrack);
		Track tTrack(*rTrackView);

		tlv.SetXYZM(tTrack.px(),tTrack.py(),tTrack.pz(),MPION);
		pilist.push_back(make_pair(itrack,tlv));

		tlv.SetXYZM(tTrack.px(),tTrack.py(),tTrack.pz(),MPROTON);
		plist.push_back(make_pair(itrack,tlv));
		
		tlv.SetXYZM(tTrack.px(),tTrack.py(),tTrack.pz(),MMUON);
		if(2 == fPsiMuons) {
			for (unsigned int im = 0; im < muonIndices.size(); ++im) {
				if(muonIndices[im] == itrack) {
					++isMuon;
					tlist1.push_back(make_pair(itrack,tlv));
					tlist2.push_back(make_pair(itrack,tlv));
				}
			}
		} else if (1 == fPsiMuons) {
			for (unsigned int im = 0; im < muonIndices.size(); ++im) {
				if (muonIndices[im] == itrack) {
					++isMuon;
					tlist1.push_back(make_pair(itrack,tlv));
				}
			}
			tlist2.push_back(make_pair(itrack,tlv));
		} else {
			tlist1.push_back(make_pair(itrack,tlv));
			tlist2.push_back(make_pair(itrack,tlv));
		}
	}
	
	if (isMuon < fPsiMuons) {
		if(fVerbose > 0) cout << "==>HFLambdas> Not enough muons found for J/Psi candidates combinatorics: isMuon = "
			<< isMuon << ", required are " << fPsiMuons << endl;
		return;
	}
	
	HFTwoParticleCombinatorics a(fVerbose);
	vector<pair<int, int> > psiList;
	a.combine(psiList, tlist1, tlist2, 2.8, 3.4, 1);
	if(fVerbose > 0) cout << "==>HFLambdas> J/Psi list size: " << psiList.size() << endl;
	vector<pair<int, int> > l0List;
	a.combine(l0List, pilist, plist, 0.2, 1.4, 0);
	if (fVerbose > 0) cout << "==>HFLambdas> Lambda_0 list size: " << l0List.size() << endl;
	
	// do the vertex fitting...
	HFSequentialVertexFit aSeq(hTracks, fTTB.product(), fPV, fVerbose);
	TLorentzVector tlvPsi, tlvMu1, tlvMu2, tlvPion, tlvProton, tlvLamdba0, tlvLambdaB;

	for (unsigned int i = 0; i < psiList.size(); i++) {
		unsigned int iMuon1 = psiList[i].first;
		unsigned int iMuon2 = psiList[i].second;
		
		TrackBaseRef mu1TrackView(hTracks, iMuon1);
		Track tMuon1(*mu1TrackView);
		if (tMuon1.pt() < fMuonPt) continue;
		tlvMu1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);
		
		TrackBaseRef mu2TrackView(hTracks, iMuon2);
		Track tMuon2(*mu2TrackView);
		if (tMuon2.pt() < fMuonPt) continue;
		tlvMu2.SetPtEtaPhiM(tMuon2.pt(),tMuon2.eta(),tMuon2.phi(), MMUON);
		
		tlvPsi = tlvMu1 + tlvMu2;
		if( (TMath::Abs(tlvPsi.M() - MJPSI) > fPsiWindow) ) continue;
		
		for(unsigned int j = 0; j < l0List.size(); j++) {
			
			unsigned int iPion = l0List[j].first;
			unsigned int iProton = l0List[j].second;
			
			if (iPion == iMuon1 || iPion == iMuon2) continue;
			if (iProton == iMuon1 || iProton == iMuon2) continue;
			
			// get lorentz vector of first pion
			TrackBaseRef pionTrackView(hTracks,iPion);
			Track tPion(*pionTrackView);
			if(tPion.pt() < fPionPt) continue;
			tlvPion.SetXYZM(tPion.px(),tPion.py(),tPion.pz(),MPION);
			if(tlvPsi.DeltaR(tlvPion) > fDeltaR) continue;
			
			// get lorentz vector of second pion
			TrackBaseRef protonTrackView(hTracks,iProton);
			Track tProton(*protonTrackView);
			if(tProton.pt() < fPionPt) continue;
			tlvProton.SetXYZM(tProton.px(),tProton.py(),tProton.pz(),MPROTON);
			if(tlvPsi.DeltaR(tlvProton) > fDeltaR) continue;
			
			tlvLamdba0 = tlvPion + tlvProton;
			if (TMath::Abs(tlvLamdba0.M() - MLAMBDA_0) > fL0Window) continue;
			
			tlvLambdaB = tlvPsi + tlvLamdba0;
			if (TMath::Abs(tlvLambdaB.M() - MLAMBDA_B) > fLbWindow) continue;

			// without J/Psi mass constraint
			HFDecayTree theTree(605122);
			
			HFDecayTreeIterator iterator = theTree.addDecayTree(600443,0); // no vertexing & no mass constraint
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(603122,1); // no mass constraint, with own vertex of Lamda0
			iterator->addTrack(iPion,211);
			iterator->addTrack(iProton,2212);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngleL0, &(*iterator->getNodeCut()),fMaxDoca)));
			
			aSeq.doFit(&theTree);
			
			// without J/Psi mass constraint, but with an own J/Psi vertex
			theTree.clear();
			theTree.particleID = 705122;
			
			iterator = theTree.addDecayTree(700443,1); // vertexing but no mass constraint...
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(703122,1); // Lambda0 with vertexing
			iterator->addTrack(iPion,211);
			iterator->addTrack(iProton,2212);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngleL0, &(*iterator->getNodeCut()), fMaxDoca)));
			
			aSeq.doFit(&theTree);

			// with J/Psi mass constraint                                                                                                    
			theTree.clear();
			theTree.particleID = 805122;
			
			iterator = theTree.addDecayTree(800443,1,MJPSI); // J/Psi needs an own vertex!!
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(803122);
			iterator->addTrack(iPion,211);
			iterator->addTrack(iProton,2212);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngleL0,&(*iterator->getNodeCut()),fMaxDoca)));
			
			aSeq.doFit(&theTree);

			// with mass constraint for J/Psi and Lambda0
			theTree.clear();
			theTree.particleID = 905122;
			
			iterator = theTree.addDecayTree(900443,1,MJPSI); 
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(903122,1,MLAMBDA_0);
			iterator->addTrack(iPion,211);
			iterator->addTrack(iProton,2212);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngleL0,&(*iterator->getNodeCut()),fMaxDoca)));
			
			aSeq.doFit(&theTree);
		}
	}
} // analyze()

// define this as a plug-in
DEFINE_FWK_MODULE(HFLambdas);
