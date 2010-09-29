#include "HFBd2JpsiKs.h"
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

class HFKshortCut : public HFMaxDocaCut {
	public:
		HFKshortCut(double maxPAngle = 3.2, HFNodeCut *ksCut = NULL, double docaCut = 1E7) :
		  HFMaxDocaCut(docaCut), fKsCut(ksCut), fMaxPAngle(maxPAngle) {}
		virtual bool operator()();
		
		void setKshorCut(HFKshortCut* newKsCut) {fKsCut = newKsCut;}
		void setMaxPAngle(double newAngle) {fMaxPAngle = newAngle;}
	protected:
		HFNodeCut *fKsCut; // cross reference for pointing angle calculation!
		double fMaxPAngle; // max pointing angle (in radians)
};

bool HFKshortCut::operator()()
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

HFBd2JpsiKs::HFBd2JpsiKs(const ParameterSet& iConfig) :
	fVerbose(iConfig.getUntrackedParameter<int>("verbose",0)),
	fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))),
	fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
	fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
	fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
	fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.4)),
	fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
	fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow",0.3)),
	fKsWindow(iConfig.getUntrackedParameter<double>("ksWindow",0.3)),
	fBdWindow(iConfig.getUntrackedParameter<double>("bdWindow",0.8)),
	fDeltaR(iConfig.getUntrackedParameter<double>("deltaR",99.0)),
	fVertexing(iConfig.getUntrackedParameter<int>("vertexing",1)),
	fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
	fPAngleKs(iConfig.getUntrackedParameter<double>("pAngleKs",0.1))
{
	using namespace std;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "--- HFBd2JpsiKs constructor" << endl;
	cout << "---  verbose:					" << fVerbose << endl;
	cout << "---  tracksLabel:                              " << fTracksLabel << endl;
	cout << "---  PrimaryVertexLabel:                       " << fPrimaryVertexLabel << endl;
	cout << "---  muonsLabel:				" << fMuonsLabel << endl;
	cout << "---  muonPt:                                   " << fMuonPt << endl;
	cout << "---  pionPt:                                   " << fPionPt << endl;
	cout << "---  psiMuons:					" << fPsiMuons << endl;
	cout << "---  psiWindow:                                " << fPsiWindow << endl;
	cout << "---  fKsWindow:                                " << fKsWindow << endl;
	cout << "---  bdWindow:                                 " << fBdWindow << endl;
	cout << "---  deltaR:                                   " << fDeltaR << endl;
	cout << "---  vertexing:                                " << fVertexing << endl;
	cout << "---  maxDoca:			       		" << fMaxDoca << endl;
	cout << "---  pAngleKs:                                 " << fPAngleKs << endl;

} // HFBd2JpsiKs()

HFBd2JpsiKs::~HFBd2JpsiKs()
{}

void HFBd2JpsiKs::beginJob() {} // beginJob()
void HFBd2JpsiKs::endJob() {} // endJob()

void HFBd2JpsiKs::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	// -- get the primary vertex
        Handle<VertexCollection> recoPrimaryVertexCollection;
	iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
	if(!recoPrimaryVertexCollection.isValid()) {
		cout << "==>HFBd2JpsiKs> No primary vertex collection found, skipping" << endl;
		return;
	}
	
	const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
	if (vertices.size() == 0) {
		cout << "==>HFBd2JpsiKs> No primary vertex found, skipping" << endl;
		return;
	}
	
	fPV = vertices[gHFEvent->fBestPV];
	if(fVerbose > 0)
		cout << "==>HFBd2JpsiKs: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
	
	// -- get the collection of muons
	Handle<MuonCollection> hMuons;
	iEvent.getByLabel(fMuonsLabel, hMuons);
	if(!hMuons.isValid()) {
		cout << "==>HFBd2JpsiKs> No valid MuonCollection with label "<< fMuonsLabel << " found, skipping" << endl;
		return;
	}
	
	// -- get the collection of tracks
	Handle<View<Track> > hTracks;
	iEvent.getByLabel(fTracksLabel, hTracks);
	if(!hTracks.isValid()) {
		cout << "==>HFBd2JpsiKs> No valid TrackCollection with label "<< fTracksLabel << " found, skipping" << endl;
		return;
	}
	
	// -- Transient tracks for vertexing
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
	if(!fTTB.isValid()) {
		cout << "==>HFBd2JpsiKs: Error: no TransientTrackBuilder found."<<endl;
		return;
	}
	
	// -- get the collection of muons and store their corresponding track indices
	vector<unsigned int> muonIndices;
	for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
		int im = muon->track().index();
		if(im >= 0) muonIndices.push_back(im);
	}
	
	if (fVerbose > 0) {
		cout << "==>HFBd2JpsiKs> nMuons = " << hMuons->size() << endl;
		cout << "==>HFBd2JpsiKs> nMuonIndices = " << muonIndices.size() << endl;
	}
	
	if(muonIndices.size() < (unsigned)fPsiMuons) return;
	
	// -- Build muon list
	TLorentzVector tlv;
	vector<pair<int,TLorentzVector> > pilist1,pilist2,tlist1,tlist2;
	int isMuon(0);
	if(2 == fPsiMuons) {
		tlist1.reserve(10);
		tlist2.reserve(10);
	} else {
		tlist1.reserve(100);
		tlist2.reserve(100);
	}
	pilist1.reserve(100);
	pilist2.reserve(100);
	
	for(unsigned int itrack = 0; itrack < hTracks->size(); ++itrack) {
		TrackBaseRef rTrackView(hTracks,itrack);
		Track tTrack(*rTrackView);
		tlv.SetXYZM(tTrack.px(),tTrack.py(),tTrack.pz(),MPION);
		pilist1.push_back(make_pair(itrack,tlv));
		pilist2.push_back(make_pair(itrack,tlv));
		
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
		if(fVerbose > 0) cout << "==>HFBd2JpsiKs> Not enough muons found for J/Psi candidates combinatorics: isMuon = "
			<< isMuon << ", required are " << fPsiMuons << endl;
		return;
	}
	
	HFTwoParticleCombinatorics a(fVerbose);
	vector<pair<int, int> > psiList;
	a.combine(psiList, tlist1, tlist2, 2.8, 3.4, 1);
	if(fVerbose > 0) cout << "==>HFBd2JpsiKs> J/Psi list size: " << psiList.size() << endl;
	vector<pair<int, int> > ksList;
	a.combine(ksList, pilist1, pilist2, 0.2, 1.0, 1);
	if (fVerbose > 0) cout << "==>HFBd2JpsiKs> Ks list size: " << ksList.size() << endl;
	
	// do the vertex fitting...
	HFSequentialVertexFit aSeq(hTracks, fTTB.product(), fPV, fVerbose);
	TLorentzVector psi,m1,m2,pi1,pi2,ks,bd;

	for (unsigned int i = 0; i < psiList.size(); i++) {
		unsigned int iMuon1 = psiList[i].first;
		unsigned int iMuon2 = psiList[i].second;
		
		TrackBaseRef mu1TrackView(hTracks, iMuon1);
		Track tMuon1(*mu1TrackView);
		if (tMuon1.pt() < fMuonPt) continue;
		m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);
		
		TrackBaseRef mu2TrackView(hTracks, iMuon2);
		Track tMuon2(*mu2TrackView);
		if (tMuon2.pt() < fMuonPt) continue;
		m2.SetPtEtaPhiM(tMuon2.pt(),tMuon2.eta(),tMuon2.phi(), MMUON);
		
		psi = m1 + m2;
		if( (TMath::Abs(psi.M() - MJPSI) > fPsiWindow) ) continue;
		
		for(unsigned int j = 0; j < ksList.size(); j++) {
			
			unsigned int iPion1 = ksList[j].first;
			unsigned int iPion2 = ksList[j].second;
			
			if (iPion1 == iMuon1 || iPion1 == iMuon2) continue;
			if (iPion2 == iMuon1 || iPion2 == iMuon2) continue;
			
			// get lorentz vector of first pion
			TrackBaseRef pi1TrackView(hTracks,iPion1);
			Track tPi1(*pi1TrackView);
			if(tPi1.pt() < fPionPt) continue;
			pi1.SetXYZM(tPi1.px(),tPi1.py(),tPi1.pz(),MPION);
			if(psi.DeltaR(pi1) > fDeltaR) continue;
			
			// get lorentz vector of second pion
			TrackBaseRef pi2TrackView(hTracks,iPion2);
			Track tPi2(*pi2TrackView);
			if(tPi2.pt() < fPionPt) continue;
			pi2.SetXYZM(tPi2.px(),tPi2.py(),tPi2.pz(),MPION);
			if(psi.DeltaR(pi2) > fDeltaR) continue;
			
			ks = pi1 + pi2;
			if (TMath::Abs(ks.M() - MKSHORT) > fKsWindow) continue;
			
			bd = psi + ks;
			if (TMath::Abs(bd.M() - MB_0) > fBdWindow) continue;

			// without J/Psi mass constraint
			HFDecayTree theTree(600511);
			
			HFDecayTreeIterator iterator = theTree.addDecayTree(600443,0); // no vertexing & no mass constraint
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(600310,1); // no mass constraint, with own vertex of Ks
			iterator->addTrack(iPion1,211);
			iterator->addTrack(iPion2,211);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			theTree.setNodeCut(RefCountedHFNodeCut(new HFKshortCut(fPAngleKs, &(*iterator->getNodeCut()),fMaxDoca)));
			
			aSeq.doFit(&theTree);
			
			// without J/Psi mass constraint, but with an own J/Psi vertex
			theTree.clear();
			theTree.particleID = 700511;
			
			iterator = theTree.addDecayTree(700443,1); // vertexing but no mass constraint...
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(700310,1); // Ks with vertexing
			iterator->addTrack(iPion1,211);
			iterator->addTrack(iPion2,211);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			theTree.setNodeCut(RefCountedHFNodeCut(new HFKshortCut(fPAngleKs, &(*iterator->getNodeCut()), fMaxDoca)));
			
			aSeq.doFit(&theTree);
			
			// with J/Psi mass constraint                                                                                                    
			theTree.clear();
			theTree.particleID = 800511;
			
			iterator = theTree.addDecayTree(800443,1,MJPSI); // J/Psi needs an own vertex!!
			iterator->addTrack(iMuon1,13);
			iterator->addTrack(iMuon2,13);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
			iterator = theTree.addDecayTree(800310);
			iterator->addTrack(iPion1,211);
			iterator->addTrack(iPion2,211);
			iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			theTree.setNodeCut(RefCountedHFNodeCut(new HFKshortCut(fPAngleKs,&(*iterator->getNodeCut()),fMaxDoca)));
			
			aSeq.doFit(&theTree);
		}
	}
} // analyze()

// define this as a plug-in
DEFINE_FWK_MODULE(HFBd2JpsiKs);
