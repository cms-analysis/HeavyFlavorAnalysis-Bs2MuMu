#include <iostream>

#include "HFB2muD0.hh"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;

#define MMUON 0.10566
#define MKAON 0.49368
#define MPION 0.139570
#define MB_0 5.27953
#define MBPLUS 5.27927
#define MD_0 1.86484
#define MD_S 2.01027


// ----------------------------------------------------------------------
HFB2muD0::HFB2muD0(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", string("goodTracks"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel",
							      string("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
  fProtonPt(iConfig.getUntrackedParameter<double>("protonPt", 1.0)),
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.4)),
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.)),
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 1.0)),
  fDeltaMD0(iConfig.getUntrackedParameter<double>("deltaMD0", 0.4)),
  fDeltaMDs(iConfig.getUntrackedParameter<double>("deltaMDs", 0.4)),
  fSlowPt(iConfig.getUntrackedParameter<double>("slowPt", 0.1)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- " << __FILE__ << endl;
  cout << "--- $Revision: 1.3 $" << endl;
  cout << "--- HFB2muD0 constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  protonPt:                 " << fProtonPt << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  slowPt:                   " << fSlowPt << endl;
  cout << "---  pionPt:                   " << fPionPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  deltaMD0:                 " << fDeltaMD0 << endl;
  cout << "---  deltaMDs:                 " << fDeltaMDs << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFB2muD0::~HFB2muD0() {
  
}


// ----------------------------------------------------------------------
void HFB2muD0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // -- get the primary vertex
  
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    if (fVerbose > 0) cout << "==>HFB2muD0> No primary vertex found, skipping" << endl;
    return;
  }
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    if (fVerbose > 0) cout << "==>HFB2muD0> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fBestPV]; 
  if (fVerbose > 1) {
    cout << "==>HFB2muD0> Taking vertex " << gHFEvent->fBestPV
	 << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of muons
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  if (!hMuons.isValid()) {
    if (fVerbose > 0)
      cout << "==>HFB2muD0> No valid MuonCollection with label " << fMuonsLabel
	   << " found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  edm::Handle<edm::View<reco::Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    if (fVerbose > 0)
      cout << "==>HFB2muD0> No valid TrackCollection with label " << fTracksLabel
	   << " found, skipping" << endl;
    return;
  }

  if (hTracks->size() > 100) {
    if (fVerbose > 0)
      cout << "==>HFB2muD0> TrackCollection with label " << fTracksLabel
	   << " has " << hTracks->size() << " entries, skipping event!!" << endl;
  }
 
  // -- Fill muons into muonTracks
  if (fVerbose > 1)
    cout << "==>HFB2muD0> nMuons = " << hMuons->size() << endl;


  // -- get the collection of muons and store their corresponding track indices
  vector<int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index();
    if (im > 0) muonIndices.push_back(im);
  }
  if (fVerbose > 1)
    cout << "==>HFB2muD0> nMuonIndices = " << muonIndices.size() << endl;

  if (muonIndices.size() < 1) return;

  // =========================================================================================================

  int nD0, nDs;
  nD0 = nDs = 0;

  // go through muon track list
  TLorentzVector vmuon, vtrack1, vtrack2, vd0, vkaon, vpion;
  TLorentzVector vslow, vds;
  for (unsigned int muon_idx = 0; muon_idx < muonIndices.size(); muon_idx++) {
    unsigned int imuon = muonIndices[muon_idx];
    if (fVerbose > 1)
      cout << "==>HFB2muD0> Checking Muon #" << muon_idx + 1 << endl;

    TrackBaseRef muTrackView(hTracks, imuon); // <-------------------------------------------------------- ???
    Track muon(*muTrackView);

    vmuon.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), MMUON);

    if (muon.pt() < fMuonPt) {
      if (fVerbose > 1)
	cout << "==>HFB2muD0> muon.pt() < " << fMuonPt << endl;
      continue;
    }
    
    vector<unsigned int> track_list;
    vector<unsigned int> slow_list;
    for (unsigned int itrack = 0; itrack < hTracks->size(); itrack++) {
      if (itrack == imuon)
	continue;
      
      TrackBaseRef tr1TrackView(hTracks, itrack);
      Track track(*tr1TrackView);

      vtrack1.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), MPION);
      
      if (vmuon.DeltaR(vtrack1) < fDeltaR) {
	if (track.pt() > fTrackPt)
	  track_list.push_back(itrack);
	if (track.pt() > fSlowPt)
	  slow_list.push_back(itrack);
      }
    }
    
    if (fVerbose > 1)
      cout << "==>HFB2muD0> track_list size = " << track_list.size() << endl;

    if (track_list.size() < 2) {
      if (fVerbose > 1)
	cout << "==>HFB2muD0> Too little tracks in muon cone!" << endl;
      continue;
    }
    
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
    if (!fTTB.isValid()) {
      if (fVerbose > 0)
	cout << " -->HFB2muD0: Error: no TransientTrackBuilder found." << endl;
      return;
    }
    
    // D0 reconstruction #####################################################################################

    TAnaTrack * k_track, * p_track, * m_track, * s_track;
    vector<unsigned int>::iterator ia, ib;
    for (ia = track_list.begin(); ia != track_list.end(); ++ia) {
      for (ib = ia + 1; ib != track_list.end(); ++ib) {
	TrackBaseRef tr1TrackView(hTracks, *ia);
	Track a(*tr1TrackView);
	TrackBaseRef tr2TrackView(hTracks, *ib);
	Track b(*tr2TrackView);

	if (a.charge() * b.charge() < 0) {
	  Track kaon, pion;
	  unsigned int ikaon, ipion;
	  if (a.charge() * muon.charge() > 0) {
	    // a is the kaon
	    kaon = a;
	    pion = b;
	    ikaon = *ia;
	    ipion = *ib;
	    vkaon.SetPtEtaPhiM(a.pt(), a.eta(), a.phi(), MKAON);
	    vpion.SetPtEtaPhiM(b.pt(), b.eta(), b.phi(), MPION);
	  } else {
	    // b is the kaon
	    kaon = b;
	    pion = a;
	    ikaon = *ib;
	    ipion = *ia;
	    vpion.SetPtEtaPhiM(a.pt(), a.eta(), a.phi(), MPION);
	    vkaon.SetPtEtaPhiM(b.pt(), b.eta(), b.phi(), MKAON);
	  }
	  
	  vd0 = vkaon + vpion;
	  if (TMath::Abs(vd0.M()- MD_0) < fDeltaMD0) {
	    if (fVerbose > 2)
	      cout << "==>HFB2muD0> D0 candidate found! Mass: " << vd0.M() << endl;
	    nD0++;

	    std::vector<reco::Track> fitTracks;
	    fitTracks.clear();
	    fitTracks.push_back(kaon);
	    fitTracks.push_back(pion);

	    //m1.SetXYZM(tMuon1.px(), tMuon1.py(), tMuon1.pz(), MMUON);
	    //m2.SetXYZM(tMuon2.px(), tMuon2.py(), tMuon2.pz(), MMUON);

	    TAnaVertex D0Vtx;
	    D0Vtx = DoVertexFit(fitTracks);
	    D0Vtx.addTrack(ikaon);
	    D0Vtx.addTrack(ipion);

 	    //vkaon.SetXYZM(fitTracks[0].px(), fitTracks[0].py(), fitTracks[0].pz(), MKAON);
 	    //vpion.SetXYZM(fitTracks[1].px(), fitTracks[1].py(), fitTracks[1].pz(), MPION);

	    int id0 = gHFEvent->nCands();
	    TAnaCand *d0_cand;
	    d0_cand        = gHFEvent->addCand();
	    d0_cand->fPlab = vd0.Vect();
	    d0_cand->fMass = vd0.M();
	    d0_cand->fType = 421;
	    d0_cand->fVtx  = D0Vtx;
	    d0_cand->fSig1 = gHFEvent->nSigTracks();
	    d0_cand->fSig2 = d0_cand->fSig1 + 2;

	    k_track            = gHFEvent->addSigTrack();
	    k_track->fMCID     = kaon.charge() * 321;
	    k_track->fGenIndex = -1;
	    k_track->fQ        = kaon.charge();
	    k_track->fPlab.SetPtEtaPhi(vkaon.Pt(), vkaon.Eta(), vkaon.Phi());
	    k_track->fIndex    = ikaon;

	    p_track            = gHFEvent->addSigTrack();
	    p_track->fMCID     = pion.charge() * 211;
	    p_track->fGenIndex = -1;
	    p_track->fQ        = pion.charge();
	    p_track->fPlab.SetPtEtaPhi(vpion.Pt(), vpion.Eta(), vpion.Phi());
	    p_track->fIndex    = ipion;

	    m_track            = gHFEvent->addSigTrack();
	    m_track->fMCID     = - muon.charge() * 13;
	    m_track->fGenIndex = -1;
	    m_track->fQ        = muon.charge();
	    m_track->fPlab.SetPtEtaPhi(vmuon.Pt(), vmuon.Eta(), vmuon.Phi());
	    m_track->fIndex    = imuon;

	    // D* reconstruction #############################################################################
	    vector<unsigned int>::iterator islow;
	    for (islow = slow_list.begin(); islow != slow_list.end(); ++islow) {
	      if (*islow == ikaon || *islow == ipion)
		continue;
	      TrackBaseRef tr1TrackView(hTracks, *islow);
	      Track slow(*tr1TrackView);
	      if (muon.charge() * slow.charge() < 0) {
		vslow.SetPtEtaPhiM(slow.pt(), slow.eta(), slow.phi(), MPION);
		vds = vslow + vkaon + vpion;

		// D* mass cut
		if (TMath::Abs(vds.M() - vd0.M() + MD_0 - MD_S) < fDeltaMDs) {
		  if (fVerbose > 2)
		    cout << "==>HFB2muD0> D* candidate found! Mass: " << vds.M() << endl;
		  nDs++;
		  
		  // Fill TAnaCand
		  std::vector<reco::Track> fitTracks;
		  fitTracks.clear();
		  fitTracks.push_back(muon);
		  fitTracks.push_back(slow);
		  // fitTracks.push_back(D0 candidate);  <------------------------------------------------ !!!
		  
		  //m1.SetXYZM(tMuon1.px(), tMuon1.py(), tMuon1.pz(), MMUON);
		  //m2.SetXYZM(tMuon2.px(), tMuon2.py(), tMuon2.pz(), MMUON);
		  
		  TAnaVertex DsVtx;
		  DsVtx = DoVertexFit(fitTracks);
		  DsVtx.addTrack(*islow);
		  DsVtx.addTrack(imuon);
		  
		  //vkaon.SetXYZM(fitTracks[0].px(), fitTracks[0].py(), fitTracks[0].pz(), MKAON);
		  //vpion.SetXYZM(fitTracks[1].px(), fitTracks[1].py(), fitTracks[1].pz(), MPION);
		  
		  TAnaCand *ds_cand;
		  ds_cand        = gHFEvent->addCand();
		  ds_cand->fPlab = vds.Vect();
		  ds_cand->fMass = vds.M();
		  ds_cand->fType = 413;
		  ds_cand->fVtx  = DsVtx;
		  ds_cand->fDau1  = id0;
		  ds_cand->fDau2  = ds_cand->fDau1;
		  ds_cand->fSig1 = gHFEvent->nSigTracks();
		  ds_cand->fSig2 = ds_cand->fSig1 + 3;

		  k_track            = gHFEvent->addSigTrack();
		  k_track->fMCID     = kaon.charge() * 321;
		  k_track->fGenIndex = -1;
		  k_track->fQ        = kaon.charge();
		  k_track->fPlab.SetPtEtaPhi(vkaon.Pt(), vkaon.Eta(), vkaon.Phi());
		  k_track->fIndex    = ikaon;
		  
		  p_track            = gHFEvent->addSigTrack();
		  p_track->fMCID     = pion.charge() * 211;
		  p_track->fGenIndex = -1;
		  p_track->fQ        = pion.charge();
		  p_track->fPlab.SetPtEtaPhi(vpion.Pt(), vpion.Eta(), vpion.Phi());
		  p_track->fIndex    = ipion;

		  s_track            = gHFEvent->addSigTrack();
		  s_track->fMCID     = slow.charge() * 211;
		  s_track->fGenIndex = -1;
		  s_track->fQ        = slow.charge();
		  s_track->fPlab.SetPtEtaPhi(vslow.Pt(), vslow.Eta(), vslow.Phi());
		  s_track->fIndex    = *islow;

		  m_track            = gHFEvent->addSigTrack();
		  m_track->fMCID     = - muon.charge() * 13;
		  m_track->fGenIndex = -1;
		  m_track->fQ        = muon.charge();
		  m_track->fPlab.SetPtEtaPhi(vmuon.Pt(), vmuon.Eta(), vmuon.Phi());
		  m_track->fIndex    = imuon;

		} // mass cut
	      } // charge cut D*
	    } // track c
	  } // mass cut
	} // charge cut D0
      } // track b (D0 reco)
    } // track a (D0 reco)
  } // muons

  if (fVerbose > 0 && nD0 > 0)
    cout << "==>HFB2muD0> " << nD0 << " D0 candidates found." << endl;
  if (fVerbose > 0 && nDs > 0)
    cout << "==>HFB2muD0> " << nDs << " D* candidates found." << endl;

  return;
}


// ----------------------------------------------------------------------
TAnaVertex HFB2muD0::DoVertexFit(std::vector<reco::Track> &Tracks){
  KalmanVertexFitter theFitter(true);
  std::vector<reco::TransientTrack> TransTracks;
  for(reco::TrackCollection::const_iterator iTrack=Tracks.begin(); iTrack!=Tracks.end(); iTrack++)
    TransTracks.push_back(fTTB->build(*iTrack));

  TransientVertex TransSecVtx = theFitter.vertex(TransTracks);
  TAnaVertex anaVt;
  if ( TransSecVtx.isValid() ) {
    if ( isnan(TransSecVtx.position().x())
	 || isnan(TransSecVtx.position().y())
	 || isnan(TransSecVtx.position().z()) ) {
      if (fVerbose > 0) cout << "==>HFB2muD0> Something went wrong! SecVtx nan - Aborting... !" << endl;
      anaVt.fPoint.SetXYZ(0,0,0);
      return anaVt;
    }
  } else {
    if (fVerbose > 0) cout << "==>HFB2muD0> KVF failed! Aborting... !" << endl;
    anaVt.fPoint.SetXYZ(0,0,0);
    return anaVt;
  }

  // -- Give back refitted tracks
  Tracks.clear();
  std::vector<reco::TransientTrack> refTT = TransSecVtx.refittedTracks();
  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); i != refTT.end(); i++) {
    const Track & ftt = i->track();
    Tracks.push_back(ftt);
  }

  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());

  // -- Build vertex for ntuple
  anaVt.setInfo(chi.value(), int(chi.degreesOfFreedom()), chi.probability(), 1, 0);
  anaVt.fPoint.SetXYZ(TransSecVtx.position().x(),
		      TransSecVtx.position().y(),
		      TransSecVtx.position().z());
 
  // -- Distance to primary vertex

  VertexDistanceXY axy;
  double dXY      = axy.distance(fPV, TransSecVtx).value();
  double dXYE     = axy.distance(fPV, TransSecVtx).error();
 
  VertexDistance3D a3d;
  double d3d      = a3d.distance(fPV, TransSecVtx).value();
  double d3dE     = a3d.distance(fPV, TransSecVtx).error();
 
  anaVt.fDxy  = dXY;
  anaVt.fDxyE = dXYE;

  anaVt.fD3d  = d3d;
  anaVt.fD3dE = d3dE;

  return anaVt;
}


// ------------ method called once each job just before starting event loop  ------------
void  HFB2muD0::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFB2muD0::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFB2muD0);
