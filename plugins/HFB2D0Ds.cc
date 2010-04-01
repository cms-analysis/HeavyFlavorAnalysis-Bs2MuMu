#include <iostream>

#include "HFB2D0Ds.hh"

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
HFB2D0Ds::HFB2D0Ds(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", string("goodTracks"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel",
							      string("offlinePrimaryVertices"))),
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.5)),
  fKaonPt(iConfig.getUntrackedParameter<double>("kaonPt", 0.5)),
  fDeltaMD0(iConfig.getUntrackedParameter<double>("deltaMD0", 0.4)),
  fDeltaMDs(iConfig.getUntrackedParameter<double>("deltaMDs", 0.4)),
  fSlowPt(iConfig.getUntrackedParameter<double>("slowPt", 0.1)),
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 3.2)),
  fVtxErr(iConfig.getUntrackedParameter<double>("vtxE", 20.0)),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxP", 0.0)),
  fVtxAngle(iConfig.getUntrackedParameter<double>("vtxA", 5.0)),
  fVtxRDist(iConfig.getUntrackedParameter<double>("vtxRD", 0.0)),
  fVtxADist(iConfig.getUntrackedParameter<double>("vtxAD", 0.0)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- " << __FILE__ << endl;
  cout << "--- $Revision: 1.1 $" << endl;
  cout << "--- HFB2D0Ds constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  slowPt:                   " << fSlowPt << endl;
  cout << "---  pionPt:                   " << fPionPt << endl;
  cout << "---  kaonPt:                   " << fKaonPt << endl;
  cout << "---  deltaMD0:                 " << fDeltaMD0 << endl;
  cout << "---  deltaMDs:                 " << fDeltaMDs << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  vtxE:                     " << fVtxErr << endl;
  cout << "---  vtxP:                     " << fVtxProb << endl;
  cout << "---  vtxA:                     " << fVtxAngle << endl;
  cout << "---  vtxRD:                    " << fVtxRDist << endl;
  cout << "---  vtxAD:                    " << fVtxADist << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFB2D0Ds::~HFB2D0Ds() {
  
}


// ----------------------------------------------------------------------
void HFB2D0Ds::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // -- get the primary vertex
  
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    if (fVerbose > 0) cout << "==>HFB2D0Ds> No primary vertex found, skipping" << endl;
    return;
  }
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    if (fVerbose > 0) cout << "==>HFB2D0Ds> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fEventTag]; 
  if (fVerbose > 1) {
    cout << "==>HFB2D0Ds> Taking vertex " << gHFEvent->fEventTag
	 << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of tracks
  edm::Handle<edm::View<reco::Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    if (fVerbose > 0)
      cout << "==>HFB2D0Ds> No valid TrackCollection with label " << fTracksLabel
	   << " found, skipping" << endl;
    return;
  }

  if (hTracks->size() > 100) {
    if (fVerbose > 0)
      cout << "==>HFB2D0Ds> TrackCollection with label " << fTracksLabel
	   << " has " << hTracks->size() << " entries, skipping event!!" << endl;
    return;
  }
 
  // =========================================================================================================

  int nD0, nDs;
  nD0 = nDs = 0;

  // go through muon track list
  TLorentzVector vtrack1, vtrack2, vd0, vkaon, vpion;
  TLorentzVector vslow, vds;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    if (fVerbose > 0)
      cout << " -->HFB2D0Ds: Error: no TransientTrackBuilder found." << endl;
    return;
  }
    
  // D0 reconstruction #####################################################################################

  TAnaTrack * k_track, * p_track, * s_track;
  unsigned int ia, ib;
  for (ia = 0; ia < hTracks->size(); ++ia) {
    for (ib = 0; ib < hTracks->size(); ++ib) {
      if (ia == ib)
	continue;
      TrackBaseRef tr1TrackView(hTracks, ia);
      Track a(*tr1TrackView);
      TrackBaseRef tr2TrackView(hTracks, ib);
      Track b(*tr2TrackView);

      if (a.charge() * b.charge() < 0) {
	Track kaon, pion;
	unsigned int ikaon, ipion;
	kaon = a;
	pion = b;
	ikaon = ia;
	ipion = ib;
	vkaon.SetPtEtaPhiM(a.pt(), a.eta(), a.phi(), MKAON);
	vpion.SetPtEtaPhiM(b.pt(), b.eta(), b.phi(), MPION);
	  
	vd0 = vkaon + vpion;

	// cuts
	if (vkaon.Pt() < fKaonPt || vpion.Pt() < fPionPt)
	  continue;

	if (vkaon.DeltaR(vpion) > fDeltaR)
	  continue;

	if (TMath::Abs(vd0.M()- MD_0) < fDeltaMD0) {

	  std::vector<reco::Track> fitTracks;
	  fitTracks.clear();
	  fitTracks.push_back(kaon);
	  fitTracks.push_back(pion);

	  //m1.SetXYZM(tMuon1.px(), tMuon1.py(), tMuon1.pz(), MMUON);
	  //m2.SetXYZM(tMuon2.px(), tMuon2.py(), tMuon2.pz(), MMUON);

	  TAnaVertex D0Vtx;
	  D0Vtx.fStatus = 0;
	  D0Vtx = DoVertexFit(fitTracks);
	  D0Vtx.addTrack(ikaon);
	  D0Vtx.addTrack(ipion);

	  if (D0Vtx.fStatus == 0)
	    continue;

	  //vkaon.SetXYZM(fitTracks[0].px(), fitTracks[0].py(), fitTracks[0].pz(), MKAON);
	  //vpion.SetXYZM(fitTracks[1].px(), fitTracks[1].py(), fitTracks[1].pz(), MPION);

	  // vertex cuts
	  if (D0Vtx.fD3dE > fVtxErr)
	    continue;
	  if (D0Vtx.fProb < fVtxProb)
	    continue;
	  if (D0Vtx.fD3d / D0Vtx.fD3dE < fVtxRDist)
	    continue;
	  if (D0Vtx.fD3d < fVtxADist)
	    continue;
	  if (vd0.Angle(D0Vtx.fPoint - gHFEvent->getPV(gHFEvent->fEventTag)->fPoint) > fVtxAngle)
	    continue;

	  if (fVerbose > 1)
	    cout << "==>HFB2D0Ds> D0 candidate found! Mass: " << vd0.M() << endl;
	  nD0++;

	  int id0 = gHFEvent->nCands();
	  TAnaCand *d0_cand;
	  d0_cand        = gHFEvent->addCand();
	  d0_cand->fPlab = vd0.Vect();
	  d0_cand->fMass = vd0.M();
	  d0_cand->fType = 421;
	  d0_cand->fVtx  = D0Vtx;
	  d0_cand->fSig1 = gHFEvent->nSigTracks();
	  d0_cand->fSig2 = d0_cand->fSig1 + 1;

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

	  // D* reconstruction #############################################################################
	  for (unsigned int islow = 0; islow < hTracks->size(); ++islow) {
	    if (islow == ikaon || islow == ipion)
	      continue;
	    TrackBaseRef tr1TrackView(hTracks, islow);
	    Track slow(*tr1TrackView);
	    if (kaon.charge() * slow.charge() < 0) {
	      vslow.SetPtEtaPhiM(slow.pt(), slow.eta(), slow.phi(), MPION);
	      vds = vslow + vkaon + vpion;

	      // cuts
	      if (vslow.Pt() < fSlowPt)
		continue;

	      // D* mass cut
	      if (TMath::Abs(vds.M() - vd0.M() + MD_0 - MD_S) < fDeltaMDs) {
		if (fVerbose > 1)
		  cout << "==>HFB2D0Ds> D* candidate found! Mass: " << vds.M() << endl;
		nDs++;
		  
		// Fill TAnaCand
		std::vector<reco::Track> fitTracks;
		fitTracks.clear();
		fitTracks.push_back(slow);
		// fitTracks.push_back(D0 candidate);
		  
		//m1.SetXYZM(tMuon1.px(), tMuon1.py(), tMuon1.pz(), MMUON);
		//m2.SetXYZM(tMuon2.px(), tMuon2.py(), tMuon2.pz(), MMUON);
		  
		TAnaVertex DsVtx;
		//DsVtx = DoVertexFit(fitTracks);
		DsVtx.addTrack(islow);
		  
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
		ds_cand->fSig2 = ds_cand->fSig1 + 2;

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
		s_track->fIndex    = islow;

	      } // mass cut
	    } // charge cut D*
	  } // track c
	} // mass cut
      } // charge cut D0
    } // track b (D0 reco)
  } // track a (D0 reco)

  if (fVerbose > 0 && nD0 > 0)
    cout << "==>HFB2D0Ds> " << nD0 << " D0 candidates found." << endl;
  if (fVerbose > 0 && nDs > 0)
    cout << "==>HFB2D0Ds> " << nDs << " D* candidates found." << endl;

  return;
}


// ----------------------------------------------------------------------
TAnaVertex HFB2D0Ds::DoVertexFit(std::vector<reco::Track> &Tracks){
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
      if (fVerbose > 1) cout << "==>HFB2D0Ds> Something went wrong! SecVtx nan - Aborting... !" << endl;
      anaVt.fPoint.SetXYZ(0,0,0);
      return anaVt;
    }
  } else {
    if (fVerbose > 1) cout << "==>HFB2D0Ds> KVF failed! Aborting... !" << endl;
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
void  HFB2D0Ds::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFB2D0Ds::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFB2D0Ds);
