#include "HFLambdas.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <TrackingTools/TransientTrack/interface/TrackTransientTrack.h>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include<vector>

extern TAna01Event *gHFEvent;

using std::cout;
using std::endl;
using namespace edm;
using namespace reco;

class HFLambdaCut : public HFMaxDocaCut
{
public:
    HFLambdaCut(double maxPAngle = 3.2, HFNodeCut *lambdaCut = NULL, double docaCut = 1E7) :
        HFMaxDocaCut(docaCut), fLambdaCut(lambdaCut), fMaxPAngle(maxPAngle) {}
    virtual bool operator()();

protected:
    HFNodeCut *fLambdaCut; // cross reference for pointing angle calculation!
    double fMaxPAngle; // max pointing angle (in radians)
};

bool HFLambdaCut::operator()()
{
    bool result = HFMaxDocaCut::operator()();

    if(result && fLambdaCut)
    {

        // check the pointing angle
        const double pAngle = fLambdaCut->fPtCand.Angle(fLambdaCut->fVtxPos - this->fVtxPos);
        result = pAngle < fMaxPAngle;
    }

    return result;
}

class HFV0Cut : public HFNodeCut
{
public:
    HFV0Cut(double maxChi2, double minD3, double maxDoca) :
        fCutChi2(maxChi2), fCutD3(minD3), fCutMaxDoca(maxDoca) {};
    virtual bool operator()();

protected:
    double fCutChi2;
    double fCutD3;
    double fCutMaxDoca;
};

bool HFV0Cut::operator()()
{
    return ((fVtxChi2 < fCutChi2) && (fVtxPos.Mag() > fCutD3) && (fMaxDoca < fCutMaxDoca));
}

HFLambdas::HFLambdas(const ParameterSet& iConfig) :
    fVerbose(iConfig.getUntrackedParameter<int>("verbose",0)),
    fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))),
    fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
    fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
    fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
    fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 0.4)),
    fProtonPt(iConfig.getUntrackedParameter<double>("protonPt", 0.4)),
    fTrackNormChi2(iConfig.getUntrackedParameter<double>("trackNormChi2", 7)),
    fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 1)),
    fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow",0.3)),
    fksWindow(iConfig.getUntrackedParameter<double>("ksWindow",0.3)),
    fL0Window(iConfig.getUntrackedParameter<double>("L0Window",0.3)),
    fLbWindow(iConfig.getUntrackedParameter<double>("LbWindow",0.8)),
    fPsiEffWindow(iConfig.getUntrackedParameter<double>("psiEffWindow",0.12)),
    fL0EffWindow(iConfig.getUntrackedParameter<double>("L0EffWindow",0.015)),
    fEffMaxChi2(iConfig.getUntrackedParameter<double>("EffMaxChi2",3.84)),
    fEffMin3d(iConfig.getUntrackedParameter<double>("EffMin3d",.5)),
    fEffMaxDoca(iConfig.getUntrackedParameter<double>("EffMaxDoca",.05)),
    fDeltaR(iConfig.getUntrackedParameter<double>("deltaR",99.0)),
    fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
    fPAngle(iConfig.getUntrackedParameter<double>("pAngle",0.1)),
    fUseV0producer(iConfig.getUntrackedParameter<bool>("useV0",false)),
    fDoVcands(iConfig.getUntrackedParameter<bool>("doVcands",false))

{
    cout << "----------------------------------------------------------------------" << endl;
    cout << "--- HFLambdas constructor" << endl;
    cout << "---  verbose:					" << fVerbose << endl;
    cout << "---  tracksLabel:                              " << fTracksLabel << endl;
    cout << "---  PrimaryVertexLabel:                       " << fPrimaryVertexLabel << endl;
    cout << "---  muonsLabel:				" << fMuonsLabel << endl;
    cout << "---  muonPt:                                   " << fMuonPt << endl;
    cout << "---  pionPt:                                   " << fPionPt << endl;
    cout << "---  protonPt:                                 " << fProtonPt << endl;
    cout << "---  fTrackNormChi2:                           " << fTrackNormChi2 << endl;
    cout << "---  psiMuons:					" << fPsiMuons << endl;
    cout << "---  psiWindow:                                " << fPsiWindow << endl;
    cout << "---  ksWindow:                                 " << fksWindow << endl;
    cout << "---  L0Window:                                 " << fL0Window << endl;
    cout << "---  LbWindow:                                 " << fLbWindow << endl;
    cout << "---  PsiEffWindow:                             " << fL0EffWindow << endl;
    cout << "---  L0EffWindow:                              " << fPsiEffWindow << endl;
    cout << "---  EffMaxChi2:                               " << fEffMaxChi2 << endl;
    cout << "---  EffMin3d:                                 " << fEffMin3d << endl;
    cout << "---  EffMaxDoca:                               " << fEffMaxDoca << endl;
    cout << "---  deltaR:                                   " << fDeltaR << endl;
    cout << "---  maxDoca:				    " << fMaxDoca << endl;
    cout << "---  pAngle:                                   " << fPAngle << endl;
    cout << "---  fUseV0producer                            " << fUseV0producer << endl;
    cout << "---  fDoVcands:                                " << fDoVcands << endl;

} // HFLambdas()

HFLambdas::~HFLambdas()
{}

void HFLambdas::beginJob() {} // beginJob()
void HFLambdas::endJob() {} // endJob()

void HFLambdas::analyze(const Event& iEvent, const EventSetup& iSetup)
{
    // -- get the magnetic field
    ESHandle<MagneticField> magfield;
    iSetup.get<IdealMagneticFieldRecord>().get(magfield);
    const MagneticField *field = magfield.product();

    // -- get the primary vertex
    Handle<VertexCollection> recoPrimaryVertexCollection;
    iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
    if(!recoPrimaryVertexCollection.isValid())
    {
        cout << "==>HFLambdas> No primary vertex collection found, skipping" << endl;
        return;
    }

    const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
    if (vertices.size() == 0)
    {
        cout << "==>HFLambdas> No primary vertex found, skipping" << endl;
        return;
    }

    fPV = vertices[gHFEvent->fBestPV];
    if(fVerbose > 0)
        cout << "==>HFLambdas: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;

    // -- get the collection of muons
    Handle<MuonCollection> hMuons;
    iEvent.getByLabel(fMuonsLabel, hMuons);
    if(!hMuons.isValid())
    {
        cout << "==>HFLambdas> No valid MuonCollection with label "<< fMuonsLabel << " found, skipping" << endl;
        return;
    }

    // -- get the collection of tracks
    Handle<View<Track> > hTracks;
    iEvent.getByLabel(fTracksLabel, hTracks);
    if(!hTracks.isValid())
    {
        cout << "==>HFLambdas> No valid TrackCollection with label "<< fTracksLabel << " found, skipping" << endl;
        return;
    }

    // -- Transient tracks for vertexing
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
    if(!fTTB.isValid())
    {
        cout << "==>HFLambdas: Error: no TransientTrackBuilder found."<<endl;
        return;
    }

    // -- get the collection of muons and store their corresponding track indices
    vector<unsigned int> muonIndices;
    for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon)
    {
        int im = muon->track().index();
        if(im >= 0) muonIndices.push_back(im);
    }

    if (fVerbose > 0)
    {
        cout << "==>HFLambdas> nMuons = " << hMuons->size() << endl;
        cout << "==>HFLambdas> nMuonIndices = " << muonIndices.size() << endl;
    }

    if(muonIndices.size() < (unsigned)fPsiMuons) return;

    // -- Build muon list
    TLorentzVector tlv;
    vector<pair<int, reco::Track> > pilist,prlist;
    vector<pair<int, TLorentzVector> > mulist1,mulist2;
    int isMuon(0);
    if(2 == fPsiMuons)
    {
        mulist1.reserve(10);
        mulist2.reserve(10);
    }
    else
    {
        mulist1.reserve(100);
        mulist2.reserve(100);
    }
    if(!fUseV0producer)
    {
        pilist.reserve(100);
        prlist.reserve(100);
    }

    for(unsigned int itrack = 0; itrack < hTracks->size(); ++itrack)
    {
        TrackBaseRef rTrackView(hTracks,itrack);
        Track tTrack(*rTrackView);
        if(!fUseV0producer && (tTrack.chi2()/tTrack.ndof()) < fTrackNormChi2)
        {
            if(tTrack.pt() > fPionPt)   pilist.push_back(make_pair(itrack,tTrack));
            if(tTrack.pt() > fProtonPt) prlist.push_back(make_pair(itrack,tTrack));
        }

        tlv.SetXYZM(tTrack.px(),tTrack.py(),tTrack.pz(),MMUON);
        if(tTrack.pt() < fMuonPt) continue;
        if(2 == fPsiMuons)
        {
            for (unsigned int im = 0; im < muonIndices.size(); ++im)
            {
                if(muonIndices[im] == itrack)
                {
                    ++isMuon;
                    mulist1.push_back(make_pair(itrack,tlv));
                    mulist2.push_back(make_pair(itrack,tlv));
                }
            }
        }
        else if (1 == fPsiMuons)
        {
            for (unsigned int im = 0; im < muonIndices.size(); ++im)
            {
                if (muonIndices[im] == itrack)
                {
                    ++isMuon;
                    mulist1.push_back(make_pair(itrack,tlv));
                }
            }
            mulist2.push_back(make_pair(itrack,tlv));
        }
        else
        {
            mulist1.push_back(make_pair(itrack,tlv));
            mulist2.push_back(make_pair(itrack,tlv));
        }
    }

    if (isMuon < fPsiMuons)
    {
        if(fVerbose > 0) cout << "==>HFLambdas> Not enough muons found for J/Psi candidates combinatorics: isMuon = "
                                  << isMuon << ", required are " << fPsiMuons << endl;
        return;
    }

    if(fVerbose > 5) cout << "-------------------------------- Run: " << iEvent.id().run()
                              << " Event: " << iEvent.id().event() << " LS: " << iEvent.luminosityBlock() << endl;

    HFTwoParticleCombinatorics a(fVerbose, fTTB.product());
    vector<pair<int, int> > psiList;
    a.combine(psiList, mulist1, mulist2, MJPSI-fPsiWindow, MJPSI+fPsiWindow, 1);
    if(fVerbose > 0) cout << "==>HFLambdas> J/Psi list size: " << psiList.size() << endl;

    vector<HFTwoParticleState> l0List;
    if(!fUseV0producer)
    {
        a.combine(l0List, pilist, MPION, prlist, MPROTON, MLAMBDA_0-fL0Window, MLAMBDA_0+fL0Window, fMaxDoca, 0);
    }
    else
    {
        HFTwoParticleState tps;
        edm::InputTag lamCollectionTag("generalV0Candidates:Lambda");
        Handle<reco::VertexCompositeCandidateCollection> lambdaCollection;
        iEvent.getByLabel(lamCollectionTag, lambdaCollection);
        if(lambdaCollection->size()==0) return;
        for(VertexCompositeCandidateCollection::const_iterator iLam = lambdaCollection->begin(); iLam != lambdaCollection->end(); iLam++)
        {
            Track cand1 = *((RecoChargedCandidate *) (iLam->daughter(0)))->track() ;
            Track cand2 = *((RecoChargedCandidate *) (iLam->daughter(1)))->track() ;
            int iProton =-1;
            int iPion = -1;

            for(unsigned int itrack = 0; itrack < hTracks->size(); ++itrack)
            {
                TrackBaseRef rTrackView(hTracks,itrack);
                Track tTrack(*rTrackView);
                if(tTrack.pt()>fProtonPt &&
                        tTrack.referencePoint() == cand1.referencePoint()
                        && tTrack.momentum() == cand1.momentum()
                        && tTrack.charge() == cand1.charge()) iProton=itrack;
                if(tTrack.pt() > fPionPt &&
                        tTrack.referencePoint() == cand2.referencePoint()
                        && tTrack.momentum() == cand2.momentum()
                        && tTrack.charge() == cand2.charge()) iPion=itrack;
            }
            if(iProton>-1 && iPion > -1 && iProton!=iPion)
            {
                tps.id1=iPion;
                tps.id2=iProton;
                math::XYZTLorentzVector p4=iLam->p4();
                tps.mother.SetXYZT(p4.x(), p4.y(), p4.z(), p4.t());
                l0List.push_back(tps);
            }
            else
            {
                if(cand2.pt()> fPionPt && cand1.pt()>fProtonPt) cout <<"V0 not found!  "<<iProton<<"  "<<iPion<<endl;
            }
        }
    }

    if (fVerbose > 0) cout << "==>HFLambdas> Lambda_0 list size: " << l0List.size() << endl;

    vector<HFTwoParticleState > K0List;
    if(!fUseV0producer)
    {
        a.combine(K0List, pilist, MPION, pilist, MPION, MKSHORT-fksWindow, MKSHORT+fksWindow, fMaxDoca, 1);
        if (fVerbose > 0) cout << "==>HFLambdas> K0 list size: " << K0List.size() << endl;
    }

    const int V0Cand = fUseV0producer ? 10000 : 0;
    // do the vertex fitting...
    HFSequentialVertexFit aSeq(hTracks, fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);
    //TLorentzVector tlvPsi, tlvMu1, tlvMu2, tlvPion, tlvProton, tlvLambda0, tlvLambdaB;

    for (unsigned int i = 0; i < psiList.size(); i++)
    {
        unsigned int iMuon1 = psiList[i].first;
        unsigned int iMuon2 = psiList[i].second;

        TrackBaseRef mu1TrackView(hTracks, iMuon1);
        Track tMuon1(*mu1TrackView);
        TLorentzVector tlvMu1;
        tlvMu1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);

        TrackBaseRef mu2TrackView(hTracks, iMuon2);
        Track tMuon2(*mu2TrackView);
        TLorentzVector tlvMu2;
        tlvMu2.SetPtEtaPhiM(tMuon2.pt(),tMuon2.eta(),tMuon2.phi(), MMUON);

        const TLorentzVector tlvPsi = tlvMu1 + tlvMu2;
        if(fDoVcands) // this is mainly for efficiency studies to have the neutral cands separately from the full fit
            // here we have to apply more stringent cuts to have the combinatorics under control
        {
            HFDecayTree theTree(500443+V0Cand);
            // make a Jpsi
            if (TMath::Abs(tlvPsi.M() - MJPSI) <= fPsiEffWindow)
            {
                theTree.addTrack(iMuon1,13);
                theTree.addTrack(iMuon2,13);
                theTree.setNodeCut(RefCountedHFNodeCut(new HFV0Cut(fEffMaxChi2, 0, fEffMaxDoca)));
                aSeq.doFit(&theTree);
            }
        }

        // create Lambda_B candidates
        for(unsigned int j = 0; j < l0List.size(); j++)
        {

            unsigned int iPion = l0List[j].id1;
            unsigned int iProton = l0List[j].id2;

            if (iPion == iMuon1 || iPion == iMuon2) continue;
            if (iProton == iMuon1 || iProton == iMuon2) continue;

            const TLorentzVector tlvLambda0 = l0List[j].mother;

            const TLorentzVector tlvLambdaB = tlvPsi + tlvLambda0;
            if (TMath::Abs(tlvLambdaB.M() - MLAMBDA_B) > fLbWindow) continue;

            if(fVerbose > 5) cout << "Making HFDecayTree's with tracks (mu:mu:pr:pi): "
                                      << iMuon1 << ":" << iMuon2 << ":" << iPion << ":" << iProton << endl;

            // without J/Psi mass constraint
            HFDecayTree theTree(605122+V0Cand);
            HFDecayTreeIterator iterator;
            iterator = theTree.addDecayTree(600443+V0Cand,0); // no vertexing & no mass constraint
            iterator->addTrack(iMuon1,13);
            iterator->addTrack(iMuon2,13);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            iterator = theTree.addDecayTree(603122+V0Cand,1); // no mass constraint, with own vertex of Lamda0
            iterator->addTrack(iPion,211);
            iterator->addTrack(iProton,2212);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle, &(*iterator->getNodeCut()),fMaxDoca)));

            aSeq.doFit(&theTree);

            // without J/Psi mass constraint, but with an own J/Psi vertex
            theTree.clear();
            theTree.particleID = 705122+V0Cand;

            iterator = theTree.addDecayTree(700443+V0Cand,1); // vertexing but no mass constraint...
            iterator->addTrack(iMuon1,13);
            iterator->addTrack(iMuon2,13);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            iterator = theTree.addDecayTree(703122+V0Cand,1); // Lambda0 with vertexing
            iterator->addTrack(iPion,211);
            iterator->addTrack(iProton,2212);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle, &(*iterator->getNodeCut()), fMaxDoca)));

            aSeq.doFit(&theTree);

            // with J/Psi mass constraint
            theTree.clear();
            theTree.particleID = 805122+V0Cand;

            iterator = theTree.addDecayTree(800443+V0Cand,1,MJPSI); // J/Psi needs an own vertex!!
            iterator->addTrack(iMuon1,13);
            iterator->addTrack(iMuon2,13);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            iterator = theTree.addDecayTree(803122+V0Cand,1);
            iterator->addTrack(iPion,211);
            iterator->addTrack(iProton,2212);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
            theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle,&(*iterator->getNodeCut()),fMaxDoca)));

            aSeq.doFit(&theTree);

            // with mass constraint for J/Psi and Lambda0
            theTree.clear();
            theTree.particleID = 905122+V0Cand;

            iterator = theTree.addDecayTree(900443+V0Cand,1,MJPSI);
            iterator->addTrack(iMuon1,13);
            iterator->addTrack(iMuon2,13);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

            iterator = theTree.addDecayTree(903122+V0Cand,1,MLAMBDA_0);
            iterator->addTrack(iPion,211);
            iterator->addTrack(iProton,2212);
            iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
            theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle,&(*iterator->getNodeCut()),fMaxDoca)));

            aSeq.doFit(&theTree);
        }

        if(!fUseV0producer)
        {
            // create B0 candidates
            for(unsigned int j = 0; j < K0List.size(); j++)
            {

                unsigned int iPion1 = K0List[j].id1;
                unsigned int iPion2 = K0List[j].id2;

                if (iPion1 == iMuon1 || iPion1 == iMuon2) continue;
                if (iPion2 == iMuon1 || iPion2 == iMuon2) continue;

                const TLorentzVector tlvK0 = K0List[j].mother;

                const TLorentzVector tlvB0 = tlvPsi + tlvK0;
                if (TMath::Abs(tlvB0.M() - MB_0) > fksWindow) continue;

                // without J/Psi mass constraint
                HFDecayTree theTree(600511);

                HFDecayTreeIterator iterator = theTree.addDecayTree(600443,0); // no vertexing & no mass constraint
                iterator->addTrack(iMuon1,13);
                iterator->addTrack(iMuon2,13);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                iterator = theTree.addDecayTree(600310,1); // no mass constraint, with own vertex of Kshort
                iterator->addTrack(iPion1,211);
                iterator->addTrack(iPion2,211);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle, &(*iterator->getNodeCut()),fMaxDoca)));

                aSeq.doFit(&theTree);

                // without J/Psi mass constraint, but with an own J/Psi vertex
                theTree.clear();
                theTree.particleID = 700511;

                iterator = theTree.addDecayTree(700443,1); // vertexing but no mass constraint...
                iterator->addTrack(iMuon1,13);
                iterator->addTrack(iMuon2,13);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                iterator = theTree.addDecayTree(700310,1); // Lambda0 with vertexing
                iterator->addTrack(iPion1,211);
                iterator->addTrack(iPion2,211);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle, &(*iterator->getNodeCut()), fMaxDoca)));

                aSeq.doFit(&theTree);

                // with J/Psi mass constraint
                theTree.clear();
                theTree.particleID = 800511;

                iterator = theTree.addDecayTree(800443,1,MJPSI); // J/Psi needs an own vertex!!
                iterator->addTrack(iMuon1,13);
                iterator->addTrack(iMuon2,13);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                iterator = theTree.addDecayTree(800310,1);
                iterator->addTrack(iPion1,211);
                iterator->addTrack(iPion2,211);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
                theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle,&(*iterator->getNodeCut()),fMaxDoca)));

                aSeq.doFit(&theTree);

                // with mass constraint for J/Psi and Lambda0
                theTree.clear();
                theTree.particleID = 900511;

                iterator = theTree.addDecayTree(900443,1,MJPSI);
                iterator->addTrack(iMuon1,13);
                iterator->addTrack(iMuon2,13);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));

                iterator = theTree.addDecayTree(900310,1,MB_0);
                iterator->addTrack(iPion1,211);
                iterator->addTrack(iPion2,211);
                iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
                theTree.setNodeCut(RefCountedHFNodeCut(new HFLambdaCut(fPAngle,&(*iterator->getNodeCut()),fMaxDoca)));

                aSeq.doFit(&theTree);
            }      // K0List
        }         // fUseV0producer
    }               // psiList

    if(fDoVcands) // this is mainly for efficiency studies to have the neutral cands separately from the full fit
        // here we have to apply more stringent cuts to have the combinatorics under control
    {
        for(unsigned int j = 0; j < l0List.size(); j++)
        {
            const unsigned int iPion = l0List[j].id1;
            const unsigned int iProton = l0List[j].id2;
            const TLorentzVector tlvLambda0 = l0List[j].mother;
            // make a Lambda
            if (TMath::Abs(tlvLambda0.M() - MLAMBDA_0) <= fL0EffWindow)
            {
                TrackBaseRef tbrPi(hTracks, iPion);
                Track tPi(*tbrPi);
		if (tPi.chi2()/tPi.ndof() > 5) continue;
                TrackBaseRef tbrPr(hTracks, iProton);
                Track tPr(*tbrPr);
		if (tPr.chi2()/tPr.ndof() > 5) continue;
                if(tPr.pt() > tPi.pt())
                {
                    HFDecayTree theTree(503122+V0Cand);
                    theTree.addTrack(iPion,211);
                    theTree.addTrack(iProton,2212);
                    theTree.setNodeCut(RefCountedHFNodeCut(new HFV0Cut(fEffMaxChi2, fEffMin3d, fEffMaxDoca)));
                    aSeq.doFit(&theTree);
                }
            }
        }
    }
} // analyze()

// define this as a plug-in
DEFINE_FWK_MODULE(HFLambdas);
