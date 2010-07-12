#include "massReader.hh"
#include <cstdlib>

using namespace std;

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),reduced_tree(NULL)
{
	fMomentumPtr = &fMomentum;
	fPVPositionPtr = &fPVPosition;
	fCandVertexPtr = &fCandVertex;
	
	fTreeName = "massReader reduced tree";
} // massReader()

massReader::~massReader()
{ } // ~massReader()

void massReader::eventProcessing()
{
	int j,nc;
	
	// Fill a reduced tree
	nc = fpEvt->nCands();
	for (j = 0; j<nc; j++) {
		
		if(loadCandidateVariables(fpEvt->getCand(j)))
			reduced_tree->Fill();
	}
} // massReader::eventProcessing()

int massReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaCand *momCand;
	TVector3 v1,v2;
	
	// Save in the tree
	fCandidate = pCand->fType;
	fMomentum = pCand->fPlab;
	fPVPosition = fpEvt->bestPV()->fPoint;
	fCandVertex = pCand->fVtx.fPoint;
	fPt = pCand->fPlab.Perp();
	fMass = pCand->fMass;
	fNbrMuons = countMuons(pCand);
	fD3 = pCand->fVtx.fD3d;
	fD3E = pCand->fVtx.fD3dE;
	fDxy = pCand->fVtx.fDxy;
	fDxyE = pCand->fVtx.fDxyE;
	fChi2 = pCand->fVtx.fChi2;
	fNdof = pCand->fVtx.fNdof;
	fMaxDoca = pCand->fMaxDoca;
	fIso = calculateIsolation(pCand,1.0);
	
	if (pCand->fMom >= 0) {
		momCand = fpEvt->getCand(pCand->fMom);
		v1 = pCand->fPlab;
		v2 = pCand->fVtx.fPoint - momCand->fVtx.fPoint;
	} else {
		v1 = pCand->fPlab;
		v2 = pCand->fVtx.fPoint - fpEvt->bestPV()->fPoint;
	}
	
	fAlpha = v1.Angle(v2);
	// project to xy plane
	v1.SetZ(0);
	v2.SetZ(0);
	fAlphaXY = v1.Angle(v2);
	
	// do this at the end so the checkTruth algorithm can use
	// all variables of this candidate.
	fTruth = checkTruth(pCand);
	
	return 1;
} // loadCandidateVariables()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T",fTreeName);
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("p","TVector3",&fMomentumPtr);
	reduced_tree->Branch("pv_position","TVector3",&fPVPositionPtr);
	reduced_tree->Branch("vertex_position","TVector3",&fCandVertexPtr);
	reduced_tree->Branch("pt",&fPt,"pt/F");
	reduced_tree->Branch("mass",&fMass,"mass/F");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
	reduced_tree->Branch("ident_muons",&fNbrMuons,"ident_muons/F");
	reduced_tree->Branch("d3",&fD3,"d3/F");
	reduced_tree->Branch("d3e",&fD3E,"d3e/F");
	reduced_tree->Branch("dxy",&fDxy,"dxy/F");
	reduced_tree->Branch("dxye",&fDxyE,"dxye/F");
	reduced_tree->Branch("alpha",&fAlpha,"alpha/F");
	reduced_tree->Branch("alpha_xy",&fAlphaXY,"alpha_xy/F");
	reduced_tree->Branch("chi2",&fChi2,"chi2/F");
	reduced_tree->Branch("Ndof",&fNdof,"Ndof/F");
	reduced_tree->Branch("max_doca",&fMaxDoca,"max_doca/F");
	reduced_tree->Branch("iso",&fIso,"iso/F");
} // massReader::bookHist()

int massReader::checkTruth(TAnaCand *cand)
{
	int result = 0;
	TAnaTrack *pTrack;
	TGenCand *truthParticle;
	TGenCand *trackParticle;
	int nGens, truth_type, j;
	
	truth_type = cand->fType % 1000;
	nGens = fpEvt->nGenCands();
	
	pTrack = fpEvt->getSigTrack(cand->fSig1);
	pTrack = fpEvt->getRecTrack(pTrack->fIndex);
	
	if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= nGens) goto bail;
	truthParticle = fpEvt->getGenCand(pTrack->fGenIndex);
	
	while (abs(truthParticle->fID) != truth_type) {
		if (truthParticle->fMom1 < 0 || truthParticle->fMom1 >= nGens) goto bail;
		truthParticle = fpEvt->getGenCand(truthParticle->fMom1);
	}
	
	// check to see if the other tracks originate from the same track...
	for (j = cand->fSig1+1; j <= cand->fSig2; j++) {
		
		pTrack = fpEvt->getSigTrack(j);
		pTrack = fpEvt->getRecTrack(pTrack->fIndex);
		if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= nGens) goto bail;
		
		trackParticle = fpEvt->getGenCand(pTrack->fGenIndex);
		while (trackParticle->fNumber != truthParticle->fNumber) {
			if (trackParticle->fMom1 < 0 || trackParticle->fMom1 >= nGens) goto bail;
			trackParticle = fpEvt->getGenCand(trackParticle->fMom1);
		}
	}
	
	result = 1;
bail:
	return result;
} // checkTruth()

// count the number of identified muons in the candidate 'cand'
int massReader::countMuons(TAnaCand *cand)
{
	TAnaTrack *track;
	unsigned nbrMu = 0;
	int j;
	
	// check for invalid candidates
	if (cand->fSig1 < 0)
		return 0;
	
	for (j = cand->fSig1; j <= cand->fSig2; j++) {
		
		track = fpEvt->getSigTrack(j); // signal track
		// check if this is supposed to be a muon!
		if (abs(track->fMCID) == 13) {
			track = fpEvt->getRecTrack(track->fIndex); // rec track
			if (track->fMuID > 0 && (track->fMuID & 6))
				nbrMu++;
		}
	}
	
	return nbrMu;
} // countMuons()

void massReader::buildDecay(TGenCand *gen, multiset<int> *particles)
{
	particles->insert(abs(gen->fID));
	for (int j = gen->fDau1; j <= gen->fDau2 && j>= 0; j++)
		buildDecay(fpEvt->getGenCand(j),particles);
} // buildDecay()

float massReader::calculateIsolation(TAnaCand *pCand, double openingAngle)
{
	double iso; // calculate in double precision and return single
	TVector3 plabB;
	int j,ntracks;
	TAnaTrack *pTrack;
	double r;
	set<int> usedTracks;
	
	plabB = pCand->fPlab;
	
	// build list of tracks in the plabB
	for (j = pCand->fSig1; j <= pCand->fSig2; j++) {
		pTrack = fpEvt->getSigTrack(j);
		usedTracks.insert(pTrack->fIndex);
	}
	
	iso = 0.0;
	ntracks = fpEvt->nRecTracks();
	for (j = 0; j < ntracks; j++) {
		
		if (usedTracks.find(j) != usedTracks.end()) continue; // this tracks belongs to the candidate
		
		pTrack = fpEvt->getRecTrack(j);
		r = plabB.DeltaR(pTrack->fPlab);
		if (r < openingAngle)
			iso += pTrack->fPlab.Perp();
	}
	
	iso = plabB.Perp() / (plabB.Perp() + iso);
	
	return (float)iso;
} // calculateIsolation()
