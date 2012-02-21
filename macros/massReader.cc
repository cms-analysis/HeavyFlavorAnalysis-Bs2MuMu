#include "massReader.hh"
#include <cstdlib>
#include <cmath>

using namespace std;

trigger_table_t g_trigger_table_signal [] = {
	{kHLT_DoubleMu3_Bs_Bit,	"HLT_DoubleMu3_Bs", std::pair<int64_t,int64_t>(160329ll,161176ll)},
	{kHLT_DoubleMu2_Bs_Bit, "HLT_DoubleMu2_Bs", std::pair<int64_t,int64_t>(161216ll,167913ll)},
	{kHLT_Dimuon4_Bs_Barrel_Bit, "HLT_Dimuon4_Bs_Barrel", std::pair<int64_t,int64_t>(170249ll,173198ll)}, // barrel trigger
	{kHLT_Dimuon6_Bs_Bit, "HLT_Dimuon6_Bs", std::pair<int64_t,int64_t>(170249ll,173198ll)}, // endcap trigger
	{kHLT_DoubleMu4_Dimuon4_Bs_Barrel_Bit, "HLT_DoubleMu4_Dimuon4_Bs_Barrel", std::pair<int64_t,int64_t>(173236ll,180252ll)}, // barrel trigger
	{kHLT_DoubleMu4_Dimuon6_Bs_Bit, "HLT_DoubleMu4_Dimuon6_Bs", std::pair<int64_t,int64_t>(173236ll,180252ll)} // endcap trigger
};

trigger_table_t g_trigger_table_norm [] = {
	{kHLT_DoubleMu3_Jpsi_Bit, "HLT_DoubleMu3_Jpsi", std::pair<int64_t,int64_t>(160329ll,163261ll)},
	{kHLT_Dimuon6p5_Jpsi_Displaced_Bit,"HLT_Dimuon6p5_Jpsi_Displaced", std::pair<int64_t,int64_t>(163269ll,163869ll)},
	{kHLT_Dimuon7_Jpsi_Displaced_Bit,"HLT_Dimuon7_Jpsi_Displaced", std::pair<int64_t,int64_t>(165088ll,167913ll)},
	{kHLT_DoubleMu3p5_Jpsi_Displaced_Bit,"HLT_DoubleMu3p5_Jpsi_Displaced", std::pair<int64_t,int64_t>(170249ll,173198ll)},
	{kHLT_DoubleMu4_Jpsi_Displaced_Bit,"HLT_DoubleMu4_Jpsi_Displaced", std::pair<int64_t,int64_t>(173236ll,180252ll)}
};

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),
	reduced_tree(NULL),
	fTruthType(0),
	fCutFileParsed(false),
	fCutTriggered(false),
	fCutCand(0),
	fCutFlight3dSign(0.0),
	fCutChi2(0.0),
	fCutPt(0.0),
	fCutAlpha(0.0),
	fCutChi2ByNdof(0.0),
	fCutTruth(0),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutMuID_mask(0),
	fCutMuID_reqall(false),
	fCutOppSign_mu(false)
{	
	fTreeName = "massReader reduced tree";
	
	// add only the ones we actually need.
	// feel free to enlarge this set.
	stableParticles.insert(13); // muon
	stableParticles.insert(321); // kaon
} // massReader()

massReader::~massReader()
{ } // ~massReader()

void massReader::eventProcessing()
{
	int j,nc;
	
	// Fill Generator Block candidates (for efficiency determination)
	nc = fpEvt->nGenCands();
	for (j = 0; j < nc; j++) {
		if (loadGeneratorVariables(fpEvt->getGenCand(j)))
			reduced_tree->Fill();
	}
	
	// Fill a reduced tree
	nc = fpEvt->nCands();
	for (j = 0; j<nc; j++) {
		
		if(loadCandidateVariables(fpEvt->getCand(j)) && (!fCutFileParsed || applyCut()))
			reduced_tree->Fill();
	}
} // massReader::eventProcessing()

void massReader::clearVariables()
{
	fCandidate = 0;
	fMass = 0.0;
	fMassConstraint = 0.0;
	fTruth = 0;
	fTruthFlags = 0;
	fEffFlags = 0;
	fPt = 0.0;
	fNbrMuons = 0;
	fD3 = 0.0;
	fD3E = 0.0;
	fDxy = 0.0;
	fDxyE = 0.0;
	fAlpha = 0.0;
	fAlphaXY = 0.0;
	fChi2 = 0.0;
	fNdof = 0.0;
	fMaxDoca = 0.0;
	fIsoMoriond12 = 0.0;
	fDoca0 = FLT_MAX;
	fNbrNearby = 0;
	fTriggers = 0;
	fTriggersError = 0;
	fTriggersFound = 0;
	fTriggeredJPsi = 0;
	fTriggeredBs = 0;
	fCtau = 0.0;
	fCtauE = 0.0;
	fEta = 0.0;
	fPMu1_Gen = 0.0;
	fPMu2_Gen = 0.0;
	fPtMu1_Gen = 0.0;
	fPtMu2_Gen = 0.0;
	fEtaMu1_Gen = 0.0; // eta of gen muon 1
	fEtaMu2_Gen = 0.0; // eta of gen muon 2
	fPtMu1 = fPtMu2 = 0.0f;
	fMuID1 = fMuID2 = 0;
	fMuTight1 = fMuTight2 = 0;
	fEtaMu1 = fEtaMu2 = 0.0f;
	fTrackQual_mu1 = fTrackQual_mu2 = 0;
	fQ_mu1 = fQ_mu2 = 0;
	fDeltaR = 0.0f; // deltaR
	fNbrPV = 0;
	fDeltaPhiMu = 0;
	fIPCand = 0.0f;
	fIPCandE = 0.0f;
	
	memset(fTracksIx,0,sizeof(fTracksIx));
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksIPE,0,sizeof(fTracksIPE));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
} // clearVariables()

int massReader::loadGeneratorVariables(TGenCand *pGen)
{
	int save = 0;
	multiset<int> particles;
	map<int,int> genStructure;
	TGenCand *dau;
	bool firstMu;
	
	if (abs(pGen->fID) != fTruthType)
		goto bail;
	
	if (trueDecay.size() == 0)
		goto bail;
	
	buildDecay(pGen,&particles);
	particles.erase(22); // remove Bremsstrahlung
	if (particles != trueDecay)
		goto bail;
	
	// reset variables
	clearVariables();
	
	// this is the right decay...
	fCandidate = 0; // use 0 as special case...
	fMass = pGen->fP.M();
	fEffFlags = loadEfficiencyFlags(pGen);
	fPt = pGen->fP.Perp();
	fEta = pGen->fP.Eta();
	fTriggers = loadTrigger(&fTriggersError,&fTriggersFound);
	fTriggeredJPsi = hasTriggeredNorm();
	fTriggeredBs = hasTriggeredSignal();
	fNbrPV = fpEvt->nPV();
	
	// save the muon pt
	firstMu = true;
	findGenStructure(pGen,&genStructure);
	for (map<int,int>::const_iterator it = genStructure.begin(); it != genStructure.end(); ++it) {
		dau = fpEvt->getGenCand(it->first);
		switch (abs(dau->fID)) {
			case 13: // muon
				if (firstMu) {
					fPtMu1_Gen = dau->fP.Perp();
					fEtaMu1_Gen = dau->fP.Eta();
					fPMu1_Gen = dau->fP.P();
					if (it->second >= 0) {
						fPtMu1 = fpEvt->getRecTrack(it->second)->fPlab.Perp();
						fEtaMu1 = fpEvt->getRecTrack(it->second)->fPlab.Eta();
					}
				} else {
					fPtMu2_Gen = dau->fP.Perp();
					fEtaMu2_Gen = dau->fP.Eta();
					fPMu2_Gen = dau->fP.P();
					if (it->second >= 0) {
						fPtMu2 = fpEvt->getRecTrack(it->second)->fPlab.Perp();
						fEtaMu2 = fpEvt->getRecTrack(it->second)->fPlab.Eta();
					}
				}
				firstMu = false;
				break;
			default:
				break;
		}
	}
	if (fPtMu2_Gen > fPtMu1_Gen) {
		swap(fPtMu1,fPtMu2);
		swap(fPtMu1_Gen,fPtMu2_Gen);
		swap(fEtaMu1,fEtaMu2);
		swap(fEtaMu1_Gen,fEtaMu2_Gen);
		swap(fPMu1_Gen,fPMu2_Gen);
	}
	
	// save the candidate...
	save = 1;
bail:
	return save;
} // loadGeneratorVariables()

int massReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaTrack *pTrack;
	TAnaTrack *sigTrack;
	TAnaTrack *recTrack;
	TAnaCand *momCand;
	TGenCand *muGen;
	TVector3 v1,v2,uVector;
	unsigned j,k;
	bool firstMu = true;
	map<int,int> aTracks;
	map<int,int> cand_tracks;
	TVector3 plabMu1,plabMu2;
	
	clearVariables();
	
	// Save in the tree
	fCandidate = pCand->fType;
	fPt = pCand->fPlab.Perp();
	fMass = pCand->fMass;
	fMassConstraint = -1.0f;
	fNbrMuons = countMuons(pCand);
	fD3 = pCand->fVtx.fD3d;
	fD3E = pCand->fVtx.fD3dE;
	fDxy = pCand->fVtx.fDxy;
	fDxyE = pCand->fVtx.fDxyE;
	fChi2 = pCand->fVtx.fChi2;
	fNdof = pCand->fVtx.fNdof;
	fMaxDoca = pCand->fMaxDoca;
	fNbrPV = fpEvt->nPV();
	fIPCand = pCand->fPvLip;
	fIPCandE = pCand->fPvLipE;	
	fIsoMoriond12 = calculateIsolation(pCand);
	fDoca0 = calculateDoca0(pCand);
	fNbrNearby = countTracksNearby(pCand);
	
	fCtau = pCand->fTau3d;
	fCtauE = pCand->fTau3dE;
	fEta = pCand->fPlab.Eta();
	
	// Load muon variables
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it!=cand_tracks.end(); ++it) {
		
		sigTrack = fpEvt->getSigTrack(it->second);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		
		switch (abs(sigTrack->fMCID)) {
			case 13: // muon
				if (firstMu) {
					plabMu1 = sigTrack->fPlab;
					fPtMu1 = sigTrack->fPlab.Perp();
					fMuID1 = recTrack->fMuID > 0 ? recTrack->fMuID : 0;
					fTrackQual_mu1 = recTrack->fTrackQuality;
					fMuTight1 = isMuonTight(sigTrack);
					fEtaMu1 = sigTrack->fPlab.Eta();
					fQ_mu1 = recTrack->fQ;
				} else {
					plabMu2 = sigTrack->fPlab;
					fPtMu2 = sigTrack->fPlab.Perp();
					fMuID2 = recTrack->fMuID > 0 ? recTrack->fMuID : 0;
					fTrackQual_mu2 = recTrack->fTrackQuality;
					fMuTight2 = isMuonTight(sigTrack);
					fEtaMu2 = sigTrack->fPlab.Eta();
					fQ_mu2 = recTrack->fQ;
				}
				firstMu = false;
				break;
			default:
				break;
		}
	}
	
	// make sure mu1 is leading and mu2 subleading
	if (fPtMu1 < fPtMu2) {
		swap(plabMu1,plabMu2);
		swap(fPtMu1,fPtMu2);
		swap(fMuID1,fMuID2);
		swap(fMuTight1,fMuTight2);
		swap(fEtaMu1,fEtaMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
	}
	
	fDeltaPhiMu = plabMu1.DeltaPhi(plabMu2);
	
	// Clean entries of nearest tracks
	for (j = 0; j < NBR_TRACKS_STORE; j++) fTracksIx[j] = -1;
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksIPE,0,sizeof(fTracksIPE));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
	
	// fill the histogramms...
	for (j = 0, k = 0; j < pCand->fNstTracks.size(); j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		uVector = pCand->fPlab.Unit();
		
		if (k < NBR_TRACKS_STORE) {
			fTracksIx[k] = pCand->fNstTracks[j].first;
			fTracksIP[k] = pCand->fNstTracks[j].second.first;
			fTracksIPE[k] = pCand->fNstTracks[j].second.second;
			fTracksPT[k] = pTrack->fPlab.Perp();
			fTracksPTRel[k] = (pTrack->fPlab - (pTrack->fPlab * uVector) * uVector).Mag();
			k++;
		}
	}
	
	if (pCand->fMom >= 0) {
		momCand = fpEvt->getCand(pCand->fMom);
		v1 = pCand->fPlab;
		v2 = pCand->fVtx.fPoint - momCand->fVtx.fPoint;
	} else {
		v1 = pCand->fPlab;
		if (0 <= pCand->fPvIdx && pCand->fPvIdx < fpEvt->nPV())
			v2 = pCand->fVtx.fPoint - fpEvt->getPV(pCand->fPvIdx)->fPoint;
		else
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
	fTruthFlags = loadTruthFlags(pCand);
	fTriggers = loadTrigger(&fTriggersError,&fTriggersFound);
	fTriggeredJPsi = hasTriggeredNorm();
	fTriggeredBs = hasTriggeredSignal();
	fEffFlags = 0;
	
	// generator pt of muons...
	firstMu = true;
	findAllTrackIndices(pCand,&aTracks);
	for (map<int,int>::const_iterator it = aTracks.begin(); it != aTracks.end(); ++it) {
		pTrack = fpEvt->getRecTrack(it->first);
		if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= fpEvt->nGenCands())
			continue;
		
		muGen = fpEvt->getGenCand(pTrack->fGenIndex);
		if (abs(muGen->fID) == 13) {
			if (firstMu) {
				fPtMu1_Gen = muGen->fP.Perp();
				fEtaMu1_Gen = muGen->fP.Eta();
				fPMu1_Gen = muGen->fP.P();
			}
			else {
				fPtMu2_Gen = muGen->fP.Perp();
				fEtaMu2_Gen = muGen->fP.Eta();
				fPMu2_Gen = muGen->fP.P();
			}
			firstMu = false;
		}
	}
	
	if (fPtMu1_Gen < fPtMu2_Gen) {
		swap(fPtMu1_Gen,fPtMu2_Gen);
		swap(fEtaMu1_Gen,fEtaMu2_Gen);
		swap(fPMu1_Gen,fPMu2_Gen);
	}
	
	if (plabMu2.Perp() > 0)
		fDeltaR = plabMu1.DeltaR(plabMu2);
	
	return 1;
} // loadCandidateVariables()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T",fTreeName);
	
	// and add the branches
	reduced_tree->Branch("run",&fRun,"run/I");
	reduced_tree->Branch("event",&fEvt,"event/I");
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("pt",&fPt,"pt/F");
	reduced_tree->Branch("mass",&fMass,"mass/F");
	reduced_tree->Branch("mass_c",&fMassConstraint,"mass_c/F");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
	reduced_tree->Branch("truth_flags",&fTruthFlags,"truth_flags/I");
	reduced_tree->Branch("eff_flags",&fEffFlags,"eff_flags/I");
	reduced_tree->Branch("ident_muons",&fNbrMuons,"ident_muons/F");
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("id_mu1",&fMuID1,"id_mu1/I");
	reduced_tree->Branch("id_mu2",&fMuID2,"id_mu2/I");
	reduced_tree->Branch("tight_mu1",&fMuTight1,"tight_mu1/I");
	reduced_tree->Branch("tight_mu2",&fMuTight2,"tight_mu2/I");
	reduced_tree->Branch("eta_mu1",&fEtaMu1,"eta_mu1/F");
	reduced_tree->Branch("eta_mu2",&fEtaMu2,"eta_mu2/F");
	reduced_tree->Branch("track_qual_mu1",&fTrackQual_mu1,"track_qual_mu1/I");
	reduced_tree->Branch("track_qual_mu2",&fTrackQual_mu2,"track_qual_mu2/I");
	reduced_tree->Branch("q_mu1",&fQ_mu1,"q_mu1/I");
	reduced_tree->Branch("q_mu2",&fQ_mu2,"q_mu2/I");
	reduced_tree->Branch("delta_phi",&fDeltaPhiMu,"delta_phi/F");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("d3",&fD3,"d3/F");
	reduced_tree->Branch("d3e",&fD3E,"d3e/F");
	reduced_tree->Branch("dxy",&fDxy,"dxy/F");
	reduced_tree->Branch("dxye",&fDxyE,"dxye/F");
	reduced_tree->Branch("alpha",&fAlpha,"alpha/F");
	reduced_tree->Branch("alpha_xy",&fAlphaXY,"alpha_xy/F");
	reduced_tree->Branch("chi2",&fChi2,"chi2/F");
	reduced_tree->Branch("Ndof",&fNdof,"Ndof/F");
	reduced_tree->Branch("max_doca",&fMaxDoca,"max_doca/F");
	reduced_tree->Branch("pt_mu1_gen",&fPtMu1_Gen,"pt_mu1_gen/F");
	reduced_tree->Branch("pt_mu2_gen",&fPtMu2_Gen,"pt_mu2_gen/F");
	reduced_tree->Branch("eta_mu1_gen",&fEtaMu1_Gen,"eta_mu1_gen/F");
	reduced_tree->Branch("eta_mu2_gen",&fEtaMu2_Gen,"eta_mu2_gen/F");
	reduced_tree->Branch("p_mu1_gen",&fPMu1_Gen,"p_mu1_gen/F");
	reduced_tree->Branch("p_mu2_gen",&fPMu2_Gen,"p_mu2_gen/F");
	reduced_tree->Branch("iso_mor12",&fIsoMoriond12,"iso_mor12/F");
	reduced_tree->Branch("doca0",&fDoca0,"doca0/F");
	reduced_tree->Branch("ntrk",&fNbrNearby,"ntrk/I");
	reduced_tree->Branch("triggers",&fTriggers,"triggers/I");
	reduced_tree->Branch("triggers_error",&fTriggersError,"triggers_error/I");
	reduced_tree->Branch("triggers_found",&fTriggersFound,"triggers_found/I");
	reduced_tree->Branch("triggered_jpsi",&fTriggeredJPsi,"triggered_jpsi/I");
	reduced_tree->Branch("triggered_bs",&fTriggeredBs,"triggered_bs/I");
	reduced_tree->Branch("ctau",&fCtau,"ctau/F");
	reduced_tree->Branch("ctaue",&fCtauE,"ctaue/F");
	reduced_tree->Branch("eta",&fEta,"eta/F");
	reduced_tree->Branch("nbr_pv",&fNbrPV,"nbr_pv/I");
	reduced_tree->Branch("ip",&fIPCand,"ip/F");
	reduced_tree->Branch("ipe",&fIPCandE,"ipe/F");
	reduced_tree->Branch("tracks_ix",fTracksIx,Form("tracks_ix[%d]/I",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ip",fTracksIP,Form("tracks_ip[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ipe",fTracksIPE,Form("tracks_ipe[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_pt",fTracksPT,Form("tracks_pt[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ptrel",fTracksPTRel,Form("tracks_ptrel[%d]/F",NBR_TRACKS_STORE));
} // massReader::bookHist()

void massReader::closeHistFile()
{
	fpHistFile = reduced_tree->GetCurrentFile();
	treeReader01::closeHistFile();
} // massReader::closeHistFile()

int massReader::checkTruth(TAnaCand *cand)
{
	int result = 0;
	TAnaTrack *pTrack;
	TGenCand *truthParticle;
	TGenCand *trackParticle;
	int theTruthType;
	int nGens, j;
	
	if (fTruthType) theTruthType = fTruthType;
	else theTruthType = cand->fType % 10000;
	
	nGens = fpEvt->nGenCands();
	
	if (cand->fSig1 < 0 || cand->fSig2 < cand->fSig1) goto bail;
	pTrack = fpEvt->getSigTrack(cand->fSig1);
	pTrack = fpEvt->getRecTrack(pTrack->fIndex);
	
	if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= nGens) goto bail;
	truthParticle = fpEvt->getGenCand(pTrack->fGenIndex);
	
	while (abs(truthParticle->fID) != theTruthType) {
		if (truthParticle->fMom1 < 0 || truthParticle->fMom1 >= nGens) goto bail;
		truthParticle = fpEvt->getGenCand(truthParticle->fMom1);
	}
	
	// check to see if the other tracks originate from the same particle...
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

int massReader::loadTruthFlags(TAnaCand *cand)
{
	TAnaTrack *pTrack;
	TGenCand *truthParticle = NULL;
	TGenCand *pGen;
	int nGens;
	int result = 0;
	set<int> bhadrons;
	
	map<int,int> tracks_indices; // map(recTrackIx,sigTrackIndex)
	map<int,int>::iterator it;
	
	// build the number of b hadrons
	bhadrons.insert(511); // B0
	bhadrons.insert(521); // B+
	bhadrons.insert(531); // Bs
	bhadrons.insert(5122); // Lambda_b
	
	// get the list of tracks of the candidate
	findAllTrackIndices(cand,&tracks_indices);
	
	// iterate through all tracks and see if they have the same origin
	nGens = fpEvt->nGenCands();
	for (it = tracks_indices.begin(); it!=tracks_indices.end(); ++it) {
		pTrack = fpEvt->getRecTrack(it->first);
		if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= nGens) // no generator info available
			goto bail;
		
		// get the originating b-hadron (if available)
		pGen = fpEvt->getGenCand(pTrack->fGenIndex);
		while (bhadrons.count(abs(pGen->fID)) == 0) {
			if (0 <= pGen->fMom1 && pGen->fMom1 < nGens)
				pGen = fpEvt->getGenCand(pGen->fMom1);
			else
				goto bail; // not found
		}
		
		// found an originating b-hadron
		if (!truthParticle) truthParticle = pGen;
		else if (truthParticle->fNumber != pGen->fNumber) goto bail;
	}
	// still here?
	result |= kTruthSameB_Bit;
	
bail:
	return result;
} // loadTruthFlags()

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
			if (track->fMuID > 0 && (track->fMuID & 6)) // needs to be a tracker or global muon
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

TAnaCand* massReader::findCandidate(int candID, map<int,int> *particles)
{
	TAnaCand *result = NULL;
	TAnaCand *pCand;
	map<int,int> candParticles;
	map<int,int>::iterator it;
	int j;
	
	for (j = 0; j < fpEvt->nCands(); j++) {
		pCand = fpEvt->getCand(j);
		
		if (pCand->fType != candID) continue;
		
		// check if the tracks are compatible...
		candParticles.clear();
		findCandStructure(pCand,&candParticles);
		
		// now, compare the particles
		if (candParticles == (*particles)) {
			result = pCand;
			break;
		}
	}
	
	return result;
} // findCandidate()

void massReader::findCandStructure(TAnaCand* pCand, map<int,int> *particles)
{
	map<int,int>::iterator it;
	TAnaTrack *sigTrack;
	
	particles->clear();
	findAllTrackIndices(pCand,particles);
	
	for (it = particles->begin(); it != particles->end(); ++it) {
		sigTrack = fpEvt->getSigTrack(it->second);
		it->second = sigTrack->fMCID;
	}
} // findCandStructure()

void massReader::findGenStructure(TGenCand *pGen, map<int,int> *particles)
{
	map<int,int>::const_iterator it;
	int j;
	
	// insert this
	for (j = 0; j < fpEvt->nRecTracks(); j++) {
		if (fpEvt->getRecTrack(j)->fGenIndex == pGen->fNumber)
			break;
	}
	particles->insert(make_pair<int,int>(pGen->fNumber, j < fpEvt->nRecTracks() ? j : -1));
	
	
	for (j = pGen->fDau1; j <= pGen->fDau2; j++)
		findGenStructure(fpEvt->getGenCand(j),particles);
} // findGenStructure()

void massReader::findAllTrackIndices(TAnaCand* pCand, map<int,int> *indices)
{
	int j;
	
	// iterate through all own tracks. has to be done first, so the duplicate signal tracks
	// won't be added in the daughter anymore
	for (j = pCand->fSig1; j <= pCand->fSig2 && j>=0; j++)
		indices->insert(make_pair(fpEvt->getSigTrack(j)->fIndex,j));
	
	for (j = pCand->fDau1; j <= pCand->fDau2 && j>=0; j++)
		findAllTrackIndices(fpEvt->getCand(j),indices);
} // findAllTrackIndices()

float massReader::calculateIsolation(TAnaCand *pCand)
{
	const double cone = 0.7;
	const double pt_thres = 0.9;
	const double maxDocaSV = 0.05; // 500 um
	double sum_pt = 0.0;
	double isolation;
	int j,ntracks;
	map<int,int> candTracks;
	set<int> nearSV;
	TAnaTrack *pTrack;
	
	// build a set of all tracks near SV
	ntracks = pCand->fNstTracks.size();
	for (j = 0; j < ntracks; j++) {
		if (pCand->fNstTracks[j].second.first < maxDocaSV)
			nearSV.insert(pCand->fNstTracks[j].first);
	}
	
	// find the candidate tracks
	findAllTrackIndices(pCand,&candTracks);
	
	ntracks = fpEvt->nRecTracks();
	for (j = 0; j < ntracks; j++) {
		if (candTracks.count(j) > 0) // this track belongs to the candidate...
			continue;
		
		pTrack = fpEvt->getRecTrack(j);
		
		// check for sufficient pt of track
		if (pTrack->fPlab.Pt() <= pt_thres)
			continue;
		
		// check if track_j is within cone
		if (pCand->fPlab.DeltaR(pTrack->fPlab) >= cone)
			continue;
		
		// check that not the same PV
		if (pCand->fPvIdx != pTrack->fPvIdx && (pTrack->fPvIdx >= 0 || nearSV.count(j) == 0))
			continue;
		
		sum_pt += pTrack->fPlab.Pt();
	}
	
	isolation = pCand->fPlab.Pt();
	isolation = isolation / (isolation + sum_pt);
	
	return (float)isolation;
} // calculateIsolation()

float massReader::calculateDoca0(TAnaCand *pCand)
{
	TAnaTrack *pTrack;
	float doca0 = FLT_MAX;
	int j,ntracks = pCand->fNstTracks.size();
	
	for (j = 0; j < ntracks; j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		
		if (pTrack->fPvIdx < 0 || pTrack->fPvIdx == pCand->fPvIdx) {
			if (pCand->fNstTracks[j].second.first < doca0)
				doca0 = pCand->fNstTracks[j].second.first;
		}
	}
	
	return doca0;
} // calculateDoca0()

int massReader::countTracksNearby(TAnaCand *pCand)
{
	int result = 0;
	int j,ntracks = pCand->fNstTracks.size();
	const double pt_thres = 0.5;
	const double maxDocaSV = 0.03;
	TAnaTrack *pTrack;
	
	for (j = 0; j < ntracks; j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		
		if (pTrack->fPlab.Pt() <= pt_thres)
			continue;
		
		if (pCand->fNstTracks[j].second.first < maxDocaSV)
			result++;
	}
	
	return result;
} // countTracksNearby()

int massReader::loadTrigger(int *errTriggerOut, int *triggersFoundOut)
{
	unsigned j,k;
	trigger_table_t trg;
	int triggers = 0;
	int triggers_found = 0;
	int triggers_err = 0;
	
	for (j = 0; j < NHLT; j++) {
		
		std::string name(fpEvt->fHLTNames[j].Data());
		
		// check the name...
		for (k = 0; k < sizeof(g_trigger_table_signal)/sizeof(trigger_table_t); k++) {
			
			trg = g_trigger_table_signal[k];
			if (name.find(trg.trigger_name.c_str()) <= name.size()) {
				triggers_found |= trg.t_bit;
				if (fpEvt->fHLTError[j]) {
					triggers_err |= trg.t_bit;
					continue;
				}
				if (fpEvt->fHLTResult[j])
					triggers |= trg.t_bit;
			}
		}
		
		for (k = 0; k < sizeof(g_trigger_table_norm)/sizeof(trigger_table_t); k++) {
			
			trg = g_trigger_table_norm[k];
			if (name.find(trg.trigger_name.c_str()) <= name.size()) {
				triggers_found |= trg.t_bit;
				if (fpEvt->fHLTError[j]) {
					triggers_err |= trg.t_bit;
					continue;
				}
				if (fpEvt->fHLTResult[j])
					triggers |= trg.t_bit;
			}
		}
	}
	
	if (errTriggerOut) *errTriggerOut = triggers_err;
	if (triggersFoundOut) *triggersFoundOut = triggers_found;
	
	return triggers;
} // loadTrigger()

int massReader::isMuonTight(TAnaTrack *sigTrack)
{
	int result;
	TAnaTrack *recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
	TAnaMuon *mu;
	
	// high purity track
	result = (recTrack->fTrackQuality & 4) != 0;
	if (!result) goto bail;
	
	// muon::GlobalMuonPromptTight()
	result = (0 <= recTrack->fMuIndex) && (recTrack->fMuIndex < fpEvt->nMuons());
	if(!result) goto bail;
	
	// check for the muon id (muon::TrackerMuonArbitrated && muon::GlobalMuonPromptTight)
	result = ((recTrack->fMuID & 80) == 80);
	if (!result) goto bail;
	
	// Additional requirements
	mu = fpEvt->getMuon(recTrack->fMuIndex);
	
	//	* at least two muon stations are matched to the track
	result = mu->fTimeNdof >= 2;
	if (!result) goto bail;
	
	//	* the transverse impact parameter with respect to the beamspot is < 0.2cm
	result = (recTrack->fBsTip < 0.2);
	if (!result) goto bail;
	
	//	* at least one pixel detector hit be part of the muon track
	result = (numberOfPixLayers(recTrack) >= 1);
	if (!result) goto bail;
	
	//	* the muon track has at least 11 (strip or pixel) hits
	result = (recTrack->fValidHits >= 11);
	
bail:
	return result;
} // isMuonTight()

int massReader::hasTriggeredNorm()
{
	return hasTriggered(fTriggers, g_trigger_table_norm, sizeof(g_trigger_table_norm)/sizeof(trigger_table_t));
} // hasTriggeredNorm()

int massReader::hasTriggeredSignal()
{
	return hasTriggered(fTriggers, g_trigger_table_signal, sizeof(g_trigger_table_signal)/sizeof(trigger_table_t));
} // hasTriggeredSignal()

int massReader::hasTriggered(int triggers,trigger_table_t *table, unsigned size)
{
	unsigned j;
	int result = 0;
	
	for (j = 0; !result && j < size; j++) {
		
		if((table[j].t_bit & triggers) != 0) {
			// we triggered this one. Hence check the run range (or in case of mc always proceed)
			result = ((table[j].run_range.first <= fRun) && (fRun <= table[j].run_range.second))
			|| (fIsMC);
		}
	}
	
	return result;
} // hasTriggered()

int massReader::loadEfficiencyFlags(TGenCand *gen)
{
	int effFlags = kGeneratorCand; // the fact we're here already implies we're a generator cand.
	TGenCand *dau;
	TAnaTrack *track;
	map<int,int> genStruct; // (genIx,recTrackIx)
	map<int,int>::const_iterator it;
	bool muon;
	
	findGenStructure(gen,&genStruct);
	
	/* Check acceptance
	 *	both muons have high purity track with pt_gen > 1 && |eta_gen| < 2.5 (pt > 2.0, |eta < 2.4| on reco)
	 *	all kaons have high purity track with pt_gen > .4 && |eta_gen| < 2.5 (pt > 0.5, |eta < 2.4| on reco)
	 */
	for (it = genStruct.begin(); it != genStruct.end(); ++it) {
		
		dau = fpEvt->getGenCand(it->first);
		if (stableParticles.count(abs(dau->fID)) == 0) // only check stable particles
			continue;
		
		muon = (abs(dau->fID) == 13);
		
		// pt (muon > 3 && kaon > 0.5)
		if (dau->fP.Perp() <= (muon ? 1. : 0.4)) // not in pt acceptance
			goto bail;
		
		// eta
		if (fabs(dau->fP.Eta()) >= 2.5) // not in eta acceptance
			goto bail;
		
		if (it->second < 0) // no associated track
			goto bail;
		
		track = fpEvt->getRecTrack(it->second);
		if ((track->fTrackQuality & 0x1<<2) == 0) // no high purity track
			goto bail;
		
		// momentum and eta of reco tracks
		if (track->fPlab.Perp() <= (muon ? 2.0 : 0.5))
			goto bail;
		if (fabs(track->fPlab.Eta()) >= 2.4)
			goto bail;
	}
	effFlags |= kAcceptance; // candidate was in acceptance.
	
	/* Check muon efficiency
	 *	both muons are global && tracker muon
	 */
	for (it = genStruct.begin(); it != genStruct.end(); ++it) {
		
		dau = fpEvt->getGenCand(it->first);
		if (abs(dau->fID) != 13)
			continue;
		
		// check muon id
		track = fpEvt->getRecTrack(it->second);
		if (track->fMuIndex < 0 || (track->fMuID & 6) != 6)
			goto bail;
	}
	effFlags |= kEffMuon;
	
bail:
	return effFlags;
} // loadEfficiencyFlags()

void massReader::readCuts(TString filename, int dump)
{
	// parse the files...
	char buffer[1024];
	char cutName[128];
	float cutLow,cutHigh;
	int ok;
	
	FILE *cutFile = fopen(filename.Data(),"r");
	if(!cutFile) goto bail;
	
	if (dump) {
		cout << "===================================================" << endl;
		cout << "==> massReader: Cut file " << filename.Data() << endl;
		cout << "---------------------------------------------------" << endl;
	}
	
	fCutFileParsed = true;
	
	while (fgets(buffer,sizeof(buffer),cutFile)) {
		
		if (buffer[strlen(buffer)-1] == '\n')
			buffer[strlen(buffer)-1] = '\0';
		
		// skip comment line
		if (buffer[0] == '#') continue;
		
		ok = sscanf(buffer, "%s %f %f", cutName, &cutLow, &cutHigh);
		if(ok < 2) {
			cout << "==> massReader: Cut not parse input line:" << endl;
			cout << buffer << endl;
			cout << "==> massReader: abort..." << endl;
			abort();
			continue;
		}
		else if(ok == 2) // no upper cut read...
			cutHigh = 0.0f;
		
		if(!parseCut(cutName, cutLow, cutHigh)) {
			cout << "==> massReader: Error parsing variable " << cutName << ". abort..." << endl;
			abort();
		}
	}
	
	if (dump) cout << "---------------------------------------------------" << endl;
	
	// close the cut file...
bail:
	if (cutFile) fclose(cutFile);
} // readCuts()

bool massReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	bool parsed;
	// parse the options...
	parsed = (strcmp(cutName,"CAND") == 0);
	if (parsed) {
		fCutCand = (int)cutLow;
		if (dump) cout << "CAND: " << fCutCand << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"D3SIGN") == 0);
	if (parsed) {
		fCutFlight3dSign = cutLow;
		if (dump) cout << "D3SIGN: " << fCutFlight3dSign << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"CHI2") == 0);
	if (parsed) {
		fCutChi2 = cutLow;
		if (dump) cout << "CHI2: " << fCutChi2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT") == 0);
	if (parsed) {
		fCutPt = cutLow;
		if (dump) cout << "PT: " << fCutPt << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"ALPHA") == 0);
	if (parsed) {
		fCutAlpha = cutLow;
		if (dump) cout << "ALPHA: " << fCutAlpha << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRUTH") == 0);
	if (parsed) {
		fCutTruth = (int)cutLow;
		if (dump) cout << "TRUTH: " << fCutTruth << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRIGGERED") == 0);
	if (parsed) {
		fCutTriggered = (cutLow != 0.0f);
		if (dump) cout << "TRIGGERED: " << fCutTriggered << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"CHI2BYNDOF") == 0);
	if (parsed) {
		fCutChi2ByNdof = cutLow;
		if (dump) cout << "CHI2BYNDOF: " << fCutChi2ByNdof << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_MU1") == 0);
	if (parsed) {
		fCutTrackQual_mu1 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_MU1: " << fCutTrackQual_mu1 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_MU2") == 0);
	if (parsed) {
		fCutTrackQual_mu2 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_MU2: " << fCutTrackQual_mu2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_MU") == 0);
	if (parsed) {
		fCutOppSign_mu = (cutLow != 0.0f);
		if (dump) cout << "OPPOSITE_SIGN_MU: " << fCutOppSign_mu << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MUID") == 0);
	if (parsed) {
		fCutMuID_mask = (int)cutLow;
		fCutMuID_reqall = (cutHigh != 0.0);
		if (dump) cout << "MUID: mask " << fCutMuID_mask << " and require all " << fCutMuID_reqall << endl;
		goto bail;
	}
	
bail:
	return parsed;
} // parseCut()

bool massReader::applyCut()
{
	bool pass = true;
	
	// check the candidate
	if (fCutCand != 0 && fCandidate != 0) { // cut on candidate requested
		pass = fCandidate == fCutCand;
		if (!pass) goto bail;
	}
	
	if (fCutFlight3dSign > 0.0) {
		pass = fD3 / fD3E > fCutFlight3dSign;
		if (!pass) goto bail;
	}
	
	if (fCutChi2 > 0.0) {
		pass = fChi2 < fCutChi2;
		if (!pass) goto bail;
	}
	
	if (fCutPt > 0.0) {
		pass = fPt > fCutPt;
		if (!pass) goto bail;
	}
	
	if (fCutAlpha > 0.0) {
		pass = fAlpha < fCutAlpha;
		if (!pass) goto bail;
	}
	
	if (fCutTruth != 0) {
		pass = (fTruth == fCutTruth); // can be used to deselect all truth matched ones
		if (!pass) goto bail;
	}
	
	if (fCutChi2ByNdof > 0.0) {
		pass = fChi2 / fNdof < fCutChi2ByNdof;
		if (!pass) goto bail;
	}
	
	if (fCutTriggered) {
		pass = fTriggeredJPsi || fTriggeredBs;
		if (!pass) goto bail;
	}
	
	// check track quality of mu1
	if (fCutTrackQual_mu1) {
		pass = (fTrackQual_mu1 & fCutTrackQual_mu1) != 0;
		if (!pass) goto bail;
	}
	
	// check track quality of mu2
	if (fCutTrackQual_mu2) {
		pass = (fTrackQual_mu2 & fCutTrackQual_mu2) != 0;
		if (!pass) goto bail;
	}
	
	// check the Muon ID
	if (fCutMuID_mask != 0) {
		if (fCutMuID_reqall)
			pass = ((fMuID1 & fCutMuID_mask) == fCutMuID_mask) && ((fMuID2 & fCutMuID_mask) == fCutMuID_mask);
		else
			pass = ((fMuID1 & fCutMuID_mask) != 0) && ((fMuID2 & fCutMuID_mask) != 0);
		
		if (!pass) goto bail;
	}
	
	// check opposite sign of muons
	if (fCutOppSign_mu) {
		pass = fQ_mu1 != fQ_mu2;
		if (!pass) goto bail;
	}	
	
bail:
	return pass;
} // applyCut()
