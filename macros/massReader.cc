#include "massReader.hh"

#include <utility>
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

static decay_t make_decay(int nbr, ...)
{
	decay_t dec;
	va_list ids;
	int id,i;
	
	va_start(ids,nbr);
	for (i = 0; i < nbr; i++) {
		id = va_arg(ids,int);
		dec.insert(id);
	}
	va_end(ids);
	
	return dec;
} // make_decay()

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),
	reduced_tree(NULL),
	fCutFileParsed(false),
	fCutTriggered(false),
	fCutCand(0),
	fCutFlight3dSign(0.0),
	fCutChi2(0.0),
	fCutPt(0.0),
	fCutAlpha(0.0),
	fCutChi2ByNdof(0.0),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutMuID_reqall(false),
	fCutOppSign_mu(false),
	fCutMass_JPsiLow(0.0),
	fCutMass_JPsiHigh(0.0),
	fCutTrackQual_kp1(0),
	fCutTrackQual_kp2(0),
	fCutPt_Kaon1(0.0),
	fCutPt_Kaon2(0.0),
	fCutOppSign_kp(false),
	fCutMass_PhiLow(0.0),
	fCutMass_PhiHigh(0.0)
{
	fTreeName = "massReader reduced tree";
	
	// add only the ones we actually need.
	// feel free to enlarge this set.
	stableParticles.insert(13); // muon
	stableParticles.insert(22); // photon
	stableParticles.insert(321); // kaon
	stableParticles.insert(211); // pion+
	stableParticles.insert(111); // pion0
	
	// build the decays we're interested in...
	decayTable.insert( pair<decay_t,int>(make_decay(3, 531, 13, 13), kDecay_BsToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(4, 531, 13, 13, 22), kDecay_BsToMuMuGa) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 531, 321, 321), kDecay_BsToKK) );	
	decayTable.insert( pair<decay_t,int>(make_decay(3, 531, 321, 211), kDecay_BsToKPi) );	
	decayTable.insert( pair<decay_t,int>(make_decay(3, 531, 211, 211), kDecay_BsToPiPi) );
	decayTable.insert( pair<decay_t,int>(make_decay(4, 531, 211, 13, 14), kDecay_BsToPiMuNu) );
	decayTable.insert( pair<decay_t,int>(make_decay(4, 531, 321, 13, 14), kDecay_BsToKMuNu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 511, 13, 13), kDecay_BdToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 511, 211, 211), kDecay_BdToPiPi) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 511, 321, 211), kDecay_BdToKPi) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 511, 321, 321), kDecay_BdToKK) );
	decayTable.insert( pair<decay_t,int>(make_decay(4, 511, 13, 13, 111), kDecay_BdToMuMuPi0) );	
	decayTable.insert( pair<decay_t,int>(make_decay(4, 511, 211, 13, 14), kDecay_BdToPiMuNu) );	
	decayTable.insert( pair<decay_t,int>(make_decay(5, 521, 13, 13, 13, 14), kDecay_BuTo3MuNu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 5122, 2212, 211), kDecay_LambdaBToPPi) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 5122, 2212, 321), kDecay_LambdaBToPK) );
	decayTable.insert( pair<decay_t,int>(make_decay(4, 5122, 2212, 13, 14), kDecay_LambdaBToPMuNu) );
	decayTable.insert( pair<decay_t,int>(make_decay(7, 531, 443, 333, 13, 13, 321, 321), kDecay_Bs2JpsiPhi) );
	decayTable.insert( pair<decay_t,int>(make_decay(5, 521, 443, 13, 13, 321), kDecay_Bu2JpsiKp) );
	decayTable.insert( pair<decay_t,int>(make_decay(7, 511, 443, 13, 13, 313, 321, 211), kDecay_Bd2JpsiKstar) );
	decayTable.insert( pair<decay_t,int>(make_decay(7, 511, 443, 13, 13, 310, 211, 211), kDecay_Bd2JpsiKs) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 443, 13, 13), kDecay_PsiToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 100443, 13, 13), kDecay_Psi2SToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 553, 13, 13), kDecay_Ups1SToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 100553, 13, 13), kDecay_Ups2SToMuMu) );
	decayTable.insert( pair<decay_t,int>(make_decay(3, 200553, 13, 13), kDecay_Ups3SToMuMu) );
	
	// mothers we're interested in
	validMothers.insert(511); // Bd
	validMothers.insert(521); // Bu
	validMothers.insert(531); // Bs
	validMothers.insert(5122); // Lambdab
} // massReader()

massReader::~massReader()
{ } // ~massReader()

void massReader::eventProcessing()
{
	int j,nc;
	
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
	fSameMother = 0;
	fTrueDecay = 0;
	fPt = 0.0;
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
	fEtaMu1 = fEtaMu2 = 0.0f;
	fMuTight1 = fMuTight2 = 0;
	fTrackQual_mu1 = fTrackQual_mu2 = 0;
	fQ_mu1 = fQ_mu2 = 0;
	fDeltaR = 0.0f; // deltaR
	fNbrPV = 0;
	fDeltaPhiMu = 0;
	fIPCand = 0.0f;
	fIPCandE = 0.0f;
	
	// jpsi
	fPtJPsi = 0.0;
	fMassJPsi = 0.0;
	fChi2Jpsi = 0.0;
	
	// phi
	fMassPhi = 0.0f;
	
	// kaons
	fPtKp1 = fPtKp2 = 0.0;
	fEtaKp1 = fEtaKp2 = 0.0;
	fPtKp_Gen1 = fPtKp_Gen2 = 0.0;
	fEtaKp_Gen1 = fEtaKp_Gen2 = 0.0;
	fTrackQual_kp1 = fTrackQual_kp2 = 0;
	fQ_kp1 = fQ_kp2 = 0;
	fDeltaR_Kaons = 0.0f;
	
	memset(fTracksIx,0,sizeof(fTracksIx));
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksIPE,0,sizeof(fTracksIPE));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
} // clearVariables()

int massReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaTrack *pTrack;
	TAnaTrack *sigTrack;
	TAnaTrack *recTrack;
	TAnaCand *momCand;
	TAnaCand *dauCand;
	TGenCand *muGen;
	TVector3 v1,v2,uVector;
	int j,k;
	bool firstMu = true;
	bool firstKaon = true;
	map<int,int> aTracks;
	map<int,int> cand_tracks;
	TVector3 plabMu1,plabMu2;
	TVector3 plabKp1,plabKp2;
	TLorentzVector pMu1,pMu2;

	int result;
	
	clearVariables();
	
	// Save in the tree
	fCandidate = pCand->fType;
	fPt = pCand->fPlab.Perp();
	fMass = pCand->fMass;
	fMassConstraint = -1.0f;
	fD3 = pCand->fVtx.fD3d;
	fD3E = pCand->fVtx.fD3dE;
	fDxy = pCand->fVtx.fDxy;
	fDxyE = pCand->fVtx.fDxyE;
	fChi2 = pCand->fVtx.fChi2;
	fNdof = pCand->fVtx.fNdof;
	fMaxDoca = pCand->fMaxDoca;
	fNbrPV = fpEvt->nPV();
	fIPCand = TMath::Sqrt(pCand->fPvLip*pCand->fPvLip + pCand->fPvTip*pCand->fPvTip);
	fIPCandE = TMath::Sqrt((pCand->fPvLip*pCand->fPvLip)/(fIPCand*fIPCand)*(pCand->fPvLipE*pCand->fPvLipE) + (pCand->fPvTip*pCand->fPvTip)/(fIPCand*fIPCand)*(pCand->fPvTipE*pCand->fPvTipE));
	fIsoMoriond12 = calculateIsolation(pCand);
	fDoca0 = calculateDoca0(pCand);
	fNbrNearby = countTracksNearby(pCand);
	
	fCtau = pCand->fTau3d;
	fCtauE = pCand->fTau3dE;
	fEta = pCand->fPlab.Eta();
	
	// set the constraint mass value
	findCandStructure(pCand, &cand_tracks);
	dauCand = findCandidate(400000 + (pCand->fType % 1000), &cand_tracks);
	if (dauCand) fMassConstraint = dauCand->fMass;
	
	// Load muon variables
	cand_tracks.clear();
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it!=cand_tracks.end(); ++it) {
		
		sigTrack = fpEvt->getSigTrack(it->second);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		
		switch (abs(sigTrack->fMCID)) {
			case 13: // muon
				if (firstMu) {
					plabMu1 = sigTrack->fPlab;
					fPtMu1 = sigTrack->fPlab.Perp();
					fTrackQual_mu1 = recTrack->fTrackQuality;
					fEtaMu1 = sigTrack->fPlab.Eta();
					fQ_mu1 = recTrack->fQ;
					fMuTight1 = isMuonTight(sigTrack);
				} else {
					plabMu2 = sigTrack->fPlab;
					fPtMu2 = sigTrack->fPlab.Perp();
					fTrackQual_mu2 = recTrack->fTrackQuality;
					fEtaMu2 = sigTrack->fPlab.Eta();
					fQ_mu2 = recTrack->fQ;
					fMuTight2 = isMuonTight(sigTrack);
				}
				firstMu = false;
				break;
			case 321: // kaon
				if (firstKaon) {
					plabKp1 = sigTrack->fPlab;
					fPtKp1 = plabKp1.Perp();
					fEtaKp1 = plabKp1.Eta();
					fQ_kp1 = recTrack->fQ;
					fTrackQual_kp1 = recTrack->fTrackQuality;					
				} else {
					plabKp2 = sigTrack->fPlab;
					fPtKp2 = plabKp2.Perp();
					fEtaKp2 = plabKp2.Eta();
					fQ_kp2 = recTrack->fQ;
					fTrackQual_kp2 = recTrack->fTrackQuality;
				}
				firstKaon = false;
				break;
			default:
				break;
		}
	}
	
	// make sure mu1 is leading and mu2 subleading
	if (fPtMu1 < fPtMu2) {
		swap(plabMu1,plabMu2);
		swap(fPtMu1,fPtMu2);
		swap(fEtaMu1,fEtaMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
		swap(fMuTight1,fMuTight2);
	}
	
	fDeltaPhiMu = plabMu1.DeltaPhi(plabMu2);
	
	// kaon ordering
	if (fPtKp1 < fPtKp2) {
		swap(plabKp1,plabKp2);
		swap(fPtKp1,fPtKp2);
		swap(fEtaKp1,fEtaKp2);
		swap(fQ_kp1,fQ_kp2);
		swap(fTrackQual_kp1,fTrackQual_kp2);
	}
	
	// Clean entries of nearest tracks
	for (j = 0; j < NBR_TRACKS_STORE; j++) fTracksIx[j] = -1;
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksIPE,0,sizeof(fTracksIPE));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
	
	// fill the histogramms...
	for (j = 0, k = 0; j < (int)pCand->fNstTracks.size(); j++) {
		
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
	fSameMother = sameMother(pCand) >= 0;
	fTrueDecay = loadDecay(pCand);
	fTriggers = loadTrigger(&fTriggersError,&fTriggersFound);
	fTriggeredJPsi = hasTriggeredNorm();
	fTriggeredBs = hasTriggeredSignal();
	
	// generator tracks...
	firstMu = true;
	firstKaon = true;
	findAllTrackIndices(pCand,&aTracks);
	for (map<int,int>::const_iterator it = aTracks.begin(); it != aTracks.end(); ++it) {
		pTrack = fpEvt->getRecTrack(it->first);
		if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= fpEvt->nGenCands())
			continue;
		
		muGen = fpEvt->getGenCand(pTrack->fGenIndex);
		switch (abs(muGen->fID)) {
			case 13: // muon
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
				break;
			case 321: // kaon
				if (firstKaon) {
					fPtKp_Gen1 = muGen->fP.Perp();
					fEtaKp_Gen1 = muGen->fP.Eta();
				} else {
					fPtKp_Gen1 = muGen->fP.Perp();
					fEtaKp_Gen1 = muGen->fP.Eta();
				}
				firstKaon = false;
				break;
			default:
				break;
		}
	}
	
	// muon generator ordering
	if (fPtMu1_Gen < fPtMu2_Gen) {
		swap(fPtMu1_Gen,fPtMu2_Gen);
		swap(fEtaMu1_Gen,fEtaMu2_Gen);
		swap(fPMu1_Gen,fPMu2_Gen);
	}
		
	// dimuon
	pMu1.SetVectM(plabMu1,MMUON);
	pMu2.SetVectM(plabMu2,MMUON);
	fPtJPsi = (pMu1 + pMu2).Perp();
	fMassJPsi = (pMu1 + pMu2).M();
	if (plabMu1.Perp() > 0 && plabMu2.Perp() > 0)
		fDeltaR = plabMu1.DeltaR(plabMu2);
	
	// dikaon
	pMu1.SetVectM(plabKp1,MKAON);
	pMu2.SetVectM(plabKp2,MKAON);
	fMassPhi = (pMu1 + pMu2).M();
	if (plabKp1.Perp() > 0 && plabKp2.Perp() > 0)
		fDeltaR_Kaons = (float)plabKp1.DeltaR(plabKp2);
	
	// daughters
	for (j = pCand->fDau1; 0 <= j && j <= pCand->fDau2; j++) {
		
		dauCand = fpEvt->getCand(j);
		if ((dauCand->fType % 1000) == 443) {
			
			fPtJPsi = dauCand->fPlab.Perp();
			fMassJPsi = dauCand->fMass;
			fChi2Jpsi = dauCand->fVtx.fChi2;
			break;
		} else if ((dauCand->fType % 1000) == 333) {
			fMassPhi = dauCand->fMass;
		}
	}
	
	if (fCandidate != 301313) // maybe run blind...
		result = 1;
	else
		result = !(BLIND && 5.2 < pCand->fMass && pCand->fMass < 5.45);
	
	return result;
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
	reduced_tree->Branch("same_mother",&fSameMother,"same_mother/I");
	reduced_tree->Branch("true_decay",&fTrueDecay,"true_decay/I");
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
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
	reduced_tree->Branch("pt_dimuon",&fPtJPsi,"pt_dimuon/F");
	reduced_tree->Branch("mass_dimuon",&fMassJPsi,"mass_dimuon/F");
	reduced_tree->Branch("pt_kp1",&fPtKp1,"pt_kp1/F");
	reduced_tree->Branch("pt_kp2",&fPtKp2,"pt_kp2/F");
	reduced_tree->Branch("pt_kp1_gen",&fPtKp_Gen1,"pt_kp1_gen/F");
	reduced_tree->Branch("pt_kp2_gen",&fPtKp_Gen2,"pt_kp2_gen/F");
	reduced_tree->Branch("eta_kp1",&fEtaKp1,"eta_kp1/F");
	reduced_tree->Branch("eta_kp2",&fEtaKp2,"eta_kp2/F");
	reduced_tree->Branch("eta_kp1_gen",&fEtaKp_Gen1,"eta_kp1_gen/F");
	reduced_tree->Branch("eta_kp2_gen",&fEtaKp_Gen2,"eta_kp2_gen/F");
	reduced_tree->Branch("track_qual_kp1",&fTrackQual_kp1,"track_qual_kp1/I");
	reduced_tree->Branch("track_qual_kp2",&fTrackQual_kp2,"track_qual_kp2/I");
	reduced_tree->Branch("q_kp1",&fQ_kp1,"q_kp1/I");
	reduced_tree->Branch("q_kp2",&fQ_kp2,"q_kp2/I");
	reduced_tree->Branch("mass_dikaon",&fMassPhi,"mass_dikaon/F");
	reduced_tree->Branch("deltaR_kaons",&fDeltaR_Kaons,"deltaR_kaons/F");
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

int massReader::sameMother(TAnaCand *cand)
{
	int result = -1;
	TAnaTrack *pTrack;
	TGenCand *motherParticle = NULL;
	TGenCand *trackParticle;
	int nGens, j;
		
	nGens = fpEvt->nGenCands();
	for (j = cand->fSig1; j >= 0 && j <= cand->fSig2; j++) {
		pTrack = fpEvt->getSigTrack(j);
		pTrack = fpEvt->getRecTrack(pTrack->fIndex);
		
		// get the generator info
		if (pTrack->fGenIndex < 0 || pTrack->fGenIndex >= nGens) goto bail;
		
		// look for the mother particle
		trackParticle = fpEvt->getGenCand(pTrack->fGenIndex);
		while (validMothers.count(abs(trackParticle->fID)) == 0) {
			if (trackParticle->fMom1 < 0 || trackParticle->fMom1 >= nGens) goto bail;
			trackParticle = fpEvt->getGenCand(trackParticle->fMom1);
		}
		
		if (motherParticle) {
			if (motherParticle->fNumber != trackParticle->fNumber)
				goto bail;
		} else {
			motherParticle = trackParticle;
		}

	}
	
	// still here? all particle from the same 'valid' mother
	if (motherParticle) result = motherParticle->fNumber;
bail:
	return result;
} // checkTruth()

int massReader::loadDecay(TAnaCand *anaCand)
{
	TGenCand *mother;
	int ix, result = 0;
	decay_t dec;
	map<decay_t,int>::const_iterator it;
	bool correct;
	TAnaTrack *sigTrack,*recTrack;
	
	// search for the decay..
	if ((ix = sameMother(anaCand)) >= 0)
		mother = fpEvt->getGenCand(ix);
	else goto bail;
	
	buildDecay(mother, &dec);
	
	it = decayTable.find(dec);
	if (it != decayTable.end()) {
		// check wether all tracks are assigned correctly
		correct = true;
		for (ix = anaCand->fSig1; correct && 0 <= ix && ix <= anaCand->fSig2; ix++) {
			sigTrack = fpEvt->getSigTrack(ix);
			recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
			correct = (abs(sigTrack->fMCID) == abs(recTrack->fMCID)); // track was assigned correctly
		}
		if (correct)
			result = it->second;
	}
bail:
	return result;
} // loadDecay()

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
		
		if ((pCand->fPvIdx == pTrack->fPvIdx) || (nearSV.count(j)>0 && pTrack->fPvIdx < 0))
			sum_pt += pTrack->fPlab.Pt();


		if ((pCand->fPvIdx != pTrack->fPvIdx) && (nearSV.count(j)==0 || pTrack->fPvIdx >= 0))
			continue;
				
		if ((pCand->fPvIdx != pTrack->fPvIdx) && (pTrack->fPvIdx >= 0 || nearSV.count(j) == 0))
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
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP1") == 0);
	if (parsed) {
		fCutTrackQual_kp1 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP1: " << fCutTrackQual_kp1 << endl;
		goto bail;
	}
	parsed = (strcmp(cutName,"TRACK_QUAL_KP2") == 0);
	if (parsed) {
		fCutTrackQual_kp2 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP2: " << fCutTrackQual_kp2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MASS_JPSI") == 0);
	if (parsed) {
		fCutMass_JPsiLow = cutLow;
		fCutMass_JPsiHigh = cutHigh;
		if (dump) cout << "MASS_JPSI: (" << fCutMass_JPsiLow << ", " << fCutMass_JPsiHigh << ")" << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT_KAON1") == 0);
	if (parsed) {
		fCutPt_Kaon1 = cutLow;
		if (dump) cout << "PT_KAON1: " << fCutPt_Kaon1 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT_KAON2") == 0);
	if (parsed) {
		fCutPt_Kaon2 = cutLow;
		if (dump) cout << "PT_KAON2: " << fCutPt_Kaon2 << endl;
		goto bail;
	}
	
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_KP") == 0);
	if (parsed) {
		fCutOppSign_kp = (cutLow != 0.0f);
		if (dump) cout << "OPPOSITE_SIGN_KP: " << fCutOppSign_kp << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MASS_PHI") == 0);
	if (parsed) {
		fCutMass_PhiLow = cutLow;
		fCutMass_PhiHigh = cutHigh;
		if (dump) cout << "MASS_PHI: (" << fCutMass_PhiLow << ", " << fCutMass_PhiHigh << ")" << endl;
		goto bail;
	}
	
bail:
	return parsed;
} // parseCut()

bool massReader::applyCut()
{
	bool pass = true;
	
	// check the candidate
	if (fCutCand != 0) {
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
	
	// check opposite sign of muons
	if (fCutOppSign_mu) {
		pass = fQ_mu1 != fQ_mu2;
		if (!pass) goto bail;
	}	
	
	// check track quality of kp
	if (fCutTrackQual_kp1) {
		pass = (fTrackQual_kp1 & fCutTrackQual_kp1) != 0;
		if (!pass) goto bail;
	}
	
	if (fCutTrackQual_kp2) {
		pass = (fTrackQual_kp2 & fCutTrackQual_kp2) != 0;
		if (!pass) goto bail;
	}
	
	if (fCutMass_JPsiLow > 0.0) {
		pass = fMassJPsi > fCutMass_JPsiLow;
		if(!pass) goto bail;
	}
	
	if (fCutMass_JPsiHigh > 0.0) {
		pass = fMassJPsi < fCutMass_JPsiHigh;
		if(!pass) goto bail;
	}
	
	if (fCutPt_Kaon1 > 0.0) {
		pass = fPtKp1 > fCutPt_Kaon1;
		if (!pass) goto bail;
	}
	
	if (fCutPt_Kaon2 > 0.0) {
		pass = fPtKp2 > fCutPt_Kaon2;
		if (!pass) goto bail;
	}
	
	if (fCutOppSign_kp) {
		pass = fQ_kp1 != fQ_kp2;
		if (!pass) goto bail;
	}
	
	if (fCutMass_PhiLow > 0.0) {
		pass = fMassPhi > fCutMass_PhiLow;
		if (!pass) goto bail;
	}
	
	if (fCutMass_PhiHigh > 0.0) {
		pass = fMassPhi < fCutMass_PhiHigh;
		if (!pass) goto bail;
	}
	
bail:
	return pass;
} // applyCut()
