#include "massReader.hh"
#include <cstdlib>

using namespace std;

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
	fCutTruth(0)
{	
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
		
		if(loadCandidateVariables(fpEvt->getCand(j)) && (!fCutFileParsed || applyCut()))
			reduced_tree->Fill();
	}
} // massReader::eventProcessing()

int massReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaTrack *pTrack;
	TAnaCand *momCand;
	TVector3 v1,v2,uVector;
	unsigned j,k;
	
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
	fIso7_pt0 = calculateIsolation(pCand, 0.7, 0.0, false);
	fIso7_pt5 = calculateIsolation(pCand, 0.7, 0.5, false);
	fIso7_pt7 = calculateIsolation(pCand, 0.7, 0.7, false);
	fIso7_pt10 = calculateIsolation(pCand, 0.7, 1.0, false);
	fIso10_pt0 = calculateIsolation(pCand, 1.0, 0.0, false);
	fIso10_pt5 = calculateIsolation(pCand, 1.0, 0.5, false);
	fIso10_pt7 = calculateIsolation(pCand, 1.0, 0.7, false);
	fIso10_pt9 = calculateIsolation(pCand, 1.0, 0.9, false);
	fIso10_pt10 = calculateIsolation(pCand, 1.0, 1.0, false);
	fIso7_pt0_pv = calculateIsolation(pCand, 0.7, 0.0, true);
	fIso7_pt5_pv = calculateIsolation(pCand, 0.7, 0.5, true);
	fIso7_pt7_pv = calculateIsolation(pCand, 0.7, 0.7, true);
	fIso7_pt10_pv = calculateIsolation(pCand, 0.7, 1.0, true);
	fIso10_pt0_pv = calculateIsolation(pCand, 1.0, 0.0, true);
	fIso10_pt5_pv = calculateIsolation(pCand, 1.0, 0.5, true);
	fIso10_pt7_pv = calculateIsolation(pCand, 1.0, 0.7, true);
	fIso10_pt9_pv = calculateIsolation(pCand, 1.0, 0.9, true);
	fIso10_pt10_pv = calculateIsolation(pCand, 1.0, 1.0, true);
	fCtau = 0.0f; // we can compute this only in a subclass, so initialize to zero
	fEta = pCand->fPlab.Eta();
	
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
		v2 = pCand->fVtx.fPoint - fpEvt->bestPV()->fPoint;
	}
	
	fAlpha = v1.Angle(v2);
	fD3_Para = v1.Dot(v2) / v1.Mag(); // projection of v2 to v1
	fD3_Perp = (v2 - (fD3_Para / v1.Mag() * v1) ).Mag(); // perpendicular part of v2 w.r.t v1
	
	// project to xy plane
	v1.SetZ(0);
	v2.SetZ(0);
	fAlphaXY = v1.Angle(v2);
	
	fDxy_Para = v1.Dot(v2) / v1.Mag(); // projection of v2 to v1
	fDxy_Para = (v2 - (fD3_Para / v1.Mag() * v1) ).Mag(); // perpendicular part of v2 w.r.t v1
	
	// do this at the end so the checkTruth algorithm can use
	// all variables of this candidate.
	fTruth = checkTruth(pCand);
	fTruthFlags = loadTruthFlags(pCand);
	fTriggers = loadTrigger(&fTriggersError,&fTriggersFound);
	
	return 1;
} // loadCandidateVariables()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T",fTreeName);
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("pt",&fPt,"pt/F");
	reduced_tree->Branch("mass",&fMass,"mass/F");
	reduced_tree->Branch("mass_c",&fMassConstraint,"mass_c/F");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
	reduced_tree->Branch("truth_flags",&fTruthFlags,"truth_flags/I");
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
	reduced_tree->Branch("iso7_pt0",&fIso7_pt0,"iso7_pt0/F");
	reduced_tree->Branch("iso7_pt5",&fIso7_pt5,"iso7_pt5/F");
	reduced_tree->Branch("iso7_pt7",&fIso7_pt7,"iso7_pt7/F");
	reduced_tree->Branch("iso7_pt10",&fIso7_pt10,"iso7_pt10/F");
	reduced_tree->Branch("iso10_pt0",&fIso10_pt0,"iso10_pt0/F");
	reduced_tree->Branch("iso10_pt5",&fIso10_pt5,"iso10_pt5/F");
	reduced_tree->Branch("iso10_pt7",&fIso10_pt7,"iso10_pt7/F");
	reduced_tree->Branch("iso10_pt9",&fIso10_pt9,"iso10_pt9/F");
	reduced_tree->Branch("iso10_pt10",&fIso10_pt10,"iso10_pt10/F");
	reduced_tree->Branch("iso7_pt0_pv",&fIso7_pt0_pv,"iso7_pt0_pv/F");
	reduced_tree->Branch("iso7_pt5_pv",&fIso7_pt5_pv,"iso7_pt5_pv/F");
	reduced_tree->Branch("iso7_pt7_pv",&fIso7_pt7_pv,"iso7_pt7_pv/F");
	reduced_tree->Branch("iso7_pt10_pv",&fIso7_pt10_pv,"iso7_pt10_pv/F");
	reduced_tree->Branch("iso10_pt0_pv",&fIso10_pt0_pv,"iso10_pt0_pv/F");
	reduced_tree->Branch("iso10_pt5_pv",&fIso10_pt5_pv,"iso10_pt5_pv/F");
	reduced_tree->Branch("iso10_pt7_pv",&fIso10_pt7_pv,"iso10_pt7_pv/F");
	reduced_tree->Branch("iso10_pt9_pv",&fIso10_pt9_pv,"iso10_pt9_pv/F");
	reduced_tree->Branch("iso10_pt10_pv",&fIso10_pt10_pv,"iso10_pt10_pv/F");
	reduced_tree->Branch("triggers",&fTriggers,"triggers/I");
	reduced_tree->Branch("triggers_error",&fTriggersError,"triggers_error/I");
	reduced_tree->Branch("triggers_found",&fTriggersFound,"triggers_found/I");
	reduced_tree->Branch("ctau",&fCtau,"ctau/F");
	reduced_tree->Branch("eta",&fEta,"eta/F");
	reduced_tree->Branch("d3_perp",&fD3_Perp,"d3_perp/F");
	reduced_tree->Branch("d3_para",&fD3_Para,"d3_para/F");
	reduced_tree->Branch("dxy_perp",&fDxy_Perp,"dxy_perp/F");
	reduced_tree->Branch("dxy_para",&fDxy_Para,"dxy_para/F");
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

float massReader::calculateIsolation(TAnaCand *pCand, double openingAngle, double minPt, bool sameVertex)
{
	double iso; // calculate in double precision and return single
	TVector3 plabB;
	int j,ntracks;
	TAnaTrack *pTrack;
	double r;
	map<int,int> usedTracks;
	
	plabB = pCand->fPlab;
	findAllTrackIndices(pCand,&usedTracks);
	
	iso = 0.0;
	ntracks = fpEvt->nRecTracks();
	for (j = 0; j < ntracks; j++) {
		
		if (usedTracks.count(j) > 0) continue; // this track belongs to the candidate
		
		pTrack = fpEvt->getRecTrack(j);
		
		// if samevertex, veto tracks not from the same PV
		if (sameVertex && (pTrack->fPvIdx != pCand->fPvIdx)) continue;
		
		// pt veto
		if (pTrack->fPlab.Perp() <= minPt) continue;
		
		// opening angle
		r = plabB.DeltaR(pTrack->fPlab);
		if (r < openingAngle)
			iso += pTrack->fPlab.Perp();
	}
	
	iso = plabB.Perp() / (plabB.Perp() + iso);
	
	return (float)iso;
} // calculateIsolation()

int massReader::loadTrigger(int *errTriggerOut, int *triggersFoundOut)
{
	unsigned j;
	int triggers = 0;
	TString string;
	int triggers_found = 0; // store the trigger we found!!
	int triggers_err = 0;
	
	for (j = 0; j < NHLT; j++) {
		if (fpEvt->fHLTNames[j] == "HLT_DoubleMu3") {
			triggers_found |= kHLT_DoubleMu3_Bit;
			
			if (fpEvt->fHLTError[j]) {
				triggers_err |= kHLT_DoubleMu3_Bit;
				continue;
			}
			
			if(fpEvt->fHLTResult[j])
				triggers |= kHLT_DoubleMu3_Bit;
		} else if (fpEvt->fHLTNames[j] == "HLT_DoubleMu0") {
			triggers_found |= kHLT_DoubleMu0_Bit;
			
			if (fpEvt->fHLTError[j]) {
				triggers_err |= kHLT_DoubleMu0_Bit;
				continue;
			}
			
			if(fpEvt->fHLTResult[j])
				triggers |= kHLT_DoubleMu0_Bit;
		} else if (fpEvt->fHLTNames[j] == "HLT_DoubleMu0_Quarkonium_v1") {
			triggers_found |= kHLT_DoubleMu0_Quarkonium_v1_Bit;
			
			if (fpEvt->fHLTError[j]) {
				triggers_err |= kHLT_DoubleMu0_Quarkonium_v1_Bit;
				continue;
			}
			
			if (fpEvt->fHLTResult[j])
				triggers |= kHLT_DoubleMu0_Quarkonium_v1_Bit;
		}
	}
	
	if (errTriggerOut)
		*errTriggerOut = triggers_err;

	if (triggersFoundOut)
		*triggersFoundOut = triggers_found;
	
	return triggers;
} // loadTrigger()

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
	
bail:
	return parsed;
} // parseCut()

bool massReader::applyCut()
{
	bool pass = true;
	
	// check the candidate
	if (fCutCand != 0) { // cut on candidate requested
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
		pass = ( (fTriggersFound & kHLT_DoubleMu0_Quarkonium_v1_Bit) && (fTriggers & kHLT_DoubleMu0_Quarkonium_v1_Bit) ) || ( ((fTriggersFound & kHLT_DoubleMu0_Quarkonium_v1_Bit) == 0) && (fTriggers & kHLT_DoubleMu0_Bit) );
		if (!pass) goto bail;
	}
	
bail:
	return pass;
} // applyCut()
