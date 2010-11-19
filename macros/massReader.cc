#include "massReader.hh"
#include <cstdlib>

using namespace std;

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),
	reduced_tree(NULL),
	fCutFileParsed(false),
	fCutCand(0),
	fCutFlight3dSign(0.0),
	fCutChi2(0.0),
	fCutPt(0.0),
	fCutAlpha(0.0),
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
	TAnaCand *momCand;
	TVector3 v1,v2;
	
	// Save in the tree
	fCandidate = pCand->fType;
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
	fIso7_pt0 = calculateIsolation(pCand, 0.7, 0.0);
	fIso7_pt5 = calculateIsolation(pCand, 0.7, 0.5);
	fIso7_pt7 = calculateIsolation(pCand, 0.7, 0.7);
	fIso7_pt10 = calculateIsolation(pCand, 0.7, 1.0);
	fIso10_pt0 = calculateIsolation(pCand, 1.0, 0.0);
	fIso10_pt5 = calculateIsolation(pCand, 1.0, 0.5);
	fIso10_pt7 = calculateIsolation(pCand, 1.0, 0.7);
	fIso10_pt10 = calculateIsolation(pCand, 1.0, 1.0);
	fCtau = 0.0; // we can compute this only in a subclass, so initialize to zero
	fEta = pCand->fPlab.Eta();
	
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
	reduced_tree->Branch("iso7_pt0",&fIso7_pt0,"iso7_pt0/F");
	reduced_tree->Branch("iso7_pt5",&fIso7_pt5,"iso7_pt5/F");
	reduced_tree->Branch("iso7_pt7",&fIso7_pt7,"iso7_pt7/F");
	reduced_tree->Branch("iso7_pt10",&fIso7_pt10,"iso7_pt10/F");
	reduced_tree->Branch("iso10_pt0",&fIso10_pt0,"iso10_pt0/F");
	reduced_tree->Branch("iso10_pt5",&fIso10_pt5,"iso10_pt5/F");
	reduced_tree->Branch("iso10_pt7",&fIso10_pt7,"iso10_pt7/F");
	reduced_tree->Branch("iso10_pt10",&fIso10_pt10,"iso10_pt10/F");
	reduced_tree->Branch("triggers",&fTriggers,"triggers/I");
	reduced_tree->Branch("triggers_error",&fTriggersError,"triggers_error/I");
	reduced_tree->Branch("triggers_found",&fTriggersFound,"triggers_found/I");
	reduced_tree->Branch("ctau",&fCtau,"ctau/F");
	reduced_tree->Branch("eta",&fEta,"eta/F");
	reduced_tree->Branch("d3_perp",&fD3_Perp,"d3_perp/F");
	reduced_tree->Branch("d3_para",&fD3_Para,"d3_para/F");
	reduced_tree->Branch("dxy_perp",&fDxy_Perp,"dxy_perp/F");
	reduced_tree->Branch("dxy_para",&fDxy_Para,"dxy_para/F");
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

float massReader::calculateIsolation(TAnaCand *pCand, double openingAngle, double minPt)
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
		if (pTrack->fPlab.Perp() <= minPt) continue;
		
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
	char buffer[256];
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
		
		if(buffer[strlen(buffer)-1] != '\n') {
			cout << "==> massReader: Cut file line is too long, skipping..." << endl;
			continue;
		}
		
		// skip comment line
		if (buffer[0] == '#') continue;
		
		ok = sscanf(buffer, "%s %f %f", cutName, &cutLow, &cutHigh);
		if(ok < 2) {
			cout << "==> massReader: Cut not parse input line:" << endl;
			cout << buffer << endl;
			cout << "==> massReader: Skipping this line..." << endl;
			continue;
		}
		else if(ok == 2) // no upper cut read...
			cutHigh = 0.0f;
		
		if(!parseCut(cutName, cutLow, cutHigh))
			cout << "==> massReader: Error parsing variable " << cutName << ". Ignoring this cut..." << endl;
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
	
bail:
	return pass;
} // applyCut()
