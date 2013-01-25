/*
 *  pixelReader.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 30.07.12.
 *
 */

#include "pixelReader.hh"

#include <utility>

const double kMuMToCM = 1E-4;

bool operator<(res_t r1, res_t r2)
{
	return r1.p < r2.p;
} // operator<()

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

pixelReader::pixelReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName), reduced_tree(NULL), fD0Resolution(0.0), fDzResolution(0.0), fPhiResolution(0.0), fCotThetaResolution(0.0), fPtResolution(0.0), fPVResolutionXY(0.0), fPVResolutionZ(0.0), fNumCands(1), fResMode(kResFile)
{
	fStableParticles.insert(11); // e
	fStableParticles.insert(13); // mu
	fStableParticles.insert(22); // gamma
	fStableParticles.insert(211); // pi+
	fStableParticles.insert(321); // K+
	fStableParticles.insert(2212); // p
	fStableParticles.insert(2112); // n
	
	// build the decays we're interested in...
	fDecayTable.insert( std::pair<decay_t,int>(make_decay(3, 531, 13, 13), kDecay_BsToMuMu) );
	
	// eta resolution bins
	fEtaBins[0] = 1.0;
	fEtaBins[1] = 1.5;
	fEtaBins[2] = 2.0;
	fEtaBins[3] = 2.5;
} // pixelReader()

pixelReader::~pixelReader()
{} // ~pixelReader()
	
void pixelReader::bookHist()
{
	using std::cout; using std::endl;
	unsigned j,k;
	res_t res;
	double *pnt = NULL;
	
	reduced_tree = new TTree("T","Pixel2012 Reduced Tree");
	
	reduced_tree->Branch("mass",&fMass,"mass/F");
	reduced_tree->Branch("pt",&fPt,"pt/F");
	reduced_tree->Branch("eta",&fEta,"eta/F");
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("eta_mu1",&fEtaMu1,"eta_mu1/F");
	reduced_tree->Branch("eta_mu2",&fEtaMu2,"eta_mu2/F");
	reduced_tree->Branch("doca",&fDoca,"doca/F");
	reduced_tree->Branch("doca_z",&fDocaZ,"doca_z/F");
	reduced_tree->Branch("doca_xy",&fDocaXY,"doca_xy/F");
	reduced_tree->Branch("d3",&fD3,"d3/F");
	reduced_tree->Branch("d3_true",&fD3Truth,"d3_true/F");
	reduced_tree->Branch("alpha",&fAlpha,"alpha/F");
	reduced_tree->Branch("iso",&fIso,"iso/F");
	reduced_tree->Branch("pv_z",&fPvZ,"pv_z/F"),
	reduced_tree->Branch("pv_xy",&fPvXY,"pv_xy/F");
	reduced_tree->Branch("true_decay",&fTrueDecay,"true_decay/I");
	
	// dump track resolutions
	cout << "	Number of candidates: " << fNumCands << endl;
	cout << "	d0 resolution: " << fD0Resolution << " um" << endl;
	cout << "	dz resolution: " << fDzResolution << " um" << endl;
	cout << "	phi resolution: " << fPhiResolution << " x 10^(-3)" << endl;
	cout << "	cot(theta) resolution: " << fCotThetaResolution << " x 10^(-3)" << endl;
	cout << "	Pt resolution: " << fPtResolution << " %" << endl;
	// dump pv resolutions
	cout << "	XY PV Resolution: " << fPVResolutionXY << " um" << endl;
	cout << "	Z  PV Resolution: " << fPVResolutionZ << " um" << endl;
	// dump resolution mode
	cout << "	Resolution Mode: " << fResMode << endl;
	
	if (fResMode != kResFile) {
		
		switch (fResMode) {
			case kResCMSStd:
				pnt = &res.std_geom;
				break;
			case kResCMSUpg:
				pnt = &res.upg_geom;
				break;
			case kResCMSLong:
				pnt = &res.lng_geom;
				break;
			default:
				std::cerr << "ERROR: Selected resolution mode not yet supported!!!" << endl;
				break;
		}
		for (j = 0; j < NBR_ETA_RESOLUTION; j++) {
			cout << Form("\t\t===== eta %f=====",fEtaBins[j]) << endl;
			for (k = 0; k < fResVecXY[j].size(); k++) {
				res = fResVecXY[j][k];
				cout << Form("\t\tp = %f, sigma(xy) = %f", res.p, *pnt) << endl;
			}
			for (k = 0; k < fResVecZ[j].size(); k++) {
				res = fResVecZ[j][k];
				cout << Form("\t\tp = %f, sigma(z) = %f", res.p, *pnt) << endl;
			}
		}
	}
} // bookHist()

void pixelReader::eventProcessing()
{
	int j,nc;
	unsigned k;
	
	nc = fpEvt->nCands();
	for (j = 0; j < nc; j++) {
		for (k = 0; k < fNumCands; k++) {
			if (loadCandidateVariables(fpEvt->getCand(j)))
				reduced_tree->Fill();
		}
	}
} // eventProcessing()

void pixelReader::closeHistFile()
{
	fpHistFile = reduced_tree->GetCurrentFile();
	treeReader01::closeHistFile();
} // closeHistFile()

void pixelReader::readCuts(TString filename, int dump)
{
	char buffer[1024];
	char name[1024];
	float value;
	FILE *file = fopen(filename.Data(), "r");
	
	if (!file)
		goto bail;
	
	// read the cuts file...
	while(fgets(buffer, sizeof(buffer), file) != NULL) {
		if (buffer[0] == '#')
			continue; // comment
		
		if (buffer[strlen(buffer)-1] == '\n')
			buffer[strlen(buffer)-1] = 0;
		
		if(sscanf(buffer, "%s%f\n", name, &value) < 2) {
			std::cerr << "Unable to parse line:" << std::endl;
			std::cerr << buffer << std::endl;
			continue;
		}
		
		if (strcmp(name, "TRK_D0") == 0)
			fD0Resolution = value;
		else if (strcmp(name, "TRK_DZ") == 0)
			fDzResolution = value;
		else if (strcmp(name, "TRK_PHI") == 0)
			fPhiResolution = value;
		else if (strcmp(name, "TRK_COTTHETA") == 0)
			fCotThetaResolution = value;
		else if (strcmp(name, "TRK_PT") == 0)
			fPtResolution = value;
		else if (strcmp(name, "PV_XY") == 0)
			fPVResolutionXY = value;
		else if (strcmp(name, "PV_Z") == 0)
			fPVResolutionZ = value;
		else if (strcmp(name, "NUM_CANDS") == 0)
			fNumCands = (unsigned)value;
		else if (strcmp(name, "RES_MODE") == 0)
			fResMode = (unsigned)value;
		else
			std::cerr << "readCuts(): Unknown variable '" << name << "'found!" << std::endl;
	}
	
	fclose(file);
	
	// read the resolution if set...
	readResolution();
bail:
	return;
} // readCuts()

bool pixelReader::loadCandidateVariables(TAnaCand *pCand)
{
	using std::cout; using std::endl;
	int j;
	TAnaTrack *trk;
	TGenCand *gen;
	bool first;
	
	TVector3 pVtx;
	
	// poca berechnung
	TVector3 q;
	TVector3 p1,p2;
	TVector3 v1,v2;
	TVector3 a1,a2;
	
	clearVariables();
	
	// get simulated primary vertex
	for (j = 0; j < fpEvt->nGenCands(); j++) {
		gen = fpEvt->getGenCand(j);
		if (gen->fMom1 == 0) {
			pVtx = gen->fV;
			break;
		}
	}
	
	// set variables
	fMass = pCand->fMass;
	fPt = pCand->fPlab.Perp();
	fEta = pCand->fPlab.Eta();
	
	first = true;
	for (j = pCand->fSig1; j <= pCand->fSig2; j++) {
		trk = fpEvt->getSigTrack(j);
		gen = fpEvt->getGenCand(trk->fGenIndex);
		if (first) {
			p1 = gen->fP.Vect();
			v1 = gen->fV;
		} else {
			p2 = gen->fP.Vect();
			v2 = gen->fV;
		}

		first = false;
	}
	
	// Compute 3D distance w/o resolution
	q = p1.Cross(p2);
	a1 = v1 - p2.Cross(q).Dot(v1 - v2) / q.Mag2() * p1;
	q = p2.Cross(p1);
	a2 = v2 - p1.Cross(q).Dot(v2 - v1) / q.Mag2() * p2;
	// compute the flight length
	q = 0.5 * (a1 + a2); // center of pocas
	fD3Truth = (q - pVtx).Mag();
	
	// resolution of track
	smearTrack(&v1, &p1);
	fPtMu1 = p1.Perp();
	fEtaMu1 = p1.Eta();
	smearTrack(&v2, &p2);
	fPtMu2 = p2.Perp();
	fEtaMu2 = p2.Eta();
	// resolution of PV
	if (fPVResolutionXY > 0) {
		pVtx = smearPhi(pVtx,fPVResolutionXY*kMuMToCM);
		pVtx = smearR(pVtx,fPVResolutionXY*kMuMToCM);
	}
	if (fPVResolutionZ)
		pVtx = smearZ(pVtx,fPVResolutionZ*kMuMToCM);
	
	// make muon1 the leading pt muon
	if (fPtMu1 < fPtMu2) {
		std::swap(fPtMu1,fPtMu2);
		std::swap(fEtaMu1,fEtaMu2);
		std::swap(v1,v2);
		std::swap(p1,p2);
	}
	
	// compute the point of closes approach on each trajectory
	q = p1.Cross(p2);
	a1 = v1 - p2.Cross(q).Dot(v1 - v2) / q.Mag2() * p1;
	q = p2.Cross(p1);
	a2 = v2 - p1.Cross(q).Dot(v2 - v1) / q.Mag2() * p2;
	
	q = a2 - a1; // distance in 3d
	fDoca = q.Mag();
	fDocaZ = q.Z();
	q.SetZ(0); // project in xy plane
	fDocaXY = q.Mag();
	
	// compute the flight length
	q = 0.5 * (a1 + a2); // center of pocas
	fD3 = (q - pVtx).Mag();
	fAlpha = pCand->fPlab.Angle(q - pVtx);
	
	fIso = (float)compIso(pCand);
	
	// compute primary vertex position
	fPvZ = pVtx.Z();
	pVtx.SetZ(0);
	fPvXY = pVtx.Mag();
	pVtx.SetZ(fPvZ);
	
	// truth info
	fTrueDecay = loadTruth(pCand);
	
	return true;
} // loadCandidateVariables()

void pixelReader::clearVariables()
{
	fMass = -99.f;
	fPt = -99.f;
	fEta = -99.f;
	fPtMu1 = -99.f;
	fPtMu2 = -99.f;
	fEtaMu1 = -99.f;
	fEtaMu2 = -99.f;
	
	fDoca = 99.f;
	fDocaZ = -99.f;
	fDocaXY = -99.f;
	fD3 = -99.f;
	fAlpha = -99.f;
	
	fIso = -99.f;
	
	fPvZ = -9999.f;
	fPvXY = -9999.f;
	
	fTrueDecay = 0;
} // clearVariables()

double pixelReader::compIso(TAnaCand *pCand)
{
	TGenCand *gen;
	const double cone = 0.7;
	const double pt_thres = 0.9;
	double sum_pt = 0.0;
	double iso;
	int ix;
	
	// sum pt
	for (ix = 0; ix < fpEvt->nGenCands(); ix++) {
		gen = fpEvt->getGenCand(ix);
		if (gen->fQ == 0)
			continue; // only charged particles have a track
		if (fStableParticles.count(TMath::Abs(gen->fID)) == 0)
			continue; // only stable particles have a track
		if (gen->fP.Vect().Perp() <= pt_thres)
			continue; // fails pt threshold
		if (gen->fP.Vect().DeltaR(pCand->fPlab) >= cone)
			continue; // not in cone around candidate
		
		// sum it up
		sum_pt += gen->fP.Vect().Perp();
	}
	
	// compute the isolation of the candidate
	iso = pCand->fPlab.Perp();
	iso = iso / (iso + sum_pt);
	
	return iso;
} // compIso()

int pixelReader::loadTruth(TAnaCand *pCand)
{
	int j;
	TAnaTrack *trk;
	TGenCand *gen;
	TGenCand *theMom = NULL;
	int result = 0;
	decay_t dec;
	std::map<decay_t,int>::const_iterator it;
	
	for (j = pCand->fSig1; 0 <= j && j <= pCand->fSig2; j++) {
		trk = fpEvt->getSigTrack(j);
		gen = fpEvt->getGenCand(trk->fGenIndex);
		
		// look for the mother of the particle
		while (TMath::Abs(gen->fID) != 531) {
			if (gen->fMom1 < 0 || gen->fMom1 >= fpEvt->nGenCands())
				goto bail;
			gen = fpEvt->getGenCand(gen->fMom1);
		}
		
		if (theMom) {
			if (theMom->fNumber != gen->fNumber) goto bail; // not the same
		} else {
			theMom = gen;
		}
	}
	
	// still here, then we found the mother
	buildDecay(theMom, &dec);
	
	it = fDecayTable.find(dec);
	if (it != fDecayTable.end())
		result = it->second;
	
bail:
	return result;
} // loadTruth()

void pixelReader::buildDecay(TGenCand *pGen, decay_t *dec)
{
	dec->insert(abs(pGen->fID));
	for (int j = pGen->fDau1; 0 <= j && j <= pGen->fDau2; j++)
		buildDecay(fpEvt->getGenCand(j), dec);
} // buildDecay()

void pixelReader::dumpGenerator()
{
	using std::cout; using std::endl;
	int j;
	TGenCand *gen;
	
	cout << "================================" << endl;
	for (j = 0; j < fpEvt->nGenCands(); j++) {
		gen = fpEvt->getGenCand(j);
		cout << j << ": " << gen->fID << Form(" at (%f,%f,%f)", gen->fV.X(), gen->fV.Y(), gen->fV.Z());
		cout << " with dau [" << gen->fDau1 << ", " << gen->fDau2 << "] and moms = [" << gen->fMom1 << ", " << gen->fMom2 << "]" << endl;
	}
	cout << "================================" << endl;
} // dumpGenerator()

TVector3 pixelReader::smearPhi(TVector3 v, double res_cm)
{
	TVector3 vec_phi = TVector3(-v.Y(),v.X(),0);
	TVector3 result = v + fRand.Gaus(0, res_cm)*vec_phi.Unit();
	
	return result;
} // smearXY()

TVector3 pixelReader::smearR(TVector3 v, double res_cm)
{
	TVector3 vec_r(v.X(), v.Y(), 0);
	TVector3 result = v + fRand.Gaus(0, res_cm)*vec_r.Unit();
	
	return result;
} // smearR()

TVector3 pixelReader::smearZ(TVector3 v, double res_cm)
{
	TVector3 result(v);
	
	result.SetZ(fRand.Gaus(v.Z(), res_cm));
	
	return result;
} // smearZ()

TVector3 pixelReader::smearD0(TVector3 v, double res_cm, TVector3 plab)
{
	TVector3 e(-plab.Y(),plab.X(),0.0); // direction to smear
	TVector3 result = v + fRand.Gaus(0, res_cm)*e.Unit();
	
	return result;
} // smearD0()

void pixelReader::smearTrack(TVector3 *v, TVector3 *p)
{
	double pt,cot_theta,phi;
	TVector3 save = *v;
	
	// smear impact parameters
	if (fResMode == kResFile) {
		*v = smearD0(*v,fD0Resolution*kMuMToCM, *p);
		*v = smearZ(*v, fDzResolution*kMuMToCM);
	} else {
		*v = smearD0(*v,findResolution(*p, fResVecXY), *p);
		*v = smearZ(*v, findResolution(*p, fResVecZ));
	}
	
	// smear momentum
	pt = p->Perp();
	cot_theta = p->Z()/TMath::Sqrt(p->X()*p->X() + p->Y()*p->Y());
	phi = p->Phi();
	
	pt = fRand.Gaus(pt, fPtResolution/100. * pt);
	cot_theta = fRand.Gaus(cot_theta, fCotThetaResolution);
	phi = fRand.Gaus(phi, fPhiResolution);
	
	// build momentum vector again
	cot_theta = TMath::Pi()/2.0 - TMath::ATan(cot_theta); // now it is theta
	pt = pt / TMath::Sin(cot_theta); // now it is p
	p->SetMagThetaPhi(pt,cot_theta,phi);
} // smearTrack()

void pixelReader::readResolution()
{
	int j;
	
	// parse xy file
	parseResFile(fResVecXY, "resolution/ip_xy_pu50_dloss.txt");
	parseResFile(fResVecZ, "resolution/ip_z_pu50_dloss.txt");
	
	// sort the vectors
	for (j = 0; j < NBR_ETA_RESOLUTION; j++) {
		std::sort(fResVecXY[j].begin(), fResVecXY[j].end());
		std::sort(fResVecZ[j].begin(), fResVecZ[j].end());
	}
} // readResolution()

void pixelReader::parseResFile(std::vector<res_t> *resVec, const char *resFileName)
{
	using namespace std;
	FILE *resFile = fopen(resFileName, "r");
	char buffer[1024];
	double lastP;
	res_t res;
	int parsed;
	
	if (!resFile)
		goto bail;
	
	lastP = -1;
	
	while(fgets(buffer, sizeof(buffer), resFile) != NULL) {
		if (buffer[0] == '#') // comment
			continue;
		
		if (buffer[0] == '\n' || buffer[0] == '\r') // empty line
			continue;
		
		if( (parsed = sscanf(buffer, "%lf %lf %lf %lf",&res.p, &res.std_geom, &res.upg_geom, &res.lng_geom)) < 3) {
			cerr << "ERROR: parseResFile(). parsed = " << parsed << " < 3" << endl;
			continue;
		}
		if (parsed == 3) // no entry for long geometry, use standard (i.e. xy)
			res.lng_geom = res.std_geom;
		
		if (res.p < lastP)
			resVec++;
		
		resVec->push_back(res);
		lastP = res.p;
	}
	
bail:
	if (resFile) fclose(resFile);
} // parseResFile()

double pixelReader::findResolution(TVector3 plab, std::vector<res_t> *resVec)
{
	unsigned j,k;
	double result = 0.0;
	double dist, minDist = 1.e30;
	res_t res;
	double *ptr;
	
	switch (fResMode) {
		case kResCMSStd:
			ptr = &res.std_geom;
			break;
		case kResCMSUpg:
			ptr = &res.upg_geom;
			break;
		case kResCMSLong:
			ptr = &res.lng_geom;
			break;
		default:
			ptr = NULL;
			break;
	}
	
	for (j = 0; j < NBR_ETA_RESOLUTION-1; j++) { // all eta > max => eta = max
		if (TMath::Abs(plab.Eta()) < fEtaBins[j])
			break;
	}
	
	for (k = 0; k < resVec[j].size(); k++) {
		dist = TMath::Abs(plab.Mag() - resVec[j][k].p);
		if (dist < minDist) {
			res = resVec[j][k];
			result = *ptr;
			minDist = dist;
		}
	}
	
	return result;
} // findResolution()
