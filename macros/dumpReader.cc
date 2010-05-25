#include "dumpReader.hh"

#include <set>
#include <stack>

static void dumpTabs(unsigned nbr)
{
	for (unsigned j = 0; j < nbr; j++) cout << '\t';
} // dumpTabs()

dumpReader::dumpReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName)
{
	// build the true decay channel
	true_channel.insert(13);
	true_channel.insert(13);
	true_channel.insert(211);
	true_channel.insert(211);
	true_channel.insert(310);
	true_channel.insert(443);
	true_channel.insert(511);
	
	unreco_counter = 0;
} // dumpReader()

dumpReader::~dumpReader()
{
	map<int,int>::const_iterator it;
	cout << "=====================" << endl;
	cout << "dumpReader destructor" << endl;
	
	for (it = evt_nbr.begin(); it != evt_nbr.end(); ++it) {
		cout << it->first << " of B0->J/Psi Ks decays in " << it->second << " events." << endl;
	}
	
	for (it = rec_nbr.begin(); it != rec_nbr.end(); ++it) {
		cout << it->first << " tracks reconstructed in " << it->second << " events." << endl;
	}
	
	for (it = muon_nbr.begin(); it != muon_nbr.end(); ++it) {
		cout << it->first << " muons out of 4 in " << it->second << " events." << endl;
	}
	
	cout << "Unrecovered Decays counted: " << unreco_counter << endl;
} // ~dumpreader()

void dumpReader::eventProcessing()
{
	unsigned j,k,nc;
	TGenCand *pGen;
	unsigned nbrTrack,nbrMuons;
	multiset<int> particles;
	set<int> recoverableGens;
	set<int>::iterator it;
	TAnaCand *pCand;
	static unsigned evt_counter = 0;
	bool muCut = true;
	bool piCut = true;
	
	decay_counter.clear();
	
	nc = fpEvt->nGenCands();
	for (j = 0;  j < nc; j++) {
		
		pGen = fpEvt->getGenCand(j);
		if(abs(pGen->fID) == 511) {
			
			particles.clear();
			nbrMuons = 0;
			
			nbrTrack = checkReconstruction(j,&particles,&nbrMuons);
			
			// remove the photons (Bremsstrahlung)
			particles.erase(22);
			decay_counter[particles]++;
			
			if (particles == true_channel) {
				rec_nbr[nbrTrack]++;
				
				if(nbrTrack == 4) {
					muon_nbr[nbrMuons]++;
					
					// dump candidates
					if (nbrMuons > 0) {
						recoverableGens.insert(j);
						// dumpGenerator(pGen);
					}
				}
			}
		}
	}
	
	// now, remove all the decays we actually have recovered from recoverableGens
	nc = fpEvt->nCands();
	for (j = 0; j < nc; j++) {
		pCand = fpEvt->getCand(j);
		
		if (pCand->fType != 600511)
			continue;
		
		k = getGenIndex(pCand);
		recoverableGens.erase(k);
	}
	
	// remove all the candidates from the muon pt cut (< 1.0 GeV => out)
	if (muCut) {
		
		set<int> muonOut;
		set<int> tracks;
		set<int>::const_iterator trackIt;
		TAnaTrack *pTrack;
		
		for (it = recoverableGens.begin(); it != recoverableGens.end(); ++it) {
			
			tracks.clear();
			buildTracks(*it,&tracks,13);
			
			for (trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {
				
				pTrack = fpEvt->getRecTrack(*trackIt);
				if (pTrack->fPlab.Pt() < 1.0)
					muonOut.insert(*it);
			}
		}
		
		// remove all the muonOut candidates
		for (trackIt = muonOut.begin(); trackIt!=muonOut.end(); ++trackIt)
			recoverableGens.erase(*trackIt);
	}
	
	// remove all the candidates from the pion pt cut (< 0.4 GeV => out)
	if (piCut) {
		
		set<int> pionOut;
		set<int> tracks;
		set<int>::const_iterator trackIt;
		TAnaTrack *pTrack;
	
		for (it = recoverableGens.begin(); it != recoverableGens.end(); ++it) {
			
			tracks.clear();
			buildTracks(*it,&tracks,211);
			
			for (trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {
				
				pTrack = fpEvt->getRecTrack(*trackIt);
				if (pTrack->fPlab.Pt() < 0.4)
					pionOut.insert(*it);
			}
		}
		
		// remove all the muonOut candidates
		for (trackIt = pionOut.begin(); trackIt!=pionOut.end(); ++trackIt)
			recoverableGens.erase(*trackIt);
	}
	
	
	if (recoverableGens.size() > 0) {
		unreco_counter++; // counted one more
		
		/*if (unreco_counter == 1)*/ {
			// Print info on the tracks!!!
			cout << "Event Number: " << evt_counter << endl;
			for (it = recoverableGens.begin(); it != recoverableGens.end(); ++it)
				dumpUnrecovered(fpEvt->getGenCand(*it));
		}
	}
	
	evt_nbr[decay_counter[true_channel]]++;
	evt_counter++;
} // eventProcessing()

int dumpReader::getGenIndex(TAnaCand *pCand)
{
	TAnaTrack *pTrack;
	TGenCand *momCand,*pGen;;
	int j;
	
	// set the momCand out of the first signal track...
	pTrack = fpEvt->getSigTrack(pCand->fSig1);
	pTrack = fpEvt->getRecTrack(pTrack->fIndex);
	if (pTrack->fGenIndex < 0) return -1;
	momCand = fpEvt->getGenCand(pTrack->fGenIndex);
	
	while (abs(momCand->fID) != 511 && momCand->fMom1 >= 0)
		momCand = fpEvt->getGenCand(momCand->fMom1);
	
	if (abs(momCand->fID) != 511) return -1;
	
	for (j = pCand->fSig1+1; j <= pCand->fSig2; j++) {
		
		pTrack = fpEvt->getSigTrack(j);
		pTrack = fpEvt->getRecTrack(pTrack->fIndex);
		
		if (pTrack->fGenIndex < 0) return -1;
		pGen = fpEvt->getGenCand(pTrack->fGenIndex);
		
		while(abs(pGen->fID) != 511 && pGen->fMom1 >= 0)
			pGen = fpEvt->getGenCand(pGen->fMom1);
		
		if((abs(pGen->fID)!=511) || (pGen->fNumber != momCand->fNumber)) return -1;
	}
	
	return momCand->fNumber;
} // getGenIndex()

void dumpReader::dumpCandidate(TAnaCand *pCand, unsigned indent)
{
	int j;
	
	dumpTabs(indent);
	cout << pCand->fType << '{' << endl;
	
	if (pCand->fDau1 >= 0 && pCand->fDau2 >= 0) {
		for (j = pCand->fDau1; j <= pCand->fDau2; j++) {
			dumpCandidate(fpEvt->getCand(j), indent+1);
		}
	}
	dumpTabs(indent);
	cout << '}' << endl;
} // dumpCandidate()

int dumpReader::checkReconstruction(int genIx, multiset<int> *particles, unsigned *nbrMuons)
{
	TAnaTrack *pTrack;
	TGenCand *pGen;
	int j,result = 0;
	
	// check the daughters
	pGen = fpEvt->getGenCand(genIx);
	
	// add this particle to the particles
	particles->insert(abs(pGen->fID));
	
	// check the daughters
	for (j = pGen->fDau1; j <= pGen->fDau2 && j >= 0; j++) {
		result += checkReconstruction(j,particles,nbrMuons);
	}
	
	// check the tracks
	for (j = 0; j < fpEvt->nRecTracks(); j++) {
		pTrack = fpEvt->getRecTrack(j);
		if (genIx == pTrack->fGenIndex) {
			result++;
			
			if(abs(pGen->fID) == 13 && pTrack->fMuID > 0 && (pTrack->fMuID & 6))
				(*nbrMuons)++;
			
			// FIXME: HIER ÜBERPRÜFEN
			break;
		}
	}
	
	return result;
} // checkReconstruction()

void dumpReader::dumpGenerator(TGenCand *pGen)
{
	set<int> trackIndices;
	set<int>::const_iterator it;
	
	// get all the track indices...
	buildTracks(pGen->fNumber,&trackIndices);
	
	cout << "Tracks for Generator Candidate " << pGen->fNumber << endl;
	for (it = trackIndices.begin(); it != trackIndices.end(); ++it) {
		cout << '\t' << *it;
	}
	cout << endl << endl;
	
} // dumpGenerator()

void dumpReader::buildTracks(int genIx, set<int> *trackIndices, int ID)
{
	TGenCand *pGen = fpEvt->getGenCand(genIx);
	TAnaTrack *pTrack;
	int j;
	
	// fill the daughter tracks
	for (j = pGen->fDau1; j <= pGen->fDau2; j++)
		buildTracks(j,trackIndices, ID);
	
	// fill the tracks pointing to this particle
	for (j = 0; j < fpEvt->nRecTracks() && abs(pGen->fID) == ID; j++) {
		pTrack = fpEvt->getRecTrack(j);
		if (genIx == pTrack->fGenIndex) {
			trackIndices->insert(j);
			break;
		}
	}
} // buildTracks()

void dumpReader::dumpUnrecovered(TGenCand *pGen)
{
	stack<int> indices;
	int j,ix,nc;
	TAnaTrack *pTrack;
	
	indices.push(pGen->fNumber);
	
	cout << "Dumping unrecovered Generator Candidate:" << endl;
	
	while (!indices.empty()) {
		
		// get the top-most element
		ix = indices.top();
		indices.pop();
		
		pGen = fpEvt->getGenCand(ix);
		for (j = pGen->fDau1; j <= pGen->fDau2; j++)
			indices.push(j);
		
		// try to find matching tracks!
		nc = fpEvt->nRecTracks();
		for (j = 0; j < nc; j++) {
			pTrack = fpEvt->getRecTrack(j);
			if (pTrack->fGenIndex == pGen->fNumber) {
				// dump the track!!!
				cout << "RecIndex = " << j << ", MCID = " << pTrack->fMCID << ", GenIndex = " << pTrack->fGenIndex
					<< ", (x,y,z) = " << '(' << pTrack->fPlab.x() << ',' << pTrack->fPlab.y() << ',' << pTrack->fPlab.z()
					<< "), MuID = " << pTrack->fMuID << ", Pt = " << pTrack->fPlab.Pt() << endl;
			}
		}
	}
	
	cout << "===================" << endl;
} // dumpUnrecovered()
