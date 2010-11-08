#include "impactReader.hh"

#include <set>
#include <stack>
#include <cmath>

using std::cout;
using std::endl;

impactReader::impactReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName)
{} // impactReader()

void impactReader::bookHist()
{
	treeReader01::bookHist();
	
	reduced_tree = new TTree("T","");
	
	// add the branches
	reduced_tree->Branch("lip_cmssw",&fLip_CMSSW,"lip_cmssw/F");
	reduced_tree->Branch("tip_cmssw",&fTip_CMSSW,"tip_cmssw/F");
	reduced_tree->Branch("lip_geom",&fLip_geom,"lip_geom/F");
	reduced_tree->Branch("tip_geom",&fTip_geom,"tip_geom/F");
	reduced_tree->Branch("ix_cmssw",&fIx_CMSSW,"ix_cmssw/I");
	reduced_tree->Branch("ix_geom",&fIx_geom,"ix_geom/I");
	reduced_tree->Branch("nbr_pv",&fNbrPV,"nbr_pv/I");
	reduced_tree->Branch("event_nbr",&fEvent,"event_nbr/I");
} // bookHist()

void impactReader::closeHistFile()
{
	fpHistFile = reduced_tree->GetCurrentFile();
	treeReader01::closeHistFile();
} // closeHistFile()

void impactReader::eventProcessing()
{
	int j;
	TAnaCand *pCand;
	TVector3 d3;
	TVector3 minDist;
	
	for (j = 0; j < fpEvt->nCands(); j++) {
		
		pCand = fpEvt->getCand(j);
		if (pCand->fType != 300521)
			continue;
		
		// set variables
		fIx_CMSSW = pCand->fPvIdx;
		fLip_CMSSW = pCand->fPvLip;
		fTip_CMSSW = pCand->fPvTip;
		fIx_geom = calculatePVIx(pCand);
		fLip_geom = calculatePVDist(pCand,fIx_geom).Z();
		fTip_geom = calculatePVDist(pCand, fIx_geom).Perp();
		fNbrPV = fpEvt->nPV();
		
		// fill the tree
		reduced_tree->Fill();
	}
} // eventProcessing()

int impactReader::calculatePVIx(TAnaCand *pCand)
{
	int ix = -1;
	int j;
	float minLip = -1.0,lip;
	
	for (j = 0; j < fpEvt->nPV(); j++) {
		
		lip = calculatePVDist(pCand, j).Z();
		if (ix < 0 || fabs(lip) < fabs(minLip)) {
			minLip = lip;
			ix = j;
		}
	}
	
	return ix;
} //  calculatePVIx()

TVector3 impactReader::calculatePVDist(TAnaCand *pCand, int pvIx)
{
	TVector3 d3 = pCand->fVtx.fPoint - fpEvt->getPV(pvIx)->fPoint;
	TVector3 result;
	double scale;
	
	scale = (d3 * pCand->fPlab)/pCand->fPlab.Mag2();
	result = d3 - scale * pCand->fPlab;
	
	return result;
} // calculatePVDist()
