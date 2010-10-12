#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"
#include <set>

const static double MMUON = 0.1057;
const static double MPION = 0.1396;
const static double MKAON = 0.4937;

// The trigger information stored.
// When adding one, be sure to update the code in massReader::loadTrigger()
enum trigger_bits
{
	kHLT_DoubleMu3_Bit = 1 << 0,
	kHLT_DoubleMu0_Bit = 1 << 1
};

class massReader : public treeReader01 {
	
	public:
		massReader(TChain *tree, TString evtClassName);
		~massReader();

		virtual void bookHist();
		virtual void eventProcessing();
		virtual void closeHistFile();

	protected:
		// For subclasses
		TTree *reduced_tree;
		std::multiset<int> trueDecay;
		
		virtual int loadCandidateVariables(TAnaCand *pCand);
		
		// loading variables
		virtual int checkTruth(TAnaCand *cand); // check if all are originating from the same particle
		virtual int countMuons(TAnaCand *cand); // count the number of identified muons
		float calculateIsolation(TAnaCand *pCand, double openingAngle); // calculate the isolation of the candidate
		int loadTrigger();

		// creates the decay of the TGenCand
		void buildDecay(TGenCand *gen, std::multiset<int> *particles);
	
	protected:
		const char *fTreeName;
		
	protected:
		int fCandidate;
		TVector3 fMomentum;
		TVector3 fPVPosition;
		TVector3 fCandVertex;
		TVector3 *fMomentumPtr;
		TVector3 *fPVPositionPtr;
		TVector3 *fCandVertexPtr;
		float fMass;
		int fTruth;		// is this background or a true candidate?
		float fPt; // pt of the top particle
		float fNbrMuons;  // number of muons in the muon list.
		float fD3;
		float fD3E;
		float fDxy;	// distance to originating vertex (2d)
		float fDxyE;	// error (2d)
		float fAlpha; // angle between momentum and dist(vertex, motherVertex)
		float fAlphaXY; // Angle between momentum and dist(vertex, motherVertex) in xy-plane
		float fChi2; // chi2 of the vertex
		float fNdof; // number of degrees of freedom of vertex
		float fMaxDoca; // max doca
		float fIso; // isolation
		int fTriggers; // store some trigger information
		float fCtau; // proper time (note can be filled only in subclasses as requires knowledge of m)
		float fD3_Perp;	// Perpendicular part of distance d3 w.r.t. momentum of candidate
		float fD3_Para; // Parallel part of distance d3 w.r.t. momentum of candidate
		float fDxy_Perp; // Perpendicular part of distance dxy w.r.t. momentum of candidate
		float fDxy_Para; // Parallel part of distance dxy w.r.t. momentum of candidate
};

#endif
