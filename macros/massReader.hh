#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"
#include <set>
#include <map>

// we store the first four tracks  
#define NBR_TRACKS_STORE 4

const static double MMUON = 0.1057;
const static double MPION = 0.1396;
const static double MKAON = 0.4937;

// The trigger information stored.
// When adding one, be sure to update the code in massReader::loadTrigger()
enum trigger_bits
{
	kHLT_DoubleMu3_Bit = 1 << 0,
	kHLT_DoubleMu0_Bit = 1 << 1,
	kHLT_DoubleMu0_Quarkonium_v1_Bit = 1 << 2
};

// Truth Flags Bits
enum truth_bits
{
	kTruthSameB_Bit = 1 << 0
};

// Efficiency Flags
enum
{
	kGeneratorCand	= 1 << 0,
	kAcceptance		= 1 << 1,
	kEffMuon		= 1 << 2,
//	kEffTrig is not needed
//	kEffCand is not needed
//	kEffAna is not needed
};

class massReader : public treeReader01 {
	
	public:
		massReader(TChain *tree, TString evtClassName);
		~massReader();

		virtual void bookHist();
		virtual void eventProcessing();
		virtual void closeHistFile();
		virtual void readCuts(TString filename, int dump = 1);

	protected:
		// For subclasses
		TTree *reduced_tree;
		std::multiset<int> trueDecay;
		std::set<int> stableParticles; // the particles in here are considered stable
		int fTruthType;
		
		virtual void clearVariables();
		virtual int loadGeneratorVariables(TGenCand *pGen);
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
		
		// loading variables
		virtual int checkTruth(TAnaCand *cand); // check if all are originating from the same particle
		virtual int loadTruthFlags(TAnaCand *cand); // set truth flags
		virtual int countMuons(TAnaCand *cand); // count the number of identified muons
		float calculateIsolation(TAnaCand *pCand, double openingAngle, double minPt, bool sameVertex);
		int loadTrigger(int *errTriggerOut = NULL, int *triggersFoundOut = NULL);
		virtual int loadEfficiencyFlags(TGenCand *gen);

		// other utility routines
		TAnaCand *findCandidate(int candID, std::map<int,int> *particles); // particles = map(recTrackIx,particleID)
		void buildDecay(TGenCand *gen, std::multiset<int> *particles);
		void findCandStructure(TAnaCand* pCand, std::map<int,int> *particles); // map(recTrackIx,particleID)
		void findGenStructure(TGenCand *pGen, std::map<int,int> *particles); // map(genIx,recTrackIx)
		void findAllTrackIndices(TAnaCand* pCand, std::map<int,int> *indices); // map(recTrackIx,sigTrackIndex)
	protected:
		const char *fTreeName;
		
	protected:
		int fCandidate;
		float fMass;
		float fMassConstraint; // mass of the constraint candidate...
		int fTruth;		// is this background or a true candidate?
		int fTruthFlags; // several partial truth flags
		int fEffFlags; // flags to calculate the efficiencies...
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
		float fPtMu1_Gen; // pt of generator muon
		float fPtMu2_Gen; // pt of generator muon
		// isolation variables. fIsoX_ptY means opening angle deltaR < X/10 and only sum over
		// tracks with pt > Y/10 GeV
		float fIso7_pt0;
		float fIso7_pt5;
		float fIso7_pt7;
		float fIso7_pt10;
		float fIso10_pt0;
		float fIso10_pt5;
		float fIso10_pt7;
		float fIso10_pt9; // default isolation of Analysis node
		float fIso10_pt10;
		float fIso7_pt0_pv;
		float fIso7_pt5_pv;
		float fIso7_pt7_pv;
		float fIso7_pt10_pv;
		float fIso10_pt0_pv;
		float fIso10_pt5_pv;
		float fIso10_pt7_pv;
		float fIso10_pt9_pv;
		float fIso10_pt10_pv;
		int fTriggers; // store some trigger information
		int fTriggersError; // error information of trigger
		int fTriggersFound; // what triggers were available
		float fCtau; // proper time (note can be filled only in subclasses as requires knowledge of m)
		float fEta; // eta of the candidate
		int fTracksIx[NBR_TRACKS_STORE];
		float fTracksIP[NBR_TRACKS_STORE];
		float fTracksIPE[NBR_TRACKS_STORE];
		float fTracksPT[NBR_TRACKS_STORE];
		float fTracksPTRel[NBR_TRACKS_STORE];
	
	// Cut variables
	protected:
		bool fCutFileParsed;
		bool fCutTriggered;
		int fCutCand; // the candidate to extract
		double fCutFlight3dSign;
		double fCutChi2;
		double fCutPt;
		double fCutAlpha;
		double fCutChi2ByNdof;
		int fCutTruth; // truth matching
};

#endif
