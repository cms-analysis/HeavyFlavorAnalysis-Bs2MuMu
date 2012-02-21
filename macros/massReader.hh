#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"
#include <cmath>
#include <set>
#include <map>
#include <stdlib.h>
#include <stdint.h>

// we store the first four tracks  
#define NBR_TRACKS_STORE 4

const static double MMUON = 0.1057;
const static double MPION = 0.1396;
const static double MKAON = 0.4937;

// The trigger information stored.
// When adding one, be sure to update the code in massReader::loadTrigger()
enum trigger_bits
{
	kHLT_DoubleMu3_Bit						= 1 << 0,
	kHLT_DoubleMu0_Bit						= 1 << 1,
	kHLT_DoubleMu0_Quarkonium_Bit			= 1 << 2,
	kHLT_DoubleMu3_Jpsi_Bit					= 1 << 3,
	kHLT_DoubleMu3_Bs_Bit					= 1 << 4,
	kHLT_DoubleMu2_Bs_Bit					= 1 << 5,
	kHLT_Dimuon6p5_Jpsi_Displaced_Bit		= 1 << 6,
	kHLT_Dimuon7_Jpsi_Displaced_Bit			= 1 << 7,
	kHLT_Dimuon6_Bs_Bit						= 1 << 8,
	kHLT_Dimuon4_Bs_Barrel_Bit				= 1 << 9,
	kHLT_DoubleMu4_Dimuon6_Bs_Bit			= 1 << 10,
	kHLT_DoubleMu4_Dimuon4_Bs_Barrel_Bit	= 1 << 11,
	kHLT_DoubleMu3p5_Jpsi_Displaced_Bit		= 1 << 12,
	kHLT_DoubleMu4_Jpsi_Displaced_Bit		= 1 << 13
};

struct trigger_table_t {
	trigger_bits t_bit;
	std::string trigger_name;
	std::pair<int64_t,int64_t> run_range;
}; // trigger_table_t

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
		float calculateIsolation(TAnaCand *pCand);
		float calculateDoca0(TAnaCand *pCand);
		int countTracksNearby(TAnaCand *pCand);
		int loadTrigger(int *errTriggerOut = NULL, int *triggersFoundOut = NULL);
		int isMuonTight(TAnaTrack *sigTrack);
		int hasTriggeredNorm();
		int hasTriggeredSignal();
		int hasTriggered(int triggers,trigger_table_t *table, unsigned size);
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
		float fPMu1_Gen; // p of generator muon 1
		float fPMu2_Gen; // p of generator muon 2
		float fPtMu1_Gen; // pt of generator muon
		float fPtMu2_Gen; // pt of generator muon
		float fEtaMu1_Gen; // eta of gen muon 1
		float fEtaMu2_Gen; // eta of gen muon 2
		int fTrackQual_mu1; // track quality of muon 1
		int fTrackQual_mu2; // track quality of muon 2
		int fQ_mu1; // charge of muon 1
		int fQ_mu2;	// charge of muon 2
		float fDeltaPhiMu; // mu1.Phi(mu2)
		float fIsoMoriond12; // isolation variable defined as for moriond 12
		float fDoca0;
		int fNbrNearby;
		// muon properties
		float fPtMu1,fPtMu2;
		int fMuID1,fMuID2;
		int fMuTight1,fMuTight2; // muon is identified as tight muond
		float fEtaMu1,fEtaMu2;
		float fDeltaR; // deltaR of the muons
		// triggers
		int	fTriggers; // store some trigger information
		int fTriggersError; // error information of trigger
		int fTriggersFound; // what triggers were available
		int fTriggeredJPsi;
		int fTriggeredBs;
		float fCtau;
		float fCtauE;
		float fEta; // eta of the candidate
		int fNbrPV; // nbr of PV in this event
		float fIPCand;
		float fIPCandE;
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
		int fCutTrackQual_mu1;
		int fCutTrackQual_mu2;
		int fCutMuID_mask;
		bool fCutMuID_reqall;
		bool fCutOppSign_mu;
};

#endif
