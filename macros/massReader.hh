#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"
#include <cmath>
#include <set>
#include <map>
#include <stdlib.h>
#include <stdint.h>

const static double MMUON = 0.1057;
const static double MPION = 0.1396;
const static double MKAON = 0.4937;

// enumeration of known decays for truth matching
// corresponds to HFTruthCandidates_cff.py
enum {
	kDecay_BsToMuMu			= 1,
	kDecay_BsToMuMuGa		= 2,
	kDecay_BsToKK			= 3,
	kDecay_BsToKPi			= 4,
	kDecay_BsToPiPi			= 5,
	kDecay_BsToPiMuNu		= 6,
	kDecay_BsToKMuNu		= 7,
	kDecay_BdToMuMu			= 8,
	kDecay_BdToPiPi			= 9,
	kDecay_BdToKPi			= 10,
	kDecay_BdToKK			= 11,
	kDecay_BdToMuMuPi0		= 12,
	kDecay_BdToPiMuNu		= 13,
	kDecay_BuTo3MuNu		= 14,
	kDecay_LambdaBToPPi		= 15,
	kDecay_LambdaBToPK		= 16,
	kDecay_LambdaBToPMuNu	= 17,
	kDecay_Bs2JpsiPhi		= 18,
	kDecay_Bu2JpsiKp		= 19,
	kDecay_Bd2JpsiKstar		= 20,
	kDecay_Bd2JpsiKs		= 21,
	kDecay_PsiToMuMu		= 22,
	kDecay_Psi2SToMuMu		= 23,
	kDecay_Ups1SToMuMu		= 24,
	kDecay_Ups2SToMuMu		= 25,
	kDecay_Ups3SToMuMu		= 26
};

// used to store decays...
typedef std::multiset<int> decay_t;

// The trigger information stored.
// When adding one, be sure to update the code in massReader::loadTrigger()
enum trigger_bits
{
	kHLT_DoubleMu3_Bit							= 1 << 0,
	kHLT_DoubleMu0_Bit							= 1 << 1,
	kHLT_DoubleMu0_Quarkonium_Bit				= 1 << 2,
	kHLT_DoubleMu3_Jpsi_Bit						= 1 << 3,
	kHLT_DoubleMu3_Bs_Bit						= 1 << 4,
	kHLT_DoubleMu2_Bs_Bit						= 1 << 5,
	kHLT_Dimuon6p5_Jpsi_Displaced_Bit			= 1 << 6,
	kHLT_Dimuon7_Jpsi_Displaced_Bit				= 1 << 7,
	kHLT_Dimuon6_Bs_Bit							= 1 << 8,
	kHLT_Dimuon4_Bs_Barrel_Bit					= 1 << 9,
	kHLT_DoubleMu4_Dimuon6_Bs_Bit				= 1 << 10,
	kHLT_DoubleMu4_Dimuon4_Bs_Barrel_Bit		= 1 << 11,
	kHLT_DoubleMu3p5_Jpsi_Displaced_Bit			= 1 << 12,
	kHLT_DoubleMu4_Jpsi_Displaced_Bit			= 1 << 13,
	kHLT_DoubleMu3_4_Dimuon5_Bs_Central_Bit		= 1 << 14,
	kHLT_DoubleMu3p5_4_Dimuon5_Bs_Central_Bit	= 1 << 15,
	kHLT_DoubleMu4_Dimuon7_Bs_Forward_Bit		= 1 << 16
};

struct trigger_table_t {
	trigger_bits t_bit;
	std::string trigger_name;
	std::pair<int64_t,int64_t> run_range;
}; // trigger_table_t

class massReader : public treeReader01 {
	
	public:
		massReader(TChain *tree, TString evtClassName);
		virtual ~massReader();

		virtual void bookHist();
		virtual void eventProcessing();
		virtual void closeHistFile();
		virtual void readCuts(TString filename, int dump = 1);

	protected:
		// For subclasses
		TTree *reduced_tree;
		std::map<decay_t,int> decayTable;
		std::set<int> stableParticles; // the particles in here are considered stable
		std::set<int> validMothers;
		
		virtual void clearVariables();
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
		
		// loading variables
		int sameMother(TAnaCand *cand); // check if all are originating from the same mother
		int loadDecay(TAnaCand *anaCand);
		float calculateIsolation(TAnaCand *pCand);
		float calculateMuonIsolation(TAnaCand *pCand, int muonIx);
		float calculateDoca0(TAnaCand *pCand);
		int countTracksNearby(TAnaCand *pCand);
		int loadTrigger(int *errTriggerOut = NULL, int *triggersFoundOut = NULL);
		int isMuonTight(TAnaTrack *sigTrack);
		int hasTriggeredNorm();
		int hasTriggeredSignal(bool barrel);
		int hasTriggered(int triggers,trigger_table_t *table, unsigned size);

		// other utility routines
		TAnaCand *findCandidate(int candID, std::map<int,int> *particles); // particles = map(recTrackIx,particleID)
		void buildDecay(TGenCand *gen, std::multiset<int> *particles);
		void findCandStructure(TAnaCand* pCand, std::map<int,int> *particles); // map(recTrackIx,particleID)
		void findGenStructure(TGenCand *pGen, std::map<int,int> *particles); // map(genIx,recTrackIx)
		void findAllTrackIndices(TAnaCand* pCand, std::map<int,int> *indices); // map(recTrackIx,sigTrackIndex)
	protected:
		const char *fTreeName;
		
	protected:
		// Candidate Variables
		int fCandidate;
		float fMass;
		float fMassConstraint; // mass of the constraint candidate...
		float fPt; // pt of the top particle
		int fSameMother; // set if the originate from the same 'valid' mother
		int fTrueDecay; // set the MC true decay
		float fMaxDoca; // max doca
		float fEta; // eta of the candidate
		float fIPCand;
		float fIPCandE;
		// SV Vertex Variables
		float fD3;
		float fD3E;
		float fCtau;
		float fCtauE;
		float fDxy;	// distance to originating vertex (2d)
		float fDxyE;	// error (2d)
		float fAlpha; // angle between momentum and dist(vertex, motherVertex)
		float fAlphaXY; // Angle between momentum and dist(vertex, motherVertex) in xy-plane
		float fChi2; // chi2 of the vertex
		float fNdof; // number of degrees of freedom of vertex
		// Isolation variables
		float fIsoMoriond12; // isolation variable defined as for moriond 12
		float fDoca0;
		float fNbrNearby;
		// muon properties
		float fPtMu1,fPtMu2;
		int fMuTight1,fMuTight2;
		float fEtaMu1,fEtaMu2;
		float fDeltaPhiMu; // mu1.Phi(mu2)
		float fDeltaR; // deltaR of the muons
		float fPMu1_Gen; // p of generator muon 1
		float fPMu2_Gen; // p of generator muon 2
		float fPtMu1_Gen; // pt of generator muon
		float fPtMu2_Gen; // pt of generator muon
		float fEtaMu1_Gen; // eta of gen muon 1
		float fEtaMu2_Gen; // eta of gen muon 2
		int fHighPur_mu1; // track 1 is high purity
		int fHighPur_mu2; // track 2 is high purity
		int fQ_mu1; // charge of muon 1
		int fQ_mu2;	// charge of muon 2	
		// J/Psi
		float fPtJPsi;
		float fMassJPsi;
		float fChi2Jpsi;
		// phi
		float fMassPhi;
		// kaons
		float fPtKp1;
		float fPtKp2;
		float fEtaKp1;
		float fEtaKp2;
		float fPtKp_Gen1;
		float fPtKp_Gen2;
		float fEtaKp_Gen1;
		float fEtaKp_Gen2;
		int fHighPur_kp1;
		int fHighPur_kp2;
		int fQ_kp1;
		int fQ_kp2;
		float fDeltaR_Kaons; // delta R of the two Kaons
		// triggers
		int	fTriggers; // store some trigger information
		int fTriggersError; // error information of trigger
		int fTriggersFound; // what triggers were available
		int fTriggeredJPsi;
		int fTriggeredBsBarrel;
		int fTriggeredBsEndcap;
		
		// PV
		int fNbrPV; // nbr of PV in this event
		float fPVz;
		float fPVTrkWeight;
		
		// muon isolation
		float fIsoMu1;
		float fIsoMu2;
		
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
		bool fCutMuID_reqall;
		bool fCutOppSign_mu;
		// Additional cut variables for j/psi
		double fCutMass_JPsiLow;
		double fCutMass_JPsiHigh;
		// Additional cut variables for kaon
		double fCutPt_Kaon1;
		double fCutPt_Kaon2;
		bool fCutOppSign_kp;
		double fCutMass_PhiLow;
		double fCutMass_PhiHigh;
};

#endif
