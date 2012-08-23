/*
 *  pixelReader.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 30.07.12.
 *
 */

#ifndef PIXELREADER_H
#define PIXELREADER_H

#include "treeReader01.hh"
#include "massReader.hh"

#include <set>

#include <TRandom3.h>

class pixelReader : public treeReader01 {
	
	public:
		pixelReader(TChain *tree, TString evtClassName);
		virtual ~pixelReader();
		
		virtual void bookHist();
		virtual void eventProcessing();
		virtual void closeHistFile();
		virtual void readCuts(TString filename, int dump = 1);
	
	private:
		void clearVariables();
		bool loadCandidateVariables(TAnaCand *pCand);
		double compIso(TAnaCand *pCand);
		int loadTruth(TAnaCand *pCand);
		void buildDecay(TGenCand *pGen, decay_t *dec);
		
		// utility function
		void dumpGenerator();
	
	private:
		std::set<int> fStableParticles;
		std::map<decay_t,int> fDecayTable;
	
	private:
		TTree *reduced_tree;
		
		/* candidate variables */
		float fMass;
		float fPt;
		float fEta;
		
		/* vertex variables */
		float fDoca;
		float fDocaZ;
		float fDocaXY;
		float fD3;
		float fD3Truth;
		float fAlpha;
		
		/* muon variables */
		float fPtMu1;
		float fPtMu2;
		float fEtaMu1;
		float fEtaMu2;
		
		/* other */
		float fIso;
		
		/* pv variables */
		float fPvZ;
		float fPvXY;
		
		/* truth info */
		int fTrueDecay;
	
	private:
		// random number generator
		TRandom3 fRand;
		// Track resolution
		double fD0Resolution;
		double fDzResolution;
		double fPhiResolution;
		double fCotThetaResolution;
		double fPtResolution;
		// primary vertex resolution
		double fPVResolutionXY;
		double fPVResolutionZ;
		
		void smearTrack(TVector3 *v, TVector3 *p);
		TVector3 smearPhi(TVector3 v, double res_cm);
		TVector3 smearZ(TVector3 v, double res_cm);
		TVector3 smearD0(TVector3 v,double res_cm, TVector3 plab);
		TVector3 smearR(TVector3 v, double res_cm);
	
	private:
		// resolution histograms
		TH2D *fHistoResIP_XY;
		TH2D *fHistoResIP_Z;
		
		void readResolution();
	
	private:
		unsigned fNumCands;
	
	private:
		// DEBUG. Remove again afterwards
		TH1D *hXres;
		TH1D *hYres;
		TH1D *hZres;
};

#endif
