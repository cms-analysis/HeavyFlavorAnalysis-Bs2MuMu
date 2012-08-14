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
		// Assumed resolution in µm
		double fResolutionXY;
		double fResolutionZ;
		
		TVector3 smearXY(TVector3 v);
		TVector3 smearZ(TVector3 v);
};

#endif
