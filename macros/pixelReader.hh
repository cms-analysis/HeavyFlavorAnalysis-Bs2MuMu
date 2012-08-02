/*
 *  pixelReader.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 30.07.12.
 *
 */

#ifndef PIXELREADER_H
#define PIXELREADER_H

#include "treeReader01.hh"

class pixelReader : public treeReader01 {
	
	public:
		pixelReader(TChain *tree, TString evtClassName);
		virtual ~pixelReader();
		
		virtual void bookHist();
		virtual void eventProcessing();
		virtual void closeHistFile();
	
	private:
		bool loadCandidateVariables(TAnaCand *pCand);
		double calcMaxDoca(TAnaCand *pCand);
	
	private:
		TTree *reduced_tree;
		
		float fMass;
		float fPt;
		float fDoca;
};

#endif
