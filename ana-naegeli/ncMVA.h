/*
 *  ncMVA.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.04.12.
 */

#ifndef __NCMVA__
#define __NCMVA__

#include "ncCut.h"

#include <set>
#include <map>
#include <string>
#include <utility>
#include <TTree.h>
#include <TCut.h>

class ncMVA {
	public:
		ncMVA();
		~ncMVA();
		
		void splitTree(TTree *tree, bool save = true);
		void prepTraining(unsigned channelIx, std::set<ncCut> *vars);
		void runTraining(unsigned split, unsigned channelIx, bool prep = true);
		void runAllTrainings();
		
		void evalFile(const char *filename, bool progressReport = true);
		void evalFile(const char *filename, unsigned split, unsigned channelIx, bool progressReport = true);
		void evalFile(const char *filename, std::set<std::pair<unsigned,unsigned> >*mvas, bool progressReport = true);
		TTree *evalTree(TTree *tree, std::set<std::pair<unsigned,unsigned> > *mvas, bool progressReport = true);
		
		std::set<ncCut> getMVAVariables();
	public:
		// Accessor / Settor functions
		std::string getMVAOpts(unsigned channelIx) { return fMVAOpts[channelIx]; }
		std::string getMCPath() {return fMCPath;}
		std::string getDataPath() {return fDataPath;}
		
		void setMVAOpts(unsigned channelIx, std::string opts) {fMVAOpts[channelIx] = opts;}
		void setMCPath(std::string path) {fMCPath = path;}
		void setDataPath(std::string path) {fDataPath = path;}
	private:
		typedef std::pair<Long64_t,Long64_t> ncPair;
	
	private:
		std::string fMVAOpts[2];
		
		std::string fMCPath;
		std::string fDataPath;
		std::string fVarPath;
		
		double fSigWeight;
		double fBkgWeight;
	
	private:
		std::string *fTrainFilename;
};

#endif
