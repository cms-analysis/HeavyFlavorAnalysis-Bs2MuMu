/*
 *  ncEvaluate.h
 *  cmssw
 *
 *  Created by Christoph on 20.02.12.
 *  Copyright 2012 PSI. All rights reserved.
 *
 */

#ifndef __NCEVALUATE__
#define __NCEVALUATE__

#include <stdint.h>
#include <set>
#include <map>
#include <string>

#include <TTree.h>
#include <TTreeFormula.h>
#include <TMVA/Reader.h>

class ncEvaluate {
	
	public:
		ncEvaluate(const char *inFile, TTree *inTree, int channelIx);
		~ncEvaluate();
		
		double eval(int64_t j);
		
		static TCut getPreselCut();
		static std::map<std::string,std::string>* getDefaultVariables();
	
	private:
		TTree *tree;
		const char *weight_file;
		TMVA::Reader *reader;
		std::string methodTitle;
		std::map<std::string,Float_t> vals;
		std::map<std::string,TTreeFormula*> formulas;
};

#endif
