/*
 *  ncEvaluate.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 20.02.12.
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
		ncEvaluate(const char *inFile, TTree *inTree, const char *title);
		~ncEvaluate();
		
		double eval(int64_t j); // load tree entry j and evaluate
		double eval(); // evaluate the current vals
		
		std::map<std::string,Float_t> *getVals() {return &vals;}
		
	private:
		TTree *tree;
		const char *weight_file;
		TMVA::Reader *reader;
		std::string methodTitle;
		std::map<std::string,Float_t> vals;
		std::map<std::string,TTreeFormula*> formulas;
};

#endif
