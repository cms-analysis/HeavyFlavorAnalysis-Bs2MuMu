/*
 *  ncEvaluate.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 20.02.12.
 *
 */

#include "ncEvaluate.h"
#include "ncMVA.h"

#include <TCut.h>

using namespace std;

ncEvaluate::ncEvaluate(const char *inFile, TTree *inTree, const char *title) : tree(inTree), weight_file(inFile), methodTitle(title)
{
	reader = new TMVA::Reader("v");
	map<string,string> *names = ncMVA::getMVAVariables();
	map<string,string>::const_iterator it;
	const char *c;
	
	for (it = names->begin(); it != names->end(); ++it) {
		
		if (it->first.compare(it->second) == 0)
			c = it->first.c_str();
		else
			c = Form("%s := %s", it->first.c_str(), it->second.c_str());
		
		reader->AddVariable(c,&vals[it->first.c_str()]); // implicitly creates the entry in vals
		
		// associated treeformula
		formulas[it->first] = new TTreeFormula("",it->second.c_str(),tree);
	}
	
	reader->BookMVA(methodTitle,weight_file);
} // ncEvaluate()

ncEvaluate::~ncEvaluate()
{
	delete reader;
	for(map<string,TTreeFormula*>::iterator it = formulas.begin(); it != formulas.end(); ++it)
		delete it->second;
} // ~ncEvalulate()

double ncEvaluate::eval(int64_t j)
{
	map<string,Float_t>::iterator it;
	map<string,TTreeFormula*>::iterator formIt;
	
	tree->GetEntry(j);
	for(it = vals.begin(); it != vals.end(); ++it)
		it->second = formulas[it->first]->EvalInstance();
	
	return eval();
} // eval()

double ncEvaluate::eval()
{
	double value = (double)reader->EvaluateMVA(methodTitle);
	return value;
} // eval()
