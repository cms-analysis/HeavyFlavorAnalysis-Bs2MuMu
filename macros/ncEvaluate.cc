/*
 *  ncEvaluate.cpp
 *  cmssw
 *
 *  Created by Christoph on 20.02.12.
 *  Copyright 2012 PSI. All rights reserved.
 *
 */

#include <TCut.h>

#include "ncEvaluate.hh"

using namespace std;

ncEvaluate::ncEvaluate(const char *inFile, TTree *inTree, int channelIx) : tree(inTree), weight_file(inFile),methodTitle(Form("MLP_%d",channelIx))
{
	reader = new TMVA::Reader("v");
	map<string,string> *names = getDefaultVariables();
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
	double value;
	map<string,Float_t>::iterator it;
	map<string,TTreeFormula*>::iterator formIt;
	
	tree->GetEntry(j);
	for(it = vals.begin(); it != vals.end(); ++it)
		it->second = formulas[it->first]->EvalInstance();
	
	value = (Float_t)reader->EvaluateMVA(methodTitle);
	return value;
} // eval()

TCut ncEvaluate::getPreselCut()
{
	return TCut("4 < pt_mu1 && pt_mu1 < 40 && 4.0 < pt_mu2 && pt_mu2 < 20 && pt < 60 && ip < 0.03 && ipe > 0 && ip/ipe < 6 && d3/d3e < 100 && d3 < 2 && alpha < 0.2 && chi2/Ndof < 4 && iso_mor12 > .6 && ntrk < 15 && doca0 < .1 && TMath::Abs(eta) < 2.5 && 0 < d3e && 4.9 < mass && mass < 5.9 && TMath::Abs(eta_mu1) < 2.4 && TMath::Abs(eta_mu2) < 2.4 && (track_qual_mu1 & 4) && (track_qual_mu2 & 4)");
} // preselCut()

map<string,string>* ncEvaluate::getDefaultVariables()
{
	static map<string,string> *vars = NULL;
	
	if (!vars) {
		vars = new map<string,string>;
		vars->insert(pair<string,string>("pt","pt"));
		vars->insert(pair<string,string>("pt_mu1","pt_mu1"));
		vars->insert(pair<string,string>("pt_mu2","pt_mu2"));
		vars->insert(pair<string,string>("eta","eta"));
		
		vars->insert(pair<string,string>("sig3d","d3 / d3e"));
		vars->insert(pair<string,string>("d3","d3"));
		vars->insert(pair<string,string>("alpha","alpha"));
		vars->insert(pair<string,string>("normChi2","chi2/Ndof"));
		vars->insert(pair<string,string>("ip","ip"));
		vars->insert(pair<string,string>("sigip","ip / ipe"));
		
		vars->insert(pair<string,string>("iso_mor12","iso_mor12"));
		vars->insert(pair<string,string>("doca0","doca0"));
		vars->insert(pair<string,string>("ntrk","ntrk"));
	}
	
	return vars;
} // getDefaultVariables()
