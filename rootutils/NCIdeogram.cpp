//
//  NCIdeogram.cpp
//
//  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 27.06.13.
//
//

#include "NCIdeogram.h"

#include <map>

#include <TMath.h>
#include <TLeaf.h>
#include <TEventList.h>

static std::map<int,NCIdeogram> fGlobals;

NCIdeogram *NCIdeogram::getDefaultIdeogram(int ix)
{
	return &fGlobals[ix];
} // getDefaultIdeogram()

double NCIdeogram::evaluate(double x)
{
	size_t j;
	double result = 0.0;
	
	for (j = 0; j < fValues.size(); j++)
		result += evalGauss(x,fValues[j].first,fValues[j].second);
	
	return result;
} // evaluate()

double NCIdeogram::evalGauss(double x, double mu, double sigma)
{
	double result = (x - mu)/sigma;
	
	// compute exponential
	result = TMath::Exp(-0.5*result*result);
	
	// normalize
	result /=  TMath::Sqrt(2.0*TMath::Pi())*sigma;
	
	return result;
} // evalGauss()

void NCIdeogram::addTree(TTree *tree, const char *mu, const char *sigma)
{
	TLeaf *leafMu = tree->FindLeaf(mu);
	TLeaf *leafSigma = tree->FindLeaf(sigma);
	TEventList *elist = tree->GetEventList();
	Long64_t j,ix,nbr;
	
	nbr = elist ? elist->GetN() : tree->GetEntries();
	for (j = 0; j < nbr; j++) {
		ix = elist ? elist->GetEntry(j) : j;
		tree->GetEntry(ix);
		addPoint(leafMu->GetValue(), leafSigma->GetValue());
	}
} // addTree()
