/*
 *  mlp_eff.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.05.12.
 *
 */

#include "ncPropErr.hh"
#include "ncEvaluate.hh"

#include <TFile.h>
#include <TEventList.h>
#include <TRandom3.h>

//#define DEBUG

using namespace std;

void uncertainty_mlp(TTree *tree, pair<double,double> range, double relInError)
{
	TFile dumpFile("error.root","recreate");
	TTree *dumpTree = new TTree("err","");
	ncEvaluate barEval("weights/ncMVA_MLP_0.weights.xml",tree,0);
//	ncEvaluate endEval("weights/ncMVA_MLP_1.weights.xml",tree,1);
	map<string,Float_t> *values;
	map<string,Float_t>::iterator it;
	Int_t k;
	Long64_t total = 0;
	Long64_t passedOrig = 0;
	Long64_t passedRand = 0;
	TEventList elist("elist");
	TRandom3 rd;
	double valOrig,valRand;
	double effOrig,effRand;
	
	// build the branches
	dumpTree->Branch("mlp",&valOrig,"mlp/D");
	dumpTree->Branch("mlpErr",&valRand,"mlpErr/D");
	
	cout << "Processing tree with " << tree->GetEntries() << " entries..." << endl;
	
	// get the entries after the preselection cut
	tree->Draw(">>elist", TCut("candidate == 301313 && d3e > 0 && ipe > 0 && true_decay == 1") && ncEvaluate::getPreselCut());
	elist.Sort();
	
	cout << elist.GetN() << " events preselected..." << endl;
	
	for (k = 0; k < elist.GetN(); k++) {
		
		total++; // increment the total process event counter
		
		// evaluate the entry
		valOrig = barEval.eval(elist.GetEntry(k));
		values = barEval.getVals();
		
		// shuffle the input value by 5 percent (for now)
		for (it = values->begin(); it != values->end(); ++it) {
//			cout << "'" << it->first << "':	" << it->second;
			it->second = rd.Gaus(it->second,relInError * (it->second));
//			cout << " -> " << it->second << endl;
		}
		
		valRand = barEval.eval(); // no reload of data
		
		// check the counters
#ifdef DEBUG
		cout << "Original value	" << valOrig << endl;
		cout << "Randomized value	" << valRand << endl;
#endif
		
		if (range.first < valOrig && valOrig < range.second)
			passedOrig++;
		
		if (range.first < valRand && valRand < range.second)
			passedRand++;
		
		dumpTree->Fill();
	}
	
	effOrig = (double)passedOrig / (double)total;
	effRand = (double)passedRand / (double)total;
	cout << "Efficiency Orig = " << effOrig << " = " << passedOrig << " / " << total << endl;
	cout << "Efficiency Rand = " << effRand << " = " << passedRand << " / " << total << endl;
	cout << "Delta eps = " << (effRand-effOrig)/effOrig << endl;
	
	dumpTree->Write(dumpTree->GetName(),TObject::kOverwrite);
	dumpFile.Close();
} // uncertainty_mlp()

// easy access
void test_uncertainty(double relErr)
{
	TFile f("/Users/cn/CMSData/Reduced/production-mix-general.root");
	TTree *tree = (TTree*)f.Get("T");
	
	uncertainty_mlp(tree,make_pair(0.99,1.5),relErr);
} // test_uncertainty()
