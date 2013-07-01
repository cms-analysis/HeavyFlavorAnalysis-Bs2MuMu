//
//  NCIdeogram.h
//
//  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 27.06.13.
//
//

#ifndef NCIdeogram_
#define NCIdeogram_

#include <vector>
#include <utility>

#include <TTree.h>

class NCIdeogram {
	
	public:
		NCIdeogram() {}
		
		double evaluate(double x);
		
		void clear() { fValues.clear(); }
		void addPoint(double mu, double sigma) { fValues.push_back(std::make_pair(mu,sigma)); }
		void addTree(TTree *tree, const char *mu, const char *sigma);
		
		size_t getN() { return fValues.size(); }
	
	public:
		static NCIdeogram *getDefaultIdeogram(int ix);
	
	private:
		double evalGauss(double x, double mu, double sigma);
	
	private:
		std::vector<std::pair<double,double> > fValues;
};

#endif
