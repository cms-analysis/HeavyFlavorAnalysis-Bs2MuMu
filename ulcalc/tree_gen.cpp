/*
 *  tree_gen.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 23.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "../rootutils/NCRootUtils.h"

// Standard headers
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

// ROOT headers
#include <TFile.h>
#include <TTree.h>

using namespace std;

struct cut_t {
	
	cut_t(string theName, float_t theMin, float_t theMax, uint32_t thePoints, bool theLowerBound) : name(theName),cut_min(theMin),cut_max(theMax),points(thePoints),lower_bound(theLowerBound) {}
	string name;
	float_t cut_min;
	float_t cut_max;
	float_t value;
	uint32_t points;
	bool lower_bound;
};

void saveCube(vector<cut_t>::iterator begin, vector<cut_t>::iterator end, TTree *tree)
{
	uint32_t j,npoints;
	float_t *val;
	float_t min,max;
	bool low;
	
	if (begin == end) {
		tree->Fill();
		goto bail;
	}
	
	npoints = begin->points;
	max = begin->cut_max;
	min = begin->cut_min;
	low = begin->lower_bound;
	val = &(begin->value);
	
	begin++; // go to next
	for (j = 0; j < npoints; j++) {
		if (npoints == 1)	*val = min; // take minimum
		else				*val = (float_t)j / (float_t)(npoints - 1) * (max - min) + min;
		
		// negative values for upper cut...
		if (!low) *val = -(*val);
		saveCube(begin,end,tree);
	}
	
bail:
	return;
} // saveCube()


static void usage()
{
	std::cout << "tree_gen <configfile>" << endl;
	abort();
} // usage()

static void build_cube(function_range_t *rg, void* param)
{
	vector<cut_t> *cube = static_cast<vector<cut_t>*> (param);
	
	if (rg->vals.size() < 3) {
		cerr << "ERROR: Couldn't parse config line for variable '" << rg->name.Data() << "'." << endl;
		goto bail;
	}
	
	cube->push_back(cut_t(string(rg->name.Data()),(float_t)rg->vals[0],(float_t)rg->vals[1],(uint32_t)rg->vals[2], (rg->options.size() > 0) ? rg->options[0] == 'l' : true));
bail:
	return;
} // build_cube()

int main(int argc, char *argv [])
{
	TFile f("grid_tree.root","recreate");
	TTree *tree = new TTree("GT","Tree with grid points");
	vector<cut_t> cube;
	
	if (argc < 2) {
		usage();
		goto bail;
	}
	
	parse_variable_config(argv[1], build_cube, &cube);
	
	for (vector<cut_t>::iterator it = cube.begin(); it != cube.end(); ++it) {
		cout << it->name.c_str() << (it->lower_bound ? " >" : " <") << " [" << it->cut_min << ", " << it->cut_max << "] with N = " << it->points << endl;
		tree->Branch(it->name.c_str(), &(it->value), Form("%s/F",it->name.c_str()));
	}
	
	saveCube(cube.begin(),cube.end(),tree);
	
	// save
	tree->Write();
	f.Close();
	
bail:
	return 0;
} // main()
