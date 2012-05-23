/*
 *  mlp_eff.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.05.12.
 *
 */

// Standard headers
#include <utility>

// ROOT headers
#include <TTree.h>

void uncertainty_mlp(TTree *tree, std::pair<double,double> range, double relInError);
void test_uncertainty(double relErr = 0.05);
