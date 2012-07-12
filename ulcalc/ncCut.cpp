/*
 *  ncCut.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 28.06.12.
 *
 */

#include "ncCut.h"

ncCut::ncCut()
{
	init("",ncCutRange(-1e30,1e30),ncCutRange(-1e30,1e30),"");
} // ncCut()

ncCut::ncCut(const char *formula, ncCutRange domain, ncCutRange cut, const char *name)
{
	init(formula, domain, cut, (name ? name : formula) );
} // ncCut()

void ncCut::init(const char *formula, ncCutRange domain, ncCutRange cut, const char *name)
{
	fName = std::string(name);
	fFormula = std::string(formula);
	fDomain = domain;	// Adapt RooFit convention for infinity
	fCut = cut;
} // init()
