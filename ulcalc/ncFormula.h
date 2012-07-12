/*
 *  ncFormula.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 02.07.12.
 *
 */

#ifndef NCFORMULA
#define NCFORMULA

// my headers
#include "ncTree.h"

// standard headers
#include <set>
#include <string>

// ROOT headers
#include <TCut.h>

class ncFormula {
	
	public:
		ncFormula();
		ncFormula(std::string name, bool op);
		
		std::string toString() { return fName; }
		
	public:
		// accessors
		void setName(std::string name) { fName = name; }
		void setOp(bool op) {fOp = op; }
		
		std::string getName() {return fName;}
		bool isOp() { return fOp; }
	private:
		std::string fName;
		bool fOp;
};

ncTree<ncFormula> *read_formula(std::string formula);
// nur simple cuts möglich!
std::set<std::string> get_cuts(ncTree<ncFormula> *tree);
std::set<std::string> get_dependencies(ncTree<ncFormula> *tree);

#endif
