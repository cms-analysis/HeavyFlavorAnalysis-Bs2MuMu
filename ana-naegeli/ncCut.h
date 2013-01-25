/*
 *  ncCut.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 28.06.12.
 *
 */

#ifndef NCCUT_H
#define NCCUT_H

// Standard headers
#include <utility>
#include <string>

class ncCut {
	
	public:
		typedef std::pair<double,double> ncCutRange;
	
	public:
		ncCut();
		ncCut(const char *formula, ncCutRange domain, ncCutRange cut, const char *name = NULL);
		
		const char *getName() const { return fName.c_str(); }
		const char *getFormula() const { return fFormula.c_str(); }
		ncCutRange getDomain() const { return fDomain; }
		ncCutRange getCut() const { return fCut; }
		
		void setName(const char *name) { fName = std::string(name); }
		void setFormula(const char *formula) { fFormula = std::string(formula); }
		void setDomain(ncCutRange domain) { fDomain = domain; }
		void setCut(ncCutRange cut) { fCut = cut; }
		
		void dump() const;
	
	private:
		
		void init(const char *formula, ncCutRange domain, ncCutRange cut, const char *name);
		
		std::string fName;
		std::string fFormula;
		ncCutRange fDomain;	// meaningful region
		ncCutRange fCut;	// cut on this variable
};

bool operator<(ncCut c1, ncCut c2);

#endif
