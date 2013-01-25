/*
 *  ncVarReader.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 14.12.12.
 *
 */

#ifndef NCVARREADER_H
#define NCVARREADER_H

#include "ncCut.h"

#include <string>
#include <set>
#include <vector>

class ncVarReader {
	
	public:
		void clear();
		void loadFile(const char *filename);
		std::set<ncCut> *getVars(size_t ix) { return (ix < fVecVars.size()) ? &(fVecVars[ix]) : NULL; }
		size_t getNbr() { return fVecVars.size(); }
		
	private:
		std::string nextAlphaNum(FILE *file);
		int nextChar(FILE *file);
		ncCut readCut(FILE *file);
		std::set<ncCut> readEntry(FILE *file);
		
	private:
		std::vector<std::set<ncCut> > fVecVars;
};

#endif
