/*
 *  ncVarReader.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 14.12.12.
 *
 */

#include "ncVarReader.h"

#include <cctype>
#include <cstdio>

void ncVarReader::clear()
{
	fVecVars.clear();
} // clear()

int ncVarReader::nextChar(FILE *file)
{
	int c;
	while ( (c = fgetc(file)) != EOF) {
		if (!isspace(c))
			break;
	}
	
	return c;
} // nextChar()

std::string ncVarReader::nextAlphaNum(FILE *file)
{
	std::string allowed("./_+-[]");
	std::string result;
	int c;
	bool b;
	
	do {
		c = nextChar(file);
		b = isalnum(c) || (std::find(allowed.begin(),allowed.end(),c) != allowed.end());
		if (b) result.push_back(c);
	} while(b);
	
	ungetc(c, file);
	
	return result;
} // nextAlphaNum()

ncCut ncVarReader::readCut(FILE *file)
{
	int c;
	ncCut cut;
	ncCut::ncCutRange rg;
	
	c = nextChar(file);
	if (c != '{') throw std::string("No opening '{' found in ncVarReader::readCut()");
	
	// read alpha_num
	cut.setName(nextAlphaNum(file).c_str());
	
	if ( (c = nextChar(file)) != ',') throw std::string("No separating ',' found in ncVarReader::readCut()");
	cut.setFormula(nextAlphaNum(file).c_str());
	
	if ( (c = nextChar(file)) != ',') throw std::string("No separating ',' found in ncVarReader::readCut()");
	rg.first = atof(nextAlphaNum(file).c_str());
	
	if ( (c = nextChar(file)) != ',') throw std::string("No separating ',' found in ncVarReader::readCut()");
	rg.second = atof(nextAlphaNum(file).c_str());
	cut.setDomain(rg);
	
	if ( (c = nextChar(file)) != ',') throw std::string("No separating ',' found in ncVarReader::readCut()");
	rg.first = atof(nextAlphaNum(file).c_str());
	
	if ( (c = nextChar(file)) != ',') throw std::string("No separating ',' found in ncVarReader::readCut()");
	rg.second = atof(nextAlphaNum(file).c_str());
	cut.setCut(rg);
	
	if ( (c = nextChar(file)) != '}') throw std::string("No closeing '}' found in ncVarReader::readCut()");
	
	return cut;
} // readCut()

std::set<ncCut> ncVarReader::readEntry(FILE *file)
{
	int c;
	std::set<ncCut> theSet;
	ncCut cut;
	
	while (true) {
		cut = readCut(file);
		theSet.insert(cut);
		c = nextChar(file);
		if (c != ',') {
			ungetc(c, file);
			break;
		}
	}
	
	return theSet;
} // readEntry()

void ncVarReader::loadFile(const char *filename)
{
	FILE *file = fopen(filename, "r");
	std::string read;
	int c;
	std::set<ncCut> entry;
	
	while ( (c = nextChar(file)) != EOF ) {
		
		if (c != '{') throw std::string("No opening '{' found in ncVarReader::loadFile()");
		entry = readEntry(file);
		c = nextChar(file);
		if (c != '}') throw std::string("No closing '}' found in ncVarReader::loadFile()");
		
		fVecVars.push_back(entry);
	}
	
	fclose(file);
} // loadFile()
