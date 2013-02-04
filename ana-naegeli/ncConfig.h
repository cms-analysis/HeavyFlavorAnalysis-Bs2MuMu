/*
 *  ncConfig.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.02.13.
 *
 */

#ifndef NCCONFIG_H
#define NCCONFIG_H

#include <string>

#include <TCut.h>

class ncConfig {
	
	public:
		ncConfig();
		
		const char *getDataFile() { return fDataFile.c_str(); }
		const char *getMCFile() { return fMCFile.c_str(); }
		const char *getVarFile() { return fVarFile.c_str(); }
		const char *getStyle() { return fStyle.c_str(); }
		TCut getPreselection() { return fCutPreselection; }
	
	private:
		std::string fDataFile;
		std::string fMCFile;
		std::string fVarFile;
		std::string fStyle;
		TCut fCutPreselection;
};

#endif
