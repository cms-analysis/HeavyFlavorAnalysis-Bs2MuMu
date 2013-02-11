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

#ifndef __CINT__
#include <libxml/parser.h>
#endif

class ncConfig {
	
	public:
		ncConfig();
		
#ifndef __CINT__
		ncConfig(xmlNodePtr node);
#endif
		
		const char *getDataFile() { return fDataFile.c_str(); }
		const char *getMCFile() { return fMCFile.c_str(); }
		const char *getVarFile() { return fVarFile.c_str(); }
		const char *getStyle() { return fStyle.c_str(); }
		TCut getPreselection() { return fCutPreselection; }
		TCut getMCSelection() { return fMCSelection; }
		TCut getDataSelection() { return fDataSelection; }
		
		void dump();
		
	private:
		std::string fDataFile;
		std::string fAccFile;
		std::string fMCFile;
		std::string fVarFile;
		std::string fStyle;
		TCut fCutPreselection;
		TCut fMCSelection;
		TCut fDataSelection;
};

#endif
