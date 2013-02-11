/*
 *  ncConfig.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.02.13.
 *
 */

#include "ncConfig.h"

#include <iostream>

ncConfig::ncConfig()
{
} // ncConfig()

ncConfig::ncConfig(xmlNodePtr node)
{
	xmlNodePtr child;
	xmlChar *c;
	
	for (child = node->xmlChildrenNode; child != NULL; child = child->next) {
		if (xmlStrcmp(child->name,(const xmlChar*)"datafile") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fDataFile = std::string((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"acceptance") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fAccFile = std::string((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"mcfile") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fMCFile = std::string((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"variables") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fVarFile = std::string((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"style") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fStyle = std::string((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"preselection") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fCutPreselection = TCut((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"mcselection") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fMCSelection = TCut((const char*)c);
			xmlFree(c);
		} else if (xmlStrcmp(child->name,(const xmlChar*)"dataselection") == 0) {
			c = xmlNodeListGetString(child->doc, child->xmlChildrenNode, 1);
			fDataSelection = TCut((const char*)c);
			xmlFree(c);
		}
	}
} // ncConfig()

void ncConfig::dump()
{
	using std::cout; using std::endl;
	cout << "ncConfig:" << endl;
	cout << "	fDataFile = " << fDataFile << endl;
	cout << "	fAccFile = " << fAccFile << endl;
	cout << "	fMCFile = " << fMCFile << endl;
	cout << "	fVarFile = " << fVarFile << endl;
	cout << "	fStyle = " << fStyle << endl;
	cout << "	fCutPreselection = " << fCutPreselection.GetTitle() << endl;
	cout << Form("	fMCSelection = '%s'",fMCSelection.GetTitle()) << endl;
	cout << Form("	fDataSelection = '%s'",fDataSelection.GetTitle()) << endl;
} // dump()
