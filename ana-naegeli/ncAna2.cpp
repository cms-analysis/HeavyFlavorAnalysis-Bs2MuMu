/*
 *  ncAna2.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.02.13.
 *
 */

#include "ncAna2.h"
#include "ncVarReader.h"

#include <iostream>
#include <libxml/parser.h>

#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>

using std::cout;
using std::cerr;
using std::endl;

ncAna2::ncAna2() : fConfigFile("configs/setup2012.xml"), fDefaultBinning(50), fWorkFile("ANA.root","update")
{
} // ncAna2()

void ncAna2::readConfig()
{
	xmlDocPtr doc = xmlParseFile(fConfigFile.c_str());
	xmlNodePtr node;
	
	// remove the old configuration
	fAnalyses.clear();
	
	if (!doc) {
		cerr << Form("Unable to open xml file '%s'",fConfigFile.c_str()) << endl;
		goto bail;
	}
	
	node = xmlDocGetRootElement(doc);
	if (xmlStrcmp(node->name, (const xmlChar*)"ana") != 0) {
		cerr << Form("No 'ana' topnode found in '%s'",fConfigFile.c_str()) << endl;
		goto bail;
	}
	
	for (node = node->xmlChildrenNode; node != NULL; node = node->next) {
		if (xmlStrcmp(node->name, (const xmlChar*)"ana-decay") == 0) {
			ncConfig conf(node);
			conf.dump();
			fAnalyses.push_back(conf);
		}
	}
	
bail:
	if (doc) xmlFreeDoc(doc);
} // readConfig()

void ncAna2::showAllVarPlots()
{
	size_t j;
		
	// reload the configuration
	readConfig();
	
	// do the variable plot for the analyses given
	for (j = 0; j < fAnalyses.size(); j++)
		showVarPlots(&fAnalyses[j]);
} // showVarPlots()

void ncAna2::showVarPlots(ncConfig *conf)
{
	std::set<ncCut> vars;
	std::set<ncCut>::const_iterator it;
	ncVarReader reader;
	TFile *dataFile = TFile::Open(conf->getDataFile());
	TFile *mcFile = TFile::Open(conf->getMCFile());
	TTree *dataTree = (TTree*)dataFile->Get("T");
	TTree *mcTree = (TTree*)mcFile->Get("T");
	TCut presel = conf->getPreselection();
	TCut mcCut = conf->getMCSelection();
	TCut dataCut = conf->getDataSelection();
	TCut cut;
	TCanvas *c;
	
	try {
		reader.loadFile(conf->getVarFile());
	}
	catch (std::string err) {
		cerr << Form("ncAna2::showVarPlots(): error reading variable file: '%s'",err.c_str()) << endl;
		goto bail;
	}
	
	if (reader.getNbr() == 0) {
		cerr << "ncAna2::showVarPlots(): no variables entry found." << endl;
		goto bail;
	}
	
	vars = *reader.getVars(0);
	
	fWorkFile.cd();
	// show the histogram
	for (it = vars.begin(); it != vars.end(); ++it) {
		cout << "======================" << endl;
		cout << "Processing:" << endl;
		it->dump();
		cout << "======================" << endl;
		
		// build the additional cut
		cut = buildCut(&vars, it->getName());
		TH1D* histoMC = new TH1D(Form("%s_mc",it->getName()),"",fDefaultBinning,it->getDomain().first,it->getDomain().second);
		TH1D* histoData = new TH1D(Form("%s_data",it->getName()),"",fDefaultBinning,it->getDomain().first,it->getDomain().second);
		
		setHistoStyle(histoMC,conf->getStyle());
		
		dataTree->Draw(Form("%s >> %s_data",it->getFormula(),it->getName()), presel && cut && dataCut);
		mcTree->Draw(Form("%s >> %s_mc",it->getFormula(),it->getName()), presel && cut && mcCut);
		
		// normalize to unity
		histoMC->Scale(1./histoMC->Integral());
		histoData->Scale(1./histoData->Integral());
		
		c = new TCanvas;
		histoMC->Draw();
		histoData->Draw("sameE1");
		
		histoMC->Write(histoMC->GetName(),TObject::kOverwrite);
		histoData->Write(histoData->GetName(), TObject::kOverwrite);
	}
	
bail:
	if (dataFile) delete dataFile;
	if (mcFile) delete mcFile;
	return;
} // showVarPlots()

void ncAna2::setHistoStyle(TH1D *h, const char *style)
{
	h->SetLineColor(kBlue);
	h->SetFillColor(kBlue);
	h->SetFillStyle(3005);
} // setHistoStyle()

TCut ncAna2::buildCut(std::set<ncCut> *vars, const char *exclude)
{
	TCut cut;
	std::set<ncCut>::const_iterator it;
	
	for (it = vars->begin(); it != vars->end(); ++it) {
		if (strcmp(it->getName(),exclude) == 0)
			continue;
		
		cut = cut && TCut(Form("%e < %s && %s < %e",it->getCut().first,it->getFormula(),it->getFormula(),it->getCut().second));
	}
	
	return cut;
} // buildCut()
