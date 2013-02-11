/*
 *  ncAna2.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.02.13.
 *
 */

#ifndef NCANA2_H
#define NCANA2_H

#include "ncConfig.h"
#include "ncCut.h"

#include <set>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TFile.h>

class ncAna2 {
	
	public:
		ncAna2();
		
		// show the Variable plots
		void showAllVarPlots();
		void showVarPlots(ncConfig *conf);
	
	public:
		// static routines
		static void setHistoStyle(TH1D *h, const char *style);
	
	private:
		std::string fConfigFile;
	
	private:
		void readConfig();
		TCut buildCut(std::set<ncCut> *vars, const char *exclude);
		std::vector<ncConfig> fAnalyses;
	
	private:
		int fDefaultBinning;
		TFile fWorkFile;
};

#endif
