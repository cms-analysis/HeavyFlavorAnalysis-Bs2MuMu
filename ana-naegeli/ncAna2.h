/*
 *  ncAna2.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.02.13.
 *
 */

#ifndef NCANA2_H
#define NCANA2_H

#include "ncConfig.h"

#include <string>
#include <vector>

#include <TH1D.h>

class ncAna2 {
	
	public:
		ncAna2();
		
		// show the Variable plots
		void showAllVarPlots();
		void showVarPlots(ncConfig *conf);
	
	private:
		std::string fConfigFile;
	
	private:
		void readConfig();
		std::vector<ncConfig> fAnalyses;
	
	private:
		int fDefaultBinning;
		
		void setHistoStyle(TH1D *h, const char *style);
};

#endif
