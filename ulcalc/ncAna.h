/*
 *  ncAna.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 28.06.12.
 *
 */

#ifndef NCANA_H
#define NCANA_H

// my headers
#include "ncCut.h"

// Standard headers
#include <stdint.h>

// STL
#include <set>
#include <string>

// ROOT headers
#include <TH1D.h>
#include <TCut.h>
#include <TFile.h>

class ncAna {
	
	public:
		ncAna();
		~ncAna();
		
		// Possible 'analysis'
		void showSidebandSubtraction(bool massConstraint = true);
		double getSystematicsEffAna(TCut cut, bool massConstraint = false);
		double getSystematicsEffAna(std::set<std::string> cuts, bool massConstraint = false);
	
	public:
		// Standard cuts
		static TCut cutTrigger(bool signal = true) { return signal ? TCut("triggered_bs") : TCut("triggered_jpsi"); }
		static TCut cutMuon() { return TCut("tight_mu1 && tight_mu2"); } 
		static TCut cutNormCand() { return TCut("candidate == 300521"); }
		static TCut cutNormTruth() { return TCut("true_decay == 19"); }
		
		// Standard histograms
		TH1D *getDefaultMassHistogram(const char *name = "mass");
	
	private:
		void initFiles();
	
	private:
		// Pointers to ROOT files containing the Reduced trees
		TFile *dataFile;
		TFile *mcFile;
	
	public:
		// Accessors to configurable parameters
		void setNbrMassBins(int32_t massBins) { fNbrMassBins = massBins; }
		void setMassRange(std::pair<double,double> massRange) { fMassRange = massRange; }
		void setNbrSidebandBins(int32_t nbrBins) { fNbrSidebandBins = nbrBins; }
		void setLowSidebandRange(std::pair<double,double> lowRange) { fLowSidebandRange = lowRange; }
		void setHighSidebandRange(std::pair<double,double> highRange) { fHighSidebandRange = highRange; }
		void setPlotDir(const char *pDir) { fPlotDir = std::string(pDir); }
		
		int32_t getNbrMassBins() { return fNbrMassBins;}
		std::pair<double,double> getMassRange() { return fMassRange; }
		int32_t getNbrSidebandBins() { return fNbrSidebandBins; }
		std::pair<double,double> getLowSidebandRange() { return fLowSidebandRange; }
		std::pair<double,double> getHighSidebandRange() { return fHighSidebandRange; }
		const char *getPlotDir() { return fPlotDir.c_str(); }
		
	private:
		// Configuration options
		int32_t fNbrMassBins;
		std::pair<double,double> fMassRange;
		int32_t fNbrSidebandBins;
		std::pair<double,double> fLowSidebandRange;
		std::pair<double,double> fHighSidebandRange;
		// saving plots
		std::string fPlotDir;
		
	private:
		// for sideband subtraction
		std::vector<ncCut> fSidebandCuts;
};

#endif
