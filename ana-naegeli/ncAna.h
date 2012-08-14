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
#include "../macros/massReader.hh"
#include "../ulcalc/external_constants.h"

// Standard headers
#include <stdint.h>

// STL
#include <set>
#include <string>

// ROOT headers
#include <TH1D.h>
#include <TCut.h>
#include <TFile.h>

// RooFit
#include <RooWorkspace.h>

class ncAna {
	
	public:
		ncAna();
		~ncAna();
		
		// Possible 'analysis'
		void showSidebandSubtraction(bool massConstraint = true);
		double getSystematicsEffAna(TCut cut, bool massConstraint = false);
		double getSystematicsEffAna(std::set<std::string> cuts, bool massConstraint = false);
		
		// ulc files
		void computeBplus(unsigned channelIx, bool reload = true);
		void computeBmm(unsigned channelIx, bool bsmm, bool reload = true);
		void computePeaking(unsigned channelIx, bool reload = true);
		void appendPeakingChannel(int trueCand, measurement_t muMis1, measurement_t muMis2, measurement_t effFilter, measurement_t bf, unsigned channelIx, bool reload);
		measurement_t fitBplus(unsigned channelIx, bool reload = true, RooWorkspace **wout = NULL);
		void writeULC(const char *ulcname, bool reload = true);
		
		// load systematics
		void loadSystematics(bool def = true);
	public:
		// Standard cuts (static)
		static TCut cutTrigger(bool signal = true) { return signal ? TCut("triggered_bs") : TCut("triggered_jpsi"); }
		static TCut cutMuon() { return TCut("tight_mu1 && tight_mu2"); } 
		static TCut cutAcceptanceMC(int nKaons = 0);
		static TCut cutAcceptanceMCHard(int nKaons = 0);
		static TCut cutAcceptanceData(int nKaons = 0);
		static TCut cutNormCandGen(bool reco) { return TCut(Form("candidate == %d", (reco ? 68 : -68))); }
		static TCut cutNormCand() { return TCut("candidate == 300521"); }
		static TCut cutNormTruth() { return TCut("true_decay == 19"); }
		static TCut cutSigCand() { return TCut("candidate == 301313"); }
		static TCut cutSigCandGen(bool bsmm, bool reco);
		static TCut cutSigTruth(bool bsmm = true) { return TCut(Form("true_decay == %d", (bsmm ? kDecay_BsToMuMu : kDecay_BdToMuMu))); }
		static TCut cutSanity() { return TCut("d3e > 0 && ipe > 0"); }
		static TCut cutMVAPresel();
		
		// non-static cuts
		TCut cutSigWindow(bool bsmm) { return (bsmm ? TCut(Form("%f < mass && mass < %f",fBsWindowRegion.first, fBsWindowRegion.second)) : TCut(Form("%f < mass && mass < %f",fBdWindowRegion.first, fBdWindowRegion.second))); }
		TCut cutHisto() {return TCut(Form("%f < mass && mass < %f", fMassRange.first, fMassRange.second));}
		TCut cutBlindRegion() {return TCut(Form("%f < mass && mass < %f", fBlindedRegion.first, fBlindedRegion.second));}
		TCut cutChannel(unsigned ix);
		TCut cutAna(unsigned ix, bool jpsi = false);
		
		// Standard histograms
		TH1D *getDefaultMassHistogram(const char *name = "mass");
	
	public:
		// Accessors to configurable parameters
		void setNbrMassBins(int32_t massBins) { fNbrMassBins = massBins; }
		void setMassRange(std::pair<double,double> massRange) { fMassRange = massRange; }
		void setNbrSidebandBins(int32_t nbrBins) { fNbrSidebandBins = nbrBins; }
		void setLowSidebandRange(std::pair<double,double> lowRange) { fLowSidebandRange = lowRange; }
		void setHighSidebandRange(std::pair<double,double> highRange) { fHighSidebandRange = highRange; }
		void setPlotDir(const char *pDir) { fPlotDir = std::string(pDir); }
		void setChannels(std::vector<std::pair<ncCut::ncCutRange,std::string> > channels) { fChannels = channels; }
		
		int32_t getNbrMassBins() { return fNbrMassBins;}
		std::pair<double,double> getMassRange() { return fMassRange; }
		int32_t getNbrSidebandBins() { return fNbrSidebandBins; }
		std::pair<double,double> getLowSidebandRange() { return fLowSidebandRange; }
		std::pair<double,double> getHighSidebandRange() { return fHighSidebandRange; }
		const char *getPlotDir() { return fPlotDir.c_str(); }
		std::vector<std::pair<ncCut::ncCutRange,std::string> > getChannels() { return fChannels; }
	
	private:
		// Private datatypes
		typedef std::map<bmm_param,measurement_t> ulc_t;
		
	private:
		// filenames
		std::string fDataFileName;
		std::string fMCFileName;
		std::string fPeakFileName;
		std::string fAccFileName;
		
		// Configuration options
		int32_t fNbrMassBins;
		std::pair<double,double> fMassRange;
		int32_t fNbrSidebandBins;
		std::pair<double,double> fLowSidebandRange;
		std::pair<double,double> fHighSidebandRange;
		// regions
		std::pair<double,double> fBlindedRegion;
		std::pair<double,double> fBsWindowRegion;
		std::pair<double,double> fBdWindowRegion;
		std::pair<double,double> fFitRangeNorm;
		// saving plots
		std::string fPlotDir;
		
		// channels
		std::vector<std::pair<ncCut::ncCutRange,std::string> > fChannels;
		
		// ulc parameters
		ulc_t fUlcBsmm;
		ulc_t fUlcBdmm;
		ulc_t fUlcBplus;
		
		// systematics
		std::map<systematics_t,double> *fSystematicsTable;
		measurement_t fMisIDKaonPion;
		measurement_t fMisIDProton;
		
	private:
		// for sideband subtraction
		std::vector<ncCut> fSidebandCuts;
};

#endif
