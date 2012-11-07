/*
 *  plotResult.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 19.09.12.
 *
 */

#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TVirtualPad.h>
#include <RooNumber.h>
#include <RooWorkspace.h>

#include <RooStats/SamplingDistribution.h>
#include <RooStats/HypoTestInverterResult.h>

// upper limit strings
#define kVDEF_ulObs			"2011:ulObs:val"
#define kVDEF_ulObsExponent	"2011:ulObs:exponent"
#define kVDEF_ulBkgMed		"2011:ulBkg:val"
#define kVDEF_ulBkgErrHi	"2011:ulBkg:errHi"
#define kVDEF_ulBkgErrLo	"2011:ulBkg:errLo"
#define kVDEF_ulBkgExponent	"2011:ulBkg:exponent"
#define kVDEF_ulSMMed		"2011:ulSM:val"
#define kVDEF_ulSMErrHi		"2011:ulSM:errHi"
#define kVDEF_ulSMErrLo		"2011:ulSM:errLo"
#define kVDEF_ulSMExponent	"2011:ulSM:exponent"

// two sided intervals strings
#define kVDEF_intObs			"2011:intObs:val"
#define kVDEF_intObsExponent	"2011:intObs:exponent"
#define kVDEF_intObsErrHi		"2011:intObs:errHi"
#define kVDEF_intObsErrLo		"2011:intObs:errLo"
#define kVDEF_intSMExp			"2011:intSM:val"
#define kVDEF_intSMErrHi		"2011:intSM:errHi"
#define kVDEF_intSMErrLo		"2011:intSM:errLo"
#define kVDEF_intSMExponent		"2011:intSM:exponent"

// for table comparison
#define kVDEF_tblCombBsmm	"2011:fitComb:val:bsmm"
#define kVDEF_tblPeakBsmm	"2011:fitPeak:val:bsmm"
#define kVDEF_tblSigBsmm	"2011:fitSig:val:bsmm"
#define kVDEF_tblCombBdmm	"2011:fitComb:val:bdmm"
#define kVDEF_tblPeakBdmm	"2011:fitPeak:val:bdmm"
#define kVDEF_tblSigBdmm	"2011:fitSig:val:bdmm"

const static float kBFBsmm = 3.2e-9;
const static float kBFBdmm = 1.0e-10;

using namespace RooStats;

class clsPlotter {
	
	public:
		explicit clsPlotter(const char *filename) : fFilename(filename), fResultName("result_mu_s"), fResultsDir(NULL), fSMBands(false), fDrawLegend(true), fPlotObs(true), fUseCLs(true), fAllDigits(true), fFirst(true), fObs(NULL), fMgSM(NULL), fG0SM(NULL), fG1SM(NULL), fG2SM(NULL), fMgBkg(NULL), fG0Bkg(NULL), fG1Bkg(NULL), fG2Bkg(NULL), fWspace(NULL) {}
		
		void plot(double cl = 0.95);
		void print(double cl = 0.95);
	
	public:
		void setSMBands(bool smBands) { fSMBands = smBands; }
		void setPlotObs(bool plotObs) { fPlotObs = plotObs; }
		void setUseCLs(bool useCLs) { fUseCLs = useCLs; }
		void setAllDigits(bool allDigits) { fAllDigits = allDigits; }
		void setResultName(std::string resultName) { fResultName = resultName; }
		void setResultsDir(const char *resultsDir);
		void setPOI(const char *poi) { fPOI = std::string(poi); }
		bool getSMBands() { return fSMBands; }
		bool getPlotObs() { return fPlotObs; }
		bool getUseCLs() { return fUseCLs; }
		bool getAllDigits() { return fAllDigits; }
		std::string getResultName() { return fResultName; }
		const char *getResultsDir() { return (fResultsDir ? fResultsDir->c_str() : NULL); }
		const char *getPOI() { return fPOI.c_str(); }
	
	private:
		void makeObs(bool reload = false);
		void makeSMBands(bool reload = false);
		void makeBkgBands(bool reload = false);
		HypoTestInverterResult *loadResult(TFile *file, const char *name);
		
		void ulFromPoints(double testSize, vector<pair<double,double> > *vals, double *outUL, double *outLL);
		
	private:
		std::string fFilename;
		std::string fResultName;
		std::string *fResultsDir;
		std::string fPOI;
		bool fSMBands;
		bool fDrawLegend;
		bool fPlotObs;
		bool fUseCLs;
		bool fAllDigits;
		bool fFirst;
		
		TGraph *fObs;
		TMultiGraph *fMgSM;
		TGraph *fG0SM;
		TGraphAsymmErrors *fG1SM;
		TGraphAsymmErrors *fG2SM;
		
		TMultiGraph *fMgBkg;
		TGraph *fG0Bkg;
		TGraphAsymmErrors *fG1Bkg;
		TGraphAsymmErrors *fG2Bkg;
		
		RooWorkspace *fWspace;
		std::map<std::string,HypoTestInverterResult*> fResults;
};

void clsPlotter::setResultsDir(const char *resultsDir)
{
	std::string tmp( (resultsDir ? resultsDir : "") );
	if (fResultsDir) delete fResultsDir;
	fResultsDir = NULL;
	if (resultsDir)
		fResultsDir = new std::string(tmp);
} // setResultsDir()

void clsPlotter::ulFromPoints(double testSize, vector<pair<double,double> > *vals, double *outUL, double *outLL)
{
	size_t j;
	double y1,y0;
	vector<double> uls;
	vector<double> lls;

	std::sort(vals->begin(), vals->end());
	
	for (j = 1; j < vals->size(); j++) {
		y1 = (*vals)[j-0].second;
		y0 = (*vals)[j-1].second;
		
		if (y1 <= testSize && testSize <= y0)
			uls.push_back(((*vals)[j].first - (*vals)[j-1].first)/((*vals)[j].second - (*vals)[j-1].second)*(testSize - (*vals)[j-1].second) + (*vals)[j-1].first);
		if (y0 <= testSize && testSize <= y1)
			lls.push_back(((*vals)[j].first - (*vals)[j-1].first)/((*vals)[j].second - (*vals)[j-1].second)*(testSize - (*vals)[j-1].second) + (*vals)[j-1].first);
	}
	
	// this should not be needed, but lets do it anyway
	std::sort(uls.begin(),uls.end());
	std::sort(lls.begin(),lls.end());
	if(outUL)
		*outUL = (uls.size() > 0 ? uls.back() : numeric_limits<double>::quiet_NaN());
	if(outLL)
		*outLL = (lls.size() > 0 ? lls.front() : numeric_limits<double>::quiet_NaN());
} // ulFromPoints()

HypoTestInverterResult *clsPlotter::loadResult(TFile *file, const char *name)
{
	HypoTestInverterResult *result = NULL;
	if (fResults.count(std::string(name)) > 0)
		result = fResults[std::string(name)];
	
	if (!result) {
		if(!fWspace)
			fWspace = (RooWorkspace*)file->Get("wspace");
		
		fResults[std::string(name)] = result = (HypoTestInverterResult*)fWspace->obj(name);
	}
	
	return result;
} // loadResult()

void clsPlotter::print(double cl)
{
	TFile *file = TFile::Open(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(file, Form("%s_%s",fResultName.c_str(),fPOI.c_str()));
	double ul, ll, smPl,smMi, bf = 0;
	double x,y;
	double testSize = 1.0 - cl;
	vector<pair<double,double> > vals;
	Int_t j;
	FILE *latexFile = NULL;
	int expo;
	double theoreticalBF = -1;
	
	if (!resultBkg)
		goto bail;
	
	if (fPOI.compare("mu_s") == 0) {
		theoreticalBF = kBFBsmm;
	} else if (fPOI.compare("mu_d") == 0) {
		theoreticalBF = kBFBdmm;
	} else {
		cerr << "clsPlotter::print() running with unknown POI: " << fPOI << endl;
		abort();
	}

	
	if (fResultsDir) {
		latexFile = fopen(Form("%s/%s", fResultsDir->c_str(), (fUseCLs ? "cls.tex" : "clsplusb.tex")), (fFirst ? "w" : "a"));
		fFirst = false;
	}
	
	resultBkg->UseCLs(fUseCLs);
	resultBkg->SetTestSize(testSize);
	
	ul = resultBkg->UpperLimit();
	ll = resultBkg->LowerLimit();
	cout << "UL(" << fPOI << ") = " << ul*theoreticalBF << "\t(" << ul << ") @ " << (int)(cl*100.) << " % CL" << endl;
	if (!fUseCLs) cout << "LL("<< fPOI << ") = " << ll*theoreticalBF << "\t(" << ll << ") @ " << (int)(cl*100.) << " % CL" << endl;
	
	if (!fUseCLs) {
		// state the observed value due to the fit
		RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
		if (wspace) {
			wspace->pdf("total_pdf")->fitTo(*wspace->data("data"), RooFit::GlobalObservables(*wspace->set("nui")), RooFit::PrintLevel(-1));
			bf = wspace->var(fPOI.c_str())->getVal();
		}
		
		// print the numbers for the table
		if (fPOI.compare("mu_s") == 0) {
			double nbr;
			int k;
			
			for (j = 0; j <= 1; j++) {
				const char *channelName = j == 0 ? "barrel" : "endcap";
				RooRealVar *var;
				
				nbr = wspace->var(Form("TauD_%d",j))->getVal() * wspace->var(Form("nu_b_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("nu_b_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("TauD_%d",k))->getVal() * wspace->var(Form("nu_b_%d",k))->getVal();
				cout << Form("Ncomb(%s,Bdmm) = ", channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblCombBdmm,channelName), nbr);
				
				nbr = wspace->var(Form("PeakBkgBd_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("PeakBkgBd_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("PeakBkgBd_%d",k))->getVal();
				cout << Form("Npeak(%s,Bdmm) = ", channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblPeakBdmm,channelName), nbr);
				
				nbr = wspace->var(Form("Pdd_%d",j))->getVal() * wspace->function(Form("NuD_%d",j))->getVal() * wspace->var("mu_d")->getVal();
				for (k = j+2; (var = wspace->var(Form("Pdd_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("Pdd_%d",k))->getVal() * wspace->function(Form("NuD_%d",k))->getVal() * wspace->var("mu_d")->getVal();
				cout << Form("Nsig(%s,Bdmm) = ",channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblSigBdmm,channelName), nbr);
				
				nbr = wspace->var(Form("TauS_%d",j))->getVal() * wspace->var(Form("nu_b_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("nu_b_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("TauS_%d",k))->getVal() * wspace->var(Form("nu_b_%d",k))->getVal();
				cout << Form("Ncomb(%s,Bsmm) = ",channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblCombBsmm,channelName), nbr);
				
				nbr = wspace->var(Form("PeakBkgBs_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("PeakBkgBs_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("PeakBkgBs_%d",k))->getVal();
				cout << Form("Npeak(%s,Bsmm) = ",channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblPeakBsmm,channelName), nbr);
				
				nbr = wspace->var(Form("Pss_%d",j))->getVal() * wspace->function(Form("NuS_%d",j))->getVal() * wspace->var("mu_s")->getVal();
				for (k = j+2; (var = wspace->var(Form("Pss_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("Pss_%d",k))->getVal() * wspace->function(Form("NuS_%d",k))->getVal() * wspace->var("mu_s")->getVal();
				cout << Form("Nsig(%s,Bsmm) = ",channelName) << nbr << endl;
				if (latexFile)
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.3f}}}\n", Form("%s:%s",kVDEF_tblSigBsmm,channelName), nbr);
			}
		}
	}
	
	makeObs(true);
	vals.clear();
	for(j = 0; j < fObs->GetN(); j++) {
		if(fObs->GetPoint(j, x, y) >= 0)
			vals.push_back(make_pair(x,y));
	}
	ulFromPoints(testSize,&vals,&ul,&ll);
	cout << "UL("<< fPOI << ") = " << ul*theoreticalBF << "\t(" << ul << ") @ " << (int)(cl * 100.) << " % CL" << endl;
	if (!fUseCLs) {
		cout << "LL(" << fPOI << ") = " << ll*theoreticalBF << "\t(" << ll << ") @ " << (int)(cl*100.) << " % CL" << endl;
		cout << Form("BF(%s) = %e + %e - %e\t(%e + %e - %e)", fPOI.c_str(), bf*theoreticalBF, (ul-bf)*theoreticalBF, (bf-ll)*theoreticalBF, bf, ul-bf, bf-ll) << endl;
	}
	
	if(latexFile) {
		if(fUseCLs) {
			// print the observed upper limit
			ul *= theoreticalBF;
			if (fAllDigits) {
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulObs,fPOI.c_str()), ul);
			} else {
				expo = floor(log10(ul));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulObs,fPOI.c_str()), ul*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_ulObsExponent,fPOI.c_str()), expo);
			}
		}
		else {
			// print the observed branching fraction
			if (isnan(ll)) ll = 0;
			ul = ul - bf; ll = bf - ll;
			bf *= theoreticalBF; ul *= theoreticalBF; ll *= theoreticalBF;
			if (fAllDigits) {
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intObs,fPOI.c_str()), bf);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intObsErrHi,fPOI.c_str()), ul);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intObsErrLo,fPOI.c_str()), ll);
			} else {
				expo = floor(log10(bf));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intObs,fPOI.c_str()), bf*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intObsErrHi,fPOI.c_str()), ul*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intObsErrLo,fPOI.c_str()), ll*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_intObsExponent,fPOI.c_str()), expo);
			}
		}
	}
	
	makeSMBands(true);
	vals.clear();
	for(j = 0; j < fG0SM->GetN(); j++) {
		if (fG0SM->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x, y));
	}
	ulFromPoints(testSize,&vals,&ul,&ll);
	
	vals.clear();
	for(j = 0; j < fG1SM->GetN(); j++) {
		if (fG1SM->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x, y + fG1SM->GetErrorYhigh(j)));
	}
	ulFromPoints(testSize,&vals,&smPl,NULL);
	
	vals.clear();
	for(j = 0; j < fG1SM->GetN(); j++) {
		if (fG1SM->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x,y - fG1SM->GetErrorYlow(j)));
	}
	ulFromPoints(testSize,&vals,&smMi,NULL);
	
	if (fUseCLs)	cout << "EXP_SM[UL(" << fPOI << ")] = " << ul*theoreticalBF << "+" << (smPl-ul)*theoreticalBF << "-" << (ul-smMi)*theoreticalBF << "\t(" << ul << "+" << smPl-ul << "-" << (ul-smMi) << ")" << endl;
	else			cout << "EXP_SM[BF(" << fPOI << ")] = " << theoreticalBF << "+" << (ul-1.)*theoreticalBF << "-" << (1.-ll)*theoreticalBF << "\t(" << 1. << "+" << ul-1. << "-" << 1.-ll << ")" << endl;
	
	if (latexFile) {
		if (fUseCLs) {
			smPl = smPl - ul; smMi = ul - smMi;
			ul *= theoreticalBF; smPl *= theoreticalBF; smMi *= theoreticalBF;
			if(fAllDigits) {
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulSMMed,fPOI.c_str()), ul);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulSMErrHi,fPOI.c_str()), smPl);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulSMErrLo,fPOI.c_str()), smMi);
			} else {
				expo = floor(log10(ul));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulSMMed,fPOI.c_str()), ul*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulSMErrHi,fPOI.c_str()), smPl*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulSMErrLo,fPOI.c_str()), smMi*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_ulSMExponent,fPOI.c_str()), expo);
			}
		} else {
			if (fAllDigits) {
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intSMExp,fPOI.c_str()), theoreticalBF);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intSMErrHi,fPOI.c_str()), (ul-1.)*theoreticalBF);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_intSMErrLo,fPOI.c_str()), (1.-ll)*theoreticalBF);
			} else {
				if (isnan(ll)) ll = 0.;
				expo = floor(log10(theoreticalBF));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intSMExp,fPOI.c_str()), theoreticalBF*pow(10., -expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intSMErrHi,fPOI.c_str()), (ul-1.)*theoreticalBF*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_intSMErrLo,fPOI.c_str()), (1.-ll)*theoreticalBF*pow(10.,-expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_intSMExponent,fPOI.c_str()), expo);
			}
		}
	}
	
	// compute bkg expectation
	makeBkgBands(true);
	vals.clear();
	for(j = 0; j < fG0Bkg->GetN(); j++) {
		if (fG0Bkg->GetPoint(j,x,y) >= 0) {
			vals.push_back(make_pair(x, y));
		}
	}
	ulFromPoints(testSize,&vals,&ul,NULL);
	
	vals.clear();
	for(j = 0; j < fG1Bkg->GetN(); j++) {
		if (fG1Bkg->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x, y + fG1Bkg->GetErrorYhigh(j)));
	}
	ulFromPoints(testSize,&vals,&smPl,NULL);
	
	vals.clear();
	for(j = 0; j < fG1Bkg->GetN(); j++) {
		if (fG1Bkg->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x,y - fG1Bkg->GetErrorYlow(j)));
	}
	ulFromPoints(testSize,&vals,&smMi,NULL);
	
	if(fUseCLs) cout << "EXP_Bkg[UL(" << fPOI << ")] = " << ul*theoreticalBF << "+" << (smPl-ul)*theoreticalBF << "-" << (ul-smMi)*theoreticalBF << "\t(" << ul << "+" << smPl-ul << "-" << (ul-smMi) << ")" << endl;
	if (latexFile && fUseCLs) {
		smPl = smPl - ul; smMi = ul - smMi;
		ul *= theoreticalBF; smPl *= theoreticalBF; smMi *= theoreticalBF;
		if (fAllDigits) {
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulBkgMed,fPOI.c_str()), ul);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulBkgErrHi,fPOI.c_str()), smPl);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_ulBkgErrLo,fPOI.c_str()), smMi);
		} else {
			expo = floor(log10(ul));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulBkgMed,fPOI.c_str()), ul*pow(10., -expo));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulBkgErrHi,fPOI.c_str()), smPl*pow(10.,-expo));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_ulBkgErrLo,fPOI.c_str()), smMi*pow(10.,-expo));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_ulBkgExponent,fPOI.c_str()), expo);
		}
	}
	
bail:
	if (latexFile) fclose(latexFile);
	delete file;
} // print()

void clsPlotter::makeObs(bool reload)
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&f,Form("%s_%s",fResultName.c_str(),fPOI.c_str()));
	int j,size = result ? result->ArraySize() : 0;
	std::vector<double> xvalues;
	std::vector<unsigned int> index(size);
	double clevel;
	
	if (!result)
		goto bail;
	
	if (fObs && reload)
		fObs = NULL;
	
	if (!fObs) {
		fObs = new TGraph;
		
		fObs->SetTitle(Form("%s observed", (fUseCLs ? "CLs" : "CLs+b")));
		fObs->SetLineWidth(2);
		
		// sequential access
		for (j = 0; j < result->ArraySize(); j++) xvalues.push_back(result->GetXValue(j));
		TMath::SortItr(xvalues.begin(), xvalues.end(), index.begin(), false);
		
		for (j = 0; j < size; j++) {
		  clevel = fUseCLs ? result->CLs(index[j]) : result->CLsplusb(index[j]);
		  if (clevel < 0) clevel = 1.0;
		  fObs->SetPoint(j, result->GetXValue(index[j]), clevel);
		}
	}
bail:
	return;
} // makeObs()

void clsPlotter::makeSMBands(bool reload)
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(&f, Form("%s_%s",fResultName.c_str(),fPOI.c_str()));
	HypoTestInverterResult *resultSM = loadResult(&f, Form("%s_%s_SM",fResultName.c_str(),fPOI.c_str()));
	vector<unsigned int> indexBkg(resultBkg ? resultBkg->ArraySize() : 0);
	vector<unsigned int> indexSM(resultSM ? resultSM->ArraySize() : 0);
	HypoTestResult *hypoBkg,*hypoSM;
	SamplingDistribution *dist;
	double p[5];
	double q[5];
	double *x;
	unsigned int j,k,ixBkg,ixSM;
	vector<double> pvalues;
	vector<double> xvaluesBkg;
	vector<double> xvaluesSM;
	
	if (!resultBkg || !resultSM)
		goto bail;
	
	if (reload) {
		if(fMgSM) fMgSM = NULL;
		if(fG0SM) fG0SM = NULL;
		if(fG1SM) fG1SM = NULL;
		if(fG2SM) fG2SM = NULL;
	}
	
	if(fMgSM) goto bail; // already done
	
	fMgSM = new TMultiGraph;
	fG0SM = new TGraph;
	fG1SM = new TGraphAsymmErrors;
	fG2SM = new TGraphAsymmErrors;
	
	fG0SM->SetTitle(Form("Expected SM %s - Median",(fUseCLs ? "CLs" : "CLs+b")));
	fG1SM->SetTitle(Form("Expected SM %s #pm 1 #sigma",(fUseCLs ? "CLs" : "CLs+b")));
	fG2SM->SetTitle(Form("Expected SM %s #pm 2 #sigma",(fUseCLs ? "CLs" : "CLs+b")));
	
	// sequential access
	for (j = 0; j < (unsigned int)resultBkg->ArraySize(); j++) xvaluesBkg.push_back(resultBkg->GetXValue(j));
	for (j = 0; j < (unsigned int)resultSM->ArraySize(); j++) xvaluesSM.push_back(resultSM->GetXValue(j));
	TMath::SortItr(xvaluesBkg.begin(), xvaluesBkg.end(), indexBkg.begin(), false);
	TMath::SortItr(xvaluesSM.begin(), xvaluesSM.end(), indexSM.begin(), false);
	
	p[0] = ROOT::Math::normal_cdf(-2);
	p[1] = ROOT::Math::normal_cdf(-1);
	p[2] = 0.5;
	p[3] = ROOT::Math::normal_cdf(1);
	p[4] = ROOT::Math::normal_cdf(2);
	
	for (j = 0; j < indexBkg.size(); j++) {
		HypoTestResult tempResult;
		
		ixBkg = indexBkg[j];
		ixSM = indexSM[j];
		
		if (resultBkg->GetXValue(ixBkg) != resultSM->GetXValue(ixSM)) {
			cerr << "TWO INVERTER RESULTS DO NOT MATCH!!!! ABORTING..." << endl;
			abort();
		}
		
		hypoBkg = resultBkg->GetResult(ixBkg);
		hypoSM = resultSM->GetResult(ixSM);
		
		tempResult.SetPValueIsRightTail( hypoBkg->GetPValueIsRightTail() );
		tempResult.SetBackgroundAsAlt( hypoBkg->GetBackGroundIsAlt() );
		tempResult.SetNullDistribution( hypoBkg->GetNullDistribution() );
		tempResult.SetAltDistribution( hypoBkg->GetAltDistribution() );
		
		// SM distribution saved in background
		dist = hypoSM->GetBackGroundIsAlt() ? hypoSM->GetAltDistribution() : hypoSM->GetNullDistribution();
		const vector<Double_t> &vec = dist->GetSamplingDistribution();
		pvalues.clear();
		for (k = 0; k < vec.size(); k++) {
			tempResult.SetTestStatisticData(vec[k]);
			pvalues.push_back(fUseCLs ? tempResult.CLs() : tempResult.CLsplusb());
		}
		
		x = const_cast<double*>(&pvalues[0]);
		TMath::Quantiles(pvalues.size(), 5, x, q, p, false);
		
		fG0SM->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		
		fG1SM->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		fG1SM->SetPointEYlow(j, q[2] - q[1]);
		fG1SM->SetPointEYhigh(j, q[3] - q[2]);
		
		fG2SM->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		fG2SM->SetPointEYlow(j, q[2] - q[0]);
		fG2SM->SetPointEYhigh(j, q[4] - q[2]);
	}
	
	fG2SM->SetFillColor(kYellow);
	fMgSM->Add(fG2SM,"3");
	fG1SM->SetFillColor(kGreen);
	fMgSM->Add(fG1SM,"3");
	fG0SM->SetLineStyle(2);
	fG0SM->SetLineWidth(2);
	fMgSM->Add(fG0SM,"L");
	
bail:
	return;
} // makeSMBands()

void clsPlotter::makeBkgBands(bool reload)
{
	TFile file(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&file, Form("%s_%s",fResultName.c_str(),fPOI.c_str()));
	SamplingDistribution *dist;
	unsigned int j,ix,nbr = (result ? result->ArraySize() : 0);
	std::vector<unsigned int> index(nbr);
	std::vector<double> xvalues;
	double p[5];
	double q[5];
	double *x;
	
	if (!result)
		goto bail;
	
	if (reload) {
		if(fMgBkg) fMgBkg = NULL;
		if(fG0Bkg) fG0Bkg = NULL;
		if(fG1Bkg) fG1Bkg = NULL;
		if(fG2Bkg) fG2Bkg = NULL;
	}
	
	if (fMgBkg) goto bail;
	
	fMgBkg = new TMultiGraph;
	fG0Bkg = new TGraph;
	fG1Bkg = new TGraphAsymmErrors;
	fG2Bkg = new TGraphAsymmErrors;
	
	fG0Bkg->SetTitle(Form("Expected %s - Median", (fUseCLs ? "CLs" : "CLs+b")));
	fG1Bkg->SetTitle(Form("Expected %s #pm 1 #sigma", (fUseCLs ? "CLs" : "CLs+b")));
	fG2Bkg->SetTitle(Form("Expected %s #pm 2 #sigma", (fUseCLs ? "CLs" : "CLs+b")));
	
	// seq access
	for (j = 0; j < (unsigned int)result->ArraySize(); j++) xvalues.push_back(result->GetXValue(j));
	TMath::SortItr(xvalues.begin(), xvalues.end(), index.begin(), false);
	p[0] = ROOT::Math::normal_cdf(-2);
	p[1] = ROOT::Math::normal_cdf(-1);
	p[2] = 0.5;
	p[3] = ROOT::Math::normal_cdf(1);
	p[4] = ROOT::Math::normal_cdf(2);
	
	result->UseCLs(fUseCLs);
	for (j = 0; j < nbr; j++) {
		ix = index[j];
		
		dist = result->GetExpectedPValueDist(ix);
		const std::vector<double> &values = dist->GetSamplingDistribution();
		x = const_cast<double*>(&values[0]);
		
		TMath::Quantiles(values.size(), 5, x, q, p, false);
		
		fG0Bkg->SetPoint(j, result->GetXValue(ix), q[2]);
		
		fG1Bkg->SetPoint(j, result->GetXValue(ix), q[2]);
		fG1Bkg->SetPointEYlow(j, q[2] - q[1]); // -1 sigma error
		fG1Bkg->SetPointEYhigh(j, q[3] - q[2]); // +1 sigma error
		
		fG2Bkg->SetPoint(j, result->GetXValue(ix), q[2]);
		fG2Bkg->SetPointEYlow(j, q[2] - q[0]);
		fG2Bkg->SetPointEYhigh(j, q[4] - q[2]);
		
		delete dist;
	}
	
	fG2Bkg->SetFillColor(kYellow);
	fMgBkg->Add(fG2Bkg,"3");
	fG1Bkg->SetFillColor(kGreen);
	fMgBkg->Add(fG1Bkg,"3");
	fG0Bkg->SetLineStyle(2);
	fG0Bkg->SetLineWidth(2);
	fMgBkg->Add(fG0Bkg,"L");
bail:
	return;
} // makeBkgBands()

void clsPlotter::plot(double cl)
{
	TLegend *leg = NULL;
	TObject *obj;
	TLine *line;
	int j,nbr;
	TMultiGraph *mg;
	TGraph *gr = NULL;
	
	makeObs();
	makeSMBands();
	makeBkgBands();
	mg = fSMBands ? fMgSM : fMgBkg;
	
	if (fDrawLegend) {
		leg = new TLegend(0.53,0.63,0.83,0.83,"","NDC");
		
		if(fPlotObs)
		  leg->AddEntry(fObs,"","L");
		
		nbr = mg->GetListOfGraphs()->GetSize();
		for (j = nbr-1; j>=0; j--)
			if ((obj = mg->GetListOfGraphs()->At(j)) != NULL) leg->AddEntry(obj,"", ((j == nbr-1) ? "L" : "F"));
		
		leg->SetFillColor(kWhite);
		leg->SetLineWidth(0);
		leg->SetLineColor(0);
		leg->SetShadowColor(kWhite);
	}
	
	// draw the graphs
	if (fPlotObs) {
		fObs->Draw("AL");
		fObs->GetXaxis()->SetTitle("BF / BF_{SM}");
		fObs->GetYaxis()->SetTitle(Form("CL_{%s}", (fUseCLs ? "s" : "s+b")));
		fObs->GetYaxis()->SetTitleOffset(1.0);
	} else {
		gr = fSMBands ? fG0SM : fG0Bkg;
		gr->Draw("AL");
		gr->GetXaxis()->SetTitle("BF / BF_{SM}");
		gr->GetYaxis()->SetTitle(Form("CL_{%s}", (fUseCLs ? "s" : "s+b")));
		gr->GetYaxis()->SetTitleOffset(1.0);
	}
	mg->Draw("");
	if (fPlotObs) {
		fObs->Draw("same");
		gr = fObs;
	}
	
	// draw confidence level line
	double alpha = 1. - cl;
	line = new TLine(gr->GetXaxis()->GetXmin(), alpha, gr->GetXaxis()->GetXmax(), alpha);
	line->SetLineColor(kRed);
	line->Draw();
	
	// draw legend
	leg->Draw();
	
	gPad->RedrawAxis();
} // plot()

void plotCLs(const char *file, const char *resultsDir = NULL, bool plotObs = true, bool allDigits = true, double cl = 0.95)
{
	TCanvas *c;
	const char *poi [] = {"mu_s","mu_d"};
	unsigned j;
	clsPlotter plot(file);
	
	for (j = 0; j < sizeof(poi)/sizeof(const char *); j++) {
		
		plot.setResultName("Hybrid_Ul");
		plot.setPOI(poi[j]);
		
		plot.setResultsDir(resultsDir);
		plot.setAllDigits(allDigits);
		plot.setUseCLs(true);
		plot.setPlotObs(plotObs);
		
		try { plot.print(cl); }
		catch (std::string e) { continue; }
		
		plot.setSMBands(false);
		c = new TCanvas;
		plot.plot(cl);
		
		if (resultsDir)
			c->SaveAs(Form("%s/cls_graph_bkg_%s.pdf",resultsDir,poi[j]));
		
		plot.setSMBands(true);
		c = new TCanvas;
		
		plot.plot(cl);
		if (resultsDir)
			c->SaveAs(Form("%s/cls_graph_sm_%s.pdf",resultsDir,poi[j]));
	}
} // plotCLs()

void plotTwoSided(const char *file, const char *resultsDir =  NULL, bool plotObs = true, bool allDigits = true, double cl = 0.68)
{
	TCanvas *c;
	const char *poi [] = {"mu_s","mu_d"};
	clsPlotter plot(file);
	unsigned j;
	
	for (j = 0; j < sizeof(poi)/sizeof(const char *); j++) {
		plot.setResultName("Hybrid_Int");
		plot.setPOI(poi[j]);
		plot.setResultsDir(resultsDir);
		plot.setAllDigits(allDigits);
		plot.setUseCLs(false);
		plot.setPlotObs(plotObs);
		
		try { plot.print(cl); }
		catch (std::string e) { continue; }
		
		plot.setSMBands(false);
		c = new TCanvas;
		plot.plot(cl);
		
		if (resultsDir)
			c->SaveAs(Form("%s/clsplusb_graph_bkg_%s.pdf",resultsDir,poi[j]));
		
		plot.setSMBands(true);
		c = new TCanvas;
		plot.plot(cl);
		
		if (resultsDir)
			c->SaveAs(Form("%s/clsplusb_graph_sm_%s.pdf",resultsDir,poi[j]));
	}
} // plotTwoSided()
