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
#include <TH1.h>
#include <RooNumber.h>
#include <RooWorkspace.h>

#include <RooStats/SamplingDistribution.h>
#include <RooStats/HypoTestInverterResult.h>

const static float kBFBsmm = 3.2e-9;
const static float kBFBdmm = 1.0e-10;

using namespace RooStats;

class clsPlotter {
	
	public:
		explicit clsPlotter(const char *filename) : fFilename(filename), fResultName("result_mu_s"), fResultsDir(NULL), fSMBands(false), fDrawLegend(true), fPlotObs(true), fUseCLs(true), fAllDigits(true), fFirst(true), fObs(NULL), fMgSM(NULL), fG0SM(NULL), fG1SM(NULL), fG2SM(NULL), fMgBkg(NULL), fG0Bkg(NULL), fG1Bkg(NULL), fG2Bkg(NULL), fWspace(NULL) {}
		
		void plot(double cl = 0.95);
                void print(double cl, TGraphAsymmErrors *graph, double xpos);
	
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

void clsPlotter::print(double cl, TGraphAsymmErrors *graph, double xpos)
{
	TFile *file = TFile::Open(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(file, Form("%s_%s",fResultName.c_str(),fPOI.c_str()));
	double ul, ll, smPl,smMi, bf = 0;
	double x,y;
	double testSize = 1.0 - cl;
	vector<pair<double,double> > vals;
	Int_t j;
	FILE *latexFile = NULL;
	int npos;
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
				
				nbr = wspace->var(Form("PeakBkgBd_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("PeakBkgBd_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("PeakBkgBd_%d",k))->getVal();
				cout << Form("Npeak(%s,Bdmm) = ", channelName) << nbr << endl;
				
				nbr = wspace->var(Form("Pdd_%d",j))->getVal() * wspace->function(Form("NuD_%d",j))->getVal() * wspace->var("mu_d")->getVal();
				for (k = j+2; (var = wspace->var(Form("Pdd_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("Pdd_%d",k))->getVal() * wspace->function(Form("NuD_%d",k))->getVal() * wspace->var("mu_d")->getVal();
				cout << Form("Nsig(%s,Bdmm) = ",channelName) << nbr << endl;
				
				nbr = wspace->var(Form("TauS_%d",j))->getVal() * wspace->var(Form("nu_b_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("nu_b_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("TauS_%d",k))->getVal() * wspace->var(Form("nu_b_%d",k))->getVal();
				cout << Form("Ncomb(%s,Bsmm) = ",channelName) << nbr << endl;
				
				nbr = wspace->var(Form("PeakBkgBs_%d",j))->getVal();
				for (k = j+2; (var = wspace->var(Form("PeakBkgBs_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("PeakBkgBs_%d",k))->getVal();
				cout << Form("Npeak(%s,Bsmm) = ",channelName) << nbr << endl;
				
				nbr = wspace->var(Form("Pss_%d",j))->getVal() * wspace->function(Form("NuS_%d",j))->getVal() * wspace->var("mu_s")->getVal();
				for (k = j+2; (var = wspace->var(Form("Pss_%d",k))) != NULL; k+=2)
					nbr += wspace->var(Form("Pss_%d",k))->getVal() * wspace->function(Form("NuS_%d",k))->getVal() * wspace->var("mu_s")->getVal();
				cout << Form("Nsig(%s,Bsmm) = ",channelName) << nbr << endl;
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
	
	// add to the graph
	npos = graph->GetN();
	graph->SetPoint(npos, xpos, ul);
	if( !isfinite(smPl) )
	  smPl = smMi;
	graph->SetPointError(npos, 0, 0, smMi, smPl);

	if (fUseCLs)	cout << "EXP_SM[UL(" << fPOI << ")] = " << ul*theoreticalBF << "+" << (smPl-ul)*theoreticalBF << "-" << (ul-smMi)*theoreticalBF << "\t(" << ul << "+" << smPl-ul << "-" << (ul-smMi) << ")" << endl;
	else			cout << "EXP_SM[BF(" << fPOI << ")] = " << theoreticalBF << "+" << (ul-1.)*theoreticalBF << "-" << (1.-ll)*theoreticalBF << "\t(" << 1. << "+" << ul-1. << "-" << 1.-ll << ")" << endl;
	
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

void addEntry(TGraphAsymmErrors *graph, double xpos, int ix, const char *filename)
{
	const char *poi [] = {"mu_s","mu_d"};
	clsPlotter plot(filename);
	
	plot.setResultName("Hybrid_Ul");
	plot.setPOI(poi[ix]);
	
	plot.setUseCLs(true);
	
	try { plot.print(0.95,graph,xpos); }
	catch (std::string e) { }
} // plotCLs()

void scanCLs()
{
  TGraphAsymmErrors *graphBsmm = new TGraphAsymmErrors;
  TGraphAsymmErrors *graphBdmm = new TGraphAsymmErrors;
  TCanvas c;

  graphBsmm->SetNameTitle("bsmm_ul","Bsmm upper limit window search");
  graphBdmm->SetNameTitle("bdmm_ul","Bdmm upper limit window search");

  addEntry(graphBsmm, 5.27, 0, "/scratch/naegelic/anaBmm.plotResults.2012-527.ulc.root");
  addEntry(graphBsmm, 5.28, 0, "/scratch/naegelic/anaBmm.plotResults.2012-528.ulc.root");
  addEntry(graphBsmm, 5.29, 0, "/scratch/naegelic/anaBmm.plotResults.2012-529.ulc.root");
  addEntry(graphBsmm, 5.30, 0, "/scratch/naegelic/anaBmm.plotResults.2012-530.ulc.root");
  addEntry(graphBsmm, 5.31, 0, "/scratch/naegelic/anaBmm.plotResults.2012-531.ulc.root");
  addEntry(graphBsmm, 5.32, 0, "/scratch/naegelic/anaBmm.plotResults.2012-532.ulc.root");
  addEntry(graphBsmm, 5.33, 0, "/scratch/naegelic/anaBmm.plotResults.2012-533.ulc.root");
  addEntry(graphBsmm, 5.34, 0, "/scratch/naegelic/anaBmm.plotResults.2012-534.ulc.root");
  addEntry(graphBsmm, 5.35, 0, "/scratch/naegelic/anaBmm.plotResults.2012-535.ulc.root");

  addEntry(graphBdmm, 5.27, 1, "/scratch/naegelic/anaBmm.plotResults.2012-527.ulc.root");
  addEntry(graphBdmm, 5.28, 1, "/scratch/naegelic/anaBmm.plotResults.2012-528.ulc.root");
  addEntry(graphBdmm, 5.29, 1, "/scratch/naegelic/anaBmm.plotResults.2012-529.ulc.root");
  addEntry(graphBdmm, 5.30, 1, "/scratch/naegelic/anaBmm.plotResults.2012-530.ulc.root");
  addEntry(graphBdmm, 5.31, 1, "/scratch/naegelic/anaBmm.plotResults.2012-531.ulc.root");
  addEntry(graphBdmm, 5.32, 1, "/scratch/naegelic/anaBmm.plotResults.2012-532.ulc.root");
  addEntry(graphBdmm, 5.33, 1, "/scratch/naegelic/anaBmm.plotResults.2012-533.ulc.root");
  addEntry(graphBdmm, 5.34, 1, "/scratch/naegelic/anaBmm.plotResults.2012-534.ulc.root");
  addEntry(graphBdmm, 5.35, 1, "/scratch/naegelic/anaBmm.plotResults.2012-535.ulc.root");

  graphBsmm->Draw("alp");  
  graphBdmm->Draw("alp");

  graphBsmm->GetHistogram()->SetXTitle("Window Cut / GeV");
  graphBsmm->GetHistogram()->SetYTitle("BF / BF_{SM}");
  graphBdmm->GetHistogram()->SetXTitle("Window Cut / GeV");
  graphBdmm->GetHistogram()->SetYTitle("BF / BF_{SM}");

  new TFile("scan.root","update");
  graphBsmm->Write();
  graphBdmm->Write();
} // scanCLs()
