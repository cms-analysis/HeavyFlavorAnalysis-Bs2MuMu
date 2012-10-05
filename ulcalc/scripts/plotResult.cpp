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
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TVirtualPad.h>
#include <RooNumber.h>
#include <RooWorkspace.h>

#include <RooStats/SamplingDistribution.h>
#include <RooStats/HypoTestInverterResult.h>

// upper limit strings
#define kVDEF_ulObs			"default-11:ulObs:val"
#define kVDEF_ulBkgMed		"default-11:ulBkg:val"
#define kVDEF_ulBkgErrHi	"default-11:ulBkg:errHi"
#define kVDEF_ulBkgErrLo	"default-11:ulBkg:errLo"
#define kVDEF_ulSMMed		"default-11:ulSM:val"
#define kVDEF_ulSMErrHi		"default-11:ulSM:errHi"
#define kVDEF_ulSMErrLo		"default-11:ulSM:errLo"

// two sided intervals strings
#define kVDEF_intObs		"default-11:intObs:val"
#define kVDEF_intObsErrHi	"default-11:intObs:errHi"
#define kVDEF_intObsErrLo	"default-11:intObs:errLo"
#define kVDEF_intSMExp		"default-11:intSM:val"
#define kVDEF_intSMErrHi	"default-11:intSM:errHi"
#define kVDEF_intSMErrLo	"default-11:intSM:errLo"

const static float kBFBsmm = 3.2e-9;

using namespace RooStats;

class clsPlotter {
	
	public:
	explicit clsPlotter(const char *filename) : fFilename(filename), fResultName("result_mu_s"), fSMBands(false), fDrawLegend(true), fPlotObs(true), fUseCLs(true), fObs(NULL), fMgSM(NULL), fG0SM(NULL), fG1SM(NULL), fG2SM(NULL), fMgBkg(NULL), fG0Bkg(NULL), fG1Bkg(NULL), fG2Bkg(NULL) {}
		
		void plot(double cl = 0.95);
		void print(double cl = 0.95, const char *latexName = NULL);
	
	public:
		void setSMBands(bool smBands) { fSMBands = smBands; }
		void setPlotObs(bool plotObs) { fPlotObs = plotObs; }
		void setUseCLs(bool useCLs) { fUseCLs = useCLs; }
		void setResultName(std::string resultName) { fResultName = resultName; }
		bool getSMBands() { return fSMBands; }
		bool getPlotObs() { return fPlotObs; }
		bool getUseCLs() { return fUseCLs; }
		std::string getResultName() { return fResultName; }
	
	private:
		void makeObs(bool reload = false);
		void makeSMBands(bool reload = false);
		void makeBkgBands(bool reload = false);
		HypoTestInverterResult *loadResult(TFile *file, const char *name);
		
		void ulFromPoints(double testSize, vector<pair<double,double> > *vals, double *outUL, double *outLL);
		
	private:
		std::string fFilename;
		std::string fResultName;
		bool fSMBands;
		bool fDrawLegend;
		bool fPlotObs;
		bool fUseCLs;
		
		TGraph *fObs;
		TMultiGraph *fMgSM;
		TGraph *fG0SM;
		TGraphAsymmErrors *fG1SM;
		TGraphAsymmErrors *fG2SM;
		
		TMultiGraph *fMgBkg;
		TGraph *fG0Bkg;
		TGraphAsymmErrors *fG1Bkg;
		TGraphAsymmErrors *fG2Bkg;
};

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
	HypoTestInverterResult *result = (HypoTestInverterResult*)file->Get(name);
	if (!result) {
		RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
		result = (HypoTestInverterResult*)wspace->obj(name);
	}
	
	result->UseCLs(fUseCLs);
	
	return result;
} // loadResult()

void clsPlotter::print(double cl, const char *latexName)
{
	TFile *file = TFile::Open(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(file, fResultName.c_str());
	double ul, ll, smPl,smMi, bf = 0;
	double x,y;
	double testSize = 1.0 - cl;
	vector<pair<double,double> > vals;
	Int_t j;
	FILE *latexFile = NULL;
	
	if (latexName) latexFile = fopen(latexName,"w");
	
	resultBkg->UseCLs(fUseCLs);
	resultBkg->SetTestSize(testSize);
	
	ul = resultBkg->UpperLimit();
	ll = resultBkg->LowerLimit();
	cout << "UL(Bs -> mumu) = " << ul*kBFBsmm << "\t(" << ul << ") @ " << (int)(cl*100.) << " % CL" << endl;
	if (!fUseCLs) cout << "LL(Bs -> mumu) = " << ll*kBFBsmm << "\t(" << ll << ") @ " << (int)(cl*100.) << " % CL" << endl;
	
	if (!fUseCLs) {
		// state the observed value due to the fit
		RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
		if (wspace) {
			wspace->pdf("total_pdf")->fitTo(*wspace->data("data"), RooFit::GlobalObservables(*wspace->set("nui")), RooFit::PrintLevel(-1));
			bf = wspace->var("mu_s")->getVal();
			cout << Form("BF(Bs -> mumu) = %e + %e - %e\t(%e + %e - %e)", bf*kBFBsmm, (ul-bf)*kBFBsmm, (bf-ll)*kBFBsmm, bf, ul-bf, bf-ll) << endl;
		}
	}	
	if(latexFile) {
		if(fUseCLs) fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulObs, ul*kBFBsmm); // print the observed upper limit
		else {
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intObs, bf*kBFBsmm);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intObsErrHi, (ul - bf)*kBFBsmm);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intObsErrLo, (bf - ll)*kBFBsmm);
		}
	}
	
	
	makeObs();
	vals.clear();
	for(j = 0; j < fObs->GetN(); j++) {
		if(fObs->GetPoint(j, x, y) >= 0)
			vals.push_back(make_pair(x,y));
	}
	ulFromPoints(testSize,&vals,&ul,&ll);
	cout << "UL(Bs -> mumu) = " << ul*kBFBsmm << "\t(" << ul << ") @ " << (int)(cl * 100.) << " % CL" << endl;
	if (!fUseCLs) cout << "LL(Bs -> mumu) = " << ll*kBFBsmm << "\t(" << ll << ") @ " << (int)(cl*100.) << " % CL" << endl;
	
	makeSMBands();
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
	
	if (fUseCLs)	cout << "EXP_SM[UL(Bs->mumu)] = " << ul*kBFBsmm << "+" << (smPl-ul)*kBFBsmm << "-" << (ul-smMi)*kBFBsmm << "\t(" << ul << "+" << smPl-ul << "-" << (ul-smMi) << ")" << endl;
	else			cout << "EXP_SM[BF(Bs->mumu)] = " << kBFBsmm << "+" << (ul-1.)*kBFBsmm << "-" << (1.-ll)*kBFBsmm << "\t(" << 1. << "+" << ul-1. << "-" << 1.-ll << ")" << endl;
	
	if (latexFile) {
		if (fUseCLs) {
			smPl = smPl - ul; smMi = ul - smMi;
			ul *= kBFBsmm; smPl *= kBFBsmm; smMi *= kBFBsmm;
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulSMMed, ul);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulSMErrHi, smPl);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulSMErrLo, smMi);
		} else {
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intSMExp, kBFBsmm);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intSMErrHi, (ul-1.)*kBFBsmm);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_intSMErrLo, (1.-ll)*kBFBsmm);
		}
	}
	
	// compute bkg expectation
	makeBkgBands();
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
	
	if(fUseCLs) cout << "EXP_Bkg[UL(Bs->mumu)] = " << ul*kBFBsmm << "+" << (smPl-ul)*kBFBsmm << "-" << (ul-smMi)*kBFBsmm << "\t(" << ul << "+" << smPl-ul << "-" << (ul-smMi) << ")" << endl;
	if (latexFile && fUseCLs) {
		smPl = smPl - ul; smMi = ul - smMi;
		ul *= kBFBsmm; smPl *= kBFBsmm; smMi *= kBFBsmm;
		fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulBkgMed, ul);
		fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulBkgErrHi, smPl);
		fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_ulBkgErrLo, smMi);
	}
	
	if (latexFile) fclose(latexFile);
	delete file;
} // print()

void clsPlotter::makeObs(bool reload)
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&f,fResultName.c_str());
	int j,size = result->ArraySize();
	std::vector<double> xvalues;
	std::vector<unsigned int> index(size);
	double clevel;
	
	if (fObs && reload) {
		delete fObs;
		fObs = NULL;
	}
	
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
} // makeObs()

void clsPlotter::makeSMBands(bool reload)
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(&f, fResultName.c_str());
	HypoTestInverterResult *resultSM = loadResult(&f, Form("%s_SM",fResultName.c_str()));
	vector<unsigned int> indexBkg(resultBkg->ArraySize());
	vector<unsigned int> indexSM(resultSM->ArraySize());
	HypoTestResult *hypoBkg,*hypoSM;
	SamplingDistribution *dist;
	double p[5];
	double q[5];
	double *x;
	unsigned int j,k,ixBkg,ixSM;
	vector<double> pvalues;
	vector<double> xvaluesBkg;
	vector<double> xvaluesSM;
	
	if (reload) {
		if(fMgSM) {
			delete fMgSM;
			fMgSM = NULL;
		}
		
		if(fG0SM) {
			delete fG0SM;
			fG0SM = NULL;
		}
		
		if(fG1SM) {
			delete fG1SM;
			fG1SM = NULL;
		}
		
		if(fG2SM) {
			delete fG2SM;
			fG2SM = NULL;
		}
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
	HypoTestInverterResult *result = loadResult(&file, fResultName.c_str());
	SamplingDistribution *dist;
	unsigned int j,ix,nbr = result->ArraySize();
	std::vector<unsigned int> index(nbr);
	std::vector<double> xvalues;
	double p[5];
	double q[5];
	double *x;
	
	if (reload) {
		if(fMgBkg) {
			delete fMgBkg;
			fMgBkg = NULL;
		}
		
		if(fG0Bkg) {
			delete fG0Bkg;
			fG0Bkg = NULL;
		}
		
		if(fG1Bkg) {
			delete fG1Bkg;
			fG1Bkg = NULL;
		}
		
		if(fG2Bkg) {
			delete fG2Bkg;
			fG2Bkg = NULL;
		}
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
	}
	mg->Draw( (fPlotObs ? "" : "A") );
	if (fPlotObs) fObs->Draw("same");
	
	// draw confidence level line
	double alpha = 1. - cl;
	line = new TLine(fObs->GetXaxis()->GetXmin(), alpha, fObs->GetXaxis()->GetXmax(), alpha);
	line->SetLineColor(kRed);
	line->Draw();
	
	// draw legend
	leg->Draw();
	
	gPad->RedrawAxis();
} // plot()

void plotCLs(const char *file, bool SMExp, const char *latexName = NULL, bool plotObs = true, double cl = 0.95)
{
	clsPlotter plot(file);
	
	plot.setResultName("Hybrid_Ul_mu_s");
	plot.setUseCLs(true);
	plot.setPlotObs(plotObs);
	plot.setSMBands(SMExp);
	plot.plot(cl);
	plot.print(cl, latexName);
} // plotCLs()

void plotTwoSided(const char *file, bool SMExp, const char *latexName =  NULL, bool plotObs = true, double cl = 0.68)
{
	clsPlotter plot(file);
	
	plot.setResultName("Hybrid_Int_mu_s");
	plot.setUseCLs(false);
	plot.setPlotObs(plotObs);
	plot.setSMBands(SMExp);
	plot.plot(cl);
	plot.print(cl, latexName);
} // plotTwoSided()
