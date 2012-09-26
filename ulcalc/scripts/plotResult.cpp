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
#include <RooNumber.h>
#include <RooWorkspace.h>

#include <RooStats/SamplingDistribution.h>
#include <RooStats/HypoTestInverterResult.h>

const static double kBFBsmm = 3.2e-9;

using namespace RooStats;

class clsPlotter {
	
	public:
	explicit clsPlotter(const char *filename) : fFilename(filename), fResultName("result_mu_s"), fSMBands(false), fDrawLegend(true), fPlotObs(true), fObs(NULL), fMgSM(NULL), fG0SM(NULL), fG1SM(NULL), fG2SM(NULL), fMgBkg(NULL), fG0Bkg(NULL), fG1Bkg(NULL), fG2Bkg(NULL) {}
		
		void plot(double cl = 0.95);
		void print(double cl = 0.95);
	
	public:
		void setSMBands(bool smBands) { fSMBands = smBands; }
		void setPlotObs(bool plotObs) { fPlotObs = plotObs; }
		bool getSMBands() { return fSMBands; }
		bool getPlotObs() { return fPlotObs; }
	
	private:
		void makeObs(bool reload = false);
		void makeSMBands(bool reload = false);
		void makeBkgBands(bool reload = false);
		HypoTestInverterResult *loadResult(TFile *file, const char *name);
		
		double ulFromPoints(double testSize, vector<pair<double,double> > *vals);
		
	private:
		std::string fFilename;
		std::string fResultName;
		bool fSMBands;
		bool fDrawLegend;
		bool fPlotObs;
		
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

double clsPlotter::ulFromPoints(double testSize, vector<pair<double,double> > *vals)
{
	size_t j;
	double y1,y0;
	sort(vals->begin(), vals->end());
	vector<double> uls;
	
	for (j = 1; j < vals->size(); j++) {
		y1 = (*vals)[j-0].second;
		y0 = (*vals)[j-1].second;
		
		if (y1 < y0) std::swap(y1,y0);
		if (y0 <= testSize && testSize <= y1)
			uls.push_back(((*vals)[j].first - (*vals)[j-1].first)/((*vals)[j].second - (*vals)[j-1].second)*(testSize - (*vals)[j-1].second) + (*vals)[j-1].first);
	}
	
	std::sort(uls.begin(),uls.end());
	
	return (uls.size() > 0 ? uls.back() : numeric_limits<double>::quiet_NaN());
} // ulFromPoints()

HypoTestInverterResult *clsPlotter::loadResult(TFile *file, const char *name)
{
	HypoTestInverterResult *result = (HypoTestInverterResult*)file->Get(name);
	if (!result) {
		RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
		result = (HypoTestInverterResult*)wspace->obj(name);
	}
	
	return result;
} // loadResult()

void clsPlotter::print(double cl)
{
	TFile *file = TFile::Open(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(file, fResultName.c_str());
	double ul,x,y;
	double smPl,smMi;
	double testSize = 1.0 - cl;
	vector<pair<double,double> > vals;
	vector<double> uls;
	Int_t j;
	
	resultBkg->UseCLs();
	resultBkg->SetTestSize(testSize);
	
	ul = resultBkg->UpperLimit();
	
	cout << "UL(Bs -> mumu) = " << ul*kBFBsmm << "\t(" << ul << ") @ " << (int)(cl*100.) << " % CL" << endl;
	
	makeObs();
	vals.clear();
	for(j = 0; j < fObs->GetN(); j++) {
		if(fObs->GetPoint(j, x, y) >= 0)
			vals.push_back(make_pair(x,y));
	}
	ul = ulFromPoints(testSize,&vals);
	cout << "UL(Bs -> mumu) = " << ul*kBFBsmm << "\t(" << ul << ") @ " << (int)(cl * 100.) << " % CL" << endl;
	
	makeSMBands();
	vals.clear();
	for(j = 0; j < fG0SM->GetN(); j++) {
		if (fG0SM->GetPoint(j,x,y) >= 0) {
			vals.push_back(make_pair(x, y));
		}
	}
	ul = ulFromPoints(testSize,&vals);
	
	vals.clear();
	for(j = 0; j < fG1SM->GetN(); j++) {
		if (fG1SM->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x, y + fG1SM->GetErrorYhigh(j)));
	}
	smPl = ulFromPoints(testSize,&vals);
	
	vals.clear();
	for(j = 0; j < fG1SM->GetN(); j++) {
		if (fG1SM->GetPoint(j,x,y) >= 0)
			vals.push_back(make_pair(x,y - fG1SM->GetErrorYlow(j)));
	}
	smMi = ulFromPoints(testSize,&vals);
	
	cout << "EXP[UL(Bs->mumu)] = " << ul << "+" << smPl-ul << "-" << ul-smMi << endl;
	
	delete file;
} // print()

void clsPlotter::makeObs(bool reload)
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&f,fResultName.c_str());
	int j,size = result->ArraySize();
	std::vector<unsigned int> index(size);
	double cls;
	
	if (fObs && reload) {
		delete fObs;
		fObs = NULL;
	}
	
	if (!fObs) {
		fObs = new TGraph;
		
		fObs->SetTitle("CLs observed");
		fObs->SetLineWidth(2);
		
		// sequential access
		TMath::SortItr(result->fXValues.begin(), result->fXValues.end(), index.begin(), false);
		
		for (j = 0; j < size; j++) {
		  cls = result->CLs(index[j]);
		  if (cls < 0) cls = 1.0;
		  fObs->SetPoint(j, result->GetXValue(index[j]), cls);
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
	vector<double> cls;
	
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
	
	fG0SM->SetTitle("Expected SM CLs - Median");
	fG1SM->SetTitle("Expected SM CLs #pm 1 #sigma");
	fG2SM->SetTitle("Expected SM CLs #pm 2 #sigma");
	
	// sequential access
	TMath::SortItr(resultBkg->fXValues.begin(), resultBkg->fXValues.end(), indexBkg.begin(), false);
	TMath::SortItr(resultSM->fXValues.begin(), resultSM->fXValues.end(), indexSM.begin(), false);
	
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
		cls.clear();
		for (k = 0; k < vec.size(); k++) {
			tempResult.SetTestStatisticData(vec[k]);
			cls.push_back(tempResult.CLs());
		}
		
		x = const_cast<double*>(&cls[0]);
		TMath::Quantiles(cls.size(), 5, x, q, p, false);
		
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
	
	fG0Bkg->SetTitle("Expected CLs - Median");
	fG1Bkg->SetTitle("Expected CLs #pm 1 #sigma");
	fG2Bkg->SetTitle("Expected CLs #pm 2 #sigma");
	
	// seq access
	TMath::SortItr(result->fXValues.begin(), result->fXValues.end(), index.begin(), false);
	p[0] = ROOT::Math::normal_cdf(-2);
	p[1] = ROOT::Math::normal_cdf(-1);
	p[2] = 0.5;
	p[3] = ROOT::Math::normal_cdf(1);
	p[4] = ROOT::Math::normal_cdf(2);
	
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
		fObs->GetYaxis()->SetTitle("CL_{s}");
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
} // plot()

void plotResult(const char *file, bool SMExp, bool plotObs = true, double cl = 0.95)
{
	clsPlotter plot(file);
	
	plot.setPlotObs(plotObs);
	plot.setSMBands(SMExp);
	plot.plot(cl);
	plot.print();
} // plotResult()
