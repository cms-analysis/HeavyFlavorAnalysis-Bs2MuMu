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

using namespace RooStats;

class clsPlotter {
	
	public:
		explicit clsPlotter(const char *filename) : fFilename(filename), fResultName("result_mu_s"), fSMBands(false), fDrawLegend(true) {}
		
		void plot(double cl);
	
	public:
		void setSMBands(bool smBands) { fSMBands = smBands; }
		
		bool getSMBands() { return fSMBands; }
	
	private:
		TGraph *makeObs();
		TMultiGraph *makeSMBands();
		TMultiGraph *makeBkgBands();
		HypoTestInverterResult *loadResult(TFile *file, const char *name);
	private:
		std::string fFilename;
		std::string fResultName;
		bool fSMBands;
		bool fDrawLegend;
};

HypoTestInverterResult *clsPlotter::loadResult(TFile *file, const char *name)
{
	HypoTestInverterResult *result = (HypoTestInverterResult*)file->Get(name);
	if (!result) {
		RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
		result = (HypoTestInverterResult*)wspace->obj(name);
	}
	
	return result;
} // loadResult()

TGraph *clsPlotter::makeObs()
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&f,fResultName.c_str());
	TGraph *obs = new TGraph;
	int j,size = result->ArraySize();
	std::vector<unsigned int> index(size);
	
	obs->SetTitle("CLs observed");
	obs->SetLineWidth(2);
	
	// sequential access
	TMath::SortItr(result->fXValues.begin(), result->fXValues.end(), index.begin(), false);
	
	for (j = 0; j < size; j++)
		obs->SetPoint(j, result->GetXValue(index[j]), result->CLs(index[j]));
	
	return obs;
} // makeObs()

TMultiGraph *clsPlotter::makeSMBands()
{
	TFile f(fFilename.c_str());
	HypoTestInverterResult *resultBkg = loadResult(&f, fResultName.c_str());
	HypoTestInverterResult *resultSM = loadResult(&f, Form("%s_SM",fResultName.c_str()));
	TMultiGraph *mg = new TMultiGraph;
	TGraph *g0 = new TGraph;
	TGraphAsymmErrors *g1 = new TGraphAsymmErrors;
	TGraphAsymmErrors *g2 = new TGraphAsymmErrors;
	vector<unsigned int> indexBkg(resultBkg->ArraySize());
	vector<unsigned int> indexSM(resultSM->ArraySize());
	HypoTestResult *hypoBkg,*hypoSM;
	SamplingDistribution *dist;
	double p[5];
	double q[5];
	double *x;
	unsigned int j,k,ixBkg,ixSM;
	vector<double> cls;
	
	g0->SetTitle("Expected SM CLs - Median");
	g1->SetTitle("Expected SM CLs #pm 1 #sigma");
	g2->SetTitle("Expected SM CLs #pm 2 #sigma");
	
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
		
		g0->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		
		g1->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		g1->SetPointEYlow(j, q[2] - q[1]);
		g1->SetPointEYhigh(j, q[3] - q[2]);
		
		g2->SetPoint(j, resultBkg->GetXValue(ixBkg), q[2]);
		g2->SetPointEYlow(j, q[2] - q[0]);
		g2->SetPointEYhigh(j, q[4] - q[2]);
	}
	
	g2->SetFillColor(kYellow);
	mg->Add(g2,"3");
	g1->SetFillColor(kGreen);
	mg->Add(g1,"3");
	g0->SetLineStyle(2);
	g0->SetLineWidth(2);
	mg->Add(g0,"L");
	
	return mg;
} // makeSMBands()

TMultiGraph *clsPlotter::makeBkgBands()
{
	TFile file(fFilename.c_str());
	HypoTestInverterResult *result = loadResult(&file, fResultName.c_str());
	SamplingDistribution *dist;
	TGraph *g0 = new TGraph;
	TGraphAsymmErrors *g1 = new TGraphAsymmErrors;
	TGraphAsymmErrors *g2 = new TGraphAsymmErrors;
	TMultiGraph *mg = new TMultiGraph;
	unsigned int j,ix,nbr = result->ArraySize();
	std::vector<unsigned int> index(nbr);
	double p[5];
	double q[5];
	double *x;
	
	g0->SetTitle("Expected CLs - Median");
	g1->SetTitle("Expected CLs #pm 1 #sigma");
	g2->SetTitle("Expected CLs #pm 2 #sigma");
	
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
		
		g0->SetPoint(j, result->GetXValue(ix), q[2]);
		
		g1->SetPoint(j, result->GetXValue(ix), q[2]);
		g1->SetPointEYlow(j, q[2] - q[1]); // -1 sigma error
		g1->SetPointEYhigh(j, q[3] - q[2]); // +1 sigma error
		
		g2->SetPoint(j, result->GetXValue(ix), q[2]);
		g2->SetPointEYlow(j, q[2] - q[0]);
		g2->SetPointEYhigh(j, q[4] - q[2]);
		
		delete dist;
	}
	
	g2->SetFillColor(kYellow);
	mg->Add(g2,"3");
	g1->SetFillColor(kGreen);
	mg->Add(g1,"3");
	g0->SetLineStyle(2);
	g0->SetLineWidth(2);
	mg->Add(g0,"L");
	
	return mg;
} // makeBkgBands()

void clsPlotter::plot(double cl)
{
	TLegend *leg = NULL;
	TGraph *gobs = makeObs();
	TObject *obj;
	TLine *line;
	TMultiGraph *mg = fSMBands ? makeSMBands() : makeBkgBands();
	int j,nbr;
	
	if (fDrawLegend) {
		leg = new TLegend(0.53,0.63,0.83,0.83,"","NDC");
		leg->AddEntry(gobs,"","L");
		
		nbr = mg->GetListOfGraphs()->GetSize();
		for (j = nbr-1; j>=0; j--)
			if ((obj = mg->GetListOfGraphs()->At(j)) != NULL) leg->AddEntry(obj,"", ((j == nbr-1) ? "L" : "F"));
		
		leg->SetFillColor(kWhite);
		leg->SetLineWidth(0);
		leg->SetLineColor(0);
		leg->SetShadowColor(kWhite);
	}
	
	// draw the graphs
	gobs->Draw("AL");
	gobs->GetXaxis()->SetTitle("BF / BF_{SM}");
	gobs->GetYaxis()->SetTitle("CL_{s}");
	gobs->GetYaxis()->SetTitleOffset(1.0);
	mg->Draw();
	gobs->Draw("same");
	
	// draw confidence level line
	double alpha = 1. - cl;
	line = new TLine(gobs->GetXaxis()->GetXmin(), alpha, gobs->GetXaxis()->GetXmax(), alpha);
	line->SetLineColor(kRed);
	line->Draw();
	
	// draw legend
	leg->Draw();
} // plot()

void plotResult(const char *file, bool SMExp, double cl = 0.95)
{
	clsPlotter plot(file);
	
	plot.setSMBands(SMExp);
	plot.plot(cl);
} // plotResult()
