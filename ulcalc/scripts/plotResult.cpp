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
#include <RooWorkspace.h>

#include <RooStats/SamplingDistribution.h>
#include <RooStats/HypoTestInverterResult.h>

using namespace RooStats;

TGraph *plotResult(const char *filename)
{
	TFile *f = new TFile(filename);
	RooWorkspace *wspace = (RooWorkspace*)f->Get("wspace");
	HypoTestInverterResult *result = (HypoTestInverterResult*)wspace->obj("result_mu_s");
	HypoTestInverterResult *resultSM = (HypoTestInverterResult*)wspace->obj("result_mu_s_SM");
	TGraph *gobs = new TGraph; // observed values
	TGraph *g0 = new TGraph;
	TGraphAsymmErrors *g1 = new TGraphAsymmErrors;
	TGraphAsymmErrors *g2 = new TGraphAsymmErrors;
	TMultiGraph *mg = new TMultiGraph;
	int nEntries = result->ArraySize();
	int j;
	double p[5];
	double q[5];
	
	gobs->SetTitle("CLs observed");
	g0->SetTitle("Expected CLs - Median");
	g1->SetTitle("Expected CLs #pm 1 #sigma");
	g2->SetTitle("Expected CLs #pm 2 #sigma");
	
	// access sequentially
	std::vector<unsigned int> index(nEntries);
	TMath::SortItr(result->fXValues.begin(), result->fXValues.end(), index.begin(), false);
	
	p[0] = ROOT::Math::normal_cdf(-2);
	p[1] = ROOT::Math::normal_cdf(-1);
	p[2] = 0.5;
	p[3] = ROOT::Math::normal_cdf(1);
	p[4] = ROOT::Math::normal_cdf(2);
	
	// fill the graph
	for (j = 0; j < nEntries; j++) {
		SamplingDistribution *dist;
		int ix = index[j];
		gobs->SetPoint(j, result->GetXValue(ix), result->CLs(ix));
		
		dist = result->GetExpectedPValueDist(ix);
		const std::vector<double> &values = dist->GetSamplingDistribution();
		
		double *x = const_cast<double*>(&values[0]);
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
	
	gobs->SetLineWidth(2);
	
	// draw graphs
	gobs->Draw("AL");
	mg->Draw();
	gobs->Draw("same");
	
	// draw confidence level line
	double alpha = 1. - result->ConfidenceLevel();
	TLine *line = new TLine(gobs->GetXaxis()->GetXmin(), alpha, gobs->GetXaxis()->GetXmax(), alpha);
	line->SetLineColor(kRed);
	line->Draw();
	
	// draw legend
	TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","NDC");
	leg->AddEntry(gobs,"","L");
	leg->AddEntry(g0,"","L");
	leg->AddEntry(g1,"","F");
	leg->AddEntry(g2,"","F");
	
	leg->Draw();
	
	delete f;
	return gobs;
} // plotResult()
