/*
 *  printCLb.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 25.09.12.
 *
 */

#include <TFile.h>
#include <RooWorkspace.h>

#include <RooStats/HypoTestResult.h>

#define kVDEF_pvalueObs			"default-11:pvalueObs:val"
#define kVDEF_signObs			"default-11:signObs:val"
#define kVDEF_pvalueSM_Med		"default-11:pvalueSM:val"
#define kVDEF_pvalueSM_ErrHi	"default-11:pvalueSM:errHi"
#define kVDEF_pvalueSM_ErrLo	"default-11:pvalueSM:errLo"
#define kVDEF_signSM_Med		"default-11:signSM:val"
#define kVDEF_signSM_ErrHi		"default-11:signSM:errHi"
#define kVDEF_signSM_ErrLo		"default-11:signSM:errLo"
#define kVDEF_prob3Sigma		"default-11:prob3Sigma:val"
#define kVDEF_prob5Sigma		"default-11:prob5Sigma:val"

using namespace RooStats;

static double sigFromP(double pvalue)
{
	return TMath::Sqrt(2.)*TMath::ErfInverse(1.-2.*pvalue);
} // sigFromP()

void printCLb(const char *filename, const char *latexOut = NULL)
{
	TFile *file = TFile::Open(filename);
	RooWorkspace *wspace = (RooWorkspace*)file->Get("wspace");
	HypoTestResult *resultBkg = (RooStats::HypoTestResult*)wspace->obj("HypoTestCalculator_result");
	HypoTestResult *resultSM = (RooStats::HypoTestResult*)wspace->obj("HypoTestCalculator_result_SM");
	size_t j;
	double p[5];
	double q[5];
	double *x;
	FILE *latexFile = NULL;
	
	if (latexOut)
		latexFile = fopen(latexOut, "w");
	
	if (resultBkg) {
		cout << "p value for background model: " << resultBkg->CLsplusb() << " corresponding to " << sigFromP(resultBkg->CLsplusb()) << " sigmas." << endl;
		if (latexFile) {
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_pvalueObs, resultBkg->CLsplusb());
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_signObs, sigFromP(resultBkg->CLsplusb()));
		}
	}
	if (resultSM) {
		HypoTestResult tempResult;
		SamplingDistribution *dist = resultSM->GetAltDistribution();
		SamplingDistribution *pValueDist;
		const vector<Double_t> sm_values = dist->GetSamplingDistribution();
		vector<Double_t> exp_values;
		
		tempResult.SetPValueIsRightTail(resultBkg->GetPValueIsRightTail());
		tempResult.SetBackgroundAsAlt(resultBkg->GetBackGroundIsAlt());
		tempResult.SetNullDistribution(resultBkg->GetNullDistribution());
		tempResult.SetAltDistribution(resultBkg->GetAltDistribution());
		
		for (j = 0; j < sm_values.size(); j++) {
			tempResult.SetTestStatisticData(sm_values[j]);
			exp_values.push_back(tempResult.CLsplusb());
		}
		
		// evaluate the median & 1 sigma band
		p[0] = ROOT::Math::normal_cdf(-1);
		p[1] = 0.5;
		p[2] = ROOT::Math::normal_cdf(1);
		p[3] = 0.5 * (1 - TMath::Erf(3.0/TMath::Sqrt(2))); // 3 sigma observation
		p[4] = 0.5 * (1 - TMath::Erf(5.0/TMath::Sqrt(2))); // 5 sigma observation
		
		x = const_cast<double*>(&exp_values[0]);
		TMath::Quantiles(exp_values.size(), 3, x, q, p, false);
		
		pValueDist = new SamplingDistribution("","",exp_values);
		q[3] = pValueDist->Integral(-RooNumber::infinity(), p[3], kTRUE, kTRUE, kTRUE);
		q[4] = pValueDist->Integral(-RooNumber::infinity(), p[4], kTRUE, kTRUE, kTRUE);
		delete pValueDist;
		
		cout << Form("Exp[p_SM] = %f+%f-%f = %f+%f-%f", q[1], q[2] - q[1], q[1] - q[0],sigFromP(q[1]),sigFromP(q[0])-sigFromP(q[1]),sigFromP(q[1])-sigFromP(q[2])) << endl;
		cout << Form("	P(3 sigma) = %f", q[3]) << endl;
		cout << Form("	P(5 sigma) = %f", q[4]) << endl;
		if (latexFile) {
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_pvalueSM_Med, q[1]);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_pvalueSM_ErrHi, q[2]-q[1]);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_pvalueSM_ErrLo, q[1]-q[0]);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_signSM_Med, sigFromP(q[1]));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_signSM_ErrHi, sigFromP(q[0])-sigFromP(q[1]));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_signSM_ErrLo, sigFromP(q[1])-sigFromP(q[2]));
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_prob3Sigma, q[3]);
			fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", kVDEF_prob5Sigma, q[4]);
		}
	}
	
	if (latexFile) fclose(latexFile);
	delete file;
} // printCLb()
