/*
 *  printCLb.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 25.09.12.
 *
 */

#include <TFile.h>
#include <TEfficiency.h>
#include <RooWorkspace.h>

#include <RooStats/HypoTestResult.h>

#define kVDEF_pvalueObs				"2011:pvalueObs:val"
#define kVDEF_pvalueObsExponent		"2011:pvalueObs:exponent"
#define kVDEF_signObs				"2011:signObs:val"
#define kVDEF_pvalueObsSM			"2011:pvalueObsSM:val"
#define kVDEF_pvalueObsSMExponent	"2011:pvalueObsSM:exponent"
#define kVDEF_signObsSM				"2011:signObsSM:val"
#define kVDEF_pvalueSM_Med			"2011:pvalueSM:val"
#define kVDEF_pvalueSM_ErrHi		"2011:pvalueSM:errHi"
#define kVDEF_pvalueSM_ErrLo		"2011:pvalueSM:errLo"
#define kVDef_pvalueSM_Exponent		"2011:pvalueSM:exponent"
#define kVDEF_signSM_Med			"2011:signSM:val"
#define kVDEF_signSM_ErrHi			"2011:signSM:errHi"
#define kVDEF_signSM_ErrLo			"2011:signSM:errLo"
#define kVDEF_prob3Sigma			"2011:prob3Sigma:val"
#define kVDEF_prob5Sigma			"2011:prob5Sigma:val"

using namespace RooStats;

static HypoTestResult *loadResult(TFile *file, const char *name)
{
	RooWorkspace *wspace;
	HypoTestResult *result = dynamic_cast<HypoTestResult*> (file->Get(name));
	if(!result) {
		wspace = dynamic_cast<RooWorkspace*> (file->Get("wspace"));
		if (wspace) result = dynamic_cast<HypoTestResult*> (wspace->obj(name));
	}
	
	return result;
} // loadResult()

void printCLb(const char *filename, const char *resultsDir = NULL, bool allDigits = true)
{
	using RooStats::PValueToSignificance;
	TFile *file = TFile::Open(filename);
	const char *poi[] = {"mu_s","mu_d"};
	size_t j,k;
	int expo;
	HypoTestResult *resultBkg;
	HypoTestResult *resultSM;
	double p[5];
	double q[5];
	double *x;
	FILE *latexFile = NULL;
	
	if (resultsDir)
		latexFile = fopen(Form("%s/clb.tex",resultsDir), "w");
	
	// Evaluate the total p value and significance
	resultBkg = loadResult(file, "Hybrid_Bkg");
	if (resultBkg) {
		TEfficiency effCalc("calc", "", 2, 0.0, 1.0);
		effCalc.SetStatisticOption(TEfficiency::kFCP);
		effCalc.SetConfidenceLevel(.68);
		effCalc.SetTotalEvents(1, resultBkg->GetNullDistribution()->GetSize());
		effCalc.SetPassedEvents(1, round(resultBkg->CLsplusb()*(double)resultBkg->GetNullDistribution()->GetSize()));
		effCalc.SetTotalEvents(2, resultBkg->GetAltDistribution()->GetSize());
		effCalc.SetPassedEvents(2, round(resultBkg->CLb()*(double)resultBkg->GetAltDistribution()->GetSize()));
		
		cout << Form("p value for BKG: (%e)+(%e)-(%e) corresponding to (%f)+(%f)-(%f) sigmas.",resultBkg->CLsplusb(),effCalc.GetEfficiencyErrorUp(1),effCalc.GetEfficiencyErrorLow(1),PValueToSignificance(resultBkg->CLsplusb()),PValueToSignificance(resultBkg->CLsplusb()-effCalc.GetEfficiencyErrorLow(1)) - PValueToSignificance(resultBkg->CLsplusb()),PValueToSignificance(resultBkg->CLsplusb()) - PValueToSignificance(resultBkg->CLsplusb() + effCalc.GetEfficiencyErrorUp(1))) << endl;
		cout << Form("p value for SM: (%e)+(%e)-(%e) corresponding to (%f)+(%f)-(%f) sigmas.",resultBkg->CLb(),effCalc.GetEfficiencyErrorUp(2),effCalc.GetEfficiencyErrorLow(2),PValueToSignificance(resultBkg->CLb()),PValueToSignificance(resultBkg->CLb()-effCalc.GetEfficiencyErrorLow(2)) - PValueToSignificance(resultBkg->CLb()),PValueToSignificance(resultBkg->CLb()) - PValueToSignificance(resultBkg->CLb() + effCalc.GetEfficiencyErrorUp(2))) << endl;
		
		if (latexFile) {
			if (allDigits) {
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_pvalueObs, resultBkg->CLsplusb());
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_signObs, PValueToSignificance(resultBkg->CLsplusb()));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_pvalueObsSM, resultBkg->CLb());
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", kVDEF_signObsSM, PValueToSignificance(resultBkg->CLb()));
			} else {
				expo = (int)floor(log10(resultBkg->CLsplusb()));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", kVDEF_pvalueObs, resultBkg->CLsplusb()*pow(10., -expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", kVDEF_pvalueObsExponent, expo);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", kVDEF_signObs, PValueToSignificance(resultBkg->CLsplusb()));
				expo = (int)floor(log10(resultBkg->CLb()));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", kVDEF_pvalueObsSM, resultBkg->CLb()*pow(10., -expo));
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", kVDEF_pvalueObsSMExponent, expo);
				fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", kVDEF_signObsSM, PValueToSignificance(resultBkg->CLb()));
			}
		}
	}
	
	for (j = 0; j < sizeof(poi)/sizeof(const char*); j++) {
		
		resultBkg = loadResult(file,Form("Hybrid_CLb_%s",poi[j]));
		resultSM = loadResult(file,Form("Hybrid_CLb_%s_SM",poi[j]));
		
		if (resultBkg) {
			TEfficiency effCalc("calc", "", 1, 0.0, 1.0);
			effCalc.SetStatisticOption(TEfficiency::kFCP);
			effCalc.SetConfidenceLevel(0.68);
			effCalc.SetTotalEvents(1, resultBkg->GetNullDistribution()->GetSize());
			effCalc.SetPassedEvents(1, round(resultBkg->CLsplusb()*(double)resultBkg->GetNullDistribution()->GetSize()));
			
			cout << Form("p value %s for background model: (%e)+(%e)-(%e) corresponding to (%f)+(%f)-(%f) sigmas.", (j == 0 ? "bsmm" : "bdmm"), resultBkg->CLsplusb(), effCalc.GetEfficiencyErrorUp(1), effCalc.GetEfficiencyErrorLow(1),PValueToSignificance(resultBkg->CLsplusb()),PValueToSignificance(resultBkg->CLsplusb()-effCalc.GetEfficiencyErrorLow(1)) - PValueToSignificance(resultBkg->CLsplusb()), PValueToSignificance(resultBkg->CLsplusb())-PValueToSignificance(resultBkg->CLsplusb() + effCalc.GetEfficiencyErrorUp(1))) << endl;
			
			if (latexFile) {
				if (allDigits) {
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_pvalueObs, poi[j]), resultBkg->CLsplusb());
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_signObs, poi[j]), PValueToSignificance(resultBkg->CLsplusb()));
				} else {
					expo = (int)floor(log10(resultBkg->CLsplusb()));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_pvalueObs,poi[j]), resultBkg->CLsplusb()*pow(10., -expo));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_signObs,poi[j]), PValueToSignificance(resultBkg->CLsplusb()));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_pvalueObsExponent,poi[j]), expo);
				}
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
			
			for (k = 0; k < sm_values.size(); k++) {
				tempResult.SetTestStatisticData(sm_values[k]);
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
			
			cout << Form("Exp[p_SM,%s] = %e+%e-%e = %f+%f-%f", poi[j], q[1], q[2] - q[1], q[1] - q[0],PValueToSignificance(q[1]),PValueToSignificance(q[0])-PValueToSignificance(q[1]),PValueToSignificance(q[1])-PValueToSignificance(q[2])) << endl;
			cout << Form("	P(3 sigma, %s) = %f", poi[j], q[3]) << endl;
			cout << Form("	P(5 sigma, %s) = %f", poi[j], q[4]) << endl;
			
			if (latexFile) {
				if (allDigits) {
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_pvalueSM_Med,poi[j]), q[1]);
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_pvalueSM_ErrHi,poi[j]), q[2]-q[1]);
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%e}}}\n", Form("%s:%s",kVDEF_pvalueSM_ErrLo,poi[j]), q[1]-q[0]);
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_signSM_Med,poi[j]), PValueToSignificance(q[1]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_signSM_ErrHi,poi[j]), PValueToSignificance(q[0])-PValueToSignificance(q[1]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_signSM_ErrLo,poi[j]), PValueToSignificance(q[1])-PValueToSignificance(q[2]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_prob3Sigma,poi[j]), q[3]);
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%f}}}\n", Form("%s:%s",kVDEF_prob5Sigma,poi[j]), q[4]);
				} else {
					expo = (int)floor(log10(q[1]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_pvalueSM_Med,poi[j]), q[1]*pow(10., -expo));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_pvalueSM_ErrHi,poi[j]), (q[2]-q[1])*pow(10.,-expo));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_pvalueSM_ErrLo,poi[j]), (q[1]-q[0])*pow(10.,-expo));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDef_pvalueSM_Exponent,poi[j]), expo);
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_signSM_Med,poi[j]), PValueToSignificance(q[1]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_signSM_ErrHi,poi[j]), PValueToSignificance(q[0])-PValueToSignificance(q[1]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%.1f}}}\n", Form("%s:%s",kVDEF_signSM_ErrLo,poi[j]), PValueToSignificance(q[1])-PValueToSignificance(q[2]));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_prob3Sigma,poi[j]), (int)round(q[3]*100.));
					fprintf(latexFile, "\\vdef{%s}	{\\ensuremath{{%d}}}\n", Form("%s:%s",kVDEF_prob5Sigma,poi[j]), (int)round(q[4]*100.));	
				}
			}
		}
	}
	
	if (latexFile) fclose(latexFile);
	delete file;
} // printCLb()
