/*
 *  printCLb.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 25.09.12.
 *
 */

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <RooWorkspace.h>
#include <TGraphAsymmErrors.h>

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

void addEntry(TGraphAsymmErrors *graph, double xpos, int ix, const char *filename)
{
	using RooStats::PValueToSignificance;
	TFile *file = TFile::Open(filename);
	const char *poi[] = {"mu_s","mu_d"};
	size_t k;
	int npoint;
	HypoTestResult *resultBkg;
	HypoTestResult *resultSM;
	double p[5];
	double q[5];
	double *x;
	double sigma,lo,hi;
	
	resultBkg = loadResult(file,Form("Hybrid_CLb_%s",poi[ix]));
	resultSM = loadResult(file,Form("Hybrid_CLb_%s_SM",poi[ix]));
	
	if (resultBkg) {
	  TEfficiency effCalc("calc", "", 1, 0.0, 1.0);
	  effCalc.SetStatisticOption(TEfficiency::kFCP);
	  effCalc.SetConfidenceLevel(0.68);
	  effCalc.SetTotalEvents(1, resultBkg->GetNullDistribution()->GetSize());
	  effCalc.SetPassedEvents(1, round(resultBkg->CLsplusb()*(double)resultBkg->GetNullDistribution()->GetSize()));
	  
	  cout << Form("p value %s for background model: (%e)+(%e)-(%e) corresponding to (%f)+(%f)-(%f) sigmas.", (ix == 0 ? "bsmm" : "bdmm"), resultBkg->CLsplusb(), effCalc.GetEfficiencyErrorUp(1), effCalc.GetEfficiencyErrorLow(1),PValueToSignificance(resultBkg->CLsplusb()),PValueToSignificance(resultBkg->CLsplusb()-effCalc.GetEfficiencyErrorLow(1)) - PValueToSignificance(resultBkg->CLsplusb()), PValueToSignificance(resultBkg->CLsplusb())-PValueToSignificance(resultBkg->CLsplusb() + effCalc.GetEfficiencyErrorUp(1))) << endl;			
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
	  
	  cout << Form("Exp[p_SM,%s] = %e+%e-%e = %f+%f-%f", poi[ix], q[1], q[2] - q[1], q[1] - q[0],PValueToSignificance(q[1]),PValueToSignificance(q[0])-PValueToSignificance(q[1]),PValueToSignificance(q[1])-PValueToSignificance(q[2])) << endl;
	  cout << Form("	P(3 sigma, %s) = %f", poi[ix], q[3]) << endl;
	  cout << Form("	P(5 sigma, %s) = %f", poi[ix], q[4]) << endl;
	  
	  npoint = graph->GetN();
	  sigma = PValueToSignificance(q[1]);
	  cout << Form("Adding point (%.2f, %.3f)",xpos, sigma) << endl;
	  graph->SetPoint(npoint, xpos, sigma);
	  lo = PValueToSignificance(q[1])-PValueToSignificance(q[2]);
	  hi = PValueToSignificance(q[0])-PValueToSignificance(q[1]);
	  
	  if ( !isfinite(hi) )
	    hi = lo; // just very high

	  cout << Form("With Error lo=%.3f, hi=%.3f", lo, hi) << endl;
	  graph->SetPointError( npoint, 0, 0, lo, hi );
	}
	
	delete file;
} // addEntry()

void scanCLb()
{
  TGraphAsymmErrors *graphBsmm = new TGraphAsymmErrors;
  TGraphAsymmErrors *graphBdmm = new TGraphAsymmErrors;
  TCanvas c;

  graphBsmm->SetNameTitle("bsmm_sign","Bsmm Significance window search");
  graphBdmm->SetNameTitle("bdmm_sign","Bdmm Significance window search");

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
  
  graphBsmm->Draw("a");
  graphBdmm->Draw("a");
  

  graphBsmm->GetXaxis()->SetTitle("Window cut / GeV");
  graphBsmm->GetYaxis()->SetTitle("#sigma");

  graphBdmm->GetXaxis()->SetTitle("Window cut / GeV");
  graphBdmm->GetYaxis()->SetTitle("#sigma");

  new TFile("scan.root","update");
  graphBsmm->Write();
  graphBdmm->Write();
} // scanCLb()
