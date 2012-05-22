/* 
 * File:   pdf_toyMC.h
 * Author: lucamartini
 *
 * Created on 26 aprile 2012, 11.03
 */

#ifndef PDF_TOYMC_H
#define	PDF_TOYMC_H

#include "pdf_analysis.h"

#include "TPaveStats.h"

#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooMCStudy.h"
#include "RooDLLSignificanceMCSModule.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/HybridCalculator.h"

class pdf_toyMC : public pdf_analysis {
public:
  
  pdf_toyMC(string input_estimates, bool print, string meth, string ch_s, string range = "all");
  void generate(int NExp, string pdf_toy);
  void mcstudy(int NExp, string pdf_toy);
  void pvalue(int NExp);
  void fit_pdf (string pdf, RooAbsData* data, int printlevel = -1);
  void fit_pulls();
  void parse_estimate();
  void unset_constant();
  void set_ws(RooWorkspace *ws);
  
  TH1D* pull_h_bs;
  TH1D* pull_h_bd;
  
private:
  bool parse(char *cutName, float cut);
  double estimate_bs;
  double estimate_bd;
  double estimate_rare;
  double estimate_comb;

  string input_estimates_;
  string pdf_toy_;
  
  RooDataSet* pull_rds_bs;
  RooDataSet* pull_rds_bd;
  RooRealVar* pull_bd;
  RooRealVar* pull_bs;

};

#endif	/* PDF_TOYMC_H */

