/* 
 * File:   pdf_toyMC.h
 * Author: lucamartini
 *
 * Created on 26 aprile 2012, 11.03
 */

#ifndef PDF_TOYMC_H
#define	PDF_TOYMC_H

#include "pdf_fitData.h"

#include "TPaveStats.h"
#include "TH2D.h"

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

class pdf_toyMC : public pdf_fitData {
public:
  
  pdf_toyMC(bool print, int inputs = 1, string input_estimates = "", string input_cuts = "", string meth = "bdt", string range = "all", bool SM = false, bool bd_constr = false, TTree *input_tree = 0, string bias = "no", bool simul = false, bool pee_ = false, bool bdt_fit = false, string ch_s = "0");

  void generate(int NExp, string pdf_toy, string test_pdf = "total");
  void mcstudy(int NExp, string pdf_toy);
  RooFitResult* fit_pdf (string pdf, RooAbsData* data, int printlevel = -1, RooWorkspace *ws = 0);
  void fit_pulls();
  void unset_constant();
  void set_ws(RooWorkspace *ws);
  
  TH1D* pull_h_bs;
  TH1D* pull_h_bd;
  
private:
  string bias_;

  string pdf_toy_;
  string pdf_test_;
  
  vector <RooDataSet*> pull_rds_bs;
  vector <RooDataSet*> pull_rds_bd;
  vector <RooRealVar*> pull_bd;
  vector <RooRealVar*> pull_bs;

  void print(string output, RooWorkspace* ws);

};

#endif	/* PDF_TOYMC_H */

