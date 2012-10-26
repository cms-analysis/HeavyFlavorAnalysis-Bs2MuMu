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
  pdf_toyMC(bool print, int inputs = 1, int inputs_bdt = 1, string input_estimates = "", string meth = "bdt", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, bool simul = false, bool simulbdt = false, bool pee_ = false, bool bdt_fit = false, string ch_s = "0", int sig = -1, bool asimov = false, bool syste = false, bool randomsyste = false, int nexp = 10, string bias = "no");
  ~pdf_toyMC();

  void generate(string pdf_toy, string test_pdf = "total");
  void mcstudy(string pdf_toy, string test_pdf = "total");
  void unset_constant();

private:
  string bias_;

  string pdf_toy_;
  string pdf_test_;
  
  vector <vector <RooDataSet*> > residual_rds_bs;
  vector <vector <RooDataSet*> > residual_rds_bd;
  vector <vector <RooDataSet*> > pull_rds_bs;
  vector <vector <RooDataSet*> > pull_rds_bd;
  vector <vector <RooDataSet*> > pull_rds_semi;
  vector <vector <RooDataSet*> > pull_rds_comb;
  vector <vector <RooRealVar*> > residual_bs;
  vector <vector <RooRealVar*> > residual_bd;
  vector <vector <RooRealVar*> > pull_bd;
  vector <vector <RooRealVar*> > pull_bs;
  vector <vector <RooRealVar*> > pull_semi;
  vector <vector <RooRealVar*> > pull_comb;
  RooRealVar* pull_BF_bs;
  RooDataSet* pull_rds_BF_bs;
  RooRealVar* pull_BF_bd;
  RooDataSet* pull_rds_BF_bd;

  void print(string output, RooWorkspace* ws);

  void fit_pulls(RooRealVar *pull, RooDataSet *rds, int i, int j);
  void print_histos(TH1D* histos, int i, int j);
  RooFitResult* fit_pdf (string pdf, RooAbsData* data, int printlevel = -1, RooWorkspace *ws = 0);
  Double_t sig_hand(RooAbsData *data, int printlevel, RooWorkspace *ws);
  void do_bias(RooWorkspace* ws);

};

#endif	/* PDF_TOYMC_H */

