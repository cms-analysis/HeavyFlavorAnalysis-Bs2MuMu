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

class pdf_toyMC : public pdf_analysis {
public:
  
  pdf_toyMC(string input_estimates, bool print, string meth, string ch_s, string range = "all");
  void generate(int NExp);
  void fit_pdf (string pdf, RooAbsData* data, int printlevel = -1);
  void fit_pulls();
  void parse_estimate();
  void unset_constant();
  
  TH1D* pull_h_bs;
  TH1D* pull_h_bd;
  
private:
  bool parse(char *cutName, float cut);
  double estimate_bs;
  double estimate_bd;
  double estimate_rare;
  double estimate_comb;

  string input_estimates_;
  
  RooDataSet* pull_rds_bs;
  RooDataSet* pull_rds_bd;
  RooRealVar* pull_bd;
  RooRealVar* pull_bs;

};

#endif	/* PDF_TOYMC_H */

