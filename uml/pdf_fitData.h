#ifndef PDF_FITDATA_H
#define PDF_FITDATA_H

#include <sstream>
#include <iostream>
#include <iomanip>

#include "pdf_analysis.h"

#include "TTree.h"

#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HybridCalculator.h"

class pdf_fitData : public pdf_analysis {
  public:
  pdf_fitData(bool print, int inputs = 1, string input_estimates = "", string input_cuts = "", string meth = "bdt", string range = "all", bool SM = false, bool bd_constr = false, TTree *input_tree = 0, bool simul = true, string ch_s = "0");

    void print();
    void print_each_channel();

    void define_channels();
    void make_dataset();
    void make_pdf_input();
    void make_pdf();
    RooDataSet* global_data;
    RooSimultaneous* simul_pdf;

    void fit_pdf(string pdf, RooAbsData* data, bool extended, bool sumw2error = true, bool hesse = true);
    void fit_pdf();
    void significance(int meth);
    void save();

  protected:
    void parse_estimate();
    bool parse(char *cutName, float cut);
    string input_cuts_;
    bool random;

  private:

    vector < TFile*> ws_file_input;
    vector < RooWorkspace*> ws_input;

    //RooCategory* channels_cat;
    void FillRooDataSet(RooDataSet* dataset, string cuts_f);
    void changeName(RooWorkspace *ws, int str);
    TTree* tree;

    vector <RooAbsPdf*> total_pdf_i;

    void sig_hand();
    void sig_plhc();
    void sig_plhts();
};

#endif // PDF_FITDATA_H
