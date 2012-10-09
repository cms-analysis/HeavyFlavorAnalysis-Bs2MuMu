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
#include "RooGamma.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"

using namespace RooStats;

class pdf_fitData : public pdf_analysis {
  public:
    pdf_fitData(bool print, int inputs = 1, int inputs_bdt = 1, string input_estimates = "", string meth = "bdt", string range = "all", bool SM = false, bool bd_constr = false, TTree *input_tree = 0, bool simul = false, bool simulbdt = false, bool pee_ = false , bool bdt_fit = false , string ch_s = "0", int sig = -1);
    ~pdf_fitData();
    void print();
    void print_each_channel();

    void make_dataset(bool cut_b, vector<double> cut_, TF1* MassRes_f, string cuts);
    void make_pdf_input();
    void make_pdf();

    void BF(string eff_filename, string numbers_filename);

    RooDataSet* global_data;
    RooSimultaneous* simul_pdf;

    void fit_pdf(bool do_not_import = false);
    void significance();
    void save();

    vector <pair <double, double> > eff_bd;
    vector <pair <double, double> > eff_bs;
    vector <pair <double, double> > eff_bu;
    vector <pair <double, double> > N_bu;
    vector <pair <double, double> > BF_bs;
    vector <pair <double, double> > BF_bd;

  protected:

    string input_estimates_;
    vector <double> estimate_bs;
    vector <double> estimate_bd;
    vector <double> estimate_rare;
    vector <double> estimate_comb;
    vector <double> estimate_channel;

    vector <vector <double> > estimate2D_bs;
    vector <vector <double> > estimate2D_bd;
    vector <vector <double> > estimate2D_rare;
    vector <vector <double> > estimate2D_comb;
    vector <vector <double> > estimate2D_channel;

    void parse_estimate();
    bool parse(char *cutName, float cut);
    string input_cuts_;
    bool random;
    int sign;

  private:

    TFile* ws_file_input;
    RooWorkspace* ws_input;

    void FillRooDataSet(RooDataSet* dataset, bool cut_b, vector<double> cut_, TF1* MassRes_f, string cuts);
    void changeName(RooWorkspace *ws, int str);
    TTree* tree;

    Double_t sig_hand();
    void sig_plhc();
    void sig_plhts();
    void sig_hybrid_plhts();
    void sig_hybrid_roplhts();
    void make_prior();
    void make_models();

    void parse_external_numbers(string filename);
    void parse_efficiency_numbers(string filename);

    pair <double, double> fs_over_fu;
    pair <double, double> Jpsi2MuMu_BF;
    pair <double, double> Bu2JpsiK_BF;


};

#endif // PDF_FITDATA_H
