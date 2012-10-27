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
    pdf_fitData(bool print, int inputs = 1, int inputs_bdt = 1, string input_estimates = "", string meth = "bdt", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, bool simul = false, bool simulbdt = false, bool pee_ = false , bool bdt_fit = false , string ch_s = "0", int sig = -1, bool asimov = false, bool syste = false, bool randomsyste = false, int nexp = 10);
    ~pdf_fitData();
    void print();
    void print_each_channel();

    void make_dataset(bool cut_b, vector<double> cut_, string cuts, TTree *tree, int offset = 0);
    void make_pdf_input();
    void make_pdf();

    void parse_systematics(string filename);

    void BF(string eff_filename, string numbers_filename);

    RooDataSet* global_data;
    RooDataHist* global_datahist;
    RooSimultaneous* simul_pdf;

    void fit_pdf(bool do_not_import = false);
    void significance();
    void save();

    double lumi;
    bool random;
    void setnewlumi();
    void setsyst();

    int proof;

    void extract_N_inRanges();

  protected:

    void randomize_constraints(RooWorkspace* ws);

    string input_estimates_;
    vector <double> estimate_bs;
    vector <double> estimate_bd;
    vector <double> estimate_semi;
    vector <double> estimate_comb;

    vector <vector <double> > estimate2D_bs;
    vector <vector <double> > estimate2D_bd;
    vector <vector <double> > estimate2D_semi;
    vector <vector <double> > estimate2D_comb;
    vector <vector <double> > estimate2D_channel;

    void parse_estimate();
    bool parse(char *cutName, float cut);
    string input_cuts_;
    bool asimov_;
    int sign;

    void addsyst();
    bool parse_sys(char *cutName, double cut);

    string input_systematics_;
    vector <double> systematics_bs;
    vector <double> systematics_bd;
    vector <double> systematics_semi;
    vector <double> systematics_comb;
    vector <double> systematics_channel;

    vector <vector <double> > systematics2D_bs;
    vector <vector <double> > systematics2D_bd;
    vector <vector <double> > systematics2D_semi;
    vector <vector <double> > systematics2D_comb;
    vector <vector <double> > systematics2D_channel;

    string pdfname;
    int NExp;

  private:

    TFile* ws_file_input;
    RooWorkspace* ws_input;

    void FillRooDataSet(RooDataSet* dataset, bool cut_b, vector<double> cut_, string cuts, TTree *tree, int offset);
    void changeName(RooWorkspace *ws, int str);

    Double_t sig_hand();
    void sig_plhc();
    void sig_plhts();
    void sig_hybrid_plhts();
    void sig_hybrid_roplhts();
    void make_prior();
    void make_models();

};

#endif // PDF_FITDATA_H
