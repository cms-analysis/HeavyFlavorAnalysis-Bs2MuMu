#ifndef PDF_ANALYSIS_H
#define PDF_ANALYSIS_H

#include <string>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "RooWorkspace.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooUniform.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooSimWSTool.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooExponential.h"
#include "RooSimultaneous.h"
#include "RooProduct.h"
#include "RooGaussModel.h"
#include "RooFFTConvPdf.h"
#include "RooSimWSTool.h"
#include "RooSuperCategory.h"
#include "RooMultiCategory.h"

using namespace std;
using namespace RooFit;

class pdf_analysis {
public:
  pdf_analysis(bool print, string meth = "bdt", string ch_s = "0", string range = "all", int BF = 0, bool SM = false, bool bd_constr = false, bool simul = false, bool simulbdt = false, bool pee_ = false, bool bdt_fit = false);
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}

  void set_rad(RooAbsData* rad) {rds_ = rad;}
  RooAbsData* get_rad() {return rds_;}
  
  void initialize();
  RooHistPdf* define_MassRes_pdf(RooDataSet *rds, string name);
  RooHistPdf* define_bdt_pdf(RooDataSet *rds, string name);

  void define_pdfs();
  void define_bs(int i, int j);
  void define_bd(int i, int j);
  void define_peaking(int i, int j);
  void define_nonpeaking(int i, int j);
  void define_comb(int i, int j);
  void define_signals(int i, int j);
  void define_rare(int i, int j);
  void define_rare2(RooDataHist *data, int i, int j);
  void define_rare3(int i, int j);
  void define_bkg_fractional(int i, int j);
  void define_bkg_extended(int i, int j);
  void define_signalsrare(int i, int j);
  void define_bf(int i, int j);
  void define_constraints(int i, int j);

  void set_SMratio(double ratio) {ratio_ = ratio;}

  double getformulaerror(RooFormulaVar* formula, int vars);
  
  string define_pdf_sum(string name, int i = 0);
  void define_total_fractional(int i, int j); // final pdf with fractional components, and also extended
  void define_total_extended(int i, int j); // final pdf with all extended components

  void define_simul(bool simulbdt);

  void getBFnumbers(string numbers_filename);

  void fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error = true, bool hesse = true);
  void print(RooAbsData *data, string output = "");
  void set_pdf_constant(string pdf);
  void set_rare_normalization(string input, bool extended = false); //set peak fraction parameter to Bu2JpsiK

  string pdf_name;

  int channel;
  int channel_bdt;
  bool SM_;
  bool bd_constr_;
  bool simul_;
  bool simul_bdt_;
  int BF_;
  RooRealVar* Mass;
  RooRealVar* MassRes;
  RooRealVar* bdt;
  RooRealVar* eta;
  RooRealVar* m1eta;
  RooRealVar* m2eta;
  RooRealVar* weight;

  RooCategory* channels_cat;
  RooCategory* bdt_cat;
  RooSuperCategory* super_cat;

  int channels;
  int channels_bdt;
  string range_;
  RooFitResult* RFR;

  int verbosity;

  void simsplit();

  double getErrorHigh(RooRealVar* var);
  double getErrorLow(RooRealVar* var);

  bool old_tree;
  bool pee;
  bool bdt_fit_;
  bool no_legend;

  void gen_and_fit(string pdfname);
  bool print_;

  const char* name(string name, int i, int j = -1);

  bool newcomb_;
  void setSBslope(RooAbsData* sb_data);

protected:
  string meth_;
  string ch_s_;
  int ch_i_;
  RooWorkspace* ws_;
  RooAbsData* rds_;

  double ratio_;

  RooDataHist getRandom_rdh();
  string get_address(string name, string pdf = "", bool channeling = true);

  RooArgSet* obs;

  vector < string > source;

  void parse_external_numbers(string filename);
  void parse_efficiency_numbers(int offset = 0);

  double Bs2MuMu_SM_BF_val;
  double Bs2MuMu_SM_BF_err;
  double Bd2MuMu_SM_BF_val;
  double Bd2MuMu_SM_BF_err;
  double fs_over_fu_val;
  double Jpsi2MuMu_BF_val;
  double Bu2JpsiK_BF_val;
  double fs_over_fu_err;
  double Jpsi2MuMu_BF_err;
  double Bu2JpsiK_BF_err;

  double K_cor_bs;
  double K_cor_bd;

  vector <double> k_unc_bs;
  vector <double> k_unc_bd;

  vector <double> eff_bd_val;
  vector <double> eff_bs_val;
  vector <double> eff_bu_val;
  vector <double> N_bu_val;
  vector <double> eff_bd_err;
  vector <double> eff_bs_err;
  vector <double> eff_bu_err;
  vector <double> N_bu_err;

  vector <double> BF_bs_val;
  vector <double> BF_bd_val;
  vector <double> BF_bs_err;
  vector <double> BF_bd_err;

private:

};

#endif // PDF_ANALYSIS_H
