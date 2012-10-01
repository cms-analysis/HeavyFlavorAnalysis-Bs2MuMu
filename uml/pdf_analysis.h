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

using namespace std;
using namespace RooFit;

class pdf_analysis {
public:
  pdf_analysis(bool print, string meth = "bdt", string ch_s = "0", string range = "all", bool SM = false, bool bd_constr = false, bool simul = false, bool pee_ = false, bool bdt_fit = false);
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}

  void set_rad(RooAbsData* rad) {rds_ = rad;}
  RooAbsData* get_rad() {return rds_;}
  
  void initialize();
  RooHistPdf* define_MassRes_pdf(RooDataSet *rds, string name);
  RooHistPdf* define_bdt_pdf(RooDataSet *rds, string name);

  void define_pdfs();
  void define_bs(int i);
  void define_bd(int i);
  void define_peaking(int i);
  void define_nonpeaking(int i);
  void define_comb(int i);
  void define_signals(int i);
  void define_rare(int i);
  void define_rare2(RooDataHist *data, int i);
  void define_rare3(int i);
  void define_bkg_fractional(int i);
  void define_bkg_extended(int i);
  void define_signalsrare(int i);

  void set_SMratio(double ratio) {ratio_ = ratio;}
  
  string define_pdf_sum(string name, int i = 0);
  void define_total_fractional(int i); // final pdf with fractional components, and also extended
  void define_total_extended(int i); // final pdf with all extended components

  void define_simul();

  void fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error = true, bool hesse = true);
  void print(RooAbsData *data, string output = "");
  void set_pdf_constant(string pdf);
  void set_rare_normalization(string input, bool extended = false); //set peak fraction parameter to Bu2JpsiK

  string pdf_name;

  int channel;
  bool SM_;
  bool bd_constr_;
  bool simul_;
  RooRealVar* Mass;
  RooRealVar* MassRes;
  RooRealVar* bdt;
  RooRealVar* eta;
  RooRealVar* m1eta;
  RooRealVar* m2eta;
  RooRealVar* weight;
  RooCategory* channels_cat;

  int channels;
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

  const char* name(string name, int i);

protected:
  string meth_;
  string ch_s_;
  int ch_i_;
  RooWorkspace* ws_;
  RooAbsData* rds_;

  double ratio_;

  string input_estimates_;
  vector <double> estimate_bs;
  vector <double> estimate_bd;
  vector <double> estimate_rare;
  vector <double> estimate_comb;
  vector <double> estimate_channel;

  RooDataHist getRandom_rdh();
  string get_address_root(string name);

  RooArgSet* obs;

  vector < string > source;

private:

};

#endif // PDF_ANALYSIS_H
