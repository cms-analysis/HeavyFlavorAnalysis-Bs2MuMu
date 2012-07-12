#ifndef PDF_ANALYSIS_H
#define PDF_ANALYSIS_H

#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

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

using namespace std;
using namespace RooFit;

class pdf_analysis {


public:
  pdf_analysis(bool print, string meth = "bdt", string ch_s = "0", string range = "all", bool SM = false, bool bd_constr = false);
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}

  void set_rad(RooAbsData* rad) {rds_ = rad;}
  RooAbsData* get_rad() {return rds_;}
  
  void define_pdfs();
  void define_bs();
  void define_bd();
  void define_peaking();
  void define_nonpeaking();
  void define_comb();
  void define_signals();
  void define_rare();
  void define_rare2(RooDataHist *data);
  void define_rare3();
  void define_bkg_fractional();
  void define_bkg_extended();
  void define_signalsrare();

  void set_SMratio(double ratio) {ratio_ = ratio;}
  
  string define_pdf_sum(string name);
  void define_total_fractional(); // final pdf with fractional components, and also extended
  void define_total_extended(); // final pdf with all extended components

  void fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error = true, bool hesse = true);
  void print(string output = "", RooWorkspace *ws = 0);
  void set_pdf_constant(string pdf);
  
  string pdf_name;

  bool SM_;
  bool bd_constr_;
  RooRealVar* Mass;

  int channels;
  string range_;
  RooFitResult* RFR;

protected:
  bool print_;
  string meth_;
  string ch_s_;
  RooWorkspace* ws_;
  RooAbsData* rds_;

  double ratio_;

  string input_estimates_;
  vector <double> estimate_bs;
  vector <double> estimate_bd;
  vector <double> estimate_rare;
  vector <double> estimate_comb;
  vector <double> estimate_channel;
  
};

#endif // PDF_ANALYSIS_H
