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

using namespace std;
using namespace RooFit;

class pdf_analysis {


public:
  pdf_analysis(bool print, string meth, string ch_s, string range = "all");
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}
  
  void define_pdfs();
  void define_bs();
  void define_bd();
  void define_peaking();
  void define_nonpeaking();
  void define_comb();
  void define_signals();
  void define_rare();
  void define_bkg();
  void define_signalsrare();

  void set_SMconstraint(double ratio) {SM_ = true; ratio_ = ratio;}
  
  string define_pdf_sum(string name);
  void define_total_fractional(); // final pdf with fractional components, and also extended
  void define_total_extended(); // final pdf with all extended components

  void fit_pdf (string pdf, RooAbsData* data, bool extended);
  void print(string output = "");
  void set_pdf_constant(string pdf);
  
  string pdf_name;
  string rdh_name;


protected:
  bool print_;
  string meth_;
  string ch_s_;
  RooWorkspace* ws_;
  RooAbsData* rds_;
  string range_;

  double ratio_;
  bool SM_;
  
};

#endif // PDF_ANALYSIS_H
