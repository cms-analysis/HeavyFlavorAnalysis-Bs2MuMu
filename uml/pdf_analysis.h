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
#include "RooPlot.h"

using namespace std;
using namespace RooFit;

class pdf_analysis {
public:
  pdf_analysis(RooWorkspace *ws, bool print, string source_name, string meth, string ch_s);
  void set_ws(RooWorkspace *ws) {ws_ = ws;}
  RooWorkspace* get_ws() {return ws_;}
  void fill_inputs(string input);

  void define_pdf(string pdf_name);
  void define_bs();
  void define_bd();
  void define_peaking();
  void define_nonpeaking();
  void define_comb();

  void fit_pdf ();

  void print();
  string pdf_name;
  string rdh_name;

private:
  bool print_;
  string input_f;
  string meth_;
  string ch_s_;
  string source_name_;
  RooWorkspace* ws_;
  vector <string> pdf_primitives_;

};

#endif // PDF_ANALYSIS_H
