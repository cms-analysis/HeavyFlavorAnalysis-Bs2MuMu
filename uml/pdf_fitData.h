#ifndef PDF_FITDATA_H
#define PDF_FITDATA_H

#include <sstream>

#include "pdf_analysis.h"

#include "TTree.h"

#include "RooCategory.h"
#include "RooSimultaneous.h"

class pdf_fitData : public pdf_analysis {
  public:
    pdf_fitData(bool print, string meth, string ch_s, string range = "all", bool SM = false, bool bd_constr = false, TTree *input_tree = 0, int inputs = 1);

    void print(string output = "");
    void print_each_channel();

    void make_dataset();
    void make_pdf_input();
    void make_pdf();
    RooDataSet* global_data;
    RooSimultaneous* simul_pdf;
    int channels;

  private:

    vector < TFile*> ws_file_input;
    vector < RooWorkspace*> ws_input;

    RooCategory* channel;
    void FillRooDataSet(TTree* tree, RooDataSet* dataset, RooRealVar *Mass, int ch_i);
    void changeName(RooWorkspace *ws, int str);
    TTree* tree;
};

#endif // PDF_FITDATA_H
