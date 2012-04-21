#include "pdf_analysis.h"

#include <string>
#include <vector>
#include <iostream>

#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

using namespace std;

/// options
static string input_name;
static string output_name;
static string pdf_source;
static string meth;
static string ch_s;
static bool print = false;

void help() {
  cout << "-s {bs, bd, signals, peak, nonpeak, rare, comb, bkg, all} \t choose which data fit (mandatory)" << endl;
  cout << "-print \t save the fits to gif" << endl;
  cout << "-i filename \t input (mandatory)" << endl;
  cout << "-o filename \t output (mandatory)" << endl;
  cout << "-meth {cnc, bdt} \t cut and count OR boosted decision tree" << endl;
  cout << "-cha {0, 1} \t barrel OR endcap" << endl;
  exit(0);
}

void parse_options(int argc, char* argv[]){
  bool pdf = false, input = false, output = false;
  vector <string> pdf_name(9);
  pdf_name[0] = "bs";
  pdf_name[1] = "bd";
  pdf_name[2] = "signals";
  pdf_name[3] = "peak";
  pdf_name[4] = "nonpeak";
  pdf_name[5] = "rare";
  pdf_name[6] = "comb";
  pdf_name[7] = "bkg";
  pdf_name[8] = "all";
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-s")) {
      pdf_source = (argv[i+1]);
      cout << "fit channel " << pdf_source << endl;
      for (vector<string>::size_type i = 0; i < 9; i++) {
        if (pdf_source == pdf_name[i]) {
          pdf = true;
          break;
        }
      }
    }
    if (!strcmp(argv[i],"-meth")) {
      if (!strcmp(argv[i+1],"cnc")) {
        meth = "cnc";
      }
      if (!strcmp(argv[i+1],"bdt")) {
        meth = "bdt";
      }
      cout << "method: " << meth << endl;
    }
    if (!strcmp(argv[i],"-cha")) {
      ch_s = argv[i+1];
      cout << "channel: " << ch_s << endl;
    }
    if (!strcmp(argv[i],"-print")) {
      cout << "print plots" << endl;
      print = true;
    }
    if (!strcmp(argv[i],"-i")) {
      input_name = argv[i+1];
      cout << "input = " << input_name << endl;
      input = true;
    }
    if (!strcmp(argv[i],"-o")) {
      output_name = argv[i+1];
      cout << "output = " << output_name << endl;
      output = true;
    }
    if (!strcmp(argv[i],"-h")) help();
  }
  if (!pdf || !input || !output) help();
}


int main(int argc, char* argv[]) {

  parse_options(argc, argv);
  RooWorkspace *ws = new RooWorkspace("ws", "workspace");
  RooRealVar* Mass = new RooRealVar("Mass", "Candidate invariant mass", 5.5, 4.9, 5.9, "GeV/c^{2}");
  ws->import(*Mass);

  pdf_analysis ana1(ws, print, pdf_source ,meth, ch_s);

  ana1.fill_inputs(input_name);
  ana1.define_pdf(pdf_source);
  ana1.fit_pdf();
  ws->Print();
  ws->SaveAs(Form("output/fit_ws_%s_%s.root", meth.c_str(), ch_s.c_str()));

  return 0;
}
