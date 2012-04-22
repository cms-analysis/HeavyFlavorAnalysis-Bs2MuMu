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
static string meth;
static string ch_s;
static bool print = false;

void help() {
  cout << "-print \t save the fits to gif" << endl;
  cout << "-i filename \t input (mandatory)" << endl;
  cout << "-o filename \t output (mandatory)" << endl;
  cout << "-meth {cnc, bdt} \t cut and count OR boosted decision tree (mandatory)" << endl;
  cout << "-cha {0, 1} \t barrel OR endcap (mandatory)" << endl;
  exit(0);
}

void parse_options(int argc, char* argv[]){
  bool input = false, output = false, method = false, channel = false;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-meth")) {
      if (!strcmp(argv[i+1],"cnc")) {
        meth = "cnc";
        method = true;
      }
      if (!strcmp(argv[i+1],"bdt")) {
        meth = "bdt";
        method = true;
      }
      cout << "method: " << meth << endl;
    }
    if (!strcmp(argv[i],"-cha")) {
      ch_s = argv[i+1];
      channel = true;
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
  if (!input || !output || !method || !channel) help();
}

int main(int argc, char* argv[]) {

  parse_options(argc, argv);
  RooWorkspace *ws = new RooWorkspace("ws", "workspace");
  RooRealVar* Mass = new RooRealVar("Mass", "Candidate invariant mass", 5.5, 4.9, 5.9, "GeV/c^{2}");
  ws->import(*Mass);

  pdf_analysis ana1(ws, print ,meth, ch_s);
  ana1.define_pdfs();
  
  /// INPUTS
  TFile* input_f = new TFile(input_name.c_str());
  TH1D* Bs_h = (TH1D*)input_f->Get(Form("Bs_%s_chan%s", meth.c_str(), ch_s.c_str()));
  TH1D* Bd_h = (TH1D*)input_f->Get(Form("Bd_%s_chan%s", meth.c_str(), ch_s.c_str()));
  TH1D* Rare_h;
  if (ch_s=="0") Rare_h = (TH1D*)input_f->Get(Form("bRare_%s", meth.c_str()));
  if (ch_s=="1") Rare_h = (TH1D*)input_f->Get(Form("eRare_%s", meth.c_str()));
  ///histo SB real data
  TH1D* data_h = (TH1D*)input_f->Get(Form("hMassWithAllCuts_%s_5_chan%s", meth.c_str(), ch_s.c_str()));
  for (int i = 1; i < Bs_h->GetNbinsX(); i++){
    if (Bs_h->GetBinCenter(i) < 4.90 || Bs_h->GetBinCenter(i) > 5.90) {
      Bs_h->SetBinContent(i, 0.);
      Bs_h->SetBinError(i, 0.);
      Bd_h->SetBinContent(i, 0.);
      Bd_h->SetBinError(i, 0.);
      Rare_h->SetBinContent(i, 0.);
      Rare_h->SetBinError(i, 0.);
      data_h->SetBinContent(i, 0.);
      data_h->SetBinError(i, 0.);
    }
  }
  //// FIXMEEEEE ERRORS!!!
  for (int i = 1; i < Bs_h->GetNbinsX(); i++) {
    Bs_h->SetBinError(i, Bs_h->GetBinContent(i)*0.15);
    Bd_h->SetBinError(i, Bd_h->GetBinContent(i)*0.15);
    Rare_h->SetBinError(i, Rare_h->GetBinContent(i)*0.20);
  }
  TH1D* signals = (TH1D*)Bs_h->Clone("signals");
  signals->Add(Bd_h, 1);

  /// FITS
  /// bs
  RooDataHist* rdh_bs = new RooDataHist("rdh_bs", "dataset template for Bs", *Mass, Bs_h);
  ana1.fit_pdf("bs", rdh_bs, false);
   
  /// bd
  RooDataHist* rdh_bd = new RooDataHist("rdh_bd", "dataset template for Bd", *Mass, Bd_h);
  ana1.fit_pdf("bd", rdh_bd, false);
  
  /// signals
  RooDataHist* rdh_signals = new RooDataHist("rdh_signals", "dataset template for signals", *Mass, signals);
  
  ana1.fit_pdf("signals", rdh_signals, false);
  
  
  /// rare
  RooDataHist* rdh_rare = new RooDataHist("rdh_rare", "dataset template for rare bkg", *Mass, Rare_h);
  ana1.fit_pdf("rare", rdh_rare, false);
  
  ws->Print();
  ws->SaveAs(Form("output/fit_ws_%s_%s.root", meth.c_str(), ch_s.c_str()));

  return 0;
}
