#include "pdf_analysis.h"

#include <string>
#include <vector>
#include <iostream>

#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "TRandom3.h"

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
  TH1D* bkg_sub_h = (TH1D*)data_h->Clone("bkg_sub_h");
  bkg_sub_h->Add(Rare_h, -1.);
  for (int i = 1; i < bkg_sub_h->GetNbinsX(); i++){
    if ( (bkg_sub_h->GetBinCenter(i) >= 5.20 && bkg_sub_h->GetBinCenter(i) <= 5.45) || bkg_sub_h->GetBinCenter(i) < 4.90 || bkg_sub_h->GetBinCenter(i) > 5.90) {
      bkg_sub_h->SetBinContent(i, 0.);
      bkg_sub_h->SetBinError(i, 0.);
      //bkg_h->SetBinContent(i, 0.);
      //bkg_h->SetBinError(i, 0.);
    }
  }


  /// MC shapes
  pdf_analysis ana1(print ,meth, ch_s);
  RooWorkspace *ws = ana1.get_ws();

  ana1.define_pdfs();
  /// FITS
  /// bs
  RooDataHist* rdh_bs = new RooDataHist("rdh_bs", "dataset template for Bs", *ws->var("Mass"), Bs_h);
  ana1.fit_pdf("bs", rdh_bs, false);
   
  /// bd
  RooDataHist* rdh_bd = new RooDataHist("rdh_bd", "dataset template for Bd", *ws->var("Mass"), Bd_h);
  ana1.fit_pdf("bd", rdh_bd, false);
  
  /// signals
  RooDataHist* rdh_signals = new RooDataHist("rdh_signals", "dataset template for signals", *ws->var("Mass"), signals);
  ana1.fit_pdf("signals", rdh_signals, false);
  
  /// rare
  RooDataHist* rdh_rare = new RooDataHist("rdh_rare", "dataset template for rare bkg", *ws->var("Mass"), Rare_h);
  ana1.fit_pdf("rare", rdh_rare, false);

  /// signals + rare
  RooDataHist* rdh_signalsrare = (RooDataHist*)rdh_signals->Clone("rdh_signalsrare");
  rdh_signalsrare->add(*rdh_rare);
  ana1.fit_pdf("signalsrare", rdh_signalsrare, false);
  
  ///
  /// SB shapes
  pdf_analysis ana1_SB(print ,meth, ch_s, "sb_lo,sb_hi");
  RooWorkspace *ws_SB = ana1_SB.get_ws();

  
  /// comb in sidebands
  ana1_SB.define_comb();
  RooDataHist* comb_SB_rdh = new RooDataHist("comb_SB_rdh", "comb_SB_rdh", *ws->var("Mass"), bkg_sub_h);

  ana1_SB.fit_pdf("comb", comb_SB_rdh, true);
  
  /// combinatorial estimation
  double X_comb_estimated = ws_SB->var("N_comb")->getVal() * (1. / 0.75); 
  int N_comb_estimated = (int)X_comb_estimated;
  if (X_comb_estimated - N_comb_estimated > 0.5) N_comb_estimated++;
  cout << "estimated comb = " << N_comb_estimated << endl;
  TH1D* uniform_histo = new TH1D("uniform_histo", "uniform_histo", signals->GetNbinsX(), 4.9, 5.9);
  for (int i = 1; i <= N_comb_estimated; i++) {
    TRandom3 rand(0);
    double n = rand.Uniform(4.9, 5.9);
    uniform_histo->Fill(n);
  }
  RooDataHist* uniform_rdh = new RooDataHist("uniform_rdh", "uniform_rdh", *ws->var("Mass"), uniform_histo);

  ///total
  RooDataHist* total_rdh = new RooDataHist("total_rdh", "total_rdh", *ws->var("Mass"));
  total_rdh->add(*rdh_bs);
  total_rdh->add(*rdh_bd);
  total_rdh->add(*rdh_rare);
  total_rdh->add(*uniform_rdh);
  ana1.fit_pdf("total", total_rdh, true);
  //ana1.fit_pdf("all", total_rdh, true);
  
  ws->SaveAs(Form("output/fit_ws_%s_%s.root", meth.c_str(), ch_s.c_str()));
  
  cout << "starting Bs " << Bs_h->Integral() << " fit Bs " << ws->var("N_bs")->getVal() << endl;
  cout << "starting Bd " << Bd_h->Integral() << " fit Bd " << ws->var("N_bd")->getVal() << endl;
  cout << "starting rare " << Rare_h->Integral() << " fit rare " << ws->var("N_rare")->getVal() << endl;
  cout << "starting comb " << uniform_histo->Integral() << " fit comb " << ws->var("N_comb")->getVal() << endl;
  
  return 0;
}
