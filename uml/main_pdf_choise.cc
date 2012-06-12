#include "CommonFun.h"
#include "pdf_analysis.h"

#include <string>
#include <vector>

#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "TRandom3.h"
#include "TH1D.h"

int main(int argc, char* argv[]) {

  parse_options(argc, argv);
  if (!input || !method || !channel) help();
  
  /// INPUTS
  TFile* input_f = new TFile(input_name.c_str());
  TH1D* Bs_h = (TH1D*)input_f->Get(Form("Bs_%s_chan%s", meth.c_str(), ch_s.c_str()));
  TH1D* Bd_h = (TH1D*)input_f->Get(Form("Bd_%s_chan%s", meth.c_str(), ch_s.c_str()));
  double SM_ratio = Bd_h->Integral() / Bs_h->Integral();
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
  TH1D* signalsrare = (TH1D*)signals->Clone("signalsrare");
  signalsrare->Add(Rare_h, 1);
  TH1D* bkg_sub_h = (TH1D*)data_h->Clone("bkg_sub_h");
  bkg_sub_h->Add(Rare_h, -1.);
  for (int i = 1; i < bkg_sub_h->GetNbinsX(); i++){
    if ( (bkg_sub_h->GetBinCenter(i) >= 5.20 && bkg_sub_h->GetBinCenter(i) <= 5.45) || bkg_sub_h->GetBinCenter(i) < 4.90 || bkg_sub_h->GetBinCenter(i) > 5.90) {
      bkg_sub_h->SetBinContent(i, 0.);
      bkg_sub_h->SetBinError(i, 0.);
    }
  }
/// combinatorial estimation
  double X_comb_estimated =bkg_sub_h->Integral() * (1. / 0.75); 
  int N_comb_estimated = (int)X_comb_estimated;
  if (X_comb_estimated - N_comb_estimated > 0.5) N_comb_estimated++;
  cout << "estimated comb = " << N_comb_estimated << endl;
  TH1D* uniform_histo = new TH1D("uniform_histo", "uniform_histo", signals->GetNbinsX(), signals->GetXaxis()->GetXmin(), signals->GetXaxis()->GetXmax());
  TRandom3 rand0(0);
  double N = rand0.Poisson(N_comb_estimated);
  for (int i = 1; i <= N; i++) {
    TRandom3 rand1(0);
    double n = rand1.Uniform(4.9, 5.9);
    uniform_histo->Fill(n);
  }
  TH1D* total = (TH1D*)signalsrare->Clone("total");
  total->Add(uniform_histo);
  
  /// MC shapes
  pdf_analysis ana1(print, meth, ch_s);
  if (SM && bd_const) {cout << "please select SM OR bd_const, not both" << endl; return (EXIT_SUCCESS);}
  if (SM) ana1.set_SMconstraint(SM_ratio);
  if (bd_const) ana1.set_bdconstraint();
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
  RooDataHist* rdh_signalsrare = new RooDataHist("rdh_signalsrare", "dataset template for signals + rare bkg", *ws->var("Mass"), signalsrare);
  ana1.fit_pdf("signalsrare", rdh_signalsrare, false);
  
  /// comb
  RooDataHist* uniform_rdh = new RooDataHist("uniform_rdh", "uniform_rdh", *ws->var("Mass"), uniform_histo);
  //ana1.fit_pdf("comb", uniform_rdh, false);
  
  ///total
  RooDataHist* total_rdh = new RooDataHist("total_rdh", "total_rdh", *ws->var("Mass")
          ,total
          );
//  total_rdh->add(*rdh_bs);
//  total_rdh->add(*rdh_bd);
//  total_rdh->add(*rdh_rare);
//  total_rdh->add(*uniform_rdh);
  ana1.fit_pdf("total", total_rdh, true);
  //ana1.fit_pdf("all", total_rdh, true);
  
  string output_s = "output/fit_ws_" + meth + "_" + ch_s;
  if (SM) output_s += "_SM";
  if (bd_const) output_s += "_BdConst";
  output_s += ".root";
  ws->SaveAs(output_s.c_str());
  
  cout << "starting Bs " << Bs_h->Integral() << " fit Bs " << ws->var("N_bs")->getVal() << endl;
  if (SM) cout << "starting Bd " << Bd_h->Integral() << " fit Bd constraint to SM" << endl;
  if (!SM) cout << "starting Bd " << Bd_h->Integral() << " fit Bd " << ws->var("N_bd")->getVal() << endl;
  cout << "starting rare " << Rare_h->Integral() << " fit rare " << ws->var("N_rare")->getVal() << endl;
  cout << "starting comb " << uniform_histo->Integral() << " fit comb " << ws->var("N_comb")->getVal() << endl;
    
  return (EXIT_SUCCESS);
}
