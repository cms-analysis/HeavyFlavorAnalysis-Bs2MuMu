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
#include "TTree.h"

#include "RooSimWSTool.h"
#include "RooCategory.h"

int main(int argc, char* argv[]) {

  parse_options(argc, argv);
  if (!input || !method || inputs == -1) help();

  /// MC shapes
  if (SM && bd_const) {cout << "please select SM OR bd_const, not both" << endl; return (EXIT_FAILURE);}
  if (pee && strcmp(input_name.c_str(), "new")) {cout << "per event error only with new trees, select -i new" << endl; return EXIT_FAILURE;}

  pdf_analysis ana1(print, meth, "0", "all", SM, bd_const);
  if (pee) ana1.pee = true;
  if (strcmp(input_name.c_str(), "new")) ana1.old_tree = true;
  ana1.initialize();
  RooWorkspace *ws = ana1.get_ws();
  ana1.channels = inputs;
  ana1.simul_ = true;
  ana1.define_pdfs();
  ana1.simsplit();

  /// INPUTS
  vector <RooAbsData*> rad_bs(inputs);
  vector <RooAbsData*> rad_bd(inputs);
  vector <RooAbsData*> rad_signals(inputs);
  vector <RooAbsData*> rad_peak(inputs);
  vector <RooAbsData*> rad_semi(inputs);
  vector <RooAbsData*> rad_rare(inputs);
  vector <RooAbsData*> rad_signalsrare(inputs);
  vector <RooAbsData*> rad_uniform(inputs);
  vector <RooAbsData*> rad_total(inputs);

  vector <RooHistPdf* > eta_rhp_bs(inputs);
  vector <RooHistPdf* > eta_rhp_bd(inputs);
  vector <RooHistPdf* > eta_rhp_signals(inputs);
  vector <RooHistPdf* > eta_rhp_peak(inputs);
  vector <RooHistPdf* > eta_rhp_semi(inputs);
  vector <RooHistPdf* > eta_rhp_rare(inputs);
  vector <RooHistPdf* > eta_rhp_signalsrare(inputs);
  vector <RooHistPdf* > eta_rhp_uniform(inputs);
  vector <RooHistPdf* > eta_rhp_total(inputs);

  if (strcmp(input_name.c_str(), "new")) { // old way
    TFile* input_f = new TFile(input_name.c_str());

    vector <TH1D*> Bs_h(inputs);
    vector <TH1D*> Bd_h(inputs);
    vector <TH1D*> Rare_h(inputs);
    vector <TH1D*> data_h(inputs);
    vector <double> SM_ratio(inputs);

    for (int j = 0; j < inputs; j++) {
      Bs_h[j] = (TH1D*)input_f->Get(Form("Bs_%s_chan%d", meth.c_str(), j));
      Bd_h[j] = (TH1D*)input_f->Get(Form("Bd_%s_chan%d", meth.c_str(), j));
      SM_ratio[j] = Bd_h[j]->Integral() / Bs_h[j]->Integral();
      ana1.set_SMratio(SM_ratio[j]);
      cout << "SM ratio for channel " << j << " = " << SM_ratio[j] << endl;
      if (j == 0) Rare_h[j] = (TH1D*)input_f->Get(Form("bRare_%s", meth.c_str()));
      if (j == 1) Rare_h[j] = (TH1D*)input_f->Get(Form("eRare_%s", meth.c_str()));
      ///histo SB real data
      data_h[j] = (TH1D*)input_f->Get(Form("hMassWithAllCuts_%s_5_chan%d", meth.c_str(), j));
      for (int i = 1; i < Bs_h[j]->GetNbinsX(); i++){
        if (Bs_h[j]->GetBinCenter(i) < 4.90 || Bs_h[j]->GetBinCenter(i) > 5.90) {
          Bs_h[j]->SetBinContent(i, 0.);
          Bs_h[j]->SetBinError(i, 0.);
          Bd_h[j]->SetBinContent(i, 0.);
          Bd_h[j]->SetBinError(i, 0.);
          Rare_h[j]->SetBinContent(i, 0.);
          Rare_h[j]->SetBinError(i, 0.);
          data_h[j]->SetBinContent(i, 0.);
          data_h[j]->SetBinError(i, 0.);
        }
      }
    }
    //// FIXMEEEEE ERRORS!!!
    for (int j = 0; j < inputs; j++) {
      for (int i = 1; i < Bs_h[j]->GetNbinsX(); i++) {
        Bs_h[j]->SetBinError(i, Bs_h[j]->GetBinContent(i)*0.15);
        Bd_h[j]->SetBinError(i, Bd_h[j]->GetBinContent(i)*0.15);
        Rare_h[j]->SetBinError(i, Rare_h[j]->GetBinContent(i)*0.20);
      }
    }
    vector <TH1D*> signals_h(inputs);
    vector <TH1D*> signalsrare(inputs);
    vector <TH1D*> bkg_sub_h(inputs);
    for (int j = 0; j < inputs; j++) {
      signals_h[j] = (TH1D*)Bs_h[j]->Clone("signals");
      signals_h[j]->Add(Bd_h[j], 1);
      signalsrare[j] = (TH1D*)signals_h[j]->Clone("signalsrare");
      signalsrare[j]->Add(Rare_h[j], 1);
      bkg_sub_h[j] = (TH1D*)data_h[j]->Clone("bkg_sub_h");
      bkg_sub_h[j]->Add(Rare_h[j], -1.);
      for (int i = 1; i < bkg_sub_h[j]->GetNbinsX(); i++){
        if ( (bkg_sub_h[j]->GetBinCenter(i) >= 5.20 && bkg_sub_h[j]->GetBinCenter(i) <= 5.45) || bkg_sub_h[j]->GetBinCenter(i) < 4.90 || bkg_sub_h[j]->GetBinCenter(i) > 5.90) {
          bkg_sub_h[j]->SetBinContent(i, 0.);
          bkg_sub_h[j]->SetBinError(i, 0.);
        }
      }
    }
    /// combinatorial estimation
    vector <double> X_comb_estimated(inputs);
    vector <double> N_comb_estimated(inputs);
    vector <TH1D*> uniform_histo(inputs);
    for (int j = 0; j < inputs; j++) {
      X_comb_estimated[j] = bkg_sub_h[j]->Integral() * (1. / 0.75);
      N_comb_estimated[j] = (int)X_comb_estimated[j];
      if (X_comb_estimated[j] - N_comb_estimated[j] > 0.5) N_comb_estimated[j]++;
      cout << "estimated comb = " << N_comb_estimated[j] << endl;
      uniform_histo[j] = new TH1D(Form("uniform_histo_%d", j), "uniform_histo", signals_h[j]->GetNbinsX(), signals_h[j]->GetXaxis()->GetXmin(), signals_h[j]->GetXaxis()->GetXmax());
      TRandom3 rand0(0);
      double N = rand0.Poisson(N_comb_estimated[j]);
      for (int i = 1; i <= N; i++) {
        TRandom3 rand1(0);
        double n = rand1.Uniform(4.9, 5.9);
        uniform_histo[j]->Fill(n);
      }
    }
    vector <TH1D*> total(inputs);
    for (int j = 0; j < inputs; j++) {
      total[j] = (TH1D*)signalsrare[j]->Clone("total");
      total[j]->Add(uniform_histo[j]);
    }
    /// FITS
    vector <RooDataHist*> rdh_bs(inputs);
    vector <RooDataHist*> rdh_bd(inputs);
    vector <RooDataHist*> rdh_signals(inputs);
    vector <RooDataHist*> rdh_rare(inputs);
    vector <RooDataHist*> rdh_signalsrare(inputs);
    vector <RooDataHist*> rdh_uniform(inputs);
    vector <RooDataHist*> rdh_total(inputs);

    for (int j = 0; j < inputs; j++) {
      rdh_bs[j] = new RooDataHist("rdh_bs", "dataset template for Bs", *ws->var("Mass"), Bs_h[j]);
      rdh_bd[j] = new RooDataHist("rdh_bd", "dataset template for Bd", *ws->var("Mass"), Bd_h[j]);
      rdh_rare[j] = new RooDataHist("rdh_rare", "dataset template for rare bkg", *ws->var("Mass"), Rare_h[j]);
      rdh_uniform[j] = new RooDataHist("rdh_uniform", "rdh_uniform", *ws->var("Mass"), uniform_histo[j]);
      rdh_total[j] = new RooDataHist("rdh_total", "rdh_total", *ws->var("Mass"), total[j]);

      rad_bs[j] = rdh_bs[j];
      rad_bd[j] = rdh_bd[j];
      rad_signals[j] = rdh_signals[j];
      rad_rare[j] = rdh_rare[j];
      rad_signalsrare[j] = rdh_signalsrare[j];
      rad_uniform[j] = rdh_uniform[j];
      rad_total[j] = rdh_total[j];
    }
  }
  else { //new small trees way

    string decays[] = {"SgMc", "BdMc", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgLb2KP", "bgLb2PiP", "bgBd2PiMuNu", "bgBs2KMuNu", "bgLb2PMuNu"};
    double weight_i[] = { 1,      1,      1./27.6,   1./1.1,     1./1.6,      1./1.0,    1./2.5,     1./2.5,      1./1.3,    1./2.2,     1./1.0,        1./1.1,       1./1.1};
    // misid:
    double pion_misid    = 0.0015;
    double kaons_misid   = 0.0017;
    double protons_misid = 0.0005;
    weight_i[2]  *= kaons_misid*kaons_misid;
    weight_i[3]  *= kaons_misid*pion_misid;
    weight_i[4]  *= pion_misid*pion_misid;
    weight_i[5]  *= kaons_misid*kaons_misid;
    weight_i[6]  *= kaons_misid*pion_misid;
    weight_i[7]  *= pion_misid*pion_misid;
    weight_i[8]  *= kaons_misid*protons_misid;
    weight_i[9]  *= pion_misid*protons_misid;
    weight_i[10] *= pion_misid;
    weight_i[11] *= kaons_misid;
    weight_i[12] *= protons_misid;

    int decays_n = sizeof(decays)/sizeof(string);

    for (int j = 0; j < inputs; j++) {

      string cut = get_cut(j);

      vector <string> decays_filename(decays_n);
      vector <string> decays_treename(decays_n);
      vector <string> decays_rdsname(decays_n);

      vector <TFile*> smalltree_f(decays_n);
      vector <TTree*> smalltree(decays_n);
      vector <RooDataSet*> rds_smalltree(decays_n);

      RooRealVar* m = ws->var("Mass");
      RooRealVar* eta = ws->var("eta");
      RooRealVar* weight = new RooRealVar("weight", "weight", 0., 100000., "");
      m->SetName("m");

      for (int i = 0; i < decays_n; i++) {
        decays_filename[i] = "input/small-" + decays[i] + ".root";
        decays_treename[i] = decays[i] + "_bdt";
        decays_rdsname[i] = decays[i] + "_rds";

        smalltree_f[i] = new TFile(decays_filename[i].c_str(), "UPDATE");
        smalltree[i] = (TTree*)smalltree_f[i]->Get(decays_treename[i].c_str());
        TTree* reduced_tree = smalltree[i]->CopyTree(cut.c_str());
        reduced_tree->Branch("weight", &weight_i[i], "weight/D");
        rds_smalltree[i] = new RooDataSet(decays_rdsname[i].c_str(), decays_rdsname[i].c_str(), RooArgSet(*m, *eta, *weight), Import(*reduced_tree), WeightVar("weight"));
        rds_smalltree[i]->changeObservableName("m", "Mass");
        cout << rds_smalltree[i]->GetName() << " done: " << reduced_tree->GetEntries() << " / " << smalltree[i]->GetEntries() << endl;
      }
      m->SetName("Mass");

      rad_bs[j] = rds_smalltree[0];
      eta_rhp_bs[j] = ana1.define_etapdf(rds_smalltree[0], Form("bs_channel_%d", j));

      rad_bd[j] = rds_smalltree[1];
      eta_rhp_bd[j] = ana1.define_etapdf(rds_smalltree[1], Form("bd_channel_%d", j));

      RooDataSet* rds_signals = (RooDataSet*)rds_smalltree[0]->Clone("rds_signals");
      rds_signals->append(*rds_smalltree[1]);
      eta_rhp_signals[j] = ana1.define_etapdf(rds_signals, Form("signals_channel_%d", j));
      rad_signals[j] = rds_signals;

      RooDataSet* rds_semi = (RooDataSet*)rds_smalltree[10]->Clone("rds_semi");
      rds_semi->append(*rds_smalltree[11]);
      rds_semi->append(*rds_smalltree[12]);
      eta_rhp_semi[j] = ana1.define_etapdf(rds_semi, Form("semi_channel_%d", j));
      rad_semi[j] = rds_semi;

      RooDataSet* rds_peak = (RooDataSet*)rds_smalltree[2]->Clone("rds_peak");
      rds_peak->append(*rds_smalltree[3]);
      rds_peak->append(*rds_smalltree[4]);
      rds_peak->append(*rds_smalltree[5]);
      rds_peak->append(*rds_smalltree[6]);
      rds_peak->append(*rds_smalltree[7]);
      rds_peak->append(*rds_smalltree[8]);
      rds_peak->append(*rds_smalltree[9]);
      eta_rhp_peak[j] = ana1.define_etapdf(rds_peak, Form("peak_channel_%d", j));
      rad_peak[j] = rds_peak;

      RooDataSet* rds_rare = (RooDataSet*)rds_peak->Clone("rds_rare");
      rds_rare->append(*rds_semi);
      eta_rhp_rare[j] = ana1.define_etapdf(rds_rare, Form("rare_channel_%d", j));
      rad_rare[j] = rds_rare;

      RooDataSet* rds_signalsrare = (RooDataSet*)rds_signals->Clone("rds_signalsrare");
      rds_signalsrare->append(*rds_rare);
      eta_rhp_signalsrare[j] = ana1.define_etapdf(rds_signalsrare, Form("signalsrare_channel_%d", j));
      rad_signalsrare[j] = rds_signalsrare;

      //no uniform without estimation

      //    rad_uniform = rdh_uniform;
      //    rad_total = rdh_total;
    }
  }
//  RooWorkspace * ws_temp = ana1.get_ws();
//  ws_temp->Print();
//  RooArgSet data = ws_temp->allPdfs();
//  TIterator* it = data.createIterator();
//  RooAbsPdf* var_Obj = 0;
//  RooPlot *frame = ws_temp->var("eta")->frame();
//  int i =0;
//  while((var_Obj = (RooAbsPdf*)it->Next())){
//    var_Obj->plotOn(frame, LineColor(kBlue+i));
//    i++;
//    cout << i << " " << var_Obj->GetName() << endl;
//  }
//  TCanvas c("c", "c", 600, 600);
//  frame->Draw();
//  c.Print("fig/all_eta.gif");
//  delete frame;
//  return 0;

  /// FITS
  for (int j = 0; j < inputs; j++) {
    ana1.channel = j;
  /// bs
    ana1.fit_pdf(Form("bs_channel_%d", j), rad_bs[j], false);

  /// bd
    ana1.fit_pdf(Form("bd_channel_%d", j), rad_bd[j], false);

    if (!strcmp(input_name.c_str(), "new")) {
  /// peak
      ana1.fit_pdf(Form("peaking_channel_%d", j), rad_peak[j], false);

  /// semi
      ana1.fit_pdf(Form("nonpeaking_channel_%d", j), rad_semi[j], false);
    }

  /// rare
    ana1.fit_pdf(Form("rare_channel_%d", j), rad_rare[j], false);

    if (strcmp(input_name.c_str(), "new")) {
  /// comb
    //ana1.fit_pdf(Form("comb_channel_%d", j), rad_uniform, false);

  /// total
    ana1.fit_pdf(Form("total_channel_%d", j), rad_total[j], true);
    }
  }

  string output_s = "output/ws_simul_" + meth;
  if (SM) output_s += "_SM";
  if (bd_const) output_s += "_BdConst";
  output_s += ".root";
  ws->SaveAs(output_s.c_str());
  ws->pdf("pdf_ext_simul")->graphVizTree(Form("sim_%d_pdf.dot", inputs));
ws->Print();
  return EXIT_SUCCESS;
}

