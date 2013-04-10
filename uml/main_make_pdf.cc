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
#include "TCanvas.h"

#include "RooSimWSTool.h"
#include "RooCategory.h"

int main(int argc, char* argv[]) {

  parse_options(argc, argv);
  if (simul && channel) help();

  /// MC shapes
  if (SM && bd_const) {cout << "please select SM OR bd_const, not both" << endl; return EXIT_FAILURE;}

  pdf_analysis ana1(print, ch_s, "all", BF, SM, bd_const, inputs, (!simul_all) ? inputs_bdt : 1, inputs_all, pee, bdt_fit);
  ana1.set_SMratio(0.09);
  if (no_legend) ana1.no_legend = true;
  RooWorkspace *ws = ana1.get_ws();

  /// INPUTS

  static string decays[] = {"SgMc", "BdMc", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgLb2KP", "bgLb2PiP", "bgBd2PiMuNu", "bgBs2KMuNu", "bgLb2PMuNu", "SgData"};
  int decays_n = sizeof(decays)/sizeof(string);

  string year_s[2] = {"2011", "2012"};
  int years, years_i;
  if (years_opt == "all") years = 2;
  else {
    years = 1;
    years_i = atoi(years_opt.c_str());
  }

  vector <double> cuts_v(2*years, -1);
  if (cuts_f_b) cuts_v = cut_bdt_file();

  vector <string> decays_filename(decays_n);
  vector <string> decays_treename(decays_n);
  vector <string> decays_rdsname(decays_n);

  vector <RooDataSet*> rds_smalltree(decays_n);

  RooRealVar* m = ws->var("Mass");
  RooRealVar* bdt = ws->var("bdt");
  RooRealVar* eta = ws->var("eta");
  RooRealVar* m1eta = ws->var("m1eta");
  RooRealVar* m2eta = ws->var("m2eta");
  RooRealVar* weight = ws->var("weight");
  RooRealVar* MassRes = ws->var("MassRes");
  RooCategory* channel_cat = ws->cat("etacat");
  RooCategory* bdt_cat = ws->cat("bdtcat");
  RooCategory* all_cat = ws->cat("allcat");

  vector < double > exp_v_0(get_singlerare_normalization("input/2011/anaBmm.plotResults.2011.tex", 0, decays_n));
  vector < double > exp_v_1(get_singlerare_normalization("input/2011/anaBmm.plotResults.2011.tex", 1, decays_n));
  vector < double > exp_v_2(get_singlerare_normalization("input/2012/anaBmm.plotResults.2012.tex", 0, decays_n));
  vector < double > exp_v_3(get_singlerare_normalization("input/2012/anaBmm.plotResults.2012.tex", 1, decays_n));

  for (int i = 0; i < decays_n; i++) {

    decays_treename[i] = decays[i] + "_bdt";
    decays_rdsname[i] = decays[i] + "_rds";
    RooArgList varlist(*m, *MassRes, *channel_cat);
    if (bdt_fit) varlist.add(*bdt);
    if (simul_bdt || simul_all) varlist.add(*bdt_cat);
    if (simul_all) varlist.add(*all_cat);
    varlist.add(*weight);
    rds_smalltree[i] = new RooDataSet(decays_rdsname[i].c_str(), decays_rdsname[i].c_str(), varlist, "weight");

    for (int yy = 0; yy < years; yy++) {
      int y = (years == 2) ? yy : years_i;
      decays_filename[i] = "input/" + year_s[y] + "/small-" + decays[i] + ".root";
      cout << decays_filename[i] << endl;
      TFile* smalltree_f = new TFile(decays_filename[i].c_str(), "UPDATE");
      TTree* smalltree = (TTree*)smalltree_f->Get(decays_treename[i].c_str());
      TTree* reduced_tree = smalltree->CopyTree(cuts.c_str()); // string cuts
      Double_t m_t, eta_t, m1eta_t, m2eta_t, bdt_t, me_t;
      reduced_tree->SetBranchAddress("m",     &m_t);
      reduced_tree->SetBranchAddress("bdt",   &bdt_t);
      reduced_tree->SetBranchAddress("eta",   &eta_t);
      reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
      reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
      reduced_tree->SetBranchAddress("me",    &me_t);
      double entries = reduced_tree->GetEntries();
      double events_0 = 0, events_1 = 0, events_2 = 0, events_3 = 0;
      for (int j = 0; j < entries; j++) {
        reduced_tree->GetEntry(j);
        if (y == 0) {
          if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) events_0++;
          else events_1++;
        }
        else {
          if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) events_2++;
          else events_3++;
        }
      }
      for (int j = 0; j < entries; j++) {
        reduced_tree->GetEntry(j);
        m->setVal(m_t);
        eta->setVal(eta_t);
        m1eta->setVal(m1eta_t);
        m2eta->setVal(m2eta_t);
        bdt->setVal(bdt_t);
        if (m_t > 5.20 && m_t < 5.45 && i == 13) continue; // skip signal windows for comb bkg
        if (m_t < 4.90 || m_t > 5.90) continue; // skip outside range
        if (me_t < 0.0 || me_t > 0.2) continue; //skip wrong mass scale
        MassRes->setVal(me_t);
        /// eta channels
        int eta_channel = -1;
        if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) {
          eta_channel = 0 + 2*yy;
          channel_cat->setIndex(eta_channel);
          if (i != 13 || newcomb) {
            if (cuts_f_b && bdt_t < cuts_v[0 + 2*yy]) continue;
          }
          else {
            if (cuts_f_b && bdt_t < 0.1) continue; // give some statistics to comb!
          }
        }
        else {
          eta_channel = 1 + 2*yy;
          channel_cat->setIndex(eta_channel);
          if (i != 13 || newcomb ) {
            if (cuts_f_b && bdt_t < cuts_v[1 + 2*yy]) continue;
          }
          else {
            if (cuts_f_b && bdt_t < 0.1) continue; // give some statistics to comb!
          }
        }
        /// bdt channels
        int bdt_channel = ana1.bdt_index(eta_channel, bdt_t);
        if (simul_bdt || simul_all) {
          if (bdt_channel == -1) continue; /// bdt < 0.1
          bdt_cat->setIndex(bdt_channel);
        }
        if (simul_all) all_cat->setIndex(ana1.super_index(eta_channel, bdt_channel));

        double weight = 1;
        if (i > 1 && i < 13) {
          if (y == 0) {
            if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) weight = exp_v_0[i] / events_0;
            else weight = exp_v_1[i] / events_1;
          }
          else {
            if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) weight = exp_v_2[i] / events_2;
            else weight = exp_v_3[i] / events_3;
          }
        }
        RooArgList varlist_tmp(*m, *MassRes, *channel_cat);
        if (bdt_fit) varlist_tmp.add(*bdt);
        if (simul_bdt || simul_all) varlist_tmp.add(*bdt_cat);
        if (simul_all) varlist_tmp.add(*all_cat);
        rds_smalltree[i]->add(varlist_tmp, weight);
      }
      cout << rds_smalltree[i]->GetName() << " done: " << rds_smalltree[i]->sumEntries() << " <--- " << smalltree->GetEntries() << endl;
    }
  }

  /// datasets
  RooAbsData* rad_bs = rds_smalltree[0];

  RooAbsData* rad_bd = rds_smalltree[1];

  RooDataSet* rds_semi = (RooDataSet*)rds_smalltree[10]->Clone("rds_semi");
  rds_semi->append(*rds_smalltree[11]);
  rds_semi->append(*rds_smalltree[12]);
  RooAbsData* rad_semi = rds_semi;

  RooDataSet* rds_peak = (RooDataSet*)rds_smalltree[2]->Clone("rds_peak");
  rds_peak->append(*rds_smalltree[3]);
  rds_peak->append(*rds_smalltree[4]);
  rds_peak->append(*rds_smalltree[5]);
  rds_peak->append(*rds_smalltree[6]);
  rds_peak->append(*rds_smalltree[7]);
  rds_peak->append(*rds_smalltree[8]);
  rds_peak->append(*rds_smalltree[9]);
  RooAbsData* rad_peak = rds_peak;

  RooAbsData* rad_comb = rds_smalltree[13];

  for (int j = 0; j < inputs; j++) {
    ana1.channel = simul ? j : ch_i;
    /// 1D
    for (int k = 0; k < ana1.bdt_index_max(j); k++) {
      ana1.channel_bdt = (simul_bdt || simul_all) ? k : ch_bdt_i;
      ana1.define_MassRes_pdf(rds_smalltree[0], "bs");
      ana1.define_MassRes_pdf(rds_smalltree[1], "bd");
      ana1.define_MassRes_pdf(rds_semi, "semi");
      ana1.define_MassRes_pdf(rds_peak, "peak");
      ana1.define_MassRes_pdf(rds_smalltree[13], "comb");
    }
    /// 2D
    if (bdt_fit) {
      ana1.define_bdt_pdf(rds_smalltree[0], "bs", 0.);
      ana1.define_bdt_pdf(rds_smalltree[1], "bd", 0.);
      ana1.define_bdt_pdf(rds_semi, "semi", 0.);
      ana1.define_bdt_pdf(rds_peak, "peak", 0.);
      ana1.define_bdt_pdf(rds_smalltree[13], "comb", 0.);
    }
  }
  ana1.define_N();

  get_rare_normalization("anaBmm.plotResults.2011.tex", "./input/2011/");
  get_rare_normalization("anaBmm.plotResults.2012.tex", "./input/2012/", 2);
  system("rm input/rare_frac.txt; cat input/2011/rare_frac.txt >> input/rare_frac.txt; cat input/2012/rare_frac.txt >> input/rare_frac.txt;");
  ana1.set_rare_normalization("input/rare_frac.txt");

  ana1.define_pdfs();
//  if (simul) ana1.define_simul();

  /// FITS
  for (int j = 0; j < inputs; j++) {
    for (int k = 0; k < ana1.bdt_index_max(j); k++) {
      ana1.channel = simul ? j : ch_i;
      ana1.channel_bdt = (simul_bdt || simul_all) ? k : ch_bdt_i;
      /// bs
      ana1.fit_pdf(ana1.name("bs", j, k), rad_bs, false, true, false);

      /// bd
      ana1.fit_pdf(ana1.name("bd", j, k), rad_bd, false, true, false);

      /// peak
      ana1.fit_pdf(ana1.name("peak", j, k), rad_peak, false, true, false);

      /// semi
      ana1.fit_pdf(ana1.name("semi", j, k), rad_semi, false, true, false);

      /// comb
      ana1.setSBslope(ana1.name("comb", j, k), rad_comb);
//        ana1.fit_pdf(ana1.name("comb", j, k), rad_comb, false, true, false, false);

//      ana1.print_ = false;
//      ana1.define_rare3(j, k);
//      ana1.fit_pdf(ana1.name("expo3", j, k), rad_semi, false, true, false);
//      ana1.print_ = print;
    }
  }

  ostringstream inputs_oss; inputs_oss << inputs;
  string output_s = "output/ws_";
  if (simul) output_s += "simul" + inputs_oss.str() + "_" + meth;
  else output_s += "pdf_" + meth + "_" + ch_s;
  if (simul_bdt) output_s += "_simulBdt";
  if (simul_all) output_s += "_simulAll";
  if (BF > 0) output_s += Form("_BF%d", BF);
  if (SM) output_s += "_SM";
  if (bd_const) output_s += "_BdConst";
  if (bdt_fit) output_s += "_2D";
  if (pee) output_s += "_PEE";
  output_s += ".root";
  TFile * output_f = new TFile(output_s.c_str(), "RECREATE");
  ws->Write();
//  MassRes_0_h->Write();
//  MassRes_2_h->Write();
  output_f->Close();
  if (!simul && !simul_bdt) ws->pdf("pdf_ext_total")->graphVizTree(Form("ext_%s.dot", ch_s.c_str()));
  ws->Print();

//  if (!simul && !SM && !bd_const) {
//    ws->var("N_bs")->setVal(5);
//    ws->var("N_bd")->setVal(2);
//    ws->var("N_peak")->setVal(1);
//    ws->var("N_semi")->setVal(1);
//    ws->var("N_comb")->setVal(25);
//    ana1.gen_and_fit("pdf_ext_total");
//  }
  cout << "workspace saved to " << output_s << endl;
  return EXIT_SUCCESS;
}

