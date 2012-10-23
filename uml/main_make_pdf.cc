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
  if (!input || !method || (simul && channel)) help();

  /// MC shapes
  if (SM && bd_const) {cout << "please select SM OR bd_const, not both" << endl; return (EXIT_FAILURE);}
  if (pee && strcmp(input_name.c_str(), "new")) {cout << "per event error only with new trees, select -i new" << endl; return EXIT_FAILURE;}

  pdf_analysis ana1(print, meth, ch_s, "all", BF, SM, bd_const, simul, simul_bdt, pee, bdt_fit);
  if (strcmp(input_name.c_str(), "new")) ana1.old_tree = true;
  if (no_legend) ana1.no_legend = true;
  ana1.set_SMratio(0.09);

  ana1.channels = inputs;
  ana1.channels_bdt = inputs_bdt;
  ana1.initialize();
  RooWorkspace *ws = ana1.get_ws();

  /// INPUTS

  static string decays[] = {"SgMc", "BdMc", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgLb2KP", "bgLb2PiP", "bgBd2PiMuNu", "bgBs2KMuNu", "bgLb2PMuNu", "SgData"};
  double weight_i[] = { 1,      1,      1./27.6,   1./1.1,     1./1.6,      1./1.0,    1./2.5,     1./2.5,      1./1.3,    1./2.2,     1./1.0,        1./1.1,       1./1.1,    1};

  int decays_n = sizeof(decays)/sizeof(string);

  string year_s[2] = {"2011", "2012"};
  int years, years_i;
  if (years_opt == "all") years = 2;
  else {
    years = 1;
    years_i = atoi(years_opt.c_str());
  }

  vector <double> cuts_v(2*years, -10);
  if (cuts_f_b) cuts_v = cut_bdt_file();
  TF1* MassRes_0_h = Fit_MassRes("input/2011/small-SgMc.root", cuts, cuts_v, 0);
  TF1* MassRes_2_h = Fit_MassRes("input/2012/small-SgMc.root", cuts, cuts_v, 1);

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

  vector < double > exp_v_0(get_singlerare_normalization("input/2011/anaBmm.plotResults.2011.tex", 0, decays_n));
  vector < double > exp_v_1(get_singlerare_normalization("input/2011/anaBmm.plotResults.2011.tex", 1, decays_n));
  vector < double > exp_v_2(get_singlerare_normalization("input/2012/anaBmm.plotResults.2012.tex", 0, decays_n));
  vector < double > exp_v_3(get_singlerare_normalization("input/2012/anaBmm.plotResults.2012.tex", 1, decays_n));

  for (int i = 0; i < decays_n; i++) {

    decays_treename[i] = decays[i] + "_bdt";
    decays_rdsname[i] = decays[i] + "_rds";
    RooArgList varlist(*m, *MassRes, /*, *eta, *m1eta, *m2eta*/ *bdt, *channel_cat/*, bdt_cat*/, *weight);
    rds_smalltree[i] = new RooDataSet(decays_rdsname[i].c_str(), decays_rdsname[i].c_str(), varlist, "weight");

    for (int yy = 0; yy < years; yy++) {
      int y = (years == 2) ? yy : years_i;
      decays_filename[i] = "input/" + year_s[y] + "/small-" + decays[i] + ".root";
      cout << decays_filename[i] << endl;
      TFile* smalltree_f = new TFile(decays_filename[i].c_str(), "UPDATE");
      TTree* smalltree = (TTree*)smalltree_f->Get(decays_treename[i].c_str());
      TTree* reduced_tree = smalltree->CopyTree(cuts.c_str());
      Double_t m_t, eta_t, m1eta_t, m2eta_t, bdt_t;
      reduced_tree->SetBranchAddress("m",     &m_t);
      reduced_tree->SetBranchAddress("bdt",   &bdt_t);
      reduced_tree->SetBranchAddress("eta",   &eta_t);
      reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
      reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
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
        /// mass resolution
        if (y == 0) MassRes->setVal(MassRes_0_h->Eval(eta_t));
        else if (y == 1) MassRes->setVal(MassRes_2_h->Eval(eta_t));
        /// eta channels
        if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) {
          channel_cat->setIndex(0 + 2*yy);
          if (i != 13 || newcomb) {
            if (cuts_f_b && bdt_t < cuts_v[0 + 2*yy]) continue;
          }
          else {
            if (cuts_f_b && bdt_t < 0.1) continue;
          }
        }
        else {
          channel_cat->setIndex(1 + 2*yy);
          if (i != 13 || newcomb ) {
            if (cuts_f_b && bdt_t < cuts_v[1 + 2*yy]) continue;
          }
          else {
            if (cuts_f_b && bdt_t < 0.1) continue;
          }
        }
        /// bdt channels
        if (bdt_t < 0.1) bdt_cat->setIndex(0);
        else if (bdt_t < 0.18) bdt_cat->setIndex(1);
        else bdt_cat->setIndex(2);

        double weight = 1;
        if (i > 2 && i < 12) {
          if (y == 0) {
            if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) weight = exp_v_0[i] / events_0;
            else weight = exp_v_1[i] / events_1;
          }
          else {
            if ( fabs(m1eta_t) < 1.4 && fabs(m2eta_t) < 1.4) weight = exp_v_2[i] / events_2;
            else weight = exp_v_3[i] / events_3;
          }
        }
        RooArgSet varlist_tmp(*m, *MassRes,/*, *eta, *m1eta, *m2eta*/ *bdt, *channel_cat/*, *bdt_cat*/);
        //weight = weight/1000;
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

  RooDataSet* rds_rare = (RooDataSet*)rds_peak->Clone("rds_rare");
  rds_rare->append(*rds_semi);
  RooAbsData* rad_rare = rds_rare;

  RooAbsData* rad_comb = rds_smalltree[13];

  for (int j = 0; j < inputs; j++) {
    for (int k = 0; k < inputs_bdt; k++) {
      ana1.channel = simul ? j : ch_i;
      ana1.channel_bdt = simul_bdt ? k : ch_bdt_i;

      ana1.define_MassRes_pdf(rds_smalltree[0], "bs");
      ana1.define_MassRes_pdf(rds_smalltree[1], "bd");
      ana1.define_MassRes_pdf(rds_semi, "semi");
      ana1.define_MassRes_pdf(rds_peak, "peak");
//      ana1.define_MassRes_pdf(rds_rare, "rare");
      ana1.define_MassRes_pdf(rds_smalltree[13], "comb");

      if (bdt_fit) {
        ana1.define_bdt_pdf(rds_smalltree[0], "bs");
        ana1.define_bdt_pdf(rds_smalltree[1], "bd");
        ana1.define_bdt_pdf(rds_semi, "semi");
        ana1.define_bdt_pdf(rds_peak, "peak");
//        ana1.define_bdt_pdf(rds_rare, "rare");
        ana1.define_bdt_pdf(rds_smalltree[13], "comb");
      }
    }
  }
  if (newcomb) {
    ana1.newcomb_ = true;
    ana1.setSBslope(rad_comb);
  }
  ana1.define_pdfs();

  if (simul) ana1.define_simul(simul_bdt);

  get_rare_normalization("anaBmm.plotResults.2011.tex", "./input/2011/");
  get_rare_normalization("anaBmm.plotResults.2012.tex", "./input/2012/", 2);
  system("rm input/rare_frac.txt; cat input/2011/rare_frac.txt >> input/rare_frac.txt; cat input/2012/rare_frac.txt >> input/rare_frac.txt;");
  ana1.set_rare_normalization("input/rare_frac.txt");

  /// FITS
  for (int j = 0; j < inputs; j++) {
    for (int k = 0; k < inputs_bdt; k++) {
      ana1.channel = simul ? j : ch_i;
      ana1.channel_bdt = simul_bdt ? k : ch_bdt_i;
      /// bs
      ana1.fit_pdf(ana1.name("bs", j, k), rad_bs, false);

      /// bd
      ana1.fit_pdf(ana1.name("bd", j, k), rad_bd, false);

      /// peak
      ana1.fit_pdf(ana1.name("peak", j, k), rad_peak, false);

      /// semi
      ana1.fit_pdf(ana1.name("semi", j, k), rad_semi, false);

      /// comb
//      ana1.fit_pdf(ana1.name("comb", j, k), rad_comb, false);

      ana1.print_ = false;
      ana1.define_rare3(j, k);
      ana1.fit_pdf(ana1.name("expo3", j, k), rad_rare, false);
      ana1.print_ = print;
    }
  }

  ostringstream inputs_oss; inputs_oss << inputs;
  string output_s = "output/ws_";
  if (simul) output_s += "simul" + inputs_oss.str() + "_" + meth;
  else output_s += "pdf_" + meth + "_" + ch_s;
  if (simul_bdt) output_s += "_simulBdt";
  if (BF) output_s += Form("_BF%d", BF);
  if (SM) output_s += "_SM";
  if (bd_const) output_s += "_BdConst";
  if (bdt_fit) output_s += "_2D";
  if (pee) output_s += "_PEE";
  output_s += ".root";
  TFile * output_f = new TFile(output_s.c_str(), "RECREATE");
  ws->Write();
  MassRes_0_h->Write();
  MassRes_2_h->Write();
  output_f->Close();
  if (simul && !simul_bdt) ws->pdf("pdf_ext_simul")->graphVizTree(Form("sim_%d_pdf.dot", inputs));
  else if (!simul && !simul_bdt) ws->pdf("pdf_ext_total")->graphVizTree(Form("ext_%s.dot", ch_s.c_str()));
  else ws->pdf("pdf_ext_simul")->graphVizTree(Form("prodext_%s.dot", ch_s.c_str()));
  ws->Print();

  if (!simul && !SM && !bd_const) {
    ws->var("N_bs")->setVal(5);
    ws->var("N_bd")->setVal(2);
    ws->var("N_rare")->setVal(18);
    ws->var("N_comb")->setVal(25);
    ana1.gen_and_fit("pdf_ext_total");
  }
  cout << "workspace saved to " << output_s << endl;
  return EXIT_SUCCESS;
}

