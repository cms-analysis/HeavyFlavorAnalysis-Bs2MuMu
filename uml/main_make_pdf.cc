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

  pdf_analysis ana1(print, meth, ch_s, "all", SM, bd_const, simul, pee, bdt_fit);
  if (strcmp(input_name.c_str(), "new")) ana1.old_tree = true;
  if (no_legend) ana1.no_legend = true;
  ana1.set_SMratio(0.09);
  ana1.initialize();
  RooWorkspace *ws = ana1.get_ws();

  ana1.channels = inputs;
  ana1.channels_bdt = inputs_bdt;

  /// INPUTS
  vector <RooAbsData*> rad_bs(inputs);
  vector <RooAbsData*> rad_bd(inputs);
  vector <RooAbsData*> rad_peak(inputs);
  vector <RooAbsData*> rad_semi(inputs);
  vector <RooAbsData*> rad_rare(inputs);

  static string decays[] = {"SgMc", "BdMc", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgLb2KP", "bgLb2PiP", "bgBd2PiMuNu", "bgBs2KMuNu", "bgLb2PMuNu", "SgData"};
  double weight_i[] = { 1,      1,      1./27.6,   1./1.1,     1./1.6,      1./1.0,    1./2.5,     1./2.5,      1./1.3,    1./2.2,     1./1.0,        1./1.1,       1./1.1,    1};
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

  vector <double> cuts_v(2, -10);
  if (cuts_f_b) cuts_v = cut_bdt_file();
  TF1* MassRes_h = Fit_MassRes("input/small-SgMc.root", cuts_b ? cuts : "", cuts_v);

  for (int j = 0; j < inputs; j++) {

    vector <string> decays_filename(decays_n);
    vector <string> decays_treename(decays_n);
    vector <string> decays_rdsname(decays_n);

    vector <TFile*> smalltree_f(decays_n);
    vector <TTree*> smalltree(decays_n);
    vector <RooDataSet*> rds_smalltree(decays_n);

    RooRealVar* m = ws->var("Mass");
    RooRealVar* bdt = ws->var("bdt");
    RooRealVar* eta = ws->var("eta");
    RooRealVar* m1eta = ws->var("m1eta");
    RooRealVar* m2eta = ws->var("m2eta");
    RooRealVar* weight = ws->var("weight");
    RooRealVar* MassRes = ws->var("MassRes");
    RooCategory* channel_cat = ws->cat("channels");
    RooCategory* bdt_cat = ws->cat("bdtcat");

    for (int i = 0; i < decays_n; i++) {
      decays_filename[i] = "input/small-" + decays[i] + ".root";
      decays_treename[i] = decays[i] + "_bdt";
      decays_rdsname[i] = decays[i] + "_rds";

      smalltree_f[i] = new TFile(decays_filename[i].c_str(), "UPDATE");
      smalltree[i] = (TTree*)smalltree_f[i]->Get(decays_treename[i].c_str());
      TTree* reduced_tree = smalltree[i]->CopyTree(cuts_b ? cuts.c_str() : "");
      Double_t m_t, eta_t, m1eta_t, m2eta_t, bdt_t;
      reduced_tree->SetBranchAddress("m",     &m_t);
      reduced_tree->SetBranchAddress("bdt",   &bdt_t);
      reduced_tree->SetBranchAddress("eta",   &eta_t);
      reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
      reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
      RooArgList varlist(*m, *MassRes, *eta, *m1eta, *m2eta, *bdt, *channel_cat, *bdt_cat, *weight);
      rds_smalltree[i] = new RooDataSet(decays_rdsname[i].c_str(), decays_rdsname[i].c_str(), varlist, "weight");
      for (int j = 0; j<reduced_tree->GetEntries(); j++) {
        reduced_tree->GetEntry(j);
        m->setVal(m_t);
        eta->setVal(eta_t);
        m1eta->setVal(m1eta_t);
        m2eta->setVal(m2eta_t);
        bdt->setVal(bdt_t);
        /// mass resolution
        MassRes->setVal(MassRes_h->Eval(eta_t));
        /// eta channels
        if (fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4) {
          channel_cat->setIndex(0);
          if (cuts_f_b && bdt_t < cuts_v[0]) continue;
        }
        else {
          channel_cat->setIndex(1);
          if (cuts_f_b && bdt_t < cuts_v[1]) continue;
        }
        /// bdt channels
        if (bdt_t < 0.1) bdt_cat->setIndex(0);
        else if (bdt_t < 0.18) bdt_cat->setIndex(1);
        else bdt_cat->setIndex(2);

        RooArgSet varlist_tmp(*m, *MassRes, *eta, *m1eta, *m2eta, *bdt, *channel_cat, *bdt_cat);
        rds_smalltree[i]->add(varlist_tmp, weight_i[i]);
      }
      cout << rds_smalltree[i]->GetName() << " done: " << rds_smalltree[i]->sumEntries() << " / " << smalltree[i]->GetEntries() << endl;
    }

    ana1.channel = simul ? j : ch_i;
    rad_bs[j] = rds_smalltree[0];
    ana1.define_MassRes_pdf(rds_smalltree[0], "bs");
    ana1.define_bdt_pdf(rds_smalltree[0], "bs");

    rad_bd[j] = rds_smalltree[1];
    ana1.define_MassRes_pdf(rds_smalltree[1], "bd");
    ana1.define_bdt_pdf(rds_smalltree[1], "bd");

    RooDataSet* rds_semi = (RooDataSet*)rds_smalltree[10]->Clone("rds_semi");
    rds_semi->append(*rds_smalltree[11]);
    rds_semi->append(*rds_smalltree[12]);
    rad_semi[j] = rds_semi;
    ana1.define_MassRes_pdf(rds_semi, "semi");
    ana1.define_bdt_pdf(rds_semi, "semi");

    RooDataSet* rds_peak = (RooDataSet*)rds_smalltree[2]->Clone("rds_peak");
    rds_peak->append(*rds_smalltree[3]);
    rds_peak->append(*rds_smalltree[4]);
    rds_peak->append(*rds_smalltree[5]);
    rds_peak->append(*rds_smalltree[6]);
    rds_peak->append(*rds_smalltree[7]);
    rds_peak->append(*rds_smalltree[8]);
    rds_peak->append(*rds_smalltree[9]);
    rad_peak[j] = rds_peak;
    ana1.define_MassRes_pdf(rds_peak, "peak");
    ana1.define_bdt_pdf(rds_peak, "peak");

    RooDataSet* rds_rare = (RooDataSet*)rds_peak->Clone("rds_rare");
    rds_rare->append(*rds_semi);
    rad_rare[j] = rds_rare;
    ana1.define_MassRes_pdf(rds_rare, "rare");
    ana1.define_bdt_pdf(rds_rare, "rare");

    ana1.define_MassRes_pdf(rds_smalltree[13], "comb");
    ana1.define_bdt_pdf(rds_smalltree[13], "comb");
  }

  ana1.define_pdfs();

  if (strcmp(rare_f.c_str(),"no")) {
    get_rare_normalization("./input/anaBmm.plotResults.default-11.tex");
    ana1.set_rare_normalization(rare_f);
  }

  if (simul) ana1.define_simul(simulbdt);

  /// FITS
  for (int j = 0; j < inputs; j++) {
    ana1.channel = simul ? j : ch_i;
  /// bs
    ana1.fit_pdf(ana1.name("bs", j), rad_bs[j], false);

  /// bd
    ana1.fit_pdf(ana1.name("bd", j), rad_bd[j], false);

    /// peak
    ana1.fit_pdf(ana1.name("peak", j), rad_peak[j], false);

  /// semi
      ana1.fit_pdf(ana1.name("semi", j), rad_semi[j], false);

  /// rare
    if (!strcmp(rare_f.c_str(),"no")) ana1.fit_pdf(ana1.name("rare", j), rad_rare[j], false);

    ana1.print_ = false;
    ana1.define_rare3(j);
    ana1.fit_pdf(ana1.name("expo3", j), rad_rare[j], false);
    ana1.print_ = print;
  }

  string output_s = "output/ws_";
  if (simul) output_s += "simul_" + meth;
  else output_s += "pdf_" + meth + "_" + ch_s;
  if (SM) output_s += "_SM";
  if (bd_const) output_s += "_BdConst";
  if (bdt_fit) output_s += "_2D";
  if (pee) output_s += "_PEE";
  output_s += ".root";
  ws->SaveAs(output_s.c_str());
  if (simul) ws->pdf("pdf_ext_simul")->graphVizTree(Form("sim_%d_pdf.dot", inputs));
  else ws->pdf("pdf_ext_total")->graphVizTree(Form("ext_%s.dot", ch_s.c_str()));
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

