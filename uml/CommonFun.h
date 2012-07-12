#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace std;

/// options
static string input_name;
static string input_estimates;
static string meth;
static string ch_s = "-1";
static string pdf_toy = "total";
static string pdf_test;
static string tree_name = "bdt";
static string bias_s = "no";
static string cuts_f = "no";
static bool print = false;
static bool simul = false;
static int NExp = 1;
static int ch_i = -1;
static int inputs = 1;
bool input = false, output = false, method = false, channel = false, estimate = false, pdf = false, roomcs = false, pvalue = false, SM = false, bd_const = false, pdf_test_b = false, bias = false, SB = false;

static string channels[5] = {"bs", "bd", "rare", "comb", "total"};

void help() {
  cout << ">>>>>>>>> options for all programs:" << endl;
  cout << "-i #filename \t input (mandatory)" << endl;
  cout << "-t treename" << endl;
  cout << endl;
  cout << ">>>>>>>>> only for pdf_choise:" << endl;
  cout << "-meth {cnc, bdt} \t cut and count OR boosted decision tree input (mandatory)" << endl;
  cout << "-cha {0, 1} \t barrel OR endcap input (mandatory)" << endl;
  cout << "-SM \t SM constraints" << endl;
  cout << "-bd_const \t Bd constrainted to Bs, over all different channels" << endl;
  cout << "-print \t save the fits to gif and pdf" << endl;
  cout << endl;
  cout << ">>>>>>>>> only for toyMC:" << endl;
  cout << "-e #filename \t estimates file (mandatory)" << endl;
  cout << "-nexp # \t number of experiments" << endl;
  cout << "-pdf {bs, bd, rare, comb, total} \t combination of pdf names, for generating (mandatory)" << endl;
  cout << "-test \t fitting pdf, if different from pdf" << endl;
  cout << "-bias [c+,c-,p+,p-]\t biasing rare pdf parameters" << endl;
  cout << "-roomcs \t toy mc with RooMCStudy" << endl;
  cout << "-pvalue \t pvalue with RooStats" << endl;
  cout << endl;
  cout << ">>>>>>>>> only for fitData:" << endl;
  cout << "-SM \t SM constraints" << endl;
  cout << "-bd_const \t Bd constrainted to Bs, over all different channels" << endl;
  cout << "-print \t save the fits to gif and pdf" << endl;
  cout << "-simul # \t simultaneous fit of # channels (default 1)" << endl;
  cout << "-cuts #filename \t file with MVA selections" << endl;
  cout << "-SB \t fit side-bands only" << endl;
  exit(0);
}

void parse_options(int argc, char* argv[]){
  
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
      ch_i = atoi(ch_s.c_str());
      channel = true;
      cout << "channel: " << ch_s << endl;
    }
    if (!strcmp(argv[i],"-print")) {
      cout << "print plots" << endl;
      print = true;
    }
    if (!strcmp(argv[i],"-nexp")) {
      NExp = atoi(argv[i+1]);
      cout << "number of experiments: " << NExp << endl;
    }
    if (!strcmp(argv[i],"-i")) {
      input_name = argv[i+1];
      cout << "input = " << input_name << endl;
      input = true;
    }
    if (!strcmp(argv[i],"-e")) {
      input_estimates = argv[i+1];
      cout << "estimate file = " << input_estimates << endl;
      estimate = true;
    }
    if (!strcmp(argv[i],"-cuts")) {
      cuts_f = argv[i+1];
      cout << "cuts file = " << cuts_f << endl;
    }
    if (!strcmp(argv[i],"-pdf")) {
      pdf_toy = argv[i+1];
      cout << "pdf = " << pdf_toy << endl;
      pdf = true;
    }
    if (!strcmp(argv[i],"-test")) {
      pdf_test = argv[i+1];
      cout << "testing pdf = " << pdf_test << endl;
      pdf_test_b = true;
    }
    if (!strcmp(argv[i],"-bias")) {
      bias_s = argv[i+1];
      cout << "biasing " << bias_s << endl;
      bias = true;
    }
    if (!strcmp(argv[i],"-roomcs")) {
      cout << "using RooMCStudy" << endl;
      roomcs = true;
    }
    if (!strcmp(argv[i],"-pvalue")) {
      cout << "evaluates pvalue" << endl;
      pvalue = true;
    }
    if (!strcmp(argv[i],"-SM")) {
      cout << "SM constraints" << endl;
      SM = true;
    }
    if (!strcmp(argv[i],"-bd_const")) {
      cout << "Bd constrainted" << endl;
      bd_const = true;
    }
    if (!strcmp(argv[i],"-simul")) {
      inputs = atoi(argv[i+1]);
      cout << inputs << " simultaneous fits of " << inputs << " channels" << endl;
      simul = true;
    }
    if (!strcmp(argv[i],"-t")) {
      tree_name = argv[i+1];
      cout << "tree name = " << tree_name << endl;
    }
    if (!strcmp(argv[i],"-SB")) {
      SB = true;
      cout << "fitting only the side-bands" << endl;
    }
    if (!strcmp(argv[i],"-h")) help();
  }
}

void parse_input (string input) {
  string a_meth_s[2] = {"cnc", "bdt"};
  string a_ch_s[2] = {"0", "1"};
  size_t found;
  for (int i = 0; i < 2; i++) {
    found = input.find(a_meth_s[i]);
    if (found!=string::npos) meth = a_meth_s[i];
    cout << "meth = " << meth << endl;
  }
  for (int i = 0; i < 2; i++) {
    found = input.find(a_ch_s[i]);
    if (found!=string::npos) ch_s = a_ch_s[i];
    cout << "channel = " << ch_s << endl;
  }
  found = input.find("SM");
  if (found!=string::npos) {
    SM = true;
    cout << "SM" << endl;
  }
  found = input.find("BdConst");
  if (found!=string::npos) {
    bd_const = true;
    cout << "bd constrained" << endl;
  }
  found = input.find("simul");
  if (found!=string::npos) {
    simul = true;
    size_t found2;
    found2 = input.find_first_of("0123456789");
    ostringstream number;
    number<< input[found2];
    inputs = atoi(number.str().c_str());
    cout << "simultaneous " << inputs << endl;
  }
}
