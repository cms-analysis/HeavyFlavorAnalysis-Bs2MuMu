#include <vector>
#include <string>
#include <iostream>

using namespace std;

/// options
static string input_name;
static string input_estimates;
static string meth;
static string ch_s = "-1";
static string pdf_toy = "total";
static bool print = false;
static bool simul = false;
static int NExp = 1;
static int ch_i = -1;
static int inputs = 1;
bool input = false, output = false, method = false, channel = false, estimate = false, pdf = false, roomcs = false, pvalue = false, SM = false, bd_const = false;

static string channels[5] = {"bs", "bd", "rare", "comb", "total"};

void help() {
  cout << ">>>>>>>>> options for all programs:" << endl;
  cout << "-i #filename \t input (mandatory)" << endl;
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
  cout << "-pdf {bs, bd, rare, comb, total} \t combination of pdf names (mandatory)" << endl;
  cout << "-roomcs \t toy mc with RooMCStudy" << endl;
  cout << "-pvalue \t pvalue with RooStats" << endl;
  cout << endl;
  cout << ">>>>>>>>> only for fitData:" << endl;
  cout << "-SM \t SM constraints" << endl;
  cout << "-bd_const \t Bd constrainted to Bs, over all different channels" << endl;
  cout << "-print \t save the fits to gif and pdf" << endl;
  cout << "-simul # \t simultaneous fit of # channels (default 1)" << endl;
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
    if (!strcmp(argv[i],"-pdf")) {
      pdf_toy = argv[i+1];
      cout << "pdf = " << pdf_toy << endl;
      pdf = true;
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
      cout << "simultaneous fit" << endl;
      inputs = atoi(argv[i+1]);
      simul = true;
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
  }
  for (int i = 0; i < 2; i++) {
    found = input.find(a_ch_s[i]);
    if (found!=string::npos) ch_s = a_ch_s[i];
  }
}
