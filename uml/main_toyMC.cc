/* 
 * File:   main_toyMC.cc
 * Author: lucamartini
 *
 * Created on 24 aprile 2012, 18.25
 */
#include "pdf_analysis.h"
#include "pdf_toyMC.h"

#include <stdio.h>
#include <stdlib.h>

#include "CommonFun.h"

#include "TFile.h"

#include "RooWorkspace.h"

/*
 * 
 */

int main(int argc, char** argv) {

  parse_options(argc, argv);
  if (!input || !estimate || !pdf) help();
  
  TFile* input_f = new TFile(input_name.c_str());
  RooWorkspace* ws = (RooWorkspace*)input_f->Get("ws");
  ws->Print();

  parse_input(input_name);

  pdf_toyMC toy1(input_estimates, false, meth, ch_s);
  toy1.set_ws(ws);
 // toy1.unset_constant();
  if (roomcs) toy1.mcstudy(NExp, pdf_toy);
  if (pvalue) toy1.pvalue(NExp);
  if (!roomcs && !pvalue) toy1.generate(NExp, pdf_toy);
  
  delete input_f;

  return (EXIT_SUCCESS);
}

