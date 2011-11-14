#include "plotIso.hh"

#include "../macros/AnalysisDistribution.hh"
#include "TMath.h"

using namespace std; 
using std::string; 

ClassImp(plotIso)

// ----------------------------------------------------------------------
plotIso::plotIso(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

}


// ----------------------------------------------------------------------
plotIso::~plotIso() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotIso::makeAll(int channels) {
  
  fMode = 0; 

}


// ----------------------------------------------------------------------
void plotIso::dataVsMc(string file1, string dir1, string file2, string dir2, const char *selection) {

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  char option1[100], option2[100],  loption1[100], loption2[100]; 
  int color1(1), fill1(1000), color2(1), fill2(1000), marker1(20), marker2(25); 

  if (string::npos != file1.find("No")) {
    color1 = kBlack; 
    fill1  = 3004; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("No")) {
    color2 = kBlue; 
    fill2  = 3004; 
    fill2  = 3356; 
  }
  if (string::npos != file1.find("Cs")) {
    color1 = kBlack; 
    fill1  = 3005; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("Cs")) {
    color2 = kRed; 
    fill2  = 3005; 
    fill2  = 3356; 
  }

  if (string::npos != file1.find("Sg")) {
    // -- data is first, then signal MC
    color1 = kBlack; 
    fill1  = 0; 
    sprintf(option1, "e"); 
    sprintf(loption1, "p");
    color2 = kBlue; 
    fill2  = 3005; 
    fill2  = 3356; 
    sprintf(option2, "hist"); 
    sprintf(loption2, "f");
  }
  
  // -- now fix overlays of two data files
  sprintf(option1, "e");
  sprintf(loption1, "p");
  
  sprintf(option2,  "hist");
  sprintf(loption2, "f");
  
  if (string::npos != file2.find("Data")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(option2, "e"); 
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 21; 
  }
  
  if (string::npos != file1.find("Mc")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 25; 
  }
 
  //  ofstream OUT("testUL.txt", ios::app);

  fF[file1]->cd(dir1.c_str()); 
  TH1D *h = (TH1D*)gDirectory->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }

  TH1D *h1(0), *h2(0);
  AnalysisDistribution a("A_pvz"); 
  //  a.fVerbose = 1; 

  TCanvas *c1;
  string cut, pdfname; 
  vector<string> doList; 

  for (unsigned int i = 0; i < doList.size(); ++i) {
    cut = Form("%s", doList[i].c_str()); 

    pdfname = Form("%s/%s_%s-%s_sbs_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), cut.c_str(), selection);
    
    cout << pdfname << endl;
    
    fF[file1]->cd(dir1.c_str()); 
    cout << "==> File1: "; gDirectory->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    // -- separate mm from rest because the former don't need sbs
    if (0 == fMode) {
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h1 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
      } else {
	if (1 == fMode) h1 = a.sbsDistribution(cut.c_str(), selection);
	if (2 == fMode) h1 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
      }
    }

    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint){
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection));
    }

    fF[file2]->cd(dir2.c_str()); 
    cout << "==> File2: pwd() = "; gDirectory->pwd(); 
    if (0 == fMode) {
      h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
      } else {
	if (1 == fMode) h2 = a.sbsDistribution(cut.c_str(), selection);
	if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
      }
      //      if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection));
    }

    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    c0->cd(); 
    c0->Clear();
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 

    h1->SetMinimum(0.1);
    h1->SetMaximum(dmax); 
    gPad->SetLogy(0); 

    setHist(h1, color1, marker1, 1.5); 
    setFilledHist(h1, color1, color1, fill1); 
    h1->SetTitle("");
    h1->Draw(option1);
    
    setHist(h2, color2, marker2, 1.5); 
    setFilledHist(h2, color2, color2, fill2); 
    h2->Draw(Form("same%s", option2));

    if (string::npos != file1.find("Sg")) {
      legg->AddEntry(h1, "Data (sideband)", loption1); 
    } else {
      legg->AddEntry(h1, "Data", loption1); 
    }
    legg->AddEntry(h2, fName[file2].c_str(), loption2); 
    legg->Draw(); 

    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  }

} 








