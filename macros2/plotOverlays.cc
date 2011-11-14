#include "plotOverlays.hh"

#include "../macros/AnalysisDistribution.hh"
#include "TMath.h"

using namespace std; 
using std::string; 

ClassImp(plotOverlays)

// ----------------------------------------------------------------------
plotOverlays::plotOverlays(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fMode = 2; // expo+gauss
  fMode = 1; // pol1+gauss
}


// ----------------------------------------------------------------------
plotOverlays::~plotOverlays() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotOverlays::makeAll(int verbose) {

  fVerbose = verbose;
  
  fMode = 0; 
  sbsDistributionOverlay("SgData", "candAnaMuMu", "SgMc", "candAnaMuMu", "Ao", "A");
  sbsDistributionOverlay("SgData", "candAnaMuMu", "SgMc", "candAnaMuMu", "Ao", "B");
  sbsDistributionOverlay("SgData", "candAnaMuMu", "SgMc", "candAnaMuMu", "Ao", "E");

  // -- For normalization sample, use pol1 + error function for bg parametrization
  fMode = 1; 
  fMode = 3; 
  fPreco = 5.1;
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "NoMc", "candAnaBu2JpsiK", "Ao", "A");
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "NoMc", "candAnaBu2JpsiK", "Ao", "B");
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "NoMc", "candAnaBu2JpsiK", "Ao", "E");

  // -- For control sample, use expo function for bg parametrization
  fMode = 2; 
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "CsMc", "candAnaBs2JpsiPhi", "Ao", "A");
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "CsMc", "candAnaBs2JpsiPhi", "Ao", "B");
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "CsMc", "candAnaBs2JpsiPhi", "Ao", "E");
  

}


// ----------------------------------------------------------------------
void plotOverlays::sbsDistributionOverlay(string file1, string dir1, string file2, string dir2, const char *selection, const char *region) {

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
  a.fVerbose = fVerbose; 
  TH1D *hcuts = (TH1D*)gDirectory->Get("hcuts"); 

  vector<string> doList; 
  map<string, string> cutMap; 
  doList.push_back("pvz"); 
  doList.push_back("pvn");
  doList.push_back("pvntrk");
  doList.push_back("muon1pt");
  doList.push_back("muon2pt"); cutMap.insert(make_pair("muon2pt", "MUPTLO"));
  doList.push_back("muonseta");
  doList.push_back("pt"); cutMap.insert(make_pair("pt", "CANDPTLO"));
  doList.push_back("eta");

  doList.push_back("fls3d");  cutMap.insert(make_pair("fls3d", "CANDFLS3D"));
  doList.push_back("flsxy");
  doList.push_back("chi2dof"); cutMap.insert(make_pair("chi2dof", "CANDVTXCHI2"));
  doList.push_back("pchi2dof");
  doList.push_back("alpha"); cutMap.insert(make_pair("alpha", "CANDCOSALPHA"));
  doList.push_back("iso");
  doList.push_back("docatrk");
  doList.push_back("isotrk");

  if (string::npos != file1.find("No")) {
    doList.push_back("kaonpt");
    doList.push_back("psipt");
  } 

  if (string::npos != file1.find("Cs")) {
    doList.push_back("kaonspt");
    doList.push_back("psipt");
    doList.push_back("phipt");
    doList.push_back("deltar");
    doList.push_back("mkk");
  }

  vector<string> leftList; 
  leftList.push_back("iso"); 

  vector<string> skipList; 
  skipList.push_back("pvz"); 
  skipList.push_back("psipt"); 
  skipList.push_back("muonseta"); 
  skipList.push_back("pchi2dof"); 


  TCanvas *c1;
  string skipregion, skipcut, cut, pdfname; 
  for (unsigned int i = 0; i < doList.size(); ++i) {
    cut = Form("%s_%s", region, doList[i].c_str()); 
    skipregion =  cut.substr(0, cut.find_first_of("_")); 

    pdfname = Form("%s/%s_%s-%s_sbs_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), cut.c_str(), selection);
    
    cout << pdfname << endl;
    
    fF[file1]->cd(dir1.c_str()); 
    cout << "==> File1: "; gDirectory->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    // -- 0 == fMode is dimuon case: separate this from rest because the former don't need sbs
    if (0 == fMode) {
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h1 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
      } else {
	if (1 == fMode) h1 = a.sbsDistribution(cut.c_str(), selection);
	if (2 == fMode) h1 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
	if (3 == fMode) h1 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection, fPreco);
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
      if (string::npos != file2.find("Mc")) {
	// -- For MC use the signal distributions directly
	h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
      } else {
	if (1 == fMode) h2 = a.sbsDistribution(cut.c_str(), selection);
	if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
	if (3 == fMode) h2 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection, fPreco);
      }
      //      if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection));
    }

    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    // -- Dump hist-based efficiencies
    double cutval = 0.75; 
    string label;
    string bla = cutMap[doList[i]];
    if (bla != "") {
      for (int j = 1; j < hcuts->GetNbinsX(); ++j) {
	label = string(hcuts->GetXaxis()->GetBinLabel(j));
	if (string::npos != label.find(bla)) {
	  bla = label.substr(label.find_last_of("::")+1);
	  cutval = atof(bla.c_str());
	  if (doList[i] == "alpha") cutval = 0.05; 
	  break;
	}
      }
    }
  
    double ntot   = h1->Integral(0, h1->GetNbinsX()+1);
    //    double ntot   = h1->Integral();
    double ucut   = h1->Integral(h1->FindBin(cutval), h1->GetNbinsX()+1);
    double ueff   = ucut/ntot; 
    double ustat  = dEff(static_cast<int>(ucut), static_cast<int>(ntot));

    double lcut   = h1->Integral(0, h1->FindBin(cutval)-1);
    double leff   = lcut/ntot; 
    double lstat  = dEff(static_cast<int>(lcut), static_cast<int>(ntot));

    fTEX << formatTex(ueff, Form("%s:%s:%s:ueff", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(ustat, Form("%s:%s:%s:ueffE", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(leff, Form("%s:%s:%s:leff", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(lstat, Form("%s:%s:%s:leffE", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;


    ntot   = h2->GetSumOfWeights();
    ucut   = h2->Integral(h2->FindBin(cutval), h2->GetNbinsX());
    ueff   = ucut/ntot; 
    ustat  = dEff(static_cast<int>(ucut), static_cast<int>(ntot));
    
    lcut   = h2->Integral(1, h2->FindBin(cutval)-1);
    leff   = lcut/ntot; 
    lstat  = dEff(static_cast<int>(lcut), static_cast<int>(ntot));

    fTEX << formatTex(ueff, Form("%s:%s:%s:ueff", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(ustat, Form("%s:%s:%s:ueffE", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(leff, Form("%s:%s:%s:leff", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(lstat, Form("%s:%s:%s:leffE", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;

    c0->cd(); 
    c0->Clear();
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 

    //     if (doList[i] == "docatrk") {
    //       h1->GetXaxis()->SetTitle("d_{ca}^{min} [cm]");
    //     }

    if (string::npos != file1.find("Sg") && doList[i] == "fls3d") {
      h1->SetMaximum(2.*dmax); 
      gPad->SetLogy(1); 
    } else {
      h1->SetMaximum(dmax); 
      gPad->SetLogy(0); 
    }
    h1->SetMinimum(0.1);
    string xtitle = h1->GetXaxis()->GetTitle();
    if (string::npos != xtitle.find("[")) {
      string unit = xtitle.substr(xtitle.find("[")+1, xtitle.find("]")-xtitle.find("[")-1); 
      cout << "%%%%%%%%% > unit: " << unit << endl;
      if (TMath::Abs(h1->GetBinWidth(1) - 1.) < 0.1) {
	setTitles(h1, h1->GetXaxis()->GetTitle(), Form("Candidates / %2.0f %s", h1->GetBinWidth(1), unit.c_str()), 0.05, 1.1, 1.6); 
      } else if (h1->GetBinWidth(1) < 0.1) {
	setTitles(h1, h1->GetXaxis()->GetTitle(), Form("Candidates / %4.3f %s", h1->GetBinWidth(1), unit.c_str()), 0.05, 1.1, 1.6); 
      }	else {
	setTitles(h1, h1->GetXaxis()->GetTitle(), Form("Candidates / %2.1f %s", h1->GetBinWidth(1), unit.c_str()), 0.05, 1.1, 1.6); 
      }
    } else {
      setTitles(h1, h1->GetXaxis()->GetTitle(), Form("Candidates"), 0.05, 1.1, 1.6); 
    }
    
    setHist(h1, color1, marker1, 1.5); 
    setFilledHist(h1, color1, color1, fill1); 
    h1->SetTitle("");
    h1->Draw(option1);
    
    setHist(h2, color2, marker2, 1.5); 
    setFilledHist(h2, color2, color2, fill2); 
    h2->Draw(Form("same%s", option2));

    if (skipList.end() != find(skipList.begin(), skipList.end(), doList[i])) {
      // do nothing 
    } else {
      if (leftList.end() != find(leftList.begin(), leftList.end(), doList[i])) {
	newLegend(0.25, 0.7, 0.50, 0.85); 
      } else {
	newLegend(0.50, 0.7, 0.75, 0.85); 
      }
      //      legg->AddEntry(h1, fName[file1].c_str(), loption1); 
      if (string::npos != file1.find("Sg")) {
	legg->AddEntry(h1, "Data (sideband)", loption1); 
      } else {
	legg->AddEntry(h1, "Data", loption1); 
      }
      legg->AddEntry(h2, fName[file2].c_str(), loption2); 
      legg->Draw(); 
    }

    //    stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  }

} 


// ----------------------------------------------------------------------
void plotOverlays::sbsDistributionOverlaySameFile(string file1, string dir1, string r1, string r2, string sel, string L1, string L2) {

}






