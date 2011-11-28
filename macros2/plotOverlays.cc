#include "plotOverlays.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"

using namespace std; 

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

  int all(0);
  
  fMode = 0; 
  sbsDistributionOverlay("SgData", "candAnaMuMu", "A", "SgMc", "candAnaMuMu", "A", "Ao");
  if (all) sbsDistributionOverlay("SgData", "candAnaMuMu", "B", "SgMc", "candAnaMuMu", "B", "Ao");
  if (all) sbsDistributionOverlay("SgData", "candAnaMuMu", "E", "SgMc", "candAnaMuMu", "E", "Ao");

  // -- For normalization sample, use pol1 + error function for bg parametrization
  fMode = 1; 
  fMode = 3; 
  fPreco = 5.1;
  // -- default/mix
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "A", "NoMc", "candAnaBu2JpsiK", "A", "Ao");
  if (all) sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "B", "NoMc", "candAnaBu2JpsiK", "B", "Ao");
  if (all) sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "E", "NoMc", "candAnaBu2JpsiK", "E", "Ao");
  // -- run ranges
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "AR3", "NoMc2e33", "candAnaBu2JpsiK", "A", "Ao");
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "AR4", "NoMc3e33", "candAnaBu2JpsiK", "A", "Ao");

  // -- For control sample, use expo function for bg parametrization
  fMode = 2; 
  // -- default/mix
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "A", "CsMc", "candAnaBs2JpsiPhi", "A", "Ao");
  if (all) sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "B", "CsMc", "candAnaBs2JpsiPhi", "B", "Ao");
  if (all) sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "E", "CsMc", "candAnaBs2JpsiPhi", "E", "Ao");
  // -- run ranges
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "AR3", "CsMc2e33", "candAnaBs2JpsiPhi", "A", "Ao");
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "AR4", "CsMc3e33", "candAnaBs2JpsiPhi", "A", "Ao");
  

}


// ----------------------------------------------------------------------
void plotOverlays::sbsDistributionOverlay(string file1, string dir1, string region1, string file2, string dir2, string region2, 
					  string selection) {

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
  vector<string> leftList, skipList; 

  map<string, string> cutMap; 
  doList.push_back("pvz");  skipList.push_back("pvz"); 
  doList.push_back("pvn");
  doList.push_back("pvavew8"); leftList.push_back("pvavew8"); 
  doList.push_back("pvntrk");
  doList.push_back("muon1pt");
  doList.push_back("muon2pt"); cutMap.insert(make_pair("muon2pt", "MUPTLO"));
  doList.push_back("muonseta");  skipList.push_back("muonseta"); 
  doList.push_back("pt"); cutMap.insert(make_pair("pt", "CANDPTLO"));
  doList.push_back("eta"); skipList.push_back("eta"); 
  doList.push_back("bdt"); skipList.push_back("bdt"); 

  doList.push_back("fl3d");  
  doList.push_back("fls3d");  cutMap.insert(make_pair("fls3d", "CANDFLS3D"));
  doList.push_back("flsxy");
  doList.push_back("chi2dof"); cutMap.insert(make_pair("chi2dof", "CANDVTXCHI2"));
  doList.push_back("pchi2dof"); leftList.push_back("pchi2dof"); 
  doList.push_back("alpha"); cutMap.insert(make_pair("alpha", "CANDCOSALPHA"));
  doList.push_back("iso");  leftList.push_back("iso"); 
  doList.push_back("docatrk");
  doList.push_back("isotrk");
  doList.push_back("closetrk");
  doList.push_back("lip"); skipList.push_back("lip");
  doList.push_back("lips"); skipList.push_back("lips");
  doList.push_back("lip2"); skipList.push_back("lip2");
  doList.push_back("lips2"); skipList.push_back("lips2");

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




  TCanvas *c1;
  string skipregion, skipcut, cut, pdfname; 
  for (unsigned int i = 0; i < doList.size(); ++i) {
    cut = Form("%s_%s", region1.c_str(), doList[i].c_str()); 
    skipregion =  cut.substr(0, cut.find_first_of("_")); 

    pdfname = Form("%s/%s_%s-%s_sbs_%s-%s_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), region1.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str(), selection.c_str());
    
    cout << pdfname << endl;
    
    fF[file1]->cd(dir1.c_str()); 
    cout << "==> File1: "; gDirectory->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    // -- 0 == fMode is dimuon case: separate this from rest because the former don't need sbs
    if (0 == fMode) {
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h1 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h1 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h1 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h1 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
    }

    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint){
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection.c_str()));
    }

    fF[file2]->cd(dir2.c_str()); 
    cut = Form("%s_%s", region2.c_str(), doList[i].c_str()); 
    cout << "==> File2: pwd() = "; gDirectory->pwd(); 
    if (0 == fMode) {
      h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file2.find("Mc")) {
	// -- For MC use the signal distributions directly
	h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h2 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h2 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
      //      if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection.c_str()));
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

    HistCutEfficiency a(h1); 
    a.eff(h1, cutval); 
    fTEX << formatTex(a.loEff, Form("%s:%s:%s:loEff", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.loErr, Form("%s:%s:%s:loEffE", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.hiEff, Form("%s:%s:%s:hiEff", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.hiErr, Form("%s:%s:%s:hiEffE", fSuffix.c_str(), file1.c_str(), cut.c_str()), 3) << endl;

    a.eff(h2, cutval); 
    fTEX << formatTex(a.loEff, Form("%s:%s:%s:loEff", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.loErr, Form("%s:%s:%s:loEffE", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.hiEff, Form("%s:%s:%s:hiEff", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;
    fTEX << formatTex(a.hiErr, Form("%s:%s:%s:hiEffE", fSuffix.c_str(), file2.c_str(), cut.c_str()), 3) << endl;

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

    //     double kprob = h1->KolmogorovTest(h2); 
    //     double xprob = h1->Chi2Test(h2, "WW"); 
    //     tl->DrawLatex(0.20, 0.92, Form("P(K): %4.3f", kprob));
    //     tl->DrawLatex(0.60, 0.92, Form("P(#chi^{2}): %4.3f", xprob));
    
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






