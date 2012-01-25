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

  fDoUseBDT = false; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 


  fNumbersFileName = fDirectory + "/anaBmm.plotOverlays." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
plotOverlays::~plotOverlays() {

  cout << "fHistFile = " << fHistFile << endl;
  if (fHistFile) fHistFile->Write();
  if (fHistFile) fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotOverlays::makeAll(int verbose) {

  fVerbose = verbose;

  int all(0);
  
  fMode = 0; 
  cout << " ########################## SG DATA/MC A #########################" << endl;
  sbsDistributionOverlay("SgData", "candAnaMuMu", "A", "SgMc", "candAnaMuMu", "A", "Ao");
  cout << " ########################## SG DATA/MCPU A #########################" << endl;
  sbsDistributionOverlay("SgData", "candAnaMuMu", "A", "SgMcPU", "candAnaMuMu", "A", "Ao");
  cout << " ########################## SG DATA/MCPU APV0 #########################" << endl;
  sbsDistributionOverlay("SgData", "candAnaMuMu", "APV0", "SgMcPU", "candAnaMuMu", "APV0", "Ao"); 
  cout << " ########################## SG DATA/MCPU APV1 #########################" << endl;
  sbsDistributionOverlay("SgData", "candAnaMuMu", "APV1", "SgMcPU", "candAnaMuMu", "APV1", "Ao"); 
  cout << " ########################## MCPU APV0/MCPU APV1 #########################" << endl;
  sbsDistributionOverlay("SgMcPU", "candAnaMuMu", "APV0", "SgMcPU", "candAnaMuMu", "APV1", "Ao"); 

  if (all) sbsDistributionOverlay("SgData", "candAnaMuMu", "B", "SgMc", "candAnaMuMu", "B", "Ao");
  if (all) sbsDistributionOverlay("SgData", "candAnaMuMu", "E", "SgMc", "candAnaMuMu", "E", "Ao");

  // -- For normalization sample, use pol1 + error function for bg parametrization
  fMode = 1; 
  fMode = 3; 
  fPreco = 5.1;
  // -- default/mix
  cout << " ########################## NO DATA/MC A #########################" << endl;
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "A", "NoMc", "candAnaBu2JpsiK", "A", "Ao");
  cout << " ########################## NO DATA/MCPU APV0 #########################" << endl;
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "APV0", "NoMcPU", "candAnaBu2JpsiK", "APV0", "Ao");
  cout << " ########################## NO DATA/MCPU APV1 #########################" << endl;
  sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "APV1", "NoMcPU", "candAnaBu2JpsiK", "APV1", "Ao");
  cout << " ########################## NO MCPU APV0/MCPU APV1 #########################" << endl;
  sbsDistributionOverlay("NoMcPU", "candAnaBu2JpsiK", "APV0", "NoMcPU", "candAnaBu2JpsiK", "APV1", "Ao"); 
  if (all) sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "B", "NoMc", "candAnaBu2JpsiK", "B", "Ao");
  if (all) sbsDistributionOverlay("NoData", "candAnaBu2JpsiK", "E", "NoMc", "candAnaBu2JpsiK", "E", "Ao");

  // -- For control sample, use expo function for bg parametrization
  fMode = 2; 
  // -- default/mix
  cout << " ########################## CS DATA/MC A #########################" << endl;
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "A", "CsMc", "candAnaBs2JpsiPhi", "A", "Ao");
  cout << " ########################## CS DATA/MCPU APV0 #########################" << endl;
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "APV0", "CsMcPU", "candAnaBs2JpsiPhi", "APV0", "Ao");
  cout << " ########################## CS DATA/MCPU APV0 #########################" << endl;
  sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "APV1", "CsMcPU", "candAnaBs2JpsiPhi", "APV1", "Ao");
  cout << " ########################## CS MCPU APV1/MCPU APV0 #########################" << endl;
  sbsDistributionOverlay("CsMcPU", "candAnaBs2JpsiPhi", "APV0", "CsMcPU", "candAnaBs2JpsiPhi", "APV1", "Ao");
  if (all) sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "B", "CsMc", "candAnaBs2JpsiPhi", "B", "Ao");
  if (all) sbsDistributionOverlay("CsData", "candAnaBs2JpsiPhi", "E", "CsMc", "candAnaBs2JpsiPhi", "E", "Ao");
  

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
  
  if ((string::npos != file1.find("Mc")) && (string::npos != file2.find("Mc"))) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    fill1   = 3356; 
    color1  = kBlue; 

    sprintf(option2, "hist"); 
    sprintf(loption2, "f");
    fill2   = 3365; 
    color2  = kRed; 
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
  map<string, double> cutVal; 
  doList.push_back("hlt");  
  doList.push_back("muonsid");  
  doList.push_back("tracksqual"); 
  doList.push_back("pvz");  skipList.push_back("pvz"); 
  doList.push_back("pvn");
  doList.push_back("pvavew8"); leftList.push_back("pvavew8"); 
  doList.push_back("pvntrk");
  doList.push_back("muon1pt");
  doList.push_back("muon2pt");  cutMap.insert(make_pair("muon2pt", "MUPTLO")); cutVal.insert(make_pair("muon2pt", fCuts[0]->m2pt));
  doList.push_back("muonseta");  skipList.push_back("muonseta"); 
  doList.push_back("pt");       cutMap.insert(make_pair("pt", "CANDPTLO"));  cutVal.insert(make_pair("pt", fCuts[0]->pt));
  doList.push_back("eta");            skipList.push_back("eta"); 
  doList.push_back("bdt");            skipList.push_back("bdt"); 

  doList.push_back("fl3d");  
  doList.push_back("fls3d");    cutMap.insert(make_pair("fls3d", "CANDFLS3D")); cutVal.insert(make_pair("fls3d", fCuts[0]->fls3d));
  doList.push_back("flsxy");
  doList.push_back("chi2dof");  cutMap.insert(make_pair("chi2dof", "CANDVTXCHI2")); cutVal.insert(make_pair("chi2dof", fCuts[0]->chi2dof));
  doList.push_back("pchi2dof");       leftList.push_back("pchi2dof"); 
  doList.push_back("alpha");    cutMap.insert(make_pair("alpha", "CANDALPHA")); cutVal.insert(make_pair("alpha", fCuts[0]->alpha));
  doList.push_back("iso");      cutMap.insert(make_pair("iso", "CANDISOLATION"));  leftList.push_back("iso");  cutVal.insert(make_pair("iso", fCuts[0]->iso));
  doList.push_back("docatrk");  cutMap.insert(make_pair("docatrk", "CANDDOCATRK")); cutVal.insert(make_pair("docatrk", fCuts[0]->docatrk));
  doList.push_back("isotrk");
  doList.push_back("closetrk"); cutMap.insert(make_pair("closetrk", "CANDCLOSETRK"));  cutVal.insert(make_pair("closetrk", fCuts[0]->closetrk));
  doList.push_back("lip");      cutMap.insert(make_pair("lip", "CANDLIP")); skipList.push_back("lip"); cutVal.insert(make_pair("lip", fCuts[0]->pvlip));
  doList.push_back("lips");     cutMap.insert(make_pair("lips", "CANDLIPS")); skipList.push_back("lips"); cutVal.insert(make_pair("lips", fCuts[0]->pvlips));
  doList.push_back("ip");       cutMap.insert(make_pair("ip", "CANDIP"));  cutVal.insert(make_pair("ip", fCuts[0]->pvip));
  doList.push_back("ips");      cutMap.insert(make_pair("ips", "CANDIPS")); cutVal.insert(make_pair("ips", fCuts[0]->pvips));
  doList.push_back("maxdoca");  cutMap.insert(make_pair("maxdoca", "MAXDOCA")); cutVal.insert(make_pair("maxdoca", fCuts[0]->maxdoca));

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

//   doList.clear(); 
//   doList.push_back("ips"); 


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
	  break;
	}
      }
    }
    
    //    cutval = cutVal[doList[i]];

    if (doList[i] == "hlt") cutval = 1; 
    if (doList[i] == "muonsid") cutval = 1; 
    if (doList[i] == "tracksqual") cutval = 1; 

    cout << "--> cutval = " << cutval << endl;

    HistCutEfficiency a(h1); 
    a.fVerbose = 1; 
    a.fIncludeOverflow = 0;
    a.eff(h1, cutval); 
    fTEX << formatTex(a.loEff, Form("%s:%s-%s:%s:loEff",  fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.loErr, Form("%s:%s-%s:%s:loEffE", fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.hiEff, Form("%s:%s-%s:%s:hiEff",  fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.hiErr, Form("%s:%s-%s:%s:hiEffE", fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;
    double lo1  = a.loEff; 
    double lo1E = a.loErr; 
    double hi1  = a.hiEff; 
    double hi1E = a.hiErr; 

    a.eff(h2, cutval); 
    fTEX << formatTex(a.loEff, Form("%s:%s-%s:%s:loEff",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.loErr, Form("%s:%s-%s:%s:loEffE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.hiEff, Form("%s:%s-%s:%s:hiEff",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(a.hiErr, Form("%s:%s-%s:%s:hiEffE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;

    double lo2  = a.loEff; 
    double lo2E = a.loErr; 
    double hi2  = a.hiEff; 
    double hi2E = a.hiErr; 

    vector<double> loDelta = computeDelta(lo1, lo1E, lo2, lo2E); 
    cout << "..> delta = " << loDelta[0] << " +/- " << loDelta[1] << endl;
    vector<double> hiDelta = computeDelta(hi1, hi1E, hi2, hi2E); 
    cout << "..> delta = " << hiDelta[0] << " +/- " << hiDelta[1] << endl;
    
    fTEX << formatTex(loDelta[0],  Form("%s:%s-%s:%s:loDelta",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3, 1) << endl;
    fTEX << formatTex(loDelta[1], Form("%s:%s-%s:%s:loDeltaE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
    fTEX << formatTex(hiDelta[0], Form("%s:%s-%s:%s:hiDelta",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3, 1) << endl;
    fTEX << formatTex(hiDelta[1], Form("%s:%s-%s:%s:hiDeltaE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;


    if (doList[i] == "lip" || doList[i] == "lips") {
      HistCutEfficiency a(h1); 
      //      a.fVerbose = 1; 
      cout << "--> cutval = " << cutval << endl;
      a.eff(h1, -cutval, cutval); 
      double lo1  = a.inEff; 
      double lo1E = a.inErr; 
      fTEX << formatTex(a.inEff, Form("%s:%s-%s:%s:inEff",  fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;
      fTEX << formatTex(a.inErr, Form("%s:%s-%s:%s:inEffE", fSuffix.c_str(), file1.c_str(), region1.c_str(), doList[i].c_str()), 3) << endl;

      a.eff(h2, -cutval, cutval); 
      fTEX << formatTex(a.inEff, Form("%s:%s-%s:%s:inEff",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
      fTEX << formatTex(a.inErr, Form("%s:%s-%s:%s:inEffE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
      double lo2  = a.inEff; 
      double lo2E = a.inErr; 
      
      vector<double> loDelta = computeDelta(lo1, lo1E, lo2, lo2E); 
      
      fTEX << formatTex(loDelta[0], Form("%s:%s-%s:%s:inDelta",  fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3, 1) << endl;
      fTEX << formatTex(loDelta[1], Form("%s:%s-%s:%s:inDeltaE", fSuffix.c_str(), file2.c_str(), region2.c_str(), doList[i].c_str()), 3) << endl;
    }      
    
    c0->cd(); 
    c0->Clear();
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 

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
      if (string::npos != file1.find("SgMcPU") && string::npos != file2.find("SgMcPU")) {
	legg->AddEntry(h1, "low PU ", loption1); 
	legg->AddEntry(h2, "high PU ", loption2); 
      } else if (string::npos != file1.find("Sg")) {
	legg->AddEntry(h1, "Data (sideband)", loption1); 
	legg->AddEntry(h2, fName[file2].c_str(), loption2); 
      } else {
	legg->AddEntry(h1, "Data", loption1); 
	legg->AddEntry(h2, fName[file2].c_str(), loption2); 
      }
      legg->Draw(); 
    }

    //    stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  }

} 


// ----------------------------------------------------------------------
vector<double> plotOverlays::computeDelta(double lo1, double lo1E, double lo2, double lo2E) {
  vector<double> result; 

  double a = 2.*(lo1-lo2)/(lo1+lo2);
  double c = lo1+lo2; c = c*c*c*c;
  double b = TMath::Sqrt(16.*lo2*lo2*lo1E*lo1E/c + 16*lo1*lo1*lo2E*lo2E/c);
  cout << "--> delta = " << a << " +/- " << b << endl;
  result.push_back(a); 
  result.push_back(b); 
  return result;
}


// ----------------------------------------------------------------------
void plotOverlays::sbsDistributionOverlaySameFile(string file1, string dir1, string r1, string r2, string sel, string L1, string L2) {

}






