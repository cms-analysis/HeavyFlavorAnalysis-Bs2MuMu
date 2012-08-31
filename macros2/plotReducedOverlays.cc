#include "plotReducedOverlays.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std; 

ClassImp(plotReducedOverlays)

// ----------------------------------------------------------------------
plotReducedOverlays::plotReducedOverlays(const char *files, const char *cuts, const char *dir) : 
plotClass(files, cuts, dir, 11) { 

  fIsMC = false; 
  fIsSignal = false; 
  fSetup = "A";
  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 

  fNumbersFileName = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  
  fDoList.push_back("bdt");

  fDoList.push_back("muon1pt");
  fDoList.push_back("muon2pt");

  fDoList.push_back("muonseta");
  fDoList.push_back("pt");
  fDoList.push_back("p");
  fDoList.push_back("eta");
  fDoList.push_back("alpha");

  fDoList.push_back("iso");
  fDoList.push_back("closetrk");
  fDoList.push_back("docatrk");

  fDoList.push_back("chi2dof");
  fDoList.push_back("pchi2dof");
  fDoList.push_back("fls3d");
  fDoList.push_back("fl3d");
  fDoList.push_back("fl3de");

  fDoList.push_back("maxdoca");
  fDoList.push_back("ip");
  fDoList.push_back("ips");
  fDoList.push_back("pvn");
  fDoList.push_back("pvavew8");


}


// ----------------------------------------------------------------------
plotReducedOverlays::~plotReducedOverlays() {
  cout << "closing fHistFile = " << fHistFile << endl;
  if (fHistFile) {
    fHistFile->Close();
  }
}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeAll(string selection) {
  vector<string> dolist; 
  dolist.push_back("B"); 
  dolist.push_back("E"); 

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  for (unsigned i = 0; i < dolist.size(); ++i) {
    fSetup = dolist[i]; 
    makeComparison("SgData", "SgMc", selection); 
    makeComparison("NoData", "NoMc", selection); 
    makeComparison("CsData", "CsMc", selection); 

    overlay("SgMc", "CsData", "HLT", "bdt"); 
    overlay("SgMc", "CsMc", "HLT", "bdt"); 
    overlay("SgMc", "NoMc", "HLT", "bdt"); 
    overlay("NoMc", "CsMc", "HLT", "bdt"); 
    overlay("NoData", "CsData", "HLT", "bdt"); 
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeComparison(string sample1, string sample2, string selection) {

  fOffset = 0; 
  makeSample(sample1, selection);   
  fOffset = 1; 
  makeSample(sample2, selection);   

  overlay(sample1, sample2, selection); 
  overlay(sample1, sample2, "HLT", "bdt"); 
}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay(string sample1, string sample2, string channel, string selection) {
  
  fSetup = channel; 
  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str());
  
  overlay(sample1, sample2, selection); 
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeSample(string mode, string selection) {

  fSample = mode;

  if (string::npos != mode.find("Mc")) {
    fIsMC = true; 
  } else {
    fIsMC = false; 
  }

  MASSMIN = 4.5; 
  MASSMAX = 6.5;
  if (string::npos != mode.find("Sg")) {
    fIsSignal = true; 
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.45;
    BGLBOXMIN = 4.80;
    BGLBOXMAX = 5.20;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 6.00;
  }

  if (string::npos != mode.find("No")) {
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.35;
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.18;
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

  if (string::npos != mode.find("Cs")) {
    SIGBOXMIN = 5.27;
    SIGBOXMAX = 5.47;
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.20;
    BGHBOXMIN = 5.50;
    BGHBOXMAX = 5.70;
  }

  bookDistributions(mode);
  TTree *t = getTree(mode); 

  setupTree(t, mode); 
  //  loopOverTree(t, mode, 1, 100000);
  loopOverTree(t, mode, 1);

  sbsDistributions(mode, selection);
  sbsDistributions(mode, "HLT", "bdt");
}



// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions(string mode) {

  fHistFile->cd();
  string name = Form("%s_%s_", fSetup.c_str(), mode.c_str());

  fpMuon1Pt[fOffset]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 30, 0., 30.); 
  fpMuon2Pt[fOffset]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.); 
  fpMuonsEta[fOffset]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", 20, -2.5, 2.5); 
  fpPt[fOffset]        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 30, 0., 60.); 
  fpP[fOffset]         = bookDistribution(Form("%sp", name.c_str()), "p(B) [GeV]", "fGoodPt", 50, 0., 100.); 
  fpEta[fOffset]       = bookDistribution(Form("%seta", name.c_str()), "#eta(B)", "fGoodEta", 20, -2.5, 2.5); 
  fpAlpha[fOffset]     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", "fGoodAlpha", 40, 0., 0.08); 
  fpIso[fOffset]       = bookDistribution(Form("%siso", name.c_str()),  "isolation", "fGoodIso", 52, 0., 1.04); 
  fpCloseTrk[fOffset]  = bookDistribution(Form("%sclosetrk", name.c_str()),  "N_{trk}^{close}", "fGoodCloseTrack", 10, 0., 10.); 
  fpDocaTrk[fOffset]   = bookDistribution(Form("%sdocatrk", name.c_str()), "d_{ca}^{0} [cm]", "fGoodDocaTrk", 50, 0., 0.25);   

  fpChi2Dof[fOffset]   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", "fGoodChi2", 40, 0., 4.);
  fpPChi2Dof[fOffset]  = bookDistribution(Form("%spchi2dof", name.c_str()),  "P(#chi^{2},dof)", "fGoodChi2", 50, 0., 1.0);    

  fpFLS3d[fOffset]     = bookDistribution(Form("%sfls3d", name.c_str()), "l_{3D}/#sigma(l_{3D})", "fGoodFLS", 50, 0., 150.);  
  fpFL3d[fOffset]      = bookDistribution(Form("%sfl3d", name.c_str()),  "l_{3D} [cm]", "fGoodFLS", 30, 0., 1.5);  
  fpFL3dE[fOffset]     = bookDistribution(Form("%sfl3de", name.c_str()), "#sigma(l_{3D}) [cm]", "fGoodFLS", 25, 0., 0.05);  

  fpMaxDoca[fOffset]   = bookDistribution(Form("%smaxdoca", name.c_str()), "d^{max} [cm]", "fGoodMaxDoca", 60, 0., 0.03);   
  fpIp[fOffset]        = bookDistribution(Form("%sip", name.c_str()), "#delta_{3D} [cm]", "fGoodIp", 50, 0., 0.015);   
  fpIpS[fOffset]       = bookDistribution(Form("%sips", name.c_str()), "#delta_{3D}/#sigma(#delta_{3D})", "fGoodIpS", 50, 0., 4);
  //  fpPvZ[fOffset]       = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 40, -20., 20.);           
  fpPvN[fOffset]       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 40, 0., 40.);           
  fpPvAveW8[fOffset]   = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", 50, 0.5, 1.);           

  fpBDT[fOffset]       = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodHLT", 40, -1.0, 1.0);   
  
}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsDistributions(string mode, string selection, string what) {

  //  string sbsControlPlotsFileName = Form("%s_%s_%d_", mode.c_str(), fSetup.c_str(), fOffset);
  string sbsControlPlotsFileName = "ad";
  fHistFile->cd();

  AnalysisDistribution a(Form("%s_%s_muon1pt", fSetup.c_str(), mode.c_str())); 
  a.fVerbose = 1; 
  //  string bla = Form("%s_%s_%d", fSetup.c_str(), mode.c_str(), fOffset);
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  int type(0); 
  if (string::npos != mode.find("SgData")) type = 0; // sidebands
  if (string::npos != mode.find("Mc"))     type = 1; // signal window
  if (string::npos != mode.find("NoData")) type = 2; // pol1+expo 
  if (string::npos != mode.find("CsData")) type = 3; // expo
  double preco(5.15);
  if (string::npos != mode.find("CsData")) {
    preco = 5.1;
    a.fMassPeak = 5.37;
    a.fMassSigma = 0.06;
  } else {
    a.fMassPeak = -1.;
    a.fMassSigma = -1.;
  }

  TH1D *h(0); 
  bool restricted = (what != ""); 
  string bla; 
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    bla =  Form("%s_%s_%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str());
    if (0 == type) { 
      cout << "=> Looking for sideband histogram " << Form("%s%s1", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s1", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  " 
	   << Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (1 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  " 
	   << Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (2 == type) {
      cout << "=> sbsDistributionPol1ErrGauss histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistributionPol1ErrGauss(bla.c_str(), selection.c_str(), preco);
    } else if (3 == type) {
      cout << "=> sbsDistributionExpoGauss histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistributionExpoGauss(bla.c_str(), selection.c_str());
    }
    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;
    h->Write();
  }


}


// ----------------------------------------------------------------------
void plotReducedOverlays::overlay(string sample1, string sample2, string selection, string what) {
  gStyle->SetOptTitle(0); 
  c0->cd();
  TH1D *h1(0), *h2(0); 
  string n1, n2; 
  bool restricted = (what != ""); 
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    n1 =  Form("sbs_%s_%s_%s%s", fSetup.c_str(), sample1.c_str(), fDoList[i].c_str(), selection.c_str());
    n2 =  Form("sbs_%s_%s_%s%s", fSetup.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str());
    cout << "n1: " << n1 << endl;
    cout << "n2: " << n2 << endl;
    h1 = (TH1D*)fHistFile->Get(n1.c_str());
    h2 = (TH1D*)fHistFile->Get(n2.c_str());
    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    setHist(h1, kBlack, 20, 1.5); 
    if (string::npos != sample2.find("CsMc")) {
      setHist(h2, kRed, 20, 1.5); 
      setFilledHist(h2, kRed, kRed, 3365); 
    } else {
      setHist(h2, kBlue, 20, 1.5); 
      setFilledHist(h2, kBlue, kBlue, 3365); 
    }

    double ymax = (h1->GetMaximum() > h2->GetMaximum()? 1.2*h1->GetMaximum() : 1.2*h2->GetMaximum());
    h1->SetMinimum(0.01);
    h1->SetMaximum(ymax);
    h1->Draw("e");
    h2->Draw("samehist");

    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/overlay_%s_%s_%s_%s_%s.pdf", 
		    fDirectory.c_str(), fSetup.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str())); 
  }


}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction(int mode) {
  if (1 == mode) loopFunction1();
}

// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction1() {

  if (fSetup == "B" && fChan != 0) return;
  if (fSetup == "E" && fChan != 1) return;
  fillDistributions();

}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions() {

  double mass = fb.cm;
  if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;

  TLorentzVector a;
  a.SetPtEtaPhiM(fb.pt,fb.eta,fb.phi,mass);

  fpMuon1Pt[fOffset]->fill(fb.m1pt, mass);
  fpMuon2Pt[fOffset]->fill(fb.m2pt, mass);

  fpMuonsEta[fOffset]->fill(fb.m1eta, mass);
  fpMuonsEta[fOffset]->fill(fb.m2eta, mass);
  fpPt[fOffset]->fill(fb.pt, mass);
  fpP[fOffset]->fill(a.P(), mass);
  fpEta[fOffset]->fill(fb.eta, mass);
  fpAlpha[fOffset]->fill(fb.alpha, mass);

  fpIso[fOffset]->fill(fb.iso, mass);
  fpCloseTrk[fOffset]->fill(fb.closetrk, mass);
  fpDocaTrk[fOffset]->fill(fb.docatrk, mass);

  fpChi2Dof[fOffset]->fill(fb.chi2/fb.dof, mass); 
  fpPChi2Dof[fOffset]->fill(fb.pchi2dof, mass); 
	      
  fpFLS3d[fOffset]->fill(fb.fls3d, mass);  
  fpFL3d[fOffset]->fill(fb.fl3d, mass); 
  fpFL3dE[fOffset]->fill(fb.fl3dE, mass); 
	      
  fpMaxDoca[fOffset]->fill(fb.maxdoca, mass); 
  fpIp[fOffset]->fill(fb.pvip, mass); 
  fpIpS[fOffset]->fill(fb.pvips, mass); 
  //  fpPvZ[fOffset]->fill(fb.pvz, mass);
  fpPvN[fOffset]->fill(fb.pvn, mass); 
  fpPvAveW8[fOffset]->fill(fb.pvw8, mass);

  fpBDT[fOffset]->fill(fBDT, mass);
	      

}




// ----------------------------------------------------------------------
AnalysisDistribution* plotReducedOverlays::bookDistribution(string hn, string ht, string hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi); 
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX); 
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX); 
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX); 
  p->setAnalysisCuts(&fAnaCuts, hc.c_str()); 
  p->setPreselCut(&fPreselection); 

  return p; 
}


