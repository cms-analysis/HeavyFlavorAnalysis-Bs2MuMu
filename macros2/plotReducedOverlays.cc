#include "plotReducedOverlays.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVirtualFitter.h"

using namespace std; 

ClassImp(plotReducedOverlays)

// ----------------------------------------------------------------------
plotReducedOverlays::plotReducedOverlays(const char *files, const char *dir, const char *cuts) : 
plotClass(files, dir, cuts, 11) { 

  TVirtualFitter::SetMaxIterations(50000);

  fNumbersFileName = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  fIsMC = false; 
  fIsSignal = false; 
  fSetup = "A";
  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 

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

  fDoList.push_back("bdt");

  //  gROOT->cd();

}


// ----------------------------------------------------------------------
plotReducedOverlays::~plotReducedOverlays() {

}



// ----------------------------------------------------------------------
void plotReducedOverlays::makeAll(string selection) {
//   vector<string> dolist; 
//   dolist.push_back("B"); 
//   dolist.push_back("E"); 

//   string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
//   cout << "fHistFile: " << hfname << endl;
//   fHistFile = TFile::Open(hfname.c_str(), "UPDATE");

//   for (unsigned i = 0; i < dolist.size(); ++i) {
//     fSetup = dolist[i]; 
//     makeSampleOverlay("NoData", "NoMc", selection); 
//     makeSampleOverlay("CsData", "CsMc", selection); 
//     makeSampleOverlay("SgData", "SgMc", selection); 

//     // -- unconventional overlays
//     overlay("SgMc", "CsData", "HLT", "bdt"); 
//     overlay("SgMc", "CsMc", "HLT", "bdt"); 
//     overlay("SgMc", "NoMc", "HLT", "bdt"); 
//     overlay("NoMc", "CsMc", "HLT", "bdt"); 
//     overlay("NoData", "CsData", "HLT", "bdt"); 
//   }
//   //  systematics("SgMc", "NoMc"); 
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeSampleOverlay(string sample1, string sample2, string channel, string selection) {

  fSetup = channel;

  fOffset = 1; 
  makeSample(sample2, selection);   
  fOffset = 0; 
  makeSample(sample1, selection);   

  overlay(sample1, sample2, selection); 
  overlay(sample1, sample2, "HLT", "bdt"); 
}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay(string sample1, string sample2, string channel, string selection) {
  
  fSetup = channel; 

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "UPDATE");

  // -- debug problems
  //   sbsDistributions(sample1, selection);
  //   sbsDistributions(sample1, "HLT", "bdt");
  //   sbsDistributions(sample2, selection);
  //   sbsDistributions(sample2, "HLT", "bdt");
  // -- debug problems

  overlay(sample1, sample2, selection); 
  overlay(sample1, sample2, "HLT", "bdt"); 

  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeSample(string mode, string selection) {


  // -- dump histograms
  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname;
  fHistFile = TFile::Open(hfname.c_str(), "UPDATE");
  cout << " opened " << endl;

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
    BGLBOXMIN = 4.80;
    BGLBOXMAX = 5.20;
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.45;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 6.00;
  }

  if (string::npos != mode.find("No")) {
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.18;
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.35;
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

  if (string::npos != mode.find("Cs")) {
    BGLBOXMIN = 5.10;
    BGLBOXMAX = 5.29;
    SIGBOXMIN = 5.30;
    SIGBOXMAX = 5.43;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 5.70;
  }

  bookDistributions(mode);
  TTree *t = getTree(mode); 
  cout << "getTree(" << mode << "): " << t << endl;
  if (0 == t) {
    cout << "tree for mode = " << mode << " not found" << endl;
    return;
  }
  setupTree(t, mode); 
  //loopOverTree(t, mode, 1, 200000);
  loopOverTree(t, mode, 1);

  sbsDistributions(mode, selection);
  sbsDistributions(mode, "HLT", "bdt");

  fHistFile->Write();
  fHistFile->Close();
  //  gROOT->cd();
}



// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions(string mode) {

  //fHistFile->cd();
  string name = Form("%s_%s_", fSetup.c_str(), mode.c_str());

  cout << "bla " << endl;
  fpMuon1Pt[fOffset]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 60, 0., 30.); 
  fpMuon2Pt[fOffset]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.); 
  fpMuonsEta[fOffset]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", 40, -2.5, 2.5); 
  fpPt[fOffset]        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 60, 0., 60.); 
  fpP[fOffset]         = bookDistribution(Form("%sp", name.c_str()), "p(B) [GeV]", "fGoodPt", 50, 0., 100.); 
  fpEta[fOffset]       = bookDistribution(Form("%seta", name.c_str()), "#eta(B)", "fGoodEta", 40, -2.5, 2.5); 
  fpAlpha[fOffset]     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", "fGoodAlpha", 50, 0., 0.1); 
  fpIso[fOffset]       = bookDistribution(Form("%siso", name.c_str()),  "isolation", "fGoodIso", 52, 0., 1.04); 
  fpCloseTrk[fOffset]  = bookDistribution(Form("%sclosetrk", name.c_str()),  "N_{trk}^{close}", "fGoodCloseTrack", 10, 0., 10.); 
  fpDocaTrk[fOffset]   = bookDistribution(Form("%sdocatrk", name.c_str()), "d_{ca}^{0} [cm]", "fGoodDocaTrk", 50, 0., 0.20);   

  fpChi2Dof[fOffset]   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", "fGoodChi2", 40, 0., 4.);
  fpPChi2Dof[fOffset]  = bookDistribution(Form("%spchi2dof", name.c_str()),  "P(#chi^{2},dof)", "fGoodChi2", 50, 0., 1.0);    

  fpFLS3d[fOffset]     = bookDistribution(Form("%sfls3d", name.c_str()), "l_{3D}/#sigma(l_{3D})", "fGoodFLS", 60, 0., 120.);  
  fpFL3d[fOffset]      = bookDistribution(Form("%sfl3d", name.c_str()),  "l_{3D} [cm]", "fGoodFLS", 60, 0., 1.5);  
  fpFL3dE[fOffset]     = bookDistribution(Form("%sfl3de", name.c_str()), "#sigma(l_{3D}) [cm]", "fGoodFLS", 50, 0., 0.05);  

  fpMaxDoca[fOffset]   = bookDistribution(Form("%smaxdoca", name.c_str()), "d^{max} [cm]", "fGoodMaxDoca", 60, 0., 0.03);   
  fpIp[fOffset]        = bookDistribution(Form("%sip", name.c_str()), "#delta_{3D} [cm]", "fGoodIp", 50, 0., 0.015);   
  fpIpS[fOffset]       = bookDistribution(Form("%sips", name.c_str()), "#delta_{3D}/#sigma(#delta_{3D})", "fGoodIpS", 50, 0., 4);
  //  fpPvZ[fOffset]       = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 40, -20., 20.);           
  fpPvN[fOffset]       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 40, 0., 40.);           
  fpPvAveW8[fOffset]   = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", 50, 0.5, 1.);           

  fpBDT[fOffset]       = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodHLT", 200, -1.0, 1.0);   
  
}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsDistributions(string mode, string selection, string what) {

  string sbsControlPlotsFileName = "ad";
  //  fHistFile->cd();

  AnalysisDistribution a(Form("%s_%s_muon1pt", fSetup.c_str(), mode.c_str())); 
  a.fVerbose = 1; 
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  int type(0); 
  if (string::npos != mode.find("SgData")) type = 0; // sidebands
  if (string::npos != mode.find("Mc"))     type = 1; // signal window
  if (string::npos != mode.find("NoData")) type = 2; // pol1+err
  if (string::npos != mode.find("CsData")) type = 3; // expo
  double preco(5.15);
  if (string::npos != mode.find("CsData")) {
    preco = 5.2;
    a.fMassPeak = 5.37;
    a.fMassSigma = 0.06;
  } else {
    a.fMassPeak = -1.;
    a.fMassSigma = -1.;
  }

  cout << "gDIRECTORY: "; gDirectory->pwd();
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
    //    h->Write();
  }


  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    bla =  Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str());
    cout << bla << "  " ; 
    h = (TH1D*)gDirectory->Get(bla.c_str());
    //    h->SetDirectory(fHistFile); 
    if (h) {
      cout << "writing to " << fHistFile->GetName() << " " << h << "  " << (h==0? " nada": h->GetName()) << endl;
      h->Write();
    } else {
      cout << "histogram " << bla << " not found" << endl;
    }
  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::allSystematics() {
  for (int i = 0; i < 2; ++i) {
    systematics("CsData", "CsMc", i); 
    systematics("NoData", "NoMc", i); 
    systematics("CsData", "SgMc", i); 
    systematics("SgMc",   "CsMc", i); 
  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::systematics(string sample1, string sample2, int chan) {
  gStyle->SetOptTitle(0); 
  c0->cd();

  string sChan = (chan == 0? "B": "E"); 

  double bdtCut = fCuts[chan]->bdt;
  cout << "bdtCut = " << bdtCut << endl;

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "");

  TH1D *h1 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_bdtHLT", sChan.c_str(), sample1.c_str()));
  for (int i = 1; i <= h1->GetNbinsX(); ++i) if (h1->GetBinContent(i) < 0) h1->SetBinContent(i, -0.0001); 
  HistCutEfficiency a1(h1, bdtCut, 0); 
  double eps1 = a1.hiEff;
  double eps1E = a1.hiErr;
  cout << "eps1 = " << eps1 << " lo eff = " << a1.loEff << endl;

  TH1D *h2 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_bdtHLT", sChan.c_str(), sample2.c_str()));
  for (int i = 1; i <= h2->GetNbinsX(); ++i) if (h2->GetBinContent(i) < 0) h2->SetBinContent(i, -0.0001); 
  HistCutEfficiency a2(h2, bdtCut, 0); 
  double eps2 = a2.hiEff;
  double eps2E = a2.hiErr;

  cout << "eps2 = " << eps2 << " lo eff = " << a2.loEff << endl;
  double deltaEps = eps1-eps2; 
  double adeltaEps = TMath::Abs(eps1-eps2); 
  double rdeltaEps = 2.*TMath::Abs(eps1-eps2)/(eps1+eps2); 

  fTEX << formatTex(deltaEps, Form("%s:deltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(adeltaEps, Form("%s:absDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(rdeltaEps, Form("%s:relDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;

  zone(1);

  if (string::npos != sample2.find("CsMc")) {
    setHist(h2, kRed, 20, 1.5); 
    setFilledHist(h2, kRed, kRed, 3365); 
  } else if (string::npos != sample2.find("SgMc")) {
    setHist(h2, kRed, 20, 1.5); 
    setFilledHist(h2, kRed, kRed, 3344); 
  } else {
    setHist(h2, kBlue, 20, 1.5); 
    setFilledHist(h2, kBlue, kBlue, 3365); 
  }

  h1->SetMinimum(0.);
  h1->SetMaximum(1.3*h1->GetMaximum());
  h1->Draw("e");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("histsame");

  newLegend(0.2, 0.77, 0.45, 0.87); 
  legg->SetTextSize(0.035);  
  string text1, text2;
  if (string::npos != sample1.find("CsData")) text1 = "B_{s}^{0} #rightarrow J/#psi #phi (data)";
  if (string::npos != sample1.find("NoData")) text1 = "B^{+} #rightarrow J/#psi K (data)";
  if (string::npos != sample1.find("SgMc"))   text1 = "B_{s}^{0} #rightarrow #mu #mu (MC)";

  if (string::npos != sample2.find("SgMc")) text2 = "B_{s}^{0} #rightarrow #mu #mu (MC)";
  if (string::npos != sample2.find("CsMc")) text2 = "B_{s}^{0} #rightarrow J/#psi #phi (MC)";
  if (string::npos != sample2.find("NoMc")) text2 = "B^{+} #rightarrow J/#psi K (MC)";


  legg->AddEntry(h1, Form("#varepsilon = %4.3f#pm%4.3f, %s", eps1, eps1E, text1.c_str()), "p");
  legg->AddEntry(h2, Form("#varepsilon = %4.3f#pm%4.3f, %s", eps2, eps2E, text2.c_str()), "f");
  legg->Draw();

  double yhi = 0.3*h1->GetMaximum();
  double ylo = 0.;
  pa->DrawArrow(bdtCut, yhi, bdtCut, ylo); 

  c0->SaveAs(Form("%s/systematics_%s_%s_chan%d.pdf", fDirectory.c_str(), sample1.c_str(), sample2.c_str(), chan)); 
  

// -- This is the old cumulative distribution approach from Gemma
//   TH1D *i1 = new TH1D("i1", "", h1->GetNbinsX(), h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()+1)); setHist(i1, kBlack); 
//   for (int i = 1; i <= h1->GetNbinsX(); ++i) i1->SetBinContent(i, h1->Integral(1, i)); 
//   i1->Scale(1./i1->GetBinContent(h1->GetNbinsX()));
//   int bin1 = i1->FindBin(bdtCut)-1; 
//   double eff1 = i1->GetBinContent(bin1);
//   double frac1 = i1->GetBinContent(bin1); 
//   cout << "eff1 = " << eff1 << endl;

//   TH1D *i2 = new TH1D("i2", "", h2->GetNbinsX(), h2->GetBinLowEdge(1), h2->GetBinLowEdge(h2->GetNbinsX()+1)); setHist(i2, kRed); 
//   for (int i = 1; i <= h2->GetNbinsX(); ++i) i2->SetBinContent(i, h2->Integral(1, i)); 
//   i2->Scale(1./i2->GetBinContent(h1->GetNbinsX()));
//   int bin2 = i2->FindFirstBinAbove(frac1)-1; 
//   double bdt2 = i2->GetBinCenter(bin2); 
//   double eff2 = i2->GetBinContent(bin2);
//   cout << "bin1 = " << bin1 << " frac1 = " << frac1 << " -> bdt2 = " << bdt2 << " frac = " << eff2 << " at bin " << bin2 << endl;
//   cout << "eff2 = " << eff2 << endl;

//   cout << "sys1 = " << 2.*((eps1-eps2)/(eps1+eps2)) << endl;
//   cout << "sys2 = " << 2.*((1-eff1) - (1-eff2))/((1-eff1) + (1-eff2)) << endl;
//   i1->Draw("hist");
//   i2->Draw("histsame");

//   c0->SaveAs(Form("%s/cumulative_%s_%s_%s_chan%d.pdf", 
// 		  fDirectory.c_str(), fSetup.c_str(), sample1.c_str(), sample2.c_str(), chan)); 

}

// ----------------------------------------------------------------------
void plotReducedOverlays::overlay(string sample1, string sample2, string selection, string what) {

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  //  string hfname  = fDirectory + "/test.root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str());

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
    h1 = (TH1D*)fHistFile->Get(n1.c_str());
    cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)fHistFile->Get(n2.c_str());
    cout << "n2: " << n2 << " -> " << h2 << endl;
    if (0 == h1 || 0 == h2) {
      cout << "  histograms not found" << endl;
      continue;
    }
    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
    
    h1->Draw();
    h2->Draw("samehist");
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

  // -- modify here the fGoodHLT to increase S/B for the BDT distribution
  fGoodHLT        = fb.hlt && fGoodMuonsID && (fb.fls3d > 5) && (fb.chi2/fb.dof < 5) && (fBDT > -1.);
  fGoodHLT        = fb.hlt && fGoodMuonsID && (fBDT > -1.);

  // -- modify even more for the control sample
  if (string::npos != fSample.find("Cs")) {
    fGoodHLT      = fb.hlt && fGoodMuonsID && fGoodMuonsPt 
      && (fb.fls3d > 5) && (fb.alpha < 0.1)&& (fb.chi2/fb.dof < 4) 
      && (fb.iso > 0.5) && (fb.closetrk < 5)
      && (fBDT > -1.)
      ;
    fGoodHLT      = fb.hlt && fGoodMuonsID && (fBDT > -1.);
  }

  // -- update ana cuts!
  fAnaCuts.update(); 

  if (fSetup == "B" && fChan != 0) return;
  if (fSetup == "E" && fChan != 1) return;
  fillDistributions();

}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions() {

  //  cout << "BDT: " << fBDT << " mass =  " << fb.m << " pt = " << fb.pt << endl;

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



