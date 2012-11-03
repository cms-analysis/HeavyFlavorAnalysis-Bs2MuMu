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
plotReducedOverlays::plotReducedOverlays(const char *files, const char *dir, const char *cuts, int mode) : 
plotClass(files, dir, cuts, mode) { 

  TVirtualFitter::SetMaxIterations(50000);

  cout << "==> plotReducedOverlays files: " << files << " dir: " << dir << " cuts: " << cuts << endl;

  fNumbersFileName = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  fOffset = 0; 
  
  fIsMC = false; 
  fIsSignal = false; 
  fSetup = "A";
  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  fSaveSmallTree = false; 

  fSel0 = false; 
  fSel1 = false; 
  fSel2 = false; 

  fDoList.push_back("muon1pt");
  fDoList.push_back("muon2pt");

  fDoList.push_back("muonseta");
  fDoList.push_back("pt");
  fDoList.push_back("p");
  fDoList.push_back("pz");
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
  fDoList.push_back("bdtsel0");
  fDoList.push_back("bdtsel1");
  fDoList.push_back("bdtsel2");

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
  makeSample(sample2, selection, channel);   
  fOffset = 0; 
  makeSample(sample1, selection, channel);   

  overlay(sample1, sample2, selection); 
  overlay(sample1, sample2, "HLT", "bdt"); 
}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay(string sample1, string sample2, string channel, string selection) {
  
  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  cout << "fHistFile: " << hfname;
  fHistFile = TFile::Open(hfname.c_str());
  cout << " opened " << endl;


  fSetup = channel; 

  sbsDistributions(sample1, selection);
  sbsDistributions(sample1, "HLT", "bdt");

  sbsDistributions(sample2, selection);
  sbsDistributions(sample2, "HLT", "bdt");

  overlay(sample1, sample2, selection); 
  overlay(sample1, sample2, "HLT", "bdt"); 


  // -- save histograms into separate file for systematics sbs_E_NoMc_bdtHLT
  TH1D *h1  = (TH1D*)gDirectory->Get(Form("sbs_%s_%s_bdtHLT", fSetup.c_str(), sample1.c_str())); 
  TH1D *h2  = (TH1D*)gDirectory->Get(Form("sbs_%s_%s_bdtHLT", fSetup.c_str(), sample2.c_str()));
  TDirectory *pD = gDirectory; 

  hfname  = fDirectory + "/anaBmm.plotReducedOverlaysSystematics." + fSuffix + ".root";
  TFile *fl = TFile::Open(hfname.c_str(), "UPDATE");
  h1->SetDirectory(fl);
  h1->Write();
  h2->SetDirectory(fl);
  h2->Write();

  // -- and now all histograms as well
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    h1 = (TH1D*)pD->Get(Form("sbs_%s_%s_%s%s", fSetup.c_str(), sample1.c_str(), fDoList[i].c_str(), selection.c_str())); 
    h1->SetDirectory(fl); 
    h1->Write();
    h2 = (TH1D*)pD->Get(Form("sbs_%s_%s_%s%s", fSetup.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str())); 
    h2->SetDirectory(fl); 
    h2->Write();
  }

  fl->Close();

  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeSample(string mode, string selection, string channel) {

  fSetup = channel;

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
    if (fSetup == "B") {
      SIGBOXMIN = 5.23;
      SIGBOXMAX = 5.33;
    } else {
      SIGBOXMIN = 5.21;
      SIGBOXMAX = 5.34;
    }
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

  if (string::npos != mode.find("Cs")) {
    BGLBOXMIN = 5.10;
    BGLBOXMAX = 5.29;
    if (fSetup == "B") {
      SIGBOXMIN = 5.34;
      SIGBOXMAX = 5.40;
    } else {
      SIGBOXMIN = 5.32;
      SIGBOXMAX = 5.41;
    }
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
  //  loopOverTree(t, mode, 1, 1000000);
  loopOverTree(t, mode, 1);

  fHistFile->Write();
  fHistFile->Close();
}



// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions(string mode) {

  //fHistFile->cd();
  string name = Form("%s_%s_", fSetup.c_str(), mode.c_str());

  cout << "fOffset: " << fOffset << " name = " << name << endl;
  fpMuon1Pt[fOffset]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 60, 0., 30.); 
  fpMuon2Pt[fOffset]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.); 
  fpMuonsEta[fOffset]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", 40, -2.5, 2.5); 
  fpPt[fOffset]        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 60, 0., 60.); 
  fpP[fOffset]         = bookDistribution(Form("%sp", name.c_str()), "p(B) [GeV]", "fGoodPt", 50, 0., 100.); 
  fpPz[fOffset]        = bookDistribution(Form("%spz", name.c_str()), "p_{z}(B) [GeV]", "fGoodPt", 50, 0., 100.); 
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

  fpBDTSel0[fOffset]   = bookSpecialDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodHLT", 200, -1.0, 1.0, &fSel0);   
  fpBDTSel1[fOffset]   = bookSpecialDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodHLT", 200, -1.0, 1.0, &fSel1);   
  fpBDTSel2[fOffset]   = bookSpecialDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodHLT", 200, -1.0, 1.0, &fSel2);   
  
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
    a.fMassPeak = 5.27;
    if (fSetup == "B") {
      a.fMassSigma = 0.02;
    } else {
      a.fMassSigma = 0.03;
    }
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


//   for (unsigned int i = 0; i < fDoList.size(); ++i) {
//     if (restricted) {
//       if (string::npos == fDoList[i].find(what)) continue;
//     }
//     bla =  Form("sbs_%s_%s_%s%s", fSetup.c_str(), mode.c_str(), fDoList[i].c_str(), selection.c_str());
//     cout << bla << "  " ; 
//     h = (TH1D*)gDirectory->Get(bla.c_str());
//     //    h->SetDirectory(fHistFile); 
//     if (h) {
//       cout << "writing to " << fHistFile->GetName() << " " << h << "  " << (h==0? " nada": h->GetName()) << endl;
//       h->Write();
//     } else {
//       cout << "histogram " << bla << " not found" << endl;
//     }
//   }

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
  zone(1);
  c0->cd();

  string sChan = (chan == 0? "B": "E"); 

  double bdtCut = fCuts[chan]->bdt;
  cout << "bdtCut = " << bdtCut << endl;

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlaysSystematics." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "");

  // -- extract here the means of the NPV distributions  
  TH1D *hpv1 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_pvnPresel", sChan.c_str(), sample1.c_str()));
  TH1D *hpv2 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_pvnPresel", sChan.c_str(), sample2.c_str()));

  if (string::npos != sample2.find("CsMc")) {
    setHist(hpv2, kRed, 20, 1.5); 
    setFilledHist(hpv2, kRed, kRed, 3365); 
  } else if (string::npos != sample2.find("SgMc")) {
    setHist(hpv2, kRed, 20, 1.5); 
    setFilledHist(hpv2, kRed, kRed, 3344); 
  } else {
    setHist(hpv2, kBlue, 20, 1.5); 
    setFilledHist(hpv2, kBlue, kBlue, 3365); 
  }

  hpv1->SetMinimum(0.); 
  hpv1->Draw();
  hpv2->Scale(hpv1->GetSumOfWeights()/hpv2->GetSumOfWeights()); 
  hpv2->Draw("samehist"); 
  tl->DrawLatex(0.15, 0.92, Form("Data: %3.2f#pm%3.2f", hpv1->GetMean(), hpv1->GetMeanError())); 
  tl->DrawLatex(0.6, 0.92, Form("MC: %3.2f#pm%3.2f", hpv2->GetMean(), hpv2->GetMeanError())); 
  c0->SaveAs(Form("%s/systematics-npv_%s_%s_chan%d.pdf", fDirectory.c_str(), sample1.c_str(), sample2.c_str(), chan)); 


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

  fTEX << formatTex(fCuts[chan]->bdt, Form("%s:sysCutOnBdtchan%i:val", fSuffix.c_str(), chan), 3) << endl;
  fTEX << formatTex(deltaEps, Form("%s:deltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(adeltaEps, Form("%s:absDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(rdeltaEps, Form("%s:relDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;

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

  tl->DrawLatex(0.2, 0.92, Form("rel difference: %4.3f", rdeltaEps)); 

  double yhi = 0.3*h1->GetMaximum();
  double ylo = 0.;
  pa->DrawArrow(bdtCut, yhi, bdtCut, ylo); 

  c0->SaveAs(Form("%s/systematics_%s_%s_chan%d.pdf", fDirectory.c_str(), sample1.c_str(), sample2.c_str(), chan)); 
  
}


// ----------------------------------------------------------------------
void plotReducedOverlays::massScale(int year) {
  TH1D *h(0); 
  vector<string> dolist, titlist; 
  dolist.push_back("B"); titlist.push_back("Barrel"); 
  dolist.push_back("E"); titlist.push_back("Endcap"); 
  dolist.push_back("APV0"); titlist.push_back("N_{PV} < 6"); 
  dolist.push_back("APV1"); titlist.push_back("N_{PV} > 10"); 

  double hm, hmE, hr, hrE; 
  double fm, fmE, fr, frE; 
  string dset; 

  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  gStyle->SetOptTitle(0); 

  // -- J/psi peak
  TH1D *h1 = new TH1D("h1", "", 100, -0.1, 0.1); 
  double mpsiPdg(3.0969); 
  double mpv0, mpv1; 
  for (int year = 2011; year <= 2012; ++year) {
    dset = (year == 2012?"NoData": "NoData2011"); 
    fF[dset]->cd("candAnaBu2JpsiK"); 
    for (int i = 0; i < dolist.size(); ++i) {
      h   = (TH1D*)gDirectory->Get(Form("%s_mpsiSi0", dolist[i].c_str())); 
      hm  = h->GetMean(); 
      hmE = h->GetMeanError(); 
      hr  = h->GetRMS(); 
      hrE = h->GetRMSError(); 
      h->Fit("gaus"); 
      fm  = h->GetFunction("gaus")->GetParameter(1); 
      fmE = h->GetFunction("gaus")->GetParError(1); 
      fr  = h->GetFunction("gaus")->GetParameter(2); 
      frE = h->GetFunction("gaus")->GetParError(2); 
      
      if ("APV0" == dolist[i]) mpv0 = hm;
      if ("APV1" == dolist[i]) mpv1 = hm;
      h1->Fill((hm-mpsiPdg)/mpsiPdg); 

      h->SetMaximum(1.4*h->GetMaximum()); 
      h->SetAxisRange(2.9, 3.21, "X"); 
      h->SetLineWidth(2); 
      h->Draw("");
      pl->SetLineColor(kRed); 
      pl->DrawLine(mpsiPdg, 0., mpsiPdg, h->GetMaximum()); 
      
      tl->SetTextSize(0.03); 
      tl->DrawLatex(0.2, 0.86, Form("hist mean: %5.4f#pm%5.4f", hm, hmE)); 
      tl->DrawLatex(0.2, 0.80, Form("hist rms:  %5.4f#pm%5.4f", hr, hrE)); 
      
      tl->DrawLatex(0.2, 0.74, Form("fit mean: %5.4f#pm%5.4f", fm, fmE)); 
      tl->DrawLatex(0.2, 0.68, Form("fit #sigma:       %5.4f#pm%5.4f", fr, frE)); 
      
      fTEX << formatTex(hm,  Form("%s:ms-mean-jpsi-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 5) << endl;
      fTEX << formatTex(hmE, Form("%s:ms-mean-jpsi-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 5) << endl;
      fTEX << formatTex(hr,  Form("%s:ms-rms-jpsi-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 5) << endl;
      fTEX << formatTex(hrE, Form("%s:ms-rms-jpsi-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 5) << endl;
      
      tl->SetTextSize(0.05); 
      tl->DrawLatex(0.7, 0.92, titlist[i].c_str()); 
      c0->SaveAs(Form("%s/syst-ms-jpsi-%s-%s.pdf", fDirectory.c_str(), dset.c_str(), dolist[i].c_str())); 
    }
    fTEX << formatTex(mpv1-mpv0, Form("%s:ms-scale-jpsi-%d:diff", fSuffix.c_str(), year), 5) << endl;
  }

  fTEX << formatTex(mpv1-mpv0, Form("%s:ms-scale-jpsi:diff", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(h1->GetRMS(), Form("%s:ms-scale-jpsi:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(h1->GetRMSError(), Form("%s:ms-scale-jpsi:err", fSuffix.c_str()), 5) << endl;
  h1->Reset();


  // -- B+ peak
  double mbuPdg(5.27925); 
  for (int year = 2011; year <= 2012; ++year) {
    dset = (year == 2012?"NoData": "NoData2011"); 
    fF[dset]->cd("candAnaBu2JpsiK"); 
    for (int i = 0; i < dolist.size(); ++i) {
      h   = (TH1D*)gDirectory->Get(Form("%s_mpsiMassAo", dolist[i].c_str())); 
      normYield2(h, 0, 5.0); 
      fm  = h->GetFunction("f1_expo_err_gauss2_landau")->GetParameter(1); 
      fmE = h->GetFunction("f1_expo_err_gauss2_landau")->GetParError(1); 
      fr  = h->GetFunction("f1_expo_err_gauss2_landau")->GetParameter(2); 
      frE = h->GetFunction("f1_expo_err_gauss2_landau")->GetParError(2); 
      
      if ("APV0" == dolist[i]) mpv0 = fm;
      if ("APV1" == dolist[i]) mpv1 = fm;

      h->SetMaximum(1.4*h->GetMaximum()); 
      h->SetAxisRange(2.9, 3.21, "X"); 
      h->SetLineWidth(2); 
      h->Draw("");
      pl->SetLineColor(kRed); 
      pl->DrawLine(mbuPdg, 0., mbuPdg, h->GetMaximum()); 

      h1->Fill((fm-mbuPdg)/mbuPdg); 
      
      tl->SetTextSize(0.04); 
      tl->DrawLatex(0.2, 0.92, Form("fit mean: %5.4f#pm%5.4f", fm, fmE)); 
      
      fTEX << formatTex(fm,  Form("%s:ms-mean-Bu-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(fmE, Form("%s:ms-mean-Bu-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(hr,  Form("%s:ms-rms-Bu-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(hrE, Form("%s:ms-rms-Bu-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      
      tl->SetTextSize(0.05); 
      tl->DrawLatex(0.7, 0.92, titlist[i].c_str()); 
      c0->SaveAs(Form("%s/syst-massscale-Bu-%s-%s.pdf", fDirectory.c_str(), dset.c_str(), dolist[i].c_str())); 
    }
    fTEX << formatTex(mpv1-mpv0, Form("%s:ms-scale-Bu-%d:diff", fSuffix.c_str(), year), 5) << endl;
  }

  fTEX << formatTex(h1->GetRMS(), Form("%s:ms-scale-bu:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(h1->GetRMSError(), Form("%s:ms-scale-bu:err", fSuffix.c_str()), 5) << endl;
  h1->Reset();

  // -- Bs peak
  double mbsPdg(5.36677); 
  for (int year = 2011; year <= 2012; ++year) {
    dset = (year == 2012?"CsData": "CsData2011"); 
    fF[dset]->cd("candAnaBs2JpsiPhi"); 

    for (int i = 0; i < dolist.size(); ++i) {
      h   = (TH1D*)gDirectory->Get(Form("%s_mpsiMassAo", dolist[i].c_str())); 
      if (dolist[i] == "E") {
	csYield2(h, 1, 5.0); 
      } else {
	csYield2(h, 0, 5.0); 
      }	
      fm  = h->GetFunction("f1_expo_err_gauss2f")->GetParameter(1); 
      fmE = h->GetFunction("f1_expo_err_gauss2f")->GetParError(1); 
      fr  = h->GetFunction("f1_expo_err_gauss2f")->GetParameter(2); 
      frE = h->GetFunction("f1_expo_err_gauss2f")->GetParError(2); 

      if ("APV0" == dolist[i]) mpv0 = fm;
      if ("APV1" == dolist[i]) mpv1 = fm;

      h1->Fill((fm-mbsPdg)/mbsPdg); 
      
      h->SetMaximum(1.4*h->GetMaximum()); 
      h->SetAxisRange(2.9, 3.21, "X"); 
      h->SetLineWidth(2); 
      h->Draw("");
      pl->SetLineColor(kRed); 
      pl->DrawLine(mbsPdg, 0., mbsPdg, h->GetMaximum()); 
      
      tl->SetTextSize(0.04); 
      tl->DrawLatex(0.2, 0.92, Form("fit mean: %5.4f#pm%5.4f", fm, fmE)); 

      fTEX << formatTex(fm,  Form("%s:ms-mean-Bs-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(fmE, Form("%s:ms-mean-Bs-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(hr,  Form("%s:ms-rms-Bs-%d-%s:val", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;
      fTEX << formatTex(hrE, Form("%s:ms-rms-Bs-%d-%s:err", fSuffix.c_str(), year, dolist[i].c_str()), 3) << endl;   

      tl->SetTextSize(0.05); 
      tl->DrawLatex(0.7, 0.92, titlist[i].c_str()); 
      c0->SaveAs(Form("%s/play2-Bs-%s-%s.pdf", fDirectory.c_str(), dset.c_str(), dolist[i].c_str())); 
    }
    fTEX << formatTex(mpv1-mpv0, Form("%s:ms-scale-Bs-%d:diff", fSuffix.c_str(), year), 5) << endl;
  }

  fTEX << formatTex(h1->GetRMS(), Form("%s:ms-scale-bs:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(h1->GetRMSError(), Form("%s:ms-scale-bs:err", fSuffix.c_str()), 5) << endl;
  h1->Reset();

}


// ----------------------------------------------------------------------
void plotReducedOverlays::overlay(string sample1, string sample2, string selection, string what) {

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlays." + fSuffix + ".root";
  //  string hfname  = fDirectory + "/test.root";
  //   cout << "fHistFile: " << hfname << endl;
  //   fHistFile = TFile::Open(hfname.c_str());


  
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
    h1 = (TH1D*)gDirectory->Get(n1.c_str());
    cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)gDirectory->Get(n2.c_str());
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


    if (string::npos != fDoList[i].find("npv")) {
      tl->DrawLatex(0.2, 0.92, Form("means: MC(%4.3f) Data(%4.3f)", h1->GetMean(), h2->GetMean()));
    }
    
    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/overlay_%s_%s_%s_%s_%s.pdf", 
		    fDirectory.c_str(), fSetup.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str())); 
  }


}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction(int function, int mode) {
  if (1 == function) loopFunction1(mode);
}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction1(int mode) {

  // -- modify here the fGoodHLT to increase S/B for the BDT distribution
  fGoodHLT        = fb.hlt && fGoodMuonsID && (fBDT > -1.);

  // -- update ana cuts!
  fAnaCuts.update(); 

  if (fSetup == "B" && fChan != 0) return;
  if (fSetup == "E" && fChan != 1) return;

  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 10) {
    fSel0 = true;
  } else {
    fSel0 = false; 
  }

  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 15) {
    fSel1 = true;
  } else {
    fSel1 = false; 
  }

  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 20) {
    fSel2 = true;
  } else {
    fSel2 = false; 
  }

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
  fpPz[fOffset]->fill(a.Pz(), mass);
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
	      
  fpBDTSel0[fOffset]->fill(fBDT, mass);
  fpBDTSel1[fOffset]->fill(fBDT, mass);
  fpBDTSel2[fOffset]->fill(fBDT, mass);
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


// ----------------------------------------------------------------------
AnalysisDistribution* plotReducedOverlays::bookSpecialDistribution(string hn, string ht, string hc, int nbins, double lo, double hi, bool *presel) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi); 
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX); 
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX); 
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX); 
  p->setAnalysisCuts(&fAnaCuts, hc.c_str()); 
  p->setPreselCut(presel); 

  return p; 
}

