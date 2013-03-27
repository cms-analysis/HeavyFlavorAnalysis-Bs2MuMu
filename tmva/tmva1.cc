#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "tmvaglob.C"
#include "TMVA/Config.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TText.h"
#include "TH2.h"
#include "TGraph.h"
#include "TPaveStats.h"

#include "tmva1.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)#
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

//2011
//#define LUMISCALE 2.01e-4

//2012  12/7405 = 0.00162
#define LUMISCALE 0.00162

ClassImp(tmva1)

using namespace std; 

// ----------------------------------------------------------------------
// -- 
// -- USAGE: a.makeAll(0, 1); > TMVA-0.log
// --
// ----------------------------------------------------------------------

tmva1::tmva1(int year) {

  cout << "tmva1 hello: setup for year = " << year << endl;

  legg = 0;
  legge = 0; 
  tl = new TLatex(); 
  tl->SetTextFont(42);
  tl->SetTextSize(0.03);
  tl->SetNDC(kTRUE); 

  fYear = year; 

  if (year == 2011) {
      fInputFiles.sname = "/scratch/ursl/bdt/v14-2011-mix-Bs2MuMu.root"; 
      fInputFiles.dname = "/scratch/ursl/bdt/v14-2011-bmmLoose.root";
  } else {
      fInputFiles.sname = "/scratch/ursl/bdt/v14-2012-cms-BsToMuMu.root"; 
      fInputFiles.dname = "/scratch/ursl/bdt/v14-2012-bmmLoose.root";
  }

  // -- BDT setup 108/109
  fBdtSetup.NTrees = 800;
  fBdtSetup.nEventsMin = 50; 
  fBdtSetup.MaxDepth = 2;
  fBdtSetup.nCuts = 20; 
  fBdtSetup.AdaBoostBeta = 1.0; 
  fBdtSetup.NNodesMax = 5;

  fApplyOn0 = false;
  fApplyOn1 = false;
  fApplyOn2 = false;
  fTrainAntiMuon = false; 
  fChannel = 0; 

}


// ----------------------------------------------------------------------
tmva1::~tmva1() {
  cout << "tmva1 good bye " << endl;
}


// ----------------------------------------------------------------------
TCanvas* tmva1::getC0() {
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
  if (0 == c0) c0 = new TCanvas("c0","--c0--",2303,0,656,700);
  return c0; 
} 



// ----------------------------------------------------------------------
void tmva1::setupTree(TTree *t, RedTreeData &b) {
  t->SetBranchAddress("evt", &b.evt);
  t->SetBranchAddress("gmuid", &b.gmuid);
  t->SetBranchAddress("hlt", &b.hlt);
  t->SetBranchAddress("m1pt", &b.m1pt);
  t->SetBranchAddress("m2pt", &b.m2pt);
  t->SetBranchAddress("m1eta", &b.m1eta);
  t->SetBranchAddress("m2eta", &b.m2eta);
  t->SetBranchAddress("pt", &b.pt);
  t->SetBranchAddress("eta", &b.eta);
  t->SetBranchAddress("pvlip", &b.pvlip);
  t->SetBranchAddress("pvlips", &b.pvlips);
  t->SetBranchAddress("fl3d", &b.fl3d);
  t->SetBranchAddress("fls3d", &b.fls3d);
  t->SetBranchAddress("flsxy", &b.flsxy);
  t->SetBranchAddress("alpha", &b.alpha);
  t->SetBranchAddress("maxdoca", &b.maxdoca);
  t->SetBranchAddress("pvip", &b.pvip);
  t->SetBranchAddress("pvips", &b.pvips);
  t->SetBranchAddress("iso", &b.iso);
  t->SetBranchAddress("docatrk", &b.docatrk);
  t->SetBranchAddress("chi2", &b.chi2);
  t->SetBranchAddress("dof", &b.dof);
  t->SetBranchAddress("closetrk", &b.closetrk);
  t->SetBranchAddress("m", &b.m);

  t->SetBranchAddress("m1iso",&b.m1iso);
  t->SetBranchAddress("m2iso",&b.m2iso);
  t->SetBranchAddress("pvdchi2",&b.pvdchi2);
  t->SetBranchAddress("m1xpdist",&b.m1xpdist);
  t->SetBranchAddress("m2xpdist",&b.m2xpdist);
}


// ----------------------------------------------------------------------
void tmva1::train(string oname, string filename) {
   // This loads the library
   TMVA::Tools::Instance();
   
   (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40; 
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 40; 
   
   // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName(Form("%s.root", oname.c_str()));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   
   TH1D *hSetup = new TH1D("hSetup", "hSetup", 100, 0., 100.); 
   int i(0); 
   i =  1; hSetup->SetBinContent(i, fTrainAntiMuon?1:0); hSetup->GetXaxis()->SetBinLabel(i, "antimuon");
   i =  3; hSetup->SetBinContent(i, fRsigma); hSetup->GetXaxis()->SetBinLabel(i, "rsigma");
   i =  5; hSetup->SetBinContent(i, fApplyOn0?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn0");
   i =  6; hSetup->SetBinContent(i, fApplyOn1?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn1");
   i =  7; hSetup->SetBinContent(i, fApplyOn2?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn2");
   i = 10; hSetup->SetBinContent(i, fBdtSetup.NTrees); hSetup->GetXaxis()->SetBinLabel(i, "NTrees");
   i = 11; hSetup->SetBinContent(i, fBdtSetup.nEventsMin); hSetup->GetXaxis()->SetBinLabel(i, "nEventsMin");
   i = 12; hSetup->SetBinContent(i, fBdtSetup.nCuts); hSetup->GetXaxis()->SetBinLabel(i, "nCuts");
   i = 13; hSetup->SetBinContent(i, fBdtSetup.AdaBoostBeta); hSetup->GetXaxis()->SetBinLabel(i, "AdaBoostBeta");

   i = 20; hSetup->SetBinContent(i, fBdtSetup.MaxDepth); hSetup->GetXaxis()->SetBinLabel(i, "MaxDepth");
   i = 21; hSetup->SetBinContent(i, fBdtSetup.NNodesMax); hSetup->GetXaxis()->SetBinLabel(i, "NNodesMax");

   cout << "----------------------------------------------------------------------" << endl;
   cout << "==> oname: " << oname << " antimuon: " << fTrainAntiMuon <<  endl;

   string optstring = "V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
   optstring        = "V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
   cout << "==> Factory: " << optstring << endl;
   TMVA::Factory *factory = new TMVA::Factory(Form("%s", oname.c_str()), outputFile,  optstring.c_str());

//    factory->AddVariable("m1pt",       'F' );
//    factory->AddVariable("m1eta",      'F' );
//    factory->AddVariable("m2pt",       'F' );
//    factory->AddVariable("m2eta",      'F' );
   factory->AddVariable("pt",         'F' );
   factory->AddVariable("eta",        'F' );
   factory->AddVariable("fls3d",      'F' );
   factory->AddVariable("alpha",      'F' );
   factory->AddVariable("maxdoca",    'F' );
   factory->AddVariable("pvip",       'F' );
   factory->AddVariable("pvips",      'F' );
   factory->AddVariable("iso",        'F' );
   factory->AddVariable("docatrk",    'F' );
   factory->AddVariable("closetrk",   'I' );
   factory->AddVariable("chi2dof := chi2/dof",    'F' );

   factory->AddVariable("m1iso",   'F' );
   factory->AddVariable("m2iso",   'F' );
   factory->AddVariable("pvdchi2", 'F' );


   factory->AddSpectator("m",  "mass", "GeV", 'F' );
   
   TFile* inFile;
   TTree *applySg(0), *trainSg(0), *testSg(0), *applyBg(0), *trainBg(0), *testBg(0); 

   inFile = TFile::Open(filename.c_str());
   if (fApplyOn0) {
     cout << "==============> Apply on events0, train on events1, test on events2" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
     cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
   } 


   if (fApplyOn1) {
     cout << "==============> Apply on events1, train on events2, test on events0" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
     cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
   } 

   if (fApplyOn2) {
     cout << "==============> Apply on events2, train on events0, test on events1" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
     cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
     cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
   } 

   i = 30; hSetup->SetBinContent(i, applySg->GetEntries()); hSetup->GetXaxis()->SetBinLabel(i, "sgcnt");
   i = 31; hSetup->SetBinContent(i, applyBg->GetEntries()); hSetup->GetXaxis()->SetBinLabel(i, "bgcnt");
   writeOut(outputFile, hSetup); 

   Double_t signalWeight      = 1.; //= LUMISCALE; // 0.000388
   Double_t rbackgroundWeight = 1.;
   Double_t cbackgroundWeight = 1.;
   Double_t tbackgroundWeight = cbackgroundWeight;

   cout << "--> signal weight:     " << signalWeight << endl;
   cout << "--> cbackground weight: " << cbackgroundWeight << endl;
   cout << "--> rbackground weight: " << rbackgroundWeight << endl;

   factory->AddTree(trainSg,     "Signal",     signalWeight,  "", "train");
   factory->AddTree(testSg,      "Signal",     signalWeight,  "", "test");
   factory->AddTree(trainBg, "Background", cbackgroundWeight, "", "train");
   factory->AddTree(testBg,  "Background", tbackgroundWeight, "", "test");

   int nSgTrain = trainSg->GetEntries();
   int nSgTest  = testSg->GetEntries();

   int nBgTrain = trainBg->GetEntries();
   int nBgTest  = testBg->GetEntries();

   //   optstring = Form("nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:V"); 
   optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V", 
		    nSgTrain, nSgTest, nBgTrain, nBgTest); 
   cout << "==> PrepareTrainingAndTestTree: " << optstring << endl;
   factory->PrepareTrainingAndTestTree("", "", optstring.c_str());
   
   if (1) {
     optstring = Form("!H:V:NTrees=%d:nEventsMin=%d", fBdtSetup.NTrees, fBdtSetup.nEventsMin);
     optstring += Form(":BoostType=AdaBoost:AdaBoostBeta=%f:SeparationType=GiniIndex:nCuts=%d:PruneMethod=NoPruning", 
		       fBdtSetup.AdaBoostBeta, fBdtSetup.nCuts);

     //      fBdtSetup.MaxDepth  = 10;
     //      fBdtSetup.NNodesMax = 10;
     optstring += Form(":MaxDepth=%d:NNodesMax=%d", fBdtSetup.MaxDepth, fBdtSetup.NNodesMax);
   } else {
     optstring = "";
   }

// -- Josh's (modified) proposal
//    optstring = "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1";
//    optstring += ":UseBaggedGrad=F:nCuts=200:MaxDepth=3:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=1000:";

   cout << "==> BookMethod: " << optstring << endl;
   factory->BookMethod( TMVA::Types::kBDT, "BDT", optstring);

   cout << "==> TrainAllMethods " << endl;
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   cout << "==> TestAllMethods " << endl;
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   cout << "==> EvaluateAllMethods " << endl;
   factory->EvaluateAllMethods();

   // Save the output
   outputFile->Close();

   cout << "==> Wrote root file: " << outputFile->GetName() << endl;
   cout << "==> TMVAClassification is done!" << endl;

   delete factory;

}



// ----------------------------------------------------------------------
void tmva1::apply(const char *fname) {

  // --- Book the MVA methods
  string methodName("BDT");
  string dir("weights");
  string XmlName = Form("%s/%s-Events0_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd)); 
  XmlName = Form("%s/%s-Events1_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd)); 
  XmlName = Form("%s/%s-Events2_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd)); 

  // -- open files
  TFile *dfile(0);
  dfile = TFile::Open(fInputFiles.dname.c_str()); 
  if (!dfile) {
    cout << "ERROR: could not open data file" << endl;
    exit(1);
  }

  TFile *sfile(0);
  sfile = TFile::Open(fInputFiles.sname.c_str()); 
  if (!sfile) {
    cout << "ERROR: could not open signal file" << endl;
    exit(1);
  }
  
  TFile f(Form("%s-combined.root", fname), "RECREATE"); 

  TTree *tt = new TTree("bdtTree", "bdtTree");
  int classID;
  double w8; 
  tt->Branch("m",       &ftd.m,     "m/D");
  tt->Branch("eta",     &ftd.eta,   "eta/D");
  tt->Branch("weight",  &w8,        "weight/D");
  tt->Branch("bdt",     &fBDT,      "bdt/D");
  tt->Branch("bdt0",    &fBDT0,     "bdt0/D");
  tt->Branch("bdt1",    &fBDT1,     "bdt1/D");
  tt->Branch("bdt2",    &fBDT2,     "bdt2/D");
  tt->Branch("classID", &classID,   "classID/I");
  tt->Branch("hlt",     &ftd.hlt,   "hlt/O");
  tt->Branch("gmuid",   &ftd.gmuid, "gmuid/O");
  tt->Branch("evt",     &ftd.evt,   "evt/I");

  TH1D *hd = new TH1D("bdtTreeData", "", 20, 0., 20.);
  TH1D *hs = new TH1D("bdtTreeSignal", "", 20, 0., 20.);
 

  // -- data processing
  TTree* t = (TTree*)dfile->Get("candAnaMuMu/events");
  setupTree(t, ftd); 
  w8 = 1.;
  classID = 1;
  cout << "--- Processing data: " << t->GetEntries() << " events" << endl;
  int nEvent = t->GetEntries();
  double lostEvents(0); 
  double totalEvents(0); 
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing data event: " << ievt << std::endl;
    t->GetEntry(ievt);
    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8; 
    if (ftd.gmuid && ftd.hlt && preselection(ftd, fChannel)) {
      calcBDT();
      tt->Fill();
    } else {
      lostEvents += w8;
      hd->Fill(TMath::Abs(ftd.evt%3), w8); 
      hd->Fill(9, w8); 
    }
  }
  hd->SetBinContent(10, lostEvents); 
  cout << "lost events: " << lostEvents << " out of " << totalEvents << " events in total" << endl;

  // -- signal MC processing
  t = (TTree*)sfile->Get("candAnaMuMu/events");
  classID = 0;
  setupTree(t, ftd); 
  w8 = LUMISCALE;
  cout << "--- Processing signal: " << t->GetEntries() << " events" << endl;
  nEvent = t->GetEntries();
  lostEvents = 0; 
  totalEvents = 0; 
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing signal event: " << ievt << std::endl;
    t->GetEntry(ievt);

    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8; 
    if (ftd.gmuid && ftd.hlt && preselection(ftd, fChannel)) {
      calcBDT();
      tt->Fill();
    } else {
      lostEvents += w8;
      hs->Fill(TMath::Abs(ftd.evt%3), w8); 
      hs->Fill(9, w8); 
    }
  }
  hs->SetBinContent(10, lostEvents); 
  cout << "lost events: " << lostEvents << " out of " << totalEvents << " events in total" << endl;

  hd->Write();
  hs->Write();
  f.Write();
  f.Close();

}


// ----------------------------------------------------------------------
void tmva1::reAnalyze(int imin, int imax) {

  string name;
  for (int i = imin; i <= imax; ++i) {
    name = Form("TMVA-%d", i);
    cout << name << endl;
    analyze(name.c_str()); 
  }

}


// ----------------------------------------------------------------------
void tmva1::analyze(const char *fname) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  int bdtBins(200); 
  double bdtMin(-1.0), bdtMax(1.0); 

  TH1D *sm = new TH1D("sm", "m(signal)", 100, 4.9, 5.9); sm->Sumw2();
  TH1D *dm = new TH1D("dm", "m(data)", 100, 4.9, 5.9); dm->Sumw2();

  TH1D *sm0 = new TH1D("sm0", "m(signal)", 100, 4.9, 5.9); sm0->Sumw2();
  TH1D *dm0 = new TH1D("dm0", "m(data)", 100, 4.9, 5.9); dm0->Sumw2();

  TH1D *sm1 = new TH1D("sm1", "m(signal)", 100, 4.9, 5.9); sm1->Sumw2();
  TH1D *dm1 = new TH1D("dm1", "m(data)", 100, 4.9, 5.9); dm1->Sumw2();

  TH1D *sm2 = new TH1D("sm2", "m(signal)", 100, 4.9, 5.9); sm2->Sumw2();
  TH1D *dm2 = new TH1D("dm2", "m(data)", 100, 4.9, 5.9); dm2->Sumw2();

  cout << "Open " << Form("%s-combined.root", fname) << endl;
  TFile *f = TFile::Open(Form("%s-combined.root", fname), "UPDATE"); 
  TH1D *hd = (TH1D*)f->Get("bdtTreeData"); 
  double dLostEvents  = hd->GetBinContent(10); 
  double dLostEvents0 = hd->GetBinContent(1); 
  double dLostEvents1 = hd->GetBinContent(2); 
  double dLostEvents2 = hd->GetBinContent(3); 
  TH1D *hs = (TH1D*)f->Get("bdtTreeSignal"); 
  double sLostEvents  = hs->GetBinContent(10); 
  double sLostEvents0 = hs->GetBinContent(1); 
  double sLostEvents1 = hs->GetBinContent(2); 
  double sLostEvents2 = hs->GetBinContent(3); 

  TH1D *h = new TH1D("s1", "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); h->Sumw2();
  setTitles(h, "b > ", "S/#sqrt{S+B}"); 
  TH1D *hs2 = new TH1D("s2", "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); hs2->Sumw2();
  setTitles(hs2, "b > ", "S/#sqrt{S+B'}"); 
  TGraph *grocop = new TGraph(1); grocop->SetMarkerStyle(34); grocop->SetMarkerSize(3.); grocop->SetMarkerColor(kBlack);
  TGraph *groc   = new TGraph(bdtBins);   groc->SetMarkerStyle(24);  groc->SetMarkerSize(1.5);    groc->SetMarkerColor(kBlack); 
  TGraph *groc0  = new TGraph(bdtBins);   groc0->SetMarkerStyle(20); groc0->SetMarkerSize(0.8);   groc0->SetMarkerColor(kMagenta); 
  TGraph *groc1  = new TGraph(bdtBins);   groc1->SetMarkerStyle(20); groc1->SetMarkerSize(0.8);   groc1->SetMarkerColor(kRed); 
  TGraph *groc2  = new TGraph(bdtBins);   groc2->SetMarkerStyle(20); groc2->SetMarkerSize(0.8);   groc2->SetMarkerColor(kBlue); 

  TH1D *hroc = new TH1D("hroc", "", 200, 0., 1.); 
  TH1D *hroc0 = new TH1D("hroc0", "", 200, 0., 1.); 
  TH1D *hroc1 = new TH1D("hroc1", "", 200, 0., 1.); 
  TH1D *hroc2 = new TH1D("hroc2", "", 200, 0., 1.); 
  TFile *f0 = TFile::Open(Form("%s-Events0.root", fname)); 
  TFile *f1 = TFile::Open(Form("%s-Events1.root", fname)); 
  TFile *f2 = TFile::Open(Form("%s-Events2.root", fname)); 

  if (!f0 || !f1 || !f2) {
    cout << "ERROR: could not open a file" << endl;
    exit(1);
  }

  TH1D *hr01 = getRanking(fname, "IdTransformation", "events0");
  hr01->SetDirectory(f); 
  TH1D *hr02 =  getRanking(fname, "BDT", "events0");
  hr02->SetDirectory(f); 

  TH1D *hr11 = getRanking(fname, "IdTransformation", "events1");
  hr11->SetDirectory(f); 
  TH1D *hr12 =  getRanking(fname, "BDT", "events1");
  hr12->SetDirectory(f); 

  TH1D *hr21 = getRanking(fname, "IdTransformation", "events2");
  hr21->SetDirectory(f); 
  TH1D *hr22 =  getRanking(fname, "BDT", "events2");
  hr22->SetDirectory(f); 

  TH1F *trainBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT0->SetLineColor(kBlack);
  TH1F *testBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT0->SetLineColor(kBlack);

  TH1F *trainBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT1->SetLineColor(kRed);
  TH1F *testBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT1->SetLineColor(kRed);

  TH1F *trainBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT2->SetLineColor(kBlue);
  TH1F *testBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT2->SetLineColor(kBlue);

  TH1D *bBDT = new TH1D("bBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *cBDT = new TH1D("cBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *h2 = new TH1D("res_ssb", "max(S/Sqrt(S+B))", 10, 0., 10.); h2->Sumw2();  h2->SetDirectory(f); 

  // -- new rebinned versions of the BDT output plots. These ARE identical to the TMVA output IF the binning IS the same!
  TH1D *ap0bgBDT = new TH1D("ap0bgBDT", "",    400, -1., 1.);
  TH1D *tr0bgBDT = new TH1D("tr0bgBDT", "",    400, -1., 1.);
  TH1D *te0bgBDT = new TH1D("te0bgBDT", "",    400, -1., 1.);
  ap0bgBDT->SetLineColor(kBlack); te0bgBDT->SetLineColor(kBlack); tr0bgBDT->SetLineColor(kBlack);
  ap0bgBDT->SetMarkerColor(kBlack); te0bgBDT->SetMarkerColor(kBlack); tr0bgBDT->SetMarkerColor(kBlack);
  
  TH1D *ap1bgBDT = new TH1D("ap1bgBDT", "",    400, -1., 1.);
  TH1D *tr1bgBDT = new TH1D("tr1bgBDT", "",    400, -1., 1.);
  TH1D *te1bgBDT = new TH1D("te1bgBDT", "",    400, -1., 1.);
  ap1bgBDT->SetLineColor(kRed); te1bgBDT->SetLineColor(kRed); tr1bgBDT->SetLineColor(kRed);
  ap1bgBDT->SetMarkerColor(kRed); te1bgBDT->SetMarkerColor(kRed); tr1bgBDT->SetMarkerColor(kRed);

  TH1D *ap2bgBDT = new TH1D("ap2bgBDT", "",    400, -1., 1.);
  TH1D *tr2bgBDT = new TH1D("tr2bgBDT", "",    400, -1., 1.);
  TH1D *te2bgBDT = new TH1D("te2bgBDT", "",    400, -1., 1.);
  ap2bgBDT->SetLineColor(kBlue); te2bgBDT->SetLineColor(kBlue); tr2bgBDT->SetLineColor(kBlue);
  ap2bgBDT->SetMarkerColor(kBlue); te2bgBDT->SetMarkerColor(kBlue); tr2bgBDT->SetMarkerColor(kBlue);

  TH1D *ap0sgBDT = new TH1D("ap0sgBDT", "",    400, -1., 1.);
  TH1D *tr0sgBDT = new TH1D("tr0sgBDT", "",    400, -1., 1.);
  TH1D *te0sgBDT = new TH1D("te0sgBDT", "",    400, -1., 1.);
  ap0sgBDT->SetLineColor(kBlack); te0sgBDT->SetLineColor(kBlack); tr0sgBDT->SetLineColor(kBlack);
  ap0sgBDT->SetMarkerColor(kBlack); te0sgBDT->SetMarkerColor(kBlack); tr0sgBDT->SetMarkerColor(kBlack);
  
  TH1D *ap1sgBDT = new TH1D("ap1sgBDT", "",    400, -1., 1.);
  TH1D *tr1sgBDT = new TH1D("tr1sgBDT", "",    400, -1., 1.);
  TH1D *te1sgBDT = new TH1D("te1sgBDT", "",    400, -1., 1.);
  ap1sgBDT->SetLineColor(kRed); te1sgBDT->SetLineColor(kRed); tr1sgBDT->SetLineColor(kRed);
  ap1sgBDT->SetMarkerColor(kRed); te1sgBDT->SetMarkerColor(kRed); tr1sgBDT->SetMarkerColor(kRed);

  TH1D *ap2sgBDT = new TH1D("ap2sgBDT", "",    400, -1., 1.);
  TH1D *tr2sgBDT = new TH1D("tr2sgBDT", "",    400, -1., 1.);
  TH1D *te2sgBDT = new TH1D("te2sgBDT", "",    400, -1., 1.);
  ap2sgBDT->SetLineColor(kBlue); te2sgBDT->SetLineColor(kBlue); tr2sgBDT->SetLineColor(kBlue);
  ap2sgBDT->SetMarkerColor(kBlue); te2sgBDT->SetMarkerColor(kBlue); tr2sgBDT->SetMarkerColor(kBlue);

  
  TTree *t = (TTree*)f->Get("bdtTree");
  double bdt, m, w8; 
  double bdt0, bdt1, bdt2; 
  int classID, evt;
  bool gmuid, hlt;
  t->SetBranchAddress("bdt", &bdt);
  t->SetBranchAddress("bdt0", &bdt0);
  t->SetBranchAddress("bdt1", &bdt1);
  t->SetBranchAddress("bdt2", &bdt2);
  t->SetBranchAddress("classID", &classID);
  t->SetBranchAddress("m", &m);
  t->SetBranchAddress("weight", &w8);
  t->SetBranchAddress("hlt", &hlt);
  t->SetBranchAddress("gmuid", &gmuid);
  t->SetBranchAddress("evt", &evt);
  
  // -- data (overall) distribution
  TCanvas *c0 = getC0();
  int nEvent(0); 
  nEvent = t->GetEntries();
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    t->GetEntry(ievt);
    if (1 == classID) {
      if (hlt && gmuid) bBDT->Fill(bdt);
      if (hlt && !gmuid) cBDT->Fill(bdt);

      if (!hlt) continue;
      if (!gmuid) continue;

      ap0bgBDT->Fill(bdt0);
      ap1bgBDT->Fill(bdt1);
      ap2bgBDT->Fill(bdt2);

      if (evt%3==0) {
	te1bgBDT->Fill(bdt1);
	tr2bgBDT->Fill(bdt2);
      }
      if (evt%3==1) {
	te2bgBDT->Fill(bdt2);
	tr0bgBDT->Fill(bdt0);
      }
      if (evt%3==2) {
	te0bgBDT->Fill(bdt0);
	tr1bgBDT->Fill(bdt1);
      }
    }

    if (0 == classID) {
      if (!hlt) continue;
      if (!gmuid) continue;

      ap0sgBDT->Fill(bdt0);
      ap1sgBDT->Fill(bdt1);
      ap2sgBDT->Fill(bdt2);

      if (evt%3==0) {
	te1sgBDT->Fill(bdt1);
	tr2sgBDT->Fill(bdt2);
      }
      if (evt%3==1) {
	te2sgBDT->Fill(bdt2);
	tr0sgBDT->Fill(bdt0);
      }
      if (evt%3==2) {
	te0sgBDT->Fill(bdt0);
	tr1sgBDT->Fill(bdt1);
      }
    }
  }

  double bgks00   = tr0bgBDT->KolmogorovTest(te0bgBDT); 
  double bgks11   = tr1bgBDT->KolmogorovTest(te1bgBDT); 
  double bgks22   = tr2bgBDT->KolmogorovTest(te2bgBDT); 

  double bgks01te = te0bgBDT->KolmogorovTest(te1bgBDT);
  double bgks01tr = tr0bgBDT->KolmogorovTest(tr1bgBDT);

  double bgks12te = te1bgBDT->KolmogorovTest(te2bgBDT);
  double bgks12tr = tr1bgBDT->KolmogorovTest(tr2bgBDT);

  double bgks20te = te2bgBDT->KolmogorovTest(te0bgBDT);
  double bgks20tr = tr2bgBDT->KolmogorovTest(tr0bgBDT);

  trainBDT0->Scale(1./trainBDT0->GetSumOfWeights());
  testBDT0->Scale(1./testBDT0->GetSumOfWeights());
  trainBDT1->Scale(1./trainBDT1->GetSumOfWeights());
  testBDT1->Scale(1./testBDT1->GetSumOfWeights());
  trainBDT2->Scale(1./trainBDT2->GetSumOfWeights());
  testBDT2->Scale(1./testBDT2->GetSumOfWeights());

  ap0bgBDT->Scale(1./ap0bgBDT->GetSumOfWeights());
  tr0bgBDT->Scale(1./tr0bgBDT->GetSumOfWeights());
  te0bgBDT->Scale(1./te0bgBDT->GetSumOfWeights());

  ap1bgBDT->Scale(1./ap1bgBDT->GetSumOfWeights());
  tr1bgBDT->Scale(1./tr1bgBDT->GetSumOfWeights());
  te1bgBDT->Scale(1./te1bgBDT->GetSumOfWeights());

  ap2bgBDT->Scale(1./ap2bgBDT->GetSumOfWeights());
  tr2bgBDT->Scale(1./tr2bgBDT->GetSumOfWeights());
  te2bgBDT->Scale(1./te2bgBDT->GetSumOfWeights());


  double sgks00   = tr0sgBDT->KolmogorovTest(te0sgBDT); 
  double sgks11   = tr1sgBDT->KolmogorovTest(te1sgBDT); 
  double sgks22   = tr2sgBDT->KolmogorovTest(te2sgBDT); 

  double sgks01te = te0sgBDT->KolmogorovTest(te1sgBDT);
  double sgks01tr = tr0sgBDT->KolmogorovTest(tr1sgBDT);

  double sgks12te = te1sgBDT->KolmogorovTest(te2sgBDT);
  double sgks12tr = tr1sgBDT->KolmogorovTest(tr2sgBDT);

  double sgks20te = te2sgBDT->KolmogorovTest(te0sgBDT);
  double sgks20tr = tr2sgBDT->KolmogorovTest(tr0sgBDT);

  ap0sgBDT->Scale(1./ap0sgBDT->GetSumOfWeights());
  tr0sgBDT->Scale(1./tr0sgBDT->GetSumOfWeights());
  te0sgBDT->Scale(1./te0sgBDT->GetSumOfWeights());

  ap1sgBDT->Scale(1./ap1sgBDT->GetSumOfWeights());
  tr1sgBDT->Scale(1./tr1sgBDT->GetSumOfWeights());
  te1sgBDT->Scale(1./te1sgBDT->GetSumOfWeights());

  ap2sgBDT->Scale(1./ap2sgBDT->GetSumOfWeights());
  tr2sgBDT->Scale(1./tr2sgBDT->GetSumOfWeights());
  te2sgBDT->Scale(1./te2sgBDT->GetSumOfWeights());

  double hmax = tr0bgBDT->GetMaximum(); 
  if (tr1bgBDT->GetMaximum() > hmax) hmax = tr1bgBDT->GetMaximum();
  if (tr2bgBDT->GetMaximum() > hmax) hmax = tr2bgBDT->GetMaximum();
  tr0bgBDT->SetMaximum(1.3*hmax);

  tr0bgBDT->Draw("p");
  te0bgBDT->Draw("same");

  tr1bgBDT->Draw("psame");
  te1bgBDT->Draw("same");

  tr2bgBDT->Draw("psame");
  te2bgBDT->Draw("same");

  c0->SaveAs(Form("plots/%s-rebinned-bg-overlays.pdf", fname)); 

  hmax = tr0sgBDT->GetMaximum(); 
  if (tr1sgBDT->GetMaximum() > hmax) hmax = tr1sgBDT->GetMaximum();
  if (tr2sgBDT->GetMaximum() > hmax) hmax = tr2sgBDT->GetMaximum();
  tr0sgBDT->SetMaximum(1.3*hmax);

  tr0sgBDT->Draw("p");
  te0sgBDT->Draw("same");

  tr1sgBDT->Draw("psame");
  te1sgBDT->Draw("same");

  tr2sgBDT->Draw("psame");
  te2sgBDT->Draw("same");

  c0->SaveAs(Form("plots/%s-rebinned-bg-overlays.pdf", fname));


  writeOut(f, ap0bgBDT); 
  writeOut(f, tr0bgBDT); 
  writeOut(f, te0bgBDT); 

  writeOut(f, ap1bgBDT); 
  writeOut(f, tr1bgBDT); 
  writeOut(f, te1bgBDT); 

  writeOut(f, ap2bgBDT); 
  writeOut(f, tr2bgBDT); 
  writeOut(f, te2bgBDT); 

  writeOut(f, ap0sgBDT); 
  writeOut(f, tr0sgBDT); 
  writeOut(f, te0sgBDT); 

  writeOut(f, ap1sgBDT); 
  writeOut(f, tr1sgBDT); 
  writeOut(f, te1sgBDT); 

  writeOut(f, ap2sgBDT); 
  writeOut(f, tr2sgBDT); 
  writeOut(f, te2sgBDT); 

  // -- compute S and B, SSB
  double bdtCut, maxSSB(-1.), maxBDT(-1.), sCnt(0.), dCnt(0.); 
  double maxSSBsimple(-1.), maxBDTsimple(-1.);
  double sCnt0(0.), dCnt0(0.), sCnt1(0.), dCnt1(0.), sCnt2(0.), dCnt2(0.);
  double seffTot(0.), seffMax(0.), deffMax(0.); 
  double seffBinD99(0.), sMax(-1.), bMax(-1);
  int ibin(0), gbin(0), dMax(-1), dhiMax(-1); 
  for (ibin = bdtBins; ibin >=0; --ibin) {
    bdtCut = bdtMin + ibin*(bdtMax-bdtMin)/bdtBins;
    if (0 == ibin%10) cout << " bin " << ibin << " cutting at bdt > " << bdtCut << endl;
    nEvent = t->GetEntries();
    sm->Reset();    dm->Reset();
    sm0->Reset();   dm0->Reset();
    sm1->Reset();   dm1->Reset();
    sm2->Reset();   dm2->Reset();
    sCnt = dCnt = 0.;
    sCnt0 = dCnt0 = sCnt1 = dCnt1 = sCnt2 = dCnt2 = 0.;
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      t->GetEntry(ievt);

      if (false == hlt) continue;
      if (false == gmuid) continue;
      if (0 == classID) {
	sCnt += w8; 
	if (0 == TMath::Abs(evt%3)) sCnt0 += w8;
	if (1 == TMath::Abs(evt%3)) sCnt1 += w8;
	if (2 == TMath::Abs(evt%3)) sCnt2 += w8;
      } else {
	dCnt += w8; 
	if (0 == TMath::Abs(evt%3)) dCnt0 += w8;
	if (1 == TMath::Abs(evt%3)) dCnt1 += w8;
	if (2 == TMath::Abs(evt%3)) dCnt2 += w8;
      }

      if (bdt <= bdtCut) continue;

      if (0 == classID) {
	sm->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) sm0->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) sm1->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) sm2->Fill(m, w8); 
      } else {
	dm->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) dm0->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) dm1->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) dm2->Fill(m, w8); 
      }
    }

    //    cout << "dcnt = " << dCnt << "/" << dCnt0+dCnt1+dCnt2 << " 0 = " << dCnt0 << " 1 = " << dCnt1 << " 2 = " << dCnt2 << endl;

    double s = sm->Integral(sm->FindBin(5.3), sm->FindBin(5.45));
    double d = dm->Integral(0, dm->GetNbinsX());
    double seff = s/sCnt;
    double bsimple = d*(5.45-5.30)/(5.9-4.9-0.25);
    double deff = bsimple/(dLostEvents + dCnt);

    double dhi = dm->Integral(dm->FindBin(5.45), dm->GetNbinsX());
    double pbg = 0.07*s;
    double b = dhi*(5.45-5.30)/(5.9-5.45);

    if ((1.-deff) > 0.999998) seffBinD99 = seff;

    groc->SetPoint(gbin, seff, 1.-deff); 
    hroc->SetBinContent(hroc->FindBin(seff), 1.-deff); 

    double s0 = sm0->Integral(sm0->FindBin(5.3), sm0->FindBin(5.45));
    double d0 = dm0->Integral(0, dm0->GetNbinsX());
    double b0 = d0*(5.45-5.30)/(5.9-4.9-0.25);
    double seff0 = s0/sCnt0;
    double deff0 = b0/(dLostEvents0 + dCnt0);

    groc0->SetPoint(gbin, seff0, 1.-deff0); 
    hroc0->SetBinContent(hroc0->FindBin(seff0), 1.-deff0); 

    double s1 = sm1->Integral(sm1->FindBin(5.3), sm1->FindBin(5.45));
    double d1 = dm1->Integral(1, dm1->GetNbinsX());
    double b1 = d1*(5.45-5.30)/(5.9-4.9-0.25);
    double seff1 = s1/sCnt1;
    double deff1 = b1/(dLostEvents1 + dCnt1);

    groc1->SetPoint(gbin, seff1, 1.-deff1); 
    hroc1->SetBinContent(hroc1->FindBin(seff1), 1.-deff1); 

    double s2 = sm2->Integral(sm2->FindBin(5.3), sm2->FindBin(5.45));
    double d2 = dm2->Integral(1, dm2->GetNbinsX());
    double b2 = d2*(5.45-5.30)/(5.9-4.9-0.25);
    double seff2 = s2/sCnt2;
    double deff2 = b2/(dLostEvents2 + dCnt2);

    groc2->SetPoint(gbin, seff2, 1.-deff2); 
    hroc2->SetBinContent(hroc2->FindBin(seff2), 1.-deff2); 
    //    cout << "roc2: " << seff2 << " " << 1.-deff2 << endl;

    ++gbin;
    if (s+b+pbg >0) {
      double ssb = s/TMath::Sqrt(s+b+pbg);
      if (ssb > maxSSB) {
	sMax = s;
	dMax = d; 
	dhiMax = dhi; 
	bMax = b+pbg; 
	maxSSB = ssb; 
	maxBDT = bdtCut;
	seffMax = seff;
	deffMax = deff;
	seffTot = s/(sLostEvents + sCnt); 
      }
      h->SetBinContent(ibin, ssb); 
    } else {
      h->SetBinContent(ibin, 0); 
    }
    if (s+bsimple >0) {
      double ssb = s/TMath::Sqrt(s+bsimple);
      if (ssb > maxSSBsimple) {
	maxSSBsimple = ssb; 
	maxBDTsimple = bdtCut;
      }
      hs2->SetBinContent(ibin, ssb); 
    } else {
      hs2->SetBinContent(ibin, 0); 
    }
    //    cout << "S = " << s << " B = " << b << " => S/sqrt(S+B) = " << s/TMath::Sqrt(s+b) << endl;
    c0->Clear();
    dm->SetTitle(Form("bdt >%f, S/B = %5.2f/%5.2f = %4.3f   D = %5.0f", bdtCut, s, b, s/b, d));
    dm->Draw("e");
    sm->Draw("samehist");
    c0->Modified();
    c0->Update();
  }

  grocop->SetPoint(0, seffMax, 1.-deffMax); 

  // -- patch hroc for empty bins
  double okVal(1.); 
  for (int i = 1; i< hroc->GetNbinsX(); ++i) {
    if (hroc->GetBinContent(i) < 1.e-5) hroc->SetBinContent(i, okVal); 
    okVal = hroc->GetBinContent(i);
  }

  okVal = 1.;
  for (int i = 1; i< hroc0->GetNbinsX(); ++i) {
    if (hroc0->GetBinContent(i) < 1.e-5) hroc0->SetBinContent(i, okVal); 
    okVal = hroc0->GetBinContent(i);
  }

  okVal = 1.;
  for (int i = 1; i< hroc1->GetNbinsX(); ++i) {
    if (hroc1->GetBinContent(i) < 1.e-5) hroc1->SetBinContent(i, okVal); 
    okVal = hroc1->GetBinContent(i);
  }

  okVal = 1.;
  for (int i = 1; i< hroc2->GetNbinsX(); ++i) {
    if (hroc2->GetBinContent(i) < 1.e-5) hroc2->SetBinContent(i, okVal); 
    okVal = hroc2->GetBinContent(i);
  }

  for (int i = hroc->GetNbinsX()-1; i > 0; --i) {
    if (!(hroc->GetBinContent(i) < hroc->GetBinContent(i-1))) {
      hroc->SetBinContent(i, 0); 
    } else {
      break;
    }
  }

  for (int i = hroc0->GetNbinsX()-1; i > 0; --i) {
    if (!(hroc0->GetBinContent(i) < hroc0->GetBinContent(i-1))) {
      hroc0->SetBinContent(i, 0); 
    } else {
      break;
    }
  }

  for (int i = hroc1->GetNbinsX()-1; i > 0; --i) {
    if (!(hroc1->GetBinContent(i) < hroc1->GetBinContent(i-1))) {
      hroc1->SetBinContent(i, 0); 
    } else {
      break;
    }
  }

  for (int i = hroc2->GetNbinsX()-1; i > 0; --i) {
    if (!(hroc2->GetBinContent(i) < hroc2->GetBinContent(i-1))) {
      hroc2->SetBinContent(i, 0); 
    } else {
      break;
    }
  }

  bBDT->Scale(1./bBDT->GetSumOfWeights());
  cBDT->Scale(1./cBDT->GetSumOfWeights());

  // -- BG overlays
  // --------------
  c0->Clear();
  c0->Divide(1,2);
  // -- linear plot
  c0->cd(1);
  setHist(trainBDT0, kRed, 24); 
  setTitles(trainBDT0, "b", "(1/N) dN/db", 0.05, 1.2, 1.5); 

  hmax = trainBDT0->GetMaximum(); 
  if (trainBDT1->GetMaximum() > hmax) hmax = trainBDT1->GetMaximum();
  if (trainBDT2->GetMaximum() > hmax) hmax = trainBDT2->GetMaximum();

  trainBDT0->SetMaximum(1.2*hmax); 
  trainBDT0->DrawCopy("e");
  setHist(testBDT0, kRed); 
  testBDT0->Draw("histsame");

  setHist(trainBDT1, kBlue, 24); 
  trainBDT1->Draw("esame");
  setHist(testBDT1, kBlue); 
  testBDT1->Draw("histsame");

  setHist(trainBDT2, kBlack, 24); 
  trainBDT2->Draw("esame");
  setHist(testBDT2, kBlack); 
  testBDT2->Draw("histsame");

  bBDT->SetLineColor(kBlack);
  bBDT->Draw("histsame");
  cBDT->SetLineColor(kBlack);
  cBDT->SetLineStyle(kDashed);
  cBDT->Draw("histsame");

  newLegend(0.55, 0.5, 0.85, 0.85, "Samples"); 
  legg->AddEntry(trainBDT0, "0 train", "p"); 
  legg->AddEntry(testBDT0, "0 test", "l"); 
  legg->AddEntry(trainBDT1, "1 train", "p"); 
  legg->AddEntry(testBDT1, "1 test", "l"); 
  legg->AddEntry(trainBDT2, "2 train", "p"); 
  legg->AddEntry(testBDT2, "2 test", "l"); 
  legg->AddEntry(bBDT, "all events", "l"); 
  legg->AddEntry(cBDT, "all events + failed #mu ID", "l"); 
  legg->Draw();

  tl->SetTextSize(0.04);
  tl->DrawLatex(0.6, 0.80, Form("KS 0/1: %5.4f/%5.4f", bgks01te, bgks01tr)); 
  tl->DrawLatex(0.6, 0.72, Form("KS 1/2: %5.4f/%5.4f", bgks12te, bgks12tr)); 
  tl->DrawLatex(0.6, 0.64, Form("KS 2/0: %5.4f/%5.4f", bgks20te, bgks20tr)); 
  tl->SetTextSize(0.03);

  // -- logarithmic plot
  c0->cd(2);
  setHist(trainBDT0, kRed, 24); 
  trainBDT0->SetMaximum(1.5*hmax); 
  trainBDT0->Draw("e");
  setHist(testBDT0, kRed); 
  testBDT0->Draw("histsame");

  setHist(trainBDT1, kBlue, 24); 
  trainBDT1->Draw("esame");
  setHist(testBDT1, kBlue); 
  testBDT1->Draw("histsame");

  setHist(trainBDT2, kBlack, 24); 
  trainBDT2->Draw("esame");
  setHist(testBDT2, kBlack); 
  testBDT2->Draw("histsame");

  bBDT->SetLineColor(kBlack);
  bBDT->Draw("histsame");
  cBDT->SetLineColor(kBlack);
  cBDT->SetLineStyle(kDashed);
  cBDT->Draw("histsame");

  gPad->SetLogy(1);
  c0->SaveAs(Form("plots/%s-tmva-overlays.pdf", fname)); 

  c0->Clear();
  TH2F* frame(0); 
  frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.65, 100, 0.99999, 1.000001);
  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
  frame->Draw();  

  gPad->SetLogy(0);
  gPad->SetLeftMargin(0.15);
  groc->Draw("p"); 
  groc->GetXaxis()->SetTitle("#epsilon_{ S}"); 
  groc->GetYaxis()->SetTitle("1 - #epsilon_{ B}"); 
  double rocInt = hroc->Integral(1, hroc->GetNbinsX())*hroc->GetBinWidth(1); 
  double rocInt2= hroc->Integral(1, hroc->FindBin(seffBinD99))*hroc->GetBinWidth(1); 
  cout << " ==> seffBinD99 " << seffBinD99 << " in bin " << hroc->FindBin(seffBinD99) << endl;
  groc->SetName("groc"); 
  groc->SetTitle(Form("integral = %5.3f", rocInt));

  grocop->SetName("grocop"); 
  groc0->SetName("groc0"); 
  groc1->SetName("groc1"); 
  groc2->SetName("groc2"); 

  groc0->Draw("p");
  groc1->Draw("p");
  groc2->Draw("p");
  grocop->Draw("p");
  
  tl->DrawLatex(0.25, 0.44, Form("D/B/S = %d/%2.1f/%2.1f", dMax, bMax, sMax)); 
  tl->DrawLatex(0.25, 0.40, Form("S/#sqrt{S+B}(MC) = %4.3f (%4.3f)", maxSSB, maxSSBsimple)); 
  tl->DrawLatex(0.25, 0.36, Form("b_{max}(MC) = %4.3f (%4.3f)", maxBDT, maxBDTsimple)); 
  tl->DrawLatex(0.25, 0.32, Form("#epsilon_{BDT} = %4.3f", seffMax)); 
  tl->DrawLatex(0.25, 0.28, Form("#epsilon_{tot} = %6.5f", seffTot)); 
  tl->DrawLatex(0.25, 0.24, Form("I_{tot} = %6.5f", rocInt)); 
  tl->DrawLatex(0.25, 0.20, Form("I_{part} = %6.5f", rocInt2)); 

  string texname = string(fname) + ".tex";
  ofstream TEX(texname.c_str());
  
  TEX << Form("\\vdef{s%s:string}       {%s}", fname, fname) << endl;
  TEX << Form("\\vdef{s%s:ssb}          {%4.3f}", fname, maxSSB) << endl;
  TEX << Form("\\vdef{s%s:maxbdt}       {%4.3f}", fname, maxBDT) << endl;
  TEX << Form("\\vdef{s%s:ssbsimple}    {%4.3f}", fname, maxSSBsimple) << endl;
  TEX << Form("\\vdef{s%s:maxbdtsimple} {%4.3f}", fname, maxBDTsimple) << endl;
  TEX << Form("\\vdef{s%s:Smc}          {%4.3f}", fname, sMax) << endl;
  TEX << Form("\\vdef{s%s:D}            {%d}", fname, dMax) << endl;
  TEX << Form("\\vdef{s%s:Dhi}          {%d}", fname, dhiMax) << endl;
  TEX << Form("\\vdef{s%s:B}            {%4.3f}", fname, bMax) << endl;
  TEX << Form("\\vdef{s%s:ipart}        {%6.5f}", fname, rocInt2) << endl;
  TEX << Form("\\vdef{s%s:itot}         {%6.5f}", fname, rocInt) << endl;
  TEX << Form("\\vdef{s%s:epstot}       {%6.5f}", fname, seffTot) << endl;
  TEX << Form("\\vdef{s%s:epsbdt}       {%6.5f}", fname, seffMax) << endl;
  TEX.close();

  newLegend(0.22, 0.47, 0.55, 0.67); 

  legg->SetTextSize(0.035);  
  legg->AddEntry(grocop, "operating point", "p"); 
  legg->AddEntry(groc,  "combined", "p"); 
  legg->AddEntry(groc0,  "BDT 0", "p"); 
  legg->AddEntry(groc1,  "BDT 1", "p"); 
  legg->AddEntry(groc2,  "BDT 2", "p"); 
  legg->Draw();

  
  c0->SaveAs(Form("plots/%s-roc0.pdf", fname)); 
  //   gPad->SetLogy(1);
  //   groc->SetMinimum(1.e-6); 
  //   groc->Draw("ap"); 
  //   c0->SaveAs(Form("plots/%s-roc1.pdf", fname)); 

  c0->Clear();
  h->Draw();
  hs2->Draw("same");
  tl->DrawLatex(0.2, 0.85, Form("SSB_{max} = %4.3f (%4.3f)", maxSSB, maxSSBsimple)); 
  tl->DrawLatex(0.2, 0.80, Form("BDT_{max} > %4.3f (%4.3f)", maxBDT, maxBDTsimple)); 
  tl->DrawLatex(0.2, 0.75, Form("ROC_{int} = %4.3f", rocInt)); 
  c0->SaveAs(Form("plots/%s-ssb.pdf", fname)); 

  cout << "Write out SSB histograms" << endl;
  cout << "  maxSSB: " << maxSSB << " at BDT > " << maxBDT << endl;
  h->SetDirectory(f);

  h2->SetBinContent(1, maxSSB);  h2->GetXaxis()->SetBinLabel(1, "maxSSB");
  h2->SetBinContent(2, maxBDT);  h2->GetXaxis()->SetBinLabel(2, "maxBDT");
  h2->SetBinContent(3, rocInt);  h2->GetXaxis()->SetBinLabel(2, "rocInt");

  f->cd();
  grocop->Write();
  groc->Write();
  groc0->Write();
  groc1->Write();
  groc2->Write();

  f->Write();
  f->Close();

} 


// ----------------------------------------------------------------------
void tmva1::mvas(string fname) { //, HistType htype, Bool_t useTMVAStyle ) {
   // set style and remove existing canvas'
   TMVAGlob::Initialize( kTRUE );


   // checks if file with name "fin" is already open, and if not opens one
   TString fin = Form("%s.root", fname.c_str()); 
   TFile* file = TFile::Open(fin.Data(), "UPDATE");  
   TH1D *h2 = new TH1D("res_ks", "KS probabilities", 10, 0., 10.); h2->Sumw2();

   // define Canvas layout here!
   const Int_t width = 600;   // size of canvas

   // this defines how many canvases we need
   TCanvas *c = 0;

   // counter variables
   Int_t countCanvas = 0;

   // search for the right histograms in full list of keys
   TIter next(file->GetListOfKeys());
   TKey *key(0);   
   while ((key = (TKey*)next())) {

      if (!TString(key->GetName()).BeginsWith("Method_")) continue;
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      TString methodName;
      TMVAGlob::GetMethodName(methodName,key);

      TDirectory* mDir = (TDirectory*)key->ReadObj();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while ((titkey = (TKey*)keyIt())) {

         if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

         TDirectory *titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);

         cout << "--- Found directory for method: " << methodName << "::" << methodTitle << flush;
         TString hname = "MVA_" + methodTitle;
         TH1* sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
         TH1* bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));


         cout << " containing " << hname << "_S/_B" << endl;
         // chop off useless stuff
	 sig->SetTitle( Form("TMVA overtraining check for classifier: %s", methodTitle.Data()) );
         
         // create new canvas
         TString ctitle = Form("TMVA comparison %s",methodTitle.Data()) ;
         
         c = new TCanvas( Form("canvas%d", countCanvas+1), ctitle, 
                          countCanvas*50+200, countCanvas*20, width, (Int_t)width*0.78 ); 
    
         // set the histogram style
         TMVAGlob::SetSignalAndBackgroundStyle( sig, bgd );
         
         // normalise both signal and background
         TMVAGlob::NormalizeHists( sig, bgd );
         
         // frame limits (choose judicuous x range)
         Float_t nrms = 10;
         cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
         cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
         Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(), 
                                               bgd->GetMean() - nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmin() );
         Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(), 
                                               bgd->GetMean() + nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmax() );
         Float_t ymin = 0;
         Float_t maxMult = 2.0;
         Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;

	 xmin = -1.;
	 xmax = 1.;
   
         // build a frame
         Int_t nb = 500;
         TString hFrameName(TString("frame") + methodTitle);
         TObject *o = gROOT->FindObject(hFrameName);
         if(o) delete o;
         TH2F* frame = new TH2F( hFrameName, sig->GetTitle(), 
                                 nb, xmin, xmax, nb, ymin, ymax );
         frame->GetXaxis()->SetTitle( methodTitle + " response"  );
         frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
         TMVAGlob::SetFrameStyle( frame );
   
         // eventually: draw the frame
         frame->Draw();  
    
         c->GetPad(0)->SetLeftMargin( 0.105 );
         frame->GetYaxis()->SetTitleOffset( 1.2 );

         // Draw legend               
         TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12, 
                                       c->GetLeftMargin() +  0.40, 1 - c->GetTopMargin() );
         legend->SetFillStyle( 1 );
         legend->AddEntry(sig,TString("Signal")     + " (test sample)", "F");
         legend->AddEntry(bgd,TString("Background") + " (test sample)", "F");
         legend->SetBorderSize(1);
         legend->SetMargin(0.2);
         legend->Draw("same");

         // overlay signal and background histograms
         sig->Draw("samehist");
         bgd->Draw("samehist");
   
	 TH1* sigOv = 0;
	 TH1* bgdOv = 0;
	 
	 TString ovname = hname += "_Train";
	 sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
	 bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));
	 
	 if (sigOv == 0 || bgdOv == 0) {
	   cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
	 }
	 else {
	   cout << "--- Found comparison histograms for overtraining check" << endl;
	   
	   TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
					  1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
	   legend2->SetFillStyle( 1 );
	   legend2->SetBorderSize(1);
	   legend2->AddEntry(sigOv,"Signal (training sample)","P");
	   legend2->AddEntry(bgdOv,"Background (training sample)","P");
	   legend2->SetMargin( 0.1 );
	   legend2->Draw("same");
	 }
	 // normalise both signal and background
	 TMVAGlob::NormalizeHists( sigOv, bgdOv );
	 
	 Int_t col = sig->GetLineColor();
	 sigOv->SetMarkerColor( col );
	 sigOv->SetMarkerSize( 0.7 );
	 sigOv->SetMarkerStyle( 20 );
	 sigOv->SetLineWidth( 1 );
	 sigOv->SetLineColor( col );
	 sigOv->Draw("e1same");
	 
	 col = bgd->GetLineColor();
	 bgdOv->SetMarkerColor( col );
	 bgdOv->SetMarkerSize( 0.7 );
	 bgdOv->SetMarkerStyle( 20 );
	 bgdOv->SetLineWidth( 1 );
	 bgdOv->SetLineColor( col );
	 bgdOv->Draw("e1same");
	 
	 ymax = TMath::Max(ymax, 
			   TMath::Max( static_cast<Float_t>(sigOv->GetMaximum()), static_cast<Float_t>(bgdOv->GetMaximum()) )*maxMult);
	 frame->GetYaxis()->SetLimits( 0, ymax );
	 
	 // for better visibility, plot thinner lines
	 sig->SetLineWidth( 1 );
	 bgd->SetLineWidth( 1 );
	 
	 // perform K-S test
	 cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
	 Double_t kolS = sig->KolmogorovTest( sigOv );
	 Double_t kolB = bgd->KolmogorovTest( bgdOv );
	 cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;
	 
	 TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
	 TText* tt = new TText( 0.12, 0.91, probatext );
	 tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad(); 
	 
	 h2->SetBinContent(1, kolS); h2->GetXaxis()->SetBinLabel(1, "KS(sg)");
	 h2->SetBinContent(2, kolB); h2->GetXaxis()->SetBinLabel(2, "KS(bg)");
	 writeOut(file, h2); 
	 // h2->Write();
	 


         // redraw axes
         frame->Draw("sameaxis");
	 
         // text for overflows
         Int_t    nbin = sig->GetNbinsX();
         Double_t dxu  = sig->GetBinWidth(0);
         Double_t dxo  = sig->GetBinWidth(nbin+1);
         TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
                                sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
                                sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
         TText* t = new TText( 0.975, 0.115, uoflow );
         t->SetNDC();
         t->SetTextSize( 0.030 );
         t->SetTextAngle( 90 );
         t->AppendPad();    
   
         // update canvas
         c->Update();

         // save canvas to file
	 c->SaveAs(Form("plots/%s-overtrain0.pdf", fname.c_str())); 
	 frame->GetYaxis()->SetLimits(1.e-5, 5.e1);
	 c->SetLogy(1);
	 c->SaveAs(Form("plots/%s-overtrain1.pdf", fname.c_str())); 

         countCanvas++;
      }
      cout << "";
   }
   cout << endl;
   file->Write();
   file->Close();
}



// ----------------------------------------------------------------------
TH1D* tmva1::getRanking(string fname, string prefix, std::string type) {
  TH1D *h1 = new TH1D(Form("rank_%s_%s", type.c_str(), prefix.c_str()), Form("rank_%s", prefix.c_str()), 100, 0., 100.);
  // -- read in variable ranking from logfile
  vector<string> allLines; 
  char  buffer[2000];
  cout << "getRanking: open file " << Form("%s.log", fname.c_str()) << endl;
  ifstream is(Form("%s.log", fname.c_str()));
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  string::size_type m1, m2;
  string varn, vars, varp; 
  int bail(0), istart(0); 
  string after;

  if (type == "events0") {
    after = "Apply on events0, train on events1, test on events2";
  } else if (type == "events1") {
    after = "Apply on events1, train on events2, test on events0";
  } else if (type == "events2") {
    after = "Apply on events2, train on events0, test on events1";
  }

  for (unsigned int i = 0; i < allLines.size(); ++i) {
    if (string::npos != allLines[i].find(after)) {
      istart = i; 
      break;
    }
  }

  cout << "start after: " << after << " at line " << istart << endl;
  for (unsigned int i = istart; i < allLines.size(); ++i) {
    // -- method unspecific classification
    if ((string::npos != allLines[i].find(Form("--- %s", prefix.c_str())))
	&& (string::npos != allLines[i].find(": Rank : Variable "))
	) {
      bail = 0; 
      for (unsigned int j = i+2; j < i+100; ++j) {
	if (string::npos != allLines[j].find(": ---------------------------------")) {
	  bail = 1;
	  cout << "  -> breaking out " << endl;
	  break;
	}
	
	m1 = allLines[j].find(":"); 
	m2 = allLines[j].find(":", m1+1);
	varn = allLines[j].substr(m1+2, m2-m1-2); 
	m1 = m2; 
	m2 = allLines[j].find(":", m1+1);
	vars = allLines[j].substr(m1+2, m2-m1-2); 
	m1 = m2; 
	m2 = allLines[j].find(":", m1+1);
	varp = allLines[j].substr(m1+2, m2-m1-2); 
	cout << varn << "-> " << vars << " -> " << varp << endl;
	//	cout << allLines[j] << endl;
	int ibin = atoi(varn.c_str()); 
	h1->GetXaxis()->SetBinLabel(ibin, vars.c_str());
	h1->SetBinContent(ibin, atof(varp.c_str()));
      }
      if (1 == bail) break;
    }
  }

  return h1;
}    


// ----------------------------------------------------------------------
void tmva1::cleanup(string fname) {
  cout << "     cleanup of " << Form("%s-Events0.root, %s-Events1.root, %s-Events2.root", fname.c_str(), fname.c_str(), fname.c_str()) << endl;
  TFile *f = TFile::Open(Form("%s-Events0.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h0 = (TH1D*)f->Get("res_ks"); 
  if (0 == h0) {
    cout << "      hist res_ks not found " << endl;
    return;
  }

  f = TFile::Open(Form("%s-Events1.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h1 = (TH1D*)f->Get("res_ks"); 
  if (0 == h1) {
    cout << "      hist res_ks not found " << endl;
    return;
  }

  f = TFile::Open(Form("%s-Events2.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h2 = (TH1D*)f->Get("res_ks"); 
  if (0 == h2) {
    cout << "      hist res_ks not found " << endl;
    return;
  }
  
  // -- change!
  f = TFile::Open(Form("%s-combined.root", fname.c_str()));
  TH1D *H = (TH1D*)f->Get("res_ssb");
  if (0 == H) {
    cout << "      hist res_ssb not found " << endl;
    return;
  }

  double kssg0 = h0->GetBinContent(1); 
  double ksbg0 = h0->GetBinContent(2); 
  double kssg1 = h1->GetBinContent(1); 
  double ksbg1 = h1->GetBinContent(2); 
  double kssg2 = h2->GetBinContent(1); 
  double ksbg2 = h2->GetBinContent(2); 
  double ssb  = H->GetBinContent(1); 
  cout << fname << " performance: kssg0 = " << kssg0 << " ksbg0 = " << ksbg0
       << " kssg1 = " << kssg1 << " ksbg1 = " << ksbg1 
       << " kssg2 = " << kssg2 << " ksbg2 = " << ksbg2 
       << " ssb = "  << ssb << endl;
  double ssbCut = 1.3; 
  if (1 == fChannel) ssbCut = 0.9; 
  if (kssg0 < 0.05 || ksbg0 < 0.05 || kssg1 < 0.05 || ksbg1 < 0.05 || kssg2 < 0.05 || ksbg2 < 0.05 || ssb < ssbCut) {
    system(Form("/bin/rm -f %s-Events0.root", fname.c_str()));
    system(Form("/bin/rm -f %s-Events1.root", fname.c_str()));
    system(Form("/bin/rm -f %s-Events2.root", fname.c_str()));
    cout << Form("===> REMOVED %s-[Events0,Events1,Events2].root, kssg0 = %f ksbg0 = %f kssg1 = %f ksbg1 = %f kssg2 = %f ksbg2 = %f SSB = %f", 
		 fname.c_str(), kssg0, ksbg0, kssg1, ksbg1, kssg2, ksbg2, ssb) << endl;
    system(Form("/bin/rm -f weights/%s-Events0_BDT.weights.xml", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events1_BDT.weights.xml", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events2_BDT.weights.xml", fname.c_str()));
    cout << Form("     REMOVED weights/%s-[Events0,Events1,Events2]_BDT.weights.xml", fname.c_str()) << endl;
    system(Form("/bin/rm -f weights/%s-Events0_BDT.class.C", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events1_BDT.class.C", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events2_BDT.class.C", fname.c_str()));
    cout << Form("     REMOVED weights/%s-[Events0,Events1,Events2]_BDT.class.C", fname.c_str()) << endl;
  } else {
    cout << Form("===> KEEP %s.root, kssg0 = %f ksbg0 = %f kssg1 = %f ksbg1 = %f kssg0 = %f ksbg0 = %f SSB = %f", 
		 fname.c_str(), kssg0, ksbg0, kssg1, ksbg1, kssg2, ksbg2, ssb) << endl;
  }
  
}


// ----------------------------------------------------------------------
void tmva1::makeAll(int offset, string filename, int clean) {
  // createInputFile(filename); 

  if (filename == "") {
    if (2011 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-2011.root"; 
    } 
    if (2012 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-2012.root"; 
    } 
  }


  make(offset, filename, 0, clean);
  make(offset, filename, 1, clean);
  make(offset, filename, 2, clean);
  
  string oname = Form("TMVA-%d", offset); 
  cout << "-->apply(...)" << endl;
  apply(oname.c_str());
  cout << "-->analyze(...)" << endl;
  analyze(oname.c_str()); 

  cout << "-->mvas(...)" << endl;
  string sEvents = oname + "-Events0";
  mvas(sEvents.c_str()); 
  sEvents = oname + "-Events1";
  mvas(sEvents.c_str()); 
  sEvents = oname + "-Events2";
  mvas(sEvents.c_str()); 

  cout << "-->cleanup(...)" << endl;
  if (clean) cleanup(oname);

}

// ----------------------------------------------------------------------
void tmva1::make(int offset, string filename, int evt, int clean) {

  if (0 == evt)  setApply0();
  if (1 == evt)  setApply1();
  if (2 == evt)  setApply2();

  string type; 
  switch (evt) {
  case 0: 
    type = "Events0"; 
    break;
  case 1: 
    type = "Events1"; 
    break;
  case 2: 
    type = "Events2"; 
    break;
  default: 
    cout << "All hell break loose" << endl;
  }

  string oname = Form("TMVA-%d-%s", offset, type.c_str()); 
  cout << "======================================================================" << endl;
  cout << "==> tmva1(" << oname << ") " << endl;
  cout << "======================================================================" << endl;

  cout << "-->train(...) with oname = " << oname << " and filename = " << filename << endl;
  train(oname, filename);
}


// ----------------------------------------------------------------------
void tmva1::createInputFile(string filename) {
  TFile *sinput = TFile::Open(fInputFiles.sname.c_str());
  TFile *dinput = TFile::Open(fInputFiles.dname.c_str());

  // -- cuts definition: this should be identical to preselection()
  //   TCut sgcut0 = "hlt&&gmuid&&pt<70&&pt>6"; 
  //   sgcut0 += "m1pt>4.0&&m1pt<50";
  //   sgcut0 += "m2pt>4.0&&m2pt<20"; 
  //   sgcut0 += "fl3d<2&&chi2/dof<10&&pvip<0.02&&!TMath::IsNaN(pvips)&&pvips<5&&pvips>0&&maxdoca<0.02";
  //   sgcut0 += "closetrk<21"; 
  //   sgcut0 += "fls3d<100&&docatrk<0.2";
  //   sgcut0 += "alpha<0.3&&iso>0.6&&chi2/dof<10"; 

  TCut sgcut = preselection().c_str(); 
  
  //  TCut sgcut = sgcut0;  

  //   cout << "original: " << endl;
  //   cout << sgcut0 << endl;

  cout << "new: " << endl;
  cout << sgcut << endl;
  
  TCut masscut = "m>4.9&&m<5.9"; 
  TCut massbg  = "!(5.2<m&&m<5.45)";
  
  cout << "==> signal input file:     " << sinput->GetName() << std::endl;
  cout << "==> background input file: " << dinput->GetName() << std::endl;
  
  TTree *signal      = (TTree*)sinput->Get("candAnaMuMu/events");
  TTree *cbackground = (TTree*)dinput->Get("candAnaMuMu/events");
  
  TFile *outFile = TFile::Open(filename.c_str(),"RECREATE");

  // -- channel selection/definition
  string chanDef[] = {"TMath::Abs(m1eta) < 1.4 && TMath::Abs(m2eta) < 1.4", 
		      "(TMath::Abs(m1eta)>1.4||TMath::Abs(m2eta)>1.4)&&TMath::Abs(m1eta)<2.4&&TMath::Abs(m2eta)<2.4"
  };

  int nchan = 2; 
  string sdir, type; 
  TTree *copyTree(0);
  TCut copyCuts; 
  TCut chanCut, typeCut;
  for (int j = 0; j < 3; ++j) {
    if (0 == j) {
      type = "Events0"; 
      typeCut = "TMath::Abs(evt%3)==0";
      //      typeCut = "3*rndm%3==0";
    } else if (1 == j) {
      type = "Events1";
      typeCut = "TMath::Abs(evt%3)==1";
      //      typeCut = "3*rndm%3==1";
    } else if (2 == j) {
      type = "Events2";
      typeCut = "TMath::Abs(evt%3)==2";
      //      typeCut = "3*rndm%3==2";
    }

    for (int i = 0; i < nchan; ++i) {
      // -- signal
      sdir = Form("signalChan%d%s", i, type.c_str()); 
      chanCut = chanDef[i].c_str();
      outFile->mkdir(sdir.c_str());
      outFile->cd(sdir.c_str());
      copyCuts = sgcut + chanCut + typeCut;
      cout << "sg copyCuts: " << copyCuts << endl;
      copyTree = signal->CopyTree(copyCuts);
      cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;

      // -- background
      sdir = Form("sidebandChan%d%s", i, type.c_str()); 
      chanCut = chanDef[i].c_str();
      outFile->mkdir(sdir.c_str());
      outFile->cd(sdir.c_str());
      copyCuts = sgcut + massbg + masscut + chanCut + typeCut;
      cout << "bg copyCuts: " << copyCuts << endl;
      copyTree = cbackground->CopyTree(copyCuts);
      cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;
    }
  
  }

  outFile->Write();
  outFile->Close();   
  
  sinput->Close();
  dinput->Close();
}



// ----------------------------------------------------------------------
void tmva1::calcBDT() {
  fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;

  frd.pt = ftd.pt; 
  frd.eta = ftd.eta; 
  frd.m1eta = ftd.m1eta; 
  frd.m2eta = ftd.m2eta; 
  frd.m1pt = ftd.m1pt; 
  frd.m2pt = ftd.m2pt;
  frd.fls3d = ftd.fls3d; 
  frd.alpha = ftd.alpha; 
  frd.maxdoca = ftd.maxdoca;
  frd.pvip = ftd.pvip; 
  frd.pvips = ftd.pvips; 
  frd.iso = ftd.iso; 
  frd.docatrk = ftd.docatrk; 
  frd.chi2dof = ftd.chi2/ftd.dof; 
  frd.closetrk = ftd.closetrk; 

  frd.m1iso = ftd.m1iso; 
  frd.m2iso = ftd.m2iso; 
  frd.pvdchi2 = ftd.pvdchi2; 
  
  frd.m  = ftd.m; 
  int ichan = 0; 
  if (TMath::Abs(ftd.evt%3) == 1) ichan = 1; 
  if (TMath::Abs(ftd.evt%3) == 2) ichan = 2; 
  fBDT   = fReader[ichan]->EvaluateMVA("BDT"); 
  fBDT0  = fReader[0]->EvaluateMVA("BDT"); 
  fBDT1  = fReader[1]->EvaluateMVA("BDT"); 
  fBDT2  = fReader[2]->EvaluateMVA("BDT"); 
  //  cout << "calcBDT: evt = " << ichan << " BDT = " <<  fBDT << endl; 
}



// ----------------------------------------------------------------------
TMVA::Reader* setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  
  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  //  TString methodName = TString(fBdt) + TString(" method");
  //  TString weightfile = dir + fBdt + "_" + methodNameprefix + TString(".weights.xml");
  TString weightfile = xmlFile;

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  cout << "setupReader, open file " << weightfile << endl;
  ifstream is(weightfile); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1); 
  string::size_type m1, m2;
  string stype; 
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	//	cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
	if (stype == "m1pt") {
	  cout << "  adding m1pt" << endl;
	  reader->AddVariable( "m1pt", &rd.m1pt);
	}
	if (stype == "m2pt") {
	  cout << "  adding m2pt" << endl;
	  reader->AddVariable( "m2pt", &rd.m2pt);
	}
	if (stype == "m1eta") {
	  cout << "  adding m1eta" << endl;
	  reader->AddVariable( "m1eta", &rd.m1eta);
	}
	if (stype == "m2eta") {
	  reader->AddVariable( "m2eta", &rd.m2eta);
	  cout << "  adding m2eta" << endl;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable( "pt", &rd.pt);
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable( "eta", &rd.eta);
	}
	if (stype == "fls3d") {
	  cout << "  adding fls3d" << endl;
	  reader->AddVariable( "fls3d", &rd.fls3d);
	}
	if (stype == "alpha") {
	  cout << "  adding alpha" << endl;
	  reader->AddVariable( "alpha", &rd.alpha);
	}
	if (stype == "maxdoca") {
	  cout << "  adding maxdoca" << endl;
	  reader->AddVariable( "maxdoca", &rd.maxdoca);
	}
	if (stype == "pvip") {
	  cout << "  adding pvip" << endl;
	  reader->AddVariable( "pvip", &rd.pvip);
	}
	if (stype == "pvips") {
	  cout << "  adding pvips" << endl;
	  reader->AddVariable( "pvips", &rd.pvips);
	}
	if (stype == "iso") {
	  cout << "  adding iso" << endl;
	  reader->AddVariable( "iso", &rd.iso);
	}
	if (stype == "docatrk") {
	  cout << "  adding docatrk" << endl;
	  reader->AddVariable( "docatrk", &rd.docatrk);
	}
	if (stype == "closetrk") {
	  cout << "  adding closetrk" << endl;
	  reader->AddVariable( "closetrk", &rd.closetrk);
	}
	if (stype == "chi2/dof") {
	  cout << "  adding chi2/dof" << endl;
	  reader->AddVariable( "chi2dof := chi2/dof", &rd.chi2dof);
	}
	if (stype == "m1iso") {
	  cout << "  adding m1iso" << endl;
	  reader->AddVariable( "m1iso", &rd.m1iso);
	}
	if (stype == "m2iso") {
	  cout << "  adding m2iso" << endl;
	  reader->AddVariable( "m2iso", &rd.m2iso);
	}
	if (stype == "pvdchi2") {
	  cout << "  adding pvdchi2" << endl;
	  reader->AddVariable( "pvdchi2", &rd.pvdchi2);
	}
      }
      break;
    }
  }
  
  nvars = -1; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
	if (stype == "m") {
	  cout << "  adding m as spectator" << endl;
	  reader->AddSpectator( "m", &rd.m);  
	}
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", weightfile); 
  return reader; 
}


// ----------------------------------------------------------------------
void tmva1::writeOut(TFile *f, TH1 *h) {
  TDirectory *pD = gDirectory; 
  f->cd();
  h->SetDirectory(f); 
  h->Write();
  pD->cd();
}


// ----------------------------------------------------------------------
void tmva1::redrawStats(double x, double y, const char *newname, int color) {

  TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
  st->SetName(newname);
  st->SetX1NDC(x); 
  st->SetY1NDC(y); 
  st->SetTextColor(color); 
  //  st->Draw("sames");
}

// ----------------------------------------------------------------------
void tmva1::newLegend(double x1, double y1, double x2, double y2, string title) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}

// ----------------------------------------------------------------------
void setTitles(TH1 *h, const char *sx, const char *sy, float size, 
	       float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}


// ----------------------------------------------------------------------
void tmva1::toyRuns(string ifilename, int nruns) {
  string oname; 
  int offset(100);
  int seed(0); 
  int nsg(20000), nbg(25000);
  vector<int> vcolor; 
  vcolor.push_back(kBlue); 
  vcolor.push_back(kCyan); 
  vcolor.push_back(kRed); 
  vcolor.push_back(kMagenta); 
  vcolor.push_back(kBlack); 
  if (1) {
    for (int i = 0; i < nruns; ++i) {
      seed = offset + i;
      oname = Form("/scratch/ursl/tmva-toy-%d.root", seed); 
      createToyData(ifilename, oname, seed, nsg, nbg); 
      trainOnToyData(oname, Form("toy-%d", seed)); 
    }
  } else {
    cout << "UNCOMMENT CODE IF YOU WANT TO RUN NEW SAMPLES" << endl;
  }

  TCanvas *c0 = getC0();
  c0->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double ymax = 0.07; 
  TH2F* frame = new TH2F("frame", "BDT output distributions", 100, -1., 1., 100, 0., ymax);
  frame->GetXaxis()->SetTitle(" b");
  frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
  frame->Draw();  

  int color(0);
  for (int i = 0; i < nruns; ++i) {
    seed = offset + i;
    oname = Form("toy-%d.root", seed);
    TFile *file = TFile::Open(oname.c_str()); 
    
    color = (i < vcolor.size()? vcolor[i] :kBlack); 
    TH1F *hTrainBg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_B"); hTrainBg->SetLineColor(color); hTrainBg->SetMarkerColor(color);
    TH1F *hTrainSg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_S"); hTrainSg->SetLineColor(color); hTrainSg->SetMarkerColor(color);
    
    hTrainBg->Scale(1./hTrainBg->GetSumOfWeights());
    hTrainSg->Scale(1./hTrainSg->GetSumOfWeights());
    
    hTrainBg->DrawCopy("histsame");
    hTrainSg->DrawCopy("histsame");
    file->Close();
  }

  c0->SaveAs(Form("toys-training-%d.pdf", nruns));

  c0->Clear();
  frame->Draw();  
  for (int i = 0; i < nruns; ++i) {
    seed = offset + i;
    oname = Form("toy-%d.root", seed);
    TFile *file = TFile::Open(oname.c_str()); 

    color = (i < vcolor.size()? vcolor[i] :kBlack); 
    TH1F *hTestBg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_B"); hTestBg->SetLineColor(color); hTestBg->SetMarkerColor(color);
    TH1F *hTestSg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_S"); hTestSg->SetLineColor(color); hTestSg->SetMarkerColor(color);
    
    hTestBg->Scale(1./hTestBg->GetSumOfWeights());
    hTestSg->Scale(1./hTestSg->GetSumOfWeights());
    
    hTestBg->DrawCopy("histsame");
    hTestSg->DrawCopy("histsame");

    file->Close();
  }

  c0->SaveAs(Form("toys-testing-%d.pdf", nruns));

}

// ----------------------------------------------------------------------
void tmva1::createToyData(string ifilename, string ofilename, int seed, int nsg, int nbg) {
  delete gRandom; 
  gRandom = new TRandom3(seed); 

  cout << "==> CREATE TOY DATA for seed = " << seed << endl;

  int channel = 0; 

  vector<string> vNames;
  vector<double> vMin, vMax; 
  vector<int> vNbins; 
  vNames.push_back("m1pt"); vMin.push_back(0.); vMax.push_back(40.); vNbins.push_back(100); 
  vNames.push_back("m2pt"); vMin.push_back(0.); vMax.push_back(20.); vNbins.push_back(100); 
  vNames.push_back("m1eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("m2eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("pt");    vMin.push_back(0.); vMax.push_back(40.); vNbins.push_back(100); 
  vNames.push_back("eta");   vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("fls3d"); vMin.push_back(0.); vMax.push_back(100.); vNbins.push_back(120); 
  vNames.push_back("alpha"); vMin.push_back(0.); vMax.push_back(0.3); vNbins.push_back(100); 
  vNames.push_back("maxdoca"); vMin.push_back(0.); vMax.push_back(0.03); vNbins.push_back(100); 
  vNames.push_back("pvip"); vMin.push_back(0.); vMax.push_back(0.05); vNbins.push_back(100); 
  vNames.push_back("pvips"); vMin.push_back(0.); vMax.push_back(5); vNbins.push_back(100); 
  vNames.push_back("iso"); vMin.push_back(0.6); vMax.push_back(1.01); vNbins.push_back(41); 
  vNames.push_back("closetrk"); vMin.push_back(0.); vMax.push_back(21); vNbins.push_back(21); 
  vNames.push_back("docatrk"); vMin.push_back(0.); vMax.push_back(0.1); vNbins.push_back(100); 
  vNames.push_back("chi2dof"); vMin.push_back(0.); vMax.push_back(5); vNbins.push_back(100); 
  vNames.push_back("m1iso"); vMin.push_back(0.); vMax.push_back(1.01); vNbins.push_back(101); 
  vNames.push_back("m2iso"); vMin.push_back(0.); vMax.push_back(1.01); vNbins.push_back(101); 
  vNames.push_back("pvdchi2"); vMin.push_back(0.); vMax.push_back(10.0); vNbins.push_back(100); 

  TFile *inFile = TFile::Open(ifilename.c_str());
  TTree *tsg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", channel));
  TTree *tbg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", channel));

  TH1D *hs, *hb; 
  string hsName, hbName; 

  TCanvas *c0 = getC0();
  c0->Clear();
  c0->Divide(4,4);
  for (int i = 0; i < vNames.size(); ++i) {
    hsName = Form("hs_%s", vNames[i].c_str());
    hs  = new TH1D(hsName.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);

    hbName = Form("hb_%s", vNames[i].c_str());
    hb  = new TH1D(hbName.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);
    
    string var = vNames[i]; 
    if (var == string("chi2dof")) var = "chi2/dof"; 
    cout << "--> " << var << endl;
    tsg->Draw(Form("%s>>%s", var.c_str(), hsName.c_str()), "", "goff");
    tbg->Draw(Form("%s>>%s", var.c_str(), hbName.c_str()), "", "goff");
    c0->cd(i+1); 
    hs->Draw("hist"); 
    hb->Draw("esame");
  }

  c0->cd(vNames.size()+1);

  TFile *outfile = TFile::Open(ofilename.c_str(), "RECREATE");  
  struct readerData rd; 

  outfile->mkdir("signalChan0Events0");
  outfile->cd("signalChan0Events0");
  TTree *osg0 = createTree(rd); 
  for (int i = 0 ; i < nsg; ++i) {
    rd.m = 5.37;
    rd.pt = ((TH1D*)inFile->Get("hs_pt"))->GetRandom();
    rd.eta = ((TH1D*)inFile->Get("hs_eta"))->GetRandom();
    rd.m1eta = ((TH1D*)inFile->Get("hs_m1eta"))->GetRandom();
    rd.m2eta = ((TH1D*)inFile->Get("hs_m2eta"))->GetRandom();
    rd.m1pt = ((TH1D*)inFile->Get("hs_m1pt"))->GetRandom();
    rd.m2pt = ((TH1D*)inFile->Get("hs_m2pt"))->GetRandom();

    rd.fls3d = ((TH1D*)inFile->Get("hs_fls3d"))->GetRandom();
    rd.alpha = ((TH1D*)inFile->Get("hs_alpha"))->GetRandom();
    rd.maxdoca = ((TH1D*)inFile->Get("hs_maxdoca"))->GetRandom();
    rd.pvip = ((TH1D*)inFile->Get("hs_pvip"))->GetRandom();
    rd.pvips = ((TH1D*)inFile->Get("hs_pvips"))->GetRandom();
    rd.iso = ((TH1D*)inFile->Get("hs_iso"))->GetRandom();
    rd.docatrk = ((TH1D*)inFile->Get("hs_docatrk"))->GetRandom();
    rd.chi2dof = ((TH1D*)inFile->Get("hs_chi2dof"))->GetRandom();
    rd.closetrk = ((TH1D*)inFile->Get("hs_closetrk"))->GetRandom();

    rd.m1iso = ((TH1D*)inFile->Get("hs_m1iso"))->GetRandom();
    rd.m2iso = ((TH1D*)inFile->Get("hs_m2iso"))->GetRandom();
    rd.pvdchi2 = ((TH1D*)inFile->Get("hs_pvdchi2"))->GetRandom();

    osg0->Fill();
  }

  outfile->mkdir("signalChan0Events1");
  outfile->cd("signalChan0Events1");
  TTree *osg1 = createTree(rd); 
  for (int i = 0 ; i < nsg; ++i) {
    rd.m = 5.37;
    rd.pt = ((TH1D*)inFile->Get("hs_pt"))->GetRandom();
    rd.eta = ((TH1D*)inFile->Get("hs_eta"))->GetRandom();
    rd.m1eta = ((TH1D*)inFile->Get("hs_m1eta"))->GetRandom();
    rd.m2eta = ((TH1D*)inFile->Get("hs_m2eta"))->GetRandom();
    rd.m1pt = ((TH1D*)inFile->Get("hs_m1pt"))->GetRandom();
    rd.m2pt = ((TH1D*)inFile->Get("hs_m2pt"))->GetRandom();

    rd.fls3d = ((TH1D*)inFile->Get("hs_fls3d"))->GetRandom();
    rd.alpha = ((TH1D*)inFile->Get("hs_alpha"))->GetRandom();
    rd.maxdoca = ((TH1D*)inFile->Get("hs_maxdoca"))->GetRandom();
    rd.pvip = ((TH1D*)inFile->Get("hs_pvip"))->GetRandom();
    rd.pvips = ((TH1D*)inFile->Get("hs_pvips"))->GetRandom();
    rd.iso = ((TH1D*)inFile->Get("hs_iso"))->GetRandom();
    rd.docatrk = ((TH1D*)inFile->Get("hs_docatrk"))->GetRandom();
    rd.chi2dof = ((TH1D*)inFile->Get("hs_chi2dof"))->GetRandom();
    rd.closetrk = ((TH1D*)inFile->Get("hs_closetrk"))->GetRandom();

    rd.m1iso = ((TH1D*)inFile->Get("hs_m1iso"))->GetRandom();
    rd.m2iso = ((TH1D*)inFile->Get("hs_m2iso"))->GetRandom();
    rd.pvdchi2 = ((TH1D*)inFile->Get("hs_pvdchi2"))->GetRandom();

    osg1->Fill();
  }

  outfile->mkdir("sidebandChan0Events0");
  outfile->cd("sidebandChan0Events0");
  TTree *obg0 = createTree(rd); 
  for (int i = 0 ; i < nbg; ++i) {
    rd.m = 5.37;
    rd.pt = ((TH1D*)inFile->Get("hb_pt"))->GetRandom();
    rd.eta = ((TH1D*)inFile->Get("hb_eta"))->GetRandom();
    rd.m1eta = ((TH1D*)inFile->Get("hb_m1eta"))->GetRandom();
    rd.m2eta = ((TH1D*)inFile->Get("hb_m2eta"))->GetRandom();
    rd.m1pt = ((TH1D*)inFile->Get("hb_m1pt"))->GetRandom();
    rd.m2pt = ((TH1D*)inFile->Get("hb_m2pt"))->GetRandom();

    rd.fls3d = ((TH1D*)inFile->Get("hb_fls3d"))->GetRandom();
    rd.alpha = ((TH1D*)inFile->Get("hb_alpha"))->GetRandom();
    rd.maxdoca = ((TH1D*)inFile->Get("hb_maxdoca"))->GetRandom();
    rd.pvip = ((TH1D*)inFile->Get("hb_pvip"))->GetRandom();
    rd.pvips = ((TH1D*)inFile->Get("hb_pvips"))->GetRandom();
    rd.iso = ((TH1D*)inFile->Get("hb_iso"))->GetRandom();
    rd.docatrk = ((TH1D*)inFile->Get("hb_docatrk"))->GetRandom();
    rd.chi2dof = ((TH1D*)inFile->Get("hb_chi2dof"))->GetRandom();
    rd.closetrk = ((TH1D*)inFile->Get("hb_closetrk"))->GetRandom();

    rd.m1iso = ((TH1D*)inFile->Get("hb_m1iso"))->GetRandom();
    rd.m2iso = ((TH1D*)inFile->Get("hb_m2iso"))->GetRandom();
    rd.pvdchi2 = ((TH1D*)inFile->Get("hb_pvdchi2"))->GetRandom();

    obg0->Fill();
  }

  outfile->mkdir("sidebandChan0Events1");
  outfile->cd("sidebandChan0Events1");
  TTree *obg1 = createTree(rd); 
  for (int i = 0 ; i < nbg; ++i) {
    rd.m = 5.37;
    rd.pt = ((TH1D*)inFile->Get("hb_pt"))->GetRandom();
    rd.eta = ((TH1D*)inFile->Get("hb_eta"))->GetRandom();
    rd.m1eta = ((TH1D*)inFile->Get("hb_m1eta"))->GetRandom();
    rd.m2eta = ((TH1D*)inFile->Get("hb_m2eta"))->GetRandom();
    rd.m1pt = ((TH1D*)inFile->Get("hb_m1pt"))->GetRandom();
    rd.m2pt = ((TH1D*)inFile->Get("hb_m2pt"))->GetRandom();

    rd.fls3d = ((TH1D*)inFile->Get("hb_fls3d"))->GetRandom();
    rd.alpha = ((TH1D*)inFile->Get("hb_alpha"))->GetRandom();
    rd.maxdoca = ((TH1D*)inFile->Get("hb_maxdoca"))->GetRandom();
    rd.pvip = ((TH1D*)inFile->Get("hb_pvip"))->GetRandom();
    rd.pvips = ((TH1D*)inFile->Get("hb_pvips"))->GetRandom();
    rd.iso = ((TH1D*)inFile->Get("hb_iso"))->GetRandom();
    rd.docatrk = ((TH1D*)inFile->Get("hb_docatrk"))->GetRandom();
    rd.chi2dof = ((TH1D*)inFile->Get("hb_chi2dof"))->GetRandom();
    rd.closetrk = ((TH1D*)inFile->Get("hb_closetrk"))->GetRandom();

    rd.m1iso = ((TH1D*)inFile->Get("hb_m1iso"))->GetRandom();
    rd.m2iso = ((TH1D*)inFile->Get("hb_m2iso"))->GetRandom();
    rd.pvdchi2 = ((TH1D*)inFile->Get("hb_pvdchi2"))->GetRandom();

    obg1->Fill();
  }

  osg0->Write();
  osg1->Write();
  obg0->Write();
  obg1->Write();

  outfile->Write();
  outfile->Close();

}


// ----------------------------------------------------------------------
void tmva1::trainOnToyData(string iname, string oname) {

  // This loads the library
  TMVA::Tools::Instance();
  
  (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40; 
  (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 100; 
  
  // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName(Form("%s.root", oname.c_str()));
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  string optstring = "V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
  optstring        = "V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
  cout << "==> Factory: " << optstring << endl;
  TMVA::Factory *factory = new TMVA::Factory(Form("%s", oname.c_str()), outputFile,  optstring.c_str());
  
  factory->AddVariable("m1pt",       'F' );
  factory->AddVariable("m2pt",       'F' );
  factory->AddVariable("m1eta",      'F' );
  factory->AddVariable("m2eta",      'F' );
  factory->AddVariable("pt",         'F' );
  factory->AddVariable("eta",        'F' );
  factory->AddVariable("fls3d",      'F' );
  factory->AddVariable("alpha",      'F' );
  factory->AddVariable("maxdoca",    'F' );
  factory->AddVariable("pvip",       'F' );
  factory->AddVariable("pvips",      'F' );
  factory->AddVariable("iso",        'F' );
  factory->AddVariable("docatrk",    'F' );
  factory->AddVariable("closetrk",   'I' );
  factory->AddVariable("chi2dof",    'F' );

  factory->AddVariable("m1iso",    'F' );
  factory->AddVariable("m2iso",    'F' );
  factory->AddVariable("pvdchi2",    'F' );
  
  factory->AddSpectator("m",  "mass", "GeV", 'F' );
  
  TTree *trainSg(0), *testSg(0), *trainBg(0), *testBg(0); 
  
  int channel(0); 
  
  TFile *inFile = TFile::Open(iname.c_str());
  cout << "==============> Apply on events0, train on events1, test on events2" << endl;
  trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", channel));
  testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", channel));
  trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", channel));
  testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", channel));
  cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
  cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
  
  Double_t signalWeight      = 1.; //= LUMISCALE; // 0.000388
  Double_t rbackgroundWeight = 1.;
  Double_t cbackgroundWeight = 1.;
  Double_t tbackgroundWeight = cbackgroundWeight;
  
  factory->AddTree(trainSg,     "Signal",     signalWeight,  "", "train");
  factory->AddTree(testSg,      "Signal",     signalWeight,  "", "test");
  factory->AddTree(trainBg, "Background", cbackgroundWeight, "", "train");
  factory->AddTree(testBg,  "Background", tbackgroundWeight, "", "test");
  
  int nSgTrain = trainSg->GetEntries();
  int nSgTest  = testSg->GetEntries();
  
  int nBgTrain = trainBg->GetEntries();
  int nBgTest  = testBg->GetEntries();
  
  optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V", 
		   nSgTrain, nSgTest, nBgTrain, nBgTest); 
  cout << "==> PrepareTrainingAndTestTree: " << optstring << endl;
  factory->PrepareTrainingAndTestTree("", "", optstring.c_str());
  
  if (1) {
    optstring = Form("!H:V:NTrees=%d:nEventsMin=%d", fBdtSetup.NTrees, fBdtSetup.nEventsMin);
    optstring += Form(":BoostType=AdaBoost:AdaBoostBeta=%f:SeparationType=GiniIndex:nCuts=%d:PruneMethod=NoPruning", 
		      fBdtSetup.AdaBoostBeta, fBdtSetup.nCuts);
    
    optstring += Form(":MaxDepth=%d:NNodesMax=%d", fBdtSetup.MaxDepth, fBdtSetup.NNodesMax);
  } else {
    optstring = "";
  }
  cout << "==> BookMethod: " << optstring << endl;
  factory->BookMethod( TMVA::Types::kBDT, "BDT", optstring);
  
  cout << "==> TrainAllMethods " << endl;
  factory->TrainAllMethods();
  
  // ---- Evaluate all MVAs using the set of test events
  cout << "==> TestAllMethods " << endl;
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  cout << "==> EvaluateAllMethods " << endl;
  factory->EvaluateAllMethods();
  
  // Save the output
  outputFile->Close();
  
  TFile *file = TFile::Open(outfileName.Data()); 
  
  TH1F *hTrainBg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_B"); hTrainBg->SetLineColor(kRed); hTrainBg->SetMarkerColor(kRed);
  TH1F *hTestBg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_B"); hTestBg->SetLineColor(kRed); hTestBg->SetMarkerColor(kRed);
  
  TH1F *hTrainSg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_S"); hTrainSg->SetLineColor(kBlue); hTrainSg->SetMarkerColor(kBlue);
  TH1F *hTestSg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_S"); hTestSg->SetLineColor(kBlue); hTestSg->SetMarkerColor(kBlue);
  
  TCanvas *c0 = getC0();

  hTrainBg->Scale(1./hTrainBg->GetSumOfWeights());
  hTestBg->Scale(1./hTestBg->GetSumOfWeights());
  hTrainSg->Scale(1./hTrainSg->GetSumOfWeights());
  hTestSg->Scale(1./hTestSg->GetSumOfWeights());

  gStyle->SetOptStat(0);
  double ymax = (hTrainBg->GetMaximum() > hTrainSg->GetMaximum()?hTrainBg->GetMaximum():hTrainSg->GetMaximum());
  TH2F* frame = new TH2F("frame", "BDT output distributions", 100, -1., 1., 100, 0., 1.3*ymax);
  frame->GetXaxis()->SetTitle(" b");
  frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
  frame->Draw();  

  hTrainBg->Draw("hist");
  hTestBg->Draw("psame");

  hTrainSg->Draw("histsame");
  hTestSg->Draw("psame");
  
  c0->SaveAs(Form("%s.pdf", oname.c_str()));

  delete factory; 
   
}


// ----------------------------------------------------------------------
TTree* tmva1::createTree(struct readerData &rd) {


  TTree *tree = new TTree("events", "events");
  tree->Branch("m", &rd.m, "m/F");
  tree->Branch("pt", &rd.pt, "pt/F");
  tree->Branch("eta", &rd.eta, "eta/F");
  tree->Branch("m1eta", &rd.m1eta, "m1eta/F");
  tree->Branch("m2eta", &rd.m2eta, "m2eta/F");
  tree->Branch("m1pt", &rd.m1pt, "m1pt/F");
  tree->Branch("m2pt", &rd.m2pt, "m2pt/F");

  tree->Branch("fls3d", &rd.fls3d, "fls3d/F");
  tree->Branch("alpha", &rd.alpha, "alpha/F");
  tree->Branch("maxdoca", &rd.maxdoca, "maxdoca/F");
  tree->Branch("pvip", &rd.pvip, "pvip/F");
  tree->Branch("pvips", &rd.pvips, "pvips/F");
  tree->Branch("iso", &rd.iso, "iso/F");
  tree->Branch("docatrk", &rd.docatrk, "docatrk/F");
  tree->Branch("chi2dof", &rd.chi2dof, "chi2dof/F");
  tree->Branch("closetrk", &rd.closetrk, "closetrk/F");

  tree->Branch("m1iso", &rd.m1iso, "m1iso/F");
  tree->Branch("m2iso", &rd.m2iso, "m2iso/F");
  tree->Branch("pvdchi2", &rd.pvdchi2, "pvdchi2/F");

  return tree; 
}
