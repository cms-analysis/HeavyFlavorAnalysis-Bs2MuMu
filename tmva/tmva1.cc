#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "tmvaglob.hh"
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
#include "TMinuit.h"
#include "TVirtualFitter.h"

#include "tmva1.hh"
#include "../macros/initFunc.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)#
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

//2011
//#define LUMISCALE 2.01e-4

//2012  12/7405 = 0.00162
//#define LUMISCALE 0.00162

ClassImp(tmva1)

using namespace std; 

// ----------------------------------------------------------------------
// -- 
// -- USAGE: a.makeAll(0, 1); > TMVA-0.log
// --
// ----------------------------------------------------------------------

tmva1::tmva1(int year, string vars) {

  cout << "tmva1 hello: setup for year = " << year << " with variables: " << vars << endl;

  fVariables = vars; 
  fBDTParameters = ""; 

  legg = 0;
  legge = 0; 
  tl = new TLatex(); 
  tl->SetTextFont(42);
  tl->SetTextSize(0.03);
  tl->SetNDC(kTRUE); 

  fYear = year; 

  if (year == 2011) {
    fLumiScale = 3.1e-4; // 4.9/16000
    fInputFiles.sname = "/scratch/ursl/bdt/v16-2011-mix-Bs2MuMu.root"; 
    fInputFiles.dname = "/scratch/ursl/bdt/v16-2011-data-bmmLoose-2.root";
  } else {
    fLumiScale = 2.8e-4; // 20/714000
    fInputFiles.sname = "/scratch/ursl/bdt/v16-2012-cms-BsToMuMu-small.root"; 
    fInputFiles.dname = "/scratch/ursl/bdt/v16-2012-data-bmmLoose-2.root";
  }

//   // -- BDT setup 108/109
//   fBdtSetup.NTrees = 800;
//   fBdtSetup.nEventsMin = 50; 
//   fBdtSetup.MaxDepth = 2;
//   //  fBdtSetup.MaxDepth = 3;
//   fBdtSetup.nCuts = 20; 
//   fBdtSetup.AdaBoostBeta = 1.0; 
//   fBdtSetup.NNodesMax = 5;
//   //  fBdtSetup.NNodesMax = 20;

  // -- TMVA default
  fBdtSetup.NTrees = 200;
  fBdtSetup.nEventsMin = 500; 
  fBdtSetup.MaxDepth = 3;
  fBdtSetup.nCuts = 20; 
  fBdtSetup.AdaBoostBeta = 1.0; 
  fBdtSetup.NNodesMax = 100000;

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
  t->SetBranchAddress("run", &b.run);
  t->SetBranchAddress("json", &b.json);
  t->SetBranchAddress("pvw8", &b.pvw8);
  t->SetBranchAddress("gmuid", &b.gmuid);
  t->SetBranchAddress("gmutmid", &b.gmutmid);
  t->SetBranchAddress("gtqual", &b.gtqual);
  t->SetBranchAddress("hlt", &b.hlt);
  t->SetBranchAddress("hltm", &b.hltm);
  t->SetBranchAddress("hltm2", &b.hltm2);
  t->SetBranchAddress("m1pt", &b.m1pt);
  t->SetBranchAddress("m1eta", &b.m1eta);
  t->SetBranchAddress("m1q", &b.m1q);
  t->SetBranchAddress("m2pt", &b.m2pt);
  t->SetBranchAddress("m2eta", &b.m2eta);
  t->SetBranchAddress("m2q", &b.m2q);
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
  t->SetBranchAddress("chi2dof", &b.chi2dof);
  t->SetBranchAddress("closetrk", &b.closetrk);
  t->SetBranchAddress("m", &b.m);

  t->SetBranchAddress("m1iso",&b.m1iso);
  t->SetBranchAddress("m2iso",&b.m2iso);
  t->SetBranchAddress("closetrks1", &b.closetrks1);
  t->SetBranchAddress("closetrks2", &b.closetrks2);
  t->SetBranchAddress("closetrks3", &b.closetrks3);
  t->SetBranchAddress("pvdchi2",&b.pvdchi2);
  t->SetBranchAddress("othervtx",&b.othervtx);
  t->SetBranchAddress("m1xpdist",&b.m1xpdist);
  t->SetBranchAddress("m2xpdist",&b.m2xpdist);

  t->SetBranchAddress("pvlips2", &b.pvlips2);
  t->SetBranchAddress("pvlip2", &b.pvlip2);
}


// ----------------------------------------------------------------------
void tmva1::train(string oname, string filename, int nsg, int nbg) {
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

   // -- parse string with all variables into a vector
   vector<string> vVar; 
   string allvar = fVariables; 
   string::size_type m0(allvar.size());
   cout << "--> allvar = " << allvar << endl;
   string var; 
   while (string::npos != m0) {
     m0 = allvar.find(":");
     var = allvar.substr(0, m0); 
     vVar.push_back(var); 
     allvar = allvar.substr(m0+1); 
     m0 = allvar.find(":");
   }
   var = allvar; 
   vVar.push_back(var); 

   for (unsigned int i = 0; i < vVar.size(); ++i) {
     cout << "  " << vVar[i] << endl;
     //      if (string::npos != vVar[i].find("closetrk")) {
     //        factory->AddVariable(vVar[i].c_str(), 'I');        
     //      } else {
     factory->AddVariable(vVar[i].c_str(), 'F');        
     //      }
   }

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
     cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
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
     cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
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
     cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
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

   if (nsg > -1) {
     nSgTrain = nsg; 
     nSgTest = nsg; 
   } 

   if (nbg > -1) {
     nBgTrain = nbg; 
     nBgTest = nbg; 
   }

   int seed = static_cast<int>(100*gRandom->Rndm()); 

   // optstring=Form("nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:V"); 
   // optstring=Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=random:SplitSeed=%d:NormMode=None:V", 
   // nSgTrain, nSgTest, nBgTrain, nBgTest, seed); 

   optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V", 
		    nSgTrain, nSgTest, nBgTrain, nBgTest); 
   cout << "==> PrepareTrainingAndTestTree: " << optstring << endl;
   factory->PrepareTrainingAndTestTree("", "", optstring.c_str());
   
   if (0) {
     optstring = Form("!H:V:NTrees=%d", fBdtSetup.NTrees);
     optstring += Form(":nCuts=%d:PruneMethod=NoPruning", fBdtSetup.nCuts);
     optstring += Form(":PruneMethod=NoPruning");
     optstring += Form(":BoostType=AdaBoost:AdaBoostBeta=%f:SeparationType=GiniIndex", fBdtSetup.AdaBoostBeta);
     optstring += Form(":MaxDepth=%d", fBdtSetup.MaxDepth);
     optstring += Form(":NNodesMax=%d", fBdtSetup.NNodesMax);
     //     optstring += Form(":nEventsMin=%d", fBdtSetup.nEventsMin);
   } else {
     optstring = "!H:V" + fBDTParameters;
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
   gROOT->Clear();  gROOT->DeleteAll();

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
  tt->Branch("me",      &ftd.me,     "me/D");
  tt->Branch("weight",  &w8,        "weight/D");
  tt->Branch("bdt",     &fBDT,      "bdt/D");
  tt->Branch("bdt0",    &fBDT0,     "bdt0/D");
  tt->Branch("bdt1",    &fBDT1,     "bdt1/D");
  tt->Branch("bdt2",    &fBDT2,     "bdt2/D");
  tt->Branch("classID", &classID,   "classID/I");
  tt->Branch("hlt",     &ftd.hlt,   "hlt/O");
  tt->Branch("hltm",    &ftd.hltm,  "hltm/O");
  tt->Branch("hltm2",   &ftd.hltm2,  "hltm2/O");
  tt->Branch("gmuid",   &ftd.gmuid, "gmuid/O");
  tt->Branch("evt",     &ftd.evt,   "evt/I");

  // -- add the usual variables
  tt->Branch("m1pt",    &ftd.m1pt,     "m1pt/D");
  tt->Branch("m2pt",    &ftd.m2pt,     "m2pt/D");
  tt->Branch("m1eta",   &ftd.m1eta,    "m1eta/D");
  tt->Branch("m2eta",   &ftd.m2eta,    "m2etat/D");
  tt->Branch("pt",      &ftd.pt,       "pt/D");
  tt->Branch("eta",     &ftd.eta,      "eta/D");

  tt->Branch("chi2dof", &ftd.chi2dof,  "chi2dof/D");
  tt->Branch("maxdoca", &ftd.maxdoca,  "maxdoca/D");
  tt->Branch("fls3d",   &ftd.fls3d,  "fls3d/D");
  tt->Branch("fl3d",    &ftd.fl3d,  "fl3d/D");
  tt->Branch("flsxy",   &ftd.flsxy,  "flsxy/D");
  tt->Branch("alpha",   &ftd.alpha,  "alpha/D");
  tt->Branch("pvip",    &ftd.pvip,  "pvip/D");
  tt->Branch("pvips",   &ftd.pvips,  "pvips/D");
  tt->Branch("pvlip",   &ftd.pvlip,  "pvlip/D");
  tt->Branch("pvlips",  &ftd.pvlips,  "pvlips/D");

  tt->Branch("iso",     &ftd.iso,  "iso/D");
  tt->Branch("m1iso",   &ftd.m1iso,  "m1iso/D");
  tt->Branch("m2iso",   &ftd.m2iso,  "m2iso/D");
  tt->Branch("docatrk", &ftd.docatrk,  "docatrk/D");
  tt->Branch("closetrk",&ftd.closetrk,  "closetrk/D");


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
  bool otherCuts(false);
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing data event: " << ievt << std::endl;
    t->GetEntry(ievt);
    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8; 
    otherCuts = (ftd.m > 4.9 && ftd.m < 5.9) && (ftd.m1q*ftd.m2q < 0) && (ftd.pvw8 > 0.7) && (ftd.gtqual);
    
    if (ftd.json && otherCuts && ftd.gmuid && ftd.hlt && ftd.hltm2 && preselection(ftd, fChannel)) {
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
  w8 = fLumiScale;
  cout << "--- Processing signal: " << t->GetEntries() << " events and weight " << w8 << endl;
  nEvent = t->GetEntries();
  lostEvents = 0; 
  totalEvents = 0; 
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing signal event: " << ievt << std::endl;
    t->GetEntry(ievt);

    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8; 
    otherCuts = (ftd.m > 4.9 && ftd.m < 5.9) && (ftd.m1q*ftd.m2q < 0) && (ftd.pvw8 > 0.7) && (ftd.gtqual);
    if (otherCuts && ftd.gmuid && ftd.hlt && ftd.hltm2 && preselection(ftd, fChannel)) {
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
  TH1D *hs = (TH1D*)f->Get("bdtTreeSignal"); 
  double sLostEvents  = hs->GetBinContent(10); 

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

  TH1F *trainSgBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT0->SetLineColor(kBlack);
  TH1F *testSgBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT0->SetLineColor(kBlack);

  TH1F *trainSgBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT1->SetLineColor(kRed);
  TH1F *testSgBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT1->SetLineColor(kRed);

  TH1F *trainSgBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT2->SetLineColor(kBlue);
  TH1F *testSgBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT2->SetLineColor(kBlue);

  TH1D *bBDT = new TH1D("bBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *cBDT = new TH1D("cBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *h2 = new TH1D("res_ssb", "max(S/Sqrt(S+B))", 10, 0., 10.); h2->Sumw2();  h2->SetDirectory(f); 

  // -- new rebinned versions of the BDT output plots. These ARE identical to the TMVA output IF the binning IS the same!
  const int BINS(100); 
  TH1D *ap0bgBDT = new TH1D("ap0bgBDT", "",    BINS, -1., 1.);
  TH1D *tr0bgBDT = new TH1D("tr0bgBDT", "",    BINS, -1., 1.);
  TH1D *te0bgBDT = new TH1D("te0bgBDT", "",    BINS, -1., 1.);
  ap0bgBDT->SetLineColor(kMagenta);   te0bgBDT->SetLineColor(kMagenta);   tr0bgBDT->SetLineColor(kMagenta);
  ap0bgBDT->SetMarkerColor(kMagenta); te0bgBDT->SetMarkerColor(kMagenta); tr0bgBDT->SetMarkerColor(kMagenta); 
  
  TH1D *ap1bgBDT = new TH1D("ap1bgBDT", "",    BINS, -1., 1.);
  TH1D *tr1bgBDT = new TH1D("tr1bgBDT", "",    BINS, -1., 1.);
  TH1D *te1bgBDT = new TH1D("te1bgBDT", "",    BINS, -1., 1.);
  ap1bgBDT->SetLineColor(kRed);   te1bgBDT->SetLineColor(kRed);   tr1bgBDT->SetLineColor(kRed);
  ap1bgBDT->SetMarkerColor(kRed); te1bgBDT->SetMarkerColor(kRed); tr1bgBDT->SetMarkerColor(kRed);

  TH1D *ap2bgBDT = new TH1D("ap2bgBDT", "",    BINS, -1., 1.);
  TH1D *tr2bgBDT = new TH1D("tr2bgBDT", "",    BINS, -1., 1.);
  TH1D *te2bgBDT = new TH1D("te2bgBDT", "",    BINS, -1., 1.);
  ap2bgBDT->SetLineColor(kBlue);   te2bgBDT->SetLineColor(kBlue);   tr2bgBDT->SetLineColor(kBlue);
  ap2bgBDT->SetMarkerColor(kBlue); te2bgBDT->SetMarkerColor(kBlue); tr2bgBDT->SetMarkerColor(kBlue);

  TH1D *ap3bgBDT = new TH1D("ap3bgBDT", "",    BINS, -1., 1.);
  ap3bgBDT->SetLineColor(kBlack);   
  ap3bgBDT->SetMarkerColor(kBlack); 

  TH1D *ap0sgBDT = new TH1D("ap0sgBDT", "",    BINS, -1., 1.);
  TH1D *tr0sgBDT = new TH1D("tr0sgBDT", "",    BINS, -1., 1.);
  TH1D *te0sgBDT = new TH1D("te0sgBDT", "",    BINS, -1., 1.);
  ap0sgBDT->SetLineColor(kMagenta);   te0sgBDT->SetLineColor(kMagenta);   tr0sgBDT->SetLineColor(kMagenta);
  ap0sgBDT->SetMarkerColor(kMagenta); te0sgBDT->SetMarkerColor(kMagenta); tr0sgBDT->SetMarkerColor(kMagenta);
  
  TH1D *ap1sgBDT = new TH1D("ap1sgBDT", "",    BINS, -1., 1.);
  TH1D *tr1sgBDT = new TH1D("tr1sgBDT", "",    BINS, -1., 1.);
  TH1D *te1sgBDT = new TH1D("te1sgBDT", "",    BINS, -1., 1.);
  ap1sgBDT->SetLineColor(kRed);   te1sgBDT->SetLineColor(kRed);   tr1sgBDT->SetLineColor(kRed);
  ap1sgBDT->SetMarkerColor(kRed); te1sgBDT->SetMarkerColor(kRed); tr1sgBDT->SetMarkerColor(kRed);

  TH1D *ap2sgBDT = new TH1D("ap2sgBDT", "",    BINS, -1., 1.);
  TH1D *tr2sgBDT = new TH1D("tr2sgBDT", "",    BINS, -1., 1.);
  TH1D *te2sgBDT = new TH1D("te2sgBDT", "",    BINS, -1., 1.);
  ap2sgBDT->SetLineColor(kBlue);   te2sgBDT->SetLineColor(kBlue);   tr2sgBDT->SetLineColor(kBlue);
  ap2sgBDT->SetMarkerColor(kBlue); te2sgBDT->SetMarkerColor(kBlue); tr2sgBDT->SetMarkerColor(kBlue);

  TH1D *ap3sgBDT = new TH1D("ap3sgBDT", "",    BINS, -1., 1.);
  ap3sgBDT->SetLineColor(kBlack);   
  ap3sgBDT->SetMarkerColor(kBlack); 
  
  TTree *t = (TTree*)f->Get("bdtTree");
  cout << "bdtTree with entries = " << t->GetEntries() << endl;
  double bdt, m, w8; 
  double bdt0, bdt1, bdt2; 
  int classID, evt;
  bool gmuid, hlt, hltm;
  t->SetBranchAddress("bdt", &bdt);
  t->SetBranchAddress("bdt0", &bdt0);
  t->SetBranchAddress("bdt1", &bdt1);
  t->SetBranchAddress("bdt2", &bdt2);
  t->SetBranchAddress("classID", &classID);
  t->SetBranchAddress("m", &m);
  t->SetBranchAddress("weight", &w8);
  t->SetBranchAddress("hlt", &hlt);
  t->SetBranchAddress("hltm", &hltm);
  t->SetBranchAddress("gmuid", &gmuid);
  t->SetBranchAddress("evt", &evt);
  
  // -- data (overall) distribution
  TCanvas *c0 = getC0();
  int nEvent(0); 
  nEvent = t->GetEntries();
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    t->GetEntry(ievt);
    if (1 == classID) {
      ap0bgBDT->Fill(bdt0);
      ap1bgBDT->Fill(bdt1);
      ap2bgBDT->Fill(bdt2);
      ap3bgBDT->Fill(bdt);

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
      ap0sgBDT->Fill(bdt0);
      ap1sgBDT->Fill(bdt1);
      ap2sgBDT->Fill(bdt2);
      ap3sgBDT->Fill(bdt);

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

  ap0sgBDT->Scale(1./ap0sgBDT->GetSumOfWeights());
  tr0sgBDT->Scale(1./tr0sgBDT->GetSumOfWeights());
  te0sgBDT->Scale(1./te0sgBDT->GetSumOfWeights());

  ap1sgBDT->Scale(1./ap1sgBDT->GetSumOfWeights());
  tr1sgBDT->Scale(1./tr1sgBDT->GetSumOfWeights());
  te1sgBDT->Scale(1./te1sgBDT->GetSumOfWeights());

  ap2sgBDT->Scale(1./ap2sgBDT->GetSumOfWeights());
  tr2sgBDT->Scale(1./tr2sgBDT->GetSumOfWeights());
  te2sgBDT->Scale(1./te2sgBDT->GetSumOfWeights());

  ap3sgBDT->Scale(1./ap3sgBDT->GetSumOfWeights());
  ap3bgBDT->Scale(1./ap3bgBDT->GetSumOfWeights());

  double hmax = tr0bgBDT->GetMaximum(); 
  if (tr1bgBDT->GetMaximum() > hmax) hmax = tr1bgBDT->GetMaximum();
  if (tr2bgBDT->GetMaximum() > hmax) hmax = tr2bgBDT->GetMaximum();
  tr0bgBDT->SetMaximum(1.3*hmax);

  /*
  tr0bgBDT->SetMarkerStyle(20); 
  tr0bgBDT->SetMarkerSize(0.7); 
  tr0bgBDT->Draw("p");
  te0bgBDT->Draw("same");

  tr1bgBDT->SetMarkerStyle(20); 
  tr1bgBDT->SetMarkerSize(0.7); 
  tr1bgBDT->Draw("psame");
  te1bgBDT->Draw("same");

  tr2bgBDT->SetMarkerStyle(20); 
  tr2bgBDT->SetMarkerSize(0.7); 
  tr2bgBDT->Draw("psame");
  te2bgBDT->Draw("same");
  */

  c0->Clear();
  shrinkPad(0.15, 0.15); 
  setTitles(ap0bgBDT, "b", "a.u.", 0.05, 1.1, 1.5);
  ap0bgBDT->Draw();
  ap1bgBDT->Draw("same");
  ap2bgBDT->Draw("same");
  ap3bgBDT->Draw("same");

  newLegend(0.60, 0.67, 0.93, 0.87); 
  legg->SetTextSize(0.035);  
  legg->SetHeader("Background events"); 
  legg->AddEntry(ap3sgBDT, "combined", "l"); 
  legg->AddEntry(ap0sgBDT, "BDT 0", "l"); 
  legg->AddEntry(ap1sgBDT, "BDT 1", "l"); 
  legg->AddEntry(ap2sgBDT, "BDT 2", "l"); 
  legg->Draw();

  c0->SaveAs(Form("plots/%s-rebinned-bg-overlays.pdf", fname)); 

  hmax = tr0sgBDT->GetMaximum(); 
  if (tr1sgBDT->GetMaximum() > hmax) hmax = tr1sgBDT->GetMaximum();
  if (tr2sgBDT->GetMaximum() > hmax) hmax = tr2sgBDT->GetMaximum();
  tr0sgBDT->SetMaximum(1.3*hmax);

  /*
  tr0sgBDT->SetMarkerStyle(20); 
  tr0sgBDT->SetMarkerSize(0.7); 
  tr0sgBDT->Draw("p");
  te0sgBDT->Draw("same");

  tr1sgBDT->SetMarkerStyle(20); 
  tr1sgBDT->SetMarkerSize(0.7); 
  tr1sgBDT->Draw("psame");
  te1sgBDT->Draw("same");

  tr2sgBDT->SetMarkerStyle(20); 
  tr2sgBDT->SetMarkerSize(0.7); 
  tr2sgBDT->Draw("psame");
  te2sgBDT->Draw("same");
  */

  shrinkPad(0.15, 0.15); 
  setTitles(ap0sgBDT, "b", "a.u.", 0.05, 1.1, 1.5);
  ap0sgBDT->Draw();
  ap1sgBDT->Draw("same");
  ap2sgBDT->Draw("same");
  ap3sgBDT->Draw("same");

  newLegend(0.20, 0.67, 0.5, 0.87); 
  
  legg->SetTextSize(0.035);  
  legg->SetHeader("Signal events"); 
  legg->AddEntry(ap3sgBDT, "combined", "l"); 
  legg->AddEntry(ap0sgBDT, "BDT 0", "l"); 
  legg->AddEntry(ap1sgBDT, "BDT 1", "l"); 
  legg->AddEntry(ap2sgBDT, "BDT 2", "l"); 
  legg->Draw();

  c0->SaveAs(Form("plots/%s-rebinned-sg-overlays.pdf", fname));


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


  // -- new world: compute S, B, and SSB for each event category separately
  gStyle->SetOptTitle(1); 
  int vDmax[4], vDhiMax[4];
  double vSeffTot[4], vSeffMax[4], vSmax[4];
  double vDeffMax[4]; 
  double vSeffBinD99[4], vBmax[4], vMaxSSB[4], vMaxBDT[4];
  double vMaxSSBsimple[4], vMaxBDTsimple[4];
  double vMaxSSBfit[4], vMaxBDTfit[4];

  TH1D *hvsg[4], *hvbg[4];
  TH1D *hvsgAll[4], *hvbgAll[4];

  TH1D *hrocs[4];
  TH1D *hssb[4], *h2ssb[4];
  TGraph *grocs[4]; 
  TGraph *grocsOp[4]; 
  for (int ie = 0; ie < 4; ++ie) {
    vDhiMax[ie] = vDmax[ie] = 0; 

    vSeffTot[ie] = vSeffMax[ie] = vSmax[ie] = 
      vDeffMax[ie] = 
      vSeffBinD99[ie] = vBmax[ie] = vMaxSSB[ie] = vMaxBDT[ie] = 
      vMaxSSBsimple[ie] = vMaxBDTsimple[ie] = 
      vMaxSSBfit[ie] = vMaxBDTfit[ie] = 0; 

    hrocs[ie] = new TH1D(Form("hroc%d", ie), Form("hroc%d", ie), 200, 0., 1.); 
    hssb[ie]  = new TH1D(Form("hssb%d", ie), Form("hssb%d", ie), bdtBins, bdtMin, bdtMax); hssb[ie]->Sumw2();
    setTitles(hssb[ie], "b > ", "S/#sqrt{S+B}"); 
    h2ssb[ie] = new TH1D(Form("h2ssb%d", ie), "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); h2ssb[ie]->Sumw2();
    setTitles(h2ssb[ie], "b > ", "S/#sqrt{S+B'}"); 

    grocs[ie]   = new TGraph(bdtBins);   
    grocsOp[ie] = new TGraph(1); grocsOp[ie]->SetMarkerStyle(kOpenCross); grocsOp[ie]->SetMarkerSize(3.); grocsOp[ie]->SetMarkerColor(kBlack);

    grocs[ie]->SetMarkerStyle(20); grocs[ie]->SetMarkerSize(0.8);   
    if (0 == ie) grocs[ie]->SetMarkerColor(kMagenta); 
    if (1 == ie) grocs[ie]->SetMarkerColor(kRed); 
    if (2 == ie) grocs[ie]->SetMarkerColor(kBlue); 
    if (3 == ie) grocs[ie]->SetMarkerColor(kBlack); 

    hvsg[ie] = new TH1D(Form("hvsg%d", ie), Form("hvsg%d", ie), 100, 4.9, 5.9); hvsg[ie]->Sumw2();
    hvbg[ie] = new TH1D(Form("hvbg%d", ie), Form("hvbg%d", ie), 100, 4.9, 5.9); hvbg[ie]->Sumw2(); 

    hvsgAll[ie] = new TH1D(Form("hvsgAll%d", ie), Form("hvsgAll%d", ie), 100, 4.9, 5.9); hvsgAll[ie]->Sumw2();
    hvbgAll[ie] = new TH1D(Form("hvbgAll%d", ie), Form("hvbgAll%d", ie), 100, 4.9, 5.9); hvbgAll[ie]->Sumw2(); 
  }
  
  int ibin(0); 
  double bdtCut; 
  for (ibin = bdtBins; ibin >=0; --ibin) {
    bdtCut = bdtMin + ibin*(bdtMax-bdtMin)/bdtBins;
    cout << "+-+-+-+-+-+-+- bin " << ibin << " cutting at bdt > " << bdtCut << endl;
    for (int ie = 0; ie < 4; ++ie) {
      hvsg[ie]->Reset(); 
      hvbg[ie]->Reset(); 
      hvsgAll[ie]->Reset(); 
      hvbgAll[ie]->Reset(); 
    }

    // -- loop over tree for BDT cut
    nEvent = t->GetEntries();
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      t->GetEntry(ievt);
      if (0 == classID) {
	hvsgAll[3]->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) hvsgAll[0]->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) hvsgAll[1]->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) hvsgAll[2]->Fill(m, w8); 
      } else {
	hvbgAll[3]->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) hvbgAll[0]->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) hvbgAll[1]->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) hvbgAll[2]->Fill(m, w8); 
      }
      if (bdt <= bdtCut) continue;
      if (0 == classID) {
	hvsg[3]->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) hvsg[0]->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) hvsg[1]->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) hvsg[2]->Fill(m, w8); 
      } else {
	hvbg[3]->Fill(m, w8); 
	if (0 == TMath::Abs(evt%3)) hvbg[0]->Fill(m, w8); 
	if (1 == TMath::Abs(evt%3)) hvbg[1]->Fill(m, w8); 
	if (2 == TMath::Abs(evt%3)) hvbg[2]->Fill(m, w8); 
      }
    }

    // -- determine per event type the characteristics
    for (int ie = 0; ie < 4; ++ie) {
      double sCnt = hvsgAll[ie]->Integral(1, hvsgAll[ie]->GetNbinsX());
      double dCnt = hvbgAll[ie]->Integral(1, hvbgAll[ie]->GetNbinsX());
      
      h = hvsg[ie]; 
      double s = h->Integral(h->FindBin(5.3), h->FindBin(5.45));
      
      h = hvbg[ie]; 
      double d = h->Integral(1, h->GetNbinsX());
      double dhi = h->Integral(h->FindBin(5.45), h->GetNbinsX());
      double bsimple = d*(5.45-5.30)/(5.9-4.9-0.25);

      double seff = s/sCnt;
      double deff = bsimple/(dLostEvents + dCnt);
      double pbg  = 0.07*s;
      
      double b = bgBlind(h, 3, 4.9, 5.9);

      if ((1.-deff) > 0.999998) vSeffBinD99[ie] = seff;
      grocs[ie]->SetPoint(ibin, seff, 1.-deff); 
      hrocs[ie]->SetBinContent(hrocs[ie]->FindBin(seff), 1.-deff); 
      
      double ssb(0.); 
      cout << "** " << ie << "** bdt> " << bdtCut << " d = " << d << " s = " << s << " b = " << b
	   << " seff = " << seff << " deff = " << deff << endl;
      if (s+b+pbg >0) {
	ssb = s/TMath::Sqrt(s+b+pbg);
	cout << "** ** ssb " << ssb << endl;
	if (ssb > vMaxSSB[ie]) {
	  vSmax[ie] = s;
	  vDmax[ie] = static_cast<int>(d); 
	  vDhiMax[ie] = static_cast<int>(dhi); 
	  vBmax[ie] = b+pbg; 
	  vMaxSSB[ie] = ssb; 
	  vMaxBDT[ie] = bdtCut;
	  vSeffMax[ie] = seff;
	  vDeffMax[ie] = deff;
	  vSeffTot[ie] = s/(sLostEvents + sCnt); 
	}
	hssb[ie]->SetBinContent(ibin, ssb); 
      } else {
	hssb[ie]->SetBinContent(ibin, 0); 
      }

      double ssbs(0.); 
      if (s+bsimple >0) {
	ssbs = s/TMath::Sqrt(s+bsimple);
	cout << "** ** ssbs " << ssbs << endl;
	if (ssbs > vMaxSSBsimple[ie]) {
	  vMaxSSBsimple[ie] = ssbs; 
	  vMaxBDTsimple[ie] = bdtCut;
	}
	h2ssb[ie]->SetBinContent(ibin, ssbs); 
      } else {
	h2ssb[ie]->SetBinContent(ibin, 0); 
      }


      //    cout << "S = " << s << " B = " << b << " => S/sqrt(S+B) = " << s/TMath::Sqrt(s+b) << endl;
      c0->Clear();
      if (3 == ie) {
	hvbg[ie]->SetTitle(Form("evt type %d, bdt >%3.2f, S/B/D = %5.2f/%5.2f/%5.0f ssb = %4.3f/%4.3f", 
				ie, bdtCut, s, b, d, ssb, ssbs));
	hvbg[ie]->Draw("e");
	hvsg[ie]->Draw("samehist");
	c0->Modified();
	c0->Update();
      }
    }
  }

  for (int ie = 0; ie < 4; ++ie) {
    grocsOp[ie]->SetPoint(0, vSeffMax[ie], 1.-vDeffMax[ie]); 
  }

  // -- patch empty bins
  double okVal(1.); 
  for (int ie = 0; ie < 4; ++ie) {
    for (int i = 1; i< hrocs[ie]->GetNbinsX(); ++i) {
      if (hrocs[ie]->GetBinContent(i) < 1.e-5) hrocs[ie]->SetBinContent(i, okVal); 
      okVal = hrocs[ie]->GetBinContent(i);
    }

    for (int i = hrocs[ie]->GetNbinsX()-1; i > 0; --i) {
      if (!(hrocs[ie]->GetBinContent(i) < hrocs[ie]->GetBinContent(i-1))) {
	hrocs[ie]->SetBinContent(i, 0); 
      } else {
	break;
      }
    }
  }


  // -- fit for the maximum ssb and bdt cut
  initFunc *pFunc  = new initFunc(); 
  gStyle->SetOptFit(0); 
  for (int ie = 0; ie < 4; ++ie) {
    double xmax(0.), xmin(0.); 
    int nbins(0); 
    double maxVal = -1.;
    for (int i = 1; i < hssb[ie]->GetNbinsX(); ++i) 
      if (hssb[ie]->GetBinContent(i) > maxVal) maxVal = hssb[ie]->GetBinContent(i);
    
    for (int i = 1; i < hssb[ie]->GetNbinsX(); ++i) {
      if (hssb[ie]->GetBinContent(i) > 0.7*maxVal) {
	xmax = hssb[ie]->GetBinCenter(i); 
	nbins = i - hssb[ie]->GetMaximumBin(); 
      }
      hssb[ie]->SetBinError(i, 0.03*hssb[ie]->GetBinContent(i)); 
    }
    xmin = hssb[ie]->GetBinCenter(hssb[ie]->GetMaximumBin() - TMath::Abs(nbins)); 
    
    cout << "maxval: " << hssb[ie]->GetMaximum() << endl;
    cout << "maxbin: " << hssb[ie]->GetMaximumBin() << endl;
    cout << "xmax: " << xmax << endl;
    cout << "xmin: " << xmin << endl;
    cout << "nbins: " << nbins << endl;
    
    TF1 *f1 = pFunc->pol2local(hssb[ie], 0.05); 
    hssb[ie]->Fit(f1, "r", "", xmin, xmax); 
    double maxfitssbX = hssb[ie]->GetFunction("iF_pol2local")->GetParameter(2); 
    double maxfitssbY = hssb[ie]->GetFunction("iF_pol2local")->GetParameter(0); 
    vMaxSSBfit[ie] = maxfitssbX; 
    vMaxBDTfit[ie] = maxfitssbY; 
  }
  delete pFunc; 

  // -- BG overlays
  // --------------
  gStyle->SetOptTitle(0); 
  c0->Clear();
  TH2F* frame(0); 
  frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.65, 100, 0.99999, 1.000001);
  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
  
  string texname = string(fname) + ".tex";
  system(Form("/bin/rm -f %s", texname.c_str())); 


  int dMaxSum(0);
  double bMaxSum(0.), sMaxSum(0.); 
  double a(0.), dMaxBDT(0.), dMaxBDT3(0.); 

  dMaxBDT = TMath::Abs(vMaxBDT[0] - vMaxBDT[1]);
  a = TMath::Abs(vMaxBDT[0] - vMaxBDT[2]);
  if (a > dMaxBDT) dMaxBDT = a; 
  a = TMath::Abs(vMaxBDT[1] - vMaxBDT[2]);
  if (a > dMaxBDT) dMaxBDT = a; 

  a = 0.; 
  for (int ie = 0; ie < 4; ++ie) {

    if (ie < 3) {
      a = TMath::Abs(vMaxBDT[3] - vMaxBDT[ie]);
      if (a > dMaxBDT3) dMaxBDT3 = a; 
    }

    c0->Clear();
    frame->Draw();  
    gPad->SetLogy(0);
    gPad->SetLeftMargin(0.15);
    grocs[ie]->Draw("p"); 
    grocs[ie]->GetXaxis()->SetTitle("#epsilon_{ S}"); 
    grocs[ie]->GetYaxis()->SetTitle("1 - #epsilon_{ B}"); 
    double rocInt = hrocs[ie]->Integral(1, hrocs[ie]->GetNbinsX())*hrocs[ie]->GetBinWidth(1); 
    double rocInt2= hrocs[ie]->Integral(1, hrocs[ie]->FindBin(vSeffBinD99[ie]))*hrocs[ie]->GetBinWidth(1); 
    cout << " ==> seffBinD99 " << vSeffBinD99[ie] << " in bin " << hrocs[ie]->FindBin(vSeffBinD99[ie]) << endl;
    grocs[ie]->SetName("groc"); 
    grocs[ie]->SetTitle(Form("integral = %5.3f", rocInt));

    grocsOp[ie]->SetName(Form("grocsOp%d", ie)); 
    grocs[ie]->SetName(Form("groc%d", ie)); 

    grocs[ie]->Draw("p");
    grocsOp[ie]->Draw("p");

    int dMax = vDmax[ie];
    int dhiMax = vDhiMax[ie];
    double bMax = vBmax[ie];
    double sMax = vSmax[ie];
    double maxSSB = vMaxSSB[ie];
    double maxBDT = vMaxBDT[ie];
    double maxBDTsimple = vMaxBDTsimple[ie];
    double maxSSBsimple = vMaxSSBsimple[ie];
    double maxBDTfit = vMaxSSBfit[ie];
    double maxSSBfit = vMaxBDTfit[ie];
    double seffMax = vSeffMax[ie];
    double seffTot = vSeffTot[ie];

    if (ie < 3) {
      dMaxSum += dMax; 
      sMaxSum += sMax; 
      bMaxSum += bMax; 
    }

    tl->DrawLatex(0.25, 0.44, Form("D/B/S = %d/%2.1f/%2.1f", dMax, bMax, sMax)); 
    tl->DrawLatex(0.25, 0.40, Form("S/#sqrt{S+B}(MC) = %4.3f (%4.3f)", maxSSB, maxSSBsimple)); 
    tl->DrawLatex(0.25, 0.36, Form("b_{max}(MC) = %4.3f (%4.3f)", maxBDT, maxBDTsimple)); 
    tl->DrawLatex(0.25, 0.32, Form("#epsilon_{BDT} = %4.3f", seffMax)); 
    tl->DrawLatex(0.25, 0.28, Form("#epsilon_{tot} = %6.5f", seffTot)); 
    tl->DrawLatex(0.25, 0.24, Form("I_{tot} = %6.5f", rocInt)); 
    tl->DrawLatex(0.25, 0.20, Form("I_{part} = %6.5f", rocInt2)); 
    
    ofstream TEX(texname.c_str(), ios::app);
    
    TEX << Form("\\vdef{s%s:ie%d:string}       {%s}", fname, ie, fname) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssb}          {%4.3f}", fname, ie, maxSSB) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdt}       {%4.3f}", fname, ie, maxBDT) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssbsimple}    {%4.3f}", fname, ie, maxSSBsimple) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdtsimple} {%4.3f}", fname, ie, maxBDTsimple) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssbfit}       {%4.3f}", fname, ie, maxSSBfit) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdtfit}    {%4.3f}", fname, ie, maxBDTfit) << endl;
    TEX << Form("\\vdef{s%s:ie%d:Smc}          {%4.3f}", fname, ie, sMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:D}            {%d}", fname, ie, dMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:Dhi}          {%d}", fname, ie, dhiMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:B}            {%4.3f}", fname, ie, bMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ipart}        {%6.5f}", fname, ie, rocInt2) << endl;
    TEX << Form("\\vdef{s%s:ie%d:itot}         {%6.5f}", fname, ie, rocInt) << endl;
    TEX << Form("\\vdef{s%s:ie%d:epstot}       {%6.5f}", fname, ie, seffTot) << endl;
    TEX << Form("\\vdef{s%s:ie%d:epsbdt}       {%6.5f}", fname, ie, seffMax) << endl;
    if (ie == 3) {
      TEX << Form("\\vdef{s%s:BDTparameters}  {%s}", fname, fBDTParameters.c_str()) << endl;
      TEX << Form("\\vdef{s%s:BDTvariables}  {%s}", fname, fVariables.c_str()) << endl;
      TEX << Form("\\vdef{s%s:sum:Smc}     {%4.3f}", fname, sMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:D}       {%d}", fname, dMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:B}       {%4.3f}", fname, bMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:ssb}     {%4.3f}", fname, sMaxSum/sqrt(sMaxSum+bMaxSum)) << endl;

      double ks0 = trainBDT0->KolmogorovTest(testBDT0);
      double ks1 = trainBDT1->KolmogorovTest(testBDT1);
      double ks2 = trainBDT2->KolmogorovTest(testBDT2);
      TEX << Form("\\vdef{s%s:ie0:ksBg}    {%4.3f}", fname, ks0) << endl;
      TEX << Form("\\vdef{s%s:ie1:ksBg}    {%4.3f}", fname, ks1) << endl;
      TEX << Form("\\vdef{s%s:ie2:ksBg}    {%4.3f}", fname, ks2) << endl;

      ks0 = trainSgBDT0->KolmogorovTest(testSgBDT0);
      ks1 = trainSgBDT1->KolmogorovTest(testSgBDT1);
      ks2 = trainSgBDT2->KolmogorovTest(testSgBDT2);
      TEX << Form("\\vdef{s%s:ie0:ksSg}    {%4.3f}", fname, ks0) << endl;
      TEX << Form("\\vdef{s%s:ie1:ksSg}    {%4.3f}", fname, ks1) << endl;
      TEX << Form("\\vdef{s%s:ie2:ksSg}    {%4.3f}", fname, ks2) << endl;
      TEX << Form("\\vdef{s%s:dMaxBDT}     {%4.3f}", fname, dMaxBDT) << endl;
      TEX << Form("\\vdef{s%s:dMaxBDT3}    {%4.3f}", fname, dMaxBDT3) << endl;
    }

    TEX.close();
    system(Form("/bin/cp %s plots", texname.c_str())); 
   
    newLegend(0.22, 0.47, 0.55, 0.67); 
    
    legg->SetTextSize(0.035);  
    legg->AddEntry(grocsOp[ie], Form("operating point b > %4.3f", maxBDT), "p"); 
    if (ie < 3) {
      legg->AddEntry(grocs[ie],  Form("BDT %d", ie), "p"); 
    } else {
      legg->AddEntry(grocs[ie],  "combined", "p"); 
    }
    legg->Draw();
    
    c0->SaveAs(Form("plots/%s-roc-ie%d.pdf", fname, ie)); 
    
    c0->Clear();
    hssb[ie]->Draw();
    h2ssb[ie]->Draw("same");
    tl->DrawLatex(0.2, 0.90, Form("event type %d", ie)); 
    tl->DrawLatex(0.2, 0.85, Form("SSB_{max} = %4.3f (%4.3f/%4.3f)", maxSSB, maxSSBsimple, maxSSBfit)); 
    tl->DrawLatex(0.2, 0.80, Form("BDT_{max} > %4.3f (%4.3f/%4.3f)", maxBDT, maxBDTsimple, maxBDTfit)); 
    tl->DrawLatex(0.2, 0.75, Form("ROC_{int} = %4.3f", rocInt)); 
    c0->SaveAs(Form("plots/%s-ssb-ie%d.pdf", fname, ie)); 
    
    cout << "Write out SSB histograms" << endl;
    cout << "  maxSSB: " << maxSSB << " at BDT > " << maxBDT << endl;
    hssb[ie]->SetDirectory(f);

    f->cd();
    grocsOp[ie]->Write();
    grocs[ie]->Write();
  }


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
                          countCanvas*50+200, countCanvas*20, width, static_cast<int>(width*0.78) ); 
    
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
	 cout << "legend->AddEntry(sig...)" << endl;
         legend->AddEntry(sig, TString("Signal")     + " (test sample)", "F");
         legend->AddEntry(bgd, TString("Background") + " (test sample)", "F");
         legend->SetBorderSize(1);
         legend->SetMargin(0.2);
         legend->Draw();

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
	   legend2->Draw();
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
	 
	 //          // text for overflows
	 //          Int_t    nbin = sig->GetNbinsX();
	 //          Double_t dxu  = sig->GetBinWidth(0);
	 //          Double_t dxo  = sig->GetBinWidth(nbin+1);
	 //          TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
	 //                                 sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
	 //                                 sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
	 TText* t = new TText( 0.975, 0.115, fname.c_str());
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
	if (0 == ibin) cout << "===> bin == 0: " << varn << endl;
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
      filename = "/scratch/ursl/bdt/tmva-trees-0-2011.root"; 
    } 
    if (2012 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-0-2012.root"; 
    } 
  }


  make(offset, filename, 0, clean);
  make(offset, filename, 1, clean);
  make(offset, filename, 2, clean);
  
  string oname = Form("TMVA-%d", offset); 
  cout << "-->apply(" << oname.c_str() << ")" << endl;
  apply(oname.c_str());
  cout << "-->analyze(" << oname.c_str() << ")" << endl;
  analyze(oname.c_str()); 

  cout << "-->mvas(...)" << endl;
  string sEvents = oname + "-Events0";
  mvas(sEvents.c_str()); 
  sEvents = oname + "-Events1";
  mvas(sEvents.c_str()); 
  sEvents = oname + "-Events2";
  mvas(sEvents.c_str()); 

  if (clean) {
    cout << "-->cleanup(...)" << endl;
    cleanup(oname);
  }

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
  cout << "==> tmva1(" << oname << ", " << filename << ") " << endl;
  cout << "======================================================================" << endl;

  cout << "-->train(...) with oname = " << oname << " and filename = " << filename << endl;
  train(oname, filename);
}


// ----------------------------------------------------------------------
void tmva1::createInputFile(string filename, int randomSeed) {
  TFile *sinput = TFile::Open(fInputFiles.sname.c_str());
  TFile *dinput = TFile::Open(fInputFiles.dname.c_str());

  TCut sgcut = preselection().c_str(); 
  
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
  if (randomSeed > -1) gRandom->SetSeed(randomSeed);
  for (int j = 0; j < 3; ++j) {
    if (0 == j) {
      type = "Events0"; 
      typeCut = "TMath::Abs(evt%3)==0";
      if (randomSeed > -1) typeCut = "3*rndm%3==0";
    } else if (1 == j) {
      type = "Events1";
      typeCut = "TMath::Abs(evt%3)==1";
      if (randomSeed > -1) typeCut = "3*rndm%3==1";
    } else if (2 == j) {
      type = "Events2";
      typeCut = "TMath::Abs(evt%3)==2";
      if (randomSeed > -1) typeCut = "3*rndm%3==2";
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
  frd.chi2dof = ftd.chi2dof;
  frd.closetrk = ftd.closetrk; 
  frd.closetrks1 = ftd.closetrks1; 
  frd.closetrks2 = ftd.closetrks2; 
  frd.closetrks3 = ftd.closetrks3; 

  frd.m1iso = ftd.m1iso; 
  frd.m2iso = ftd.m2iso; 
  frd.pvdchi2 = ftd.pvdchi2; 
  frd.othervtx = ftd.othervtx; 

  frd.pvlips2 = ftd.pvlips2; 
  frd.pvlip2 = ftd.pvlip2; 
  
  frd.m  = ftd.m; 
  int ichan = 0; 

  if (ftd.evt < 0) {
    cout << "XXXXXXXXXXXXXXXXXXXXX event number still smaller than zero!!!!!" << endl;
  }
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
	if (stype == "closetrks1") {
	  cout << "  adding closetrks1" << endl;
	  reader->AddVariable( "closetrks1", &rd.closetrks1);
	}
	if (stype == "closetrks2") {
	  cout << "  adding closetrks2" << endl;
	  reader->AddVariable( "closetrks2", &rd.closetrks2);
	}
	if (stype == "closetrks3") {
	  cout << "  adding closetrks3" << endl;
	  reader->AddVariable( "closetrks3", &rd.closetrks3);
	}
	if (stype == "chi2dof") {
	  cout << "  adding chi2dof" << endl;
	  reader->AddVariable( "chi2dof", &rd.chi2dof);
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
	if (stype == "othervtx") {
	  cout << "  adding othervtx" << endl;
	  reader->AddVariable( "othervtx", &rd.othervtx);
	}
	if (stype == "pvlips2") {
	  cout << "  adding pvlips2" << endl;
	  reader->AddVariable( "pvlips2", &rd.pvlips2);
	}
	if (stype == "pvlip2") {
	  cout << "  adding pvlip2" << endl;
	  reader->AddVariable( "pvlip2", &rd.pvlip2);
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
void shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}


// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(static_cast<Width_t>(width));
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
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
  tree->Branch("closetrks1", &rd.closetrks1, "closetrks1/I");
  tree->Branch("closetrks2", &rd.closetrks2, "closetrks2/I");
  tree->Branch("closetrks3", &rd.closetrks3, "closetrks3/I");

  tree->Branch("m1iso", &rd.m1iso, "m1iso/F");
  tree->Branch("m2iso", &rd.m2iso, "m2iso/F");
  tree->Branch("pvdchi2", &rd.pvdchi2, "pvdchi2/F");
  tree->Branch("othervtx", &rd.othervtx, "othervtx/F");
  tree->Branch("pvlips2", &rd.pvlips2, "pvlips2/F");
  tree->Branch("pvlip2", &rd.pvlip2, "pvlip2/F");

  return tree; 
}


// ----------------------------------------------------------------------
double tmva1::bgBlind(TH1 *h, int mode, double lo, double hi) {
  
  TVirtualFitter::SetMaxIterations(20000);

  if (0 == h) { 
    cout << "tmva1::bgBlind(...): No histogram passed! mode = " << mode << endl;
    return -1.;
  }
  
  TF1 *lF1(0), *lF2(0);

  initFunc *pFunc  = new initFunc(); 

  double BgLo = 4.9; 
  double BgHi = 5.9; 
  double histCount = h->Integral(h->FindBin(BgLo+0.0001), h->FindBin(BgHi-0.0001)); 
  cout << "bgBlind: histCount = " << histCount 
       << " starting at " << BgLo+0.0001 << " bin(" << h->FindBin(BgLo+0.0001) << ")"
       << " to " << BgHi-0.0001 << " bin(" << h->FindBin(BgHi-0.0001) << ")" 
       << " mode: " << mode 
       << endl;

  double BgHist  = histCount; 
  double BgHistE  = TMath::Sqrt(histCount); 
  if (histCount > 0) {
    BgHistE = TMath::Sqrt(histCount)/histCount*BgHist;
  } else {
    BgHistE  = 0.2; // FIXME?!
    return 0.;
  }

  if (3 == mode) {
    pFunc->resetLimits(); 
    lF1 = pFunc->pol1BsBlind(h); 
    lF2 = pFunc->pol1(4.9, 5.9); 
  } else {
    cout << " implement missing code!" << endl;
  }
  
  lF2->SetLineStyle(kDashed);
  h->Fit(lF1, "rl", "", lo, hi); 
  h->DrawCopy();
  lF2->SetLineColor(kBlue);
  lF2->Draw("same");
  lF2->SetParameters(lF1->GetParameters());
  lF2->SetParErrors(lF1->GetParErrors());

  double BsBgExp; 
  if (!strcmp(gMinuit->fCstatu.Data(), "CONVERGED ")) {
    lF2->Update();
    double integral = lF2->Integral(5.30, 5.45, static_cast<const Double_t*>(0), 1.e-15);
    BsBgExp  = integral/h->GetBinWidth(1); 
  } else {
    BsBgExp  = (5.45-5.30)/(5.9-4.9-0.25)*BgHist;
    cout << "+++ Fit did not converge, take flat bg interpretation, fCstatu = ->" << gMinuit->fCstatu.Data() << "<-" 
	 << ", BsBgExp = " << BsBgExp
	 << endl;
  }

  delete lF1; 
  delete lF2; 
  delete pFunc; 

  return BsBgExp;
}



// ----------------------------------------------------------------------
void tmva1::createToyData(string sgfilename, string bgfilename, string ofilename, int seed, int nsg, int nbg) {
  delete gRandom; 
  gRandom = new TRandom3(seed); 

  cout << "==> CREATE TOY DATA for seed = " << seed << endl;

  int channel = 0; 

  // -- define variables
  vector<string> vNames;
  vector<double> vMin, vMax; 
  vector<int> vNbins; 
  vNames.push_back("m1pt"); vMin.push_back(0.); vMax.push_back(60.); vNbins.push_back(120); 
  vNames.push_back("m2pt"); vMin.push_back(0.); vMax.push_back(40.); vNbins.push_back(80); 
  vNames.push_back("m1eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("m2eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("pt");    vMin.push_back(0.); vMax.push_back(100.); vNbins.push_back(200); 
  vNames.push_back("eta");   vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); 
  vNames.push_back("fls3d"); vMin.push_back(0.); vMax.push_back(150.); vNbins.push_back(150); 
  vNames.push_back("alpha"); vMin.push_back(0.); vMax.push_back(1.); vNbins.push_back(100); 
  vNames.push_back("maxdoca"); vMin.push_back(0.); vMax.push_back(0.1); vNbins.push_back(200); 
  vNames.push_back("pvip"); vMin.push_back(0.); vMax.push_back(0.1); vNbins.push_back(100); 
  vNames.push_back("pvips"); vMin.push_back(0.); vMax.push_back(5); vNbins.push_back(100); 
  vNames.push_back("iso"); vMin.push_back(0.); vMax.push_back(1.01); vNbins.push_back(101); 
  vNames.push_back("m1iso"); vMin.push_back(0.); vMax.push_back(1.01); vNbins.push_back(101); 
  vNames.push_back("m2iso"); vMin.push_back(0.); vMax.push_back(1.01); vNbins.push_back(101); 
  vNames.push_back("closetrk"); vMin.push_back(0.); vMax.push_back(21); vNbins.push_back(21); 
  vNames.push_back("docatrk"); vMin.push_back(0.); vMax.push_back(0.5); vNbins.push_back(100); 
  vNames.push_back("chi2dof"); vMin.push_back(0.); vMax.push_back(10); vNbins.push_back(100); 

  TFile *sgFile = TFile::Open(sgfilename.c_str());
  TTree *tsg = (TTree*)sgFile->Get("candAnaMuMu/events");
  TFile *bgFile = TFile::Open(bgfilename.c_str());
  TTree *tbg = (TTree*)bgFile->Get("candAnaMuMu/events");

  // -- create histograms
  TH1D *hs, *hb; 
  string hsName, hbName; 
  string presel = preselection();
  TCanvas *c0 = getC0();
  c0->Clear();
  c0->Divide(6,6);
  sgFile->cd();
  cout << "Preselection: " 
       << endl
       << presel
       << endl;
  for (unsigned int i = 0; i < vNames.size(); ++i) {
    hsName = Form("hs_%s", vNames[i].c_str());
    hs  = new TH1D(hsName.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);
    setHist(hs, kBlue); 

    hbName = Form("hb_%s", vNames[i].c_str());
    hb  = new TH1D(hbName.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);
    setHist(hb, kRed); 
    
    string var = vNames[i]; 
    cout << "--> " << var <<  endl;
    tsg->Draw(Form("%s>>%s", var.c_str(), hsName.c_str()), presel.c_str(), "goff");
    tbg->Draw(Form("%s>>%s", var.c_str(), hbName.c_str()), presel.c_str(), "goff");
    c0->cd(2*i+1); 
    hs->Draw("hist"); 
    c0->cd(2*i+2); 
    hb->Draw("hist");
    c0->Modified(); 
    c0->Update();
  }

  c0->SaveAs(Form("toyData-%d-%d-%d.pdf", seed, nsg, nbg)); 
  
  c0->cd(vNames.size()+1);

  // -- fill and dump toy trees
  TFile *outfile = TFile::Open(ofilename.c_str(), "RECREATE");  
  struct readerData rd; 

  TTree *ts[3], *tb[3]; 
  for (int ie = 0; ie < 3; ++ie) {
    outfile->mkdir(Form("signalChan0Events%d", ie));
    outfile->cd(Form("signalChan0Events%d", ie));
    ts[ie] = createTree(rd); 
    for (int i = 0 ; i < nsg; ++i) {
      rd.m = 5.37;
      rd.pt = ((TH1D*)sgFile->Get("hs_pt"))->GetRandom();
      rd.eta = ((TH1D*)sgFile->Get("hs_eta"))->GetRandom();
      rd.m1eta = ((TH1D*)sgFile->Get("hs_m1eta"))->GetRandom();
      rd.m2eta = ((TH1D*)sgFile->Get("hs_m2eta"))->GetRandom();
      rd.m1pt = ((TH1D*)sgFile->Get("hs_m1pt"))->GetRandom();
      rd.m2pt = ((TH1D*)sgFile->Get("hs_m2pt"))->GetRandom();
      
      rd.fls3d = ((TH1D*)sgFile->Get("hs_fls3d"))->GetRandom();
      rd.alpha = ((TH1D*)sgFile->Get("hs_alpha"))->GetRandom();
      rd.maxdoca = ((TH1D*)sgFile->Get("hs_maxdoca"))->GetRandom();
      rd.pvip = ((TH1D*)sgFile->Get("hs_pvip"))->GetRandom();
      rd.pvips = ((TH1D*)sgFile->Get("hs_pvips"))->GetRandom();
      rd.iso = ((TH1D*)sgFile->Get("hs_iso"))->GetRandom();
      rd.m1iso = ((TH1D*)sgFile->Get("hs_m1iso"))->GetRandom();
      rd.m2iso = ((TH1D*)sgFile->Get("hs_m2iso"))->GetRandom();
      rd.docatrk = ((TH1D*)sgFile->Get("hs_docatrk"))->GetRandom();
      rd.chi2dof = ((TH1D*)sgFile->Get("hs_chi2dof"))->GetRandom();
      rd.closetrk = ((TH1D*)sgFile->Get("hs_closetrk"))->GetRandom();
      ts[ie]->Fill();
    }
  }

  for (int ie = 0; ie < 3; ++ie) {
    outfile->mkdir(Form("sidebandChan0Events%d", ie));
    outfile->cd(Form("sidebandChan0Events%d", ie));

    tb[ie] = createTree(rd); 
    for (int i = 0 ; i < nbg; ++i) {
      rd.m = 5.37;
      rd.pt = ((TH1D*)sgFile->Get("hb_pt"))->GetRandom();
      rd.eta = ((TH1D*)sgFile->Get("hb_eta"))->GetRandom();
      rd.m1eta = ((TH1D*)sgFile->Get("hb_m1eta"))->GetRandom();
      rd.m2eta = ((TH1D*)sgFile->Get("hb_m2eta"))->GetRandom();
      rd.m1pt = ((TH1D*)sgFile->Get("hb_m1pt"))->GetRandom();
      rd.m2pt = ((TH1D*)sgFile->Get("hb_m2pt"))->GetRandom();
      
      rd.fls3d = ((TH1D*)sgFile->Get("hb_fls3d"))->GetRandom();
      rd.alpha = ((TH1D*)sgFile->Get("hb_alpha"))->GetRandom();
      rd.maxdoca = ((TH1D*)sgFile->Get("hb_maxdoca"))->GetRandom();
      rd.pvip = ((TH1D*)sgFile->Get("hb_pvip"))->GetRandom();
      rd.pvips = ((TH1D*)sgFile->Get("hb_pvips"))->GetRandom();
      rd.iso = ((TH1D*)sgFile->Get("hb_iso"))->GetRandom();
      rd.m1iso = ((TH1D*)sgFile->Get("hb_m1iso"))->GetRandom();
      rd.m2iso = ((TH1D*)sgFile->Get("hb_m2iso"))->GetRandom();
      rd.docatrk = ((TH1D*)sgFile->Get("hb_docatrk"))->GetRandom();
      rd.chi2dof = ((TH1D*)sgFile->Get("hb_chi2dof"))->GetRandom();
      rd.closetrk = ((TH1D*)sgFile->Get("hb_closetrk"))->GetRandom();
      tb[ie]->Fill();
    }
  }

  for (int ie = 0; ie < 3; ++ie) {
    ts[ie]->Write();
    tb[ie]->Write();
  }

  outfile->Write();
  outfile->Close();

}


// ----------------------------------------------------------------------
void tmva1::toyRun(string modifier, string vars, string bdtpars, int seed, int nsg, int nbg) {
  string oname; 

  string sgfilename = "/scratch/ursl/bdt/v16-2011-mix-Bs2MuMu.root"; 
  string bgfilename = "/scratch/ursl/bdt/v16-2011-data-bmmLoose-1.root";
  
  if (!strcmp(vars.c_str(), "")) {
    fVariables = "pt:eta:alpha:chi2dof:maxdoca:fls3d:pvip:pvips:closetrk:iso:m1iso:m2iso"; 
  } else {
    fVariables = vars; 
  }

  if (!strcmp(bdtpars.c_str(), "")) {
    fBDTParameters = ":NTrees=200:nCuts=20:BoostType=AdaBoost:AdaBoostBeta=1.0:MaxDepth=3:NNodesMax=100000:nEventsMin=50"; 
  } else {
    fBDTParameters = bdtpars; 
  }

  oname = Form("/scratch/ursl/toys/toy-%d.root", seed); 

  ifstream ifile(oname.c_str());
  if (!ifile) {
    createToyData(sgfilename, bgfilename, oname, seed, nsg, nbg); 
  } else {
    ifile.close(); 
  }
  

  if (1) {
    setApply0(); train(Form("toy-%d-0", seed), oname, nsg, nbg); 
    setApply1(); train(Form("toy-%d-1", seed), oname, nsg, nbg); 
    setApply2(); train(Form("toy-%d-2", seed), oname, nsg, nbg); 
  }

  TCanvas *c0 = getC0();
  c0->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double ymax = 0.5; 
  TH2F* frame = new TH2F("frame", "BDT output distributions", 100, -1., 1., 100, 0., ymax);
  frame->GetXaxis()->SetTitle(" b");
  frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
  frame->Draw();  

  int color(0);
  for (int ie = 0; ie < 3; ++ie) {
    oname = Form("toy-%d-%d.root", seed, ie);
    TFile *file = TFile::Open(oname.c_str()); 
    
    if (0 == ie) color = kMagenta;
    if (1 == ie) color = kRed;
    if (2 == ie) color = kBlue;
    TH1F *hTrainSg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_S"); 
    hTrainSg->SetLineColor(color); 
    hTrainSg->SetMarkerColor(color);
    hTrainSg->Scale(1./hTrainSg->GetSumOfWeights());
    c0->cd(); 
    hTrainSg->DrawCopy("histesame");

    TH1F *hTrainBg = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_B"); 
    hTrainBg->SetLineColor(color); 
    hTrainBg->SetMarkerColor(color);
    hTrainBg->Scale(1./hTrainBg->GetSumOfWeights());
    hTrainBg->DrawCopy("histesame");

    TPad *t = new TPad(Form("pad%d", ie), Form("pad%d", ie), 0.2+0.2*ie, 0.4, 0.2+0.2*ie+0.2, 0.6); 
    t->Draw();
    t->cd(); 
    TH1F *hs0 = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_S"); hs0->SetMarkerColor(color); 
    TH1F *hs1 = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_S"); hs1->SetLineColor(color);
    hs0->DrawNormalized("e");
    hs1->DrawNormalized("histsame");
    Double_t kolS = hs1->KolmogorovTest(hs0);

    TH1F *hb0 = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_Train_B"); hb0->SetMarkerColor(color); 
    TH1F *hb1 = (TH1F*)file->Get("Method_BDT/BDT/MVA_BDT_B"); hb1->SetLineColor(color);
    hb0->DrawNormalized("esame");
    hb1->DrawNormalized("histsame");
    Double_t kolB = hb1->KolmogorovTest(hb0);


    tl->SetTextSize(0.1);
    tl->DrawLatex(0.2, 0.92, Form("%4.3f", kolB)); 
    tl->DrawLatex(0.6, 0.92, Form("%4.3f", kolS)); 
    
    file->Close();
  }
  
  c0->cd();
  tl->SetTextSize(0.02); 
  tl->DrawLatex(0.13, 0.92, Form("toys-output-%d-%s", seed, modifier.c_str())); 
  tl->SetTextSize(0.02); 
  tl->DrawLatex(0.13, 0.80, Form("Nsg/Nbg %d/%d", nsg, nbg)); 
  tl->DrawLatex(0.13, 0.76, fVariables.c_str());  
  tl->SetTextSize(0.014); 
  tl->DrawLatex(0.13, 0.72, fBDTParameters.c_str());  
  c0->SaveAs(Form("toys-output-%d-%s.pdf", seed, modifier.c_str()));
}


// ----------------------------------------------------------------------
void tmva1::analyzeTexFiles(std::string dir, int start, int end, string what) {

  int ivar(0); 
  if (string::npos != what.find("ssbfit")) {
    ivar = 1; 
  } else   if (string::npos != what.find("ssbsimple")) {
    ivar = 2; 
  } else if (string::npos != what.find("ssb")) {
    ivar = 3; 
  }

  vector<string> lines; 
  double best_x(-1.); 
  int best_idx(-1);

  const double KSCUT(0.1); 
  int combD, combDhi, sumD, sumDhi;
  double x, ssbfit, ssb, ssbs, combBg, combSg, sumBg, sumSg, sssb; 
  double bks0, bks1, bks2, sks0, sks1, sks2; 
  double dMaxBDT, dMaxBDT3; 
  string setup;
  for (int i = start; i <= end; ++i) {
    lines.clear(); 
    readTexFile(Form("%s/tmp-%d/TMVA-%d.tex", dir.c_str(), i, i), lines); 
    sumD = sumDhi = combD = combDhi = -99;
    x = ssbfit = ssb = ssbs = combBg = combSg = sumBg = sumSg = sssb = -99.;
    sks0 = sks1 = sks2 = bks0 = bks1 = bks2 = -99.;
    dMaxBDT = dMaxBDT3 = 99.; 
    setup = ""; 

    for (unsigned int j = 0; j < lines.size(); ++j) {
      if (string::npos != lines[j].find("ie0:ksBg}")) {bks0 = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie1:ksBg}")) {bks1 = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie2:ksBg}")) {bks2 = parseTexLine(lines[j]); continue;}

      if (string::npos != lines[j].find("ie0:ksSg}")) {sks0 = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie1:ksSg}")) {sks1 = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie2:ksSg}")) {sks2 = parseTexLine(lines[j]); continue;}
      
      if (string::npos != lines[j].find("sum:ssb}")) {sssb = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("sum:Smc}")) {sumSg = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("sum:D}"))   {sumD = static_cast<int>(parseTexLine(lines[j])); continue;}
      if (string::npos != lines[j].find("sum:Dhi}")) {sumDhi = static_cast<int>(parseTexLine(lines[j])); continue;}
      if (string::npos != lines[j].find("sum:B}"))   {sumBg = parseTexLine(lines[j]); continue;}

      if (string::npos != lines[j].find("ie3:Smc}")) {combSg = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie3:Dhi}")) {combDhi = static_cast<int>(parseTexLine(lines[j])); continue;}
      if (string::npos != lines[j].find("ie3:D}"))   {combD = static_cast<int>(parseTexLine(lines[j])); continue;}
      if (string::npos != lines[j].find("ie3:B}"))   {combBg = parseTexLine(lines[j]); continue;}

      if (string::npos != lines[j].find("dMaxBDT}"))   {dMaxBDT = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("dMaxBDT3}"))   {dMaxBDT3 = parseTexLine(lines[j]); continue;}

      if (string::npos != lines[j].find("ie3:ssb}")) {ssb = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie3:ssbfit}")) {ssbfit = parseTexLine(lines[j]); continue;}
      if (string::npos != lines[j].find("ie3:ssbsimple}")) {ssbs = parseTexLine(lines[j]); continue;}

      if (string::npos != lines[j].find("BDTparameters}")) {
	setup += lines[j].substr(lines[j].find("BDTparameters} ")+string("BDTparameters} ").length()); 
	continue;
      }
      if (string::npos != lines[j].find("BDTvariables}")) {
	setup += lines[j].substr(lines[j].find("BDTvariables} ")+string("BDTvariables} ").length()); 
	continue;
      }
  
    }


    if (bks0 < KSCUT || bks1 < KSCUT || bks2 < KSCUT || sks0 < KSCUT || sks1 < KSCUT || sks2 < KSCUT) {
      //       cout << "skipping " << i << " because of low KS probs: " 
      // 	   << bks0 << " " << bks1 << " " << bks2 << " " << sks0 << " " << sks1 << " " << sks2 
      // 	   << endl;
      continue;
    }

    if (1 == ivar) {
      x = ssbfit;
    } else if (2 == ivar) {
      x = ssbs;
    } else if (3 == ivar) {
      x = ssb;
    }

    if (x > best_x) {
      cout << "at " << i << " new better " << what << ": " << ssb << "/" << ssbs << "/" << ssbfit 
	   << "S/B/D = " << combSg << "/" << combBg << "/" << combD
	   << " dmaxBDT = " << dMaxBDT << "/" << dMaxBDT3 << " " << setup << endl;
      best_x = x; 
      best_idx = i; 
    } else if (x > 0.99*best_x) {
      cout << "   at " << i << " similar " << what << ": " << ssb << "/" << ssbs << "/" << ssbfit 
	   << "S/B/D = " << combSg << "/" << combBg << "/" << combD
	   << " dmaxBDT = " << dMaxBDT << "/" << dMaxBDT3 << " " << setup << endl;
    }
    
  }
    
}

// ----------------------------------------------------------------------
float tmva1::parseTexLine(string line) {
  //  cout << "0:" << line << endl;
  string::size_type m1 = line.rfind("{"); 
  string::size_type m2 = line.rfind("}"); 
  string stype = line.substr(m1+1, m2-m1-1); 
  //  cout << "1:" << stype << endl;
  return atof(stype.c_str());
}

// ----------------------------------------------------------------------
void tmva1::readTexFile(string filename, vector<string> &lines) {
  ifstream is(filename.c_str());
  char  buffer[200];
  while (is.getline(buffer, 200, '\n')) {
    lines.push_back(string(buffer));
  }
}
