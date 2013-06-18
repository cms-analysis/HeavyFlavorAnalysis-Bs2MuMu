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

#include "tmva10.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)#
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

ClassImp(tmva10)

using namespace std; 

// ----------------------------------------------------------------------
// -- 
// -- USAGE: a.makeAll(0, 1); > TMVA-0.log
// --
// ----------------------------------------------------------------------

tmva10::tmva10(string vars) {

  cout << "tmva10 hello: setup with variables: " << vars << endl;

  fVariables = vars; 
  fBDTParameters = ""; 

  legg = 0;
  legge = 0; 
  tl = new TLatex(); 
  tl->SetTextFont(42);
  tl->SetTextSize(0.03);
  tl->SetNDC(kTRUE); 

  // -- TMVA default
  fBdtSetup.NTrees = 200;
  fBdtSetup.nEventsMin = 500; 
  fBdtSetup.MaxDepth = 3;
  fBdtSetup.nCuts = 20; 
  fBdtSetup.AdaBoostBeta = 1.0; 
  fBdtSetup.NNodesMax = 100000;

}


// ----------------------------------------------------------------------
tmva10::~tmva10() {
  cout << "tmva10 good bye " << endl;
}


// ----------------------------------------------------------------------
TCanvas* tmva10::getC0() {
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
  if (0 == c0) c0 = new TCanvas("c0","--c0--",2303,0,656,700);
  return c0; 
} 

// ----------------------------------------------------------------------
void tmva10::train(string oname, string iname) {
   // This loads the library
   TMVA::Tools::Instance();
   
   (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40; 
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 40; 
   
   // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName(Form("%s.root", oname.c_str()));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TFile *iFile = TFile::Open(iname.c_str());
   TTree *trainSg(0), *testSg(0), *trainBg(0), *testBg(0); 
   
   cout << "==============> Train on events0, test on events1" << endl;
   trainSg = (TTree*)iFile->Get("sg0/muonidtree");
   testSg  = (TTree*)iFile->Get("sg1/muonidtree");
   trainBg = (TTree*)iFile->Get("bg0/muonidtree");
   testBg  = (TTree*)iFile->Get("bg1/muonidtree");

   
   TH1D *hSetup = new TH1D("hSetup", "hSetup", 100, 0., 100.); 
   int i(0); 
   i = 10; hSetup->SetBinContent(i, fBdtSetup.NTrees); hSetup->GetXaxis()->SetBinLabel(i, "NTrees");
   i = 11; hSetup->SetBinContent(i, fBdtSetup.nEventsMin); hSetup->GetXaxis()->SetBinLabel(i, "nEventsMin");
   i = 12; hSetup->SetBinContent(i, fBdtSetup.nCuts); hSetup->GetXaxis()->SetBinLabel(i, "nCuts");
   i = 13; hSetup->SetBinContent(i, fBdtSetup.AdaBoostBeta); hSetup->GetXaxis()->SetBinLabel(i, "AdaBoostBeta");

   i = 20; hSetup->SetBinContent(i, fBdtSetup.MaxDepth); hSetup->GetXaxis()->SetBinLabel(i, "MaxDepth");
   i = 21; hSetup->SetBinContent(i, fBdtSetup.NNodesMax); hSetup->GetXaxis()->SetBinLabel(i, "NNodesMax");

   cout << "----------------------------------------------------------------------" << endl;
   cout << "==> oname: " << oname << " antimuon: " << endl;

   string optstring = "V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
   optstring        = "V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
   cout << "==> Factory: " << optstring << endl;
   TMVA::Factory *factory = new TMVA::Factory(Form("%s", oname.c_str()), outputFile,  optstring.c_str());

   // -- parse string with all variables into a vector
   vector<string> vVar; 
   if (fVariables == "all") {
     fVariables = "pt:eta";
     fVariables += ":intnmatchedstations:intvalidpixelhits:inttrklayerswithhits";
     fVariables += ":gchi2:itrkvalidfraction:segcomp:chi2lmom:chi2lpos:gtrkprob:ntrkvhits:ntrkehitsout";
     fVariables += ":kink";
     fVariables += ":dptrel:deta:dphi";
   }
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
     if (string::npos != vVar[i].find("int")) {
       factory->AddVariable(vVar[i].c_str(), 'I');        
     } else {
       factory->AddVariable(vVar[i].c_str(), 'F');        
     }
   }

   writeOut(outputFile, hSetup); 

   Double_t signalWeight      = 1.; //= LUMISCALE; // 0.000388
   Double_t backgroundWeight = 1.;

   cout << "--> signal weight:     " << signalWeight << endl;
   cout << "--> background weight: " << backgroundWeight << endl;

   factory->AddTree(trainSg,     "Signal",     signalWeight,  "", "train");
   factory->AddTree(testSg,      "Signal",     signalWeight,  "", "test");
   factory->AddTree(trainBg, "Background", backgroundWeight, "", "train");
   factory->AddTree(testBg,  "Background", backgroundWeight, "", "test");

   int nSgTrain = trainSg->GetEntries();
   int nSgTest  = testSg->GetEntries();

   //    nSgTrain = 20000;
   //    nSgTest  = 20000;

   int nBgTrain = trainBg->GetEntries();
   int nBgTest  = testBg->GetEntries();

   nSgTrain = 0; 
   nSgTest  = 0; 

   nBgTrain = 0; 
   nBgTest  = 0; 

   optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V", 
		    nSgTrain, nSgTest, nBgTrain, nBgTest); 
   cout << "==> PrepareTrainingAndTestTree: " << optstring << endl;
   factory->PrepareTrainingAndTestTree("", "", optstring.c_str());
   
   optstring = "!H:V" + fBDTParameters;

   cout << "==> BookMethod: " << optstring << endl;
   factory->BookMethod( TMVA::Types::kBDT, "BDT", optstring.c_str());
   //   factory->BookMethod( TMVA::Types::kMLP, "MLP_ANN");

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

   mvas(oname); 

   delete factory;

}



// ----------------------------------------------------------------------
void tmva10::mvas(string fname) { //, HistType htype, Bool_t useTMVAStyle ) {
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
TH1D* tmva10::getRanking(string fname, string prefix, std::string type) {
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
    after = "Train on events1, test on events2";
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
void tmva10::createInputFile(string ofile, string sgfile, string bgfile) {
  TFile *sinput = TFile::Open(sgfile.c_str());
  TFile *binput = TFile::Open(bgfile.c_str());

  cout << "==> signal input file:     " << sinput->GetName() << std::endl;
  cout << "==> background input file: " << binput->GetName() << std::endl;
  
  TTree *signal      = (TTree*)sinput->Get("candAnaMuMu/muonidtree");
  TTree *cbackground = (TTree*)binput->Get("candAnaMuMu/muonidtree");
  
  TFile *outFile = TFile::Open(ofile.c_str(),"RECREATE");

  string sdir, type; 
  TTree *copyTree(0);
  TCut copyCuts; 
  TCut typeCut;
  for (int j = 0; j < 2; ++j) {
    type = Form("Events%d", j); 
    typeCut = Form("%s==%d", "2*rndm%2", j);
    
    // -- signal
    sdir = Form("sg%d", j); 
    outFile->mkdir(sdir.c_str());
    outFile->cd(sdir.c_str());
    copyCuts = typeCut;
    cout << "sg copyCuts: " << copyCuts << endl;
    copyTree = signal->CopyTree(copyCuts);
    cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;
    
    // -- background
    sdir = Form("bg%d", j); 
    outFile->mkdir(sdir.c_str());
    outFile->cd(sdir.c_str());
    copyCuts = typeCut;
    cout << "bg copyCuts: " << copyCuts << endl;
    copyTree = cbackground->CopyTree(copyCuts);
    cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;
  
  }

  outFile->Write();
  outFile->Close();   
  
  sinput->Close();
  binput->Close();
}

// ----------------------------------------------------------------------
void tmva10::writeOut(TFile *f, TH1 *h) {
  TDirectory *pD = gDirectory; 
  f->cd();
  h->SetDirectory(f); 
  h->Write();
  pD->cd();
}


// ----------------------------------------------------------------------
void tmva10::redrawStats(double x, double y, const char *newname, int color) {

  TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
  st->SetName(newname);
  st->SetX1NDC(x); 
  st->SetY1NDC(y); 
  st->SetTextColor(color); 
  //  st->Draw("sames");
}

// ----------------------------------------------------------------------
void tmva10::newLegend(double x1, double y1, double x2, double y2, string title) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}

