{

  string version = gSystem->Getenv("VERSION");
  cout << "=> Local rootlogon.C <=" << endl;
  cout << "Loading libPhysics.so" << endl;
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");

  cout << "Loading libMinuit.so" << endl;
  gSystem->Load("libMinuit.so");

  cout << "Loading libTMVA.so" << endl;
  gSystem->Load("libTMVA.so");

  cout << "Loading libGpad.so" << endl;
  gSystem->Load("libGpad.so");

  cout << "Loading lib/libTmvaClasses.so" << endl;
  gSystem->Load("lib/libTmvaClasses.so");

  gSystem->Load("libRooFit.so");
  using namespace RooFit;

  //  gROOT->Macro("cms-tdr.C");
  //  gROOT->ForceStyle();

  //  gStyle->SetTitleBorderSize(0);  // no border around histogram title (font size can't be calculated anyways ...)
  //  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetStatStyle(0);      // for a completely transparent stat box
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 

  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextFont(42); 

  TLine *pl = new TLine();  

  // --- Cleanup if this is not the first call to rootlogon.C
  TCanvas *c = 0;
  c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
  p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
  // --- Create a new canvas.
  //  TCanvas c0("c0","--c0--",815,0,656,700);
  TCanvas c0("c0","--c0--",680,0,600,700);
  c0->ToggleEventStatus();
}



