{

  cout << "Loading libPhysics.so" << endl;
  gSystem->Load("libPhysics.so");

  cout << "Loading libAna00.so" << endl;
  gSystem->Load("../rootio/lib/libAna00.so");

  cout << "Loading libAnaClasses.so" << endl;
  gSystem->Load("../rootio/lib/libAnaClasses.so");

  gROOT->SetStyle("Plain");

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetStatStyle(0);      // for a completely transparent stat box
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.8);
  gStyle->SetMarkerColor(1);

  gStyle->SetNdivisions(505, "x");
  gStyle->SetTickLength(-0.02, "x");
  gStyle->SetTickLength(-0.02, "y");
  gStyle->SetLabelOffset(0.02, "x");
  gStyle->SetLabelOffset(0.02, "y");
  gStyle->SetPadLeftMargin(0.15); 
  gStyle->SetPadBottomMargin(015); 

  gStyle->SetTitleBorderSize(0);  // no border around histogram title (font size can't be calculated anyways ...)
  gStyle->SetStatFont(132); 
  gStyle->SetTextFont(132); 
  gStyle->SetLabelFont(132, "X"); 
  gStyle->SetLabelFont(132, "Y"); 
  gStyle->SetTitleFont(132); 

//  gROOT->ForceStyle();

  TLatex *tl = new TLatex();
  TLine *pl = new TLine();  

  // --- Cleanup if this is not the first call to rootlogon.C
  TCanvas *c = 0;
  c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
  p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
  // --- Create a new canvas.
  //  TCanvas c0("c0","--c0--",356,0,656,700);
  TCanvas c0("c0","--c0--",615,0,656,700);
  c0->ToggleEventStatus();
  c0->SetBottomMargin(0.15);
}



