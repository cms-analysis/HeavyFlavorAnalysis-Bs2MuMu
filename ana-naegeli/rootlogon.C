{
  // Load my Personal utility functions
  gSystem->AddDynamicPath("/Users/cn/Documents/PSI/Projects/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/rootutils/lib");
  cout << "Loading utility functions..." << flush;
  gSystem->Load("libNCRootUtils.so");
  cout << "ok" << endl;

  gSystem->AddDynamicPath("/Users/cn/Documents/PSI/Projects/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/ana-naegeli/lib");
  cout << "Loading NCAna library..." << flush;
  gSystem->Load("libNCAna.so");
  cout << "ok" << endl;

  cout << "Setting CMS histogram default style" << endl;
  gROOT->Macro("/Users/cn/Documents/PSI/Projects/macros/cms-tdr.C");
  gROOT->ForceStyle();

  gStyle->SetPalette(1);
}
