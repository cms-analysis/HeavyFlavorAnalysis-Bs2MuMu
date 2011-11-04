#include "plotClass.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "../macros/AnalysisDistribution.hh"

#include "TF1.h"
#include "THStack.h"
#include "TKey.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

#include <iomanip>
#include <string>
#include <list>

using namespace std; 
using std::string; 

ClassImp(plotClass)

// ----------------------------------------------------------------------
double fa_expo(double *x, double *par) {
  return par[0]*TMath::Exp(x[0]*par[1]);
} 

// ----------------------------------------------------------------------
double fa_err(double *x, double *par) {
  // from DK: TMath::Erf((a1-x)/a2))+a3
  return par[3]*(TMath::Erf((par[0]-x[0])/par[1])+par[2]); 
} 

// ----------------------------------------------------------------------
double fa_pol1(double *x, double *par) {
  return par[0] + par[1]*x[0]; 
}

// ----------------------------------------------------------------------
// expo and err and gauss 
double fa_expo_err(double *x, double *par) {
  //   par[0] = norm
  //   par[1] = exp
  //   par[2] = par[0] of err
  //   par[3] = par[1] of err
  //   par[4] = par[2] of err
  //   par[5] = par[3] of err
  return  (fa_err(x, &par[2]) + fa_expo(x, &par[0]));
}

// ----------------------------------------------------------------------
// expo and err and gauss 
double fa_pol1_err(double *x, double *par) {
  //   par[0] = const 
  //   par[1] = slope
  //   par[2] = par[0] of err
  //   par[3] = par[1] of err
  //   par[4] = par[2] of err
  //   par[5] = par[3] of err
  return  (fa_err(x, &par[2]) + fa_pol1(x, &par[0]));
}


// ----------------------------------------------------------------------
plotClass::plotClass(const char *files, const char *cuts, const char *dir, int mode) { 

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  init(files, cuts, dir, mode);
}

// ----------------------------------------------------------------------
plotClass::~plotClass() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotClass::init(const char *files, const char *cuts, const char *dir, int mode) {

  gStyle->SetHatchesSpacing(2);

  fDoPrint = true; // create output

  legg = 0; 

  fSize = 0.05; 
  
  fpFunc  = new initFunc(); 

  fMassLo = 4.5; 
  fMassHi = 6.5;

  fSgLo = 5.27;
  fSgHi = 5.47;

  fNoLo = 5.10;
  fNoHi = 5.40;

  fCsLo = 5.27;
  fCsHi = 5.47;

  fBgLo = 4.9;
  fBgHi = 5.9;

  // -- initialize cuts
  cout << "Reading cuts from " << Form("plotClass.%s.cuts", cuts) << endl;
  readCuts(Form("plotClass.%s.cuts", cuts)); 

  printCuts(cout); 

  fFont = 42; 
  fMode = mode;  

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(fSize);
  tl->SetTextFont(fFont);

  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  f0 = new TF1("f0", f_p1, 0., 7., 2); 

  f1 = new TF1("f1", f_expo, 0., 7.0, 2);

  f2 = new TF1("f2", f_p1aG, 0., 7., 5);
  f2->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f3 = new TF1("f3", f_eaG, 0., 7., 5);
  f3->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f4 = new TF1("f4", f_ea2G, 0., 7., 8);
  f4->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 

  fDirectory = dir; 
  fSuffix    = cuts; 
  if (fMode > 0)  fSuffix += Form("-%d", fMode); 
  cout << "--> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".txt";
  // FIXME
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fOUT.open(fNumbersFileName.c_str(), ios::app);
  // FIXME
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);
  
  loadFiles(files);
  string hfname  = fDirectory + "/anaBmm." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  printCuts(fOUT); 

  dumpSamples();
  fF["SgMc"]->cd();
  dumpCutNames("candAnaMuMu/hcuts");

}


// ----------------------------------------------------------------------
void plotClass::loadFiles(const char *files) {

  char buffer[1000];
  ifstream is(files);
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    
    string sbuffer = string(buffer); 

    string::size_type m1 = sbuffer.find("dset="); 
    string stype = sbuffer.substr(5, m1-6); 

    string::size_type m2 = sbuffer.find("lumi="); 
    string sdset = sbuffer.substr(m1+5, m2-m1-6); 
    m1 = m2; 

    m2 = sbuffer.find("file="); 
    string slumi = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 
    string sname; 
    
    // -- DATA
    TFile *pF(0); 
    if (string::npos != sdset.find("data")) {
      pF = loadFile(sfile); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgData"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	sname = "NoData";
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsData"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }
    } else {
      string sfilter = sdset; 
      replaceAll(sfilter, "mc,", ""); 
      double effFilter = atof(sfilter.c_str());
      // -- MC
      pF = loadFile(sfile); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgMc"; 
	fF.insert(make_pair(sname, pF)); 	
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("sg")) {
	sname = "SgMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC 2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("sg")) {
	sname = "SgMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC 3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	sname = "BdMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("bd")) {
	sname = "BdMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC 2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("bd")) {
	sname = "BdMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC 3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	sname = "NoMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("no")) {
	sname = "NoMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (MC 2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("no")) {
	sname = "NoMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (MC 3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("2e33") && string::npos != stype.find("cs")) {
	sname = "CsMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (MC 2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("3e33") && string::npos != stype.find("cs")) {
	sname = "CsMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (MC 3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,82")) {
	sname = "bg82"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,83")) {
	sname = "bg83"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,84")) {
	sname = "bg84"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,86")) {
	sname = "bg86"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{-}#mu^{+}#nu")); 
      }	

      if (string::npos != stype.find("bg,95")) {
	sname = "bg95"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{-}#mu^{+}#nu")); 
      }	
      if (string::npos != stype.find("bg,93")) {
	sname = "bg93"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,92")) {
	sname = "bg92"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,91")) {
	sname = "bg91"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,60")) {
	sname = "bg60"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p K^{-}")); 
      }	
      if (string::npos != stype.find("bg,61")) {
	sname = "bg61"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p #pi^{-}")); 
      }	
      cout << "open MC file "  << sfile  << " as " << sname << " (" << stype << ") with lumi = " << slumi << endl;
    }
  }


}


// ----------------------------------------------------------------------
TFile* plotClass::loadFile(string file) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}


// ----------------------------------------------------------------------
void plotClass::makeAll(int channel) {


}


// ----------------------------------------------------------------------
void plotClass::dumpSamples() {
  // FIXME
  std::map<string, double> ngen;
  ngen.insert(make_pair("bg60", 20.e6)); 
  ngen.insert(make_pair("bg61", 20.e6)); 

  ngen.insert(make_pair("bg86", 20.e6)); 
  ngen.insert(make_pair("bg84", 20.e6)); 
  ngen.insert(make_pair("bg83", 20.e6)); 
  ngen.insert(make_pair("bg82", 20.e6)); 

  ngen.insert(make_pair("bg91", 20.e6)); 
  ngen.insert(make_pair("bg92", 20.e6)); 
  ngen.insert(make_pair("bg93", 20.e6)); 
  ngen.insert(make_pair("bg95", 16.8e6)); 

  ngen.insert(make_pair("SgMcAcc", 100.e6)); 

  ngen.insert(make_pair("BdMcAcc", 40.e6)); 

  ngen.insert(make_pair("NoMcAcc", 160.e6)); 
  //  ngen.insert(make_pair("CsMcAcc", 40.e6)); 
  ngen.insert(make_pair("CsMcAcc", 154.e6)); 
  
  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);

  fTEX << "% ----------------------------------------------------------------------" << endl;
  string name; 
  double lumi(0), n(0), f(0); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;
    name = imap->second; 
    lumi = fLumi[imap->first];
    n = ngen[imap->first];
    //FIXME    f = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1);
    replaceAll(name, "#", "\\"); 
    //    cout <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    fTEX <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    //    cout <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    fTEX <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    if (n>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
    if (f>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), f) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
  }
  fTEX.flush();

}



// ----------------------------------------------------------------------
void plotClass::dumpCutNames(const char *hname) {

  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  
  TH1D *h = (TH1D*)gFile->Get(hname); 
  string cut, empty(""), cutstring, cutline, cutvalue; 
  double value; 
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i));
    value = h->GetBinContent(i);
    if (cut == empty) continue;
    string::size_type n1 = cut.find_first_of(":"); 
    string::size_type n2 = cut.find_last_of(":")-1; 
    
    cutstring = cut.substr(0, n1-1); 
    cutline   = cut.substr(n1+3, n2-n1-4);
    replaceAll(cutline, "#", "\\"); 
    cutvalue  = cut.substr(n2+2, cut.size()-n2);

    //    cout << cut << " cutstring==>" << cutstring << "<== cutvalue==>" << cutvalue << " cutline==>" << cutline << "<==" << endl;
    if (string::npos != cutstring.find("JSON")) {
      replaceAll(cutvalue, "_", "\\_"); 
      replaceAll(cutvalue, "cuts/", ""); 
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    } else if (string::npos != cutstring.find("CANDCOSALPHA")) {
      fTEX <<  Form("\\vdef{%s:CANDALPHA:cutLine}   {\\ensuremath{{\\alpha } } }", fSuffix.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:CANDALPHA:cutValue}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), TMath::ACos(atof(cutvalue.c_str()))) 
	   << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    }

  }
  fTEX.flush();
}



// ----------------------------------------------------------------------
string plotClass::scientificTex(double n, double nE, std::string name, double base, int digits) {

  char line[200]; 
  double a1 = n/base; 
  double a2 = nE/base; 
  int  expo = TMath::Log10(base); 

  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name.c_str(), static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%2.1f \\pm %2.1f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%3.2f \\pm %3.2f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%4.3f \\pm %4.3f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.4f \\pm %5.4f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%6.5f \\pm %6.5f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%f \\pm %f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  }

  return string(line); 
}


// ----------------------------------------------------------------------
string plotClass::formatTex(double n, std::string name, int digits) {

  char line[200]; 
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name.c_str(), static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.1f } } }", name.c_str(), n);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.2f } } }", name.c_str(), n);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.3f } } }", name.c_str(), n);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.4f } } }", name.c_str(), n);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%6.5f } } }", name.c_str(), n);
  } else if (6 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%7.6f } } }", name.c_str(), n);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%f } } }", name.c_str(), n);
  }

  return string(line); 
}



// ----------------------------------------------------------------------
string plotClass::formatTex(double n, string name, string tag) {
  
  char line[200]; 
  
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str(), tag.c_str());
    //   } else if ( n > 1.e10) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
    //   } else if ( n > 1.e4) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
  } else if ( n > 100. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%6.0f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.2f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-1) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-3) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
    //   } else if ( n > 1.e-9 ) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
    //   } else if ( n > 1.e-19 ){
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
  } else {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{0.0 } } }", name.c_str(), tag.c_str());
  }
  
  string result(line); 
  return result;
}



// ----------------------------------------------------------------------
void plotClass::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}



// ----------------------------------------------------------------------
void plotClass::replaceAll(std::string &s, std::string a, std::string b) {
  
  TString ts(s.c_str()); 
  ts.ReplaceAll(a.c_str(), b.c_str()); 
  s = ts.Data(); 

}


// ----------------------------------------------------------------------
void plotClass::setErrors(TH1D *h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinLowEdge(i) > 5.1 && h->GetBinLowEdge(i+1) <= 5.5) continue;
    if (h->GetBinContent(i) < 1) {
      h->SetBinError(i, 1.); 
    }
  }
}


// ----------------------------------------------------------------------
void plotClass::newLegend(double x1, double y1, double x2, double y2) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2);
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}


// ----------------------------------------------------------------------
void plotClass::readCuts(const char *filename) {
  cout << "==> plotClass: Reading " << filename << " for cut settings" << endl;
  vector<string> cutLines; 
  char  buffer[200];
  ifstream is(filename);
  while (is.getline(buffer, 200, '\n')) {
    cutLines.push_back(string(buffer));
  }

  char CutName[100];
  float CutValue;
  int dump(0), ok(0);

  cuts *a = 0;
  
  fCuts.clear();

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "index")) {
      ok = 1;
      if (dump) cout << "index:            " << CutValue << endl;
      if (a) fCuts.push_back(a); 
      a = new cuts; 
      a->index = static_cast<int>(CutValue); 
    }
    
    if (!strcmp(CutName, "mBdLo")) {
      a->mBdLo = CutValue; ok = 1;
      if (dump) cout << "mBdLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBdHi")) {
      a->mBdHi = CutValue; ok = 1;
      if (dump) cout << "mBdHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsLo")) {
      a->mBsLo = CutValue; ok = 1;
      if (dump) cout << "mBsLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsHi")) {
      a->mBsHi = CutValue; ok = 1;
      if (dump) cout << "mBsHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMin")) {
      a->etaMin = CutValue; ok = 1;
      if (dump) cout << "etaMin:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMax")) {
      a->etaMax = CutValue; ok = 1;
      if (dump) cout << "etaMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "pt")) {
      a->pt = CutValue; ok = 1;
      if (dump) cout << "pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1pt")) {
      a->m1pt = CutValue; ok = 1;
      if (dump) cout << "m1pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2pt")) {
      a->m2pt = CutValue; ok = 1;
      if (dump) cout << "m2pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1eta")) {
      a->m1eta = CutValue; ok = 1;
      if (dump) cout << "m1eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2eta")) {
      a->m2eta = CutValue; ok = 1;
      if (dump) cout << "m2eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "iso")) {
      a->iso = CutValue; ok = 1;
      if (dump) cout << "iso:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "chi2dof")) {
      a->chi2dof = CutValue; ok = 1;
      if (dump) cout << "chi2dof:             " << CutValue << endl;
    }

    if (!strcmp(CutName, "alpha")) {
      a->alpha = CutValue; ok = 1;
      if (dump) cout << "alpha:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "fls3d")) {
      a->fls3d = CutValue; ok = 1;
      if (dump) cout << "fls3d:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "docatrk")) {
      a->docatrk = CutValue; ok = 1;
      if (dump) cout << "docatrk:               " << CutValue << endl;
    }

  }

  if (a) fCuts.push_back(a); 

  if (!ok) cout << "==> what about " << CutName << endl;
  

}


// ----------------------------------------------------------------------
void plotClass::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    cuts *a = fCuts[i]; 
    OUT << "# -- channel " << a->index << endl;
    OUT << "index   " << a->index << endl;
    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% -- Cuts for channel " << a->index << endl;

    OUT << "mBdLo   " << Form("%4.3f", a->mBdLo) << endl;
    OUT << "mBdHi   " << Form("%4.3f", a->mBdHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBdLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBdHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdHi) << endl;

    OUT << "mBsLo   " << Form("%4.3f", a->mBsLo) << endl;
    OUT << "mBsHi   " << Form("%4.3f", a->mBsHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBsLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBsHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsHi) << endl;

    OUT << "etaMin  " << Form("%3.1f", a->etaMin) << endl;
    OUT << "etaMax  " << Form("%3.1f", a->etaMax) << endl;
    fTEX <<  Form("\\vdef{%s:etaMin:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMin) << endl;
    fTEX <<  Form("\\vdef{%s:etaMax:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMax) << endl;

    OUT << "pt      " << Form("%3.1f", a->pt) << endl;
    fTEX <<  Form("\\vdef{%s:pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->pt) << endl;
    OUT << "m1pt    " << a->m1pt << endl;
    fTEX <<  Form("\\vdef{%s:m1pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1pt) << endl;
    OUT << "m2pt    " << a->m2pt << endl;
    fTEX <<  Form("\\vdef{%s:m2pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2pt) << endl;
    OUT << "m1eta   " << a->m1eta << endl;
    fTEX <<  Form("\\vdef{%s:m1eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1eta) << endl;
    OUT << "m2eta   " << a->m2eta << endl;
    fTEX <<  Form("\\vdef{%s:m2eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2eta) << endl;

    OUT << "iso    " << a->iso << endl;
    fTEX <<  Form("\\vdef{%s:iso:%d}   {\\ensuremath{{%3.2f } } }", fSuffix.c_str(), a->index, a->iso) << endl;
    OUT << "chi2dof " << a->chi2dof << endl;
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->chi2dof) << endl;
    OUT << "alpha   " << a->alpha << endl;
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->alpha) << endl;
    OUT << "fls3d   " << a->fls3d << endl;
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->fls3d) << endl;
    OUT << "docatrk   " << a->docatrk << endl;
    fTEX <<  Form("\\vdef{%s:docatrk:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->docatrk) << endl;
  }
  OUT.flush();
}



// ----------------------------------------------------------------------
void plotClass::stamp(double x1, string text1, double x2, string text2) {
  tl->SetTextSize(fSize); 
  tl->DrawLatex(x1, 0.91, text1.c_str());   
  tl->DrawLatex(x2, 0.91, text2.c_str()); 
  //  tl->DrawLatex(x2, 0.85, text2.c_str()); 
}



// ----------------------------------------------------------------------
void plotClass::drawArrow(double height, int mode, int color) {

  double ylo(0.01); 
  pl->SetLineWidth(3.); 
  
  double d(0.08), y(0.80), x(5.25); 
  
  if (1 == mode) {
    pl->SetLineColor(kBlue); 
    pl->SetLineColor(kBlue); 
    pl->SetLineStyle(kSolid); 
    pl->DrawLine(fCuts[0]->mBsLo, height, fCuts[0]->mBsHi, height); 
    pl->SetLineWidth(2.); 
    pl->DrawLine(fCuts[0]->mBsLo, height+d, fCuts[0]->mBsLo, height-d); 
    pl->DrawLine(fCuts[0]->mBsHi, height+d, fCuts[0]->mBsHi, height-d); 

    if (1) {
      y = 1.75;
      x = 5.25; 
      d = 0.05;
      pl->SetLineWidth(3.); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineWidth(2.); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-0.05, "B_{s}^{0} signal window");
    }

  } else if (2 == mode) {
    pl->SetLineColor(kRed); 
    pl->SetLineColor(kRed); 
    pl->SetLineStyle(kDashed); 
    pl->DrawLine(fCuts[0]->mBdLo, height, fCuts[0]->mBdHi, height); 
    pl->SetLineStyle(kSolid); 
    pl->SetLineWidth(2.); 
    pl->DrawLine(fCuts[0]->mBdLo, height+d, fCuts[0]->mBdLo, height-d); 
    pl->DrawLine(fCuts[0]->mBdHi, height+d, fCuts[0]->mBdHi, height-d); 

    if (1) {
      x = 5.25; 
      y = 1.55;
      d = 0.05;
      pl->SetLineWidth(3.); 
      pl->SetLineStyle(kDashed); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineStyle(kSolid); 
      pl->SetLineWidth(2.); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-0.05, "B^{0} signal window");
    }

  } else if (3 == mode) {
    pa->SetLineColor(kBlack); 
    pa->SetFillColor(kBlack); 
    pa->DrawArrow(4.9, 0.5, 4.9, ylo); 
    pa->DrawArrow(5.2, 0.5, 5.2, ylo); 
    pl->SetLineColor(kBlack); 
    pl->DrawLine(4.9, 0.5, 5.2, 0.5);

    pa->DrawArrow(5.45, 0.5, 5.45, ylo); 
    pa->DrawArrow(5.90, 0.5, 5.90, ylo); 
    pl->SetLineColor(kBlack); 
    pl->DrawLine(5.45, 0.5, 5.90, 0.5);
  }

  
 
}


// ----------------------------------------------------------------------
void plotClass::drawBox(int mode, double hi, int ylo) {

  TBox *b = new TBox; 

  if (1 == mode) {
    b->SetFillColor(kBlue); 
    b->SetLineColor(kBlue); 
    b->SetFillStyle(3004); 
    b->SetFillStyle(3356); 
    b->DrawBox(fCuts[0]->mBsLo, ylo, fCuts[0]->mBsHi, hi); 
  } else if (2 == mode) {
    b->SetLineColor(kRed); 
    b->SetFillColor(kRed); 
    b->SetFillStyle(3005); 
    b->SetFillStyle(3365); 
    b->DrawBox(fCuts[0]->mBdLo, ylo, fCuts[0]->mBdHi, hi); 
  } else if (3 == mode) {
    b->SetLineColor(kBlack); 
    b->SetFillColor(kBlack); 
    b->DrawBox(4.9, ylo, 5.2, hi); 
    b->DrawBox(5.45, ylo, 5.90, hi); 
  }

  
 
}
