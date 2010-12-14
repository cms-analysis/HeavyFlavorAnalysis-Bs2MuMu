#include "anaBmm.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "AnalysisDistribution.hh"

#include "TF1.h"
#include "TKey.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TCanvas.h"

#include <iomanip>
#include <string>

// Usage: root[0] anaBmm(

using namespace std; 
using std::string; 

ClassImp(anaBmm)

// ----------------------------------------------------------------------
anaBmm::anaBmm(const char *files, const char *dir, int mode) { 
  init(files, dir, mode);
}

// ----------------------------------------------------------------------
void anaBmm::init(const char *files, const char *dir, int mode) {

  fNData = fNMc = 0; 
  fSgData = fSgMc = fNoData = fNoMc = 0; 

  fFont = 42; 
  fMode = mode;  

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(0.07);
  tl->SetTextFont(fFont);

  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  f0 = new TF1("f0", f_p1, 0., 6., 2); 

  f1 = new TF1("f1", f_expo, 0., 6.0, 2);

  f2 = new TF1("f2", f_p1aG, 0., 6., 5);
  f2->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f3 = new TF1("f3", f_eaG, 0., 6., 5);
  f3->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f4 = new TF1("f4", f_ea2G, 0., 6., 8);
  f4->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 

  fDirectory = dir; 
  fSuffix    = dir; 
  if (fMode > 0)  fSuffix += Form("%d", fMode); 
  cout << "--> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  loadFiles(files);

}


// ----------------------------------------------------------------------
void anaBmm::loadFiles(const char *files) {

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
    
    if (string::npos != sdset.find("data")) {
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) fSgData = fNData;
      if (string::npos != stype.find("default") && string::npos != stype.find("no")) fNoData = fNData;
      fpData[fNData] = loadFile(sfile, stype); 
      cout << "open data " << sfile << " as " << stype << endl;
      ++fNData;
    } else {
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) fSgMc = fNMc;
      if (string::npos != stype.find("default") && string::npos != stype.find("no")) fNoMc = fNMc;
      fpMc[fNMc] = loadFile(sfile, stype); 
      cout << "open MC " << sfile << " as " << stype << endl;
      ++fNMc;
    }
  }
}


// ----------------------------------------------------------------------
TFile* anaBmm::loadFile(string file, string type) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}


// ----------------------------------------------------------------------
void anaBmm::makeAll(int channel) {
  fpMc[fSgMc]->cd();
  dumpCutNames();

//   fpData[fNoData]->cd(); 
//   effTable("NoData");
//   return; 

  // -- BMM
  fpMc[fSgMc]->cd(); 
  effTable("SgMc");
  fpData[fSgData]->cd(); 
  effTable("SgData");

  // -- B+
  fpData[fNoData]->cd(); 
  effTable("NoData");
  fpMc[fNoMc]->cd(); 
  effTable("NoMc");

}


// ----------------------------------------------------------------------
void anaBmm::dumpCutNames() {

  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  
  TH1D *h = (TH1D*)gFile->Get("hcuts"); 
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
    OUT <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
    OUT <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
  }
}


// ----------------------------------------------------------------------
void anaBmm::effTable(string smode) {
  TH1D *h = (TH1D*)gFile->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }
  
  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << " % --- " << smode << endl;
  int mode(0); 
  if (string::npos != smode.find("SgMc")) {
    cout << "==> This is signal MC, counting events in histogram" << endl;
    mode = 0; 
  }
  if (string::npos != smode.find("SgData")) {
    cout << "==> This is signal Data, counting events in blind signal box from sidebands" << endl;
    mode = 1; 
  }
  
  if (string::npos != smode.find("No")) {
    if (0 == fMode) {
      cout << "==> This is normalization, fitting mode = 11 = expo+Gaus" << endl;
      mode = 11; 
    } else {
      cout << "==> This is normalization, fitting mode = " << fMode << endl;
      mode = fMode; 
    }
  }
  
  string cut, pdfname;
  double n(0.), nE(0.), nprev(0.), nprevE(0.); 
  double norm(0.),  normE(0.), relEff(0.), relEffE(0.), cumEff(0.), cumEffE(0.);
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 
    
    if (cut == string("")) {
      cout << "empty string found, break ..." << endl;
      break;
    }
    pdfname = Form("%s/%s_%s.pdf", fDirectory.c_str(), smode.c_str(), cut.c_str());
    cout << "AD for " << cut << " results in " << pdfname << endl;
    AnalysisDistribution a(cut.c_str());
    a.fMassLo = 5.0; 
    a.fMassHi = 5.5; 
    a.hMassCu->SetMinimum(0.);
    n = a.fitMass(a.hMassCu, nE, mode); 
    
    if ((string::npos != cut.find("tracks")) || (string::npos != cut.find("muons"))) {
      n *=0.5; 
      nE *=0.5; 
    }

    if ("hlt" == cut) {
      cout << "initialize normalization numbers" << endl;
      norm  = n; 
      normE = nE; 
      nprev = norm; 
      nprevE= normE; 
    }

    relEff = n/nprev;
    relEffE= dEff(n, nE, nprev, nprevE);
    cumEff = n/norm;
    cumEffE= dEff(n, nE, norm, normE);
    
    cout << cut << " n = " << n << "+/-" << nE 
	 << "rel eff = " << relEff << "+/-" << relEffE
	 << "cum eff = " << cumEff << "+/-" << cumEffE
	 << endl;
    pdfname = Form("%s/%s_%s_%s_hMassCu.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str(), cut.c_str());
    c0->SaveAs(pdfname.c_str(), "Portrait");

    OUT << Form("%s", (formatTex(n, fSuffix+":"+cut+"N", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(nE, fSuffix+":"+cut+"NE", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(relEff, fSuffix+":"+cut+"eRel", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(relEffE, fSuffix+":"+cut+"eRelE", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(cumEff, fSuffix+":"+cut+"eCum", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(cumEffE, fSuffix+":"+cut+"eCumE", smode)).c_str()) << endl;

    nprev = n; 
    nprevE= nE; 
  }
}




// ----------------------------------------------------------------------
string anaBmm::formatTex(double n, string name, string tag) {
  
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
void anaBmm::makeCanvas(int i) {
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
int anaBmm::wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
}


// ----------------------------------------------------------------------
void anaBmm::replaceAll(std::string &s, std::string a, std::string b) {
  
  TString ts(s.c_str()); 
  ts.ReplaceAll(a.c_str(), b.c_str()); 
  s = ts.Data(); 

}
