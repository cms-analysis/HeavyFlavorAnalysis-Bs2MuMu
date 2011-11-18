#include "plotIso.hh"

#include <algorithm>

#include "../macros/AnalysisDistribution.hh"
#include "TMath.h"

using namespace std; 
using std::string; 

ClassImp(plotIso)

// ----------------------------------------------------------------------
plotIso::plotIso(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

}


// ----------------------------------------------------------------------
plotIso::~plotIso() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotIso::makeAll(double isoCut) {
  
  fIsoCut = isoCut; 
  fMode = 3; 
  fPreco = 5.1;
  dataVsMc("NoData", "candAnaBu2JpsiK",   "NoMc", "candAnaBu2JpsiK",   "Ao"); 

  fMode = 2; 
  dataVsMc("CsData", "candAnaBs2JpsiPhi", "CsMc", "candAnaBs2JpsiPhi", "Ao"); 

  fMode = 0; 
  dataVsMc("SgData", "candAnaMuMu", "SgMc", "candAnaMuMu", "Ao"); 

  anaTextFiles(static_cast<int>(100*isoCut)); 

}

// ----------------------------------------------------------------------
bool largerRatio(const data a, const data b) {
  return (a.ratio > b.ratio); 
}

// ----------------------------------------------------------------------
bool closerToOne(const data a, const data b) {
  return (TMath::Abs(a.ratio - 1) < TMath::Abs(b.ratio - 1)); 
}

// ----------------------------------------------------------------------
void plotIso::anaTextFiles(int cutval) {

  vector<data> mm, no, cs; 
  
  readFile(Form("plotIso-%d-SgData-SgMc.txt", cutval), mm);
  readFile(Form("plotIso-%d-NoData-NoMc.txt", cutval), no);
  readFile(Form("plotIso-%d-CsData-CsMc.txt", cutval), cs);

  sort(mm.begin(), mm.end(), largerRatio); 
  double maxDev(1.03); 

  for (unsigned int i = 0; i < mm.size(); ++i) {
    for (unsigned int j = 0; j < no.size(); ++j) {
      for (unsigned int k = 0; k < cs.size(); ++k) {
	if (no[j].name == mm[i].name && (cs[k].name == mm[i].name)) {
	  if (no[j].ratio < maxDev &&cs[k].ratio < maxDev) 
	    cout << mm[i].name << ": " << mm[i].ratio << " r = " << mm[i].r << " pt = " << mm[i].pt << " doca = " << mm[i].doca 
		 << " and " << no[j].name << " ratio = " << no[j].ratio 
		 << " and " << cs[k].name << " ratio = " << cs[k].ratio 
		 << endl;
	}
      }
    }
  }

  for (unsigned int i = 0; i < mm.size(); ++i) {
    for (unsigned int j = 0; j < no.size(); ++j) {
      for (unsigned int k = 0; k < cs.size(); ++k) {
	if ((mm[i].name == "Doca05isor010pt09") && (no[j].name == mm[i].name) && (cs[k].name == mm[i].name)) {
	  cout << "* " << mm[i].name << ": " << mm[i].ratio << " r = " << mm[i].r << " pt = " << mm[i].pt << " doca = " << mm[i].doca 
	       << " and " << no[j].name << " ratio = " << no[j].ratio 
	       << " and " << cs[k].name << " ratio = " << cs[k].ratio 
	       << endl;
	}
      }
    }
  }
}


// ----------------------------------------------------------------------
void plotIso::readFile(string fname, vector<data> &v) {

  int i0, i1, i2; 
  float f0, f1, f2; 
  char buffer[1000]; 
  ifstream is(fname.c_str()); 
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "Doca%disor%dpt%d: %f %f %f", &i0, &i1, &i2, &f0, &f1, &f2); 
    //    cout << "doca: " << i0*0.01 << " r = " << i1*0.1 << " pt = " << i2*0.1 << endl;
    data a; 
    a.name = Form("Doca0%disor0%dpt0%d", i0, i1, i2);
    a.doca = i0*0.01; 
    a.r    = i1*0.1; 
    a.pt   = i2*0.1;
    a.eff1 = f0; 
    a.eff2 = f1; 
    a.ratio= f2; 
    v.push_back(a); 
  }
  is.close(); 

}


// ----------------------------------------------------------------------
void plotIso::dataVsMc(string file1, string dir1, string file2, string dir2, string selection) {

  double cutval = fIsoCut; 
  ofstream OUT(Form("plotIso-%d-%s-%s.txt", static_cast<int>(cutval*100), file1.c_str(), file2.c_str()));

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  char option1[100], option2[100],  loption1[100], loption2[100]; 
  int color1(1), fill1(1000), color2(1), fill2(1000), marker1(20), marker2(25); 

  if (string::npos != file1.find("No")) {
    color1 = kBlack; 
    fill1  = 3004; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("No")) {
    color2 = kBlue; 
    fill2  = 3004; 
    fill2  = 3356; 
  }
  if (string::npos != file1.find("Cs")) {
    color1 = kBlack; 
    fill1  = 3005; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("Cs")) {
    color2 = kRed; 
    fill2  = 3005; 
    fill2  = 3356; 
  }

  if (string::npos != file1.find("Sg")) {
    // -- data is first, then signal MC
    color1 = kBlack; 
    fill1  = 0; 
    sprintf(option1, "e"); 
    sprintf(loption1, "p");
    color2 = kBlue; 
    fill2  = 3005; 
    fill2  = 3356; 
    sprintf(option2, "hist"); 
    sprintf(loption2, "f");
  }
  
  // -- now fix overlays of two data files
  sprintf(option1, "e");
  sprintf(loption1, "p");
  
  sprintf(option2,  "hist");
  sprintf(loption2, "f");
  
  if (string::npos != file2.find("Data")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(option2, "e"); 
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 21; 
  }
  
  if (string::npos != file1.find("Mc")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 25; 
  }
 
  //  ofstream OUT("testUL.txt", ios::app);

  fF[file1]->cd(Form("%s", dir1.c_str())); 
  TH1D *h(0), *h1(0), *h2(0);

  AnalysisDistribution a("A_pvz"); 

  fF[file1]->cd(Form("%s/%s", dir1.c_str(), "Isolation")); 

  TCanvas *c1;
  string cut, pdfname; 
  vector<string> doList; 
  TString hname; 

  // -- Suck in all AD's
  TObject *key;
  TIter next(gDirectory->GetListOfKeys());                           
  while ((key = (TObject*)next())) {
    hname = key->GetName();
    if (!hname.Contains("MassSi")) continue;
    hname.ReplaceAll("MassSi", ""); 
    doList.push_back(hname.Data()); 
  }



  for (unsigned int i = 0; i < doList.size(); ++i) {
    cut = Form("%s", doList[i].c_str()); 

    pdfname = Form("%s/%s_%s-%s_sbs-%d_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), static_cast<int>(100*cutval), cut.c_str(), selection.c_str());
    
    cout << pdfname << endl;
    
    fF[file1]->cd(Form("%s/%s", dir1.c_str(), "Isolation")); 
    cout << "==> File1: "; gDirectory->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    // -- separate mm from rest because the former don't need sbs
    if (0 == fMode) {
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h1 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h1 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h1 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h1 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
    }

    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint){
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection.c_str()));
    }

    fF[file2]->cd(Form("%s/%s", dir2.c_str(), "Isolation")); 
    cout << "==> File2: pwd() = "; gDirectory->pwd(); 
    if (0 == fMode) {
      h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h2 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h2 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
      //      if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection.c_str()));
    }

    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    c0->cd(); 
    c0->Clear();
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 

    h1->SetMinimum(0.1);
    h1->SetMaximum(dmax); 
    gPad->SetLogy(0); 

    setHist(h1, color1, marker1, 1.5); 
    setFilledHist(h1, color1, color1, fill1); 
    h1->SetTitle("");
    h1->Draw(option1);
    
    setHist(h2, color2, marker2, 1.5); 
    setFilledHist(h2, color2, color2, fill2); 
    h2->Draw(Form("same%s", option2));

    if (string::npos != cut.find("iso")) {
      double ntot   = h1->Integral(0, h1->GetNbinsX()+1);
      double ucut   = h1->Integral(h1->FindBin(cutval), h1->GetNbinsX()+1);
      double eps1   = ucut/ntot; 
      cout << "====> " << cutval << " at bin " << h1->FindBin(cutval) << endl;
      
      ntot          = h2->Integral(0, h2->GetNbinsX()+1);
      ucut          = h2->Integral(h2->FindBin(cutval), h2->GetNbinsX()+1);
      double eps2   = ucut/ntot; 
      
      tl->DrawLatex(0.2, 0.92, cut.c_str());
      tl->DrawLatex(0.2, 0.80, Form("#epsilon_{D}  = %4.3f", eps1));
      tl->DrawLatex(0.2, 0.72, Form("#epsilon_{M}  = %4.3f", eps2));
      tl->DrawLatex(0.2, 0.64, Form("r_{M/D} = %4.3f", eps2/eps1));
      
      OUT << Form("%s: %4.3f %4.3f %4.3f", cut.c_str(), eps1, eps2, eps2/eps1) << endl;
    }

    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  }

} 








