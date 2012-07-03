#include "plotPU.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TPaveStats.h"

using namespace std; 
using std::string; 

ClassImp(plotPU)

// ----------------------------------------------------------------------
plotPU::plotPU(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fDoUseBDT = false; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 

  fNumbersFileName = fDirectory + "/anaBmm.plotPU." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}

// ----------------------------------------------------------------------
plotPU::~plotPU() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotPU::makeAll(int channels) {
  
  fDoPrint = true; 
  string candName("candAnaBu2JpsiK"); 
  fFile = "NoData"; 
  cd(fFile.c_str());  
  effVsNpv("alpha",   -0.05,   "#alpha < 0.05",      "A", candName.c_str(), "Ao"); 
  effVsNpv("fls3d",    15.0,   "l_{3D}/#sigma(l_{3D}) > 15", "A", candName.c_str(), "Ao"); 
  effVsNpv("iso",      0.80,   "I > 0.80",           "A", candName.c_str(), "Ao"); 
  effVsNpv("chi2dof", -2.00,   "#chi^{2}/dof < 2.0", "A", candName.c_str(), "Ao"); 
  effVsNpv("pchi2dof", 0.10,   "prob > 0.10",        "A", candName.c_str(), "Ao"); 
  effVsNpv("flsxy",    15.0,   "l_{xy}/#sigma > 15", "A", candName.c_str(), "Ao"); 
  effVsNpv("docatrk",  0.015,  "d^{0}_{ca} > 0.015 cm", "A", candName.c_str(), "Ao"); 
  effVsNpv("closetrk", -2,     "N^{close}_{trk} < 2",        "A", candName.c_str(), "Ao"); 
  effVsNpv("lip",      -0.008, "l_{z} < 0.008 cm",      "A", candName.c_str(), "Ao"); 
  effVsNpv("lips",     -2.0,   "l_{z}/#sigma(l_{z}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("ip",       -0.008, "#delta_{3D} < 0.008 cm", "A", candName.c_str(), "Ao"); 
  effVsNpv("ips",      -2.0,   "#delta_{3D}/#sigma(#delta_{3D}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("pvavew8",  0.6,    "< w_{trk} > 0.6",     "A", candName.c_str(), "Ao"); 

	/////////////// this part was commented out ////////////////
//   effVsNpv("iso0",     0.75, "#epsilon(I0>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso1",     0.75, "#epsilon(I1>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso2",     0.75, "#epsilon(I2>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso3",     0.75, "#epsilon(I3>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso4",     0.75, "#epsilon(I4>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso5",     0.75, "#epsilon(I5>0.75)",           "A", candName.c_str(), "Ao"); 
	////////////////////
  candName = "candAnaBs2JpsiPhi"; 
  fFile = "CsData"; 
  cd(fFile.c_str());  
  effVsNpv("alpha",   -0.05,   "#alpha < 0.05",      "A", candName.c_str(), "Ao"); 
  effVsNpv("fls3d",    15.0,   "l_{3D}/#sigma(l_{3D}) > 15", "A", candName.c_str(), "Ao"); 
  effVsNpv("iso",      0.80,   "I > 0.80",           "A", candName.c_str(), "Ao"); 
  effVsNpv("chi2dof", -2.0,    "#chi^{2}/dof < 2.0", "A", candName.c_str(), "Ao"); 
  effVsNpv("pchi2dof", 0.10,   "prob > 0.10",        "A", candName.c_str(), "Ao"); 
  effVsNpv("flsxy",    15.0,   "l_{xy}/#sigma > 15", "A", candName.c_str(), "Ao"); 
  effVsNpv("docatrk", 0.015,   "d^{0}_{ca} > 0.015 cm", "A", candName.c_str(), "Ao"); 
  effVsNpv("closetrk", -2,     "N^{close}_{trk} < 2",        "A", candName.c_str(), "Ao"); 
  effVsNpv("lip",      -0.008, "l_{z} < 0.008 cm",      "A", candName.c_str(), "Ao"); 
  effVsNpv("lips",     -2.0,   "l_{z}/#sigma(l_{z}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("ip",       -0.008, "#delta_{3D} < 0.008 cm", "A", candName.c_str(), "Ao"); 
  effVsNpv("ips",      -2.0,   "#delta_{3D}/#sigma(#delta_{3D}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("pvavew8",  0.6,    "< w_{trk} > 0.6",     "A", candName.c_str(), "Ao"); 
}




// ----------------------------------------------------------------------
void plotPU::effVsNpv(const char *var, double cut, const char *ylabel, const char *chan, const char *dir, const char *selection) {

  TH1D *h[25];
  
  bool larger(true); 
  if (cut < 0) {
    larger = false; 
    cut = -cut; 
  }

  gFile->cd(dir); 
  AnalysisDistribution a("A_pvz"); 
  cout << "==> gDirectory = "; gDirectory->pwd(); 
  c0->Clear();
  c0->Divide(5,5);
  for (int i = 0; i < 25; ++i) {
    c0->cd(i+1);
    // A_Npv0/A_npv0_iso
	  cout << "   directory: " << Form("%s/Npv%d", dir, i) << endl;
    gFile->cd(Form("%s/%s_Npv%d", dir, chan, i));
    cout << "     sbsDistribution(" << Form("%s_npv%d_%s", chan, i, var) << ", " << selection << ") " 
	  << var << (larger?" > ":" < ") << cut 
	  << endl;
    h[i] = a.sbsDistribution(Form("%s_npv%d_%s", chan, i, var), selection);
    // -- check that the bins are not negative. If they are, reset to zero
    for (int ix = 1; ix <= h[i]->GetNbinsX(); ++ix) {
      if (h[i]->GetBinContent(ix) < 0) {
	h[i]->SetBinContent(ix, 0.); 
	h[i]->SetBinError(ix, 0.); 
      }	
    }
	  cout << "Number of entries in the " << i <<"th hist = " << h[i]->GetEntries() << endl;
    h[i]->Draw("hist");
  }

  if (fDoPrint)  
    c0->SaveAs(Form("%s/effVsNpv0-%s-%s-%s-%s-0_%d.pdf", fDirectory.c_str(), fFile.c_str(), dir, chan, var, static_cast<int>(100.*cut)));


	TH1D *h_pass = new TH1D("h_pass","pass",25,0,50);
	TH1D *h_tot = new TH1D("h_tot","tot",25,0,50);
	TEfficiency* pEff = 0;
  double eff(0.); 
  double stat(0.); 
	double ntot(0.);
	double ncut(0.);
  for (int i = 0; i < 25; ++i) {
    HistCutEfficiency a(h[i]); 
    a.fVerbose = 1; 
    a.fIncludeOverflow = 0;
    a.eff(h[i], cut); 

    if (larger) {
      eff  = a.hiEff; 
      stat = a.hiErr; 
		ntot = a.nTot;
		ncut = a.lCut;
    } else {
      eff  = a.loEff; 
      stat = a.loErr; 
		ntot = a.nTot;
		ncut = a.uCut;
    }
    double effE = TMath::Sqrt(stat*stat + 0.02*eff*0.02*eff);
    effE = stat;
    // -- can go above 1, if single bins with too large negative entries are present
    if (eff > 1) {
      eff = 0.; 
      effE = 0.;
    }
    cout << "h[i] : " << h[i]->GetName() << endl;
    cout << "h[i]->FindBin(cut) = " << h[i]->FindBin(cut) << endl;
    cout << "h[i]->GetNbinsX()  = " << h[i]->GetNbinsX() << endl;
    cout << "eff = " << eff << " +- (effE) = " << effE << endl;
	cout << "ntot = " << ntot << ", ncut = " << ncut << endl;
	h_tot->SetBinContent(i+1,ntot);
	h_pass->SetBinContent(i+1,ncut);
  }
	if (TEfficiency::CheckConsistency(*h_pass,*h_tot)) {
		pEff = new TEfficiency(*h_pass,*h_tot);
	}
	pEff->Draw();
//	cout << "****Made eff histo, now will make the pdf..." << endl;

  c0->Clear();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
	c0->Update();
//	gStyle->SetOptFit(1111); 

	TF1* f1 = new TF1("f1","pol0",0,50);
	f1->SetLineWidth(3);

	pEff->Fit(f1);
	pEff->SetStatisticOption(TEfficiency:: kFCP);
//	TPaveStats *st = (TPaveStats*)pEff->GetPaintedGraph()->FindObject("f1");
//	st->SetOptStat(0); 
//	st->SetOptFit(0); 

	pEff->Draw();
	gPad->Update();
//	cout << "Made Efficiency" << endl;
	pEff->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
	pEff->GetPaintedGraph()->SetMinimum(0.75);
	pEff->GetPaintedGraph()->SetMaximum(1.05);
	pEff->GetPaintedGraph()->GetXaxis()->SetTitle("N_{PV}");
	pEff->SetMarkerStyle(20);
	pEff->SetMarkerSize(1.4);
//	cout << "Making the text" << endl;
	pEff->Write();
	pEff->Draw("ap");

    tl->SetTextSize(0.04); 
	cout << "after setting text size" << endl;
	tl->DrawLatex(0.18, 0.92, Form("Prob = %4.3f", 
								   f1->GetProb()));
	
	tl->DrawLatex(0.45, 0.92, ylabel);//Form("Cut = %5.4f", cut));

    tl->SetTextSize(0.05); 
    tl->DrawLatex(0.25, 0.82, Form("#chi^{2}/dof = %3.1f/%i", 
								  f1->GetChisquare(),
								  f1->GetNDF()));

	stamp(0.2, "CMS, 4.9 fb^{-1}", 0.7, "#sqrt{s} = 7 TeV"); 
	
//	TPave *pave = new TPave(25.11976,1.055334,37.65518,1.086642,4,"br");
//	pave->SetFillColor(0);
//	pave->SetLineColor(0);
//	pave->SetLineWidth(6);
//	pave->SetShadowColor(0);
//	pave->Draw();
	
  if (fDoPrint)  
    c0->SaveAs(Form("%s/effVsNpv-%s-%s-%s-%s-0_%d.pdf", fDirectory.c_str(), fFile.c_str(), dir, chan, var, static_cast<int>(100.*cut)));

}
