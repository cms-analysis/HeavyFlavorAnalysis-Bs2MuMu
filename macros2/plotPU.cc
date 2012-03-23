#include "plotPU.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"

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
  effVsNpv("lip",      -0.008, "l_{z}<0.008 cm",      "A", candName.c_str(), "Ao"); 
  effVsNpv("lips",     -2.0,   "l_{z}/#sigma(l_{z}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("ip",       -0.008, "#delta_{3D} < 0.008 cm", "A", candName.c_str(), "Ao"); 
  effVsNpv("ips",      -2.0,   "#delta_{3D}/#sigma(#delta_{3D}) < 2",  "A", candName.c_str(), "Ao"); 
  effVsNpv("pvavew8",  0.6,    "<w_{trk} > 0.6",     "A", candName.c_str(), "Ao"); 

//   effVsNpv("iso0",     0.75, "#epsilon(I0>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso1",     0.75, "#epsilon(I1>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso2",     0.75, "#epsilon(I2>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso3",     0.75, "#epsilon(I3>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso4",     0.75, "#epsilon(I4>0.75)",           "A", candName.c_str(), "Ao"); 
//   effVsNpv("iso5",     0.75, "#epsilon(I5>0.75)",           "A", candName.c_str(), "Ao"); 

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
  effVsNpv("pvavew8",  0.6,    "<w_{trk} > 0.6",     "A", candName.c_str(), "Ao"); 
}




// ----------------------------------------------------------------------
void plotPU::effVsNpv(const char *var, double cut, const char *ylabel, const char *chan, const char *dir, const char *selection) {

  TH1D *h[15];
  
  bool larger(true); 
  if (cut < 0) {
    larger = false; 
    cut = -cut; 
  }

  gFile->cd(dir); 
  AnalysisDistribution a("A_pvz"); 
  cout << "==> gDirectory = "; gDirectory->pwd(); 
  c0->Clear();
  c0->Divide(4,4);
  for (int i = 0; i < 15; ++i) {
    c0->cd(i+1);
    // A_Npv0/A_npv0_iso
    //    cout << "   directory: " << Form("%s/Npv%d", dir, i) << endl;
    gFile->cd(Form("%s/%s_Npv%d", dir, chan, i));
    //     cout << "     sbsDistribution(" << Form("%s_npv%d_%s", chan, i, var) << ", " << selection << ") " 
    // 	 << var << (larger?" > ":" < ") << cut 
    // 	  << endl;
    h[i] = a.sbsDistribution(Form("%s_npv%d_%s", chan, i, var), selection);
    // -- check that the bins are not negative. If they are, reset to zero
    for (int ix = 1; ix <= h[i]->GetNbinsX(); ++ix) {
      if (h[i]->GetBinContent(ix) < 0) {
	h[i]->SetBinContent(ix, 0.); 
	h[i]->SetBinError(ix, 0.); 
      }	
    }

    h[i]->Draw("hist");
  }

  if (fDoPrint)  
    c0->SaveAs(Form("%s/effVsNpv0-%s-%s-%s-%s-0_%d.pdf", fDirectory.c_str(), fFile.c_str(), dir, chan, var, static_cast<int>(100.*cut)));


  TH1D *heff = new TH1D("heff", "", 15, 0., 30.);
  double eff(0.); 
  double stat(0.); 
  for (int i = 0; i < 15; ++i) {
    HistCutEfficiency a(h[i]); 
    a.fVerbose = 1; 
    a.fIncludeOverflow = 0;
    a.eff(h[i], cut); 

    if (larger) {
      eff  = a.hiEff; 
      stat = a.hiErr; 
    } else {
      eff  = a.loEff; 
      stat = a.loErr; 
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
    //    cout << "ntot = " << ntot << endl;
    //     cout << "ncut = " << ncut << endl;
    cout << "eff = " << eff << endl;
    cout << "effE = " << effE << endl;
    heff->SetBinContent(i+1, eff); 
    heff->SetBinError(i+1, effE); 
  }

  c0->Clear();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  //  setTitles(heff, "N_{PV}", Form("%s", ylabel), 0.05); 
  setTitles(heff, "N_{PV}", "efficiency", 0.06, 1.1, 1.4); 
  heff->SetMaximum(1.05);
  if (heff->GetMinimum() > 0.5) {
    heff->SetMinimum(0.75);
  } else {
    heff->SetMinimum(0.);
  }
  heff->SetMarkerSize(2); 
  heff->Fit("pol0", "FM");
  heff->GetFunction("pol0")->SetLineWidth(3);
  tl->SetTextSize(0.07); 
  tl->DrawLatex(0.25, 0.2, ylabel); 
  tl->SetTextSize(0.05); 
//   tl->DrawLatex(0.55, 0.91, Form("#epsilon = %4.3f#pm%4.3f", 
// 				 heff->GetFunction("pol0")->GetParameter(0), 
// 				 heff->GetFunction("pol0")->GetParError(0))); 

//   tl->DrawLatex(0.25, 0.82, Form("#chi^{2}/dof = %3.1f/%i", 
// 				 heff->GetFunction("pol0")->GetChisquare(), 
// 				 heff->GetFunction("pol0")->GetNDF())); 

  stamp(0.18, "CMS, 5 fb^{-1}", 0.65, "#sqrt{s} = 7 TeV"); 
  //  stamp(0.18, "CMS, 5.3 fb^{-1}", 0.65, "#sqrt{s} = 7 TeV"); 
  
  if (fDoPrint)  
    c0->SaveAs(Form("%s/effVsNpv-%s-%s-%s-%s-0_%d.pdf", fDirectory.c_str(), fFile.c_str(), dir, chan, var, static_cast<int>(100.*cut)));


  heff->Fit("pol1", "FM");
  heff->GetFunction("pol1")->SetLineWidth(3);
  tl->SetTextSize(0.04); 
  tl->DrawLatex(0.16, 0.91, Form("p0 = %4.3f#pm%4.3f", 
				 heff->GetFunction("pol1")->GetParameter(0), 
				 heff->GetFunction("pol1")->GetParError(0))); 

  tl->DrawLatex(0.50, 0.91, Form("p1 = %5.4f#pm%5.4f", 
				 heff->GetFunction("pol1")->GetParameter(1), 
				 heff->GetFunction("pol1")->GetParError(1))); 

  tl->SetTextSize(0.05); 
  tl->DrawLatex(0.25, 0.82, Form("#chi^{2}/dof = %3.1f/%i", 
				 heff->GetFunction("pol1")->GetChisquare(), 
				 heff->GetFunction("pol1")->GetNDF())); 

  //  stamp(0.2, "CMS, 4.9 fb^{-1}", 0.7, "#sqrt{s} = 7 TeV"); 
  
  if (fDoPrint)  
    c0->SaveAs(Form("%s/effVsNpv-%s-%s-%s-%s-1_%d.pdf", fDirectory.c_str(), fFile.c_str(), dir, chan, var, static_cast<int>(100.*cut)));

}
