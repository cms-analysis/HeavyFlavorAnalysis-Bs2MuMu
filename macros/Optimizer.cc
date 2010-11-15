#include "Optimizer.hh"

#include <iostream>
#include <iomanip>


ClassImp(Optimizer)

// ----------------------------------------------------------------------
Optimizer::Optimizer(TH1D *s, TH1D *b) {

  fpSignal = s; 
  fpBackground = b; 

}


// ----------------------------------------------------------------------
Optimizer::~Optimizer() {

}


// ----------------------------------------------------------------------
TH1D* Optimizer::SoverB(int upperCut) {
  TH1D *h = new TH1D(*fpSignal); h->SetName(Form("SoverB-%s", fpSignal->GetName())); h->Reset();
  
  double s(0.), b(0.);
  if (upperCut) {
    for (int i = 0; i <= fpSignal->GetNbinsX(); ++i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (b>0) {
	h->SetBinContent(i, s/b);
      } else {
	h->SetBinContent(i, 0.);
      }	
    }
  } else {
    for (int i = fpSignal->GetNbinsX(); i >=0 ; --i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (b > 0.) {
	h->SetBinContent(i, s/b);
      } else {
	h->SetBinContent(i, 0.);
      }
    }
  }

  return h;
}

// ----------------------------------------------------------------------
TH1D* Optimizer::EffS2overEffB(int upperCut) {
  TH1D *h = new TH1D(*fpSignal); h->SetName(Form("SoverB-%s", fpSignal->GetName())); h->Reset();
  
  double s(0.), b(0.);
  double S = fpSignal->GetSumOfWeights(), 
    B = fpBackground->GetSumOfWeights(); 
  
  if (upperCut) {
    for (int i = 0; i <= fpSignal->GetNbinsX(); ++i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (S > 0 && B > 0 && b > 0) {
	h->SetBinContent(i,  ((s/S)*(s/S)) / (b/B));
      } else {
	h->SetBinContent(i,  0.);
      }
    }
  } else {
    for (int i = fpSignal->GetNbinsX(); i >=0 ; --i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (S > 0 && B > 0 && b > 0) {
	h->SetBinContent(i, ((s/S)*(s/S)) / (b/B));
      } else {
	h->SetBinContent(i, 0.);
      }
    }
  }

  return h;
}

// ----------------------------------------------------------------------
TH1D* Optimizer::S2overSplusB(int upperCut) {
  TH1D *h = new TH1D(*fpSignal); h->SetName(Form("SoverB-%s", fpSignal->GetName())); h->Reset();
  
  double s(0.), b(0.);

  if (upperCut) {
    for (int i = 0; i <= fpSignal->GetNbinsX(); ++i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (s+b > 0.) {
	h->SetBinContent(i, s*s/(s+b));
      } else {
	h->SetBinContent(i, 0.);
      }	
    }
  } else {
    for (int i = fpSignal->GetNbinsX(); i >=0 ; --i) {
      s += fpSignal->GetBinContent(i);
      b += fpBackground->GetBinContent(i);
      if (s+b > 0.) {
	h->SetBinContent(i, s*s/(s+b));
      } else {
	h->SetBinContent(i, 0.);
      }
    }
  }
  
  return h;
}
  
