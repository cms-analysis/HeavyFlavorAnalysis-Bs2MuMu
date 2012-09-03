#include "preselection.hh"
#include "TString.h"

#define PTMIN 6.0
#define PTMAX 70.0

#define M1PTMIN 4.0
#define M1PTMAX 50.0

#define M2PTMIN 4.0
#define M2PTMAX 20.0

#define FL3DMAX 1.5
#define CHI2DOFMAX 8.0

#define PVIPMAX 0.04
#define PVIPSMAX 10.0

#define MAXDOCAMAX 0.03
#define CLOSETRKMAX 21

#define FLS3DMAX 100.0
#define DOCATRKMAX 0.15

#define ISOMIN 0.6
#define ALPHAMAX 0.3    

// ----------------------------------------------------------------------
std::string preselection() {

//   TCut sgcut = "hlt&&gmuid&&pt<70&&pt>6"; 
//   sgcut += "m1pt>4.0&&m1pt<50";
//   sgcut += "m2pt>4.0&&m2pt<20"; 
//   sgcut += "fl3d<2&&chi2/dof<10&&pvip<0.02&&!TMath::IsNaN(pvips)&&pvips<5&&pvips>0&&maxdoca<0.02";
//   sgcut += "closetrk<21"; 
//   sgcut += "fls3d<100&&docatrk<0.2";
//   sgcut += "alpha<0.3&&iso>0.6&&chi2/dof<10"; 

  std::string cut = Form("hlt && gmuid && (%3.2f<pt)&&(pt<%3.2f) && (%3.2f<m1pt)&&(m1pt<%3.2f) && (%3.2f<m2pt)&&(m2pt<%3.2f)", 
			 PTMIN, PTMAX, M1PTMIN, M1PTMAX, M2PTMIN, M2PTMAX);
  cut += std::string(Form(" && (fl3d<%3.2f) && (pvip < %3.2f) && !(TMath::IsNaN(pvips)) && (pvips>0) && (pvips<%3.2f)", 
			  FL3DMAX, PVIPMAX, PVIPSMAX));
  cut += std::string(Form(" && (closetrk<%d) && (fls3d<%3.2f) && (docatrk<%3.2f) && (maxdoca<%f)",
			  CLOSETRKMAX, FLS3DMAX, DOCATRKMAX, MAXDOCAMAX));
  cut += std::string(Form(" && (chi2/dof<%3.2f) && (iso>%3.2f) && (alpha<%3.2f)", 
			  CHI2DOFMAX, ISOMIN, ALPHAMAX)); 
  return cut; 
}

// ----------------------------------------------------------------------
bool preselection(RedTreeData &b, int channel, bool rejectInvIso) {
  //   if (fTrainAntiMuon && b.gmuid) return false;
  //   if (!fTrainAntiMuon && !b.gmuid) return false;
  if (!b.hlt) return false;
  if (b.pt > PTMAX) return false;
  if (b.pt < PTMIN) return false;

  if (b.m1pt < M1PTMIN) return false;
  if (b.m1pt > M1PTMAX) return false;
  if (b.m2pt < M2PTMIN) return false;
  if (b.m2pt > M2PTMAX) return false;

  if (b.fl3d > FL3DMAX) return false;
  if (b.chi2/b.dof > CHI2DOFMAX) return false;
  if (b.pvip > PVIPMAX) return false;
  if (b.pvips < 0) return false;
  if (b.pvips > PVIPSMAX) return false;
  if (TMath::IsNaN(b.pvips)) return false;
  if (b.maxdoca > MAXDOCAMAX) return false;

  if (b.closetrk > CLOSETRKMAX) return false;

  if (b.fls3d > FLS3DMAX) return false;
  if (b.docatrk > DOCATRKMAX) return false;

  if (rejectInvIso && 5.2 < b.m && b.m < 5.45 && b.iso < 0.7) return false;

  if (b.m < 4.9) return false;
  if (b.m > 5.9) return false;

  // -- physics preselection: reduce background by factor 7, signal efficiency >90%
  if (b.chi2/b.dof > 10) return false;
  if (b.iso<0.6) return false; 
  if (b.alpha > 0.3) return false; 

  if (0 == channel) {
    if (TMath::Abs(b.m1eta) > 1.4) return false;
    if (TMath::Abs(b.m2eta) > 1.4) return false;
  } else if (1 == channel) {
    if (TMath::Abs(b.m1eta) < 1.4 && TMath::Abs(b.m2eta) < 1.4) return false;
    if (TMath::Abs(b.m1eta) > 2.4 || TMath::Abs(b.m2eta) > 2.4) return false;
  }

  return true;
}

