#ifndef GUARD_LAMBDAEFFREADER_H
#define GUARD_LAMBDAEFFREADER_H

#include <iostream>
#include <vector>
#include <utility>
#include <sstream>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaJet.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "treeReader01.hh"


class lambdaEffReader : public treeReader01
{

public:
    lambdaEffReader(TChain *tree, TString evtClassName);
    ~lambdaEffReader();

    void bookHist();
    void bookReducedTree();
    void startAnalysis();
    void endAnalysis();
    void eventProcessing();
    void fillHist();
    void readCuts(TString filename, int dump = 1);
    void initVariables();
    void doL1stuff();
    void doHLTstuff();

    // some functions to look for relatives
    int getMom(TGenCand *cand, int id, int pos = 1);
    int getDau(TGenCand *cand, int id, int pos = 1);
    int getNDau(TGenCand *cand);
    int getCandId(const TAnaCand *tac, const int id, int pos = 1);
    int getSigId(const TAnaCand *tac, const int id, int pos = 1);
    int getCheckedGenIndex(const TAnaCand *tac, const int id, int pos = 1);
    int getCheckedGenIndex(const int sigid);

    // helper functions
    std::string stripOff(std::string instring, char symbol);
    void setKinematicVariables(TLorentzVector tlv, double &pt, double &eta, double &phi, double &rapidity);

    /*! Converts simple types to strings
      \param i variable to be converted
      \return string
      */
    template <typename T> std::string toString(T i)
    {
	std::ostringstream oss;
	oss << i;
	return oss.str();
    };

    template <typename T> void setCut(T &var, std::string value)
    {
        std::istringstream iss(value);
        iss >> var;
    };

    template <typename T> void setCut(T &var, std::string value, TH1D *hcuts, int ibin, std::string title)
    {
        std::istringstream iss(value);
        iss >> var;
        hcuts->SetBinContent(ibin, var);
        hcuts->GetXaxis()->SetBinLabel(ibin, title.c_str());
    };

    // -- Cut values
    int CUTLbCandidate;

    // -- Variables
    TAnaCand    *fpCand1, *fpCand2;

    // reduced tree
    bool fMom5122; // true if the mom of a 3122/443 is a 5122
    bool fJpsiL0; // true if Lb has at least L0 and Jpsi
    bool fJpsiL0excl; // true if Lb decays in exactly L0 and JPsi
    bool fL0PrPi; // true if L0 decays in pr+pi
    bool fL0PrPiexcl; // true if L0 decays exactly in pr+pi
    bool fJpsiMuMu; // true if the JPsi decays in two muons
    bool fJpsiMuMuexcl; // true if the JPsi decays exactly in two muons

    bool fL0matched; // true if there exists a valid MC match

    // vertex of L0
    double fgenL0vtxR, fgenL0vtxZ;

    // kinematic variables for the particles
    double fgenprpt, fgenpreta, fgenprphi, fgenpry;
    double fgenpipt, fgenpieta, fgenpiphi, fgenpiy;
    double fgenL0pt, fgenL0eta, fgenL0phi, fgenL0y;
    double frecprpt, frecpreta, frecprphi, frecpry;
    double frecpipt, frecpieta, frecpiphi, frecpiy;
    double frecL0pt, frecL0eta, frecL0phi, frecL0y;


    // for trigger stuff
    bool fL1TDMu0, fL1TDMu3;
    bool fL1TMuBH, fL1TMu0, fL1TMu3, fL1TMu5, fL1TMu7, fL1TMu10, fL1TMu14, fL1TMu20;

    bool fHLTqrk, fHLTqrkLS;
    bool fHLTDMuOp, fHLTDMu0, fHLTDMu3, fHLTDMu5;
    bool fHLTMu3t3jp, fHLTMu3t5jp, fHLTMu5t0jp;
    bool fHLTMu0jp, fHLTMu0jpT, fHLTMu3jp, fHLTMu3jpT, fHLTMu5jp, fHLTMu5jpT;
    bool fHLTL1DMu0, fHLTL2DMu0, fHLTL2Mu0;
};

#endif
