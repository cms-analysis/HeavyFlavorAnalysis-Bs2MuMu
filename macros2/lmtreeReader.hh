#ifndef LMTREEREADER_HH
#define LMTREEREADER_HH

#include <utility>
#include <iostream>
#include <cmath>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "../macros/treeReader01.hh"

const static double MUON_MASS = 0.1057;
const static double KAON_MASS = 0.4937;

const static int array_size = 256;
const static int Muon_array_size = 1024;

class lmtreeReader : public treeReader01 {
  public:
    lmtreeReader(TChain *tree, TString evtClassName);
    virtual ~lmtreeReader();

    virtual void readCuts(TString filename, int dump = 1);
    virtual bool parseCut(char *cutName, float cut, int dump = 1);
    virtual void bookHist();
    virtual void eventProcessing();
    virtual void fillGen();
    virtual void fillTrigger();
    virtual void fillRecoCand();
    virtual void fillB2MMCand(TAnaCand* Candi);
    virtual void fillB2JpsiKCand(TAnaCand* Candi);
    virtual void fillB2JpsiPhiCand(TAnaCand* Candi);
    virtual void fillJpsiCand(TAnaCand* Candi);
    virtual void fillMuons();
    virtual bool selCand(TAnaCand* Candi);
    virtual bool selTrack(TAnaTrack* Track);
    virtual bool selMuon(TAnaTrack* Muon);

    virtual double isoClassicWithDOCA(TAnaCand*, double dca, double r = 1.0, double ptmin = 0.9);
    virtual int nCloseTracks(TAnaCand*, double dca, double pt = 0.5);


    virtual void clearVariables();
    virtual void closeHistFile();

  protected:
    //cuts
    int SAVECANDS, NCANDS, B2MMTYPE, B2JpsiKTYPE, B2JpsiPhiTYPE, CANDCLOSETRK;
    double PVAVEW8, PVLIP, PVLIPS;
    double CANDMINMASS, CANDMAXMASS, CANDPT, CANDETA;
    double MUONPT, MUONETA;
    // options
    int SAVEALWAYS, MC;

    std::vector<std::string> trigger_name;

    //output tree
    TTree* T1;
    /// RECO
    int Reco_NPV;



    TClonesArray* Reco_B2MM_4mom;
    int Reco_B2MM_size;
    double Reco_B2MM_Dxyz[array_size];
    double Reco_B2MM_DxyzE[array_size];
    double Reco_B2MM_Dxy[array_size];
    double Reco_B2MM_DxyE[array_size];
    double Reco_B2MM_Alpha[array_size];
    double Reco_B2MM_Alphaxy[array_size];
    double Reco_B2MM_VxRChi2[array_size];
    double Reco_B2MM_DOCAtrk[array_size];
    int Reco_B2MM_Type[array_size];
    double Reco_B2MM_Iso1_pt09[array_size];
    bool Reco_B2MM_MCTruth[array_size];
    int Reco_B2MM_Mu1_index[array_size];
    int Reco_B2MM_Mu2_index[array_size];

    TClonesArray* Reco_B2JpsiK_4mom;
    TClonesArray* Reco_B2JpsiK_Jpsi_4mom;
    TClonesArray* Reco_B2JpsiK_Tk1_4mom;
    int Reco_B2JpsiK_size;
    double Reco_B2JpsiK_Dxyz[array_size];
    double Reco_B2JpsiK_DxyzE[array_size];
    double Reco_B2JpsiK_Dxy[array_size];
    double Reco_B2JpsiK_DxyE[array_size];
    double Reco_B2JpsiK_Alpha[array_size];
    double Reco_B2JpsiK_Alphaxy[array_size];
    double Reco_B2JpsiK_VxRChi2[array_size];
    double Reco_B2JpsiK_DOCAtrk[array_size];
    int Reco_B2JpsiK_Type[array_size];
    double Reco_B2JpsiK_Iso1_pt09[array_size];
    bool Reco_B2JpsiK_MCTruth[array_size];
    int Reco_B2JpsiK_Mu1_index[array_size];
    int Reco_B2JpsiK_Mu2_index[array_size];
    int Reco_B2JpsiK_Tk1_Charge[array_size];

    TClonesArray* Reco_B2JpsiPhi_4mom;
    TClonesArray* Reco_B2JpsiPhi_Jpsi_4mom;
    TClonesArray* Reco_B2JpsiPhi_Phi_4mom;
    TClonesArray* Reco_B2JpsiPhi_Tk1_4mom;
    TClonesArray* Reco_B2JpsiPhi_Tk2_4mom;
    int Reco_B2JpsiPhi_size;
    double Reco_B2JpsiPhi_Dxyz[array_size];
    double Reco_B2JpsiPhi_DxyzE[array_size];
    double Reco_B2JpsiPhi_Dxy[array_size];
    double Reco_B2JpsiPhi_DxyE[array_size];
    double Reco_B2JpsiPhi_Alpha[array_size];
    double Reco_B2JpsiPhi_Alphaxy[array_size];
    double Reco_B2JpsiPhi_VxRChi2[array_size];
    double Reco_B2JpsiPhi_DOCAtrk[array_size];
    int Reco_B2JpsiPhi_Type[array_size];
    double Reco_B2JpsiPhi_Iso1_pt09[array_size];
    bool Reco_B2JpsiPhi_MCTruth[array_size];
    int Reco_B2JpsiPhi_Mu1_index[array_size];
    int Reco_B2JpsiPhi_Mu2_index[array_size];
    int Reco_B2JpsiPhi_Tk1_Charge[array_size];
    int Reco_B2JpsiPhi_Tk2_Charge[array_size];

    TClonesArray* Reco_Jpsi_4mom;
    int Reco_Jpsi_size;
    double Reco_Jpsi_Dxyz[array_size];
    double Reco_Jpsi_DxyzE[array_size];
    double Reco_Jpsi_Dxy[array_size];
    double Reco_Jpsi_DxyE[array_size];
    double Reco_Jpsi_Alpha[array_size];
    double Reco_Jpsi_Alphaxy[array_size];
    double Reco_Jpsi_VxRChi2[array_size];
    double Reco_Jpsi_DOCAtrk[array_size];
    double Reco_Jpsi_Iso1_pt09[array_size];
    bool Reco_Jpsi_MCTruth[array_size];
    int Reco_Jpsi_Mu1_index[array_size];
    int Reco_Jpsi_Mu2_index[array_size];

    int Reco_Muon_size;
    TClonesArray* Reco_Muon_4mom;
    TClonesArray* Reco_Muon_PosM2;
    int Reco_Muon_Charge[Muon_array_size];
    double Reco_Muon_Chi2AntiKink[Muon_array_size];
    int Reco_Muon_ID[Muon_array_size];
    int Reco_Muon_ValidHits[Muon_array_size];
    int Reco_Muon_PxValidHits[Muon_array_size];
    int Reco_Muon_pdgId[Muon_array_size];
    int Reco_Muon_Index[Muon_array_size];

    /// GEN
    int Gen_Cand_size;
    TClonesArray* Gen_Cand_4mom;
    TClonesArray* Gen_Cand_Tk1_4mom;
    TClonesArray* Gen_Cand_Tk2_4mom;
    TClonesArray* Gen_Cand_Jpsi_4mom;
    TClonesArray* Gen_Cand_Phi_4mom;
    int Gen_Cand_pdgId[array_size];
    int Gen_Cand_Mu1_index[array_size];
    int Gen_Cand_Mu2_index[array_size];
    int Gen_Cand_proc[array_size];

    TClonesArray* Gen_Muon_4mom;
    int Gen_Muon_size;
    int Gen_Muon_Charge[Muon_array_size];
    int Gen_Muon_Index[Muon_array_size];

    std::vector<std::pair <std::string, bool> > *HLT;

};

#endif // LMTREEREADER_HH
