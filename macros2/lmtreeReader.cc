#include "lmtreeReader.hh"



using std::cout;
using std::endl;
using std::string;

lmtreeReader::lmtreeReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> lmtreeReader: constructor..." << endl;

  fForceJson = false;
  SAVECANDS = 0;

  Reco_B2MM_size = 0;
  Reco_B2JpsiK_size = 0;
  Reco_B2JpsiPhi_size = 0;
  Reco_Muon_size = 0;

  Gen_Cand_size = 0;
  Gen_Muon_size = 0;

  trigger_name.push_back("HLT_DoubleMu3_Jpsi");
  trigger_name.push_back("HLT_Dimuon6p5_Jpsi");
  trigger_name.push_back("HLT_Dimuon6p5_Jpsi_Displaced");
  trigger_name.push_back("HLT_Dimuon7_Jpsi_Displaced");
  trigger_name.push_back("HLT_DoubleMu3p5_Jpsi_Displaced");
  trigger_name.push_back("HLT_DoubleMu4_Jpsi_Displaced");

  trigger_name.push_back("HLT_DoubleMu3_Bs");
  trigger_name.push_back("HLT_DoubleMu2_Bs");
  trigger_name.push_back("HLT_Dimuon4_Bs_Barrel");
  trigger_name.push_back("HLT_Dimuon6_Bs");
  trigger_name.push_back("HLT_DoubleMu4_Dimuon6_Bs");
  trigger_name.push_back("HLT_DoubleMu4_Dimuon4_Bs_Barrel");

  HLT = new std::vector <std::pair <std::string,bool> >;
  HLT->resize(trigger_name.size());

}

lmtreeReader::~lmtreeReader() {
  cout << "==> lmtreeReader: destructor..." << endl;
}

void lmtreeReader::readCuts(TString filename, int dump) {
  char buffer[1024];
  char cutName[128];
  float cut;
  int ok;
  FILE *cutFile = fopen(filename.Data(), "r");
  if (!cutFile) {
    cout << "cut file not existing: " << filename.Data() << " aborting..." << endl;
    abort();
  }
  if (dump) {
    cout << "===================================================" << endl;
    cout << "==> lmtreeReader: Cut file " << filename.Data() << endl;
    cout << "---------------------------------------------------" << endl;
  }

  while (fgets(buffer, sizeof(buffer), cutFile)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
	// skip comment line
	if (buffer[0] == '#') continue;
	ok = sscanf(buffer, "%s %f", cutName, &cut);
	if(ok < 1) {
	  cout << "==> massReader: Cut not parse input line:" << endl;
	  cout << buffer << endl;
	  cout << "==> massReader: abort..." << endl;
	  abort();
	}
	if (!parseCut(cutName, cut)) {
	  cout << "==> lmtreeReader: Error parsing variable " << cutName << ". abort..." << endl;
	  abort();
	}
  }
  if (dump) cout << "---------------------------------------------------" << endl;
  if (cutFile) fclose(cutFile);
} // readCuts()

bool lmtreeReader::parseCut(char *cutName, float cut, int dump) {
  // parse the options...
  if (!strcmp(cutName,"SAVEALWAYS")) {
    SAVEALWAYS = (int)cut;
    if (dump) cout << "SAVEALWAYS: " << SAVEALWAYS << endl;
    return true;
  }
  if (!strcmp(cutName,"SAVECANDS")) {
    SAVECANDS = (int)cut;
    if (dump) cout << "SAVECANDS: " << SAVECANDS << endl;
    return true;
  }
  if (!strcmp(cutName,"MC")) {
    MC = (int)cut;
    if (dump) cout << "MC: " << MC << endl;
    return true;
  }
  if (!strcmp(cutName,"B2MMTYPE")) {
    B2MMTYPE = (int)cut;
    if (dump) cout << "B2MMTYPE: " << B2MMTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"B2JpsiKTYPE")) {
    B2JpsiKTYPE = (int)cut;
    if (dump) cout << "B2JpsiKTYPE: " << B2JpsiKTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"B2JpsiPhiTYPE")) {
    B2JpsiPhiTYPE = (int)cut;
    if (dump) cout << "B2JpsiPhiTYPE: " << B2JpsiPhiTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"NCANDS")) {
    NCANDS = (int)cut;
    if (dump) cout << "NCANDS: " << NCANDS << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDPT")) {
    CANDPT = cut;
    if (dump) cout << "CANDPT: " << CANDPT << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDETA")) {
    CANDETA = cut;
    if (dump) cout << "CANDETA: " << CANDETA << endl;
    return true;
  }
  if (!strcmp(cutName,"MUONPT")) {
    MUONPT = cut;
    if (dump) cout << "MUONPT: " << MUONPT << endl;
    return true;
  }
  if (!strcmp(cutName,"MUONETA")) {
    MUONETA = cut;
    if (dump) cout << "MUONETA: " << MUONETA << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDMINMASS")) {
    CANDMINMASS = cut;
    if (dump) cout << "CANDMINMASS: " << CANDMINMASS << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDMAXMASS")) {
    CANDMAXMASS = cut;
    if (dump) cout << "CANDMAXMASS: " << CANDMAXMASS << endl;
    return true;
  }
  return false;
} // parseCut()

void lmtreeReader::bookHist() {

  // create the tree
  T1 = new TTree("T1", "Reduced tree for Bs2MuMu studies by lm");

  T1->Branch("run", &fRun, "run/I");
  T1->Branch("LS", &fLS, "LS/I");
  T1->Branch("event", &fEvt, "event/I");
  T1->Branch("Reco_NPV", &Reco_NPV, "Reco_NPV/I");

  for (unsigned int i = 0; i < trigger_name.size(); ++i){
    T1->Branch(trigger_name[i].c_str(), &(HLT->at(i).second), (trigger_name[i] + "/O").c_str());
  }

  Reco_Muon_4mom  = new TClonesArray("TLorentzVector");
  Reco_Muon_PosM2 = new TClonesArray("TVector3");
  T1->Branch("Reco_Muon_size", &Reco_Muon_size, "Reco_Muon_size/I");
  T1->Branch("Reco_Muon_4mom", "TClonesArray", &Reco_Muon_4mom, 32000, 0);
  T1->Branch("Reco_Muon_PosM2", "TClonesArray", &Reco_Muon_PosM2, 32000, 0);
  T1->Branch("Reco_Muon_Charge", Reco_Muon_Charge, "Reco_Muon_Charge[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_Chi2AntiKink", Reco_Muon_Chi2AntiKink, "Reco_Muon_Chi2AntiKink[Reco_Muon_size]/D");
  T1->Branch("Reco_Muon_ID", Reco_Muon_ID, "Reco_Muon_ID[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_pdgId", Reco_Muon_pdgId, "Reco_Muon_pdgId[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_ValidHits", Reco_Muon_ValidHits, "Reco_Muon_ValidHits[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_Index", Reco_Muon_Index, "Reco_Muon_Index[Reco_Muon_size]/I");

  if (SAVECANDS) {
    Reco_B2MM_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2MM_size", &Reco_B2MM_size, "Reco_B2MM_size/I");
    T1->Branch("Reco_B2MM_MCTruth", Reco_B2MM_MCTruth, "Reco_B2MM_MCTruth[Reco_B2MM_size]/O");
    T1->Branch("Reco_B2MM_4mom", "TClonesArray", &Reco_B2MM_4mom, 32000, 0);
 //   T1->Branch("Reco_B2MM_Type", Reco_B2MM_Type, "Reco_B2MM_Type[Reco_B2MM_size]/I");
    T1->Branch("Reco_B2MM_Dxy", Reco_B2MM_Dxy, "Reco_B2MM_Dxy[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DxyE", Reco_B2MM_DxyE, "Reco_B2MM_Dxyz[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Dxyz", Reco_B2MM_Dxyz, "Reco_B2MM_Dxyz[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DxyzE", Reco_B2MM_DxyzE, "Reco_B2MM_DxyzE[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Alpha", Reco_B2MM_Alpha, "Reco_B2MM_Alpha[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Alphaxy", Reco_B2MM_Alphaxy, "Reco_B2MM_Alphaxy[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_VxChi2", Reco_B2MM_VxChi2, "Reco_B2MM_VxChi2[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DOCA", Reco_B2MM_DOCA, "Reco_B2MM_DOCA[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_NDOF", Reco_B2MM_NDOF, "Reco_B2MM_NDOF[Reco_B2MM_size]/I");
    T1->Branch("Reco_B2MM_Iso1_pt09", Reco_B2MM_Iso1_pt09, "Reco_B2MM_Iso1_pt09[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Mu1_index", Reco_B2MM_Mu1_index, "Reco_B2MM_Mu1_index[Reco_B2MM_size]/I");
    T1->Branch("Reco_B2MM_Mu2_index", Reco_B2MM_Mu2_index, "Reco_B2MM_Mu2_index[Reco_B2MM_size]/I");

    Reco_B2JpsiK_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiK_Tk1_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiK_Jpsi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2JpsiK_size", &Reco_B2JpsiK_size, "Reco_B2JpsiK_size/I");
    T1->Branch("Reco_B2JpsiK_MCTruth", Reco_B2JpsiK_MCTruth, "Reco_B2JpsiK_MCTruth[Reco_B2JpsiK_size]/O");
    T1->Branch("Reco_B2JpsiK_4mom", "TClonesArray", &Reco_B2JpsiK_4mom, 32000, 0);
 //   T1->Branch("Reco_B2JpsiK_Type", Reco_B2JpsiK_Type, "Reco_B2JpsiK_Type[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Dxy", Reco_B2JpsiK_Dxy, "Reco_B2JpsiK_Dxy[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DxyE", Reco_B2JpsiK_DxyE, "Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Dxyz", Reco_B2JpsiK_Dxyz, "Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DxyzE", Reco_B2JpsiK_DxyzE, "Reco_B2JpsiK_DxyzE[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Alpha", Reco_B2JpsiK_Alpha, "Reco_B2JpsiK_Alpha[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Alphaxy", Reco_B2JpsiK_Alphaxy, "Reco_B2JpsiK_Alphaxy[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_VxChi2", Reco_B2JpsiK_VxChi2, "Reco_B2JpsiK_VxChi2[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DOCA", Reco_B2JpsiK_DOCA, "Reco_B2JpsiK_DOCA[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_NDOF", Reco_B2JpsiK_NDOF, "Reco_B2JpsiK_NDOF[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Iso1_pt09", Reco_B2JpsiK_Iso1_pt09, "Reco_B2JpsiK_Iso1_pt09[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Tk1_4mom", "TClonesArray", &Reco_B2JpsiK_Tk1_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiK_Tk1_Charge", Reco_B2JpsiK_Tk1_Charge, "Reco_B2JpsiK_Tk1_Charge[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Jpsi_4mom", "TClonesArray", &Reco_B2JpsiK_Jpsi_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiK_Mu1_index", Reco_B2JpsiK_Mu1_index, "Reco_B2JpsiK_Mu1_index[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Mu2_index", Reco_B2JpsiK_Mu2_index, "Reco_B2JpsiK_Mu2_index[Reco_B2JpsiK_size]/I");

    Reco_B2JpsiPhi_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Tk1_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Tk2_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Jpsi_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Phi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2JpsiPhi_size", &Reco_B2JpsiPhi_size, "Reco_B2JpsiPhi_size/I");
    T1->Branch("Reco_B2JpsiPhi_MCTruth", Reco_B2JpsiPhi_MCTruth, "Reco_B2JpsiPhi_MCTruth[Reco_B2JpsiPhi_size]/O");
    T1->Branch("Reco_B2JpsiPhi_4mom", "TClonesArray", &Reco_B2JpsiPhi_4mom, 32000, 0);
 //   T1->Branch("Reco_B2JpsiPhi_Type", Reco_B2JpsiPhi_Type, "Reco_B2JpsiPhi_Type[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Dxy", Reco_B2JpsiPhi_Dxy, "Reco_B2JpsiPhi_Dxy[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DxyE", Reco_B2JpsiPhi_DxyE, "Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Dxyz", Reco_B2JpsiPhi_Dxyz, "Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DxyzE", Reco_B2JpsiPhi_DxyzE, "Reco_B2JpsiPhi_DxyzE[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Alpha", Reco_B2JpsiPhi_Alpha, "Reco_B2JpsiPhi_Alpha[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Alphaxy", Reco_B2JpsiPhi_Alphaxy, "Reco_B2JpsiPhi_Alphaxy[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_VxChi2", Reco_B2JpsiPhi_VxChi2, "Reco_B2JpsiPhi_VxChi2[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DOCA", Reco_B2JpsiPhi_DOCA, "Reco_B2JpsiPhi_DOCA[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_NDOF", Reco_B2JpsiPhi_NDOF, "Reco_B2JpsiPhi_NDOF[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Iso1_pt09", Reco_B2JpsiPhi_Iso1_pt09, "Reco_B2JpsiPhi_Iso1_pt09[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Mu1_index", Reco_B2JpsiPhi_Mu1_index, "Reco_B2JpsiPhi_Mu1_index[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Mu2_index", Reco_B2JpsiPhi_Mu2_index, "Reco_B2JpsiPhi_Mu2_index[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Tk1_4mom", "TClonesArray", &Reco_B2JpsiPhi_Tk1_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Tk1_Charge", Reco_B2JpsiPhi_Tk1_Charge, "Reco_B2JpsiPhi_Tk1_Charge[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Tk2_4mom", "TClonesArray", &Reco_B2JpsiPhi_Tk2_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Tk2_Charge", Reco_B2JpsiPhi_Tk2_Charge, "Reco_B2JpsiPhi_Tk2_Charge[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Jpsi_4mom", "TClonesArray", &Reco_B2JpsiPhi_Jpsi_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Phi_4mom", "TClonesArray", &Reco_B2JpsiPhi_Phi_4mom, 32000, 0);

    Reco_Jpsi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_Jpsi_size", &Reco_Jpsi_size, "Reco_Jpsi_size/I");
    T1->Branch("Reco_Jpsi_MCTruth", Reco_Jpsi_MCTruth, "Reco_Jpsi_MCTruth[Reco_Jpsi_size]/O");
    T1->Branch("Reco_Jpsi_4mom", "TClonesArray", &Reco_Jpsi_4mom, 32000, 0);
 //   T1->Branch("Reco_Jpsi_Type", Reco_Jpsi_Type, "Reco_Jpsi_Type[Reco_Jpsi_size]/I");
    T1->Branch("Reco_Jpsi_Dxy", Reco_Jpsi_Dxy, "Reco_Jpsi_Dxy[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DxyE", Reco_Jpsi_DxyE, "Reco_Jpsi_Dxyz[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Dxyz", Reco_Jpsi_Dxyz, "Reco_Jpsi_Dxyz[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DxyzE", Reco_Jpsi_DxyzE, "Reco_Jpsi_DxyzE[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Alpha", Reco_Jpsi_Alpha, "Reco_Jpsi_Alpha[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Alphaxy", Reco_Jpsi_Alphaxy, "Reco_Jpsi_Alphaxy[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_VxChi2", Reco_Jpsi_VxChi2, "Reco_Jpsi_VxChi2[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DOCA", Reco_Jpsi_DOCA, "Reco_Jpsi_DOCA[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_NDOF", Reco_Jpsi_NDOF, "Reco_Jpsi_NDOF[Reco_Jpsi_size]/I");
    T1->Branch("Reco_Jpsi_Iso1_pt09", Reco_Jpsi_Iso1_pt09, "Reco_Jpsi_Iso1_pt09[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Mu1_index", Reco_Jpsi_Mu1_index, "Reco_Jpsi_Mu1_index[Reco_Jpsi_size]/I");
    T1->Branch("Reco_Jpsi_Mu2_index", Reco_Jpsi_Mu2_index, "Reco_Jpsi_Mu2_index[Reco_Jpsi_size]/I");

  }
  if (MC) {
    Gen_Cand_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Tk1_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Tk2_4mom  = new TClonesArray("TLorentzVector");
    T1->Branch("Gen_Cand_size", &Gen_Cand_size, "Gen_Cand_size/I");
    T1->Branch("Gen_Cand_4mom", "TClonesArray", &Gen_Cand_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Tk1_4mom", "TClonesArray", &Gen_Cand_Tk1_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Tk2_4mom", "TClonesArray", &Gen_Cand_Tk2_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Mu1_index", Gen_Cand_Mu1_index, "Gen_Cand_Mu1_index[Gen_Cand_size]/I");
    T1->Branch("Gen_Cand_Mu2_index", Gen_Cand_Mu2_index, "Gen_Cand_Mu2_index[Gen_Cand_size]/I");
    T1->Branch("Gen_Cand_proc", Gen_Cand_proc, "Gen_Cand_proc[Gen_Cand_size]/I");

    Gen_Muon_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Jpsi_4mom = new TClonesArray("TLorentzVector");
    Gen_Cand_Phi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Gen_Muon_size", &Gen_Muon_size, "Gen_Muon_size/I");
    T1->Branch("Gen_Muon_4mom", "TClonesArray", &Gen_Muon_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Jpsi_4mom", "TClonesArray", &Gen_Cand_Jpsi_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Phi_4mom", "TClonesArray", &Gen_Cand_Phi_4mom, 32000, 0);
    T1->Branch("Gen_Muon_Charge", Gen_Muon_Charge, "Gen_Muon_Charge[Gen_Muon_size]/I");
    T1->Branch("Gen_Muon_Index", Gen_Muon_Index, "Gen_Muon_Index[Gen_Muon_size]/I");
  }

} // massReader::bookHist()

void lmtreeReader::eventProcessing() {

  // Fill Generator Block candidates
  if (MC) fillGen();

  // fill trigger infos
  fillTrigger();

  // Fill Reconstructed Block candidates
  fillMuons();

  // Fill only if there are enough cands
  if (SAVECANDS) {
    fillRecoCand();
  }

  if (SAVEALWAYS || Reco_B2MM_size >= NCANDS || Reco_B2JpsiK_size  >= NCANDS || Reco_B2JpsiPhi_size >= NCANDS) T1->Fill();

  clearVariables();
} // massReader::eventProcessing()

void lmtreeReader::fillGen() {
  int nc = fpEvt->nGenCands();
  for (int j = 0; j < nc; j++) {
    TGenCand* Cand = fpEvt->getGenCand(j);
    if (abs(Cand->fID) == 531 || abs(Cand->fID) == 521) {
      if (abs(Cand->fID) == 531) {
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = k+1; l <= Cand->fDau2; l++) {
            TGenCand* Mu1 = fpEvt->getGenCand(k);
            TGenCand* Mu2 = fpEvt->getGenCand(l);
            if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
              new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
              Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
              Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
              Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
              new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector();
              Gen_Cand_proc[Gen_Cand_size] = 5311313;
              Gen_Cand_size++;
            }
          }
        }
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = Cand->fDau1; l <= Cand->fDau2; l++) {
            TGenCand* jpsi = fpEvt->getGenCand(k);
            TGenCand* phi = fpEvt->getGenCand(l);
            if (abs(jpsi->fID) == 443 && abs(phi->fID) == 333) {
              for (int m = jpsi->fDau1; m <= jpsi->fDau2; m++) {
                for (int n = m+1; n <= jpsi->fDau2; n++) {
                  TGenCand* Mu1 = fpEvt->getGenCand(m);
                  TGenCand* Mu2 = fpEvt->getGenCand(n);
                  if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
                    for (int o = phi->fDau1; o <= phi->fDau2; o++) {
                      for (int p = o+1; p <= phi->fDau2; p++) {
                        TGenCand* ka1 = fpEvt->getGenCand(o);
                        TGenCand* ka2 = fpEvt->getGenCand(p);
                        if ( (ka1->fID == 321 && ka2->fID == -321) || (ka1->fID == -321 && ka2->fID == 321)) {
                          new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
                          Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
                          Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
                          Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
                          new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector(ka1->fP);
                          new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector(ka2->fP);
                          new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector(jpsi->fP);
                          new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector(phi->fP);
                          Gen_Cand_proc[Gen_Cand_size] = 531443333;
                          Gen_Cand_size++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (abs(Cand->fID) == 521) {
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = Cand->fDau1; l <= Cand->fDau2; l++) {
            TGenCand* jpsi = fpEvt->getGenCand(k);
            TGenCand* kaon = fpEvt->getGenCand(l);
            if (abs(jpsi->fID) == 443 && abs(kaon->fID) == 321) {
              for (int m = jpsi->fDau1; m <= jpsi->fDau2; m++) {
                for (int n = m+1; n <= jpsi->fDau2; n++) {
                  TGenCand* Mu1 = fpEvt->getGenCand(m);
                  TGenCand* Mu2 = fpEvt->getGenCand(n);
                  if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
                    new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
                    Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
                    Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
                    Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
                    new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector(kaon->fP);
                    new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector();
                    new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector(jpsi->fP);
                    new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector();
                    Gen_Cand_proc[Gen_Cand_size] = 521443321;
                    Gen_Cand_size++;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (abs(Cand->fID) == 13) {
      new ((*Gen_Muon_4mom)[Gen_Muon_size])TLorentzVector(Cand->fP);
      Gen_Muon_Index[Gen_Muon_size] = Cand->fNumber;
      if (Cand->fNumber != j) cout << "non e' lo stesso numero: fNumber = " << Cand->fNumber << " e j = " << j << endl;
      Gen_Muon_size++;
    }
  }
}

void lmtreeReader::fillTrigger() {

  size_t found;
  bool exists;
  int lenght = sizeof(fpEvt->fHLTNames) / sizeof (TString);
  for (unsigned int k = 0; k < trigger_name.size(); k++) {
    exists = false;
    for (int j = 0; j < lenght; j++) {
      string htl_name(fpEvt->fHLTNames[j].Data());
      found = htl_name.find(trigger_name[k]);
      if (found != string::npos) {
        std::pair < std::string, bool> fire(trigger_name[k], fpEvt->fHLTResult[j]);
        HLT->at(k) = fire;
        exists = true;
        break;
      }
    }
    if (exists == false) {
      std::pair < std::string, bool> fire(trigger_name[k], false);
      HLT->at(k) = fire;
    }
  }
}

void lmtreeReader::fillRecoCand() {

  Reco_NPV = fpEvt->nPV();
  int nc = fpEvt->nCands();

  for (int j = 0; j < nc; j++) {
    TAnaCand* Cand = fpEvt->getCand(j);
    //if ( j != Cand->fIndex) cout << j << " = " << Cand->fIndex << " ?" << endl;
    if (selCand(Cand)) {
      if (Cand->fType == B2MMTYPE)  fillB2MMCand(Cand);
      if (Cand->fType == B2JpsiKTYPE) fillB2JpsiKCand(Cand);
      if (Cand->fType == B2JpsiPhiTYPE) fillB2JpsiPhiCand(Cand);

      if (Reco_B2MM_size >= array_size) {
        cout << "event "<< fEvt <<" B2MM array too small: "<< array_size << " nCands = " << Reco_B2MM_size << endl;
        abort();
      }
    }

    //if (Cand->fType == 300443) {
    //  fillJpsiCand(Cand);
    //  jpsi++;
    //}
    //cout << "reco jpsi = " << jpsi << endl;

  }
  return;
}

void lmtreeReader::fillB2MMCand(TAnaCand* Cand) {

  TAnaTrack *mu0;
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);
  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    mu0 = fpEvt->getSigTrack(it);
    if (TMath::Abs(mu0->fMCID) == 13) {
      if (mu1s == 0) mu1s = mu0;
      else mu2s = mu0;
    }
  }
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2MM_4mom)[Reco_B2MM_size])TLorentzVector(quadriMom);
  //Reco_B2MM_Type[Reco_B2MM_size]   = Cand->fType;
  Reco_B2MM_Dxy[Reco_B2MM_size]    = Cand->fVtx.fDxy;
  Reco_B2MM_DxyE[Reco_B2MM_size]   = Cand->fVtx.fDxyE;
  Reco_B2MM_Dxyz[Reco_B2MM_size]   = Cand->fVtx.fD3d;
  Reco_B2MM_DxyzE[Reco_B2MM_size]  = Cand->fVtx.fD3dE;
  Reco_B2MM_VxChi2[Reco_B2MM_size] = Cand->fVtx.fChi2;
  Reco_B2MM_NDOF[Reco_B2MM_size]   = Cand->fVtx.fNdof;
  Reco_B2MM_DOCA[Reco_B2MM_size]   = Cand->fMaxDoca;
  Reco_B2MM_Mu1_index[Reco_B2MM_size] = mu1_index;
  Reco_B2MM_Mu2_index[Reco_B2MM_size] = mu2_index;

  Reco_B2MM_size++;

  return;
}

void lmtreeReader::fillB2JpsiKCand(TAnaCand* Cand) {

  int counter = 0;
  int zerocounter = 0;
  TAnaTrack *sigtrack0;
  TAnaTrack *kappas(0), *kappar(0);
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);

  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    zerocounter++;
    sigtrack0 = fpEvt->getSigTrack(it);
    if (abs(sigtrack0->fMCID) == 321) {
      kappas = sigtrack0;
      counter++;
    }
    if (abs(sigtrack0->fMCID) == 13) {
      if (mu1s == 0) mu1s = sigtrack0;
      else mu2s = sigtrack0;
      counter++;
    }
  }
  if (counter != 3) {cout << "counter " << counter << "zer0counter " << counter << endl; abort();}
  kappar = fpEvt->getRecTrack(kappas->fIndex);
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);

  if (!selTrack(kappar)) return;
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2JpsiK_4mom)[Reco_B2JpsiK_size])TLorentzVector(quadriMom);
  //Reco_B2JpsiK_Type[Reco_B2JpsiK_size]   = Cand->fType;
  Reco_B2JpsiK_Dxy[Reco_B2JpsiK_size]    = Cand->fVtx.fDxy;
  Reco_B2JpsiK_DxyE[Reco_B2JpsiK_size]   = Cand->fVtx.fDxyE;
  Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]   = Cand->fVtx.fD3d;
  Reco_B2JpsiK_DxyzE[Reco_B2JpsiK_size]  = Cand->fVtx.fD3dE;
  Reco_B2JpsiK_VxChi2[Reco_B2JpsiK_size] = Cand->fVtx.fChi2;
  Reco_B2JpsiK_NDOF[Reco_B2JpsiK_size]   = Cand->fVtx.fNdof;
  Reco_B2JpsiK_DOCA[Reco_B2JpsiK_size]   = Cand->fMaxDoca;
  Reco_B2JpsiK_Mu1_index[Reco_B2JpsiK_size] = mu1_index;
  Reco_B2JpsiK_Mu2_index[Reco_B2JpsiK_size] = mu2_index;

  pt = kappar->fPlab.Pt();
  eta = kappar->fPlab.Eta();
  phi = kappar->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiK_Tk1_4mom)[Reco_B2JpsiK_size])TLorentzVector(quadriMom);
  Reco_B2JpsiK_Tk1_Charge[Reco_B2JpsiK_size] = kappar->fQ;

  pt   = mu1r->fPlab.Pt();
  eta  = mu1r->fPlab.Eta();
  phi  = mu1r->fPlab.Phi();
  TLorentzVector m1_4mom;
  m1_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  pt   = mu2r->fPlab.Pt();
  eta  = mu2r->fPlab.Eta();
  phi  = mu2r->fPlab.Phi();
  TLorentzVector m2_4mom;
  m2_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  TLorentzVector jpsi;
  jpsi = m1_4mom + m2_4mom;
  new ((*Reco_B2JpsiK_Jpsi_4mom)[Reco_B2JpsiK_size])TLorentzVector(jpsi);

  Reco_B2JpsiK_size++;

  return;
}

void lmtreeReader::fillB2JpsiPhiCand(TAnaCand* Cand) {

  int counter = 0;
  int zerocounter = 0;
  TAnaTrack *sigtrack0;
  TAnaTrack *kappa1s(0), *kappa1r(0);
  TAnaTrack *kappa2s(0), *kappa2r(0);

  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);

  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    zerocounter++;
    sigtrack0 = fpEvt->getSigTrack(it);
    if (abs(sigtrack0->fMCID) == 321) {
      if (kappa1s == 0) kappa1s = sigtrack0;
      else kappa2s = sigtrack0;
      counter++;
    }
    if (abs(sigtrack0->fMCID) == 13) {
      if (mu1s == 0) mu1s = sigtrack0;
      else mu2s = sigtrack0;
      counter++;
    }
  }
  if (counter != 4) {cout << "Jpsi Phi counter " << counter << "zer0counter " << counter << endl; abort();}
  kappa1r = fpEvt->getRecTrack(kappa1s->fIndex);
  kappa2r = fpEvt->getRecTrack(kappa2s->fIndex);
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);

  if (!selTrack(kappa1r)) return;
  if (!selTrack(kappa2r)) return;
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2JpsiPhi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  //Reco_B2JpsiPhi_Type[Reco_B2JpsiPhi_size]   = Cand->fType;
  Reco_B2JpsiPhi_Dxy[Reco_B2JpsiPhi_size]    = Cand->fVtx.fDxy;
  Reco_B2JpsiPhi_DxyE[Reco_B2JpsiPhi_size]   = Cand->fVtx.fDxyE;
  Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]   = Cand->fVtx.fD3d;
  Reco_B2JpsiPhi_DxyzE[Reco_B2JpsiPhi_size]  = Cand->fVtx.fD3dE;
  Reco_B2JpsiPhi_VxChi2[Reco_B2JpsiPhi_size] = Cand->fVtx.fChi2;
  Reco_B2JpsiPhi_NDOF[Reco_B2JpsiPhi_size]   = Cand->fVtx.fNdof;
  Reco_B2JpsiPhi_DOCA[Reco_B2JpsiPhi_size]   = Cand->fMaxDoca;
  Reco_B2JpsiPhi_Mu1_index[Reco_B2JpsiPhi_size] = mu1_index;
  Reco_B2JpsiPhi_Mu2_index[Reco_B2JpsiPhi_size] = mu2_index;

  pt = kappa1r->fPlab.Pt();
  eta = kappa1r->fPlab.Eta();
  phi = kappa1r->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiPhi_Tk1_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  Reco_B2JpsiPhi_Tk1_Charge[Reco_B2JpsiPhi_size] = kappa1r->fQ;
  pt = kappa2r->fPlab.Pt();
  eta = kappa2r->fPlab.Eta();
  phi = kappa2r->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiPhi_Tk2_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  Reco_B2JpsiPhi_Tk2_Charge[Reco_B2JpsiPhi_size] = kappa2r->fQ;

  pt   = mu1r->fPlab.Pt();
  eta  = mu1r->fPlab.Eta();
  phi  = mu1r->fPlab.Phi();
  TLorentzVector m1_4mom;
  m1_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  pt   = mu2r->fPlab.Pt();
  eta  = mu2r->fPlab.Eta();
  phi  = mu2r->fPlab.Phi();
  TLorentzVector m2_4mom;
  m2_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  TLorentzVector jpsi;
  jpsi = m1_4mom + m2_4mom;
  new ((*Reco_B2JpsiPhi_Jpsi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(jpsi);
  pt   = kappa1r->fPlab.Pt();
  eta  = kappa1r->fPlab.Eta();
  phi  = kappa1r->fPlab.Phi();
  TLorentzVector kappa1_4mom;
  kappa1_4mom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  pt   = kappa2r->fPlab.Pt();
  eta  = kappa2r->fPlab.Eta();
  phi  = kappa2r->fPlab.Phi();
  TLorentzVector kappa2_4mom;
  kappa2_4mom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  TLorentzVector Phi;
  Phi = kappa1_4mom + kappa2_4mom;
  new ((*Reco_B2JpsiPhi_Phi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(Phi);

  Reco_B2JpsiPhi_size++;

}

void lmtreeReader::fillJpsiCand(TAnaCand* Cand) {

  TAnaTrack *mu0;
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);
  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    mu0 = fpEvt->getSigTrack(it);
    if (TMath::Abs(mu0->fMCID) == 13) {
      if (mu1s == 0) mu1s = mu0;
      else mu2s = mu0;
    }
  }
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();// cout << pt << endl;
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_Jpsi_4mom)[Reco_Jpsi_size])TLorentzVector(quadriMom);
  //Reco_Jpsi_Type[Reco_Jpsi_size]   = Cand->fType;
  Reco_Jpsi_Dxy[Reco_Jpsi_size]    = Cand->fVtx.fDxy;
  Reco_Jpsi_DxyE[Reco_Jpsi_size]   = Cand->fVtx.fDxyE;
  Reco_Jpsi_Dxyz[Reco_Jpsi_size]   = Cand->fVtx.fD3d;
  Reco_Jpsi_DxyzE[Reco_Jpsi_size]  = Cand->fVtx.fD3dE;
  Reco_Jpsi_VxChi2[Reco_Jpsi_size] = Cand->fVtx.fChi2;
  Reco_Jpsi_NDOF[Reco_Jpsi_size]   = Cand->fVtx.fNdof;
  Reco_Jpsi_DOCA[Reco_Jpsi_size]   = Cand->fMaxDoca;
  Reco_Jpsi_Mu1_index[Reco_Jpsi_size] = mu1_index;
  Reco_Jpsi_Mu2_index[Reco_Jpsi_size] = mu2_index;

  Reco_Jpsi_size++;

  return;
}

void lmtreeReader::fillMuons() {

  int nm = fpEvt->nRecTracks();
  for (int j = 0; j < nm; j++) {
    TAnaTrack* Track = fpEvt->getRecTrack(j);
    if (selMuon(Track)) {
      Reco_Muon_Index[Reco_Muon_size] = j;
      Reco_Muon_ValidHits[Reco_Muon_size] = Track->fValidHits;
      double pt   = Track->fPlab.Pt();
      double eta  = Track->fPlab.Eta();
      double phi  = Track->fPlab.Phi();
      TLorentzVector quadriMom;
      quadriMom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
      new ((*Reco_Muon_4mom)[Reco_Muon_size])TLorentzVector(quadriMom);

      Reco_Muon_ID[Reco_Muon_size] = Track->fMuID;
      Reco_Muon_Charge[Reco_Muon_size] = Track->fQ;
      Reco_Muon_pdgId[Reco_Muon_size] = Track->fMCID;

      TAnaMuon* Muon = fpEvt->getMuon(Track->fMuIndex);
      double x = Muon->fPositionAtM2.X();
      double y = Muon->fPositionAtM2.Y();
      double z = Muon->fPositionAtM2.Z();
      TVector3 triVec;
      triVec.SetXYZ(x, y, z);
      new ((*Reco_Muon_PosM2)[Reco_Muon_size])TVector3(triVec);
      Reco_Muon_Chi2AntiKink[Reco_Muon_size] = Muon->fMuonChi2;
      Reco_Muon_size++;
      if (Reco_Muon_size >= Muon_array_size) {
        cout << "muon array too small: "<< Reco_Muon_size << " nMuons = " << nm << endl;
        abort();
      }
    }
  }
  return;
}

bool lmtreeReader::selCand(TAnaCand* Candi) {
  if (! ( (Candi->fType == B2MMTYPE) || (Candi->fType == B2JpsiKTYPE) || (Candi->fType == B2JpsiPhiTYPE))) return false;
  if (! (Candi->fPlab.Pt() > CANDPT)) return false;
  if (! (fabs(Candi->fPlab.Eta()) < CANDETA)) return false;
  if (! (Candi->fMass > CANDMINMASS && Candi->fMass < CANDMAXMASS) ) return false;
  return true;
}

bool lmtreeReader::selMuon(TAnaTrack* Muon) {
  if (! ((Muon->fTrackQuality & 0x1<<2) == 0x1<<2)) return false;  // track is HP
  if (! (Muon->fMuIndex >= 0 && Muon->fMuID >= 0)) return false;  // track is a muon too
  //if (! ((Muon->fMuID & 6) == 6)) return false;  // muon is global and tracker
  if (! (Muon->fPlab.Pt() > MUONPT)) return false;
  if (! (fabs(Muon->fPlab.Eta()) < MUONETA)) return false;
  return true;
}

bool lmtreeReader::selTrack(TAnaTrack* Track) {

  // check the candidate
  if (! ((Track->fTrackQuality & 0x1<<2) == 0x1<<2)) return false;  // track is HP

  return true;
}

void lmtreeReader::clearVariables() {

  Reco_NPV = 0;

  if (MC) {
    if (Gen_Cand_size > 0) {
      Gen_Cand_4mom->Clear();
      Gen_Cand_Phi_4mom->Clear();
      Gen_Cand_Jpsi_4mom->Clear();
      Gen_Cand_Tk1_4mom->Clear();
      Gen_Cand_Tk2_4mom->Clear();
    }
    if (Gen_Muon_size > 0) Gen_Muon_4mom->Clear();
    Gen_Cand_size = 0;
    Gen_Muon_size = 0;
  }

  if (SAVECANDS) {
    if (Reco_B2MM_size > 0) Reco_B2MM_4mom->Clear();
    if (Reco_B2JpsiK_size > 0) {
      Reco_B2JpsiK_4mom->Clear();
      Reco_B2JpsiK_Tk1_4mom->Clear();
      Reco_B2JpsiK_Jpsi_4mom->Clear();
    }
    if (Reco_B2JpsiPhi_size > 0) {
      Reco_B2JpsiPhi_4mom->Clear();
      Reco_B2JpsiPhi_Tk1_4mom->Clear();
      Reco_B2JpsiPhi_Tk2_4mom->Clear();
      Reco_B2JpsiPhi_Jpsi_4mom->Clear();
      Reco_B2JpsiPhi_Phi_4mom->Clear();
    }
    if (Reco_Jpsi_size > 0) Reco_Jpsi_4mom->Clear();
    Reco_B2MM_size = 0;
    Reco_B2JpsiK_size = 0;
    Reco_B2JpsiPhi_size = 0;
    Reco_Jpsi_size = 0;
  }

  if (Reco_Muon_size > 0) {
    Reco_Muon_4mom->Clear();
    Reco_Muon_PosM2->Clear();
  }
  Reco_Muon_size = 0;

}

void lmtreeReader::closeHistFile() {

  fpHistFile = T1->GetCurrentFile();
  treeReader01::closeHistFile();

} // massReader::closeHistFile()
