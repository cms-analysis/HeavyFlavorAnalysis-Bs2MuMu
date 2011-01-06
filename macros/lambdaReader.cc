#include "lambdaReader.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/lambdaReader.default.cuts
//           ./runTreeReaders -f test.root
// ----------------------------------------------------------------------

/*! Converts simple types to strings
  \param i variable to be converted
  \return string
  */
template <typename T>
std::string toString(T i)
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

/*! Strips off trailing characters from a string
  /param instring String where characters should be cut off
  /param symbol Symbol after which text should be cut off, including the symbol
  /return remaining string
  */
std::string stripOff(std::string instring, char symbol)
{
    std::string::iterator iter = std::find(instring.begin(), instring.end(), symbol);
    std::string outstring;
    std::copy(instring.begin(), iter, std::back_inserter(outstring));
    return outstring;
}

//==========================================================================
decayMap::decayMap() { };

decayMap::pos_t decayMap::getPos(std::string key)
{
    pos_t ret;
    map_t::iterator pos = myMap.find(key);
    if(pos!=myMap.end())
    {
        ret=pos->second;
    }
    else
    {
        ret=myMap.size();
        myMap[key]=ret;
    }
    return ret;
}

void decayMap::printMap()
{
    createRevMap();
    cout << "map contains " << myMap.size() << " entries" << endl;
    cout << "revmap contains " << myRevMap.size() << " entries" << endl;
    for(revmap_t::const_iterator it = myRevMap.begin(); it!=myRevMap.end(); it++)
    {
        cout << it->first << " - " << it->second << endl;
    }
}

void decayMap::createRevMap()
{
    myRevMap.empty();
    for(map_t::const_iterator it = myMap.begin(); it!=myMap.end(); it++)
    {
        myRevMap[it->second]=it->first;
    }
}

unsigned int decayMap::readFile(std::string filename)
{
    std::ifstream infile(filename.c_str());
    if(!infile.good())
    {
        throw std::invalid_argument("File " + filename + " not found.");
    }
    std::string s;
    unsigned int counter(0);
    while (getline(infile,s))
    {
        getPos(s);
        counter++;
    }
    return counter;
}

//==========================================================================
class findDaughters
{
public:
    findDaughters() : nDaughters(0) {};
    findDaughters(TGenCand* curParticle, TAna01Event* fpEvt) : nDaughters(0)
    {
        strDaughters = "";
        strDauGrdau = toString(abs(curParticle->fID)) + " -> ";
        //strDauGrdau = "";
        const int nGenCands(fpEvt->nGenCands());
        mother = curParticle;
        if(mother->fDau1 > -1 && mother->fDau1 < nGenCands && mother->fDau2 > -1 && mother->fDau2 < nGenCands && mother->fDau1 < mother->fDau2)
        {
            nDaughters = 0;
            nDauGrdau = 0;
            for(int dauit=mother->fDau1; dauit<=mother->fDau2; dauit++)
            {
                // <= is ok because fDau2 is id of last daughter in list
                TGenCand *gcDau = fpEvt->getGenCand(dauit);
                const unsigned int id = abs(gcDau->fID);
                // check if these particles come from the right mother
                if(gcDau->fMom1==mother->fNumber&&gcDau->fMom2==mother->fNumber)
                {
                    nDaughters++;
                    strDaughters += toString(id) + " ";
                    findDaughters grdau(gcDau,fpEvt);
                    strDauGrdau += toString(id) + " (" + grdau.getDaughters() + ") " ;
                    strDauGrdaus += toString(id) + " (" + grdau.getDauGrdau() + ") " ;
                    nDauGrdau += grdau.getNDaughters();
                }
                // else: particle comes from hadronic interactions or decays in flight
                // discard this for the moment
            }
        }
    };
    std::string getDaughters()
    {
        return strDaughters;
    };
    std::string getDauGrdau()
    {
        return strDauGrdau;
    };
    std::string getDauGrdaus()
    {
        return strDauGrdaus;
    };
    unsigned int getNDaughters()
    {
        return nDaughters;
    };
    unsigned int getNDauGrdau()
    {
        return nDauGrdau;
    };

private:
    TGenCand* mother;
    std::string strDaughters;
    std::string strDauGrdau;
    std::string strDauGrdaus;
    unsigned int nDaughters; // no. of direct daughters
    unsigned int nDauGrdau; // no. of daughters of the direct daughters
};

//==========================================================================

lambdaReader::lambdaReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName)
{
    cout << "==> lambdaReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
lambdaReader::~lambdaReader()
{
    cout << "==> lambdaReader: destructor..." << endl;
}

// ----------------------------------------------------------------------
void lambdaReader::startAnalysis()
{
    cout << "==> lambdaReader: Starting analysis..." << endl;
    // fill the map with some anchor values
    if(CUTReadDecayMaps)
    {
        cout << "Reading " << CUTDecayMap1 << " - " << myDecayMap1Gen.readFile(CUTDecayMap1.c_str()) << " entries read." << endl;
        cout << "Reading " << CUTDecayMap2 << " - " << myDecayMap2Gen.readFile(CUTDecayMap2.c_str()) << " entries read." << endl;
    }
    // book reduced trees
    bookReducedTree();
}

// ----------------------------------------------------------------------
void lambdaReader::endAnalysis()
{
    if (CUTPrintDecayMaps)
    {
        cout << "Decay map 1Gen: " << endl;
        myDecayMap1Gen.printMap();
        cout << "Decay map 2Gen: " << endl;
        myDecayMap2Gen.printMap();
    }
}

// ----------------------------------------------------------------------
void lambdaReader::eventProcessing()
{

    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());
    ((TH1D*)fpHistFile->Get("h2"))->Fill(fpEvt->nCands());

    TAnaCand *pCand;

    // initialize a flag to mark in the genTree if there exists a corresponding candidate from signal
    fghasCand = 0;

    // MC only
    if (fIsMC) doGenLevelStuff();

    // do some trigger stuff
    doL1stuff();
    doHLTstuff();

    // main loop
    if (CUTLbCandidateFit == 0 || CUTLbCandidateFit == CUTLbCandidate)
    {
        // we fill everything from the same candidate
        int candcounter(0);
        for (int iC = 0; iC < fpEvt->nCands(); ++iC)
        {
            pCand = fpEvt->getCand(iC);
            ((TH1D*)fpHistFile->Get("h3"))->Fill(pCand->fType);

            if (CUTLbCandidate == pCand->fType)
            {
                CheckedLbCand curCand = getCheckedLbCand(pCand);
                if (curCand.isOk)
                {
                    if (doCandStuff(curCand) && doCandFitStuff(curCand)) fTree->Fill();
                    candcounter++;
                }
            }
        }
        ((TH1D*)fpHistFile->Get("h4"))->Fill(candcounter);
    }
    else
    {
        // we take data from different candidates, needs matching
        // populate the list of candidates
        std::vector<CheckedLbCand> candVec, candVecFit;
        for (int iC = 0; iC < fpEvt->nCands(); ++iC)
        {
            pCand = fpEvt->getCand(iC);
            ((TH1D*)fpHistFile->Get("h3"))->Fill(pCand->fType);
            if (CUTLbCandidate == pCand->fType)
            {
                CheckedLbCand curCand = getCheckedLbCand(pCand);
                if (curCand.isOk) candVec.push_back(curCand);
            }
            if (CUTLbCandidateFit == pCand->fType)
            {
                CheckedLbCand curCand = getCheckedLbCand(pCand);
                if (curCand.isOk) candVecFit.push_back(curCand);
            }
        }
        std::vector<std::pair<CheckedLbCand, CheckedLbCand> > candPairVec;
        for(std::vector<CheckedLbCand>::iterator it1 = candVec.begin(); it1!=candVec.end(); it1++)
            for(std::vector<CheckedLbCand>::iterator it2 = candVecFit.begin(); it2!=candVecFit.end(); it2++)
            {
                if(compareCands( *it1,*it2 ))
                    candPairVec.push_back(std::make_pair((*it1),(*it2)));
            }
        ((TH1D*)fpHistFile->Get("h4"))->Fill(candVec.size());
        ((TH1D*)fpHistFile->Get("h5"))->Fill(candVecFit.size());
        ((TH1D*)fpHistFile->Get("h6"))->Fill(candPairVec.size());
        for(std::vector<std::pair<CheckedLbCand,CheckedLbCand> >::iterator it = candPairVec.begin();
                it != candPairVec.end(); it++)
        {
            if (doCandStuff((*it).first) && doCandFitStuff((*it).second)) fTree->Fill();
        }
    }

    if (fIsMC && fIsSig) fGenTree->Fill();
    return;
}

// ----------------------------------------------------------------------
int lambdaReader::getSigId(const TAnaCand *tac, const int id, int pos)
{
    int ret(-1);
    //cout << "getSigId: " << tac->fSig1 << "-" << tac->fSig2 << " ";
    for(int i = tac->fSig1; i <= tac->fSig2; i++)
    {
        const TAnaTrack* tat = fpEvt->getSigTrack(i);
        if(TMath::Abs(tat->fMCID) == id) pos--;
        if(0 == pos)
        {
            ret = i;
            break;
        }
    }
    //cout << "found: " << ret << endl;
    return ret;
}

// ----------------------------------------------------------------------
int lambdaReader::getCandId(const TAnaCand *tac, const int id, int pos)
{
    int ret(-1);
    //cout << "getCandId: " << tac->fDau1 << "-" << tac->fDau2 << " ";
    for(int i = tac->fDau1; i <= tac->fDau2; i++)
    {
        const TAnaCand* tacDau = fpEvt->getCand(i);
        const int tacDauId = tacDau->fType - int(tacDau->fType / 10000) * 10000;
        if(TMath::Abs(tacDauId) == id) pos--;
        if(0 == pos)
        {
            ret = i;
            break;
        }
    }
    //cout << "found: " << id << " at " << ret << endl;
    return ret;
}

// ----------------------------------------------------------------------
// check if two candidates point to the same two sig tracks for each daughter
bool lambdaReader::compareCands(const CheckedLbCand clc1, const CheckedLbCand clc2)
{
    const int clc1mu1 = fpEvt->getSigTrack(clc1.mu1)->fIndex;
    const int clc1mu2 = fpEvt->getSigTrack(clc1.mu2)->fIndex;
    const int clc1pi = fpEvt->getSigTrack(clc1.pi)->fIndex;
    const int clc1pr = fpEvt->getSigTrack(clc1.pr)->fIndex;

    const int clc2mu1 = fpEvt->getSigTrack(clc2.mu1)->fIndex;
    const int clc2mu2 = fpEvt->getSigTrack(clc2.mu2)->fIndex;
    const int clc2pi = fpEvt->getSigTrack(clc2.pi)->fIndex;
    const int clc2pr = fpEvt->getSigTrack(clc2.pr)->fIndex;

    //cout << "compareCands: " << clc1mu1 << " " << clc1mu2 << " " << clc1pi << " " << clc1pr << " : "
    // << clc2mu1 << " " << clc2mu2 << " " << clc2pi << " " << clc2pr << endl;

    return ((clc1mu1 == clc2mu1)
            && (clc1mu2 == clc2mu2)
            && (clc1pi == clc2pi)
            && (clc1pr == clc2pr));
}

// ----------------------------------------------------------------------
// check if the candidates are complete
lambdaReader::CheckedLbCand lambdaReader::getCheckedLbCand(const TAnaCand *tac)
{
    CheckedLbCand ret;
    ret.isOk = false;
    ret.lb = tac->fIndex;
    // check if daughters point to something
    if (tac->fDau1 < 0 || tac->fDau2 < 0)
    {
        cout << "Problem: tac->fDau1: " << tac->fDau1 << " tac->fDau2: " << tac->fDau2
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }

    // now let's see if we have a jp and a l0
    ret.jp = getCandId(tac,443,1);
    ret.l0 = getCandId(tac,3122,1);
    if (ret.jp < 0 || ret.l0 < 0)
    {
        cout << "Problem finding jp/l0: ret.jp: " << ret.jp << " ret.l0: " << ret.l0
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }

    // check if daughters signal tracks exist
    TAnaCand *tacDau1 = fpEvt->getCand(tac->fDau1);
    if (tacDau1->fSig1 < 0 || tacDau1->fSig2 < 0)
    {
        cout << "Problem: tacDau1->fSig1: " << tacDau1->fSig1 << " tacCur->fDaui->fSig2: " << tacDau1->fSig2
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }
    const TAnaTrack *tatDau1Sig1 = fpEvt->getSigTrack(tacDau1->fSig1);
    const TAnaTrack *tatDau1Sig2 = fpEvt->getSigTrack(tacDau1->fSig2);
    if (tatDau1Sig1->fIndex < 0 || tatDau1Sig2->fIndex < 0)
    {
        cout << "Problem: tatDau1Sig1->fIndex: " << tatDau1Sig1->fIndex << " tatDau1Sig2->fIndex: " << tatDau1Sig2->fIndex
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }
    // and the same for the 2nd daughter
    TAnaCand *tacDau2 = fpEvt->getCand(tac->fDau2);
    if (tacDau2->fSig1 < 0 || tacDau2->fSig2 < 0)
    {
        cout << "Problem: tacDau2->fSig1: " << tacDau2->fSig1 << " tacCur->fDaui->fSig2: " << tacDau2->fSig2
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }
    const TAnaTrack *tatDau2Sig1 = fpEvt->getSigTrack(tacDau2->fSig1);
    const TAnaTrack *tatDau2Sig2 = fpEvt->getSigTrack(tacDau2->fSig2);
    if (tatDau2Sig1->fIndex < 0 || tatDau2Sig2->fIndex < 0)
    {
        cout << "Problem: tatDau2Sig1->fIndex: " << tatDau2Sig1->fIndex << " tatDau2Sig2->fIndex: " << tatDau2Sig2->fIndex
             << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
        return ret;
    }

    // ok, everything seems to point to something meaningful. Now tests we we can get our decay signals properly
    ret.mu1 = getSigId(fpEvt->getCand(ret.jp),13,1);
    ret.mu2 = getSigId(fpEvt->getCand(ret.jp),13,2);
    ret.pi = getSigId(fpEvt->getCand(ret.l0),211,1);;
    ret.pr = getSigId(fpEvt->getCand(ret.l0),2212,1);;
    // cout << ret.jp << " " << ret.l0 << " " << ret.mu1 << " " << ret.mu2 << " " << ret.pi << " " << ret.pr << endl;
    ret.isOk = (ret.mu1 > -1) && (ret.mu2 > -1) && (ret.pi > -1) && (ret.pr > -1);
    return ret;
}

// ----------------------------------------------------------------------
bool lambdaReader::doCandStuff(const CheckedLbCand &clc)
{
    const TAnaCand *tacCur = fpEvt->getCand(clc.lb);
    // ---------------------------------------------------------
    // First we get some handier objects for later calculations
    const TAnaCand *tacJpsi = fpEvt->getCand(clc.jp);
    const TAnaCand *tacLambda0 = fpEvt->getCand(clc.l0);

    // get the muons
    const TAnaTrack *tatSigMu1 = fpEvt->getSigTrack(clc.mu1);
    const TAnaTrack *tatSigMu2 = fpEvt->getSigTrack(clc.mu2);
    const TAnaTrack *tatRecMu1 = fpEvt->getRecTrack(tatSigMu1->fIndex);
    const TAnaTrack *tatRecMu2 = fpEvt->getRecTrack(tatSigMu2->fIndex);

    // require the muons to be of at least some quality
    if(tatRecMu1->fMuID==-1||(tatRecMu1->fMuID & CUTMuId1)!=CUTMuId1) return false;
    if(tatRecMu2->fMuID==-1||(tatRecMu2->fMuID & CUTMuId2)!=CUTMuId2) return false;

    // get the lambda daughters
    const TAnaTrack *tatSigPi = fpEvt->getSigTrack(clc.pi);
    const TAnaTrack *tatSigPr = fpEvt->getSigTrack(clc.pr);
    const TAnaTrack *tatRecPi = fpEvt->getRecTrack(tatSigPi->fIndex);
    const TAnaTrack *tatRecPr = fpEvt->getRecTrack(tatSigPr->fIndex);

    // get some TLorentzVector for eta and phi (and reuse it later)
    TLorentzVector tlvRecMu1, tlvRecMu2, tlvRecPr, tlvRecPi;
    tlvRecMu1.SetVectM(tatRecMu1->fPlab, MMUON);
    tlvRecMu2.SetVectM(tatRecMu2->fPlab, MMUON);
    tlvRecPr.SetVectM(tatRecPr->fPlab, MPROTON);
    tlvRecPi.SetVectM(tatRecPi->fPlab, MPION);

    TLorentzVector tlvLambdaB, tlvLambda0, tlvJpsi;
    tlvLambdaB.SetVectM(tacCur->fPlab, MLAMBDA_B);
    tlvLambda0.SetVectM(tacLambda0->fPlab, MLAMBDA_0);
    tlvJpsi.SetVectM(tacJpsi->fPlab, MJPSI);

    // sigtracks
    TLorentzVector tlvSigMu1, tlvSigMu2, tlvSigPr, tlvSigPi;
    tlvSigMu1.SetVectM(tatSigMu1->fPlab, MMUON);
    tlvSigMu2.SetVectM(tatSigMu2->fPlab, MMUON);
    tlvSigPr.SetVectM(tatSigPr->fPlab,MPROTON);
    tlvSigPi.SetVectM(tatSigPi->fPlab,MPION);

    // ---------------------------------------------------------
    // now we do some calculations

    // K_s hypothesis
    TLorentzVector tlvPrAsPiHypo;
    tlvPrAsPiHypo.SetVectM(tatRecPr->fPlab, MPION);
    const TLorentzVector tlvKsHypoth = tlvRecPi+tlvPrAsPiHypo;

    fghasCand = 1; // obviously we have a candidate, so we want to mark this

    // ---------------------------------------------------------
    // now comes the boring part: filling of the reduced tree

    // candidate data from fit
    fml0=tacLambda0->fMass;
    fptl0=tacLambda0->fPlab.Perp();
    fpl0=tacLambda0->fPlab.Mag();
    fetal0=tlvLambda0.Eta();
    fphil0=tlvLambda0.Phi();

    fmjp=tacJpsi->fMass;
    fptjp=tacJpsi->fPlab.Perp();
    fpjp=tacJpsi->fPlab.Mag();
    fetajp=tlvJpsi.Eta();
    fphijp=tlvJpsi.Phi();

    // reco tracks
    frpt1m=tatRecMu1->fPlab.Perp();
    freta1m=tlvRecMu1.Eta();
    frphi1m=tlvRecMu1.Phi();
    frq1m=tatRecMu1->fQ;
    if(tatRecMu1->fMuID==-1) frid1m=0;
    else frid1m=tatRecMu1->fMuID;
    fchi21m=tatRecMu1->fChi2;
    fprob1m=TMath::Prob(tatRecMu1->fChi2,tatRecMu1->fDof);
    fndof1m=tatRecMu1->fDof;
    fqual1m=tatRecMu1->fTrackQuality;

    frpt2m=tatRecMu2->fPlab.Perp();
    freta2m=tlvRecMu2.Eta();
    frphi2m=tlvRecMu2.Phi();
    frq2m=tatRecMu2->fQ;
    if(tatRecMu2->fMuID==-1) frid2m=0;
    else frid2m=tatRecMu2->fMuID;
    fchi22m=tatRecMu2->fChi2;
    fprob2m=TMath::Prob(tatRecMu2->fChi2,tatRecMu2->fDof);
    fndof2m=tatRecMu2->fDof;
    fqual2m=tatRecMu2->fTrackQuality;

    frptpr=tatRecPr->fPlab.Perp();
    fretapr=tlvRecPr.Eta();
    frphipr=tlvRecPr.Phi();
    frqpr=tatRecPr->fQ;
    fchi2pr=tatRecPr->fChi2;
    fprobpr=TMath::Prob(tatRecPr->fChi2,tatRecPr->fDof);
    fndofpr=tatRecPr->fDof;
    fqualpr=tatRecPr->fTrackQuality;

    frptpi=tatRecPi->fPlab.Perp();
    fretapi=tlvRecPi.Eta();
    frphipi=tlvRecPi.Phi();
    frqpi=tatRecPi->fQ;
    fchi2pi=tatRecPi->fChi2;
    fprobpi=TMath::Prob(tatRecPi->fChi2,tatRecPi->fDof);
    fndofpi=tatRecPi->fDof;
    fqualpi=tatRecPi->fTrackQuality;

    // qualities of vertex fits
    fchi2l0=tacLambda0->fVtx.fChi2;
    fchi2jp=tacJpsi->fVtx.fChi2;
    fndofl0=tacLambda0->fVtx.fNdof;
    fndofjp=tacJpsi->fVtx.fNdof;
    fprobl0=TMath::Prob(fchi2l0,fndofl0);
    fprobjp=TMath::Prob(fchi2jp,fndofjp);

    // doca
    fmaxdocalb=tacCur->fMaxDoca;
    fmaxdocal0=tacLambda0->fMaxDoca;
    fmaxdocajp=tacJpsi->fMaxDoca;

    // dR and angles
    fdRprpi=tlvRecPr.DeltaR(tlvRecPi);
    fdRmumu=tlvRecMu1.DeltaR(tlvRecMu2);

    // K_s hypothesis
    fKshypo=tlvKsHypoth.M();

    return true;
}

// ----------------------------------------------------------------------
bool lambdaReader::doCandFitStuff(const CheckedLbCand &clc)
{
    const TAnaCand *tacCur = fpEvt->getCand(clc.lb);
    // ---------------------------------------------------------
    // First we get some handier objects for later calculations
    const TAnaCand *tacJpsi = fpEvt->getCand(clc.jp);
    const TAnaCand *tacLambda0 = fpEvt->getCand(clc.l0);

    // get the muons
    const TAnaTrack *tatSigMu1 = fpEvt->getSigTrack(clc.mu1);
    const TAnaTrack *tatSigMu2 = fpEvt->getSigTrack(clc.mu2);
    const TAnaTrack *tatRecMu1 = fpEvt->getRecTrack(tatSigMu1->fIndex);
    const TAnaTrack *tatRecMu2 = fpEvt->getRecTrack(tatSigMu2->fIndex);

    // require the muons to be of at least some quality
    if(tatRecMu1->fMuID==-1||(tatRecMu1->fMuID & CUTMuId1)!=CUTMuId1) return false;
    if(tatRecMu2->fMuID==-1||(tatRecMu2->fMuID & CUTMuId2)!=CUTMuId2) return false;

    // get the lambda daughters
    const TAnaTrack *tatSigPi = fpEvt->getSigTrack(clc.pi);
    const TAnaTrack *tatSigPr = fpEvt->getSigTrack(clc.pr);
    const TAnaTrack *tatRecPi = fpEvt->getRecTrack(tatSigPi->fIndex);
    const TAnaTrack *tatRecPr = fpEvt->getRecTrack(tatSigPr->fIndex);

    // get some TLorentzVector for eta and phi (and reuse it later)
    TLorentzVector tlvRecMu1, tlvRecMu2, tlvRecPr, tlvRecPi;
    tlvRecMu1.SetVectM(tatRecMu1->fPlab, MMUON);
    tlvRecMu2.SetVectM(tatRecMu2->fPlab, MMUON);
    tlvRecPr.SetVectM(tatRecPr->fPlab, MPROTON);
    tlvRecPi.SetVectM(tatRecPi->fPlab, MPION);

    TLorentzVector tlvLambdaB, tlvLambda0, tlvJpsi;
    tlvLambdaB.SetVectM(tacCur->fPlab, MLAMBDA_B);
    tlvLambda0.SetVectM(tacLambda0->fPlab, MLAMBDA_0);
    tlvJpsi.SetVectM(tacJpsi->fPlab, MJPSI);

    // sigtracks
    TLorentzVector tlvSigMu1, tlvSigMu2, tlvSigPr, tlvSigPi;
    tlvSigMu1.SetVectM(tatSigMu1->fPlab, MMUON);
    tlvSigMu2.SetVectM(tatSigMu2->fPlab, MMUON);
    tlvSigPr.SetVectM(tatSigPr->fPlab,MPROTON);
    tlvSigPi.SetVectM(tatSigPi->fPlab,MPION);

    // ---------------------------------------------------------
    // now we do some calculations

    // calculate alpha (pointing angles in rads)
    const TVector3 vecJpsiL0 = tacLambda0->fVtx.fPoint - tacJpsi->fVtx.fPoint;
    const TVector3 vecPV = fpEvt->getPV(tacCur->fPvIdx)->fPoint;
    const TVector3 vecPVLb = tacJpsi->fVtx.fPoint - vecPV;
    const double alphal0 = vecJpsiL0.Angle(tacLambda0->fPlab);
    const double alphalb = vecPVLb.Angle(tacCur->fPlab);

    // calculate beta factor as 3d vector
    const TVector3 vecBetaL0 = vecJpsiL0 * (1. / tlvLambda0.E());
    const TVector3 vecBetaLb = vecPVLb   * (1. / tlvLambdaB.E());

    fghasCand = 1; // obviously we have a candidate, so we want to mark this

    // ---------------------------------------------------------
    // now comes the boring part: filling of the reduced tree

    // candidate data from fit
    fmlb=tacCur->fMass;
    fptlb=tacCur->fPlab.Perp();
    fplb=tacCur->fPlab.Mag();
    fetalb=tlvLambdaB.Eta();
    fphilb=tlvLambdaB.Phi();

    // refitted tracks, "signal" tracks
    fSpt1m=tlvSigMu1.Perp();
    fSeta1m=tlvSigMu1.Eta();
    fSphi1m=tlvSigMu1.Phi();

    fSpt2m=tlvSigMu2.Perp();
    fSeta2m=tlvSigMu2.Eta();
    fSphi2m=tlvSigMu2.Phi();

    fSptpr=tlvSigPr.Perp();
    fSetapr=tlvSigPr.Eta();
    fSphipr=tlvSigPr.Phi();

    fSptpi=tlvSigPi.Perp();
    fSetapi=tlvSigPi.Eta();
    fSphipi=tlvSigPi.Phi();

    // flight lengths and their errors
    fd3lb=tacCur->fVtx.fD3d;
    fd3l0=tacLambda0->fVtx.fD3d;
    fd3jp=tacJpsi->fVtx.fD3d;
    fd3Elb=tacCur->fVtx.fD3dE;
    fd3El0=tacLambda0->fVtx.fD3dE;
    fd3Ejp=tacJpsi->fVtx.fD3dE;

    fdxylb=tacCur->fVtx.fDxy;
    fdxyl0=tacLambda0->fVtx.fDxy;
    fdxyjp=tacJpsi->fVtx.fDxy;
    fdxyElb=tacCur->fVtx.fDxyE;
    fdxyEl0=tacLambda0->fVtx.fDxyE;
    fdxyEjp=tacJpsi->fVtx.fDxyE;

    // qualities of vertex fits
    fchi2lb=tacCur->fVtx.fChi2;
    fndoflb=tacCur->fVtx.fNdof;
    fproblb=TMath::Prob(fchi2lb,fndoflb);

    // ctau - written as dist/p*m as this is numerically more stable than using
    // the TLorentzVector::Gamma() functions to write beta*gamma
    fctlb = tacCur->fVtx.fD3d / tacCur->fPlab.Mag() * MLAMBDA_B;
    fctl0 = tacLambda0->fVtx.fD3d / tacLambda0->fPlab.Mag() * MLAMBDA_0;

    // vertex positions
    fvxlb = tacCur->fVtx.fPoint.x();
    fvylb = tacCur->fVtx.fPoint.y();
    fvzlb = tacCur->fVtx.fPoint.z();
    fvrlb = tacCur->fVtx.fPoint.XYvector().Mod();
    fvxl0 = tacLambda0->fVtx.fPoint.x();
    fvyl0 = tacLambda0->fVtx.fPoint.y();
    fvzl0 = tacLambda0->fVtx.fPoint.z();
    fvrl0 = tacLambda0->fVtx.fPoint.XYvector().Mod();

    // beta vector components
    fbtlbx=vecBetaL0.x();
    fbtlby=vecBetaL0.y();
    fbtlbz=vecBetaL0.z();
    fbtl0x=vecBetaLb.x();
    fbtl0y=vecBetaLb.y();
    fbtl0z=vecBetaLb.z();

    // other cut values
    // pointing angles
    falphalb=alphalb;
    falphal0=alphal0;

    // dR and angles
    fdRl0jp=tlvLambda0.DeltaR(tlvLambdaB);
    fanglbl0=tlvLambda0.Angle(tlvLambdaB.Vect());

    // Do truth matching
    if (fIsMC)
    {
	// get the muons
	const TAnaTrack *tatSigMu1 = fpEvt->getSigTrack(clc.mu1);
	const TAnaTrack *tatSigMu2 = fpEvt->getSigTrack(clc.mu2);
	const TAnaTrack *tatRecMu1 = fpEvt->getRecTrack(tatSigMu1->fIndex);
	const TAnaTrack *tatRecMu2 = fpEvt->getRecTrack(tatSigMu2->fIndex);
	// get the lambda daughters
	const TAnaTrack *tatSigPi = fpEvt->getSigTrack(clc.pi);
	const TAnaTrack *tatSigPr = fpEvt->getSigTrack(clc.pr);
	const TAnaTrack *tatRecPi = fpEvt->getRecTrack(tatSigPi->fIndex);
	const TAnaTrack *tatRecPr = fpEvt->getRecTrack(tatSigPr->fIndex);

	fIsMCmatch = doTruthMatchingLb(tatRecMu1, tatRecMu2, tatRecPi, tatRecPr, tacLambda0->fVtx.fPoint);
    }

    return true;
}

// ----------------------------------------------------------------------
bool lambdaReader::doTruthMatchingLb(const TAnaTrack *Mu1, const TAnaTrack *Mu2, const TAnaTrack *Pi, const TAnaTrack *Pr, const TVector3 &vtx)
{
    bool ret(false);
    TLorentzVector tlvMu1, tlvMu2, tlvPi, tlvPr;
    tlvMu1.SetVectM(Mu1->fPlab, MMUON);
    tlvMu2.SetVectM(Mu2->fPlab, MMUON);
    tlvPi.SetVectM(Pi->fPlab, MPION);
    tlvPr.SetVectM(Pr->fPlab, MPROTON);

    const int giMu1(fpEvt->getGenIndexWithDeltaR(tlvMu1,Mu1->fQ));
    const int giMu2(fpEvt->getGenIndexWithDeltaR(tlvMu2,Mu2->fQ));
    const int giPi(fpEvt->getGenIndexWithDeltaR(tlvPi,vtx,Pi->fQ));
    const int giPr(fpEvt->getGenIndexWithDeltaR(tlvPr,vtx,Pr->fQ));
    const int nGenCands(fpEvt->nGenCands());
    if( giMu1 < nGenCands && giMu2 < nGenCands && giPi < nGenCands && giPr < nGenCands &&
            giMu1 > -1 && giMu2 > -1 && giPi > -1 && giPr > -1 )
    {
        TGenCand *gcMu1 = fpEvt->getGenCand(giMu1);
	TGenCand *gcMu1Mom = fpEvt->getGenCand(gcMu1->fMom1);
        TGenCand *gcMu2 = fpEvt->getGenCand(giMu2);
	TGenCand *gcMu2Mom = fpEvt->getGenCand(gcMu2->fMom1);
        TGenCand *gcPion = fpEvt->getGenCand(giPi);
	TGenCand *gcPionMom = fpEvt->getGenCand(gcPion->fMom1);
        TGenCand *gcProton = fpEvt->getGenCand(giPr);
	TGenCand *gcProtonMom = fpEvt->getGenCand(gcProton->fMom1);

        ret = (TMath::Abs(gcMu1->fID) == 13
	       && TMath::Abs(gcMu1Mom->fID) == 443
               && TMath::Abs(gcMu2->fID) == 13
	       && TMath::Abs(gcMu2Mom->fID) == 443
               && TMath::Abs(gcPion->fID) == 211
               && TMath::Abs(gcPionMom->fID) == 3122
               && TMath::Abs(gcProton->fID) == 2212
	       && TMath::Abs(gcProtonMom->fID) == 3122);

        if (fVerbose > 5)
        {
            cout << "Truth matching: Mu1: " << gcMu1->fID
                 << " Mu2: " << gcMu2->fID
                 << " Pion: " << gcPion->fID
                 << " Proton: " << gcProton->fID
                 << " isMatch: " << fIsMCmatch << " isSig: " << fIsSig << endl;
        }
    }
    return ret;
}

// ----------------------------------------------------------------------
void lambdaReader::doGenLevelStuff()
{
    // Do some generator level stuff
    const int nGenCands(fpEvt->nGenCands());
    for(int gcit=0; gcit!=nGenCands; gcit++)
    {
        TGenCand *gc5122 = fpEvt->getGenCand(gcit);
        if (511==abs(gc5122->fID) || 521==abs(gc5122->fID) || 513==abs(gc5122->fID) || 523==abs(gc5122->fID))
        {
            findDaughters fd(gc5122,fpEvt);
        }
        if(5122==abs(gc5122->fID))
        {
            // reset tree variables
            fgmlb = fgmlbsw = fgml0 = fgml0sw = 9999;
            fgptpr = fgptpi = fgptmu1 = fgptmu2 = 9999;
            fgetapr = fgetapi = fgetamu1 = fgetamu2 = 9999;
            fgphipr = fgphipi = fgphimu1 = fgphimu2 = 9999;
            fgpmu1 = fgpmu2 = fgppr = fgppi = 9999;
            fgptl0 = fgpl0 = 9999;
            fgdRprpi = fgdRmumu = fgdRl0lb = 9999;
            fganprpi = fganmumu = fganl0jp = 9999;
            fgnDaughters = fgnGrandDaughters = fgmapRef1Gen = -1;
            // analyse the decay and extract some data to store in the tree
            findDaughters fd(gc5122,fpEvt);
            fIsSig = ("5122 -> 443 (13 13 ) 3122 (211 2212 ) "==fd.getDauGrdau()) ? 1 : 0;
            fnDaughters = fgnDaughters = fd.getNDaughters();
            fnGrandDaughters = fgnGrandDaughters = fd.getNDauGrdau();
            fmapRef1Gen = fgmapRef1Gen = myDecayMap1Gen.getPos(fd.getDaughters());
            fmapRef2Gen = fgmapRef2Gen = myDecayMap2Gen.getPos(fd.getDauGrdau());
            // now extract some kinematic data
            TLorentzVector tlvGenJp, tlvGenL0, tlvGenPr, tlvGenPi, tlvGenMu1, tlvGenMu2, tlvGenLambda0, tlvGenLambdaB;
            TVector3 tv3PV, tv3LbV, tv3L0V;
            if (fIsSig)
            {
                tv3PV = gc5122->fV; // vertices in candblock denote where the particle was produced
                for(int gcDauit=gc5122->fDau1; gcDauit<=gc5122->fDau2; gcDauit++)
                {
                    TGenCand* gcDau=fpEvt->getGenCand(gcDauit);
                    if(443==abs(gcDau->fID))
                    {
                        tv3LbV = gcDau->fV;
                        fgvxlb = tv3LbV.x();
                        fgvylb = tv3LbV.y();
                        fgvzlb = tv3LbV.z();
                        fgvrlb = tv3LbV.XYvector().Mod();
                        // ct. Observe that Lb is in gc5122
                        fgctlb = (tv3LbV-tv3PV).Mag() / gc5122->fP.Vect().Mag() * MLAMBDA_B;
                        unsigned int mucounter(0);
                        tlvGenJp=gcDau->fP;
                        for(int gcGrDauit=gcDau->fDau1; gcGrDauit<=gcDau->fDau2; gcGrDauit++)
                        {
                            TGenCand* gcGrDau=fpEvt->getGenCand(gcGrDauit);
                            if(13==abs(gcGrDau->fID)&&gcGrDau->fMom1==gcDau->fNumber&&gcGrDau->fMom2==gcDau->fNumber)
                            {
                                mucounter++;
                                if (mucounter==1) tlvGenMu1=gcGrDau->fP;
                                if (mucounter==2) tlvGenMu2=gcGrDau->fP;
                            }
                        }
                        if(mucounter>=2)
                        {
                            fgptmu1 = tlvGenMu1.Pt();
                            fgptmu2 = tlvGenMu2.Pt();
                            fgpmu1 = tlvGenMu1.P();
                            fgpmu2 = tlvGenMu2.P();
                            fgetamu1 = tlvGenMu1.Eta();
                            fgetamu2 = tlvGenMu2.Eta();
                            fgphimu1 = tlvGenMu1.Phi();
                            fgphimu2 = tlvGenMu2.Phi();
                            fgdRmumu = tlvGenMu1.DeltaR(tlvGenMu2);
                            fganmumu = tlvGenMu1.Angle(tlvGenMu2.Vect());
                        }
                    }
                    if(3122==abs(gcDau->fID))
                    {
                        tlvGenL0=gcDau->fP;
                        bool prFound(false), piFound(false);
                        for(int gcGrDauit=gcDau->fDau1; gcGrDauit<=gcDau->fDau2; gcGrDauit++)
                        {
                            TGenCand* gcGrDau=fpEvt->getGenCand(gcGrDauit);
                            if(!prFound&&2212==abs(gcGrDau->fID)&&gcGrDau->fMom1==gcDau->fNumber&&gcGrDau->fMom2==gcDau->fNumber)
                            {	// in case we have more than one pr we get the first
                                tlvGenPr=gcGrDau->fP;
                                prFound=true;

                                // fill in data for the L0 decay vertex
                                tv3L0V = gcGrDau->fV;
                                fgvxl0 = tv3L0V.x();
                                fgvyl0 = tv3L0V.y();
                                fgvzl0 = tv3L0V.z();
                                fgvrl0 = tv3L0V.XYvector().Mod();
                                // ct. Observe that L0 is in gcDau
                                fgctl0 = (tv3L0V-tv3PV).Mag() / gcDau->fP.Vect().Mag() * MLAMBDA_0;
                            }
                            if(!piFound&&211==abs(gcGrDau->fID)&&gcGrDau->fMom1==gcDau->fNumber&&gcGrDau->fMom2==gcDau->fNumber)
                            {	// in case we have more than one pi we get the first
                                tlvGenPi=gcGrDau->fP;
                                piFound=true;
                            }
                        }
                        if(prFound&&piFound) // we have a pr and a pi
                        {
                            fgptpr = tlvGenPr.Pt();
                            fgppr = tlvGenPr.P();
                            fgetapr = tlvGenPr.Eta();
                            fgphipr = tlvGenPr.Phi();

                            fgptpi = tlvGenPi.Pt();
                            fgppi = tlvGenPi.P();
                            fgetapi = tlvGenPi.Eta();
                            fgphipi = tlvGenPi.Phi();

                            fgdRprpi = tlvGenPr.DeltaR(tlvGenPi);
                            fganprpi = tlvGenPr.Angle(tlvGenPi.Vect());
                            // make a L0
                            tlvGenLambda0 = tlvGenPr + tlvGenPi;
                            fgml0 = tlvGenLambda0.M();
                            fgptl0 = tlvGenLambda0.Pt();
                            fgpl0 = tlvGenLambda0.P();
                            // make a L0 with swapped mass assumptions for pr and pi
                            TLorentzVector tlvGenPrAsPi, tlvGenPiAsPr;
                            tlvGenPrAsPi.SetVectM(tlvGenPr.Vect(),MPION);
                            tlvGenPiAsPr.SetVectM(tlvGenPi.Vect(),MPROTON);
                            const TLorentzVector tlvGenLambdaSwapped=tlvGenPrAsPi+tlvGenPiAsPr;
                            fgml0sw = tlvGenLambdaSwapped.M();
                            // make a Lb
                            tlvGenLambdaB = tlvGenLambda0+tlvGenJp;
                            fgmlb = tlvGenLambdaB.M();
                            const TLorentzVector tlvGenLambdaBSwapped = tlvGenLambdaSwapped+tlvGenJp;
                            fgmlbsw = tlvGenLambdaBSwapped.M();

                            //
                            //if (fgml0>2) fpEvt->dumpGenBlock();
                        }
                    }
                    fgdRl0lb = tlvGenJp.DeltaR(tlvGenLambda0);
                    fganl0jp = tlvGenJp.Angle(tlvGenLambda0.Vect());
                    fganl0lb = tlvGenLambdaB.Angle(tlvGenLambda0.Vect());
                    const double anl0mu1 = tlvGenLambda0.Angle(tlvGenMu1.Vect());
                    const double anl0mu2 = tlvGenLambda0.Angle(tlvGenMu2.Vect());
                    fganl0mumin = anl0mu1<anl0mu2 ? anl0mu1 : anl0mu2;
                    fganl0muPt = tlvGenMu1.Pt() > tlvGenMu2.Pt() ? anl0mu1 : anl0mu2;
                }
            }
        }
    }
}

void lambdaReader::doL1stuff()
{
    fL1TDMu0 = fL1TDMu3 = false;
    fL1TMuBH = fL1TMu0 = fL1TMu3 = fL1TMu5 = fL1TMu7 = fL1TMu10 = fL1TMu14 = fL1TMu20 = false;
    for(int i=0; i!=NL1T; i++)
    {
	const std::string name = fpEvt->fL1TNames[i].Data();
	//const double prescale = fpEvt->fL1TPrescale[i];
	const bool result = fpEvt->fL1TResult[i];
	const bool error = fpEvt->fL1TError[i];
	//cout << name << " - prescale: " << prescale << " result: " << result << " error: " << error << endl;
	if (result && !error)
	{
	    if ("L1_DoubleMuOpen"== name) fL1TDMu0 = true;
	    if ("L1_DoubleMu3"== name) fL1TDMu3 = true;
	    if ("L1_SingleMuBeamHalo"== name) fL1TMuBH = true;
	    if ("L1_SingleMu0"== name) fL1TMu0 = true;
	    if ("L1_SingleMu3"== name) fL1TMu3 = true;
	    if ("L1_SingleMu5"== name) fL1TMu5 = true;
	    if ("L1_SingleMu7"== name) fL1TMu7 = true;
	    if ("L1_SingleMu10"== name) fL1TMu10 = true;
	    if ("L1_SingleMu14"== name) fL1TMu14 = true;
	    if ("L1_SingleMu20"== name) fL1TMu20 = true;
	}
    }
}

void lambdaReader::doHLTstuff()
{
    fHLTqrk = fHLTqrkLS = false;
    fHLTDMuOp = fHLTDMu0 = fHLTDMu3 = fHLTDMu5 = false;
    fHLTMu3t3jp = fHLTMu3t5jp = fHLTMu5t0jp = false;
    fHLTMu0jp = fHLTMu0jpT = fHLTMu3jp = fHLTMu3jpT = fHLTMu5jp = fHLTMu5jpT = false;
    fHLTL1DMu0 = fHLTL2DMu0 = fHLTL2Mu0 = false;
    for(int i=0; i!=NHLT; i++)
    {
	const std::string name = fpEvt->fHLTNames[i].Data();
	if (name == "") break; // if name is empty we reached the end of the list (according to Urs)
	const double prescale = fpEvt->fHLTPrescale[i];
	const bool wasrun = fpEvt->fHLTWasRun[i];
	const bool result = fpEvt->fHLTResult[i];
	const bool error = fpEvt->fHLTError[i];
	//cout << name << " - wasrun: " << wasrun << " prescale: " << prescale << " result: " << result << " error: " << error << endl;
	if (result && wasrun &&!error)
	{
	    // Quarkonium triggers
	    if ("HLT_DoubleMu0_Quarkonium_v1" == name)    fHLTqrk = true;
	    if ("HLT_DoubleMu0_Quarkonium_LS_v1" == name) fHLTqrkLS = true;
	    // DoubleMu triggers (usually prescaled)
	    if ("HLT_DoubleMuOpen"       == name) fHLTDMuOp = true;
	    if ("HLT_DoubleMu0"          == name) fHLTDMu0 = true;
	    if ("HLT_DoubleMu3_v2"       == name) fHLTDMu3 = true;
	    if ("HLT_DoubleMu5_v1"       == name) fHLTDMu5 = true;
	    // Jpsi triggers
	    if ("HLT_Mu3_Track3_Jpsi"    == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track3_Jpsi_v2" == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track3_Jpsi_v3" == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v1" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v2" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v3" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu5_Track0_Jpsi"    == name) fHLTMu5t0jp = true;
	    if ("HLT_Mu5_Track0_Jpsi_v2" == name) fHLTMu5t0jp = true;
	    // Jpsi with OST
	    if ("HLT_Mu0_TkMu0_OST_Jpsi" == name) fHLTMu0jp = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi" == name) fHLTMu3jp = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu3jpT = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3" == name) fHLTMu3jpT = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi" == name) fHLTMu5jp = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1" == name) fHLTMu5jpT = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu5jpT = true;
	    // L1 passthrough triggers
	    if ("HLT_L1DoubleMuOpen"     == name) fHLTL1DMu0 = true;
	    if ("HLT_L2DoubleMu0"        == name) fHLTL2DMu0 = true;
	    if ("HLT_Mu5_L2Mu0"          == name) fHLTL2Mu0 = true;
	}
    }
}

// ----------------------------------------------------------------------
void lambdaReader::initVariables()
{

}

// ----------------------------------------------------------------------
void lambdaReader::fillHist()
{

    cout << "fillHist()" << endl;
    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());

    if (0 != fpCand1)
    {
        ((TH1D*)fpHistFile->Get("h10"))->Fill(fpCand1->fPlab.Perp());
        ((TH1D*)fpHistFile->Get("h11"))->Fill(fpCand1->fMass);
        ((TH1D*)fpHistFile->Get("h12"))->Fill(fpCand1->fVtx.fChi2);

        fTree->Fill();
    }

}

// ----------------------------------------------------------------------
void lambdaReader::bookHist()
{
    cout << "==> lambdaReader: bookHist " << endl;

    TH1D *h;
    h = new TH1D("h1", "Ntrk", 500, 0., 1000.);
    h = new TH1D("h2", "NCand", 20, 0., 20.);
    h = new TH1D("h3", "cand ID", 1000100, -100., 1000000.);

    h = new TH1D("h4", ("NCand " + toString(CUTLbCandidate)).c_str(), 20, 0., 20.);
    h = new TH1D("h5", ("NCand fit " + toString(CUTLbCandidateFit)).c_str(), 20, 0., 20.);
    h = new TH1D("h6", ("NCand matched " + toString(CUTLbCandidate) + " to " + toString(CUTLbCandidateFit)).c_str(), 20, 0., 20.);

    h = new TH1D("h10", "pT", 40, 0., 20.);
    h = new TH1D("h11", "mass", 50, 1.6, 2.1);
    h = new TH1D("h12", "chi2", 50, 0., 10.);

    return;
}


// ----------------------------------------------------------------------
void lambdaReader::bookReducedTree()
{
    cout << "==> lambdaReader: bookReducedTree" << endl;

    // create the events tree ======================================
    fTree = new TTree("events", "events");

    // run info
    fTree->Branch("run",     &fRun,     "run/I");
    fTree->Branch("event",   &fEvt,   "event/I");
    fTree->Branch("LS",      &fLS,      "LS/I");

    // candidates
    fTree->Branch("mlb",     &fmlb,     "mlb/D");
    fTree->Branch("ml0",     &fml0,     "ml0/D");
    fTree->Branch("mjp",     &fmjp,     "mjp/D");
    fTree->Branch("ptlb",    &fptlb,    "ptlb/D");
    fTree->Branch("ptl0",    &fptl0,    "ptl0/D");
    fTree->Branch("pjp",     &fpjp,     "pjp/D");
    fTree->Branch("plb",     &fplb,     "plb/D");
    fTree->Branch("pl0",     &fpl0,     "pl0/D");
    fTree->Branch("ptjp",    &fptjp,    "ptjp/D");
    fTree->Branch("etalb",   &fetalb,   "etalb/D");
    fTree->Branch("etal0",   &fetal0,   "etal0/D");
    fTree->Branch("etajp",   &fetajp,   "etajp/D");
    fTree->Branch("philb",   &fphilb,   "philb/D");
    fTree->Branch("phil0",   &fphil0,   "phil0/D");
    fTree->Branch("phijp",   &fphijp,   "phijp/D");

    // signal tracks
    fTree->Branch("rpt1m",   &frpt1m,   "rpt1m/D");
    fTree->Branch("rpt2m",   &frpt2m,   "rpt2m/D");
    fTree->Branch("rptpr",   &frptpr,   "rptpr/D");
    fTree->Branch("rptpi",   &frptpi,   "rptpi/D");
    fTree->Branch("reta1m",  &freta1m,  "reta1m/D");
    fTree->Branch("reta2m",  &freta2m,  "reta2m/D");
    fTree->Branch("retapr",  &fretapr,  "retapr/D");
    fTree->Branch("retapi",  &fretapi,  "retapi/D");
    fTree->Branch("rphi1m",  &frphi1m,  "rphi1m/D");
    fTree->Branch("rphi2m",  &frphi2m,  "rphi2m/D");
    fTree->Branch("rphipr",  &frphipr,  "rphipr/D");
    fTree->Branch("rphipi",  &frphipi,  "rphipi/D");
    fTree->Branch("rid1m",   &frid1m,   "rid1m/I");
    fTree->Branch("rid2m",   &frid2m,   "rid2m/I");
    fTree->Branch("rq1m",    &frq1m,    "rq1m/I");
    fTree->Branch("rq2m",    &frq2m,    "rq2m/I");
    fTree->Branch("rqpr",    &frqpr,    "rqpr/I");
    fTree->Branch("rqpi",    &frqpi,    "rqpi/I");
    fTree->Branch("chi21m",  &fchi21m,  "chi21m/D");
    fTree->Branch("chi22m",  &fchi22m,  "chi22m/D");
    fTree->Branch("chi2pr",  &fchi2pr,  "chi2pr/D");
    fTree->Branch("chi2pi",  &fchi2pi,  "chi2pi/D");
    fTree->Branch("prob1m",  &fprob1m,  "prob1m/D");
    fTree->Branch("prob2m",  &fprob2m,  "prob2m/D");
    fTree->Branch("probpr",  &fprobpr,  "probpr/D");
    fTree->Branch("probpi",  &fprobpi,  "probpi/D");
    fTree->Branch("ndof1m",  &fndof1m,  "ndof1m/I");
    fTree->Branch("ndof2m",  &fndof2m,  "ndof2m/I");
    fTree->Branch("ndofpr",  &fndofpr,  "ndofpr/I");
    fTree->Branch("ndofpi",  &fndofpi,  "ndofpi/I");
    fTree->Branch("qual1m",  &fqual1m,  "qual1m/I");
    fTree->Branch("qual2m",  &fqual2m,  "qual2m/I");
    fTree->Branch("qualpr",  &fqualpr,  "qualpr/I");
    fTree->Branch("qualpi",  &fqualpi,  "qualpi/I");

    fTree->Branch("Kshypo",  &fKshypo,  "Kshypo/D");
    fTree->Branch("alphalb", &falphalb, "alphalb/D");
    fTree->Branch("alphal0", &falphal0, "alphal0/D");
    fTree->Branch("maxdocalb",&fmaxdocalb,    "maxdocalb/D");
    fTree->Branch("maxdocal0",&fmaxdocal0,    "maxdocal0/D");
    fTree->Branch("maxdocajp",&fmaxdocajp,    "maxdocajp/D");
    fTree->Branch("d3lb",    &fd3lb,    "d3lb/D");
    fTree->Branch("d3l0",    &fd3l0,    "d3l0/D");
    fTree->Branch("d3jp",    &fd3jp,    "d3jp/D");
    fTree->Branch("d3Elb",   &fd3Elb,   "d3Elb/D");
    fTree->Branch("d3El0",   &fd3El0,   "d3El0/D");
    fTree->Branch("d3Ejp",   &fd3Ejp,   "d3Ejp/D");
    fTree->Branch("dxylb",   &fdxylb,   "dxylb/D");
    fTree->Branch("dxyl0",   &fdxyl0,   "dxyl0/D");
    fTree->Branch("dxyjp",   &fdxyjp,   "dxyjp/D");
    fTree->Branch("dxyElb",  &fdxyElb,  "dxyElb/D");
    fTree->Branch("dxyEl0",  &fdxyEl0,  "dxyEl0/D");
    fTree->Branch("dxyEjp",  &fdxyEjp,  "dxyEjp/D");

    fTree->Branch("ctlb",    &fctlb,    "ctlb/D");
    fTree->Branch("ctl0",    &fctl0,    "ctl0/D");

    fTree->Branch("btlbx",   &fbtlbx,   "btlbx/D");
    fTree->Branch("btlby",   &fbtlby,   "btlby/D");
    fTree->Branch("btlbz",   &fbtlbz,   "btlbz/D");
    fTree->Branch("btl0x",   &fbtl0x,   "btl0x/D");
    fTree->Branch("btl0y",   &fbtl0y,   "btl0y/D");
    fTree->Branch("btl0z",   &fbtl0z,   "btl0z/D");

    fTree->Branch("vxlb",    &fvxlb,    "vxlb/D");
    fTree->Branch("vylb",    &fvylb,    "vylb/D");
    fTree->Branch("vzlb",    &fvzlb,    "vzlb/D");
    fTree->Branch("vrlb",    &fvrlb,    "vrlb/D");
    fTree->Branch("vxl0",    &fvxl0,    "vxl0/D");
    fTree->Branch("vyl0",    &fvyl0,    "vyl0/D");
    fTree->Branch("vzl0",    &fvzl0,    "vzl0/D");
    fTree->Branch("vrl0",    &fvrl0,    "vrl0/D");

    fTree->Branch("chi2lb",  &fchi2lb,  "chi2lb/D");
    fTree->Branch("chi2l0",  &fchi2l0,  "chi2l0/D");
    fTree->Branch("chi2jp",  &fchi2jp,  "chi2jp/D");
    fTree->Branch("ndoflb",  &fndoflb,  "ndoflb/D");
    fTree->Branch("ndofl0",  &fndofl0,  "ndofl0/D");
    fTree->Branch("ndofjp",  &fndofjp,  "ndofjp/D");
    fTree->Branch("problb",  &fproblb,  "problb/D");
    fTree->Branch("probl0",  &fprobl0,  "probl0/D");
    fTree->Branch("probjp",  &fprobjp,  "probjp/D");

    fTree->Branch("dRprpi",  &fdRprpi,  "dRprpi/D");
    fTree->Branch("dRmumu",  &fdRmumu,  "dRmumu/D");
    fTree->Branch("dRl0jp",  &fdRl0jp,  "dRl0jp/D");

    fTree->Branch("anglbl0",&fanglbl0,"anglbl0/D");
    // gen info to main tree
    fTree->Branch("nDau",    &fnDaughters,"nDau/I");
    fTree->Branch("nGDau",   &fnGrandDaughters,"nGDau/I");
    fTree->Branch("nRef1G",  &fmapRef1Gen,"nRef1G/I");
    fTree->Branch("nRef2G",  &fmapRef2Gen,"nRef2G/I");
    fTree->Branch("isSig",   &fIsSig,     "isSig/I");
    fTree->Branch("isMCmatch",&fIsMCmatch, "isMCmatch/I");

    // sigtrack info
    fTree->Branch("Spt1m",   &fSpt1m,   "Spt1m/D");
    fTree->Branch("Seta1m",  &fSeta1m,  "Seta1m/D");
    fTree->Branch("Sphi1m",  &fSphi1m,  "Sphi1m/D");
    fTree->Branch("Spt2m",   &fSpt2m,   "Spt2m/D");
    fTree->Branch("Seta2m",  &fSeta2m,  "Seta2m/D");
    fTree->Branch("Sphi2m",  &fSphi2m,  "Sphi2m/D");
    fTree->Branch("Sptpr",   &fSptpr,   "Sptpr/D");
    fTree->Branch("Setapr",  &fSetapr,  "Setapr/D");
    fTree->Branch("Sphipr",  &fSphipr,  "Sphipr/D");
    fTree->Branch("Sptpi",   &fSptpi,   "Sptpi/D");
    fTree->Branch("Setapi",  &fSetapi,  "Setapi/D");
    fTree->Branch("Sphipi",  &fSphipi,  "Sphipi/D");

    // Trigger info
    fTree->Branch("L1TDMu0", &fL1TDMu0, "L1TDMu0/O");
    fTree->Branch("L1TDMu3", &fL1TDMu3, "L1TDMu3/O");
    fTree->Branch("L1TMuBH", &fL1TMuBH, "L1TMuBH/O");
    fTree->Branch("L1TMu0",  &fL1TMu0,  "L1TMu0/O");
    fTree->Branch("L1TMu3",  &fL1TMu3,  "L1TMu3/O");
    fTree->Branch("L1TMu5",  &fL1TMu5,  "L1TMu5/O");
    fTree->Branch("L1TMu7",  &fL1TMu7,  "L1TMu7/O");
    fTree->Branch("L1TMu10", &fL1TMu10, "L1TMu10/O");
    fTree->Branch("L1TMu14", &fL1TMu14, "L1TMu14/O");
    fTree->Branch("L1TMu20", &fL1TMu20, "L1TMu20/O");

    fTree->Branch("HLTqrk", &fHLTqrk, "HLTqrk/O");
    fTree->Branch("HLTqrkLS", &fHLTqrkLS, "HLTqrkLS/O");
    fTree->Branch("HLTDMuOp", &fHLTDMuOp, "HLTDMuOp/O");
    fTree->Branch("HLTDMu0", &fHLTDMu0, "HLTDMu0/O");
    fTree->Branch("HLTDMu3", &fHLTDMu3, "HLTDMu3/O");
    fTree->Branch("HLTDMu5", &fHLTDMu5, "HLTDMu5/O");
    fTree->Branch("HLTMu3t3jp", &fHLTMu3t3jp, "HLTMu3t3jp/O");
    fTree->Branch("HLTMu3t5jp", &fHLTMu3t5jp, "HLTMu3t5jp/O");
    fTree->Branch("HLTMu5t0jp", &fHLTMu5t0jp, "HLTMu5t0jp/O");
    fTree->Branch("HLTMu0jp", &fHLTMu0jp, "HLTMu0jp/O");
    fTree->Branch("HLTMu0jpT", &fHLTMu0jpT, "HLTMu0jpT/O");
    fTree->Branch("HLTMu3jp", &fHLTMu3jp, "HLTMu3jp/O");
    fTree->Branch("HLTMu3jpT", &fHLTMu3jpT, "HLTMu3jpT/O");
    fTree->Branch("HLTMu5jp", &fHLTMu5jp, "HLTMu5jp/O");
    fTree->Branch("HLTMu5jpT", &fHLTMu5jpT, "HLTMu5jpT/O");
    fTree->Branch("HLTL1DMu0", &fHLTL1DMu0, "HLTL1DMu0/O");
    fTree->Branch("HLTL2DMu0", &fHLTL2DMu0, "HLTL2DMu0/O");
    fTree->Branch("HLTL2Mu0", &fHLTL2Mu0, "HLTL2Mu0/O");

    // Generator tree ==========================================
    fGenTree = new TTree("genevents", "genevents");
    // run info
    fGenTree->Branch("run",     &fRun,     "run/I");
    fGenTree->Branch("event",   &fEvt,   "event/I");
    fGenTree->Branch("LS",      &fLS,      "LS/I");

    // MC info
    fGenTree->Branch("nDau",    &fgnDaughters,"nDau/I");
    fGenTree->Branch("nGDau",   &fgnGrandDaughters,"nGDau/I");
    fGenTree->Branch("nRef1G",  &fgmapRef1Gen,"nRef1G/I");
    fGenTree->Branch("nRef2G",  &fgmapRef2Gen,"nRef2G/I");
    fGenTree->Branch("hasCand", &fghasCand,"hasCand/I");

    // Lambda_b
    fGenTree->Branch("mlb",     &fgmlb,    "gmlb/D");
    fGenTree->Branch("mlbSwap", &fgmlbsw,  "gmlbsw/D");
    fGenTree->Branch("ml0",     &fgml0,    "gml0/D");
    fGenTree->Branch("ml0Swap", &fgml0sw,  "gml0sw/D");

    fGenTree->Branch("vxlb",    &fgvxlb,   "vxlb/D");
    fGenTree->Branch("vylb",    &fgvylb,   "vylb/D");
    fGenTree->Branch("vzlb",    &fgvzlb,   "vzlb/D");
    fGenTree->Branch("vrlb",    &fgvrlb,   "vrlb/D");
    fGenTree->Branch("ctlb",    &fgctlb,   "ctlb/D");

    fGenTree->Branch("vxl0",    &fgvxl0,   "vxl0/D");
    fGenTree->Branch("vyl0",    &fgvyl0,   "vyl0/D");
    fGenTree->Branch("vzl0",    &fgvzl0,   "vzl0/D");
    fGenTree->Branch("vrl0",    &fgvrl0,   "vrl0/D");
    fGenTree->Branch("ctl0",    &fgctl0,   "ctl0/D");

    // Signal tracks
    fGenTree->Branch("ptmu1",   &fgptmu1,  "ptmu1/D");
    fGenTree->Branch("pmu1",    &fgpmu1,  "pmu1/D");
    fGenTree->Branch("etamu1",  &fgetamu1, "etamu1/D");
    fGenTree->Branch("phimu1",  &fgphimu1, "phimu1/D");

    fGenTree->Branch("ptmu2",   &fgptmu2,  "ptmu2/D");
    fGenTree->Branch("pmu2",    &fgpmu2,  "pmu2/D");
    fGenTree->Branch("etamu2",  &fgetamu2, "etamu2/D");
    fGenTree->Branch("phimu2",  &fgphimu2, "phimu2/D");

    fGenTree->Branch("ptpr",    &fgptpr,   "ptpr/D");
    fGenTree->Branch("ppr",     &fgppr,    "ppr/D");
    fGenTree->Branch("etapr",   &fgetapr,  "etapr/D");
    fGenTree->Branch("phipr",   &fgphipr,  "phipr/D");

    fGenTree->Branch("ptpi",    &fgptpi,   "ptpi/D");
    fGenTree->Branch("ppi",     &fgppi,    "ppi/D");
    fGenTree->Branch("etapi",   &fgetapi,  "etapi/D");
    fGenTree->Branch("phipi",   &fgphipi,  "phipi/D");

    // Lambda_b direct daughters
    fGenTree->Branch("ptl0",    &fgptl0,   "ptl0/D");
    fGenTree->Branch("pl0",     &fgpl0,    "pl0/D");

    fGenTree->Branch("dRprpi",  &fgdRprpi, "dRprpi/D");
    fGenTree->Branch("dRmumu",  &fgdRmumu, "dRmumu/D");
    fGenTree->Branch("dRl0lb",  &fgdRl0lb, "dRl0lb/D");

    // Angles between tracks
    fGenTree->Branch("anprpi",  &fganprpi, "anprpi/D");
    fGenTree->Branch("anmumu",  &fganmumu, "anmumu/D");
    fGenTree->Branch("anl0jp",  &fganl0jp, "anl0jp/D");
    fGenTree->Branch("anl0lb",  &fganl0lb, "anl0lb/D");
    fGenTree->Branch("anl0mumin",&fganl0mumin, "anl0mumin/D");
    fGenTree->Branch("anl0muPt",&fganl0muPt,"anl0muPt/D");

}

// ----------------------------------------------------------------------
void lambdaReader::readCuts(TString filename, int dump)
{
    // set some default cuts
    CUTLbCandidateFit = 0; // indicates non set value
    // read cut file
    fCutFile = filename;
    if (dump) cout << "==> lambdaReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
    ifstream infile(fCutFile.Data());

    TString fn(fCutFile.Data());

    if (dump)
    {
        cout << "====================================" << endl;
        cout << "==> lambdaReader: Cut file  " << fCutFile.Data() << endl;
        cout << "------------------------------------" << endl;
    }

    TH1D *hcuts = new TH1D("hcuts", "", 40, 0., 40.);
    hcuts->GetXaxis()->SetBinLabel(1, fn.Data());

    // loop over lines of input file
    std::string strLine;
    while (std::getline(infile,strLine))
    {
        // strip off everything after a hash or a slash
        strLine = stripOff(strLine, '#');
        strLine = stripOff(strLine, '/');
        std::istringstream iss(strLine, std::istringstream::in);
        // intermediate storage of the entries
        std::string key, value;
        if (iss >> key >> value)
        {
            if (dump) cout << key << ": " << value << endl;
            // a switch statement doesn't work as C++ does not allow to take strings as switch argument
            // http://www.codeguru.com/cpp/cpp/cpp_mfc/article.php/c4067
            if("CUTLbCandidate" == key) // The candidate to retrieve signal values
            {
                setCut(CUTLbCandidate, value);
                //setCut(CUTLbCandidate, value, hcuts, 1, "Candidate");
                continue;
            }
            if("CUTLbCandidateFit" == key) // The candidate to retrieve kinematic vertex fit data
            {
                setCut(CUTLbCandidateFit, value);
                continue;
            }
            if("CUTMuId1" == key)
            {
                setCut(CUTMuId1, value);
                continue;
            }
            if("CUTMuId2" == key)
            {
                setCut(CUTMuId2, value);
                continue;
            }
            if("CUTMjpMin" == key)
            {
                setCut(CUTMjpMin, value, hcuts, 1, "m^{min}(J/#psi)");
                continue;
            }
            if("CUTMjpMax" == key)
            {
                setCut(CUTMjpMax, value, hcuts, 2, "m^{max}(J/#psi)");
                continue;
            }
            if("CUTAlphal0Max" == key)
            {
                setCut(CUTAlphal0Max, value, hcuts, 3, "#alpha^{max}(#Lambda^{b})");
                continue;
            }
            if("CUTMl0Max" == key)
            {
                setCut(CUTMl0Max, value, hcuts, 4, "m^{max}(#Lambda^{0})");
                continue;
            }
            if("CUTD3dl0Min" == key)
            {
                setCut(CUTD3dl0Min, value, hcuts, 5, "d3^{min}(#Lambda^{0})");
                continue;
            }
            if("CUTPtl0Min" == key)
            {
                setCut(CUTPtl0Min, value, hcuts, 6, "p_{#perp}^{min}(#Lambda^{0})");
                continue;
            }
            if("CUTReadDecayMaps" == key)
            {
                setCut(CUTReadDecayMaps, value);
                continue;
            }
            if("CUTPrintDecayMaps" == key)
            {
                setCut(CUTPrintDecayMaps, value);
                continue;
            }
            if("CUTDecayMap1" == key)
            {
                setCut(CUTDecayMap1, value);
                continue;
            }
            if("CUTDecayMap2" == key)
            {
                setCut(CUTDecayMap2, value);
                continue;
            }

            cout << "==> lambdaReader: ERROR in cutfile: Don't know what to do with " << key << "=" << value << endl;
        }
    }

    if (dump)  cout << "------------------------------------" << endl;

}
