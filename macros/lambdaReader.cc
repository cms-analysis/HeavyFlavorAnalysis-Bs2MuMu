#include "lambdaReader.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>
#include <string>
#include <map>
#include <sstream>

using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/lambdaReader.default.cuts
//           ./runTreeReaders -f test.root
// ----------------------------------------------------------------------

template <typename T>
std::string toString(T i)
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

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
    std::cout << "map contains " << myMap.size() << " entries" << std::endl;
    std::cout << "revmap contains " << myRevMap.size() << " entries" << std::endl;
    for(revmap_t::const_iterator it = myRevMap.begin(); it!=myRevMap.end(); it++)
    {
        std::cout << it->first << " - " << it->second << std::endl;
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

// ----------------------------------------------------------------------
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
    myDecayMap1Gen.getPos("443 3122 "); // my daughters
    myDecayMap1Gen.getPos("443 3122 211 211 ");
    myDecayMap1Gen.getPos("443 3122 211 211 22 ");
    myDecayMap1Gen.getPos("443 3122 211 211 22 22 ");
}

void lambdaReader::endAnalysis()
{
    std::cout << "Decay map 1Gen: " << std::endl;
    myDecayMap1Gen.printMap();
    //std::cout << "Decay map 2Gen: " << std::endl;
    //myDecayMap2Gen.printMap();
}

// ----------------------------------------------------------------------
void lambdaReader::eventProcessing()
{

    // fpEvt->dumpGenBlock();
    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());
    ((TH1D*)fpHistFile->Get("h2"))->Fill(fpEvt->nCands());

    TAnaCand *pCand;
    TH1D *h;
    int n1313(0), n1300(0);



    for (int iC = 0; iC < fpEvt->nCands(); ++iC)
    {
        pCand = fpEvt->getCand(iC);
        ((TH1D*)fpHistFile->Get("h3"))->Fill(pCand->fType);
        if (pCand->fType < 100)
        {
            //      cout << Form("%6i %i", iC, pCand->fType) << endl;
        }

        if (h = (TH1D*)fpHistFile->Get(Form("m%d", pCand->fType)))
        {
            if (pCand->fType > 20000 && pCand->fType < 21000)
            {
                h->Fill(pCand->fVar1);
            }
            else
            {
                h->Fill(pCand->fMass);
            }
        }
        else
        {
            //      cout << "Unknown candidate " << pCand->fType << endl;
        }

        //if (100521 == pCand->fType) doBplus(pCand);
        //if (20010 == pCand->fType) doDzero(pCand);
        //if (1300 == pCand->fType) {
        //++n1300;
        //}
        //if (1313 == pCand->fType) {
        if (600443 == pCand->fType)
        {
            TAnaTrack *mu1, *mu2;
            ((TH1D*)fpHistFile->Get("m443_0"))->Fill(pCand->fMass);
            mu1 = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
            mu2 = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex);
            if( (mu1->fMuID & 6)==6)
            {
                ((TH1D*)fpHistFile->Get("m443_1"))->Fill(pCand->fMass);
            }
            if( (mu2->fMuID & 6)==6)
            {
                ((TH1D*)fpHistFile->Get("m443_2"))->Fill(pCand->fMass);
            }
            if( (mu1->fMuID & 6)==6 && (mu2->fMuID & 6)==6)
            {
                ((TH1D*)fpHistFile->Get("m443_3"))->Fill(pCand->fMass);
            }

            //++n1313;
            //doUpsilon(pCand);
        }

        if (cutLbCandidate == pCand->fType)
        {
            TAnaCand *jpsi, *lambda0;
            jpsi = fpEvt->getCand(pCand->fDau1);
            lambda0 = fpEvt->getCand(pCand->fDau2);

            // calculate alpha
            const TVector3 vecVtxLambda0 = lambda0->fVtx.fPoint;
            const TVector3 vecVtxJpsi = jpsi->fVtx.fPoint;
            const TVector3 tmpVect = vecVtxLambda0-vecVtxJpsi;
            //const double alpha = vecVtxLambda0.Angle(vecVtxJpsi);
            const double alphal0 = tmpVect.Angle(lambda0->fPlab);
            const double alphalb = vecVtxLambda0.Angle(vecVtxJpsi);
            //const double vtxDistLambda0(vecVtxLambda0.Mag());
            //const double vtxDistLambdaB(vecVtxJpsi.Mag());

            // get the muons
            TAnaTrack *mu1, *mu2;
            mu1 = fpEvt->getRecTrack(fpEvt->getSigTrack(jpsi->fSig1)->fIndex);
            mu2 = fpEvt->getRecTrack(fpEvt->getSigTrack(jpsi->fSig2)->fIndex);
	    // require the muons to be of at least some quality
	    const int muQual=4;
	    if(mu1->fMuID==-1||(mu1->fMuID & muQual)!=muQual) continue;
	    if(mu2->fMuID==-1||(mu2->fMuID & muQual)!=muQual) continue;

            // get the lambda daughters
            TAnaTrack *pion, *proton;
            pion = fpEvt->getRecTrack(fpEvt->getSigTrack(lambda0->fSig1)->fIndex);
            proton = fpEvt->getRecTrack(fpEvt->getSigTrack(lambda0->fSig2)->fIndex);
	    TAnaTrack *tatSigPi, *tatSigPr;
	    tatSigPi = fpEvt->getSigTrack(lambda0->fSig1);
	    tatSigPr = fpEvt->getSigTrack(lambda0->fSig2);

            const TVector3 tv3L0 = pion->fPlab+proton->fPlab;
            if(lambda0->fMass < 1.12
                //    && TMath::Abs(1-tv3L0.Perp()/lambda0->fPlab.Perp())>.2
              )
            {
                if (TMath::Abs(1-tv3L0.Perp()/lambda0->fPlab.Perp())>.2)
		    std::cout << "------- large discrepancy in pT ";
                std::cout << "----------------------------------------" << std::endl;
                std::cout << pCand->fVtx.fChi2 << "/" << pCand->fVtx.fNdof << "  "
                          << lambda0->fVtx.fChi2 << "/" << lambda0->fVtx.fNdof << "  "
                          << jpsi->fVtx.fChi2 << "/" << jpsi->fVtx.fNdof << std::endl;
                TAnaCand* mypCand;
                std::cout << iC << " "
                          << "fDau: " << pCand->fDau1 << " " << pCand->fDau2 << " "
                          << "tv3L0.Perp(): " << tv3L0.Perp() << " "
                          << "lambda0->fPlab.Perp(): " << lambda0->fPlab.Perp()
                          << std::endl;
                for (int IC = 0; IC < fpEvt->nCands(); ++IC)
                {
                    mypCand = fpEvt->getCand(IC);
                    //if(mypCand->fVtx.fChi2/mypCand->fVtx.fNdof < 10)
                    {
                        std::cout << IC << " " << mypCand->fType << " mypCand->fDau " << mypCand->fDau1 << " " << mypCand->fDau2 << " " ;
                        if(mypCand->fDau2>-1)
                        {
                            std::cout << "fDau1)->fSigN: " << fpEvt->getCand(mypCand->fDau1)->fSig1 << " " << fpEvt->getCand(mypCand->fDau1)->fSig2 << " ";
                            std::cout << "fDau2)->fSigN: " << fpEvt->getCand(mypCand->fDau2)->fSig1 << " " << fpEvt->getCand(mypCand->fDau2)->fSig2 << " ";
                            if(fpEvt->getCand(mypCand->fDau2)->fSig1 > -1)
                                std::cout << "fSig1)->fIndex: " << fpEvt->getSigTrack(fpEvt->getCand(mypCand->fDau2)->fSig1)->fIndex << " ";
                            if(fpEvt->getCand(mypCand->fDau2)->fSig2 > -1)
                                std::cout << "fSig2)->fIndex: " << fpEvt->getSigTrack(fpEvt->getCand(mypCand->fDau2)->fSig2)->fIndex << " ";
                        }
                        else
                        {
                            std::cout << "mypCand->fSig: " << mypCand->fSig1 << " " << mypCand->fSig2 << " ";
                        }
                        std::cout << " : ";
                        //mypCand->dump();
                    }
                }
                std::cout << "                -----" << std::endl;
                for (int iRT = 0; iRT< fpEvt->nRecTracks(); ++iRT)
                {
                    std::cout << iRT << " ";
                    //fpEvt->getRecTrack(iRT)->dump();
                }
            }

            // Histos before cuts
            // -----------------------------------------------------

            // alpha
            ((TH1D*)fpHistFile->Get("alpha3122vor"))->Fill(alphal0);

            // pt Histos
            ((TH1D*)fpHistFile->Get("pt_3122"))->Fill(lambda0->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_5122"))->Fill(pCand->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_443"))->Fill(jpsi->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_2212"))->Fill(proton->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_211"))->Fill(pion->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_13_1"))->Fill(mu1->fPlab.Perp());
            ((TH1D*)fpHistFile->Get("pt_13_2"))->Fill(mu2->fPlab.Perp());

            // maxDoca Histos
            ((TH1D*)fpHistFile->Get("maxDoca_5122"))->Fill(pCand->fMaxDoca);
            ((TH1D*)fpHistFile->Get("maxDoca_3122"))->Fill(lambda0->fMaxDoca);
            ((TH1D*)fpHistFile->Get("maxDoca_443"))->Fill(jpsi->fMaxDoca);

            // d3 Histos
            ((TH1D*)fpHistFile->Get("d3_5122"))->Fill(pCand->fVtx.fD3d);
            ((TH1D*)fpHistFile->Get("d3E_5122"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
            ((TH1D*)fpHistFile->Get("dxy_5122"))->Fill(pCand->fVtx.fDxy);
            ((TH1D*)fpHistFile->Get("dxyE_5122"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
            ((TH1D*)fpHistFile->Get("d3_3122"))->Fill(lambda0->fVtx.fD3d);
            ((TH1D*)fpHistFile->Get("d3E_3122"))->Fill(lambda0->fVtx.fD3d/lambda0->fVtx.fD3dE);
            ((TH1D*)fpHistFile->Get("dxy_3122"))->Fill(lambda0->fVtx.fDxy);
            ((TH1D*)fpHistFile->Get("dxyE_3122"))->Fill(lambda0->fVtx.fDxy/lambda0->fVtx.fDxyE);
            ((TH1D*)fpHistFile->Get("d3_443"))->Fill(jpsi->fVtx.fD3d);
            ((TH1D*)fpHistFile->Get("d3E_443"))->Fill(jpsi->fVtx.fD3d/jpsi->fVtx.fD3dE);
            ((TH1D*)fpHistFile->Get("dxy_443"))->Fill(jpsi->fVtx.fDxy);
            ((TH1D*)fpHistFile->Get("dxyE_443"))->Fill(jpsi->fVtx.fDxy/jpsi->fVtx.fDxyE);

            // vtx chi2 histos
            ((TH1D*)fpHistFile->Get("chi2ndof_5122"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
            ((TH1D*)fpHistFile->Get("chi2ndof_3122"))->Fill(lambda0->fVtx.fChi2/lambda0->fVtx.fNdof);
            ((TH1D*)fpHistFile->Get("chi2ndof_443"))->Fill(jpsi->fVtx.fChi2/jpsi->fVtx.fNdof);
            ((TH1D*)fpHistFile->Get("probChi2_5122"))->Fill(TMath::Prob(pCand->fVtx.fChi2,pCand->fVtx.fNdof));
            ((TH1D*)fpHistFile->Get("probChi2_3122"))->Fill(TMath::Prob(lambda0->fVtx.fChi2,lambda0->fVtx.fNdof));
            ((TH1D*)fpHistFile->Get("probChi2_443"))->Fill(TMath::Prob(jpsi->fVtx.fChi2,jpsi->fVtx.fNdof));

            // get some TLorentzVector for eta and phi (and reuse it later)
            TLorentzVector tlvMu1, tlvMu2, tlvPr, tlvPi;
            const TVector3 tv3Mu1 = mu1->fPlab;
            tlvMu1.SetXYZM(tv3Mu1.x(), tv3Mu1.y(), tv3Mu1.z(), MMUON);
            const TVector3 tv3Mu2 = mu2->fPlab;
            tlvMu2.SetXYZM(tv3Mu2.x(), tv3Mu2.y(), tv3Mu2.z(), MMUON);
            const TVector3 tv3Pr = proton->fPlab;
            tlvPr.SetXYZM(tv3Pr.x(), tv3Pr.y(), tv3Pr.z(), MPROTON);
            const TVector3 tv3Pi = pion->fPlab;
            tlvPi.SetXYZM(tv3Pi.x(), tv3Pi.y(), tv3Pi.z(), MPION);

            TLorentzVector tlvLambda0, tlvLambdaB, tlvJpsi;
            const TVector3 tv3Lambda0 = lambda0->fPlab;
            tlvLambda0.SetXYZM(tv3Lambda0.x(), tv3Lambda0.y(), tv3Lambda0.z(), MLAMBDA_0);
            const TVector3 tv3LambdaB = pCand->fPlab;
            tlvLambdaB.SetXYZM(tv3LambdaB.x(), tv3LambdaB.y(), tv3LambdaB.z(), MLAMBDA_B);
            const TVector3 tv3Jpsi = jpsi->fPlab;
            tlvJpsi.SetXYZM(tv3Jpsi.x(), tv3Jpsi.y(), tv3Jpsi.z(), MJPSI);

	    // sigtracks
	    TLorentzVector tlvSigPr, tlvSigPi, tlvSigMu1, tlvSigMu2;
	    tlvSigPr.SetVectM(tatSigPr->fPlab,MPROTON);
	    tlvSigPi.SetVectM(tatSigPi->fPlab,MPION);

            // Set values for reduced tree
            //fRun
            fmlb=pCand->fMass;
            fml0=lambda0->fMass;
            fmjp=jpsi->fMass;
            fptlb=pCand->fPlab.Perp();
            fptl0=lambda0->fPlab.Perp();
            fptjp=jpsi->fPlab.Perp();
            fetalb=tlvLambdaB.Eta();
            fetal0=tlvLambda0.Eta();
            fetajp=tlvJpsi.Eta();
            fphilb=tlvLambdaB.Phi();
            fphil0=tlvLambda0.Phi();
            fphijp=tlvJpsi.Phi();

            fpt1m=mu1->fPlab.Perp();
            fpt2m=mu2->fPlab.Perp();
            fptpr=proton->fPlab.Perp();
            fptpi=pion->fPlab.Perp();
            feta1m=tlvMu1.Eta();
            feta2m=tlvMu2.Eta();
            fetapr=tlvPr.Eta();
            fetapi=tlvPi.Eta();
            fphi1m=tlvMu1.Phi();
            fphi2m=tlvMu2.Phi();
            fphipr=tlvPr.Phi();
            fphipi=tlvPi.Phi();
            fq1m=mu1->fQ;
            fq2m=mu2->fQ;
            fqpr=proton->fQ;
            fqpi=pion->fQ;
            if(mu1->fMuID==-1) fid1m=0;
            else fid1m=mu1->fMuID;
            if(mu2->fMuID==-1) fid2m=0;
            else fid2m=mu2->fMuID;
            falphalb=alphalb;
            falphal0=alphal0;
            fmaxdocalb=pCand->fMaxDoca;
            fmaxdocal0=lambda0->fMaxDoca;
            fmaxdocajp=jpsi->fMaxDoca;

            fd3lb=pCand->fVtx.fD3d;
            fd3l0=lambda0->fVtx.fD3d;
            fd3jp=jpsi->fVtx.fD3d;
            fd3Elb=pCand->fVtx.fD3dE;
            fd3El0=lambda0->fVtx.fD3dE;
            fd3Ejp=jpsi->fVtx.fD3dE;
            fdxylb=pCand->fVtx.fDxy;
            fdxyl0=lambda0->fVtx.fDxy;
            fdxyjp=jpsi->fVtx.fDxy;
            fdxyElb=pCand->fVtx.fDxyE;
            fdxyEl0=lambda0->fVtx.fDxyE;
            fdxyEjp=jpsi->fVtx.fDxyE;

            fchi2lb=pCand->fVtx.fChi2;
            fchi2l0=lambda0->fVtx.fChi2;
            fchi2jp=jpsi->fVtx.fChi2;
            fndoflb=pCand->fVtx.fNdof;
            fndofl0=lambda0->fVtx.fNdof;
            fndofjp=jpsi->fVtx.fNdof;

            fdRprpi=tlvPr.DeltaR(tlvPi);
            fdRmumu=tlvMu1.DeltaR(tlvMu2);
            fdRl0jp=tlvLambda0.DeltaR(tlvLambdaB);

            fanglbl0=tlvLambda0.Angle(tlvLambdaB.Vect());

            // K_s hypothesis
            TLorentzVector tlvPasPiHypo, tlvKsHypoth;
            tlvPasPiHypo.SetXYZM(tv3Pr.x(), tv3Pr.y(), tv3Pr.z(), MPION);
            tlvKsHypoth = tlvPi+tlvPasPiHypo;
            fKshypo=tlvKsHypoth.M(); // fill
            ((TH1D*)fpHistFile->Get("mKsHypo"))->Fill(tlvKsHypoth.M());
            ((TH1D*)fpHistFile->Get("mKsHypo_wide"))->Fill(tlvKsHypoth.M());

	    // sigtracks
	    fSgptpr=tlvSigPr.Perp();
	    fSgetapr=tlvSigPr.Eta();
	    fSgphipr=tlvSigPr.Phi();
	    fSgptpi=tlvSigPi.Perp();
	    fSgetapi=tlvSigPi.Eta();
	    fSgphipi=tlvSigPi.Phi();

            // Do some generator level stuff
            {
                const int nGenCands(fpEvt->nGenCands());
                for(int gcit=0; gcit!=nGenCands; gcit++)
                {
                    TGenCand *gc5122 = fpEvt->getGenCand(gcit);
		    if (511==abs(gc5122->fID) || 521==abs(gc5122->fID) || 513==abs(gc5122->fID) || 523==abs(gc5122->fID))
		    {
			std::cout << "B found: ";
			findDaughters fd(gc5122,fpEvt);
			std::cout << fd.getDauGrdau() << " (" << fd.getNDaughters() << " daughters)" << std::endl;
		    }
		    if(5122==abs(gc5122->fID))
                    {
                        if(gc5122 == gcPrev) continue; // the info should be persistent for the main red. tree
                        gcPrev=gc5122;
                        // reset tree variables
                        fgmlb = fgmlbsw = fgml0 = fgml0sw = 9999;
                        fgptpr = fgptpi = fgppr = fgppi = 9999;
                        fgptl0 = fgpl0 = 9999;
                        fgdRprpi = 9999;
                        fgnDaughters = fgnGrandDaughters = fgmapRef1Gen = -1;
                        //fgnDaughters = fgnGrandDaughters = fgmapRef1Gen = fgmapRef2Gen = -1;
                        // analyse the decay and extract some data to store in the tree
                        findDaughters fd(gc5122,fpEvt);
			std::cout << "L_b found: " << fd.getDauGrdau() << " (" << fd.getNDaughters() << " daughters)" << std::endl;
                        fIsSig = ("443 (13 13 ) 3122 (211 2212 ) "==fd.getDauGrdau()) ? 1 : 0;
                        fnDaughters = fgnDaughters = fd.getNDaughters();
                        fnGrandDaughters = fgnGrandDaughters = fd.getNDauGrdau();
                        fmapRef1Gen = fgmapRef1Gen = myDecayMap1Gen.getPos(fd.getDaughters());
                        //fmapRef2Gen = fgmapRef2Gen = myDecayMap2Gen.getPos(fd.getDauGrdau());
                        // now extract some kinematic data
                        TLorentzVector tlvGenJp, tlvGenL0, tlvGenPr, tlvGenPi;
                        for(int gcDauit=gc5122->fDau1; gcDauit<=gc5122->fDau2; gcDauit++)
                        {
                            TGenCand* gcDau=fpEvt->getGenCand(gcDauit);
                            if(443==abs(gcDau->fID))
                            {
                                tlvGenJp=gcDau->fP;
                            }
                            if(3122==abs(gcDau->fID))
                            {
                                tlvGenL0=gcDau->fP;
                                bool prFound(false), piFound(false);
                                for(int gcGrDauit=gcDau->fDau1; gcGrDauit<=gcDau->fDau2; gcGrDauit++)
                                {
                                    TGenCand* gcGrDau=fpEvt->getGenCand(gcGrDauit);
                                    if(!prFound&&2212==abs(gcGrDau->fID)) // in case we have more than one pr we get the first
                                    {
                                        tlvGenPr=gcGrDau->fP;
                                        prFound=true;
                                    }
                                    if(!piFound&&211==abs(gcGrDau->fID)) // in case we have more than one pi we get the first
                                    {
                                        tlvGenPi=gcGrDau->fP;
                                        piFound=true;
                                    }
                                }
                                if(prFound&&piFound) // we have a pr and a pi
                                {
                                    fgptpr = tlvGenPr.Pt();
                                    fgppr = tlvGenPr.P();
                                    fgptpi = tlvGenPi.Pt();
                                    fgppi = tlvGenPi.P();
                                    fgdRprpi = tlvGenPr.DeltaR(tlvGenPi);
                                    // make a L0
                                    const TLorentzVector tlvGenLambda0 = tlvGenPr + tlvGenPi;
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
                                    const TLorentzVector tlvGenLambdaB = tlvGenLambda0+tlvGenJp;
                                    fgmlb = tlvGenLambdaB.M();
                                    const TLorentzVector tlvGenLambdaBSwapped = tlvGenLambdaSwapped+tlvGenJp;
                                    fgmlbsw = tlvGenLambdaBSwapped.M();

                                    //
                                    //if (fgml0>2) fpEvt->dumpGenBlock();
                                }
                            }
                        }
                        // fill the genTree
                        fGenTree->Fill();
                    }
                }
                //if(found5122==2) fpEvt->dumpGenBlock();
            }

            // Fill the trees
            fTree->Fill();

            // Apply cuts
            if ((mu1->fMuID & cutMuId1)==cutMuId1
                    && (mu2->fMuID & cutMuId2)==cutMuId2
                    && jpsi->fMass > cutMjpMin
                    && jpsi->fMass < cutMjpMax
                    && alphal0 < cutAlphal0Max
                    && lambda0->fMass < cutMl0Max
                    && lambda0->fVtx.fD3d > cutD3dl0Min
                    && lambda0->fPlab.Perp() > cutPtl0Min )
            {
                ((TH1D*)fpHistFile->Get("m5122"))->Fill(pCand->fMass);
                ((TH1D*)fpHistFile->Get("m5122_wide"))->Fill(pCand->fMass);
                ((TH1D*)fpHistFile->Get("m3122"))->Fill(lambda0->fMass);
                ((TH1D*)fpHistFile->Get("m443"))->Fill(jpsi->fMass);
                ((TH1D*)fpHistFile->Get("alpha3122"))->Fill(alphal0);
                // Truth checking

                // getting tlv's
                TLorentzVector tlvLamda0Add;
                tlvLamda0Add=tlvPi+tlvPr;
                ((TH1D*)fpHistFile->Get("tm_m3122"))->Fill(lambda0->fMass);
                ((TH1D*)fpHistFile->Get("tm_m3122_tlvadd"))->Fill(tlvLamda0Add.M());


                const int giMu1(mu1->fGenIndex);
                const int giMu2(mu2->fGenIndex);
                const int giPion(pion->fGenIndex);
                const int giProton(proton->fGenIndex);
                const int nGenCands(fpEvt->nGenCands());
                // do truth matching
                /*
                if( giMu1 < nGenCands && giMu2 < nGenCands && giPion < nGenCands && giProton < nGenCands &&
                    giMu1 > -1 && giMu2 > -1 && giPion > -1 && giProton > -1 )
                {
                    TGenCand *gcMu1 = fpEvt->getGenCand(giMu1);
                    TGenCand *gcMu2 = fpEvt->getGenCand(giMu2);
                    TGenCand *gcPion = fpEvt->getGenCand(giPion);
                    TGenCand *gcProton = fpEvt->getGenCand(giProton);

                    std::cout << "Mu1: " << gcMu1->fID
                	<< " Mu2: " << gcMu2->fID
                	<< " Pion: " << gcPion->fID << " " << pion->fPlab.Perp() << " " << gcPion->fMom1
                	<< " Proton: " << gcProton->fID << " " << proton->fPlab.Perp() << " " << gcProton->fMom1;
                    const bool truthMatched = abs(gcMu1->fID) == 13 && abs(gcMu2->fID) == 13 && abs(gcPion->fID) == 211 && abs(gcProton->fID) == 2212;
                    std::cout << " -- " << truthMatched << std::endl;
                    if(truthMatched)
                    {
                	TLorentzVector tlvGenLamdba0Add, tlvGenJpsiAdd;
                	tlvGenLamdba0Add=gcPion->fP+gcProton->fP;
                	tlvGenJpsiAdd=gcMu1->fP+gcMu2->fP;

                    	((TH1D*)fpHistFile->Get("tm_alpha"))->Fill(alpha);
                	((TH1D*)fpHistFile->Get("tm_m3122"))->Fill(lambda0->fMass);
                	((TH1D*)fpHistFile->Get("tm_m3122_tlvadd"))->Fill(tlvLamda0Add.M());
                	((TH1D*)fpHistFile->Get("gen_m3122_tlvadd"))->Fill(tlvGenLamdba0Add.M());
                	((TH1D*)fpHistFile->Get("gen_m443_tlvadd"))->Fill(tlvGenJpsiAdd.M());
                	((TH1D*)fpHistFile->Get("tm_d3_3122"))->Fill(lambda0->fVtx.fD3d);
                	((TH1D*)fpHistFile->Get("tm_d3E_3122"))->Fill(lambda0->fVtx.fD3d/lambda0->fVtx.fD3dE);
                	((TH1D*)fpHistFile->Get("tm_dxy_3122"))->Fill(lambda0->fVtx.fDxy);
                	((TH1D*)fpHistFile->Get("tm_dxyE_3122"))->Fill(lambda0->fVtx.fDxy/lambda0->fVtx.fDxyE);
                	((TH1D*)fpHistFile->Get("tm_pt_3122"))->Fill(lambda0->fPlab.Perp());
                    }
                }
                */
            }
        }
    }

    ((TH1D*)fpHistFile->Get("h100"))->Fill(n1313);
    ((TH1D*)fpHistFile->Get("h101"))->Fill(n1300);

    return;

    // -- initialize all variables
    initVariables();

    if (!goodRun())
    {
        //cout << "not a good run: " << fRun << endl;
        return;
    }

    // -- track selection for all candidates
    trackSelection();

    // -- Select a candidate
    candidateSelection(0);

//   if (0 != fpCand) {
//     fillHist();
//   }

}


// ----------------------------------------------------------------------
void lambdaReader::initVariables()
{



}



// ----------------------------------------------------------------------
void lambdaReader::doBplus(TAnaCand *pCand )
{

    ((TH1D*)fpHistFile->Get("m100521h0"))->Fill(pCand->fMass);

    if (pCand->fVtx.fChi2 > 5.) return;
    ((TH1D*)fpHistFile->Get("m100521h1"))->Fill(pCand->fMass);

    ((TH1D*)fpHistFile->Get("m100521h10"))->Fill(pCand->fPlab.Perp());

    if (pCand->fPlab.Perp() < 3.) return;
    ((TH1D*)fpHistFile->Get("m100521h2"))->Fill(pCand->fMass);

    if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 1) return;
    ((TH1D*)fpHistFile->Get("m100521h3"))->Fill(pCand->fMass);


}


// ----------------------------------------------------------------------
void lambdaReader::doDzero(TAnaCand *pCand )
{

    TAnaTrack *pt1, *pt2, *pt3;
    //  cout << pCand->fSig1 << " .. " << pCand->fSig2 << endl;
    pt1 = fpEvt->getSigTrack(pCand->fSig1);
    pt2 = fpEvt->getSigTrack(pCand->fSig1+1);
    pt3 = fpEvt->getSigTrack(pCand->fSig2);

    ((TH1D*)fpHistFile->Get("m20010h0"))->Fill(pCand->fVar1);
    if (pCand->fMaxDoca > 0.01) return;
    ((TH1D*)fpHistFile->Get("m20010h1"))->Fill(pCand->fVar1);

    if (pCand->fPlab.Perp() < 4) return;
    ((TH1D*)fpHistFile->Get("m20010h2"))->Fill(pCand->fVar1);
}



// ----------------------------------------------------------------------
void lambdaReader::candidateSelection(int mode)
{

    TAnaCand *pCand(0);
    vector<int> lCands;
    for (int iC = 0; iC < fpEvt->nCands(); ++iC)
    {
        pCand = fpEvt->getCand(iC);
        if (MCTYPE == pCand->fType)
        {
            fillTMCand(pCand, 50);
        }

        if (TYPE != pCand->fType) continue;

        // -- check that candidate fulfilled track selection
        TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1);
        TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1);
        if (0 == pc1->fInt1) continue;
        if (0 == pc2->fInt1) continue;

        if (pCand->fVtx.fChi2 > VTXCHI2) continue;
        if (pCand->fPlab.Perp() < CHARMPTLO) continue;

        lCands.push_back(iC);
    }

    int nc(lCands.size());

    ((TH1D*)fpHistFile->Get("h2"))->Fill(nc);
    if (0 == nc) return;

    int best(0);
    if (nc > 1)
    {
        double ptMax(-1), pt(0.);
        for (unsigned int iC = 0; iC < lCands.size(); ++iC)
        {
            pCand = fpEvt->getCand(lCands[iC]);

            pt = pCand->fPlab.Perp();

            if (0 == mode)
            {
                if (pt > ptMax)
                {
                    ptMax = pt;
                    best = lCands[iC];
                }
            }
        }
    }

    if (best > -1)
    {
        fpCand1 = fpEvt->getCand(best);

        //     TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1);
        //     TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1);

        //     TLorentzVector a1, a2, a0;
        //     a1.SetVectM(pc1->fPlab, MKAON);
        //     a2.SetVectM(pc2->fPlab, MPION);
        //     a0 = a1 + a2;

        //     fCandMass = a0.M();

    }
}


// ----------------------------------------------------------------------
void lambdaReader::fillTMCand(TAnaCand *pCand, int type)
{

    // -- fill simple TM mass
    ((TH1D*)fpHistFile->Get(Form("h%i", 1000+type)))->Fill(pCand->fMass);

    // -- now try to find an additional muon
    TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1);
    TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1);

    int pgI = pc1->fGenIndex;
    TGenCand *pB;
    int foundB(0);
    if (pgI < fpEvt->nGenCands())
    {
        TGenCand *pg = fpEvt->getGenCand(pgI);
        int momI = pg->fMom1;
        TGenCand *pD0 = fpEvt->getGenCand(momI);

        pB = pD0;
        int id(999), cnt(0);
        while (id > 100)
        {
            momI = pB->fMom1;
            pB = fpEvt->getGenCand(momI);
            id  = TMath::Abs(pB->fID%1000);
            ++cnt;
            if (id > 499 && id < 599)
            {
                foundB = 1;
                break;
            }
        }

        if (1 == foundB)
        {
            cout << "========> Found a B" << endl;
            pB->dump();
            TGenCand *pG;
            for (int ig = pB->fDau1; ig <= pB->fDau2; ++ig)
            {
                pG = fpEvt->getGenCand(ig);
                pG->dump();
                if (13 == TMath::Abs(pG->fID))
                {
                    cout << "++++++++++++++++++++++++++++ with a MUON" << endl;
                }
            }
        }

    }

}


// ----------------------------------------------------------------------
void lambdaReader::MCKinematics()
{

}

// ----------------------------------------------------------------------
void lambdaReader::L1TSelection()
{

}

// ----------------------------------------------------------------------
void lambdaReader::HLTSelection()
{

}

// ----------------------------------------------------------------------
void lambdaReader::trackSelection()
{

    TAnaCand *pCand;
    TAnaTrack *pt, *ps[2];

    TLorentzVector pb, pm1, pm2;

    for (int iC = 0; iC < fpEvt->nCands(); ++iC)
    {
        pCand = fpEvt->getCand(iC);
        if (TYPE != pCand->fType) continue;

        // -- Get the 2 signal tracks
        ps[0] = fpEvt->getSigTrack(pCand->fSig1);
        ps[1] = fpEvt->getSigTrack(pCand->fSig2);

        for (int i = 0; i < 2; ++i)
        {

            // -- Get the corresponding RecTrack
            pt = fpEvt->getRecTrack(ps[i]->fIndex);

            // -- Use spare variables in the RecTrack as bookmark whether it passed the signal selection
            ps[i]->fInt1 = 1;

            if (0 == i && pt->fPlab.Perp() < KAPTLO) ps[i]->fInt1 = 0;
            if (1 == i && pt->fPlab.Perp() < PIPTLO) ps[i]->fInt1 = 0;
        }

    }

}


// ----------------------------------------------------------------------
void lambdaReader::muonSelection()
{

}


// ----------------------------------------------------------------------
void lambdaReader::fillHist()
{

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
    h = new TH1D("h100", "1313 multiplicity", 20, 0., 20.);
    h = new TH1D("h101", "1300 multiplicity", 100, 0., 100.);
    h = new TH1D("h10", "pT", 40, 0., 20.);
    h = new TH1D("h11", "mass", 50, 1.6, 2.1);
    h = new TH1D("h12", "chi2", 50, 0., 10.);

    h = new TH1D("h1050", "TM D0->Kpi", 50, 1.6, 2.1);
    h = new TH1D("h2050", "TM mu D0->Kpi", 50, 1.6, 2.1);

    h = new TH1D("m1300", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1313", "mass 1313", 40, 2.7, 3.5);
    h = new TH1D("m20010", "mass 20010", 70, 1.7, 2.4);
    h = new TH1D("m20020", "mass 20020", 70, 1.7, 2.4);
    h = new TH1D("m20030", "mass 20030", 70, 1.7, 2.4);
    h = new TH1D("m20040", "mass 20040", 70, 1.7, 2.4);
    h = new TH1D("m20050", "mass 20050", 70, 1.7, 2.4);
    h = new TH1D("m20060", "mass 20060", 70, 1.7, 2.4);

    h = new TH1D("ups1313h0", "mass 1313", 40, 9.0, 11.0);
    h = new TH1D("ups1313h1", "mass 1313", 40, 9.0, 11.0);
    h = new TH1D("ups1313h2", "mass 1313", 40, 9.0, 11.0);

    // Frank plots

    h = new TH1D("m443_0", "mass 443 - alle Muonen", 40, 2.7, 3.5);
    h = new TH1D("m443_1", "mass 443 - muon1==6", 40, 2.7, 3.5);
    h = new TH1D("m443_2", "mass 443 - muon2==6", 40, 2.7, 3.5);
    h = new TH1D("m443_3", "mass 443 - muonen==6", 40, 2.7, 3.5);

    h = new TH1D("m443", "mass 443", 40, 2.7, 3.5);
    h = new TH1D("m5122", "mass 5122 #Lambda_{b}", 40, 5.4, 5.8);
    h = new TH1D("m5122_wide", "mass 5122 #Lambda_{b}", 200, 4.0, 6.0);
    //h = new TH1D("d5122", "vertex distance 5122 #Lambda_{b}", 40, 0, 10);
    h = new TH1D("alpha3122", "alpha 3122 #Lambda", 80, 0.0, 0.002);
    h = new TH1D("alpha3122vor", "alpha 3122 #Lambda", 80, 0.0, 0.02);
    h = new TH1D("m3122", "mass 3122 #Lambda", 40, 1.0, 1.3);
    //h = new TH1D("d3122", "vertex distance 3122 #Lambda", 40, 0, 10);

    h = new TH1D("tm_alpha", "truth matched -- alpha 3122 #Lambda", 40, 0, 0.1);
    h = new TH1D("tm_m3122", "truth matched -- mass 3122 #Lambda", 40, 1.0, 1.5);
    h = new TH1D("tm_m3122_tlvadd", "truth matched -- mass 3122 #Lambda from tlv addition", 40, 1.0, 1.5);
    h = new TH1D("gen_m3122_tlvadd", "generator -- mass 3122 #Lambda from tlv addition", 40, 1.0, 1.5);
    h = new TH1D("gen_m443_tlvadd", "generator -- mass 443 J/#psi from tlv addition", 40, 2.8, 3.4);

    h = new TH1D("mKsHypo", "K_{s} hypotheses on p#pi", 80, 0.45, 0.55);
    h = new TH1D("mKsHypo_wide", "K_{s} hypotheses on p#pi", 80, 0.2, 2.0);

    h = new TH1D("pt_3122", "p_{T} 3122 #Lambda", 40, 0.0, 20.);
    h = new TH1D("pt_5122", "p_{T} 3122 #Lambda_{b}", 40, 0.0, 20.);
    h = new TH1D("pt_443", "p_{T} 3122 J/#psi", 40, 0.0, 20.);
    h = new TH1D("pt_2212", "p_{T} 3122 p", 40, 0.0, 20.);
    h = new TH1D("pt_211", "p_{T} 3122 #pi", 40, 0.0, 20.);
    h = new TH1D("pt_13_1", "p_{T} 3122 #mu_{1}", 40, 0.0, 20.);
    h = new TH1D("pt_13_2", "p_{T} 3122 #mu_{2}", 40, 0.0, 20.);

    h = new TH1D("maxDoca_5122", "maxDoca 5122 #Lambda_{b}", 40, 0.0, 0.2);
    h = new TH1D("maxDoca_3122", "maxDoca 3122 #Lambda", 40, 0.0, 0.2);
    h = new TH1D("maxDoca_443", "maxDoca 443 J/#psi", 40, 0.0, 0.2);

    h = new TH1D("d3_5122", "d3d 5122 #Lambda_{b}", 80, 0.0, 120.0);
    h = new TH1D("d3E_5122", "d3d/d3dE 5122 #Lambda_{b}", 40, 0.0, 200);
    h = new TH1D("dxy_5122", "dxy 5122 #Lambda_{b}", 80, 0.0, 40.0);
    h = new TH1D("dxyE_5122", "dxy/dxyE 5122 #Lambda_{b}", 40, 0.0, 200);
    h = new TH1D("d3_3122", "d3d 3122 #Lambda", 80, 0.0, 120.0);
    h = new TH1D("d3E_3122", "d3d/d3dE 3122 #Lambda", 40, 0.0, 200);
    h = new TH1D("dxy_3122", "dxy 3122 #Lambda", 80, 0.0, 40.0);
    h = new TH1D("dxyE_3122", "dxy/dxyE 3122 #Lambda", 40, 0.0, 200);
    h = new TH1D("d3_443", "d3d 443 J/#psi", 80, 0.0, 120.0);
    h = new TH1D("d3E_443", "d3d/d3dE 443 J/#psi", 40, 0.0, 200);
    h = new TH1D("dxy_443", "dxy 443 J/#psi", 80, 0.0, 40.0);
    h = new TH1D("dxyE_443", "dxy/dxyE 443 J/#psi", 40, 0.0, 200);

    h = new TH1D("chi2ndof_5122", "#chi^{2}/ndof 5122 #Lambda_{b}", 40, 0.0, 20);
    h = new TH1D("chi2ndof_3122", "#chi^{2}/ndof 3122 #Lambda", 40, 0.0, 20);
    h = new TH1D("chi2ndof_443", "#chi^{2}/ndof 443 J/#psi", 40, 0.0, 20);
    h = new TH1D("probChi2_5122", "Prob(#chi^{2},ndof) 5122 #Lambda_{b}", 40, 0.0, 1.0);
    h = new TH1D("probChi2_3122", "Prob(#chi^{2},ndof) 3122 #Lambda", 40, 0.0, 1.0);
    h = new TH1D("probChi2_443", "Prob(#chi^{2},ndof) 443 J/#psi", 40, 0.0, 1.0);

    //h = new TH1D("tm_m5122", "truth matched -- mass 3122 #Lambda_{b}", 40, 5.0, 6.2);

    // Ende Frankplots

    h = new TH1D("m100521h0",  "mass 100521", 50, 5.0, 5.5);
    h = new TH1D("m100521h1",  "mass 100521", 60, 5.0, 5.6);
    h = new TH1D("m100521h2",  "mass 100521", 60, 5.0, 5.6);
    h = new TH1D("m100521h3",  "mass 100521", 60, 5.0, 5.6);

    h = new TH1D("m100521h10",  "pT 100521", 40, 0., 10.);

    h = new TH1D("m20010h0",  "mass 100521", 60, 1.7, 2.0);
    h = new TH1D("m20010h1",  "mass 100521", 60, 1.7, 2.0);
    h = new TH1D("m20010h2",  "mass 100521", 60, 1.7, 2.0);
    h = new TH1D("m20010h3",  "mass 100521", 60, 1.7, 2.0);

    h = new TH1D("m1300h0", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1300h1", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1300h2", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1300h3", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1300h4", "mass 1300", 40, 2.7, 3.5);
    h = new TH1D("m1300h5", "mass 1300", 40, 2.7, 3.5);


    // -- Reduced Tree
    fTree = new TTree("events", "events");
    fTree->Branch("run",     &fRun,     "run/I");
    // candidates
    fTree->Branch("mlb",     &fmlb,     "mlb/D");
    fTree->Branch("ml0",     &fml0,     "ml0/D");
    fTree->Branch("mjp",     &fmjp,     "mjp/D");
    fTree->Branch("ptlb",    &fptlb,    "ptlb/D");
    fTree->Branch("ptl0",    &fptl0,    "ptl0/D");
    fTree->Branch("ptjp",    &fptjp,    "ptjp/D");
    fTree->Branch("etalb",   &fetalb,   "etalb/D");
    fTree->Branch("etal0",   &fetal0,   "etal0/D");
    fTree->Branch("etajp",   &fetajp,   "etajp/D");
    fTree->Branch("philb",   &fphilb,   "philb/D");
    fTree->Branch("phil0",   &fphil0,   "phil0/D");
    fTree->Branch("phijp",   &fphijp,   "phijp/D");
    // signal tracks
    fTree->Branch("pt1m",    &fpt1m,    "pt1m/D");
    fTree->Branch("pt2m",    &fpt2m,    "pt2m/D");
    fTree->Branch("ptpr",    &fptpr,    "ptpr/D");
    fTree->Branch("ptpi",    &fptpi,    "ptpi/D");
    fTree->Branch("eta1m",   &feta1m,   "eta1m/D");
    fTree->Branch("eta2m",   &feta2m,   "eta2m/D");
    fTree->Branch("etapr",   &fetapr,   "etapr/D");
    fTree->Branch("etapi",   &fetapi,   "etapi/D");
    fTree->Branch("phi1m",   &fphi1m,   "phi1m/D");
    fTree->Branch("phi2m",   &fphi2m,   "phi2m/D");
    fTree->Branch("phipr",   &fphipr,   "phipr/D");
    fTree->Branch("phipi",   &fphipi,   "phipi/D");
    fTree->Branch("id1m",    &fid1m,    "id1m/I");
    fTree->Branch("id2m",    &fid2m,    "id2m/I");
    fTree->Branch("q1m",     &fq1m,     "q1m/I");
    fTree->Branch("q2m",     &fq2m,     "q2m/I");
    fTree->Branch("qpr",     &fqpr,     "qpr/I");
    fTree->Branch("qpi",     &fqpi,     "qpi/I");

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
    fTree->Branch("chi2lb",  &fchi2lb,  "chi2lb/D");
    fTree->Branch("chi2l0",  &fchi2l0,  "chi2l0/D");
    fTree->Branch("chi2jp",  &fchi2jp,  "chi2jp/D");
    fTree->Branch("ndoflb",  &fndoflb,  "ndoflb/D");
    fTree->Branch("ndofl0",  &fndofl0,  "ndofl0/D");
    fTree->Branch("ndofjp",  &fndofjp,  "ndofjp/D");

    fTree->Branch("dRprpi",  &fdRprpi,  "dRprpi/D");
    fTree->Branch("dRmumu",  &fdRmumu,  "dRmumu/D");
    fTree->Branch("dRl0jp",  &fdRl0jp,  "dRl0jp/D");

    fTree->Branch("anglbl0",&fanglbl0,"anglbl0/D");
    // gen info to main tree
    fTree->Branch("nDau",    &fnDaughters,"nDau/I");
    fTree->Branch("nGDau",   &fnGrandDaughters,"nGDau/I");
    fTree->Branch("nRef1G",  &fmapRef1Gen,"nRef1G/I");
    //fTree->Branch("nRef2G",  &fmapRef2Gen,"nRef2G/I");
    fTree->Branch("isSig",   &fIsSig,     "isSig/I");

    // sigtrack info
    fTree->Branch("Sgptpr",  &fSgptpr,  "Sgptpr/D");
    fTree->Branch("Sgetapr", &fSgetapr, "Sgetapr/D");
    fTree->Branch("Sgphipr", &fSgphipr, "Sgphipr/D");
    fTree->Branch("Sgptpi",  &fSgptpi,  "Sgptpi/D");
    fTree->Branch("Sgetapi", &fSgetapi, "Sgetapi/D");
    fTree->Branch("Sgphipi", &fSgphipi, "Sgphipi/D");

    // Generator tree
    fGenTree = new TTree("genevents", "genevents");
    fGenTree->Branch("nDau",    &fgnDaughters,"nDau/I");
    fGenTree->Branch("nGDau",   &fgnGrandDaughters,"nGDau/I");
    fGenTree->Branch("nRef1G",  &fgmapRef1Gen,"nRef1G/I");
    //fGenTree->Branch("nRef2G",  &fgmapRef2Gen,"nRef2G/I");
    fGenTree->Branch("mlb",     &fgmlb,    "gmlb/D");
    fGenTree->Branch("mlbSwap", &fgmlbsw,  "gmlbsw/D");
    fGenTree->Branch("ml0",     &fgml0,    "gml0/D");
    fGenTree->Branch("ml0Swap", &fgml0sw,  "gml0sw/D");
    fGenTree->Branch("ptpr",    &fgptpr,   "ptpr/D");
    fGenTree->Branch("ptpi",    &fgptpi,   "ptpi/D");
    fGenTree->Branch("ppr",     &fgppr,    "ppr/D");
    fGenTree->Branch("ppi",     &fgppi,    "ppi/D");
    fGenTree->Branch("ptl0",    &fgptl0,   "ptl0/D");
    fGenTree->Branch("pl0",     &fgpl0,    "pl0/D");
    fGenTree->Branch("dRprpi",  &fgdRprpi, "dRprpi/D");
}

// ----------------------------------------------------------------------
void lambdaReader::readCuts(TString filename, int dump)
{
    char  buffer[200];
    fCutFile = filename;
    if (dump) cout << "==> lambdaReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
    sprintf(buffer, "%s", fCutFile.Data());
    ifstream is(buffer);
    char CutName[100];
    float CutValue;
    int ok(0);

    TString fn(fCutFile.Data());

    if (dump)
    {
        cout << "====================================" << endl;
        cout << "==> lambdaReader: Cut file  " << fCutFile.Data() << endl;
        cout << "------------------------------------" << endl;
    }

    TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
    hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
    int ibin;
    while (is.getline(buffer, 200, '\n'))
    {
        ok = 0;
        if (buffer[0] == '#')
        {
            continue;
        }
        if (buffer[0] == '/')
        {
            continue;
        }
        sscanf(buffer, "%s %f", CutName, &CutValue);

        // Frank cuts
        if (!strcmp(CutName, "cutLbCandidate"))
        {
            cutLbCandidate = int(CutValue);
            ok = 1;
            if (dump) cout << "cutLbCandidate:    " << cutLbCandidate << endl;
        }

        if (!strcmp(CutName, "cutMuId1"))
        {
            cutMuId1 = int(CutValue);
            ok = 1;
            if (dump) cout << "cutMuId1:         " << cutMuId1 << endl;
        }

        if (!strcmp(CutName, "cutMuId2"))
        {
            cutMuId2 = int(CutValue);
            ok = 1;
            if (dump) cout << "cutMuId2:         " << cutMuId2 << endl;
        }

        if (!strcmp(CutName, "cutMjpMin"))
        {
            cutMjpMin = CutValue;
            ok = 1;
            if (dump) cout << "cutMjpMin:        " << cutMjpMin << endl;
        }

        if (!strcmp(CutName, "cutMjpMax"))
        {
            cutMjpMax = CutValue;
            ok = 1;
            if (dump) cout << "cutMjpMax:        " << cutMjpMax << endl;
        }

        if (!strcmp(CutName, "cutAlphal0Max"))
        {
            cutAlphal0Max = CutValue;
            ok = 1;
            if (dump) cout << "cutAlphal0Max:    " << cutAlphal0Max << endl;
        }

        if (!strcmp(CutName, "cutMl0Max"))
        {
            cutMl0Max = CutValue;
            ok = 1;
            if (dump) cout << "cutMl0Max:        " << cutMl0Max << endl;
        }

        if (!strcmp(CutName, "cutD3dl0Min"))
        {
            cutD3dl0Min = CutValue;
            ok = 1;
            if (dump) cout << "cutD3dl0Min:      " << cutD3dl0Min << endl;
        }

        if (!strcmp(CutName, "cutPtl0Min"))
        {
            cutPtl0Min = CutValue;
            ok = 1;
            if (dump) cout << "cutPtl0Min:       " << cutPtl0Min << endl;
        }

        // alte Cuts --------------------------------------------
        if (!strcmp(CutName, "TYPE"))
        {
            TYPE = int(CutValue);
            ok = 1;
            if (dump) cout << "TYPE:           " << TYPE << endl;
        }

        if (!strcmp(CutName, "MCTYPE"))
        {
            MCTYPE = int(CutValue);
            ok = 1;
            if (dump) cout << "MCTYPE:           " << MCTYPE << endl;
        }

        if (!strcmp(CutName, "CHARMPTLO"))
        {
            CHARMPTLO = CutValue;
            ok = 1;
            if (dump) cout << "CHARMPTLO:           " << CHARMPTLO << " GeV" << endl;
            ibin = 11;
            hcuts->SetBinContent(ibin, CHARMPTLO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(B_{s}) [GeV]");
        }

        if (!strcmp(CutName, "VTXCHI2"))
        {
            VTXCHI2 = CutValue;
            ok = 1;
            if (dump) cout << "VTXCHI2:           " << VTXCHI2 << " " << endl;
            ibin = 11;
            hcuts->SetBinContent(ibin, VTXCHI2);
            hcuts->GetXaxis()->SetBinLabel(ibin, "#chi^{2}");
        }

        if (!strcmp(CutName, "CHARMETALO"))
        {
            CHARMETALO = CutValue;
            ok = 1;
            if (dump) cout << "CHARMETALO:           " << CHARMETALO << endl;
            ibin = 12;
            hcuts->SetBinContent(ibin, CHARMETALO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(B_{s})");
        }

        if (!strcmp(CutName, "CHARMETAHI"))
        {
            CHARMETAHI = CutValue;
            ok = 1;
            if (dump) cout << "CHARMETAHI:           " << CHARMETAHI << endl;
            ibin = 13;
            hcuts->SetBinContent(ibin, CHARMETAHI);
            hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{max}(B_{s})");
        }

        if (!strcmp(CutName, "KAPTLO"))
        {
            KAPTLO = CutValue;
            ok = 1;
            if (dump) cout << "KAPTLO:           " << KAPTLO << " GeV" << endl;
            ibin = 21;
            hcuts->SetBinContent(ibin, KAPTLO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(K)");
        }

        if (!strcmp(CutName, "PIPTLO"))
        {
            PIPTLO = CutValue;
            ok = 1;
            if (dump) cout << "PIPTLO:           " << PIPTLO << " GeV" << endl;
            ibin = 21;
            hcuts->SetBinContent(ibin, PIPTLO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(#pi)");
        }
        if (!strcmp(CutName, "MUPTLO"))
        {
            MUPTLO = CutValue;
            ok = 1;
            if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
            ibin = 21;
            hcuts->SetBinContent(ibin, MUPTLO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(#mu)");
        }

        if (!strcmp(CutName, "MUPTHI"))
        {
            MUPTHI = CutValue;
            ok = 1;
            if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
            ibin = 22;
            hcuts->SetBinContent(ibin, MUPTHI);
            hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(#mu)");
        }

        if (!strcmp(CutName, "MUETALO"))
        {
            MUETALO = CutValue;
            ok = 1;
            if (dump) cout << "MUETALO:           " << MUETALO << endl;
            ibin = 23;
            hcuts->SetBinContent(ibin, MUETALO);
            hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(#mu)");
        }

        if (!strcmp(CutName, "MUETAHI"))
        {
            MUETAHI = CutValue;
            ok = 1;
            if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
            ibin = 24;
            hcuts->SetBinContent(ibin, MUETAHI);
            hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{max}(#mu)");
        }


        if (!ok) cout << "==> lambdaReader: ERROR: Don't know about variable " << CutName << endl;
    }

    if (dump)  cout << "------------------------------------" << endl;

}
