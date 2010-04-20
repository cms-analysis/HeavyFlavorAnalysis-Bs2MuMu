#include "dpReader.hh"
#include <cstdlib>

#define require_true(COND,LABEL) if( !(COND) ) goto LABEL
#define MMUON 0.10566
#define MKAON 0.49368
#define MPION 0.139570
#define MDPLUS 1.86962
#define MB_0 5.27953
#define PI 3.141592654

// test a <= x < b
static inline int in_interval(int x, int a, int b)
{
	return (a <= x) && (x < b);
} // in_interval()

dpReader::dpReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),reduced_tree(NULL)
{
  fMomentumPtr = &fMomentum;
} // dpReader()

dpReader::~dpReader()
{ } // ~dpReader()

void dpReader::eventProcessing()
{
	int j,nc;
	TAnaCand *pCand;
	TAnaTrack *pi1Trk, *pi2Trk, *kaTrk, *muTrk;


	// Fill a reduced tree with the following types
	//	- momentum
	//	- mass
	//	- type
	nc = fpEvt->nCands();
        runNr = fpEvt->fRunNumber;
	goodrun = CheckGoodRun(runNr);

	for (j=0; j<nc; j++) {
	    pCand = fpEvt->getCand(j);
		
            // Save in the tree
	    fCandidate = pCand->fType;
	    fMass = pCand->fMass;
	    if(fCandidate==10030 && fabs(fMass-MDPLUS)<0.3 ){

	          dppt= dpm= dpdxy= dpdxyE= dpd3d= dpd3dE= -9.;
	          bpt= bm= bdxy= bdxyE= bd3d= bd3dE= -9.;
	          fdistPV= fdistDp= fdeltaR= fdeltaRMu= fdeltaRMax =-9.;
		  chi2vtx= dalpha= -9.;
		  nPV=0;

	          fMomentum = pCand->fPlab;
		  
                  dppt   = fMomentum.Perp();
                  dpm    = pCand->fMass;

                  dpdxy  = pCand->fVtx.fDxy;
                  dpdxyE = pCand->fVtx.fDxyE;
                  dpd3d  = pCand->fVtx.fD3d;
                  dpd3dE = pCand->fVtx.fD3dE;

		  if((pCand->fSig2-pCand->fSig1) != 2) {
		    cout <<  "Error detected for 10030: nTracks="<<pCand->fSig2-pCand->fSig1<<endl; 
		    continue;
		  }
		  kaTrk = fpEvt->getSigTrack(pCand->fSig1);
		  pi1Trk = fpEvt->getSigTrack(pCand->fSig2-1);
		  pi2Trk = fpEvt->getSigTrack(pCand->fSig2);


		  fdeltaR   = fMomentum.DeltaR(kaTrk->fPlab);
		  double dR =  fMomentum.DeltaR(pi1Trk->fPlab);
		  fdeltaR = dR<fdeltaR ? dR : fdeltaR;
		  dR =  fMomentum.DeltaR(pi2Trk->fPlab);
		  fdeltaR = dR<fdeltaR ? dR : fdeltaR;

		  fdeltaRMax   = fMomentum.DeltaR(kaTrk->fPlab);
		  dR =  fMomentum.DeltaR(pi1Trk->fPlab);
		  fdeltaRMax = dR>fdeltaRMax ? dR : fdeltaRMax;
		  dR =  fMomentum.DeltaR(pi2Trk->fPlab);
		  fdeltaRMax = dR>fdeltaRMax ? dR : fdeltaRMax;

		  ptmin=kaTrk->fPlab.Perp();
		  etamax=fabs(kaTrk->fPlab.Eta());

		  double pt_tmp=pi1Trk->fPlab.Perp();
		  double eta_tmp=fabs(pi1Trk->fPlab.Eta());
		  if(pt_tmp<ptmin) ptmin=pt_tmp;
		  if(eta_tmp>etamax) etamax=eta_tmp;
		  pt_tmp=pi2Trk->fPlab.Perp();
		  eta_tmp=fabs(pi2Trk->fPlab.Eta());
		  if(pt_tmp<ptmin) ptmin=pt_tmp;
		  if(eta_tmp>etamax) etamax=eta_tmp;

		  dcamin = pCand->fVar1;
		  dcamax = pCand->fVar2;

		  chi2vtx = pCand->fVtx.fChi2;
		  nPV=fpEvt->nPV();
		  if(nPV>0){
		    TVector3 pVtx=fpEvt->getPV(0)->fPoint;
		    TVector3 svtx = pCand->fVtx.fPoint-pVtx;
		    dalpha = fMomentum.Angle(svtx);
		    //if (dalpha>PI/2) dalpha-=PI;
		  }

        	  fTruth = checkTruth(pCand, TYPE);
  		  if(dcamax<DCAMAX && dpdxy/dpdxyE>SIGVTX) reduced_tree->Fill();
	     } else if (fCandidate==20030 && fabs(fMass-MB_0)<0.5){
	          dppt= dpm= dpdxy= dpdxyE= dpd3d= dpd3dE= -9.;
	          bpt= bm= bdxy= bdxyE= bd3d= bd3dE= -9.;
	          fdistPV= fdistDp= fdeltaR= fdeltaRMu= fdeltaRMax =-9.;
		  chi2vtx= dalpha= -9.;
		  nPV=0;
	       
	          fMomentum = pCand->fPlab;

                  bpt   = fMomentum.Perp();
                  bm    = pCand->fMass;
                  bdxy  = pCand->fVtx.fDxy;
                  bdxyE = pCand->fVtx.fDxyE;
                  bd3d  = pCand->fVtx.fD3d;
                  bd3dE = pCand->fVtx.fD3dE;

		  if((pCand->fSig2-pCand->fSig1) != 3) {
		    cout <<  "Error detected for 20030: nTracks="<<pCand->fSig2-pCand->fSig1<<endl;
		    continue;
		  }
		  kaTrk = fpEvt->getSigTrack(pCand->fSig1);
		  pi1Trk = fpEvt->getSigTrack(pCand->fSig1+1);
		  pi2Trk = fpEvt->getSigTrack(pCand->fSig1+2);
		  muTrk = fpEvt->getSigTrack(pCand->fSig2);

		  fdeltaR   = fMomentum.DeltaR(kaTrk->fPlab);
		  double dR =  fMomentum.DeltaR(pi1Trk->fPlab);
		  fdeltaR = dR<fdeltaR ? dR : fdeltaR;
		  dR =  fMomentum.DeltaR(pi2Trk->fPlab);
		  fdeltaR = dR<fdeltaR ? dR : fdeltaR;

		  fdeltaRMax   = fMomentum.DeltaR(kaTrk->fPlab);
		  dR =  fMomentum.DeltaR(pi1Trk->fPlab);
		  fdeltaRMax = dR>fdeltaRMax ? dR : fdeltaRMax;
		  dR =  fMomentum.DeltaR(pi2Trk->fPlab);
		  fdeltaRMax = dR>fdeltaRMax ? dR : fdeltaRMax;

		  fdeltaRMu= muTrk->fPlab.DeltaR(kaTrk->fPlab);
		  dR =  muTrk->fPlab.DeltaR(pi1Trk->fPlab);
		  fdeltaRMu = dR<fdeltaRMu ? dR : fdeltaRMu;
		  dR =  muTrk->fPlab.DeltaR(pi2Trk->fPlab);
		  fdeltaRMu = dR<fdeltaRMu ? dR : fdeltaRMu;

		  ptmin=kaTrk->fPlab.Perp();
		  etamax=fabs(kaTrk->fPlab.Eta());

		  double pt_tmp=pi1Trk->fPlab.Perp();
		  double eta_tmp=fabs(pi1Trk->fPlab.Eta());
		  if(pt_tmp<ptmin) ptmin=pt_tmp;
		  if(eta_tmp>etamax) etamax=eta_tmp;
		  pt_tmp=pi2Trk->fPlab.Perp();
		  eta_tmp=fabs(pi2Trk->fPlab.Eta());
		  if(pt_tmp<ptmin) ptmin=pt_tmp;
		  if(eta_tmp>etamax) etamax=eta_tmp;

		  ptmu=muTrk->fPlab.Perp();
		  etamu=fabs(muTrk->fPlab.Eta());

 		  dcamin=pCand->fVar1;
		  dcamax=pCand->fVar2;

		  chi2vtx = pCand->fVtx.fChi2;
		  nPV=fpEvt->nPV();
		  if(nPV>0){
		    TVector3 pVtx=fpEvt->getPV(0)->fPoint;
		    TVector3 svtx = pCand->fVtx.fPoint-pVtx;
		    dalpha = fMomentum.Angle(svtx);
		    //if (dalpha>PI/2) dalpha-=PI;
		  }

	          fTruth = checkTruth(pCand, TYPE);
		  if(dcamax<DCAMAX && bdxy/bdxyE>SIGVTX) reduced_tree->Fill();
	    }
	}
} // dplusReader::eventProcessing()

void dpReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T","Candidate Mass / Eta / Charge");
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("p","TVector3",&fMomentumPtr);
	reduced_tree->Branch("mass",&fMass,"mass/D");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
	reduced_tree->Branch("goodrun",&goodrun,"goodrun/I");
	reduced_tree->Branch("nPV",     &nPV ,    "nPV/I");
        reduced_tree->Branch("run",     &runNr ,    "runNr/I");

	reduced_tree->Branch("chi2vtx",&chi2vtx,"chi2vtx/D");
	reduced_tree->Branch("dalpha",&dalpha,"dalpha/D");

        reduced_tree->Branch("dppt",   &dppt,   "dppt/D");
        reduced_tree->Branch("dpm",    &dpm,    "dpm/D");
        reduced_tree->Branch("dpdxy",  &dpdxy,  "dpdxy/D");
        reduced_tree->Branch("dpdxyE", &dpdxyE, "dpdxyE/D");
        reduced_tree->Branch("dpd3d",  &dpd3d,  "dpd3d/D");
        reduced_tree->Branch("dpd3dE", &dpd3dE, "dpd3dE/D");

        reduced_tree->Branch("bpt",   &bpt,   "bpt/D");
        reduced_tree->Branch("bm",    &bm,    "bm/D");
        reduced_tree->Branch("bdxy",  &bdxy,  "bdxy/D");
        reduced_tree->Branch("bdxyE", &bdxyE, "bdxyE/D");
        reduced_tree->Branch("bd3d",  &bd3d,  "bd3d/D");
        reduced_tree->Branch("bd3dE", &bd3dE, "bd3dE/D");

	reduced_tree->Branch("dcamin", &dcamin, "dcamin/D");
	reduced_tree->Branch("dcamax", &dcamax, "dcamax/D");
	reduced_tree->Branch("ptmin", &ptmin, "ptmin/D"); 
	reduced_tree->Branch("etamax", &etamax, "etamax/D");
	reduced_tree->Branch("ptmu", &ptmu, "ptmu/D"); 
	reduced_tree->Branch("etamu", &etamu, "etamu/D");  

        reduced_tree->Branch("deltaRMin",     &fdeltaR ,    "deltaRMin/D");
        reduced_tree->Branch("deltaRMax",     &fdeltaRMax ,    "deltaRMax/D");
        reduced_tree->Branch("deltaRMu",     &fdeltaRMu ,    "deltaRMu/D");

} // dpReader::bookHist()



int dpReader::checkTruth(TAnaCand *cand, int truth_type)
{
	TAnaTrack *sgTrack;
	TAnaTrack *recTrack;
	TGenCand *truthParticle;
	TGenCand *trackParticle;
	int j,succes = 0;
	int nSigs,nRecs,nGens;
	
	
	nSigs = fpEvt->nSigTracks();
	nRecs = fpEvt->nRecTracks();
	nGens = fpEvt->nGenCands();
	
	// get the first track and get the originating particle of type 'truth_type' => truthParticle
	require_true(in_interval(cand->fSig1, 0, nSigs),bail);
	sgTrack = fpEvt->getSigTrack(cand->fSig1);
	
	require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
	recTrack = fpEvt->getRecTrack(sgTrack->fIndex);
	
	require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
	truthParticle = fpEvt->getGenCand(recTrack->fGenIndex);
	
	while (abs(truthParticle->fID) != truth_type && in_interval(truthParticle->fMom1, 0, nGens)){
	  truthParticle = fpEvt->getGenCand(truthParticle->fMom1);
	  //	  if(abs(truthParticle->fID)==411) cout << "Truth type = "<<truthParticle->fID<<"   "<<truth_type<<endl;
	}

	// check if our original is a valid.
	require_true(abs(truthParticle->fID) == truth_type,bail);
	
	// reconstruct the other tracks
	for (j = cand->fSig1+1; j <=cand->fSig2; j++) {
		
		require_true(in_interval(j, 0, nSigs),bail);
		sgTrack = fpEvt->getSigTrack(j);
		
		require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
		recTrack = fpEvt->getRecTrack(sgTrack->fIndex);
		
		require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
		trackParticle = fpEvt->getGenCand(recTrack->fGenIndex);
		
		while (abs(trackParticle->fID)!=truth_type && in_interval(trackParticle->fMom1, 0, nGens))
			trackParticle = fpEvt->getGenCand(trackParticle->fMom1);
		
		// check the particle
		require_true(trackParticle->fNumber==truthParticle->fNumber,bail);
	}
	
	// still here? then every track originated from the right particle
	succes = 1;
	
bail:
	return succes;
} // dpReader::checkTruth()


int dpReader::CheckGoodRun(int run){
  for (vector<int>::iterator it=goodRuns.begin(); it!=goodRuns.end(); it++){
    if(run==*it) return 1;
  }
  return 0;
}


void dpReader::startAnalysis(){
  ifstream infile("goodruns.txt",ifstream::in);
  int tmp_run;
  while(!infile.eof()){
    infile >> tmp_run;
    goodRuns.push_back(tmp_run);
    cout <<"Good tracker run: "<<tmp_run<<endl;
  }
  infile.close();
}


void dpReader::readCuts(TString filename, int dump) {
  dump=1;
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> dpReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> dpReader: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue);
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!strcmp(CutName, "SIGVTX")) {
      SIGVTX = CutValue;
      if (dump) cout << "SIGVTX:           " << SIGVTX << endl;
    }

    if (!strcmp(CutName, "DCAMAX")) {
      DCAMAX = CutValue;
      if (dump) cout << "DCAMAX:           " << DCAMAX << endl;
    }
  }
}
