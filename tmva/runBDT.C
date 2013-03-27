// ----------------------------------------------------------------------
void runAll() {
  runBDT(0, 0); 
  runBDT(1, 1); 
}

// ----------------------------------------------------------------------
void runBDT(int seed, string filename, int iChannel = 0, int year = 2011) {
  gSystem->Load("libTMVA.so");
  gRandom->SetSeed(seed);

  int nvars(6);
  float loVal[] = {100, 100, 2,  20, 0.05, 1e4};
  float hiVal[] = {1e3, 1e3, 7, 100, 0.9, 1e6};

  float cut(0); 
  int ntrees, nevts, maxdepth, ncuts, nnodesmax; 
  float beta;
  for (int i = 0; i < nvars; ++i) {
    cut = gRandom->Rndm()*(hiVal[i]-loVal[i]) + loVal[i];
    if (0 == i) ntrees = static_cast<int>(static_cast<int>(cut)/10.)*10; 
    if (1 == i) nevts  = static_cast<int>(static_cast<int>(cut)/10.)*10;
    if (2 == i) maxdepth = static_cast<int>(static_cast<int>(cut)/1.)*1; 
    if (3 == i) ncuts = static_cast<int>(static_cast<int>(cut)/1.)*1;
    if (4 == i) beta = static_cast<float>(static_cast<int>(20*cut))/20.;
    if (5 == i) nnodesmax = static_cast<int>(cut/10000)*10000.;
  }

  // -- 205
  ntrees = 100; 
  ncuts = 50; 
  nevts = 200; 
  beta = 1.0; 
  maxdepth = 3; 
  nnodesmax = 4e5; 

//   // -- 108
//   ntrees = 100; 
//   ncuts = 50; 
//   nevts = 200; 
//   beta = 0.3; 
//   maxdepth = 4; 
//   nnodesmax = 4e5; 


  ntrees = 800;
  nevts = 50; 
  ncuts = 20; 
  maxdepth = 2;
  nnodesmax = 5;
  beta = 1.0; 


  cout << "seed: " << seed << " ntrees = " << ntrees << " nevts = " << nevts << " maxdepth = " 
       << maxdepth << " ncuts = " << ncuts << " beta = " << beta << " nnodesmax = " << nnodesmax
       << endl;

  tmva1 aT(year);

  aT.setChannel(iChannel);

  aT.setNTrees(ntrees);
  aT.setnEventsMin(nevts); 
  aT.setnCuts(ncuts); 
  aT.setAdaBoostBeta(beta); 

  aT.setMaxDepth(maxdepth);
  aT.setNNodesMax(nnodesmax);

  aT.makeAll(seed, filename, 0);
}
