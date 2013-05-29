// ----------------------------------------------------------------------
void runAll() {
  runBDT(0, 0); 
  runBDT(1, 1); 
}

// ----------------------------------------------------------------------
void runBDT(int seed, string filename, string vars, string pars, int iChannel = 0, int year = 2011) {
  gSystem->Load("libTMVA.so");
  gRandom->SetSeed(seed);

  int npars(5);
  float loVal[] = {200, 100, 2,  20, 0.05};
  float hiVal[] = {2e3, 1e3, 6, 100, 1.20};

  float cut(0); 
  int ntrees, nevts, maxdepth, ncuts, nnodesmax; 
  float beta;
  for (int i = 0; i < npars; ++i) {
    cut = gRandom->Rndm()*(hiVal[i]-loVal[i]) + loVal[i];
    if (0 == i) ntrees = static_cast<int>(static_cast<int>(cut)/10.)*10; 
    if (1 == i) nevts  = static_cast<int>(static_cast<int>(cut)/10.)*10;
    if (2 == i) {
      maxdepth = static_cast<int>(static_cast<int>(cut)/1.)*1; 
      nnodesmax = 100*maxdepth;
    }
    if (3 == i) ncuts = static_cast<int>(static_cast<int>(cut)/1.)*1;
    if (4 == i) beta = static_cast<float>(static_cast<int>(20*cut))/20.;
  }

  string optstring(pars); 
  if (pars == "random") {
    optstring = "!H:V";
    optstring += Form(":NTrees=%d", ntrees);
    optstring += Form(":nEventsMin=%d", nevts);
    optstring += Form(":nCuts=%d:PruneMethod=NoPruning", ncuts);
    optstring += Form(":BoostType=AdaBoost:AdaBoostBeta=%f:SeparationType=GiniIndex", beta);
    optstring += Form(":MaxDepth=%d", maxdepth);
    optstring += Form(":NNodesMax=%d", nnodesmax);
  }

  tmva1 aT(year, vars);
  aT.setBDTParameters(optstring);
  aT.setChannel(iChannel);


  aT.makeAll(seed, filename, 0);
}
