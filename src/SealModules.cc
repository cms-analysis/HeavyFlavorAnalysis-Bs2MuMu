#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTree.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpStuff.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpGenerator.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpTracks.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpMuons.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpCandidates.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpSignal.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/L1TrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HLTrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BmmDumpTrigger.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HFTree);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpStuff);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpGenerator);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpTracks);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpMuons);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpSignal);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpCandidates);
DEFINE_ANOTHER_FWK_MODULE(L1TrigReport);
DEFINE_ANOTHER_FWK_MODULE(HLTrigReport);
DEFINE_ANOTHER_FWK_MODULE(BmmDumpTrigger);

