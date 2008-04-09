#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTree.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpStuff.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpGenerator.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpTracks.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpMuons.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpCandidates.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpSignal.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/L1TrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HLTrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDumpTrigger.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HFTree);
DEFINE_ANOTHER_FWK_MODULE(HFDumpStuff);
DEFINE_ANOTHER_FWK_MODULE(HFDumpGenerator);
DEFINE_ANOTHER_FWK_MODULE(HFDumpTracks);
DEFINE_ANOTHER_FWK_MODULE(HFDumpMuons);
DEFINE_ANOTHER_FWK_MODULE(HFDumpSignal);
DEFINE_ANOTHER_FWK_MODULE(HFDumpCandidates);
DEFINE_ANOTHER_FWK_MODULE(L1TrigReport);
DEFINE_ANOTHER_FWK_MODULE(HLTrigReport);
DEFINE_ANOTHER_FWK_MODULE(HFDumpTrigger);

