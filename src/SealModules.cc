#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSTree.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpStuff.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpGenerator.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpTracks.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpMuons.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpCandidates.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpSignal.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/L1TrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HLTrigReport.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/BSDumpTrigger.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(BSTree);
DEFINE_ANOTHER_FWK_MODULE(BSDumpStuff);
DEFINE_ANOTHER_FWK_MODULE(BSDumpGenerator);
DEFINE_ANOTHER_FWK_MODULE(BSDumpTracks);
DEFINE_ANOTHER_FWK_MODULE(BSDumpMuons);
DEFINE_ANOTHER_FWK_MODULE(BSDumpSignal);
DEFINE_ANOTHER_FWK_MODULE(BSDumpCandidates);
DEFINE_ANOTHER_FWK_MODULE(L1TrigReport);
DEFINE_ANOTHER_FWK_MODULE(HLTrigReport);
DEFINE_ANOTHER_FWK_MODULE(BSDumpTrigger);

