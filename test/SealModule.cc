#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoParticleFlow/PFTracking/test/TestPreId.h"
#include "RecoParticleFlow/PFTracking/test/FakeEval.h"
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TestPreId);
DEFINE_ANOTHER_FWK_MODULE(FakeEval);

