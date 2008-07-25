import FWCore.ParameterSet.Config as cms

process = cms.Process("BLOCK")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("RecoTracker.IterativeTracking.iterativeTk_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/2008/2/11/RelVal-RelValTTbar-1202721222/0000/06D6FF94-F0D8-DC11-8AC0-000423D9939C.root')
)

process.MessageLogger = cms.Service("MessageLogger",
    pippo = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('pippo')
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.GEN-SIM-DIGI-RECO = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep recoPFRecHits_*_*_*', 
        'keep recoPFClusters_*_*_*', 
        'keep recoPFRecTracks_*_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFSimParticles_*_*_*', 
        'keep recoTracks_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep CaloTowersSorted_*_*_*', 
        'keep edmHepMCProduct_*_*_*'),
    fileName = cms.untracked.string('blocks.root')
)

process.newiterativeTk = cms.Sequence(process.newTracking*process.iterTracking)
process.p = cms.Path(process.newiterativeTk*process.caloTowersPFRec*process.particleFlowCluster*process.particleFlowTrack*process.particleFlowSimParticle*process.particleFlowBlock)
process.outpath = cms.EndPath(process.GEN-SIM-DIGI-RECO)
process.particleFlowCluster.verbose = True
process.particleFlowBlock.verbose = True
process.particleFlowSimParticle.verbose = True
process.elecpreid.TkColList = ['highPurityTracks:', 'secStep:', 'thStep:']


