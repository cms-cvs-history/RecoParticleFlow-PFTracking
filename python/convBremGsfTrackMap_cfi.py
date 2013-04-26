import FWCore.ParameterSet.Config as cms

convBremGsfTrackMap = cms.EDProducer("ConvBremGsfTrackMapProducer",
                                     TrackQuality = cms.string('highPurity'),
                                     UseQuality = cms.bool(True),
                                     GsfTrackModuleLabel = cms.InputTag("electronGsfTracks"),
                                     TrackCollection = cms.InputTag("generalTracks"),
                                     PFEcalClusters  = cms.InputTag("particleFlowClusterECAL"),                 
                                     PrimaryVertexLabel = cms.InputTag("offlinePrimaryVertices"),
                                     useConvBremFinder  = cms.bool(True),
                                     pf_convBremFinderID_mvaCut =  cms.double(0.2),
                                     pf_convBremFinderID_mvaWeightFile = cms.string('RecoParticleFlow/PFTracking/data/MVAnalysis_BDT.weights_convBremFinder_19Apr_IntToFloat.txt'),   
                                     )


