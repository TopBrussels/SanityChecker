import FWCore.ParameterSet.Config as cms

jetmet = cms.EDAnalyzer('JetMetChecker',
	jetsName     = cms.InputTag('selectedLayer1Jets'),
	vertexName   = cms.InputTag('offlinePrimaryVertices'),
#	muonsName    = cms.InputTag('selectedLayer1Muons'),
        metsName     = cms.InputTag('selectedLayer1METs'),
        matchingAlgo    = cms.int32(2),
        useMaxDist      = cms.bool(True),
        useDeltaR       = cms.bool(True),
        maxDist         = cms.double(0.3)
) 
