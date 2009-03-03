import FWCore.ParameterSet.Config as cms

jetmet = cms.EDAnalyzer('JetMetChecker',
	jetsName     = cms.InputTag('selectedLayer1Jets'),
	vertexName   = cms.InputTag('offlinePrimaryVertices'),
#	muonsName    = cms.InputTag('selectedLayer1Muons'),
        metsName     = cms.InputTag('selectedLayer1METs')

) 
