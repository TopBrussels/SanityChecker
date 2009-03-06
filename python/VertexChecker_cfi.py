import FWCore.ParameterSet.Config as cms

vertex = cms.EDAnalyzer('VertexChecker',
	vertexName   = cms.InputTag('offlinePrimaryVertices'),

) 
