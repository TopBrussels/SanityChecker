import FWCore.ParameterSet.Config as cms

kinematics = cms.EDAnalyzer('KinematicsChecker',
	jetsName        = cms.InputTag('selectedLayer1Jets'),
	muonsName       = cms.InputTag('selectedLayer1Muons'),
        metsName        = cms.InputTag('selectedLayer1METs'),
        matchingAlgo    = cms.int32(0),
        useMaxDist      = cms.bool(True),
        useDeltaR       = cms.bool(False),
        maxDist         = cms.double(0.5),
        jetsAcceptance  = cms.vdouble(2.4,15),                 #acceptance region for jets in jet collection   (eta,pt)
        muonsAcceptance  = cms.vdouble(2.5,15)                 #acceptance region for muons in muon collection (eta,pt)
) 
