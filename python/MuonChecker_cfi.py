import FWCore.ParameterSet.Config as cms

muonchecker = cms.EDFilter("MuonChecker",
	muonTag   = cms.InputTag("selectedLayer1Muons")

)
