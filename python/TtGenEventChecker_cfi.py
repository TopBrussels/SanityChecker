import FWCore.ParameterSet.Config as cms


TtGenEventChecker = cms.EDAnalyzer('TtGenEventChecker',
      genEvtCollectionName=cms.InputTag('genEvt'), 
      doit_all = cms.untracked.bool(True),
      doit_semileptonic = cms.untracked.bool(False),
      doit_semileptonicMuon = cms.untracked.bool(False)
)

	
