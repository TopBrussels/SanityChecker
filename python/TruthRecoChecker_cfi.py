import FWCore.ParameterSet.Config as cms

TruthReco = cms.EDAnalyzer('TruthReco',
	genEventCollectionName	= cms.InputTag('genEvt'),
	patMuons								= cms.InputTag("selectedLayer1Muons"),
	patJets									= cms.InputTag("selectedLayer1Jets"),
	patMET									= cms.InputTag("selectedLayer1METs"),
	useMatchingFromPAT			= cms.bool(False),	
  matchingAlgo    				= cms.int32(2),
  useMaxDist      				= cms.bool(True),
  useDeltaR       				= cms.bool(True),
  maxDist         				= cms.double(0.7),
	doMuonIsolation					=	cms.bool(True),		
	muonRelIso							= cms.double(0.),
	muonECALVetoConeEt			= cms.double(4.),
	muonHCALVetoConeEt			= cms.double(6.),
	muonPt									= cms.double(15),
	jetPt										= cms.double(15),
	muonEta									= cms.double(2.5),
	jetEta									= cms.double(2.5),
	muonMinDR								= cms.double(0.5),
	doNonOverlappingJets		= cms.bool(False),	#not yet implemented
	nonOverlappingJetsMinDR	= cms.double(1.0),
	NJets										= cms.int32(4),	
	calculateJECinbins			= cms.bool(True),	
	etabinValues						= cms.vdouble(0, 0.17, 0.35, 0.5, 0.7, 0.9, 1.15, 1.4, 1.7, 2.1, 2.5),
	pTbinValues							= cms.vdouble(20,30,40,50,60,70,80,90,110,130,200),
	minDRunmatchedParton		=	cms.double(0.5)
)