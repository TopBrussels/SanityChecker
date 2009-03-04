import FWCore.ParameterSet.Config as cms
  
process = cms.Process("SanityCheck")

#Message Logger (more info see: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMessageLogger)
process.MessageLogger = cms.Service("MessageLogger",
   
    destinations = cms.untracked.vstring("Warning_SC", "ErrorSummary_SC", "InfoSummary_SC"), # 3 files for 3 different type of output

    debugModules = cms.untracked.vstring('*'), #for all modules

    categories = cms.untracked.vstring("inputChain","decayChain","NoDataFound","LinkBroken","SummaryError","MainResults"), # list of categories
                                                                                                                           # inputChain and decayChain come from the TopDecaySubsetModule and need to be suppressedin warning- and info-summary

    Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string("DEBUG"), 
                                     NoDataFound = cms.untracked.PSet(reportEvery = cms.untracked.int32(10)), 
                                     LinkBroken = cms.untracked.PSet(reportEvery = cms.untracked.int32(2)),
                                     inputChain = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                     decayChain = cms.untracked.PSet(limit = cms.untracked.int32(0))
				     ),
    ErrorSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("ERROR"), categories = cms.untracked.vstring("SummaryError") ),
    InfoSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("INFO"), categories = cms.untracked.vstring("MainResults"),
                                         inputChain = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                         decayChain = cms.untracked.PSet(limit = cms.untracked.int32(0)))
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
    )
)

## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.load("TopBrussels.SanityChecker.TtGenEventChecker_cfi")
process.load("TopBrussels.SanityChecker.ResolutionChecker_cfi")
process.load("TopBrussels.SanityChecker.KinematicsChecker_cfi")


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('SanityChecker.root')
	)
	

process.p = cms.Path(process.makeGenEvt * (process.TtGenEventChecker + process.Resolutions_lJets + process.Resolutions_bJets + process.kinematics + process.jetmet))
