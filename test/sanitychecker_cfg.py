import FWCore.ParameterSet.Config as cms
  
process = cms.Process("SanityCheck")

#Message Logger

process.MessageLogger = cms.Service("MessageLogger",
    # 3 files for 3 different type of output
    destinations = cms.untracked.vstring("Warning_SC", "ErrorSummary_SC", "InfoSummary_SC"),
    debugModules = cms.untracked.vstring('*'),#for all modules
    categories = cms.untracked.vstring("NoDataFound","LinkBroken","SummaryError","MainResults"), # list of categories
    #Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'), categories = cms.untracked.vstring("cat") ),
    Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string("DEBUG"), 
                                     NoDataFound = cms.untracked.PSet(reportEvery = cms.untracked.int32(10)), 
                                     LinkBroken = cms.untracked.PSet(reportEvery = cms.untracked.int32(2)) 
				     ),
    ErrorSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("ERROR"), categories = cms.untracked.vstring("SummaryError") ),
    InfoSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("INFO"), categories = cms.untracked.vstring("MainResults") )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    )
)

process.load("TopBrussels.SanityChecker.TtGenEventChecker_cfi")


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('SanityChecker.root')
	)
	

process.p = cms.Path(process.TtGenEventChecker)
