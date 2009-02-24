import FWCore.ParameterSet.Config as cms
  
process = cms.Process("Demo")

#Message Logger
#process.MessageLogger = cms.Service("MessageLogger",
#    # 3 files for 3 different type of output
#    destinations = cms.untracked.vstring("Warning_SC", "ErrorSummary_SC", "InfoSummary_SC"),
#    debugModules = cms.untracked.vstring('Resolutions_bJets','Resolutions_lJets'),#for all modules
#    categories = cms.untracked.vstring("NoDataFound","LinkBroken","SummaryError","MainResults"), # list of categories
#    #Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'), categories = cms.untracked.vstring("cat") ),
#    Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string("DEBUG"), 
#                                     NoDataFound = cms.untracked.PSet(reportEvery = cms.untracked.int32(10)), 
#                                     LinkBroken = cms.untracked.PSet(reportEvery = cms.untracked.int32(2)) 
#				     ),
#    ErrorSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("ERROR"), categories = cms.untracked.vstring("SummaryError") ),
#    InfoSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("INFO"), categories = cms.untracked.vstring("MainResults") )
#)

process.load("TopBrussels.SanityChecker.PATLayer1_Ttjets_MG_NoSel_input_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
process.decaySubset.pdgId = cms.uint32(999)

process.Resolutions_electrons = cms.EDAnalyzer("ResolutionChecker",
	object    	= cms.string('electron'), 
	label  		= cms.string('selectedLayer1Electrons'),
	minMatchingDR	= cms.double(0.1),
	etabinValues 	= cms.vdouble(0,0.17,0.35,0.5,0.7,0.9,1.15,1.3,1.6,1.9,2.5),
	pTbinValues  	= cms.vdouble(15,22,28,35,41,49,57,68,81,104,200)
)
   
process.Resolutions_muons  = cms.EDAnalyzer("ResolutionChecker",
	object    	= cms.string('muon'),
	label  		= cms.string('selectedLayer1Muons'),
	minMatchingDR	= cms.double(0.1),
	etabinValues 	= cms.vdouble(0,0.17,0.35,0.5,0.7,0.9,1.15,1.3,1.6,1.9,2.5),
	pTbinValues  	= cms.vdouble(15,22,28,35,41,49,57,68,81,104,200)
)
     
process.Resolutions_lJets  = cms.EDAnalyzer("ResolutionChecker",
	object    	= cms.string('lJets'), 
	label  		= cms.string('selectedLayer1Jets'),
	minMatchingDR	= cms.double(0.3),
	etabinValues 	= cms.vdouble(0, 0.17, 0.35, 0.5, 0.7, 0.9, 1.15, 1.4, 1.7, 2.1, 2.5),
	pTbinValues  	= cms.vdouble(20,30,40,50,60,70,80,90,110,130,200)
)
 
process.Resolutions_bJets  = cms.EDAnalyzer("ResolutionChecker",
	object    	= cms.string('bJets'),
	label  		= cms.string('selectedLayer1Jets'),
	minMatchingDR	= cms.double(0.3),
	etabinValues 	= cms.vdouble(0, 0.17, 0.35, 0.5, 0.7, 0.9, 1.15, 1.4, 1.7, 2.1, 2.5),
	pTbinValues  	= cms.vdouble(20,30,40,50,60,70,80,90,110,130,200)
)

process.Resolutions_met  = cms.EDAnalyzer("ResolutionChecker",
	object    	= cms.string('met'),
	label  		= cms.string('selectedLayer1METs'),
	minMatchingDR	= cms.double(1000.),
	pTbinValues  	= cms.vdouble(20,29,37,44,51,59,69,80,96,122,200)	
) 

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('DummyChecker_testresolutions.root')
	)
	

process.p = cms.Path(process.makeGenEvt *(process.Resolutions_lJets + process.Resolutions_bJets))
