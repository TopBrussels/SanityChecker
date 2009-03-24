import FWCore.ParameterSet.Config as cms
  
process = cms.Process("SanityCheck")

#Message Logger (more info see: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMessageLogger)
process.MessageLogger = cms.Service("MessageLogger",
   
    destinations = cms.untracked.vstring("Warning_SC", "ErrorSummary_SC", "InfoSummary_SC"), # 3 files for 3 different type of output
    statistics = cms.untracked.vstring("Statistics_SC"),
    debugModules = cms.untracked.vstring('*'), #for all modules

    categories = cms.untracked.vstring("inputChain","decayChain","NoDataFound",
                                       "NoDataFound_NotSemiMu","NoDataFound_LessThen4SelJets","NoDataFound_LessThen1SelMuon",
                                       "NoDataFound_LessThen1IsoMuon","NoDataFound_MuonUnmatched","NoDataFound_NoRadiation",
                                       "NoDataFound_NoRadiation2","NoDataFound_usePATMatching","LinkBroken_noJetsFound",
                                       "LinkBroken_noMuonsFound","LinkBroken_noMetsFound","LinkBroken_noGenEvtFound",
                                       "NoDataFound_MC_noMuon","NoDataFound_MC_noGenEvt","LinkBroken_MC_noGenMuon",
                                       "LinkBroken_ISRJets","LinkBroken_TopRadJets","SummaryError","MainResults",
				       #TtGenEventChecker
				       "NoDataFound_top","NoDataFound_b","NoDataFound_w","NoDataFound_Not2B","NoDataFound_Not1LepFromW",
				       "NoDataFound_Lepton","NoDataFound_Neutrino","NoDataFound_HadTop","NoDataFound_HadW",
				       "NoDataFound_HadB","NoDataFound_HadQ","NoDataFound_HadQBar","NoDataFound_LepB",
				       "NoAcess_ISR","NoAccess_LepRad","NoAccess_HadRad"
				       ), # list of categories
                                    # inputChain and decayChain come from the TopDecaySubsetModule and need to be suppressedin warning- and info-summary
                                    
    Warning_SC = cms.untracked.PSet( threshold = cms.untracked.string("DEBUG"), 
                                     NoDataFound = cms.untracked.PSet(reportEvery = cms.untracked.int32(10)), 
                                     #LinkBroken = cms.untracked.PSet(reportEvery = cms.untracked.int32(2)),
                                     LinkBroken = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     LinkBroken_noGenEvt = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     LinkBroken_ISRJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     LinkBroken_TopRadJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_NotSemiMu = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_LessThen4SelJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_LessThen1SelMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_LessThen1IsoMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_MuonUnmatched = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_NoRadiation = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_NoRadiation2 = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_usePATMatching = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     LinkBroken_noMuonsFound = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     NoDataFound_MC_noGenEvt = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                     LinkBroken_MC_noGenMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     #TtGenEventChecker
				     NoDataFound_top = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_b = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_w = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_Not2B = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_Not1LepFromW = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_Lepton = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_Neutrino = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_HadTop = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_HadW = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_HadB = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_HadQ = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_HadQBar = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoDataFound_LepB = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoAcess_ISR = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoAccess_LepRad = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     NoAccess_HadRad = cms.untracked.PSet(limit = cms.untracked.int32(5)),
				     inputChain = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                     decayChain = cms.untracked.PSet(limit = cms.untracked.int32(0))
				     ),
    ErrorSummary_SC = cms.untracked.PSet( threshold = cms.untracked.string("ERROR"), categories = cms.untracked.vstring("SummaryError") ),
    InfoSummary_SC = cms.untracked.PSet( 	threshold = cms.untracked.string("INFO"), categories = cms.untracked.vstring("MainResults"),
                                                NoDataFound = cms.untracked.PSet(reportEvery = cms.untracked.int32(10)),
                                                LinkBroken = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                         	LinkBroken_noGenEvt = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                         	LinkBroken_ISRJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                         	LinkBroken_TopRadJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_NotSemiMu = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_LessThen4SelJets = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_LessThen1SelMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_LessThen1IsoMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_MuonUnmatched = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_NoRadiation = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_NoRadiation2 = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_usePATMatching = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                LinkBroken_noMuonsFound = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                NoDataFound_MC_noGenEvt = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                                LinkBroken_MC_noGenMuon = cms.untracked.PSet(limit = cms.untracked.int32(5)),
                                         	inputChain = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                         	decayChain = cms.untracked.PSet(limit = cms.untracked.int32(0))),
    #Statistics_SC = cms.untracked.PSet(threshold = cms.untracked.string("DEBUG"))
)

#process.load("TopBrussels.SanityChecker.PATLayer1_Ttjets_MG_input_cfi")
#process.load("TopBrussels.SanityChecker.PATLayer1_Ttjets_MG_NoSel_input_cfi")
#process.load("TopBrussels.SanityChecker.PATLayer1_R1_QCD100to250_MG_input_cfi")
#process.load("TopBrussels.SanityChecker.PATLayer1_R1_TTJets_MG_input_cfi")
process.load("TopBrussels.SanityChecker.PATLayer1_R1_TauolaTTbar_input_cfi")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.load("TopBrussels.SanityChecker.TtGenEventChecker_cfi")
process.load("TopBrussels.SanityChecker.ResolutionChecker_cfi")
process.load("TopBrussels.SanityChecker.KinematicsChecker_cfi")
process.load("TopBrussels.SanityChecker.JetMetChecker_cfi")
process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
process.ttDecaySelection.channel_1 = [0, 1, 0]
process.load("TopBrussels.SanityChecker.TruthRecoChecker_cfi")
process.load("TopBrussels.SanityChecker.MuonChecker_cfi")
process.load("TopBrussels.SanityChecker.VertexChecker_cfi")

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('SanityChecker_test.root')
	)
	
process.p = cms.Path(
    process.makeGenEvt  
    + process.TtGenEventChecker 				#this line should be commented when not running over ttbar
    + process.Resolutions_lJets 				#this line should be commented when not running over ttbar
    + process.Resolutions_bJets 				#this line should be commented when not running over ttbar
		+ process.Resolutions_muons  				#this line should be commented when not running over ttbar  
		+ process.Resolutions_electrons    	#this line should be commented when not running over ttbar
		+ process.Resolutions_met   				#this line should be commented when not running over ttbar 
		+ process.jetmet 
    + process.muonchecker 
    + process.vertex 
    + process.jetmet 
    + process.kinematics 
    + (process.ttDecaySelection + process.TruthReco) # This line should be commented in case of not running on ttbar events
                                                     # because the ttDecaySelection will throw an edm exception thrown if no genevent is there
)

