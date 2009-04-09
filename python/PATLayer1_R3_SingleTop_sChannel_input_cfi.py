import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
       '/store/user/pvmulder/SingleTop_sChannel/CMSSW227_R3_SingleTop_sChannel/20cbf087f60b951e3b7f50103ab2bf6b/PATLayer1_1.root',
        '/store/user/pvmulder/SingleTop_sChannel/CMSSW227_R3_SingleTop_sChannel/20cbf087f60b951e3b7f50103ab2bf6b/PATLayer1_2.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
)
