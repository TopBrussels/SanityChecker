import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_9.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_8.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_7.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_6.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_5.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_4.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_3.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_2.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_13.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_12.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_11.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_10.root',
       '/store/user/ghammad/SUSY_LM1-sftsht/CMSSW227_R3_SUSY_LM1-sftsht_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_1.root' ] );


secFiles.extend( [
               ] )
