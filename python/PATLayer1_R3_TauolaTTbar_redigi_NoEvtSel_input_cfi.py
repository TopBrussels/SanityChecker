import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
 				'/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_1.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_10.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_11.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_2.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_3.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_4.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_5.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_6.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_7.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_8.root',
        '/store/user/ghammad/TauolaTTbar/CMSSW227_R3_TauolaTTbar_redigi_nocuts/6d522cadf5d449331aa73a468b4e5e1f/PATLayer1_9.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
)
