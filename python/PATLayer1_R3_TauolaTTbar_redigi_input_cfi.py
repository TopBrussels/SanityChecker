import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_1.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_10.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_11.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_2.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_3.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_4.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_5.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_6.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_7.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_8.root',
        '/store/user/pvmulder/TauolaTTbar/CMSSW227_R3_TauolaTTbar-new/f3f718855883784a6d2302631bb7cd61/PATLayer1_9.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
)
